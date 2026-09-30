"""
High-tier margin ranking loss (disease-program-aligned; replaces the variance floor).

Why
---
Per-cell fine-tier `exact` is noise-limited (achievable ceiling ~0.28; split-half
test-retest kappa ~0.04 on the fine/adjacent tiers). But the HIGH tier (tier 4 =
high expression) IS reproducible (split-half half-vs-full ~0.71). So instead of
ordering ALL cells (a full pairwise ranking would fit the adjacent-tier NOISE),
this loss pushes ONLY the reliable distinction: cells that are truly high-
expression (tier 4) should be predicted clearly ABOVE cells that are truly low
(tier 1/2) within the same celltype×gene. This targets tier-4 recall / high-vs-
low separation — the signal that drives disease programs & drug targets — NOT
global exact.

Three design choices that matter (see the design doc):
  * **detach(low):** gradient pushes the HIGH group UP only; the LOW group is a
    fixed (detached) reference. Pushing true-nonzero low cells DOWN would send
    them toward 0 and REVIVE the nonzero->0 leak we fight with `zero_pos_weight`.
    So low is a *baseline*, not a target.
  * **running low reference:** at the small per-rank batch (32) a given
    celltype×gene rarely has enough high AND low cells in the SAME batch, so the
    low baseline is an EMA accumulated across batches (per celltype×gene,
    detached). Any batch with enough tier-4 cells is pushed above the stored low
    level. (This is the coverage fix learned from the variance floor's n_grp=2.)
  * **tier 3 is IGNORED:** it is the ambiguous boundary in split-half; only
    tier4 (reliable) vs tier1/2 (low reference) are used.

Off (`enabled=False` or `lambda_high_margin=0`) ⇒ no term added ⇒ byte-identical.
Non-loss fields are detached scalars.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class HighMarginRankConfig:
    enabled: bool = False
    high_tier: int = 4
    low_tiers: Tuple[int, ...] = (1, 2)
    margin: float = 1.0
    lambda_high_margin: float = 0.0
    min_high_cells: int = 2          # min true-tier4 cells of a (celltype,gene) in the batch
    min_low_ref_count: int = 8       # running low-ref must have >= this many batch-updates before use
    ema_momentum: float = 0.9
    t4_floor_target: float = 0.0     # absolute floor on mean E[tier] of true-tier4 cells (0 = off)
    lambda_t4_floor: float = 0.0     # weight of the absolute tier4 floor (gradient to high cells)
    lambda_tier4_focal: float = 0.0  # tier4 binary-focal: P[tier4] UP on true-4, DOWN on true-low(1/2) → precision
    lambda_zero_focal: float = 0.0   # stable-zero focal-negative: P[tier4] DOWN on true-zeros, weighted (1-pi_tech) so only STABLE bio zeros penalized (tech-zeros spared)
    focal_gamma: float = 2.0         # focal down-weighting of easy examples
    tier4_focal_pos_weight: float = 1.0  # weight on POSITIVE focal term (P[tier4] up on true-4 = recall)
    tier4_focal_neg_weight: float = 1.0  # weight on NEGATIVE focal term (P[tier4] down on true-low = precision)
    eps: float = 1e-6


@dataclass
class HighMarginOutput:
    loss: torch.Tensor             # ranking (margin) loss; 0 when no qualifying group
    loss_t4_floor: torch.Tensor    # absolute tier4 floor loss (grad path); 0 when off / no tier4
    margin_gap: torch.Tensor       # mean(mean_high - low_ref) over active groups (detached)
    tier4_pred_level: torch.Tensor # mean predicted E[tier] of true-tier4 cells (→ should rise toward 4)
    low_pred_level: torch.Tensor   # mean low reference used
    tier4_recall: torch.Tensor     # frac of true-tier4 positions whose argmax == tier4
    active_frac: torch.Tensor      # frac of evaluated groups with deficit > 0
    n_groups: torch.Tensor         # active CELLTYPES (>=1 qualifying gene) — NOT (celltype,gene) pairs
    active_pairs: torch.Tensor       # total active (celltype,gene) pairs across celltypes
    active_high_level: torch.Tensor  # mean E[tier] of tier4 cells WITHIN active groups (vs global tier4_pred_level)
    t4_fp_on_low: torch.Tensor       # frac of true low-tier(1/2) cells predicted == tier4 (precision guardrail)
    loss_tier4_focal: torch.Tensor   # tier4 binary-focal loss (grad path; P[tier4] up on 4 / down on low); 0 when off
    loss_zero_focal: torch.Tensor    # stable-zero focal-negative loss (grad path; P[tier4] down on stable bio zeros); 0 when off
    t4_fp_on_zero: torch.Tensor      # frac of true-ZERO positions predicted == tier4 (the DOMINANT false-positive source)


class RunningLowTierScore:
    """EMA of per-(celltype, gene) mean predicted score (E[tier]) among LOW-tier
    cells. Detached reference for the high-margin loss; lets tier-4 cells in any
    batch be pushed above a stable low baseline (coverage at batch=32). Per-rank
    under DDP; not in the model state_dict (resets on resume, re-warms)."""

    def __init__(self, n_celltypes: int, n_genes: int, device, momentum: float = 0.9):
        self.n_celltypes = int(n_celltypes)
        self.momentum = float(momentum)
        shape = (self.n_celltypes, int(n_genes))
        self.ema = torch.zeros(shape, device=device, dtype=torch.float32)
        self.updates = torch.zeros(shape, device=device, dtype=torch.float32)

    @torch.no_grad()
    def update(self, score: torch.Tensor, low_mask: torch.Tensor, celltype_id: torch.Tensor) -> None:
        """score [B,G] (E[tier]); low_mask [B,G] float (1 at true low-tier cells)."""
        s_low = score.float() * low_mask
        valid = (celltype_id >= 0) & (celltype_id < self.n_celltypes)
        if not bool(valid.any()):
            return
        idx = celltype_id[valid].long()
        sum_s = torch.zeros_like(self.ema).index_add_(0, idx, s_low[valid])
        cnt = torch.zeros_like(self.ema).index_add_(0, idx, low_mask[valid].float())
        has = cnt > 0
        batch_mean = sum_s / cnt.clamp_min(1.0)
        ema_upd = self.momentum * self.ema + (1.0 - self.momentum) * batch_mean
        # EMA where already initialised & has data; init to batch_mean where fresh; keep otherwise
        self.ema = torch.where(has & (self.updates > 0), ema_upd,
                               torch.where(has, batch_mean, self.ema))
        self.updates = self.updates + has.to(torch.float32)

    def reference(self) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.ema, self.updates


def high_margin_loss(
    score: torch.Tensor,        # [B,G] predicted E[tier] (carries gradient)
    y_ord: torch.Tensor,        # [B,G] long
    celltype_id: torch.Tensor,  # [B] long
    low_ref: torch.Tensor,      # [C,G] detached running low baseline
    low_ref_updates: torch.Tensor,  # [C,G]
    probs: torch.Tensor,        # [B,G,K] (for the tier4_recall diagnostic)
    cfg: HighMarginRankConfig,
    pi_tech: Optional[torch.Tensor] = None,  # [B,G] tech-dropout prob; (1-pi_tech) weights stable-zero negatives
) -> HighMarginOutput:
    dev = score.device
    z = torch.zeros((), device=dev, dtype=torch.float32)
    sc = score.float()
    high = y_ord == int(cfg.high_tier)              # [B,G] bool

    terms = []
    gap_sum = z.clone(); lo_sum = z.clone(); act_sum = z.clone()
    hi_sum = z.clone(); pair_sum = z.clone()
    ng = 0
    for c in torch.unique(celltype_id):
        ci = int(c)
        if ci < 0:
            continue
        m = celltype_id == c
        hmask = high[m].to(torch.float32)            # [n, G]
        cnt_hi = hmask.sum(0)                          # [G]
        gsel = (cnt_hi >= float(cfg.min_high_cells)) & \
               (low_ref_updates[ci] >= float(cfg.min_low_ref_count))
        if not bool(gsel.any()):
            continue
        sc_c = sc[m]                                   # [n, G]
        # mean predicted score of the tier-4 cells (gradient flows ONLY here)
        mean_high = (sc_c * hmask).sum(0) / cnt_hi.clamp_min(1.0)   # [G]
        lref = low_ref[ci].detach()                    # [G] detached baseline (no grad to low)
        mh = mean_high[gsel]
        lr = lref[gsel]
        deficit = F.relu(float(cfg.margin) - (mh - lr))
        terms.append(deficit.pow(2).mean())
        gap_sum = gap_sum + (mh.detach() - lr).mean()
        lo_sum = lo_sum + lr.mean()
        hi_sum = hi_sum + mh.detach().mean()
        act_sum = act_sum + (deficit.detach() > 0).to(torch.float32).mean()
        pair_sum = pair_sum + gsel.sum().to(torch.float32)
        ng += 1

    # absolute tier4 floor (direct upward pull on true-tier4 E[tier]) + diagnostics over
    # ALL true-tier4 positions. Coverage-free (global), so it works even when no
    # per-(c,g) ranking group qualifies.
    pred_tier = probs.detach().argmax(dim=-1)
    if bool(high.any()):
        t4_mean = sc[high].mean()                              # WITH gradient
        t4_level = t4_mean.detach()
        t4_recall = (pred_tier[high] == int(cfg.high_tier)).to(torch.float32).mean()
        if float(cfg.lambda_t4_floor) > 0.0:
            loss_t4_floor = F.relu(float(cfg.t4_floor_target) - t4_mean).pow(2)   # grad -> high
        else:
            loss_t4_floor = z.clone()
    else:
        t4_level = z
        t4_recall = z
        loss_t4_floor = z.clone()

    # precision guardrail: fraction of true LOW-tier (1/2) cells wrongly predicted == tier4
    low_any = torch.zeros_like(y_ord, dtype=torch.bool)
    for lt in cfg.low_tiers:
        low_any = low_any | (y_ord == int(lt))
    if bool(low_any.any()):
        t4_fp_on_low = (pred_tier[low_any] == int(cfg.high_tier)).to(torch.float32).mean()
    else:
        t4_fp_on_low = z

    # tier4 binary-focal (P[tier4] UP on true-4, DOWN on true-low 1/2) + stable-zero
    # focal-negative (P[tier4] DOWN on true-ZEROS, weighted by (1-pi_tech) so only
    # confident BIOLOGICAL zeros are penalized, tech-zeros spared). The negative terms
    # are what the margin/floor lack -> they raise PRECISION. tier3 is always ignored.
    zero_mask = (y_ord == 0)
    if bool(zero_mask.any()):
        t4_fp_on_zero = (pred_tier[zero_mask] == int(cfg.high_tier)).to(torch.float32).mean()
    else:
        t4_fp_on_zero = z
    loss_tier4_focal = z.clone()
    loss_zero_focal = z.clone()
    if float(cfg.lambda_tier4_focal) > 0.0 or float(cfg.lambda_zero_focal) > 0.0:
        p4 = probs[..., int(cfg.high_tier)].float().clamp(float(cfg.eps), 1.0 - float(cfg.eps))
        gfoc = float(cfg.focal_gamma)
        if float(cfg.lambda_tier4_focal) > 0.0:
            pos_l = z.clone()
            neg_l = z.clone()
            if bool(high.any()):
                pp = p4[high]
                pos_l = (-((1.0 - pp) ** gfoc) * torch.log(pp)).mean()
            if bool(low_any.any()):
                pn = p4[low_any]
                neg_l = (-(pn ** gfoc) * torch.log1p(-pn)).mean()
            loss_tier4_focal = (float(cfg.tier4_focal_pos_weight) * pos_l
                                + float(cfg.tier4_focal_neg_weight) * neg_l)
        if float(cfg.lambda_zero_focal) > 0.0 and pi_tech is not None and bool(zero_mask.any()):
            w = (1.0 - pi_tech.float()).clamp(0.0, 1.0)   # stable-zero weight: low pi_tech -> ~1
            pz = p4[zero_mask]
            wz = w[zero_mask]
            denom = wz.sum().clamp_min(1.0)
            loss_zero_focal = (wz * (pz ** gfoc) * (-torch.log1p(-pz))).sum() / denom

    if ng == 0:
        return HighMarginOutput(z.clone(), loss_t4_floor, z, t4_level, z, t4_recall, z, z,
                                z, z, t4_fp_on_low, loss_tier4_focal, loss_zero_focal, t4_fp_on_zero)

    ngf = float(ng)
    return HighMarginOutput(
        loss=torch.stack(terms).mean(),
        loss_t4_floor=loss_t4_floor,
        margin_gap=gap_sum / ngf,
        tier4_pred_level=t4_level,
        low_pred_level=lo_sum / ngf,
        tier4_recall=t4_recall,
        active_frac=act_sum / ngf,
        n_groups=torch.tensor(ngf, device=dev),
        active_pairs=pair_sum,
        active_high_level=hi_sum / ngf,
        t4_fp_on_low=t4_fp_on_low,
        loss_tier4_focal=loss_tier4_focal,
        loss_zero_focal=loss_zero_focal,
        t4_fp_on_zero=t4_fp_on_zero,
    )
