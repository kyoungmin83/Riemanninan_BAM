"""
Within-celltype state-variance floor (Phase 2).

Why
---
The binding constraint on `nonzero_exact` is the per-cell dynamic range
(R ≈ 0.38, saturated): within a celltype the model collapses its per-cell
predictions onto the celltype baseline, so cells of the same (celltype, gene)
get nearly identical tiers and exact stays near the majority-tier floor. The
asymmetric loss (Phase 1) supplies DIRECTION (push high tiers up) but not SPREAD.

This module supplies the spread: a floor that penalises the model when its
within-celltype *predicted*-tier variance falls below the within-celltype *true*-
tier variance, so it must use the per-cell latent to differentiate cells. To keep
that spread (1) in the right direction it pairs with asym, and (2) biological
rather than technical it adds a depth-decorrelation guardrail.

Target modes (see doc/within_celltype_variance_floor_design_2026-06-02.md §3.3):
  * MVP (`target_mode="in_batch"`, `nonzero_only=False`): Var_true is the in-batch
    variance of y_ord (INCLUDING zeros). Fast, but the target is noisy at small
    per-rank batches and is dominated by the zero/nonzero split rather than the
    fine nonzero 1/2/3/4 spread.
  * v1.1 (`target_mode="running"`, `nonzero_only=True`): Var_true is a STABLE
    streaming estimate (`RunningNonzeroTierVariance`, sufficient statistics over
    training batches) computed over NONZERO cells only — directly targeting the
    within-nonzero spread that drives `nonzero_exact`. Var_pred is likewise the
    nonzero-only within-celltype variance.

`group = celltype only` (no disease split — no label leakage; disease emerges as
spread). Homogeneous genes (Var_true ≈ 0) get target ≈ 0 ⇒ no penalty. Restricted
to genes with enough nonzero support per group and to celltypes with enough cells.

The trainer adds ``lambda_var * loss_var + lambda_depth_decorr * loss_depth_decorr``
only when ``enabled``; ``enabled=False`` ⇒ the block is skipped ⇒ byte-identical.
All non-loss fields are detached scalars (no GPU sync in the loss hot path).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class VarianceFloorConfig:
    enabled: bool = False
    beta: float = 0.5                 # require Var_pred >= beta * Var_true
    lambda_var: float = 0.0           # weight of the floor in the total loss
    lambda_depth_decorr: float = 0.0  # weight of the depth-decorrelation guardrail
    min_cells_per_group: int = 16     # skip celltypes with fewer cells in the batch
    min_nonzero_frac: float = 0.1     # (in_batch/total mode) genes with >= this nonzero frac
    n_bins: int = 5                   # number of ordinal tiers (0..n_bins-1); informational
    # --- v1.1 ---
    nonzero_only: bool = False        # compute Var over nonzero (y>0) cells only
    target_mode: str = "in_batch"     # "in_batch" | "running" (stable streaming target)
    min_target_count: int = 50        # (running) min accumulated nonzero cells before a (c,g) target is used
    min_nonzero_cells: int = 4        # (nonzero_only) min nonzero cells in the batch to form Var_pred
    eps: float = 1e-6


@dataclass
class VarianceFloorOutput:
    loss_var: torch.Tensor            # floor loss (grad path); 0 when no qualifying group
    loss_depth_decorr: torch.Tensor   # guardrail loss (grad path); 0 when no qualifying group
    # detached scalar diagnostics:
    var_pred_mean: torch.Tensor       # mean within-celltype predicted-tier variance (over penalised g)
    var_true_mean: torch.Tensor       # mean within-celltype true-tier variance (the target)
    var_gap: torch.Tensor             # var_pred_mean / (var_true_mean + eps): →1 as the floor closes
    active_frac: torch.Tensor         # fraction of evaluated (c,g) where the floor is still active
    state_depth_corr: torch.Tensor    # mean |corr(state magnitude, log depth)| over groups
    n_groups: torch.Tensor            # number of qualifying celltypes this batch


class RunningNonzeroTierVariance:
    """Streaming per-(celltype, gene) variance of the NONZERO ordinal tier.

    Sufficient statistics (sum, sumsq, count over y>0 cells) accumulated across
    training batches with no gradient; ``variance()`` returns the population
    within-celltype nonzero-tier variance, which converges to the true value as
    batches accumulate (a per-(c,g) target warms up once its count crosses
    ``min_target_count``). float64 buffers avoid precision loss at high counts.

    Per-rank under DDP (each rank accumulates its own and converges independently;
    no all-reduce). Not in the model state_dict ⇒ resets on resume and re-warms.
    """

    def __init__(self, n_celltypes: int, n_genes: int, device: torch.device):
        self.n_celltypes = int(n_celltypes)
        shape = (self.n_celltypes, int(n_genes))
        self.sum = torch.zeros(shape, device=device, dtype=torch.float64)
        self.sumsq = torch.zeros(shape, device=device, dtype=torch.float64)
        self.count = torch.zeros(shape, device=device, dtype=torch.float64)

    @torch.no_grad()
    def update(self, y_ord: torch.Tensor, celltype_id: torch.Tensor) -> None:
        yf = y_ord.to(torch.float64)
        nz = (y_ord > 0).to(torch.float64)               # [B, G]
        yv = yf * nz                                     # y on nonzero, 0 elsewhere
        valid = (celltype_id >= 0) & (celltype_id < self.n_celltypes)
        if not bool(valid.any()):
            return
        idx = celltype_id[valid].to(torch.long)
        self.sum.index_add_(0, idx, yv[valid])
        self.sumsq.index_add_(0, idx, (yv * yv)[valid])
        self.count.index_add_(0, idx, nz[valid])

    def variance(self) -> Tuple[torch.Tensor, torch.Tensor]:
        """Return (var [C,G] float32, count [C,G] float32). var=0 where count=0."""
        cnt = self.count.clamp_min(1.0)
        mean = self.sum / cnt
        var = (self.sumsq / cnt - mean * mean).clamp_min(0.0)
        return var.to(torch.float32), self.count.to(torch.float32)


def compute_variance_floor(
    probs: torch.Tensor,         # [B, G, K] tier probabilities
    y_ord: torch.Tensor,         # [B, G] long, true ordinal tier
    state_score: torch.Tensor,   # [B, G] per-cell state contribution to the score
    celltype_id: torch.Tensor,   # [B] long
    log_depth: torch.Tensor,     # [B] log library/depth proxy (for the guardrail)
    cfg: VarianceFloorConfig,
    var_true_target: Optional[torch.Tensor] = None,   # [C, G] stable target (running/precomputed); None ⇒ in-batch
    target_count: Optional[torch.Tensor] = None,      # [C, G] accumulated nonzero counts (for warmup masking)
) -> VarianceFloorOutput:
    """Within-celltype predicted-tier variance floor + depth-decorrelation guardrail.

    Computed in float32 (robust under bf16 autocast). Iterates unique celltypes
    (≈20), skipping any with < ``min_cells_per_group`` cells in the batch and any
    gene with insufficient (nonzero) support. Returns zeros (no gradient) when no
    group qualifies. With ``var_true_target`` provided the floor uses that stable
    target (and masks genes whose ``target_count`` < ``min_target_count``);
    otherwise it falls back to the in-batch true-tier variance."""
    dev = probs.device
    z = torch.zeros((), device=dev, dtype=torch.float32)
    eps = float(cfg.eps)
    nonzero_only = bool(cfg.nonzero_only)

    K = probs.shape[-1]
    levels = torch.arange(K, device=dev, dtype=torch.float32)
    e_pred = (probs.float() * levels).sum(dim=-1)        # [B, G] expected tier
    y_f = y_ord.to(torch.float32)                        # [B, G]
    state_mag = state_score.float().abs().mean(dim=1)    # [B] per-cell state magnitude
    ld = log_depth.to(torch.float32)                     # [B]

    var_terms = []
    dec_terms = []
    vp_sum = z.clone()
    vt_sum = z.clone()
    active_sum = z.clone()
    sdc_sum = z.clone()
    ng = 0

    for c in torch.unique(celltype_id):
        ci = int(c)
        if ci < 0:                                       # skip "unknown" celltype id
            continue
        m = celltype_id == c
        n = int(m.sum())
        if n < int(cfg.min_cells_per_group):
            continue
        ep_c = e_pred[m]                                 # [n, G]
        yt_c = y_f[m]                                    # [n, G]
        nzmask = (y_ord[m] > 0).to(torch.float32)        # [n, G]

        if nonzero_only:
            cnt = nzmask.sum(dim=0)                       # [G] nonzero cells per gene
            denom = cnt.clamp_min(1.0)
            mean_p = (ep_c * nzmask).sum(dim=0) / denom
            var_pred = (((ep_c - mean_p) ** 2) * nzmask).sum(dim=0) / denom   # [G]
            gsel = cnt >= float(cfg.min_nonzero_cells)
        else:
            var_pred = ep_c.var(dim=0, unbiased=False)   # [G]
            gsel = nzmask.mean(dim=0) >= float(cfg.min_nonzero_frac)

        # true-tier variance target
        if var_true_target is not None and ci < int(var_true_target.shape[0]):
            var_true = var_true_target[ci].to(torch.float32)           # [G] stable, detached
            if target_count is not None:
                gsel = gsel & (target_count[ci] >= float(cfg.min_target_count))
        elif nonzero_only:
            cnt = nzmask.sum(dim=0)
            denom = cnt.clamp_min(1.0)
            mean_t = (yt_c * nzmask).sum(dim=0) / denom
            var_true = (((yt_c - mean_t) ** 2) * nzmask).sum(dim=0) / denom
        else:
            var_true = yt_c.var(dim=0, unbiased=False)

        if not bool(gsel.any()):
            continue
        vp = var_pred[gsel]
        vt = var_true[gsel]
        deficit = F.relu(float(cfg.beta) * vt - vp)
        var_terms.append(deficit.pow(2).mean())

        # depth-decorrelation guardrail within this celltype
        sm = state_mag[m]
        sm_z = sm - sm.mean()
        ld_z = ld[m] - ld[m].mean()
        denom_d = (sm_z.pow(2).sum().clamp_min(eps) * ld_z.pow(2).sum().clamp_min(eps)).sqrt()
        corr = (sm_z * ld_z).sum() / denom_d.clamp_min(eps)
        dec_terms.append(corr.pow(2))

        # detached diagnostics
        vp_sum = vp_sum + vp.detach().mean()
        vt_sum = vt_sum + vt.detach().mean()
        active_sum = active_sum + (deficit.detach() > 0).to(torch.float32).mean()
        sdc_sum = sdc_sum + corr.detach().abs()
        ng += 1

    if ng == 0:
        return VarianceFloorOutput(z.clone(), z.clone(), z, z, z, z, z, z)

    ngf = float(ng)
    return VarianceFloorOutput(
        loss_var=torch.stack(var_terms).mean(),
        loss_depth_decorr=torch.stack(dec_terms).mean(),
        var_pred_mean=vp_sum / ngf,
        var_true_mean=vt_sum / ngf,
        var_gap=vp_sum / (vt_sum + eps),
        active_frac=active_sum / ngf,
        state_depth_corr=sdc_sum / ngf,
        n_groups=torch.tensor(ngf, device=dev),
    )
