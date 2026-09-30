# KMLEE-BAM CP-BAM: Consensus-Private decomposition head (CORE v1)
#
# Splits the disease latent z_perp into:
#   z_common  -> the cross-donor COMMON disease fingerprint (predicts pathology, donor-invariant
#                within a pathology x celltype context, pulled toward a leave-one-donor consensus)
#   z_private -> the donor-PRIVATE residual (shrunk, kept for precision-medicine fingerprints)
#
# Faithful CORE of the CP-BAM design (sections 3-6, 9). Stage-2 decoder-aware
# losses are configured here but attached inside model/system.py so decoder
# parameters remain DDP-visible.
#
# DEFERRED to later stages (NOT here):
#   - disease-programme mixture D = sum_k pi_k B_k (low-rank)
#   - consensus prior on the BAM *attention* (q_d -> p_pop^{-d})
#   - group-DRO / split-half-as-loss
# This module is self-contained and GATED (cfg.enabled); it does NOT touch the encoder, decoder,
# or the existing pathology_aux. It reads z_perp and returns its own loss.
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


DEFAULT_PATHOLOGY_AXES: Tuple[str, ...] = (
    "axis_braak_high",
    "axis_thal_high",
    "axis_late_positive",
    "axis_lewy_positive",
)


@dataclass
class CPBAMConfig:
    enabled: bool = False
    d_z: int = 32
    hidden_dim: int = 0  # 0 => linear projections/heads
    axes: Tuple[str, ...] = DEFAULT_PATHOLOGY_AXES
    n_donor: int = 0          # required (number of donors) for the donor adversary
    n_celltype: int = 0       # required for the celltype embedding
    min_group_cells: int = 2  # (donor,celltype) groups smaller than this are masked out
    tau_y: float = 1.0        # pathology-similarity kernel temperature
    ema_decay: float = 0.99   # EMA decay for the per-(donor,celltype) z_common consensus bank
    donor_adv_strength: float = 1.0   # GRL strength (ramped externally if desired)
    celltype_embed_dim: int = 8
    # loss weights
    lambda_path: float = 1.0
    lambda_consensus: float = 0.3
    lambda_donor_adv: float = 0.05
    lambda_orth: float = 0.1
    lambda_shrink: float = 0.01
    lambda_latent_rec: float = 1.0
    # Stage-2 decoder-aware decomposition weights (default OFF).
    #
    # lambda_common_decode:
    #   decode(z_common) -> y_ord NLL. Small values make z_common a decoder-native
    #   population programme instead of a pure auxiliary projection.
    # lambda_split_full_decode:
    #   decode(z_common + z_private) -> y_ord NLL. Ensures the split remains
    #   directly usable by the existing ordinal decoder.
    # lambda_private_resid_score:
    #   decode(stopgrad(z_common) + z_private) should match the original full
    #   decoder score, forcing z_private to fill residual variation rather than
    #   rewrite the common branch.
    # lambda_ref_zero:
    #   on reference-clean cells, z_common/z_private and their decoder state
    #   scores are pulled toward zero, aligning the split with normal reference.
    lambda_common_decode: float = 0.0
    lambda_split_full_decode: float = 0.0
    lambda_private_resid_score: float = 0.0
    lambda_ref_zero: float = 0.0


@dataclass
class CPBAMOutput:
    loss: torch.Tensor
    details: Dict[str, float] = field(default_factory=dict)
    z_common: Optional[torch.Tensor] = None
    z_private: Optional[torch.Tensor] = None


class _GradReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, strength):
        ctx.strength = float(strength)
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad):
        return -ctx.strength * grad, None


def grad_reverse(x: torch.Tensor, strength: float) -> torch.Tensor:
    return _GradReverse.apply(x, strength)


def _mlp(d_in: int, d_out: int, hidden: int) -> nn.Module:
    if hidden and hidden > 0:
        return nn.Sequential(nn.Linear(d_in, hidden), nn.GELU(), nn.Linear(hidden, d_out))
    return nn.Linear(d_in, d_out)


class CPBAMHead(nn.Module):
    """Consensus-Private decomposition of z_perp. See module docstring."""

    def __init__(self, config: CPBAMConfig):
        super().__init__()
        self.cfg = config
        d_z = int(config.d_z)
        h = int(config.hidden_dim)
        self.n_axes = len(config.axes)
        self.axes = tuple(config.axes)

        # --- the split: two parallel projections of z_perp -------------------------------------
        self.common_proj = _mlp(d_z, d_z, h)
        self.private_proj = _mlp(d_z, d_z, h)

        # --- pathology read-out from z_common (group-level) ------------------------------------
        self.path_head = _mlp(d_z, self.n_axes, h)

        # --- conditional donor adversary on z_common (GRL + celltype + pathology context) ------
        self.ct_embed = nn.Embedding(max(1, int(config.n_celltype)), int(config.celltype_embed_dim))
        adv_in = d_z + int(config.celltype_embed_dim) + self.n_axes
        adv_hidden = h if h and h > 0 else max(64, 2 * d_z)
        self.donor_adv = nn.Sequential(
            nn.Linear(adv_in, adv_hidden), nn.GELU(),
            nn.Linear(adv_hidden, max(1, int(config.n_donor))),
        )

        # --- consensus EMA bank: per (donor, celltype) z_common aggregate + pathology vector ----
        # detached running state (NOT synced across DDP ranks in CORE v1 -> each rank's bank
        # converges to cover the donors it sees; stage-2 could all-reduce it).
        self._n_donor_bank = max(1, int(config.n_donor))
        self._n_celltype_bank = max(1, int(config.n_celltype))
        self.register_buffer("bank_zc", torch.zeros(self._n_donor_bank, self._n_celltype_bank, d_z), persistent=False)
        self.register_buffer("bank_y", torch.zeros(self._n_donor_bank, self._n_celltype_bank, self.n_axes), persistent=False)
        self.register_buffer("bank_seen", torch.zeros(self._n_donor_bank, self._n_celltype_bank), persistent=False)

    # -------------------------------------------------------------------------------------------
    @staticmethod
    def _group(donor_id: torch.Tensor, celltype_id: torch.Tensor):
        """Composite (donor, celltype) grouping. Returns (inv, G, cnt)."""
        key = donor_id.reshape(-1).long() * 100000 + celltype_id.reshape(-1).long()
        uniq, inv = torch.unique(key, return_inverse=True)
        G = int(uniq.numel())
        return inv, G

    def forward(
        self,
        z_perp: torch.Tensor,                 # [B, d_z]  (the existing disease latent)
        donor_id: torch.Tensor,               # [B]
        celltype_id: torch.Tensor,            # [B]
        labels: Dict[str, torch.Tensor],      # axis_* -> [B] (donor-level, repeated per cell)
    ) -> CPBAMOutput:
        cfg = self.cfg
        B = z_perp.shape[0]
        dev = z_perp.device
        details: Dict[str, float] = {}

        # ===== 1. split =====================================================================
        z_common = self.common_proj(z_perp)
        z_private = self.private_proj(z_perp)
        # latent reconstruction: keep the split lossless w.r.t. z_perp (no decoder change).
        loss_latent = ((z_perp - z_common - z_private) ** 2).mean()

        # per-cell pathology vector y (float) — donor-level labels repeated per cell.
        y = torch.stack([labels[a].reshape(-1).float() for a in self.axes], dim=1)  # [B, n_axes]

        # ===== group by (donor, celltype) ===================================================
        inv, G = self._group(donor_id, celltype_id)
        ones = z_common.new_ones(B)
        cnt = z_common.new_zeros(G).index_add_(0, inv, ones)                      # [G]

        def gmean(t):  # differentiable group mean
            s = t.new_zeros(G, t.shape[1]).index_add_(0, inv, t)
            return s / cnt.clamp_min(1.0).unsqueeze(1)

        zc_g = gmean(z_common)                                                    # [G, d_z]
        y_g = gmean(y)                                                            # [G, n_axes]
        # donor / celltype per group (constant within a group)
        donor_g = donor_id.new_zeros(G); donor_g[inv] = donor_id.long()
        ct_g = celltype_id.new_zeros(G); ct_g[inv] = celltype_id.long()
        keep = (cnt >= float(cfg.min_group_cells))                               # [G] bool
        n_keep = int(keep.sum().item())

        # ===== 2. pathology from z_common (group level) =====================================
        path_logits = self.path_head(zc_g)                                       # [G, n_axes]
        if n_keep > 0:
            loss_path = F.binary_cross_entropy_with_logits(
                path_logits[keep], y_g[keep].clamp(0.0, 1.0)
            )
        else:
            loss_path = z_perp.new_zeros(())

        # ===== 3. leave-one-donor-out consensus via EMA bank (latent z_common) ===============
        # Pull each (donor,celltype) z_common aggregate toward the pathology-similarity-weighted
        # consensus of OTHER donors of the SAME celltype, read from the EMA bank (ALL donors, not
        # just this batch). The consensus target is stop-grad (bank is no_grad); the CURRENT donor
        # is excluded. Grad flows only through zc_g -> z_common -> common_proj -> z_perp.
        n_c = self._n_celltype_bank
        ctg = ct_g.long().clamp(0, n_c - 1)
        dg = donor_g.long().clamp(0, self._n_donor_bank - 1)
        with torch.no_grad():
            bz = self.bank_zc[:, ctg, :]                  # [n_donor, G, d_z]
            by = self.bank_y[:, ctg, :]                   # [n_donor, G, n_axes]
            bs = self.bank_seen[:, ctg]                   # [n_donor, G]
            d2 = ((y_g.unsqueeze(0) - by) ** 2).sum(-1)   # [n_donor, G]
            Kw = torch.exp(-d2 / max(cfg.tau_y, 1e-6))    # [n_donor, G]
            dvec = torch.arange(self._n_donor_bank, device=dev).unsqueeze(1)
            diff = (dvec != dg.unsqueeze(0)).float()      # [n_donor, G]  exclude current donor
            W = Kw * bs * diff                            # [n_donor, G]
            row = W.sum(0)                                # [G]
            consensus = (W.unsqueeze(-1) * bz).sum(0) / row.clamp_min(1e-6).unsqueeze(-1)  # [G, d_z]
            valid = (row > 0) & keep                      # groups that have bank peers
        n_valid = int(valid.sum().item())
        if n_valid > 0:
            loss_consensus = ((zc_g[valid] - consensus[valid]) ** 2).mean()   # consensus is no_grad
        else:
            loss_consensus = z_perp.new_zeros(())
        # Update the EMA bank ONLY on a real TRAINING step, AFTER using it (consensus uses the
        # PREVIOUS state). self.training + grad-enabled gates out eval/validation/test forwards so
        # validation donors NEVER enter the train consensus bank — essential for the donor-disjoint
        # generalisation CP-BAM is built for.
        if self.training and torch.is_grad_enabled():
            with torch.no_grad():
                decay = float(cfg.ema_decay)
                g_flat = (dg * n_c + ctg).clamp(0, self._n_donor_bank * n_c - 1)
                bzc = self.bank_zc.view(-1, z_common.shape[1])
                byv = self.bank_y.view(-1, self.n_axes)
                bsv = self.bank_seen.view(-1)
                seen_now = bsv[g_flat] > 0
                new = zc_g.detach().to(bzc.dtype)  # bank is float32; z_common may be bfloat16 (autocast)
                bzc[g_flat] = torch.where(
                    seen_now.unsqueeze(1), decay * bzc[g_flat] + (1.0 - decay) * new, new
                )
                byv[g_flat] = y_g.detach().to(byv.dtype)
                bsv[g_flat] = 1.0

        # ===== 4. conditional donor adversary on z_common (per cell) ========================
        # GRL: z_common is pushed to NOT reveal donor given (celltype, pathology) context.
        z_grl = grad_reverse(z_common, cfg.donor_adv_strength)
        # NOT clamp_(): celltype_id is already long, so .long() returns the SAME batch tensor and an
        # in-place clamp_ would mutate batch["celltype_id"] (downstream code reads it). Use clamp().
        ct_idx = celltype_id.long().clamp(0, self.ct_embed.num_embeddings - 1)
        ct_e = self.ct_embed(ct_idx)
        adv_in = torch.cat([z_grl, ct_e, y], dim=1)
        donor_logits = self.donor_adv(adv_in)                                    # [B, n_donor]
        # balanced per-batch weighting (rare donors not ignored)
        with torch.no_grad():
            uid, ucnt = torch.unique(donor_id.long(), return_counts=True)
            # float32: under autocast donor_logits is bfloat16; a bfloat16 wtab indexed by a float32
            # source raises an index_put dtype error, and CE wants a float weight. Keep both float32.
            wtab = torch.zeros(donor_logits.shape[1], device=donor_logits.device, dtype=torch.float32)
            wtab[uid.clamp(0, wtab.numel() - 1)] = (B / (uid.numel() * ucnt.float())).clamp_max(50.0)
        loss_donor_adv = F.cross_entropy(donor_logits.float(), donor_id.long(), weight=wtab)
        with torch.no_grad():
            donor_acc = (donor_logits.argmax(1) == donor_id.long()).float().mean()

        # ===== 5. orthogonality: z_common <-> z_private cross-covariance =====================
        zc = z_common - z_common.mean(0, keepdim=True)
        zp = z_private - z_private.mean(0, keepdim=True)
        cov = (zc.transpose(0, 1) @ zp) / max(B, 1)                              # [d_z, d_z]
        loss_orth = (cov ** 2).mean()

        # ===== 6. private shrinkage (L1) ====================================================
        loss_shrink = z_private.abs().mean()

        # ===== total ========================================================================
        total = (
            cfg.lambda_path * loss_path
            + cfg.lambda_consensus * loss_consensus
            + cfg.lambda_donor_adv * loss_donor_adv
            + cfg.lambda_orth * loss_orth
            + cfg.lambda_shrink * loss_shrink
            + cfg.lambda_latent_rec * loss_latent
        )

        def _f(t):
            try:
                return float(t.detach().cpu())
            except Exception:
                return float("nan")

        details = {
            "loss/cpbam_path": _f(loss_path),
            "loss/cpbam_consensus": _f(loss_consensus),
            "loss/cpbam_donor_adv": _f(loss_donor_adv),
            "loss/cpbam_orth": _f(loss_orth),
            "loss/cpbam_shrink": _f(loss_shrink),
            "loss/cpbam_latent_rec": _f(loss_latent),
            "loss/cpbam_total": _f(total),
            "metric/cpbam_donor_adv_acc": _f(donor_acc),       # want LOW (donor not decodable)
            "metric/cpbam_zcommon_norm": _f(z_common.norm(dim=1).mean()),
            "metric/cpbam_zprivate_norm": _f(z_private.norm(dim=1).mean()),
            "metric/cpbam_n_groups": float(G),
            "metric/cpbam_n_kept": float(n_keep),
            "metric/cpbam_valid_consensus": float(n_valid),            # groups that found bank peers
            "metric/cpbam_bank_seen_frac": _f(self.bank_seen.mean()),  # bank coverage (rank-local, #2)
        }
        return CPBAMOutput(loss=total, details=details, z_common=z_common, z_private=z_private)
