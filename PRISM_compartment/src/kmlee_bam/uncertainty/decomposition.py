"""
u_total = sqrt(u_cell² + u_donor²)  — composite uncertainty target.

This is the value BAM head is trained to predict. The two components are:

    u_cell  = (noise_score - bank.mean[d, t]) / bank.std[d, t]
              -> "how far is this cell from its donor's typical noise level"

    u_donor = ancova.residual_mean[d, t] / ancova.residual_std[t]
              -> "how unusually noisy is this donor (after pathology control)"

Combining them as L2 norm captures total deviation magnitude without losing
either component. Cell-level outlier and donor-level outlier both raise u_total.

See: model_v7a_implementation_plan.md §4.2.6, model_v7_design...md §5.3
"""

from __future__ import annotations

from typing import Optional, Tuple

import torch

from kmlee_bam.uncertainty.ancova import ANCOVAFit
from kmlee_bam.uncertainty.bank import NoiseScoreBank


@torch.no_grad()
def compute_u_total(
    noise_scores: torch.Tensor,    # [B]
    donors: torch.Tensor,          # [B] long
    celltypes: torch.Tensor,       # [B] long
    bank: NoiseScoreBank,
    ancova: ANCOVAFit,
    *,
    bank_std_floor: float = 0.05,
    residual_std_floor: float = 0.05,
    min_bank_n: float = 5.0,
    u_total_cap: Optional[float] = None,
    return_components: bool = False,
):
    """
    Returns u_total tensor [B]. If `return_components` is True, also returns
    (u_cell, u_donor) for diagnostics.

    Behavior under partial state:
        - For (d, t) where bank effective count `n` is below `min_bank_n`,
          u_cell is set to 0. This prevents the alignment target from being
          poisoned by undersampled bank stats (e.g. singletons where var ~ 0
          gives huge u_cell magnitudes that the BAM head can't track).
        - If bank is not initialized for (d, t), same: u_cell = 0.
        - If ANCOVA is not fitted for celltype t (e.g. too few donors),
          u_donor = 0.
        - `bank_std_floor` (default 0.05) guards against tiny EMA variances
          from single-cell batches.
    """
    device = noise_scores.device
    ns = noise_scores.to(device).float()
    donors = donors.to(device).long()
    celltypes = celltypes.to(device).long()

    bank_mean = bank.get_mean_tensor().to(device)             # [D, T]
    bank_std = bank.get_std_tensor().to(device).clamp_min(float(bank_std_floor))
    bank_init = bank.get_initialized_mask().to(device)        # [D, T] bool
    bank_n = bank.n.to(device)                                # [D, T]

    res_mean = ancova.residual_mean.to(device)                # [D, T]
    res_std = ancova.residual_std.to(device).clamp_min(float(residual_std_floor))
    is_fit = ancova.is_fitted.to(device)                      # [T] bool

    # Per-cell stats lookup.
    cell_bank_mean = bank_mean[donors, celltypes]
    cell_bank_std = bank_std[donors, celltypes]
    cell_bank_init = bank_init[donors, celltypes]             # [B] bool
    cell_bank_n = bank_n[donors, celltypes]                   # [B]

    # Gate u_cell on having enough effective samples.
    valid_cell = cell_bank_init & (cell_bank_n >= float(min_bank_n))

    u_cell = torch.zeros_like(ns)
    if bool(valid_cell.any()):
        u_cell[valid_cell] = (
            (ns[valid_cell] - cell_bank_mean[valid_cell]) / cell_bank_std[valid_cell]
        )

    # Donor effect.
    cell_res_mean = res_mean[donors, celltypes]
    cell_res_std = res_std[celltypes]
    cell_is_fit = is_fit[celltypes]                            # [B] bool

    u_donor = torch.zeros_like(ns)
    if bool(cell_is_fit.any()):
        u_donor[cell_is_fit] = (
            cell_res_mean[cell_is_fit] / cell_res_std[cell_is_fit]
        )

    u_total = (u_cell.pow(2) + u_donor.pow(2)).clamp_min(1e-12).sqrt()

    if u_total_cap is not None:
        u_total = u_total.clamp(max=float(u_total_cap))

    if return_components:
        return u_total, u_cell, u_donor
    return u_total
