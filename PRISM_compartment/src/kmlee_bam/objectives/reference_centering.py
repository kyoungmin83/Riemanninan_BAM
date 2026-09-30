from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch


def _zero_like(x: torch.Tensor) -> torch.Tensor:
    return torch.zeros((), dtype=x.dtype, device=x.device)


@dataclass(frozen=True)
class DonorBalancedReferenceCenterConfig:
    min_cells_per_donor: int = 2
    min_donors_per_celltype: int = 2


@dataclass(frozen=True)
class DonorBalancedReferenceCenterStats:
    loss: torch.Tensor
    active_fraction: torch.Tensor
    n_reference_celltypes: int
    n_active_celltypes: int


def donor_balanced_reference_center_loss(
    z_perp: Optional[torch.Tensor],
    celltype_id: torch.Tensor,
    donor_id: torch.Tensor,
    is_reference: Optional[torch.Tensor],
    *,
    config: DonorBalancedReferenceCenterConfig = DonorBalancedReferenceCenterConfig(),
) -> torch.Tensor:
    """
    Penalise non-zero reference origin using donor means, not cell means.

    For each cell type c:

        mean_c = mean_d mean_{i: c_i=c, donor_i=d, reference_i} z_perp_i

    The loss is mean_c ||mean_c||^2.

    This prevents donors with many cells from dominating the reference-normal
    origin and is the first v3 fix for train-reference origin transfer.
    """
    return donor_balanced_reference_center_stats(
        z_perp=z_perp,
        celltype_id=celltype_id,
        donor_id=donor_id,
        is_reference=is_reference,
        config=config,
    ).loss


def donor_balanced_reference_center_stats(
    z_perp: Optional[torch.Tensor],
    celltype_id: torch.Tensor,
    donor_id: torch.Tensor,
    is_reference: Optional[torch.Tensor],
    *,
    config: DonorBalancedReferenceCenterConfig = DonorBalancedReferenceCenterConfig(),
) -> DonorBalancedReferenceCenterStats:
    """
    Donor-balanced reference center loss plus minibatch activation diagnostics.

    `active_fraction` is the fraction of reference-containing cell types in the
    microbatch that had enough donor-balanced support to contribute a non-zero
    loss term. With small microbatches this is as important as the loss value.
    """
    if z_perp is None:
        zero = _zero_like(celltype_id.float())
        return DonorBalancedReferenceCenterStats(zero, zero, 0, 0)
    if is_reference is None:
        zero = _zero_like(z_perp)
        return DonorBalancedReferenceCenterStats(zero, zero, 0, 0)
    if z_perp.ndim != 2:
        raise ValueError(f"z_perp must have shape [B, d_z], got {tuple(z_perp.shape)}.")
    if celltype_id.ndim != 1 or donor_id.ndim != 1:
        raise ValueError("celltype_id and donor_id must have shape [B].")
    if celltype_id.shape[0] != z_perp.shape[0] or donor_id.shape[0] != z_perp.shape[0]:
        raise ValueError("Batch size mismatch among z_perp, celltype_id, and donor_id.")

    is_reference = is_reference.bool()
    if is_reference.ndim != 1 or is_reference.shape[0] != z_perp.shape[0]:
        raise ValueError("is_reference must have shape [B].")
    if int(is_reference.sum().item()) == 0:
        zero = _zero_like(z_perp)
        return DonorBalancedReferenceCenterStats(zero, zero, 0, 0)

    losses = []
    reference_celltypes = torch.unique(celltype_id[is_reference])
    n_reference_celltypes = int(reference_celltypes.numel())
    n_active_celltypes = 0

    for c in reference_celltypes:
        donor_means = []
        c_mask = is_reference & (celltype_id == c)
        for d in torch.unique(donor_id[c_mask]):
            m = c_mask & (donor_id == d)
            if int(m.sum().item()) < int(config.min_cells_per_donor):
                continue
            donor_means.append(z_perp[m].mean(dim=0))

        if len(donor_means) < int(config.min_donors_per_celltype):
            continue

        n_active_celltypes += 1
        mean_c = torch.stack(donor_means, dim=0).mean(dim=0)
        losses.append(mean_c.pow(2).mean())

    if not losses:
        loss = _zero_like(z_perp)
    else:
        loss = torch.stack(losses).mean()

    active_fraction = torch.tensor(
        float(n_active_celltypes) / max(float(n_reference_celltypes), 1.0),
        dtype=z_perp.dtype,
        device=z_perp.device,
    )
    return DonorBalancedReferenceCenterStats(
        loss=loss,
        active_fraction=active_fraction,
        n_reference_celltypes=n_reference_celltypes,
        n_active_celltypes=n_active_celltypes,
    )



def group_relative_uncertainty(
    uncertainty: torch.Tensor,
    group_id: torch.Tensor,
    *,
    min_group_size: int = 4,
    robust: bool = False,
    eps: float = 1e-6,
) -> torch.Tensor:
    """
    Convert absolute uncertainty into within-group relative uncertainty.

    The group can be celltype, celltype x tech, or celltype x tech x depth_bin.
    For small groups, the original batch-wide normalisation is used.

    This is intentionally simple and batch-local. For long training, prefer an
    EMA table keyed by celltype/tech/depth once the audit confirms which
    confounders dominate raw BAM entropy.
    """
    if uncertainty.ndim != 1 or group_id.ndim != 1:
        raise ValueError("uncertainty and group_id must both have shape [B].")
    if uncertainty.shape[0] != group_id.shape[0]:
        raise ValueError("Batch size mismatch between uncertainty and group_id.")

    out = torch.empty_like(uncertainty)
    global_center = uncertainty.median() if robust else uncertainty.mean()
    if robust:
        global_scale = (uncertainty - global_center).abs().median().clamp_min(eps)
    else:
        global_scale = uncertainty.std(unbiased=False).clamp_min(eps)

    for g in torch.unique(group_id):
        m = group_id == g
        if int(m.sum().item()) < int(min_group_size):
            out[m] = (uncertainty[m] - global_center) / global_scale
            continue
        center = uncertainty[m].median() if robust else uncertainty[m].mean()
        if robust:
            scale = (uncertainty[m] - center).abs().median().clamp_min(eps)
        else:
            scale = uncertainty[m].std(unbiased=False).clamp_min(eps)
        out[m] = (uncertainty[m] - center) / scale

    return out
