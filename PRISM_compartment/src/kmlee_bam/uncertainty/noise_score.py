"""
noise_score — magnitude-only cell-level uncertainty target.

This is the value that BAM head is trained to predict (via alignment loss in
v7a). It is built by:

    1. Per-gene robust z-score using the (gene, bin, celltype) difficulty
       profile: z = (rec_err - mu) / (mad * 1.4826), then winsorized.
    2. Cell-level aggregation via trimmed-mean of |z| across genes.

The magnitude-only aggregation is *intentional* — pathology direction is
handled by separate z_state heads in v7b, not by this score.

See: model_v7a_implementation_plan.md §4.2.3, model_v7_design...md §5.0.2
"""

from __future__ import annotations

from typing import Optional

import torch

from kmlee_bam.uncertainty.difficulty_profile import DifficultyProfile


@torch.no_grad()
def normalize_rec_err(
    rec_err: torch.Tensor,           # [B, G]
    true_bins: torch.Tensor,         # [B, G] int
    celltypes: torch.Tensor,         # [B] int
    profile: DifficultyProfile,
    *,
    winsorize_z: float = 5.0,
) -> torch.Tensor:
    """Per-(cell, gene) robust z-score against the difficulty profile."""
    if rec_err.ndim != 2:
        raise ValueError(f"rec_err must be [B, G], got {tuple(rec_err.shape)}")
    if true_bins.shape != rec_err.shape:
        raise ValueError("true_bins must match rec_err shape")
    if celltypes.shape[0] != rec_err.shape[0]:
        raise ValueError("celltypes must have batch dim matching rec_err")

    device = rec_err.device
    B, G = rec_err.shape

    # Build broadcastable indices.
    gene_idx = torch.arange(G, device=device, dtype=torch.long).unsqueeze(0).expand(B, G)
    bin_idx = true_bins.long()
    celltype_idx = celltypes.long().unsqueeze(-1).expand(B, G)

    mu, mad, _ = profile.lookup_vectorized(gene_idx, bin_idx, celltype_idx)

    sigma = (mad * float(profile.config.mad_to_std)).clamp_min(
        float(profile.config.mad_floor)
    )
    z = (rec_err.float() - mu) / sigma
    z = z.clamp(-float(winsorize_z), float(winsorize_z))
    return z


def _trimmed_mean(values: torch.Tensor, *, fraction: float, dim: int = -1) -> torch.Tensor:
    """
    Symmetric trimmed mean along `dim`. Trims `fraction` from each tail.

    values: real-valued tensor
    fraction: fraction to trim from each side, e.g. 0.10 means drop 10% on
              each side, keep middle 80%.
    """
    if fraction <= 0.0:
        return values.mean(dim=dim)
    n = int(values.shape[dim])
    if n == 0:
        return values.new_zeros(values.shape[:dim] + values.shape[dim + 1:])
    k = max(0, int(round(n * float(fraction))))
    if k * 2 >= n:
        # Fall back to median when trim is too aggressive
        return values.median(dim=dim).values
    sorted_, _ = values.sort(dim=dim)
    # slice off k from each end
    sliced = sorted_.narrow(dim, k, n - 2 * k)
    return sliced.mean(dim=dim)


@torch.no_grad()
def compute_noise_score(
    rec_err: torch.Tensor,           # [B, G]
    true_bins: torch.Tensor,         # [B, G]
    celltypes: torch.Tensor,         # [B]
    profile: DifficultyProfile,
    *,
    winsorize_z: float = 5.0,
    trim_fraction: float = 0.10,
) -> torch.Tensor:
    """
    Returns: [B] noise_score (>= 0, magnitude only).
    """
    z = normalize_rec_err(
        rec_err, true_bins, celltypes, profile,
        winsorize_z=winsorize_z,
    )
    abs_z = z.abs()
    return _trimmed_mean(abs_z, fraction=trim_fraction, dim=-1)
