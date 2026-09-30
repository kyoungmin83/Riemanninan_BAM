"""
Stratified cell sampling for difficulty profile bootstrap.

Goal: pick ~50k cells, balanced across celltypes, so that the per-(gene, bin,
celltype) difficulty profile is well-estimated even for rare celltypes. The
full SEA-AD DLPFC train split (~1M cells) is too large to hold all per-gene
reconstruction errors in memory; sampling gives statistical accuracy within
0.5% (median/MAD CI) at a small compute cost.

See: model_v7a_implementation_plan.md §4.2.1
"""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
import torch


def stratified_cell_sample(
    celltype_ids: Sequence[int] | np.ndarray | torch.Tensor,
    *,
    sample_size: int = 50_000,
    cells_per_celltype: Optional[int] = None,
    seed: int = 42,
) -> np.ndarray:
    """
    Sample cell indices stratified by celltype.

    Args:
        celltype_ids: per-cell celltype integer labels, length `n_cells`.
        sample_size: total cells to draw. If `cells_per_celltype` is None,
            this is divided evenly across celltypes.
        cells_per_celltype: per-celltype quota. Overrides `sample_size`-based
            quota if given. Rare celltypes contribute all available cells.
        seed: RNG seed for reproducibility.

    Returns:
        np.ndarray of cell indices (1D, int64), with one row per sampled cell.

    Notes:
        - Rare celltypes (cell count < quota) contribute all their cells.
          The leftover budget is *not* redistributed — sample_size is a target,
          not a hard floor.
        - Order within returned indices is *not* random; caller should shuffle
          if order matters.
    """
    if isinstance(celltype_ids, torch.Tensor):
        celltype_ids = celltype_ids.detach().cpu().numpy()
    elif not isinstance(celltype_ids, np.ndarray):
        celltype_ids = np.asarray(celltype_ids)
    celltype_ids = celltype_ids.astype(np.int64, copy=False)

    if celltype_ids.ndim != 1:
        raise ValueError(f"celltype_ids must be 1D, got shape {celltype_ids.shape}")

    n_cells = int(celltype_ids.shape[0])
    if n_cells == 0:
        return np.zeros((0,), dtype=np.int64)

    unique_celltypes, inverse = np.unique(celltype_ids, return_inverse=True)
    n_celltypes = int(unique_celltypes.shape[0])

    if cells_per_celltype is None:
        quota = max(1, int(sample_size // max(n_celltypes, 1)))
    else:
        quota = int(cells_per_celltype)

    rng = np.random.default_rng(seed)
    chosen: List[np.ndarray] = []
    for t_idx in range(n_celltypes):
        cell_idx = np.flatnonzero(inverse == t_idx)
        take = min(quota, cell_idx.size)
        if take == cell_idx.size:
            chosen.append(cell_idx)
        else:
            chosen.append(rng.choice(cell_idx, size=take, replace=False))

    return np.concatenate(chosen).astype(np.int64, copy=False)
