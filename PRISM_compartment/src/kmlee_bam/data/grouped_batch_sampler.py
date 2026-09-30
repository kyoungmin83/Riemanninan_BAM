# KMLEE-BAM grouped (donor x celltype) block batch sampler
from __future__ import annotations

from collections import defaultdict
from typing import Iterator

import math
import numpy as np
from torch.utils.data import Sampler


class GroupedBatchSampler(Sampler[list[int]]):
    """Batch sampler that builds each batch from a few DENSE (donor, celltype) blocks.

    Each batch = ``n_blocks`` distinct (celltype, donor) groups x ``block_size`` cells
    (n_blocks * block_size == batch_size). This makes the v32 stochastic-consistency loss
    meaningful: with random batching ~89% of (donor,celltype) groups in a batch are singletons,
    so the donor x celltype group-mean is a single cell (per-cell noise, filtered out). Dense
    blocks give the consistency real multi-cell aggregates (e.g. 16 groups x 8 cells).

    Design
    ------
    - Yields dataset-LOCAL indices (like ReferenceAnchorBatchSampler), not zarr rows.
    - Group selection is weighted by group size and a fixed ``block_size`` is drawn from each
      chosen group -> P(cell per batch) ~ size * (block_size / size) = const => UNIFORM per-cell
      coverage (same expected coverage as random sampling), while keeping blocks dense.
    - min_group_cells defaults to block_size: every kept group has >= block_size cells, so each
      block is block_size REAL distinct cells (no replacement padding). Lowering it admits small
      groups (a < block_size group is then padded with replacement). Groups below the threshold
      are excluded from the grouped runs entirely (both treatment+control, so the pair stays matched).
    - DDP: each rank uses a rank-specific RNG and the SAME number of batches per epoch
      (so the backward all_reduce stays in lockstep).
    - Safe with v31a: no BatchNorm in the model; ordinal_balance / thinning / subgroup_cvar /
      donor_balanced_ref_center are loss-level; reference_anchor_sampler is OFF.
    """

    def __init__(
        self,
        dataset,
        *,
        batch_size: int,
        block_size: int = 8,
        num_replicas: int = 1,
        rank: int = 0,
        seed: int = 42,
        drop_last: bool = False,
        min_group_cells: int | None = None,
    ) -> None:
        if batch_size <= 0:
            raise ValueError("batch_size must be positive.")
        if block_size <= 0:
            raise ValueError("block_size must be positive.")
        if batch_size % block_size != 0:
            raise ValueError(
                f"batch_size ({batch_size}) must be divisible by block_size ({block_size}) "
                "so every batch is exactly n_blocks * block_size cells (no silent shrink)."
            )
        if num_replicas <= 0:
            raise ValueError("num_replicas must be positive.")
        if not (0 <= rank < num_replicas):
            raise ValueError(f"rank must be in [0, {num_replicas}), got {rank}.")

        required_attrs = ["row_idx", "celltype_ids", "donor_ids"]
        missing = [a for a in required_attrs if not hasattr(dataset, a)]
        if missing:
            raise AttributeError(
                f"GroupedBatchSampler requires dataset attributes {required_attrs}. "
                f"Missing: {missing}"
            )

        self.dataset = dataset
        self.batch_size = int(batch_size)
        self.block_size = int(block_size)
        self.n_blocks = max(1, self.batch_size // self.block_size)
        self.effective_batch = self.n_blocks * self.block_size  # == batch_size if divisible
        self.num_replicas = int(num_replicas)
        self.rank = int(rank)
        self.seed = int(seed)
        self.drop_last = bool(drop_last)
        # default min_group_cells = block_size -> every kept group has >= block_size cells, so each
        # block is block_size REAL distinct cells (no replacement padding), which is what the
        # consistency wants. Lower it to admit small groups (then padded with replacement).
        self.min_group_cells = int(block_size if min_group_cells is None else min_group_cells)
        self.epoch = 0

        self.n = len(dataset)

        # DDP-like epoch length (identical on every rank -> lockstep backward).
        local_n = self.n / float(self.num_replicas)
        if self.drop_last:
            self.num_batches = max(1, int(math.floor(local_n / self.batch_size)))
        else:
            self.num_batches = max(1, int(math.ceil(local_n / self.batch_size)))

        self._build_group_index()

    def _build_group_index(self) -> None:
        """group_cells[g] = dataset-local indices of cells sharing a (celltype, donor)."""
        row_idx = np.asarray(self.dataset.row_idx, dtype=np.int64)
        celltype_ids = np.asarray(self.dataset.celltype_ids, dtype=np.int64)
        donor_ids = np.asarray(self.dataset.donor_ids, dtype=np.int64)

        groups: dict[tuple[int, int], list[int]] = defaultdict(list)
        for local_i, row in enumerate(row_idx):
            r = int(row)
            groups[(int(celltype_ids[r]), int(donor_ids[r]))].append(int(local_i))

        kept = [
            np.asarray(v, dtype=np.int64)
            for v in groups.values()
            if len(v) >= self.min_group_cells
        ]
        self.n_groups_total = int(len(groups))
        self.n_groups_excluded = int(self.n_groups_total - len(kept))
        self.n_cells_eligible = int(sum(len(group) for group in kept))
        self.n_cells_excluded = int(self.n - self.n_cells_eligible)
        if len(kept) < self.n_blocks:
            raise ValueError(
                f"GroupedBatchSampler: only {len(kept)} (celltype,donor) groups have "
                f">= {self.min_group_cells} cells, need >= n_blocks={self.n_blocks}."
            )
        self.group_cells = kept
        sizes = np.asarray([len(g) for g in kept], dtype=np.float64)
        self.group_probs = sizes / sizes.sum()  # weight by size => uniform per-cell coverage
        self.n_groups = len(kept)

    def set_epoch(self, epoch: int) -> None:
        self.epoch = int(epoch)

    def __len__(self) -> int:
        return int(self.num_batches)

    def __iter__(self) -> Iterator[list[int]]:
        rng = np.random.default_rng(
            self.seed + 100_003 * int(self.epoch) + 9_973 * int(self.rank)
        )
        group_idx = np.arange(self.n_groups)
        for _ in range(self.num_batches):
            chosen = rng.choice(
                group_idx, size=self.n_blocks, replace=False, p=self.group_probs
            )
            batch: list[int] = []
            for gi in chosen.tolist():
                pool = self.group_cells[gi]
                cells = rng.choice(
                    pool, size=self.block_size, replace=(len(pool) < self.block_size)
                )
                batch.extend(int(x) for x in cells.tolist())
            yield batch
