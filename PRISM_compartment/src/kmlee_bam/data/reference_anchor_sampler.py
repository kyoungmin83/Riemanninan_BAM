# KMLEE-BAM reference-anchor sampler
from __future__ import annotations

from collections import defaultdict
from typing import Iterator, Optional

import math
import numpy as np
import torch
from torch.utils.data import Sampler


class ReferenceAnchorBatchSampler(Sampler[list[int]]):
    """
    Batch sampler that forces each training microbatch to contain a small
    donor-balanced reference anchor when possible.

    Intended purpose
    ----------------
    Make donor-balanced reference-centering active inside local DDP microbatches.

    Each anchored batch tries to include:

        same cell type c
        + reference cells only
        + ref_donors_per_batch distinct donors
        + ref_cells_per_donor cells from each donor

    The rest of the batch is filled with ordinary random cells.

    Notes
    -----
    - The sampler yields dataset-local indices, not global zarr row indices.
    - It samples with replacement. This is deliberate: its purpose is not exact
      one-pass epoch coverage, but stable activation of the v3 ref-center loss.
    - It supports DDP by giving every rank the same number of batches per epoch.
    """

    def __init__(
        self,
        dataset,
        *,
        batch_size: int,
        num_replicas: int = 1,
        rank: int = 0,
        seed: int = 42,
        drop_last: bool = False,
        anchor_probability: float = 1.0,
        ref_cells_per_donor: int = 1,
        ref_donors_per_batch: int = 2,
        fill_from_reference: bool = False,
    ) -> None:
        if batch_size <= 0:
            raise ValueError("batch_size must be positive.")
        if num_replicas <= 0:
            raise ValueError("num_replicas must be positive.")
        if not (0 <= rank < num_replicas):
            raise ValueError(f"rank must be in [0, {num_replicas}), got {rank}.")
        if not (0.0 <= anchor_probability <= 1.0):
            raise ValueError("anchor_probability must be in [0, 1].")
        if ref_cells_per_donor <= 0:
            raise ValueError("ref_cells_per_donor must be positive.")
        if ref_donors_per_batch <= 0:
            raise ValueError("ref_donors_per_batch must be positive.")

        anchor_size = int(ref_cells_per_donor) * int(ref_donors_per_batch)
        if anchor_size > batch_size:
            raise ValueError(
                "Reference anchor is larger than batch_size: "
                f"ref_cells_per_donor * ref_donors_per_batch = {anchor_size}, "
                f"batch_size = {batch_size}."
            )

        required_attrs = ["row_idx", "celltype_ids", "is_reference_origin", "donor_ids"]
        missing = [a for a in required_attrs if not hasattr(dataset, a)]
        if missing:
            raise AttributeError(
                "ReferenceAnchorBatchSampler requires dataset attributes: "
                f"{required_attrs}. Missing: {missing}"
            )

        self.dataset = dataset
        self.batch_size = int(batch_size)
        self.num_replicas = int(num_replicas)
        self.rank = int(rank)
        self.seed = int(seed)
        self.drop_last = bool(drop_last)
        self.anchor_probability = float(anchor_probability)
        self.ref_cells_per_donor = int(ref_cells_per_donor)
        self.ref_donors_per_batch = int(ref_donors_per_batch)
        self.fill_from_reference = bool(fill_from_reference)
        self.epoch = 0

        self.n = len(dataset)
        self.all_local_indices = np.arange(self.n, dtype=np.int64)

        # DDP-like epoch length.
        local_n = self.n / float(self.num_replicas)
        if self.drop_last:
            self.num_batches = max(1, int(math.floor(local_n / self.batch_size)))
        else:
            self.num_batches = max(1, int(math.ceil(local_n / self.batch_size)))

        self._build_reference_index()

    def _build_reference_index(self) -> None:
        """
        Build:
            ref_by_celltype_donor[celltype][donor] = [dataset-local indices]
        """
        row_idx = np.asarray(self.dataset.row_idx, dtype=np.int64)

        celltype_ids = np.asarray(self.dataset.celltype_ids, dtype=np.int64)
        is_reference = np.asarray(self.dataset.is_reference_origin, dtype=bool)
        donor_ids = np.asarray(self.dataset.donor_ids, dtype=np.int64)

        ref_by_ct_donor: dict[int, dict[int, list[int]]] = defaultdict(
            lambda: defaultdict(list)
        )

        all_ref_local = []

        for local_i, row in enumerate(row_idx):
            r = int(row)
            if not bool(is_reference[r]):
                continue

            ct = int(celltype_ids[r])
            donor = int(donor_ids[r])

            ref_by_ct_donor[ct][donor].append(int(local_i))
            all_ref_local.append(int(local_i))

        valid_celltypes = []
        for ct, donor_map in ref_by_ct_donor.items():
            valid_donors = [
                d
                for d, idxs in donor_map.items()
                if len(idxs) >= self.ref_cells_per_donor
            ]
            if len(valid_donors) >= self.ref_donors_per_batch:
                valid_celltypes.append(int(ct))

        self.ref_by_ct_donor = ref_by_ct_donor
        self.valid_celltypes = np.asarray(sorted(valid_celltypes), dtype=np.int64)
        self.all_ref_local_indices = np.asarray(all_ref_local, dtype=np.int64)

        if self.fill_from_reference and len(self.all_ref_local_indices) > 0:
            self.fill_pool = self.all_ref_local_indices
        else:
            self.fill_pool = self.all_local_indices

    def set_epoch(self, epoch: int) -> None:
        self.epoch = int(epoch)

    def __len__(self) -> int:
        return int(self.num_batches)

    def _sample_reference_anchor(self, rng: np.random.Generator) -> list[int]:
        if len(self.valid_celltypes) == 0:
            return []

        ct = int(rng.choice(self.valid_celltypes))
        donor_map = self.ref_by_ct_donor[ct]

        valid_donors = [
            d
            for d, idxs in donor_map.items()
            if len(idxs) >= self.ref_cells_per_donor
        ]

        if len(valid_donors) < self.ref_donors_per_batch:
            return []

        chosen_donors = rng.choice(
            np.asarray(valid_donors, dtype=np.int64),
            size=self.ref_donors_per_batch,
            replace=False,
        )

        anchor: list[int] = []
        for d in chosen_donors.tolist():
            pool = np.asarray(donor_map[int(d)], dtype=np.int64)
            chosen_cells = rng.choice(
                pool,
                size=self.ref_cells_per_donor,
                replace=(len(pool) < self.ref_cells_per_donor),
            )
            anchor.extend(int(x) for x in chosen_cells.tolist())

        return anchor

    def __iter__(self) -> Iterator[list[int]]:
        rng = np.random.default_rng(
            self.seed + 100_003 * int(self.epoch) + 9_973 * int(self.rank)
        )

        for _ in range(self.num_batches):
            batch: list[int] = []

            use_anchor = (
                self.anchor_probability > 0.0
                and rng.random() < self.anchor_probability
            )

            if use_anchor:
                batch.extend(self._sample_reference_anchor(rng))

            n_fill = self.batch_size - len(batch)
            if n_fill > 0:
                fill = rng.choice(
                    self.fill_pool,
                    size=n_fill,
                    replace=(len(self.fill_pool) < n_fill),
                )
                batch.extend(int(x) for x in fill.tolist())

            if len(batch) < self.batch_size:
                extra = rng.choice(
                    self.all_local_indices,
                    size=self.batch_size - len(batch),
                    replace=True,
                )
                batch.extend(int(x) for x in extra.tolist())

            rng.shuffle(batch)
            yield batch[: self.batch_size]