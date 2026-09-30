"""Leakage-safe donor split used by PRISM architecture search.

The official validation and test donors are never touched here.  This module
splits the donors already assigned to the training split into:

* ``weight_train`` donors used for ordinary model/module-rescue updates; and
* ``architecture_train`` donors used only for generator-gate optimisation.

The split is deterministic and balances donor-level pathology, missingness,
sex, age, region support, cell-type support, and cell count.  ``DonorIndexView``
is deliberately a transparent dataset proxy rather than ``torch.Subset``:
legacy trainers inspect attributes such as ``row_idx``, ``donor_ids``, and
``spec`` on ``loader.dataset``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any, Iterable

import numpy as np
import torch
from torch.utils.data import Dataset


@dataclass(frozen=True)
class ArchitectureDonorSplit:
    """One immutable train-only donor partition and its audit manifest."""

    weight_donor_ids: tuple[int, ...]
    architecture_donor_ids: tuple[int, ...]
    manifest: dict[str, Any]


def _string_sequence_sha256(values: Iterable[str]) -> str:
    payload = json.dumps(
        [str(value) for value in values],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _donor_constant(
    values: np.ndarray,
    rows: np.ndarray,
    *,
    valid: np.ndarray | None = None,
) -> float:
    selected = np.asarray(values)[rows]
    if valid is not None:
        selected = selected[np.asarray(valid)[rows].astype(bool)]
    if selected.ndim > 1:
        raise ValueError("_donor_constant expects one scalar column")
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return float("nan")
    if float(np.ptp(selected.astype(np.float64))) > 1e-6:
        raise ValueError("a donor-level covariate varies within donor")
    return float(selected[0])


def _build_donor_feature_matrix(
    dataset: Any,
    donor_ids: np.ndarray,
) -> tuple[np.ndarray, tuple[str, ...]]:
    row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
    global_donor = np.asarray(dataset.donor_ids, dtype=np.int64)
    global_celltype = np.asarray(dataset.celltype_ids, dtype=np.int64)
    global_region = np.asarray(dataset.region_ids, dtype=np.int64)
    pathology = np.asarray(dataset.prism_pathology, dtype=np.float64)
    pathology_valid = np.asarray(
        dataset.prism_pathology_valid, dtype=bool
    )
    n_celltype = len(dataset.spec.celltype_vocab)
    n_region = int(getattr(dataset, "n_region", int(global_region.max()) + 1))

    split_mask = np.zeros(global_donor.shape[0], dtype=bool)
    split_mask[row_idx] = True
    feature_rows: list[list[float]] = []
    names = (
        tuple(f"pathology_{axis}" for axis in range(pathology.shape[1]))
        + tuple(f"pathology_valid_{axis}" for axis in range(pathology.shape[1]))
        + ("sex", "age", "log_cell_count")
        + tuple(f"region_fraction_{region}" for region in range(n_region))
        + tuple(f"celltype_present_{celltype}" for celltype in range(n_celltype))
    )

    sex_source = (
        dataset.sex_ids
        if hasattr(dataset, "sex_ids")
        else np.full(global_donor.shape, -1.0, dtype=np.float64)
    )
    age_source = (
        dataset.age_years
        if hasattr(dataset, "age_years")
        else np.full(global_donor.shape, np.nan, dtype=np.float64)
    )
    sex_values = np.asarray(sex_source, dtype=np.float64)
    age_values = np.asarray(age_source, dtype=np.float64)
    age_valid = np.asarray(
        getattr(dataset, "age_valid", np.isfinite(age_values)),
        dtype=bool,
    )

    for donor in donor_ids.tolist():
        rows = np.flatnonzero(split_mask & (global_donor == int(donor)))
        if rows.size == 0:
            raise ValueError(f"train donor {donor} has no selected cells")

        p_values: list[float] = []
        p_valid: list[float] = []
        for axis in range(pathology.shape[1]):
            valid = pathology_valid[rows, axis]
            p_valid.append(float(bool(valid.any())))
            finite = pathology[rows[valid], axis]
            finite = finite[np.isfinite(finite)]
            p_values.append(
                float(np.median(finite)) if finite.size else float("nan")
            )

        sex = _donor_constant(sex_values, rows)
        age = _donor_constant(age_values, rows, valid=age_valid)
        region_count = np.bincount(
            global_region[rows], minlength=n_region
        ).astype(np.float64)
        region_fraction = region_count / max(float(region_count.sum()), 1.0)
        celltype_present = (
            np.bincount(global_celltype[rows], minlength=n_celltype) > 0
        ).astype(np.float64)
        feature_rows.append(
            p_values
            + p_valid
            + [sex, age, float(np.log1p(rows.size))]
            + region_fraction.tolist()
            + celltype_present.tolist()
        )

    matrix = np.asarray(feature_rows, dtype=np.float64)
    if matrix.shape != (len(donor_ids), len(names)):
        raise RuntimeError(
            f"donor feature shape mismatch: {matrix.shape} vs "
            f"{(len(donor_ids), len(names))}"
        )
    median = np.nanmedian(matrix, axis=0)
    median = np.where(np.isfinite(median), median, 0.0)
    matrix = np.where(np.isfinite(matrix), matrix, median[None, :])
    scale = np.std(matrix, axis=0)
    keep = scale > 1e-8
    if not bool(keep.any()):
        raise ValueError("all donor-balancing features are constant")
    matrix = (matrix[:, keep] - matrix[:, keep].mean(axis=0)) / scale[keep]
    kept_names = tuple(name for name, flag in zip(names, keep) if bool(flag))
    return matrix.astype(np.float64), kept_names


def _partition_score(
    features: np.ndarray,
    architecture_index: np.ndarray,
) -> float:
    """Balance both means and variances against the complete train cohort."""

    select = np.zeros(features.shape[0], dtype=bool)
    select[architecture_index] = True
    architecture = features[select]
    weight = features[~select]
    full_mean = features.mean(axis=0)
    full_var = features.var(axis=0)
    mean_error = np.mean(
        np.square(architecture.mean(axis=0) - full_mean)
        + np.square(weight.mean(axis=0) - full_mean)
    )
    variance_error = np.mean(
        np.square(architecture.var(axis=0) - full_var)
        + np.square(weight.var(axis=0) - full_var)
    )
    return float(mean_error + 0.20 * variance_error)


def make_architecture_donor_split(
    dataset: Any,
    *,
    architecture_donor_count: int = 12,
    seed: int = 420730,
    search_trials: int = 20_000,
) -> ArchitectureDonorSplit:
    """Return a deterministic balanced partition of the dataset's train donors."""

    if int(architecture_donor_count) <= 0:
        raise ValueError("architecture_donor_count must be positive")
    if int(search_trials) <= 0:
        raise ValueError("search_trials must be positive")
    row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
    donor_ids = np.unique(np.asarray(dataset.donor_ids)[row_idx]).astype(
        np.int64
    )
    n_donor = int(donor_ids.size)
    n_arch = int(architecture_donor_count)
    if n_arch >= n_donor:
        raise ValueError(
            "architecture_donor_count must be smaller than train donor count"
        )
    features, feature_names = _build_donor_feature_matrix(dataset, donor_ids)

    rng = np.random.default_rng(int(seed))
    best_index: np.ndarray | None = None
    best_score = float("inf")
    for _ in range(int(search_trials)):
        candidate = np.sort(
            rng.choice(n_donor, size=n_arch, replace=False).astype(np.int64)
        )
        score = _partition_score(features, candidate)
        if score < best_score:
            best_score = score
            best_index = candidate
    if best_index is None:
        raise RuntimeError("architecture donor search produced no candidate")

    # Deterministic one-swap local improvement.
    improved = True
    while improved:
        improved = False
        selected = set(int(value) for value in best_index.tolist())
        for outgoing in sorted(selected):
            for incoming in range(n_donor):
                if incoming in selected:
                    continue
                candidate_set = (selected - {outgoing}) | {incoming}
                candidate = np.asarray(sorted(candidate_set), dtype=np.int64)
                score = _partition_score(features, candidate)
                if score + 1e-12 < best_score:
                    best_score = score
                    best_index = candidate
                    improved = True
                    break
            if improved:
                break

    architecture_ids = tuple(
        int(value) for value in donor_ids[best_index].tolist()
    )
    architecture_set = set(architecture_ids)
    weight_ids = tuple(
        int(value)
        for value in donor_ids.tolist()
        if int(value) not in architecture_set
    )
    donor_names = np.asarray(dataset.donor_vocab, dtype=object).astype(str)
    architecture_names = tuple(donor_names[list(architecture_ids)].tolist())
    weight_names = tuple(donor_names[list(weight_ids)].tolist())
    manifest: dict[str, Any] = {
        "schema_version": "kmlee_bam.architecture_donor_split.v1",
        "source_split": "train_only",
        "seed": int(seed),
        "search_trials": int(search_trials),
        "balance_score": float(best_score),
        "feature_names": list(feature_names),
        "train_donor_count": n_donor,
        "weight_donor_count": len(weight_ids),
        "architecture_donor_count": len(architecture_ids),
        "weight_donor_ids": list(weight_ids),
        "architecture_donor_ids": list(architecture_ids),
        "weight_donor_names": list(weight_names),
        "architecture_donor_names": list(architecture_names),
        "weight_donor_sha256": _string_sequence_sha256(weight_names),
        "architecture_donor_sha256": _string_sequence_sha256(
            architecture_names
        ),
        "official_validation_used": False,
        "official_test_used": False,
    }
    return ArchitectureDonorSplit(
        weight_donor_ids=weight_ids,
        architecture_donor_ids=architecture_ids,
        manifest=manifest,
    )


class DonorIndexView(Dataset):
    """Dataset proxy restricted to a fixed donor set.

    The view keeps the original global donor vocabulary/ids.  Only ``row_idx``
    and per-cell donor-balance weights change.
    """

    def __init__(self, dataset: Dataset, donor_ids: Iterable[int]) -> None:
        super().__init__()
        self.base_dataset = dataset
        selected_donors = tuple(sorted({int(value) for value in donor_ids}))
        if not selected_donors:
            raise ValueError("DonorIndexView requires at least one donor")
        self.selected_donor_ids = selected_donors

        base_row = np.asarray(getattr(dataset, "row_idx"), dtype=np.int64)
        donor_by_row = np.asarray(getattr(dataset, "donor_ids"), dtype=np.int64)
        keep = np.isin(
            donor_by_row[base_row],
            np.asarray(selected_donors, dtype=np.int64),
        )
        self.local_indices = np.flatnonzero(keep).astype(np.int64)
        self.row_idx = base_row[self.local_indices]
        if self.row_idx.size == 0:
            raise ValueError("DonorIndexView selected no cells")

        selected_cell_donor = donor_by_row[self.row_idx]
        count = np.bincount(
            selected_cell_donor,
            minlength=int(getattr(dataset, "n_donor")),
        ).astype(np.float64)
        present = count > 0
        per_donor = np.zeros_like(count, dtype=np.float64)
        per_donor[present] = len(self.row_idx) / (
            int(present.sum()) * count[present]
        )
        self._balance_weight = per_donor[selected_cell_donor].astype(
            np.float32
        )

    def __len__(self) -> int:
        return int(self.local_indices.size)

    def __getitem__(self, index: int) -> dict[str, Any]:
        source_index = int(self.local_indices[int(index)])
        item = self.base_dataset[source_index]
        item["donor_balance_weight"] = torch.tensor(
            float(self._balance_weight[int(index)]),
            dtype=torch.float32,
        )
        return item

    def __getattr__(self, name: str) -> Any:
        if name in {
            "base_dataset",
            "selected_donor_ids",
            "local_indices",
            "row_idx",
            "_balance_weight",
        }:
            raise AttributeError(name)
        return getattr(self.base_dataset, name)
