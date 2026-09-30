"""Cross-fitted donor-context data for PRISM Stage-2 shared pathology.

The prediction target removes measured nuisance effects but deliberately keeps
pathology-associated signal.  The personal encoder sees a stricter residual
from which a train-only fixed pathology ridge prediction has also been
removed.  Thus pathology can reach the learned shared branch explicitly,
while the personal branch receives a disease-firewalled source summary.

Every inner-training donor is transformed out-of-fold.  Inner-validation
donors are transformed by a model fitted on all inner-training donors.  The
actual target context's cell count and technology fraction are used only to
define the nuisance-adjusted outcome; they are never returned as model input.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import torch
from torch.utils.data import Dataset

from kmlee_bam.analysis.conditional_pathology_axes import balanced_donor_folds
from kmlee_bam.data.donor_context_dataset import (
    MASK_MODES,
    DonorContextTable,
    _fit_ridge,
    _safe_scale,
    pathology_features,
)


@dataclass
class PreparedSharedPathologyTable:
    raw: DonorContextTable
    personal_source_module: np.ndarray
    personal_source_latent: np.ndarray
    target_module: np.ndarray
    target_latent: np.ndarray
    ridge_disease_module: np.ndarray
    ridge_disease_latent: np.ndarray
    inner_train_donors: np.ndarray
    inner_validation_donors: np.ndarray
    crossfit_fold_assignment: np.ndarray
    nuisance_ridge: float
    pathology_ridge: float


def _nuisance_design(
    table: DonorContextTable,
    context: int,
    fit: np.ndarray,
    eps: float,
) -> np.ndarray:
    age_mean = float(np.nanmean(table.age[fit]))
    age_scale = float(np.nanstd(table.age[fit]))
    if not np.isfinite(age_scale) or age_scale < eps:
        age_scale = 1.0
    age = np.nan_to_num((table.age - age_mean) / age_scale)
    sex_mean = float(np.nanmean(table.sex[fit]))
    sex = np.nan_to_num(table.sex, nan=sex_mean)
    tech = table.tech_fraction[:, context]
    tech_one = tech[:, 1] if tech.shape[1] > 1 else tech[:, 0]
    log_count = np.log1p(table.cell_count[:, context].astype(np.float64))
    count_mean = float(np.mean(log_count[fit]))
    count_scale = float(np.std(log_count[fit]))
    if not np.isfinite(count_scale) or count_scale < eps:
        count_scale = 1.0
    return np.column_stack(
        [
            age,
            sex,
            tech_one,
            (log_count - count_mean) / count_scale,
            table.pathology_missing[:, 3].astype(np.float64),
        ]
    ).astype(np.float64)


def _transform_values(
    values: np.ndarray,
    nuisance: np.ndarray,
    pathology: np.ndarray,
    fit: np.ndarray,
    apply: np.ndarray,
    *,
    nuisance_ridge: float,
    pathology_ridge: float,
    eps: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = values.astype(np.float64)
    nuisance_beta, nuisance_mean, outcome_mean, _ = _fit_ridge(
        nuisance[fit], y[fit], nuisance_ridge
    )
    nuisance_prediction = (nuisance - nuisance_mean) @ nuisance_beta + outcome_mean
    nuisance_residual = y - nuisance_prediction
    target_scale = _safe_scale(nuisance_residual[fit], eps)
    target = nuisance_residual / target_scale

    pathology_beta, pathology_mean, pathology_outcome_mean, _ = _fit_ridge(
        pathology[fit], target[fit], pathology_ridge
    )
    disease_prediction = (
        (pathology - pathology_mean) @ pathology_beta + pathology_outcome_mean
    )
    personal_residual = target - disease_prediction
    personal_scale = _safe_scale(personal_residual[fit], eps)
    personal_source = personal_residual / personal_scale
    return (
        target[apply].astype(np.float32),
        personal_source[apply].astype(np.float32),
        disease_prediction[apply].astype(np.float32),
    )


def prepare_shared_pathology_crossfit(
    table: DonorContextTable,
    inner_train_donors: Sequence[int],
    inner_validation_donors: Sequence[int],
    *,
    nuisance_ridge: float = 100.0,
    pathology_ridge: float = 100.0,
    n_crossfit_folds: int = 5,
    fold_candidates: int = 2000,
    seed: int = 20260716,
    eps: float = 1e-5,
) -> PreparedSharedPathologyTable:
    """Create out-of-fold training residuals and held-out validation residuals."""
    inner_train = np.asarray(inner_train_donors, dtype=np.int64)
    inner_validation = np.asarray(inner_validation_donors, dtype=np.int64)
    folds = balanced_donor_folds(
        table.pathology,
        table.age,
        table.sex,
        inner_train,
        n_folds=n_crossfit_folds,
        candidates=fold_candidates,
        seed=seed,
    )
    shapes = {
        "module": table.module.shape,
        "latent": table.latent.shape,
    }
    target_module = np.zeros(shapes["module"], dtype=np.float32)
    target_latent = np.zeros(shapes["latent"], dtype=np.float32)
    personal_module = np.zeros(shapes["module"], dtype=np.float32)
    personal_latent = np.zeros(shapes["latent"], dtype=np.float32)
    ridge_module = np.zeros(shapes["module"], dtype=np.float32)
    ridge_latent = np.zeros(shapes["latent"], dtype=np.float32)
    pathology = pathology_features(table.pathology)

    transformations: list[tuple[np.ndarray, np.ndarray]] = []
    for fold in range(n_crossfit_folds):
        transformations.append(
            (inner_train[folds[inner_train] != fold], inner_train[folds[inner_train] == fold])
        )
    transformations.append((inner_train, inner_validation))

    for fit_donors, apply_donors in transformations:
        fit_donor_mask = np.zeros(len(table.donor_names), dtype=bool)
        fit_donor_mask[fit_donors] = True
        apply_donor_mask = np.zeros(len(table.donor_names), dtype=bool)
        apply_donor_mask[apply_donors] = True
        for context in range(len(table.context_names)):
            fit = fit_donor_mask & table.observed[:, context]
            apply = apply_donor_mask & table.observed[:, context]
            if not np.any(apply):
                continue
            if int(fit.sum()) < 12:
                raise ValueError(
                    f"context {table.context_names[context]} has only {int(fit.sum())} fit donors"
                )
            nuisance = _nuisance_design(table, context, fit, eps)
            for raw_values, target_out, personal_out, ridge_out in (
                (table.module[:, context], target_module, personal_module, ridge_module),
                (table.latent[:, context], target_latent, personal_latent, ridge_latent),
            ):
                target, personal, fixed = _transform_values(
                    raw_values,
                    nuisance,
                    pathology,
                    fit,
                    apply,
                    nuisance_ridge=nuisance_ridge,
                    pathology_ridge=pathology_ridge,
                    eps=eps,
                )
                rows = np.flatnonzero(apply)
                target_out[rows, context] = target
                personal_out[rows, context] = personal
                ridge_out[rows, context] = fixed

    retained = np.concatenate([inner_train, inner_validation])
    for name, values in (
        ("target_module", target_module),
        ("target_latent", target_latent),
        ("personal_source_module", personal_module),
        ("personal_source_latent", personal_latent),
        ("ridge_disease_module", ridge_module),
        ("ridge_disease_latent", ridge_latent),
    ):
        observed = table.observed[retained]
        expanded = observed[..., None]
        if not np.isfinite(values[retained][expanded.repeat(values.shape[-1], axis=-1)]).all():
            raise ValueError(f"non-finite values in {name}")

    return PreparedSharedPathologyTable(
        raw=table,
        personal_source_module=personal_module,
        personal_source_latent=personal_latent,
        target_module=target_module,
        target_latent=target_latent,
        ridge_disease_module=ridge_module,
        ridge_disease_latent=ridge_latent,
        inner_train_donors=inner_train,
        inner_validation_donors=inner_validation,
        crossfit_fold_assignment=folds,
        nuisance_ridge=float(nuisance_ridge),
        pathology_ridge=float(pathology_ridge),
    )


def joint_pathology_permutation(
    table: DonorContextTable,
    donor_groups: Sequence[Sequence[int]],
    *,
    seed: int,
    window_size: int = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Jointly shuffle all axes within sex/nearby-age windows.

    Keeping each five-axis row intact preserves Thal/Braak/CERAD/LATE/Lewy
    correlation.  Separate shuffling inside each supplied split prevents
    information from crossing train/validation boundaries.
    """
    pathology = table.pathology.copy()
    missing = table.pathology_missing.copy()
    source_donor = np.arange(len(table.donor_names), dtype=np.int64)
    rng = np.random.default_rng(seed)
    for raw_group in donor_groups:
        group = np.asarray(raw_group, dtype=np.int64)
        for sex_value in np.unique(table.sex[group][np.isfinite(table.sex[group])]):
            same_sex = group[np.isclose(table.sex[group], sex_value)]
            order = same_sex[np.argsort(table.age[same_sex], kind="stable")]
            for start in range(0, len(order), window_size):
                window = order[start : start + window_size]
                if len(window) < 2:
                    # Attach a singleton to the preceding window when possible.
                    previous = order[max(0, start - window_size) : start + len(window)]
                    window = previous if len(previous) > 1 else window
                if len(window) < 2:
                    continue
                shift = int(rng.integers(1, len(window)))
                sources = np.roll(window, shift)
                pathology[window] = table.pathology[sources]
                missing[window] = table.pathology_missing[sources]
                source_donor[window] = sources
    return pathology, missing, source_donor


class SharedPathologyMaskedDataset(Dataset[dict[str, torch.Tensor]]):
    """Masked context examples for the Stage-2 A/B/R/P/E/S arms."""

    def __init__(
        self,
        table: PreparedSharedPathologyTable,
        donors: Iterable[int],
        *,
        mask_modes: Sequence[str] = MASK_MODES,
        reliability_cap: int = 100,
        min_source_contexts: int = 2,
        pathology_override: np.ndarray | None = None,
        pathology_missing_override: np.ndarray | None = None,
    ) -> None:
        self.table = table
        self.donors = np.asarray(list(donors), dtype=np.int64)
        self.mask_modes = tuple(str(mode) for mode in mask_modes)
        invalid = sorted(set(self.mask_modes).difference({"context", "celltype", "region"}))
        if invalid:
            raise ValueError(f"unknown mask modes: {invalid}")
        self.reliability_cap = int(reliability_cap)
        self.min_source_contexts = int(min_source_contexts)
        raw = table.raw
        self.pathology = raw.pathology if pathology_override is None else pathology_override
        self.pathology_missing = (
            raw.pathology_missing
            if pathology_missing_override is None
            else pathology_missing_override
        )
        self.samples: list[tuple[int, int, str]] = []
        for donor in self.donors:
            for target in np.flatnonzero(raw.observed[int(donor)]):
                for mode in self.mask_modes:
                    source = self.source_mask(int(donor), int(target), mode)
                    if int(source.sum()) >= self.min_source_contexts:
                        self.samples.append((int(donor), int(target), mode))
        if not self.samples:
            raise ValueError("masked dataset has no usable samples")
        per_donor = {int(donor): 0 for donor in self.donors}
        for donor, _, _ in self.samples:
            per_donor[donor] += 1
        self.sample_weight = np.asarray(
            [1.0 / per_donor[donor] for donor, _, _ in self.samples], dtype=np.float32
        )

    def source_mask(self, donor: int, target: int, mode: str) -> np.ndarray:
        raw = self.table.raw
        source = raw.observed[donor].copy()
        if mode == "context":
            source[target] = False
        elif mode == "celltype":
            source[raw.context_celltype == raw.context_celltype[target]] = False
        elif mode == "region":
            source[raw.context_region == raw.context_region[target]] = False
        else:
            raise ValueError(f"unknown mask mode: {mode}")
        return source

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        donor, target, mode = self.samples[index]
        table = self.table
        raw = table.raw
        source = self.source_mask(donor, target, mode)
        reliability = np.sqrt(
            np.minimum(raw.cell_count[donor], self.reliability_cap)
            / float(self.reliability_cap)
        ).astype(np.float32)
        reliability *= source.astype(np.float32)
        mode_id = {"context": 0, "celltype": 1, "region": 2}[mode]
        return {
            "source_module": torch.from_numpy(table.personal_source_module[donor]),
            "source_latent": torch.from_numpy(table.personal_source_latent[donor]),
            "source_mask": torch.from_numpy(source),
            "source_reliability": torch.from_numpy(reliability),
            "target_module": torch.from_numpy(table.target_module[donor, target]),
            "target_latent": torch.from_numpy(table.target_latent[donor, target]),
            "personal_target_module": torch.from_numpy(
                table.target_module[donor, target]
                - table.ridge_disease_module[donor, target]
            ),
            "personal_target_latent": torch.from_numpy(
                table.target_latent[donor, target]
                - table.ridge_disease_latent[donor, target]
            ),
            "ridge_disease_module": torch.from_numpy(
                table.ridge_disease_module[donor, target]
            ),
            "ridge_disease_latent": torch.from_numpy(
                table.ridge_disease_latent[donor, target]
            ),
            "target_context": torch.tensor(target, dtype=torch.long),
            "target_region": torch.tensor(raw.context_region[target], dtype=torch.long),
            "target_celltype": torch.tensor(raw.context_celltype[target], dtype=torch.long),
            "pathology": torch.from_numpy(self.pathology[donor]),
            "pathology_missing": torch.from_numpy(self.pathology_missing[donor]),
            "donor": torch.tensor(donor, dtype=torch.long),
            "mask_mode": torch.tensor(mode_id, dtype=torch.long),
            "sample_weight": torch.tensor(self.sample_weight[index], dtype=torch.float32),
        }
