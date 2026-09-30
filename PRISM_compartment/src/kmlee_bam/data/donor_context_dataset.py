"""Leakage-resistant donor-context tables for PRISM-BAM Stage 1.

The source artifact contains one frozen summary per
``donor x cell type x region``.  This module turns those rows into a dense
donor-context table, fits *training-donor-only* multi-pathology/nuisance
baselines, and constructs masked prediction examples.

ADNC is intentionally report-only.  The conditioning coordinates are the
named pathologies Thal, Braak, CERAD, LATE-NC, and Lewy pathology.  LATE and
Lewy are co-pathologies, not components of ADNC.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import torch
from torch.utils.data import Dataset


PATHOLOGY_NAMES = ("thal", "braak", "cerad", "late", "lewy")
PATHOLOGY_LABELS = (
    "amyloid_Thal",
    "tau_Braak",
    "neuritic_plaque_CERAD",
    "LATE_NC_TDP43",
    "Lewy_body",
)
PATHOLOGY_SCALES = np.asarray([5.0, 6.0, 3.0, 3.0, 1.0], dtype=np.float64)
MASK_MODES = ("celltype", "region")


def pathology_features(pathology: np.ndarray) -> np.ndarray:
    """Fixed low-degree monotone hinge features, three per named axis."""
    p = np.asarray(pathology, dtype=np.float64)
    if p.shape[-1] != len(PATHOLOGY_NAMES):
        raise ValueError("pathology must end in five named axes")
    return np.concatenate(
        [p, np.maximum(p - 1.0 / 3.0, 0.0), np.maximum(p - 2.0 / 3.0, 0.0)],
        axis=-1,
    )


def _fit_ridge(
    x: np.ndarray,
    y: np.ndarray,
    ridge: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit a centred multivariate ridge with an unpenalized intercept."""
    x_mean = np.mean(x, axis=0)
    y_mean = np.mean(y, axis=0)
    xc = x - x_mean
    yc = y - y_mean
    gram = xc.T @ xc + float(ridge) * np.eye(x.shape[1], dtype=np.float64)
    beta = np.linalg.solve(gram, xc.T @ yc)
    return beta, x_mean, y_mean, (x - x_mean) @ beta + y_mean


def _safe_scale(values: np.ndarray, eps: float) -> np.ndarray:
    scale = np.std(values, axis=0)
    return np.where(np.isfinite(scale) & (scale > eps), scale, 1.0)


@dataclass
class DonorContextTable:
    module: np.ndarray
    latent: np.ndarray
    observed: np.ndarray
    cell_count: np.ndarray
    tech_fraction: np.ndarray
    donor_names: np.ndarray
    context_names: np.ndarray
    celltype_names: np.ndarray
    context_celltype: np.ndarray
    context_region: np.ndarray
    region_names: np.ndarray
    donor_split: np.ndarray
    pathology: np.ndarray
    pathology_missing: np.ndarray
    adnc_report: np.ndarray
    age: np.ndarray
    sex: np.ndarray
    module_names: np.ndarray
    source_path: str
    source_cells: int

    @classmethod
    def from_npz(cls, path: str | Path, min_cells: int = 10) -> "DonorContextTable":
        source = Path(path)
        archive = np.load(source, allow_pickle=True)
        required = {
            "z", "module_activity", "donor", "celltype", "region", "split",
            "age", "sex", "cell_count", "tech_fraction", "celltype_vocab",
            "donor_vocab", "module_names", "thal", "braak", "cerad", "late",
            "adnc", "lewy",
        }
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(f"missing arrays in {source}: {missing}")

        # ``NpzFile`` decompresses lazily on every ``archive[key]`` access.
        # Cache each array once; row-wise access to the archive would otherwise
        # decompress the 4108 x 414 module matrix thousands of times.
        data = {key: archive[key] for key in required}
        n_source_cells = int(
            np.asarray(archive["n_source_cells"]).item()
            if "n_source_cells" in archive.files else 0
        )
        support_donor_scope = str(
            np.asarray(archive["support_donor_scope"]).item()
            if "support_donor_scope" in archive.files
            else "all_available_splits"
        )
        allow_unavailable_donors = support_donor_scope == "train64"
        archive.close()

        donor_names = np.asarray(data["donor_vocab"], dtype=object).astype(str)
        celltype_names = np.asarray(data["celltype_vocab"], dtype=object).astype(str)
        region_names = np.asarray(sorted({str(x) for x in data["region"]}), dtype=object)
        context_names = np.asarray(
            [f"{ct}|{region}" for ct in celltype_names for region in region_names],
            dtype=object,
        )
        context_lookup = {str(name): i for i, name in enumerate(context_names)}
        n_donors = len(donor_names)
        n_contexts = len(context_names)
        n_modules = int(data["module_activity"].shape[1])
        latent_dim = int(data["z"].shape[1])

        module = np.zeros((n_donors, n_contexts, n_modules), dtype=np.float32)
        latent = np.zeros((n_donors, n_contexts, latent_dim), dtype=np.float32)
        observed = np.zeros((n_donors, n_contexts), dtype=bool)
        cell_count = np.zeros((n_donors, n_contexts), dtype=np.int64)
        n_tech = int(data["tech_fraction"].shape[1])
        tech_fraction = np.zeros((n_donors, n_contexts, n_tech), dtype=np.float32)

        row_donor = np.asarray(data["donor"], dtype=np.int64)
        row_celltype = np.asarray(data["celltype"], dtype=np.int64)
        row_region = np.asarray(data["region"], dtype=object).astype(str)
        row_count = np.asarray(data["cell_count"], dtype=np.int64)
        for row, (donor, celltype, region, count) in enumerate(
            zip(row_donor, row_celltype, row_region, row_count)
        ):
            context = context_lookup[f"{celltype_names[int(celltype)]}|{region}"]
            if cell_count[int(donor), context] != 0:
                raise ValueError(f"duplicate donor/context row: {donor}/{context}")
            cell_count[int(donor), context] = int(count)
            tech_fraction[int(donor), context] = data["tech_fraction"][row]
            if int(count) < int(min_cells):
                continue
            values_module = np.asarray(data["module_activity"][row], dtype=np.float32)
            values_latent = np.asarray(data["z"][row], dtype=np.float32)
            if not np.isfinite(values_module).all() or not np.isfinite(values_latent).all():
                continue
            module[int(donor), context] = values_module
            latent[int(donor), context] = values_latent
            observed[int(donor), context] = True

        donor_split = np.empty(n_donors, dtype=object)
        age = np.full(n_donors, np.nan, dtype=np.float64)
        sex = np.full(n_donors, np.nan, dtype=np.float64)
        pathology = np.full((n_donors, len(PATHOLOGY_NAMES)), np.nan, dtype=np.float64)
        adnc_report = np.full(n_donors, np.nan, dtype=np.float64)
        for donor in range(n_donors):
            rows = np.flatnonzero(row_donor == donor)
            if rows.size == 0:
                if allow_unavailable_donors:
                    donor_split[donor] = "unavailable"
                    continue
                raise ValueError(f"donor {donor} has no rows")

            def constant(key: str) -> float:
                values = np.asarray(data[key][rows], dtype=np.float64)
                finite = values[np.isfinite(values)]
                if finite.size == 0:
                    return float("nan")
                if np.ptp(finite) > 1e-8:
                    raise ValueError(f"{key} varies within donor {donor}")
                return float(finite[0])

            splits = {str(value) for value in data["split"][rows]}
            if len(splits) != 1:
                raise ValueError(f"split varies within donor {donor}: {splits}")
            donor_split[donor] = splits.pop()
            age[donor] = constant("age")
            sex[donor] = constant("sex")
            pathology[donor] = [constant(name) for name in PATHOLOGY_NAMES]
            adnc_report[donor] = constant("adnc")

        pathology /= PATHOLOGY_SCALES[None, :]
        pathology_missing = ~np.isfinite(pathology)
        train = donor_split.astype(str) == "train"
        impute = np.nanmedian(pathology[train], axis=0)
        pathology = np.where(pathology_missing, impute[None, :], pathology)
        pathology = np.clip(pathology, 0.0, 1.0)

        region_lookup = {str(name): i for i, name in enumerate(region_names)}
        context_celltype = np.repeat(np.arange(len(celltype_names)), len(region_names))
        context_region = np.tile(
            np.asarray([region_lookup[str(x)] for x in region_names], dtype=np.int64),
            len(celltype_names),
        )
        return cls(
            module=module,
            latent=latent,
            observed=observed,
            cell_count=cell_count,
            tech_fraction=tech_fraction,
            donor_names=donor_names,
            context_names=context_names,
            celltype_names=celltype_names,
            context_celltype=context_celltype.astype(np.int64),
            context_region=context_region,
            region_names=region_names,
            donor_split=donor_split.astype(str),
            pathology=pathology.astype(np.float32),
            pathology_missing=pathology_missing,
            adnc_report=adnc_report.astype(np.float32),
            age=age.astype(np.float32),
            sex=sex.astype(np.float32),
            module_names=np.asarray(data["module_names"], dtype=object).astype(str),
            source_path=str(source),
            source_cells=n_source_cells,
        )

    def prepare(
        self,
        fit_donors: Sequence[int],
        ridge: float = 10.0,
        eps: float = 1e-5,
    ) -> "PreparedDonorContextTable":
        """Fit common named-pathology/nuisance effects on fit donors only."""
        fit_donors = np.asarray(fit_donors, dtype=np.int64)
        fit_mask = np.zeros(len(self.donor_names), dtype=bool)
        fit_mask[fit_donors] = True
        p_features = pathology_features(self.pathology)
        n_path_features = p_features.shape[1]

        age_mean = float(np.nanmean(self.age[fit_mask]))
        age_scale = float(np.nanstd(self.age[fit_mask]))
        if not np.isfinite(age_scale) or age_scale < eps:
            age_scale = 1.0
        age_z = np.nan_to_num((self.age - age_mean) / age_scale)

        module_residual = np.zeros_like(self.module, dtype=np.float32)
        latent_residual = np.zeros_like(self.latent, dtype=np.float32)
        module_scale = np.ones((len(self.context_names), self.module.shape[-1]), dtype=np.float32)
        latent_scale = np.ones((len(self.context_names), self.latent.shape[-1]), dtype=np.float32)
        common_module_coef = np.zeros(
            (len(self.context_names), n_path_features, self.module.shape[-1]), dtype=np.float32
        )
        common_latent_coef = np.zeros(
            (len(self.context_names), n_path_features, self.latent.shape[-1]), dtype=np.float32
        )
        axis_direction_module = np.zeros(
            (len(self.context_names), len(PATHOLOGY_NAMES), self.module.shape[-1]),
            dtype=np.float32,
        )
        axis_direction_latent = np.zeros(
            (len(self.context_names), len(PATHOLOGY_NAMES), self.latent.shape[-1]),
            dtype=np.float32,
        )
        n_fit_context = np.zeros(len(self.context_names), dtype=np.int64)

        for context in range(len(self.context_names)):
            observed = self.observed[:, context]
            fit = fit_mask & observed
            n_fit_context[context] = int(fit.sum())
            if int(fit.sum()) < 12:
                raise ValueError(
                    f"context {self.context_names[context]} has only {int(fit.sum())} fit donors"
                )
            tech = self.tech_fraction[:, context]
            tech_one = tech[:, 1] if tech.shape[1] > 1 else tech[:, 0]
            log_count = np.log1p(self.cell_count[:, context].astype(np.float64))
            log_mean = float(np.mean(log_count[fit]))
            log_scale = float(np.std(log_count[fit]))
            if log_scale < eps:
                log_scale = 1.0
            nuisance = np.column_stack(
                [
                    age_z,
                    np.nan_to_num(self.sex, nan=float(np.nanmean(self.sex[fit_mask]))),
                    tech_one,
                    (log_count - log_mean) / log_scale,
                    self.pathology_missing[:, 3].astype(np.float64),
                ]
            )
            design = np.column_stack([p_features, nuisance])

            for values, residual_out, scale_out, coef_out, direction_out in (
                (
                    self.module, module_residual, module_scale,
                    common_module_coef, axis_direction_module,
                ),
                (
                    self.latent, latent_residual, latent_scale,
                    common_latent_coef, axis_direction_latent,
                ),
            ):
                beta, x_mean, y_mean, _ = _fit_ridge(
                    design[fit], values[fit, context].astype(np.float64), ridge
                )
                prediction = (design - x_mean) @ beta + y_mean
                residual = values[:, context].astype(np.float64) - prediction
                scale = _safe_scale(residual[fit], eps)
                standardized = residual / scale
                standardized[~observed] = 0.0
                residual_out[:, context] = standardized.astype(np.float32)
                scale_out[context] = scale.astype(np.float32)
                coef_out[context] = beta[:n_path_features].astype(np.float32)

                reference = np.nanmedian(self.pathology[fit], axis=0)
                for axis in range(len(PATHOLOGY_NAMES)):
                    low = reference.copy()
                    high = reference.copy()
                    low[axis] = 0.0
                    high[axis] = 1.0
                    delta = pathology_features(high[None])[0] - pathology_features(low[None])[0]
                    direction_out[context, axis] = ((delta @ beta[:n_path_features]) / scale).astype(
                        np.float32
                    )

        return PreparedDonorContextTable(
            raw=self,
            module_residual=module_residual,
            latent_residual=latent_residual,
            module_scale=module_scale,
            latent_scale=latent_scale,
            common_module_coef=common_module_coef,
            common_latent_coef=common_latent_coef,
            axis_direction_module=axis_direction_module,
            axis_direction_latent=axis_direction_latent,
            n_fit_context=n_fit_context,
            fit_donors=fit_donors,
            ridge=float(ridge),
        )


@dataclass
class PreparedDonorContextTable:
    raw: DonorContextTable
    module_residual: np.ndarray
    latent_residual: np.ndarray
    module_scale: np.ndarray
    latent_scale: np.ndarray
    common_module_coef: np.ndarray
    common_latent_coef: np.ndarray
    axis_direction_module: np.ndarray
    axis_direction_latent: np.ndarray
    n_fit_context: np.ndarray
    fit_donors: np.ndarray
    ridge: float


class MaskedDonorContextDataset(Dataset[dict[str, torch.Tensor]]):
    """One target context with either its cell type or its region hidden."""

    def __init__(
        self,
        table: PreparedDonorContextTable,
        donors: Iterable[int],
        mask_modes: Sequence[str] = MASK_MODES,
        reliability_cap: int = 100,
        min_source_contexts: int = 2,
    ) -> None:
        self.table = table
        self.donors = np.asarray(list(donors), dtype=np.int64)
        self.mask_modes = tuple(str(mode) for mode in mask_modes)
        invalid = sorted(set(self.mask_modes).difference({"context", "celltype", "region"}))
        if invalid:
            raise ValueError(f"unknown mask modes: {invalid}")
        self.reliability_cap = int(reliability_cap)
        self.min_source_contexts = int(min_source_contexts)
        self.samples: list[tuple[int, int, str]] = []
        raw = table.raw
        for donor in self.donors:
            for target in np.flatnonzero(raw.observed[int(donor)]):
                for mode in self.mask_modes:
                    source = self.source_mask(int(donor), int(target), mode)
                    if int(source.sum()) >= self.min_source_contexts:
                        self.samples.append((int(donor), int(target), mode))
        if not self.samples:
            raise ValueError("masked dataset has no usable samples")
        per_donor = {int(d): 0 for d in self.donors}
        for donor, _, _ in self.samples:
            per_donor[donor] += 1
        self.sample_weight = np.asarray(
            [1.0 / per_donor[donor] for donor, _, _ in self.samples], dtype=np.float32
        )

    def source_mask(self, donor: int, target: int, mode: str) -> np.ndarray:
        raw = self.table.raw
        source = raw.observed[int(donor)].copy()
        if mode == "context":
            source[int(target)] = False
        elif mode == "celltype":
            source[raw.context_celltype == raw.context_celltype[int(target)]] = False
        elif mode == "region":
            source[raw.context_region == raw.context_region[int(target)]] = False
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
            "source_module": torch.from_numpy(table.module_residual[donor]),
            "source_latent": torch.from_numpy(table.latent_residual[donor]),
            "source_mask": torch.from_numpy(source),
            "source_reliability": torch.from_numpy(reliability),
            "target_module": torch.from_numpy(table.module_residual[donor, target]),
            "target_latent": torch.from_numpy(table.latent_residual[donor, target]),
            "target_context": torch.tensor(target, dtype=torch.long),
            "target_region": torch.tensor(raw.context_region[target], dtype=torch.long),
            "target_celltype": torch.tensor(raw.context_celltype[target], dtype=torch.long),
            "pathology": torch.from_numpy(raw.pathology[donor]),
            "pathology_missing": torch.from_numpy(raw.pathology_missing[donor]),
            "donor": torch.tensor(donor, dtype=torch.long),
            "mask_mode": torch.tensor(mode_id, dtype=torch.long),
            "sample_weight": torch.tensor(self.sample_weight[index], dtype=torch.float32),
        }
