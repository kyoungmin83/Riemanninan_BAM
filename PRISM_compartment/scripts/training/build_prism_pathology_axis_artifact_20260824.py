#!/usr/bin/env python3
"""Augment the frozen train64 module target with train-only pathology axes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


AXES = ("braak", "thal", "cerad", "late", "lewy")
SEED = 420824


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fold(name: str) -> int:
    payload = f"{SEED}:{name}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % 2


def _fit_predict_ridge(x_train, y_train, x_test, ridge: float = 1.0):
    mean = np.nanmedian(x_train, axis=0)
    x_train = np.where(np.isfinite(x_train), x_train, mean)
    x_test = np.where(np.isfinite(x_test), x_test, mean)
    scale = np.nanstd(x_train, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    x_train = (x_train - mean) / scale
    x_test = (x_test - mean) / scale
    design = np.column_stack((np.ones(x_train.shape[0]), x_train))
    test_design = np.column_stack((np.ones(x_test.shape[0]), x_test))
    penalty = np.eye(design.shape[1]) * float(ridge)
    penalty[0, 0] = 0.0
    coef = np.linalg.solve(design.T @ design + penalty, design.T @ y_train)
    return test_design @ coef


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    with np.load(args.base, allow_pickle=False) as archive:
        base = {key: np.asarray(archive[key]) for key in archive.files}
    with np.load(args.context, allow_pickle=True) as archive:
        split = np.asarray(archive["split"]).astype(str)
        train_row = split == "train"
        donor_row = np.asarray(archive["donor"], dtype=np.int64)[train_row]
        context_donors = np.asarray(archive["donor_vocab"]).astype(str)
        context_celltypes = np.asarray(archive["celltype_vocab"]).astype(str)
        context_modules = np.asarray(archive["module_names"]).astype(str)
        age_row = np.asarray(archive["age"], dtype=np.float64)[train_row]
        sex_row = np.asarray(archive["sex"], dtype=np.float64)[train_row]
        pathology_rows = np.column_stack(
            [np.asarray(archive[name], dtype=np.float64)[train_row] for name in AXES]
        )

    donor_names = np.asarray(base["donor_names"]).astype(str)
    celltype_names = np.asarray(base["celltype_names"]).astype(str)
    module_names = np.asarray(base["module_names"]).astype(str)
    if not np.array_equal(donor_names, context_donors):
        raise ValueError("context/base donor vocabulary mismatch")
    if not np.array_equal(celltype_names, context_celltypes):
        raise ValueError("context/base celltype vocabulary mismatch")
    if not np.array_equal(module_names, context_modules):
        raise ValueError("context/base module vocabulary mismatch")

    n_donor = donor_names.size
    n_axis = len(AXES)
    donor_path = np.full((n_donor, n_axis), np.nan, dtype=np.float64)
    donor_age = np.full(n_donor, np.nan, dtype=np.float64)
    donor_sex = np.full(n_donor, np.nan, dtype=np.float64)
    train_donor = np.zeros(n_donor, dtype=bool)
    for donor in np.unique(donor_row):
        rows = donor_row == donor
        train_donor[donor] = True
        donor_age[donor] = np.nanmedian(age_row[rows])
        donor_sex[donor] = np.nanmedian(sex_row[rows])
        for axis in range(n_axis):
            values = pathology_rows[rows, axis]
            finite = values[np.isfinite(values)]
            if finite.size:
                donor_path[donor, axis] = np.median(finite)
    if int(train_donor.sum()) != 64:
        raise ValueError(f"expected 64 train donors, found {int(train_donor.sum())}")

    folds = np.asarray([_fold(name) for name in donor_names], dtype=np.int64)
    residual = np.zeros((n_donor, n_axis), dtype=np.float32)
    residual_valid = np.zeros((n_donor, n_axis), dtype=bool)
    for axis in range(n_axis):
        covariates = np.column_stack(
            (
                np.delete(donor_path, axis, axis=1),
                donor_age,
                donor_sex,
            )
        )
        valid_target = train_donor & np.isfinite(donor_path[:, axis])
        axis_residual = np.full(n_donor, np.nan, dtype=np.float64)
        for held_fold in (0, 1):
            fit = valid_target & (folds != held_fold)
            held = valid_target & (folds == held_fold)
            if int(fit.sum()) < 16 or int(held.sum()) < 8:
                raise ValueError(f"axis {AXES[axis]} has an invalid cross-fit fold")
            prediction = _fit_predict_ridge(
                covariates[fit], donor_path[fit, axis], covariates[held]
            )
            axis_residual[held] = donor_path[held, axis] - prediction
        scale = np.nanstd(axis_residual[valid_target])
        if not np.isfinite(scale) or scale <= 1.0e-8:
            raise ValueError(f"axis {AXES[axis]} residual scale collapsed")
        axis_residual[valid_target] = (
            axis_residual[valid_target] - np.nanmean(axis_residual[valid_target])
        ) / scale
        residual[valid_target, axis] = axis_residual[valid_target].astype(np.float32)
        residual_valid[valid_target, axis] = True

    group_module = np.asarray(base["group_module"], dtype=np.float64)
    group_count = np.asarray(base["group_count"], dtype=np.int64)
    celltype_scale = np.asarray(base["celltype_scale"], dtype=np.float64)
    reliable = np.zeros(
        (celltype_names.size, n_axis, module_names.size), dtype=bool
    )
    reliability_rows = []
    for celltype in range(celltype_names.size):
        y_all = group_module[:, celltype, :] / np.maximum(
            celltype_scale[celltype][None, :], 0.05
        )
        for axis in range(n_axis):
            eligible = (
                train_donor
                & residual_valid[:, axis]
                & (group_count[:, celltype] >= 10)
                & np.isfinite(y_all).all(axis=1)
            )
            coefficients = []
            for fold in (0, 1):
                rows = eligible & (folds == fold)
                if int(rows.sum()) < 8:
                    raise ValueError(
                        f"too few reliable donors for {celltype_names[celltype]} {AXES[axis]}"
                    )
                x = residual[rows, axis].astype(np.float64)
                x = x - x.mean()
                y = y_all[rows] - y_all[rows].mean(axis=0, keepdims=True)
                coefficients.append((x[:, None] * y).sum(axis=0) / max((x * x).sum(), 1.0e-8))
            first, second = coefficients
            same_sign = (first * second) > 0.0
            score = np.minimum(np.abs(first), np.abs(second))
            candidates = np.flatnonzero(same_sign & np.isfinite(score))
            if candidates.size < 32:
                candidates = np.flatnonzero(np.isfinite(score))
            order = candidates[np.argsort(score[candidates])[::-1]]
            keep_n = max(32, int(np.ceil(0.25 * order.size)))
            chosen = order[: min(keep_n, order.size)]
            reliable[celltype, axis, chosen] = True
            reliability_rows.append(
                {
                    "celltype": str(celltype_names[celltype]),
                    "axis": AXES[axis],
                    "eligible_train_donors": int(eligible.sum()),
                    "same_sign_modules": int(same_sign.sum()),
                    "selected_modules": int(chosen.size),
                }
            )

    metadata = json.loads(str(np.asarray(base["metadata_json"]).item()))
    metadata.update(
        {
            "pathology_axis_schema": "kmlee_bam.prism_pathology_axis_rescue.v1",
            "pathology_axis_names": list(AXES),
            "pathology_axis_source_split": "train_only",
            "pathology_axis_official_validation_used": False,
            "pathology_axis_official_test_used": False,
            "pathology_axis_fit_donor_count": int(train_donor.sum()),
            "pathology_axis_crossfit_folds": 2,
            "pathology_axis_crossfit_seed": SEED,
            "pathology_axis_covariates": ["other_pathology_axes", "age", "sex"],
            "pathology_axis_region_policy": "region_balanced_combined_group_module_target",
            "pathology_axis_reliability": "deterministic_train_donor_split_half_same_sign_top_quartile_min32",
            "base_artifact_sha256": _sha256(args.base),
            "context_artifact_sha256": _sha256(args.context),
        }
    )
    base["metadata_json"] = np.asarray(json.dumps(metadata, sort_keys=True))
    base["pathology_residual"] = residual
    base["pathology_residual_valid"] = residual_valid
    base["pathology_reliable_mask"] = reliable
    base["pathology_axis_names"] = np.asarray(AXES)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **base)
    report = {
        "output": str(args.output),
        "sha256": _sha256(args.output),
        "train_donors": int(train_donor.sum()),
        "validation_used": False,
        "test_used": False,
        "reliability": reliability_rows,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
