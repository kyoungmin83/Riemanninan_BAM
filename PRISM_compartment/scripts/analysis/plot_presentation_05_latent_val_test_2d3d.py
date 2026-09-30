#!/usr/bin/env python3
"""Validation-fitted 32-D PRISM latent maps with frozen test projection.

This is deliberately a donor-level analysis.  Cells are first averaged within
donor x cell type, cell types are centred on validation, and the 24 centred
cell-type means are then equally averaged per donor.  For each named pathology
axis, a ridge direction is learned from validation donors only.  Validation
performance is nested leave-one-donor-out; test donors are projected without
refitting.  The 3-D height is a descriptive validation-donor KDE, not a
geodesic, Waddington potential, or empirical cell-level likelihood.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/prism_latent_val_test_mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/prism_latent_val_test_cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge


EXPECTED_CHECKPOINT_SHA256 = (
    "2cdc91ae005e6183cd7aa9f620ab463e0c348db9752247adf320ecc261a99bd0"
)
PATHOLOGY_LABELS = {
    "thal": "Thal",
    "braak": "Braak",
    "cerad": "CERAD",
    "late": "LATE",
    "lewy": "Lewy",
}
ALPHAS = np.logspace(-4.0, 4.0, 17)
VAL_COLOR = "#1769aa"
TEST_COLOR = "#e67e22"
GREEN = "#008f6b"
ORANGE = "#e66b20"
NAVY = "#123b72"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def set_style() -> None:
    korean_font = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    font_family = "DejaVu Sans"
    if korean_font.exists():
        font_manager.fontManager.addfont(str(korean_font))
        font_family = font_manager.FontProperties(fname=str(korean_font)).get_name()
    plt.rcParams.update(
        {
            "font.family": font_family,
            "font.size": 8.2,
            "axes.titlesize": 9.0,
            "axes.labelsize": 7.6,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "figure.facecolor": "white",
            "axes.facecolor": "#fbfcfe",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


@dataclass
class Extract:
    path: Path
    z: np.ndarray
    donor: np.ndarray
    celltype: np.ndarray
    pathology: np.ndarray
    pathology_valid: np.ndarray
    split: np.ndarray
    donor_vocab: np.ndarray
    celltype_vocab: np.ndarray
    pathology_names: np.ndarray
    checkpoint_sha256: str


def load_extract(path: Path) -> Extract:
    data = np.load(path, allow_pickle=True)
    required = {
        "z",
        "donor",
        "celltype",
        "pathology",
        "pathology_valid",
        "split",
        "donor_vocab",
        "celltype_vocab",
        "pathology_names",
        "source_checkpoint_sha256",
    }
    missing = required.difference(data.files)
    if missing:
        raise RuntimeError(f"{path} is missing {sorted(missing)}")
    checkpoint_sha256 = str(data["source_checkpoint_sha256"].item())
    if checkpoint_sha256 != EXPECTED_CHECKPOINT_SHA256:
        raise RuntimeError(
            f"checkpoint mismatch for {path}: {checkpoint_sha256} != "
            f"{EXPECTED_CHECKPOINT_SHA256}"
        )
    return Extract(
        path=path,
        z=np.asarray(data["z"], dtype=np.float64),
        donor=np.asarray(data["donor"], dtype=np.int64),
        celltype=np.asarray(data["celltype"], dtype=np.int64),
        pathology=np.asarray(data["pathology"], dtype=np.float64),
        pathology_valid=np.asarray(data["pathology_valid"], dtype=bool),
        split=np.asarray(data["split"], dtype=object).astype(str),
        donor_vocab=np.asarray(data["donor_vocab"], dtype=object).astype(str),
        celltype_vocab=np.asarray(data["celltype_vocab"], dtype=object).astype(str),
        pathology_names=np.asarray(data["pathology_names"], dtype=object).astype(str),
        checkpoint_sha256=checkpoint_sha256,
    )


@dataclass
class Groups:
    donor: np.ndarray
    donor_name: np.ndarray
    celltype: np.ndarray
    z: np.ndarray
    pathology: np.ndarray
    valid: np.ndarray
    n_cells: np.ndarray


def aggregate_donor_celltype(extract: Extract, split: str, min_cells: int) -> Groups:
    selected = np.flatnonzero(extract.split == split)
    if not len(selected):
        raise RuntimeError(f"split {split!r} is absent from {extract.path}")
    keys = np.column_stack((extract.donor[selected], extract.celltype[selected]))
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    z_sum = np.zeros((len(unique), extract.z.shape[1]), dtype=np.float64)
    count = np.zeros(len(unique), dtype=np.int64)
    np.add.at(z_sum, inverse, extract.z[selected])
    np.add.at(count, inverse, 1)
    keep = count >= min_cells
    unique = unique[keep]
    z_mean = z_sum[keep] / count[keep, None]
    count = count[keep]

    pathology = np.full((len(unique), extract.pathology.shape[1]), np.nan)
    valid = np.zeros_like(pathology, dtype=bool)
    for row, (donor, celltype) in enumerate(unique):
        mask = selected[
            (extract.donor[selected] == donor)
            & (extract.celltype[selected] == celltype)
        ]
        for axis in range(pathology.shape[1]):
            ok = extract.pathology_valid[mask, axis] & np.isfinite(
                extract.pathology[mask, axis]
            )
            if np.any(ok):
                pathology[row, axis] = float(np.mean(extract.pathology[mask, axis][ok]))
                valid[row, axis] = True
    return Groups(
        donor=unique[:, 0].astype(np.int64),
        donor_name=extract.donor_vocab[unique[:, 0]],
        celltype=unique[:, 1].astype(np.int64),
        z=z_mean,
        pathology=pathology,
        valid=valid,
        n_cells=count,
    )


@dataclass
class Donors:
    name: np.ndarray
    z: np.ndarray
    pathology: np.ndarray
    valid: np.ndarray
    n_celltypes: np.ndarray
    n_cells: np.ndarray


def validation_celltype_standardization(
    validation: Groups,
    test: Groups,
    n_celltypes: int,
) -> tuple[Groups, Groups, np.ndarray, np.ndarray]:
    d_z = validation.z.shape[1]
    centers = np.zeros((n_celltypes, d_z), dtype=np.float64)
    for celltype in range(n_celltypes):
        mask = validation.celltype == celltype
        if not np.any(mask):
            raise RuntimeError(f"validation lacks cell type id {celltype}")
        centers[celltype] = validation.z[mask].mean(axis=0)
    val_centered = validation.z - centers[validation.celltype]
    test_centered = test.z - centers[test.celltype]
    scale = val_centered.std(axis=0, ddof=1)
    positive = scale[scale > 1e-8]
    floor = max(float(np.median(positive)) * 0.05, 1e-6)
    scale = np.maximum(scale, floor)

    def replace(groups: Groups, values: np.ndarray) -> Groups:
        return Groups(
            donor=groups.donor,
            donor_name=groups.donor_name,
            celltype=groups.celltype,
            z=values / scale,
            pathology=groups.pathology,
            valid=groups.valid,
            n_cells=groups.n_cells,
        )

    return replace(validation, val_centered), replace(test, test_centered), centers, scale


def equal_celltype_donor_mean(groups: Groups, min_celltypes: int) -> Donors:
    names = np.unique(groups.donor_name)
    rows = []
    for name in names:
        mask = groups.donor_name == name
        n_celltypes = len(np.unique(groups.celltype[mask]))
        if n_celltypes < min_celltypes:
            continue
        axis_values = np.full(groups.pathology.shape[1], np.nan)
        axis_valid = np.zeros(groups.pathology.shape[1], dtype=bool)
        for axis in range(groups.pathology.shape[1]):
            ok = groups.valid[mask, axis] & np.isfinite(groups.pathology[mask, axis])
            if np.any(ok):
                axis_values[axis] = float(np.mean(groups.pathology[mask, axis][ok]))
                axis_valid[axis] = True
        rows.append(
            (
                name,
                groups.z[mask].mean(axis=0),
                axis_values,
                axis_valid,
                n_celltypes,
                int(groups.n_cells[mask].sum()),
            )
        )
    if not rows:
        raise RuntimeError("no donors survive the cell-type completeness gate")
    return Donors(
        name=np.asarray([row[0] for row in rows], dtype=object),
        z=np.stack([row[1] for row in rows]),
        pathology=np.stack([row[2] for row in rows]),
        valid=np.stack([row[3] for row in rows]),
        n_celltypes=np.asarray([row[4] for row in rows], dtype=np.int64),
        n_cells=np.asarray([row[5] for row in rows], dtype=np.int64),
    )


def ridge_predictions(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, alpha: float):
    model = Ridge(alpha=float(alpha), fit_intercept=True)
    model.fit(x_train, y_train)
    return model.predict(x_test), model


def balanced_rank_folds(y: np.ndarray, n_splits: int) -> np.ndarray:
    """Deterministic severity-balanced folds without random/test-informed tuning."""
    if n_splits < 2 or len(y) < 2 * n_splits:
        raise ValueError(f"cannot form {n_splits} folds from {len(y)} donors")
    order = np.argsort(y, kind="mergesort")
    pattern = np.concatenate((np.arange(n_splits), np.arange(n_splits - 1, -1, -1)))
    folds = np.empty(len(y), dtype=np.int64)
    folds[order] = np.resize(pattern, len(y))
    return folds


def crossfit_predictions(
    x: np.ndarray,
    y: np.ndarray,
    alpha: float,
    *,
    n_splits: int,
) -> np.ndarray:
    folds = balanced_rank_folds(y, n_splits)
    predictions = np.empty(len(y), dtype=np.float64)
    for fold in range(n_splits):
        test = folds == fold
        train = ~test
        predictions[test] = ridge_predictions(
            x[train], y[train], x[test], alpha
        )[0]
    return predictions


def select_alpha(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_splits: int,
) -> tuple[float, np.ndarray, list[dict]]:
    audits = []
    best = None
    for alpha in ALPHAS:
        prediction = crossfit_predictions(
            x, y, float(alpha), n_splits=n_splits
        )
        mse = float(np.mean(np.square(prediction - y)))
        audits.append(
            {
                "alpha": float(alpha),
                "crossfit_mse": mse,
                "n_splits": int(n_splits),
            }
        )
        key = (mse, float(alpha))
        if best is None or key < best[0]:
            best = (key, float(alpha), prediction)
    assert best is not None
    return best[1], best[2], audits


def nested_balanced_predictions(
    x: np.ndarray,
    y: np.ndarray,
    *,
    outer_splits: int = 4,
    inner_splits: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    folds = balanced_rank_folds(y, outer_splits)
    prediction = np.empty(len(y), dtype=np.float64)
    alpha_used = np.empty(outer_splits, dtype=np.float64)
    for fold in range(outer_splits):
        test = folds == fold
        train = ~test
        alpha, _, _ = select_alpha(
            x[train], y[train], n_splits=inner_splits
        )
        prediction[test] = ridge_predictions(
            x[train], y[train], x[test], alpha
        )[0]
        alpha_used[fold] = alpha
    return prediction, alpha_used


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return float("nan")
    return float(spearmanr(x, y).statistic)


def bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_boot: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        index = rng.integers(0, len(y), size=len(y))
        value = safe_spearman(x[index], y[index])
        if np.isfinite(value):
            values.append(value)
    if len(values) < max(50, n_boot // 20):
        return float("nan"), float("nan")
    return tuple(np.quantile(values, (0.025, 0.975)).tolist())


def permutation_pvalue(x: np.ndarray, y: np.ndarray, *, n_perm: int, seed: int) -> float:
    observed = abs(safe_spearman(x, y))
    if not np.isfinite(observed):
        return float("nan")
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_perm):
        exceed += abs(safe_spearman(x, rng.permutation(y))) >= observed
    return float((exceed + 1) / (n_perm + 1))


def orthogonal_plane(x_val: np.ndarray, x_test: np.ndarray, coefficient: np.ndarray):
    e1 = np.asarray(coefficient, dtype=np.float64)
    if np.linalg.norm(e1) < 1e-12:
        raise RuntimeError("ridge coefficient is numerically zero")
    e1 /= np.linalg.norm(e1)
    val_x = x_val @ e1
    residual = x_val - np.outer(val_x, e1)
    _, _, vh = np.linalg.svd(residual - residual.mean(axis=0), full_matrices=False)
    e2 = vh[0]
    e2 -= float(e2 @ e1) * e1
    e2 /= max(np.linalg.norm(e2), 1e-12)
    if e2[np.argmax(np.abs(e2))] < 0:
        e2 *= -1.0

    def project(value: np.ndarray) -> np.ndarray:
        return np.column_stack((value @ e1, value @ e2))

    return project(x_val), project(x_test), e1, e2


def density_surface(points: np.ndarray, severity: np.ndarray, test_points: np.ndarray):
    all_points = np.vstack((points, test_points))
    span = np.ptp(all_points, axis=0)
    padding = np.maximum(0.24 * span, 0.35)
    lower = all_points.min(axis=0) - padding
    upper = all_points.max(axis=0) + padding
    gx = np.linspace(lower[0], upper[0], 78)
    gy = np.linspace(lower[1], upper[1], 66)
    xx, yy = np.meshgrid(gx, gy, indexing="xy")
    grid = np.stack((xx, yy), axis=-1)
    std = np.std(points, axis=0, ddof=1)
    bandwidth = np.maximum(0.58 * std, 0.18 * np.maximum(span, 1e-6))
    delta = (grid[:, :, None, :] - points[None, None, :, :]) / bandwidth
    weight = np.exp(-0.5 * np.sum(np.square(delta), axis=-1))
    density = weight.mean(axis=-1)
    density /= max(float(density.max()), 1e-12)
    severity_field = np.sum(weight * severity[None, None, :], axis=-1) / np.maximum(
        weight.sum(axis=-1), 1e-12
    )
    return gx, gy, xx, yy, density, np.clip(severity_field, 0.0, 1.0), bandwidth


def lifted_height(gx: np.ndarray, gy: np.ndarray, density: np.ndarray, points: np.ndarray):
    interpolator = RegularGridInterpolator(
        (gy, gx), density, bounds_error=False, fill_value=0.0
    )
    return np.asarray(interpolator(points[:, [1, 0]]), dtype=np.float64)


def extreme_centres(points: np.ndarray, severity: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    k = max(2, min(4, len(severity) // 3))
    order = np.argsort(severity)
    return points[order[:k]].mean(axis=0), points[order[-k:]].mean(axis=0), k


def draw_2d(
    ax,
    name: str,
    val_points: np.ndarray,
    val_y: np.ndarray,
    test_points: np.ndarray,
    test_y: np.ndarray,
    surface: tuple,
    metrics: dict,
    cmap,
    norm,
    show_legend: bool,
):
    gx, gy, xx, yy, density, severity_field, _ = surface
    mask = density >= 0.08
    shown = np.ma.masked_where(~mask, severity_field)
    ax.contourf(xx, yy, shown, levels=np.linspace(0.0, 1.0, 13), cmap=cmap, norm=norm, alpha=0.30)
    ax.contour(xx, yy, density, levels=(0.18, 0.35, 0.55, 0.75), colors="#65758b", linewidths=(0.7, 1.0, 1.3, 1.7), alpha=0.72)
    ax.scatter(
        val_points[:, 0], val_points[:, 1], s=58, marker="o",
        facecolors="white", edgecolors=cmap(norm(val_y)), linewidths=1.7,
        zorder=8,
    )
    ax.scatter(
        test_points[:, 0], test_points[:, 1], s=62, marker="D",
        c=test_y, cmap=cmap, norm=norm, edgecolors="#13243b", linewidths=0.75,
        zorder=9,
    )
    val_low, val_high, _ = extreme_centres(val_points, val_y)
    test_low, test_high, _ = extreme_centres(test_points, test_y)
    ax.annotate("", xy=val_high, xytext=val_low, arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2.5))
    ax.annotate("", xy=test_high, xytext=test_low, arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.2, ls="--"))
    ax.set_xlabel("validation-fitted latent pathology score")
    ax.set_ylabel("orthogonal latent variation")
    val_ci = metrics["validation_nested_4fold_rho_ci"]
    test_ci = metrics["test_rho_ci"]
    ax.set_title(
        f"{name} | val nested 4-fold ρ={metrics['validation_nested_4fold_rho']:+.2f} "
        f"[{val_ci[0]:+.2f},{val_ci[1]:+.2f}]\n"
        f"test ρ={metrics['test_rho']:+.2f} [{test_ci[0]:+.2f},{test_ci[1]:+.2f}]  "
        f"cos={metrics['val_test_direction_cosine']:+.2f}",
        loc="left", fontweight="bold",
    )
    ax.grid(color="#d5dde8", lw=0.55, alpha=0.65)
    if show_legend:
        ax.legend(
            handles=(
                Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=VAL_COLOR, markeredgewidth=1.7, label="Validation donor (16명)"),
                Line2D([0], [0], marker="D", color="none", markerfacecolor=TEST_COLOR, markeredgecolor="#13243b", label="Test donor (9명, 무재학습)"),
                Line2D([0], [0], color=GREEN, lw=2.5, label="Validation 저→고병리 중심"),
                Line2D([0], [0], color=ORANGE, lw=2.2, ls="--", label="Test 저→고병리 중심"),
            ),
            loc="best", fontsize=5.8, framealpha=0.94,
        )


def draw_3d(
    ax,
    name: str,
    val_points: np.ndarray,
    val_y: np.ndarray,
    test_points: np.ndarray,
    test_y: np.ndarray,
    surface: tuple,
    cmap,
    norm,
):
    gx, gy, xx, yy, density, severity_field, _ = surface
    colors = cmap(norm(severity_field))
    colors[..., 3] = 0.78
    ax.plot_surface(
        xx[::2, ::2], yy[::2, ::2], density[::2, ::2],
        facecolors=colors[::2, ::2], rstride=1, cstride=1,
        linewidth=0.12, edgecolor=(0.2, 0.2, 0.2, 0.18),
        shade=False, antialiased=True,
    )
    val_z = lifted_height(gx, gy, density, val_points) + 0.055
    test_z = lifted_height(gx, gy, density, test_points) + 0.085
    ax.scatter(
        val_points[:, 0], val_points[:, 1], val_z, s=34, marker="o",
        facecolors="white", edgecolors=cmap(norm(val_y)), linewidths=1.15,
        depthshade=False,
    )
    ax.scatter(
        test_points[:, 0], test_points[:, 1], test_z, s=42, marker="D",
        c=test_y, cmap=cmap, norm=norm, edgecolors="#13243b", linewidths=0.55,
        depthshade=False,
    )
    for points, severity, color, linestyle in (
        (val_points, val_y, GREEN, "-"),
        (test_points, test_y, ORANGE, "--"),
    ):
        low, high, _ = extreme_centres(points, severity)
        line = np.linspace(low, high, 80)
        line_z = lifted_height(gx, gy, density, line) + 0.10
        ax.plot(line[:, 0], line[:, 1], line_z, color=color, ls=linestyle, lw=2.4)
    ax.set_xlabel("latent pathology score", labelpad=1)
    ax.set_ylabel("orthogonal", labelpad=1)
    ax.set_zlabel("validation donor density", labelpad=1)
    ax.set_zlim(0.0, 1.18)
    ax.set_zticks((0.0, 0.5, 1.0))
    ax.view_init(elev=31, azim=-61)
    ax.set_title(f"{name} | 같은 donor 지도를 3-D 밀도로 올림", fontweight="bold", pad=1)
    ax.grid(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-extract", type=Path, required=True)
    parser.add_argument("--test-extract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-cells-group", type=int, default=10)
    parser.add_argument("--min-celltypes-donor", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=4000)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()

    set_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    validation_extract = load_extract(args.validation_extract)
    test_extract = load_extract(args.test_extract)
    if not np.array_equal(validation_extract.celltype_vocab, test_extract.celltype_vocab):
        raise RuntimeError("cell-type vocabularies do not match")
    if not np.array_equal(validation_extract.pathology_names, test_extract.pathology_names):
        raise RuntimeError("pathology vocabularies do not match")

    val_groups_raw = aggregate_donor_celltype(validation_extract, "val", args.min_cells_group)
    test_groups_raw = aggregate_donor_celltype(test_extract, "test", args.min_cells_group)
    val_groups, test_groups, celltype_centers, latent_scale = validation_celltype_standardization(
        val_groups_raw,
        test_groups_raw,
        len(validation_extract.celltype_vocab),
    )
    val_donors = equal_celltype_donor_mean(val_groups, args.min_celltypes_donor)
    test_donors = equal_celltype_donor_mean(test_groups, args.min_celltypes_donor)
    if len(val_donors.name) != 16 or len(test_donors.name) != 9:
        raise RuntimeError(
            f"expected 16 validation and 9 test donors, got "
            f"{len(val_donors.name)} and {len(test_donors.name)}"
        )

    cmap = matplotlib.colormaps["coolwarm"]
    norm = Normalize(vmin=0.0, vmax=1.0)
    fig = plt.figure(figsize=(22.4, 10.8))
    metrics_rows = []
    coordinate_rows = []
    exact = {
        "validation_donor": val_donors.name,
        "test_donor": test_donors.name,
        "validation_z_equal_celltype_mean": val_donors.z.astype(np.float32),
        "test_z_equal_celltype_mean": test_donors.z.astype(np.float32),
        "pathology_names": validation_extract.pathology_names,
        "celltype_centers": celltype_centers.astype(np.float32),
        "latent_scale": latent_scale.astype(np.float32),
    }
    all_metrics = []

    for axis, raw_name in enumerate(validation_extract.pathology_names):
        display_name = PATHOLOGY_LABELS.get(raw_name.lower(), raw_name)
        val_ok = val_donors.valid[:, axis] & np.isfinite(val_donors.pathology[:, axis])
        test_ok = test_donors.valid[:, axis] & np.isfinite(test_donors.pathology[:, axis])
        if not np.all(val_ok) or not np.all(test_ok):
            raise RuntimeError(f"missing {display_name} pathology values")
        x_val = val_donors.z[val_ok]
        y_val = val_donors.pathology[val_ok, axis]
        x_test = test_donors.z[test_ok]
        y_test = test_donors.pathology[test_ok, axis]

        nested_prediction, outer_alpha = nested_balanced_predictions(
            x_val, y_val, outer_splits=4, inner_splits=3
        )
        full_alpha, ordinary_crossfit, alpha_audit = select_alpha(
            x_val, y_val, n_splits=4
        )
        test_prediction, model = ridge_predictions(x_val, y_val, x_test, full_alpha)
        val_rho = safe_spearman(nested_prediction, y_val)
        test_rho = safe_spearman(test_prediction, y_test)
        val_ci = bootstrap_spearman(
            nested_prediction, y_val, n_boot=args.bootstrap, seed=args.seed + 100 * axis
        )
        test_ci = bootstrap_spearman(
            test_prediction, y_test, n_boot=args.bootstrap, seed=args.seed + 100 * axis + 1
        )
        test_p = permutation_pvalue(
            test_prediction, y_test, n_perm=args.permutations, seed=args.seed + 100 * axis + 2
        )
        val_points, test_points, e1, e2 = orthogonal_plane(
            x_val, x_test, model.coef_
        )
        val_low, val_high, val_k = extreme_centres(x_val, y_val)
        test_low, test_high, test_k = extreme_centres(x_test, y_test)
        val_direction = val_high - val_low
        test_direction = test_high - test_low
        cosine = float(
            np.dot(val_direction, test_direction)
            / max(np.linalg.norm(val_direction) * np.linalg.norm(test_direction), 1e-12)
        )
        surface = density_surface(val_points, y_val, test_points)
        metric = {
            "pathology": display_name,
            "n_validation_donors": int(len(y_val)),
            "n_test_donors": int(len(y_test)),
            "ridge_alpha_selected_validation_4fold": float(full_alpha),
            "validation_nested_4fold_rho": val_rho,
            "validation_nested_4fold_rho_ci": val_ci,
            "test_rho": test_rho,
            "test_rho_ci": test_ci,
            "test_permutation_p_two_sided": test_p,
            "val_test_direction_cosine": cosine,
            "validation_extreme_k": val_k,
            "test_extreme_k": test_k,
            "outer_fold_alphas": outer_alpha.tolist(),
            "alpha_audit_full_validation": alpha_audit,
            "warning": "Named pathology scores are correlated; this is marginal association, not a unique-axis or causal effect.",
        }
        all_metrics.append(metric)
        metrics_rows.append(metric)

        exact[f"{raw_name}_ridge_coefficient"] = np.asarray(model.coef_, dtype=np.float32)
        exact[f"{raw_name}_plane_e1"] = e1.astype(np.float32)
        exact[f"{raw_name}_plane_e2"] = e2.astype(np.float32)
        exact[f"{raw_name}_validation_xy"] = val_points.astype(np.float32)
        exact[f"{raw_name}_test_xy"] = test_points.astype(np.float32)
        exact[f"{raw_name}_validation_pathology"] = y_val.astype(np.float32)
        exact[f"{raw_name}_test_pathology"] = y_test.astype(np.float32)
        exact[f"{raw_name}_validation_nested_prediction"] = nested_prediction.astype(np.float32)
        exact[f"{raw_name}_test_prediction"] = test_prediction.astype(np.float32)

        for split, donors, points, observed, predicted in (
            ("validation", val_donors.name[val_ok], val_points, y_val, nested_prediction),
            ("test", test_donors.name[test_ok], test_points, y_test, test_prediction),
        ):
            for donor, point, obs, pred in zip(donors, points, observed, predicted):
                coordinate_rows.append(
                    {
                        "pathology": display_name,
                        "split": split,
                        "donor": str(donor),
                        "latent_pathology_score": float(point[0]),
                        "orthogonal_latent_variation": float(point[1]),
                        "observed_pathology": float(obs),
                        "cross_validated_or_frozen_prediction": float(pred),
                    }
                )

        ax2 = fig.add_subplot(2, 5, axis + 1)
        draw_2d(
            ax2,
            display_name,
            val_points,
            y_val,
            test_points,
            y_test,
            surface,
            metric,
            cmap,
            norm,
            show_legend=(axis == 0),
        )
        if axis > 0:
            ax2.set_ylabel("")
            ax2.set_yticklabels([])
        ax3 = fig.add_subplot(2, 5, axis + 6, projection="3d")
        draw_3d(
            ax3,
            display_name,
            val_points,
            y_val,
            test_points,
            y_test,
            surface,
            cmap,
            norm,
        )

    positive_test = int(sum(metric["test_rho"] > 0 for metric in all_metrics))
    test_ci_excludes_zero = int(
        sum(
            metric["test_rho_ci"][0] > 0.0 or metric["test_rho_ci"][1] < 0.0
            for metric in all_metrics
        )
    )
    median_test = float(np.nanmedian([metric["test_rho"] for metric in all_metrics]))
    median_cosine = float(
        np.nanmedian([metric["val_test_direction_cosine"] for metric in all_metrics])
    )
    fig.suptitle(
        "PRISM 32차원 zclean에서 validation이 정한 병리 방향을 test donor에 고정 투영했다\n"
        "donor×세포타입 평균 → 세포타입 동일비중 donor 지문 | 위: 실제 donor 좌표 | 아래: 같은 validation donor 밀도의 3-D lift",
        fontsize=15.2,
        fontweight="bold",
        color="#102b55",
        y=0.987,
    )
    color_axis = fig.add_axes((0.970, 0.57, 0.008, 0.255))
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    colorbar = fig.colorbar(sm, cax=color_axis)
    colorbar.set_label("관찰 병리 중증도 (0–1)", fontsize=7.2)
    fig.text(
        0.5,
        0.050,
        f"핵심 | test 중앙 ρ={median_test:+.2f}, 95% CI가 0을 제외한 축={test_ci_excludes_zero}/5, "
        f"validation–test 방향 중앙 cosine={median_cosine:+.2f}: 명명된 병리방향의 재현성은 아직 확정되지 않았다.",
        ha="center",
        va="center",
        fontsize=10.3,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.65", facecolor=NAVY, edgecolor=NAVY),
    )
    fig.text(
        0.5,
        0.015,
        "주의: 3-D 높이=validation donor KDE(설명용)이며 geodesic·시간경로·Waddington potential이 아니다. "
        "다섯 병리축은 서로 상관되어 있으므로 각 패널을 고유 병리 인과효과로 해석하지 않는다.",
        ha="center",
        fontsize=8.0,
        color="#5b6574",
    )
    fig.subplots_adjust(left=0.035, right=0.955, top=0.875, bottom=0.105, wspace=0.24, hspace=0.30)

    stem = args.output_dir / "figure_05_latent_val_test_2d3d"
    fig.savefig(stem.with_suffix(".png"), dpi=190, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    import csv

    with (args.output_dir / "figure_05_axis_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "pathology",
            "n_validation_donors",
            "n_test_donors",
            "ridge_alpha_selected_validation_4fold",
            "validation_nested_4fold_rho",
            "validation_nested_4fold_rho_ci_low",
            "validation_nested_4fold_rho_ci_high",
            "test_rho",
            "test_rho_ci_low",
            "test_rho_ci_high",
            "test_permutation_p_two_sided",
            "val_test_direction_cosine",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for metric in metrics_rows:
            writer.writerow(
                {
                    "pathology": metric["pathology"],
                    "n_validation_donors": metric["n_validation_donors"],
                    "n_test_donors": metric["n_test_donors"],
                    "ridge_alpha_selected_validation_4fold": metric["ridge_alpha_selected_validation_4fold"],
                    "validation_nested_4fold_rho": metric["validation_nested_4fold_rho"],
                    "validation_nested_4fold_rho_ci_low": metric["validation_nested_4fold_rho_ci"][0],
                    "validation_nested_4fold_rho_ci_high": metric["validation_nested_4fold_rho_ci"][1],
                    "test_rho": metric["test_rho"],
                    "test_rho_ci_low": metric["test_rho_ci"][0],
                    "test_rho_ci_high": metric["test_rho_ci"][1],
                    "test_permutation_p_two_sided": metric["test_permutation_p_two_sided"],
                    "val_test_direction_cosine": metric["val_test_direction_cosine"],
                }
            )
    with (args.output_dir / "figure_05_donor_coordinates.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(coordinate_rows[0])
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(coordinate_rows)
    np.savez_compressed(args.output_dir / "figure_05_exact_arrays.npz", **exact)

    pathology_correlation = {
        "validation": np.corrcoef(val_donors.pathology.T),
        "test": np.corrcoef(test_donors.pathology.T),
    }
    manifest = {
        "schema_version": "prism.presentation.latent_val_test_2d3d.v1",
        "model": "frozen PRISM personal-rank2 epoch 20",
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "validation_extract": {
            "path": str(args.validation_extract.resolve()),
            "sha256": sha256_file(args.validation_extract),
        },
        "test_extract": {
            "path": str(args.test_extract.resolve()),
            "sha256": sha256_file(args.test_extract),
        },
        "analysis_contract": {
            "analysis_level": "independent donor",
            "cell_aggregation": "mean within donor x cell type",
            "celltype_aggregation": "validation-centred, dimension-scaled, equal-weight mean across cell types",
            "minimum_cells_per_donor_celltype": args.min_cells_group,
            "minimum_celltypes_per_donor": args.min_celltypes_donor,
            "pathology_direction": "ridge fit on validation only",
            "alpha_selection": "validation severity-balanced 4-fold over fixed logspace(-4,4,17)",
            "validation_metric": "nested severity-balanced 4-fold donor Spearman; 3-fold inner selection",
            "test_metric": "frozen validation fit projected to test without refit",
            "three_dimensional_height": "descriptive Gaussian KDE of 16 validation donor coordinates",
            "test_status": "exploratory posthoc; not prespecified before earlier test access",
            "unique_axis_adjustment": "none; named pathology scores are correlated",
        },
        "counts": {
            "validation_cells": int(np.sum(validation_extract.split == "val")),
            "test_cells": int(np.sum(test_extract.split == "test")),
            "validation_donor_celltype_groups": int(len(val_groups_raw.donor)),
            "test_donor_celltype_groups": int(len(test_groups_raw.donor)),
            "validation_donors": int(len(val_donors.name)),
            "test_donors": int(len(test_donors.name)),
        },
        "metrics": all_metrics,
        "pathology_correlation": pathology_correlation,
        "summary": {
            "positive_test_axes": positive_test,
            "test_axes_with_bootstrap_ci_excluding_zero": test_ci_excludes_zero,
            "total_axes": len(all_metrics),
            "median_test_spearman": median_test,
            "median_val_test_direction_cosine": median_cosine,
        },
        "outputs": {
            "png": stem.with_suffix(".png").name,
            "pdf": stem.with_suffix(".pdf").name,
            "svg": stem.with_suffix(".svg").name,
            "metrics_tsv": "figure_05_axis_metrics.tsv",
            "coordinates_tsv": "figure_05_donor_coordinates.tsv",
            "exact_arrays": "figure_05_exact_arrays.npz",
        },
    }
    with (args.output_dir / "figure_05_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(json_safe(manifest), handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps(json_safe(manifest["summary"]), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
