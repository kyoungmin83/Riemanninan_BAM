#!/usr/bin/env python3
"""Validate frozen common-only pathology directions against observed module slopes.

The primary comparison deliberately uses one common output unit:

1. Decode the frozen common-only normal and one-axis disease endpoints to
   expected ordinal gene tiers.
2. Project the endpoint difference back to the exact tokenizer module
   coordinate with the frozen activity-weight matrix.
3. Estimate observed donor-level module slopes for all five pathology axes in
   one simultaneous OLS model within each cell type.
4. Compare the two 404-module directions with Spearman correlation.

The 9-donor test analysis is descriptive/exploratory: the locked-test plan
pre-specified that per-axis donor-disjoint inference was underpowered below 12
donors.  This script never trains or changes the frozen checkpoint.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch
from scipy.stats import rankdata


PRISM_ROOT = Path(__file__).resolve().parents[2]
KMLEE_ROOT = PRISM_ROOT.parent / "kmlee_bam"
PRESENTATION_ROOT = PRISM_ROOT / "analysis_outputs" / "presentation_results_20260819"

DEFAULT_GEOMETRY = (
    KMLEE_ROOT
    / "analysis_outputs"
    / "prism_rank2_figure11_13_3d_20260811"
    / "geometry"
    / "prism_fingerprints_exact.npz"
)
DEFAULT_COMMON = (
    KMLEE_ROOT
    / "analysis_outputs"
    / "prism_rank2_figure11_13_3d_20260811"
    / "common_paths"
    / "common_rgf_exact_arrays.npz"
)
DEFAULT_VALIDATION = (
    KMLEE_ROOT
    / "analysis_outputs"
    / "prism_rank2_rank4_geometry_reaudit_20260810"
    / "stage0c"
    / "rank2_e20_fingerprint.npz"
)
DEFAULT_TEST = (
    PRESENTATION_ROOT
    / "00_locked_test_rank2_e20"
    / "rank2_e20_test_fingerprint.npz"
)
DEFAULT_WEIGHT = (
    PRESENTATION_ROOT
    / "04_common_axis_observed_validation"
    / "inputs"
    / "activity_weight_kme_or_l2_membership.npz"
)
DEFAULT_OUTPUT = PRESENTATION_ROOT / "04_common_axis_observed_validation"

AXES = ("Thal", "Braak", "CERAD", "LATE", "Lewy")
AXIS_COLORS = {
    "Thal": "#168A83",
    "Braak": "#E67E32",
    "CERAD": "#D6A119",
    "LATE": "#7551A8",
    "Lewy": "#B23A7A",
}
COLORS = {
    "navy": "#123B73",
    "blue": "#1D66B2",
    "teal": "#15958F",
    "orange": "#E87522",
    "red": "#C53B32",
    "dark": "#172033",
    "muted": "#667085",
    "grid": "#D7E0EA",
    "pale_blue": "#EEF5FB",
    "pale_orange": "#FFF5EC",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    font_paths = [
        Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"),
        Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"),
    ]
    for font_path in font_paths:
        if font_path.exists():
            font_manager.fontManager.addfont(font_path)
    family = (
        font_manager.FontProperties(fname=font_paths[0]).get_name()
        if font_paths[0].exists()
        else "DejaVu Sans"
    )
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [family, "DejaVu Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def index_exact(source: np.ndarray, target: np.ndarray, label: str) -> np.ndarray:
    source = np.asarray(source, dtype=str)
    target = np.asarray(target, dtype=str)
    if len(set(source)) != len(source):
        raise ValueError(f"duplicate {label} names in source")
    lookup = {value: index for index, value in enumerate(source)}
    missing = [value for value in target if value not in lookup]
    if missing:
        raise ValueError(f"missing {label} names: {missing[:5]}")
    return np.asarray([lookup[value] for value in target], dtype=np.int64)


def standardize(matrix: np.ndarray) -> np.ndarray | None:
    matrix = np.asarray(matrix, dtype=np.float64)
    mean = np.nanmean(matrix, axis=0)
    sd = np.nanstd(matrix, axis=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 1e-12):
        return None
    result = (matrix - mean) / sd
    if np.any(~np.isfinite(result)):
        return None
    return result


def fit_partial_slopes(
    observed: np.ndarray,
    pathology: np.ndarray,
    donor_indices: np.ndarray | None = None,
    extra_covariates: np.ndarray | None = None,
    ridge: float = 0.0,
) -> np.ndarray | None:
    """Return simultaneous pathology slopes [cell type, axis, module]."""
    if donor_indices is None:
        donor_indices = np.arange(len(pathology), dtype=np.int64)
    donor_indices = np.asarray(donor_indices, dtype=np.int64)
    design_raw = pathology[donor_indices]
    if extra_covariates is not None:
        design_raw = np.column_stack([design_raw, extra_covariates[donor_indices]])
    design = standardize(design_raw)
    if design is None:
        return None
    n_axes = pathology.shape[1]
    slopes = np.full(
        (observed.shape[1], n_axes, observed.shape[2]), np.nan, dtype=np.float64
    )
    for celltype in range(observed.shape[1]):
        response = observed[donor_indices, celltype]
        finite = np.all(np.isfinite(response), axis=1)
        if int(finite.sum()) < n_axes + 1:
            continue
        x = design[finite]
        y = response[finite] - np.nanmean(response[finite], axis=0, keepdims=True)
        gram = x.T @ x + float(ridge) * np.eye(x.shape[1])
        beta = np.linalg.pinv(gram) @ x.T @ y
        slopes[celltype] = beta[:n_axes]
    return slopes


def spearman_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    if left.shape != right.shape or left.ndim != 3:
        raise ValueError("Spearman inputs must have equal [cell type, axis, module] shape")
    complete = np.all(np.isfinite(left) & np.isfinite(right), axis=-1)
    safe_left = np.where(np.isfinite(left), left, 0.0)
    safe_right = np.where(np.isfinite(right), right, 0.0)
    left_rank = rankdata(safe_left, axis=-1, method="average")
    right_rank = rankdata(safe_right, axis=-1, method="average")
    left_rank = left_rank - np.mean(left_rank, axis=-1, keepdims=True)
    right_rank = right_rank - np.mean(right_rank, axis=-1, keepdims=True)
    numerator = np.sum(left_rank * right_rank, axis=-1)
    denominator = np.sqrt(
        np.sum(np.square(left_rank), axis=-1)
        * np.sum(np.square(right_rank), axis=-1)
    )
    result = numerator / np.where(denominator > 0, denominator, np.nan)
    return np.where(complete, result, np.nan)


def sign_agreement_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    finite = np.isfinite(left) & np.isfinite(right) & (left != 0) & (right != 0)
    agree = finite & (np.sign(left) == np.sign(right))
    count = finite.sum(axis=-1)
    return agree.sum(axis=-1) / np.where(count > 0, count, np.nan)


def axis_median(values: np.ndarray) -> np.ndarray:
    return np.nanmedian(values, axis=0)


def bootstrap_model_observed(
    observed: np.ndarray,
    pathology: np.ndarray,
    model_direction: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    output = np.full((n_bootstrap, pathology.shape[1]), np.nan, dtype=np.float64)
    accepted = 0
    attempts = 0
    max_attempts = n_bootstrap * 40
    while accepted < n_bootstrap and attempts < max_attempts:
        attempts += 1
        sample = rng.integers(0, len(pathology), size=len(pathology))
        slopes = fit_partial_slopes(observed, pathology, sample)
        if slopes is None or np.all(~np.isfinite(slopes)):
            continue
        output[accepted] = axis_median(spearman_rows(model_direction, slopes))
        accepted += 1
    if accepted < n_bootstrap:
        raise RuntimeError(f"only {accepted}/{n_bootstrap} model-observed bootstraps succeeded")
    return output


def bootstrap_cell_half(
    observed_a: np.ndarray,
    observed_b: np.ndarray,
    pathology: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    output = np.full((n_bootstrap, pathology.shape[1]), np.nan, dtype=np.float64)
    accepted = 0
    attempts = 0
    max_attempts = n_bootstrap * 40
    while accepted < n_bootstrap and attempts < max_attempts:
        attempts += 1
        sample = rng.integers(0, len(pathology), size=len(pathology))
        slope_a = fit_partial_slopes(observed_a, pathology, sample)
        slope_b = fit_partial_slopes(observed_b, pathology, sample)
        if slope_a is None or slope_b is None:
            continue
        output[accepted] = axis_median(spearman_rows(slope_a, slope_b))
        accepted += 1
    if accepted < n_bootstrap:
        raise RuntimeError(f"only {accepted}/{n_bootstrap} cell-half bootstraps succeeded")
    return output


def bootstrap_cross_split(
    validation: np.ndarray,
    validation_pathology: np.ndarray,
    test: np.ndarray,
    test_pathology: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    output = np.full((n_bootstrap, validation_pathology.shape[1]), np.nan, dtype=np.float64)
    accepted = 0
    attempts = 0
    max_attempts = n_bootstrap * 60
    while accepted < n_bootstrap and attempts < max_attempts:
        attempts += 1
        val_sample = rng.integers(0, len(validation_pathology), size=len(validation_pathology))
        test_sample = rng.integers(0, len(test_pathology), size=len(test_pathology))
        val_slope = fit_partial_slopes(validation, validation_pathology, val_sample)
        test_slope = fit_partial_slopes(test, test_pathology, test_sample)
        if val_slope is None or test_slope is None:
            continue
        output[accepted] = axis_median(spearman_rows(val_slope, test_slope))
        accepted += 1
    if accepted < n_bootstrap:
        raise RuntimeError(f"only {accepted}/{n_bootstrap} cross-split bootstraps succeeded")
    return output


def interval(bootstrap: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.nanquantile(bootstrap, 0.025, axis=0),
        np.nanquantile(bootstrap, 0.975, axis=0),
    )


def draw_interval_points(
    axis: plt.Axes,
    y: np.ndarray,
    point: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    color: str,
    label: str,
    filled: bool,
    offset: float,
) -> None:
    yy = y + offset
    axis.hlines(yy, low, high, color=color, linewidth=2.1, alpha=0.9, zorder=2)
    face = color if filled else "white"
    axis.scatter(
        point,
        yy,
        s=72,
        facecolor=face,
        edgecolor=color,
        linewidth=2.0,
        label=label,
        zorder=3,
    )


def draw_figure(
    summary: pd.DataFrame,
    output: Path,
    n_validation: int,
    n_test: int,
    n_modules: int,
    max_corr_validation: float,
    max_corr_test: float,
) -> None:
    configure_style()
    figure = plt.figure(figsize=(16, 9), facecolor="white")
    grid = figure.add_gridspec(
        1, 2, left=0.075, right=0.965, bottom=0.235, top=0.72, width_ratios=[1.12, 1.0], wspace=0.18
    )
    ax_left = figure.add_subplot(grid[0, 0])
    ax_right = figure.add_subplot(grid[0, 1])

    figure.text(
        0.055,
        0.935,
        "개별 병리축 common 지문은 현재 donor 수에서 독립적으로 재현되지 않았다",
        fontsize=25,
        fontweight=900,
        color=COLORS["navy"],
        va="top",
    )
    figure.text(
        0.055,
        0.865,
        (
            "Frozen personal-rank2 epoch 20 | common-only 0→1 방향 vs 관측 five-axis partial slope | "
            f"validation {n_validation}명 · test {n_test}명 · {n_modules}개 모듈"
        ),
        fontsize=13.2,
        color=COLORS["muted"],
        va="top",
    )
    figure.text(
        0.945,
        0.91,
        "TEST PER-AXIS = EXPLORATORY",
        ha="right",
        va="center",
        fontsize=11.2,
        fontweight=800,
        color=COLORS["orange"],
        bbox={"boxstyle": "round,pad=0.34,rounding_size=0.12", "fc": "#FFF8F1", "ec": "#E9A56D", "lw": 1.1},
    )

    y = np.arange(len(AXES), dtype=float)
    validation = summary.loc[summary["metric"].eq("model_observed_validation")].set_index("axis").loc[list(AXES)]
    test = summary.loc[summary["metric"].eq("model_observed_test")].set_index("axis").loc[list(AXES)]
    draw_interval_points(
        ax_left,
        y,
        validation["point"].to_numpy(),
        validation["ci_low"].to_numpy(),
        validation["ci_high"].to_numpy(),
        COLORS["blue"],
        f"Validation · {n_validation} donors",
        False,
        -0.11,
    )
    draw_interval_points(
        ax_left,
        y,
        test["point"].to_numpy(),
        test["ci_low"].to_numpy(),
        test["ci_high"].to_numpy(),
        COLORS["orange"],
        f"Held-out test · {n_test} donors",
        True,
        0.11,
    )
    ax_left.axvline(0.0, color="#8796A8", linewidth=1.5)
    ax_left.set_xlim(-1.03, 1.03)
    ax_left.set_ylim(len(AXES) - 0.45, -0.55)
    ax_left.set_yticks(y)
    ax_left.set_yticklabels(AXES, fontsize=12.5, fontweight=800)
    for label, name in zip(ax_left.get_yticklabels(), AXES, strict=True):
        label.set_color(AXIS_COLORS[name])
    ax_left.set_xticks(np.linspace(-1, 1, 5))
    ax_left.grid(axis="x", color=COLORS["grid"], linewidth=1.0)
    ax_left.set_axisbelow(True)
    ax_left.set_xlabel("모델 common 방향 ↔ 관측 partial module slope의 Spearman ρ", fontsize=11.5, color=COLORS["muted"], labelpad=10)
    ax_left.set_title("A. common-only 지문과 관측 변화의 직접 일치도", loc="left", fontsize=16, fontweight=900, color=COLORS["dark"], pad=15)
    ax_left.text(
        0.0,
        1.01,
        "점 = 24개 세포타입 중앙값 · 선 = donor-bootstrap 95% CI",
        transform=ax_left.transAxes,
        fontsize=10.5,
        color=COLORS["muted"],
        va="bottom",
    )
    handles, _ = ax_left.get_legend_handles_labels()
    ax_left.legend(
        handles,
        [f"Validation ({n_validation})", f"Test ({n_test}, exploratory)"],
        frameon=False,
        loc="upper center",
        fontsize=9.8,
        ncol=2,
        bbox_to_anchor=(0.5, -0.095),
    )

    val_half = summary.loc[summary["metric"].eq("observed_cellhalf_validation")].set_index("axis").loc[list(AXES)]
    test_half = summary.loc[summary["metric"].eq("observed_cellhalf_test")].set_index("axis").loc[list(AXES)]
    cross = summary.loc[summary["metric"].eq("observed_validation_test")].set_index("axis").loc[list(AXES)]
    draw_interval_points(
        ax_right,
        y,
        val_half["point"].to_numpy(),
        val_half["ci_low"].to_numpy(),
        val_half["ci_high"].to_numpy(),
        COLORS["teal"],
        "같은 validation donor · 세포 절반 A↔B",
        False,
        -0.16,
    )
    draw_interval_points(
        ax_right,
        y,
        test_half["point"].to_numpy(),
        test_half["ci_low"].to_numpy(),
        test_half["ci_high"].to_numpy(),
        COLORS["navy"],
        "같은 test donor · 세포 절반 A↔B",
        True,
        0.0,
    )
    draw_interval_points(
        ax_right,
        y,
        cross["point"].to_numpy(),
        cross["ci_low"].to_numpy(),
        cross["ci_high"].to_numpy(),
        COLORS["red"],
        "Validation donor ↔ test donor",
        True,
        0.16,
    )
    ax_right.axvline(0.0, color="#8796A8", linewidth=1.5)
    ax_right.set_xlim(-1.03, 1.03)
    ax_right.set_ylim(len(AXES) - 0.45, -0.55)
    ax_right.set_yticks(y)
    ax_right.set_yticklabels([])
    ax_right.set_xticks(np.linspace(-1, 1, 5))
    ax_right.grid(axis="x", color=COLORS["grid"], linewidth=1.0)
    ax_right.set_axisbelow(True)
    ax_right.set_xlabel("관측 partial module slope끼리의 Spearman ρ", fontsize=11.5, color=COLORS["muted"], labelpad=10)
    ax_right.set_title("B. 관측 병리축 방향 자체의 재현성", loc="left", fontsize=16, fontweight=900, color=COLORS["dark"], pad=15)
    ax_right.text(
        0.0,
        1.01,
        "같은 donor의 독립 세포 절반은 안정적이지만 donor split이 바뀌면 방향이 반전",
        transform=ax_right.transAxes,
        fontsize=10.5,
        color=COLORS["muted"],
        va="bottom",
    )
    handles, _ = ax_right.get_legend_handles_labels()
    ax_right.legend(
        handles,
        ["Validation cell-half", "Test cell-half", "Validation ↔ test"],
        frameon=False,
        loc="upper center",
        fontsize=8.7,
        ncol=3,
        bbox_to_anchor=(0.5, -0.095),
        columnspacing=1.0,
        handletextpad=0.35,
    )

    for axis in (ax_left, ax_right):
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.tick_params(axis="x", colors=COLORS["muted"], labelsize=10.5, length=0)
        axis.tick_params(axis="y", length=0)

    warning = FancyBboxPatch(
        (0.055, 0.075),
        0.89,
        0.085,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        transform=figure.transFigure,
        facecolor="#FFF7EA",
        edgecolor="#F0C681",
        linewidth=1.2,
    )
    figure.add_artist(warning)
    figure.text(0.073, 0.119, "!", fontsize=22, fontweight=900, color=COLORS["orange"], ha="center", va="center")
    figure.text(
        0.095,
        0.119,
        (
            "해석 | 세포 수 부족이 아니라 donor 수준의 축 분리 문제이다. "
            f"병리축 최대 상관은 validation {max_corr_validation:.2f}, test {max_corr_test:.2f}였고, "
            "test는 9명뿐이므로 개별 common 축의 관측 검증을 주장할 수 없다."
        ),
        fontsize=12.0,
        color=COLORS["dark"],
        va="center",
    )
    figure.text(
        0.055,
        0.025,
        (
            "Primary: age-standardized, region/technology/known-sex balanced common-only expected-tier endpoint; "
            "observed donor×cell-type module pseudobulk; five pathologies fit simultaneously. "
            "Test per-axis analysis was not a pre-specified confirmatory endpoint."
        ),
        fontsize=8.5,
        color=COLORS["muted"],
    )

    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        kwargs = {"dpi": 220} if extension == "png" else {}
        figure.savefig(output / f"figure_04_common_axis_observed_validation.{extension}", bbox_inches="tight", **kwargs)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--common", type=Path, default=DEFAULT_COMMON)
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--test", type=Path, default=DEFAULT_TEST)
    parser.add_argument("--activity-weight", type=Path, default=DEFAULT_WEIGHT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()

    for path in (args.geometry, args.common, args.validation, args.test, args.activity_weight):
        if not path.exists():
            raise FileNotFoundError(path)

    geometry = np.load(args.geometry, allow_pickle=True)
    common = np.load(args.common, allow_pickle=True)
    weight_payload = np.load(args.activity_weight, allow_pickle=True)
    validation = np.load(args.validation, allow_pickle=True)
    test = np.load(args.test, allow_pickle=True)

    pathology_names = np.asarray(common["pathology_names"], dtype=str)
    if tuple(pathology_names) != AXES:
        raise ValueError(f"unexpected pathology order: {pathology_names.tolist()}")
    common_celltypes = np.asarray(common["celltype_names"], dtype=str)
    validation_celltypes = np.asarray(validation["celltypes"], dtype=str)
    test_celltypes = np.asarray(test["celltypes"], dtype=str)
    if not np.array_equal(validation_celltypes, test_celltypes):
        raise ValueError("validation/test cell-type orders differ")
    celltype_index = index_exact(common_celltypes, validation_celltypes, "cell type")

    validation_modules = np.asarray(validation["modules"], dtype=str)
    test_modules = np.asarray(test["modules"], dtype=str)
    if not np.array_equal(validation_modules, test_modules):
        raise ValueError("validation/test module orders differ")
    weight_modules = np.asarray(weight_payload["module_names"], dtype=str)
    module_index = index_exact(weight_modules, validation_modules, "module")
    activity_weight = np.asarray(weight_payload["activity_weight"], dtype=np.float64)[module_index]
    if activity_weight.shape[1] != common["normal_expected_tier"].shape[1]:
        raise ValueError("activity-weight gene count differs from common endpoint gene count")

    expected_delta = (
        np.asarray(common["disease_expected_tier"], dtype=np.float64)
        - np.asarray(common["normal_expected_tier"], dtype=np.float64)[:, None, :]
    )[celltype_index]
    model_direction = np.einsum(
        "ckg,mg->ckm", expected_delta, activity_weight, optimize=True
    )

    geometry_donors = np.asarray(geometry["donor_names"], dtype=str)
    geometry_pathology_names = np.asarray(geometry["pathology_names"], dtype=str)
    pathology_index = index_exact(geometry_pathology_names, pathology_names, "pathology")

    def split_pathology(payload: np.lib.npyio.NpzFile) -> np.ndarray:
        donors = np.asarray(payload["donors"], dtype=str)
        donor_index = index_exact(geometry_donors, donors, "donor")
        valid = np.asarray(geometry["pathology_valid"], dtype=bool)[donor_index][:, pathology_index]
        if not np.all(valid):
            raise ValueError("a selected donor has a missing named-pathology value")
        return np.asarray(geometry["pathology"], dtype=np.float64)[donor_index][:, pathology_index]

    val_pathology = split_pathology(validation)
    test_pathology = split_pathology(test)
    val_observed = np.asarray(validation["obs_fp"], dtype=np.float64)
    test_observed = np.asarray(test["obs_fp"], dtype=np.float64)
    val_slope = fit_partial_slopes(val_observed, val_pathology)
    test_slope = fit_partial_slopes(test_observed, test_pathology)
    if val_slope is None or test_slope is None:
        raise RuntimeError("primary partial slopes could not be fit")
    model_val = spearman_rows(model_direction, val_slope)
    model_test = spearman_rows(model_direction, test_slope)
    model_val_sign = sign_agreement_rows(model_direction, val_slope)
    model_test_sign = sign_agreement_rows(model_direction, test_slope)

    val_a = fit_partial_slopes(np.asarray(validation["obs_fp_A"], dtype=np.float64), val_pathology)
    val_b = fit_partial_slopes(np.asarray(validation["obs_fp_B"], dtype=np.float64), val_pathology)
    test_a = fit_partial_slopes(np.asarray(test["obs_fp_A"], dtype=np.float64), test_pathology)
    test_b = fit_partial_slopes(np.asarray(test["obs_fp_B"], dtype=np.float64), test_pathology)
    if any(value is None for value in (val_a, val_b, test_a, test_b)):
        raise RuntimeError("cell-half slopes could not be fit")
    val_half = spearman_rows(val_a, val_b)
    test_half = spearman_rows(test_a, test_b)
    val_test = spearman_rows(val_slope, test_slope)

    seed_sequence = np.random.SeedSequence(args.seed)
    child_seeds = seed_sequence.spawn(5)
    boot_model_val = bootstrap_model_observed(
        val_observed, val_pathology, model_direction, args.bootstraps, np.random.default_rng(child_seeds[0])
    )
    boot_model_test = bootstrap_model_observed(
        test_observed, test_pathology, model_direction, args.bootstraps, np.random.default_rng(child_seeds[1])
    )
    boot_val_half = bootstrap_cell_half(
        np.asarray(validation["obs_fp_A"], dtype=np.float64),
        np.asarray(validation["obs_fp_B"], dtype=np.float64),
        val_pathology,
        args.bootstraps,
        np.random.default_rng(child_seeds[2]),
    )
    boot_test_half = bootstrap_cell_half(
        np.asarray(test["obs_fp_A"], dtype=np.float64),
        np.asarray(test["obs_fp_B"], dtype=np.float64),
        test_pathology,
        args.bootstraps,
        np.random.default_rng(child_seeds[3]),
    )
    boot_val_test = bootstrap_cross_split(
        val_observed,
        val_pathology,
        test_observed,
        test_pathology,
        args.bootstraps,
        np.random.default_rng(child_seeds[4]),
    )

    metrics = {
        "model_observed_validation": (axis_median(model_val), boot_model_val),
        "model_observed_test": (axis_median(model_test), boot_model_test),
        "observed_cellhalf_validation": (axis_median(val_half), boot_val_half),
        "observed_cellhalf_test": (axis_median(test_half), boot_test_half),
        "observed_validation_test": (axis_median(val_test), boot_val_test),
    }
    summary_rows: list[dict[str, object]] = []
    for metric, (point, bootstrap) in metrics.items():
        low, high = interval(bootstrap)
        for axis_index, axis_name in enumerate(AXES):
            summary_rows.append(
                {
                    "metric": metric,
                    "axis": axis_name,
                    "point": float(point[axis_index]),
                    "ci_low": float(low[axis_index]),
                    "ci_high": float(high[axis_index]),
                    "n_bootstrap": int(args.bootstraps),
                }
            )
    summary = pd.DataFrame(summary_rows)

    per_celltype_rows: list[dict[str, object]] = []
    for celltype_index_value, celltype_name in enumerate(validation_celltypes):
        for axis_index, axis_name in enumerate(AXES):
            per_celltype_rows.append(
                {
                    "celltype": celltype_name,
                    "axis": axis_name,
                    "model_observed_validation_rho": float(model_val[celltype_index_value, axis_index]),
                    "model_observed_test_rho": float(model_test[celltype_index_value, axis_index]),
                    "model_observed_validation_sign_agreement": float(model_val_sign[celltype_index_value, axis_index]),
                    "model_observed_test_sign_agreement": float(model_test_sign[celltype_index_value, axis_index]),
                    "observed_cellhalf_validation_rho": float(val_half[celltype_index_value, axis_index]),
                    "observed_cellhalf_test_rho": float(test_half[celltype_index_value, axis_index]),
                    "observed_validation_test_rho": float(val_test[celltype_index_value, axis_index]),
                }
            )
    per_celltype = pd.DataFrame(per_celltype_rows)

    # Sensitivity analyses are descriptive and never used to choose the primary method.
    sensitivity_rows: list[dict[str, object]] = []
    for split_name, payload, pathology, observed in (
        ("validation", validation, val_pathology, val_observed),
        ("test", test, test_pathology, test_observed),
    ):
        covariates = np.column_stack(
            [
                np.asarray(payload["sex"], dtype=np.float64),
                np.asarray(payload["age"], dtype=np.float64),
                np.log1p(np.asarray(payload["n_cells"], dtype=np.float64)),
                np.log1p(np.asarray(payload["mean_depth"], dtype=np.float64)),
                np.asarray(payload["dom_batch_frac"], dtype=np.float64),
            ]
        )
        specifications = (
            ("five-axis OLS", None, 0.0),
            ("five-axis ridge1", None, 1.0),
            ("five-axis + sex + age ridge1", covariates[:, :2], 1.0),
            ("five-axis + sex + age + technical ridge1", covariates, 1.0),
        )
        for specification, extra, ridge in specifications:
            slope = fit_partial_slopes(observed, pathology, extra_covariates=extra, ridge=ridge)
            if slope is None:
                continue
            rho = axis_median(spearman_rows(model_direction, slope))
            for axis_index, axis_name in enumerate(AXES):
                sensitivity_rows.append(
                    {
                        "split": split_name,
                        "specification": specification,
                        "axis": axis_name,
                        "median_rho": float(rho[axis_index]),
                    }
                )
    sensitivity = pd.DataFrame(sensitivity_rows)

    def maximum_axis_correlation(pathology: np.ndarray) -> float:
        correlation = np.corrcoef(pathology, rowvar=False)
        upper = np.abs(correlation[np.triu_indices_from(correlation, k=1)])
        return float(np.max(upper))

    max_corr_validation = maximum_axis_correlation(val_pathology)
    max_corr_test = maximum_axis_correlation(test_pathology)
    draw_figure(
        summary,
        args.output,
        n_validation=len(val_pathology),
        n_test=len(test_pathology),
        n_modules=len(validation_modules),
        max_corr_validation=max_corr_validation,
        max_corr_test=max_corr_test,
    )

    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / "figure_04_axis_summary.tsv", sep="\t", index=False)
    per_celltype.to_csv(args.output / "figure_04_per_celltype.tsv", sep="\t", index=False)
    sensitivity.to_csv(args.output / "figure_04_sensitivity.tsv", sep="\t", index=False)
    np.savez_compressed(
        args.output / "figure_04_exact_arrays.npz",
        axes=np.asarray(AXES, dtype=object),
        celltypes=validation_celltypes.astype(object),
        modules=validation_modules.astype(object),
        model_common_module_direction=model_direction.astype(np.float32),
        observed_validation_partial_slope=val_slope.astype(np.float32),
        observed_test_partial_slope=test_slope.astype(np.float32),
        model_observed_validation_rho=model_val.astype(np.float32),
        model_observed_test_rho=model_test.astype(np.float32),
        observed_cellhalf_validation_rho=val_half.astype(np.float32),
        observed_cellhalf_test_rho=test_half.astype(np.float32),
        observed_validation_test_rho=val_test.astype(np.float32),
        bootstrap_model_validation=boot_model_val.astype(np.float32),
        bootstrap_model_test=boot_model_test.astype(np.float32),
        bootstrap_cellhalf_validation=boot_val_half.astype(np.float32),
        bootstrap_cellhalf_test=boot_test_half.astype(np.float32),
        bootstrap_validation_test=boot_val_test.astype(np.float32),
    )

    manifest = {
        "schema_version": "prism.presentation.figure04.v1",
        "model": "PRISM personal-rank2 epoch 20",
        "analysis": "frozen common-only named-pathology direction versus observed module partial slope",
        "primary_model_direction": (
            "region/technology/known-sex balanced, age-standardized common-only expected-tier endpoint "
            "difference projected to the exact tokenizer module coordinate"
        ),
        "primary_observed_direction": (
            "donor-by-cell-type observed module pseudobulk; five named pathologies standardized and fit simultaneously by OLS"
        ),
        "comparison": "Spearman across 404 modules within cell type; median across 24 cell types",
        "uncertainty": f"{args.bootstraps} donor-bootstrap resamples",
        "validation_donors": int(len(val_pathology)),
        "test_donors": int(len(test_pathology)),
        "test_status": (
            "exploratory for the per-axis endpoint; original locked-test plan excluded test-only per-axis donor-disjoint inference below 12 donors"
        ),
        "maximum_absolute_pathology_correlation": {
            "validation": max_corr_validation,
            "test": max_corr_test,
        },
        "sources": {
            "geometry": {"path": str(args.geometry.resolve()), "sha256": sha256(args.geometry)},
            "common": {"path": str(args.common.resolve()), "sha256": sha256(args.common)},
            "validation": {"path": str(args.validation.resolve()), "sha256": sha256(args.validation)},
            "test": {"path": str(args.test.resolve()), "sha256": sha256(args.test)},
            "activity_weight": {"path": str(args.activity_weight.resolve()), "sha256": sha256(args.activity_weight)},
        },
        "outputs": {
            "figure_png": "figure_04_common_axis_observed_validation.png",
            "figure_pdf": "figure_04_common_axis_observed_validation.pdf",
            "figure_svg": "figure_04_common_axis_observed_validation.svg",
            "summary": "figure_04_axis_summary.tsv",
            "per_celltype": "figure_04_per_celltype.tsv",
            "sensitivity": "figure_04_sensitivity.tsv",
            "exact_arrays": "figure_04_exact_arrays.npz",
        },
    }
    (args.output / "figure_04_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
