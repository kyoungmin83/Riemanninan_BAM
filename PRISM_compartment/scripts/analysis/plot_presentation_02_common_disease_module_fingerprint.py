#!/usr/bin/env python3
"""Render a readable five-axis common-disease module fingerprint atlas."""

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


PRISM_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    PRISM_ROOT.parent
    / "kmlee_bam"
    / "analysis_outputs"
    / "prism_rank2_figure11_13_3d_20260811"
    / "geometry"
    / "prism_fingerprints_exact.npz"
)
DEFAULT_VALIDATION_RECOVERY = (
    PRISM_ROOT.parent
    / "kmlee_bam"
    / "analysis_outputs"
    / "ct64_rank2_rank4_k2_comparison_20260811"
    / "celltype_module_recovery_long.tsv"
)
DEFAULT_TEST_RECOVERY = (
    PRISM_ROOT
    / "analysis_outputs"
    / "presentation_results_20260819"
    / "00_locked_test_rank2_e20"
    / "rank2_e20_module_disease_test.json"
)
DEFAULT_OUTPUT = (
    PRISM_ROOT
    / "analysis_outputs"
    / "presentation_results_20260819"
    / "02_common_disease_module_fingerprint"
)

CELLTYPE_ORDER = [
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
    "Chandelier",
    "L2/3 IT",
    "L4 IT",
    "L5 ET",
    "L5 IT",
    "L5/6 NP",
    "L6 CT",
    "L6 IT",
    "L6 IT Car3",
    "L6b",
    "Lamp5",
    "Lamp5 Lhx6",
    "Pax6",
    "Pvalb",
    "Sncg",
    "Sst",
    "Sst Chodl",
    "Vip",
]

AXIS_COLORS = {
    "Thal": "#168A83",
    "Braak": "#E67E32",
    "CERAD": "#D6A119",
    "LATE": "#7551A8",
    "Lewy": "#B23A7A",
}

COLORS = {
    "navy": "#123B73",
    "dark": "#172033",
    "muted": "#667085",
    "green": "#2E8B57",
    "grid": "#DCE3EC",
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
    font_family = (
        font_manager.FontProperties(fname=font_paths[0]).get_name()
        if font_paths[0].exists()
        else "DejaVu Sans"
    )
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font_family, "DejaVu Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def short_celltype(name: str) -> str:
    return name.replace("Microglia-PVM", "Microglia").replace("Oligodendrocyte", "Oligodendro.")


def short_module_name(name: str) -> str:
    if name.startswith("HALLMARK::HALLMARK_"):
        pathway = name.split("HALLMARK::HALLMARK_", 1)[1]
        aliases = {
            "OXIDATIVE_PHOSPHORYLATION": "OxPhos",
            "CHOLESTEROL_HOMEOSTASIS": "Cholesterol",
            "DNA_REPAIR": "DNA repair",
            "MYC_TARGETS_V1": "MYC targets",
        }
        return f"Hallmark\n{aliases.get(pathway, pathway.replace('_', ' ').title())}"
    if name.startswith("hdWGCNA_"):
        prefix, celltype, color = name.split("::", 2)
        region = prefix.replace("hdWGCNA_", "")
        clean_celltype = celltype.replace("_", " ").replace("Microglia-PVM", "Microglia")
        clean_celltype = clean_celltype.replace("Oligodendrocyte", "Oligo.")
        return f"{region} {clean_celltype}\n{color}"
    return name.replace("::", "\n")[:32]


def load_common(source: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(source, allow_pickle=True) as payload:
        common = payload["common_high"].astype(np.float64)
        celltypes = payload["celltype_names"].astype(str)
        axes = payload["pathology_names"].astype(str)
        modules = payload["module_names"].astype(str)

    if common.shape != (24, 2, 5, 414):
        raise ValueError(f"Unexpected common_high shape: {common.shape}")
    if set(celltypes) != set(CELLTYPE_ORDER):
        raise ValueError("Cell-type vocabulary does not match the frozen presentation order")
    order = np.asarray([int(np.where(celltypes == name)[0][0]) for name in CELLTYPE_ORDER])
    common_region_mean = common[order].mean(axis=1)
    return common_region_mean, celltypes[order], axes, modules


def load_recovery_summary(validation_path: Path, test_path: Path) -> dict[str, float | int]:
    validation = pd.read_csv(validation_path, sep="\t")
    validation = validation.loc[validation["model"].eq("rank2_e20")].copy()
    if set(validation["celltype"]) != set(CELLTYPE_ORDER):
        raise ValueError("Validation recovery table does not contain the expected 24 cell types")
    with test_path.open(encoding="utf-8") as handle:
        test_json = json.load(handle)
    test_values = np.asarray(
        [test_json["per_celltype"][name]["AD_module_spearman"] for name in CELLTYPE_ORDER],
        dtype=float,
    )
    validation_values = validation["AD_module"].to_numpy(float)
    return {
        "validation_median_AD_module": float(np.nanmedian(validation_values)),
        "test_median_AD_module": float(np.nanmedian(test_values)),
        "validation_positive_celltypes": int(np.sum(validation_values > 0)),
        "test_positive_celltypes": int(np.sum(test_values > 0)),
        "n_celltypes": int(len(CELLTYPE_ORDER)),
        "test_donors": int(test_json.get("module_to_ADNC", {}).get("n_donors", len(test_json.get("per_donor", [])))),
    }


def draw_figure(
    common: np.ndarray,
    celltypes: np.ndarray,
    axes_names: np.ndarray,
    modules: np.ndarray,
    output: Path,
    top_per_axis: int,
    recovery: dict[str, float | int],
) -> tuple[pd.DataFrame, float]:
    configure_style()

    selected: list[np.ndarray] = []
    records: list[dict[str, object]] = []
    for axis_index, axis_name in enumerate(axes_names):
        rms = np.sqrt(np.mean(np.square(common[:, axis_index, :]), axis=0))
        indices = np.argsort(-rms, kind="stable")[:top_per_axis]
        selected.append(indices)
        for rank, module_index in enumerate(indices, start=1):
            records.append(
                {
                    "pathology": str(axis_name),
                    "rank_within_pathology": rank,
                    "module_index": int(module_index),
                    "module": str(modules[module_index]),
                    "cross_celltype_rms": float(rms[module_index]),
                }
            )

    shown_values = np.concatenate(
        [common[:, axis_index, indices].ravel() for axis_index, indices in enumerate(selected)]
    )
    color_limit = float(np.quantile(np.abs(shown_values[np.isfinite(shown_values)]), 0.98))
    color_limit = max(color_limit, 1e-8)

    figure = plt.figure(figsize=(16, 9), facecolor="white")
    grid = figure.add_gridspec(
        nrows=1,
        ncols=5,
        left=0.165,
        right=0.935,
        bottom=0.205,
        top=0.715,
        wspace=0.10,
    )
    image = None
    for axis_index, (axis_name, indices) in enumerate(zip(axes_names, selected, strict=True)):
        axis = figure.add_subplot(grid[0, axis_index])
        matrix = common[:, axis_index, indices]
        x = np.arange(top_per_axis + 1)
        y = np.arange(len(celltypes) + 1)
        image = axis.pcolormesh(
            x,
            y,
            matrix,
            cmap="RdBu_r",
            vmin=-color_limit,
            vmax=color_limit,
            shading="flat",
            edgecolors="white",
            linewidth=0.65,
            rasterized=False,
        )
        axis.set_xlim(0, top_per_axis)
        axis.set_ylim(len(celltypes), 0)
        axis.set_xticks(np.arange(top_per_axis) + 0.5)
        axis.set_xticklabels(
            [short_module_name(str(modules[index])) for index in indices],
            rotation=58,
            ha="right",
            va="top",
            rotation_mode="anchor",
            fontsize=8.4,
            color=COLORS["dark"],
        )
        axis.tick_params(axis="x", length=0, pad=5)
        axis.set_yticks(np.arange(len(celltypes)) + 0.5)
        if axis_index == 0:
            axis.set_yticklabels([short_celltype(str(name)) for name in celltypes], fontsize=10.2)
            for index, label in enumerate(axis.get_yticklabels()):
                label.set_color(COLORS["green"] if index < 6 else COLORS["dark"])
                label.set_fontweight(700 if index < 6 else 500)
        else:
            axis.set_yticklabels([])
        axis.tick_params(axis="y", length=0, pad=6)
        axis.axhline(6, color="#8496AA", linewidth=1.7)
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.set_title(
            str(axis_name),
            fontsize=15,
            color="white",
            fontweight=800,
            pad=11,
            bbox={
                "boxstyle": "round,pad=0.45,rounding_size=0.18",
                "facecolor": AXIS_COLORS[str(axis_name)],
                "edgecolor": AXIS_COLORS[str(axis_name)],
                "linewidth": 0.0,
            },
        )

    assert image is not None
    color_axis = figure.add_axes([0.950, 0.275, 0.014, 0.365])
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("Common module coefficient", fontsize=10.5, color=COLORS["muted"], labelpad=7)
    colorbar.ax.tick_params(labelsize=9.5, colors=COLORS["muted"], length=3)
    colorbar.outline.set_edgecolor("#AAB5C2")

    figure.text(
        0.045,
        0.955,
        "PRISM의 공통 병리효과는 세포타입마다 서로 다른 모듈 조합으로 나타난다",
        ha="left",
        va="top",
        fontsize=25,
        fontweight=800,
        color=COLORS["navy"],
    )
    figure.text(
        0.047,
        0.907,
        "PRISM personal-rank2 · epoch 20  |  Personal code = 0 · 개인×병리반응 = 0 · DLPFC/MTG 평균",
        ha="left",
        va="top",
        fontsize=13.2,
        color=COLORS["muted"],
    )
    figure.text(
        0.957,
        0.912,
        "FROZEN MODEL · MODEL-DERIVED POSTHOC",
        ha="right",
        va="center",
        fontsize=10.2,
        color=COLORS["navy"],
        fontweight=800,
        bbox={
            "boxstyle": "round,pad=0.45,rounding_size=0.2",
            "facecolor": "#EDF3FB",
            "edgecolor": "#B9CBE2",
            "linewidth": 1.0,
        },
    )
    figure.text(
        0.165,
        0.790,
        f"각 병리축에서 세포타입 간 RMS가 큰 상위 {top_per_axis}개 모듈",
        ha="left",
        va="center",
        fontsize=13.5,
        color=COLORS["dark"],
        fontweight=700,
    )
    figure.text(
        0.932,
        0.790,
        "관측 복원 확인  |  validation ρ = "
        f"{recovery['validation_median_AD_module']:.3f}  →  locked test ρ = {recovery['test_median_AD_module']:.3f}\n"
        f"Test: {recovery['test_positive_celltypes']}/{recovery['n_celltypes']} 세포타입에서 양의 AD 모듈 방향",
        ha="right",
        va="center",
        fontsize=10.2,
        color="#C7621F",
        fontweight=700,
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.42,rounding_size=0.18",
            "facecolor": "#FFF8F1",
            "edgecolor": "#E67E32",
            "linewidth": 1.1,
        },
    )
    figure.text(
        0.058,
        0.645,
        "비신경\n6종",
        ha="center",
        va="center",
        fontsize=11.5,
        color=COLORS["green"],
        fontweight=800,
    )
    figure.text(
        0.058,
        0.410,
        "신경\n18종",
        ha="center",
        va="center",
        fontsize=11.5,
        color=COLORS["navy"],
        fontweight=800,
    )

    banner = FancyBboxPatch(
        (0.040, 0.060),
        0.920,
        0.055,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=figure.transFigure,
        facecolor=COLORS["navy"],
        edgecolor=COLORS["navy"],
        linewidth=0.0,
        zorder=2,
    )
    figure.add_artist(banner)
    figure.text(
        0.500,
        0.087,
        "핵심  |  하나의 공통 질병효과를 모든 세포에 복사하지 않고, 병리축과 세포타입에 따라 모듈 방향과 크기를 다르게 추정한다.",
        ha="center",
        va="center",
        fontsize=12.4,
        color="white",
        fontweight=700,
        zorder=3,
    )
    figure.text(
        0.047,
        0.024,
        "빨강 = 해당 모듈의 양의 common coefficient, 파랑 = 음의 coefficient. 모듈은 축별 cross-celltype RMS로 선택했다. "
        "이는 frozen decoder의 모델 내부 지문이며 관측 differential expression 또는 인과효과 자체를 뜻하지 않는다.",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color=COLORS["muted"],
    )

    output.mkdir(parents=True, exist_ok=True)
    base = output / "figure_02_common_disease_module_fingerprint"
    figure.savefig(base.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(base.with_suffix(".pdf"), facecolor="white")
    figure.savefig(base.with_suffix(".svg"), facecolor="white")
    plt.close(figure)
    return pd.DataFrame.from_records(records), color_limit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--validation-recovery", type=Path, default=DEFAULT_VALIDATION_RECOVERY)
    parser.add_argument("--test-recovery", type=Path, default=DEFAULT_TEST_RECOVERY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-per-axis", type=int, default=4)
    args = parser.parse_args()

    source = args.source.resolve()
    validation_recovery = args.validation_recovery.resolve()
    test_recovery = args.test_recovery.resolve()
    output = args.output.resolve()
    common, celltypes, axes_names, modules = load_common(source)
    recovery = load_recovery_summary(validation_recovery, test_recovery)
    selected, color_limit = draw_figure(
        common, celltypes, axes_names, modules, output, args.top_per_axis, recovery
    )
    selected_path = output / "figure_02_selected_modules.tsv"
    selected.to_csv(selected_path, sep="\t", index=False)
    manifest = {
        "schema_version": "prism.presentation.figure02.v1",
        "model": "PRISM personal-rank2 epoch 20",
        "analysis": "common-only disease module coefficients",
        "personal_code": 0,
        "personal_response": 0,
        "region_summary": "mean of DLPFC and MTG common coefficients",
        "source": str(source),
        "source_sha256": sha256(source),
        "common_array": "common_high",
        "common_array_shape": [24, 2, 5, 414],
        "top_modules_per_axis": int(args.top_per_axis),
        "module_selection": "largest RMS across 24 cell types within each pathology axis",
        "shared_color_limit_abs_q98": color_limit,
        "recovery_validation_source": str(validation_recovery),
        "recovery_validation_source_sha256": sha256(validation_recovery),
        "recovery_test_source": str(test_recovery),
        "recovery_test_source_sha256": sha256(test_recovery),
        "recovery_summary": recovery,
        "test_data_used": "summary badge only; coefficient heatmap is split-invariant",
        "outputs": {
            "png": "figure_02_common_disease_module_fingerprint.png",
            "pdf": "figure_02_common_disease_module_fingerprint.pdf",
            "svg": "figure_02_common_disease_module_fingerprint.svg",
            "selected_modules": selected_path.name,
        },
    }
    (output / "figure_02_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
