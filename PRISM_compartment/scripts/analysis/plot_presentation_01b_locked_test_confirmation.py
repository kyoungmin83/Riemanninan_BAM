#!/usr/bin/env python3
"""Compare validation and locked-test module recovery for frozen rank2 epoch 20."""

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
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


PRISM_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VALIDATION = (
    PRISM_ROOT.parent
    / "kmlee_bam"
    / "analysis_outputs"
    / "ct64_rank2_rank4_k2_comparison_20260811"
    / "celltype_module_recovery_long.tsv"
)
DEFAULT_TEST = (
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
    / "01b_locked_test_confirmation"
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

COLORS = {
    "navy": "#123B73",
    "teal": "#168A83",
    "orange": "#E67E32",
    "green": "#2E8B57",
    "dark": "#172033",
    "muted": "#667085",
    "grid": "#D9E1EA",
    "connector": "#B8C2CE",
    "neuron_bg": "#F4F7FC",
    "non_neuron_bg": "#F1F9F5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    regular = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    bold = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc")
    for font_path in (regular, bold):
        if font_path.exists():
            font_manager.fontManager.addfont(font_path)
    family = (
        font_manager.FontProperties(fname=regular).get_name()
        if regular.exists()
        else "DejaVu Sans"
    )
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [family, "DejaVu Sans"],
            "axes.unicode_minus": False,
            "axes.titlesize": 17,
            "axes.titleweight": 700,
            "axes.labelsize": 13,
            "xtick.labelsize": 11.5,
            "ytick.labelsize": 12.0,
            "legend.fontsize": 11.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_data(validation_path: Path, test_path: Path) -> tuple[pd.DataFrame, dict]:
    validation_all = pd.read_csv(validation_path, sep="\t")
    validation = validation_all.loc[validation_all["model"].eq("rank2_e20")].copy()
    if set(validation["celltype"]) != set(CELLTYPE_ORDER):
        raise ValueError("Validation table does not contain the expected 24 cell types")

    with test_path.open(encoding="utf-8") as handle:
        test_json = json.load(handle)
    test_per_celltype = test_json.get("per_celltype", {})
    if set(test_per_celltype) != set(CELLTYPE_ORDER):
        absent = sorted(set(CELLTYPE_ORDER).difference(test_per_celltype))
        extra = sorted(set(test_per_celltype).difference(CELLTYPE_ORDER))
        raise ValueError(f"Unexpected test cell types; absent={absent}, extra={extra}")

    validation = validation.set_index("celltype")
    rows = []
    for celltype in CELLTYPE_ORDER:
        val = validation.loc[celltype]
        test = test_per_celltype[celltype]
        rows.append(
            {
                "celltype": celltype,
                "class": val["class"],
                "validation_n_donors": int(val["n_donors"]),
                "test_n_donors": int(test["n_donors"]),
                "validation_blind_centered": float(val["blind_centered"]),
                "test_blind_centered": float(test["blind_centered"]),
                "validation_AD_module": float(val["AD_module"]),
                "test_AD_module": float(test["AD_module_spearman"]),
            }
        )
    return pd.DataFrame(rows), test_json


def style_axis(axis: plt.Axes, n_rows: int, non_neuronal_count: int) -> None:
    axis.set_xlim(-0.1, 1.0)
    axis.set_xticks(np.arange(0.0, 1.01, 0.2))
    axis.set_ylim(n_rows - 0.5, -0.5)
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.9, zorder=0)
    axis.axvline(0.0, color="#91A0B2", linewidth=1.2, zorder=1)
    axis.axhspan(-0.5, non_neuronal_count - 0.5, color=COLORS["non_neuron_bg"], zorder=-3)
    axis.axhspan(non_neuronal_count - 0.5, n_rows - 0.5, color=COLORS["neuron_bg"], zorder=-3)
    axis.axhline(non_neuronal_count - 0.5, color="#B7C4D2", linewidth=1.2, zorder=1)
    axis.tick_params(axis="x", colors=COLORS["muted"], length=0, pad=7)
    axis.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#AAB5C2")


def plot_pair(
    axis: plt.Axes,
    y: np.ndarray,
    validation: np.ndarray,
    test: np.ndarray,
    color: str,
) -> None:
    for row, val_value, test_value in zip(y, validation, test, strict=True):
        axis.plot(
            [val_value, test_value],
            [row, row],
            color=COLORS["connector"],
            linewidth=2.0,
            solid_capstyle="round",
            zorder=2,
        )
    axis.scatter(
        validation,
        y,
        s=72,
        facecolor="white",
        edgecolor=color,
        linewidth=1.8,
        zorder=4,
    )
    axis.scatter(
        test,
        y,
        s=82,
        facecolor=color,
        edgecolor="white",
        linewidth=1.1,
        zorder=5,
    )


def add_metric_chip(
    figure: plt.Figure,
    x: float,
    y: float,
    label: str,
    validation_median: float,
    test_median: float,
    color: str,
) -> None:
    figure.text(
        x,
        y,
        f"{label}  validation {validation_median:.3f}  →  test {test_median:.3f}",
        ha="left",
        va="center",
        fontsize=11.4,
        color=color,
        fontweight=700,
        bbox={
            "boxstyle": "round,pad=0.42,rounding_size=0.18",
            "facecolor": "white",
            "edgecolor": color,
            "linewidth": 1.2,
        },
    )


def draw_figure(data: pd.DataFrame, test_json: dict, output: Path) -> dict:
    configure_style()
    figure = plt.figure(figsize=(16, 9), facecolor="white")
    grid = figure.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.16,
        right=0.975,
        bottom=0.195,
        top=0.725,
        width_ratios=(1.0, 1.0),
        wspace=0.15,
    )
    axis_blind = figure.add_subplot(grid[0, 0])
    axis_ad = figure.add_subplot(grid[0, 1], sharey=axis_blind)

    y = np.arange(len(data))
    non_neuronal_count = int(data["class"].eq("non-neuronal").sum())
    for axis in (axis_blind, axis_ad):
        style_axis(axis, len(data), non_neuronal_count)

    val_blind = data["validation_blind_centered"].to_numpy(float)
    test_blind = data["test_blind_centered"].to_numpy(float)
    val_ad = data["validation_AD_module"].to_numpy(float)
    test_ad = data["test_AD_module"].to_numpy(float)

    plot_pair(axis_blind, y, val_blind, test_blind, COLORS["teal"])
    plot_pair(axis_ad, y, val_ad, test_ad, COLORS["orange"])

    axis_blind.set_yticks(y)
    axis_blind.set_yticklabels(data["celltype"])
    for index, label in enumerate(axis_blind.get_yticklabels()):
        label.set_color(COLORS["green"] if index < non_neuronal_count else COLORS["dark"])
        label.set_fontweight(700 if index < non_neuronal_count else 500)
    axis_ad.tick_params(labelleft=False)
    axis_blind.set_xlabel("관측–예측 일반 모듈 효과의 Spearman ρ", color=COLORS["muted"], labelpad=10)
    axis_ad.set_xlabel("관측–예측 AD 모듈 방향의 Spearman ρ", color=COLORS["muted"], labelpad=10)

    provenance = test_json.get("evaluation_provenance", {})
    selected_cells = int(provenance.get("selected_cells", 0))
    test_donors = int(
        test_json.get("module_to_ADNC", {}).get(
            "n_donors", len(test_json.get("per_donor", []))
        )
    )
    medians = {
        "validation_blind_centered": float(np.median(val_blind)),
        "test_blind_centered": float(np.median(test_blind)),
        "validation_AD_module": float(np.median(val_ad)),
        "test_AD_module": float(np.median(test_ad)),
    }

    figure.text(
        0.045,
        0.968,
        f"Frozen PRISM을 잠가 둔 {test_donors}명 test donor에서 독립 평가했다",
        ha="left",
        va="top",
        fontsize=24,
        fontweight=800,
        color=COLORS["navy"],
    )
    figure.text(
        0.047,
        0.898,
        f"PRISM personal-rank2 · epoch 20  |  validation 16명 → locked test {test_donors}명 · {selected_cells:,}개 test 세포",
        ha="left",
        va="top",
        fontsize=13.2,
        color=COLORS["muted"],
    )
    figure.text(
        0.957,
        0.901,
        "LOCKED TEST · NO RETUNING",
        ha="right",
        va="center",
        fontsize=10.5,
        color=COLORS["navy"],
        fontweight=800,
        bbox={
            "boxstyle": "round,pad=0.45,rounding_size=0.2",
            "facecolor": "#EDF3FB",
            "edgecolor": "#B9CBE2",
            "linewidth": 1.0,
        },
    )

    figure.text(0.160, 0.820, "A. 일반 모듈 모양 복원", ha="left", va="center", fontsize=16, color=COLORS["dark"], fontweight=800)
    figure.text(0.600, 0.820, "B. AD 관련 모듈 방향 복원", ha="left", va="center", fontsize=16, color=COLORS["dark"], fontweight=800)
    add_metric_chip(figure, 0.160, 0.770, "중앙값 ρ", medians["validation_blind_centered"], medians["test_blind_centered"], COLORS["teal"])
    add_metric_chip(figure, 0.600, 0.770, "중앙값 ρ", medians["validation_AD_module"], medians["test_AD_module"], COLORS["orange"])

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=COLORS["navy"], markeredgewidth=1.8, markersize=8, label="Validation"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["navy"], markeredgecolor="white", markeredgewidth=1.0, markersize=9, label="Locked test"),
    ]
    figure.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.750), ncol=2, frameon=False, handletextpad=0.5, columnspacing=1.8)

    figure.text(0.032, 0.660, "비신경\n6종", ha="center", va="center", fontsize=11.5, color=COLORS["green"], fontweight=800)
    figure.text(0.032, 0.390, "신경\n18종", ha="center", va="center", fontsize=11.5, color=COLORS["navy"], fontweight=800)

    banner = FancyBboxPatch(
        (0.040, 0.060),
        0.920,
        0.052,
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
        0.086,
        "핵심  |  test 중앙값은 일반 모듈 ρ = "
        f"{medians['test_blind_centered']:.3f}, AD 모듈 ρ = {medians['test_AD_module']:.3f}였다.",
        ha="center",
        va="center",
        fontsize=12.7,
        color="white",
        fontweight=700,
        zorder=3,
    )
    figure.text(
        0.047,
        0.027,
        f"빈 점 = validation, 채운 점 = 사전에 잠근 {test_donors}명 test donor. 세포타입별 유효 donor는 "
        f"{int(data['test_n_donors'].min())}–{int(data['test_n_donors'].max())}명이다. Test는 checkpoint 선택·재학습·threshold 조정에 사용하지 않았다.\n"
        f"AD-module은 ADNC·LATE·Lewy·성별·나이를 조정한다. Test n = {test_donors}이므로 donor-disjoint 절반분할은 수행하지 않았다. 점은 세포타입별 점추정이다.",
        ha="left",
        va="center",
        fontsize=8.8,
        color=COLORS["muted"],
        linespacing=1.35,
    )

    output.mkdir(parents=True, exist_ok=True)
    base = output / "figure_01b_locked_test_confirmation"
    figure.savefig(base.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(base.with_suffix(".pdf"), facecolor="white")
    figure.savefig(base.with_suffix(".svg"), facecolor="white")
    plt.close(figure)
    return {
        "medians": medians,
        "selected_test_cells": selected_cells,
        "test_donors_min": int(data["test_n_donors"].min()),
        "test_donors_max": int(data["test_n_donors"].max()),
        "test_donors_overall": test_donors,
        "all_test_blind_positive": bool((test_blind > 0).all()),
        "all_test_AD_module_positive": bool((test_ad > 0).all()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--test", type=Path, default=DEFAULT_TEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    validation_path = args.validation.resolve()
    test_path = args.test.resolve()
    output = args.output.resolve()
    data, test_json = load_data(validation_path, test_path)
    summary = draw_figure(data, test_json, output)

    data_path = output / "figure_01b_validation_vs_test_celltype_values.tsv"
    data.to_csv(data_path, sep="\t", index=False)
    manifest = {
        "schema_version": "prism.presentation.figure01b.v1",
        "model": "PRISM personal-rank2 epoch 20",
        "comparison": "validation versus locked test",
        "validation_source": str(validation_path),
        "validation_source_sha256": sha256(validation_path),
        "test_source": str(test_path),
        "test_source_sha256": sha256(test_path),
        "n_celltypes": int(len(data)),
        **summary,
        "outputs": {
            "png": "figure_01b_locked_test_confirmation.png",
            "pdf": "figure_01b_locked_test_confirmation.pdf",
            "svg": "figure_01b_locked_test_confirmation.svg",
            "data": data_path.name,
        },
    }
    (output / "figure_01b_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
