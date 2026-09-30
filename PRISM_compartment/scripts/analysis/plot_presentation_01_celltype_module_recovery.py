#!/usr/bin/env python3
"""Render the slide-ready PRISM rank2 cell-type module-recovery figure."""

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
    / "ct64_rank2_rank4_k2_comparison_20260811"
    / "celltype_module_recovery_long.tsv"
)
DEFAULT_OUTPUT = (
    PRISM_ROOT
    / "analysis_outputs"
    / "presentation_results_20260819"
    / "01_celltype_module_recovery"
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
    "blue": "#1F67B1",
    "teal": "#168A83",
    "orange": "#E67E32",
    "green": "#2E8B57",
    "dark": "#172033",
    "muted": "#667085",
    "grid": "#D9E1EA",
    "neuron_bg": "#F4F7FC",
    "non_neuron_bg": "#F1F9F5",
    "connector": "#B8C2CE",
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
            "axes.titlesize": 17,
            "axes.titleweight": 700,
            "axes.labelsize": 13,
            "xtick.labelsize": 11.5,
            "ytick.labelsize": 12.0,
            "legend.fontsize": 11.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_rank2(source: Path) -> pd.DataFrame:
    data = pd.read_csv(source, sep="\t")
    required = {
        "celltype",
        "class",
        "model",
        "n_donors",
        "blind_centered",
        "AD_module",
        "cross_disjoint",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    rank2 = data.loc[data["model"].eq("rank2_e20")].copy()
    if set(rank2["celltype"]) != set(CELLTYPE_ORDER):
        absent = sorted(set(CELLTYPE_ORDER).difference(rank2["celltype"]))
        extra = sorted(set(rank2["celltype"]).difference(CELLTYPE_ORDER))
        raise ValueError(f"Unexpected cell types; absent={absent}, extra={extra}")

    rank2["celltype"] = pd.Categorical(
        rank2["celltype"], categories=CELLTYPE_ORDER, ordered=True
    )
    rank2 = rank2.sort_values("celltype").reset_index(drop=True)
    return rank2


def style_axis(axis: plt.Axes) -> None:
    axis.set_xlim(0.0, 1.0)
    axis.set_xticks(np.arange(0.0, 1.01, 0.2))
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.9, zorder=0)
    axis.tick_params(axis="x", colors=COLORS["muted"], length=0, pad=7)
    axis.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#AAB5C2")
    axis.spines["bottom"].set_linewidth(1.0)


def shade_cell_classes(axis: plt.Axes, non_neuronal_count: int) -> None:
    axis.axhspan(-0.5, non_neuronal_count - 0.5, color=COLORS["non_neuron_bg"], zorder=-3)
    axis.axhspan(
        non_neuronal_count - 0.5,
        len(CELLTYPE_ORDER) - 0.5,
        color=COLORS["neuron_bg"],
        zorder=-3,
    )
    axis.axhline(
        non_neuronal_count - 0.5,
        color="#B7C4D2",
        linewidth=1.2,
        zorder=1,
    )


def add_metric_chip(
    figure: plt.Figure,
    x: float,
    y: float,
    label: str,
    value: float,
    color: str,
) -> None:
    figure.text(
        x,
        y,
        f"{label}  중앙값 ρ = {value:.3f}",
        ha="left",
        va="center",
        fontsize=11.5,
        color=color,
        fontweight=700,
        bbox={
            "boxstyle": "round,pad=0.42,rounding_size=0.18",
            "facecolor": "white",
            "edgecolor": color,
            "linewidth": 1.2,
        },
    )


def draw_figure(rank2: pd.DataFrame, output: Path) -> dict[str, float]:
    configure_style()

    figure = plt.figure(figsize=(16, 9), facecolor="white")
    grid = figure.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.16,
        right=0.975,
        bottom=0.190,
        top=0.745,
        width_ratios=(0.88, 1.20),
        wspace=0.16,
    )
    axis_blind = figure.add_subplot(grid[0, 0])
    axis_disease = figure.add_subplot(grid[0, 1], sharey=axis_blind)

    y = np.arange(len(rank2))
    non_neuronal_count = int(rank2["class"].eq("non-neuronal").sum())

    for axis in (axis_blind, axis_disease):
        style_axis(axis)
        shade_cell_classes(axis, non_neuronal_count)
        axis.set_ylim(len(rank2) - 0.5, -0.5)

    blind = rank2["blind_centered"].to_numpy(float)
    ad_module = rank2["AD_module"].to_numpy(float)
    cross = rank2["cross_disjoint"].to_numpy(float)

    axis_blind.hlines(
        y,
        0.0,
        blind,
        color="#B8C8D8",
        linewidth=2.6,
        zorder=2,
    )
    axis_blind.scatter(
        blind,
        y,
        s=72,
        color=COLORS["teal"],
        edgecolor="white",
        linewidth=1.2,
        zorder=4,
    )

    for row, ad_value, cross_value in zip(y, ad_module, cross, strict=True):
        axis_disease.plot(
            [cross_value, ad_value],
            [row, row],
            color=COLORS["connector"],
            linewidth=2.1,
            solid_capstyle="round",
            zorder=2,
        )
    axis_disease.scatter(
        ad_module,
        y,
        s=70,
        color=COLORS["orange"],
        edgecolor="white",
        linewidth=1.1,
        label="같은 validation donor 기반 AD 복원",
        zorder=4,
    )
    axis_disease.scatter(
        cross,
        y,
        s=76,
        color=COLORS["navy"],
        edgecolor="white",
        linewidth=1.1,
        label="Donor-disjoint 질병 방향 복원",
        zorder=5,
    )

    axis_blind.set_yticks(y)
    axis_blind.set_yticklabels(rank2["celltype"].astype(str))
    for index, label in enumerate(axis_blind.get_yticklabels()):
        label.set_color(
            COLORS["green"] if index < non_neuronal_count else COLORS["dark"]
        )
        label.set_fontweight(700 if index < non_neuronal_count else 500)
    axis_disease.tick_params(labelleft=False)

    axis_blind.set_xlabel("관측–예측 모듈 효과의 Spearman ρ", color=COLORS["muted"], labelpad=10)
    axis_disease.set_xlabel("관측–예측 질병 방향의 Spearman ρ", color=COLORS["muted"], labelpad=10)

    figure.text(
        0.045,
        0.955,
        "PRISM은 세포타입별 모듈 변화를 복원하고, 질병 방향을 새 donor로 전달한다",
        ha="left",
        va="top",
        fontsize=25,
        fontweight=800,
        color=COLORS["navy"],
    )
    figure.text(
        0.047,
        0.907,
        "PRISM personal-rank2 · epoch 20  |  24개 세포타입 · 404개 평가 모듈 · 13–16명 validation donor",
        ha="left",
        va="top",
        fontsize=13.2,
        color=COLORS["muted"],
    )
    figure.text(
        0.957,
        0.912,
        "VALIDATION ONLY · TEST UNUSED",
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

    medians = {
        "blind_centered": float(np.median(blind)),
        "AD_module": float(np.median(ad_module)),
        "cross_disjoint": float(np.median(cross)),
    }
    figure.text(
        0.160,
        0.840,
        "A. 세포타입별 일반 모듈 모양 복원",
        ha="left",
        va="center",
        fontsize=17,
        color=COLORS["dark"],
        fontweight=800,
    )
    figure.text(
        0.571,
        0.840,
        "B. 세포타입별 질병 모듈 방향 복원",
        ha="left",
        va="center",
        fontsize=17,
        color=COLORS["dark"],
        fontweight=800,
    )
    add_metric_chip(figure, 0.160, 0.790, "● 일반 모듈", medians["blind_centered"], COLORS["teal"])
    add_metric_chip(figure, 0.571, 0.790, "● 같은-donor AD", medians["AD_module"], COLORS["orange"])
    add_metric_chip(figure, 0.766, 0.790, "● Donor-disjoint", medians["cross_disjoint"], COLORS["navy"])

    figure.text(
        0.032,
        0.660,
        "비신경\n6종",
        ha="center",
        va="center",
        fontsize=11.5,
        color=COLORS["green"],
        fontweight=800,
    )
    figure.text(
        0.032,
        0.390,
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
        "핵심  |  Donor-disjoint 질병 방향 복원은 24개 세포타입 모두에서 양의 상관을 보였으며, 중앙값은 ρ = "
        f"{medians['cross_disjoint']:.3f}였다.",
        ha="center",
        va="center",
        fontsize=12.7,
        color="white",
        fontweight=700,
        zorder=3,
    )
    figure.text(
        0.047,
        0.024,
        "Blind = donor별 평균을 제거한 전체 모듈 모양 복원.  AD-module = ADNC·LATE·Lewy·성별·나이 조정.  "
        "Donor-disjoint = 겹치지 않는 donor 절반에서 추정·평가한 질병 방향. 점은 세포타입별 점추정이며 독립 유의성을 뜻하지 않는다.",
        ha="left",
        va="bottom",
        fontsize=9.6,
        color=COLORS["muted"],
    )

    output.mkdir(parents=True, exist_ok=True)
    base = output / "figure_01_celltype_module_and_disease_recovery"
    figure.savefig(base.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(base.with_suffix(".pdf"), facecolor="white")
    figure.savefig(base.with_suffix(".svg"), facecolor="white")
    plt.close(figure)
    return medians


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    source = args.source.resolve()
    output = args.output.resolve()
    rank2 = load_rank2(source)
    medians = draw_figure(rank2, output)

    table_path = output / "figure_01_rank2_e20_celltype_values.tsv"
    rank2.to_csv(table_path, sep="\t", index=False)
    manifest = {
        "schema_version": "prism.presentation.figure01.v1",
        "model": "PRISM personal-rank2 epoch 20",
        "evaluation_split": "validation only",
        "test_used": False,
        "source": str(source),
        "source_sha256": sha256(source),
        "n_celltypes": int(len(rank2)),
        "n_evaluated_modules": 404,
        "n_donors_min": int(rank2["n_donors"].min()),
        "n_donors_max": int(rank2["n_donors"].max()),
        "medians": medians,
        "all_cross_disjoint_positive": bool((rank2["cross_disjoint"] > 0).all()),
        "outputs": {
            "png": "figure_01_celltype_module_and_disease_recovery.png",
            "pdf": "figure_01_celltype_module_and_disease_recovery.pdf",
            "svg": "figure_01_celltype_module_and_disease_recovery.svg",
            "data": table_path.name,
        },
    }
    (output / "figure_01_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
