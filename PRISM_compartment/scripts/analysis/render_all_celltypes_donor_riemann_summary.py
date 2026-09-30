#!/usr/bin/env python3
"""Presentation renderer for the all-cell-type donor Riemann audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


PATHOLOGY = ("Thal", "Braak", "CERAD", "LATE", "Lewy")
NAVY = "#123b72"


def set_style() -> None:
    korean_font = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    family = "DejaVu Sans"
    if korean_font.exists():
        font_manager.fontManager.addfont(str(korean_font))
        family = font_manager.FontProperties(fname=str(korean_font)).get_name()
    plt.rcParams.update(
        {
            "font.family": family,
            "font.size": 8.2,
            "axes.facecolor": "#fbfcfe",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
        }
    )


def heatmap(
    ax,
    values: np.ndarray,
    row_labels: list[str],
    title: str,
    cmap: str,
    *,
    vmin: float,
    vmax: float,
    fmt: str,
    colorbar_label: str,
) -> None:
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(5), PATHOLOGY, rotation=25, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)), row_labels, fontsize=7.4)
    ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold", color="#18243b")
    midpoint = 0.5 * (vmin + vmax)
    for row in range(len(row_labels)):
        for col in range(5):
            value = values[row, col]
            if np.isfinite(value):
                color = "white" if value > midpoint else "#152238"
                ax.text(col, row, format(value, fmt), ha="center", va="center", fontsize=5.8, color=color)
    bar = plt.colorbar(image, ax=ax, fraction=0.027, pad=0.012)
    bar.set_label(colorbar_label, fontsize=7.3)
    bar.ax.tick_params(labelsize=6.7)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arrays", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    set_style()
    data = np.load(args.arrays, allow_pickle=True)
    celltypes = [str(value) for value in data["celltype_names"]]
    test_similarity = np.asarray(data["test_endpoint_metric_similarity"], dtype=float)
    half_similarity = np.asarray(data["split_half_metric_similarity"], dtype=float)
    specificity = np.asarray(data["test_axis_specificity"], dtype=float)
    shortening = np.asarray(data["test_fixed_full5_shortening_pct"], dtype=float)
    shorter_fraction = np.asarray(data["test_fixed_full5_shorter_fraction"], dtype=float)
    response = np.asarray(data["test_response_log_metric_deformation"], dtype=float)
    test_mismatch = 1e6 * (1.0 - test_similarity)
    half_mismatch = 1e6 * (1.0 - half_similarity)

    fig, axes = plt.subplots(2, 3, figsize=(20.5, 13.2), constrained_layout=True)
    heatmap(
        axes[0, 0],
        test_mismatch,
        celltypes,
        "A. locked test → validation 지문 불일치",
        "YlOrRd",
        vmin=0.0,
        vmax=max(50.0, float(np.nanpercentile(test_mismatch, 99))),
        fmt=".1f",
        colorbar_label="불일치 = (1−ρ) × 10⁶; 낮을수록 안정",
    )
    heatmap(
        axes[0, 1],
        half_mismatch,
        celltypes,
        "B. 같은 donor split-half 기술적 불일치",
        "YlGnBu",
        vmin=0.0,
        vmax=max(2.5, float(np.nanpercentile(half_mismatch, 99))),
        fmt=".2f",
        colorbar_label="불일치 = (1−ρ) × 10⁶; 기술적 바닥",
    )
    heatmap(
        axes[0, 2],
        specificity,
        celltypes,
        "C. 같은 질병축 − 다른 질병축 유사도",
        "GnBu",
        vmin=0.0,
        vmax=max(0.06, float(np.nanmax(specificity))),
        fmt=".3f",
        colorbar_label="양수일수록 질병축별 shape 구분",
    )
    heatmap(
        axes[1, 0],
        shortening,
        celltypes,
        "D. test에서 고정된 Figure 13 5D 경로 단축률",
        "magma",
        vmin=0.0,
        vmax=max(14.0, float(np.nanmax(shortening))),
        fmt=".1f",
        colorbar_label="직선 대비 Fisher 길이 단축률 (%)",
    )
    heatmap(
        axes[1, 1],
        shorter_fraction,
        celltypes,
        "E. 고정 5D 경로가 더 짧은 test donor 비율",
        "Blues",
        vmin=0.0,
        vmax=1.0,
        fmt=".2f",
        colorbar_label="1.00 = 해당 test donor 모두에서 더 짧음",
    )
    heatmap(
        axes[1, 2],
        response,
        celltypes,
        "F. 탐색적 personal×병리반응의 기하 변형",
        "Purples",
        vmin=0.0,
        vmax=max(1.25, float(np.nanmax(response))),
        fmt=".2f",
        colorbar_label="정규화 metric의 log-Euclidean 변화",
    )
    fig.suptitle(
        "PRISM 공통 Riemann 지문은 24개 세포타입의 locked test donor 문맥에서도 유지됐다",
        fontsize=20.5,
        fontweight="bold",
        color=NAVY,
    )
    median_similarity = float(np.nanmedian(test_similarity))
    median_specificity = float(np.nanmedian(specificity))
    median_shortening = float(np.nanmedian(shortening))
    positive_specificity = int(np.sum(specificity > 0))
    fixed_shorter = int(np.sum(shorter_fraction == 1.0))
    total = int(np.isfinite(shorter_fraction).sum())
    fig.text(
        0.5,
        -0.014,
        f"요약 | test 지문 ρ 중앙값={median_similarity:.6f} · 같은 질병축 우위={positive_specificity}/120 "
        f"(중앙값 Δ={median_specificity:.3f}) · 고정 5D 경로 test 전원 단축={fixed_shorter}/{total} "
        f"(단축률 중앙값={median_shortening:.2f}%)",
        ha="center",
        fontsize=11.2,
        color=NAVY,
        fontweight="bold",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
