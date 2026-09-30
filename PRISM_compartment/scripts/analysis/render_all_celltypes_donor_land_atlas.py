#!/usr/bin/env python3
"""Render separate validation and locked-test Figure-13-style 2D+3D atlases."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import RegularGridInterpolator


BLUE = "#1464c0"
RED = "#d82e4b"
GREEN = "#009b72"
ORANGE = "#ef6c23"
GREY = "#4f5965"
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
            "font.size": 8.5,
            "axes.facecolor": "#fbfcfe",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
        }
    )


def safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def normalized(value: np.ndarray) -> np.ndarray:
    return value / max(float(np.nanmax(value)), 1e-300)


def mixture_height(normal: np.ndarray, disease: np.ndarray) -> np.ndarray:
    return 0.5 * (normalized(normal) + normalized(disease))


def path_height(x: np.ndarray, y: np.ndarray, height: np.ndarray, path: np.ndarray, offset: float) -> np.ndarray:
    interpolator = RegularGridInterpolator((y, x), height, bounds_error=False, fill_value=0.0)
    return interpolator(np.column_stack((path[:, 1], path[:, 0]))) + offset


def facecolors(normal: np.ndarray, disease: np.ndarray) -> np.ndarray:
    n = normalized(normal)
    d = normalized(disease)
    mix = d / np.maximum(n + d, 1e-15)
    strength = np.maximum(n, d)
    color = plt.get_cmap("coolwarm")(mix)
    color[..., :3] = 0.18 + 0.82 * color[..., :3]
    color[..., 3] = 0.28 + 0.68 * np.sqrt(strength)
    return color


def contour_levels(value: np.ndarray) -> np.ndarray:
    peak = float(np.nanmax(value))
    return np.linspace(0.16, 0.91, 7) * peak


def render_page(
    celltype: str,
    pathology: list[str],
    split_label: str,
    donor_count: int,
    x: np.ndarray,
    y: np.ndarray,
    normal: np.ndarray,
    disease: np.ndarray,
    geodesic: np.ndarray,
    direct: np.ndarray,
    shortening: np.ndarray,
    shortening_q10: np.ndarray,
    shortening_q90: np.ndarray,
    shorter_fraction: np.ndarray,
    capture: np.ndarray,
) -> plt.Figure:
    fig = plt.figure(figsize=(21.5, 10.6), facecolor="white")
    grid = fig.add_gridspec(2, 5, left=0.035, right=0.985, top=0.83, bottom=0.105, hspace=0.28, wspace=0.18)
    xx, yy = np.meshgrid(x, y)
    for axis, name in enumerate(pathology):
        path = geodesic[axis]
        straight = direct[axis]
        max_off = float(np.max(path[:, 1]))
        zoom_y = max(0.16, min(float(y[-1]), 2.15 * max_off + 0.035))
        keep_y = y <= zoom_y + 1e-8
        if np.sum(keep_y) < 5:
            keep_y[:5] = True
        y_view = y[keep_y]
        n_view = normal[axis, keep_y]
        d_view = disease[axis, keep_y]
        xx_view, yy_view = np.meshgrid(x, y_view)

        ax = fig.add_subplot(grid[0, axis])
        ax.contour(xx_view, yy_view, n_view, levels=contour_levels(n_view), colors=BLUE, linewidths=1.45)
        ax.contour(
            xx_view,
            yy_view,
            d_view,
            levels=contour_levels(d_view),
            colors=RED,
            linewidths=1.45,
            linestyles="--",
        )
        ax.plot(straight[:, 0], straight[:, 1], color=ORANGE, lw=2.5, ls="--", zorder=8)
        ax.plot(path[:, 0], path[:, 1], color=GREEN, lw=3.6, zorder=9)
        peak = int(np.argmax(path[:, 1]))
        ax.scatter(path[peak, 0], path[peak, 1], s=34, color=GREEN, edgecolor="white", lw=0.7, zorder=10)
        ax.scatter((0, 1), (0, 0), c=(BLUE, RED), s=28, zorder=11)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, zoom_y)
        ax.grid(alpha=0.14)
        ax.set_xlabel("named pathology severity")
        if axis == 0:
            ax.set_ylabel("off-axis coordinate (확대)")
        ax.set_title(
            f"{name}\n5D short={shortening[axis]:.2f}% "
            f"[{shortening_q10[axis]:.2f},{shortening_q90[axis]:.2f}]\n"
            f"off={max_off:.3f} · path capture={100*capture[axis]:.0f}%",
            fontsize=10.3,
            fontweight="bold",
        )

        ax3 = fig.add_subplot(grid[1, axis], projection="3d")
        height = mixture_height(normal[axis], disease[axis])
        h_view = height[keep_y]
        colors = facecolors(normal[axis], disease[axis])[keep_y]
        ax3.plot_surface(
            xx_view,
            yy_view,
            h_view,
            facecolors=colors,
            rstride=1,
            cstride=1,
            linewidth=0,
            antialiased=True,
            shade=False,
        )
        z_geo = path_height(x, y, height, path, 0.075)
        z_euclid = path_height(x, y, height, straight, 0.035)
        ax3.plot(path[:, 0], path[:, 1], z_geo, color=GREEN, lw=4.2, zorder=15)
        ax3.plot(straight[:, 0], straight[:, 1], z_euclid, color=GREY, lw=2.4, ls=":", zorder=14)
        chord_z = np.linspace(z_geo[0], z_geo[-1], len(straight)) + 0.055
        ax3.plot(straight[:, 0], straight[:, 1], chord_z, color=ORANGE, lw=2.8, ls="--", zorder=16)
        ax3.scatter(path[peak, 0], path[peak, 1], z_geo[peak], color=GREEN, s=28, edgecolor="white", lw=0.5)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, zoom_y)
        ax3.set_zlim(0, 1.22)
        ax3.view_init(elev=28, azim=-57)
        ax3.set_xlabel("severity", labelpad=-1)
        if axis == 0:
            ax3.set_ylabel("off-axis", labelpad=-1)
        ax3.set_zlabel("LAND", labelpad=-3)
        ax3.tick_params(labelsize=6.8, pad=0)
        ax3.set_title(
            f"shorter donor={100*shorter_fraction[axis]:.0f}%",
            fontsize=9.5,
            pad=0,
        )

    locked = " · path/model 재조정 없음" if split_label == "LOCKED TEST" else ""
    fig.suptitle(
        f"{celltype} | {split_label} donor-conditioned common Riemann fingerprint",
        fontsize=21,
        fontweight="bold",
        color=NAVY,
        y=0.965,
    )
    fig.text(
        0.5,
        0.895,
        f"donor {donor_count}명{locked} | green=저장된 실제 5D Fisher path의 투영 | orange=병리축 직선",
        ha="center",
        fontsize=12.2,
        color="#536077",
    )
    fig.legend(
        handles=(
            Line2D([0], [0], color=BLUE, lw=2, label="normal LAND"),
            Line2D([0], [0], color=RED, lw=2, ls="--", label="disease LAND"),
            Line2D([0], [0], color=GREEN, lw=4, label="saved full-5D Fisher geodesic (projected)"),
            Line2D([0], [0], color=GREY, lw=2.4, ls=":", label="coordinate-Euclidean path (surface-lifted)"),
            Line2D([0], [0], color=ORANGE, lw=2.8, ls="--", label="literal straight 3-D chord"),
        ),
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=10.2,
        bbox_to_anchor=(0.5, 0.018),
    )
    fig.text(
        0.5,
        0.003,
        "3D 높이는 LAND density의 scalar lift이며 5D manifold의 등거리 embedding이나 시간적 질병 진행경로가 아니다.",
        ha="center",
        fontsize=8.5,
        color="#697386",
    )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arrays", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    set_style()
    data = np.load(args.arrays, allow_pickle=True)
    celltypes = [str(value) for value in data["celltype_names"]]
    pathology = [str(value) for value in data["pathology_names"]]
    split_names = [str(value) for value in data["split_names"]]
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    normal = np.asarray(data["normal_density_template"], dtype=float)
    disease = np.asarray(data["disease_density_template"], dtype=float)
    geodesic = np.asarray(data["projected_saved_full5_path"], dtype=float)
    direct = np.asarray(data["direct_coordinate_path"], dtype=float)
    shortening = np.asarray(data["shortening_median_pct"], dtype=float)
    shortening_q10 = np.asarray(data["shortening_q10_pct"], dtype=float)
    shortening_q90 = np.asarray(data["shortening_q90_pct"], dtype=float)
    shorter_fraction = np.asarray(data["shorter_fraction"], dtype=float)
    capture = np.asarray(data["projection_capture"], dtype=float)
    donor_counts = (
        np.asarray(data["validation_donor_count"], dtype=int),
        np.asarray(data["locked_test_donor_count"], dtype=int),
    )
    split_labels = ("VALIDATION", "LOCKED TEST")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for split_index, split_name in enumerate(split_names):
        split_dir = args.out_dir / split_name
        split_dir.mkdir(parents=True, exist_ok=True)
        pdf_path = args.out_dir / f"figure_07_{split_name}_all24_2d3d_atlas.pdf"
        with PdfPages(pdf_path) as pdf:
            for celltype, name in enumerate(celltypes):
                fig = render_page(
                    name,
                    pathology,
                    split_labels[split_index],
                    int(donor_counts[split_index][celltype]),
                    x,
                    y,
                    normal[split_index, celltype],
                    disease[split_index, celltype],
                    geodesic[celltype],
                    direct[celltype],
                    shortening[split_index, celltype],
                    shortening_q10[split_index, celltype],
                    shortening_q90[split_index, celltype],
                    shorter_fraction[split_index, celltype],
                    capture[celltype],
                )
                pdf.savefig(fig, bbox_inches="tight")
                page = split_dir / f"figure_07_{split_name}_page_{celltype + 1:02d}_{safe_name(name)}.png"
                fig.savefig(page, dpi=185, bbox_inches="tight")
                plt.close(fig)
                print(f"[render-land-atlas] {split_name} {celltype + 1}/{len(celltypes)} {name}", flush=True)
        print(f"[render-land-atlas] wrote {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
