#!/usr/bin/env python3
"""Render donor overlays for common-only PRISM endpoint fingerprints.

Each donor's frozen ordinal-Fisher metric has already been projected into the
fixed Figure-13 named-axis/off-axis chart.  This script converts each 2x2 SPD
metric into a closed local-Fisher HDR ellipse.  Validation donors define the
consensus and its 10--90% radial band; locked-test donors are overlaid without
refitting or retuning.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


NAVY = "#173f76"
VAL_INDIVIDUAL = "#84a9d8"
TEST_INDIVIDUAL = "#e88b50"
TEST_MEDIAN = "#e46723"
NORMAL_CONSENSUS = "#075ca8"
DISEASE_CONSENSUS = "#b21f3d"
GRID = "#d7dee8"
TEXT = "#4d5b70"
PATHOLOGY_LABELS = {
    "thal": "Thal",
    "braak": "Braak",
    "cerad": "CERAD",
    "late": "LATE",
    "lewy": "Lewy",
}


def safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def register_font() -> None:
    candidates = (
        Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"),
        Path("/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc"),
    )
    for path in candidates:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(path)).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False


def regularize_metric(metric: np.ndarray) -> np.ndarray:
    metric = 0.5 * (metric + metric.T)
    values, vectors = np.linalg.eigh(metric)
    floor = max(float(np.max(values)) * 1e-7, 1e-10)
    values = np.maximum(values, floor)
    return (vectors * values) @ vectors.T


def metric_radius(metric: np.ndarray, sigma: float, unit: np.ndarray, mass: float) -> np.ndarray:
    """Radius of a constant-mass contour for exp(-delta' G delta / 2 sigma^2)."""
    metric = regularize_metric(np.asarray(metric, dtype=np.float64))
    quadratic = np.einsum("ti,ij,tj->t", unit, metric, unit)
    threshold = -2.0 * np.log1p(-mass)
    return float(sigma) * np.sqrt(threshold / np.maximum(quadratic, 1e-12))


def area_normalize(radius: np.ndarray) -> np.ndarray:
    scale = np.sqrt(np.mean(np.square(radius), axis=-1, keepdims=True))
    return radius / np.maximum(scale, 1e-12)


def radial_overlap(radius_a: np.ndarray, radius_b: np.ndarray) -> float:
    """Intersection-over-union for two centred star-shaped contours."""
    a2 = np.square(area_normalize(np.asarray(radius_a, dtype=np.float64)))
    b2 = np.square(area_normalize(np.asarray(radius_b, dtype=np.float64)))
    return float(np.sum(np.minimum(a2, b2)) / np.maximum(np.sum(np.maximum(a2, b2)), 1e-12))


def prepare(data: np.lib.npyio.NpzFile, mass: float, angle_count: int):
    donors = [str(value) for value in data["donor_names"]]
    splits = np.asarray([str(value) for value in data["donor_split"]], dtype=object)
    celltypes = [str(value) for value in data["celltype_names"]]
    pathologies = [str(value) for value in data["pathology_names"]]
    eligible = np.asarray(data["eligible"], dtype=bool)
    sigma = np.asarray(data["sigma_fisher"], dtype=np.float64)
    theta = np.linspace(0.0, 2.0 * np.pi, angle_count, endpoint=True)
    unit = np.column_stack((np.cos(theta), np.sin(theta)))

    endpoint_metrics = {
        "normal": np.asarray(data["normal_metric_g2"], dtype=np.float64),
        "disease": np.asarray(data["disease_metric_g2"], dtype=np.float64),
    }
    records = {}
    overlap = np.full((2, len(donors), len(celltypes), len(pathologies)), np.nan, dtype=np.float32)
    scale_ratio = np.full_like(overlap, np.nan)
    validation_consensus = np.full(
        (2, len(celltypes), len(pathologies), angle_count), np.nan, dtype=np.float32
    )
    test_consensus = np.full_like(validation_consensus, np.nan)

    for endpoint_index, endpoint in enumerate(("normal", "disease")):
        metrics = endpoint_metrics[endpoint]
        for celltype in range(len(celltypes)):
            val_mask = eligible[:, celltype] & (splits == "val")
            test_mask = eligible[:, celltype] & (splits == "test")
            for axis in range(len(pathologies)):
                radii = np.full((len(donors), angle_count), np.nan, dtype=np.float64)
                for donor in np.flatnonzero(val_mask | test_mask):
                    radii[donor] = metric_radius(metrics[donor, celltype, axis], sigma[celltype, axis], unit, mass)
                val_radii = radii[val_mask]
                test_radii = radii[test_mask]
                val_median = np.median(val_radii, axis=0)
                test_median = np.median(test_radii, axis=0)
                validation_consensus[endpoint_index, celltype, axis] = val_median
                test_consensus[endpoint_index, celltype, axis] = test_median
                val_scale = float(np.sqrt(np.mean(np.square(val_median))))
                for donor in np.flatnonzero(val_mask | test_mask):
                    overlap[endpoint_index, donor, celltype, axis] = radial_overlap(radii[donor], val_median)
                    donor_scale = float(np.sqrt(np.mean(np.square(radii[donor]))))
                    scale_ratio[endpoint_index, donor, celltype, axis] = donor_scale / max(val_scale, 1e-12)
                records[(endpoint, celltype, axis)] = {
                    "radii": radii,
                    "val_mask": val_mask,
                    "test_mask": test_mask,
                    "val_median": val_median,
                    "test_median": test_median,
                    "val_low": np.quantile(val_radii, 0.10, axis=0),
                    "val_high": np.quantile(val_radii, 0.90, axis=0),
                }
    return donors, splits, celltypes, pathologies, theta, unit, records, overlap, scale_ratio, validation_consensus, test_consensus


def common_bounds(records: dict, endpoint: str, celltype_count: int, axis: int, unit: np.ndarray):
    maximum = 0.0
    for celltype in range(celltype_count):
        record = records[(endpoint, celltype, axis)]
        selected = record["radii"][record["val_mask"] | record["test_mask"]]
        xy = selected[:, :, None] * unit[None, :, :]
        maximum = max(maximum, float(np.nanmax(np.abs(xy))))
    bound = max(maximum * 1.10, 1e-6)
    return (-bound, bound)


def draw_panel(ax, record: dict, unit: np.ndarray, endpoint: str, title: str, bound, panel_overlap: float):
    val_radii = record["radii"][record["val_mask"]]
    test_radii = record["radii"][record["test_mask"]]
    for radius in val_radii:
        xy = radius[:, None] * unit
        ax.plot(xy[:, 0], xy[:, 1], color=VAL_INDIVIDUAL, lw=0.55, alpha=0.42, zorder=2)
    for radius in test_radii:
        xy = radius[:, None] * unit
        ax.plot(xy[:, 0], xy[:, 1], color=TEST_INDIVIDUAL, lw=0.65, alpha=0.58, ls="--", zorder=3)

    low_xy = record["val_low"][:, None] * unit
    high_xy = record["val_high"][:, None] * unit
    band_x = np.concatenate((low_xy[:, 0], high_xy[::-1, 0]))
    band_y = np.concatenate((low_xy[:, 1], high_xy[::-1, 1]))
    consensus_color = NORMAL_CONSENSUS if endpoint == "normal" else DISEASE_CONSENSUS
    ax.fill(band_x, band_y, color=consensus_color, alpha=0.10, linewidth=0.0, zorder=1)
    val_xy = record["val_median"][:, None] * unit
    test_xy = record["test_median"][:, None] * unit
    ax.plot(val_xy[:, 0], val_xy[:, 1], color=consensus_color, lw=2.2, zorder=5)
    ax.plot(test_xy[:, 0], test_xy[:, 1], color=TEST_MEDIAN, lw=1.8, ls="--", zorder=6)
    ax.scatter((0.0,), (0.0,), s=10, color=NAVY, edgecolor="white", linewidth=0.4, zorder=8)
    ax.axhline(0.0, color=GRID, lw=0.45, zorder=0)
    ax.axvline(0.0, color=GRID, lw=0.45, zorder=0)
    ax.set_xlim(bound)
    ax.set_ylim(bound)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(title, fontsize=9.0, fontweight="bold", color=NAVY, pad=2)
    ax.text(
        0.97,
        0.04,
        f"test overlap {panel_overlap:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        color=TEXT,
    )
    ax.text(
        0.03,
        0.04,
        f"val n={val_radii.shape[0]} | test n={test_radii.shape[0]}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.2,
        color=TEXT,
    )
    for spine in ax.spines.values():
        spine.set_color("#aab4c3")
        spine.set_linewidth(0.55)


def legend_handles(endpoint: str):
    consensus_color = NORMAL_CONSENSUS if endpoint == "normal" else DISEASE_CONSENSUS
    return (
        Line2D((0,), (0,), color=VAL_INDIVIDUAL, lw=1.3, label="Validation donor (개별)"),
        Line2D((0,), (0,), color=TEST_INDIVIDUAL, lw=1.3, ls="--", label="Locked-test donor (개별)"),
        Line2D((0,), (0,), color=consensus_color, lw=3.0, label="Validation 공통 윤곽"),
        Line2D((0,), (0,), color=TEST_MEDIAN, lw=2.2, ls="--", label="Locked-test 중앙 윤곽 (평가만)"),
        Patch(facecolor=consensus_color, alpha=0.12, edgecolor="none", label="Validation donor 10–90% 범위"),
    )


def render_atlases(out_dir: Path, celltypes, pathologies, unit, records, overlap, splits):
    test_index = splits == "test"
    for endpoint_index, endpoint in enumerate(("normal", "disease")):
        korean = "정상 끝점" if endpoint == "normal" else "질병 끝점"
        pdf_path = out_dir / f"figure_10_{endpoint}_donor_contour_overlay_atlas.pdf"
        with PdfPages(pdf_path) as writer:
            for axis, pathology in enumerate(pathologies):
                bound = common_bounds(records, endpoint, len(celltypes), axis, unit)
                fig, axes = plt.subplots(4, 6, figsize=(17.5, 11.8))
                for celltype, ax in enumerate(axes.ravel()):
                    values = overlap[endpoint_index, test_index, celltype, axis]
                    panel_overlap = float(np.nanmedian(values))
                    draw_panel(
                        ax,
                        records[(endpoint, celltype, axis)],
                        unit,
                        endpoint,
                        celltypes[celltype],
                        bound,
                        panel_overlap,
                    )
                label = PATHOLOGY_LABELS.get(pathology.lower(), pathology)
                fig.suptitle(
                    f"{label}: 여러 donor에서 반복되는 세포타입별 common-only {korean} 지문",
                    fontsize=19.5,
                    fontweight="bold",
                    color=NAVY,
                    y=0.985,
                )
                fig.text(
                    0.5,
                    0.948,
                    "explicit personal·personal×병리·성별·기술 항 = 0 | donor zclean baseline 문맥 유지 | locked test 재학습 없음",
                    ha="center",
                    fontsize=10.2,
                    color=TEXT,
                )
                fig.legend(handles=legend_handles(endpoint), loc="lower center", ncol=5, frameon=False, fontsize=9.2)
                fig.text(
                    0.5,
                    0.026,
                    "overlap은 크기를 맞춘 닫힌 윤곽의 교집합/합집합이며 1에 가까울수록 모양이 같다. 모든 패널은 해당 병리축에서 동일한 축척이다.",
                    ha="center",
                    fontsize=8.8,
                    color="#667085",
                )
                fig.subplots_adjust(left=0.025, right=0.985, top=0.902, bottom=0.088, wspace=0.12, hspace=0.25)
                png_path = out_dir / f"figure_10_{endpoint}_{axis + 1:02d}_{safe_name(label)}_donor_overlay.png"
                fig.savefig(png_path, dpi=210, bbox_inches="tight")
                writer.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                print(f"[donor-overlay] rendered {endpoint} {label}", flush=True)


def render_representative(out_dir: Path, celltypes, pathologies, unit, records, overlap, splits, chosen_celltype: str):
    celltype = celltypes.index(chosen_celltype) if chosen_celltype in celltypes else 0
    test_index = splits == "test"
    fig, axes = plt.subplots(2, len(pathologies), figsize=(18.2, 7.8))
    for endpoint_index, endpoint in enumerate(("normal", "disease")):
        for axis, pathology in enumerate(pathologies):
            selected = records[(endpoint, celltype, axis)]["radii"]
            selected = selected[records[(endpoint, celltype, axis)]["val_mask"] | records[(endpoint, celltype, axis)]["test_mask"]]
            maximum = float(np.nanmax(np.abs(selected[:, :, None] * unit[None, :, :])))
            bound = (-1.10 * maximum, 1.10 * maximum)
            panel_overlap = float(np.nanmedian(overlap[endpoint_index, test_index, celltype, axis]))
            draw_panel(
                axes[endpoint_index, axis],
                records[(endpoint, celltype, axis)],
                unit,
                endpoint,
                PATHOLOGY_LABELS.get(pathology.lower(), pathology),
                bound,
                panel_overlap,
            )
            if axis == 0:
                axes[endpoint_index, axis].set_ylabel(
                    "정상 끝점" if endpoint == "normal" else "질병 끝점",
                    fontsize=11,
                    fontweight="bold",
                    color=NAVY,
                    labelpad=10,
                )
    fig.suptitle(
        f"{chosen_celltype}: validation 공통 지문이 여러 locked-test donor에서 반복된다",
        fontsize=20,
        fontweight="bold",
        color=NAVY,
        y=0.982,
    )
    fig.text(0.5, 0.936, "위: 정상 끝점 | 아래: 질병 끝점 | explicit personal 항은 0, donor zclean baseline 문맥은 유지", ha="center", fontsize=10.5, color=TEXT)
    fig.legend(handles=legend_handles("normal"), loc="lower center", ncol=5, frameon=False, fontsize=9.2)
    fig.subplots_adjust(left=0.045, right=0.99, top=0.875, bottom=0.105, wspace=0.15, hspace=0.25)
    for suffix in ("png", "pdf"):
        fig.savefig(out_dir / f"figure_10_representative_{safe_name(chosen_celltype)}_normal_disease_donor_overlay.{suffix}", dpi=220 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def render_summary(out_dir: Path, celltypes, pathologies, overlap, scale_ratio, splits):
    test_index = splits == "test"
    median_overlap = np.nanmedian(overlap[:, test_index], axis=1)
    median_scale = np.nanmedian(scale_ratio[:, test_index], axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 10.5), constrained_layout=True)
    endpoint_labels = ("정상 끝점", "질병 끝점")
    images = []
    for endpoint_index, ax in enumerate(axes):
        image = ax.imshow(median_overlap[endpoint_index], vmin=0.98, vmax=1.0, cmap="YlGnBu", aspect="auto")
        images.append(image)
        ax.set_xticks(np.arange(len(pathologies)), [PATHOLOGY_LABELS.get(value.lower(), value) for value in pathologies], rotation=35, ha="right")
        ax.set_yticks(np.arange(len(celltypes)), celltypes)
        ax.set_title(endpoint_labels[endpoint_index], fontsize=15, fontweight="bold", color=NAVY, pad=10)
        for row in range(len(celltypes)):
            for column in range(len(pathologies)):
                value = median_overlap[endpoint_index, row, column]
                ax.text(column, row, f"{value:.3f}", ha="center", va="center", fontsize=6.8, color="white" if value > 0.994 else "#132f52")
        ax.set_xlabel("병리축", color=TEXT)
    colorbar = fig.colorbar(images[0], ax=axes, shrink=0.72, pad=0.02)
    colorbar.set_label("Locked-test → validation 공통 윤곽의 shape overlap", color=TEXT)
    fig.suptitle("Explicit personal 항 제거 후 common-only 지문의 locked-test 재현성", fontsize=20, fontweight="bold", color=NAVY)
    for suffix in ("png", "pdf"):
        fig.savefig(out_dir / f"figure_10_donor_shape_overlap_summary.{suffix}", dpi=220 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return median_overlap, median_scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--mass", type=float, default=0.68)
    parser.add_argument("--angle-count", type=int, default=361)
    parser.add_argument("--representative-celltype", default="L4 IT")
    args = parser.parse_args()
    register_font()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data = np.load(args.input, allow_pickle=True)
    prepared = prepare(data, args.mass, args.angle_count)
    donors, splits, celltypes, pathologies, theta, unit, records, overlap, scale_ratio, validation_consensus, test_consensus = prepared
    render_atlases(args.out_dir, celltypes, pathologies, unit, records, overlap, splits)
    render_representative(args.out_dir, celltypes, pathologies, unit, records, overlap, splits, args.representative_celltype)
    median_overlap, median_scale = render_summary(args.out_dir, celltypes, pathologies, overlap, scale_ratio, splits)

    np.savez_compressed(
        args.out_dir / "donor_contour_shape_statistics.npz",
        theta=theta.astype(np.float32),
        donor_names=np.asarray(donors, dtype=object),
        donor_split=splits,
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(pathologies, dtype=object),
        validation_consensus_radius=validation_consensus,
        locked_test_median_radius=test_consensus,
        donor_to_validation_shape_overlap=overlap,
        donor_to_validation_scale_ratio=scale_ratio,
        locked_test_median_shape_overlap=median_overlap,
        locked_test_median_scale_ratio=median_scale,
    )
    val_index = splits == "val"
    test_index = splits == "test"
    summary = {
        "schema_version": "prism.common_only_donor_contour_overlay.v1",
        "input": str(args.input.resolve()),
        "checkpoint_sha256": str(data["checkpoint_sha256"].item()),
        "checkpoint_epoch": int(data["checkpoint_epoch"].item()),
        "contour_mass": float(args.mass),
        "test_refit": False,
        "personal_baseline": 0,
        "personal_response": 0,
        "sex_score": 0,
        "technology_score": 0,
        "validation_donor_overlap": {
            endpoint: {
                "median": float(np.nanmedian(overlap[index, val_index])),
                "q10": float(np.nanquantile(overlap[index, val_index], 0.10)),
                "q90": float(np.nanquantile(overlap[index, val_index], 0.90)),
            }
            for index, endpoint in enumerate(("normal", "disease"))
        },
        "locked_test_donor_overlap": {
            endpoint: {
                "median": float(np.nanmedian(overlap[index, test_index])),
                "q10": float(np.nanquantile(overlap[index, test_index], 0.10)),
                "q90": float(np.nanquantile(overlap[index, test_index], 0.90)),
                "minimum_celltype_pathology_median": float(np.nanmin(median_overlap[index])),
                "maximum_celltype_pathology_median": float(np.nanmax(median_overlap[index])),
            }
            for index, endpoint in enumerate(("normal", "disease"))
        },
        "locked_test_scale_ratio": {
            endpoint: {
                "median": float(np.nanmedian(scale_ratio[index, test_index])),
                "q10": float(np.nanquantile(scale_ratio[index, test_index], 0.10)),
                "q90": float(np.nanquantile(scale_ratio[index, test_index], 0.90)),
            }
            for index, endpoint in enumerate(("normal", "disease"))
        },
        "interpretation": "Shape overlap is area-normalized radial intersection-over-union; scale ratio is donor equivalent radius divided by validation-consensus equivalent radius.",
    }
    with (args.out_dir / "donor_contour_shape_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
