#!/usr/bin/env python3
"""Overlay donor-specific multi-level local-Fisher fingerprint contours.

The saved common 5-D path and display chart are fixed for every donor.  Each
donor contributes its own ordinal-Fisher metric along that path.  The resulting
anchored Gaussian mixture is the same Figure-11-style local fingerprint used in
the presentation, not an endpoint-ellipse approximation.
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
import numpy as np


NAVY = "#173f76"
NORMAL = "#0067c5"
DISEASE = "#c81932"
TEXT = "#53637a"
GRID = "#cfd7e3"
HDR_MASSES = (0.86, 0.74, 0.62, 0.50, 0.38, 0.26, 0.16)
PATHOLOGY_LABELS = {"thal": "Thal", "braak": "Braak", "cerad": "CERAD", "late": "LATE", "lewy": "Lewy"}
LINESTYLES = ("-", "--", ":", "-.", "-", "--", ":", "-.", "--")


def safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def register_font() -> None:
    path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if path.exists():
        font_manager.fontManager.addfont(str(path))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(path)).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42


def hdr_threshold(density: np.ndarray, mass: float) -> float:
    values = np.asarray(density, dtype=np.float64).ravel()
    order = np.argsort(values)[::-1]
    cumulative = np.cumsum(values[order])
    cumulative /= max(float(cumulative[-1]), 1e-15)
    index = min(int(np.searchsorted(cumulative, mass)), len(order) - 1)
    return float(values[order[index]])


def regularize_metric(metric: np.ndarray) -> np.ndarray:
    metric = 0.5 * (metric + np.swapaxes(metric, -1, -2))
    values, vectors = np.linalg.eigh(metric)
    floor = np.maximum(np.max(values, axis=-1) * 1e-8, 1e-12)
    values = np.maximum(values, floor[..., None])
    return np.einsum("...ik,...k,...jk->...ij", vectors, values, vectors)


def densify(path5: np.ndarray, path_xy: np.ndarray, metric: np.ndarray, count: int):
    segment = np.linalg.norm(np.diff(path5, axis=0), axis=1)
    parameter = np.concatenate(([0.0], np.cumsum(segment)))
    keep = np.concatenate(([True], np.diff(parameter) > 1e-12))
    parameter = parameter[keep]
    path_xy = path_xy[keep]
    metric = metric[keep]
    if parameter[-1] <= 1e-12:
        parameter = np.linspace(0.0, 1.0, len(parameter))
    else:
        parameter /= parameter[-1]
    dense = np.linspace(0.0, 1.0, count)
    dense_xy = np.column_stack([np.interp(dense, parameter, path_xy[:, dim]) for dim in range(2)])
    dense_metric = np.empty((count, 2, 2), dtype=np.float64)
    for row in range(2):
        for column in range(2):
            dense_metric[:, row, column] = np.interp(dense, parameter, metric[:, row, column])
    return dense, dense_xy, regularize_metric(dense_metric)


def display_bounds(path_xy: np.ndarray, path_metric: np.ndarray, sigma: float, scale: float):
    covariance = sigma**2 * np.linalg.inv(path_metric)
    std_x = np.sqrt(np.maximum(covariance[:, 0, 0], 1e-12))
    std_y = np.sqrt(np.maximum(covariance[:, 1, 1], 1e-12))
    left = max(0.15, min(0.38, 2.35 * float(std_x[0])))
    right = max(0.15, min(0.38, 2.35 * float(std_x[-1])))
    below = max(0.12, min(0.52, 2.35 * float(max(std_y[0], std_y[-1]))))
    above = max(0.15, min(0.58, 2.15 * float(np.percentile(std_y, 90))))
    raw_x = (-left, 1.0 + right)
    raw_y = (-below, float(np.max(path_xy[:, 1])) + above)
    x_mid, y_mid = 0.5 * sum(raw_x), 0.5 * sum(raw_y)
    x_half, y_half = 0.5 * (raw_x[1] - raw_x[0]) * scale, 0.5 * (raw_y[1] - raw_y[0]) * scale
    return (x_mid - x_half, x_mid + x_half), (y_mid - y_half, y_mid + y_half)


def densities_on_grid(grid_x, grid_y, path_xy, metric, parameter, sigma):
    xx, yy = np.meshgrid(grid_x, grid_y, indexing="xy")
    delta = np.stack((xx[..., None] - path_xy[:, 0], yy[..., None] - path_xy[:, 1]), axis=-1)
    local_square = np.einsum("...ni,nij,...nj->...n", delta, metric, delta, optimize=True)
    kernel = np.exp(-0.5 * local_square / max(float(sigma) ** 2, 1e-30))
    normal_weight = np.exp(-parameter / 0.26)
    disease_weight = np.exp(-(1.0 - parameter) / 0.26)
    normal_weight[0] += 0.80 * normal_weight.sum()
    disease_weight[-1] += 0.80 * disease_weight.sum()
    normal = np.einsum("...n,n->...", kernel, normal_weight, optimize=True)
    disease = np.einsum("...n,n->...", kernel, disease_weight, optimize=True)
    normal /= max(float(normal.max()), 1e-300)
    disease /= max(float(disease.max()), 1e-300)
    return normal.astype(np.float32), disease.astype(np.float32)


def fingerprint_vector(density: np.ndarray) -> np.ndarray:
    values = np.log(np.maximum(np.asarray(density, dtype=np.float64).ravel(), 1e-8))
    values -= values.mean()
    return values / max(float(np.linalg.norm(values)), 1e-30)


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = fingerprint_vector(left)
    right = fingerprint_vector(right)
    return float(left @ right)


def is_closed(density: np.ndarray) -> bool:
    outer = hdr_threshold(density, HDR_MASSES[0])
    mask = density >= outer
    return not bool(np.any(np.concatenate((mask[0], mask[-1], mask[:, 0], mask[:, -1]))))


def prepare_record(data, celltype: int, axis: int, test_global: np.ndarray, nx: int, ny: int, dense_nodes: int, base_scale: float):
    eligible = np.asarray(data["eligible"], dtype=bool)[:, celltype]
    split = np.asarray([str(value) for value in data["donor_split"]], dtype=object)
    val = np.flatnonzero(eligible & (split == "val"))
    test = np.asarray([donor for donor in test_global if eligible[donor]], dtype=int)
    path5 = np.asarray(data["path5"][celltype, axis], dtype=np.float64)
    path_xy = np.asarray(data["path_xy"][celltype, axis], dtype=np.float64)
    path_metric = np.asarray(data["path_metric_g2"][:, celltype, axis], dtype=np.float64)
    sigma = float(data["sigma_fisher"][celltype, axis])
    consensus_raw = np.nanmedian(path_metric[val], axis=0)
    parameter, dense_xy, consensus_metric = densify(path5, path_xy, consensus_raw, dense_nodes)

    for scale in (base_scale, base_scale * 1.30, base_scale * 1.65):
        x_bounds, y_bounds = display_bounds(dense_xy, consensus_metric, sigma, scale)
        grid_x = np.linspace(*x_bounds, nx)
        grid_y = np.linspace(*y_bounds, ny)
        consensus_normal, consensus_disease = densities_on_grid(
            grid_x, grid_y, dense_xy, consensus_metric, parameter, sigma
        )
        test_normal, test_disease = [], []
        for donor in test:
            _, _, donor_metric = densify(path5, path_xy, path_metric[donor], dense_nodes)
            normal, disease = densities_on_grid(grid_x, grid_y, dense_xy, donor_metric, parameter, sigma)
            test_normal.append(normal)
            test_disease.append(disease)
        test_normal = np.stack(test_normal)
        test_disease = np.stack(test_disease)
        all_closed = all(
            is_closed(value)
            for value in (consensus_normal, consensus_disease, *test_normal, *test_disease)
        )
        if all_closed:
            break
    return {
        "grid_x": grid_x,
        "grid_y": grid_y,
        "consensus_normal": consensus_normal,
        "consensus_disease": consensus_disease,
        "test_normal": test_normal,
        "test_disease": test_disease,
        "test_donors": test,
        "closed": all_closed,
        "bound_scale": scale,
        "normal_correlation": np.asarray([correlation(value, consensus_normal) for value in test_normal]),
        "disease_correlation": np.asarray([correlation(value, consensus_disease) for value in test_disease]),
    }


def contour_levels(density: np.ndarray):
    return sorted({hdr_threshold(density, mass) for mass in HDR_MASSES})


def page_bounds(records: list[dict], endpoint: str):
    anchor = 0.0 if endpoint == "normal" else 1.0
    x_min, x_max, y_min, y_max = np.inf, -np.inf, np.inf, -np.inf
    for record in records:
        densities = [record[f"consensus_{endpoint}"], *record[f"test_{endpoint}"]]
        for density in densities:
            outer = min(contour_levels(density))
            rows, columns = np.where(density >= outer)
            x = record["grid_x"][columns] - anchor
            y = record["grid_y"][rows]
            x_min, x_max = min(x_min, float(x.min())), max(x_max, float(x.max()))
            y_min, y_max = min(y_min, float(y.min())), max(y_max, float(y.max()))
    x_pad = 0.08 * max(x_max - x_min, 1e-6)
    y_pad = 0.08 * max(y_max - y_min, 1e-6)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def draw_panel(ax, record: dict, endpoint: str, celltype_name: str, colors, test_order, bounds):
    anchor = 0.0 if endpoint == "normal" else 1.0
    consensus = record[f"consensus_{endpoint}"]
    tests = record[f"test_{endpoint}"]
    correlations = record[f"{endpoint}_correlation"]
    xx, yy = np.meshgrid(record["grid_x"] - anchor, record["grid_y"], indexing="xy")
    endpoint_color = NORMAL if endpoint == "normal" else DISEASE
    levels = contour_levels(consensus)
    fill_levels = [levels[0], *levels[2::2], 1.01]
    fill_levels = sorted(set(fill_levels))
    ax.contourf(xx, yy, consensus, levels=fill_levels, cmap="Blues" if endpoint == "normal" else "Reds", alpha=0.44, zorder=0)
    ax.contour(xx, yy, consensus, levels=levels, colors=endpoint_color, linewidths=np.linspace(0.9, 1.9, len(levels)), zorder=2)
    for local_index, donor in enumerate(record["test_donors"]):
        global_test_index = int(np.flatnonzero(test_order == donor)[0])
        donor_levels = contour_levels(tests[local_index])
        ax.contour(
            xx,
            yy,
            tests[local_index],
            levels=donor_levels,
            colors=[colors[global_test_index]],
            linewidths=np.linspace(0.48, 0.92, len(donor_levels)),
            linestyles=LINESTYLES[global_test_index % len(LINESTYLES)],
            alpha=0.72,
            zorder=4 + global_test_index,
        )
    ax.scatter((0.0,), (0.0,), s=13, color=endpoint_color, edgecolor="white", linewidth=0.55, zorder=20)
    ax.axhline(0.0, color=GRID, lw=0.45, zorder=-1)
    ax.axvline(0.0, color=GRID, lw=0.45, zorder=-1)
    ax.set_xlim(*bounds[0])
    ax.set_ylim(*bounds[1])
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(celltype_name, fontsize=8.7, fontweight="bold", color=NAVY, pad=2)
    ax.text(0.97, 0.04, f"test median ρ={np.nanmedian(correlations):.5f}", transform=ax.transAxes, ha="right", fontsize=6.2, color=TEXT)
    ax.text(0.03, 0.04, f"n={len(tests)}", transform=ax.transAxes, ha="left", fontsize=6.2, color=TEXT)
    for spine in ax.spines.values():
        spine.set_color("#aab4c3")
        spine.set_linewidth(0.55)


def render_page(out_dir, writer, endpoint, axis, pathology, celltypes, records, colors, test_order, donor_names):
    bounds = page_bounds(records, endpoint)
    fig, axes = plt.subplots(4, 6, figsize=(18.2, 12.2))
    for celltype, ax in enumerate(axes.ravel()):
        draw_panel(ax, records[celltype], endpoint, celltypes[celltype], colors, test_order, bounds)
    label = PATHOLOGY_LABELS.get(pathology.lower(), pathology)
    korean = "정상" if endpoint == "normal" else "질병"
    fig.suptitle(f"{label}: locked-test donor마다 반복되는 세포타입별 {korean} 다중등고선 지문", fontsize=19.5, fontweight="bold", color=NAVY, y=0.987)
    fig.text(0.5, 0.953, "Figure-11 방식 | 같은 저장된 5-D 공통 경로 | donor별 경로상 ordinal-Fisher metric | test 재적합 없음", ha="center", fontsize=10.0, color=TEXT)
    handles = [Line2D((0,), (0,), color=NORMAL if endpoint == "normal" else DISEASE, lw=3, label="Validation consensus")]
    handles += [Line2D((0,), (0,), color=colors[index], lw=1.5, ls=LINESTYLES[index], label=f"T{index + 1}") for index in range(len(test_order))]
    fig.legend(handles=handles, loc="lower center", ncol=10, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.052))
    mapping = "  |  ".join(f"T{index + 1}={donor_names[donor]}" for index, donor in enumerate(test_order))
    fig.text(0.5, 0.031, mapping, ha="center", fontsize=7.1, color="#667085")
    fig.text(0.5, 0.014, "각 donor의 7개 HDR 등고선을 직접 중첩했다. ρ는 전체 local density shape의 validation consensus 상관이다.", ha="center", fontsize=8.3, color=TEXT)
    fig.subplots_adjust(left=0.022, right=0.99, top=0.91, bottom=0.088, wspace=0.10, hspace=0.23)
    png = out_dir / f"figure_11_{endpoint}_{axis + 1:02d}_{safe_name(label)}_locked_test_donor_multicontour.png"
    fig.savefig(png, dpi=215, bbox_inches="tight")
    writer.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    return png


def render_representative(out_dir, endpoint, celltype_name, celltype_index, pathologies, representative_records, colors, test_order, donor_names):
    fig, axes = plt.subplots(1, len(pathologies), figsize=(18.5, 4.5))
    for axis, ax in enumerate(axes):
        record = representative_records[axis]
        bounds = page_bounds([record], endpoint)
        draw_panel(ax, record, endpoint, PATHOLOGY_LABELS.get(pathologies[axis].lower(), pathologies[axis]), colors, test_order, bounds)
    korean = "정상" if endpoint == "normal" else "질병"
    fig.suptitle(f"{celltype_name}: 9명 locked-test donor에서 반복되는 {korean} 다중등고선 지문", fontsize=19.0, fontweight="bold", color=NAVY, y=0.985)
    fig.text(0.5, 0.913, "색·선형 = 서로 다른 test donor | 붉거나 푸른 배경 = validation consensus", ha="center", fontsize=10.0, color=TEXT)
    handles = [Line2D((0,), (0,), color=NORMAL if endpoint == "normal" else DISEASE, lw=3, label="Validation consensus")]
    handles += [Line2D((0,), (0,), color=colors[index], lw=1.5, ls=LINESTYLES[index], label=f"T{index + 1}") for index in range(len(test_order))]
    fig.legend(handles=handles, loc="lower center", ncol=10, frameon=False, fontsize=8.5)
    fig.subplots_adjust(left=0.025, right=0.99, top=0.82, bottom=0.20, wspace=0.16)
    for suffix in ("png", "pdf"):
        fig.savefig(out_dir / f"figure_11_representative_{safe_name(celltype_name)}_{endpoint}_donor_multicontour.{suffix}", dpi=220 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--grid-x", type=int, default=130)
    parser.add_argument("--grid-y", type=int, default=100)
    parser.add_argument("--dense-nodes", type=int, default=81)
    parser.add_argument("--bound-scale", type=float, default=2.4)
    parser.add_argument("--representative-celltype", default="L4 IT")
    args = parser.parse_args()
    register_font()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data = np.load(args.input, allow_pickle=True)
    donor_names = [str(value) for value in data["donor_names"]]
    splits = np.asarray([str(value) for value in data["donor_split"]], dtype=object)
    celltypes = [str(value) for value in data["celltype_names"]]
    pathologies = [str(value) for value in data["pathology_names"]]
    test_order = np.flatnonzero(splits == "test")
    colors = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(test_order)))
    representative_index = celltypes.index(args.representative_celltype)
    representative_records = {}
    shape_median = np.full((2, len(celltypes), len(pathologies)), np.nan, dtype=np.float32)
    minimum_shape = np.full_like(shape_median, np.nan)
    closed = np.zeros((2, len(celltypes), len(pathologies)), dtype=bool)
    bound_scale_used = np.full((len(celltypes), len(pathologies)), np.nan, dtype=np.float32)

    normal_pdf = PdfPages(args.out_dir / "figure_11_normal_locked_test_donor_multicontour_atlas.pdf")
    disease_pdf = PdfPages(args.out_dir / "figure_11_disease_locked_test_donor_multicontour_atlas.pdf")
    try:
        for axis, pathology in enumerate(pathologies):
            records = []
            for celltype, name in enumerate(celltypes):
                record = prepare_record(data, celltype, axis, test_order, args.grid_x, args.grid_y, args.dense_nodes, args.bound_scale)
                records.append(record)
                shape_median[0, celltype, axis] = float(np.nanmedian(record["normal_correlation"]))
                shape_median[1, celltype, axis] = float(np.nanmedian(record["disease_correlation"]))
                minimum_shape[0, celltype, axis] = float(np.nanmin(record["normal_correlation"]))
                minimum_shape[1, celltype, axis] = float(np.nanmin(record["disease_correlation"]))
                closed[:, celltype, axis] = bool(record["closed"])
                bound_scale_used[celltype, axis] = float(record["bound_scale"])
                if celltype == representative_index:
                    representative_records[axis] = record
            render_page(args.out_dir, normal_pdf, "normal", axis, pathology, celltypes, records, colors, test_order, donor_names)
            render_page(args.out_dir, disease_pdf, "disease", axis, pathology, celltypes, records, colors, test_order, donor_names)
            print(f"[multicontour] rendered {axis + 1}/{len(pathologies)} {pathology}", flush=True)
    finally:
        normal_pdf.close()
        disease_pdf.close()

    render_representative(args.out_dir, "normal", args.representative_celltype, representative_index, pathologies, representative_records, colors, test_order, donor_names)
    render_representative(args.out_dir, "disease", args.representative_celltype, representative_index, pathologies, representative_records, colors, test_order, donor_names)
    np.savez_compressed(
        args.out_dir / "donor_multicontour_shape_statistics.npz",
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(pathologies, dtype=object),
        endpoint_names=np.asarray(("normal", "disease"), dtype=object),
        locked_test_donor_names=np.asarray([donor_names[index] for index in test_order], dtype=object),
        locked_test_shape_correlation_median=shape_median,
        locked_test_shape_correlation_minimum=minimum_shape,
        all_requested_hdr_closed=closed,
        bound_scale_used=bound_scale_used,
    )
    summary = {
        "schema_version": "prism.donor_multicontour_local_fingerprint.v1",
        "input": str(args.input.resolve()),
        "checkpoint_sha256": str(data["checkpoint_sha256"].item()),
        "checkpoint_epoch": int(data["checkpoint_epoch"].item()),
        "locked_test_donors": [donor_names[index] for index in test_order],
        "test_refit": False,
        "path": "same saved common 5-D path for every donor",
        "metric": "donor-specific ordinal-Fisher metric evaluated along the saved path",
        "fingerprint": "Figure-11-style anchored mixture of local Fisher Gaussian kernels",
        "normal_shape_correlation": {
            "median_across_celltype_pathology": float(np.nanmedian(shape_median[0])),
            "minimum_celltype_pathology_median": float(np.nanmin(shape_median[0])),
        },
        "disease_shape_correlation": {
            "median_across_celltype_pathology": float(np.nanmedian(shape_median[1])),
            "minimum_celltype_pathology_median": float(np.nanmin(shape_median[1])),
        },
        "all_contours_closed": bool(np.all(closed)),
        "warning": "High similarity partly reflects the deliberately shared common path; independent biological validation still requires held-out target-cell empirical comparison and matched-versus-shuffled personal-code tests.",
    }
    (args.out_dir / "donor_multicontour_shape_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
