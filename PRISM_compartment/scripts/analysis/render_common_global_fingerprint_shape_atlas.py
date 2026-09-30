#!/usr/bin/env python3
"""Render path-free common whole-5D Riemann fingerprint shape atlases.

The detailed fingerprint is reconstructed from the saved Figure-13 whole-cube
distance and Riemann-volume fields.  Personal and personal-by-pathology terms
are absent from those saved common fields.  The atlas deliberately removes all
geodesic/Euclidean paths so that only the normal/disease LAND shape remains.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import qmc


BLUE = "#0067c5"
RED = "#d7193f"
NAVY = "#173f76"


def safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def load_figure13_module(path: Path):
    spec = importlib.util.spec_from_file_location("prism_figure13_renderer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import Figure-13 renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def correlation(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=np.float64).ravel()
    second = np.asarray(second, dtype=np.float64).ravel()
    first = first - first.mean()
    second = second - second.mean()
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.dot(first, second) / max(denominator, 1e-30))


def prepare_shapes(args, figure13):
    data = np.load(args.global_npz, allow_pickle=True)
    source_summary = json.load(args.global_json.open(encoding="utf-8"))
    slice_data = np.load(args.slice_npz, allow_pickle=True)
    donor = np.load(args.donor_summary_npz, allow_pickle=True)
    celltypes = [str(value) for value in data["celltype_names"]]
    pathologies = [str(value) for value in data["pathology_names"]]
    if celltypes != [str(value) for value in donor["celltype_names"]]:
        raise RuntimeError("donor-summary cell-type vocabulary mismatch")
    if pathologies != [str(value) for value in donor["pathology_names"]]:
        raise RuntimeError("donor-summary pathology vocabulary mismatch")

    coordinate = np.asarray(data["grid_coordinate"], dtype=np.float64)
    distance = np.asarray(data["global_distance"], dtype=np.float64)
    volume = np.asarray(data["riemann_volume"], dtype=np.float64)
    selected_lengths = np.asarray(data["selected_path_length"], dtype=np.float64)
    shared_sigma = 0.34 * float(np.median(selected_lengths))
    sobol = qmc.Sobol(d=5, scramble=True, seed=20260810)
    points = sobol.random_base2(m=args.sobol_power)
    row_lookup = {
        (row["celltype"], row["pathology"]): row for row in source_summary["rows"]
    }
    direction_lookup = {
        (int(celltype), int(axis)): row
        for row, (celltype, axis) in enumerate(
            zip(slice_data["pair_celltype"], slice_data["pair_axis"])
        )
    }

    shape = None
    normal_rel = None
    disease_rel = None
    y_limit = np.zeros((len(celltypes), len(pathologies)), dtype=np.float64)
    overlap = np.zeros_like(y_limit)
    for celltype_index, name in enumerate(celltypes):
        interpolation_grid = (coordinate,) * 5
        volume_sample = RegularGridInterpolator(
            interpolation_grid,
            volume[celltype_index].reshape((len(coordinate),) * 5),
            bounds_error=False,
            fill_value=None,
        )(points)
        distance_sample = np.stack(
            [
                RegularGridInterpolator(
                    interpolation_grid,
                    distance[celltype_index, source].reshape((len(coordinate),) * 5),
                    bounds_error=False,
                    fill_value=None,
                )(points)
                for source in range(6)
            ]
        )
        volume_weight = np.clip(volume_sample, 1e-300, None) / float(len(points))
        normal_density_5d, normal_mass = figure13.global_land(
            distance_sample[0], volume_weight, shared_sigma
        )
        log_volume = np.log10(np.clip(volume_sample, 1e-300, None))
        for axis, pathology in enumerate(pathologies):
            disease_density_5d, disease_mass = figure13.global_land(
                distance_sample[axis + 1], volume_weight, shared_sigma
            )
            pair_index = celltype_index * len(pathologies) + axis
            count = int(data["selected_path_node_count"][pair_index])
            path5 = np.asarray(data["selected_path"][pair_index, :count], dtype=np.float64)
            direction_row = direction_lookup[(celltype_index, axis)]
            direction = np.asarray(slice_data["slice_direction"][direction_row], dtype=np.float64)
            result = figure13.projection_fields(
                points,
                direction,
                axis,
                normal_mass,
                disease_mass,
                log_volume,
                path5,
            )
            if shape is None:
                grid_shape = result["surface"].shape
                shape = np.empty((len(celltypes), len(pathologies), *grid_shape), dtype=np.float32)
                normal_rel = np.empty_like(shape)
                disease_rel = np.empty_like(shape)
            shape[celltype_index, axis] = result["surface"]
            normal_rel[celltype_index, axis] = result["normal_rel"]
            disease_rel[celltype_index, axis] = result["disease_rel"]
            y_limit[celltype_index, axis] = float(result["y"][-1])
            overlap[celltype_index, axis] = float(
                np.sum(np.minimum(normal_density_5d, disease_density_5d) * volume_weight)
            )
        print(f"[shape-atlas] prepared {celltype_index + 1}/{len(celltypes)} {name}", flush=True)

    assert shape is not None and normal_rel is not None and disease_rel is not None
    feature = np.concatenate(
        (normal_rel.reshape(len(celltypes), len(pathologies), -1),
         disease_rel.reshape(len(celltypes), len(pathologies), -1)),
        axis=-1,
    )
    pairwise = np.empty((len(pathologies), len(celltypes), len(celltypes)), dtype=np.float64)
    consensus_similarity = np.empty((len(celltypes), len(pathologies)), dtype=np.float64)
    for axis in range(len(pathologies)):
        for first in range(len(celltypes)):
            for second in range(len(celltypes)):
                pairwise[axis, first, second] = correlation(feature[first, axis], feature[second, axis])
            leave_one_out = np.mean(np.delete(feature[:, axis], first, axis=0), axis=0)
            consensus_similarity[first, axis] = correlation(feature[first, axis], leave_one_out)

    return {
        "celltypes": celltypes,
        "pathologies": pathologies,
        "normal_rel": normal_rel,
        "disease_rel": disease_rel,
        "surface": shape,
        "y_limit": y_limit,
        "global_land_overlap": overlap,
        "pairwise_shape_similarity": pairwise,
        "cross_celltype_consensus_similarity": consensus_similarity,
        "donor_test_validation_similarity": np.asarray(
            donor["test_endpoint_metric_similarity"], dtype=np.float64
        ),
        "shared_sigma": shared_sigma,
    }


def draw_shape_2d(ax, normal: np.ndarray, disease: np.ndarray, title: str) -> None:
    y = np.linspace(0.0, 1.0, normal.shape[0])
    x = np.linspace(0.0, 1.0, normal.shape[1])
    xx, yy = np.meshgrid(x, y)
    fill_levels = (0.10, 0.25, 0.45, 0.65, 0.85, 1.01)
    line_levels = (0.20, 0.40, 0.60, 0.80)
    ax.contourf(xx, yy, normal, levels=fill_levels, cmap="Blues", alpha=0.66)
    ax.contourf(xx, yy, disease, levels=fill_levels, cmap="Reds", alpha=0.50)
    ax.contour(xx, yy, normal, levels=line_levels, colors=BLUE, linewidths=(0.55, 0.8, 1.1, 1.45))
    ax.contour(xx, yy, disease, levels=line_levels, colors=RED, linewidths=(0.55, 0.8, 1.1, 1.45), linestyles="--")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks((0.0, 0.5, 1.0))
    ax.set_yticks((0.0, 0.5, 1.0))
    ax.tick_params(labelsize=6, length=2)
    ax.set_title(title, fontsize=9.2, fontweight="bold", color=NAVY, pad=3)
    ax.set_facecolor("white")


def draw_shape_3d(ax, normal: np.ndarray, disease: np.ndarray, figure13, title: str) -> None:
    y = np.linspace(0.0, 1.0, normal.shape[0])
    x = np.linspace(0.0, 1.0, normal.shape[1])
    xx, yy = np.meshgrid(x, y)
    surface = 0.5 * (normal + disease)
    surface /= max(float(surface.max()), 1e-30)
    stride = 3
    ax.plot_surface(
        xx[::stride, ::stride],
        yy[::stride, ::stride],
        surface[::stride, ::stride],
        facecolors=figure13.surface_facecolors(
            normal[::stride, ::stride], disease[::stride, ::stride]
        ),
        rstride=1,
        cstride=1,
        linewidth=0.16,
        edgecolor=(0.20, 0.20, 0.20, 0.20),
        antialiased=True,
        shade=False,
        alpha=0.92,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_zlim(0.0, 1.05)
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_zticks(())
    ax.set_title(title, fontsize=8.8, fontweight="bold", color=NAVY, pad=0)
    ax.view_init(elev=29, azim=-58)
    ax.grid(False)


def render_shape_pages(result: dict, out_dir: Path, figure13) -> None:
    celltypes = result["celltypes"]
    pathologies = result["pathologies"]
    normal = result["normal_rel"]
    disease = result["disease_rel"]
    donor_similarity = result["donor_test_validation_similarity"]
    out_dir.mkdir(parents=True, exist_ok=True)
    legend = (
        Line2D((0,), (0,), color=BLUE, lw=3.0, label="정상 공통 LAND"),
        Line2D((0,), (0,), color=RED, lw=3.0, ls="--", label="질병 공통 LAND"),
    )

    pdf2d = out_dir / "figure_08a_pathology_by_celltype_common_shape_2d_atlas.pdf"
    pdf3d = out_dir / "figure_08b_pathology_by_celltype_common_shape_3d_atlas.pdf"
    with PdfPages(pdf2d) as writer2d, PdfPages(pdf3d) as writer3d:
        for axis, pathology in enumerate(pathologies):
            fig2, axes2 = plt.subplots(4, 6, figsize=(17.5, 11.8), constrained_layout=False)
            for index, ax in enumerate(axes2.ravel()):
                draw_shape_2d(ax, normal[index, axis], disease[index, axis], celltypes[index])
                if index % 6 != 0:
                    ax.set_yticklabels([])
                if index < 18:
                    ax.set_xticklabels([])
            median_donor = float(np.median(donor_similarity[:, axis]))
            fig2.suptitle(
                f"{pathology}: 개인 성분을 제거한 세포타입별 공통 Riemann 지문 모양",
                fontsize=20,
                fontweight="bold",
                color=NAVY,
                y=0.985,
            )
            fig2.text(
                0.5,
                0.948,
                f"경로와 metric 배경을 제거하고 정상·질병 LAND 윤곽만 표시 | locked test↔validation donor 유사도 중앙값={median_donor:.6f}",
                ha="center",
                fontsize=11,
                color="#53637a",
            )
            fig2.legend(handles=legend, loc="lower center", ncol=2, frameon=False, fontsize=10)
            fig2.text(0.5, 0.030, "가로=병리 심도 | 세로=정규화된 dominant off-axis | 세포타입 간 세로 눈금은 모양 비교용으로 정규화", ha="center", fontsize=9, color="#667085")
            fig2.subplots_adjust(left=0.045, right=0.985, top=0.905, bottom=0.085, wspace=0.11, hspace=0.24)
            path2 = out_dir / f"figure_08a_{axis + 1:02d}_{safe_name(pathology)}_common_shape_2d.png"
            fig2.savefig(path2, dpi=210, bbox_inches="tight")
            writer2d.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)

            fig3 = plt.figure(figsize=(18.2, 12.2))
            for index, celltype in enumerate(celltypes):
                ax = fig3.add_subplot(4, 6, index + 1, projection="3d")
                draw_shape_3d(ax, normal[index, axis], disease[index, axis], figure13, celltype)
            fig3.suptitle(
                f"{pathology}: 세포타입별 공통 Riemann 지문의 3D scalar lift",
                fontsize=20,
                fontweight="bold",
                color=NAVY,
                y=0.982,
            )
            fig3.text(
                0.5,
                0.947,
                "파랑=정상 공통 LAND | 빨강=질병 공통 LAND | 경로를 제거하여 지형의 shape만 비교",
                ha="center",
                fontsize=11,
                color="#53637a",
            )
            fig3.text(0.5, 0.022, "높이는 LAND 확률밀도의 scalar lift이며 5차원 다양체의 등거리 embedding이 아니다.", ha="center", fontsize=9, color="#667085")
            fig3.subplots_adjust(left=0.02, right=0.985, top=0.91, bottom=0.05, wspace=0.02, hspace=0.03)
            path3 = out_dir / f"figure_08b_{axis + 1:02d}_{safe_name(pathology)}_common_shape_3d.png"
            fig3.savefig(path3, dpi=205, bbox_inches="tight")
            writer3d.savefig(fig3, bbox_inches="tight")
            plt.close(fig3)
            print(f"[shape-atlas] rendered {pathology}", flush=True)


def render_similarity_summary(result: dict, out_dir: Path) -> None:
    donor = result["donor_test_validation_similarity"]
    cross = result["cross_celltype_consensus_similarity"]
    celltypes = result["celltypes"]
    pathologies = result["pathologies"]
    fig, axes = plt.subplots(1, 2, figsize=(15.8, 10.5), gridspec_kw={"width_ratios": (1.0, 1.0)})
    donor_error = 1.0 - donor
    image0 = axes[0].imshow(donor_error, aspect="auto", cmap="magma_r")
    axes[0].set_title("A. 같은 세포타입·같은 병리: donor가 바뀐 차이", fontweight="bold", color=NAVY)
    axes[0].set_xticks(range(len(pathologies)), pathologies)
    axes[0].set_yticks(range(len(celltypes)), celltypes)
    axes[0].set_xlabel("값이 0에 가까울수록 validation과 test shape가 같음")
    fig.colorbar(image0, ax=axes[0], fraction=0.047, pad=0.03, label="1 - validation/test correlation")

    image1 = axes[1].imshow(cross, aspect="auto", cmap="coolwarm", vmin=-1.0, vmax=1.0)
    axes[1].set_title("B. 서로 다른 세포타입: 병리축별 평균 shape와의 상관", fontweight="bold", color=NAVY)
    axes[1].set_xticks(range(len(pathologies)), pathologies)
    axes[1].set_yticks(range(len(celltypes)), celltypes)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("높을수록 여러 세포타입이 비슷한 모양")
    fig.colorbar(image1, ax=axes[1], fraction=0.047, pad=0.03, label="leave-one-celltype-out shape correlation")
    fig.suptitle(
        "공통 지문은 donor 사이에서는 매우 안정적이지만, 세포타입 사이에서는 동일한 모양일 필요가 없다",
        fontsize=18,
        fontweight="bold",
        color=NAVY,
        y=0.995,
    )
    fig.text(
        0.5,
        0.018,
        "왼쪽은 frozen validation endpoint-metric template과 locked-test donor의 유사도이며, 오른쪽은 전역 LAND 정상·질병 shape의 세포타입 간 비교이다.",
        ha="center",
        fontsize=9.5,
        color="#667085",
    )
    fig.subplots_adjust(left=0.13, right=0.96, top=0.94, bottom=0.065, wspace=0.24)
    fig.savefig(out_dir / "figure_08c_donor_vs_celltype_shape_similarity.png", dpi=210, bbox_inches="tight")
    fig.savefig(out_dir / "figure_08c_donor_vs_celltype_shape_similarity.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure13-renderer", required=True, type=Path)
    parser.add_argument("--global-npz", required=True, type=Path)
    parser.add_argument("--global-json", required=True, type=Path)
    parser.add_argument("--slice-npz", required=True, type=Path)
    parser.add_argument("--donor-summary-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sobol-power", type=int, default=19)
    parser.add_argument("--reuse-prepared", action="store_true")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    figure13 = load_figure13_module(args.figure13_renderer)
    figure13.set_style()
    korean_font_path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if korean_font_path.exists():
        font_manager.fontManager.addfont(korean_font_path)
        korean_family = font_manager.FontProperties(fname=korean_font_path).get_name()
    else:
        korean_family = "DejaVu Sans"
    plt.rcParams.update({"font.family": korean_family, "axes.unicode_minus": False})
    prepared_path = args.out_dir / "common_global_fingerprint_shape_arrays.npz"
    if args.reuse_prepared and prepared_path.exists():
        prepared = np.load(prepared_path, allow_pickle=True)
        result = {
            "celltypes": [str(value) for value in prepared["celltype_names"]],
            "pathologies": [str(value) for value in prepared["pathology_names"]],
            "normal_rel": np.asarray(prepared["normal_rel"]),
            "disease_rel": np.asarray(prepared["disease_rel"]),
            "surface": np.asarray(prepared["surface"]),
            "y_limit": np.asarray(prepared["y_limit"]),
            "global_land_overlap": np.asarray(prepared["global_land_overlap"]),
            "donor_test_validation_similarity": np.asarray(prepared["donor_test_validation_similarity"]),
            "cross_celltype_consensus_similarity": np.asarray(prepared["cross_celltype_consensus_similarity"]),
            "pairwise_shape_similarity": np.asarray(prepared["pairwise_shape_similarity"]),
            "shared_sigma": float(prepared["shared_fisher_sigma"]),
        }
        print(f"[shape-atlas] reused {prepared_path}", flush=True)
    else:
        result = prepare_shapes(args, figure13)
    np.savez_compressed(
        prepared_path,
        schema_version=np.asarray("prism.common_global_fingerprint_shape.v1", dtype=object),
        celltype_names=np.asarray(result["celltypes"], dtype=object),
        pathology_names=np.asarray(result["pathologies"], dtype=object),
        normal_rel=result["normal_rel"],
        disease_rel=result["disease_rel"],
        surface=result["surface"],
        y_limit=result["y_limit"],
        global_land_overlap=result["global_land_overlap"],
        donor_test_validation_similarity=result["donor_test_validation_similarity"],
        cross_celltype_consensus_similarity=result["cross_celltype_consensus_similarity"],
        pairwise_shape_similarity=result["pairwise_shape_similarity"],
        shared_fisher_sigma=np.asarray(result["shared_sigma"]),
    )
    render_shape_pages(result, args.out_dir, figure13)
    render_similarity_summary(result, args.out_dir)
    pairwise = result["pairwise_shape_similarity"]
    off_diagonal = np.concatenate(
        [pairwise[axis][~np.eye(len(result["celltypes"]), dtype=bool)] for axis in range(len(result["pathologies"]))]
    )
    summary = {
        "schema_version": "prism.common_global_fingerprint_shape.summary.v1",
        "celltypes": result["celltypes"],
        "pathologies": result["pathologies"],
        "personal_terms": "off in the saved common Figure-13 fields",
        "donor_test_validation_similarity_median": float(np.median(result["donor_test_validation_similarity"])),
        "donor_test_validation_similarity_range": [
            float(np.min(result["donor_test_validation_similarity"])),
            float(np.max(result["donor_test_validation_similarity"])),
        ],
        "cross_celltype_consensus_similarity_median": float(np.median(result["cross_celltype_consensus_similarity"])),
        "cross_celltype_consensus_similarity_range": [
            float(np.min(result["cross_celltype_consensus_similarity"])),
            float(np.max(result["cross_celltype_consensus_similarity"])),
        ],
        "pairwise_cross_celltype_shape_similarity_median": float(np.median(off_diagonal)),
        "pairwise_cross_celltype_shape_similarity_q10_q90": [
            float(np.percentile(off_diagonal, 10)),
            float(np.percentile(off_diagonal, 90)),
        ],
        "shape_definition": "concatenated normal and disease whole-5D LAND relative-density maps in normalized named/off-axis display coordinates",
        "donor_similarity_definition": "locked-test endpoint Fisher metric fingerprint versus frozen validation template",
        "warning": "donor stability within a celltype-pathology pair does not imply that different celltypes share the same fingerprint shape",
    }
    (args.out_dir / "common_global_fingerprint_shape_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[shape-atlas] wrote {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
