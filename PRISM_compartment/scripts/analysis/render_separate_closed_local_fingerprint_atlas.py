#!/usr/bin/env python3
"""Render normal and disease common local-Fisher fingerprints separately.

Unlike the whole-cube Figure-13 LAND marginal, these endpoint-anchored local
Fisher kernels are evaluated beyond the pathology-cube boundary solely in the
local display chart.  This makes every requested HDR contour closed without
pretending that out-of-domain values are a global LAND estimate.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
import numpy as np


BLUE = "#0067c5"
RED = "#d7193f"
NAVY = "#173f76"
HDR_MASSES = (0.95, 0.88, 0.78, 0.66, 0.53, 0.40, 0.27)


def safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def load_module(path: Path):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("prism_local_fingerprint", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import local-fingerprint renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_records(args, local):
    metric = np.load(args.metric_arrays, allow_pickle=True)
    common = np.load(args.common_arrays, allow_pickle=True)
    celltypes = [str(value) for value in metric["celltype_names"]]
    pathologies = [str(value) for value in metric["pathology_names"]]
    metric_lookup = local.pair_lookup(metric)
    common_lookup = local.pair_lookup(common)
    records: dict[tuple[int, int], dict] = {}
    for celltype, name in enumerate(celltypes):
        for axis, pathology in enumerate(pathologies):
            index = metric_lookup[(celltype, axis)]
            if index != common_lookup[(celltype, axis)]:
                raise RuntimeError("metric/common pair ordering differs")
            result = local.fingerprint_fields(
                metric,
                common,
                index,
                nx=args.grid_x,
                ny=args.grid_y,
                bound_scale=args.bound_scale,
            )
            for endpoint, key, anchor_x in (
                ("normal", "density_normal", 0.0),
                ("disease", "density_disease", 1.0),
            ):
                density = np.asarray(result[key], dtype=np.float64)
                density /= max(float(density.max()), 1e-300)
                levels = sorted(
                    {float(local.hdr_threshold(density, mass)) for mass in HDR_MASSES}
                )
                outer = levels[0]
                mask = density >= outer
                boundary = np.concatenate(
                    (mask[0], mask[-1], mask[:, 0], mask[:, -1])
                )
                records[(celltype, axis, endpoint)] = {
                    "density": density,
                    "levels": levels,
                    "x": np.asarray(result["grid_x"], dtype=np.float64) - anchor_x,
                    "y": np.asarray(result["grid_y"], dtype=np.float64),
                    "closed": not bool(np.any(boundary)),
                    "outer": outer,
                    "anchor_x": anchor_x,
                }
        print(f"[closed-atlas] prepared {celltype + 1}/{len(celltypes)} {name}", flush=True)
    return celltypes, pathologies, records


def shared_bounds(records, celltype_count: int, axis: int, endpoint: str):
    x_min, x_max, y_min, y_max = np.inf, -np.inf, np.inf, -np.inf
    for celltype in range(celltype_count):
        record = records[(celltype, axis, endpoint)]
        mask = record["density"] >= record["outer"]
        rows, columns = np.where(mask)
        x_min = min(x_min, float(record["x"][columns].min()))
        x_max = max(x_max, float(record["x"][columns].max()))
        y_min = min(y_min, float(record["y"][rows].min()))
        y_max = max(y_max, float(record["y"][rows].max()))
    x_pad = 0.10 * max(x_max - x_min, 1e-6)
    y_pad = 0.10 * max(y_max - y_min, 1e-6)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def draw_2d(ax, record: dict, endpoint: str, title: str, bounds) -> None:
    xx, yy = np.meshgrid(record["x"], record["y"], indexing="xy")
    color = BLUE if endpoint == "normal" else RED
    cmap = "Blues" if endpoint == "normal" else "Reds"
    ax.contourf(xx, yy, record["density"], levels=[record["outer"], 0.2, 0.4, 0.6, 0.8, 1.01], cmap=cmap, alpha=0.88)
    ax.contour(
        xx,
        yy,
        record["density"],
        levels=record["levels"],
        colors=color,
        linewidths=np.linspace(0.65, 1.8, len(record["levels"])),
        linestyles="-" if endpoint == "normal" else "--",
    )
    ax.scatter((0.0,), (0.0,), s=22, color=color, edgecolor="white", linewidth=0.7, zorder=10)
    ax.axhline(0.0, color="#a8b0bd", lw=0.45, zorder=0)
    ax.axvline(0.0, color="#a8b0bd", lw=0.45, zorder=0)
    ax.set_xlim(*bounds[0])
    ax.set_ylim(*bounds[1])
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(title, fontsize=9.2, fontweight="bold", color=NAVY, pad=3)
    for spine in ax.spines.values():
        spine.set_color("#98a2b3")
        spine.set_linewidth(0.55)


def draw_3d(ax, record: dict, endpoint: str, title: str, bounds) -> None:
    xx, yy = np.meshgrid(record["x"], record["y"], indexing="xy")
    density = record["density"]
    base_color = np.asarray(to_rgb(BLUE if endpoint == "normal" else RED))
    strength = density[..., None]
    facecolors = 1.0 - (1.0 - base_color) * (0.10 + 0.90 * strength)
    stride = 3
    ax.plot_surface(
        xx[::stride, ::stride],
        yy[::stride, ::stride],
        density[::stride, ::stride],
        facecolors=facecolors[::stride, ::stride],
        rstride=1,
        cstride=1,
        linewidth=0.15,
        edgecolor=(0.20, 0.20, 0.20, 0.20),
        antialiased=True,
        shade=False,
        alpha=0.94,
    )
    ax.set_xlim(*bounds[0])
    ax.set_ylim(*bounds[1])
    ax.set_zlim(0.0, 1.05)
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_zticks(())
    ax.set_title(title, fontsize=8.8, fontweight="bold", color=NAVY, pad=0)
    ax.view_init(elev=30, azim=-60)
    ax.grid(False)


def render_endpoint(endpoint: str, celltypes, pathologies, records, out_dir: Path) -> None:
    korean = "정상" if endpoint == "normal" else "질병"
    color = BLUE if endpoint == "normal" else RED
    pdf2 = out_dir / f"figure_09_{endpoint}_closed_local_fingerprint_2d_atlas.pdf"
    pdf3 = out_dir / f"figure_09_{endpoint}_closed_local_fingerprint_3d_atlas.pdf"
    with PdfPages(pdf2) as writer2, PdfPages(pdf3) as writer3:
        for axis, pathology in enumerate(pathologies):
            bounds = shared_bounds(records, len(celltypes), axis, endpoint)
            fig2, axes = plt.subplots(4, 6, figsize=(17.5, 11.8))
            for celltype, ax in enumerate(axes.ravel()):
                draw_2d(ax, records[(celltype, axis, endpoint)], endpoint, celltypes[celltype], bounds)
            fig2.suptitle(
                f"{pathology}: 개인 성분을 제거한 세포타입별 {korean} local-Fisher 지문",
                fontsize=20,
                fontweight="bold",
                color=NAVY,
                y=0.985,
            )
            fig2.text(
                0.5,
                0.946,
                "끝점을 중심으로 좌표를 이동하고 모든 HDR 등고선이 닫히도록 local chart를 확장 | 세포타입 간 동일 축척",
                ha="center",
                fontsize=10.5,
                color="#53637a",
            )
            fig2.legend(
                handles=(Line2D((0,), (0,), color=color, lw=3, ls="-" if endpoint == "normal" else "--", label=f"{korean} common local-Fisher HDR"),),
                loc="lower center",
                frameon=False,
                fontsize=10,
            )
            fig2.text(0.5, 0.027, "이 닫힌 윤곽은 끝점 주변의 local-Fisher neighbourhood이며 전역 5차원 LAND의 경계 밖 확률질량이 아니다.", ha="center", fontsize=9, color="#667085")
            fig2.subplots_adjust(left=0.025, right=0.985, top=0.90, bottom=0.075, wspace=0.10, hspace=0.24)
            image2 = out_dir / f"figure_09_{endpoint}_{axis + 1:02d}_{safe_name(pathology)}_closed_2d.png"
            fig2.savefig(image2, dpi=210, bbox_inches="tight")
            writer2.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)

            fig3 = plt.figure(figsize=(18.2, 12.2))
            for celltype, name in enumerate(celltypes):
                ax = fig3.add_subplot(4, 6, celltype + 1, projection="3d")
                draw_3d(ax, records[(celltype, axis, endpoint)], endpoint, name, bounds)
            fig3.suptitle(
                f"{pathology}: 세포타입별 {korean} local-Fisher 지문의 3D scalar lift",
                fontsize=20,
                fontweight="bold",
                color=NAVY,
                y=0.982,
            )
            fig3.text(0.5, 0.946, f"{korean} 지문만 분리 | 경로·반대편 지문·전역 metric 배경 제거", ha="center", fontsize=10.5, color="#53637a")
            fig3.text(0.5, 0.022, "높이는 local-Fisher display density이며 원래 5차원 다양체의 등거리 embedding이 아니다.", ha="center", fontsize=9, color="#667085")
            fig3.subplots_adjust(left=0.015, right=0.99, top=0.905, bottom=0.045, wspace=0.01, hspace=0.02)
            image3 = out_dir / f"figure_09_{endpoint}_{axis + 1:02d}_{safe_name(pathology)}_closed_3d.png"
            fig3.savefig(image3, dpi=205, bbox_inches="tight")
            writer3.savefig(fig3, bbox_inches="tight")
            plt.close(fig3)
            print(f"[closed-atlas] rendered {endpoint} {pathology}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--local-renderer", required=True, type=Path)
    parser.add_argument("--metric-arrays", required=True, type=Path)
    parser.add_argument("--common-arrays", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--grid-x", type=int, default=175)
    parser.add_argument("--grid-y", type=int, default=135)
    parser.add_argument("--bound-scale", type=float, default=1.65)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    font_path = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    if font_path.exists():
        font_manager.fontManager.addfont(font_path)
        family = font_manager.FontProperties(fname=font_path).get_name()
    else:
        family = "DejaVu Sans"
    plt.rcParams.update({"font.family": family, "axes.unicode_minus": False, "pdf.fonttype": 42})
    local = load_module(args.local_renderer)
    celltypes, pathologies, records = prepare_records(args, local)
    closed = {
        endpoint: np.asarray(
            [[records[(celltype, axis, endpoint)]["closed"] for axis in range(len(pathologies))] for celltype in range(len(celltypes))],
            dtype=bool,
        )
        for endpoint in ("normal", "disease")
    }
    if not all(bool(np.all(value)) for value in closed.values()):
        failures = {
            endpoint: [
                (celltypes[celltype], pathologies[axis])
                for celltype, axis in zip(*np.where(~value))
            ]
            for endpoint, value in closed.items()
        }
        raise RuntimeError(f"outer HDR contour touched display boundary: {failures}")
    render_endpoint("normal", celltypes, pathologies, records, args.out_dir)
    render_endpoint("disease", celltypes, pathologies, records, args.out_dir)
    summary = {
        "schema_version": "prism.separate_closed_local_fingerprint_atlas.v1",
        "celltypes": celltypes,
        "pathologies": pathologies,
        "normal_panel_count": len(celltypes) * len(pathologies),
        "disease_panel_count": len(celltypes) * len(pathologies),
        "all_requested_hdr_contours_closed": True,
        "bound_scale": args.bound_scale,
        "personal_terms": "off in the saved common Figure-11 fields",
        "interpretation": "endpoint-anchored local-Fisher neighbourhood shape; not whole-cube global LAND outside-domain mass",
        "warning": "Use Figure 13 for global N-to-D geometry and this atlas only for fully closed local endpoint fingerprints.",
    }
    (args.out_dir / "separate_closed_local_fingerprint_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[closed-atlas] wrote {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
