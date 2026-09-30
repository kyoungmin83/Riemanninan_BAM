#!/usr/bin/env python3
"""Donor-conditioned common-only PRISM Riemann fingerprint audit.

The explicit personal baseline and personal-by-pathology response are zeroed in
the primary audit.  Sex and technology scores are also zeroed.  Donor-specific
clean latent centroids remain only as decoder conditioning contexts, allowing a
non-circular test of whether the same explicit common pathology map induces a
similar ordinal-Fisher geometry across validation and locked-test donors.

An exploratory secondary audit adds only the saved personal-response slope.
No parameter is fitted and no checkpoint is changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import torch
import torch.nn.functional as F
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disease_program_eval import build_pooled_dataset
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


PATHOLOGY = ("Thal", "Braak", "CERAD", "LATE", "Lewy")
BLUE = "#1f62c4"
RED = "#d62f4b"
GREEN = "#009b72"
ORANGE = "#ed6a24"
PURPLE = "#762a83"


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hinge(value: torch.Tensor) -> torch.Tensor:
    return torch.stack(
        (value, F.relu(value - 1.0 / 3.0), F.relu(value - 2.0 / 3.0)), dim=-1
    )


def hinge_derivative(value: torch.Tensor) -> torch.Tensor:
    return torch.stack(
        (
            torch.ones_like(value),
            (value > 1.0 / 3.0).to(value.dtype),
            (value > 2.0 / 3.0).to(value.dtype),
        ),
        dim=-1,
    )


def fisher_info(decoder, score: torch.Tensor, thresholds: torch.Tensor) -> torch.Tensor:
    cumulative, probability = decoder._score_to_probs(score, thresholds)
    sd = cumulative * (1.0 - cumulative)
    derivative = torch.cat(
        (-sd[..., :1], sd[..., :-1] - sd[..., 1:], sd[..., -1:]), dim=-1
    )
    return (derivative.square() / probability.clamp_min(1e-8)).sum(dim=-1)


def regularize_metric(metric: np.ndarray) -> np.ndarray:
    metric = 0.5 * (metric + np.swapaxes(metric, -1, -2))
    trace = np.trace(metric, axis1=-2, axis2=-1)
    scale = max(float(np.nanmedian(trace[trace > 0])) if np.any(trace > 0) else 1.0, 1e-12)
    ridge = 1e-7 * np.maximum(trace / 2.0, scale * 1e-8)
    return metric + ridge[..., None, None] * np.eye(2)


def metric_graph(metric: np.ndarray, x: np.ndarray, y: np.ndarray) -> csr_matrix:
    ny, nx = metric.shape[:2]
    node = np.arange(ny * nx, dtype=np.int64).reshape(ny, nx)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for oy, ox in ((0, 1), (1, 0), (1, 1), (1, -1)):
        y0 = slice(0, ny - oy)
        y1 = slice(oy, ny)
        if ox >= 0:
            x0, x1 = slice(0, nx - ox), slice(ox, nx)
        else:
            x0, x1 = slice(-ox, nx), slice(0, nx + ox)
        source = node[y0, x0].ravel()
        target = node[y1, x1].ravel()
        local = 0.5 * (metric[y0, x0] + metric[y1, x1])
        delta = np.asarray((float(x[abs(ox)] - x[0]) * np.sign(ox), float(y[oy] - y[0])))
        square = np.einsum("i,...ij,j->...", delta, local, delta)
        weight = np.sqrt(np.maximum(square, 1e-30)).ravel()
        rows.extend((source, target))
        cols.extend((target, source))
        weights.extend((weight, weight))
    return coo_matrix(
        (np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))),
        shape=(ny * nx, ny * nx),
    ).tocsr()


def reconstruct_path(predecessor: np.ndarray, source: int, target: int, points: np.ndarray) -> np.ndarray:
    nodes = [int(target)]
    cursor = int(target)
    while cursor != source and len(nodes) <= len(predecessor) + 1:
        cursor = int(predecessor[cursor])
        if cursor < 0:
            return np.vstack((points[source], points[target]))
        nodes.append(cursor)
    nodes.reverse()
    return points[np.asarray(nodes, dtype=np.int64)]


def path_length_axis(metric: np.ndarray, x: np.ndarray) -> float:
    total = 0.0
    for index in range(len(x) - 1):
        local = 0.5 * (metric[0, index] + metric[0, index + 1])
        delta = np.asarray((x[index + 1] - x[index], 0.0))
        total += math.sqrt(max(float(delta @ local @ delta), 1e-30))
    return total


def land_density(distance: np.ndarray, volume: np.ndarray, sigma: float, cell_area: float) -> np.ndarray:
    kernel = np.exp(-0.5 * np.square(distance / max(sigma, 1e-12)))
    normalizer = float(np.sum(kernel * volume) * cell_area)
    return kernel / max(normalizer, 1e-300)


def fingerprint_vector(normal: np.ndarray, disease: np.ndarray) -> np.ndarray:
    values = np.concatenate((normal.ravel(), disease.ravel()))
    values = np.log(np.maximum(values / max(float(values.max()), 1e-300), 1e-8))
    values -= values.mean()
    scale = float(np.linalg.norm(values))
    return values / max(scale, 1e-30)


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=np.float64).ravel()
    right = np.asarray(right, dtype=np.float64).ravel()
    left = left - left.mean()
    right = right - right.mean()
    den = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(left @ right / den) if den > 1e-30 else float("nan")


def pairwise_median(vectors: np.ndarray) -> float:
    values = [correlation(vectors[i], vectors[j]) for i in range(len(vectors)) for j in range(i + 1, len(vectors))]
    return float(np.nanmedian(values)) if values else float("nan")


def resample_path(path: np.ndarray, size: int = 60) -> np.ndarray:
    delta = np.diff(path, axis=0)
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(delta, axis=1))))
    if arc[-1] <= 1e-12:
        return np.repeat(path[:1], size, axis=0)
    target = np.linspace(0.0, arc[-1], size)
    return np.column_stack([np.interp(target, arc, path[:, dim]) for dim in range(path.shape[1])])


def bootstrap_ci(values: np.ndarray, seed: int, n_boot: int = 5000) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(n_boot, len(values)), replace=True)
    stat = np.median(draws, axis=1)
    return float(np.median(values)), float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5))


@torch.no_grad()
def donor_base_scores(system, z_values: np.ndarray, celltype: int, device: torch.device) -> np.ndarray:
    valid = np.all(np.isfinite(z_values), axis=1)
    output = np.full((len(z_values), system.decoder.n_genes), np.nan, dtype=np.float32)
    if not np.any(valid):
        return output
    z = torch.as_tensor(z_values[valid], dtype=torch.float32, device=device)
    ct = torch.full((len(z),), celltype, dtype=torch.long, device=device)
    base_score = system.decoder.celltype_baseline(ct)
    state_score = system.decoder.state_score_from_z_and_base(
        z,
        base_score,
        celltype_id=ct,
        generator_gate=system._generator_gate_for_forward(),
    )
    zero_tech = torch.zeros_like(base_score)
    if system.decoder.score_mixer is None:
        score = base_score + state_score
    else:
        score, _ = system.decoder.score_mixer(z, base_score, zero_tech, state_score)
    output[valid] = score.float().cpu().numpy()
    return output


@torch.no_grad()
def evaluate_metric_grid(
    decoder,
    thresholds: torch.Tensor,
    base_score: np.ndarray,
    score_dictionary: torch.Tensor,
    chart: torch.Tensor,
    pathology_points: torch.Tensor,
    response_score: torch.Tensor | None,
    batch_size: int,
) -> np.ndarray:
    n_point = len(pathology_points)
    output = np.empty((n_point, 2, 2), dtype=np.float64)
    base_tensor = torch.as_tensor(base_score, dtype=torch.float32, device=pathology_points.device)
    for start in range(0, n_point, batch_size):
        stop = min(start + batch_size, n_point)
        point = pathology_points[start:stop]
        feature = hinge(point)
        score = base_tensor.unsqueeze(0) + torch.einsum("bkf,kfg->bg", feature, score_dictionary)
        derivative = hinge_derivative(point)
        jacobian5 = torch.einsum("bkf,kfg->bkg", derivative, score_dictionary)
        if response_score is not None:
            score = score + point @ response_score
            jacobian5 = jacobian5 + response_score.unsqueeze(0)
        jacobian2 = torch.einsum("ka,bkg->bag", chart, jacobian5)
        info = fisher_info(decoder, score, thresholds)
        metric = torch.einsum("bag,bg,bcg->bac", jacobian2, info, jacobian2)
        output[start:stop] = (metric / float(decoder.n_genes)).double().cpu().numpy()
    return output


def analyse_surface(metric_flat: np.ndarray, x: np.ndarray, y: np.ndarray, sigma: float) -> dict[str, Any]:
    metric = regularize_metric(metric_flat.reshape(len(y), len(x), 2, 2))
    graph = metric_graph(metric, x, y)
    source = 0
    target = len(x) - 1
    distance, predecessor = dijkstra(graph, directed=False, indices=(source, target), return_predecessors=True)
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    path = reconstruct_path(predecessor[0], source, target, points)
    sign, logdet = np.linalg.slogdet(metric)
    if not np.all(sign > 0):
        raise RuntimeError("non-positive regularized chart metric")
    volume = np.exp(0.5 * logdet)
    area = float((x[1] - x[0]) * (y[1] - y[0]))
    normal = land_density(distance[0].reshape(len(y), len(x)), volume, sigma, area)
    disease = land_density(distance[1].reshape(len(y), len(x)), volume, sigma, area)
    axis_length = path_length_axis(metric, x)
    geodesic_length = float(distance[0, target])
    return {
        "metric": metric,
        "normal": normal,
        "disease": disease,
        "path": path,
        "axis_length": axis_length,
        "geodesic_length": geodesic_length,
        "shortening": 100.0 * (axis_length - geodesic_length) / max(axis_length, 1e-30),
        "off_axis": float(np.max(path[:, 1])),
    }


def template_density(values: np.ndarray) -> np.ndarray:
    result = np.nanmean(values, axis=0)
    return result / max(float(np.nanmax(result)), 1e-300)


def surface_height(normal: np.ndarray, disease: np.ndarray) -> np.ndarray:
    n = normal / max(float(normal.max()), 1e-300)
    d = disease / max(float(disease.max()), 1e-300)
    return 0.5 * (n + d)


def path_height(x: np.ndarray, y: np.ndarray, height: np.ndarray, path: np.ndarray) -> np.ndarray:
    interp = RegularGridInterpolator((y, x), height, bounds_error=False, fill_value=0.0)
    return interp(np.column_stack((path[:, 1], path[:, 0]))) + 0.04


def render_common(
    out: Path,
    celltype_name: str,
    x: np.ndarray,
    y: np.ndarray,
    normal_template: np.ndarray,
    disease_template: np.ndarray,
    path_template: np.ndarray,
    test_paths: np.ndarray,
    rows: list[dict[str, Any]],
) -> None:
    fig = plt.figure(figsize=(21, 10.8), facecolor="white")
    grid = fig.add_gridspec(2, 5, hspace=0.28, wspace=0.18, top=0.84, bottom=0.10)
    xx, yy = np.meshgrid(x, y)
    for axis, name in enumerate(PATHOLOGY):
        ax = fig.add_subplot(grid[0, axis])
        n = normal_template[axis]
        d = disease_template[axis]
        levels_n = np.linspace(0.18, 0.92, 6) * float(n.max())
        levels_d = np.linspace(0.18, 0.92, 6) * float(d.max())
        ax.contour(xx, yy, n, levels=levels_n, colors=BLUE, linewidths=1.5)
        ax.contour(xx, yy, d, levels=levels_d, colors=RED, linewidths=1.5, linestyles="--")
        for path in test_paths[:, axis]:
            ax.plot(path[:, 0], path[:, 1], color="#7851a9", alpha=0.22, lw=1.0)
        ax.plot(path_template[axis, :, 0], path_template[axis, :, 1], color=GREEN, lw=3.0)
        ax.plot((0, 1), (0, 0), color=ORANGE, lw=2.0, ls="--")
        ax.scatter((0, 1), (0, 0), c=(BLUE, RED), s=26, zorder=8)
        one = rows[axis]
        ax.set_title(
            f"{name}\ntest→val shape ρ={one['test_template_shape_median']:.2f}\n"
            f"split-half ρ={one['split_half_median']:.2f}",
            fontsize=11,
            fontweight="bold",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, y[-1])
        ax.set_xlabel("named pathology severity")
        if axis == 0:
            ax.set_ylabel("fixed off-axis coordinate")
        ax.grid(alpha=0.15)

        ax3 = fig.add_subplot(grid[1, axis], projection="3d")
        height = surface_height(n, d)
        color_mix = d / np.maximum(n + d, 1e-30)
        face = plt.get_cmap("coolwarm")(color_mix)
        face[..., 3] = 0.82
        ax3.plot_surface(xx, yy, height, facecolors=face, rstride=1, cstride=1, linewidth=0, antialiased=True)
        path = path_template[axis]
        ax3.plot(path[:, 0], path[:, 1], path_height(x, y, height, path), color=GREEN, lw=3.2)
        direct = np.column_stack((np.linspace(0, 1, 80), np.zeros(80)))
        ax3.plot(direct[:, 0], direct[:, 1], path_height(x, y, height, direct), color=ORANGE, lw=2.0, ls="--")
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, y[-1])
        ax3.set_zlim(0, 1.15)
        ax3.view_init(elev=27, azim=-58)
        ax3.set_xlabel("severity", labelpad=-1)
        if axis == 0:
            ax3.set_ylabel("off-axis", labelpad=-1)
        ax3.set_zlabel("LAND", labelpad=-2)
        ax3.tick_params(labelsize=7, pad=0)
        ax3.set_title(f"test shortening={one['test_shortening_median']:.2f}%", fontsize=10, pad=0)

    fig.suptitle(
        f"{celltype_name}: the validation common-only Riemann fingerprint was frozen and tested in 9 donors",
        fontsize=22,
        fontweight="bold",
        color="#173f76",
        y=0.965,
    )
    fig.text(
        0.5,
        0.905,
        "personal and personal×pathology terms = 0 | sex/technology score = 0 | donor clean-state centroid retained only as decoder context",
        ha="center",
        fontsize=12.5,
        color="#4d5c73",
    )
    fig.legend(
        handles=(
            Line2D([0], [0], color=BLUE, lw=2, label="normal LAND"),
            Line2D([0], [0], color=RED, lw=2, ls="--", label="disease LAND"),
            Line2D([0], [0], color=GREEN, lw=3, label="validation-template Fisher path"),
            Line2D([0], [0], color=PURPLE, lw=1.5, alpha=0.5, label="locked-test donor paths"),
            Line2D([0], [0], color=ORANGE, lw=2, ls="--", label="coordinate-straight path"),
        ),
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=10.5,
        bbox_to_anchor=(0.5, 0.015),
    )
    fig.savefig(out.with_suffix(".png"), dpi=190, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def render_response(
    out: Path,
    celltype_name: str,
    x: np.ndarray,
    y: np.ndarray,
    common_normal: np.ndarray,
    common_disease: np.ndarray,
    response_normal: np.ndarray,
    response_disease: np.ndarray,
    common_path: np.ndarray,
    response_path: np.ndarray,
    donor_names: np.ndarray,
    representative_index: np.ndarray,
    response_rows: list[dict[str, Any]],
) -> None:
    fig = plt.figure(figsize=(21, 10.8), facecolor="white")
    grid = fig.add_gridspec(2, 5, hspace=0.26, wspace=0.18, top=0.84, bottom=0.10)
    xx, yy = np.meshgrid(x, y)
    for axis, name in enumerate(PATHOLOGY):
        donor = int(representative_index[axis])
        n0, d0 = common_normal[donor, axis], common_disease[donor, axis]
        nr, dr = response_normal[donor, axis], response_disease[donor, axis]
        ax = fig.add_subplot(grid[0, axis])
        ax.contour(xx, yy, d0, levels=np.linspace(0.2, 0.9, 6) * d0.max(), colors=GREEN, linewidths=1.8)
        ax.contour(xx, yy, dr, levels=np.linspace(0.2, 0.9, 6) * dr.max(), colors=PURPLE, linewidths=1.8, linestyles="--")
        ax.plot(common_path[donor, axis, :, 0], common_path[donor, axis, :, 1], color=GREEN, lw=3)
        ax.plot(response_path[donor, axis, :, 0], response_path[donor, axis, :, 1], color=PURPLE, lw=3)
        ax.plot((0, 1), (0, 0), color=ORANGE, lw=1.7, ls="--")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, y[-1])
        ax.grid(alpha=0.15)
        ax.set_xlabel("named pathology severity")
        if axis == 0:
            ax.set_ylabel("fixed off-axis coordinate")
        ax.set_title(
            f"{name} | {donor_names[donor]}\nresponse shape Δ={response_rows[axis]['representative_shape_change']:.2f}",
            fontsize=10.5,
            fontweight="bold",
        )

        ax3 = fig.add_subplot(grid[1, axis], projection="3d")
        h0 = surface_height(n0, d0)
        hr = surface_height(nr, dr)
        ax3.plot_wireframe(xx, yy, h0, color=GREEN, alpha=0.38, rstride=2, cstride=2, linewidth=0.7)
        ax3.plot_surface(xx, yy, hr, color=PURPLE, alpha=0.48, linewidth=0, antialiased=True)
        path0 = common_path[donor, axis]
        pathr = response_path[donor, axis]
        ax3.plot(path0[:, 0], path0[:, 1], path_height(x, y, h0, path0), color=GREEN, lw=3)
        ax3.plot(pathr[:, 0], pathr[:, 1], path_height(x, y, hr, pathr), color=PURPLE, lw=3)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, y[-1])
        ax3.set_zlim(0, 1.15)
        ax3.view_init(elev=27, azim=-58)
        ax3.set_xlabel("severity", labelpad=-1)
        if axis == 0:
            ax3.set_ylabel("off-axis", labelpad=-1)
        ax3.set_zlabel("LAND", labelpad=-2)
        ax3.tick_params(labelsize=7, pad=0)
        ax3.set_title(
            f"test median response Δ={response_rows[axis]['test_shape_change_median']:.2f}", fontsize=9.5, pad=0
        )

    fig.suptitle(
        f"{celltype_name}: exploratory model-estimated personal×pathology response deforms the common fingerprint",
        fontsize=21,
        fontweight="bold",
        color="#173f76",
        y=0.965,
    )
    fig.text(
        0.5,
        0.905,
        "green = same donor context with common branch only | purple = personal baseline still zero, response slope added",
        ha="center",
        fontsize=12.5,
        color="#4d5c73",
    )
    fig.legend(
        handles=(
            Line2D([0], [0], color=GREEN, lw=3, label="common-only"),
            Line2D([0], [0], color=PURPLE, lw=3, label="response-on"),
            Line2D([0], [0], color=ORANGE, lw=2, ls="--", label="coordinate-straight path"),
        ),
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=11,
        bbox_to_anchor=(0.5, 0.02),
    )
    fig.savefig(out.with_suffix(".png"), dpi=190, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--figure13-summary", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--celltype", default="L4 IT")
    parser.add_argument("--grid-x", type=int, default=37)
    parser.add_argument("--grid-y", type=int, default=21)
    parser.add_argument("--y-max", type=float, default=0.90)
    parser.add_argument("--batch-size", type=int, default=192)
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(20260820)
    torch.manual_seed(20260820)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    expected_sha = str(centroid["checkpoint_sha256"].item())
    actual_sha = sha256(args.checkpoint)
    if expected_sha != actual_sha:
        raise RuntimeError(f"checkpoint SHA mismatch: {expected_sha} versus {actual_sha}")
    celltype_names = [str(value) for value in centroid["celltype_vocab"]]
    exact_celltypes = [str(value) for value in exact["celltype_names"]]
    if celltype_names != exact_celltypes:
        raise RuntimeError("centroid and exact-array cell-type vocabularies differ")
    if args.celltype not in celltype_names:
        raise ValueError(f"unknown cell type {args.celltype!r}")
    celltype = celltype_names.index(args.celltype)

    raw = json.load(args.config.open())
    current.set_settings(raw)
    current.install_hooks()
    cfg = base.load_config(args.config)
    ds, obs, _, _ = build_pooled_dataset(cfg)
    device = torch.device(args.device)
    system = build_posthoc_system(cfg, ds, obs)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = payload.get("system_state_dict", payload) if isinstance(payload, dict) else payload
    missing, unexpected = system.load_state_dict(state, strict=False)
    compatibility = validate_posthoc_checkpoint_compatibility(missing, unexpected)
    print(f"[donor-rgf] checkpoint compatibility={compatibility}", flush=True)
    system.to(device).eval()
    for parameter in system.parameters():
        parameter.requires_grad_(False)

    split_all = np.asarray(centroid["split"], dtype=str)
    donors = np.flatnonzero(np.isin(split_all, ("val", "test")) & (centroid["count"][:, celltype] >= 20))
    split = split_all[donors]
    val_index = np.flatnonzero(split == "val")
    test_index = np.flatnonzero(split == "test")
    if len(val_index) != 16 or len(test_index) != 9:
        raise RuntimeError(f"expected 16 validation and 9 test donors, got {len(val_index)} and {len(test_index)}")
    z_full = np.asarray(centroid["z_centroid"], dtype=np.float32)[donors, celltype]
    z_half = np.asarray(centroid["z_half_centroid"], dtype=np.float32)[:, donors, celltype]
    contexts = np.concatenate((z_full[None], z_half), axis=0)
    base_score = np.stack(
        [donor_base_scores(system, context, celltype, device) for context in contexts], axis=0
    )

    common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
    region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
    lift = system.module_tokenizer.activity_weight.detach().float().to(device)
    score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift
    response_module = torch.as_tensor(exact["response_unit"][donors, celltype], dtype=torch.float32, device=device)
    response_score = torch.einsum("dkm,mg->dkg", response_module, lift)
    thresholds = system.decoder._compute_thresholds().detach().float().to(device)

    figure13 = json.load(args.figure13_summary.open())
    row_lookup = {(row["celltype"], row["pathology"]): row for row in figure13["rows"]}
    directions = np.stack(
        [np.asarray(row_lookup[(args.celltype, name)]["slice_direction"], dtype=np.float32) for name in PATHOLOGY]
    )
    shared_sigma = float(figure13["shared_fisher_sigma"])
    x = np.linspace(0.0, 1.0, args.grid_x, dtype=np.float32)
    y = np.linspace(0.0, args.y_max, args.grid_y, dtype=np.float32)
    xx, yy = np.meshgrid(x, y)
    chart_points = np.column_stack((xx.ravel(), yy.ravel())).astype(np.float32)

    shape = (3, len(donors), 5, args.grid_y, args.grid_x)
    normal = np.full(shape, np.nan, dtype=np.float32)
    disease = np.full(shape, np.nan, dtype=np.float32)
    paths = np.full((3, len(donors), 5, 60, 2), np.nan, dtype=np.float32)
    shortening = np.full((3, len(donors), 5), np.nan, dtype=np.float32)
    off_axis = np.full_like(shortening, np.nan)
    response_normal = np.full(shape[1:], np.nan, dtype=np.float32)
    response_disease = np.full(shape[1:], np.nan, dtype=np.float32)
    response_paths = np.full((len(donors), 5, 60, 2), np.nan, dtype=np.float32)
    response_shortening = np.full((len(donors), 5), np.nan, dtype=np.float32)

    for axis in range(5):
        e_axis = np.zeros(5, dtype=np.float32)
        e_axis[axis] = 1.0
        chart_np = np.column_stack((e_axis, directions[axis])).astype(np.float32)
        pathology_points_np = chart_points @ chart_np.T
        if np.min(pathology_points_np) < -1e-6 or np.max(pathology_points_np) > 1.0 + 1e-6:
            raise RuntimeError(f"chart for {PATHOLOGY[axis]} left [0,1]^5")
        chart = torch.as_tensor(chart_np, dtype=torch.float32, device=device)
        pathology_points = torch.as_tensor(pathology_points_np, dtype=torch.float32, device=device)
        for context_id in range(3):
            for donor_index in range(len(donors)):
                metric = evaluate_metric_grid(
                    system.decoder,
                    thresholds,
                    base_score[context_id, donor_index],
                    score_dictionary,
                    chart,
                    pathology_points,
                    None,
                    args.batch_size,
                )
                result = analyse_surface(metric, x, y, shared_sigma)
                normal[context_id, donor_index, axis] = result["normal"]
                disease[context_id, donor_index, axis] = result["disease"]
                paths[context_id, donor_index, axis] = resample_path(result["path"])
                shortening[context_id, donor_index, axis] = result["shortening"]
                off_axis[context_id, donor_index, axis] = result["off_axis"]
        for donor_index in range(len(donors)):
            metric = evaluate_metric_grid(
                system.decoder,
                thresholds,
                base_score[0, donor_index],
                score_dictionary,
                chart,
                pathology_points,
                response_score[donor_index],
                args.batch_size,
            )
            result = analyse_surface(metric, x, y, shared_sigma)
            response_normal[donor_index, axis] = result["normal"]
            response_disease[donor_index, axis] = result["disease"]
            response_paths[donor_index, axis] = resample_path(result["path"])
            response_shortening[donor_index, axis] = result["shortening"]
        print(f"[donor-rgf] {args.celltype} {PATHOLOGY[axis]} complete", flush=True)

    vector = np.empty((3, len(donors), 5, 2 * args.grid_y * args.grid_x), dtype=np.float32)
    response_vector = np.empty((len(donors), 5, 2 * args.grid_y * args.grid_x), dtype=np.float32)
    for context_id in range(3):
        for donor_index in range(len(donors)):
            for axis in range(5):
                vector[context_id, donor_index, axis] = fingerprint_vector(
                    normal[context_id, donor_index, axis], disease[context_id, donor_index, axis]
                )
    for donor_index in range(len(donors)):
        for axis in range(5):
            response_vector[donor_index, axis] = fingerprint_vector(
                response_normal[donor_index, axis], response_disease[donor_index, axis]
            )

    val_template_vector = vector[0, val_index].mean(axis=0)
    val_template_vector /= np.linalg.norm(val_template_vector, axis=1, keepdims=True).clip(min=1e-30)
    val_template_normal = np.stack([template_density(normal[0, val_index, axis]) for axis in range(5)])
    val_template_disease = np.stack([template_density(disease[0, val_index, axis]) for axis in range(5)])
    val_template_path = np.nanmean(paths[0, val_index], axis=0)
    rows: list[dict[str, Any]] = []
    response_rows: list[dict[str, Any]] = []
    representative = np.zeros(5, dtype=np.int16)
    for axis, name in enumerate(PATHOLOGY):
        test_shape = np.asarray([correlation(vector[0, donor, axis], val_template_vector[axis]) for donor in test_index])
        val_shape = np.asarray([correlation(vector[0, donor, axis], val_template_vector[axis]) for donor in val_index])
        split_half = np.asarray([correlation(vector[1, donor, axis], vector[2, donor, axis]) for donor in range(len(donors))])
        test_different = np.asarray(
            [np.nanmedian([correlation(vector[0, donor, axis], val_template_vector[other]) for other in range(5) if other != axis]) for donor in test_index]
        )
        path_rms = np.sqrt(np.mean(np.square(paths[0, test_index, axis] - val_template_path[axis]), axis=(1, 2)))
        response_change = np.asarray(
            [1.0 - correlation(response_vector[donor, axis], vector[0, donor, axis]) for donor in range(len(donors))]
        )
        test_response_change = response_change[test_index]
        representative[axis] = int(test_index[np.nanargmax(test_response_change)])
        test_ci = bootstrap_ci(test_shape, 20260820 + axis)
        half_ci = bootstrap_ci(split_half, 20260920 + axis)
        short_ci = bootstrap_ci(shortening[0, test_index, axis], 20261020 + axis)
        rows.append(
            {
                "pathology": name,
                "validation_pairwise_shape_median": pairwise_median(vector[0, val_index, axis]),
                "test_template_shape_median": test_ci[0],
                "test_template_shape_ci95": test_ci[1:],
                "split_half_median": half_ci[0],
                "split_half_ci95": half_ci[1:],
                "test_same_minus_different_median": float(np.nanmedian(test_shape - test_different)),
                "test_path_rms_median": float(np.nanmedian(path_rms)),
                "validation_shortening_median": float(np.nanmedian(shortening[0, val_index, axis])),
                "test_shortening_median": short_ci[0],
                "test_shortening_ci95": short_ci[1:],
                "validation_off_axis_median": float(np.nanmedian(off_axis[0, val_index, axis])),
                "test_off_axis_median": float(np.nanmedian(off_axis[0, test_index, axis])),
            }
        )
        response_rows.append(
            {
                "pathology": name,
                "validation_shape_change_median": float(np.nanmedian(response_change[val_index])),
                "test_shape_change_median": float(np.nanmedian(test_response_change)),
                "test_path_shift_rms_median": float(
                    np.nanmedian(
                        np.sqrt(np.mean(np.square(response_paths[test_index, axis] - paths[0, test_index, axis]), axis=(1, 2)))
                    )
                ),
                "test_shortening_change_median": float(
                    np.nanmedian(response_shortening[test_index, axis] - shortening[0, test_index, axis])
                ),
                "representative_donor": str(centroid["donor_vocab"][donors[representative[axis]]]),
                "representative_shape_change": float(response_change[representative[axis]]),
            }
        )

    exact_out = args.out_dir / "donor_conditioned_riemann_exact_arrays.npz"
    np.savez_compressed(
        exact_out,
        schema_version=np.asarray("prism.donor_conditioned_riemann.v1", dtype=object),
        checkpoint_sha256=np.asarray(actual_sha, dtype=object),
        checkpoint_epoch=np.asarray(int(payload.get("epoch", -1))),
        celltype=np.asarray(args.celltype, dtype=object),
        pathology_names=np.asarray(PATHOLOGY, dtype=object),
        donor_ids=donors,
        donor_names=np.asarray(centroid["donor_vocab"], dtype=object)[donors],
        split=split,
        count=np.asarray(centroid["count"])[donors, celltype],
        reference_fraction=np.asarray(centroid["reference_fraction"])[donors],
        x=x,
        y=y,
        slice_direction=directions,
        shared_fisher_sigma=np.asarray(shared_sigma),
        normal_density=normal,
        disease_density=disease,
        path=paths,
        shortening_pct=shortening,
        max_off_axis=off_axis,
        response_normal_density=response_normal,
        response_disease_density=response_disease,
        response_path=response_paths,
        response_shortening_pct=response_shortening,
        validation_template_normal=val_template_normal,
        validation_template_disease=val_template_disease,
        validation_template_path=val_template_path,
        representative_donor_index=representative,
    )

    render_common(
        args.out_dir / "figure_06a_l4it_common_fingerprint_validation_test",
        args.celltype,
        x,
        y,
        val_template_normal,
        val_template_disease,
        val_template_path,
        paths[0, test_index],
        rows,
    )
    render_response(
        args.out_dir / "figure_06b_l4it_personal_response_exploratory",
        args.celltype,
        x,
        y,
        normal[0],
        disease[0],
        response_normal,
        response_disease,
        paths[0],
        response_paths,
        np.asarray(centroid["donor_vocab"], dtype=object)[donors],
        representative,
        response_rows,
    )

    summary = {
        "schema_version": "prism.donor_conditioned_riemann.summary.v1",
        "provenance": {
            "checkpoint": str(args.checkpoint.resolve()),
            "checkpoint_sha256": actual_sha,
            "checkpoint_epoch": int(payload.get("epoch", -1)),
            "config": str(args.config.resolve()),
            "centroids": str(args.centroids.resolve()),
            "exact_npz": str(args.exact_npz.resolve()),
            "compatibility": compatibility,
            "locked_test_used_for_refitting": False,
        },
        "estimand": {
            "primary": "common explicit pathology geometry conditional on each donor clean-state centroid",
            "personal_baseline": 0,
            "personal_by_pathology_response": 0,
            "sex_score": 0,
            "technology_score": 0,
            "age": 0,
            "region": "equal-weight mean of DLPFC and MTG explicit gates",
            "donor_context": "mean frozen zclean for donor-by-celltype; not an explicit personal term",
            "metric": "ordinal-Fisher pullback in a fixed validation/common-template 2-D chart",
            "display_warning": "3-D height is a scalar LAND lift, not an isometric embedding or temporal trajectory",
        },
        "data": {
            "celltype": args.celltype,
            "validation_donors": int(len(val_index)),
            "locked_test_donors": int(len(test_index)),
            "cells_per_donor_celltype": np.asarray(centroid["count"])[donors, celltype].tolist(),
            "strict_reference_validation_donors": int(np.sum(np.asarray(centroid["reference_fraction"])[donors[val_index]] > 0.5)),
            "strict_reference_test_donors": int(np.sum(np.asarray(centroid["reference_fraction"])[donors[test_index]] > 0.5)),
        },
        "common_only": rows,
        "personal_response_exploratory": response_rows,
        "claim_rule": {
            "noise_ceiling": "split-half fingerprint correlation",
            "cross_donor": "locked-test donor correlation with frozen validation template",
            "axis_specificity": "same-axis correlation minus median different-axis correlation",
            "normal_replication_limitation": "strict-reference normal has no locked-test donors, so normal biology is not independently test-replicated",
            "response_limitation": "response branch is model-estimated and remains exploratory until direct branch ablation supports it",
        },
        "files": {
            "exact_arrays": exact_out.name,
            "common_figure_png": "figure_06a_l4it_common_fingerprint_validation_test.png",
            "common_figure_pdf": "figure_06a_l4it_common_fingerprint_validation_test.pdf",
            "response_figure_png": "figure_06b_l4it_personal_response_exploratory.png",
            "response_figure_pdf": "figure_06b_l4it_personal_response_exploratory.pdf",
        },
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (args.out_dir / "README_KO.md").write_text(
        "# Donor-conditioned PRISM Riemann fingerprint\n\n"
        "동결된 personal-rank2 epoch 20만 사용한 사후분석이다. 주 분석은 explicit personal과 response, 성별, 기술 점수를 0으로 두고 donor별 zclean 평균만 decoder 문맥으로 남긴다. validation 16명에서 평균 지문을 고정한 뒤 locked test 9명을 재학습 없이 투영했다. 정상점의 엄격한 reference donor는 validation 1명, test 0명이므로 정상 생물학의 독립 test 재현으로 해석하면 안 된다. response 그림은 branch ablation 전까지 탐색적 모델 추정치이다.\n",
        encoding="utf-8",
    )
    print(f"[donor-rgf] wrote {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
