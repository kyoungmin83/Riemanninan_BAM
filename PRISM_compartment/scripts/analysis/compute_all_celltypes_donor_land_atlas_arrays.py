#!/usr/bin/env python3
"""Compute validation/test donor-conditioned LAND surfaces for all cell types.

The display chart is fixed from the original Figure 13 common template.  The
green path saved below is the actual saved full-5D Figure-13 path projected to
that chart; it is not a geodesic re-optimized in the 2-D display slice.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_donor_conditioned_riemann_fingerprint import (
    PATHOLOGY,
    analyse_surface,
    donor_base_scores,
    evaluate_metric_grid,
    json_safe,
    sha256,
)
from disease_program_eval import build_pooled_dataset
from evaluate_donor_conditioned_full5_paths import path_length
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


def resample_path(path: np.ndarray, size: int = 121) -> np.ndarray:
    delta = np.diff(path, axis=0)
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(delta, axis=1))))
    if arc[-1] <= 1e-12:
        return np.repeat(path[:1], size, axis=0)
    target = np.linspace(0.0, arc[-1], size)
    return np.column_stack([np.interp(target, arc, path[:, dim]) for dim in range(path.shape[1])])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--global5d-npz", required=True, type=Path)
    parser.add_argument("--figure13-summary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--minimum-cells", type=int, default=20)
    parser.add_argument("--grid-x", type=int, default=45)
    parser.add_argument("--grid-y", type=int, default=25)
    parser.add_argument("--y-max", type=float, default=0.90)
    parser.add_argument("--batch-size", type=int, default=192)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    global5 = np.load(args.global5d_npz, allow_pickle=True)
    figure13 = json.load(args.figure13_summary.open())
    actual_sha = sha256(args.checkpoint)
    if actual_sha != str(centroid["checkpoint_sha256"].item()):
        raise RuntimeError("checkpoint SHA mismatch")
    celltypes = [str(value) for value in centroid["celltype_vocab"]]
    if celltypes != [str(value) for value in exact["celltype_names"]]:
        raise RuntimeError("exact-array cell-type vocabulary mismatch")
    if celltypes != [str(value) for value in global5["celltype_names"]]:
        raise RuntimeError("global-5D cell-type vocabulary mismatch")
    row_lookup = {(row["celltype"], row["pathology"]): row for row in figure13["rows"]}

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
    system.to(device).eval()
    for parameter in system.parameters():
        parameter.requires_grad_(False)

    lift = system.module_tokenizer.activity_weight.detach().float().to(device)
    thresholds = system.decoder._compute_thresholds().detach().float().to(device)
    split_all = np.asarray(centroid["split"], dtype=str)
    count_all = np.asarray(centroid["count"])
    z_all = np.asarray(centroid["z_centroid"], dtype=np.float32)
    path_saved = np.asarray(global5["selected_path"], dtype=np.float32)
    path_count = np.asarray(global5["selected_path_node_count"], dtype=int)
    sigma = float(figure13["shared_fisher_sigma"])
    x = np.linspace(0.0, 1.0, args.grid_x, dtype=np.float32)
    y = np.linspace(0.0, args.y_max, args.grid_y, dtype=np.float32)
    xx, yy = np.meshgrid(x, y)
    display_points = np.column_stack((xx.ravel(), yy.ravel())).astype(np.float32)

    c_count = len(celltypes)
    template_shape = (2, c_count, 5, args.grid_y, args.grid_x)
    normal_template = np.full(template_shape, np.nan, dtype=np.float32)
    disease_template = np.full(template_shape, np.nan, dtype=np.float32)
    projected_path = np.full((c_count, 5, 121, 2), np.nan, dtype=np.float32)
    direct_path = np.zeros((c_count, 5, 121, 2), dtype=np.float32)
    path5d_resampled = np.full((c_count, 5, 121, 5), np.nan, dtype=np.float32)
    shortening_median = np.full((2, c_count, 5), np.nan, dtype=np.float32)
    shortening_q10 = np.full_like(shortening_median, np.nan)
    shortening_q90 = np.full_like(shortening_median, np.nan)
    shorter_fraction = np.full_like(shortening_median, np.nan)
    val_count = np.zeros(c_count, dtype=np.int16)
    test_count = np.zeros(c_count, dtype=np.int16)
    directions = np.empty((c_count, 5, 5), dtype=np.float32)
    projection_capture = np.empty((c_count, 5), dtype=np.float32)
    summary_rows: list[dict[str, object]] = []

    for celltype, name in enumerate(celltypes):
        donors = np.flatnonzero(
            np.isin(split_all, ("val", "test")) & (count_all[:, celltype] >= args.minimum_cells)
        )
        split = split_all[donors]
        val = np.flatnonzero(split == "val")
        test = np.flatnonzero(split == "test")
        val_count[celltype], test_count[celltype] = len(val), len(test)
        if len(val) < 5 or len(test) < 3:
            print(f"[land-atlas] skip {name}: val={len(val)} test={len(test)}", flush=True)
            continue
        base_score = donor_base_scores(system, z_all[donors, celltype], celltype, device)
        common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
        region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
        score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift

        for axis, pathology_name in enumerate(PATHOLOGY):
            source_row = row_lookup[(name, pathology_name)]
            direction = np.asarray(source_row["slice_direction"], dtype=np.float32)
            directions[celltype, axis] = direction
            projection_capture[celltype, axis] = float(source_row["dominant_off_axis_projection_captured"])
            named = np.zeros(5, dtype=np.float32)
            named[axis] = 1.0
            chart_np = np.column_stack((named, direction)).astype(np.float32)
            pathology_points_np = display_points @ chart_np.T
            if np.min(pathology_points_np) < -1e-6 or np.max(pathology_points_np) > 1.0 + 1e-6:
                raise RuntimeError(f"display chart for {name}|{pathology_name} left [0,1]^5")
            chart = torch.as_tensor(chart_np, dtype=torch.float32, device=device)
            pathology_points = torch.as_tensor(pathology_points_np, dtype=torch.float32, device=device)
            normal_donor = np.empty((len(donors), args.grid_y, args.grid_x), dtype=np.float32)
            disease_donor = np.empty_like(normal_donor)
            for donor in range(len(donors)):
                metric = evaluate_metric_grid(
                    system.decoder,
                    thresholds,
                    base_score[donor],
                    score_dictionary,
                    chart,
                    pathology_points,
                    None,
                    args.batch_size,
                )
                result = analyse_surface(metric, x, y, sigma)
                normal_donor[donor] = result["normal"]
                disease_donor[donor] = result["disease"]
            for split_index, member in enumerate((val, test)):
                normal_template[split_index, celltype, axis] = np.mean(normal_donor[member], axis=0)
                disease_template[split_index, celltype, axis] = np.mean(disease_donor[member], axis=0)

            row = celltype * 5 + axis
            fixed5_raw = path_saved[row, : path_count[row]]
            # Preserve the exact saved polyline for every length calculation.
            # Independent-coordinate resampling can cut across sharp optimized
            # vertices and would otherwise make the displayed path look
            # artificially shorter.  Resampling is display-only.
            fixed5_display = resample_path(fixed5_raw)
            path5d_resampled[celltype, axis] = fixed5_display
            projected = np.column_stack(
                (fixed5_display[:, axis], fixed5_display @ direction)
            )
            projected_path[celltype, axis] = projected
            direct_path[celltype, axis, :, 0] = np.linspace(0.0, 1.0, 121)
            if np.max(projected[:, 1]) > args.y_max + 1e-6 or np.min(projected[:, 1]) < -1e-6:
                raise RuntimeError(
                    f"projected Figure-13 path for {name}|{pathology_name} is outside display y-range: "
                    f"[{projected[:, 1].min():.3f},{projected[:, 1].max():.3f}]"
                )
            direct5 = np.linspace(0.0, 1.0, 73, dtype=np.float32)[:, None] * named[None]
            donor_shortening = np.empty(len(donors), dtype=np.float64)
            for donor in range(len(donors)):
                direct_length = path_length(
                    system.decoder, thresholds, base_score[donor], score_dictionary, direct5, None
                )
                fixed_length = path_length(
                    system.decoder,
                    thresholds,
                    base_score[donor],
                    score_dictionary,
                    fixed5_raw,
                    None,
                )
                donor_shortening[donor] = 100.0 * (direct_length - fixed_length) / max(direct_length, 1e-30)
            for split_index, member in enumerate((val, test)):
                values = donor_shortening[member]
                shortening_median[split_index, celltype, axis] = float(np.median(values))
                shortening_q10[split_index, celltype, axis] = float(np.percentile(values, 10))
                shortening_q90[split_index, celltype, axis] = float(np.percentile(values, 90))
                shorter_fraction[split_index, celltype, axis] = float(np.mean(values > 0))
            summary_rows.append(
                {
                    "celltype": name,
                    "pathology": pathology_name,
                    "validation_donors": len(val),
                    "locked_test_donors": len(test),
                    "validation_shortening_median_pct": float(shortening_median[0, celltype, axis]),
                    "locked_test_shortening_median_pct": float(shortening_median[1, celltype, axis]),
                    "validation_shorter_fraction": float(shorter_fraction[0, celltype, axis]),
                    "locked_test_shorter_fraction": float(shorter_fraction[1, celltype, axis]),
                    "projected_max_off_axis": float(np.max(projected[:, 1])),
                    "projection_capture": float(projection_capture[celltype, axis]),
                }
            )
        print(f"[land-atlas] {celltype + 1}/{c_count} {name} complete", flush=True)

    np.savez_compressed(
        args.out,
        schema_version=np.asarray("prism.all_celltypes_donor_land_atlas.v1", dtype=object),
        checkpoint_sha256=np.asarray(actual_sha, dtype=object),
        checkpoint_epoch=np.asarray(int(payload.get("epoch", -1))),
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(PATHOLOGY, dtype=object),
        split_names=np.asarray(("validation", "locked_test"), dtype=object),
        validation_donor_count=val_count,
        locked_test_donor_count=test_count,
        x=x,
        y=y,
        shared_fisher_sigma=np.asarray(sigma),
        normal_density_template=normal_template,
        disease_density_template=disease_template,
        slice_direction=directions,
        projection_capture=projection_capture,
        projected_saved_full5_path=projected_path,
        direct_coordinate_path=direct_path,
        saved_full5_path=path5d_resampled,
        shortening_median_pct=shortening_median,
        shortening_q10_pct=shortening_q10,
        shortening_q90_pct=shortening_q90,
        shorter_fraction=shorter_fraction,
    )
    summary = {
        "schema_version": "prism.all_celltypes_donor_land_atlas.summary.v1",
        "checkpoint_sha256": actual_sha,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "compatibility": compatibility,
        "output_npz": str(args.out.resolve()),
        "celltypes": celltypes,
        "pathologies": list(PATHOLOGY),
        "validation_template": "mean of individually normalized donor-conditioned LAND densities",
        "locked_test_template": "same calculation with locked-test donors; no model or path refit",
        "green_path": "actual saved full-5D Figure-13 common path projected into the fixed Figure-13 display chart",
        "orange_path": "one-axis coordinate-Euclidean path",
        "three_dimensional_height": "scalar LAND mixture lift; not an isometric embedding or temporal trajectory",
        "personal_baseline": 0,
        "personal_response": 0,
        "sex_score": 0,
        "technology_score": 0,
        "rows": summary_rows,
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[land-atlas] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
