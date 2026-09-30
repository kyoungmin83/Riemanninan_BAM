#!/usr/bin/env python3
"""Compute donor-specific local Fisher metrics along saved common 5-D paths.

This is a frozen-checkpoint posthoc calculation.  It does not optimize a path,
fit a donor template, or update the model.  The same saved common path and the
same fixed Figure-13 display chart are used for validation and locked test.
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
    donor_base_scores,
    evaluate_metric_grid,
    json_safe,
    sha256,
)
from disease_program_eval import build_pooled_dataset
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--common-path-npz", required=True, type=Path)
    parser.add_argument("--slice-npz", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--minimum-cells", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=192)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    common = np.load(args.common_path_npz, allow_pickle=True)
    slice_data = np.load(args.slice_npz, allow_pickle=True)
    actual_sha = sha256(args.checkpoint)
    if actual_sha != str(centroid["checkpoint_sha256"].item()):
        raise RuntimeError("checkpoint SHA mismatch")
    celltypes = [str(value) for value in centroid["celltype_vocab"]]
    pathologies = [str(value) for value in common["pathology_names"]]
    if celltypes != [str(value) for value in exact["celltype_names"]]:
        raise RuntimeError("centroid/exact cell-type vocabulary mismatch")
    if celltypes != [str(value) for value in common["celltype_names"]]:
        raise RuntimeError("centroid/common cell-type vocabulary mismatch")
    pair_lookup = {
        (int(celltype), int(axis)): row
        for row, (celltype, axis) in enumerate(zip(common["pair_celltype"], common["pair_axis"]))
    }
    slice_lookup = {
        (int(celltype), int(axis)): row
        for row, (celltype, axis) in enumerate(zip(slice_data["pair_celltype"], slice_data["pair_axis"]))
    }

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
    split = np.asarray(centroid["split"], dtype=str)
    count = np.asarray(centroid["count"])
    z = np.asarray(centroid["z_centroid"], dtype=np.float32)
    eligible = np.isin(split, ("val", "test"))[:, None] & (count >= args.minimum_cells)
    donor_count = len(split)
    celltype_count = len(celltypes)
    pathology_count = len(pathologies)
    node_count = int(common["path_geodesic"].shape[1])

    path5 = np.full((celltype_count, pathology_count, node_count, 5), np.nan, dtype=np.float32)
    path_xy = np.full((celltype_count, pathology_count, node_count, 2), np.nan, dtype=np.float32)
    path_metric = np.full(
        (donor_count, celltype_count, pathology_count, node_count, 2, 2),
        np.nan,
        dtype=np.float32,
    )
    directions = np.full((celltype_count, pathology_count, 5), np.nan, dtype=np.float32)
    sigma = np.full((celltype_count, pathology_count), np.nan, dtype=np.float32)

    for celltype, name in enumerate(celltypes):
        donor_indices = np.flatnonzero(eligible[:, celltype])
        if len(donor_indices) == 0:
            continue
        base_score = donor_base_scores(system, z[donor_indices, celltype], celltype, device)
        common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
        region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
        score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift
        for axis in range(pathology_count):
            row = pair_lookup[(celltype, axis)]
            slice_row = slice_lookup[(celltype, axis)]
            saved_path = np.asarray(common["path_geodesic"][row], dtype=np.float32)
            direction = np.asarray(slice_data["slice_direction"][slice_row], dtype=np.float32)
            named = np.zeros(5, dtype=np.float32)
            named[axis] = 1.0
            chart_np = np.column_stack((named, direction)).astype(np.float32)
            chart = torch.as_tensor(chart_np, dtype=torch.float32, device=device)
            points = torch.as_tensor(saved_path, dtype=torch.float32, device=device)
            path5[celltype, axis] = saved_path
            path_xy[celltype, axis, :, 0] = saved_path[:, axis]
            path_xy[celltype, axis, :, 1] = saved_path @ direction
            directions[celltype, axis] = direction
            sigma[celltype, axis] = max(0.085 * float(slice_data["full5_geodesic_length"][slice_row]), 1e-6)
            for local_index, donor in enumerate(donor_indices):
                path_metric[donor, celltype, axis] = evaluate_metric_grid(
                    system.decoder,
                    thresholds,
                    base_score[local_index],
                    score_dictionary,
                    chart,
                    points,
                    None,
                    args.batch_size,
                )
        print(f"[donor-local-path] {celltype + 1}/{celltype_count} {name} donors={len(donor_indices)}", flush=True)

    np.savez_compressed(
        args.out,
        schema_version=np.asarray("prism.donor_local_fisher_path_metrics.v1", dtype=object),
        checkpoint_sha256=np.asarray(actual_sha, dtype=object),
        checkpoint_epoch=np.asarray(int(payload.get("epoch", -1))),
        donor_names=np.asarray(centroid["donor_vocab"], dtype=object),
        donor_split=split.astype(object),
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(pathologies, dtype=object),
        eligible=eligible,
        cell_count=count,
        path5=path5,
        path_xy=path_xy,
        path_metric_g2=path_metric,
        slice_direction=directions,
        sigma_fisher=sigma,
    )
    summary = {
        "schema_version": "prism.donor_local_fisher_path_metrics.summary.v1",
        "checkpoint_sha256": actual_sha,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "compatibility": compatibility,
        "output": str(args.out.resolve()),
        "minimum_cells": args.minimum_cells,
        "path": "same saved 5-D common Fisher geodesic for validation and locked test",
        "chart": "same fixed Figure-13 named-axis/dominant-off-axis chart",
        "personal_baseline": 0,
        "personal_response": 0,
        "sex_score": 0,
        "technology_score": 0,
        "donor_context": "frozen donor-by-celltype zclean centroid retained in decoder baseline",
        "test_refit": False,
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[donor-local-path] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
