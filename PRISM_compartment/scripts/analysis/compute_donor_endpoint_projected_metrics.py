#!/usr/bin/env python3
"""Compute donor-specific 2-D endpoint Fisher metrics for contour overlays."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_donor_conditioned_riemann_fingerprint import donor_base_scores, json_safe, sha256
from disease_program_eval import build_pooled_dataset
from evaluate_all_celltypes_donor_riemann_summary import endpoint_metric
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--slice-npz", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--minimum-cells", type=int, default=20)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    slice_data = np.load(args.slice_npz, allow_pickle=True)
    actual_sha = sha256(args.checkpoint)
    if actual_sha != str(centroid["checkpoint_sha256"].item()):
        raise RuntimeError("checkpoint SHA mismatch")
    celltypes = [str(value) for value in centroid["celltype_vocab"]]
    pathologies = [str(value) for value in centroid["pathology_names"]]
    if celltypes != [str(value) for value in exact["celltype_names"]]:
        raise RuntimeError("exact-array cell-type vocabulary mismatch")
    if celltypes != [str(value) for value in slice_data["celltype_names"]]:
        raise RuntimeError("slice-array cell-type vocabulary mismatch")
    direction_lookup = {
        (int(celltype), int(axis)): row
        for row, (celltype, axis) in enumerate(
            zip(slice_data["pair_celltype"], slice_data["pair_axis"])
        )
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
    donor_count = len(split)
    projected_normal = np.full((donor_count, len(celltypes), 5, 2, 2), np.nan, dtype=np.float32)
    projected_disease = np.full_like(projected_normal, np.nan)
    eligible = np.isin(split, ("val", "test"))[:, None] & (count >= args.minimum_cells)
    directions = np.empty((len(celltypes), 5, 5), dtype=np.float32)
    sigma = np.empty((len(celltypes), 5), dtype=np.float32)

    for celltype, name in enumerate(celltypes):
        donors = np.flatnonzero(eligible[:, celltype])
        if len(donors) == 0:
            continue
        base_score = donor_base_scores(system, z[donors, celltype], celltype, device)
        common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
        region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
        score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift
        normal5 = endpoint_metric(
            system.decoder,
            thresholds,
            base_score,
            score_dictionary,
            np.zeros(5, dtype=np.float32),
            None,
        )
        for axis, _ in enumerate(pathologies):
            direction_row = direction_lookup[(celltype, axis)]
            direction = np.asarray(slice_data["slice_direction"][direction_row], dtype=np.float64)
            directions[celltype, axis] = direction
            sigma[celltype, axis] = max(
                0.085 * float(slice_data["full5_geodesic_length"][direction_row]), 1e-6
            )
            named = np.zeros(5, dtype=np.float64)
            named[axis] = 1.0
            chart = np.column_stack((named, direction))
            disease_point = named.astype(np.float32)
            disease5 = endpoint_metric(
                system.decoder,
                thresholds,
                base_score,
                score_dictionary,
                disease_point,
                None,
            )
            projected_normal[donors, celltype, axis] = np.einsum(
                "ai,dab,bj->dij", chart, normal5, chart, optimize=True
            )
            projected_disease[donors, celltype, axis] = np.einsum(
                "ai,dab,bj->dij", chart, disease5, chart, optimize=True
            )
        print(f"[donor-endpoint] {celltype + 1}/{len(celltypes)} {name} donors={len(donors)}", flush=True)

    np.savez_compressed(
        args.out,
        schema_version=np.asarray("prism.donor_endpoint_projected_metrics.v1", dtype=object),
        checkpoint_sha256=np.asarray(actual_sha, dtype=object),
        checkpoint_epoch=np.asarray(int(payload.get("epoch", -1))),
        donor_names=np.asarray(centroid["donor_vocab"], dtype=object),
        donor_split=split.astype(object),
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(pathologies, dtype=object),
        eligible=eligible,
        cell_count=count,
        slice_direction=directions,
        sigma_fisher=sigma,
        normal_metric_g2=projected_normal,
        disease_metric_g2=projected_disease,
    )
    summary = {
        "schema_version": "prism.donor_endpoint_projected_metrics.summary.v1",
        "checkpoint_sha256": actual_sha,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "compatibility": compatibility,
        "output": str(args.out.resolve()),
        "minimum_cells": args.minimum_cells,
        "validation_donors_by_celltype": {
            name: int(np.sum(eligible[:, index] & (split == "val"))) for index, name in enumerate(celltypes)
        },
        "locked_test_donors_by_celltype": {
            name: int(np.sum(eligible[:, index] & (split == "test"))) for index, name in enumerate(celltypes)
        },
        "personal_baseline": 0,
        "personal_response": 0,
        "sex_score": 0,
        "technology_score": 0,
        "region": "equal-weight mean of DLPFC and MTG explicit gates",
        "test_refit": False,
        "metric": "ordinal-Fisher endpoint metric projected into the fixed Figure-13 named/off-axis chart",
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[donor-endpoint] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
