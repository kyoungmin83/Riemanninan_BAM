#!/usr/bin/env python3
"""Evaluate frozen Figure-13 full-5D paths in donor-conditioned score contexts."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_donor_conditioned_riemann_fingerprint import (
    PATHOLOGY,
    bootstrap_ci,
    donor_base_scores,
    fisher_info,
    hinge,
    hinge_derivative,
    json_safe,
    sha256,
)
from disease_program_eval import build_pooled_dataset
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


@torch.no_grad()
def path_length(
    decoder,
    thresholds: torch.Tensor,
    base_score: np.ndarray,
    score_dictionary: torch.Tensor,
    path: np.ndarray,
    response_score: torch.Tensor | None,
) -> float:
    point = torch.as_tensor(path, dtype=torch.float32, device=score_dictionary.device)
    midpoint = 0.5 * (point[:-1] + point[1:])
    delta = point[1:] - point[:-1]
    feature = hinge(midpoint)
    score = torch.as_tensor(base_score, dtype=torch.float32, device=point.device).unsqueeze(0)
    score = score + torch.einsum("skf,kfg->sg", feature, score_dictionary)
    jacobian = torch.einsum("skf,kfg->skg", hinge_derivative(midpoint), score_dictionary)
    if response_score is not None:
        score = score + midpoint @ response_score
        jacobian = jacobian + response_score.unsqueeze(0)
    info = fisher_info(decoder, score, thresholds)
    metric = torch.einsum("skg,sg,slg->skl", jacobian, info, jacobian) / float(decoder.n_genes)
    square = torch.einsum("si,sij,sj->s", delta, metric, delta).clamp_min(1e-30)
    return float(torch.sqrt(square).sum().double().cpu())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--global5d-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--celltype", default="L4 IT")
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    global5 = np.load(args.global5d_npz, allow_pickle=True)
    actual_sha = sha256(args.checkpoint)
    if actual_sha != str(centroid["checkpoint_sha256"].item()):
        raise RuntimeError("checkpoint SHA mismatch")
    celltypes = [str(value) for value in centroid["celltype_vocab"]]
    global_celltypes = [str(value) for value in global5["celltype_names"]]
    if args.celltype not in celltypes or args.celltype not in global_celltypes:
        raise ValueError(f"unknown cell type {args.celltype!r}")
    celltype = celltypes.index(args.celltype)
    global_celltype = global_celltypes.index(args.celltype)

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

    split_all = np.asarray(centroid["split"], dtype=str)
    count_all = np.asarray(centroid["count"])
    donors = np.flatnonzero(np.isin(split_all, ("val", "test")) & (count_all[:, celltype] >= 20))
    split = split_all[donors]
    val = np.flatnonzero(split == "val")
    test = np.flatnonzero(split == "test")
    z = np.asarray(centroid["z_centroid"], dtype=np.float32)[donors, celltype]
    base_score = donor_base_scores(system, z, celltype, device)

    lift = system.module_tokenizer.activity_weight.detach().float().to(device)
    common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
    region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
    score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift
    response_module = torch.as_tensor(exact["response_unit"][donors, celltype], dtype=torch.float32, device=device)
    response_score = torch.einsum("dkm,mg->dkg", response_module, lift)
    thresholds = system.decoder._compute_thresholds().detach().float().to(device)

    direct_length = np.empty((len(donors), 5), dtype=np.float64)
    candidate_length = np.empty_like(direct_length)
    response_direct_length = np.empty_like(direct_length)
    response_candidate_length = np.empty_like(direct_length)
    path_count = np.asarray(global5["selected_path_node_count"], dtype=int)
    path_saved = np.asarray(global5["selected_path"], dtype=np.float32)
    fixed_paths: list[np.ndarray] = []
    direct_paths: list[np.ndarray] = []
    for axis in range(5):
        row = global_celltype * 5 + axis
        fixed = path_saved[row, : path_count[row]]
        target = np.zeros(5, dtype=np.float32)
        target[axis] = 1.0
        direct = np.linspace(0.0, 1.0, 73, dtype=np.float32)[:, None] * target[None]
        fixed_paths.append(fixed)
        direct_paths.append(direct)
        for donor in range(len(donors)):
            direct_length[donor, axis] = path_length(
                system.decoder, thresholds, base_score[donor], score_dictionary, direct, None
            )
            candidate_length[donor, axis] = path_length(
                system.decoder, thresholds, base_score[donor], score_dictionary, fixed, None
            )
            response_direct_length[donor, axis] = path_length(
                system.decoder, thresholds, base_score[donor], score_dictionary, direct, response_score[donor]
            )
            response_candidate_length[donor, axis] = path_length(
                system.decoder, thresholds, base_score[donor], score_dictionary, fixed, response_score[donor]
            )
        print(f"[full5-robustness] {PATHOLOGY[axis]} complete", flush=True)

    shortening = 100.0 * (direct_length - candidate_length) / np.maximum(direct_length, 1e-30)
    response_shortening = 100.0 * (response_direct_length - response_candidate_length) / np.maximum(
        response_direct_length, 1e-30
    )
    rows = []
    for axis, name in enumerate(PATHOLOGY):
        val_ci = bootstrap_ci(shortening[val, axis], 20260820 + axis)
        test_ci = bootstrap_ci(shortening[test, axis], 20260920 + axis)
        response_ci = bootstrap_ci(response_shortening[test, axis], 20261020 + axis)
        rows.append(
            {
                "pathology": name,
                "validation_fixed_path_shortening_median_pct": val_ci[0],
                "validation_fixed_path_shortening_ci95_pct": val_ci[1:],
                "test_fixed_path_shortening_median_pct": test_ci[0],
                "test_fixed_path_shortening_ci95_pct": test_ci[1:],
                "test_fixed_path_shorter_fraction": float(np.mean(shortening[test, axis] > 0)),
                "test_response_on_shortening_median_pct": response_ci[0],
                "test_response_on_shortening_ci95_pct": response_ci[1:],
                "common_template_max_off_axis": float(np.nanmax(np.abs(fixed_paths[axis][:, np.arange(5) != axis]))),
            }
        )

    np.savez_compressed(
        args.out_dir / "full5_path_robustness_exact_arrays.npz",
        schema_version=np.asarray("prism.donor_conditioned_full5_path_robustness.v1", dtype=object),
        donor_ids=donors,
        donor_names=np.asarray(centroid["donor_vocab"], dtype=object)[donors],
        split=split,
        pathology_names=np.asarray(PATHOLOGY, dtype=object),
        direct_length=direct_length,
        fixed_path_length=candidate_length,
        fixed_path_shortening_pct=shortening,
        response_direct_length=response_direct_length,
        response_fixed_path_length=response_candidate_length,
        response_fixed_path_shortening_pct=response_shortening,
    )

    fig, axes = plt.subplots(1, 5, figsize=(19, 4.8), sharey=True)
    for axis, name in enumerate(PATHOLOGY):
        ax = axes[axis]
        rng = np.random.default_rng(20260820 + axis)
        for donor_set, xpos, color, label in (
            (val, 0.0, "#2f78bd", "validation"),
            (test, 1.0, "#7b3294", "locked test"),
        ):
            jitter = rng.uniform(-0.12, 0.12, len(donor_set))
            ax.scatter(xpos + jitter, shortening[donor_set, axis], s=32, color=color, alpha=0.78, edgecolor="white", lw=0.4)
            median = float(np.median(shortening[donor_set, axis]))
            ax.plot((xpos - 0.22, xpos + 0.22), (median, median), color=color, lw=3)
        ax.axhline(0, color="0.2", lw=1, ls="--")
        ax.set_xticks((0, 1), ("Val", "Test"))
        ax.set_xlim(-0.45, 1.45)
        ax.grid(axis="y", alpha=0.2)
        ax.set_title(
            f"{name}\nTest median={rows[axis]['test_fixed_path_shortening_median_pct']:.2f}%\n"
            f"shorter in {100*rows[axis]['test_fixed_path_shorter_fraction']:.0f}%",
            fontsize=11,
            fontweight="bold",
        )
        if axis == 0:
            ax.set_ylabel("fixed Figure-13 path shortening vs one-axis route (%)")
    fig.suptitle(
        f"{args.celltype}: the frozen full-5D Figure-13 detour was re-scored in each donor context",
        fontsize=18,
        fontweight="bold",
        color="#173f76",
        y=1.03,
    )
    fig.text(
        0.5,
        -0.02,
        "Positive = the pre-specified Figure-13 path remains shorter; no test-donor path was re-optimized.",
        ha="center",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.out_dir / "figure_06c_full5_path_robustness.png", dpi=190, bbox_inches="tight")
    fig.savefig(args.out_dir / "figure_06c_full5_path_robustness.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "schema_version": "prism.donor_conditioned_full5_path_robustness.summary.v1",
        "checkpoint_sha256": actual_sha,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "compatibility": compatibility,
        "celltype": args.celltype,
        "validation_donors": int(len(val)),
        "locked_test_donors": int(len(test)),
        "test_path_refit": False,
        "path_definition": "saved common-only whole-5D Figure-13 path fixed before donor-conditioned evaluation",
        "personal_baseline": 0,
        "personal_response_primary": 0,
        "rows": rows,
        "warning": "The fixed path is a model-space shortest-path candidate, not a temporal or causal disease trajectory.",
    }
    (args.out_dir / "full5_path_robustness_summary.json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[full5-robustness] wrote {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
