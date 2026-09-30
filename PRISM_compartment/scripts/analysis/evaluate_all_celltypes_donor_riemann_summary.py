#!/usr/bin/env python3
"""All-cell-type donor-conditioned endpoint and full-5D path audit."""

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
    correlation,
    donor_base_scores,
    fisher_info,
    hinge,
    hinge_derivative,
    json_safe,
    sha256,
)
from disease_program_eval import build_pooled_dataset
from evaluate_donor_conditioned_full5_paths import path_length
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


def matrix_log(metric: np.ndarray) -> np.ndarray:
    metric = 0.5 * (metric + metric.T)
    trace = max(float(np.trace(metric)), 1e-30)
    metric = metric / trace
    value, vector = np.linalg.eigh(metric + 1e-7 * np.eye(len(metric)))
    return (vector * np.log(np.clip(value, 1e-30, None))) @ vector.T


def metric_vector(normal: np.ndarray, disease: np.ndarray) -> np.ndarray:
    upper = np.triu_indices(5)
    value = np.concatenate((matrix_log(normal)[upper], matrix_log(disease)[upper]))
    value -= value.mean()
    return value / max(float(np.linalg.norm(value)), 1e-30)


def log_metric_distance(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(matrix_log(left) - matrix_log(right)))


@torch.no_grad()
def endpoint_metric(
    decoder,
    thresholds: torch.Tensor,
    base_score: np.ndarray,
    score_dictionary: torch.Tensor,
    pathology: np.ndarray,
    response_score: torch.Tensor | None,
) -> np.ndarray:
    p = torch.as_tensor(pathology, dtype=torch.float32, device=score_dictionary.device).reshape(1, 5)
    feature = hinge(p)[0]
    derivative = hinge_derivative(p)[0]
    common_score = torch.einsum("kf,kfg->g", feature, score_dictionary)
    common_jacobian = torch.einsum("kf,kfg->kg", derivative, score_dictionary)
    base = torch.as_tensor(base_score, dtype=torch.float32, device=p.device)
    score = base + common_score.unsqueeze(0)
    jacobian = common_jacobian.unsqueeze(0).expand(len(base), -1, -1)
    if response_score is not None:
        score = score + torch.einsum("k,dkg->dg", p[0], response_score)
        jacobian = jacobian + response_score
    info = fisher_info(decoder, score, thresholds)
    metric = torch.einsum("dkg,dg,dlg->dkl", jacobian, info, jacobian) / float(decoder.n_genes)
    return metric.double().cpu().numpy()


def heatmap(ax, values: np.ndarray, rows: list[str], title: str, cmap: str, vmin=None, vmax=None, fmt=".3f") -> None:
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(5), PATHOLOGY, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(rows)), rows, fontsize=8)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    for i in range(len(rows)):
        for j in range(5):
            if np.isfinite(values[i, j]):
                ax.text(j, i, format(values[i, j], fmt), ha="center", va="center", fontsize=6.1,
                        color="white" if abs(values[i, j] - np.nanmean(values)) > 0.65 * np.nanstd(values) else "black")
    plt.colorbar(image, ax=ax, fraction=0.028, pad=0.015)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--centroids", required=True, type=Path)
    parser.add_argument("--exact-npz", required=True, type=Path)
    parser.add_argument("--global5d-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--minimum-cells", type=int, default=20)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    centroid = np.load(args.centroids, allow_pickle=True)
    exact = np.load(args.exact_npz, allow_pickle=True)
    global5 = np.load(args.global5d_npz, allow_pickle=True)
    actual_sha = sha256(args.checkpoint)
    if actual_sha != str(centroid["checkpoint_sha256"].item()):
        raise RuntimeError("checkpoint SHA mismatch")
    celltypes = [str(value) for value in centroid["celltype_vocab"]]
    if celltypes != [str(value) for value in global5["celltype_names"]]:
        raise RuntimeError("global-5D and centroid cell-type vocabularies differ")

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
    z_half_all = np.asarray(centroid["z_half_centroid"], dtype=np.float32)
    path_saved = np.asarray(global5["selected_path"], dtype=np.float32)
    path_count = np.asarray(global5["selected_path_node_count"], dtype=int)

    dimensions = (len(celltypes), 5)
    test_similarity = np.full(dimensions, np.nan)
    split_half_similarity = np.full(dimensions, np.nan)
    axis_specificity = np.full(dimensions, np.nan)
    full5_shortening = np.full(dimensions, np.nan)
    full5_shorter_fraction = np.full(dimensions, np.nan)
    response_metric_deformation = np.full(dimensions, np.nan)
    response_shortening = np.full(dimensions, np.nan)
    val_n = np.zeros(len(celltypes), dtype=int)
    test_n = np.zeros(len(celltypes), dtype=int)
    summary_rows: list[dict[str, object]] = []

    for celltype, name in enumerate(celltypes):
        donors = np.flatnonzero(
            np.isin(split_all, ("val", "test")) & (count_all[:, celltype] >= args.minimum_cells)
        )
        split = split_all[donors]
        val = np.flatnonzero(split == "val")
        test = np.flatnonzero(split == "test")
        val_n[celltype], test_n[celltype] = len(val), len(test)
        if len(val) < 5 or len(test) < 3:
            print(f"[all-ct-rgf] skip {name}: val={len(val)} test={len(test)}", flush=True)
            continue
        base_full = donor_base_scores(system, z_all[donors, celltype], celltype, device)
        base_half = np.stack(
            [donor_base_scores(system, z_half_all[half, donors, celltype], celltype, device) for half in (0, 1)]
        )
        common_basis = torch.as_tensor(exact["common_basis"][celltype], dtype=torch.float32, device=device)
        region_gate = torch.as_tensor(exact["region_gate"], dtype=torch.float32, device=device).mean(dim=0)
        score_dictionary = (common_basis * region_gate.unsqueeze(-1)) @ lift
        response_module = torch.as_tensor(exact["response_unit"][donors, celltype], dtype=torch.float32, device=device)
        response_score = torch.einsum("dkm,mg->dkg", response_module, lift)

        normal_full = endpoint_metric(
            system.decoder, thresholds, base_full, score_dictionary, np.zeros(5, dtype=np.float32), None
        )
        normal_half = np.stack(
            [endpoint_metric(system.decoder, thresholds, base_half[half], score_dictionary, np.zeros(5, dtype=np.float32), None) for half in (0, 1)]
        )
        normal_response = endpoint_metric(
            system.decoder, thresholds, base_full, score_dictionary, np.zeros(5, dtype=np.float32), response_score
        )
        fingerprints = np.empty((len(donors), 5, 30), dtype=np.float64)
        half_fingerprints = np.empty((2, len(donors), 5, 30), dtype=np.float64)
        response_fingerprints = np.empty((len(donors), 5, 30), dtype=np.float64)
        shortening_donor = np.empty((len(donors), 5), dtype=np.float64)
        response_shortening_donor = np.empty_like(shortening_donor)
        response_distance_donor = np.empty_like(shortening_donor)

        for axis in range(5):
            endpoint = np.zeros(5, dtype=np.float32)
            endpoint[axis] = 1.0
            disease_full = endpoint_metric(system.decoder, thresholds, base_full, score_dictionary, endpoint, None)
            disease_half = np.stack(
                [endpoint_metric(system.decoder, thresholds, base_half[half], score_dictionary, endpoint, None) for half in (0, 1)]
            )
            disease_response = endpoint_metric(
                system.decoder, thresholds, base_full, score_dictionary, endpoint, response_score
            )
            for donor in range(len(donors)):
                fingerprints[donor, axis] = metric_vector(normal_full[donor], disease_full[donor])
                response_fingerprints[donor, axis] = metric_vector(normal_response[donor], disease_response[donor])
                response_distance_donor[donor, axis] = 0.5 * (
                    log_metric_distance(normal_full[donor], normal_response[donor])
                    + log_metric_distance(disease_full[donor], disease_response[donor])
                )
                for half in (0, 1):
                    half_fingerprints[half, donor, axis] = metric_vector(
                        normal_half[half, donor], disease_half[half, donor]
                    )

            row = celltype * 5 + axis
            fixed = path_saved[row, : path_count[row]]
            direct = np.linspace(0.0, 1.0, 73, dtype=np.float32)[:, None] * endpoint[None]
            for donor in range(len(donors)):
                direct_value = path_length(
                    system.decoder, thresholds, base_full[donor], score_dictionary, direct, None
                )
                fixed_value = path_length(
                    system.decoder, thresholds, base_full[donor], score_dictionary, fixed, None
                )
                shortening_donor[donor, axis] = 100.0 * (direct_value - fixed_value) / max(direct_value, 1e-30)
                response_direct = path_length(
                    system.decoder, thresholds, base_full[donor], score_dictionary, direct, response_score[donor]
                )
                response_fixed = path_length(
                    system.decoder, thresholds, base_full[donor], score_dictionary, fixed, response_score[donor]
                )
                response_shortening_donor[donor, axis] = 100.0 * (response_direct - response_fixed) / max(response_direct, 1e-30)

        template = fingerprints[val].mean(axis=0)
        template /= np.linalg.norm(template, axis=1, keepdims=True).clip(min=1e-30)
        for axis in range(5):
            test_same = np.asarray([correlation(fingerprints[donor, axis], template[axis]) for donor in test])
            test_other = np.asarray(
                [np.median([correlation(fingerprints[donor, axis], template[other]) for other in range(5) if other != axis]) for donor in test]
            )
            half_value = np.asarray(
                [correlation(half_fingerprints[0, donor, axis], half_fingerprints[1, donor, axis]) for donor in range(len(donors))]
            )
            test_similarity[celltype, axis] = np.median(test_same)
            split_half_similarity[celltype, axis] = np.median(half_value)
            axis_specificity[celltype, axis] = np.median(test_same - test_other)
            full5_shortening[celltype, axis] = np.median(shortening_donor[test, axis])
            full5_shorter_fraction[celltype, axis] = np.mean(shortening_donor[test, axis] > 0)
            response_metric_deformation[celltype, axis] = np.median(response_distance_donor[test, axis])
            response_shortening[celltype, axis] = np.median(response_shortening_donor[test, axis])
            ci = bootstrap_ci(test_same, 20260820 + 100 * celltype + axis)
            summary_rows.append(
                {
                    "celltype": name,
                    "pathology": PATHOLOGY[axis],
                    "validation_donors": len(val),
                    "locked_test_donors": len(test),
                    "test_endpoint_metric_similarity_median": ci[0],
                    "test_endpoint_metric_similarity_ci95": ci[1:],
                    "split_half_metric_similarity_median": split_half_similarity[celltype, axis],
                    "test_same_minus_different_axis_median": axis_specificity[celltype, axis],
                    "test_fixed_full5_path_shortening_median_pct": full5_shortening[celltype, axis],
                    "test_fixed_full5_path_shorter_fraction": full5_shorter_fraction[celltype, axis],
                    "test_response_log_metric_deformation_median": response_metric_deformation[celltype, axis],
                    "test_response_on_fixed_path_shortening_median_pct": response_shortening[celltype, axis],
                }
            )
        print(f"[all-ct-rgf] {celltype + 1}/{len(celltypes)} {name} complete", flush=True)

    np.savez_compressed(
        args.out_dir / "all_celltypes_donor_riemann_summary_arrays.npz",
        schema_version=np.asarray("prism.all_celltypes_donor_riemann_summary.v1", dtype=object),
        celltype_names=np.asarray(celltypes, dtype=object),
        pathology_names=np.asarray(PATHOLOGY, dtype=object),
        validation_donor_count=val_n,
        test_donor_count=test_n,
        test_endpoint_metric_similarity=test_similarity,
        split_half_metric_similarity=split_half_similarity,
        test_axis_specificity=axis_specificity,
        test_fixed_full5_shortening_pct=full5_shortening,
        test_fixed_full5_shorter_fraction=full5_shorter_fraction,
        test_response_log_metric_deformation=response_metric_deformation,
        test_response_fixed_path_shortening_pct=response_shortening,
    )

    fig, axes = plt.subplots(2, 3, figsize=(22, 14), constrained_layout=True)
    heatmap(axes[0, 0], test_similarity, celltypes, "A. Test → frozen validation endpoint-metric similarity", "viridis", vmin=0.90, vmax=1.0, fmt=".3f")
    heatmap(axes[0, 1], split_half_similarity, celltypes, "B. Within-donor split-half reliability", "viridis", vmin=0.90, vmax=1.0, fmt=".3f")
    heatmap(axes[0, 2], axis_specificity, celltypes, "C. Same-axis minus different-axis similarity", "YlGnBu", vmin=0.0, fmt=".3f")
    heatmap(axes[1, 0], full5_shortening, celltypes, "D. Fixed Figure-13 full-5D path shortening in test (%)", "magma", vmin=0.0, fmt=".1f")
    heatmap(axes[1, 1], full5_shorter_fraction, celltypes, "E. Test donors where fixed path is shorter", "Blues", vmin=0.0, vmax=1.0, fmt=".2f")
    heatmap(axes[1, 2], response_metric_deformation, celltypes, "F. Exploratory personal-response metric deformation", "Purples", vmin=0.0, fmt=".2f")
    fig.suptitle(
        "Frozen PRISM personal-rank2: donor-conditioned Riemann fingerprint audit across all 24 cell types",
        fontsize=21,
        fontweight="bold",
        color="#173f76",
    )
    fig.savefig(args.out_dir / "figure_06d_all_celltypes_donor_riemann_summary.png", dpi=190, bbox_inches="tight")
    fig.savefig(args.out_dir / "figure_06d_all_celltypes_donor_riemann_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "schema_version": "prism.all_celltypes_donor_riemann_summary.v1",
        "checkpoint_sha256": actual_sha,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "compatibility": compatibility,
        "minimum_cells_per_donor_celltype": args.minimum_cells,
        "validation_donor_count_by_celltype": dict(zip(celltypes, val_n.tolist())),
        "locked_test_donor_count_by_celltype": dict(zip(celltypes, test_n.tolist())),
        "test_refit": False,
        "primary_terms_zeroed": ["explicit personal baseline", "personal-by-pathology response", "sex score", "technology score"],
        "donor_context_retained": "frozen donor-by-celltype zclean centroid",
        "normal_replication_warning": "strict reference donors are absent from locked test; endpoint normal geometry is a donor-context robustness audit, not independent healthy-donor replication",
        "response_warning": "personal-response metrics are exploratory until direct branch ablation supports the response branch",
        "rows": summary_rows,
        "files": {
            "arrays": "all_celltypes_donor_riemann_summary_arrays.npz",
            "figure_png": "figure_06d_all_celltypes_donor_riemann_summary.png",
            "figure_pdf": "figure_06d_all_celltypes_donor_riemann_summary.pdf",
        },
    }
    (args.out_dir / "all_celltypes_summary.json").write_text(
        json.dumps(json_safe(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[all-ct-rgf] wrote {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
