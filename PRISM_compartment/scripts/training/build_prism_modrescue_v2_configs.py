#!/usr/bin/env python3
"""Build PRISM projection-aware/module-rescue-v2 configs without launching.

The Huber-only 2A arm is always emitted.  The CCC 2B arm is emitted only when
``--ccc-weight`` is supplied, because that coefficient must come from the
train-only calibration pilot rather than validation or test performance.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "sv6": ROOT
    / "configs/train_config_kmlee_bam_dlpfc_mtg_bins4_"
    "prism_modrescue_recovery_v1_nointeraction_s42_sv6.json",
    "sv7": ROOT
    / "configs/train_config_kmlee_bam_dlpfc_mtg_bins4_"
    "prism_modrescue_recovery_v1_interaction_s42_sv7.json",
}

# Frozen paired-input provenance, measured independently on sv6 and sv7 on
# 2026-08-04.  Both servers must contain byte-identical copies before launch.
WARM_START_SHA256 = (
    "548b0a0fae8beb5c3f2e312384d0ea781a91e13af8927f1375044186d1f9c3e0"
)
MODULE_RESCUE_STATS_SHA256 = (
    "c1d70a9dc88ee78fad047b2a4ebbf723898ebc3cf5dafb639a35b3cbc1d08147"
)


def _build(
    source: Path,
    *,
    server: str,
    arm: str,
    ccc_weight: float,
    ccc_calibration_sha256: str | None,
) -> dict:
    config = copy.deepcopy(json.loads(source.read_text(encoding="utf-8")))
    stage_description = (
        "Huber-only Stage 2A package"
        if float(ccc_weight) == 0.0
        else "train-only calibrated Huber+CCC Stage 2B package"
    )
    config["_doc_pointer"] = (
        "doc/prism_projection_whitening_module_rescue_v2_design_20260804.md "
        "§§3-6; Stage 1 is discharged by the analytic/numerical gradient-"
        f"equivalence gate, and this config is the {stage_description}."
    )
    # Only train.out_dir is consumed.  Remove the inherited legacy duplicate so
    # provenance tools cannot mistake a 20260729 arm for this experiment.
    config.pop("out_dir", None)
    server_root = (
        "/home/kmlee/project/riemann_bam"
        if server == "sv7"
        else "/home/kmlee/project/kmlee_bam"
    )
    config["train"].update(
        {
            "epochs": 12,
            "save_every": 1,
            "save_best": False,
            "early_stopping": False,
            "reload_best_before_test": False,
            "skip_final_test": True,
            "out_dir": (
                f"{server_root}/outputs/SEA_AD_OUTPUTS/train_runs/"
                f"kmlee_bam_prism_modrescue_v2{arm}_nointeraction_"
                "s42_20260804"
            ),
        }
    )
    config["precision_medicine"]["interaction_enabled"] = False
    config["warm_start"].update(
        {
            "expected_sha256": WARM_START_SHA256,
            "require_precision_head": True,
        }
    )
    config["nuisance_projection"]["projection_aware_whitening"] = True
    config["prism_module_rescue"].update(
        {
            "enabled": True,
            "expected_stats_sha256": MODULE_RESCUE_STATS_SHA256,
            "lambda_module": 0.10,
            "branch_fraction": 0.25,
            "huber_delta": 1.0,
            "ccc_weight": float(ccc_weight),
            "steps_per_epoch": 24 if server == "sv7" else 16,
            "draws_per_celltype_per_epoch": 4,
            "cells_per_donor_celltype": 2,
            "minimum_donors_per_celltype": 8,
            "minimum_target_cells": 10,
            "warmup_epochs": 2,
            "nonoverlap_sampling": True,
            "deterministic_forward": True,
            "restrict_parameter_updates": True,
            "contract_version": "v2",
            "global_blocks_per_optimizer_step": 12,
            "equal_donor_rescue_gradient": True,
            "ccc_target_variance_floor": 1e-4,
            "ccc_calibration_sha256": (
                ccc_calibration_sha256
            ),
            # Hold the recovery-v1 sampler seed fixed.  The scheduler changes,
            # but an unnecessary random-seed change must not become a confound.
            "seed": 420803,
        }
    )
    config["learned_generator_count"]["enabled"] = False
    config["generator_budget_search"]["enabled"] = False
    config["experiment_manifest"] = {
        "schema_version": f"kmlee_bam.prism_modrescue_v2{arm}",
        "server": server,
        "projection_aware_whitening": True,
        "effective_latent_dimension_when_sex_rank4": 28,
        "generator_bank_all_active": True,
        "generator_count_search": False,
        "pathology_interaction": False,
        "module_rescue_global_blocks_per_epoch": 96,
        "module_rescue_draws_per_celltype_per_epoch": 4,
        "module_rescue_cells_per_donor_per_draw": 2,
        "module_rescue_nonoverlap_sampling": True,
        "module_rescue_deterministic_posterior_mean": True,
        "module_rescue_restricted_parameter_updates": True,
        "module_rescue_global_blocks_per_optimizer_step": 12,
        "module_rescue_optimizer_updates_per_epoch": 8,
        "module_rescue_equal_donor_gradient": True,
        "module_rescue_parameter_scope": (
            "decoder state heads plus explicit common/personal/response/"
            "interaction output dictionaries; context encoder and age/region "
            "baselines frozen during rescue-only steps"
        ),
        "module_rescue_branch_fraction": 0.25,
        "module_rescue_full_fraction": 0.75,
        "module_rescue_ccc_weight": float(ccc_weight),
        "warm_start_sha256": WARM_START_SHA256,
        "module_rescue_stats_sha256": MODULE_RESCUE_STATS_SHA256,
        "whitening_stage1_validation": "analytic_plus_bf16_regression_smoke",
        "whitening_expected_identity_loss_floor_removed": 4.0 / (32.0**2),
        "runtime_manifest_required": True,
        "official_test_opened": False,
        "checkpoint_selection": "validation_donor_centered_module_recovery",
    }
    return config


def _smoke(production: dict, *, server: str, arm: str) -> dict:
    config = copy.deepcopy(production)
    config["train"].update(
        {
            "epochs": 1,
            "max_train_steps_per_epoch": 1,
            "max_eval_steps": 1,
            "progress_log_every": 1,
            "out_dir": f"/tmp/prism_modrescue_v2{arm}_{server}_ddp_smoke_r2",
        }
    )
    config["loader"].update({"batch_size": 16, "eval_batch_size": 16})
    config["prism_module_rescue"].update(
        {
            # Production-topology smoke: one global draw for all 24 cell types.
            # sv6 uses 4 local slots x 6 ranks and sv7 uses 6 x 4; both execute
            # two identical 12-global-block optimizer updates.
            "steps_per_epoch": 6 if server == "sv7" else 4,
            "draws_per_celltype_per_epoch": 1,
            "contract_version": "v2",
            "global_blocks_per_optimizer_step": 12,
            "warmup_epochs": 1,
        }
    )
    config["experiment_manifest"].update(
        {
            "schema_version": f"kmlee_bam.prism_modrescue_v2{arm}.smoke",
            "smoke_only": True,
            "smoke_world_size": 4 if server == "sv7" else 6,
            "smoke_kind": "production_topology_ddp",
            "module_rescue_global_blocks_per_epoch": 24,
            "module_rescue_draws_per_celltype_per_epoch": 1,
            "module_rescue_optimizer_updates_per_epoch": 2,
        }
    )
    return config


def _calibration(production: dict, *, server: str) -> dict:
    """Single-GPU, train-only input for CCC gradient calibration."""

    config = copy.deepcopy(production)
    config["train"].update(
        {
            "epochs": 1,
            "max_train_steps_per_epoch": 1,
            "max_eval_steps": 0,
            "out_dir": f"/tmp/prism_modrescue_v2_ccc_calibration_{server}_r3",
        }
    )
    # The warm-start immediately restores decoder thresholds, and calibration
    # consumes only the module-rescue objective.  Avoid spending minutes on an
    # 8,192-cell cold-start threshold estimate that is discarded before the
    # first calibration block.
    config["decoder"]["bin_marginal_n_sample"] = 128
    config["v4"]["ordinal_balance"]["max_cells_for_bin_stats"] = 128
    config["prism_module_rescue"].update(
        {
            "ccc_weight": 0.0,
            "ccc_calibration_sha256": None,
            # Calibration measures exactly one deterministic train-only block
            # for each of the 24 eligible cell types on one GPU.  Do not retain
            # the four-GPU production draw contract here.
            "steps_per_epoch": 24,
            "draws_per_celltype_per_epoch": 1,
        }
    )
    config["experiment_manifest"].update(
        {
            "schema_version": (
                "kmlee_bam.prism_modrescue_v2.ccc_gradient_calibration"
            ),
            "calibration_only": True,
            "calibration_source_split": "train_only",
            "calibration_validation_or_test_used": False,
            "calibration_blocks": 24,
            "calibration_world_size": 1,
        }
    )
    return config


def _write(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ccc-weight",
        type=float,
        default=None,
        help=(
            "Train-only calibrated CCC coefficient. Omit to emit only the "
            "Huber-only 2A arm."
        ),
    )
    parser.add_argument(
        "--ccc-calibration-sha256",
        default=None,
        help="SHA256 of the train-only CCC gradient calibration record.",
    )
    args = parser.parse_args()
    if args.ccc_weight is not None and (
        not math.isfinite(args.ccc_weight) or args.ccc_weight <= 0.0
    ):
        parser.error("--ccc-weight must be finite and positive")
    if args.ccc_weight is not None:
        digest = str(args.ccc_calibration_sha256 or "").lower()
        if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            parser.error(
                "--ccc-calibration-sha256 must be a 64-character hex digest "
                "when --ccc-weight is supplied"
            )

    arms = [("a_huber", 0.0)]
    if args.ccc_weight is not None:
        arms.append(("b_ccc", float(args.ccc_weight)))
    for server, source in SOURCES.items():
        for arm, ccc_weight in arms:
            production = _build(
                source,
                server=server,
                arm=arm,
                ccc_weight=ccc_weight,
                ccc_calibration_sha256=(
                    None
                    if float(ccc_weight) == 0.0
                    else str(args.ccc_calibration_sha256).lower()
                ),
            )
            stem = (
                "train_config_kmlee_bam_dlpfc_mtg_bins4_"
                f"prism_modrescue_v2{arm}_nointeraction_s42_{server}"
            )
            output = ROOT / "configs" / f"{stem}.json"
            smoke = ROOT / "configs" / f"_smoke_{stem}.json"
            _write(output, production)
            _write(smoke, _smoke(production, server=server, arm=arm))
            if arm == "a_huber":
                calibration = ROOT / "configs" / (
                    "_calibration_train_config_kmlee_bam_dlpfc_mtg_bins4_"
                    f"prism_modrescue_v2_ccc_nointeraction_s42_{server}.json"
                )
                _write(calibration, _calibration(production, server=server))


if __name__ == "__main__":
    main()
