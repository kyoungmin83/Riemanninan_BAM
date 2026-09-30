#!/usr/bin/env python3
"""Build the authorized SV6 integrated rank-first compartmental PRISM config."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = (
    ROOT
    / "configs/final/train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
)
BLUEPRINT = (
    ROOT
    / "configs/final/prism_integrated_compartmental_nonlinearity_techzero_v2_blueprint_20260827.json"
)
RUN_NAME = "kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_20260828"
RUN_DIR = (
    "/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
    + RUN_NAME
)
RUNTIME = "/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828"
OUTPUT_CAP_CHECKPOINT_SHA256 = (
    "2cdc91ae005e6183cd7aa9f620ab463e0c348db9752247adf320ecc261a99bd0"
)
OUTPUT_CAP_CONFIG_SHA256 = (
    "b9428e1f1a71304fd26b90ae5527ca6b2eb9d0686bde49a7853082e65924b420"
)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: str, expected_sha256: str | None = None) -> str:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    observed = sha256_file(source)
    if expected_sha256 is not None and observed != expected_sha256:
        raise RuntimeError(f"artifact SHA256 mismatch for {source}: {observed}")
    return observed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--reliability", required=True)
    parser.add_argument("--reliability-sha256", required=True)
    parser.add_argument("--compartment-graph", required=True)
    parser.add_argument("--compartment-graph-sha256", required=True)
    parser.add_argument("--output-cap", required=True)
    parser.add_argument("--output-cap-sha256", required=True)
    parser.add_argument("--module-rescue", required=True)
    parser.add_argument("--module-rescue-sha256", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = copy.deepcopy(json.loads(BASE.read_text(encoding="utf-8")))
    blueprint = json.loads(BLUEPRINT.read_text(encoding="utf-8"))
    artifact_hashes = {
        "module_local_reliability_sha256": require_file(
            args.reliability, args.reliability_sha256
        ),
        "module_local_compartment_graph_sha256": require_file(
            args.compartment_graph, args.compartment_graph_sha256
        ),
        "module_local_output_cap_sha256": require_file(
            args.output_cap, args.output_cap_sha256
        ),
        "module_rescue_stats_sha256": require_file(
            args.module_rescue, args.module_rescue_sha256
        ),
    }

    config["train"].update(
        {
            "epochs": 55,
            "out_dir": RUN_DIR,
            "early_stopping": False,
            "early_stopping_min_epochs": 55,
            "require_validation_for_early_stopping": False,
            "reload_best_before_test": False,
            "skip_final_test": True,
            "resume_checkpoint": None,
            "resume_optimizer": True,
            "resume_history": True,
            "log_every": 100,
            "progress_log_every": 500,
            "eval_log_every": 2000,
            "grad_accum_steps": 1,
        }
    )
    config["loader"].update({"batch_size": 64, "eval_batch_size": 32})
    config["optim"]["lr"] = 1.0e-4
    config["decoder"].update(
        {
            "use_pathology_decoder": True,
            "n_pathology_axes": 4,
            "pathology_pairwise": False,
            "use_pathology_interactions": False,
            "pathology_rank": 0,
        }
    )
    config["learned_pathology_rank"] = {
        "enabled": True,
        "initial_keep_probability": 0.995,
        "warmup_fixed_rank": 8,
        "warmup_end_epoch": 24,
        "initial_extra_keep_probability": 0.5,
        "require_canonical_phase1_rank8": True,
        "require_active_pathology_route_during_search": True,
        "temperature_start": 2.0,
        "temperature_end": 0.35,
        "soft_start_epoch": 25,
        "hard_start_epoch": 33,
        "freeze_epoch": 37,
        "sparsity_start_epoch": 29,
        "sparsity_ramp_epochs": 4,
        "cardinality_loss_fraction": 0.0025,
        "learning_rate": 3.0e-5,
        "weight_decay": 0.0,
        "hard_threshold": 0.5,
        "require_stable_freeze": True,
        "freeze_threshold_low": 0.45,
        "freeze_threshold_high": 0.55,
        "maximum_freeze_count_spread": 4,
        "maximum_freeze_uncertain_fraction": 0.15,
    }

    generator = config["learned_generator_count"]
    generator.update(
        {
            "enabled": True,
            "mode": "joint",
            "initial_keep_probability": 0.995,
            "minimum_active_generators": None,
            "protect_singletons": True,
            "protect_unique_gene_coverage": True,
            "start_epoch": 37,
            "shadow_end_epoch": 40,
            "soft_end_epoch": 44,
            "base_end_epoch": 55,
            "safe_mask_confirmations": 3,
            "rollback_to_last_safe_mask": True,
            "require_frozen_pathology_rank_before_search": True,
        }
    )
    config["generator_budget_search"] = {"enabled": False}

    curriculum = config["integrated_phase_curriculum"]
    curriculum.update(
        {
            "enabled": True,
            "phase1_end_epoch": 12,
            "phase2_start_epoch": 13,
            "pathology_crossfade_end_epoch": 16,
            "phase1_learning_rate": 1.0e-4,
            "phase2_learning_rate": 2.0e-5,
            "phase1_grad_accum_steps": 1,
            "phase2_grad_accum_steps": 2,
            "phase1_batch_size": 64,
            "phase2_batch_size": 32,
            "phase1_pathology_scale": 1.0,
            "phase2_pathology_scale": 1.0,
            "skip_dormant_precision_forward": True,
            "require_exact_phase1_contract": True,
        }
    )

    precision = config["precision_medicine"]
    precision.update(
        {
            "context_npz_path": (
                "/home/kmlee/project_sv6/kmlee_bam/outputs/"
                "population_fingerprint_20260714/"
                "ctemb64_epoch12_full_region_pseudobulk.npz"
            ),
            "enabled": True,
            "personal_rank": 2,
            "loss_start_epoch": 13,
            "output_start_epoch": 13,
            "loss_ramp_epochs": 3,
            "output_ramp_epochs": 3,
            "scale_regularizers_with_loss_ramp": True,
            "module_local_enabled": True,
            "module_local_reliability_npz_path": str(Path(args.reliability).resolve()),
            "module_local_rank": 8,
            "module_local_start_epoch": 17,
            "module_local_ramp_epochs": 8,
            "module_local_rescue_start_epoch": 25,
            "module_local_output_cap_npz_path": str(Path(args.output_cap).resolve()),
            "module_local_output_cap_source_checkpoint_sha256": OUTPUT_CAP_CHECKPOINT_SHA256,
            "module_local_output_cap_source_config_sha256": OUTPUT_CAP_CONFIG_SHA256,
            "module_local_output_cap_quantile": 0.995,
            "module_local_input_clip_quantile": 0.995,
            "module_local_projection_ridge": 1.0e-4,
            "module_local_size_ratio_cap": 1.0,
            "module_local_lr_multiplier": 10.0,
            "lambda_module_local_center": 2.0e-3,
            "lambda_module_local_hierarchy": 1.0e-4,
            "lambda_module_local_size": 2.0e-3,
            "lambda_module_local_pathology_leak": 1.0e-2,
            "lambda_module_local_age_leak": 5.0e-3,
            "module_local_nonlinear_enabled": True,
            "module_local_nonlinear_variant": "compartmental_threshold",
            "module_local_compartment_graph_npz_path": str(
                Path(args.compartment_graph).resolve()
            ),
            "module_local_nonlinear_start_epoch": 21,
            "module_local_nonlinear_ramp_epochs": 4,
            "module_local_compartment_mix_max": 0.5,
            "module_local_threshold_min": 0.5,
            "module_local_threshold_max": 2.0,
            "module_local_slope_min": 1.0,
            "module_local_slope_max": 8.0,
            "module_local_gain_max": 1.0,
        }
    )

    rescue = config["prism_module_rescue"]
    rescue.update(
        {
            "enabled": True,
            "stats_path": str(Path(args.module_rescue).resolve()),
            "expected_stats_sha256": args.module_rescue_sha256,
            "start_epoch": 25,
            "warmup_epochs": 4,
            "allow_state_readout_updates": False,
            "pathology_axis_enabled": True,
            "pathology_axis_start_epoch": 25,
            "celltype_tail_enabled": True,
            "celltype_tail_start_epoch": 25,
        }
    )

    pi_blueprint = blueprint["pi_tech"]
    config["pi_tech"] = {
        "enabled": True,
        "mode": "same_celltype_donor_bank",
        "k": 16,
        "in_batch_k": 8,
        "bank_start_epoch": 13,
        "bank_ramp_epochs": 4,
        "bank_refresh_epochs": 1,
        "bank_per_donor_celltype": 2,
        "bank_batch_size": 48,
        "bank_num_workers": 0,
        "bank_min_distinct_donors": 4,
        "bank_detect_shrinkage_donors": 8.0,
        "bank_seed": 20260827,
        "bank_exclude_query_donor": True,
        "bank_invalid_fallback": "zero",
        "bank_match_region": True,
        "diagnostic_min_positions": 256,
        "diagnostic_min_cells": 8,
        "diagnostic_min_distinct_donors": 2,
        "sex_linked_gene_names": pi_blueprint["sex_linked_gene_names"],
        "allow_missing_sex_linked_genes": False,
        "sex_linked_expected_resolution": pi_blueprint[
            "sex_linked_expected_resolution"
        ],
        "pi_cap": 0.7,
        "depth_proxy": "n_detected",
        "depth_low_pct": 0.25,
        "depth_low_temp": 0.05,
        "detect_floor": 0.05,
        "neighbor_on_hi": 0.5,
        "safety_neighbor_off": 0.2,
        "w_measure": 0.5,
        "w_neighbor": 0.5,
        "pi_soft_target_zero": 0.5,
        "pi_downweight_zero": 0.7,
        "eps": 1.0e-6,
    }

    critical = {
        "precision_context_npz_sha256": require_file(
            precision["context_npz_path"]
        ),
        "module_registry_json_sha256": require_file(
            config["module_tokenizer"]["registry_json_path"]
        ),
        "activity_weight_npz_sha256": require_file(
            config["module_tokenizer"]["activity_weight_path"]
        ),
        **artifact_hashes,
    }
    config["experiment_manifest"] = {
        "schema_version": "kmlee_bam.prism_integrated_compartmental_pathrank.v1",
        "run_name": RUN_NAME,
        "runtime_root": RUNTIME,
        "host": "SV6",
        "training_form": "one_architecture_one_optimizer_one_checkpoint_lineage",
        "initialization": "scratch",
        "canonical_phase1_source": str(BASE),
        "design_blueprint": str(BLUEPRINT),
        "phase1_exact_epochs": [1, 12],
        "phase1_effective_pathology_rank": 8,
        "pathology_rank_allocation": "automatic_algebraic_capacity_96",
        "pathology_rank_learning": "soft25_hard33_freeze37",
        "generator_cardinality_learning": "shadow37_40_soft41_44_hard45_55",
        "rank_generator_cardinality_overlap_forbidden": True,
        "pathology_route_scale_all_epochs": 1.0,
        "personal_rank": 2,
        "module_registry_count": 414,
        "paper_source": blueprint["paper_translation"]["source_pdf"],
        "paper_translation": "computational_analogy_not_observed_dendritic_biophysics",
        "paper_arm": "compartmental_threshold",
        "official_legacy_test_donors": "previously_consumed_not_model_selection_eligible",
        "official_test_dataset_opened_by_training": False,
        "confirmatory_holdout_status": "required_before_final_release_not_available_at_training_launch",
        "critical_model_inputs": critical,
    }
    config["_launch_guard"] = {
        "launch_allowed": True,
        "authorized_design_answer": "integrated",
        "explicit_launch_approval": "2026-08-28 user requested tmux training on SV6",
        "phase1_epoch12_parity_preflight_required": True,
        "resume_equivalence_preflight_required": True,
        "merged_curriculum_smoke_required": True,
        "exact_budget_audit_required": True,
        "official_test_dataset_open_forbidden": True,
        "rank_commit_before_generator_required": True,
        "new_confirmatory_holdout_required_before_final_release": True,
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"wrote {output} sha256={sha256_file(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
