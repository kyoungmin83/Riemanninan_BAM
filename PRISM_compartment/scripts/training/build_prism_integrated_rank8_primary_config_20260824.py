#!/usr/bin/env python3
"""Build the executable one-lineage Phase-I -> Phase-II rank8 config."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE1 = ROOT / "configs/phase1/train_config_kmlee_bam_dlpfc_mtg_bins4_v31a_cons_h8_ct64_gen414_singletonr2_s42.json"
PHASE2 = ROOT / "configs/final/resolved_config_epoch20.json"
BLUEPRINT = ROOT / "configs/final/prism_integrated_celltype_ad_module_curriculum_v2_blueprint_20260824.json"
OUTPUT = ROOT / "configs/final/train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"

REMOTE_RUNTIME = "/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_20260824"
RUN_NAME = "kmlee_bam_prism_integrated_rank8_primary_s42_sv6_retry1_20260824"
RUN_DIR = f"/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/{RUN_NAME}"
ARTIFACT = f"{REMOTE_RUNTIME}/artifacts/prism_module_rescue_stats_train64_pathaxis_v3_20260824.npz"
ARTIFACT_SHA256 = "b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c"


def main() -> None:
    phase1 = json.loads(PHASE1.read_text(encoding="utf-8"))
    phase2 = json.loads(PHASE2.read_text(encoding="utf-8"))
    blueprint = json.loads(BLUEPRINT.read_text(encoding="utf-8"))
    config = copy.deepcopy(phase1)

    # The epoch-1 surface remains canonical Phase I.
    config["loader"]["batch_size"] = 64
    config["optim"]["lr"] = 1.0e-4
    config["train"].update(
        {
            "epochs": 55,
            "out_dir": RUN_DIR,
            "grad_accum_steps": 1,
            "early_stopping": False,
            "early_stopping_min_epochs": 55,
            "require_validation_for_early_stopping": False,
            "reload_best_before_test": False,
            "skip_final_test": True,
            "ddp_skip_initial_sync": True,
            "resume_checkpoint": None,
            "resume_optimizer": True,
            "resume_history": True,
        }
    )
    config["nuisance_projection"]["projection_aware_whitening"] = False

    # Phase II modules exist in the checkpoint from epoch 1, but their forward
    # path is skipped until global epoch 13 by the integrated controller.
    precision = copy.deepcopy(phase2["precision_medicine"])
    precision.update(
        {
            "enabled": True,
            "personal_rank": 2,
            "loss_start_epoch": 13,
            "output_start_epoch": 13,
            "loss_ramp_epochs": 3,
            "output_ramp_epochs": 3,
            "scale_regularizers_with_loss_ramp": True,
        }
    )
    config["precision_medicine"] = precision

    config["v7a"] = {
        "pathology_aware_uncertainty": {
            "enabled": True,
            "alignment": {
                "enabled": True,
                "lambda_alignment": 0.005,
                "detach_target": True,
                "loss_form": "huber",
                "huber_delta": 1.0,
                "alignment_start_epoch": 17,
                "ramp_epochs": 4,
                "u_total_cap": 5.0,
                "lambda_relative_rec": 0.5,
                "tau_relative": 0.25,
                "relative_start_epoch": 25,
                "relative_ramp_epochs": 5,
                "relative_weight_min": 0.8,
                "relative_weight_max": 1.25,
                "relative_weight_source": "target",
                "detach_relative_source": True,
                "normalize_relative_within_celltype": True,
                "minimum_relative_ess_fraction": 0.9,
                "eps": 1.0e-6,
            },
        }
    }

    rescue = copy.deepcopy(phase2["prism_module_rescue"])
    rescue.update(
        {
            "enabled": True,
            "stats_path": ARTIFACT,
            "expected_stats_sha256": ARTIFACT_SHA256,
            "start_epoch": 25,
            "warmup_epochs": 4,
            "non_neuronal_pseudobulk_mode": "aggregate",
            "non_neuronal_aggregate_draws": 4,
            # The primary architecture uses CLS pooling.  The rescue
            # state-readout allowlist targets persistent AGP pooling modules
            # and must remain disabled for CLS checkpoints.
            "allow_state_readout_updates": False,
            "state_readout_start_epoch": 29,
            "pathology_axis_enabled": True,
            "pathology_axis_start_epoch": 25,
            "lambda_pathology_axis": 0.05,
            "pathology_axis_magnitude_weight": 0.1,
            "celltype_tail_enabled": True,
            "celltype_tail_start_epoch": 25,
            "celltype_deficit_ema_decay": 0.8,
            "celltype_weight_min": 0.5,
            "celltype_weight_max": 2.0,
            "celltype_tail_fraction": 0.25,
            "ad_recovery_baseline_by_celltype": blueprint[
                "per_celltype_ad_module_baseline"
            ],
        }
    )
    config["prism_module_rescue"] = rescue

    gate = copy.deepcopy(phase2["learned_generator_count"])
    gate.update(
        {
            "enabled": True,
            "mode": "joint",
            "minimum_active_generators": None,
            "protect_singletons": True,
            "protect_unique_gene_coverage": True,
            # Preserve the exact Phase-I forward through epoch 12, then begin
            # cardinality learning immediately.  Shadow search cannot alter
            # the live decoder; soft gates begin after the cross-fade, and
            # irreversible hard masks remain protected until epoch 25.
            "start_epoch": 13,
            "shadow_end_epoch": 16,
            "soft_end_epoch": 24,
            "base_end_epoch": 55,
            "safe_mask_confirmations": 3,
            "rollback_to_last_safe_mask": True,
        }
    )
    config["learned_generator_count"] = gate
    config["generator_budget_search"] = {"enabled": False}

    config["integrated_phase_curriculum"] = {
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
        "phase2_pathology_scale": 0.0,
        "phase1_lambda_state_abs": 0.08,
        "phase1_lambda_state_fraction": 0.01,
        "phase2_lambda_state_abs": 0.0,
        "phase2_lambda_state_fraction": 0.0,
        "phase1_projection_aware_whitening": False,
        "phase2_projection_aware_whitening": True,
        "phase1_phu_alignment_lambda": 0.001,
        "phase1_phu_alignment_start_epoch": 1,
        "phase1_phu_alignment_ramp_epochs": 5,
        "phase2_phu_alignment_lambda": 0.005,
        "phase2_phu_alignment_start_epoch": 17,
        "phase2_phu_alignment_ramp_epochs": 4,
        "phase1_phu_relative_lambda": 0.0,
        "phase2_phu_relative_lambda": 0.5,
        "phase2_phu_relative_start_epoch": 25,
        "phase2_phu_relative_ramp_epochs": 5,
        "skip_dormant_precision_forward": True,
        "require_exact_phase1_contract": True,
    }

    config["experiment_manifest"] = {
        "schema_version": "kmlee_bam.prism_integrated_rank8_primary.v1",
        "run_name": RUN_NAME,
        "runtime_root": REMOTE_RUNTIME,
        "training_form": "one_architecture_one_optimizer_one_checkpoint_lineage",
        "initialization": "scratch",
        "retry_of": "kmlee_bam_prism_integrated_rank8_primary_s42_sv6_20260824",
        "retry_reason": "pre_optimizer_FrozenInstanceError_in_curriculum_apply_epoch",
        "phase1_source": str(PHASE1),
        "phase2_source": str(PHASE2),
        "design_blueprint": str(BLUEPRINT),
        "phase1_exact_epochs": [1, 12],
        "phase2_start_epoch": 13,
        "pathology_rank": 8,
        "personal_rank": 2,
        "module_registry_count": 414,
        "matched_separated_budget_checkpoint_epoch": 32,
        "extension_end_epoch": 55,
        "official_test_donors": "sealed_not_opened",
        "pathology_axis_artifact_sha256": ARTIFACT_SHA256,
    }
    config["_launch_guard"] = {
        "launch_allowed": True,
        "authorized_design_answer": "integrated",
        "phase1_parity_preflight_required": True,
        "resume_preflight_required": True,
        "compressed_stage_smoke_required": True,
        "official_test_dataset_open_forbidden": True,
    }
    OUTPUT.write_text(
        json.dumps(config, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
