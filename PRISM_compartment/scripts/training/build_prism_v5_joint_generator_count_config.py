#!/usr/bin/env python3
"""Build the sv6 PRISM-v5 online learned-generator-count full config.

The source is the audited v4 recovery configuration.  This builder changes
only the run identity/training horizon and the learned-count controller.  It
does not create W44/R10/C10 partitions and never emits a candidate-K grid.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path


DEFAULT_SOURCE = Path(
    "configs/train_config_prism_v4_sv6_lineage_ccc_fasttrack_20260805.json"
)
DEFAULT_OUTPUT = Path(
    "configs/train_config_prism_v5_sv6_joint_generator_count_s42_20260805.json"
)
DEFAULT_SMOKE_OUTPUT = Path(
    "configs/_smoke_train_config_prism_v5_sv6_joint_generator_count_s42_20260805.json"
)


def build(source: Path, output: Path, smoke_output: Path) -> None:
    config = json.loads(source.read_text(encoding="utf-8"))

    train = config["train"]
    train["out_dir"] = (
        "/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
        "kmlee_bam_prism_v5_sv6_joint_l0_generator_count_s42_20260805"
    )
    train["epochs"] = 36
    train["early_stopping"] = True
    train["early_stopping_patience"] = 8
    train["early_stopping_min_epochs"] = 18
    train["require_validation_for_early_stopping"] = True
    train["skip_final_test"] = True

    criteria = list(train.get("early_stopping_criteria", []))
    criteria.extend(
        [
            {
                "name": "metric/generator_candidate_full_nll",
                "mode": "min",
                "min_delta": 0.0005,
                "weight": 0.5,
                "short_name": "gate_full",
            },
            {
                "name": "metric/generator_candidate_isolated_nll",
                "mode": "min",
                "min_delta": 0.0005,
                "weight": 0.5,
                "short_name": "gate_isolated",
            },
            {
                "name": "metric/generator_expected_active",
                "mode": "min",
                "min_delta": 0.5,
                "weight": 0.25,
                "short_name": "gate_count",
            },
        ]
    )
    train["early_stopping_criteria"] = criteria

    config["generator_budget_search"]["enabled"] = False
    learned = config["learned_generator_count"]
    learned.update(
        {
            "enabled": True,
            "mode": "joint",
            "initial_keep_probability": 0.995,
            "temperature_start": 2.0,
            "temperature_end": 0.3,
            "hard_threshold": 0.5,
            "minimum_active_generators": 32,
            "protect_singletons": True,
            "protect_unique_gene_coverage": False,
            "augmented_rho": 2.0,
            "dual_learning_rate": 0.0001,
            "constraint_ema": 0.99,
            "start_epoch": 3,
            "base_end_epoch": 18,
            "learning_rate": 0.0002,
            "weight_decay": 0.0,
            "joint_objective_weight": 0.05,
            "joint_full_nll_relative_margin": 0.01,
            "joint_isolated_nll_relative_margin": 0.015,
            "joint_constraint_scale_floor": 0.001,
        }
    )
    # These keys belong exclusively to the retired gate-only W44/R10/C10
    # controller.  Joint mode ignores them in code, but removing them from the
    # reviewed artifact makes it impossible to mistake this run for an outer
    # candidate-K search.
    for legacy_key in (
        "coalition_samples",
        "max_extension_epochs",
        "architecture_donor_count",
        "architecture_split_seed",
        "architecture_split_trials",
        "gate_updates_per_epoch",
        "minimum_gate_updates",
        "steps_per_epoch",
        "architecture_draws_per_celltype",
        "cells_per_donor_celltype",
        "minimum_donors_per_celltype",
        "bootstrap_replicates",
        "noninferiority_se_multiplier",
        "constraint_scale_floor",
        "stability_window",
        "stable_count_tolerance",
        "stable_mask_jaccard",
        "maximum_uncertain_fraction",
        "minimum_search_epochs",
        "minimum_pruned_generators",
        "seed",
    ):
        learned.pop(legacy_key, None)

    design_path = Path(
        "doc/prism_v5_online_joint_generator_count_design_20260805.md"
    )
    runtime_sources = (
        "src/kmlee_bam/model/latent_nuisance_projection.py",
        "src/kmlee_bam/model/lie_ordinal_decoder.py",
        "src/kmlee_bam/model/system.py",
        "src/kmlee_bam/objectives/celltype_alignment.py",
        "src/kmlee_bam/training/core_trainer.py",
        "src/kmlee_bam/training/learned_generator_count.py",
        "src/kmlee_bam/training/prism_module_rescue_training.py",
        "src/kmlee_bam/training/run_current.py",
        "src/kmlee_bam/training/runner_base.py",
    )

    config["_doc_pointer"] = {
        "design": str(design_path),
        "selection": "online exact-hard Binary-Concrete; no candidate-K grid",
        "base_checkpoint": (
            "PRISM v2a validation-selected epoch 11, "
            "SHA256 9bee5aba99e38453499d6cd2333831509b811945a33e77c2b381eb3ea655c555"
        ),
    }
    manifest = config["experiment_manifest"]
    manifest.update(
        {
            "schema_version": "kmlee_bam.prism_v5_joint_generator_count.v1",
            "purpose": "online_joint_optimization_of_decoder_generator_count",
            "design_document": str(design_path),
            "design_document_sha256": hashlib.sha256(
                design_path.read_bytes()
            ).hexdigest(),
            "scientific_changes_from_selected_v2a": [
                "module_rescue_ccc_weight_0_to_1",
                "celltype_loss_aggregation_uniform_to_lineage_equal_mass",
                "online_joint_exact_hard_generator_gates",
                "paired_full_and_generator_isolated_noninferiority",
                "module_rescue_straight_through_feedback_to_generator_gate",
            ],
            "training_control": {
                "maximum_epochs": 36,
                "minimum_epochs_before_stop": 18,
                "early_stopping_patience": 8,
                "monitor_split": "val",
                "rank0_exact_validation": True,
                "save_every_epoch": True,
                "test_used": False,
            },
            "generator_count_contract": {
                "mode": "joint",
                "candidate_count": 414,
                "candidate_k_grid_used": False,
                "train_donor_count": 64,
                "warmup_all_on_epochs": 2,
                "minimum_active_generators": 32,
                "temperature_start": 2.0,
                "temperature_end": 0.3,
                "anneal_end_epoch": 18,
                "full_nll_relative_margin": 0.01,
                "isolated_nll_relative_margin": 0.015,
                "protected_singleton_count": 10,
                "unique_gene_coverage_blanket_protection": False,
                "checkpoint_selection_requires_val_full_violation_le_zero": True,
                "checkpoint_selection_requires_val_isolated_violation_le_zero": True,
                "no_feasible_pruned_epoch_fallback": (
                    "last_best_feasible_all_on_or_larger_K"
                ),
            },
            "ccc_calibration_scope": {
                "evidence_sha256": config["prism_module_rescue"][
                    "ccc_calibration_sha256"
                ],
                "calibrated_weight": 1.0,
                "inherited_from_v2_allowlist": True,
                "calibration_excluded_new_joint_gate_parameter": True,
                "raw_v2_weight_exceeded_clip_ceiling": True,
                "claim": (
                    "The inherited train-only calibration justifies the clipped "
                    "CCC coefficient for the original rescue parameters; it is "
                    "not presented as a recalibration of the new gate logit."
                ),
            },
            "runtime_environment_contract": {
                "WORLD_SIZE": 6,
                "KMLEE_GROUPED_BLOCK_SIZE": 8,
                "KMLEE_PATH_CONSISTENCY": 0,
                "KMLEE_PATH_CONS_MINCELLS": 2,
            },
            "runtime_source_sha256": {
                relative: hashlib.sha256(Path(relative).read_bytes()).hexdigest()
                for relative in runtime_sources
            },
            "full_training_authorized": False,
            "requires_explicit_user_approval": True,
            "launch_performed": False,
        }
    )
    config["_launch_guard"] = {
        "launch_allowed": False,
        "reason": (
            "review-only full config; the delegated-approval launcher creates "
            "a hash-bound authorized derivative only after smoke passes"
        ),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    smoke = copy.deepcopy(config)
    smoke["train"].update(
        {
            "out_dir": "/tmp/prism_v5_joint_generator_count_smoke_sv6_20260805",
            "epochs": 3,
            "early_stopping": False,
            "early_stopping_min_epochs": 3,
            "early_stopping_patience": 3,
            "require_validation_for_early_stopping": False,
            "save_best": False,
            "save_every": 1,
            "max_train_steps_per_epoch": 2,
            "max_eval_steps": 1,
            "skip_final_test": True,
        }
    )
    smoke["learned_generator_count"].update(
        {"start_epoch": 2, "base_end_epoch": 3}
    )
    smoke["_launch_guard"] = {
        "launch_allowed": True,
        "reason": "production-topology smoke authorized by delegated user approval",
    }
    smoke["experiment_manifest"].update(
        {
            "smoke_only": True,
            "training_control": {
                "maximum_epochs": 3,
                "minimum_epochs_before_stop": 3,
                "early_stopping_patience": 3,
                "monitor_split": "val",
                "rank0_exact_validation": True,
                "save_every_epoch": True,
                "test_used": False,
            },
            "full_training_authorized": False,
            "requires_explicit_user_approval": True,
            "launch_performed": False,
        }
    )
    smoke_output.parent.mkdir(parents=True, exist_ok=True)
    smoke_output.write_text(
        json.dumps(smoke, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--smoke-output", type=Path, default=DEFAULT_SMOKE_OUTPUT
    )
    args = parser.parse_args()
    build(args.source, args.output, args.smoke_output)


if __name__ == "__main__":
    main()
