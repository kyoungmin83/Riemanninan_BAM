#!/usr/bin/env python3
"""Build the fail-closed SV6 continuation config after the CLS/AGP guard stop."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "configs/final/train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
)
OUTPUT = (
    ROOT
    / "configs/final/train_config_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826.json"
)

SOURCE_RUN = (
    "/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
    "kmlee_bam_prism_integrated_rank8_primary_s42_sv6_retry1_20260824"
)
SOURCE_CHECKPOINT = f"{SOURCE_RUN}/checkpoint_epoch_028.pt"
SOURCE_CHECKPOINT_SHA256 = (
    "8ebb2f0644da8bdd985fa11f66dba6719ad7fd20ee6141f54dac71ae064eb108"
)
RUN_NAME = (
    "kmlee_bam_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826"
)
RUNTIME_ROOT = (
    "/home/kmlee/project_sv6/"
    "kmlee_bam_integrated_rank8_primary_resume_e28_20260826"
)
RUN_DIR = (
    "/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
    f"{RUN_NAME}"
)


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    config = copy.deepcopy(source)

    if config["encoder"]["pooling"] != "cls":
        raise ValueError("the primary continuation must retain CLS pooling")
    if config["prism_module_rescue"]["allow_state_readout_updates"] is not False:
        raise ValueError("CLS continuation must disable AGP-only state-readout rescue")
    if config["train"]["skip_final_test"] is not True:
        raise ValueError("official test donors must remain sealed")

    config["train"].update(
        {
            "out_dir": RUN_DIR,
            "resume_checkpoint": SOURCE_CHECKPOINT,
            "resume_optimizer": True,
            "resume_history": True,
            "early_stopping": False,
            "reload_best_before_test": False,
            "skip_final_test": True,
        }
    )
    config["experiment_manifest"].update(
        {
            "run_name": RUN_NAME,
            "runtime_root": RUNTIME_ROOT,
            "continuation_of": config["experiment_manifest"]["run_name"],
            "resume_source_checkpoint": SOURCE_CHECKPOINT,
            "resume_source_checkpoint_sha256": SOURCE_CHECKPOINT_SHA256,
            "resume_source_epoch": 28,
            "resume_next_epoch": 29,
            "retry_reason": (
                "epoch29 fail-closed guard rejected AGP-only state-readout "
                "rescue updates in the primary CLS architecture"
            ),
            "corrective_change": (
                "prism_module_rescue.allow_state_readout_updates=false; "
                "all learned model and optimizer states restored strictly"
            ),
            "official_test_donors": "sealed_not_opened",
        }
    )
    config["_launch_guard"].update(
        {
            "launch_allowed": True,
            "resume_preflight_required": True,
            "resume_checkpoint_sha256_required": SOURCE_CHECKPOINT_SHA256,
            "resume_failure_fix": "disable_agp_only_state_readout_rescue_for_cls",
            "explicit_user_reapproval_required": False,
            "explicit_user_reapproval_received": (
                "2026-08-26 KST: user confirmed continuation as integrated PRISM"
            ),
        }
    )

    OUTPUT.write_text(
        json.dumps(config, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
