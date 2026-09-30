#!/usr/bin/env python3
"""Build the same-lineage SV7 config with the invalid CLS/AGP rescue flag fixed."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "configs/final/train_config_prism_integrated_pathranklearn_s42_sv7_pre_recovery_20260826.json"
)
OUTPUT = (
    ROOT
    / "configs/final/train_config_prism_integrated_pathranklearn_s42_sv7_cls_recovery_20260826.json"
)
SOURCE_RUN = (
    "/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
    "kmlee_bam_prism_integrated_pathranklearn_s42_sv7_20260824"
)
SOURCE_CHECKPOINT = f"{SOURCE_RUN}/checkpoint_epoch_016.pt"
SOURCE_CHECKPOINT_SHA256 = (
    "c92bb3bbe0258bfba3d30626248a0b83663e0142746e9b0ad24f9c3603004aef"
)
RECOVERY_RUN_NAME = (
    "kmlee_bam_prism_integrated_pathranklearn_s42_sv7_resume_e16_20260826"
)
RECOVERY_RUN = (
    "/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/"
    f"{RECOVERY_RUN_NAME}"
)


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    config = copy.deepcopy(source)

    if config["encoder"]["pooling"] != "cls":
        raise ValueError("the SV7 integrated model must retain CLS pooling")
    if config["prism_module_rescue"]["allow_state_readout_updates"] is not True:
        raise ValueError("the source does not contain the diagnosed CLS/AGP conflict")
    if config["learned_pathology_rank"]["enabled"] is not True:
        raise ValueError("learned pathology rank must remain enabled")
    if int(config["decoder"]["pathology_rank"]) != 0:
        raise ValueError("learned-rank AUTO sentinel changed")
    if config["train"]["skip_final_test"] is not True:
        raise ValueError("official test donors must remain sealed")

    config["prism_module_rescue"]["allow_state_readout_updates"] = False
    config["train"].update(
        {
            "out_dir": RECOVERY_RUN,
            "resume_checkpoint": SOURCE_CHECKPOINT,
            "resume_optimizer": True,
            "resume_history": True,
            "skip_final_test": True,
        }
    )
    config["experiment_manifest"].update(
        {
            "run_name": RECOVERY_RUN_NAME,
            "continuation_of": source["experiment_manifest"]["run_name"],
            "same_lineage_corrective_resume": True,
            "corrective_resume_source_epoch": 16,
            "corrective_resume_next_epoch": 17,
            "corrective_resume_source_checkpoint": SOURCE_CHECKPOINT,
            "corrective_resume_source_checkpoint_sha256": SOURCE_CHECKPOINT_SHA256,
            "corrective_change": (
                "disable AGP-only state-readout rescue updates in the unchanged "
                "CLS architecture"
            ),
            "corrective_change_scope": (
                "no model, optimizer, scheduler, curriculum, learned-rank, data, "
                "or test-split change"
            ),
            "explicit_user_reapproval_received": (
                "2026-08-26 KST: user confirmed continued integrated PRISM"
            ),
        }
    )
    config["_launch_guard"].update(
        {
            "launch_allowed": True,
            "explicit_user_reapproval_required": False,
            "cls_agp_conflict_corrected": True,
            "resume_checkpoint_sha256_required": SOURCE_CHECKPOINT_SHA256,
        }
    )

    OUTPUT.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
