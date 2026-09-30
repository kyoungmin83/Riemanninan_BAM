#!/usr/bin/env python3
"""Build a fail-closed one-update epoch-12 to epoch-13 BF16 resume smoke."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import torch


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--checkpoint-sha256", required=True)
    parser.add_argument("--fixed-source", required=True, type=Path)
    parser.add_argument("--fixed-source-sha256", required=True)
    parser.add_argument("--output-config", required=True, type=Path)
    parser.add_argument("--output-run", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    for output in (args.output_config, args.audit, args.output_run):
        if output.exists():
            raise FileExistsError(f"refusing to overwrite smoke output: {output}")

    observed_checkpoint_sha = _sha256(args.checkpoint)
    if observed_checkpoint_sha != args.checkpoint_sha256.lower():
        raise RuntimeError("epoch-12 checkpoint SHA-256 mismatch")
    observed_source_sha = _sha256(args.fixed_source)
    if observed_source_sha != args.fixed_source_sha256.lower():
        raise RuntimeError("BF16 projection-fix source SHA-256 mismatch")

    source = json.loads(args.source_config.read_text(encoding="utf-8"))
    if source["train"]["amp"] is not True:
        raise RuntimeError("production config does not enable AMP")
    if source["train"]["amp_dtype"] != "bfloat16":
        raise RuntimeError("production config does not use BF16")
    if source["integrated_phase_curriculum"]["phase1_end_epoch"] != 12:
        raise RuntimeError("Phase-I boundary is not epoch 12")
    if source["integrated_phase_curriculum"]["phase2_start_epoch"] != 13:
        raise RuntimeError("Phase-II boundary is not epoch 13")
    if source["train"]["skip_final_test"] is not True:
        raise RuntimeError("official test split is not sealed")

    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    required = {
        "system_state_dict",
        "criterion_state_dict",
        "optimizer_state_dict",
        "scheduler_state_dict",
        "amp_scaler_state_dict",
        "history",
        "v7a_state",
        "v8_state",
        "curriculum_state",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise RuntimeError(f"epoch-12 checkpoint is incomplete: {missing}")
    if int(payload.get("epoch", 0)) != 12:
        raise RuntimeError("resume source is not epoch 12")
    if payload.get("step_in_epoch") is not None:
        raise RuntimeError("resume source is not an epoch-boundary checkpoint")
    if len(payload["history"]) != 12:
        raise RuntimeError("resume source does not contain 12 history records")

    smoke = copy.deepcopy(source)
    smoke["train"].update(
        {
            "epochs": 13,
            "out_dir": str(args.output_run),
            "resume_checkpoint": str(args.checkpoint),
            "resume_optimizer": True,
            "resume_history": True,
            "allow_in_place_resume": False,
            "early_stopping": False,
            "save_best": False,
            "save_every": 1,
            "max_train_steps_per_epoch": 2,
            "max_eval_steps": 1,
            "log_every": 1,
            "progress_log_every": 1,
            "eval_log_every": 1,
            "reload_best_before_test": False,
            "skip_final_test": True,
        }
    )
    smoke["data"].update(
        {"num_workers": 0, "persistent_workers": False, "pin_memory": False}
    )
    ordinal = smoke.setdefault("v4", {}).setdefault("ordinal_balance", {})
    ordinal["max_cells_for_bin_stats"] = 64
    ordinal["bin_stats_log_every"] = 65
    smoke["experiment_manifest"] = {
        **smoke["experiment_manifest"],
        "run_name": (
            smoke["experiment_manifest"]["run_name"]
            + "_e12_e13_bf16_recovery_smoke"
        ),
        "smoke_only": True,
        "resume_source_checkpoint": str(args.checkpoint),
        "resume_source_checkpoint_sha256": observed_checkpoint_sha,
        "bf16_projection_fix_source_sha256": observed_source_sha,
        "official_test_dataset_opened_by_training": False,
    }
    smoke["_launch_guard"] = {
        **smoke.get("_launch_guard", {}),
        "launch_allowed": True,
        "smoke_only": True,
        "official_test_dataset_open_forbidden": True,
    }

    args.output_config.parent.mkdir(parents=True, exist_ok=True)
    args.output_config.write_text(
        json.dumps(smoke, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    audit = {
        "schema_version": "kmlee_bam.prism_compartmental_e12_e13_bf16_smoke.v1",
        "status": "READY",
        "source_epoch": 12,
        "next_epoch": 13,
        "train_microbatches": 2,
        "expected_optimizer_updates": 1,
        "world_size_required": 6,
        "amp_dtype": "bfloat16",
        "checkpoint_sha256": observed_checkpoint_sha,
        "fixed_source_sha256": observed_source_sha,
        "optimizer_state_present": True,
        "amp_state_present": True,
        "curriculum_state_present": True,
        "official_test_dataset_opened": False,
        "production_run_modified": False,
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
