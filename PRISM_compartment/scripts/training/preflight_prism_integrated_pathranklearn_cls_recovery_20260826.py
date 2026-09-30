#!/usr/bin/env python3
"""Fail-closed audit for the SV7 epoch-16 CLS corrective continuation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import torch

from kmlee_bam.training.runner_base import load_config


ALLOWED_DIFFERENCES = {
    "prism_module_rescue.allow_state_readout_updates",
    "train.out_dir",
    "train.resume_checkpoint",
}


def _flatten(value: object, prefix: str = "") -> dict[str, object]:
    if not isinstance(value, dict):
        return {prefix: value}
    result: dict[str, object] = {}
    for key, child in value.items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if path.startswith("experiment_manifest") or path.startswith("_launch_guard"):
            continue
        result.update(_flatten(child, path))
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-config", required=True, type=Path)
    parser.add_argument("--recovery-config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    source = json.loads(args.source_config.read_text(encoding="utf-8"))
    recovery = json.loads(args.recovery_config.read_text(encoding="utf-8"))
    parsed = load_config(args.recovery_config)
    source_flat = _flatten(source)
    recovery_flat = _flatten(recovery)
    differing = {
        key
        for key in set(source_flat) | set(recovery_flat)
        if source_flat.get(key) != recovery_flat.get(key)
    }
    unexpected = sorted(differing - ALLOWED_DIFFERENCES)
    if unexpected:
        raise RuntimeError(f"unexpected corrective config changes: {unexpected}")
    if parsed.encoder.pooling != "cls":
        raise RuntimeError("corrective continuation changed CLS pooling")
    if recovery["prism_module_rescue"]["allow_state_readout_updates"] is not False:
        raise RuntimeError("AGP-only state-readout rescue remains enabled")
    if recovery["learned_pathology_rank"]["enabled"] is not True:
        raise RuntimeError("learned pathology rank was disabled")
    if int(recovery["decoder"]["pathology_rank"]) != 0:
        raise RuntimeError("learned-rank AUTO sentinel changed")
    if recovery["train"]["skip_final_test"] is not True:
        raise RuntimeError("official test donors are not sealed")
    if recovery["train"]["out_dir"] == source["train"]["out_dir"]:
        raise RuntimeError("corrective continuation must preserve the source output")
    if Path(recovery["train"]["resume_checkpoint"]).resolve() != args.checkpoint.resolve():
        raise RuntimeError("configured resume checkpoint differs from the audited source")

    checkpoint_sha256 = _sha256(args.checkpoint)
    payload = torch.load(
        args.checkpoint,
        map_location="cpu",
        weights_only=False,
        mmap=True,
    )
    if int(payload.get("epoch", 0)) != 16 or payload.get("step_in_epoch") is not None:
        raise RuntimeError("source checkpoint is not the complete epoch-16 boundary")
    required = {
        "system_state_dict",
        "criterion_state_dict",
        "optimizer_state_dict",
        "scheduler_state_dict",
        "history",
        "v7a_state",
        "v8_state",
        "curriculum_state",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise RuntimeError(f"epoch-16 checkpoint is incomplete: {missing}")
    if len(payload["history"]) != 16:
        raise RuntimeError("epoch-16 checkpoint history length is not 16")
    state = payload["system_state_dict"]
    rank_keys = {
        "decoder.pathology_rank_gate.log_alpha",
        "decoder.pathology_rank_gate.epoch_state",
        "decoder.pathology_rank_gate.finalized_state",
        "decoder.pathology_rank_gate.frozen_mask",
    }
    missing_rank = sorted(rank_keys - set(state))
    if missing_rank:
        raise RuntimeError(f"learned-rank state is incomplete: {missing_rank}")
    if bool(state["decoder.pathology_rank_gate.finalized_state"].item()):
        raise RuntimeError("learned pathology rank finalized before epoch 17")

    report = {
        "schema_version": "kmlee_bam.prism_pathrank_cls_recovery_preflight.v1",
        "status": "PASS",
        "source_epoch": 16,
        "next_epoch": 17,
        "checkpoint": str(args.checkpoint.resolve()),
        "checkpoint_sha256": checkpoint_sha256,
        "optimizer_restored": True,
        "scheduler_restored": True,
        "history_restored": True,
        "curriculum_state_present": True,
        "learned_rank_state_present": True,
        "learned_rank_finalized_before_resume": False,
        "pooling": "cls",
        "agp_only_state_readout_rescue_enabled": False,
        "official_test_dataset_opened": False,
        "source_output_immutable": True,
        "allowed_config_differences": sorted(differing),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
