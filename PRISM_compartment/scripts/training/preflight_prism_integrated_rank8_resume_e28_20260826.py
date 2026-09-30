#!/usr/bin/env python3
"""Fail-closed audit for the epoch-28 CLS continuation configuration."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import torch

from kmlee_bam.training.runner_base import load_config


ALLOWED_DIFFERENCES = {
    "train.out_dir",
    "train.resume_checkpoint",
    "train.resume_optimizer",
    "train.resume_history",
    "prism_module_rescue.allow_state_readout_updates",
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
    parser.add_argument("--resume-config", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    source = json.loads(args.source_config.read_text(encoding="utf-8"))
    resume = json.loads(args.resume_config.read_text(encoding="utf-8"))
    parsed_resume = load_config(args.resume_config)
    source_flat = _flatten(source)
    resume_flat = _flatten(resume)
    differing = {
        key
        for key in set(source_flat) | set(resume_flat)
        if source_flat.get(key) != resume_flat.get(key)
    }
    unexpected = sorted(differing - ALLOWED_DIFFERENCES)
    if unexpected:
        raise RuntimeError(f"unexpected continuation config changes: {unexpected}")

    if resume["encoder"]["pooling"] != "cls":
        raise RuntimeError("resume pooling is not CLS")
    if parsed_resume.encoder.pooling != "cls":
        raise RuntimeError("parsed resume pooling is not CLS")
    if resume["prism_module_rescue"]["allow_state_readout_updates"] is not False:
        raise RuntimeError("CLS resume still enables AGP-only state-readout rescue")
    if resume["train"]["skip_final_test"] is not True:
        raise RuntimeError("official test split is not sealed")
    if resume["train"]["resume_optimizer"] is not True:
        raise RuntimeError("optimizer state would not be restored")
    if resume["train"]["resume_history"] is not True:
        raise RuntimeError("history would not be restored")
    if resume["train"]["out_dir"] == source["train"]["out_dir"]:
        raise RuntimeError("resume output must differ from the source run")

    checkpoint = Path(resume["train"]["resume_checkpoint"])
    expected_sha256 = resume["experiment_manifest"][
        "resume_source_checkpoint_sha256"
    ]
    observed_sha256 = _sha256(checkpoint)
    if observed_sha256 != expected_sha256:
        raise RuntimeError(
            f"resume checkpoint SHA-256 differs: {observed_sha256} != {expected_sha256}"
        )
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    if int(payload.get("epoch", 0)) != 28 or payload.get("step_in_epoch") is not None:
        raise RuntimeError("resume source is not the complete epoch-28 boundary")
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
        raise RuntimeError(f"resume checkpoint is incomplete: {missing}")
    if len(payload["history"]) != 28:
        raise RuntimeError("resume checkpoint history does not contain 28 epochs")

    report = {
        "schema_version": "kmlee_bam.prism_integrated_rank8_resume_preflight.v1",
        "status": "PASS",
        "source_epoch": 28,
        "next_epoch": 29,
        "checkpoint_sha256": observed_sha256,
        "optimizer_restored": True,
        "history_restored": True,
        "curriculum_state_present": True,
        "pooling": "cls",
        "agp_only_state_readout_rescue_enabled": False,
        "source_run_immutable": True,
        "official_test_dataset_opened": False,
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
