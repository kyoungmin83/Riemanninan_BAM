#!/usr/bin/env python3
"""Create an audited continuation checkpoint with singleton-only gate protection.

This migration deliberately preserves model parameters, gate logits, optimizer,
scheduler, history, and curriculum state.  It changes only the persistent gate
protection buffer and controller state whose meaning depends on that protection
set.  The source checkpoint is never modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import torch


KEYS = {
    "log_alpha": "generator_count_gate.log_alpha",
    "protected": "generator_count_gate.protected_mask",
    "last_safe": "generator_count_gate.last_safe_log_alpha",
    "streak": "generator_count_gate.safe_audit_streak",
    "has_safe": "generator_count_gate.has_confirmed_safe_mask",
    "dual": "generator_count_objective.dual",
    "violation_ema": "generator_count_objective.violation_ema",
    "ema_initialized": "generator_count_objective.ema_initialized",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-checkpoint", type=Path, required=True)
    parser.add_argument("--joint-count-json", type=Path, required=True)
    parser.add_argument("--output-checkpoint", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--expected-epoch", type=int, required=True)
    parser.add_argument("--expected-old-protected", type=int, default=326)
    parser.add_argument("--expected-singletons", type=int, default=10)
    return parser.parse_args()


def scalar_value(value: torch.Tensor) -> Any:
    if value.numel() != 1:
        raise ValueError(f"expected scalar tensor, got shape={tuple(value.shape)}")
    return value.item()


def main() -> None:
    args = parse_args()
    source = args.source_checkpoint.expanduser().resolve()
    joint_path = args.joint_count_json.expanduser().resolve()
    output = args.output_checkpoint.expanduser().resolve()
    receipt_path = args.receipt.expanduser().resolve()

    if not source.is_file():
        raise FileNotFoundError(source)
    if not joint_path.is_file():
        raise FileNotFoundError(joint_path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output checkpoint: {output}")
    if receipt_path.exists():
        raise FileExistsError(f"refusing to overwrite receipt: {receipt_path}")

    joint = json.loads(joint_path.read_text(encoding="utf-8"))
    if int(joint["epoch"]) != int(args.expected_epoch):
        raise ValueError(
            f"joint-count epoch={joint['epoch']} != expected={args.expected_epoch}"
        )
    local_ids = joint.get("selected_local_generator_ids")
    names = joint.get("selected_registry_module_names")
    if not isinstance(local_ids, list) or not isinstance(names, list):
        raise ValueError("joint-count artifact lacks selected generator mapping")
    if len(local_ids) != len(names):
        raise ValueError("selected generator ids/names have different lengths")
    singleton_ids = sorted(
        int(local_id)
        for local_id, name in zip(local_ids, names, strict=True)
        if str(name).startswith("SINGLETON::")
    )
    if len(singleton_ids) != int(args.expected_singletons):
        raise ValueError(
            f"singleton count={len(singleton_ids)} != expected={args.expected_singletons}"
        )

    source_sha = sha256(source)
    payload = torch.load(source, map_location="cpu", weights_only=False)
    if int(payload.get("epoch", -1)) != int(args.expected_epoch):
        raise ValueError(
            f"checkpoint epoch={payload.get('epoch')} != expected={args.expected_epoch}"
        )
    if payload.get("step_in_epoch") is not None:
        raise ValueError("migration requires an epoch-boundary checkpoint")
    state = payload.get("system_state_dict")
    if not isinstance(state, dict):
        raise ValueError("checkpoint lacks system_state_dict")
    missing = [key for key in KEYS.values() if key not in state]
    if missing:
        raise KeyError(f"checkpoint lacks generator-controller keys: {missing}")

    protected = state[KEYS["protected"]]
    log_alpha = state[KEYS["log_alpha"]]
    if protected.dtype != torch.bool or protected.ndim != 1:
        raise ValueError("protected_mask must be a one-dimensional bool tensor")
    if log_alpha.shape != protected.shape:
        raise ValueError("log_alpha/protected_mask shapes differ")
    old_count = int(protected.sum().item())
    if old_count != int(args.expected_old_protected):
        raise ValueError(
            f"old protected count={old_count} != expected={args.expected_old_protected}"
        )
    if not all(bool(protected[index]) for index in singleton_ids):
        raise ValueError("one or more singleton generators were not protected in source")

    old_log_alpha = log_alpha.detach().clone()
    new_protected = torch.zeros_like(protected)
    new_protected[singleton_ids] = True
    state[KEYS["protected"]] = new_protected
    state[KEYS["last_safe"]] = log_alpha.detach().clone()
    state[KEYS["streak"]] = torch.zeros_like(state[KEYS["streak"]])
    state[KEYS["has_safe"]] = torch.zeros_like(state[KEYS["has_safe"]])
    state[KEYS["dual"]] = torch.zeros_like(state[KEYS["dual"]])
    state[KEYS["violation_ema"]] = torch.zeros_like(state[KEYS["violation_ema"]])
    state[KEYS["ema_initialized"]] = torch.zeros_like(
        state[KEYS["ema_initialized"]]
    )

    if not torch.equal(state[KEYS["log_alpha"]], old_log_alpha):
        raise AssertionError("gate log_alpha changed during migration")
    if int(state[KEYS["protected"]].sum().item()) != len(singleton_ids):
        raise AssertionError("new protected count is incorrect")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}")
    torch.save(payload, temporary)
    os.replace(temporary, output)
    output_sha = sha256(output)

    reloaded = torch.load(output, map_location="cpu", weights_only=False)
    reloaded_state = reloaded["system_state_dict"]
    checks = {
        "epoch_preserved": int(reloaded["epoch"]) == int(args.expected_epoch),
        "global_step_preserved": reloaded.get("global_step")
        == payload.get("global_step"),
        "epoch_boundary": reloaded.get("step_in_epoch") is None,
        "log_alpha_bitwise_preserved": torch.equal(
            reloaded_state[KEYS["log_alpha"]], old_log_alpha
        ),
        "singleton_only_protection": int(
            reloaded_state[KEYS["protected"]].sum().item()
        )
        == len(singleton_ids),
        "safe_audit_streak_reset": scalar_value(
            reloaded_state[KEYS["streak"]]
        )
        == 0,
        "confirmed_safe_mask_reset": not bool(
            scalar_value(reloaded_state[KEYS["has_safe"]])
        ),
        "constraint_dual_reset": bool(
            torch.count_nonzero(reloaded_state[KEYS["dual"]]) == 0
        ),
        "constraint_ema_reset": bool(
            torch.count_nonzero(reloaded_state[KEYS["violation_ema"]]) == 0
        ),
        "constraint_ema_initialization_reset": not bool(
            scalar_value(reloaded_state[KEYS["ema_initialized"]])
        ),
        "optimizer_state_present": "optimizer_state_dict" in reloaded,
        "scheduler_state_present": "scheduler_state_dict" in reloaded,
        "curriculum_state_present": "curriculum_state" in reloaded,
        "official_test_not_marked_used": not bool(
            reloaded.get("official_test_used", False)
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"post-migration verification failed: {checks}")

    receipt = {
        "schema_version": "kmlee_bam.generator_protection_migration.v1",
        "status": "PASS",
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_sha,
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": output_sha,
        "source_epoch": int(args.expected_epoch),
        "global_step": int(payload.get("global_step", -1)),
        "old_protected_count": old_count,
        "new_protected_count": len(singleton_ids),
        "new_protected_local_generator_ids": singleton_ids,
        "preserved_state": [
            "model_parameters",
            "generator_gate_log_alpha",
            "optimizer",
            "scheduler",
            "history",
            "curriculum",
            "PHU_and_adaptive_state",
        ],
        "reset_state": [
            "generator_count_gate.protected_mask",
            "generator_count_gate.last_safe_log_alpha",
            "generator_count_gate.safe_audit_streak",
            "generator_count_gate.has_confirmed_safe_mask",
            "generator_count_objective.dual",
            "generator_count_objective.violation_ema",
            "generator_count_objective.ema_initialized",
        ],
        "checks": checks,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
