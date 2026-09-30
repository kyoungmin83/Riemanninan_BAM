#!/usr/bin/env python3
"""Reactivate learned pathology rank in an audited epoch-boundary checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import torch


RANK_KEYS = {
    "log_alpha": "decoder.pathology_rank_gate.log_alpha",
    "epoch": "decoder.pathology_rank_gate.epoch_state",
    "finalized": "decoder.pathology_rank_gate.finalized_state",
    "frozen_mask": "decoder.pathology_rank_gate.frozen_mask",
    "pathology_scale": "decoder.pathology_curriculum_scale",
}
GENERATOR_PROTECTED_KEY = "generator_count_gate.protected_mask"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-checkpoint", type=Path, required=True)
    parser.add_argument("--output-checkpoint", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--expected-epoch", type=int, required=True)
    parser.add_argument("--phase2-pathology-scale", type=float, required=True)
    parser.add_argument("--expected-generator-protected", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.source_checkpoint.expanduser().resolve()
    output = args.output_checkpoint.expanduser().resolve()
    receipt_path = args.receipt.expanduser().resolve()
    scale = float(args.phase2_pathology_scale)
    if not 0.0 < scale <= 1.0:
        raise ValueError("phase2 pathology scale must lie in (0, 1]")
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite rank migration output")

    source_sha = sha256(source)
    payload = torch.load(source, map_location="cpu", weights_only=False)
    if int(payload.get("epoch", -1)) != int(args.expected_epoch):
        raise ValueError("source epoch mismatch")
    if payload.get("step_in_epoch") is not None:
        raise ValueError("rank migration requires an epoch-boundary checkpoint")
    state = payload["system_state_dict"]
    missing = [key for key in RANK_KEYS.values() if key not in state]
    if missing:
        raise KeyError(f"checkpoint lacks rank state: {missing}")
    if GENERATOR_PROTECTED_KEY not in state:
        raise KeyError("checkpoint lacks generator protected mask")
    protected_count = int(state[GENERATOR_PROTECTED_KEY].sum().item())
    if protected_count != int(args.expected_generator_protected):
        raise ValueError(
            f"generator protected count={protected_count} != "
            f"expected={args.expected_generator_protected}"
        )

    log_alpha = state[RANK_KEYS["log_alpha"]].detach().clone()
    old_expected = float(torch.sigmoid(log_alpha).sum().item())
    old_hard = int((torch.sigmoid(log_alpha) >= 0.5).sum().item())
    old_finalized = bool(state[RANK_KEYS["finalized"]].item())
    if not old_finalized:
        raise ValueError("source rank gate is not finalized")

    state[RANK_KEYS["finalized"]] = torch.zeros_like(
        state[RANK_KEYS["finalized"]]
    )
    # Keep the learned logits and old hard mask bitwise intact.  The mask is
    # ignored while active and remains useful provenance until finalization.
    state[RANK_KEYS["pathology_scale"]] = torch.full_like(
        state[RANK_KEYS["pathology_scale"]], scale
    )
    curriculum = payload.get("curriculum_state")
    if not isinstance(curriculum, dict):
        raise ValueError("checkpoint lacks curriculum_state")
    controller = curriculum.get("integrated_phase_controller")
    if not isinstance(controller, dict):
        raise ValueError("checkpoint lacks integrated phase controller state")
    if int(controller.get("epoch", -1)) != int(args.expected_epoch):
        raise ValueError("curriculum epoch mismatch")
    controller["pathology_scale"] = scale

    if not torch.equal(state[RANK_KEYS["log_alpha"]], log_alpha):
        raise AssertionError("rank log_alpha changed during migration")

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
        "rank_log_alpha_bitwise_preserved": torch.equal(
            reloaded_state[RANK_KEYS["log_alpha"]], log_alpha
        ),
        "rank_gate_reactivated": not bool(
            reloaded_state[RANK_KEYS["finalized"]].item()
        ),
        "phase2_pathology_scale_installed": abs(
            float(reloaded_state[RANK_KEYS["pathology_scale"]].item()) - scale
        )
        <= 1.0e-7,
        "curriculum_scale_migrated": abs(
            float(
                reloaded["curriculum_state"]["integrated_phase_controller"][
                    "pathology_scale"
                ]
            )
            - scale
        )
        <= 1.0e-12,
        "generator_singleton_only_protection_preserved": int(
            reloaded_state[GENERATOR_PROTECTED_KEY].sum().item()
        )
        == int(args.expected_generator_protected),
        "optimizer_state_present": "optimizer_state_dict" in reloaded,
        "scheduler_state_present": "scheduler_state_dict" in reloaded,
        "history_present": "history" in reloaded,
        "v7a_state_present": "v7a_state" in reloaded,
        "v8_state_present": "v8_state" in reloaded,
    }
    if not all(checks.values()):
        raise AssertionError(f"rank migration verification failed: {checks}")

    receipt = {
        "schema_version": "kmlee_bam.pathology_rank_learning_migration.v1",
        "status": "PASS",
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_sha,
        "output_checkpoint": str(output),
        "output_checkpoint_sha256": output_sha,
        "source_epoch": int(args.expected_epoch),
        "global_step": int(payload.get("global_step", -1)),
        "rank_capacity": int(log_alpha.numel()),
        "rank_expected_at_migration": old_expected,
        "rank_hard_at_migration": old_hard,
        "rank_finalized_before": old_finalized,
        "rank_finalized_after": False,
        "phase2_pathology_scale": scale,
        "generator_protected_count": protected_count,
        "checks": checks,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
