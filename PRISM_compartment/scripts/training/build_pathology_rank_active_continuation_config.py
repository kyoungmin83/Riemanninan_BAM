#!/usr/bin/env python3
"""Build an audited SV7 continuation with active pathology-rank learning."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-config", type=Path, required=True)
    parser.add_argument("--output-config", type=Path, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--resume-checkpoint", required=True)
    parser.add_argument("--resume-checkpoint-sha256", required=True)
    parser.add_argument("--migration-receipt", required=True)
    parser.add_argument("--source-epoch", type=int, required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--phase2-pathology-scale", type=float, default=0.1)
    parser.add_argument("--freeze-epoch", type=int, default=64)
    parser.add_argument("--audit", type=Path, required=True)
    return parser.parse_args()


def contract_without_authorized_changes(config: dict[str, Any]) -> dict[str, Any]:
    value = copy.deepcopy(config)
    train = value["train"]
    for key in (
        "out_dir",
        "resume_checkpoint",
        "resume_optimizer",
        "resume_history",
        "allow_in_place_resume",
    ):
        train.pop(key, None)
    value["learned_pathology_rank"].pop("freeze_epoch", None)
    value["integrated_phase_curriculum"].pop("phase2_pathology_scale", None)
    value.pop("_launch_guard", None)
    value.pop("experiment_manifest", None)
    return value


def canonical_sha(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    args = parse_args()
    source_path = args.source_config.expanduser().resolve()
    output_path = args.output_config.expanduser().resolve()
    audit_path = args.audit.expanduser().resolve()
    if output_path.exists() or audit_path.exists():
        raise FileExistsError("refusing to overwrite generated rank config")
    config = json.loads(source_path.read_text(encoding="utf-8"))
    source = copy.deepcopy(config)
    scale = float(args.phase2_pathology_scale)
    if not 0.0 < scale <= 1.0:
        raise ValueError("phase2 pathology scale must lie in (0, 1]")
    if int(args.freeze_epoch) <= int(config["train"]["epochs"]):
        raise ValueError("freeze_epoch must remain after the nominal training run")
    if not config["learned_pathology_rank"]["enabled"]:
        raise ValueError("learned pathology rank is disabled")
    if int(config["decoder"]["pathology_rank"]) != 0:
        raise ValueError("learned-rank AUTO sentinel changed")
    gate = config["learned_generator_count"]
    if gate.get("protect_unique_gene_coverage", True):
        raise ValueError("generator blanket protection reappeared")

    train = config["train"]
    train["out_dir"] = args.output_dir
    train["resume_checkpoint"] = args.resume_checkpoint
    train["resume_optimizer"] = True
    train["resume_history"] = True
    if "allow_in_place_resume" in train:
        train["allow_in_place_resume"] = False
    config["learned_pathology_rank"]["freeze_epoch"] = int(args.freeze_epoch)
    config["integrated_phase_curriculum"]["phase2_pathology_scale"] = scale

    guard = config.setdefault("_launch_guard", {})
    guard.update(
        {
            "launch_allowed": True,
            "authorized_design_answer": "integrated",
            "official_test_dataset_open_forbidden": True,
            "resume_checkpoint_sha256_required": args.resume_checkpoint_sha256,
            "pathology_rank_learning_required": True,
            "pathology_rank_migration_receipt": args.migration_receipt,
            "pathology_rank_expected_finalized_at_resume": False,
            "pathology_rank_freeze_epoch": int(args.freeze_epoch),
            "phase2_pathology_scale": scale,
            "explicit_user_reapproval_required": False,
            "explicit_user_reapproval_received": (
                "2026-08-26 KST: user explicitly ordered continued pathology-rank "
                "cardinality learning and rank-count training logs"
            ),
        }
    )
    manifest = config.setdefault("experiment_manifest", {})
    previous = manifest.get("run_name")
    manifest.update(
        {
            "run_name": args.run_name,
            "continuation_of": previous,
            "resume_source_checkpoint": args.resume_checkpoint,
            "resume_source_checkpoint_sha256": args.resume_checkpoint_sha256,
            "resume_source_epoch": int(args.source_epoch),
            "resume_next_epoch": int(args.source_epoch) + 1,
            "pathology_rank_learning": "active_hard_straight_through",
            "pathology_rank_target": None,
            "pathology_rank_minimum": 0,
            "pathology_rank_freeze_epoch": int(args.freeze_epoch),
            "phase2_pathology_scale": scale,
            "pathology_rank_epoch_logging": "start_and_end_main_training_log",
            "official_test_donors": "sealed_not_opened",
        }
    )

    source_contract = contract_without_authorized_changes(source)
    output_contract = contract_without_authorized_changes(config)
    if source_contract != output_contract:
        raise AssertionError("unauthorized scientific configuration change")
    checks = {
        "integrated_curriculum_preserved_except_pathology_floor": True,
        "personal_rank2_preserved": config["precision_medicine"]["personal_rank"]
        == 2,
        "learned_pathology_rank_enabled": bool(
            config["learned_pathology_rank"]["enabled"]
        ),
        "rank_freeze_after_nominal_training": int(
            config["learned_pathology_rank"]["freeze_epoch"]
        )
        > int(train["epochs"]),
        "phase2_pathology_gradient_active": float(
            config["integrated_phase_curriculum"]["phase2_pathology_scale"]
        )
        > 0.0,
        "rank_has_no_target": manifest["pathology_rank_target"] is None,
        "rank_has_no_minimum": int(manifest["pathology_rank_minimum"]) == 0,
        "generator_blanket_protection_disabled": not bool(
            gate["protect_unique_gene_coverage"]
        ),
        "official_test_sealed": bool(train["skip_final_test"]),
        "optimizer_resume_enabled": bool(train["resume_optimizer"]),
        "history_resume_enabled": bool(train["resume_history"]),
        "authorized_contract_only": source_contract == output_contract,
    }
    if not all(checks.values()):
        raise AssertionError(f"rank continuation config checks failed: {checks}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    audit = {
        "schema_version": "kmlee_bam.pathology_rank_active_continuation.v1",
        "status": "PASS",
        "source_config": str(source_path),
        "output_config": str(output_path),
        "source_epoch": int(args.source_epoch),
        "next_epoch": int(args.source_epoch) + 1,
        "phase2_pathology_scale": scale,
        "rank_freeze_epoch": int(args.freeze_epoch),
        "source_contract_sha256": canonical_sha(source_contract),
        "output_contract_sha256": canonical_sha(output_contract),
        "checks": checks,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
