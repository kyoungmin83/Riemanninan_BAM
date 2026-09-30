#!/usr/bin/env python3
"""Build an audited integrated continuation config with free generator count."""

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
    parser.add_argument("--host", required=True)
    parser.add_argument("--authorization-note", required=True)
    parser.add_argument("--audit", type=Path, required=True)
    return parser.parse_args()


def canonical_sha(value: Any) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(data).hexdigest()


def scientific_contract(config: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(config)
    train = result["train"]
    for key in (
        "out_dir",
        "resume_checkpoint",
        "resume_optimizer",
        "resume_history",
        "allow_in_place_resume",
    ):
        train.pop(key, None)
    result["learned_generator_count"].pop(
        "protect_unique_gene_coverage", None
    )
    result.pop("_launch_guard", None)
    result.pop("experiment_manifest", None)
    return result


def main() -> None:
    args = parse_args()
    source_path = args.source_config.expanduser().resolve()
    output_path = args.output_config.expanduser().resolve()
    audit_path = args.audit.expanduser().resolve()
    if output_path.exists() or audit_path.exists():
        raise FileExistsError("refusing to overwrite generated config or audit")
    source = json.loads(source_path.read_text(encoding="utf-8"))
    config = copy.deepcopy(source)

    if config["precision_medicine"]["personal_rank"] != 2:
        raise ValueError("continuation is not personal-rank2")
    if not config["learned_generator_count"]["enabled"]:
        raise ValueError("generator-count learning is not enabled")
    if not config["learned_generator_count"].get(
        "protect_unique_gene_coverage", False
    ):
        raise ValueError("source config did not contain the erroneous protection")
    if not config["train"].get("skip_final_test", False):
        raise ValueError("official-test seal is not active")

    train = config["train"]
    train["out_dir"] = args.output_dir
    train["resume_checkpoint"] = args.resume_checkpoint
    train["resume_optimizer"] = True
    train["resume_history"] = True
    # Older frozen runtimes predate this optional TrainConfig field.  Preserve
    # schema compatibility by changing it only when the source runtime knows it.
    if "allow_in_place_resume" in train:
        train["allow_in_place_resume"] = False

    gate = config["learned_generator_count"]
    gate["protect_unique_gene_coverage"] = False
    gate["protect_singletons"] = True
    gate["minimum_active_generators"] = None

    guard = config.setdefault("_launch_guard", {})
    guard.update(
        {
            "launch_allowed": True,
            "authorized_design_answer": "integrated",
            "official_test_dataset_open_forbidden": True,
            "resume_checkpoint_sha256_required": args.resume_checkpoint_sha256,
            "explicit_user_reapproval_required": False,
            "explicit_user_reapproval_received": args.authorization_note,
            "generator_protection_migration_required": True,
            "generator_protection_migration_receipt": args.migration_receipt,
            "expected_generator_protected_count": 10,
        }
    )

    manifest = config.setdefault("experiment_manifest", {})
    previous_run_name = manifest.get("run_name")
    manifest.update(
        {
            "run_name": args.run_name,
            "host": args.host,
            "training_form": "one_architecture_one_optimizer_one_checkpoint_lineage",
            "continuation_of": previous_run_name,
            "resume_source_checkpoint": args.resume_checkpoint,
            "resume_source_checkpoint_sha256": args.resume_checkpoint_sha256,
            "resume_source_epoch": args.source_epoch,
            "resume_next_epoch": args.source_epoch + 1,
            "generator_cardinality_correction": (
                "remove unrequested blanket unique-gene-coverage protection; "
                "retain only 10 singleton generators"
            ),
            "generator_count_target": None,
            "generator_minimum_active": None,
            "generator_protected_singletons": 10,
            "official_test_donors": "sealed_not_opened",
        }
    )

    source_contract = scientific_contract(source)
    output_contract = scientific_contract(config)
    if source_contract != output_contract:
        raise AssertionError("a non-generator scientific setting changed")

    checks = {
        "personal_rank2_preserved": config["precision_medicine"]["personal_rank"]
        == 2,
        "integrated_curriculum_preserved": config.get(
            "integrated_phase_curriculum"
        )
        == source.get("integrated_phase_curriculum"),
        "unique_coverage_protection_disabled": not gate[
            "protect_unique_gene_coverage"
        ],
        "singleton_protection_retained": bool(gate["protect_singletons"]),
        "no_generator_count_floor": gate["minimum_active_generators"] is None,
        "optimizer_resume_enabled": bool(train["resume_optimizer"]),
        "history_resume_enabled": bool(train["resume_history"]),
        "new_output_directory_required": not bool(
            train.get("allow_in_place_resume", False)
        ),
        "official_test_sealed": bool(train["skip_final_test"]),
        "scientific_contract_identical_except_protection": source_contract
        == output_contract,
    }
    if not all(checks.values()):
        raise AssertionError(f"config verification failed: {checks}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    audit = {
        "schema_version": "kmlee_bam.generator_unprotected_continuation_config.v1",
        "status": "PASS",
        "source_config": str(source_path),
        "output_config": str(output_path),
        "source_epoch": args.source_epoch,
        "next_epoch": args.source_epoch + 1,
        "authorization_note": args.authorization_note,
        "source_scientific_contract_sha256": canonical_sha(source_contract),
        "output_scientific_contract_sha256": canonical_sha(output_contract),
        "checks": checks,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
