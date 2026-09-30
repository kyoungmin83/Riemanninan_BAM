#!/usr/bin/env python3
"""Verify every prelaunch receipt and seal the exact executable config."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--preflight", required=True)
    parser.add_argument("--smoke-checkpoint", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    preflight_path = Path(args.preflight).resolve()
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    smoke_checkpoint = Path(args.smoke_checkpoint).resolve()
    if preflight.get("status") != "PASS":
        raise RuntimeError("full-system preflight did not pass")
    if preflight.get("official_test_dataset_opened") is not False:
        raise RuntimeError("preflight opened the test dataset")
    if not smoke_checkpoint.is_file():
        raise FileNotFoundError("real DDP resume smoke checkpoint is missing")
    if config["_launch_guard"].get("launch_allowed") is not True:
        raise RuntimeError("config is not explicitly authorized")
    if config["train"].get("skip_final_test") is not True:
        raise RuntimeError("final config can open the official test split")
    if config["decoder"].get("pathology_rank") != 0:
        raise RuntimeError("learned pathology rank AUTO sentinel is missing")
    rank = config["learned_pathology_rank"]
    generator = config["learned_generator_count"]
    if (rank["warmup_fixed_rank"], rank["soft_start_epoch"], rank["freeze_epoch"]) != (
        8,
        25,
        37,
    ):
        raise RuntimeError("rank-first schedule changed")
    if (
        generator["start_epoch"],
        generator["shadow_end_epoch"],
        generator["soft_end_epoch"],
    ) != (37, 40, 44):
        raise RuntimeError("generator schedule changed")
    if generator.get("require_frozen_pathology_rank_before_search") is not True:
        raise RuntimeError("rank-to-generator ordering guard is disabled")
    if config["integrated_phase_curriculum"].get("phase2_pathology_scale") != 1.0:
        raise RuntimeError("pathology route becomes dormant")
    precision = config["precision_medicine"]
    if (
        precision.get("module_local_nonlinear_enabled") is not True
        or precision.get("module_local_nonlinear_variant")
        != "compartmental_threshold"
    ):
        raise RuntimeError("paper-inspired primary arm is not active")
    if config["pi_tech"].get("bank_invalid_fallback") != "zero":
        raise RuntimeError("technical-zero bank does not fail closed")
    critical = config["experiment_manifest"]["critical_model_inputs"]
    path_map = {
        "precision_context_npz_sha256": precision["context_npz_path"],
        "module_registry_json_sha256": config["module_tokenizer"][
            "registry_json_path"
        ],
        "activity_weight_npz_sha256": config["module_tokenizer"][
            "activity_weight_path"
        ],
        "module_local_reliability_sha256": precision[
            "module_local_reliability_npz_path"
        ],
        "module_local_compartment_graph_sha256": precision[
            "module_local_compartment_graph_npz_path"
        ],
        "module_local_output_cap_sha256": precision[
            "module_local_output_cap_npz_path"
        ],
        "module_rescue_stats_sha256": config["prism_module_rescue"]["stats_path"],
    }
    observed = {}
    for key, path in path_map.items():
        observed[key] = sha256_file(path)
        if observed[key] != critical[key]:
            raise RuntimeError(f"critical input changed after config build: {key}")
    receipt = {
        "schema_version": "kmlee_bam.prism_integrated_compartmental_pathrank_launch.v1",
        "status": "PASS",
        "config_path": str(config_path),
        "config_sha256": sha256_file(config_path),
        "preflight_path": str(preflight_path),
        "preflight_sha256": sha256_file(preflight_path),
        "ddp_resume_smoke_checkpoint": str(smoke_checkpoint),
        "ddp_resume_smoke_checkpoint_sha256": sha256_file(smoke_checkpoint),
        "critical_model_inputs": observed,
        "official_test_dataset_opened": False,
        "legacy_test_used_for_model_selection": False,
        "confirmatory_holdout_required_before_final_release": True,
        "explicit_user_launch_approval": config["_launch_guard"][
            "explicit_launch_approval"
        ],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
