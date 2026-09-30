"""Production adapter for one sealed fixed-generator training job.

The adapter turns an exact-evaluation training job into a normal
``run_current`` configuration.  It does not implement a second trainer: the
existing runner remains responsible for data loading, optimisation, DDP and
checkpointing.  This module only validates immutable inputs, installs the
sealed donor/mask contract, and audits the produced checkpoint.

``--validate-only`` is deliberately non-executing and is the mode used before
operator approval.  Omitting it launches the ordinary trainer and can be
expensive.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, Sequence

import torch

from kmlee_bam.training.generator_exact_confirmation import (
    canonical_json_sha256,
    generator_mask_sha256,
)
from kmlee_bam.training.generator_exact_evaluation import JOB_SCHEMA
from kmlee_bam.training.generator_exact_executor import TRAINING_OUTPUT_SCHEMA


BACKEND_PREFLIGHT_SCHEMA = "kmlee_bam.generator_exact_training_preflight.v1"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json_object(path: str | Path, *, label: str) -> dict[str, Any]:
    source = Path(path).resolve()
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"failed to read {label}: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return value


def _positive_int(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return int(value)


def validate_fixed_mask_training_job(job: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete immutable job before config generation."""

    result = dict(job)
    if result.get("schema_version") != JOB_SCHEMA:
        raise ValueError("fixed-mask training job schema mismatch")
    supplied = result.get("job_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "job_sha256"}
    )
    if supplied != expected:
        raise ValueError("fixed-mask training job SHA256 mismatch")
    partition = result.get("donor_partition")
    expected_kind = {
        "ranking": "fixed_mask_adaptation_w44",
        "confirmation": "paired_fixed_mask_retrain_w_plus_r54",
    }.get(partition)
    if result.get("job_kind") != expected_kind:
        raise ValueError("fixed-mask training job kind/partition mismatch")
    candidate_count = _positive_int(
        result.get("candidate_count"), label="candidate_count"
    )
    selected = result.get("selected_local_generator_ids")
    if not isinstance(selected, list) or not selected:
        raise ValueError("selected_local_generator_ids must be non-empty")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in selected):
        raise ValueError("selected_local_generator_ids must contain integers")
    selected_ids = tuple(int(value) for value in selected)
    if len(set(selected_ids)) != len(selected_ids) or any(
        value < 0 or value >= candidate_count for value in selected_ids
    ):
        raise ValueError("selected generator ids are duplicate or out of range")
    if result.get("generator_count") != len(selected_ids):
        raise ValueError("generator_count differs from selected-id count")
    if result.get("mask_sha256") != generator_mask_sha256(
        selected_ids, candidate_count=candidate_count
    ):
        raise ValueError("fixed-mask job binary mask SHA256 mismatch")

    donors = result.get("training_donor_names")
    expected_count = 44 if partition == "ranking" else 54
    if not isinstance(donors, list) or len(donors) != expected_count:
        raise ValueError(f"training job requires exactly {expected_count} donors")
    if len({str(value) for value in donors}) != expected_count:
        raise ValueError("training donor names must be unique")
    if result.get("training_donor_names_sha256") != canonical_json_sha256(donors):
        raise ValueError("training donor order SHA256 mismatch")

    contract = result.get("training_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("fixed-mask training job lacks training_contract")
    for field in (
        "fresh_optimizer",
        "fixed_mask_for_entire_run",
        "same_seed_schedule_across_k",
        "initial_weights_identical_source_checkpoint",
    ):
        if contract.get(field) is not True:
            raise ValueError(f"training contract requires {field}=true")
    for field in (
        "mask_trainable",
        "acute_mask_ablation",
        "learned_generator_count_enabled",
        "legacy_generator_budget_search_enabled",
        "official_validation_evaluation",
        "official_test_evaluation",
    ):
        if contract.get(field) is not False:
            raise ValueError(f"training contract requires {field}=false")
    _positive_int(contract.get("epochs"), label="training_contract.epochs")
    seed = contract.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("training_contract.seed must be a non-negative integer")

    paths = result.get("source_paths")
    hashes = result.get("source_sha256")
    if not isinstance(paths, Mapping) or not isinstance(hashes, Mapping):
        raise ValueError("training job requires source_paths/source_sha256")
    required_paths = {
        "config",
        "checkpoint",
        "registry",
        "scaler",
        "module_stats",
        "donor_context",
        "age_normalization",
        "donor_split_manifest",
    }
    missing = sorted(required_paths - set(paths))
    if missing:
        raise ValueError(f"training job lacks source paths: {missing}")
    for name, raw_path in paths.items():
        source = Path(str(raw_path)).resolve()
        if not source.is_file():
            raise ValueError(f"training source {name} does not exist: {source}")
        if hashes.get(name) != sha256_file(source):
            raise ValueError(f"training source {name} SHA256 mismatch")
    return result


def build_run_current_config(
    job: Mapping[str, Any],
    *,
    job_path: str | Path,
    work_dir: str | Path,
) -> dict[str, Any]:
    """Create the exact normal-runner config without launching it."""

    validated = validate_fixed_mask_training_job(job)
    source_paths = validated["source_paths"]
    raw = read_json_object(source_paths["config"], label="source training config")
    config = deepcopy(raw)
    contract = validated["training_contract"]
    output = Path(work_dir).resolve()

    config["_doc_pointer"] = (
        "generated fixed-mask exact-evaluation retrain; do not edit or reuse "
        "outside the sealed job"
    )
    config["warm_start"] = {
        "init_weights_path": str(Path(source_paths["checkpoint"]).resolve()),
        "expected_sha256": validated["source_sha256"]["checkpoint"],
        "require_precision_head": True,
    }
    config["generator_exact_training"] = {
        "enabled": True,
        "job_path": str(Path(job_path).resolve()),
        "expected_job_sha256": validated["job_sha256"],
        "donor_split_manifest_path": str(
            Path(source_paths["donor_split_manifest"]).resolve()
        ),
        "donor_split_manifest_sha256": validated["source_sha256"][
            "donor_split_manifest"
        ],
        "fresh_optimizer": True,
        "resume_forbidden": True,
    }
    config.setdefault("learned_generator_count", {})["enabled"] = False
    config.setdefault("generator_budget_search", {})["enabled"] = False
    config.setdefault("module_tokenizer", {})["registry_json_path"] = str(
        Path(source_paths["registry"]).resolve()
    )
    config.setdefault("data", {})["reference_scaler_path"] = str(
        Path(source_paths["scaler"]).resolve()
    )
    config.setdefault("precision_medicine", {})["context_npz_path"] = str(
        Path(source_paths["donor_context"]).resolve()
    )
    rescue = config.setdefault("prism_module_rescue", {})
    if bool(rescue.get("enabled", False)):
        rescue["stats_path"] = str(Path(source_paths["module_stats"]).resolve())
        rescue["expected_stats_sha256"] = validated["source_sha256"][
            "module_stats"
        ]
    decoder_count = config.setdefault("decoder", {}).get("n_generators")
    if decoder_count != validated["candidate_count"]:
        raise ValueError(
            "source config decoder.n_generators differs from candidate_count"
        )
    train = config.setdefault("train", {})
    train.update(
        {
            "out_dir": str(output),
            "epochs": int(contract["epochs"]),
            "seed": int(contract["seed"]),
            "early_stopping": False,
            "require_validation_for_early_stopping": False,
            "save_best": False,
            "save_split_latents": False,
            "reload_best_before_test": False,
            "skip_final_test": True,
            "use_rank0_eval": False,
            "save_every": 0,
            "save_every_steps": 0,
        }
    )
    return config


def build_training_preflight(
    job: Mapping[str, Any],
    *,
    job_path: str | Path,
    work_dir: str | Path,
) -> dict[str, Any]:
    config = build_run_current_config(job, job_path=job_path, work_dir=work_dir)
    validated = validate_fixed_mask_training_job(job)
    result: dict[str, Any] = {
        "schema_version": BACKEND_PREFLIGHT_SCHEMA,
        "status": "validated_not_executed",
        "job_sha256": validated["job_sha256"],
        "donor_partition": validated["donor_partition"],
        "training_donor_count": len(validated["training_donor_names"]),
        "generator_count": validated["generator_count"],
        "mask_sha256": validated["mask_sha256"],
        "epochs": validated["training_contract"]["epochs"],
        "seed": validated["training_contract"]["seed"],
        "run_current_config_sha256": canonical_json_sha256(config),
        "training_started": False,
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["preflight_sha256"] = canonical_json_sha256(result)
    return result


def audit_training_checkpoint(
    checkpoint_path: str | Path,
    job: Mapping[str, Any],
) -> dict[str, Any]:
    validated = validate_fixed_mask_training_job(job)
    checkpoint = Path(checkpoint_path).resolve()
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    if not isinstance(payload, Mapping):
        raise ValueError("fixed-mask checkpoint must contain a mapping")
    provenance = payload.get("fixed_generator_result_provenance")
    integrity = payload.get("fixed_generator_integrity")
    if not isinstance(provenance, Mapping) or not isinstance(integrity, Mapping):
        raise ValueError("checkpoint lacks fixed-mask provenance/integrity audit")
    required = {
        "exact_training_job_sha256": validated["job_sha256"],
        "candidate_count": validated["candidate_count"],
        "hard_active": validated["generator_count"],
        "selected_local_generator_ids": validated[
            "selected_local_generator_ids"
        ],
        "mask_sha256": validated["mask_sha256"],
        "row_level_donor_allowlist_enforced": True,
        "fresh_optimizer": True,
        "official_validation_used": False,
        "official_test_used": False,
    }
    for field, expected in required.items():
        if provenance.get(field) != expected:
            raise ValueError(f"checkpoint fixed-mask provenance {field} mismatch")
    if integrity.get("mask_sha256") != validated["mask_sha256"]:
        raise ValueError("checkpoint integrity mask SHA256 mismatch")
    state = payload.get("system_state_dict")
    if not isinstance(state, Mapping):
        raise ValueError("checkpoint lacks system_state_dict")
    mask_tensors = [
        value
        for name, value in state.items()
        if str(name).endswith("decoder.generator_active_mask")
    ]
    if len(mask_tensors) != 1:
        raise ValueError("checkpoint must contain exactly one generator mask buffer")
    mask = torch.as_tensor(mask_tensors[0]).detach().cpu()
    expected_mask = torch.zeros(validated["candidate_count"], dtype=torch.float32)
    expected_mask[validated["selected_local_generator_ids"]] = 1.0
    if mask.shape != expected_mask.shape or not torch.equal(
        mask.to(dtype=torch.float32), expected_mask
    ):
        raise ValueError("checkpoint-persisted generator mask differs from job")
    return dict(provenance)


def _write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(dict(value), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def run_backend(
    *,
    job_path: str | Path,
    out_path: str | Path,
    work_dir: str | Path,
    validate_only: bool,
) -> dict[str, Any]:
    job_source = Path(job_path).resolve()
    output = Path(out_path).resolve()
    work = Path(work_dir).resolve()
    if output.exists() or work.exists():
        raise FileExistsError("training backend output/work path already exists")
    job = read_json_object(job_source, label="fixed-mask training job")
    preflight = build_training_preflight(
        job, job_path=job_source, work_dir=work / "model"
    )
    if validate_only:
        _write_json_exclusive(output, preflight)
        return preflight

    work.mkdir(parents=True, exist_ok=False)
    model_dir = work / "model"
    config = build_run_current_config(
        job, job_path=job_source, work_dir=model_dir
    )
    config_path = work / "generated_train_config.json"
    _write_json_exclusive(config_path, config)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "kmlee_bam.training.run_current",
            "--config",
            str(config_path),
        ],
        check=False,
        stdin=subprocess.DEVNULL,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"run_current fixed-mask training exited with {completed.returncode}"
        )
    checkpoint = model_dir / "checkpoint_last.pt"
    audit_training_checkpoint(checkpoint, job)
    validated = validate_fixed_mask_training_job(job)
    result: dict[str, Any] = {
        "schema_version": TRAINING_OUTPUT_SCHEMA,
        "job_sha256": validated["job_sha256"],
        "donor_partition": validated["donor_partition"],
        "candidate_count": validated["candidate_count"],
        "generator_count": validated["generator_count"],
        "selected_local_generator_ids": validated[
            "selected_local_generator_ids"
        ],
        "mask_sha256": validated["mask_sha256"],
        "training_donor_names_sha256": validated[
            "training_donor_names_sha256"
        ],
        "source_sha256": validated["source_sha256"],
        "fixed_mask_for_entire_run": True,
        "fixed_mask_persisted_in_checkpoint": True,
        "fresh_optimizer": True,
        "same_seed_schedule_across_k": True,
        "row_level_donor_allowlist_enforced": True,
        "official_validation_used": False,
        "official_test_used": False,
        "checkpoint_path": str(checkpoint.resolve()),
        "checkpoint_sha256": sha256_file(checkpoint),
        "generated_config_path": str(config_path.resolve()),
        "generated_config_sha256": sha256_file(config_path),
    }
    result["output_sha256"] = canonical_json_sha256(result)
    _write_json_exclusive(output, result)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate or run one sealed fixed-mask training job"
    )
    parser.add_argument("--job-json", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    run_backend(
        job_path=args.job_json,
        out_path=args.out_json,
        work_dir=args.work_dir,
        validate_only=bool(args.validate_only),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
