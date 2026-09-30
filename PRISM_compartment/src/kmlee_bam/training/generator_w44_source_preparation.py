"""Fail-closed source preparation for canonical PRISM-v4 generator search.

This module deliberately does not pretend that a path plus a ``fit_partition``
label proves which donors were used.  A preparation bundle is non-executable
until every producer can consume the sealed W44/R10/C10 manifest as a row-level
allow-list and embed the *observed* fitting donor names in its output artifact.

All five consumers have opt-in manifest plumbing, but a newly written bundle is
still non-executable until concrete configs are materialized, reviewed, and its
outputs pass :func:`audit_completed_w44_sources`.  The bundle therefore keeps
``executor_ready=false`` and has no launch command.  No training is started by
this module.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from kmlee_bam.training.generator_exact_confirmation import (
    canonical_json_sha256,
    validate_canonical_v4_donor_split_manifest,
)
from kmlee_bam.training.sealed_donor_allowlist import (
    load_sealed_donor_allowlist,
    validate_w44_source_checkpoint,
)


SOURCE_PLAN_SCHEMA = "kmlee_bam.generator_w44_source_plan.v4"
SOURCE_AUDIT_SCHEMA = "kmlee_bam.generator_w44_source_audit.v4"
ALL64_ARTIFACT_SCHEMA = "kmlee_bam.generator_all64_artifacts.v4"

W44_ARTIFACTS = (
    "source_checkpoint",
    "age_normalization",
    "reference_scaler",
    "module_stats",
    "donor_context",
)

CONSUMER_CONNECTIONS: dict[str, str] = {
    "source_checkpoint": (
        "run_current canonical_w44_source restricts actual optimizer rows and "
        "embeds checkpoint provenance"
    ),
    "age_normalization": (
        "the same source hook refits age moments on observed W44 ages and "
        "writes a hash-bound JSON artifact"
    ),
    "reference_scaler": (
        "build_celltype_control_scaler canonical opt-in masks moment rows to W44"
    ),
    "module_stats": (
        "module-stat locations/scales fit W44 while raw R10/C10 targets remain "
        "evaluation-only"
    ),
    "donor_context": (
        "context extraction validates an attested W44 checkpoint; residualizer "
        "fit is W44 and support rows are train64-only"
    ),
}

CONSUMER_OPT_IN: dict[str, dict[str, Any]] = {
    "source_checkpoint": {
        "config_section": "canonical_w44_source",
        "enabled": True,
        "required_fields": (
            "donor_split_manifest_path",
            "age_normalization_out_path",
        ),
        "required_training_controls": {
            "train.early_stopping": False,
            "train.save_best": False,
            "train.skip_final_test": True,
        },
    },
    "age_normalization": {
        "producer": "source_checkpoint_job",
        "output_field": "canonical_w44_source.age_normalization_out_path",
    },
    "reference_scaler": {
        "config_fields": {
            "donor_split_manifest_path": "sealed_manifest_path",
            "canonical_fit_partition": "W44",
        }
    },
    "module_stats": {
        "cli_flags": (
            "--donor-split-manifest",
            "sealed_manifest_path",
            "--canonical-fit-partition",
            "W44",
        )
    },
    "donor_context": {
        "cli_flags": (
            "--donor-split-manifest",
            "sealed_manifest_path",
            "--canonical-fit-partition",
            "W44",
            "--source-checkpoint-expected-sha256",
            "attested_W44_encoder_checkpoint_sha256",
        )
    },
}


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_sha256(value: Any) -> bool:
    text = str(value).lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _load_json_object(path: str | Path, *, label: str) -> dict[str, Any]:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"failed to read {label}: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def _exact_names(
    observed: Sequence[Any],
    expected: Sequence[str],
    *,
    label: str,
) -> tuple[str, ...]:
    names = tuple(str(value) for value in observed)
    expected_names = tuple(str(value) for value in expected)
    if names != expected_names:
        raise ValueError(
            f"{label} must equal the sealed W44 donor names in manifest order"
        )
    if len(names) != 44 or len(set(names)) != 44:
        raise ValueError(f"{label} must contain exactly 44 unique donors")
    return names


def _npz_scalar(archive: Any, key: str) -> Any:
    if key not in archive.files:
        raise ValueError(f"artifact is missing embedded field {key!r}")
    value = np.asarray(archive[key])
    if value.size != 1:
        raise ValueError(f"embedded field {key!r} must be scalar")
    return value.reshape(()).item()


def build_w44_source_preparation_plan(
    *,
    base_config_path: str | Path,
    donor_split_manifest_path: str | Path,
    output_root: str | Path,
) -> dict[str, Any]:
    """Build a deterministic, non-runnable W44 source-production plan."""

    config_path = Path(base_config_path).resolve()
    split_path = Path(donor_split_manifest_path).resolve()
    output = Path(output_root).resolve()
    base_config = _load_json_object(config_path, label="base config")
    split = validate_canonical_v4_donor_split_manifest(
        _load_json_object(split_path, label="donor split manifest")
    )
    weight_names = tuple(str(value) for value in split["weight_donor_names"])
    _exact_names(weight_names, weight_names, label="weight_donor_names")

    required_sources: dict[str, dict[str, str]] = {}
    path_specs = (
        ("zarr", base_config.get("data", {}).get("zarr_path")),
        ("ordinal_spec", base_config.get("data", {}).get("spec_path")),
        (
            "module_registry",
            base_config.get("module_tokenizer", {}).get("registry_json_path"),
        ),
        (
            "activity_dictionary",
            base_config.get("module_tokenizer", {}).get("activity_weight_path"),
        ),
    )
    for name, raw_path in path_specs:
        if not isinstance(raw_path, str) or not raw_path:
            raise ValueError(f"base config is missing required source path {name}")
        source_path = Path(raw_path).expanduser().resolve()
        # Zarr is a directory and cannot be represented by one file digest here.
        # Its immutable dataset/spec identity must therefore be supplied by the
        # future row-filtering producer.  Ordinary files are bound immediately.
        entry = {"path": str(source_path)}
        if source_path.is_file():
            entry["sha256"] = _sha256_file(source_path)
        else:
            entry["sha256"] = "producer_must_emit_dataset_content_sha256"
        required_sources[name] = entry

    artifact_paths = {
        "source_checkpoint": output / "artifacts" / "w44_supernet_checkpoint.pt",
        "age_normalization": output / "artifacts" / "w44_age_normalization.json",
        "reference_scaler": output / "artifacts" / "w44_reference_scaler.npz",
        "module_stats": output / "artifacts" / "w44_module_rescue_stats.npz",
        "donor_context": output / "artifacts" / "w44_donor_context.npz",
    }
    jobs: list[dict[str, Any]] = []
    for artifact in W44_ARTIFACTS:
        jobs.append(
            {
                "artifact": artifact,
                "fit_partition": "W44",
                "sealed_donor_split_manifest_path": str(split_path),
                "sealed_donor_split_manifest_sha256": split["manifest_sha256"],
                "required_observed_fit_donor_names": list(weight_names),
                "required_observed_fit_donor_sha256": split[
                    "weight_donor_sha256"
                ],
                "output_path": str(artifact_paths[artifact]),
                "row_level_allowlist_required": True,
                "consumer_connected": True,
                "consumer_connection": CONSUMER_CONNECTIONS[artifact],
                "required_opt_in": CONSUMER_OPT_IN[artifact],
                "executor_ready": False,
                "launch_command": None,
                "blocker": (
                    "concrete artifact/config has not been built and audited; "
                    "user review is required before launch"
                ),
            }
        )

    plan: dict[str, Any] = {
        "schema_version": SOURCE_PLAN_SCHEMA,
        "purpose": "canonical_v4_generator_selection_W44_sources",
        "base_config_path": str(config_path),
        "base_config_sha256": _sha256_file(config_path),
        "donor_split_manifest_path": str(split_path),
        "donor_split_manifest_sha256": split["manifest_sha256"],
        "weight_donor_names": list(weight_names),
        "weight_donor_sha256": split["weight_donor_sha256"],
        "ranking_donor_names": list(split["ranking_donor_names"]),
        "confirmation_donor_names": list(split["confirmation_donor_names"]),
        "required_sources": required_sources,
        "jobs": jobs,
        "executor_ready": False,
        "launch_allowed": False,
        "full_training_started": False,
        "blocked_artifacts": list(W44_ARTIFACTS),
        "unconnected_consumers": [],
        "consumer_connections": dict(CONSUMER_CONNECTIONS),
        "review_required_before_any_full_training": True,
    }
    plan["manifest_sha256"] = canonical_json_sha256(plan)
    return plan


def write_w44_source_preparation_bundle(
    *,
    base_config_path: str | Path,
    donor_split_manifest_path: str | Path,
    output_root: str | Path,
    force: bool = False,
) -> tuple[Path, ...]:
    """Write only review artifacts; never execute a producer or trainer."""

    output = Path(output_root).resolve()
    plan = build_w44_source_preparation_plan(
        base_config_path=base_config_path,
        donor_split_manifest_path=donor_split_manifest_path,
        output_root=output,
    )
    paths = [output / "w44_source_plan.json"] + [
        output / "jobs" / f"{job['artifact']}.job.json" for job in plan["jobs"]
    ]
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        raise FileExistsError(
            "refusing to overwrite W44 preparation bundle: "
            + ", ".join(str(path) for path in existing)
        )
    (output / "jobs").mkdir(parents=True, exist_ok=True)
    paths[0].write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for path, job in zip(paths[1:], plan["jobs"]):
        payload = {
            "schema_version": "kmlee_bam.generator_w44_source_job.v4",
            **job,
            "source_plan_sha256": plan["manifest_sha256"],
        }
        payload["manifest_sha256"] = canonical_json_sha256(payload)
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return tuple(paths)


def _validate_embedded_common(
    payload: Mapping[str, Any],
    *,
    split: Mapping[str, Any],
    donor_field: str,
    artifact: str,
) -> None:
    _exact_names(
        payload.get(donor_field, ()),
        split["weight_donor_names"],
        label=f"{artifact}.{donor_field}",
    )
    if payload.get("donor_split_manifest_sha256") != split["manifest_sha256"]:
        raise ValueError(f"{artifact} split-manifest SHA256 mismatch")
    if payload.get("official_validation_used") is not False:
        raise ValueError(f"{artifact} must embed official_validation_used=false")
    if payload.get("official_test_used") is not False:
        raise ValueError(f"{artifact} must embed official_test_used=false")
    if payload.get("row_level_allowlist_enforced") is not True:
        raise ValueError(f"{artifact} lacks observed row-level allow-list evidence")


def audit_completed_w44_sources(
    *,
    donor_split_manifest_path: str | Path,
    source_checkpoint_path: str | Path,
    age_normalization_path: str | Path,
    reference_scaler_path: str | Path,
    module_stats_path: str | Path,
    donor_context_path: str | Path,
) -> dict[str, Any]:
    """Open and audit completed sources; current legacy artifacts fail closed."""

    split = validate_canonical_v4_donor_split_manifest(
        _load_json_object(donor_split_manifest_path, label="donor split manifest")
    )
    checkpoint_path = Path(source_checkpoint_path).resolve()
    allowlist = load_sealed_donor_allowlist(
        donor_split_manifest_path,
        partition="W44",
    )
    checkpoint_sha, _ = validate_w44_source_checkpoint(
        checkpoint_path, allowlist
    )

    age_path = Path(age_normalization_path).resolve()
    age = _load_json_object(age_path, label="W44 age normalization")
    _exact_names(
        age.get("allowed_fit_donor_names", ()),
        split["weight_donor_names"],
        label="age_normalization.allowed_fit_donor_names",
    )
    if age.get("donor_split_manifest_sha256") != split["manifest_sha256"]:
        raise ValueError("age_normalization split-manifest SHA256 mismatch")
    if age.get("row_level_allowlist_enforced") is not True:
        raise ValueError("age_normalization lacks row-level allow-list evidence")
    if age.get("official_validation_used") is not False or age.get(
        "official_test_used"
    ) is not False:
        raise ValueError("age_normalization used an official held-out split")
    observed_age_names = tuple(
        str(value) for value in age.get("observed_age_donor_names", ())
    )
    if (
        len(observed_age_names) < 2
        or not set(observed_age_names).issubset(split["weight_donor_names"])
    ):
        raise ValueError("age moments were not an observed W44-only subset")
    supplied_age_sha = age.get("manifest_sha256")
    expected_age_sha = canonical_json_sha256(
        {key: value for key, value in age.items() if key != "manifest_sha256"}
    )
    if supplied_age_sha != expected_age_sha:
        raise ValueError("age normalization manifest SHA256 mismatch")
    if not np.isfinite(float(age.get("age_mean", np.nan))) or not np.isfinite(
        float(age.get("age_standard_deviation", np.nan))
    ):
        raise ValueError("W44 age normalization contains non-finite moments")
    if float(age["age_standard_deviation"]) <= 0.0:
        raise ValueError("W44 age standard deviation must be positive")

    npz_specs = (
        (
            "reference_scaler",
            Path(reference_scaler_path).resolve(),
            "observed_moment_donor_names",
            False,
        ),
        (
            "module_stats",
            Path(module_stats_path).resolve(),
            "observed_statistic_donor_names",
            True,
        ),
        (
            "donor_context",
            Path(donor_context_path).resolve(),
            "observed_residualizer_fit_donor_names",
            True,
        ),
    )
    artifact_sha = {
        "source_checkpoint": checkpoint_sha,
        "age_normalization": _sha256_file(age_path),
    }
    context_extractor_checkpoint_sha256 = None
    context_extractor_provenance_sha256 = None
    for artifact, path, donor_field, require_all in npz_specs:
        with np.load(path, allow_pickle=True) as archive:
            allowed_names = (
                np.asarray(archive["allowed_fit_donor_names"], dtype=object)
                .astype(str)
                .tolist()
                if "allowed_fit_donor_names" in archive.files
                else ()
            )
            _exact_names(
                allowed_names,
                split["weight_donor_names"],
                label=f"{artifact}.allowed_fit_donor_names",
            )
            embedded = {
                donor_field: np.asarray(archive[donor_field], dtype=object)
                .astype(str)
                .tolist()
                if donor_field in archive.files
                else (),
                "donor_split_manifest_sha256": str(
                    _npz_scalar(archive, "donor_split_manifest_sha256")
                ),
                "official_validation_used": bool(
                    _npz_scalar(archive, "official_validation_used")
                ),
                "official_test_used": bool(
                    _npz_scalar(archive, "official_test_used")
                ),
                "row_level_allowlist_enforced": bool(
                    _npz_scalar(archive, "row_level_allowlist_enforced")
                ),
            }
            if require_all:
                _validate_embedded_common(
                    embedded,
                    split=split,
                    donor_field=donor_field,
                    artifact=artifact,
                )
            else:
                observed_names = tuple(
                    str(value) for value in embedded[donor_field]
                )
                if not observed_names or not set(observed_names).issubset(
                    split["weight_donor_names"]
                ):
                    raise ValueError(
                        "reference_scaler moment donors must be a non-empty "
                        "W44-only subset"
                    )
                for flag, expected in (
                    ("donor_split_manifest_sha256", split["manifest_sha256"]),
                    ("official_validation_used", False),
                    ("official_test_used", False),
                    ("row_level_allowlist_enforced", True),
                ):
                    if embedded.get(flag) != expected:
                        raise ValueError(
                            f"reference_scaler requires {flag}={expected!r}"
                        )
            if artifact == "donor_context":
                observed_checkpoint_sha = str(
                    _npz_scalar(archive, "extractor_checkpoint_sha256")
                )
                observed_checkpoint_provenance_sha = str(
                    _npz_scalar(
                        archive,
                        "extractor_checkpoint_provenance_sha256",
                    )
                )
                if not _is_sha256(observed_checkpoint_sha) or not _is_sha256(
                    observed_checkpoint_provenance_sha
                ):
                    raise ValueError(
                        "donor_context lacks a hash-bound attested W44 extractor "
                        "checkpoint"
                    )
                context_extractor_checkpoint_sha256 = observed_checkpoint_sha
                context_extractor_provenance_sha256 = (
                    observed_checkpoint_provenance_sha
                )
                source_fit_names = (
                    np.asarray(
                        archive["observed_source_model_fit_donor_names"],
                        dtype=object,
                    )
                    .astype(str)
                    .tolist()
                    if "observed_source_model_fit_donor_names" in archive.files
                    else ()
                )
                _exact_names(
                    source_fit_names,
                    split["weight_donor_names"],
                    label="donor_context.observed_source_model_fit_donor_names",
                )
        artifact_sha[artifact] = _sha256_file(path)

    audit: dict[str, Any] = {
        "schema_version": SOURCE_AUDIT_SCHEMA,
        "donor_split_manifest_sha256": split["manifest_sha256"],
        "observed_weight_donor_names": list(split["weight_donor_names"]),
        "observed_weight_donor_sha256": split["weight_donor_sha256"],
        "artifact_sha256": artifact_sha,
        "context_extractor_checkpoint_sha256": (
            context_extractor_checkpoint_sha256
        ),
        "context_extractor_checkpoint_provenance_sha256": (
            context_extractor_provenance_sha256
        ),
        "executor_ready": True,
        "launch_allowed": False,
        "reason_launch_is_still_blocked": (
            "artifact audit does not replace user review/approval"
        ),
    }
    audit["manifest_sha256"] = canonical_json_sha256(audit)
    return audit


def validate_all64_artifact_manifest(
    manifest: Mapping[str, Any],
    *,
    split_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate final-retrain artifacts against the complete train64 set."""

    result = dict(manifest)
    if result.get("schema_version") != ALL64_ARTIFACT_SCHEMA:
        raise ValueError("all64 artifact manifest schema mismatch")
    split = validate_canonical_v4_donor_split_manifest(split_manifest)
    expected = {
        str(value)
        for prefix in ("weight", "ranking", "confirmation")
        for value in split[f"{prefix}_donor_names"]
    }
    if len(expected) != 64:
        raise RuntimeError("sealed donor split does not contain 64 unique donors")
    for artifact in ("reference_scaler", "module_stats", "donor_context"):
        entry = result.get(artifact)
        if not isinstance(entry, Mapping):
            raise ValueError(f"all64 manifest lacks {artifact}")
        path = Path(str(entry.get("path", ""))).resolve()
        digest = str(entry.get("sha256", "")).lower()
        if not _is_sha256(digest) or _sha256_file(path) != digest:
            raise ValueError(f"all64 {artifact} SHA256 mismatch")
        # Do not accept donor names written only in this JSON sidecar.  Open
        # the artifact and require the producer's row-level observation to be
        # embedded in the hashed NPZ itself.
        with np.load(path, allow_pickle=True) as archive:
            if "observed_fit_donor_names" not in archive.files:
                raise ValueError(
                    f"all64 {artifact} lacks embedded fit-donor evidence"
                )
            embedded_names = tuple(
                np.asarray(
                    archive["observed_fit_donor_names"], dtype=object
                )
                .astype(str)
                .tolist()
            )
            embedded_enforced = bool(
                _npz_scalar(archive, "row_level_allowlist_enforced")
            )
            embedded_validation = bool(
                _npz_scalar(archive, "official_validation_used")
            )
            embedded_test = bool(
                _npz_scalar(archive, "official_test_used")
            )
        names = tuple(str(value) for value in entry.get("observed_fit_donor_names", ()))
        if names != embedded_names:
            raise ValueError(
                f"all64 {artifact} sidecar donor names differ from the artifact"
            )
        if len(embedded_names) != 64 or set(embedded_names) != expected:
            raise ValueError(f"all64 {artifact} was not fit on exactly train64")
        if (
            entry.get("row_level_allowlist_enforced") is not True
            or embedded_enforced is not True
        ):
            raise ValueError(f"all64 {artifact} lacks row-level fit evidence")
        if embedded_validation or embedded_test:
            raise ValueError(
                f"all64 {artifact} embedded official held-out usage"
            )
    if result.get("official_validation_used") is not False:
        raise ValueError("all64 artifacts must not use official validation")
    if result.get("official_test_used") is not False:
        raise ValueError("all64 artifacts must not use official test")
    supplied = result.get("manifest_sha256")
    expected_sha = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "manifest_sha256"}
    )
    if supplied != expected_sha:
        raise ValueError("all64 artifact manifest SHA256 mismatch")
    return result


def build_all64_fixed_mask_retrain_config(
    *,
    base_config: Mapping[str, Any],
    donor_split_manifest: Mapping[str, Any],
    all64_artifact_manifest: Mapping[str, Any],
    source_checkpoint_path: str | Path,
    source_checkpoint_sha256: str,
    fixed_generator_result_path: str | Path,
    fixed_generator_result_sha256: str,
    out_dir: str | Path,
) -> dict[str, Any]:
    """Return final-train config semantics after all64 artifacts are audited."""

    artifacts = validate_all64_artifact_manifest(
        all64_artifact_manifest,
        split_manifest=donor_split_manifest,
    )
    checkpoint_path = Path(source_checkpoint_path).resolve()
    result_path = Path(fixed_generator_result_path).resolve()
    for label, path, digest in (
        ("source checkpoint", checkpoint_path, source_checkpoint_sha256),
        ("fixed generator result", result_path, fixed_generator_result_sha256),
    ):
        if not _is_sha256(digest) or _sha256_file(path) != str(digest).lower():
            raise ValueError(f"{label} SHA256 mismatch")

    config = deepcopy(dict(base_config))
    config.setdefault("data", {})["reference_scaler_path"] = str(
        Path(artifacts["reference_scaler"]["path"]).resolve()
    )
    config.setdefault("precision_medicine", {})["context_npz_path"] = str(
        Path(artifacts["donor_context"]["path"]).resolve()
    )
    config.setdefault("prism_module_rescue", {})["stats_path"] = str(
        Path(artifacts["module_stats"]["path"]).resolve()
    )
    config["prism_module_rescue"]["expected_stats_sha256"] = artifacts[
        "module_stats"
    ]["sha256"]
    config.setdefault("learned_generator_count", {})["enabled"] = False
    config.setdefault("generator_budget_search", {})["enabled"] = False
    config["warm_start"] = {
        "init_weights_path": str(checkpoint_path),
        "expected_sha256": str(source_checkpoint_sha256).lower(),
        "require_precision_head": True,
        "fixed_generator_result_path": str(result_path),
        "fixed_generator_result_expected_sha256": str(
            fixed_generator_result_sha256
        ).lower(),
        # build_system has already refit the context residualizer on train64;
        # these exact buffers must survive loading the W44 checkpoint.
        "rebuild_precision_context_buffers": True,
        "rebuilt_precision_context_fit_partition": "train64",
        "rebuilt_precision_context_path": artifacts["donor_context"]["path"],
        "rebuilt_precision_context_expected_sha256": artifacts[
            "donor_context"
        ]["sha256"],
    }
    config.setdefault("train", {})["out_dir"] = str(Path(out_dir).resolve())
    config["_canonical_v4_final_retrain"] = {
        "schema_version": "kmlee_bam.generator_fixed_mask_retrain.v4",
        "donor_split_manifest_sha256": donor_split_manifest["manifest_sha256"],
        "all64_artifact_manifest_sha256": artifacts["manifest_sha256"],
        "context_buffers_built_before_warm_start": True,
        "w44_context_buffers_must_not_load": True,
        "fresh_optimizer": True,
        "user_approval_required_before_launch": True,
    }
    config["_launch_guard"] = {
        "launch_allowed": False,
        "reason": (
            "review-only final fixed-mask config; explicit user approval is "
            "required before full training"
        ),
    }
    return config
