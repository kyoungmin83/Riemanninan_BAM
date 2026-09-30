"""Frozen-additive metric backend for exact generator confirmation.

The model-facing collector produces one deterministic donor feature vector
from observed, full-decoder and generator-isolated module pseudobulks.  Every
scientific choice (centres, disease directions, worst-quartile coordinates and
leakage probes) must be fitted on R10 and encoded as an immutable affine metric
program.  Consequently C10 only applies frozen arithmetic and never fits a
probe or chooses a coordinate.

The current repository has no producer for that R10 program.  Preflight thus
fails closed on historical placeholder ``frozen_signature`` files rather than
inventing definitions, especially for the two under-specified leakage metrics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from kmlee_bam.training.generator_exact_confirmation import (
    CANONICAL_V4_METRIC_NAMES,
    canonical_json_sha256,
)
from kmlee_bam.training.generator_exact_evaluation import (
    JOB_SCHEMA,
    OUTPUT_SCHEMA,
    metric_contract,
)
from kmlee_bam.training.generator_exact_executor import (
    RUNTIME_JOB_SCHEMA,
    TRAINING_OUTPUT_SCHEMA,
)


FROZEN_PROGRAM_SCHEMA = "kmlee_bam.generator_frozen_additive_program.v1"
DENOMINATOR_SCHEMA = "kmlee_bam.generator_metric_denominator.v1"
METRIC_BACKEND_PREFLIGHT_SCHEMA = (
    "kmlee_bam.generator_exact_metric_backend_preflight.v1"
)
SUPPORTED_METRIC_KINDS = frozenset({"huber", "direct", "binary_logloss"})


class ProductionMetricBlocker(RuntimeError):
    """A missing scientific definition, not a transient runtime failure."""


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: str | Path, *, label: str) -> dict[str, Any]:
    source = Path(path).resolve()
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"failed to read {label}: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return value


def _require_digest(value: Any, *, label: str) -> str:
    digest = str(value).lower()
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise ValueError(f"{label} must be a lowercase SHA256 digest")
    return digest


def validate_runtime_job(runtime: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(runtime)
    if result.get("schema_version") != RUNTIME_JOB_SCHEMA:
        raise ValueError("metric backend runtime-job schema mismatch")
    supplied = result.get("runtime_job_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "runtime_job_sha256"}
    )
    if supplied != expected:
        raise ValueError("metric backend runtime-job SHA256 mismatch")
    job = result.get("evaluation_job")
    training = result.get("training_output")
    if not isinstance(job, Mapping) or job.get("schema_version") != JOB_SCHEMA:
        raise ValueError("runtime lacks a valid evaluation job")
    if job.get("job_sha256") != canonical_json_sha256(
        {key: value for key, value in job.items() if key != "job_sha256"}
    ):
        raise ValueError("runtime evaluation-job SHA256 mismatch")
    if not isinstance(training, Mapping) or training.get(
        "schema_version"
    ) != TRAINING_OUTPUT_SCHEMA:
        raise ValueError("runtime lacks a fixed-mask training output")
    if training.get("job_sha256") != job.get("upstream_training_job_sha256"):
        raise ValueError("runtime training/evaluation dependency mismatch")
    for field in (
        "donor_partition",
        "candidate_count",
        "generator_count",
        "selected_local_generator_ids",
        "mask_sha256",
        "source_sha256",
    ):
        if training.get(field) != job.get(field):
            raise ValueError(f"runtime training/evaluation {field} mismatch")
    checkpoint = Path(str(result.get("evaluated_checkpoint_path", ""))).resolve()
    if not checkpoint.is_file():
        raise ValueError("runtime evaluated checkpoint does not exist")
    checkpoint_sha = _sha256_file(checkpoint)
    if checkpoint_sha != result.get("evaluated_checkpoint_sha256") or (
        checkpoint_sha != training.get("checkpoint_sha256")
    ):
        raise ValueError("runtime evaluated checkpoint SHA256 mismatch")
    if result.get("upstream_training_manifest_sha256") != training.get(
        "output_sha256"
    ):
        raise ValueError("runtime upstream training manifest mismatch")
    for name, raw_path in job.get("source_paths", {}).items():
        source = Path(str(raw_path)).resolve()
        if not source.is_file() or _sha256_file(source) != job.get(
            "source_sha256", {}
        ).get(name):
            raise ValueError(f"runtime source {name} path/SHA256 mismatch")
    return result


def load_metric_denominators(
    path: str | Path,
    *,
    expected_sha256: str,
    donor_count: int,
) -> dict[str, int]:
    source = Path(path).resolve()
    if _sha256_file(source) != _require_digest(
        expected_sha256, label="metric_denominator_sha256"
    ):
        raise ValueError("metric denominator raw SHA256 mismatch")
    payload = _read_json(source, label="metric denominator")
    if payload.get("schema_version") != DENOMINATOR_SCHEMA:
        raise ProductionMetricBlocker(
            "metric denominator artifact has no canonical v1 schema"
        )
    if payload.get("fit_partition") != "R10":
        raise ValueError("metric denominators must be frozen on R10")
    if payload.get("official_validation_used") is not False or payload.get(
        "official_test_used"
    ) is not False:
        raise ValueError("metric denominator used a forbidden official split")
    values = payload.get("denominators")
    if not isinstance(values, Mapping) or set(values) != CANONICAL_V4_METRIC_NAMES:
        raise ValueError("metric denominator names differ from canonical registry")
    parsed: dict[str, int] = {}
    for name, raw in values.items():
        if isinstance(raw, bool) or not isinstance(raw, int) or not 1 <= raw <= donor_count:
            raise ValueError(f"invalid metric denominator for {name}")
        parsed[str(name)] = int(raw)
    supplied = payload.get("manifest_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != "manifest_sha256"}
    )
    if supplied != expected:
        raise ValueError("metric denominator manifest SHA256 mismatch")
    return parsed


def load_frozen_metric_program(
    path: str | Path,
    *,
    expected_sha256: str,
    expected_metric_registry_sha256: str,
) -> dict[str, Any]:
    source = Path(path).resolve()
    if _sha256_file(source) != _require_digest(
        expected_sha256, label="frozen_signature_sha256"
    ):
        raise ValueError("frozen metric program raw SHA256 mismatch")
    try:
        archive_context = np.load(source, allow_pickle=False)
    except Exception as exc:
        raise ProductionMetricBlocker(
            "frozen_signature is not a canonical NPZ metric program"
        ) from exc
    with archive_context as archive:
        required = {
            "schema_version",
            "metadata_json",
            "metric_names",
            "feature_names",
            "metric_kind",
            "prediction_weight",
            "target_weight",
            "prediction_bias",
            "target_bias",
            "huber_delta",
            "support_feature_index",
            "support_minimum",
        }
        missing = sorted(required - set(archive.files))
        if missing:
            raise ProductionMetricBlocker(
                "frozen metric program lacks arrays: " + ", ".join(missing)
            )
        schema = str(np.asarray(archive["schema_version"]).reshape(()).item())
        if schema != FROZEN_PROGRAM_SCHEMA:
            raise ProductionMetricBlocker(
                f"frozen metric program schema is {schema!r}, expected "
                f"{FROZEN_PROGRAM_SCHEMA!r}"
            )
        try:
            metadata = json.loads(
                str(np.asarray(archive["metadata_json"]).reshape(()).item())
            )
        except Exception as exc:
            raise ValueError("frozen metric metadata_json is invalid") from exc
        metric_names = tuple(np.asarray(archive["metric_names"]).astype(str).tolist())
        feature_names = tuple(
            np.asarray(archive["feature_names"]).astype(str).tolist()
        )
        metric_kind = tuple(np.asarray(archive["metric_kind"]).astype(str).tolist())
        arrays = {
            name: np.asarray(archive[name])
            for name in required
            if name
            not in {"schema_version", "metadata_json", "metric_names", "feature_names", "metric_kind"}
        }
    canonical_order = tuple(sorted(CANONICAL_V4_METRIC_NAMES))
    if metric_names != canonical_order:
        raise ValueError("frozen metric program metric order is not canonical")
    if len(set(feature_names)) != len(feature_names) or not feature_names:
        raise ValueError("frozen metric feature names must be unique/non-empty")
    n_metric, n_feature = len(metric_names), len(feature_names)
    if metric_kind.__len__() != n_metric or not set(metric_kind).issubset(
        SUPPORTED_METRIC_KINDS
    ):
        raise ValueError("frozen metric kinds are invalid")
    expected_shapes = {
        "prediction_weight": (n_metric, n_feature),
        "target_weight": (n_metric, n_feature),
        "prediction_bias": (n_metric,),
        "target_bias": (n_metric,),
        "huber_delta": (n_metric,),
        "support_feature_index": (n_metric,),
        "support_minimum": (n_metric,),
    }
    for name, shape in expected_shapes.items():
        if arrays[name].shape != shape or not np.isfinite(
            arrays[name].astype(np.float64)
        ).all():
            raise ValueError(f"frozen metric array {name} shape/finite mismatch")
    support_index = arrays["support_feature_index"].astype(np.int64)
    if np.any((support_index < -1) | (support_index >= n_feature)):
        raise ValueError("frozen metric support_feature_index is out of range")
    if metadata.get("schema_version") != FROZEN_PROGRAM_SCHEMA:
        raise ValueError("frozen metric metadata schema mismatch")
    if metadata.get("fit_partition") != "R10":
        raise ValueError("frozen metric program must be fitted on R10")
    if metadata.get("metric_registry_sha256") != expected_metric_registry_sha256:
        raise ValueError("frozen metric program registry SHA256 mismatch")
    if metadata.get("feature_names_sha256") != canonical_json_sha256(
        list(feature_names)
    ):
        raise ValueError("frozen metric program feature-name SHA256 mismatch")
    for field in ("official_validation_used", "official_test_used"):
        if metadata.get(field) is not False:
            raise ValueError(f"frozen metric program requires {field}=false")
    return {
        "metadata": metadata,
        "metric_names": metric_names,
        "feature_names": feature_names,
        "metric_kind": metric_kind,
        **arrays,
    }


def evaluate_frozen_additive_program(
    donor_features: np.ndarray,
    *,
    feature_names: Sequence[str],
    program: Mapping[str, Any],
    denominators: Mapping[str, int],
) -> tuple[dict[str, list[float]], dict[str, int]]:
    """Apply an R10-frozen program; no fitting or selection occurs here."""

    features = np.asarray(donor_features, dtype=np.float64)
    if features.ndim != 2 or not np.isfinite(features).all():
        raise ValueError("donor feature matrix must be finite and two-dimensional")
    if tuple(str(value) for value in feature_names) != tuple(program["feature_names"]):
        raise ValueError("runtime donor features differ from frozen program order")
    n_donor = int(features.shape[0])
    prediction = features @ np.asarray(
        program["prediction_weight"], dtype=np.float64
    ).T + np.asarray(program["prediction_bias"], dtype=np.float64)[None, :]
    target = features @ np.asarray(
        program["target_weight"], dtype=np.float64
    ).T + np.asarray(program["target_bias"], dtype=np.float64)[None, :]
    support_index = np.asarray(program["support_feature_index"], dtype=np.int64)
    support_minimum = np.asarray(program["support_minimum"], dtype=np.float64)
    output: dict[str, list[float]] = {}
    support_counts: dict[str, int] = {}
    for metric_index, name in enumerate(program["metric_names"]):
        if support_index[metric_index] < 0:
            support = np.ones(n_donor, dtype=bool)
        else:
            support = (
                features[:, support_index[metric_index]]
                >= support_minimum[metric_index]
            )
        observed_support = int(support.sum())
        denominator = int(denominators[name])
        if observed_support != denominator:
            raise ValueError(
                f"runtime support for {name} differs from frozen denominator: "
                f"{observed_support} != {denominator}"
            )
        kind = program["metric_kind"][metric_index]
        residual = prediction[:, metric_index] - target[:, metric_index]
        if kind == "direct":
            loss = prediction[:, metric_index]
            if np.any(loss[support] < 0):
                raise ValueError(f"direct loss {name} contains negative values")
        elif kind == "binary_logloss":
            labels = target[:, metric_index]
            if np.any((labels[support] < 0) | (labels[support] > 1)):
                raise ValueError(f"binary leakage target for {name} is outside [0,1]")
            logits = prediction[:, metric_index]
            loss = np.maximum(logits, 0) - logits * labels + np.log1p(
                np.exp(-np.abs(logits))
            )
        else:
            delta = float(np.asarray(program["huber_delta"])[metric_index])
            if not math.isfinite(delta) or delta <= 0:
                raise ValueError(f"Huber delta for {name} must be positive")
            absolute = np.abs(residual)
            loss = np.where(
                absolute <= delta,
                0.5 * np.square(residual),
                delta * (absolute - 0.5 * delta),
            )
        # Exact confirmation takes an unweighted mean over ten donor entries.
        # Zero unsupported contributions and rescale by the frozen denominator
        # so that mean(output) equals the supported-donor mean loss.
        contribution = np.where(support, loss * n_donor / denominator, 0.0)
        if not np.isfinite(contribution).all():
            raise ValueError(f"additive contribution {name} is non-finite")
        output[name] = [float(value) for value in contribution.tolist()]
        support_counts[name] = observed_support
    return output, support_counts


def preflight_metric_backend(runtime: Mapping[str, Any]) -> dict[str, Any]:
    validated = validate_runtime_job(runtime)
    job = validated["evaluation_job"]
    blockers: list[str] = []
    if job.get("job_kind") == "ranking_evaluation_r10":
        blockers.append(
            "R10 job does not define the leakage target/probe family or a "
            "canonical producer for its frozen 24-metric program"
        )
    elif job.get("job_kind") == "one_shot_confirmation_evaluation_c10":
        contract = job.get("additive_confirmation_contract")
        if not isinstance(contract, Mapping):
            raise ValueError("C10 evaluation lacks additive confirmation contract")
        try:
            program = load_frozen_metric_program(
                job["source_paths"]["frozen_signature"],
                expected_sha256=contract["frozen_signature_sha256"],
                expected_metric_registry_sha256=job["metric_registry_sha256"],
            )
            load_metric_denominators(
                job["source_paths"]["metric_denominator"],
                expected_sha256=contract["metric_denominator_sha256"],
                donor_count=len(job["donor_names"]),
            )
            if not {
                "leakage.full.additive_probe_loss",
                "leakage.isolated.additive_probe_loss",
            }.issubset(program["metric_names"]):
                blockers.append("frozen program omits a canonical leakage probe")
        except ProductionMetricBlocker as exc:
            blockers.append(str(exc))
    else:
        raise ValueError("unsupported evaluation job kind")
    result: dict[str, Any] = {
        "schema_version": METRIC_BACKEND_PREFLIGHT_SCHEMA,
        "status": "ready" if not blockers else "blocked_before_data_or_gpu",
        "job_sha256": job["job_sha256"],
        "donor_partition": job["donor_partition"],
        "metric_registry_sha256": job["metric_registry_sha256"],
        "metric_count": len(metric_contract()),
        "blockers": blockers,
        "data_opened": False,
        "gpu_opened": False,
        "inference_started": False,
    }
    result["preflight_sha256"] = canonical_json_sha256(result)
    return result


def run_backend(
    *,
    runtime_job_path: str | Path,
    output_path: str | Path,
    preflight_only: bool,
) -> dict[str, Any]:
    runtime = _read_json(runtime_job_path, label="metric runtime job")
    preflight = preflight_metric_backend(runtime)
    if preflight_only:
        output = Path(output_path).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as handle:
            json.dump(preflight, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        return preflight
    if preflight["blockers"]:
        raise ProductionMetricBlocker("; ".join(preflight["blockers"]))
    # Deliberately unreachable until an audited R10 program producer exists.
    # Keeping this hard error is safer than silently substituting the legacy
    # pooled posthoc metrics, whose nonlinear aggregation violates the exact
    # donor-additive confirmation contract.
    raise ProductionMetricBlocker(
        "canonical model-to-donor feature collection is not released until "
        "the R10 frozen-program producer and leakage definition are approved"
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Preflight or run exact frozen-additive metric inference"
    )
    parser.add_argument("--runtime-job-json", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args(argv)
    run_backend(
        runtime_job_path=args.runtime_job_json,
        output_path=args.out_json,
        preflight_only=bool(args.preflight_only),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
