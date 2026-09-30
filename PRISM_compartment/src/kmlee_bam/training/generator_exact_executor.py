"""Fail-closed subprocess executor for exact fixed-mask evaluation plans.

The executor owns orchestration and integrity checks, while model-specific
training/inference backends own tensor computation.  A backend cannot change
the donor partition, generator mask, or upstream checkpoint unnoticed: every
input and output is hash-bound to the sealed plan.  Confirmation execution is
impossible until the C10 ledger atomically consumes the claimed partition.

Importing or dry-running this module never launches training or inference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Mapping, Sequence

from kmlee_bam.training.generator_c10_consumption_ledger import (
    begin_c10_execution_once,
    claim_exact_evaluation_plan,
    record_c10_confirmation_failure,
    verify_claimed_exact_evaluation_plan,
)
from kmlee_bam.training.generator_exact_confirmation import canonical_json_sha256
from kmlee_bam.training.generator_exact_evaluation import (
    JOB_SCHEMA,
    OUTPUT_SCHEMA,
    PLAN_SCHEMA,
    validate_exact_evaluation_output,
)


TRAINING_OUTPUT_SCHEMA = "kmlee_bam.generator_exact_training_output.v1"
RUNTIME_JOB_SCHEMA = "kmlee_bam.generator_exact_runtime_job.v1"
EXECUTION_MANIFEST_SCHEMA = "kmlee_bam.generator_exact_execution.v1"


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
        raise ValueError(f"{label} must be a JSON object")
    return value


def _write_json_exclusive(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(dict(payload), handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return path


def _validate_embedded_digest(
    payload: Mapping[str, Any], *, field: str, label: str
) -> str:
    supplied = str(payload.get(field, "")).lower()
    if len(supplied) != 64 or any(char not in "0123456789abcdef" for char in supplied):
        raise ValueError(f"{label}.{field} must be a SHA256 digest")
    expected = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != field}
    )
    if supplied != expected:
        raise ValueError(f"{label} {field} mismatch")
    return supplied


def validate_execution_plan(plan: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the complete dependency graph without touching donor data."""

    result = dict(plan)
    if result.get("schema_version") != PLAN_SCHEMA:
        raise ValueError("exact evaluation plan schema mismatch")
    _validate_embedded_digest(result, field="plan_sha256", label="plan")
    if result.get("source_split") != "train_only":
        raise ValueError("executor accepts train_only plans exclusively")
    for field in ("official_validation_used", "official_test_used"):
        if result.get(field) is not False:
            raise ValueError(f"executor requires {field}=false")
    partition = result.get("donor_partition")
    if partition not in {"ranking", "confirmation"}:
        raise ValueError("executor plan donor_partition is invalid")
    training = result.get("training_jobs")
    evaluation = result.get("evaluation_jobs")
    if not isinstance(training, list) or not isinstance(evaluation, list):
        raise ValueError("executor plan must contain training/evaluation jobs")
    if not training or len(training) != len(evaluation):
        raise ValueError("every evaluation job requires exactly one training job")

    training_by_sha: dict[str, dict[str, Any]] = {}
    for job in training:
        if not isinstance(job, Mapping) or job.get("schema_version") != JOB_SCHEMA:
            raise ValueError("training job schema mismatch")
        job = dict(job)
        digest = _validate_embedded_digest(
            job, field="job_sha256", label="training job"
        )
        if digest in training_by_sha:
            raise ValueError("duplicate training job SHA256")
        contract = job.get("training_contract")
        if not isinstance(contract, Mapping):
            raise ValueError("training job lacks training_contract")
        required_contract = {
            "fresh_optimizer": True,
            "fixed_mask_for_entire_run": True,
            "mask_trainable": False,
            "same_seed_schedule_across_k": True,
            "acute_mask_ablation": False,
            "learned_generator_count_enabled": False,
            "legacy_generator_budget_search_enabled": False,
            "official_validation_evaluation": False,
            "official_test_evaluation": False,
        }
        for name, expected in required_contract.items():
            if contract.get(name) != expected:
                raise ValueError(f"training contract requires {name}={expected!r}")
        donors = job.get("training_donor_names")
        expected_count = 44 if partition == "ranking" else 54
        if not isinstance(donors, list) or len(donors) != expected_count:
            raise ValueError(f"training job must contain exactly {expected_count} donors")
        if len(set(str(value) for value in donors)) != expected_count:
            raise ValueError("training donors must be unique")
        if job.get("training_donor_names_sha256") != canonical_json_sha256(donors):
            raise ValueError("training donor order SHA256 mismatch")
        training_by_sha[digest] = job

    seen_upstream: set[str] = set()
    for job in evaluation:
        if not isinstance(job, Mapping) or job.get("schema_version") != JOB_SCHEMA:
            raise ValueError("evaluation job schema mismatch")
        job = dict(job)
        _validate_embedded_digest(job, field="job_sha256", label="evaluation job")
        upstream = job.get("upstream_training_job_sha256")
        if upstream not in training_by_sha:
            raise ValueError("evaluation job references an unknown training job")
        if upstream in seen_upstream:
            raise ValueError("multiple evaluation jobs reuse one training job")
        seen_upstream.add(str(upstream))
        train_job = training_by_sha[str(upstream)]
        for field in (
            "donor_partition",
            "candidate_count",
            "generator_count",
            "selected_local_generator_ids",
            "mask_sha256",
            "source_sha256",
        ):
            if job.get(field) != train_job.get(field):
                raise ValueError(f"training/evaluation job {field} mismatch")
        expected_donors = 10
        donors = job.get("donor_names")
        if not isinstance(donors, list) or len(donors) != expected_donors:
            raise ValueError("evaluation job must contain exactly ten donors")
        if len(set(str(value) for value in donors)) != expected_donors:
            raise ValueError("evaluation donors must be unique")
        if set(donors) & set(train_job["training_donor_names"]):
            raise ValueError("training and evaluation donor partitions overlap")
        if job.get("donor_names_sha256") != canonical_json_sha256(donors):
            raise ValueError("evaluation donor order SHA256 mismatch")
    if seen_upstream != set(training_by_sha):
        raise ValueError("one or more training jobs have no evaluation consumer")
    if partition == "confirmation":
        candidates = [job for job in evaluation if job.get("all_on_reference") is False]
        references = [job for job in evaluation if job.get("all_on_reference") is True]
        if len(candidates) != 1 or len(references) != 1:
            raise ValueError("confirmation execution requires one K and one all-on job")
    return result


def validate_training_output(
    output: Mapping[str, Any], job: Mapping[str, Any]
) -> dict[str, Any]:
    """Validate a backend's fixed-mask checkpoint before inference can start."""

    result = dict(output)
    expected = dict(job)
    if result.get("schema_version") != TRAINING_OUTPUT_SCHEMA:
        raise ValueError("training backend output schema mismatch")
    if result.get("job_sha256") != expected.get("job_sha256"):
        raise ValueError("training output job SHA256 mismatch")
    for field in (
        "donor_partition",
        "candidate_count",
        "generator_count",
        "selected_local_generator_ids",
        "mask_sha256",
        "training_donor_names_sha256",
        "source_sha256",
    ):
        if result.get(field) != expected.get(field):
            raise ValueError(f"training output {field} differs from job")
    for field in (
        "fixed_mask_for_entire_run",
        "fixed_mask_persisted_in_checkpoint",
        "fresh_optimizer",
        "same_seed_schedule_across_k",
        "row_level_donor_allowlist_enforced",
    ):
        if result.get(field) is not True:
            raise ValueError(f"training output requires {field}=true")
    for field in ("official_validation_used", "official_test_used"):
        if result.get(field) is not False:
            raise ValueError(f"training output requires {field}=false")
    checkpoint = Path(str(result.get("checkpoint_path", ""))).resolve()
    if not checkpoint.is_file():
        raise ValueError("training output checkpoint does not exist")
    observed_checkpoint_sha = _sha256_file(checkpoint)
    if result.get("checkpoint_sha256") != observed_checkpoint_sha:
        raise ValueError("training output checkpoint SHA256 mismatch")
    supplied = result.get("output_sha256")
    expected_output_sha = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "output_sha256"}
    )
    if supplied != expected_output_sha:
        raise ValueError("training output payload SHA256 mismatch")
    return result


def _backend_command(backend: str | Path, arguments: Sequence[str]) -> list[str]:
    path = Path(backend).resolve()
    if not path.is_file():
        raise ValueError(f"backend does not exist: {path}")
    if path.suffix == ".py":
        return [sys.executable, str(path), *arguments]
    return [str(path), *arguments]


def _run_backend(command: Sequence[str], *, label: str) -> None:
    completed = subprocess.run(
        list(command),
        check=False,
        stdin=subprocess.DEVNULL,
        stdout=None,
        stderr=None,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{label} backend exited with code {completed.returncode}")


def dry_run_execution(
    plan: Mapping[str, Any],
    *,
    training_backend: str | Path,
    inference_backend: str | Path,
    work_dir: str | Path,
) -> dict[str, Any]:
    """Validate commands and output paths without claiming C10 or executing."""

    validated = validate_execution_plan(plan)
    train_backend = Path(training_backend).resolve()
    infer_backend = Path(inference_backend).resolve()
    if not train_backend.is_file() or not infer_backend.is_file():
        raise ValueError("training and inference backend files must exist")
    root = Path(work_dir).resolve()
    if root.exists():
        raise FileExistsError("executor work directory must not already exist")
    commands: list[dict[str, Any]] = []
    for job in validated["training_jobs"]:
        stem = job["job_sha256"]
        commands.append(
            {
                "stage": "training",
                "job_sha256": stem,
                "command": _backend_command(
                    train_backend,
                    (
                        "--job-json",
                        str(root / "jobs" / f"train_{stem}.json"),
                        "--out-json",
                        str(root / "outputs" / f"train_{stem}.json"),
                        "--work-dir",
                        str(root / "runs" / stem),
                    ),
                ),
            }
        )
    for job in validated["evaluation_jobs"]:
        stem = job["job_sha256"]
        commands.append(
            {
                "stage": "inference",
                "job_sha256": stem,
                "command": _backend_command(
                    infer_backend,
                    (
                        "--runtime-job-json",
                        str(root / "jobs" / f"eval_{stem}.json"),
                        "--out-json",
                        str(root / "outputs" / f"eval_{stem}.json"),
                    ),
                ),
            }
        )
    result: dict[str, Any] = {
        "schema_version": EXECUTION_MANIFEST_SCHEMA,
        "status": "dry_run_validated_not_executed",
        "plan_sha256": validated["plan_sha256"],
        "donor_partition": validated["donor_partition"],
        "c10_claimed": False,
        "canonical_release_api": (
            "build_canonical_generator_result_with_c10_ledger"
            if validated["donor_partition"] == "confirmation"
            else None
        ),
        "training_or_inference_started": False,
        "training_backend": {
            "path": str(train_backend),
            "sha256": _sha256_file(train_backend),
            "protocol": "job_json_to_fixed_mask_checkpoint_v1",
        },
        "inference_backend": {
            "path": str(infer_backend),
            "sha256": _sha256_file(infer_backend),
            "protocol": "runtime_job_json_to_exact_metrics_v1",
        },
        "commands": commands,
    }
    result["manifest_sha256"] = canonical_json_sha256(result)
    return result


def execute_plan(
    plan: Mapping[str, Any],
    *,
    training_backend: str | Path,
    inference_backend: str | Path,
    work_dir: str | Path,
    ledger_root: str | Path | None = None,
) -> dict[str, Any]:
    """Execute all fixed-mask training jobs and their bound inference jobs."""

    validated = validate_execution_plan(plan)
    preflight = dry_run_execution(
        validated,
        training_backend=training_backend,
        inference_backend=inference_backend,
        work_dir=work_dir,
    )
    if validated["donor_partition"] == "confirmation":
        if ledger_root is None:
            raise ValueError("confirmation execution requires ledger_root")
        if "c10_consumption" not in validated:
            raise ValueError(
                "confirmation plan must be atomically claimed before execution"
            )
        verify_claimed_exact_evaluation_plan(
            validated, ledger_root, require_active=True
        )
        begin_c10_execution_once(
            validated,
            ledger_root,
            execution_input_sha256=preflight["manifest_sha256"],
            executor_work_id=str(Path(work_dir).resolve()),
        )
    elif "c10_consumption" in validated:
        raise ValueError("ranking plan must not carry a C10 claim")

    root = Path(work_dir).resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / "jobs").mkdir()
    (root / "outputs").mkdir()
    (root / "runs").mkdir()
    training_outputs: dict[str, dict[str, Any]] = {}
    for job in validated["training_jobs"]:
        digest = str(job["job_sha256"])
        job_path = _write_json_exclusive(root / "jobs" / f"train_{digest}.json", job)
        output_path = root / "outputs" / f"train_{digest}.json"
        run_dir = root / "runs" / digest
        command = _backend_command(
            training_backend,
            (
                "--job-json",
                str(job_path),
                "--out-json",
                str(output_path),
                "--work-dir",
                str(run_dir),
            ),
        )
        if _sha256_file(training_backend) != preflight["training_backend"]["sha256"]:
            raise RuntimeError("training backend changed after preflight")
        _run_backend(command, label=f"training job {digest}")
        output = validate_training_output(
            _read_json(output_path, label="training backend output"), job
        )
        training_outputs[digest] = output

    evaluation_outputs: list[dict[str, Any]] = []
    for job in validated["evaluation_jobs"]:
        upstream_sha = str(job["upstream_training_job_sha256"])
        training_output = training_outputs[upstream_sha]
        runtime: dict[str, Any] = {
            "schema_version": RUNTIME_JOB_SCHEMA,
            "evaluation_job": job,
            "training_output": training_output,
            "evaluated_checkpoint_path": training_output["checkpoint_path"],
            "evaluated_checkpoint_sha256": training_output["checkpoint_sha256"],
            "upstream_training_manifest_sha256": training_output["output_sha256"],
        }
        runtime["runtime_job_sha256"] = canonical_json_sha256(runtime)
        digest = str(job["job_sha256"])
        runtime_path = _write_json_exclusive(
            root / "jobs" / f"eval_{digest}.json", runtime
        )
        output_path = root / "outputs" / f"eval_{digest}.json"
        command = _backend_command(
            inference_backend,
            (
                "--runtime-job-json",
                str(runtime_path),
                "--out-json",
                str(output_path),
            ),
        )
        if _sha256_file(inference_backend) != preflight["inference_backend"]["sha256"]:
            raise RuntimeError("inference backend changed after preflight")
        _run_backend(command, label=f"inference job {digest}")
        raw_output = _read_json(output_path, label="inference backend output")
        validated_output = validate_exact_evaluation_output(raw_output, job)
        if raw_output.get("evaluated_checkpoint_sha256") != training_output[
            "checkpoint_sha256"
        ]:
            raise ValueError("inference output used the wrong adapted checkpoint")
        if raw_output.get("upstream_training_manifest_sha256") != training_output[
            "output_sha256"
        ]:
            raise ValueError("inference output used the wrong training manifest")
        evaluation_outputs.append(
            {
                "job_sha256": digest,
                "output_path": str(output_path),
                "output_sha256": validated_output["output_sha256"],
            }
        )

    result: dict[str, Any] = {
        "schema_version": EXECUTION_MANIFEST_SCHEMA,
        "status": "fixed_mask_training_and_inference_completed",
        "plan_sha256": validated["plan_sha256"],
        "donor_partition": validated["donor_partition"],
        "canonical_mask_allowed": False,
        "training_backend_sha256": preflight["training_backend"]["sha256"],
        "inference_backend_sha256": preflight["inference_backend"]["sha256"],
        "training_outputs": [
            {
                "job_sha256": digest,
                "output_sha256": output["output_sha256"],
                "checkpoint_sha256": output["checkpoint_sha256"],
            }
            for digest, output in training_outputs.items()
        ],
        "evaluation_outputs": evaluation_outputs,
        "c10_claim_sha256": (
            validated.get("c10_consumption", {})
            .get("claim", {})
            .get("claim_sha256")
        ),
        "canonical_release_api": (
            "build_canonical_generator_result_with_c10_ledger"
            if validated["donor_partition"] == "confirmation"
            else None
        ),
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["manifest_sha256"] = canonical_json_sha256(result)
    _write_json_exclusive(root / "execution_manifest.json", result)
    return result


def _claim_for_execution(
    plan: Mapping[str, Any],
    *,
    ledger_root: str | Path,
    donor_split_path: str | Path,
    preselection_path: str | Path,
) -> dict[str, Any]:
    if plan.get("donor_partition") != "confirmation":
        return dict(plan)
    if "c10_consumption" in plan:
        verify_claimed_exact_evaluation_plan(plan, ledger_root, require_active=True)
        return dict(plan)
    return claim_exact_evaluation_plan(
        plan,
        ledger_root,
        donor_split_manifest=_read_json(donor_split_path, label="donor split"),
        preselection_manifest=_read_json(preselection_path, label="preselection"),
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Execute or dry-run a sealed generator exact-evaluation plan"
    )
    parser.add_argument("--plan", required=True)
    parser.add_argument("--training-backend", required=True)
    parser.add_argument("--inference-backend", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--ledger-root")
    parser.add_argument("--donor-split-manifest")
    parser.add_argument("--preselection-manifest")
    args = parser.parse_args(argv)
    out_manifest = Path(args.out_manifest).resolve()
    if out_manifest.exists():
        raise FileExistsError(f"refusing to overwrite output manifest: {out_manifest}")
    plan = _read_json(args.plan, label="evaluation plan")
    if not args.execute:
        result = dry_run_execution(
            plan,
            training_backend=args.training_backend,
            inference_backend=args.inference_backend,
            work_dir=args.work_dir,
        )
    else:
        # Exhaust every read-only/path-level failure before an irreversible C10
        # claim is published.  This preflight never creates ``work_dir`` and
        # never opens donor data.
        dry_run_execution(
            plan,
            training_backend=args.training_backend,
            inference_backend=args.inference_backend,
            work_dir=args.work_dir,
        )
        if plan.get("donor_partition") == "confirmation":
            if not all(
                (
                    args.ledger_root,
                    args.donor_split_manifest,
                    args.preselection_manifest,
                )
            ) and "c10_consumption" not in plan:
                raise ValueError(
                    "unclaimed confirmation execution requires ledger, split, "
                    "and preselection paths"
                )
            plan = _claim_for_execution(
                plan,
                ledger_root=args.ledger_root,
                donor_split_path=args.donor_split_manifest,
                preselection_path=args.preselection_manifest,
            )
        try:
            result = execute_plan(
                plan,
                training_backend=args.training_backend,
                inference_backend=args.inference_backend,
                work_dir=args.work_dir,
                ledger_root=args.ledger_root,
            )
        except Exception as exc:
            if plan.get("donor_partition") == "confirmation" and plan.get(
                "c10_consumption"
            ) is not None:
                failure_evidence = canonical_json_sha256(
                    {
                        "schema_version": "kmlee_bam.generator_exact_execution_failure.v1",
                        "plan_sha256": plan.get("plan_sha256"),
                        "error_type": type(exc).__name__,
                        "error_message": str(exc),
                    }
                )
                record_c10_confirmation_failure(
                    plan,
                    args.ledger_root,
                    failure_evidence_sha256=failure_evidence,
                    detail=(
                        "fixed-mask training/inference executor failed after "
                        f"C10 claim: {type(exc).__name__}: {exc}"
                    ),
                )
            raise
    _write_json_exclusive(out_manifest, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
