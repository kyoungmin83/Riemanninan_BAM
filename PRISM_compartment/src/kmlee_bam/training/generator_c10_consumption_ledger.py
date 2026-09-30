"""Best-effort append-only, one-shot C10 consumption ledger.

This module reserves a canonical generator-confirmation donor partition before
any C10 execution is allowed.  The partition identity is derived from the
*unordered* set of confirmation donor names, so merely reordering the same ten
donors cannot create a second attempt.  An immutable claim binds that partition
to the sealed donor-split manifest, the R10 preselection manifest, the selected
mask/K, and the exact constraint-registry hash.

The implementation uses a fully written temporary file followed by an atomic
same-filesystem hard link to publish each claim or terminal outcome without an
overwrite path.  A claim is consumption: even if execution crashes, is
abandoned, or fails, another claim for the same C10 partition is refused.
The only API state transitions are ``claimed -> completed`` and
``claimed -> failed``; both are terminal and each record seals plan/job plus
confirmation-input/canonical-output hashes when those artifacts exist.
Executors must also atomically publish the sole execution-start marker; a
read-only claim verification is deliberately insufficient to authorize a run.

This is a fail-closed operational control, not a cryptographic transparency
service.  A user able to copy, delete, replace, restore, or roll back the ledger
filesystem can erase history.  Read-only modes, hashes, and append-only API
semantics improve local evidence preservation but cannot prevent such rollback.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import socket
from typing import Any, Iterable, Mapping
import uuid


CLAIM_SCHEMA = "kmlee_bam.generator_c10_consumption_claim.v1"
OUTCOME_SCHEMA = "kmlee_bam.generator_c10_consumption_outcome.v1"
EXECUTION_START_SCHEMA = "kmlee_bam.generator_c10_execution_start.v1"
EVENT_SCHEMA = "kmlee_bam.generator_c10_consumption_event.v1"
PARTITION_SCHEMA = "kmlee_bam.generator_c10_partition_identity.v1"
PLAN_BINDING_SCHEMA = "kmlee_bam.generator_c10_claimed_plan.v1"

DONOR_SPLIT_SCHEMA = "kmlee_bam.generator_selection_donor_split.v1"
PRESELECTION_SCHEMA = "kmlee_bam.generator_K_preselection.v4"
EVALUATION_PLAN_SCHEMA = "kmlee_bam.generator_exact_evaluation_plan.v1"
CANONICAL_RESULT_SCHEMA = "kmlee_bam.learned_generator_result.v1"
EVALUATION_JOB_SCHEMA = "kmlee_bam.generator_exact_evaluation_job.v1"
TERMINAL_STATES = frozenset({"completed", "failed"})

FILESYSTEM_LIMITATION = (
    "Best-effort local-filesystem evidence only: filesystem copy, deletion, "
    "replacement, snapshot restore, and rollback cannot be cryptographically "
    "prevented without an external append-only witness."
)


class C10LedgerError(RuntimeError):
    """Base class for ledger publication and validation failures."""


class C10AlreadyConsumedError(C10LedgerError):
    """Raised whenever a C10 partition already has a claim."""


class C10LedgerIntegrityError(C10LedgerError):
    """Raised for corrupt, missing, or inconsistent ledger evidence."""


def canonical_json_sha256(value: Any) -> str:
    """Hash a JSON-compatible value with the project's canonical encoding."""

    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            dict(value),
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _require_sha256(value: Any, *, label: str) -> str:
    digest = str(value).lower()
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise ValueError(f"{label} must be a lowercase SHA256 digest")
    return digest


def _require_positive_int(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return int(value)


def _validate_sha256_fields(
    value: Mapping[str, Any],
    *,
    required_fields: Iterable[str],
    label: str,
) -> dict[str, str]:
    expected = tuple(required_fields)
    if set(value) != set(expected):
        raise ValueError(
            f"{label} fields must be exactly {sorted(expected)}"
        )
    return {
        field: _require_sha256(value[field], label=f"{label}.{field}")
        for field in expected
    }


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace(
        "+00:00", "Z"
    )


def _mask_sha256(selected_ids: Iterable[int], *, candidate_count: int) -> str:
    count = _require_positive_int(candidate_count, label="candidate_count")
    raw = tuple(selected_ids)
    if any(isinstance(value, bool) or not isinstance(value, int) for value in raw):
        raise ValueError("selected_local_generator_ids must contain integers")
    selected = tuple(sorted(int(value) for value in raw))
    if len(set(selected)) != len(selected) or any(
        value < 0 or value >= count for value in selected
    ):
        raise ValueError(
            "selected_local_generator_ids must be unique and within candidate_count"
        )
    mask = bytearray(count)
    for index in selected:
        mask[index] = 1
    return hashlib.sha256(bytes(mask)).hexdigest()


def _validated_manifest_sha256(
    manifest: Mapping[str, Any],
    *,
    field: str,
    label: str,
) -> str:
    supplied = _require_sha256(manifest.get(field), label=f"{label}.{field}")
    expected = canonical_json_sha256(
        {key: value for key, value in manifest.items() if key != field}
    )
    if supplied != expected:
        raise ValueError(f"{label} {field} mismatch")
    return supplied


def build_c10_consumption_binding(
    *,
    donor_split_manifest: Mapping[str, Any],
    preselection_manifest: Mapping[str, Any],
    constraint_registry_sha256: str,
) -> dict[str, Any]:
    """Validate and bind the three sealed inputs of one C10 attempt.

    The returned ``c10_partition_sha256`` deliberately hashes sorted donor
    names.  Donor order is separately sealed for paired evaluation, while the
    set hash prevents a reordered C10 list from evading one-shot consumption.
    """

    split = dict(donor_split_manifest)
    if split.get("schema_version") != DONOR_SPLIT_SCHEMA:
        raise ValueError("donor split manifest schema mismatch")
    for field, expected in (
        ("source_split", "train_only"),
        ("disjoint", True),
        ("complete", True),
        ("official_validation_used", False),
        ("official_test_used", False),
        ("train_donor_count", 64),
        ("weight_donor_count", 44),
        ("ranking_donor_count", 10),
        ("confirmation_donor_count", 10),
        ("hard_support_passed", True),
    ):
        if split.get(field) != expected:
            raise ValueError(f"donor split manifest requires {field}={expected!r}")
    split_sha256 = _validated_manifest_sha256(
        split,
        field="manifest_sha256",
        label="donor split manifest",
    )
    support = split.get("hard_support_audit")
    if not isinstance(support, Mapping):
        raise ValueError("donor split manifest lacks hard-support audit evidence")
    support = dict(support)
    if support.get("schema_version") != "kmlee_bam.generator_split_hard_support.v4":
        raise ValueError("donor split hard-support audit schema mismatch")
    if support.get("passed") is not True:
        raise ValueError("donor split hard-support audit did not pass")
    partition_support = support.get("partitions")
    if not isinstance(partition_support, Mapping) or set(partition_support) != {
        "W44",
        "R10",
        "C10",
    }:
        raise ValueError("donor split hard-support partitions are incomplete")
    if any(
        not isinstance(partition_support[name], Mapping)
        or partition_support[name].get("passed") is not True
        for name in ("W44", "R10", "C10")
    ):
        raise ValueError("a donor split partition failed hard support")
    supplied_support_sha256 = _require_sha256(
        support.get("audit_sha256"), label="hard_support_audit.audit_sha256"
    )
    expected_support_sha256 = canonical_json_sha256(
        {key: value for key, value in support.items() if key != "audit_sha256"}
    )
    if supplied_support_sha256 != expected_support_sha256:
        raise ValueError("donor split hard-support audit SHA256 mismatch")
    donor_names_raw = split.get("confirmation_donor_names")
    if not isinstance(donor_names_raw, list):
        raise ValueError("confirmation_donor_names must be a list")
    donor_names = tuple(str(value) for value in donor_names_raw)
    if len(donor_names) != 10 or any(not value for value in donor_names):
        raise ValueError("C10 must contain exactly ten non-empty donor names")
    if len(set(donor_names)) != len(donor_names):
        raise ValueError("C10 donor names must be unique")
    donor_order_sha256 = canonical_json_sha256(list(donor_names))
    if split.get("confirmation_donor_sha256") != donor_order_sha256:
        raise ValueError("confirmation donor-order SHA256 mismatch")
    all_partition_names: list[str] = []
    all_partition_ids: list[int] = []
    for prefix, expected_count in (
        ("weight", 44),
        ("ranking", 10),
        ("confirmation", 10),
    ):
        values = split.get(f"{prefix}_donor_names")
        if not isinstance(values, list) or len(values) != expected_count:
            raise ValueError(
                f"{prefix}_donor_names must contain exactly {expected_count} names"
            )
        parsed = [str(value) for value in values]
        if any(not value for value in parsed) or len(set(parsed)) != expected_count:
            raise ValueError(f"{prefix}_donor_names must be unique and non-empty")
        if split.get(f"{prefix}_donor_sha256") != canonical_json_sha256(parsed):
            raise ValueError(f"{prefix} donor-order SHA256 mismatch")
        all_partition_names.extend(parsed)
        ids = split.get(f"{prefix}_donor_ids")
        if not isinstance(ids, list) or len(ids) != expected_count or any(
            isinstance(value, bool) or not isinstance(value, int) for value in ids
        ):
            raise ValueError(
                f"{prefix}_donor_ids must contain exactly {expected_count} integers"
            )
        if len(set(ids)) != expected_count:
            raise ValueError(f"{prefix}_donor_ids must be unique")
        all_partition_ids.extend(ids)
    if len(set(all_partition_names)) != 64:
        raise ValueError("W44, R10, and C10 donor names must be mutually disjoint")
    if len(set(all_partition_ids)) != 64:
        raise ValueError("W44, R10, and C10 donor ids must be mutually disjoint")
    if split.get("train_donor_sha256") != canonical_json_sha256(
        all_partition_names
    ):
        raise ValueError("train donor-order SHA256 mismatch")
    partition_sha256 = canonical_json_sha256(
        {
            "schema_version": PARTITION_SCHEMA,
            "confirmation_donor_names_sorted": sorted(donor_names),
        }
    )

    preselection = dict(preselection_manifest)
    if preselection.get("schema_version") != PRESELECTION_SCHEMA:
        raise ValueError("preselection manifest schema mismatch")
    for field, expected in (
        ("selection_partition", "R10"),
        ("confirmation_partition_opened", False),
        ("prior_confirmation_attempts", 0),
        ("failure_policy", "abort_without_trying_another_K"),
        ("official_validation_used", False),
        ("official_test_used", False),
    ):
        if preselection.get(field) != expected:
            raise ValueError(f"preselection manifest requires {field}={expected!r}")
    preselection_sha256 = _validated_manifest_sha256(
        preselection,
        field="manifest_sha256",
        label="preselection manifest",
    )
    if preselection.get("donor_split_manifest_sha256") != split_sha256:
        raise ValueError("preselection donor-split SHA256 mismatch")
    ranking_result_sha256 = _require_sha256(
        preselection.get("ranking_result_sha256"),
        label="preselection.ranking_result_sha256",
    )
    candidate_count = _require_positive_int(
        preselection.get("candidate_count"), label="preselection.candidate_count"
    )
    preselected_k = _require_positive_int(
        preselection.get("preselected_k"), label="preselection.preselected_k"
    )
    if preselected_k >= candidate_count:
        raise ValueError("canonical C10 preselection must be non-all-on")
    selected_raw = preselection.get("selected_local_generator_ids")
    if not isinstance(selected_raw, list):
        raise ValueError("preselection selected_local_generator_ids must be a list")
    if len(selected_raw) != preselected_k:
        raise ValueError("preselection K differs from selected generator ids")
    selected = tuple(sorted(int(value) for value in selected_raw)) if all(
        isinstance(value, int) and not isinstance(value, bool) for value in selected_raw
    ) else ()
    if list(selected) != selected_raw:
        raise ValueError("preselection generator ids must be sorted unique integers")
    expected_mask_sha256 = _mask_sha256(selected, candidate_count=candidate_count)
    mask_sha256 = _require_sha256(
        preselection.get("mask_sha256"), label="preselection.mask_sha256"
    )
    if mask_sha256 != expected_mask_sha256:
        raise ValueError("preselection mask SHA256 mismatch")
    registry_sha256 = _require_sha256(
        constraint_registry_sha256, label="constraint_registry_sha256"
    )

    binding: dict[str, Any] = {
        "c10_partition_sha256": partition_sha256,
        "confirmation_donor_order_sha256": donor_order_sha256,
        "confirmation_donor_count": len(donor_names),
        "donor_split_manifest_sha256": split_sha256,
        "preselection_manifest_sha256": preselection_sha256,
        "constraint_registry_sha256": registry_sha256,
        "ranking_result_sha256": ranking_result_sha256,
        "candidate_count": candidate_count,
        "preselected_k": preselected_k,
        "selected_local_generator_ids": list(selected),
        "mask_sha256": mask_sha256,
    }
    binding["binding_sha256"] = canonical_json_sha256(binding)
    return binding


def _ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    if path.is_symlink() or not path.is_dir():
        raise C10LedgerIntegrityError(f"ledger path is not a real directory: {path}")


def _ledger_layout(root: str | Path) -> tuple[Path, Path, Path, Path, Path]:
    ledger_root = Path(root).absolute()
    _ensure_directory(ledger_root)
    claims = ledger_root / "claims"
    execution_starts = ledger_root / "execution_starts"
    outcomes = ledger_root / "outcomes"
    events = ledger_root / "events"
    for path in (claims, execution_starts, outcomes, events):
        _ensure_directory(path)
    return ledger_root, claims, execution_starts, outcomes, events


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_exclusive_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Publish a complete immutable file without an overwrite operation."""

    if path.parent.is_symlink() or not path.parent.is_dir():
        raise C10LedgerIntegrityError(
            f"ledger record parent is not a real directory: {path.parent}"
        )
    temporary = path.parent / f".{path.name}.tmp-{uuid.uuid4().hex}"
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        try:
            data = _json_bytes(payload)
            view = memoryview(data)
            while view:
                written = os.write(descriptor, view)
                if written <= 0:
                    raise OSError("short write while publishing ledger record")
                view = view[written:]
            os.fsync(descriptor)
            os.fchmod(descriptor, 0o444)
        finally:
            os.close(descriptor)
        # link(2) is atomic and, unlike replace(2), has no overwrite mode.
        os.link(temporary, path, follow_symlinks=False)
        _fsync_directory(path.parent)
    finally:
        try:
            temporary.unlink()
            _fsync_directory(path.parent)
        except FileNotFoundError:
            pass


def _read_json_record(
    path: Path,
    *,
    schema: str,
    digest_field: str,
) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise C10LedgerIntegrityError(f"ledger record is missing/not regular: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise C10LedgerIntegrityError(f"cannot read ledger record: {path}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != schema:
        raise C10LedgerIntegrityError(f"ledger record schema mismatch: {path}")
    supplied = payload.get(digest_field)
    expected = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != digest_field}
    )
    if supplied != expected:
        raise C10LedgerIntegrityError(f"ledger record SHA256 mismatch: {path}")
    return payload


def _claim_path(claims: Path, partition_sha256: str) -> Path:
    return claims / f"{_require_sha256(partition_sha256, label='c10_partition_sha256')}.json"


def _outcome_path(outcomes: Path, partition_sha256: str) -> Path:
    return outcomes / f"{_require_sha256(partition_sha256, label='c10_partition_sha256')}.json"


def _execution_start_path(
    execution_starts: Path,
    partition_sha256: str,
) -> Path:
    return execution_starts / (
        f"{_require_sha256(partition_sha256, label='c10_partition_sha256')}.json"
    )


def _append_rejection_event(
    events_root: Path,
    *,
    existing_claim: Mapping[str, Any],
    proposed_binding: Mapping[str, Any],
) -> None:
    """Best-effort evidence for refused duplicate claims; refusal never depends on it."""

    try:
        partition = str(existing_claim["binding"]["c10_partition_sha256"])
        directory = events_root / partition
        _ensure_directory(directory)
        event: dict[str, Any] = {
            "schema_version": EVENT_SCHEMA,
            "event_type": "duplicate_claim_refused",
            "recorded_at_utc": _utc_now(),
            "writer_pid": os.getpid(),
            "writer_host": socket.gethostname(),
            "c10_partition_sha256": partition,
            "existing_claim_sha256": existing_claim["claim_sha256"],
            "existing_binding_sha256": existing_claim["binding"]["binding_sha256"],
            "proposed_binding_sha256": proposed_binding["binding_sha256"],
            "proposed_preselected_k": proposed_binding["preselected_k"],
            "filesystem_limitation": FILESYSTEM_LIMITATION,
        }
        event["event_sha256"] = canonical_json_sha256(event)
        name = (
            f"{event['recorded_at_utc'].replace(':', '').replace('-', '')}-"
            f"{uuid.uuid4().hex}-{event['event_sha256'][:12]}.json"
        )
        _publish_exclusive_json(directory / name, event)
    except (OSError, ValueError, C10LedgerError, KeyError, TypeError):
        # The claim itself is authoritative.  A failure in advisory rejection
        # logging must never turn a refusal into a permitted second attempt.
        return


def _append_execution_rejection_event(
    events_root: Path,
    *,
    claim: Mapping[str, Any],
    existing_start: Mapping[str, Any],
    execution_input_sha256: str,
    executor_work_id: str,
) -> None:
    """Best-effort append-only evidence for a refused backend rerun."""

    try:
        partition = str(claim["binding"]["c10_partition_sha256"])
        directory = events_root / partition
        _ensure_directory(directory)
        event: dict[str, Any] = {
            "schema_version": EVENT_SCHEMA,
            "event_type": "duplicate_execution_start_refused",
            "recorded_at_utc": _utc_now(),
            "writer_pid": os.getpid(),
            "writer_host": socket.gethostname(),
            "c10_partition_sha256": partition,
            "claim_sha256": claim["claim_sha256"],
            "existing_execution_start_sha256": existing_start[
                "execution_start_sha256"
            ],
            "proposed_execution_input_sha256": execution_input_sha256,
            "proposed_executor_work_id": executor_work_id,
            "rerun_refused": True,
            "filesystem_limitation": FILESYSTEM_LIMITATION,
        }
        event["event_sha256"] = canonical_json_sha256(event)
        name = (
            f"{event['recorded_at_utc'].replace(':', '').replace('-', '')}-"
            f"{uuid.uuid4().hex}-{event['event_sha256'][:12]}.json"
        )
        _publish_exclusive_json(directory / name, event)
    except (OSError, ValueError, C10LedgerError, KeyError, TypeError):
        return


def claim_c10_confirmation(
    ledger_root: str | Path,
    *,
    donor_split_manifest: Mapping[str, Any],
    preselection_manifest: Mapping[str, Any],
    constraint_registry_sha256: str,
    plan_intent_sha256: str | None = None,
    plan_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Atomically consume C10 for exactly one sealed preselected candidate.

    There is intentionally no idempotent success case.  Calling this function a
    second time is another attempt and is refused, even when every input is
    byte-for-byte identical to the original claim.
    """

    binding = build_c10_consumption_binding(
        donor_split_manifest=donor_split_manifest,
        preselection_manifest=preselection_manifest,
        constraint_registry_sha256=constraint_registry_sha256,
    )
    intent_sha256 = (
        None
        if plan_intent_sha256 is None
        else _require_sha256(plan_intent_sha256, label="plan_intent_sha256")
    )
    sealed_plan_evidence: dict[str, str] | None = None
    if plan_evidence is not None:
        sealed_plan_evidence = _validate_sha256_fields(
            plan_evidence,
            required_fields=(
                "plan_sha256",
                "candidate_training_job_sha256",
                "all_on_training_job_sha256",
                "candidate_evaluation_job_sha256",
                "all_on_evaluation_job_sha256",
            ),
            label="plan_evidence",
        )
        if intent_sha256 != sealed_plan_evidence["plan_sha256"]:
            raise ValueError("plan evidence differs from plan_intent_sha256")
    _, claims, _, _, events = _ledger_layout(ledger_root)
    target = _claim_path(claims, binding["c10_partition_sha256"])
    claim: dict[str, Any] = {
        "schema_version": CLAIM_SCHEMA,
        "record_type": "claim",
        "state": "claimed",
        "claimed_at_utc": _utc_now(),
        "attempt_id": str(uuid.uuid4()),
        "writer_pid": os.getpid(),
        "writer_host": socket.gethostname(),
        "binding": binding,
        "plan_sha256": intent_sha256,
        "plan_intent_sha256": intent_sha256,
        "plan_evidence": sealed_plan_evidence,
        "one_shot_policy": {
            "claim_consumes_partition": True,
            "retry_after_claim_forbidden": True,
            "different_k_after_claim_forbidden": True,
            "atomic_execution_start_marker_required": True,
            "same_claim_backend_rerun_forbidden": True,
            "retry_after_failed_forbidden": True,
            "retry_after_completed_forbidden": True,
            "allowed_state_transitions": ["claimed->completed", "claimed->failed"],
        },
        "publication": "complete-temp-file_then_atomic_hard-link_no-overwrite",
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    claim["claim_sha256"] = canonical_json_sha256(claim)
    try:
        _publish_exclusive_json(target, claim)
    except FileExistsError:
        existing = _read_json_record(
            target, schema=CLAIM_SCHEMA, digest_field="claim_sha256"
        )
        _append_rejection_event(
            events,
            existing_claim=existing,
            proposed_binding=binding,
        )
        existing_binding = existing.get("binding", {})
        existing_k = existing_binding.get("preselected_k", "unknown")
        existing_sha = existing.get("claim_sha256", "unknown")
        raise C10AlreadyConsumedError(
            "C10 partition is already consumed by claim "
            f"{existing_sha} (preselected K={existing_k}); no second attempt is allowed"
        ) from None
    return claim


def read_c10_consumption(
    ledger_root: str | Path,
    *,
    c10_partition_sha256: str,
) -> dict[str, Any]:
    """Read and integrity-check the claim plus its optional terminal outcome."""

    _, claims, execution_starts, outcomes, _ = _ledger_layout(ledger_root)
    partition = _require_sha256(
        c10_partition_sha256, label="c10_partition_sha256"
    )
    claim = _read_json_record(
        _claim_path(claims, partition),
        schema=CLAIM_SCHEMA,
        digest_field="claim_sha256",
    )
    binding = claim.get("binding")
    if not isinstance(binding, dict):
        raise C10LedgerIntegrityError("claim binding is missing")
    binding_sha256 = binding.get("binding_sha256")
    if binding_sha256 != canonical_json_sha256(
        {key: value for key, value in binding.items() if key != "binding_sha256"}
    ):
        raise C10LedgerIntegrityError("claim binding SHA256 mismatch")
    if binding.get("c10_partition_sha256") != partition:
        raise C10LedgerIntegrityError("claim is stored under the wrong C10 partition")

    execution_target = _execution_start_path(execution_starts, partition)
    execution_start: dict[str, Any] | None = None
    if execution_target.exists() or execution_target.is_symlink():
        execution_start = _read_json_record(
            execution_target,
            schema=EXECUTION_START_SCHEMA,
            digest_field="execution_start_sha256",
        )
        if execution_start.get("c10_partition_sha256") != partition:
            raise C10LedgerIntegrityError(
                "execution start is stored under the wrong partition"
            )
        if execution_start.get("claim_sha256") != claim.get("claim_sha256"):
            raise C10LedgerIntegrityError(
                "execution start does not extend the immutable claim"
            )
        if execution_start.get("binding_sha256") != binding.get("binding_sha256"):
            raise C10LedgerIntegrityError("execution start binding SHA256 mismatch")
        if execution_start.get("plan_sha256") != claim.get("plan_sha256"):
            raise C10LedgerIntegrityError("execution start plan SHA256 mismatch")
        if execution_start.get("plan_evidence") != claim.get("plan_evidence"):
            raise C10LedgerIntegrityError("execution start plan/job evidence mismatch")
        if execution_start.get("previous_record_sha256") != claim.get(
            "claim_sha256"
        ):
            raise C10LedgerIntegrityError("execution start hash chain is broken")

    outcome_target = _outcome_path(outcomes, partition)
    outcome: dict[str, Any] | None = None
    if outcome_target.exists() or outcome_target.is_symlink():
        outcome = _read_json_record(
            outcome_target,
            schema=OUTCOME_SCHEMA,
            digest_field="outcome_sha256",
        )
        if outcome.get("c10_partition_sha256") != partition:
            raise C10LedgerIntegrityError("outcome is stored under the wrong partition")
        if outcome.get("claim_sha256") != claim.get("claim_sha256"):
            raise C10LedgerIntegrityError("outcome does not extend the immutable claim")
        if outcome.get("binding_sha256") != binding.get("binding_sha256"):
            raise C10LedgerIntegrityError("outcome binding SHA256 mismatch")
        expected_execution_sha256 = (
            None
            if execution_start is None
            else execution_start["execution_start_sha256"]
        )
        if outcome.get("execution_start_sha256") != expected_execution_sha256:
            raise C10LedgerIntegrityError(
                "outcome execution-start evidence differs from the ledger"
            )
        expected_previous_sha256 = (
            claim.get("claim_sha256")
            if execution_start is None
            else expected_execution_sha256
        )
        if outcome.get("previous_record_sha256") != expected_previous_sha256:
            raise C10LedgerIntegrityError("outcome hash chain is broken")
        if outcome.get("state") not in TERMINAL_STATES or outcome.get(
            "outcome"
        ) != outcome.get("state"):
            raise C10LedgerIntegrityError("outcome has an invalid terminal state")
        if outcome["state"] == "completed":
            if execution_start is None:
                raise C10LedgerIntegrityError(
                    "completed outcome has no one-shot execution-start marker"
                )
            try:
                input_sha256 = _require_sha256(
                    outcome.get("confirmation_input_sha256"),
                    label="outcome.confirmation_input_sha256",
                )
                output_sha256 = _require_sha256(
                    outcome.get("canonical_output_sha256"),
                    label="outcome.canonical_output_sha256",
                )
            except ValueError as exc:
                raise C10LedgerIntegrityError(
                    "completed outcome lacks input/output SHA256 evidence"
                ) from exc
            if outcome.get("evidence_sha256") != output_sha256 or not input_sha256:
                raise C10LedgerIntegrityError(
                    "completed outcome output/evidence SHA256 mismatch"
                )
        elif outcome.get("canonical_output_sha256") is not None:
            raise C10LedgerIntegrityError(
                "failed outcome cannot carry canonical output evidence"
            )
    state = "claimed" if outcome is None else str(outcome["state"])
    return {
        "state": state,
        "execution_started": execution_start is not None,
        "claim": claim,
        "execution_start": execution_start,
        "outcome": outcome,
    }


def verify_c10_claim(
    ledger_root: str | Path,
    claim: Mapping[str, Any],
    *,
    require_active: bool = True,
) -> dict[str, Any]:
    """Verify immutable evidence; this read-only check never authorizes execution.

    ``require_active`` means only "no terminal outcome exists."  Executors must
    separately win :func:`begin_c10_execution_once` before touching C10.
    """

    supplied = dict(claim)
    binding = supplied.get("binding")
    if not isinstance(binding, Mapping):
        raise C10LedgerIntegrityError("supplied claim lacks a binding")
    partition = _require_sha256(
        binding.get("c10_partition_sha256"), label="c10_partition_sha256"
    )
    entry = read_c10_consumption(
        ledger_root, c10_partition_sha256=partition
    )
    if entry["claim"] != supplied:
        raise C10LedgerIntegrityError("supplied claim differs from on-disk claim")
    if require_active and entry["outcome"] is not None:
        raise C10AlreadyConsumedError(
            f"C10 claim already has terminal outcome {entry['state']!r}"
        )
    return entry


def record_c10_outcome(
    ledger_root: str | Path,
    claim: Mapping[str, Any],
    *,
    outcome: str,
    evidence_sha256: str,
    detail: str,
    confirmation_input_sha256: str | None = None,
    canonical_output_sha256: str | None = None,
) -> dict[str, Any]:
    """Append the sole terminal ``failed`` or ``completed`` record for a claim."""

    if outcome not in TERMINAL_STATES:
        raise ValueError("outcome must be 'failed' or 'completed'")
    if not isinstance(detail, str) or not detail.strip():
        raise ValueError("outcome detail must be a non-empty string")
    evidence = _require_sha256(evidence_sha256, label="evidence_sha256")
    input_sha256 = (
        None
        if confirmation_input_sha256 is None
        else _require_sha256(
            confirmation_input_sha256, label="confirmation_input_sha256"
        )
    )
    output_sha256 = (
        None
        if canonical_output_sha256 is None
        else _require_sha256(
            canonical_output_sha256, label="canonical_output_sha256"
        )
    )
    if outcome == "completed" and (
        input_sha256 is None or output_sha256 is None
    ):
        raise ValueError(
            "completed outcome requires confirmation_input_sha256 and "
            "canonical_output_sha256"
        )
    if outcome == "failed" and output_sha256 is not None:
        raise ValueError("failed outcome cannot carry a canonical_output_sha256")
    if output_sha256 is not None and evidence != output_sha256:
        raise ValueError("completed evidence_sha256 must equal canonical_output_sha256")
    entry = verify_c10_claim(ledger_root, claim, require_active=True)
    immutable_claim = entry["claim"]
    binding = immutable_claim["binding"]
    execution_start = entry["execution_start"]
    if outcome == "completed" and execution_start is None:
        raise C10LedgerIntegrityError(
            "completed outcome requires the durable one-shot execution-start marker"
        )
    previous_record_sha256 = (
        immutable_claim["claim_sha256"]
        if execution_start is None
        else execution_start["execution_start_sha256"]
    )
    _, _, _, outcomes, _ = _ledger_layout(ledger_root)
    target = _outcome_path(outcomes, binding["c10_partition_sha256"])
    result: dict[str, Any] = {
        "schema_version": OUTCOME_SCHEMA,
        "record_type": "terminal_outcome",
        "state": outcome,
        "outcome": outcome,
        "recorded_at_utc": _utc_now(),
        "writer_pid": os.getpid(),
        "writer_host": socket.gethostname(),
        "c10_partition_sha256": binding["c10_partition_sha256"],
        "binding_sha256": binding["binding_sha256"],
        "claim_sha256": immutable_claim["claim_sha256"],
        "execution_start_sha256": (
            None
            if execution_start is None
            else execution_start["execution_start_sha256"]
        ),
        "previous_record_sha256": previous_record_sha256,
        "evidence_sha256": evidence,
        "confirmation_input_sha256": input_sha256,
        "canonical_output_sha256": output_sha256,
        "detail": detail.strip(),
        "retry_forbidden": True,
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    result["outcome_sha256"] = canonical_json_sha256(result)
    try:
        _publish_exclusive_json(target, result)
    except FileExistsError:
        existing = _read_json_record(
            target, schema=OUTCOME_SCHEMA, digest_field="outcome_sha256"
        )
        raise C10AlreadyConsumedError(
            "C10 claim already has terminal outcome "
            f"{existing.get('outcome')!r}; terminal evidence cannot be replaced"
        ) from None
    return result


def _validated_plan_hash(plan: Mapping[str, Any]) -> str:
    supplied = _require_sha256(plan.get("plan_sha256"), label="plan.plan_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in plan.items() if key != "plan_sha256"}
    )
    if supplied != expected:
        raise ValueError("evaluation plan SHA256 mismatch")
    return supplied


def _validated_job_sha256(job: Mapping[str, Any], *, label: str) -> str:
    if job.get("schema_version") != EVALUATION_JOB_SCHEMA:
        raise ValueError(f"{label} schema mismatch")
    supplied = _require_sha256(job.get("job_sha256"), label=f"{label}.job_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in job.items() if key != "job_sha256"}
    )
    if supplied != expected:
        raise ValueError(f"{label} job SHA256 mismatch")
    return supplied


def _one_candidate_reference_pair(
    jobs: Any,
    *,
    label: str,
) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    if not isinstance(jobs, list) or len(jobs) != 2 or any(
        not isinstance(job, Mapping) for job in jobs
    ):
        raise ValueError(f"C10 plan must contain exactly two {label} jobs")
    candidates = [job for job in jobs if job.get("all_on_reference") is False]
    references = [job for job in jobs if job.get("all_on_reference") is True]
    if len(candidates) != 1 or len(references) != 1:
        raise ValueError(
            f"C10 plan must contain exactly one preselected K and one all-on {label} job"
        )
    return candidates[0], references[0]


def _validate_exact_job_pair(
    plan: Mapping[str, Any],
    binding: Mapping[str, Any],
) -> dict[str, str]:
    training_candidate, training_reference = _one_candidate_reference_pair(
        plan.get("training_jobs"), label="training"
    )
    evaluation_candidate, evaluation_reference = _one_candidate_reference_pair(
        plan.get("evaluation_jobs"), label="evaluation"
    )
    jobs = (
        (training_candidate, "candidate training", "paired_fixed_mask_retrain_w_plus_r54"),
        (training_reference, "all-on training", "paired_fixed_mask_retrain_w_plus_r54"),
        (evaluation_candidate, "candidate evaluation", "one_shot_confirmation_evaluation_c10"),
        (evaluation_reference, "all-on evaluation", "one_shot_confirmation_evaluation_c10"),
    )
    job_hashes: dict[str, str] = {}
    for job, label, expected_kind in jobs:
        if job.get("donor_partition") != "confirmation":
            raise ValueError(f"{label} job is not on the confirmation partition")
        if job.get("job_kind") != expected_kind:
            raise ValueError(f"{label} job kind mismatch")
        for field, expected in (
            ("source_split", "train_only"),
            ("official_validation_used", False),
            ("official_test_used", False),
        ):
            if job.get(field) != expected:
                raise ValueError(f"{label} requires {field}={expected!r}")
        job_hashes[label] = _validated_job_sha256(job, label=label)

    for candidate, label in (
        (training_candidate, "candidate training"),
        (evaluation_candidate, "candidate evaluation"),
    ):
        for field in (
            "candidate_count",
            "selected_local_generator_ids",
            "mask_sha256",
        ):
            if candidate.get(field) != binding[field]:
                raise ValueError(f"{label} {field} differs from preselection")
        if candidate.get("generator_count") != binding["preselected_k"]:
            raise ValueError(f"{label} K differs from preselection")

    expected_all_on_ids = list(range(int(binding["candidate_count"])))
    expected_all_on_mask = _mask_sha256(
        expected_all_on_ids, candidate_count=int(binding["candidate_count"])
    )
    for reference, label in (
        (training_reference, "all-on training"),
        (evaluation_reference, "all-on evaluation"),
    ):
        expected = {
            "candidate_count": binding["candidate_count"],
            "generator_count": binding["candidate_count"],
            "selected_local_generator_ids": expected_all_on_ids,
            "mask_sha256": expected_all_on_mask,
        }
        for field, value in expected.items():
            if reference.get(field) != value:
                raise ValueError(f"{label} {field} is not the exact all-on reference")

    for evaluation, label in (
        (evaluation_candidate, "candidate evaluation"),
        (evaluation_reference, "all-on evaluation"),
    ):
        if evaluation.get("donor_names_sha256") != binding[
            "confirmation_donor_order_sha256"
        ]:
            raise ValueError(f"{label} donor order differs from the sealed split")
        additive_contract = evaluation.get("additive_confirmation_contract")
        if not isinstance(additive_contract, Mapping) or additive_contract.get(
            "constraint_registry_sha256"
        ) != binding["constraint_registry_sha256"]:
            raise ValueError(f"{label} constraint registry differs from the binding")

    split = _donor_split_from_plan(plan)
    expected_training_names = list(split["weight_donor_names"]) + list(
        split["ranking_donor_names"]
    )
    expected_confirmation_names = list(split["confirmation_donor_names"])
    for training, label in (
        (training_candidate, "candidate training"),
        (training_reference, "all-on training"),
    ):
        if training.get("training_donor_names") != expected_training_names:
            raise ValueError(f"{label} donors are not exactly sealed W44+R10")
        if training.get("training_donor_names_sha256") != canonical_json_sha256(
            expected_training_names
        ):
            raise ValueError(f"{label} donor-name SHA256 mismatch")
    for evaluation, label in (
        (evaluation_candidate, "candidate evaluation"),
        (evaluation_reference, "all-on evaluation"),
    ):
        if evaluation.get("donor_names") != expected_confirmation_names:
            raise ValueError(f"{label} donors are not exactly sealed C10")

    if evaluation_candidate.get("upstream_training_job_sha256") != job_hashes[
        "candidate training"
    ]:
        raise ValueError("candidate evaluation is not bound to candidate training")
    if evaluation_reference.get("upstream_training_job_sha256") != job_hashes[
        "all-on training"
    ]:
        raise ValueError("all-on evaluation is not bound to all-on training")
    source_values = [job.get("source_sha256") for job, _, _ in jobs]
    if any(not isinstance(value, Mapping) or not value for value in source_values):
        raise ValueError("C10 jobs must embed non-empty source SHA256 mappings")
    for value in source_values:
        for name, digest in value.items():
            if not isinstance(name, str) or not name:
                raise ValueError("C10 source SHA256 names must be non-empty strings")
            _require_sha256(digest, label=f"source_sha256.{name}")
    if len({canonical_json_sha256(value) for value in source_values}) != 1:
        raise ValueError("C10 K/all-on jobs differ in source_sha256")
    if any(
        job.get("paired_retrain_group") != "frozen_final_k_vs_all_on_w_plus_r54"
        for job, _, _ in jobs
    ):
        raise ValueError("C10 K/all-on jobs have an invalid paired retrain group")

    return {
        "plan_sha256": _require_sha256(
            plan.get("plan_sha256"), label="plan.plan_sha256"
        ),
        "candidate_training_job_sha256": job_hashes["candidate training"],
        "all_on_training_job_sha256": job_hashes["all-on training"],
        "candidate_evaluation_job_sha256": job_hashes["candidate evaluation"],
        "all_on_evaluation_job_sha256": job_hashes["all-on evaluation"],
    }


def claim_exact_evaluation_plan(
    plan: Mapping[str, Any],
    ledger_root: str | Path,
    *,
    donor_split_manifest: Mapping[str, Any],
    preselection_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    """Claim C10 and return a hash-sealed exact-evaluation plan copy.

    This is the narrow integration hook for ``generator_exact_evaluation``.  It
    accepts only a validated confirmation plan and does not execute its training
    or evaluation jobs.
    """

    payload = dict(plan)
    if payload.get("schema_version") != EVALUATION_PLAN_SCHEMA:
        raise ValueError("evaluation plan schema mismatch")
    if payload.get("donor_partition") != "confirmation":
        raise ValueError("only a C10 confirmation plan can consume C10")
    for field, expected in (
        ("source_split", "train_only"),
        ("official_validation_used", False),
        ("official_test_used", False),
        ("canonical_mask_allowed", False),
    ):
        if payload.get(field) != expected:
            raise ValueError(f"C10 evaluation plan requires {field}={expected!r}")
    if "c10_consumption" in payload:
        raise C10AlreadyConsumedError("evaluation plan already carries a C10 claim")
    _validated_plan_hash(payload)
    embedded_split = payload.get("donor_split_manifest")
    if embedded_split is not None and dict(embedded_split) != dict(
        donor_split_manifest
    ):
        raise ValueError("evaluation plan embeds a different donor split manifest")
    payload["donor_split_manifest"] = dict(donor_split_manifest)
    payload["plan_sha256"] = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != "plan_sha256"}
    )
    plan_intent_sha256 = payload["plan_sha256"]
    constraint_sha256 = _require_sha256(
        payload.get("exact_constraint_registry_sha256"),
        label="plan.exact_constraint_registry_sha256",
    )
    binding = build_c10_consumption_binding(
        donor_split_manifest=donor_split_manifest,
        preselection_manifest=preselection_manifest,
        constraint_registry_sha256=constraint_sha256,
    )
    if payload.get("donor_split_manifest_sha256") != binding[
        "donor_split_manifest_sha256"
    ]:
        raise ValueError("evaluation plan donor-split SHA256 differs from binding")
    if payload.get("requested_generator_counts") != [
        binding["preselected_k"],
        binding["candidate_count"],
    ]:
        raise ValueError("C10 plan must request only preselected K and all-on")
    selection_policy = payload.get("selection_policy")
    if not isinstance(selection_policy, Mapping) or any(
        selection_policy.get(field) != expected
        for field, expected in (
            ("stage", "C10_one_shot_confirmation"),
            ("next_k_after_failure_forbidden", True),
            ("failure_is_terminal_no_canonical_mask", True),
        )
    ):
        raise ValueError("C10 evaluation plan one-shot selection policy is invalid")
    plan_evidence = _validate_exact_job_pair(payload, binding)

    claim = claim_c10_confirmation(
        ledger_root,
        donor_split_manifest=donor_split_manifest,
        preselection_manifest=preselection_manifest,
        constraint_registry_sha256=constraint_sha256,
        plan_intent_sha256=plan_intent_sha256,
        plan_evidence=plan_evidence,
    )
    result = dict(payload)
    result["unclaimed_plan_sha256"] = plan_intent_sha256
    result["preselection_manifest"] = dict(preselection_manifest)
    result["c10_consumption"] = {
        "schema_version": PLAN_BINDING_SCHEMA,
        "claim_required_before_any_c10_access": True,
        "claim_is_consumption": True,
        "atomic_execution_start_marker_required": True,
        "read_only_claim_verification_never_authorizes_execution": True,
        "claim": claim,
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    result["plan_sha256"] = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "plan_sha256"}
    )
    return result


def verify_claimed_exact_evaluation_plan(
    plan: Mapping[str, Any],
    ledger_root: str | Path,
    *,
    require_active: bool = True,
) -> dict[str, Any]:
    """Validate a claimed plan; this read-only check never authorizes execution."""

    payload = dict(plan)
    if payload.get("schema_version") != EVALUATION_PLAN_SCHEMA:
        raise C10LedgerIntegrityError("claimed evaluation plan schema mismatch")
    try:
        _validated_plan_hash(payload)
    except ValueError as exc:
        raise C10LedgerIntegrityError(str(exc)) from exc
    consumption = payload.get("c10_consumption")
    if not isinstance(consumption, Mapping) or consumption.get(
        "schema_version"
    ) != PLAN_BINDING_SCHEMA:
        raise C10LedgerIntegrityError("evaluation plan lacks a C10 claim binding")
    claim = consumption.get("claim")
    if not isinstance(claim, Mapping):
        raise C10LedgerIntegrityError("evaluation plan lacks the immutable C10 claim")
    entry = verify_c10_claim(ledger_root, claim, require_active=require_active)
    intent_sha256 = payload.get("unclaimed_plan_sha256")
    if claim.get("plan_intent_sha256") != intent_sha256:
        raise C10LedgerIntegrityError("claim does not bind the unclaimed plan hash")
    unclaimed = {
        key: value
        for key, value in payload.items()
        if key not in {"c10_consumption", "preselection_manifest", "unclaimed_plan_sha256"}
    }
    unclaimed["plan_sha256"] = intent_sha256
    try:
        _validated_plan_hash(unclaimed)
    except ValueError as exc:
        raise C10LedgerIntegrityError("unclaimed plan intent cannot be reconstructed") from exc
    binding = claim["binding"]
    preselection = payload.get("preselection_manifest")
    if not isinstance(preselection, Mapping):
        raise C10LedgerIntegrityError("claimed plan lacks its preselection manifest")
    try:
        rebuilt = build_c10_consumption_binding(
            donor_split_manifest=_donor_split_from_plan(payload),
            preselection_manifest=preselection,
            constraint_registry_sha256=payload.get(
                "exact_constraint_registry_sha256"
            ),
        )
    except ValueError as exc:
        raise C10LedgerIntegrityError(f"claimed plan binding is invalid: {exc}") from exc
    if rebuilt != binding:
        raise C10LedgerIntegrityError("claimed plan inputs differ from the on-disk binding")
    try:
        plan_evidence = _validate_exact_job_pair(unclaimed, binding)
    except ValueError as exc:
        raise C10LedgerIntegrityError(f"claimed plan job pair is invalid: {exc}") from exc
    if claim.get("plan_evidence") != plan_evidence:
        raise C10LedgerIntegrityError("claimed plan/job SHA evidence differs from claim")
    return entry


def begin_c10_execution_once(
    claimed_plan: Mapping[str, Any],
    ledger_root: str | Path,
    *,
    execution_input_sha256: str,
    executor_work_id: str,
) -> dict[str, Any]:
    """Atomically publish the sole permitted execution-start marker.

    Merely verifying a claim is read-only and never authorizes a backend run.
    Every executor must successfully call this function immediately before its
    first C10-affecting action.  A crash after publication remains consumed, and
    every subsequent call is refused even for the same plan, work id, or input.
    """

    if not isinstance(executor_work_id, str) or not executor_work_id.strip():
        raise ValueError("executor_work_id must be a non-empty string")
    input_sha256 = _require_sha256(
        execution_input_sha256, label="execution_input_sha256"
    )
    entry = verify_claimed_exact_evaluation_plan(
        claimed_plan, ledger_root, require_active=True
    )
    claim = entry["claim"]
    partition = claim["binding"]["c10_partition_sha256"]
    _, _, execution_starts, _, events = _ledger_layout(ledger_root)
    target = _execution_start_path(execution_starts, partition)
    marker: dict[str, Any] = {
        "schema_version": EXECUTION_START_SCHEMA,
        "record_type": "one_shot_execution_start",
        "execution_started": True,
        "started_at_utc": _utc_now(),
        "writer_pid": os.getpid(),
        "writer_host": socket.gethostname(),
        "executor_work_id": executor_work_id.strip(),
        "execution_input_sha256": input_sha256,
        "c10_partition_sha256": partition,
        "claim_sha256": claim["claim_sha256"],
        "binding_sha256": claim["binding"]["binding_sha256"],
        "plan_sha256": claim["plan_sha256"],
        "plan_evidence": claim["plan_evidence"],
        "previous_record_sha256": claim["claim_sha256"],
        "second_execution_forbidden": True,
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    marker["execution_start_sha256"] = canonical_json_sha256(marker)
    try:
        _publish_exclusive_json(target, marker)
    except FileExistsError:
        existing = _read_json_record(
            target,
            schema=EXECUTION_START_SCHEMA,
            digest_field="execution_start_sha256",
        )
        _append_execution_rejection_event(
            events,
            claim=claim,
            existing_start=existing,
            execution_input_sha256=input_sha256,
            executor_work_id=executor_work_id.strip(),
        )
        raise C10AlreadyConsumedError(
            "C10 execution was already started by marker "
            f"{existing.get('execution_start_sha256')}; backend rerun is forbidden"
        ) from None
    return marker


def _donor_split_from_plan(plan: Mapping[str, Any]) -> Mapping[str, Any]:
    split = plan.get("donor_split_manifest")
    if not isinstance(split, Mapping):
        raise ValueError(
            "claimed plan must embed donor_split_manifest for offline ledger verification"
        )
    return split


def attach_c10_canonical_release_evidence(
    result: Mapping[str, Any],
    claimed_plan: Mapping[str, Any],
    ledger_root: str | Path,
) -> dict[str, Any]:
    """Verify the active claim, terminally record success, and bind the result.

    This is the narrow canonical-release integration hook.  The input must be a
    result already validated by ``build_canonical_generator_result``.  No model
    training or evaluation is performed here.  A successful call terminally
    records ``completed``; it is not a preview/dry-run API.
    """

    entry = verify_claimed_exact_evaluation_plan(
        claimed_plan, ledger_root, require_active=True
    )
    if entry["execution_start"] is None:
        raise C10LedgerIntegrityError(
            "canonical release requires the durable one-shot execution-start marker"
        )
    claim = entry["claim"]
    binding = claim["binding"]
    payload = dict(result)
    if payload.get("schema_version") != CANONICAL_RESULT_SCHEMA:
        raise ValueError("canonical result schema mismatch")
    if payload.get("canonical_mask_allowed") is not True:
        raise ValueError("canonical result is not releasable")
    if payload.get("selection_stage") != "exact_confirmation":
        raise ValueError("canonical result is not exact confirmation")
    for field, expected in (
        ("stable", True),
        ("exact_confirmation_passed", True),
        ("every_confirmation_constraint_passed", True),
        ("official_validation_used", False),
        ("official_test_used", False),
    ):
        if payload.get(field) is not expected:
            raise ValueError(f"canonical result requires {field}={expected!r}")
    for field in (
        "constraint_registry_sha256",
        "mask_sha256",
        "selected_local_generator_ids",
    ):
        if payload.get(field) != binding[field]:
            raise ValueError(f"canonical result {field} differs from C10 claim")
    if payload.get("hard_active") != binding["preselected_k"]:
        raise ValueError("canonical result K differs from C10 claim")
    if payload.get("candidate_count") != binding["candidate_count"] or payload.get(
        "n_generators"
    ) != binding["candidate_count"]:
        raise ValueError("canonical result candidate count differs from C10 claim")
    if payload.get("donor_split_manifest") != claimed_plan.get(
        "donor_split_manifest"
    ):
        raise ValueError("canonical result embeds a different donor split manifest")
    if payload.get("preselection_manifest") != claimed_plan.get(
        "preselection_manifest"
    ):
        raise ValueError("canonical result embeds a different preselection manifest")
    source_hashes = payload.get("source_hashes")
    if not isinstance(source_hashes, Mapping):
        raise ValueError("canonical result lacks source_hashes")
    for result_field, binding_field in (
        ("donor_split_manifest_sha256", "donor_split_manifest_sha256"),
        ("preselection_manifest_sha256", "preselection_manifest_sha256"),
        ("ranking_result_sha256", "ranking_result_sha256"),
    ):
        if source_hashes.get(result_field) != binding[binding_field]:
            raise ValueError(f"canonical result {result_field} differs from C10 claim")
    confirmation_input_sha256 = _require_sha256(
        source_hashes.get("confirmation_input_sha256"),
        label="canonical result confirmation_input_sha256",
    )
    supplied_result_sha256 = payload.get("canonical_payload_sha256")
    expected_result_sha256 = canonical_json_sha256(
        {
            key: value
            for key, value in payload.items()
            if key != "canonical_payload_sha256"
        }
    )
    if supplied_result_sha256 != expected_result_sha256:
        raise ValueError("canonical result payload SHA256 mismatch")

    payload["c10_consumption_evidence"] = {
        "claim_sha256": claim["claim_sha256"],
        "binding_sha256": binding["binding_sha256"],
        "c10_partition_sha256": binding["c10_partition_sha256"],
        "plan_sha256": claim["plan_sha256"],
        "plan_evidence": claim["plan_evidence"],
        "confirmation_input_sha256": confirmation_input_sha256,
        "one_shot_claim_verified": True,
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    payload["canonical_payload_sha256"] = canonical_json_sha256(
        {
            key: value
            for key, value in payload.items()
            if key != "canonical_payload_sha256"
        }
    )
    record_c10_outcome(
        ledger_root,
        claim,
        outcome="completed",
        evidence_sha256=payload["canonical_payload_sha256"],
        confirmation_input_sha256=confirmation_input_sha256,
        canonical_output_sha256=payload["canonical_payload_sha256"],
        detail="canonical generator result passed and was bound to this one-shot C10 claim",
    )
    return payload


def verify_c10_canonical_release_evidence(
    result: Mapping[str, Any],
    claimed_plan: Mapping[str, Any],
    ledger_root: str | Path,
) -> dict[str, Any]:
    """Audit a persisted canonical result against its completed ledger chain."""

    entry = verify_claimed_exact_evaluation_plan(
        claimed_plan, ledger_root, require_active=False
    )
    if entry["state"] != "completed" or entry["execution_start"] is None:
        raise C10LedgerIntegrityError(
            "canonical release evidence requires a completed, execution-started claim"
        )
    claim = entry["claim"]
    binding = claim["binding"]
    outcome = entry["outcome"]
    payload = dict(result)
    if payload.get("schema_version") != CANONICAL_RESULT_SCHEMA:
        raise C10LedgerIntegrityError("canonical release schema mismatch")
    for field, expected in (
        ("stable", True),
        ("selection_stage", "exact_confirmation"),
        ("exact_confirmation_passed", True),
        ("every_confirmation_constraint_passed", True),
        ("canonical_mask_allowed", True),
        ("official_validation_used", False),
        ("official_test_used", False),
        ("candidate_count", binding["candidate_count"]),
        ("n_generators", binding["candidate_count"]),
        ("hard_active", binding["preselected_k"]),
        ("selected_local_generator_ids", binding["selected_local_generator_ids"]),
        ("mask_sha256", binding["mask_sha256"]),
        ("constraint_registry_sha256", binding["constraint_registry_sha256"]),
    ):
        if payload.get(field) != expected:
            raise C10LedgerIntegrityError(
                f"canonical release {field} differs from completed C10 claim"
            )
    if payload.get("donor_split_manifest") != claimed_plan.get(
        "donor_split_manifest"
    ) or payload.get("preselection_manifest") != claimed_plan.get(
        "preselection_manifest"
    ):
        raise C10LedgerIntegrityError(
            "canonical release embedded split/preselection differs from claimed plan"
        )
    source_hashes = payload.get("source_hashes")
    if not isinstance(source_hashes, Mapping):
        raise C10LedgerIntegrityError("canonical release lacks source hashes")
    for field in (
        "donor_split_manifest_sha256",
        "preselection_manifest_sha256",
        "ranking_result_sha256",
    ):
        if source_hashes.get(field) != binding[field]:
            raise C10LedgerIntegrityError(
                f"canonical release source {field} differs from C10 claim"
            )
    try:
        input_sha256 = _require_sha256(
            source_hashes.get("confirmation_input_sha256"),
            label="canonical release confirmation_input_sha256",
        )
    except ValueError as exc:
        raise C10LedgerIntegrityError(str(exc)) from exc
    expected_release_evidence = {
        "claim_sha256": claim["claim_sha256"],
        "binding_sha256": binding["binding_sha256"],
        "c10_partition_sha256": binding["c10_partition_sha256"],
        "plan_sha256": claim["plan_sha256"],
        "plan_evidence": claim["plan_evidence"],
        "confirmation_input_sha256": input_sha256,
        "one_shot_claim_verified": True,
        "filesystem_limitation": FILESYSTEM_LIMITATION,
    }
    if payload.get("c10_consumption_evidence") != expected_release_evidence:
        raise C10LedgerIntegrityError(
            "canonical release C10 evidence differs from immutable claim"
        )
    supplied_output_sha256 = payload.get("canonical_payload_sha256")
    expected_output_sha256 = canonical_json_sha256(
        {
            key: value
            for key, value in payload.items()
            if key != "canonical_payload_sha256"
        }
    )
    if supplied_output_sha256 != expected_output_sha256:
        raise C10LedgerIntegrityError("canonical release output SHA256 mismatch")
    if (
        outcome.get("confirmation_input_sha256") != input_sha256
        or outcome.get("canonical_output_sha256") != supplied_output_sha256
        or outcome.get("evidence_sha256") != supplied_output_sha256
        or outcome.get("execution_start_sha256")
        != entry["execution_start"]["execution_start_sha256"]
    ):
        raise C10LedgerIntegrityError(
            "canonical release input/output evidence differs from terminal outcome"
        )
    return entry


def record_c10_confirmation_failure(
    claimed_plan: Mapping[str, Any],
    ledger_root: str | Path,
    *,
    failure_evidence_sha256: str,
    detail: str,
    confirmation_input_sha256: str | None = None,
) -> dict[str, Any]:
    """Append terminal failure evidence for a claimed evaluation plan."""

    entry = verify_claimed_exact_evaluation_plan(
        claimed_plan, ledger_root, require_active=True
    )
    return record_c10_outcome(
        ledger_root,
        entry["claim"],
        outcome="failed",
        evidence_sha256=failure_evidence_sha256,
        confirmation_input_sha256=confirmation_input_sha256,
        detail=detail,
    )


def build_claimed_exact_evaluation_plan(
    *,
    c10_ledger_root: str | Path,
    **plan_kwargs: Any,
) -> dict[str, Any]:
    """Build (but do not execute) a confirmation plan and atomically claim C10.

    This convenience entry point derives the canonical preselection manifest
    from the validated exact-evaluation plan.  It intentionally rejects ranking
    plans because R10 exploration does not consume C10.
    """

    if plan_kwargs.get("donor_partition") != "confirmation":
        raise ValueError("claimed plan builder requires donor_partition='confirmation'")
    from kmlee_bam.training.generator_exact_confirmation import (
        build_generator_k_preselection_manifest,
    )
    from kmlee_bam.training.generator_exact_evaluation import (
        build_exact_evaluation_plan,
    )

    plan = build_exact_evaluation_plan(**plan_kwargs)
    split_path = Path(plan_kwargs["donor_split_manifest_path"])
    try:
        split = json.loads(split_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read donor split manifest: {split_path}") from exc
    if not isinstance(split, dict):
        raise ValueError("donor split manifest must contain a JSON object")
    jobs = plan.get("evaluation_jobs", [])
    candidates = [job for job in jobs if job.get("all_on_reference") is False]
    if len(candidates) != 1:
        raise ValueError("confirmation plan must have exactly one non-all-on candidate")
    candidate = candidates[0]
    preselection = build_generator_k_preselection_manifest(
        selected_generator_ids=candidate["selected_local_generator_ids"],
        candidate_count=candidate["candidate_count"],
        ranking_result_sha256=plan["ranking_result_sha256"],
        donor_split_manifest_sha256=plan["donor_split_manifest_sha256"],
    )
    plan_with_split = dict(plan)
    # Embedding the already hash-sealed manifest makes later release verification
    # independent of mutable source paths.
    plan_with_split["donor_split_manifest"] = split
    plan_with_split["plan_sha256"] = canonical_json_sha256(
        {
            key: value
            for key, value in plan_with_split.items()
            if key != "plan_sha256"
        }
    )
    return claim_exact_evaluation_plan(
        plan_with_split,
        c10_ledger_root,
        donor_split_manifest=split,
        preselection_manifest=preselection,
    )


def build_canonical_generator_result_with_c10_ledger(
    selected: Any,
    *,
    c10_ledger_root: str | Path,
    claimed_evaluation_plan: Mapping[str, Any],
    **release_kwargs: Any,
) -> dict[str, Any]:
    """Build and terminally complete a canonical result through the ledger gate."""

    from kmlee_bam.training.generator_exact_confirmation import (
        build_canonical_generator_result,
    )

    verify_claimed_exact_evaluation_plan(
        claimed_evaluation_plan, c10_ledger_root, require_active=True
    )
    result = build_canonical_generator_result(selected, **release_kwargs)
    return attach_c10_canonical_release_evidence(
        result,
        claimed_evaluation_plan,
        c10_ledger_root,
    )
