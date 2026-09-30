"""Fail-closed orchestration contract for exact generator-mask evaluation.

This module intentionally does not perform model inference.  It builds frozen,
train-only GPU job manifests from a Hard-Concrete *ranking* and validates the
per-donor results emitted by an inference backend.  Separating orchestration
from inference prevents an evaluator from silently opening validation/test
donors, changing donor order between masks, or substituting an aggregate metric
that cannot be paired by donor.

Nonlinear CCC/correlation/slope metrics are intentionally *not* emitted as
pre-aggregated donor scalars.  R10 jobs preserve raw sufficient statistics;
C10 jobs use only R10-frozen additive surrogate losses from the exact canonical
registry.  The actual callback/GPU executor remains external and the generated
plan therefore cannot itself release a canonical mask.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from kmlee_bam.training.generator_exact_confirmation import (
    CANONICAL_V4_METRIC_REGISTRY,
    CONSENSUS_SCHEMA,
    RESULT_SCHEMA,
    SPLIT_SCHEMA,
    canonical_json_sha256,
    generator_mask_sha256,
    nested_topk_generator_ids,
    validate_canonical_v4_donor_split_manifest,
)


PLAN_SCHEMA = "kmlee_bam.generator_exact_evaluation_plan.v1"
JOB_SCHEMA = "kmlee_bam.generator_exact_evaluation_job.v1"
OUTPUT_SCHEMA = "kmlee_bam.generator_exact_evaluation_output.v1"
RESAMPLING_SCHEMA = "kmlee_bam.generator_exact_resampling_sufficient_stats.v1"
R10_SELECTION_SCHEMA = "kmlee_bam.generator_r10_selection.v1"
W44_SOURCE_SCHEMA = "kmlee_bam.generator_weight_only_provenance.v4"

NON_NEURONAL_CELLTYPES = (
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
)
PATHOLOGY_AXES = ("Braak", "Thal", "CERAD", "LATE", "Lewy")


def _raw_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_sha256(value: Any, *, label: str) -> str:
    digest = str(value).lower()
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise ValueError(f"{label} must be a lowercase SHA256 digest")
    return digest


def _read_json_object(path: str | Path, *, label: str) -> tuple[Path, dict[str, Any]]:
    source = Path(path).resolve()
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {label} JSON: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return source, value


def metric_contract() -> tuple[dict[str, str], ...]:
    """Return exact-core's immutable 24-metric additive registry."""

    return tuple(
        {
            "name": name,
            "category": category,
            "direction": direction,
            "estimator": "frozen_r10_defined_additive_donor_loss",
        }
        for name, (category, direction) in sorted(
            CANONICAL_V4_METRIC_REGISTRY.items()
        )
    )


def _validate_split_manifest(split: Mapping[str, Any]) -> dict[str, Any]:
    result = validate_canonical_v4_donor_split_manifest(split)
    if result.get("schema_version") != SPLIT_SCHEMA:
        raise ValueError("donor split manifest schema mismatch")
    if result.get("source_split") != "train_only":
        raise ValueError("generator evaluation requires source_split='train_only'")
    if result.get("disjoint") is not True or result.get("complete") is not True:
        raise ValueError("donor split manifest must be disjoint and complete")
    if result.get("official_validation_used") is not False:
        raise ValueError("official validation donors are forbidden")
    if result.get("official_test_used") is not False:
        raise ValueError("official test donors are forbidden")
    supplied = _require_sha256(result.get("manifest_sha256"), label="manifest_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "manifest_sha256"}
    )
    if supplied != expected:
        raise ValueError("donor split manifest SHA256 mismatch")
    expected_counts = {
        "weight_donor_names": 44,
        "ranking_donor_names": 10,
        "confirmation_donor_names": 10,
    }
    for field, expected_count in expected_counts.items():
        names = result.get(field)
        if not isinstance(names, list) or len(names) != expected_count:
            raise ValueError(
                f"{field} must contain exactly {expected_count} donor names"
            )
        parsed = tuple(str(name) for name in names)
        if any(not name for name in parsed) or len(set(parsed)) != len(parsed):
            raise ValueError(f"{field} must contain unique non-empty names")
    partitions = [set(result[field]) for field in expected_counts]
    if any(
        partitions[left] & partitions[right]
        for left in range(3)
        for right in range(left + 1, 3)
    ):
        raise ValueError("weight, ranking, and confirmation donors overlap")
    return result


def _validate_ranking(ranking: Mapping[str, Any]) -> tuple[np.ndarray, tuple[int, ...], tuple[str, ...], tuple[int, ...]]:
    if ranking.get("schema_version") != RESULT_SCHEMA:
        raise ValueError("ranking result schema mismatch")
    if ranking.get("selection_stage") != "ranking_only":
        raise ValueError("exact evaluation requires a ranking_only result")
    if ranking.get("consensus_schema_version") != CONSENSUS_SCHEMA:
        raise ValueError("exact evaluation requires a two-seed consensus ranking")
    if ranking.get("ranking_method") != (
        "median_normalized_rank_two_independent_seeds"
    ):
        raise ValueError("exact evaluation ranking method is not canonical consensus")
    if ranking.get("stability_gate_passed") is not True:
        raise ValueError("exact evaluation requires the ranking stability gate")
    if ranking.get("canonical_mask_allowed") is not False:
        raise ValueError("ranking result must not claim a canonical mask")
    if ranking.get("exact_confirmation_required") is not True:
        raise ValueError("ranking result must require exact confirmation")
    if ranking.get("official_validation_used") is not False:
        raise ValueError("ranking used official validation donors")
    if ranking.get("official_test_used") is not False:
        raise ValueError("ranking used official test donors")
    count = ranking.get("candidate_count")
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError("ranking candidate_count must be a positive integer")
    scores = np.asarray(ranking.get("ranking_scores"), dtype=np.float64)
    if scores.shape != (count,) or not bool(np.isfinite(scores).all()):
        raise ValueError("ranking_scores must be finite with candidate_count entries")
    if ranking.get("ranking_score_direction") != "higher_is_more_important":
        raise ValueError("ranking score direction must be higher_is_more_important")
    indices_raw = ranking.get("candidate_registry_indices")
    names_raw = ranking.get("candidate_registry_module_names")
    if not isinstance(indices_raw, list) or not isinstance(names_raw, list):
        raise ValueError("ranking must preserve the complete candidate registry mapping")
    if len(indices_raw) != count or len(names_raw) != count:
        raise ValueError("candidate registry mapping length mismatch")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in indices_raw):
        raise ValueError("candidate_registry_indices must contain integers")
    indices = tuple(int(value) for value in indices_raw)
    names = tuple(str(value) for value in names_raw)
    if len(set(indices)) != count or any(not value for value in names):
        raise ValueError("candidate registry mapping is not unique/non-empty")
    protected_raw = ranking.get("protected_local_generator_ids", [])
    if not isinstance(protected_raw, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in protected_raw
    ):
        raise ValueError("protected_local_generator_ids must contain integers")
    protected = tuple(int(value) for value in protected_raw)
    if len(set(protected)) != len(protected) or any(
        value < 0 or value >= count for value in protected
    ):
        raise ValueError("protected generator ids are duplicate/out of range")
    return scores, indices, names, protected


def build_exact_evaluation_plan(
    *,
    config_path: str | Path,
    checkpoint_path: str | Path,
    registry_path: str | Path,
    scaler_path: str | Path,
    module_stats_path: str | Path,
    donor_context_path: str | Path,
    age_normalization_path: str | Path,
    weight_only_source_manifest_path: str | Path,
    frozen_signature_path: str | Path | None = None,
    metric_denominator_path: str | Path | None = None,
    ranking_result_path: str | Path,
    donor_split_manifest_path: str | Path,
    requested_counts: Iterable[int],
    donor_partition: str,
    r10_selection_result_path: str | Path | None = None,
    adaptation_epochs: int = 4,
    confirmation_retrain_epochs: int = 6,
    training_seed: int = 420807,
) -> dict[str, Any]:
    """Build fixed-mask adaptation/retrain plus evaluation jobs.

    ``ranking`` builds one W44 adaptation and one R10 evaluation for every K.
    ``confirmation`` accepts exactly one frozen non-all-on K and builds only
    that K and all-on as paired W+R54 retrains followed by one-shot C10
    evaluation.  Acute decoder-mask ablation is never an allowed evaluation.
    """

    if donor_partition not in {"ranking", "confirmation"}:
        raise ValueError("donor_partition must be 'ranking' or 'confirmation'")
    for label, value in (
        ("adaptation_epochs", adaptation_epochs),
        ("confirmation_retrain_epochs", confirmation_retrain_epochs),
    ):
        if isinstance(value, bool) or int(value) <= 0:
            raise ValueError(f"{label} must be a positive integer")
    if isinstance(training_seed, bool) or int(training_seed) < 0:
        raise ValueError("training_seed must be a non-negative integer")
    paths: dict[str, Path] = {}
    documents: dict[str, dict[str, Any]] = {}
    for label, path in (
        ("config", config_path),
        ("ranking_result", ranking_result_path),
        ("donor_split_manifest", donor_split_manifest_path),
    ):
        paths[label], documents[label] = _read_json_object(path, label=label)
    paths["weight_only_source_manifest"], documents[
        "weight_only_source_manifest"
    ] = _read_json_object(
        weight_only_source_manifest_path,
        label="weight_only_source_manifest",
    )
    if donor_partition == "confirmation":
        if r10_selection_result_path is None:
            raise ValueError(
                "confirmation requires a frozen R10 selection result"
            )
        paths["r10_selection_result"], documents["r10_selection_result"] = (
            _read_json_object(
                r10_selection_result_path,
                label="r10_selection_result",
            )
        )
        if frozen_signature_path is None or metric_denominator_path is None:
            raise ValueError(
                "confirmation requires frozen R10 signature and metric denominator artifacts"
            )
        paths["frozen_signature"] = Path(frozen_signature_path).resolve()
        paths["metric_denominator"] = Path(metric_denominator_path).resolve()
    paths["checkpoint"] = Path(checkpoint_path).resolve()
    paths["registry"] = Path(registry_path).resolve()
    paths["scaler"] = Path(scaler_path).resolve()
    paths["module_stats"] = Path(module_stats_path).resolve()
    paths["donor_context"] = Path(donor_context_path).resolve()
    paths["age_normalization"] = Path(age_normalization_path).resolve()
    for label in (
        "checkpoint",
        "registry",
        "scaler",
        "module_stats",
        "donor_context",
        "age_normalization",
    ):
        if not paths[label].is_file():
            raise ValueError(f"{label} file does not exist: {paths[label]}")
    for label in ("frozen_signature", "metric_denominator"):
        if label in paths and not paths[label].is_file():
            raise ValueError(f"{label} file does not exist: {paths[label]}")

    split = _validate_split_manifest(documents["donor_split_manifest"])
    scores, registry_indices, registry_names, protected = _validate_ranking(
        documents["ranking_result"]
    )
    requested = tuple(sorted({int(value) for value in requested_counts}))
    if not requested:
        raise ValueError("requested_counts must be non-empty")
    if donor_partition == "confirmation":
        non_all_on = tuple(value for value in requested if value != int(scores.size))
        if len(non_all_on) != 1 or len(requested) != len(non_all_on):
            raise ValueError(
                "confirmation requires exactly one frozen non-all-on K; "
                "all-on is added automatically"
            )
        counts = (non_all_on[0], int(scores.size))
    else:
        counts = tuple(sorted(set(requested) | {int(scores.size)}))
    masks = nested_topk_generator_ids(
        scores,
        counts,
        protected_generator_ids=protected,
    )
    evaluation_donors = tuple(split[f"{donor_partition}_donor_names"])
    weight_donors = tuple(split["weight_donor_names"])
    ranking_donors = tuple(split["ranking_donor_names"])
    if donor_partition == "ranking":
        training_donors = weight_donors
        training_kind = "fixed_mask_adaptation_w44"
        evaluation_kind = "ranking_evaluation_r10"
    else:
        training_donors = weight_donors + ranking_donors
        training_kind = "paired_fixed_mask_retrain_w_plus_r54"
        evaluation_kind = "one_shot_confirmation_evaluation_c10"
    source_hashes = {
        label: _raw_sha256(path)
        for label, path in paths.items()
    }
    source_manifest = documents["weight_only_source_manifest"]
    if source_manifest.get("schema_version") != W44_SOURCE_SCHEMA:
        raise ValueError("W44 source artifact manifest schema mismatch")
    provenance_expected = {
        "checkpoint_sha256": source_hashes["checkpoint"],
        "module_stats_sha256": source_hashes["module_stats"],
        "source_config_sha256": source_hashes["config"],
        "registry_sha256": source_hashes["registry"],
        "scaler_sha256": source_hashes["scaler"],
        "donor_context_sha256": source_hashes["donor_context"],
        "age_normalization_sha256": source_hashes["age_normalization"],
        "donor_split_manifest_sha256": split["manifest_sha256"],
        "weight_donor_sha256": split["weight_donor_sha256"],
        "weight_donor_names_sha256": split["weight_donor_sha256"],
        "weight_donors_only": True,
        "model_weight_fit_partition": "W44",
        "context_fit_partition": "W44",
        "module_stats_fit_partition": "W44",
        "official_validation_used": False,
        "official_test_used": False,
    }
    for field, expected in provenance_expected.items():
        if source_manifest.get(field) != expected:
            raise ValueError(f"W44 source provenance {field} mismatch")
    source_manifest_sha = canonical_json_sha256(
        {
            key: value
            for key, value in source_manifest.items()
            if key != "manifest_sha256"
        }
    )
    if source_manifest.get("manifest_sha256") != source_manifest_sha:
        raise ValueError("W44 source artifact manifest SHA256 mismatch")
    ranking_declared_checkpoint = _require_sha256(
        documents["ranking_result"].get("source_checkpoint_sha256"),
        label="source_checkpoint_sha256",
    )
    if ranking_declared_checkpoint != source_hashes["checkpoint"]:
        raise ValueError("ranking source checkpoint SHA256 mismatch")
    ranking_declared_registry = _require_sha256(
        documents["ranking_result"].get("registry_sha256"),
        label="registry_sha256",
    )
    if ranking_declared_registry != source_hashes["registry"]:
        raise ValueError("ranking registry SHA256 mismatch")
    artifact_provenance = documents["ranking_result"].get("artifact_provenance")
    if not isinstance(artifact_provenance, Mapping):
        raise ValueError("ranking result lacks W44 artifact_provenance")
    if dict(artifact_provenance) != source_manifest:
        raise ValueError(
            "ranking artifact_provenance differs from the audited W44 source manifest"
        )
    if donor_partition == "confirmation":
        selection = documents["r10_selection_result"]
        if selection.get("schema_version") != R10_SELECTION_SCHEMA:
            raise ValueError("R10 selection result schema mismatch")
        if selection.get("selection_stage") != "R10_grid_selection":
            raise ValueError("confirmation K was not selected on R10")
        if selection.get("selection_locked") is not True:
            raise ValueError("R10 final K is not frozen")
        if selection.get("official_confirmation_used") is not False:
            raise ValueError("R10 selection opened confirmation donors")
        if selection.get("official_validation_used") is not False:
            raise ValueError("R10 selection opened official validation")
        if selection.get("official_test_used") is not False:
            raise ValueError("R10 selection opened official test")
        if selection.get("next_k_after_confirmation_failure_forbidden") is not True:
            raise ValueError("R10 result does not forbid post-confirmation K fallback")
        frozen_k = selection.get("selected_generator_count")
        if frozen_k != counts[0]:
            raise ValueError("requested confirmation K differs from frozen R10 K")
        if selection.get("candidate_count") != int(scores.size):
            raise ValueError("R10 selection candidate bank size mismatch")
        if selection.get("selected_local_generator_ids") != list(masks[counts[0]]):
            raise ValueError("R10 selected generator mask differs from frozen ranking")
        if selection.get("mask_sha256") != generator_mask_sha256(
            masks[counts[0]], candidate_count=int(scores.size)
        ):
            raise ValueError("R10 selected mask SHA256 mismatch")
        if selection.get("hard_concrete_ranking_result_sha256") != source_hashes[
            "ranking_result"
        ]:
            raise ValueError("R10 selection ranking-result SHA256 mismatch")
        if selection.get("donor_split_manifest_sha256") != split["manifest_sha256"]:
            raise ValueError("R10 selection donor-split SHA256 mismatch")
        for field, expected in (
            ("source_checkpoint_sha256", source_hashes["checkpoint"]),
            ("registry_sha256", source_hashes["registry"]),
        ):
            if selection.get(field) != expected:
                raise ValueError(f"R10 selection {field} mismatch")
        for field, source_name in (
            ("frozen_signature_sha256", "frozen_signature"),
            ("metric_denominator_sha256", "metric_denominator"),
        ):
            if selection.get(field) != source_hashes[source_name]:
                raise ValueError(f"R10 selection {field} mismatch")
        _require_sha256(
            selection.get("constraint_registry_sha256"),
            label="constraint_registry_sha256",
        )
        output_hashes = selection.get("r10_evaluation_output_sha256_by_k")
        if not isinstance(output_hashes, Mapping):
            raise ValueError("R10 selection lacks evaluation output hashes")
        for count in (counts[0], int(scores.size)):
            _require_sha256(
                output_hashes.get(str(count)),
                label=f"R10 evaluation output SHA256 for K={count}",
            )

    registry_document = documents.get("registry")
    if registry_document is None:
        _, registry_document = _read_json_object(paths["registry"], label="registry")
    module_names = registry_document.get("module_names")
    if not isinstance(module_names, list) or not module_names:
        raise ValueError("registry must contain a non-empty module_names list")
    for local_id, (registry_index, expected_name) in enumerate(
        zip(registry_indices, registry_names)
    ):
        if registry_index < 0 or registry_index >= len(module_names):
            raise ValueError(
                f"candidate {local_id} registry index is out of range"
            )
        if str(module_names[registry_index]) != expected_name:
            raise ValueError(
                f"candidate {local_id} registry name/index mapping mismatch"
            )

    contract = list(metric_contract())
    metric_registry_sha256 = canonical_json_sha256(contract)
    exact_constraint_registry_sha256 = (
        documents["r10_selection_result"]["constraint_registry_sha256"]
        if donor_partition == "confirmation"
        else None
    )
    training_jobs: list[dict[str, Any]] = []
    evaluation_jobs: list[dict[str, Any]] = []
    for count in counts:
        selected = masks[count]
        common = {
            "source_split": "train_only",
            "official_validation_used": False,
            "official_test_used": False,
            "candidate_count": int(scores.size),
            "generator_count": count,
            "all_on_reference": bool(count == int(scores.size)),
            "selected_local_generator_ids": list(selected),
            "selected_registry_indices": [registry_indices[index] for index in selected],
            "selected_registry_module_names": [registry_names[index] for index in selected],
            "mask_sha256": generator_mask_sha256(
                selected, candidate_count=int(scores.size)
            ),
            "source_paths": {key: str(value) for key, value in sorted(paths.items())},
            "source_sha256": dict(sorted(source_hashes.items())),
        }
        training_job: dict[str, Any] = {
            "schema_version": JOB_SCHEMA,
            "job_kind": training_kind,
            "donor_partition": donor_partition,
            **common,
            "training_donor_names": list(training_donors),
            "training_donor_names_sha256": canonical_json_sha256(
                list(training_donors)
            ),
            "training_contract": {
                "epochs": int(
                    adaptation_epochs
                    if donor_partition == "ranking"
                    else confirmation_retrain_epochs
                ),
                "seed": int(training_seed),
                "fresh_optimizer": True,
                "fixed_mask_for_entire_run": True,
                "mask_trainable": False,
                "initial_weights_identical_source_checkpoint": True,
                "same_seed_schedule_across_k": True,
                "acute_mask_ablation": False,
                "learned_generator_count_enabled": False,
                "legacy_generator_budget_search_enabled": False,
                "official_validation_evaluation": False,
                "official_test_evaluation": False,
                "output": "fixed_mask_adapted_checkpoint_and_manifest",
            },
        }
        if donor_partition == "confirmation":
            training_job["paired_retrain_group"] = (
                "frozen_final_k_vs_all_on_w_plus_r54"
            )
        training_job["job_sha256"] = canonical_json_sha256(training_job)
        training_jobs.append(training_job)

        evaluation_job: dict[str, Any] = {
            "schema_version": JOB_SCHEMA,
            "job_kind": evaluation_kind,
            "donor_partition": donor_partition,
            **common,
            "upstream_training_job_sha256": training_job["job_sha256"],
            "donor_names": list(evaluation_donors),
            "donor_names_sha256": canonical_json_sha256(
                list(evaluation_donors)
            ),
            "registry_module_count": len(module_names),
            "required_nonneuronal_celltypes": list(NON_NEURONAL_CELLTYPES),
            "inference_contract": {
                "model_weights_frozen": True,
                "sample_latent": False,
                "shuffle": False,
                "donor_order": "exactly_as_listed",
                "mask_application": (
                    "none_during_evaluation;mask_already_fixed_for_upstream_training"
                ),
                "acute_mask_ablation": False,
                "full_view": "upstream_fixed_mask_trained_decoder",
                "isolated_view": (
                    "base_plus_tech_plus_sex_plus_masked_generator_only;"
                    "mixer_direct_pathology_interaction_PRISM_precision_excluded"
                ),
                "bootstrap_unit": "donor",
                "nonlinear_metric_policy": (
                    "recompute_from_raw_sufficient_statistics_each_resample"
                ),
                "precomputed_donor_scalar_for_ccc_correlation_slope": False,
            },
            "metric_contract": contract,
            "metric_registry_sha256": metric_registry_sha256,
        }
        if donor_partition == "ranking":
            evaluation_job["resampling_contract"] = {
                "schema_version": RESAMPLING_SCHEMA,
                "callback": "kmlee_bam.generator_exact_metric_callback.v1",
                "shared_resample_indices_for_candidate_and_all_on": True,
                "raw_module_sums_required": True,
                "raw_pathology_and_covariates_required": True,
                "additive_nll_sums_required": True,
                "additive_leakage_sums_required": True,
            }
        else:
            evaluation_job["additive_confirmation_contract"] = {
                "metric_registry": "CANONICAL_V4_METRIC_REGISTRY",
                "constraint_registry_sha256": exact_constraint_registry_sha256,
                "frozen_signature_sha256": source_hashes["frozen_signature"],
                "metric_denominator_sha256": source_hashes[
                    "metric_denominator"
                ],
                "donor_scalar_mean_for_nonlinear_metric": False,
                "all_metrics_are_frozen_additive_losses": True,
                "full_and_isolated_nll_in_exact_registry": True,
            }
        if donor_partition == "confirmation":
            evaluation_job["inference_contract"]["nonlinear_metric_policy"] = (
                "frozen_R10_defined_additive_losses_only"
            )
            evaluation_job["paired_retrain_group"] = (
                "frozen_final_k_vs_all_on_w_plus_r54"
            )
        evaluation_job["job_sha256"] = canonical_json_sha256(evaluation_job)
        evaluation_jobs.append(evaluation_job)

    plan: dict[str, Any] = {
        "schema_version": PLAN_SCHEMA,
        "status": "training_and_evaluation_jobs_built_not_run",
        "executor": "required_external_fixed_mask_training_and_gpu_backend",
        "executor_ready": False,
        "canonical_mask_allowed": False,
        "blocked_until": (
            "audited nonlinear donor-resampling callback and GPU executor are implemented"
        ),
        "donor_partition": donor_partition,
        "source_split": "train_only",
        "official_validation_used": False,
        "official_test_used": False,
        "ranking_result_sha256": source_hashes["ranking_result"],
        "metric_registry_sha256": metric_registry_sha256,
        "exact_constraint_registry_sha256": exact_constraint_registry_sha256,
        "donor_split_manifest_sha256": split["manifest_sha256"],
        "requested_generator_counts": list(counts),
        "training_schedule": {
            "adaptation_epochs": int(adaptation_epochs),
            "confirmation_retrain_epochs": int(confirmation_retrain_epochs),
            "seed": int(training_seed),
            "same_across_compared_k": True,
        },
        "training_jobs": training_jobs,
        "evaluation_jobs": evaluation_jobs,
        "selection_policy": (
            {
                "stage": "R10_grid_selection",
                "smallest_feasible_k_selected_once": True,
                "confirmation_donors_opened": False,
            }
            if donor_partition == "ranking"
            else {
                "stage": "C10_one_shot_confirmation",
                "compared_models": "frozen_final_k_vs_all_on_only",
                "next_k_after_failure_forbidden": True,
                "failure_is_terminal_no_canonical_mask": True,
            }
        ),
    }
    plan["plan_sha256"] = canonical_json_sha256(plan)
    return plan


def validate_exact_evaluation_output(
    output: Mapping[str, Any],
    job: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate one GPU result and return paired-bootstrap-ready arrays."""

    result = dict(output)
    expected_job = dict(job)
    if result.get("schema_version") != OUTPUT_SCHEMA:
        raise ValueError("evaluation output schema mismatch")
    if expected_job.get("schema_version") != JOB_SCHEMA:
        raise ValueError("evaluation job schema mismatch")
    if expected_job.get("job_kind") not in {
        "ranking_evaluation_r10",
        "one_shot_confirmation_evaluation_c10",
    }:
        raise ValueError("supplied job is not an evaluation job")
    supplied_job_sha = expected_job.get("job_sha256")
    if supplied_job_sha != canonical_json_sha256(
        {key: value for key, value in expected_job.items() if key != "job_sha256"}
    ):
        raise ValueError("evaluation job SHA256 mismatch")
    if result.get("job_sha256") != supplied_job_sha:
        raise ValueError("output was not produced for the supplied job")
    for flag in ("official_validation_used", "official_test_used"):
        if result.get(flag) is not False:
            raise ValueError(f"evaluation output has unsafe {flag}")
    if result.get("source_split") != "train_only":
        raise ValueError("evaluation output must declare train_only")
    if result.get("deterministic_inference") is not True:
        raise ValueError("evaluation output must use deterministic inference")
    if result.get("model_weights_frozen") is not True:
        raise ValueError("evaluation output must freeze model weights")
    if result.get("acute_mask_ablation") is not False:
        raise ValueError("acute-mask evaluation is forbidden")
    if result.get("source_sha256") != expected_job.get("source_sha256"):
        raise ValueError("evaluation source hashes differ from the job")
    for field in (
        "donor_partition",
        "candidate_count",
        "generator_count",
        "selected_local_generator_ids",
        "mask_sha256",
    ):
        if result.get(field) != expected_job.get(field):
            raise ValueError(f"evaluation output {field} differs from the job")
    donor_ids = tuple(str(value) for value in result.get("donor_names", ()))
    if list(donor_ids) != expected_job.get("donor_names"):
        raise ValueError("evaluation donor names/order differ from the job")
    if len(donor_ids) < 2 or len(set(donor_ids)) != len(donor_ids):
        raise ValueError("evaluation donor names must be unique")
    if result.get("donor_names_sha256") != expected_job.get("donor_names_sha256"):
        raise ValueError("evaluation donor-name SHA256 mismatch")
    if result.get("inference_contract") != expected_job.get("inference_contract"):
        raise ValueError("evaluation inference contract differs from the job")

    if result.get("upstream_training_job_sha256") != expected_job.get(
        "upstream_training_job_sha256"
    ):
        raise ValueError("output is not from the required fixed-mask training job")
    _require_sha256(
        result.get("upstream_training_manifest_sha256"),
        label="upstream_training_manifest_sha256",
    )
    _require_sha256(
        result.get("evaluated_checkpoint_sha256"),
        label="evaluated_checkpoint_sha256",
    )

    expected_contract = expected_job.get("metric_contract")
    if expected_contract != list(metric_contract()):
        raise ValueError("evaluation job metric contract is not canonical")
    expected_metric_hash = canonical_json_sha256(expected_contract)
    if expected_job.get("metric_registry_sha256") != expected_metric_hash:
        raise ValueError("evaluation job metric registry SHA256 mismatch")
    if result.get("metric_registry_sha256") != expected_metric_hash:
        raise ValueError("evaluation output metric registry SHA256 mismatch")
    if expected_job["job_kind"] == "one_shot_confirmation_evaluation_c10":
        additive_contract = expected_job.get("additive_confirmation_contract")
        if result.get("additive_confirmation_contract") != additive_contract:
            raise ValueError("frozen additive confirmation contract differs")
        metrics = result.get("additive_metrics")
        expected_names = {entry["name"] for entry in expected_contract}
        if not isinstance(metrics, Mapping) or set(metrics) != expected_names:
            raise ValueError("additive confirmation metric names differ from registry")
        validated_metrics: dict[str, list[float]] = {}
        for name in sorted(expected_names):
            values = np.asarray(metrics[name], dtype=np.float64)
            if values.shape != (len(donor_ids),) or not bool(
                np.isfinite(values).all()
            ):
                raise ValueError(
                    f"additive metric {name!r} must be finite with one value per donor"
                )
            validated_metrics[name] = [float(value) for value in values.tolist()]
        return {
            "donor_ids": donor_ids,
            "metrics": validated_metrics,
            "generator_count": int(result["generator_count"]),
            "selected_generator_ids": tuple(
                int(value) for value in result["selected_local_generator_ids"]
            ),
            "output_sha256": canonical_json_sha256(result),
        }
    if "metrics" in result:
        raise ValueError(
            "pre-aggregated donor metrics are forbidden for nonlinear confirmation"
        )
    if result.get("resampling_contract") != expected_job.get("resampling_contract"):
        raise ValueError("resampling callback contract differs from the job")
    statistics = result.get("donor_sufficient_statistics")
    if not isinstance(statistics, list) or len(statistics) != len(donor_ids):
        raise ValueError("one sufficient-statistics row is required per donor")
    module_count = int(expected_job["registry_module_count"])
    required_ct = set(expected_job["required_nonneuronal_celltypes"])
    validated_statistics: list[dict[str, Any]] = []
    for expected_donor, row in zip(donor_ids, statistics):
        if not isinstance(row, Mapping) or row.get("donor_name") != expected_donor:
            raise ValueError("sufficient-statistics donor names/order differ")
        for name in ("full_nll_sum", "isolated_nll_sum", "leakage_advantage_sum"):
            value = row.get(name)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        n_cells = row.get("n_cells")
        if isinstance(n_cells, bool) or not isinstance(n_cells, int) or n_cells <= 0:
            raise ValueError("n_cells must be a positive integer")
        celltypes = row.get("celltypes")
        if not isinstance(celltypes, Mapping) or not required_ct.issubset(celltypes):
            raise ValueError("raw module sums omit a fixed non-neuronal cell type")
        for celltype, values in celltypes.items():
            if not isinstance(values, Mapping):
                raise ValueError(f"celltype {celltype!r} statistics must be an object")
            group_count = values.get("n_cells")
            if isinstance(group_count, bool) or not isinstance(group_count, int) or group_count < 0:
                raise ValueError("celltype n_cells must be a non-negative integer")
            for field in (
                "observed_module_sum",
                "full_predicted_module_sum",
                "isolated_predicted_module_sum",
            ):
                array = np.asarray(values.get(field), dtype=np.float64)
                if array.shape != (module_count,) or not bool(np.isfinite(array).all()):
                    raise ValueError(
                        f"{celltype}.{field} must be finite with registry_module_count entries"
                    )
        pathology = row.get("pathology")
        if not isinstance(pathology, Mapping) or set(pathology) != set(PATHOLOGY_AXES):
            raise ValueError("pathology sufficient statistics are incomplete")
        for axis, observation in pathology.items():
            if not isinstance(observation, Mapping) or not isinstance(
                observation.get("observed"), bool
            ):
                raise ValueError(f"pathology {axis} observation is malformed")
            value = observation.get("value")
            if observation["observed"]:
                if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                    raise ValueError(f"observed pathology {axis} must be finite")
            elif value is not None:
                raise ValueError(f"missing pathology {axis} must have value=null")
        validated_statistics.append(dict(row))
    return {
        "donor_ids": donor_ids,
        "donor_sufficient_statistics": validated_statistics,
        "resampling_contract": dict(result["resampling_contract"]),
        "generator_count": int(result["generator_count"]),
        "selected_generator_ids": tuple(
            int(value) for value in result["selected_local_generator_ids"]
        ),
        "output_sha256": canonical_json_sha256(result),
    }


def prepare_paired_confirmation_inputs(
    *,
    candidate_output: Mapping[str, Any],
    candidate_job: Mapping[str, Any],
    all_on_output: Mapping[str, Any],
    all_on_job: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind two validated outputs into exact-core keyword arguments."""

    candidate = validate_exact_evaluation_output(candidate_output, candidate_job)
    reference = validate_exact_evaluation_output(all_on_output, all_on_job)
    if candidate_job.get("job_kind") != "one_shot_confirmation_evaluation_c10":
        raise ValueError("canonical confirmation cannot use ranking donors")
    if all_on_job.get("job_kind") != "one_shot_confirmation_evaluation_c10":
        raise ValueError("all-on reference cannot use ranking donors")
    if all_on_job.get("all_on_reference") is not True:
        raise ValueError("paired reference job is not the all-on mask")
    if candidate_job.get("all_on_reference") is True:
        raise ValueError("candidate job must not be the all-on reference")
    if candidate_job.get("candidate_count") != all_on_job.get("candidate_count"):
        raise ValueError("paired jobs use different candidate banks")
    if candidate_job.get("source_sha256") != all_on_job.get("source_sha256"):
        raise ValueError("paired jobs use different frozen source artifacts")
    if candidate_job.get("additive_confirmation_contract") != all_on_job.get(
        "additive_confirmation_contract"
    ):
        raise ValueError("paired jobs use different frozen additive signatures")
    if candidate["donor_ids"] != reference["donor_ids"]:
        raise ValueError("paired outputs use different donor order")
    if set(candidate["metrics"]) != set(reference["metrics"]):
        raise ValueError("paired outputs use different additive metric names")
    return {
        "generator_count": candidate["generator_count"],
        "selected_generator_ids": candidate["selected_generator_ids"],
        "candidate_donor_ids": candidate["donor_ids"],
        "all_on_donor_ids": reference["donor_ids"],
        "candidate_metrics": candidate["metrics"],
        "all_on_metrics": reference["metrics"],
        "constraint_registry_sha256": candidate_job[
            "additive_confirmation_contract"
        ]["constraint_registry_sha256"],
        "candidate_output_sha256": candidate["output_sha256"],
        "all_on_output_sha256": reference["output_sha256"],
    }


def write_evaluation_plan(path: str | Path, plan: Mapping[str, Any]) -> Path:
    payload = dict(plan)
    if payload.get("schema_version") != PLAN_SCHEMA:
        raise ValueError("evaluation plan schema mismatch")
    if payload.get("executor_ready") is not False:
        raise ValueError("evaluation plan must remain executor_ready=false")
    if payload.get("canonical_mask_allowed") is not False:
        raise ValueError("incomplete evaluation plan cannot allow a canonical mask")
    supplied = payload.get("plan_sha256")
    expected = canonical_json_sha256(
        {key: value for key, value in payload.items() if key != "plan_sha256"}
    )
    if supplied != expected:
        raise ValueError("evaluation plan SHA256 mismatch")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(target)
    return target
