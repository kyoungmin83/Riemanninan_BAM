"""Exact, train-only confirmation for PRISM generator-count selection.

The differentiable Hard-Concrete stage is deliberately treated as a ranking
method only.  This module turns that ranking into nested fixed masks, compares
each mask with the all-on supernet on held-out *training* donors, and releases a
canonical mask only when every registered non-inferiority constraint passes.

No model inference is performed here.  Callers are responsible for producing
one paired metric value per confirmation donor for every candidate and the
all-on reference.  Keeping inference outside this module makes the statistical
and provenance contract independently testable.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

from kmlee_bam.training.architecture_donor_split import (
    _build_donor_feature_matrix,
)


RESULT_SCHEMA = "kmlee_bam.learned_generator_result.v1"
SPLIT_SCHEMA = "kmlee_bam.generator_selection_donor_split.v1"
CONFIRMATION_SCHEMA = "kmlee_bam.generator_exact_confirmation.v1"

REQUIRED_CONSTRAINT_CATEGORIES = frozenset(
    {"overall", "nonneuronal", "celltype", "pathology", "leakage"}
)
REQUIRED_SOURCE_HASHES = frozenset(
    {
        "ranking_result_sha256",
        "confirmation_input_sha256",
        "source_checkpoint_sha256",
        "donor_split_manifest_sha256",
        "registry_sha256",
        "preselection_manifest_sha256",
    }
)
PRESELECTION_SCHEMA = "kmlee_bam.generator_K_preselection.v4"
CONSENSUS_SCHEMA = "kmlee_bam.generator_ranking_consensus.v4"
NO_COMPRESSION_SCHEMA = "kmlee_bam.generator_no_compression_certificate.v4"

# Canonical-v4 uses donor-additive losses only.  In particular, a mean of
# donor-level CCC/Pearson values is not a bootstrap replicate of the pooled
# CCC/Pearson statistic.  Non-linear summaries are therefore represented by
# frozen, ranking-donor-defined additive substitutes on confirmation donors.
CANONICAL_V4_NON_NEURONAL = (
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
)
CANONICAL_V4_PATHOLOGY_AXES = ("THAL", "BRAAK", "CERAD", "LATE", "LEWY")
CANONICAL_V4_METRIC_REGISTRY: dict[str, tuple[str, str]] = {
    "overall.centered_huber": ("overall", "lower"),
    "decoder.full_nll": ("overall", "lower"),
    "decoder.isolated_nll": ("overall", "lower"),
    "neuronal.frozen_mean_centered_huber": ("celltype", "lower"),
    "worst_quartile.frozen_mean_centered_huber": ("celltype", "lower"),
    **{
        f"nonneuronal.{celltype}.centered_huber": ("nonneuronal", "lower")
        for celltype in CANONICAL_V4_NON_NEURONAL
    },
    **{
        f"nonneuronal.{celltype}.ad_huber": ("nonneuronal", "lower")
        for celltype in CANONICAL_V4_NON_NEURONAL
    },
    **{
        f"pathology.{axis}.frozen_additive_slope_huber": ("pathology", "lower")
        for axis in CANONICAL_V4_PATHOLOGY_AXES
    },
    "leakage.full.additive_probe_loss": ("leakage", "lower"),
    "leakage.isolated.additive_probe_loss": ("leakage", "lower"),
}
CANONICAL_V4_METRIC_NAMES = frozenset(CANONICAL_V4_METRIC_REGISTRY)


def canonical_v4_minimum_support(metric_name: str) -> int:
    """Return the predeclared minimum number of informative C10 donors."""

    if metric_name not in CANONICAL_V4_METRIC_NAMES:
        raise ValueError(f"unknown canonical_v4 metric: {metric_name}")
    if metric_name.startswith("nonneuronal.") or metric_name.startswith(
        "pathology."
    ):
        return 6
    return 8


def _is_sha256(value: Any) -> bool:
    digest = str(value).lower()
    return len(digest) == 64 and all(char in "0123456789abcdef" for char in digest)


def canonical_v4_constraint_registry_sha256(
    constraints: Sequence["MetricConstraint"],
    *,
    denominator_sha256: str,
) -> str:
    """Hash the complete predeclared metric/margin/denominator contract."""

    if not _is_sha256(denominator_sha256):
        raise ValueError("canonical_v4 denominator_sha256 must be a SHA256 digest")
    specs = tuple(constraints)
    names = {spec.name for spec in specs}
    if names != CANONICAL_V4_METRIC_NAMES:
        raise ValueError(
            "canonical_v4 constraint metric names must exactly equal the "
            "predeclared registry; missing="
            f"{sorted(CANONICAL_V4_METRIC_NAMES - names)} extra="
            f"{sorted(names - CANONICAL_V4_METRIC_NAMES)}"
        )
    if len(specs) != len(names):
        raise ValueError("canonical_v4 constraint names must be unique")
    rows: list[dict[str, Any]] = []
    for spec in sorted(specs, key=lambda item: item.name):
        expected_category, expected_direction = CANONICAL_V4_METRIC_REGISTRY[
            spec.name
        ]
        if (spec.category, spec.direction) != (
            expected_category,
            expected_direction,
        ):
            raise ValueError(
                f"canonical_v4 constraint metadata mismatch for {spec.name!r}"
            )
        rows.append(
            {
                "name": spec.name,
                "category": spec.category,
                "direction": spec.direction,
                "margin": float(spec.margin),
                "minimum_support": canonical_v4_minimum_support(spec.name),
                "aggregation": "additive_donor_mean",
            }
        )
    return canonical_json_sha256(
        {
            "schema_version": "kmlee_bam.generator_constraint_registry.v4",
            "denominator_sha256": str(denominator_sha256).lower(),
            "constraints": rows,
        }
    )


def canonical_json_sha256(value: Any) -> str:
    """Hash a JSON-compatible value using one stable serialization."""

    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _names_sha256(values: Iterable[str]) -> str:
    return canonical_json_sha256([str(value) for value in values])


@dataclass(frozen=True)
class GeneratorSelectionDonorSplit:
    """Disjoint train-only donors for weights, ranking, and confirmation."""

    weight_donor_ids: tuple[int, ...]
    ranking_donor_ids: tuple[int, ...]
    confirmation_donor_ids: tuple[int, ...]
    manifest: dict[str, Any]


def _three_way_balance_score(
    features: np.ndarray,
    labels: np.ndarray,
) -> float:
    """Balance each partition's means and variances against all train donors."""

    full_mean = features.mean(axis=0)
    full_var = features.var(axis=0)
    losses: list[float] = []
    for label in range(3):
        group = features[labels == label]
        if group.shape[0] < 2:
            return float("inf")
        mean_error = np.mean(np.square(group.mean(axis=0) - full_mean))
        variance_error = np.mean(np.square(group.var(axis=0) - full_var))
        losses.append(float(mean_error + 0.20 * variance_error))
    return float(np.mean(losses))


def _canonical_v4_hard_support_audit(
    dataset: Any,
    partitions: Mapping[str, tuple[int, ...]],
) -> dict[str, Any]:
    """Audit support required before a 44/10/10 split can be sealed."""

    row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
    donor_by_row = np.asarray(dataset.donor_ids, dtype=np.int64)
    pathology = np.asarray(dataset.prism_pathology, dtype=np.float64)
    pathology_valid = np.asarray(dataset.prism_pathology_valid, dtype=bool)
    celltype = np.asarray(dataset.celltype_ids, dtype=np.int64)
    region = np.asarray(dataset.region_ids, dtype=np.int64)
    sex = np.asarray(dataset.sex_ids, dtype=np.int64)
    age = np.asarray(dataset.age_years, dtype=np.float64)
    celltype_names = tuple(str(value) for value in dataset.spec.celltype_vocab)
    glia_ids = tuple(
        celltype_names.index(name) for name in CANONICAL_V4_NON_NEURONAL
    )
    all_donors = tuple(
        int(value) for value in np.unique(donor_by_row[row_idx]).tolist()
    )
    donor_rows = {
        donor: row_idx[donor_by_row[row_idx] == donor] for donor in all_donors
    }
    donor_pathology = np.full((len(all_donors), pathology.shape[1]), np.nan)
    donor_age = np.full(len(all_donors), np.nan)
    donor_sex = np.full(len(all_donors), -1, dtype=np.int64)
    donor_regions: list[set[int]] = []
    donor_region_fraction = np.zeros(
        (len(all_donors), int(np.max(region[row_idx])) + 1), dtype=np.float64
    )
    donor_celltypes: list[set[int]] = []
    position = {donor: index for index, donor in enumerate(all_donors)}
    for donor, rows in donor_rows.items():
        index = position[donor]
        for axis in range(pathology.shape[1]):
            valid = pathology_valid[rows, axis] & np.isfinite(pathology[rows, axis])
            if bool(valid.any()):
                donor_pathology[index, axis] = float(
                    np.median(pathology[rows[valid], axis])
                )
        finite_age = age[rows][np.isfinite(age[rows])]
        if finite_age.size:
            donor_age[index] = float(finite_age[0])
        donor_sex[index] = int(sex[rows][0])
        donor_regions.append(set(int(value) for value in np.unique(region[rows])))
        counts = np.bincount(
            region[rows], minlength=donor_region_fraction.shape[1]
        ).astype(np.float64)
        donor_region_fraction[index] = counts / max(float(counts.sum()), 1.0)
        donor_celltypes.append(set(int(value) for value in np.unique(celltype[rows])))

    full_region = set().union(*donor_regions)
    full_sex = set(int(value) for value in donor_sex if int(value) >= 0)
    pathology_q33 = np.nanquantile(donor_pathology, 1.0 / 3.0, axis=0)
    pathology_q67 = np.nanquantile(donor_pathology, 2.0 / 3.0, axis=0)
    pathology_median = np.nanmedian(donor_pathology, axis=0)
    age_median = float(np.nanmedian(donor_age))
    partition_audits: dict[str, Any] = {}
    all_pass = True
    for label, donor_ids in partitions.items():
        indices = np.asarray([position[int(value)] for value in donor_ids], dtype=np.int64)
        pathology_audit: list[dict[str, Any]] = []
        for axis in range(pathology.shape[1]):
            values = donor_pathology[indices, axis]
            finite = values[np.isfinite(values)]
            axis_pass = bool(
                finite.size >= 2
                and np.unique(finite).size >= 2
                and (finite <= pathology_q33[axis]).any()
                and (finite >= pathology_q67[axis]).any()
            )
            pathology_audit.append(
                {
                    "axis": int(axis),
                    "observed_count": int(finite.size),
                    "unique_count": int(np.unique(finite).size),
                    "low_tier_present": bool((finite <= pathology_q33[axis]).any()),
                    "high_tier_present": bool((finite >= pathology_q67[axis]).any()),
                    "passed": axis_pass,
                }
            )
        sex_classes = set(int(value) for value in donor_sex[indices] if int(value) >= 0)
        region_classes = set().union(*(donor_regions[index] for index in indices))
        glia_support = {
            celltype_names[ct]: int(
                sum(ct in donor_celltypes[index] for index in indices)
            )
            for ct in glia_ids
        }
        design_pathology = np.where(
            np.isfinite(donor_pathology[indices]),
            donor_pathology[indices],
            pathology_median[None, :],
        )
        design_age = np.where(
            np.isfinite(donor_age[indices]), donor_age[indices], age_median
        )
        design = np.column_stack(
            (
                np.ones(indices.size),
                donor_sex[indices].astype(np.float64),
                design_age,
                design_pathology,
                donor_region_fraction[indices, 1:],
            )
        )
        scale = np.std(design[:, 1:], axis=0)
        scaled = np.column_stack(
            (
                design[:, 0],
                (design[:, 1:] - design[:, 1:].mean(axis=0))
                / np.where(scale > 1e-12, scale, 1.0),
            )
        )
        rank = int(np.linalg.matrix_rank(scaled))
        condition = float(np.linalg.cond(scaled))
        design_pass = bool(rank == scaled.shape[1] and condition <= 1.0e6)
        passed = bool(
            all(item["passed"] for item in pathology_audit)
            and (len(full_sex) < 2 or sex_classes == full_sex)
            and region_classes == full_region
            and all(count >= 6 for count in glia_support.values())
            and design_pass
        )
        partition_audits[label] = {
            "donor_count": int(indices.size),
            "pathology": pathology_audit,
            "sex_classes": sorted(sex_classes),
            "region_classes": sorted(region_classes),
            "glia_donor_support": glia_support,
            "design_column_count": int(scaled.shape[1]),
            "design_rank": rank,
            "design_condition_number": condition,
            "passed": passed,
        }
        all_pass = all_pass and passed
    audit: dict[str, Any] = {
        "schema_version": "kmlee_bam.generator_split_hard_support.v4",
        "requirements": {
            "pathology_each_axis_observed_low_high_and_variable": True,
            "all_observed_sex_and_region_classes": True,
            "six_non_neuronal_celltypes_at_least_six_donors": True,
            "design_full_column_rank": True,
            "design_condition_number_max": 1.0e6,
        },
        "partitions": partition_audits,
        "passed": bool(all_pass),
    }
    audit["audit_sha256"] = canonical_json_sha256(audit)
    return audit


def make_generator_selection_donor_split(
    dataset: Any,
    *,
    weight_donor_count: int = 44,
    ranking_donor_count: int = 10,
    confirmation_donor_count: int = 10,
    seed: int = 420805,
    search_trials: int = 20_000,
) -> GeneratorSelectionDonorSplit:
    """Make a deterministic, balanced 3-way partition of train donors.

    The defaults intentionally sum to the current 64-donor training cohort.
    For another cohort size the caller must state all three counts explicitly;
    silently dropping or duplicating donors is prohibited.
    """

    counts = (
        int(weight_donor_count),
        int(ranking_donor_count),
        int(confirmation_donor_count),
    )
    if any(count < 2 for count in counts):
        raise ValueError("every donor partition must contain at least two donors")
    if int(search_trials) <= 0:
        raise ValueError("search_trials must be positive")

    row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
    donor_by_row = np.asarray(dataset.donor_ids, dtype=np.int64)
    donor_ids = np.unique(donor_by_row[row_idx]).astype(np.int64)
    n_donor = int(donor_ids.size)
    if sum(counts) != n_donor:
        raise ValueError(
            "weight/ranking/confirmation counts must sum exactly to the "
            f"train donor count ({sum(counts)} != {n_donor})"
        )

    features, feature_names = _build_donor_feature_matrix(dataset, donor_ids)
    rng = np.random.default_rng(int(seed))
    n_weight, n_ranking, _ = counts
    canonical_v4 = counts == (44, 10, 10) and n_donor == 64
    best_labels: np.ndarray | None = None
    best_support_audit: dict[str, Any] | None = None
    best_score = float("inf")

    def labels_to_partitions(labels: np.ndarray) -> tuple[tuple[int, ...], ...]:
        return tuple(
            tuple(int(value) for value in donor_ids[labels == label].tolist())
            for label in range(3)
        )

    def support_audit(labels: np.ndarray) -> dict[str, Any] | None:
        if not canonical_v4:
            return None
        weight, ranking, confirmation = labels_to_partitions(labels)
        return _canonical_v4_hard_support_audit(
            dataset,
            {"W44": weight, "R10": ranking, "C10": confirmation},
        )

    for _ in range(int(search_trials)):
        permutation = rng.permutation(n_donor)
        labels = np.full(n_donor, 2, dtype=np.int8)
        labels[permutation[:n_weight]] = 0
        labels[permutation[n_weight : n_weight + n_ranking]] = 1
        score = _three_way_balance_score(features, labels)
        if score + 1e-15 < best_score:
            candidate_support = support_audit(labels)
            if (
                candidate_support is not None
                and candidate_support["passed"] is not True
            ):
                continue
            best_score = score
            best_labels = labels.copy()
            best_support_audit = candidate_support
    if best_labels is None:
        if canonical_v4:
            raise RuntimeError(
                "three-way donor search produced no candidate satisfying the "
                "canonical_v4 hard-support contract"
            )
        raise RuntimeError("three-way donor search produced no candidate")

    # Deterministic pairwise swaps preserve exact partition sizes.
    improved = True
    while improved:
        improved = False
        for left_label in range(3):
            left = np.flatnonzero(best_labels == left_label)
            for right_label in range(left_label + 1, 3):
                right = np.flatnonzero(best_labels == right_label)
                for left_index in left.tolist():
                    for right_index in right.tolist():
                        candidate = best_labels.copy()
                        candidate[left_index] = right_label
                        candidate[right_index] = left_label
                        score = _three_way_balance_score(features, candidate)
                        if score + 1e-12 < best_score:
                            candidate_support = support_audit(candidate)
                            if (
                                candidate_support is not None
                                and candidate_support["passed"] is not True
                            ):
                                continue
                            best_score = score
                            best_labels = candidate
                            best_support_audit = candidate_support
                            improved = True
                            break
                    if improved:
                        break
                if improved:
                    break
            if improved:
                break

    partitions = labels_to_partitions(best_labels)
    weight_ids, ranking_ids, confirmation_ids = partitions
    combined = weight_ids + ranking_ids + confirmation_ids
    if len(combined) != n_donor or set(combined) != set(donor_ids.tolist()):
        raise RuntimeError("internal error: donor partition is not complete")
    if len(set(combined)) != n_donor:
        raise RuntimeError("internal error: donor partitions overlap")

    donor_names = np.asarray(dataset.donor_vocab, dtype=object).astype(str)

    def names(ids: tuple[int, ...]) -> tuple[str, ...]:
        return tuple(donor_names[list(ids)].tolist())

    weight_names = names(weight_ids)
    ranking_names = names(ranking_ids)
    confirmation_names = names(confirmation_ids)
    train_names = names(tuple(int(value) for value in donor_ids.tolist()))
    manifest: dict[str, Any] = {
        "schema_version": SPLIT_SCHEMA,
        "source_split": "train_only",
        "seed": int(seed),
        "search_trials": int(search_trials),
        "balance_score": float(best_score),
        "balance_feature_names": list(feature_names),
        "train_donor_count": n_donor,
        "weight_donor_count": len(weight_ids),
        "ranking_donor_count": len(ranking_ids),
        "confirmation_donor_count": len(confirmation_ids),
        "weight_donor_ids": list(weight_ids),
        "ranking_donor_ids": list(ranking_ids),
        "confirmation_donor_ids": list(confirmation_ids),
        "weight_donor_names": list(weight_names),
        "ranking_donor_names": list(ranking_names),
        "confirmation_donor_names": list(confirmation_names),
        "train_donor_sha256": _names_sha256(train_names),
        "weight_donor_sha256": _names_sha256(weight_names),
        "ranking_donor_sha256": _names_sha256(ranking_names),
        "confirmation_donor_sha256": _names_sha256(confirmation_names),
        "disjoint": True,
        "complete": True,
        "official_validation_used": False,
        "official_test_used": False,
    }
    if canonical_v4:
        if best_support_audit is None or best_support_audit["passed"] is not True:
            raise RuntimeError("canonical_v4 search lost its hard-support evidence")
        manifest["hard_support_passed"] = True
        manifest["hard_support_audit"] = best_support_audit
    manifest["manifest_sha256"] = canonical_json_sha256(manifest)
    return GeneratorSelectionDonorSplit(
        weight_donor_ids=weight_ids,
        ranking_donor_ids=ranking_ids,
        confirmation_donor_ids=confirmation_ids,
        manifest=manifest,
    )


def validate_canonical_v4_donor_split_manifest(
    manifest: Mapping[str, Any],
    *,
    expected_train_donor_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Validate the immutable 44/10/10 W/R/C train-only split contract."""

    result = dict(manifest)
    if result.get("schema_version") != SPLIT_SCHEMA:
        raise ValueError("canonical_v4 donor split manifest schema mismatch")
    expected_flags = {
        "source_split": "train_only",
        "disjoint": True,
        "complete": True,
        "official_validation_used": False,
        "official_test_used": False,
        "train_donor_count": 64,
        "weight_donor_count": 44,
        "ranking_donor_count": 10,
        "confirmation_donor_count": 10,
    }
    for name, expected in expected_flags.items():
        if result.get(name) != expected:
            raise ValueError(
                f"canonical_v4 donor split requires {name}={expected!r}"
            )
    if result.get("hard_support_passed") is not True:
        raise ValueError("canonical_v4 donor split lacks passed hard-support audit")
    support = result.get("hard_support_audit")
    if not isinstance(support, Mapping):
        raise ValueError("canonical_v4 donor split lacks hard-support details")
    support = dict(support)
    if support.get("schema_version") != "kmlee_bam.generator_split_hard_support.v4":
        raise ValueError("canonical_v4 hard-support audit schema mismatch")
    if support.get("passed") is not True:
        raise ValueError("canonical_v4 hard-support audit did not pass")
    required_support_contract = {
        "pathology_each_axis_observed_low_high_and_variable": True,
        "all_observed_sex_and_region_classes": True,
        "six_non_neuronal_celltypes_at_least_six_donors": True,
        "design_full_column_rank": True,
        "design_condition_number_max": 1.0e6,
    }
    if support.get("requirements") != required_support_contract:
        raise ValueError("canonical_v4 hard-support requirements mismatch")
    partition_support = support.get("partitions")
    if not isinstance(partition_support, Mapping) or set(partition_support) != {
        "W44",
        "R10",
        "C10",
    }:
        raise ValueError("canonical_v4 hard-support partitions are incomplete")
    if any(
        not isinstance(partition_support[name], Mapping)
        or partition_support[name].get("passed") is not True
        for name in ("W44", "R10", "C10")
    ):
        raise ValueError("canonical_v4 partition hard support did not pass")
    expected_support_sha = canonical_json_sha256(
        {key: value for key, value in support.items() if key != "audit_sha256"}
    )
    if support.get("audit_sha256") != expected_support_sha:
        raise ValueError("canonical_v4 hard-support audit SHA256 mismatch")
    partitions: list[tuple[int, ...]] = []
    name_partitions: list[tuple[str, ...]] = []
    for prefix, expected_count in (("weight", 44), ("ranking", 10), ("confirmation", 10)):
        ids = tuple(result.get(f"{prefix}_donor_ids", ()))
        names = tuple(str(value) for value in result.get(f"{prefix}_donor_names", ()))
        if len(ids) != expected_count or len(names) != expected_count:
            raise ValueError(f"canonical_v4 {prefix} partition has wrong size")
        if any(isinstance(value, bool) or not isinstance(value, int) for value in ids):
            raise ValueError(f"canonical_v4 {prefix} donor ids must be integers")
        if len(set(ids)) != expected_count or len(set(names)) != expected_count:
            raise ValueError(f"canonical_v4 {prefix} partition contains duplicates")
        if result.get(f"{prefix}_donor_sha256") != _names_sha256(names):
            raise ValueError(f"canonical_v4 {prefix} donor SHA256 mismatch")
        partitions.append(ids)
        name_partitions.append(names)
    combined_ids = tuple(value for group in partitions for value in group)
    combined_names = tuple(value for group in name_partitions for value in group)
    if len(set(combined_ids)) != 64 or len(set(combined_names)) != 64:
        raise ValueError("canonical_v4 W/R/C donor partitions overlap")
    if expected_train_donor_names is not None:
        expected_names = tuple(str(value) for value in expected_train_donor_names)
        if len(expected_names) != 64 or set(expected_names) != set(combined_names):
            raise ValueError("canonical_v4 split does not match the train donors")
        # Generator split manifests hash train names in dataset donor-id order.
        if result.get("train_donor_sha256") != _names_sha256(expected_names):
            raise ValueError("canonical_v4 train donor SHA256 mismatch")
    elif not _is_sha256(result.get("train_donor_sha256")):
        raise ValueError("canonical_v4 train donor SHA256 is missing or invalid")
    supplied_sha = result.get("manifest_sha256")
    expected_sha = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "manifest_sha256"}
    )
    if supplied_sha != expected_sha:
        raise ValueError("canonical_v4 donor split manifest SHA256 mismatch")
    return result


def _validate_generator_ids(
    values: Iterable[int],
    *,
    n_generators: int,
    label: str,
) -> tuple[int, ...]:
    raw = tuple(values)
    if any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        for value in raw
    ):
        raise ValueError(f"{label} must contain integer generator ids")
    parsed = tuple(int(value) for value in raw)
    if len(set(parsed)) != len(parsed):
        raise ValueError(f"{label} must be unique")
    if any(value < 0 or value >= int(n_generators) for value in parsed):
        raise ValueError(f"{label} contains an out-of-range generator id")
    return parsed


def nested_topk_generator_ids(
    ranking_scores: Sequence[float] | np.ndarray,
    requested_counts: Iterable[int],
    *,
    protected_generator_ids: Iterable[int] = (),
) -> dict[int, tuple[int, ...]]:
    """Return deterministic, nested, exact-size top-K generator sets."""

    scores = np.asarray(ranking_scores, dtype=np.float64)
    if scores.ndim != 1 or scores.size == 0:
        raise ValueError("ranking_scores must be a non-empty one-dimensional array")
    if not bool(np.isfinite(scores).all()):
        raise ValueError("ranking_scores must be finite")
    n_generators = int(scores.size)
    protected = _validate_generator_ids(
        protected_generator_ids,
        n_generators=n_generators,
        label="protected_generator_ids",
    )
    counts = tuple(sorted({int(value) for value in requested_counts}))
    if not counts:
        raise ValueError("requested_counts must be non-empty")
    if any(value < len(protected) or value <= 0 for value in counts):
        raise ValueError("every requested K must include all protected generators")
    if any(value > n_generators for value in counts):
        raise ValueError("requested K exceeds the generator count")

    protected_set = set(protected)
    local_ids = np.arange(n_generators, dtype=np.int64)
    nonprotected = local_ids[
        np.asarray(
            [value not in protected_set for value in local_ids.tolist()],
            dtype=bool,
        )
    ]
    # lexsort's last key is primary: descending score, then ascending local id.
    order = np.lexsort((nonprotected, -scores[nonprotected]))
    ranked_nonprotected = tuple(int(value) for value in nonprotected[order])

    result: dict[int, tuple[int, ...]] = {}
    previous: set[int] = set()
    for count in counts:
        selected = set(protected)
        selected.update(ranked_nonprotected[: count - len(protected)])
        if len(selected) != count:
            raise RuntimeError("internal error: top-K mask has the wrong size")
        if previous and not previous.issubset(selected):
            raise RuntimeError("internal error: top-K masks are not nested")
        result[count] = tuple(sorted(selected))
        previous = selected
    return result


def build_consensus_generator_ranking(
    rankings: Sequence[Mapping[str, Any]],
    *,
    ranking_result_sha256: Sequence[str],
    topk_counts: Iterable[int] = (128, 192, 256, 320, 384),
    minimum_spearman: float = 0.80,
    minimum_topk_jaccard: float = 0.60,
) -> dict[str, Any]:
    """Build a fail-closed median normalized-rank consensus from two seeds."""

    items = tuple(dict(item) for item in rankings)
    digests = tuple(str(value).lower() for value in ranking_result_sha256)
    if len(items) != 2 or len(digests) != 2:
        raise ValueError("canonical_v4 consensus requires exactly two rankings")
    if any(not _is_sha256(value) for value in digests) or digests[0] == digests[1]:
        raise ValueError("ranking result SHA256 values must be distinct and valid")
    for item in items:
        required = {
            "schema_version": RESULT_SCHEMA,
            "selection_stage": "ranking_only",
            "canonical_mask_allowed": False,
            "exact_confirmation_required": True,
            "canonical_v4": True,
            "stable": True,
            "ranking_donor_partition": "R10",
            "lineage_aggregation": "neuronal_50_non_neuronal_50",
            "official_validation_used": False,
            "official_test_used": False,
        }
        for name, expected in required.items():
            if item.get(name) != expected:
                raise ValueError(f"consensus ranking requires {name}={expected!r}")
    split_hashes = {item.get("donor_split_manifest_sha256") for item in items}
    if len(split_hashes) != 1 or not _is_sha256(next(iter(split_hashes))):
        raise ValueError("rankings do not use the same sealed donor split")
    checkpoint_hashes = {item.get("source_checkpoint_sha256") for item in items}
    registry_hashes = {item.get("registry_sha256") for item in items}
    if len(checkpoint_hashes) != 1 or not _is_sha256(next(iter(checkpoint_hashes))):
        raise ValueError("rankings do not use the same source checkpoint")
    if len(registry_hashes) != 1 or not _is_sha256(next(iter(registry_hashes))):
        raise ValueError("rankings do not use the same generator registry")
    artifact_provenance = items[0].get("artifact_provenance")
    if not isinstance(artifact_provenance, Mapping) or dict(
        items[1].get("artifact_provenance", {})
    ) != dict(artifact_provenance):
        raise ValueError("rankings do not use identical sealed W44 artifacts")
    seeds = tuple(item.get("ranking_seed") for item in items)
    if any(isinstance(seed, bool) or not isinstance(seed, int) for seed in seeds):
        raise ValueError("rankings must declare integer ranking_seed")
    if seeds[0] == seeds[1]:
        raise ValueError("consensus rankings must use independent seeds")
    candidate_counts = {item.get("candidate_count") for item in items}
    if len(candidate_counts) != 1:
        raise ValueError("ranking candidate counts differ")
    n_generator = int(next(iter(candidate_counts)))
    if n_generator <= 1:
        raise ValueError("consensus requires at least two generators")
    registry_indices = tuple(items[0].get("candidate_registry_indices", ()))
    registry_names = tuple(items[0].get("candidate_registry_module_names", ()))
    protected = tuple(items[0].get("protected_local_generator_ids", ()))
    for item in items[1:]:
        if tuple(item.get("candidate_registry_indices", ())) != registry_indices:
            raise ValueError("ranking registry indices differ")
        if tuple(item.get("candidate_registry_module_names", ())) != registry_names:
            raise ValueError("ranking registry names differ")
        if tuple(item.get("protected_local_generator_ids", ())) != protected:
            raise ValueError("ranking protected-generator sets differ")
    if len(registry_indices) != n_generator or len(registry_names) != n_generator:
        raise ValueError("ranking registry mapping is incomplete")

    normalized_ranks: list[np.ndarray] = []
    ordered_ids: list[np.ndarray] = []
    for item in items:
        scores = np.asarray(item.get("ranking_scores"), dtype=np.float64)
        if scores.shape != (n_generator,) or not bool(np.isfinite(scores).all()):
            raise ValueError("ranking scores are incomplete/non-finite")
        order = np.lexsort((np.arange(n_generator), -scores))
        rank = np.empty(n_generator, dtype=np.float64)
        rank[order] = np.arange(n_generator, dtype=np.float64)
        normalized_ranks.append(1.0 - rank / float(n_generator - 1))
        ordered_ids.append(order)
    spearman = float(np.corrcoef(normalized_ranks[0], normalized_ranks[1])[0, 1])
    if not np.isfinite(spearman):
        raise ValueError("ranking Spearman correlation is non-finite")
    counts = tuple(
        sorted(
            {
                int(value)
                for value in topk_counts
                if 0 < int(value) < n_generator
            }
        )
    )
    if not counts:
        raise ValueError("topk_counts contains no proper generator subset")
    jaccard: dict[str, float] = {}
    for count in counts:
        left = set(int(value) for value in ordered_ids[0][:count])
        right = set(int(value) for value in ordered_ids[1][:count])
        jaccard[str(count)] = len(left & right) / len(left | right)
    stable = bool(
        spearman >= float(minimum_spearman)
        and all(
            value >= float(minimum_topk_jaccard) for value in jaccard.values()
        )
    )
    consensus = np.median(np.stack(normalized_ranks, axis=0), axis=0)
    result: dict[str, Any] = {
        "schema_version": RESULT_SCHEMA,
        "consensus_schema_version": CONSENSUS_SCHEMA,
        "selection_stage": "ranking_only",
        "stable": stable,
        "canonical_mask_allowed": False,
        "exact_confirmation_required": True,
        "candidate_count": n_generator,
        "ranking_method": "median_normalized_rank_two_independent_seeds",
        "ranking_score_direction": "higher_is_more_important",
        "ranking_scores": [float(value) for value in consensus.tolist()],
        "ranking_seeds": list(seeds),
        "ranking_result_sha256": list(digests),
        "donor_split_manifest_sha256": next(iter(split_hashes)),
        "source_checkpoint_sha256": next(iter(checkpoint_hashes)),
        "registry_sha256": next(iter(registry_hashes)),
        "artifact_provenance": dict(artifact_provenance),
        "candidate_registry_indices": list(registry_indices),
        "candidate_registry_module_names": list(registry_names),
        "protected_local_generator_ids": list(protected),
        "spearman": spearman,
        "minimum_spearman": float(minimum_spearman),
        "topk_jaccard": jaccard,
        "minimum_topk_jaccard": float(minimum_topk_jaccard),
        "stability_gate_passed": stable,
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["manifest_sha256"] = canonical_json_sha256(result)
    return result


def build_no_compression_certificate(
    consensus_ranking: Mapping[str, Any],
    *,
    r10_feasibility_by_k: Mapping[int, bool],
    r10_evaluation_sha256_by_k: Mapping[int, str],
) -> dict[str, Any]:
    """Certify K=N only when the exact nested K=N-1 boundary fails on R10."""

    consensus = dict(consensus_ranking)
    if consensus.get("schema_version") != RESULT_SCHEMA or consensus.get(
        "consensus_schema_version"
    ) != CONSENSUS_SCHEMA:
        raise ValueError("no-compression certificate requires consensus schema")
    if consensus.get("stability_gate_passed") is not True:
        raise ValueError("no-compression certificate requires stable consensus")
    supplied_consensus_sha = consensus.get("manifest_sha256")
    expected_consensus_sha = canonical_json_sha256(
        {key: value for key, value in consensus.items() if key != "manifest_sha256"}
    )
    if supplied_consensus_sha != expected_consensus_sha:
        raise ValueError("consensus manifest SHA256 mismatch")
    candidate_count = int(consensus["candidate_count"])
    feasibility = {int(key): bool(value) for key, value in r10_feasibility_by_k.items()}
    if feasibility.get(candidate_count) is not True or feasibility.get(
        candidate_count - 1
    ) is not False:
        raise ValueError(
            "K=N no-compression requires exact R10 evidence that K=N passes "
            "and nested K=N-1 fails"
        )
    hashes = {int(key): str(value).lower() for key, value in r10_evaluation_sha256_by_k.items()}
    if set(hashes) != {candidate_count - 1, candidate_count} or any(
        not _is_sha256(value) for value in hashes.values()
    ):
        raise ValueError("no-compression boundary evaluation SHA256 set is invalid")
    result: dict[str, Any] = {
        "schema_version": NO_COMPRESSION_SCHEMA,
        "status": "no_compression_certified",
        "candidate_count": candidate_count,
        "selected_generator_count": candidate_count,
        "no_compression_certificate_present": True,
        "boundary_rule": "nested_K_N_minus_1_fails_and_K_N_passes_on_R10",
        "r10_feasibility_by_k": {
            str(key): feasibility[key] for key in sorted(feasibility)
        },
        "r10_boundary_evaluation_sha256_by_k": {
            str(key): hashes[key] for key in sorted(hashes)
        },
        "consensus_manifest_sha256": supplied_consensus_sha,
        "donor_split_manifest_sha256": consensus["donor_split_manifest_sha256"],
        "confirmation_donors_opened": False,
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["certificate_sha256"] = canonical_json_sha256(result)
    return result


def coarse_generator_counts(
    n_generators: int,
    *,
    protected_count: int = 0,
    requested: Iterable[int] = (128, 192, 256, 320, 384, 414),
) -> tuple[int, ...]:
    """Clip and deduplicate the standard coarse K grid."""

    n_generators = int(n_generators)
    minimum = max(1, int(protected_count))
    if n_generators <= 0 or minimum > n_generators:
        raise ValueError("invalid generator/protected count")
    raw = tuple(int(value) for value in requested)
    if not raw:
        raise ValueError("requested coarse grid must be non-empty")
    clipped = {min(n_generators, max(minimum, value)) for value in raw}
    clipped.add(n_generators)
    return tuple(sorted(clipped))


@dataclass(frozen=True)
class BoundarySearchResult:
    selected_k: int
    feasibility_by_k: dict[int, bool]
    stages: tuple[dict[str, Any], ...]
    monotonicity_assumption: str


def _assert_observed_monotonic(feasibility: Mapping[int, bool]) -> None:
    seen_feasible = False
    for count in sorted(feasibility):
        if bool(feasibility[count]):
            seen_feasible = True
        elif seen_feasible:
            raise ValueError(
                "observed feasibility is non-monotonic across nested masks; "
                "automatic boundary refinement is unsafe"
            )


def refine_smallest_feasible_count(
    evaluate_feasibility: Callable[[int], bool],
    *,
    n_generators: int,
    protected_count: int = 0,
    coarse_counts: Iterable[int] | None = None,
    refinement_steps: Iterable[int] = (16, 4, 1),
) -> BoundarySearchResult:
    """Evaluate a coarse grid, then refine its feasible boundary to one K.

    The refinement is valid only under the audited assumption that feasibility
    is monotone non-decreasing as generators are added to the nested masks.
    Any observed counterexample raises instead of releasing a boundary.
    """

    minimum = max(1, int(protected_count))
    if minimum > int(n_generators):
        raise ValueError("protected_count exceeds n_generators")
    coarse = (
        coarse_generator_counts(
            n_generators,
            protected_count=protected_count,
        )
        if coarse_counts is None
        else tuple(
            sorted(
                {
                    min(int(n_generators), max(minimum, int(value)))
                    for value in coarse_counts
                }
                | {int(n_generators)}
            )
        )
    )
    steps = tuple(int(value) for value in refinement_steps)
    if not steps or any(value <= 0 for value in steps):
        raise ValueError("refinement_steps must be positive")
    if any(left <= right for left, right in zip(steps, steps[1:])):
        raise ValueError("refinement_steps must be strictly decreasing")
    if steps[-1] != 1:
        raise ValueError("the final refinement step must be one")

    cache: dict[int, bool] = {}
    stages: list[dict[str, Any]] = []

    def evaluate(counts: Iterable[int], stage: str) -> None:
        evaluated: list[int] = []
        for count in sorted(set(int(value) for value in counts)):
            if count < minimum or count > int(n_generators) or count in cache:
                continue
            cache[count] = bool(evaluate_feasibility(count))
            evaluated.append(count)
        _assert_observed_monotonic(cache)
        stages.append(
            {
                "stage": stage,
                "counts": evaluated,
                "feasible": [count for count in evaluated if cache[count]],
            }
        )

    evaluate(coarse, "coarse")
    if not any(cache.values()):
        raise RuntimeError("no coarse generator count passed exact confirmation")

    for step in steps:
        upper = min(count for count, feasible in cache.items() if feasible)
        lower_candidates = [
            count for count, feasible in cache.items() if count < upper and not feasible
        ]
        lower = max(lower_candidates, default=minimum - 1)
        # Anchor to the known feasible upper bound so every stage tightens the
        # same bracket even when it is not divisible by ``step``.
        candidates: list[int] = []
        value = upper - step
        while value > lower:
            candidates.append(value)
            value -= step
        evaluate(candidates, f"refine_step_{step}")

    selected = min(count for count, feasible in cache.items() if feasible)
    # The step-1 stage must leave an explicitly observed infeasible predecessor
    # (or reach the allowed minimum); otherwise the boundary is not exact.
    if selected > minimum and cache.get(selected - 1) is not False:
        raise RuntimeError("unit refinement did not certify the exact K boundary")
    return BoundarySearchResult(
        selected_k=selected,
        feasibility_by_k=dict(sorted(cache.items())),
        stages=tuple(stages),
        monotonicity_assumption=(
            "For deterministic nested top-K masks, exact-confirmation feasibility "
            "is monotone non-decreasing in K; observed violations fail closed."
        ),
    )


@dataclass(frozen=True)
class MetricConstraint:
    """One allowable degradation relative to the all-on reference."""

    name: str
    category: str
    direction: str
    margin: float

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("constraint name must be non-empty")
        if self.category not in REQUIRED_CONSTRAINT_CATEGORIES:
            raise ValueError(f"unsupported constraint category: {self.category}")
        if self.direction not in {"higher", "lower"}:
            raise ValueError("constraint direction must be 'higher' or 'lower'")
        if not np.isfinite(self.margin) or float(self.margin) < 0.0:
            raise ValueError("constraint margin must be finite and non-negative")


@dataclass(frozen=True)
class ConstraintConfirmation:
    name: str
    category: str
    direction: str
    margin: float
    candidate_mean: float
    all_on_mean: float
    observed_harm: float
    bootstrap_upper_95: float
    bootstrap_standard_error: float
    passed: bool
    support_count: int | None = None
    minimum_support: int | None = None
    sufficient_power: bool = True

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "category": self.category,
            "direction": self.direction,
            "margin": self.margin,
            "candidate_mean": self.candidate_mean,
            "all_on_mean": self.all_on_mean,
            "observed_harm": self.observed_harm,
            "bootstrap_upper_95": self.bootstrap_upper_95,
            "bootstrap_standard_error": self.bootstrap_standard_error,
            "passed": self.passed,
            "support_count": self.support_count,
            "minimum_support": self.minimum_support,
            "sufficient_power": self.sufficient_power,
        }


@dataclass(frozen=True)
class ExactCandidateConfirmation:
    generator_count: int
    selected_generator_ids: tuple[int, ...]
    constraints: tuple[ConstraintConfirmation, ...]
    primary_score_mean: float
    primary_score_standard_error: float
    primary_score_direction: str
    donor_ids: tuple[str, ...]
    bootstrap_seed: int
    bootstrap_replicates: int
    bootstrap_indices_sha256: str
    constraint_registry_sha256: str | None = None
    denominator_sha256: str | None = None
    metric_aggregation: str = "legacy_unspecified"
    confirmation_role: str = "discovery_only"

    @property
    def insufficient_power_metrics(self) -> tuple[str, ...]:
        return tuple(
            constraint.name
            for constraint in self.constraints
            if not constraint.sufficient_power
        )

    @property
    def status(self) -> str:
        if self.insufficient_power_metrics:
            return "insufficient_power"
        return "passed" if self.feasible else "constraint_failed"

    @property
    def feasible(self) -> bool:
        categories = {constraint.category for constraint in self.constraints}
        return (
            categories.issuperset(REQUIRED_CONSTRAINT_CATEGORIES)
            and bool(self.constraints)
            and all(constraint.passed for constraint in self.constraints)
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "generator_count": self.generator_count,
            "selected_generator_ids": list(self.selected_generator_ids),
            "feasible": self.feasible,
            "status": self.status,
            "insufficient_power_metrics": list(self.insufficient_power_metrics),
            "primary_score_mean": self.primary_score_mean,
            "primary_score_standard_error": self.primary_score_standard_error,
            "primary_score_direction": self.primary_score_direction,
            "confirmation_donor_ids": list(self.donor_ids),
            "bootstrap_seed": self.bootstrap_seed,
            "bootstrap_replicates": self.bootstrap_replicates,
            "bootstrap_indices_sha256": self.bootstrap_indices_sha256,
            "constraint_registry_sha256": self.constraint_registry_sha256,
            "denominator_sha256": self.denominator_sha256,
            "metric_aggregation": self.metric_aggregation,
            "confirmation_role": self.confirmation_role,
            "constraints": [constraint.as_dict() for constraint in self.constraints],
        }


def _aligned_metric_arrays(
    metrics: Mapping[str, Sequence[float] | np.ndarray],
    *,
    expected_names: set[str],
    n_donor: int,
    label: str,
) -> dict[str, np.ndarray]:
    names = set(metrics)
    if names != expected_names:
        raise ValueError(
            f"{label} metric names differ from constraints: "
            f"missing={sorted(expected_names - names)} "
            f"extra={sorted(names - expected_names)}"
        )
    result: dict[str, np.ndarray] = {}
    for name, values in metrics.items():
        array = np.asarray(values, dtype=np.float64)
        if array.shape != (n_donor,):
            raise ValueError(
                f"{label} metric {name!r} must have shape ({n_donor},)"
            )
        if not bool(np.isfinite(array).all()):
            raise ValueError(f"{label} metric {name!r} contains non-finite values")
        result[name] = array
    return result


def confirm_candidate_with_paired_bootstrap(
    *,
    generator_count: int,
    selected_generator_ids: Iterable[int],
    candidate_donor_ids: Sequence[str],
    all_on_donor_ids: Sequence[str],
    candidate_metrics: Mapping[str, Sequence[float] | np.ndarray],
    all_on_metrics: Mapping[str, Sequence[float] | np.ndarray],
    constraints: Sequence[MetricConstraint],
    candidate_primary_scores: Sequence[float] | np.ndarray,
    primary_score_direction: str = "higher",
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 420806,
    canonical_v4: bool = False,
    denominator_sha256: str | None = None,
    confirmation_role: str = "discovery_only",
    metric_support_counts: Mapping[str, int] | None = None,
) -> ExactCandidateConfirmation:
    """Run a shared-index bootstrap over donor-additive metrics.

    ``canonical_v4=True`` deliberately rejects arbitrary/non-linear metric
    collections.  The caller must provide the exact frozen additive registry
    and its denominator manifest hash.  Raw-data callback bootstrapping can be
    added later; this one-shot array API must never claim pooled CCC validity.
    """

    candidate_ids = tuple(str(value) for value in candidate_donor_ids)
    reference_ids = tuple(str(value) for value in all_on_donor_ids)
    if len(candidate_ids) < 2 or candidate_ids != reference_ids:
        raise ValueError(
            "candidate and all-on donor ids must contain at least two donors "
            "and be identically ordered"
        )
    if len(set(candidate_ids)) != len(candidate_ids):
        raise ValueError("confirmation donor ids must be unique")
    if int(bootstrap_replicates) < 100:
        raise ValueError("bootstrap_replicates must be at least 100")
    if primary_score_direction not in {"higher", "lower"}:
        raise ValueError("primary_score_direction must be 'higher' or 'lower'")

    specs = tuple(constraints)
    if not specs:
        raise ValueError("at least one exact confirmation constraint is required")
    if len({spec.name for spec in specs}) != len(specs):
        raise ValueError("constraint names must be unique")
    registry_sha256: str | None = None
    if canonical_v4:
        if confirmation_role not in {"preselected_canonical", "discovery_only"}:
            raise ValueError("unsupported canonical_v4 confirmation_role")
        registry_sha256 = canonical_v4_constraint_registry_sha256(
            specs,
            denominator_sha256=str(denominator_sha256),
        )
        if metric_support_counts is None or set(metric_support_counts) != {
            spec.name for spec in specs
        }:
            raise ValueError(
                "canonical_v4 metric_support_counts must exactly match the registry"
            )
    expected_names = {spec.name for spec in specs}
    n_donor = len(candidate_ids)
    candidate = _aligned_metric_arrays(
        candidate_metrics,
        expected_names=expected_names,
        n_donor=n_donor,
        label="candidate",
    )
    all_on = _aligned_metric_arrays(
        all_on_metrics,
        expected_names=expected_names,
        n_donor=n_donor,
        label="all-on",
    )
    primary = np.asarray(candidate_primary_scores, dtype=np.float64)
    if primary.shape != (n_donor,) or not bool(np.isfinite(primary).all()):
        raise ValueError(
            f"candidate_primary_scores must be finite with shape ({n_donor},)"
        )

    selected_raw = tuple(selected_generator_ids)
    if any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        for value in selected_raw
    ):
        raise ValueError("selected_generator_ids must contain integers")
    selected = tuple(sorted(int(value) for value in selected_raw))
    if len(set(selected)) != len(selected):
        raise ValueError("selected_generator_ids must be unique")
    if len(selected) != int(generator_count) or int(generator_count) <= 0:
        raise ValueError("selected_generator_ids must have exactly generator_count ids")
    if selected[0] < 0:
        raise ValueError("selected generator ids must be non-negative")

    rng = np.random.default_rng(int(bootstrap_seed))
    indices = rng.integers(
        0,
        n_donor,
        size=(int(bootstrap_replicates), n_donor),
        endpoint=False,
        dtype=np.int64,
    )
    index_sha256 = hashlib.sha256(indices.tobytes(order="C")).hexdigest()
    confirmations: list[ConstraintConfirmation] = []
    for spec in specs:
        candidate_values = candidate[spec.name]
        reference_values = all_on[spec.name]
        if spec.direction == "higher":
            paired_harm = reference_values - candidate_values
        else:
            paired_harm = candidate_values - reference_values
        bootstrap_harm = paired_harm[indices].mean(axis=1)
        upper = float(np.quantile(bootstrap_harm, 0.95))
        standard_error = float(np.std(bootstrap_harm, ddof=1))
        support_count: int | None = None
        minimum_support: int | None = None
        sufficient_power = True
        if canonical_v4:
            raw_support = metric_support_counts[spec.name]
            if isinstance(raw_support, (bool, np.bool_)) or not isinstance(
                raw_support, (int, np.integer)
            ):
                raise ValueError("canonical_v4 metric support counts must be integers")
            support_count = int(raw_support)
            if support_count < 0 or support_count > n_donor:
                raise ValueError("canonical_v4 metric support count is out of range")
            minimum_support = canonical_v4_minimum_support(spec.name)
            sufficient_power = support_count >= minimum_support
        confirmations.append(
            ConstraintConfirmation(
                name=spec.name,
                category=spec.category,
                direction=spec.direction,
                margin=float(spec.margin),
                candidate_mean=float(candidate_values.mean()),
                all_on_mean=float(reference_values.mean()),
                observed_harm=float(paired_harm.mean()),
                bootstrap_upper_95=upper,
                bootstrap_standard_error=standard_error,
                passed=bool(
                    upper <= float(spec.margin) and sufficient_power
                ),
                support_count=support_count,
                minimum_support=minimum_support,
                sufficient_power=sufficient_power,
            )
        )

    primary_se = (
        float(np.std(primary, ddof=1) / np.sqrt(n_donor))
        if n_donor > 1
        else 0.0
    )
    return ExactCandidateConfirmation(
        generator_count=int(generator_count),
        selected_generator_ids=selected,
        constraints=tuple(confirmations),
        primary_score_mean=float(primary.mean()),
        primary_score_standard_error=primary_se,
        primary_score_direction=primary_score_direction,
        donor_ids=candidate_ids,
        bootstrap_seed=int(bootstrap_seed),
        bootstrap_replicates=int(bootstrap_replicates),
        bootstrap_indices_sha256=index_sha256,
        constraint_registry_sha256=registry_sha256,
        denominator_sha256=(
            str(denominator_sha256).lower() if canonical_v4 else None
        ),
        metric_aggregation=(
            "additive_donor_mean" if canonical_v4 else "legacy_unspecified"
        ),
        confirmation_role=confirmation_role,
    )


def choose_smallest_confirmed_candidate(
    candidates: Sequence[ExactCandidateConfirmation],
    *,
    one_standard_error_guard: bool = False,
) -> ExactCandidateConfirmation:
    """Discovery-only retrospective selector; never authorizes canonical use.

    Comparing several K values on the same confirmation donors adapts to C and
    invalidates their role as untouched confirmation data.  This helper is
    retained for diagnostics/legacy reports only.
    """

    candidates = tuple(candidates)
    if not candidates:
        raise ValueError("no exact candidate confirmations were supplied")
    counts = [candidate.generator_count for candidate in candidates]
    if len(set(counts)) != len(counts):
        raise ValueError("candidate generator counts must be unique")
    directions = {candidate.primary_score_direction for candidate in candidates}
    if len(directions) != 1:
        raise ValueError("all candidates must use the same primary-score direction")
    feasible = [candidate for candidate in candidates if candidate.feasible]
    if not feasible:
        raise RuntimeError("no candidate passed every exact confirmation constraint")

    eligible = feasible
    if one_standard_error_guard:
        direction = next(iter(directions))
        if direction == "higher":
            best = max(feasible, key=lambda item: item.primary_score_mean)
            threshold = best.primary_score_mean - best.primary_score_standard_error
            eligible = [item for item in feasible if item.primary_score_mean >= threshold]
        else:
            best = min(feasible, key=lambda item: item.primary_score_mean)
            threshold = best.primary_score_mean + best.primary_score_standard_error
            eligible = [item for item in feasible if item.primary_score_mean <= threshold]
    return min(eligible, key=lambda item: item.generator_count)


def _validate_sha256_mapping(source_hashes: Mapping[str, str]) -> dict[str, str]:
    missing = REQUIRED_SOURCE_HASHES - set(source_hashes)
    if missing:
        raise ValueError(f"missing required source hashes: {sorted(missing)}")
    validated: dict[str, str] = {}
    for name, digest in source_hashes.items():
        digest = str(digest).lower()
        if len(digest) != 64 or any(value not in "0123456789abcdef" for value in digest):
            raise ValueError(f"invalid SHA256 digest for {name}")
        validated[str(name)] = digest
    return dict(sorted(validated.items()))


def generator_mask_sha256(
    selected_generator_ids: Iterable[int],
    *,
    candidate_count: int,
) -> str:
    """SHA256 of one uint8 0/1 byte per local candidate, in local-id order."""

    candidate_count = int(candidate_count)
    if candidate_count <= 0:
        raise ValueError("candidate_count must be positive")
    selected = _validate_generator_ids(
        selected_generator_ids,
        n_generators=candidate_count,
        label="selected_generator_ids",
    )
    mask = np.zeros(candidate_count, dtype=np.uint8)
    mask[list(selected)] = 1
    return hashlib.sha256(mask.tobytes(order="C")).hexdigest()


def build_generator_k_preselection_manifest(
    *,
    selected_generator_ids: Iterable[int],
    candidate_count: int,
    ranking_result_sha256: str,
    donor_split_manifest_sha256: str,
) -> dict[str, Any]:
    """Freeze exactly one R-selected K before confirmation donors are opened."""

    if not _is_sha256(ranking_result_sha256) or not _is_sha256(
        donor_split_manifest_sha256
    ):
        raise ValueError("preselection source SHA256 is invalid")
    selected = _validate_generator_ids(
        selected_generator_ids,
        n_generators=int(candidate_count),
        label="selected_generator_ids",
    )
    result: dict[str, Any] = {
        "schema_version": PRESELECTION_SCHEMA,
        "selection_partition": "R10",
        "confirmation_partition_opened": False,
        "prior_confirmation_attempts": 0,
        "failure_policy": "abort_without_trying_another_K",
        "candidate_count": int(candidate_count),
        "preselected_k": len(selected),
        "selected_local_generator_ids": list(selected),
        "mask_sha256": generator_mask_sha256(
            selected, candidate_count=int(candidate_count)
        ),
        "ranking_result_sha256": str(ranking_result_sha256).lower(),
        "donor_split_manifest_sha256": str(
            donor_split_manifest_sha256
        ).lower(),
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["manifest_sha256"] = canonical_json_sha256(result)
    return result


def validate_generator_k_preselection_manifest(
    manifest: Mapping[str, Any],
    *,
    selected_generator_ids: Iterable[int],
    candidate_count: int,
    ranking_result_sha256: str,
    donor_split_manifest_sha256: str,
) -> dict[str, Any]:
    expected = build_generator_k_preselection_manifest(
        selected_generator_ids=selected_generator_ids,
        candidate_count=candidate_count,
        ranking_result_sha256=ranking_result_sha256,
        donor_split_manifest_sha256=donor_split_manifest_sha256,
    )
    supplied = dict(manifest)
    if supplied != expected:
        raise ValueError(
            "preselection manifest mismatch: K/mask must be frozen on R10 "
            "before C10 is opened"
        )
    return supplied


def build_canonical_generator_result(
    selected: ExactCandidateConfirmation,
    *,
    n_generators: int,
    protected_generator_ids: Iterable[int],
    source_hashes: Mapping[str, str],
    donor_split_manifest: Mapping[str, Any],
    preselection_manifest: Mapping[str, Any],
    evaluated_candidates: Sequence[ExactCandidateConfirmation],
    registry_indices: Sequence[int],
    registry_names: Sequence[str],
    ranking_method: str = "hard_concrete_ranking_only",
    one_standard_error_guard: bool = False,
) -> dict[str, Any]:
    """Build a warm-start-compatible result, or fail closed."""

    if selected.metric_aggregation == "additive_donor_mean":
        # Validate the exact registry before the broader legacy category audit
        # so a category-complete but metric-incomplete payload is unambiguous.
        canonical_v4_constraint_registry_sha256(
            tuple(
                MetricConstraint(
                    constraint.name,
                    constraint.category,
                    constraint.direction,
                    constraint.margin,
                )
                for constraint in selected.constraints
            ),
            denominator_sha256=str(selected.denominator_sha256),
        )
    if not selected.feasible:
        raise ValueError("selected candidate did not pass every constraint category")
    # K must be frozen using W/R before C is opened.  Testing several K values
    # on C and choosing the first pass is selection on confirmation data.
    if len(evaluated_candidates) != 1:
        raise ValueError(
            "canonical release requires exactly one preselected K evaluated "
            "once on confirmation donors; failure forbids trying the next K"
        )
    if selected.confirmation_role != "preselected_canonical":
        raise ValueError("canonical release requires a preselected_canonical candidate")
    if selected.metric_aggregation != "additive_donor_mean":
        raise ValueError("canonical release requires additive donor metrics")
    if not _is_sha256(selected.denominator_sha256):
        raise ValueError("canonical release requires denominator_sha256")
    expected_registry_sha = canonical_v4_constraint_registry_sha256(
        tuple(
            MetricConstraint(
                constraint.name,
                constraint.category,
                constraint.direction,
                constraint.margin,
            )
            for constraint in selected.constraints
        ),
        denominator_sha256=str(selected.denominator_sha256),
    )
    if selected.constraint_registry_sha256 != expected_registry_sha:
        raise ValueError("canonical constraint registry SHA256 mismatch")
    evaluated_by_count = {
        candidate.generator_count: candidate for candidate in evaluated_candidates
    }
    if len(evaluated_by_count) != len(evaluated_candidates):
        raise ValueError("evaluated candidate generator counts must be unique")
    if evaluated_by_count.get(selected.generator_count) != selected:
        raise ValueError("selected candidate is absent from evaluated_candidates")
    if one_standard_error_guard:
        raise ValueError(
            "the 1-SE rule is sensitivity-only and cannot govern canonical release"
        )
    n_generators = int(n_generators)
    selected_ids = _validate_generator_ids(
        selected.selected_generator_ids,
        n_generators=n_generators,
        label="selected_generator_ids",
    )
    protected_ids = _validate_generator_ids(
        protected_generator_ids,
        n_generators=n_generators,
        label="protected_generator_ids",
    )
    if not set(protected_ids).issubset(selected_ids):
        raise ValueError("selected candidate omits a protected generator")
    if len(selected_ids) != selected.generator_count:
        raise ValueError("selected generator count does not match selected ids")
    hashes = _validate_sha256_mapping(source_hashes)

    if len(registry_indices) != n_generators or len(registry_names) != n_generators:
        raise ValueError(
            "registry_indices and registry_names must have candidate_count entries"
        )
    if any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        for value in registry_indices
    ):
        raise ValueError("registry_indices must contain integers")
    parsed_registry_indices = tuple(int(value) for value in registry_indices)
    parsed_registry_names = tuple(str(value) for value in registry_names)
    if len(set(parsed_registry_indices)) != n_generators:
        raise ValueError("registry_indices must be unique")
    if any(not value for value in parsed_registry_names):
        raise ValueError("registry_names must be non-empty")

    split = validate_canonical_v4_donor_split_manifest(donor_split_manifest)
    if split.get("disjoint") is not True or split.get("complete") is not True:
        raise ValueError("donor split manifest is not disjoint and complete")
    if split.get("official_validation_used") is not False:
        raise ValueError("donor split manifest used official validation donors")
    if split.get("official_test_used") is not False:
        raise ValueError("donor split manifest used official test donors")
    if hashes["donor_split_manifest_sha256"] != canonical_json_sha256(
        {key: value for key, value in split.items() if key != "manifest_sha256"}
    ):
        raise ValueError("donor split manifest SHA256 mismatch")
    if tuple(split.get("confirmation_donor_names", ())) != selected.donor_ids:
        raise ValueError("confirmation donors differ from the split manifest")
    preselection = validate_generator_k_preselection_manifest(
        preselection_manifest,
        selected_generator_ids=selected_ids,
        candidate_count=n_generators,
        ranking_result_sha256=hashes["ranking_result_sha256"],
        donor_split_manifest_sha256=hashes["donor_split_manifest_sha256"],
    )
    if hashes["preselection_manifest_sha256"] != preselection["manifest_sha256"]:
        raise ValueError("preselection manifest SHA256 mismatch")

    mask_sha256 = generator_mask_sha256(
        selected_ids,
        candidate_count=n_generators,
    )
    selected_registry_indices = [
        parsed_registry_indices[index] for index in selected_ids
    ]
    selected_registry_names = [parsed_registry_names[index] for index in selected_ids]
    sensitivity_selected = choose_smallest_confirmed_candidate(
        tuple(evaluated_by_count.values()),
        one_standard_error_guard=True,
    )

    result: dict[str, Any] = {
        "schema_version": RESULT_SCHEMA,
        "confirmation_schema_version": CONFIRMATION_SCHEMA,
        "status": "exact_confirmation_passed",
        "stable": True,
        "selection_stage": "exact_confirmation",
        "exact_confirmation_passed": True,
        "every_confirmation_constraint_passed": True,
        "ranking_method": str(ranking_method),
        "canonical_selection_rule": "one_preselected_K_single_C_evaluation",
        "confirmation_failure_policy": "abort_without_trying_another_K",
        "one_standard_error_sensitivity": {
            "used_for_canonical_selection": False,
            "selected_generator_count_if_applied": (
                sensitivity_selected.generator_count
            ),
            "differs_from_canonical": bool(sensitivity_selected != selected),
        },
        "nested_mask_feasibility_assumption": (
            "Feasibility is assumed monotone non-decreasing in K for the "
            "deterministic nested top-K masks; any observed violation fails closed."
        ),
        "exact_confirmation_required": False,
        "canonical_mask_allowed": True,
        "candidate_count": n_generators,
        "hard_active": len(selected_ids),
        "pruned_generators": n_generators - len(selected_ids),
        "n_generators": n_generators,
        "selected_local_generator_ids": list(selected_ids),
        "protected_local_generator_ids": list(protected_ids),
        "selected_registry_indices": selected_registry_indices,
        "selected_registry_names": selected_registry_names,
        # Keep the discovery-stage spelling for downstream report compatibility.
        "selected_registry_module_names": selected_registry_names,
        "mask_sha256": mask_sha256,
        "mask_sha256_definition": (
            "SHA256 of candidate_count uint8 bytes in local-generator-id order; "
            "byte is 1 when selected and 0 otherwise"
        ),
        "selected_confirmation": selected.as_dict(),
        "evaluated_generator_counts": sorted(evaluated_by_count),
        "evaluated_candidates": [
            evaluated_by_count[count].as_dict() for count in sorted(evaluated_by_count)
        ],
        "source_hashes": hashes,
        "source_hash_definitions": {
            "ranking_result_sha256": "raw-byte SHA256 of ranking result JSON",
            "confirmation_input_sha256": (
                "raw-byte SHA256 of paired per-donor confirmation input"
            ),
            "source_checkpoint_sha256": "raw-byte SHA256 of source checkpoint",
            "registry_sha256": "raw-byte SHA256 of source registry file",
            "donor_split_manifest_sha256": (
                "canonical JSON SHA256 of donor split manifest excluding its "
                "manifest_sha256 field"
            ),
        },
        "source_checkpoint_sha256": hashes["source_checkpoint_sha256"],
        "registry_sha256": hashes["registry_sha256"],
        "ranking_result_sha256": hashes["ranking_result_sha256"],
        "donor_split_manifest": split,
        "preselection_manifest": preselection,
        "constraint_registry_sha256": expected_registry_sha,
        "denominator_sha256": selected.denominator_sha256,
        "final_fixed_structure_retraining_required": True,
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["canonical_payload_sha256"] = canonical_json_sha256(result)
    return result


def write_canonical_generator_result(
    path: str | Path,
    result: Mapping[str, Any],
) -> Path:
    """Atomically write a previously validated canonical result payload."""

    payload = dict(result)
    if payload.get("schema_version") != RESULT_SCHEMA:
        raise ValueError("canonical result schema mismatch")
    if payload.get("stable") is not True:
        raise ValueError("canonical result requires stable=true")
    if payload.get("canonical_mask_allowed") is not True:
        raise ValueError("canonical result requires canonical_mask_allowed=true")
    if payload.get("selection_stage") != "exact_confirmation":
        raise ValueError("canonical result requires exact_confirmation stage")
    if payload.get("exact_confirmation_passed") is not True:
        raise ValueError("canonical result requires exact confirmation to pass")
    if payload.get("every_confirmation_constraint_passed") is not True:
        raise ValueError("canonical result requires every constraint to pass")
    if payload.get("official_validation_used") is not False:
        raise ValueError("canonical result must not use official validation")
    if payload.get("official_test_used") is not False:
        raise ValueError("canonical result must not use official test")
    candidate_count = payload.get("candidate_count")
    hard_active = payload.get("hard_active")
    selected_ids = payload.get("selected_local_generator_ids")
    if (
        isinstance(candidate_count, bool)
        or not isinstance(candidate_count, int)
        or isinstance(hard_active, bool)
        or not isinstance(hard_active, int)
        or not isinstance(selected_ids, list)
    ):
        raise ValueError("canonical result has invalid mask count fields")
    parsed_ids = _validate_generator_ids(
        selected_ids,
        n_generators=candidate_count,
        label="selected_local_generator_ids",
    )
    if len(parsed_ids) != hard_active:
        raise ValueError("canonical result hard_active does not match selected ids")
    expected_mask_sha256 = generator_mask_sha256(
        parsed_ids,
        candidate_count=candidate_count,
    )
    if payload.get("mask_sha256") != expected_mask_sha256:
        raise ValueError("canonical result mask_sha256 mismatch")
    supplied_payload_sha256 = payload.get("canonical_payload_sha256")
    expected_payload_sha256 = canonical_json_sha256(
        {
            key: value
            for key, value in payload.items()
            if key != "canonical_payload_sha256"
        }
    )
    if supplied_payload_sha256 != expected_payload_sha256:
        raise ValueError("canonical result payload SHA256 mismatch")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(target)
    return target
