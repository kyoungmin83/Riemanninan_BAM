"""One fail-closed donor allow-list shared by canonical PRISM-v4 producers."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import torch

from kmlee_bam.training.architecture_donor_split import DonorIndexView
from kmlee_bam.training.generator_exact_confirmation import (
    canonical_json_sha256,
    validate_canonical_v4_donor_split_manifest,
)


ALLOWLIST_SCHEMA = "kmlee_bam.sealed_donor_allowlist.v4"
AGE_NORMALIZATION_SCHEMA = "kmlee_bam.weight_only_age_normalization.v4"

_PARTITION_PREFIX = {"W44": "weight", "R10": "ranking", "C10": "confirmation"}


@dataclass(frozen=True)
class SealedDonorAllowlist:
    partition: str
    donor_ids: tuple[int, ...]
    donor_names: tuple[str, ...]
    donor_names_sha256: str
    manifest_path: str
    manifest_sha256: str
    manifest: dict[str, Any]

    def evidence(self) -> dict[str, Any]:
        value: dict[str, Any] = {
            "schema_version": ALLOWLIST_SCHEMA,
            "partition": self.partition,
            "donor_ids": list(self.donor_ids),
            "donor_names": list(self.donor_names),
            "donor_names_sha256": self.donor_names_sha256,
            "donor_split_manifest_path": self.manifest_path,
            "donor_split_manifest_sha256": self.manifest_sha256,
            "row_level_allowlist_enforced": True,
            "official_validation_used": False,
            "official_test_used": False,
        }
        value["evidence_sha256"] = canonical_json_sha256(value)
        return value


def _load_json(path: str | Path) -> dict[str, Any]:
    source = Path(path).resolve()
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"failed to read sealed donor manifest: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError("sealed donor manifest must be a JSON object")
    return value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_sealed_donor_allowlist(
    manifest_path: str | Path,
    *,
    partition: str = "W44",
    expected_train_donor_names: Sequence[str] | None = None,
) -> SealedDonorAllowlist:
    """Load one immutable partition; a label without a valid manifest fails."""

    partition = str(partition).upper()
    if partition not in _PARTITION_PREFIX:
        raise ValueError("sealed donor partition must be W44, R10, or C10")
    source = Path(manifest_path).resolve()
    manifest = validate_canonical_v4_donor_split_manifest(
        _load_json(source),
        expected_train_donor_names=expected_train_donor_names,
    )
    prefix = _PARTITION_PREFIX[partition]
    ids = tuple(int(value) for value in manifest[f"{prefix}_donor_ids"])
    names = tuple(str(value) for value in manifest[f"{prefix}_donor_names"])
    expected_count = {"W44": 44, "R10": 10, "C10": 10}[partition]
    if len(ids) != expected_count or len(names) != expected_count:
        raise ValueError(f"sealed {partition} partition has the wrong size")
    return SealedDonorAllowlist(
        partition=partition,
        donor_ids=ids,
        donor_names=names,
        donor_names_sha256=str(manifest[f"{prefix}_donor_sha256"]),
        manifest_path=str(source),
        manifest_sha256=str(manifest["manifest_sha256"]),
        manifest=manifest,
    )


def dataset_train_donor_names(dataset: Any) -> tuple[str, ...]:
    """Return actual selected donor names in global donor-id order."""

    rows = np.asarray(dataset.row_idx, dtype=np.int64)
    donor_by_row = np.asarray(dataset.donor_ids, dtype=np.int64)
    vocab = np.asarray(dataset.donor_vocab, dtype=object).astype(str)
    ids = np.unique(donor_by_row[rows]).astype(np.int64)
    if ids.size == 0 or int(ids.min()) < 0 or int(ids.max()) >= len(vocab):
        raise ValueError("dataset has invalid or empty donor membership")
    return tuple(str(value) for value in vocab[ids].tolist())


def bind_allowlist_to_dataset(
    dataset: Any,
    allowlist: SealedDonorAllowlist,
) -> tuple[int, ...]:
    """Map the sealed names to this dataset and verify every selected row."""

    train_names = dataset_train_donor_names(dataset)
    # Revalidate the split against the actual dataset order, not just itself.
    validate_canonical_v4_donor_split_manifest(
        allowlist.manifest,
        expected_train_donor_names=train_names,
    )
    vocab = np.asarray(dataset.donor_vocab, dtype=object).astype(str)
    mapped_ids: list[int] = []
    for manifest_id, name in zip(allowlist.donor_ids, allowlist.donor_names):
        if manifest_id < 0 or manifest_id >= len(vocab):
            raise ValueError(f"sealed donor id {manifest_id} is outside dataset vocab")
        if str(vocab[manifest_id]) != name:
            raise ValueError(
                "sealed donor id/name pair differs from the dataset: "
                f"{manifest_id}/{name!r} != {vocab[manifest_id]!r}"
            )
        mapped_ids.append(int(manifest_id))
    rows = np.asarray(dataset.row_idx, dtype=np.int64)
    donor_by_row = np.asarray(dataset.donor_ids, dtype=np.int64)
    observed = tuple(
        int(value)
        for value in np.unique(
            donor_by_row[rows][
                np.isin(donor_by_row[rows], np.asarray(mapped_ids, dtype=np.int64))
            ]
        ).tolist()
    )
    if observed != tuple(sorted(mapped_ids)):
        raise ValueError(
            f"dataset has no selected rows for one or more {allowlist.partition} donors"
        )
    return tuple(mapped_ids)


def restrict_dataset_to_allowlist(
    dataset: Any,
    allowlist: SealedDonorAllowlist,
) -> DonorIndexView:
    donor_ids = bind_allowlist_to_dataset(dataset, allowlist)
    view = DonorIndexView(dataset, donor_ids)
    observed_names = dataset_train_donor_names(view)
    if observed_names != allowlist.donor_names:
        raise RuntimeError(
            "row-filtered dataset donor names differ from the sealed partition"
        )
    view.sealed_donor_allowlist = allowlist.evidence()
    view.canonical_w44_source_enabled = allowlist.partition == "W44"
    return view


def refit_age_normalization_for_allowlist(
    dataset: Any,
    allowlist: SealedDonorAllowlist,
) -> dict[str, Any]:
    """Refit age moments on observed ages from only the sealed donor rows."""

    donor_ids = bind_allowlist_to_dataset(dataset, allowlist)
    donor_age = np.asarray(dataset.donor_age_years, dtype=np.float64)
    if max(donor_ids) >= donor_age.size:
        raise ValueError("donor_age_years does not cover the sealed donor ids")
    selected = donor_age[np.asarray(donor_ids, dtype=np.int64)]
    finite = np.isfinite(selected)
    observed = selected[finite]
    observed_names = tuple(
        name for name, keep in zip(allowlist.donor_names, finite.tolist()) if keep
    )
    if observed.size < 2:
        raise ValueError("sealed donor set has fewer than two observed ages")
    mean = float(np.mean(observed))
    standard_deviation = float(np.std(observed))
    if not np.isfinite(standard_deviation) or standard_deviation < 1e-6:
        raise ValueError("sealed donor age standard deviation is degenerate")
    dataset.age_train_mean = mean
    dataset.age_train_std = standard_deviation
    dataset.age_z = np.where(
        np.asarray(dataset.age_valid, dtype=bool),
        (np.asarray(dataset.age_years, dtype=np.float64) - mean)
        / standard_deviation,
        0.0,
    ).astype(np.float32)
    result: dict[str, Any] = {
        "schema_version": AGE_NORMALIZATION_SCHEMA,
        "fit_partition": allowlist.partition,
        "allowed_fit_donor_names": list(allowlist.donor_names),
        "allowed_fit_donor_sha256": allowlist.donor_names_sha256,
        "observed_age_donor_names": list(observed_names),
        "observed_age_donor_count": int(observed.size),
        "age_mean": mean,
        "age_standard_deviation": standard_deviation,
        "donor_split_manifest_sha256": allowlist.manifest_sha256,
        "row_level_allowlist_enforced": True,
        "official_validation_used": False,
        "official_test_used": False,
    }
    result["manifest_sha256"] = canonical_json_sha256(result)
    return result


def validate_train_name_universe(
    donor_values: Iterable[Any],
    allowlist: SealedDonorAllowlist,
) -> tuple[str, ...]:
    """Prove an obs table's train donors equal the sealed W+R+C universe."""

    observed = {str(value) for value in donor_values}
    expected = {
        str(value)
        for prefix in ("weight", "ranking", "confirmation")
        for value in allowlist.manifest[f"{prefix}_donor_names"]
    }
    if len(expected) != 64 or observed != expected:
        raise ValueError(
            "observed train-donor universe differs from the sealed 64 donors"
        )
    return tuple(sorted(observed))


def donor_name_row_mask(
    donor_values: Sequence[Any] | np.ndarray,
    allowlist: SealedDonorAllowlist,
) -> np.ndarray:
    values = np.asarray(donor_values, dtype=object).astype(str)
    mask = np.isin(values, np.asarray(allowlist.donor_names, dtype=str))
    observed = {str(value) for value in values[mask].tolist()}
    if observed != set(allowlist.donor_names):
        raise ValueError(
            f"row table lacks one or more sealed {allowlist.partition} donors"
        )
    return mask


def validate_w44_source_checkpoint(
    checkpoint_path: str | Path,
    allowlist: SealedDonorAllowlist,
) -> tuple[str, dict[str, Any]]:
    """Open a checkpoint and bind its actual optimizer-row evidence to W44."""

    if allowlist.partition != "W44":
        raise ValueError("source checkpoint validation requires the W44 partition")
    source = Path(checkpoint_path).resolve()
    digest = sha256_file(source)
    payload = torch.load(source, map_location="cpu", weights_only=False)
    if not isinstance(payload, Mapping):
        raise ValueError("canonical W44 source checkpoint must be a mapping")
    provenance = payload.get("canonical_w44_source_provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError(
            "checkpoint lacks canonical_w44_source_provenance"
        )
    result = dict(provenance)
    required = {
        "schema_version": "kmlee_bam.canonical_w44_source_checkpoint.v4",
        "fit_partition": "W44",
        "observed_optimizer_donor_names": list(allowlist.donor_names),
        "observed_optimizer_donor_count": 44,
        "allowed_fit_donor_names": list(allowlist.donor_names),
        "allowed_fit_donor_sha256": allowlist.donor_names_sha256,
        "donor_split_manifest_sha256": allowlist.manifest_sha256,
        "row_level_allowlist_enforced": True,
        "official_validation_used": False,
        "official_test_used": False,
    }
    for name, expected in required.items():
        if result.get(name) != expected:
            raise ValueError(
                f"canonical W44 checkpoint requires {name}={expected!r}"
            )
    age_sha = str(result.get("age_normalization_manifest_sha256", "")).lower()
    if len(age_sha) != 64 or any(
        character not in "0123456789abcdef" for character in age_sha
    ):
        raise ValueError(
            "canonical W44 checkpoint lacks age-normalization provenance"
        )
    supplied = result.get("manifest_sha256")
    expected_manifest_sha = canonical_json_sha256(
        {key: value for key, value in result.items() if key != "manifest_sha256"}
    )
    if supplied != expected_manifest_sha:
        raise ValueError("canonical W44 checkpoint provenance SHA256 mismatch")
    return digest, result
