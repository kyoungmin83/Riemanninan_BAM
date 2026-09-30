"""Leakage-audited reliability artifacts for PRISM module-local personal paths.

The artifact is deliberately small and immutable: it contains one split-half
reliability value per source context and biological module, plus enough
provenance to prove that only the sealed fitting-donor partition was used.
Loading is fail-closed because a context/module order mismatch would silently
route reliability weights to the wrong biological programme.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np


MODULE_LOCAL_RELIABILITY_SCHEMA = "kmlee_bam.module_local_reliability.v1"
MODULE_LOCAL_OUTPUT_CAP_SCHEMA = "kmlee_bam.module_local_output_cap.v1"
MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA = (
    "kmlee_bam.module_local_compartment_graph.v1"
)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_string_sequence(values: Iterable[str]) -> str:
    payload = json.dumps(
        [str(value) for value in values],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _is_sha256(value: object) -> bool:
    text = str(value)
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _scalar(archive: np.lib.npyio.NpzFile, key: str) -> object:
    value = np.asarray(archive[key])
    if value.size != 1:
        raise ValueError(f"module-local artifact field {key!r} must be scalar")
    return value.reshape(-1)[0].item()


@dataclass(frozen=True)
class ModuleLocalReliabilityArtifact:
    reliability: np.ndarray
    reliable_mask: np.ndarray
    context_names: tuple[str, ...]
    module_names: tuple[str, ...]
    train_donor_names: tuple[str, ...]
    split_seed: int
    split_rule: str
    source_context_sha256: str
    registry_sha256: str
    activity_dictionary_sha256: str
    train_donor_allowlist_sha256: str
    artifact_sha256: str
    source_path: str

    @property
    def effective_reliability(self) -> np.ndarray:
        return (
            np.asarray(self.reliability, dtype=np.float32)
            * np.asarray(self.reliable_mask, dtype=np.float32)
        )


@dataclass(frozen=True)
class ModuleLocalOutputCapArtifact:
    output_cap: float
    absolute_quantile: float
    personal_rank: int
    personal_coefficient_count: int
    train_donor_names: tuple[str, ...]
    checkpoint_sha256: str
    source_config_sha256: str
    source_context_sha256: str
    train_donor_allowlist_sha256: str
    coefficient_origin: str
    artifact_sha256: str
    source_path: str


@dataclass(frozen=True)
class ModuleLocalCompartmentGraphArtifact:
    """Immutable sparse module-overlap graph for the nonlinear local path."""

    adjacency: np.ndarray
    module_names: tuple[str, ...]
    registry_sha256: str
    topk: int
    minimum_jaccard: float
    artifact_sha256: str
    source_path: str


def load_module_local_compartment_graph_artifact(
    path: str | Path,
    *,
    expected_module_names: Sequence[str],
    expected_registry_sha256: Optional[str] = None,
) -> ModuleLocalCompartmentGraphArtifact:
    """Load a train/test-label-free, row-normalized module-overlap graph.

    The diagonal must be exactly zero.  A zero row is allowed for a singleton
    or otherwise isolated module; the runtime then leaves that compartment
    independent instead of inventing a connection.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(
            f"module-local compartment graph does not exist: {source}"
        )
    required = {
        "schema_version",
        "adjacency",
        "module_names",
        "registry_sha256",
        "topk",
        "minimum_jaccard",
        "validation_donors_used",
        "test_donors_used",
    }
    with np.load(source, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(
                "module-local compartment graph is missing fields: "
                f"{missing}"
            )
        schema = str(_scalar(archive, "schema_version"))
        if schema != MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA:
            raise ValueError(
                f"unsupported module-local compartment graph schema: {schema!r}"
            )
        if bool(_scalar(archive, "validation_donors_used")):
            raise ValueError(
                "validation donors were used to build the compartment graph"
            )
        if bool(_scalar(archive, "test_donors_used")):
            raise ValueError("test donors were used to build the compartment graph")
        adjacency = np.asarray(archive["adjacency"], dtype=np.float32).copy()
        module_names = tuple(
            str(value) for value in np.asarray(archive["module_names"]).tolist()
        )
        registry_sha256 = str(_scalar(archive, "registry_sha256"))
        topk = int(_scalar(archive, "topk"))
        minimum_jaccard = float(_scalar(archive, "minimum_jaccard"))

    expected_modules = tuple(str(value) for value in expected_module_names)
    if module_names != expected_modules:
        raise ValueError(
            "module-local compartment graph module names/order do not match"
        )
    expected_shape = (len(expected_modules), len(expected_modules))
    if adjacency.shape != expected_shape:
        raise ValueError(
            "module-local compartment adjacency shape mismatch: "
            f"{adjacency.shape} != {expected_shape}"
        )
    if not np.isfinite(adjacency).all():
        raise ValueError("module-local compartment adjacency is non-finite")
    if np.any(adjacency < 0.0):
        raise ValueError("module-local compartment adjacency must be non-negative")
    if not np.array_equal(np.diag(adjacency), np.zeros(len(expected_modules))):
        raise ValueError("module-local compartment adjacency diagonal must be zero")
    row_sum = adjacency.sum(axis=1)
    nonzero = row_sum > 0.0
    if np.any(np.abs(row_sum[nonzero] - 1.0) > 1.0e-5):
        raise ValueError(
            "non-empty module-local compartment adjacency rows must sum to one"
        )
    if topk <= 0:
        raise ValueError("module-local compartment graph topk must be positive")
    if bool(((adjacency > 0.0).sum(axis=1) > topk).any()):
        raise ValueError(
            "module-local compartment graph has more nonzero edges than topk"
        )
    if not np.isfinite(minimum_jaccard) or not 0.0 <= minimum_jaccard <= 1.0:
        raise ValueError(
            "module-local compartment minimum_jaccard must lie in [0,1]"
        )
    if not _is_sha256(registry_sha256):
        raise ValueError("module-local compartment registry SHA-256 is invalid")
    if (
        expected_registry_sha256 is not None
        and registry_sha256 != str(expected_registry_sha256)
    ):
        raise ValueError(
            "module-local compartment registry SHA-256 mismatch: "
            f"{registry_sha256} != {expected_registry_sha256}"
        )

    return ModuleLocalCompartmentGraphArtifact(
        adjacency=adjacency,
        module_names=module_names,
        registry_sha256=registry_sha256,
        topk=topk,
        minimum_jaccard=minimum_jaccard,
        artifact_sha256=sha256_file(source),
        source_path=str(source),
    )


def load_module_local_reliability_artifact(
    path: str | Path,
    *,
    expected_context_names: Sequence[str],
    expected_module_names: Sequence[str],
    expected_train_donor_names: Sequence[str],
    expected_source_context_sha256: Optional[str] = None,
    expected_registry_sha256: Optional[str] = None,
    expected_activity_dictionary_sha256: Optional[str] = None,
) -> ModuleLocalReliabilityArtifact:
    """Load and fully validate a train-only context-by-module artifact."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"module-local reliability artifact does not exist: {source}")
    required = {
        "schema_version",
        "context_module_reliability",
        "context_module_reliable_mask",
        "context_names",
        "module_names",
        "train_donor_allowlist",
        "train_donor_allowlist_sha256",
        "split_seed",
        "split_rule",
        "source_context_sha256",
        "registry_sha256",
        "activity_dictionary_sha256",
        "validation_donors_used",
        "test_donors_used",
    }
    with np.load(source, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(
                f"module-local reliability artifact is missing fields: {missing}"
            )
        schema = str(_scalar(archive, "schema_version"))
        if schema != MODULE_LOCAL_RELIABILITY_SCHEMA:
            raise ValueError(
                f"unsupported module-local reliability schema: {schema!r}"
            )
        if bool(_scalar(archive, "validation_donors_used")):
            raise ValueError("validation donors were used to fit module reliability")
        if bool(_scalar(archive, "test_donors_used")):
            raise ValueError("test donors were used to fit module reliability")

        reliability = np.asarray(
            archive["context_module_reliability"], dtype=np.float32
        ).copy()
        reliable_mask = np.asarray(
            archive["context_module_reliable_mask"], dtype=bool
        ).copy()
        context_names = tuple(
            str(value) for value in np.asarray(archive["context_names"]).tolist()
        )
        module_names = tuple(
            str(value) for value in np.asarray(archive["module_names"]).tolist()
        )
        train_donor_names = tuple(
            str(value)
            for value in np.asarray(archive["train_donor_allowlist"]).tolist()
        )
        split_seed = int(_scalar(archive, "split_seed"))
        split_rule = str(_scalar(archive, "split_rule"))
        source_context_sha256 = str(_scalar(archive, "source_context_sha256"))
        registry_sha256 = str(_scalar(archive, "registry_sha256"))
        activity_dictionary_sha256 = str(
            _scalar(archive, "activity_dictionary_sha256")
        )
        train_hash = str(_scalar(archive, "train_donor_allowlist_sha256"))

    expected_context = tuple(str(value) for value in expected_context_names)
    expected_modules = tuple(str(value) for value in expected_module_names)
    expected_donors = tuple(str(value) for value in expected_train_donor_names)
    if context_names != expected_context:
        raise ValueError("module-local artifact context names/order do not match")
    if module_names != expected_modules:
        raise ValueError("module-local artifact module names/order do not match")
    if train_donor_names != expected_donors:
        raise ValueError(
            "module-local artifact fitting donors differ from the runtime fit allow-list"
        )
    expected_shape = (len(expected_context), len(expected_modules))
    if reliability.shape != expected_shape or reliable_mask.shape != expected_shape:
        raise ValueError(
            "module-local reliability/mask shape mismatch: "
            f"{reliability.shape}/{reliable_mask.shape} != {expected_shape}"
        )
    if not np.isfinite(reliability).all():
        raise ValueError("module-local reliability contains non-finite values")
    if np.any(reliability < 0.0) or np.any(reliability > 1.0):
        raise ValueError("module-local reliability must lie in [0,1]")
    if not np.any(reliable_mask):
        raise ValueError("module-local artifact has no reliable context-module entries")
    if split_seed < 0 or not split_rule.strip():
        raise ValueError("module-local split provenance is incomplete")

    observed_train_hash = sha256_string_sequence(train_donor_names)
    if train_hash != observed_train_hash:
        raise ValueError("module-local train-donor allow-list hash is invalid")
    for label, observed in (
        ("source_context_sha256", source_context_sha256),
        ("registry_sha256", registry_sha256),
        ("activity_dictionary_sha256", activity_dictionary_sha256),
        ("train_donor_allowlist_sha256", train_hash),
    ):
        if not _is_sha256(observed):
            raise ValueError(f"module-local {label} is not a SHA-256 digest")

    expected_hashes = (
        ("source context", source_context_sha256, expected_source_context_sha256),
        ("registry", registry_sha256, expected_registry_sha256),
        (
            "activity dictionary",
            activity_dictionary_sha256,
            expected_activity_dictionary_sha256,
        ),
    )
    for label, observed, expected in expected_hashes:
        if expected is not None and observed != str(expected):
            raise ValueError(
                f"module-local {label} SHA-256 mismatch: {observed} != {expected}"
            )

    return ModuleLocalReliabilityArtifact(
        reliability=reliability,
        reliable_mask=reliable_mask,
        context_names=context_names,
        module_names=module_names,
        train_donor_names=train_donor_names,
        split_seed=split_seed,
        split_rule=split_rule,
        source_context_sha256=source_context_sha256,
        registry_sha256=registry_sha256,
        activity_dictionary_sha256=activity_dictionary_sha256,
        train_donor_allowlist_sha256=train_hash,
        artifact_sha256=sha256_file(source),
        source_path=str(source),
    )


def load_module_local_output_cap_artifact(
    path: str | Path,
    *,
    expected_train_donor_names: Sequence[str],
    expected_checkpoint_sha256: str,
    expected_source_config_sha256: str,
    expected_source_context_sha256: str,
    expected_absolute_quantile: float,
    expected_personal_rank: int = 2,
) -> ModuleLocalOutputCapArtifact:
    """Load a cap calibrated only from a sealed frozen rank-2 comparator."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(
            f"module-local output-cap artifact does not exist: {source}"
        )
    required = {
        "schema_version",
        "output_cap",
        "absolute_quantile",
        "personal_rank",
        "personal_coefficient_count",
        "train_donor_allowlist",
        "train_donor_allowlist_sha256",
        "checkpoint_sha256",
        "source_config_sha256",
        "source_context_sha256",
        "coefficient_origin",
        "validation_donors_used",
        "test_donors_used",
    }
    with np.load(source, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(
                f"module-local output-cap artifact is missing fields: {missing}"
            )
        schema = str(_scalar(archive, "schema_version"))
        if schema != MODULE_LOCAL_OUTPUT_CAP_SCHEMA:
            raise ValueError(
                f"unsupported module-local output-cap schema: {schema!r}"
            )
        if bool(_scalar(archive, "validation_donors_used")):
            raise ValueError("validation donors were used to calibrate output cap")
        if bool(_scalar(archive, "test_donors_used")):
            raise ValueError("test donors were used to calibrate output cap")
        output_cap = float(_scalar(archive, "output_cap"))
        absolute_quantile = float(_scalar(archive, "absolute_quantile"))
        personal_rank = int(_scalar(archive, "personal_rank"))
        coefficient_count = int(
            _scalar(archive, "personal_coefficient_count")
        )
        train_donor_names = tuple(
            str(value)
            for value in np.asarray(archive["train_donor_allowlist"]).tolist()
        )
        train_hash = str(_scalar(archive, "train_donor_allowlist_sha256"))
        checkpoint_sha256 = str(_scalar(archive, "checkpoint_sha256"))
        source_config_sha256 = str(_scalar(archive, "source_config_sha256"))
        source_context_sha256 = str(_scalar(archive, "source_context_sha256"))
        coefficient_origin = str(_scalar(archive, "coefficient_origin"))

    if not np.isfinite(output_cap) or output_cap <= 0.0:
        raise ValueError("module-local output cap must be finite and positive")
    if not 0.5 < absolute_quantile < 1.0:
        raise ValueError("module-local output-cap quantile must lie in (0.5,1)")
    if not np.isclose(
        absolute_quantile,
        float(expected_absolute_quantile),
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError(
            "module-local output-cap quantile differs from the configured value"
        )
    if personal_rank != int(expected_personal_rank):
        raise ValueError(
            "module-local output cap was not calibrated from the expected rank"
        )
    if coefficient_count <= 0:
        raise ValueError("module-local output-cap artifact has no coefficients")
    expected_donors = tuple(
        str(value) for value in expected_train_donor_names
    )
    if train_donor_names != expected_donors:
        raise ValueError(
            "module-local output-cap donors differ from the runtime fit allow-list"
        )
    if train_hash != sha256_string_sequence(train_donor_names):
        raise ValueError("module-local output-cap donor hash is invalid")
    if checkpoint_sha256 != str(expected_checkpoint_sha256):
        raise ValueError(
            "module-local output-cap checkpoint SHA-256 mismatch"
        )
    if source_config_sha256 != str(expected_source_config_sha256):
        raise ValueError(
            "module-local output-cap source-config SHA-256 mismatch"
        )
    if source_context_sha256 != str(expected_source_context_sha256):
        raise ValueError(
            "module-local output-cap context SHA-256 mismatch"
        )
    for label, digest in (
        ("checkpoint_sha256", checkpoint_sha256),
        ("source_config_sha256", source_config_sha256),
        ("source_context_sha256", source_context_sha256),
        ("train_donor_allowlist_sha256", train_hash),
    ):
        if not _is_sha256(digest):
            raise ValueError(f"module-local output-cap {label} is invalid")
    expected_origin = (
        "frozen_rank2_train_donor_personal_coeff_absolute_quantile"
    )
    if coefficient_origin != expected_origin:
        raise ValueError("module-local output-cap coefficient origin is invalid")

    return ModuleLocalOutputCapArtifact(
        output_cap=output_cap,
        absolute_quantile=absolute_quantile,
        personal_rank=personal_rank,
        personal_coefficient_count=coefficient_count,
        train_donor_names=train_donor_names,
        checkpoint_sha256=checkpoint_sha256,
        source_config_sha256=source_config_sha256,
        source_context_sha256=source_context_sha256,
        train_donor_allowlist_sha256=train_hash,
        coefficient_origin=coefficient_origin,
        artifact_sha256=sha256_file(source),
        source_path=str(source),
    )
