#!/usr/bin/env python3
"""Build the sealed train-only split-half reliability artifact for PRISM.

The input can contain either per-cell module activities or already aggregated
half-pseudobulks.  In per-cell mode the split is deterministic from
``seed|cell_id`` and never consults validation/test outcomes.  Reliability is
the non-negative train-donor Pearson reproducibility for each
context-by-module coordinate; low-support/low-correlation coordinates are
masked rather than selected by downstream validation performance.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import sys
import tempfile

import numpy as np


SCHEMA = "kmlee_bam.module_local_reliability.v1"
SPLIT_RULE = "blake2b64(seed|cell_id)_least_significant_bit"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_string_sequence(values: list[str]) -> str:
    import json

    payload = json.dumps(
        [str(value) for value in values],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def deterministic_half(cell_ids: np.ndarray, seed: int) -> np.ndarray:
    result = np.empty(len(cell_ids), dtype=np.int8)
    for index, value in enumerate(cell_ids):
        payload = f"{int(seed)}|{str(value)}".encode("utf-8")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        result[index] = int.from_bytes(digest, "little") & 1
    return result


def aggregate_per_cell(
    activity: np.ndarray,
    donor: np.ndarray,
    context: np.ndarray,
    cell_id: np.ndarray,
    *,
    n_donors: int,
    n_contexts: int,
    seed: int,
    minimum_cells_per_half: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if activity.ndim != 2:
        raise ValueError("cell_module_activity must be [cells,modules]")
    n_cells, n_modules = activity.shape
    if donor.shape != (n_cells,) or context.shape != (n_cells,) or cell_id.shape != (
        n_cells,
    ):
        raise ValueError("per-cell donor/context/cell_id arrays have wrong shapes")
    if not np.isfinite(activity).all():
        raise ValueError("cell_module_activity contains non-finite values")
    if np.any(donor < 0) or np.any(donor >= n_donors):
        raise ValueError("cell_donor contains an out-of-range index")
    if np.any(context < 0) or np.any(context >= n_contexts):
        raise ValueError("cell_context contains an out-of-range index")

    half = deterministic_half(cell_id, seed)
    sums = [
        np.zeros((n_donors, n_contexts, n_modules), dtype=np.float64)
        for _ in range(2)
    ]
    counts = [
        np.zeros((n_donors, n_contexts), dtype=np.int64) for _ in range(2)
    ]
    for side in (0, 1):
        selected = np.flatnonzero(half == side)
        np.add.at(sums[side], (donor[selected], context[selected]), activity[selected])
        np.add.at(counts[side], (donor[selected], context[selected]), 1)
        sums[side] /= np.maximum(counts[side][..., None], 1)
    observed = [count >= int(minimum_cells_per_half) for count in counts]
    return (
        sums[0].astype(np.float32),
        sums[1].astype(np.float32),
        observed[0],
        observed[1],
    )


def split_half_reliability(
    half_a: np.ndarray,
    half_b: np.ndarray,
    observed_a: np.ndarray,
    observed_b: np.ndarray,
    train_donor: np.ndarray,
    *,
    minimum_donors: int,
    minimum_reliability: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if half_a.shape != half_b.shape or half_a.ndim != 3:
        raise ValueError("half module arrays must share shape [donor,context,module]")
    if observed_a.shape != half_a.shape[:2] or observed_b.shape != half_a.shape[:2]:
        raise ValueError("half observation masks have the wrong shape")
    if train_donor.shape != (half_a.shape[0],):
        raise ValueError("train donor mask has the wrong shape")
    n_contexts, n_modules = half_a.shape[1:]
    reliability = np.zeros((n_contexts, n_modules), dtype=np.float32)
    donor_count = np.zeros((n_contexts, n_modules), dtype=np.int32)

    for context in range(n_contexts):
        base_valid = train_donor & observed_a[:, context] & observed_b[:, context]
        a_context = half_a[:, context]
        b_context = half_b[:, context]
        for module in range(n_modules):
            valid = (
                base_valid
                & np.isfinite(a_context[:, module])
                & np.isfinite(b_context[:, module])
            )
            count = int(valid.sum())
            donor_count[context, module] = count
            if count < int(minimum_donors):
                continue
            left = a_context[valid, module].astype(np.float64)
            right = b_context[valid, module].astype(np.float64)
            left -= left.mean()
            right -= right.mean()
            denominator = float(np.sqrt(np.dot(left, left) * np.dot(right, right)))
            if not np.isfinite(denominator) or denominator <= 1.0e-12:
                continue
            correlation = float(np.dot(left, right) / denominator)
            reliability[context, module] = np.float32(
                np.clip(correlation, 0.0, 1.0)
            )
    reliable_mask = (
        (donor_count >= int(minimum_donors))
        & (reliability >= float(minimum_reliability))
    )
    return reliability, reliable_mask, donor_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Per-cell or half-pseudobulk NPZ")
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-context", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--activity-dictionary", default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--minimum-cells-per-half", type=int, default=5)
    parser.add_argument("--minimum-donors", type=int, default=12)
    parser.add_argument("--minimum-reliability", type=float, default=0.20)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.seed < 0:
        raise ValueError("seed must be non-negative")
    if args.minimum_cells_per_half <= 0 or args.minimum_donors < 3:
        raise ValueError("minimum cell/donor support is invalid")
    if not 0.0 <= args.minimum_reliability <= 1.0:
        raise ValueError("minimum_reliability must lie in [0,1]")

    input_path = Path(args.input).expanduser().resolve()
    with np.load(input_path, allow_pickle=False) as archive:
        required_metadata = {
            "donor_names",
            "donor_split",
            "context_names",
            "module_names",
        }
        missing = sorted(required_metadata.difference(archive.files))
        if missing:
            raise ValueError(f"input archive is missing metadata: {missing}")
        donor_names = np.asarray(archive["donor_names"]).astype(str)
        donor_split = np.asarray(archive["donor_split"]).astype(str)
        context_names = np.asarray(archive["context_names"]).astype(str)
        module_names = np.asarray(archive["module_names"]).astype(str)
        if donor_split.shape != donor_names.shape:
            raise ValueError("donor_split and donor_names disagree")
        train_donor = donor_split == "train"
        if not np.any(train_donor):
            raise ValueError("input archive has no training donors")

        per_cell_keys = {
            "cell_module_activity",
            "cell_donor",
            "cell_context",
            "cell_id",
        }
        half_keys = {
            "module_half_a",
            "module_half_b",
            "observed_half_a",
            "observed_half_b",
        }
        if per_cell_keys.issubset(archive.files):
            half_a, half_b, observed_a, observed_b = aggregate_per_cell(
                np.asarray(archive["cell_module_activity"], dtype=np.float32),
                np.asarray(archive["cell_donor"], dtype=np.int64),
                np.asarray(archive["cell_context"], dtype=np.int64),
                np.asarray(archive["cell_id"]).astype(str),
                n_donors=len(donor_names),
                n_contexts=len(context_names),
                seed=int(args.seed),
                minimum_cells_per_half=int(args.minimum_cells_per_half),
            )
        elif half_keys.issubset(archive.files):
            if "split_seed" not in archive.files or int(
                np.asarray(archive["split_seed"]).item()
            ) != int(args.seed):
                raise ValueError("pre-aggregated halves lack the requested split_seed")
            half_a = np.asarray(archive["module_half_a"], dtype=np.float32)
            half_b = np.asarray(archive["module_half_b"], dtype=np.float32)
            observed_a = np.asarray(archive["observed_half_a"], dtype=bool)
            observed_b = np.asarray(archive["observed_half_b"], dtype=bool)
        else:
            raise ValueError(
                "input must provide either per-cell module values or both half pseudobulks"
            )

    expected_shape = (len(donor_names), len(context_names), len(module_names))
    if half_a.shape != expected_shape or half_b.shape != expected_shape:
        raise ValueError(
            f"half-pseudobulk shape mismatch: {half_a.shape}/{half_b.shape} != {expected_shape}"
        )
    reliability, reliable_mask, donor_count = split_half_reliability(
        half_a,
        half_b,
        observed_a,
        observed_b,
        train_donor,
        minimum_donors=int(args.minimum_donors),
        minimum_reliability=float(args.minimum_reliability),
    )
    if not np.any(reliable_mask):
        raise RuntimeError("no context-module coordinate passed the fixed reliability rule")

    train_names = [str(value) for value in donor_names[train_donor].tolist()]
    registry_path = Path(args.registry).expanduser().resolve()
    activity_path = (
        Path(args.activity_dictionary).expanduser().resolve()
        if args.activity_dictionary
        else registry_path
    )
    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": np.asarray(SCHEMA),
        "context_module_reliability": reliability,
        "context_module_reliable_mask": reliable_mask,
        "context_module_train_donor_count": donor_count,
        "context_names": context_names.astype(str),
        "module_names": module_names.astype(str),
        "train_donor_allowlist": np.asarray(train_names, dtype=str),
        "train_donor_allowlist_sha256": np.asarray(
            sha256_string_sequence(train_names)
        ),
        "split_seed": np.asarray(int(args.seed), dtype=np.int64),
        "split_rule": np.asarray(SPLIT_RULE),
        "source_context_sha256": np.asarray(sha256_file(args.source_context)),
        "registry_sha256": np.asarray(sha256_file(registry_path)),
        "activity_dictionary_sha256": np.asarray(sha256_file(activity_path)),
        "validation_donors_used": np.asarray(False),
        "test_donors_used": np.asarray(False),
        "minimum_cells_per_half": np.asarray(
            int(args.minimum_cells_per_half), dtype=np.int64
        ),
        "minimum_donors": np.asarray(int(args.minimum_donors), dtype=np.int64),
        "minimum_reliability": np.asarray(
            float(args.minimum_reliability), dtype=np.float32
        ),
        "input_source_sha256": np.asarray(sha256_file(input_path)),
    }
    with tempfile.NamedTemporaryFile(
        prefix=output_path.name + ".",
        suffix=".tmp.npz",
        dir=output_path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
    try:
        np.savez_compressed(temporary_path, **payload)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()
    print(
        f"wrote {output_path} sha256={sha256_file(output_path)} "
        f"train_donors={len(train_names)} reliable="
        f"{int(reliable_mask.sum())}/{int(reliable_mask.size)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FATAL: {error}", file=sys.stderr, flush=True)
        raise
