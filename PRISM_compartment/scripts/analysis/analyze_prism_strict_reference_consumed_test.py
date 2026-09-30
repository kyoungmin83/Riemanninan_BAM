#!/usr/bin/env python3
"""Strict-reference nuisance probes on an already-consumed test cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import zarr


ANALYSIS_SCRIPTS = os.environ.get(
    "PRISM_ANALYSIS_SCRIPTS", "/home/kmlee/project_sv6/kmlee_bam/scripts"
)
sys.path.insert(0, ANALYSIS_SCRIPTS)

from analyze_prism_leakage import classify_train_to_split, json_safe  # noqa: E402
from kmlee_bam.data.ordinal_dataset import read_zarr_dataframe  # noqa: E402
from kmlee_bam.data.pathology_masks import build_pathology_masks  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope-lock", required=True)
    parser.add_argument("--expected-lock-sha256", required=True)
    parser.add_argument("--out-json", required=True)
    return parser.parse_args()


def verify_scope(path: str, expected_sha256: str) -> tuple[dict, str]:
    with open(path, "rb") as handle:
        payload = handle.read()
    actual_sha256 = hashlib.sha256(payload).hexdigest()
    if actual_sha256 != expected_sha256:
        raise RuntimeError(
            f"scope lock hash mismatch: expected={expected_sha256} actual={actual_sha256}"
        )
    scope = json.loads(payload.decode("utf-8"))
    if scope.get("official_test_status_at_scope_lock") != "consumed_by_prior_sv6_locked_audit":
        raise RuntimeError("scope does not declare prior test consumption")
    if scope.get("test_reuse_may_change_validation_roles") is not False:
        raise RuntimeError("scope permits post-test checkpoint reselection")
    if scope.get("test_reuse_may_authorize_baseline_replacement") is not False:
        raise RuntimeError("scope permits baseline replacement")
    return scope, actual_sha256


def main() -> None:
    args = parse_args()
    scope, lock_sha256 = verify_scope(args.scope_lock, args.expected_lock_sha256)
    data = np.load(args.extract, allow_pickle=True)
    with open(args.config, encoding="utf-8") as handle:
        config = json.load(handle)

    split = data["split"].astype(str)
    present_splits = set(split.tolist())
    if present_splits != {"train", "test"}:
        raise RuntimeError(
            "consumed-test strict-reference probe requires exactly train+test rows; "
            f"present={sorted(present_splits)}"
        )

    root = zarr.open(config["data"]["zarr_path"], mode="r")
    obs = read_zarr_dataframe(root, "obs")
    row_index = data["row_index"].astype(np.int64)
    strict_reference = build_pathology_masks(obs)["is_reference_origin"][row_index]
    selected = strict_reference & np.isin(split, ("train", "test"))

    features = data["z"][selected].astype(float)
    donor = data["donor"][selected].astype(int)
    selected_split = split[selected]
    targets = {
        "celltype": data["celltype"][selected].astype(int),
        "sex": data["sex"][selected].astype(int),
        "technology": data["tech"][selected].astype(int),
        "region": data["region_id"][selected].astype(int),
    }

    result = {
        "schema_version": "kmlee_bam.prism_strict_reference_consumed_test.v1",
        "scope": "strict pathology-clean reference cells; fit train donors; descriptively evaluate already-consumed test donors",
        "test_used": True,
        "official_test_status": "consumed_reuse",
        "independent_holdout": False,
        "validation_roles_locked_before_sv7_test_reuse": True,
        "test_reuse_may_change_validation_roles": False,
        "scope_lock": os.path.abspath(args.scope_lock),
        "scope_lock_sha256": lock_sha256,
        "validation_primary_id": scope["validation_primary_id"],
        "extract": os.path.abspath(args.extract),
        "config": os.path.abspath(args.config),
        "target_split": "test",
        "counts": {
            name: int(np.sum(selected_split == name)) for name in ("train", "test")
        },
        "donor_counts": {
            name: int(np.unique(donor[selected_split == name]).size)
            for name in ("train", "test")
        },
        "celltype_counts": {
            name: int(np.unique(targets["celltype"][selected_split == name]).size)
            for name in ("train", "test")
        },
        "probes": {},
    }
    for index, (name, target) in enumerate(targets.items()):
        result["probes"][name] = classify_train_to_split(
            features,
            target,
            donor,
            selected_split,
            "test",
            seed=20260828 + index,
        )

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(json_safe(result), handle, indent=2, allow_nan=False)
        handle.write("\n")
    scores = []
    for name, probe in result["probes"].items():
        value = probe.get("balanced_accuracy")
        scores.append(f"z->{name}={value:.3f}" if isinstance(value, (int, float)) else f"z->{name}=NA")
    print("[strict-reference-consumed-test] " + " ".join(scores), flush=True)
    print(f"[strict-reference-consumed-test] wrote {args.out_json}", flush=True)


if __name__ == "__main__":
    main()
