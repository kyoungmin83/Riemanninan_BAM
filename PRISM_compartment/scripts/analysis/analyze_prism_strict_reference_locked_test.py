#!/usr/bin/env python3
"""Probe strict-reference nuisance leakage on an explicitly authorized locked test set."""

from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np
import zarr

from analyze_prism_leakage import classify_train_to_split, json_safe
from kmlee_bam.data.ordinal_dataset import read_zarr_dataframe
from kmlee_bam.data.pathology_masks import build_pathology_masks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--pretest-lock", required=True)
    parser.add_argument("--expected-lock-sha256", required=True)
    parser.add_argument("--out-json", required=True)
    return parser.parse_args()


def verify_lock(path: str, expected_sha256: str) -> tuple[dict, str]:
    with open(path, "rb") as handle:
        payload = handle.read()
    actual_sha256 = hashlib.sha256(payload).hexdigest()
    if actual_sha256 != expected_sha256:
        raise RuntimeError(
            f"pre-test lock hash mismatch: expected={expected_sha256} actual={actual_sha256}"
        )
    lock = json.loads(payload.decode("utf-8"))
    if lock.get("official_test_previously_used") is not False:
        raise RuntimeError("pre-test lock does not certify a previously sealed test set")
    if lock.get("test_audit_may_change_selection") is not False:
        raise RuntimeError("test audit must not be allowed to change the locked selection")
    return lock, actual_sha256


def main() -> None:
    args = parse_args()
    lock, lock_sha256 = verify_lock(args.pretest_lock, args.expected_lock_sha256)
    data = np.load(args.extract, allow_pickle=True)
    with open(args.config, encoding="utf-8") as handle:
        config = json.load(handle)

    split = data["split"].astype(str)
    present_splits = set(split.tolist())
    if present_splits != {"train", "test"}:
        raise RuntimeError(
            "locked-test strict-reference probe requires exactly train+test rows; "
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
        "schema_version": "kmlee_bam.prism_strict_reference_locked_test.v1",
        "scope": "strict pathology-clean reference cells; fit train donors; evaluate locked test donors",
        "test_used": True,
        "selection_locked_before_test": True,
        "pretest_lock": os.path.abspath(args.pretest_lock),
        "pretest_lock_sha256": lock_sha256,
        "locked_primary_id": next(
            record["id"]
            for record in lock["locked_roles"]
            if record["role"] == "primary_validation_candidate"
        ),
        "extract": os.path.abspath(args.extract),
        "config": os.path.abspath(args.config),
        "target_split": "test",
        "counts": {
            name: int(np.sum(selected_split == name))
            for name in ("train", "test")
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
    print(
        "[strict-reference-locked-test] "
        + " ".join(
            f"z->{name}={probe.get('balanced_accuracy', float('nan')):.3f}"
            for name, probe in result["probes"].items()
        ),
        flush=True,
    )
    print(f"[strict-reference-locked-test] wrote {args.out_json}", flush=True)


if __name__ == "__main__":
    main()
