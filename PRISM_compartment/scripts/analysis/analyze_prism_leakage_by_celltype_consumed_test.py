#!/usr/bin/env python3
"""Train-to-consumed-test leakage probes by cell type under a fixed SV7 scope."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys

import numpy as np


ANALYSIS_SCRIPTS = os.environ.get(
    "PRISM_ANALYSIS_SCRIPTS", "/home/kmlee/project_sv6/kmlee_bam/scripts"
)
sys.path.insert(0, ANALYSIS_SCRIPTS)

from analyze_prism_leakage import (  # noqa: E402
    aggregate_donor_celltype,
    classify_train_to_split,
    donor_target,
    regress_train_to_split,
)


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


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
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", required=True)
    parser.add_argument("--scope-lock", required=True)
    parser.add_argument("--expected-lock-sha256", required=True)
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args()

    scope, lock_sha256 = verify_scope(args.scope_lock, args.expected_lock_sha256)
    data = np.load(args.extract, allow_pickle=True)
    split_cell = data["split"].astype(str)
    present = set(split_cell.tolist())
    if present != {"train", "test"}:
        raise ValueError(
            "consumed-test cell-type audit requires exactly train+test rows; "
            f"present={sorted(present)}"
        )

    z_dc, donor, celltype, split = aggregate_donor_celltype(
        data, data["z"].astype(float)
    )
    personal_dc, donor_p, celltype_p, split_p = aggregate_donor_celltype(
        data, data["personal_code"].astype(float)
    )
    if not (
        np.array_equal(donor, donor_p)
        and np.array_equal(celltype, celltype_p)
        and np.array_equal(split, split_p)
    ):
        raise RuntimeError("z and personal-code aggregation order differs")

    support_count, donor_s, celltype_s, split_s = aggregate_donor_celltype(
        data, data["support_count"].astype(float)[:, None]
    )
    support_reliability, donor_r, celltype_r, split_r = aggregate_donor_celltype(
        data, data["support_reliability"].astype(float)[:, None]
    )
    if not (
        np.array_equal(donor, donor_s)
        and np.array_equal(donor, donor_r)
        and np.array_equal(celltype, celltype_s)
        and np.array_equal(celltype, celltype_r)
        and np.array_equal(split, split_s)
        and np.array_equal(split, split_r)
    ):
        raise RuntimeError("support aggregation order differs")

    celltype_vocab = np.asarray(data["celltype_vocab"], dtype=object)
    pathology_names = [str(name) for name in data["pathology_names"]]
    targets_continuous = {
        "age": donor_target(data, donor, "age"),
        "ADNC_report_only": donor_target(data, donor, "adnc"),
        "support_context_count": support_count[:, 0],
        "support_context_reliability": support_reliability[:, 0],
    }
    for axis, name in enumerate(pathology_names):
        targets_continuous[name] = donor_target(data, donor, "pathology", axis)
    targets_categorical = {
        "sex": np.rint(donor_target(data, donor, "sex")).astype(int),
        "technology_majority": np.rint(donor_target(data, donor, "tech")).astype(int),
    }

    per_celltype = {}
    for ct in sorted(np.unique(celltype).tolist()):
        mask = celltype == ct
        row = {
            "n_train_donors": int(len(np.unique(donor[mask & (split == "train")]))),
            "n_test_donors": int(len(np.unique(donor[mask & (split == "test")]))),
            "z": {},
            "personal_code": {},
        }
        for target, values in targets_categorical.items():
            row["z"][target] = classify_train_to_split(
                z_dc[mask],
                values[mask],
                donor[mask],
                split[mask],
                "test",
                seed=1200 + int(ct),
            )
            row["personal_code"][target] = classify_train_to_split(
                personal_dc[mask],
                values[mask],
                donor[mask],
                split[mask],
                "test",
                seed=2200 + int(ct),
            )
        for target, values in targets_continuous.items():
            row["z"][target] = regress_train_to_split(
                z_dc[mask],
                values[mask],
                donor[mask],
                split[mask],
                "test",
                seed=3200 + int(ct),
            )
            row["personal_code"][target] = regress_train_to_split(
                personal_dc[mask],
                values[mask],
                donor[mask],
                split[mask],
                "test",
                seed=4200 + int(ct),
            )
        per_celltype[str(celltype_vocab[int(ct)])] = row

    output = {
        "schema_version": "kmlee_bam.prism_leakage_by_celltype_consumed_test.v1",
        "scope": "train-donor fit -> already-consumed test-donor descriptive evaluation",
        "test_used": True,
        "official_test_status": "consumed_reuse",
        "independent_holdout": False,
        "validation_roles_locked_before_sv7_test_reuse": True,
        "test_reuse_may_change_validation_roles": False,
        "scope_lock": os.path.abspath(args.scope_lock),
        "scope_lock_sha256": lock_sha256,
        "validation_primary_id": scope["validation_primary_id"],
        "extract": os.path.abspath(args.extract),
        "target_split": "test",
        "per_celltype": per_celltype,
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(json_safe(output), handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    print(
        f"[leakage-by-celltype-consumed-test] wrote {args.out_json}; "
        f"celltypes={len(per_celltype)}; lock={lock_sha256}",
        flush=True,
    )


if __name__ == "__main__":
    main()
