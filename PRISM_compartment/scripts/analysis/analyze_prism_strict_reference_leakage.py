#!/usr/bin/env python3
"""Probe nuisance leakage on the strict pathology-clean reference subset."""

from __future__ import annotations

import argparse
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
    parser.add_argument("--target-split", default="val", choices=("val",))
    parser.add_argument("--out-json", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = np.load(args.extract, allow_pickle=True)
    with open(args.config, encoding="utf-8") as handle:
        config = json.load(handle)

    split = data["split"].astype(str)
    present_splits = set(split.tolist())
    if "test" in present_splits:
        raise RuntimeError("strict-reference validation probe refuses extracts containing test")

    root = zarr.open(config["data"]["zarr_path"], mode="r")
    obs = read_zarr_dataframe(root, "obs")
    row_index = data["row_index"].astype(np.int64)
    strict_reference = build_pathology_masks(obs)["is_reference_origin"][row_index]
    selected = strict_reference & np.isin(split, ("train", args.target_split))

    X = data["z"][selected].astype(float)
    donor = data["donor"][selected].astype(int)
    selected_split = split[selected]
    targets = {
        "celltype": data["celltype"][selected].astype(int),
        "sex": data["sex"][selected].astype(int),
        "technology": data["tech"][selected].astype(int),
        "region": data["region_id"][selected].astype(int),
    }

    result = {
        "schema_version": "kmlee_bam.prism_strict_reference_leakage.v1",
        "scope": "strict pathology-clean reference cells; fit train donors; evaluate validation donors",
        "test_used": False,
        "extract": os.path.abspath(args.extract),
        "config": os.path.abspath(args.config),
        "target_split": args.target_split,
        "counts": {
            name: int(np.sum(selected_split == name))
            for name in ("train", args.target_split)
        },
        "donor_counts": {
            name: int(np.unique(donor[selected_split == name]).size)
            for name in ("train", args.target_split)
        },
        "celltype_counts": {
            name: int(np.unique(targets["celltype"][selected_split == name]).size)
            for name in ("train", args.target_split)
        },
        "probes": {},
    }
    for index, (name, target) in enumerate(targets.items()):
        result["probes"][name] = classify_train_to_split(
            X,
            target,
            donor,
            selected_split,
            args.target_split,
            seed=20260826 + index,
        )

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(json_safe(result), handle, indent=2, allow_nan=False)
    print(
        "[strict-reference-leakage] "
        + " ".join(
            f"z->{name}={probe.get('balanced_accuracy', float('nan')):.3f}"
            for name, probe in result["probes"].items()
        ),
        flush=True,
    )
    print(f"[strict-reference-leakage] wrote {args.out_json}", flush=True)


if __name__ == "__main__":
    main()
