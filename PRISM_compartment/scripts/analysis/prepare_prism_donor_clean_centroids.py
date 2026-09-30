#!/usr/bin/env python3
"""Aggregate frozen PRISM cell extracts into donor-by-cell-type clean centroids.

The output is intentionally small enough to move between hosts.  It contains
the full centroid and two deterministic cell halves for a split-half noise
ceiling.  No model is fitted and the locked-test labels are never used to tune
parameters.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def scalar_text(data: np.lib.npyio.NpzFile, key: str) -> str:
    value = data[key]
    return str(value.item() if value.shape == () else value.reshape(-1)[0])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train-val", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    sources = [args.train_val, args.test]
    opened = [np.load(path, allow_pickle=True) for path in sources]
    checkpoint_sha = {scalar_text(data, "source_checkpoint_sha256") for data in opened}
    if len(checkpoint_sha) != 1:
        raise RuntimeError(f"checkpoint mismatch across extracts: {checkpoint_sha}")
    celltype_vocab = np.asarray(opened[0]["celltype_vocab"], dtype=object)
    donor_vocab = np.asarray(opened[0]["donor_vocab"], dtype=object)
    region_vocab = np.asarray(opened[0]["region_vocab"], dtype=object)
    for data in opened[1:]:
        if not np.array_equal(celltype_vocab, data["celltype_vocab"]):
            raise RuntimeError("cell-type vocabulary mismatch")
        if not np.array_equal(donor_vocab, data["donor_vocab"]):
            raise RuntimeError("donor vocabulary mismatch")

    d_count = len(donor_vocab)
    c_count = len(celltype_vocab)
    r_count = len(region_vocab)
    z_dim = int(opened[0]["z"].shape[1])
    p_count = int(opened[0]["pathology"].shape[1])
    centroid = np.full((d_count, c_count, z_dim), np.nan, dtype=np.float32)
    half_centroid = np.full((2, d_count, c_count, z_dim), np.nan, dtype=np.float32)
    count = np.zeros((d_count, c_count), dtype=np.int32)
    half_count = np.zeros((2, d_count, c_count), dtype=np.int32)
    region_fraction = np.zeros((d_count, c_count, r_count), dtype=np.float32)
    split = np.full(d_count, "", dtype=object)
    pathology = np.full((d_count, p_count), np.nan, dtype=np.float32)
    pathology_valid = np.zeros((d_count, p_count), dtype=bool)
    reference_fraction = np.full(d_count, np.nan, dtype=np.float32)

    donor_z: dict[tuple[int, int], list[np.ndarray]] = {}
    donor_half_z: dict[tuple[int, int, int], list[np.ndarray]] = {}
    donor_region: dict[tuple[int, int], list[np.ndarray]] = {}
    donor_path: dict[int, list[np.ndarray]] = {}
    donor_path_valid: dict[int, list[np.ndarray]] = {}
    donor_reference: dict[int, list[np.ndarray]] = {}

    for data in opened:
        row = np.asarray(data["row_index"], dtype=np.int64)
        z = np.asarray(data["z"], dtype=np.float32)
        donor = np.asarray(data["donor"], dtype=np.int64)
        celltype = np.asarray(data["celltype"], dtype=np.int64)
        region = np.asarray(data["region_id"], dtype=np.int64)
        split_cell = np.asarray(data["split"], dtype=str)
        path = np.asarray(data["pathology"], dtype=np.float32)
        path_valid = np.asarray(data["pathology_valid"], dtype=bool)
        is_reference = (
            np.asarray(data["isref"], dtype=np.float32)
            if "isref" in data.files
            else np.zeros(len(row), dtype=np.float32)
        )
        keep = np.isin(split_cell, np.asarray(("val", "test")))
        for donor_id in np.unique(donor[keep]):
            donor_mask = keep & (donor == donor_id)
            donor_split = np.unique(split_cell[donor_mask])
            if len(donor_split) != 1:
                raise RuntimeError(f"donor {donor_id} spans splits {donor_split.tolist()}")
            if split[donor_id] and split[donor_id] != donor_split[0]:
                raise RuntimeError(f"donor {donor_id} split changed")
            split[donor_id] = donor_split[0]
            donor_path.setdefault(int(donor_id), []).append(path[donor_mask])
            donor_path_valid.setdefault(int(donor_id), []).append(path_valid[donor_mask])
            donor_reference.setdefault(int(donor_id), []).append(is_reference[donor_mask])
            for celltype_id in np.unique(celltype[donor_mask]):
                mask = donor_mask & (celltype == celltype_id)
                key = (int(donor_id), int(celltype_id))
                donor_z.setdefault(key, []).append(z[mask])
                donor_region.setdefault(key, []).append(region[mask])
                # Stable pseudo-random split from the immutable global row id.
                half = ((row[mask] * np.int64(1103515245) + np.int64(12345)) & 1).astype(np.int8)
                values = z[mask]
                for half_id in (0, 1):
                    donor_half_z.setdefault((half_id, *key), []).append(values[half == half_id])

    for (donor_id, celltype_id), blocks in donor_z.items():
        values = np.concatenate(blocks, axis=0)
        centroid[donor_id, celltype_id] = values.mean(axis=0)
        count[donor_id, celltype_id] = len(values)
        regions = np.concatenate(donor_region[(donor_id, celltype_id)])
        region_fraction[donor_id, celltype_id] = np.bincount(
            regions, minlength=r_count
        ) / max(len(regions), 1)
        for half_id in (0, 1):
            half_blocks = donor_half_z[(half_id, donor_id, celltype_id)]
            nonempty = [block for block in half_blocks if len(block)]
            if nonempty:
                half_values = np.concatenate(nonempty, axis=0)
                half_centroid[half_id, donor_id, celltype_id] = half_values.mean(axis=0)
                half_count[half_id, donor_id, celltype_id] = len(half_values)

    for donor_id in np.flatnonzero(split != ""):
        values = np.concatenate(donor_path[int(donor_id)], axis=0)
        valid = np.concatenate(donor_path_valid[int(donor_id)], axis=0)
        for axis in range(p_count):
            good = valid[:, axis] & np.isfinite(values[:, axis])
            if np.any(good):
                pathology[donor_id, axis] = float(np.median(values[good, axis]))
                pathology_valid[donor_id, axis] = True
        reference_fraction[donor_id] = float(
            np.mean(np.concatenate(donor_reference[int(donor_id)], axis=0))
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        schema_version=np.asarray("prism.donor_clean_centroids.v1", dtype=object),
        z_centroid=centroid,
        z_half_centroid=half_centroid,
        count=count,
        half_count=half_count,
        region_fraction=region_fraction,
        split=split,
        pathology=pathology,
        pathology_valid=pathology_valid,
        reference_fraction=reference_fraction,
        donor_vocab=donor_vocab,
        celltype_vocab=celltype_vocab,
        region_vocab=region_vocab,
        pathology_names=np.asarray(opened[0]["pathology_names"], dtype=object),
        checkpoint_sha256=np.asarray(next(iter(checkpoint_sha)), dtype=object),
        source_paths=np.asarray([str(path.resolve()) for path in sources], dtype=object),
        source_sha256=np.asarray([sha256(path) for path in sources], dtype=object),
    )
    selected = np.flatnonzero(split != "")
    manifest = {
        "schema_version": "prism.donor_clean_centroids.manifest.v1",
        "output": str(args.out.resolve()),
        "donors": int(len(selected)),
        "validation_donors": int(np.sum(split[selected] == "val")),
        "test_donors": int(np.sum(split[selected] == "test")),
        "celltypes": int(c_count),
        "latent_dim": int(z_dim),
        "minimum_nonzero_group_count": int(count[count > 0].min()),
        "checkpoint_sha256": next(iter(checkpoint_sha)),
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
