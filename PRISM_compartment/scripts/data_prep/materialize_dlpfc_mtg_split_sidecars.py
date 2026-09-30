#!/usr/bin/env python
"""Attach region/source-index columns to DLPFC+MTG pathology split sidecars."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_INPUT_DIR = (
    "/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/"
    "module_mapping_outputs/dlpfc_mtg/pathology_aware_split_inputs"
)
DEFAULT_SPLIT_DIR = (
    "/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/"
    "module_mapping_outputs/dlpfc_mtg/pathology_aware_split"
)


def log(msg: str) -> None:
    print(f"[phase1-sidecars] {msg}", flush=True)


def read_table(path_stem: Path) -> pd.DataFrame:
    parquet_path = path_stem.with_suffix(".parquet")
    csv_path = path_stem.with_suffix(".csv")
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if csv_path.exists():
        return pd.read_csv(csv_path, index_col=0)
    raise FileNotFoundError(f"Missing {parquet_path} or {csv_path}")


def write_table(df: pd.DataFrame, path_stem: Path) -> Path:
    parquet_path = path_stem.with_suffix(".parquet")
    csv_path = path_stem.with_suffix(".csv")
    try:
        df.to_parquet(parquet_path, index=True)
        return parquet_path
    except Exception:
        df.to_csv(csv_path, index=True)
        return csv_path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=DEFAULT_INPUT_DIR)
    ap.add_argument("--split-dir", default=DEFAULT_SPLIT_DIR)
    args = ap.parse_args()

    input_dir = Path(args.input_dir)
    split_dir = Path(args.split_dir)

    combined_obs_stem = input_dir / "combined_obs_for_split"
    split_sidecar_stem = split_dir / "obs_with_pathology_split"

    log(f"reading raw combined obs: {combined_obs_stem}.[parquet|csv]")
    raw = read_table(combined_obs_stem)
    raw = raw.loc[
        :,
        [
            "brain_region",
            "celltype_region",
            "source_obs_index",
            "source_zarr",
            "Number of UMIs",
            "Genes detected",
            "Fraction mitochrondrial UMIs",
        ],
    ]
    raw.index = raw.index.astype(str)

    log(f"reading split sidecar: {split_sidecar_stem}.[parquet|csv]")
    side = read_table(split_sidecar_stem)
    side.index = side.index.astype(str)

    missing = int((~side.index.isin(raw.index)).sum())
    if missing:
        raise ValueError(f"{missing:,} split sidecar rows are absent from combined raw obs.")

    side = side.join(raw, how="left")
    missing_split = int((side["split"] == "__MISSING_DONOR_SPLIT__").sum()) if "split" in side else -1
    if missing_split:
        raise ValueError(f"{missing_split:,} cells have missing donor split.")

    combined_out = split_dir / "obs_with_pathology_split_region"
    dlpfc_out = split_dir / "obs_with_pathology_split_dlpfc"
    mtg_out = split_dir / "obs_with_pathology_split_mtg"
    manifest_out = split_dir / "region_sidecar_manifest.json"

    log(f"writing combined enriched sidecar: {combined_out}")
    combined_written = write_table(side, combined_out)

    region_outputs = {}
    for region, path in [("DLPFC", dlpfc_out), ("MTG", mtg_out)]:
        sub = side.loc[side["brain_region"] == region].copy()
        sub.index = sub["source_obs_index"].astype(str)
        written = write_table(sub, path)
        region_outputs[region] = {"path": str(written), "n_cells": int(len(sub))}
        log(f"wrote {region} sidecar: {written} ({len(sub):,} cells)")

    donor_region_split = (
        side.groupby(["donor_id", "brain_region"], observed=True)["split"]
        .nunique()
        .reset_index(name="n_unique_splits")
    )
    inconsistent = donor_region_split.loc[donor_region_split["n_unique_splits"] != 1]
    if not inconsistent.empty:
        raise ValueError("Found donor/region with multiple split labels.")

    donor_split_consistency = side.groupby("donor_id", observed=True)["split"].nunique()
    donors_multi_split = donor_split_consistency.loc[donor_split_consistency != 1]
    if not donors_multi_split.empty:
        raise ValueError("Found donors assigned to multiple splits across regions.")

    split_counts = side["split"].value_counts().sort_index().to_dict()
    region_split_counts = (
        side.groupby(["brain_region", "split"], observed=True)
        .size()
        .unstack(fill_value=0)
        .astype(int)
    )

    manifest = {
        "input_dir": str(input_dir),
        "split_dir": str(split_dir),
        "n_cells": int(len(side)),
        "n_donors": int(side["donor_id"].nunique()),
        "missing_split_cells": int(missing_split),
        "split_counts_cells": {str(k): int(v) for k, v in split_counts.items()},
        "region_split_counts_cells": {
            str(region): {str(k): int(v) for k, v in row.items()}
            for region, row in region_split_counts.iterrows()
        },
        "outputs": {
            "combined_region_sidecar": str(combined_written),
            **region_outputs,
        },
    }
    with open(manifest_out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    log(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
