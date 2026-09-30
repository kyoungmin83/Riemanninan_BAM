#!/usr/bin/env python
"""Prepare DLPFC+MTG obs and donor table for donor-level pathology split.

This is Phase 1 of the DLPFC+MTG integration pipeline. It does not touch the
expression matrices. It only builds small metadata artifacts used by
make_pathology_aware_donor_split.py:

  - combined_obs_for_split.csv
  - donor_pathology_table.csv
  - donor_value_conflicts.csv
  - phase1_input_manifest.json

All large outputs should live on the gstorage-backed project path, not sv7
local scratch.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr


DEFAULT_DLPFC_ZARR = (
    "/home/kmlee/project/kmlee_bam/data/scRNA/SEA_AD/whole_taxnomy/"
    "sea_ad_dlpfc_rechunked.zarr"
)
DEFAULT_MTG_ZARR = (
    "/home/kmlee/project/kmlee_bam/data/scRNA/SEA_AD/whole_taxnomy/"
    "sea_ad_mtg_rechunked.zarr"
)
DEFAULT_OUT_DIR = (
    "/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/"
    "module_mapping_outputs/dlpfc_mtg/pathology_aware_split_inputs"
)

DONOR_KEY = "donor_id"
CELLTYPE_KEY = "Subclass"
TECH_SOURCE_COLS = ["assay", "assay_ontology_term_id", "PMI"]
CELL_QC_COLS = ["Number of UMIs", "Genes detected", "Fraction mitochrondrial UMIs"]
DONOR_LEVEL_COLS = [
    "disease",
    "Neurotypical reference",
    "ADNC",
    "Braak stage",
    "Thal phase",
    "CERAD score",
    "LATE-NC stage",
    "Lewy body disease pathology",
    "Microinfarct pathology",
    "APOE4 status",
    "Cognitive status",
    "Age at death",
    "sex",
]
CELLTYPE_REGION_KEY = "celltype_region"
KEEP_COLS = [DONOR_KEY, CELLTYPE_KEY] + TECH_SOURCE_COLS + CELL_QC_COLS + DONOR_LEVEL_COLS


def log(msg: str) -> None:
    print(f"[phase1-inputs] {msg}", flush=True)


def read_zarr_obs(path: str) -> pd.DataFrame:
    root = zarr.open(path, mode="r")
    try:
        from anndata.experimental import read_elem  # type: ignore
    except Exception:
        from anndata.io import read_elem  # type: ignore
    obs = read_elem(root["obs"])
    if not isinstance(obs, pd.DataFrame):
        raise TypeError(f"obs in {path} is not a pandas DataFrame.")
    return obs.copy()


def normalise_value(x: Any) -> Any:
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "na", "n/a"}:
        return np.nan
    return x


def collapse_mode_nonnull(values: pd.Series) -> tuple[Any, int, str]:
    vals = values.map(normalise_value).dropna()
    if vals.empty:
        return np.nan, 0, ""
    counts = vals.astype(str).value_counts(dropna=True)
    chosen = counts.index[0]
    return chosen, int(len(counts)), "|".join([f"{k}:{int(v)}" for k, v in counts.items()])


def add_region(obs: pd.DataFrame, region: str, source_zarr: str) -> pd.DataFrame:
    missing = [c for c in KEEP_COLS if c not in obs.columns]
    if missing:
        raise KeyError(f"{region} obs is missing required columns: {missing}")

    out = obs.loc[:, KEEP_COLS].copy()
    out.index = out.index.astype(str)
    out["source_obs_index"] = out.index
    out["brain_region"] = region
    out[CELLTYPE_REGION_KEY] = out[CELLTYPE_KEY].astype(str) + "__" + region
    out["source_zarr"] = source_zarr
    out.index = [f"{region}::{idx}" for idx in out.index]
    return out


def build_donor_table(combined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []

    for donor, sub in combined.groupby(DONOR_KEY, sort=True, observed=True):
        row: dict[str, Any] = {
            "donor_id": str(donor),
            "n_cells": int(len(sub)),
            "n_cells_DLPFC": int((sub["brain_region"] == "DLPFC").sum()),
            "n_cells_MTG": int((sub["brain_region"] == "MTG").sum()),
            "n_regions": int(sub["brain_region"].nunique(dropna=True)),
            "regions": ",".join(sorted(sub["brain_region"].dropna().astype(str).unique().tolist())),
        }
        for col in DONOR_LEVEL_COLS:
            chosen, n_unique, detail = collapse_mode_nonnull(sub[col])
            row[col] = chosen
            if n_unique > 1:
                conflicts.append(
                    {
                        "donor_id": str(donor),
                        "column": col,
                        "n_unique": n_unique,
                        "values": detail,
                    }
                )
        rows.append(row)

    donor_table = pd.DataFrame(rows).set_index("donor_id")
    conflict_table = pd.DataFrame(conflicts)
    if not conflict_table.empty:
        conflict_table = conflict_table.sort_values(["donor_id", "column"]).reset_index(drop=True)
    return donor_table, conflict_table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dlpfc-zarr", default=DEFAULT_DLPFC_ZARR)
    ap.add_argument("--mtg-zarr", default=DEFAULT_MTG_ZARR)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"reading DLPFC obs: {args.dlpfc_zarr}")
    dlpfc = add_region(read_zarr_obs(args.dlpfc_zarr), "DLPFC", args.dlpfc_zarr)
    log(f"DLPFC cells={len(dlpfc):,}")

    log(f"reading MTG obs: {args.mtg_zarr}")
    mtg = add_region(read_zarr_obs(args.mtg_zarr), "MTG", args.mtg_zarr)
    log(f"MTG cells={len(mtg):,}")

    combined = pd.concat([dlpfc, mtg], axis=0, copy=False)
    if not combined.index.is_unique:
        raise ValueError("Combined obs index is not unique after region prefixing.")

    donor_table, conflict_table = build_donor_table(combined)

    combined_path = out_dir / "combined_obs_for_split.csv"
    donor_path = out_dir / "donor_pathology_table.csv"
    conflicts_path = out_dir / "donor_value_conflicts.csv"
    manifest_path = out_dir / "phase1_input_manifest.json"

    log(f"writing combined obs: {combined_path}")
    combined.to_csv(combined_path, index=True)
    donor_table.to_csv(donor_path, index=True)
    conflict_table.to_csv(conflicts_path, index=False)

    region_counts = combined["brain_region"].value_counts().sort_index().to_dict()
    donor_region_counts = (
        combined.groupby([DONOR_KEY, "brain_region"], observed=True)
        .size()
        .unstack(fill_value=0)
        .astype(int)
    )
    mtg_only = donor_region_counts.index[
        (donor_region_counts.get("DLPFC", 0) == 0) & (donor_region_counts.get("MTG", 0) > 0)
    ].astype(str).tolist()

    manifest = {
        "dlpfc_zarr": args.dlpfc_zarr,
        "mtg_zarr": args.mtg_zarr,
        "out_dir": str(out_dir),
        "n_cells": int(len(combined)),
        "n_cells_by_region": {str(k): int(v) for k, v in region_counts.items()},
        "n_donors": int(len(donor_table)),
        "n_donors_with_both_regions": int((donor_region_counts.gt(0).sum(axis=1) == 2).sum()),
        "mtg_only_donors": mtg_only,
        "n_donor_value_conflicts": int(len(conflict_table)),
        "outputs": {
            "combined_obs_for_split": str(combined_path),
            "donor_pathology_table": str(donor_path),
            "donor_value_conflicts": str(conflicts_path),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    log(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
