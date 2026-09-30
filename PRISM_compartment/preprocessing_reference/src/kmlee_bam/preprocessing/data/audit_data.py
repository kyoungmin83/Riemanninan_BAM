from __future__ import annotations

"""
Donor-level pathology audit for SEA-AD / neurodegeneration scRNA-seq.

Purpose
-------
Before fitting final ordinal bins or training the module-token model, check whether
current donor splits and pathology labels can support biological/pathological-axis
stratification.

This script builds a donor-level table from obs and audits:

    - donor-level pathology-label consistency
    - donor/cell counts per split
    - distributions of Braak, Thal, CERAD, ADNC, LATE-NC, Lewy, Microinfarct, APOE4
    - cross-tabs such as Braak x LATE-NC and ADNC x LATE-NC
    - train/val/test coverage for pathology levels
    - donor-level multi-pathology combinations

It only reads obs metadata. It does NOT read the expression matrix.

Typical use
-----------
python donor_pathology_audit.py --config donor_pathology_audit_cfg.json
"""

import argparse
import json
import re
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd


# =====================================================================
# Config
# =====================================================================
@dataclass
class DonorPathologyAuditConfig:
    # Exactly one input source is required.
    obs_h5ad_path: Optional[str] = None
    zarr_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_parquet_path: Optional[str] = None

    out_dir: str = "donor_pathology_audit"

    donor_key: str = "donor_id"
    donor_key_fallbacks: list[str] = field(
        default_factory=lambda: ["subject_id", "donor_id", "Specimen ID"]
    )
    split_key: str = "split"
    celltype_key: str = "Subclass"
    train_split_name: str = "train"

    tech_columns: list[str] = field(
        default_factory=lambda: [
            "tech_id_model",
            "assay",
            "assay_ontology_term_id",
            "PMI",
            "Number of UMIs",
            "Genes detected",
            "Fraction mitochrondrial UMIs",
            "Specimen ID",
        ]
    )

    pathology_columns: list[str] = field(
        default_factory=lambda: [
            "disease",
            "Neurotypical reference",
            "Cognitive status",
            "ADNC",
            "Braak stage",
            "Thal phase",
            "CERAD score",
            "Lewy body disease pathology",
            "LATE-NC stage",
            "Microinfarct pathology",
            "APOE4 status",
            "Age at death",
            "Years of education",
            "sex",
        ]
    )

    combination_columns: list[str] = field(
        default_factory=lambda: [
            "ADNC",
            "Braak stage",
            "Thal phase",
            "CERAD score",
            "LATE-NC stage",
            "Lewy body disease pathology",
            "Microinfarct pathology",
        ]
    )

    crosstab_pairs: list[list[str]] = field(
        default_factory=lambda: [
            ["Braak stage", "LATE-NC stage"],
            ["ADNC", "LATE-NC stage"],
            ["Thal phase", "LATE-NC stage"],
            ["CERAD score", "LATE-NC stage"],
            ["Lewy body disease pathology", "ADNC"],
            ["Lewy body disease pathology", "LATE-NC stage"],
            ["Microinfarct pathology", "Cognitive status"],
            ["APOE4 status", "ADNC"],
            ["APOE4 status", "LATE-NC stage"],
        ]
    )

    min_train_donors_per_level: int = 3
    save_parquet: bool = True
    save_csv: bool = True
    drop_missing_donor: bool = True


# =====================================================================
# IO
# =====================================================================
def _log(msg: str) -> None:
    print(f"[donor-pathology-audit] {msg}", flush=True)


def load_config(path: str | Path) -> DonorPathologyAuditConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return DonorPathologyAuditConfig(**raw)


def _read_zarr_obs(zarr_path: str | Path) -> pd.DataFrame:
    import zarr

    root = zarr.open(str(zarr_path), mode="r")
    try:
        from anndata.experimental import read_elem  # type: ignore
        obs = read_elem(root["obs"])
        if not isinstance(obs, pd.DataFrame):
            raise TypeError("Zarr obs element is not a pandas DataFrame.")
        return obs.copy()
    except Exception:
        import anndata as ad  # type: ignore
        return ad.read_zarr(str(zarr_path)).obs.copy()


def load_obs(cfg: DonorPathologyAuditConfig) -> pd.DataFrame:
    n_sources = sum(
        x is not None
        for x in [cfg.obs_h5ad_path, cfg.zarr_path, cfg.obs_csv_path, cfg.obs_parquet_path]
    )
    if n_sources != 1:
        raise ValueError(
            "Provide exactly one of obs_h5ad_path, zarr_path, obs_csv_path, obs_parquet_path."
        )

    if cfg.obs_h5ad_path is not None:
        import anndata as ad
        _log(f"loading obs from h5ad: {cfg.obs_h5ad_path}")
        return ad.read_h5ad(cfg.obs_h5ad_path, backed="r").obs.copy()

    if cfg.zarr_path is not None:
        _log(f"loading obs from zarr: {cfg.zarr_path}")
        return _read_zarr_obs(cfg.zarr_path)

    if cfg.obs_csv_path is not None:
        _log(f"loading obs from csv: {cfg.obs_csv_path}")
        return pd.read_csv(cfg.obs_csv_path, index_col=0)

    if cfg.obs_parquet_path is not None:
        _log(f"loading obs from parquet: {cfg.obs_parquet_path}")
        return pd.read_parquet(cfg.obs_parquet_path)

    raise RuntimeError("unreachable")


# =====================================================================
# Helpers
# =====================================================================
def safe_filename(x: Any) -> str:
    x = str(x)
    x = re.sub(r"[^\w.\-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def choose_existing_key(obs: pd.DataFrame, preferred: str, fallbacks: list[str], *, role: str) -> str:
    if preferred in obs.columns:
        return preferred
    for key in fallbacks:
        if key in obs.columns:
            warnings.warn(
                f"{role}: preferred key {preferred!r} not found; using fallback {key!r}.",
                RuntimeWarning,
            )
            return key
    raise KeyError(
        f"No valid {role} key found. Tried preferred={preferred!r}, fallbacks={fallbacks!r}."
    )


def normalise_value(x: Any) -> str:
    if pd.isna(x):
        return "__NA__"
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "na", "n/a"}:
        return "__NA__"
    return s


def collapse_unique_values(s: pd.Series) -> tuple[str, int, str]:
    vals = [normalise_value(x) for x in s.tolist()]
    vals_non_na = sorted({v for v in vals if v != "__NA__"})

    if len(vals_non_na) == 0:
        return "__NA__", 0, ""
    if len(vals_non_na) == 1:
        return vals_non_na[0], 1, vals_non_na[0]
    return "__MULTI__", len(vals_non_na), " | ".join(vals_non_na)


def majority_value(s: pd.Series) -> str:
    vals = pd.Series([normalise_value(x) for x in s.tolist()])
    vals = vals[vals != "__NA__"]
    if len(vals) == 0:
        return "__NA__"
    return str(vals.value_counts().index[0])


def numeric_summary(s: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(s, errors="coerce")
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"mean": np.nan, "median": np.nan, "min": np.nan, "max": np.nan}
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def save_df(df: pd.DataFrame, path_stem: Path, cfg: DonorPathologyAuditConfig, *, index: bool = True) -> None:
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    if cfg.save_csv:
        df.to_csv(path_stem.with_suffix(".csv"), index=index)
    if cfg.save_parquet:
        try:
            df.to_parquet(path_stem.with_suffix(".parquet"), index=index)
        except Exception as e:
            warnings.warn(f"Could not save parquet for {path_stem.name}: {repr(e)}")


def crosstab_with_totals(df: pd.DataFrame, row: str, col: str) -> pd.DataFrame:
    tab = pd.crosstab(df[row].astype(str), df[col].astype(str), dropna=False)
    tab.loc["__TOTAL__"] = tab.sum(axis=0)
    tab["__TOTAL__"] = tab.sum(axis=1)
    return tab


# =====================================================================
# Main audit pieces
# =====================================================================
def build_donor_table(
    obs: pd.DataFrame,
    *,
    donor_key: str,
    split_key: str,
    celltype_key: str,
    pathology_columns: list[str],
    tech_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    available_path = [c for c in pathology_columns if c in obs.columns]
    available_tech = [c for c in tech_columns if c in obs.columns]
    available = list(dict.fromkeys(available_path + available_tech))

    rows = []
    conflicts = []

    for donor, g in obs.groupby(donor_key, observed=True, sort=True):
        row: dict[str, Any] = {
            "donor_id_audit": str(donor),
            "n_cells": int(len(g)),
        }

        if split_key in g.columns:
            val, nuniq, detail = collapse_unique_values(g[split_key])
            row[split_key] = val
            row[f"{split_key}__n_unique"] = nuniq
            if nuniq > 1:
                conflicts.append({"donor_id_audit": str(donor), "column": split_key, "n_unique": nuniq, "values": detail})
        else:
            row[split_key] = "__MISSING_COLUMN__"
            row[f"{split_key}__n_unique"] = 0

        if celltype_key in g.columns:
            row["n_celltypes"] = int(g[celltype_key].astype(str).nunique(dropna=False))
            row["celltypes_seen"] = " | ".join(sorted(g[celltype_key].astype(str).unique()))
            row["majority_celltype"] = majority_value(g[celltype_key])
        else:
            row["n_celltypes"] = 0
            row["celltypes_seen"] = ""
            row["majority_celltype"] = "__MISSING_COLUMN__"

        for col in available:
            val, nuniq, detail = collapse_unique_values(g[col])
            row[col] = val
            row[f"{col}__n_unique"] = nuniq
            if nuniq > 1:
                conflicts.append({"donor_id_audit": str(donor), "column": col, "n_unique": nuniq, "values": detail})

        for col in ["Number of UMIs", "Genes detected", "Fraction mitochrondrial UMIs", "PMI"]:
            if col in g.columns:
                summ = numeric_summary(g[col])
                for k, v in summ.items():
                    row[f"{col}__{k}"] = v

        rows.append(row)

    donor_table = pd.DataFrame(rows).set_index("donor_id_audit")
    conflict_table = pd.DataFrame(conflicts)
    if len(conflict_table) > 0:
        conflict_table = conflict_table.set_index(["donor_id_audit", "column"]).sort_index()
    return donor_table, conflict_table


def build_split_summary(donor_table: pd.DataFrame, split_key: str) -> pd.DataFrame:
    if split_key not in donor_table.columns:
        return pd.DataFrame()
    rows = []
    for split, g in donor_table.groupby(split_key, observed=True, sort=True):
        rows.append({"split": split, "n_donors": int(len(g)), "n_cells": int(g["n_cells"].sum())})
    return pd.DataFrame(rows).set_index("split")


def build_pathology_level_counts(
    donor_table: pd.DataFrame,
    *,
    pathology_columns: list[str],
    split_key: str,
) -> pd.DataFrame:
    rows = []
    cols = [c for c in pathology_columns if c in donor_table.columns]
    for col in cols:
        vc = donor_table[col].astype(str).value_counts(dropna=False)
        for level, n in vc.items():
            rows.append({"column": col, "level": str(level), "split": "__ALL__", "n_donors": int(n)})
        if split_key in donor_table.columns:
            for split, g in donor_table.groupby(split_key, observed=True, sort=True):
                vc_s = g[col].astype(str).value_counts(dropna=False)
                for level, n in vc_s.items():
                    rows.append({"column": col, "level": str(level), "split": str(split), "n_donors": int(n)})
    return pd.DataFrame(rows)


def build_train_level_warnings(level_counts: pd.DataFrame, *, train_split_name: str, min_train_donors_per_level: int) -> pd.DataFrame:
    if level_counts.empty:
        return pd.DataFrame()
    train = level_counts[level_counts["split"].astype(str) == str(train_split_name)].copy()
    train = train[~train["level"].isin(["__NA__", "__MULTI__", "__MISSING_COLUMN__"])]
    warn = train[train["n_donors"] < int(min_train_donors_per_level)].copy()
    warn["warning"] = f"fewer than {min_train_donors_per_level} train donors for this pathology level"
    return warn.sort_values(["column", "n_donors", "level"])


def build_combination_counts(donor_table: pd.DataFrame, *, combination_columns: list[str], split_key: str) -> pd.DataFrame:
    cols = [c for c in combination_columns if c in donor_table.columns]
    if len(cols) == 0:
        return pd.DataFrame()
    use_cols = []
    if split_key in donor_table.columns:
        use_cols.append(split_key)
    use_cols.extend(cols)
    return (
        donor_table.assign(n_donors=1)
        .groupby(use_cols, observed=True, dropna=False)["n_donors"]
        .sum()
        .reset_index()
        .sort_values("n_donors", ascending=False)
    )


# =====================================================================
# Main
# =====================================================================
def run_audit(cfg: DonorPathologyAuditConfig) -> dict[str, Any]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    obs = load_obs(cfg)
    obs.index = obs.index.astype(str)

    donor_key = choose_existing_key(obs, cfg.donor_key, cfg.donor_key_fallbacks, role="donor")

    if cfg.drop_missing_donor:
        donor_vals = obs[donor_key].map(normalise_value)
        keep = donor_vals != "__NA__"
        if (~keep).sum() > 0:
            _log(f"dropping cells with missing donor: {(~keep).sum():,}")
        obs = obs.loc[keep].copy()

    if cfg.split_key not in obs.columns:
        warnings.warn(f"split_key={cfg.split_key!r} not found. Split-specific summaries will be skipped.")
    if cfg.celltype_key not in obs.columns:
        warnings.warn(f"celltype_key={cfg.celltype_key!r} not found. Celltype summaries will be limited.")

    _log(f"obs cells: {len(obs):,}")
    _log(f"donor key: {donor_key}")
    _log(f"unique donors: {obs[donor_key].astype(str).nunique():,}")

    donor_table, conflict_table = build_donor_table(
        obs,
        donor_key=donor_key,
        split_key=cfg.split_key,
        celltype_key=cfg.celltype_key,
        pathology_columns=cfg.pathology_columns,
        tech_columns=cfg.tech_columns,
    )
    save_df(donor_table, out_dir / "donor_pathology_table", cfg)

    if conflict_table.empty:
        conflict_table = pd.DataFrame(columns=["n_unique", "values"])
    save_df(conflict_table, out_dir / "donor_value_conflicts", cfg)

    split_summary = build_split_summary(donor_table, cfg.split_key)
    if not split_summary.empty:
        save_df(split_summary, out_dir / "split_summary", cfg)

    level_counts = build_pathology_level_counts(donor_table, pathology_columns=cfg.pathology_columns, split_key=cfg.split_key)
    if not level_counts.empty:
        save_df(level_counts, out_dir / "pathology_level_counts", cfg, index=False)

    train_warnings = build_train_level_warnings(
        level_counts,
        train_split_name=cfg.train_split_name,
        min_train_donors_per_level=cfg.min_train_donors_per_level,
    )
    if not train_warnings.empty:
        save_df(train_warnings, out_dir / "low_train_donor_pathology_levels", cfg, index=False)

    combo_counts = build_combination_counts(donor_table, combination_columns=cfg.combination_columns, split_key=cfg.split_key)
    if not combo_counts.empty:
        save_df(combo_counts, out_dir / "multi_pathology_combination_counts", cfg, index=False)

    crosstab_dir = out_dir / "crosstabs"
    crosstab_dir.mkdir(parents=True, exist_ok=True)
    crosstab_written = []
    for pair in cfg.crosstab_pairs:
        if len(pair) != 2:
            continue
        a, b = pair
        if a not in donor_table.columns or b not in donor_table.columns:
            continue
        tab = crosstab_with_totals(donor_table, a, b)
        fname = f"{safe_filename(a)}__x__{safe_filename(b)}.csv"
        tab.to_csv(crosstab_dir / fname)
        crosstab_written.append(fname)
        if cfg.split_key in donor_table.columns:
            for split, g in donor_table.groupby(cfg.split_key, observed=True, sort=True):
                tab_s = crosstab_with_totals(g, a, b)
                fname_s = f"{safe_filename(a)}__x__{safe_filename(b)}__split_{safe_filename(split)}.csv"
                tab_s.to_csv(crosstab_dir / fname_s)
                crosstab_written.append(fname_s)

    if cfg.celltype_key in obs.columns:
        base_cols = [donor_key, cfg.celltype_key]
        if cfg.split_key in obs.columns:
            base_cols.append(cfg.split_key)
        tmp = obs[base_cols].copy()
        tmp["_n_cells"] = 1
        group_cols = [cfg.celltype_key]
        if cfg.split_key in obs.columns:
            group_cols = [cfg.split_key, cfg.celltype_key]
        celltype_summary = (
            tmp.groupby(group_cols, observed=True)
            .agg(n_cells=("_n_cells", "sum"), n_donors=(donor_key, lambda x: x.astype(str).nunique()))
            .reset_index()
            .sort_values(["n_donors", "n_cells"], ascending=False)
        )
        save_df(celltype_summary, out_dir / "celltype_split_summary", cfg, index=False)

    summary = {
        "config": asdict(cfg),
        "n_cells": int(len(obs)),
        "n_donors": int(len(donor_table)),
        "donor_key_used": donor_key,
        "available_pathology_columns": [c for c in cfg.pathology_columns if c in donor_table.columns],
        "missing_pathology_columns": [c for c in cfg.pathology_columns if c not in donor_table.columns],
        "n_conflict_records": int(len(conflict_table)),
        "split_summary": split_summary.to_dict(orient="index") if not split_summary.empty else {},
        "n_low_train_donor_pathology_levels": int(len(train_warnings)),
        "crosstabs_written": crosstab_written,
        "outputs": {
            "donor_pathology_table": str(out_dir / "donor_pathology_table.csv"),
            "donor_value_conflicts": str(out_dir / "donor_value_conflicts.csv"),
            "pathology_level_counts": str(out_dir / "pathology_level_counts.csv"),
            "multi_pathology_combination_counts": str(out_dir / "multi_pathology_combination_counts.csv"),
            "crosstab_dir": str(crosstab_dir),
        },
    }

    with open(out_dir / "audit_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    _log("done")
    _log(f"main output: {out_dir / 'donor_pathology_table.csv'}")
    _log(f"summary: {out_dir / 'audit_summary.json'}")
    return summary


# =====================================================================
# CLI
# =====================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_audit(cfg)


if __name__ == "__main__":
    main()