from __future__ import annotations

import argparse
import gc
import json
import re
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# =====================================================================
# Config
# =====================================================================
@dataclass
class ResidualHVGConfig:
    h5ad_path: str
    out_dir: str

    celltype_key: str = "Subclass"
    donor_key_preferred: str = "subject_id"
    donor_key_fallback: str = "donor_id"
    counts_layer: str = "counts"

    n_top_genes: int = 4000
    max_cells_for_hvg: Optional[int] = 100000
    min_cells_per_gene: int = 5

    # Always subclass-specific HVG.
    # Only choose whether to use all subclasses or a filtered subset.
    use_all_subclasses: bool = True
    min_cells_per_subclass: int = 2000
    min_donors_per_subclass: int = 20

    save_obs_tables: bool = True
    obs_table_format: str = "parquet"  # parquet or csv
    save_full_subset_h5ad: bool = False

    seed: int = 42


# =====================================================================
# Logging / timing
# =====================================================================
def _log(msg: str) -> None:
    print(f"[residual-hvg] {msg}", flush=True)


class StageTimer:
    def __init__(self, name: str) -> None:
        self.name = name
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = time.time()
        _log(f"START: {self.name}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.time() - self.t0
        tag = "DONE" if exc_type is None else "FAIL"
        _log(f"{tag}: {self.name} ({dt:.1f}s)")
        return False


# =====================================================================
# Helpers
# =====================================================================
def safe_name(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"[^\w\-.]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def normalize_label(x: str) -> str:
    return str(x).strip()


def load_config(path: str | Path) -> ResidualHVGConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return ResidualHVGConfig(**raw)


def get_donor_key(obs: pd.DataFrame, preferred: str, fallback: str) -> Optional[str]:
    if preferred in obs.columns:
        return preferred
    if fallback in obs.columns:
        return fallback
    return None


def save_table(df: pd.DataFrame, path_base: Path, fmt: str) -> Path:
    fmt = str(fmt).strip().lower()
    if fmt == "parquet":
        out_path = path_base.with_suffix(".parquet")
        try:
            df.to_parquet(out_path, index=True)
            return out_path
        except Exception as e:
            warnings.warn(
                f"Parquet save failed for {out_path.name}; falling back to csv. reason={repr(e)}"
            )
    out_path = path_base.with_suffix(".csv")
    df.to_csv(out_path, index=True)
    return out_path


def backed_slice_to_memory(adata_b: ad.AnnData, row_idx: np.ndarray) -> ad.AnnData:
    view = adata_b[row_idx, :]
    if hasattr(view, "to_memory"):
        return view.to_memory()
    return view.copy()


def ensure_counts_layer(adata_sub: ad.AnnData, counts_layer: str = "counts") -> None:
    if counts_layer not in adata_sub.layers:
        adata_sub.layers[counts_layer] = adata_sub.X.copy()


def donor_balanced_subsample_local_positions(
    obs_sub: pd.DataFrame,
    max_cells: Optional[int],
    donor_key: Optional[str],
    seed: int = 42,
) -> np.ndarray:
    n = len(obs_sub)
    if max_cells is None or n <= max_cells:
        return np.arange(n, dtype=np.int64)

    rng = np.random.default_rng(seed)

    if donor_key is None or donor_key not in obs_sub.columns:
        return np.sort(rng.choice(n, size=max_cells, replace=False)).astype(np.int64)

    groups = obs_sub.groupby(donor_key, observed=True).indices
    total = float(n)

    selected = []
    for _, positions in groups.items():
        positions = np.asarray(positions, dtype=np.int64)
        k = max(1, int(round(len(positions) / total * max_cells)))
        k = min(k, len(positions))
        take = rng.choice(positions, size=k, replace=False)
        selected.append(take)

    selected = np.concatenate(selected) if len(selected) > 0 else np.array([], dtype=np.int64)
    selected = np.unique(selected)

    if len(selected) > max_cells:
        selected = np.sort(rng.choice(selected, size=max_cells, replace=False)).astype(np.int64)
    elif len(selected) < max_cells:
        remaining = np.setdiff1d(np.arange(n, dtype=np.int64), selected, assume_unique=False)
        extra_k = min(max_cells - len(selected), len(remaining))
        if extra_k > 0:
            extra = rng.choice(remaining, size=extra_k, replace=False)
            selected = np.sort(np.concatenate([selected, extra])).astype(np.int64)

    return np.sort(selected).astype(np.int64)


def prefilter_for_hvg(
    adata_sub: ad.AnnData,
    *,
    counts_layer: str = "counts",
    min_cells_per_gene: int = 5,
) -> ad.AnnData:
    X = adata_sub.layers[counts_layer] if counts_layer in adata_sub.layers else adata_sub.X

    if sparse.issparse(X):
        cell_sums = np.asarray(X.sum(axis=1)).ravel()
        gene_nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        cell_sums = np.asarray(X.sum(axis=1)).ravel()
        gene_nnz = np.asarray((X > 0).sum(axis=0)).ravel()

    cell_keep = np.isfinite(cell_sums) & (cell_sums > 0)
    gene_keep = np.isfinite(gene_nnz) & (gene_nnz >= min_cells_per_gene)

    return adata_sub[cell_keep, gene_keep].copy()


def compute_hvg_seurat_v3(
    adata_sub: ad.AnnData,
    *,
    n_top_genes: int = 4000,
    counts_layer: str = "counts",
) -> ad.AnnData:
    """
    HVG is always computed within subclass.
    We intentionally use batch_key=None
    """
    tmp = adata_sub.copy()
    ensure_counts_layer(tmp, counts_layer=counts_layer)

    sc.pp.highly_variable_genes(
        tmp,
        layer=counts_layer,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=None,
        subset=False,
        inplace=True,
    )
    tmp.uns["hvg_method_used"] = "seurat_v3"
    tmp.uns["hvg_batch_key"] = None
    return tmp


def build_summary_and_indices(
    obs_all: pd.DataFrame,
    *,
    celltype_key: str,
    donor_key: Optional[str],
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    if celltype_key not in obs_all.columns:
        raise KeyError(f"'{celltype_key}' not found in obs.")

    labels = obs_all[celltype_key].astype(str).map(normalize_label)

    work = pd.DataFrame(index=obs_all.index)
    work["subclass"] = labels.values

    if donor_key is not None and donor_key in obs_all.columns:
        work["donor"] = obs_all[donor_key].astype(str).values
        summary = (
            work.groupby("subclass", observed=True)
            .agg(
                n_cells=("subclass", "size"),
                n_donors=("donor", "nunique"),
            )
            .sort_values(["n_cells", "n_donors"], ascending=[False, False])
            .reset_index()
        )
    else:
        summary = (
            work.groupby("subclass", observed=True)
            .agg(n_cells=("subclass", "size"))
            .sort_values(["n_cells"], ascending=[False])
            .reset_index()
        )
        summary["n_donors"] = np.nan

    grouped = work.groupby("subclass", observed=True).indices
    subclass_to_idx = {
        str(k): np.asarray(v, dtype=np.int64)
        for k, v in grouped.items()
    }
    return summary, subclass_to_idx


def select_target_subclasses(summary: pd.DataFrame, cfg: ResidualHVGConfig) -> list[str]:
    if cfg.use_all_subclasses:
        targets = summary["subclass"].astype(str).tolist()
    else:
        targets = summary.loc[
            (summary["n_cells"] >= cfg.min_cells_per_subclass)
            & (summary["n_donors"].fillna(0) >= cfg.min_donors_per_subclass),
            "subclass",
        ].astype(str).tolist()
    return targets


def export_one_subclass_hvg_from_open_backed(
    adata_b: ad.AnnData,
    obs_all: pd.DataFrame,
    subclass_to_idx: dict[str, np.ndarray],
    subclass_name: str,
    out_dir: Path,
    cfg: ResidualHVGConfig,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    subclass_safe = safe_name(subclass_name)

    global_idx = subclass_to_idx.get(str(subclass_name), None)
    if global_idx is None or len(global_idx) == 0:
        raise ValueError(f"No cells found for subclass='{subclass_name}'.")

    obs_sub_full = obs_all.iloc[global_idx].copy()
    donor_key = get_donor_key(obs_sub_full, cfg.donor_key_preferred, cfg.donor_key_fallback)

    if cfg.save_obs_tables:
        obs_path = save_table(obs_sub_full, out_dir / f"{subclass_safe}_obs", cfg.obs_table_format)
        _log(f"[{subclass_name}] obs table saved: {obs_path.name}")

    local_keep = donor_balanced_subsample_local_positions(
        obs_sub=obs_sub_full,
        max_cells=cfg.max_cells_for_hvg,
        donor_key=donor_key,
        seed=cfg.seed,
    )
    global_keep = global_idx[local_keep]

    _log(f"[{subclass_name}] total cells = {len(global_idx):,}")
    _log(f"[{subclass_name}] cells used for HVG = {len(global_keep):,}")
    if donor_key is not None:
        _log(f"[{subclass_name}] donor key = {donor_key}, donors = {obs_sub_full[donor_key].astype(str).nunique():,}")

    with StageTimer(f"subset counts into memory for '{subclass_name}'"):
        adata_sub = backed_slice_to_memory(adata_b, global_keep)
        ensure_counts_layer(adata_sub, counts_layer=cfg.counts_layer)

    with StageTimer(f"prefilter for HVG '{subclass_name}'"):
        n_obs_before, n_vars_before = adata_sub.n_obs, adata_sub.n_vars
        adata_sub = prefilter_for_hvg(
            adata_sub,
            counts_layer=cfg.counts_layer,
            min_cells_per_gene=cfg.min_cells_per_gene,
        )
        _log(
            f"[{subclass_name}] prefilter: cells {n_obs_before:,} -> {adata_sub.n_obs:,}, "
            f"genes {n_vars_before:,} -> {adata_sub.n_vars:,}"
        )

    with StageTimer(f"compute HVG '{subclass_name}'"):
        tmp_hvg = compute_hvg_seurat_v3(
            adata_sub,
            n_top_genes=cfg.n_top_genes,
            counts_layer=cfg.counts_layer,
        )

    hvg_cols = [
        c
        for c in [
            "highly_variable",
            "highly_variable_rank",
            "means",
            "variances",
            "variances_norm",
            "highly_variable_nbatches",
            "highly_variable_intersection",
        ]
        if c in tmp_hvg.var.columns
    ]

    adata_sub.var = adata_sub.var.copy()
    for c in hvg_cols:
        adata_sub.var[c] = tmp_hvg.var[c].values

    adata_sub.uns["hvg_method_used"] = tmp_hvg.uns.get("hvg_method_used", "unknown")
    adata_sub.uns["hvg_batch_key"] = tmp_hvg.uns.get("hvg_batch_key", None)
    adata_sub.uns["hvg_n_top_genes"] = int(cfg.n_top_genes)
    adata_sub.uns["source_subclass"] = str(subclass_name)

    var_df = adata_sub.var.copy()
    var_df["gene_name"] = adata_sub.var_names.astype(str)
    if "highly_variable" not in var_df.columns:
        raise RuntimeError("HVG annotation was not created correctly.")

    hvg_df = var_df[var_df["highly_variable"].astype(bool)].copy()
    sort_cols = [c for c in ["highly_variable_rank", "variances_norm", "variances"] if c in hvg_df.columns]
    if len(sort_cols) > 0:
        ascending = [True] + [False] * (len(sort_cols) - 1)
        hvg_df = hvg_df.sort_values(sort_cols, ascending=ascending)

    hvg_table_path = out_dir / f"{subclass_safe}_hvg_table.csv"
    hvg_df.to_csv(hvg_table_path, index=True)

    hvg_list_path = out_dir / f"{subclass_safe}_hvg_list.json"
    with open(hvg_list_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "subclass": str(subclass_name),
                "n_cells_total": int(len(global_idx)),
                "n_cells_used_for_hvg": int(len(global_keep)),
                "n_cells_after_prefilter": int(adata_sub.n_obs),
                "n_genes_after_prefilter": int(adata_sub.n_vars),
                "n_hvg": int(hvg_df.shape[0]),
                "hvg_method_used": adata_sub.uns["hvg_method_used"],
                "hvg_batch_key": adata_sub.uns["hvg_batch_key"],
                "genes": hvg_df["gene_name"].astype(str).tolist(),
            },
            f,
            ensure_ascii=False,
            indent=2,
        )

    with StageTimer(f"write HVG-only h5ad '{subclass_name}'"):
        hvg_mask = adata_sub.var["highly_variable"].astype(bool).to_numpy()
        adata_hvg = adata_sub[:, hvg_mask].copy()
        ensure_counts_layer(adata_hvg, counts_layer=cfg.counts_layer)
        hvg_h5ad_path = out_dir / f"{subclass_safe}_hvg{cfg.n_top_genes}.h5ad"
        adata_hvg.write_h5ad(hvg_h5ad_path)

    _log(f"[{subclass_name}] HVG-only h5ad saved: {hvg_h5ad_path}")
    _log(f"[{subclass_name}] shape = {adata_hvg.shape}")

    full_subset_path = None
    if cfg.save_full_subset_h5ad:
        with StageTimer(f"write full filtered subset h5ad '{subclass_name}'"):
            full_subset_path = out_dir / f"{subclass_safe}_full_subset.h5ad"
            adata_sub.write_h5ad(full_subset_path)
        _log(f"[{subclass_name}] full subset h5ad saved: {full_subset_path}")

    result = {
        "subclass": str(subclass_name),
        "subclass_safe": subclass_safe,
        "n_cells_total": int(len(global_idx)),
        "n_cells_used_for_hvg": int(len(global_keep)),
        "n_cells_after_prefilter": int(adata_sub.n_obs),
        "n_genes_after_prefilter": int(adata_sub.n_vars),
        "n_hvg": int(hvg_df.shape[0]),
        "hvg_method_used": adata_sub.uns["hvg_method_used"],
        "obs_table": None if not cfg.save_obs_tables else str((out_dir / f"{subclass_safe}_obs").with_suffix(".parquet" if cfg.obs_table_format.lower() == "parquet" else ".csv")),
        "hvg_table": str(hvg_table_path),
        "hvg_list_json": str(hvg_list_path),
        "hvg_h5ad": str(hvg_h5ad_path),
        "full_subset_h5ad": None if full_subset_path is None else str(full_subset_path),
    }

    del tmp_hvg, adata_sub, adata_hvg, obs_sub_full
    gc.collect()
    return result


# =====================================================================
# Main
# =====================================================================
def run(cfg: ResidualHVGConfig) -> None:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    obs_dir = out_dir / "subclass_obs"
    hvg_dir = out_dir / "subclass_hvg"
    obs_dir.mkdir(parents=True, exist_ok=True)
    hvg_dir.mkdir(parents=True, exist_ok=True)

    _log(f"output directory: {out_dir}")
    _log(f"subclass mode: {'all subclasses' if cfg.use_all_subclasses else 'filtered subclasses only'}")

    adata_b: Optional[ad.AnnData] = None
    results: list[dict] = []

    try:
        with StageTimer("open residual h5ad once in backed mode"):
            adata_b = ad.read_h5ad(cfg.h5ad_path, backed="r")
            obs_all = adata_b.obs.copy()
            _log(f"obs shape = {obs_all.shape}")

        donor_key = get_donor_key(obs_all, cfg.donor_key_preferred, cfg.donor_key_fallback)
        _log(f"preferred donor schema = {donor_key}")

        with StageTimer("build subclass summary and indices once"):
            summary, subclass_to_idx = build_summary_and_indices(
                obs_all,
                celltype_key=cfg.celltype_key,
                donor_key=donor_key,
            )
            summary_path = out_dir / "subclass_summary.csv"
            summary.to_csv(summary_path, index=False)
            _log(f"subclass summary saved: {summary_path}")

        targets = select_target_subclasses(summary, cfg)
        selected_df = summary[summary["subclass"].isin(targets)].copy()
        selected_path = out_dir / "selected_subclasses.csv"
        selected_df.to_csv(selected_path, index=False)
        _log(f"selected subclasses = {len(targets):,}")
        _log(f"selected subclass table saved: {selected_path}")

        if len(targets) == 0:
            raise ValueError("No subclasses were selected. Relax thresholds or set use_all_subclasses=True.")

        for i, subclass_name in enumerate(targets, start=1):
            _log("-" * 72)
            _log(f"[{i}/{len(targets)}] processing subclass = {subclass_name}")
            res = export_one_subclass_hvg_from_open_backed(
                adata_b=adata_b,
                obs_all=obs_all,
                subclass_to_idx=subclass_to_idx,
                subclass_name=subclass_name,
                out_dir=hvg_dir / safe_name(subclass_name),
                cfg=cfg,
            )
            results.append(res)

        results_df = pd.DataFrame(results)
        results_csv = out_dir / "results_manifest.csv"
        results_df.to_csv(results_csv, index=False)
        _log(f"results manifest saved: {results_csv}")

        manifest = {
            "config": asdict(cfg),
            "n_total_cells": int(len(obs_all)),
            "n_subclasses_total": int(summary.shape[0]),
            "n_subclasses_selected": int(len(targets)),
            "selected_subclasses": targets,
            "summary_csv": str(summary_path),
            "selected_subclasses_csv": str(selected_path),
            "results_manifest_csv": str(results_csv),
        }
        manifest_path = out_dir / "manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)
        _log(f"manifest saved: {manifest_path}")

        _log("ALL DONE")

    finally:
        if adata_b is not None:
            try:
                adata_b.file.close()
                _log("closed backed h5ad handle")
            except Exception:
                pass


# =====================================================================
# CLI
# =====================================================================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Extract subclass-specific HVGs from residual_train.h5ad. "
            "HVGs are always computed within each subclass. "
            "You may choose either all subclasses or only subclasses that pass "
            "cell/donor thresholds."
        )
    )
    p.add_argument("--config", type=str, required=True, help="Path to JSON config.")
    return p


def main() -> None:
    args = build_argparser().parse_args()
    cfg = load_config(args.config)
    run(cfg)


if __name__ == "__main__":
    main()