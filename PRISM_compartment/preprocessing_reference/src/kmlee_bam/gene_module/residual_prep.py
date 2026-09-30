from __future__ import annotations

"""
Residual-gene preparation for hdWGCNA.

What this script does
---------------------
1) Reads the residual-gene list from an existing registry JSON.
2) Uses module_h5ad.obs as the authoritative cell universe for this stage.
3) Maps those cells back to the source zarr robustly by barcode (not just by order).
4) Audits donor × disease confounding and screens candidate tech covariates.
5) Performs donor-level train/val/test split using donor_split from preprocess_scRNA_01.py.
6) Emits subject_id / subject_int_id / tech_id_model schema.
7) Writes residual_train.h5ad (and optionally val/test) as raw-count sparse matrices.

Design principles
-----------------
- Residual genes are a feature universe, not a cell filter.
- Only residual_train.h5ad is allowed to feed hdWGCNA module discovery.
- donor_id is never used as tech_id_model.
- Auto-picking a tech covariate is conservative: subject-like columns such as
  sample_id are blocked by default, and suspicious high-cardinality keys are
  rejected from auto-pick.
- X always stores raw counts; layers['counts'] is optional and controlled by split.
- Optional obs sidecars can provide an externally fixed donor split. This is
  used when module discovery must share the exact same train/val/test donors as
  a downstream model dataset.
"""

import argparse
import inspect
import json
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr


# =====================================================================
# Config
# =====================================================================
@dataclass
class ResidualPrepConfig:
    # ---- inputs -----------------------------------------------------
    zarr_path: str
    registry_json_path: str
    module_h5ad_path: str

    # ---- outputs ----------------------------------------------------
    out_dir: str

    # ---- optional externally prepared obs / split sidecar -----------
    obs_parquet_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_h5ad_path: Optional[str] = None
    merge_obs_sidecar: bool = True

    # ---- obs column names ------------------------------------------
    donor_key: str = "donor_id"
    disease_key: str = "disease"
    celltype_key: str = "Subclass"
    primary_label_key: Optional[str] = None

    # ---- split parameters ------------------------------------------
    seed: int = 42
    train_frac: float = 0.85
    val_frac: float = 0.10
    split_mode: str = "pretrain"  # "pretrain" or "benchmark"
    min_donors_per_stratum: int = 3
    drop_mixed_in_benchmark: bool = False
    use_existing_split: bool = False
    existing_split_key: str = "split"

    # ---- extraction -------------------------------------------------
    max_gap: int = 3_000_000
    chunk_rows_fallback: int = 4096
    skip_val_test_extraction: bool = True
    write_counts_layer_train: bool = True
    write_counts_layer_val_test: bool = False
    write_counts_layer: Optional[bool] = None  # deprecated global override

    # ---- batch/tech audit ------------------------------------------
    batch_candidate_keys: List[str] = field(
        default_factory=lambda: [
            "Library prep",
            "library_prep",
            "Single Nucleus Sequencing batch",
            "sequencing_batch",
            "Batch",
            "batch",
            "10x chemistry",
            "chemistry",
            "assay",
            "platform",
            "sample_id",
        ]
    )
    tech_mixing_threshold: float = 0.5
    tech_min_median_donors_per_level: float = 2.0
    tech_max_levels_fraction_of_donors: float = 0.50
    auto_pick_priority: List[str] = field(
        default_factory=lambda: [
            "Library prep",
            "library_prep",
            "Single Nucleus Sequencing batch",
            "sequencing_batch",
            "Batch",
            "batch",
            "10x chemistry",
            "chemistry",
            "assay",
            "platform",
            "sample_id",
        ]
    )
    auto_pick_blocklist: List[str] = field(default_factory=lambda: ["sample_id"])
    fail_if_no_safe_tech: bool = False

    # ---- imports ----------------------------------------------------
    preprocess_module_dir: Optional[str] = None


# =====================================================================
# Logging
# =====================================================================
def _log(msg: str) -> None:
    print(f"[residual-prep] {msg}", flush=True)


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
# Small helpers
# =====================================================================
def _unique_split_or_multi(s: pd.Series) -> str:
    vals = pd.unique(s.astype(str))
    return str(vals[0]) if len(vals) == 1 else "__MULTI__"


def _sorted_non_null_str_tuple(s: pd.Series) -> tuple[str, ...]:
    return tuple(sorted({str(v) for v in s.dropna()}))


def _sorted_str_tuple(s: pd.Series) -> tuple[str, ...]:
    return tuple(sorted(set(map(str, s))))


def load_obs_sidecar(cfg: ResidualPrepConfig) -> Optional[pd.DataFrame]:
    n_sources = sum(x is not None for x in [cfg.obs_parquet_path, cfg.obs_csv_path, cfg.obs_h5ad_path])
    if n_sources > 1:
        raise ValueError("Provide at most one of obs_parquet_path, obs_csv_path, obs_h5ad_path.")
    if cfg.obs_parquet_path is not None:
        _log(f"loading obs sidecar parquet: {cfg.obs_parquet_path}")
        return pd.read_parquet(cfg.obs_parquet_path)
    if cfg.obs_csv_path is not None:
        _log(f"loading obs sidecar csv: {cfg.obs_csv_path}")
        return pd.read_csv(cfg.obs_csv_path, index_col=0)
    if cfg.obs_h5ad_path is not None:
        _log(f"loading obs sidecar h5ad: {cfg.obs_h5ad_path}")
        return ad.read_h5ad(cfg.obs_h5ad_path, backed="r").obs.copy()
    return None


def merge_obs_sidecar(obs: pd.DataFrame, sidecar: Optional[pd.DataFrame], *, merge: bool) -> pd.DataFrame:
    obs = obs.copy()
    obs.index = obs.index.astype(str)
    if sidecar is None:
        return obs

    sidecar = sidecar.copy()
    sidecar.index = sidecar.index.astype(str)
    if sidecar.index.has_duplicates:
        dup = sidecar.index[sidecar.index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"obs sidecar has duplicated barcodes. Examples: {dup}")

    missing = pd.Index(obs.index).difference(pd.Index(sidecar.index))
    if len(missing) > 0:
        raise ValueError(f"module_h5ad cells missing from obs sidecar. Examples: {missing[:5].tolist()}")

    aligned = sidecar.reindex(obs.index)
    if not merge:
        return aligned

    out = obs.copy()
    for c in aligned.columns:
        out[c] = aligned[c].values
    return out


# =====================================================================
# Helpers: import / registry / zarr
# =====================================================================
def import_preprocess_helpers(
    preprocess_module_dir: Optional[str],
) -> tuple[Callable[..., Any], Optional[Callable[..., sp.csr_matrix]]]:
    if preprocess_module_dir:
        sys.path.append(str(preprocess_module_dir))
    try:
        import scRNA_preprocess_01 as prep  # type: ignore
    except ImportError as e:
        raise ImportError(
            "Could not import preprocess_scRNA_01. Set preprocess_module_dir or PYTHONPATH."
        ) from e

    donor_split = getattr(prep, "donor_split", None)
    if donor_split is None:
        raise ImportError("preprocess_scRNA_01.py does not expose donor_split.")

    extract_fn = getattr(prep, "extract_csr_subset_from_zarr", None)
    return donor_split, extract_fn


def load_registry_residual_genes(
    registry_json_path: str | Path,
) -> tuple[List[str], List[int], int]:
    with open(registry_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    gene_names = data["gene_names"]
    residual_indices = data.get("residual_gene_indices", [])
    if len(residual_indices) == 0:
        raise ValueError(
            "No residual_gene_indices found in registry. "
            "Did the build script run with add_residual_module=True?"
        )

    residual_names = [gene_names[i] for i in residual_indices]
    return residual_names, residual_indices, len(gene_names)


def _normalize_ensembl_index(idx: pd.Index) -> pd.Index:
    return idx.astype(str).str.replace(r"\.\d+$", "", regex=True).str.upper()


def _read_zarr_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError("Could not read AnnData-Zarr metadata element.") from e


def _read_zarr_dataframe(root, key_path: str) -> pd.DataFrame:
    group = root
    for part in key_path.split("/"):
        if part:
            group = group[part]
    out = _read_zarr_elem(group)
    if not isinstance(out, pd.DataFrame):
        raise TypeError(f"Zarr element {key_path!r} is not a pandas DataFrame.")
    return out.copy()


def _read_obs_var_names_from_zarr(zarr_path: str | Path) -> tuple[pd.DataFrame, pd.Index, bool]:
    root = zarr.open(str(zarr_path), mode="r")
    obs = _read_zarr_dataframe(root, "obs")
    if "raw" in root and "var" in root["raw"]:
        var = _read_zarr_dataframe(root, "raw/var")
        use_raw = True
    else:
        var = _read_zarr_dataframe(root, "var")
        use_raw = False
    return obs, pd.Index(var.index.astype(str)), use_raw


def _choose_matrix_and_var(zarr_path: str | Path):
    return _read_obs_var_names_from_zarr(zarr_path)


# =====================================================================
# Alignment
# =====================================================================
def build_module_obs_to_zarr_positions(
    zarr_path: str | Path,
    module_h5ad_obs: pd.DataFrame,
) -> np.ndarray:
    """
    Return zarr row positions aligned to module_h5ad_obs rows.

    Fast path:
      - exact same length and exact same barcode order
    Fallback:
      - explicit barcode lookup into zarr obs index
    """
    zarr_obs, _, _ = _choose_matrix_and_var(zarr_path)
    zarr_index = pd.Index(zarr_obs.index.astype(str))
    module_index = pd.Index(module_h5ad_obs.index.astype(str))

    if len(zarr_index) == len(module_index) and np.array_equal(
        zarr_index.to_numpy(), module_index.to_numpy()
    ):
        _log("alignment: exact 1-to-1 barcode order match")
        return np.arange(len(module_index), dtype=np.int64)

    if module_index.has_duplicates:
        dup = module_index[module_index.duplicated()].unique().tolist()[:5]
        raise RuntimeError(f"module_h5ad has duplicated barcodes. Examples: {dup}")

    row_positions = zarr_index.get_indexer(module_index)
    missing_mask = row_positions < 0
    if missing_mask.any():
        missing = module_index[missing_mask].tolist()[:5]
        raise RuntimeError(
            "Some module_h5ad barcodes were not found in zarr. "
            f"Examples: {missing}"
        )

    _log(
        "alignment: explicit barcode mapping used "
        f"({len(module_index):,} module cells mapped into {len(zarr_index):,} zarr rows)"
    )
    return row_positions.astype(np.int64)


# =====================================================================
# Confounding / tech selection
# =====================================================================
def diagnose_batch_confounding(
    obs: pd.DataFrame,
    *,
    donor_key: str,
    disease_key: str,
    candidate_batch_keys: Sequence[str],
    mixing_threshold: float,
    min_median_donors_per_level: float,
    max_levels_fraction_of_donors: float,
    donor_key_forbid: Optional[str] = None,
) -> dict:
    report: dict = {}

    if donor_key not in obs.columns or disease_key not in obs.columns:
        raise KeyError(f"Missing required columns: {donor_key}, {disease_key}")

    donor_disease = (
        obs.groupby(donor_key, observed=True)[disease_key]
        .apply(_sorted_non_null_str_tuple)
        .rename("disease_set")
        .reset_index()
    )
    donor_disease["n_diseases"] = donor_disease["disease_set"].map(len)
    n_donors = int(len(donor_disease))
    n_donors_multi = int(donor_disease["n_diseases"].gt(1).sum())
    donor_purity = float((donor_disease["n_diseases"] == 1).mean()) if n_donors > 0 else 0.0

    report["n_donors"] = n_donors
    report["n_donors_with_multiple_diseases"] = n_donors_multi
    report["donor_disease_purity"] = donor_purity
    report["donor_perfectly_confounded_with_disease"] = donor_purity >= 0.999
    report["disease_counts"] = {
        str(k): int(v) for k, v in obs[disease_key].astype(str).value_counts().to_dict().items()
    }

    donor_disease_long = (
        donor_disease[[donor_key, "disease_set"]]
        .explode("disease_set")
        .rename(columns={"disease_set": disease_key})
    )
    donor_disease_long[disease_key] = donor_disease_long[disease_key].fillna("__NA__").astype(str)

    cand_report: Dict[str, dict] = {}
    for key in candidate_batch_keys:
        entry: dict = {"present": key in obs.columns}
        if not entry["present"]:
            cand_report[key] = entry
            continue

        if donor_key_forbid is not None and key == donor_key_forbid:
            entry["reasonable_tech_covariate"] = False
            entry["reason"] = "forbidden: donor key"
            cand_report[key] = entry
            continue

        levels = obs[key].astype(str)
        entry["n_levels"] = int(levels.nunique(dropna=False))
        if entry["n_levels"] <= 1:
            entry["reasonable_tech_covariate"] = False
            entry["reason"] = "only 1 level"
            cand_report[key] = entry
            continue

        # 1) cell-level disease mixing
        ct_dis_cell = pd.crosstab(levels, obs[disease_key].astype(str))
        levels_multi_disease_cell = int(ct_dis_cell.gt(0).sum(axis=1).gt(1).sum())
        frac_cell_multi_disease = float(levels_multi_disease_cell / max(entry["n_levels"], 1))

        # 2) donor-level disease mixing
        level_donor = (
            obs[[key, donor_key]]
            .drop_duplicates()
            .merge(donor_disease_long, on=donor_key, how="left")
        )
        ct_dis_donor = pd.crosstab(
            level_donor[key].astype(str),
            level_donor[disease_key].astype(str),
        )
        levels_multi_disease_donor = int(ct_dis_donor.gt(0).sum(axis=1).gt(1).sum())
        frac_donor_multi_disease = float(levels_multi_disease_donor / max(entry["n_levels"], 1))

        # 3) donor replication / subject-like cardinality
        ct_donor = pd.crosstab(levels, obs[donor_key].astype(str))
        donors_per_level = ct_donor.gt(0).sum(axis=1)
        median_donors_per_level = float(np.median(donors_per_level.to_numpy()))
        mean_donors_per_level = float(np.mean(donors_per_level.to_numpy()))
        levels_fraction_of_donors = float(entry["n_levels"] / max(n_donors, 1))

        passes_cell_mixing = frac_cell_multi_disease >= mixing_threshold
        passes_donor_mixing = frac_donor_multi_disease >= mixing_threshold
        passes_mixing = passes_cell_mixing and passes_donor_mixing
        passes_replication = median_donors_per_level >= min_median_donors_per_level
        passes_cardinality = levels_fraction_of_donors <= max_levels_fraction_of_donors

        entry.update(
            {
                # backward compatibility
                "n_levels_covering_multiple_diseases": levels_multi_disease_cell,
                "fraction_levels_multi_disease": frac_cell_multi_disease,
                # detailed diagnostics
                "n_levels_covering_multiple_diseases_cell": levels_multi_disease_cell,
                "fraction_levels_multi_disease_cell": frac_cell_multi_disease,
                "n_levels_covering_multiple_diseases_donor": levels_multi_disease_donor,
                "fraction_levels_multi_disease_donor": frac_donor_multi_disease,
                "median_donors_per_level": median_donors_per_level,
                "mean_donors_per_level": mean_donors_per_level,
                "levels_fraction_of_donors": levels_fraction_of_donors,
                "levels_preview": sorted(map(str, pd.unique(levels)))[:6],
                "passes_cell_mixing": passes_cell_mixing,
                "passes_donor_mixing": passes_donor_mixing,
                "passes_mixing": passes_mixing,
                "passes_replication": passes_replication,
                "passes_cardinality": passes_cardinality,
            }
        )
        entry["reasonable_tech_covariate"] = (
            passes_mixing and passes_replication and passes_cardinality
        )

        reasons: List[str] = []
        if not passes_cell_mixing:
            reasons.append("poor cell-level disease mixing")
        if not passes_donor_mixing:
            reasons.append("poor donor-level disease mixing")
        if not passes_replication:
            reasons.append("too few donors per level")
        if not passes_cardinality:
            reasons.append("too many levels relative to donors (subject-like)")
        entry["reason"] = "; ".join(reasons) if reasons else "ok"
        cand_report[key] = entry

    report["batch_candidates"] = cand_report
    usable = [k for k, v in cand_report.items() if v.get("reasonable_tech_covariate", False)]
    report["recommended_tech_covariates"] = usable
    report["recommendation"] = (
        "Use a recommended tech covariate for tech_id_model; never use donor_id."
        if usable
        else "No safe tech covariate was found. Use '__GLOBAL__' and disable or heavily shrink the tech term."
    )
    return report


def print_diagnosis(report: dict) -> None:
    _log(f"n_donors = {report['n_donors']}")
    _log(f"donor disease purity = {report['donor_disease_purity']:.4f}")
    _log(
        "donor perfectly confounded with disease: "
        f"{report['donor_perfectly_confounded_with_disease']}"
    )
    _log(f"disease counts: {report['disease_counts']}")
    _log("candidate tech covariates:")
    for k, v in report["batch_candidates"].items():
        if not v.get("present", False):
            _log(f"  {k:35s} : (not present)")
            continue
        status = "OK" if v.get("reasonable_tech_covariate") else "NO"
        _log(
            f"  {k:35s} : [{status}] n_levels={v.get('n_levels')}, "
            f"cell-mix={v.get('fraction_levels_multi_disease_cell', v.get('fraction_levels_multi_disease', 0)):.2f}, "
            f"donor-mix={v.get('fraction_levels_multi_disease_donor', 0):.2f}, "
            f"median donors/level={v.get('median_donors_per_level', 0):.1f}, "
            f"levels/donors={v.get('levels_fraction_of_donors', 0):.2f}, "
            f"reason={v.get('reason')}"
        )
    _log(f"recommended tech covariates: {report['recommended_tech_covariates']}")
    _log(f"note: {report['recommendation']}")


# =====================================================================
# Tech covariate selection
# =====================================================================
def choose_tech_covariate(
    obs: pd.DataFrame,
    report: dict,
    *,
    donor_key: str,
    preferred_tech_covariate: Optional[str],
    auto_pick_priority: Sequence[str],
    auto_pick_blocklist: Sequence[str],
    fail_if_no_safe_tech: bool,
    allow_unsafe_preferred: bool = False,
) -> tuple[pd.DataFrame, Optional[str]]:
    obs = obs.copy()
    obs["subject_id"] = obs[donor_key].astype(str)
    obs["subject_int_id"] = obs["subject_id"].astype("category").cat.codes.astype(np.int64)

    if preferred_tech_covariate is not None:
        if preferred_tech_covariate == donor_key:
            raise ValueError("preferred tech covariate cannot be donor_key.")
        if preferred_tech_covariate not in obs.columns:
            raise ValueError(f"preferred tech covariate '{preferred_tech_covariate}' not found in obs.")

        cand = report["batch_candidates"].get(preferred_tech_covariate, {})
        safe = cand.get("reasonable_tech_covariate", False)
        if (not safe) and (not allow_unsafe_preferred):
            raise ValueError(
                f"preferred tech covariate '{preferred_tech_covariate}' failed safety screening: {cand.get('reason')}"
            )
        obs["tech_id_model"] = obs[preferred_tech_covariate].astype(str)
        return obs, preferred_tech_covariate

    usable = set(report.get("recommended_tech_covariates", []))
    blocked = set(auto_pick_blocklist) | {donor_key}
    chosen: Optional[str] = None
    for key in auto_pick_priority:
        if key in blocked:
            continue
        if key in usable and key in obs.columns:
            chosen = key
            break

    if chosen is None:
        remaining = [k for k in usable if k not in blocked and k in obs.columns]
        if remaining:
            chosen = sorted(remaining)[0]

    if chosen is None:
        if fail_if_no_safe_tech:
            raise RuntimeError(
                "No safe tech covariate was found. Inspect obs columns manually or set fail_if_no_safe_tech=False."
            )
        obs["tech_id_model"] = "__GLOBAL__"
        return obs, None

    obs["tech_id_model"] = obs[chosen].astype(str)
    return obs, chosen


# =====================================================================
# donor_split wrapper
# =====================================================================
def call_donor_split(
    donor_split_fn: Callable[..., Any],
    adata_obs: ad.AnnData,
    cfg: ResidualPrepConfig,
) -> ad.AnnData:
    kwargs = {
        "donor_key": cfg.donor_key,
        "disease_key": cfg.disease_key,
        "primary_label_key": cfg.primary_label_key,
        "mode": cfg.split_mode,
        "seed": cfg.seed,
        "train_frac": cfg.train_frac,
        "val_frac": cfg.val_frac,
        "min_donors_per_stratum": cfg.min_donors_per_stratum,
        "drop_mixed_in_benchmark": cfg.drop_mixed_in_benchmark,
    }
    sig = inspect.signature(donor_split_fn)
    filtered = {k: v for k, v in kwargs.items() if k in sig.parameters and v is not None}
    return donor_split_fn(adata_obs, **filtered)


# =====================================================================
# Counts-layer policy
# =====================================================================
def resolve_write_counts_layer_for_split(
    cfg: ResidualPrepConfig,
    split_name: str,
) -> bool:
    """
    Decide whether to explicitly write layers['counts'] for a given split.

    Priority
    --------
    1) Deprecated global override:
         cfg.write_counts_layer is not None
       -> use that value for every split.

    2) Split-specific policy:
         - train -> cfg.write_counts_layer_train
         - val/test -> cfg.write_counts_layer_val_test
    """
    split = str(split_name).strip().lower()

    if cfg.write_counts_layer is not None:
        return bool(cfg.write_counts_layer)
    if split == "train":
        return bool(cfg.write_counts_layer_train)
    if split in {"val", "test"}:
        return bool(cfg.write_counts_layer_val_test)

    raise ValueError(
        f"Unknown split_name='{split_name}'. Expected one of: 'train', 'val', 'test'."
    )


# =====================================================================
# Extraction
# =====================================================================
def _fallback_extract_rows(
    z_X,
    orig_row_positions: np.ndarray,
    n_vars: int,
    chunk_rows: int,
) -> sp.csr_matrix:
    if len(orig_row_positions) == 0:
        return sp.csr_matrix((0, n_vars), dtype=np.float32)

    sort_idx = np.argsort(orig_row_positions)
    sorted_rows = orig_row_positions[sort_idx]
    inv = np.argsort(sort_idx)

    data_arr = z_X["data"]
    indices_arr = z_X["indices"]
    indptr_arr = z_X["indptr"]

    out_blocks: List[sp.csr_matrix] = []
    n = len(sorted_rows)
    for start in range(0, n, chunk_rows):
        stop = min(n, start + chunk_rows)
        block_rows = sorted_rows[start:stop]
        row_datas, row_indices, row_lengths = [], [], []
        for r in block_rows:
            p0 = int(indptr_arr[r])
            p1 = int(indptr_arr[r + 1])
            row_datas.append(np.asarray(data_arr[p0:p1], dtype=np.float32))
            row_indices.append(np.asarray(indices_arr[p0:p1], dtype=np.int64))
            row_lengths.append(p1 - p0)
        data = np.concatenate(row_datas) if row_datas else np.array([], dtype=np.float32)
        indices = np.concatenate(row_indices) if row_indices else np.array([], dtype=np.int64)
        indptr = np.zeros(len(block_rows) + 1, dtype=np.int64)
        np.cumsum(row_lengths, out=indptr[1:])
        out_blocks.append(
            sp.csr_matrix((data, indices, indptr), shape=(len(block_rows), n_vars), dtype=np.float32)
        )
    return sp.vstack(out_blocks, format="csr")[inv]


def extract_residual_counts_for_cells(
    *,
    zarr_path: str | Path,
    residual_gene_names: Sequence[str],
    row_positions: np.ndarray,
    extract_csr_subset_from_zarr: Optional[Callable[..., sp.csr_matrix]],
    max_gap: int,
    chunk_rows_fallback: int,
) -> sp.csr_matrix:
    z = zarr.open(str(zarr_path), mode="r")
    _, var_names, use_raw = _choose_matrix_and_var(zarr_path)
    xgrp = z["raw"]["X"] if use_raw else z["X"]

    var_norm = _normalize_ensembl_index(var_names)
    var_to_idx = {str(v): i for i, v in enumerate(var_norm.tolist())}

    col_idx, missing = [], []
    for g in residual_gene_names:
        raw_key = str(g).strip()
        key = raw_key.split(".")[0].upper() if raw_key.upper().startswith("ENSG") else raw_key.upper()
        if key in var_to_idx:
            col_idx.append(var_to_idx[key])
        else:
            missing.append(g)
    if missing:
        raise KeyError(
            f"{len(missing)} residual genes not found in zarr var_names. Examples: {missing[:5]}"
        )
    col_idx = np.asarray(col_idx, dtype=np.int64)

    if extract_csr_subset_from_zarr is not None:
        X_rows = extract_csr_subset_from_zarr(
            xgrp,
            row_positions.astype(np.int64),
            len(var_names),
            max_gap=max_gap,
        )
    else:
        X_rows = _fallback_extract_rows(
            xgrp,
            row_positions.astype(np.int64),
            len(var_names),
            chunk_rows_fallback,
        )

    return X_rows[:, col_idx].tocsr()


# =====================================================================
# Reporting
# =====================================================================
def write_text_report(
    out_path: Path,
    *,
    cfg: ResidualPrepConfig,
    report: dict,
    residual_n: int,
    registry_n: int,
    n_cells_total: int,
    split_counts: Dict[str, int],
    tech_key_chosen: Optional[str],
    counts_layer_policy: Dict[str, bool],
) -> None:
    lines: List[str] = []
    lines.append("Residual hdWGCNA preparation report")
    lines.append("=" * 72)
    lines.append("")
    lines.append("[inputs]")
    lines.append(f"  zarr            : {cfg.zarr_path}")
    lines.append(f"  registry        : {cfg.registry_json_path}")
    lines.append(f"  module h5ad     : {cfg.module_h5ad_path}")
    lines.append("")
    lines.append("[gene universe]")
    lines.append(f"  total registry genes : {registry_n:,}")
    lines.append(f"  residual genes       : {residual_n:,}")
    lines.append("  note: residual genes are a feature universe, not a cell filter.")
    lines.append("")
    lines.append("[cells]")
    lines.append(f"  total cells : {n_cells_total:,}")
    for s in ("train", "val", "test"):
        lines.append(f"  {s:5s}     : {split_counts.get(s, 0):,}")
    lines.append("")
    lines.append("[confounding audit]")
    lines.append(f"  n_donors             : {report['n_donors']}")
    lines.append(f"  donor disease purity : {report['donor_disease_purity']:.4f}")
    lines.append(f"  donor ≡ disease      : {report['donor_perfectly_confounded_with_disease']}")
    lines.append("")
    lines.append("[tech covariate selection]")
    lines.append(f"  recommended : {report['recommended_tech_covariates']}")
    lines.append(f"  chosen      : {tech_key_chosen if tech_key_chosen else '__GLOBAL__'}")
    lines.append(f"  note        : {report['recommendation']}")
    lines.append("")
    lines.append("[schema emitted in output h5ad obs]")
    lines.append("  subject_id      : donor_id as string (biology / subject)")
    lines.append("  subject_int_id  : integer code of subject_id")
    lines.append("  tech_id_model   : chosen safe tech covariate, or '__GLOBAL__'")
    lines.append("  split           : train / val / test")
    lines.append("")
    lines.append("[counts storage policy]")
    lines.append(f"  train counts layer : {counts_layer_policy.get('train', False)}")
    lines.append(f"  val counts layer   : {counts_layer_policy.get('val', False)}")
    lines.append(f"  test counts layer  : {counts_layer_policy.get('test', False)}")
    lines.append("  note: X always stores raw counts; layers['counts'] is written")
    lines.append("        explicitly only when the policy for that split is True.")
    lines.append("")
    lines.append("[next stages]")
    lines.append("  Stage 3 : hdWGCNA on residual_train.h5ad ONLY")
    lines.append("  Stage 4 : integrate discovered sub-modules into registry")
    lines.append("  Stage 5 : rebuild module_encoder_input for all splits")
    out_path.write_text("\n".join(lines), encoding="utf-8")


# =====================================================================
# Main driver
# =====================================================================
def prepare_residual_for_hdwgcna(
    cfg: ResidualPrepConfig,
    *,
    preferred_tech_covariate: Optional[str] = None,
    allow_unsafe_preferred_tech: bool = False,
) -> None:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"output directory: {out_dir}")

    donor_split_fn, extract_csr_subset_from_zarr = import_preprocess_helpers(cfg.preprocess_module_dir)

    with StageTimer("load registry residual genes"):
        residual_names, residual_indices, registry_n = load_registry_residual_genes(cfg.registry_json_path)
        _log(f"total registry genes = {registry_n:,}")
        _log(f"residual genes       = {len(residual_names):,}")
        with open(out_dir / "residual_gene_list.json", "w", encoding="utf-8") as f:
            json.dump(
                {
                    "n_total_registry_genes": registry_n,
                    "n_residual_genes": len(residual_names),
                    "residual_gene_names": residual_names,
                    "registry_indices": residual_indices,
                },
                f,
                ensure_ascii=False,
                indent=2,
            )

    with StageTimer("load cell obs from module h5ad"):
        adata_meta = ad.read_h5ad(cfg.module_h5ad_path, backed="r")
        obs = adata_meta.obs.copy()
        adata_meta.file.close()
        obs = merge_obs_sidecar(obs, load_obs_sidecar(cfg), merge=cfg.merge_obs_sidecar)
        _log(f"n_cells = {len(obs):,}")
        for k in [cfg.donor_key, cfg.disease_key, cfg.celltype_key]:
            if k not in obs.columns:
                raise KeyError(f"'{k}' not found in module_h5ad.obs")

    with StageTimer("build zarr row-position mapping"):
        zarr_row_positions = build_module_obs_to_zarr_positions(cfg.zarr_path, obs)

    with StageTimer("confounding audit"):
        report = diagnose_batch_confounding(
            obs,
            donor_key=cfg.donor_key,
            disease_key=cfg.disease_key,
            candidate_batch_keys=cfg.batch_candidate_keys,
            mixing_threshold=cfg.tech_mixing_threshold,
            min_median_donors_per_level=cfg.tech_min_median_donors_per_level,
            max_levels_fraction_of_donors=cfg.tech_max_levels_fraction_of_donors,
            donor_key_forbid=cfg.donor_key,
        )
        print_diagnosis(report)
        with open(out_dir / "confounding_diagnosis.json", "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    with StageTimer("donor-level split"):
        if cfg.use_existing_split:
            if cfg.existing_split_key not in obs.columns:
                raise KeyError(
                    f"use_existing_split=true but {cfg.existing_split_key!r} not found in obs. "
                    "Provide an obs sidecar with a split column."
                )
            obs_with_split = obs.copy()
            if cfg.existing_split_key != "split":
                obs_with_split["split"] = obs_with_split[cfg.existing_split_key].astype(str).values
            else:
                obs_with_split["split"] = obs_with_split["split"].astype(str).values
            _log(f"using existing split column: {cfg.existing_split_key}")
        else:
            adata_obs = ad.AnnData(obs=obs.copy())
            adata_obs = call_donor_split(donor_split_fn, adata_obs, cfg)
            obs_with_split = adata_obs.obs.copy()
        split_counts = {str(k): int(v) for k, v in obs_with_split["split"].value_counts().items()}
        _log(f"split cell counts: {split_counts}")

    with StageTimer("attach subject_id / tech_id_model"):
        obs_with_split, tech_key_chosen = choose_tech_covariate(
            obs_with_split,
            report,
            donor_key=cfg.donor_key,
            preferred_tech_covariate=preferred_tech_covariate,
            auto_pick_priority=cfg.auto_pick_priority,
            auto_pick_blocklist=cfg.auto_pick_blocklist,
            fail_if_no_safe_tech=cfg.fail_if_no_safe_tech,
            allow_unsafe_preferred=allow_unsafe_preferred_tech,
        )
        _log(f"subject_id    <- {cfg.donor_key}")
        _log(f"tech_id_model <- {tech_key_chosen if tech_key_chosen else '__GLOBAL__'}")

        keep_cols = [
            c
            for c in [
                cfg.donor_key,
                cfg.disease_key,
                cfg.celltype_key,
                "split",
                "subject_id",
                "subject_int_id",
                "tech_id_model",
            ]
            if c in obs_with_split.columns
        ]
        obs_with_split[keep_cols].to_csv(out_dir / "residual_cells_split.csv", index=True)

        donor_summary = (
            obs_with_split.groupby(cfg.donor_key, observed=True)
            .agg(
                split=("split", _unique_split_or_multi),
                disease_set=(cfg.disease_key, _sorted_non_null_str_tuple),
                n_cells=(cfg.donor_key, "size"),
                tech_levels=("tech_id_model", _sorted_str_tuple),
            )
            .reset_index()
        )
        donor_summary.to_csv(out_dir / "donor_split_summary.csv", index=False)

    splits_to_emit = ["train"]
    if not cfg.skip_val_test_extraction:
        splits_to_emit += ["val", "test"]

    for split_name in splits_to_emit:
        mask = obs_with_split["split"].astype(str).values == split_name
        if not np.any(mask):
            _log(f"skip {split_name}: no cells")
            continue
        row_positions = zarr_row_positions[mask].astype(np.int64)

        with StageTimer(
            f"extract residual counts for split='{split_name}' "
            f"({int(mask.sum()):,} cells × {len(residual_names):,} genes)"
        ):
            X = extract_residual_counts_for_cells(
                zarr_path=cfg.zarr_path,
                residual_gene_names=residual_names,
                row_positions=row_positions,
                extract_csr_subset_from_zarr=extract_csr_subset_from_zarr,
                max_gap=cfg.max_gap,
                chunk_rows_fallback=cfg.chunk_rows_fallback,
            )

        split_obs = obs_with_split.loc[mask].copy()
        split_var = pd.DataFrame(index=pd.Index(residual_names, name="ensembl_id"))
        adata_out = ad.AnnData(X=X, obs=split_obs, var=split_var)

        write_counts_layer_this_split = resolve_write_counts_layer_for_split(cfg, split_name)
        if write_counts_layer_this_split:
            adata_out.layers["counts"] = X.copy()

        adata_out.uns["description"] = (
            f"Residual-gene subset of SEA-AD (split='{split_name}'). "
            "Residual genes are a feature universe, not a cell filter. "
            "Cell selection is driven purely by the donor-level split."
        )
        adata_out.uns["registry_path"] = str(cfg.registry_json_path)
        adata_out.uns["split"] = split_name
        adata_out.uns["subject_key"] = "subject_id"
        adata_out.uns["tech_key_model"] = "tech_id_model"
        adata_out.uns["tech_key_source"] = tech_key_chosen if tech_key_chosen else "__GLOBAL__"
        adata_out.uns["X_semantics"] = "raw_counts"
        adata_out.uns["counts_stored_in_X"] = True
        adata_out.uns["counts_layer_written"] = write_counts_layer_this_split
        if write_counts_layer_this_split:
            adata_out.uns["counts_layer_semantics"] = "raw_counts"
        adata_out.uns["hdwgcna_input_eligible"] = split_name == "train"

        out_path = out_dir / f"residual_{split_name}.h5ad"
        adata_out.write_h5ad(str(out_path))
        _log(f"saved {split_name}: {out_path} (shape={adata_out.shape})")

    counts_layer_policy = {
        "train": resolve_write_counts_layer_for_split(cfg, "train"),
        "val": resolve_write_counts_layer_for_split(cfg, "val"),
        "test": resolve_write_counts_layer_for_split(cfg, "test"),
    }

    manifest = {
        "config": asdict(cfg),
        "registry": {
            "n_total_registry_genes": registry_n,
            "n_residual_genes": len(residual_names),
        },
        "cells": {
            "n_total": int(len(obs_with_split)),
            "split_counts": split_counts,
        },
        "audit": {
            "donor_disease_purity": report["donor_disease_purity"],
            "donor_perfectly_confounded_with_disease": report["donor_perfectly_confounded_with_disease"],
            "recommended_tech_covariates": report["recommended_tech_covariates"],
            "tech_key_chosen_for_tech_id_model": tech_key_chosen if tech_key_chosen else "__GLOBAL__",
        },
        "schema_emitted": {
            "subject_id": f"{cfg.donor_key} as string (biological / subject)",
            "subject_int_id": "int64 code for subject_id",
            "tech_id_model": "chosen safe tech covariate, or '__GLOBAL__'",
            "split": "train / val / test",
        },
        "counts_layer_policy": counts_layer_policy,
        "outputs": {
            "residual_gene_list_json": str(out_dir / "residual_gene_list.json"),
            "confounding_diagnosis_json": str(out_dir / "confounding_diagnosis.json"),
            "residual_cells_split_csv": str(out_dir / "residual_cells_split.csv"),
            "donor_split_summary_csv": str(out_dir / "donor_split_summary.csv"),
            "residual_train_h5ad": str(out_dir / "residual_train.h5ad"),
            "residual_val_h5ad": None if cfg.skip_val_test_extraction else str(out_dir / "residual_val.h5ad"),
            "residual_test_h5ad": None if cfg.skip_val_test_extraction else str(out_dir / "residual_test.h5ad"),
            "report_txt": str(out_dir / "report.txt"),
        },
    }
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=str)

    write_text_report(
        out_dir / "report.txt",
        cfg=cfg,
        report=report,
        residual_n=len(residual_names),
        registry_n=registry_n,
        n_cells_total=int(len(obs_with_split)),
        split_counts=split_counts,
        tech_key_chosen=tech_key_chosen,
        counts_layer_policy=counts_layer_policy,
    )
    _log("ALL DONE")


# =====================================================================
# CLI
# =====================================================================
def load_cfg(path: str | Path) -> ResidualPrepConfig:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    return ResidualPrepConfig(**d)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Prepare residual-gene subsets (train/val/test) for hdWGCNA. "
            "Also emits subject_id / tech_id_model schema and a donor × disease confounding audit."
        )
    )
    p.add_argument("--config", type=str, required=True, help="Path to JSON with ResidualPrepConfig fields.")
    p.add_argument("--preferred-tech", type=str, default=None, help="Safely prefer a specific tech covariate column.")
    p.add_argument(
        "--allow-unsafe-preferred-tech",
        action="store_true",
        help="Allow --preferred-tech even if it fails safety screening (not recommended).",
    )
    return p


def main() -> None:
    args = build_argparser().parse_args()
    cfg = load_cfg(args.config)
    if args.preferred_tech == cfg.donor_key:
        raise ValueError(
            f"--preferred-tech cannot equal donor_key ('{cfg.donor_key}'). donor_id must never be used as tech_id_model."
        )
    prepare_residual_for_hdwgcna(
        cfg,
        preferred_tech_covariate=args.preferred_tech,
        allow_unsafe_preferred_tech=args.allow_unsafe_preferred_tech,
    )


if __name__ == "__main__":
    main()
