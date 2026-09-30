from __future__ import annotations

"""
build_celltype_control_scaler_zarr_fast.py
==========================================

Fast one-pass builder for celltype × control-train gene scaler.

Run after:
  1. pathology-aware donor split
       -> obs_with_pathology_split.parquet
  2. make_ordinal_data_zarr.py
       -> final_gene_union_ordinal_spec.npz

It builds:

    x_gene_scalar[i, g]
      = (log1p(count[i, g]) - mean_control_train[celltype_i, g])
        / (std_control_train[celltype_i, g] + eps)

for the final module-registry gene union.

Why this version is faster
--------------------------
Older versions computed:

    global moments once
    + celltype moments once per celltype

which re-read the same Zarr CSR matrix many times.

This version reads train-control rows once and accumulates:

    global moments
    all celltype-specific moments

in the same pass.
"""

import argparse
import json
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import zarr


@dataclass
class ControlFilter:
    column: str
    values: list[Any]


@dataclass
class CelltypeControlScalerFastConfig:
    zarr_path: str
    spec_path: str
    gene_indices_path: str
    out_npz_path: str
    out_json_path: Optional[str] = None

    obs_parquet_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_h5ad_path: Optional[str] = None

    matrix_path: str = "raw/X"
    counts_layer: str = "counts"
    prefer_raw: bool = True

    split_key: str = "split"
    train_split_name: str = "train"
    celltype_key: str = "Subclass"
    donor_key: str = "donor_id"

    control_filters: list[dict[str, Any]] = field(
        default_factory=lambda: [
            {"column": "disease", "values": ["normal"]},
            {"column": "ADNC", "values": ["Reference", "Not AD", "Low"]},
            {"column": "LATE-NC stage", "values": ["Reference", "Not Identified"]},
            {"column": "Lewy body disease pathology", "values": ["Reference", "Not Identified"]},
            {"column": "Microinfarct pathology", "values": ["Reference", "0 to 3 microinfarcts"]},
        ]
    )
    case_insensitive_filter: bool = True

    donor_balanced: bool = True
    chunk_size: int = 2048
    std_floor: float = 1e-3

    max_csr_data_span: int = 50_000_000
    max_span_to_selected_nnz_ratio: float = 20.0

    use_shrinkage: bool = True
    shrinkage_k_donors: float = 5.0
    shrinkage_k_cells: Optional[float] = None

    warn_min_control_cells: int = 30
    warn_min_control_donors: int = 3

    check_counts_sanity: bool = True
    counts_sanity_n_cells: int = 64
    counts_integer_atol: float = 1e-3
    counts_integer_fraction_threshold: float = 0.95

    scaler_version: str = "celltype_control_scaler_zarr_fast_v1_one_pass"
    description: str = (
        "fast one-pass celltype x low-pathology control-train donor-balanced "
        "log1p scaler for final module-registry gene union"
    )


def _log(msg: str) -> None:
    print(f"[celltype-control-scaler-fast] {msg}", flush=True)


def load_config(path: str | Path) -> CelltypeControlScalerFastConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if "control_filters" not in raw and "control_key" in raw:
        raw["control_filters"] = [
            {"column": raw.pop("control_key"), "values": raw.pop("control_values")}
        ]
    return CelltypeControlScalerFastConfig(**raw)


@dataclass
class MinimalOrdinalSpec:
    gene_names: np.ndarray
    celltype_vocab: np.ndarray

    @classmethod
    def load(cls, npz_path: str | Path) -> "MinimalOrdinalSpec":
        obj = np.load(npz_path, allow_pickle=True)
        for key in ["gene_names", "celltype_vocab"]:
            if key not in obj:
                raise KeyError(f"ordinal spec must contain {key!r}.")
        return cls(
            gene_names=np.asarray(obj["gene_names"], dtype=object),
            celltype_vocab=np.asarray(obj["celltype_vocab"], dtype=object),
        )


def _read_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError(
            "Could not read AnnData-Zarr metadata with anndata.experimental.read_elem."
        ) from e


def read_zarr_dataframe(root, key_path: str) -> pd.DataFrame:
    group = root
    for part in key_path.split("/"):
        if part:
            group = group[part]
    out = _read_elem(group)
    if not isinstance(out, pd.DataFrame):
        raise TypeError(f"Zarr element {key_path!r} is not a pandas DataFrame.")
    return out.copy()


def _zarr_get(root, path: str):
    obj = root
    for part in path.split("/"):
        if part:
            obj = obj[part]
    return obj


def choose_zarr_matrix(root, cfg: CelltypeControlScalerFastConfig):
    requested = str(cfg.matrix_path).strip()
    if requested and requested != "auto":
        matrix = _zarr_get(root, requested)
        var_key = "raw/var" if requested.startswith("raw/") else "var"
        return matrix, requested, var_key

    layer_path = f"layers/{cfg.counts_layer}"
    if "layers" in root and cfg.counts_layer in root["layers"]:
        return root["layers"][cfg.counts_layer], layer_path, "var"
    if cfg.prefer_raw and "raw" in root and "X" in root["raw"]:
        return root["raw"]["X"], "raw/X", "raw/var"
    if "X" in root:
        return root["X"], "X", "var"
    raise KeyError("Could not find matrix path: no layers/counts, raw/X, or X.")


def load_obs_source(cfg: CelltypeControlScalerFastConfig, zarr_obs: pd.DataFrame) -> pd.DataFrame:
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
        import anndata as ad
        _log(f"loading obs sidecar h5ad: {cfg.obs_h5ad_path}")
        return ad.read_h5ad(cfg.obs_h5ad_path, backed="r").obs.copy()

    _log("using zarr obs directly")
    return zarr_obs.copy()


def align_obs_to_zarr(obs_source: pd.DataFrame, zarr_obs: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    obs_source = obs_source.copy()
    obs_source.index = obs_source.index.astype(str)

    zarr_index = pd.Index(zarr_obs.index.astype(str))
    source_index = pd.Index(obs_source.index.astype(str))

    if len(source_index) == len(zarr_index) and np.array_equal(source_index.to_numpy(), zarr_index.to_numpy()):
        return obs_source, np.arange(len(zarr_index), dtype=np.int64)

    if source_index.has_duplicates:
        dup = source_index[source_index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"obs source has duplicated cell barcodes. Examples: {dup}")

    pos = zarr_index.get_indexer(source_index)
    missing = pos < 0
    if missing.any():
        examples = source_index[missing].tolist()[:5]
        raise ValueError(f"Some obs sidecar barcodes were not found in zarr obs. Examples: {examples}")

    return obs_source, pos.astype(np.int64)


def _normalise_allowed(values: list[Any], *, case_insensitive: bool) -> set[str]:
    out = {str(v) for v in values}
    return {v.lower() for v in out} if case_insensitive else out


def _series_as_str(series: pd.Series, *, case_insensitive: bool) -> pd.Series:
    s = series.astype(str).fillna("__NA__")
    return s.str.lower() if case_insensitive else s


def build_control_mask(obs: pd.DataFrame, cfg: CelltypeControlScalerFastConfig) -> np.ndarray:
    if not cfg.control_filters:
        raise ValueError("At least one control filter must be provided.")

    mask = np.ones(len(obs), dtype=bool)
    for filt_raw in cfg.control_filters:
        filt = ControlFilter(**filt_raw)
        if filt.column not in obs.columns:
            raise KeyError(f"control filter column {filt.column!r} not found in obs source.")
        allowed = _normalise_allowed(list(filt.values), case_insensitive=cfg.case_insensitive_filter)
        vals = _series_as_str(obs[filt.column], case_insensitive=cfg.case_insensitive_filter)
        mask &= vals.isin(allowed).to_numpy()
    return mask


class ZarrSelectedGeneReader:
    def __init__(self, matrix, *, n_cols_full: int, gene_indices: np.ndarray, name: str) -> None:
        self.matrix = matrix
        self.n_cols_full = int(n_cols_full)
        self.gene_indices = np.asarray(gene_indices, dtype=np.int64)
        self.name = str(name)

        if self.gene_indices.ndim != 1:
            raise ValueError("gene_indices must be 1D.")
        if len(self.gene_indices) == 0:
            raise ValueError("gene_indices is empty.")
        if self.gene_indices.min() < 0 or self.gene_indices.max() >= self.n_cols_full:
            raise ValueError(
                f"gene_indices out of bounds for n_cols_full={self.n_cols_full}: "
                f"min={self.gene_indices.min()}, max={self.gene_indices.max()}."
            )
        if len(np.unique(self.gene_indices)) != len(self.gene_indices):
            raise ValueError("gene_indices contains duplicates.")

        self.n_genes = int(len(self.gene_indices))
        self.is_csr = hasattr(matrix, "keys") and all(k in matrix for k in ("data", "indices", "indptr"))

        if self.is_csr:
            self.data = matrix["data"]
            self.indices = matrix["indices"]
            self.indptr = matrix["indptr"][:]  # load once; cheap and much faster
            self.shape = (len(self.indptr) - 1, self.n_cols_full)
            self.full_to_selected = np.full(self.n_cols_full, -1, dtype=np.int32)
            self.full_to_selected[self.gene_indices] = np.arange(self.n_genes, dtype=np.int32)
        else:
            self.shape = tuple(matrix.shape)
            if len(self.shape) != 2:
                raise ValueError(f"Matrix {name!r} must be 2D or CSR group, got shape={self.shape}.")
            if self.shape[1] != self.n_cols_full:
                raise ValueError(
                    f"Matrix {name!r} gene dimension mismatch: shape[1]={self.shape[1]}, "
                    f"n_cols_full={self.n_cols_full}."
                )

    def rows_to_dense(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        if len(rows) == 0:
            return np.zeros((0, self.n_genes), dtype=np.float32)
        if rows.min() < 0 or rows.max() >= self.shape[0]:
            raise IndexError("row index out of bounds for zarr matrix.")

        if self.is_csr:
            out = np.zeros((len(rows), self.n_genes), dtype=np.float32)
            for i, r in enumerate(rows):
                s = int(self.indptr[int(r)])
                e = int(self.indptr[int(r) + 1])
                if e <= s:
                    continue
                idx = np.asarray(self.indices[s:e], dtype=np.int64)
                dat = np.asarray(self.data[s:e], dtype=np.float32)
                mapped = self.full_to_selected[idx]
                keep = mapped >= 0
                if np.any(keep):
                    out[i, mapped[keep]] = dat[keep]
            return out

        try:
            block = self.matrix.get_orthogonal_selection((rows, self.gene_indices))
        except Exception:
            try:
                block = self.matrix[np.ix_(rows, self.gene_indices)]
            except Exception:
                block = np.vstack([np.asarray(self.matrix[int(r), :])[self.gene_indices] for r in rows])
        return np.asarray(block, dtype=np.float32)


def warn_if_not_raw_counts(reader: ZarrSelectedGeneReader, candidate_rows: np.ndarray, *, cfg: CelltypeControlScalerFastConfig) -> None:
    if not cfg.check_counts_sanity or cfg.counts_sanity_n_cells <= 0 or len(candidate_rows) == 0:
        return

    rows = np.asarray(candidate_rows, dtype=np.int64)
    take_n = min(cfg.counts_sanity_n_cells, len(rows))
    take_pos = np.linspace(0, len(rows) - 1, num=take_n, dtype=np.int64)
    sample_rows = np.unique(rows[take_pos])

    block = reader.rows_to_dense(sample_rows).ravel()
    finite = np.isfinite(block)
    if not finite.all():
        warnings.warn("Non-finite values found in count sanity sample.", RuntimeWarning)
        block = block[finite]
    if block.size == 0:
        return
    if np.any(block < 0):
        warnings.warn("Negative values found; raw non-negative counts expected.", RuntimeWarning)

    nonzero = block[np.abs(block) > 0]
    vals = nonzero if nonzero.size > 0 else block
    frac_int_like = float(np.mean(np.isclose(vals, np.round(vals), atol=cfg.counts_integer_atol, rtol=0.0)))
    if frac_int_like < cfg.counts_integer_fraction_threshold:
        warnings.warn(
            f"Matrix {reader.name} may not be raw counts: integer-like fraction={frac_int_like:.3f}.",
            RuntimeWarning,
        )


def donor_balanced_weights(obs_sub: pd.DataFrame, donor_key: str) -> tuple[np.ndarray, int]:
    if donor_key not in obs_sub.columns:
        raise KeyError(f"donor_key {donor_key!r} not found in selected obs.")

    donors = obs_sub[donor_key].astype(str).fillna("__NA__").to_numpy()
    unique, inverse, counts = np.unique(donors, return_inverse=True, return_counts=True)
    n_donors = len(unique)
    if n_donors == 0:
        return np.zeros(0, dtype=np.float64), 0

    weights = 1.0 / (float(n_donors) * counts[inverse].astype(np.float64))
    return weights.astype(np.float64), int(n_donors)


def uniform_weights(n: int) -> tuple[np.ndarray, int]:
    if n <= 0:
        return np.zeros(0, dtype=np.float64), 0
    return np.full(n, 1.0 / float(n), dtype=np.float64), int(n)


def count_unique_donors(obs_sub: pd.DataFrame, donor_key: str) -> int:
    if donor_key not in obs_sub.columns:
        return 0
    return int(obs_sub[donor_key].astype(str).fillna("__NA__").nunique())


def _accumulate_row(
    *,
    sum_x: np.ndarray,
    sum_x2: np.ndarray,
    out_row: int,
    mapped: np.ndarray,
    data: np.ndarray,
    weight: float,
) -> None:
    if weight <= 0.0 or mapped.size == 0:
        return
    x = np.log1p(data.astype(np.float32, copy=False)).astype(np.float64, copy=False)
    m = mapped.astype(np.int64, copy=False)
    sum_x[out_row, m] += weight * x
    sum_x2[out_row, m] += weight * x * x


def compute_all_celltype_moments_one_pass(
    *,
    reader: ZarrSelectedGeneReader,
    rows: np.ndarray,
    global_weights: np.ndarray,
    celltype_row_ids: np.ndarray,
    celltype_weights: np.ndarray,
    n_celltypes: int,
    chunk_size: int,
    std_floor: float,
    max_csr_data_span: int,
    max_span_to_selected_nnz_ratio: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute all celltype moments plus global moments in one matrix pass."""
    rows = np.asarray(rows, dtype=np.int64)
    global_weights = np.asarray(global_weights, dtype=np.float64)
    celltype_row_ids = np.asarray(celltype_row_ids, dtype=np.int64)
    celltype_weights = np.asarray(celltype_weights, dtype=np.float64)

    if len(rows) == 0:
        raise ValueError("Cannot compute moments for an empty row set.")
    if not (len(rows) == len(global_weights) == len(celltype_row_ids) == len(celltype_weights)):
        raise ValueError("rows / weights / celltype ids length mismatch.")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive.")

    wg_sum = float(global_weights.sum())
    if not np.isfinite(wg_sum) or wg_sum <= 0:
        raise ValueError("global_weights must have positive finite sum.")
    global_weights = global_weights / wg_sum

    n_rows_out = n_celltypes + 1
    global_row = n_celltypes
    n_genes = reader.n_genes

    sum_x = np.zeros((n_rows_out, n_genes), dtype=np.float64)
    sum_x2 = np.zeros((n_rows_out, n_genes), dtype=np.float64)

    order = np.argsort(rows)
    rows_s = rows[order]
    wg_s = global_weights[order]
    ct_s = celltype_row_ids[order]
    wc_s = celltype_weights[order]

    if not reader.is_csr:
        for start in range(0, len(rows_s), chunk_size):
         # rows_s: zarr row _index in train-control cell
            end = min(start + chunk_size, len(rows_s))
            rr = rows_s[start:end] # zarr row indices for this chunk
            block = reader.rows_to_dense(rr)
            np.log1p(block, out=block)

            for j in range(len(rr)):
                wg = wg_s[start + j]
                if wg > 0:
                    sum_x[global_row] += wg * block[j]
                    sum_x2[global_row] += wg * block[j] * block[j]

                ct = int(ct_s[start + j])
                wc = float(wc_s[start + j])
                if 0 <= ct < n_celltypes and wc > 0:
                    sum_x[ct] += wc * block[j]
                    sum_x2[ct] += wc * block[j] * block[j]

            if end == len(rows_s) or (end % (chunk_size * 20) == 0):
                _log(f"  one-pass dense moment rows processed: {end:,}/{len(rows_s):,}")

    else:
        indptr = reader.indptr

        for start in range(0, len(rows_s), chunk_size):
        # rows_s: zarr row _index in train-control cell
            end = min(start + chunk_size, len(rows_s))
            rr = rows_s[start:end] # zarr row indices for this chunk
	    
	     # indptr[r], indptr[r+1] tells us starting and end points
    	     # of nonzero data of row r in CSR. 
            selected_nnz = int(np.sum(indptr[rr + 1] - indptr[rr])) # the number of non-zeros
            
            # p0~p1 is continuos interval that cover from 1st and last row
            p0 = int(indptr[int(rr[0])])
            p1 = int(indptr[int(rr[-1]) + 1])
            span = max(0, p1 - p0)
	    
	    # use_span=True then read all at once in this chunk rather than not reading every data in every row 
            use_span = (
                span > 0
                and span <= int(max_csr_data_span)
                and span <= float(max_span_to_selected_nnz_ratio) * max(selected_nnz, 1)
            )

            if use_span:
                idx_span = np.asarray(reader.indices[p0:p1], dtype=np.int64)
                dat_span = np.asarray(reader.data[p0:p1], dtype=np.float32)

                for j, r in enumerate(rr):
                    # convert start and end position of row r to local index of p0 	 
                    s = int(indptr[int(r)]) - p0
                    e = int(indptr[int(r) + 1]) - p0
                    if e <= s:
                        continue

                    idx = idx_span[s:e]
                    dat = dat_span[s:e]

                    mapped = reader.full_to_selected[idx]
                    keep = mapped >= 0
                    if not np.any(keep):
                        continue
		    # m:gene-index in final gene-union
		    # xdat: raw count in that gene   	
                    m = mapped[keep]
                    xdat = dat[keep]

                    _accumulate_row(
                        sum_x=sum_x,
                        sum_x2=sum_x2,
                        out_row=global_row,
                        mapped=m,
                        data=xdat,
                        weight=float(wg_s[start + j]),
                    )
                    
                    # celltype-specific moment
                    # ct: cell type row id of the cell

                    ct = int(ct_s[start + j])
                    if 0 <= ct < n_celltypes:
                        _accumulate_row(
                            sum_x=sum_x,
                            sum_x2=sum_x2,
                            out_row=ct,
                            mapped=m,
                            data=xdat,
                            weight=float(wc_s[start + j]),
                        )

            else:
                for j, r in enumerate(rr):
                    s = int(indptr[int(r)])
                    e = int(indptr[int(r) + 1])
                    if e <= s:
                        continue

                    idx = np.asarray(reader.indices[s:e], dtype=np.int64)
                    dat = np.asarray(reader.data[s:e], dtype=np.float32)

                    mapped = reader.full_to_selected[idx]
                    keep = mapped >= 0
                    if not np.any(keep):
                        continue

                    m = mapped[keep]
                    xdat = dat[keep]

                    _accumulate_row(
                        sum_x=sum_x,
                        sum_x2=sum_x2,
                        out_row=global_row,
                        mapped=m,
                        data=xdat,
                        weight=float(wg_s[start + j]),
                    )

                    ct = int(ct_s[start + j])
                    if 0 <= ct < n_celltypes:
                        _accumulate_row(
                            sum_x=sum_x,
                            sum_x2=sum_x2,
                            out_row=ct,
                            mapped=m,
                            data=xdat,
                            weight=float(wc_s[start + j]),
                        )

            if end == len(rows_s) or (end % (chunk_size * 20) == 0):
                _log(f"  one-pass CSR moment rows processed: {end:,}/{len(rows_s):,}")

    var = sum_x2 - np.square(sum_x)
    var = np.maximum(var, 0.0)
    var = np.maximum(var, float(std_floor) ** 2)
    std = np.sqrt(var)

    return sum_x.astype(np.float32), std.astype(np.float32)


def shrink_moments(
    mean_ct: np.ndarray,
    std_ct: np.ndarray,
    mean_global: np.ndarray,
    std_global: np.ndarray,
    *,
    n_effective: int,
    cfg: CelltypeControlScalerFastConfig,
) -> tuple[np.ndarray, np.ndarray, float]:
    if not cfg.use_shrinkage:
        return mean_ct.astype(np.float32), std_ct.astype(np.float32), 0.0

    k = cfg.shrinkage_k_donors if cfg.donor_balanced or cfg.shrinkage_k_cells is None else cfg.shrinkage_k_cells
    k = float(k)
    if k < 0:
        raise ValueError("shrinkage k must be non-negative.")

    lam = k / (k + max(float(n_effective), 0.0))
    lam = min(max(lam, 0.0), 1.0)

    mean = (1.0 - lam) * mean_ct + lam * mean_global
    var = (1.0 - lam) * np.square(std_ct) + lam * np.square(std_global)
    std = np.sqrt(np.maximum(var, float(cfg.std_floor) ** 2))

    return mean.astype(np.float32), std.astype(np.float32), float(lam)


def build_celltype_control_scaler(cfg: CelltypeControlScalerFastConfig) -> dict[str, Any]:
    out_npz = Path(cfg.out_npz_path)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    if cfg.out_json_path is not None:
        Path(cfg.out_json_path).parent.mkdir(parents=True, exist_ok=True)

    if cfg.chunk_size <= 0:
        raise ValueError("chunk_size must be positive.")
    if cfg.std_floor <= 0:
        raise ValueError("std_floor must be positive.")

    _log(f"loading ordinal spec: {cfg.spec_path}")
    spec = MinimalOrdinalSpec.load(cfg.spec_path)

    _log(f"loading gene indices: {cfg.gene_indices_path}")
    gene_indices = np.load(cfg.gene_indices_path).astype(np.int64)

    _log(f"opening zarr: {cfg.zarr_path}")
    root = zarr.open(cfg.zarr_path, mode="r")

    _log("loading zarr obs/var metadata")
    zarr_obs = read_zarr_dataframe(root, "obs")
    matrix, matrix_source, var_key = choose_zarr_matrix(root, cfg)
    var = read_zarr_dataframe(root, var_key)
    var_names_full = np.asarray(var.index.astype(str), dtype=object)

    if len(gene_indices) != len(spec.gene_names):
        raise ValueError(
            f"gene_indices length ({len(gene_indices)}) != spec gene_names length ({len(spec.gene_names)})."
        )
    if gene_indices.min() < 0 or gene_indices.max() >= len(var_names_full):
        raise ValueError("gene_indices out of bounds for zarr var_names.")

    gene_names_from_zarr = np.asarray(var_names_full[gene_indices], dtype=object)
    if not np.array_equal(gene_names_from_zarr, np.asarray(spec.gene_names, dtype=object)):
        raise ValueError(
            "Gene order mismatch: zarr var_names[gene_indices] != spec.gene_names. "
            "Check gene_union_original_indices.npy and matrix_path/raw var."
        )

    obs_source = load_obs_source(cfg, zarr_obs)
    obs, zarr_row_positions = align_obs_to_zarr(obs_source, zarr_obs)

    for key in [cfg.split_key, cfg.celltype_key, cfg.donor_key]:
        if key not in obs.columns:
            raise KeyError(f"{key!r} not found in obs source.")

    reader = ZarrSelectedGeneReader(
        matrix,
        n_cols_full=len(var_names_full),
        gene_indices=gene_indices,
        name=matrix_source,
    )

    train_mask_local = obs[cfg.split_key].astype(str).to_numpy() == str(cfg.train_split_name)
    control_mask_local = build_control_mask(obs, cfg)
    selected_mask_local = train_mask_local & control_mask_local

    if selected_mask_local.sum() == 0:
        raise ValueError("No train-control cells found. Check split/control filters.")

    control_rows_local = np.where(selected_mask_local)[0].astype(np.int64)
    control_rows_zarr = zarr_row_positions[control_rows_local].astype(np.int64)
    obs_control = obs.iloc[control_rows_local].copy()

    _log(
        f"matrix={matrix_source}; final genes={len(gene_indices):,}; "
        f"obs cells={len(obs):,}; train cells={int(train_mask_local.sum()):,}; "
        f"train-control cells={len(control_rows_zarr):,}"
    )

    warn_if_not_raw_counts(reader, control_rows_zarr, cfg=cfg)

    spec_celltypes = [str(x) for x in spec.celltype_vocab.tolist()]
    global_name = "__GLOBAL__"
    output_celltypes = list(spec_celltypes)
    if global_name not in output_celltypes:
        output_celltypes.append(global_name)

    n_celltypes = len(spec_celltypes)
    global_row = n_celltypes
    n_rows = n_celltypes + 1
    n_genes = len(spec.gene_names)

    global_actual_donors = count_unique_donors(obs_control, cfg.donor_key)
    if cfg.donor_balanced:
        global_weights, global_effective_n = donor_balanced_weights(obs_control, cfg.donor_key)
    else:
        global_weights, global_effective_n = uniform_weights(len(obs_control))

    celltype_to_row = {ct: i for i, ct in enumerate(spec_celltypes)}
    celltype_values_control = (
        obs_control[cfg.celltype_key]
        .astype(str)
        .fillna("__NA__")
        .to_numpy()
    )
    celltype_row_ids = np.asarray(
        [celltype_to_row.get(ct, -1) for ct in celltype_values_control],
        dtype=np.int64,
    )

    unknown_ct = int(np.sum(celltype_row_ids < 0))
    if unknown_ct > 0:
        warnings.warn(
            f"{unknown_ct:,} train-control cells have celltype not present in spec.celltype_vocab. "
            "They contribute only to __GLOBAL__ moments.",
            RuntimeWarning,
        )

    celltype_weights = np.zeros(len(obs_control), dtype=np.float64)

    n_control_cells = np.zeros(n_rows, dtype=np.int64)
    n_control_donors = np.zeros(n_rows, dtype=np.int64)
    shrinkage_effective_n = np.zeros(n_rows, dtype=np.int64)

    for ct_id, ct in enumerate(spec_celltypes):
        mask = celltype_row_ids == ct_id
        n_control_cells[ct_id] = int(mask.sum())

        if not np.any(mask):
            continue

        obs_ct = obs_control.iloc[np.where(mask)[0]]
        n_control_donors[ct_id] = count_unique_donors(obs_ct, cfg.donor_key)

        if cfg.donor_balanced:
            w_ct, eff_ct = donor_balanced_weights(obs_ct, cfg.donor_key)
        else:
            w_ct, eff_ct = uniform_weights(int(mask.sum()))

        celltype_weights[mask] = w_ct
        shrinkage_effective_n[ct_id] = int(eff_ct)

    n_control_cells[global_row] = len(control_rows_zarr)
    n_control_donors[global_row] = global_actual_donors
    shrinkage_effective_n[global_row] = int(global_effective_n)

    _log(
        f"computing moments in one pass: cells={len(control_rows_zarr):,}, "
        f"global_donors={global_actual_donors}, global_effective_n={global_effective_n}"
    )

    mean_raw, std_raw = compute_all_celltype_moments_one_pass(
        reader=reader,
        rows=control_rows_zarr,
        global_weights=global_weights,
        celltype_row_ids=celltype_row_ids,
        celltype_weights=celltype_weights,
        n_celltypes=n_celltypes,
        chunk_size=cfg.chunk_size,
        std_floor=cfg.std_floor,
        max_csr_data_span=cfg.max_csr_data_span,
        max_span_to_selected_nnz_ratio=cfg.max_span_to_selected_nnz_ratio,
    )

    mean_global = mean_raw[global_row]
    std_global = std_raw[global_row]

    mean = np.zeros((n_rows, n_genes), dtype=np.float32)
    std = np.zeros((n_rows, n_genes), dtype=np.float32)
    shrinkage_lambda = np.zeros(n_rows, dtype=np.float32)
    used_global_fallback = np.zeros(n_rows, dtype=bool)

    for row_id, ct in enumerate(output_celltypes):
        if ct == global_name:
            mean[row_id] = mean_global
            std[row_id] = std_global
            continue

        if n_control_cells[row_id] == 0:
            warnings.warn(
                f"No train-control cells for celltype {ct!r}. Using __GLOBAL__ scaler row.",
                RuntimeWarning,
            )
            mean[row_id] = mean_global
            std[row_id] = std_global
            shrinkage_effective_n[row_id] = 0
            shrinkage_lambda[row_id] = 1.0
            used_global_fallback[row_id] = True
            continue

        if n_control_cells[row_id] < cfg.warn_min_control_cells or n_control_donors[row_id] < cfg.warn_min_control_donors:
            warnings.warn(
                f"Celltype {ct!r} has only {int(n_control_cells[row_id])} train-control cells and "
                f"{int(n_control_donors[row_id])} train-control donors. Shrinkage may dominate.",
                RuntimeWarning,
            )

        mean[row_id], std[row_id], lam = shrink_moments(
            mean_raw[row_id],
            std_raw[row_id],
            mean_global,
            std_global,
            n_effective=int(shrinkage_effective_n[row_id]),
            cfg=cfg,
        )
        shrinkage_lambda[row_id] = lam

    metadata = {
        "scaler_version": cfg.scaler_version,
        "description": cfg.description,
        "build_date": datetime.now(timezone.utc).isoformat(),
        "zarr_path": str(cfg.zarr_path),
        "spec_path": str(cfg.spec_path),
        "gene_indices_path": str(cfg.gene_indices_path),
        "obs_source": (
            str(cfg.obs_parquet_path)
            if cfg.obs_parquet_path is not None
            else str(cfg.obs_csv_path)
            if cfg.obs_csv_path is not None
            else str(cfg.obs_h5ad_path)
            if cfg.obs_h5ad_path is not None
            else "zarr_obs"
        ),
        "matrix_source": matrix_source,
        "matrix_path": cfg.matrix_path,
        "split_key": cfg.split_key,
        "train_split_name": cfg.train_split_name,
        "celltype_key": cfg.celltype_key,
        "donor_key": cfg.donor_key,
        "donor_balanced": bool(cfg.donor_balanced),
        "control_filters": cfg.control_filters,
        "case_insensitive_filter": bool(cfg.case_insensitive_filter),
        "n_obs_source": int(len(obs)),
        "n_train_cells": int(train_mask_local.sum()),
        "n_control_cells_total": int(selected_mask_local.sum()),
        "n_control_donors_global": int(global_actual_donors),
        "shrinkage_effective_n_global": int(global_effective_n),
        "n_final_genes": int(n_genes),
        "use_shrinkage": bool(cfg.use_shrinkage),
        "shrinkage_k_donors": float(cfg.shrinkage_k_donors),
        "shrinkage_k_cells": None if cfg.shrinkage_k_cells is None else float(cfg.shrinkage_k_cells),
        "std_floor": float(cfg.std_floor),
        "chunk_size": int(cfg.chunk_size),
        "max_csr_data_span": int(cfg.max_csr_data_span),
        "max_span_to_selected_nnz_ratio": float(cfg.max_span_to_selected_nnz_ratio),
        "moment_algorithm": "one_pass_all_celltypes_float64_clamped_zarr_selected_gene_streaming",
        "check_counts_sanity": bool(cfg.check_counts_sanity),
    }

    metadata_json = json.dumps(metadata, ensure_ascii=False, sort_keys=True)

    _log(f"saving scaler npz: {out_npz}")
    np.savez_compressed(
        out_npz,
        gene_names=np.asarray(spec.gene_names, dtype=object),
        celltype_vocab=np.asarray(output_celltypes, dtype=object),
        mean_control_train=mean,
        std_control_train=std,
        n_control=n_control_cells,
        n_control_cells=n_control_cells,
        n_control_donors=n_control_donors,
        shrinkage_effective_n=shrinkage_effective_n,
        shrinkage_lambda=shrinkage_lambda,
        used_global_fallback=used_global_fallback,
        metadata_json=np.asarray(metadata_json),
        control_definition=np.asarray(json.dumps(cfg.control_filters, ensure_ascii=False)),
        donor_balanced=np.asarray(bool(cfg.donor_balanced)),
        donor_key=np.asarray(str(cfg.donor_key)),
        shrinkage_used=np.asarray(bool(cfg.use_shrinkage)),
        n_train_cells=np.asarray(int(train_mask_local.sum())),
        n_control_cells_total=np.asarray(int(selected_mask_local.sum())),
        scaler_version=np.asarray(str(cfg.scaler_version)),
        build_date=np.asarray(metadata["build_date"]),
    )

    summary = {
        "out_npz_path": str(out_npz),
        "metadata": metadata,
        "celltype_summary": [
            {
                "celltype": str(ct),
                "n_control_cells": int(n_control_cells[i]),
                "n_control_donors": int(n_control_donors[i]),
                "shrinkage_effective_n": int(shrinkage_effective_n[i]),
                "shrinkage_lambda": float(shrinkage_lambda[i]),
                "used_global_fallback": bool(used_global_fallback[i]),
            }
            for i, ct in enumerate(output_celltypes)
        ],
    }

    if cfg.out_json_path is not None:
        _log(f"saving summary json: {cfg.out_json_path}")
        with open(cfg.out_json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)

    _log("done")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    build_celltype_control_scaler(cfg)


if __name__ == "__main__":
    main()