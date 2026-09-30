from __future__ import annotations

"""
build_gene_union_subset_zarr_hybrid.py
======================================

Hybrid fast/safe builder for a final-gene-union AnnData-Zarr store.

This script combines two desirable properties:

1. Fast extraction logic:
   - Read the original full-gene CSR matrix once.
   - Build chunk-level subset CSR matrices.
   - Append subset data/indices/indptr to a new Zarr store.
   - Verify sampled rows against the original source matrix.

2. Pipeline safety:
   - Merge pathology-aware obs sidecar into the output obs.
   - Preserve split / fold / pathology / tech columns.
   - Write an identity gene-index file for downstream subset-zarr configs:
         gene_union_subset_indices.npy = np.arange(G_final)
   - Write a manifest with exact source/output paths.

The output is NOT a cell x module matrix. It is still gene-level raw counts:

    cells x final_gene_union_genes

Downstream after this script:
    zarr_path         = out_zarr_path
    matrix_path       = "X"
    gene_indices_path = gene_union_subset_indices.npy
"""

import argparse
import json
import shutil
import time
import warnings
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr


# =============================================================================
# Config
# =============================================================================
@dataclass
class GeneUnionSubsetZarrHybridConfig:
    # Input source. `original_zarr_path` is accepted as an alias in load_config.
    zarr_path: str
    gene_indices_path: str
    out_zarr_path: str

    # Source matrix inside original zarr.
    matrix_path: str = "raw/X"
    counts_layer: str = "counts"
    prefer_raw: bool = True

    # Optional pathology-aware obs sidecar.
    obs_parquet_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_h5ad_path: Optional[str] = None
    merge_obs_sidecar: bool = True

    # Optional sanity check against final_registry/gene_union_ensembl.json.
    check_gene_names_json: Optional[str] = None

    # Output side files.
    out_gene_indices_path: Optional[str] = None
    out_gene_names_json: Optional[str] = None
    out_manifest_path: Optional[str] = None

    # Chunking / memory.
    row_chunk_size: int = 8192
    data_chunk_size: int = 2_000_000
    indptr_chunk_size: int = 1_000_000

    # Span heuristic. Because chunks here are contiguous rows, use_span is
    # usually safe and fast, but we keep bounds to avoid pathological memory.
    max_csr_data_span: int = 200_000_000
    max_span_to_selected_nnz_ratio: float = 20.0

    # Optional raw alias. Usually false; downstream should use matrix_path="X".
    write_raw_alias: bool = False

    # Safety / verification.
    overwrite: bool = False
    verify_n_rows: int = 64

    # Compression. None means copy source compressor when available.
    compression_level: Optional[int] = None


# =============================================================================
# Logging / config
# =============================================================================
def _log(msg: str) -> None:
    print(f"[gene-union-subset-hybrid] {msg}", flush=True)


class Timer:
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


def load_config(path: str | Path) -> GeneUnionSubsetZarrHybridConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    # Alias for compatibility with the other proposal.
    if "zarr_path" not in raw and "original_zarr_path" in raw:
        raw["zarr_path"] = raw.pop("original_zarr_path")

    # Alias for compatibility with the other proposal.
    if "data_chunk_size" in raw and "zarr_data_chunk" not in raw:
        raw["data_chunk_size"] = raw.pop("data_chunk_size")
    if "indptr_chunk_size" in raw and "zarr_indptr_chunk" not in raw:
        raw["indptr_chunk_size"] = raw.pop("indptr_chunk_size")

    return GeneUnionSubsetZarrHybridConfig(**raw)


# =============================================================================
# Zarr helpers
# =============================================================================
def _read_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError(
            "Could not read AnnData-Zarr metadata with anndata.experimental.read_elem."
        ) from e


def _write_elem(root, key: str, value) -> None:
    try:
        from anndata.experimental import write_elem  # type: ignore
        write_elem(root, key, value)
    except Exception as e:
        raise RuntimeError(
            f"Could not write AnnData-Zarr element {key!r} with anndata.experimental.write_elem."
        ) from e


def _zarr_get(root, path: str):
    obj = root
    for part in path.split("/"):
        if part:
            if part not in obj:
                raise KeyError(f"Zarr path {path!r} missing component {part!r}.")
            obj = obj[part]
    return obj


def read_zarr_dataframe(root, key_path: str) -> pd.DataFrame:
    g = _zarr_get(root, key_path)
    out = _read_elem(g)
    if not isinstance(out, pd.DataFrame):
        raise TypeError(f"Zarr element {key_path!r} is not a pandas DataFrame.")
    return out.copy()


def choose_zarr_matrix(root, cfg: GeneUnionSubsetZarrHybridConfig):
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
    raise KeyError("Could not find matrix in original zarr.")


def is_csr_group(g) -> bool:
    return hasattr(g, "keys") and all(k in g for k in ("data", "indices", "indptr"))


def load_obs_sidecar(cfg: GeneUnionSubsetZarrHybridConfig) -> Optional[pd.DataFrame]:
    n_src = sum(x is not None for x in [cfg.obs_parquet_path, cfg.obs_csv_path, cfg.obs_h5ad_path])
    if n_src > 1:
        raise ValueError("Provide at most one obs sidecar source.")

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

    return None


def merge_obs_sidecar(zarr_obs: pd.DataFrame, sidecar: Optional[pd.DataFrame], *, merge: bool) -> pd.DataFrame:
    zarr_obs = zarr_obs.copy()
    zarr_obs.index = zarr_obs.index.astype(str)

    if sidecar is None:
        return zarr_obs

    sidecar = sidecar.copy()
    sidecar.index = sidecar.index.astype(str)

    zidx = pd.Index(zarr_obs.index.astype(str))
    sidx = pd.Index(sidecar.index.astype(str))

    if sidx.has_duplicates:
        dup = sidx[sidx.duplicated()].unique().tolist()[:5]
        raise ValueError(f"obs sidecar has duplicated barcodes. Examples: {dup}")

    missing = sidx.difference(zidx)
    if len(missing) > 0:
        raise ValueError(f"obs sidecar barcodes not found in zarr obs. Examples: {missing[:5].tolist()}")

    aligned = sidecar.reindex(zidx)

    if not merge:
        if len(sidecar) != len(zarr_obs):
            raise ValueError("merge_obs_sidecar=false requires sidecar to contain all zarr cells.")
        return aligned

    out = zarr_obs.copy()
    # Sidecar columns overwrite original columns when names overlap.
    for c in aligned.columns:
        out[c] = aligned[c].values
    return out


def maybe_get_compressor(src_arr, cfg: GeneUnionSubsetZarrHybridConfig):
    comp = getattr(src_arr, "compressor", None)
    if cfg.compression_level is None or comp is None:
        return comp

    try:
        from numcodecs import Blosc
        if isinstance(comp, Blosc):
            return Blosc(cname=comp.cname, clevel=int(cfg.compression_level), shuffle=comp.shuffle)
    except Exception:
        pass

    return comp


def prepare_output_root(out_path: Path, *, overwrite: bool):
    if out_path.exists():
        if not overwrite:
            raise FileExistsError(f"{out_path} already exists. Set overwrite=true to replace it.")
        _log(f"removing existing output: {out_path}")
        if out_path.is_dir():
            shutil.rmtree(out_path)
        else:
            out_path.unlink()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    root = zarr.open(str(out_path), mode="w")
    root.attrs["encoding-type"] = "anndata"
    root.attrs["encoding-version"] = "0.1.0"
    return root


def create_csr_group(
    root,
    group_name: str,
    *,
    n_rows: int,
    n_cols: int,
    data_dtype,
    indices_dtype,
    indptr_dtype,
    data_chunk_size: int,
    indptr_chunk_size: int,
    compressor,
):
    if group_name in root:
        del root[group_name]

    g = root.create_group(group_name)
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = [int(n_rows), int(n_cols)]

    data_z = g.create_dataset(
        "data",
        shape=(0,),
        chunks=(int(data_chunk_size),),
        dtype=data_dtype,
        compressor=compressor,
    )
    indices_z = g.create_dataset(
        "indices",
        shape=(0,),
        chunks=(int(data_chunk_size),),
        dtype=indices_dtype,
        compressor=compressor,
    )
    indptr_z = g.create_dataset(
        "indptr",
        shape=(n_rows + 1,),
        chunks=(min(int(indptr_chunk_size), n_rows + 1),),
        dtype=indptr_dtype,
        compressor=compressor,
    )
    return g, data_z, indices_z, indptr_z


def append_1d(arr, values: np.ndarray) -> None:
    values = np.asarray(values)
    if values.size == 0:
        return
    old = int(arr.shape[0])
    new = old + int(values.shape[0])
    arr.resize((new,))
    arr[old:new] = values


# =============================================================================
# Source CSR reader / chunk builder
# =============================================================================
class SourceCSRReader:
    def __init__(self, matrix, *, n_cols_full: int, name: str) -> None:
        if not is_csr_group(matrix):
            raise TypeError(
                f"Matrix at {name!r} is not a CSR group. This script currently supports CSR sources only."
            )

        self.name = str(name)
        self.data = matrix["data"]
        self.indices = matrix["indices"]
        self.indptr = matrix["indptr"][:].astype(np.int64)
        self.n_cols_full = int(n_cols_full)
        self.shape = (len(self.indptr) - 1, self.n_cols_full)


def build_subset_csr_for_chunk(
    *,
    reader: SourceCSRReader,
    rr: np.ndarray,
    full_to_selected: np.ndarray,
    n_final_genes: int,
    max_csr_data_span: int,
    max_span_to_selected_nnz_ratio: float,
) -> sp.csr_matrix:
    """
    Build [len(rr), n_final_genes] CSR chunk from source CSR rows.

    The fast path reads one contiguous CSR data span. The robust mapping then
    ensures only entries belonging to the requested rows rr are retained. For
    contiguous rr this is equivalent to a fast cumulative-count method; for
    future non-contiguous rr it remains correct.
    """
    indptr = reader.indptr
    rr = np.asarray(rr, dtype=np.int64)

    if len(rr) == 0:
        return sp.csr_matrix((0, n_final_genes), dtype=np.float32)

    row_starts_full = indptr[rr].astype(np.int64)
    row_ends_full = indptr[rr + 1].astype(np.int64)
    selected_nnz = int((row_ends_full - row_starts_full).sum())

    p0 = int(row_starts_full[0])
    p1 = int(row_ends_full[-1])
    span = max(0, p1 - p0)

    use_span = (
        span > 0
        and span <= int(max_csr_data_span)
        and span <= float(max_span_to_selected_nnz_ratio) * max(selected_nnz, 1)
    )

    if use_span:
        idx_span = np.asarray(reader.indices[p0:p1], dtype=np.int64)
        dat_span = np.asarray(reader.data[p0:p1], dtype=np.float32)

        mapped_span = full_to_selected[idx_span]
        gene_keep = mapped_span >= 0

        if not np.any(gene_keep):
            return sp.csr_matrix((len(rr), n_final_genes), dtype=np.float32)

        kept_pos = np.flatnonzero(gene_keep).astype(np.int64)

        row_starts_local = row_starts_full - p0
        row_ends_local = row_ends_full - p0

        # Fast path for contiguous row chunks.
        if int(rr[-1]) - int(rr[0]) + 1 == len(rr):
            keep_i64 = gene_keep.astype(np.int64)
            keep_cum = np.empty(len(keep_i64) + 1, dtype=np.int64)
            keep_cum[0] = 0
            np.cumsum(keep_i64, out=keep_cum[1:])

            kept_per_row = keep_cum[row_ends_local] - keep_cum[row_starts_local]
            chunk_indptr = np.empty(len(rr) + 1, dtype=np.int64)
            chunk_indptr[0] = 0
            np.cumsum(kept_per_row, out=chunk_indptr[1:])

            mapped_kept = mapped_span[gene_keep].astype(np.int32)
            data_kept = dat_span[gene_keep].astype(np.float32, copy=False)

            return sp.csr_matrix(
                (data_kept, mapped_kept, chunk_indptr),
                shape=(len(rr), n_final_genes),
                dtype=np.float32,
            )

        # Robust path for non-contiguous row selections.
        candidate = np.searchsorted(row_ends_local, kept_pos, side="right")
        in_range = candidate < len(rr)
        if not np.any(in_range):
            return sp.csr_matrix((len(rr), n_final_genes), dtype=np.float32)

        candidate2 = candidate[in_range]
        kept_pos2 = kept_pos[in_range]
        valid = kept_pos2 >= row_starts_local[candidate2]
        if not np.any(valid):
            return sp.csr_matrix((len(rr), n_final_genes), dtype=np.float32)

        row_idx = candidate2[valid].astype(np.int32)
        pos_final = kept_pos2[valid]
        col_idx = mapped_span[pos_final].astype(np.int32)
        data_kept = dat_span[pos_final].astype(np.float32, copy=False)

        return sp.csr_matrix(
            (data_kept, (row_idx, col_idx)),
            shape=(len(rr), n_final_genes),
            dtype=np.float32,
        )

    # Fallback: per-row reads, only when span is too large/wasteful.
    rows_data: list[np.ndarray] = []
    rows_idx: list[np.ndarray] = []
    rows_ptr = [0]

    for r in rr:
        s = int(indptr[int(r)])
        e = int(indptr[int(r) + 1])
        if e > s:
            idx = np.asarray(reader.indices[s:e], dtype=np.int64)
            dat = np.asarray(reader.data[s:e], dtype=np.float32)
            mapped = full_to_selected[idx]
            keep = mapped >= 0
            if np.any(keep):
                rows_data.append(dat[keep])
                rows_idx.append(mapped[keep].astype(np.int32))
                rows_ptr.append(rows_ptr[-1] + int(keep.sum()))
                continue
        rows_ptr.append(rows_ptr[-1])

    if rows_ptr[-1] == 0:
        return sp.csr_matrix((len(rr), n_final_genes), dtype=np.float32)

    return sp.csr_matrix(
        (
            np.concatenate(rows_data),
            np.concatenate(rows_idx),
            np.asarray(rows_ptr, dtype=np.int64),
        ),
        shape=(len(rr), n_final_genes),
        dtype=np.float32,
    )


# =============================================================================
# Verification
# =============================================================================
def verify_subset(
    *,
    reader: SourceCSRReader,
    out_zarr_path: str,
    full_to_selected: np.ndarray,
    n_final_genes: int,
    n_check: int,
    seed: int = 0,
) -> dict[str, Any]:
    if n_check <= 0:
        return {"n_rows_checked": 0, "mismatches": 0}

    _log(f"verifying {n_check} random rows against source")
    rng = np.random.default_rng(seed)
    n_cells = reader.shape[0]
    rows = np.sort(rng.choice(n_cells, size=min(n_check, n_cells), replace=False))

    dst = zarr.open(str(out_zarr_path), mode="r")
    dst_indptr = dst["X"]["indptr"][:]
    dst_data = dst["X"]["data"]
    dst_indices = dst["X"]["indices"]

    mismatches = 0
    first_mismatch: Optional[int] = None

    for r in rows:
        s = int(reader.indptr[int(r)])
        e = int(reader.indptr[int(r) + 1])
        src_idx = np.asarray(reader.indices[s:e], dtype=np.int64)
        src_dat = np.asarray(reader.data[s:e], dtype=np.float32)

        mapped = full_to_selected[src_idx]
        keep = mapped >= 0
        ref_idx = mapped[keep].astype(np.int32)
        ref_dat = src_dat[keep]

        order = np.argsort(ref_idx, kind="stable")
        ref_idx = ref_idx[order]
        ref_dat = ref_dat[order]

        ds = int(dst_indptr[int(r)])
        de = int(dst_indptr[int(r) + 1])
        got_idx = np.asarray(dst_indices[ds:de], dtype=np.int32)
        got_dat = np.asarray(dst_data[ds:de], dtype=np.float32)

        order2 = np.argsort(got_idx, kind="stable")
        got_idx = got_idx[order2]
        got_dat = got_dat[order2]

        ok = (
            got_idx.shape == ref_idx.shape
            and np.array_equal(got_idx, ref_idx)
            and np.allclose(got_dat, ref_dat, rtol=0, atol=0)
        )
        if not ok:
            mismatches += 1
            if first_mismatch is None:
                first_mismatch = int(r)

    if mismatches > 0:
        warnings.warn(
            f"verification failed: {mismatches}/{len(rows)} sampled rows mismatched. "
            f"first_mismatch={first_mismatch}",
            RuntimeWarning,
        )
    else:
        _log(f"verification passed on {len(rows)} rows")

    return {
        "n_rows_checked": int(len(rows)),
        "mismatches": int(mismatches),
        "first_mismatch": first_mismatch,
    }


# =============================================================================
# Main
# =============================================================================
def build_subset_zarr(cfg: GeneUnionSubsetZarrHybridConfig) -> dict[str, Any]:
    out_path = Path(cfg.out_zarr_path)

    if cfg.row_chunk_size <= 0:
        raise ValueError("row_chunk_size must be positive.")
    if cfg.data_chunk_size <= 0:
        raise ValueError("data_chunk_size must be positive.")
    if cfg.indptr_chunk_size <= 0:
        raise ValueError("indptr_chunk_size must be positive.")
    if cfg.max_csr_data_span <= 0:
        raise ValueError("max_csr_data_span must be positive.")

    with Timer("open source zarr and metadata"):
        src = zarr.open(cfg.zarr_path, mode="r")
        src_obs = read_zarr_dataframe(src, "obs")
        matrix, matrix_source, var_key = choose_zarr_matrix(src, cfg)
        if not is_csr_group(matrix):
            raise TypeError(f"Source matrix {matrix_source!r} is not CSR.")

        src_var = read_zarr_dataframe(src, var_key)
        n_cells = len(src_obs)
        n_full_genes = len(src_var)

        gene_indices = np.load(cfg.gene_indices_path).astype(np.int64)
        if gene_indices.ndim != 1 or len(gene_indices) == 0:
            raise ValueError("gene_indices must be a non-empty 1D array.")
        if gene_indices.min() < 0 or gene_indices.max() >= n_full_genes:
            raise ValueError("gene_indices out of bounds.")
        if len(np.unique(gene_indices)) != len(gene_indices):
            raise ValueError("gene_indices contains duplicates.")

        n_final_genes = int(len(gene_indices))
        src_var_subset = src_var.iloc[gene_indices].copy()
        gene_names = np.asarray(src_var_subset.index.astype(str), dtype=object)

        if cfg.check_gene_names_json is not None:
            with open(cfg.check_gene_names_json, "r", encoding="utf-8") as f:
                expected = np.asarray(json.load(f), dtype=object)
            if not np.array_equal(expected, gene_names):
                raise ValueError(
                    "Selected source gene names do not match check_gene_names_json. "
                    "Check matrix_path and gene_indices_path."
                )

        sidecar = load_obs_sidecar(cfg)
        obs_out = merge_obs_sidecar(src_obs, sidecar, merge=cfg.merge_obs_sidecar)

        if len(obs_out) != n_cells:
            raise ValueError(
                "This builder creates an all-cell subset zarr. obs_out length must equal source n_cells."
            )

        reader = SourceCSRReader(matrix, n_cols_full=n_full_genes, name=matrix_source)

    _log(f"source matrix={matrix_source}; cells={n_cells:,}; full genes={n_full_genes:,}")
    _log(f"selected final genes={n_final_genes:,}")

    with Timer("create output zarr and metadata"):
        root = prepare_output_root(out_path, overwrite=cfg.overwrite)
        _write_elem(root, "obs", obs_out)
        _write_elem(root, "var", src_var_subset)

        compressor = maybe_get_compressor(matrix["data"], cfg)
        _, data_z, indices_z, indptr_z = create_csr_group(
            root,
            "X",
            n_rows=n_cells,
            n_cols=n_final_genes,
            data_dtype=matrix["data"].dtype,
            indices_dtype=np.int32,
            indptr_dtype=np.int64,
            data_chunk_size=cfg.data_chunk_size,
            indptr_chunk_size=cfg.indptr_chunk_size,
            compressor=compressor,
        )

    full_to_selected = np.full(n_full_genes, -1, dtype=np.int32)
    full_to_selected[gene_indices] = np.arange(n_final_genes, dtype=np.int32)

    out_indptr = np.empty(n_cells + 1, dtype=np.int64)
    out_indptr[0] = 0
    nnz_out = 0

    with Timer("extract selected-gene CSR matrix"):
        row = 0
        while row < n_cells:
            stop = min(row + int(cfg.row_chunk_size), n_cells)

            # Shrink row block if the source CSR span is too large.
            while stop > row + 1:
                p0_test = int(reader.indptr[row])
                p1_test = int(reader.indptr[stop])
                if (p1_test - p0_test) <= int(cfg.max_csr_data_span):
                    break
                stop = row + max(1, (stop - row) // 2)

            rr = np.arange(row, stop, dtype=np.int64)

            chunk_csr = build_subset_csr_for_chunk(
                reader=reader,
                rr=rr,
                full_to_selected=full_to_selected,
                n_final_genes=n_final_genes,
                max_csr_data_span=cfg.max_csr_data_span,
                max_span_to_selected_nnz_ratio=cfg.max_span_to_selected_nnz_ratio,
            )

            chunk_data = chunk_csr.data.astype(matrix["data"].dtype, copy=False)
            chunk_indices = chunk_csr.indices.astype(np.int32, copy=False)

            append_1d(data_z, chunk_data)
            append_1d(indices_z, chunk_indices)

            out_indptr[row + 1:stop + 1] = nnz_out + chunk_csr.indptr[1:].astype(np.int64, copy=False)
            nnz_out = int(np.int64(nnz_out) + np.int64(chunk_csr.nnz))

            row = stop
            if row == n_cells or (row % (int(cfg.row_chunk_size) * 10) == 0):
                _log(f"  rows processed: {row:,}/{n_cells:,}; subset nnz={nnz_out:,}")

    with Timer("write indptr and optional raw alias"):
        indptr_z[:] = out_indptr
        root.attrs["shape"] = [int(n_cells), int(n_final_genes)]

        if cfg.write_raw_alias:
            _log("write_raw_alias=true: duplicating X into raw/X. This doubles matrix storage.")
            if "raw" in root:
                del root["raw"]
            raw = root.create_group("raw")
            raw.attrs["encoding-type"] = "raw"
            raw.attrs["encoding-version"] = "0.1.0"
            _write_elem(raw, "var", src_var_subset)

            _, raw_data, raw_indices, raw_indptr = create_csr_group(
                raw,
                "X",
                n_rows=n_cells,
                n_cols=n_final_genes,
                data_dtype=data_z.dtype,
                indices_dtype=indices_z.dtype,
                indptr_dtype=indptr_z.dtype,
                data_chunk_size=cfg.data_chunk_size,
                indptr_chunk_size=cfg.indptr_chunk_size,
                compressor=maybe_get_compressor(matrix["data"], cfg),
            )
            raw_data.resize(data_z.shape)
            raw_indices.resize(indices_z.shape)
            raw_data[:] = data_z[:]
            raw_indices[:] = indices_z[:]
            raw_indptr[:] = indptr_z[:]

    with Timer("verification"):
        verification = verify_subset(
            reader=reader,
            out_zarr_path=str(out_path),
            full_to_selected=full_to_selected,
            n_final_genes=n_final_genes,
            n_check=int(cfg.verify_n_rows),
        )

    out_dir = out_path.parent
    out_gene_indices_path = Path(cfg.out_gene_indices_path) if cfg.out_gene_indices_path else out_dir / "gene_union_subset_indices.npy"
    out_gene_names_json = Path(cfg.out_gene_names_json) if cfg.out_gene_names_json else out_dir / "gene_union_subset_gene_names.json"
    out_manifest_path = Path(cfg.out_manifest_path) if cfg.out_manifest_path else out_dir / "gene_union_subset_manifest.json"

    with Timer("write side files"):
        out_gene_indices_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(out_gene_indices_path, np.arange(n_final_genes, dtype=np.int64))

        with open(out_gene_names_json, "w", encoding="utf-8") as f:
            json.dump([str(x) for x in gene_names.tolist()], f, ensure_ascii=False, indent=2)

        manifest = {
            "config": asdict(cfg),
            "build_date": datetime.now(timezone.utc).isoformat(),
            "source_zarr_path": str(cfg.zarr_path),
            "source_matrix_path": str(cfg.matrix_path),
            "source_matrix_resolved": str(matrix_source),
            "source_var_key": str(var_key),
            "out_zarr_path": str(out_path),
            "out_matrix_path": "X",
            "n_obs": int(n_cells),
            "n_full_genes": int(n_full_genes),
            "n_final_genes": int(n_final_genes),
            "nnz_out": int(nnz_out),
            "density": float(nnz_out / max(n_cells * n_final_genes, 1)),
            "gene_indices_path_original": str(cfg.gene_indices_path),
            "gene_indices_path_for_subset": str(out_gene_indices_path),
            "gene_names_json_for_subset": str(out_gene_names_json),
            "verification": verification,
            "obs_columns": list(map(str, obs_out.columns)),
            "downstream_config_note": {
                "zarr_path": str(out_path),
                "matrix_path": "X",
                "gene_indices_path": str(out_gene_indices_path),
            },
        }

        with open(out_manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)

    _log("done")
    _log(f"subset zarr: {out_path}")
    _log(f"subset gene indices: {out_gene_indices_path}")
    return manifest


# =============================================================================
# CLI
# =============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Build final gene-union subset AnnData-Zarr.")
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()

    cfg = load_config(args.config)
    build_subset_zarr(cfg)


if __name__ == "__main__":
    main()