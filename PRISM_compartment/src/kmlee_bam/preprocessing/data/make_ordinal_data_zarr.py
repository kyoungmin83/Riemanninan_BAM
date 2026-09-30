from __future__ import annotations

"""
Make ordinal specification directly from an AnnData Zarr store.

This is the Zarr-native replacement for make_ordinal_data_02.py.
It creates the same NPZ schema:
    gene_names, edges, celltype_vocab, batch_vocab, n_bins
but uses a final registry gene-union index file instead of an HVG h5ad.
"""

import argparse
import json
import os
import tempfile
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import zarr


@dataclass
class OrdinalZarrConfig:
    zarr_path: str
    gene_indices_path: str
    out_npz_path: str
    out_json_path: str

    # Matrix source: "raw/X", "X", "layers/counts", or "auto".
    matrix_path: str = "raw/X"
    counts_layer: str = "counts"
    prefer_raw: bool = True

    # Optional obs source. Use this when original zarr obs lacks split/tech_id_model.
    obs_h5ad_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_parquet_path: Optional[str] = None

    split_key: str = "split"
    train_split_name: str = "train"
    celltype_key: str = "Subclass"
    batch_key: str = "tech_id_model"

    n_bins: int = 10
    max_train_cells_for_edges: Optional[int] = 100_000
    seed: int = 42
    zero_aware_min_frac: float = 0.05
    fallback_eps: float = 1e-3

    row_chunk_size: int = 512
    tmp_dir: Optional[str] = None
    use_memmap: bool = True

    check_gene_names_json: Optional[str] = None
    fail_on_negative_counts: bool = True
    warn_if_noninteger_counts: bool = True
    integer_atol: float = 1e-3


def _log(msg: str) -> None:
    print(f"[ordinal-zarr] {msg}", flush=True)


def load_config(path: str | Path) -> OrdinalZarrConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return OrdinalZarrConfig(**raw)


def _strictly_increasing(edges: np.ndarray) -> np.ndarray:
    edges = edges.astype(np.float32, copy=True)
    for i in range(1, len(edges)):
        if edges[i] <= edges[i - 1]:
            edges[i] = np.nextafter(edges[i - 1], np.float32(np.inf))
    return edges


def _compute_gene_edges(
    values: np.ndarray,
    n_bins: int,
    *,
    fallback_eps: float = 1e-3,
    zero_aware_min_frac: float = 0.05,
) -> np.ndarray:
    """Hybrid zero-aware ordinal binning on log1p counts."""
    if values.ndim != 1:
        raise ValueError("values must be 1D.")
    if n_bins < 2:
        raise ValueError("n_bins must be >= 2.")

    values = values.astype(np.float32, copy=False)
    nonzero_mask = values > 0.0
    n = values.size
    n_nonzero = int(nonzero_mask.sum())
    zero_frac = 1.0 - (n_nonzero / max(n, 1))

    def _standard_quantile_edges(x: np.ndarray) -> np.ndarray:
        qs = np.linspace(1.0 / n_bins, (n_bins - 1) / n_bins, n_bins - 1)
        edges = np.quantile(x, qs).astype(np.float32)
        if np.unique(edges).size < edges.size:
            xmin = float(x.min())
            xmax = float(x.max())
            if xmax <= xmin:
                xmax = xmin + fallback_eps
            edges = np.linspace(xmin, xmax, n_bins + 1, dtype=np.float32)[1:-1]
        return _strictly_increasing(edges)

    if zero_frac < zero_aware_min_frac:
        return _standard_quantile_edges(values)

    if n_nonzero < 2:
        vmin = 0.0
        vmax = float(values.max()) if n_nonzero == 1 else fallback_eps
        edges = np.linspace(vmin, vmax, n_bins + 1, dtype=np.float32)[1:-1]
        return _strictly_increasing(edges)

    nz_vals = values[nonzero_mask]
    nz_min = float(nz_vals.min())
    first_edge = np.float32(nz_min / 2.0)

    n_inner = n_bins - 2
    if n_inner < 1:
        return _strictly_increasing(np.array([first_edge], dtype=np.float32))

    qs = np.linspace(1.0 / (n_inner + 1), n_inner / (n_inner + 1), n_inner)
    inner_edges = np.quantile(nz_vals, qs).astype(np.float32)

    if np.unique(inner_edges).size < inner_edges.size:
        nz_max = float(nz_vals.max())
        if nz_max <= nz_min:
            nz_max = nz_min + fallback_eps
        inner_edges = np.linspace(nz_min, nz_max, n_inner + 2, dtype=np.float32)[1:-1]

    edges = np.concatenate([[first_edge], inner_edges]).astype(np.float32)
    return _strictly_increasing(edges)


def _unique_str_values(series: pd.Series) -> np.ndarray:
    vals = series.astype(str).fillna("__NA__")
    return np.array(sorted(pd.unique(vals)), dtype=object)


def _read_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError(
            "Could not read AnnData-Zarr metadata via anndata.experimental.read_elem. "
            "Install a recent anndata version."
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


def _get_zarr_group(root, path: str):
    group = root
    for part in path.split("/"):
        if part:
            if part not in group:
                raise KeyError(f"Zarr path {path!r} not found; missing {part!r}.")
            group = group[part]
    return group


def _is_csr_group(g: Any) -> bool:
    return all(k in g for k in ["data", "indices", "indptr"])


def choose_zarr_matrix(root, cfg: OrdinalZarrConfig):
    if cfg.matrix_path != "auto":
        z_X = _get_zarr_group(root, cfg.matrix_path)
        if not _is_csr_group(z_X):
            raise ValueError(f"Matrix path {cfg.matrix_path!r} is not a CSR Zarr group.")
        var_key = "raw/var" if cfg.matrix_path.startswith("raw/") else "var"
        return z_X, cfg.matrix_path, var_key

    layer_path = f"layers/{cfg.counts_layer}"
    if "layers" in root and cfg.counts_layer in root["layers"] and _is_csr_group(root["layers"][cfg.counts_layer]):
        return root["layers"][cfg.counts_layer], layer_path, "var"
    if cfg.prefer_raw and "raw" in root and "X" in root["raw"] and _is_csr_group(root["raw"]["X"]):
        return root["raw"]["X"], "raw/X", "raw/var"
    if "X" in root and _is_csr_group(root["X"]):
        return root["X"], "X", "var"
    raise ValueError("Could not find CSR matrix in layers/counts, raw/X, or X.")


def load_obs_source(cfg: OrdinalZarrConfig, zarr_obs: pd.DataFrame) -> pd.DataFrame:
    n_sources = sum(x is not None for x in [cfg.obs_h5ad_path, cfg.obs_csv_path, cfg.obs_parquet_path])
    if n_sources > 1:
        raise ValueError("Provide only one of obs_h5ad_path, obs_csv_path, obs_parquet_path.")

    if cfg.obs_h5ad_path is not None:
        import anndata as ad  # local import
        _log(f"loading external obs from h5ad: {cfg.obs_h5ad_path}")
        return ad.read_h5ad(cfg.obs_h5ad_path, backed="r").obs.copy()
    if cfg.obs_csv_path is not None:
        _log(f"loading external obs from csv: {cfg.obs_csv_path}")
        return pd.read_csv(cfg.obs_csv_path, index_col=0)
    if cfg.obs_parquet_path is not None:
        _log(f"loading external obs from parquet: {cfg.obs_parquet_path}")
        return pd.read_parquet(cfg.obs_parquet_path)
    return zarr_obs.copy()


def align_obs_to_zarr(obs_source: pd.DataFrame, zarr_obs: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    zarr_index = pd.Index(zarr_obs.index.astype(str))
    source_index = pd.Index(obs_source.index.astype(str))

    if len(source_index) == len(zarr_index) and np.array_equal(source_index.to_numpy(), zarr_index.to_numpy()):
        return obs_source.copy(), np.arange(len(zarr_index), dtype=np.int64)

    if source_index.has_duplicates:
        dups = source_index[source_index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"obs source contains duplicated barcodes. Examples: {dups}")

    pos = zarr_index.get_indexer(source_index)
    missing = pos < 0
    if missing.any():
        examples = source_index[missing].tolist()[:5]
        raise ValueError(f"Some obs barcodes not found in zarr obs. Examples: {examples}")
    return obs_source.copy(), pos.astype(np.int64)


def read_var_names(root, var_key: str) -> np.ndarray:
    var = read_zarr_dataframe(root, var_key)
    return np.asarray(var.index.astype(str), dtype=object)


def infer_csr_shape(z_X, n_obs: int, n_vars: int) -> tuple[int, int]:
    shape = z_X.attrs.get("shape", None)
    if shape is not None:
        return int(shape[0]), int(shape[1])
    return n_obs, n_vars


def fill_log1p_selected_matrix(
    *,
    z_X,
    rows: np.ndarray,
    gene_indices: np.ndarray,
    n_vars: int,
    out: np.ndarray,
    row_chunk_size: int,
    fail_on_negative_counts: bool,
    warn_if_noninteger_counts: bool,
    integer_atol: float,
) -> None:
    """Fill out [len(rows), len(gene_indices)] with log1p selected raw counts."""
    if out.shape != (len(rows), len(gene_indices)):
        raise ValueError("Output matrix shape mismatch.")

    full_to_sel = np.full(n_vars, -1, dtype=np.int32)
    full_to_sel[gene_indices.astype(np.int64)] = np.arange(len(gene_indices), dtype=np.int32)

    data_arr = z_X["data"]
    indices_arr = z_X["indices"]
    indptr = z_X["indptr"][:]

    order = np.argsort(rows)
    rows_sorted = rows[order]
    warned_noninteger = False

    for start in range(0, len(rows_sorted), row_chunk_size):
        end = min(start + row_chunk_size, len(rows_sorted))
        rr = rows_sorted[start:end]
        chunk = np.zeros((len(rr), len(gene_indices)), dtype=np.float32)

        for local_i, r in enumerate(rr):
            s = int(indptr[r])
            e = int(indptr[r + 1])
            if e <= s:
                continue
            idx = indices_arr[s:e]
            val = data_arr[s:e].astype(np.float32, copy=False)

            if fail_on_negative_counts and np.any(val < 0):
                raise ValueError(f"Negative count values found in row {r}.")
            if warn_if_noninteger_counts and not warned_noninteger:
                nz = val[np.abs(val) > 0]
                if nz.size > 0:
                    frac = float(np.mean(np.isclose(nz, np.round(nz), atol=integer_atol, rtol=0.0)))
                    if frac < 0.95:
                        warnings.warn(
                            f"Row {r} contains non-integer-like count values; check matrix_path.",
                            RuntimeWarning,
                        )
                        warned_noninteger = True

            mapped = full_to_sel[idx.astype(np.int64, copy=False)]
            keep = mapped >= 0
            if np.any(keep):
                chunk[local_i, mapped[keep]] = val[keep]

        np.log1p(chunk, out=chunk)
        out[order[start:end], :] = chunk

        if (end % (row_chunk_size * 10) == 0) or end == len(rows_sorted):
            _log(f"  extracted log1p rows: {end:,}/{len(rows_sorted):,}")


def fit_ordinal_spec_zarr(cfg: OrdinalZarrConfig) -> dict[str, Any]:
    if cfg.n_bins < 2:
        raise ValueError("n_bins must be >= 2.")
    if cfg.max_train_cells_for_edges is not None and cfg.max_train_cells_for_edges <= 0:
        raise ValueError("max_train_cells_for_edges must be positive or null.")
    if cfg.row_chunk_size <= 0:
        raise ValueError("row_chunk_size must be positive.")

    out_npz = Path(cfg.out_npz_path)
    out_json = Path(cfg.out_json_path)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    _log(f"opening zarr: {cfg.zarr_path}")
    root = zarr.open(cfg.zarr_path, mode="r")

    _log("reading zarr metadata")
    zarr_obs = read_zarr_dataframe(root, "obs")
    z_X, matrix_source, var_key = choose_zarr_matrix(root, cfg)
    var_names = read_var_names(root, var_key)
    _log(f"using matrix={matrix_source}; var={var_key}")

    n_obs, n_vars = infer_csr_shape(z_X, n_obs=len(zarr_obs), n_vars=len(var_names))
    if n_vars != len(var_names):
        raise ValueError(f"Matrix n_vars={n_vars} but var_names length={len(var_names)}.")

    obs_source = load_obs_source(cfg, zarr_obs)
    obs, zarr_row_positions = align_obs_to_zarr(obs_source, zarr_obs)

    for key in [cfg.split_key, cfg.celltype_key, cfg.batch_key]:
        if key not in obs.columns:
            raise KeyError(f"'{key}' not found in obs source.")

    gene_indices = np.load(cfg.gene_indices_path).astype(np.int64)
    if gene_indices.ndim != 1 or len(gene_indices) == 0:
        raise ValueError("gene_indices_path must contain a non-empty 1D numpy array.")
    if gene_indices.min() < 0 or gene_indices.max() >= n_vars:
        raise ValueError(
            f"gene_indices out of bounds for matrix with n_vars={n_vars}: "
            f"min={gene_indices.min()}, max={gene_indices.max()}."
        )
    if len(np.unique(gene_indices)) != len(gene_indices):
        raise ValueError("gene_indices contains duplicates.")

    gene_names = np.asarray(var_names[gene_indices], dtype=object)

    if cfg.check_gene_names_json is not None:
        with open(cfg.check_gene_names_json, "r", encoding="utf-8") as f:
            expected = np.asarray(json.load(f), dtype=object)
        if not np.array_equal(expected, gene_names):
            raise ValueError(
                "Gene names from zarr/gene_indices do not match check_gene_names_json. "
                "Check matrix_path and gene_union_original_indices.npy."
            )

    split_vals = obs[cfg.split_key].astype(str).to_numpy()
    train_local = np.where(split_vals == str(cfg.train_split_name))[0]
    if len(train_local) == 0:
        raise ValueError(f"No cells found for split={cfg.train_split_name!r}.")

    rng = np.random.default_rng(cfg.seed)
    n_train_available = len(train_local)
    if cfg.max_train_cells_for_edges is not None and len(train_local) > cfg.max_train_cells_for_edges:
        train_local = np.sort(rng.choice(train_local, size=int(cfg.max_train_cells_for_edges), replace=False))

    train_zarr_rows = zarr_row_positions[train_local].astype(np.int64)
    n_train = len(train_zarr_rows)
    n_genes = len(gene_indices)
    _log(f"fitting ordinal spec from {n_train:,} train cells x {n_genes:,} final genes")

    tmp_path = None
    if cfg.use_memmap:
        tmp_dir = cfg.tmp_dir if cfg.tmp_dir is not None else tempfile.gettempdir()
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        tmp_path = Path(tmp_dir) / f"ordinal_zarr_log1p_{os.getpid()}_{out_npz.stem}.dat"
        _log(f"using memmap: {tmp_path}")
        X_log = np.memmap(tmp_path, mode="w+", dtype=np.float32, shape=(n_train, n_genes))
    else:
        X_log = np.zeros((n_train, n_genes), dtype=np.float32)

    fill_log1p_selected_matrix(
        z_X=z_X,
        rows=train_zarr_rows,
        gene_indices=gene_indices,
        n_vars=n_vars,
        out=X_log,
        row_chunk_size=cfg.row_chunk_size,
        fail_on_negative_counts=cfg.fail_on_negative_counts,
        warn_if_noninteger_counts=cfg.warn_if_noninteger_counts,
        integer_atol=cfg.integer_atol,
    )

    edges = np.zeros((n_genes, cfg.n_bins - 1), dtype=np.float32)
    for g in range(n_genes):
        values = np.asarray(X_log[:, g], dtype=np.float32)
        edges[g] = _compute_gene_edges(
            values,
            n_bins=cfg.n_bins,
            fallback_eps=cfg.fallback_eps,
            zero_aware_min_frac=cfg.zero_aware_min_frac,
        )
        if (g + 1) % 500 == 0 or (g + 1) == n_genes:
            _log(f"  processed genes: {g + 1:,}/{n_genes:,}")

    if isinstance(X_log, np.memmap):
        X_log.flush()
        if tmp_path is not None:
            try:
                Path(tmp_path).unlink(missing_ok=True)
            except Exception:
                pass

    celltype_vocab = _unique_str_values(obs[cfg.celltype_key])
    batch_vocab = _unique_str_values(obs[cfg.batch_key])

    _log(f"saving ordinal spec: {out_npz}")
    np.savez_compressed(
        out_npz,
        gene_names=gene_names,
        edges=edges,
        celltype_vocab=celltype_vocab,
        batch_vocab=batch_vocab,
        n_bins=np.array([cfg.n_bins], dtype=np.int32),
    )

    manifest = {
        "config": asdict(cfg),
        "zarr_path": cfg.zarr_path,
        "matrix_source": matrix_source,
        "var_key": var_key,
        "gene_indices_path": cfg.gene_indices_path,
        "n_obs_zarr": int(n_obs),
        "n_obs_source": int(len(obs)),
        "n_genes_matrix": int(n_vars),
        "n_genes_final": int(n_genes),
        "n_train_cells_available": int(n_train_available),
        "n_train_cells_used_for_edges": int(n_train),
        "split_key": cfg.split_key,
        "train_split_name": cfg.train_split_name,
        "celltype_key": cfg.celltype_key,
        "batch_key": cfg.batch_key,
        "n_bins": int(cfg.n_bins),
    }

    _log(f"saving manifest: {out_json}")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    _log("done")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()

    cfg = load_config(args.config)
    fit_ordinal_spec_zarr(cfg)


if __name__ == "__main__":
    main()