from __future__ import annotations

"""
Build a CP10K-normalized AnnData Zarr from a raw-count AnnData Zarr.

This is the canonical preprocessing step that turns the gene-union raw-count
store produced by build_gene_union_subset_zarr.py into the CP10K expression
store consumed by training configs as data.zarr_path.

It intentionally does only one thing:

    raw counts / obs["Number of UMIs"] * 10_000 -> CP10K sparse X

Ordinal spec generation belongs in make_ordinal_data_zarr.py.
Reference scaler generation belongs in build_celltype_control_scaler.py.
"""

import argparse
import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import zarr


@dataclass
class BuildCp10kZarrConfig:
    src_zarr_path: str
    out_zarr_path: str

    # Source sparse matrix inside src_zarr_path. For the current final registry
    # raw-count companion this is "X".
    source_matrix: str = "X"

    # If null, choose from common UMI columns in obs.
    umi_col: Optional[str] = "Number of UMIs"

    scale_factor: float = 10_000.0
    row_block_size: int = 4096
    overwrite: bool = False

    # Optional manifest with exact provenance. If null, write next to out_zarr.
    out_manifest_path: Optional[str] = None


def _log(msg: str) -> None:
    print(f"[build-cp10k-zarr] {msg}", flush=True)


def load_config(path: str | Path) -> BuildCp10kZarrConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return BuildCp10kZarrConfig(**raw)


def _get_group(root, path: str):
    obj = root
    for part in path.split("/"):
        if part:
            obj = obj[part]
    return obj


def _read_obs_var(zarr_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    import anndata as ad

    adata = ad.read_zarr(zarr_path)
    return adata.obs.copy(), adata.var.copy()


def _choose_umi(obs: pd.DataFrame, umi_col: Optional[str]) -> tuple[str, np.ndarray]:
    if umi_col is not None:
        if umi_col not in obs.columns:
            raise KeyError(f"umi_col={umi_col!r} not found in obs.")
        return umi_col, obs[umi_col].to_numpy(dtype=np.float64)

    candidates = ["Number of UMIs", "total_counts", "n_counts", "UMI_count"]
    for col in candidates:
        if col in obs.columns:
            _log(f"Using obs[{col!r}] as total UMI column.")
            return col, obs[col].to_numpy(dtype=np.float64)

    raise KeyError(
        "Could not find a total UMI column in obs. "
        "Pass umi_col explicitly in the config."
    )


def _copy_metadata(src_root, dst_root) -> None:
    dst_root.attrs.update(dict(src_root.attrs))
    for key in src_root.keys():
        if key == "X":
            continue
        _log(f"Copying metadata group/key: {key}")
        zarr.copy(src_root[key], dst_root, name=key, if_exists="replace")


def _validate_csr_group(X, matrix_path: str) -> None:
    missing = [name for name in ("data", "indices", "indptr") if name not in X]
    if missing:
        raise KeyError(
            f"Source matrix {matrix_path!r} does not look like CSR zarr; "
            f"missing {missing}."
        )


def _copy_1d_array(src, dst, *, chunk_size: int, label: str) -> None:
    n = int(src.shape[0])
    done = 0
    while done < n:
        stop = min(done + int(chunk_size), n)
        dst[done:stop] = src[done:stop]
        done = stop
        if done == n or (done % (int(chunk_size) * 10) == 0):
            _log(f"  copied {label}: {done:,}/{n:,}")


def _build_cp10k_x(
    src_root,
    dst_root,
    *,
    source_matrix: str,
    umi: np.ndarray,
    scale_factor: float,
    row_block_size: int,
) -> dict[str, int | float | str]:
    X = _get_group(src_root, source_matrix)
    _validate_csr_group(X, source_matrix)

    data = X["data"]
    indices = X["indices"]
    indptr = X["indptr"]

    n_rows = len(indptr) - 1
    n_nnz = len(data)
    if len(umi) != n_rows:
        raise ValueError(f"UMI length mismatch: umi={len(umi)}, n_rows={n_rows}")
    if row_block_size <= 0:
        raise ValueError("row_block_size must be positive.")
    if scale_factor <= 0:
        raise ValueError("scale_factor must be positive.")

    _log(f"Creating CP10K CSR X: n_rows={n_rows:,}, nnz={n_nnz:,}")
    dst_X = dst_root.create_group("X", overwrite=True)
    dst_X.attrs.update(dict(X.attrs))

    chunk_nnz = int(getattr(data, "chunks", (10_000_000,))[0])
    chunk_indptr = min(n_rows + 1, 1_000_000)

    dst_data = dst_X.create_dataset(
        "data",
        shape=data.shape,
        chunks=(chunk_nnz,),
        dtype="float32",
        compressor=getattr(data, "compressor", None),
        overwrite=True,
    )
    dst_indices = dst_X.create_dataset(
        "indices",
        shape=indices.shape,
        chunks=getattr(indices, "chunks", (chunk_nnz,)),
        dtype=indices.dtype,
        compressor=getattr(indices, "compressor", None),
        overwrite=True,
    )
    dst_X.create_dataset(
        "indptr",
        data=indptr[:],
        chunks=(chunk_indptr,),
        dtype=indptr.dtype,
        compressor=getattr(indptr, "compressor", None),
        overwrite=True,
    )
    _copy_1d_array(indices, dst_indices, chunk_size=chunk_nnz, label="indices")

    indptr_np = indptr[:]
    umi_safe = np.asarray(umi, dtype=np.float64)
    bad_umi = int((~np.isfinite(umi_safe) | (umi_safe <= 0)).sum())
    umi_safe[~np.isfinite(umi_safe)] = 0.0
    umi_safe = np.maximum(umi_safe, 1.0)
    factors = scale_factor / umi_safe

    for r0 in range(0, n_rows, row_block_size):
        r1 = min(r0 + row_block_size, n_rows)
        start = int(indptr_np[r0])
        end = int(indptr_np[r1])

        if end > start:
            row_nnz = np.diff(indptr_np[r0 : r1 + 1])
            row_factors = np.repeat(factors[r0:r1], row_nnz).astype(np.float32)
            dst_data[start:end] = data[start:end].astype(np.float32) * row_factors

        if (r0 % (row_block_size * 20) == 0) or r1 == n_rows:
            _log(f"  rows {r0:,}-{r1:,} / {n_rows:,}")

    _log("Finished CP10K X.")
    return {
        "source_matrix": source_matrix,
        "n_rows": int(n_rows),
        "n_nnz": int(n_nnz),
        "scale_factor": float(scale_factor),
        "bad_or_nonpositive_umi_rows": int(bad_umi),
    }


def build_cp10k_zarr(cfg: BuildCp10kZarrConfig) -> dict:
    src_zarr = Path(cfg.src_zarr_path)
    out_zarr = Path(cfg.out_zarr_path)
    if not src_zarr.exists():
        raise FileNotFoundError(f"src_zarr_path does not exist: {src_zarr}")

    if out_zarr.exists():
        if not cfg.overwrite:
            raise FileExistsError(f"{out_zarr} exists. Set overwrite=true to replace it.")
        _log(f"Removing existing output zarr: {out_zarr}")
        shutil.rmtree(out_zarr)

    out_zarr.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = (
        Path(cfg.out_manifest_path)
        if cfg.out_manifest_path is not None
        else out_zarr.with_suffix(out_zarr.suffix + ".manifest.json")
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    _log("Reading obs/var")
    obs, var = _read_obs_var(str(src_zarr))
    umi_name, umi = _choose_umi(obs, cfg.umi_col)

    _log("Opening zarr stores")
    src_root = zarr.open(str(src_zarr), mode="r")
    dst_root = zarr.open(str(out_zarr), mode="w")

    _copy_metadata(src_root, dst_root)
    matrix_meta = _build_cp10k_x(
        src_root,
        dst_root,
        source_matrix=cfg.source_matrix,
        umi=umi,
        scale_factor=cfg.scale_factor,
        row_block_size=cfg.row_block_size,
    )

    manifest = {
        "purpose": "raw_count_to_cp10k_zarr",
        "config": asdict(cfg),
        "src_zarr_path": str(src_zarr),
        "out_zarr_path": str(out_zarr),
        "umi_col": umi_name,
        "n_obs": int(len(obs)),
        "n_vars": int(len(var)),
        "matrix": matrix_meta,
        "notes": [
            "Output X is sparse CSR CP10K, not log1p.",
            "Ordinal thresholds should be fitted downstream with make_ordinal_data_zarr.py.",
            "Reference scalers should be fitted downstream with build_celltype_control_scaler.py.",
        ],
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    _log(f"Wrote manifest: {manifest_path}")
    _log("done")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()
    cfg = load_config(args.config)
    build_cp10k_zarr(cfg)


if __name__ == "__main__":
    main()
