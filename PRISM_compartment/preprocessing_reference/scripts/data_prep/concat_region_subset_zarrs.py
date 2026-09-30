#!/usr/bin/env python
"""Concatenate region-level gene-union AnnData-Zarr stores.

This script is intentionally narrow: it concatenates DLPFC and MTG raw
gene-union subset zarrs that already share the same gene order. It copies CSR
matrix arrays in chunks, so the full data/indices arrays are never loaded into
memory at once.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zarr


@dataclass
class RegionInput:
    name: str
    zarr_path: str
    region_label: str


@dataclass
class ConcatRegionSubsetConfig:
    inputs: list[dict[str, str]]
    out_zarr_path: str
    out_gene_indices_path: str
    out_gene_names_json: str
    out_manifest_path: str
    data_chunk_size: int = 2_000_000
    indptr_chunk_size: int = 1_000_000
    copy_chunk_size: int = 20_000_000
    prefix_obs_index: bool = True
    overwrite: bool = False


def _log(msg: str) -> None:
    print(f"[concat-region-zarr] {msg}", flush=True)


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
        _log(f"{'DONE' if exc_type is None else 'FAIL'}: {self.name} ({dt:.1f}s)")
        return False


def _read_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError("Could not read AnnData-Zarr metadata.") from e


def _write_elem(root, key: str, value) -> None:
    try:
        from anndata.experimental import write_elem  # type: ignore
        write_elem(root, key, value)
    except Exception as e:
        raise RuntimeError(f"Could not write AnnData-Zarr element {key!r}.") from e


def load_config(path: str | Path) -> ConcatRegionSubsetConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return ConcatRegionSubsetConfig(**raw)


def _csr_shape(root) -> tuple[int, int]:
    shape = root["X"].attrs.get("shape")
    if shape is None:
        raise ValueError(f"Missing X.attrs['shape'] in {root.store}")
    return int(shape[0]), int(shape[1])


def _copy_1d(src, dst, *, dst_offset: int, chunk: int, label: str) -> None:
    n = int(src.shape[0])
    done = 0
    while done < n:
        stop = min(done + int(chunk), n)
        dst[dst_offset + done: dst_offset + stop] = src[done:stop]
        done = stop
        _log(f"  copied {label}: {done:,}/{n:,}")


def concat_region_subset_zarrs(cfg: ConcatRegionSubsetConfig) -> dict[str, Any]:
    inputs = [RegionInput(**x) for x in cfg.inputs]
    if len(inputs) < 2:
        raise ValueError("At least two inputs are required.")

    out_path = Path(cfg.out_zarr_path)
    if out_path.exists():
        if not cfg.overwrite:
            raise FileExistsError(f"{out_path} already exists. Set overwrite=true to replace it.")
        _log(f"removing existing output: {out_path}")
        shutil.rmtree(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    roots = []
    obs_frames = []
    shapes: list[tuple[int, int]] = []
    nnzs: list[int] = []

    with Timer("open inputs and metadata"):
        for inp in inputs:
            root = zarr.open(inp.zarr_path, mode="r")
            if "X" not in root or not all(k in root["X"] for k in ("data", "indices", "indptr")):
                raise ValueError(f"Input {inp.zarr_path} does not contain CSR X.")
            if "obs" not in root or "var" not in root:
                raise ValueError(f"Input {inp.zarr_path} must contain obs and var.")

            obs = _read_elem(root["obs"])
            if not isinstance(obs, pd.DataFrame):
                raise TypeError(f"obs in {inp.zarr_path} is not a DataFrame.")
            obs = obs.copy()
            obs.index = obs.index.astype(str)
            obs["concat_region"] = inp.region_label
            obs["concat_source_name"] = inp.name
            obs["concat_source_obs_index"] = obs.index.astype(str)
            if cfg.prefix_obs_index:
                obs.index = pd.Index([f"{inp.region_label}::{x}" for x in obs.index], name=obs.index.name)

            shape = _csr_shape(root)
            nnz = int(root["X"]["data"].shape[0])
            if len(obs) != shape[0]:
                raise ValueError(f"obs length mismatch for {inp.name}: {len(obs)} vs {shape[0]}")

            roots.append(root)
            obs_frames.append(obs)
            shapes.append(shape)
            nnzs.append(nnz)
            _log(f"input {inp.name}: cells={shape[0]:,}, genes={shape[1]:,}, nnz={nnz:,}")

        n_genes = shapes[0][1]
        for inp, shape in zip(inputs, shapes):
            if shape[1] != n_genes:
                raise ValueError(f"Gene count mismatch for {inp.name}: {shape[1]} vs {n_genes}")

        var0 = _read_elem(roots[0]["var"]).copy()
        for inp, root in zip(inputs[1:], roots[1:]):
            var_i = _read_elem(root["var"])
            if not var0.index.astype(str).equals(var_i.index.astype(str)):
                raise ValueError(f"var gene order mismatch for {inp.name}")

        obs_out = pd.concat(obs_frames, axis=0, join="outer", copy=False)
        if not obs_out.index.is_unique:
            dup = obs_out.index[obs_out.index.duplicated()].unique().tolist()[:5]
            raise ValueError(f"Concatenated obs index still has duplicates: {dup}")

    n_cells_total = int(sum(x[0] for x in shapes))
    nnz_total = int(sum(nnzs))

    with Timer("create output zarr and metadata"):
        out = zarr.open(str(out_path), mode="w")
        out.attrs["encoding-type"] = "anndata"
        out.attrs["encoding-version"] = "0.1.0"
        out.attrs["shape"] = [n_cells_total, n_genes]
        _write_elem(out, "obs", obs_out)
        _write_elem(out, "var", var0)

        comp = getattr(roots[0]["X"]["data"], "compressor", None)
        x = out.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs["shape"] = [n_cells_total, n_genes]
        data_z = x.create_dataset(
            "data",
            shape=(nnz_total,),
            chunks=(int(cfg.data_chunk_size),),
            dtype=roots[0]["X"]["data"].dtype,
            compressor=comp,
        )
        indices_z = x.create_dataset(
            "indices",
            shape=(nnz_total,),
            chunks=(int(cfg.data_chunk_size),),
            dtype=np.int32,
            compressor=comp,
        )
        indptr_z = x.create_dataset(
            "indptr",
            shape=(n_cells_total + 1,),
            chunks=(min(int(cfg.indptr_chunk_size), n_cells_total + 1),),
            dtype=np.int64,
            compressor=comp,
        )

    with Timer("copy CSR arrays"):
        cell_offset = 0
        nnz_offset = 0
        indptr_z[0] = 0
        for inp, root, shape, nnz in zip(inputs, roots, shapes, nnzs):
            _log(f"copy matrix for {inp.name}")
            _copy_1d(
                root["X"]["data"],
                data_z,
                dst_offset=nnz_offset,
                chunk=int(cfg.copy_chunk_size),
                label=f"{inp.name}/data",
            )
            _copy_1d(
                root["X"]["indices"],
                indices_z,
                dst_offset=nnz_offset,
                chunk=int(cfg.copy_chunk_size),
                label=f"{inp.name}/indices",
            )

            in_indptr = root["X"]["indptr"]
            rows = int(shape[0])
            indptr_z[cell_offset + 1: cell_offset + rows + 1] = (
                np.asarray(in_indptr[1:], dtype=np.int64) + np.int64(nnz_offset)
            )
            cell_offset += rows
            nnz_offset += nnz
            _log(f"  finished {inp.name}: cell_offset={cell_offset:,}, nnz_offset={nnz_offset:,}")

    out_gene_indices_path = Path(cfg.out_gene_indices_path)
    out_gene_names_json = Path(cfg.out_gene_names_json)
    out_manifest_path = Path(cfg.out_manifest_path)
    with Timer("write side files"):
        out_gene_indices_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(out_gene_indices_path, np.arange(n_genes, dtype=np.int64))
        with open(out_gene_names_json, "w", encoding="utf-8") as f:
            json.dump([str(x) for x in var0.index.astype(str).tolist()], f, ensure_ascii=False, indent=2)

        manifest = {
            "config": asdict(cfg),
            "build_date": datetime.now(timezone.utc).isoformat(),
            "inputs": [
                {
                    "name": inp.name,
                    "region_label": inp.region_label,
                    "zarr_path": inp.zarr_path,
                    "n_obs": shape[0],
                    "n_genes": shape[1],
                    "nnz": nnz,
                }
                for inp, shape, nnz in zip(inputs, shapes, nnzs)
            ],
            "out_zarr_path": str(out_path),
            "out_matrix_path": "X",
            "n_obs": n_cells_total,
            "n_final_genes": n_genes,
            "nnz_out": nnz_total,
            "density": float(nnz_total / max(n_cells_total * n_genes, 1)),
            "obs_columns": list(map(str, obs_out.columns)),
            "gene_indices_path_for_subset": str(out_gene_indices_path),
            "gene_names_json_for_subset": str(out_gene_names_json),
            "downstream_config_note": {
                "zarr_path": str(out_path),
                "matrix_path": "X",
                "gene_indices_path": str(out_gene_indices_path),
            },
        }
        with open(out_manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)

    _log("done")
    _log(f"combined zarr: {out_path}")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Concatenate region gene-union subset zarrs.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)
    concat_region_subset_zarrs(cfg)


if __name__ == "__main__":
    main()
