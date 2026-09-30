from __future__ import annotations

"""
Build module-first artifacts from the original full-gene zarr.

What this script produces
-------------------------
1. ENSG <-> gene symbol mapping table
2. GeneModuleRegistry JSON (ENSG as internal key)
3. LieAction decoder module masks
4. Module-level encoder input h5ad [cells x modules]
5. Gene-union metadata files for the future gene-level decoder branch

This version addresses four issues explicitly:
- symbol->ENSG mapping policy is configurable (default: unambiguous)
- raw-count assumption is explicit and checked
- residual inclusion in decoder gene-union is configurable
- dead split_key API is removed
"""

from pathlib import Path
import argparse
from typing import Dict, Iterable, List, Literal, Optional, Sequence
import json
import warnings
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import zarr

try:
    from kmlee_bam.gene_module.gene_module_registry import (
        GeneModuleRegistry,
        ModuleCollectionSpec,
        build_symbol_to_gene_ids_map,
    )
except ImportError:  # direct-file execution fallback
    from gene_module_registry import (
        GeneModuleRegistry,
        ModuleCollectionSpec,
        build_symbol_to_gene_ids_map,
    )

MappingStrategy = Literal["all", "first", "unambiguous"]


def _log(msg: str) -> None:
    print(f"[module-build] {msg}", flush=True)


class StageTimer:
    def __init__(self, name: str) -> None:
        self.name = name
        self.t0 = None

    def __enter__(self):
        self.t0 = time.time()
        _log(f"START: {self.name}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.time() - self.t0
        if exc_type is None:
            _log(f"DONE : {self.name} ({dt:.1f}s)")
        else:
            _log(f"FAIL : {self.name} after {dt:.1f}s")
        return False


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


def _read_obs_var_names_from_zarr(zarr_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.Index, bool]:
    root = zarr.open(str(zarr_path), mode="r")
    obs = _read_zarr_dataframe(root, "obs")
    if "raw" in root and "var" in root["raw"]:
        var = _read_zarr_dataframe(root, "raw/var")
        use_raw = True
    else:
        var = _read_zarr_dataframe(root, "var")
        use_raw = False
    var_names = pd.Index(var.index.astype(str))
    return obs, var, var_names, use_raw

# ---------------------------------------------------------------------
# Mapping helpers
# ---------------------------------------------------------------------
def build_ensembl_symbol_table_from_zarr(
    zarr_path: str | Path,
    *,
    symbol_col_candidates: tuple[str, ...] = (
        "gene_symbol",
        "symbol",
        "feature_name",
        "gene_name",
        "external_gene_name",
    ),
) -> pd.DataFrame:
    _, var, var_names, _ = _read_obs_var_names_from_zarr(zarr_path)

    ensembl_id = var_names.str.replace(r"\.\d+$", "", regex=True).str.upper()

    symbol_col = None
    for c in symbol_col_candidates:
        if c in var.columns:
            symbol_col = c
            break
    if symbol_col is None:
        raise ValueError(f"Could not find a gene symbol column in var. Tried: {symbol_col_candidates}")

    gene_symbol = var[symbol_col].astype(str).fillna("").str.strip()

    table = pd.DataFrame(
        {
            "ensembl_id": ensembl_id.to_numpy(),
            "gene_symbol": gene_symbol.to_numpy(),
        }
    )
    table["display_name"] = table["gene_symbol"]
    table.loc[table["display_name"].eq(""), "display_name"] = table["ensembl_id"]
    table = table.drop_duplicates(subset=["ensembl_id"], keep="first").reset_index(drop=True)
    return table


# ---------------------------------------------------------------------
# GMT / Reactome helpers
# ---------------------------------------------------------------------
def _read_gmt(path: str | Path) -> List[tuple[str, List[str]]]:
    out: List[tuple[str, List[str]]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = [g.strip() for g in parts[2:] if g.strip()]
            out.append((name, genes))
    return out


def normalize_set_name(name: str) -> str:
    x = name.upper()
    x = x.replace("REACTOME_", "")
    x = x.replace("HALLMARK_", "")
    x = x.replace("__", "_")
    return x


REACTOME_KEYWORDS = [
    # immune / glia
    "IMMUNE", "INNATE", "ADAPTIVE", "INTERFERON", "INTERLEUKIN", "CYTOKINE",
    "COMPLEMENT", "ANTIGEN", "MICROGLIA", "TLR", "TOLL", "NFKB",
    # AD / neurodegeneration relevant
    "AMYLOID", "APP", "TAU", "MAPT", "APOPTOSIS", "AUTOPHAGY", "LYSOSOME",
    "ENDOCYTOSIS", "VESICLE", "SYNAPT", "NEURON", "NEURONAL", "AXON",
    # metabolism / mitochondria / stress
    "MITOCHON", "RESPIRATORY", "OXIDATIVE", "ELECTRON TRANSPORT",
    "CHOLESTEROL", "LIPID", "GLYCOLYSIS", "ROS", "PROTEIN FOLDING", "UPR",
    # signaling
    "PI3K", "AKT", "MTOR", "NOTCH", "WNT", "TGF", "MAPK", "JAK", "STAT",
]


def keyword_score(name: str, keywords: list[str]) -> int:
    norm = normalize_set_name(name)
    return sum(int(kw in norm) for kw in keywords)


def build_reactome_candidate_table(
    reactome_gmt: str | Path,
    symbol_to_gene_ids: Dict[str, Sequence[str]],
    *,
    min_size: int = 8,
    max_size: int = 250,
    keywords: Optional[list[str]] = None,
) -> pd.DataFrame:
    if keywords is None:
        keywords = REACTOME_KEYWORDS

    rows = []
    for name, genes in _read_gmt(reactome_gmt):
        n = 0
        seen = set()
        for g in genes:
            hits = symbol_to_gene_ids.get(g.upper(), [])
            for h in hits:
                if h not in seen:
                    seen.add(h)
                    n += 1
        if n < min_size or n > max_size:
            continue

        score = keyword_score(name, keywords)
        if score <= 0:
            continue

        rows.append(
            {
                "name": name,
                "n_genes_after_map": n,
                "keyword_score": score,
                "norm_name": normalize_set_name(name),
            }
        )

    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df
    return df.sort_values(["keyword_score", "n_genes_after_map", "name"], ascending=[False, True, True]).reset_index(drop=True)


# ---------------------------------------------------------------------
# Registry builder
# ---------------------------------------------------------------------
def build_registry_from_full_zarr(
    *,
    zarr_path: str | Path,
    hallmark_gmt: str | Path,
    reactome_gmt: Optional[str | Path] = None,
    singleton_genes: Optional[Sequence[str]] = None,
    reactome_top_k: int = 25,
    reactome_manual_keep: Optional[Sequence[str]] = None,
    reactome_manual_drop: Optional[Sequence[str]] = None,
    min_module_size: int = 8,
    max_module_size: int = 250,
    prune_redundancy_jaccard: float = 0.92,
    symbol_map_strategy: MappingStrategy = "unambiguous",
    out_dir: Optional[str | Path] = None,
) -> tuple[GeneModuleRegistry, pd.DataFrame, pd.DataFrame, Dict[str, List[str]]]:
    _log("Preparing registry build from full zarr")

    with StageTimer("build ENSG <-> symbol mapping table"):
        table = build_ensembl_symbol_table_from_zarr(zarr_path)
        gene_names_ensg = table["ensembl_id"].tolist()
        gene_symbols = table["gene_symbol"].tolist()
        symbol_to_gene_ids = build_symbol_to_gene_ids_map(
            gene_ids=gene_names_ensg,
            gene_symbols=gene_symbols,
            strategy=symbol_map_strategy,
        )
        _log(f"#genes in zarr = {len(gene_names_ensg):,}")
        _log(f"#symbols mapped = {(table['gene_symbol'] != '').sum():,}")
        _log(f"#symbol->gene ids entries = {len(symbol_to_gene_ids):,}")


    specs: List[ModuleCollectionSpec] = [
        ModuleCollectionSpec(
            path=str(hallmark_gmt),
            source="hallmark",
            token_type="regular",
            min_size=min_module_size,
            max_size=max_module_size,
            prefix="HALLMARK::",
        )
    ]

    reactome_candidate_df = pd.DataFrame()
    if reactome_gmt is not None:
        with StageTimer("build Reactome candidate table"):
            reactome_candidate_df = build_reactome_candidate_table(
                reactome_gmt,
                symbol_to_gene_ids,
                min_size=min_module_size,
                max_size=max_module_size,
            )
            _log(f"#Reactome candidates = {len(reactome_candidate_df):,}")

        auto_keep = reactome_candidate_df["name"].head(reactome_top_k).tolist() if len(reactome_candidate_df) > 0 else []
        selected = sorted((set(auto_keep) | set(reactome_manual_keep or [])) - set(reactome_manual_drop or []))
        _log(f"#Reactome selected = {len(selected):,}")

        specs.append(
            ModuleCollectionSpec(
                path=str(reactome_gmt),
                source="reactome",
                token_type="regular",
                min_size=min_module_size,
                max_size=max_module_size,
                prefix="REACTOME::",
                keep_names=selected if len(selected) > 0 else None,
            )
        )

    with StageTimer("construct GeneModuleRegistry"):
        registry = GeneModuleRegistry.from_collection_specs(
            gene_names=gene_names_ensg,
            collection_specs=specs,
            singleton_genes=singleton_genes,
            add_residual_module=True,
            residual_name="RESIDUAL_LEFTOVER",
            prune_redundancy_jaccard=prune_redundancy_jaccard,
            source_priority={"hallmark": 0, "reactome": 1},
            symbol_to_gene_ids=symbol_to_gene_ids,
            allow_direct_gene_name_match=True,
        )
        _log(f"#modules total = {registry.n_modules:,}")
        _log(f"#regular modules = {len(registry.regular_module_ids):,}")
        _log(f"#residual modules = {len(registry.residual_module_ids):,}")
        _log(f"#singleton modules = {len(registry.singleton_module_ids):,}")


    if out_dir is not None:
        with StageTimer("save registry artifacts"):
            out_dir = Path(out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            _log(f"Saving registry artifacts to: {out_dir}")
            table.to_csv(out_dir / "ensembl_symbol_table.csv", index=False)
            with open(out_dir / "symbol_to_gene_ids.json", "w", encoding="utf-8") as f:
                json.dump(symbol_to_gene_ids, f, ensure_ascii=False)
            registry.save(out_dir / "gene_module_registry.json")
            pd.DataFrame(
                {
                    "module_name": registry.module_names,
                    "source": registry.module_sources,
                    "token_type": registry.module_token_types,
                    "size": [m.size for m in registry.modules],
                }
            ).to_csv(out_dir / "module_table.csv", index=False)
            if len(reactome_candidate_df) > 0:
                reactome_candidate_df.to_csv(out_dir / "reactome_candidates.csv", index=False)
            torch.save(
                registry.make_decoder_module_masks(32, strategy="proportional"),
                out_dir / "lie_decoder_module_masks.pt",
            )

    return registry, table, reactome_candidate_df, symbol_to_gene_ids


def _build_decoder_gene_union(
    registry: GeneModuleRegistry,
    *,
    include_regular: bool = True,
    include_residual: bool = False,
    include_singletons: bool = True,
) -> List[str]:
    ids: List[int] = []
    if include_regular:
        ids.extend(registry.regular_module_ids)
    if include_residual:
        ids.extend(registry.residual_module_ids)
    if include_singletons:
        ids.extend(registry.singleton_module_ids)

    ids = sorted(set(ids))
    if len(ids) == 0:
        return []

    union_mask = registry.membership_binary[ids].sum(dim=0) > 0
    gene_union_idx = union_mask.nonzero(as_tuple=True)[0].tolist()
    return [registry.gene_names[i] for i in gene_union_idx]


# ---------------------------------------------------------------------
# Chunked zarr -> module matrix aggregation
# ---------------------------------------------------------------------
def _choose_matrix_and_var(zarr_path: str | Path):
    obs, _, var_names, use_raw = _read_obs_var_names_from_zarr(zarr_path)
    return obs, var_names, use_raw


def _check_raw_count_assumption(data_sample: np.ndarray) -> None:
    if (data_sample < 0).any():
        raise ValueError("Negative values found in expression matrix; this does not look like raw counts.")
    if np.any(np.abs(data_sample - np.round(data_sample)) > 1e-6):
        warnings.warn(
            "Expression values are non-integer. Ensure this matrix truly represents raw counts before log1p.",
            RuntimeWarning,
        )


def aggregate_full_zarr_to_module_h5ad(
    *,
    zarr_path: str | Path,
    registry: GeneModuleRegistry,
    out_h5ad_path: str | Path,
    chunk_rows: int = 4096,
    assume_raw_counts: bool = True,
    print_every: int = 25,
    output_dtype: np.dtype = np.float32,
) -> None:
    """
    Build a module-level h5ad [cells x modules] from the original full-gene zarr.

    The matrix X in the output h5ad is continuous module values computed as
    weighted-mean expression over registry genes.
    """
    _log("Starting full-zarr -> module-h5ad aggregation")
    _log(f"zarr_path = {zarr_path}")
    _log(f"out_h5ad_path = {out_h5ad_path}")
    _log(f"chunk_rows = {chunk_rows}")
    _log(f"assume_raw_counts = {assume_raw_counts}")

    obs, var_names, use_raw = _choose_matrix_and_var(zarr_path)
    _log(f"use_raw = {use_raw}")
    _log(f"#cells = {len(obs):,}")
    _log(f"#total zarr genes = {len(var_names):,}")
    _log(f"#registry genes = {len(registry.gene_names):,}")
    _log(f"#modules = {registry.n_modules:,}")

    z = zarr.open(str(zarr_path), mode="r")
    xgrp = z["raw"]["X"] if use_raw else z["X"]

    var_names = var_names.str.replace(r"\.\d+$", "", regex=True).str.upper()
    var_to_idx = {str(v): i for i, v in enumerate(var_names.tolist())}

    col_idx = []
    for g in registry.gene_names:
        key = g.upper()
        if key not in var_to_idx:
            raise KeyError(f"Registry gene '{g}' not found in zarr var_names.")
        col_idx.append(var_to_idx[key])
    col_idx = np.asarray(col_idx, dtype=np.int64)

    n_obs = len(obs)
    n_vars_total = len(var_names)
    col_idx_is_identity = len(col_idx) == n_vars_total and np.array_equal(col_idx, np.arange(n_vars_total))
    if col_idx_is_identity:
        _log("registry gene order covers all zarr genes; skipping per-chunk column slicing")
    W_T = registry.membership_weight.t().cpu().numpy().astype(output_dtype)  # [G, M]

    chunks: List[np.ndarray] = []

    data_arr = xgrp["data"]
    indices_arr = xgrp["indices"]
    indptr_arr = xgrp["indptr"]

    # raw-count sanity check on a small slice
    if assume_raw_counts:
        _log("Checking raw-count assumption on a small sample")
        sample_n = min(2048, len(data_arr))
        if sample_n > 0:
            sample = np.asarray(data_arr[:sample_n], dtype=np.float32)
            _check_raw_count_assumption(sample)

    n_chunks = (n_obs + chunk_rows - 1) // chunk_rows
    for chunk_id, start in enumerate(range(0, n_obs, chunk_rows), start=1):
        stop = min(n_obs, start + chunk_rows)
        should_print = print_every > 0 and (chunk_id == 1 or chunk_id % print_every == 0 or chunk_id == n_chunks)
        if should_print:
            print(f"[module-agg] START chunk {chunk_id:04d}/{n_chunks:04d} rows {start:,}:{stop:,}", flush=True)

        p0 = int(indptr_arr[start])
        p1 = int(indptr_arr[stop])

        data = np.asarray(data_arr[p0:p1], dtype=output_dtype)
        indices = np.asarray(indices_arr[p0:p1], dtype=np.int64)
        indptr = np.asarray(indptr_arr[start : stop + 1], dtype=np.int64) - p0

        X_chunk = sp.csr_matrix((data, indices, indptr), shape=(stop - start, n_vars_total), dtype=output_dtype)
        if col_idx_is_identity:
            X_sel = X_chunk
        else:
            X_sel = X_chunk[:, col_idx].tocsr(copy=True)

        if assume_raw_counts:
            X_sel.data = np.log1p(X_sel.data)

        mod_chunk = (X_sel @ W_T).astype(output_dtype)   # [B_chunk, M]
        chunks.append(mod_chunk)

        if should_print:
            print(f"[module-agg] DONE  chunk {chunk_id:04d}/{n_chunks:04d} rows {start:,}:{stop:,}", flush=True)

    _log("Combining chunk outputs into one module matrix")
    X_mod = np.vstack(chunks)

    var = pd.DataFrame(
        {
            "module_name": registry.module_names,
            "source": registry.module_sources,
            "token_type": registry.module_token_types,
            "size": [m.size for m in registry.modules],
        },
        index=pd.Index(registry.module_names, name="module_name"),
    )

    adata_out = ad.AnnData(X=X_mod, obs=obs.copy(), var=var)
    adata_out.uns["registry_gene_names"] = registry.gene_names
    adata_out.uns["module_type_ids"] = registry.module_type_ids.tolist()
    adata_out.write_h5ad(str(out_h5ad_path))
    
    _log("Module-level h5ad build completed")


def save_decoder_gene_union_artifacts(
    registry: GeneModuleRegistry,
    *,
    out_dir: str | Path,
    include_regular: bool = True,
    include_residual: bool = False,
    include_singletons: bool = True,
) -> None:
    _log("Saving decoder gene-union artifacts")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_union_ensg = _build_decoder_gene_union(
        registry,
        include_regular=include_regular,
        include_residual=include_residual,
        include_singletons=include_singletons,
    )
    with open(out_dir / "decoder_gene_union_ensembl.json", "w", encoding="utf-8") as f:
        json.dump(gene_union_ensg, f, ensure_ascii=False)

    policy = {
        "include_regular": include_regular,
        "include_residual": include_residual,
        "include_singletons": include_singletons,
        "n_decoder_genes": len(gene_union_ensg),
    }
    with open(out_dir / "decoder_gene_union_policy.json", "w", encoding="utf-8") as f:
        json.dump(policy, f, ensure_ascii=False, indent=2)
    
    _log("Saving decoder gene-union artifacts")


def load_build_config(config_path: str | Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build module-first artifacts from the original full-gene zarr."
    )
    p.add_argument("--config", type=str, required=True, help="Path to build config JSON.")
    return p


def main() -> None:
    args = build_argparser().parse_args()
    cfg = load_build_config(args.config)
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    _log(f"Loaded config: {args.config}")
    _log(f"Output directory: {out_dir}")

    
    with StageTimer("build registry from full zarr"):
        registry, table, reactome_candidate_df, symbol_to_gene_ids = build_registry_from_full_zarr(
            zarr_path=cfg["zarr_path"],
            hallmark_gmt=cfg["hallmark_gmt"],
            reactome_gmt=cfg.get("reactome_gmt"),
            singleton_genes=cfg.get("singleton_genes"),
            reactome_top_k=cfg.get("reactome_top_k", 25),
            reactome_manual_keep=cfg.get("reactome_manual_keep"),
            reactome_manual_drop=cfg.get("reactome_manual_drop"),
            min_module_size=cfg.get("min_module_size", 8),
            max_module_size=cfg.get("max_module_size", 250),
            prune_redundancy_jaccard=cfg.get("prune_redundancy_jaccard", 0.92),
            symbol_map_strategy=cfg.get("symbol_map_strategy", "unambiguous"),
            out_dir=out_dir,
        )

    with StageTimer("save decoder gene-union"):
        save_decoder_gene_union_artifacts(
            registry,
            out_dir=out_dir,
            include_regular=cfg.get("include_regular_in_decoder", True),
            include_residual=cfg.get("include_residual_in_decoder", False),
            include_singletons=cfg.get("include_singletons_in_decoder", True),
        )

    with StageTimer("aggregate full zarr to module h5ad"):
        aggregate_full_zarr_to_module_h5ad(
            zarr_path=cfg["zarr_path"],
            registry=registry,
            out_h5ad_path=cfg.get(
                "module_h5ad_path",
                str(out_dir / "module_encoder_input.h5ad"),
            ),
            chunk_rows=cfg.get("chunk_rows", 4096),
            assume_raw_counts=cfg.get("assume_raw_counts", True),
            print_every=cfg.get("print_every", 25),
        )
  
    _log(f"ALL DONE. saved artifacts to: {out_dir}")


if __name__ == "__main__":
    main()
