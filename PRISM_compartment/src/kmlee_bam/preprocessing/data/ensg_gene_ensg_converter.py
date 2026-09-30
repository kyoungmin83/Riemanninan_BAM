from __future__ import annotations

from pathlib import Path
from typing import Optional
import json

import anndata as ad
import numpy as np
import pandas as pd


def load_spec_gene_ids(npz_path: str) -> list[str]:
    obj = np.load(npz_path, allow_pickle=True)
    gene_ids = [str(x) for x in obj["gene_names"].tolist()]
    return gene_ids


def build_ensembl_symbol_table_from_h5ad(
    h5ad_path: str,
    *,
    ensembl_col: Optional[str] = None,
    symbol_col_candidates: tuple[str, ...] = (
        "gene_symbol",
        "symbol",
        "feature_name",
        "gene_name",
        "external_gene_name",
    ),
) -> pd.DataFrame:
    """
    Build ENSG <-> symbol table from adata.var.

    Priority:
    1) If var_names are ENSG-like, use var_names as ensembl_id.
    2) Else use `ensembl_col` if provided.
    3) Symbol column is searched from candidates.
    """
    adata = ad.read_h5ad(h5ad_path)
    var = adata.var.copy()

    var_names = pd.Index(var.index.astype(str))
    looks_like_ensg = var_names.str.startswith("ENSG").mean() > 0.8

    if looks_like_ensg:
        ensembl_id = var_names
    elif ensembl_col is not None and ensembl_col in var.columns:
        ensembl_id = var[ensembl_col].astype(str)
    else:
        raise ValueError(
            "Could not determine Ensembl IDs. "
            "Either var_names must be ENSG-like or ensembl_col must be provided."
        )

    symbol_col = None
    for c in symbol_col_candidates:
        if c in var.columns:
            symbol_col = c
            break

    if symbol_col is None:
        raise ValueError(
            f"Could not find a symbol column in adata.var. "
            f"Tried: {symbol_col_candidates}"
        )

    symbol = var[symbol_col].astype(str).fillna("")

    table = pd.DataFrame(
        {
            "ensembl_id": ensembl_id.to_numpy(),
            "gene_symbol": symbol.to_numpy(),
        }
    )

    # clean
    table["ensembl_id"] = table["ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
    table["gene_symbol"] = table["gene_symbol"].str.strip()

    # drop empty ENSG
    table = table[table["ensembl_id"] != ""].copy()

    # keep first non-empty symbol per ENSG
    table["symbol_is_empty"] = table["gene_symbol"].eq("")
    table = table.sort_values(["ensembl_id", "symbol_is_empty"])
    table = table.drop_duplicates(subset=["ensembl_id"], keep="first")
    table = table.drop(columns=["symbol_is_empty"]).reset_index(drop=True)

    return table


def restrict_table_to_spec_genes(
    spec_gene_ids: list[str],
    ensg_symbol_table: pd.DataFrame,
) -> pd.DataFrame:
    spec_df = pd.DataFrame({"ensembl_id": [g.split(".")[0] for g in spec_gene_ids]})
    out = spec_df.merge(ensg_symbol_table, on="ensembl_id", how="left")

    # fallback: if symbol missing, keep ENSG itself as display label
    out["display_name"] = out["gene_symbol"].fillna("")
    out.loc[out["display_name"].eq(""), "display_name"] = out["ensembl_id"]

    return out


def build_symbol_to_ensembl_map(
    table: pd.DataFrame,
    *,
    uppercase_symbol: bool = True,
) -> dict[str, list[str]]:
    """
    Symbol -> possibly multiple ENSG IDs.
    Keep list because symbol collisions happen.
    """
    tmp = table.copy()
    sym = tmp["gene_symbol"].fillna("").astype(str).str.strip()
    if uppercase_symbol:
        sym = sym.str.upper()
    tmp["gene_symbol_norm"] = sym
    tmp = tmp[tmp["gene_symbol_norm"] != ""].copy()

    grouped = tmp.groupby("gene_symbol_norm")["ensembl_id"].apply(list)
    return grouped.to_dict()


def save_mapping_artifacts(
    table: pd.DataFrame,
    out_dir: str,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    table.to_csv(out / "ensembl_symbol_table.csv", index=False)

    ensg_to_symbol = dict(zip(table["ensembl_id"], table["display_name"]))
    with open(out / "ensembl_to_symbol.json", "w", encoding="utf-8") as f:
        json.dump(ensg_to_symbol, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    
    npz_path = "/home/kmlee/project_local/kmlee_bam/outputs/SEA_AD/dlpfc_rechunked_sample_3000hvg_ordinal_spec.npz"
    h5ad_path = "/home/kmlee/project_local/kmlee_bam/outputs/SEA_AD/dlpfc_rechunked_sample_3000hvg.h5ad"
    out_dir = "/home/kmlee/project_local/kmlee_bam/outputs/SEA_AD/module_mapping_outputs"

    spec_gene_ids = load_spec_gene_ids(npz_path)
    table0 = build_ensembl_symbol_table_from_h5ad(h5ad_path)
    table = restrict_table_to_spec_genes(spec_gene_ids, table0)

    save_mapping_artifacts(table, out_dir)

    sym2ensg = build_symbol_to_ensembl_map(table)
    print(table.head())
    print(f"#spec genes: {len(spec_gene_ids)}")
    print(f"#mapped symbols: {(table['gene_symbol'].fillna('') != '').sum()}")
    print(f"#unique symbols: {table['gene_symbol'].fillna('').nunique()}")