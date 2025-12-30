# %%
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import anndata as ad


#%%
@dataclass
class QCRules:
    # common rules
    min_genes: int = 400
    max_genes: int = 8000
    min_counts: int = 2500
    max_mt: float = 10.0

    # optional for Smart-seq
    max_ercc: Optional[float] = None

    # override rules per assay
    per_assay: Optional[Dict[str, Dict[str, float]]] = None

#%%
def _get_value(rule_dict: Dict[str, float], key: str, default):
    return rule_dict.get(key, default)


def build_keep_mask(
    adata: ad.AnnData,
    rules: QCRules,
    assay_key: str = "assay",
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Returns:
      keep_mask: boolean array of length n_obs
      reasons: per-cell boolean dataframe for each rule (True = pass)
    """
    obs = adata.obs

    required = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    for c in required:
        if c not in obs.columns:
            raise ValueError(f"Missing required QC column in adata.obs: '{c}'")

    # pass/fail
    pass_min_genes = (obs["n_genes_by_counts"] >= rules.min_genes).astype(bool)
    pass_max_genes = (obs["n_genes_by_counts"] <= rules.max_genes).astype(bool)
    pass_min_counts = (obs["total_counts"] >= rules.min_counts).astype(bool)
    pass_max_mt = (obs["pct_counts_mt"] <= rules.max_mt).astype(bool)

    
    pass_max_counts = pd.Series(True, index=obs.index).astype(bool)


    has_ercc = "pct_counts_ercc" in obs.columns
    if rules.max_ercc is not None and has_ercc:
        pass_max_ercc = (obs["pct_counts_ercc"] <= rules.max_ercc).astype(bool)
    else:
        pass_max_ercc = pd.Series(True, index=obs.index, dtype=bool)

    # applying override per assay
    if rules.per_assay is not None:
        if assay_key not in obs.columns:
            raise ValueError(f"assay_key='{assay_key}' not in adata.obs.columns")

        assays = obs[assay_key].astype(str)

        for assay_name, override in rules.per_assay.items():
            m = assays == assay_name
            if not np.any(m):
                continue

            # override keys: min_genes, max_genes, min_counts, max_counts, max_mt, max_ercc
            min_genes = _get_value(override, "min_genes", rules.min_genes)
            max_genes = _get_value(override, "max_genes", rules.max_genes)
            min_counts = _get_value(override, "min_counts", rules.min_counts)
            max_counts = _get_value(override, "max_counts", None)
            max_mt = _get_value(override, "max_mt", rules.max_mt)
            max_ercc = _get_value(override, "max_ercc", rules.max_ercc)

            pass_min_genes[m] = (obs.loc[m, "n_genes_by_counts"] >= min_genes).values
            pass_max_genes[m] = (obs.loc[m, "n_genes_by_counts"] <= max_genes).values
            pass_min_counts[m] = (obs.loc[m, "total_counts"] >= min_counts).values
            pass_max_mt[m] = (obs.loc[m, "pct_counts_mt"] <= max_mt).values

            if max_counts is not None:
                pass_max_counts[m] = (obs.loc[m, "total_counts"] <= max_counts).values

            if max_ercc is not None and has_ercc:
                pass_max_ercc[m] = (obs.loc[m, "pct_counts_ercc"] <= max_ercc).values

    reasons = pd.DataFrame(
        {
            "min_genes": pass_min_genes.values,
            "max_genes": pass_max_genes.values,
            "min_counts": pass_min_counts.values,
            "max_counts": pass_max_counts.values,
            "max_mt": pass_max_mt.values,
            "max_ercc": pass_max_ercc.values,
        },
        index=obs.index,
    )

    keep = reasons.all(axis=1).values
    return keep, reasons


def filter_adata(
    adata: ad.AnnData,
    rules: QCRules,
    assay_key: str = "assay",
    copy: bool = True,
) -> Tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      adata_qc: filtered AnnData
      retention_by_assay: per-assay before/after counts
      drop_reasons: counts of failures per rule (overall and by assay)
    """
    keep, reasons = build_keep_mask(adata, rules, assay_key=assay_key)

    # retention table
    obs = adata.obs
    if assay_key in obs.columns:
        assay = obs[assay_key].astype(str)
        ret = (
            pd.DataFrame({"assay": assay, "keep": keep})
            .groupby("assay", dropna=False)["keep"]
            .agg(n_before="size", n_after="sum")
            .reset_index()
        )
        ret["retention_rate"] = ret["n_after"] / ret["n_before"]
    else:
        ret = pd.DataFrame(
            {"assay": ["ALL"], "n_before": [adata.n_obs], "n_after": [int(keep.sum())], "retention_rate": [keep.mean()]}
        )

    # drop reason summary (overall)
    fail = ~reasons.astype(bool)
    overall = fail.sum(axis=0).sort_values(ascending=False).to_frame("n_failed")

    # drop reason by assay (optional)
    if assay_key in obs.columns:
        tmp = fail.copy()
        tmp["assay"] = obs[assay_key].astype(str).values
        by_assay = tmp.groupby("assay").sum(numeric_only=True)
    else:
        by_assay = pd.DataFrame()

    adata_qc = adata[keep].copy() if copy else adata[keep]
    # record
    adata_qc.uns["qc_rules"] = {
        "min_genes": rules.min_genes,
        "max_genes": rules.max_genes,
        "min_counts": rules.min_counts,
        "max_mt": rules.max_mt,
        "max_ercc": rules.max_ercc,
        "assay_key": assay_key,
        "per_assay": rules.per_assay,
    }
    return adata_qc, ret, (overall, by_assay)

