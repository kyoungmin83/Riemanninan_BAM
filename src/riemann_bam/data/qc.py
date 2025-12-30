from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import anndata as ad

# Criteria list for QC
DEFAULT_QC_COLS = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ercc",
]

@dataclass
class QCSummary:
    groupby: Optional[str] #allow only str and None
    n_cells: int
    table: pd.DataFrame

def _quantiles(series: pd.Series, qs: Iterable[float]) -> pd.Series:
    series = pd.to_numeric(series, errors="coerce")
    return series.quantile(list(qs))

def summarize_qc(adata: ad.AnnData, qc_cols: Iterable[str] = DEFAULT_QC_COLS,
                  groupby: Optional[str] = None, 
                  quantiles: Iterable[float] = (0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999),)-> QCSummary:
    """Summarize QC metrics for an AnnData object."""
    obs = adata.obs.copy()
    qc_cols = [c for c in qc_cols if c in obs.columns]

    if len(qc_cols) == 0:
        raise ValueError("No QC columns found in adata.obs. Check qc_cols / adata.obs.columns.")
    
    qs = tuple(float(q) for q in quantiles)

    if groupby is None:
        rows =[]
        for c in qc_cols:
            qv = _quantiles(obs[c], qs)
            rows.append(pd.DataFrame({"metric":c, "q":qv.index, "value":qv.values}))
        out = pd.concat(rows, ignore_index=True).pivot(index="metric", columns="q", values="value")
        out.columns = [f"q{str(col).replace('0.', '.').replace('.', '')}" if col < 1 else f"q{int(col)}" for col in out.columns]
        out.insert(0, "n_nonnull", obs[qc_cols].notna().sum().min())
        return QCSummary(groupby=None, n_cells=adata.n_obs, table=out)
    
    if groupby not in obs.columns:
        raise ValueError(f"groupby='{groupby}' not in adata.obs.columns")

    rows = []
    grouped = obs.groupby(groupby, dropna=False)
    for g, df in grouped:
        for c in qc_cols:
            qv = _quantiles(df[c], qs)
            tmp = pd.DataFrame(
                {"group": str(g), "metric": c, "q": qv.index, "value": qv.values, "n_cells": len(df)}
            )
            rows.append(tmp)

    out = pd.concat(rows, ignore_index=True)
    pivot = out.pivot_table(index=["group", "metric"], columns="q", values="value", aggfunc="first")
    pivot.columns = [f"q{str(col).replace('0.', '.').replace('.', '')}" if col < 1 else f"q{int(col)}" for col in pivot.columns]
    
    sizes = out.groupby("group")["n_cells"].first()
    pivot.insert(0, "n_cells", pivot.index.get_level_values("group").map(sizes))
    return QCSummary(groupby=groupby, n_cells=adata.n_obs, table=pivot)

def pretty_print(summary: QCSummary, max_groups: int = 30) -> None:

    print("=" * 80)
    if summary.groupby is None:
        print(f"[QC SUMMARY] all cells | n_cells={summary.n_cells}")
        print(summary.table.round(3))
        return

    print(f"[QC SUMMARY] groupby='{summary.groupby}' | n_cells={summary.n_cells}")
    table = summary.table
    groups = table.index.get_level_values("group").unique()

    if len(groups) > max_groups:
        print(f"NOTE: groups={len(groups)} > max_groups={max_groups}. Showing first {max_groups} groups.")
        groups = groups[:max_groups]

  
    for g in groups:
        sub = table.loc[g]
        print("-" * 80)
        print(f"GROUP: {g} | n_cells={int(sub['n_cells'].iloc[0]) if 'n_cells' in sub.columns else 'NA'}")
        sub2 = sub.drop(columns=["n_cells"], errors="ignore")
        print(sub2.round(3))

