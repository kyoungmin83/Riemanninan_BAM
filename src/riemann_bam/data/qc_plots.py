import os
import matplotlib.pyplot as plt
import pandas as pd

QC_METRICS = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ercc"]

def plot_qc_histograms(adata, out_dir:str, bins: int = 60):
    os.makedirs(out_dir, exist_ok=True)
    obs = adata.obs

    for m in QC_METRICS:
        if m not in obs.columns:
            continue
        x = pd.to_numeric(obs[m], errors="coerce").dropna()

        plt.figure()
        plt.hist(x, bins=bins)
        plt.title(f"QC Distribution: {m} (n={len(x)})")
        plt.xlabel(m)
        plt.ylabel("count")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"qc_hist_{m}.png"),dpi=200)

        plt.close()

def plot_qc_by_group_box(adata, groupby: str, metric: str, out_dir: str, max_groups: int = 30):
    os.makedirs(out_dir, exist_ok=True)
    obs = adata.obs
    if groupby not in obs.columns or metric not in obs.columns:
        return
    df = obs[[groupby, metric]].copy()
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    df = df.dropna()

    # Ordered by size pf group
    top = df[groupby].value_counts().head(max_groups).index
    df = df[df[groupby].isin(top)]

    groups = [g for g in top]
    data = [df[df[groupby] == g][metric].values for g in groups]

    plt.figure(figsize=(max(8, len(groups)*0.4), 4))
    plt.boxplot(data, labels=groups, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.title(f"{metric} by {groupby} (top {len(groups)})")
    plt.ylabel(metric)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"qc_box_{metric}_by_{groupby}.png"), dpi=200)
    plt.close()
