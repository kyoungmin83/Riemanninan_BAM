#%%
import os
import scanpy as sc
#from riemann_bam.data.preprocess import load_counts
from riemann_bam.data.qc_filter import QCRules, filter_adata

# %%
INPUT = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/train_raw.h5ad")
OUTPUT = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/TS_Blood_qc.h5ad")
OUT_DIR = os.path.expanduser("~/project_local/riemann_bam/data/tabula_blood_normal/outputs/qc_filter")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    adata = sc.read_h5ad(INPUT)

    rules = QCRules(
        min_genes=400,
        max_genes=7000,
        min_counts=500,
        max_mt=10.0,
        max_ercc = None,
        per_assay={
            "10x 3' v3": {"max_counts": 50000},
            "10x 5' v2": {"max_counts": 40000},
            "Smart-seq2": {"max_counts": 6500000, "max_mt": 15.0, "max_ercc":40.0, "min_genes":500},
        }
    )

    adata_qc, retention, (overall_fail, fail_by_assay) =  filter_adata(adata, rules, assay_key="assay")

    adata_qc.write_h5ad(OUTPUT)

    retention.to_csv(os.path.join(OUT_DIR, "retention_by_assay.csv"), index=False)
    overall_fail.to_csv(os.path.join(OUT_DIR, "fail_counts_overall.csv"))
    if  len(fail_by_assay) > 0:
        fail_by_assay.to_csv(os.path.join(OUT_DIR, "fail_counts_by_assay.csv"))

    print("Saved:", OUTPUT)
    print("\nRetention by assay:")
    print(retention)
    print("\nOverall fail counts (which rules drop most cells):")
    print(overall_fail.head(10))

if __name__ == "__main__":
    main()
