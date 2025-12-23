import os
import scanpy as sc

from riemann_bam.data.preprocess import load_counts
from riemann_bam.data.qc import summarize_qc, pretty_print


INPUT_PATH = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/TS_Blood.h5ad")
OUT_DIR = os.path.expanduser("~/project/riemann_bam/outputs/qc_reports")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    adata = load_counts(INPUT_PATH)
    print(adata)


    s_all = summarize_qc(adata, groupby=None)
    pretty_print(s_all)

    s_all.table.to_csv(os.path.join(OUT_DIR, "qc_summary_all.csv"))

  
    if "assay" in adata.obs.columns:
        s_assay = summarize_qc(adata, groupby="assay")
        pretty_print(s_assay, max_groups=50)
        s_assay.table.to_csv(os.path.join(OUT_DIR, "qc_summary_by_assay.csv"))
    else:
        print("WARNING: 'assay' not found in adata.obs. Skipping assay group summary.")

    print(f"\nSaved QC reports to: {OUT_DIR}")


if __name__ == "__main__":
    main()
