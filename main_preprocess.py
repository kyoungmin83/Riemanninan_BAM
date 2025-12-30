import os
from riemann_bam.data.preprocess import load_counts

INPUT_PATH = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/TS_Blood.h5ad")
OUTPUT_PATH = os.path.expanduser("~/project/riemann_bam/data/tabula_blood_normal/train_raw.h5ad")

def main():
    adata = load_counts(INPUT_PATH)
    print(adata)
    print("count_source:", adata.uns["count_source"])
    adata.write_h5ad(OUTPUT_PATH)

if __name__ == "__main__":
    main()