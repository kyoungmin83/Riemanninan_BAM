import scanpy as sc
from scipy import sparse

def load_counts(path):
    adata = sc.read_h5ad(path)

    if "decountXcounts" in adata.layers:
        adata.X = adata.layers["decontXcounts"].copy()
        adata.uns["counts_source"] = "decontXcounts"
    else:
        adata.uns["count_source"] = "X"

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)

    return adata