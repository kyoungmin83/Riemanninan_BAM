"""
Per-gene marginal bin frequency helper.

Used to initialise the LieActionOrdinalDecoder's thresholds at data-derived
quantile positions so that, at score=0, the decoder already reproduces the
empirical per-gene bin distribution.

See `doc/threshold_quantile_init.md` for the philosophy.
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np

from kmlee_bam.data.ordinal_dataset import OrdinalScDataset


def compute_bin_marginal_freq(
    dataset: OrdinalScDataset,
    *,
    n_sample_cells: int = 4096,
    seed: int = 42,
    smoothing: float = 1.0,
) -> np.ndarray:
    """
    Estimate per-gene marginal ordinal-bin frequency [G, K].

    Parameters
    ----------
    dataset :
        OrdinalScDataset built on the same zarr + spec as training.
    n_sample_cells :
        Number of cells to sample for the estimate. 4096 is enough to
        stabilise the per-gene distribution for K=7 bins.
    seed :
        Random seed for sampling.
    smoothing :
        Laplace smoothing added to every (gene, bin) count before
        normalisation. Prevents degenerate logit when a bin has 0
        observed cells.

    Returns
    -------
    freq : np.ndarray, shape [G, K], dtype float32
        Per-gene marginal probability over bins, summing to 1 along
        the bin axis.

    Notes
    -----
    Uses the chunked CSR reader directly to avoid the per-call overhead
    of `__getitem__` (which also builds pathology mask tensors we don't
    need). This makes a 4096-cell pass on SEA-AD ~30 seconds.
    """
    if n_sample_cells <= 0:
        raise ValueError("n_sample_cells must be positive.")
    if smoothing < 0:
        raise ValueError("smoothing must be non-negative.")

    n_total = len(dataset)
    if n_total == 0:
        raise RuntimeError("dataset is empty; cannot estimate bin marginals.")
    n_sample = int(min(n_sample_cells, n_total))

    rng = np.random.default_rng(int(seed))
    # `dataset.row_idx` indexes into the underlying zarr matrix.
    sample_ds_idx = rng.choice(n_total, size=n_sample, replace=False)
    sample_row_idx = np.asarray(dataset.row_idx, dtype=np.int64)[sample_ds_idx]

    g = len(dataset.spec.gene_names)
    k = int(dataset.spec.n_bins)
    edges = np.asarray(dataset.edges, dtype=np.float32)

    counts = np.full((g, k), float(smoothing), dtype=np.float64)

    # Read in chunks. ZarrMatrixReader can handle large row arrays; we
    # still split to keep per-call memory bounded.
    chunk = 256
    for start in range(0, n_sample, chunk):
        rows = np.sort(sample_row_idx[start : start + chunk])
        x_counts = dataset.X.rows_to_dense(rows)
        x_log = np.log1p(x_counts.astype(np.float32, copy=False))
        # bin = number of edges below value (same as OrdinalScDataset._bin_one_cell_from_log)
        # vectorised over [chunk, G]
        bins = (x_log[:, :, None] > edges[None, :, :]).sum(axis=-1).astype(np.int64)
        # Accumulate per-(gene, bin)
        for cell_bin_row in bins:
            np.add.at(counts, (np.arange(g), cell_bin_row), 1.0)

    freq = counts / counts.sum(axis=1, keepdims=True).clip(min=1e-12)
    return freq.astype(np.float32)


def save_bin_marginal_freq(
    freq: np.ndarray,
    gene_names: np.ndarray,
    path: str,
) -> None:
    """Persist a computed marginal so subsequent runs can skip the pass."""
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    np.savez_compressed(
        path,
        freq=freq.astype(np.float32),
        gene_names=np.asarray(gene_names, dtype=object),
    )


def load_bin_marginal_freq(path: str) -> tuple[np.ndarray, np.ndarray]:
    obj = np.load(path, allow_pickle=True)
    freq = obj["freq"].astype(np.float32)
    gene_names = np.asarray(obj["gene_names"], dtype=object)
    return freq, gene_names


def maybe_load_or_compute_bin_marginal(
    dataset: OrdinalScDataset,
    *,
    cache_path: Optional[str] = None,
    n_sample_cells: int = 4096,
    seed: int = 42,
    smoothing: float = 1.0,
) -> np.ndarray:
    """
    Convenience: load from cache_path if it exists and gene_names match,
    else compute and (if path given) save.
    """
    if cache_path is not None:
        try:
            freq, cached_names = load_bin_marginal_freq(cache_path)
            cur_names = np.asarray(dataset.spec.gene_names, dtype=object)
            if np.array_equal(cached_names, cur_names) and freq.shape == (
                len(dataset.spec.gene_names),
                dataset.spec.n_bins,
            ):
                print(
                    f"[bin_marginal] loaded cached marginal from {cache_path}",
                    flush=True,
                )
                return freq
            print(
                f"[bin_marginal] cache at {cache_path} mismatched, recomputing",
                flush=True,
            )
        except FileNotFoundError:
            pass
        except Exception as e:
            print(
                f"[bin_marginal] cache load failed ({e}); recomputing",
                flush=True,
            )

    print(
        f"[bin_marginal] computing from {n_sample_cells} sampled cells "
        f"(seed={seed}, smoothing={smoothing})...",
        flush=True,
    )
    freq = compute_bin_marginal_freq(
        dataset,
        n_sample_cells=n_sample_cells,
        seed=seed,
        smoothing=smoothing,
    )
    if cache_path is not None:
        save_bin_marginal_freq(freq, np.asarray(dataset.spec.gene_names, dtype=object), cache_path)
        print(f"[bin_marginal] cached marginal to {cache_path}", flush=True)
    return freq
