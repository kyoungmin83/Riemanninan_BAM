from __future__ import annotations

"""
Zarr-backed OrdinalScDataset for ordinal single-cell expression modelling.

This version avoids loading the full h5ad matrix into memory. It reads obs/var
metadata from an AnnData Zarr store and streams rows from X, raw/X, or
layers/<counts_layer>.

Existing outputs:
    y_ord, celltype_id, batch_id, row_index, optional x_log1p

New optional output:
    x_gene_scalar [G]

x_gene_scalar is intended for GeneModuleTokenizer.module_activity, not for
ordinal binning itself. The recommended scaler is celltype x control-train
reference standardisation.
"""

import json
import math
import warnings
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
import zarr
try:
    from kmlee_bam.preprocessing.data.pathology_masks import build_pathology_masks
except ImportError:
    from pathology_masks import build_pathology_masks


# =====================================================================
# Ordinal spec
# =====================================================================
@dataclass
class OrdinalSpec:
    gene_names: np.ndarray
    edges: np.ndarray          # [G, K-1]
    celltype_vocab: np.ndarray
    batch_vocab: np.ndarray
    n_bins: int

    @classmethod
    def load(cls, npz_path: str) -> "OrdinalSpec":
        obj = np.load(npz_path, allow_pickle=True)
        return cls(
            gene_names=np.asarray(obj["gene_names"], dtype=object),
            edges=obj["edges"].astype(np.float32),
            celltype_vocab=np.asarray(obj["celltype_vocab"], dtype=object),
            batch_vocab=np.asarray(obj["batch_vocab"], dtype=object),
            n_bins=int(np.asarray(obj["n_bins"]).ravel()[0]),
        )


# =====================================================================
# Reference scaler
# =====================================================================
@dataclass
class ReferenceGeneScaler:
    gene_names: np.ndarray
    mean: np.ndarray
    std: np.ndarray
    celltype_vocab: Optional[np.ndarray] = None
    is_global: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)
    n_control: Optional[np.ndarray] = None

    @classmethod
    def load(cls, npz_path: str) -> "ReferenceGeneScaler":
        obj = np.load(npz_path, allow_pickle=True)

        if "gene_names" not in obj:
            raise KeyError("reference scaler npz must contain 'gene_names'.")
        gene_names = np.asarray(obj["gene_names"], dtype=object)

        mean_key = _first_present_key(obj, ["mean_control_train", "mean"])
        std_key = _first_present_key(obj, ["std_control_train", "std"])
        if mean_key is None:
            raise KeyError("reference scaler npz must contain 'mean_control_train' or 'mean'.")
        if std_key is None:
            raise KeyError("reference scaler npz must contain 'std_control_train' or 'std'.")

        mean = obj[mean_key].astype(np.float32)
        std = obj[std_key].astype(np.float32)

        if mean.shape != std.shape:
            raise ValueError(f"reference scaler mean/std shape mismatch: mean={mean.shape}, std={std.shape}.")
        if mean.ndim not in {1, 2}:
            raise ValueError("reference scaler mean/std must have shape [G] or [C, G].")
        if mean.shape[-1] != len(gene_names):
            raise ValueError(
                "reference scaler gene dimension mismatch: "
                f"mean.shape[-1]={mean.shape[-1]}, len(gene_names)={len(gene_names)}."
            )
        if not np.isfinite(mean).all():
            raise ValueError("reference scaler mean contains non-finite values.")
        if not np.isfinite(std).all():
            raise ValueError("reference scaler std contains non-finite values.")
        if (std < 0).any():
            raise ValueError("reference scaler std contains negative values.")

        is_global = mean.ndim == 1

        if "celltype_vocab" in obj:
            celltype_vocab = np.asarray(obj["celltype_vocab"], dtype=object)
        elif is_global:
            celltype_vocab = np.asarray(["__GLOBAL__"], dtype=object)
        else:
            raise KeyError("celltype-aware reference scaler must contain 'celltype_vocab'.")

        if not is_global and len(celltype_vocab) != mean.shape[0]:
            raise ValueError(
                "reference scaler celltype dimension mismatch: "
                f"len(celltype_vocab)={len(celltype_vocab)}, mean.shape[0]={mean.shape[0]}."
            )

        metadata = _load_scaler_metadata(obj)
        n_control = None
        n_key = _first_present_key(obj, ["n_control", "n_control_cells"])
        if n_key is not None:
            n_control = np.asarray(obj[n_key])

        return cls(
            gene_names=gene_names,
            mean=mean,
            std=std,
            celltype_vocab=celltype_vocab,
            is_global=is_global,
            metadata=metadata,
            n_control=n_control,
        )


def _first_present_key(npz_obj: Any, candidates: list[str]) -> Optional[str]:
    for key in candidates:
        if key in npz_obj:
            return key
    return None


def _to_python_scalar_or_list(x: np.ndarray) -> Any:
    arr = np.asarray(x)
    if arr.ndim == 0:
        return arr.item()
    return arr.tolist()


def _load_scaler_metadata(npz_obj: Any) -> dict[str, Any]:
    metadata: dict[str, Any] = {}

    if "metadata_json" in npz_obj:
        raw = np.asarray(npz_obj["metadata_json"]).item()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        if isinstance(raw, str) and raw.strip():
            try:
                parsed = json.loads(raw)
                if isinstance(parsed, dict):
                    metadata.update(parsed)
            except json.JSONDecodeError:
                metadata["metadata_json_parse_error"] = str(raw)

    metadata_keys = [
        "control_definition", "control_key", "control_values",
        "donor_balanced", "donor_key", "shrinkage_used", "shrinkage_lambda",
        "n_train_donors", "n_train_cells", "n_control_donors", "n_control_cells_total",
        "scaler_version", "build_date",
    ]
    for key in metadata_keys:
        if key in npz_obj:
            metadata[key] = _to_python_scalar_or_list(npz_obj[key])
    return metadata


def _build_celltype_to_scaler_row(
    spec_celltype_vocab: np.ndarray,
    scaler_celltype_vocab: np.ndarray,
) -> tuple[np.ndarray, list[str]]:
    scaler_map = {str(v): i for i, v in enumerate(scaler_celltype_vocab.tolist())}
    global_idx = scaler_map.get("__GLOBAL__")

    out = np.empty(len(spec_celltype_vocab), dtype=np.int64)
    missing: list[str] = []
    fallback_used: list[str] = []

    for i, ct in enumerate(spec_celltype_vocab.tolist()):
        key = str(ct)
        if key in scaler_map:
            out[i] = scaler_map[key]
        elif global_idx is not None:
            out[i] = int(global_idx)
            fallback_used.append(key)
        else:
            missing.append(key)

    if missing:
        preview = missing[:10]
        suffix = "..." if len(missing) > 10 else ""
        raise KeyError(
            "Reference scaler is missing celltype(s) and has no '__GLOBAL__' fallback: "
            f"{preview}{suffix}"
        )

    return out, fallback_used


# =====================================================================
# Zarr metadata / matrix helpers
# =====================================================================
def _read_elem(group):
    try:
        from anndata.experimental import read_elem  # type: ignore
        return read_elem(group)
    except Exception as e:
        raise RuntimeError(
            "Could not read Zarr metadata via anndata.experimental.read_elem. "
            "Install a recent anndata version."
        ) from e


def read_zarr_dataframe(root, key_path: str) -> pd.DataFrame:
    group = root
    for part in key_path.split("/"):
        if part == "":
            continue
        group = group[part]
    out = _read_elem(group)
    if not isinstance(out, pd.DataFrame):
        raise TypeError(f"Zarr element '{key_path}' is not a pandas DataFrame.")
    return out.copy()


class ZarrMatrixReader:
    """Read one or more rows from a dense or CSR sparse AnnData Zarr matrix."""

    def __init__(
        self,
        matrix,
        *,
        n_cols: int,
        name: str,
        max_csr_data_span: int = 50_000_000,
        max_span_to_selected_nnz_ratio: float = 20.0,
    ) -> None:
        self.matrix = matrix
        self.n_cols = int(n_cols)
        self.name = name
        self.max_csr_data_span = int(max_csr_data_span)
        self.max_span_to_selected_nnz_ratio = float(max_span_to_selected_nnz_ratio)

        self.is_csr = hasattr(matrix, "keys") and all(k in matrix for k in ("data", "indices", "indptr"))
        if self.is_csr:
            self.data = matrix["data"]
            self.indices = matrix["indices"]
            self.indptr = matrix["indptr"]
            self.shape = (len(self.indptr) - 1, self.n_cols)
        else:
            self.shape = tuple(matrix.shape)
            if len(self.shape) != 2:
                raise ValueError(f"Matrix '{name}' must be 2D or CSR group, got shape={self.shape}.")
            if self.shape[1] != self.n_cols:
                raise ValueError(f"Matrix '{name}' gene dimension mismatch: shape[1]={self.shape[1]}, n_cols={self.n_cols}.")

    def row_to_dense(self, row: int) -> np.ndarray:
        return self.rows_to_dense(np.asarray([int(row)], dtype=np.int64))[0]

    def rows_to_dense(self, rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        if rows.ndim != 1:
            raise ValueError("rows must be 1D.")
        if len(rows) == 0:
            return np.zeros((0, self.n_cols), dtype=np.float32)
        if (rows < 0).any() or (rows >= self.shape[0]).any():
            raise IndexError("row index out of bounds for zarr matrix.")
        if self.is_csr:
            return self._csr_rows_to_dense(rows)
        return self._dense_rows_to_dense(rows)

    def _dense_rows_to_dense(self, rows: np.ndarray) -> np.ndarray:
        try:
            block = self.matrix.get_orthogonal_selection((rows, slice(None)))
        except Exception:
            try:
                block = self.matrix[rows, :]
            except Exception:
                block = np.vstack([np.asarray(self.matrix[int(r), :]) for r in rows])
        return np.asarray(block, dtype=np.float32)

    def _csr_rows_to_dense(self, rows: np.ndarray) -> np.ndarray:
        order = np.argsort(rows)
        rows_sorted = rows[order]
        base = int(rows_sorted[0])
        end_row = int(rows_sorted[-1])

        indptr_slice = np.asarray(self.indptr[base : end_row + 2], dtype=np.int64)
        starts = indptr_slice[rows_sorted - base]
        stops = indptr_slice[rows_sorted - base + 1]
        selected_nnz = int(np.sum(stops - starts))
        span_start = int(starts[0])
        span_end = int(stops[-1])
        span_len = max(0, span_end - span_start)

        out_sorted = np.zeros((len(rows_sorted), self.n_cols), dtype=np.float32)
        use_per_row = False
        if selected_nnz > 0:
            use_per_row = (
                span_len > self.max_csr_data_span
                and span_len > self.max_span_to_selected_nnz_ratio * selected_nnz
            )

        if selected_nnz == 0:
            pass
        elif use_per_row:
            for j, (s, e) in enumerate(zip(starts, stops)):
                if e <= s:
                    continue
                idx = np.asarray(self.indices[int(s):int(e)], dtype=np.int64)
                dat = np.asarray(self.data[int(s):int(e)], dtype=np.float32)
                out_sorted[j, idx] = dat
        else:
            idx_span = np.asarray(self.indices[span_start:span_end], dtype=np.int64)
            dat_span = np.asarray(self.data[span_start:span_end], dtype=np.float32)
            for j, (s, e) in enumerate(zip(starts, stops)):
                if e <= s:
                    continue
                a = int(s - span_start)
                b = int(e - span_start)
                out_sorted[j, idx_span[a:b]] = dat_span[a:b]

        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        return out_sorted[inv]


def _zarr_get(root, path: str):
    obj = root
    for part in path.split("/"):
        if part == "":
            continue
        obj = obj[part]
    return obj


def choose_zarr_matrix(
    root,
    *,
    matrix_path: str = "auto",
    counts_layer: str = "counts",
    prefer_raw: bool = True,
) -> tuple[object, str, str]:
    requested = str(matrix_path).strip()
    if requested and requested != "auto":
        matrix = _zarr_get(root, requested)
        var_key = "raw/var" if requested.startswith("raw/") else "var"
        return matrix, requested, var_key

    layer_path = f"layers/{counts_layer}"
    if "layers" in root and counts_layer in root["layers"]:
        return root["layers"][counts_layer], layer_path, "var"

    if prefer_raw and "raw" in root and "X" in root["raw"]:
        return root["raw"]["X"], "raw/X", "raw/var"

    if "X" not in root:
        raise KeyError("Could not find matrix path: no layers/counts, raw/X, or X in zarr store.")
    return root["X"], "X", "var"


# =====================================================================
# Dataset
# =====================================================================
class OrdinalScDataset(Dataset):
    def __init__(
        self,
        zarr_path: str,
        spec_path: str,
        *,
        split: str,
        matrix_path: str = "auto",
        counts_layer: str = "counts",
        prefer_raw: bool = True,
        split_key: str = "split",
        celltype_key: str = "Subclass",
        batch_key: str = "donor_id",
        return_log_input: bool = True,
        reference_scaler_path: Optional[str] = None,
        return_x_gene_scalar: bool = False,
        x_gene_scalar_clip: Optional[float] = 10.0,
        x_gene_scalar_eps: float = 1e-4,
        warn_min_control_cells: Optional[int] = 30,
        max_csr_data_span: int = 50_000_000,
        max_span_to_selected_nnz_ratio: float = 20.0,
    ):
        self.zarr_path = str(zarr_path)
        self.root = zarr.open(self.zarr_path, mode="r")
        self.spec = OrdinalSpec.load(spec_path)
        obs = read_zarr_dataframe(self.root, "obs")
        matrix, self.matrix_source, var_key = choose_zarr_matrix(
            self.root,
            matrix_path=matrix_path,
            counts_layer=counts_layer,
            prefer_raw=prefer_raw,
        )
        var = read_zarr_dataframe(self.root, var_key)
        var_names = np.asarray(var.index.astype(str), dtype=object)

        for key in [split_key, celltype_key, batch_key]:
            if key not in obs.columns:
                raise KeyError(f"'{key}' not found in zarr obs")
            
        # Pathology / reference masks aligned to full obs rows.
        self.pathology_masks = build_pathology_masks(obs)

        self.is_reference_origin = self.pathology_masks["is_reference_origin"].astype(np.bool_)
        self.is_abnormal_pathology = self.pathology_masks["is_abnormal_pathology"].astype(np.bool_)
            
        self.axis_late_positive = self.pathology_masks["axis_late_positive"].astype(np.bool_)
        self.axis_late_advanced = self.pathology_masks["axis_late_advanced"].astype(np.bool_)
        self.axis_braak_mid = self.pathology_masks["axis_braak_mid"].astype(np.bool_)
        self.axis_braak_high = self.pathology_masks["axis_braak_high"].astype(np.bool_)
        self.axis_thal_positive = self.pathology_masks["axis_thal_positive"].astype(np.bool_)
        self.axis_thal_high = self.pathology_masks["axis_thal_high"].astype(np.bool_)
        self.axis_cerad_positive = self.pathology_masks["axis_cerad_positive"].astype(np.bool_)
        self.axis_cerad_moderate_or_frequent = self.pathology_masks["axis_cerad_moderate_or_frequent"].astype(np.bool_)
        self.axis_cerad_frequent = self.pathology_masks["axis_cerad_frequent"].astype(np.bool_)
        self.axis_adnc_intermediate_high = self.pathology_masks["axis_adnc_intermediate_high"].astype(np.bool_)
        self.axis_adnc_high = self.pathology_masks["axis_adnc_high"].astype(np.bool_)
        self.axis_lewy_positive = self.pathology_masks["axis_lewy_positive"].astype(np.bool_)
        self.axis_microinfarct_positive = self.pathology_masks["axis_microinfarct_positive"].astype(np.bool_)
        self.axis_apoe4_carrier = self.pathology_masks["axis_apoe4_carrier"].astype(np.bool_)
    

        if len(var_names) != len(self.spec.gene_names):
            raise ValueError("Gene dimension mismatch between zarr var and spec")
        if not np.array_equal(var_names, np.asarray(self.spec.gene_names, dtype=object)):
            raise ValueError("Gene order mismatch between zarr var and spec")

        self.X = ZarrMatrixReader(
            matrix,
            n_cols=len(self.spec.gene_names),
            name=self.matrix_source,
            max_csr_data_span=max_csr_data_span,
            max_span_to_selected_nnz_ratio=max_span_to_selected_nnz_ratio,
        )

        self.row_idx = np.where(obs[split_key].astype(str).values == split)[0]
        if len(self.row_idx) == 0:
            raise ValueError(f"No rows found for split='{split}'")

        if x_gene_scalar_eps <= 0:
            raise ValueError("x_gene_scalar_eps must be positive.")
        if x_gene_scalar_clip is not None and x_gene_scalar_clip <= 0:
            raise ValueError("x_gene_scalar_clip must be positive when not None.")

        self.return_log_input = bool(return_log_input)
        self.return_x_gene_scalar = bool(return_x_gene_scalar)
        self.x_gene_scalar_clip = x_gene_scalar_clip
        self.x_gene_scalar_eps = float(x_gene_scalar_eps)
        self.warn_min_control_cells = warn_min_control_cells
        self.edges = self.spec.edges.astype(np.float32, copy=False)

        celltype_to_id = {str(v): i for i, v in enumerate(self.spec.celltype_vocab.tolist())}
        batch_to_id = {str(v): i for i, v in enumerate(self.spec.batch_vocab.tolist())}

        self.celltype_ids = np.array(
            [celltype_to_id[str(v)] for v in obs[celltype_key].astype(str).values],
            dtype=np.int64,
        )
        self.batch_ids = np.array(
            [batch_to_id[str(v)] for v in obs[batch_key].astype(str).values],
            dtype=np.int64,
        )

        self.reference_scaler: Optional[ReferenceGeneScaler] = None
        self.celltype_to_scaler_row: Optional[np.ndarray] = None
        self.fallback_celltypes: list[str] = []
        self._scaler_mean: Optional[np.ndarray] = None
        self._scaler_std_eps: Optional[np.ndarray] = None

        if reference_scaler_path is not None:
            self.reference_scaler = ReferenceGeneScaler.load(reference_scaler_path)
            self._prepare_reference_scaler_cache()

        if self.return_x_gene_scalar and self.reference_scaler is None:
            raise ValueError("return_x_gene_scalar=True requires reference_scaler_path to be provided.")

    def _prepare_reference_scaler_cache(self) -> None:
        assert self.reference_scaler is not None

        if not np.array_equal(
            np.array(self.reference_scaler.gene_names, dtype=object),
            np.array(self.spec.gene_names, dtype=object),
        ):
            raise ValueError("Gene order mismatch between reference scaler and ordinal spec.")

        scaler = self.reference_scaler
        if scaler.is_global:
            mean_by_ct = np.broadcast_to(
                scaler.mean.reshape(1, -1),
                (len(self.spec.celltype_vocab), len(self.spec.gene_names)),
            )
            std_by_ct = np.broadcast_to(
                scaler.std.reshape(1, -1),
                (len(self.spec.celltype_vocab), len(self.spec.gene_names)),
            )
            self.celltype_to_scaler_row = np.zeros(len(self.spec.celltype_vocab), dtype=np.int64)
            self.fallback_celltypes = []
        else:
            assert scaler.celltype_vocab is not None
            self.celltype_to_scaler_row, self.fallback_celltypes = _build_celltype_to_scaler_row(
                self.spec.celltype_vocab,
                scaler.celltype_vocab,
            )
            if self.fallback_celltypes:
                preview = self.fallback_celltypes[:10]
                suffix = "..." if len(self.fallback_celltypes) > 10 else ""
                warnings.warn(
                    "Reference scaler used '__GLOBAL__' fallback for "
                    f"{len(self.fallback_celltypes)} celltype(s): {preview}{suffix}",
                    RuntimeWarning,
                )
            mean_by_ct = scaler.mean[self.celltype_to_scaler_row]
            std_by_ct = scaler.std[self.celltype_to_scaler_row]

        self._scaler_mean = np.ascontiguousarray(mean_by_ct, dtype=np.float32)
        self._scaler_std_eps = np.ascontiguousarray(std_by_ct + self.x_gene_scalar_eps, dtype=np.float32)
        self._warn_if_low_control_counts()

    def _warn_if_low_control_counts(self) -> None:
        if self.reference_scaler is None or self.reference_scaler.n_control is None:
            return
        if self.warn_min_control_cells is None:
            return

        n_control = np.asarray(self.reference_scaler.n_control)
        if n_control.ndim == 0:
            return
        threshold = int(self.warn_min_control_cells)
        if threshold <= 0:
            return

        if n_control.ndim == 2:
            n_by_row = np.min(n_control, axis=1)
        elif n_control.ndim == 1:
            n_by_row = n_control
        else:
            return

        if self.reference_scaler.is_global:
            low = np.where(n_by_row < threshold)[0]
            if len(low) > 0:
                warnings.warn("Global reference scaler has low n_control; x_gene_scalar may be noisy.", RuntimeWarning)
            return

        if self.celltype_to_scaler_row is None:
            return

        mapped_counts = n_by_row[self.celltype_to_scaler_row]
        low_spec_rows = np.where(mapped_counts < threshold)[0]
        if len(low_spec_rows) > 0:
            low_names = [str(self.spec.celltype_vocab[i]) for i in low_spec_rows[:10]]
            suffix = "..." if len(low_spec_rows) > 10 else ""
            warnings.warn(
                f"{len(low_spec_rows)} celltype(s) have n_control < {threshold} in scaler: "
                f"{low_names}{suffix}. x_gene_scalar may be noisy.",
                RuntimeWarning,
            )

    def __len__(self) -> int:
        return len(self.row_idx)

    def _get_row_dense(self, r: int) -> np.ndarray:
        return self.X.row_to_dense(int(r)).astype(np.float32, copy=False)

    def _bin_one_cell_from_log(self, x_log: np.ndarray) -> np.ndarray:
        return np.sum(x_log[:, None] > self.edges, axis=1).astype(np.int64)

    def _make_x_gene_scalar(self, x_log: np.ndarray, celltype_id: int) -> np.ndarray:
        if self._scaler_mean is None or self._scaler_std_eps is None:
            raise RuntimeError("Reference scaler cache is not initialised.")

        x_scalar = ((x_log - self._scaler_mean[celltype_id]) / self._scaler_std_eps[celltype_id]).astype(
            np.float32,
            copy=False,
        )
        if self.x_gene_scalar_clip is not None:
            clip = float(self.x_gene_scalar_clip)
            np.clip(x_scalar, -clip, clip, out=x_scalar)
        return x_scalar

    def __getitem__(self, idx: int):
        r = int(self.row_idx[idx])
        x_counts = self._get_row_dense(r)
        x_log = np.log1p(x_counts).astype(np.float32, copy=False)
        y_ord = self._bin_one_cell_from_log(x_log)
        celltype_id = int(self.celltype_ids[r])

        item = {
            "y_ord": torch.from_numpy(y_ord),
            "celltype_id": torch.tensor(celltype_id, dtype=torch.long),
            "batch_id": torch.tensor(self.batch_ids[r], dtype=torch.long),
            "row_index": torch.tensor(r, dtype=torch.long),
        }

        item.update(
            {
                "is_reference": torch.tensor(bool(self.is_reference_origin[r]), dtype=torch.bool),
                "is_abnormal": torch.tensor(bool(self.is_abnormal_pathology[r]), dtype=torch.bool),

                "axis_late_positive": torch.tensor(bool(self.axis_late_positive[r]), dtype=torch.bool),
                "axis_late_advanced": torch.tensor(bool(self.axis_late_advanced[r]), dtype=torch.bool),
                "axis_braak_mid": torch.tensor(bool(self.axis_braak_mid[r]), dtype=torch.bool),
                "axis_braak_high": torch.tensor(bool(self.axis_braak_high[r]), dtype=torch.bool),
                "axis_thal_positive": torch.tensor(bool(self.axis_thal_positive[r]), dtype=torch.bool),
                "axis_thal_high": torch.tensor(bool(self.axis_thal_high[r]), dtype=torch.bool),
                "axis_cerad_positive": torch.tensor(bool(self.axis_cerad_positive[r]), dtype=torch.bool),
                "axis_cerad_moderate_or_frequent": torch.tensor(bool(self.axis_cerad_moderate_or_frequent[r]), dtype=torch.bool),
                "axis_cerad_frequent": torch.tensor(bool(self.axis_cerad_frequent[r]), dtype=torch.bool),
                "axis_adnc_intermediate_high": torch.tensor(bool(self.axis_adnc_intermediate_high[r]), dtype=torch.bool),
                "axis_adnc_high": torch.tensor(bool(self.axis_adnc_high[r]), dtype=torch.bool),
                "axis_lewy_positive": torch.tensor(bool(self.axis_lewy_positive[r]), dtype=torch.bool),
                "axis_microinfarct_positive": torch.tensor(bool(self.axis_microinfarct_positive[r]), dtype=torch.bool),
                "axis_apoe4_carrier": torch.tensor(bool(self.axis_apoe4_carrier[r]), dtype=torch.bool),
            }
        )
        if self.return_log_input:
            item["x_log1p"] = torch.from_numpy(x_log)
        if self.return_x_gene_scalar:
            item["x_gene_scalar"] = torch.from_numpy(self._make_x_gene_scalar(x_log, celltype_id))
        return item


if __name__ == "__main__":
    print("This module defines Zarr-backed OrdinalScDataset.")