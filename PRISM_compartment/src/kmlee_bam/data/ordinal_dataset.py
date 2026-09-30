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
from kmlee_bam.data.pathology_masks import build_pathology_masks, compute_residualized_targets


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


def _read_sex_male_ids(obs: pd.DataFrame) -> np.ndarray:
    """Per-cell sex id (0=female, 1=male; -1=unknown) aligned to full obs rows.

    Prefers the existing 0/1 ``axis_sex_male`` obs column (per the data-build
    convention: female=0, male=1). Falls back to the free-text ``sex`` column
    (``m*``→1, ``f*``→0) when the axis column is absent. Any cell whose value is
    not a clean 0/1 (missing / unparseable) gets the sentinel ``-1`` so the
    sex-alignment loss can simply treat it as its own ignored group rather than
    crash — matching how the dataset never breaks on a missing covariate. Output
    is ``int64`` with values in ``{-1, 0, 1}``.
    """
    n = len(obs)
    if "axis_sex_male" in obs.columns:
        vals = pd.to_numeric(obs["axis_sex_male"], errors="coerce")
        out = np.full(n, -1, dtype=np.int64)
        arr = vals.to_numpy()
        finite = np.isfinite(arr)
        rounded = np.where(finite, np.rint(arr), -1).astype(np.int64)
        valid = finite & np.isin(rounded, (0, 1))
        out[valid] = rounded[valid]
        return out
    if "sex" in obs.columns:
        s = obs["sex"].astype(str).str.strip().str.lower()
        out = np.full(n, -1, dtype=np.int64)
        out[s.str.startswith("m").to_numpy()] = 1
        out[s.str.startswith("f").to_numpy()] = 0
        return out
    return np.full(n, -1, dtype=np.int64)


def _read_age_years(obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Read age-at-death without silently treating a missing age as zero.

    SEA-AD carries an exact ``development_stage`` string for most donors and a
    top-coded ``Age at death`` category for the oldest donors.  We prefer the
    exact year, but deliberately replace the generic ``80 and over`` ontology
    label with the midpoint-like value 95 for the explicit ``90+`` category.
    The returned validity mask remains authoritative for any unresolved row.
    """

    n = len(obs)
    out = np.full(n, np.nan, dtype=np.float32)

    if "development_stage" in obs.columns:
        stage = obs["development_stage"].fillna("").astype(str).str.strip()
        exact = pd.to_numeric(
            stage.str.extract(r"(?i)(\d+(?:\.\d+)?)\s*-?\s*year-old", expand=False),
            errors="coerce",
        )
        top_coded = stage.str.lower().str.contains("and over", regex=False)
        exact[top_coded] = np.nan
        out = exact.to_numpy(dtype=np.float32)

    if "Age at death" in obs.columns:
        category = obs["Age at death"].fillna("").astype(str).str.strip().str.lower()
        fallback = np.full(n, np.nan, dtype=np.float32)
        interval = category.str.extract(
            r"(\d+(?:\.\d+)?)\s*to\s*(\d+(?:\.\d+)?)", expand=True
        )
        lo = pd.to_numeric(interval[0], errors="coerce").to_numpy(dtype=np.float32)
        hi = pd.to_numeric(interval[1], errors="coerce").to_numpy(dtype=np.float32)
        has_interval = np.isfinite(lo) & np.isfinite(hi)
        fallback[has_interval] = 0.5 * (lo[has_interval] + hi[has_interval])
        lower = pd.to_numeric(
            category.str.extract(r"(\d+(?:\.\d+)?)", expand=False), errors="coerce"
        ).to_numpy(dtype=np.float32)
        less = category.str.contains("less than", regex=False).to_numpy()
        plus = category.str.contains("+", regex=False).to_numpy()
        fallback[less & np.isfinite(lower)] = lower[less & np.isfinite(lower)] - 5.0
        fallback[plus & np.isfinite(lower)] = lower[plus & np.isfinite(lower)] + 5.0
        out = np.where(np.isfinite(out), out, fallback).astype(np.float32)

    valid = np.isfinite(out)
    return np.where(valid, out, 0.0).astype(np.float32), valid.astype(np.bool_)


# PRISM end-to-end named pathology contract.  ADNC is intentionally absent:
# it is a report-only composite of AD-core measures, not an independent input.
PRISM_PATHOLOGY_NAMES: tuple[str, ...] = ("thal", "braak", "cerad", "late", "lewy")


def _read_prism_pathology(obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return five ordinal pathology coordinates and a separate validity mask.

    Coordinates are normalised to [0,1] as Thal/5, Braak/6, CERAD/3,
    LATE/3, and binary Lewy.  A strict reference is a valid biological zero.
    Unstageable/missing observations are invalid and are never silently turned
    into normal; callers must carry the returned mask.
    """

    n = len(obs)
    raw = np.full((n, len(PRISM_PATHOLOGY_NAMES)), np.nan, dtype=np.float32)

    def text(column: str) -> pd.Series:
        if column not in obs.columns:
            return pd.Series([""] * n, index=obs.index, dtype=object)
        return obs[column].fillna("").astype(str).str.strip()

    thal = text("Thal phase")
    thal_num = pd.to_numeric(thal.str.extract(r"(?i)thal\s*([0-5])", expand=False), errors="coerce")
    thal_num[thal.str.lower().eq("reference")] = 0.0
    raw[:, 0] = thal_num.to_numpy(dtype=np.float32)

    braak = text("Braak stage")
    braak_map = {
        "braak 0": 0.0,
        "braak i": 1.0,
        "braak ii": 2.0,
        "braak iii": 3.0,
        "braak iv": 4.0,
        "braak v": 5.0,
        "braak vi": 6.0,
        "reference": 0.0,
    }
    raw[:, 1] = braak.str.lower().map(braak_map).to_numpy(dtype=np.float32)

    cerad = text("CERAD score")
    cerad_map = {
        "absent": 0.0,
        "sparse": 1.0,
        "moderate": 2.0,
        "frequent": 3.0,
        "reference": 0.0,
    }
    raw[:, 2] = cerad.str.lower().map(cerad_map).to_numpy(dtype=np.float32)

    late = text("LATE-NC stage")
    late_lower = late.str.lower()
    late_num = pd.to_numeric(
        late.str.extract(r"(?i)late\s*stage\s*([0-3])", expand=False), errors="coerce"
    )
    late_num[late_lower.eq("not identified") | late_lower.eq("reference")] = 0.0
    # Strings such as "Staging Precluded ..." deliberately remain NaN/invalid.
    raw[:, 3] = late_num.to_numpy(dtype=np.float32)

    lewy = text("Lewy body disease pathology")
    lewy_lower = lewy.str.lower()
    lewy_num = np.full(n, np.nan, dtype=np.float32)
    known_zero = lewy_lower.eq("reference") | lewy_lower.str.startswith("not identified")
    known_positive = (
        lewy_lower.str.contains("olfactory bulb only", regex=False)
        | lewy_lower.str.contains("brainstem-predominant", regex=False)
        | lewy_lower.str.contains("amygdala-predominant", regex=False)
        | lewy_lower.str.contains("limbic (transitional)", regex=False)
        | lewy_lower.str.contains("neocortical (diffuse)", regex=False)
    )
    lewy_num[known_zero.to_numpy()] = 0.0
    lewy_num[known_positive.to_numpy()] = 1.0
    raw[:, 4] = lewy_num

    valid = np.isfinite(raw)
    scale = np.asarray([5.0, 6.0, 3.0, 3.0, 1.0], dtype=np.float32)
    normalised = raw / scale[None, :]
    # Stored value for invalid coordinates is a harmless placeholder; validity
    # is the source of truth and the model masks the value before use.
    normalised = np.where(valid, np.clip(normalised, 0.0, 1.0), 0.0).astype(np.float32)
    return normalised, valid.astype(np.bool_)


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
        # Thinning (zero-origin work, 2026-06): when True, each item also carries
        # a *thinned* view — raw counts binomially down-sampled (each molecule
        # kept with prob `thin_rho`) and re-binned with the SAME edges/log1p
        # code. Used by V8Trainer for the depth-invariance consistency loss and
        # for proven technical-zero labels. Off ⇒ items are byte-identical.
        # See doc/zero_origin_and_capacity_design_2026-06-01.md.
        return_thinned: bool = False,
        thin_rho: float = 0.5,
        # Faithful thinning needs RAW integer counts. When the training matrix
        # is CP10K (fractional), point `thin_counts_zarr_path` at a row/gene-
        # aligned raw-count companion store; the thinned view is then rebuilt as
        # raw → binomial thin → CP10K renorm by the *full* library size
        # (`thin_library_size_key` in obs) → log1p → re-bin, matching the
        # log1p(CP10K) training scale. If left None, the legacy path thins the
        # main matrix directly (valid only when it is raw counts on log1p(raw)).
        # See doc/thinning_stageA_data_verification_2026-06-02.md.
        thin_counts_zarr_path: Optional[str] = None,
        thin_counts_matrix_path: str = "X",
        thin_library_size_key: str = "Number of UMIs",
    ):
        self.zarr_path = str(zarr_path)
        self.split = str(split)
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
        self.gene_symbols = np.asarray(
            (
                var["feature_name"].astype(str).values
                if "feature_name" in var.columns
                else var_names
            ),
            dtype=object,
        )

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
        self.axis_cognitive_impaired = self.pathology_masks["axis_cognitive_impaired"].astype(np.bool_)

        # PRISM E2E: five named, ordinal pathology coordinates.  These are
        # separate from the historical high/positive binary masks above.
        # ADNC remains report-only and is intentionally not emitted here.
        self.prism_pathology, self.prism_pathology_valid = _read_prism_pathology(obs)

        # v23: donor-level RESIDUALIZED pathology targets (axis MINUS what overall
        # severity predicts). Regression fit on TRAIN donors ONLY (no val/test label
        # leak), standardized per axis, missing-masked. Donor-constant, aligned to FULL
        # obs rows; the path_aux head applies a MASKED MSE on (donor x celltype)
        # AGGREGATES -> forces axis-SPECIFIC signal, blocks the 1-D severity collapse.
        # Unused unless pathology_aux.lambda_resid_axis > 0.
        try:
            _rax = np.stack([self.axis_braak_high, self.axis_thal_high,
                             self.axis_late_positive, self.axis_lewy_positive], axis=1).astype(float)
            _donor_str = obs["donor_id"].astype(str).values
            _is_train = obs[split_key].astype(str).values == "train"      # fit on TRAIN donors only
            _ud, _invd = np.unique(_donor_str, return_inverse=True)
            _nD = len(_ud)
            _cnt = np.bincount(_invd, minlength=_nD).clip(1, None)
            _Ad = (np.stack([np.bincount(_invd, weights=_rax[:, j], minlength=_nD) / _cnt
                             for j in range(4)], axis=1) >= 0.5).astype(float)
            _trd = np.zeros(_nD, bool); np.maximum.at(_trd, _invd, _is_train)
            _resid_d, _rval_d = compute_residualized_targets(_Ad, _trd, np.ones_like(_Ad, bool))
            self.axis_resid = np.nan_to_num(_resid_d[_invd], nan=0.0).astype(np.float32)   # [n_obs, 4]
            self.axis_resid_valid = _rval_d[_invd]                                          # [n_obs, 4]
        except Exception:
            _n = len(self.axis_braak_high)
            self.axis_resid = np.zeros((_n, 4), np.float32)
            self.axis_resid_valid = np.zeros((_n, 4), bool)

        # Per-cell SEX id (0=female, 1=male) for the OPTIONAL sex-alignment loss
        # term (grouped by sex_id, reference-only, same EMA covariance engine as
        # celltype). Read straight from the original zarr obs column
        # ``axis_sex_male`` (0/1; in the SEA-AD DLPFC+MTG obs: female=0, male=1).
        # Cells lacking a usable 0/1 value get the sentinel -1, which the loss's
        # grouping simply treats as its own (ignored) group rather than crashing —
        # this matches how the dataset never breaks on missing covariates. Aligned
        # to FULL obs rows (indexed by row_index), exactly like the axes above.
        self.sex_ids = _read_sex_male_ids(obs)
        self.age_years, self.age_valid = _read_age_years(obs)

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

        self.return_thinned = bool(return_thinned)
        if not (0.0 < float(thin_rho) < 1.0):
            raise ValueError("thin_rho must be in (0, 1).")
        self.thin_rho = float(thin_rho)
        # Lazily-created per-worker RNG.  It is seeded from PyTorch's worker
        # seed instead of fresh OS entropy, so a fixed training seed,
        # world-size, sampler order and worker count reproduce the exact
        # thinning stream across runs while persistent workers still produce
        # a new draw for successive examples/epochs.
        self._thin_rng = None

        # Raw-count companion for faithful CP10K-scale thinning (see signature).
        self.thin_counts_zarr_path = (
            str(thin_counts_zarr_path) if thin_counts_zarr_path is not None else None
        )
        self.thin_library_size_key = str(thin_library_size_key)
        self.X_thin_counts: Optional[ZarrMatrixReader] = None
        self.thin_library_size: Optional[np.ndarray] = None

        # Binomial thinning is only valid on integer molecule counts.
        #   - companion given: thin the row-aligned raw companion, then CP10K-
        #     renormalize by the full library size to match the training scale.
        #   - companion None: legacy path — thin the main matrix directly, which
        #     requires it to be integer (the `auto` resolver can fall back to a
        #     cp10k `X`; np.rint of such values fabricates counts). Fail loud.
        if self.return_thinned:
            if self.thin_counts_zarr_path is not None:
                self._setup_thin_counts_companion(
                    obs=obs,
                    thin_counts_matrix_path=str(thin_counts_matrix_path),
                    max_csr_data_span=max_csr_data_span,
                    max_span_to_selected_nnz_ratio=max_span_to_selected_nnz_ratio,
                )
            else:
                self._assert_integer_counts_for_thinning()

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

        # Per-cell DONOR index (v20 donor-balanced pathology aggregation). Distinct
        # from batch_key (= technical batch tech_id_model): donors are the effective
        # sample unit for the donor-level pathology labels. Built over the FULL obs
        # (global, split-independent) so train/val share one donor->id map.
        _donor_vals = obs["donor_id"].astype(str).values
        self.donor_vocab = np.asarray(sorted(set(_donor_vals.tolist())), dtype=object)
        _donor_to_id = {str(d): i for i, d in enumerate(self.donor_vocab.tolist())}
        self.donor_ids = np.array(
            [_donor_to_id[str(v)] for v in _donor_vals], dtype=np.int64
        )
        self.n_donor = len(_donor_to_id)

        # Age is normalised using each TRAIN donor exactly once.  Computing the
        # mean over cells would let donors with more captured cells dominate the
        # nuisance reference and would leak held-out donor ages into the scale.
        donor_age = np.full(self.n_donor, np.nan, dtype=np.float64)
        for donor in range(self.n_donor):
            rows = np.flatnonzero((self.donor_ids == donor) & self.age_valid)
            if rows.size:
                values = self.age_years[rows].astype(np.float64)
                if float(np.ptp(values)) > 1e-6:
                    raise ValueError(f"age varies within donor_id={self.donor_vocab[donor]}")
                donor_age[donor] = float(values[0])
        self.donor_age_years = donor_age.astype(np.float32)
        train_rows_all = obs[split_key].astype(str).values == "train"
        train_donor_ids = np.unique(self.donor_ids[train_rows_all])
        train_age = donor_age[train_donor_ids]
        train_age = train_age[np.isfinite(train_age)]
        if train_age.size < 2:
            raise ValueError("PRISM E2E requires age for at least two training donors")
        self.age_train_mean = float(np.mean(train_age))
        self.age_train_std = float(np.std(train_age))
        if not np.isfinite(self.age_train_std) or self.age_train_std < 1e-6:
            self.age_train_std = 1.0
        self.age_z = np.where(
            self.age_valid,
            (self.age_years - self.age_train_mean) / self.age_train_std,
            0.0,
        ).astype(np.float32)

        # Region is an explicit biological context, not folded into a technical
        # batch.  Prefer the harmonised brain_region label used by the combined
        # DLPFC/MTG artifact.
        _region_key = "brain_region" if "brain_region" in obs.columns else "concat_region"
        if _region_key not in obs.columns:
            raise KeyError("PRISM E2E requires 'brain_region' or 'concat_region' in obs")
        _region_vals = obs[_region_key].astype(str).values
        self.region_vocab = np.asarray(sorted(set(_region_vals.tolist())), dtype=object)
        _region_to_id = {str(v): i for i, v in enumerate(self.region_vocab.tolist())}
        self.region_ids = np.asarray(
            [_region_to_id[str(v)] for v in _region_vals], dtype=np.int64
        )
        self.n_region = len(self.region_vocab)

        # A per-cell inverse-frequency weight makes any scalar marked as
        # donor-balanced give every donor equal total mass inside this split.
        split_donor = self.donor_ids[self.row_idx]
        split_counts = np.bincount(split_donor, minlength=self.n_donor).astype(np.float64)
        n_split_donors = int(np.sum(split_counts > 0))
        self.donor_balance_weight = np.zeros(len(obs), dtype=np.float32)
        if n_split_donors > 0:
            present = split_counts > 0
            per_donor = np.zeros_like(split_counts)
            per_donor[present] = len(self.row_idx) / (
                float(n_split_donors) * split_counts[present]
            )
            self.donor_balance_weight[self.row_idx] = per_donor[split_donor].astype(np.float32)

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

    def _assert_integer_counts_for_thinning(
        self, n_probe: int = 16, reader: Optional["ZarrMatrixReader"] = None
    ) -> None:
        """Reject non-integer matrices when thinning is requested.

        Faithful binomial thinning (``np.rint`` → binomial) requires raw molecule
        counts. Probe a few rows of ``reader`` (the raw companion when given, else
        the main matrix) and raise if any value is negative or fractional, so a
        missing/normalized counts source surfaces at construction instead of
        silently corrupting the thinned view.
        """
        rdr = reader if reader is not None else self.X
        src = rdr.name
        n = len(self.row_idx)
        if n == 0:
            return
        k = int(min(n_probe, n))
        probe = np.unique(self.row_idx[np.linspace(0, n - 1, num=k).astype(np.int64)])
        max_frac = 0.0
        min_val = 0.0
        for r in probe:
            x = rdr.row_to_dense(int(r)).astype(np.float32, copy=False)
            if x.size == 0:
                continue
            min_val = min(min_val, float(x.min()))
            max_frac = max(max_frac, float(np.abs(x - np.rint(x)).max()))
        tol = 1e-4
        if min_val < -tol or max_frac > tol:
            raise ValueError(
                f"return_thinned=True requires integer raw counts, but the resolved "
                f"matrix '{src}' is not integer-valued "
                f"(min={min_val:.4g}, max|x-rint(x)|={max_frac:.4g} over {k} probed "
                f"rows). Binomial thinning of non-count (e.g. cp10k-normalized) "
                f"values is invalid. Point the dataset at an integer counts layer / "
                f"raw matrix (counts_layer, prefer_raw, or an explicit matrix_path / "
                f"thin_counts_zarr_path), or set return_thinned=False."
            )

    def _setup_thin_counts_companion(
        self,
        *,
        obs: pd.DataFrame,
        thin_counts_matrix_path: str,
        max_csr_data_span: int,
        max_span_to_selected_nnz_ratio: float,
    ) -> None:
        """Open the row/gene-aligned raw-count companion store used to rebuild
        the thinned view on the training (log1p CP10K) scale, and cache the
        per-cell full library size used as the renormalization denominator."""
        comp_root = zarr.open(self.thin_counts_zarr_path, mode="r")
        matrix, comp_src, comp_var_key = choose_zarr_matrix(
            comp_root, matrix_path=thin_counts_matrix_path, prefer_raw=True
        )
        comp_var = read_zarr_dataframe(comp_root, comp_var_key)
        comp_var_names = np.asarray(comp_var.index.astype(str), dtype=object)
        if not np.array_equal(comp_var_names, np.asarray(self.spec.gene_names, dtype=object)):
            raise ValueError(
                "Thinning counts companion gene order/identity does not match the "
                f"ordinal spec (companion '{self.thin_counts_zarr_path}')."
            )
        reader = ZarrMatrixReader(
            matrix,
            n_cols=len(self.spec.gene_names),
            name=f"thin_counts:{comp_src}",
            max_csr_data_span=max_csr_data_span,
            max_span_to_selected_nnz_ratio=max_span_to_selected_nnz_ratio,
        )
        if reader.shape[0] != self.X.shape[0]:
            raise ValueError(
                f"Thinning companion row count {reader.shape[0]} != main matrix row "
                f"count {self.X.shape[0]}; the companion must be row-aligned (same "
                f"cells in the same order)."
            )
        # Hard row-IDENTITY guard. Same row COUNT is not enough: if the cell ORDER
        # differs, thinning silently pairs the wrong cells and corrupts everything.
        # Require the companion obs index to be identical AND in the same order.
        comp_obs = read_zarr_dataframe(comp_root, "obs")
        comp_index = np.asarray(comp_obs.index.astype(str), dtype=object)
        main_index = np.asarray(obs.index.astype(str), dtype=object)
        if comp_index.shape[0] != main_index.shape[0] or not np.array_equal(comp_index, main_index):
            raise ValueError(
                "Thinning companion obs index does not match the main obs index in "
                "identity/order; the companion must contain the SAME cells in the "
                f"SAME order (companion '{self.thin_counts_zarr_path}'). A matching "
                "row count is NOT sufficient — mismatched order silently corrupts "
                "the thinned view."
            )
        self.X_thin_counts = reader
        # Companion must be integer molecule counts.
        self._assert_integer_counts_for_thinning(reader=reader)

        # Per-cell full library size = the CP10K normalization denominator the
        # training matrix used (verified: cp10k = raw / NumUMIs * 1e4). NOT the
        # subset-gene sum. Used to renormalize the thinned view on the same basis.
        key = self.thin_library_size_key
        if key not in obs.columns:
            raise KeyError(
                f"thin_library_size_key '{key}' not found in obs; it is required as "
                f"the full-library CP10K renorm denominator for the thinned view."
            )
        lib = np.asarray(pd.to_numeric(obs[key], errors="coerce").to_numpy(), dtype=np.float64)
        if lib.shape[0] != self.X.shape[0]:
            raise ValueError(
                f"library size column '{key}' length {lib.shape[0]} != matrix rows "
                f"{self.X.shape[0]}."
            )
        self.thin_library_size = lib

        # Hard guard (early): the full library size must be finite and >= the
        # subset raw sum for every cell (subset genes ⊆ full library). A violation
        # signals a library-size / alignment / source mismatch and must NOT be
        # silently repaired. Probe a few rows at construction to fail fast.
        n = len(self.row_idx)
        if n > 0:
            k = int(min(16, n))
            probe = np.unique(self.row_idx[np.linspace(0, n - 1, num=k).astype(np.int64)])
            for r in probe:
                subset_sum = float(np.rint(np.clip(reader.row_to_dense(int(r)), 0.0, None)).sum())
                lv = float(lib[int(r)])
                if not np.isfinite(lv) or lv < subset_sum - 1e-6:
                    raise ValueError(
                        f"thin_library_size_key '{key}'={lv} at row {int(r)} is non-finite "
                        f"or below the subset raw sum {subset_sum}; this indicates a "
                        f"library-size / alignment / source mismatch. Refusing to "
                        f"proceed (the thinned-view CP10K renorm denominator would be wrong)."
                    )

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
            "donor_id": torch.tensor(int(self.donor_ids[r]), dtype=torch.long),
            "region_id": torch.tensor(int(self.region_ids[r]), dtype=torch.long),
            "donor_balance_weight": torch.tensor(
                float(self.donor_balance_weight[r]), dtype=torch.float32
            ),
            "age_z": torch.tensor(float(self.age_z[r]), dtype=torch.float32),
            "age_valid": torch.tensor(bool(self.age_valid[r]), dtype=torch.bool),
            "row_index": torch.tensor(r, dtype=torch.long),
        }

        item.update(
            {
                "is_reference": torch.tensor(bool(self.is_reference_origin[r]), dtype=torch.bool),
                "is_abnormal": torch.tensor(bool(self.is_abnormal_pathology[r]), dtype=torch.bool),

                # Per-cell sex id (0/1; -1 sentinel = unknown, ignored by the loss)
                # for the optional reference-only sex-alignment term.
                "sex_id": torch.tensor(int(self.sex_ids[r]), dtype=torch.long),

                "axis_late_positive": torch.tensor(bool(self.axis_late_positive[r]), dtype=torch.bool),
                "axis_cognitive_impaired": torch.tensor(bool(self.axis_cognitive_impaired[r]), dtype=torch.bool),
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
                # v23 residualized pathology targets (donor-constant) + validity mask
                "axis_resid": torch.from_numpy(self.axis_resid[r].copy()),                 # [4] float32
                "axis_resid_valid": torch.from_numpy(self.axis_resid_valid[r].copy()),     # [4] bool
                # PRISM E2E named pathology contract, order:
                # Thal, Braak, CERAD, LATE, Lewy.  Missing is represented only
                # by the validity mask; ADNC is never an input.
                "prism_pathology": torch.from_numpy(self.prism_pathology[r].copy()),       # [5]
                "prism_pathology_valid": torch.from_numpy(
                    self.prism_pathology_valid[r].copy()
                ),                                                                           # [5]
            }
        )
        if self.return_log_input:
            item["x_log1p"] = torch.from_numpy(x_log)
        if self.return_x_gene_scalar:
            item["x_gene_scalar"] = torch.from_numpy(self._make_x_gene_scalar(x_log, celltype_id))

        if self.return_thinned:
            if self._thin_rng is None:
                # DataLoader deterministically assigns torch.initial_seed() to
                # each worker. In the num_workers=0 path this is the process
                # seed set by the runner. Mask to uint64 for numpy portability.
                self._thin_rng = np.random.default_rng(
                    int(torch.initial_seed()) & ((1 << 64) - 1)
                )
            rng = self._thin_rng

            if self.X_thin_counts is not None:
                # Faithful CP10K-scale thinning: thin the row-aligned RAW counts,
                # then renormalize by the *thinned full library size* (same basis
                # as training: cp10k = raw / NumUMIs * 1e4), log1p, and re-bin
                # with the same edges. A gene nonzero in full but thinned to 0 is
                # a proven technical zero. See thinning_stageA doc.
                raw_full = self.X_thin_counts.row_to_dense(r)
                raw_int = np.rint(np.clip(raw_full, 0.0, None)).astype(np.int64)
                raw_thin = rng.binomial(raw_int, self.thin_rho)
                subset_sum = int(raw_int.sum())

                numi = (
                    float(self.thin_library_size[r])
                    if self.thin_library_size is not None
                    else float(subset_sum)
                )
                # Hard fail (no silent repair): full library < subset raw sum (or
                # non-finite) is impossible for consistent data (subset ⊆ full) and
                # signals an alignment / metadata / source mismatch. Repairing it
                # would feed a wrong CP10K renorm denominator and corrupt results.
                if not math.isfinite(numi) or numi < subset_sum - 1e-6:
                    raise ValueError(
                        f"thinning: library size ({numi}) is non-finite or below the "
                        f"subset raw sum ({subset_sum}) at row {r}; indicates a "
                        f"library-size / alignment / source mismatch."
                    )
                # Thin the non-subset molecules too, so the denominator is the
                # thinned FULL library (subset thinned + non-subset thinned).
                other_full = max(int(round(numi)) - subset_sum, 0)
                other_thin = int(rng.binomial(other_full, self.thin_rho)) if other_full > 0 else 0
                numi_thin = float(int(raw_thin.sum()) + other_thin)
                denom = numi_thin if numi_thin > 0.0 else 1.0

                cp10k_thin = (raw_thin.astype(np.float32) / denom) * 1e4
                x_log_thin = np.log1p(cp10k_thin).astype(np.float32, copy=False)
            else:
                # Legacy path: thin the main matrix directly (valid only when it
                # is raw counts and the spec is on log1p(raw)). No CP10K renorm.
                counts_int = np.rint(np.clip(x_counts, 0.0, None)).astype(np.int64)
                counts_thin = rng.binomial(counts_int, self.thin_rho).astype(np.float32)
                x_log_thin = np.log1p(counts_thin).astype(np.float32, copy=False)

            y_ord_thin = self._bin_one_cell_from_log(x_log_thin)
            item["y_ord_thin"] = torch.from_numpy(y_ord_thin)
            if self.return_log_input:
                item["x_log1p_thin"] = torch.from_numpy(x_log_thin)
            if self.return_x_gene_scalar:
                item["x_gene_scalar_thin"] = torch.from_numpy(
                    self._make_x_gene_scalar(x_log_thin, celltype_id)
                )
        return item


if __name__ == "__main__":
    print("This module defines Zarr-backed OrdinalScDataset.")
