# src/data/pathology_masks.py
from __future__ import annotations

from typing import Any
import numpy as np
import pandas as pd


# =============================================================================
# Small helpers
# =============================================================================
def _s(x: Any) -> str:
    if pd.isna(x):
        return "__NA__"
    return str(x).strip()


def _lower(x: Any) -> str:
    return _s(x).lower()


def _boolish_true(x: Any) -> bool:
    s = _lower(x)
    return s in {"true", "yes", "y", "1", "reference"}


def _col(obs: pd.DataFrame, name: str, default: str = "__NA__") -> pd.Series:
    if name in obs.columns:
        return obs[name].astype(str)
    return pd.Series([default] * len(obs), index=obs.index, dtype="object")


def _axis_zero(obs: pd.DataFrame, name: str) -> np.ndarray:
    """
    Return True where an existing 0/1 axis is explicitly zero.

    For strict reference, missing or non-numeric values are treated as NOT clean.
    """
    if name not in obs.columns:
        return np.zeros(len(obs), dtype=bool)

    x = pd.to_numeric(obs[name], errors="coerce")
    return x.fillna(1).astype(int).eq(0).to_numpy()


def _axis_one(obs: pd.DataFrame, name: str) -> np.ndarray:
    """
    Return True where an existing 0/1 axis is explicitly one.

    Missing values are treated as negative for positive pathology axes.
    """
    if name not in obs.columns:
        return np.zeros(len(obs), dtype=bool)

    x = pd.to_numeric(obs[name], errors="coerce")
    return x.fillna(0).astype(int).eq(1).to_numpy()


def _isin(obs: pd.DataFrame, name: str, allowed: set[str]) -> np.ndarray:
    if name not in obs.columns:
        return np.zeros(len(obs), dtype=bool)
    return obs[name].astype(str).isin(allowed).to_numpy()


def _not_in_known_clean(obs: pd.DataFrame, name: str, clean: set[str], positive: set[str]) -> np.ndarray:
    """
    Positive only if value is explicitly in the positive set.

    Unknown/unclassifiable labels are not treated as positive, but they also
    should not pass strict reference-clean checks.
    """
    if name not in obs.columns:
        return np.zeros(len(obs), dtype=bool)
    return obs[name].astype(str).isin(positive).to_numpy()


# =============================================================================
# Main mask builder
# =============================================================================
def build_pathology_masks(
    obs: pd.DataFrame,
    *,
    exclude_apoe4_from_reference: bool = False,
) -> dict[str, np.ndarray]:
    """
    Build precision-medicine-aware pathology masks.

    Main outputs
    ------------
    is_reference_origin:
        Strict pathology-clean reference cells used to anchor each Subclass
        latent universe at zero.

    is_abnormal_pathology:
        Broad pathology/disease-positive cells used for residual/pathology
        analysis.

    Philosophy
    ----------
    Reference origin should be stricter than "clinical normal".
    In particular, CERAD Sparse, LATE positive, Lewy positive, and
    microinfarct positive cells should not define the normal origin.

    APOE4 is treated as a risk / precision-medicine axis, not necessarily
    as current tissue pathology. It is therefore not excluded from reference
    by default, but can be excluded via exclude_apoe4_from_reference=True.
    """

    n = len(obs)

    # -------------------------------------------------------------------------
    # Clinical / disease labels
    # -------------------------------------------------------------------------
    disease = _col(obs, "disease").str.lower()
    disease_normal = disease.eq("normal").to_numpy()
    disease_abnormal = (~disease.eq("normal")).to_numpy() if "disease" in obs.columns else np.zeros(n, dtype=bool)

    if "Neurotypical reference" in obs.columns:
        neuro_ref = obs["Neurotypical reference"].map(_boolish_true).to_numpy()
    else:
        neuro_ref = np.zeros(n, dtype=bool)

    # Require at least one clinical/reference-like indicator.
    # Do NOT use axis_control_like alone; it is intentionally broad.
    clinical_reference = disease_normal | neuro_ref

    # -------------------------------------------------------------------------
    # Raw pathology clean definitions
    # -------------------------------------------------------------------------
    adnc_clean = _isin(obs, "ADNC", {"Not AD", "Reference"})

    braak_clean = _isin(obs, "Braak stage", {"Braak 0", "Reference", "0"})
    thal_clean = _isin(obs, "Thal phase", {"Thal 0", "Reference", "0"})

    # IMPORTANT:
    # CERAD Sparse is NOT reference-clean. It is low plaque burden, but still
    # pathology-positive for our reference-origin purpose.
    cerad_clean = _isin(obs, "CERAD score", {"Absent", "Reference"})

    late_clean = _isin(obs, "LATE-NC stage", {"Not Identified", "Reference"})

    # -------------------------------------------------------------------------
    # Axis-clean definitions, using explicit axis columns when available
    # -------------------------------------------------------------------------
    thal_axis_clean = (
        _axis_zero(obs, "axis_Thal_positive")
        if "axis_Thal_positive" in obs.columns
        else thal_clean
    )

    braak_axis_clean = (
        _axis_zero(obs, "axis_Braak_mid")
        if "axis_Braak_mid" in obs.columns
        else braak_clean
    )

    late_axis_clean = (
        _axis_zero(obs, "axis_LATE_positive")
        if "axis_LATE_positive" in obs.columns
        else late_clean
    )

    lewy_axis_clean = (
        _axis_zero(obs, "axis_Lewy_positive")
        if "axis_Lewy_positive" in obs.columns
        else np.ones(n, dtype=bool)
    )

    micro_axis_clean = (
        _axis_zero(obs, "axis_Microinfarct_positive")
        if "axis_Microinfarct_positive" in obs.columns
        else np.ones(n, dtype=bool)
    )

    # Your current obs may not contain CERAD axis columns.
    # If not present, raw CERAD labels control reference cleanliness.
    if "axis_CERAD_sparse_or_higher" in obs.columns:
        cerad_axis_clean = _axis_zero(obs, "axis_CERAD_sparse_or_higher")
    elif "axis_CERAD_moderate_or_frequent" in obs.columns:
        # This is weaker because Sparse would not be removed.
        # Combine with raw cerad_clean below, so Sparse is still excluded.
        cerad_axis_clean = _axis_zero(obs, "axis_CERAD_moderate_or_frequent")
    else:
        cerad_axis_clean = cerad_clean

    if "axis_ADNC_intermediate_or_high_report" in obs.columns:
        adnc_axis_clean = _axis_zero(obs, "axis_ADNC_intermediate_or_high_report")
    else:
        adnc_axis_clean = adnc_clean

    # -------------------------------------------------------------------------
    # Strict reference-origin mask
    # -------------------------------------------------------------------------
    reference = (
        clinical_reference
        & adnc_clean
        & adnc_axis_clean
        & braak_clean
        & braak_axis_clean
        & thal_clean
        & thal_axis_clean
        & cerad_clean
        & cerad_axis_clean
        & late_clean
        & late_axis_clean
        & lewy_axis_clean
        & micro_axis_clean
    )

    if exclude_apoe4_from_reference and "axis_APOE4_carrier" in obs.columns:
        reference &= _axis_zero(obs, "axis_APOE4_carrier")

    # -------------------------------------------------------------------------
    # Positive pathology axes
    # -------------------------------------------------------------------------
    thal_pos = (
        _axis_one(obs, "axis_Thal_positive")
        if "axis_Thal_positive" in obs.columns
        else _isin(obs, "Thal phase", {"Thal 1", "Thal 2", "Thal 3", "Thal 4", "Thal 5"})
    )
    thal_high = (
        _axis_one(obs, "axis_Thal_high")
        if "axis_Thal_high" in obs.columns
        else _isin(obs, "Thal phase", {"Thal 4", "Thal 5"})
    )

    braak_mid = (
        _axis_one(obs, "axis_Braak_mid")
        if "axis_Braak_mid" in obs.columns
        else _isin(obs, "Braak stage", {"Braak III", "Braak IV", "Braak V", "Braak VI"})
    )
    braak_high = (
        _axis_one(obs, "axis_Braak_high")
        if "axis_Braak_high" in obs.columns
        else _isin(obs, "Braak stage", {"Braak V", "Braak VI"})
    )

    # CERAD:
    # - positive: Sparse / Moderate / Frequent
    # - moderate_or_frequent: Moderate / Frequent
    # - frequent: Frequent
    cerad_positive = _isin(obs, "CERAD score", {"Sparse", "Moderate", "Frequent"})
    cerad_moderate_or_frequent = _isin(obs, "CERAD score", {"Moderate", "Frequent"})
    cerad_frequent = _isin(obs, "CERAD score", {"Frequent"})

    if "axis_CERAD_sparse_or_higher" in obs.columns:
        cerad_positive = _axis_one(obs, "axis_CERAD_sparse_or_higher")

    if "axis_CERAD_moderate_or_frequent" in obs.columns:
        cerad_moderate_or_frequent = _axis_one(obs, "axis_CERAD_moderate_or_frequent")

    if "axis_CERAD_frequent" in obs.columns:
        cerad_frequent = _axis_one(obs, "axis_CERAD_frequent")

    adnc_intermediate_high = (
        _axis_one(obs, "axis_ADNC_intermediate_or_high_report")
        | _isin(obs, "ADNC", {"Intermediate", "High"})
    )
    adnc_high = (
        _axis_one(obs, "axis_ADNC_high_report")
        | _isin(obs, "ADNC", {"High"})
    )

    late_pos = (
        _axis_one(obs, "axis_LATE_positive")
        | _isin(obs, "LATE-NC stage", {"LATE Stage 1", "LATE Stage 2", "LATE Stage 3"})
    )
    late_advanced = (
        _axis_one(obs, "axis_LATE_advanced")
        | _isin(obs, "LATE-NC stage", {"LATE Stage 2", "LATE Stage 3"})
    )

    lewy_pos = _axis_one(obs, "axis_Lewy_positive")
    micro_pos = _axis_one(obs, "axis_Microinfarct_positive")
    apoe4 = _axis_one(obs, "axis_APOE4_carrier")

    abnormal_pathology = (
        disease_abnormal
        | thal_pos
        | braak_mid
        | cerad_positive
        | adnc_intermediate_high
        | late_pos
        | lewy_pos
        | micro_pos
    )

    return {
        # Main masks
        "is_reference_origin": reference.astype(bool),
        "is_abnormal_pathology": abnormal_pathology.astype(bool),

        # Amyloid / tau / plaque axes
        "axis_thal_positive": thal_pos.astype(bool),
        "axis_thal_high": thal_high.astype(bool),
        "axis_braak_mid": braak_mid.astype(bool),
        "axis_braak_high": braak_high.astype(bool),

        # CERAD axes
        "axis_cerad_positive": cerad_positive.astype(bool),
        "axis_cerad_moderate_or_frequent": cerad_moderate_or_frequent.astype(bool),
        "axis_cerad_high": cerad_moderate_or_frequent.astype(bool),  # backward-compatible alias
        "axis_cerad_frequent": cerad_frequent.astype(bool),

        # ADNC report axes
        "axis_adnc_intermediate_high": adnc_intermediate_high.astype(bool),
        "axis_adnc_high": adnc_high.astype(bool),

        # Precision-medicine-relevant co-pathology axes
        "axis_late_positive": late_pos.astype(bool),
        "axis_late_advanced": late_advanced.astype(bool),
        "axis_lewy_positive": lewy_pos.astype(bool),
        "axis_microinfarct_positive": micro_pos.astype(bool),
        "axis_apoe4_carrier": apoe4.astype(bool),
    }