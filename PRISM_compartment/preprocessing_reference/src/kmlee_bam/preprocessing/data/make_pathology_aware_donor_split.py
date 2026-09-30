from __future__ import annotations

"""
make_pathology_aware_donor_split_v4.py
======================================

Hybrid donor-level split generator for SEA-AD-like sc/snRNA cohorts.

This v4 combines the strongest parts of the previous proposals:

Taken from the stricter v2-style proposal
-----------------------------------------
1. Explicit ordinal mappings for Thal, Braak, CERAD, ADNC, LATE.
2. Correct APOE4 parser for both Y/N and genotype-like strings such as 3/4.
3. Iterative multi-label stratification when `iterative-stratification` exists.
4. Locked final hold-out test set plus K-fold validation folds.
5. Optional amyloid-axis choice: Thal, CERAD, or both.

Taken from the practical v3-style proposal
------------------------------------------
1. ADNC and Cognitive status are report-only, not split-objective axes.
2. Braak/Thal/CERAD/LATE/Lewy/Microinfarct are the primary pathology axes.
3. APOE4 is treated as a strong precision-medicine/risk stratification axis.
4. Cell-count balance is considered when choosing the hold-out fold.
5. Obs sidecar is written in Parquet by default and keeps donor-level pathology
   columns for downstream probing.

Key design decision
-------------------
ADNC is NOT used for stratification if Thal/Braak/CERAD are available, because
ADNC is a composite of those components. Cognitive status is also NOT used for
stratification because it is a clinical phenotype/outcome, not a primary
biological pathology axis.

Outputs
-------
out_dir/
  donor_split_table.csv/parquet
  donor_label_matrix.csv/parquet
  obs_with_pathology_split.parquet
  pathology_balance_report.csv
  pathology_group_counts_by_split.csv
  ordinal_mapping_used.json
  split_manifest.json

Typical use
-----------
python make_pathology_aware_donor_split_v4.py --config pathology_split_v4_cfg.json
"""

import argparse
import json
import re
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd


# =============================================================================
# 1. Explicit SEA-AD-compatible mappings
# =============================================================================
THAL_ORDINAL: dict[str, int] = {
    "Reference": 0,
    "Thal 0": 0,
    "Thal 1": 1,
    "Thal 2": 2,
    "Thal 3": 3,
    "Thal 4": 4,
    "Thal 5": 5,
    "0": 0,
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5": 5,
}

BRAAK_ORDINAL: dict[str, int] = {
    "Reference": 0,
    "Braak 0": 0,
    "Braak I": 1,
    "Braak II": 2,
    "Braak III": 3,
    "Braak IV": 4,
    "Braak V": 5,
    "Braak VI": 6,
    "0": 0,
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
}

CERAD_ORDINAL: dict[str, int] = {
    "Reference": 0,
    "Absent": 0,
    "Sparse": 1,
    "Moderate": 2,
    "Frequent": 3,
}

ADNC_ORDINAL: dict[str, int] = {
    "Reference": 0,
    "Not AD": 0,
    "Low": 1,
    "Intermediate": 2,
    "High": 3,
}

LATE_ORDINAL: dict[str, int] = {
    "Reference": 0,
    "Not Identified": 0,
    "LATE Stage 1": 1,
    "LATE Stage 2": 2,
    "LATE Stage 3": 3,
    "Stage 1": 1,
    "Stage 2": 2,
    "Stage 3": 3,
    # Unclassifiable is not negative; it is unknown. We parse it as None below.
}

LEWY_NEGATIVE_PATTERNS: tuple[str, ...] = (
    "not identified",
    "reference",
)

LEWY_POSITIVE_PATTERNS: tuple[str, ...] = (
    "amygdala",
    "olfactory bulb only",
    "brainstem",
    "limbic",
    "transitional",
    "neocortical",
    "diffuse",
)

MICROINFARCT_NEGATIVE_PATTERNS: tuple[str, ...] = (
    "reference",
    "0 to 3 microinfarcts",
    "0-3 microinfarcts",
    "0 to 3",
    "not identified",
    "absent",
)

MICROINFARCT_POSITIVE_PATTERNS: tuple[str, ...] = (
    "4 to 6 microinfarcts",
    "7 to 10 microinfarcts",
    "4-6 microinfarcts",
    "7-10 microinfarcts",
    "4 to 6",
    "7 to 10",
    "present",
)


# =============================================================================
# 2. Config
# =============================================================================
@dataclass
class PathologyAwareSplitConfigV4:
    donor_pathology_table_path: str

    # Exactly one obs source is required. zarr_path is usually preferred.
    zarr_path: Optional[str] = None
    obs_h5ad_path: Optional[str] = None
    obs_csv_path: Optional[str] = None
    obs_parquet_path: Optional[str] = None

    out_dir: str = "pathology_aware_split_v4"

    # Core schema.
    donor_key: str = "donor_id"
    split_key: str = "split"
    fold_key: str = "fold_id"
    celltype_key: str = "Subclass"

    # Donor pathology table columns.
    col_thal: str = "Thal phase"
    col_braak: str = "Braak stage"
    col_cerad: str = "CERAD score"
    col_adnc: str = "ADNC"
    col_late: str = "LATE-NC stage"
    col_lewy: str = "Lewy body disease pathology"
    col_microinfarct: str = "Microinfarct pathology"
    col_apoe4: str = "APOE4 status"
    col_disease: str = "disease"
    col_cognitive: str = "Cognitive status"
    col_neurotypical: str = "Neurotypical reference"
    col_sex: str = "sex"
    col_age: str = "Age at death"

    # Split structure.
    n_folds: int = 5
    final_holdout_frac: float = 0.10
    val_fold_index: int = 0
    seed: int = 42

    # Stratification policy.
    # "thal" is usually preferred because Thal is an amyloid-distribution axis.
    # "cerad" uses neuritic plaque burden instead.
    # "both" is allowed but more constraining with only ~83 donors.
    amyloid_axis_choice: str = "thal"  # thal | cerad | both

    # Thresholds.
    thal_positive_threshold: int = 1
    thal_high_threshold: int = 4
    braak_mid_threshold: int = 3
    braak_high_threshold: int = 5
    cerad_moderate_threshold: int = 2
    cerad_frequent_threshold: int = 3
    late_positive_threshold: int = 1
    late_advanced_threshold: int = 2
    age_band_threshold: float = 85.0

    # Parser behaviour.
    strict_ordinal: bool = False
    use_demographics_as_weak_balance: bool = True

    # Cell count balance.
    # The iterative stratifier balances donors/labels. This coefficient is used
    # only when selecting which candidate fold becomes the final hold-out.
    holdout_cell_count_weight: float = 0.25

    # Technical covariate for downstream model batch/tech id.
    tech_id_model_key: str = "tech_id_model"
    tech_source_candidates: list[str] = field(
        default_factory=lambda: ["assay", "assay_ontology_term_id", "PMI"]
    )
    tech_missing_value: str = "__GLOBAL__"

    # Output.
    save_csv: bool = True
    save_parquet: bool = True
    save_obs_csv: bool = False
    save_obs_parquet: bool = True


# =============================================================================
# 3. IO
# =============================================================================
def _log(msg: str) -> None:
    print(f"[pathology-split-v4] {msg}", flush=True)


def load_config(path: str | Path) -> PathologyAwareSplitConfigV4:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return PathologyAwareSplitConfigV4(**raw)


def _read_zarr_obs(zarr_path: str | Path) -> pd.DataFrame:
    import zarr

    root = zarr.open(str(zarr_path), mode="r")
    try:
        from anndata.experimental import read_elem  # type: ignore
        obs = read_elem(root["obs"])
        if not isinstance(obs, pd.DataFrame):
            raise TypeError("Zarr obs is not a pandas DataFrame.")
        return obs.copy()
    except Exception:
        import anndata as ad
        return ad.read_zarr(str(zarr_path)).obs.copy()


def load_obs(cfg: PathologyAwareSplitConfigV4) -> pd.DataFrame:
    n_src = sum(
        x is not None
        for x in [cfg.zarr_path, cfg.obs_h5ad_path, cfg.obs_csv_path, cfg.obs_parquet_path]
    )
    if n_src != 1:
        raise ValueError(
            "Provide exactly one obs source: zarr_path, obs_h5ad_path, obs_csv_path, or obs_parquet_path."
        )

    if cfg.zarr_path is not None:
        _log(f"loading obs from zarr: {cfg.zarr_path}")
        return _read_zarr_obs(cfg.zarr_path)

    if cfg.obs_h5ad_path is not None:
        import anndata as ad
        _log(f"loading obs from h5ad: {cfg.obs_h5ad_path}")
        return ad.read_h5ad(cfg.obs_h5ad_path, backed="r").obs.copy()

    if cfg.obs_csv_path is not None:
        _log(f"loading obs from csv: {cfg.obs_csv_path}")
        return pd.read_csv(cfg.obs_csv_path, index_col=0)

    if cfg.obs_parquet_path is not None:
        _log(f"loading obs from parquet: {cfg.obs_parquet_path}")
        return pd.read_parquet(cfg.obs_parquet_path)

    raise RuntimeError("unreachable")


def save_df(df: pd.DataFrame, path_stem: Path, cfg: PathologyAwareSplitConfigV4, *, index: bool = True) -> None:
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    if cfg.save_csv:
        df.to_csv(path_stem.with_suffix(".csv"), index=index)
    if cfg.save_parquet:
        try:
            df.to_parquet(path_stem.with_suffix(".parquet"), index=index)
        except Exception as e:
            warnings.warn(f"Could not save parquet for {path_stem.name}: {repr(e)}")


def save_obs(df: pd.DataFrame, path_stem: Path, cfg: PathologyAwareSplitConfigV4) -> None:
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    if cfg.save_obs_parquet:
        try:
            df.to_parquet(path_stem.with_suffix(".parquet"), index=True)
        except Exception as e:
            warnings.warn(f"Could not save obs parquet; falling back to csv. reason={repr(e)}")
            df.to_csv(path_stem.with_suffix(".csv"), index=True)
    if cfg.save_obs_csv:
        df.to_csv(path_stem.with_suffix(".csv"), index=True)


# =============================================================================
# 4. Parsing helpers
# =============================================================================
def normalise_level(x: Any) -> str:
    if pd.isna(x):
        return "__NA__"
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "na", "n/a"}:
        return "__NA__"
    return s


def level_lower(x: Any) -> str:
    return normalise_level(x).lower()


def ordinal_lookup(value: Any, table: dict[str, int], *, axis_name: str, strict: bool) -> Optional[int]:
    s = normalise_level(value)
    if s == "__NA__":
        return None

    if "unclassifiable" in s.lower() or "precluded" in s.lower():
        return None

    if s in table:
        return table[s]

    for k, v in table.items():
        if k.lower() == s.lower():
            return v

    msg = f"Unknown level for {axis_name}: {s!r}. Add it to the mapping if valid."
    if strict:
        raise KeyError(msg)
    warnings.warn(msg, RuntimeWarning)
    return None


def pattern_binary(value: Any, *, positive: tuple[str, ...], negative: tuple[str, ...]) -> Optional[int]:
    s = level_lower(value)
    if s == "__na__":
        return None

    for pat in negative:
        if pat.lower() in s:
            return 0
    for pat in positive:
        if pat.lower() in s:
            return 1
    return None


def apoe4_carrier(value: Any) -> Optional[int]:
    s = normalise_level(value)
    if s == "__NA__":
        return None

    low = s.lower()
    if low in {"reference", "n", "no", "negative", "0"}:
        return 0
    if low in {"y", "yes", "positive", "carrier", "1"}:
        return 1

    # Genotype-like strings: 3/4, 4/4, ε3/ε4, e3/e4.
    cleaned = (
        s.replace("ε", "")
        .replace("E", "")
        .replace("e", "")
        .replace("\\", "/")
        .replace("|", "/")
    )
    parts = [p.strip() for p in cleaned.split("/") if p.strip()]
    if len(parts) == 2:
        try:
            a, b = int(parts[0]), int(parts[1])
            return int(a == 4 or b == 4)
        except ValueError:
            pass

    if "apoe4" in low or "e4" in low:
        return 1

    return None


def age_high(value: Any, threshold: float) -> Optional[int]:
    s = normalise_level(value)
    if s == "__NA__":
        return None

    low = s.lower()
    if "less than" in low:
        nums = re.findall(r"\d+", low)
        if nums:
            return int(float(nums[0]) >= threshold)
        return 0

    if "+" in low:
        nums = re.findall(r"\d+", low)
        if nums:
            return int(float(nums[0]) >= threshold)
        return None

    nums = re.findall(r"\d+", low)
    if len(nums) == 1:
        return int(float(nums[0]) >= threshold)
    if len(nums) >= 2:
        lo, hi = float(nums[0]), float(nums[1])
        mid = (lo + hi) / 2.0
        return int(mid >= threshold)

    return None


def sex_male(value: Any) -> Optional[int]:
    s = level_lower(value)
    if s == "__na__":
        return None
    if s.startswith("m"):
        return 1
    if s.startswith("f"):
        return 0
    return None


def control_like(row: pd.Series, cfg: PathologyAwareSplitConfigV4) -> int:
    disease = level_lower(row.get(cfg.col_disease, "__NA__"))
    neuro = row.get(cfg.col_neurotypical, "__NA__")
    neuro_s = level_lower(neuro)
    adnc = level_lower(row.get(cfg.col_adnc, "__NA__"))

    if disease == "normal":
        return 1
    if neuro_s in {"true", "yes", "1"}:
        return 1
    if "reference" in neuro_s:
        return 1
    if "not ad" in adnc:
        return 1
    return 0


def cognitive_impaired(value: Any) -> Optional[int]:
    """Derived report flag, not a measured cognitive score or dementia-only target.

    SEA-AD 'No dementia' excludes a recorded dementia diagnosis, not all possible
    cognitive impairment. Keep the original label for phenotype selection.
    Historical generated stores can retain the old substring bug even when this
    parser is correct; see repair_cognitive_report_20260921.py for audited repair.
    """
    s = level_lower(value)
    if s == "__na__":
        return None
    # Check explicit negative/reference labels first.  ``"no dementia"``
    # contains the positive substring ``"dementia"`` and would otherwise be
    # misclassified as impaired.
    if "no dementia" in s or "reference" in s:
        return 0
    if "dementia" in s or "mci" in s or "impaired" in s:
        return 1
    return None


# =============================================================================
# 5. Label construction
# =============================================================================
def build_labels(donor_table: pd.DataFrame, cfg: PathologyAwareSplitConfigV4) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    rows = []

    for donor, row in donor_table.iterrows():
        thal = ordinal_lookup(row.get(cfg.col_thal), THAL_ORDINAL, axis_name="Thal", strict=cfg.strict_ordinal)
        braak = ordinal_lookup(row.get(cfg.col_braak), BRAAK_ORDINAL, axis_name="Braak", strict=cfg.strict_ordinal)
        cerad = ordinal_lookup(row.get(cfg.col_cerad), CERAD_ORDINAL, axis_name="CERAD", strict=cfg.strict_ordinal)
        adnc = ordinal_lookup(row.get(cfg.col_adnc), ADNC_ORDINAL, axis_name="ADNC", strict=cfg.strict_ordinal)
        late = ordinal_lookup(row.get(cfg.col_late), LATE_ORDINAL, axis_name="LATE", strict=cfg.strict_ordinal)

        lewy = pattern_binary(
            row.get(cfg.col_lewy),
            positive=LEWY_POSITIVE_PATTERNS,
            negative=LEWY_NEGATIVE_PATTERNS,
        )
        micro = pattern_binary(
            row.get(cfg.col_microinfarct),
            positive=MICROINFARCT_POSITIVE_PATTERNS,
            negative=MICROINFARCT_NEGATIVE_PATTERNS,
        )
        apoe = apoe4_carrier(row.get(cfg.col_apoe4))
        sex_m = sex_male(row.get(cfg.col_sex))
        age_hi = age_high(row.get(cfg.col_age), cfg.age_band_threshold)
        cog_imp = cognitive_impaired(row.get(cfg.col_cognitive))

        d = {
            "donor_id": str(donor),

            # Primary pathology axes.
            "axis_Thal_positive": int(thal is not None and thal >= cfg.thal_positive_threshold),
            "axis_Thal_high": int(thal is not None and thal >= cfg.thal_high_threshold),
            "axis_CERAD_moderate_or_frequent": int(cerad is not None and cerad >= cfg.cerad_moderate_threshold),
            "axis_CERAD_frequent": int(cerad is not None and cerad >= cfg.cerad_frequent_threshold),
            "axis_Braak_mid": int(braak is not None and braak >= cfg.braak_mid_threshold),
            "axis_Braak_high": int(braak is not None and braak >= cfg.braak_high_threshold),
            "axis_LATE_positive": int(late is not None and late >= cfg.late_positive_threshold),
            "axis_LATE_advanced": int(late is not None and late >= cfg.late_advanced_threshold),
            "axis_Lewy_positive": 0 if lewy is None else int(lewy),
            "axis_Microinfarct_positive": 0 if micro is None else int(micro),
            "axis_APOE4_carrier": 0 if apoe is None else int(apoe),

            # Weak/demographic axes.
            "axis_control_like": control_like(row, cfg),
            "axis_sex_male": 0 if sex_m is None else int(sex_m),
            "axis_age_high": 0 if age_hi is None else int(age_hi),

            # Report-only axes.
            "axis_ADNC_high_report": int(adnc is not None and adnc >= 3),
            "axis_ADNC_intermediate_or_high_report": int(adnc is not None and adnc >= 2),
            "axis_ADNC_low_or_notAD_report": int(adnc is not None and adnc <= 1),
            "axis_Cognitive_impaired_report": 0 if cog_imp is None else int(cog_imp),

            # Numeric ordinals for later analysis.
            "ordinal_Thal": -1 if thal is None else int(thal),
            "ordinal_Braak": -1 if braak is None else int(braak),
            "ordinal_CERAD": -1 if cerad is None else int(cerad),
            "ordinal_ADNC": -1 if adnc is None else int(adnc),
            "ordinal_LATE": -1 if late is None else int(late),
        }
        rows.append(d)

    labels = pd.DataFrame(rows).set_index("donor_id")

    primary: list[str] = []
    if cfg.amyloid_axis_choice == "thal":
        primary.extend(["axis_Thal_positive", "axis_Thal_high"])
    elif cfg.amyloid_axis_choice == "cerad":
        primary.extend(["axis_CERAD_moderate_or_frequent", "axis_CERAD_frequent"])
    elif cfg.amyloid_axis_choice == "both":
        primary.extend([
            "axis_Thal_positive",
            "axis_Thal_high",
            "axis_CERAD_moderate_or_frequent",
            "axis_CERAD_frequent",
        ])
    else:
        raise ValueError("amyloid_axis_choice must be one of: thal, cerad, both.")

    primary.extend([
        "axis_Braak_mid",
        "axis_Braak_high",
        "axis_LATE_positive",
        "axis_LATE_advanced",
        "axis_Lewy_positive",
        "axis_Microinfarct_positive",
        "axis_APOE4_carrier",
    ])

    weak = ["axis_control_like"]
    if cfg.use_demographics_as_weak_balance:
        weak.extend(["axis_sex_male", "axis_age_high"])

    report = [
        "axis_ADNC_high_report",
        "axis_ADNC_intermediate_or_high_report",
        "axis_ADNC_low_or_notAD_report",
        "axis_Cognitive_impaired_report",
    ]

    return labels, primary, weak, report


def report_group(row: pd.Series) -> str:
    if int(row.get("axis_control_like", 0)) == 1:
        return "control_like"

    amy = int(row.get("axis_Thal_high", 0)) or int(row.get("axis_CERAD_frequent", 0))
    tau = int(row.get("axis_Braak_high", 0))
    late = int(row.get("axis_LATE_positive", 0))
    late_adv = int(row.get("axis_LATE_advanced", 0))
    lewy = int(row.get("axis_Lewy_positive", 0))
    micro = int(row.get("axis_Microinfarct_positive", 0))

    if (amy or tau) and not late:
        return "AD_pathology_without_LATE"
    if (amy or tau) and late:
        return "AD_pathology_with_LATE"
    if late_adv and not (amy or tau):
        return "LATE_advanced_without_high_AD_pathology"
    if lewy:
        return "Lewy_mixed"
    if micro:
        return "vascular_mixed"
    return "other_mixed_or_unclassified"


# =============================================================================
# 6. Iterative stratification + fallback
# =============================================================================
def iterative_kfold_indices(Y: np.ndarray, n_splits: int, seed: int) -> list[np.ndarray]:
    try:
        from iterstrat.ml_stratifiers import MultilabelStratifiedKFold  # type: ignore
        mskf = MultilabelStratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        X_dummy = np.zeros((Y.shape[0], 1), dtype=np.int8)
        return [test_idx for _, test_idx in mskf.split(X_dummy, Y)]
    except ImportError:
        warnings.warn(
            "Package 'iterative-stratification' is not installed. "
            "Falling back to a deterministic greedy multi-label splitter.",
            RuntimeWarning,
        )
        return greedy_kfold_indices(Y, n_splits=n_splits, seed=seed)


def greedy_kfold_indices(Y: np.ndarray, n_splits: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    n, m = Y.shape

    target_pos = Y.sum(axis=0) / float(n_splits)
    fold_pos = np.zeros((n_splits, m), dtype=float)
    fold_size = np.zeros(n_splits, dtype=int)
    target_size = n / float(n_splits)

    label_sum = Y.sum(axis=1)
    order = np.lexsort((rng.random(n), -label_sum))

    assignment = np.full(n, -1, dtype=int)
    for idx in order:
        active = Y[idx].astype(bool)
        if active.any():
            deficit_score = (target_pos - fold_pos)[:, active].sum(axis=1)
        else:
            deficit_score = np.zeros(n_splits)
        size_penalty = 0.05 * np.abs((fold_size + 1) - target_size)
        score = deficit_score - size_penalty
        k = int(np.argmax(score))
        assignment[idx] = k
        fold_pos[k] += Y[idx]
        fold_size[k] += 1

    return [np.where(assignment == k)[0] for k in range(n_splits)]


def fold_objective(
    idx: np.ndarray,
    *,
    Y: np.ndarray,
    n_target: float,
    cell_counts: np.ndarray,
    cell_target: float,
    total_positive: np.ndarray,
    cell_weight: float,
) -> float:
    if len(idx) == 0:
        return float("inf")

    y = Y[idx].sum(axis=0)
    y_target = total_positive * (len(idx) / max(len(Y), 1))
    label_err = float(np.sum(np.square(y - y_target)))
    donor_err = float((len(idx) - n_target) ** 2)

    cell_sum = float(cell_counts[idx].sum())
    cell_err = float(((cell_sum - cell_target) / max(cell_target, 1.0)) ** 2)

    return label_err + donor_err + cell_weight * cell_err


def choose_holdout_fold(
    candidate_folds: list[np.ndarray],
    *,
    Y: np.ndarray,
    cell_counts: np.ndarray,
    target_frac: float,
    cell_weight: float,
) -> np.ndarray:
    n = Y.shape[0]
    n_target = n * target_frac
    cell_target = cell_counts.sum() * target_frac
    total_positive = Y.sum(axis=0)

    scores = [
        fold_objective(
            idx,
            Y=Y,
            n_target=n_target,
            cell_counts=cell_counts,
            cell_target=cell_target,
            total_positive=total_positive,
            cell_weight=cell_weight,
        )
        for idx in candidate_folds
    ]
    return candidate_folds[int(np.argmin(scores))]


def assign_splits(
    donor_df: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    primary_cols: list[str],
    weak_cols: list[str],
    cfg: PathologyAwareSplitConfigV4,
) -> pd.DataFrame:
    if not (0.0 < cfg.final_holdout_frac < 0.5):
        raise ValueError("final_holdout_frac must be in (0, 0.5).")
    if cfg.n_folds < 2:
        raise ValueError("n_folds must be >= 2.")
    if not (0 <= cfg.val_fold_index < cfg.n_folds):
        raise ValueError("val_fold_index must be within [0, n_folds).")

    donor_ids = donor_df.index.astype(str).to_numpy()
    strat_cols = primary_cols + weak_cols
    Y = labels.loc[donor_ids, strat_cols].fillna(0).astype(np.int8).to_numpy()

    cell_counts = donor_df["n_cells"].to_numpy(dtype=float) if "n_cells" in donor_df else np.ones(len(donor_ids))

    # Stage 1: locked hold-out.
    pseudo_folds = max(2, int(round(1.0 / cfg.final_holdout_frac))) #??
    candidate_holdouts = iterative_kfold_indices(Y, n_splits=pseudo_folds, seed=cfg.seed)
    holdout_idx = choose_holdout_fold(
        candidate_holdouts,
        Y=Y,
        cell_counts=cell_counts,
        target_frac=cfg.final_holdout_frac,
        cell_weight=cfg.holdout_cell_count_weight,
    )
    holdout_mask = np.zeros(len(donor_ids), dtype=bool)
    holdout_mask[holdout_idx] = True

    # Stage 2: K-fold on the remaining donor pool.
    remain_idx = np.where(~holdout_mask)[0]
    Y_remain = Y[remain_idx]
    cv_folds_local = iterative_kfold_indices(Y_remain, n_splits=cfg.n_folds, seed=cfg.seed + 1)

    split = np.empty(len(donor_ids), dtype=object)
    fold_id = np.full(len(donor_ids), -1, dtype=int)

    split[holdout_idx] = "test"

    for fold, local_idx in enumerate(cv_folds_local):
        global_idx = remain_idx[local_idx]
        fold_id[global_idx] = fold
        split[global_idx] = "val" if fold == cfg.val_fold_index else "train"

    out = donor_df.copy()
    out[cfg.split_key] = split
    out[cfg.fold_key] = fold_id
    out["pathology_group"] = [report_group(labels.loc[d]) for d in donor_ids]

    for c in labels.columns:
        out[f"label__{c}"] = labels.loc[donor_ids, c].to_numpy()

    return out


# =============================================================================
# 7. Reporting
# =============================================================================
def balance_report(
    donor_split: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    primary_cols: list[str],
    weak_cols: list[str],
    report_cols: list[str],
    cfg: PathologyAwareSplitConfigV4,
) -> pd.DataFrame:
    roles = (
        {c: "primary_stratification_axis" for c in primary_cols}
        | {c: "weak_balance_axis" for c in weak_cols}
        | {c: "report_only_axis" for c in report_cols}
    )

    rows = []
    for c, role in roles.items():
        col = f"label__{c}"
        if col not in donor_split.columns:
            continue

        for split, g in donor_split.groupby(cfg.split_key, observed=True):
            vals = g[col].astype(float)
            rows.append(
                {
                    "axis": c,
                    "role": role,
                    "stratum_kind": "split",
                    "stratum": str(split),
                    "n_donors": int(len(g)),
                    "n_positive": int(vals.sum()),
                    "positive_rate": float(vals.mean()) if len(vals) else np.nan,
                }
            )

        for fold, g in donor_split[donor_split[cfg.fold_key] >= 0].groupby(cfg.fold_key, observed=True):
            vals = g[col].astype(float)
            rows.append(
                {
                    "axis": c,
                    "role": role,
                    "stratum_kind": "fold",
                    "stratum": int(fold),
                    "n_donors": int(len(g)),
                    "n_positive": int(vals.sum()),
                    "positive_rate": float(vals.mean()) if len(vals) else np.nan,
                }
            )

    return pd.DataFrame(rows)


def group_counts_by_split(donor_split: pd.DataFrame, cfg: PathologyAwareSplitConfigV4) -> pd.DataFrame:
    tab = pd.pivot_table(
        donor_split.assign(_n=1),
        index="pathology_group",
        columns=cfg.split_key,
        values="_n",
        aggfunc="sum",
        fill_value=0,
        observed=True,
    )
    tab["__TOTAL__"] = tab.sum(axis=1)
    return tab.sort_values("__TOTAL__", ascending=False)


# =============================================================================
# 8. Obs sidecar
# =============================================================================
def build_obs_sidecar(
    obs: pd.DataFrame,
    donor_split: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    primary_cols: list[str],
    weak_cols: list[str],
    report_cols: list[str],
    cfg: PathologyAwareSplitConfigV4,
) -> pd.DataFrame:
    obs = obs.copy()
    obs.index = obs.index.astype(str)

    if cfg.donor_key not in obs.columns:
        raise KeyError(f"donor_key={cfg.donor_key!r} not found in obs.")

    donor = obs[cfg.donor_key].astype(str)

    side = pd.DataFrame(index=obs.index)
    side[cfg.donor_key] = donor
    side[cfg.split_key] = donor.map(donor_split[cfg.split_key].to_dict()).fillna("__MISSING_DONOR_SPLIT__")
    side[cfg.fold_key] = donor.map(donor_split[cfg.fold_key].to_dict()).fillna(-1).astype(np.int64)
    side["pathology_group"] = donor.map(donor_split["pathology_group"].to_dict()).fillna("__MISSING_DONOR_GROUP__")

    # Required by make_ordinal_data_zarr.py and training dataset.
    side[cfg.celltype_key] = obs[cfg.celltype_key].astype(str) if cfg.celltype_key in obs.columns else "__UNKNOWN_CELLTYPE__"

    if cfg.tech_id_model_key in obs.columns:
        side[cfg.tech_id_model_key] = obs[cfg.tech_id_model_key].astype(str).fillna(cfg.tech_missing_value)
    else:
        chosen = next((c for c in cfg.tech_source_candidates if c in obs.columns), None)
        if chosen is None:
            warnings.warn(
                f"No tech source found among {cfg.tech_source_candidates}; using {cfg.tech_missing_value!r}.",
                RuntimeWarning,
            )
            side[cfg.tech_id_model_key] = cfg.tech_missing_value
        else:
            _log(f"creating {cfg.tech_id_model_key} from obs column {chosen!r}")
            side[cfg.tech_id_model_key] = obs[chosen].astype(str).fillna(cfg.tech_missing_value)

    # Preserve donor-level pathology metadata in cell-level sidecar.
    donor_level_cols = [
        cfg.col_disease,
        cfg.col_neurotypical,
        cfg.col_adnc,
        cfg.col_braak,
        cfg.col_thal,
        cfg.col_cerad,
        cfg.col_late,
        cfg.col_lewy,
        cfg.col_microinfarct,
        cfg.col_apoe4,
        cfg.col_cognitive,
        cfg.col_age,
        cfg.col_sex,
    ]
    for c in donor_level_cols:
        if c in donor_split.columns:
            side[c] = donor.map(donor_split[c].to_dict()).fillna("__NA__")

    # Attach labels.
    for c in primary_cols + weak_cols + report_cols:
        if c in labels.columns:
            side[c] = donor.map(labels[c].to_dict()).fillna(0).astype(np.int8)

    missing = int((side[cfg.split_key] == "__MISSING_DONOR_SPLIT__").sum())
    if missing > 0:
        warnings.warn(f"{missing:,} cells have donor_id absent from donor_split table.", RuntimeWarning)

    return side


# =============================================================================
# 9. Main
# =============================================================================
def make_split(cfg: PathologyAwareSplitConfigV4) -> dict[str, Any]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _log(f"loading donor table: {cfg.donor_pathology_table_path}")
    donor_table = pd.read_csv(cfg.donor_pathology_table_path, index_col=0)
    donor_table.index = donor_table.index.astype(str)

    if "n_cells" not in donor_table.columns:
        donor_table["n_cells"] = 1

    _log(f"donors: {len(donor_table):,}")

    labels, primary_cols, weak_cols, report_cols = build_labels(donor_table, cfg)

    _log(f"primary axes: {primary_cols}")
    _log(f"weak axes: {weak_cols}")
    _log(f"report-only axes: {report_cols}")
    _log(f"positive counts primary: {labels[primary_cols].sum().to_dict()}")

    donor_split = assign_splits(
        donor_table,
        labels,
        primary_cols=primary_cols,
        weak_cols=weak_cols,
        cfg=cfg,
    )

    _log("loading obs")
    obs = load_obs(cfg)
    sidecar = build_obs_sidecar(
        obs,
        donor_split,
        labels,
        primary_cols=primary_cols,
        weak_cols=weak_cols,
        report_cols=report_cols,
        cfg=cfg,
    )

    save_df(donor_split, out_dir / "donor_split_table", cfg, index=True)
    save_df(labels, out_dir / "donor_label_matrix", cfg, index=True)
    save_obs(sidecar, out_dir / "obs_with_pathology_split", cfg)

    balance = balance_report(
        donor_split,
        labels,
        primary_cols=primary_cols,
        weak_cols=weak_cols,
        report_cols=report_cols,
        cfg=cfg,
    )
    balance.to_csv(out_dir / "pathology_balance_report.csv", index=False)

    group_counts = group_counts_by_split(donor_split, cfg)
    save_df(group_counts, out_dir / "pathology_group_counts_by_split", cfg, index=True)

    split_counts = donor_split[cfg.split_key].value_counts().to_dict()
    fold_counts = donor_split[cfg.fold_key].value_counts().sort_index().to_dict()

    mapping_used = {
        "THAL_ORDINAL": THAL_ORDINAL,
        "BRAAK_ORDINAL": BRAAK_ORDINAL,
        "CERAD_ORDINAL": CERAD_ORDINAL,
        "ADNC_ORDINAL": ADNC_ORDINAL,
        "LATE_ORDINAL": LATE_ORDINAL,
        "LEWY_POSITIVE_PATTERNS": list(LEWY_POSITIVE_PATTERNS),
        "LEWY_NEGATIVE_PATTERNS": list(LEWY_NEGATIVE_PATTERNS),
        "MICROINFARCT_POSITIVE_PATTERNS": list(MICROINFARCT_POSITIVE_PATTERNS),
        "MICROINFARCT_NEGATIVE_PATTERNS": list(MICROINFARCT_NEGATIVE_PATTERNS),
    }
    with open(out_dir / "ordinal_mapping_used.json", "w", encoding="utf-8") as f:
        json.dump(mapping_used, f, ensure_ascii=False, indent=2)

    manifest = {
        "config": asdict(cfg),
        "n_donors": int(len(donor_split)),
        "n_cells": int(len(sidecar)),
        "primary_axes": primary_cols,
        "weak_axes": weak_cols,
        "report_only_axes": report_cols,
        "positive_counts_primary": {k: int(v) for k, v in labels[primary_cols].sum().to_dict().items()},
        "positive_counts_weak": {k: int(v) for k, v in labels[weak_cols].sum().to_dict().items()},
        "positive_counts_report": {k: int(v) for k, v in labels[report_cols].sum().to_dict().items()},
        "donor_split_counts": {str(k): int(v) for k, v in split_counts.items()},
        "donor_fold_counts": {str(k): int(v) for k, v in fold_counts.items()},
        "outputs": {
            "donor_split_table_csv": str(out_dir / "donor_split_table.csv"),
            "donor_split_table_parquet": str(out_dir / "donor_split_table.parquet"),
            "donor_label_matrix_csv": str(out_dir / "donor_label_matrix.csv"),
            "obs_with_pathology_split_parquet": str(out_dir / "obs_with_pathology_split.parquet"),
            "pathology_balance_report_csv": str(out_dir / "pathology_balance_report.csv"),
            "pathology_group_counts_by_split_csv": str(out_dir / "pathology_group_counts_by_split.csv"),
            "ordinal_mapping_used_json": str(out_dir / "ordinal_mapping_used.json"),
        },
    }

    with open(out_dir / "split_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    _log("done")
    _log(f"obs sidecar: {out_dir / 'obs_with_pathology_split.parquet'}")
    return manifest


# =============================================================================
# 10. CLI
# =============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Pathology-aware donor split v4.")
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()

    cfg = load_config(args.config)
    make_split(cfg)


if __name__ == "__main__":
    main()
