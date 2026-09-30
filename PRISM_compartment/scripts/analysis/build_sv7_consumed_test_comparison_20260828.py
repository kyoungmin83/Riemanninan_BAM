#!/usr/bin/env python3
"""Build descriptive SV7 test-reuse tables without changing validation roles."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
OUT_ROOT = ROOT / "analysis_outputs" / "prism_sv7_pathrank_consumed_test_posthoc_20260828"
BASE_ROOT = ROOT / "analysis_outputs" / "prism_locked_test_audit_20260828"
VAL_ROOT = ROOT / "analysis_outputs" / "prism_integrated_pathranklearn_sv7_candidate_comparison_20260828"
SCOPE_PATH = OUT_ROOT / "SV7_CONSUMED_TEST_SCOPE_LOCK.json"
SCOPE_SHA256 = "5105a96025c69f90d7b19128e3a1e52ffc80608d5d546fd574c81f6a87989790"

NON_NEURONAL_CELLTYPES = {
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
}

SOURCES = [
    {
        "id": "rank2_e20",
        "validation_id": "rank2_e20",
        "label": "Frozen rank2 e20",
        "role": "frozen baseline",
        "directory": BASE_ROOT / "rank2_e20",
        "sv7": False,
    },
    {
        "id": "sv6_e33_protected",
        "validation_id": "sv6_fixed_e33",
        "label": "SV6 fixed-rank8 e33 protected",
        "role": "prior integrated protected control",
        "directory": BASE_ROOT / "sv6_e33_protected",
        "sv7": False,
    },
    {
        "id": "sv6_e42",
        "validation_id": "sv6_fixed_e42",
        "label": "SV6 fixed-rank8 e42 free",
        "role": "prior validation-selected integrated primary",
        "directory": BASE_ROOT / "sv6_e42",
        "sv7": False,
    },
    {
        "id": "sv7_pathrank_e33",
        "validation_id": "sv7_pathrank_e33",
        "label": "SV7 learned-path-rank e33",
        "role": "balanced validation primary",
        "directory": OUT_ROOT / "sv7_pathrank_e33",
        "sv7": True,
    },
    {
        "id": "sv7_pathrank_e36",
        "validation_id": "sv7_pathrank_e36",
        "label": "SV7 learned-path-rank e36",
        "role": "aggregate-biology validation secondary",
        "directory": OUT_ROOT / "sv7_pathrank_e36",
        "sv7": True,
    },
    {
        "id": "sv7_pathrank_e39",
        "validation_id": "sv7_pathrank_e39",
        "label": "SV7 learned-path-rank e39",
        "role": "disease/tail secondary; validation-primary rejected for sex leakage",
        "directory": OUT_ROOT / "sv7_pathrank_e39",
        "sv7": True,
    },
]


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def nested(data: Any, *keys: str) -> Any:
    value = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def unwrap_split(record: Any, split: str = "test") -> dict[str, Any] | None:
    if not isinstance(record, dict):
        return None
    if isinstance(record.get(split), dict):
        return record[split]
    return record


def finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def numeric(value: Any) -> float | None:
    if value in (None, ""):
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def mean_or_none(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def max_or_none(values: list[float]) -> float | None:
    return max(values) if values else None


def difference(left: Any, right: Any) -> float | None:
    if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
        return None
    return float(left) - float(right)


def ci(record: dict[str, Any] | None, index: int) -> Any:
    if not record:
        return None
    values = record.get("donor_bootstrap_95ci") or record.get(
        "donor_bootstrap_spearman_95ci"
    )
    return values[index] if isinstance(values, list) and len(values) == 2 else None


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: "" if finite(row.get(key)) is None else finite(row.get(key))
                    for key in fields
                }
            )


def main() -> None:
    scope = load_json(SCOPE_PATH)
    assert scope["official_test_status_at_scope_lock"] == "consumed_by_prior_sv6_locked_audit"
    assert scope["test_reuse_may_change_validation_roles"] is False
    assert scope["test_reuse_may_authorize_baseline_replacement"] is False
    locked_roles = {item["id"]: item for item in scope["locked_roles"]}
    expected_locked_roles = {
        "sv7_pathrank_e33": "balanced_validation_primary",
        "sv7_pathrank_e36": "aggregate_biology_secondary",
        "sv7_pathrank_e39": "disease_and_non_neuronal_tail_secondary_rejected_as_primary_for_validation_sex_leakage",
    }

    with (VAL_ROOT / "global_candidate_comparison.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        validation = {row["id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    data: dict[str, dict[str, Any]] = {}
    for source in SOURCES:
        directory = source["directory"]
        module = load_json(directory / "module_disease_test_bal40.json")
        leakage = load_json(directory / "leakage_test.json")
        by_celltype = load_json(directory / "leakage_by_celltype_test.json")
        strict = load_json(directory / "strict_reference_leakage_test.json")
        assert module["evaluation_provenance"]["split"] == "test"
        assert module["evaluation_provenance"]["selected_cells"] == 8478
        assert len(module["per_donor"]) == 9
        assert len(module["per_celltype"]) == 24
        assert leakage["test_used"] is True
        assert leakage["cell_counts"]["test"] == 6411
        assert leakage["cell_counts"].get("val", 0) == 0
        assert by_celltype["test_used"] is True
        assert len(by_celltype["per_celltype"]) == 24
        assert strict["test_used"] is True
        if source["sv7"]:
            assert (directory / "CONSUMED_TEST_POSTHOC_COMPLETE").is_file()
            assert by_celltype["official_test_status"] == "consumed_reuse"
            assert by_celltype["scope_lock_sha256"] == SCOPE_SHA256
            assert strict["scope_lock_sha256"] == SCOPE_SHA256
            assert locked_roles[source["id"]]["role"] == expected_locked_roles[source["id"]]
        data[source["id"]] = {
            "module": module,
            "leakage": leakage,
            "by_celltype": by_celltype,
            "strict": strict,
        }

    rank2_modules = data["rank2_e20"]["module"]["per_celltype"]
    global_rows: list[dict[str, Any]] = []
    for source in SOURCES:
        item = data[source["id"]]
        module = item["module"]
        leakage = item["leakage"]
        by_celltype = item["by_celltype"]["per_celltype"]
        strict = item["strict"]
        val = validation[source["validation_id"]]
        module_by_celltype = module["per_celltype"]
        module_values = [row["AD_module_spearman"] for row in module_by_celltype.values()]
        non_neuronal = [
            module_by_celltype[name]["AD_module_spearman"]
            for name in sorted(NON_NEURONAL_CELLTYPES)
        ]
        per_donor = module["per_donor"]

        sex = unwrap_split(nested(leakage, "cell_level_nuisance", "sex"))
        tech = unwrap_split(nested(leakage, "cell_level_nuisance", "technology"))
        celltype = unwrap_split(nested(leakage, "cell_level_nuisance", "celltype"))
        region = unwrap_split(nested(leakage, "cell_level_nuisance", "region"))
        umi = unwrap_split(nested(leakage, "cell_level_continuous_nuisance", "log1p_UMI"))
        detected = unwrap_split(
            nested(leakage, "cell_level_continuous_nuisance", "log1p_detected_genes")
        )
        donor_z_sex = unwrap_split(nested(leakage, "donor_celltype_targets", "sex", "z"))
        donor_p_sex = unwrap_split(
            nested(leakage, "donor_celltype_targets", "sex", "personal_code")
        )
        donor_z_adnc = unwrap_split(
            nested(leakage, "donor_celltype_targets", "ADNC_report_only", "z")
        )

        z_sex_records = [nested(row, "z", "sex") for row in by_celltype.values()]
        z_sex_records = [
            row for row in z_sex_records if isinstance(row, dict) and row.get("status") == "ok"
        ]
        p_sex_records = [nested(row, "personal_code", "sex") for row in by_celltype.values()]
        p_sex_records = [
            row for row in p_sex_records if isinstance(row, dict) and row.get("status") == "ok"
        ]
        z_sex_values = [float(row["balanced_accuracy"]) for row in z_sex_records]
        p_sex_values = [float(row["balanced_accuracy"]) for row in p_sex_records]
        z_sex_gaps = [
            float(row["balanced_accuracy"] - row["donor_label_shuffle_bacc"])
            for row in z_sex_records
            if isinstance(row.get("donor_label_shuffle_bacc"), (int, float))
        ]
        p_sex_gaps = [
            float(row["balanced_accuracy"] - row["donor_label_shuffle_bacc"])
            for row in p_sex_records
            if isinstance(row.get("donor_label_shuffle_bacc"), (int, float))
        ]

        test_ad = module["median_AD_module_spearman"]
        test_pooled = module["pooled_spearman"]
        test_blind = module["blind_module_centered_median"]
        test_disease = statistics.mean(row["recon_disease"] for row in per_donor)
        full_nll = nested(leakage, "sample_nll", "test", "full")
        global_rows.append(
            {
                "id": source["id"],
                "label": source["label"],
                "locked_validation_role": source["role"],
                "epoch": numeric(val.get("epoch")),
                "personal_rank": numeric(val.get("personal_rank")),
                "pathology_hard": numeric(val.get("pathology_hard")),
                "pathology_expected": numeric(val.get("pathology_expected")),
                "pathology_finalized": val.get("pathology_finalized"),
                "generator_hard": numeric(val.get("generator_hard")),
                "generator_expected": numeric(val.get("generator_expected")),
                "validation_ad_module": numeric(val.get("ad_module")),
                "test_ad_module": test_ad,
                "test_minus_validation_ad_module": difference(test_ad, numeric(val.get("ad_module"))),
                "validation_pooled": numeric(val.get("pooled")),
                "test_pooled": test_pooled,
                "test_minus_validation_pooled": difference(test_pooled, numeric(val.get("pooled"))),
                "validation_blind_centered": numeric(val.get("blind_centered")),
                "test_blind_centered": test_blind,
                "test_minus_validation_blind_centered": difference(test_blind, numeric(val.get("blind_centered"))),
                "validation_disease_recon": numeric(val.get("disease_recon")),
                "test_disease_recon": test_disease,
                "test_minus_validation_disease_recon": difference(test_disease, numeric(val.get("disease_recon"))),
                "test_celltypes_better_than_rank2": sum(
                    row["AD_module_spearman"] > rank2_modules[name]["AD_module_spearman"]
                    for name, row in module_by_celltype.items()
                ),
                "test_ad_module_celltype_mean": statistics.mean(module_values),
                "test_ad_module_celltype_min": min(module_values),
                "test_non_neuronal_mean": statistics.mean(non_neuronal),
                "test_non_neuronal_min": min(non_neuronal),
                "validation_full_nll": numeric(val.get("val_rec_nll")) or numeric(val.get("posthoc_full_nll")),
                "test_full_nll": full_nll,
                "test_branch_nll": nested(leakage, "sample_nll", "test", "target_latent_free_branch"),
                "sex_bacc": nested(sex, "balanced_accuracy"),
                "sex_ci_low": ci(sex, 0),
                "sex_ci_high": ci(sex, 1),
                "technology_bacc": nested(tech, "balanced_accuracy"),
                "technology_ci_low": ci(tech, 0),
                "technology_ci_high": ci(tech, 1),
                "celltype_bacc": nested(celltype, "balanced_accuracy"),
                "region_bacc": nested(region, "balanced_accuracy"),
                "umi_spearman": nested(umi, "spearman"),
                "detected_genes_spearman": nested(detected, "spearman"),
                "donor_celltype_z_sex_bacc": nested(donor_z_sex, "balanced_accuracy"),
                "donor_celltype_personal_sex_bacc": nested(donor_p_sex, "balanced_accuracy"),
                "donor_celltype_z_adnc_spearman": nested(donor_z_adnc, "spearman"),
                "z_donor_retrieval_top1": nested(leakage, "donor_retrieval", "test", "z_cross_celltype", "top1_accuracy"),
                "personal_donor_retrieval_top1": nested(leakage, "donor_retrieval", "test", "personal_code_cross_celltype", "top1_accuracy"),
                "celltype_z_sex_evaluable_count": len(z_sex_values),
                "celltype_z_sex_median": median_or_none(z_sex_values),
                "celltype_z_sex_max": max_or_none(z_sex_values),
                "celltype_z_sex_shuffle_gap_median": median_or_none(z_sex_gaps),
                "celltype_personal_sex_evaluable_count": len(p_sex_values),
                "celltype_personal_sex_median": median_or_none(p_sex_values),
                "celltype_personal_sex_max": max_or_none(p_sex_values),
                "celltype_personal_sex_shuffle_gap_median": median_or_none(p_sex_gaps),
                "strict_test_donors": nested(strict, "donor_counts", "test"),
                "test_donors": len(per_donor),
                "module_test_cells": nested(module, "evaluation_provenance", "selected_cells"),
                "leak_train_cells": nested(leakage, "cell_counts", "train"),
                "leak_test_cells": nested(leakage, "cell_counts", "test"),
                "official_test_status": "consumed_reuse" if source["sv7"] else "consumed_initial_audit",
                "independent_sv7_holdout": False,
                "validation_roles_may_change": False,
            }
        )

    write_tsv(OUT_ROOT / "global_consumed_test_comparison.tsv", global_rows)
    with (OUT_ROOT / "global_consumed_test_comparison.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(global_rows, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    celltypes = list(rank2_modules)
    module_rows: list[dict[str, Any]] = []
    sex_rows: list[dict[str, Any]] = []
    for celltype_name in celltypes:
        module_row: dict[str, Any] = {"celltype": celltype_name}
        sex_row: dict[str, Any] = {"celltype": celltype_name}
        for source in SOURCES:
            record = data[source["id"]]["module"]["per_celltype"][celltype_name]
            leak_record = data[source["id"]]["by_celltype"]["per_celltype"][celltype_name]
            module_row[f"{source['id']}_ad_module"] = record.get("AD_module_spearman")
            module_row[f"{source['id']}_blind_centered"] = record.get("blind_centered")
            module_row[f"{source['id']}_delta_ad_vs_rank2"] = difference(
                record.get("AD_module_spearman"), rank2_modules[celltype_name]["AD_module_spearman"]
            )
            sex_row[f"{source['id']}_n_test_donors"] = leak_record.get("n_test_donors")
            sex_row[f"{source['id']}_z_sex_status"] = nested(leak_record, "z", "sex", "status")
            sex_row[f"{source['id']}_z_sex_bacc"] = nested(leak_record, "z", "sex", "balanced_accuracy")
            sex_row[f"{source['id']}_personal_sex_status"] = nested(
                leak_record, "personal_code", "sex", "status"
            )
            sex_row[f"{source['id']}_personal_sex_bacc"] = nested(
                leak_record, "personal_code", "sex", "balanced_accuracy"
            )
        module_rows.append(module_row)
        sex_rows.append(sex_row)
    write_tsv(OUT_ROOT / "celltype_consumed_test_module_comparison.tsv", module_rows)
    write_tsv(OUT_ROOT / "celltype_consumed_test_sex_leakage_comparison.tsv", sex_rows)

    sv7_rows = [row for row in global_rows if row["id"].startswith("sv7_")]
    metric_directions = {
        "test_ad_module": "max",
        "test_pooled": "max",
        "test_blind_centered": "max",
        "test_disease_recon": "max",
        "test_full_nll": "min",
        "test_non_neuronal_mean": "max",
        "test_non_neuronal_min": "max",
    }
    descriptive_leaders = {}
    for metric, direction in metric_directions.items():
        usable = [row for row in sv7_rows if isinstance(row.get(metric), (int, float))]
        key = (lambda row: row[metric])
        winner = min(usable, key=key) if direction == "min" else max(usable, key=key)
        descriptive_leaders[metric] = {
            "id": winner["id"],
            "value": winner[metric],
            "direction": direction,
        }

    conclusion = {
        "schema_version": "kmlee_bam.sv7_consumed_test_posthoc_conclusion.v1",
        "official_test_status": "consumed_reuse",
        "independent_holdout": False,
        "validation_primary_before_test_reuse": scope["validation_primary_id"],
        "validation_primary_after_test_reuse": scope["validation_primary_id"],
        "baseline_replacement_authorized": False,
        "test_based_epoch_reselection_authorized": False,
        "descriptive_sv7_metric_leaders": descriptive_leaders,
        "confirmatory_next_step": "new untouched external cohort, new holdout, or nested cross-validation",
    }
    with (OUT_ROOT / "POSTHOC_CONCLUSION.json").open("w", encoding="utf-8") as handle:
        json.dump(conclusion, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    manifest = {
        "schema_version": "kmlee_bam.sv7_consumed_test_posthoc_manifest.v1",
        "scope_lock": str(SCOPE_PATH),
        "scope_lock_sha256": SCOPE_SHA256,
        "official_test_status": "consumed_reuse",
        "independent_holdout": False,
        "validation_roles_may_change": False,
        "module_protocol": {
            "split": "test",
            "donors": 9,
            "cells": 8478,
            "celltypes": 24,
            "modules": 404,
            "per_donor_celltype": 40,
            "seed": 20260814,
        },
        "leakage_protocol": {
            "fit": "train",
            "evaluate": "test",
            "train_cells": 44960,
            "test_cells": 6411,
            "per_donor_celltype": 30,
            "seed": 20260814,
        },
        "sources": [
            {
                "id": source["id"],
                "label": source["label"],
                "role": source["role"],
                "directory": str(source["directory"]),
            }
            for source in SOURCES
        ],
        "outputs": [
            "global_consumed_test_comparison.tsv",
            "global_consumed_test_comparison.json",
            "celltype_consumed_test_module_comparison.tsv",
            "celltype_consumed_test_sex_leakage_comparison.tsv",
            "POSTHOC_CONCLUSION.json",
            "comparison_manifest.json",
        ],
    }
    with (OUT_ROOT / "comparison_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
