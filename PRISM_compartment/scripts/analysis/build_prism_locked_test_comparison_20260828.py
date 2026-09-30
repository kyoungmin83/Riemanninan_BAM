#!/usr/bin/env python3
"""Build a post-selection comparison for the explicitly consumed PRISM test set."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
TEST_ROOT = ROOT / "analysis_outputs" / "prism_locked_test_audit_20260828"
VALIDATION_ROOT = (
    ROOT
    / "analysis_outputs"
    / "prism_integrated_rank8_primary_final_candidate_selection_20260828"
)
LOCK_PATH = VALIDATION_ROOT / "PRE_TEST_SELECTION_LOCK.json"
LOCK_SHA256 = "f88f29426e625ceeb2fa3467d000db2e5d6a44a478f71afc21285019354d96ff"

SOURCES = [
    {
        "id": "rank2_e20",
        "validation_id": "rank2_e20",
        "label": "Frozen rank2 e20",
        "epoch": 20,
        "generators": 97,
        "role": "frozen baseline",
    },
    {
        "id": "sv6_e33_protected",
        "validation_id": "sv6_e33",
        "label": "SV6 e33 protected",
        "epoch": 33,
        "generators": 326,
        "role": "protected-generator control",
    },
    {
        "id": "sv6_e42",
        "validation_id": "sv6_e42",
        "label": "SV6 e42 primary",
        "epoch": 42,
        "generators": 144,
        "role": "locked validation primary",
    },
    {
        "id": "sv6_e47",
        "validation_id": "sv6_e47",
        "label": "SV6 e47 aggregate",
        "epoch": 47,
        "generators": 130,
        "role": "aggregate AD-module secondary",
    },
    {
        "id": "sv6_e50",
        "validation_id": "sv6_e50",
        "label": "SV6 e50 tail",
        "epoch": 50,
        "generators": 123,
        "role": "non-neuronal tail secondary",
    },
    {
        "id": "sv6_e55",
        "validation_id": "sv6_e55",
        "label": "SV6 e55 compression",
        "epoch": 55,
        "generators": 114,
        "role": "maximum-compression ablation",
    },
]

NON_NEURONAL_CELLTYPES = {
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
}


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def nested(data: dict[str, Any] | None, *keys: str) -> Any:
    value: Any = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def unwrap(record: Any) -> dict[str, Any] | None:
    if not isinstance(record, dict):
        return None
    split_records = [record[key] for key in ("test", "val") if isinstance(record.get(key), dict)]
    if len(split_records) == 1:
        return split_records[0]
    return record


def ci(record: dict[str, Any] | None, index: int) -> Any:
    if not record:
        return None
    values = record.get("donor_bootstrap_95ci") or record.get("donor_bootstrap_spearman_95ci")
    return values[index] if isinstance(values, list) and len(values) == 2 else None


def finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def numeric(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def difference(test_value: Any, validation_value: Any) -> float | None:
    if not isinstance(test_value, (int, float)) or not isinstance(validation_value, (int, float)):
        return None
    return float(test_value) - float(validation_value)


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def max_or_none(values: list[float]) -> float | None:
    return max(values) if values else None


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
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
    lock = load_json(LOCK_PATH)
    if lock.get("test_audit_may_change_selection") is not False:
        raise AssertionError("pre-test lock permits post-test reselection")

    with (VALIDATION_ROOT / "global_candidate_comparison.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        validation_rows = {row["id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    data: dict[str, dict[str, Any]] = {}
    for source in SOURCES:
        directory = TEST_ROOT / source["id"]
        module = load_json(directory / "module_disease_test_bal40.json")
        leakage = load_json(directory / "leakage_test.json")
        by_celltype = load_json(directory / "leakage_by_celltype_test.json")
        strict = load_json(directory / "strict_reference_leakage_test.json")
        if module["evaluation_provenance"]["split"] != "test":
            raise AssertionError(f"module split is not test: {source['id']}")
        if len(module.get("per_donor", [])) != 9:
            raise AssertionError(f"test donor count is not 9: {source['id']}")
        for name, result in (
            ("leakage", leakage),
            ("celltype leakage", by_celltype),
            ("strict leakage", strict),
        ):
            if result.get("test_used") is not True:
                raise AssertionError(f"{name} did not declare test use: {source['id']}")
        if strict.get("pretest_lock_sha256") != LOCK_SHA256:
            raise AssertionError(f"strict result lock mismatch: {source['id']}")
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
        validation = validation_rows[source["validation_id"]]
        module_by_celltype = module["per_celltype"]
        module_values = [record["AD_module_spearman"] for record in module_by_celltype.values()]
        non_neuronal_values = [
            module_by_celltype[celltype]["AD_module_spearman"]
            for celltype in NON_NEURONAL_CELLTYPES
        ]
        per_donor = module["per_donor"]

        sex = unwrap(nested(leakage, "cell_level_nuisance", "sex"))
        tech = unwrap(nested(leakage, "cell_level_nuisance", "technology"))
        celltype = unwrap(nested(leakage, "cell_level_nuisance", "celltype"))
        region = unwrap(nested(leakage, "cell_level_nuisance", "region"))
        umi = unwrap(nested(leakage, "cell_level_continuous_nuisance", "log1p_UMI"))
        detected = unwrap(
            nested(leakage, "cell_level_continuous_nuisance", "log1p_detected_genes")
        )
        donor_z_sex = unwrap(nested(leakage, "donor_celltype_targets", "sex", "z"))
        donor_p_sex = unwrap(
            nested(leakage, "donor_celltype_targets", "sex", "personal_code")
        )
        donor_z_adnc = unwrap(
            nested(leakage, "donor_celltype_targets", "ADNC_report_only", "z")
        )

        z_sex = [nested(record, "z", "sex") for record in by_celltype.values()]
        z_sex = [record for record in z_sex if isinstance(record, dict) and record.get("status") == "ok"]
        p_sex_values = [
            nested(record, "personal_code", "sex", "balanced_accuracy")
            for record in by_celltype.values()
        ]
        p_sex_values = [
            value
            for value in p_sex_values
            if isinstance(value, (int, float)) and math.isfinite(value)
        ]
        z_sex_values = [record["balanced_accuracy"] for record in z_sex]
        z_sex_auc_values = [record["roc_auc"] for record in z_sex]
        z_sex_gaps = [
            record["balanced_accuracy"] - record["donor_label_shuffle_bacc"]
            for record in z_sex
        ]

        test_ad = module["median_AD_module_spearman"]
        test_pooled = module["pooled_spearman"]
        test_blind = module["blind_module_centered_median"]
        test_disease = statistics.mean(record["recon_disease"] for record in per_donor)
        test_full_nll = nested(leakage, "sample_nll", "test", "full")
        global_rows.append(
            {
                "id": source["id"],
                "label": source["label"],
                "locked_role": source["role"],
                "epoch": source["epoch"],
                "generators": source["generators"],
                "test_module_cells": nested(module, "evaluation_provenance", "selected_cells"),
                "test_donors": len(per_donor),
                "validation_ad_module": numeric(validation.get("ad_module_median")),
                "test_ad_module": test_ad,
                "test_minus_validation_ad_module": difference(
                    test_ad, numeric(validation.get("ad_module_median"))
                ),
                "validation_pooled": numeric(validation.get("pooled_spearman")),
                "test_pooled": test_pooled,
                "test_minus_validation_pooled": difference(
                    test_pooled, numeric(validation.get("pooled_spearman"))
                ),
                "validation_blind_centered": numeric(
                    validation.get("blind_centered_median")
                ),
                "test_blind_centered": test_blind,
                "test_minus_validation_blind_centered": difference(
                    test_blind, numeric(validation.get("blind_centered_median"))
                ),
                "validation_disease_recon": numeric(validation.get("disease_recon_mean")),
                "test_disease_recon": test_disease,
                "test_minus_validation_disease_recon": difference(
                    test_disease, numeric(validation.get("disease_recon_mean"))
                ),
                "test_celltypes_better_than_rank2": sum(
                    record["AD_module_spearman"]
                    > rank2_modules[celltype_name]["AD_module_spearman"]
                    for celltype_name, record in module_by_celltype.items()
                ),
                "test_ad_module_celltype_min": min(module_values),
                "test_non_neuronal_mean": statistics.mean(non_neuronal_values),
                "test_non_neuronal_min": min(non_neuronal_values),
                "validation_posthoc_full_nll": numeric(validation.get("posthoc_full_nll")),
                "test_posthoc_full_nll": test_full_nll,
                "test_minus_validation_full_nll": difference(
                    test_full_nll, numeric(validation.get("posthoc_full_nll"))
                ),
                "test_posthoc_branch_nll": nested(
                    leakage, "sample_nll", "test", "target_latent_free_branch"
                ),
                "leak_train_cells": nested(leakage, "cell_counts", "train"),
                "leak_test_cells": nested(leakage, "cell_counts", "test"),
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
                "donor_celltype_personal_sex_bacc": nested(
                    donor_p_sex, "balanced_accuracy"
                ),
                "donor_celltype_z_adnc_spearman": nested(donor_z_adnc, "spearman"),
                "celltype_z_sex_evaluable_count": len(z_sex_values),
                "celltype_z_sex_median": median_or_none(z_sex_values),
                "celltype_z_sex_auc_median": median_or_none(z_sex_auc_values),
                "celltype_z_sex_max": max_or_none(z_sex_values),
                "celltype_z_sex_shuffle_gap_median": median_or_none(z_sex_gaps),
                "celltype_z_sex_above_shuffle_count": sum(gap > 0 for gap in z_sex_gaps),
                "celltype_z_sex_shuffle_gap_ge_010_count": sum(
                    gap >= 0.1 for gap in z_sex_gaps
                ),
                "celltype_personal_sex_median": median_or_none(p_sex_values),
                "z_donor_retrieval_top1": nested(
                    leakage, "donor_retrieval", "test", "z_cross_celltype", "top1_accuracy"
                ),
                "personal_donor_retrieval_top1": nested(
                    leakage,
                    "donor_retrieval",
                    "test",
                    "personal_code_cross_celltype",
                    "top1_accuracy",
                ),
                "strict_test_donors": nested(strict, "donor_counts", "test"),
                "strict_celltype_bacc": nested(
                    strict, "probes", "celltype", "balanced_accuracy"
                ),
                "strict_sex_bacc": nested(strict, "probes", "sex", "balanced_accuracy"),
                "strict_technology_bacc": nested(
                    strict, "probes", "technology", "balanced_accuracy"
                ),
                "strict_region_bacc": nested(strict, "probes", "region", "balanced_accuracy"),
                "test_used": True,
                "selection_locked_before_test": True,
            }
        )

    write_tsv(TEST_ROOT / "global_locked_test_comparison.tsv", global_rows)
    with (TEST_ROOT / "global_locked_test_comparison.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(global_rows, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")

    celltypes = list(data["rank2_e20"]["module"]["per_celltype"])
    module_rows: list[dict[str, Any]] = []
    leakage_rows: list[dict[str, Any]] = []
    for celltype_name in celltypes:
        module_row: dict[str, Any] = {"celltype": celltype_name}
        leakage_row: dict[str, Any] = {"celltype": celltype_name}
        for source in SOURCES:
            item = data[source["id"]]
            module_record = item["module"]["per_celltype"][celltype_name]
            leakage_record = item["by_celltype"]["per_celltype"].get(celltype_name, {})
            module_row[f"{source['id']}_ad_module"] = module_record.get(
                "AD_module_spearman"
            )
            module_row[f"{source['id']}_blind_centered"] = module_record.get(
                "blind_centered"
            )
            leakage_row[f"{source['id']}_z_sex_bacc"] = nested(
                leakage_record, "z", "sex", "balanced_accuracy"
            )
            leakage_row[f"{source['id']}_personal_sex_bacc"] = nested(
                leakage_record, "personal_code", "sex", "balanced_accuracy"
            )
        module_rows.append(module_row)
        leakage_rows.append(leakage_row)
    write_tsv(TEST_ROOT / "celltype_locked_test_module_comparison.tsv", module_rows)
    write_tsv(TEST_ROOT / "celltype_locked_test_sex_leakage_comparison.tsv", leakage_rows)

    manifest = {
        "schema_version": "kmlee_bam.prism_locked_test_audit.v1",
        "official_test_status": "consumed",
        "official_test_first_open_date": "2026-08-28",
        "selection_locked_before_test": True,
        "test_results_may_change_locked_selection": False,
        "pretest_lock": str(LOCK_PATH),
        "pretest_lock_sha256": LOCK_SHA256,
        "primary_id_locked_before_test": "sv6_e42",
        "excluded_before_test": ["sv6_e46"],
        "module_protocol": {
            "split": "test",
            "donors": 9,
            "celltypes": 24,
            "modules": 404,
            "selected_cells": 8478,
            "per_donor_celltype": 40,
            "seed": 20260814,
        },
        "leakage_protocol": {
            "fit": "train",
            "evaluate": "test",
            "per_donor_celltype": 30,
            "seed": 20260814,
        },
        "sources": SOURCES,
        "outputs": [
            "global_locked_test_comparison.tsv",
            "global_locked_test_comparison.json",
            "celltype_locked_test_module_comparison.tsv",
            "celltype_locked_test_sex_leakage_comparison.tsv",
            "locked_test_manifest.json",
        ],
    }
    with (TEST_ROOT / "locked_test_manifest.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
