#!/usr/bin/env python3
"""Build validation-only comparison tables for completed SV6 rank-8 candidates."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_final_candidate_selection_20260828"

SOURCES = [
    {
        "id": "rank2_e20",
        "label": "Frozen rank2 e20",
        "epoch": 20,
        "directory": ROOT / "analysis_outputs" / "prism_comprehensive_checkpoint_comparison_20260827" / "raw" / "rank2_e20",
        "module_name": "module.json",
        "leakage_name": "leakage.json",
        "celltype_leakage_name": "leakage_by_celltype.json",
        "history_name": None,
        "generator_hard": 97,
        "generator_expected": 97.0,
        "generator_policy": "fixed mask",
    },
    {
        "id": "sv6_e33",
        "label": "SV6 e33 protected",
        "epoch": 33,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e33_interim_posthoc_20260826",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "interim_summary.json",
        "generator_hard": 326,
        "generator_expected": 326.0633850097656,
        "generator_policy": "326 unique-coverage protected",
    },
    {
        "id": "sv6_e42",
        "label": "SV6 e42 free",
        "epoch": 42,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e42_interim_posthoc_20260827",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "history/interim_summary.json",
        "generator_hard": 144,
        "generator_expected": 152.21578979492188,
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv6_e46",
        "label": "SV6 e46 candidate",
        "epoch": 46,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e46_candidate_posthoc_20260828",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "history/interim_summary.json",
        "generator_json": "joint_generator_count_epoch_046.json",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv6_e47",
        "label": "SV6 e47 free",
        "epoch": 47,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e47_interim_posthoc_20260827",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "history/interim_summary.json",
        "generator_json": "joint_generator_count_epoch_047.json",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv6_e50",
        "label": "SV6 e50 candidate",
        "epoch": 50,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e50_candidate_posthoc_20260828",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "history/interim_summary.json",
        "generator_json": "joint_generator_count_epoch_050.json",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv6_e55",
        "label": "SV6 e55 final",
        "epoch": 55,
        "directory": ROOT / "analysis_outputs" / "prism_integrated_rank8_primary_e55_candidate_posthoc_20260828",
        "module_name": "module_disease_val_bal40.json",
        "leakage_name": "leakage_val.json",
        "celltype_leakage_name": "leakage_by_celltype_val.json",
        "history_name": "history/interim_summary.json",
        "generator_json": "joint_generator_count_epoch_055.json",
        "generator_policy": "free, singleton-only protection",
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


def load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
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
    if isinstance(record.get("val"), dict):
        return record["val"]
    return record


def finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def ci(record: dict[str, Any] | None, index: int) -> Any:
    if not record:
        return None
    values = record.get("donor_bootstrap_95ci") or record.get("donor_bootstrap_spearman_95ci")
    return values[index] if isinstance(values, list) and len(values) == 2 else None


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if finite(row.get(key)) is None else finite(row.get(key)) for key in fields})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data: dict[str, dict[str, Any]] = {}
    for source in SOURCES:
        directory = source["directory"]
        module = load_json(directory / source["module_name"])
        leakage = load_json(directory / source["leakage_name"])
        celltype_leakage = load_json(directory / source["celltype_leakage_name"])
        strict_leakage = load_json(directory / "strict_reference_leakage_val.json")
        history = load_json(directory / source["history_name"]) if source.get("history_name") else None
        if module is None or leakage is None or celltype_leakage is None:
            raise FileNotFoundError(f"incomplete candidate artifacts: {source['id']} at {directory}")
        if leakage.get("test_used") is not False or celltype_leakage.get("test_used") is not False:
            raise AssertionError(f"test split was used for {source['id']}")
        if strict_leakage is not None and strict_leakage.get("test_used") is not False:
            raise AssertionError(f"test split was used for strict probe {source['id']}")
        generator = load_json(directory / source["generator_json"]) if source.get("generator_json") else None
        data[source["id"]] = {
            "module": module,
            "leakage": leakage,
            "celltype_leakage": celltype_leakage,
            "strict_leakage": strict_leakage,
            "history": history,
            "generator": generator,
        }

    global_rows: list[dict[str, Any]] = []
    for source in SOURCES:
        item = data[source["id"]]
        module = item["module"]
        leakage = item["leakage"]
        strict_leakage = item["strict_leakage"]
        history = item["history"]
        generator = item["generator"]
        per_donor = module.get("per_donor", [])
        sex = unwrap(nested(leakage, "cell_level_nuisance", "sex"))
        tech = unwrap(nested(leakage, "cell_level_nuisance", "technology"))
        celltype = unwrap(nested(leakage, "cell_level_nuisance", "celltype"))
        region = unwrap(nested(leakage, "cell_level_nuisance", "region"))
        umi = unwrap(nested(leakage, "cell_level_continuous_nuisance", "log1p_UMI"))
        detected = unwrap(nested(leakage, "cell_level_continuous_nuisance", "log1p_detected_genes"))
        donor_z_sex = unwrap(nested(leakage, "donor_celltype_targets", "sex", "z"))
        donor_p_sex = unwrap(nested(leakage, "donor_celltype_targets", "sex", "personal_code"))
        donor_z_adnc = unwrap(nested(leakage, "donor_celltype_targets", "ADNC_report_only", "z"))
        by_celltype = item["celltype_leakage"].get("per_celltype", {})
        z_sex_values = [
            nested(record, "z", "sex", "balanced_accuracy")
            for record in by_celltype.values()
        ]
        z_sex_values = [value for value in z_sex_values if isinstance(value, (int, float)) and math.isfinite(value)]
        p_sex_values = [
            nested(record, "personal_code", "sex", "balanced_accuracy")
            for record in by_celltype.values()
        ]
        p_sex_values = [value for value in p_sex_values if isinstance(value, (int, float)) and math.isfinite(value)]
        z_sex_auc_values = [
            nested(record, "z", "sex", "roc_auc")
            for record in by_celltype.values()
        ]
        z_sex_auc_values = [value for value in z_sex_auc_values if isinstance(value, (int, float)) and math.isfinite(value)]
        z_sex_shuffle_gaps = [
            nested(record, "z", "sex", "balanced_accuracy")
            - nested(record, "z", "sex", "donor_label_shuffle_bacc")
            for record in by_celltype.values()
            if isinstance(nested(record, "z", "sex", "balanced_accuracy"), (int, float))
            and isinstance(nested(record, "z", "sex", "donor_label_shuffle_bacc"), (int, float))
        ]
        module_by_celltype = module.get("per_celltype", {})
        rank2_module_by_celltype = data["rank2_e20"]["module"].get("per_celltype", {})
        module_values = [
            record.get("AD_module_spearman")
            for record in module_by_celltype.values()
            if isinstance(record.get("AD_module_spearman"), (int, float))
        ]
        non_neuronal_values = [
            nested(module_by_celltype, celltype_name, "AD_module_spearman")
            for celltype_name in NON_NEURONAL_CELLTYPES
        ]
        non_neuronal_values = [value for value in non_neuronal_values if isinstance(value, (int, float))]
        latest = history.get("latest", {}) if history else {}
        global_rows.append(
            {
                "id": source["id"],
                "label": source["label"],
                "epoch": source["epoch"],
                "generator_policy": source["generator_policy"],
                "generator_hard": generator.get("hard_active") if generator else source.get("generator_hard"),
                "generator_expected": generator.get("expected_active") if generator else source.get("generator_expected"),
                "module_cells": nested(module, "evaluation_provenance", "selected_cells"),
                "ad_module_median": module.get("median_AD_module_spearman"),
                "pooled_spearman": module.get("pooled_spearman"),
                "blind_centered_median": module.get("blind_module_centered_median"),
                "module_to_adnc_lam100": nested(module, "module_to_ADNC", "model_lam100"),
                "disease_recon_mean": statistics.mean(row["recon_disease"] for row in per_donor),
                "full_recon_mean": statistics.mean(row["recon_full"] for row in per_donor),
                "per_module_bias_mean": module.get("per_module_bias_mean"),
                "celltypes_better_than_rank2": sum(
                    nested(module_by_celltype, celltype_name, "AD_module_spearman")
                    > nested(rank2_module_by_celltype, celltype_name, "AD_module_spearman")
                    for celltype_name in module_by_celltype
                ),
                "ad_module_celltype_min": min(module_values),
                "non_neuronal_ad_module_mean": statistics.mean(non_neuronal_values),
                "non_neuronal_ad_module_min": min(non_neuronal_values),
                "val_rec_nll": latest.get("val:loss/rec"),
                "val_branch_nll": latest.get("val:loss/prism_branch_nll"),
                "val_nonzero_acc": latest.get("val:metric/ordinal_nonzero_acc"),
                "val_within1": latest.get("val:metric/ordinal_nonzero_within1"),
                "val_balanced_recall": latest.get("val:metric/ordinal_balanced_recall"),
                "val_z_participation": latest.get("val:metric/v20_z_participation"),
                "posthoc_full_nll": nested(leakage, "sample_nll", "val", "full"),
                "posthoc_branch_nll": nested(leakage, "sample_nll", "val", "target_latent_free_branch"),
                "leak_train_cells": nested(leakage, "cell_counts", "train"),
                "leak_val_cells": nested(leakage, "cell_counts", "val"),
                "leak_test_cells": nested(leakage, "cell_counts", "test"),
                "sex_bacc": nested(sex, "balanced_accuracy"),
                "sex_ci_low": ci(sex, 0),
                "sex_ci_high": ci(sex, 1),
                "technology_bacc": nested(tech, "balanced_accuracy"),
                "technology_ci_low": ci(tech, 0),
                "technology_ci_high": ci(tech, 1),
                "celltype_bacc": nested(celltype, "balanced_accuracy"),
                "region_bacc": nested(region, "balanced_accuracy"),
                "region_ci_low": ci(region, 0),
                "region_ci_high": ci(region, 1),
                "umi_spearman": nested(umi, "spearman"),
                "detected_genes_spearman": nested(detected, "spearman"),
                "donor_celltype_z_sex_bacc": nested(donor_z_sex, "balanced_accuracy"),
                "donor_celltype_personal_sex_bacc": nested(donor_p_sex, "balanced_accuracy"),
                "donor_celltype_z_adnc_spearman": nested(donor_z_adnc, "spearman"),
                "celltype_z_sex_median": statistics.median(z_sex_values),
                "celltype_z_sex_auc_median": statistics.median(z_sex_auc_values),
                "celltype_z_sex_max": max(z_sex_values),
                "celltype_z_sex_ge_075_count": sum(value >= 0.75 for value in z_sex_values),
                "celltype_z_sex_shuffle_gap_median": statistics.median(z_sex_shuffle_gaps),
                "celltype_z_sex_above_shuffle_count": sum(value > 0 for value in z_sex_shuffle_gaps),
                "celltype_z_sex_shuffle_gap_ge_010_count": sum(value >= 0.1 for value in z_sex_shuffle_gaps),
                "celltype_personal_sex_median": statistics.median(p_sex_values),
                "z_donor_retrieval_top1": nested(leakage, "donor_retrieval", "val", "z_cross_celltype", "top1_accuracy"),
                "personal_donor_retrieval_top1": nested(leakage, "donor_retrieval", "val", "personal_code_cross_celltype", "top1_accuracy"),
                "strict_val_donors": nested(strict_leakage, "donor_counts", "val"),
                "strict_celltype_bacc": nested(strict_leakage, "probes", "celltype", "balanced_accuracy"),
                "strict_region_bacc": nested(strict_leakage, "probes", "region", "balanced_accuracy"),
            }
        )
    write_tsv(OUT / "global_candidate_comparison.tsv", global_rows)

    celltypes = list(data["rank2_e20"]["module"]["per_celltype"])
    module_rows: list[dict[str, Any]] = []
    leakage_rows: list[dict[str, Any]] = []
    for celltype_name in celltypes:
        module_row: dict[str, Any] = {"celltype": celltype_name}
        leakage_row: dict[str, Any] = {"celltype": celltype_name}
        for source in SOURCES:
            item = data[source["id"]]
            module_record = nested(item["module"], "per_celltype", celltype_name) or {}
            leakage_record = nested(item["celltype_leakage"], "per_celltype", celltype_name) or {}
            module_row[f"{source['id']}_ad_module"] = module_record.get("AD_module_spearman")
            module_row[f"{source['id']}_blind_centered"] = module_record.get("blind_centered")
            leakage_row[f"{source['id']}_z_sex_bacc"] = nested(leakage_record, "z", "sex", "balanced_accuracy")
            leakage_row[f"{source['id']}_personal_sex_bacc"] = nested(leakage_record, "personal_code", "sex", "balanced_accuracy")
        module_rows.append(module_row)
        leakage_rows.append(leakage_row)
    write_tsv(OUT / "celltype_module_comparison.tsv", module_rows)
    write_tsv(OUT / "celltype_sex_leakage_comparison.tsv", leakage_rows)

    manifest = {
        "schema_version": "prism.sv6_rank8_final_candidate_comparison.v1",
        "official_test_used": False,
        "selection_scope": "validation-only completed-run checkpoint selection",
        "module_protocol": {
            "split": "validation",
            "donors": 16,
            "celltypes": 24,
            "modules": 404,
            "per_donor_celltype": 40,
            "seed": 20260814,
        },
        "leakage_protocol": {
            "fit": "train",
            "evaluate": "validation",
            "recent_train_cells": 44960,
            "recent_validation_cells": 11317,
            "test_cells": 0,
        },
        "fixed_architecture": {"personal_rank": 2, "pathology_rank": 8},
        "sources": [
            {key: str(value) if isinstance(value, Path) else value for key, value in source.items()}
            for source in SOURCES
        ],
        "outputs": [
            "global_candidate_comparison.tsv",
            "celltype_module_comparison.tsv",
            "celltype_sex_leakage_comparison.tsv",
            "comparison_manifest.json",
        ],
    }
    with (OUT / "comparison_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
