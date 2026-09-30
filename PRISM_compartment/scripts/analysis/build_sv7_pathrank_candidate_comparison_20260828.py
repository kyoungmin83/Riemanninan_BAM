#!/usr/bin/env python3
"""Build the validation-only comparison for the active SV7 path-rank run."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE_RAW = ROOT / "analysis_outputs" / "prism_comprehensive_checkpoint_comparison_20260827" / "raw"
OUT = ROOT / "analysis_outputs" / "prism_integrated_pathranklearn_sv7_candidate_comparison_20260828"

NON_NEURONAL = [
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
]

CHECKPOINTS: list[dict[str, Any]] = [
    {
        "id": "rank2_e20",
        "label": "Frozen rank2 e20",
        "epoch": 20,
        "source": BASE_RAW / "rank2_e20",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "pathology_finalized": True,
        "generator_hard": 97,
        "generator_expected": 97.0,
        "generator_policy": "fixed mask",
    },
    {
        "id": "sv7_pathrank_e14",
        "label": "SV7 path-rank e14",
        "epoch": 14,
        "source": BASE_RAW / "sv7_pathrank_e14",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "pathology_hard": 96,
        "pathology_expected": 73.95681762695312,
        "pathology_finalized": False,
        "generator_hard": 414,
        "generator_expected": 408.3619384765625,
        "generator_policy": "326 unique-coverage protected",
    },
    {
        "id": "sv7_pathrank_e26",
        "label": "SV7 path-rank e26",
        "epoch": 26,
        "source": BASE_RAW / "sv7_pathrank_e26",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "pathology_hard": 72,
        "pathology_expected": 50.01966094970703,
        "pathology_finalized": False,
        "generator_hard": 26,
        "generator_expected": 62.802459716796875,
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv7_pathrank_e30",
        "label": "SV7 path-rank e30",
        "epoch": 30,
        "source": BASE_RAW / "sv7_pathrank_e30",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "pathology_hard": 50,
        "pathology_expected": 47.55648422241211,
        "pathology_finalized": False,
        "generator_hard": 44,
        "generator_expected": 62.41193389892578,
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv7_pathrank_e33",
        "label": "SV7 path-rank e33",
        "epoch": 33,
        "source": ROOT / "analysis_outputs" / "prism_integrated_pathranklearn_sv7_e33_interim_posthoc_20260828",
        "source_kind": "candidate",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv7_pathrank_e36",
        "label": "SV7 path-rank e36",
        "epoch": 36,
        "source": ROOT / "analysis_outputs" / "prism_integrated_pathranklearn_sv7_e36_interim_posthoc_20260828",
        "source_kind": "candidate",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv7_pathrank_e39",
        "label": "SV7 path-rank e39",
        "epoch": 39,
        "source": ROOT / "analysis_outputs" / "prism_integrated_pathranklearn_sv7_e39_interim_posthoc_20260828",
        "source_kind": "candidate",
        "personal_rank": 2,
        "pathology_rank_policy": "learned",
        "generator_policy": "free, singleton-only protection",
    },
    {
        "id": "sv6_fixed_e33",
        "label": "SV6 fixed-rank8 e33 protected",
        "epoch": 33,
        "source": BASE_RAW / "sv6_e33",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "pathology_finalized": True,
        "generator_hard": 326,
        "generator_expected": 326.0633850097656,
        "generator_policy": "326 unique-coverage protected",
    },
    {
        "id": "sv6_fixed_e42",
        "label": "SV6 fixed-rank8 e42 free",
        "epoch": 42,
        "source": BASE_RAW / "sv6_e42",
        "source_kind": "raw",
        "personal_rank": 2,
        "pathology_rank_policy": "fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "pathology_finalized": True,
        "generator_hard": 144,
        "generator_expected": 152.21578979492188,
        "generator_policy": "free, singleton-only protection",
    },
]


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_rank_record(path: Path, epoch: int) -> dict[str, Any] | None:
    if not path.exists():
        return None
    match = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("epoch") == epoch:
                match = record
    return match


def nested(data: dict[str, Any] | None, *keys: str) -> Any:
    value: Any = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def ci_value(record: dict[str, Any] | None, index: int) -> Any:
    if not record:
        return None
    ci = record.get("donor_bootstrap_95ci") or record.get("donor_bootstrap_spearman_95ci")
    return ci[index] if isinstance(ci, list) and len(ci) == 2 else None


def leakage_record(leak: dict[str, Any] | None, group: str, target: str) -> dict[str, Any] | None:
    record = nested(leak, group, target)
    if isinstance(record, dict) and "val" in record:
        record = record["val"]
    return record if isinstance(record, dict) else None


def donor_target_record(leak: dict[str, Any] | None, target: str, representation: str) -> dict[str, Any] | None:
    record = nested(leak, "donor_celltype_targets", target, representation)
    if isinstance(record, dict) and "val" in record:
        record = record["val"]
    return record if isinstance(record, dict) else None


def load_source(item: dict[str, Any]) -> dict[str, Any]:
    source = item["source"]
    if item["source_kind"] == "raw":
        names = {
            "module": "module.json",
            "leakage": "leakage.json",
            "bycell": "leakage_by_celltype.json",
            "history": "history_summary.json",
            "strict": "strict_reference.json",
        }
    else:
        names = {
            "module": "module_disease_val_bal40.json",
            "leakage": "leakage_val.json",
            "bycell": "leakage_by_celltype_val.json",
            "history": "history/interim_summary.json",
            "strict": "strict_reference_leakage_val.json",
        }
    data = {key: read_json(source / filename) for key, filename in names.items()}
    if data["module"] is None:
        raise FileNotFoundError(source / names["module"])
    if nested(data["module"], "evaluation_provenance", "split") != "val":
        raise RuntimeError(f"non-validation module source: {item['id']}")
    if data["leakage"] is not None and data["leakage"].get("test_used") is not False:
        raise RuntimeError(f"test leakage source detected: {item['id']}")
    return data


def complete_dynamic_metadata(item: dict[str, Any]) -> None:
    if item["source_kind"] != "candidate":
        return
    epoch = item["epoch"]
    generator = read_json(item["source"] / f"joint_generator_count_epoch_{epoch:03d}.json")
    rank = read_rank_record(item["source"] / "pathology_correction_rank_snapshot.jsonl", epoch)
    if generator is None or rank is None:
        raise FileNotFoundError(f"missing rank/generator evidence for {item['id']}")
    learned = rank["learned_rank"]
    item["generator_hard"] = generator["hard_active"]
    item["generator_expected"] = generator["expected_active"]
    item["pathology_hard"] = learned["hard_rank"]
    item["pathology_expected"] = learned["expected_rank"]
    item["pathology_finalized"] = learned["finalized"]
    item["pathology_uncertain_fraction"] = learned["uncertain_fraction"]
    item["pathology_rank95_axes"] = [axis["rank_95_energy"] for axis in rank["axes"]]
    item["pathology_participation_axes"] = [axis["participation_rank"] for axis in rank["axes"]]


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if finite(row.get(key)) is None else finite(row.get(key)) for key in fields})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for item in CHECKPOINTS:
        complete_dynamic_metadata(item)
    data = {item["id"]: load_source(item) for item in CHECKPOINTS}
    baseline_per_cell = data["rank2_e20"]["module"]["per_celltype"]

    global_rows: list[dict[str, Any]] = []
    for item in CHECKPOINTS:
        bundle = data[item["id"]]
        module = bundle["module"]
        leakage = bundle["leakage"]
        bycell = bundle["bycell"] or {}
        history = bundle["history"] or {}
        latest = history.get("latest", {})
        per_cell = module["per_celltype"]
        per_donor = module.get("per_donor", [])
        ad_values = {name: record.get("AD_module_spearman") for name, record in per_cell.items()}
        non_neuronal = [ad_values[name] for name in NON_NEURONAL]
        z_sex_values = [
            nested(record, "z", "sex", "balanced_accuracy")
            for record in bycell.get("per_celltype", {}).values()
        ]
        z_sex_values = [value for value in z_sex_values if value is not None]
        p_sex_values = [
            nested(record, "personal_code", "sex", "balanced_accuracy")
            for record in bycell.get("per_celltype", {}).values()
        ]
        p_sex_values = [value for value in p_sex_values if value is not None]
        z_sex_shuffle_gaps = [
            nested(record, "z", "sex", "balanced_accuracy")
            - nested(record, "z", "sex", "donor_label_shuffle_bacc")
            for record in bycell.get("per_celltype", {}).values()
            if nested(record, "z", "sex", "balanced_accuracy") is not None
            and nested(record, "z", "sex", "donor_label_shuffle_bacc") is not None
        ]
        p_sex_shuffle_gaps = [
            nested(record, "personal_code", "sex", "balanced_accuracy")
            - nested(record, "personal_code", "sex", "donor_label_shuffle_bacc")
            for record in bycell.get("per_celltype", {}).values()
            if nested(record, "personal_code", "sex", "balanced_accuracy") is not None
            and nested(record, "personal_code", "sex", "donor_label_shuffle_bacc") is not None
        ]
        sex = leakage_record(leakage, "cell_level_nuisance", "sex")
        technology = leakage_record(leakage, "cell_level_nuisance", "technology")
        celltype = leakage_record(leakage, "cell_level_nuisance", "celltype")
        region = leakage_record(leakage, "cell_level_nuisance", "region")
        umi = leakage_record(leakage, "cell_level_continuous_nuisance", "log1p_UMI")
        detected = leakage_record(leakage, "cell_level_continuous_nuisance", "log1p_detected_genes")
        z_sex = donor_target_record(leakage, "sex", "z")
        p_sex = donor_target_record(leakage, "sex", "personal_code")
        z_adnc = donor_target_record(leakage, "ADNC_report_only", "z")
        row = {
            **{key: value for key, value in item.items() if key not in {"source", "source_kind"}},
            "module_cells": nested(module, "evaluation_provenance", "selected_cells"),
            "module_celltypes": module.get("celltypes"),
            "module_count": module.get("modules"),
            "ad_module": module.get("median_AD_module_spearman"),
            "pooled": module.get("pooled_spearman"),
            "blind_centered": module.get("blind_module_centered_median"),
            "module_to_adnc_lam100": nested(module, "module_to_ADNC", "model_lam100"),
            "disease_recon": statistics.mean(record["recon_disease"] for record in per_donor) if per_donor else None,
            "celltypes_better_than_rank2": sum(
                value is not None and value > baseline_per_cell[name]["AD_module_spearman"]
                for name, value in ad_values.items()
            ),
            "celltype_ad_min": min(ad_values.values()),
            "non_neuronal_mean": statistics.mean(non_neuronal),
            "non_neuronal_min": min(non_neuronal),
            "val_rec_nll": latest.get("val:loss/rec"),
            "val_branch_nll": latest.get("val:loss/prism_branch_nll"),
            "val_nonzero_acc": latest.get("val:metric/ordinal_nonzero_acc"),
            "val_within1": latest.get("val:metric/ordinal_nonzero_within1"),
            "val_balanced_recall": latest.get("val:metric/ordinal_balanced_recall"),
            "val_z_participation": latest.get("val:metric/v20_z_participation"),
            "leak_train_cells": nested(leakage, "cell_counts", "train"),
            "leak_val_cells": nested(leakage, "cell_counts", "val"),
            "leak_test_cells": nested(leakage, "cell_counts", "test"),
            "test_used": nested(leakage, "test_used"),
            "sex_bacc": nested(sex, "balanced_accuracy"),
            "sex_ci_low": ci_value(sex, 0),
            "sex_ci_high": ci_value(sex, 1),
            "technology_bacc": nested(technology, "balanced_accuracy"),
            "technology_ci_low": ci_value(technology, 0),
            "technology_ci_high": ci_value(technology, 1),
            "celltype_bacc": nested(celltype, "balanced_accuracy"),
            "region_bacc": nested(region, "balanced_accuracy"),
            "umi_spearman": nested(umi, "spearman"),
            "detected_genes_spearman": nested(detected, "spearman"),
            "donor_celltype_z_sex_bacc": nested(z_sex, "balanced_accuracy"),
            "donor_celltype_personal_sex_bacc": nested(p_sex, "balanced_accuracy"),
            "donor_celltype_z_adnc_spearman": nested(z_adnc, "spearman"),
            "celltype_z_sex_median": statistics.median(z_sex_values) if z_sex_values else None,
            "celltype_z_sex_max": max(z_sex_values) if z_sex_values else None,
            "celltype_personal_sex_median": statistics.median(p_sex_values) if p_sex_values else None,
            "celltype_z_sex_shuffle_gap_median": statistics.median(z_sex_shuffle_gaps) if z_sex_shuffle_gaps else None,
            "celltype_z_sex_shuffle_gap_positive_count": sum(value > 0 for value in z_sex_shuffle_gaps),
            "celltype_z_sex_shuffle_gap_ge_010_count": sum(value >= 0.1 for value in z_sex_shuffle_gaps),
            "celltype_personal_sex_shuffle_gap_median": statistics.median(p_sex_shuffle_gaps) if p_sex_shuffle_gaps else None,
            "celltype_personal_sex_shuffle_gap_positive_count": sum(value > 0 for value in p_sex_shuffle_gaps),
            "celltype_personal_sex_shuffle_gap_ge_010_count": sum(value >= 0.1 for value in p_sex_shuffle_gaps),
            "z_representation_participation": nested(leakage, "representation", "z", "participation_ratio"),
            "personal_representation_participation": nested(leakage, "representation", "personal_code", "participation_ratio"),
            "z_donor_retrieval_top1": nested(leakage, "donor_retrieval", "val", "z_cross_celltype", "top1_accuracy"),
            "personal_donor_retrieval_top1": nested(leakage, "donor_retrieval", "val", "personal_code_cross_celltype", "top1_accuracy"),
            "strict_val_donors": nested(bundle["strict"], "donor_counts", "val"),
        }
        if item.get("pathology_participation_axes"):
            row["pathology_participation_min"] = min(item["pathology_participation_axes"])
            row["pathology_participation_max"] = max(item["pathology_participation_axes"])
        global_rows.append(row)

    global_fields = list(global_rows[0].keys())
    for row in global_rows[1:]:
        for key in row:
            if key not in global_fields:
                global_fields.append(key)
    write_tsv(OUT / "global_candidate_comparison.tsv", global_rows, global_fields)

    celltypes = sorted(baseline_per_cell)
    module_rows: list[dict[str, Any]] = []
    leakage_rows: list[dict[str, Any]] = []
    for name in celltypes:
        module_row: dict[str, Any] = {"celltype": name}
        leak_row: dict[str, Any] = {"celltype": name}
        for item in CHECKPOINTS:
            source = data[item["id"]]
            module_record = nested(source["module"], "per_celltype", name) or {}
            module_row[f"{item['id']}_ad_module"] = module_record.get("AD_module_spearman")
            module_row[f"{item['id']}_blind_centered"] = module_record.get("blind_centered")
            cell_leak = nested(source["bycell"], "per_celltype", name) or {}
            leak_row[f"{item['id']}_z_sex_bacc"] = nested(cell_leak, "z", "sex", "balanced_accuracy")
            leak_row[f"{item['id']}_personal_sex_bacc"] = nested(cell_leak, "personal_code", "sex", "balanced_accuracy")
        module_rows.append(module_row)
        leakage_rows.append(leak_row)
    write_tsv(OUT / "celltype_module_comparison.tsv", module_rows, list(module_rows[0].keys()))
    write_tsv(OUT / "celltype_sex_leakage_comparison.tsv", leakage_rows, list(leakage_rows[0].keys()))

    manifest = {
        "schema_version": "kmlee_bam.sv7_pathrank_validation_candidate_comparison.v1",
        "official_test_used": False,
        "official_test_status": "consumed_elsewhere_not_reused",
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
            "recent_candidates_per_donor_celltype": 30,
            "seed": 20260814,
        },
        "checkpoints": [
            {key: value for key, value in item.items() if key not in {"source", "source_kind"}}
            for item in CHECKPOINTS
        ],
        "outputs": [
            "global_candidate_comparison.tsv",
            "celltype_module_comparison.tsv",
            "celltype_sex_leakage_comparison.tsv",
            "comparison_manifest.json",
            "ANALYSIS_CUTOFF_LOCK.json",
            "SV7_PATHRANK_INTERIM_POSTHOC_KO.md",
            "INTERIM_SELECTION_DECISION.json",
            "pathology_rank_threshold_sensitivity.json",
        ],
    }
    with (OUT / "comparison_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
