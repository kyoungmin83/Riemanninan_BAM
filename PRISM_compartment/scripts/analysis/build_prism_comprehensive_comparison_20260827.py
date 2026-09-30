#!/usr/bin/env python3
"""Build the validation-only PRISM checkpoint comparison tables."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis_outputs" / "prism_comprehensive_checkpoint_comparison_20260827"
RAW = OUT / "raw"

CHECKPOINTS = [
    {
        "id": "rank2_e20",
        "label": "Frozen rank2 e20",
        "server": "baseline",
        "epoch": 20,
        "lineage": "frozen CLS personal-rank2",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "fixed mask",
        "generator_hard": 97,
        "generator_expected": 97.0,
        "generator_protected": 97,
    },
    {
        "id": "sv7_ctemb128_e14",
        "label": "SV7 ctemb128 e14",
        "server": "SV7",
        "epoch": 14,
        "lineage": "Phase-II ctemb128",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 99,
        "generator_expected": 100.77603912353516,
        "generator_protected": 10,
    },
    {
        "id": "sv7_ctemb128_e19",
        "label": "SV7 ctemb128 e19",
        "server": "SV7",
        "epoch": 19,
        "lineage": "Phase-II ctemb128",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 94,
        "generator_expected": 94.16015625,
        "generator_protected": 10,
    },
    {
        "id": "sv7_ctemb128_e26",
        "label": "SV7 ctemb128 e26",
        "server": "SV7",
        "epoch": 26,
        "lineage": "Phase-II ctemb128",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 90,
        "generator_expected": 90.20425415039062,
        "generator_protected": 10,
    },
    {
        "id": "sv7_pathrank_e14",
        "label": "SV7 path-rank e14",
        "server": "SV7",
        "epoch": 14,
        "lineage": "integrated learnable pathology-rank",
        "personal_rank": "2 fixed",
        "pathology_rank": "learned, not finalized",
        "pathology_hard": 96,
        "pathology_expected": 73.95681762695312,
        "generator_policy": "326 unique-coverage protected",
        "generator_hard": 414,
        "generator_expected": 408.3619384765625,
        "generator_protected": 326,
    },
    {
        "id": "sv7_pathrank_e26",
        "label": "SV7 path-rank e26",
        "server": "SV7",
        "epoch": 26,
        "lineage": "integrated learnable pathology-rank",
        "personal_rank": "2 fixed",
        "pathology_rank": "learned, not finalized",
        "pathology_hard": 72,
        "pathology_expected": 50.01966094970703,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 26,
        "generator_expected": 62.802459716796875,
        "generator_protected": 10,
    },
    {
        "id": "sv7_pathrank_e30",
        "label": "SV7 path-rank e30",
        "server": "SV7",
        "epoch": 30,
        "lineage": "integrated learnable pathology-rank",
        "personal_rank": "2 fixed",
        "pathology_rank": "learned, not finalized",
        "pathology_hard": 50,
        "pathology_expected": 47.55648422241211,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 44,
        "generator_expected": 62.41193389892578,
        "generator_protected": 10,
    },
    {
        "id": "sv6_e33",
        "label": "SV6 protected e33",
        "server": "SV6",
        "epoch": 33,
        "lineage": "integrated fixed pathology-rank8",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "326 unique-coverage protected",
        "generator_hard": 326,
        "generator_expected": 326.0633850097656,
        "generator_protected": 326,
    },
    {
        "id": "sv6_e42",
        "label": "SV6 free e42",
        "server": "SV6",
        "epoch": 42,
        "lineage": "integrated fixed pathology-rank8",
        "personal_rank": "2 fixed",
        "pathology_rank": "8 fixed",
        "pathology_hard": 8,
        "pathology_expected": 8.0,
        "generator_policy": "free, singleton-only protection",
        "generator_hard": 144,
        "generator_expected": 152.21578979492188,
        "generator_protected": 10,
    },
]


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
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


def donor_target_record(
    leak: dict[str, Any] | None, target: str, representation: str
) -> dict[str, Any] | None:
    record = nested(leak, "donor_celltype_targets", target, representation)
    if isinstance(record, dict) and "val" in record:
        record = record["val"]
    return record if isinstance(record, dict) else None


def history_value(history: dict[str, Any] | None, key: str) -> Any:
    return finite(nested(history, "latest", key))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if finite(row.get(key)) is None else finite(row.get(key)) for key in fields})


def fmt(value: Any, digits: int = 4) -> str:
    value = finite(value)
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def markdown_table(headers: list[str], rows: list[list[Any]], digits: int = 4) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] + ["---:"] * (len(headers) - 1)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(value, digits) for value in row) + " |")
    return "\n".join(lines)


def main() -> None:
    data: dict[str, dict[str, Any]] = {}
    for checkpoint in CHECKPOINTS:
        raw = RAW / checkpoint["id"]
        data[checkpoint["id"]] = {
            "module": read_json(raw / "module.json"),
            "leakage": read_json(raw / "leakage.json"),
            "leakage_by_celltype": read_json(raw / "leakage_by_celltype.json"),
            "history": read_json(raw / "history_summary.json"),
        }

    global_rows: list[dict[str, Any]] = []
    for checkpoint in CHECKPOINTS:
        item = data[checkpoint["id"]]
        module = item["module"]
        leak = item["leakage"]
        history = item["history"]
        sex = leakage_record(leak, "cell_level_nuisance", "sex")
        tech = leakage_record(leak, "cell_level_nuisance", "technology")
        celltype = leakage_record(leak, "cell_level_nuisance", "celltype")
        region = leakage_record(leak, "cell_level_nuisance", "region")
        umi = leakage_record(leak, "cell_level_continuous_nuisance", "log1p_UMI")
        detected = leakage_record(leak, "cell_level_continuous_nuisance", "log1p_detected_genes")
        z_sex = donor_target_record(leak, "sex", "z")
        p_sex = donor_target_record(leak, "sex", "personal_code")
        z_age = donor_target_record(leak, "age", "z")
        p_age = donor_target_record(leak, "age", "personal_code")
        z_adnc = donor_target_record(leak, "ADNC_report_only", "z")
        p_adnc = donor_target_record(leak, "ADNC_report_only", "personal_code")
        per_donor = module.get("per_donor", []) if module else []
        row = dict(checkpoint)
        row.update(
            {
                "module_selected_cells": nested(module, "evaluation_provenance", "selected_cells"),
                "module_celltypes": nested(module, "celltypes"),
                "module_count": nested(module, "modules"),
                "ad_module_median": nested(module, "median_AD_module_spearman"),
                "pooled_spearman": nested(module, "pooled_spearman"),
                "blind_raw_median": nested(module, "blind_module_raw_median"),
                "blind_centered_median": nested(module, "blind_module_centered_median"),
                "module_to_adnc_lam30": nested(module, "module_to_ADNC", "model_lam30"),
                "module_to_adnc_lam100": nested(module, "module_to_ADNC", "model_lam100"),
                "module_to_adnc_lam300": nested(module, "module_to_ADNC", "model_lam300"),
                "disease_recon_mean": statistics.mean(x["recon_disease"] for x in per_donor) if per_donor else None,
                "full_recon_mean": statistics.mean(x["recon_full"] for x in per_donor) if per_donor else None,
                "per_module_bias_mean": nested(module, "per_module_bias_mean"),
                "module_sexdis_tested": nested(module, "sexdis_summary", "n_tested"),
                "module_sexdis_significant": nested(module, "sexdis_summary", "n_signif_p05"),
                "val_rec_nll": history_value(history, "val:loss/rec"),
                "val_branch_nll": history_value(history, "val:loss/prism_branch_nll"),
                "val_nonzero_acc": history_value(history, "val:metric/ordinal_nonzero_acc"),
                "val_within1": history_value(history, "val:metric/ordinal_nonzero_within1"),
                "val_balanced_recall": history_value(history, "val:metric/ordinal_balanced_recall"),
                "val_z_participation": history_value(history, "val:metric/v20_z_participation"),
                "leak_train_cells": nested(leak, "cell_counts", "train"),
                "leak_val_cells": nested(leak, "cell_counts", "val"),
                "leak_test_cells": nested(leak, "cell_counts", "test"),
                "test_used": nested(leak, "test_used"),
                "sex_bacc": nested(sex, "balanced_accuracy"),
                "sex_ci_low": ci_value(sex, 0),
                "sex_ci_high": ci_value(sex, 1),
                "technology_bacc": nested(tech, "balanced_accuracy"),
                "technology_ci_low": ci_value(tech, 0),
                "technology_ci_high": ci_value(tech, 1),
                "celltype_bacc": nested(celltype, "balanced_accuracy"),
                "celltype_ci_low": ci_value(celltype, 0),
                "celltype_ci_high": ci_value(celltype, 1),
                "region_bacc": nested(region, "balanced_accuracy"),
                "region_ci_low": ci_value(region, 0),
                "region_ci_high": ci_value(region, 1),
                "umi_spearman": nested(umi, "spearman"),
                "detected_genes_spearman": nested(detected, "spearman"),
                "donor_celltype_z_sex_bacc": nested(z_sex, "balanced_accuracy"),
                "donor_celltype_z_sex_ci_low": ci_value(z_sex, 0),
                "donor_celltype_z_sex_ci_high": ci_value(z_sex, 1),
                "donor_celltype_personal_sex_bacc": nested(p_sex, "balanced_accuracy"),
                "donor_celltype_z_age_spearman": nested(z_age, "spearman"),
                "donor_celltype_personal_age_spearman": nested(p_age, "spearman"),
                "donor_celltype_z_adnc_spearman": nested(z_adnc, "spearman"),
                "donor_celltype_personal_adnc_spearman": nested(p_adnc, "spearman"),
                "z_donor_retrieval_top1": nested(leak, "donor_retrieval", "val", "z_cross_celltype", "top1_accuracy"),
                "personal_donor_retrieval_top1": nested(leak, "donor_retrieval", "val", "personal_code_cross_celltype", "top1_accuracy"),
            }
        )
        global_rows.append(row)

    global_fields = list(global_rows[0].keys())
    write_tsv(OUT / "checkpoint_global_comparison.tsv", global_rows, global_fields)

    celltypes = sorted(
        {
            celltype
            for checkpoint in CHECKPOINTS
            for celltype in (data[checkpoint["id"]]["module"] or {}).get("per_celltype", {})
        }
    )
    module_wide_rows: list[dict[str, Any]] = []
    module_long_rows: list[dict[str, Any]] = []
    for celltype_name in celltypes:
        wide: dict[str, Any] = {"celltype": celltype_name}
        baseline_value = nested(data["rank2_e20"]["module"], "per_celltype", celltype_name, "AD_module_spearman")
        for checkpoint in CHECKPOINTS:
            record = nested(data[checkpoint["id"]]["module"], "per_celltype", celltype_name) or {}
            ad_value = record.get("AD_module_spearman")
            wide[checkpoint["id"]] = ad_value
            wide[f"{checkpoint['id']}_minus_rank2"] = ad_value - baseline_value if ad_value is not None and baseline_value is not None else None
            sexdis = record.get("sexdis", {})
            module_long_rows.append(
                {
                    "checkpoint_id": checkpoint["id"],
                    "checkpoint_label": checkpoint["label"],
                    "celltype": celltype_name,
                    "n_donors": record.get("n_donors"),
                    "ad_module_spearman": ad_value,
                    "blind_raw": record.get("blind_raw"),
                    "blind_centered": record.get("blind_centered"),
                    "cross_disjoint": record.get("cross_disjoint"),
                    "obs_ceiling": record.get("obs_ceiling"),
                    "model_ceiling": record.get("model_ceiling"),
                    "recovery_shared": record.get("recovery_shared"),
                    "sexdis_n_male": sexdis.get("n_M"),
                    "sexdis_n_female": sexdis.get("n_F"),
                    "sexdis_obs_mvf_spearman": sexdis.get("obs_MvF_spearman"),
                    "sexdis_perm_null_median": sexdis.get("perm_null_median"),
                    "sexdis_pval": sexdis.get("pval"),
                    "sexdis_model_mvf_spearman": sexdis.get("model_MvF_spearman"),
                }
            )
        module_wide_rows.append(wide)

    module_wide_fields = ["celltype"]
    for checkpoint in CHECKPOINTS:
        module_wide_fields.extend([checkpoint["id"], f"{checkpoint['id']}_minus_rank2"])
    write_tsv(OUT / "celltype_ad_module_recovery_wide.tsv", module_wide_rows, module_wide_fields)
    write_tsv(OUT / "celltype_module_metrics_long.tsv", module_long_rows, list(module_long_rows[0].keys()))

    leakage_wide_rows: list[dict[str, Any]] = []
    leakage_long_rows: list[dict[str, Any]] = []
    continuous_targets = ["age", "ADNC_report_only", "thal", "braak", "cerad", "late", "lewy", "support_context_count", "support_context_reliability"]
    for celltype_name in celltypes:
        wide: dict[str, Any] = {"celltype": celltype_name}
        for checkpoint in CHECKPOINTS:
            source = nested(data[checkpoint["id"]]["leakage_by_celltype"], "per_celltype", celltype_name) or {}
            z_sex_record = nested(source, "z", "sex") or {}
            p_sex_record = nested(source, "personal_code", "sex") or {}
            wide[f"{checkpoint['id']}_z_sex_bacc"] = z_sex_record.get("balanced_accuracy")
            wide[f"{checkpoint['id']}_personal_sex_bacc"] = p_sex_record.get("balanced_accuracy")
            row: dict[str, Any] = {
                "checkpoint_id": checkpoint["id"],
                "checkpoint_label": checkpoint["label"],
                "celltype": celltype_name,
                "n_train_donors": source.get("n_train_donors"),
                "n_val_donors": source.get("n_val_donors"),
                "z_sex_bacc": z_sex_record.get("balanced_accuracy"),
                "z_sex_ci_low": ci_value(z_sex_record, 0),
                "z_sex_ci_high": ci_value(z_sex_record, 1),
                "z_sex_auc": z_sex_record.get("roc_auc"),
                "personal_sex_bacc": p_sex_record.get("balanced_accuracy"),
                "personal_sex_ci_low": ci_value(p_sex_record, 0),
                "personal_sex_ci_high": ci_value(p_sex_record, 1),
                "personal_sex_auc": p_sex_record.get("roc_auc"),
            }
            for representation in ("z", "personal_code"):
                for target in continuous_targets:
                    record = nested(source, representation, target) or {}
                    prefix = f"{representation}_{target.lower()}"
                    row[f"{prefix}_spearman"] = record.get("spearman")
                    row[f"{prefix}_ci_low"] = ci_value(record, 0)
                    row[f"{prefix}_ci_high"] = ci_value(record, 1)
            leakage_long_rows.append(row)
        leakage_wide_rows.append(wide)

    leakage_wide_fields = ["celltype"]
    for checkpoint in CHECKPOINTS:
        leakage_wide_fields.extend([f"{checkpoint['id']}_z_sex_bacc", f"{checkpoint['id']}_personal_sex_bacc"])
    write_tsv(OUT / "celltype_sex_leakage_wide.tsv", leakage_wide_rows, leakage_wide_fields)
    write_tsv(OUT / "celltype_leakage_metrics_long.tsv", leakage_long_rows, list(leakage_long_rows[0].keys()))

    manifest = {
        "schema_version": "prism.comprehensive_validation_checkpoint_comparison.v1",
        "official_test_used": False,
        "module_protocol": {
            "split": "validation",
            "donors": 16,
            "celltypes": 24,
            "modules": 404,
            "per_donor_celltype": 40,
            "sample_seed": 20260814,
        },
        "leakage_protocol_note": "Historical rank2/ctemb128 extracts use cap15 (22,702 train, 5,688 val); recent integrated extracts use cap30 (44,960 train, 11,317 val). Missing leakage was not imputed.",
        "checkpoints": CHECKPOINTS,
        "outputs": [
            "checkpoint_global_comparison.tsv",
            "celltype_ad_module_recovery_wide.tsv",
            "celltype_module_metrics_long.tsv",
            "celltype_sex_leakage_wide.tsv",
            "celltype_leakage_metrics_long.tsv",
            "COMPREHENSIVE_COMPARISON_KO.md",
        ],
    }
    with (OUT / "comparison_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    global_by_id = {row["id"]: row for row in global_rows}
    global_headers = ["Checkpoint", "Gen hard/expected", "Path hard/expected", "AD median", "Pooled", "Blind-centered", "Module→ADNC (ridge 100)", "Disease recon"]
    global_table = [
        [
            checkpoint["label"],
            f"{global_by_id[checkpoint['id']]['generator_hard']}/{global_by_id[checkpoint['id']]['generator_expected']:.2f}",
            f"{global_by_id[checkpoint['id']]['pathology_hard']}/{global_by_id[checkpoint['id']]['pathology_expected']:.2f}",
            global_by_id[checkpoint["id"]]["ad_module_median"],
            global_by_id[checkpoint["id"]]["pooled_spearman"],
            global_by_id[checkpoint["id"]]["blind_centered_median"],
            global_by_id[checkpoint["id"]]["module_to_adnc_lam100"],
            global_by_id[checkpoint["id"]]["disease_recon_mean"],
        ]
        for checkpoint in CHECKPOINTS
    ]
    leak_headers = ["Checkpoint", "Train/val cells", "Sex", "Technology", "Cell type", "Region", "DC z→sex", "CT-sex median"]
    leak_table: list[list[Any]] = []
    for checkpoint in CHECKPOINTS:
        row = global_by_id[checkpoint["id"]]
        byct = data[checkpoint["id"]]["leakage_by_celltype"]
        z_sex_values = [
            nested(record, "z", "sex", "balanced_accuracy")
            for record in (byct or {}).get("per_celltype", {}).values()
        ]
        z_sex_values = [value for value in z_sex_values if value is not None]
        leak_table.append(
            [
                checkpoint["label"],
                f"{row['leak_train_cells']}/{row['leak_val_cells']}" if row["leak_train_cells"] is not None else "NA",
                row["sex_bacc"],
                row["technology_bacc"],
                row["celltype_bacc"],
                row["region_bacc"],
                row["donor_celltype_z_sex_bacc"],
                statistics.median(z_sex_values) if z_sex_values else None,
            ]
        )

    module_celltype_headers = ["Cell type"] + [checkpoint["label"] for checkpoint in CHECKPOINTS]
    module_celltype_table = [
        [row["celltype"]] + [row[checkpoint["id"]] for checkpoint in CHECKPOINTS]
        for row in module_wide_rows
    ]
    module_winner_counts = {checkpoint["id"]: 0 for checkpoint in CHECKPOINTS}
    for row in module_wide_rows:
        best = max(row[checkpoint["id"]] for checkpoint in CHECKPOINTS)
        for checkpoint in CHECKPOINTS:
            if row[checkpoint["id"]] == best:
                module_winner_counts[checkpoint["id"]] += 1
    module_winner_summary = ", ".join(
        f"{checkpoint['label']} {module_winner_counts[checkpoint['id']]}개"
        for checkpoint in sorted(CHECKPOINTS, key=lambda item: module_winner_counts[item["id"]], reverse=True)
        if module_winner_counts[checkpoint["id"]] > 0
    )
    leakage_celltype_headers = ["Cell type"] + [checkpoint["label"] for checkpoint in CHECKPOINTS]
    z_leakage_celltype_table = [
        [row["celltype"]] + [row[f"{checkpoint['id']}_z_sex_bacc"] for checkpoint in CHECKPOINTS]
        for row in leakage_wide_rows
    ]
    personal_leakage_celltype_table = [
        [row["celltype"]] + [row[f"{checkpoint['id']}_personal_sex_bacc"] for checkpoint in CHECKPOINTS]
        for row in leakage_wide_rows
    ]

    report = f"""# PRISM 전체 validation 비교표

- 작성일: 2026-08-27 KST
- 공식 test donor 사용: 없음
- 모듈 비교: {len(CHECKPOINTS)}개 체크포인트 모두 validation 16 donors, 24 cell types, 404 modules, donor-celltype당 최대 40 cells, 동일 seed
- 누출 비교: 과거 rank2/ctemb128은 cap15, 최신 integrated는 cap30이다. 표본 수를 병기했으며 미측정값은 `NA`다.
- `rank2`는 personal rank를 뜻한다. 현재 SV6·SV7도 personal rank는 2이며, SV7에서 학습 중인 rank는 pathology correction rank다.

## 전체 AD 관련 모듈 복원

{markdown_table(global_headers, global_table)}

핵심 판정:

- AD-module 중앙값은 SV6 e42가 가장 높고, module→ADNC는 SV7 ctemb128 e14가 가장 높다.
- 24개 세포유형별 AD-module 1위 수: {module_winner_summary}.
- SV7 path-rank e30은 e26의 blind-centered 및 donor disease-isolated 붕괴에서 크게 회복했고 pooled 복원은 이 표에서 가장 높다.
- e30은 frozen rank2 e20보다 pooled, blind-centered, module→ADNC, disease-isolated 복원이 높지만, path-rank e14의 blind-centered 및 disease-isolated 수준에는 아직 못 미친다.
- Frozen rank2 e20은 생물학과 nuisance의 균형 기준선이다.

## 전체 누출 비교

{markdown_table(leak_headers, leak_table)}

- SV7 path-rank e30은 global sex가 우연 수준이고 세포유형별 z→sex 중앙값도 e26보다 낮아졌지만, UMI·detected-gene depth 및 donor identity 누출은 커졌다.
- SV7 path-rank e14와 ctemb128 e26은 독립 누출 재추론이 없어 `NA`다. 학습 adversary 수치로 대체하지 않았다.
- 누출 CI, UMI, detected genes, age, ADNC, donor retrieval까지의 전체 수치는 `checkpoint_global_comparison.tsv`에 있다.

## 24개 세포유형별 AD 모듈 복원도

{markdown_table(module_celltype_headers, module_celltype_table)}

## 24개 세포유형별 z→sex 누출

{markdown_table(leakage_celltype_headers, z_leakage_celltype_table)}

## 24개 세포유형별 personal-code→sex 누출

{markdown_table(leakage_celltype_headers, personal_leakage_celltype_table)}

## 세부 원자료 표

- `checkpoint_global_comparison.tsv`: 전체 모듈, reconstruction, generator/rank, nuisance, CI, donor-level 질병 probe
- `celltype_ad_module_recovery_wide.tsv`: 세포유형별 AD 모듈 복원 및 rank2 대비 변화
- `celltype_module_metrics_long.tsv`: 세포유형별 AD, blind, cross-disjoint, ceiling, sex-discordance 전체 지표
- `celltype_sex_leakage_wide.tsv`: 세포유형별 z/personal sex 누출
- `celltype_leakage_metrics_long.tsv`: 세포유형별 sex CI/AUC, age, ADNC, 병리축 및 support-context probe
- `comparison_manifest.json`: 체크포인트 계약과 평가 규약
"""
    with (OUT / "COMPREHENSIVE_COMPARISON_KO.md").open("w", encoding="utf-8") as handle:
        handle.write(report)


if __name__ == "__main__":
    main()
