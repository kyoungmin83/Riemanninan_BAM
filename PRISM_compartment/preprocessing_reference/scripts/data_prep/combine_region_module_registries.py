#!/usr/bin/env python
"""Combine region-specific final module registries into one training registry.

Common pathway modules (Hallmark/Reactome/singletons) are kept once.
Region-specific hdWGCNA residual modules are kept separately, e.g.
hdWGCNA_DLPFC::* and hdWGCNA_MTG::*.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


def _load_json(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _is_hdwgcna(name: str, source: str) -> bool:
    return name.lower().startswith("hdwgcna") or source.lower().startswith("hdwgcna")


def _row_indices(row: list[float]) -> set[int]:
    return {i for i, v in enumerate(row) if float(v) > 0.0}


def _l2_binary_weight(indices: set[int], n_genes: int) -> np.ndarray:
    w = np.zeros(n_genes, dtype=np.float32)
    if indices:
        val = 1.0 / float(np.sqrt(len(indices)))
        for i in indices:
            w[i] = val
    return w


def _map_weight(
    weight_row: np.ndarray,
    src_gene_names: list[str],
    dst_gene_to_idx: dict[str, int],
    n_dst: int,
) -> np.ndarray:
    out = np.zeros(n_dst, dtype=np.float32)
    for old_i, val in enumerate(weight_row):
        if val == 0:
            continue
        new_i = dst_gene_to_idx.get(src_gene_names[old_i])
        if new_i is not None:
            out[new_i] = float(val)
    norm = float(np.sqrt(np.sum(out * out)))
    if norm > 0:
        out /= norm
    return out


def combine(cfg: dict[str, Any]) -> dict[str, Any]:
    regions = cfg["regions"]
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    loaded: list[dict[str, Any]] = []
    gene_union: set[str] = set()
    for item in regions:
        reg = _load_json(item["registry_json"])
        act = np.load(item["activity_weight_npz"], allow_pickle=True)["activity_weight"]
        if act.shape != (len(reg["module_names"]), len(reg["gene_names"])):
            raise ValueError(
                f"activity_weight shape mismatch for {item['region']}: "
                f"{act.shape} vs registry {(len(reg['module_names']), len(reg['gene_names']))}"
            )
        loaded.append({"meta": item, "registry": reg, "activity": act})
        gene_union.update(str(g) for g in reg["gene_names"])

    gene_names = sorted(gene_union)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    common: dict[str, dict[str, Any]] = {}
    common_order: list[str] = []
    region_modules: list[dict[str, Any]] = []

    for blob in loaded:
        region = str(blob["meta"]["region"])
        reg = blob["registry"]
        act = blob["activity"]
        src_genes = [str(g) for g in reg["gene_names"]]

        for mi, (name, source, token_type, row) in enumerate(
            zip(
                reg["module_names"],
                reg["module_sources"],
                reg["module_token_types"],
                reg["membership_binary"],
            )
        ):
            name = str(name)
            source = str(source)
            token_type = str(token_type)
            old_idx = _row_indices(row)
            new_idx = {gene_to_idx[src_genes[i]] for i in old_idx}

            if _is_hdwgcna(name, source):
                region_modules.append(
                    {
                        "name": name,
                        "source": source,
                        "token_type": token_type,
                        "gene_indices": new_idx,
                        "activity": _map_weight(act[mi], src_genes, gene_to_idx, n_genes),
                        "regions": [region],
                    }
                )
                continue

            if name not in common:
                common[name] = {
                    "name": name,
                    "source": source,
                    "token_type": token_type,
                    "gene_indices": set(),
                    "regions": [],
                }
                common_order.append(name)
            entry = common[name]
            if entry["source"] != source or entry["token_type"] != token_type:
                raise ValueError(f"Common module metadata mismatch for {name}")
            entry["gene_indices"].update(new_idx)
            if region not in entry["regions"]:
                entry["regions"].append(region)

    modules: list[dict[str, Any]] = []
    for name in common_order:
        entry = common[name]
        entry["activity"] = _l2_binary_weight(entry["gene_indices"], n_genes)
        modules.append(entry)

    modules.extend(sorted(region_modules, key=lambda m: (m["source"], m["name"])))

    module_names = [m["name"] for m in modules]
    module_sources = [m["source"] for m in modules]
    module_token_types = [m["token_type"] for m in modules]
    membership_binary: list[list[float]] = []
    activity = np.zeros((len(modules), n_genes), dtype=np.float32)

    for mi, mod in enumerate(modules):
        row = [0.0] * n_genes
        for gi in sorted(mod["gene_indices"]):
            row[int(gi)] = 1.0
        membership_binary.append(row)
        activity[mi] = mod["activity"]

    singleton_gene_indices = sorted(
        {
            int(gi)
            for mod in modules
            if mod["token_type"] == "singleton"
            for gi in mod["gene_indices"]
        }
    )
    residual_gene_indices = sorted(
        {
            gene_to_idx[blob["registry"]["gene_names"][int(i)]]
            for blob in loaded
            for i in blob["registry"].get("residual_gene_indices", [])
            if int(i) < len(blob["registry"]["gene_names"])
        }
    )

    registry = {
        "gene_names": gene_names,
        "module_names": module_names,
        "module_sources": module_sources,
        "module_token_types": module_token_types,
        "membership_binary": membership_binary,
        "singleton_gene_indices": singleton_gene_indices,
        "residual_gene_indices": residual_gene_indices,
    }

    compact_path = out_dir / "final_gene_module_registry_compact.json"
    with open(compact_path, "w", encoding="utf-8") as f:
        json.dump(registry, f, ensure_ascii=False, indent=2)

    np.savez_compressed(
        out_dir / "activity_weight_kme_or_l2_membership.npz",
        activity_weight=activity,
        gene_names=np.array(gene_names, dtype=object),
        module_names=np.array(module_names, dtype=object),
        registry_json_path=str(compact_path),
        note=(
            "DLPFC+MTG combined registry: common pathway modules merged once; "
            "region-specific hdWGCNA modules retained with signed kME weights."
        ),
    )
    with open(out_dir / "gene_union_ensembl.json", "w", encoding="utf-8") as f:
        json.dump(gene_names, f, ensure_ascii=False, indent=2)

    with open(out_dir / "final_module_table.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["module_name", "module_source", "token_type", "n_genes", "regions"],
        )
        writer.writeheader()
        for mod in modules:
            writer.writerow(
                {
                    "module_name": mod["name"],
                    "module_source": mod["source"],
                    "token_type": mod["token_type"],
                    "n_genes": len(mod["gene_indices"]),
                    "regions": "|".join(mod["regions"]),
                }
            )

    summary = {
        "regions": [str(x["region"]) for x in regions],
        "input_registry_jsons": [str(x["registry_json"]) for x in regions],
        "n_genes": n_genes,
        "n_modules": len(modules),
        "n_common_modules": len(common_order),
        "n_region_hdwgcna_modules": len(region_modules),
        "n_singleton_gene_indices": len(singleton_gene_indices),
        "n_residual_gene_indices": len(residual_gene_indices),
        "compact_registry_path": str(compact_path),
        "activity_weight_path": str(out_dir / "activity_weight_kme_or_l2_membership.npz"),
    }
    with open(out_dir / "combined_registry_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = _load_json(args.config)
    summary = combine(cfg)
    print("[combine-region-registry] done", flush=True)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
