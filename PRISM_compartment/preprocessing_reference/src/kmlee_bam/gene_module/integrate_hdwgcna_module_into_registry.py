from __future__ import annotations

"""
Integrate subclass-specific hdWGCNA residual modules into a final
GeneModuleRegistry JSON.

Why this script exists
----------------------
After build_full_zarr_module_artifacts.py, the initial registry contains:

    regular modules  : Hallmark / Reactome
    singleton modules: selected individual genes
    residual module  : a large RESIDUAL_LEFTOVER pseudo-module

You then ran hdWGCNA on residual HVG subsets per subclass. Those hdWGCNA
modules should replace the large residual pseudo-module as biologically more
specific residual submodules.

This script performs that integration.

Important
---------
This script does NOT create a [cell x module] scalar matrix. For the final
module-token model, the training input should still be gene-level expression
over the final registry gene union. GeneModuleTokenizer then lifts gene tokens
to module tokens during training.

Typical usage
-------------
python integrate_hdwgcna_modules_into_registry.py --config integrate_hdwgcna_cfg.json
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


# =====================================================================
# Config
# =====================================================================
@dataclass
class IntegrateHDWGCNAConfig:
    # ---- inputs / outputs -------------------------------------------
    base_registry_json: str
    hdwgcna_root_dir: str
    out_dir: str

    # Optional symbol map saved by build_full_zarr_module_artifacts.py.
    # Use it if hdWGCNA CSVs contain gene symbols while registry uses ENSG.
    symbol_to_gene_ids_json: Optional[str] = None

    # ---- hdWGCNA CSV discovery --------------------------------------
    modules_glob: str = "**/*_modules_with_kME_self.csv"
    fallback_modules_glob: str = "**/*_modules.csv"

    # Set explicitly if auto-detection fails.
    gene_column: Optional[str] = None
    module_column: Optional[str] = None
    kme_column: str = "kME_self"

    # ---- filtering ---------------------------------------------------
    drop_grey: bool = True
    grey_names: list[str] = field(default_factory=lambda: ["grey", "gray"])
    min_module_size: int = 8
    max_module_size: Optional[int] = 500

    # ---- module naming ----------------------------------------------
    hdwgcna_name_prefix: str = "hdWGCNA"
    hdwgcna_source_prefix: str = "hdwgcna"

    # ---- base registry handling -------------------------------------
    # Recommended: drop the old giant residual pseudo-module.
    drop_existing_residual_modules: bool = True
    drop_existing_hdwgcna_modules: bool = True

    # Usually false. If true, remaining old residual genes not covered by
    # regular/singleton/hdWGCNA modules are added as one fallback module.
    add_leftover_residual_module: bool = False
    leftover_module_name: str = "RESIDUAL_LEFTOVER_AFTER_HDWGCNA"

    # ---- redundancy among imported hdWGCNA modules -------------------
    prune_hdwgcna_redundancy_jaccard: Optional[float] = None

    # ---- gene matching -----------------------------------------------
    allow_direct_gene_name_match: bool = True

    # ---- outputs -----------------------------------------------------
    write_full_registry: bool = True
    write_compact_registry: bool = True
    write_gene_union: bool = True
    write_activity_weight: bool = True

    # ---- safety ------------------------------------------------------
    fail_if_no_hdwgcna_modules: bool = True


# =====================================================================
# Small utilities
# =====================================================================
def _log(msg: str) -> None:
    print(f"[integrate-hdwgcna] {msg}", flush=True)


def load_config(path: str | Path) -> IntegrateHDWGCNAConfig:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return IntegrateHDWGCNAConfig(**raw)


def _normalize_gene_id(x: Any) -> str:
    s = str(x).strip()
    if s.upper().startswith("ENSG"):
        return s.split(".", 1)[0].upper()
    return s.upper()


def _safe_name(x: Any) -> str:
    s = str(x).strip()
    s = re.sub(r"[^\w.\-]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


def _infer_subclass_from_path(path: Path, root: Path) -> str:
    try:
        rel = path.relative_to(root)
        if len(rel.parts) >= 2:
            return rel.parts[0]
    except Exception:
        pass
    return path.parent.name


def _jaccard(a: list[int], b: list[int]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / max(len(sa | sb), 1)


# =====================================================================
# Registry JSON
# =====================================================================
def load_registry_json(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    required = [
        "gene_names",
        "module_names",
        "module_sources",
        "module_token_types",
        "membership_binary",
    ]
    missing = [k for k in required if k not in data]
    if missing:
        raise KeyError(f"Base registry missing required keys: {missing}")

    n_modules = len(data["module_names"])
    if not (
        len(data["module_sources"]) == n_modules
        and len(data["module_token_types"]) == n_modules
        and len(data["membership_binary"]) == n_modules
    ):
        raise ValueError("Base registry module arrays have inconsistent lengths.")

    n_genes = len(data["gene_names"])
    for i, row in enumerate(data["membership_binary"]):
        if len(row) != n_genes:
            raise ValueError(
                f"membership_binary row {i} has length {len(row)}; expected {n_genes}."
            )

    return data


def registry_json_to_modules(data: dict[str, Any]) -> list[dict[str, Any]]:
    modules: list[dict[str, Any]] = []
    for i, (name, source, token_type, row) in enumerate(
        zip(
            data["module_names"],
            data["module_sources"],
            data["module_token_types"],
            data["membership_binary"],
        )
    ):
        gene_indices = [j for j, v in enumerate(row) if float(v) > 0.0]
        modules.append(
            {
                "old_module_id": i,
                "name": str(name),
                "source": str(source),
                "token_type": str(token_type),
                "gene_indices": gene_indices,
            }
        )
    return modules


def write_registry_json(
    path: str | Path,
    *,
    gene_names: list[str],
    modules: list[dict[str, Any]],
    singleton_gene_indices: list[int],
    residual_gene_indices: list[int],
) -> None:
    n_genes = len(gene_names)

    module_names = []
    module_sources = []
    module_token_types = []
    membership_binary: list[list[float]] = []

    for mod in modules:
        row = [0.0] * n_genes
        for g in sorted(set(map(int, mod["gene_indices"]))):
            if g < 0 or g >= n_genes:
                raise ValueError(f"Gene index {g} out of range in module {mod['name']}.")
            row[g] = 1.0

        module_names.append(str(mod["name"]))
        module_sources.append(str(mod["source"]))
        module_token_types.append(str(mod["token_type"]))
        membership_binary.append(row)

    out = {
        "gene_names": gene_names,
        "module_names": module_names,
        "module_sources": module_sources,
        "module_token_types": module_token_types,
        "membership_binary": membership_binary,
        "singleton_gene_indices": sorted({int(x) for x in singleton_gene_indices}),
        "residual_gene_indices": sorted({int(x) for x in residual_gene_indices}),
    }

    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)


# =====================================================================
# Gene resolver
# =====================================================================
def load_symbol_to_gene_ids(path: Optional[str]) -> Optional[dict[str, list[str]]]:
    if path is None:
        return None
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return {str(k).upper(): [str(v) for v in vals] for k, vals in raw.items()}


def build_gene_resolver(
    gene_names: list[str],
    *,
    symbol_to_gene_ids: Optional[dict[str, list[str]]],
    allow_direct: bool,
):
    norm_to_idx = {_normalize_gene_id(g): i for i, g in enumerate(gene_names)}

    def resolve(ref: Any) -> list[int]:
        hits: list[int] = []

        key = _normalize_gene_id(ref)
        if allow_direct and key in norm_to_idx:
            hits.append(norm_to_idx[key])

        if symbol_to_gene_ids is not None:
            sym_key = str(ref).strip().upper()
            for gid in symbol_to_gene_ids.get(sym_key, []):
                gid_key = _normalize_gene_id(gid)
                if gid_key in norm_to_idx:
                    hits.append(norm_to_idx[gid_key])

        return sorted(set(hits))

    return resolve


# =====================================================================
# hdWGCNA CSV parsing
# =====================================================================
GENE_COL_CANDIDATES = [
    "gene_name",
    "gene",
    "genes",
    "feature",
    "feature_name",
    "gene_symbol",
    "symbol",
    "ensembl_id",
    "gene_id",
    "name",
]

MODULE_COL_CANDIDATES = [
    "module",
    "Module",
    "module_color",
    "color",
    "colour",
    "wgcna_module",
]


def discover_module_csvs(root: str | Path, cfg: IntegrateHDWGCNAConfig) -> list[Path]:
    root = Path(root).expanduser().resolve()

    if not root.exists():
        raise FileNotFoundError(f"hdwgcna_root_dir does not exist: {root}")

    # 1) Prefer files with kME_self
    paths = sorted(root.rglob("*_modules_with_kME_self.csv"))

    # 2) Fallback to ordinary module files
    if not paths:
        paths = sorted(root.rglob("*_modules.csv"))

    out: list[Path] = []
    for p in paths:
        lower = p.name.lower()

        # Keep real module assignment files only
        if "modules_with_kme_self" in lower or lower.endswith("_modules.csv"):
            pass
        else:
            continue

        # Exclude auxiliary files
        if any(bad in lower for bad in ["top20", "hub_genes", "power_table", "summary", "soft_powers"]):
            continue

        if any(part.lower() == "tom" for part in p.parts):
            continue

        out.append(p)

    if not out:
        example = list(root.rglob("*.csv"))[:10]
        raise FileNotFoundError(
            f"No hdWGCNA module CSVs found under {root}.\n"
            f"Example CSVs seen: {[str(x) for x in example]}"
        )

    print(f"[discover_module_csvs] found {len(out)} module CSV files")
    for p in out[:5]:
        print(f"  - {p}")

    return out
    

def choose_column(df: pd.DataFrame, explicit: Optional[str], candidates: list[str], kind: str) -> str:
    if explicit is not None:
        if explicit not in df.columns:
            raise KeyError(f"Explicit {kind}_column={explicit!r} not found. Columns={list(df.columns)}")
        return explicit

    exact = {str(c): c for c in df.columns}
    lower = {str(c).lower(): c for c in df.columns}

    for cand in candidates:
        if cand in exact:
            return exact[cand]
        if cand.lower() in lower:
            return lower[cand.lower()]

    if kind == "gene":
        # Common case: gene names were saved as an unnamed first index column.
        first = str(df.columns[0])
        if first.startswith("Unnamed"):
            return df.columns[0]

        for c in df.columns:
            cl = str(c).lower()
            if "gene" in cl or "symbol" in cl or "feature" in cl:
                return c

    raise KeyError(
        f"Could not infer {kind} column. Columns={list(df.columns)}. "
        f"Pass {kind}_column explicitly in config."
    )


def parse_hdwgcna_csv(
    path: Path,
    *,
    root: Path,
    cfg: IntegrateHDWGCNAConfig,
    resolve_gene,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    df = pd.read_csv(path)

    gene_col = choose_column(df, cfg.gene_column, GENE_COL_CANDIDATES, "gene")
    module_col = choose_column(df, cfg.module_column, MODULE_COL_CANDIDATES, "module")

    subclass = _safe_name(_infer_subclass_from_path(path, root))
    grey = {str(x).lower() for x in cfg.grey_names}

    modules: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []

    for module_value, gdf in df.groupby(module_col, observed=True, sort=True):
        module_label = str(module_value).strip()
        if cfg.drop_grey and module_label.lower() in grey:
            continue

        raw_gene_refs = gdf[gene_col].astype(str).tolist()
        mapped: list[int] = [] # set of index number in a module
        n_unmapped = 0

        for ref in raw_gene_refs:
            hits = resolve_gene(ref)
            if not hits:
                n_unmapped += 1
            mapped.extend(hits)

        idx = sorted(set(mapped))
        status = "kept"
        if len(idx) < cfg.min_module_size:
            status = "drop_small"
        elif cfg.max_module_size is not None and len(idx) > cfg.max_module_size:
            status = "drop_large"

        name = f"{cfg.hdwgcna_name_prefix}::{subclass}::{module_label}"
        source = f"{cfg.hdwgcna_source_prefix}::{subclass}"

        audit.append(
            {
                "csv_path": str(path),
                "subclass": subclass,
                "raw_module": module_label,
                "module_name": name,
                "source": source,
                "gene_column": str(gene_col),
                "module_column": str(module_col),
                "raw_unique_genes": int(len(set(raw_gene_refs))),
                "mapped_genes": int(len(idx)),
                "unmapped_gene_refs": int(n_unmapped),
                "status": status,
            }
        )

        if status != "kept":
            continue

        mod: dict[str, Any] = {
            "name": name,
            "source": source,
            "token_type": "residual",
            "gene_indices": idx,
            "subclass": subclass,
            "raw_module": module_label,
            "csv_path": str(path),
        }

        if cfg.kme_column in gdf.columns:
            kme_by_gene: dict[int, float] = {}
            for ref, val in zip(gdf[gene_col].astype(str).tolist(), gdf[cfg.kme_column].tolist()):
                if pd.isna(val):
                    continue
                for gi in resolve_gene(ref):
                    kme_by_gene[int(gi)] = float(val)
            if kme_by_gene:
                mod["kme_by_gene_index"] = kme_by_gene

        modules.append(mod)

    return modules, audit


def prune_redundant_modules(
    modules: list[dict[str, Any]],
    threshold: Optional[float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if threshold is None:
        return modules, []

    kept: list[dict[str, Any]] = []
    dropped: list[dict[str, Any]] = []

    for mod in sorted(modules, key=lambda m: (-len(m["gene_indices"]), m["name"])):
        redundant_with = None
        for prev in kept:
            if _jaccard(mod["gene_indices"], prev["gene_indices"]) >= float(threshold):
                redundant_with = prev["name"]
                break

        if redundant_with is None:
            kept.append(mod)
        else:
            d = dict(mod)
            d["status"] = "drop_redundant"
            d["redundant_with"] = redundant_with
            dropped.append(d)

    return kept, dropped


# =====================================================================
# Compact registry and weights
# =====================================================================
def make_compact_registry(
    *,
    gene_names_full: list[str],
    modules_full: list[dict[str, Any]],
    singleton_gene_indices_full: list[int],
    residual_gene_indices_full: list[int],
) -> tuple[list[str], list[dict[str, Any]], list[int], list[int], list[int]]:
    covered = sorted({int(g) for m in modules_full for g in m["gene_indices"]})
    old_to_new = {old: new for new, old in enumerate(covered)}

    gene_names = [gene_names_full[i] for i in covered]
    modules: list[dict[str, Any]] = []

    for mod in modules_full:
        new_idx = sorted({old_to_new[g] for g in mod["gene_indices"] if g in old_to_new})
        if not new_idx:
            continue

        new_mod = {
            "name": mod["name"],
            "source": mod["source"],
            "token_type": mod["token_type"],
            "gene_indices": new_idx,
        }

        if "kme_by_gene_index" in mod:
            new_mod["kme_by_gene_index"] = {
                old_to_new[int(k)]: float(v)
                for k, v in mod["kme_by_gene_index"].items()
                if int(k) in old_to_new
            }

        modules.append(new_mod)

    singleton_new = sorted({old_to_new[g] for g in singleton_gene_indices_full if g in old_to_new})
    residual_new = sorted({old_to_new[g] for g in residual_gene_indices_full if g in old_to_new})

    return gene_names, modules, singleton_new, residual_new, covered


def build_activity_weight(modules: list[dict[str, Any]], n_genes: int) -> np.ndarray:
    """
    Build [M, G] activity_weight.

    - hdWGCNA modules with kME_by_gene_index use signed kME values.
    - all other modules use positive L2 membership fallback.
    - every row is L2-normalised.
    """
    W = np.zeros((len(modules), n_genes), dtype=np.float32)

    for m, mod in enumerate(modules):
        kme = mod.get("kme_by_gene_index", {})
        if kme:
            for g, v in kme.items():
                gi = int(g)
                if 0 <= gi < n_genes:
                    W[m, gi] = float(v)

        if not np.any(W[m] != 0):
            for g in mod["gene_indices"]:
                gi = int(g)
                if 0 <= gi < n_genes:
                    W[m, gi] = 1.0

        norm = float(np.sqrt(np.sum(W[m] ** 2)))
        if norm > 0:
            W[m] /= norm

    return W


def strip_helper_fields(modules: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for m in modules:
        out.append(
            {
                "name": m["name"],
                "source": m["source"],
                "token_type": m["token_type"],
                "gene_indices": sorted(set(map(int, m["gene_indices"]))),
            }
        )
    return out


# =====================================================================
# Main
# =====================================================================
def integrate_hdwgcna_modules(cfg: IntegrateHDWGCNAConfig) -> dict[str, Any]:
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = load_registry_json(cfg.base_registry_json)
    gene_names_full = [str(x) for x in base["gene_names"]]
    base_modules = registry_json_to_modules(base)

    symbol_map = load_symbol_to_gene_ids(cfg.symbol_to_gene_ids_json)
    resolve_gene = build_gene_resolver(
        gene_names_full,
        symbol_to_gene_ids=symbol_map,
        allow_direct=cfg.allow_direct_gene_name_match,
    )

    kept_base: list[dict[str, Any]] = []
    dropped_base: list[dict[str, Any]] = []

    for mod in base_modules:
        source_lower = mod["source"].lower()
        is_existing_hdw = source_lower.startswith(cfg.hdwgcna_source_prefix.lower())
        drop = False

        if cfg.drop_existing_residual_modules and mod["token_type"] == "residual":
            drop = True
        if cfg.drop_existing_hdwgcna_modules and is_existing_hdw:
            drop = True

        if drop:
            dropped_base.append(mod)
        else:
            kept_base.append(
                {
                    "name": mod["name"],
                    "source": mod["source"],
                    "token_type": mod["token_type"],
                    "gene_indices": sorted(set(map(int, mod["gene_indices"]))),
                }
            )

    _log(
        f"base modules={len(base_modules):,}; kept={len(kept_base):,}; "
        f"dropped={len(dropped_base):,}"
    )

    root = Path(cfg.hdwgcna_root_dir)
    csvs = discover_module_csvs(root, cfg)
    _log(f"hdWGCNA CSV files discovered: {len(csvs):,}")

    imported_raw: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []

    for p in csvs:
        mods, rows = parse_hdwgcna_csv(p, root=root, cfg=cfg, resolve_gene=resolve_gene)
        imported_raw.extend(mods)
        audit.extend(rows)

    imported, redundant = prune_redundant_modules(
        imported_raw,
        cfg.prune_hdwgcna_redundancy_jaccard,
    )

    for d in redundant:
        audit.append(
            {
                "csv_path": d.get("csv_path", ""),
                "subclass": d.get("subclass", ""),
                "raw_module": d.get("raw_module", ""),
                "module_name": d["name"],
                "source": d["source"],
                "gene_column": "",
                "module_column": "",
                "raw_unique_genes": len(d["gene_indices"]),
                "mapped_genes": len(d["gene_indices"]),
                "unmapped_gene_refs": 0,
                "status": "drop_redundant",
                "redundant_with": d.get("redundant_with", ""),
            }
        )

    if cfg.fail_if_no_hdwgcna_modules and not imported:
        raise RuntimeError(
            "No hdWGCNA modules survived filtering. Check module CSV path, "
            "gene mapping, module column, gene column, and min/max size."
        )

    _log(
        f"hdWGCNA modules imported={len(imported):,} "
        f"(raw kept before redundancy={len(imported_raw):,}, redundant dropped={len(redundant):,})"
    )

    modules_full_with_helpers: list[dict[str, Any]] = kept_base + imported

    old_residual_full = [int(x) for x in base.get("residual_gene_indices", [])]
    singleton_full = [int(x) for x in base.get("singleton_gene_indices", [])]

    if cfg.add_leftover_residual_module and old_residual_full:
        covered = {g for m in modules_full_with_helpers for g in m["gene_indices"]}
        leftover = sorted(set(old_residual_full) - covered)
        if leftover:
            modules_full_with_helpers.append(
                {
                    "name": cfg.leftover_module_name,
                    "source": "derived::residual_leftover_after_hdwgcna",
                    "token_type": "residual",
                    "gene_indices": leftover,
                }
            )
            _log(f"leftover residual module added: {len(leftover):,} genes")
        else:
            _log("no leftover residual genes remain")

    residual_full = sorted({
        int(g)
        for m in modules_full_with_helpers
        if m["token_type"] == "residual"
        for g in m["gene_indices"]
    })

    modules_full_json = strip_helper_fields(modules_full_with_helpers)

    full_registry_path = out_dir / "final_gene_module_registry_full.json"
    if cfg.write_full_registry:
        _log(f"writing full registry: {full_registry_path}")
        write_registry_json(
            full_registry_path,
            gene_names=gene_names_full,
            modules=modules_full_json,
            singleton_gene_indices=singleton_full,
            residual_gene_indices=residual_full,
        )

    compact_registry_path = out_dir / "final_gene_module_registry_compact.json"
    compact_gene_names: list[str] = []
    compact_modules_with_helpers: list[dict[str, Any]] = []
    compact_modules_json: list[dict[str, Any]] = []
    compact_singletons: list[int] = []
    compact_residuals: list[int] = []
    original_indices: list[int] = []

    if cfg.write_compact_registry:
        (
            compact_gene_names,
            compact_modules_with_helpers,
            compact_singletons,
            compact_residuals,
            original_indices,
        ) = make_compact_registry(
            gene_names_full=gene_names_full,
            modules_full=modules_full_with_helpers,
            singleton_gene_indices_full=singleton_full,
            residual_gene_indices_full=residual_full,
        )

        compact_modules_json = strip_helper_fields(compact_modules_with_helpers)

        _log(
            f"writing compact registry: {compact_registry_path} "
            f"({len(compact_gene_names):,} genes, {len(compact_modules_json):,} modules)"
        )
        write_registry_json(
            compact_registry_path,
            gene_names=compact_gene_names,
            modules=compact_modules_json,
            singleton_gene_indices=compact_singletons,
            residual_gene_indices=compact_residuals,
        )

        if cfg.write_gene_union:
            with open(out_dir / "gene_union_ensembl.json", "w", encoding="utf-8") as f:
                json.dump(compact_gene_names, f, ensure_ascii=False, indent=2)
            np.save(out_dir / "gene_union_original_indices.npy", np.asarray(original_indices, dtype=np.int64))

    if cfg.write_activity_weight:
        if cfg.write_compact_registry:
            W = build_activity_weight(compact_modules_with_helpers, len(compact_gene_names))
            weight_gene_names = compact_gene_names
            weight_module_names = [m["name"] for m in compact_modules_with_helpers]
            weight_registry = str(compact_registry_path)
        else:
            W = build_activity_weight(modules_full_with_helpers, len(gene_names_full))
            weight_gene_names = gene_names_full
            weight_module_names = [m["name"] for m in modules_full_with_helpers]
            weight_registry = str(full_registry_path)

        np.savez_compressed(
            out_dir / "activity_weight_kme_or_l2_membership.npz",
            activity_weight=W,
            gene_names=np.asarray(weight_gene_names, dtype=object),
            module_names=np.asarray(weight_module_names, dtype=object),
            registry_json_path=np.asarray(weight_registry),
            note=np.asarray(
                "hdWGCNA residual modules use signed kME_self when available; "
                "all other modules use positive L2 membership fallback."
            ),
        )

    if cfg.write_compact_registry:
        module_table = pd.DataFrame(
            {
                "module_name": [m["name"] for m in compact_modules_json],
                "source": [m["source"] for m in compact_modules_json],
                "token_type": [m["token_type"] for m in compact_modules_json],
                "size": [len(m["gene_indices"]) for m in compact_modules_json],
            }
        )
    else:
        module_table = pd.DataFrame(
            {
                "module_name": [m["name"] for m in modules_full_json],
                "source": [m["source"] for m in modules_full_json],
                "token_type": [m["token_type"] for m in modules_full_json],
                "size": [len(m["gene_indices"]) for m in modules_full_json],
            }
        )

    module_table.to_csv(out_dir / "final_module_table.csv", index=False)
    pd.DataFrame(audit).to_csv(out_dir / "hdwgcna_module_import_audit.csv", index=False)

    summary = {
        "config": asdict(cfg),
        "base_registry_json": str(cfg.base_registry_json),
        "n_base_genes": int(len(gene_names_full)),
        "n_base_modules": int(len(base_modules)),
        "n_base_modules_kept": int(len(kept_base)),
        "n_base_modules_dropped": int(len(dropped_base)),
        "n_hdwgcna_csvs": int(len(csvs)),
        "n_hdwgcna_modules_imported": int(len(imported)),
        "n_final_full_modules": int(len(modules_full_json)),
        "n_final_full_genes": int(len(gene_names_full)),
        "n_final_compact_modules": int(len(compact_modules_json)) if cfg.write_compact_registry else None,
        "n_final_compact_genes": int(len(compact_gene_names)) if cfg.write_compact_registry else None,
        "full_registry_path": str(full_registry_path) if cfg.write_full_registry else None,
        "compact_registry_path": str(compact_registry_path) if cfg.write_compact_registry else None,
        "gene_union_ensembl_path": str(out_dir / "gene_union_ensembl.json")
        if cfg.write_gene_union and cfg.write_compact_registry else None,
        "gene_union_original_indices_path": str(out_dir / "gene_union_original_indices.npy")
        if cfg.write_gene_union and cfg.write_compact_registry else None,
        "activity_weight_path": str(out_dir / "activity_weight_kme_or_l2_membership.npz")
        if cfg.write_activity_weight else None,
    }

    with open(out_dir / "hdwgcna_integration_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    _log("done")
    return summary


# =====================================================================
# CLI
# =====================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str)
    args = parser.parse_args()

    cfg = load_config(args.config)
    integrate_hdwgcna_modules(cfg)


if __name__ == "__main__":
    main()