from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Sequence

import torch


TokenType = Literal["regular", "residual", "singleton"]
MappingStrategy = Literal["all", "first", "unambiguous"]


# =====================================================================
# Dataclasses
# =====================================================================
@dataclass(frozen=True)
class ModuleCollectionSpec:
    """
    Specification for one external gene-set collection.

    Parameters
    ----------
    path:
        Local path to a GMT file.
    source:
        Human-readable source label such as "reactome", "go_bp", "regulon",
        "kegg", or "hallmark".
    token_type:
        Usually "regular". Residual and singleton tokens are constructed by the
        registry itself and should not usually be supplied as collections.
    min_size / max_size:
        Filter on the number of genes that survive intersection with the current
        model vocabulary.
    prefix:
        Optional prefix added to module names, e.g. "REACTOME::".
    max_modules:
        Optional cap after filtering and pruning. When used, larger retained
        modules are prioritized within the already size-filtered range.
    keep_names:
        Optional allow-list of set names.
    drop_names:
        Optional deny-list of set names.
    """

    path: str
    source: str
    token_type: TokenType = "regular"
    min_size: int = 10
    max_size: int = 300
    prefix: Optional[str] = None
    max_modules: Optional[int] = None
    keep_names: Optional[Sequence[str]] = None
    drop_names: Optional[Sequence[str]] = None


@dataclass
class ModuleRecord:
    module_id: int
    name: str
    source: str
    token_type: TokenType
    gene_indices: List[int]

    @property
    def size(self) -> int:
        return len(self.gene_indices)


@dataclass
class RegistryState:
    gene_names: List[str]
    module_names: List[str]
    module_sources: List[str]
    module_token_types: List[str]
    membership_binary: List[List[float]]
    singleton_gene_indices: List[int]
    residual_gene_indices: List[int]


# =====================================================================
# Helpers
# =====================================================================
def _as_upper(x: str) -> str:
    return str(x).strip().upper()


def _normalize_id(x: str) -> str:
    """Normalize identifiers like ENSG000001234.5 -> ENSG000001234."""
    s = str(x).strip()
    if s.upper().startswith("ENSG"):
        return s.split(".", 1)[0].upper()
    return _as_upper(s)


def _read_gmt(path: str) -> List[tuple[str, str, List[str]]]:
    out: List[tuple[str, str, List[str]]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            desc = parts[1]
            genes = [g for g in parts[2:] if g]
            out.append((name, desc, genes))
    return out


def _jaccard(a: Sequence[int], b: Sequence[int]) -> float:
    sa = set(a)
    sb = set(b)
    if not sa and not sb:
        return 1.0
    inter = len(sa & sb)
    union = len(sa | sb)
    return inter / max(union, 1)


def build_symbol_to_gene_ids_map(
    gene_ids: Sequence[str],
    gene_symbols: Sequence[str],
    *,
    strategy: MappingStrategy = "unambiguous",
) -> Dict[str, List[str]]:
    """
    Build symbol -> internal-gene-id map.

    Parameters
    ----------
    gene_ids:
        Usually ENSG ids in the model vocabulary.
    gene_symbols:
        Display symbols aligned to gene_ids.
    strategy:
        - "all"          : keep all mappings
        - "first"        : keep first mapping only
        - "unambiguous"  : keep only symbols mapping to exactly one gene id
    """
    if len(gene_ids) != len(gene_symbols):
        raise ValueError("gene_ids and gene_symbols must have the same length.")

    grouped: Dict[str, List[str]] = {}
    for gid, sym in zip(gene_ids, gene_symbols):
        sym = _as_upper(sym)
        gid = _normalize_id(gid)
        if sym == "":
            continue
        grouped.setdefault(sym, []).append(gid)

    if strategy == "all":
        return grouped
    if strategy == "first":
        return {k: [v[0]] for k, v in grouped.items() if len(v) > 0}
    if strategy == "unambiguous":
        return {k: v for k, v in grouped.items() if len(set(v)) == 1}
    raise ValueError(f"Unknown strategy '{strategy}'.")


# =====================================================================
# Main registry
# =====================================================================
class GeneModuleRegistry:
    """
    External gene-set registry for pathway / regulon / WGCNA-style modules.

    Philosophy
    ----------
    1. The registry is static infrastructure.
    2. It stores external gene-module definitions.
    3. It supports overlapping memberships, which is crucial for pathway and
       regulon databases.
    4. It can emit both:
         - weighted membership matrices for encoder tokenization
         - binary masks for decoder-side Lie generator restriction

    Typical usage
    -------------
    - Internal keys: ENSG (recommended)
    - External GMTs: symbols
    - Bridge: symbol_to_gene_ids mapping
    """

    TOKEN_TYPE_TO_ID = {"regular": 0, "residual": 1, "singleton": 2}

    def __init__(
        self,
        gene_names: Sequence[str],
        modules: List[ModuleRecord],
        membership_binary: torch.Tensor,
        singleton_gene_indices: Sequence[int],
        residual_gene_indices: Sequence[int],
    ) -> None:
        if membership_binary.ndim != 2:
            raise ValueError("membership_binary must have shape [M, G].")
        n_modules, n_genes = membership_binary.shape
        if len(gene_names) != n_genes:
            raise ValueError("gene_names length must match membership matrix width.")
        if len(modules) != n_modules:
            raise ValueError("modules length must match membership matrix height.")

        self.gene_names = list(gene_names)
        self.n_genes = n_genes
        self.modules = modules
        self.n_modules = n_modules
        self.gene_to_idx = {_normalize_id(g): i for i, g in enumerate(self.gene_names)}

        self.membership_binary = membership_binary.to(dtype=torch.float32)
        row_sums = self.membership_binary.sum(dim=1, keepdim=True).clamp_min(1.0)
        self.membership_weight = self.membership_binary / row_sums

        self.singleton_gene_indices = sorted({int(x) for x in singleton_gene_indices})
        self.residual_gene_indices = sorted({int(x) for x in residual_gene_indices})

        self.module_names = [m.name for m in self.modules]
        self.module_sources = [m.source for m in self.modules]
        self.module_token_types = [m.token_type for m in self.modules]
        self.module_type_ids = torch.tensor(
            [self.TOKEN_TYPE_TO_ID[m.token_type] for m in self.modules],
            dtype=torch.long,
        )

    # ==================================================================
    # Constructors
    # ==================================================================
    @classmethod
    def from_collection_specs(
        cls,
        gene_names: Sequence[str],
        collection_specs: Sequence[ModuleCollectionSpec],
        *,
        singleton_genes: Optional[Sequence[str | int]] = None,
        add_residual_module: bool = True,
        residual_name: str = "RESIDUAL_LEFTOVER",
        prune_redundancy_jaccard: Optional[float] = 0.95,
        source_priority: Optional[Dict[str, int]] = None,
        symbol_to_gene_ids: Optional[Dict[str, Sequence[str]]] = None,
        allow_direct_gene_name_match: bool = True,
    ) -> "GeneModuleRegistry":
        """
        Build a registry from multiple local GMT collections.

        Important
        ---------
        If gene_names are ENSG ids and GMT files use symbols, pass
        symbol_to_gene_ids.
        """
        if len(gene_names) == 0:
            raise ValueError("gene_names must not be empty.")

        gene_names = list(gene_names)
        gene_to_idx = {_normalize_id(g): i for i, g in enumerate(gene_names)}

        singleton_idx = cls._resolve_gene_refs(
            refs=singleton_genes or [],
            gene_to_idx=gene_to_idx,
            n_genes=len(gene_names),
            symbol_to_gene_ids=symbol_to_gene_ids,
            allow_direct_gene_name_match=allow_direct_gene_name_match,
        )

        if source_priority is None:
            source_priority = {}

        candidates: List[ModuleRecord] = []
        for spec in collection_specs:
            raw_sets = _read_gmt(spec.path)
            keep_names = None if spec.keep_names is None else {_as_upper(x) for x in spec.keep_names}
            drop_names = None if spec.drop_names is None else {_as_upper(x) for x in spec.drop_names}

            local_records: List[ModuleRecord] = []
            for name, _desc, genes in raw_sets:
                up_name = _as_upper(name)
                if keep_names is not None and up_name not in keep_names:
                    continue
                if drop_names is not None and up_name in drop_names:
                    continue

                idx = cls._resolve_gene_refs(
                    refs=genes,
                    gene_to_idx=gene_to_idx,
                    n_genes=len(gene_names),
                    symbol_to_gene_ids=symbol_to_gene_ids,
                    allow_direct_gene_name_match=allow_direct_gene_name_match,
                )
                if len(idx) < spec.min_size or len(idx) > spec.max_size:
                    continue

                full_name = f"{spec.prefix}{name}" if spec.prefix else name
                local_records.append(
                    ModuleRecord(
                        module_id=-1,
                        name=full_name,
                        source=spec.source,
                        token_type=spec.token_type,
                        gene_indices=idx,
                    )
                )

            if spec.max_modules is not None and len(local_records) > spec.max_modules:
                # Keep the largest modules within the already size-filtered band.
                local_records = sorted(local_records, key=lambda m: (-m.size, m.name))[: spec.max_modules]

            candidates.extend(local_records)

        # Greedy redundancy pruning.
        if prune_redundancy_jaccard is not None:
            kept: List[ModuleRecord] = []

            def _rank_key(rec: ModuleRecord) -> tuple:
                pr = source_priority.get(rec.source, 999)
                return (pr, -rec.size, rec.name)

            for rec in sorted(candidates, key=_rank_key):
                redundant = False
                for prev in kept:
                    if _jaccard(rec.gene_indices, prev.gene_indices) >= prune_redundancy_jaccard:
                        redundant = True
                        break
                if not redundant:
                    kept.append(rec)
            candidates = kept

        # Reassign module ids after pruning.
        modules: List[ModuleRecord] = []
        next_mid = 0
        for rec in candidates:
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=rec.name,
                    source=rec.source,
                    token_type=rec.token_type,
                    gene_indices=list(rec.gene_indices),
                )
            )
            next_mid += 1

        # Residual pseudo-module.
        covered = torch.zeros(len(gene_names), dtype=torch.bool)
        for mod in modules:
            covered[torch.tensor(mod.gene_indices, dtype=torch.long)] = True

        residual_gene_indices = (~covered).nonzero(as_tuple=True)[0].tolist()
        if add_residual_module and len(residual_gene_indices) > 0:
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=residual_name,
                    source="derived::residual",
                    token_type="residual",
                    gene_indices=list(residual_gene_indices),
                )
            )
            next_mid += 1

        # Singleton modules are additional tokens, not replacements.
        for g in sorted(singleton_idx):
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=f"SINGLETON::{gene_names[g]}",
                    source="derived::singleton",
                    token_type="singleton",
                    gene_indices=[g],
                )
            )
            next_mid += 1

        membership_binary = torch.zeros(len(modules), len(gene_names), dtype=torch.float32)
        for mod in modules:
            membership_binary[mod.module_id, mod.gene_indices] = 1.0

        return cls(
            gene_names=gene_names,
            modules=modules,
            membership_binary=membership_binary,
            singleton_gene_indices=singleton_idx,
            residual_gene_indices=residual_gene_indices,
        )

    @classmethod
    def from_module_sets(
        cls,
        gene_names: Sequence[str],
        module_sets: Dict[str, Sequence[str | int]],
        *,
        source: str = "custom",
        singleton_genes: Optional[Sequence[str | int]] = None,
        add_residual_module: bool = True,
        residual_name: str = "RESIDUAL_LEFTOVER",
        min_size: int = 10,
        max_size: int = 300,
        prune_redundancy_jaccard: Optional[float] = None,
        symbol_to_gene_ids: Optional[Dict[str, Sequence[str]]] = None,
        allow_direct_gene_name_match: bool = True,
    ) -> "GeneModuleRegistry":
        gene_names = list(gene_names)
        gene_to_idx = {_normalize_id(g): i for i, g in enumerate(gene_names)}
        singleton_idx = cls._resolve_gene_refs(
            refs=singleton_genes or [],
            gene_to_idx=gene_to_idx,
            n_genes=len(gene_names),
            symbol_to_gene_ids=symbol_to_gene_ids,
            allow_direct_gene_name_match=allow_direct_gene_name_match,
        )

        modules: List[ModuleRecord] = []
        next_mid = 0
        for name, refs in module_sets.items():
            idx = sorted(
                cls._resolve_gene_refs(
                    refs,
                    gene_to_idx,
                    len(gene_names),
                    symbol_to_gene_ids=symbol_to_gene_ids,
                    allow_direct_gene_name_match=allow_direct_gene_name_match,
                )
            )
            if len(idx) < min_size or len(idx) > max_size:
                continue
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=name,
                    source=source,
                    token_type="regular",
                    gene_indices=idx,
                )
            )
            next_mid += 1

        if prune_redundancy_jaccard is not None:
            kept: List[ModuleRecord] = []
            for rec in sorted(modules, key=lambda m: (-m.size, m.name)):
                redundant = False
                for prev in kept:
                    if _jaccard(rec.gene_indices, prev.gene_indices) >= prune_redundancy_jaccard:
                        redundant = True
                        break
                if not redundant:
                    kept.append(rec)
            modules = [
                ModuleRecord(
                    module_id=i,
                    name=m.name,
                    source=m.source,
                    token_type=m.token_type,
                    gene_indices=m.gene_indices,
                )
                for i, m in enumerate(kept)
            ]
            next_mid = len(modules)

        covered = torch.zeros(len(gene_names), dtype=torch.bool)
        for mod in modules:
            covered[torch.tensor(mod.gene_indices, dtype=torch.long)] = True
        residual_gene_indices = (~covered).nonzero(as_tuple=True)[0].tolist()

        if add_residual_module and len(residual_gene_indices) > 0:
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=residual_name,
                    source="derived::residual",
                    token_type="residual",
                    gene_indices=residual_gene_indices,
                )
            )
            next_mid += 1

        for g in sorted(singleton_idx):
            modules.append(
                ModuleRecord(
                    module_id=next_mid,
                    name=f"SINGLETON::{gene_names[g]}",
                    source="derived::singleton",
                    token_type="singleton",
                    gene_indices=[g],
                )
            )
            next_mid += 1

        membership_binary = torch.zeros(len(modules), len(gene_names), dtype=torch.float32)
        for mod in modules:
            membership_binary[mod.module_id, mod.gene_indices] = 1.0

        return cls(
            gene_names=gene_names,
            modules=modules,
            membership_binary=membership_binary,
            singleton_gene_indices=singleton_idx,
            residual_gene_indices=residual_gene_indices,
        )

    @staticmethod
    def _resolve_gene_refs(
        refs: Sequence[str | int],
        gene_to_idx: Dict[str, int],
        n_genes: int,
        *,
        symbol_to_gene_ids: Optional[Dict[str, Sequence[str]]] = None,
        allow_direct_gene_name_match: bool = True,
    ) -> List[int]:
        out: List[int] = []
        for x in refs:
            if isinstance(x, int):
                if x < 0 or x >= n_genes:
                    raise ValueError(f"Gene index {x} out of range [0, {n_genes}).")
                out.append(int(x))
                continue

            key = _normalize_id(x)

            # 1) direct match to internal ids/names
            if allow_direct_gene_name_match and key in gene_to_idx:
                out.append(gene_to_idx[key])
                continue

            # 2) symbol -> internal ids mapping
            if symbol_to_gene_ids is not None:
                symbol_key = _as_upper(x)
                hits = symbol_to_gene_ids.get(symbol_key, [])
                for gid in hits:
                    gid_norm = _normalize_id(gid)
                    if gid_norm in gene_to_idx:
                        out.append(gene_to_idx[gid_norm])

        return sorted(set(out))

    # ==================================================================
    # Queries / summaries
    # ==================================================================
    @property
    def regular_module_ids(self) -> List[int]:
        return [m.module_id for m in self.modules if m.token_type == "regular"]

    @property
    def residual_module_ids(self) -> List[int]:
        return [m.module_id for m in self.modules if m.token_type == "residual"]

    @property
    def singleton_module_ids(self) -> List[int]:
        return [m.module_id for m in self.modules if m.token_type == "singleton"]

    def regular_module_mask(self) -> torch.Tensor:
        return torch.tensor([m.token_type == "regular" for m in self.modules], dtype=torch.bool)

    def residual_module_mask(self) -> torch.Tensor:
        return torch.tensor([m.token_type == "residual" for m in self.modules], dtype=torch.bool)

    def singleton_module_mask(self) -> torch.Tensor:
        return torch.tensor([m.token_type == "singleton" for m in self.modules], dtype=torch.bool)

    def module_sizes(self) -> torch.Tensor:
        return self.membership_binary.sum(dim=1).long()

    def gene_cover_counts(self) -> torch.Tensor:
        return self.membership_binary.sum(dim=0).long()

    def get_module_genes(self, module_id: int) -> List[int]:
        if module_id < 0 or module_id >= self.n_modules:
            raise ValueError(f"module_id {module_id} out of range [0, {self.n_modules}).")
        return list(self.modules[module_id].gene_indices)

    def summary(self, max_lines: int = 20) -> str:
        lines = [
            "GeneModuleRegistry",
            f"  n_genes          : {self.n_genes}",
            f"  n_modules        : {self.n_modules}",
            f"    regular        : {len(self.regular_module_ids)}",
            f"    residual       : {len(self.residual_module_ids)}",
            f"    singleton      : {len(self.singleton_module_ids)}",
            f"  uncovered genes  : {(self.gene_cover_counts() == 0).sum().item()}",
            f"  multi-covered    : {(self.gene_cover_counts() > 1).sum().item()}",
        ]
        lines.append("  modules:")
        for mod in self.modules[:max_lines]:
            lines.append(
                f"    [{mod.module_id:04d}] {mod.name:40s}  size={mod.size:4d}  "
                f"type={mod.token_type:9s}  source={mod.source}"
            )
        if len(self.modules) > max_lines:
            lines.append(f"    ... ({len(self.modules) - max_lines} more)")
        return "\n".join(lines)

    # ==================================================================
    # Decoder mask generation
    # ==================================================================
    def make_decoder_module_masks(
        self,
        n_generators: int,
        *,
        strategy: Literal["one_to_one", "proportional", "cyclic", "full"] = "proportional",
        restrict_to_regular_modules: bool = False,
        device: torch.device = torch.device("cpu"),
    ) -> torch.Tensor:
        if n_generators <= 0:
            raise ValueError("n_generators must be positive.")

        if strategy == "full":
            return torch.ones(n_generators, self.n_genes, device=device)

        if restrict_to_regular_modules:
            base_ids = self.regular_module_ids
        else:
            base_ids = list(range(self.n_modules))
        if len(base_ids) == 0:
            return torch.ones(n_generators, self.n_genes, device=device)

        masks = torch.zeros(n_generators, self.n_genes, device=device)

        if strategy == "one_to_one":
            ranked = sorted(base_ids, key=lambda mid: self.modules[mid].size, reverse=True)
            for gen_idx in range(n_generators):
                if gen_idx < len(ranked):
                    masks[gen_idx] = self.membership_binary[ranked[gen_idx]].to(device)
                else:
                    masks[gen_idx] = 1.0
            return masks

        if strategy == "cyclic":
            ranked = sorted(base_ids, key=lambda mid: self.modules[mid].size, reverse=True)
            for gen_idx in range(n_generators):
                mid = ranked[gen_idx % len(ranked)]
                masks[gen_idx] = self.membership_binary[mid].to(device)
            return masks

        if strategy == "proportional":
            sizes = torch.tensor([self.modules[mid].size for mid in base_ids], dtype=torch.float32)
            probs = sizes / sizes.sum().clamp_min(1.0)
            raw = torch.floor(probs * n_generators).to(torch.long)
            alloc = raw.tolist()
            leftover = n_generators - sum(alloc)
            rank = torch.argsort(sizes, descending=True).tolist()
            for i in range(leftover):
                alloc[rank[i % len(rank)]] += 1

            gen_idx = 0
            for local_idx, mid in enumerate(base_ids):
                for _ in range(alloc[local_idx]):
                    if gen_idx >= n_generators:
                        break
                    masks[gen_idx] = self.membership_binary[mid].to(device)
                    gen_idx += 1
            while gen_idx < n_generators:
                masks[gen_idx] = 1.0
                gen_idx += 1
            return masks

        raise ValueError(f"Unknown strategy '{strategy}'.")

    # ==================================================================
    # Sparse selection helpers
    # ==================================================================
    def select_genes_from_modules(
        self,
        active_modules: torch.Tensor,
        *,
        include_residual: bool = True,
        include_singletons: bool = True,
    ) -> torch.Tensor:
        if active_modules.shape != (self.n_modules,):
            raise ValueError(
                f"active_modules must have shape ({self.n_modules},), got {tuple(active_modules.shape)}."
            )
        gene_mask = torch.zeros(self.n_genes, dtype=torch.bool, device=active_modules.device)
        chosen = active_modules.clone()
        if include_residual and len(self.residual_module_ids) > 0:
            chosen[self.residual_module_ids] = True
        if include_singletons and len(self.singleton_module_ids) > 0:
            chosen[self.singleton_module_ids] = True
        chosen_idx = chosen.nonzero(as_tuple=True)[0]
        if len(chosen_idx) > 0:
            gene_mask = self.membership_binary[chosen_idx].sum(dim=0).to(device=active_modules.device) > 0
        return gene_mask

    def exploration_gene_mask(
        self,
        active_modules: torch.Tensor,
        *,
        explore_fraction: float = 0.05,
        include_residual: bool = True,
        include_singletons: bool = True,
        generator: Optional[torch.Generator] = None,
    ) -> torch.Tensor:
        if not (0.0 <= explore_fraction <= 1.0):
            raise ValueError("explore_fraction must be in [0, 1].")
        base = self.select_genes_from_modules(
            active_modules,
            include_residual=include_residual,
            include_singletons=include_singletons,
        )
        if explore_fraction == 0.0:
            return base
        inactive = (~base).nonzero(as_tuple=True)[0]
        n_inactive = inactive.shape[0]
        if n_inactive == 0:
            return base
        n_pick = max(1, int(n_inactive * explore_fraction))
        perm = torch.randperm(n_inactive, generator=generator, device=base.device)
        picked = inactive[perm[:n_pick]]
        out = base.clone()
        out[picked] = True
        return out

    def random_gene_subsample_mask(
        self,
        n_genes_per_step: int,
        *,
        ensure_singletons: bool = True,
        ensure_residual: bool = True,
        generator: Optional[torch.Generator] = None,
        device: torch.device = torch.device("cpu"),
    ) -> torch.Tensor:
        if n_genes_per_step <= 0:
            raise ValueError("n_genes_per_step must be positive.")
        if n_genes_per_step >= self.n_genes:
            return torch.ones(self.n_genes, dtype=torch.bool, device=device)

        mask = torch.zeros(self.n_genes, dtype=torch.bool, device=device)
        forced: List[int] = []
        if ensure_residual:
            forced.extend(self.residual_gene_indices[:1])
        if ensure_singletons:
            forced.extend(self.singleton_gene_indices)
        forced = sorted(set([g for g in forced if 0 <= g < self.n_genes]))
        if len(forced) > 0:
            mask[torch.tensor(forced, dtype=torch.long, device=device)] = True

        budget = max(0, n_genes_per_step - int(mask.sum().item()))
        if budget > 0:
            cand = (~mask).nonzero(as_tuple=True)[0]
            n_cand = cand.shape[0]
            n_pick = min(budget, n_cand)
            perm = torch.randperm(n_cand, generator=generator, device=device)
            mask[cand[perm[:n_pick]]] = True
        return mask

    # ==================================================================
    # Serialization
    # ==================================================================
    def state(self) -> RegistryState:
        return RegistryState(
            gene_names=self.gene_names,
            module_names=[m.name for m in self.modules],
            module_sources=[m.source for m in self.modules],
            module_token_types=[m.token_type for m in self.modules],
            membership_binary=self.membership_binary.tolist(),
            singleton_gene_indices=self.singleton_gene_indices,
            residual_gene_indices=self.residual_gene_indices,
        )

    def save(self, path: str | Path) -> None:
        st = self.state()
        data = {
            "gene_names": st.gene_names,
            "module_names": st.module_names,
            "module_sources": st.module_sources,
            "module_token_types": st.module_token_types,
            "membership_binary": st.membership_binary,
            "singleton_gene_indices": st.singleton_gene_indices,
            "residual_gene_indices": st.residual_gene_indices,
        }
        with open(Path(path), "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: str | Path) -> "GeneModuleRegistry":
        with open(Path(path), "r", encoding="utf-8") as f:
            data = json.load(f)

        membership_binary = torch.tensor(data["membership_binary"], dtype=torch.float32)
        modules: List[ModuleRecord] = []
        for mid, (name, source, token_type) in enumerate(
            zip(data["module_names"], data["module_sources"], data["module_token_types"])
        ):
            idx = (membership_binary[mid] > 0).nonzero(as_tuple=True)[0].tolist()
            modules.append(
                ModuleRecord(
                    module_id=mid,
                    name=name,
                    source=source,
                    token_type=token_type,
                    gene_indices=idx,
                )
            )

        return cls(
            gene_names=data["gene_names"],
            modules=modules,
            membership_binary=membership_binary,
            singleton_gene_indices=data.get("singleton_gene_indices", []),
            residual_gene_indices=data.get("residual_gene_indices", []),
        )

    def __repr__(self) -> str:
        return (
            f"GeneModuleRegistry(n_genes={self.n_genes}, n_modules={self.n_modules}, "
            f"regular={len(self.regular_module_ids)}, residual={len(self.residual_module_ids)}, "
            f"singletons={len(self.singleton_module_ids)})"
        )


if __name__ == "__main__":
    gene_names = [f"ENSG{i:06d}" for i in range(20)]
    gene_symbols = [f"G{i}" for i in range(20)]
    s2g = build_symbol_to_gene_ids_map(gene_names, gene_symbols, strategy="unambiguous")

    module_sets = {
        "reactome_synapse": ["G0", "G1", "G2", "G3", "G4", "G5"],
        "go_stress": ["G4", "G5", "G6", "G7", "G8", "G9"],
        "regulon_tf1": ["G8", "G9", "G10", "G11", "G12", "G13"],
    }
    reg = GeneModuleRegistry.from_module_sets(
        gene_names=gene_names,
        module_sets=module_sets,
        singleton_genes=["G1", "G17"],
        min_size=2,
        max_size=20,
        add_residual_module=True,
        symbol_to_gene_ids=s2g,
    )
    print(reg.summary())
    masks = reg.make_decoder_module_masks(6, strategy="proportional")
    print("decoder masks:", tuple(masks.shape))