"""
Latent-kNN π_tech: a disease-aware estimate of P(an observed zero is technical).

Problem
-------
In scRNA-seq a measured ``0`` for gene ``g`` in cell ``i`` is either
  * a *biological* zero — the gene is genuinely off in that cell's state, or
  * a *technical* zero — the gene was expressed but the measurement missed it
    (low capture / shallow sequencing).
The L1 zero/nonzero loss (`hierarchical_ordinal.compute_zero_nonzero_loss`)
currently treats *every* observed zero as a hard biological-off label, which
punishes the model for suspecting an expressed-but-dropped gene and caps
nonzero learning.

We cannot label a single zero with certainty, so we estimate a probability
``π_tech[i,g] = P(this observed zero is technical)`` and use it to *relax* the
L1 penalty on suspicious zeros (soft target + down-weight). Two signals:

  1. Measurement quality (biology-agnostic, safe everywhere): a zero is more
     likely technical when the *cell* was shallowly sequenced AND the *gene* is
     generally detectable. Neither part assumes the gene "should" be on, so it
     is valid in disease cells too.
  2. Neighborhood (disease-aware): do cells biologically *similar* to ``i`` —
     its nearest neighbors in the latent state ``z_perp`` — detect gene ``g``?
     Neighbors are taken in ``z_perp`` (the state-deviation latent), NOT by
     clinical pathology labels, so "similar" means "similar disease state".
       * neighbors detect it, cell i is 0, depth low  → likely technical.
       * neighbors are ALSO 0                         → genuine off (e.g. a real
         disease down-regulation) → π_tech forced low (SAFETY rule), so the
         biological signal is protected, not imputed away.

Identifiability
---------------
π_tech here is a *fixed, non-learned* function of depth, gene detection rate,
z_perp neighbors, and sealed biological-context metadata. It has no trainable parameters the model could game, hence
no collapse risk (contrast the tech adversary, which collapsed to the majority
class). It enters the loss only as a detached soft-target / weight. The natural
place to *calibrate* it later is binomial-thinning supervision (see
`objectives/thinning.py`), which yields proven technical-zero labels.

See doc/zero_origin_and_capacity_design_2026-06-01.md for the full rationale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence

import numpy as np

import torch


@dataclass(frozen=True)
class LatentKNNPiTechConfig:
    """Configuration for latent-kNN π_tech estimation."""

    enabled: bool = False

    # "in_batch": nearest neighbors are found within the current minibatch.
    # "same_celltype_donor_bank": query a detached, train-only epoch bank;
    # neighbors share cell type (and region by default), and each donor
    # contributes at most one.
    mode: str = "in_batch"
    k: int = 32                       # bank neighbors per cell
    # Optional legacy in-batch K.  Keeping this at the canonical value while
    # increasing the bank K preserves Phase-I trajectory parity.
    in_batch_k: Optional[int] = None

    # Train-only same-cell-type, distinct-donor latent bank.  The bank is
    # rebuilt from deterministic donor×cell-type slots at epoch boundaries.
    bank_start_epoch: int = 13
    bank_ramp_epochs: int = 4
    bank_refresh_epochs: int = 1
    bank_per_donor_celltype: int = 2
    bank_batch_size: int = 48
    bank_num_workers: int = 0
    bank_min_distinct_donors: int = 4
    bank_detect_shrinkage_donors: float = 8.0
    bank_seed: int = 20260827
    bank_exclude_query_donor: bool = True
    # A query without the required distinct-donor support must fail closed by
    # default.  ``in_batch`` remains available only as an explicit legacy
    # ablation during the bank ramp.
    bank_invalid_fallback: str = "zero"
    # Region is known biological context and is modeled outside z_perp by the
    # PRISM decoder.  Match it explicitly; technology remains unrestricted.
    bank_match_region: bool = True
    # Indices are resolved from Ensembl IDs/gene symbols by a sealed builder
    # before launch.  These genes use same-sex neighbors; unknown-sex queries
    # receive π=0 at these positions.
    sex_linked_gene_indices: tuple[int, ...] = ()

    # Epoch-level thinning diagnostics are emitted only when both classes have
    # adequate support.  Counts are always logged so an omitted metric means
    # unsupported/NA rather than a numerical zero.
    diagnostic_min_positions: int = 256
    diagnostic_min_cells: int = 8
    diagnostic_min_distinct_donors: int = 2

    # Hard upper bound on π_tech: no single zero may be discounted more than this.
    pi_cap: float = 0.6

    # Depth proxy used for the measurement-quality term.
    #   "depth_value" - batch["depth_value"] (library size) if present
    #   "n_detected"  - number of detected genes per cell  (default; always available)
    #   "expm1_sum"   - sum of expm1(x_log1p) per cell (approx total counts)
    depth_proxy: str = "n_detected"
    depth_low_pct: float = 0.25       # cells below this depth-quantile count as "low depth"
    depth_low_temp: float = 0.05      # softness of the low-depth gate

    # Gene global detection-rate gate: a gene is "detectable" between
    # detect_floor and neighbor_on_hi.
    detect_floor: float = 0.05

    # Neighborhood gates.
    neighbor_on_hi: float = 0.5       # neighbor detection rate above which gene is "expected ON"
    safety_neighbor_off: float = 0.2  # at/below this neighbor detection rate -> force π_tech low

    # Convex-ish combination of the two terms.
    w_measure: float = 0.5
    w_neighbor: float = 0.5

    eps: float = 1e-6

    def __post_init__(self) -> None:
        if self.mode not in {"in_batch", "same_celltype_donor_bank"}:
            raise ValueError(f"unsupported pi_tech mode: {self.mode!r}")
        if self.bank_invalid_fallback not in {"zero", "in_batch"}:
            raise ValueError(
                "bank_invalid_fallback must be 'zero' or 'in_batch'"
            )
        if int(self.k) <= 0:
            raise ValueError("pi_tech k must be positive")
        if self.in_batch_k is not None and int(self.in_batch_k) <= 0:
            raise ValueError("pi_tech in_batch_k must be positive when set")
        if not 0.0 <= float(self.pi_cap) <= 1.0:
            raise ValueError("pi_tech pi_cap must be in [0,1]")
        if float(self.neighbor_on_hi) <= float(self.safety_neighbor_off):
            raise ValueError("neighbor_on_hi must exceed safety_neighbor_off")
        if float(self.neighbor_on_hi) <= float(self.detect_floor):
            raise ValueError("neighbor_on_hi must exceed detect_floor")
        if self.mode == "same_celltype_donor_bank":
            positive = {
                "bank_start_epoch": self.bank_start_epoch,
                "bank_refresh_epochs": self.bank_refresh_epochs,
                "bank_per_donor_celltype": self.bank_per_donor_celltype,
                "bank_batch_size": self.bank_batch_size,
                "bank_min_distinct_donors": self.bank_min_distinct_donors,
            }
            for name, value in positive.items():
                if int(value) <= 0:
                    raise ValueError(f"{name} must be positive")
            if int(self.bank_ramp_epochs) < 0 or int(self.bank_num_workers) < 0:
                raise ValueError("bank ramp epochs/workers must be non-negative")
            if int(self.bank_min_distinct_donors) > int(self.k):
                raise ValueError("bank_min_distinct_donors must not exceed k")
            if float(self.bank_detect_shrinkage_donors) < 0.0:
                raise ValueError("bank_detect_shrinkage_donors must be non-negative")
        for name in (
            "diagnostic_min_positions",
            "diagnostic_min_cells",
            "diagnostic_min_distinct_donors",
        ):
            if int(getattr(self, name)) <= 0:
                raise ValueError(f"{name} must be positive")
        indices = tuple(int(index) for index in self.sex_linked_gene_indices)
        if any(index < 0 for index in indices) or len(indices) != len(set(indices)):
            raise ValueError("sex_linked_gene_indices must be unique and non-negative")


def _rank_quantile(x: torch.Tensor) -> torch.Tensor:
    """Map values to their rank-based quantile in [0, 1] (0 = smallest)."""
    n = int(x.shape[0])
    if n <= 1:
        return torch.zeros_like(x)
    order = torch.argsort(x)
    ranks = torch.empty_like(order)
    ranks[order] = torch.arange(n, device=x.device)
    return ranks.to(dtype=x.dtype) / float(n - 1)


def depth_proxy_from_batch(
    batch: Dict[str, torch.Tensor],
    y_ord: torch.Tensor,
    config: LatentKNNPiTechConfig,
) -> torch.Tensor:
    """Per-cell depth proxy [B], robust to the missing-depth dataset."""
    if config.depth_proxy == "depth_value":
        dv = batch.get("depth_value", None)
        if dv is not None:
            return dv.to(dtype=torch.float32).view(-1)
    if config.depth_proxy == "expm1_sum":
        x_log1p = batch.get("x_log1p", None)
        if x_log1p is not None:
            return torch.expm1(x_log1p.to(dtype=torch.float32)).sum(dim=1)
    # default / fallback: number of detected genes (always derivable from y_ord)
    return (y_ord > 0).to(dtype=torch.float32).sum(dim=1)


def _pi_from_neighbor_rate(
    y_ord: torch.Tensor,
    depth: torch.Tensor,
    gene_detect_rate: torch.Tensor,
    neighbor_rate: torch.Tensor,
    config: LatentKNNPiTechConfig,
    *,
    depth_quantile: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """Shared detached score transform for in-batch and bank neighbors."""

    depth_q = (
        _rank_quantile(depth.to(dtype=torch.float32))
        if depth_quantile is None
        else depth_quantile.to(device=y_ord.device, dtype=torch.float32)
    )
    if depth_q.shape != (y_ord.shape[0],):
        raise ValueError("depth_quantile must have shape [B]")
    low_depth = torch.sigmoid(
        (float(config.depth_low_pct) - depth_q)
        / max(float(config.depth_low_temp), 1.0e-6)
    )
    span_g = max(
        float(config.neighbor_on_hi) - float(config.detect_floor),
        float(config.eps),
    )
    gdr = gene_detect_rate.to(device=y_ord.device, dtype=torch.float32)
    if gdr.ndim == 1:
        gdr = gdr.unsqueeze(0)
    if gdr.shape not in {(1, y_ord.shape[1]), tuple(y_ord.shape)}:
        raise ValueError(
            "gene_detect_rate must have shape [G] or [B,G], got "
            f"{tuple(gdr.shape)}"
        )
    detectable = ((gdr - float(config.detect_floor)) / span_g).clamp(0.0, 1.0)
    m_measure = low_depth.unsqueeze(1) * detectable
    span_n = max(
        float(config.neighbor_on_hi) - float(config.safety_neighbor_off),
        float(config.eps),
    )
    m_neighbor = (
        (neighbor_rate.to(torch.float32) - float(config.safety_neighbor_off))
        / span_n
    ).clamp(0.0, 1.0)
    pi_raw = (
        float(config.w_measure) * m_measure
        + float(config.w_neighbor) * m_neighbor
    )
    pi = torch.where(
        neighbor_rate <= float(config.safety_neighbor_off),
        torch.zeros_like(pi_raw),
        pi_raw,
    )
    pi = pi.clamp(min=0.0, max=float(config.pi_cap))
    return (pi * (y_ord == 0).to(pi.dtype)).detach()


@torch.no_grad()
def compute_pi_tech_in_batch(
    z_perp: torch.Tensor,           # [B, d_z]
    y_ord: torch.Tensor,            # [B, G] long
    depth: torch.Tensor,            # [B] float
    gene_detect_rate: torch.Tensor,  # [G] in [0, 1], global EMA detection rate
    config: LatentKNNPiTechConfig,
) -> torch.Tensor:
    """
    Estimate π_tech[i,g] in [0, pi_cap], nonzero only where y_ord[i,g] == 0.

    Returns a detached [B, G] tensor (π_tech is used as a target/weight, never
    backprop'd through neighbor selection).
    """
    B, G = int(y_ord.shape[0]), int(y_ord.shape[1])
    device = y_ord.device
    detect = (y_ord > 0).to(dtype=torch.float32)                     # [B, G]

    # ---- neighborhood term: neighbor detection rate per gene ----
    if B >= 2:
        z = z_perp.detach().to(dtype=torch.float32)
        requested_k = (
            int(config.in_batch_k)
            if config.in_batch_k is not None
            else int(config.k)
        )
        k = max(1, min(requested_k, B - 1))
        d2 = torch.cdist(z, z)                                       # [B, B]
        # exclude self without in-place diagonal fill (keeps it simple/safe)
        d2 = d2 + torch.eye(B, device=device, dtype=d2.dtype) * 1.0e9
        nbr = torch.topk(d2, k, dim=1, largest=False).indices        # [B, k]
        # membership matrix @ detect avoids a [B, k, G] gather (memory-friendly)
        member = torch.zeros(B, B, device=device, dtype=torch.float32)
        member.scatter_(1, nbr, 1.0)
        nb_rate = (member @ detect) / float(k)                       # [B, G]
    else:
        nb_rate = torch.zeros_like(detect)

    return _pi_from_neighbor_rate(
        y_ord, depth, gene_detect_rate, nb_rate, config
    )


def select_donor_celltype_bank_indices(
    *,
    row_indices: Sequence[int],
    celltype_ids: Sequence[int],
    donor_ids: Sequence[int],
    region_ids: Optional[Sequence[int]] = None,
    per_group: int,
    seed: int,
) -> np.ndarray:
    """Select deterministic local indices per donor×cell type, region-stratified."""

    if int(per_group) <= 0:
        raise ValueError("per_group must be positive")
    rows = np.asarray(row_indices, dtype=np.int64)
    celltypes = np.asarray(celltype_ids)
    donors = np.asarray(donor_ids)
    regions = None if region_ids is None else np.asarray(region_ids)
    if rows.ndim != 1 or celltypes.ndim != 1 or donors.ndim != 1:
        raise ValueError("bank row/celltype/donor arrays must be one-dimensional")
    if celltypes.shape != donors.shape:
        raise ValueError("bank celltype/donor arrays must share the full-row shape")
    if regions is not None and regions.shape != celltypes.shape:
        raise ValueError("bank region array must share the full-row shape")
    if rows.size == 0 or bool((rows < 0).any()) or bool(
        (rows >= celltypes.shape[0]).any()
    ):
        raise ValueError("bank row indices are empty or outside the full-row arrays")
    if bool((donors[rows] < 0).any()):
        raise ValueError("bank donor ids must be non-negative")
    groups: dict[tuple[int, int], list[int]] = {}
    for local_index, row in enumerate(rows):
        groups.setdefault(
            (int(celltypes[row]), int(donors[row])), []
        ).append(int(local_index))
    rng = np.random.default_rng(int(seed))
    selected: list[np.ndarray] = []
    for key in sorted(groups):
        candidates = np.asarray(groups[key], dtype=np.int64)
        take = min(int(per_group), int(candidates.size))
        if regions is None or take <= 1:
            selected.append(rng.choice(candidates, size=take, replace=False))
            continue
        # Guarantee region coverage before filling remaining slots.  This
        # keeps a two-slot donor×cell-type bank from accidentally selecting
        # both cells from only DLPFC or only MTG when both are available.
        candidate_regions = regions[rows[candidates]]
        region_values = np.asarray(
            sorted(set(int(value) for value in candidate_regions.tolist())),
            dtype=np.int64,
        )
        region_values = rng.permutation(region_values)
        chosen: list[int] = []
        for region in region_values[:take]:
            pool = candidates[candidate_regions == int(region)]
            chosen.append(int(rng.choice(pool)))
        remaining_slots = take - len(chosen)
        if remaining_slots > 0:
            remaining = np.asarray(
                [value for value in candidates.tolist() if value not in chosen],
                dtype=np.int64,
            )
            chosen.extend(
                int(value)
                for value in rng.choice(
                    remaining, size=remaining_slots, replace=False
                ).tolist()
            )
        selected.append(np.asarray(chosen, dtype=np.int64))
    if not selected:
        raise ValueError("train-only bank selection found no donor×cell-type groups")
    return np.sort(np.concatenate(selected)).astype(np.int64)


def _pack_detection(detect: torch.Tensor) -> tuple[torch.Tensor, int]:
    detect = detect.detach().to(device="cpu", dtype=torch.uint8)
    n_genes = int(detect.shape[1])
    pad = (-n_genes) % 8
    if pad:
        detect = torch.cat(
            [detect, torch.zeros(detect.shape[0], pad, dtype=torch.uint8)], dim=1
        )
    weights = (2 ** torch.arange(8, dtype=torch.int16)).view(1, 1, 8)
    packed = (
        detect.reshape(detect.shape[0], -1, 8).to(torch.int16) * weights
    ).sum(dim=2).to(torch.uint8)
    return packed, n_genes


def _unpack_detection(packed: torch.Tensor, n_genes: int) -> torch.Tensor:
    shifts = torch.arange(8, dtype=torch.uint8).view(1, 1, 8)
    bits = ((packed.to(torch.uint8).unsqueeze(-1) >> shifts) & 1).reshape(
        packed.shape[0], -1
    )
    return bits[:, : int(n_genes)].contiguous()


@dataclass
class TrainOnlyPiTechBank:
    """Detached context-matched latent bank with one vote per distinct donor."""

    z: torch.Tensor
    detect: torch.Tensor
    row: torch.Tensor
    celltype: torch.Tensor
    region: torch.Tensor
    donor: torch.Tensor
    sex: torch.Tensor
    epoch: int

    def __post_init__(self) -> None:
        n = int(self.z.shape[0])
        if self.z.ndim != 2 or self.detect.ndim != 2:
            raise ValueError("bank z/detect must have shapes [N,D] and [N,G]")
        if int(self.detect.shape[0]) != n:
            raise ValueError("bank z/detect row count mismatch")
        if int(self.detect.shape[1]) == 0:
            raise ValueError("technical-zero bank must contain at least one gene")
        for name in ("row", "celltype", "region", "donor", "sex"):
            value = getattr(self, name)
            if value.shape != (n,):
                raise ValueError(f"bank {name} must have shape [N]")
        if n == 0:
            raise ValueError("technical-zero bank must not be empty")
        if not bool(torch.isfinite(self.z.float()).all()):
            raise ValueError("technical-zero bank z contains non-finite values")
        if bool(((self.detect != 0) & (self.detect != 1)).any()):
            raise ValueError("technical-zero bank detect must be binary")
        if bool((self.donor < 0).any()):
            raise ValueError("technical-zero bank donor ids must be non-negative")
        # Query-only sufficient statistics are deterministic functions of the
        # epoch bank.  They are rebuilt lazily after refresh/restore and are
        # intentionally excluded from checkpoints.
        self._query_statistics_cache: dict[
            tuple[bool, tuple[int, ...]], dict[str, object]
        ] = {}
        self._query_statistics_cache_builds = 0

    @property
    def n_genes(self) -> int:
        return int(self.detect.shape[1])

    def to(self, device: torch.device | str) -> "TrainOnlyPiTechBank":
        return TrainOnlyPiTechBank(
            z=self.z.to(device=device, dtype=torch.float32),
            detect=self.detect.to(device=device, dtype=torch.uint8),
            row=self.row.to(device=device, dtype=torch.long),
            celltype=self.celltype.to(device=device, dtype=torch.long),
            region=self.region.to(device=device, dtype=torch.long),
            donor=self.donor.to(device=device, dtype=torch.long),
            sex=self.sex.to(device=device, dtype=torch.long),
            epoch=int(self.epoch),
        )

    def state_dict(self) -> dict[str, object]:
        packed, n_genes = _pack_detection(self.detect)
        return {
            "schema_version": "kmlee_bam.pi_tech_train_bank.v2",
            "z": self.z.detach().cpu().float(),
            "detect_packed": packed,
            "n_genes": int(n_genes),
            "row": self.row.detach().cpu().long(),
            "celltype": self.celltype.detach().cpu().long(),
            "region": self.region.detach().cpu().long(),
            "donor": self.donor.detach().cpu().long(),
            "sex": self.sex.detach().cpu().long(),
            "epoch": int(self.epoch),
        }

    @classmethod
    def from_state_dict(
        cls, state: dict[str, object], *, device: torch.device | str
    ) -> "TrainOnlyPiTechBank":
        if state.get("schema_version") != "kmlee_bam.pi_tech_train_bank.v2":
            raise ValueError("unsupported technical-zero bank checkpoint schema")
        packed = torch.as_tensor(state["detect_packed"], dtype=torch.uint8)
        bank = cls(
            z=torch.as_tensor(state["z"], dtype=torch.float32),
            detect=_unpack_detection(packed, int(state["n_genes"])),
            row=torch.as_tensor(state["row"], dtype=torch.long),
            celltype=torch.as_tensor(state["celltype"], dtype=torch.long),
            region=torch.as_tensor(state["region"], dtype=torch.long),
            donor=torch.as_tensor(state["donor"], dtype=torch.long),
            sex=torch.as_tensor(state["sex"], dtype=torch.long),
            epoch=int(state["epoch"]),
        )
        return bank.to(device)

    def _neighbor_rate(
        self,
        query_z: torch.Tensor,
        query_celltype: torch.Tensor,
        query_row: torch.Tensor,
        *,
        k: int,
        minimum_distinct_donors: int,
        query_sex: Optional[torch.Tensor] = None,
        query_donor: Optional[torch.Tensor] = None,
        query_region: Optional[torch.Tensor] = None,
        gene_indices: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        device = query_z.device
        bank = self if self.z.device == device else self.to(device)
        distance = torch.cdist(query_z.detach().float(), bank.z.float())
        allowed = query_celltype.long()[:, None].eq(bank.celltype[None, :])
        allowed &= query_row.long()[:, None].ne(bank.row[None, :])
        if query_region is not None:
            allowed &= query_region.long()[:, None].eq(bank.region[None, :])
        if query_donor is not None:
            allowed &= query_donor.long()[:, None].ne(bank.donor[None, :])
        if query_sex is not None:
            known = query_sex.long() >= 0
            allowed &= known[:, None]
            allowed &= query_sex.long()[:, None].eq(bank.sex[None, :])
        distance = distance.masked_fill(~allowed, torch.inf)
        n_donors = int(bank.donor.max().item()) + 1
        donor_distance = torch.full(
            (query_z.shape[0], n_donors),
            torch.inf,
            device=device,
            dtype=distance.dtype,
        )
        donor_distance.scatter_reduce_(
            1,
            bank.donor.unsqueeze(0).expand(query_z.shape[0], -1),
            distance,
            reduce="amin",
            include_self=True,
        )
        resolved_k = max(1, min(int(k), n_donors))
        donor_value, chosen_donor = torch.topk(
            donor_distance, resolved_k, dim=1, largest=False
        )
        valid_vote = torch.isfinite(donor_value)
        candidate = chosen_donor[:, :, None].eq(bank.donor[None, None, :])
        candidate &= allowed[:, None, :]
        chosen_distance = distance[:, None, :].masked_fill(
            ~candidate, torch.inf
        )
        neighbor_index = chosen_distance.argmin(dim=2)
        detect_source = (
            bank.detect
            if gene_indices is None
            else bank.detect[:, gene_indices.to(device=device, dtype=torch.long)]
        )
        gathered = detect_source[neighbor_index].float()
        gathered = gathered * valid_vote[:, :, None].to(gathered.dtype)
        effective_k = valid_vote.sum(dim=1)
        neighbor_rate = gathered.sum(dim=1) / effective_k.clamp_min(1).to(
            gathered.dtype
        ).unsqueeze(1)
        valid_query = effective_k >= int(minimum_distinct_donors)
        return neighbor_rate, valid_query, effective_k

    @staticmethod
    def _donor_balanced_mean(
        values: torch.Tensor,
        donor: torch.Tensor,
    ) -> tuple[torch.Tensor, int]:
        """Average rows within donor first, then give every donor equal mass."""

        unique_donors = torch.unique(donor.long(), sorted=True)
        if int(unique_donors.numel()) == 0:
            raise ValueError("donor-balanced mean requires at least one donor")
        donor_means = torch.stack(
            [
                values[donor.long() == int(donor_id)].float().mean(dim=0)
                for donor_id in unique_donors.tolist()
            ],
            dim=0,
        )
        return donor_means.mean(dim=0), int(unique_donors.numel())

    def _query_gene_detect_rate_reference(
        self,
        global_rate: torch.Tensor,
        query_celltype: torch.Tensor,
        query_region: torch.Tensor,
        query_sex: torch.Tensor,
        *,
        shrinkage_donors: float,
        sex_linked_gene_indices: Sequence[int],
        match_region: bool,
    ) -> torch.Tensor:
        bank = self if self.z.device == global_rate.device else self.to(global_rate.device)
        global_rate = global_rate.float()
        result = global_rate.unsqueeze(0).expand(query_celltype.shape[0], -1).clone()
        for celltype in torch.unique(query_celltype.long()).tolist():
            cell_query = query_celltype == int(celltype)
            regions: list[Optional[int]] = (
                [
                    int(value)
                    for value in torch.unique(
                        query_region[cell_query].long(), sorted=True
                    ).tolist()
                ]
                if bool(match_region)
                else [None]
            )
            for region in regions:
                query_mask = cell_query
                bank_mask = bank.celltype == int(celltype)
                if region is not None:
                    query_mask = query_mask & (query_region == int(region))
                    bank_mask = bank_mask & (bank.region == int(region))
                if not bool(query_mask.any()) or not bool(bank_mask.any()):
                    continue
                rate, distinct = self._donor_balanced_mean(
                    bank.detect[bank_mask], bank.donor[bank_mask]
                )
                weight = distinct / (distinct + float(shrinkage_donors))
                result[query_mask] = (
                    weight * rate + (1.0 - weight) * global_rate
                )
        indices = [
            int(index)
            for index in sex_linked_gene_indices
            if 0 <= int(index) < self.n_genes
        ]
        if indices:
            gene_index = torch.tensor(indices, device=result.device, dtype=torch.long)
            for celltype in torch.unique(query_celltype.long()).tolist():
                cell_query = query_celltype == int(celltype)
                regions = (
                    [
                        int(value)
                        for value in torch.unique(
                            query_region[cell_query].long(), sorted=True
                        ).tolist()
                    ]
                    if bool(match_region)
                    else [None]
                )
                for region in regions:
                    for sex in (0, 1):
                        query_mask = cell_query & (query_sex == int(sex))
                        bank_mask = (bank.celltype == int(celltype)) & (
                            bank.sex == int(sex)
                        )
                        if region is not None:
                            query_mask = query_mask & (
                                query_region == int(region)
                            )
                            bank_mask = bank_mask & (bank.region == int(region))
                        if not bool(query_mask.any()) or not bool(bank_mask.any()):
                            continue
                        rate, distinct = self._donor_balanced_mean(
                            bank.detect[bank_mask][:, gene_index],
                            bank.donor[bank_mask],
                        )
                        weight = distinct / (
                            distinct + float(shrinkage_donors)
                        )
                        rows = torch.nonzero(
                            query_mask, as_tuple=False
                        ).flatten()
                        result[rows[:, None], gene_index[None, :]] = (
                            weight * rate
                            + (1.0 - weight) * global_rate[gene_index]
                        ).unsqueeze(0)
        return result

    def _query_depth_quantile_reference(
        self,
        query_depth: torch.Tensor,
        query_celltype: torch.Tensor,
        query_region: torch.Tensor,
        *,
        match_region: bool,
    ) -> torch.Tensor:
        """Midrank depth percentile within a donor-balanced cell-type bank."""

        bank = self if self.z.device == query_depth.device else self.to(query_depth.device)
        result = _rank_quantile(query_depth.float())
        bank_depth = bank.detect.sum(dim=1).float()
        for celltype in torch.unique(query_celltype.long()).tolist():
            cell_query = query_celltype.long() == int(celltype)
            regions: list[Optional[int]] = (
                [
                    int(value)
                    for value in torch.unique(
                        query_region[cell_query].long(), sorted=True
                    ).tolist()
                ]
                if bool(match_region)
                else [None]
            )
            for region in regions:
                query_mask = cell_query
                bank_mask = bank.celltype == int(celltype)
                if region is not None:
                    query_mask = query_mask & (query_region == int(region))
                    bank_mask = bank_mask & (bank.region == int(region))
                reference = bank_depth[bank_mask]
                if int(reference.numel()) < 2 or not bool(query_mask.any()):
                    continue
                query_value = query_depth[query_mask].float().unsqueeze(1)
                reference_donor = bank.donor[bank_mask]
                donor_less = []
                donor_equal = []
                for donor_id in torch.unique(
                    reference_donor.long(), sorted=True
                ).tolist():
                    donor_reference = reference[
                        reference_donor.long() == int(donor_id)
                    ]
                    donor_less.append(
                        (donor_reference.unsqueeze(0) < query_value)
                        .float()
                        .mean(dim=1)
                    )
                    donor_equal.append(
                        (donor_reference.unsqueeze(0) == query_value)
                        .float()
                        .mean(dim=1)
                    )
                less = torch.stack(donor_less, dim=0).mean(dim=0)
                equal = torch.stack(donor_equal, dim=0).mean(dim=0)
                result[query_mask] = (less + 0.5 * equal).clamp(0.0, 1.0)
        return result

    def _build_query_statistics(
        self,
        *,
        match_region: bool,
        sex_linked_gene_indices: Sequence[int],
    ) -> dict[str, object]:
        """Cache epoch-frozen bank statistics, never the live EMA shrinkage."""

        sex_indices = tuple(
            int(index)
            for index in sex_linked_gene_indices
            if 0 <= int(index) < self.n_genes
        )
        cache_key = (bool(match_region), sex_indices)
        cached = self._query_statistics_cache.get(cache_key)
        if cached is not None:
            return cached

        device = self.z.device
        celltype_limit = int(self.celltype.max().item()) + 1
        if bool(match_region):
            region_limit = max(1, int(self.region.max().item()) + 2)
            region_slot = self.region.long() + 1
        else:
            region_limit = 1
            region_slot = torch.zeros_like(self.region, dtype=torch.long)
        group_lookup = torch.full(
            (celltype_limit, region_limit),
            -1,
            device=device,
            dtype=torch.long,
        )
        group_key_tensor = torch.stack(
            (self.celltype.long(), region_slot), dim=1
        )
        group_keys = torch.unique(group_key_tensor, dim=0, sorted=True)
        bank_depth = self.detect.sum(dim=1).long()
        rate_rows = []
        donor_counts = []
        depth_rows = []
        depth_valid = []
        sex_rate_rows = []
        sex_donor_rows = []
        gene_index = torch.tensor(
            sex_indices, device=device, dtype=torch.long
        )

        for group_index, key in enumerate(group_keys.tolist()):
            celltype, region_index = int(key[0]), int(key[1])
            mask = (self.celltype.long() == celltype) & (
                region_slot == region_index
            )
            rate, distinct = self._donor_balanced_mean(
                self.detect[mask], self.donor[mask]
            )
            rate_rows.append(rate)
            donor_counts.append(float(distinct))
            group_lookup[celltype, region_index] = int(group_index)

            reference_depth = bank_depth[mask]
            if int(reference_depth.numel()) >= 2:
                donor_depth_rows = []
                reference_donor = self.donor[mask].long()
                for donor_id in torch.unique(
                    reference_donor, sorted=True
                ).tolist():
                    donor_depth = reference_depth[
                        reference_donor == int(donor_id)
                    ]
                    histogram = torch.bincount(
                        donor_depth,
                        minlength=self.n_genes + 1,
                    ).float()
                    less = torch.cumsum(histogram, dim=0) - histogram
                    donor_depth_rows.append(
                        (less + 0.5 * histogram)
                        / float(donor_depth.numel())
                    )
                depth_rows.append(torch.stack(donor_depth_rows).mean(dim=0))
                depth_valid.append(True)
            else:
                depth_rows.append(
                    torch.zeros(
                        self.n_genes + 1,
                        device=device,
                        dtype=torch.float32,
                    )
                )
                depth_valid.append(False)

            group_sex_rate = torch.zeros(
                2, len(sex_indices), device=device, dtype=torch.float32
            )
            group_sex_donors = torch.zeros(
                2, device=device, dtype=torch.float32
            )
            if sex_indices:
                for sex in (0, 1):
                    sex_mask = mask & (self.sex.long() == int(sex))
                    if bool(sex_mask.any()):
                        sex_rate, sex_distinct = self._donor_balanced_mean(
                            self.detect[sex_mask][:, gene_index],
                            self.donor[sex_mask],
                        )
                        group_sex_rate[sex] = sex_rate
                        group_sex_donors[sex] = float(sex_distinct)
            sex_rate_rows.append(group_sex_rate)
            sex_donor_rows.append(group_sex_donors)

        cached = {
            "match_region": bool(match_region),
            "group_lookup": group_lookup,
            "rate": torch.stack(rate_rows).float(),
            "donor_count": torch.tensor(
                donor_counts, device=device, dtype=torch.float32
            ),
            "depth_midrank": torch.stack(depth_rows).float(),
            "depth_valid": torch.tensor(
                depth_valid, device=device, dtype=torch.bool
            ),
            "sex_rate": torch.stack(sex_rate_rows).float(),
            "sex_donor_count": torch.stack(sex_donor_rows).float(),
            "sex_indices": sex_indices,
        }
        self._query_statistics_cache[cache_key] = cached
        self._query_statistics_cache_builds += 1
        return cached

    @staticmethod
    def _query_group_index(
        cache: dict[str, object],
        query_celltype: torch.Tensor,
        query_region: torch.Tensor,
    ) -> torch.Tensor:
        lookup = cache["group_lookup"]
        assert isinstance(lookup, torch.Tensor)
        celltype = query_celltype.long()
        region = (
            query_region.long() + 1
            if bool(cache["match_region"])
            else torch.zeros_like(celltype)
        )
        valid = (
            (celltype >= 0)
            & (celltype < int(lookup.shape[0]))
            & (region >= 0)
            & (region < int(lookup.shape[1]))
        )
        safe_celltype = celltype.clamp(0, int(lookup.shape[0]) - 1)
        safe_region = region.clamp(0, int(lookup.shape[1]) - 1)
        resolved = lookup[safe_celltype, safe_region]
        return torch.where(valid, resolved, torch.full_like(resolved, -1))

    def query_gene_detect_rate(
        self,
        global_rate: torch.Tensor,
        query_celltype: torch.Tensor,
        query_region: torch.Tensor,
        query_sex: torch.Tensor,
        *,
        shrinkage_donors: float,
        sex_linked_gene_indices: Sequence[int],
        match_region: bool,
    ) -> torch.Tensor:
        """Vectorized cached query with step-live global EMA shrinkage."""

        bank = self if self.z.device == global_rate.device else self.to(global_rate.device)
        global_rate = global_rate.float()
        if global_rate.shape != (bank.n_genes,):
            raise ValueError("global_rate must have shape [G]")
        cache = bank._build_query_statistics(
            match_region=bool(match_region),
            sex_linked_gene_indices=sex_linked_gene_indices,
        )
        group_index = bank._query_group_index(
            cache, query_celltype, query_region
        )
        group_valid = group_index >= 0
        safe_group = group_index.clamp_min(0)
        rate = cache["rate"]
        donor_count = cache["donor_count"]
        assert isinstance(rate, torch.Tensor)
        assert isinstance(donor_count, torch.Tensor)
        selected_rate = rate[safe_group]
        selected_donors = donor_count[safe_group]
        denominator = (
            selected_donors + float(shrinkage_donors)
        ).clamp_min(1.0e-12)
        weight = selected_donors / denominator
        shrunk = (
            weight[:, None] * selected_rate
            + (1.0 - weight[:, None]) * global_rate[None, :]
        )
        result = torch.where(
            group_valid[:, None],
            shrunk,
            global_rate[None, :].expand(query_celltype.shape[0], -1),
        ).clone()

        sex_indices = cache["sex_indices"]
        assert isinstance(sex_indices, tuple)
        if sex_indices:
            gene_index = torch.tensor(
                sex_indices, device=result.device, dtype=torch.long
            )
            sex_rate = cache["sex_rate"]
            sex_donor_count = cache["sex_donor_count"]
            assert isinstance(sex_rate, torch.Tensor)
            assert isinstance(sex_donor_count, torch.Tensor)
            sex = query_sex.long()
            known_sex = (sex == 0) | (sex == 1)
            safe_sex = sex.clamp(0, 1)
            selected_sex_rate = sex_rate[safe_group, safe_sex]
            selected_sex_donors = sex_donor_count[safe_group, safe_sex]
            sex_denominator = (
                selected_sex_donors + float(shrinkage_donors)
            ).clamp_min(1.0e-12)
            sex_weight = selected_sex_donors / sex_denominator
            sex_shrunk = (
                sex_weight[:, None] * selected_sex_rate
                + (1.0 - sex_weight[:, None])
                * global_rate[gene_index][None, :]
            )
            sex_valid = (
                group_valid & known_sex & (selected_sex_donors > 0.0)
            )
            result[:, gene_index] = torch.where(
                sex_valid[:, None],
                sex_shrunk,
                result[:, gene_index],
            )
        return result

    def query_depth_quantile(
        self,
        query_depth: torch.Tensor,
        query_celltype: torch.Tensor,
        query_region: torch.Tensor,
        *,
        match_region: bool,
        sex_linked_gene_indices: Sequence[int] = (),
    ) -> torch.Tensor:
        """Vectorized exact lookup for integer detected-gene depth."""

        bank = self if self.z.device == query_depth.device else self.to(query_depth.device)
        cache = bank._build_query_statistics(
            match_region=bool(match_region),
            sex_linked_gene_indices=sex_linked_gene_indices,
        )
        group_index = bank._query_group_index(
            cache, query_celltype, query_region
        )
        group_valid = group_index >= 0
        safe_group = group_index.clamp_min(0)
        depth = query_depth.float()
        depth_index = depth.round().long().clamp(0, bank.n_genes)
        integer_depth = (depth - depth_index.float()).abs() <= 1.0e-6
        depth_midrank = cache["depth_midrank"]
        depth_valid = cache["depth_valid"]
        assert isinstance(depth_midrank, torch.Tensor)
        assert isinstance(depth_valid, torch.Tensor)
        cached_value = depth_midrank[safe_group, depth_index]
        use_cache = group_valid & depth_valid[safe_group] & integer_depth
        return torch.where(
            use_cache,
            cached_value,
            _rank_quantile(depth),
        ).clamp(0.0, 1.0)


@torch.no_grad()
def compute_pi_tech_from_bank(
    z_perp: torch.Tensor,
    y_ord: torch.Tensor,
    depth: torch.Tensor,
    gene_detect_rate: torch.Tensor,
    celltype_id: torch.Tensor,
    region_id: torch.Tensor,
    donor_id: torch.Tensor,
    row_index: torch.Tensor,
    sex_id: torch.Tensor,
    bank: TrainOnlyPiTechBank,
    config: LatentKNNPiTechConfig,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    """Estimate π from a cell-type/region bank with distinct-donor votes."""

    bank = bank if bank.z.device == y_ord.device else bank.to(y_ord.device)
    neighbor_rate, valid, effective_k = bank._neighbor_rate(
        z_perp,
        celltype_id,
        row_index,
        k=int(config.k),
        minimum_distinct_donors=int(config.bank_min_distinct_donors),
        query_donor=(
            donor_id if bool(config.bank_exclude_query_donor) else None
        ),
        query_region=(region_id if bool(config.bank_match_region) else None),
    )
    sex_indices = [
        int(index)
        for index in config.sex_linked_gene_indices
        if 0 <= int(index) < y_ord.shape[1]
    ]
    sex_valid = valid.clone()
    if sex_indices:
        gene_index = torch.tensor(
            sex_indices, device=y_ord.device, dtype=torch.long
        )
        sex_rate, sex_valid, _ = bank._neighbor_rate(
            z_perp,
            celltype_id,
            row_index,
            k=int(config.k),
            minimum_distinct_donors=int(config.bank_min_distinct_donors),
            query_sex=sex_id,
            query_donor=(
                donor_id if bool(config.bank_exclude_query_donor) else None
            ),
            query_region=(
                region_id if bool(config.bank_match_region) else None
            ),
            gene_indices=gene_index,
        )
        neighbor_rate[:, gene_index] = sex_rate
    query_rate = bank.query_gene_detect_rate(
        gene_detect_rate,
        celltype_id,
        region_id,
        sex_id,
        shrinkage_donors=float(config.bank_detect_shrinkage_donors),
        sex_linked_gene_indices=sex_indices,
        match_region=bool(config.bank_match_region),
    )
    depth_quantile = (
        bank.query_depth_quantile(
            depth,
            celltype_id,
            region_id,
            match_region=bool(config.bank_match_region),
            sex_linked_gene_indices=sex_indices,
        )
        if str(config.depth_proxy) == "n_detected"
        else None
    )
    pi = _pi_from_neighbor_rate(
        y_ord,
        depth,
        query_rate,
        neighbor_rate,
        config,
        depth_quantile=depth_quantile,
    )
    pi = torch.where(valid[:, None], pi, torch.zeros_like(pi))
    if sex_indices:
        gene_index = torch.tensor(
            sex_indices, device=y_ord.device, dtype=torch.long
        )
        safe = sex_valid & (sex_id >= 0)
        pi[:, gene_index] = torch.where(
            safe[:, None], pi[:, gene_index], torch.zeros_like(pi[:, gene_index])
        )
    return pi.detach(), {
        "valid_query": valid,
        "sex_valid_query": sex_valid,
        "effective_k": effective_k,
        "neighbor_rate": neighbor_rate,
        "depth_quantile": (
            depth_quantile
            if depth_quantile is not None
            else _rank_quantile(depth.float())
        ),
    }
