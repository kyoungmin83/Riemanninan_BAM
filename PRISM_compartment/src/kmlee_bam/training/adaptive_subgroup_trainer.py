"""
V8 trainer — Zero-Origin + Capacity.

Sits at the top of the trainer MRO (V8 → V7a → V6 → V4 → V3 → core) and wires
three default-OFF mechanisms on top of the fully-assembled v7a/v6/hierarchical
loss; everything below is untouched when v8 is disabled. (Per-rank coefficients
live in the decoder itself and are controlled by `decoder.per_rank_coeff`, so
they need no trainer hook.)

  * π_tech soft-zero — computes a per-(cell,gene) technical-dropout probability
    from depth + latent-kNN neighbors and stashes it on `_pending_pi_tech` just
    before `super()._compute_loss`, so V4's single hierarchical-loss call relaxes
    the L1 zero penalty on suspicious zeros.
  * Thinning — runs a second forward on the dataset's thinned view and adds a
    `z_perp` depth-invariance consistency loss (+ optional proven-tech-zero
    supervision). Training-only, every Nth step.
  * Adaptive subgroup CVaR — monitors per-(celltype×tech×depth) metrics each
    epoch and auto-engages bounded, ramped loss reweighting on collapsed groups.

See doc/zero_origin_and_capacity_design_2026-06-01.md.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import torch
import torch.distributed as dist
from torch.utils.data import DataLoader, Subset

try:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.pathology_aware_uncertainty_trainer import V7aTrainer
    from kmlee_bam.training.core_trainer import ModelForwardOutput, _resolve_tech_id
    from kmlee_bam.objectives.latent_knn_pi_tech import (
        LatentKNNPiTechConfig,
        TrainOnlyPiTechBank,
        compute_pi_tech_from_bank,
        compute_pi_tech_in_batch,
        depth_proxy_from_batch,
        select_donor_celltype_bank_indices,
    )
    from kmlee_bam.objectives.thinning import (
        ThinningConfig,
        ZeroReliabilityConfig,
        thinning_consistency_loss,
        technical_zero_mask,
        technical_zero_supervision_loss,
        build_thinned_batch,
        zero_reliability_weights,
    )
    from kmlee_bam.objectives.subgroup_cvar import (
        SubgroupCVaRConfig,
        AdaptiveSubgroupCVaR,
        per_cell_metric,
    )
    from kmlee_bam.objectives.variance_floor import (
        VarianceFloorConfig,
        compute_variance_floor,
        RunningNonzeroTierVariance,
    )
    from kmlee_bam.objectives.ranking_loss import (
        HighMarginRankConfig,
        high_margin_loss,
        RunningLowTierScore,
    )
except ImportError:  # pragma: no cover - mirror v7a's import fallback
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.pathology_aware_uncertainty_trainer import V7aTrainer
    from kmlee_bam.training.core_trainer import ModelForwardOutput, _resolve_tech_id
    from kmlee_bam.objectives.latent_knn_pi_tech import (
        LatentKNNPiTechConfig,
        TrainOnlyPiTechBank,
        compute_pi_tech_from_bank,
        compute_pi_tech_in_batch,
        depth_proxy_from_batch,
        select_donor_celltype_bank_indices,
    )
    from kmlee_bam.objectives.thinning import (
        ThinningConfig,
        ZeroReliabilityConfig,
        thinning_consistency_loss,
        technical_zero_mask,
        technical_zero_supervision_loss,
        build_thinned_batch,
        zero_reliability_weights,
    )
    from kmlee_bam.objectives.subgroup_cvar import (
        SubgroupCVaRConfig,
        AdaptiveSubgroupCVaR,
        per_cell_metric,
    )
    from kmlee_bam.objectives.variance_floor import (
        VarianceFloorConfig,
        compute_variance_floor,
        RunningNonzeroTierVariance,
    )
    from kmlee_bam.objectives.ranking_loss import (
        HighMarginRankConfig,
        high_margin_loss,
        RunningLowTierScore,
    )

from kmlee_bam.training.tech_invariance import (
    TechInvarianceConfig,
    TechInvarianceAugmenter,
    z_consistency_loss,
    tech_risk_bce,
    rec_consistency_loss,
)


class V8Trainer(V7aTrainer):
    def __init__(
        self,
        *args,
        v8_enabled: bool = False,
        pi_tech_config: Optional[LatentKNNPiTechConfig] = None,
        thinning_config: Optional[ThinningConfig] = None,
        subgroup_cvar_config: Optional[SubgroupCVaRConfig] = None,
        zero_reliability_config: Optional[ZeroReliabilityConfig] = None,
        variance_floor_config: Optional[VarianceFloorConfig] = None,
        high_margin_config: Optional[HighMarginRankConfig] = None,
        tech_invariance_config: Optional[TechInvarianceConfig] = None,
        n_tech_max: Optional[int] = None,
        n_celltypes: Optional[int] = None,
        n_genes: Optional[int] = None,
        pi_tech_bank_dataset=None,
        gene_detect_ema_momentum: float = 0.95,
        **kwargs,
    ):
        super().__init__(*args, n_celltypes=n_celltypes, n_genes=n_genes, **kwargs)

        self.v8_enabled = bool(v8_enabled)
        self.pi_tech_config = pi_tech_config or LatentKNNPiTechConfig(enabled=False)
        self.thinning_config = thinning_config or ThinningConfig(enabled=False)
        self.subgroup_cvar_config = subgroup_cvar_config or SubgroupCVaRConfig(enabled=False)
        self.zero_reliability_config = zero_reliability_config or ZeroReliabilityConfig(enabled=False)
        self.variance_floor_config = variance_floor_config or VarianceFloorConfig(enabled=False)
        # Stable streaming target for the variance floor (v1.1, target_mode="running");
        # None ⇒ in-batch target (MVP). Created below once the device is known.
        self._var_running = None
        self.high_margin_config = high_margin_config or HighMarginRankConfig(enabled=False)
        self.tech_inv_config = tech_invariance_config or TechInvarianceConfig(enabled=False)
        self.tech_augmenter = None  # lazily built from the train dataset's reference scaler
        self._disease_dir_tensor = None  # optional [C,G] disease direction for the negative control
        self._low_ref = None  # RunningLowTierScore; created below once device is known
        # Stage D conduit: per-(cell,gene) reconstruction-reliability weights read
        # by core_trainer's criterion call. None ⇒ reconstruction byte-identical.
        self._pending_rec_weight_per_gene = None

        # Conduit read by V4Trainer._compute_loss (None ⇒ hierarchical loss unchanged).
        self._pending_pi_tech = None
        self.gene_detect_rate: Optional[torch.Tensor] = None
        self.cvar: Optional[AdaptiveSubgroupCVaR] = None
        self._gene_detect_ema_momentum = float(gene_detect_ema_momentum)
        self._v8_step_counter = 0
        self._v8_sync_group = None
        self._pi_tech_bank_dataset = pi_tech_bank_dataset
        self._pi_tech_bank: Optional[TrainOnlyPiTechBank] = None
        self._pi_tech_bank_local_indices: Optional[np.ndarray] = None
        self._pi_tech_epoch = 0
        self._pi_tech_epoch_diagnostics: Optional[
            Dict[str, torch.Tensor]
        ] = None
        self._pi_tech_diag_n_celltypes = int(n_celltypes or 0)
        self._pi_tech_diag_n_regions = 1
        self._pi_tech_diag_n_donors = 1

        if not self.v8_enabled:
            return

        device = self._infer_device()

        if (
            self.variance_floor_config.enabled
            and str(self.variance_floor_config.target_mode) == "running"
            and n_celltypes is not None
            and n_genes is not None
        ):
            self._var_running = RunningNonzeroTierVariance(
                int(n_celltypes), int(n_genes), device
            )

        if (
            self.high_margin_config.enabled
            and n_celltypes is not None
            and n_genes is not None
        ):
            self._low_ref = RunningLowTierScore(
                int(n_celltypes), int(n_genes), device,
                momentum=self.high_margin_config.ema_momentum,
            )

        if self.pi_tech_config.enabled:
            if n_genes is None:
                raise ValueError("v8 pi_tech enabled requires n_genes.")
            mode = str(self.pi_tech_config.mode)
            if mode not in {"in_batch", "same_celltype_donor_bank"}:
                raise ValueError(
                    "pi_tech.mode must be 'in_batch' or "
                    "'same_celltype_donor_bank'"
                )
            # Neutral prior detection rate until the EMA warms up.
            self.gene_detect_rate = torch.full(
                (int(n_genes),), 0.1, device=device, dtype=torch.float32
            )
            if mode == "same_celltype_donor_bank":
                if self._pi_tech_bank_dataset is None:
                    raise ValueError(
                        "same_celltype_donor_bank requires the train dataset"
                    )
                dataset = self._pi_tech_bank_dataset
                if str(getattr(dataset, "split", "")) != "train":
                    raise ValueError(
                        "technical-zero bank dataset must be the train split"
                    )
                required_attributes = ["row_idx", "celltype_ids", "donor_ids"]
                if bool(self.pi_tech_config.bank_match_region):
                    required_attributes.append("region_ids")
                for attr in required_attributes:
                    if getattr(dataset, attr, None) is None:
                        raise AttributeError(
                            "pi_tech train-bank dataset must expose "
                            f"{attr}"
                        )
                self._pi_tech_bank_local_indices = (
                    select_donor_celltype_bank_indices(
                        row_indices=np.asarray(dataset.row_idx, dtype=np.int64),
                        celltype_ids=np.asarray(dataset.celltype_ids),
                        donor_ids=np.asarray(dataset.donor_ids),
                        region_ids=(
                            np.asarray(dataset.region_ids)
                            if bool(self.pi_tech_config.bank_match_region)
                            else None
                        ),
                        per_group=int(
                            self.pi_tech_config.bank_per_donor_celltype
                        ),
                        seed=int(self.pi_tech_config.bank_seed),
                    )
                )
                selected_rows = np.asarray(dataset.row_idx, dtype=np.int64)
                if self._pi_tech_diag_n_celltypes <= 0:
                    self._pi_tech_diag_n_celltypes = int(
                        np.asarray(dataset.celltype_ids)[selected_rows].max()
                    ) + 1
                if bool(self.pi_tech_config.bank_match_region):
                    self._pi_tech_diag_n_regions = int(
                        np.asarray(dataset.region_ids)[selected_rows].max()
                    ) + 1
                self._pi_tech_diag_n_donors = int(
                    np.asarray(dataset.donor_ids)[selected_rows].max()
                ) + 1

        if self.subgroup_cvar_config.enabled:
            if n_celltypes is None or n_tech_max is None:
                # best-effort inference of n_tech from the decoder
                if n_tech_max is None:
                    dec = getattr(self._unwrapped_system(), "decoder", None)
                    n_tech_max = int(getattr(dec, "n_tech", 0)) or None
                if n_celltypes is None or n_tech_max is None:
                    raise ValueError(
                        "v8 subgroup_cvar enabled requires n_celltypes and n_tech_max."
                    )
            self.cvar = AdaptiveSubgroupCVaR(
                n_celltypes=int(n_celltypes),
                n_tech_max=int(n_tech_max),
                config=self.subgroup_cvar_config,
                device=device,
            )

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def _unwrapped_system(self):
        sys = self.system
        return getattr(sys, "module", sys)

    def set_ddp_control_group(self, group) -> None:
        self._v8_sync_group = group
        if self.cvar is not None:
            self.cvar.set_sync_process_group(group)
        # Preserve V6's reference-prior group wiring.
        parent = getattr(super(), "set_ddp_control_group", None)
        if parent is not None:
            parent(group)

    @torch.no_grad()
    def _update_gene_detect_rate(self, y_ord: torch.Tensor) -> None:
        if self.gene_detect_rate is None:
            return
        batch_rate = (y_ord > 0).to(dtype=torch.float32).mean(dim=0)  # [G]
        m = self._gene_detect_ema_momentum
        self.gene_detect_rate.mul_(m).add_(batch_rate, alpha=1.0 - m)

    @torch.no_grad()
    def _sync_gene_detect_rate(self) -> None:
        if self.gene_detect_rate is None:
            return
        if dist.is_available() and dist.is_initialized():
            dist.all_reduce(self.gene_detect_rate, op=dist.ReduceOp.SUM, group=self._v8_sync_group)
            self.gene_detect_rate /= float(dist.get_world_size())

    def _reset_pi_tech_epoch_diagnostics(self) -> None:
        """Allocate exact position/cell/donor denominators for one epoch."""

        if not self.pi_tech_config.enabled:
            self._pi_tech_epoch_diagnostics = None
            return
        n_groups = max(
            1,
            int(self._pi_tech_diag_n_celltypes)
            * int(self._pi_tech_diag_n_regions),
        )
        device = self._infer_device()

        def _zeros(size) -> torch.Tensor:
            return torch.zeros(size, device=device, dtype=torch.float64)

        self._pi_tech_epoch_diagnostics = {
            "bank_query": _zeros(n_groups),
            "bank_valid": _zeros(n_groups),
            "bank_effective_k": _zeros(n_groups),
            "bank_group_batches": _zeros(n_groups),
            "bank_donor_seen": _zeros(
                (n_groups, int(self._pi_tech_diag_n_donors))
            ),
            "proven_sum": _zeros(n_groups),
            "proven_count": _zeros(n_groups),
            "stable_sum": _zeros(n_groups),
            "stable_count": _zeros(n_groups),
            "thin_cells": _zeros(n_groups),
            "thin_group_batches": _zeros(n_groups),
            "thin_donor_seen": _zeros(
                (n_groups, int(self._pi_tech_diag_n_donors))
            ),
            "thin_probe_steps": _zeros(1),
        }

    def _pi_tech_group_index(
        self, batch: Dict[str, torch.Tensor]
    ) -> tuple[torch.Tensor, torch.Tensor]:
        celltype = batch["celltype_id"].long()
        if bool(self.pi_tech_config.bank_match_region):
            region = batch["region_id"].long()
        else:
            region = torch.zeros_like(celltype)
        valid = (
            (celltype >= 0)
            & (celltype < int(self._pi_tech_diag_n_celltypes))
            & (region >= 0)
            & (region < int(self._pi_tech_diag_n_regions))
        )
        group = celltype * int(self._pi_tech_diag_n_regions) + region
        return group, valid

    @torch.no_grad()
    def _accumulate_pi_tech_bank_diagnostics(
        self,
        batch: Dict[str, torch.Tensor],
        valid_query: torch.Tensor,
        effective_k: torch.Tensor,
    ) -> None:
        state = self._pi_tech_epoch_diagnostics
        if state is None:
            return
        group, valid_group = self._pi_tech_group_index(batch)
        group = group[valid_group]
        if int(group.numel()) == 0:
            return
        ones = torch.ones_like(group, dtype=torch.float64)
        state["bank_query"].index_add_(0, group, ones)
        state["bank_valid"].index_add_(
            0, group, valid_query[valid_group].to(torch.float64)
        )
        state["bank_effective_k"].index_add_(
            0, group, effective_k[valid_group].to(torch.float64)
        )
        unique_group = torch.unique(group)
        state["bank_group_batches"][unique_group] += 1.0
        donor = batch["donor_id"].long()[valid_group]
        valid_donor = (
            (donor >= 0) & (donor < int(self._pi_tech_diag_n_donors))
        )
        state["bank_donor_seen"][
            group[valid_donor], donor[valid_donor]
        ] = 1.0

    @torch.no_grad()
    def _accumulate_pi_tech_thinning_diagnostics(
        self,
        batch: Dict[str, torch.Tensor],
        pi_thin: torch.Tensor,
        proven: torch.Tensor,
        stable: torch.Tensor,
    ) -> None:
        state = self._pi_tech_epoch_diagnostics
        if state is None:
            return
        group, valid_group = self._pi_tech_group_index(batch)
        group = group[valid_group]
        if int(group.numel()) == 0:
            return
        proven_float = proven.to(torch.float32)
        stable_float = stable.to(torch.float32)
        per_cell_proven_count = proven.sum(dim=1).to(torch.float64)[valid_group]
        per_cell_stable_count = stable.sum(dim=1).to(torch.float64)[valid_group]
        per_cell_proven_sum = (
            pi_thin.to(torch.float32) * proven_float
        ).sum(dim=1).to(torch.float64)[valid_group]
        per_cell_stable_sum = (
            pi_thin.to(torch.float32) * stable_float
        ).sum(dim=1).to(torch.float64)[valid_group]
        for name, value in (
            ("proven_count", per_cell_proven_count),
            ("stable_count", per_cell_stable_count),
            ("proven_sum", per_cell_proven_sum),
            ("stable_sum", per_cell_stable_sum),
        ):
            state[name].index_add_(0, group, value)
        state["thin_cells"].index_add_(
            0, group, torch.ones_like(group, dtype=torch.float64)
        )
        unique_group = torch.unique(group)
        state["thin_group_batches"][unique_group] += 1.0
        state["thin_probe_steps"] += 1.0
        donor = batch["donor_id"].long()[valid_group]
        valid_donor = (
            (donor >= 0) & (donor < int(self._pi_tech_diag_n_donors))
        )
        state["thin_donor_seen"][
            group[valid_donor], donor[valid_donor]
        ] = 1.0

    @torch.no_grad()
    def _finalize_pi_tech_epoch_diagnostics(self) -> Dict[str, float]:
        state = self._pi_tech_epoch_diagnostics
        if state is None:
            return {}
        reduced = {name: value.clone() for name, value in state.items()}
        if dist.is_available() and dist.is_initialized():
            for value in reduced.values():
                dist.all_reduce(
                    value, op=dist.ReduceOp.SUM, group=self._v8_sync_group
                )
        metrics: Dict[str, float] = {}
        cfg = self.pi_tech_config

        def _count(name: str, value: torch.Tensor) -> None:
            metrics[name] = float(value.detach().item())

        proven_count = reduced["proven_count"].sum()
        stable_count = reduced["stable_count"].sum()
        thin_cells = reduced["thin_cells"].sum()
        thin_donors = (reduced["thin_donor_seen"].sum(dim=0) > 0).sum()
        _count("count/pi_tech_proven_positions", proven_count)
        _count("count/pi_tech_stable_positions", stable_count)
        _count("count/pi_tech_thinning_cells", thin_cells)
        _count("count/pi_tech_thinning_probe_steps", reduced["thin_probe_steps"].sum())
        _count("count/pi_tech_thinning_distinct_donors", thin_donors)
        global_supported = (
            float(proven_count) >= int(cfg.diagnostic_min_positions)
            and float(stable_count) >= int(cfg.diagnostic_min_positions)
            and float(thin_cells) >= int(cfg.diagnostic_min_cells)
            and float(thin_donors)
            >= int(cfg.diagnostic_min_distinct_donors)
        )
        if global_supported:
            proven_mean = reduced["proven_sum"].sum() / proven_count
            stable_mean = reduced["stable_sum"].sum() / stable_count
            metrics["metric/pi_tech_at_proven_dropout"] = float(proven_mean)
            metrics["metric/pi_tech_at_stable_zero"] = float(stable_mean)
            metrics["metric/pi_tech_proven_stable_gap"] = float(
                proven_mean - stable_mean
            )

        bank_query_total = reduced["bank_query"].sum()
        if float(bank_query_total) > 0.0:
            metrics["count/pi_tech_bank_queries"] = float(bank_query_total)
            metrics["metric/pi_tech_bank_valid_query_fraction"] = float(
                reduced["bank_valid"].sum() / bank_query_total
            )
            metrics["metric/pi_tech_bank_effective_k"] = float(
                reduced["bank_effective_k"].sum() / bank_query_total
            )
            metrics["count/pi_tech_bank_query_distinct_donors"] = float(
                (reduced["bank_donor_seen"].sum(dim=0) > 0).sum()
            )

        n_regions = int(self._pi_tech_diag_n_regions)
        n_groups = int(reduced["proven_count"].numel())
        for group_index in range(n_groups):
            celltype = group_index // n_regions
            region = group_index % n_regions
            suffix = f"ct{celltype}_region{region}"
            bank_query = reduced["bank_query"][group_index]
            if float(bank_query) > 0.0:
                metrics[
                    f"count/pi_tech_bank_queries_{suffix}"
                ] = float(bank_query)
                metrics[
                    f"metric/pi_tech_bank_valid_query_fraction_{suffix}"
                ] = float(
                    reduced["bank_valid"][group_index] / bank_query
                )
                metrics[
                    f"metric/pi_tech_bank_effective_k_{suffix}"
                ] = float(
                    reduced["bank_effective_k"][group_index] / bank_query
                )
                metrics[
                    f"count/pi_tech_bank_query_distinct_donors_{suffix}"
                ] = float(
                    (
                        reduced["bank_donor_seen"][group_index] > 0
                    ).sum()
                )

            group_proven = reduced["proven_count"][group_index]
            group_stable = reduced["stable_count"][group_index]
            group_cells = reduced["thin_cells"][group_index]
            group_donors = (
                reduced["thin_donor_seen"][group_index] > 0
            ).sum()
            if float(group_cells) <= 0.0:
                continue
            metrics[f"count/pi_tech_proven_positions_{suffix}"] = float(
                group_proven
            )
            metrics[f"count/pi_tech_stable_positions_{suffix}"] = float(
                group_stable
            )
            metrics[f"count/pi_tech_thinning_cells_{suffix}"] = float(
                group_cells
            )
            metrics[
                f"count/pi_tech_thinning_batches_{suffix}"
            ] = float(reduced["thin_group_batches"][group_index])
            metrics[
                f"count/pi_tech_thinning_distinct_donors_{suffix}"
            ] = float(group_donors)
            supported = (
                float(group_proven) >= int(cfg.diagnostic_min_positions)
                and float(group_stable) >= int(cfg.diagnostic_min_positions)
                and float(group_cells) >= int(cfg.diagnostic_min_cells)
                and float(group_donors)
                >= int(cfg.diagnostic_min_distinct_donors)
            )
            if not supported:
                continue
            proven_mean = (
                reduced["proven_sum"][group_index] / group_proven
            )
            stable_mean = (
                reduced["stable_sum"][group_index] / group_stable
            )
            metrics[f"metric/pi_tech_at_proven_dropout_{suffix}"] = float(
                proven_mean
            )
            metrics[f"metric/pi_tech_at_stable_zero_{suffix}"] = float(
                stable_mean
            )
            metrics[f"metric/pi_tech_proven_stable_gap_{suffix}"] = float(
                proven_mean - stable_mean
            )
        return metrics

    def _pi_tech_bank_blend(self) -> float:
        cfg = self.pi_tech_config
        epoch = int(self._pi_tech_epoch)
        start = int(cfg.bank_start_epoch)
        if epoch < start:
            return 0.0
        ramp = int(cfg.bank_ramp_epochs)
        if ramp <= 0:
            return 1.0
        return min(1.0, max(0.0, float(epoch - start + 1) / float(ramp)))

    def _validate_pi_tech_bank_against_train_dataset(
        self, bank: TrainOnlyPiTechBank
    ) -> None:
        """Prove a refreshed/restored bank is exactly the configured train view."""

        dataset = self._pi_tech_bank_dataset
        local_indices = self._pi_tech_bank_local_indices
        if dataset is None or local_indices is None:
            raise RuntimeError(
                "technical-zero bank cannot be validated without train selection"
            )
        train_rows = np.asarray(dataset.row_idx, dtype=np.int64)
        expected_rows = train_rows[np.asarray(local_indices, dtype=np.int64)]
        observed_rows = bank.row.detach().cpu().numpy().astype(np.int64)
        if not np.array_equal(np.sort(observed_rows), np.sort(expected_rows)):
            raise RuntimeError(
                "technical-zero checkpoint bank rows differ from the deterministic "
                "train-only selection"
            )
        if len(set(int(value) for value in observed_rows.tolist())) != len(
            observed_rows
        ):
            raise RuntimeError("technical-zero checkpoint bank rows are duplicated")

        order = np.argsort(observed_rows)
        rows = observed_rows[order]
        observed = {
            "celltype": bank.celltype.detach().cpu().numpy()[order],
            "donor": bank.donor.detach().cpu().numpy()[order],
        }
        expected = {
            "celltype": np.asarray(dataset.celltype_ids)[rows],
            "donor": np.asarray(dataset.donor_ids)[rows],
        }
        if bool(self.pi_tech_config.bank_match_region):
            observed["region"] = bank.region.detach().cpu().numpy()[order]
            expected["region"] = np.asarray(dataset.region_ids)[rows]
        else:
            if not bool((bank.region.detach().cpu() == -1).all()):
                raise RuntimeError(
                    "non-region-matched technical-zero bank must store region=-1"
                )
        dataset_sex = getattr(dataset, "sex_ids", None)
        observed["sex"] = bank.sex.detach().cpu().numpy()[order]
        expected["sex"] = (
            np.asarray(dataset_sex)[rows]
            if dataset_sex is not None
            else np.full(rows.shape, -1, dtype=np.int64)
        )
        mismatched = [
            name
            for name in observed
            if not np.array_equal(
                np.asarray(observed[name], dtype=np.int64),
                np.asarray(expected[name], dtype=np.int64),
            )
        ]
        if mismatched:
            raise RuntimeError(
                "technical-zero checkpoint bank metadata mismatch: "
                + ", ".join(mismatched)
            )
        if self.gene_detect_rate is not None and bank.n_genes != int(
            self.gene_detect_rate.numel()
        ):
            raise RuntimeError(
                "technical-zero checkpoint bank gene count differs from decoder"
            )

    @torch.no_grad()
    def _refresh_pi_tech_bank(self, epoch_index: int) -> None:
        cfg = self.pi_tech_config
        if (
            not self.v8_enabled
            or not cfg.enabled
            or str(cfg.mode) != "same_celltype_donor_bank"
            or int(epoch_index) < int(cfg.bank_start_epoch)
        ):
            return
        refresh_every = max(1, int(cfg.bank_refresh_epochs))
        if (
            self._pi_tech_bank is not None
            and int(epoch_index) != int(cfg.bank_start_epoch)
            and (int(epoch_index) - int(cfg.bank_start_epoch)) % refresh_every != 0
        ):
            return
        dataset = self._pi_tech_bank_dataset
        indices = self._pi_tech_bank_local_indices
        if dataset is None or indices is None or int(indices.size) == 0:
            raise RuntimeError("technical-zero train bank has no selected rows")

        system = self._unwrapped_system()
        was_training = bool(system.training)
        # The training dataset may generate a stochastic binomial-thinning
        # view in __getitem__.  The bank needs only the untouched y_ord, and it
        # must not advance that persistent thinning stream at an epoch boundary.
        thinning_owner = dataset
        while "base_dataset" in vars(thinning_owner):
            thinning_owner = vars(thinning_owner)["base_dataset"]
        original_return_thinned = getattr(
            thinning_owner, "return_thinned", None
        )
        data_generator = torch.Generator()
        data_generator.manual_seed(
            int(cfg.bank_seed) + 1_000_003 * int(epoch_index)
        )
        z_parts = []
        detect_parts = []
        row_parts = []
        celltype_parts = []
        region_parts = []
        donor_parts = []
        sex_parts = []
        try:
            if isinstance(original_return_thinned, bool):
                thinning_owner.return_thinned = False
            loader = DataLoader(
                Subset(dataset, indices.tolist()),
                batch_size=max(1, int(cfg.bank_batch_size)),
                shuffle=False,
                num_workers=max(0, int(cfg.bank_num_workers)),
                pin_memory=self._infer_device().type == "cuda",
                drop_last=False,
                generator=data_generator,
            )
            system.eval()
            for raw_batch in loader:
                bank_batch = self._move_batch_to_device(raw_batch)
                with self._autocast_context():
                    # The bank consumes only z_perp.  Skipping the ordinal
                    # decoder and auxiliary heads avoids a full [B,G] decode
                    # at every epoch boundary and cannot alter the bank value.
                    bank_out = system(
                        bank_batch,
                        sample_latent=False,
                        compute_decoder=False,
                    )
                z_parts.append(bank_out.z_perp.detach().float().cpu())
                detect_parts.append(
                    (bank_batch["y_ord"] > 0).detach().to("cpu", torch.uint8)
                )
                row_parts.append(bank_batch["row_index"].detach().long().cpu())
                celltype_parts.append(
                    bank_batch["celltype_id"].detach().long().cpu()
                )
                region = bank_batch.get("region_id")
                if bool(cfg.bank_match_region) and region is None:
                    raise KeyError(
                        "region-matched technical-zero bank requires region_id"
                    )
                if not bool(cfg.bank_match_region) or region is None:
                    region = torch.full_like(bank_batch["donor_id"], -1)
                region_parts.append(region.detach().long().cpu())
                donor_parts.append(bank_batch["donor_id"].detach().long().cpu())
                sex = bank_batch.get("sex_id")
                if sex is None:
                    sex = torch.full_like(bank_batch["donor_id"], -1)
                sex_parts.append(sex.detach().long().cpu())
        finally:
            system.train(was_training)
            if isinstance(original_return_thinned, bool):
                thinning_owner.return_thinned = original_return_thinned

        bank_cpu = TrainOnlyPiTechBank(
            z=torch.cat(z_parts, dim=0),
            detect=torch.cat(detect_parts, dim=0),
            row=torch.cat(row_parts, dim=0),
            celltype=torch.cat(celltype_parts, dim=0),
            region=torch.cat(region_parts, dim=0),
            donor=torch.cat(donor_parts, dim=0),
            sex=torch.cat(sex_parts, dim=0),
            epoch=int(epoch_index),
        )
        self._validate_pi_tech_bank_against_train_dataset(bank_cpu)
        if self.gene_detect_rate is not None and (
            bank_cpu.n_genes != int(self.gene_detect_rate.numel())
        ):
            raise RuntimeError(
                "technical-zero bank gene count does not match decoder: "
                f"{bank_cpu.n_genes} vs {int(self.gene_detect_rate.numel())}"
            )
        self._pi_tech_bank = bank_cpu.to(self._infer_device())
        self._pi_tech_bank._build_query_statistics(
            match_region=bool(cfg.bank_match_region),
            sex_linked_gene_indices=tuple(cfg.sex_linked_gene_indices),
        )
        if self._is_rank0():
            print(
                "[pi-tech-bank] "
                f"epoch={int(epoch_index)} cells={int(bank_cpu.z.shape[0])} "
                f"celltypes={int(torch.unique(bank_cpu.celltype).numel())} "
                f"regions={int(torch.unique(bank_cpu.region).numel())} "
                f"donors={int(torch.unique(bank_cpu.donor).numel())}",
                flush=True,
            )

    @torch.no_grad()
    def _compute_pi_tech(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
        y_ord: torch.Tensor,
        gene_detect_rate: torch.Tensor,
        *,
        accumulate_epoch: bool = False,
    ) -> tuple[torch.Tensor, Dict[str, float]]:
        cfg = self.pi_tech_config
        depth = depth_proxy_from_batch(batch, y_ord, cfg)
        in_batch = compute_pi_tech_in_batch(
            model_out.z_perp, y_ord, depth, gene_detect_rate, cfg
        )
        sex_indices = tuple(
            int(index)
            for index in cfg.sex_linked_gene_indices
            if 0 <= int(index) < int(y_ord.shape[1])
        )
        diagnostics: Dict[str, float] = {
            "metric/pi_tech_bank_blend": 0.0,
        }
        blend = self._pi_tech_bank_blend()
        if (
            str(cfg.mode) != "same_celltype_donor_bank"
            or self._pi_tech_bank is None
            or blend <= 0.0
        ):
            return in_batch, diagnostics
        required = ["celltype_id", "donor_id", "row_index"]
        if bool(cfg.bank_match_region):
            required.append("region_id")
        missing = [key for key in required if key not in batch]
        if missing:
            raise KeyError(
                "same_celltype_donor_bank batch is missing " + ", ".join(missing)
            )
        sex = batch.get("sex_id")
        if sex is None:
            sex = torch.full_like(batch["donor_id"], -1)
        bank_pi, bank_diag = compute_pi_tech_from_bank(
            model_out.z_perp,
            y_ord,
            depth,
            gene_detect_rate,
            batch["celltype_id"],
            (
                batch["region_id"]
                if bool(cfg.bank_match_region)
                else torch.full_like(batch["celltype_id"], -1)
            ),
            batch["donor_id"],
            batch["row_index"],
            sex,
            self._pi_tech_bank,
            cfg,
        )
        valid = bank_diag["valid_query"].bool()
        invalid_fallback = (
            torch.zeros_like(in_batch)
            if str(cfg.bank_invalid_fallback) == "zero"
            else in_batch
        )
        safe_bank = torch.where(valid[:, None], bank_pi, invalid_fallback)
        if sex_indices:
            gene_index = torch.tensor(
                sex_indices, device=safe_bank.device, dtype=torch.long
            )
            sex_safe = bank_diag["sex_valid_query"].bool() & (sex >= 0)
            safe_bank[:, gene_index] = torch.where(
                sex_safe[:, None],
                bank_pi[:, gene_index],
                torch.zeros_like(bank_pi[:, gene_index]),
            )
        pi = in_batch.lerp(safe_bank, float(blend)).detach()
        effective_k = bank_diag["effective_k"].float()
        diagnostics.update(
            {
                "metric/pi_tech_bank_blend": float(blend),
                "metric/pi_tech_bank_valid_query_fraction": float(
                    valid.float().mean().detach()
                ),
                "metric/pi_tech_bank_effective_k": float(
                    effective_k.mean().detach()
                ),
                "metric/pi_tech_bank_epoch": float(self._pi_tech_bank.epoch),
                "metric/pi_tech_bank_invalid_zero_fraction": (
                    float((~valid).float().mean().detach())
                    if str(cfg.bank_invalid_fallback) == "zero"
                    else 0.0
                ),
            }
        )
        if tuple(cfg.sex_linked_gene_indices):
            diagnostics["metric/pi_tech_bank_sex_valid_query_fraction"] = float(
                bank_diag["sex_valid_query"].float().mean().detach()
            )
        if bool(accumulate_epoch):
            self._accumulate_pi_tech_bank_diagnostics(
                batch, valid, effective_k
            )
        return pi, diagnostics

    # ------------------------------------------------------------------ #
    # Epoch hooks
    # ------------------------------------------------------------------ #
    def train_epoch(self, loader, *args, epoch_index: Optional[int] = None, **kwargs):
        self._pi_tech_epoch = int(epoch_index or 1)
        if self.v8_enabled:
            self._reset_pi_tech_epoch_diagnostics()
            self._refresh_pi_tech_bank(self._pi_tech_epoch)
        if self.v8_enabled and self.cvar is not None and self.cvar.cfg.enabled:
            self.cvar.begin_epoch(epoch_index)
        out = super().train_epoch(loader, *args, epoch_index=epoch_index, **kwargs)
        if self.v8_enabled:
            if self.cvar is not None and self.cvar.cfg.enabled:
                self.cvar.end_epoch_sync_and_update(epoch_index)
                if self._is_rank0():
                    pool = self.cvar.ema_count > 0
                    worst = float(self.cvar.ema_metric[pool].max().item()) if bool(pool.any()) else 0.0
                    print(
                        f"[v8-cvar] epoch={epoch_index} "
                        f"triggered={int(self.cvar.triggered.sum().item())}/{self.cvar.n_groups} "
                        f"active={self.cvar.is_active()} worst_group_metric={worst:.4f}",
                        flush=True,
                    )
            self._sync_gene_detect_rate()
            out.metrics.update(self._finalize_pi_tech_epoch_diagnostics())
        return out

    # ------------------------------------------------------------------ #
    # Loss composition
    # ------------------------------------------------------------------ #
    def _build_tech_augmenter(self, device):
        """Lazily build the tech-invariance augmenter from the train dataset's reference scaler
        (held by run_current.V7A_TRAIN_DATASET). disease_dir (optional [C,G]) is read from
        self._disease_dir_tensor if precomputed; otherwise the disease negative control is skipped."""
        try:
            sm = getattr(self, "_tech_scaler_mean", None)
            ss = getattr(self, "_tech_scaler_std_eps", None)
            clip = getattr(self, "_tech_scaler_clip", None)
            if sm is None or ss is None:  # fallback: the run_current global
                from kmlee_bam.training import run_current as _rc
                ds = getattr(_rc, "V7A_TRAIN_DATASET", None)
                sm = getattr(ds, "_scaler_mean", None)
                ss = getattr(ds, "_scaler_std_eps", None)
                clip = getattr(ds, "x_gene_scalar_clip", None)
            if sm is None or ss is None:
                print("[tech-inv] no reference-scaler cache (trainer + global both empty) -> tech-invariance OFF", flush=True)
                return None
            tic = self.tech_inv_config
            aug = TechInvarianceAugmenter(
                torch.as_tensor(sm), torch.as_tensor(ss), clip,
                tech_kinds=list(tic.tech_kinds), dropout_p=tic.dropout_p, gain_lo=tic.gain_lo,
                gain_hi=tic.gain_hi, ambient_a=tic.ambient_a, disease_dir=self._disease_dir_tensor,
            ).to(device)
            print(f"[tech-inv] augmenter ready (kinds={tic.tech_kinds}, "
                  f"disease_ctrl={'on' if self._disease_dir_tensor is not None else 'off'})", flush=True)
            return aug
        except Exception as e:  # never let augmentation crash training
            print(f"[tech-inv] augmenter build failed: {e!r} -> tech-invariance OFF", flush=True)
            return None

    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        if not self.v8_enabled:
            return super()._compute_loss(batch, model_out)

        is_training = bool(self.system.training and torch.is_grad_enabled())
        y_ord = batch["y_ord"].long()

        # ---- π_tech: compute + stash BEFORE super() so V4's L1 picks it up ----
        self._pending_pi_tech = None
        pi_tech_local = None
        pi_tech_diagnostics: Dict[str, float] = {}
        if self.pi_tech_config.enabled:
            if is_training:
                self._update_gene_detect_rate(y_ord)
            gdr = self.gene_detect_rate
            if gdr is None:
                gdr = (y_ord > 0).to(dtype=torch.float32).mean(dim=0)
            self._pending_pi_tech, pi_tech_diagnostics = self._compute_pi_tech(
                batch,
                model_out,
                y_ord,
                gdr,
                accumulate_epoch=is_training,
            )
            pi_tech_local = self._pending_pi_tech  # keep after the conduit is cleared

        # ---- Stage D: per-gene reconstruction reliability. Downweight suspicious
        #      observed zeros (high π_tech) in the MAIN reconstruction so technical
        #      zeros do not pull the latent toward a hard biological off. Stash
        #      BEFORE super() so core_trainer's criterion picks it up. Off ⇒ None
        #      ⇒ reconstruction byte-identical.
        self._pending_rec_weight_per_gene = None
        rec_w = None
        if self.zero_reliability_config.enabled and pi_tech_local is not None:
            rec_w = zero_reliability_weights(
                y_ord,
                pi_tech_local,
                alpha=float(self.zero_reliability_config.alpha),
                floor=float(self.zero_reliability_config.floor),
            )
            self._pending_rec_weight_per_gene = rec_w

        out = super()._compute_loss(batch, model_out)
        self._pending_pi_tech = None  # clear the conduit
        self._pending_rec_weight_per_gene = None  # clear the conduit

        if is_training:
            self._v8_step_counter += 1

        total = out.total
        details = dict(out.details)
        details.update(pi_tech_diagnostics)
        # DDP-stable: keep tech_risk_head params in the autograd graph EVERY step via a
        # zero-scaled touch, so they always receive a (zero) grad even when the tech-inv block
        # is skipped this step. Prevents the "unused parameter" reducer crash.
        _trl = getattr(getattr(model_out, "encoder_out", None), "tech_risk_logit", None)
        if _trl is not None:
            total = total + 0.0 * _trl.sum()
        if rec_w is not None:
            details["metric/zero_reliability_mean"] = float(rec_w.mean().detach())
            details["metric/zero_reliability_min"] = float(rec_w.min().detach())
            details["metric/zero_reliability_zero_frac_affected"] = float(
                (rec_w < 1.0).to(torch.float32).mean().detach()
            )

        # ---- thinning consistency (+ optional supervision), training only ----
        if (
            self.thinning_config.enabled
            and is_training
            and (self._v8_step_counter % max(1, int(self.thinning_config.every_n_steps)) == 0)
        ):
            thinned = build_thinned_batch(batch)
            if thinned is not None:
                with self._autocast_context():
                    thin_out = self.system(thinned, sample_latent=False)
                cons = thinning_consistency_loss(
                    model_out.z_perp,
                    thin_out.z_perp,
                    metric=str(self.thinning_config.consistency_metric),
                    detach_full=bool(self.thinning_config.detach_full),
                )
                total = total + float(self.thinning_config.lambda_consistency) * cons
                details["loss/thin_consistency"] = float(cons.detach())
                # Stage E: validate the π_tech "judge" against proven labels.
                # IMPORTANT: π_tech is forced 0 at observed-NONZERO positions, and a
                # proven dropout is NONZERO in the FULL view — so full-view π_tech is
                # ~0 there (useless). Recompute π_tech on the THINNED view (where the
                # proven dropout appears as a zero) and compare proven dropouts vs
                # zeros that stayed zero. A working judge ⇒ pi@proven > pi@stable
                # (positive gap). Diagnostic only (no gradient); runs every_n_steps.
                if self.pi_tech_config.enabled and "y_ord_thin" in batch:
                    yt = batch["y_ord_thin"].long()
                    gdr = self.gene_detect_rate
                    if gdr is None:
                        gdr = (yt > 0).to(torch.float32).mean(dim=0)
                    with torch.no_grad():
                        pi_thin, _ = self._compute_pi_tech(
                            thinned, thin_out, yt, gdr
                        )
                    proven = (y_ord > 0) & (yt == 0)
                    stable0 = (y_ord == 0) & (yt == 0)
                    self._accumulate_pi_tech_thinning_diagnostics(
                        batch, pi_thin, proven, stable0
                    )
                if float(self.thinning_config.lambda_tech_zero_sup) > 0.0 and "y_ord_thin" in batch:
                    mask = technical_zero_mask(y_ord, batch["y_ord_thin"].long())
                    sup = technical_zero_supervision_loss(thin_out.decoder_out.probs, mask)
                    total = total + float(self.thinning_config.lambda_tech_zero_sup) * sup
                    details["loss/thin_tech_zero_sup"] = float(sup.detach())
                    details["metric/thin_tech_zero_frac"] = float(mask.to(torch.float32).mean().detach())

        # ---- tech-invariance: BAM as supervised tech-OOD detector + z robustness ----
        #  Build a technically-corrupted VIEW of the batch (dropout/gain/ambient) and require:
        #   (1) z invariant to it, (2) the risk head flags it (clean=0/tech=1), (3) the prediction
        #   stays consistent with the CLEAN prediction. Optional disease-view negative control keeps
        #   the risk head from firing on real biology. CORE = (1)+(2)+(3); risk-weighted recon deferred.
        tic = self.tech_inv_config
        if (
            tic.enabled
            and is_training
            and (self._v8_step_counter % max(1, int(tic.every_n_steps)) == 0)
        ):
            if self.tech_augmenter is None:
                self.tech_augmenter = self._build_tech_augmenter(model_out.z_perp.device)
            aug = self.tech_augmenter
            if aug is not None:
                ramp = 1.0 if tic.warmup_steps <= 0 else min(1.0, self._v8_step_counter / float(tic.warmup_steps))
                tech_view = aug.make_tech_view(batch)
                with self._autocast_context():
                    tech_out = self.system(tech_view, sample_latent=False)
                rc = getattr(model_out.encoder_out, "tech_risk_logit", None)
                rt = getattr(tech_out.encoder_out, "tech_risk_logit", None)
                rd = None
                if tic.lambda_disease_ctrl > 0.0 and self._disease_dir_tensor is not None:
                    dis_view = aug.make_disease_view(batch, gamma=tic.disease_gamma)
                    if dis_view is not None:
                        with self._autocast_context():
                            dis_out = self.system(dis_view, sample_latent=False)
                        rd = getattr(dis_out.encoder_out, "tech_risk_logit", None)
                # (1) z-consistency: tech corruption must not move the biological state ESTIMATE.
                #     Compare the DETERMINISTIC posterior mean mu_q (NOT the sampled z_perp) so that
                #     (a) posterior sampling noise stays out of the target, (b) we don't accidentally
                #     shrink the posterior variance, and (c) pooled_state stays tech-aware for the
                #     risk head while mu_q is the invariant (resolves the invariant-vs-sensitive tension).
                if tic.lambda_z_cons > 0.0:
                    zc = z_consistency_loss(
                        model_out.encoder_out.mu_q, tech_out.encoder_out.mu_q,
                        detach_clean=tic.detach_clean,
                    )
                    total = total + ramp * float(tic.lambda_z_cons) * zc
                    details["loss/tech_z_cons"] = float(zc.detach())
                # (2) supervised tech-risk detector (+ disease negative control)
                if tic.lambda_bam_tech > 0.0 and rc is not None and rt is not None:
                    bce = tech_risk_bce(rc, rt, rd, disease_weight=float(tic.lambda_disease_ctrl))
                    total = total + ramp * float(tic.lambda_bam_tech) * bce
                    details["loss/tech_risk_bce"] = float(bce.detach())
                    details["metric/risk_clean"] = float(torch.sigmoid(rc).mean().detach())
                    details["metric/risk_tech"] = float(torch.sigmoid(rt).mean().detach())
                    if rd is not None:
                        details["metric/risk_disease"] = float(torch.sigmoid(rd).mean().detach())
                elif tic.lambda_bam_tech > 0.0 and rc is None:
                    details["metric/tech_risk_head_missing"] = 1.0  # encoder.tech_risk_head not enabled
                # (3) prediction consistency to the CLEAN prediction (never the corrupted target)
                if tic.lambda_rec_cons > 0.0:
                    rcl = rec_consistency_loss(
                        model_out.decoder_out.probs, tech_out.decoder_out.probs,
                        n_bins=model_out.decoder_out.probs.shape[-1], metric=str(tic.rec_metric),
                    )
                    total = total + ramp * float(tic.lambda_rec_cons) * rcl
                    details["loss/tech_rec_cons"] = float(rcl.detach())

        # ---- adaptive subgroup CVaR: monitor + (maybe) intervene ----
        if self.cvar is not None and self.cvar.cfg.enabled:
            depth_c = depth_proxy_from_batch(batch, y_ord, self.pi_tech_config)
            gid = self.cvar.group_ids(batch["celltype_id"], _resolve_tech_id(batch), depth_c)
            metric = per_cell_metric(
                model_out.decoder_out.nll_per_cell,
                model_out.decoder_out.probs,
                y_ord,
                metric=str(self.cvar.cfg.metric),
            )
            if is_training:
                self.cvar.observe_batch(gid, metric)

            if self.cvar.is_active():
                w_cell = self.cvar.weights_for(gid)                      # mean ≈ 1, no-grad
                nll = model_out.decoder_out.nll_per_cell
                # Reweight-only delta. Gradient is mean((w-1)·∇nll): exactly 0
                # when w_cell ≡ 1, and it never re-adds the base reconstruction
                # gradient (base `total` already contains nll). The earlier form
                # `(w*nll).mean() - nll.mean().detach()` had VALUE 0 at w≡1 but a
                # non-zero gradient equal to a full extra recon push of magnitude
                # lambda_cvar·∇nll.mean() throughout the active phase. w_cell is
                # detached (weights_for is @torch.no_grad) so only nll backprops.
                delta = ((w_cell.detach() - 1.0) * nll).mean()
                total = total + float(self.cvar.cfg.lambda_cvar) * delta
                details["metric/cvar_w_mean"] = float(w_cell.mean().detach())
                details["metric/cvar_w_max"] = float(w_cell.max().detach())

            # always-on monitoring (visible from epoch 1)
            details["metric/cvar_n_triggered"] = float(self.cvar.triggered.sum().item())
            details["metric/cvar_active"] = float(self.cvar.is_active())
            pool = self.cvar.ema_count > 0
            details["metric/cvar_worst_group_metric"] = (
                float(self.cvar.ema_metric[pool].max().item()) if bool(pool.any()) else 0.0
            )

        # ---- within-celltype variance floor (Phase 2), training only ----
        # Penalise predicted-tier variance below the in-batch true-tier variance
        # within each celltype (raise per-cell dynamic range R), plus an active
        # depth-decorrelation guardrail so the new spread is biological. Off ⇒
        # block skipped ⇒ byte-identical.
        if self.variance_floor_config.enabled and is_training:
            vf_target = None
            vf_count = None
            if self._var_running is not None:
                self._var_running.update(y_ord, batch["celltype_id"])
                vf_target, vf_count = self._var_running.variance()
            log_depth_vf = torch.log(
                depth_proxy_from_batch(batch, y_ord, self.pi_tech_config).clamp_min(1.0)
            )
            vf = compute_variance_floor(
                model_out.decoder_out.probs,
                y_ord,
                model_out.decoder_out.state_score,
                batch["celltype_id"],
                log_depth_vf,
                self.variance_floor_config,
                var_true_target=vf_target,
                target_count=vf_count,
            )
            total = (
                total
                + float(self.variance_floor_config.lambda_var) * vf.loss_var
                + float(self.variance_floor_config.lambda_depth_decorr) * vf.loss_depth_decorr
            )
            details["loss/var_floor"] = float(vf.loss_var.detach())
            details["loss/var_depth_decorr"] = float(vf.loss_depth_decorr.detach())
            details["metric/var_pred_mean"] = float(vf.var_pred_mean)
            details["metric/var_true_mean"] = float(vf.var_true_mean)
            details["metric/var_gap"] = float(vf.var_gap)
            details["metric/var_floor_active_frac"] = float(vf.active_frac)
            details["metric/var_state_depth_corr"] = float(vf.state_depth_corr)
            details["metric/var_floor_n_groups"] = float(vf.n_groups)

        # ---- high-tier margin ranking (disease-aligned), training only ----
        # Push true-tier4 cells' predicted E[tier] above a DETACHED running low
        # reference (tier1/2), per celltype×gene. Gradient flows only to the high
        # cells (low is not pushed down -> protects the zero-leak goal). Off ⇒
        # block skipped ⇒ byte-identical.
        if self.high_margin_config.enabled and is_training and self._low_ref is not None:
            probs = model_out.decoder_out.probs
            levels = torch.arange(probs.shape[-1], device=probs.device, dtype=torch.float32)
            e_tier = (probs.float() * levels).sum(dim=-1)            # [B,G] E[tier] (grad)
            low_mask = torch.zeros_like(y_ord, dtype=torch.float32)
            for lt in self.high_margin_config.low_tiers:
                low_mask = low_mask + (y_ord == int(lt)).to(torch.float32)
            self._low_ref.update(e_tier.detach(), low_mask, batch["celltype_id"])
            low_ref, low_upd = self._low_ref.reference()
            hm = high_margin_loss(
                e_tier, y_ord, batch["celltype_id"], low_ref, low_upd, probs,
                self.high_margin_config, pi_tech=pi_tech_local,
            )
            total = (
                total
                + float(self.high_margin_config.lambda_high_margin) * hm.loss
                + float(self.high_margin_config.lambda_t4_floor) * hm.loss_t4_floor
                + float(self.high_margin_config.lambda_tier4_focal) * hm.loss_tier4_focal
                + float(self.high_margin_config.lambda_zero_focal) * hm.loss_zero_focal
            )
            details["loss/high_margin"] = float(hm.loss.detach())
            details["loss/t4_floor"] = float(hm.loss_t4_floor.detach())
            details["loss/tier4_focal"] = float(hm.loss_tier4_focal.detach())
            details["loss/zero_focal"] = float(hm.loss_zero_focal.detach())
            details["metric/tier4_fp_on_zero"] = float(hm.t4_fp_on_zero)
            details["metric/high_margin_gap"] = float(hm.margin_gap)
            details["metric/tier4_pred_level"] = float(hm.tier4_pred_level)
            details["metric/low_pred_level"] = float(hm.low_pred_level)
            details["metric/tier4_recall"] = float(hm.tier4_recall)
            details["metric/high_margin_active_frac"] = float(hm.active_frac)
            details["metric/high_margin_n_groups"] = float(hm.n_groups)
            details["metric/high_margin_active_pairs"] = float(hm.active_pairs)
            details["metric/high_margin_active_high_level"] = float(hm.active_high_level)
            details["metric/tier4_fp_on_low"] = float(hm.t4_fp_on_low)

        # ---- BAM/decoder confidence diagnostic (calibrated readout, 2026-06-04) ----
        # The usable per-cell confidence is the DECODER predictive uncertainty, NOT the BAM
        # attention entropy (which is anti-calibrated, partial nzmis|depth -0.26). pe_pnz =
        # P(nonzero)-weighted predictive entropy (deployable; ranks error, partial nzmis|depth
        # ~0.66). Detached, training-neutral. Also log the (weak) BAM attn entropy for contrast.
        with torch.no_grad():
            _p = model_out.decoder_out.probs.float()
            _ent = -(_p * torch.log(_p.clamp_min(1e-8))).sum(dim=-1)             # [B,G]
            _wnz = (1.0 - _p[..., 0]).clamp(0.0, 1.0)                            # P(nonzero)
            _pe_pnz = (_ent * _wnz).sum(dim=1) / _wnz.sum(dim=1).clamp_min(1e-6)  # [B]
            details["metric/pred_entropy_pnz"] = float(_pe_pnz.mean())
            _cu = getattr(getattr(model_out, "encoder_out", None), "cell_uncertainty", None)
            if _cu is not None:
                details["metric/bam_attn_entropy"] = float(_cu.detach().float().mean())

        details["loss/total"] = float(total.detach())
        out.total = total
        out.details = details
        return out

    # ------------------------------------------------------------------ #
    # State dict (checkpoint persistence)
    # ------------------------------------------------------------------ #
    def v8_state_dict(self) -> dict:
        if not self.v8_enabled:
            return {}
        return {
            "cvar": self.cvar.state_dict() if self.cvar is not None else None,
            "gene_detect_rate": (
                self.gene_detect_rate.detach().cpu() if self.gene_detect_rate is not None else None
            ),
            "step_counter": int(self._v8_step_counter),
            "pi_tech_epoch": int(self._pi_tech_epoch),
            "pi_tech_bank": (
                self._pi_tech_bank.state_dict()
                if self._pi_tech_bank is not None
                else None
            ),
        }

    def v8_load_state_dict(self, state: dict) -> None:
        if not self.v8_enabled or not state:
            return
        if state.get("cvar") is not None and self.cvar is not None:
            self.cvar.load_state_dict(state["cvar"])
        if state.get("gene_detect_rate") is not None and self.gene_detect_rate is not None:
            self.gene_detect_rate = state["gene_detect_rate"].to(self._infer_device())
        self._v8_step_counter = int(state.get("step_counter", 0))
        self._pi_tech_epoch = int(state.get("pi_tech_epoch", 0))
        if (
            state.get("pi_tech_bank") is not None
            and str(self.pi_tech_config.mode) == "same_celltype_donor_bank"
        ):
            self._pi_tech_bank = TrainOnlyPiTechBank.from_state_dict(
                state["pi_tech_bank"], device=self._infer_device()
            )
            self._validate_pi_tech_bank_against_train_dataset(
                self._pi_tech_bank
            )
            self._pi_tech_bank._build_query_statistics(
                match_region=bool(self.pi_tech_config.bank_match_region),
                sex_linked_gene_indices=tuple(
                    self.pi_tech_config.sex_linked_gene_indices
                ),
            )
