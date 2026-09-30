"""Runtime hard birth/death search for Lie decoder generators.

This module is intentionally separate from the ordinary training objective.
It never adds an L1/L0 term and never learns a soft architecture gate.  At
prespecified epoch boundaries it:

1. builds one of two deterministic, donor-disjoint validation panels;
2. caches the encoder/PRISM outputs with all model weights frozen;
3. ranks currently active, unprotected generators by the first-order harm of
   removing their *joint* Lie-linear plus affine-translation contribution;
4. evaluates the proposed hard block removal with exact decoder forwards;
5. asks :class:`HardGeneratorBudgetController` to enforce paired
   non-inferiority and two-panel hysteresis.

The full view retains Direct and PRISM.  The generator-isolated view removes
both, so a proposal cannot be accepted merely because an additive bypass has
already copied the generator's work.  Search metrics are computed on validation
only; the held-out test split is never opened here.
"""

from __future__ import annotations

import json
import math
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

from kmlee_bam.training.generator_budget_controller import (
    HardGeneratorBudgetController,
    MetricRule,
)


@dataclass
class GeneratorBudgetSearchConfig:
    enabled: bool = False
    start_epoch: int = 4
    interval_epochs: int = 1
    block_size: int = 16
    minimum_active_generators: int = 128
    protect_singletons: bool = True
    required_confirmations: int = 2
    cooldown_checks: int = 1
    cells_per_donor_celltype: int = 4
    minimum_donors_per_celltype: int = 3
    panel_count: int = 2
    seed: int = 420729
    eval_batch_size: int = 32
    full_module_noninferiority: float = 0.015
    isolated_module_noninferiority: float = 0.020
    full_nll_noninferiority: float = 0.010
    isolated_nll_noninferiority: float = 0.020
    rank_isolated_weight: float = 1.0
    max_cached_cells: int = 4096

    def validate(self) -> None:
        integer_positive = (
            "start_epoch",
            "interval_epochs",
            "block_size",
            "minimum_active_generators",
            "required_confirmations",
            "cells_per_donor_celltype",
            "minimum_donors_per_celltype",
            "panel_count",
            "eval_batch_size",
            "max_cached_cells",
        )
        for name in integer_positive:
            value = getattr(self, name)
            if isinstance(value, bool) or int(value) != value or int(value) <= 0:
                raise ValueError(f"generator_budget_search.{name} must be positive.")
        if int(self.panel_count) < int(self.required_confirmations):
            raise ValueError(
                "generator_budget_search.panel_count must be at least "
                "required_confirmations."
            )
        if int(self.cooldown_checks) < 0:
            raise ValueError(
                "generator_budget_search.cooldown_checks must be non-negative."
            )
        for name in (
            "full_module_noninferiority",
            "isolated_module_noninferiority",
            "full_nll_noninferiority",
            "isolated_nll_noninferiority",
            "rank_isolated_weight",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(
                    f"generator_budget_search.{name} must be finite and non-negative."
                )


@dataclass
class _CachedDecoderBatch:
    z: torch.Tensor
    celltype: torch.Tensor
    tech: torch.Tensor
    sex: Optional[torch.Tensor]
    precision: Optional[torch.Tensor]
    target: torch.Tensor
    donor: torch.Tensor


@dataclass
class _PanelCache:
    batches: list[_CachedDecoderBatch]
    target_module: np.ndarray
    group_celltype: np.ndarray
    module_keep: np.ndarray
    n_cells: int
    n_groups: int
    n_celltypes: int


def _unwrap(module: nn.Module) -> nn.Module:
    if hasattr(module, "module"):
        return module.module
    if hasattr(module, "_ddp_model"):
        ddp = getattr(module, "_ddp_model")
        if hasattr(ddp, "module"):
            return ddp.module
    return module


@contextmanager
def _frozen_parameters(module: nn.Module):
    parameters = tuple(module.parameters())
    flags = tuple(parameter.requires_grad for parameter in parameters)
    try:
        for parameter in parameters:
            parameter.requires_grad_(False)
        yield
    finally:
        for parameter, flag in zip(parameters, flags):
            parameter.requires_grad_(flag)


@contextmanager
def _preserve_training_mode(module: nn.Module):
    was_training = bool(module.training)
    module.eval()
    try:
        yield
    finally:
        module.train(was_training)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Stable average ranks without a scipy dependency."""

    flat = np.asarray(values, dtype=np.float64).reshape(-1)
    order = np.argsort(flat, kind="mergesort")
    sorted_values = flat[order]
    ranks = np.empty(len(flat), dtype=np.float64)
    start = 0
    while start < len(flat):
        end = start + 1
        while end < len(flat) and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1.0
        start = end
    return ranks


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3 or y.size != x.size:
        return float("nan")
    rx = _average_ranks(x)
    ry = _average_ranks(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not math.isfinite(denominator) or denominator <= 1e-12:
        return float("nan")
    return float(np.dot(rx, ry) / denominator)


class GeneratorBudgetSearch:
    """Rank-0 search driver; the caller broadcasts the resulting hard mask."""

    def __init__(
        self,
        *,
        config: GeneratorBudgetSearchConfig,
        system: nn.Module,
        validation_dataset: Any,
        device: torch.device,
        output_dir: str | Path,
    ) -> None:
        config.validate()
        self.config = config
        self.system = _unwrap(system)
        self.decoder = self.system.decoder
        self.validation_dataset = validation_dataset
        self.device = device
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.event_path = self.output_dir / "generator_budget_events.jsonl"
        self.state_path = self.output_dir / "generator_budget_state.json"
        self.manifest_path = self.output_dir / "generator_budget_manifest.json"
        self.check_index = 0
        self.pending_block: Optional[tuple[int, ...]] = None
        self.panel_baselines: dict[int, dict[str, float]] = {}

        if not hasattr(self.decoder, "generator_active_mask"):
            raise RuntimeError(
                "decoder lacks generator_active_mask; hard search cannot run."
            )
        n_generators = int(self.decoder.n_generators)
        if int(config.minimum_active_generators) >= n_generators:
            raise ValueError(
                "minimum_active_generators must be smaller than the candidate bank."
            )
        module_sizes = (
            (self.decoder.module_masks.detach().cpu() > 0).sum(dim=1).numpy()
        )
        protected = (
            tuple(int(i) for i in np.flatnonzero(module_sizes == 1))
            if bool(config.protect_singletons)
            else ()
        )
        fresh_controller = HardGeneratorBudgetController(
            num_generators=n_generators,
            metric_rules=(
                MetricRule(
                    "module_recovery",
                    "maximize",
                    noninferiority_margin=float(
                        config.full_module_noninferiority
                    ),
                ),
                MetricRule(
                    "isolated_module_recovery",
                    "maximize",
                    noninferiority_margin=float(
                        config.isolated_module_noninferiority
                    ),
                ),
                MetricRule(
                    "full_nll",
                    "minimize",
                    noninferiority_margin=float(config.full_nll_noninferiority),
                ),
                MetricRule(
                    "isolated_nll",
                    "minimize",
                    noninferiority_margin=float(
                        config.isolated_nll_noninferiority
                    ),
                ),
            ),
            protected_ids=protected,
            minimum_active_generators=int(config.minimum_active_generators),
            maximum_active_generators=n_generators,
            required_confirmations=int(config.required_confirmations),
            cooldown_checks=int(config.cooldown_checks),
        )
        if self.state_path.exists():
            # runner_base currently supports fresh warm starts, not an exact
            # optimizer/scheduler/epoch resume.  Restoring only this controller
            # onto the original warm-start weights would bind old architecture
            # decisions to the wrong model state.  Refuse that unsafe partial
            # resume until an atomic whole-training resume contract exists.
            raise RuntimeError(
                "generator-budget state already exists in the output "
                "directory. Partial controller-only resume is forbidden; "
                "use a fresh output directory."
            )
        self.controller = fresh_controller
        self.decoder.set_generator_active_mask(
            torch.tensor(
                self.controller.active_mask,
                dtype=torch.bool,
                device=self.device,
            )
        )
        self._write_manifest()

    @property
    def active_count(self) -> int:
        return len(self.controller.active_ids)

    def due(self, epoch: int) -> bool:
        return bool(
            self.config.enabled
            and int(epoch) >= int(self.config.start_epoch)
            and (
                int(epoch) - int(self.config.start_epoch)
            )
            % int(self.config.interval_epochs)
            == 0
            and self.active_count > int(self.config.minimum_active_generators)
        )

    def _panel_local_indices(self, panel_index: int) -> np.ndarray:
        dataset = self.validation_dataset
        row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
        donor = np.asarray(dataset.donor_ids[row_idx], dtype=np.int64)
        celltype = np.asarray(dataset.celltype_ids[row_idx], dtype=np.int64)
        donor_values = np.unique(donor)
        rng = np.random.default_rng(int(self.config.seed))
        donor_order = donor_values.copy()
        rng.shuffle(donor_order)
        panel_donors = donor_order[
            int(panel_index) :: int(self.config.panel_count)
        ]
        panel_donor_set = set(int(value) for value in panel_donors.tolist())

        selected: list[int] = []
        cells_per_group = int(self.config.cells_per_donor_celltype)
        for ct in sorted(int(value) for value in np.unique(celltype).tolist()):
            for dn in sorted(panel_donor_set):
                candidates = np.flatnonzero((celltype == ct) & (donor == dn))
                if len(candidates) < cells_per_group:
                    continue
                group_seed = (
                    int(self.config.seed)
                    + 104729 * int(panel_index)
                    + 1009 * int(ct)
                    + 9176 * int(dn)
                )
                group_rng = np.random.default_rng(group_seed)
                picked = group_rng.choice(
                    candidates,
                    size=cells_per_group,
                    replace=False,
                )
                selected.extend(int(value) for value in picked.tolist())
        selected = sorted(set(selected))
        if len(selected) > int(self.config.max_cached_cells):
            cap_rng = np.random.default_rng(
                int(self.config.seed) + 65537 * int(panel_index)
            )
            selected = sorted(
                int(value)
                for value in cap_rng.choice(
                    np.asarray(selected, dtype=np.int64),
                    size=int(self.config.max_cached_cells),
                    replace=False,
                ).tolist()
            )
        if not selected:
            raise RuntimeError(
                f"generator search panel {panel_index} has no eligible cells."
            )
        return np.asarray(selected, dtype=np.int64)

    def _module_projection(
        self,
        *,
        gene_sum: torch.Tensor,
        count: torch.Tensor,
        group_celltype: torch.Tensor,
    ) -> tuple[np.ndarray, np.ndarray]:
        valid = count > 0
        means = gene_sum[valid] / count[valid, None].to(gene_sum.dtype)
        weight = self.system.module_tokenizer.activity_weight.to(
            device=means.device,
            dtype=means.dtype,
        )
        module = means @ weight.transpose(0, 1)
        return (
            module.float().cpu().numpy(),
            group_celltype[valid].cpu().numpy(),
        )

    def _build_panel_cache(self, panel_index: int) -> _PanelCache:
        local_indices = self._panel_local_indices(panel_index)
        subset = Subset(
            self.validation_dataset,
            [int(value) for value in local_indices.tolist()],
        )
        loader = DataLoader(
            subset,
            batch_size=int(self.config.eval_batch_size),
            shuffle=False,
            num_workers=0,
            pin_memory=False,
            drop_last=False,
        )
        n_donor = int(self.validation_dataset.n_donor)
        n_celltype = len(self.validation_dataset.spec.celltype_vocab)
        n_group = n_donor * n_celltype
        n_gene = int(self.decoder.n_genes)
        target_sum = torch.zeros(
            n_group, n_gene, dtype=torch.float32, device=self.device
        )
        count = torch.zeros(n_group, dtype=torch.long, device=self.device)
        group_ct = (
            torch.arange(n_group, device=self.device, dtype=torch.long)
            % n_celltype
        )
        cached: list[_CachedDecoderBatch] = []

        with _preserve_training_mode(self.system), torch.no_grad():
            for raw_batch in loader:
                batch = {
                    key: (
                        value.to(self.device, non_blocking=False)
                        if isinstance(value, torch.Tensor)
                        else value
                    )
                    for key, value in raw_batch.items()
                }
                out = self.system(batch, sample_latent=False)
                if out.decoder_out is None:
                    raise RuntimeError("generator search requires decoder output.")
                donor = batch["donor_id"].long()
                celltype = batch["celltype_id"].long()
                group = donor * n_celltype + celltype
                target = batch["y_ord"].float()
                target_sum.index_add_(0, group, target)
                count.index_add_(
                    0,
                    group,
                    torch.ones_like(group, dtype=torch.long),
                )
                cached.append(
                    _CachedDecoderBatch(
                        z=out.z_perp.detach(),
                        celltype=celltype.detach(),
                        tech=batch.get("tech_id", batch["batch_id"]).long().detach(),
                        sex=(
                            batch["sex_id"].long().detach()
                            if "sex_id" in batch
                            else None
                        ),
                        precision=(
                            out.decoder_out.precision_score.detach()
                            if out.decoder_out.precision_score is not None
                            else None
                        ),
                        target=batch["y_ord"].long().detach(),
                        donor=donor.detach(),
                    )
                )

        target_module, valid_group_ct = self._module_projection(
            gene_sum=target_sum,
            count=count,
            group_celltype=group_ct,
        )
        module_keep = (
            (self.decoder.module_masks.detach().cpu() > 0)
            .sum(dim=1)
            .numpy()
            > 1
        )
        if int(target_module.shape[1]) != int(module_keep.shape[0]):
            raise RuntimeError(
                "generator-budget module recovery requires a one-to-one "
                "candidate bank and tokenizer registry; got "
                f"{module_keep.shape[0]} generator rows and "
                f"{target_module.shape[1]} tokenizer modules."
            )
        return _PanelCache(
            batches=cached,
            target_module=target_module,
            group_celltype=valid_group_ct,
            module_keep=module_keep,
            n_cells=int(sum(batch.target.shape[0] for batch in cached)),
            n_groups=int(len(valid_group_ct)),
            n_celltypes=n_celltype,
        )

    def _decoder_score(
        self,
        batch: _CachedDecoderBatch,
        *,
        gate: torch.Tensor,
        isolated: bool,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        decoder = self.decoder
        base = decoder.celltype_baseline(batch.celltype)
        tech = decoder.tech_baseline(batch.tech)
        state = decoder.state_score_from_z_and_base(
            batch.z,
            base,
            celltype_id=batch.celltype,
            generator_gate=gate,
            include_direct=not isolated,
            include_pathology=not isolated,
            include_interactions=not isolated,
        )
        if decoder.score_mixer is not None:
            score, _ = decoder.score_mixer(
                z_perp=batch.z,
                base_score=base,
                tech_score=tech,
                state_score=state,
            )
        else:
            score = base + tech + state
        if decoder.sex_baseline is not None and batch.sex is not None:
            score = score + decoder.sex_baseline(decoder._sex_index(batch.sex))
        if not isolated and batch.precision is not None:
            score = score + batch.precision
        thresholds = decoder._compute_thresholds()
        _, probabilities = decoder._score_to_probs(score, thresholds)
        _, nll = decoder._nll_from_probs(probabilities, batch.target)
        bins = torch.arange(
            decoder.n_bins,
            dtype=probabilities.dtype,
            device=probabilities.device,
        )
        expected = (probabilities * bins.view(1, 1, -1)).sum(dim=-1)
        return nll, expected

    def _module_recovery(
        self,
        prediction_module: np.ndarray,
        cache: _PanelCache,
    ) -> float:
        values: list[float] = []
        keep = cache.module_keep
        for celltype in range(cache.n_celltypes):
            rows = cache.group_celltype == celltype
            if int(rows.sum()) < int(self.config.minimum_donors_per_celltype):
                continue
            prediction = prediction_module[rows][:, keep]
            target = cache.target_module[rows][:, keep]
            prediction = prediction - prediction.mean(axis=0, keepdims=True)
            target = target - target.mean(axis=0, keepdims=True)
            rho = _spearman(prediction.reshape(-1), target.reshape(-1))
            if math.isfinite(rho):
                values.append(rho)
        if not values:
            return float("nan")
        return float(np.median(np.asarray(values, dtype=np.float64)))

    def _exact_metrics(
        self,
        cache: _PanelCache,
        *,
        active_mask: torch.Tensor,
    ) -> dict[str, float]:
        n_donor = int(self.validation_dataset.n_donor)
        n_group = n_donor * cache.n_celltypes
        n_gene = int(self.decoder.n_genes)
        sums = {
            False: torch.zeros(
                n_group, n_gene, dtype=torch.float32, device=self.device
            ),
            True: torch.zeros(
                n_group, n_gene, dtype=torch.float32, device=self.device
            ),
        }
        counts = torch.zeros(n_group, dtype=torch.long, device=self.device)
        group_ct = (
            torch.arange(n_group, device=self.device, dtype=torch.long)
            % cache.n_celltypes
        )
        nll_sum = {False: 0.0, True: 0.0}
        n_cell = 0
        gate = active_mask.to(device=self.device, dtype=torch.float32)
        with _preserve_training_mode(self.decoder), torch.no_grad():
            old_stash = getattr(self.decoder, "_stash_diag", True)
            self.decoder._stash_diag = False
            try:
                for batch in cache.batches:
                    group = batch.donor * cache.n_celltypes + batch.celltype
                    counts.index_add_(
                        0,
                        group,
                        torch.ones_like(group, dtype=torch.long),
                    )
                    for isolated in (False, True):
                        nll, expected = self._decoder_score(
                            batch,
                            gate=gate,
                            isolated=isolated,
                        )
                        sums[isolated].index_add_(
                            0, group, expected.float()
                        )
                        nll_sum[isolated] += float(nll.double().sum().item())
                    n_cell += int(batch.target.shape[0])
            finally:
                self.decoder._stash_diag = old_stash

        result: dict[str, float] = {}
        for isolated, prefix in ((False, ""), (True, "isolated_")):
            module, valid_ct = self._module_projection(
                gene_sum=sums[isolated],
                count=counts,
                group_celltype=group_ct,
            )
            if not np.array_equal(valid_ct, cache.group_celltype):
                raise RuntimeError(
                    "prediction and target donor×celltype group order differs."
                )
            result[f"{prefix}module_recovery"] = self._module_recovery(
                module, cache
            )
            result[f"{'isolated_' if isolated else 'full_'}nll"] = (
                nll_sum[isolated] / max(n_cell, 1)
            )
        return result

    def _rank_removal_harm(
        self,
        cache: _PanelCache,
        *,
        active_mask: torch.Tensor,
    ) -> torch.Tensor:
        total_full = torch.zeros_like(active_mask, dtype=torch.float64)
        total_isolated = torch.zeros_like(active_mask, dtype=torch.float64)
        total_cells = 0
        with _preserve_training_mode(self.decoder), _frozen_parameters(
            self.decoder
        ):
            old_stash = getattr(self.decoder, "_stash_diag", True)
            self.decoder._stash_diag = False
            try:
                for batch in cache.batches:
                    gate = active_mask.detach().to(
                        device=self.device, dtype=torch.float32
                    )
                    gate.requires_grad_(True)
                    full_nll, _ = self._decoder_score(
                        batch, gate=gate, isolated=False
                    )
                    grad_full = torch.autograd.grad(
                        full_nll.mean(), gate, retain_graph=False
                    )[0]

                    gate_iso = active_mask.detach().to(
                        device=self.device, dtype=torch.float32
                    )
                    gate_iso.requires_grad_(True)
                    isolated_nll, _ = self._decoder_score(
                        batch, gate=gate_iso, isolated=True
                    )
                    grad_isolated = torch.autograd.grad(
                        isolated_nll.mean(), gate_iso, retain_graph=False
                    )[0]
                    cells = int(batch.target.shape[0])
                    total_full += (-grad_full.detach()).double() * cells
                    total_isolated += (
                        -grad_isolated.detach()
                    ).double() * cells
                    total_cells += cells
            finally:
                self.decoder._stash_diag = old_stash

        full_harm = total_full / max(total_cells, 1)
        isolated_harm = total_isolated / max(total_cells, 1)
        active = active_mask.bool()
        full_scale = full_harm[active].abs().median().clamp_min(1e-12)
        isolated_scale = (
            isolated_harm[active].abs().median().clamp_min(1e-12)
        )
        return torch.maximum(
            full_harm / full_scale,
            float(self.config.rank_isolated_weight)
            * isolated_harm
            / isolated_scale,
        )

    def _write_manifest(self) -> None:
        payload = self.controller.manifest()
        payload.update(
            {
                "search_config": asdict(self.config),
                "selector_views": [
                    "full_with_fixed_direct_and_prism",
                    "generator_isolated_without_direct_pathology_prism",
                ],
                "ranking": "paired_first_order_gate_harm_then_exact_block_ablation",
                "runtime_search": "hard_backward_elimination_from_all_active",
                "runtime_proposal_kinds": ["death"],
                "validation_panels": int(self.config.panel_count),
                "fixed_all_active_panel_baselines": {
                    str(panel): metrics
                    for panel, metrics in sorted(self.panel_baselines.items())
                },
                "test_split_opened": False,
            }
        )
        self.manifest_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)
            + "\n",
            encoding="utf-8",
        )

    def run_check(self, *, epoch: int) -> dict[str, Any]:
        if not self.due(epoch):
            return {
                "due": False,
                "epoch": int(epoch),
                "active_count": self.active_count,
            }
        panel_index = self.check_index % int(self.config.panel_count)
        cache = self._build_panel_cache(panel_index)
        active_mask = torch.tensor(
            self.controller.active_mask,
            dtype=torch.bool,
            device=self.device,
        )
        incumbent = self._exact_metrics(
            cache,
            active_mask=active_mask,
        )
        if panel_index not in self.panel_baselines:
            if not all(math.isfinite(float(value)) for value in incumbent.values()):
                raise RuntimeError(
                    "cannot establish a finite all-active generator baseline "
                    f"for validation panel {panel_index}: {incumbent}"
                )
            self.panel_baselines[panel_index] = {
                key: float(value) for key, value in incumbent.items()
            }

        if self.pending_block is None:
            harm = self._rank_removal_harm(
                cache,
                active_mask=active_mask,
            )
            eligible = [
                index
                for index in self.controller.active_ids
                if index not in set(self.controller.protected_ids)
            ]
            remaining_capacity = (
                self.active_count
                - int(self.config.minimum_active_generators)
            )
            take = min(
                int(self.config.block_size),
                len(eligible),
                remaining_capacity,
            )
            if take <= 0:
                return {
                    "due": False,
                    "epoch": int(epoch),
                    "active_count": self.active_count,
                    "reason": "minimum_active_reached",
                }
            block = tuple(
                sorted(
                    sorted(eligible, key=lambda index: (float(harm[index]), index))[
                        :take
                    ]
                )
            )
            ranking_values = {
                str(index): float(harm[index]) for index in block
            }
        else:
            block = self.pending_block
            ranking_values = {}

        candidate_mask = active_mask.clone()
        candidate_mask[list(block)] = False
        candidate = self._exact_metrics(
            cache,
            active_mask=candidate_mask,
        )
        fixed_baseline = self.panel_baselines[panel_index]
        fixed_baseline_reasons: list[str] = []
        for rule in self.controller.metric_rules:
            baseline_value = float(fixed_baseline[rule.name])
            candidate_value = float(candidate[rule.name])
            gain = rule.oriented_gain(baseline_value, candidate_value)
            if not (
                math.isfinite(baseline_value)
                and math.isfinite(candidate_value)
                and math.isfinite(gain)
            ):
                fixed_baseline_reasons.append(
                    f"fixed_baseline_nonfinite:{rule.name}"
                )
            elif gain + rule.atol < -rule.noninferiority_margin:
                fixed_baseline_reasons.append(
                    f"fixed_baseline_noninferiority:{rule.name}"
                )
        decision = self.controller.consider(
            step=self.check_index,
            kind="death",
            generator_ids=block,
            incumbent_metrics=incumbent,
            candidate_metrics=candidate,
            external_reason_codes=fixed_baseline_reasons,
        )
        if decision.applied:
            self.pending_block = None
        elif decision.eligible:
            self.pending_block = block
        else:
            self.pending_block = None
        self.decoder.set_generator_active_mask(
            torch.tensor(
                self.controller.active_mask,
                dtype=torch.bool,
                device=self.device,
            )
        )
        event = {
            "due": True,
            "epoch": int(epoch),
            "check_index": int(self.check_index),
            "panel_index": int(panel_index),
            "n_cells": int(cache.n_cells),
            "n_groups": int(cache.n_groups),
            "proposal": list(block),
            "proposal_rank_harm": ranking_values,
            "incumbent_metrics": incumbent,
            "candidate_metrics": candidate,
            "fixed_all_active_baseline": fixed_baseline,
            "fixed_baseline_reason_codes": fixed_baseline_reasons,
            "decision": decision.to_dict(),
            "active_count": self.active_count,
            "active_ids": list(self.controller.active_ids),
        }
        def _strict_json(value):
            if isinstance(value, float) and not math.isfinite(value):
                return None
            if isinstance(value, dict):
                return {str(key): _strict_json(item) for key, item in value.items()}
            if isinstance(value, (list, tuple)):
                return [_strict_json(item) for item in value]
            return value

        serializable_event = _strict_json(event)
        if not isinstance(serializable_event, dict):
            raise RuntimeError("generator-budget event serialization failed.")
        with self.event_path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(
                    serializable_event,
                    sort_keys=True,
                    allow_nan=False,
                )
                + "\n"
            )
        self.controller.write_state_json(self.state_path)
        self.check_index += 1
        self._write_manifest()
        return event
