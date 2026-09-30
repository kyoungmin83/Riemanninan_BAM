"""Learnable hard gates for constrained Lie-generator cardinality search.

This module contains only the architecture variables and the mathematical
non-inferiority objective.  It does not modify a decoder, split donors, open a
validation/test set, or run an optimiser.  The training integration must pass
``GateSample.straight_through`` to the decoder's existing
``generator_gate=...`` argument and must construct leakage-safe baseline and
candidate losses.

The forward gate is exactly binary.  Gradients use a stretched Binary-Concrete
surrogate, so a gate may turn off and later reactivate.  Cardinality is
minimised under explicit non-inferiority constraints rather than by adding an
L1 penalty to generator weights.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import torch
from torch import nn


@dataclass(frozen=True)
class LearnedGeneratorCountConfig:
    """Compatibility-safe configuration for learned generator gates.

    ``enabled=False`` is authoritative.  Merely constructing this dataclass
    must not alter a historical model or checkpoint.
    """

    enabled: bool = False
    # ``gate_only`` preserves the historical two-stage search runner.  In
    # ``joint`` mode the same exact-forward hard gates are part of the ordinary
    # model graph, so generator weights, PRISM and the architecture logits are
    # learned together in one full training run.  No candidate-K grid is used.
    mode: str = "gate_only"
    initial_keep_probability: float = 0.95
    temperature_start: float = 1.5
    temperature_end: float = 0.3
    stretch_low: float = -0.1
    stretch_high: float = 1.1
    hard_threshold: float = 0.5
    # None means no user-chosen cardinality floor.  The only lower bound then
    # emerges from protected singleton/unique-coverage routes and learned
    # non-inferiority constraints.
    minimum_active_generators: int | None = 1
    protect_singletons: bool = True
    coalition_samples: int = 4
    augmented_rho: float = 2.0
    dual_learning_rate: float = 0.05
    constraint_ema: float = 0.9
    # Runtime architecture-search contract.  Epochs 1..start_epoch-1 train an
    # all-on reference model; the separate gate-search runner labels its first
    # search pass ``start_epoch``.
    start_epoch: int = 13
    # Optional integrated schedule.  Zero preserves the historical exact-hard
    # behavior from start_epoch.  A non-zero shadow window learns counterfactual
    # importance while the live decoder remains all-on; a following soft window
    # lets weights adapt before reversible exact-hard masks begin.
    shadow_end_epoch: int = 0
    soft_end_epoch: int = 0
    base_end_epoch: int = 40
    max_extension_epochs: int = 8
    architecture_donor_count: int = 12
    architecture_split_seed: int = 420730
    architecture_split_trials: int = 20_000
    learning_rate: float = 0.03
    weight_decay: float = 0.0
    # A gate update is an actual optimiser step.  ``steps_per_epoch`` is the
    # number of complete architecture blocks accumulated into each update.
    # Keeping the two controls separate prevents a long block loop from being
    # mistaken for many opportunities to move the hard-gate logits.
    gate_updates_per_epoch: int = 8
    minimum_gate_updates: int = 96
    steps_per_epoch: int = 4
    # Cache several deterministic, non-overlapping cell draws for every
    # eligible cell type and rotate through them during gate updates.
    architecture_draws_per_celltype: int = 4
    cells_per_donor_celltype: int = 2
    minimum_donors_per_celltype: int = 8
    bootstrap_replicates: int = 200
    noninferiority_se_multiplier: float = 1.0
    constraint_scale_floor: float = 1e-4
    stability_window: int = 3
    stable_count_tolerance: int = 2
    stable_mask_jaccard: float = 0.98
    maximum_uncertain_fraction: float = 0.05
    minimum_search_epochs: int = 8
    # A saturated all-on mask is not evidence that 414 generators are
    # necessary.  At least this many generators must be removed before a mask
    # can be released as canonical; a future exact no-compression certificate
    # may provide a separate release path.
    minimum_pruned_generators: int = 1
    seed: int = 420732
    # Joint-training controls.  The paired all-on counterfactual is evaluated
    # on the same minibatch and with the same current weights.  The gate
    # objective minimises expected cardinality subject to the gated full and
    # generator-isolated NLLs remaining within these relative margins.
    joint_objective_weight: float = 0.05
    joint_full_nll_relative_margin: float = 0.01
    joint_isolated_nll_relative_margin: float = 0.015
    joint_constraint_scale_floor: float = 1e-3
    protect_unique_gene_coverage: bool = False
    safe_mask_confirmations: int = 3
    rollback_to_last_safe_mask: bool = True
    # New integrated learned-rank runs can require the pathology-rank mask to
    # be committed before the generator-count objective is allowed to start.
    # The cross-config epoch ordering is validated in runner_base.  Default
    # false preserves historical manifests.
    require_frozen_pathology_rank_before_search: bool = False

    def validate(self) -> None:
        if str(self.mode) not in {"gate_only", "joint"}:
            raise ValueError("mode must be either 'gate_only' or 'joint'.")
        if not (0.5 < float(self.initial_keep_probability) < 1.0):
            raise ValueError(
                "initial_keep_probability must lie strictly between 0.5 and 1."
            )
        if not (
            math.isfinite(float(self.temperature_start))
            and math.isfinite(float(self.temperature_end))
            and float(self.temperature_start) > 0.0
            and float(self.temperature_end) > 0.0
        ):
            raise ValueError("gate temperatures must be finite and positive.")
        if float(self.temperature_end) > float(self.temperature_start):
            raise ValueError("temperature_end cannot exceed temperature_start.")
        if not (
            math.isfinite(float(self.stretch_low))
            and math.isfinite(float(self.stretch_high))
            and float(self.stretch_low) < 0.0
            and float(self.stretch_high) > 1.0
        ):
            raise ValueError(
                "stretched Binary-Concrete bounds must satisfy low < 0 < 1 < high."
            )
        if not 0.0 < float(self.hard_threshold) < 1.0:
            raise ValueError("hard_threshold must lie strictly inside (0, 1).")
        if self.minimum_active_generators is not None and (
            isinstance(self.minimum_active_generators, bool)
            or int(self.minimum_active_generators) <= 0
        ):
            raise ValueError("minimum_active_generators must be a positive integer.")
        if (
            isinstance(self.coalition_samples, bool)
            or int(self.coalition_samples) <= 0
        ):
            raise ValueError("coalition_samples must be a positive integer.")
        if not math.isfinite(float(self.augmented_rho)) or float(
            self.augmented_rho
        ) < 0.0:
            raise ValueError("augmented_rho must be finite and non-negative.")
        if not math.isfinite(float(self.dual_learning_rate)) or float(
            self.dual_learning_rate
        ) <= 0.0:
            raise ValueError("dual_learning_rate must be finite and positive.")
        if not 0.0 <= float(self.constraint_ema) < 1.0:
            raise ValueError("constraint_ema must lie in [0, 1).")
        for name in (
            "start_epoch",
            "base_end_epoch",
            "architecture_donor_count",
            "architecture_split_trials",
            "gate_updates_per_epoch",
            "minimum_gate_updates",
            "steps_per_epoch",
            "architecture_draws_per_celltype",
            "cells_per_donor_celltype",
            "minimum_donors_per_celltype",
            "bootstrap_replicates",
            "stability_window",
            "minimum_search_epochs",
            "minimum_pruned_generators",
            "safe_mask_confirmations",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or int(value) <= 0:
                raise ValueError(f"{name} must be a positive integer.")
        if int(self.base_end_epoch) < int(self.start_epoch):
            raise ValueError("base_end_epoch cannot precede start_epoch.")
        if int(self.shadow_end_epoch) < 0 or int(self.soft_end_epoch) < 0:
            raise ValueError("shadow_end_epoch and soft_end_epoch must be non-negative")
        if int(self.shadow_end_epoch) > 0:
            if int(self.shadow_end_epoch) < int(self.start_epoch):
                raise ValueError("shadow_end_epoch cannot precede start_epoch")
            if int(self.soft_end_epoch) < int(self.shadow_end_epoch):
                raise ValueError("soft_end_epoch cannot precede shadow_end_epoch")
        if int(self.max_extension_epochs) < 0:
            raise ValueError("max_extension_epochs must be non-negative.")
        if int(self.stable_count_tolerance) < 0:
            raise ValueError("stable_count_tolerance must be non-negative.")
        for name in ("learning_rate", "constraint_scale_floor"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if not math.isfinite(float(self.weight_decay)) or float(
            self.weight_decay
        ) < 0.0:
            raise ValueError("weight_decay must be finite and non-negative.")
        if not math.isfinite(float(self.joint_objective_weight)) or float(
            self.joint_objective_weight
        ) <= 0.0:
            raise ValueError(
                "joint_objective_weight must be finite and positive."
            )
        for name in (
            "joint_full_nll_relative_margin",
            "joint_isolated_nll_relative_margin",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if not math.isfinite(float(self.joint_constraint_scale_floor)) or float(
            self.joint_constraint_scale_floor
        ) <= 0.0:
            raise ValueError(
                "joint_constraint_scale_floor must be finite and positive."
            )
        if not math.isfinite(float(self.noninferiority_se_multiplier)) or float(
            self.noninferiority_se_multiplier
        ) < 0.0:
            raise ValueError(
                "noninferiority_se_multiplier must be finite/non-negative."
            )
        if not 0.0 <= float(self.stable_mask_jaccard) <= 1.0:
            raise ValueError("stable_mask_jaccard must lie in [0, 1].")
        if not 0.0 <= float(self.maximum_uncertain_fraction) <= 1.0:
            raise ValueError(
                "maximum_uncertain_fraction must lie in [0, 1]."
            )


@dataclass(frozen=True)
class GateSample:
    """One stochastic coalition and its differentiable surrogate."""

    hard: torch.Tensor
    straight_through: torch.Tensor
    soft: torch.Tensor
    keep_probability: torch.Tensor

    @property
    def active_count(self) -> int:
        return int(self.hard.detach().sum().item())


class HardBinaryConcreteGeneratorGate(nn.Module):
    """One learnable exact-forward switch per candidate generator.

    Protected generators are always one in both forward and backward paths.
    When an optional stochastic floor is configured and a draw falls below it,
    the generators with the largest keep probabilities are forced on.
    The soft backward path is retained for every unprotected generator.
    """

    def __init__(
        self,
        num_generators: int,
        *,
        config: LearnedGeneratorCountConfig,
        protected_mask: torch.Tensor | None = None,
    ) -> None:
        super().__init__()
        config.validate()
        if isinstance(num_generators, bool) or int(num_generators) <= 0:
            raise ValueError("num_generators must be a positive integer.")
        self.num_generators = int(num_generators)
        if (
            config.minimum_active_generators is not None
            and int(config.minimum_active_generators) > self.num_generators
        ):
            raise ValueError(
                "minimum_active_generators cannot exceed num_generators."
            )
        self.config = config

        if protected_mask is None:
            protected = torch.zeros(self.num_generators, dtype=torch.bool)
        else:
            protected = torch.as_tensor(protected_mask, dtype=torch.bool)
            if protected.shape != (self.num_generators,):
                raise ValueError(
                    "protected_mask must have shape "
                    f"({self.num_generators},), got {tuple(protected.shape)}."
                )
        if int(protected.sum()) > self.num_generators:
            raise ValueError("protected_mask contains too many protected entries.")
        self.register_buffer("protected_mask", protected.clone(), persistent=True)

        probability = float(config.initial_keep_probability)
        threshold_relaxed = self._threshold_in_relaxed_space()
        threshold_logit = math.log(
            threshold_relaxed / (1.0 - threshold_relaxed)
        )
        # P(hard=1) = sigmoid(log_alpha - tau*logit(threshold_relaxed)).
        log_alpha = (
            math.log(probability / (1.0 - probability))
            + float(config.temperature_start) * threshold_logit
        )
        self.log_alpha = nn.Parameter(
            torch.full((self.num_generators,), log_alpha, dtype=torch.float32)
        )
        self.register_buffer(
            "last_safe_log_alpha",
            self.log_alpha.detach().clone(),
            persistent=True,
        )
        self.register_buffer(
            "safe_audit_streak", torch.tensor(0, dtype=torch.int64), persistent=True
        )
        self.register_buffer(
            "has_confirmed_safe_mask", torch.tensor(False), persistent=True
        )

    @torch.no_grad()
    def apply_safety_audit(self, *, feasible: bool) -> dict[str, int | bool]:
        """Confirm a safe state or roll back logits after a failed audit."""

        rolled_back = False
        if bool(feasible):
            self.safe_audit_streak.add_(1)
            if int(self.safe_audit_streak.item()) >= int(
                self.config.safe_mask_confirmations
            ):
                self.last_safe_log_alpha.copy_(self.log_alpha.detach())
                self.has_confirmed_safe_mask.fill_(True)
        else:
            self.safe_audit_streak.zero_()
            if bool(self.config.rollback_to_last_safe_mask) and bool(
                self.has_confirmed_safe_mask.item()
            ):
                self.log_alpha.copy_(self.last_safe_log_alpha)
                rolled_back = True
        return {
            "feasible": bool(feasible),
            "streak": int(self.safe_audit_streak.item()),
            "has_confirmed_safe_mask": bool(
                self.has_confirmed_safe_mask.item()
            ),
            "rolled_back": bool(rolled_back),
        }

    def _threshold_in_relaxed_space(self) -> float:
        low = float(self.config.stretch_low)
        high = float(self.config.stretch_high)
        value = (float(self.config.hard_threshold) - low) / (high - low)
        if not 0.0 < value < 1.0:
            raise ValueError(
                "hard_threshold maps outside the unstretched Binary-Concrete range."
            )
        return value

    def temperature(self, progress: float) -> float:
        """Log-linearly anneal temperature for progress in ``[0, 1]``."""

        if not math.isfinite(float(progress)):
            raise ValueError("progress must be finite.")
        progress = min(max(float(progress), 0.0), 1.0)
        start = float(self.config.temperature_start)
        end = float(self.config.temperature_end)
        return float(math.exp((1.0 - progress) * math.log(start) + progress * math.log(end)))

    def keep_probability(self, *, temperature: float) -> torch.Tensor:
        """Analytic probability that the exact hard forward switch is one."""

        if not math.isfinite(float(temperature)) or float(temperature) <= 0.0:
            raise ValueError("temperature must be finite and positive.")
        with torch.autocast(
            device_type=self.log_alpha.device.type, enabled=False
        ):
            threshold = self._threshold_in_relaxed_space()
            threshold_logit = math.log(threshold / (1.0 - threshold))
            probability = torch.sigmoid(
                self.log_alpha.float() - float(temperature) * threshold_logit
            )
            return torch.where(
                self.protected_mask,
                torch.ones_like(probability),
                probability,
            )

    def _enforce_hard_contract(
        self,
        hard: torch.Tensor,
        probability: torch.Tensor,
    ) -> torch.Tensor:
        hard = torch.where(
            self.protected_mask,
            torch.ones_like(hard),
            hard,
        )
        if self.config.minimum_active_generators is None:
            return hard
        minimum = int(self.config.minimum_active_generators)
        if int(hard.detach().sum().item()) >= minimum:
            return hard
        top = torch.topk(probability, k=minimum, largest=True, sorted=False).indices
        forced = hard.clone()
        forced[top] = 1.0
        return forced

    def sample(
        self,
        *,
        temperature: float,
        uniform: torch.Tensor | None = None,
    ) -> GateSample:
        """Sample an exact 0/1 coalition with a straight-through gradient."""

        with torch.autocast(
            device_type=self.log_alpha.device.type, enabled=False
        ):
            probability = self.keep_probability(temperature=temperature)
            if uniform is None:
                uniform = torch.rand_like(self.log_alpha)
            else:
                uniform = torch.as_tensor(
                    uniform,
                    dtype=torch.float32,
                    device=self.log_alpha.device,
                )
                if uniform.shape != self.log_alpha.shape:
                    raise ValueError(
                        f"uniform must have shape {tuple(self.log_alpha.shape)}."
                    )
            eps = torch.finfo(torch.float32).eps
            uniform = uniform.float().clamp(min=eps, max=1.0 - eps)
            logistic = torch.log(uniform) - torch.log1p(-uniform)
            relaxed = torch.sigmoid(
                (logistic + self.log_alpha.float()) / float(temperature)
            )
            stretched = (
                relaxed
                * (float(self.config.stretch_high) - float(self.config.stretch_low))
                + float(self.config.stretch_low)
            )
            soft = stretched.clamp(0.0, 1.0)
            hard = (soft >= float(self.config.hard_threshold)).to(soft.dtype)
            hard = self._enforce_hard_contract(hard, probability)

            protected_float = self.protected_mask.to(dtype=soft.dtype)
            soft = torch.where(self.protected_mask, protected_float, soft)
            # Exact hard values in forward; d/dalpha follows ``soft``.
            straight_through = hard.detach() - soft.detach() + soft
            return GateSample(
                hard=hard.detach(),
                straight_through=straight_through,
                soft=soft,
                keep_probability=probability,
            )

    def deterministic_mask(self, *, temperature: float) -> torch.Tensor:
        """Return the reproducible hard architecture implied by current logits."""

        probability = self.keep_probability(temperature=temperature)
        hard = (
            probability >= float(self.config.hard_threshold)
        ).to(probability.dtype)
        return self._enforce_hard_contract(hard, probability).detach()

    def deterministic_straight_through(
        self, *, temperature: float
    ) -> GateSample:
        """Return a repeatable exact hard mask with probability gradients.

        This is used only by deterministic grouped module-rescue updates: the
        same architecture is evaluated on every rank, while recovery gradients
        can still tell the gate logits which routes must remain available.
        """

        probability = self.keep_probability(temperature=temperature)
        hard = (
            probability >= float(self.config.hard_threshold)
        ).to(probability.dtype)
        hard = self._enforce_hard_contract(hard, probability).detach()
        straight_through = hard - probability.detach() + probability
        return GateSample(
            hard=hard,
            straight_through=straight_through,
            soft=probability,
            keep_probability=probability,
        )

    def raw_expected_active_count(self, *, temperature: float) -> torch.Tensor:
        """Return the unconstrained sum of analytic keep probabilities."""

        return self.keep_probability(temperature=temperature).sum()

    def expected_active_count(self, *, temperature: float) -> torch.Tensor:
        """Return expected cardinality consistent with the hard forward floor.

        The analytic expectation is clamped only to constraints that also bind
        the hard forward: protected routes and, when non-null, the configured
        technical floor.
        """

        raw = self.raw_expected_active_count(temperature=temperature)
        configured_floor = (
            0
            if self.config.minimum_active_generators is None
            else int(self.config.minimum_active_generators)
        )
        floor = max(configured_floor, int(self.protected_mask.sum().item()))
        return raw.clamp_min(float(floor))


@dataclass(frozen=True)
class ConstraintBatch:
    """Loss-space non-inferiority inputs.

    Every metric must be oriented so that *smaller is better*.  ``margin`` and
    ``scale`` are fixed train-architecture-donor quantities.  A positive
    normalized violation means the gated model is worse than allowed.
    """

    baseline: Mapping[str, torch.Tensor]
    candidate: Mapping[str, torch.Tensor]
    margin: Mapping[str, float]
    scale: Mapping[str, float]


@dataclass(frozen=True)
class ConstrainedCountOutput:
    loss: torch.Tensor
    cardinality_fraction: torch.Tensor
    violations: Mapping[str, torch.Tensor]
    positive_violations: Mapping[str, torch.Tensor]
    dual_values: Mapping[str, torch.Tensor]


class ConstrainedGeneratorCountObjective(nn.Module):
    """Augmented-Lagrangian objective for minimum feasible cardinality."""

    def __init__(
        self,
        metric_names: Sequence[str],
        *,
        config: LearnedGeneratorCountConfig,
    ) -> None:
        super().__init__()
        config.validate()
        names = tuple(str(name) for name in metric_names)
        if not names or any(not name for name in names):
            raise ValueError("metric_names must contain non-empty names.")
        if len(set(names)) != len(names):
            raise ValueError("metric_names must be unique.")
        self.metric_names = names
        self.config = config
        self.register_buffer(
            "dual",
            torch.zeros(len(names), dtype=torch.float32),
            persistent=True,
        )
        self.register_buffer(
            "violation_ema",
            torch.zeros(len(names), dtype=torch.float32),
            persistent=True,
        )
        self.register_buffer(
            "ema_initialized",
            torch.tensor(False, dtype=torch.bool),
            persistent=True,
        )

    def _validate_batch(self, batch: ConstraintBatch) -> None:
        expected = set(self.metric_names)
        for label, values in (
            ("baseline", batch.baseline),
            ("candidate", batch.candidate),
            ("margin", batch.margin),
            ("scale", batch.scale),
        ):
            if set(values) != expected:
                raise ValueError(
                    f"{label} keys must be exactly {sorted(expected)}, "
                    f"got {sorted(values)}."
                )
        for name in self.metric_names:
            baseline = batch.baseline[name]
            candidate = batch.candidate[name]
            if not isinstance(baseline, torch.Tensor) or baseline.numel() != 1:
                raise ValueError(f"baseline[{name!r}] must be a scalar tensor.")
            if not isinstance(candidate, torch.Tensor) or candidate.numel() != 1:
                raise ValueError(f"candidate[{name!r}] must be a scalar tensor.")
            if not bool(torch.isfinite(baseline).all()) or not bool(
                torch.isfinite(candidate).all()
            ):
                raise ValueError(f"constraint metric {name!r} is non-finite.")
            margin = float(batch.margin[name])
            scale = float(batch.scale[name])
            if not math.isfinite(margin) or margin < 0.0:
                raise ValueError(f"margin[{name!r}] must be finite and non-negative.")
            if not math.isfinite(scale) or scale <= 0.0:
                raise ValueError(f"scale[{name!r}] must be finite and positive.")

    def forward(
        self,
        *,
        expected_active_count: torch.Tensor,
        num_generators: int,
        constraints: ConstraintBatch,
    ) -> ConstrainedCountOutput:
        self._validate_batch(constraints)
        if not isinstance(expected_active_count, torch.Tensor) or (
            expected_active_count.numel() != 1
        ):
            raise ValueError("expected_active_count must be a scalar tensor.")
        if isinstance(num_generators, bool) or int(num_generators) <= 0:
            raise ValueError("num_generators must be a positive integer.")
        cardinality = expected_active_count.reshape(()) / float(num_generators)
        violations: dict[str, torch.Tensor] = {}
        positive: dict[str, torch.Tensor] = {}
        loss = cardinality
        rho = float(self.config.augmented_rho)
        for index, name in enumerate(self.metric_names):
            violation = (
                constraints.candidate[name]
                - constraints.baseline[name].detach()
                - float(constraints.margin[name])
            ) / float(constraints.scale[name])
            violation = violation.reshape(())
            positive_violation = torch.relu(violation)
            violations[name] = violation
            positive[name] = positive_violation
            loss = (
                loss
                + self.dual[index].detach() * positive_violation
                + 0.5 * rho * positive_violation.square()
            )
        return ConstrainedCountOutput(
            loss=loss,
            cardinality_fraction=cardinality,
            violations=violations,
            positive_violations=positive,
            dual_values={
                name: self.dual[index].detach()
                for index, name in enumerate(self.metric_names)
            },
        )

    @torch.no_grad()
    def update_dual(self, violations: Mapping[str, torch.Tensor]) -> None:
        if set(violations) != set(self.metric_names):
            raise ValueError("violations keys do not match metric_names.")
        current = torch.stack(
            [
                torch.as_tensor(
                    violations[name],
                    dtype=self.dual.dtype,
                    device=self.dual.device,
                ).reshape(())
                for name in self.metric_names
            ]
        )
        if not bool(torch.isfinite(current).all()):
            raise ValueError("dual update received non-finite violations.")
        decay = float(self.config.constraint_ema)
        if bool(self.ema_initialized):
            self.violation_ema.mul_(decay).add_(current, alpha=1.0 - decay)
        else:
            self.violation_ema.copy_(current)
            self.ema_initialized.fill_(True)
        self.dual.add_(
            float(self.config.dual_learning_rate) * self.violation_ema
        ).clamp_(min=0.0)


__all__ = [
    "ConstrainedCountOutput",
    "ConstrainedGeneratorCountObjective",
    "ConstraintBatch",
    "GateSample",
    "HardBinaryConcreteGeneratorGate",
    "LearnedGeneratorCountConfig",
]
