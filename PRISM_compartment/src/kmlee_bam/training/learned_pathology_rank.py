"""Learn the effective rank of the pathology-conditioned decoder correction.

The tensor allocation uses the algebraic rank capacity implied by the data
shape.  It is not a selected rank.  A shared, structured gate learns how many
of those basis components are actually used, without a hand-set minimum or a
target rank.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn as nn


def algebraic_pathology_rank_capacity(
    n_celltypes: int,
    n_pathology_axes: int,
    n_genes: int,
) -> int:
    """Return the largest possible rank of the stacked correction matrix."""

    values = (int(n_celltypes), int(n_pathology_axes), int(n_genes))
    if any(value <= 0 for value in values):
        raise ValueError("celltypes, pathology axes, and genes must be positive")
    return min(values[0] * values[1], values[2])


@dataclass(frozen=True)
class LearnedPathologyRankConfig:
    enabled: bool = False
    initial_keep_probability: float = 0.995
    # Optional one-architecture canonical warm-up.  The full algebraic basis is
    # allocated before epoch one, but only the first ``warmup_fixed_rank``
    # components participate through ``warmup_end_epoch``.  This makes an
    # integrated learned-rank run functionally equivalent to the historical
    # fixed-rank decoder during its biological warm-up without checkpoint
    # surgery or a second optimizer lineage.
    warmup_fixed_rank: int = 0
    warmup_end_epoch: int = 0
    # ``None`` preserves the historical all-components initialization.  A
    # separate probability lets the additional components enter the rank-
    # learning stage at a neutral value while the canonical warm-up components
    # retain their established initialization.  There is no cross-fade: at
    # ``soft_start_epoch`` the full learned mask becomes active immediately.
    initial_extra_keep_probability: float | None = None
    # These flags opt a new primary run into stricter cross-section validation
    # in runner_base while leaving historical manifests loadable.
    require_canonical_phase1_rank8: bool = False
    require_active_pathology_route_during_search: bool = False
    temperature_start: float = 2.0
    temperature_end: float = 0.35
    soft_start_epoch: int = 1
    hard_start_epoch: int = 13
    freeze_epoch: int = 17
    sparsity_start_epoch: int = 3
    sparsity_ramp_epochs: int = 6
    cardinality_loss_fraction: float = 0.0025
    learning_rate: float = 3.0e-5
    weight_decay: float = 0.0
    hard_threshold: float = 0.5
    # Optional fail-closed structural commit gate.  It prevents a threshold-
    # sensitive rank mask from being frozen merely because the calendar
    # reached ``freeze_epoch``.  Historical configurations keep it disabled.
    require_stable_freeze: bool = False
    freeze_threshold_low: float = 0.45
    freeze_threshold_high: float = 0.55
    maximum_freeze_count_spread: int = 4
    maximum_freeze_uncertain_fraction: float = 0.15

    def validate(self) -> None:
        if not bool(self.enabled):
            return
        if not 0.5 < float(self.initial_keep_probability) < 1.0:
            raise ValueError("initial_keep_probability must lie in (0.5, 1)")
        if isinstance(self.warmup_fixed_rank, bool) or int(
            self.warmup_fixed_rank
        ) < 0:
            raise ValueError("warmup_fixed_rank must be a non-negative integer")
        if isinstance(self.warmup_end_epoch, bool) or int(
            self.warmup_end_epoch
        ) < 0:
            raise ValueError("warmup_end_epoch must be a non-negative integer")
        warmup_rank = int(self.warmup_fixed_rank)
        warmup_end = int(self.warmup_end_epoch)
        if (warmup_rank == 0) != (warmup_end == 0):
            raise ValueError(
                "warmup_fixed_rank and warmup_end_epoch must both be zero or "
                "both be positive"
            )
        if warmup_rank > 0 and int(self.soft_start_epoch) <= warmup_end:
            raise ValueError(
                "soft_start_epoch must follow the fixed-rank warm-up"
            )
        if bool(self.require_canonical_phase1_rank8) and warmup_rank != 8:
            raise ValueError(
                "require_canonical_phase1_rank8 requires warmup_fixed_rank=8"
            )
        if self.initial_extra_keep_probability is not None and not (
            0.0 < float(self.initial_extra_keep_probability) < 1.0
        ):
            raise ValueError(
                "initial_extra_keep_probability must lie in (0, 1)"
            )
        if float(self.temperature_start) <= 0.0 or float(self.temperature_end) <= 0.0:
            raise ValueError("pathology-rank temperatures must be positive")
        if int(self.soft_start_epoch) < 1:
            raise ValueError("soft_start_epoch must be positive")
        if int(self.hard_start_epoch) <= int(self.soft_start_epoch):
            raise ValueError("hard_start_epoch must follow soft_start_epoch")
        if int(self.freeze_epoch) <= int(self.hard_start_epoch):
            raise ValueError("freeze_epoch must follow hard_start_epoch")
        if int(self.sparsity_start_epoch) < int(self.soft_start_epoch):
            raise ValueError("sparsity cannot begin before soft rank learning")
        if int(self.sparsity_ramp_epochs) < 1:
            raise ValueError("sparsity_ramp_epochs must be positive")
        if not 0.0 < float(self.cardinality_loss_fraction) < 0.1:
            raise ValueError("cardinality_loss_fraction must lie in (0, 0.1)")
        if float(self.learning_rate) <= 0.0:
            raise ValueError("pathology-rank learning_rate must be positive")
        if float(self.weight_decay) < 0.0:
            raise ValueError("pathology-rank weight_decay cannot be negative")
        if not 0.0 < float(self.hard_threshold) < 1.0:
            raise ValueError("hard_threshold must lie in (0, 1)")
        low = float(self.freeze_threshold_low)
        high = float(self.freeze_threshold_high)
        if not 0.0 < low < float(self.hard_threshold) < high < 1.0:
            raise ValueError(
                "freeze thresholds must satisfy 0 < low < hard < high < 1"
            )
        if isinstance(self.maximum_freeze_count_spread, bool) or int(
            self.maximum_freeze_count_spread
        ) < 0:
            raise ValueError(
                "maximum_freeze_count_spread must be a non-negative integer"
            )
        if not 0.0 <= float(self.maximum_freeze_uncertain_fraction) <= 1.0:
            raise ValueError(
                "maximum_freeze_uncertain_fraction must lie in [0, 1]"
            )


class LearnedPathologyRankGate(nn.Module):
    """Target-free structured rank selection with a checkpointed hard result.

    Every algebraically possible component is allocated at initialization.
    An optional prefix mask keeps the canonical fixed-rank decoder exact during
    biological warm-up.  At ``soft_start_epoch`` the full differentiable mask
    is activated directly; epochs from ``hard_start_epoch`` use an exact hard
    forward with a straight-through gradient.  At
    ``freeze_epoch`` the learned mask is checkpointed and cannot collapse
    merely because the pathology route is subsequently dormant.
    """

    def __init__(
        self,
        capacity: int,
        *,
        config: LearnedPathologyRankConfig,
    ) -> None:
        super().__init__()
        config.validate()
        if int(capacity) <= 0:
            raise ValueError("pathology rank capacity must be positive")
        self.capacity = int(capacity)
        self.config = config
        warmup_rank = int(config.warmup_fixed_rank)
        if warmup_rank > self.capacity:
            raise ValueError(
                "warmup_fixed_rank cannot exceed pathology rank capacity"
            )
        probability = float(config.initial_keep_probability)
        extra_probability = (
            probability
            if config.initial_extra_keep_probability is None
            else float(config.initial_extra_keep_probability)
        )
        initial_probability = torch.full(
            (self.capacity,), extra_probability, dtype=torch.float32
        )
        if warmup_rank > 0:
            initial_probability[:warmup_rank] = probability
        initial_logit = torch.log(initial_probability) - torch.log1p(
            -initial_probability
        )
        self.log_alpha = nn.Parameter(initial_logit)
        self.register_buffer(
            "epoch_state",
            torch.tensor(1, dtype=torch.int64),
            persistent=True,
        )
        self.register_buffer(
            "finalized_state",
            torch.tensor(False, dtype=torch.bool),
            persistent=True,
        )
        self.register_buffer(
            "frozen_mask",
            torch.ones(self.capacity, dtype=torch.float32),
            persistent=True,
        )
        warmup_mask = torch.ones(self.capacity, dtype=torch.float32)
        if warmup_rank > 0:
            warmup_mask.zero_()
            warmup_mask[:warmup_rank] = 1.0
        self.register_buffer(
            "warmup_mask",
            warmup_mask,
            # Deterministically reconstructed from the immutable config.  Do
            # not add a new checkpoint key: historical learned-rank checkpoints
            # must remain strict-load compatible.
            persistent=False,
        )

    @property
    def epoch(self) -> int:
        return int(self.epoch_state.item())

    def temperature(self, epoch: int | None = None) -> float:
        epoch_value = self.epoch if epoch is None else int(epoch)
        start = int(self.config.soft_start_epoch)
        end = max(int(self.config.hard_start_epoch) - 1, start)
        progress = min(max((epoch_value - start) / max(end - start, 1), 0.0), 1.0)
        log_temperature = (
            math.log(float(self.config.temperature_start))
            + progress
            * (
                math.log(float(self.config.temperature_end))
                - math.log(float(self.config.temperature_start))
            )
        )
        return float(math.exp(log_temperature))

    def keep_probability(self) -> torch.Tensor:
        return torch.sigmoid(self.log_alpha)

    def deterministic_hard_mask(self) -> torch.Tensor:
        return (
            self.keep_probability() >= float(self.config.hard_threshold)
        ).to(dtype=self.log_alpha.dtype)

    @torch.no_grad()
    def freeze_readiness(self) -> dict[str, float | int | bool]:
        """Return threshold-sensitivity evidence for the rank commit."""

        probability = self.keep_probability().detach().float()
        low = float(self.config.freeze_threshold_low)
        high = float(self.config.freeze_threshold_high)
        count_low = int((probability >= low).sum().item())
        count_mid = int(
            (probability >= float(self.config.hard_threshold)).sum().item()
        )
        count_high = int((probability >= high).sum().item())
        spread = count_low - count_high
        uncertain_fraction = float(
            ((probability >= low) & (probability < high))
            .float()
            .mean()
            .cpu()
        )
        ready = bool(
            spread <= int(self.config.maximum_freeze_count_spread)
            and uncertain_fraction
            <= float(self.config.maximum_freeze_uncertain_fraction)
        )
        return {
            "ready": ready,
            "count_at_low_threshold": count_low,
            "count_at_hard_threshold": count_mid,
            "count_at_high_threshold": count_high,
            "count_spread": spread,
            "near_threshold_uncertain_fraction": uncertain_fraction,
        }

    @torch.no_grad()
    def set_epoch(self, epoch: int) -> None:
        epoch_value = int(epoch)
        if epoch_value < 1:
            raise ValueError("pathology-rank epoch must be positive")
        self.epoch_state.fill_(epoch_value)
        if epoch_value >= int(self.config.freeze_epoch) and not bool(
            self.finalized_state.item()
        ):
            readiness = self.freeze_readiness()
            if bool(self.config.require_stable_freeze) and not bool(
                readiness["ready"]
            ):
                raise RuntimeError(
                    "pathology-rank freeze rejected: "
                    f"count_spread={readiness['count_spread']} "
                    f"uncertain_fraction="
                    f"{readiness['near_threshold_uncertain_fraction']:.6f}"
                )
            self.frozen_mask.copy_(self.deterministic_hard_mask())
            self.finalized_state.fill_(True)

    def mode(self) -> str:
        if bool(self.finalized_state.item()):
            return "frozen_hard"
        if (
            int(self.config.warmup_fixed_rank) > 0
            and self.epoch <= int(self.config.warmup_end_epoch)
        ):
            return "fixed_rank_warmup"
        if self.epoch < int(self.config.soft_start_epoch):
            return "all_open"
        if self.epoch < int(self.config.hard_start_epoch):
            return "soft_rank_learning"
        return "hard_straight_through"

    def forward(self) -> torch.Tensor:
        mode = self.mode()
        if mode == "fixed_rank_warmup":
            # Exact fixed-rank forward with an explicit zero derivative keeps
            # the gate parameter in DDP's graph while preventing warm-up drift.
            return self.warmup_mask.to(dtype=self.log_alpha.dtype) + 0.0 * self.log_alpha
        if mode == "all_open":
            return torch.ones_like(self.log_alpha) + 0.0 * self.log_alpha
        if mode == "frozen_hard":
            # Keep the parameter in the DDP graph with an exactly zero
            # derivative; the frozen architecture itself is immutable.
            return self.frozen_mask.to(dtype=self.log_alpha.dtype) + 0.0 * self.log_alpha
        if mode == "soft_rank_learning":
            return self.keep_probability()
        temperature = self.temperature()
        soft = torch.sigmoid(self.log_alpha / temperature)
        hard = self.deterministic_hard_mask()
        return hard + soft - soft.detach()

    def sparsity_multiplier(self) -> float:
        epoch = self.epoch
        start = int(self.config.sparsity_start_epoch)
        if epoch < start or bool(self.finalized_state.item()):
            return 0.0
        ramp = int(self.config.sparsity_ramp_epochs)
        return float(min(1.0, (epoch - start + 1) / float(ramp)))

    def cardinality_penalty(self) -> torch.Tensor:
        """Normalized expected cardinality; no target or lower-bound term."""

        return self.keep_probability().mean()

    @torch.no_grad()
    def diagnostics(self) -> dict[str, float | int | str | bool]:
        probability = self.keep_probability().detach().float()
        mode = self.mode()
        search_hard = self.deterministic_hard_mask().detach().float()
        if mode == "fixed_rank_warmup":
            hard = self.warmup_mask.detach().float()
            expected_rank = float(hard.sum().cpu())
        elif mode == "all_open":
            hard = torch.ones_like(search_hard)
            expected_rank = float(self.capacity)
        elif mode == "frozen_hard":
            hard = self.frozen_mask.detach().float()
            expected_rank = float(hard.sum().cpu())
        else:
            hard = search_hard
            expected_rank = float(probability.sum().cpu())
        readiness = self.freeze_readiness()
        return {
            "capacity": int(self.capacity),
            "expected_rank": expected_rank,
            "hard_rank": int(hard.sum().item()),
            "search_expected_rank": float(probability.sum().cpu()),
            "search_hard_rank": int(search_hard.sum().item()),
            "uncertain_fraction": float(
                ((probability > 0.1) & (probability < 0.9)).float().mean().cpu()
            ),
            "temperature": float(self.temperature()),
            "mode": mode,
            "finalized": bool(self.finalized_state.item()),
            "warmup_fixed_rank": int(self.config.warmup_fixed_rank),
            "warmup_end_epoch": int(self.config.warmup_end_epoch),
            "freeze_ready": bool(readiness["ready"]),
            "freeze_count_at_low_threshold": int(
                readiness["count_at_low_threshold"]
            ),
            "freeze_count_at_hard_threshold": int(
                readiness["count_at_hard_threshold"]
            ),
            "freeze_count_at_high_threshold": int(
                readiness["count_at_high_threshold"]
            ),
            "freeze_count_spread": int(readiness["count_spread"]),
            "freeze_near_threshold_uncertain_fraction": float(
                readiness["near_threshold_uncertain_fraction"]
            ),
            "minimum_rank": 0,
            "target_rank": "none",
        }
