"""Leakage-safe pairwise pathology interaction features for PRISM.

The five named pathology main effects already use flexible monotone hinge
features.  A raw product such as ``Thal * Lewy`` is therefore not, by itself,
an identifiable interaction: part of that product can be explained by the two
main effects.  This module fits a donor-level, train-only residualizer and
returns only the part of each pair product that cannot be represented by the
two axes' zero-at-normal hinge bases.  No intercept is used, so the interaction
is exactly zero when both pathologies are zero and cannot contaminate the
common-normal component.

The fitting unit is a donor, never a cell.  Missing pathology is not treated as
normal: a pair feature is valid only when both named axes were observed.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import torch
from torch import nn
import torch.nn.functional as F


PATHOLOGY_NAMES: tuple[str, ...] = ("thal", "braak", "cerad", "late", "lewy")
PATHOLOGY_PAIRS: tuple[tuple[int, int], ...] = tuple(
    combinations(range(len(PATHOLOGY_NAMES)), 2)
)
PATHOLOGY_PAIR_NAMES: tuple[str, ...] = tuple(
    f"{PATHOLOGY_NAMES[left]}__{PATHOLOGY_NAMES[right]}"
    for left, right in PATHOLOGY_PAIRS
)


def _hinge_features(pathology: torch.Tensor) -> torch.Tensor:
    if pathology.ndim != 2 or int(pathology.shape[1]) != len(PATHOLOGY_NAMES):
        raise ValueError(
            f"pathology must have shape [N,{len(PATHOLOGY_NAMES)}], "
            f"got {tuple(pathology.shape)}."
        )
    return torch.stack(
        (
            pathology,
            F.relu(pathology - (1.0 / 3.0)),
            F.relu(pathology - (2.0 / 3.0)),
        ),
        dim=-1,
    )


@dataclass(frozen=True)
class PairwisePathologyResidualizerArtifact:
    """Train-donor statistics defining the ten residualized pair features."""

    projection: torch.Tensor       # [Q, 6]: 3 left + 3 right
    normalization_scale: torch.Tensor  # [Q], raw pair-product SD
    residual_rms: torch.Tensor         # [Q], after main-effect residualization
    residual_fraction: torch.Tensor    # [Q], residual_rms / normalization_scale
    eligible: torch.Tensor         # [Q] bool
    complete_donors: torch.Tensor  # [Q] long
    pairs: torch.Tensor            # [Q, 2] long

    def validate(self) -> None:
        q = len(PATHOLOGY_PAIRS)
        if self.projection.shape != (q, 6):
            raise ValueError(
                f"projection must have shape ({q},6), "
                f"got {tuple(self.projection.shape)}."
            )
        for name, tensor in (
            ("normalization_scale", self.normalization_scale),
            ("residual_rms", self.residual_rms),
            ("residual_fraction", self.residual_fraction),
            ("eligible", self.eligible),
            ("complete_donors", self.complete_donors),
        ):
            if tensor.shape != (q,):
                raise ValueError(f"{name} must have shape ({q},).")
        if self.pairs.shape != (q, 2):
            raise ValueError(f"pairs must have shape ({q},2).")
        if not bool(torch.isfinite(self.projection).all()):
            raise ValueError("projection contains non-finite values.")
        for name, tensor in (
            ("normalization_scale", self.normalization_scale),
            ("residual_rms", self.residual_rms),
            ("residual_fraction", self.residual_fraction),
        ):
            if not bool(torch.isfinite(tensor).all()):
                raise ValueError(f"{name} contains non-finite values.")
        if bool((self.normalization_scale <= 0).any()):
            raise ValueError("normalization_scale must be strictly positive.")
        if bool((self.residual_rms < 0).any()) or bool(
            (self.residual_fraction < 0).any()
        ):
            raise ValueError("residual diagnostics must be non-negative.")
        expected_pairs = torch.tensor(PATHOLOGY_PAIRS, dtype=torch.long)
        if not torch.equal(self.pairs.detach().cpu().long(), expected_pairs):
            raise ValueError("pairs do not match the canonical pathology pair order.")


def fit_pairwise_pathology_residualizer(
    pathology: torch.Tensor,
    pathology_valid: torch.Tensor,
    *,
    ridge: float = 1.0e-3,
    minimum_complete_donors: int = 24,
    scale_floor: float = 1.0e-4,
) -> PairwisePathologyResidualizerArtifact:
    """Fit pair residualizers on one row per training donor.

    For pair ``(a,b)`` the raw interaction is ``p_a * p_b``.  It is regressed
    on ``[hinge(p_a), hinge(p_b)]`` using complete donors only.  The
    standardized residual is consequently separated from the two named main
    effects on the fitting donors while remaining exactly zero at
    ``p_a=p_b=0``.
    """

    pathology = torch.as_tensor(pathology, dtype=torch.float64)
    pathology_valid = torch.as_tensor(pathology_valid, dtype=torch.bool)
    if pathology.shape != pathology_valid.shape:
        raise ValueError("pathology and pathology_valid must have the same shape.")
    if pathology.ndim != 2 or int(pathology.shape[1]) != len(PATHOLOGY_NAMES):
        raise ValueError(
            f"pathology must have shape [D,{len(PATHOLOGY_NAMES)}]."
        )
    if int(pathology.shape[0]) <= 0:
        raise ValueError("at least one training donor is required.")
    if not bool(torch.isfinite(pathology).all()):
        raise ValueError("pathology contains non-finite values.")
    if bool((pathology < 0).any()) or bool((pathology > 1).any()):
        raise ValueError("pathology must be normalized to [0,1].")
    if not torch.isfinite(torch.tensor(float(ridge))) or float(ridge) < 0.0:
        raise ValueError("ridge must be finite and non-negative.")
    if (
        isinstance(minimum_complete_donors, bool)
        or int(minimum_complete_donors) < 7
    ):
        raise ValueError("minimum_complete_donors must be an integer >= 7.")
    if (
        not torch.isfinite(torch.tensor(float(scale_floor)))
        or float(scale_floor) <= 0.0
    ):
        raise ValueError("scale_floor must be finite and positive.")

    hinges = _hinge_features(pathology)
    projection = torch.zeros(len(PATHOLOGY_PAIRS), 6, dtype=torch.float64)
    normalization_scale = torch.ones(len(PATHOLOGY_PAIRS), dtype=torch.float64)
    residual_rms = torch.zeros(len(PATHOLOGY_PAIRS), dtype=torch.float64)
    residual_fraction = torch.zeros(len(PATHOLOGY_PAIRS), dtype=torch.float64)
    eligible = torch.zeros(len(PATHOLOGY_PAIRS), dtype=torch.bool)
    complete_donors = torch.zeros(len(PATHOLOGY_PAIRS), dtype=torch.long)

    for pair_index, (left, right) in enumerate(PATHOLOGY_PAIRS):
        complete = pathology_valid[:, left] & pathology_valid[:, right]
        count = int(complete.sum().item())
        complete_donors[pair_index] = count
        if count < int(minimum_complete_donors):
            continue

        x = torch.cat(
            (
                hinges[complete, left],
                hinges[complete, right],
            ),
            dim=1,
        )
        target = pathology[complete, left] * pathology[complete, right]
        raw_scale = (target - target.mean()).square().mean().sqrt()
        if not bool(torch.isfinite(raw_scale)) or float(raw_scale) < float(scale_floor):
            continue
        penalty = torch.eye(6, dtype=pathology.dtype) * float(ridge)
        gram = x.T @ x + penalty
        rhs = x.T @ target
        try:
            beta = torch.linalg.solve(gram, rhs)
        except RuntimeError:
            # A zero-ridge diagnostic fit may be rank deficient when ordinal
            # pathology levels do not occupy every hinge segment.
            beta = torch.linalg.lstsq(gram, rhs.unsqueeze(-1)).solution.squeeze(-1)
        residual = target - x @ beta
        residual_scale = residual.square().mean().sqrt()
        if not bool(torch.isfinite(residual_scale)):
            continue

        projection[pair_index] = beta
        normalization_scale[pair_index] = raw_scale
        residual_rms[pair_index] = residual_scale
        residual_fraction[pair_index] = residual_scale / raw_scale
        eligible[pair_index] = True

    artifact = PairwisePathologyResidualizerArtifact(
        projection=projection.float(),
        normalization_scale=normalization_scale.float(),
        residual_rms=residual_rms.float(),
        residual_fraction=residual_fraction.float(),
        eligible=eligible,
        complete_donors=complete_donors,
        pairs=torch.tensor(PATHOLOGY_PAIRS, dtype=torch.long),
    )
    artifact.validate()
    return artifact


class PairwisePathologyInteractionBasis(nn.Module):
    """Apply a fixed train-donor residualizer to cell-level pathology rows."""

    def __init__(
        self,
        artifact: PairwisePathologyResidualizerArtifact,
        *,
        pair_names: Sequence[str] = PATHOLOGY_PAIR_NAMES,
    ) -> None:
        super().__init__()
        artifact.validate()
        if tuple(str(name) for name in pair_names) != PATHOLOGY_PAIR_NAMES:
            raise ValueError("pair_names must follow the canonical pathology pair order.")
        self.pair_names = tuple(str(name) for name in pair_names)
        self.register_buffer(
            "projection", artifact.projection.float().clone(), persistent=True
        )
        self.register_buffer(
            "normalization_scale",
            artifact.normalization_scale.float().clone(),
            persistent=True,
        )
        self.register_buffer(
            "residual_rms",
            artifact.residual_rms.float().clone(),
            persistent=True,
        )
        self.register_buffer(
            "residual_fraction",
            artifact.residual_fraction.float().clone(),
            persistent=True,
        )
        self.register_buffer(
            "eligible", artifact.eligible.bool().clone(), persistent=True
        )
        self.register_buffer(
            "complete_donors",
            artifact.complete_donors.long().clone(),
            persistent=True,
        )
        self.register_buffer("pairs", artifact.pairs.long().clone(), persistent=True)

    @property
    def n_pairs(self) -> int:
        return int(self.pairs.shape[0])

    def forward(
        self,
        pathology: torch.Tensor,
        pathology_valid: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        pathology = pathology.float()
        pathology_valid = pathology_valid.bool()
        if pathology.shape != pathology_valid.shape:
            raise ValueError("pathology and pathology_valid must have the same shape.")
        if pathology.ndim != 2 or int(pathology.shape[1]) != len(PATHOLOGY_NAMES):
            raise ValueError(
                f"pathology must have shape [B,{len(PATHOLOGY_NAMES)}]."
            )
        if not bool(torch.isfinite(pathology).all()):
            raise ValueError("pathology contains non-finite values.")

        hinges = _hinge_features(pathology)
        values = []
        valid_pairs = []
        for pair_index in range(self.n_pairs):
            left = int(self.pairs[pair_index, 0].item())
            right = int(self.pairs[pair_index, 1].item())
            design = torch.cat(
                (
                    hinges[:, left],
                    hinges[:, right],
                ),
                dim=1,
            )
            raw = pathology[:, left] * pathology[:, right]
            residual = raw - design @ self.projection[pair_index].to(
                device=pathology.device, dtype=pathology.dtype
            )
            # Normalize by the raw pair-product variation, not by the often
            # tiny residual variation.  Near-collinear pathology pairs then
            # remain weak instead of having measurement noise amplified to
            # unit variance.
            residual = residual / self.normalization_scale[pair_index].to(
                device=pathology.device, dtype=pathology.dtype
            )
            valid = (
                pathology_valid[:, left]
                & pathology_valid[:, right]
                & self.eligible[pair_index]
            )
            values.append(torch.where(valid, residual, torch.zeros_like(residual)))
            valid_pairs.append(valid)
        return torch.stack(values, dim=1), torch.stack(valid_pairs, dim=1)
