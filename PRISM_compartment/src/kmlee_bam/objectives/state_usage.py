from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class V5StateUsageConfig:
    enabled: bool = True
    lambda_state_abs: float = 0.0
    target_state_abs: float = 0.0
    lambda_state_fraction: float = 0.0
    target_state_fraction: float = 0.0
    lambda_tech_to_state: float = 0.0
    max_tech_to_state: float = 0.0
    eps: float = 1e-6


@dataclass(frozen=True)
class V5StateUsageOutput:
    loss_state_abs: torch.Tensor
    loss_state_fraction: torch.Tensor
    loss_tech_to_state: torch.Tensor
    state_abs: torch.Tensor
    base_abs: torch.Tensor
    tech_abs: torch.Tensor
    state_fraction: torch.Tensor
    tech_to_state: torch.Tensor
    # v19 biological-covariate factorization: the decoder score is now
    # base + tech + SEX + state, so the score "budget" includes a sex term.
    # These are 0 (and state_fraction is byte-identical to the pre-v19 value)
    # whenever the decoder produces no sex_score (n_sex unset / sex_id absent).
    sex_abs: torch.Tensor
    sex_to_state: torch.Tensor
    precision_abs: torch.Tensor
    # The 3-term state_fraction (base+tech+state only), kept for continuity with
    # pre-v19 logs / docs (§12.7) so the value can be compared across versions.
    state_fraction_no_sex: torch.Tensor


def v5_state_usage_losses(
    *,
    decoder_out,
    config: V5StateUsageConfig,
) -> V5StateUsageOutput:
    state_score = decoder_out.state_score
    base_score = decoder_out.base_score
    tech_score = decoder_out.tech_score
    # v19: include the SEX baseline in the score "budget" so state_fraction is the
    # share of the FULL decoder score, not just base+tech+state. Backward-compatible:
    # when the decoder has no sex term, sex_score is None ⇒ sex_abs = 0 and every
    # value below is byte-identical to the pre-v19 behaviour.
    sex_score = getattr(decoder_out, "sex_score", None)
    precision_score = getattr(decoder_out, "precision_score", None)

    state_abs = state_score.abs().mean()
    base_abs = base_score.abs().mean()
    tech_abs = tech_score.abs().mean()
    if sex_score is not None:
        sex_abs = sex_score.abs().mean()
    else:
        sex_abs = torch.zeros((), dtype=state_abs.dtype, device=state_abs.device)
    if precision_score is not None:
        precision_abs = precision_score.abs().mean()
    else:
        precision_abs = torch.zeros((), dtype=state_abs.dtype, device=state_abs.device)
    eps = float(config.eps)
    # Full score-budget denominator.  PRISM is an explicit decoder component,
    # so omitting it would overstate the residual state's share after the new
    # architecture is attached.  state_fraction_no_sex remains the historical
    # base+tech+state diagnostic for cross-version comparisons.
    denom = state_abs + base_abs + tech_abs + sex_abs + precision_abs + eps
    denom_no_sex = state_abs + base_abs + tech_abs + eps
    state_fraction = state_abs / denom
    state_fraction_no_sex = state_abs / denom_no_sex
    tech_to_state = tech_abs / (state_abs + eps)
    sex_to_state = sex_abs / (state_abs + eps)

    loss_state_abs = F.relu(float(config.target_state_abs) - state_abs).pow(2)
    loss_state_fraction = F.relu(
        float(config.target_state_fraction) - state_fraction
    ).pow(2)
    if float(config.max_tech_to_state) > 0:
        loss_tech_to_state = F.relu(
            tech_to_state - float(config.max_tech_to_state)
        ).pow(2)
    else:
        loss_tech_to_state = torch.zeros((), dtype=state_abs.dtype, device=state_abs.device)

    return V5StateUsageOutput(
        loss_state_abs=loss_state_abs,
        loss_state_fraction=loss_state_fraction,
        loss_tech_to_state=loss_tech_to_state,
        state_abs=state_abs.detach(),
        base_abs=base_abs.detach(),
        tech_abs=tech_abs.detach(),
        state_fraction=state_fraction.detach(),
        tech_to_state=tech_to_state.detach(),
        sex_abs=sex_abs.detach(),
        sex_to_state=sex_to_state.detach(),
        precision_abs=precision_abs.detach(),
        state_fraction_no_sex=state_fraction_no_sex.detach(),
    )
