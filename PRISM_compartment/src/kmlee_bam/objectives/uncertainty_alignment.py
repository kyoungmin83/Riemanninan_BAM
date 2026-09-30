"""
Uncertainty alignment loss — BAM head ↔ noise-score-derived u_total.

The model's BAM head (encoder.cell_uncertainty) is trained to regress against
the bank/ANCOVA-derived u_total. This is the supervision signal that turns
the otherwise-collapsing BAM head into a per-cell calibrated uncertainty
estimator.

Target is detached by default — we treat u_total as a *fixed-for-this-step*
supervision target, not a back-prop path into the bank machinery.

See: model_v7a_implementation_plan.md §4.2.7
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn.functional as F


@dataclass
class UncertaintyAlignmentConfig:
    enabled: bool = False
    lambda_alignment: float = 0.05
    detach_target: bool = True
    loss_form: str = "mse"            # "mse" | "huber" | "log_ratio"
    huber_delta: float = 1.0
    alignment_start_epoch: int = 1
    ramp_epochs: int = 3              # linear warmup from 0 -> 1 over these epochs
    u_total_cap: Optional[float] = None  # clamp u_total to this max; None = no cap
    # Optional pathology-aware reconstruction reweighting.  This uses the
    # detached PHU target by default, rather than an uncalibrated BAM output,
    # so the model cannot lower its loss merely by changing its own confidence.
    lambda_relative_rec: float = 0.0
    tau_relative: float = 0.25
    relative_start_epoch: int = 6
    relative_ramp_epochs: int = 5
    relative_weight_min: float = 0.67
    relative_weight_max: float = 1.50
    relative_weight_source: str = "target"  # "target" | "bam"
    detach_relative_source: bool = True
    normalize_relative_within_celltype: bool = False
    minimum_relative_ess_fraction: float = 0.0
    eps: float = 1e-6


@dataclass
class RelativeReconstructionOutput:
    delta: torch.Tensor
    weighted_rec: torch.Tensor
    unweighted_rec: torch.Tensor
    ramp: float
    weight_min: torch.Tensor
    weight_mean: torch.Tensor
    weight_max: torch.Tensor
    effective_sample_fraction: torch.Tensor
    bam_target_corr: torch.Tensor


def _compute_ramp(epoch: Optional[int], ramp_epochs: int) -> float:
    if ramp_epochs <= 0 or epoch is None:
        return 1.0
    return float(min(1.0, max(0.0, float(epoch) / float(ramp_epochs))))


def _compute_delayed_ramp(
    epoch: Optional[int],
    *,
    start_epoch: int,
    ramp_epochs: int,
) -> float:
    """Return 0 before ``start_epoch`` and then linearly ramp to 1."""
    if epoch is None:
        return 1.0
    start = max(1, int(start_epoch))
    if int(epoch) < start:
        return 0.0
    if int(ramp_epochs) <= 0:
        return 1.0
    return float(min(1.0, (int(epoch) - start + 1) / float(ramp_epochs)))


def _signed_corr(x: torch.Tensor, y: torch.Tensor, eps: float) -> torch.Tensor:
    x = x.reshape(-1).float()
    y = y.reshape(-1).float()
    if x.numel() < 2 or y.numel() != x.numel():
        return x.new_zeros(())
    xz = (x - x.mean()) / x.std(unbiased=False).clamp_min(float(eps))
    yz = (y - y.mean()) / y.std(unbiased=False).clamp_min(float(eps))
    return (xz * yz).mean()


def compute_pathology_aware_relative_reconstruction(
    bam_head_output: torch.Tensor,
    u_total: torch.Tensor,
    rec_per_cell: torch.Tensor,
    config: UncertaintyAlignmentConfig,
    *,
    epoch: Optional[int] = None,
    group_id: Optional[torch.Tensor] = None,
) -> RelativeReconstructionOutput:
    """Build bounded, mean-one PHU weights and a centred reconstruction term.

    ``delta`` is ``weighted_rec - unweighted_rec``.  Adding
    ``lambda_relative_rec * delta`` changes the relative contribution of cells
    without adding a second copy of the batch-average reconstruction loss.
    """
    rec = rec_per_cell.reshape(-1).float()
    target = u_total.reshape(-1).float()
    bam = bam_head_output.reshape(-1).float()
    if rec.numel() != target.numel() or rec.numel() != bam.numel():
        raise ValueError(
            "rec_per_cell, u_total, and bam_head_output must contain the same "
            f"number of cells; got {rec.numel()}, {target.numel()}, {bam.numel()}"
        )

    source_name = str(config.relative_weight_source).lower()
    if source_name == "target":
        source = target
    elif source_name == "bam":
        source = bam
    else:
        raise ValueError(
            "relative_weight_source must be 'target' or 'bam', got "
            f"{config.relative_weight_source!r}"
        )
    if bool(config.detach_relative_source):
        source = source.detach()

    weights = torch.exp(-float(config.tau_relative) * source)

    def normalize(value: torch.Tensor) -> torch.Tensor:
        if not bool(config.normalize_relative_within_celltype):
            return value / value.mean().clamp_min(float(config.eps))
        if group_id is None:
            raise ValueError(
                "normalize_relative_within_celltype requires group_id"
            )
        group = group_id.reshape(-1).long().to(value.device)
        if group.numel() != value.numel():
            raise ValueError("relative reconstruction group_id shape mismatch")
        result = value.clone()
        for current in torch.unique(group, sorted=True):
            rows = group == current
            result[rows] = value[rows] / value[rows].mean().clamp_min(
                float(config.eps)
            )
        return result

    weights = normalize(weights)
    lo = float(config.relative_weight_min)
    hi = float(config.relative_weight_max)
    if lo <= 0.0 or hi < lo:
        raise ValueError(
            "relative weight bounds must satisfy 0 < min <= max; got "
            f"{lo}, {hi}"
        )
    weights = weights.clamp(min=lo, max=hi)
    weights = normalize(weights)
    minimum_ess = float(config.minimum_relative_ess_fraction)
    if not 0.0 <= minimum_ess <= 1.0:
        raise ValueError("minimum_relative_ess_fraction must lie in [0, 1]")
    if minimum_ess > 0.0:
        ess = weights.mean().pow(2) / weights.pow(2).mean().clamp_min(
            float(config.eps)
        )
        if float(ess.detach()) < minimum_ess:
            low, high = 0.0, 1.0
            for _ in range(24):
                middle = 0.5 * (low + high)
                candidate = 1.0 + middle * (weights - 1.0)
                candidate_ess = candidate.mean().pow(2) / candidate.pow(2).mean().clamp_min(
                    float(config.eps)
                )
                if float(candidate_ess.detach()) >= minimum_ess:
                    low = middle
                else:
                    high = middle
            weights = 1.0 + low * (weights - 1.0)

    unweighted_rec = rec.mean()
    weighted_rec = (weights * rec).mean()
    delta = weighted_rec - unweighted_rec
    ess_fraction = weights.mean().pow(2) / weights.pow(2).mean().clamp_min(
        float(config.eps)
    )
    ramp = _compute_delayed_ramp(
        epoch,
        start_epoch=int(config.relative_start_epoch),
        ramp_epochs=int(config.relative_ramp_epochs),
    )
    return RelativeReconstructionOutput(
        delta=delta,
        weighted_rec=weighted_rec,
        unweighted_rec=unweighted_rec,
        ramp=ramp,
        weight_min=weights.detach().min(),
        weight_mean=weights.detach().mean(),
        weight_max=weights.detach().max(),
        effective_sample_fraction=ess_fraction.detach(),
        bam_target_corr=_signed_corr(bam.detach(), target.detach(), float(config.eps)).detach(),
    )


def compute_uncertainty_alignment_loss(
    bam_head_output: torch.Tensor,    # [B] — model's predicted per-cell uncertainty
    u_total: torch.Tensor,             # [B] — bank/ANCOVA-derived target
    config: UncertaintyAlignmentConfig,
    *,
    epoch: Optional[int] = None,
) -> torch.Tensor:
    """
    Returns a scalar weighted alignment loss.

    The caller is responsible for the lambda weighting *unless* it wants the
    ramp_scale applied here. To keep this function composable with the rest
    of total_loss.py (which already aggregates λ-weighted scalars), we
    return the unweighted-by-lambda loss; the trainer multiplies by lambda
    when adding into total.
    """
    if not config.enabled:
        return torch.zeros((), dtype=bam_head_output.dtype, device=bam_head_output.device)

    target = u_total.detach() if config.detach_target else u_total

    pred = bam_head_output
    if pred.shape != target.shape:
        # Best-effort reshape: BAM head typically returns [B] or [B, 1]; we
        # squeeze trailing singleton dims to align.
        if pred.dim() == target.dim() + 1 and pred.shape[-1] == 1:
            pred = pred.squeeze(-1)
        elif target.dim() == pred.dim() + 1 and target.shape[-1] == 1:
            target = target.squeeze(-1)
        else:
            raise ValueError(
                f"Shape mismatch: pred {tuple(pred.shape)} vs target "
                f"{tuple(target.shape)}"
            )

    form = str(config.loss_form).lower()
    if form == "mse":
        loss = F.mse_loss(pred, target)
    elif form == "huber":
        loss = F.huber_loss(pred, target, delta=float(config.huber_delta))
    elif form == "log_ratio":
        # log(target / pred) — both must be positive
        eps = 1e-6
        log_ratio = (
            torch.log(pred.clamp_min(eps)) - torch.log(target.clamp_min(eps))
        )
        loss = (log_ratio ** 2).mean()
    else:
        raise ValueError(f"Unknown loss_form: {config.loss_form}")

    ramp = _compute_delayed_ramp(
        epoch,
        start_epoch=int(config.alignment_start_epoch),
        ramp_epochs=int(config.ramp_epochs),
    )
    return ramp * loss
