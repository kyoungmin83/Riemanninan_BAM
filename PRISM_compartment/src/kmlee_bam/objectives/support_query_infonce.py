"""Matched support-to-query contrastive objective for PRISM personal codes."""

from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class SupportQueryInfoNCEOutput:
    loss: torch.Tensor
    contrastive: torch.Tensor
    anchor: torch.Tensor
    valid_fraction: torch.Tensor
    own_nll: torch.Tensor
    negative_nll: torch.Tensor


def anchored_support_query_infonce(
    candidate_nll: torch.Tensor,
    negative_valid: torch.Tensor,
    common_nll: torch.Tensor,
    *,
    temperature: float,
    anchor_weight: float,
    cell_weight: torch.Tensor | None = None,
) -> SupportQueryInfoNCEOutput:
    """Prefer the own-donor support code without sacrificing common NLL.

    Candidate column zero is the query donor's own support code.  Remaining
    columns are covariate/support-composition matched other-donor codes.  The
    common-only anchor prevents the contrastive term from winning merely by
    making every personal prediction worse.
    """

    if candidate_nll.ndim != 2 or candidate_nll.shape[1] < 2:
        raise ValueError("candidate_nll must have shape [B,1+K] with K>=1")
    batch, n_candidates = candidate_nll.shape
    if negative_valid.shape != (batch, n_candidates - 1):
        raise ValueError("negative_valid shape does not match candidate_nll")
    if common_nll.shape != (batch,):
        raise ValueError("common_nll must have shape [B]")
    if not torch.isfinite(candidate_nll).all() or not torch.isfinite(common_nll).all():
        raise ValueError("support-query NLL inputs must be finite")
    if float(temperature) <= 0.0:
        raise ValueError("temperature must be positive")
    if float(anchor_weight) < 0.0:
        raise ValueError("anchor_weight must be non-negative")

    row_valid = negative_valid.any(dim=1)
    logits = -candidate_nll.float() / float(temperature)
    logits[:, 1:] = logits[:, 1:].masked_fill(~negative_valid, -torch.inf)
    target = torch.zeros(batch, dtype=torch.long, device=candidate_nll.device)
    per_row = F.cross_entropy(logits, target, reduction="none")
    anchor_per_row = F.relu(candidate_nll[:, 0].float() - common_nll.float())

    if cell_weight is None:
        weight = torch.ones(batch, device=candidate_nll.device, dtype=torch.float32)
    else:
        if cell_weight.shape != (batch,):
            raise ValueError("cell_weight must have shape [B]")
        weight = cell_weight.float().clamp_min(0.0)
    valid_weight = weight * row_valid.to(weight.dtype)
    denom = valid_weight.sum().clamp_min(1.0e-8)
    contrastive = (per_row * valid_weight).sum() / denom
    anchor = (anchor_per_row * valid_weight).sum() / denom
    loss = contrastive + float(anchor_weight) * anchor

    negative_values = candidate_nll[:, 1:]
    negative_weight = negative_valid.to(negative_values.dtype)
    negative_mean = (
        (negative_values * negative_weight).sum()
        / negative_weight.sum().clamp_min(1.0)
    )
    own_mean = (
        (candidate_nll[:, 0] * valid_weight.to(candidate_nll.dtype)).sum()
        / denom.to(candidate_nll.dtype)
    )
    return SupportQueryInfoNCEOutput(
        loss=loss,
        contrastive=contrastive,
        anchor=anchor,
        valid_fraction=row_valid.float().mean(),
        own_nll=own_mean,
        negative_nll=negative_mean,
    )
