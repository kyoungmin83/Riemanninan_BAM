from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn.functional as F


def _zero_like(x: torch.Tensor) -> torch.Tensor:
    return torch.zeros((), dtype=x.dtype, device=x.device)


@dataclass(frozen=True)
class V5UncertaintyResidualConfig:
    enabled: bool = True
    lambda_relative_rec: float = 0.0
    lambda_rec_alignment: float = 0.0
    lambda_depth_corr: float = 0.0
    lambda_saturation: float = 0.0
    tau_relative: float = 0.1
    group_by: str = "celltype_tech_depth"
    n_depth_bins: int = 4
    min_group_size: int = 8
    clamp_z: float = 2.5
    raw_std_floor: float = 0.0
    detach_weight_uncertainty: bool = True
    eps: float = 1e-6


@dataclass(frozen=True)
class V5UncertaintyResidualOutput:
    relative_rec: torch.Tensor
    rec_alignment: torch.Tensor
    depth_corr: torch.Tensor
    saturation: torch.Tensor
    rel_unc_std: torch.Tensor
    raw_unc_std: torch.Tensor
    rel_w_min: torch.Tensor
    rel_w_mean: torch.Tensor
    rel_w_max: torch.Tensor


def _depth_bins(depth: torch.Tensor, n_bins: int) -> torch.Tensor:
    n_bins = max(1, int(n_bins))
    if n_bins == 1 or int(depth.numel()) == 0:
        return torch.zeros_like(depth, dtype=torch.long)
    order = torch.argsort(depth)
    bins = torch.empty_like(order)
    ranks = torch.arange(depth.numel(), device=depth.device)
    bins[order] = torch.clamp((ranks * n_bins) // max(int(depth.numel()), 1), max=n_bins - 1)
    return bins.long()


def _group_id(
    *,
    celltype_id: torch.Tensor,
    tech_id: torch.Tensor,
    depth_value: Optional[torch.Tensor],
    config: V5UncertaintyResidualConfig,
) -> torch.Tensor:
    mode = str(config.group_by).lower()
    if mode == "celltype":
        return celltype_id.long()
    if mode == "celltype_tech":
        return celltype_id.long() * (int(tech_id.max().item()) + 1) + tech_id.long()
    if mode == "celltype_tech_depth" and depth_value is not None:
        depth_bin = _depth_bins(depth_value.float(), int(config.n_depth_bins))
        n_tech = int(tech_id.max().item()) + 1
        return (celltype_id.long() * n_tech + tech_id.long()) * int(config.n_depth_bins) + depth_bin
    return celltype_id.long()


def _group_zscore(
    values: torch.Tensor,
    group_id: torch.Tensor,
    *,
    min_group_size: int,
    eps: float,
) -> torch.Tensor:
    global_center = values.mean()
    global_scale = values.std(unbiased=False).clamp_min(float(eps))
    out = torch.empty_like(values)
    for g in torch.unique(group_id):
        m = group_id == g
        if int(m.sum().item()) < int(min_group_size):
            out[m] = (values[m] - global_center) / global_scale
            continue
        center = values[m].mean()
        scale = values[m].std(unbiased=False).clamp_min(float(eps))
        out[m] = (values[m] - center) / scale
    return out


def _corr_squared(x: torch.Tensor, y: torch.Tensor, eps: float) -> torch.Tensor:
    if x.numel() < 2:
        return _zero_like(x)
    xz = (x - x.mean()) / x.std(unbiased=False).clamp_min(float(eps))
    yz = (y - y.mean()) / y.std(unbiased=False).clamp_min(float(eps))
    return (xz * yz).mean().pow(2)


def v5_uncertainty_residual_losses(
    *,
    raw_uncertainty: Optional[torch.Tensor],
    rec_per_cell: torch.Tensor,
    celltype_id: torch.Tensor,
    tech_id: torch.Tensor,
    depth_value: Optional[torch.Tensor],
    config: V5UncertaintyResidualConfig,
) -> V5UncertaintyResidualOutput:
    if raw_uncertainty is None:
        z = _zero_like(rec_per_cell)
        one = torch.ones((), dtype=rec_per_cell.dtype, device=rec_per_cell.device)
        return V5UncertaintyResidualOutput(z, z, z, z, z, z, one, one, one)

    raw = raw_uncertainty.view(-1).float()
    rec = rec_per_cell.view(-1).float()
    group_id = _group_id(
        celltype_id=celltype_id,
        tech_id=tech_id,
        depth_value=depth_value,
        config=config,
    )
    rel_unc = _group_zscore(
        raw,
        group_id,
        min_group_size=int(config.min_group_size),
        eps=float(config.eps),
    ).clamp(min=-float(config.clamp_z), max=float(config.clamp_z))

    rec_z = _group_zscore(
        rec.detach(),
        group_id,
        min_group_size=int(config.min_group_size),
        eps=float(config.eps),
    ).clamp(min=-float(config.clamp_z), max=float(config.clamp_z))
    rec_alignment = F.mse_loss(rel_unc, rec_z)

    rel_for_weight = rel_unc.detach() if bool(config.detach_weight_uncertainty) else rel_unc
    rel_weights = torch.exp(-float(config.tau_relative) * rel_for_weight)
    rel_weights = rel_weights / rel_weights.mean().clamp_min(float(config.eps))
    relative_rec = (rel_weights * rec).mean()

    if depth_value is not None:
        depth_z = _group_zscore(
            depth_value.float().view(-1).detach(),
            group_id,
            min_group_size=int(config.min_group_size),
            eps=float(config.eps),
        )
        depth_corr = _corr_squared(rel_unc, depth_z, float(config.eps))
    else:
        depth_corr = _zero_like(rec)

    raw_unc_std = raw.std(unbiased=False)
    saturation = F.relu(float(config.raw_std_floor) - raw_unc_std).pow(2)

    return V5UncertaintyResidualOutput(
        relative_rec=relative_rec,
        rec_alignment=rec_alignment,
        depth_corr=depth_corr,
        saturation=saturation,
        rel_unc_std=rel_unc.detach().std(unbiased=False),
        raw_unc_std=raw_unc_std.detach(),
        rel_w_min=rel_weights.detach().min(),
        rel_w_mean=rel_weights.detach().mean(),
        rel_w_max=rel_weights.detach().max(),
    )
