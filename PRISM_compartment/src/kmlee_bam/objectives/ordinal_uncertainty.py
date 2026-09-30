from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn.functional as F


def _zero_like(x: torch.Tensor) -> torch.Tensor:
    return torch.zeros((), dtype=x.dtype, device=x.device)


@dataclass(frozen=True)
class OrdinalBalanceConfig:
    enabled: bool = True
    lambda_balanced: float = 0.0
    lambda_nonzero: float = 0.0
    class_weight_mode: str = "sqrt_inverse"
    class_weight_cap: float = 6.0
    zero_bin_weight_floor: float = 0.2
    eps: float = 1e-6


@dataclass(frozen=True)
class OrdinalBalanceOutput:
    loss_balanced: torch.Tensor
    loss_nonzero: torch.Tensor
    balanced_recall: torch.Tensor
    nonzero_accuracy: torch.Tensor
    nonzero_within1: torch.Tensor
    true_zero_fraction: torch.Tensor
    pred_zero_fraction: torch.Tensor


def make_ordinal_class_weights(
    counts: torch.Tensor,
    *,
    mode: str = "sqrt_inverse",
    cap: float = 6.0,
    zero_bin_weight_floor: float = 0.2,
    eps: float = 1e-6,
) -> torch.Tensor:
    """Build stable ordinal-bin weights from train-level bin counts."""
    if counts.ndim != 1:
        raise ValueError(f"counts must have shape [K], got {tuple(counts.shape)}.")
    counts = counts.to(dtype=torch.float32)
    if int(counts.numel()) < 2:
        raise ValueError("ordinal class weights require at least two bins.")

    smoothed = counts.clamp_min(float(eps))
    weights = smoothed.sum() / (float(counts.numel()) * smoothed)
    mode_l = str(mode).lower()
    if mode_l == "sqrt_inverse":
        weights = torch.sqrt(weights)
    elif mode_l == "inverse":
        pass
    else:
        raise ValueError("class_weight_mode must be 'inverse' or 'sqrt_inverse'.")

    if cap is not None and float(cap) > 0:
        weights = weights.clamp(max=float(cap))
    if float(zero_bin_weight_floor) > 0:
        weights[0] = weights[0].clamp_min(float(zero_bin_weight_floor))

    # Keep the expected weight close to one under the observed train
    # distribution so the auxiliary loss is easy to tune.
    expected = (weights * (counts / counts.sum().clamp_min(1.0))).sum()
    return weights / expected.clamp_min(float(eps))


def ordinal_balance_losses(
    *,
    nll_per_gene: Optional[torch.Tensor],
    probs: Optional[torch.Tensor],
    y_ord: torch.Tensor,
    class_weights: torch.Tensor,
    config: OrdinalBalanceConfig,
) -> OrdinalBalanceOutput:
    """Auxiliary losses/diagnostics that do not let bin 0 dominate."""
    if nll_per_gene is None:
        z = _zero_like(y_ord.float())
        return OrdinalBalanceOutput(z, z, z, z, z, z, z)
    if y_ord.ndim != 2:
        raise ValueError(f"y_ord must have shape [B, G], got {tuple(y_ord.shape)}.")
    if nll_per_gene.shape != y_ord.shape:
        raise ValueError("nll_per_gene and y_ord must have matching [B, G] shapes.")

    weights = class_weights.to(device=nll_per_gene.device, dtype=nll_per_gene.dtype)
    gene_weights = weights[y_ord.long()]
    loss_balanced = (nll_per_gene * gene_weights).sum() / gene_weights.sum().clamp_min(
        float(config.eps)
    )

    nonzero_mask = y_ord > 0
    if bool(nonzero_mask.any()):
        loss_nonzero = nll_per_gene[nonzero_mask].mean()
    else:
        loss_nonzero = _zero_like(nll_per_gene)

    if probs is None:
        z = _zero_like(nll_per_gene)
        return OrdinalBalanceOutput(
            loss_balanced,
            loss_nonzero,
            z,
            z,
            z,
            z,
            z,
        )

    pred = probs.argmax(dim=-1)
    true = y_ord.long()
    n_bins = int(weights.numel())
    recalls = []
    for b in range(n_bins):
        mask = true == b
        if bool(mask.any()):
            recalls.append((pred[mask] == b).float().mean())
    balanced_recall = (
        torch.stack(recalls).mean()
        if recalls
        else _zero_like(nll_per_gene)
    )
    if bool(nonzero_mask.any()):
        nonzero_accuracy = (pred[nonzero_mask] == true[nonzero_mask]).float().mean()
        nonzero_within1 = (
            (pred[nonzero_mask] - true[nonzero_mask]).abs() <= 1
        ).float().mean()
    else:
        nonzero_accuracy = _zero_like(nll_per_gene)
        nonzero_within1 = _zero_like(nll_per_gene)

    true_zero_fraction = (true == 0).float().mean()
    pred_zero_fraction = (pred == 0).float().mean()
    return OrdinalBalanceOutput(
        loss_balanced=loss_balanced,
        loss_nonzero=loss_nonzero,
        balanced_recall=balanced_recall,
        nonzero_accuracy=nonzero_accuracy,
        nonzero_within1=nonzero_within1,
        true_zero_fraction=true_zero_fraction,
        pred_zero_fraction=pred_zero_fraction,
    )


@dataclass(frozen=True)
class GlobalReferenceBankConfig:
    enabled: bool = True
    lambda_global: float = 0.0
    ema_momentum: float = 0.95
    shrinkage_k: float = 8.0
    min_reference_cells: int = 1
    update_during_eval: bool = False


@dataclass(frozen=True)
class GlobalReferenceBankOutput:
    loss: torch.Tensor
    active_fraction: torch.Tensor
    center_norm_mean: torch.Tensor
    center_norm_max: torch.Tensor
    reliability_mean: torch.Tensor
    n_active_celltypes: int


class GlobalReferenceCenterBank:
    """
    EMA memory bank for reference-normal z_perp centers.

    The bank is not a model parameter. It is a slowly updated calibration table
    keyed by cell type. The training gradient still comes from the current
    batch, but the target is the accumulated reference origin instead of a
    noisy one-batch estimate.
    """

    def __init__(
        self,
        *,
        n_celltypes: int,
        d_z: int,
        config: GlobalReferenceBankConfig,
        device: torch.device,
        dtype: torch.dtype = torch.float32,
    ) -> None:
        self.n_celltypes = int(n_celltypes)
        self.d_z = int(d_z)
        self.config = config
        self.centers = torch.zeros(self.n_celltypes, self.d_z, device=device, dtype=dtype)
        self.support = torch.zeros(self.n_celltypes, device=device, dtype=dtype)
        self.initialized = torch.zeros(self.n_celltypes, device=device, dtype=torch.bool)

    @property
    def device(self) -> torch.device:
        return self.centers.device

    def to(self, device: torch.device, dtype: Optional[torch.dtype] = None) -> None:
        dtype = self.centers.dtype if dtype is None else dtype
        self.centers = self.centers.to(device=device, dtype=dtype)
        self.support = self.support.to(device=device, dtype=dtype)
        self.initialized = self.initialized.to(device=device)

    def shrunk_centers(self) -> tuple[torch.Tensor, torch.Tensor]:
        support = self.support.clamp_min(0.0)
        reliability = support / (support + float(self.config.shrinkage_k))
        return self.centers * reliability[:, None], reliability

    @torch.no_grad()
    def update(
        self,
        *,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        is_reference: Optional[torch.Tensor],
    ) -> None:
        if is_reference is None:
            return
        ref = is_reference.bool()
        if not bool(ref.any()):
            return
        z_det = z_perp.detach()
        for c in torch.unique(celltype_id[ref]):
            c_int = int(c.detach().cpu())
            mask = ref & (celltype_id == c)
            if int(mask.sum().item()) < int(self.config.min_reference_cells):
                continue
            mean_c = z_det[mask].mean(dim=0)
            if bool(self.initialized[c_int]):
                m = float(self.config.ema_momentum)
                self.centers[c_int].mul_(m).add_(mean_c, alpha=1.0 - m)
            else:
                self.centers[c_int].copy_(mean_c)
                self.initialized[c_int] = True
            self.support[c_int] += 1.0

    def loss(
        self,
        *,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        is_reference: Optional[torch.Tensor],
    ) -> GlobalReferenceBankOutput:
        if is_reference is None:
            z = _zero_like(z_perp)
            return GlobalReferenceBankOutput(z, z, z, z, z, 0)
        ref = is_reference.bool()
        if not bool(ref.any()):
            z = _zero_like(z_perp)
            return GlobalReferenceBankOutput(z, z, z, z, z, 0)

        centers, reliability = self.shrunk_centers()
        targets = centers[celltype_id.long()].detach()
        rel = reliability[celltype_id.long()].detach()
        active = ref & (rel > 0)
        if not bool(active.any()):
            z = _zero_like(z_perp)
            return GlobalReferenceBankOutput(z, z, z, z, z, 0)

        per_cell = (z_perp[active] - targets[active]).pow(2).mean(dim=-1)
        loss = (per_cell * rel[active]).sum() / rel[active].sum().clamp_min(1e-6)

        active_celltypes = torch.unique(celltype_id[active])
        active_centers = centers[active_celltypes.long()]
        norms = active_centers.norm(dim=-1)
        active_reliability = reliability[active_celltypes.long()]
        return GlobalReferenceBankOutput(
            loss=loss,
            active_fraction=active.float().mean(),
            center_norm_mean=norms.mean(),
            center_norm_max=norms.max(),
            reliability_mean=active_reliability.mean(),
            n_active_celltypes=int(active_celltypes.numel()),
        )


@dataclass(frozen=True)
class UncertaintyCalibrationConfig:
    enabled: bool = True
    lambda_rec_alignment: float = 0.0
    lambda_relative_rec: float = 0.0
    tau_relative: float = 0.1
    group_by: str = "celltype_tech_depth"
    n_depth_bins: int = 4
    min_group_size: int = 4
    clamp_z: float = 3.0
    eps: float = 1e-6


@dataclass(frozen=True)
class UncertaintyCalibrationOutput:
    rec_alignment: torch.Tensor
    relative_rec: torch.Tensor
    rel_unc_std: torch.Tensor
    rel_w_min: torch.Tensor
    rel_w_mean: torch.Tensor
    rel_w_max: torch.Tensor


def _depth_bins(depth: torch.Tensor, n_bins: int) -> torch.Tensor:
    if depth.ndim != 1:
        raise ValueError("depth must have shape [B].")
    if int(depth.numel()) == 0:
        return depth.long()
    n_bins = max(1, int(n_bins))
    if n_bins == 1:
        return torch.zeros_like(depth, dtype=torch.long)
    order = torch.argsort(depth)
    bins = torch.empty_like(order)
    ranks = torch.arange(depth.numel(), device=depth.device)
    bins[order] = torch.clamp((ranks * n_bins) // max(int(depth.numel()), 1), max=n_bins - 1)
    return bins.long()


def make_uncertainty_group_id(
    *,
    celltype_id: torch.Tensor,
    tech_id: torch.Tensor,
    depth_value: Optional[torch.Tensor],
    config: UncertaintyCalibrationConfig,
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


def group_zscore(
    values: torch.Tensor,
    group_id: torch.Tensor,
    *,
    min_group_size: int,
    eps: float,
) -> torch.Tensor:
    if values.ndim != 1 or group_id.ndim != 1:
        raise ValueError("values and group_id must have shape [B].")
    global_center = values.mean()
    global_scale = values.std(unbiased=False).clamp_min(float(eps))
    out = torch.empty_like(values)
    for g in torch.unique(group_id):
        mask = group_id == g
        if int(mask.sum().item()) < int(min_group_size):
            out[mask] = (values[mask] - global_center) / global_scale
            continue
        center = values[mask].mean()
        scale = values[mask].std(unbiased=False).clamp_min(float(eps))
        out[mask] = (values[mask] - center) / scale
    return out


def uncertainty_calibration_losses(
    *,
    raw_uncertainty: Optional[torch.Tensor],
    rec_per_cell: torch.Tensor,
    celltype_id: torch.Tensor,
    tech_id: torch.Tensor,
    depth_value: Optional[torch.Tensor],
    config: UncertaintyCalibrationConfig,
) -> UncertaintyCalibrationOutput:
    if raw_uncertainty is None:
        z = _zero_like(rec_per_cell)
        one = torch.ones((), dtype=rec_per_cell.dtype, device=rec_per_cell.device)
        return UncertaintyCalibrationOutput(z, z, z, one, one, one)

    group_id = make_uncertainty_group_id(
        celltype_id=celltype_id,
        tech_id=tech_id,
        depth_value=depth_value,
        config=config,
    )
    rel_unc = group_zscore(
        raw_uncertainty,
        group_id,
        min_group_size=int(config.min_group_size),
        eps=float(config.eps),
    ).clamp(min=-float(config.clamp_z), max=float(config.clamp_z))

    rec_target = group_zscore(
        rec_per_cell.detach(),
        group_id,
        min_group_size=int(config.min_group_size),
        eps=float(config.eps),
    ).clamp(min=-float(config.clamp_z), max=float(config.clamp_z))
    rec_alignment = F.mse_loss(rel_unc, rec_target)

    rel_weights = torch.exp(-float(config.tau_relative) * rel_unc.detach())
    rel_weights = rel_weights / rel_weights.mean().clamp_min(float(config.eps))
    relative_rec = (rel_weights * rec_per_cell).mean()
    return UncertaintyCalibrationOutput(
        rec_alignment=rec_alignment,
        relative_rec=relative_rec,
        rel_unc_std=rel_unc.detach().std(unbiased=False),
        rel_w_min=rel_weights.detach().min(),
        rel_w_mean=rel_weights.detach().mean(),
        rel_w_max=rel_weights.detach().max(),
    )

