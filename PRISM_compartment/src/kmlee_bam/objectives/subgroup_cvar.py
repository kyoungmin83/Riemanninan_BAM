"""
Adaptive subgroup CVaR: monitor-first, then auto-intervene.

Motivation
----------
Once we start relaxing zeros (π_tech soft-zero), the danger is *over-imputation*
in weak subgroups — rare cell types, shallow cells, certain batches — where the
average metric stays fine while those groups quietly degrade. CVaR ("care about
the worst, not the average") is the seatbelt.

Lifecycle (all DDP-safe; the only collective is one all_reduce per epoch):
  1. **Monitor (always, from epoch 1).** Each epoch, accumulate per-subgroup
     mean of a chosen per-cell metric (reconstruction NLL by default) and log it
     as `metric/cvar_*`. No effect on the loss.
  2. **Trigger.** After `warmup_epochs`, a subgroup whose EMA metric exceeds
     `median + mad_k·MAD` of the cohort is flagged "collapsed", with two-sided
     hysteresis (`persist_epochs` to engage, `release_epochs` to disengage) so
     interventions don't flap.
  3. **Intervene.** Engaged subgroups get a bounded, ramped per-cell loss
     multiplier in [1, w_max], mean-normalised to ≈1 (so total gradient scale is
     preserved). Applied by the trainer in *delta form*
     `λ·(mean(w·nll) − stop_grad(mean(nll)))` so it is exactly 0 when nothing is
     engaged and never double-counts the base reconstruction.

Subgroups reuse the existing `(celltype × tech × depth_bin)` scheme, but with a
FIXED group count (`n_celltypes · n_tech_max · n_depth_bins`) so the weight
table is identical across ranks and epochs.

See doc/zero_origin_and_capacity_design_2026-06-01.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.distributed as dist

from kmlee_bam.objectives.uncertainty_residual import _depth_bins


@dataclass(frozen=True)
class SubgroupCVaRConfig:
    enabled: bool = False
    # grouping
    group_by: str = "celltype_tech_depth"   # celltype | celltype_tech | celltype_tech_depth
    n_depth_bins: int = 4
    metric: str = "nll"                      # "nll" | "nonzero_err"
    # monitor
    ema_momentum: float = 0.9
    min_group_cells: int = 50
    warmup_epochs: int = 3
    # trigger
    mad_k: float = 3.0
    worst_alpha: float = 0.10                # CVaR tail fraction (mode="cvar")
    persist_epochs: int = 2
    release_epochs: int = 3
    # intervention
    mode: str = "reweight"                   # "reweight" (default) | "cvar"
    w_max: float = 4.0
    ramp_epochs: int = 3
    lambda_cvar: float = 1.0                 # weight on the delta term in total loss
    eps: float = 1e-6


def group_id_global(
    celltype_id: torch.Tensor,
    tech_id: torch.Tensor,
    depth: Optional[torch.Tensor],
    *,
    n_tech_max: int,
    n_depth_bins: int,
    group_by: str,
) -> torch.Tensor:
    """Composite subgroup id with a FIXED layout (constant `n_tech_max` /
    `n_depth_bins`), so ids are comparable across batches and ranks — unlike
    `uncertainty_residual._group_id`, which derives `n_tech` per batch."""
    ct = celltype_id.long()
    tech = tech_id.long().clamp(min=0, max=int(n_tech_max) - 1)
    mode = str(group_by).lower()
    if mode == "celltype":
        return ct
    if mode == "celltype_tech":
        return ct * int(n_tech_max) + tech
    if mode == "celltype_tech_depth" and depth is not None:
        db = _depth_bins(depth.float(), int(n_depth_bins))           # [B] in [0, n_depth_bins)
        return (ct * int(n_tech_max) + tech) * int(n_depth_bins) + db
    # fallback (no depth available): collapse depth axis to bin 0
    return (ct * int(n_tech_max) + tech) * int(n_depth_bins)


def per_cell_metric(
    nll_per_cell: torch.Tensor,
    probs: torch.Tensor,
    y_ord: torch.Tensor,
    *,
    metric: str,
) -> torch.Tensor:
    """Per-cell scalar [B], higher = worse (same orientation for both modes)."""
    if metric == "nonzero_err":
        p_nz = 1.0 - probs[..., 0]
        pred_nz = (p_nz > 0.5)
        true_nz = (y_ord > 0)
        return (pred_nz != true_nz).to(dtype=torch.float32).mean(dim=1)
    return nll_per_cell.detach().to(dtype=torch.float32)


class AdaptiveSubgroupCVaR:
    """Stateful per-epoch monitor + bounded intervention (see module docstring)."""

    def __init__(
        self,
        *,
        n_celltypes: int,
        n_tech_max: int,
        config: SubgroupCVaRConfig,
        device: torch.device,
    ) -> None:
        self.cfg = config
        self.n_celltypes = int(n_celltypes)
        self.n_tech_max = int(n_tech_max)
        self.n_depth_bins = max(1, int(config.n_depth_bins))
        if str(config.group_by).lower() == "celltype":
            self.n_groups = self.n_celltypes
        elif str(config.group_by).lower() == "celltype_tech":
            self.n_groups = self.n_celltypes * self.n_tech_max
        else:
            self.n_groups = self.n_celltypes * self.n_tech_max * self.n_depth_bins
        self.device = device
        self.sync_group = None
        self._epoch = 0

        z = lambda dt=torch.float32: torch.zeros(self.n_groups, device=device, dtype=dt)
        self.ema_metric = z()
        self.ema_count = z()
        self.weight_table = torch.ones(self.n_groups, device=device, dtype=torch.float32)
        self.consec_bad = z(torch.long)
        self.consec_good = z(torch.long)
        self.triggered = torch.zeros(self.n_groups, device=device, dtype=torch.bool)
        self.trigger_epoch = z(torch.long)
        self._epoch_sum = z()
        self._epoch_cnt = z()

    # --------------------------------------------------------------
    def set_sync_process_group(self, group) -> None:
        self.sync_group = group

    def group_ids(self, celltype_id, tech_id, depth):
        return group_id_global(
            celltype_id, tech_id, depth,
            n_tech_max=self.n_tech_max,
            n_depth_bins=self.n_depth_bins,
            group_by=self.cfg.group_by,
        ).clamp(min=0, max=self.n_groups - 1)

    def begin_epoch(self, epoch_index: Optional[int]) -> None:
        self._epoch = int(epoch_index) if epoch_index is not None else (self._epoch + 1)
        self._epoch_sum.zero_()
        self._epoch_cnt.zero_()

    @torch.no_grad()
    def observe_batch(self, group_id: torch.Tensor, metric_per_cell: torch.Tensor) -> None:
        gid = group_id.long().clamp(0, self.n_groups - 1)
        self._epoch_sum.scatter_add_(0, gid, metric_per_cell.detach().to(torch.float32))
        self._epoch_cnt.scatter_add_(0, gid, torch.ones_like(gid, dtype=torch.float32))

    def is_active(self) -> bool:
        return bool(self.cfg.enabled) and self._epoch > int(self.cfg.warmup_epochs) and bool(self.triggered.any())

    @torch.no_grad()
    def weights_for(self, group_id: torch.Tensor) -> torch.Tensor:
        """Per-cell multiplier, mean-normalised to ≈1 over the batch."""
        gid = group_id.long().clamp(0, self.n_groups - 1)
        w = self.weight_table.index_select(0, gid)
        return w / w.mean().clamp_min(float(self.cfg.eps))

    # --------------------------------------------------------------
    @torch.no_grad()
    def end_epoch_sync_and_update(self, epoch_index: Optional[int]) -> None:
        epoch = int(epoch_index) if epoch_index is not None else self._epoch
        if dist.is_available() and dist.is_initialized():
            stacked = torch.stack([self._epoch_sum, self._epoch_cnt], dim=0)
            dist.all_reduce(stacked, op=dist.ReduceOp.SUM, group=self.sync_group)
            self._epoch_sum, self._epoch_cnt = stacked[0].clone(), stacked[1].clone()

        valid = self._epoch_cnt >= float(self.cfg.min_group_cells)
        epoch_mean = self._epoch_sum / self._epoch_cnt.clamp_min(1.0)
        m = float(self.cfg.ema_momentum)
        fresh = valid & (self.ema_count == 0)
        ema = valid & (self.ema_count > 0)
        self.ema_metric[fresh] = epoch_mean[fresh]
        self.ema_metric[ema] = m * self.ema_metric[ema] + (1.0 - m) * epoch_mean[ema]
        self.ema_count[valid] += self._epoch_cnt[valid]

        p_g = self._epoch_cnt / self._epoch_cnt.sum().clamp_min(1.0)
        self._rebuild_weight_table(epoch, valid, p_g)

    @torch.no_grad()
    def _rebuild_weight_table(self, epoch, valid, p_g) -> None:
        L = self.ema_metric
        pool = valid & (self.ema_count > 0)
        if not bool(pool.any()):
            self.weight_table = torch.ones_like(self.weight_table)
            return

        med = torch.median(L[pool])
        mad = torch.median((L[pool] - med).abs()).clamp_min(float(self.cfg.eps))
        severity = ((L - med) / mad).clamp_min(0.0)
        collapsed = pool & (severity > float(self.cfg.mad_k))

        self.consec_bad = torch.where(collapsed, self.consec_bad + 1, torch.zeros_like(self.consec_bad))
        self.consec_good = torch.where(pool & ~collapsed, self.consec_good + 1, torch.zeros_like(self.consec_good))
        turn_on = (~self.triggered) & (self.consec_bad >= int(self.cfg.persist_epochs))
        turn_off = self.triggered & (self.consec_good >= int(self.cfg.release_epochs))
        self.triggered = (self.triggered | turn_on) & ~turn_off
        self.trigger_epoch = torch.where(turn_on, torch.full_like(self.trigger_epoch, int(epoch)), self.trigger_epoch)

        if epoch <= int(self.cfg.warmup_epochs):
            self.weight_table = torch.ones_like(self.weight_table)
            return

        if str(self.cfg.mode).lower() == "cvar":
            w_raw = self._cvar_weights(L, p_g, pool)
        else:
            w_raw = 1.0 + (float(self.cfg.w_max) - 1.0) * (severity / float(self.cfg.mad_k)).clamp(0.0, 1.0)
        w_raw = w_raw.clamp(max=float(self.cfg.w_max))

        ramp = ((float(epoch) - self.trigger_epoch.float()).clamp(min=0.0) / max(1, int(self.cfg.ramp_epochs))).clamp(max=1.0)
        self.weight_table = torch.where(
            self.triggered,
            1.0 + ramp * (w_raw - 1.0),
            torch.ones_like(self.weight_table),
        )

    @torch.no_grad()
    def _cvar_weights(self, L, p_g, pool) -> torch.Tensor:
        """Hard CVaR tail weights: groups in the worst-α mass get p_g/α, else 1."""
        w = torch.ones_like(L)
        idx = torch.nonzero(pool, as_tuple=False).view(-1)
        if idx.numel() == 0:
            return w
        order = idx[torch.argsort(L[idx], descending=True)]
        cum = torch.cumsum(p_g[order], dim=0)
        in_tail = cum <= float(self.cfg.worst_alpha)
        if not bool(in_tail.any()):
            in_tail[0] = True  # always include the single worst group
        tail = order[in_tail]
        w[tail] = (p_g[tail] / max(float(self.cfg.worst_alpha), float(self.cfg.eps))).clamp(min=1.0)
        return w

    # --------------------------------------------------------------
    def state_dict(self) -> dict:
        return {
            "ema_metric": self.ema_metric.detach().cpu(),
            "ema_count": self.ema_count.detach().cpu(),
            "weight_table": self.weight_table.detach().cpu(),
            "consec_bad": self.consec_bad.detach().cpu(),
            "consec_good": self.consec_good.detach().cpu(),
            "triggered": self.triggered.detach().cpu(),
            "trigger_epoch": self.trigger_epoch.detach().cpu(),
            "epoch": int(self._epoch),
        }

    def load_state_dict(self, state: dict) -> None:
        if not state:
            return
        dev = self.device
        self.ema_metric = state["ema_metric"].to(dev)
        self.ema_count = state["ema_count"].to(dev)
        self.weight_table = state["weight_table"].to(dev)
        self.consec_bad = state["consec_bad"].to(dev)
        self.consec_good = state["consec_good"].to(dev)
        self.triggered = state["triggered"].to(dev)
        self.trigger_epoch = state["trigger_epoch"].to(dev)
        self._epoch = int(state.get("epoch", 0))
