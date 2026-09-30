"""
NoiseScoreBank — per (donor, celltype) EMA running stats of noise_score.

This is the L1 reference for cell-level uncertainty:
    u_cell(c) = (noise_score(c) - bank.mean[d, t]) / bank.std[d, t]

Updates happen per-batch via EMA; an `all_reduce_` call once per epoch keeps
DDP ranks in sync.

See: model_v7a_implementation_plan.md §4.2.4
"""

from __future__ import annotations

from typing import Optional, Tuple

import torch
import torch.distributed as dist


class NoiseScoreBank:
    """
    Running EMA stats per (donor, celltype).

    Storage:
        mean: [n_donors, n_celltypes]
        m2  : [n_donors, n_celltypes]  (sum of squared deviations / EMA proxy)
        n   : [n_donors, n_celltypes]  (EMA "effective count")
    """

    def __init__(
        self,
        n_donors: int,
        n_celltypes: int,
        *,
        momentum: float = 0.95,
        device: Optional[torch.device] = None,
        eps: float = 1e-6,
    ):
        self.n_donors = int(n_donors)
        self.n_celltypes = int(n_celltypes)
        self.momentum = float(momentum)
        self.eps = float(eps)
        dev = device if device is not None else torch.device("cpu")

        # mean and var are tracked separately via EMA over per-batch group
        # statistics. n is an exponential effective-count proxy used only
        # for downstream gating / fallbacks (e.g. ANCOVA min_donors).
        self.mean = torch.zeros(n_donors, n_celltypes, dtype=torch.float32, device=dev)
        self.var = torch.zeros(n_donors, n_celltypes, dtype=torch.float32, device=dev)
        self.n = torch.zeros(n_donors, n_celltypes, dtype=torch.float32, device=dev)
        self._initialized = torch.zeros(
            n_donors, n_celltypes, dtype=torch.bool, device=dev
        )

    # ------------------------------------------------------------------ #
    # Updates
    # ------------------------------------------------------------------ #

    @torch.no_grad()
    def update_from_batch(
        self,
        noise_scores: torch.Tensor,     # [B]
        donors: torch.Tensor,           # [B] long
        celltypes: torch.Tensor,        # [B] long
    ) -> None:
        """
        Update EMA state from a single batch, **with collective sync built in**.

        Correctness note (rev 2 design): each rank builds per-(d, t) sufficient
        stats (count, sum, sum_sq) for its local batch, then all_reduces them
        across ranks so every rank sees the same *global* batch stats. EMA is
        then applied to the global stats on every rank in lockstep, keeping
        mean/var/n bit-exactly identical across ranks. This replaces the
        previous "EMA-locally, average-occasionally" approach which produced
        biased estimates whenever a (d, t) was seen only on a subset of ranks
        (it would be scaled by 1/world_size at sync time).
        """
        if noise_scores.shape[0] == 0:
            self._cross_rank_empty_step()
            return
        device = self.mean.device
        ns = noise_scores.detach().to(device).float()
        donors = donors.to(device).long()
        celltypes = celltypes.to(device).long()

        D = self.n_donors
        T = self.n_celltypes

        # ---- 1) Local per-(d, t) sufficient stats via scatter_add ---- #
        pair_id = donors * T + celltypes                        # [B]
        flat_count = torch.zeros(D * T, dtype=torch.float32, device=device)
        flat_sum = torch.zeros(D * T, dtype=torch.float32, device=device)
        flat_sum_sq = torch.zeros(D * T, dtype=torch.float32, device=device)
        ones = torch.ones_like(ns)
        flat_count.scatter_add_(0, pair_id, ones)
        flat_sum.scatter_add_(0, pair_id, ns)
        flat_sum_sq.scatter_add_(0, pair_id, ns * ns)

        # ---- 2) All-reduce sufficient stats across ranks ---- #
        if dist.is_available() and dist.is_initialized():
            # Stack into a single tensor to halve communication cost.
            stacked = torch.stack([flat_count, flat_sum, flat_sum_sq], dim=0)
            dist.all_reduce(stacked, op=dist.ReduceOp.SUM)
            flat_count = stacked[0]
            flat_sum = stacked[1]
            flat_sum_sq = stacked[2]

        count = flat_count.reshape(D, T)
        sum_ = flat_sum.reshape(D, T)
        sum_sq = flat_sum_sq.reshape(D, T)

        valid = count > 0
        if not bool(valid.any()):
            return

        # ---- 3) Global batch mean / var from synced sufficient stats ---- #
        safe_count = count.clamp_min(1.0)
        batch_mean = sum_ / safe_count
        # Var via E[x^2] - E[x]^2; floor at 0 for numerical safety.
        batch_var = (sum_sq / safe_count - batch_mean.pow(2)).clamp_min(0.0)

        # Single-cell groups (count == 1) have batch_var = 0 by definition.
        # We replace those with a tiny floor so the EMA-tracked var does not
        # collapse to 0 from rare singletons (downstream we clamp_min anyway,
        # but keeping the floor inside the EMA target avoids a slow drift).
        single_cell_floor = self.eps * self.eps
        batch_var = torch.where(
            count >= 2.0,
            batch_var,
            torch.full_like(batch_var, single_cell_floor),
        )

        # ---- 4) EMA update in lockstep across ranks ---- #
        # mean and var are tracked with EMA (path-dependent smoothing).
        # n is tracked as a *cumulative* count of observations contributed
        # to the bank since init, which is what downstream gating logic
        # (`min_bank_n` in decomposition) wants. Using EMA on n would yield
        # a steady-state ~= batch_size, never reflecting how often a (d, t)
        # has been seen — singletons would stay singletons forever.
        mom = float(self.momentum)
        not_initialized = valid & (~self._initialized)
        needs_ema = valid & self._initialized

        if bool(not_initialized.any()):
            self.mean[not_initialized] = batch_mean[not_initialized]
            self.var[not_initialized] = batch_var[not_initialized]
            self.n[not_initialized] = count[not_initialized]
            self._initialized[not_initialized] = True

        if bool(needs_ema.any()):
            self.mean[needs_ema] = (
                mom * self.mean[needs_ema] + (1.0 - mom) * batch_mean[needs_ema]
            )
            self.var[needs_ema] = (
                mom * self.var[needs_ema] + (1.0 - mom) * batch_var[needs_ema]
            )
            # n: cumulative, not EMA.
            self.n[needs_ema] = self.n[needs_ema] + count[needs_ema]

    @torch.no_grad()
    def _cross_rank_empty_step(self) -> None:
        """Empty local batch: still participate in collective so other ranks
        don't deadlock waiting for our all_reduce."""
        if not (dist.is_available() and dist.is_initialized()):
            return
        device = self.mean.device
        empty = torch.zeros(
            3, self.n_donors * self.n_celltypes, dtype=torch.float32, device=device
        )
        dist.all_reduce(empty, op=dist.ReduceOp.SUM)
        flat_count = empty[0]
        flat_sum = empty[1]
        flat_sum_sq = empty[2]
        # Other ranks may have contributed; run the same logic as the main path
        # but skip the scatter_add (we have nothing local).
        count = flat_count.reshape(self.n_donors, self.n_celltypes)
        sum_ = flat_sum.reshape(self.n_donors, self.n_celltypes)
        sum_sq = flat_sum_sq.reshape(self.n_donors, self.n_celltypes)
        valid = count > 0
        if not bool(valid.any()):
            return
        safe_count = count.clamp_min(1.0)
        batch_mean = sum_ / safe_count
        batch_var = (sum_sq / safe_count - batch_mean.pow(2)).clamp_min(0.0)
        single_cell_floor = self.eps * self.eps
        batch_var = torch.where(
            count >= 2.0,
            batch_var,
            torch.full_like(batch_var, single_cell_floor),
        )
        mom = float(self.momentum)
        not_initialized = valid & (~self._initialized)
        needs_ema = valid & self._initialized
        if bool(not_initialized.any()):
            self.mean[not_initialized] = batch_mean[not_initialized]
            self.var[not_initialized] = batch_var[not_initialized]
            self.n[not_initialized] = count[not_initialized]
            self._initialized[not_initialized] = True
        if bool(needs_ema.any()):
            self.mean[needs_ema] = (
                mom * self.mean[needs_ema] + (1.0 - mom) * batch_mean[needs_ema]
            )
            self.var[needs_ema] = (
                mom * self.var[needs_ema] + (1.0 - mom) * batch_var[needs_ema]
            )
            # n: cumulative (see main update path).
            self.n[needs_ema] = self.n[needs_ema] + count[needs_ema]

    @torch.no_grad()
    def all_reduce_(self) -> None:
        """No-op kept for API compatibility.

        With the per-batch sufficient-stats reduce in `update_from_batch`,
        bank state is already bit-exact across ranks at all times. Calling
        this method is harmless but unnecessary.
        """
        return

    # ------------------------------------------------------------------ #
    # Accessors
    # ------------------------------------------------------------------ #

    def get_stats(self, donor: int, celltype: int) -> Tuple[float, float, float]:
        """Returns (mean, std, n) for one (donor, celltype) cell."""
        mean = float(self.mean[donor, celltype])
        std = float(self.var[donor, celltype].clamp_min(0.0).sqrt())
        n = float(self.n[donor, celltype])
        return mean, std, n

    def get_mean_tensor(self) -> torch.Tensor:
        return self.mean

    def get_std_tensor(self) -> torch.Tensor:
        return self.var.clamp_min(self.eps * self.eps).sqrt()

    def get_initialized_mask(self) -> torch.Tensor:
        return self._initialized

    # ------------------------------------------------------------------ #
    # State dict
    # ------------------------------------------------------------------ #

    def state_dict(self) -> dict:
        return {
            "mean": self.mean.detach().cpu(),
            "var": self.var.detach().cpu(),
            "n": self.n.detach().cpu(),
            "initialized": self._initialized.detach().cpu(),
        }

    def load_state_dict(self, state: dict) -> None:
        device = self.mean.device
        self.mean = state["mean"].to(device)
        self.var = state["var"].to(device)
        self.n = state["n"].to(device)
        self._initialized = state["initialized"].to(device)
