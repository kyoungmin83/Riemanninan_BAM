from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import math
import numpy as np
import torch
from torch.utils.data import DataLoader, Subset


@dataclass(frozen=True)
class ReferenceAnchoredPriorConfig:
    """
    Configuration for the v6 reference-anchored prior mean.

    `mu_p(c)` is treated as the reference-normal origin of cell type `c`.
    It is therefore refreshed from reference cells only, not learned as a
    generic trainable parameter over normal + disease cells.
    """

    enabled: bool = True
    freeze_mu: bool = True
    refresh_before_epoch: bool = True
    refresh_every_epochs: int = 1
    refresh_batch_size: int = 64
    refresh_num_workers: int = 0
    max_refresh_cells: int = 131072
    ema_momentum: float = 0.8
    min_cells_per_donor: int = 1
    min_donors_per_celltype: int = 2
    min_reference_cells: int = 4
    sync_ddp: bool = True
    eps: float = 1e-6


@dataclass(frozen=True)
class ReferencePriorRefreshSummary:
    n_active_celltypes: int
    n_total_reference_cells: int
    n_total_reference_donors: int
    mean_update_norm: float
    mean_prior_norm: float
    max_prior_norm: float


class ReferencePriorAnchor:
    """
    Maintain `prior.mu_embedding.weight` as a reference-normal origin.

    The update uses donor-balanced reference means:

        center_c = mean_d mean_i mu_q(i), i in reference cells of cell type c

    This makes `z_perp = (z_s - mu_p(c)) / sigma_p(c)` a deviation from the
    reference-normal local origin, rather than a deviation from a mixed
    normal+disease cell-type average.
    """

    def __init__(
        self,
        *,
        system: Any,
        n_celltypes: int,
        n_donors: int,
        d_z: int,
        config: ReferenceAnchoredPriorConfig,
    ) -> None:
        self.system = system
        self.n_celltypes = int(n_celltypes)
        self.n_donors = int(n_donors)
        self.d_z = int(d_z)
        self.config = config
        self.sync_process_group = None
        self.initialized = torch.zeros(self.n_celltypes, dtype=torch.bool)
        self.support_cells = torch.zeros(self.n_celltypes, dtype=torch.float32)
        self.support_donors = torch.zeros(self.n_celltypes, dtype=torch.float32)
        self.last_summary: Optional[ReferencePriorRefreshSummary] = None

    def set_sync_process_group(self, group) -> None:
        self.sync_process_group = group

    def prior(self):
        system = self.system
        if hasattr(system, "module"):
            system = system.module
        if not hasattr(system, "prior"):
            raise AttributeError("System does not expose a CellTypePrior as `.prior`.")
        return system.prior

    def freeze_prior_mu(self) -> None:
        if not bool(self.config.freeze_mu):
            return
        self.prior().mu_embedding.weight.requires_grad_(False)

    @staticmethod
    def _distributed_rank_world() -> tuple[int, int]:
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            return torch.distributed.get_rank(), torch.distributed.get_world_size()
        return 0, 1

    def _local_refresh_indices(self, n: int) -> list[int]:
        rank, world_size = self._distributed_rank_world()
        local = np.arange(rank, int(n), world_size, dtype=np.int64)

        max_cells = int(self.config.max_refresh_cells)
        if max_cells > 0:
            local_limit = max(1, int(math.ceil(max_cells / max(world_size, 1))))
            if local.size > local_limit:
                pos = np.linspace(0, local.size - 1, num=local_limit, dtype=np.int64)
                local = local[pos]

        return local.astype(np.int64).tolist()

    @torch.no_grad()
    def refresh_from_loader(
        self,
        loader,
        trainer,
        *,
        epoch_index: Optional[int] = None,
    ) -> Optional[ReferencePriorRefreshSummary]:
        if not bool(self.config.enabled):
            return None
        dataset = getattr(loader, "dataset", None)
        if dataset is None:
            return None

        local_indices = self._local_refresh_indices(len(dataset))
        refresh_ds = Subset(dataset, local_indices)
        refresh_loader = DataLoader(
            refresh_ds,
            batch_size=max(1, int(self.config.refresh_batch_size)),
            shuffle=False,
            num_workers=max(0, int(self.config.refresh_num_workers)),
            pin_memory=False,
            drop_last=False,
        )

        sums = torch.zeros(
            self.n_celltypes * self.n_donors,
            self.d_z,
            dtype=torch.float32,
            device="cpu",
        )
        counts = torch.zeros(
            self.n_celltypes * self.n_donors,
            dtype=torch.float32,
            device="cpu",
        )

        was_training = bool(getattr(trainer.system, "training", False))
        trainer.system.eval()

        for batch in refresh_loader:
            batch = trainer._move_batch_to_device(batch)
            with trainer._autocast_context():
                model_out = trainer.system(
                    batch,
                    sample_latent=False,
                    return_all_hidden_states=False,
                    return_attn_diagnostics=False,
                )

            ref = batch.get("is_reference", None)
            donor_id = batch.get("donor_id", None)
            if ref is None or donor_id is None:
                continue
            ref = ref.bool().view(-1)
            if int(ref.sum().item()) == 0:
                continue

            ct = batch["celltype_id"][ref].long().detach().cpu().clamp(
                min=0, max=self.n_celltypes - 1
            )
            dn = donor_id[ref].long().detach().cpu().clamp(
                min=0, max=self.n_donors - 1
            )
            mu_q = model_out.encoder_out.mu_q[ref].detach().float().cpu()
            flat = ct * self.n_donors + dn
            sums.index_add_(0, flat, mu_q)
            counts.index_add_(0, flat, torch.ones_like(flat, dtype=torch.float32))

        if (
            bool(self.config.sync_ddp)
            and torch.distributed.is_available()
            and torch.distributed.is_initialized()
        ):
            group = self.sync_process_group
            torch.distributed.all_reduce(sums, op=torch.distributed.ReduceOp.SUM, group=group)
            torch.distributed.all_reduce(counts, op=torch.distributed.ReduceOp.SUM, group=group)

        summary = self._apply_dense_stats(sums=sums, counts=counts)
        self.last_summary = summary

        if was_training:
            trainer.system.train()

        return summary

    @torch.no_grad()
    def _apply_dense_stats(
        self,
        *,
        sums: torch.Tensor,
        counts: torch.Tensor,
    ) -> ReferencePriorRefreshSummary:
        sums = sums.view(self.n_celltypes, self.n_donors, self.d_z)
        counts = counts.view(self.n_celltypes, self.n_donors)

        prior = self.prior()
        weight = prior.mu_embedding.weight
        device = weight.device
        dtype = weight.dtype

        new_centers = torch.zeros(self.n_celltypes, self.d_z, dtype=torch.float32)
        active = torch.zeros(self.n_celltypes, dtype=torch.bool)
        support_cells = torch.zeros(self.n_celltypes, dtype=torch.float32)
        support_donors = torch.zeros(self.n_celltypes, dtype=torch.float32)

        min_cells_per_donor = int(self.config.min_cells_per_donor)
        min_donors = int(self.config.min_donors_per_celltype)
        min_ref_cells = int(self.config.min_reference_cells)

        for c in range(self.n_celltypes):
            valid = counts[c] >= float(min_cells_per_donor)
            n_donors = int(valid.sum().item())
            n_cells = int(counts[c, valid].sum().item()) if n_donors > 0 else 0
            if n_donors < min_donors or n_cells < min_ref_cells:
                continue
            donor_means = sums[c, valid] / counts[c, valid].unsqueeze(-1).clamp_min(
                float(self.config.eps)
            )
            new_centers[c] = donor_means.mean(dim=0)
            active[c] = True
            support_cells[c] = float(n_cells)
            support_donors[c] = float(n_donors)

        old_weight = weight.detach().float().cpu()
        updated = old_weight.clone()
        momentum = float(self.config.ema_momentum)
        momentum = min(0.9999, max(0.0, momentum))

        active_idx = active.nonzero(as_tuple=True)[0]
        if int(active_idx.numel()) > 0:
            not_initialized = ~self.initialized[active_idx]
            if bool(not_initialized.any()):
                fresh_idx = active_idx[not_initialized]
                updated[fresh_idx] = new_centers[fresh_idx]
            if bool((~not_initialized).any()):
                ema_idx = active_idx[~not_initialized]
                updated[ema_idx] = (
                    momentum * old_weight[ema_idx]
                    + (1.0 - momentum) * new_centers[ema_idx]
                )

        update_norm = (updated - old_weight).norm(dim=-1)
        with torch.no_grad():
            weight.copy_(updated.to(device=device, dtype=dtype))

        self.initialized[active] = True
        self.support_cells[active] = support_cells[active]
        self.support_donors[active] = support_donors[active]

        initialized = self.initialized
        if bool(initialized.any()):
            prior_norm = updated[initialized].norm(dim=-1)
            mean_prior_norm = float(prior_norm.mean().item())
            max_prior_norm = float(prior_norm.max().item())
        else:
            mean_prior_norm = 0.0
            max_prior_norm = 0.0

        if bool(active.any()):
            mean_update_norm = float(update_norm[active].mean().item())
        else:
            mean_update_norm = 0.0

        return ReferencePriorRefreshSummary(
            n_active_celltypes=int(active.sum().item()),
            n_total_reference_cells=int(support_cells.sum().item()),
            n_total_reference_donors=int(support_donors.sum().item()),
            mean_update_norm=mean_update_norm,
            mean_prior_norm=mean_prior_norm,
            max_prior_norm=max_prior_norm,
        )
