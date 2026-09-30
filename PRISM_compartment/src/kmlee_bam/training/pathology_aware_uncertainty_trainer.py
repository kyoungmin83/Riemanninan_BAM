"""
V7a trainer — Pathology-Aware Hierarchical Uncertainty (PHU), uncertainty-only.

Wires the uncertainty/ package (DifficultyProfile + NoiseScoreBank + ANCOVAFit
+ alignment loss) into the existing V6 trainer. Adds:

  * Profile bootstrap at the start of training (or from a previous run's
    checkpoint).
  * Bank EMA update per batch.
  * Bank DDP all_reduce + ANCOVA refit at each epoch start.
  * Alignment loss = lambda * f(BAM_head_out, u_total) added to total.

The v6 reference-anchored prior, hierarchical ordinal loss, etc. are untouched.

See: model_v7a_implementation_plan.md §4.3
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import math
import os

import numpy as np
import torch
import torch.distributed as dist

try:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.reference_anchored_trainer import V6Trainer
    from kmlee_bam.training.core_trainer import ModelForwardOutput
    from kmlee_bam.uncertainty import (
        DifficultyProfile,
        DifficultyProfileConfig,
        NoiseScoreBank,
        ANCOVAConfig,
        ANCOVAFit,
        compute_noise_score,
        compute_u_total,
        stratified_cell_sample,
    )
    from kmlee_bam.objectives.uncertainty_alignment import (
        UncertaintyAlignmentConfig,
        compute_uncertainty_alignment_loss,
        compute_pathology_aware_relative_reconstruction,
    )
except ImportError:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.reference_anchored_trainer import V6Trainer
    from kmlee_bam.training.core_trainer import ModelForwardOutput
    from kmlee_bam.uncertainty import (
        DifficultyProfile,
        DifficultyProfileConfig,
        NoiseScoreBank,
        ANCOVAConfig,
        ANCOVAFit,
        compute_noise_score,
        compute_u_total,
        stratified_cell_sample,
    )
    from kmlee_bam.objectives.uncertainty_alignment import (
        UncertaintyAlignmentConfig,
        compute_uncertainty_alignment_loss,
        compute_pathology_aware_relative_reconstruction,
    )


class V7aTrainer(V6Trainer):
    """
    PHU uncertainty-only trainer on top of V6.

    Caller is responsible for providing per-donor pathology covariates and
    the donor x celltype presence mask:

        donor_pathology       : [n_donors, n_axes] float32 (e.g. Braak/Thal/CERAD)
        donor_celltype_mask   : [n_donors, n_celltypes] bool
        donor_pathology_valid : [n_donors] bool (donors with usable labels)
    """

    def __init__(
        self,
        *args,
        v7a_enabled: bool = False,
        v7a_profile_config: Optional[DifficultyProfileConfig] = None,
        v7a_bank_momentum: float = 0.95,
        v7a_ancova_config: Optional[ANCOVAConfig] = None,
        v7a_alignment_config: Optional[UncertaintyAlignmentConfig] = None,
        n_donors: Optional[int] = None,
        n_celltypes: Optional[int] = None,
        n_genes: Optional[int] = None,
        n_bins: Optional[int] = None,
        v7a_donor_pathology: Optional[torch.Tensor] = None,           # [D, P]
        v7a_donor_celltype_mask: Optional[torch.Tensor] = None,       # [D, T]
        v7a_donor_pathology_valid: Optional[torch.Tensor] = None,     # [D]
        v7a_sample_size: int = 50_000,
        v7a_profile_init_seed: int = 42,
        v7a_diag_print_every: int = 100,
        **kwargs,
    ):
        super().__init__(*args, n_celltypes=n_celltypes, n_donors=n_donors, **kwargs)

        self.v7a_enabled = bool(v7a_enabled)

        if not self.v7a_enabled:
            # Cheap no-op trainer; nothing else is constructed.
            self.v7a_profile = None
            self.v7a_bank = None
            self.v7a_ancova = None
            self.v7a_alignment_config = UncertaintyAlignmentConfig(enabled=False)
            return

        if n_donors is None or n_celltypes is None or n_genes is None or n_bins is None:
            raise ValueError(
                "v7a_enabled=True requires n_donors, n_celltypes, n_genes, n_bins."
            )

        device = self._infer_device()

        self.v7a_profile_config = v7a_profile_config or DifficultyProfileConfig()
        self.v7a_profile = DifficultyProfile(
            self.v7a_profile_config,
            n_genes=int(n_genes),
            n_bins=int(n_bins),
            n_celltypes=int(n_celltypes),
            device=device,
        )

        self.v7a_bank = NoiseScoreBank(
            n_donors=int(n_donors),
            n_celltypes=int(n_celltypes),
            momentum=float(v7a_bank_momentum),
            device=device,
        )

        self.v7a_ancova_config = v7a_ancova_config or ANCOVAConfig()
        n_axes = len(self.v7a_ancova_config.pathology_axes)
        self.v7a_ancova = ANCOVAFit(
            self.v7a_ancova_config,
            n_donors=int(n_donors),
            n_celltypes=int(n_celltypes),
            n_axes=int(n_axes),
            device=device,
        )

        self.v7a_alignment_config = v7a_alignment_config or UncertaintyAlignmentConfig(
            enabled=True
        )

        # Donor-level metadata. v7a *requires* these for ANCOVA. Fail fast
        # at construction rather than silently disabling u_donor and then
        # crashing later at epoch 2 ANCOVA refit.
        if v7a_donor_pathology is None or v7a_donor_celltype_mask is None:
            raise RuntimeError(
                "v7a enabled but donor metadata (donor_pathology / "
                "donor_celltype_mask) was not provided. Check the "
                "[v7a-meta] log in run_current.py — most likely the "
                "obs/zarr-side pathology extraction failed silently."
            )
        self.v7a_donor_pathology = self._maybe_to_device(v7a_donor_pathology, device)
        self.v7a_donor_celltype_mask = self._maybe_to_device(v7a_donor_celltype_mask, device)
        self.v7a_donor_pathology_valid = self._maybe_to_device(
            v7a_donor_pathology_valid, device
        )

        self.v7a_sample_size = int(v7a_sample_size)
        self.v7a_profile_init_seed = int(v7a_profile_init_seed)
        self.v7a_diag_print_every = max(0, int(v7a_diag_print_every))
        self._v7a_step_counter = 0

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #

    def _infer_device(self) -> torch.device:
        # Trainer.device is typically a torch.device; fall back to cuda:0 if
        # the attribute isn't present for whatever reason.
        dev = getattr(self, "device", None)
        if isinstance(dev, torch.device):
            return dev
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")

    @staticmethod
    def _maybe_to_device(
        x: Optional[torch.Tensor], device: torch.device
    ) -> Optional[torch.Tensor]:
        if x is None:
            return None
        return x.to(device)

    def _is_rank0(self) -> bool:
        return int(os.environ.get("RANK", "0")) == 0

    def attach_donor_metadata(
        self,
        donor_pathology: torch.Tensor,
        donor_celltype_mask: torch.Tensor,
        donor_pathology_valid: Optional[torch.Tensor] = None,
    ) -> None:
        """Set donor-level covariates if they were not provided at __init__."""
        if not self.v7a_enabled:
            return
        device = self._infer_device()
        self.v7a_donor_pathology = donor_pathology.to(device)
        self.v7a_donor_celltype_mask = donor_celltype_mask.to(device).bool()
        if donor_pathology_valid is None:
            n_donors = self.v7a_donor_pathology.shape[0]
            self.v7a_donor_pathology_valid = torch.ones(
                n_donors, dtype=torch.bool, device=device
            )
        else:
            self.v7a_donor_pathology_valid = donor_pathology_valid.to(device).bool()

    # ------------------------------------------------------------------ #
    # Profile bootstrap from sampled cells
    # ------------------------------------------------------------------ #

    @torch.no_grad()
    def initialize_difficulty_profile(self, loader, dataset) -> None:
        """
        Build the difficulty profile on rank 0 only, then broadcast to all
        DDP ranks.

        Why rank-0 only: the training loader is wrapped by DistributedSampler
        so each rank sees a different shard. Letting every rank build a
        profile from its own shard yields per-rank inconsistent z-score bases,
        which silently corrupts the bank / alignment supervision.

        We build a plain (non-distributed) loader on rank 0 from the same
        dataset, run forward on the first ~sample_size cells, compute the
        profile, then broadcast all profile tensors to the other ranks.
        """
        if not self.v7a_enabled:
            return

        is_dist = dist.is_available() and dist.is_initialized()
        rank = dist.get_rank() if is_dist else 0

        # WARM-START SKIP (2026-06-04): if the difficulty profile was already restored
        # from a checkpoint (run_current._warm_start_trainer_state_from_checkpoint loads
        # it IDENTICALLY on ALL ranks), skip the bootstrap entirely. This avoids re-running
        # the rank-0 forward sweep AND the post-bootstrap DDP broadcast that has hung in
        # production. Every rank sees _initialized=True -> all skip consistently -> no
        # collective is issued, so there is no rank desync.
        if self.v7a_profile is not None and getattr(self.v7a_profile, "_initialized", False):
            if rank == 0:
                print(
                    "[v7a-profile-init] profile already initialized (warm-start restore) "
                    "-> skipping bootstrap + broadcast",
                    flush=True,
                )
            return

        success = False
        error_msg = ""

        if rank == 0:
            try:
                target_cells = int(self.v7a_sample_size)
                rec_err_chunks: List[torch.Tensor] = []
                true_bins_chunks: List[torch.Tensor] = []
                celltype_chunks: List[torch.Tensor] = []
                collected = 0

                self.system.eval()
                plain_loader = self._build_plain_loader_for_profile(dataset, loader)
                for batch in plain_loader:
                    batch = self._move_batch_to_device(batch)
                    with self._autocast_context():
                        model_out = self.system(
                            batch,
                            sample_latent=False,
                            return_all_hidden_states=False,
                            return_attn_diagnostics=False,
                        )
                    nll_per_gene = model_out.decoder_out.nll_per_gene
                    if nll_per_gene is None:
                        raise RuntimeError(
                            "decoder_out.nll_per_gene is None. Profile bootstrap "
                            "requires per-gene reconstruction error from the decoder."
                        )
                    rec_err_chunks.append(nll_per_gene.detach().float().cpu())
                    # int8 storage: bin values 0..6, celltype ids 0..~30 both fit
                    # in int8 (-128..127). Saves ~4× CPU RAM vs int64 — the main
                    # CPU peak driver for production sample_size=50k.
                    true_bins_chunks.append(
                        batch["y_ord"].detach().to(torch.int8).cpu()
                    )
                    celltype_chunks.append(
                        batch["celltype_id"].detach().to(torch.int8).cpu()
                    )

                    collected += int(batch["y_ord"].shape[0])
                    if collected >= target_cells:
                        break

                rec_err = torch.cat(rec_err_chunks, dim=0)[:target_cells]
                true_bins = torch.cat(true_bins_chunks, dim=0)[:target_cells]
                celltypes = torch.cat(celltype_chunks, dim=0)[:target_cells]

                stats = self.v7a_profile.compute_from_sample(
                    rec_err, true_bins, celltypes
                )
                print(
                    f"[v7a-profile-init] sample_cells={stats['n_sample_cells']} "
                    f"gbt_valid={stats.get('gene_bin_celltype/cells_valid', 0)} "
                    f"gb_valid={stats.get('gene_bin/cells_valid', 0)} "
                    f"bt_valid={stats.get('bin_celltype/cells_valid', 0)}",
                    flush=True,
                )
                success = True
            except Exception as exc:  # noqa: BLE001 — we re-raise after broadcast
                success = False
                error_msg = repr(exc)
                print(
                    f"[v7a-profile-init | rank=0] FAILED: {error_msg}",
                    flush=True,
                )

        # Broadcast success/failure. EVERY rank participates in this collective
        # so a rank-0 failure does not leave the others hanging at barrier.
        if is_dist:
            flag = torch.tensor(
                [1.0 if success else 0.0],
                device=self._infer_device(),
                dtype=torch.float32,
            )
            dist.broadcast(flag, src=0)
            success = float(flag.item()) >= 0.5

        if not success:
            # All ranks see the same failure mode → torchrun exits cleanly.
            raise RuntimeError(
                f"v7a profile bootstrap failed on rank 0 (rank={rank}): {error_msg or 'see rank0 logs'}"
            )

        # Broadcast: rank 0 computed -> all ranks now hold identical profile.
        if is_dist:
            self.v7a_profile.broadcast_(src_rank=0)
            if rank != 0:
                # Receiver rank just got the tensors; log a confirmation.
                print(
                    f"[v7a-profile-init | rank={rank}] received broadcast",
                    flush=True,
                )

    def _build_plain_loader_for_profile(self, dataset, ddp_loader):
        """
        Construct a non-distributed DataLoader for profile bootstrap.

        We reuse the DDP loader's batch_size + num_workers settings but drop
        the DistributedSampler so rank 0 sees the full (sequential) dataset.
        """
        from torch.utils.data import DataLoader

        batch_size = getattr(ddp_loader, "batch_size", None) or 64
        num_workers = getattr(ddp_loader, "num_workers", 0) or 0
        pin_memory = getattr(ddp_loader, "pin_memory", False)

        return DataLoader(
            dataset,
            batch_size=int(batch_size),
            shuffle=False,
            num_workers=int(num_workers),
            pin_memory=bool(pin_memory),
            drop_last=False,
        )

    @staticmethod
    def _extract_celltype_ids_from_dataset(dataset, n_cells: int) -> np.ndarray:
        """Best-effort extract per-cell celltype ids without iterating full dataset."""
        for attr in ("celltype_ids", "celltype_id", "Subclass"):
            if hasattr(dataset, attr):
                ids = getattr(dataset, attr)
                ids = np.asarray(ids).astype(np.int64)
                if ids.shape[0] == n_cells:
                    return ids
        # Fallback: iterate the dataset (slow, but correct)
        ids = np.empty(n_cells, dtype=np.int64)
        for i in range(n_cells):
            item = dataset[i]
            ids[i] = int(item["celltype_id"])
        return ids

    def _collate_indices(self, dataset, idxs: Sequence[int]) -> Dict[str, torch.Tensor]:
        """Lightweight collate: stack tensors from dataset[idx] items."""
        items = [dataset[int(i)] for i in idxs]
        keys = items[0].keys()
        out: Dict[str, torch.Tensor] = {}
        for k in keys:
            vals = [item[k] for item in items]
            if torch.is_tensor(vals[0]):
                out[k] = torch.stack(vals, dim=0)
            else:
                out[k] = torch.as_tensor(vals)
        return out

    # ------------------------------------------------------------------ #
    # Epoch-level hooks
    # ------------------------------------------------------------------ #

    def train_epoch(self, loader, *args, epoch_index: Optional[int] = None, **kwargs):
        # Profile init at the start of epoch 1 (or first call) if v7a is
        # enabled but profile hasn't been built yet. We use loader.dataset
        # to obtain stratified samples.
        if (
            self.v7a_enabled
            and self.v7a_profile is not None
            and not self.v7a_profile.initialized
        ):
            dataset = getattr(loader, "dataset", None)
            if dataset is None:
                # batch_sampler-wrapped loaders also expose .dataset via the
                # sampler; try that too.
                sampler = getattr(loader, "batch_sampler", None)
                dataset = getattr(sampler, "dataset", None) if sampler is not None else None
            if dataset is None:
                if self._is_rank0():
                    print(
                        "[v7a-profile-init] WARNING: could not resolve dataset "
                        "from loader; skipping profile init this epoch.",
                        flush=True,
                    )
            else:
                if self._is_rank0():
                    print(
                        f"[v7a-profile-init] starting profile bootstrap at "
                        f"epoch={epoch_index} sample_size={self.v7a_sample_size}",
                        flush=True,
                    )
                self.initialize_difficulty_profile(loader, dataset)

        # ANCOVA refit at the start of each epoch (after the first; on
        # epoch 1 the bank is still being warmed up).
        # Bank itself stays in sync across DDP ranks via per-batch
        # sufficient-stats reduce inside update_from_batch (see bank.py).
        if (
            self.v7a_enabled
            and self.v7a_profile is not None
            and self.v7a_profile.initialized
            and epoch_index is not None
            and int(epoch_index) > 1
        ):
            fit_stats = self.v7a_ancova.fit(
                self.v7a_bank,
                self.v7a_donor_pathology,
                self.v7a_donor_celltype_mask,
                donor_pathology_valid=self.v7a_donor_pathology_valid,
            )
            if self._is_rank0():
                print(
                    f"[v7a-ancova] epoch={epoch_index} "
                    f"fitted={fit_stats['n_celltypes_fitted']}/"
                    f"{fit_stats['n_celltypes_total']} "
                    f"beta_abs_mean={fit_stats['beta_abs_mean']:.4f} "
                    f"residual_std_mean={fit_stats['residual_std_mean']:.4f}",
                    flush=True,
                )

        return super().train_epoch(loader, *args, epoch_index=epoch_index, **kwargs)

    # ------------------------------------------------------------------ #
    # Loss composition
    # ------------------------------------------------------------------ #

    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        out = super()._compute_loss(batch, model_out)

        if not self.v7a_enabled:
            return out

        if self.v7a_profile is None or not self.v7a_profile.initialized:
            # Profile not yet initialized; alignment is zero this batch.
            return out

        # Required tensors.
        nll_per_gene = model_out.decoder_out.nll_per_gene
        if nll_per_gene is None:
            return out  # Decoder did not compute per-gene NLL; skip silently.

        if "donor_id" not in batch:
            return out
        donors = batch["donor_id"].long()
        celltypes = batch["celltype_id"].long()
        true_bins = batch["y_ord"].long()

        # Compute noise_score (detached internally — no grad through bank).
        noise_scores = compute_noise_score(
            rec_err=nll_per_gene.detach().float(),
            true_bins=true_bins,
            celltypes=celltypes,
            profile=self.v7a_profile,
            winsorize_z=float(self.v7a_profile.config.winsorize_z),
            trim_fraction=0.10,
        )

        # Update bank ONLY during training. Two reasons:
        #   1) Val/test cells should not contaminate the train-side noise
        #      statistics that drive alignment supervision.
        #   2) `bank.update_from_batch` performs a collective all_reduce. If
        #      validation runs in rank0-only mode, only rank 0 would enter the
        #      collective and the other ranks would deadlock.
        # `system.training` flips off during evaluate_epoch; `is_grad_enabled`
        # is False under @torch.no_grad context. Both checks together are
        # robust.
        is_training_step = (
            self.system.training and torch.is_grad_enabled()
        )
        if is_training_step:
            self.v7a_bank.update_from_batch(noise_scores, donors, celltypes)

        # Compute alignment loss.
        bam_head = model_out.encoder_out.cell_uncertainty
        if bam_head is None or not self.v7a_alignment_config.enabled:
            return out

        # u_total: bank/ANCOVA-derived target (no grad).
        u_total = compute_u_total(
            noise_scores=noise_scores,
            donors=donors,
            celltypes=celltypes,
            bank=self.v7a_bank,
            ancova=self.v7a_ancova,
            u_total_cap=self.v7a_alignment_config.u_total_cap,
        )

        # Decide epoch for ramp (best-effort via attribute).
        epoch_for_ramp = getattr(self, "current_epoch_index", None)
        align_loss = compute_uncertainty_alignment_loss(
            bam_head_output=bam_head,
            u_total=u_total,
            config=self.v7a_alignment_config,
            epoch=epoch_for_ramp,
        )

        # Surface the *actual* delayed curriculum scale separately from the
        # configured lambda.  The old console only printed the alignment loss,
        # which made an intentionally inactive pre-start epoch look like a
        # failed uncertainty model.  Keep this calculation local so logging is
        # explicit without coupling the trainer to a private objective helper.
        if epoch_for_ramp is None:
            alignment_ramp = 1.0
        elif int(epoch_for_ramp) < int(self.v7a_alignment_config.alignment_start_epoch):
            alignment_ramp = 0.0
        elif int(self.v7a_alignment_config.ramp_epochs) <= 0:
            alignment_ramp = 1.0
        else:
            alignment_ramp = float(
                min(
                    1.0,
                    (
                        int(epoch_for_ramp)
                        - int(self.v7a_alignment_config.alignment_start_epoch)
                        + 1
                    )
                    / float(self.v7a_alignment_config.ramp_epochs),
                )
            )

        # Add weighted alignment loss to total during training only.
        # Validation/test often contains held-out donors that never populate the
        # train-side bank; their u_total target can be all-zero by construction.
        # Logging the diagnostic is useful, but adding that invalid target to
        # val_total corrupts checkpoint selection.
        lambda_align_config = float(self.v7a_alignment_config.lambda_alignment)
        lambda_align = lambda_align_config if is_training_step else 0.0
        align_weighted = lambda_align * align_loss

        # Optional centred PHU reweighting.  The default source is detached
        # u_total, which is gene/bin/celltype-normalised and donor-pathology
        # adjusted.  Validation remains diagnostic-only because its donor bank
        # entries are intentionally not populated.
        relative_out = compute_pathology_aware_relative_reconstruction(
            bam_head_output=bam_head,
            u_total=u_total,
            rec_per_cell=model_out.decoder_out.nll_per_cell,
            config=self.v7a_alignment_config,
            epoch=epoch_for_ramp,
            group_id=celltypes,
        )
        lambda_relative_config = float(
            self.v7a_alignment_config.lambda_relative_rec
        )
        lambda_relative = lambda_relative_config if is_training_step else 0.0
        relative_weighted = (
            lambda_relative * float(relative_out.ramp) * relative_out.delta
        )

        new_total = out.total + align_weighted + relative_weighted

        # Update LossOutput fields.
        details = dict(out.details)
        details["loss/uncertainty_alignment"] = float(align_loss.detach())
        details["weight/lambda_alignment"] = lambda_align
        details["weight/phu_alignment_ramp"] = float(alignment_ramp)
        details["weight/lambda_alignment_effective"] = float(
            lambda_align * alignment_ramp
        )
        details["diag/u_total_mean"] = float(u_total.mean().detach())
        details["diag/u_total_std"] = float(u_total.std().detach())
        details["diag/noise_score_mean"] = float(noise_scores.mean().detach())
        details["diag/noise_score_std"] = float(noise_scores.std().detach())
        details["loss/phu_relative_rec_delta"] = float(relative_out.delta.detach())
        details["loss/phu_relative_rec_weighted"] = float(
            relative_out.weighted_rec.detach()
        )
        details["weight/lambda_phu_relative_rec"] = float(lambda_relative)
        details["weight/phu_relative_ramp"] = float(relative_out.ramp)
        details["weight/phu_rel_w_min"] = float(relative_out.weight_min)
        details["weight/phu_rel_w_mean"] = float(relative_out.weight_mean)
        details["weight/phu_rel_w_max"] = float(relative_out.weight_max)
        details["metric/phu_rel_ess_fraction"] = float(
            relative_out.effective_sample_fraction
        )
        details["metric/phu_bam_target_corr"] = float(relative_out.bam_target_corr)
        details["loss/total"] = float(new_total.detach())

        # (Separate [v7a-step] print removed — diagnostics flow through the
        # unified [train step] log via details dict.)

        # Reconstruct the LossOutput with the updated total.
        # Note: the v7a stats above (loss/uncertainty_alignment, diag/u_total_*,
        # diag/noise_score_*, weight/lambda_alignment) flow through into the
        # main [train step] log via the trainer's _format_step_log integration
        # so we do NOT emit a separate [v7a-step] line. The previous separate
        # log fired at a different cadence and produced duplicate / noisy
        # output; consolidation makes the per-step log a single coherent block.
        out.total = new_total
        out.details = details
        return out

    # ------------------------------------------------------------------ #
    # State dict (for checkpoint persistence)
    # ------------------------------------------------------------------ #

    def v7a_state_dict(self) -> dict:
        if not self.v7a_enabled:
            return {}
        return {
            "profile": self.v7a_profile.state_dict() if self.v7a_profile else None,
            "bank": self.v7a_bank.state_dict() if self.v7a_bank else None,
            "ancova": self.v7a_ancova.state_dict() if self.v7a_ancova else None,
        }

    def v7a_load_state_dict(self, state: dict) -> None:
        if not self.v7a_enabled:
            return
        if "profile" in state and state["profile"] is not None:
            self.v7a_profile.load_state_dict(state["profile"])
        if "bank" in state and state["bank"] is not None:
            self.v7a_bank.load_state_dict(state["bank"])
        if "ancova" in state and state["ancova"] is not None:
            self.v7a_ancova.load_state_dict(state["ancova"])
