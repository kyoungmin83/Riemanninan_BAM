from __future__ import annotations

from typing import Dict, Optional

import torch

try:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.objectives.reference_centering import (
        DonorBalancedReferenceCenterConfig,
        donor_balanced_reference_center_stats,
    )
    from kmlee_bam.objectives.aux_reference_adversary import (
        AuxReferenceAdversaryConfig,
        AuxReferenceAdversaryRunner,
        _get_celltype_adversary,
        _unwrap_system,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput, Trainer
except ImportError:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.objectives.reference_centering import (
        DonorBalancedReferenceCenterConfig,
        donor_balanced_reference_center_stats,
    )
    from kmlee_bam.objectives.aux_reference_adversary import (
        AuxReferenceAdversaryConfig,
        AuxReferenceAdversaryRunner,
        _get_celltype_adversary,
        _unwrap_system,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput, Trainer


class V3Trainer(Trainer):
    """
    Trainer extension for scoped v3 losses.

    v3 keeps the v2 reconstruction / state KL / BAM objective intact, but adds:

    1. donor-balanced reference centering
       - replaces ordinary cell-count reference centering
       - requires batch["donor_id"]

    2. weak warm-started conditional tech adversary
       - expects model_out.tech_adversary_out
       - adversary itself should apply GRL only to z_dev

    3. explicit v3 diagnostics
       - printed directly as [kmlee-loss:adversary]
       - useful because the base compact logger may not print details dict
    """

    def __init__(
        self,
        *args,
        lambda_ref_center_db: float = 0.0,
        ref_center_min_cells_per_donor: int = 2,
        ref_center_min_donors_per_celltype: int = 2,
        lambda_tech_adv: float = 0.0,
        tech_adv_warmup_epochs: int = 5,
        tech_adv_ramp_epochs: int = 0,
        tech_adv_max_lambda: Optional[float] = None,
        lambda_celltype_adv: float = 0.0,
        celltype_adv_warmup_epochs: int = 5,
        celltype_adv_ramp_epochs: int = 0,
        celltype_adv_max_lambda: Optional[float] = None,
        lambda_sex_adv: float = 0.0,
        sex_adv_warmup_epochs: int = 5,
        sex_adv_ramp_epochs: int = 0,
        sex_adv_max_lambda: Optional[float] = None,
        aux_celltype_adv_config: Optional[AuxReferenceAdversaryConfig] = None,
        v3_diag_print_every: int = 100,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)

        # --------------------------------------------------------------
        # Donor-balanced reference centering
        # --------------------------------------------------------------
        self.lambda_ref_center_db = float(lambda_ref_center_db)
        self.ref_center_db_config = DonorBalancedReferenceCenterConfig(
            min_cells_per_donor=int(ref_center_min_cells_per_donor),
            min_donors_per_celltype=int(ref_center_min_donors_per_celltype),
        )

        # --------------------------------------------------------------
        # Conditional tech adversary schedule
        # --------------------------------------------------------------
        self.lambda_tech_adv = float(lambda_tech_adv)
        self.tech_adv_warmup_epochs = int(tech_adv_warmup_epochs)
        self.tech_adv_ramp_epochs = max(0, int(tech_adv_ramp_epochs))
        self.tech_adv_max_lambda = (
            float(tech_adv_max_lambda)
            if tech_adv_max_lambda is not None
            else float(lambda_tech_adv)
        )

        # --------------------------------------------------------------
        # Reference-only celltype adversary schedule (mirrors tech adv)
        # --------------------------------------------------------------
        self.lambda_celltype_adv = float(lambda_celltype_adv)
        self.celltype_adv_warmup_epochs = int(celltype_adv_warmup_epochs)
        self.celltype_adv_ramp_epochs = max(0, int(celltype_adv_ramp_epochs))
        self.celltype_adv_max_lambda = (
            float(celltype_adv_max_lambda)
            if celltype_adv_max_lambda is not None
            else float(lambda_celltype_adv)
        )

        # --------------------------------------------------------------
        # v31 sex adversary schedule (soft sex erasure on z_clean; mirrors above)
        # --------------------------------------------------------------
        self.lambda_sex_adv = float(lambda_sex_adv)
        self.sex_adv_warmup_epochs = int(sex_adv_warmup_epochs)
        self.sex_adv_ramp_epochs = max(0, int(sex_adv_ramp_epochs))
        self.sex_adv_max_lambda = (
            float(sex_adv_max_lambda)
            if sex_adv_max_lambda is not None
            else float(lambda_sex_adv)
        )

        # --------------------------------------------------------------
        # v15 auxiliary reference batch for the celltype adversary.
        #
        # OFF by default: the runner is built only when the config is enabled
        # AND the model actually has a celltype adversary attached. Otherwise
        # self._aux_celltype_adv_runner stays None and train_epoch falls back to
        # the v14 path (byte-identical). The runner is lazily constructed in
        # train_epoch (first call) because it needs the train loader's dataset.
        # --------------------------------------------------------------
        self.aux_celltype_adv_config = aux_celltype_adv_config or AuxReferenceAdversaryConfig()
        self._aux_celltype_adv_runner: Optional[AuxReferenceAdversaryRunner] = None
        self._aux_celltype_adv_built = False
        # Last aux-step stats, surfaced in the details dict / diagnostic print so
        # the (now meaningful, n≈aux_batch_size) balanced-acc is visible.
        self._aux_celltype_adv_last: Optional[Dict[str, float]] = None

        # Current epoch is set by train_epoch().
        # Assumption: epoch_index is 1-based, as in logs epoch 001, 002, ...
        self.current_epoch_index: Optional[int] = None

        # Explicit v3 diagnostic print counter.
        self.v3_diag_print_every = max(0, int(v3_diag_print_every))
        self._v3_loss_print_counter = 0

    @staticmethod
    def _tech_collapse_summary(
        true_counts: Optional[torch.Tensor],
        pred_counts: Optional[torch.Tensor],
        recall_per_class: Optional[torch.Tensor],
    ) -> tuple[str, Dict[str, float]]:
        """Summarize majority-class collapse risk for the tech adversary."""
        if true_counts is None or pred_counts is None or recall_per_class is None:
            return "NA", {}
        with torch.no_grad():
            true_f = true_counts.detach().float()
            pred_f = pred_counts.detach().float()
            recall_f = recall_per_class.detach().float()
            present = true_f > 0
            n_present = int(present.sum().item())
            pred_total = pred_f.sum().clamp_min(1.0)
            pred_major_frac_t = pred_f.max() / pred_total
            if present.any():
                min_recall_t = recall_f[present].min()
            else:
                min_recall_t = torch.zeros((), dtype=recall_f.dtype, device=recall_f.device)

            pred_major_frac = float(pred_major_frac_t.detach().cpu())
            min_recall = float(min_recall_t.detach().cpu())

        if n_present < 2:
            status = "single-class batch"
        elif pred_major_frac >= 0.95 and min_recall <= 0.10:
            status = "SEVERE majority-class collapse"
        elif pred_major_frac >= 0.85 or min_recall <= 0.20:
            status = "WARN skew/collapse risk"
        else:
            status = "ok"

        msg = (
            f"{status}  pred_major_frac={pred_major_frac:.3f}  "
            f"min_recall={min_recall:.3f}"
        )
        return msg, {
            "tech_adv_pred_major_frac": pred_major_frac,
            "tech_adv_min_recall_present": min_recall,
            "tech_adv_present_classes": float(n_present),
        }

    # ------------------------------------------------------------------
    # Epoch bookkeeping
    # ------------------------------------------------------------------
    def _maybe_build_aux_celltype_adv_runner(self, loader) -> None:
        """Lazily build the v15 aux-reference-adversary runner from the train
        loader's dataset. No-op unless the aux config is enabled AND the model
        actually carries a celltype adversary (so disabled / off-by-default runs
        never touch the dataset and stay byte-identical)."""
        if self._aux_celltype_adv_built:
            return
        self._aux_celltype_adv_built = True

        cfg = self.aux_celltype_adv_config
        if not bool(getattr(cfg, "enabled", False)) or int(getattr(cfg, "aux_batch_size", 0)) <= 0:
            return
        if _get_celltype_adversary(self.system) is None:
            if self._is_rank0():
                print(
                    "[kmlee-aux-adv] aux reference batch is enabled but the model has no "
                    "celltype adversary; aux path disabled.",
                    flush=True,
                )
            return
        dataset = getattr(loader, "dataset", None)
        if dataset is None:
            return
        try:
            runner = AuxReferenceAdversaryRunner(dataset=dataset, config=cfg)
        except Exception as exc:  # noqa: BLE001
            if self._is_rank0():
                print(
                    f"[kmlee-aux-adv] WARNING: failed to build aux reference bank: {exc}; "
                    "aux path disabled.",
                    flush=True,
                )
            return
        if not runner.is_active():
            if self._is_rank0():
                print(
                    "[kmlee-aux-adv] aux reference bank has no usable reference cells; "
                    "aux path disabled.",
                    flush=True,
                )
            return
        self._aux_celltype_adv_runner = runner
        if self._is_rank0():
            print(
                "[kmlee-aux-adv] auxiliary reference adversary ACTIVE "
                f"aux_batch_size={cfg.aux_batch_size} aux_every={cfg.aux_every} "
                f"classifier_pretrain_epochs={cfg.classifier_pretrain_epochs} "
                f"celltype_balanced={cfg.celltype_balanced} "
                f"n_reference_cells={runner.bank.n_reference} "
                f"n_celltypes_with_ref={int(runner.bank.celltypes.size)} "
                f"grl_pretrain_strength={cfg.grl_pretrain_strength}",
                flush=True,
            )

    def _run_aux_celltype_adv_step(
        self, step_idx: int, *, at_optimizer_boundary: bool = True
    ) -> None:
        """One self-contained adversary update on a dedicated reference batch.

        Runs OUTSIDE the main batch's autograd graph: zero_grad → fresh
        encoder-only forward on aux_batch_size reference cells (on the UNWRAPPED
        module, ``compute_decoder=False`` — see ``AuxReferenceAdversaryRunner``)
        → adversary CE/GRL loss → backward → manual all-reduce of the aux grads
        across DDP ranks → optimizer step. Only the encoder (via GRL) and the
        adversary head receive gradients; the decoder / prior / classifier are
        untouched because their params have no grad after zero_grad.

        DDP sync (the v15 crash fix, part 2): because the aux forward runs on the
        UNWRAPPED module (so it never engages DDP's reducer — part 1), each rank's
        aux gradient is LOCAL (ranks sample different reference cells). Without
        syncing, the 4 replicas would diverge and DDP would break. So AFTER a
        successful backward and BEFORE the optimizer step we average the aux grads
        across ranks with an explicit all_reduce (which also realizes the intended
        world_size×aux_batch_size effective batch and works fine under
        NCCL_P2P_DISABLE=1). Placed before grad-clip/step so clipping sees the
        synced grad.

        Grad-accum boundary guard (refinement B): this step calls
        ``optimizer.zero_grad()``, which would WIPE the partially-accumulated MAIN
        gradients when ``grad_accum_steps>1`` and the main optimizer has not
        stepped yet this micro-step. So skip unless ``at_optimizer_boundary`` —
        i.e. the callback fired right after a real main optimizer step (detected
        by ``step_out.grad_norm is not None``). With ``grad_accum_steps=1`` every
        micro-step is a boundary, so behavior is unchanged.
        """
        runner = self._aux_celltype_adv_runner
        if runner is None or not runner.should_run_this_step(step_idx):
            return
        if not at_optimizer_boundary:
            return

        was_training = bool(getattr(self.system, "training", False))
        self.system.train()
        self.optimizer.zero_grad(set_to_none=True)

        result = runner.aux_loss_and_stats(
            self,
            epoch_index=self.current_epoch_index,
            step_idx=int(step_idx),
            sample_latent=True,
        )
        if result is None:
            if not was_training:
                self.system.eval()
            return

        aux_loss, stats = result
        backward_ok = self._backward(aux_loss)
        if backward_ok:
            # DDP sync (part 2): the aux forward ran on the UNWRAPPED module, so
            # each rank's grad is local (ranks draw different reference cells).
            # Average the aux grads across ranks BEFORE grad-clip/step so the 4
            # replicas stay identical and clipping sees the synced grad. After
            # zero_grad + the encoder-only aux backward, only the encoder and the
            # adversary head carry grad (p.grad is not None), so this touches
            # exactly those params. all_reduce works under NCCL_P2P_DISABLE=1.
            if (
                torch.distributed.is_available()
                and torch.distributed.is_initialized()
                and torch.distributed.get_world_size() > 1
            ):
                ws = torch.distributed.get_world_size()
                unwrapped = _unwrap_system(self.system)
                for p in unwrapped.parameters():
                    if p.grad is not None:
                        torch.distributed.all_reduce(
                            p.grad, op=torch.distributed.ReduceOp.SUM
                        )
                        p.grad /= ws
            self._finish_optimizer_step()
        else:
            self.optimizer.zero_grad(set_to_none=True)

        self._aux_celltype_adv_last = {
            "loss": float(stats.loss),
            "balanced_accuracy": float(stats.balanced_accuracy),
            "accuracy": float(stats.accuracy),
            "n_reference": float(stats.n_reference),
            "pretrain": 1.0 if stats.pretrain else 0.0,
            "grl_strength": float(stats.grl_strength),
        }
        if not was_training:
            self.system.eval()

    def train_epoch(self, *args, epoch_index: Optional[int] = None, **kwargs):
        self.current_epoch_index = epoch_index

        loader = args[0] if args else kwargs.get("loader")
        if loader is not None:
            self._maybe_build_aux_celltype_adv_runner(loader)

        # Wrap the runner's per-step aux adversary update around the existing
        # step_callback. The aux update runs after the main optimizer step
        # (the callback fires post-step), keeping it outside the main batch's
        # gradient-accumulation / no_sync window. No runner ⇒ callback passes
        # through untouched (byte-identical).
        if self._aux_celltype_adv_runner is not None:
            user_callback = kwargs.get("step_callback")

            def _aux_wrapped_callback(*, step_idx, n_total_steps, step_out, meters, epoch_index=None):
                # core_trainer sets step_out.grad_norm to a float ONLY on a true
                # optimizer-step boundary (None on accumulating micro-steps); use
                # it to gate the aux zero_grad/step so grad_accum_steps>1 cannot
                # wipe accumulating main gradients (refinement B).
                at_boundary = getattr(step_out, "grad_norm", None) is not None
                self._run_aux_celltype_adv_step(
                    int(step_idx), at_optimizer_boundary=at_boundary
                )
                if user_callback is not None:
                    user_callback(
                        step_idx=step_idx,
                        n_total_steps=n_total_steps,
                        step_out=step_out,
                        meters=meters,
                        epoch_index=epoch_index,
                    )

            kwargs["step_callback"] = _aux_wrapped_callback

        try:
            return super().train_epoch(*args, epoch_index=epoch_index, **kwargs)
        finally:
            self.current_epoch_index = None
            self._aux_celltype_adv_last = None

    # ------------------------------------------------------------------
    # Tech adversary weight schedule
    # ------------------------------------------------------------------
    def _tech_adv_weight(self) -> float:
        """
        Schedule, assuming 1-based epoch_index:

            epochs 1..warmup_epochs:
                0

            first epoch after warmup:
                lambda_tech_adv

            then linearly ramp to tech_adv_max_lambda over ramp_epochs.

        Example:
            warmup_epochs = 5
            lambda_tech_adv = 0.001
            tech_adv_max_lambda = 0.005
            tech_adv_ramp_epochs = 10

            epoch 1-5  : 0
            epoch 6    : 0.001
            epoch 15   : 0.005
            epoch 16+  : 0.005
        """
        if self.lambda_tech_adv <= 0.0:
            return 0.0

        if self.current_epoch_index is None:
            return 0.0

        epoch = int(self.current_epoch_index)

        if epoch <= self.tech_adv_warmup_epochs:
            return 0.0

        if self.tech_adv_ramp_epochs <= 0:
            return min(self.lambda_tech_adv, self.tech_adv_max_lambda)

        if self.tech_adv_max_lambda <= self.lambda_tech_adv:
            return min(self.lambda_tech_adv, self.tech_adv_max_lambda)

        epochs_after_warmup = epoch - self.tech_adv_warmup_epochs

        if self.tech_adv_ramp_epochs == 1:
            alpha = 1.0
        else:
            alpha = float(epochs_after_warmup - 1) / float(self.tech_adv_ramp_epochs - 1)
            alpha = min(1.0, max(0.0, alpha))

        return self.lambda_tech_adv + alpha * (
            self.tech_adv_max_lambda - self.lambda_tech_adv
        )

    @staticmethod
    def _adv_weight_schedule(
        epoch: Optional[int],
        lambda_base: float,
        warmup_epochs: int,
        ramp_epochs: int,
        max_lambda: float,
    ) -> float:
        """Generic warmup+linear-ramp adversary weight (1-based epoch_index).

        Same shape as :meth:`_tech_adv_weight`: 0 during warmup, ``lambda_base``
        on the first post-warmup epoch, then a linear ramp to ``max_lambda`` over
        ``ramp_epochs``.
        """
        if lambda_base <= 0.0 or epoch is None:
            return 0.0
        epoch = int(epoch)
        if epoch <= warmup_epochs:
            return 0.0
        if ramp_epochs <= 0 or max_lambda <= lambda_base:
            return min(lambda_base, max_lambda)
        epochs_after_warmup = epoch - warmup_epochs
        if ramp_epochs == 1:
            alpha = 1.0
        else:
            alpha = float(epochs_after_warmup - 1) / float(ramp_epochs - 1)
            alpha = min(1.0, max(0.0, alpha))
        return lambda_base + alpha * (max_lambda - lambda_base)

    def _celltype_adv_weight(self) -> float:
        return self._adv_weight_schedule(
            self.current_epoch_index,
            self.lambda_celltype_adv,
            self.celltype_adv_warmup_epochs,
            self.celltype_adv_ramp_epochs,
            self.celltype_adv_max_lambda,
        )

    def _sex_adv_weight(self) -> float:
        return self._adv_weight_schedule(
            self.current_epoch_index,
            self.lambda_sex_adv,
            self.sex_adv_warmup_epochs,
            self.sex_adv_ramp_epochs,
            self.sex_adv_max_lambda,
        )

    # ------------------------------------------------------------------
    # Rank helper
    # ------------------------------------------------------------------
    @staticmethod
    def _is_rank0() -> bool:
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            return torch.distributed.get_rank() == 0
        return True

    # ------------------------------------------------------------------
    # Main v3 loss extension
    # ------------------------------------------------------------------
    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        # First compute the ordinary v2 loss.
        out = super()._compute_loss(batch, model_out)

        total = out.total
        details = dict(out.details)

        # --------------------------------------------------------------
        # 1. Donor-balanced reference centering
        # --------------------------------------------------------------
        ref_center_db = torch.zeros((), dtype=total.dtype, device=total.device)
        ref_center_db_active_fraction = torch.zeros(
            (), dtype=total.dtype, device=total.device
        )
        ref_center_db_n_ref_celltypes = 0
        ref_center_db_n_active_celltypes = 0

        if self.lambda_ref_center_db > 0.0 and "donor_id" in batch:
            ref_stats = donor_balanced_reference_center_stats(
                z_perp=model_out.z_perp,
                celltype_id=batch["celltype_id"],
                donor_id=batch["donor_id"],
                is_reference=batch.get("is_reference", None),
                config=self.ref_center_db_config,
            )

            ref_center_db = ref_stats.loss
            ref_center_db_active_fraction = ref_stats.active_fraction.to(
                dtype=total.dtype
            )
            ref_center_db_n_ref_celltypes = int(ref_stats.n_reference_celltypes)
            ref_center_db_n_active_celltypes = int(ref_stats.n_active_celltypes)

            total = total + self.lambda_ref_center_db * ref_center_db

        # --------------------------------------------------------------
        # 2. Conditional tech adversary
        # --------------------------------------------------------------
        tech_adv_loss = torch.zeros((), dtype=total.dtype, device=total.device)
        tech_adv_acc = torch.zeros((), dtype=total.dtype, device=total.device)
        tech_adv_bal_acc = torch.zeros((), dtype=total.dtype, device=total.device)
        tech_adv_true_counts = None
        tech_adv_pred_counts = None
        tech_adv_correct_counts = None
        tech_adv_recall_per_class = None

        tech_adv_weight = self._tech_adv_weight()
        tech_adv_out = getattr(model_out, "tech_adversary_out", None)

        if tech_adv_out is not None:
            tech_adv_loss = tech_adv_out.loss
            tech_adv_acc = tech_adv_out.accuracy.to(dtype=total.dtype)
            tech_adv_bal_acc = tech_adv_out.balanced_accuracy.to(dtype=total.dtype)
            tech_adv_true_counts = tech_adv_out.true_counts.detach()
            tech_adv_pred_counts = tech_adv_out.pred_counts.detach()
            tech_adv_correct_counts = tech_adv_out.correct_counts.detach()
            tech_adv_recall_per_class = tech_adv_out.recall_per_class.detach()

            total = total + tech_adv_weight * tech_adv_loss

        # --------------------------------------------------------------
        # 2b. Reference-only celltype adversary
        #     Pushes z_perp to NOT encode celltype identity on control cells.
        #     The adversary already masks the CE / balanced-acc to reference
        #     cells (and returns zero loss when a batch has none), so here we
        #     just apply the warmup-scheduled weight.
        # --------------------------------------------------------------
        celltype_adv_loss = torch.zeros((), dtype=total.dtype, device=total.device)
        celltype_adv_acc = torch.zeros((), dtype=total.dtype, device=total.device)
        celltype_adv_bal_acc = torch.zeros((), dtype=total.dtype, device=total.device)
        celltype_adv_n_ref = torch.zeros((), dtype=total.dtype, device=total.device)

        celltype_adv_weight = self._celltype_adv_weight()
        celltype_adv_out = getattr(model_out, "celltype_adversary_out", None)

        if celltype_adv_out is not None:
            celltype_adv_loss = celltype_adv_out.loss
            celltype_adv_acc = celltype_adv_out.accuracy.to(dtype=total.dtype)
            celltype_adv_bal_acc = celltype_adv_out.balanced_accuracy.to(dtype=total.dtype)
            celltype_adv_n_ref = celltype_adv_out.n_reference.to(dtype=total.dtype)

            total = total + celltype_adv_weight * celltype_adv_loss

        # --------------------------------------------------------------
        # 2c. v31 sex adversary — SOFT sex erasure on z_clean (post-projection).
        #     GRL pushes the encoder to hide sex anywhere in z, catching the refill
        #     the hard projection can't. Warmup-scheduled, small lambda.
        # --------------------------------------------------------------
        sex_adv_loss = torch.zeros((), dtype=total.dtype, device=total.device)
        sex_adv_bal_acc = torch.zeros((), dtype=total.dtype, device=total.device)
        sex_adv_weight = self._sex_adv_weight()
        sex_adv_out = getattr(model_out, "sex_adversary_out", None)
        if sex_adv_out is not None:
            sex_adv_loss = sex_adv_out.loss
            sex_adv_bal_acc = sex_adv_out.balanced_accuracy.to(dtype=total.dtype)
            total = total + sex_adv_weight * sex_adv_loss

        # --------------------------------------------------------------
        # 3. Details dict for checkpoint / aggregated logging
        # --------------------------------------------------------------
        details["loss/ref_center_donor_balanced"] = float(
            ref_center_db.detach().cpu()
        )
        details["metric/ref_center_db_active_fraction"] = float(
            ref_center_db_active_fraction.detach().cpu()
        )
        details["metric/ref_center_db_n_reference_celltypes"] = float(
            ref_center_db_n_ref_celltypes
        )
        details["metric/ref_center_db_n_active_celltypes"] = float(
            ref_center_db_n_active_celltypes
        )

        details["loss/tech_adv"] = float(tech_adv_loss.detach().cpu())
        details["metric/tech_adv_acc"] = float(tech_adv_acc.detach().cpu())
        details["metric/tech_adv_bal_acc"] = float(tech_adv_bal_acc.detach().cpu())
        if tech_adv_true_counts is not None:
            for i in range(int(tech_adv_true_counts.numel())):
                details[f"metric/tech_adv_true_n_{i}"] = float(
                    tech_adv_true_counts[i].detach().cpu()
                )
                details[f"metric/tech_adv_pred_n_{i}"] = float(
                    tech_adv_pred_counts[i].detach().cpu()
                )
                details[f"metric/tech_adv_recall_{i}"] = float(
                    tech_adv_recall_per_class[i].detach().cpu()
                )
            _, collapse_metrics = self._tech_collapse_summary(
                tech_adv_true_counts,
                tech_adv_pred_counts,
                tech_adv_recall_per_class,
            )
            for key, value in collapse_metrics.items():
                details[f"metric/{key}"] = value

        details["loss/sex_adv"] = float(sex_adv_loss.detach().cpu())
        details["metric/sex_adv_bal_acc"] = float(sex_adv_bal_acc.detach().cpu())
        details["weight/lambda_sex_adv_active"] = float(sex_adv_weight)
        details["weight/lambda_sex_adv_target"] = float(self.lambda_sex_adv)

        details["loss/celltype_adv"] = float(celltype_adv_loss.detach().cpu())
        details["metric/celltype_adv_acc"] = float(celltype_adv_acc.detach().cpu())
        details["metric/celltype_adv_balanced_acc"] = float(
            celltype_adv_bal_acc.detach().cpu()
        )
        details["metric/celltype_adv_n_reference"] = float(
            celltype_adv_n_ref.detach().cpu()
        )
        details["weight/lambda_celltype_adv_active"] = float(celltype_adv_weight)
        details["weight/lambda_celltype_adv_target"] = float(self.lambda_celltype_adv)
        details["weight/lambda_celltype_adv_max"] = float(self.celltype_adv_max_lambda)
        details["weight/celltype_adv_warmup_epochs"] = float(
            self.celltype_adv_warmup_epochs
        )
        details["weight/celltype_adv_ramp_epochs"] = float(self.celltype_adv_ramp_epochs)

        # v15 auxiliary reference batch stats (meaningful balanced-acc on n≈
        # aux_batch_size reference cells; absent ⇒ 0/aux disabled). These are the
        # numbers the module gate cares about — n_ref here is the aux batch size,
        # not the ~3 reference cells in the main batch.
        aux_last = self._aux_celltype_adv_last or {}
        details["loss/celltype_adv_aux"] = float(aux_last.get("loss", 0.0))
        details["metric/celltype_adv_aux_balanced_acc"] = float(
            aux_last.get("balanced_accuracy", 0.0)
        )
        details["metric/celltype_adv_aux_acc"] = float(aux_last.get("accuracy", 0.0))
        details["metric/celltype_adv_aux_n_reference"] = float(
            aux_last.get("n_reference", 0.0)
        )
        details["metric/celltype_adv_aux_pretrain"] = float(aux_last.get("pretrain", 0.0))
        details["weight/celltype_adv_aux_grl_strength"] = float(
            aux_last.get("grl_strength", 0.0)
        )
        details["weight/celltype_adv_aux_enabled"] = float(
            1.0 if self._aux_celltype_adv_runner is not None else 0.0
        )

        details["weight/lambda_ref_center_db"] = float(self.lambda_ref_center_db)
        details["weight/lambda_tech_adv_active"] = float(tech_adv_weight)
        details["weight/lambda_tech_adv_target"] = float(self.lambda_tech_adv)
        details["weight/lambda_tech_adv_max"] = float(self.tech_adv_max_lambda)
        details["weight/tech_adv_warmup_epochs"] = float(
            self.tech_adv_warmup_epochs
        )
        details["weight/tech_adv_ramp_epochs"] = float(self.tech_adv_ramp_epochs)

        # --------------------------------------------------------------
        # 4. Explicit v3 diagnostic print
        # --------------------------------------------------------------
        is_train = bool(
            getattr(self.system, "training", False) and torch.is_grad_enabled()
        )
        # This counter controls DDP collectives below, so it must advance on
        # exactly the same calls on every rank.  In rank-0-only validation the
        # main rank executes additional loss forwards while the other ranks
        # wait.  Counting those eval forwards shifts rank 0's modulo cadence;
        # at the next training epoch rank 0 can then enter the diagnostic
        # all-reduce while the other ranks enter a different collective (for
        # example the PHU donor-by-celltype bank update), deadlocking the job.
        # Validation is diagnostic-only and must never affect this distributed
        # training clock.
        if is_train:
            self._v3_loss_print_counter += 1

        should_collect_for_print = (
            is_train
            and self.v3_diag_print_every > 0
            and self._v3_loss_print_counter % self.v3_diag_print_every == 0
            and self.current_epoch_index is not None
            and tech_adv_true_counts is not None
        )

        print_true_counts = tech_adv_true_counts
        print_pred_counts = tech_adv_pred_counts
        print_correct_counts = tech_adv_correct_counts
        print_recall_per_class = tech_adv_recall_per_class
        print_bal_acc = tech_adv_bal_acc.detach()

        if should_collect_for_print:
            if torch.distributed.is_available() and torch.distributed.is_initialized():
                print_true_counts = tech_adv_true_counts.to(
                    device=total.device,
                    dtype=total.dtype,
                ).clone()
                print_pred_counts = tech_adv_pred_counts.to(
                    device=total.device,
                    dtype=total.dtype,
                ).clone()
                print_correct_counts = tech_adv_correct_counts.to(
                    device=total.device,
                    dtype=total.dtype,
                ).clone()
                torch.distributed.all_reduce(print_true_counts, op=torch.distributed.ReduceOp.SUM)
                torch.distributed.all_reduce(print_pred_counts, op=torch.distributed.ReduceOp.SUM)
                torch.distributed.all_reduce(print_correct_counts, op=torch.distributed.ReduceOp.SUM)
                present = print_true_counts > 0
                print_recall_per_class = torch.zeros_like(print_true_counts)
                print_recall_per_class[present] = (
                    print_correct_counts[present] / print_true_counts[present]
                )
                print_bal_acc = (
                    print_recall_per_class[present].mean()
                    if present.any()
                    else torch.zeros((), dtype=total.dtype, device=total.device)
                )

        should_print = (
            is_train
            and self.v3_diag_print_every > 0
            and self._v3_loss_print_counter % self.v3_diag_print_every == 0
            and self._is_rank0()
        )

        if should_print:
            if print_true_counts is not None:
                true_counts_msg = ",".join(
                    str(int(x)) for x in print_true_counts.detach().cpu().tolist()
                )
                pred_counts_msg = ",".join(
                    str(int(x)) for x in print_pred_counts.detach().cpu().tolist()
                )
                recall_msg = ",".join(
                    f"{float(x):.3f}"
                    for x in print_recall_per_class.detach().cpu().tolist()
                )
                bal_acc_value = float(print_bal_acc.detach().cpu())
                collapse_msg, _ = self._tech_collapse_summary(
                    print_true_counts,
                    print_pred_counts,
                    print_recall_per_class,
                )
            else:
                true_counts_msg = "NA"
                pred_counts_msg = "NA"
                recall_msg = "NA"
                bal_acc_value = float(tech_adv_bal_acc.detach().cpu())
                collapse_msg = "NA"

            # Reference-center DB is off when its lambda is 0 (config
            # lambda_ref_center=0); show "disabled" instead of all-zero noise.
            if float(self.lambda_ref_center_db) > 0.0:
                ref_db_msg = (
                    f"loss={float(ref_center_db.detach().cpu()):.6f}  "
                    f"active_fraction={float(ref_center_db_active_fraction.detach().cpu()):.3f}  "
                    f"ref_celltypes={ref_center_db_n_ref_celltypes}  "
                    f"active_celltypes={ref_center_db_n_active_celltypes}  "
                    f"lambda={self.lambda_ref_center_db:.6f}"
                )
            else:
                ref_db_msg = "disabled"
            # v15 aux reference batch line: the celltype-adversary balanced-acc
            # and n_ref on the DEDICATED reference batch (n≈aux_batch_size, NOT
            # the ~3 reference cells in the main batch). "disabled" ⇒ v14 path.
            if self._aux_celltype_adv_runner is not None:
                aux_last = self._aux_celltype_adv_last or {}
                aux_mode = "pretrain" if float(aux_last.get("pretrain", 0.0)) > 0.5 else "adv"
                aux_adv_msg = (
                    f"mode={aux_mode}  "
                    f"loss={float(aux_last.get('loss', 0.0)):.6f}  "
                    f"balanced_acc={float(aux_last.get('balanced_accuracy', 0.0)):.3f}  "
                    f"n_ref={float(aux_last.get('n_reference', 0.0)):.0f}  "
                    f"grl={float(aux_last.get('grl_strength', 0.0)):.2f}  "
                    f"aux_batch_size={self.aux_celltype_adv_config.aux_batch_size}  "
                    f"aux_every={self.aux_celltype_adv_config.aux_every}"
                )
            else:
                aux_adv_msg = "disabled"
            print(
                f"[diagnostics | adversary | epoch {self.current_epoch_index}]\n"
                f"  Reference DB  : {ref_db_msg}\n"
                f"  Tech adversary: weight={tech_adv_weight:.6f}  "
                f"loss={float(tech_adv_loss.detach().cpu()):.6f}  "
                f"accuracy={float(tech_adv_acc.detach().cpu()):.3f}  "
                f"balanced_accuracy={bal_acc_value:.3f}  "
                f"target={self.lambda_tech_adv:.6f}  "
                f"max={self.tech_adv_max_lambda:.6f}  "
                f"ramp_epochs={self.tech_adv_ramp_epochs}\n"
                f"  Tech classes  : true_counts=[{true_counts_msg}]  "
                f"pred_counts=[{pred_counts_msg}]  "
                f"recall=[{recall_msg}]\n"
                f"  Tech collapse : {collapse_msg}\n"
                f"  CT adversary  : weight={celltype_adv_weight:.6f}  "
                f"loss={float(celltype_adv_loss.detach().cpu()):.6f}  "
                f"ref_balanced_acc={float(celltype_adv_bal_acc.detach().cpu()):.3f}  "
                f"n_reference={float(celltype_adv_n_ref.detach().cpu()):.0f}  "
                f"target={self.lambda_celltype_adv:.6f}  "
                f"max={self.celltype_adv_max_lambda:.6f}\n"
                f"  CT aux batch  : {aux_adv_msg}",
                flush=True,
            )

        # --------------------------------------------------------------
        # 5. Return same LossOutput type as v2, with v3 total/details
        # --------------------------------------------------------------
        return LossOutput(
            total=total,
            main=out.main,
            rec=out.rec,
            rec_weighted=out.rec_weighted,
            kl_state=out.kl_state,
            bam_kl=out.bam_kl,
            reg_tech=out.reg_tech,
            reg_gauge=out.reg_gauge,
            cls=out.cls,
            ref_center=ref_center_db.detach(),
            ref_state=out.ref_state,
            weights_per_cell=out.weights_per_cell,
            raw_uncertainty=out.raw_uncertainty,
            clipped_uncertainty=out.clipped_uncertainty,
            rec_per_cell=out.rec_per_cell,
            rec_weighted_per_cell=out.rec_weighted_per_cell,
            kl_state_per_cell=out.kl_state_per_cell,
            tech_score_mag=out.tech_score_mag,
            align=out.align,
            white=out.white,
            align_sex=out.align_sex,
            details=details,
        )
