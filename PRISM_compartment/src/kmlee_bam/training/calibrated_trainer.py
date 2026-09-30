from __future__ import annotations

from typing import Dict, Optional

import torch

try:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.adversarial_trainer import V3Trainer
    from kmlee_bam.objectives.ordinal_uncertainty import (
        GlobalReferenceBankConfig,
        GlobalReferenceCenterBank,
        OrdinalBalanceConfig,
        UncertaintyCalibrationConfig,
        ordinal_balance_losses,
        uncertainty_calibration_losses,
    )
    from kmlee_bam.objectives.hierarchical_ordinal import (
        HierarchicalOrdinalConfig,
        compute_hierarchical_ordinal_loss,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput
except ImportError:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.adversarial_trainer import V3Trainer
    from kmlee_bam.objectives.ordinal_uncertainty import (
        GlobalReferenceBankConfig,
        GlobalReferenceCenterBank,
        OrdinalBalanceConfig,
        UncertaintyCalibrationConfig,
        ordinal_balance_losses,
        uncertainty_calibration_losses,
    )
    from kmlee_bam.objectives.hierarchical_ordinal import (
        HierarchicalOrdinalConfig,
        compute_hierarchical_ordinal_loss,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput


def _resolve_tech_id(batch: Dict[str, torch.Tensor]) -> torch.Tensor:
    if "tech_id" in batch:
        return batch["tech_id"]
    if "batch_id" in batch:
        return batch["batch_id"]
    raise KeyError("Batch must contain either 'tech_id' or 'batch_id'.")


class V4Trainer(V3Trainer):
    """
    v4 trainer: v3 stability fixes plus ordinal/depth/reference calibration.

    The base v2/v3 loss remains intact. v4 adds auxiliary terms that directly
    address the post-hoc failures found in validation:
    - sparse ordinal bins were hidden by bin-0 majority accuracy
    - reference origin was too batch-local
    - BAM entropy was high but not calibrated as reliability
    """

    def __init__(
        self,
        *args,
        ordinal_class_weights: Optional[torch.Tensor] = None,
        ordinal_balance_config: Optional[OrdinalBalanceConfig] = None,
        global_ref_config: Optional[GlobalReferenceBankConfig] = None,
        uncertainty_calibration_config: Optional[UncertaintyCalibrationConfig] = None,
        hierarchical_ordinal_config: Optional[HierarchicalOrdinalConfig] = None,
        n_celltypes: Optional[int] = None,
        d_z: Optional[int] = None,
        v4_diag_print_every: int = 100,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)

        self.ordinal_balance_config = ordinal_balance_config or OrdinalBalanceConfig(
            enabled=False
        )
        # Hierarchical ordinal (L1 zero/nonzero + L2 group | nonzero + L3 EMD).
        # See doc/hierarchical_ordinal_loss.md.
        self.hierarchical_ordinal_config = (
            hierarchical_ordinal_config or HierarchicalOrdinalConfig(enabled=False)
        )
        if ordinal_class_weights is None:
            ordinal_class_weights = torch.ones(1, dtype=torch.float32)
        self.ordinal_class_weights = ordinal_class_weights.detach().float()

        self.global_ref_config = global_ref_config or GlobalReferenceBankConfig(
            enabled=False
        )
        self.global_reference_bank: Optional[GlobalReferenceCenterBank] = None
        if (
            self.global_ref_config.enabled
            and float(self.global_ref_config.lambda_global) > 0.0
            and n_celltypes is not None
            and d_z is not None
        ):
            self.global_reference_bank = GlobalReferenceCenterBank(
                n_celltypes=int(n_celltypes),
                d_z=int(d_z),
                config=self.global_ref_config,
                device=self.device,
                dtype=torch.float32,
            )

        self.uncertainty_calibration_config = (
            uncertainty_calibration_config
            or UncertaintyCalibrationConfig(enabled=False)
        )
        self.v4_diag_print_every = max(0, int(v4_diag_print_every))
        self._v4_loss_print_counter = 0

    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        out = super()._compute_loss(batch, model_out)
        total = out.total
        details = dict(out.details)

        # Keep auxiliary tensors on the same dtype/device as the active loss.
        dtype = total.dtype
        device = total.device
        y_ord = batch["y_ord"].long()
        tech_id = _resolve_tech_id(batch)

        # --------------------------------------------------------------
        # 1. Ordinal-bin balanced reconstruction pressure
        # --------------------------------------------------------------
        ordinal_out = ordinal_balance_losses(
            nll_per_gene=model_out.decoder_out.nll_per_gene,
            probs=model_out.decoder_out.probs,
            y_ord=y_ord,
            class_weights=self.ordinal_class_weights.to(device=device),
            config=self.ordinal_balance_config,
        )
        if self.ordinal_balance_config.enabled:
            total = (
                total
                + float(self.ordinal_balance_config.lambda_balanced)
                * ordinal_out.loss_balanced
                + float(self.ordinal_balance_config.lambda_nonzero)
                * ordinal_out.loss_nonzero
            )

        # --------------------------------------------------------------
        # 1b. Hierarchical ordinal loss (L1 zero/nonzero + L2 group | nonzero
        #     + L3 EMD).  See doc/hierarchical_ordinal_loss.md.
        # --------------------------------------------------------------
        epoch_hint = getattr(self, "current_epoch_index", None)
        # π_tech soft-zero conduit: a top-of-MRO trainer (V8Trainer) may stash a
        # per-(cell,gene) technical-dropout probability on `_pending_pi_tech`
        # just before calling super()._compute_loss. Absent (every non-V8
        # trainer) ⇒ None ⇒ the hierarchical loss is byte-identical.
        pending_pi_tech = getattr(self, "_pending_pi_tech", None)
        hier_out = compute_hierarchical_ordinal_loss(
            probs=model_out.decoder_out.probs,
            y_true=y_ord,
            config=self.hierarchical_ordinal_config,
            epoch=epoch_hint,
            pi_tech=pending_pi_tech,
        )
        if self.hierarchical_ordinal_config.enabled:
            total = total + hier_out.weighted_total.to(dtype=dtype)

        # --------------------------------------------------------------
        # 2. Global reference-origin memory bank
        # --------------------------------------------------------------
        global_ref_loss = torch.zeros((), dtype=dtype, device=device)
        global_ref_active_fraction = torch.zeros((), dtype=dtype, device=device)
        global_ref_center_norm_mean = torch.zeros((), dtype=dtype, device=device)
        global_ref_center_norm_max = torch.zeros((), dtype=dtype, device=device)
        global_ref_reliability_mean = torch.zeros((), dtype=dtype, device=device)
        global_ref_n_active = 0

        if self.global_reference_bank is not None:
            self.global_reference_bank.to(device=device, dtype=torch.float32)
            ref_out = self.global_reference_bank.loss(
                z_perp=model_out.z_perp.float(),
                celltype_id=batch["celltype_id"],
                is_reference=batch.get("is_reference", None),
            )
            global_ref_loss = ref_out.loss.to(dtype=dtype)
            global_ref_active_fraction = ref_out.active_fraction.to(dtype=dtype)
            global_ref_center_norm_mean = ref_out.center_norm_mean.to(dtype=dtype)
            global_ref_center_norm_max = ref_out.center_norm_max.to(dtype=dtype)
            global_ref_reliability_mean = ref_out.reliability_mean.to(dtype=dtype)
            global_ref_n_active = int(ref_out.n_active_celltypes)
            total = total + float(self.global_ref_config.lambda_global) * global_ref_loss

            is_train = bool(getattr(self.system, "training", False))
            if is_train:
                self.global_reference_bank.update(
                    z_perp=model_out.z_perp.float(),
                    celltype_id=batch["celltype_id"],
                    is_reference=batch.get("is_reference", None),
                )

        # --------------------------------------------------------------
        # 3. BAM uncertainty calibration
        # --------------------------------------------------------------
        rec_per_cell_for_grad = model_out.decoder_out.nll_per_cell
        unc_out = uncertainty_calibration_losses(
            raw_uncertainty=model_out.encoder_out.cell_uncertainty,
            rec_per_cell=rec_per_cell_for_grad,
            celltype_id=batch["celltype_id"],
            tech_id=tech_id,
            depth_value=batch.get("depth_value", None),
            config=self.uncertainty_calibration_config,
        )
        if self.uncertainty_calibration_config.enabled:
            total = (
                total
                + float(self.uncertainty_calibration_config.lambda_rec_alignment)
                * unc_out.rec_alignment
                + float(self.uncertainty_calibration_config.lambda_relative_rec)
                * unc_out.relative_rec
            )

        # --------------------------------------------------------------
        # 4. Details and explicit diagnostics
        # --------------------------------------------------------------
        details["loss/ordinal_balanced"] = float(
            ordinal_out.loss_balanced.detach().cpu()
        )
        details["loss/ordinal_nonzero"] = float(
            ordinal_out.loss_nonzero.detach().cpu()
        )
        details["metric/ordinal_balanced_recall"] = float(
            ordinal_out.balanced_recall.detach().cpu()
        )
        details["metric/ordinal_nonzero_acc"] = float(
            ordinal_out.nonzero_accuracy.detach().cpu()
        )
        details["metric/ordinal_nonzero_within1"] = float(
            ordinal_out.nonzero_within1.detach().cpu()
        )
        details["metric/ordinal_true_zero_frac"] = float(
            ordinal_out.true_zero_fraction.detach().cpu()
        )
        details["metric/ordinal_pred_zero_frac"] = float(
            ordinal_out.pred_zero_fraction.detach().cpu()
        )
        details["weight/lambda_ordinal_balanced"] = float(
            self.ordinal_balance_config.lambda_balanced
        )
        details["weight/lambda_ordinal_nonzero"] = float(
            self.ordinal_balance_config.lambda_nonzero
        )

        # Hierarchical ordinal diagnostics (always logged; weights tell whether they were applied)
        details["loss/hier_zero"] = float(hier_out.loss_zero.detach().cpu())
        details["loss/hier_group"] = float(hier_out.loss_group.detach().cpu())
        details["loss/hier_emd"] = float(hier_out.loss_emd.detach().cpu())
        details["loss/hier_within_group"] = float(hier_out.loss_within_group.detach().cpu())
        details["loss/hier_weighted_total"] = float(
            hier_out.weighted_total.detach().cpu()
        )
        details["metric/hier_pred_zero_frac"] = float(
            hier_out.pred_zero_frac.detach().cpu()
        )
        details["metric/hier_n_nonzero_pairs"] = float(hier_out.n_nonzero_pairs)
        details["metric/hier_ramp_scale"] = float(hier_out.ramp_scale_used)
        # Asymmetric under-prediction term + its required monitoring (gap,
        # high-tier exp_level, nonzero→zero leak, zero_pred vs zero_true).
        # Present only when lambda_under_asym > 0.
        if hier_out.loss_under_asym is not None:
            details["loss/hier_under_asym"] = float(hier_out.loss_under_asym.detach().cpu())
        if hier_out.asym_diag:
            # Diagnostics arrive as detached tensors; convert them in ONE batched
            # GPU→CPU sync (not per-key float(), which would sync 6× per step).
            _ak = list(hier_out.asym_diag.keys())
            _av = torch.stack([hier_out.asym_diag[k].float().reshape(()) for k in _ak]).cpu().tolist()
            details.update(dict(zip(_ak, _av)))
        for gi in range(len(self.hierarchical_ordinal_config.group_definition)):
            details[f"metric/hier_pred_group{gi}_mean"] = float(
                hier_out.pred_group_dist[gi].detach().cpu()
            )
            details[f"metric/hier_true_group{gi}_frac"] = float(
                hier_out.true_group_dist[gi].detach().cpu()
            )
        details["weight/lambda_hier_zero"] = float(
            self.hierarchical_ordinal_config.lambda_zero
        )
        details["weight/lambda_hier_group"] = float(
            self.hierarchical_ordinal_config.lambda_group
        )
        details["weight/lambda_hier_emd"] = float(
            self.hierarchical_ordinal_config.lambda_emd
        )
        details["weight/lambda_hier_within_group"] = float(
            self.hierarchical_ordinal_config.lambda_within_group
        )
        details["weight/hier_zero_pos_weight"] = float(
            self.hierarchical_ordinal_config.zero_pos_weight
        )

        details["loss/ref_center_global_bank"] = float(
            global_ref_loss.detach().cpu()
        )
        details["metric/ref_global_active_fraction"] = float(
            global_ref_active_fraction.detach().cpu()
        )
        details["metric/ref_global_center_norm_mean"] = float(
            global_ref_center_norm_mean.detach().cpu()
        )
        details["metric/ref_global_center_norm_max"] = float(
            global_ref_center_norm_max.detach().cpu()
        )
        details["metric/ref_global_reliability_mean"] = float(
            global_ref_reliability_mean.detach().cpu()
        )
        details["metric/ref_global_n_active_celltypes"] = float(global_ref_n_active)
        details["weight/lambda_ref_center_global_bank"] = float(
            self.global_ref_config.lambda_global
        )

        details["loss/unc_rec_alignment"] = float(
            unc_out.rec_alignment.detach().cpu()
        )
        details["loss/unc_relative_rec"] = float(unc_out.relative_rec.detach().cpu())
        details["unc/relative_std"] = float(unc_out.rel_unc_std.detach().cpu())
        details["weight/bam_rel_w_min"] = float(unc_out.rel_w_min.detach().cpu())
        details["weight/bam_rel_w_mean"] = float(unc_out.rel_w_mean.detach().cpu())
        details["weight/bam_rel_w_max"] = float(unc_out.rel_w_max.detach().cpu())
        details["weight/lambda_unc_rec_alignment"] = float(
            self.uncertainty_calibration_config.lambda_rec_alignment
        )
        details["weight/lambda_unc_relative_rec"] = float(
            self.uncertainty_calibration_config.lambda_relative_rec
        )
        details["loss/total"] = float(total.detach().cpu())

        is_train = bool(getattr(self.system, "training", False))

        self._v4_loss_print_counter += 1
        if (
            is_train
            and self.v4_diag_print_every > 0
            and self._v4_loss_print_counter % self.v4_diag_print_every == 0
            and self._is_rank0()
        ):
            pred_group_msg = ", ".join(
                f"{details.get(f'metric/hier_pred_group{gi}_mean', float('nan')):.3f}"
                for gi in range(len(self.hierarchical_ordinal_config.group_definition))
            )
            true_group_msg = ", ".join(
                f"{details.get(f'metric/hier_true_group{gi}_frac', float('nan')):.3f}"
                for gi in range(len(self.hierarchical_ordinal_config.group_definition))
            )
            # Global reference bank is off in this run (enabled=false/lambda=0 ⇒
            # global_reference_bank is None). Show "disabled" instead of zeros.
            if self.global_reference_bank is not None:
                ref_bank_msg = (
                    f"loss={float(global_ref_loss.detach().cpu()):.6f}  "
                    f"reliability={float(global_ref_reliability_mean.detach().cpu()):.3f}  "
                    f"center_norm={float(global_ref_center_norm_mean.detach().cpu()):.3f}"
                )
            else:
                ref_bank_msg = "disabled"
            print(
                f"[diagnostics | ordinal | epoch {self.current_epoch_index}]\n"
                f"  Ordinal metrics : balanced_recall={float(ordinal_out.balanced_recall.detach().cpu()):.3f}  "
                f"nonzero_exact={float(ordinal_out.nonzero_accuracy.detach().cpu()):.3f}  "
                f"nonzero_within1={details.get('metric/ordinal_nonzero_within1', float('nan')):.3f}  "
                f"zero_true={float(ordinal_out.true_zero_fraction.detach().cpu()):.3f}  "
                f"zero_pred={float(ordinal_out.pred_zero_fraction.detach().cpu()):.3f}\n"
                f"  Ordinal losses  : balanced={float(ordinal_out.loss_balanced.detach().cpu()):.6f}  "
                f"nonzero={float(ordinal_out.loss_nonzero.detach().cpu()):.6f}  "
                f"lambda_balanced={self.ordinal_balance_config.lambda_balanced:.4f}  "
                f"lambda_nonzero={self.ordinal_balance_config.lambda_nonzero:.4f}\n"
                f"  Hierarchical    : enabled={int(self.hierarchical_ordinal_config.enabled)}  "
                f"weighted={details['loss/hier_weighted_total']:.6f}  "
                f"zero={details['loss/hier_zero']:.6f}  "
                f"group={details['loss/hier_group']:.6f}  "
                f"EMD={details['loss/hier_emd']:.6f}  "
                f"within_group={details['loss/hier_within_group']:.6f}  "
                f"p_zero={details['metric/hier_pred_zero_frac']:.3f}  "
                f"ramp={details['metric/hier_ramp_scale']:.2f}\n"
                f"  Hier groups     : pred_group=[{pred_group_msg}]  "
                f"true_group=[{true_group_msg}]  "
                f"nonzero_pairs={details['metric/hier_n_nonzero_pairs']:.0f}\n"
                f"  Global ref bank : {ref_bank_msg}\n"
                f"  Uncertainty     : rec_alignment={float(unc_out.rec_alignment.detach().cpu()):.6f}  "
                f"relative_std={float(unc_out.rel_unc_std.detach().cpu()):.3f}  "
                f"relative_weight=[{float(unc_out.rel_w_min.detach().cpu()):.3f}, "
                f"{float(unc_out.rel_w_mean.detach().cpu()):.3f}, "
                f"{float(unc_out.rel_w_max.detach().cpu()):.3f}]",
                flush=True,
            )

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
            ref_center=out.ref_center,
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
