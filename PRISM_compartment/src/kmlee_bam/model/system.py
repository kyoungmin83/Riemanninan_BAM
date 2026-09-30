from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn

try:
    from kmlee_bam.training.core_trainer import ModelForwardOutput
    from kmlee_bam.training.core_trainer import OrdinalBAMSystem
    from kmlee_bam.objectives.tech_adversary import ConditionalTechAdversary, TechAdversaryOutput
    from kmlee_bam.objectives.reference_celltype_adversary import (
        ReferenceCelltypeAdversary,
        ReferenceCelltypeAdversaryOutput,
    )
except ImportError:
    from kmlee_bam.training.core_trainer import ModelForwardOutput
    from kmlee_bam.training.core_trainer import OrdinalBAMSystem
    from kmlee_bam.objectives.tech_adversary import ConditionalTechAdversary, TechAdversaryOutput
    from kmlee_bam.objectives.reference_celltype_adversary import (
        ReferenceCelltypeAdversary,
        ReferenceCelltypeAdversaryOutput,
    )

from kmlee_bam.objectives.pathology_aux import PathologyAuxOutput
from kmlee_bam.objectives.sex_adversary import SexAdversary, SexAdversaryOutput
from kmlee_bam.objectives.metric_contrast import MetricContrastOutput
from kmlee_bam.objectives.cp_bam import CPBAMOutput


@dataclass
class V3ModelForwardOutput(ModelForwardOutput):
    """
    Forward output with v3 tensors declared as dataclass fields.

    DDP's unused-parameter traversal only sees registered container fields, so
    avoid attaching the adversary output dynamically after the base forward.
    """

    tech_adversary_out: Optional[TechAdversaryOutput] = None
    celltype_adversary_out: Optional[ReferenceCelltypeAdversaryOutput] = None
    pathology_aux_out: Optional[PathologyAuxOutput] = None
    cp_bam_out: Optional[CPBAMOutput] = None
    sex_adversary_out: Optional[SexAdversaryOutput] = None
    metric_contrast_out: Optional[MetricContrastOutput] = None


class V3OrdinalBAMSystem(OrdinalBAMSystem):
    """
    v3 system wrapper that keeps the v2 forward path and optionally attaches a
    conditional tech adversary to the DDP-visible module graph.

    It can also optionally attach a reference-only celltype adversary that pushes
    z_perp to NOT encode cell-type identity on reference/control cells (see
    :class:`ReferenceCelltypeAdversary`). Both adversaries are off unless
    explicitly constructed and passed in (None ⇒ byte-identical forward).
    """

    def __init__(
        self,
        *,
        base_system: OrdinalBAMSystem,
        tech_adversary: Optional[ConditionalTechAdversary] = None,
        celltype_adversary: Optional[ReferenceCelltypeAdversary] = None,
        sex_adversary: Optional[SexAdversary] = None,
    ) -> None:
        nn.Module.__init__(self)
        self.gene_embedding = base_system.gene_embedding
        self.module_tokenizer = base_system.module_tokenizer
        self.state_encoder = base_system.state_encoder
        self.prior = base_system.prior
        self.decoder = base_system.decoder
        self.classifier_head = base_system.classifier_head
        self.tech_adversary = tech_adversary
        self.celltype_adversary = celltype_adversary
        self.sex_adversary = sex_adversary
        # v20 multi-pathology auxiliary head — attached to base_system in build_system.
        # Must be carried onto the V3 wrapper so it lands in the optimizer (V3 is what's
        # optimized) AND is reachable as ``self.system.pathology_aux_head`` from the
        # trainer's loss path. Without this the head is silently dropped (never called,
        # never trained). None when v20 is disabled.
        self.pathology_aux_head = getattr(base_system, "pathology_aux_head", None)
        # CP-BAM head — same explicit-copy caveat as the aux head; carried so its params land in
        # the optimizer/DDP and it is reachable as self.system.cp_bam_head from the trainer.
        self.cp_bam_head = getattr(base_system, "cp_bam_head", None)
        self.metric_contrast_head = getattr(base_system, "metric_contrast_head", None)
        # PRISM E2E precision branch.  The wrapper is the object optimized and
        # wrapped by DDP, so an explicit copy is required just like pathology_aux.
        self.precision_head = getattr(base_system, "precision_head", None)
        # Joint generator-count components must remain registered on the V3
        # wrapper because this wrapper, not ``base_system``, is optimized and
        # wrapped by DDP.  They are absent for all legacy runs.
        self.generator_count_gate = getattr(
            base_system, "generator_count_gate", None
        )
        self.generator_count_objective = getattr(
            base_system, "generator_count_objective", None
        )
        self.generator_count_config = getattr(
            base_system, "generator_count_config", None
        )
        generator_epoch_state = getattr(
            base_system, "generator_count_epoch_state", None
        )
        if generator_epoch_state is not None:
            self.register_buffer(
                "generator_count_epoch_state",
                generator_epoch_state.detach().clone(),
                persistent=True,
            )
        self.generator_count_joint_training = bool(
            getattr(base_system, "generator_count_joint_training", False)
        )
        # v25 nuisance projector — same explicit-copy caveat as the aux head: the V3 wrapper
        # rebuilds submodules by hand, so carry it here or the base forward's projection step
        # (and its EMA buffers) would be silently dropped. None when v25 is disabled.
        self.nuisance_projector = getattr(base_system, "nuisance_projector", None)
        self._validate_component_compatibility()

    @staticmethod
    def _float_detail(x) -> float:
        try:
            return float(x.detach().cpu())
        except Exception:
            return float("nan")

    @staticmethod
    def _save_decoder_transient_state(decoder):
        names = (
            "_stash_diag",
            "_last_mixer_gates",
            "_last_pathology_p",
            "_last_pathdec_absmean",
            "_last_synergy_absmean",
            "_last_gate_absmean",
            "_last_gate_strengths",
        )
        sentinel = object()
        return names, {name: getattr(decoder, name, sentinel) for name in names}, sentinel

    @staticmethod
    def _restore_decoder_transient_state(decoder, names, saved, sentinel) -> None:
        for name in names:
            val = saved.get(name, sentinel)
            if val is sentinel:
                try:
                    delattr(decoder, name)
                except AttributeError:
                    pass
            else:
                setattr(decoder, name, val)

    def _cpbam_score_from_z(self, z, batch, decoder_out):
        tech_id = batch.get("tech_id", batch.get("batch_id"))
        if tech_id is None:
            raise KeyError("CP-BAM decoder losses require 'tech_id' or 'batch_id' in batch.")
        score, _base_score, _tech_score, state_score = self.decoder._compute_score_components(
            z_perp=z,
            celltype_id=batch["celltype_id"],
            tech_id=tech_id,
        )
        if getattr(decoder_out, "sex_score", None) is not None:
            score = score + decoder_out.sex_score
        return score, state_score

    def _attach_cp_bam_decoder_losses(self, cp_bam_out, base_out, batch):
        """Stage-2 CP-BAM decoder-aware losses.

        These losses make the common/private split decoder-native instead of a
        pure auxiliary projection. They are attached inside forward so all
        participating parameters remain visible to DDP. All weights default to
        zero, so existing runs are unchanged unless explicitly configured.
        """
        if cp_bam_out is None or self.cp_bam_head is None:
            return cp_bam_out
        cfg = getattr(self.cp_bam_head, "cfg", None)
        if cfg is None:
            return cp_bam_out
        lam_common = float(getattr(cfg, "lambda_common_decode", 0.0))
        lam_split = float(getattr(cfg, "lambda_split_full_decode", 0.0))
        lam_priv = float(getattr(cfg, "lambda_private_resid_score", 0.0))
        lam_ref = float(getattr(cfg, "lambda_ref_zero", 0.0))
        if lam_common <= 0.0 and lam_split <= 0.0 and lam_priv <= 0.0 and lam_ref <= 0.0:
            return cp_bam_out

        decoder_out = getattr(base_out, "decoder_out", None)
        z_common = getattr(cp_bam_out, "z_common", None)
        z_private = getattr(cp_bam_out, "z_private", None)
        if decoder_out is None or z_common is None or z_private is None:
            return cp_bam_out

        y_ord = batch.get("y_ord", None)
        total_extra = cp_bam_out.loss.new_zeros(())
        details = cp_bam_out.details
        names, saved, sentinel = self._save_decoder_transient_state(self.decoder)
        try:
            # Extra decode passes should not overwrite the diagnostics from the
            # main decoder forward that other loss terms/logs expect.
            self.decoder._stash_diag = False

            def nll_from_score(score):
                if y_ord is None:
                    return None
                _, probs = self.decoder._score_to_probs(score, decoder_out.thresholds)
                _, nll_cell = self.decoder._nll_from_probs(probs, y_ord.long())
                return nll_cell.mean()

            if lam_common > 0.0:
                score_common, state_common = self._cpbam_score_from_z(z_common, batch, decoder_out)
                loss_common = nll_from_score(score_common)
                if loss_common is not None:
                    total_extra = total_extra + lam_common * loss_common.to(dtype=total_extra.dtype)
                    details["loss/cpbam_common_decode_nll"] = self._float_detail(loss_common)
                    details["weight/lambda_cpbam_common_decode"] = float(lam_common)
                details["metric/cpbam_common_state_abs"] = self._float_detail(state_common.float().abs().mean())
            else:
                score_common = None
                state_common = None

            if lam_split > 0.0:
                score_split, state_split = self._cpbam_score_from_z(z_common + z_private, batch, decoder_out)
                loss_split = nll_from_score(score_split)
                if loss_split is not None:
                    total_extra = total_extra + lam_split * loss_split.to(dtype=total_extra.dtype)
                    details["loss/cpbam_split_full_decode_nll"] = self._float_detail(loss_split)
                    details["weight/lambda_cpbam_split_full_decode"] = float(lam_split)
                details["metric/cpbam_split_state_abs"] = self._float_detail(state_split.float().abs().mean())

            if lam_priv > 0.0:
                # Freeze the common branch target path here: z_private should fill
                # the residual needed to reproduce the original full score.
                score_priv, state_priv = self._cpbam_score_from_z(
                    z_common.detach() + z_private,
                    batch,
                    decoder_out,
                )
                target = decoder_out.score.detach()
                loss_priv = (score_priv.float() - target.float()).pow(2).mean()
                total_extra = total_extra + lam_priv * loss_priv.to(dtype=total_extra.dtype)
                details["loss/cpbam_private_resid_score"] = self._float_detail(loss_priv)
                details["weight/lambda_cpbam_private_resid_score"] = float(lam_priv)
                details["metric/cpbam_private_state_abs"] = self._float_detail(state_priv.float().abs().mean())

            if lam_ref > 0.0:
                is_reference = batch.get("is_reference", None)
                if is_reference is not None:
                    ref = is_reference.bool().view(-1)
                    n_ref = int(ref.sum().item())
                    if n_ref > 0:
                        score_c_ref, state_c_ref = self._cpbam_score_from_z(z_common, batch, decoder_out)
                        score_p_ref, state_p_ref = self._cpbam_score_from_z(z_private, batch, decoder_out)
                        loss_ref = (
                            z_common[ref].float().pow(2).mean()
                            + z_private[ref].float().pow(2).mean()
                            + state_c_ref[ref].float().pow(2).mean()
                            + state_p_ref[ref].float().pow(2).mean()
                        )
                        total_extra = total_extra + lam_ref * loss_ref.to(dtype=total_extra.dtype)
                        details["loss/cpbam_ref_zero"] = self._float_detail(loss_ref)
                        details["weight/lambda_cpbam_ref_zero"] = float(lam_ref)
                        details["metric/cpbam_ref_zero_n_ref"] = float(n_ref)
        finally:
            self._restore_decoder_transient_state(self.decoder, names, saved, sentinel)

        if total_extra.detach().abs().item() != 0.0:
            cp_bam_out.loss = cp_bam_out.loss + total_extra
            details["loss/cpbam_stage2_decoder_total"] = self._float_detail(total_extra)
            details["loss/cpbam_total"] = self._float_detail(cp_bam_out.loss)
        return cp_bam_out

    def forward(self, batch, *, compute_decoder: bool = True, **kwargs):
        out = super().forward(batch, compute_decoder=compute_decoder, **kwargs)
        # Encoder-only path (``compute_decoder=False``): used by the v15 aux
        # celltype-adversary step, which only needs encoder_out.z_s / mu_p /
        # sigma_p and recomputes its OWN σ-detached adversary call downstream
        # (see AuxReferenceAdversaryRunner). Attaching the tech/celltype
        # adversaries here on the σ-ATTACHED z_perp would be pure waste (and the
        # aux step never reads these outputs), so skip them too — making the aux
        # forward genuinely encoder-only. The normal path (compute_decoder=True)
        # is byte-identical.
        tech_adversary_out = None
        celltype_adversary_out = None
        pathology_aux_out = None
        cp_bam_out = None
        sex_adversary_out = None
        metric_contrast_out = None
        if compute_decoder:
            if self.tech_adversary is not None:
                tech_id = batch.get("tech_id", batch.get("batch_id"))
                if tech_id is None:
                    raise KeyError("V3 tech adversary requires 'tech_id' or 'batch_id' in batch.")
                tech_adversary_out = self.tech_adversary(
                    out.z_perp,
                    batch["celltype_id"],
                    tech_id,
                )
            if self.celltype_adversary is not None:
                if "is_reference" not in batch:
                    raise KeyError(
                        "reference_celltype_adversary is enabled but the batch has no "
                        "'is_reference' mask; the reference-only adversary loss requires it. "
                        "Disable the adversary or provide is_reference in the dataset."
                    )
                celltype_adversary_out = self.celltype_adversary(
                    out.z_perp,
                    batch["celltype_id"],
                    batch["is_reference"],
                )
            if self.sex_adversary is not None:
                # SOFT sex erasure on the CLEAN (post-projection) z. sex_id: 0=F,1=M,-1=missing.
                sex_adversary_out = self.sex_adversary(out.z_perp, batch.get("sex_id"))
            # v20 multi-pathology auxiliary head — called INSIDE forward so its params
            # are tracked by the DDP reducer (calling it from the trainer's loss path
            # raised "mark a variable ready only once" on every rank). Returns its loss
            # in the output; the trainer only scales + adds it.
            if (
                self.pathology_aux_head is not None
                and getattr(self.pathology_aux_head, "cfg", None) is not None
                and self.pathology_aux_head.cfg.enabled
                and bool(getattr(self, "pathology_aux_forward_enabled", True))
                and "donor_id" in batch
                and out.z_perp is not None
            ):
                _axes = self.pathology_aux_head.axes
                if all(ax in batch for ax in _axes):
                    pathology_aux_out = self.pathology_aux_head(
                        z_perp=out.z_perp,
                        donor_id=batch["donor_id"],
                        celltype_id=batch["celltype_id"],
                        labels={ax: batch[ax] for ax in _axes},
                        resid_target=batch.get("axis_resid"),       # v23 (None unless emitted)
                        resid_valid=batch.get("axis_resid_valid"),
                    )
            # CP-BAM consensus-private decomposition — called INSIDE forward (like pathology_aux)
            # so its params are tracked by the DDP reducer; the trainer only adds cp_bam_out.loss.
            if (
                self.cp_bam_head is not None
                and getattr(self.cp_bam_head, "cfg", None) is not None
                and self.cp_bam_head.cfg.enabled
                and out.z_perp is not None
                and "donor_id" in batch
            ):
                _cax = self.cp_bam_head.axes
                if all(ax in batch for ax in _cax):
                    cp_bam_out = self.cp_bam_head(
                        out.z_perp,
                        batch["donor_id"],
                        batch["celltype_id"],
                        {ax: batch[ax] for ax in _cax},
                    )
                    cp_bam_out = self._attach_cp_bam_decoder_losses(cp_bam_out, out, batch)
            # v35 metric-contrast — pressure for amyloid-aligned curvature. Called INSIDE forward
            # (like pathology_aux) so the gate/decoder grads are DDP-tracked. Default off => skipped.
            if (
                self.metric_contrast_head is not None
                and getattr(self.metric_contrast_head, "cfg", None) is not None
                and self.metric_contrast_head.cfg.enabled
                and out.z_perp is not None
                and "donor_id" in batch
                and self.metric_contrast_head.cfg.amyloid_label_key in batch
            ):
                _mtid = batch.get("tech_id", batch.get("batch_id"))
                if _mtid is not None:
                    metric_contrast_out = self.metric_contrast_head(
                        decoder=self.decoder,
                        z_perp=out.z_perp,
                        celltype_id=batch["celltype_id"],
                        tech_id=_mtid,
                        donor_id=batch["donor_id"],
                        amyloid_label=batch[self.metric_contrast_head.cfg.amyloid_label_key],
                    )
        return V3ModelForwardOutput(
            tokens=out.tokens,
            token_padding_mask=out.token_padding_mask,
            gene_tokens=out.gene_tokens,
            module_tokens=out.module_tokens,
            pooled_gene_state=out.pooled_gene_state,
            module_activity=out.module_activity,
            module_attention=out.module_attention,
            encoder_out=out.encoder_out,
            mu_p=out.mu_p,
            logvar_p=out.logvar_p,
            sigma_p=out.sigma_p,
            z_perp=out.z_perp,
            decoder_out=out.decoder_out,
            cls_logits=out.cls_logits,
            precision_out=out.precision_out,
            tech_adversary_out=tech_adversary_out,
            celltype_adversary_out=celltype_adversary_out,
            pathology_aux_out=pathology_aux_out,
            cp_bam_out=cp_bam_out,
            sex_adversary_out=sex_adversary_out,
            metric_contrast_out=metric_contrast_out,
        )
