from __future__ import annotations

from typing import Dict, Optional

import torch

try:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.calibrated_trainer import V4Trainer, _resolve_tech_id
    from kmlee_bam.objectives.state_usage import (
        V5StateUsageConfig,
        v5_state_usage_losses,
    )
    from kmlee_bam.objectives.uncertainty_residual import (
        V5UncertaintyResidualConfig,
        v5_uncertainty_residual_losses,
    )
    from kmlee_bam.model.reference_anchored_prior import (
        ReferenceAnchoredPriorConfig,
        ReferencePriorAnchor,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput
except ImportError:
    from kmlee_bam.objectives.total_loss import LossOutput
    from kmlee_bam.training.calibrated_trainer import V4Trainer, _resolve_tech_id
    from kmlee_bam.objectives.state_usage import (
        V5StateUsageConfig,
        v5_state_usage_losses,
    )
    from kmlee_bam.objectives.uncertainty_residual import (
        V5UncertaintyResidualConfig,
        v5_uncertainty_residual_losses,
    )
    from kmlee_bam.model.reference_anchored_prior import (
        ReferenceAnchoredPriorConfig,
        ReferencePriorAnchor,
    )
    from kmlee_bam.training.core_trainer import ModelForwardOutput


class V6Trainer(V4Trainer):
    """
    v6 trainer: v4/v5 diagnostics with a reference-anchored prior mean.

    The key difference from v5 is that `mu_p(c)` is not corrected through an
    extra loss on z_perp. It is directly refreshed from reference-normal cells
    using donor-balanced EMA of encoder posterior means. This makes the local
    origin a definition, not a tug-of-war between KL and a reference penalty.
    """

    def __init__(
        self,
        *args,
        v6_reference_prior_config: Optional[ReferenceAnchoredPriorConfig] = None,
        v6_state_usage_config: Optional[V5StateUsageConfig] = None,
        v6_uncertainty_config: Optional[V5UncertaintyResidualConfig] = None,
        # v6_min_mixer additions (see doc/v6_min_mixer_design_2026-05-22.md).
        v6_lambda_gate_balance: float = 0.0,
        v6_lambda_gate_ref_safe: float = 0.0,
        v6_lambda_coeff_zero_mean: float = 0.0,
        n_celltypes: Optional[int] = None,
        n_donors: Optional[int] = None,
        d_z: Optional[int] = None,
        v6_diag_print_every: int = 100,
        **kwargs,
    ) -> None:
        super().__init__(*args, n_celltypes=n_celltypes, d_z=d_z, **kwargs)

        self.v6_reference_prior_config = (
            v6_reference_prior_config or ReferenceAnchoredPriorConfig(enabled=False)
        )
        self.v6_reference_prior: Optional[ReferencePriorAnchor] = None
        if (
            self.v6_reference_prior_config.enabled
            and n_celltypes is not None
            and n_donors is not None
            and d_z is not None
        ):
            self.v6_reference_prior = ReferencePriorAnchor(
                system=self.system,
                n_celltypes=int(n_celltypes),
                n_donors=int(n_donors),
                d_z=int(d_z),
                config=self.v6_reference_prior_config,
            )
            self.v6_reference_prior.freeze_prior_mu()

        self.v6_state_usage_config = v6_state_usage_config or V5StateUsageConfig(
            enabled=False
        )
        self.v6_uncertainty_config = (
            v6_uncertainty_config or V5UncertaintyResidualConfig(enabled=False)
        )
        # v6_min_mixer: safety regularizers on the score residual mixer and
        # on the Lie coefficient head. All default 0 so v5/v6_minimal runs
        # stay byte-for-byte identical when the mixer is disabled.
        self.v6_lambda_gate_balance = float(v6_lambda_gate_balance)
        self.v6_lambda_gate_ref_safe = float(v6_lambda_gate_ref_safe)
        self.v6_lambda_coeff_zero_mean = float(v6_lambda_coeff_zero_mean)
        self.v6_diag_print_every = max(0, int(v6_diag_print_every))
        self._v6_loss_print_counter = 0

    def set_ddp_control_group(self, group) -> None:
        if self.v6_reference_prior is not None:
            self.v6_reference_prior.set_sync_process_group(group)

    def train_epoch(self, loader, *args, epoch_index: Optional[int] = None, **kwargs):
        refresh_every = int(self.v6_reference_prior_config.refresh_every_epochs)
        should_refresh = (
            self.v6_reference_prior is not None
            and bool(self.v6_reference_prior_config.refresh_before_epoch)
            and refresh_every > 0
            and epoch_index is not None
            and (int(epoch_index) == 1 or (int(epoch_index) - 1) % refresh_every == 0)
        )
        if should_refresh:
            summary = self.v6_reference_prior.refresh_from_loader(
                loader,
                self,
                epoch_index=epoch_index,
            )
            if summary is not None and self._is_rank0():
                print(
                    "[kmlee-ref-refresh] "
                    f"epoch={epoch_index} "
                    f"active_ct={summary.n_active_celltypes} "
                    f"ref_cells={summary.n_total_reference_cells} "
                    f"ref_donors={summary.n_total_reference_donors} "
                    f"update_norm={summary.mean_update_norm:.6f} "
                    f"prior_norm_mean={summary.mean_prior_norm:.6f} "
                    f"prior_norm_max={summary.max_prior_norm:.6f}",
                    flush=True,
                )

        return super().train_epoch(loader, *args, epoch_index=epoch_index, **kwargs)

    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        out = super()._compute_loss(batch, model_out)
        total = out.total
        details = dict(out.details)

        # --------------------------------------------------------------
        # 1. State branch usage diagnostics/guards
        # --------------------------------------------------------------
        state_out = v5_state_usage_losses(
            decoder_out=model_out.decoder_out,
            config=self.v6_state_usage_config,
        )
        if self.v6_state_usage_config.enabled:
            total = (
                total
                + float(self.v6_state_usage_config.lambda_state_abs)
                * state_out.loss_state_abs
                + float(self.v6_state_usage_config.lambda_state_fraction)
                * state_out.loss_state_fraction
                + float(self.v6_state_usage_config.lambda_tech_to_state)
                * state_out.loss_tech_to_state
            )

        # --------------------------------------------------------------
        # 1b. v20 multi-pathology auxiliary head (donor-balanced aggregate)
        #     Gives z explicit, supervised reasons to represent the partially-
        #     independent disease axes (Braak/Thal/LATE/Lewy) so it stops collapsing
        #     to the single AD-severity direction. Aggregated to (donor,celltype) via
        #     the head's EMA bank to avoid donor-label pseudo-replication.
        # --------------------------------------------------------------
        path_head = getattr(self.system, "pathology_aux_head", None)
        path_out = getattr(model_out, "pathology_aux_out", None)
        path_curriculum_scale = float(
            getattr(self, "pathology_aux_curriculum_scale", 1.0)
        )
        if (
            path_out is not None
            and path_head is not None
            and float(path_head.cfg.lambda_path_aux) > 0.0
        ):
            # The head was CALLED inside the model forward (V3 system) so its params are
            # DDP-tracked; here we only scale + add the returned loss + log diagnostics.
            total = total + (
                path_curriculum_scale
                * float(path_head.cfg.lambda_path_aux)
                * path_out.loss
            )
            details["weight/pathology_curriculum_scale"] = path_curriculum_scale
            details["loss/v20_path_aux"] = float(path_out.loss.detach())
            details["metric/v20_path_groups_sup"] = float(path_out.n_groups_supervised)
            details["metric/v20_path_axes_active"] = float(path_out.n_axes_active)
            details["metric/v20_specificity"] = float(getattr(path_out, "specificity_loss", 0.0))
            details["metric/v20_residual_frac"] = float(getattr(path_out, "residual_frac", 0.0))
            if float(getattr(path_head, "lambda_resid_axis", 0.0)) > 0.0:   # v23 resid-axis: log ONLY when ON
                details["metric/v20_resid_axis"] = float(getattr(path_out, "resid_axis_loss", 0.0))   # else 0.0≠"잘 맞힘"
            for _ax, _l in path_out.per_axis_loss.items():
                details[f"loss/v20_{_ax}"] = _l
            for _ax, _pf in path_out.per_axis_pos_frac.items():
                details[f"metric/v20_posfrac_{_ax}"] = _pf

        # v23 INTERACTION usage: L2 norm of the (zero-init) decoder interaction patterns.
        # Grows from 0 as the decoder learns to use the pathology-axis PRODUCTS — i.e. how
        # much curvature is actually being exercised. None/absent when interactions are off.
        try:
            _ip = getattr(getattr(self.system, "decoder", None), "interaction_patterns", None)
            if _ip is not None:
                details["metric/v23_interaction_norm"] = float(_ip.detach().norm())
        except Exception:
            pass

        # v23 (B) SOFT-TIE: cosine-align the decoder's free interaction coords to the
        # pathology head's SUPERVISED axis DIRECTIONS, so the interaction products become
        # named-pathology (Braak/Thal/..) and curvature lands in DISEASE directions (not
        # free/noise). Off (lambda 0) or interactions-off => skipped. Direction-only
        # (cosine) so the coord MAGNITUDE stays free (absorbed by interaction_patterns).
        try:
            _lam_tie = path_curriculum_scale * (float(getattr(path_head.cfg, "lambda_interaction_tie", 0.0))
                        if path_head is not None and getattr(path_head, "cfg", None) is not None else 0.0)
            if _lam_tie > 0.0:
                _proj = getattr(getattr(self.system, "decoder", None), "interaction_proj", None)
                _hh = getattr(path_head, "head", None)
                if _proj is not None and isinstance(_hh, torch.nn.Linear):
                    _m = min(int(_proj.weight.shape[0]), int(_hh.weight.shape[0]))
                    _cos = torch.nn.functional.cosine_similarity(
                        _proj.weight[:_m], _hh.weight[:_m].detach(), dim=1)
                    _tie = (1.0 - _cos).mean()
                    total = total + _lam_tie * _tie
                    details["metric/v23_interaction_tie"] = float(_tie.detach())
        except Exception:
            pass

        # v35 DISEASE-GATE SPARSITY: L1 on the decoder disease-gate per-axis strengths s_k so the
        # model SELF-SELECTS which pathology axes get curvature (dead axes shrink to 0). Logs each
        # s_k (watch: only the signal axis should grow; a dead axis growing = memorization => raise
        # lambda_gate_sparse) + the gate magnitude. lambda 0 / gate off => skipped. Fully guarded.
        try:
            _lam_gs = path_curriculum_scale * (float(getattr(path_head.cfg, "lambda_gate_sparse", 0.0))
                       if path_head is not None and getattr(path_head, "cfg", None) is not None else 0.0)
            _dec_g = getattr(self.system, "decoder", None)
            if getattr(_dec_g, "use_disease_gate", False):
                # review finding #1: penalize the EFFECTIVE per-axis magnitude s_k·‖disease_base_k‖
                # (group-lasso), NOT s_k alone — else the model evades L1 by (s↓, base↑) rescaling and
                # the sparsity is toothless. This drives WHOLE dead axes (s_k AND base_k) to zero =
                # genuine self-selection.
                _s = torch.nn.functional.softplus(_dec_g.gate_strength_raw)             # [na] ≥0
                _bn = _dec_g.disease_base.norm(dim=1)                                   # [na] ‖base_k‖
                if getattr(_dec_g, "gate_celltype_conditioned", False) and getattr(_dec_g, "gate_rank", 0) > 0:
                    _bn = _bn + _dec_g.disease_U.norm(dim=(0, 2)) * _dec_g.disease_V.norm()
                _eff = _s * _bn                                                          # [na] EFFECTIVE strength
                if _lam_gs > 0.0:
                    total = total + _lam_gs * _eff.sum()                                # scale-invariant group penalty
                # log even when lambda==0 so the gate can be monitored during a #1-only run
                details["metric/v35_gate_absmean"] = float(getattr(_dec_g, "_last_gate_absmean", 0.0) or 0.0)
                details["metric/v35_gate_norm"] = float(_dec_g.disease_base.detach().norm())
                _axn = getattr(path_head, "axes", None)
                for _i in range(int(_eff.shape[0])):
                    _nm = _axn[_i] if (_axn is not None and _i < len(_axn)) else str(_i)
                    details[f"metric/v35_gate_eff_{_nm}"] = float(_eff[_i].detach())     # TRUE self-selection signal
                    details[f"metric/v35_gate_s_{_nm}"] = float(_s[_i].detach())
        except Exception:
            pass

        # v35 METRIC-CONTRAST: amyloid-aligned curvature pressure (donor-balanced, perm-controlled,
        # nuisance-guarded), computed in the model forward (grads DDP-tracked through the decoder
        # score). Here scale + log: contrast_real should grow >0, while contrast_perm / nuis stay ~0
        # (the donor-perm + random-direction controls). lambda 0 / off / too-few-donors => skipped.
        try:
            _mc = getattr(model_out, "metric_contrast_out", None)
            _mch = getattr(self.system, "metric_contrast_head", None)
            if _mc is not None and getattr(_mc, "active", False) and _mch is not None:
                _lam_mc = float(getattr(_mch.cfg, "lambda_metric_contrast", 0.0))
                if _lam_mc > 0.0:
                    total = total + _lam_mc * _mc.loss
                details["loss/v35_metric_contrast"] = float(_mc.loss.detach())
                details["metric/v35_contrast_real"] = float(_mc.contrast_real)
                details["metric/v35_contrast_perm"] = float(_mc.contrast_perm)
                details["metric/v35_contrast_nuis"] = float(_mc.contrast_nuisance)
                details["metric/v35_q_hi"] = float(_mc.q_hi)
                details["metric/v35_q_lo"] = float(_mc.q_lo)
        except Exception:
            pass

        # v26/v27 PATHOLOGY-TIE: cosine-align the decoder's pathology_head (which produces the
        # named-axis projection p that drives the celltype-conditioned + pairwise patterns) to
        # the SUPERVISED pathology_aux head directions. THIS is what makes the new decoder term
        # *pathology-tied* rather than free coords (v24's failure): p_k responds to the same
        # z-direction as the supervised Braak/Thal/LATE/Lewy axis, so its curvature lands in
        # DISEASE directions. Direction-only (magnitude stays free, absorbed by the patterns).
        # lambda 0 or pathology decoder off => skipped. Fully guarded.
        try:
            _lam_pt = path_curriculum_scale * (float(getattr(path_head.cfg, "lambda_pathology_tie", 0.0))
                       if path_head is not None and getattr(path_head, "cfg", None) is not None else 0.0)
            _ph = getattr(getattr(self.system, "decoder", None), "pathology_head", None)
            _hh2 = getattr(path_head, "head", None) if path_head is not None else None
            if _lam_pt > 0.0 and isinstance(_ph, torch.nn.Linear) and isinstance(_hh2, torch.nn.Linear):
                _m2 = min(int(_ph.weight.shape[0]), int(_hh2.weight.shape[0]))
                _cos2 = torch.nn.functional.cosine_similarity(
                    _ph.weight[:_m2], _hh2.weight[:_m2].detach(), dim=1)
                _tie2 = (1.0 - _cos2).mean()
                total = total + _lam_pt * _tie2
                details["metric/v26_pathology_tie"] = float(_tie2.detach())
        except Exception:
            pass

        # v26/v27 pathology DECODER magnitude — how much the celltype-conditioned (+v27 pairwise
        # synergy) term shifts the score. Stashed by the decoder each forward; surfaced here so the
        # progress log shows it growing from zero-init (the new decoder actually "turning on").
        try:
            _dec = getattr(self.system, "decoder", None)
            if _dec is not None and getattr(_dec, "use_pathology_decoder", False):
                _pm = getattr(_dec, "_last_pathdec_absmean", None)
                if _pm is not None:
                    details["metric/v26_pathdec_absmean"] = float(_pm)
                if getattr(_dec, "pathology_pairwise", False):
                    _sm = getattr(_dec, "_last_synergy_absmean", None)
                    if _sm is not None:
                        details["metric/v27_synergy_absmean"] = float(_sm)
        except Exception:
            pass

        # v20: LIVE z effective rank (participation ratio) — the headline "is z
        # spreading across axes?" signal, watched in the log instead of via a post-hoc
        # probe. TWO versions, because the σ-attached SAMPLED z_perp over-states the
        # rank whenever z is collapsed:
        #   - HEADLINE  metric/v20_z_participation  : on the posterior MEAN z_perp
        #     = (mu_q - mu_p)/sigma_p (NO sampling noise) → matches the post-hoc probe
        #     (which uses sample_latent=False). This is the number to trust.
        #   - metric/v20_z_participation_sampled    : on the sampled z_perp. The
        #     reparam noise is ~isotropic, so when the real (mean) z is collapsed it
        #     DOMINATES and inflates this toward d_z (v21: sampled 5.07 vs mean 1.14;
        #     v20b spread: sampled 3.1 ≈ mean 3.15). A large sampled≫mean GAP is itself
        #     a collapse alarm. EMA-smoothed (a 128-cell batch estimate is noisy).
        # Fully guarded: a diagnostics error must NEVER break the training step.
        try:
            def _participation(_z):
                if _z is None or int(_z.shape[0]) < 4:
                    return None
                # detach + float32 + autocast OFF: linalg.eigvalsh is not autocast-
                # eligible and raises under the bf16 autocast region (silently caught
                # before → the metric never appeared in v20's first run).
                with torch.autocast(device_type=_z.device.type, enabled=False):
                    _zc = _z.detach().float()
                    _zc = _zc - _zc.mean(0, keepdim=True)
                    _cov = (_zc.t() @ _zc) / max(int(_zc.shape[0]) - 1, 1)
                    _ev = torch.linalg.eigvalsh(_cov).clamp_min(0.0)
                _s1 = float(_ev.sum())
                _s2 = float((_ev * _ev).sum())
                return (_s1 * _s1) / _s2 if _s2 > 0.0 else None

            # MEAN z_perp from the posterior mean (no reparam noise) — the honest rank.
            _muq = getattr(model_out.encoder_out, "mu_q", None)
            _mup = getattr(model_out, "mu_p", None)
            _sigp = getattr(model_out, "sigma_p", None)
            _zp_mean_raw = (
                (_muq - _mup) / (_sigp + 1e-5)
                if (_muq is not None and _mup is not None and _sigp is not None)
                else None
            )
            # HEADLINE on the CLEAN z (after the same sex/celltype projection the decoder sees), so the
            # reported rank reflects the z the model actually USES. The RAW (pre-projection) rank looks
            # low only because sex/celltype variance dominates it (sex-rank sweep: raw 3.2 -> clean 7.6).
            _zp_mean = _zp_mean_raw
            _nproj = getattr(self.system, "nuisance_projector", None)
            # Snapshot the BATCH (real decoder forward, sampled-z) removed-energy BEFORE the mean-z
            # re-projection just below overwrites projector._diag. rank/eigvals depend only on the EMA
            # buffers (input-independent) ⇒ only removed_energy_frac needs snapshotting here.
            if _nproj is not None:
                try:
                    self._proj_batch_removed = _nproj.diagnostics().get("removed_energy_frac")
                except Exception:
                    self._proj_batch_removed = None
            if (_zp_mean_raw is not None and _nproj is not None
                    and bool(getattr(getattr(_nproj, "cfg", None), "enabled", False))):
                try:
                    with torch.no_grad():
                        _zp_mean = _nproj(_zp_mean_raw, update=False)
                except Exception:
                    _zp_mean = _zp_mean_raw
            _pr_mean = _participation(_zp_mean)
            if _pr_mean is not None:
                _prev = getattr(self, "_v20_pr_ema", None)
                _pr_mean = _pr_mean if _prev is None else 0.9 * float(_prev) + 0.1 * _pr_mean
                self._v20_pr_ema = _pr_mean
                details["metric/v20_z_participation"] = _pr_mean
            _pr_raw = _participation(_zp_mean_raw)   # pre-projection (nuisance-laden) rank, for contrast
            if _pr_raw is not None:
                details["metric/v20_z_participation_raw"] = _pr_raw

            _pr_samp = _participation(getattr(model_out, "z_perp", None))
            if _pr_samp is not None:
                _prevs = getattr(self, "_v20_pr_ema_sampled", None)
                _pr_samp = _pr_samp if _prevs is None else 0.9 * float(_prevs) + 0.1 * _pr_samp
                self._v20_pr_ema_sampled = _pr_samp
                details["metric/v20_z_participation_sampled"] = _pr_samp
        except Exception:
            pass

        # v28 nuisance-PROJECTOR self-diagnostics (multi-rank sex). removed_energy + top sex eigenvalue
        # are the cheap in-train REFILL proxies: if z keeps re-encoding sex, the top sex eigenvalue
        # stays high and removed_energy stays large epoch after epoch. (Definitive refill = the
        # periodic post-hoc z_clean->sex probe on checkpoints.)
        try:
            _proj = getattr(self.system, "nuisance_projector", None)
            if _proj is not None and bool(getattr(getattr(_proj, "cfg", None), "enabled", False)):
                _pd = _proj.diagnostics()
                _re = getattr(self, "_proj_batch_removed", None)   # batch sampled-z value (pre-overwrite)
                if _re is None:
                    _re = _pd.get("removed_energy_frac")           # fallback: whatever _diag currently holds
                if _re is not None:
                    details["metric/proj_removed_energy"] = float(_re)
                if "n_sex_dirs" in _pd:
                    details["metric/proj_sex_rank"] = float(_pd["n_sex_dirs"])
                if "n_ct_dirs" in _pd:
                    details["metric/proj_ct_rank"] = float(_pd["n_ct_dirs"])
                if "n_eff_dirs" in _pd:                     # TRUE rank removed after sex/ct dedup (SVD)
                    details["metric/proj_eff_rank"] = float(_pd["n_eff_dirs"])
                _ev = _pd.get("sex_eigvals_top")
                if _ev:
                    details["metric/proj_sex_eig_top"] = float(_ev[0])
        except Exception:
            pass

        # --------------------------------------------------------------
        # 2. BAM residual uncertainty calibration
        # --------------------------------------------------------------
        tech_id = _resolve_tech_id(batch)
        unc_out = v5_uncertainty_residual_losses(
            raw_uncertainty=model_out.encoder_out.cell_uncertainty,
            rec_per_cell=model_out.decoder_out.nll_per_cell,
            celltype_id=batch["celltype_id"],
            tech_id=tech_id,
            depth_value=batch.get("depth_value", None),
            config=self.v6_uncertainty_config,
        )
        if self.v6_uncertainty_config.enabled:
            total = (
                total
                + float(self.v6_uncertainty_config.lambda_relative_rec)
                * unc_out.relative_rec
                + float(self.v6_uncertainty_config.lambda_rec_alignment)
                * unc_out.rec_alignment
                + float(self.v6_uncertainty_config.lambda_depth_corr)
                * unc_out.depth_corr
                + float(self.v6_uncertainty_config.lambda_saturation)
                * unc_out.saturation
            )

        summary = (
            self.v6_reference_prior.last_summary
            if self.v6_reference_prior is not None
            else None
        )
        details["metric/v6_ref_prior_active_celltypes"] = float(
            0 if summary is None else summary.n_active_celltypes
        )
        details["metric/v6_ref_prior_cells"] = float(
            0 if summary is None else summary.n_total_reference_cells
        )
        details["metric/v6_ref_prior_donors"] = float(
            0 if summary is None else summary.n_total_reference_donors
        )
        details["metric/v6_ref_prior_update_norm"] = float(
            0.0 if summary is None else summary.mean_update_norm
        )
        details["metric/v6_ref_prior_norm_mean"] = float(
            0.0 if summary is None else summary.mean_prior_norm
        )
        details["metric/v6_ref_prior_norm_max"] = float(
            0.0 if summary is None else summary.max_prior_norm
        )
        details["weight/v6_ref_prior_ema_momentum"] = float(
            self.v6_reference_prior_config.ema_momentum
        )
        details["metric/v6_prior_mu_trainable"] = float(
            bool(getattr(model_out.mu_p, "requires_grad", False))
        )

        details["loss/v6_state_abs_floor"] = float(
            state_out.loss_state_abs.detach().cpu()
        )
        details["loss/v6_state_fraction_floor"] = float(
            state_out.loss_state_fraction.detach().cpu()
        )
        details["loss/v6_tech_to_state"] = float(
            state_out.loss_tech_to_state.detach().cpu()
        )
        details["metric/v6_state_abs"] = float(state_out.state_abs.detach().cpu())
        details["metric/v6_base_abs"] = float(state_out.base_abs.detach().cpu())
        details["metric/v6_tech_abs"] = float(state_out.tech_abs.detach().cpu())
        details["metric/v6_state_fraction"] = float(
            state_out.state_fraction.detach().cpu()
        )
        details["metric/v6_tech_to_state"] = float(
            state_out.tech_to_state.detach().cpu()
        )
        # v19 biological-covariate factorization: the sex baseline now shares the
        # decoder score budget. sex_abs/sex_to_state are 0 and state_fraction_no_sex
        # == state_fraction whenever the decoder has no sex term (pre-v19 identical).
        details["metric/v6_sex_abs"] = float(state_out.sex_abs.detach().cpu())
        details["metric/v6_sex_to_state"] = float(
            state_out.sex_to_state.detach().cpu()
        )
        details["metric/v6_precision_abs"] = float(
            state_out.precision_abs.detach().cpu()
        )
        details["metric/v6_state_fraction_no_sex"] = float(
            state_out.state_fraction_no_sex.detach().cpu()
        )
        details["weight/v6_lambda_state_abs"] = float(
            self.v6_state_usage_config.lambda_state_abs
        )
        details["weight/v6_lambda_state_fraction"] = float(
            self.v6_state_usage_config.lambda_state_fraction
        )
        details["weight/v6_lambda_tech_to_state"] = float(
            self.v6_state_usage_config.lambda_tech_to_state
        )

        details["loss/v6_unc_relative_rec"] = float(
            unc_out.relative_rec.detach().cpu()
        )
        details["loss/v6_unc_rec_alignment"] = float(
            unc_out.rec_alignment.detach().cpu()
        )
        details["loss/v6_unc_depth_corr"] = float(unc_out.depth_corr.detach().cpu())
        details["loss/v6_unc_saturation"] = float(unc_out.saturation.detach().cpu())
        # 1=firewall loss is ADDED to total (enabled), 0=computed-but-NOT-applied (diagnostic only).
        # The log uses this to show "firewall <v>" vs "firewall OFF·진단만 <v>" — never a misleading active-looking 0.
        details["metric/v6_firewall_active"] = 1.0 if self.v6_uncertainty_config.enabled else 0.0
        details["unc/v6_rel_std"] = float(unc_out.rel_unc_std.detach().cpu())
        details["unc/v6_raw_std"] = float(unc_out.raw_unc_std.detach().cpu())
        details["weight/v6_bam_rel_w_min"] = float(unc_out.rel_w_min.detach().cpu())
        details["weight/v6_bam_rel_w_mean"] = float(unc_out.rel_w_mean.detach().cpu())
        details["weight/v6_bam_rel_w_max"] = float(unc_out.rel_w_max.detach().cpu())

        # --------------------------------------------------------------
        # v6_min_mixer safety regularizers (only fire when mixer is on).
        # --------------------------------------------------------------
        decoder = self.system.decoder
        mixer_gates = getattr(decoder, "_last_mixer_gates", None)

        # 1) Gate balance: keep the per-branch gate mean near 1.0 so no
        #    single branch (base, tech, state) is silenced. Identity init
        #    sits at 1.0 for every gate; we don't want any drift to be
        #    free.
        if (
            mixer_gates is not None
            and self.v6_lambda_gate_balance > 0.0
        ):
            gate_mean = mixer_gates.mean(dim=0)                    # [3]
            gate_balance = (gate_mean - 1.0).pow(2).sum()
            total = total + float(self.v6_lambda_gate_balance) * gate_balance.to(dtype=total.dtype)
            details["loss/v6_mixer_gate_balance"] = float(gate_balance.detach().cpu())
            details["metric/v6_mixer_gate_base_mean"] = float(gate_mean[0].detach().cpu())
            details["metric/v6_mixer_gate_tech_mean"] = float(gate_mean[1].detach().cpu())
            details["metric/v6_mixer_gate_state_mean"] = float(gate_mean[2].detach().cpu())
            details["weight/v6_lambda_gate_balance"] = float(self.v6_lambda_gate_balance)
        elif mixer_gates is not None:
            # mixer on but regularizer off: log gate values anyway.
            with torch.no_grad():
                gm = mixer_gates.mean(dim=0)
                details["metric/v6_mixer_gate_base_mean"] = float(gm[0].cpu())
                details["metric/v6_mixer_gate_tech_mean"] = float(gm[1].cpu())
                details["metric/v6_mixer_gate_state_mean"] = float(gm[2].cpu())

        # 2) Reference-safe gate: in reference-normal cells the `state`
        #    branch should contribute little (those cells are pathology-
        #    clean). Penalise `gate_state` away from 1.0 on those cells.
        is_reference = batch.get("is_reference", None)
        if (
            mixer_gates is not None
            and self.v6_lambda_gate_ref_safe > 0.0
            and is_reference is not None
        ):
            ref_mask = is_reference.bool().view(-1)
            if bool(ref_mask.any()):
                # Want gate_state to be small (i.e., near 0, not 1) on
                # reference cells. Penalise (gate_state - 0)^2 on refs.
                ref_state_gates = mixer_gates[ref_mask, 2]
                gate_ref_safe = ref_state_gates.pow(2).mean()
                total = total + float(self.v6_lambda_gate_ref_safe) * gate_ref_safe.to(dtype=total.dtype)
                details["loss/v6_mixer_gate_ref_safe"] = float(gate_ref_safe.detach().cpu())
                details["weight/v6_lambda_gate_ref_safe"] = float(self.v6_lambda_gate_ref_safe)

        # 3) Coefficient zero mean: encourage `E[gamma | reference, c] ≈ 0`
        #    so reference cells don't accumulate non-zero Lie generator
        #    activation. Helps keep the state branch a true deviation.
        if (
            self.v6_lambda_coeff_zero_mean > 0.0
            and hasattr(decoder, "coeff_zero_mean_penalty")
            and is_reference is not None
        ):
            ref_mask = is_reference.bool().view(-1)
            if bool(ref_mask.any()):
                ref_indices = ref_mask.nonzero(as_tuple=True)[0]
                if ref_indices.numel() >= 4:
                    z_ref = model_out.z_perp[ref_indices]
                    ct_ref = batch["celltype_id"][ref_indices]
                    coeff_zm = decoder.coeff_zero_mean_penalty(z_ref, ct_ref)
                    total = total + float(self.v6_lambda_coeff_zero_mean) * coeff_zm.to(dtype=total.dtype)
                    details["loss/v6_coeff_zero_mean"] = float(coeff_zm.detach().cpu())
                    details["weight/v6_lambda_coeff_zero_mean"] = float(self.v6_lambda_coeff_zero_mean)

        # learnable-rank group-lasso (2026-07, doc/lie_rank_learnable_design_v3.md): add to the
        # BACKWARD loss like the other regularizers, warmup 0->1 over N epochs. Kept ~1-3% of loss
        # (calibrated) so the composite selection is only negligibly perturbed. Readout (PR /
        # energy-rank / channel shares / mixer gates) logged as metric/* for post-hoc analysis.
        try:
            _decL = getattr(self.system, "decoder", None)
            _lamL = float(getattr(_decL, "lambda_lie_rank_sparse", 0.0)) if _decL is not None else 0.0
            _penL = getattr(_decL, "_lie_rank_pen", None)
            if _lamL > 0.0 and _penL is not None and bool(getattr(self.system, "training", False)):
                _wu = int(getattr(_decL, "lie_rank_warmup_epochs", 0) or 0)
                # current_epoch_index starts at 1 (runner_base: `for epoch in range(1, epochs+1)`),
                # so epoch/warmup gives the intended ramp: ep1->1/3, ep2->2/3, ep3->full (warmup=3).
                _warm = 1.0 if _wu <= 0 else min(1.0, float(self.current_epoch_index) / float(_wu))
                # SELF-SCALE: interpret lambda_lie_rank_sparse as a TARGET FRACTION of the loss.
                # penalty value = warm * frac * loss; its gradient is a positively-scaled group-lasso
                # gradient. Removes fragile absolute-lambda calibration (pen magnitude drifts in
                # training); keeps pruning pressure a stable % of the objective. reviewer-suggested.
                # NOTE: this add is TRAIN-only (guard above); val loss stays penalty-free => the
                # composite/early-stop selection (monitored on val) is not contaminated. loss_total_task
                # below is the penalty-free value for a clean vs-fixed-rank comparison on the train side.
                # SAFEGUARD (self-scale gradient = frac*(total/pen)*grad(pen) blows up as pen->0):
                # floor pen at 10% of its EMA so `scale` can't exceed ~10x its running level. Bounds
                # the penalty gradient without a magic constant; if pen is pathologically small we
                # simply under-penalize (safe direction). Watch metric/lie_grad_scale for spikes.
                _pv = float(_penL.detach())
                _ema = getattr(self, "_lie_pen_ema", None)
                _ema = _pv if (_ema is None or _ema <= 0.0) else (0.98 * _ema + 0.02 * _pv)
                self._lie_pen_ema = _ema
                _pfloor = max(_pv, 0.1 * _ema)
                _scale = (total.detach() / (_pfloor + 1e-12)).to(dtype=total.dtype)
                details["metric/loss_total_task"] = float(total.detach().cpu())   # penalty-FREE
                details["metric/lie_grad_scale"] = float(_scale)                  # watch: blow-up guard
                total = total + (_warm * _lamL) * _scale * _penL.to(dtype=total.dtype)
                details["metric/lie_penalty"] = float(_penL.detach())
                details["metric/lie_penalty_frac"] = float(_warm * _lamL)      # exact loss-fraction
            _c = getattr(_decL, "_last_rank_c", None)                              # [M,R] logged even at lambda=0 (calib)
            if _c is not None:
                _s1 = _c.sum(dim=1); _s2 = _c.pow(2).sum(dim=1)
                _pr = (_s1 * _s1) / (_s2 + 1e-12)                                  # [M] participation ratio
                details["metric/lie_rank_PR_mean"] = float(_pr.mean())
                _cs = _c.sort(dim=1, descending=True).values
                _cum = _cs.cumsum(dim=1) / (_cs.sum(dim=1, keepdim=True) + 1e-12)
                details["metric/lie_rank_k95_mean"] = float(((_cum < 0.95).float().sum(dim=1) + 1.0).mean())
                details["metric/lie_rank_sumc"] = float(_c.sum())
            _cr = getattr(_decL, "_chan_rms", None)                               # channel leak check
            if isinstance(_cr, dict):
                _tot = float(sum(_cr.values())) + 1e-8
                for _k, _v in _cr.items():
                    details[f"metric/lie_chanfrac_{_k}"] = float(_v) / _tot
            _mg = getattr(_decL, "_last_mixer_gates", None)
            if _mg is not None:
                _mgm = _mg.detach().float().mean(dim=0)
                if _mgm.numel() >= 3:
                    details["metric/lie_gate_base"] = float(_mgm[0])
                    details["metric/lie_gate_tech"] = float(_mgm[1])
                    details["metric/lie_gate_state"] = float(_mgm[2])
        except Exception:
            pass

        # Pathology-correction rank is learned from the full algebraic basis.
        # The normalized expected cardinality has no target and no minimum;
        # reconstruction/pathology gradients decide which components justify
        # their cost.  Validation remains penalty-free.
        pathology_rank_gate = getattr(decoder, "pathology_rank_gate", None)
        if pathology_rank_gate is not None:
            rank_diag = pathology_rank_gate.diagnostics()
            details["metric/pathology_rank_capacity"] = float(
                rank_diag["capacity"]
            )
            details["metric/pathology_rank_expected"] = float(
                rank_diag["expected_rank"]
            )
            details["metric/pathology_rank_hard"] = float(
                rank_diag["hard_rank"]
            )
            details["metric/pathology_rank_search_expected"] = float(
                rank_diag["search_expected_rank"]
            )
            details["metric/pathology_rank_search_hard"] = float(
                rank_diag["search_hard_rank"]
            )
            details["metric/pathology_rank_uncertain_fraction"] = float(
                rank_diag["uncertain_fraction"]
            )
            details["metric/pathology_rank_temperature"] = float(
                rank_diag["temperature"]
            )
            details["metric/pathology_rank_finalized"] = float(
                bool(rank_diag["finalized"])
            )
            details["metric/pathology_rank_freeze_ready"] = float(
                bool(rank_diag["freeze_ready"])
            )
            details["metric/pathology_rank_freeze_count_spread"] = float(
                rank_diag["freeze_count_spread"]
            )
            details[
                "metric/pathology_rank_freeze_near_threshold_uncertain_fraction"
            ] = float(
                rank_diag["freeze_near_threshold_uncertain_fraction"]
            )
            details["metric/pathology_rank_freeze_count_low"] = float(
                rank_diag["freeze_count_at_low_threshold"]
            )
            details["metric/pathology_rank_freeze_count_mid"] = float(
                rank_diag["freeze_count_at_hard_threshold"]
            )
            details["metric/pathology_rank_freeze_count_high"] = float(
                rank_diag["freeze_count_at_high_threshold"]
            )
            rank_mode = str(rank_diag["mode"])
            for mode_name in (
                "fixed_rank_warmup",
                "all_open",
                "soft_rank_learning",
                "hard_straight_through",
                "frozen_hard",
            ):
                details[f"metric/pathology_rank_mode_{mode_name}"] = float(
                    rank_mode == mode_name
                )
            sparsity_multiplier = float(
                pathology_rank_gate.sparsity_multiplier()
            )
            details["weight/pathology_rank_sparsity_multiplier"] = (
                sparsity_multiplier
            )
            if bool(getattr(self.system, "training", False)) and sparsity_multiplier > 0.0:
                penalty = pathology_rank_gate.cardinality_penalty()
                penalty_value = float(penalty.detach().cpu())
                fraction = float(
                    pathology_rank_gate.config.cardinality_loss_fraction
                )
                scale = total.detach() / penalty.detach().clamp_min(1.0e-6)
                total = total + (
                    sparsity_multiplier
                    * fraction
                    * scale.to(dtype=total.dtype)
                    * penalty.to(dtype=total.dtype)
                )
                details["loss/pathology_rank_cardinality"] = penalty_value
                details["weight/pathology_rank_cardinality_fraction"] = (
                    sparsity_multiplier * fraction
                )

        details["loss/total"] = float(total.detach().cpu())

        is_train = bool(getattr(self.system, "training", False))
        self._v6_loss_print_counter += 1
        should_print = (
            is_train
            and self.v6_diag_print_every > 0
            and self._v6_loss_print_counter % self.v6_diag_print_every == 0
            and self._is_rank0()
        )
        if should_print:
            # v6_min_mixer: surface gate means + coeff_zero_mean penalty
            # so unattended runs can spot mixer drift without tailing the
            # full details dict. Defaults to "-" when the mixer is off.
            gate_b = details.get("metric/v6_mixer_gate_base_mean", None)
            gate_t = details.get("metric/v6_mixer_gate_tech_mean", None)
            gate_s = details.get("metric/v6_mixer_gate_state_mean", None)
            czm = details.get("loss/v6_coeff_zero_mean", None)
            if gate_b is not None and gate_t is not None and gate_s is not None:
                mixer_line = (
                    f"\n  Mixer gates     : base={float(gate_b):.3f}  "
                    f"tech={float(gate_t):.3f}  "
                    f"state={float(gate_s):.3f}"
                )
                if czm is not None:
                    mixer_line += f"  coeff_zm={float(czm):.4f}"
            else:
                mixer_line = ""
            print(
                f"[diagnostics | reference/state | epoch {self.current_epoch_index}]\n"
                f"  Reference prior : active_celltypes={details['metric/v6_ref_prior_active_celltypes']:.0f}  "
                f"update_norm={details['metric/v6_ref_prior_update_norm']:.6f}  "
                f"prior_norm_mean={details['metric/v6_ref_prior_norm_mean']:.3f}  "
                f"mu_trainable={details['metric/v6_prior_mu_trainable']:.0f}\n"
                f"  State usage     : state_abs={float(state_out.state_abs.detach().cpu()):.3f}  "
                f"state_fraction={float(state_out.state_fraction.detach().cpu()):.3f}  "
                f"tech_to_state={float(state_out.tech_to_state.detach().cpu()):.3f}\n"
                f"  BAM-attn diag   : relative_std={float(unc_out.rel_unc_std.detach().cpu()):.3f}  "
                f"raw_std={float(unc_out.raw_unc_std.detach().cpu()):.6f}  "
                f"depth_corr={float(unc_out.depth_corr.detach().cpu()):.6f}  "
                f"relative_weight=[{float(unc_out.rel_w_min.detach().cpu()):.3f}, "
                f"{float(unc_out.rel_w_mean.detach().cpu()):.3f}, "
                f"{float(unc_out.rel_w_max.detach().cpu()):.3f}]"
                + mixer_line,
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
