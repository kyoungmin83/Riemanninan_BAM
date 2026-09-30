from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


try:
    from kmlee_bam.model.celltype_prior import CellTypePrior
except ImportError:
    try:
        from kmlee_bam.model.celltype_prior import CellTypePrior
    except ImportError:
        from celltype_prior import CellTypePrior

try:
    from kmlee_bam.objectives.celltype_alignment import (
        compute_celltype_alignment_loss,
        EMAAlignmentState,
    )
except ImportError:
    from celltype_alignment import compute_celltype_alignment_loss, EMAAlignmentState

# Decoder is intentionally typed as Any: this loss supports both
# legacy additive decoders and the current LieActionOrdinalDecoder.


# ======================================================================
# Output container
# ======================================================================
@dataclass
class LossOutput:
    """
    Rich container for scalar losses and per-cell diagnostics.

    Attributes
    ----------
    total : scalar tensor
        Final scalar used for backward().
    main : scalar tensor
        Weighted reconstruction + beta_state * state KL.
    rec : scalar tensor
        Mean raw reconstruction NLL across cells.
    rec_weighted : scalar tensor
        Mean BAM-weighted reconstruction NLL across cells.
    kl_state : scalar tensor
        Mean state KL across cells.
    bam_kl : scalar tensor
        BAM attention KL term.
    reg_tech : scalar tensor
        Technical baseline regularizer.
    reg_gauge : scalar tensor
        Donor-balanced gauge regularizer.
    cls : scalar tensor
        Optional auxiliary classification loss.
    weights_per_cell : [B]
        BAM-derived reconstruction weights.
    raw_uncertainty : [B] or None
        Original cell-level BAM uncertainty.
    clipped_uncertainty : [B] or None
        Clipped uncertainty actually used for weighting.
    rec_per_cell : [B]
        Raw reconstruction NLL per cell.
    rec_weighted_per_cell : [B]
        Weighted reconstruction term per cell.
    kl_state_per_cell : [B]
        State KL per cell.
    details : dict[str, float]
        Flat scalar diagnostics for logging.
    """

    total: torch.Tensor
    main: torch.Tensor
    rec: torch.Tensor
    rec_weighted: torch.Tensor
    kl_state: torch.Tensor
    bam_kl: torch.Tensor
    reg_tech: torch.Tensor
    reg_gauge: torch.Tensor
    cls: torch.Tensor
    ref_center: torch.Tensor
    ref_state: torch.Tensor
    weights_per_cell: torch.Tensor
    raw_uncertainty: Optional[torch.Tensor]
    clipped_uncertainty: Optional[torch.Tensor]
    rec_per_cell: torch.Tensor
    rec_weighted_per_cell: torch.Tensor
    kl_state_per_cell: torch.Tensor
    tech_score_mag: torch.Tensor
    # Celltype-alignment + whitening on z_perp. Optional with a default so
    # existing LossOutput(...) constructions in subclass trainers keep working;
    # TotalLoss always populates them (graph-connected zeros when disabled), and
    # the propagating trainers forward them explicitly via align=out.align.
    align: Optional[torch.Tensor] = None
    white: Optional[torch.Tensor] = None
    align_sex: Optional[torch.Tensor] = None
    details: Dict[str, float] = field(default_factory=dict)


# ======================================================================
# Helper functions
# ======================================================================
def zero_scalar_like(x: torch.Tensor) -> torch.Tensor:
    return torch.zeros((), dtype=x.dtype, device=x.device)


@torch.no_grad()
def empirical_group_weights(group_id: torch.Tensor, n_groups: int) -> torch.Tensor:
    """
    Estimate empirical group weights from a minibatch.

    Parameters
    ----------
    group_id : [B] long
    n_groups : int

    Returns
    -------
    weights : [n_groups]
        Non-negative weights summing to 1.
    """
    if group_id.ndim != 1:
        raise ValueError(f"group_id must have shape [B], got {tuple(group_id.shape)}.")
    if group_id.dtype != torch.long:
        raise TypeError(f"group_id must be torch.long, got {group_id.dtype}.")
    if n_groups <= 0:
        raise ValueError("n_groups must be positive.")

    counts = torch.bincount(group_id, minlength=n_groups).to(dtype=torch.float32)
    return counts / counts.sum().clamp_min(1.0)



def build_bam_reconstruction_weights(
    cell_uncertainty: Optional[torch.Tensor],   # [B], raw entropy
    *,
    n_valid_tokens: int,                        # effective T (e.g. number of genes, excluding CLS)
    tau: float = 2.0,
    detach_uncertainty: bool = True,
    normalize_weights: bool = True,
    w_min: Optional[float] = None,
    w_max: Optional[float] = None,
    eps: float = 1e-8,
    u_min: Optional[float] = None,
    u_max: Optional[float] = None,
) -> tuple[Optional[torch.Tensor], Optional[torch.Tensor]]:
    """
    Absolute-normalized BAM reconstruction weights.

    u_i = r_i / log(T)
    w_i = exp(-tau * u_i)
    optional clamp
    optional mean normalization

    Returns
    -------
    weights : [B] or None
    normalized_uncertainty : [B] or None
    """
    if cell_uncertainty is None:
        return None, None

    if cell_uncertainty.ndim != 1:
        raise ValueError(
            f"cell_uncertainty must have shape [B], got {tuple(cell_uncertainty.shape)}."
        )
    if n_valid_tokens <= 1:
        raise ValueError("n_valid_tokens must be > 1.")
    if tau < 0:
        raise ValueError("tau must be non-negative.")

    r = cell_uncertainty.detach() if detach_uncertainty else cell_uncertainty

    entropy_scale = torch.log(
        torch.tensor(float(n_valid_tokens), device=r.device, dtype=r.dtype)
    ).clamp_min(eps)

    u = r / entropy_scale
    
    u_used =u
    if u_min is not None or u_max is not None:
    	lo = -float("inf") if u_min is None else float(u_min)
    	hi = float("inf") if u_max is None else float(u_max)
    	u_used = u.clamp(min=lo, max=hi)
    	

    weights = torch.exp(-tau * u_used)

    if w_min is not None or w_max is not None:
        lo = 0.0 if w_min is None else float(w_min)
        hi = float("inf") if w_max is None else float(w_max)
        weights = weights.clamp(min=lo, max=hi)

    if normalize_weights:
        weights = weights / weights.mean().clamp_min(eps)

    return weights, u_used

def reference_center_loss(
    z_perp: Optional[torch.Tensor],
    celltype_id: torch.Tensor,
    is_reference: Optional[torch.Tensor],
    *,
    min_cells_per_class: int = 2,
) -> torch.Tensor:
    """
    Penalise non-zero reference mean within each cell type.

        L = mean_c || mean_{i: ref_i=1, c_i=c} z_perp_i ||^2

    This anchors each Subclass-specific latent universe so that pathology-clean
    reference cells sit near the local origin.
    """
    if z_perp is None:
        return zero_scalar_like(celltype_id.float())

    if is_reference is None:
        return zero_scalar_like(z_perp)

    if z_perp.ndim != 2:
        raise ValueError(f"z_perp must have shape [B, d_z], got {tuple(z_perp.shape)}.")
    if celltype_id.ndim != 1:
        raise ValueError(f"celltype_id must have shape [B], got {tuple(celltype_id.shape)}.")

    is_reference = is_reference.bool()
    if is_reference.ndim != 1:
        raise ValueError(f"is_reference must have shape [B], got {tuple(is_reference.shape)}.")
    if is_reference.shape[0] != z_perp.shape[0]:
        raise ValueError("Batch size mismatch between z_perp and is_reference.")

    if int(is_reference.sum().item()) == 0:
        return zero_scalar_like(z_perp)

    losses = []
    classes = torch.unique(celltype_id[is_reference])

    for c in classes:
        m = is_reference & (celltype_id == c)
        if int(m.sum().item()) < int(min_cells_per_class):
            continue
        mean_c = z_perp[m].mean(dim=0)
        losses.append(mean_c.pow(2).mean())

    if len(losses) == 0:
        return zero_scalar_like(z_perp)

    return torch.stack(losses).mean()


def reference_state_penalty(
    state_score: Optional[torch.Tensor],
    is_reference: Optional[torch.Tensor],
) -> torch.Tensor:
    """
    Optional: discourage the state branch from doing too much on reference cells.

    Keep lambda_ref_state=0.0 in the first experiment.
    """
    if state_score is None:
        raise ValueError("state_score is None; decoder_out must expose state_score.")

    if is_reference is None:
        return zero_scalar_like(state_score)

    is_reference = is_reference.bool()
    if int(is_reference.sum().item()) == 0:
        return zero_scalar_like(state_score)

    return state_score[is_reference].abs().mean()

# ======================================================================
# Main loss module
# ======================================================================
class TotalLoss(nn.Module):
    """
    Hybrid project-ready loss module.

    Philosophy
    ----------
    1) Keep the project-specific mathematical structure:

           main = mean_i[w_i * rec_i] + beta_state * mean_i[KL_i]
           total = main
                 + lambda_tech  * R_tech
                 + lambda_gauge * R_gauge
                 + lambda_bam   * L_BAM_KL
                 + lambda_cls   * L_cls

       BAM uncertainty weights only the reconstruction term.
       It does NOT weaken the state KL.

    2) Keep the cleaner class-style interface and rich diagnostics.

    Expected inputs
    ---------------
    - prior       : CellTypePrior
    - decoder     : decoder object exposing n_genes, n_tech, tech_regularizer, gauge_regularizer_db_c
    - encoder_out : object with attributes
          mu_q, logvar_q, attn_kl, cell_uncertainty
    - decoder_out : decoder output with nll_per_cell and state_score
    """

    def __init__(
        self,
        *,
        beta_state: float = 1.0,
        lambda_tech: float = 0.0,
        lambda_tech_mag: float = 0.0,
        lambda_gauge: float = 0.0,
        lambda_bam: float = 1.0,
        lambda_cls: float = 0.0,
        lambda_ref_center: float = 0.0,
        lambda_ref_state: float = 0.0,
        lambda_align: float = 0.0,
        lambda_white: float = 0.0,
        lambda_align_sex: float = 0.0,
        align_sex_all_cells: bool = False,
        align_use_covariance: bool = True,
        align_reference_only: bool = True,
        align_whiten_all_cells: bool = True,
        align_cov_use_ema: bool = True,
        align_ema_decay: float = 0.99,
        ref_center_min_cells: int = 2,
        tau: float = 1.0,
        r_min: float = 0.0,
        r_max: float = 5.0,
        detach_bam_uncertainty: bool = True,
        normalize_bam_weights: bool = True,
        normalize_gene_regularizers: bool = True,
    ) -> None:
        super().__init__()

        if beta_state < 0:
            raise ValueError(f"beta_state must be non-negative, got {beta_state}.")
        if lambda_tech < 0:
            raise ValueError(f"lambda_tech must be non-negative, got {lambda_tech}.")
        if lambda_tech_mag < 0:
            raise ValueError(f"lambda_tech_mag must be non-negative, got {lambda_tech_mag}.")
        if lambda_gauge < 0:
            raise ValueError(f"lambda_gauge must be non-negative, got {lambda_gauge}.")
        if lambda_bam < 0:
            raise ValueError(f"lambda_bam must be non-negative, got {lambda_bam}.")
        if lambda_cls < 0:
            raise ValueError(f"lambda_cls must be non-negative, got {lambda_cls}.")
        if lambda_ref_center < 0:
            raise ValueError(f"lambda_ref_center must be non-negative, got {lambda_ref_center}.")
        if lambda_ref_state < 0:
            raise ValueError(f"lambda_ref_state must be non-negative, got {lambda_ref_state}.")
        if lambda_align < 0:
            raise ValueError(f"lambda_align must be non-negative, got {lambda_align}.")
        if lambda_white < 0:
            raise ValueError(f"lambda_white must be non-negative, got {lambda_white}.")
        if lambda_align_sex < 0:
            raise ValueError(f"lambda_align_sex must be non-negative, got {lambda_align_sex}.")
        if not (0.0 <= align_ema_decay < 1.0):
            raise ValueError(f"align_ema_decay must be in [0, 1), got {align_ema_decay}.")
        if ref_center_min_cells < 1:
            raise ValueError(f"ref_center_min_cells must be >= 1, got {ref_center_min_cells}.")
        if tau < 0:
            raise ValueError(f"tau must be non-negative, got {tau}.")
        if r_max < r_min:
            raise ValueError(f"r_max must be >= r_min, got r_min={r_min}, r_max={r_max}.")

        self.beta_state = beta_state
        self.lambda_tech = lambda_tech
        self.lambda_tech_mag = lambda_tech_mag
        self.lambda_gauge = lambda_gauge
        self.lambda_bam = lambda_bam
        self.lambda_cls = lambda_cls
        self.lambda_ref_center = lambda_ref_center
        self.lambda_ref_state = lambda_ref_state
        # Celltype-ALIGNMENT + WHITENING on z_perp (DDP-trivial replacement for
        # the GRL aux-adversary). Both lambdas==0 ⇒ the loss call is skipped
        # entirely below, so behaviour is byte-identical when disabled.
        self.lambda_align = lambda_align
        self.lambda_white = lambda_white
        # Optional per-SEX alignment (same mechanism as celltype). DEFAULT OFF.
        # The DLPFC+MTG training batch does NOT currently carry a sex label, so
        # enabling this without plumbing a per-cell sex id into the batch raises
        # a clear error from compute_celltype_alignment_loss (see forward()).
        self.lambda_align_sex = lambda_align_sex
        # Align the sex-mean over ALL cells (not reference-only). Reference cells
        # are ≈ single-gender, so reference-only sex alignment is a near no-op.
        self.align_sex_all_cells = bool(align_sex_all_cells)
        self.align_use_covariance = bool(align_use_covariance)
        self.align_reference_only = bool(align_reference_only)
        self.align_whiten_all_cells = bool(align_whiten_all_cells)
        # EMA-accumulated covariance alignment: activates the per-celltype
        # covariance term despite sparse per-batch reference cells by carrying
        # detached running statistics across steps. When False the loss uses the
        # original per-batch path (byte-identical to the pre-EMA behaviour).
        self.align_cov_use_ema = bool(align_cov_use_ema)
        self.align_ema_decay = float(align_ema_decay)
        # The EMA accumulators are registered submodules so their (non-persistent)
        # buffers move with .to(device)/.cuda() and there is exactly one running
        # state per rank. They are lazily sized on first use (d_z / n_celltypes
        # inferred from the data), so TotalLoss needs neither at construction.
        # Only created when the EMA path can actually be used, so the disabled
        # path adds NO new module state.
        if self.align_cov_use_ema and (self.lambda_align != 0.0 or self.lambda_align_sex != 0.0):
            self.ema_align_state = EMAAlignmentState(decay=self.align_ema_decay)
            self.ema_sex_state = EMAAlignmentState(decay=self.align_ema_decay)
        else:
            self.ema_align_state = None
            self.ema_sex_state = None
        self.ref_center_min_cells = ref_center_min_cells
        self.tau = tau
        self.r_min = r_min
        self.r_max = r_max
        self.detach_bam_uncertainty = detach_bam_uncertainty
        self.normalize_bam_weights = normalize_bam_weights
        self.normalize_gene_regularizers = normalize_gene_regularizers

    def forward(
        self,
        *,
        prior: CellTypePrior,
        decoder: Any,
        encoder_out: Any,
        decoder_out: Any,
        celltype_id: torch.Tensor,
        tech_id: torch.Tensor,
        cls_logits: Optional[torch.Tensor] = None,
        cls_targets: Optional[torch.Tensor] = None,
        tech_group_weights: Optional[torch.Tensor] = None,
        n_valid_tokens_for_uncertainty: Optional[int] = None,
        z_perp: Optional[torch.Tensor] = None,
        is_reference: Optional[torch.Tensor] = None,
        reconstruction_weight_per_gene: Optional[torch.Tensor] = None,
        sex_id: Optional[torch.Tensor] = None,
        whitening_target: Optional[torch.Tensor] = None,
    ) -> LossOutput:
        """
        Compute total training loss and diagnostics.
        """
        if decoder_out.nll_per_cell is None:
            raise ValueError(
                "decoder_out.nll_per_cell is None. Call the decoder with y_ord before the loss."
            )
        if celltype_id.dtype != torch.long:
            raise TypeError(f"celltype_id must be torch.long, got {celltype_id.dtype}.")
        if tech_id.dtype != torch.long:
            raise TypeError(f"tech_id must be torch.long, got {tech_id.dtype}.")
        if celltype_id.ndim != 1 or tech_id.ndim != 1:
            raise ValueError("celltype_id and tech_id must both have shape [B].")
        if decoder_out.nll_per_cell.ndim != 1:
            raise ValueError(
                f"decoder_out.nll_per_cell must have shape [B], got {tuple(decoder_out.nll_per_cell.shape)}."
            )
        B = decoder_out.nll_per_cell.shape[0]
        if celltype_id.shape[0] != B or tech_id.shape[0] != B:
            raise ValueError("Batch size mismatch among ids and decoder output.")

        for attr in ("mu_q", "logvar_q", "attn_kl", "cell_uncertainty"):
            if not hasattr(encoder_out, attr):
                raise AttributeError(f"encoder_out must have attribute '{attr}'.")

        # Reconstruction NLL per cell. Default: uniform mean over genes
        # (== decoder_out.nll_per_cell). Stage D: when a per-gene reliability
        # weight is supplied (suspicious observed zeros downweighted), use the
        # weighted gene-mean instead. With all-ones weights this is identical to
        # nll_per_cell, so the path is byte-identical when the weight is absent.
        if reconstruction_weight_per_gene is not None:
            if decoder_out.nll_per_gene is None:
                raise ValueError(
                    "reconstruction_weight_per_gene was provided but decoder_out.nll_per_gene "
                    "is None; per-gene reconstruction weighting needs the per-gene NLL."
                )
            nll_g = decoder_out.nll_per_gene
            w_g = reconstruction_weight_per_gene.to(dtype=nll_g.dtype, device=nll_g.device)
            if w_g.shape != nll_g.shape:
                raise ValueError(
                    f"reconstruction_weight_per_gene shape {tuple(w_g.shape)} != nll_per_gene "
                    f"shape {tuple(nll_g.shape)}."
                )
            rec_per_cell = (nll_g * w_g).sum(dim=1) / w_g.sum(dim=1).clamp_min(1e-8)
        else:
            rec_per_cell = decoder_out.nll_per_cell  # [B]
        raw_uncertainty = encoder_out.cell_uncertainty

        # --------------------------------------------------------------
        # 1) BAM-derived reconstruction weighting
        # --------------------------------------------------------------
      
        if n_valid_tokens_for_uncertainty is None:
            # Backward-compatible fallback. For the module-token model,
            # Trainer should pass system.n_effective_encoder_tokens_for_uncertainty
            # so that entropy is normalised by log(number_of_module_tokens),
            # not by log(number_of_genes).
            n_valid_tokens = int(decoder.n_genes)
        else:
            n_valid_tokens = int(n_valid_tokens_for_uncertainty)

        weights_per_cell, used_uncertainty = build_bam_reconstruction_weights(
            encoder_out.cell_uncertainty,
            n_valid_tokens=n_valid_tokens,
            tau=self.tau,
            detach_uncertainty=self.detach_bam_uncertainty,
            normalize_weights=self.normalize_bam_weights,
            w_min=None,
            w_max=None,
            eps=1e-8,
            u_min=self.r_min,
            u_max=self.r_max,
        )
        
        raw_uncertainty = encoder_out.cell_uncertainty
        clipped_uncertainty = used_uncertainty

        if weights_per_cell is None:
            weights_per_cell = torch.ones_like(rec_per_cell)
            rec_weighted_per_cell = rec_per_cell
        else:
            rec_weighted_per_cell = weights_per_cell * rec_per_cell

        rec = rec_per_cell.mean()
        rec_weighted = rec_weighted_per_cell.mean()

        # --------------------------------------------------------------
        # 2) State KL against cell-type prior
        # --------------------------------------------------------------
        mu_p, logvar_p, _ = prior(celltype_id)
        kl_state_per_cell = prior.kl_divergence(
            mu_q=encoder_out.mu_q,
            logvar_q=encoder_out.logvar_q,
            mu_p=mu_p,
            logvar_p=logvar_p,
        )
        kl_state = kl_state_per_cell.mean()

        # --------------------------------------------------------------
        # 3) Decoder regularizers
        # --------------------------------------------------------------
        # reg_tech is a *centering* penalty: pulls the weighted-mean tech
        # baseline toward 0. It does NOT constrain per-cell tech_score
        # magnitude, so tech can still absorb a large class-prior offset by
        # spreading nonzero baselines that cancel on average. The magnitude
        # penalty below is the complement.
        reg_tech = decoder.tech_regularizer(tech_group_weights)
        reg_gauge = decoder.gauge_regularizer_db_c(
            state_score=decoder_out.state_score,
            celltype_id=celltype_id,
            tech_id=tech_id,
        )
        if self.normalize_gene_regularizers:
            gene_scale = float(decoder.n_genes)
            reg_tech = reg_tech / gene_scale
            reg_gauge = reg_gauge / gene_scale

        # --------------------------------------------------------------
        # 3b) Tech-score magnitude penalty
        # --------------------------------------------------------------
        # See doc/tech_score_magnitude_penalty.md.
        # Tech is allowed to absorb a small per-batch offset, but it must
        # not become a free zero-class shortcut. We squeeze the per-cell
        # tech_score toward 0, leaving room for state and base to carry
        # the ordinal score.
        if hasattr(decoder_out, "tech_score") and decoder_out.tech_score is not None:
            tech_score_mag = decoder_out.tech_score.pow(2).mean()
        else:
            tech_score_mag = zero_scalar_like(rec)

        # --------------------------------------------------------------
        # 4) BAM KL
        # --------------------------------------------------------------
        bam_kl = encoder_out.attn_kl
        if bam_kl.ndim != 0:
            bam_kl = bam_kl.mean()

        # --------------------------------------------------------------
        # 5) Optional classifier
        # --------------------------------------------------------------
        cls = zero_scalar_like(rec)
        if self.lambda_cls != 0.0 or cls_logits is not None or cls_targets is not None:
            if cls_logits is None or cls_targets is None:
                raise ValueError(
                    "cls_logits and cls_targets must both be provided when using auxiliary classification."
                )
            if cls_targets.dtype != torch.long:
                raise TypeError(f"cls_targets must be torch.long, got {cls_targets.dtype}.")
            cls = F.cross_entropy(cls_logits, cls_targets)

        
        # --------------------------------------------------------------
        # 6) Reference-origin regularizers
        # --------------------------------------------------------------
        ref_center = reference_center_loss(
            z_perp=z_perp,
            celltype_id=celltype_id,
            is_reference=is_reference,
            min_cells_per_class=self.ref_center_min_cells,
        )

        ref_state = reference_state_penalty(
            state_score=getattr(decoder_out, "state_score", None),
            is_reference=is_reference,
        )

        # --------------------------------------------------------------
        # 6b) Celltype-ALIGNMENT + WHITENING on z_perp
        # --------------------------------------------------------------
        # Pure functions of the main forward's z_perp (no torch.distributed, no
        # aux batch, no separate backward): the existing main backward + DDP
        # reducer average their gradients exactly like the reconstruction loss.
        # When BOTH lambdas are 0 we skip the call ENTIRELY (truly byte-identical
        # to the pre-alignment path) and emit graph-connected zeros, mirroring
        # how the function itself handles its edge cases.
        align_stats: Dict[str, float] = {
            "n_reference": 0.0,
            "n_celltypes_in_ref": 0.0,
            "l_mean": 0.0,
            "l_cov": 0.0,
            "l_white": 0.0,
            "l_mean_sex": 0.0,
            "l_cov_sex": 0.0,
            "n_sex_in_ref": 0.0,
            "white_active_error": 0.0,
            "white_removed_leakage": 0.0,
            "white_cross_leakage": 0.0,
            "white_target_rank": float(z_perp.shape[1]) if z_perp is not None else 0.0,
            "white_removed_rank": 0.0,
            "white_identity_target_loss": 0.0,
            "white_target_loss_delta": 0.0,
        }
        use_sex = self.lambda_align_sex != 0.0
        any_align = (
            self.lambda_align != 0.0 or self.lambda_white != 0.0 or use_sex
        )
        if any_align and z_perp is not None:
            align, white, align_sex, align_stats = compute_celltype_alignment_loss(
                z_perp,
                celltype_id,
                is_reference,
                use_covariance=self.align_use_covariance,
                reference_only=self.align_reference_only,
                whiten_all_cells=self.align_whiten_all_cells,
                ema_state=self.ema_align_state if self.align_cov_use_ema else None,
                ema_decay=self.align_ema_decay,
                update_ema=self.training,
                sex_id=sex_id,
                sex_ema_state=self.ema_sex_state if self.align_cov_use_ema else None,
                use_sex=use_sex,
                sex_all_cells=self.align_sex_all_cells,
                whitening_target=whitening_target,
            )
        elif z_perp is not None:
            align = z_perp.sum() * 0.0
            white = z_perp.sum() * 0.0
            align_sex = z_perp.sum() * 0.0
        else:
            align = zero_scalar_like(rec)
            white = zero_scalar_like(rec)
            align_sex = zero_scalar_like(rec)

        # --------------------------------------------------------------
        # 7) Final assembly
        # --------------------------------------------------------------
        main = rec_weighted + self.beta_state * kl_state
        total = (
            main
            + self.lambda_tech * reg_tech
            + self.lambda_tech_mag * tech_score_mag
            + self.lambda_gauge * reg_gauge
            + self.lambda_bam * bam_kl
            + self.lambda_cls * cls
            + self.lambda_ref_center * ref_center
            + self.lambda_ref_state * ref_state
            + self.lambda_align * align
            + self.lambda_white * white
            + self.lambda_align_sex * align_sex
        )

        details = {
            "loss/total": float(total.detach().cpu()),
            "loss/main": float(main.detach().cpu()),
            "loss/rec": float(rec.detach().cpu()),
            "loss/rec_weighted": float(rec_weighted.detach().cpu()),
            "loss/kl_state": float(kl_state.detach().cpu()),
            "loss/bam_kl": float(bam_kl.detach().cpu()),
            "loss/reg_tech": float(reg_tech.detach().cpu()),
            "loss/tech_score_mag": float(tech_score_mag.detach().cpu()),
            "weight/lambda_tech_mag": float(self.lambda_tech_mag),
            "loss/reg_gauge": float(reg_gauge.detach().cpu()),
            "loss/cls": float(cls.detach().cpu()),
            "loss/ref_center": float(ref_center.detach().cpu()),
            "loss/ref_state": float(ref_state.detach().cpu()),
            "weight/lambda_ref_center": float(self.lambda_ref_center),
            "weight/lambda_ref_state": float(self.lambda_ref_state),
            "loss/align": float(align.detach().cpu()),
            "loss/white": float(white.detach().cpu()),
            "loss/align_sex": float(align_sex.detach().cpu()),
            "loss/align_l_mean": float(align_stats["l_mean"]),
            "loss/align_l_cov": float(align_stats["l_cov"]),
            "loss/align_l_white": float(align_stats["l_white"]),
            "loss/align_l_mean_sex": float(align_stats["l_mean_sex"]),
            "loss/align_l_cov_sex": float(align_stats["l_cov_sex"]),
            "metric/align_n_reference": float(align_stats["n_reference"]),
            "metric/align_n_celltypes_in_ref": float(align_stats["n_celltypes_in_ref"]),
            "metric/align_n_sex_in_ref": float(align_stats["n_sex_in_ref"]),
            "metric/white_active_error": float(align_stats["white_active_error"]),
            "metric/white_removed_leakage": float(align_stats["white_removed_leakage"]),
            "metric/white_cross_leakage": float(align_stats["white_cross_leakage"]),
            "metric/white_target_rank": float(align_stats["white_target_rank"]),
            "metric/white_removed_rank": float(align_stats["white_removed_rank"]),
            "metric/white_identity_target_loss": float(
                align_stats["white_identity_target_loss"]
            ),
            "metric/white_target_loss_delta": float(
                align_stats["white_target_loss_delta"]
            ),
            "weight/lambda_align": float(self.lambda_align),
            "weight/lambda_white": float(self.lambda_white),
            "weight/lambda_align_sex": float(self.lambda_align_sex),
            "weight/align_cov_use_ema": float(self.align_cov_use_ema),
            "weight/align_ema_decay": float(self.align_ema_decay),
            "weight/bam_w_mean": float(weights_per_cell.detach().mean().cpu()),
            "weight/bam_w_min": float(weights_per_cell.detach().min().cpu()),
            "weight/bam_w_max": float(weights_per_cell.detach().max().cpu()),
            "unc/n_valid_tokens": float(n_valid_tokens),
        }
        if raw_uncertainty is not None:
            details["unc/raw_mean"] = float(raw_uncertainty.detach().mean().cpu())
            details["unc/raw_min"] = float(raw_uncertainty.detach().min().cpu())
            details["unc/raw_max"] = float(raw_uncertainty.detach().max().cpu())
        if clipped_uncertainty is not None:
            details["unc/clip_mean"] = float(clipped_uncertainty.detach().mean().cpu())
            details["unc/clip_min"] = float(clipped_uncertainty.detach().min().cpu())
            details["unc/clip_max"] = float(clipped_uncertainty.detach().max().cpu())

        return LossOutput(
            total=total,
            main=main,
            rec=rec.detach(),
            rec_weighted=rec_weighted.detach(),
            kl_state=kl_state.detach(),
            bam_kl=bam_kl.detach(),
            reg_tech=reg_tech.detach(),
            reg_gauge=reg_gauge.detach(),
            cls=cls.detach(),
            ref_center=ref_center.detach(),
            ref_state=ref_state.detach(),
            align=align.detach(),
            white=white.detach(),
            align_sex=align_sex.detach(),
            weights_per_cell=weights_per_cell.detach(),
            raw_uncertainty=None if raw_uncertainty is None else raw_uncertainty.detach(),
            clipped_uncertainty=None if clipped_uncertainty is None else clipped_uncertainty.detach(),
            rec_per_cell=rec_per_cell.detach(),
            rec_weighted_per_cell=rec_weighted_per_cell.detach(),
            kl_state_per_cell=kl_state_per_cell.detach(),
            tech_score_mag=tech_score_mag.detach(),
            details=details,
        )
