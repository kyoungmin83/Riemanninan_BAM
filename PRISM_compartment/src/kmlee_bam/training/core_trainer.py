from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Optional
from contextlib import nullcontext

from pathlib import Path
import json
import math
import os
import time
import sys

# ----------------------------------------------------------------------
# Import path bootstrap
# ----------------------------------------------------------------------
# This file is part of the integrated KMLEE-BAM package.  The public imports
# below use role-based module names rather than historical version folders.
_THIS_DIR = Path(__file__).resolve().parent
_SRC_DIR = _THIS_DIR.parents[1]
_PROJECT_ROOT = _SRC_DIR.parent

for _p in (_PROJECT_ROOT, _SRC_DIR):
    _s = str(_p)
    if _s not in sys.path:
        sys.path.insert(0, _s)

import torch
import torch.distributed as dist
import torch.nn as nn
from torch.nn.utils import clip_grad_norm_

from kmlee_bam.modules.gene_embedding import GeneExpressionEmbedding
from kmlee_bam.model.state_encoder import StateEncoder, StateEncoderOutput
from kmlee_bam.modules.transformer_encoder import TransformerEncoder
from kmlee_bam.model.celltype_prior import CellTypePrior
from kmlee_bam.objectives.total_loss import TotalLoss, LossOutput, empirical_group_weights

from kmlee_bam.modules.module_tokenizer import (
    GeneModuleTokenizer,
    GeneModuleTokenizerOutput,
)
from kmlee_bam.model.lie_ordinal_decoder import (
    LieActionOrdinalDecoder,
    LieActionOrdinalDecoderOutput,
)
from kmlee_bam.model.precision_medicine import (
    PrecisionMedicineOutput,
    delayed_epoch_ramp,
    donor_balanced_gradient,
)
from kmlee_bam.objectives.support_query_infonce import (
    anchored_support_query_infonce,
)
from kmlee_bam.model.pathology_interactions import PATHOLOGY_PAIR_NAMES
from kmlee_bam.training.learned_generator_count import ConstraintBatch


PRISM_NLL_NUM_PER_EXAMPLE_KEY = "stat/prism_branch_nll_num_per_example"
PRISM_NLL_DEN_PER_EXAMPLE_KEY = "stat/prism_branch_nll_den_per_example"
PRISM_RAW_NLL_NUM_PER_EXAMPLE_KEY = "stat/prism_branch_raw_nll_num_per_example"
PRISM_BALANCED_NLL_NUM_PER_EXAMPLE_KEY = (
    "stat/prism_branch_balanced_nll_num_per_example"
)
PRISM_NONZERO_NLL_NUM_PER_EXAMPLE_KEY = (
    "stat/prism_branch_nonzero_nll_num_per_example"
)


def finalize_prism_epoch_metrics(metrics: Dict[str, float]) -> Dict[str, float]:
    """Replace batch-ratio averaging by the exact epoch-wide weighted ratio."""

    numerator = metrics.get(PRISM_NLL_NUM_PER_EXAMPLE_KEY)
    denominator = metrics.get(PRISM_NLL_DEN_PER_EXAMPLE_KEY)
    if numerator is not None and denominator is not None and float(denominator) > 0.0:
        metrics["loss/prism_branch_nll"] = float(numerator) / float(denominator)
        for output_key, numerator_key in (
            ("loss/prism_branch_raw_nll", PRISM_RAW_NLL_NUM_PER_EXAMPLE_KEY),
            ("loss/prism_branch_balanced_nll", PRISM_BALANCED_NLL_NUM_PER_EXAMPLE_KEY),
            ("loss/prism_branch_nonzero_nll", PRISM_NONZERO_NLL_NUM_PER_EXAMPLE_KEY),
        ):
            component_numerator = metrics.get(numerator_key)
            if component_numerator is not None:
                metrics[output_key] = float(component_numerator) / float(denominator)
    return metrics

def _resolve_amp_dtype(amp_dtype: Any) -> Optional[torch.dtype]:
    if amp_dtype is None:
        return None

    if isinstance(amp_dtype, torch.dtype):
        return amp_dtype

    name = str(amp_dtype).lower().strip()

    if name in {"bf16", "bfloat16", "torch.bfloat16"}:
        return torch.bfloat16
    if name in {"fp16", "float16", "half", "torch.float16"}:
        return torch.float16
    if name in {"fp32", "float32", "none", "false", "off"}:
        return None

    raise ValueError(f"Unknown amp_dtype: {amp_dtype}")
# ======================================================================
# Compatibility patch
# ======================================================================
def _patch_transformer_encoder_for_state_encoder_compat() -> None:
    """
    Some local versions of StateEncoder call TransformerEncoder.forward with
    `return_last_attn_uncertainty`, while some local versions of
    TransformerEncoder do not yet expose that argument.

    This patch makes trainer.py robust to that mild interface drift without
    changing the mathematical behavior.
    """
    import inspect

    sig = inspect.signature(TransformerEncoder.forward)
    if "return_last_attn_uncertainty" in sig.parameters:
        return

    original_forward = TransformerEncoder.forward

    def wrapped_forward(
        self,
        x: torch.Tensor,
        *,
        attn_mask: Optional[torch.Tensor] = None,
        key_padding_mask: Optional[torch.Tensor] = None,
        n_samples: int = 1,
        return_all_hidden_states: bool = False,
        return_attn_diagnostics: bool = False,
        return_last_attn_uncertainty: bool = False,
    ):
        need_attn_diag = return_attn_diagnostics or return_last_attn_uncertainty
        out = original_forward(
            self,
            x,
            attn_mask=attn_mask,
            key_padding_mask=key_padding_mask,
            n_samples=n_samples,
            return_all_hidden_states=return_all_hidden_states,
            return_attn_diagnostics=need_attn_diag,
        )

        last_unc = None
        if need_attn_diag and getattr(out, "all_attn_uncertainty", None) is not None:
            all_unc = out.all_attn_uncertainty
            if len(all_unc) > 0:
                last_unc = all_unc[-1]

        setattr(out, "last_attn_uncertainty", last_unc)

        if return_last_attn_uncertainty and not return_attn_diagnostics:
            out.all_attn_weights = None
            out.all_posterior_mean_attn = None
            out.all_prior_mean_attn = None
            out.all_attn_uncertainty = None

        return out

    TransformerEncoder.forward = wrapped_forward


_patch_transformer_encoder_for_state_encoder_compat()

# ======================================================================
# Dataclasses
# ======================================================================
@dataclass
class ModelForwardOutput:
    """
    Rich forward container for the module-token ordinal scRNA model.

    Important convention
    --------------------
    tokens is the token sequence actually passed into StateEncoder. In this
    version that means [CLS, module tokens] when preserve_cls=True, not raw
    gene-level tokens. gene_tokens are kept separately for diagnostics.
    """

    # Encoder input after module tokenization.
    tokens: torch.Tensor
    token_padding_mask: Optional[torch.Tensor]

    # Diagnostics around the gene -> module lift.
    gene_tokens: torch.Tensor
    module_tokens: torch.Tensor
    pooled_gene_state: torch.Tensor
    module_activity: Optional[torch.Tensor]
    module_attention: Optional[torch.Tensor]

    # Latent-state path.
    encoder_out: StateEncoderOutput
    mu_p: torch.Tensor
    logvar_p: torch.Tensor
    sigma_p: torch.Tensor
    z_perp: torch.Tensor

    # Gene-level ordinal reconstruction path. ``None`` only on an encoder-only
    # forward (``OrdinalBAMSystem.forward(..., compute_decoder=False)``), used by
    # the v15 aux celltype-adversary step; always present on the normal path.
    decoder_out: Optional[LieActionOrdinalDecoderOutput]
    cls_logits: Optional[torch.Tensor]
    precision_out: Optional[PrecisionMedicineOutput] = None


@dataclass
class StepOutput:
    """
    Output of one train/eval step.

    Attributes
    ----------
    loss : LossOutput
        Full loss breakdown.
    grad_norm : float or None
        Gradient norm after clipping (train only).
    lr : float
        Current learning rate from the first optimizer group.
    batch_size : int
        Number of cells in the batch.
    """

    loss: LossOutput
    grad_norm: Optional[float]
    lr: float
    batch_size: int


@dataclass
class EpochOutput:
    """
    Aggregated epoch metrics.
    """

    metrics: Dict[str, float] = field(default_factory=dict)
    n_examples: int = 0
    elapsed_sec: float = 0.0


# ======================================================================
# Full system wrapper
# ======================================================================
class OrdinalBAMSystem(nn.Module):
    """
    Full module-token BAM system.

    Composition
    -----------
        GeneExpressionEmbedding
        -> GeneModuleTokenizer
        -> StateEncoder
        -> CellTypePrior.residualize
        -> LieActionOrdinalDecoder

    The encoder now sees module tokens, while the decoder still reconstructs
    gene-level ordinal observations. This preserves gene-level reconstruction
    pressure while letting the latent state be inferred from structured modules.
    """

    def __init__(
        self,
        gene_embedding: GeneExpressionEmbedding,
        module_tokenizer: GeneModuleTokenizer,
        state_encoder: StateEncoder,
        prior: CellTypePrior,
        decoder: LieActionOrdinalDecoder,
        *,
        classifier_head: Optional[nn.Module] = None,
    ) -> None:
        super().__init__()
        self.gene_embedding = gene_embedding
        self.module_tokenizer = module_tokenizer
        self.state_encoder = state_encoder
        self.prior = prior
        self.decoder = decoder
        self.classifier_head = classifier_head

        self._validate_component_compatibility()

    def _validate_component_compatibility(self) -> None:
        if self.state_encoder.condition_on_celltype and self.state_encoder.celltype_embedding is None:
            raise ValueError("StateEncoder is configured to condition on cell type but has no celltype embedding.")

        if self.gene_embedding.d_model != self.module_tokenizer.d_model:
            raise ValueError(
                "d_model mismatch between GeneExpressionEmbedding and GeneModuleTokenizer."
            )

        if self.state_encoder.d_model != self.module_tokenizer.d_model:
            raise ValueError(
                "d_model mismatch between GeneModuleTokenizer and StateEncoder."
            )

        if self.module_tokenizer.n_genes != self.gene_embedding.n_genes:
            raise ValueError(
                "Gene dimension mismatch between GeneExpressionEmbedding and GeneModuleTokenizer. "
                "The registry gene-union order must match the dataset/spec gene order."
            )

        if self.decoder.n_genes != self.gene_embedding.n_genes:
            raise ValueError(
                "Gene dimension mismatch between decoder and GeneExpressionEmbedding."
            )

        if self.gene_embedding.n_bins != self.decoder.n_bins:
            raise ValueError(
                "Ordinal bin count mismatch between GeneExpressionEmbedding and decoder."
            )

        if self.decoder.d_z != self.prior.d_z or self.state_encoder.d_z != self.prior.d_z:
            raise ValueError(
                "Latent dimension mismatch among StateEncoder, CellTypePrior, and decoder."
            )

        if self.module_tokenizer.preserve_cls and not self.gene_embedding.use_cls:
            raise ValueError(
                "GeneModuleTokenizer.preserve_cls=True requires GeneExpressionEmbedding.use_cls=True."
            )

        if self.state_encoder.pooling == "cls" and not self.module_tokenizer.preserve_cls:
            raise ValueError(
                "StateEncoder uses CLS pooling, but GeneModuleTokenizer.preserve_cls=False."
            )

    @property
    def n_effective_encoder_tokens_for_uncertainty(self) -> int:
        """
        Effective token count for BAM entropy normalisation.

        This is not consumed by the current TotalLoss yet, but it is exposed so
        losses.py can use module-token length rather than decoder.n_genes.
        """
        n = int(self.module_tokenizer.n_modules)
        if self.module_tokenizer.preserve_cls and not (
            self.state_encoder.compute_cell_uncertainty
            and self.state_encoder.uncertainty_exclude_cls
        ):
            n += 1
        return n

    def _generator_gate_for_forward(self) -> Optional[torch.Tensor]:
        """Resolve the live joint-count gate for one ordinary model forward.

        Training uses exact 0/1 Binary-Concrete samples with a straight-through
        gradient.  Evaluation uses the deterministic hard mask.  Before the
        configured start epoch every generator is exactly on, while a zero
        dependency on ``log_alpha`` keeps DDP's parameter-usage contract
        intact without moving the logits during warm-up.
        """

        gate = getattr(self, "generator_count_gate", None)
        cfg = getattr(self, "generator_count_config", None)
        if gate is None or cfg is None or not bool(getattr(cfg, "enabled", False)):
            return None
        if str(getattr(cfg, "mode", "gate_only")) != "joint":
            return None

        epoch_state = getattr(self, "generator_count_epoch_state", None)
        epoch = 1 if epoch_state is None else int(epoch_state.item())
        start_epoch = int(cfg.start_epoch)
        end_epoch = max(int(cfg.base_end_epoch), start_epoch)
        progress = (
            0.0
            if epoch <= start_epoch
            else min(
                1.0,
                (epoch - start_epoch) / max(end_epoch - start_epoch, 1),
            )
        )
        temperature = float(gate.temperature(progress))

        shadow_end = int(getattr(cfg, "shadow_end_epoch", 0) or 0)
        soft_end = int(getattr(cfg, "soft_end_epoch", 0) or 0)
        if epoch < start_epoch:
            gate_mode = "warmup_all_on"
            hard = torch.ones_like(gate.log_alpha).detach()
            forward_gate = torch.ones_like(gate.log_alpha) + 0.0 * gate.log_alpha
            live_gate = forward_gate
            expected_active = hard.sum()
        elif bool(
            getattr(
                self,
                "generator_count_force_deterministic_straight_through",
                False,
            )
        ) and torch.is_grad_enabled():
            gate_mode = "rescue_deterministic_candidate"
            sample = gate.deterministic_straight_through(
                temperature=temperature
            )
            hard = sample.hard
            forward_gate = sample.straight_through
            live_gate = forward_gate
            expected_active = gate.expected_active_count(
                temperature=temperature
            )
        elif shadow_end >= start_epoch and epoch <= shadow_end:
            gate_mode = "shadow_all_on"
            if self.training and torch.is_grad_enabled():
                sample = gate.sample(temperature=temperature)
                hard = sample.hard
                forward_gate = sample.straight_through
            else:
                hard = gate.deterministic_mask(temperature=temperature)
                probability = gate.keep_probability(temperature=temperature)
                forward_gate = hard - probability.detach() + probability
            live_gate = torch.ones_like(gate.log_alpha) + 0.0 * gate.log_alpha
            expected_active = gate.expected_active_count(temperature=temperature)
        elif soft_end > 0 and epoch <= soft_end:
            gate_mode = "soft_adaptation"
            probability = gate.keep_probability(temperature=temperature)
            hard = gate.deterministic_mask(temperature=temperature)
            forward_gate = probability
            live_gate = probability
            expected_active = gate.expected_active_count(temperature=temperature)
        elif (
            self.training
            and torch.is_grad_enabled()
            and not bool(
                getattr(self, "generator_count_force_deterministic", False)
            )
        ):
            gate_mode = "protected_hard"
            sample = gate.sample(temperature=temperature)
            hard = sample.hard
            forward_gate = sample.straight_through
            live_gate = forward_gate
            expected_active = gate.expected_active_count(
                temperature=temperature
            )
        else:
            gate_mode = "protected_hard"
            hard = gate.deterministic_mask(temperature=temperature)
            forward_gate = hard
            live_gate = forward_gate
            expected_active = gate.expected_active_count(
                temperature=temperature
            )

        if not bool(((hard == 0) | (hard == 1)).all()):
            raise RuntimeError("joint generator forward mask is not exactly binary")
        self._last_generator_gate = forward_gate
        self._last_generator_hard_mask = hard.detach()
        self._last_generator_expected_active = expected_active
        self._last_generator_temperature = temperature
        self._last_generator_gate_mode = gate_mode
        return live_gate

    def forward(
        self,
        batch: Dict[str, torch.Tensor],
        *,
        sample_latent: bool = True,
        return_all_hidden_states: bool = False,
        return_attn_diagnostics: bool = False,
        compute_decoder: bool = True,
    ) -> ModelForwardOutput:
        # ``compute_decoder=False`` runs the ENCODER->PRIOR path only and SKIPS
        # the gene-level ordinal decoder (the single most expensive component:
        # an ordinal head over all G genes), returning ``decoder_out=None``. This
        # is used by the v15 auxiliary celltype-adversary step, which only needs
        # ``encoder_out.z_s`` / ``mu_p`` / ``sigma_p`` to form the σ-detached
        # z_perp + feed the adversary head; the decoder would be pure waste every
        # aux step. The decoder call is cleanly isolated below, so the rest of
        # the forward (and ``compute_decoder=True``) is byte-identical.
        y_ord = batch["y_ord"]
        celltype_id = batch["celltype_id"]
        tech_id = _resolve_tech_id(batch)
        x_log1p = batch.get("x_log1p")
        x_gene_scalar = batch.get("x_gene_scalar")
        gene_ids = batch.get("gene_ids")

        # Prefer module-level masks if they exist. Gene-level masks are not
        # valid after gene -> module tokenisation because the sequence length
        # has changed from G(+CLS) to M(+CLS).
        module_attn_mask = batch.get("module_attn_mask", batch.get("attn_mask"))

        gene_tokens = self.gene_embedding(
            y_ord=y_ord,
            x_log1p=x_log1p,
            gene_ids=gene_ids,
        )

        module_out: GeneModuleTokenizerOutput = self.module_tokenizer(
            gene_tokens=gene_tokens,
            x_gene_scalar=x_gene_scalar,
            return_attention=return_attn_diagnostics,
        )

        tokens = module_out.tokens
        token_padding_mask = _resolve_module_token_padding_mask(
            batch=batch,
            B=tokens.shape[0],
            T=tokens.shape[1],
            device=tokens.device,
        )

        enc_out = self.state_encoder(
            tokens,
            celltype_id=celltype_id if self.state_encoder.condition_on_celltype else None,
            attn_mask=module_attn_mask,
            key_padding_mask=token_padding_mask,
            sample_latent=sample_latent,
            return_all_hidden_states=return_all_hidden_states,
            return_attn_diagnostics=return_attn_diagnostics,
        )

        mu_p, logvar_p, sigma_p = self.prior(celltype_id)
        z_perp = self.prior.residualize(enc_out.z_s, celltype_id)
        # v25: STRUCTURAL nuisance projection-out (sex / celltype) — default-off identity.
        # Applied at the source so EVERY downstream user (decoder state, path_aux head, and
        # the trainer's align/whiten losses on model_out.z_perp) sees the CLEANED z_perp.
        # The decoder still receives sex_id below ⇒ the sex BASELINE carries the main effect
        # while z_perp is structurally barred from re-encoding it (⟂ whitening, no tug-of-war).
        _nproj = getattr(self, "nuisance_projector", None)
        if _nproj is not None:
            z_perp = _nproj(
                z_perp,
                sex_id=batch.get("sex_id"),
                celltype_id=celltype_id,
                donor_id=batch.get("donor_id"),
                update=self.training,
            )

        decoder_out = None
        precision_out = None
        if compute_decoder:
            generator_gate = self._generator_gate_for_forward()
            precision_score = None
            precision_score_for_full = None
            precision_head = getattr(self, "precision_head", None)
            if (
                precision_head is not None
                and getattr(precision_head.cfg, "enabled", False)
                and bool(getattr(self, "precision_forward_enabled", True))
            ):
                required = (
                    "donor_id",
                    "celltype_id",
                    "region_id",
                    "age_z",
                    "age_valid",
                    "prism_pathology",
                    "prism_pathology_valid",
                )
                missing = [key for key in required if key not in batch]
                if missing:
                    raise KeyError(f"PRISM E2E batch is missing required fields: {missing}")
                precision_out = precision_head(
                    donor_id=batch["donor_id"],
                    celltype_id=celltype_id,
                    region_id=batch["region_id"],
                    age_z=batch["age_z"],
                    age_valid=batch["age_valid"],
                    pathology=batch["prism_pathology"],
                    pathology_valid=batch["prism_pathology_valid"],
                    state_latent=z_perp,
                    cell_weight=batch.get("donor_balance_weight"),
                )
                lift = self.module_tokenizer.activity_weight
                if precision_out.total_module_coeff.shape[1] != lift.shape[0]:
                    raise RuntimeError(
                        "PRISM module coefficient/dictionary mismatch: "
                        f"{precision_out.total_module_coeff.shape[1]} vs {lift.shape[0]}"
                    )
                # One linear lift for the sum.  Individual common/personal/
                # response module coefficients remain explicit for logging and
                # fingerprint extraction; linearity makes their gene-score sum exact.
                precision_score = precision_out.total_module_coeff @ lift.to(
                    device=precision_out.total_module_coeff.device,
                    dtype=precision_out.total_module_coeff.dtype,
                )
                precision_out.gene_score = precision_score
                if bool(getattr(precision_head.cfg, "module_local_enabled", False)):
                    local_scale = precision_head.explicit_scale.to(
                        device=precision_out.module_local_coeff.device,
                        dtype=precision_out.module_local_coeff.dtype,
                    )
                    precision_out.module_local_gene_score = (
                        precision_out.module_local_coeff * local_scale
                    ) @ lift.to(
                        device=precision_out.module_local_coeff.device,
                        dtype=precision_out.module_local_coeff.dtype,
                    )
                    if bool(
                        getattr(
                            precision_head.cfg,
                            "module_local_nonlinear_enabled",
                            False,
                        )
                    ):
                        precision_out.module_local_nonlinear_gene_score = (
                            precision_out.module_local_nonlinear_coeff
                            * local_scale
                        ) @ lift.to(
                            device=(
                                precision_out.module_local_nonlinear_coeff.device
                            ),
                            dtype=(
                                precision_out.module_local_nonlinear_coeff.dtype
                            ),
                        )
                donor_weight = batch.get("donor_balance_weight")
                if donor_weight is None:
                    donor_weight = torch.ones(
                        precision_score.shape[0],
                        device=precision_score.device,
                        dtype=precision_score.dtype,
                    )
                # Forward values are identical to ``precision_score``.  Only
                # the full-reconstruction gradient entering the PRISM head is
                # scaled so every donor has equal total influence over an epoch.
                precision_score_for_full = donor_balanced_gradient(
                    precision_score, donor_weight
                )

            decoder_out = self.decoder(
                z_perp=z_perp,
                celltype_id=celltype_id,
                tech_id=tech_id,
                y_ord=y_ord,
                # Biological-covariate factorization: the decoder adds an explicit
                # SEX baseline when built with n_sex (else this is ignored). None
                # when the batch carries no sex label ⇒ no-op.
                sex_id=batch.get("sex_id", None),
                precision_score=precision_score_for_full,
                generator_gate=generator_gate,
            )
            if precision_out is not None:
                local_gene_score = precision_out.module_local_gene_score
                if local_gene_score is not None:
                    with torch.no_grad():
                        (
                            _,
                            precision_out.module_local_off_full_nll_per_cell,
                        ) = self.decoder.nll_from_external_score(
                            decoder_out.score.detach()
                            - local_gene_score.detach(),
                            y_ord,
                            detach_thresholds=True,
                        )
                nonlinear_gene_score = (
                    precision_out.module_local_nonlinear_gene_score
                )
                if nonlinear_gene_score is not None:
                    with torch.no_grad():
                        (
                            _,
                            precision_out.module_local_nonlinear_off_full_nll_per_cell,
                        ) = self.decoder.nll_from_external_score(
                            decoder_out.score.detach()
                            - nonlinear_gene_score.detach(),
                            y_ord,
                            detach_thresholds=True,
                        )
                # Target-latent-free likelihood: the target cell's encoder z,
                # and legacy state decoder are absent.  Known technical and sex
                # baselines are included as detached nuisance adjustments.  The
                # remaining biology beyond the common cell-type normal is
                # inferred from other cell types of the same donor plus the five
                # named pathology axes.  This prevents the unrestricted
                # autoencoder route from making the explicit branches decorative.
                branch_score = (
                    self.decoder.celltype_baseline(celltype_id).detach()
                    + self.decoder.tech_baseline(tech_id).detach()
                    + precision_score
                )
                if self.decoder.sex_baseline is not None and batch.get("sex_id") is not None:
                    branch_score = branch_score + self.decoder.sex_baseline(
                        self.decoder._sex_index(batch["sex_id"])
                    ).detach()
                (
                    precision_out.branch_nll_per_gene,
                    precision_out.branch_nll_per_cell,
                ) = self.decoder.nll_from_external_score(
                    branch_score, y_ord, detach_thresholds=True
                )
                if local_gene_score is not None:
                    with torch.no_grad():
                        (
                            _,
                            precision_out.module_local_off_branch_nll_per_cell,
                        ) = self.decoder.nll_from_external_score(
                            branch_score.detach() - local_gene_score.detach(),
                            y_ord,
                            detach_thresholds=True,
                        )
                if nonlinear_gene_score is not None:
                    with torch.no_grad():
                        (
                            _,
                            precision_out.module_local_nonlinear_off_branch_nll_per_cell,
                        ) = self.decoder.nll_from_external_score(
                            branch_score.detach()
                            - nonlinear_gene_score.detach(),
                            y_ord,
                            detach_thresholds=True,
                        )

                # Optional support -> query discrimination.  It is evaluated
                # inside the DDP-wrapped model forward so every repeated use of
                # the personal context encoder remains visible to DDP.  The
                # expensive counterfactuals are training-only; validation keeps
                # the ordinary held-out reconstruction objective unchanged.
                pm_cfg = precision_head.cfg
                if (
                    self.training
                    and float(pm_cfg.lambda_support_query_infonce) > 0.0
                ):
                    negative_donor, negative_valid, match_distance = (
                        precision_head.select_matched_negative_donors(
                            batch["donor_id"], celltype_id
                        )
                    )
                    batch_size, n_negative = negative_donor.shape
                    flat_celltype = celltype_id.unsqueeze(1).expand(
                        batch_size, n_negative
                    ).reshape(-1)
                    negative_code, _, _ = precision_head.infer_personal_code(
                        negative_donor.reshape(-1), flat_celltype
                    )
                    flat_pathology = precision_out.pathology.unsqueeze(1).expand(
                        batch_size, n_negative, -1
                    ).reshape(-1, precision_head.n_pathology)
                    flat_pathology_valid = (
                        precision_out.pathology_valid.unsqueeze(1)
                        .expand(batch_size, n_negative, -1)
                        .reshape(-1, precision_head.n_pathology)
                    )
                    negative_personal, negative_response = (
                        precision_head.personal_terms_from_code(
                            negative_code,
                            flat_celltype,
                            flat_pathology,
                            flat_pathology_valid,
                        )
                    )
                    if bool(
                        getattr(pm_cfg, "module_local_enabled", False)
                    ):
                        (
                            negative_local_raw,
                            _,
                            _,
                            _,
                        ) = precision_head.infer_module_local(
                            negative_donor.reshape(-1), flat_celltype
                        )
                        negative_local_scale = precision_head.module_local_scale.to(
                            device=negative_local_raw.device,
                            dtype=negative_local_raw.dtype,
                        )
                        negative_personal = negative_personal + (
                            negative_local_raw * negative_local_scale
                        )
                    negative_personal = negative_personal.reshape(
                        batch_size, n_negative, -1
                    )
                    negative_response = negative_response.reshape(
                        batch_size,
                        n_negative,
                        precision_head.n_pathology,
                        -1,
                    )
                    common_module = (
                        precision_out.total_module_coeff
                        - precision_out.personal_coeff
                        - precision_out.response_axis_coeff.sum(dim=1)
                    )
                    negative_module = (
                        common_module.unsqueeze(1)
                        + negative_personal
                        + negative_response.sum(dim=2)
                    )
                    candidate_module = torch.cat(
                        (
                            precision_out.total_module_coeff.unsqueeze(1),
                            negative_module,
                        ),
                        dim=1,
                    )
                    candidate_gene = candidate_module @ lift.to(
                        device=candidate_module.device,
                        dtype=candidate_module.dtype,
                    )
                    nuisance_score = (
                        self.decoder.celltype_baseline(celltype_id).detach()
                        + self.decoder.tech_baseline(tech_id).detach()
                    )
                    if (
                        self.decoder.sex_baseline is not None
                        and batch.get("sex_id") is not None
                    ):
                        nuisance_score = nuisance_score + self.decoder.sex_baseline(
                            self.decoder._sex_index(batch["sex_id"])
                        ).detach()
                    n_candidates = 1 + n_negative
                    flat_score = (
                        nuisance_score.unsqueeze(1) + candidate_gene
                    ).reshape(batch_size * n_candidates, -1)
                    flat_y = y_ord.unsqueeze(1).expand(
                        batch_size, n_candidates, -1
                    ).reshape(batch_size * n_candidates, -1)
                    _, candidate_nll = self.decoder.nll_from_external_score(
                        flat_score,
                        flat_y,
                        detach_thresholds=True,
                    )
                    common_gene = common_module @ lift.to(
                        device=common_module.device,
                        dtype=common_module.dtype,
                    )
                    _, common_nll = self.decoder.nll_from_external_score(
                        nuisance_score + common_gene,
                        y_ord,
                        detach_thresholds=True,
                    )
                    precision_out.support_query_candidate_nll = candidate_nll.reshape(
                        batch_size, n_candidates
                    )
                    precision_out.support_query_negative_valid = negative_valid
                    precision_out.support_query_common_nll = common_nll
                    precision_out.support_query_negative_donor = negative_donor
                    precision_out.support_query_match_distance = match_distance

        cls_logits = None
        if self.classifier_head is not None:
            cls_logits = self.classifier_head(enc_out.pooled_state)

        return ModelForwardOutput(
            tokens=tokens,
            token_padding_mask=token_padding_mask,
            gene_tokens=gene_tokens,
            module_tokens=module_out.module_tokens,
            pooled_gene_state=module_out.pooled_gene_state,
            module_activity=module_out.module_activity,
            module_attention=module_out.module_attention,
            encoder_out=enc_out,
            mu_p=mu_p,
            logvar_p=logvar_p,
            sigma_p=sigma_p,
            z_perp=z_perp,
            decoder_out=decoder_out,
            cls_logits=cls_logits,
            precision_out=precision_out,
        )


# ======================================================================
# Trainer
# ======================================================================
class Trainer:
    """
    Project-ready trainer for the ordinal BAM-VAE.

    Design choices
    --------------
    1) The model system stays modular: embedding / encoder / prior / decoder are
       still exposed and can be inspected independently.
    2) `TotalLoss` remains the single source of truth for the mathematical
       objective.
    3) Mixed precision is supported but disabled automatically on CPU.
    4) BAM uncertainty weights only the reconstruction term, exactly as encoded
       in the hybrid loss module.
    5) Gradient accumulation is handled at the epoch loop level so that DDP, AMP,
       clipping, and scheduler stepping all occur on true optimizer updates.
    """

    def __init__(
        self,
        system: OrdinalBAMSystem,
        criterion: TotalLoss,
        optimizer: torch.optim.Optimizer,
        *,
        device: Optional[torch.device | str] = None,
        scheduler: Optional[Any] = None,
        grad_clip_norm: Optional[float] = 1.0,
        amp: bool = False,
        amp_dtype: torch.dtype = torch.float16,
        tech_weight_mode: str = "empirical_batch",
        scheduler_step_on: str = "epoch",
        eval_sample_latent: bool = False,
        grad_accum_steps: int = 1,
        debug_z_sensitivity_every: Optional[int] = None,
        debug_z_sensitivity_path: Optional[str | Path] = None,
    ) -> None:
        self.system = system
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.grad_clip_norm = grad_clip_norm
        self.tech_weight_mode = tech_weight_mode
        self.scheduler_step_on = scheduler_step_on
        self.eval_sample_latent = eval_sample_latent
        self.debug_z_sensitivity_every = debug_z_sensitivity_every
        self.debug_z_sensitivity_path = (
            Path(debug_z_sensitivity_path)
            if debug_z_sensitivity_path is not None
            else None
        )

        if tech_weight_mode not in {"uniform", "empirical_batch"}:
            raise ValueError(
                f"tech_weight_mode must be 'uniform' or 'empirical_batch', got '{tech_weight_mode}'."
            )
        if scheduler_step_on not in {"step", "epoch", "none"}:
            raise ValueError(
                f"scheduler_step_on must be 'step', 'epoch', or 'none', got '{scheduler_step_on}'."
            )
        if not isinstance(grad_accum_steps, int) or grad_accum_steps < 1:
            raise ValueError("grad_accum_steps must be an integer >= 1.")

        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.device = torch.device(device)
        self.system.to(self.device)
        self.amp_dtype = _resolve_amp_dtype(amp_dtype)
        self.amp_enabled = bool(amp 
                                and self.device.type == "cuda"
                                and self.amp_dtype is not None)
        
        self.grad_accum_steps = grad_accum_steps
        # v32 stochastic-consistency (R-Drop style) regulariser. OFF unless env KMLEE_PATH_CONSISTENCY>0.
        # Two stochastic forward passes (attention + z sampling) must agree on the donor x celltype
        # disease projection -> a disease route that survives the perturbation. See
        # _maybe_add_path_consistency. KMLEE_PATH_CONS_MINCELLS drops small (donor,celltype) groups
        # (a singleton group == noisy per-cell consistency; pathology is donor-level).
        self.lambda_path_consistency = float(os.environ.get("KMLEE_PATH_CONSISTENCY", "0.0"))
        self.path_cons_min_cells = int(os.environ.get("KMLEE_PATH_CONS_MINCELLS", "2"))
        # Precision-medicine decoder anchor. OFF unless explicitly enabled.
        #
        # The latent prior already anchors z_perp around reference cells, but the
        # decoder celltype_baseline is still a learned mixed-cohort intercept.
        # These losses push the decoder interpretation toward:
        #   base+tech(+sex) = normal/reference reconstruction path
        #   state_score     = deviation path
        # without changing any default run.
        self.lambda_decoder_ref_state = float(os.environ.get("KMLEE_DEC_REF_STATE", "0.0"))
        self.lambda_decoder_ref_base_nll = float(os.environ.get("KMLEE_DEC_REF_BASE_NLL", "0.0"))
        self.decoder_ref_min_cells = int(os.environ.get("KMLEE_DEC_REF_MINCELLS", "4"))
        self.scaler = torch.amp.GradScaler(
            "cuda",
            enabled=(self.amp_enabled and self.amp_dtype is torch.float16),
        )

        self._last_step_was_skipped = False

    @staticmethod
    def _is_rank0() -> bool:
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            return torch.distributed.get_rank() == 0
        return True

    # ------------------------------------------------------------------
    # Public step methods
    # ------------------------------------------------------------------
    def _maybe_add_path_consistency(self, batch, model_out, loss_out):
        """v32 stochastic-consistency regulariser (R-Drop style). OFF unless KMLEE_PATH_CONSISTENCY>0.

        A 2nd encoder-only STOCHASTIC forward (fresh attention + z sample) is run; its donor x
        celltype disease projection (the pathology head) must match the first -> a disease route
        that survives the stochastic perturbation. NOT attention-only: sample_latent mixes attention
        AND z-posterior noise, so this is 'stochastic representation' consistency.
          - compute_decoder=False : skips decoder + adversaries + the pathology_aux EMA-bank update.
          - nuisance_projector.eval() : EMA scatter is NOT double-updated (its _update is gated by
            self.training) -> no running-state corruption from the 2nd pass.
          - 2nd pass under no_grad : mean-teacher, so its graph is never built (saves memory/time).
          - group at donor x celltype, keeping only groups with >= path_cons_min_cells cells; a
            singleton group == per-cell consistency (noisy, since pathology is donor-level).
          - d1=head(z1) carries grad; d2=head(z2) is a frozen target (mean-teacher).
        COLLAPSE WATCH: head is trainable via d1, so a too-large lambda can shrink the head logit
        scale -> monitor loss/path_consistency vs the pathology loss (use lambda 0.03-0.1 first)."""
        lam = self.lambda_path_consistency
        head_mod = getattr(self.system, "pathology_aux_head", None)
        if (lam <= 0.0
                or getattr(model_out, "z_perp", None) is None
                or head_mod is None or getattr(head_mod, "head", None) is None
                or "donor_id" not in batch or "celltype_id" not in batch):
            return loss_out
        proj = getattr(self.system, "nuisance_projector", None)
        proj_was_training = bool(proj.training) if proj is not None else False
        if proj is not None:
            proj.eval()  # gate off the EMA double-update on the 2nd pass
        try:
            with torch.no_grad():  # mean-teacher: 2nd pass is a frozen target, never built into a graph
                model_out2 = self.system(
                    batch,
                    sample_latent=True,
                    compute_decoder=False,
                    return_all_hidden_states=False,
                    return_attn_diagnostics=False,
                )
        finally:
            if proj is not None and proj_was_training:
                proj.train()
        if getattr(model_out2, "z_perp", None) is None:
            return loss_out
        # (donor, celltype) grouping is identical for both passes (same batch ids) -> compute once.
        z = model_out.z_perp
        key = batch["donor_id"].reshape(-1).long() * 100000 + batch["celltype_id"].reshape(-1).long()
        uniq, inv = torch.unique(key, return_inverse=True)
        n = int(uniq.numel())
        cnt = z.new_zeros(n).index_add_(0, inv, z.new_ones(z.shape[0]))

        def gmean(zz):
            s = zz.new_zeros(n, zz.shape[1]).index_add_(0, inv, zz)
            return s / cnt.clamp_min(1.0).unsqueeze(1)

        keep = cnt >= float(self.path_cons_min_cells)  # drop small groups -> avoid noisy per-cell consistency
        nkept = int(keep.sum().item())
        try:  # diagnostics logged even when the loss is skipped (reveals batch composition)
            loss_out.details["metric/path_cons_n_groups"] = float(n)
            loss_out.details["metric/path_cons_n_kept"] = float(nkept)
            loss_out.details["metric/path_cons_mean_group_cells"] = float(cnt.mean().item())
            loss_out.details["metric/path_cons_singleton_frac"] = float((cnt < 2).float().mean().item())
        except Exception:
            pass
        if nkept == 0:
            return loss_out
        head = head_mod.head
        z1 = gmean(model_out.z_perp)[keep]
        d1 = head(z1)
        with torch.no_grad():  # frozen target: no grad through z2 OR the head
            z2 = gmean(model_out2.z_perp)[keep]
            d2 = head(z2)
        cons = ((d1 - d2) ** 2).mean()
        loss_out.total = loss_out.total + lam * cons
        try:
            loss_out.details["loss/path_consistency"] = float(cons.detach().cpu())
            loss_out.details["weight/lambda_path_consistency"] = float(lam)
            loss_out.details["loss/total"] = float(loss_out.total.detach().cpu())  # keep log/finite in sync (#1)
        except Exception:
            pass
        return loss_out

    def _maybe_add_cp_bam(self, batch, model_out, loss_out):
        """Add the CP-BAM loss. The cp_bam_head is called INSIDE system.forward (so its params are
        DDP-tracked, like pathology_aux); here we ONLY read the precomputed cp_bam_out and add its
        loss. Calling a parametered head from the trainer loss path raises DDP 'mark ready once'."""
        cpb = getattr(model_out, "cp_bam_out", None)
        if cpb is None or getattr(cpb, "loss", None) is None:
            return loss_out
        loss_out.total = loss_out.total + cpb.loss
        try:
            loss_out.details.update(cpb.details)
            loss_out.details["loss/total"] = float(loss_out.total.detach().cpu())  # keep log/finite in sync
        except Exception:
            pass
        return loss_out

    def _maybe_add_precision_medicine(self, batch, model_out, loss_out):
        """Add the donor-balanced target-latent-free PRISM likelihood and guards."""

        precision = getattr(model_out, "precision_out", None)
        head = getattr(self.system, "precision_head", None)
        if precision is None or head is None or not getattr(head.cfg, "enabled", False):
            return loss_out
        if (
            precision.branch_nll_per_cell is None
            or precision.branch_nll_per_gene is None
            or precision.gene_score is None
        ):
            raise RuntimeError("PRISM E2E output is missing its branch likelihood/score")

        cfg = head.cfg
        weight = batch.get("donor_balance_weight")
        if weight is None:
            weight = torch.ones_like(precision.branch_nll_per_cell)
        weight = weight.to(precision.branch_nll_per_cell.dtype).reshape(-1)
        branch_raw_values = precision.branch_nll_per_cell.reshape(-1)
        branch_gene_nll = precision.branch_nll_per_gene
        y_ord = batch["y_ord"].long()
        ordinal_weights = getattr(self, "ordinal_class_weights", None)
        if ordinal_weights is None or int(ordinal_weights.numel()) <= int(y_ord.max()):
            ordinal_weights = torch.ones(
                int(y_ord.max().item()) + 1,
                device=branch_gene_nll.device,
                dtype=branch_gene_nll.dtype,
            )
        else:
            ordinal_weights = ordinal_weights.to(
                device=branch_gene_nll.device, dtype=branch_gene_nll.dtype
            )
        branch_balanced_values = (
            branch_gene_nll * ordinal_weights[y_ord]
        ).mean(dim=1)
        nonzero = y_ord > 0
        branch_nonzero_values = (
            (branch_gene_nll * nonzero.to(branch_gene_nll.dtype)).sum(dim=1)
            / nonzero.sum(dim=1).clamp_min(1).to(branch_gene_nll.dtype)
        )
        branch_values = (
            branch_raw_values
            + float(cfg.branch_balanced_weight) * branch_balanced_values
            + float(cfg.branch_nonzero_weight) * branch_nonzero_values
        )
        weighted_numerator = (branch_values * weight).sum()
        weighted_denominator = weight.sum().clamp_min(1e-8)
        branch_nll = weighted_numerator / weighted_denominator
        # This is the optimisation term.  Dataset weights have split-wide mean
        # one, so a plain per-cell mean gives an unbiased donor-balanced epoch
        # objective even when a batch contains mostly one donor.
        branch_objective = (branch_values * weight).mean()

        current_epoch = getattr(self, "current_epoch_index", None)
        if current_epoch is not None:
            self._prism_last_epoch = int(current_epoch)
        epoch_for_ramp = int(getattr(self, "_prism_last_epoch", 1))
        ramp = delayed_epoch_ramp(
            epoch_for_ramp,
            start_epoch=int(getattr(cfg, "loss_start_epoch", 1) or 1),
            ramp_epochs=int(getattr(cfg, "loss_ramp_epochs", 0) or 0),
        )

        penalties = head.parameter_penalties()
        code = precision.personal_code
        code_l2 = (code.square().mean(dim=1) * weight).mean()
        weighted_code_mean = (code * weight.unsqueeze(1)).sum(dim=0) / weight.sum().clamp_min(
            1e-8
        )
        code_center = weighted_code_mean.square().mean()
        weighted_response_mean = (
            precision.response_axis_coeff * weight[:, None, None]
        ).sum(dim=0) / weight.sum().clamp_min(1e-8)
        response_mean = weighted_response_mean.square().mean()

        branch_term = float(cfg.lambda_branch_nll) * ramp * branch_objective
        is_train = bool(getattr(self.system, "training", False))
        support_query = None
        support_query_added = branch_term.detach() * 0.0
        if is_train and float(cfg.lambda_support_query_infonce) > 0.0:
            if (
                precision.support_query_candidate_nll is None
                or precision.support_query_negative_valid is None
                or precision.support_query_common_nll is None
            ):
                raise RuntimeError(
                    "support-query InfoNCE is enabled but counterfactual NLLs are missing"
                )
            support_query = anchored_support_query_infonce(
                precision.support_query_candidate_nll,
                precision.support_query_negative_valid,
                precision.support_query_common_nll,
                temperature=float(cfg.support_query_temperature),
                anchor_weight=float(cfg.support_query_anchor_weight),
                cell_weight=weight,
            )
            support_query_added = (
                float(cfg.lambda_support_query_infonce) * support_query.loss
            )
        regularizer = (
            float(cfg.lambda_normal_context_l2) * penalties["normal_context"]
            + float(cfg.lambda_common_l2) * penalties["common"]
            + float(cfg.lambda_personal_l2) * penalties["personal"]
            + float(cfg.lambda_response_l2) * penalties["response"]
            + float(cfg.lambda_interaction_l2) * penalties["interaction"]
            + float(cfg.lambda_response_zero_mean) * response_mean
            + float(cfg.lambda_code_l2) * code_l2
            + float(cfg.lambda_code_center) * code_center
            + float(cfg.lambda_pathology_leak) * precision.pathology_leak_loss
            + float(cfg.lambda_age_leak) * precision.age_leak_loss
            + float(cfg.lambda_state_pathology_leak)
            * precision.state_pathology_leak_loss
            + float(cfg.lambda_module_local_center)
            * precision.module_local_center_loss
            + float(cfg.lambda_module_local_hierarchy)
            * penalties["module_local_hierarchy"]
            + float(cfg.lambda_module_local_size)
            * precision.module_local_size_loss
            + float(cfg.lambda_module_local_pathology_leak)
            * precision.module_local_pathology_leak_loss
            + float(cfg.lambda_module_local_age_leak)
            * precision.module_local_age_leak_loss
        )
        # Parameter/code guards train the representation but are excluded from
        # validation loss.  The raw donor-balanced branch NLL remains logged on
        # both splits and is the precision early-stopping criterion.
        regularizer_scale = (
            ramp
            if bool(getattr(cfg, "scale_regularizers_with_loss_ramp", False))
            else 1.0
        )
        added = (
            branch_term
            + support_query_added
            + regularizer_scale
            * (regularizer if is_train else regularizer.detach() * 0.0)
        )
        loss_out.total = loss_out.total + added.to(loss_out.total.dtype)

        def scalar(value) -> float:
            return float(value.detach().float().cpu())

        details = loss_out.details
        details["loss/prism_branch_nll"] = scalar(branch_nll)
        details["loss/prism_branch_raw_nll"] = scalar(
            (branch_raw_values * weight).sum() / weighted_denominator
        )
        details["loss/prism_branch_balanced_nll"] = scalar(
            (branch_balanced_values * weight).sum() / weighted_denominator
        )
        details["loss/prism_branch_nonzero_nll"] = scalar(
            (branch_nonzero_values * weight).sum() / weighted_denominator
        )
        batch_size = max(int(branch_values.numel()), 1)
        details[PRISM_NLL_NUM_PER_EXAMPLE_KEY] = scalar(
            weighted_numerator / float(batch_size)
        )
        details[PRISM_RAW_NLL_NUM_PER_EXAMPLE_KEY] = scalar(
            (branch_raw_values * weight).sum() / float(batch_size)
        )
        details[PRISM_BALANCED_NLL_NUM_PER_EXAMPLE_KEY] = scalar(
            (branch_balanced_values * weight).sum() / float(batch_size)
        )
        details[PRISM_NONZERO_NLL_NUM_PER_EXAMPLE_KEY] = scalar(
            (branch_nonzero_values * weight).sum() / float(batch_size)
        )
        details[PRISM_NLL_DEN_PER_EXAMPLE_KEY] = scalar(
            weighted_denominator / float(batch_size)
        )
        details["loss/prism_regularizer"] = scalar(regularizer)
        details["loss/prism_pathology_leak"] = scalar(precision.pathology_leak_loss)
        details["loss/prism_age_leak"] = scalar(precision.age_leak_loss)
        details["loss/prism_state_pathology_leak"] = scalar(
            precision.state_pathology_leak_loss
        )
        details["loss/prism_module_local_center"] = scalar(
            precision.module_local_center_loss
        )
        details["loss/prism_module_local_hierarchy"] = scalar(
            penalties["module_local_hierarchy"]
        )
        details["loss/prism_module_local_size"] = scalar(
            precision.module_local_size_loss
        )
        details["loss/prism_module_local_pathology_leak"] = scalar(
            precision.module_local_pathology_leak_loss
        )
        details["loss/prism_module_local_age_leak"] = scalar(
            precision.module_local_age_leak_loss
        )
        details["loss/prism_added"] = scalar(added)
        details["loss/prism_support_query_infonce"] = scalar(
            support_query.loss if support_query is not None else support_query_added
        )
        details["loss/prism_support_query_contrastive"] = scalar(
            support_query.contrastive
            if support_query is not None
            else support_query_added
        )
        details["loss/prism_support_query_anchor"] = scalar(
            support_query.anchor if support_query is not None else support_query_added
        )
        details["metric/prism_support_query_valid_fraction"] = scalar(
            support_query.valid_fraction
            if support_query is not None
            else support_query_added
        )
        details["metric/prism_support_query_own_nll"] = scalar(
            support_query.own_nll if support_query is not None else support_query_added
        )
        details["metric/prism_support_query_negative_nll"] = scalar(
            support_query.negative_nll
            if support_query is not None
            else support_query_added
        )
        if support_query is not None and precision.support_query_match_distance is not None:
            finite_match = precision.support_query_match_distance[
                precision.support_query_negative_valid
            ]
            details["metric/prism_support_query_match_distance"] = scalar(
                finite_match.mean()
                if finite_match.numel() > 0
                else support_query_added
            )
        else:
            details["metric/prism_support_query_match_distance"] = 0.0
        details["weight/prism_ramp"] = float(ramp)
        details["weight/prism_regularizer_ramp"] = float(regularizer_scale)
        details["weight/prism_output_ramp"] = float(
            getattr(head, "explicit_scale", torch.ones(())).detach().float().cpu()
        )
        details["metric/prism_normal_region_module_rms"] = scalar(
            precision.normal_region_coeff.float().square().mean().sqrt()
        )
        details["metric/prism_age_module_rms"] = scalar(
            precision.age_coeff.float().square().mean().sqrt()
        )
        details["metric/prism_common_module_rms"] = scalar(
            precision.common_axis_coeff.float().square().mean().sqrt()
        )
        details["metric/prism_personal_module_rms"] = scalar(
            precision.personal_coeff.float().square().mean().sqrt()
        )
        local_module_rms = precision.module_local_coeff.float().square().mean().sqrt()
        local_raw_rms = (
            precision.module_local_raw_coeff.float().square().mean().sqrt()
        )
        rank2_module_rms = (
            precision.personal_rank2_coeff.float().square().mean().sqrt()
        )
        details["metric/prism_module_local_module_rms"] = scalar(
            local_module_rms
        )
        details["metric/prism_module_local_raw_module_rms"] = scalar(
            local_raw_rms
        )
        nonlinear_local_rms = (
            precision.module_local_nonlinear_coeff.float()
            .square()
            .mean()
            .sqrt()
        )
        details["metric/prism_module_local_nonlinear_rms"] = scalar(
            nonlinear_local_rms
        )
        details["metric/prism_module_local_nonlinear_to_local_ratio"] = scalar(
            nonlinear_local_rms / local_module_rms.clamp_min(1.0e-8)
        )
        details["metric/prism_module_local_threshold_crossing_fraction"] = scalar(
            precision.module_local_threshold_crossing_fraction.float().mean()
        )
        details["metric/prism_rank2_personal_module_rms"] = scalar(
            rank2_module_rms
        )
        details["metric/prism_module_local_to_rank2_rms_ratio"] = scalar(
            local_module_rms / rank2_module_rms.clamp_min(1.0e-8)
        )
        details["metric/prism_module_local_support_coverage"] = scalar(
            precision.module_local_support_coverage.float().mean()
        )
        details["metric/prism_module_local_zero_coverage_fraction"] = scalar(
            precision.module_local_zero_coverage_fraction.float().mean()
        )
        details["metric/prism_module_local_personal_span_overlap_mean"] = scalar(
            precision.module_local_personal_span_overlap.float().mean()
        )
        details["metric/prism_module_local_personal_span_overlap_max"] = scalar(
            precision.module_local_personal_span_overlap.float().max()
        )
        local_scale_tensor = getattr(head, "module_local_scale", None)
        details["weight/prism_module_local_ramp"] = (
            scalar(local_scale_tensor)
            if isinstance(local_scale_tensor, torch.Tensor)
            else 0.0
        )
        nonlinear_scale_tensor = getattr(
            head, "module_local_nonlinear_scale", None
        )
        details["weight/prism_module_local_nonlinear_ramp"] = (
            scalar(nonlinear_scale_tensor)
            if isinstance(nonlinear_scale_tensor, torch.Tensor)
            else 0.0
        )
        paper_branch_enabled = bool(
            getattr(head, "module_local_nonlinear_enabled", False)
        )
        paper_branch_variant = str(
            getattr(head, "module_local_nonlinear_variant", "disabled")
        )
        details["metric/prism_paper_branch_enabled"] = float(
            paper_branch_enabled
        )
        for variant_name in (
            "compartmental_threshold",
            "graph_linear_control",
        ):
            details[
                f"metric/prism_paper_branch_variant_{variant_name}"
            ] = float(
                paper_branch_enabled and paper_branch_variant == variant_name
            )
        parameter_diagnostics = getattr(
            head, "module_local_nonlinear_parameter_diagnostics", None
        )
        if callable(parameter_diagnostics):
            for name, value in parameter_diagnostics().items():
                details[f"metric/prism_module_local_{name}"] = scalar(value)

        local_code = precision.module_local_code.detach().float()
        if local_code.ndim == 2 and local_code.shape[1] > 0:
            centered_local_code = local_code - local_code.mean(
                dim=0, keepdim=True
            )
            singular_values = torch.linalg.svdvals(centered_local_code)
            singular_power = singular_values.square()
            power_sum = singular_power.sum()
            participation_rank = power_sum.square() / singular_power.square().sum().clamp_min(
                1.0e-12
            )
            if float(power_sum.cpu()) > 0.0:
                cumulative_power = torch.cumsum(singular_power, dim=0) / power_sum
                rank95 = int(
                    torch.searchsorted(
                        cumulative_power,
                        cumulative_power.new_tensor(0.95),
                    ).item()
                ) + 1
            else:
                rank95 = 0
            details["metric/prism_module_local_code_participation_rank"] = scalar(
                participation_rank
            )
            details["metric/prism_module_local_code_rank95"] = float(rank95)
            for index, singular_value in enumerate(singular_values):
                details[
                    f"metric/prism_module_local_code_singular_value_{index + 1}"
                ] = scalar(singular_value)
        else:
            details["metric/prism_module_local_code_participation_rank"] = 0.0
            details["metric/prism_module_local_code_rank95"] = 0.0

        if precision.module_local_off_branch_nll_per_cell is not None:
            local_off_branch = (
                precision.module_local_off_branch_nll_per_cell.float().reshape(-1)
            )
            local_on_branch = branch_raw_values.float()
            details["metric/prism_module_local_off_branch_nll"] = scalar(
                (local_off_branch * weight).sum() / weighted_denominator
            )
            details["metric/prism_module_local_branch_nll_gain"] = scalar(
                ((local_off_branch - local_on_branch) * weight).sum()
                / weighted_denominator
            )
        else:
            details["metric/prism_module_local_off_branch_nll"] = 0.0
            details["metric/prism_module_local_branch_nll_gain"] = 0.0
        full_on_nll = getattr(model_out.decoder_out, "nll_per_cell", None)
        if (
            precision.module_local_off_full_nll_per_cell is not None
            and full_on_nll is not None
        ):
            local_off_full = (
                precision.module_local_off_full_nll_per_cell.float().reshape(-1)
            )
            full_on = full_on_nll.float().reshape(-1)
            details["metric/prism_module_local_off_full_nll"] = scalar(
                (local_off_full * weight).sum() / weighted_denominator
            )
            details["metric/prism_module_local_full_nll_gain"] = scalar(
                ((local_off_full - full_on) * weight).sum()
                / weighted_denominator
            )
        else:
            details["metric/prism_module_local_off_full_nll"] = 0.0
            details["metric/prism_module_local_full_nll_gain"] = 0.0
        if precision.module_local_nonlinear_off_branch_nll_per_cell is not None:
            nonlinear_off_branch = (
                precision.module_local_nonlinear_off_branch_nll_per_cell
                .float()
                .reshape(-1)
            )
            nonlinear_on_branch = branch_raw_values.float()
            details[
                "metric/prism_module_local_nonlinear_off_branch_nll"
            ] = scalar(
                (nonlinear_off_branch * weight).sum() / weighted_denominator
            )
            details[
                "metric/prism_module_local_nonlinear_branch_nll_gain"
            ] = scalar(
                ((nonlinear_off_branch - nonlinear_on_branch) * weight).sum()
                / weighted_denominator
            )
        else:
            details[
                "metric/prism_module_local_nonlinear_off_branch_nll"
            ] = 0.0
            details[
                "metric/prism_module_local_nonlinear_branch_nll_gain"
            ] = 0.0
        if (
            precision.module_local_nonlinear_off_full_nll_per_cell is not None
            and full_on_nll is not None
        ):
            nonlinear_off_full = (
                precision.module_local_nonlinear_off_full_nll_per_cell
                .float()
                .reshape(-1)
            )
            full_on = full_on_nll.float().reshape(-1)
            details[
                "metric/prism_module_local_nonlinear_off_full_nll"
            ] = scalar(
                (nonlinear_off_full * weight).sum() / weighted_denominator
            )
            details[
                "metric/prism_module_local_nonlinear_full_nll_gain"
            ] = scalar(
                ((nonlinear_off_full - full_on) * weight).sum()
                / weighted_denominator
            )
        else:
            details[
                "metric/prism_module_local_nonlinear_off_full_nll"
            ] = 0.0
            details[
                "metric/prism_module_local_nonlinear_full_nll_gain"
            ] = 0.0
        details["metric/prism_response_module_rms"] = scalar(
            precision.response_axis_coeff.float().square().mean().sqrt()
        )
        if int(precision.interaction_pair_coeff.shape[1]) > 0:
            details["metric/prism_interaction_module_rms"] = scalar(
                precision.interaction_pair_coeff.float().square().mean().sqrt()
            )
            details["metric/prism_interaction_valid_frac"] = scalar(
                precision.interaction_valid.float().mean()
            )
            details["weight/prism_interaction_ramp"] = scalar(
                head.interaction_scale
            )
            for pair_index, pair_name in enumerate(PATHOLOGY_PAIR_NAMES):
                details[f"metric/prism_interaction_{pair_name}_rms"] = scalar(
                    precision.interaction_pair_coeff[:, pair_index]
                    .float()
                    .square()
                    .mean()
                    .sqrt()
                )
        else:
            details["metric/prism_interaction_module_rms"] = 0.0
            details["metric/prism_interaction_valid_frac"] = 0.0
            details["weight/prism_interaction_ramp"] = 0.0
        details["metric/prism_gene_score_rms"] = scalar(
            precision.gene_score.float().square().mean().sqrt()
        )
        details["metric/prism_code_rms"] = scalar(code.float().square().mean().sqrt())
        details["metric/prism_code_center"] = scalar(code_center.sqrt())
        details["metric/prism_support_contexts"] = scalar(
            precision.support_count.float().mean()
        )
        details["metric/prism_support_reliability"] = scalar(
            precision.support_reliability.float().mean()
        )
        details["metric/prism_pathology_valid_frac"] = scalar(
            precision.pathology_valid.float().mean()
        )
        axis_names = ("thal", "braak", "cerad", "late", "lewy")
        for axis, name in enumerate(axis_names):
            details[f"metric/prism_common_{name}_rms"] = scalar(
                precision.common_axis_coeff[:, axis].float().square().mean().sqrt()
            )
            details[f"metric/prism_response_{name}_rms"] = scalar(
                precision.response_axis_coeff[:, axis].float().square().mean().sqrt()
            )
        details["loss/total"] = scalar(loss_out.total)
        return loss_out

    @staticmethod
    def _add_agp_diagnostics(model_out, loss_out):
        """Add optional AGP-v2 anti-collapse loss and record diagnostics."""

        encoder_out = getattr(model_out, "encoder_out", None)
        weights = getattr(encoder_out, "pooling_attention", None)
        if weights is None:
            return loss_out
        if weights.ndim != 3:
            raise RuntimeError(
                "AGP pooling weights must have shape [batch,tokens,heads], got "
                f"{tuple(weights.shape)}"
            )

        auxiliary = getattr(encoder_out, "pooling_auxiliary", None)
        if auxiliary is not None:
            weighted_total = auxiliary.get("weighted_total")
            if weighted_total is None or weighted_total.ndim != 0:
                raise RuntimeError("AGP auxiliary weighted_total must be scalar")
            loss_out.total = loss_out.total + weighted_total.to(
                dtype=loss_out.total.dtype
            )
            details = loss_out.details
            details["loss/agp_auxiliary"] = float(weighted_total.detach().cpu())
            details["loss/agp_diversity"] = float(
                auxiliary["diversity"].detach().cpu()
            )
            details["loss/agp_entropy_band"] = float(
                auxiliary["entropy_band"].detach().cpu()
            )
            details["loss/agp_query_orthogonality"] = float(
                auxiliary["query_orthogonality"].detach().cpu()
            )
            details["metric/agp_centered_head_overlap"] = float(
                auxiliary["centred_overlap"].detach().cpu()
            )
            details["metric/agp_query_effective_rank"] = float(
                auxiliary["query_effective_rank"].detach().cpu()
            )
            details["metric/agp_score_effective_rank"] = float(
                auxiliary["score_effective_rank"].detach().cpu()
            )
            details["metric/agp_temperature"] = float(
                auxiliary["temperature_mean"].detach().cpu()
            )
            details["metric/agp_mean_residual_gate"] = float(
                auxiliary["mean_residual_gate"].detach().cpu()
            )
        with torch.no_grad():
            w = weights.detach().float().clamp_min(1e-12)
            entropy = -(w * w.log()).sum(dim=1)
            effective_tokens = entropy.exp()
            top1 = w.amax(dim=1)
            topk = min(5, int(w.shape[1]))
            top5_mass = w.topk(topk, dim=1).values.sum(dim=1)

            # Cosine overlap between pooling heads.  Zero means the heads read
            # different tokens; one means they are effectively duplicates.
            by_head = w.transpose(1, 2)
            by_head = by_head / by_head.norm(dim=2, keepdim=True).clamp_min(1e-12)
            overlap = torch.bmm(by_head, by_head.transpose(1, 2))
            n_heads = int(overlap.shape[1])
            if n_heads > 1:
                eye = torch.eye(n_heads, device=overlap.device, dtype=torch.bool)
                head_overlap = overlap[:, ~eye].mean()
            else:
                head_overlap = overlap.new_zeros(())

            details = loss_out.details
            details["metric/agp_effective_tokens"] = float(
                effective_tokens.mean().cpu()
            )
            details["metric/agp_top1_weight"] = float(top1.mean().cpu())
            details["metric/agp_top5_mass"] = float(top5_mass.mean().cpu())
            details["metric/agp_head_overlap"] = float(head_overlap.cpu())
            details["loss/total"] = float(loss_out.total.detach().cpu())
        return loss_out

    def _maybe_add_decoder_reference_anchor(self, batch, model_out, loss_out):
        """Reference-normal decoder anchoring. OFF unless KMLEE_DEC_REF_* lambdas are > 0.

        This is the first precision-medicine structural step:

          reference-clean cells should be reconstructable from base+nuisance
          (celltype baseline + tech + optional sex) with little contribution
          from the state branch.

        It deliberately does NOT alter decoder.forward; it only adds training
        pressure on reference cells, so existing v31a/cons runs remain unchanged
        when the env gates are zero.
        """
        lam_state = float(self.lambda_decoder_ref_state)
        lam_base = float(self.lambda_decoder_ref_base_nll)
        if lam_state <= 0.0 and lam_base <= 0.0:
            return loss_out

        dec_out = getattr(model_out, "decoder_out", None)
        decoder = getattr(self.system, "decoder", None)
        is_reference = batch.get("is_reference", None)
        if dec_out is None or decoder is None or is_reference is None:
            return loss_out

        ref_mask = is_reference.bool().view(-1)
        n_ref = int(ref_mask.sum().item())
        if n_ref < int(self.decoder_ref_min_cells):
            try:
                loss_out.details["metric/decoder_ref_anchor_n_ref"] = float(n_ref)
            except Exception:
                pass
            return loss_out

        total = loss_out.total
        dtype = total.dtype

        if lam_state > 0.0:
            state_ref = dec_out.state_score[ref_mask].float()
            ref_state_l2 = state_ref.pow(2).mean()
            total = total + lam_state * ref_state_l2.to(dtype=dtype)
            try:
                loss_out.details["loss/decoder_ref_state_l2"] = float(ref_state_l2.detach().cpu())
                loss_out.details["weight/lambda_decoder_ref_state"] = float(lam_state)
            except Exception:
                pass

        if lam_base > 0.0:
            y_ord = batch.get("y_ord", None)
            if y_ord is not None:
                base_score = dec_out.base_score + dec_out.tech_score
                if getattr(dec_out, "sex_score", None) is not None:
                    base_score = base_score + dec_out.sex_score
                _, base_probs = decoder._score_to_probs(base_score[ref_mask], dec_out.thresholds)
                _, base_nll_per_cell = decoder._nll_from_probs(base_probs, y_ord[ref_mask].long())
                ref_base_nll = base_nll_per_cell.mean()
                total = total + lam_base * ref_base_nll.to(dtype=dtype)
                try:
                    loss_out.details["loss/decoder_ref_base_nll"] = float(ref_base_nll.detach().cpu())
                    loss_out.details["weight/lambda_decoder_ref_base_nll"] = float(lam_base)
                    if dec_out.nll_per_cell is not None:
                        full_ref_nll = dec_out.nll_per_cell[ref_mask].float().mean()
                        loss_out.details["metric/decoder_ref_full_nll"] = float(full_ref_nll.detach().cpu())
                        loss_out.details["metric/decoder_ref_base_minus_full_nll"] = float(
                            (ref_base_nll.detach() - full_ref_nll.detach()).cpu()
                        )
                except Exception:
                    pass

        loss_out.total = total
        try:
            loss_out.details["metric/decoder_ref_anchor_n_ref"] = float(n_ref)
            loss_out.details["loss/total"] = float(loss_out.total.detach().cpu())
        except Exception:
            pass
        return loss_out

    def _maybe_add_joint_generator_count(
        self,
        batch,
        model_out,
        loss_out,
        *,
        optimize: bool,
    ):
        """Add the online constrained-L0 objective for generator count.

        The candidate is the exact hard gate already used by the ordinary
        forward.  Its paired reference is an all-on counterfactual evaluated
        with the same latent, minibatch and current weights.  A second
        generator-isolated comparison removes Direct/PRISM/mixer bypasses.
        Thus ``K`` is learned by gradients during this run; no post-hoc K grid
        or R10 frozen signature participates.
        """

        system = getattr(self.system, "module", self.system)
        gate = getattr(system, "generator_count_gate", None)
        objective = getattr(system, "generator_count_objective", None)
        cfg = getattr(system, "generator_count_config", None)
        decoder_out = getattr(model_out, "decoder_out", None)
        if (
            gate is None
            or objective is None
            or cfg is None
            or not bool(getattr(cfg, "enabled", False))
            or str(getattr(cfg, "mode", "gate_only")) != "joint"
            or decoder_out is None
        ):
            return loss_out

        hard = getattr(system, "_last_generator_hard_mask", None)
        forward_gate = getattr(system, "_last_generator_gate", None)
        expected_active = getattr(
            system, "_last_generator_expected_active", None
        )
        temperature = float(
            getattr(system, "_last_generator_temperature", math.nan)
        )
        if hard is None or forward_gate is None or expected_active is None:
            raise RuntimeError("joint generator gate state is missing after forward")
        if not bool(((hard == 0) | (hard == 1)).all()):
            raise RuntimeError("joint generator hard mask lost exact binarity")

        hard_active = int(hard.detach().sum().item())
        probability = gate.keep_probability(temperature=temperature)
        details = loss_out.details
        details["metric/generator_candidate_count"] = float(
            gate.num_generators
        )
        details["metric/generator_hard_active"] = float(hard_active)
        details["metric/generator_expected_active"] = float(
            expected_active.detach().float().cpu()
        )
        details["metric/generator_raw_expected_active"] = float(
            probability.detach().float().sum().cpu()
        )
        details["metric/generator_temperature"] = temperature
        gate_mode = str(
            getattr(system, "_last_generator_gate_mode", "protected_hard")
        )
        details["metric/generator_mode_shadow"] = float(
            gate_mode == "shadow_all_on"
        )
        details["metric/generator_mode_soft"] = float(
            gate_mode == "soft_adaptation"
        )
        details["metric/generator_mode_hard"] = float(
            gate_mode in {"protected_hard", "rescue_deterministic_candidate"}
        )
        details["metric/generator_uncertain_fraction"] = float(
            (((probability > 0.1) & (probability < 0.9)).float().mean())
            .detach()
            .cpu()
        )

        epoch_state = getattr(system, "generator_count_epoch_state", None)
        epoch = 1 if epoch_state is None else int(epoch_state.item())
        in_warmup = epoch < int(cfg.start_epoch)
        details["metric/generator_warmup_all_on"] = float(in_warmup)

        decoder = system.decoder
        y_ord = batch["y_ord"].long()
        celltype_id = batch["celltype_id"].long()
        tech_id = batch.get("tech_id", batch.get("batch_id")).long()
        sex_id = batch.get("sex_id")
        precision = getattr(model_out, "precision_out", None)
        precision_score = (
            None
            if precision is None or precision.gene_score is None
            else precision.gene_score.detach()
        )
        ones = torch.ones(
            int(decoder.n_generators),
            dtype=forward_gate.dtype,
            device=forward_gate.device,
        )

        stash = getattr(decoder, "_stash_diag", True)
        decoder._stash_diag = False
        try:
            with torch.no_grad():
                all_on = decoder(
                    z_perp=model_out.z_perp.detach(),
                    celltype_id=celltype_id,
                    tech_id=tech_id,
                    y_ord=y_ord,
                    sex_id=sex_id,
                    precision_score=precision_score,
                    generator_gate=ones,
                )
                if all_on.nll_per_cell is None:
                    raise RuntimeError("all-on counterfactual did not return NLL")
                baseline_full_nll = all_on.nll_per_cell.float().mean()

            candidate_decoder_out = decoder_out
            if gate_mode == "shadow_all_on":
                candidate_decoder_out = decoder(
                    z_perp=model_out.z_perp,
                    celltype_id=celltype_id,
                    tech_id=tech_id,
                    y_ord=y_ord,
                    sex_id=sex_id,
                    precision_score=precision_score,
                    generator_gate=forward_gate,
                )
            candidate_full_nll = candidate_decoder_out.nll_per_cell
            if candidate_full_nll is None:
                raise RuntimeError("gated ordinary forward did not return NLL")
            candidate_full_nll = candidate_full_nll.float().mean()

            base_score = decoder.celltype_baseline(celltype_id)
            tech_score = decoder.tech_baseline(tech_id)

            def isolated_nll(
                selected_gate: torch.Tensor,
                *,
                detach_inputs: bool,
            ) -> torch.Tensor:
                z = model_out.z_perp.detach() if detach_inputs else model_out.z_perp
                base = base_score.detach() if detach_inputs else base_score
                tech = tech_score.detach() if detach_inputs else tech_score
                state = decoder.state_score_from_z_and_base(
                    z,
                    base,
                    celltype_id=celltype_id,
                    generator_gate=selected_gate,
                    include_direct=False,
                    include_pathology=False,
                    include_interactions=False,
                )
                score = base + tech + state
                if decoder.sex_baseline is not None and sex_id is not None:
                    sex_score = decoder.sex_baseline(decoder._sex_index(sex_id))
                    score = score + (
                        sex_score.detach() if detach_inputs else sex_score
                    )
                _, probabilities = decoder._score_to_probs(
                    score, decoder._compute_thresholds()
                )
                _, nll_per_cell = decoder._nll_from_probs(probabilities, y_ord)
                return nll_per_cell.float().mean()

            if in_warmup:
                # Log paired all-on metrics from epoch 1 onward so composite
                # checkpoint statistics never treat the first gated epoch as
                # if two missing observations had value zero.  No gate/count
                # objective is optimized during this exact-all-on warm-up.
                with torch.no_grad():
                    candidate_isolated_nll = isolated_nll(
                        forward_gate, detach_inputs=True
                    )
            else:
                candidate_isolated_nll = isolated_nll(
                    forward_gate, detach_inputs=False
                )
            with torch.no_grad():
                baseline_isolated_nll = isolated_nll(
                    ones, detach_inputs=True
                )
        finally:
            decoder._stash_diag = stash

        full_margin = float(
            baseline_full_nll.detach().cpu()
        ) * float(cfg.joint_full_nll_relative_margin)
        isolated_margin = float(
            baseline_isolated_nll.detach().cpu()
        ) * float(cfg.joint_isolated_nll_relative_margin)
        scale_floor = float(cfg.joint_constraint_scale_floor)
        constraints = ConstraintBatch(
            baseline={
                "full_nll": baseline_full_nll,
                "isolated_nll": baseline_isolated_nll,
            },
            candidate={
                "full_nll": candidate_full_nll,
                "isolated_nll": candidate_isolated_nll,
            },
            margin={
                "full_nll": full_margin,
                "isolated_nll": isolated_margin,
            },
            scale={
                "full_nll": max(full_margin, scale_floor),
                "isolated_nll": max(isolated_margin, scale_floor),
            },
        )
        count_out = objective(
            expected_active_count=expected_active,
            num_generators=int(decoder.n_generators),
            constraints=constraints,
        )

        if optimize and not in_warmup:
            weighted = float(cfg.joint_objective_weight) * count_out.loss
            loss_out.total = loss_out.total + weighted.to(loss_out.total.dtype)
            # Every rank contributes its paired violation.  Explicit averaging
            # keeps the persistent dual buffers identical because DDP buffer
            # broadcasting is intentionally disabled in this project.
            dual_violation = torch.stack(
                [
                    count_out.violations["full_nll"].detach().float(),
                    count_out.violations["isolated_nll"].detach().float(),
                ]
            )
            if dist.is_available() and dist.is_initialized():
                dist.all_reduce(dual_violation, op=dist.ReduceOp.SUM)
                dual_violation.div_(float(dist.get_world_size()))
            pending_sum = getattr(self, "_joint_dual_violation_sum", None)
            if pending_sum is None:
                self._joint_dual_violation_sum = dual_violation.clone()
                self._joint_dual_violation_count = 1
            else:
                self._joint_dual_violation_sum = pending_sum + dual_violation
                self._joint_dual_violation_count = int(
                    getattr(self, "_joint_dual_violation_count", 0)
                ) + 1
            details["loss/generator_count"] = float(weighted.detach().cpu())
            details["loss/total"] = float(loss_out.total.detach().cpu())
        else:
            details["loss/generator_count"] = 0.0

        details["metric/generator_baseline_full_nll"] = float(
            baseline_full_nll.detach().cpu()
        )
        details["metric/generator_candidate_full_nll"] = float(
            candidate_full_nll.detach().cpu()
        )
        details["metric/generator_baseline_isolated_nll"] = float(
            baseline_isolated_nll.detach().cpu()
        )
        details["metric/generator_candidate_isolated_nll"] = float(
            candidate_isolated_nll.detach().cpu()
        )
        for name in ("full_nll", "isolated_nll"):
            details[f"metric/generator_violation_{name}"] = float(
                count_out.violations[name].detach().cpu()
            )
        return loss_out

    def train_step(self, batch: Dict[str, torch.Tensor]) -> StepOutput:
        """
        Single optimizer update on one batch.

        This method intentionally ignores `grad_accum_steps` and performs one
        full update. Gradient accumulation is handled in `train_epoch`, where the
        trainer can correctly manage optimizer-step boundaries.
        """
        self.system.train()
        batch = self._move_batch_to_device(batch)

        self.optimizer.zero_grad(set_to_none=True)

        with self._autocast_context():
            model_out = self.system(
                batch,
                sample_latent=True,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )
            loss_out = self._compute_loss(batch, model_out)
            loss_out = self._add_agp_diagnostics(model_out, loss_out)
            loss_out = self._maybe_add_precision_medicine(batch, model_out, loss_out)
            loss_out = self._maybe_add_joint_generator_count(
                batch, model_out, loss_out, optimize=True
            )
            loss_out = self._maybe_add_path_consistency(batch, model_out, loss_out)
            loss_out = self._maybe_add_cp_bam(batch, model_out, loss_out)
            loss_out = self._maybe_add_decoder_reference_anchor(batch, model_out, loss_out)

        backward_ok = self._backward(loss_out.total, loss_out=loss_out)

        if backward_ok:
            grad_norm = self._finish_optimizer_step()
        else:
            self.optimizer.zero_grad(set_to_none=True)
            self._clear_joint_generator_dual_pending()
            grad_norm = None
       

        return StepOutput(
            loss=loss_out,
            grad_norm=grad_norm,
            lr=self.current_lr,
            batch_size=int(batch["y_ord"].shape[0]),
        )

    @torch.no_grad()
    def eval_step(self, batch: Dict[str, torch.Tensor]) -> StepOutput:
        self.system.eval()
        batch = self._move_batch_to_device(batch)

        with self._autocast_context():
            model_out = self.system(
                batch,
                sample_latent=self.eval_sample_latent,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )
            loss_out = self._compute_loss(batch, model_out)
            loss_out = self._add_agp_diagnostics(model_out, loss_out)
            loss_out = self._maybe_add_precision_medicine(batch, model_out, loss_out)
            loss_out = self._maybe_add_joint_generator_count(
                batch, model_out, loss_out, optimize=False
            )

        return StepOutput(
            loss=loss_out,
            grad_norm=None,
            lr=self.current_lr,
            batch_size=int(batch["y_ord"].shape[0]),
        )

    def curriculum_state_dict(self) -> dict:
        """Persist deterministic curriculum position and live branch scales."""

        system = getattr(self.system, "module", self.system)
        head = getattr(system, "precision_head", None)
        generator_epoch = getattr(system, "generator_count_epoch_state", None)
        state = {
            "schema_version": "kmlee_bam.integrated_curriculum.v1",
            "epoch": int(getattr(self, "_prism_last_epoch", 0) or 0),
            "precision_explicit_scale": (
                None
                if head is None or not hasattr(head, "explicit_scale")
                else float(head.explicit_scale.detach().cpu())
            ),
            "precision_interaction_scale": (
                None
                if head is None or not hasattr(head, "interaction_scale")
                else float(head.interaction_scale.detach().cpu())
            ),
            "generator_epoch": (
                None
                if generator_epoch is None
                else int(generator_epoch.detach().cpu().item())
            ),
        }
        controller = getattr(self, "integrated_phase_controller", None)
        if controller is not None:
            state["integrated_phase_controller"] = controller.state_dict()
        updater = getattr(self, "module_rescue_updater", None)
        if updater is not None and hasattr(updater, "state_dict"):
            state["module_rescue_state"] = updater.state_dict()
        return state

    @torch.no_grad()
    def curriculum_load_state_dict(self, state: dict) -> None:
        """Restore a checkpointed integrated curriculum state."""

        if not state:
            return
        epoch = int(state.get("epoch", 0) or 0)
        if epoch > 0:
            self._prism_last_epoch = epoch
        system = getattr(self.system, "module", self.system)
        head = getattr(system, "precision_head", None)
        if head is not None:
            explicit = state.get("precision_explicit_scale")
            if explicit is not None and hasattr(head, "explicit_scale"):
                head.explicit_scale.fill_(float(explicit))
            interaction = state.get("precision_interaction_scale")
            if interaction is not None and hasattr(head, "interaction_scale"):
                head.interaction_scale.fill_(float(interaction))
        generator_epoch = state.get("generator_epoch")
        epoch_state = getattr(system, "generator_count_epoch_state", None)
        if generator_epoch is not None and epoch_state is not None:
            epoch_state.fill_(int(generator_epoch))
        controller_state = state.get("integrated_phase_controller")
        controller = getattr(self, "integrated_phase_controller", None)
        if controller_state is not None:
            if controller is None:
                raise RuntimeError(
                    "checkpoint contains integrated Phase curriculum state but "
                    "the current trainer has no controller"
                )
            controller.load_state_dict(controller_state)
        module_rescue_state = state.get("module_rescue_state")
        if module_rescue_state is not None:
            updater = getattr(self, "module_rescue_updater", None)
            if updater is None:
                self._pending_module_rescue_state = module_rescue_state
            else:
                updater.load_state_dict(module_rescue_state)

    def train_epoch(
        self,
        loader: Iterable[Dict[str, torch.Tensor]],
        *,
        log_every: Optional[int] = None,
        step_callback: Optional[Any] = None,
        epoch_index: Optional[int] = None,
        max_steps: Optional[int] = None,
    ) -> EpochOutput:
        self.system.train()
        controller = getattr(self, "integrated_phase_controller", None)
        if controller is not None and epoch_index is not None:
            phase_state = controller.apply_epoch(self, int(epoch_index))
            self._integrated_phase_live_state = dict(phase_state)
        if epoch_index is not None:
            joint_system = getattr(self.system, "module", self.system)
            pathology_decoder = getattr(joint_system, "decoder", None)
            pathology_rank_gate = getattr(
                pathology_decoder, "pathology_rank_gate", None
            )
            if pathology_rank_gate is not None:
                pathology_rank_gate.set_epoch(int(epoch_index))
            if bool(getattr(joint_system, "generator_count_joint_training", False)):
                epoch_state = getattr(
                    joint_system, "generator_count_epoch_state", None
                )
                if epoch_state is None:
                    raise RuntimeError(
                        "joint generator count lacks persistent epoch state"
                    )
                epoch_state.fill_(int(epoch_index))
        precision_head = getattr(self.system, "precision_head", None)
        if precision_head is not None and epoch_index is not None:
            if hasattr(precision_head, "set_curriculum_epoch"):
                precision_head.set_curriculum_epoch(int(epoch_index))
            if hasattr(precision_head, "set_interaction_epoch"):
                precision_head.set_interaction_epoch(int(epoch_index))
            if hasattr(precision_head, "set_module_local_epoch"):
                precision_head.set_module_local_epoch(int(epoch_index))
            if hasattr(precision_head, "set_module_local_nonlinear_epoch"):
                precision_head.set_module_local_nonlinear_epoch(int(epoch_index))
        meters = RunningAverages()
        epoch_start = time.perf_counter()
        step_debug = os.environ.get("KMLEE_TRAIN_STEP_DEBUG", "0") == "1"
        debug_rank = int(os.environ.get("RANK", "0"))

        def _step_log(message: str) -> None:
            if step_debug:
                print(f"[step-debug rank={debug_rank}] {message}", flush=True)

        self.optimizer.zero_grad(set_to_none=True)

        try:
            n_total_steps = len(loader)
        except TypeError:
            n_total_steps = None

        allow_no_sync = (
            self.grad_accum_steps > 1
            and n_total_steps is not None
        )

        pending_micro_steps = 0

        for step_idx, batch in enumerate(loader, start=1):
            _step_log(f"step={step_idx} batch:loaded")
            batch = self._move_batch_to_device(batch)
            _step_log(f"step={step_idx} batch:on_device")

            is_boundary = (step_idx % self.grad_accum_steps == 0)
            is_known_last = (n_total_steps is not None and step_idx == n_total_steps)
            should_step = is_boundary or is_known_last

            sync_ctx = self._ddp_no_sync_context(
                use_no_sync=(allow_no_sync and not should_step)
            )

            with sync_ctx:
                with self._autocast_context():
                    _step_log(f"step={step_idx} forward:start boundary={should_step}")
                    model_out = self.system(
                        batch,
                        sample_latent=True,
                        return_all_hidden_states=False,
                        return_attn_diagnostics=False,
                    )
                    _step_log(f"step={step_idx} forward:done")
                    _step_log(f"step={step_idx} loss:start")
                    loss_out = self._compute_loss(batch, model_out)
                    loss_out = self._add_agp_diagnostics(model_out, loss_out)
                    loss_out = self._maybe_add_precision_medicine(batch, model_out, loss_out)
                    loss_out = self._maybe_add_joint_generator_count(
                        batch, model_out, loss_out, optimize=True
                    )
                    loss_out = self._maybe_add_path_consistency(batch, model_out, loss_out)
                    loss_out = self._maybe_add_cp_bam(batch, model_out, loss_out)
                    loss_out = self._maybe_add_decoder_reference_anchor(batch, model_out, loss_out)
                    _step_log(f"step={step_idx} loss:done")

                if (
                    self.debug_z_sensitivity_every is not None
                    and self.debug_z_sensitivity_every > 0
                    and step_idx % self.debug_z_sensitivity_every == 0
                ):
                    self._debug_z_sensitivity(batch, model_out, step_idx)

                loss_for_backward = loss_out.total / float(self.grad_accum_steps)
                _step_log(f"step={step_idx} backward:start")
                backward_ok = self._backward(
                    loss_for_backward,
                    loss_out=loss_out,
                    step_idx=step_idx,
                )
                _step_log(f"step={step_idx} backward:done ok={backward_ok}")

                if not backward_ok:
                    self.optimizer.zero_grad(set_to_none=True)
                    self._clear_joint_generator_dual_pending()
                    pending_micro_steps = 0
                    continue

            pending_micro_steps += 1
            grad_norm = None

            if should_step:
                _step_log(f"step={step_idx} optimizer:start")
                grad_norm = self._finish_optimizer_step()
                _step_log(f"step={step_idx} optimizer:done grad_norm={grad_norm}")
                pending_micro_steps = 0

            step_out = StepOutput(
                loss=loss_out,
                grad_norm=grad_norm,
                lr=self.current_lr,
                batch_size=int(batch["y_ord"].shape[0]),
            )
            meters.update_from_step(step_out)

            if step_callback is not None:
                step_callback(
                    step_idx=step_idx,
                    n_total_steps=n_total_steps,
                    step_out=step_out,
                    meters=meters,
                    epoch_index=epoch_index,
                )

            if log_every is not None and step_idx % log_every == 0:
                print(self._format_step_log(prefix="train", step_idx=step_idx, meters=meters))

            if max_steps is not None and int(max_steps) > 0 and step_idx >= int(max_steps):
                break

        if pending_micro_steps > 0:
            # This path matters when `loader` does not expose __len__ and the
            # final partial accumulation chunk would otherwise never be stepped.
            self._finish_optimizer_step()

        if self.scheduler is not None and self.scheduler_step_on == "epoch":
            self.scheduler.step()

        return EpochOutput(
            metrics=meters.compute(),
            n_examples=meters.n_examples,
            elapsed_sec=time.perf_counter() - epoch_start,
        )

    @torch.no_grad()
    def evaluate_epoch(
        self,
        loader: Iterable[Dict[str, torch.Tensor]],
        *,
        log_every: Optional[int] = None,
        max_steps: Optional[int] = None,
    ) -> EpochOutput:
        meters = RunningAverages()
        epoch_start = time.perf_counter()

        for step_idx, batch in enumerate(loader, start=1):
            step_out = self.eval_step(batch)
            meters.update_from_step(step_out)

            if log_every is not None and step_idx % log_every == 0:
                print(self._format_step_log(prefix="eval", step_idx=step_idx, meters=meters))

            if max_steps is not None and int(max_steps) > 0 and step_idx >= int(max_steps):
                break

        return EpochOutput(
            metrics=meters.compute(),
            n_examples=meters.n_examples,
            elapsed_sec=time.perf_counter() - epoch_start,
        )

    def fit(
        self,
        train_loader: Iterable[Dict[str, torch.Tensor]],
        *,
        epochs: int,
        val_loader: Optional[Iterable[Dict[str, torch.Tensor]]] = None,
        log_every: Optional[int] = None,
    ) -> list[Dict[str, Dict[str, float]]]:
        if epochs <= 0:
            raise ValueError("epochs must be positive.")

        history: list[Dict[str, Dict[str, float]]] = []
        for epoch in range(1, epochs + 1):
            train_out = self.train_epoch(
                train_loader, log_every=log_every, epoch_index=epoch
            )
            epoch_record: Dict[str, Dict[str, float]] = {"train": train_out.metrics}

            if val_loader is not None:
                val_out = self.evaluate_epoch(val_loader, log_every=log_every)
                epoch_record["val"] = val_out.metrics
                print(
                    f"[epoch {epoch:03d}] "
                    f"train_total={train_out.metrics.get('loss/total', float('nan')):.6f}  "
                    f"val_total={val_out.metrics.get('loss/total', float('nan')):.6f}"
                )
            else:
                print(
                    f"[epoch {epoch:03d}] "
                    f"train_total={train_out.metrics.get('loss/total', float('nan')):.6f}"
                )

            history.append(epoch_record)

        return history

    @torch.no_grad()
    def collect_latents(
        self,
        loader: Iterable[Dict[str, torch.Tensor]],
    ) -> Dict[str, torch.Tensor]:
        """
        Collect posterior parameters, residualized latents, metadata, and
        optional cell-level BAM uncertainty for downstream analysis.
        """
        self.system.eval()

        mu_q_list = []
        logvar_q_list = []
        z_s_list = []
        z_perp_list = []
        celltype_list = []
        tech_list = []
        row_index_list = []
        cell_uncertainty_list = []

        for batch in loader:
            batch = self._move_batch_to_device(batch)
            model_out = self.system(
                batch,
                sample_latent=self.eval_sample_latent,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )

            mu_q_list.append(model_out.encoder_out.mu_q.detach().cpu())
            logvar_q_list.append(model_out.encoder_out.logvar_q.detach().cpu())
            z_s_list.append(model_out.encoder_out.z_s.detach().cpu())
            z_perp_list.append(model_out.z_perp.detach().cpu())
            celltype_list.append(batch["celltype_id"].detach().cpu())
            tech_list.append(_resolve_tech_id(batch).detach().cpu())

            if "row_index" in batch:
                row_index_list.append(batch["row_index"].detach().cpu())

            if model_out.encoder_out.cell_uncertainty is not None:
                cell_uncertainty_list.append(
                    model_out.encoder_out.cell_uncertainty.detach().cpu()
                )

        out = {
            "mu_q": torch.cat(mu_q_list, dim=0),
            "logvar_q": torch.cat(logvar_q_list, dim=0),
            "z_s": torch.cat(z_s_list, dim=0),
            "z_perp": torch.cat(z_perp_list, dim=0),
            "celltype_id": torch.cat(celltype_list, dim=0),
            "tech_id": torch.cat(tech_list, dim=0),
        }
        if row_index_list:
            out["row_index"] = torch.cat(row_index_list, dim=0)
        if cell_uncertainty_list:
            out["cell_uncertainty"] = torch.cat(cell_uncertainty_list, dim=0)
        return out

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------
    def _compute_loss(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
    ) -> LossOutput:
        celltype_id = batch["celltype_id"]
        tech_id = _resolve_tech_id(batch)
        cls_targets = batch.get("cls_target")

        tech_group_weights = None
        if self.tech_weight_mode == "empirical_batch":
            tech_group_weights = empirical_group_weights(tech_id, n_groups=self.system.decoder.n_tech)
            tech_group_weights = tech_group_weights.to(device=tech_id.device, dtype=model_out.decoder_out.score.dtype)

        whitening_target = None
        if (
            model_out.z_perp is not None
            and float(getattr(self.criterion, "lambda_white", 0.0)) != 0.0
        ):
            nuisance_projector = getattr(
                self.system,
                "nuisance_projector",
                None,
            )
            if (
                nuisance_projector is not None
                and bool(
                    getattr(
                        getattr(nuisance_projector, "cfg", None),
                        "enabled",
                        False,
                    )
                )
                and bool(
                    getattr(
                        getattr(nuisance_projector, "cfg", None),
                        "projection_aware_whitening",
                        False,
                    )
                )
            ):
                whitening_target = nuisance_projector.whitening_target(
                    device=model_out.z_perp.device,
                    dtype=torch.float32,
                )

        return self.criterion(
            prior=self.system.prior,
            decoder=self.system.decoder,
            encoder_out=model_out.encoder_out,
            decoder_out=model_out.decoder_out,
            celltype_id=celltype_id,
            tech_id=tech_id,
            cls_logits=model_out.cls_logits,
            cls_targets=cls_targets,
            tech_group_weights=tech_group_weights,
            n_valid_tokens_for_uncertainty=self.system.n_effective_encoder_tokens_for_uncertainty,
            z_perp=model_out.z_perp,
            is_reference=batch.get("is_reference", None),
            # Optional per-SEX alignment label. The DLPFC+MTG dataset does NOT
            # currently emit a sex id, so this is None and the sex term stays off
            # (lambda_align_sex defaults to 0). If a per-cell `sex_id` is added to
            # the batch later, it flows through here with no further plumbing.
            sex_id=batch.get("sex_id", None),
            whitening_target=whitening_target,
            # Stage D conduit: a top-of-MRO trainer (V8) may stash per-(cell,gene)
            # reconstruction-reliability weights just before super()._compute_loss.
            # Absent ⇒ None ⇒ reconstruction is byte-identical.
            reconstruction_weight_per_gene=getattr(self, "_pending_rec_weight_per_gene", None),
        )

    def _move_batch_to_device(self, batch: Dict[str, Any]) -> Dict[str, torch.Tensor]:
        moved: Dict[str, torch.Tensor] = {}
        for key, value in batch.items():
            if torch.is_tensor(value):
                moved[key] = value.to(self.device, non_blocking=True)
            else:
                moved[key] = value
        return moved

    def _autocast_context(self):
        if self.device.type == "cuda" and self.amp_enabled:
            return torch.autocast(
                device_type="cuda", 
                dtype=self.amp_dtype, 
                enabled=True,
                )
        return nullcontext()

    def _all_ranks_finite(self, local_finite: bool) -> bool:
        """
        Collective: returns True iff ALL ranks report finite loss.

        DDP requires all ranks to make the same skip decision; otherwise the
        skipping rank does not enter the gradient AllReduce while other ranks
        wait, causing NCCL watchdog timeout.
        """
        flag = torch.tensor(
            [1.0 if local_finite else 0.0],
            device=self.device,
            dtype=torch.float32,
        )
        if dist.is_available() and dist.is_initialized():
            dist.all_reduce(flag, op=dist.ReduceOp.MIN)
        return float(flag.item()) >= 0.5

    def _dump_non_finite_components(
        self,
        loss_out: Optional[LossOutput],
        *,
        step_idx: Optional[int] = None,
    ) -> None:
        rank = int(os.environ.get("RANK", "0"))
        label = f"step={step_idx}" if step_idx is not None else ""
        if loss_out is None or not getattr(loss_out, "details", None):
            print(f"[non-finite components | rank={rank} {label}]: <no details available>", flush=True)
            return
        bad: Dict[str, float] = {}
        for name, value in loss_out.details.items():
            try:
                v = float(value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(v):
                bad[name] = v
        if not bad:
            total_v = loss_out.details.get("loss/total")
            print(
                f"[non-finite components | rank={rank} {label}]: "
                f"<total non-finite={total_v} but no individual component flagged>",
                flush=True,
            )
        else:
            print(f"[non-finite components | rank={rank} {label}]: {bad}", flush=True)

    def _backward(
        self,
        loss: torch.Tensor,
        *,
        loss_out: Optional[LossOutput] = None,
        step_idx: Optional[int] = None,
    ) -> bool:
        """
        Backward with DDP-safe collective skip.

        If ANY rank reports non-finite loss, ALL ranks skip backward + AllReduce
        together. This prevents the NCCL watchdog timeout we observed in the
        L2.5 run where rank-local skip desynchronized collective ops.

        On the rank that is locally non-finite, we also dump which loss
        component(s) produced the NaN/Inf for downstream diagnosis.
        """
        local_finite = bool(torch.isfinite(loss).item())
        if not local_finite:
            print(
                f"[skip] non-finite loss before backward: {float(loss.detach().cpu())}",
                flush=True,
            )
            self._dump_non_finite_components(loss_out, step_idx=step_idx)

        all_finite = self._all_ranks_finite(local_finite)
        if not all_finite:
            if local_finite:
                # This rank's loss was finite but another rank reported NaN/Inf;
                # the operator looking at rank 0 logs needs a breadcrumb that
                # tells them *why* the optimizer step was skipped, since the
                # component dump only shows up on the offending rank.
                rank = int(os.environ.get("RANK", "0"))
                label = f"step={step_idx}" if step_idx is not None else ""
                print(
                    f"[skip] global non-finite loss; at least one other rank "
                    f"reported NaN/Inf (rank={rank} {label})",
                    flush=True,
                )
            return False

        if self.scaler.is_enabled():
            self.scaler.scale(loss).backward()
        else:
            loss.backward()
        return True
    
    def _clear_joint_generator_dual_pending(self) -> None:
        self._joint_dual_violation_sum = None
        self._joint_dual_violation_count = 0

    @torch.no_grad()
    def _commit_joint_generator_dual_pending(self) -> None:
        pending = getattr(self, "_joint_dual_violation_sum", None)
        count = int(getattr(self, "_joint_dual_violation_count", 0))
        if pending is None or count <= 0:
            return
        system = getattr(self.system, "module", self.system)
        objective = getattr(system, "generator_count_objective", None)
        if objective is None:
            raise RuntimeError("pending generator dual update lacks objective")
        mean = pending / float(count)
        objective.update_dual(
            {"full_nll": mean[0], "isolated_nll": mean[1]}
        )
        self._clear_joint_generator_dual_pending()

    def _finish_optimizer_step(
        self,
        *,
        advance_scheduler: bool = True,
    ) -> Optional[float]:
        opt_debug = os.environ.get("KMLEE_OPT_STEP_DEBUG", "0") == "1"
        debug_rank = int(os.environ.get("RANK", "0"))

        def _opt_log(message: str) -> None:
            if opt_debug:
                print(f"[opt-debug rank={debug_rank}] {message}", flush=True)

        self._last_step_was_skipped = False
        grad_norm: Optional[float] = None

        # fp16 GradScaler를 쓰는 경우, clipping 전에 반드시 unscale.
        if self.scaler.is_enabled():
            _opt_log("unscale:start")
            self.scaler.unscale_(self.optimizer)
            _opt_log("unscale:done")

        if self.grad_clip_norm is not None:
            _opt_log("clip:start")
            grad_norm_tensor = self._manual_clip_grad_norm(float(self.grad_clip_norm))
            _opt_log("clip:done")
            grad_norm = float(grad_norm_tensor.detach().cpu())

            if not torch.isfinite(grad_norm_tensor):
                print(f"[skip] non-finite grad norm: {grad_norm}")
                if os.environ.get("KMLEE_NONFINITE_GRAD_DEBUG", "0") == "1":
                    self._dump_non_finite_gradients()
                self.optimizer.zero_grad(set_to_none=True)

                # fp16 scaler 상태는 업데이트해줘야 다음 scale 조정이 이어진다.
                if self.scaler.is_enabled():
                    self.scaler.update()

                self._last_step_was_skipped = True
                self._clear_joint_generator_dual_pending()
                return grad_norm

        if self.scaler.is_enabled():
            old_scale = self.scaler.get_scale()
            _opt_log("scaler_step:start")
            self.scaler.step(self.optimizer)
            _opt_log("scaler_step:done")
            _opt_log("scaler_update:start")
            self.scaler.update()
            _opt_log("scaler_update:done")
            new_scale = self.scaler.get_scale()

            if new_scale < old_scale:
                print(f"[AMP] overflow detected: scale {old_scale} -> {new_scale}")
                self._last_step_was_skipped = True
        else:
            _opt_log("optimizer_step:start")
            self.optimizer.step()
            _opt_log("optimizer_step:done")

        _opt_log("zero_grad:start")
        self.optimizer.zero_grad(set_to_none=True)
        _opt_log("zero_grad:done")

        if self._last_step_was_skipped:
            self._clear_joint_generator_dual_pending()
        else:
            self._commit_joint_generator_dual_pending()

        # optimizer step을 실제로 한 경우에만 scheduler step.
        if (
            bool(advance_scheduler)
            and self.scheduler is not None
            and self.scheduler_step_on == "step"
        ):
            _opt_log("scheduler_step:start")
            self.scheduler.step()
            _opt_log("scheduler_step:done")

        return grad_norm

    def _manual_clip_grad_norm(self, max_norm: float) -> torch.Tensor:
        """Local gradient clipping that avoids torch's multi-tensor helper path."""
        grads = []
        for group in self.optimizer.param_groups:
            for p in group.get("params", []):
                if p.grad is not None:
                    grads.append(p.grad.detach())

        if not grads:
            return torch.zeros((), dtype=torch.float32, device=self.device)

        total = torch.zeros((), dtype=torch.float32, device=grads[0].device)
        for grad in grads:
            total = total + grad.float().pow(2).sum()
        total_norm = total.sqrt()

        clip_coef = float(max_norm) / (total_norm + 1e-6)
        if bool(torch.isfinite(total_norm).item()) and float(clip_coef.detach().cpu()) < 1.0:
            coef = clip_coef.to(device=grads[0].device, dtype=grads[0].dtype)
            for grad in grads:
                grad.mul_(coef.to(device=grad.device, dtype=grad.dtype))
        return total_norm

    def _dump_non_finite_gradients(self) -> None:
        """Print parameter-level NaN/Inf evidence for an explicitly debugged run."""

        rank = int(os.environ.get("RANK", "0"))
        names: Dict[int, str] = {}
        for prefix, module in (
            ("system", self.system),
            ("criterion", self.criterion),
        ):
            for name, parameter in module.named_parameters():
                names.setdefault(id(parameter), f"{prefix}.{name}")

        records = []
        for group_index, group in enumerate(self.optimizer.param_groups):
            for parameter_index, parameter in enumerate(group.get("params", [])):
                grad = parameter.grad
                if grad is None:
                    continue
                finite = torch.isfinite(grad)
                if bool(finite.all().item()):
                    continue
                detached = grad.detach()
                finite_values = detached[finite]
                records.append(
                    {
                        "name": names.get(
                            id(parameter),
                            f"optimizer_group{group_index}.parameter{parameter_index}",
                        ),
                        "shape": tuple(parameter.shape),
                        "dtype": str(detached.dtype),
                        "nan": int(torch.isnan(detached).sum().item()),
                        "posinf": int(torch.isposinf(detached).sum().item()),
                        "neginf": int(torch.isneginf(detached).sum().item()),
                        "finite_abs_max": (
                            float(finite_values.abs().max().item())
                            if int(finite_values.numel()) > 0
                            else None
                        ),
                    }
                )
        print(
            f"[non-finite gradients | rank={rank}] count={len(records)} "
            f"records={records[:64]}",
            flush=True,
        )


    def _ddp_no_sync_context(self, *, use_no_sync: bool):
        if not use_no_sync:
            return nullcontext()

        ddp_model = getattr(self.system, "_ddp_model", None)
        if ddp_model is not None and hasattr(ddp_model, "no_sync"):
            return ddp_model.no_sync()

        no_sync = getattr(self.system, "no_sync", None)
        if callable(no_sync):
            return no_sync()

        return nullcontext()
    
   
    @torch.no_grad()
    def _debug_z_sensitivity(
        self,
        batch: Dict[str, torch.Tensor],
        model_out: ModelForwardOutput,
        step_idx: int,
    ) -> None:
        was_training = self.system.training
        self.system.eval()

        # deterministic posterior mean path
        det_out = self.system(
            batch,
            sample_latent=False,
            return_all_hidden_states=False,
            return_attn_diagnostics=False,
        )

        z = det_out.z_perp.detach()
        z_perm = z[torch.randperm(z.shape[0], device=z.device)]
        z_zero = torch.zeros_like(z)

        celltype_id = batch["celltype_id"]
        tech_id = _resolve_tech_id(batch)
        y_ord = batch["y_ord"]

        dec_real = self.system.decoder(z, celltype_id, tech_id, y_ord)
        dec_perm = self.system.decoder(z_perm, celltype_id, tech_id, y_ord)
        dec_zero = self.system.decoder(z_zero, celltype_id, tech_id, y_ord)

        gamma = self.system.decoder.coeff_from_z(z)

        rec_real = dec_real.nll_per_cell.mean()
        rec_perm = dec_perm.nll_per_cell.mean()
        rec_zero = dec_zero.nll_per_cell.mean()
        z_norm = z.norm(dim=-1)
        sigma_p = det_out.sigma_p.detach()
        mu_delta = det_out.encoder_out.mu_q.detach() - det_out.mu_p.detach()

        residual_clamp = float(getattr(self.system.prior, "residual_clamp", 0.0))
        if residual_clamp > 0.0:
            clamp_eps = max(1e-6, residual_clamp * 1e-6)
            clamp_hit_rate = (z.abs() >= (residual_clamp - clamp_eps)).float().mean()
        else:
            clamp_hit_rate = torch.zeros((), dtype=z.dtype, device=z.device)

        record = {
            "event": "z_debug_det",
            "epoch": getattr(self, "current_epoch_index", None),
            "step_in_epoch": int(step_idx),
            "batch_size": int(z.shape[0]),
            "rank": int(torch.distributed.get_rank())
            if torch.distributed.is_available() and torch.distributed.is_initialized()
            else 0,
            "rec_real": float(rec_real.detach().cpu()),
            "rec_perm": float(rec_perm.detach().cpu()),
            "rec_zero": float(rec_zero.detach().cpu()),
            "rec_perm_minus_real": float((rec_perm - rec_real).detach().cpu()),
            "rec_zero_minus_real": float((rec_zero - rec_real).detach().cpu()),
            "state_abs": float(det_out.decoder_out.state_score.abs().mean().detach().cpu()),
            "base_abs": float(det_out.decoder_out.base_score.abs().mean().detach().cpu()),
            "tech_abs": float(det_out.decoder_out.tech_score.abs().mean().detach().cpu()),
            "gamma_abs": float(gamma.abs().mean().detach().cpu()),
            "z_norm_mean": float(z_norm.mean().detach().cpu()),
            "z_norm_min": float(z_norm.min().detach().cpu()),
            "z_norm_max": float(z_norm.max().detach().cpu()),
            "sigma_p_mean": float(sigma_p.mean().detach().cpu()),
            "sigma_p_min": float(sigma_p.min().detach().cpu()),
            "sigma_p_max": float(sigma_p.max().detach().cpu()),
            "mu_q_minus_mu_p_norm_mean": float(mu_delta.norm(dim=-1).mean().detach().cpu()),
            "residual_clamp": residual_clamp,
            "residual_clamp_hit_rate": float(clamp_hit_rate.detach().cpu()),
        }

        if self._is_rank0():
            print(
                f"[z-sensitivity | step {step_idx:05d}]\n"
                f"  Reconstruction : real={record['rec_real']:.6f}  "
                f"permuted_z={record['rec_perm']:.6f}  "
                f"zero_z={record['rec_zero']:.6f}  "
                f"perm_minus_real={record['rec_perm_minus_real']:.6f}  "
                f"zero_minus_real={record['rec_zero_minus_real']:.6f}\n"
                f"  Decoder scores : state_abs={record['state_abs']:.6e}  "
                f"base_abs={record['base_abs']:.6e}  "
                f"tech_abs={record['tech_abs']:.6e}  "
                f"gamma_abs={record['gamma_abs']:.6e}\n"
                f"  Latent prior   : z_norm={record['z_norm_mean']:.6f}  "
                f"sigma_p={record['sigma_p_mean']:.6f}  "
                f"mu_q_minus_mu_p_norm={record['mu_q_minus_mu_p_norm_mean']:.6f}  "
                f"clamp_hit={record['residual_clamp_hit_rate']:.4f}"
            )

            if self.debug_z_sensitivity_path is not None:
                self.debug_z_sensitivity_path.parent.mkdir(parents=True, exist_ok=True)
                with open(self.debug_z_sensitivity_path, "a", encoding="utf-8") as f:
                    f.write(json.dumps(record, ensure_ascii=False, sort_keys=True))
                    f.write("\n")

        if was_training:
            self.system.train()

    def _format_step_log(self, *, prefix: str, step_idx: int, meters: "RunningAverages") -> str:
        m = meters.compute()
        nan = float("nan")

        def _finite_value(key: str, default: float = nan) -> float:
            try:
                v = float(m.get(key, default))
            except (TypeError, ValueError):
                return default
            return v if math.isfinite(v) else default

        def _has(key: str) -> bool:
            return math.isfinite(_finite_value(key))

        def _fmt(key: str, fmt: str = ".3f", default: str = "-") -> str:
            v = _finite_value(key)
            if not math.isfinite(v):
                return default
            return format(v, fmt)

        def _item(label: str, value: str) -> str:
            return f"   {label:<10} {value}"

        def _section(title: str) -> None:
            lines.append("")
            lines.append(title)

        def _verdict(text: str) -> None:
            lines.append(_item("verdict", text))

        loss_total = _finite_value("loss/total")
        rec = _finite_value("loss/rec")
        kl_state = _finite_value("loss/kl_state")
        bam_kl = _finite_value("loss/bam_kl")
        grad = _finite_value("optim/grad_norm")
        lr = _finite_value("optim/lr")
        stable = (
            "BAD"
            if not math.isfinite(loss_total)
            else "WATCH"
            if math.isfinite(grad) and grad > 10.0
            else "STABLE"
        )

        console_style = os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower()
        if console_style in {"prism_experiment", "prism_informative"}:
            # A middle ground between the old multi-page diagnostic and the
            # over-compressed four-line experiment log.  Only live mechanisms
            # are shown; every raw metric is still retained in history.json.
            previous = getattr(self, "_prev_diag_metrics", {})

            def _delta(key: str) -> float:
                current = _finite_value(key)
                try:
                    old = float(previous[key])
                except (KeyError, TypeError, ValueError):
                    return nan
                return current - old if math.isfinite(current) and math.isfinite(old) else nan

            def _signed(value: float, digits: int = 4) -> str:
                return f"{value:+.{digits}f}" if math.isfinite(value) else "-"

            lines = [
                "═" * 92,
                (
                    f"[{prefix} step {step_idx:05d}] {stable} | "
                    f"total {_fmt('loss/total', '.4f')} | "
                    f"PRISM {_fmt('loss/prism_branch_nll', '.4f')} | "
                    f"grad {_fmt('optim/grad_norm', '.2f')} | lr {_fmt('optim/lr', '.2e')}"
                ),
            ]

            reconstruction = (
                "  [복원] "
                f"nonzero 정확도 {_fmt('metric/ordinal_nonzero_acc', '.3f')} · "
                f"±1칸 {_fmt('metric/ordinal_nonzero_within1', '.3f')} · "
                f"balanced recall {_fmt('metric/ordinal_balanced_recall', '.3f')} · "
                f"raw/nonzero NLL {_fmt('loss/prism_branch_raw_nll', '.4f')}/"
                f"{_fmt('loss/prism_branch_nonzero_nll', '.4f')}"
            )
            if _has("metric/asym_nonzero_to_zero_leak"):
                reconstruction += (
                    " · nonzero→0 누락 "
                    f"{_fmt('metric/asym_nonzero_to_zero_leak', '.3f')}"
                )
            lines.append(reconstruction)

            if _has("loss/prism_branch_nll"):
                lines.append(
                    "  [PRISM] "
                    "정상지역/나이 "
                    f"{_fmt('metric/prism_normal_region_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_age_module_rms', '.4f')} · "
                    "공통/개인/반응 "
                    f"{_fmt('metric/prism_common_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_personal_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_response_module_rms', '.4f')} · "
                    f"support {_fmt('metric/prism_support_contexts', '.1f')} "
                    f"(신뢰 {_fmt('metric/prism_support_reliability', '.2f')})"
                )

            warnings_now = []
            support_query_enabled = False
            try:
                live_system = getattr(self.system, "module", self.system)
                live_head = getattr(live_system, "precision_head", None)
                live_cfg = getattr(live_head, "cfg", None)
                support_query_enabled = (
                    live_cfg is not None
                    and float(getattr(live_cfg, "lambda_support_query_infonce", 0.0)) > 0.0
                )
            except Exception:
                support_query_enabled = False

            if support_query_enabled and _has("loss/prism_support_query_infonce"):
                own = _finite_value("metric/prism_support_query_own_nll")
                negative = _finite_value("metric/prism_support_query_negative_nll")
                gap = negative - own
                own_delta = _delta("metric/prism_support_query_own_nll")
                negative_delta = _delta("metric/prism_support_query_negative_nll")
                if not (math.isfinite(own_delta) and math.isfinite(negative_delta)):
                    infonce_status = "기준점 수집중"
                elif gap <= 0.0:
                    infonce_status = "분리 안됨"
                    warnings_now.append("InfoNCE gap≤0")
                elif own_delta >= -1.0e-4 and negative_delta > 1.0e-4:
                    infonce_status = "주의: 남의 코드만 악화"
                    warnings_now.append("InfoNCE negative-only")
                elif own_delta < -1.0e-4:
                    infonce_status = "자기 donor 개선 포함"
                else:
                    infonce_status = "혼합 변화"
                lines.append(
                    "  [개인화·InfoNCE] "
                    f"loss {_fmt('loss/prism_support_query_infonce', '.4f')} · "
                    f"자기 {own:.4f} (Δ{_signed(own_delta)}) · "
                    f"매칭 타인 {negative:.4f} (Δ{_signed(negative_delta)}) · "
                    f"격차 {gap:+.4f} · valid {_fmt('metric/prism_support_query_valid_fraction', '.2f')} · "
                    f"anchor {_fmt('loss/prism_support_query_anchor', '.4f')} | {infonce_status}"
                )

            if _has("metric/agp_effective_tokens"):
                effective = _finite_value("metric/agp_effective_tokens")
                top1 = _finite_value("metric/agp_top1_weight")
                overlap = _finite_value("metric/agp_head_overlap")
                if top1 < 0.01 and overlap > 0.90:
                    agp_status = "주의: 거의 균등·head 중복"
                    warnings_now.append("AGP uniform pooling")
                elif overlap > 0.90:
                    agp_status = "주의: head 중복"
                    warnings_now.append("AGP duplicate heads")
                else:
                    agp_status = "선택적 pooling 작동"
                lines.append(
                    "  [AGP] "
                    f"실질 토큰 {effective:.1f} · top1 {top1:.3f} · "
                    f"top5 {_fmt('metric/agp_top5_mass', '.3f')} · "
                    f"head overlap {overlap:.3f} | {agp_status}"
                )

            if _has("loss/prism_pathology_leak"):
                nuisance = (
                    "  [누출·교란] "
                    "병리/나이/state "
                    f"{_fmt('loss/prism_pathology_leak', '.4f')}/"
                    f"{_fmt('loss/prism_age_leak', '.4f')}/"
                    f"{_fmt('loss/prism_state_pathology_leak', '.4f')}"
                )
                if _has("metric/tech_adv_bal_acc"):
                    nuisance += (
                        " · tech adv bal.acc "
                        f"{_fmt('metric/tech_adv_bal_acc', '.3f')} "
                        f"(λ {_fmt('weight/lambda_tech_adv_active', '.4f')})"
                    )
                if _has("metric/tech_adv_pred_major_frac"):
                    major = _finite_value("metric/tech_adv_pred_major_frac")
                    minority = _finite_value("metric/tech_adv_min_recall_present")
                    nuisance += f" · pred-major {major:.3f}/min-recall {minority:.3f}"
                    if major >= 0.95 and minority <= 0.10:
                        warnings_now.append("tech adversary class-collapse")
                if _has("metric/sex_adv_bal_acc"):
                    nuisance += (
                        " · sex adv bal.acc "
                        f"{_fmt('metric/sex_adv_bal_acc', '.3f')}"
                    )
                lines.append(nuisance)

            uncertainty_items = []
            for key, label, fmt in (
                ("unc/raw_mean", "entropy", ".3f"),
                ("unc/clip_mean", "u", ".3f"),
                ("weight/bam_w_min", "w_min", ".3f"),
                ("weight/bam_w_max", "w_max", ".3f"),
                ("loss/v6_unc_rec_alignment", "rec-align", ".3f"),
                ("loss/v6_unc_depth_corr", "depth-r²", ".3f"),
            ):
                if _has(key):
                    uncertainty_items.append(f"{label} {_fmt(key, fmt)}")
            if uncertainty_items:
                lines.append("  [불확실성] " + " · ".join(uncertainty_items))

            if _has("metric/pathology_rank_capacity"):
                rank_modes = (
                    ("fixed_rank_warmup", "rank8-warmup"),
                    ("all_open", "all-open"),
                    ("soft_rank_learning", "rank-soft"),
                    ("hard_straight_through", "rank-hard-ST"),
                    ("frozen_hard", "rank-frozen"),
                )
                rank_mode = next(
                    (
                        label
                        for key, label in rank_modes
                        if _finite_value(
                            f"metric/pathology_rank_mode_{key}", 0.0
                        )
                        >= 0.5
                    ),
                    "-",
                )
                generator_mode = (
                    "generator-warmup"
                    if _finite_value("metric/generator_warmup_all_on", 0.0)
                    >= 0.5
                    else "generator-shadow"
                    if _finite_value("metric/generator_mode_shadow", 0.0)
                    >= 0.5
                    else "generator-soft"
                    if _finite_value("metric/generator_mode_soft", 0.0) >= 0.5
                    else "generator-hard"
                )
                lines.append(
                    "  [구조용량] "
                    f"path-rank hard/E={_fmt('metric/pathology_rank_hard', '.1f')}/"
                    f"{_fmt('metric/pathology_rank_expected', '.1f')} of "
                    f"{_fmt('metric/pathology_rank_capacity', '.0f')} "
                    f"({rank_mode}) · search E="
                    f"{_fmt('metric/pathology_rank_search_expected', '.1f')} · "
                    f"freeze spread/near="
                    f"{_fmt('metric/pathology_rank_freeze_count_spread', '.1f')}/"
                    f"{_fmt('metric/pathology_rank_freeze_near_threshold_uncertain_fraction', '.2f')} · "
                    f"gen {_fmt('metric/generator_hard_active', '.1f')}/"
                    f"{_fmt('metric/generator_candidate_count', '.0f')} "
                    f"({generator_mode})"
                )

            if math.isfinite(grad) and grad > 10.0:
                warnings_now.append(f"large grad={grad:.2f}")
            lines.append(
                "  [판정] "
                + ("주의 — " + "; ".join(warnings_now) if warnings_now else "현재 집계창에서 뚜렷한 이상 없음")
            )
            lines.append("═" * 92)
            try:
                self._prev_diag_metrics = {
                    key: float(value)
                    for key, value in m.items()
                    if isinstance(value, (int, float)) and math.isfinite(float(value))
                }
            except Exception:
                self._prev_diag_metrics = {}
            return "\n".join(lines)

        # ==================================================================
        # v19 COMPACT diagnostic.  One tight block per category (header with a
        # status word + 1-2 data lines), a BLANK LINE between categories.  The
        # full 8-section detail further below is emitted ONLY when
        # ``self._verbose_diag`` is set (default off) — nothing is lost, but the
        # default log is readable.  All values go through the safe ``_fmt`` /
        # ``_finite_value`` helpers (missing key -> "-"), so this never raises.
        # ==================================================================
        prev = getattr(self, "_prev_diag_metrics", {})

        def _tr(key: str) -> str:
            """Trend arrow ↑/↓/→ vs the previous diagnostic (blank if unknown)."""
            cur = _finite_value(key)
            try:
                p = float(prev.get(key))
            except (TypeError, ValueError):
                return ""
            if not (math.isfinite(cur) and math.isfinite(p)):
                return ""
            d = cur - p
            return " →" if abs(d) < 1e-4 else (" ↑" if d > 0 else " ↓")

        def _cat(title: str, status: str) -> None:
            lines_c.append("")
            lines_c.append(f"  {title:<22} {status}")

        def _group(title: str) -> None:
            """Top-level section divider grouping related categories (visual hierarchy)."""
            lines_c.append("")
            lines_c.append(f"  ━━ {title} " + "━" * 40)

        mode = "fp32"
        if getattr(self, "amp_enabled", False):
            mode = "bf16" if getattr(self, "amp_dtype", None) is torch.bfloat16 else "fp16"
        head_r = f"{mode} · lr {lr:.1e}" if math.isfinite(lr) else mode
        rule = "═" * 80

        lam_align = _finite_value("weight/lambda_align", 0.0)
        lam_sex = _finite_value("weight/lambda_align_sex", 0.0)
        zero_gap = _finite_value("metric/ordinal_pred_zero_frac") - _finite_value("metric/ordinal_true_zero_frac")

        # status words (same thresholds as the verbose verdicts further below)
        if not math.isfinite(loss_total):
            st_loss = "✗ 안됨 (loss 비유한!)"
        elif math.isfinite(grad) and grad > 10.0:
            st_loss = "⚠ 주의 (grad 큼)"
        else:
            st_loss = "✓ 잘됨 (안정)"
        st_recon = "⚠ 주의 (zero gap 큼)" if (math.isfinite(zero_gap) and abs(zero_gap) > 0.10) else "⚠ 정상한계 (노이즈천장)"
        st_state = "✓ 잘됨 (collapse 없음)" if _has("metric/v6_state_fraction") else "ℹ 진단 없음"
        st_clean = "→ 진행중 (sex 제거중)" if lam_sex > 0.0 else ("→ 진행중 (celltype 정렬)" if lam_align > 0.0 else "ℹ 정렬 loss 없음")
        st_dis = "ℹ 현재 run은 epoch-end module rescue와 validation에서 측정"
        st_unc = "ℹ 진단용 (최종 신뢰도 아님)"

        schedule = getattr(self, "integrated_curriculum_schedule", {}) or {}
        epoch_now = int(getattr(self, "current_epoch_index", 0) or 0)
        total_epochs = int(schedule.get("total_epochs", 0) or 0)

        def _curriculum_state(name: str, ramp_key: str) -> str:
            start = int(schedule.get(name, 1) or 1)
            if epoch_now < start:
                return f"대기(E{start})"
            ramp_value = _finite_value(ramp_key, 1.0)
            if math.isfinite(ramp_value) and ramp_value < 0.995:
                return f"ramp {ramp_value:.2f}"
            return "ON"

        if console_style == "prism_integrated":
            epoch_text = (
                f"epoch {epoch_now:02d}/{total_epochs:02d}"
                if total_epochs > 0
                else f"epoch {epoch_now:02d}"
            )
            lines_c = [
                rule,
                f"  통합 PRISM · {epoch_text} · {prefix} step {step_idx:05d} · {head_r} · {stable}",
                rule,
            ]
            lines_c.append(
                "  ▸ curriculum  "
                "기본복원 ON · "
                f"PHU정렬 {_curriculum_state('phu_alignment_start', 'weight/phu_alignment_ramp')} · "
                f"PRISM분해 {_curriculum_state('explicit_start', 'weight/prism_output_ramp')} · "
                f"모듈복원 {'대기(E' + str(int(schedule.get('module_rescue_start', 1))) + ')' if epoch_now < int(schedule.get('module_rescue_start', 1)) else 'epoch-end ON'} · "
                f"PHU가중 {_curriculum_state('phu_relative_start', 'weight/phu_relative_ramp')} · "
                f"Generator {'대기(E' + str(int(schedule.get('generator_start', 1))) + ')' if epoch_now < int(schedule.get('generator_start', 1)) else 'ON'}"
            )
        else:
            lines_c = [rule, f"  {prefix} · step {step_idx:05d}    {head_r}    (✓ 좋음 / ⚠ 주의 / ✗ 문제 / → 진행중)", rule]
        lines_c.append(
            f"  ▸ 한눈에   {st_loss.split()[0]} 전반 "
            f"{'양호' if st_loss.startswith('✓') else ('주의' if st_loss.startswith('⚠') else '문제')}"
            f" · z 유효차원·BAM↔PHU·Generator 안전도를 아래에서 동시 추적 · test 봉인"
        )

        lines_c.append("")
        lines_c.append(f"  ■ 재구성  —  세포의 유전자 발현을 모델이 다시 만들어냄        {st_recon}")
        lines_c.append(
            f"     전체 손실 (작을수록 좋음) ......... {_fmt('loss/total')}{_tr('loss/total')}    "
            f"잠재정보량 KLz {_fmt('loss/kl_state', '.1f')} · BAM신뢰 {_fmt('loss/bam_kl', '.1f')} · 기울기 {_fmt('optim/grad_norm', '.2f')}(>10이면 불안정)"
        )
        lines_c.append(
            f"     세포 상태를 정확히 맞힘 .......... {_fmt('metric/ordinal_nonzero_acc')}    "
            f"한 칸 이내로 맞힘 {_fmt('metric/ordinal_nonzero_within1')}    "
            f"소수 중증단계 회수 {_fmt('metric/ordinal_balanced_recall')}"
        )
        recon2 = (
            f"     병변 없는 세포 맞힘 .............. {_fmt('metric/ordinal_true_zero_frac')}    "
            f"과하게 0으로 예측 {_fmt('metric/asym_nonzero_to_zero_leak')}"
        )
        if math.isfinite(zero_gap):
            recon2 += f"    0예측 과부족 {zero_gap:+.3f}"
        lines_c.append(recon2)
        lines_c.append("     ↳ 세포 하나하나는 노이즈가 커서 정확도가 낮은 게 정상 — 진짜 성적표는 아래 '질병 포착'")
        aux = []
        for key, label in (("loss/ordinal_nonzero", "순서"), ("loss/hier_weighted_total", "계층"), ("loss/hier_group", "그룹"), ("loss/hier_emd", "분포거리"), ("loss/hier_under_asym", "과소벌점")):
            if _has(key):
                aux.append(f"{label} {_fmt(key, '.2f')}")
        if aux:
            lines_c.append("     ┄ 보조손실(아는 사람용)  " + " · ".join(aux))

        lines_c.append("")
        lines_c.append(f"  ■ z 잠재변수 (사영 후 = 모델이 실제 쓰는 z)  —  잡음(성별·세포종류) 빼고 '질병 신호'만 남기기   {st_clean}")
        _z_collapse_text = (
            "아래 평균-z 유효차원에서 판정"
            if console_style == "prism_integrated"
            else ("정상 ✓" if _has("metric/v6_state_fraction") else "진단없음")
        )
        lines_c.append(
            f"     z 안 무너졌나 (붕괴 여부) ........ {_z_collapse_text}    "
            f"상태가중 {_fmt('metric/v6_state_fraction')}{_tr('metric/v6_state_fraction')} · "
            f"섞임게이트 베이스 {_fmt('metric/v6_mixer_gate_base_mean', '.2f')}/상태 {_fmt('metric/v6_mixer_gate_state_mean', '.2f')}/기술 {_fmt('metric/v6_mixer_gate_tech_mean', '.2f')}"
        )
        # sex/celltype HARD-removal status. removed = fraction of z energy projected out; sex 강도 = top
        # sex eigenvalue (REFILL proxy — if it doesn't shrink, z keeps re-encoding sex). post-projection.
        if _has("metric/proj_removed_energy"):
            lines_c.append(
                f"     성별(sex) 빼낸 양 ............... {_fmt('metric/proj_removed_energy', '.2f')} 제거    "
                f"sex 강도 {_fmt('metric/proj_sex_eig_top', '.3f')}{_tr('metric/proj_sex_eig_top')}(줄면 좋음=누출안참) · "
                f"실제 제거축 {_fmt('metric/proj_eff_rank', '.0f')}개(sex{_fmt('metric/proj_sex_rank', '.0f')}+ct{_fmt('metric/proj_ct_rank', '.0f')})"
            )
        if (
            _has("metric/sex_adv_bal_acc")
            and (
                _finite_value("weight/lambda_sex_adv_target", 0.0) > 0.0
                or _finite_value("weight/lambda_sex_adv_active", 0.0) > 0.0
            )
        ):                                                  # v31 soft sex-erasure adversary
            lines_c.append(
                f"     sex adversary (소프트 제거) ...... 맞히는정도 {_fmt('metric/sex_adv_bal_acc', '.2f')}"
                f"(0.5로 내려갈수록 encoder가 sex 숨김=좋음) · λ {_fmt('weight/lambda_sex_adv_active', '.3f')}"
            )
        _refp = (f" · 참조은행 {_fmt('metric/v6_ref_prior_active_celltypes', '.0f')}셀타입 학습중"
                 if _has("metric/v6_ref_prior_active_celltypes") else "")
        lines_c.append(
            f"     세포종류 정렬 (작을수록 좋음) ..... {_fmt('loss/align', '.3f')}{_tr('loss/align')}    "
            f"sex 정렬 {_fmt('loss/align_sex', '.3f')}(λ{lam_sex:.2g}) · z 펴기 {_fmt('loss/white', '.3f')} · 원점앵커 {_fmt('loss/ref_state', '.3f')}{_refp}"
        )

        # PRISM is the new explicit precision-medicine route.  Keep its core
        # health indicators in the default compact renderer: a full epoch can
        # take hours, so an operator attached to tmux must not have to wait for
        # the epoch-end summary to learn whether these branches are alive.
        if _has("loss/prism_branch_nll"):
            _prism_obj = _finite_value("loss/prism_branch_nll")
            _prism_ramp = _finite_value("weight/prism_ramp", 0.0)
            _prism_status = (
                "✗ 문제 (비유한)"
                if not math.isfinite(_prism_obj)
                else "→ 준비중 (ramp)"
                if _prism_ramp < 0.99
                else "✓ 작동중"
            )
            lines_c.append("")
            lines_c.append(
                "  ■ PRISM 정밀의학 분해  —  공통질병·개인기본·개인반응을 따로 학습        "
                + _prism_status
            )
            lines_c.append(
                f"     다른 세포로 표적세포 복원 ........ {_fmt('loss/prism_branch_nll', '.4f')}"
                f"{_tr('loss/prism_branch_nll')}    raw {_fmt('loss/prism_branch_raw_nll', '.4f')}"
                f" · 균형 {_fmt('loss/prism_branch_balanced_nll', '.4f')}"
                f" · 질병단계>0 {_fmt('loss/prism_branch_nonzero_nll', '.4f')}"
            )
            lines_c.append(
                f"     항 크기 정상지역/나이 ........... {_fmt('metric/prism_normal_region_module_rms', '.4f')}"
                f"/{_fmt('metric/prism_age_module_rms', '.4f')}    공통/개인/반응 "
                f"{_fmt('metric/prism_common_module_rms', '.4f')}/"
                f"{_fmt('metric/prism_personal_module_rms', '.4f')}/"
                f"{_fmt('metric/prism_response_module_rms', '.4f')}"
                f" · 누출 병리/나이/state {_fmt('loss/prism_pathology_leak', '.4f')}/"
                f"{_fmt('loss/prism_age_leak', '.4f')}/"
                f"{_fmt('loss/prism_state_pathology_leak', '.4f')}"
                f" · support {_fmt('metric/prism_support_contexts', '.1f')}"
                f"(신뢰 {_fmt('metric/prism_support_reliability', '.2f')}) · ramp {_fmt('weight/prism_ramp', '.2f')}"
            )

        if (
            console_style == "prism_integrated"
            and _finite_value("metric/prism_paper_branch_enabled", 0.0) >= 0.5
        ):
            _paper_threshold = (
                _finite_value(
                    "metric/prism_paper_branch_variant_compartmental_threshold",
                    0.0,
                )
                >= 0.5
            )
            _paper_variant = (
                "paper_compartmental_threshold"
                if _paper_threshold
                else "graph_linear_control"
            )
            _paper_ramp = _finite_value(
                "weight/prism_module_local_nonlinear_ramp", 0.0
            )
            _paper_status = (
                "→ 대기/램프중"
                if _paper_ramp < 0.995
                else "✓ ON (논문 가설 검증중)"
            )
            lines_c.append("")
            lines_c.append(
                "  ■ ref 논문 적용부  —  dendritic compartment·문턱·초선형성의 계산적 비유        "
                + _paper_status
            )
            lines_c.append(
                "     arm/활성도 ........................ "
                f"{_paper_variant} · ramp "
                f"{_fmt('weight/prism_module_local_nonlinear_ramp', '.2f')} "
                "(실제 수상돌기·NMDA 측정이라는 뜻은 아님)"
            )
            lines_c.append(
                "     비선형/기본 local RMS ............. "
                f"{_fmt('metric/prism_module_local_nonlinear_rms', '.5f')}/"
                f"{_fmt('metric/prism_module_local_module_rms', '.5f')} "
                f"· 비율 {_fmt('metric/prism_module_local_nonlinear_to_local_ratio', '.3f')} "
                f"· 문턱통과 {_fmt('metric/prism_module_local_threshold_crossing_fraction', '.2%')}"
            )
            lines_c.append(
                "     이 경로를 끄면 생기는 ΔNLL ....... branch/full "
                f"{_fmt('metric/prism_module_local_nonlinear_branch_nll_gain', '+.5f')}/"
                f"{_fmt('metric/prism_module_local_nonlinear_full_nll_gain', '+.5f')} "
                "(양수면 논문 적용 경로가 도움)"
            )
            if _paper_threshold:
                lines_c.append(
                    "     학습된 mix/threshold/slope/gain .. "
                    f"{_fmt('metric/prism_module_local_mix_mean', '.3f')}/"
                    f"{_fmt('metric/prism_module_local_threshold_mean', '.3f')}/"
                    f"{_fmt('metric/prism_module_local_slope_mean', '.3f')}/"
                    f"{_fmt('metric/prism_module_local_gain_mean', '.3f')}"
                )
            else:
                lines_c.append(
                    "     학습된 mix/linear gain ........... "
                    f"{_fmt('metric/prism_module_local_mix_mean', '.3f')}/"
                    f"{_fmt('metric/prism_module_local_linear_scale_mean', '.3f')}"
                )

        if console_style == "prism_integrated":
            module_start = int(schedule.get("module_rescue_start", 1) or 1)
            previous_rescue = getattr(
                self, "_integrated_last_module_rescue_metrics", None
            )
            lines_c.append("")
            lines_c.append(
                "  ■ 모듈 복원  —  donor·세포타입별 모듈 방향을 직접 맞춤        "
                + (
                    f"→ 대기(E{module_start}부터)"
                    if epoch_now < module_start
                    else "→ 현재 epoch 종료 후 donor-balanced 갱신"
                )
            )
            if previous_rescue:
                lines_c.append(
                    "     직전 epoch branch/full CCC ...... "
                    f"{float(previous_rescue.get('metric/prism_module_rescue_branch_flattened_ccc', float('nan'))):.3f}/"
                    f"{float(previous_rescue.get('metric/prism_module_rescue_full_flattened_ccc', float('nan'))):.3f}"
                    " · Pearson "
                    f"{float(previous_rescue.get('metric/prism_module_rescue_branch_pearson', float('nan'))):.3f}/"
                    f"{float(previous_rescue.get('metric/prism_module_rescue_full_pearson', float('nan'))):.3f}"
                    " · ramp "
                    f"{float(previous_rescue.get('weight/prism_module_rescue_ramp', float('nan'))):.2f}"
                    " · donor그룹 "
                    f"{float(previous_rescue.get('metric/prism_module_rescue_groups', float('nan'))):.1f}"
                )
            else:
                lines_c.append(
                    "     CCC/Pearson ..................... 아직 없음 · 활성화 epoch부터 매 epoch 종료 시 기록"
                )

        if _has("metric/pathology_rank_capacity"):
            _rank_modes = (
                ("fixed_rank_warmup", "fixed-rank8 warmup"),
                ("all_open", "all-open"),
                ("soft_rank_learning", "soft rank learning"),
                ("hard_straight_through", "hard straight-through"),
                ("frozen_hard", "frozen hard"),
            )
            _rank_mode = next(
                (
                    label
                    for key, label in _rank_modes
                    if _finite_value(
                        f"metric/pathology_rank_mode_{key}", 0.0
                    )
                    >= 0.5
                ),
                "unknown",
            )
            _freeze_ready = (
                _finite_value("metric/pathology_rank_freeze_ready", 0.0)
                >= 0.5
            )
            _rank_sparse = _finite_value(
                "weight/pathology_rank_sparsity_multiplier", 0.0
            )
            _generator_objective = (
                _has("metric/generator_warmup_all_on")
                and _finite_value("metric/generator_warmup_all_on", 1.0) < 0.5
            )
            _capacity_overlap = _rank_sparse > 0.0 and _generator_objective
            lines_c.append("")
            lines_c.append(
                "  ■ 구조 용량 자동학습  —  pathology rank를 확정한 뒤 generator를 줄임        "
                + ("✗ 개수학습 중첩" if _capacity_overlap else "✓ 역할분리")
            )
            lines_c.append(
                "     Pathology rank live hard/E ....... "
                f"{_fmt('metric/pathology_rank_hard', '.1f')}/"
                f"{_fmt('metric/pathology_rank_expected', '.1f')} / "
                f"{_fmt('metric/pathology_rank_capacity', '.0f')}    "
                f"search hard/E {_fmt('metric/pathology_rank_search_hard', '.1f')}/"
                f"{_fmt('metric/pathology_rank_search_expected', '.1f')} · {_rank_mode}"
            )
            lines_c.append(
                "     rank gate 온도/희소화 ............. "
                f"{_fmt('metric/pathology_rank_temperature', '.3f')}/"
                f"{_fmt('weight/pathology_rank_sparsity_multiplier', '.2f')}    "
                f"pathology scale {_fmt('weight/pathology_curriculum_scale', '.2f')} "
                f"· 실제 경로크기 {_fmt('metric/v26_pathdec_absmean', '.5f')}"
            )
            lines_c.append(
                "     rank 확정 low/mid/high ........... "
                f"{_fmt('metric/pathology_rank_freeze_count_low', '.1f')}/"
                f"{_fmt('metric/pathology_rank_freeze_count_mid', '.1f')}/"
                f"{_fmt('metric/pathology_rank_freeze_count_high', '.1f')} "
                f"· spread {_fmt('metric/pathology_rank_freeze_count_spread', '.1f')} "
                f"· 경계후보 {_fmt('metric/pathology_rank_freeze_near_threshold_uncertain_fraction', '.2%')} "
                f"· ready {'YES' if _freeze_ready else 'NO'} "
                f"· frozen {_fmt('metric/pathology_rank_finalized', '.0f')}"
            )
            lines_c.append(
                "     개수학습 분리 검사 ............... "
                f"rank penalty {'ON' if _rank_sparse > 0.0 else 'OFF'} · "
                f"generator objective {'ON' if _generator_objective else 'OFF'} · "
                f"overlap {'YES' if _capacity_overlap else 'NO'}"
            )

        # Online generator-count selection.  The live hard K is a running mean
        # of exact binary masks sampled by training batches; the deterministic
        # architecture is committed and reported separately at epoch end.
        if _has("metric/generator_hard_active"):
            _gen_warmup = (
                _finite_value("metric/generator_warmup_all_on", 0.0) >= 0.5
            )
            _gen_full_violation = _finite_value(
                "metric/generator_violation_full_nll"
            )
            _gen_isolated_violation = _finite_value(
                "metric/generator_violation_isolated_nll"
            )
            if _gen_warmup:
                _gen_status = "\u2192 \uc900\ube44\uc911 (warm-up\u00b7\uc804\uccb4 ON)"
            elif (
                math.isfinite(_gen_full_violation)
                and math.isfinite(_gen_isolated_violation)
                and _gen_full_violation <= 0.0
                and _gen_isolated_violation <= 0.0
            ):
                _gen_status = "\u2713 \uc131\ub2a5 \uc548\uc804\uc870\uac74 \ud1b5\uacfc"
            else:
                _gen_status = (
                    "\u26a0 \uc131\ub2a5 \uc548\uc804\uc870\uac74 "
                    "\uc704\ubc18(\ud559\uc2b5\uc911)"
                )

            _gen_total = _finite_value("metric/generator_candidate_count")
            _gen_total_text = (
                f"{_gen_total:.0f}" if math.isfinite(_gen_total) else "-"
            )
            _gen_uncertain = _finite_value(
                "metric/generator_uncertain_fraction"
            )
            _gen_uncertain_text = (
                f"{100.0 * _gen_uncertain:.1f}%"
                if math.isfinite(_gen_uncertain)
                else "-"
            )

            def _generator_delta(candidate_key: str, baseline_key: str) -> str:
                candidate = _finite_value(candidate_key)
                baseline = _finite_value(baseline_key)
                if not (math.isfinite(candidate) and math.isfinite(baseline)):
                    return "-"
                return f"{candidate - baseline:+.4f}"

            lines_c.append("")
            lines_c.append(
                "  \u25a0 Generator \uc790\ub3d9\uc120\ud0dd  \u2014  "
                "\uc131\ub2a5\uc744 \uc9c0\ud0a4\uba70 \ud544\uc694\ud55c "
                "\uc0dd\uc131\uc790 \uc218\ub97c \ud559\uc2b5        "
                + _gen_status
            )
            lines_c.append(
                "     \ubc30\uce58\ubcc4 \uc2e4\uc81c ON \ud3c9\uade0 ........ "
                f"{_fmt('metric/generator_hard_active', '.1f')} / "
                f"{_gen_total_text}    \ud655\ub960\uc0c1 \uc608\uc0c1 "
                f"{_fmt('metric/generator_expected_active', '.1f')}"
                f" \u00b7 \uc628\ub3c4 {_fmt('metric/generator_temperature', '.3f')}"
                f" \u00b7 \ubd88\ud655\uc2e4 gate {_gen_uncertain_text}"
            )
            lines_c.append(
                "     \uc804\uccb4 \ubcf5\uc6d0 \uc601\ud5a5(NLL) ........ "
                f"{_fmt('metric/generator_baseline_full_nll', '.4f')} "
                f"\u2192 {_fmt('metric/generator_candidate_full_nll', '.4f')}"
                " \u00b7 \ubcc0\ud654 "
                f"{_generator_delta('metric/generator_candidate_full_nll', 'metric/generator_baseline_full_nll')}"
                " \u00b7 \uc548\uc804\ub3c4 "
                f"{_fmt('metric/generator_violation_full_nll', '+.3f')}"
                "(\u22640\uc774\uba74 \ud1b5\uacfc)"
            )
            lines_c.append(
                "     Generator \ub2e8\ub3c5 \uc601\ud5a5(NLL)  "
                f"{_fmt('metric/generator_baseline_isolated_nll', '.4f')} "
                f"\u2192 {_fmt('metric/generator_candidate_isolated_nll', '.4f')}"
                " \u00b7 \ubcc0\ud654 "
                f"{_generator_delta('metric/generator_candidate_isolated_nll', 'metric/generator_baseline_isolated_nll')}"
                " \u00b7 \uc548\uc804\ub3c4 "
                f"{_fmt('metric/generator_violation_isolated_nll', '+.3f')}"
                "(\u22640\uc774\uba74 \ud1b5\uacfc)"
            )
            lines_c.append(
                "     \u21b3 \uc2e4\uc2dc\uac04\uc740 batch hard-mask "
                "\ud3c9\uade0\u00b7\ud655\uc815 K\ub294 epoch \uc885\ub8cc "
                "\ud6c4 \ubcc4\ub3c4 \uae30\ub85d"
            )

        # z effective rank (clean = POST-projection = the z the decoder USES) + per-pathology-axis spread.
        # Guarded so a log-format error never breaks the step.
        try:
            if _has("metric/v20_z_participation") or _has("loss/v20_path_aux"):
                if console_style == "prism_integrated":
                    lines_c.append("")
                    lines_c.append(
                        "  ■ z 유효차원  —  posterior 평균이 32개 축을 실제로 얼마나 쓰는가        ℹ sampled z와 구분"
                    )
                _pr = _finite_value("metric/v20_z_participation") if _has("metric/v20_z_participation") else 0.0
                if console_style == "prism_integrated" and epoch_now <= 2:
                    _st_spread = "→ 초기 형성중(E3부터 붕괴 판정)"
                else:
                    _st_spread = ("✓ 다축(≥3)" if _pr >= 3.0 else ("→ 펴지는 중" if _pr >= 1.5 else "✗ 1축뿐(붕괴주의)"))
                _pr_s = _finite_value("metric/v20_z_participation_sampled") if _has("metric/v20_z_participation_sampled") else None
                _samp_note = (f" · 노이즈포함 {_pr_s:.1f}{' ⚠간극큼' if (_pr_s - _pr) >= 1.0 else ''}" if _pr_s is not None else "")
                _pr_raw_v = _finite_value("metric/v20_z_participation_raw") if _has("metric/v20_z_participation_raw") else None
                _raw_note = f" · 사영전 {_pr_raw_v:.1f}" if _pr_raw_v is not None else ""
                lines_c.append(
                    f"     z가 실제 쓰는 차원 수 ........... {_pr:.1f} / 32{_tr('metric/v20_z_participation')}    "
                    f"{_st_spread}{_raw_note}{_samp_note}  (질병 정보가 퍼진 축 개수·사영후)"
                )
                _axb = []
                for _ax in ("axis_braak_high", "axis_thal_high", "axis_late_positive", "axis_lewy_positive"):
                    if _has(f"loss/v20_{_ax}"):
                        _sn = _ax.replace("axis_", "").replace("_positive", "").replace("_high", "")
                        _axb.append(f"{_sn} {_finite_value(f'loss/v20_{_ax}'):.2f}")
                _path = f"병리축 학습 {_fmt('loss/v20_path_aux', '.2f')}{_tr('loss/v20_path_aux')}"
                lines_c.append(
                    f"     여러 병리축을 따로 잡나 ......... {_path}"
                    + (f"    축별 " + " ".join(_axb) if _axb else "")
                )
                if _has("metric/v20_specificity") or _has("metric/v20_residual_frac"):
                    _spec = _finite_value("metric/v20_specificity") if _has("metric/v20_specificity") else 0.0
                    _resid = _finite_value("metric/v20_residual_frac") if _has("metric/v20_residual_frac") else 0.0
                    lines_c.append(
                        f"     ┄ 축 구분도 {_spec:.3f}(낮을수록 축끼리 안 겹침) · 4축 너머 추가구조 {_resid:.2f}(z>4축 전엔 ~0이 정상)"
                    )
        except Exception:
            pass

        lines_c.append("")
        lines_c.append("  ■ 질병 포착  —  z가 알츠하이머 병리를 설명하나        ℹ 학습 중엔 간접측정")
        # pathology-tied decoder magnitude (re-added in v30). zero-init ⇒ grows from 0 as it learns
        # celltype-specific pathology gene patterns = the channel that routes disease INTO z.
        if _has("metric/v26_pathdec_absmean"):
            _syn = (f" · 시너지 {_fmt('metric/v27_synergy_absmean', '.4f')}" if _has("metric/v27_synergy_absmean") else "")
            lines_c.append(
                f"     병리→발현 디코더 작동 .......... 크기 {_fmt('metric/v26_pathdec_absmean', '.4f')}{_tr('metric/v26_pathdec_absmean')}    "
                f"(0→↑ = 병리가 celltype별 유전자패턴 구동 = 질병을 z로 끌어들이는 통로){_syn}"
            )
        if console_style == "prism_integrated":
            lines_c.append(
                "     질병 회복도 ..................... train 중 PRISM branch NLL·모듈 CCC로 추적 → validation checkpoint에서 확정"
            )
        else:
            lines_c.append(
                "     질병 회복도 (module_cross) ...... 학습 중엔 못 봄 → 체크포인트(posthoc)서 측정 · 최근 .65 (목표 ≥.60) · region 누출 낮게"
            )
        if _has("loss/decoder_ref_state_l2") or _has("loss/decoder_ref_base_nll"):
            lines_c.append(
                f"     정상앵커 decoder ............... state {_fmt('loss/decoder_ref_state_l2', '.4f')}"
                f" · baseNLL {_fmt('loss/decoder_ref_base_nll', '.3f')}"
                f" · base-full {_fmt('metric/decoder_ref_base_minus_full_nll', '.3f')}"
                f" · n_ref {_fmt('metric/decoder_ref_anchor_n_ref', '.0f')}"
                f"    (reference cells: base+tech(+sex)로 복원, state는 작게)"
            )
        # v32 stochastic-consistency (R-Drop): shown only when KMLEE_PATH_CONSISTENCY>0.
        if _has("loss/path_consistency"):
            lines_c.append(
                f"     병리일관성(R-Drop) ............. {_fmt('loss/path_consistency', '.4f')}{_tr('loss/path_consistency')}"
                f"  ·  λ{_fmt('weight/lambda_path_consistency', '.2f')}"
                f"  ·  그룹 {_fmt('metric/path_cons_n_kept', '.0f')}/{_fmt('metric/path_cons_n_groups', '.0f')}"
                f"  ·  단일셀 {_fmt('metric/path_cons_singleton_frac', '.2f')}"
                f"    (2번 forward 병리투영 일치도; ↓수렴=좋음 / 0급락+병리손상=head붕괴주의)"
            )

        # learnable Lie rank (2026-07): shown only when learnable_lie_rank is on. Interpret at LATE epochs.
        if _has("metric/lie_rank_PR_mean"):
            lines_c.append("")
            lines_c.append("  ■ 학습 rank  —  Lie 생성자 rank를 group-lasso로 자기선택 (rank-6 과잉공급→가지치기)   ℹ 후반 epoch만 신뢰")
            lines_c.append(
                f"     유효 rank (PR / 강한k95) ....... PR {_fmt('metric/lie_rank_PR_mean', '.2f')}{_tr('metric/lie_rank_PR_mean')}"
                f"  ·  k95 {_fmt('metric/lie_rank_k95_mean', '.2f')} / 6    (낮을수록 적은 rank 사용; 후반 spectrum/gap로만 해석)"
            )
            lines_c.append(
                f"     penalty (손실대비 목표0.02) .... frac {_fmt('metric/lie_penalty_frac', '.3f')}"
                f"  ·  pen {_fmt('metric/lie_penalty', '.2e')}  ·  grad_scale {_fmt('metric/lie_grad_scale', '.1e')}    (grad_scale 급증시 blow-up 주의)"
            )
            lines_c.append(
                f"     누출 검사 (채널 비중) .......... Lie {_fmt('metric/lie_chanfrac_linear', '.2f')}"
                f" · direct {_fmt('metric/lie_chanfrac_direct', '.2f')} · path {_fmt('metric/lie_chanfrac_pathology', '.2f')}"
                f" · trans {_fmt('metric/lie_chanfrac_translation', '.2f')}    (Lie비중↓+다른채널↑ = 변형이 새는중=오염)"
            )

        # CP-BAM (consensus-private decomposition): shown only when cp_bam.enabled.
        if _has("loss/cpbam_consensus"):
            lines_c.append("")
            lines_c.append("  ■ CP-BAM  —  donor 효과를 공통질병지문 / 개인지문으로 분해        ℹ 학습 중 간접측정")
            lines_c.append(
                f"     공통지문 합의 ......... {_fmt('loss/cpbam_consensus', '.4f')}{_tr('loss/cpbam_consensus')} (↓수렴=비슷한병리 donor끼리 공통지문 공유)"
                f"  ·  donor제거 acc {_fmt('metric/cpbam_donor_adv_acc', '.3f')} (↓좋음; 단 celltype+병리도 봐서 순수 z_common 누출은 posthoc서)"
                f"  ·  병리예측 {_fmt('loss/cpbam_path', '.3f')}"
            )
            lines_c.append(
                f"     분해 상태 ............. z_common {_fmt('metric/cpbam_zcommon_norm', '.2f')} / z_private {_fmt('metric/cpbam_zprivate_norm', '.2f')}"
                f" (z_private→0 = 개인부분 붕괴, stage2 필요)  ·  무손실 {_fmt('loss/cpbam_latent_rec', '.3f')}{_tr('loss/cpbam_latent_rec')}"
                f"  ·  bank {_fmt('metric/cpbam_bank_seen_frac', '.2f')}·valid그룹 {_fmt('metric/cpbam_valid_consensus', '.0f')}"
            )
            if _has("loss/cpbam_stage2_decoder_total"):
                lines_c.append(
                    f"     stage2 decoder ........ commonNLL {_fmt('loss/cpbam_common_decode_nll', '.3f')}"
                    f" · splitNLL {_fmt('loss/cpbam_split_full_decode_nll', '.3f')}"
                    f" · private잔차 {_fmt('loss/cpbam_private_resid_score', '.4f')}"
                    f" · ref0 {_fmt('loss/cpbam_ref_zero', '.4f')}"
                    f" · 합 {_fmt('loss/cpbam_stage2_decoder_total', '.4f')}"
                )

        lines_c.append("")
        lines_c.append("  ■ 신뢰도  —  모델이 세포별로 자기 예측을 얼마나 믿나        ℹ 진단용")
        cb = []
        if _has("loss/uncertainty_alignment"):
            cb.append(
                f"PHU정렬 loss {_fmt('loss/uncertainty_alignment', '.3f')}"
                f"{_tr('loss/uncertainty_alignment')}"
                f"·ramp {_fmt('weight/phu_alignment_ramp', '.2f')}"
            )
        if _has("weight/phu_relative_ramp"):
            cb.append(
                f"PHU상대가중[{_fmt('weight/phu_rel_w_min', '.2f')}·"
                f"{_fmt('weight/phu_rel_w_mean', '.2f')}·"
                f"{_fmt('weight/phu_rel_w_max', '.2f')}] "
                f"ramp={_fmt('weight/phu_relative_ramp', '.2f')} "
                f"ESS={_fmt('metric/phu_rel_ess_fraction', '.2f')} "
                f"BAM↔실제noise={_fmt('metric/phu_bam_target_corr', '+.2f')}"
            )
        if _has("loss/v6_unc_rec_alignment"):
            _fw_off = _has("metric/v6_firewall_active") and _finite_value("metric/v6_firewall_active") < 0.5
            cb.append((f"firewall OFF·진단만 {_fmt('loss/v6_unc_rec_alignment', '.3f')}") if _fw_off
                      else f"firewall {_fmt('loss/v6_unc_rec_alignment', '.3f')}{_tr('loss/v6_unc_rec_alignment')}")
        if _has("loss/v6_unc_depth_corr"):
            cb.append(f"깊이추적 {_finite_value('loss/v6_unc_depth_corr'):+.2f}(↑좋음)")
        if _has("weight/v6_bam_rel_w_mean"):
            cb.append(f"세포별신뢰가중[{_fmt('weight/v6_bam_rel_w_min', '.2f')}·{_fmt('weight/v6_bam_rel_w_mean', '.2f')}·{_fmt('weight/v6_bam_rel_w_max', '.2f')}]")
        if cb:
            lines_c.append("     불확실성 보정 ................... " + "  ·  ".join(cb))
        if console_style == "prism_integrated" and _has("diag/u_total_mean"):
            lines_c.append(
                "     BAM/PHU 분포 .................... "
                f"BAM entropy {_fmt('unc/raw_mean', '.3f')} · "
                f"PHU {_fmt('diag/u_total_mean', '.3f')}±{_fmt('diag/u_total_std', '.3f')} · "
                f"세포 noise {_fmt('diag/noise_score_mean', '.3f')}±{_fmt('diag/noise_score_std', '.3f')}"
            )
        rb = []
        if _has("metric/zero_reliability_mean"):
            rb.append(f"제로신뢰 {_fmt('metric/zero_reliability_mean', '.2f')}")
        if _has("metric/pi_tech_at_proven_dropout"):
            rb.append(f"기술dropout의심 {_fmt('metric/pi_tech_at_proven_dropout', '.2f')}")
        if _has("metric/tech_adv_bal_acc"):
            rb.append(f"기술잡음 못맞힘 {_fmt('metric/tech_adv_bal_acc', '.2f')}(0.5=batch정보 안 샘=좋음)")
        if _has("metric/bam_attn_entropy"):
            rb.append(f"attn {_fmt('metric/bam_attn_entropy', '.2f')}")
        if rb:
            lines_c.append("     측정 믿을만한가 ................. " + "  ·  ".join(rb))

        lines_c.append("")
        lines_c.append("─" * 80)
        if console_style == "prism_integrated":
            lines_c.append(
                "  ▸ epoch 판정   validation 복원·모듈 CCC/Pearson·z 유효차원·BAM↔noise·Generator 안전제약을 함께 확인"
            )
        else:
            lines_c.append("  ▸ 다음 확인   posthoc로 z→sex(누출 빠졌나)·z→disease(질병 잡나)·module_cross(≥.60) · ref z→celltype↓")
        lines_c.append(rule)

        # remember this step's metrics so the next diagnostic can show ↑/↓ trends
        try:
            self._prev_diag_metrics = {k: float(v) for k, v in m.items() if isinstance(v, (int, float))}
        except Exception:
            self._prev_diag_metrics = {}

        # Default: compact only.  The full 8-section detail below runs only under
        # the verbose flag (kept for deep dives / continuity with §12.7 docs).
        if not getattr(self, "_verbose_diag", False):
            return "\n".join(lines_c)

        lines = list(lines_c)
        lines.append("")
        lines.append("── verbose detail ──")

        # ------------------------------------------------------------------
        # 1. Reconstruction / ordinal diagnostics.
        # ------------------------------------------------------------------
        _section("1) RECONSTRUCTION")
        if "metric/ordinal_nonzero_acc" in m:
            zero_true = _finite_value("metric/ordinal_true_zero_frac")
            zero_pred = _finite_value("metric/ordinal_pred_zero_frac")
            zero_gap = zero_pred - zero_true if math.isfinite(zero_pred) and math.isfinite(zero_true) else nan
            nz_leak = _finite_value("metric/asym_nonzero_to_zero_leak")
            recon_bits = [
                f"nz_exact {_fmt('metric/ordinal_nonzero_acc')}",
                f"within1 {_fmt('metric/ordinal_nonzero_within1')}",
                f"bal_rec {_fmt('metric/ordinal_balanced_recall')}",
            ]
            zero_bits = [
                f"zero {_fmt('metric/ordinal_true_zero_frac')}/{_fmt('metric/ordinal_pred_zero_frac')}",
                f"gap {zero_gap:+.3f}" if math.isfinite(zero_gap) else "gap -",
            ]
            if math.isfinite(nz_leak):
                zero_bits.append(f"nz->0 {nz_leak:.3f}")
            if _has("metric/pred_entropy_pnz"):
                zero_bits.append(f"conf {_fmt('metric/pred_entropy_pnz')}")
            lines.append(_item("metrics", " | ".join(recon_bits + zero_bits)))

            verdict = "OK: per-cell is noise-limited; module/posthoc is the real gate"
            if math.isfinite(zero_gap) and abs(zero_gap) > 0.10:
                verdict = "WATCH: zero gap is large"
        else:
            verdict = "PENDING: no ordinal metrics accumulated yet"

        aux = []
        for key, label, fmt in (
            ("loss/ordinal_balanced", "ord_bal", ".4f"),
            ("loss/ordinal_nonzero", "ord_nz", ".4f"),
            ("loss/hier_weighted_total", "hier", ".4f"),
            ("loss/hier_zero", "zero", ".4f"),
            ("loss/hier_group", "group", ".4f"),
            ("loss/hier_emd", "emd", ".4f"),
            ("metric/hier_ramp_scale", "ramp", ".2f"),
        ):
            if _has(key):
                aux.append(f"{label} {format(_finite_value(key), fmt)}")
        if aux:
            lines.append(_item("losses", " | ".join(aux)))

        _verdict(verdict)

        # ------------------------------------------------------------------
        # 2. z usage / mixer.
        # ------------------------------------------------------------------
        _section("2) STATE / MIXER")
        usage_bits = []
        for key, label in (
            ("metric/v6_state_fraction", "state_frac"),
            ("metric/v6_state_abs", "state_abs"),
            ("metric/v6_base_abs", "base_abs"),
            ("metric/v6_tech_abs", "tech_abs"),
            ("metric/v6_tech_to_state", "tech/state"),
        ):
            if _has(key):
                usage_bits.append(f"{label} {_fmt(key)}")
        if usage_bits:
            lines.append(_item("usage", " | ".join(usage_bits)))
        gate_bits = []
        for key, label in (
            ("metric/v6_mixer_gate_base_mean", "base"),
            ("metric/v6_mixer_gate_tech_mean", "tech"),
            ("metric/v6_mixer_gate_state_mean", "state"),
        ):
            if _has(key):
                gate_bits.append(f"{label} {_fmt(key)}")
        if gate_bits:
            lines.append(_item("gates", " | ".join(gate_bits)))
        if _has("metric/v6_state_fraction") or _has("metric/v6_state_abs"):
            _verdict("OK: z/state path is active; no posterior-collapse signal in train log")
        else:
            _verdict("PENDING: no state-usage diagnostics in this window")

        # ------------------------------------------------------------------
        # 3. Reference alignment / whitening.
        # ------------------------------------------------------------------
        _section("3) REFERENCE z ALIGNMENT")
        origin_bits = []
        for key, label, fmt in (
            ("loss/ref_state", "ref_state", ".4f"),
            ("loss/ref_center", "ref_center", ".4f"),
            ("loss/ref_center_donor_balanced", "ref_db", ".4f"),
            ("loss/v6_coeff_zero_mean", "coeff_zm", ".4f"),
        ):
            if _has(key):
                origin_bits.append(f"{label} {format(_finite_value(key), fmt)}")
        prior_bits = []
        for key, label, fmt in (
            ("metric/v6_ref_prior_active_celltypes", "prior_ct", ".0f"),
            ("metric/v6_ref_prior_update_norm", "update", ".3f"),
            ("metric/v6_ref_prior_norm_mean", "norm", ".3f"),
        ):
            if _has(key):
                prior_bits.append(f"{label} {format(_finite_value(key), fmt)}")
        ref_bits = origin_bits + prior_bits
        if ref_bits:
            lines.append(_item("origin/prior", " | ".join(ref_bits)))

        ct_adv_active = (
            _finite_value("weight/lambda_celltype_adv_active", 0.0) > 0.0
            or _finite_value("weight/lambda_celltype_adv_max", 0.0) > 0.0
            or _finite_value("loss/celltype_adv", 0.0) > 0.0
        )
        if ct_adv_active and (_has("metric/celltype_adv_balanced_acc") or _has("loss/celltype_adv")):
            lines.append(
                _item(
                    "main_ct_adv",
                    f"w {_fmt('weight/lambda_celltype_adv_active', '.4f')}/{_fmt('weight/lambda_celltype_adv_max', '.4f')} | "
                    f"bacc {_fmt('metric/celltype_adv_balanced_acc')} | "
                    f"n_ref {_fmt('metric/celltype_adv_n_reference', '.0f')} | "
                    f"loss {_fmt('loss/celltype_adv', '.4f')}",
                )
            )
        if _finite_value("weight/celltype_adv_aux_enabled", 0.0) > 0.5:
            lines.append(
                _item(
                    "aux_ct_adv",
                    f"bacc {_fmt('metric/celltype_adv_aux_balanced_acc')} | "
                    f"n_ref {_fmt('metric/celltype_adv_aux_n_reference', '.0f')} | "
                    f"loss {_fmt('loss/celltype_adv_aux', '.4f')} | "
                    f"pretrain {_fmt('metric/celltype_adv_aux_pretrain', '.0f')}",
                )
            )

        lambda_align = _finite_value("weight/lambda_align", 0.0)
        lambda_white = _finite_value("weight/lambda_white", 0.0)
        lambda_align_sex = _finite_value("weight/lambda_align_sex", 0.0)
        if lambda_align > 0.0 or lambda_white > 0.0:
            lines.append(
                _item(
                    "align/white",
                    f"align {_fmt('loss/align', '.4f')} "
                    f"(mean {_fmt('loss/align_l_mean', '.4f')}, cov {_fmt('loss/align_l_cov', '.4f')}, EMA {_fmt('weight/align_cov_use_ema', '.0f')}) | "
                    f"white {_fmt('loss/white', '.4f')} | n_ref/ct {_fmt('metric/align_n_reference', '.0f')}/{_fmt('metric/align_n_celltypes_in_ref', '.0f')}",
                )
            )
        if lambda_align_sex > 0.0 or _has("loss/align_sex"):
            lines.append(
                _item(
                    "sex",
                    f"{_fmt('loss/align_sex', '.4f')} "
                    f"(mean {_fmt('loss/align_l_mean_sex', '.4f')}, cov {_fmt('loss/align_l_cov_sex', '.4f')}) | "
                    f"n_sex {_fmt('metric/align_n_sex_in_ref', '.0f')} | lambda {lambda_align_sex:.3g}",
                )
            )

        align_cov = _finite_value("loss/align_l_cov")
        if lambda_align_sex > 0.0:
            _verdict("WATCH: cov and sex terms must be active; posthoc decides z->celltype/sex and disease preservation")
        elif lambda_align > 0.0 or lambda_white > 0.0:
            if math.isfinite(align_cov) and align_cov <= 1e-8:
                _verdict("WATCH: mean may improve, but cov is inactive/dead; posthoc decides if this is enough")
            else:
                _verdict("WATCH: posthoc probe decides if ref z->celltype drops from +0.743")
        elif _finite_value("weight/lambda_celltype_adv_active", 0.0) > 0.0:
            ref_bacc = _finite_value("metric/celltype_adv_balanced_acc")
            if math.isfinite(ref_bacc) and ref_bacc <= 0.35:
                _verdict("OK: reference-celltype adversary is suppressing leakage")
            else:
                _verdict("WATCH: reference z still exposes celltype; should fall after ramp")
        else:
            _verdict("INFO: no active z-celltype erasure loss in this run")

        # ------------------------------------------------------------------
        # 4. Disease signal gates.
        # ------------------------------------------------------------------
        _section("4) DISEASE SIGNAL")
        if _has("loss/prism_branch_nll"):
            lines.append(
                _item(
                    "PRISM E2E",
                    f"donor-balanced masked objective {_fmt('loss/prism_branch_nll', '.4f')} "
                    f"(raw {_fmt('loss/prism_branch_raw_nll', '.4f')}, "
                    f"nz {_fmt('loss/prism_branch_nonzero_nll', '.4f')}) | "
                    f"normal-region/age "
                    f"{_fmt('metric/prism_normal_region_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_age_module_rms', '.4f')} | "
                    f"module RMS common/personal/response "
                    f"{_fmt('metric/prism_common_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_personal_module_rms', '.4f')}/"
                    f"{_fmt('metric/prism_response_module_rms', '.4f')} | "
                    f"gene-score {_fmt('metric/prism_gene_score_rms', '.4f')} | "
                    f"u {_fmt('metric/prism_code_rms', '.3f')} | "
                    f"support {_fmt('metric/prism_support_contexts', '.1f')} | "
                    f"ramp {_fmt('weight/prism_ramp', '.2f')}",
                )
            )
            axis_bits = []
            for _name in ("thal", "braak", "cerad", "late", "lewy"):
                axis_bits.append(
                    f"{_name} C{_fmt(f'metric/prism_common_{_name}_rms', '.3f')}"
                    f"/R{_fmt(f'metric/prism_response_{_name}_rms', '.3f')}"
                )
            lines.append(_item("5 axes", " | ".join(axis_bits) + " | ADNC input=NO"))
            _verdict("ACTIVE: explicit common/personal/response branches are in ordinal reconstruction")
        else:
            lines.append(_item("gate", "posthoc only: withinCT z->ADNC stable | module_cross >=0.60 (last 0.65) | z->region low"))
            _verdict("UNKNOWN in-train: require posthoc after checkpoint")

        # ------------------------------------------------------------------
        # 5. Tech invariance and adversaries.
        # ------------------------------------------------------------------
        _section("5) TECH INVARIANCE")
        tech_bits = []
        for key, label, fmt in (
            ("loss/tech_z_cons", "z_cons", ".4f"),
            ("loss/tech_rec_cons", "rec_cons", ".4f"),
        ):
            if _has(key):
                tech_bits.append(f"{label} {format(_finite_value(key), fmt)}")
        tech_summary = list(tech_bits)
        if _has("metric/risk_clean") or _has("metric/risk_tech"):
            tech_summary.append(f"risk {_fmt('metric/risk_clean')}/{_fmt('metric/risk_tech')}")
        if _has("metric/risk_disease"):
            tech_summary.append(f"risk_disease {_fmt('metric/risk_disease')}")
        if _has("metric/tech_adv_bal_acc") or _has("loss/tech_adv"):
            pred_major = _finite_value("metric/tech_adv_pred_major_frac")
            min_recall = _finite_value("metric/tech_adv_min_recall_present")
            tech_summary.append(
                f"adv_bacc {_fmt('metric/tech_adv_bal_acc')} | collapse {_fmt('metric/tech_adv_pred_major_frac')}/{_fmt('metric/tech_adv_min_recall_present')}"
            )
            lines.append(_item("tech", " | ".join(tech_summary) if tech_summary else "-"))
            if math.isfinite(pred_major) and (pred_major >= 0.85 or (math.isfinite(min_recall) and min_recall <= 0.20)):
                _verdict("WATCH: tech classifier skew/collapse; interpret assay adversary cautiously")
            else:
                _verdict("OK-ish: tech-invariance active; tech adversary not collapsed in this window")
        else:
            if tech_summary:
                lines.append(_item("tech", " | ".join(tech_summary)))
            _verdict("PENDING: no tech adversary metrics; tech-invariance losses may still be active")

        # ------------------------------------------------------------------
        # 6. Zero/dropout helper diagnostics.
        # ------------------------------------------------------------------
        _section("6) ZERO / DROPOUT")
        pi_bits = []
        thin_bits = []
        for key, label, fmt in (
            ("metric/pi_tech_at_proven_dropout", "pi@proven", ".3f"),
            ("metric/pi_tech_at_stable_zero", "pi@stable", ".3f"),
            ("metric/zero_reliability_mean", "zero_rel", ".3f"),
        ):
            if _has(key):
                pi_bits.append(f"{label} {format(_finite_value(key), fmt)}")
        if _has("metric/pi_tech_at_proven_dropout") and _has("metric/pi_tech_at_stable_zero"):
            pi_gap = _finite_value("metric/pi_tech_at_proven_dropout") - _finite_value("metric/pi_tech_at_stable_zero")
            pi_bits.append(f"gap {pi_gap:+.3f}")
        for key, label, fmt in (
            ("loss/thin_consistency", "thin_cons", ".4f"),
            ("loss/thin_tech_zero_sup", "tech_zero_sup", ".4f"),
            ("metric/thin_tech_zero_frac", "thin_zero_frac", ".3f"),
        ):
            if _has(key):
                thin_bits.append(f"{label} {format(_finite_value(key), fmt)}")
        zero_summary = []
        if pi_bits:
            zero_summary.append("pi " + " | ".join(pi_bits))
        if thin_bits:
            zero_summary.append("thin " + " | ".join(thin_bits))
        if zero_summary:
            lines.append(_item("zero", " || ".join(zero_summary)))
        lever_bits = []
        for key, label, fmt in (
            ("loss/hier_under_asym", "under_asym", ".4f"),
            ("metric/asym_nonzero_exp_gap", "exp_gap", "+.3f"),
            ("metric/asym_mean_exp_true3", "exp3", ".2f"),
            ("loss/var_floor", "var_floor", ".4f"),
            ("metric/var_gap", "var_gap", ".3f"),
            ("metric/var_state_depth_corr", "var_depth", ".3f"),
            ("loss/high_margin", "hi_margin", ".4f"),
            ("metric/tier4_recall", "hi_recall", ".3f"),
            ("metric/tier4_fp_on_zero", "hi_fp0", ".3f"),
        ):
            if _has(key):
                lever_bits.append(f"{label}={format(_finite_value(key), fmt)}")
        if lever_bits:
            lines.append(_item("levers", " | ".join(lever_bits)))
        _verdict("INFO: per-cell zero leak is not the main optimization target")

        # ------------------------------------------------------------------
        # 7. BAM / PHU uncertainty.
        # ------------------------------------------------------------------
        _section("7) UNCERTAINTY / BAM / PHU")
        bam_bits = []
        if "unc/raw_mean" in m:
            raw_min = _finite_value("unc/raw_min")
            raw_max = _finite_value("unc/raw_max")
            raw_gap = raw_max - raw_min if math.isfinite(raw_min) and math.isfinite(raw_max) else nan
            bam_bits.append(f"raw {_fmt('unc/raw_mean')} [{raw_min:.3f},{raw_max:.3f}] gap {raw_gap:.3f}")
        if "unc/clip_mean" in m:
            bam_bits.append(f"clip {_fmt('unc/clip_mean')} [{_fmt('unc/clip_min')},{_fmt('unc/clip_max')}]")
        if _has("loss/v6_unc_depth_corr"):
            bam_bits.append(f"depth_corr {_finite_value('loss/v6_unc_depth_corr'):+.3f}")
        if _has("metric/bam_attn_entropy"):
            bam_bits.append(f"attn_ent {_fmt('metric/bam_attn_entropy')}")
        if bam_bits:
            lines.append(_item("BAM", " | ".join(bam_bits)))
        if "loss/uncertainty_alignment" in m:
            lines.append(
                _item(
                    "PHU",
                    f"align {_fmt('loss/uncertainty_alignment', '.4f')} | "
                    f"u {_fmt('diag/u_total_mean')}±{_fmt('diag/u_total_std')} | "
                    f"noise {_fmt('diag/noise_score_mean')}±{_fmt('diag/noise_score_std')} | "
                    f"conf {_fmt('metric/pred_entropy_pnz')}",
                )
            )
        if "weight/bam_w_min" in m:
            lines.append(
                _item(
                    "BAM weights",
                    f"min {_fmt('weight/bam_w_min', '.4f')} | mean {_fmt('weight/bam_w_mean', '.4f')} | max {_fmt('weight/bam_w_max', '.4f')}",
                )
            )
        depth_corr = _finite_value("loss/v6_unc_depth_corr")
        if math.isfinite(depth_corr) and abs(depth_corr) > 0.50:
            _verdict("WATCH: BAM-attn remains depth-coupled; use it as diagnostic, not confidence")
        else:
            _verdict("DIAGNOSTIC ONLY: raw BAM/PHU are not final confidence; use decoder/MC posthoc gates")

        # ------------------------------------------------------------------
        # 8. Donor / region SNR and final gates.
        # ------------------------------------------------------------------
        _section("8) RUN GATES")
        snr_keys = [
            ("metric/donor_snr_state", "donor_SNR_state"),
            ("metric/donor_snr_decoder", "donor_SNR_decoder"),
            ("metric/region_cos_state", "region_cos_state"),
            ("metric/region_leak_probe", "region_leak_probe"),
        ]
        any_snr = False
        for key, label in snr_keys:
            if _has(key):
                lines.append(_item(label, _fmt(key)))
                any_snr = True
        if not any_snr:
            lines.append(_item("audit", "donor/region SNR not scheduled in this step"))
        lines.append(_item("continue", "ref z->celltype down + z->ADNC/module_cross stable"))
        lines.append(_item("posthoc", "ref z->celltype/sex, withinCT z->ADNC, module_cross; stop if module_cross<0.60 or z rank collapse"))
        if _has("metric/tech_adv_pred_major_frac"):
            lines.append(_item("[INFO ]", "tech adv collapse WARN is diagnostic, not a stop by itself"))

        if "diag/skip_ratio" in m and _finite_value("diag/skip_ratio", 0.0) > 0.0:
            lines.append(
                _item(
                    "DDP safety",
                    f"skip_ratio={m['diag/skip_ratio']:.4%}  "
                    f"skipped_batches={int(m.get('diag/skip_batches', 0))}  "
                    f"skipped_examples={int(m.get('diag/skip_examples', 0))}",
                )
            )
        return "\n".join(lines)

    @property
    def current_lr(self) -> float:
        return float(self.optimizer.param_groups[0]["lr"])


# ======================================================================
# Running averages
# ======================================================================
class RunningAverages:
    def __init__(self) -> None:
        self.weighted_sums: Dict[str, float] = {}
        # Some diagnostics are intentionally emitted only on thinning or
        # scheduled probe steps.  Track a denominator per key so a sparse
        # metric is averaged over the examples on which it was measured,
        # rather than being diluted by unrelated batches where it was absent.
        self.weight_sums: Dict[str, int] = {}
        self.n_examples: int = 0
        # Non-finite batches are excluded from the average so that a single
        # NaN/Inf batch does not poison `loss/total` (and therefore the
        # composite early-stopping metric / best checkpoint selection).
        self.n_finite_batches: int = 0
        self.n_skipped_batches: int = 0
        self.n_skipped_examples: int = 0

    def update_from_step(self, step_out: StepOutput) -> None:
        bsz = int(step_out.batch_size)

        total_raw = step_out.loss.details.get("loss/total")
        try:
            total_val = float(total_raw) if total_raw is not None else float("nan")
        except (TypeError, ValueError):
            total_val = float("nan")

        if not math.isfinite(total_val):
            self.n_skipped_batches += 1
            self.n_skipped_examples += bsz
            return

        self.n_finite_batches += 1
        self.n_examples += bsz

        for key, value in step_out.loss.details.items():
            try:
                v = float(value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(v):
                # Per-key non-finite values are also skipped so that this
                # batch's contribution to that specific key is dropped while
                # other finite keys are still accumulated.
                continue
            self.weighted_sums[key] = self.weighted_sums.get(key, 0.0) + v * bsz
            self.weight_sums[key] = self.weight_sums.get(key, 0) + bsz

        self.weighted_sums["optim/lr"] = self.weighted_sums.get("optim/lr", 0.0) + step_out.lr * bsz
        self.weight_sums["optim/lr"] = self.weight_sums.get("optim/lr", 0) + bsz
        if step_out.grad_norm is not None:
            try:
                gn = float(step_out.grad_norm)
            except (TypeError, ValueError):
                gn = float("nan")
            if math.isfinite(gn):
                self.weighted_sums["optim/grad_norm"] = (
                    self.weighted_sums.get("optim/grad_norm", 0.0) + gn * bsz
                )
                self.weight_sums["optim/grad_norm"] = (
                    self.weight_sums.get("optim/grad_norm", 0) + bsz
                )

    def compute(self) -> Dict[str, float]:
        if self.n_examples == 0:
            return {}
        out = {
            key: value / float(self.weight_sums[key])
            for key, value in self.weighted_sums.items()
            if self.weight_sums.get(key, 0) > 0
        }
        finalize_prism_epoch_metrics(out)
        if self.n_skipped_batches > 0:
            total_batches = self.n_skipped_batches + self.n_finite_batches
            out["diag/skip_batches"] = float(self.n_skipped_batches)
            out["diag/skip_examples"] = float(self.n_skipped_examples)
            out["diag/skip_ratio"] = (
                float(self.n_skipped_batches) / float(max(total_batches, 1))
            )
        return out


# ======================================================================
# Helpers
# ======================================================================
def _resolve_tech_id(batch: Dict[str, torch.Tensor]) -> torch.Tensor:
    if "tech_id" in batch:
        return batch["tech_id"]
    if "batch_id" in batch:
        return batch["batch_id"]
    raise KeyError("Batch must contain either 'tech_id' or 'batch_id'.")


def _resolve_module_token_padding_mask(
    batch: Dict[str, torch.Tensor],
    *,
    B: int,
    T: int,
    device: torch.device,
) -> Optional[torch.Tensor]:
    """
    Resolve a padding mask for the module-token sequence.

    Gene-level padding masks cannot be lifted mechanically after
    GeneModuleTokenizer, because the encoder sequence length is now M(+CLS),
    not G(+CLS). If masking is ever needed, the batch should supply either
    `module_token_padding_mask` or `token_padding_mask` with shape [B, T].
    """
    mask = batch.get("module_token_padding_mask")
    if mask is None:
        mask = batch.get("token_padding_mask")
    if mask is None:
        return None

    if mask.ndim != 2 or mask.shape != (B, T):
        raise ValueError(
            f"module token padding mask must have shape {(B, T)}, got {tuple(mask.shape)}."
        )
    if mask.dtype != torch.bool:
        raise TypeError(f"module token padding mask must be torch.bool, got {mask.dtype}.")

    return mask.to(device=device, non_blocking=True)


# ======================================================================
# Smoke test
# ======================================================================
if __name__ == "__main__":
    torch.manual_seed(42)

    B = 2
    G = 8
    M = 4
    K = 5
    d_model = 16
    d_z = 4
    n_celltypes = 3
    n_tech = 2

    gene_embedding = GeneExpressionEmbedding(
        n_genes=G,
        d_model=d_model,
        n_bins=K,
        use_continuous=True,
        use_cls=True,
        dropout=0.1,
    )

    membership = torch.zeros(M, G)
    membership[0, [0, 1, 2]] = 1.0
    membership[1, [2, 3, 4]] = 1.0
    membership[2, [4, 5, 6]] = 1.0
    membership[3, [1, 6, 7]] = 1.0

    module_tokenizer = GeneModuleTokenizer(
        membership_weight=membership,
        d_model=d_model,
        activity_weight=None,
        activity_default="l2_membership",
        pooling="mean",
        preserve_cls=True,
        use_module_id_embedding=True,
        use_activity_projection=True,
        dropout=0.1,
    )

    state_encoder = StateEncoder(
        d_model=d_model,
        n_heads=4,
        n_layers=1,
        d_z=d_z,
        n_celltypes=n_celltypes,
        condition_on_celltype=True,
        pooling="cls",
        use_cls_token=True,
        compute_cell_uncertainty=False,
        stochastic_attention=False,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
    )

    prior = CellTypePrior(n_celltypes=n_celltypes, d_z=d_z)

    decoder = LieActionOrdinalDecoder(
        n_genes=G,
        n_celltypes=n_celltypes,
        n_tech=n_tech,
        d_z=d_z,
        n_bins=K,
        n_generators=4,
        coeff_hidden_dim=8,
    )

    criterion = TotalLoss(
        beta_state=1.0,
        lambda_tech=1e-2,
        lambda_gauge=1e-2,
        lambda_bam=1.0,
        tau=1.0,
        r_min=0.0,
        r_max=5.0,
    )

    system = OrdinalBAMSystem(
        gene_embedding=gene_embedding,
        module_tokenizer=module_tokenizer,
        state_encoder=state_encoder,
        prior=prior,
        decoder=decoder,
    )

    optimizer = torch.optim.Adam(system.parameters(), lr=1e-3)
    trainer = Trainer(
        system=system,
        criterion=criterion,
        optimizer=optimizer,
        grad_clip_norm=1.0,
        amp=False,
    )

    batch = {
        "y_ord": torch.randint(0, K, (B, G), dtype=torch.long),
        "x_log1p": torch.rand(B, G),
        "x_gene_scalar": torch.randn(B, G),
        "celltype_id": torch.randint(0, n_celltypes, (B,), dtype=torch.long),
        "batch_id": torch.randint(0, n_tech, (B,), dtype=torch.long),
        "row_index": torch.arange(B, dtype=torch.long),
    }

    train_out = trainer.train_step(batch)
    eval_out = trainer.eval_step(batch)

    print("=" * 60)
    print("trainer.py — module-token + Lie decoder smoke test")
    print("=" * 60)
    print(f"train total = {train_out.loss.total.item():.6f}")
    print(f"eval total  = {eval_out.loss.total.item():.6f}")
    print(f"grad norm   = {train_out.grad_norm}")
    print(f"lr          = {train_out.lr:.6f}")

    assert torch.isfinite(train_out.loss.total)
    assert torch.isfinite(eval_out.loss.total)
    assert train_out.loss.weights_per_cell.shape == (B,)
    assert train_out.loss.rec_per_cell.shape == (B,)
    assert train_out.loss.kl_state_per_cell.shape == (B,)

    batch_dev = trainer._move_batch_to_device(batch)

    with torch.no_grad():
        model_out = system(batch_dev, sample_latent=False)
    assert model_out.gene_tokens.shape == (B, G + 1, d_model)
    assert model_out.module_tokens.shape == (B, M, d_model)
    assert model_out.tokens.shape == (B, M + 1, d_model)
    assert model_out.module_activity is not None
    assert model_out.module_activity.shape == (B, M)

    latents = trainer.collect_latents([batch])
    assert latents["mu_q"].shape == (B, d_z)
    assert latents["z_perp"].shape == (B, d_z)
    print("[OK] trainer.py module-token + Lie decoder smoke test passed.")
