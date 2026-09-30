"""End-to-end precision-medicine score decomposition for PRISM-BAM.

The module is deliberately *module first*.  It predicts additive amplitudes in
the same 414 biological-module coordinates used by the main tokenizer, then the
system lifts their sum into gene-score space with the tokeniser's fixed,
L2-normalised activity dictionary.  The resulting score participates directly
in the ordinal reconstruction likelihood, so this is part of the trained model
rather than a post-hoc adapter.

The personal code is amortised from pathology/nuisance-residualised donor
contexts.  Donor identifiers are used only to retrieve the observed support
table; there is no learned donor embedding.  All contexts belonging to the
target cell type are masked before the code is inferred, which prevents the
target expression from being copied into its own prediction and permits honest
use on donor-disjoint validation/test donors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Optional, Sequence

import torch
from torch import nn
import torch.nn.functional as F

from kmlee_bam.model.pathology_interactions import (
    PairwisePathologyInteractionBasis,
)


PATHOLOGY_NAMES: tuple[str, ...] = ("thal", "braak", "cerad", "late", "lewy")


@dataclass
class PrecisionMedicineConfig:
    """Configuration for the explicit normal/disease/personal decomposition."""

    enabled: bool = False
    context_npz_path: Optional[str] = None
    pathology_names: Sequence[str] = field(default_factory=lambda: PATHOLOGY_NAMES)
    personal_rank: int = 2
    context_hidden_dim: int = 64
    module_token_dim: int = 32
    latent_token_dim: int = 16
    use_source_latent: bool = False
    context_embedding_dim: int = 16
    context_n_heads: int = 4
    context_n_layers: int = 2
    context_dropout: float = 0.10
    context_input_dropout: float = 0.10
    context_reliability_cap: int = 100
    n_pathology_features: int = 3
    region_gate_scale: float = 0.20
    support_mask_mode: str = "celltype"
    # Loss weights.  The branch-only NLL is ramped by the trainer.
    lambda_branch_nll: float = 0.12
    branch_balanced_weight: float = 0.10
    branch_nonzero_weight: float = 0.14
    lambda_normal_context_l2: float = 2.0e-5
    lambda_common_l2: float = 1.0e-5
    lambda_personal_l2: float = 2.0e-5
    lambda_response_l2: float = 8.0e-5
    lambda_response_zero_mean: float = 2.0e-3
    lambda_code_l2: float = 2.0e-4
    lambda_code_center: float = 2.0e-3
    lambda_pathology_leak: float = 1.0e-2
    lambda_age_leak: float = 5.0e-3
    lambda_state_pathology_leak: float = 2.0e-2
    # Integrated-run curriculum. Historical configs retain their exact
    # epoch-1 behaviour through these defaults.
    loss_start_epoch: int = 1
    loss_ramp_epochs: int = 3
    output_start_epoch: int = 1
    output_ramp_epochs: int = 0
    scale_regularizers_with_loss_ramp: bool = False
    # Optimiser multiplier relative to optim.lr for precision_head.*.
    lr_multiplier: float = 25.0
    context_ridge: float = 10.0
    min_context_cells: int = 10
    # Optional common pathology-axis interactions.  The interaction feature
    # artifact is fitted on training donors only and is exactly zero at normal.
    interaction_enabled: bool = False
    interaction_rank: int = 2
    interaction_ridge: float = 1.0e-3
    interaction_min_complete_donors: int = 24
    interaction_scale_floor: float = 1.0e-4
    interaction_warmup_epochs: int = 8
    interaction_ramp_epochs: int = 4
    lambda_interaction_l2: float = 5.0e-5
    # Optional support -> query personal-code discrimination.  Other-donor
    # negatives are selected only from the training-donor pool, require the
    # exact same observed support-context set after target-celltype exclusion,
    # and are nearest-matched on disease, age, sex, and capped support depth.
    lambda_support_query_infonce: float = 0.0
    support_query_n_negatives: int = 4
    support_query_temperature: float = 0.01
    support_query_anchor_weight: float = 1.0
    # Target-celltype-excluded module-local personal residual.  This complete
    # feature family is opt-in; disabled historical runs instantiate no new
    # parameters or persistent buffers.
    module_local_enabled: bool = False
    module_local_reliability_npz_path: Optional[str] = None
    module_local_rank: int = 8
    module_local_start_epoch: int = 17
    module_local_ramp_epochs: int = 8
    module_local_rescue_start_epoch: int = 25
    module_local_output_cap_npz_path: Optional[str] = None
    module_local_output_cap_source_checkpoint_sha256: Optional[str] = None
    module_local_output_cap_source_config_sha256: Optional[str] = None
    module_local_output_cap_quantile: float = 0.995
    module_local_input_clip_quantile: float = 0.995
    module_local_projection_ridge: float = 1.0e-4
    module_local_size_ratio_cap: float = 1.0
    module_local_lr_multiplier: float = 10.0
    lambda_module_local_center: float = 2.0e-3
    lambda_module_local_hierarchy: float = 1.0e-4
    lambda_module_local_size: float = 2.0e-3
    lambda_module_local_pathology_leak: float = 1.0e-2
    lambda_module_local_age_leak: float = 5.0e-3
    # Paper-inspired compartmental threshold nonlinearity.  This remains a
    # separate default-off feature so the audited v1 path and old checkpoints
    # instantiate no extra parameters or buffers.
    module_local_nonlinear_enabled: bool = False
    # ``compartmental_threshold`` is the paper-inspired arm.
    # ``graph_linear_control`` is a function-class control with the same fixed
    # graph, curriculum, and shared mix/gain initialisation.  Threshold/slope
    # tensors remain in the checkpoint layout but are frozen and unused; this
    # arm is deliberately not described as capacity- or parameter-matched.
    module_local_nonlinear_variant: str = "compartmental_threshold"
    module_local_compartment_graph_npz_path: Optional[str] = None
    module_local_nonlinear_start_epoch: int = 21
    module_local_nonlinear_ramp_epochs: int = 4
    module_local_compartment_mix_max: float = 0.50
    module_local_threshold_min: float = 0.50
    module_local_threshold_max: float = 2.00
    module_local_slope_min: float = 1.00
    module_local_slope_max: float = 8.00
    module_local_gain_max: float = 1.00


@dataclass
class PrecisionMedicineOutput:
    """Tensors required by the likelihood, diagnostics, and later fingerprints."""

    normal_region_coeff: torch.Tensor      # [B, M]
    age_coeff: torch.Tensor                # [B, M]
    common_axis_coeff: torch.Tensor       # [B, K, M]
    personal_coeff: torch.Tensor          # [B, M], rank-2 + module-local
    personal_rank2_coeff: torch.Tensor    # [B, M]
    module_local_coeff: torch.Tensor      # [B, M], curriculum-scaled
    module_local_raw_coeff: torch.Tensor  # [B, M], before curriculum scale
    module_local_code: torch.Tensor       # [B, R], or [B,0] when disabled
    module_local_personal_span_overlap: torch.Tensor  # [B], max row cosine
    module_local_support_coverage: torch.Tensor  # [B]
    module_local_zero_coverage_fraction: torch.Tensor  # [B]
    module_local_nonlinear_coeff: torch.Tensor  # [B, M], curriculum-scaled
    module_local_nonlinear_raw_coeff: torch.Tensor  # [B, M], pre-curriculum
    module_local_threshold_crossing_fraction: torch.Tensor  # [B]
    response_axis_coeff: torch.Tensor     # [B, K, M]
    interaction_pair_coeff: torch.Tensor  # [B, Q, M], Q=0 when disabled
    interaction_feature: torch.Tensor     # [B, Q]
    interaction_valid: torch.Tensor       # [B, Q]
    total_module_coeff: torch.Tensor      # [B, M]
    personal_code: torch.Tensor           # [B, P]
    support_count: torch.Tensor           # [B]
    support_reliability: torch.Tensor     # [B]
    pathology: torch.Tensor               # [B, K], normalised to [0, 1]
    pathology_valid: torch.Tensor         # [B, K]
    pathology_leak_loss: torch.Tensor
    age_leak_loss: torch.Tensor
    state_pathology_leak_loss: torch.Tensor
    module_local_center_loss: torch.Tensor
    module_local_size_loss: torch.Tensor
    module_local_pathology_leak_loss: torch.Tensor
    module_local_age_leak_loss: torch.Tensor
    # Filled by OrdinalBAMSystem after the fixed module -> gene lift.
    gene_score: Optional[torch.Tensor] = None
    branch_nll_per_cell: Optional[torch.Tensor] = None
    branch_nll_per_gene: Optional[torch.Tensor] = None
    module_local_gene_score: Optional[torch.Tensor] = None
    module_local_nonlinear_gene_score: Optional[torch.Tensor] = None
    module_local_nonlinear_off_branch_nll_per_cell: Optional[torch.Tensor] = None
    module_local_nonlinear_off_full_nll_per_cell: Optional[torch.Tensor] = None
    module_local_off_branch_nll_per_cell: Optional[torch.Tensor] = None
    module_local_off_full_nll_per_cell: Optional[torch.Tensor] = None
    support_query_candidate_nll: Optional[torch.Tensor] = None
    support_query_negative_valid: Optional[torch.Tensor] = None
    support_query_common_nll: Optional[torch.Tensor] = None
    support_query_negative_donor: Optional[torch.Tensor] = None
    support_query_match_distance: Optional[torch.Tensor] = None


def delayed_epoch_ramp(epoch: int, *, start_epoch: int, ramp_epochs: int) -> float:
    """Return a deterministic 1-based delayed linear curriculum scale."""

    if isinstance(epoch, bool) or int(epoch) <= 0:
        raise ValueError("epoch must be a positive 1-based integer")
    start = max(1, int(start_epoch))
    if int(epoch) < start:
        return 0.0
    if int(ramp_epochs) <= 0:
        return 1.0
    return float(min(1.0, (int(epoch) - start + 1) / float(ramp_epochs)))
class _GradReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, value: torch.Tensor, strength: float) -> torch.Tensor:
        ctx.strength = float(strength)
        return value.view_as(value)

    @staticmethod
    def backward(ctx, grad: torch.Tensor):
        return -ctx.strength * grad, None


def _grad_reverse(value: torch.Tensor, strength: float = 1.0) -> torch.Tensor:
    return _GradReverse.apply(value, float(strength))


class _ScaleGradientByCell(torch.autograd.Function):
    """Keep a score unchanged while donor-balancing only its backward path."""

    @staticmethod
    def forward(ctx, value: torch.Tensor, cell_weight: torch.Tensor) -> torch.Tensor:
        if value.ndim < 2 or cell_weight.ndim != 1 or value.shape[0] != cell_weight.shape[0]:
            raise ValueError("cell_weight must be [B] for a score whose first dimension is B")
        ctx.save_for_backward(cell_weight)
        return value.view_as(value)

    @staticmethod
    def backward(ctx, grad: torch.Tensor):
        (cell_weight,) = ctx.saved_tensors
        shape = (cell_weight.shape[0],) + (1,) * (grad.ndim - 1)
        return grad * cell_weight.reshape(shape).to(grad.dtype), None


def donor_balanced_gradient(
    value: torch.Tensor,
    cell_weight: torch.Tensor,
) -> torch.Tensor:
    """Identity in the forward pass; inverse-donor weighted in backward.

    The ordinary full-cell reconstruction remains numerically unchanged, but
    its gradient into the explicit PRISM score no longer lets donors with many
    captured cells dominate the common/personal/response parameters.
    """

    return _ScaleGradientByCell.apply(value, cell_weight.float().reshape(-1))


def pathology_hinge_features(pathology: torch.Tensor) -> torch.Tensor:
    """Three monotone, zero-at-normal features per named pathology axis."""

    if pathology.ndim != 2 or pathology.shape[1] != len(PATHOLOGY_NAMES):
        raise ValueError(
            f"pathology must be [B,{len(PATHOLOGY_NAMES)}], got {tuple(pathology.shape)}"
        )
    return torch.stack(
        (
            pathology,
            F.relu(pathology - (1.0 / 3.0)),
            F.relu(pathology - (2.0 / 3.0)),
        ),
        dim=-1,
    )


class PrecisionMedicineHead(nn.Module):
    """Infer personal context and explicit module-level pathology components."""

    @property
    def module_local_nonlinear_trainable_families(self) -> tuple[str, ...]:
        """Raw parameter families that participate in this function class."""

        if self.module_local_nonlinear_variant == "graph_linear_control":
            return ("mix", "gain")
        return ("mix", "threshold", "slope", "gain")

    def __init__(
        self,
        *,
        config: PrecisionMedicineConfig,
        source_module: torch.Tensor,
        source_latent: torch.Tensor,
        source_observed: torch.Tensor,
        source_reliability: torch.Tensor,
        context_celltype: torch.Tensor,
        context_region: torch.Tensor,
        n_celltypes: int,
        n_regions: int,
        interaction_basis: Optional[PairwisePathologyInteractionBasis] = None,
        donor_pathology: Optional[torch.Tensor] = None,
        donor_pathology_valid: Optional[torch.Tensor] = None,
        donor_age_z: Optional[torch.Tensor] = None,
        donor_sex: Optional[torch.Tensor] = None,
        donor_adnc: Optional[torch.Tensor] = None,
        fit_donor_mask: Optional[torch.Tensor] = None,
        module_local_reliability: Optional[torch.Tensor] = None,
        module_local_input_clip: Optional[float] = None,
        module_local_output_cap: Optional[float] = None,
        module_local_compartment_adjacency: Optional[torch.Tensor] = None,
    ) -> None:
        super().__init__()
        self.cfg = config
        self.n_pathology = len(PATHOLOGY_NAMES)
        if tuple(str(x).lower() for x in config.pathology_names) != PATHOLOGY_NAMES:
            raise ValueError(
                "precision pathology_names must be exactly "
                f"{PATHOLOGY_NAMES}; ADNC is intentionally not an input"
            )
        if int(config.n_pathology_features) != 3:
            raise ValueError("the audited PRISM basis uses exactly three hinge features")
        if config.support_mask_mode != "celltype":
            raise ValueError("support_mask_mode must be 'celltype' for target-celltype exclusion")
        for name in ("loss_start_epoch", "output_start_epoch"):
            value = getattr(config, name)
            if isinstance(value, bool) or int(value) <= 0:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("loss_ramp_epochs", "output_ramp_epochs"):
            value = getattr(config, name)
            if isinstance(value, bool) or int(value) < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if source_module.ndim != 3 or source_latent.ndim != 3:
            raise ValueError("source context tensors must have shape [donor,context,feature]")
        if source_module.shape[:2] != source_latent.shape[:2]:
            raise ValueError("module/latent context tables disagree")
        if source_observed.shape != source_module.shape[:2]:
            raise ValueError("source_observed has the wrong shape")
        if source_reliability.shape != source_module.shape[:2]:
            raise ValueError("source_reliability has the wrong shape")
        if context_celltype.shape != source_module.shape[1:2]:
            raise ValueError("context_celltype has the wrong shape")
        if context_region.shape != source_module.shape[1:2]:
            raise ValueError("context_region has the wrong shape")

        self.n_donors = int(source_module.shape[0])
        self.n_contexts = int(source_module.shape[1])
        self.n_modules = int(source_module.shape[2])
        self.latent_dim = int(source_latent.shape[2])
        self.n_celltypes = int(n_celltypes)
        self.n_regions = int(n_regions)
        self.personal_rank = int(config.personal_rank)
        self.module_local_enabled = bool(config.module_local_enabled)
        self.module_local_nonlinear_enabled = bool(
            config.module_local_nonlinear_enabled
        )
        self.module_local_nonlinear_variant = str(
            config.module_local_nonlinear_variant
        )
        if self.module_local_nonlinear_variant not in {
            "compartmental_threshold",
            "graph_linear_control",
        }:
            raise ValueError(
                "module_local_nonlinear_variant must be "
                "'compartmental_threshold' or 'graph_linear_control'"
            )
        if self.module_local_nonlinear_enabled and not self.module_local_enabled:
            raise ValueError(
                "module_local_nonlinear_enabled requires module_local_enabled"
            )
        if isinstance(config.module_local_rank, bool) or int(
            config.module_local_rank
        ) <= 0:
            raise ValueError("module_local_rank must be a positive integer")
        if int(config.module_local_rank) > self.n_modules:
            raise ValueError("module_local_rank cannot exceed the module count")
        for name in ("module_local_start_epoch", "module_local_rescue_start_epoch"):
            value = getattr(config, name)
            if isinstance(value, bool) or int(value) <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if isinstance(config.module_local_ramp_epochs, bool) or int(
            config.module_local_ramp_epochs
        ) < 0:
            raise ValueError("module_local_ramp_epochs must be non-negative")
        if isinstance(config.module_local_nonlinear_start_epoch, bool) or int(
            config.module_local_nonlinear_start_epoch
        ) <= 0:
            raise ValueError(
                "module_local_nonlinear_start_epoch must be a positive integer"
            )
        if isinstance(config.module_local_nonlinear_ramp_epochs, bool) or int(
            config.module_local_nonlinear_ramp_epochs
        ) < 0:
            raise ValueError(
                "module_local_nonlinear_ramp_epochs must be non-negative"
            )
        for name in (
            "module_local_projection_ridge",
            "module_local_size_ratio_cap",
            "module_local_lr_multiplier",
        ):
            value = float(getattr(config, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for lower_name, upper_name in (
            ("module_local_threshold_min", "module_local_threshold_max"),
            ("module_local_slope_min", "module_local_slope_max"),
        ):
            lower = float(getattr(config, lower_name))
            upper = float(getattr(config, upper_name))
            if (
                not math.isfinite(lower)
                or not math.isfinite(upper)
                or lower <= 0.0
                or upper <= lower
            ):
                raise ValueError(
                    f"{lower_name}/{upper_name} must be finite, positive, and ordered"
                )
        for name in (
            "module_local_compartment_mix_max",
            "module_local_gain_max",
        ):
            value = float(getattr(config, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        clip_quantile = float(config.module_local_input_clip_quantile)
        if not math.isfinite(clip_quantile) or not 0.5 < clip_quantile < 1.0:
            raise ValueError("module_local_input_clip_quantile must lie in (0.5,1)")
        cap_quantile = float(config.module_local_output_cap_quantile)
        if not math.isfinite(cap_quantile) or not 0.5 < cap_quantile < 1.0:
            raise ValueError("module_local_output_cap_quantile must lie in (0.5,1)")
        for name in (
            "lambda_module_local_center",
            "lambda_module_local_hierarchy",
            "lambda_module_local_size",
            "lambda_module_local_pathology_leak",
            "lambda_module_local_age_leak",
        ):
            value = float(getattr(config, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        if self.module_local_enabled:
            if module_local_reliability is None:
                raise ValueError(
                    "module_local_enabled requires a sealed reliability artifact"
                )
            if module_local_reliability.shape != (
                self.n_contexts,
                self.n_modules,
            ):
                raise ValueError(
                    "module_local_reliability must have shape "
                    f"[{self.n_contexts},{self.n_modules}]"
                )
            reliability_float = module_local_reliability.float()
            if not bool(torch.isfinite(reliability_float).all()):
                raise ValueError("module_local_reliability contains non-finite values")
            if bool((reliability_float < 0.0).any()) or bool(
                (reliability_float > 1.0).any()
            ):
                raise ValueError("module_local_reliability must lie in [0,1]")
            if not bool((reliability_float > 0.0).any()):
                raise ValueError("module_local_reliability has no positive support")
            if module_local_input_clip is None:
                raise ValueError(
                    "module_local_enabled requires a train-only input clip"
                )
            resolved_clip = float(module_local_input_clip)
            if not math.isfinite(resolved_clip) or resolved_clip <= 0.0:
                raise ValueError("module_local_input_clip must be finite and positive")
            if module_local_output_cap is None:
                raise ValueError(
                    "module_local_enabled requires a frozen-comparator output cap"
                )
            resolved_output_cap = float(module_local_output_cap)
            if not math.isfinite(resolved_output_cap) or resolved_output_cap <= 0.0:
                raise ValueError(
                    "module_local_output_cap must be finite and positive"
                )
            if self.module_local_nonlinear_enabled:
                if module_local_compartment_adjacency is None:
                    raise ValueError(
                        "nonlinear module-local requires a sealed compartment graph"
                    )
                adjacency = module_local_compartment_adjacency.float()
                if adjacency.shape != (self.n_modules, self.n_modules):
                    raise ValueError(
                        "module-local compartment adjacency must have shape "
                        f"[{self.n_modules},{self.n_modules}]"
                    )
                if not bool(torch.isfinite(adjacency).all()):
                    raise ValueError(
                        "module-local compartment adjacency contains non-finite values"
                    )
                if bool((adjacency < 0.0).any()):
                    raise ValueError(
                        "module-local compartment adjacency must be non-negative"
                    )
                if not bool(
                    torch.equal(
                        torch.diagonal(adjacency),
                        torch.zeros_like(torch.diagonal(adjacency)),
                    )
                ):
                    raise ValueError(
                        "module-local compartment adjacency diagonal must be zero"
                    )
                row_sum = adjacency.sum(dim=1)
                nonempty = row_sum > 0.0
                if bool(
                    (torch.abs(row_sum[nonempty] - 1.0) > 1.0e-5).any()
                ):
                    raise ValueError(
                        "non-empty compartment adjacency rows must sum to one"
                    )
            elif module_local_compartment_adjacency is not None:
                raise ValueError(
                    "compartment adjacency was supplied while nonlinear local is disabled"
                )
        elif (
            module_local_reliability is not None
            or module_local_input_clip is not None
            or module_local_output_cap is not None
            or module_local_compartment_adjacency is not None
        ):
            raise ValueError(
                "module-local artifact inputs were supplied while the feature is disabled"
            )
        self.interaction_enabled = bool(config.interaction_enabled)
        if isinstance(config.interaction_rank, bool) or int(config.interaction_rank) <= 0:
            raise ValueError("interaction_rank must be a positive integer")
        if int(config.interaction_warmup_epochs) < 0:
            raise ValueError("interaction_warmup_epochs must be non-negative")
        if int(config.interaction_ramp_epochs) < 0:
            raise ValueError("interaction_ramp_epochs must be non-negative")
        if float(config.lambda_interaction_l2) < 0.0:
            raise ValueError("lambda_interaction_l2 must be non-negative")
        if float(config.lambda_support_query_infonce) < 0.0:
            raise ValueError("lambda_support_query_infonce must be non-negative")
        if isinstance(config.support_query_n_negatives, bool) or int(
            config.support_query_n_negatives
        ) <= 0:
            raise ValueError("support_query_n_negatives must be positive")
        if float(config.support_query_temperature) <= 0.0:
            raise ValueError("support_query_temperature must be positive")
        if float(config.support_query_anchor_weight) < 0.0:
            raise ValueError("support_query_anchor_weight must be non-negative")
        if self.interaction_enabled and interaction_basis is None:
            raise ValueError(
                "interaction_enabled requires a train-donor-fitted interaction_basis"
            )
        if not self.interaction_enabled and interaction_basis is not None:
            raise ValueError(
                "interaction_basis was supplied while interaction_enabled is false"
            )
        self.interaction_basis = interaction_basis

        # The context table is an input artifact, not a trainable donor lookup.
        self.register_buffer("source_module", source_module.float().clone(), persistent=True)
        self.register_buffer("source_latent", source_latent.float().clone(), persistent=True)
        self.register_buffer("source_observed", source_observed.bool().clone(), persistent=True)
        self.register_buffer(
            "source_reliability",
            source_reliability.float().clamp(0.0, 1.0).clone(),
            persistent=True,
        )
        self.register_buffer("context_celltype", context_celltype.long().clone(), persistent=True)
        self.register_buffer("context_region", context_region.long().clone(), persistent=True)

        if self.module_local_enabled:
            assert module_local_reliability is not None
            assert module_local_input_clip is not None
            self.register_buffer(
                "module_local_reliability",
                module_local_reliability.float().clamp(0.0, 1.0).clone(),
                persistent=True,
            )
            self.register_buffer(
                "module_local_input_clip",
                torch.tensor(float(module_local_input_clip), dtype=torch.float32),
                persistent=True,
            )
            self.register_buffer(
                "module_local_output_cap",
                torch.tensor(float(module_local_output_cap), dtype=torch.float32),
                persistent=True,
            )
            self.register_buffer(
                "module_local_scale",
                torch.tensor(0.0, dtype=torch.float32),
                persistent=True,
            )
            if self.module_local_nonlinear_enabled:
                assert module_local_compartment_adjacency is not None
                self.register_buffer(
                    "module_local_compartment_adjacency",
                    module_local_compartment_adjacency.float().clone(),
                    persistent=True,
                )
                self.register_buffer(
                    "module_local_nonlinear_scale",
                    torch.tensor(0.0, dtype=torch.float32),
                    persistent=True,
                )
                # Python mirror avoids a GPU -> CPU synchronisation in every
                # module-local forward.  It is refreshed by the epoch setter
                # and after checkpoint loading.
                self._module_local_nonlinear_scale_value = 0.0
            else:
                self.module_local_compartment_adjacency = None
                self.module_local_nonlinear_scale = None
                self._module_local_nonlinear_scale_value = 0.0
        else:
            self.module_local_reliability = None
            self.module_local_input_clip = None
            self.module_local_output_cap = None
            self.module_local_scale = None
            self.module_local_compartment_adjacency = None
            self.module_local_nonlinear_scale = None
            self._module_local_nonlinear_scale_value = 0.0

        self.support_query_enabled = float(config.lambda_support_query_infonce) > 0.0
        if self.support_query_enabled:
            required_matching = {
                "donor_pathology": donor_pathology,
                "donor_pathology_valid": donor_pathology_valid,
                "donor_age_z": donor_age_z,
                "donor_sex": donor_sex,
                "donor_adnc": donor_adnc,
                "fit_donor_mask": fit_donor_mask,
            }
            missing_matching = [
                name for name, value in required_matching.items() if value is None
            ]
            if missing_matching:
                raise ValueError(
                    "support-query InfoNCE requires donor matching metadata: "
                    f"{missing_matching}"
                )
            assert donor_pathology is not None
            assert donor_pathology_valid is not None
            assert donor_age_z is not None
            assert donor_sex is not None
            assert donor_adnc is not None
            assert fit_donor_mask is not None
            if donor_pathology.shape != (self.n_donors, self.n_pathology):
                raise ValueError("donor_pathology must be [n_donors,5]")
            if donor_pathology_valid.shape != donor_pathology.shape:
                raise ValueError("donor_pathology_valid shape mismatch")
            for name, value in (
                ("donor_age_z", donor_age_z),
                ("donor_sex", donor_sex),
                ("donor_adnc", donor_adnc),
                ("fit_donor_mask", fit_donor_mask),
            ):
                if value.shape != (self.n_donors,):
                    raise ValueError(f"{name} must have shape [n_donors]")
            # These are reconstructed deterministically from the sealed context
            # artifact at every run and therefore need not bloat checkpoints.
            self.register_buffer(
                "match_pathology", donor_pathology.float().clone(), persistent=False
            )
            self.register_buffer(
                "match_pathology_valid",
                donor_pathology_valid.bool().clone(),
                persistent=False,
            )
            self.register_buffer(
                "match_age_z", donor_age_z.float().clone(), persistent=False
            )
            self.register_buffer(
                "match_sex", donor_sex.float().clone(), persistent=False
            )
            self.register_buffer(
                "match_adnc", donor_adnc.float().clone(), persistent=False
            )
            self.register_buffer(
                "match_fit_donor", fit_donor_mask.bool().clone(), persistent=False
            )

        self.module_projection = nn.Sequential(
            nn.Linear(self.n_modules, int(config.module_token_dim)),
            nn.LayerNorm(int(config.module_token_dim)),
            nn.GELU(),
        )
        self.use_source_latent = bool(config.use_source_latent)
        if self.use_source_latent:
            self.latent_projection = nn.Sequential(
                nn.Linear(self.latent_dim, int(config.latent_token_dim)),
                nn.LayerNorm(int(config.latent_token_dim)),
                nn.GELU(),
            )
        else:
            self.latent_projection = None
        emb = int(config.context_embedding_dim)
        self.context_celltype_embedding = nn.Embedding(self.n_celltypes, emb)
        self.context_region_embedding = nn.Embedding(self.n_regions, emb)
        token_in = int(config.module_token_dim) + 2 * emb
        if self.use_source_latent:
            token_in += int(config.latent_token_dim)
        hidden = int(config.context_hidden_dim)
        self.context_token_projection = nn.Linear(token_in, hidden)
        layer = nn.TransformerEncoderLayer(
            d_model=hidden,
            nhead=int(config.context_n_heads),
            dim_feedforward=4 * hidden,
            dropout=float(config.context_dropout),
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.context_encoder = nn.TransformerEncoder(layer, num_layers=int(config.context_n_layers))
        self.personal_posterior = nn.Sequential(
            nn.LayerNorm(hidden),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Dropout(float(config.context_dropout)),
            nn.Linear(hidden, self.personal_rank),
            nn.Tanh(),
        )

        k, f, m = self.n_pathology, 3, self.n_modules
        # Explicit non-disease context terms.  Without these, ordinary ageing
        # and the DLPFC/MTG normal offset can be spuriously absorbed by the
        # named pathology or personal branches.
        self.normal_region_delta = nn.Parameter(
            torch.zeros(self.n_celltypes, self.n_regions, m)
        )
        self.age_global = nn.Parameter(torch.zeros(m))
        self.age_celltype_delta = nn.Parameter(torch.zeros(self.n_celltypes, m))
        # Global disease programme plus zero-centred cell-type deviation.
        self.common_global = nn.Parameter(torch.zeros(k, f, m))
        self.common_celltype_delta = nn.Parameter(torch.zeros(self.n_celltypes, k, f, m))
        self.common_region_gate = nn.Parameter(torch.zeros(self.n_regions, k, f))
        # Baseline-like personal state and low-capacity personal response.
        self.personal_basis = nn.Parameter(
            torch.zeros(self.n_celltypes, self.personal_rank, m)
        )
        self.response_basis = nn.Parameter(
            torch.zeros(self.n_celltypes, k, self.personal_rank, m)
        )
        if self.interaction_enabled:
            assert self.interaction_basis is not None
            q = int(self.interaction_basis.n_pairs)
            rank = int(config.interaction_rank)
            self.interaction_global = nn.Parameter(torch.zeros(q, m))
            self.interaction_module_basis = nn.Parameter(
                torch.empty(q, rank, m)
            )
            self.interaction_celltype_loading = nn.Parameter(
                torch.zeros(self.n_celltypes, q, rank)
            )
            # A tiny nonzero dictionary gives the zero-initialized cell-type
            # loadings a gradient when the delayed interaction ramp opens.
            nn.init.normal_(
                self.interaction_module_basis, mean=0.0, std=1.0e-4
            )
        else:
            self.register_parameter("interaction_global", None)
            self.register_parameter("interaction_module_basis", None)
            self.register_parameter("interaction_celltype_loading", None)
        initial_interaction_scale = (
            0.0
            if self.interaction_enabled and int(config.interaction_warmup_epochs) > 0
            else 1.0
        )
        self.register_buffer(
            "interaction_scale",
            torch.tensor(initial_interaction_scale, dtype=torch.float32),
            persistent=True,
        )
        initial_explicit_scale = (
            0.0
            if int(config.output_start_epoch) > 1
            else (
                1.0
                if int(config.output_ramp_epochs) <= 0
                else 1.0 / float(config.output_ramp_epochs)
            )
        )
        self.register_buffer(
            "explicit_scale",
            torch.tensor(initial_explicit_scale, dtype=torch.float32),
            persistent=True,
        )
        # Tiny non-zero personal/response dictionaries let the context encoder
        # receive a real gradient on the very first batch while keeping the
        # warm-started CT64 score perturbation negligible.
        nn.init.normal_(self.personal_basis, mean=0.0, std=1.0e-4)
        nn.init.normal_(self.response_basis, mean=0.0, std=1.0e-4)
        # Weak adversary: code should not simply copy the five pathology labels.
        self.pathology_adversary = nn.Linear(self.personal_rank, k)
        nn.init.zeros_(self.pathology_adversary.weight)
        nn.init.zeros_(self.pathology_adversary.bias)
        self.age_adversary = nn.Linear(self.personal_rank, 1)
        nn.init.zeros_(self.age_adversary.weight)
        nn.init.zeros_(self.age_adversary.bias)
        state_adv_embed = 8
        self.state_adv_celltype_embedding = nn.Embedding(self.n_celltypes, state_adv_embed)
        self.state_pathology_adversary = nn.Sequential(
            nn.Linear(self.latent_dim + state_adv_embed, 32),
            nn.GELU(),
            nn.Linear(32, k),
        )
        nn.init.zeros_(self.state_pathology_adversary[-1].weight)
        nn.init.zeros_(self.state_pathology_adversary[-1].bias)

        if self.module_local_enabled:
            local_rank = int(config.module_local_rank)
            module_index = torch.arange(m, dtype=torch.float32).unsqueeze(0)
            rank_index = torch.arange(
                1, local_rank + 1, dtype=torch.float32
            ).unsqueeze(1)
            # Deterministic DCT rows avoid advancing the historical model RNG
            # stream while providing a semi-orthogonal shared read dictionary.
            read_init = torch.cos(
                math.pi * (module_index + 0.5) * rank_index / float(m)
            )
            read_init = F.normalize(read_init, p=2.0, dim=1)
            self.module_local_read = nn.Parameter(read_init)
            self.module_local_write_global = nn.Parameter(
                read_init.transpose(0, 1).contiguous() * 1.0e-3
            )
            self.module_local_write_celltype_delta = nn.Parameter(
                torch.zeros(self.n_celltypes, m, local_rank)
            )
            self.module_local_diagonal_global = nn.Parameter(torch.ones(m))
            self.module_local_diagonal_celltype_delta = nn.Parameter(
                torch.zeros(self.n_celltypes, m)
            )
            # Exact zero gives checkpoint-preserving no-op output.  It is a
            # signed coefficient, not a probability or stochastic hard gate.
            self.module_local_output_gate = nn.Parameter(
                torch.zeros(self.n_celltypes, m)
            )
            if self.module_local_nonlinear_enabled:
                def _logit(probability: float) -> float:
                    probability = min(max(float(probability), 1.0e-6), 1.0 - 1.0e-6)
                    return math.log(probability / (1.0 - probability))

                mix_fraction = min(
                    0.10 / float(config.module_local_compartment_mix_max),
                    0.90,
                )
                threshold_fraction = (
                    1.0 - float(config.module_local_threshold_min)
                ) / (
                    float(config.module_local_threshold_max)
                    - float(config.module_local_threshold_min)
                )
                slope_fraction = (
                    4.0 - float(config.module_local_slope_min)
                ) / (
                    float(config.module_local_slope_max)
                    - float(config.module_local_slope_min)
                )
                gain_fraction = min(
                    0.10 / float(config.module_local_gain_max), 0.90
                )
                initial_raw = {
                    "mix": _logit(mix_fraction),
                    "threshold": _logit(threshold_fraction),
                    "slope": _logit(slope_fraction),
                    "gain": _logit(gain_fraction),
                }
                trainable_families = set(
                    self.module_local_nonlinear_trainable_families
                )
                for family, value in initial_raw.items():
                    requires_grad = family in trainable_families
                    setattr(
                        self,
                        f"module_local_nonlinear_{family}_global",
                        nn.Parameter(
                            torch.full((m,), float(value)),
                            requires_grad=requires_grad,
                        ),
                    )
                    setattr(
                        self,
                        f"module_local_nonlinear_{family}_celltype_delta",
                        nn.Parameter(
                            torch.zeros(self.n_celltypes, m),
                            requires_grad=requires_grad,
                        ),
                    )
            else:
                for family in ("mix", "threshold", "slope", "gain"):
                    self.register_parameter(
                        f"module_local_nonlinear_{family}_global", None
                    )
                    self.register_parameter(
                        f"module_local_nonlinear_{family}_celltype_delta", None
                    )
            self.module_local_pathology_adversary = nn.Linear(local_rank, k)
            nn.init.zeros_(self.module_local_pathology_adversary.weight)
            nn.init.zeros_(self.module_local_pathology_adversary.bias)
            self.module_local_age_adversary = nn.Linear(local_rank, 1)
            nn.init.zeros_(self.module_local_age_adversary.weight)
            nn.init.zeros_(self.module_local_age_adversary.bias)
        else:
            for name in (
                "module_local_read",
                "module_local_write_global",
                "module_local_write_celltype_delta",
                "module_local_diagonal_global",
                "module_local_diagonal_celltype_delta",
                "module_local_output_gate",
            ):
                self.register_parameter(name, None)
            for family in ("mix", "threshold", "slope", "gain"):
                self.register_parameter(
                    f"module_local_nonlinear_{family}_global", None
                )
                self.register_parameter(
                    f"module_local_nonlinear_{family}_celltype_delta", None
                )
            self.module_local_pathology_adversary = None
            self.module_local_age_adversary = None

    def _source_mask(self, donor: torch.Tensor, target_celltype: torch.Tensor) -> torch.Tensor:
        observed = self.source_observed[donor]
        target_free = self.context_celltype.unsqueeze(0) != target_celltype.unsqueeze(1)
        return observed & target_free

    def _drop_contexts(self, mask: torch.Tensor) -> torch.Tensor:
        probability = float(self.cfg.context_input_dropout)
        if not self.training or probability <= 0.0:
            return mask
        keep = torch.rand(mask.shape, device=mask.device) >= probability
        dropped = mask & keep
        empty = ~dropped.any(dim=1)
        if bool(empty.any()):
            for row in torch.nonzero(empty, as_tuple=False).flatten():
                available = torch.nonzero(mask[row], as_tuple=False).flatten()
                if available.numel() > 0:
                    dropped[row, available[0]] = True
        return dropped

    def _encode_unique(
        self,
        donor: torch.Tensor,
        target_celltype: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        donor = donor.long().clamp(0, self.n_donors - 1)
        target_celltype = target_celltype.long().clamp(0, self.n_celltypes - 1)
        mask = self._drop_contexts(self._source_mask(donor, target_celltype))
        if bool((~mask.any(dim=1)).any()):
            bad = torch.nonzero(~mask.any(dim=1), as_tuple=False).flatten().tolist()
            raise RuntimeError(f"personal context is empty after target-celltype masking: rows={bad[:8]}")

        source_module = self.source_module[donor]
        source_latent = self.source_latent[donor]
        reliability = self.source_reliability[donor] * mask.to(
            self.source_reliability.dtype
        )
        module_token = self.module_projection(source_module)
        latent_token = (
            self.latent_projection(source_latent)
            if self.latent_projection is not None
            else None
        )
        ct_token = self.context_celltype_embedding(self.context_celltype)[None].expand(
            donor.shape[0], -1, -1
        )
        region_token = self.context_region_embedding(self.context_region)[None].expand(
            donor.shape[0], -1, -1
        )
        token_parts = [module_token]
        if latent_token is not None:
            token_parts.append(latent_token)
        token_parts.extend((ct_token, region_token))
        token = self.context_token_projection(torch.cat(token_parts, dim=-1))
        token = token * reliability.clamp_min(0.0).sqrt().unsqueeze(-1)
        encoded = self.context_encoder(token, src_key_padding_mask=~mask)
        weight = reliability.to(encoded.dtype)
        pooled = (encoded * weight.unsqueeze(-1)).sum(dim=1) / weight.sum(
            dim=1, keepdim=True
        ).clamp_min(1.0)
        code = self.personal_posterior(pooled)
        return code, mask.sum(dim=1), reliability.sum(dim=1)

    def infer_personal_code(
        self,
        donor_id: torch.Tensor,
        celltype_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Infer once per unique donor×target-celltype and expand to cells."""

        key = donor_id.long() * self.n_celltypes + celltype_id.long()
        unique, inverse = torch.unique(key, sorted=True, return_inverse=True)
        donor_unique = torch.div(unique, self.n_celltypes, rounding_mode="floor")
        celltype_unique = unique.remainder(self.n_celltypes)
        code_unique, count_unique, reliability_unique = self._encode_unique(
            donor_unique, celltype_unique
        )
        return (
            code_unique[inverse],
            count_unique[inverse],
            reliability_unique[inverse],
        )

    def _module_local_pool_unique(
        self,
        donor: torch.Tensor,
        target_celltype: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Region-balanced, modulewise pooling with an exact target firewall."""

        if not self.module_local_enabled or self.module_local_reliability is None:
            raise RuntimeError("module-local pooling requested while disabled")
        donor = donor.long().clamp(0, self.n_donors - 1)
        target_celltype = target_celltype.long().clamp(
            0, self.n_celltypes - 1
        )
        # _source_mask excludes every region belonging to the target cell type.
        mask = self._source_mask(donor, target_celltype)
        if bool((~mask.any(dim=1)).any()):
            bad = torch.nonzero(~mask.any(dim=1), as_tuple=False).flatten().tolist()
            raise RuntimeError(
                "module-local context is empty after target-celltype masking: "
                f"rows={bad[:8]}"
            )
        source = self.source_module[donor]
        assert self.module_local_input_clip is not None
        clip = self.module_local_input_clip.to(
            device=source.device, dtype=source.dtype
        )
        source = source.clamp(min=-clip, max=clip)
        weight = (
            self.source_reliability[donor].unsqueeze(-1)
            * self.module_local_reliability.unsqueeze(0).to(
                device=source.device, dtype=source.dtype
            )
            * mask.unsqueeze(-1).to(source.dtype)
        )

        region_summaries = []
        region_available = []
        for region in range(self.n_regions):
            in_region = (self.context_region == int(region)).to(source.dtype)
            region_weight = weight * in_region.view(1, -1, 1)
            denominator = region_weight.sum(dim=1)
            available = denominator > 0.0
            summary = (source * region_weight).sum(dim=1) / denominator.clamp_min(
                1.0e-8
            )
            region_summaries.append(summary)
            region_available.append(available)
        stacked_summary = torch.stack(region_summaries, dim=1)
        stacked_available = torch.stack(region_available, dim=1)
        available_float = stacked_available.to(stacked_summary.dtype)
        pooled = (
            stacked_summary * available_float
        ).sum(dim=1) / available_float.sum(dim=1).clamp_min(1.0)
        coverage = stacked_available.any(dim=1)
        pooled = torch.where(coverage, pooled, torch.zeros_like(pooled))
        return pooled, coverage

    def _project_module_local_from_personal_span(
        self,
        value: torch.Tensor,
        celltype_id: torch.Tensor,
        coverage: torch.Tensor,
    ) -> torch.Tensor:
        """Remove the detached rank-P personal span without a dense MxM matrix."""

        ct = celltype_id.long().clamp(0, self.n_celltypes - 1)
        # Keep the small normal-equation solve out of mixed precision.  Under
        # bf16 autocast, einsum can produce a bf16 RHS while the regularized
        # Gram matrix remains float32; torch.linalg.solve then rejects the
        # mismatched dtypes.  Solving in float32 is also materially safer for
        # this projection than solving in bf16.  Preserve float64 explicitly
        # for numerical/Jacobian tests and cast only the returned residual.
        solve_dtype = (
            torch.float64 if value.dtype == torch.float64 else torch.float32
        )
        with torch.autocast(device_type=value.device.type, enabled=False):
            value_solve = value.to(dtype=solve_dtype)
            basis = self.personal_basis[ct].detach().to(
                device=value.device, dtype=solve_dtype
            )
            basis = basis * coverage.unsqueeze(1).to(dtype=solve_dtype)
            gram = torch.einsum("bpm,bqm->bpq", basis, basis)
            relative_scale = gram.diagonal(dim1=1, dim2=2).mean(dim=1).clamp_min(
                1.0e-12
            )
            identity = torch.eye(
                self.personal_rank, device=value.device, dtype=solve_dtype
            ).unsqueeze(0)
            gram = gram + (
                float(self.cfg.module_local_projection_ridge)
                * relative_scale[:, None, None]
                * identity
            )
            rhs = torch.einsum(
                "bpm,bm->bp", basis, value_solve
            ).unsqueeze(-1)
            coefficient = torch.linalg.solve(gram, rhs).squeeze(-1)
            projection = torch.einsum("bp,bpm->bm", coefficient, basis)
            residual = (value_solve - projection) * coverage.to(
                dtype=solve_dtype
            )
        return residual.to(dtype=value.dtype)

    def _module_local_nonlinear_transform(
        self,
        transformed: torch.Tensor,
        celltype_id: torch.Tensor,
        coverage: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return a signed, zero-tangent threshold innovation and crossings.

        The subtraction of the gate value at zero removes the first-order
        linear term near the origin.  The branch therefore has to express
        curvature/threshold behaviour instead of merely duplicating the v1
        linear adapter with another gain vector.
        """

        if not self.module_local_nonlinear_enabled:
            return torch.zeros_like(transformed), transformed.new_zeros(
                transformed.shape[0]
            )
        assert self.module_local_compartment_adjacency is not None
        raw_values: dict[str, torch.Tensor] = {}
        for family in self.module_local_nonlinear_trainable_families:
            global_parameter = getattr(
                self, f"module_local_nonlinear_{family}_global"
            )
            delta_parameter = getattr(
                self, f"module_local_nonlinear_{family}_celltype_delta"
            )
            assert global_parameter is not None
            assert delta_parameter is not None
            centered_delta = delta_parameter - delta_parameter.mean(
                dim=0, keepdim=True
            )
            raw_values[family] = global_parameter.unsqueeze(0) + centered_delta[
                celltype_id
            ]

        unit = {family: torch.sigmoid(value) for family, value in raw_values.items()}
        mix = float(self.cfg.module_local_compartment_mix_max) * unit["mix"]
        gain = float(self.cfg.module_local_gain_max) * unit["gain"]
        graph_drive = transformed @ self.module_local_compartment_adjacency.T.to(
            device=transformed.device, dtype=transformed.dtype
        )
        drive = transformed + mix * graph_drive
        if self.module_local_nonlinear_variant == "graph_linear_control":
            # Honest graph-linear function-class control.  Only mix and gain
            # are trainable, so no redundant threshold/slope path can alter
            # optimisation geometry or impose a hidden raw-logit prior.
            innovation = gain * drive
            innovation = innovation * coverage.to(innovation.dtype)
            return innovation, transformed.new_zeros(transformed.shape[0])
        threshold = float(self.cfg.module_local_threshold_min) + (
            float(self.cfg.module_local_threshold_max)
            - float(self.cfg.module_local_threshold_min)
        ) * unit["threshold"]
        slope = float(self.cfg.module_local_slope_min) + (
            float(self.cfg.module_local_slope_max)
            - float(self.cfg.module_local_slope_min)
        ) * unit["slope"]
        gate = torch.sigmoid(slope * (drive.abs() - threshold))
        zero_gate = torch.sigmoid(-slope * threshold)
        innovation = gain * drive * (gate - zero_gate)
        innovation = innovation * coverage.to(innovation.dtype)
        crossing = (
            (drive.abs() >= threshold) & coverage
        ).to(torch.float32).sum(dim=1) / coverage.to(torch.float32).sum(
            dim=1
        ).clamp_min(1.0)
        return innovation, crossing

    def _module_local_unique(
        self,
        donor: torch.Tensor,
        target_celltype: torch.Tensor,
    ) -> tuple[
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
    ]:
        if not self.module_local_enabled:
            raise RuntimeError("module-local adapter requested while disabled")
        assert self.module_local_read is not None
        assert self.module_local_write_global is not None
        assert self.module_local_write_celltype_delta is not None
        assert self.module_local_diagonal_global is not None
        assert self.module_local_diagonal_celltype_delta is not None
        assert self.module_local_output_gate is not None
        assert self.module_local_output_cap is not None

        summary, coverage = self._module_local_pool_unique(
            donor, target_celltype
        )
        ct = target_celltype.long().clamp(0, self.n_celltypes - 1)
        code = summary @ self.module_local_read.transpose(0, 1)
        write_delta = self.module_local_write_celltype_delta - (
            self.module_local_write_celltype_delta.mean(dim=0, keepdim=True)
        )
        write = self.module_local_write_global.unsqueeze(0) + write_delta[ct]
        diagonal_delta = self.module_local_diagonal_celltype_delta - (
            self.module_local_diagonal_celltype_delta.mean(dim=0, keepdim=True)
        )
        diagonal = self.module_local_diagonal_global.unsqueeze(0) + (
            diagonal_delta[ct]
        )
        transformed = diagonal * summary + torch.einsum(
            "br,bmr->bm", code, write
        )
        nonlinear_scale = (
            self._module_local_nonlinear_scale_value
            if self.module_local_nonlinear_enabled
            else 0.0
        )
        if nonlinear_scale == 0.0:
            nonlinear_innovation = torch.zeros_like(transformed)
            threshold_crossing = transformed.new_zeros(transformed.shape[0])
            combined = transformed
            if self.module_local_nonlinear_enabled:
                # DDP uses find_unused_parameters=False.  Before the nonlinear
                # ramp opens, retain a zero-gradient path to every dormant
                # nonlinear parameter without evaluating the graph transform
                # or changing the exact v1 forward values.
                parameter_touch = transformed.new_zeros(())
                for family in self.module_local_nonlinear_trainable_families:
                    for suffix in ("global", "celltype_delta"):
                        parameter = getattr(
                            self,
                            f"module_local_nonlinear_{family}_{suffix}",
                        )
                        assert parameter is not None
                        parameter_touch = parameter_touch + 0.0 * parameter.sum()
                combined = combined + parameter_touch.to(
                    device=combined.device, dtype=combined.dtype
                )
        else:
            nonlinear_innovation, threshold_crossing = (
                self._module_local_nonlinear_transform(
                    transformed, ct, coverage
                )
            )
            combined = transformed + self.module_local_nonlinear_scale.to(
                device=transformed.device, dtype=transformed.dtype
            ) * nonlinear_innovation
        bounded = self.module_local_output_cap.to(
            device=summary.device, dtype=summary.dtype
        ) * torch.tanh(self.module_local_output_gate[ct] * combined)
        bounded = bounded * coverage.to(bounded.dtype)
        residual = self._project_module_local_from_personal_span(
            bounded, ct, coverage
        )
        if nonlinear_scale == 0.0:
            nonlinear_residual = torch.zeros_like(residual)
        else:
            linear_bounded = self.module_local_output_cap.to(
                device=summary.device, dtype=summary.dtype
            ) * torch.tanh(self.module_local_output_gate[ct] * transformed)
            linear_bounded = linear_bounded * coverage.to(linear_bounded.dtype)
            linear_residual = self._project_module_local_from_personal_span(
                linear_bounded, ct, coverage
            )
            nonlinear_residual = residual - linear_residual
        support_coverage = coverage.float().mean(dim=1)
        zero_coverage_fraction = 1.0 - support_coverage
        return (
            residual,
            code,
            support_coverage,
            zero_coverage_fraction,
            nonlinear_residual,
            threshold_crossing,
        )

    def infer_module_local(
        self,
        donor_id: torch.Tensor,
        celltype_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Infer module-local terms once per donor x target-celltype pair."""

        if not self.module_local_enabled:
            raise RuntimeError("module-local inference requested while disabled")
        key = donor_id.long() * self.n_celltypes + celltype_id.long()
        unique, inverse = torch.unique(key, sorted=True, return_inverse=True)
        donor_unique = torch.div(unique, self.n_celltypes, rounding_mode="floor")
        celltype_unique = unique.remainder(self.n_celltypes)
        raw, code, coverage, zero_fraction, _, _ = self._module_local_unique(
            donor_unique, celltype_unique
        )
        return (
            raw[inverse],
            code[inverse],
            coverage[inverse],
            zero_fraction[inverse],
        )

    def infer_module_local_with_diagnostics(
        self,
        donor_id: torch.Tensor,
        celltype_id: torch.Tensor,
    ) -> tuple[
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
    ]:
        """Infer the public v1 outputs plus nonlinear innovation diagnostics."""

        if not self.module_local_enabled:
            raise RuntimeError("module-local inference requested while disabled")
        key = donor_id.long() * self.n_celltypes + celltype_id.long()
        unique, inverse = torch.unique(key, sorted=True, return_inverse=True)
        donor_unique = torch.div(unique, self.n_celltypes, rounding_mode="floor")
        celltype_unique = unique.remainder(self.n_celltypes)
        values = self._module_local_unique(donor_unique, celltype_unique)
        return tuple(value[inverse] for value in values)  # type: ignore[return-value]

    def personal_terms_from_code(
        self,
        code: torch.Tensor,
        celltype_id: torch.Tensor,
        pathology: torch.Tensor,
        pathology_valid: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Map an own- or other-donor code into personal module terms."""

        if code.ndim != 2 or code.shape[1] != self.personal_rank:
            raise ValueError("code must have shape [B,personal_rank]")
        ct = celltype_id.long().clamp(0, self.n_celltypes - 1)
        p = pathology.float().clamp(0.0, 1.0)
        valid = pathology_valid.bool()
        if p.shape != (code.shape[0], self.n_pathology) or valid.shape != p.shape:
            raise ValueError("pathology/value mask must both be [B,5]")
        p = torch.where(valid, p, torch.zeros_like(p))
        personal = torch.einsum("bp,bpm->bm", code, self.personal_basis[ct])
        response = torch.einsum(
            "bk,bp,bkpm->bkm", p, code, self.response_basis[ct]
        )
        response = response * valid.unsqueeze(-1).to(response.dtype)
        return personal, response

    @torch.no_grad()
    def select_matched_negative_donors(
        self,
        donor_id: torch.Tensor,
        target_celltype: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Choose exact-support-set, nearest-covariate training donors."""

        if not self.support_query_enabled:
            raise RuntimeError("support-query matching requested while disabled")
        donor = donor_id.long().clamp(0, self.n_donors - 1)
        target = target_celltype.long().clamp(0, self.n_celltypes - 1)
        batch = donor.shape[0]
        candidates = torch.arange(self.n_donors, device=donor.device)

        target_free = self.context_celltype.unsqueeze(0) != target.unsqueeze(1)
        own_support = self.source_observed[donor] & target_free
        candidate_support = self.source_observed.unsqueeze(0) & target_free.unsqueeze(1)
        exact_support = (candidate_support == own_support.unsqueeze(1)).all(dim=2)

        valid = exact_support & self.match_fit_donor.unsqueeze(0)
        valid &= candidates.unsqueeze(0) != donor.unsqueeze(1)
        own_validity = self.match_pathology_valid[donor]
        valid &= (
            self.match_pathology_valid.unsqueeze(0)
            == own_validity.unsqueeze(1)
        ).all(dim=2)
        own_sex = self.match_sex[donor]
        sex_known = torch.isfinite(own_sex).unsqueeze(1) & torch.isfinite(
            self.match_sex
        ).unsqueeze(0)
        same_sex = (~sex_known) | (
            (self.match_sex.unsqueeze(0) - own_sex.unsqueeze(1)).abs() < 0.5
        )
        valid &= same_sex

        path_delta = (
            self.match_pathology.unsqueeze(0)
            - self.match_pathology[donor].unsqueeze(1)
        ).abs().mean(dim=2)
        age_delta = (
            self.match_age_z.unsqueeze(0) - self.match_age_z[donor].unsqueeze(1)
        ).abs()
        adnc_delta = (
            self.match_adnc.unsqueeze(0) - self.match_adnc[donor].unsqueeze(1)
        ).abs()
        own_rel = self.source_reliability[donor] * own_support.to(
            self.source_reliability.dtype
        )
        candidate_rel = self.source_reliability.unsqueeze(0) * candidate_support.to(
            self.source_reliability.dtype
        )
        support_denom = own_support.sum(dim=1, keepdim=True).clamp_min(1)
        support_delta = (
            (candidate_rel - own_rel.unsqueeze(1)).abs()
            * own_support.unsqueeze(1).to(candidate_rel.dtype)
        ).sum(dim=2) / support_denom.to(candidate_rel.dtype)

        distance = path_delta + 0.5 * adnc_delta + 0.25 * age_delta + support_delta
        distance = distance.masked_fill(~valid, torch.inf)
        k = int(self.cfg.support_query_n_negatives)
        top_distance, top_donor = torch.topk(
            distance, k=min(k, self.n_donors), dim=1, largest=False
        )
        top_valid = torch.isfinite(top_distance)
        if top_donor.shape[1] < k:
            pad = k - top_donor.shape[1]
            top_donor = torch.cat(
                [top_donor, donor.unsqueeze(1).expand(batch, pad)], dim=1
            )
            top_distance = torch.cat(
                [top_distance, torch.full_like(top_distance[:, :pad], torch.inf)],
                dim=1,
            )
            top_valid = torch.cat(
                [top_valid, torch.zeros_like(top_valid[:, :pad])], dim=1
            )
        top_donor = torch.where(top_valid, top_donor, donor.unsqueeze(1))
        return top_donor, top_valid, top_distance

    @torch.no_grad()
    def set_interaction_epoch(self, epoch: int) -> float:
        """Set the deterministic interaction warmup/ramp for a 1-based epoch."""

        if isinstance(epoch, bool) or int(epoch) <= 0:
            raise ValueError("epoch must be a positive 1-based integer")
        if not self.interaction_enabled:
            self.interaction_scale.fill_(1.0)
            return 1.0
        warmup = int(self.cfg.interaction_warmup_epochs)
        ramp_epochs = int(self.cfg.interaction_ramp_epochs)
        if int(epoch) <= warmup:
            scale = 0.0
        elif ramp_epochs <= 0:
            scale = 1.0
        else:
            scale = min(1.0, (int(epoch) - warmup) / float(ramp_epochs))
        self.interaction_scale.fill_(float(scale))
        return float(scale)

    @torch.no_grad()
    def set_curriculum_epoch(self, epoch: int) -> float:
        """Set the explicit branch's deterministic output scale."""

        scale = delayed_epoch_ramp(
            int(epoch),
            start_epoch=int(self.cfg.output_start_epoch),
            ramp_epochs=int(self.cfg.output_ramp_epochs),
        )
        self.explicit_scale.fill_(float(scale))
        return float(scale)

    @torch.no_grad()
    def set_module_local_epoch(self, epoch: int) -> float:
        """Set and persist the module-local output curriculum for one epoch."""

        if not self.module_local_enabled:
            return 0.0
        assert self.module_local_scale is not None
        scale = delayed_epoch_ramp(
            int(epoch),
            start_epoch=int(self.cfg.module_local_start_epoch),
            ramp_epochs=int(self.cfg.module_local_ramp_epochs),
        )
        self.module_local_scale.fill_(float(scale))
        return float(scale)

    @torch.no_grad()
    def set_module_local_nonlinear_epoch(self, epoch: int) -> float:
        """Set and persist the delayed compartment-nonlinearity scale."""

        if not self.module_local_nonlinear_enabled:
            return 0.0
        assert self.module_local_nonlinear_scale is not None
        scale = delayed_epoch_ramp(
            int(epoch),
            start_epoch=int(self.cfg.module_local_nonlinear_start_epoch),
            ramp_epochs=int(self.cfg.module_local_nonlinear_ramp_epochs),
        )
        self.module_local_nonlinear_scale.fill_(float(scale))
        self._module_local_nonlinear_scale_value = float(scale)
        return float(scale)

    def _load_from_state_dict(
        self,
        state_dict,
        prefix,
        local_metadata,
        strict,
        missing_keys,
        unexpected_keys,
        error_msgs,
    ) -> None:
        """Restore the Python curriculum mirror once per checkpoint load."""

        super()._load_from_state_dict(
            state_dict,
            prefix,
            local_metadata,
            strict,
            missing_keys,
            unexpected_keys,
            error_msgs,
        )
        if (
            self.module_local_nonlinear_enabled
            and self.module_local_nonlinear_scale is not None
        ):
            self._module_local_nonlinear_scale_value = float(
                self.module_local_nonlinear_scale.detach().item()
            )

    @torch.no_grad()
    def module_local_nonlinear_parameter_diagnostics(
        self,
    ) -> dict[str, torch.Tensor]:
        """Return bounded parameter summaries without changing optimisation."""

        if not self.module_local_nonlinear_enabled:
            return {}

        def _unit(family: str) -> torch.Tensor:
            global_parameter = getattr(
                self, f"module_local_nonlinear_{family}_global"
            )
            delta_parameter = getattr(
                self, f"module_local_nonlinear_{family}_celltype_delta"
            )
            assert global_parameter is not None
            assert delta_parameter is not None
            centered = delta_parameter - delta_parameter.mean(
                dim=0, keepdim=True
            )
            return torch.sigmoid(
                global_parameter.unsqueeze(0) + centered
            )

        gain = float(self.cfg.module_local_gain_max) * _unit("gain")
        mix = float(self.cfg.module_local_compartment_mix_max) * _unit("mix")
        diagnostics = {
            "gain_mean": gain.mean(),
            "gain_rms": gain.square().mean().sqrt(),
            "gain_min": gain.min(),
            "gain_max": gain.max(),
            "mix_mean": mix.mean(),
        }
        if self.module_local_nonlinear_variant == "graph_linear_control":
            diagnostics.update(
                {
                    "linear_scale_mean": gain.mean(),
                    "linear_scale_rms": gain.square().mean().sqrt(),
                    "linear_scale_min": gain.min(),
                    "linear_scale_max": gain.max(),
                }
            )
        else:
            threshold = float(self.cfg.module_local_threshold_min) + (
                float(self.cfg.module_local_threshold_max)
                - float(self.cfg.module_local_threshold_min)
            ) * _unit("threshold")
            slope = float(self.cfg.module_local_slope_min) + (
                float(self.cfg.module_local_slope_max)
                - float(self.cfg.module_local_slope_min)
            ) * _unit("slope")
            diagnostics.update(
                {
                    "threshold_mean": threshold.mean(),
                    "threshold_min": threshold.min(),
                    "threshold_max": threshold.max(),
                    "slope_mean": slope.mean(),
                    "slope_min": slope.min(),
                    "slope_max": slope.max(),
                }
            )
        return diagnostics

    def forward(
        self,
        *,
        donor_id: torch.Tensor,
        celltype_id: torch.Tensor,
        region_id: torch.Tensor,
        age_z: torch.Tensor,
        age_valid: torch.Tensor,
        pathology: torch.Tensor,
        pathology_valid: torch.Tensor,
        state_latent: Optional[torch.Tensor] = None,
        cell_weight: Optional[torch.Tensor] = None,
    ) -> PrecisionMedicineOutput:
        pathology = pathology.float().clamp(0.0, 1.0)
        pathology_valid = pathology_valid.bool()
        if pathology.shape != pathology_valid.shape or pathology.shape[1] != self.n_pathology:
            raise ValueError("precision pathology/value mask must both be [B,5]")
        # Missing is *not* normal.  It makes no contribution and remains marked
        # invalid for the leakage loss and downstream fingerprint audit.
        p = torch.where(pathology_valid, pathology, torch.zeros_like(pathology))
        features = pathology_hinge_features(p) * pathology_valid.unsqueeze(-1).to(p.dtype)

        code, support_count, support_reliability = self.infer_personal_code(
            donor_id, celltype_id
        )
        ct = celltype_id.long().clamp(0, self.n_celltypes - 1)
        region = region_id.long().clamp(0, self.n_regions - 1)

        normal_region_basis = self.normal_region_delta - self.normal_region_delta.mean(
            dim=1, keepdim=True
        )
        normal_region_coeff = normal_region_basis[ct, region]
        age_value = torch.where(
            age_valid.bool().reshape(-1),
            age_z.float().reshape(-1),
            torch.zeros_like(age_z.float().reshape(-1)),
        )
        age_basis = self.age_global.unsqueeze(0) + (
            self.age_celltype_delta
            - self.age_celltype_delta.mean(dim=0, keepdim=True)
        )[ct]
        age_coeff = age_value.unsqueeze(-1) * age_basis

        ct_delta = self.common_celltype_delta - self.common_celltype_delta.mean(
            dim=0, keepdim=True
        )
        common_basis = self.common_global.unsqueeze(0) + ct_delta[ct]
        region_gate = self.common_region_gate - self.common_region_gate.mean(
            dim=0, keepdim=True
        )
        gate = 1.0 + float(self.cfg.region_gate_scale) * torch.tanh(region_gate[region])
        common_axis = features.unsqueeze(-1) * common_basis * gate.unsqueeze(-1)
        common_axis_coeff = common_axis.sum(dim=2)

        personal_rank2_coeff, response_axis_coeff = self.personal_terms_from_code(
            code,
            ct,
            p,
            pathology_valid,
        )
        if self.module_local_enabled:
            (
                module_local_raw_coeff,
                module_local_code,
                module_local_support_coverage,
                module_local_zero_coverage_fraction,
                module_local_nonlinear_raw_coeff,
                module_local_threshold_crossing_fraction,
            ) = self.infer_module_local_with_diagnostics(donor_id, ct)
            assert self.module_local_scale is not None
            module_local_coeff = module_local_raw_coeff * self.module_local_scale.to(
                device=module_local_raw_coeff.device,
                dtype=module_local_raw_coeff.dtype,
            )
            module_local_nonlinear_coeff = (
                module_local_nonlinear_raw_coeff
                * self.module_local_scale.to(
                    device=module_local_nonlinear_raw_coeff.device,
                    dtype=module_local_nonlinear_raw_coeff.dtype,
                )
            )
            detached_personal_basis = self.personal_basis[ct].detach()
            span_overlap = torch.einsum(
                "bpm,bm->bp",
                detached_personal_basis,
                module_local_raw_coeff,
            ).abs()
            span_denominator = (
                detached_personal_basis.norm(dim=2)
                * module_local_raw_coeff.norm(dim=1, keepdim=True)
            )
            normalized_overlap = span_overlap / span_denominator.clamp_min(
                1.0e-12
            )
            normalized_overlap = torch.where(
                span_denominator > 0.0,
                normalized_overlap,
                torch.zeros_like(normalized_overlap),
            )
            module_local_personal_span_overlap = normalized_overlap.max(dim=1).values
        else:
            module_local_raw_coeff = personal_rank2_coeff.new_zeros(
                personal_rank2_coeff.shape
            )
            module_local_coeff = module_local_raw_coeff
            module_local_code = personal_rank2_coeff.new_zeros(
                personal_rank2_coeff.shape[0], 0
            )
            module_local_personal_span_overlap = personal_rank2_coeff.new_zeros(
                personal_rank2_coeff.shape[0]
            )
            module_local_support_coverage = personal_rank2_coeff.new_zeros(
                personal_rank2_coeff.shape[0]
            )
            module_local_zero_coverage_fraction = personal_rank2_coeff.new_ones(
                personal_rank2_coeff.shape[0]
            )
            module_local_nonlinear_raw_coeff = personal_rank2_coeff.new_zeros(
                personal_rank2_coeff.shape
            )
            module_local_nonlinear_coeff = module_local_nonlinear_raw_coeff
            module_local_threshold_crossing_fraction = (
                personal_rank2_coeff.new_zeros(personal_rank2_coeff.shape[0])
            )
        personal_coeff = personal_rank2_coeff + module_local_coeff
        if self.interaction_enabled:
            assert self.interaction_basis is not None
            assert self.interaction_global is not None
            assert self.interaction_module_basis is not None
            assert self.interaction_celltype_loading is not None
            interaction_feature, interaction_valid = self.interaction_basis(
                p, pathology_valid
            )
            loading = self.interaction_celltype_loading - (
                self.interaction_celltype_loading.mean(dim=0, keepdim=True)
            )
            celltype_interaction = torch.einsum(
                "bqr,qrm->bqm",
                loading[ct],
                self.interaction_module_basis,
            )
            interaction_dictionary = (
                self.interaction_global.unsqueeze(0) + celltype_interaction
            )
            interaction_pair_coeff = (
                self.interaction_scale.to(
                    device=interaction_feature.device,
                    dtype=interaction_feature.dtype,
                )
                * interaction_feature.unsqueeze(-1)
                * interaction_dictionary
            )
        else:
            interaction_pair_coeff = p.new_zeros(
                p.shape[0], 0, self.n_modules
            )
            interaction_feature = p.new_zeros(p.shape[0], 0)
            interaction_valid = pathology_valid.new_zeros(
                pathology_valid.shape[0], 0
            )
        total_module_coeff_unscaled = (
            normal_region_coeff
            + age_coeff
            + common_axis_coeff.sum(dim=1)
            + personal_coeff
            + response_axis_coeff.sum(dim=1)
            + interaction_pair_coeff.sum(dim=1)
        )
        total_module_coeff = total_module_coeff_unscaled * self.explicit_scale.to(
            device=total_module_coeff_unscaled.device,
            dtype=total_module_coeff_unscaled.dtype,
        )

        if cell_weight is None:
            cell_weight = torch.ones(
                pathology.shape[0], device=pathology.device, dtype=pathology.dtype
            )
        cell_weight = cell_weight.float().reshape(-1)
        if cell_weight.shape[0] != pathology.shape[0]:
            raise ValueError("cell_weight must have shape [B]")

        leak_prediction = torch.sigmoid(self.pathology_adversary(_grad_reverse(code)))
        valid_float = pathology_valid.to(leak_prediction.dtype)
        pathology_leak_per_cell = (
            ((leak_prediction - pathology).square() * valid_float).sum(dim=1)
            / valid_float.sum(dim=1).clamp_min(1.0)
        )
        pathology_leak_loss = (pathology_leak_per_cell * cell_weight).mean()
        age_prediction = self.age_adversary(_grad_reverse(code)).squeeze(-1)
        age_valid_float = age_valid.bool().reshape(-1).to(age_prediction.dtype)
        age_leak_loss = (
            (age_prediction - age_z.float().reshape(-1)).square()
            * age_valid_float
            * cell_weight
        ).mean()

        if state_latent is None:
            state_pathology_leak_loss = pathology_leak_loss.detach() * 0.0
        else:
            if state_latent.ndim != 2 or state_latent.shape != (
                pathology.shape[0], self.latent_dim
            ):
                raise ValueError(
                    f"state_latent must be [B,{self.latent_dim}], got "
                    f"{tuple(state_latent.shape)}"
                )
            state_predictor_input = torch.cat(
                (
                    _grad_reverse(state_latent),
                    self.state_adv_celltype_embedding(ct),
                ),
                dim=1,
            )
            state_prediction = torch.sigmoid(
                self.state_pathology_adversary(state_predictor_input)
            )
            state_pathology_leak_per_cell = (
                ((state_prediction - pathology).square() * valid_float).sum(dim=1)
                / valid_float.sum(dim=1).clamp_min(1.0)
            )
            state_pathology_leak_loss = (
                state_pathology_leak_per_cell * cell_weight
            ).mean()

        local_zero = pathology_leak_loss.detach() * 0.0
        if self.module_local_enabled:
            assert self.module_local_pathology_adversary is not None
            assert self.module_local_age_adversary is not None
            assert self.module_local_scale is not None
            local_guard_scale = self.module_local_scale.to(
                device=module_local_code.device, dtype=module_local_code.dtype
            )
            local_pathology_prediction = torch.sigmoid(
                self.module_local_pathology_adversary(
                    _grad_reverse(module_local_code)
                )
            )
            local_pathology_per_cell = (
                (
                    (local_pathology_prediction - pathology).square()
                    * valid_float
                ).sum(dim=1)
                / valid_float.sum(dim=1).clamp_min(1.0)
            )
            module_local_pathology_leak_loss = local_guard_scale * (
                local_pathology_per_cell * cell_weight
            ).mean()
            local_age_prediction = self.module_local_age_adversary(
                _grad_reverse(module_local_code)
            ).squeeze(-1)
            module_local_age_leak_loss = local_guard_scale * (
                (local_age_prediction - age_z.float().reshape(-1)).square()
                * age_valid_float
                * cell_weight
            ).mean()

            center_terms = []
            for target in torch.unique(ct, sorted=True):
                selected = ct == target
                selected_weight = cell_weight[selected]
                weighted_mean = (
                    module_local_coeff[selected]
                    * selected_weight.unsqueeze(1)
                ).sum(dim=0) / selected_weight.sum().clamp_min(1.0e-8)
                center_terms.append(weighted_mean.square().mean())
            module_local_center_loss = (
                torch.stack(center_terms).mean() if center_terms else local_zero
            )
            denominator = cell_weight.sum().clamp_min(1.0e-8)
            # The module-local path is exactly zero before its curriculum
            # opens and while its output gate is initialized at zero.  A raw
            # sqrt at that point has an infinite derivative; multiplying the
            # resulting inactive size penalty by zero can still seed NaN
            # gradients throughout the dormant adapter.  Clamp the two powers
            # before sqrt so the inactive branch has a finite zero gradient.
            local_power = (
                module_local_coeff.square().mean(dim=1) * cell_weight
            ).sum() / denominator
            rank2_power = (
                personal_rank2_coeff.square().mean(dim=1) * cell_weight
            ).sum() / denominator
            local_rms = local_power.clamp_min(1.0e-12).sqrt()
            rank2_rms = rank2_power.clamp_min(1.0e-12).sqrt()
            local_to_rank2_ratio = local_rms / rank2_rms.clamp_min(1.0e-8)
            module_local_size_loss = F.relu(
                local_to_rank2_ratio
                - float(self.cfg.module_local_size_ratio_cap)
            ).square()
        else:
            module_local_center_loss = local_zero
            module_local_size_loss = local_zero
            module_local_pathology_leak_loss = local_zero
            module_local_age_leak_loss = local_zero

        return PrecisionMedicineOutput(
            normal_region_coeff=normal_region_coeff,
            age_coeff=age_coeff,
            common_axis_coeff=common_axis_coeff,
            personal_coeff=personal_coeff,
            personal_rank2_coeff=personal_rank2_coeff,
            module_local_coeff=module_local_coeff,
            module_local_raw_coeff=module_local_raw_coeff,
            module_local_code=module_local_code,
            module_local_personal_span_overlap=(
                module_local_personal_span_overlap
            ),
            module_local_support_coverage=module_local_support_coverage,
            module_local_zero_coverage_fraction=(
                module_local_zero_coverage_fraction
            ),
            module_local_nonlinear_coeff=module_local_nonlinear_coeff,
            module_local_nonlinear_raw_coeff=(
                module_local_nonlinear_raw_coeff
            ),
            module_local_threshold_crossing_fraction=(
                module_local_threshold_crossing_fraction
            ),
            response_axis_coeff=response_axis_coeff,
            interaction_pair_coeff=interaction_pair_coeff,
            interaction_feature=interaction_feature,
            interaction_valid=interaction_valid,
            total_module_coeff=total_module_coeff,
            personal_code=code,
            support_count=support_count,
            support_reliability=support_reliability,
            pathology=p,
            pathology_valid=pathology_valid,
            pathology_leak_loss=pathology_leak_loss,
            age_leak_loss=age_leak_loss,
            state_pathology_leak_loss=state_pathology_leak_loss,
            module_local_center_loss=module_local_center_loss,
            module_local_size_loss=module_local_size_loss,
            module_local_pathology_leak_loss=(
                module_local_pathology_leak_loss
            ),
            module_local_age_leak_loss=module_local_age_leak_loss,
        )

    def parameter_penalties(self) -> dict[str, torch.Tensor]:
        """Scale-stable regularisers used by the trainer."""

        normal_context = (
            self.normal_region_delta
            - self.normal_region_delta.mean(dim=1, keepdim=True)
        ).square().mean()
        normal_context = normal_context + self.age_global.square().mean() + (
            self.age_celltype_delta
            - self.age_celltype_delta.mean(dim=0, keepdim=True)
        ).square().mean()
        common = self.common_global.square().mean() + (
            self.common_celltype_delta - self.common_celltype_delta.mean(dim=0, keepdim=True)
        ).square().mean()
        common = common + self.common_region_gate.square().mean()
        if self.interaction_enabled:
            assert self.interaction_global is not None
            assert self.interaction_module_basis is not None
            assert self.interaction_celltype_loading is not None
            interaction = self.interaction_global.square().mean()
            interaction = interaction + self.interaction_module_basis.square().mean()
            interaction = interaction + (
                self.interaction_celltype_loading
                - self.interaction_celltype_loading.mean(dim=0, keepdim=True)
            ).square().mean()
        else:
            interaction = common.detach() * 0.0
        if self.module_local_enabled:
            assert self.module_local_write_celltype_delta is not None
            assert self.module_local_diagonal_celltype_delta is not None
            write_delta = self.module_local_write_celltype_delta - (
                self.module_local_write_celltype_delta.mean(dim=0, keepdim=True)
            )
            diagonal_delta = self.module_local_diagonal_celltype_delta - (
                self.module_local_diagonal_celltype_delta.mean(
                    dim=0, keepdim=True
                )
            )
            module_local_hierarchy = (
                write_delta.square().mean() + diagonal_delta.square().mean()
            )
            if self.module_local_nonlinear_enabled:
                for family in self.module_local_nonlinear_trainable_families:
                    delta = getattr(
                        self,
                        f"module_local_nonlinear_{family}_celltype_delta",
                    )
                    assert delta is not None
                    centered = delta - delta.mean(dim=0, keepdim=True)
                    module_local_hierarchy = (
                        module_local_hierarchy + centered.square().mean()
                    )
        else:
            module_local_hierarchy = common.detach() * 0.0
        return {
            "normal_context": normal_context,
            "common": common,
            "personal": self.personal_basis.square().mean(),
            "response": self.response_basis.square().mean(),
            "interaction": interaction,
            "module_local_hierarchy": module_local_hierarchy,
        }

    @torch.no_grad()
    def branch_parameter_count(self) -> dict[str, int]:
        groups = {
            "context": 0,
            "normal_context": 0,
            "common": 0,
            "personal": 0,
            "response": 0,
            "interaction": 0,
            "module_local": 0,
            "adversary": 0,
        }
        for name, parameter in self.named_parameters():
            if name.startswith((
                "pathology_adversary",
                "age_adversary",
                "state_pathology_adversary",
                "state_adv_celltype_embedding",
                "module_local_pathology_adversary",
                "module_local_age_adversary",
            )):
                group = "adversary"
            elif name.startswith("module_local_"):
                group = "module_local"
            elif name.startswith((
                "normal_region_",
                "age_global",
                "age_celltype_delta",
            )):
                group = "normal_context"
            elif name.startswith("common_"):
                group = "common"
            elif name.startswith("personal_basis"):
                group = "personal"
            elif name.startswith("response_basis"):
                group = "response"
            elif name.startswith("interaction_"):
                group = "interaction"
            else:
                group = "context"
            groups[group] += int(parameter.numel())
        groups["total"] = int(sum(groups.values()))
        return groups
