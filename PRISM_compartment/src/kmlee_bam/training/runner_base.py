from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional
from contextlib import contextmanager
from datetime import timedelta

import numpy as np
import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler

from kmlee_bam.data.ordinal_dataset import OrdinalScDataset
from kmlee_bam.modules.gene_embedding import GeneExpressionEmbedding
from kmlee_bam.model.state_encoder import StateEncoder
from kmlee_bam.model.celltype_prior import CellTypePrior
from kmlee_bam.modules.module_tokenizer import GeneModuleTokenizer
from kmlee_bam.model.lie_ordinal_decoder import LieActionOrdinalDecoder
from kmlee_bam.objectives.total_loss import TotalLoss
from kmlee_bam.objectives.pathology_aux import PathologyAuxConfig, PathologyAuxHead
from kmlee_bam.objectives.cp_bam import CPBAMConfig, CPBAMHead
from kmlee_bam.objectives.metric_contrast import MetricContrastConfig, MetricContrastHead
from kmlee_bam.model.latent_nuisance_projection import (
    NuisanceProjectionConfig,
    LatentNuisanceProjector,
)
from kmlee_bam.model.precision_medicine import PrecisionMedicineConfig, PrecisionMedicineHead
from kmlee_bam.model.pathology_interactions import (
    PairwisePathologyInteractionBasis,
    fit_pairwise_pathology_residualizer,
)
from kmlee_bam.data.donor_context_dataset import DonorContextTable
from kmlee_bam.data.module_local_reliability import (
    load_module_local_compartment_graph_artifact,
    load_module_local_reliability_artifact,
    load_module_local_output_cap_artifact,
    sha256_file as module_local_sha256_file,
)
from kmlee_bam.training.generator_budget_search import (
    GeneratorBudgetSearch,
    GeneratorBudgetSearchConfig,
)
from kmlee_bam.training.learned_generator_count import (
    ConstrainedGeneratorCountObjective,
    HardBinaryConcreteGeneratorGate,
    LearnedGeneratorCountConfig,
)
from kmlee_bam.training.learned_pathology_rank import (
    LearnedPathologyRankConfig,
    LearnedPathologyRankGate,
    algebraic_pathology_rank_capacity,
)
from kmlee_bam.training.architecture_capacity_logging import (
    build_architecture_capacity_record,
    format_architecture_capacity_record,
)
from kmlee_bam.training.architecture_donor_split import (
    DonorIndexView,
    make_architecture_donor_split,
)
from kmlee_bam.training.prism_module_rescue_training import (
    PrismModuleRescueTrainingConfig,
    PrismModuleRescueUpdater,
    write_module_rescue_manifest,
)
from kmlee_bam.training.integrated_phase_curriculum import (
    IntegratedPhaseCurriculumConfig,
    IntegratedPhaseCurriculumController,
    pathology_correction_rank_diagnostics,
)
from kmlee_bam.training.core_trainer import (
    OrdinalBAMSystem,
    Trainer,
    EpochOutput,
    finalize_prism_epoch_metrics,
)


# =====================================================================
# Config dataclasses
# =====================================================================
@dataclass
class DataConfig:
    # Zarr-backed final gene-union dataset.
    zarr_path: str
    spec_path: str
    matrix_path: str = "X"          # "X", "raw/X", "layers/counts", or "auto"
    counts_layer: str = "counts"
    prefer_raw: bool = True

    # obs schema
    split_key: str = "split"
    celltype_key: str = "Subclass"
    batch_key: str = "tech_id_model"

    # returned tensors
    return_log_input: bool = True
    reference_scaler_path: Optional[str] = None
    return_x_gene_scalar: bool = True
    x_gene_scalar_clip: Optional[float] = 10.0
    x_gene_scalar_eps: float = 1e-4

    # sparse-row extraction controls
    max_csr_data_span: int = 50_000_000
    max_span_to_selected_nnz_ratio: float = 20.0

    # loader controls
    num_workers: int = 0
    pin_memory: bool = True
    persistent_workers: bool = False
    drop_last_train: bool = False


@dataclass
class LoaderConfig:
    batch_size: int = 64
    eval_batch_size: Optional[int] = None
    shuffle_train: bool = True


@dataclass
class EmbeddingConfig:
    d_model: int = 128
    dropout: float = 0.1
    use_continuous: bool = True
    use_cls: bool = True
    cont_hidden_dim: Optional[int] = None
    layer_norm_eps: float = 1e-5
    init_std: float = 0.02


@dataclass
class ModuleTokenizerConfig:
    registry_json_path: str

    # Optional external activity weights. If null, tokenizer uses
    # L2-normalised positive membership by default.
    activity_weight_path: Optional[str] = None
    activity_weight_key: str = "activity_weight"

    pooling: str = "mean"               # "mean" or "attention"
    preserve_cls: bool = True
    use_module_id_embedding: bool = True
    module_id_init_scale: float = 0.01
    module_id_scale_min: float = 1e-4
    module_id_scale_max: float = 10.0
    use_activity_projection: bool = True
    activity_weight_normalization: str = "l2"
    activity_default: str = "l2_membership"
    dropout: float = 0.1
    layer_norm_eps: float = 1e-5
    init_std: float = 0.02

    # Decoder generator masks are usually the same membership rows, but
    # this can be too large if many modules exist. If decoder.n_generators is
    # null, all modules are used as Lie generators.
    decoder_mask_strategy: str = "all"   # "all", "first", or "largest"


@dataclass
class EncoderConfig:
    n_heads: int = 8
    n_layers: int = 4
    d_z: int = 32
    d_ff: Optional[int] = None
    condition_on_celltype: bool = True
    celltype_embed_dim: Optional[int] = None
    posterior_input_norm: Optional[bool] = None
    pooling: str = "cls"
    use_cls_token: bool = True
    mean_exclude_cls: bool = True
    agp_heads: int = 8
    agp_exclude_cls: bool = True
    # "legacy" preserves the original scalar-score + 8D->D combine readout.
    # "diverse_slots_v2" uses learned query/key/value slots, direct D->D
    # concatenation, bounded mean residual, and anti-collapse penalties.
    agp_variant: str = "legacy"
    agp_temperature_init: float = 0.75
    agp_temperature_min: float = 0.25
    agp_temperature_max: float = 2.0
    agp_mean_residual_init: float = 0.25
    agp_mean_residual_max: float = 0.50
    agp_diversity_weight: float = 0.0
    agp_diversity_margin: float = 0.50
    agp_query_orthogonality_weight: float = 0.0
    agp_entropy_band_weight: float = 0.0
    agp_min_effective_tokens: float = 1.0
    agp_max_effective_tokens: float = 1.0e9
    # Persistent-AGP keeps specialist slots across tokenizer/L1/L2/L3 reads.
    # ``agp_heads`` remains the specialist-slot count for this pooling mode.
    persistent_agp_attention_heads: int = 4
    persistent_agp_ffn_dim: Optional[int] = None
    persistent_agp_update_gate_init: float = 0.10
    persistent_agp_update_gate_max: float = 1.0
    compute_cell_uncertainty: bool = True
    uncertainty_exclude_cls: bool = True
    tech_risk_head: bool = False
    posterior_hidden_dim: Optional[int] = None
    posterior_dropout: float = 0.1
    logvar_min: float = -8.0
    logvar_max: float = 8.0
    stochastic_attention: bool = True
    distribution: str = "lognormal"
    sigma_mode: str = "global"
    sigma: float = 0.3
    sigma_min: float = 1e-4
    sigma_max: float = 1.5
    weibull_k: float = 20.0
    prior_d_mid: int = 16
    attn_dropout: float = 0.1
    proj_dropout: float = 0.1
    ffn_dropout: float = 0.1
    qkv_bias: bool = False
    layer_norm_eps: float = 1e-5
    init_std: float = 0.02
    final_norm: bool = True


@dataclass
class PriorConfig:
    sigma_init: float = 0.8
    sigma_min: float = 1e-4
    residual_eps: float = 1e-5
    residual_clamp: float = 10.0


@dataclass
class DecoderConfig:
    # If null, uses all selected module masks as Lie generators.
    n_generators: Optional[int] = None
    # Biological-covariate factorization: number of rows for the optional SEX
    # baseline embedding. None/0 disables the term (byte-identical). Use 3 for
    # {unknown, F, M} (raw sex_id -1/0/1 is shifted +1 inside the decoder).
    n_sex: Optional[int] = None
    coeff_hidden_dim: Optional[int] = None
    coeff_dropout: float = 0.1
    gamma_scale: float = 0.1
    generator_init_scale: float = 1e-3
    use_module_masks: bool = True
    # Fail closed unless every registry module has exactly one decoder
    # generator.  This prevents a mutable/larger registry from being silently
    # truncated by n_generators.
    require_generator_per_module: bool = False
    # A size-1 module has operator rank at most one. None retains historical
    # uniform-rank behaviour.
    singleton_lie_rank: Optional[int] = None
    # "none" (historical) or "size_overlap" (input size and output overlap
    # variance stabilization).
    module_generator_normalization: str = "none"

    center_thresholds: bool = True
    threshold_init_spacing: float = 1.0
    prob_eps: float = 1e-8
    init_std: float = 0.02
    translation_init_scale: float = 1e-2
    use_affine_translation: bool = True
    use_direct_state_head: bool = False
    direct_state_init_scale: float = 1e-2
    # Threshold initialisation. "uniform" (default) uses
    # threshold_init_spacing as before. "data_quantile" computes per-gene
    # marginal bin frequencies from training data and initialises
    # thresholds so the decoder reproduces them at score=0. When set,
    # center_thresholds is auto-disabled inside the decoder.
    threshold_init_mode: str = "uniform"
    bin_marginal_n_sample: int = 4096
    bin_marginal_seed: int = 42
    bin_marginal_smoothing: float = 1.0
    bin_marginal_cache_path: Optional[str] = None
    bin_marginal_cum_eps: float = 1e-3

    # v6_minimal — rank-R Lie generators (sum of R outer products per
    # generator). Default `lie_max_rank=1` reduces exactly to the original
    # rank-1 decoder. See doc/v6_minimal_design_2026-05-21.md §M1.
    lie_max_rank: int = 1

    # zero-origin / capacity work (2026-06) — per-rank coefficients. When True
    # (and lie_max_rank > 1) the coefficient head emits one scalar per
    # (generator, rank) so the R rank pieces are independently steerable per
    # cell. Default False = shared-coefficient decoder (byte-identical).
    # See doc/zero_origin_and_capacity_design_2026-06-01.md.
    per_rank_coeff: bool = False

    # v6_min_mixer — ScoreResidualMixer. Default `use_score_residual_mixer=False`
    # keeps the v6_minimal decoder exactly. When True, the mixer is identity-
    # initialised so behaviour at step 0 matches the no-mixer path; deviation
    # requires reconstruction gradient + safety regularizers (gate_balance,
    # gate_ref_safe). See doc/v6_min_mixer_design_2026-05-22.md.
    use_score_residual_mixer: bool = False
    mixer_hidden_dim: int = 64
    mixer_gate_clamp: float = 3.0  # softmax identity gate == gate_clamp/3

    # v23 — pathology-axis interaction terms (explicit products of a learned k-dim
    # pathology projection of z). Default OFF + zero-init => byte-identical. Adds
    # curvature (z-dependent Jacobian) + synergy. interaction_max_order: 2=pairwise,
    # 3=+triples, 4=+quad (co-occurrence count shows 3/4-way are estimable).
    use_pathology_interactions: bool = False
    n_interaction_axes: int = 4
    interaction_max_order: int = 2
    # v26/v27 — pathology-TIED, celltype-conditioned decoder (the *principled* curvature).
    # use_pathology_decoder turns it on; pathology_pairwise False=v26 (single-axis only),
    # True=v27 (+ pairwise p_k·p_l synergy = the LATE-relevant curvature). pathology_rank =
    # low-rank dim of the celltype-specific pattern correction. zero-init => no-op until trained.
    use_pathology_decoder: bool = False
    n_pathology_axes: int = 4
    pathology_pairwise: bool = False
    pathology_rank: int = 8

    # v35 — disease-GATED decoder basis. score += (s ⊙ sigmoid(gate_slope·p)) @ disease_base,
    # p = pathology_head(z). The sigmoid gate makes ∂score/∂z depend on z => G_disease ≠ G_normal
    # (real curvature) concentrated where pathology is high. Per-axis strength s is L1-penalized
    # (pathology_aux.lambda_gate_sparse) so the model SELF-SELECTS which axes curve (dead axes
    # shrink to 0). zero-init disease_base => no-op until trained; default OFF => byte-identical.
    # Requires use_pathology_decoder. gate_rank>0 + gate_celltype_conditioned adds a celltype-
    # specific low-rank disease pattern. See doc/v35_disease_gate_design.md.
    use_disease_gate: bool = False
    gate_slope: float = 1.0
    gate_rank: int = 0
    gate_celltype_conditioned: bool = False
    # learnable Lie rank (2026-07): over-provision lie_max_rank + group-lasso prune. OFF => byte-identical.
    learnable_lie_rank: bool = False
    lambda_lie_rank_sparse: float = 0.0
    lie_rank_warmup_epochs: int = 2


@dataclass
class LossConfig:
    beta_state: float = 1.0
    lambda_tech: float = 1e-2
    lambda_tech_mag: float = 0.0
    lambda_gauge: float = 1e-2
    lambda_bam: float = 1.0
    lambda_cls: float = 0.0
    lambda_ref_center: float = 0.0
    lambda_ref_state: float = 0.0
    lambda_align: float = 0.0
    lambda_white: float = 0.0
    lambda_align_sex: float = 0.0
    align_sex_all_cells: bool = False
    align_use_covariance: bool = True
    align_reference_only: bool = True
    align_whiten_all_cells: bool = True
    align_cov_use_ema: bool = True
    align_ema_decay: float = 0.99
    ref_center_min_cells: int = 2
    tau: float = 1.0
    r_min: float = 0.0
    r_max: float = 5.0
    detach_bam_uncertainty: bool = True
    normalize_bam_weights: bool = True
    normalize_gene_regularizers: bool = True


@dataclass
class OptimConfig:
    lr: float = 1e-3
    weight_decay: float = 1e-4
    betas: tuple[float, float] = (0.9, 0.999)
    eps: float = 1e-8


@dataclass
class SchedulerConfig:
    name: str = "none"  # none | cosine | step
    t_max: int = 50
    eta_min: float = 1e-5
    step_size: int = 10
    gamma: float = 0.5

@dataclass
class TrainConfig:
    epochs: int = 10
    seed: int = 42
    device: str = "cuda"
    amp: bool = True
    amp_dtype: str = "float16"
    grad_clip_norm: Optional[float] = 1.0
    tech_weight_mode: str = "empirical_batch"
    scheduler_step_on: str = "epoch"
    eval_sample_latent: bool = False
    log_every: Optional[int] = 50
    out_dir: str = "outputs/train_run"
    save_every: int = 1
    save_best: bool = True
    save_split_latents: bool = True
    metric_name: str = "loss/total"
    metric_mode: str = "min"

    # early stopping
    early_stopping: bool = False
    early_stopping_patience: int = 5
    early_stopping_min_delta: float = 0.0
    early_stopping_min_epochs: int = 5
    # Fail closed for experiments whose stopping/checkpoint decision is
    # explicitly validation-only.  Legacy runs keep the historical fallback
    # behaviour unless this flag is enabled.
    require_validation_for_early_stopping: bool = False
    # Multi-criteria early stopping. If non-empty, the patience counter only
    # advances when NONE of the listed criteria improves over their own best
    # by their own min_delta. The primary (metric_name) is automatically
    # included as the first criterion.
    # Each entry: {"name", "mode", "min_delta", "weight"(optional), "short_name"(optional)}.
    # See doc/early_stopping_multi_criteria.md.
    early_stopping_criteria: list = field(default_factory=list)
    # How to pick `checkpoint_best.pt`:
    #   "primary"   - original behaviour: improvement of `metric_name` only.
    #   "composite" - weighted z-score sum across all early_stopping_criteria.
    #                 Composite is only used from `early_stopping_min_epochs`
    #                 onward; before that, falls back to primary tracking.
    save_best_mode: str = "primary"
    # Additionally save `checkpoint_best_<short_name>.pt` for each criterion
    # whenever that criterion improves. Independent of save_best_mode.
    save_best_per_criterion: bool = False
    # Numerical safety for composite z-score: std floor.
    composite_std_eps: float = 1e-6
    reload_best_before_test: bool = True
    # Software smoke tests can exercise train+validation without opening the
    # held-out test split. Production runs leave this false.
    skip_final_test: bool = False

    # DDP runtime
    ddp_backend: Optional[str] = None
    ddp_find_unused_parameters: bool = False
    ddp_skip_initial_sync: bool = False
    ddp_timeout_minutes: int = 120
    control_timeout_minutes: int = 720

    # Evaluation mode
    # Default: normal DDP eval. Every rank evaluates its own validation shard,
    # then metrics are reduced by sync_epoch_output().
    # Optional: rank0-only eval for special cases.
    use_rank0_eval: bool = False
    eval_log_every: Optional[int] = None
    rank0_eval_batch_size: Optional[int] = None
    rank0_eval_num_workers: int = 0
    rank0_eval_pin_memory: bool = False
    rank0_eval_persistent_workers: bool = False
    rank0_eval_prefetch_factor: Optional[int] = None

    # Evaluate every N epochs when validation split exists.
    eval_every: int = 1

    # gradient accumulation
    grad_accum_steps: int = 1

    # step-level checkpoint / progress reporting
    save_every_steps: Optional[int] = None
    keep_last_step_checkpoints: int = 3
    progress_log_every: Optional[int] = None

    # diagnostic z-sensitivity check inside Trainer.train_epoch
    debug_z_sensitivity_every: int | None = None
    max_train_steps_per_epoch: Optional[int] = None
    max_eval_steps: Optional[int] = None

    # Epoch-boundary continuation. Unlike the top-level ``warm_start``
    # facility, this restores the complete optimization state and continues
    # epoch numbering from the source checkpoint. Partial, mid-epoch resume is
    # deliberately rejected because sampler/RNG state is not checkpointed.
    resume_checkpoint: Optional[str] = None
    resume_optimizer: bool = True
    resume_history: bool = True
    # Unattended process recovery may continue inside the same immutable run
    # directory.  Disabled by default; enabled only by an explicit launch
    # guard and an epoch-boundary checkpoint.
    allow_in_place_resume: bool = False


@dataclass
class RunConfig:
    data: DataConfig
    loader: LoaderConfig = field(default_factory=LoaderConfig)
    embedding: EmbeddingConfig = field(default_factory=EmbeddingConfig)
    module_tokenizer: ModuleTokenizerConfig = field(default_factory=ModuleTokenizerConfig)
    encoder: EncoderConfig = field(default_factory=EncoderConfig)
    prior: PriorConfig = field(default_factory=PriorConfig)
    decoder: DecoderConfig = field(default_factory=DecoderConfig)
    loss: LossConfig = field(default_factory=LossConfig)
    optim: OptimConfig = field(default_factory=OptimConfig)
    scheduler: SchedulerConfig = field(default_factory=SchedulerConfig)
    train: TrainConfig = field(default_factory=TrainConfig)
    pathology_aux: PathologyAuxConfig = field(default_factory=PathologyAuxConfig)
    metric_contrast: MetricContrastConfig = field(default_factory=MetricContrastConfig)
    nuisance_projection: NuisanceProjectionConfig = field(default_factory=NuisanceProjectionConfig)
    cp_bam: CPBAMConfig = field(default_factory=CPBAMConfig)
    precision_medicine: PrecisionMedicineConfig = field(default_factory=PrecisionMedicineConfig)
    generator_budget_search: GeneratorBudgetSearchConfig = field(
        default_factory=GeneratorBudgetSearchConfig
    )
    learned_generator_count: LearnedGeneratorCountConfig = field(
        default_factory=LearnedGeneratorCountConfig
    )
    learned_pathology_rank: LearnedPathologyRankConfig = field(
        default_factory=LearnedPathologyRankConfig
    )
    prism_module_rescue: PrismModuleRescueTrainingConfig = field(
        default_factory=PrismModuleRescueTrainingConfig
    )
    integrated_phase_curriculum: IntegratedPhaseCurriculumConfig = field(
        default_factory=IntegratedPhaseCurriculumConfig
    )


@dataclass
class RuntimeContext:
    distributed: bool
    rank: int
    local_rank: int
    world_size: int
    is_main_process: bool
    device: torch.device
    backend: Optional[str]
    control_group: Optional[Any]
    control_backend: Optional[str]


# =====================================================================
# Config helpers
# =====================================================================
def _merge_dataclass(dc_cls, user_dict: Optional[Dict[str, Any]]) -> Any:
    base = dc_cls()
    if user_dict is None:
        return base
    for key, value in user_dict.items():
        if not hasattr(base, key):
            raise KeyError(f"Unknown config key '{key}' for {dc_cls.__name__}.")
        setattr(base, key, value)
    return base


def load_config(config_path: str) -> RunConfig:
    with open(config_path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    if "data" not in raw:
        raise KeyError("Config must contain a top-level 'data' section.")

    data = DataConfig(**raw["data"])
    loader = _merge_dataclass(LoaderConfig, raw.get("loader"))
    embedding = _merge_dataclass(EmbeddingConfig, raw.get("embedding"))
    if "module_tokenizer" not in raw:
        raise KeyError("Config must contain a top-level 'module_tokenizer' section.")
    module_tokenizer = ModuleTokenizerConfig(**raw["module_tokenizer"])
    encoder = _merge_dataclass(EncoderConfig, raw.get("encoder"))
    prior = _merge_dataclass(PriorConfig, raw.get("prior"))
    decoder = _merge_dataclass(DecoderConfig, raw.get("decoder"))
    loss = _merge_dataclass(LossConfig, raw.get("loss"))
    optim = _merge_dataclass(OptimConfig, raw.get("optim"))
    scheduler = _merge_dataclass(SchedulerConfig, raw.get("scheduler"))
    train = _merge_dataclass(TrainConfig, raw.get("train"))
    pathology_aux = _merge_dataclass(PathologyAuxConfig, raw.get("pathology_aux"))
    metric_contrast = _merge_dataclass(MetricContrastConfig, raw.get("metric_contrast"))
    nuisance_projection = _merge_dataclass(NuisanceProjectionConfig, raw.get("nuisance_projection"))
    cp_bam = _merge_dataclass(CPBAMConfig, raw.get("cp_bam"))
    precision_medicine = _merge_dataclass(
        PrecisionMedicineConfig, raw.get("precision_medicine")
    )
    generator_budget_search = _merge_dataclass(
        GeneratorBudgetSearchConfig, raw.get("generator_budget_search")
    )
    generator_budget_search.validate()
    learned_generator_count = LearnedGeneratorCountConfig(
        **dict(raw.get("learned_generator_count", {}))
    )
    learned_generator_count.validate()
    learned_pathology_rank = LearnedPathologyRankConfig(
        **dict(raw.get("learned_pathology_rank", {}))
    )
    learned_pathology_rank.validate()
    if bool(learned_pathology_rank.enabled):
        if not bool(decoder.use_pathology_decoder):
            raise ValueError(
                "learned_pathology_rank requires decoder.use_pathology_decoder"
            )
        if int(decoder.pathology_rank) != 0:
            raise ValueError(
                "learned pathology rank uses decoder.pathology_rank=0 as the "
                "automatic algebraic-capacity sentinel; a fixed rank is not allowed"
            )
    if (
        bool(learned_generator_count.enabled)
        and str(learned_generator_count.mode) == "joint"
        and bool(generator_budget_search.enabled)
    ):
        raise ValueError(
            "joint learned_generator_count cannot be combined with the "
            "external generator_budget_search controller"
        )
    if bool(learned_generator_count.enabled) and str(
        learned_generator_count.mode
    ) == "joint":
        if bool(cp_bam.enabled) or bool(metric_contrast.enabled):
            raise ValueError(
                "joint learned_generator_count currently requires cp_bam and "
                "metric_contrast to remain disabled because their auxiliary "
                "decoder re-forwards do not yet consume the sampled gate"
            )
    prism_module_rescue = _merge_dataclass(
        PrismModuleRescueTrainingConfig,
        raw.get("prism_module_rescue"),
    )
    prism_module_rescue.validate()
    integrated_phase_curriculum = _merge_dataclass(
        IntegratedPhaseCurriculumConfig,
        raw.get("integrated_phase_curriculum"),
    )
    integrated_phase_curriculum.validate()

    if bool(integrated_phase_curriculum.enabled):
        pathology_rank_contract = (
            int(decoder.pathology_rank) == 0
            if bool(learned_pathology_rank.enabled)
            else int(decoder.pathology_rank) == 8
        )
        required_phase1 = {
            "decoder.use_pathology_decoder": bool(decoder.use_pathology_decoder),
            "decoder.pathology_rank_contract": pathology_rank_contract,
            "decoder.pathology_pairwise=false": not bool(decoder.pathology_pairwise),
            "decoder.use_pathology_interactions=false": not bool(decoder.use_pathology_interactions),
            "pathology_aux.enabled": bool(pathology_aux.enabled),
            "pathology_aux.lambda_path_aux=0.3": abs(float(pathology_aux.lambda_path_aux) - 0.3) <= 1.0e-12,
            "pathology_aux.lambda_axis_specificity=0.05": abs(float(pathology_aux.lambda_axis_specificity) - 0.05) <= 1.0e-12,
            "pathology_aux.lambda_pathology_tie=0.1": abs(float(pathology_aux.lambda_pathology_tie) - 0.1) <= 1.0e-12,
            "celltype_projection_disabled": not bool(nuisance_projection.project_celltype),
            "personal_rank=2": int(precision_medicine.personal_rank) == 2,
        }
        if bool(
            getattr(
                learned_pathology_rank,
                "require_canonical_phase1_rank8",
                False,
            )
        ):
            required_phase1.update(
                {
                    "learned_rank.warmup_fixed_rank=8": int(
                        learned_pathology_rank.warmup_fixed_rank
                    )
                    == 8,
                    "learned_rank.warmup_covers_phase1": int(
                        learned_pathology_rank.warmup_end_epoch
                    )
                    >= int(integrated_phase_curriculum.phase1_end_epoch),
                    "learned_rank.search_after_warmup": int(
                        learned_pathology_rank.soft_start_epoch
                    )
                    > int(learned_pathology_rank.warmup_end_epoch),
                }
            )
        if bool(
            getattr(
                learned_pathology_rank,
                "require_active_pathology_route_during_search",
                False,
            )
        ):
            required_phase1[
                "learned_rank.phase2_pathology_route_active"
            ] = float(integrated_phase_curriculum.phase2_pathology_scale) > 0.0
        failed = [name for name, passed in required_phase1.items() if not passed]
        if failed:
            raise ValueError(
                "integrated Phase-I parity contract failed: " + ", ".join(failed)
            )

    if bool(
        getattr(
            learned_generator_count,
            "require_frozen_pathology_rank_before_search",
            False,
        )
    ):
        if not bool(learned_pathology_rank.enabled):
            raise ValueError(
                "generator search requires an enabled learned pathology rank"
            )
        if int(learned_generator_count.start_epoch) < int(
            learned_pathology_rank.freeze_epoch
        ):
            raise ValueError(
                "generator search cannot start before pathology-rank freeze"
            )

    return RunConfig(
        data=data,
        loader=loader,
        embedding=embedding,
        module_tokenizer=module_tokenizer,
        encoder=encoder,
        prior=prior,
        decoder=decoder,
        loss=loss,
        optim=optim,
        scheduler=scheduler,
        train=train,
        pathology_aux=pathology_aux,
        metric_contrast=metric_contrast,
        nuisance_projection=nuisance_projection,
        cp_bam=cp_bam,
        precision_medicine=precision_medicine,
        generator_budget_search=generator_budget_search,
        learned_generator_count=learned_generator_count,
        learned_pathology_rank=learned_pathology_rank,
        prism_module_rescue=prism_module_rescue,
        integrated_phase_curriculum=integrated_phase_curriculum,
    )


# =====================================================================
# Distributed helpers
# =====================================================================
def setup_runtime(cfg: RunConfig) -> RuntimeContext:
    def _runtime_log(message: str) -> None:
        if os.environ.get("KMLEE_TRAIN_BOOT_DEBUG", "0") == "1":
            rank_env = os.environ.get("RANK", "?")
            print(f"[runtime-debug rank={rank_env}] {message}", flush=True)

    world_size = int(os.environ.get("WORLD_SIZE", "1"))
    rank = int(os.environ.get("RANK", "0"))
    local_rank = int(os.environ.get("LOCAL_RANK", "0"))
    distributed = world_size > 1
    _runtime_log(
        f"start world_size={world_size} rank={rank} local_rank={local_rank} distributed={distributed}"
    )

    use_cuda = cfg.train.device.startswith("cuda") and torch.cuda.is_available()
    backend = cfg.train.ddp_backend
    if backend is None and distributed:
        backend = "nccl" if use_cuda else "gloo"
    _runtime_log(f"device_select use_cuda={use_cuda} backend={backend}")

    control_group = None
    control_backend = None

    if distributed:
        if use_cuda:
            torch.cuda.set_device(local_rank)
            device = torch.device(f"cuda:{local_rank}")
        else:
            device = torch.device("cpu")

        dist.init_process_group(
            backend=backend,
            timeout=timedelta(minutes=cfg.train.ddp_timeout_minutes),
        )
        _runtime_log("init_process_group:done")

        if backend == "nccl":
            try:
                _runtime_log("control_group:create:start")
                control_group = dist.new_group(
                    backend="gloo",
                    timeout=timedelta(minutes=cfg.train.control_timeout_minutes),
                )
                control_backend = "gloo"
                _runtime_log("control_group:create:done")
            except Exception:
                control_group = None
                control_backend = None
                _runtime_log("control_group:create:failed")
        else:
            control_group = dist.group.WORLD
            control_backend = backend
    else:
        if use_cuda:
            device = torch.device(cfg.train.device)
        else:
            device = torch.device("cpu")
        backend = None

    return RuntimeContext(
        distributed=distributed,
        rank=rank,
        local_rank=local_rank,
        world_size=world_size,
        is_main_process=(rank == 0),
        device=device,
        backend=backend,
        control_group=control_group,
        control_backend=control_backend,
    )


def cleanup_runtime(rt: RuntimeContext) -> None:
    if not rt.distributed or not dist.is_initialized():
        return

    try:
        control_barrier(rt)
    except Exception:
        pass

    try:
        if rt.control_group is not None and rt.control_group is not dist.group.WORLD:
            dist.destroy_process_group(rt.control_group)
    except Exception:
        pass

    try:
        dist.destroy_process_group()
    except Exception:
        pass


def barrier(rt: RuntimeContext) -> None:
    if rt.distributed and dist.is_initialized():
        if rt.control_group is not None:
            dist.barrier(group=rt.control_group)
            return
        if rt.backend == "nccl" and rt.device.type == "cuda":
            dist.barrier(device_ids=[rt.local_rank])
        else:
            dist.barrier()


def control_barrier(rt: RuntimeContext) -> None:
    if not rt.distributed or not dist.is_initialized():
        return
    if rt.control_group is None:
        return
    dist.barrier(group=rt.control_group)


def broadcast_object_list_control(rt: RuntimeContext, obj_list: list[Any], *, src: int = 0) -> None:
    if not rt.distributed or not dist.is_initialized():
        return
    group = rt.control_group if rt.control_group is not None else None
    dist.broadcast_object_list(obj_list, src=src, group=group)


def unwrap_system(system: Any) -> OrdinalBAMSystem:
    if hasattr(system, "module") and isinstance(system.module, nn.Module):
        return system.module
    return system


def run_system_integrity_checks(
    system: Any,
    *,
    stage: str,
    distributed: bool = False,
) -> Optional[Dict[str, Any]]:
    """Run optional model-specific integrity checks without affecting defaults."""

    target = unwrap_system(system)
    hook = getattr(target, "_fixed_generator_integrity_hook", None)
    if hook is None:
        return None
    audit = hook(stage=stage, distributed=distributed)
    if audit is not None and not isinstance(audit, dict):
        raise RuntimeError("model integrity hook must return a dict or None")
    return audit


@contextmanager
def ddp_initial_sync_context(skip: bool):
    """Optionally bypass DDP's constructor-time parameter verification/sync.

    PyTorch 2.5 on the SV7 runtime used for this project can hang inside the
    DDP constructor's initial NCCL verification/sync. Newer PyTorch versions
    expose an official `init_sync=False` switch, but this environment does not.
    When enabled by config, this context mirrors that narrow behaviour and
    leaves normal gradient all-reduce untouched.
    """
    if not skip:
        yield
        return

    import torch.nn.parallel.distributed as ddp_distributed

    old_verify = ddp_distributed._verify_param_shape_across_processes
    old_sync = ddp_distributed._sync_module_states

    def _skip_verify_param_shape_across_processes(*args, **kwargs):
        return None

    def _skip_sync_module_states(*args, **kwargs):
        return None

    ddp_distributed._verify_param_shape_across_processes = (
        _skip_verify_param_shape_across_processes
    )
    ddp_distributed._sync_module_states = _skip_sync_module_states
    try:
        yield
    finally:
        ddp_distributed._verify_param_shape_across_processes = old_verify
        ddp_distributed._sync_module_states = old_sync


class DDPSystemProxy:
    """
    Preserve the attribute interface expected by Trainer while executing the
    forward pass through DistributedDataParallel.
    """

    def __init__(self, ddp_model: DDP) -> None:
        self._ddp_model = ddp_model

    def __call__(self, *args, **kwargs):
        return self._ddp_model(*args, **kwargs)

    def train(self, mode: bool = True):
        self._ddp_model.train(mode)
        return self

    def eval(self):
        self._ddp_model.eval()
        return self

    def parameters(self, recurse: bool = True):
        return self._ddp_model.parameters(recurse=recurse)

    def named_parameters(self, prefix: str = "", recurse: bool = True):
        return self._ddp_model.named_parameters(prefix=prefix, recurse=recurse)

    def state_dict(self, *args, **kwargs):
        return self._ddp_model.module.state_dict(*args, **kwargs)

    def load_state_dict(self, *args, **kwargs):
        return self._ddp_model.module.load_state_dict(*args, **kwargs)

    def to(self, *args, **kwargs):
        self._ddp_model.module.to(*args, **kwargs)
        return self

    @property
    def module(self):
        return self._ddp_model.module

    def __getattr__(self, name: str):
        return getattr(self._ddp_model.module, name)


def sync_epoch_output(epoch_out: EpochOutput, rt: RuntimeContext) -> EpochOutput:
    """
    Reduce epoch metrics across ranks using the control-plane group when available.

    Object collectives over the NCCL default group can be fragile in long runs;
    a small Gloo control group is safer for Python dictionaries and scalar metrics.
    """
    if not rt.distributed:
        return epoch_out

    payload = {
        "metrics": epoch_out.metrics,
        "n_examples": epoch_out.n_examples,
    }
    gathered = [None for _ in range(rt.world_size)]
    group = rt.control_group if rt.control_group is not None else None
    dist.all_gather_object(gathered, payload, group=group)

    total_examples = int(sum(int(item["n_examples"]) for item in gathered))
    if total_examples == 0:
        return EpochOutput(metrics={}, n_examples=0)

    all_keys = set()
    for item in gathered:
        all_keys.update(item["metrics"].keys())

    reduced_metrics: Dict[str, float] = {}
    for key in sorted(all_keys):
        weighted_sum = 0.0
        for item in gathered:
            n = int(item["n_examples"])
            weighted_sum += float(item["metrics"].get(key, 0.0)) * n
        reduced_metrics[key] = weighted_sum / total_examples

    # The PRISM masked NLL is a global weighted numerator/denominator ratio.
    # Averaging rank-local ratios would reintroduce cell-count bias, so compute
    # the ratio again only after both sufficient statistics have been reduced.
    finalize_prism_epoch_metrics(reduced_metrics)

    return EpochOutput(metrics=reduced_metrics, n_examples=total_examples)


# =====================================================================
# Reproducibility
# =====================================================================
def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


# =====================================================================
# Build objects
# =====================================================================
def build_datasets(cfg: RunConfig) -> tuple[OrdinalScDataset, Optional[OrdinalScDataset], Optional[OrdinalScDataset]]:
    common = dict(
        zarr_path=cfg.data.zarr_path,
        spec_path=cfg.data.spec_path,
        matrix_path=cfg.data.matrix_path,
        counts_layer=cfg.data.counts_layer,
        prefer_raw=cfg.data.prefer_raw,
        split_key=cfg.data.split_key,
        celltype_key=cfg.data.celltype_key,
        batch_key=cfg.data.batch_key,
        return_log_input=cfg.data.return_log_input,
        reference_scaler_path=cfg.data.reference_scaler_path,
        return_x_gene_scalar=cfg.data.return_x_gene_scalar,
        x_gene_scalar_clip=cfg.data.x_gene_scalar_clip,
        x_gene_scalar_eps=cfg.data.x_gene_scalar_eps,
        max_csr_data_span=cfg.data.max_csr_data_span,
        max_span_to_selected_nnz_ratio=cfg.data.max_span_to_selected_nnz_ratio,
    )

    train_ds = OrdinalScDataset(split="train", **common)

    val_ds = None
    test_ds = None
    try:
        val_ds = OrdinalScDataset(split="val", **common)
    except Exception as exc:
        if bool(cfg.train.require_validation_for_early_stopping):
            raise RuntimeError(
                "validation-only training control requires a readable val split"
            ) from exc
        val_ds = None

    if not bool(cfg.train.skip_final_test):
        try:
            test_ds = OrdinalScDataset(split="test", **common)
        except Exception:
            test_ds = None

    return train_ds, val_ds, test_ds


def build_loaders(
    cfg: RunConfig,
    train_ds: OrdinalScDataset,
    val_ds: Optional[OrdinalScDataset],
    test_ds: Optional[OrdinalScDataset],
    *,
    distributed: bool,
) -> tuple[
    DataLoader,
    Optional[DataLoader],
    Optional[DataLoader],
    Optional[DistributedSampler],
    Optional[DistributedSampler],
    Optional[DistributedSampler],
]:
    """
    Build train/val/test loaders.

    Default behaviour:
      - train: DistributedSampler under DDP
      - val/test: DistributedSampler under DDP, no shuffle, no drop_last

    This restores all-rank DDP evaluation, which was the stable behaviour in
    earlier runs. Rank0-only evaluation is handled separately and only used
    when cfg.train.use_rank0_eval=True.
    """
    eval_batch_size = cfg.loader.eval_batch_size or cfg.loader.batch_size

    train_sampler = None
    val_sampler = None
    test_sampler = None

    if distributed:
        world_size = int(os.environ.get("WORLD_SIZE", "1"))
        rank = int(os.environ.get("RANK", "0"))
        train_sampler = DistributedSampler(
            train_ds,
            num_replicas=world_size,
            rank=rank,
            shuffle=cfg.loader.shuffle_train,
            drop_last=cfg.data.drop_last_train,
        )

        if val_ds is not None:
            val_sampler = DistributedSampler(
                val_ds,
                num_replicas=world_size,
                rank=rank,
                shuffle=False,
                drop_last=False,
            )

        if test_ds is not None:
            test_sampler = DistributedSampler(
                test_ds,
                num_replicas=world_size,
                rank=rank,
                shuffle=False,
                drop_last=False,
            )

    train_loader = DataLoader(
        train_ds,
        batch_size=cfg.loader.batch_size,
        shuffle=(cfg.loader.shuffle_train and train_sampler is None),
        sampler=train_sampler,
        num_workers=cfg.data.num_workers,
        pin_memory=cfg.data.pin_memory,
        persistent_workers=(
            cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
        ),
        drop_last=cfg.data.drop_last_train,
    )

    val_loader = None
    if val_ds is not None:
        val_loader = DataLoader(
            val_ds,
            batch_size=eval_batch_size,
            shuffle=False,
            sampler=val_sampler,
            num_workers=cfg.data.num_workers,
            pin_memory=cfg.data.pin_memory,
            persistent_workers=(
                cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
            ),
            drop_last=False,
        )

    test_loader = None
    if test_ds is not None:
        test_loader = DataLoader(
            test_ds,
            batch_size=eval_batch_size,
            shuffle=False,
            sampler=test_sampler,
            num_workers=cfg.data.num_workers,
            pin_memory=cfg.data.pin_memory,
            persistent_workers=(
                cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
            ),
            drop_last=False,
        )

    return train_loader, val_loader, test_loader, train_sampler, val_sampler, test_sampler


def build_plain_eval_loader(
    cfg: RunConfig,
    ds: OrdinalScDataset,
    *,
    batch_size: Optional[int] = None,
    shuffle: bool = False,
) -> DataLoader:
    """
    Plain single-process eval loader.

    Use cases:
      1. Optional rank0-only eval when cfg.train.use_rank0_eval=True.
      2. Final test / latent export on rank0 after training.

    This is NOT used for normal DDP validation when use_rank0_eval=False.
    """
    actual_batch_size = (
        batch_size
        or cfg.train.rank0_eval_batch_size
        or cfg.loader.eval_batch_size
        or cfg.loader.batch_size
    )

    num_workers = max(0, int(cfg.train.rank0_eval_num_workers))
    pin_memory = bool(cfg.train.rank0_eval_pin_memory)

    kwargs: Dict[str, Any] = dict(
        dataset=ds,
        batch_size=int(actual_batch_size),
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=pin_memory,
        persistent_workers=(
            num_workers > 0 and bool(cfg.train.rank0_eval_persistent_workers)
        ),
        drop_last=False,
    )

    if num_workers > 0 and cfg.train.rank0_eval_prefetch_factor is not None:
        kwargs["prefetch_factor"] = int(cfg.train.rank0_eval_prefetch_factor)

    return DataLoader(**kwargs)


@contextmanager
def use_unwrapped_system_for_rank0_only_ops(trainer: Trainer):
    """
    Temporarily bypass DDP for rank-0-only evaluation / latent export.

    Running forward passes through DistributedDataParallel on only one rank can
    still trigger collectives (e.g. buffer broadcasts). For single-rank eval we
    therefore swap in the raw underlying module.
    """
    original_system = trainer.system
    trainer.system = unwrap_system(original_system)
    try:
        yield
    finally:
        trainer.system = original_system


@torch.no_grad()
def evaluate_epoch_rank0_only(
    trainer: Trainer,
    loader: Optional[DataLoader],
    runtime: RuntimeContext,
    *,
    log_every: Optional[int] = None,
    max_steps: Optional[int] = None,
) -> Optional[EpochOutput]:
    """
    Distributed-safe rank-0-only evaluation.

    Critical invariant:
        In DDP, EVERY rank must call this function in the same order.
        Only rank 0 actually iterates over the validation/test loader.
        Non-rank0 ranks wait in the control-plane broadcast.

    This prevents the failure mode where rank 0 is still validating while ranks
    1..N enter the next epoch's DDP all-reduce, causing NCCL timeout.
    """
    if not runtime.distributed:
        if loader is None:
            return None
        return trainer.evaluate_epoch(loader, log_every=log_every, max_steps=max_steps)

    result = None
    if runtime.is_main_process:
        if loader is None:
            raise RuntimeError(
                "Rank 0 must receive a non-None loader for rank-0-only evaluation."
            )
        with use_unwrapped_system_for_rank0_only_ops(trainer):
            result = trainer.evaluate_epoch(
                loader,
                log_every=log_every,
                max_steps=max_steps,
            )

    shared = [
        None
        if result is None
        else {
            "metrics": result.metrics,
            "n_examples": result.n_examples,
            "elapsed_sec": result.elapsed_sec,
        }
    ]

    # All ranks must participate here. Prefer the Gloo control group for object
    # broadcast; avoid NCCL object collectives unless no control group exists.
    broadcast_object_list_control(runtime, shared, src=0)
    payload = shared[0]

    if payload is None:
        return None

    return EpochOutput(
        metrics=dict(payload["metrics"]),
        n_examples=int(payload["n_examples"]),
        elapsed_sec=float(payload["elapsed_sec"]),
    )


# =====================================================================
# Module registry helpers
# =====================================================================
def _load_registry_json(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    required = ["gene_names", "module_names", "membership_binary"]
    missing = [k for k in required if k not in data]
    if missing:
        raise KeyError(f"Registry JSON is missing required key(s): {missing}")
    return data


def _registry_membership_tensor(
    registry: Dict[str, Any],
    expected_gene_names: np.ndarray,
) -> torch.Tensor:
    reg_gene_names = np.asarray(registry["gene_names"], dtype=object)
    expected_gene_names = np.asarray(expected_gene_names, dtype=object)
    if len(reg_gene_names) != len(expected_gene_names) or not np.array_equal(reg_gene_names, expected_gene_names):
        raise ValueError(
            "Gene order mismatch between registry and ordinal spec. "
            "The final registry gene union must be exactly the dataset/spec gene order."
        )

    membership = torch.tensor(registry["membership_binary"], dtype=torch.float32)
    if membership.ndim != 2:
        raise ValueError(f"membership_binary must have shape [M, G], got {tuple(membership.shape)}.")
    if membership.shape[1] != len(expected_gene_names):
        raise ValueError(
            f"membership gene dimension mismatch: {membership.shape[1]} vs {len(expected_gene_names)}."
        )
    if not torch.isfinite(membership).all():
        raise ValueError("membership_binary contains non-finite values.")
    if (membership < 0).any():
        raise ValueError("membership_binary must be non-negative.")
    if (membership.sum(dim=1) <= 0).any():
        bad = torch.where(membership.sum(dim=1) <= 0)[0].tolist()[:10]
        raise ValueError(f"Registry contains empty module rows. Examples: {bad}")
    return membership


def _load_optional_activity_weight(
    cfg: RunConfig,
    registry: Dict[str, Any],
    membership: torch.Tensor,
) -> Optional[torch.Tensor]:
    path = cfg.module_tokenizer.activity_weight_path
    key = cfg.module_tokenizer.activity_weight_key

    if path is not None:
        p = Path(path)
        if p.suffix.lower() == ".npz":
            obj = np.load(p, allow_pickle=True)
            if key not in obj:
                raise KeyError(f"Activity weight npz does not contain key '{key}'.")
            arr = obj[key]
        else:
            arr = np.load(p, allow_pickle=True)
        activity = torch.tensor(arr, dtype=torch.float32)
    elif key in registry:
        activity = torch.tensor(registry[key], dtype=torch.float32)
    elif "activity_weight" in registry:
        activity = torch.tensor(registry["activity_weight"], dtype=torch.float32)
    else:
        return None

    if activity.shape != membership.shape:
        raise ValueError(
            f"activity_weight shape mismatch: expected {tuple(membership.shape)}, got {tuple(activity.shape)}."
        )
    return activity


def _select_decoder_module_indices(
    membership: torch.Tensor,
    *,
    requested_n_generators: Optional[int],
    strategy: str,
    require_generator_per_module: bool = False,
) -> torch.Tensor:
    """Select registry rows used by decoder generators.

    ``largest_protect_singletons`` starts from the historical ``largest``
    selection, then forces every size-one module into the fixed-size budget.
    When rows must be displaced, it removes selected non-singletons one at a
    time.  Each greedy step minimises the number of genes that would lose their
    final selected-module cover; ties prefer the smaller module and then the
    lower registry index.  Retained historical rows keep their historical
    order, while newly forced singleton rows are appended in registry order.

    Returning indices separately keeps the generator-to-registry mapping
    available for checkpoint manifests without changing the long-standing
    mask-returning helper below.
    """

    M, _ = membership.shape
    strategy = str(strategy).lower()
    if require_generator_per_module:
        if strategy != "all":
            raise ValueError(
                "require_generator_per_module=True requires "
                "module_tokenizer.decoder_mask_strategy='all'."
            )
        if (
            requested_n_generators is not None
            and int(requested_n_generators) != M
        ):
            raise ValueError(
                "require_generator_per_module=True requires "
                f"decoder.n_generators to equal registry modules ({M}), "
                f"got {requested_n_generators}."
            )
    if requested_n_generators is None:
        n = M
    else:
        n = int(requested_n_generators)
        if n <= 0:
            raise ValueError("decoder.n_generators must be positive when provided.")
        n = min(n, M)

    if strategy == "all":
        if requested_n_generators is not None and n < M:
            # With an explicit cap, "all" means first n after preserving order.
            idx = torch.arange(n)
        else:
            idx = torch.arange(M)
    elif strategy == "first":
        idx = torch.arange(n)
    elif strategy == "largest":
        sizes = membership.sum(dim=1)
        idx = torch.argsort(sizes, descending=True)[:n]
    elif strategy == "largest_protect_singletons":
        support = membership > 0
        support_sizes = support.sum(dim=1)
        singleton_idx = torch.where(support_sizes == 1)[0]
        if singleton_idx.numel() > n:
            raise ValueError(
                "largest_protect_singletons cannot protect every singleton "
                f"module: found {singleton_idx.numel()} singletons but only "
                f"{n} decoder generator slots."
            )

        # Start from the exact historical selection/order.
        sizes = membership.sum(dim=1)
        legacy_idx = torch.argsort(sizes, descending=True)[:n]
        legacy_set = set(int(i) for i in legacy_idx.detach().cpu().tolist())
        missing_singletons = sorted(
            int(i)
            for i in singleton_idx.detach().cpu().tolist()
            if int(i) not in legacy_set
        )

        if not missing_singletons:
            idx = legacy_idx
        else:
            selected = legacy_idx.detach().cpu().tolist() + missing_singletons
            coverage = support[selected].sum(dim=0)
            candidates = [
                int(i)
                for i in legacy_idx.detach().cpu().tolist()
                if int(support_sizes[int(i)]) > 1
            ]
            removed: set[int] = set()

            for _ in range(len(missing_singletons)):
                if not candidates:
                    raise ValueError(
                        "largest_protect_singletons has no selected "
                        "non-singleton module left to displace."
                    )
                uniquely_covered = coverage == 1
                removal_key = []
                for registry_idx in candidates:
                    coverage_loss = int(
                        (support[registry_idx] & uniquely_covered).sum().item()
                    )
                    removal_key.append(
                        (
                            coverage_loss,
                            int(support_sizes[registry_idx].item()),
                            registry_idx,
                        )
                    )
                _, _, drop_idx = min(removal_key)
                removed.add(drop_idx)
                coverage = coverage - support[drop_idx].to(coverage.dtype)
                candidates.remove(drop_idx)

            retained = [
                int(i)
                for i in legacy_idx.detach().cpu().tolist()
                if int(i) not in removed
            ]
            final_indices = retained + missing_singletons
            if len(final_indices) != n:
                raise RuntimeError(
                    "largest_protect_singletons internal size mismatch: "
                    f"expected {n}, got {len(final_indices)}."
                )
            idx = torch.tensor(
                final_indices,
                dtype=torch.long,
                device=membership.device,
            )
    elif strategy == "largest_append_singletons":
        support_sizes = (membership > 0).sum(dim=1)
        singleton_indices = torch.where(support_sizes == 1)[0]
        n_singletons = int(singleton_indices.numel())
        if requested_n_generators is not None and int(
            requested_n_generators
        ) > M:
            raise ValueError(
                "largest_append_singletons requested more decoder generators "
                f"({requested_n_generators}) than registry modules ({M})."
            )
        if n_singletons > n:
            raise ValueError(
                "largest_append_singletons cannot append every singleton "
                f"module: found {n_singletons} singletons but only {n} "
                "decoder generator slots."
            )

        non_singleton_budget = n - n_singletons
        non_singleton_indices = torch.where(support_sizes > 1)[0]
        if non_singleton_budget > int(non_singleton_indices.numel()):
            raise ValueError(
                "largest_append_singletons does not have enough "
                "non-singleton modules to fill the requested budget: "
                f"need {non_singleton_budget}, found "
                f"{non_singleton_indices.numel()}."
            )

        # Deterministic order: larger membership first, then lower registry
        # index. On the canonical registry, the first 256 rows are exactly the
        # historical largest-256 set; all singleton rows follow in registry
        # order.
        membership_sizes = membership.sum(dim=1)
        ordered_non_singletons = sorted(
            (int(i) for i in non_singleton_indices.detach().cpu().tolist()),
            key=lambda registry_idx: (
                -float(membership_sizes[registry_idx].item()),
                registry_idx,
            ),
        )
        final_indices = (
            ordered_non_singletons[:non_singleton_budget]
            + sorted(
                int(i)
                for i in singleton_indices.detach().cpu().tolist()
            )
        )
        if len(final_indices) != n:
            raise RuntimeError(
                "largest_append_singletons internal size mismatch: "
                f"expected {n}, got {len(final_indices)}."
            )
        idx = torch.tensor(
            final_indices,
            dtype=torch.long,
            device=membership.device,
        )
    else:
        raise ValueError(
            "module_tokenizer.decoder_mask_strategy must be 'all', 'first', "
            "'largest', 'largest_protect_singletons', or "
            "'largest_append_singletons'."
        )

    return idx


def _select_decoder_module_masks(
    membership: torch.Tensor,
    *,
    requested_n_generators: Optional[int],
    strategy: str,
    require_generator_per_module: bool = False,
) -> torch.Tensor:
    idx = _select_decoder_module_indices(
        membership,
        requested_n_generators=requested_n_generators,
        strategy=strategy,
        require_generator_per_module=require_generator_per_module,
    )

    masks = membership[idx].clone().float()
    if (masks.sum(dim=1) <= 0).any():
        raise ValueError("Selected decoder module masks contain empty rows.")
    return masks


def _sha256_int64_sequence(values: tuple[int, ...]) -> str:
    """Hash an ordered integer sequence with a platform-independent encoding."""

    payload = np.asarray(values, dtype="<i8").tobytes(order="C")
    return hashlib.sha256(payload).hexdigest()


def _sha256_string_sequence(values: tuple[str, ...]) -> str:
    """Hash an ordered string sequence without delimiter ambiguity."""

    payload = json.dumps(
        list(values),
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _attach_decoder_architecture_metadata(
    decoder: LieActionOrdinalDecoder,
    *,
    registry_path: str | Path,
    selection_strategy: str,
    membership: torch.Tensor,
    module_names: list[Any] | tuple[Any, ...] | np.ndarray,
    gene_names: list[Any] | tuple[Any, ...] | np.ndarray,
    selected_indices: Optional[torch.Tensor],
) -> None:
    """Attach plain, non-state-dict provenance to a constructed decoder."""

    registry_module_names = tuple(str(v) for v in module_names)
    registry_gene_names = tuple(str(v) for v in gene_names)
    registry_module_count, registry_gene_count = map(int, membership.shape)
    if len(registry_module_names) != registry_module_count:
        raise ValueError(
            "Registry module_names length does not match membership rows: "
            f"{len(registry_module_names)} vs {registry_module_count}."
        )
    if len(registry_gene_names) != registry_gene_count:
        raise ValueError(
            "Registry gene_names length does not match membership columns: "
            f"{len(registry_gene_names)} vs {registry_gene_count}."
        )

    if selected_indices is None:
        selected_registry_indices: tuple[int, ...] = ()
        registry_mapping_available = False
    else:
        selected_registry_indices = tuple(
            int(v) for v in selected_indices.detach().cpu().tolist()
        )
        registry_mapping_available = True
        if len(selected_registry_indices) != int(decoder.n_generators):
            raise ValueError(
                "Selected registry index count does not match decoder generators: "
                f"{len(selected_registry_indices)} vs {decoder.n_generators}."
            )
        if len(set(selected_registry_indices)) != len(selected_registry_indices):
            raise ValueError("Selected decoder registry indices must be unique.")
        if any(
            index < 0 or index >= registry_module_count
            for index in selected_registry_indices
        ):
            raise ValueError("Selected decoder registry index is out of bounds.")

    selected_registry_module_names = tuple(
        registry_module_names[index] for index in selected_registry_indices
    )
    registry_support = membership.detach().cpu() > 0
    registry_module_sizes = tuple(
        int(v) for v in registry_support.sum(dim=1).tolist()
    )
    registry_singleton_indices = tuple(
        int(v)
        for v in torch.where(registry_support.sum(dim=1) == 1)[0].tolist()
    )

    metadata: Dict[str, Any] = {
        "registry_path": str(registry_path),
        "selection_strategy": str(selection_strategy).lower(),
        "registry_mapping_available": bool(registry_mapping_available),
        "registry_module_count": registry_module_count,
        "registry_gene_count": registry_gene_count,
        "registry_module_sizes": registry_module_sizes,
        "registry_singleton_indices": registry_singleton_indices,
        "selected_registry_indices": selected_registry_indices,
        "selected_registry_module_names": selected_registry_module_names,
        "selection_sha256": _sha256_int64_sequence(
            selected_registry_indices
        ),
        "registry_module_order_sha256": _sha256_string_sequence(
            registry_module_names
        ),
        "registry_gene_order_sha256": _sha256_string_sequence(
            registry_gene_names
        ),
    }

    # These are deliberately ordinary Python attributes, not parameters or
    # buffers. They therefore make the decoder auditable without altering old
    # checkpoint/state_dict contracts.
    decoder.selected_registry_module_indices = selected_registry_indices
    decoder.selected_registry_module_names = selected_registry_module_names
    decoder.decoder_architecture_metadata = metadata


def _copy_decoder_architecture_metadata_to_system(system: nn.Module) -> None:
    """Expose decoder provenance on the final base or V3 system wrapper."""

    target = unwrap_system(system)
    decoder = getattr(target, "decoder", None)
    metadata = getattr(decoder, "decoder_architecture_metadata", None)
    if metadata is None:
        raise RuntimeError(
            "Decoder architecture metadata is missing after build_system."
        )
    target.selected_registry_module_indices = (
        decoder.selected_registry_module_indices
    )
    target.selected_registry_module_names = (
        decoder.selected_registry_module_names
    )
    target.decoder_architecture_metadata = metadata


def _build_resolved_decoder_architecture_manifest(
    cfg: RunConfig,
    system: nn.Module,
) -> Dict[str, Any]:
    """Resolve the exact decoder architecture actually constructed in memory."""

    target = unwrap_system(system)
    decoder = getattr(target, "decoder", None)
    metadata = getattr(decoder, "decoder_architecture_metadata", None)
    if decoder is None or metadata is None:
        raise RuntimeError(
            "Cannot resolve decoder architecture before metadata is attached."
        )

    support = (decoder.module_masks.detach().cpu() > 0)
    module_sizes_tensor = support.sum(dim=1)
    degree_tensor = support.sum(dim=0)
    rank_vector_tensor = (
        decoder.generator_rank_mask.detach().cpu().sum(dim=1)
    )

    module_sizes = [int(v) for v in module_sizes_tensor.tolist()]
    gene_membership_degree = [int(v) for v in degree_tensor.tolist()]
    generator_rank_vector = [int(v) for v in rank_vector_tensor.tolist()]
    singleton_generator_positions = [
        int(v) for v in torch.where(module_sizes_tensor == 1)[0].tolist()
    ]
    singleton_generator_vector = [
        bool(size == 1) for size in module_sizes
    ]
    selected_registry_indices = list(metadata["selected_registry_indices"])
    selected_names = list(metadata["selected_registry_module_names"])
    selected_singleton_registry_indices = [
        selected_registry_indices[position]
        for position in singleton_generator_positions
        if position < len(selected_registry_indices)
    ]

    covered_gene_count = int((degree_tensor > 0).sum().item())
    registry_gene_count = int(metadata["registry_gene_count"])
    if int(degree_tensor.numel()) != registry_gene_count:
        raise ValueError(
            "Decoder mask gene count no longer matches the registry contract: "
            f"{degree_tensor.numel()} vs {registry_gene_count}."
        )
    positive_degree = degree_tensor[degree_tensor > 0]
    positive_degree_mean = (
        float(positive_degree.float().mean().item())
        if positive_degree.numel()
        else 0.0
    )

    return {
        "schema_version": "kmlee_bam.resolved_decoder_architecture.v1",
        "decoder_mask_strategy": metadata["selection_strategy"],
        "use_module_masks": bool(cfg.decoder.use_module_masks),
        "requested_decoder_generator_count": (
            None
            if cfg.decoder.n_generators is None
            else int(cfg.decoder.n_generators)
        ),
        "require_generator_per_module": bool(
            getattr(cfg.decoder, "require_generator_per_module", False)
        ),
        "registry_mapping_available": metadata[
            "registry_mapping_available"
        ],
        "registry_path": metadata["registry_path"],
        "registry_module_count": int(metadata["registry_module_count"]),
        "module_tokenizer_module_count": int(
            getattr(target.module_tokenizer, "n_modules", 0)
        ),
        "decoder_generator_count": int(decoder.n_generators),
        "hard_generator_mask_supported": bool(
            hasattr(decoder, "generator_active_mask")
        ),
        "initial_active_generator_count": int(
            (
                decoder.generator_active_mask.detach().cpu() > 0
            ).sum().item()
        )
        if hasattr(decoder, "generator_active_mask")
        else int(decoder.n_generators),
        "registry_gene_count": registry_gene_count,
        "decoder_gene_count": int(decoder.n_genes),
        "selected_registry_indices": selected_registry_indices,
        "selected_module_names": selected_names,
        "registry_module_sizes": list(metadata["registry_module_sizes"]),
        "module_sizes": module_sizes,
        "selected_module_sizes": module_sizes,
        "registry_singleton_indices": list(
            metadata["registry_singleton_indices"]
        ),
        "selected_singleton_registry_indices": (
            selected_singleton_registry_indices
        ),
        "singleton_generator_positions": singleton_generator_positions,
        "singleton_generator_vector": singleton_generator_vector,
        "singleton_generator_count": len(singleton_generator_positions),
        "lie_max_rank": int(decoder.lie_max_rank),
        "singleton_lie_rank": (
            None
            if decoder.singleton_lie_rank is None
            else int(decoder.singleton_lie_rank)
        ),
        "generator_rank_vector": generator_rank_vector,
        "active_rank_component_count": int(sum(generator_rank_vector)),
        "module_generator_normalization": str(
            getattr(decoder, "module_generator_normalization", "none")
        ),
        "selected_membership_nnz": int(degree_tensor.sum().item()),
        "covered_gene_count": covered_gene_count,
        "uncovered_gene_count": registry_gene_count - covered_gene_count,
        "gene_coverage_fraction": (
            covered_gene_count / registry_gene_count
            if registry_gene_count
            else 0.0
        ),
        "gene_membership_degree": gene_membership_degree,
        "degree_summary": {
            "min": int(degree_tensor.min().item())
            if degree_tensor.numel()
            else 0,
            "mean": float(degree_tensor.float().mean().item())
            if degree_tensor.numel()
            else 0.0,
            "median": float(degree_tensor.float().median().item())
            if degree_tensor.numel()
            else 0.0,
            "max": int(degree_tensor.max().item())
            if degree_tensor.numel()
            else 0,
            "positive_mean": positive_degree_mean,
        },
        "selection_sha256": metadata["selection_sha256"],
        "selection_hash_encoding": "ordered little-endian int64",
        "registry_module_order_sha256": metadata[
            "registry_module_order_sha256"
        ],
        "registry_gene_order_sha256": metadata[
            "registry_gene_order_sha256"
        ],
        "trainable_parameter_count": int(
            sum(
                parameter.numel()
                for parameter in target.parameters()
                if parameter.requires_grad
            )
        ),
        "total_parameter_count": int(
            sum(parameter.numel() for parameter in target.parameters())
        ),
    }


def build_system(cfg: RunConfig, train_ds: OrdinalScDataset) -> OrdinalBAMSystem:
    n_genes = len(train_ds.spec.gene_names)
    n_bins = int(train_ds.spec.n_bins)
    n_celltypes = len(train_ds.spec.celltype_vocab)
    n_tech = len(train_ds.spec.batch_vocab)

    registry = _load_registry_json(cfg.module_tokenizer.registry_json_path)
    membership = _registry_membership_tensor(registry, train_ds.spec.gene_names)
    activity_weight = _load_optional_activity_weight(cfg, registry, membership)

    gene_embedding = GeneExpressionEmbedding(
        n_genes=n_genes,
        d_model=cfg.embedding.d_model,
        n_bins=n_bins,
        use_continuous=cfg.embedding.use_continuous,
        use_cls=cfg.embedding.use_cls,
        dropout=cfg.embedding.dropout,
        cont_hidden_dim=cfg.embedding.cont_hidden_dim,
        layer_norm_eps=cfg.embedding.layer_norm_eps,
        init_std=cfg.embedding.init_std,
    )

    module_tokenizer = GeneModuleTokenizer(
        membership_weight=membership,
        d_model=cfg.embedding.d_model,
        activity_weight=activity_weight,
        activity_weight_normalization=cfg.module_tokenizer.activity_weight_normalization,
        activity_default=cfg.module_tokenizer.activity_default,
        pooling=cfg.module_tokenizer.pooling,
        preserve_cls=cfg.module_tokenizer.preserve_cls,
        use_module_id_embedding=cfg.module_tokenizer.use_module_id_embedding,
        module_id_init_scale=cfg.module_tokenizer.module_id_init_scale,
        module_id_scale_min=cfg.module_tokenizer.module_id_scale_min,
        module_id_scale_max=cfg.module_tokenizer.module_id_scale_max,
        use_activity_projection=cfg.module_tokenizer.use_activity_projection,
        dropout=cfg.module_tokenizer.dropout,
        layer_norm_eps=cfg.module_tokenizer.layer_norm_eps,
        init_std=cfg.module_tokenizer.init_std,
    )

    state_encoder = StateEncoder(
        d_model=cfg.embedding.d_model,
        n_heads=cfg.encoder.n_heads,
        n_layers=cfg.encoder.n_layers,
        d_z=cfg.encoder.d_z,
        condition_on_celltype=cfg.encoder.condition_on_celltype,
        n_celltypes=n_celltypes,
        celltype_embed_dim=cfg.encoder.celltype_embed_dim,
        posterior_input_norm=cfg.encoder.posterior_input_norm,
        pooling=cfg.encoder.pooling,
        use_cls_token=cfg.encoder.use_cls_token,
        mean_exclude_cls=cfg.encoder.mean_exclude_cls,
        agp_heads=cfg.encoder.agp_heads,
        agp_exclude_cls=cfg.encoder.agp_exclude_cls,
        agp_variant=cfg.encoder.agp_variant,
        agp_temperature_init=cfg.encoder.agp_temperature_init,
        agp_temperature_min=cfg.encoder.agp_temperature_min,
        agp_temperature_max=cfg.encoder.agp_temperature_max,
        agp_mean_residual_init=cfg.encoder.agp_mean_residual_init,
        agp_mean_residual_max=cfg.encoder.agp_mean_residual_max,
        agp_diversity_weight=cfg.encoder.agp_diversity_weight,
        agp_diversity_margin=cfg.encoder.agp_diversity_margin,
        agp_query_orthogonality_weight=cfg.encoder.agp_query_orthogonality_weight,
        agp_entropy_band_weight=cfg.encoder.agp_entropy_band_weight,
        agp_min_effective_tokens=cfg.encoder.agp_min_effective_tokens,
        agp_max_effective_tokens=cfg.encoder.agp_max_effective_tokens,
        persistent_agp_attention_heads=cfg.encoder.persistent_agp_attention_heads,
        persistent_agp_ffn_dim=cfg.encoder.persistent_agp_ffn_dim,
        persistent_agp_update_gate_init=cfg.encoder.persistent_agp_update_gate_init,
        persistent_agp_update_gate_max=cfg.encoder.persistent_agp_update_gate_max,
        compute_cell_uncertainty=cfg.encoder.compute_cell_uncertainty,
        uncertainty_exclude_cls=cfg.encoder.uncertainty_exclude_cls,
        tech_risk_head=cfg.encoder.tech_risk_head,
        posterior_hidden_dim=cfg.encoder.posterior_hidden_dim,
        posterior_dropout=cfg.encoder.posterior_dropout,
        logvar_min=cfg.encoder.logvar_min,
        logvar_max=cfg.encoder.logvar_max,
        d_ff=cfg.encoder.d_ff,
        stochastic_attention=cfg.encoder.stochastic_attention,
        distribution=cfg.encoder.distribution,
        sigma_mode=cfg.encoder.sigma_mode,
        sigma=cfg.encoder.sigma,
        sigma_min=cfg.encoder.sigma_min,
        sigma_max=cfg.encoder.sigma_max,
        weibull_k=cfg.encoder.weibull_k,
        prior_d_mid=cfg.encoder.prior_d_mid,
        attn_dropout=cfg.encoder.attn_dropout,
        proj_dropout=cfg.encoder.proj_dropout,
        ffn_dropout=cfg.encoder.ffn_dropout,
        qkv_bias=cfg.encoder.qkv_bias,
        layer_norm_eps=cfg.encoder.layer_norm_eps,
        init_std=cfg.encoder.init_std,
        final_norm=cfg.encoder.final_norm,
    )

    prior = CellTypePrior(
        n_celltypes=n_celltypes,
        d_z=cfg.encoder.d_z,
        sigma_init=cfg.prior.sigma_init,
        sigma_min=cfg.prior.sigma_min,
        residual_eps=cfg.prior.residual_eps,
        residual_clamp=cfg.prior.residual_clamp,
    )

    module_masks = None
    selected_module_indices = None
    n_generators = cfg.decoder.n_generators
    if cfg.decoder.use_module_masks:
        selected_module_indices = _select_decoder_module_indices(
            membership,
            requested_n_generators=n_generators,
            strategy=cfg.module_tokenizer.decoder_mask_strategy,
            require_generator_per_module=bool(
                getattr(
                    cfg.decoder,
                    "require_generator_per_module",
                    False,
                )
            ),
        )
        module_masks = membership[selected_module_indices].clone().float()
        if (module_masks.sum(dim=1) <= 0).any():
            raise ValueError("Selected decoder module masks contain empty rows.")
        n_generators = int(module_masks.shape[0])
    elif n_generators is None:
        n_generators = 32

    bin_marginal_freq_tensor = None
    init_mode = str(cfg.decoder.threshold_init_mode).lower()
    if init_mode == "data_quantile":
        from kmlee_bam.data.bin_marginal import maybe_load_or_compute_bin_marginal
        bin_marginal_freq_np = maybe_load_or_compute_bin_marginal(
            train_ds,
            cache_path=cfg.decoder.bin_marginal_cache_path,
            n_sample_cells=int(cfg.decoder.bin_marginal_n_sample),
            seed=int(cfg.decoder.bin_marginal_seed),
            smoothing=float(cfg.decoder.bin_marginal_smoothing),
        )
        bin_marginal_freq_tensor = torch.from_numpy(bin_marginal_freq_np)
    elif init_mode != "uniform":
        raise ValueError(
            f"decoder.threshold_init_mode must be 'uniform' or 'data_quantile', "
            f"got {cfg.decoder.threshold_init_mode!r}."
        )

    resolved_pathology_rank = int(getattr(cfg.decoder, "pathology_rank", 8))
    learned_pathology_rank_cfg = getattr(cfg, "learned_pathology_rank", None)
    if learned_pathology_rank_cfg is not None and bool(
        learned_pathology_rank_cfg.enabled
    ):
        resolved_pathology_rank = algebraic_pathology_rank_capacity(
            n_celltypes=n_celltypes,
            n_pathology_axes=int(getattr(cfg.decoder, "n_pathology_axes", 4)),
            n_genes=n_genes,
        )

    decoder = LieActionOrdinalDecoder(
        n_genes=n_genes,
        n_celltypes=n_celltypes,
        n_tech=n_tech,
        d_z=cfg.encoder.d_z,
        n_bins=n_bins,
        n_sex=getattr(cfg.decoder, "n_sex", None),
        n_generators=n_generators,
        coeff_hidden_dim=cfg.decoder.coeff_hidden_dim,
        coeff_dropout=cfg.decoder.coeff_dropout,
        gamma_scale=cfg.decoder.gamma_scale,
        generator_init_scale=cfg.decoder.generator_init_scale,
        translation_init_scale=cfg.decoder.translation_init_scale,
        use_affine_translation=cfg.decoder.use_affine_translation,
        use_direct_state_head=cfg.decoder.use_direct_state_head,
        direct_state_init_scale=cfg.decoder.direct_state_init_scale,
        module_masks=module_masks,
        singleton_lie_rank=getattr(
            cfg.decoder, "singleton_lie_rank", None
        ),
        module_generator_normalization=str(
            getattr(
                cfg.decoder,
                "module_generator_normalization",
                "none",
            )
        ),
        center_thresholds=cfg.decoder.center_thresholds,
        threshold_init_spacing=cfg.decoder.threshold_init_spacing,
        prob_eps=cfg.decoder.prob_eps,
        init_std=cfg.decoder.init_std,
        bin_marginal_freq=bin_marginal_freq_tensor,
        bin_marginal_cum_eps=float(cfg.decoder.bin_marginal_cum_eps),
        # v6_minimal — see doc/v6_minimal_design_2026-05-21.md §M1.
        lie_max_rank=int(getattr(cfg.decoder, "lie_max_rank", 1)),
        # zero-origin / capacity — see doc/zero_origin_and_capacity_design_2026-06-01.md.
        per_rank_coeff=bool(getattr(cfg.decoder, "per_rank_coeff", False)),
        # v6_min_mixer — see doc/v6_min_mixer_design_2026-05-22.md.
        use_score_residual_mixer=bool(getattr(cfg.decoder, "use_score_residual_mixer", False)),
        mixer_hidden_dim=int(getattr(cfg.decoder, "mixer_hidden_dim", 64)),
        mixer_gate_clamp=float(getattr(cfg.decoder, "mixer_gate_clamp", 3.0)),
        # v23 pathology-axis interactions (default OFF => byte-identical)
        use_pathology_interactions=bool(getattr(cfg.decoder, "use_pathology_interactions", False)),
        n_interaction_axes=int(getattr(cfg.decoder, "n_interaction_axes", 4)),
        interaction_max_order=int(getattr(cfg.decoder, "interaction_max_order", 2)),
        # v26/v27 pathology-tied, celltype-conditioned decoder (default OFF => byte-identical)
        use_pathology_decoder=bool(getattr(cfg.decoder, "use_pathology_decoder", False)),
        n_pathology_axes=int(getattr(cfg.decoder, "n_pathology_axes", 4)),
        pathology_pairwise=bool(getattr(cfg.decoder, "pathology_pairwise", False)),
        pathology_rank=resolved_pathology_rank,
        pathology_rank_rng_compat_prefix=(
            int(learned_pathology_rank_cfg.warmup_fixed_rank)
            if learned_pathology_rank_cfg is not None
            and bool(learned_pathology_rank_cfg.enabled)
            and bool(learned_pathology_rank_cfg.require_canonical_phase1_rank8)
            else 0
        ),
        # v35 disease-gated basis (default OFF => byte-identical). See doc/v35_disease_gate_design.md.
        use_disease_gate=bool(getattr(cfg.decoder, "use_disease_gate", False)),
        gate_slope=float(getattr(cfg.decoder, "gate_slope", 1.0)),
        gate_rank=int(getattr(cfg.decoder, "gate_rank", 0)),
        gate_celltype_conditioned=bool(getattr(cfg.decoder, "gate_celltype_conditioned", False)),
        # learnable Lie rank (2026-07) — see doc/lie_rank_learnable_design_v3.md. OFF => byte-identical.
        learnable_lie_rank=bool(getattr(cfg.decoder, "learnable_lie_rank", False)),
        lambda_lie_rank_sparse=float(getattr(cfg.decoder, "lambda_lie_rank_sparse", 0.0)),
        lie_rank_warmup_epochs=int(getattr(cfg.decoder, "lie_rank_warmup_epochs", 2)),
    )
    if learned_pathology_rank_cfg is not None and bool(
        learned_pathology_rank_cfg.enabled
    ):
        decoder.pathology_rank_gate = LearnedPathologyRankGate(
            resolved_pathology_rank,
            config=learned_pathology_rank_cfg,
        )
        decoder.pathology_rank_capacity = int(resolved_pathology_rank)
        if int(os.environ.get("RANK", "0")) == 0:
            print(
                "[learned-pathology-rank] "
                f"capacity=auto:{resolved_pathology_rank} "
                f"celltypes={n_celltypes} "
                f"axes={int(getattr(cfg.decoder, 'n_pathology_axes', 4))} "
                "fixed_rank=false target_rank=none minimum_rank=0 "
                f"warmup_fixed_rank={learned_pathology_rank_cfg.warmup_fixed_rank} "
                f"warmup_end_epoch={learned_pathology_rank_cfg.warmup_end_epoch} "
                f"soft_epoch={learned_pathology_rank_cfg.soft_start_epoch} "
                f"hard_epoch={learned_pathology_rank_cfg.hard_start_epoch} "
                f"freeze_epoch={learned_pathology_rank_cfg.freeze_epoch}",
                flush=True,
            )
    _attach_decoder_architecture_metadata(
        decoder,
        registry_path=cfg.module_tokenizer.registry_json_path,
        selection_strategy=(
            cfg.module_tokenizer.decoder_mask_strategy
            if cfg.decoder.use_module_masks
            else "unmasked"
        ),
        membership=membership,
        module_names=registry["module_names"],
        gene_names=registry["gene_names"],
        selected_indices=selected_module_indices,
    )

    if int(os.environ.get("RANK", "0")) == 0:
        degree = decoder.module_overlap_degree.detach().cpu()
        positive_degree = degree[degree > 0]
        degree_min = (
            float(positive_degree.min()) if len(positive_degree) else 0.0
        )
        degree_mean = (
            float(positive_degree.float().mean())
            if len(positive_degree)
            else 0.0
        )
        degree_max = (
            float(positive_degree.max()) if len(positive_degree) else 0.0
        )
        print(
            "[decoder-module-contract] "
            f"registry_modules={int(membership.shape[0])} "
            f"generators={decoder.n_generators} "
            f"strategy={cfg.module_tokenizer.decoder_mask_strategy} "
            f"require_1to1={bool(getattr(cfg.decoder, 'require_generator_per_module', False))} "
            f"covered_genes={decoder.n_module_covered_genes}/{n_genes} "
            f"singleton_generators={decoder.n_singleton_generators} "
            f"singleton_rank={decoder.singleton_lie_rank} "
            f"regular_rank={decoder.lie_max_rank} "
            f"normalization={decoder.module_generator_normalization} "
            f"overlap_degree[min/mean/max]="
            f"{degree_min:.0f}/{degree_mean:.2f}/{degree_max:.0f}",
            flush=True,
        )
     

    system = OrdinalBAMSystem(
        gene_embedding=gene_embedding,
        module_tokenizer=module_tokenizer,
        state_encoder=state_encoder,
        prior=prior,
        decoder=decoder,
        classifier_head=None,
    )
    _copy_decoder_architecture_metadata_to_system(system)

    # Joint learned-count mode is part of the ordinary training graph.  Unlike
    # the legacy gate-only runner, no model parameters are frozen and no donor
    # is removed from the training set.  One exact-forward hard gate controls
    # both the Lie action and its affine translation through the decoder's
    # existing ``generator_gate`` contract.
    learned_count_cfg = getattr(cfg, "learned_generator_count", None)
    if (
        learned_count_cfg is not None
        and bool(learned_count_cfg.enabled)
        and str(learned_count_cfg.mode) == "joint"
    ):
        if int(system.decoder.n_generators) != int(
            system.module_tokenizer.n_modules
        ):
            raise ValueError(
                "joint generator count requires one decoder generator per "
                "registered module"
            )
        support = system.decoder.module_masks.detach() > 0
        protected = torch.zeros(
            int(system.decoder.n_generators),
            dtype=torch.bool,
            device=support.device,
        )
        if bool(learned_count_cfg.protect_singletons):
            protected |= support.sum(dim=1) == 1
        if bool(learned_count_cfg.protect_unique_gene_coverage):
            uniquely_covered_genes = support.sum(dim=0) == 1
            if bool(uniquely_covered_genes.any()):
                protected |= support[:, uniquely_covered_genes].any(dim=1)
        system.generator_count_gate = HardBinaryConcreteGeneratorGate(
            int(system.decoder.n_generators),
            config=learned_count_cfg,
            protected_mask=protected,
        )
        system.generator_count_objective = ConstrainedGeneratorCountObjective(
            ("full_nll", "isolated_nll"),
            config=learned_count_cfg,
        )
        system.generator_count_config = learned_count_cfg
        system.register_buffer(
            "generator_count_epoch_state",
            torch.tensor(1, dtype=torch.int64),
            persistent=True,
        )
        system.generator_count_joint_training = True

    # v20 multi-pathology auxiliary head (donor-balanced aggregate supervision on
    # z_perp). Registered as a system submodule so its params join the optimizer and
    # the checkpoint; the EMA bank buffers are non-persistent (per-rank). The loss is
    # computed in the trainer (it needs the per-batch donor/celltype grouping). Built
    # only when enabled so default runs are byte-identical.
    pa_cfg = getattr(cfg, "pathology_aux", None)
    if pa_cfg is not None and pa_cfg.enabled:
        system.pathology_aux_head = PathologyAuxHead(
            d_z=cfg.encoder.d_z,
            n_donor=int(getattr(train_ds, "n_donor", 0)),
            n_celltype=n_celltypes,
            config=pa_cfg,
        )

    # CP-BAM consensus-private decomposition head (z_perp -> z_common/z_private). Mirrors
    # pathology_aux: registered submodule so its params join the optimizer/DDP/checkpoint; the
    # loss is computed in the trainer (needs per-batch donor x celltype grouping). d_z/n_donor/
    # n_celltype are filled from the run context. Built only when enabled => default runs identical.
    cpb_cfg = getattr(cfg, "cp_bam", None)
    if cpb_cfg is not None and cpb_cfg.enabled:
        cpb_cfg.d_z = cfg.encoder.d_z
        cpb_cfg.n_donor = int(getattr(train_ds, "n_donor", 0))
        cpb_cfg.n_celltype = n_celltypes
        system.cp_bam_head = CPBAMHead(cpb_cfg)

    # v35 metric-contrast head — parameter-free; bends the disease-gate metric toward amyloid via
    # grads through the decoder score. Attached to base system; V3 wrapper carries it via getattr.
    # Built only when enabled ⇒ default runs never call it ⇒ byte-identical.
    mc_cfg = getattr(cfg, "metric_contrast", None)
    if mc_cfg is not None and mc_cfg.enabled:
        system.metric_contrast_head = MetricContrastHead(config=mc_cfg)

    # v25 nuisance projector — attached to the base system; the V3 wrapper carries it via
    # getattr (see system.py). Its EMA buffers are persistent. DDP is configured
    # with broadcast_buffers=False; rank consistency instead comes from the
    # projector's explicit all-reduce of sufficient statistics before each EMA
    # update. Built only when enabled ⇒ default runs are byte-identical.
    np_cfg = getattr(cfg, "nuisance_projection", None)
    if np_cfg is not None and np_cfg.enabled:
        system.nuisance_projector = LatentNuisanceProjector(
            d_z=cfg.encoder.d_z,
            n_celltype=n_celltypes,
            config=np_cfg,
        )

    # PRISM E2E precision-medicine score decomposition.  Context regression is
    # fitted on train64 donors only; held-out donor contexts are transformed by
    # those fixed coefficients.  The target cell type is masked later inside
    # PrecisionMedicineHead, so validation/test target expression is never an
    # input to its own personal code.
    pm_cfg = getattr(cfg, "precision_medicine", None)
    if pm_cfg is not None and pm_cfg.enabled:
        if not pm_cfg.context_npz_path:
            raise ValueError("precision_medicine.context_npz_path is required")
        table = DonorContextTable.from_npz(
            pm_cfg.context_npz_path,
            min_cells=int(pm_cfg.min_context_cells),
        )
        dataset_donors = np.asarray(getattr(train_ds, "donor_vocab", []), dtype=object).astype(str)
        dataset_celltypes = np.asarray(train_ds.spec.celltype_vocab, dtype=object).astype(str)
        dataset_regions = np.asarray(getattr(train_ds, "region_vocab", []), dtype=object).astype(str)
        if not np.array_equal(table.donor_names.astype(str), dataset_donors):
            raise ValueError("PRISM context donor order does not match the cell dataset")
        if not np.array_equal(table.celltype_names.astype(str), dataset_celltypes):
            raise ValueError("PRISM context celltype order does not match the ordinal spec")
        if not np.array_equal(table.region_names.astype(str), dataset_regions):
            raise ValueError("PRISM context region order does not match the cell dataset")
        if int(table.module.shape[-1]) != int(module_tokenizer.n_modules):
            raise ValueError("PRISM context module count does not match the tokenizer")
        registry_module_names = np.asarray(registry["module_names"], dtype=object).astype(str)
        if not np.array_equal(table.module_names.astype(str), registry_module_names):
            raise ValueError("PRISM context module names/order do not match the registry")
        if int(table.latent.shape[-1]) != int(cfg.encoder.d_z):
            raise ValueError("PRISM context latent dimension does not match encoder.d_z")
        dataset_donor_age = np.asarray(
            getattr(train_ds, "donor_age_years", []), dtype=np.float32
        )
        if dataset_donor_age.shape != table.age.shape:
            raise ValueError("PRISM context age vector does not match the cell dataset")
        # The original epoch-12 extractor retained only coarse category lower
        # bounds.  Use the audited exact/top-coded age definition from the raw
        # zarr when residualising the frozen support table, without touching its
        # expression summaries.
        table.age = dataset_donor_age.copy()
        table_train_donors = np.flatnonzero(
            table.donor_split.astype(str) == "train"
        )
        fit_donors = np.unique(
            np.asarray(train_ds.donor_ids, dtype=np.int64)[
                np.asarray(train_ds.row_idx, dtype=np.int64)
            ]
        )
        if not set(int(x) for x in fit_donors.tolist()).issubset(
            set(int(x) for x in table_train_donors.tolist())
        ):
            raise ValueError(
                "PRISM fit donors contain a validation/test donor"
            )
        legacy_gate_only_search = bool(
            cfg.learned_generator_count.enabled
            and str(cfg.learned_generator_count.mode) == "gate_only"
        )
        expected_fit_count = (
            int(cfg.learned_generator_count.architecture_donor_count)
            if legacy_gate_only_search
            else 0
        )
        if legacy_gate_only_search:
            expected_weight_count = int(table_train_donors.size) - expected_fit_count
            if int(fit_donors.size) != expected_weight_count:
                raise ValueError(
                    "PRISM weight-train donor count mismatch: "
                    f"expected {expected_weight_count}, found "
                    f"{fit_donors.size}"
                )
        elif isinstance(
            getattr(train_ds, "generator_exact_training_allowlist", None),
            dict,
        ):
            exact_allowlist = train_ds.generator_exact_training_allowlist
            expected_count = int(
                exact_allowlist.get("observed_optimizer_donor_count", -1)
            )
            if (
                exact_allowlist.get("row_level_donor_allowlist_enforced")
                is not True
                or expected_count not in {44, 54}
                or int(fit_donors.size) != expected_count
            ):
                raise ValueError(
                    "generator exact context residualizer lacks a sealed "
                    "44/54-donor row allow-list"
                )
            observed_names = tuple(
                str(value) for value in table.donor_names[fit_donors].tolist()
            )
            declared_names = tuple(
                str(value)
                for value in exact_allowlist.get("declared_donor_names", ())
            )
            if set(observed_names) != set(declared_names) or len(
                declared_names
            ) != expected_count:
                raise ValueError(
                    "generator exact context fit donors differ from the "
                    "sealed job allow-list"
                )
        elif bool(
            getattr(train_ds, "canonical_w44_source_enabled", False)
        ):
            allowlist_evidence = getattr(
                train_ds, "sealed_donor_allowlist", None
            )
            if not isinstance(allowlist_evidence, dict):
                raise ValueError(
                    "canonical W44 context fit lacks sealed allow-list evidence"
                )
            if (
                allowlist_evidence.get("partition") != "W44"
                or allowlist_evidence.get("row_level_allowlist_enforced")
                is not True
                or int(fit_donors.size) != 44
            ):
                raise ValueError(
                    "canonical W44 context residualizer must fit exactly the "
                    "sealed 44-donor dataset view"
                )
            observed_names = tuple(
                str(value) for value in table.donor_names[fit_donors].tolist()
            )
            if observed_names != tuple(
                str(value)
                for value in allowlist_evidence.get("donor_names", ())
            ):
                raise ValueError(
                    "canonical W44 context fit donor names differ from the "
                    "sealed allow-list"
                )
        elif fit_donors.size != table_train_donors.size:
            raise ValueError(
                "PRISM non-architecture run must fit every train donor: "
                f"{fit_donors.size} vs {table_train_donors.size}"
            )
        prepared = table.prepare(fit_donors, ridge=float(pm_cfg.context_ridge))
        # Runtime-observed fit membership is resolved before constructing the
        # head because the optional module-local artifact is sealed against this
        # exact ordered donor allow-list.  Validation/test donors can therefore
        # never contribute to either reliability or the input clip.
        fit_donor_names = tuple(
            str(value) for value in table.donor_names[fit_donors].tolist()
        )
        module_local_reliability = None
        module_local_input_clip = None
        module_local_output_cap = None
        module_local_artifact = None
        module_local_cap_artifact = None
        module_local_graph_artifact = None
        module_local_compartment_adjacency = None
        if bool(pm_cfg.module_local_enabled):
            reliability_path = pm_cfg.module_local_reliability_npz_path
            if not reliability_path:
                raise ValueError(
                    "module_local_enabled requires "
                    "module_local_reliability_npz_path"
                )
            output_cap_path = pm_cfg.module_local_output_cap_npz_path
            output_cap_checkpoint_sha256 = (
                pm_cfg.module_local_output_cap_source_checkpoint_sha256
            )
            output_cap_config_sha256 = (
                pm_cfg.module_local_output_cap_source_config_sha256
            )
            if (
                not output_cap_path
                or not output_cap_checkpoint_sha256
                or not output_cap_config_sha256
            ):
                raise ValueError(
                    "module_local_enabled requires a sealed output-cap artifact "
                    "and frozen-comparator checkpoint/config SHA-256"
                )
            registry_path = Path(
                cfg.module_tokenizer.registry_json_path
            ).expanduser().resolve()
            context_path = Path(pm_cfg.context_npz_path).expanduser().resolve()
            activity_path_raw = cfg.module_tokenizer.activity_weight_path
            activity_path = (
                Path(activity_path_raw).expanduser().resolve()
                if activity_path_raw
                else registry_path
            )
            module_local_artifact = load_module_local_reliability_artifact(
                reliability_path,
                expected_context_names=tuple(
                    str(value) for value in table.context_names.tolist()
                ),
                expected_module_names=tuple(
                    str(value) for value in table.module_names.tolist()
                ),
                expected_train_donor_names=fit_donor_names,
                expected_source_context_sha256=module_local_sha256_file(
                    context_path
                ),
                expected_registry_sha256=module_local_sha256_file(
                    registry_path
                ),
                expected_activity_dictionary_sha256=(
                    module_local_sha256_file(activity_path)
                ),
            )
            context_sha256 = module_local_sha256_file(context_path)
            module_local_cap_artifact = (
                load_module_local_output_cap_artifact(
                    output_cap_path,
                    expected_train_donor_names=fit_donor_names,
                    expected_checkpoint_sha256=(
                        output_cap_checkpoint_sha256
                    ),
                    expected_source_config_sha256=output_cap_config_sha256,
                    expected_source_context_sha256=context_sha256,
                    expected_absolute_quantile=float(
                        pm_cfg.module_local_output_cap_quantile
                    ),
                    expected_personal_rank=2,
                )
            )
            module_local_output_cap = float(
                module_local_cap_artifact.output_cap
            )
            module_local_reliability = torch.from_numpy(
                module_local_artifact.effective_reliability
            )
            if bool(pm_cfg.module_local_nonlinear_enabled):
                graph_path = pm_cfg.module_local_compartment_graph_npz_path
                if not graph_path:
                    raise ValueError(
                        "module_local_nonlinear_enabled requires "
                        "module_local_compartment_graph_npz_path"
                    )
                module_local_graph_artifact = (
                    load_module_local_compartment_graph_artifact(
                        graph_path,
                        expected_module_names=tuple(
                            str(value) for value in table.module_names.tolist()
                        ),
                        expected_registry_sha256=module_local_sha256_file(
                            registry_path
                        ),
                    )
                )
                module_local_compartment_adjacency = torch.from_numpy(
                    module_local_graph_artifact.adjacency
                )
            fit_observed = table.observed[fit_donors].astype(bool)
            fit_values = np.abs(
                prepared.module_residual[fit_donors][fit_observed]
            ).astype(np.float64, copy=False)
            if fit_values.size == 0 or not np.isfinite(fit_values).all():
                raise ValueError(
                    "module-local train-only residuals are empty or non-finite"
                )
            module_local_input_clip = float(
                np.quantile(
                    fit_values,
                    float(pm_cfg.module_local_input_clip_quantile),
                )
            )
            if not math.isfinite(module_local_input_clip) or module_local_input_clip <= 0.0:
                raise ValueError(
                    "module-local train-only input clip must be finite and positive"
                )
        interaction_basis = None
        if bool(pm_cfg.interaction_enabled):
            interaction_artifact = fit_pairwise_pathology_residualizer(
                torch.from_numpy(table.pathology[fit_donors]),
                torch.from_numpy(~table.pathology_missing[fit_donors]),
                ridge=float(pm_cfg.interaction_ridge),
                minimum_complete_donors=int(
                    pm_cfg.interaction_min_complete_donors
                ),
                scale_floor=float(pm_cfg.interaction_scale_floor),
            )
            interaction_basis = PairwisePathologyInteractionBasis(
                interaction_artifact
            )
            if not bool(interaction_artifact.eligible.any()):
                raise ValueError(
                    "PRISM pathology interactions are enabled but no pair has "
                    "sufficient complete training-donor support"
                )
        # Dormant Phase-II construction must not advance the Phase-I RNG
        # stream.  This preserves the canonical initialization of later
        # Phase-I components (notably the sex adversary) under seed 42.
        _phase2_rng_state = torch.random.get_rng_state()
        torch.manual_seed(int(cfg.train.seed) + 2_608_241)
        system.precision_head = PrecisionMedicineHead(
            config=pm_cfg,
            source_module=torch.from_numpy(prepared.module_residual),
            source_latent=torch.from_numpy(prepared.latent_residual),
            source_observed=torch.from_numpy(table.observed),
            source_reliability=torch.from_numpy(
                np.minimum(
                    table.cell_count.astype(np.float32)
                    / float(max(int(pm_cfg.context_reliability_cap), 1)),
                    1.0,
                )
                * table.observed.astype(np.float32)
            ),
            context_celltype=torch.from_numpy(table.context_celltype),
            context_region=torch.from_numpy(table.context_region),
            n_celltypes=n_celltypes,
            n_regions=len(table.region_names),
            interaction_basis=interaction_basis,
            donor_pathology=torch.from_numpy(table.pathology),
            donor_pathology_valid=torch.from_numpy(~table.pathology_missing),
            donor_age_z=torch.from_numpy(
                np.nan_to_num(
                    (
                        table.age
                        - float(np.nanmean(table.age[fit_donors]))
                    )
                    / max(float(np.nanstd(table.age[fit_donors])), 1.0e-6)
                ).astype(np.float32)
            ),
            donor_sex=torch.from_numpy(table.sex.astype(np.float32)),
            donor_adnc=torch.from_numpy(
                np.nan_to_num(
                    (
                        table.adnc_report
                        - float(np.nanmean(table.adnc_report[fit_donors]))
                    )
                    / max(
                        float(np.nanstd(table.adnc_report[fit_donors])),
                        1.0e-6,
                    )
                ).astype(np.float32)
            ),
            fit_donor_mask=torch.from_numpy(
                np.isin(
                    np.arange(len(table.donor_names), dtype=np.int64),
                    fit_donors,
                )
            ),
            module_local_reliability=module_local_reliability,
            module_local_input_clip=module_local_input_clip,
            module_local_output_cap=module_local_output_cap,
            module_local_compartment_adjacency=(
                module_local_compartment_adjacency
            ),
        )
        torch.random.set_rng_state(_phase2_rng_state)
        # Runtime-observed fit membership, recorded from the actual dataset
        # view that constructed the context residualizer.  Final fixed-mask
        # retraining uses this evidence to prove that freshly built train64
        # buffers are retained when W44 model weights are warm-started.
        exact_allowlist = getattr(
            train_ds, "generator_exact_training_allowlist", None
        )
        context_fit_partition = (
            str(exact_allowlist["fit_partition"])
            if isinstance(exact_allowlist, dict)
            else
            "train64"
            if int(fit_donors.size) == 64
            and int(table_train_donors.size) == 64
            and np.array_equal(
                np.sort(fit_donors), np.sort(table_train_donors)
            )
            else "weight_subset"
        )
        precision_context_fit_provenance = {
            "fit_partition": context_fit_partition,
            "observed_fit_donor_ids": [
                int(value) for value in fit_donors.tolist()
            ],
            "observed_fit_donor_names": list(fit_donor_names),
            "observed_fit_donor_count": int(fit_donors.size),
            "table_train_donor_count": int(table_train_donors.size),
            "context_npz_path": str(Path(pm_cfg.context_npz_path).resolve()),
            "row_level_dataset_view_observed": True,
        }
        system.precision_context_fit_provenance = (
            precision_context_fit_provenance
        )
        system.precision_head.context_fit_provenance = (
            precision_context_fit_provenance
        )
        if module_local_artifact is not None:
            precision_module_local_provenance = {
                "artifact_path": module_local_artifact.source_path,
                "artifact_sha256": module_local_artifact.artifact_sha256,
                "schema_version": "kmlee_bam.module_local_reliability.v1",
                "split_seed": int(module_local_artifact.split_seed),
                "split_rule": module_local_artifact.split_rule,
                "source_context_sha256": (
                    module_local_artifact.source_context_sha256
                ),
                "registry_sha256": module_local_artifact.registry_sha256,
                "activity_dictionary_sha256": (
                    module_local_artifact.activity_dictionary_sha256
                ),
                "train_donor_allowlist_sha256": (
                    module_local_artifact.train_donor_allowlist_sha256
                ),
                "fit_donor_count": len(fit_donor_names),
                "input_clip_quantile": float(
                    pm_cfg.module_local_input_clip_quantile
                ),
                "input_clip": float(module_local_input_clip),
                "output_cap_artifact_path": (
                    module_local_cap_artifact.source_path
                ),
                "output_cap_artifact_sha256": (
                    module_local_cap_artifact.artifact_sha256
                ),
                "output_cap_source_checkpoint_sha256": (
                    module_local_cap_artifact.checkpoint_sha256
                ),
                "output_cap_quantile": float(
                    module_local_cap_artifact.absolute_quantile
                ),
                "output_cap": float(module_local_output_cap),
            }
            if module_local_graph_artifact is not None:
                precision_module_local_provenance.update(
                    {
                        "compartment_graph_path": (
                            module_local_graph_artifact.source_path
                        ),
                        "compartment_graph_sha256": (
                            module_local_graph_artifact.artifact_sha256
                        ),
                        "compartment_graph_schema_version": (
                            "kmlee_bam.module_local_compartment_graph.v1"
                        ),
                        "compartment_graph_topk": int(
                            module_local_graph_artifact.topk
                        ),
                        "compartment_graph_minimum_jaccard": float(
                            module_local_graph_artifact.minimum_jaccard
                        ),
                        "compartment_branch_variant": str(
                            pm_cfg.module_local_nonlinear_variant
                        ),
                    }
                )
            system.precision_module_local_provenance = (
                precision_module_local_provenance
            )
            system.precision_head.module_local_provenance = (
                precision_module_local_provenance
            )
        if int(os.environ.get("RANK", "0")) == 0:
            counts = system.precision_head.branch_parameter_count()
            print(
                "[PRISM-E2E] attached explicit score decomposition "
                f"axes=Thal,Braak,CERAD,LATE,Lewy ADNC_INPUT=false "
                f"contexts={len(table.context_names)} source_cells={table.source_cells:,} "
                f"params={counts}",
                flush=True,
            )
            if interaction_basis is not None:
                supported = [
                    name
                    for name, keep in zip(
                        interaction_basis.pair_names,
                        interaction_basis.eligible.detach().cpu().tolist(),
                    )
                    if bool(keep)
                ]
                print(
                    "[PRISM-E2E] pathology interactions "
                    f"supported={supported} "
                    f"complete_donors={interaction_basis.complete_donors.tolist()} "
                    "residual_fraction="
                    f"{[round(float(x), 4) for x in interaction_basis.residual_fraction.tolist()]} "
                    f"rank={int(pm_cfg.interaction_rank)}",
                    flush=True,
                )
            if module_local_artifact is not None:
                print(
                    "[PRISM module-local] "
                    f"artifact={module_local_artifact.artifact_sha256} "
                    f"rank={int(pm_cfg.module_local_rank)} "
                    f"fit_donors={len(fit_donor_names)} "
                    f"input_clip={float(module_local_input_clip):.6g} "
                    f"output_cap={float(module_local_output_cap):.6g} "
                    f"cap_checkpoint={module_local_cap_artifact.checkpoint_sha256}",
                    flush=True,
                )
                if module_local_graph_artifact is not None:
                    print(
                        "[PRISM paper branch] "
                        "source='ref/Dendritic morphology and synaptic nonlinearities "
                        "enhancefunctional complexity in human cortical neurons.pdf' "
                        "interpretation=computational_analogy_not_observed_dendrite_or_NMDA "
                        f"graph={module_local_graph_artifact.artifact_sha256} "
                        f"topk={module_local_graph_artifact.topk} "
                        f"variant={pm_cfg.module_local_nonlinear_variant} "
                        f"start={int(pm_cfg.module_local_nonlinear_start_epoch)} "
                        f"ramp={int(pm_cfg.module_local_nonlinear_ramp_epochs)}",
                        flush=True,
                    )

    return system


def build_criterion(cfg: RunConfig) -> TotalLoss:
    return TotalLoss(
        beta_state=cfg.loss.beta_state,
        lambda_tech=cfg.loss.lambda_tech,
        lambda_tech_mag=cfg.loss.lambda_tech_mag,
        lambda_gauge=cfg.loss.lambda_gauge,
        lambda_bam=cfg.loss.lambda_bam,
        lambda_cls=cfg.loss.lambda_cls,
        lambda_ref_center=cfg.loss.lambda_ref_center,
        lambda_ref_state=cfg.loss.lambda_ref_state,
        lambda_align=cfg.loss.lambda_align,
        lambda_white=cfg.loss.lambda_white,
        lambda_align_sex=cfg.loss.lambda_align_sex,
        align_sex_all_cells=cfg.loss.align_sex_all_cells,
        align_use_covariance=cfg.loss.align_use_covariance,
        align_reference_only=cfg.loss.align_reference_only,
        align_whiten_all_cells=cfg.loss.align_whiten_all_cells,
        align_cov_use_ema=cfg.loss.align_cov_use_ema,
        align_ema_decay=cfg.loss.align_ema_decay,
        ref_center_min_cells=cfg.loss.ref_center_min_cells,
        tau=cfg.loss.tau,
        r_min=cfg.loss.r_min,
        r_max=cfg.loss.r_max,
        detach_bam_uncertainty=cfg.loss.detach_bam_uncertainty,
        normalize_bam_weights=cfg.loss.normalize_bam_weights,
        normalize_gene_regularizers=cfg.loss.normalize_gene_regularizers,
    )


def build_optimizer(cfg: RunConfig, system: nn.Module) -> torch.optim.Optimizer:
    return torch.optim.AdamW(
        system.parameters(),
        lr=cfg.optim.lr,
        weight_decay=cfg.optim.weight_decay,
        betas=cfg.optim.betas,
        eps=cfg.optim.eps,
    )


def build_scheduler(cfg: RunConfig, optimizer: torch.optim.Optimizer):
    name = cfg.scheduler.name.lower()
    if name == "none":
        return None
    if name == "cosine":
        return torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer,
            T_max=cfg.scheduler.t_max,
            eta_min=cfg.scheduler.eta_min,
        )
    if name == "step":
        return torch.optim.lr_scheduler.StepLR(
            optimizer,
            step_size=cfg.scheduler.step_size,
            gamma=cfg.scheduler.gamma,
        )
    raise ValueError(f"Unknown scheduler name '{cfg.scheduler.name}'.")


# =====================================================================
# Checkpoint helpers
# =====================================================================
def metric_better(
    new_val: float,
    best_val: Optional[float],
    mode: str,
    *,
    min_delta: float = 0.0,
) -> bool:
    if best_val is None:
        return True
    if min_delta < 0:
        raise ValueError(f"min_delta must be non-negative, got {min_delta}.")
    if mode == "min":
        return new_val < (best_val - min_delta)
    if mode == "max":
        return new_val > (best_val + min_delta)
    raise ValueError(f"metric_mode must be 'min' or 'max', got '{mode}'.")


def joint_generator_selection_feasible(
    cfg: Any,
    metrics: Mapping[str, float],
) -> tuple[bool, dict[str, float]]:
    """Require paired validation non-inferiority before a joint-count win.

    The augmented-Lagrangian encourages feasibility but does not guarantee it
    at every epoch.  A smaller architecture may therefore compete for
    ``checkpoint_best.pt`` only when both normalized validation violations are
    finite and non-positive.  Legacy and gate-only runs are unchanged.
    """

    learned = getattr(cfg, "learned_generator_count", None)
    enabled = bool(
        learned is not None
        and getattr(learned, "enabled", False)
        and str(getattr(learned, "mode", "gate_only")) == "joint"
    )
    if not enabled:
        return True, {}
    values = {
        name: float(
            metrics.get(f"metric/generator_violation_{name}", math.nan)
        )
        for name in ("full_nll", "isolated_nll")
    }
    feasible = all(
        math.isfinite(value) and value <= 0.0 for value in values.values()
    )
    return feasible, values


def _criterion_short_name(name: str) -> str:
    """
    Auto-derive a short, filename-safe label from a metric key.

    Strips common namespace prefixes ('loss/', 'metric/', 'weight/', 'unc/')
    and replaces '/' and '.' with '_'. Example:
        "loss/total"                       -> "total"
        "metric/ordinal_nonzero_acc"       -> "ordinal_nonzero_acc"
        "loss/rec"                         -> "rec"
    """
    for prefix in ("loss/", "metric/", "weight/", "unc/"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    return name.replace("/", "_").replace(".", "_")


def normalise_early_stopping_criteria(
    *,
    primary_name: str,
    primary_mode: str,
    primary_min_delta: float,
    extra: list,
) -> list[dict]:
    """
    Return the full list of early-stopping criteria, with the primary first.
    Duplicate names (same metric listed both as primary and in extras) are
    deduplicated so each metric is tracked exactly once. Each entry is a dict
    with normalised keys: name, mode, min_delta, weight, short_name.

    `weight` defaults to 1.0 (used by composite save_best_mode).
    `short_name` is auto-derived from name if not provided (filename label).
    """
    out: list[dict] = []
    seen: set[str] = set()

    def _push(
        name: str,
        mode: str,
        min_delta: float,
        weight: float = 1.0,
        short_name: Optional[str] = None,
    ) -> None:
        if name in seen:
            return
        if mode not in ("min", "max"):
            raise ValueError(
                f"early stopping criterion '{name}' has bad mode '{mode}'; "
                "must be 'min' or 'max'."
            )
        if float(min_delta) < 0:
            raise ValueError(
                f"early stopping criterion '{name}' has negative min_delta."
            )
        if float(weight) < 0:
            raise ValueError(
                f"early stopping criterion '{name}' has negative weight."
            )
        short = short_name if short_name else _criterion_short_name(name)
        out.append({
            "name": str(name),
            "mode": str(mode),
            "min_delta": float(min_delta),
            "weight": float(weight),
            "short_name": str(short),
        })
        seen.add(name)

    # Primary may also appear in `extra` to override its weight / short_name.
    # If so, use those values when registering the primary so the user can
    # treat primary like any other criterion in the JSON list.
    primary_override: Optional[dict] = None
    for item in extra or []:
        if isinstance(item, dict) and item.get("name") == primary_name:
            primary_override = item
            break

    if primary_override is not None:
        _push(
            name=primary_name,
            mode=str(primary_override.get("mode", primary_mode)),
            min_delta=float(primary_override.get("min_delta", primary_min_delta)),
            weight=float(primary_override.get("weight", 1.0)),
            short_name=primary_override.get("short_name"),
        )
    else:
        _push(primary_name, primary_mode, primary_min_delta)

    for item in extra or []:
        if not isinstance(item, dict):
            raise ValueError(
                "early_stopping_criteria items must be dicts with keys "
                "'name', 'mode', 'min_delta'."
            )
        if "name" not in item:
            raise ValueError("early_stopping_criteria item missing 'name'.")
        if item.get("name") == primary_name:
            continue  # already pushed (with override if provided)
        _push(
            name=item["name"],
            mode=str(item.get("mode", "min")),
            min_delta=float(item.get("min_delta", 0.0)),
            weight=float(item.get("weight", 1.0)),
            short_name=item.get("short_name"),
        )
    return out


class _WelfordStats:
    """
    Online mean/variance tracker (Welford's algorithm) per metric name.

    Provides leak-free z-scores: callers should compute z using the current
    state, *then* call update() with this epoch's values. This guarantees
    the current epoch is not in its own normalisation statistics.
    """

    def __init__(self) -> None:
        # ``n`` remains the number of epoch updates for diagnostics/backward
        # compatibility.  Normalisation must use a separate finite-observation
        # count per metric: newly introduced or conditionally absent metrics
        # must not treat earlier NaNs as zero-valued observations.
        self.n: int = 0
        self.counts: Dict[str, int] = {}
        self.means: Dict[str, float] = {}
        self.M2s: Dict[str, float] = {}

    def std(self, name: str, *, eps: float) -> float:
        count = int(self.counts.get(name, 0))
        if count < 2:
            return float(eps)
        var = self.M2s.get(name, 0.0) / float(count - 1)
        return float(max(math.sqrt(max(var, 0.0)), eps))

    def z(self, name: str, value: float, *, eps: float) -> float:
        """z-score using state BEFORE this value is observed."""
        if math.isnan(value):
            return 0.0
        if int(self.counts.get(name, 0)) < 2:
            return 0.0
        mean = self.means.get(name, value)
        return (value - mean) / self.std(name, eps=eps)

    def update(self, values: Dict[str, float]) -> None:
        """Welford update; ignore NaN values for that metric only."""
        self.n += 1
        for name, val in values.items():
            if math.isnan(val):
                continue
            count = int(self.counts.get(name, 0)) + 1
            old_mean = self.means.get(name, 0.0)
            new_mean = old_mean + (val - old_mean) / float(count)
            self.counts[name] = count
            self.means[name] = new_mean
            self.M2s[name] = self.M2s.get(name, 0.0) + (val - old_mean) * (val - new_mean)


def compute_composite_score(
    *,
    criteria: list[dict],
    values: Dict[str, float],
    stats: _WelfordStats,
    eps: float,
) -> tuple[float, Dict[str, float], Dict[str, float]]:
    """
    Compute weighted z-score composite (higher = better).

    Returns
    -------
    composite_score : float
    z_per_criterion : dict {name: z_score}
    contrib_per_criterion : dict {name: weight * sign * z}

    Sign is +1 for mode="max" criteria and -1 for "min" so that improvement
    always increases the composite.
    """
    total = 0.0
    z_map: Dict[str, float] = {}
    contrib_map: Dict[str, float] = {}
    for crit in criteria:
        name = crit["name"]
        val = float(values.get(name, math.nan))
        if math.isnan(val):
            z_map[name] = 0.0
            contrib_map[name] = 0.0
            continue
        z = stats.z(name, val, eps=eps)
        sign = -1.0 if crit["mode"] == "min" else 1.0
        weight = float(crit.get("weight", 1.0))
        contrib = weight * sign * z
        z_map[name] = z
        contrib_map[name] = contrib
        total += contrib
    return total, z_map, contrib_map


def save_checkpoint(
    out_path: str | Path,
    *,
    epoch: int,
    cfg: RunConfig,
    system: Any,
    criterion: TotalLoss,
    optimizer: torch.optim.Optimizer,
    scheduler: Optional[Any],
    history: list[Dict[str, Dict[str, float]]],
    best_metric: Optional[float],
    best_epoch: Optional[int],
    stop_reason: Optional[str] = None,
    global_step: Optional[int] = None,
    step_in_epoch: Optional[int] = None,
    trainer: Optional[Any] = None,
) -> None:
    """Save a complete training checkpoint.

    `global_step` and `step_in_epoch` are optional so the same function can
    serve epoch checkpoints, step checkpoints, and interrupt checkpoints.

    If `trainer` exposes a `v7a_state_dict()` method (i.e. it's a V7aTrainer
    with PHU enabled), its profile / bank / ANCOVA state is also persisted
    under the `v7a_state` key. This lets resume preserve the difficulty
    profile and bank EMA across runs.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    base_system = unwrap_system(system)
    fixed_generator_integrity = run_system_integrity_checks(
        base_system,
        stage=f"checkpoint_save:{out_path.name}",
        distributed=False,
    )
    payload = {
        "epoch": epoch,
        "global_step": global_step,
        "step_in_epoch": step_in_epoch,
        "config": asdict(cfg),
        "system_state_dict": base_system.state_dict(),
        "criterion_state_dict": criterion.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": None if scheduler is None else scheduler.state_dict(),
        "history": history,
        "best_metric": best_metric,
        "best_epoch": best_epoch,
        "stop_reason": stop_reason,
    }
    if trainer is not None and hasattr(trainer, "scaler"):
        try:
            payload["amp_scaler_state_dict"] = trainer.scaler.state_dict()
        except Exception as exc:
            print(
                f"[checkpoint] WARNING: failed to save AMP scaler state: {exc}",
                flush=True,
            )
    fixed_generator_provenance = getattr(
        base_system,
        "fixed_generator_result_provenance",
        None,
    )
    if fixed_generator_provenance is not None:
        payload["fixed_generator_result_provenance"] = dict(
            fixed_generator_provenance
        )
        payload["fixed_generator_integrity"] = fixed_generator_integrity
    canonical_w44_source_provenance = getattr(
        base_system,
        "canonical_w44_source_provenance",
        None,
    )
    if canonical_w44_source_provenance is not None:
        if not isinstance(canonical_w44_source_provenance, dict):
            raise RuntimeError(
                "canonical W44 source checkpoint provenance must be an object"
            )
        payload["canonical_w44_source_provenance"] = dict(
            canonical_w44_source_provenance
        )

    # Optional v7a state (profile, bank, ANCOVA) for PHU runs.
    if trainer is not None and hasattr(trainer, "v7a_state_dict"):
        try:
            v7a_state = trainer.v7a_state_dict()
            if v7a_state:
                payload["v7a_state"] = v7a_state
        except Exception as exc:
            print(
                f"[checkpoint] WARNING: failed to save v7a_state: {exc}",
                flush=True,
            )

    # Optional v8 state (adaptive CVaR table + gene detection rate) for
    # zero-origin runs. Preserves the CVaR trigger/hysteresis state on resume.
    if trainer is not None and hasattr(trainer, "v8_state_dict"):
        try:
            v8_state = trainer.v8_state_dict()
            if v8_state:
                payload["v8_state"] = v8_state
        except Exception as exc:
            print(
                f"[checkpoint] WARNING: failed to save v8_state: {exc}",
                flush=True,
            )

    # Integrated curriculum position and live output scales. Although the
    # schedule is deterministic from epoch, persisting the resolved state
    # makes continuation auditable and fail-safe.
    if trainer is not None and hasattr(trainer, "curriculum_state_dict"):
        curriculum_state = trainer.curriculum_state_dict()
        if curriculum_state:
            payload["curriculum_state"] = curriculum_state

    torch.save(payload, out_path)


def _restore_checkpoint_for_training(
    ckpt_path: str | Path,
    *,
    cfg: RunConfig,
    trainer: Any,
    out_dir: Path,
    compute_sha256: bool = True,
) -> Dict[str, Any]:
    """Restore an epoch-boundary checkpoint for continued optimization.

    This is intentionally stricter than a model-only warm start: system,
    criterion, optimizer, scheduler, PHU/CVaR sidecars, epoch numbering,
    global step, and optional history are restored together. The old output
    directory remains immutable; continuation must write to a new directory.
    """

    source = Path(ckpt_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"resume checkpoint does not exist: {source}")
    if source.parent.resolve() == out_dir.resolve() and not bool(
        cfg.train.allow_in_place_resume
    ):
        raise ValueError(
            "resume continuation must use a new train.out_dir so the source "
            "run remains immutable"
        )
    if bool(cfg.train.early_stopping):
        raise ValueError(
            "resume continuation currently requires train.early_stopping=false "
            "because per-criterion/composite stopping state is not checkpointed"
        )

    payload = torch.load(source, map_location="cpu", weights_only=False)
    source_epoch = int(payload.get("epoch", 0))
    if source_epoch < 1:
        raise ValueError(f"resume checkpoint has invalid epoch: {source_epoch}")
    if payload.get("step_in_epoch") is not None:
        raise ValueError(
            "only epoch-boundary checkpoints can be resumed; "
            f"step_in_epoch={payload.get('step_in_epoch')}"
        )
    if int(cfg.train.epochs) <= source_epoch:
        raise ValueError(
            f"train.epochs={cfg.train.epochs} must exceed source epoch "
            f"{source_epoch}"
        )

    base_system = unwrap_system(trainer.system)
    base_system.load_state_dict(payload["system_state_dict"], strict=True)
    trainer.criterion.load_state_dict(payload["criterion_state_dict"], strict=True)

    if bool(cfg.train.resume_optimizer):
        optimizer_state = payload.get("optimizer_state_dict")
        if optimizer_state is None:
            raise KeyError("resume checkpoint lacks optimizer_state_dict")
        trainer.optimizer.load_state_dict(optimizer_state)
        scheduler_state = payload.get("scheduler_state_dict")
        if trainer.scheduler is None:
            if scheduler_state is not None:
                raise ValueError(
                    "resume checkpoint contains scheduler state but current "
                    "configuration has no scheduler"
                )
        elif scheduler_state is None:
            raise ValueError(
                "current configuration has a scheduler but resume checkpoint "
                "does not contain scheduler state"
            )
        else:
            trainer.scheduler.load_state_dict(scheduler_state)

    scaler_state = payload.get("amp_scaler_state_dict")
    if scaler_state is not None and hasattr(trainer, "scaler"):
        trainer.scaler.load_state_dict(scaler_state)

    if (
        hasattr(trainer, "v7a_load_state_dict")
        and payload.get("v7a_state")
    ):
        trainer.v7a_load_state_dict(payload["v7a_state"])
    if (
        hasattr(trainer, "v8_load_state_dict")
        and payload.get("v8_state")
    ):
        trainer.v8_load_state_dict(payload["v8_state"])
    if (
        hasattr(trainer, "curriculum_load_state_dict")
        and payload.get("curriculum_state")
    ):
        trainer.curriculum_load_state_dict(payload["curriculum_state"])

    run_system_integrity_checks(
        base_system,
        stage=f"training_resume:{source.name}",
        distributed=False,
    )

    source_sha256 = None
    if compute_sha256:
        digest = hashlib.sha256()
        with source.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        source_sha256 = digest.hexdigest()
    receipt = {
        "schema_version": "kmlee_bam.training_resume.v1",
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_sha256,
        "source_epoch": source_epoch,
        "next_epoch": source_epoch + 1,
        "maximum_epoch": int(cfg.train.epochs),
        "source_global_step": int(payload.get("global_step") or 0),
        "optimizer_restored": bool(cfg.train.resume_optimizer),
        "history_restored": bool(cfg.train.resume_history),
        "v7a_state_restored": bool(payload.get("v7a_state")),
        "v8_state_restored": bool(payload.get("v8_state")),
        "curriculum_state_restored": bool(payload.get("curriculum_state")),
        "official_test_used": False,
    }
    return {"payload": payload, "receipt": receipt}


def save_step_checkpoint(
    *,
    out_dir: Path,
    epoch: int,
    step_in_epoch: int,
    global_step: int,
    cfg: RunConfig,
    trainer: Trainer,
    history: list[Dict[str, Dict[str, float]]],
    best_metric: Optional[float],
    best_epoch: Optional[int],
    keep_last: int = 3,
) -> Path:
    """Save a rolling step checkpoint and keep only the most recent N."""
    ckpt_dir = out_dir / "step_checkpoints"
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    ckpt_path = ckpt_dir / f"checkpoint_step_{global_step:09d}.pt"
    save_checkpoint(
        ckpt_path,
        epoch=epoch,
        global_step=global_step,
        step_in_epoch=step_in_epoch,
        cfg=cfg,
        system=trainer.system,
        criterion=trainer.criterion,
        optimizer=trainer.optimizer,
        scheduler=trainer.scheduler,
        trainer=trainer,
        history=history,
        best_metric=best_metric,
        best_epoch=best_epoch,
        stop_reason=None,
    )

    if keep_last is not None and keep_last > 0:
        ckpts = sorted(ckpt_dir.glob("checkpoint_step_*.pt"))
        excess = len(ckpts) - int(keep_last)
        if excess > 0:
            for old in ckpts[:excess]:
                try:
                    old.unlink()
                except OSError:
                    pass

    return ckpt_path


def _load_checkpoint_for_analysis(
    ckpt_path: str | Path,
    *,
    system: Any,
    criterion: TotalLoss,
    device: str | torch.device,
    trainer: Optional[Any] = None,
) -> Dict[str, Any]:
    payload = torch.load(Path(ckpt_path), map_location=device)
    base_system = unwrap_system(system)
    base_system.load_state_dict(payload["system_state_dict"], strict=True)
    run_system_integrity_checks(
        base_system,
        stage=f"checkpoint_reload:{Path(ckpt_path).name}",
        distributed=False,
    )
    criterion.load_state_dict(payload["criterion_state_dict"])
    base_system.to(device)
    criterion.to(device)
    base_system.eval()
    criterion.eval()

    # Optional v7a state restoration. We only attempt this if the trainer
    # has v7a_load_state_dict (i.e. it's a V7aTrainer) and the checkpoint
    # actually contains a saved v7a_state payload.
    if (
        trainer is not None
        and hasattr(trainer, "v7a_load_state_dict")
        and "v7a_state" in payload
        and payload["v7a_state"]
    ):
        try:
            trainer.v7a_load_state_dict(payload["v7a_state"])
            print("[checkpoint] v7a state (profile/bank/ANCOVA) restored from checkpoint.")
        except Exception as exc:
            print(
                f"[checkpoint] WARNING: failed to restore v7a_state: {exc}",
                flush=True,
            )

    # Optional v8 state restoration (adaptive CVaR + gene detection rate).
    if (
        trainer is not None
        and hasattr(trainer, "v8_load_state_dict")
        and "v8_state" in payload
        and payload["v8_state"]
    ):
        try:
            trainer.v8_load_state_dict(payload["v8_state"])
            print("[checkpoint] v8 state (CVaR/gene-detect-rate) restored from checkpoint.")
        except Exception as exc:
            print(
                f"[checkpoint] WARNING: failed to restore v8_state: {exc}",
                flush=True,
            )

    if (
        trainer is not None
        and hasattr(trainer, "curriculum_load_state_dict")
        and payload.get("curriculum_state")
    ):
        trainer.curriculum_load_state_dict(payload["curriculum_state"])
        print("[checkpoint] integrated curriculum state restored from checkpoint.")

    return payload


def _save_split_latents(
    *,
    trainer: Trainer,
    out_dir: Path,
    train_loader: DataLoader,
    val_loader: Optional[DataLoader],
    test_loader: Optional[DataLoader],
    prefix: str,
) -> None:
    split_to_loader = {
        "train": train_loader,
        "val": val_loader,
        "test": test_loader,
    }
    for split_name, loader in split_to_loader.items():
        if loader is None:
            continue
        latents = trainer.collect_latents(loader)
        torch.save(latents, out_dir / f"latents_{prefix}_{split_name}.pt")


def _write_training_summary(
    *,
    out_dir: Path,
    stopped_epoch: int,
    best_epoch: Optional[int],
    best_metric: Optional[float],
    monitored_split: str,
    monitored_name: str,
    stop_reason: str,
    used_best_checkpoint_for_test: bool,
    runtime: RuntimeContext,
    selection_mode: str = "primary",
    primary_best_epoch: Optional[int] = None,
    primary_best_metric: Optional[float] = None,
    best_per_criterion: Optional[Dict[str, Optional[float]]] = None,
    bad_epochs_at_stop: int = 0,
    patience: int = 0,
) -> None:
    summary = {
        "stopped_epoch": stopped_epoch,
        "best_epoch": best_epoch,
        "best_metric": best_metric,
        "monitored_split": monitored_split,
        "monitored_name": monitored_name,
        "stop_reason": stop_reason,
        "used_best_checkpoint_for_test": used_best_checkpoint_for_test,
        "selection_mode": selection_mode,
        "selected_best_epoch": best_epoch,
        "selected_best_score": best_metric,
        "primary_best_epoch": primary_best_epoch,
        "primary_best_metric": primary_best_metric,
        "best_per_criterion": best_per_criterion or {},
        "bad_epochs_at_stop": int(bad_epochs_at_stop),
        "patience": int(patience),
        "distributed": runtime.distributed,
        "world_size": runtime.world_size,
        "backend": runtime.backend,
    }
    with open(out_dir / "training_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)




def _format_seconds(sec: float) -> str:
    sec = float(sec)
    if sec < 60:
        return f"{sec:.1f}s"
    minutes, seconds = divmod(sec, 60)
    if minutes < 60:
        return f"{int(minutes)}m{int(seconds)}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes)}m{int(seconds)}s"


def _format_optional_float(value: Optional[float], digits: int = 6) -> str:
    """Format an optional metric without crashing before its first update."""

    numeric = float("nan") if value is None else float(value)
    return f"{numeric:.{int(digits)}f}"


# =====================================================================
# Main run
# =====================================================================
def validate_training_control_contract(cfg: RunConfig) -> None:
    """Validate opt-in, validation-only stopping before runtime setup."""

    if not bool(cfg.train.require_validation_for_early_stopping):
        return
    if not bool(cfg.train.early_stopping):
        raise ValueError(
            "require_validation_for_early_stopping=true requires "
            "early_stopping=true"
        )
    if cfg.train.eval_every is None or int(cfg.train.eval_every) != 1:
        raise ValueError(
            "validation-only early stopping requires eval_every=1"
        )


def run_training(cfg: RunConfig) -> None:
    """Run training with DDP-safe rank-0 validation, step checkpoints, ETA, and interrupt saving."""
    validate_training_control_contract(cfg)
    runtime = setup_runtime(cfg)
    boot_debug = os.environ.get("KMLEE_TRAIN_BOOT_DEBUG", "0") == "1"

    def _boot_log(message: str) -> None:
        if boot_debug:
            print(f"[boot-debug rank={runtime.rank}] {message}", flush=True)

    trainer: Optional[Trainer] = None
    history: list[Dict[str, Dict[str, float]]] = []
    best_metric: Optional[float] = None
    best_epoch: Optional[int] = None
    bad_epochs: int = 0
    stopped_epoch: int = 0
    monitored_split: str = "train"
    stop_reason: str = "completed_all_epochs"
    global_step: int = 0
    next_step_checkpoint: Optional[int] = cfg.train.save_every_steps
    start_epoch: int = 1

    # Multi-criteria early stopping bookkeeping. Each criterion gets its own
    # rolling best so that the patience counter only advances when NONE
    # improves. Falls back to single-metric behaviour when the extra list
    # is empty.
    es_criteria = normalise_early_stopping_criteria(
        primary_name=cfg.train.metric_name,
        primary_mode=cfg.train.metric_mode,
        primary_min_delta=cfg.train.early_stopping_min_delta,
        extra=list(cfg.train.early_stopping_criteria or []),
    )
    es_best: dict[str, Optional[float]] = {c["name"]: None for c in es_criteria}

    # Composite save_best bookkeeping. Welford stats track per-metric mean/std
    # across epochs; z-scores at epoch t use stats from epochs [0..t-1] (leak
    # free). Composite is only used to update `checkpoint_best.pt` once the
    # `early_stopping_min_epochs` warm-up has passed.
    composite_stats = _WelfordStats()
    best_composite: Optional[float] = None
    best_composite_epoch: Optional[int] = None

    try:
        _boot_log("set_seed:start")
        set_seed(cfg.train.seed)
        _boot_log("set_seed:done")

        out_dir = Path(cfg.train.out_dir)
        if runtime.is_main_process:
            _boot_log("out_dir_write:start")
            out_dir.mkdir(parents=True, exist_ok=True)
            with open(out_dir / "resolved_config.json", "w", encoding="utf-8") as f:
                json.dump(asdict(cfg), f, indent=2, ensure_ascii=False)
            _boot_log("out_dir_write:done")
        _boot_log("barrier_after_config:start")
        barrier(runtime)
        _boot_log("barrier_after_config:done")

        _boot_log("build_datasets:start")
        train_ds, val_ds, test_ds = build_datasets(cfg)
        _boot_log("build_datasets:done")
        if (
            bool(cfg.train.require_validation_for_early_stopping)
            and val_ds is None
        ):
            raise RuntimeError(
                "validation-only training control resolved no validation dataset"
            )
        architecture_train_ds = None
        architecture_split = None
        if bool(
            cfg.learned_generator_count.enabled
            and str(cfg.learned_generator_count.mode) == "gate_only"
        ):
            if bool(cfg.generator_budget_search.enabled):
                raise ValueError(
                    "learned_generator_count and generator_budget_search "
                    "cannot be enabled in the same run"
                )
            full_train_ds = train_ds
            architecture_split = make_architecture_donor_split(
                full_train_ds,
                architecture_donor_count=int(
                    cfg.learned_generator_count.architecture_donor_count
                ),
                seed=int(
                    cfg.learned_generator_count.architecture_split_seed
                ),
                search_trials=int(
                    cfg.learned_generator_count.architecture_split_trials
                ),
            )
            train_ds = DonorIndexView(
                full_train_ds,
                architecture_split.weight_donor_ids,
            )
            architecture_train_ds = DonorIndexView(
                full_train_ds,
                architecture_split.architecture_donor_ids,
            )
            if runtime.is_main_process:
                split_path = out_dir / "architecture_donor_split.json"
                with split_path.open("w", encoding="utf-8") as handle:
                    json.dump(
                        architecture_split.manifest,
                        handle,
                        indent=2,
                        ensure_ascii=False,
                        sort_keys=True,
                    )
                    handle.write("\n")
                print(
                    "[learned-generator split] "
                    f"weight_train={len(architecture_split.weight_donor_ids)} "
                    f"architecture_train="
                    f"{len(architecture_split.architecture_donor_ids)} "
                    f"weight_cells={len(train_ds):,} "
                    f"architecture_cells={len(architecture_train_ds):,} "
                    f"balance_score="
                    f"{architecture_split.manifest['balance_score']:.6f} "
                    f"manifest={split_path}",
                    flush=True,
                )
            barrier(runtime)
        _boot_log("build_loaders:start")
        (
            train_loader,
            val_loader,
            test_loader,
            train_sampler,
            _val_sampler,
            _test_sampler,
        ) = build_loaders(
            cfg,
            train_ds,
            val_ds,
            test_ds,
            distributed=runtime.distributed,
        )
        phase2_train_loader = None
        phase2_train_sampler = None
        phase_cfg = cfg.integrated_phase_curriculum
        if (
            bool(phase_cfg.enabled)
            and int(phase_cfg.phase2_batch_size) != int(cfg.loader.batch_size)
        ):
            original_batch_size = int(cfg.loader.batch_size)
            try:
                cfg.loader.batch_size = int(phase_cfg.phase2_batch_size)
                (
                    phase2_train_loader,
                    _phase2_val_loader,
                    _phase2_test_loader,
                    phase2_train_sampler,
                    _phase2_val_sampler,
                    _phase2_test_sampler,
                ) = build_loaders(
                    cfg,
                    train_ds,
                    val_ds,
                    test_ds,
                    distributed=runtime.distributed,
                )
            finally:
                cfg.loader.batch_size = original_batch_size
        _boot_log("build_loaders:done")

        rank0_val_loader = None
        rank0_test_loader = None

        if bool(cfg.train.use_rank0_eval):
            if val_ds is not None and (runtime.is_main_process or not runtime.distributed):
                rank0_val_loader = build_plain_eval_loader(cfg, val_ds, shuffle=False)

            if test_ds is not None and (runtime.is_main_process or not runtime.distributed):
                rank0_test_loader = build_plain_eval_loader(cfg, test_ds, shuffle=False)


        _boot_log("build_system:start")
        system = build_system(cfg, train_ds)
        _boot_log("build_system:done")
        _copy_decoder_architecture_metadata_to_system(system)
        if runtime.is_main_process:
            _boot_log("decoder_architecture_manifest:start")
            architecture_manifest = (
                _build_resolved_decoder_architecture_manifest(cfg, system)
            )
            with open(
                out_dir / "resolved_decoder_architecture.json",
                "w",
                encoding="utf-8",
            ) as f:
                json.dump(
                    architecture_manifest,
                    f,
                    indent=2,
                    ensure_ascii=False,
                )
            print(
                "[decoder-architecture] wrote "
                f"{out_dir / 'resolved_decoder_architecture.json'} "
                f"strategy="
                f"{architecture_manifest['decoder_mask_strategy']} "
                f"registry_modules="
                f"{architecture_manifest['registry_module_count']} "
                f"generators="
                f"{architecture_manifest['decoder_generator_count']} "
                f"singletons="
                f"{architecture_manifest['singleton_generator_count']} "
                f"covered_genes="
                f"{architecture_manifest['covered_gene_count']}/"
                f"{architecture_manifest['registry_gene_count']} "
                f"nnz="
                f"{architecture_manifest['selected_membership_nnz']} "
                f"active_rank_components="
                f"{architecture_manifest['active_rank_component_count']} "
                f"selection_sha256="
                f"{architecture_manifest['selection_sha256']} "
                f"trainable_params="
                f"{architecture_manifest['trainable_parameter_count']:,}",
                flush=True,
            )
            _boot_log("decoder_architecture_manifest:done")
        _boot_log("build_criterion:start")
        criterion = build_criterion(cfg)
        _boot_log("build_criterion:done")

        _boot_log("build_optimizer:start")
        optimizer = build_optimizer(cfg, system)
        _boot_log("build_optimizer:done")
        _boot_log("build_scheduler:start")
        scheduler = build_scheduler(cfg, optimizer)
        _boot_log("build_scheduler:done")

        _boot_log("trainer_init:start")
        trainer = Trainer(
            system=system,
            criterion=criterion,
            optimizer=optimizer,
            device=runtime.device,
            scheduler=scheduler,
            grad_clip_norm=cfg.train.grad_clip_norm,
            amp=cfg.train.amp,
            amp_dtype=cfg.train.amp_dtype,
            tech_weight_mode=cfg.train.tech_weight_mode,
            scheduler_step_on=cfg.train.scheduler_step_on,
            eval_sample_latent=cfg.train.eval_sample_latent,
            grad_accum_steps=cfg.train.grad_accum_steps,
            debug_z_sensitivity_every=cfg.train.debug_z_sensitivity_every,
            debug_z_sensitivity_path=(
                out_dir / "z_debug.jsonl"
                if runtime.is_main_process and cfg.train.debug_z_sensitivity_every
                else None
            ),
        )
        if bool(cfg.integrated_phase_curriculum.enabled):
            trainer.integrated_phase_controller = (
                IntegratedPhaseCurriculumController(
                    cfg.integrated_phase_curriculum
                )
            )
            trainer.pathology_aux_curriculum_scale = float(
                cfg.integrated_phase_curriculum.phase1_pathology_scale
            )
            unwrap_system(trainer.system).precision_forward_enabled = False
        # Read-only schedule metadata for the integrated-PRISM console.  This
        # avoids hard-coded epoch numbers in the renderer and makes resumed or
        # compressed smoke curricula describe themselves correctly.
        _schedule_system = unwrap_system(trainer.system)
        _schedule_precision = getattr(
            getattr(_schedule_system, "precision_head", None), "cfg", None
        )
        _schedule_phu = getattr(trainer, "v7a_alignment_config", None)
        _schedule_generator = getattr(
            _schedule_system, "generator_count_config", None
        )
        trainer.integrated_curriculum_schedule = {
            "total_epochs": int(cfg.train.epochs),
            "phu_alignment_start": int(
                getattr(_schedule_phu, "alignment_start_epoch", 1)
            ),
            "explicit_start": int(
                getattr(_schedule_precision, "output_start_epoch", 1)
            ),
            "module_rescue_start": int(cfg.prism_module_rescue.start_epoch),
            "phu_relative_start": int(
                getattr(_schedule_phu, "relative_start_epoch", 1)
            ),
            "generator_start": int(
                getattr(_schedule_generator, "start_epoch", 1)
            ),
        }
        trainer._integrated_last_module_rescue_metrics = None
        _boot_log("trainer_init:done")

        if cfg.train.resume_checkpoint:
            _boot_log("training_resume:start")
            restored = _restore_checkpoint_for_training(
                cfg.train.resume_checkpoint,
                cfg=cfg,
                trainer=trainer,
                out_dir=out_dir,
                compute_sha256=runtime.is_main_process,
            )
            resume_payload = restored["payload"]
            start_epoch = int(resume_payload["epoch"]) + 1
            global_step = int(resume_payload.get("global_step") or 0)
            best_metric = resume_payload.get("best_metric")
            best_epoch = resume_payload.get("best_epoch")
            if bool(cfg.train.resume_history):
                history = list(resume_payload.get("history") or [])
            if (
                cfg.train.save_every_steps is not None
                and int(cfg.train.save_every_steps) > 0
            ):
                interval = int(cfg.train.save_every_steps)
                next_step_checkpoint = (global_step // interval + 1) * interval
            else:
                next_step_checkpoint = cfg.train.save_every_steps
            if runtime.is_main_process:
                receipt_path = out_dir / "resume_receipt.json"
                with receipt_path.open("w", encoding="utf-8") as handle:
                    json.dump(
                        restored["receipt"],
                        handle,
                        indent=2,
                        ensure_ascii=False,
                        sort_keys=True,
                    )
                    handle.write("\n")
                print(
                    "[training-resume] restored complete epoch-boundary state "
                    f"source_epoch={start_epoch - 1} next_epoch={start_epoch} "
                    f"global_step={global_step:,} "
                    f"optimizer={bool(cfg.train.resume_optimizer)} "
                    f"history={len(history)} receipt={receipt_path}",
                    flush=True,
                )
            barrier(runtime)
            _boot_log("training_resume:done")

        # The fixed-mask callback is absent for all ordinary and discovery
        # runs.  Exact fixed-architecture runs verify rank agreement after the
        # model has moved to its rank-local device and before DDP construction.
        run_system_integrity_checks(
            trainer.system,
            stage="pre_ddp",
            distributed=runtime.distributed,
        )

        # Expose the reference-scaler cache to tech-invariance (train_ds is in scope here;
        # the run_current global is unreliable across DDP ranks).
        if hasattr(trainer, "tech_inv_config"):
            try:
                trainer._tech_scaler_mean = torch.as_tensor(train_ds._scaler_mean)
                trainer._tech_scaler_std_eps = torch.as_tensor(train_ds._scaler_std_eps)
                trainer._tech_scaler_clip = getattr(train_ds, "x_gene_scalar_clip", None)
            except Exception as _e:
                print(f"[tech-inv] could not stash scaler on trainer: {_e!r}", flush=True)

        if runtime.distributed:
            skip_initial_sync = bool(cfg.train.ddp_skip_initial_sync)
            if runtime.is_main_process and skip_initial_sync:
                print(
                    "[ddp] constructor initial parameter verification/sync is disabled by config.",
                    flush=True,
                )
            _boot_log("ddp_wrap:start")
            with ddp_initial_sync_context(skip_initial_sync):
                ddp_model = DDP(
                    unwrap_system(trainer.system),
                    device_ids=[runtime.local_rank] if runtime.device.type == "cuda" else None,
                    output_device=runtime.local_rank if runtime.device.type == "cuda" else None,
                    find_unused_parameters=cfg.train.ddp_find_unused_parameters,
                    broadcast_buffers=False,
                )
            trainer.system = DDPSystemProxy(ddp_model)
            _boot_log("ddp_wrap:done")

        if runtime.distributed and hasattr(trainer, "set_ddp_control_group"):
            _boot_log("set_ddp_control_group:start")
            trainer.set_ddp_control_group(runtime.control_group)
            _boot_log("set_ddp_control_group:done")

        module_rescue_updater: Optional[PrismModuleRescueUpdater] = None
        if bool(cfg.prism_module_rescue.enabled):
            if not bool(cfg.precision_medicine.enabled):
                raise ValueError(
                    "prism_module_rescue requires precision_medicine.enabled"
                )
            weight_donor_ids = np.unique(
                np.asarray(train_ds.donor_ids, dtype=np.int64)[
                    np.asarray(train_ds.row_idx, dtype=np.int64)
                ]
            )
            registry_for_rescue = _load_registry_json(
                cfg.module_tokenizer.registry_json_path
            )
            module_rescue_updater = PrismModuleRescueUpdater(
                config=cfg.prism_module_rescue,
                dataset=train_ds,
                donor_ids=weight_donor_ids,
                module_names=registry_for_rescue["module_names"],
                spec_path=cfg.data.spec_path,
                registry_path=cfg.module_tokenizer.registry_json_path,
                activity_path=cfg.module_tokenizer.activity_weight_path,
                activity_normalization=(
                    cfg.module_tokenizer.activity_weight_normalization
                ),
            )
            trainer.module_rescue_updater = module_rescue_updater
            pending_rescue_state = getattr(
                trainer, "_pending_module_rescue_state", None
            )
            if pending_rescue_state is not None:
                module_rescue_updater.load_state_dict(pending_rescue_state)
                del trainer._pending_module_rescue_state
            module_rescue_updater.resolve_parameter_allowlist(trainer.system)
            # Validate the complete schedule on every rank before rank-0 I/O;
            # otherwise a rank-0-only validation exception can strand peers at
            # the following barrier.
            resolved_rescue_steps = module_rescue_updater.steps_per_rank(
                world_size=runtime.world_size
            )
            resolved_rescue_local_accum = (
                module_rescue_updater.local_blocks_per_optimizer_update(
                    world_size=runtime.world_size
                )
            )
            resolved_rescue_updates = (
                module_rescue_updater.optimizer_updates_per_epoch(
                    world_size=runtime.world_size
                )
            )
            if runtime.is_main_process:
                manifest_path = out_dir / "prism_module_rescue_manifest.json"
                write_module_rescue_manifest(
                    manifest_path,
                    config=cfg.prism_module_rescue,
                    updater=module_rescue_updater,
                    world_size=runtime.world_size,
                )
                print(
                    "[PRISM module-rescue] grouped auxiliary updates enabled "
                    f"weight_donors={len(weight_donor_ids)} "
                    f"celltypes="
                    f"{len(module_rescue_updater.sampler.eligible_celltypes)} "
                    f"steps/rank/epoch="
                    f"{resolved_rescue_steps} "
                    f"blocks/local-update="
                    f"{resolved_rescue_local_accum} "
                    f"updates/epoch="
                    f"{resolved_rescue_updates} "
                    f"cells/donor="
                    f"{cfg.prism_module_rescue.cells_per_donor_celltype} "
                    f"lambda={cfg.prism_module_rescue.lambda_module} "
                    f"ccc={cfg.prism_module_rescue.ccc_weight} "
                    f"deterministic={cfg.prism_module_rescue.deterministic_forward} "
                    f"restricted={cfg.prism_module_rescue.restrict_parameter_updates} "
                    "ordinary_minibatch_centering=false",
                    flush=True,
                )
            barrier(runtime)

        ccc_calibration_output = os.environ.get(
            "KMLEE_PRISM_CCC_CALIBRATION_OUT"
        )
        if ccc_calibration_output:
            if module_rescue_updater is None:
                raise RuntimeError(
                    "CCC calibration requires prism_module_rescue.enabled"
                )
            if runtime.world_size != 1:
                raise RuntimeError(
                    "CCC calibration is a deterministic single-GPU preflight"
                )
            report = module_rescue_updater.calibrate_ccc_gradient_weight(
                trainer,
                calibration_epoch=1,
                calibration_blocks=len(
                    module_rescue_updater.sampler.eligible_celltypes
                ),
            )
            destination = Path(ccc_calibration_output)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("w", encoding="utf-8") as handle:
                json.dump(
                    report,
                    handle,
                    indent=2,
                    ensure_ascii=False,
                    sort_keys=True,
                )
                handle.write("\n")
            print(
                "[PRISM module-rescue] CCC train-only calibration complete "
                f"weight={report['calibrated_ccc_weight']:.9g} "
                f"huber_grad_rms={report['huber_gradient_rms']:.9g} "
                f"ccc_grad_rms={report['ccc_gradient_rms']:.9g} "
                f"output={destination}",
                flush=True,
            )
            return

        generator_search: Optional[GeneratorBudgetSearch] = None
        if bool(cfg.generator_budget_search.enabled):
            if val_ds is None:
                raise RuntimeError(
                    "generator_budget_search requires a validation split."
                )
            if runtime.is_main_process:
                _boot_log("generator_budget_search:init:start")
                generator_search = GeneratorBudgetSearch(
                    config=cfg.generator_budget_search,
                    system=trainer.system,
                    validation_dataset=val_ds,
                    device=runtime.device,
                    output_dir=out_dir,
                )
                print(
                    "[generator-budget] hard backward-elimination search enabled "
                    f"candidates={generator_search.controller.num_generators} "
                    f"protected={len(generator_search.controller.protected_ids)} "
                    f"start_epoch={cfg.generator_budget_search.start_epoch} "
                    f"block={cfg.generator_budget_search.block_size} "
                    f"min_active="
                    f"{cfg.generator_budget_search.minimum_active_generators} "
                    "soft_gate=false l1_l0=false "
                    "views=full+generator-isolated",
                    flush=True,
                )
                _boot_log("generator_budget_search:init:done")
            if runtime.distributed:
                # This is also required on resume: rank 0 restores the
                # controller sidecar, then every worker must receive exactly
                # that hard architecture before the next training step.
                decoder = unwrap_system(trainer.system).decoder
                dist.broadcast(
                    decoder.generator_active_mask,
                    src=0,
                )
            barrier(runtime)

        steps_per_epoch = len(train_loader)
        updates_per_epoch = math.ceil(steps_per_epoch / max(1, cfg.train.grad_accum_steps))
        phase2_steps_per_epoch = (
            len(phase2_train_loader)
            if phase2_train_loader is not None
            else steps_per_epoch
        )
        phase2_updates_per_epoch = math.ceil(
            phase2_steps_per_epoch
            / max(1, int(cfg.integrated_phase_curriculum.phase2_grad_accum_steps))
        )

        if runtime.is_main_process:
            informative_console = (
                os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower()
                in {"prism_experiment", "prism_informative"}
            )
            if informative_console:
                test_text = (
                    "SEALED"
                    if cfg.train.skip_final_test
                    else f"{0 if test_ds is None else len(test_ds):,}"
                )
                effective_batch = (
                    int(cfg.loader.batch_size)
                    * int(cfg.train.grad_accum_steps)
                    * int(runtime.world_size)
                )
                print("═" * 92)
                print("KMLEE-BAM · PRISM training start")
                print(
                    "  [데이터] "
                    f"train {len(train_ds):,} · val {0 if val_ds is None else len(val_ds):,} · "
                    f"test {test_text} · genes {len(train_ds.spec.gene_names):,} · "
                    f"celltypes {len(train_ds.spec.celltype_vocab):,} · bins {train_ds.spec.n_bins}"
                )
                print(
                    "  [연산] "
                    f"GPU {runtime.world_size} · {runtime.backend} · {trainer.amp_dtype} · "
                    f"batch/GPU {cfg.loader.batch_size} × accum {cfg.train.grad_accum_steps} "
                    f"= effective {effective_batch}"
                )
                print(
                    "  [일정] "
                    f"epochs {cfg.train.epochs} · steps/epoch {steps_per_epoch:,} · "
                    f"updates/epoch {updates_per_epoch:,} · log every {cfg.train.log_every} · "
                    f"progress every {cfg.train.progress_log_every}"
                )
                print(
                    "  [검증] "
                    f"every {cfg.train.eval_every} epoch · split=val · "
                    f"metric={cfg.train.metric_name} · official test "
                    f"{'sealed' if cfg.train.skip_final_test else 'enabled'}"
                )
                print(f"  [출력] {out_dir}")
                print("═" * 92)
            else:
                print("=" * 72)
                print("KMLEE-BAM integrated training start")
                print("=" * 72)
                print(f"train cells       : {len(train_ds):,}")
                print(f"val cells         : {0 if val_ds is None else len(val_ds):,}")
                print(f"test cells        : {0 if test_ds is None else len(test_ds):,}")
                print(f"genes             : {len(train_ds.spec.gene_names):,}")
                print(f"bins              : {train_ds.spec.n_bins}")
                print(f"cell types        : {len(train_ds.spec.celltype_vocab):,}")
                print(f"tech groups       : {len(train_ds.spec.batch_vocab):,}")
                print(f"batch size        : {cfg.loader.batch_size}")
                print(f"grad accum        : {cfg.train.grad_accum_steps}")
                print(f"steps/epoch       : {steps_per_epoch:,}")
                print(f"updates/epoch     : {updates_per_epoch:,}")
                print(f"save_every_steps  : {cfg.train.save_every_steps}")
                print(f"progress_log_every: {cfg.train.progress_log_every}")
                print(f"debug_z_every     : {cfg.train.debug_z_sensitivity_every}")
                print(
                    "debug_z_path      : "
                    f"{out_dir / 'z_debug.jsonl' if cfg.train.debug_z_sensitivity_every else None}"
                )
                print(f"max_train_steps   : {cfg.train.max_train_steps_per_epoch}")
                print(f"max_eval_steps    : {cfg.train.max_eval_steps}")
                print(f"device            : {trainer.device}")
                print(f"amp enabled       : {trainer.amp_enabled}")
                print(f"amp dtype         : {trainer.amp_dtype}")
                print(f"distributed       : {runtime.distributed}")
                print(f"world_size        : {runtime.world_size}")
                print(f"backend           : {runtime.backend}")
                print(f"control_pg        : {runtime.control_backend}")
                print(f"use rank0 eval   : {cfg.train.use_rank0_eval}")
                if cfg.train.use_rank0_eval:
                    print(
                        "rank0 eval bs    : "
                        f"{cfg.train.rank0_eval_batch_size or cfg.loader.eval_batch_size or cfg.loader.batch_size}"
                    )
                    print(f"rank0 eval workers: {cfg.train.rank0_eval_num_workers}")
                    print(f"rank0 eval pin    : {cfg.train.rank0_eval_pin_memory}")
                    print(f"rank0 eval persist: {cfg.train.rank0_eval_persistent_workers}")
                    print(f"rank0 eval prefetch: {cfg.train.rank0_eval_prefetch_factor}")
                else:
                    print(f"eval bs           : {cfg.loader.eval_batch_size or cfg.loader.batch_size}")
                    print("rank0 eval        : disabled; using DDP/all-rank eval")
                print(f"eval_every        : {cfg.train.eval_every}")
                print("=" * 72)
            
        run_start_time = time.perf_counter()

        def make_step_callback(epoch: int, epoch_start_time: float):
            def _callback(*, step_idx, n_total_steps, step_out, meters, epoch_index=None):
                nonlocal global_step, next_step_checkpoint, best_metric, best_epoch, history
                global_step += 1

                # ------------------------------------------------------
                # Progress / ETA logging
                # ------------------------------------------------------
                should_progress_log = (
                    runtime.is_main_process
                    and cfg.train.progress_log_every is not None
                    and cfg.train.progress_log_every > 0
                    and step_idx % cfg.train.progress_log_every == 0
                )
                if should_progress_log:
                    elapsed_epoch = time.perf_counter() - epoch_start_time
                    sec_per_step = elapsed_epoch / max(step_idx, 1)
                    if n_total_steps is not None:
                        steps_left = max(0, int(n_total_steps) - int(step_idx))
                        eta_epoch = sec_per_step * steps_left
                        pct = 100.0 * step_idx / max(int(n_total_steps), 1)
                        elapsed_run = time.perf_counter() - run_start_time
                        print(
                            f"[progress epoch {epoch:03d}] "
                            f"step={step_idx:,}/{int(n_total_steps):,} "
                            f"({pct:.1f}%)  "
                            f"left={steps_left:,}  "
                            f"sec/step={sec_per_step:.3f}  "
                            f"eta_epoch={_format_seconds(eta_epoch)}  "
                            f"elapsed={_format_seconds(elapsed_run)}"
                        )
                    else:
                        print(
                            f"[progress epoch {epoch:03d}] "
                            f"step={step_idx:,}  sec/step={sec_per_step:.3f}"
                        )

                # ------------------------------------------------------
                # Step checkpoint: save on optimizer-update boundary only.
                # If save_every_steps is not divisible by grad_accum_steps,
                # save at the first optimizer boundary after the requested step.
                # ------------------------------------------------------
                should_save_step = (
                    runtime.is_main_process
                    and cfg.train.save_every_steps is not None
                    and cfg.train.save_every_steps > 0
                    and next_step_checkpoint is not None
                    and global_step >= next_step_checkpoint
                    and step_out.grad_norm is not None
                )
                if should_save_step:
                    ckpt_path = save_step_checkpoint(
                        out_dir=out_dir,
                        epoch=epoch,
                        step_in_epoch=step_idx,
                        global_step=global_step,
                        cfg=cfg,
                        trainer=trainer,
                        history=history,
                        best_metric=best_metric,
                        best_epoch=best_epoch,
                        keep_last=cfg.train.keep_last_step_checkpoints,
                    )
                    print(
                        f"[checkpoint] saved {ckpt_path} "
                        f"epoch={epoch} step={step_idx:,} global_step={global_step:,}"
                    )
                    while next_step_checkpoint is not None and global_step >= next_step_checkpoint:
                        next_step_checkpoint += int(cfg.train.save_every_steps)

            return _callback

        for epoch in range(start_epoch, cfg.train.epochs + 1):
            epoch_start_time = time.perf_counter()
            stopped_epoch = epoch
            use_phase2_loader = bool(
                cfg.integrated_phase_curriculum.enabled
                and epoch >= int(cfg.integrated_phase_curriculum.phase2_start_epoch)
                and phase2_train_loader is not None
            )
            active_train_loader = (
                phase2_train_loader if use_phase2_loader else train_loader
            )
            active_train_sampler = (
                phase2_train_sampler if use_phase2_loader else train_sampler
            )
            if active_train_sampler is not None:
                active_train_sampler.set_epoch(epoch)

            # Install the rank epoch on every DDP process before any rank-zero
            # preview logging.  A fail-closed structural freeze audit must
            # raise synchronously; raising only on rank zero would strand the
            # other workers at their first collective.
            rank_system = unwrap_system(trainer.system)
            rank_decoder = getattr(rank_system, "decoder", None)
            rank_gate = getattr(rank_decoder, "pathology_rank_gate", None)
            if rank_gate is not None:
                rank_gate.set_epoch(int(epoch))

            if runtime.is_main_process and bool(cfg.integrated_phase_curriculum.enabled):
                preview = trainer.integrated_phase_controller.resolve(epoch)
                print(
                    "[integrated-phase] "
                    f"epoch={epoch} stage={preview['stage']} "
                    f"batch={preview['batch_size']} accum={preview['grad_accum_steps']} "
                    f"lr={preview['learning_rate']:.3g} "
                    f"pathology_scale={preview['pathology_scale']:.3f} "
                    f"precision_forward={preview['precision_forward_enabled']}",
                    flush=True,
                )
                if rank_gate is not None:
                    rank_diag = rank_gate.diagnostics()
                    print(
                        "[pathology-rank] "
                        f"epoch={epoch} point=start "
                        f"capacity={rank_diag['capacity']} "
                        f"expected={rank_diag['expected_rank']:.6f} "
                        f"hard={rank_diag['hard_rank']} "
                        f"mode={rank_diag['mode']} "
                        f"finalized={rank_diag['finalized']} "
                        f"temperature={rank_diag['temperature']:.6f} "
                        f"uncertain={rank_diag['uncertain_fraction']:.6f} "
                        f"freeze_ready={rank_diag['freeze_ready']} "
                        f"freeze_spread={rank_diag['freeze_count_spread']} "
                        f"freeze_near={rank_diag['freeze_near_threshold_uncertain_fraction']:.6f}",
                        flush=True,
                    )

            _boot_log(f"train_epoch:start epoch={epoch}")
            train_out_local = trainer.train_epoch(
                active_train_loader,
                log_every=(cfg.train.log_every if runtime.is_main_process else None),
                step_callback=make_step_callback(epoch, epoch_start_time),
                epoch_index=epoch,
                max_steps=cfg.train.max_train_steps_per_epoch,
            )
            if module_rescue_updater is not None:
                _boot_log(
                    f"prism_module_rescue:start epoch={epoch}"
                )
                rescue_system = unwrap_system(trainer.system)
                joint_rescue_gate = bool(
                    getattr(
                        rescue_system,
                        "generator_count_joint_training",
                        False,
                    )
                )
                if joint_rescue_gate:
                    rescue_system.generator_count_force_deterministic_straight_through = True
                try:
                    rescue_metrics = module_rescue_updater.run_epoch(
                        trainer,
                        epoch=epoch,
                        rank=runtime.rank,
                        world_size=runtime.world_size,
                    )
                    trainer._integrated_last_module_rescue_metrics = dict(
                        rescue_metrics
                    )
                finally:
                    if joint_rescue_gate:
                        rescue_system.generator_count_force_deterministic_straight_through = False
                train_out_local.metrics.update(rescue_metrics)
                _boot_log(
                    f"prism_module_rescue:done epoch={epoch}"
                )
            _boot_log(f"train_epoch:done epoch={epoch}")
            train_out = sync_epoch_output(train_out_local, runtime)
            if runtime.is_main_process and bool(
                cfg.integrated_phase_curriculum.enabled
            ):
                rank_report = pathology_correction_rank_diagnostics(
                    trainer.system
                )
                rank_report.update(
                    {
                        "schema_version": "kmlee_bam.pathology_correction_rank.v2",
                        "epoch": int(epoch),
                        "optimization_effect": (
                            "learned_gate_and_cardinality"
                            if rank_report.get("rank_selection") == "learned"
                            else "none_diagnostic_only"
                        ),
                    }
                )
                with (out_dir / "pathology_correction_rank.jsonl").open(
                    "a", encoding="utf-8"
                ) as handle:
                    handle.write(json.dumps(rank_report, sort_keys=True))
                    handle.write("\n")
                participation_values = [
                    float(row["participation_rank"])
                    for row in rank_report["axes"]
                ]
                train_out.metrics[
                    "metric/pathology_correction_participation_mean"
                ] = float(np.mean(participation_values))
                train_out.metrics[
                    "metric/pathology_correction_rank95_max"
                ] = float(
                    max(
                        int(row["rank_95_energy"])
                        for row in rank_report["axes"]
                    )
                )
                learned_rank_report = rank_report.get("learned_rank")
                if isinstance(learned_rank_report, dict):
                    train_out.metrics[
                        "metric/pathology_rank_expected_epoch"
                    ] = float(learned_rank_report["expected_rank"])
                    train_out.metrics[
                        "metric/pathology_rank_hard_epoch"
                    ] = float(learned_rank_report["hard_rank"])
                    print(
                        "[pathology-rank] "
                        f"epoch={epoch} point=end "
                        f"capacity={learned_rank_report['capacity']} "
                        f"expected={learned_rank_report['expected_rank']:.6f} "
                        f"hard={learned_rank_report['hard_rank']} "
                        f"mode={learned_rank_report['mode']} "
                        f"finalized={learned_rank_report['finalized']} "
                        f"temperature={learned_rank_report['temperature']:.6f} "
                        f"uncertain={learned_rank_report['uncertain_fraction']:.6f} "
                        f"freeze_ready={learned_rank_report['freeze_ready']} "
                        f"freeze_spread={learned_rank_report['freeze_count_spread']} "
                        f"freeze_near={learned_rank_report['freeze_near_threshold_uncertain_fraction']:.6f}",
                        flush=True,
                    )
            if (
                runtime.is_main_process
                and bool(cfg.learned_generator_count.enabled)
                and str(cfg.learned_generator_count.mode) == "joint"
            ):
                joint_system = unwrap_system(trainer.system)
                gate = joint_system.generator_count_gate
                start = int(cfg.learned_generator_count.start_epoch)
                end = max(int(cfg.learned_generator_count.base_end_epoch), start)
                progress = (
                    0.0
                    if epoch <= start
                    else min(1.0, (epoch - start) / max(end - start, 1))
                )
                temperature = float(gate.temperature(progress))
                deterministic_mask = gate.deterministic_mask(
                    temperature=temperature
                ).bool()
                keep_probability = gate.keep_probability(
                    temperature=temperature
                ).detach().float()
                selected_ids = torch.where(deterministic_mask)[0].cpu().tolist()
                registry_indices = tuple(
                    int(value)
                    for value in getattr(
                        joint_system.decoder,
                        "selected_registry_module_indices",
                        tuple(range(int(gate.num_generators))),
                    )
                )
                registry_names = tuple(
                    str(value)
                    for value in getattr(
                        joint_system.decoder,
                        "selected_registry_module_names",
                        tuple(str(value) for value in registry_indices),
                    )
                )
                mask_sha256 = hashlib.sha256(
                    deterministic_mask.detach().cpu().numpy().tobytes()
                ).hexdigest()
                joint_payload = {
                    "schema_version": "kmlee_bam.joint_generator_count.v1",
                    "epoch": int(epoch),
                    "candidate_count": int(gate.num_generators),
                    "hard_active": int(deterministic_mask.sum().item()),
                    "expected_active": float(
                        gate.expected_active_count(
                            temperature=temperature
                        ).detach().cpu()
                    ),
                    "raw_expected_active": float(
                        keep_probability.sum().item()
                    ),
                    "temperature": temperature,
                    "selected_local_generator_ids": selected_ids,
                    "selected_registry_indices": [
                        registry_indices[index] for index in selected_ids
                    ],
                    "selected_registry_module_names": [
                        registry_names[index] for index in selected_ids
                    ],
                    "protected_local_generator_ids": torch.where(
                        gate.protected_mask
                    )[0].cpu().tolist(),
                    "mask_sha256": mask_sha256,
                    "selection_method": "online_joint_hard_concrete",
                    "candidate_k_grid_used": False,
                    "official_test_used_for_count": False,
                }
                for joint_path in (
                    out_dir / "joint_generator_count_latest.json",
                    out_dir / f"joint_generator_count_epoch_{epoch:03d}.json",
                ):
                    with open(joint_path, "w", encoding="utf-8") as handle:
                        json.dump(
                            joint_payload, handle, indent=2, sort_keys=True
                        )
                        handle.write("\n")
                if os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower() == "prism_integrated":
                    print(
                        "═" * 80
                        + f"\n  [epoch {epoch:03d} · Generator 확정] "
                        f"K={joint_payload['hard_active']}/{gate.num_generators} · "
                        f"확률상 E[K]={joint_payload['expected_active']:.1f} · 온도={temperature:.3f}"
                        + "\n     복원 안전제약 full/isolated .... "
                        f"{train_out.metrics.get('metric/generator_violation_full_nll', float('nan')):+.3f}/"
                        f"{train_out.metrics.get('metric/generator_violation_isolated_nll', float('nan')):+.3f} "
                        "(둘 다 0 이하이면 통과)"
                        + f"\n     mask SHA256 .................... {mask_sha256}"
                        + "\n" + "═" * 80,
                        flush=True,
                    )
                else:
                    print(
                        f"[joint-generator epoch {epoch:03d}] "
                        f"K={joint_payload['hard_active']}/{gate.num_generators} "
                        f"E[K]={joint_payload['expected_active']:.1f} "
                        f"tau={temperature:.3f} "
                        f"train_full_violation="
                        f"{train_out.metrics.get('metric/generator_violation_full_nll', float('nan')):.3f} "
                        f"train_isolated_violation="
                        f"{train_out.metrics.get('metric/generator_violation_isolated_nll', float('nan')):.3f} "
                        f"mask_sha256={mask_sha256}",
                        flush=True,
                    )
            if runtime.is_main_process and module_rescue_updater is not None:
                if os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower() == "prism_integrated":
                    print(
                        "═" * 80
                        + f"\n  [epoch {epoch:03d} · 모듈 복원] "
                        f"ramp={train_out.metrics.get('weight/prism_module_rescue_ramp', float('nan')):.2f} · "
                        f"donor그룹={train_out.metrics.get('metric/prism_module_rescue_groups', float('nan')):.1f}"
                        + "\n     branch/full CCC ................. "
                        f"{train_out.metrics.get('metric/prism_module_rescue_branch_flattened_ccc', float('nan')):.3f}/"
                        f"{train_out.metrics.get('metric/prism_module_rescue_full_flattened_ccc', float('nan')):.3f}"
                        + "\n     branch/full Pearson ............. "
                        f"{train_out.metrics.get('metric/prism_module_rescue_branch_pearson', float('nan')):.3f}/"
                        f"{train_out.metrics.get('metric/prism_module_rescue_full_pearson', float('nan')):.3f}"
                        + "\n     weighted loss/grad .............. "
                        f"{train_out.metrics.get('loss/prism_module_rescue', float('nan')):.5f}/"
                        f"{train_out.metrics.get('metric/prism_module_rescue_grad_norm', float('nan')):.3f}"
                        + "\n     celltype draw min-max ........... "
                        f"{train_out.metrics.get('metric/prism_module_rescue_celltype_draw_min', float('nan')):.0f}-"
                        f"{train_out.metrics.get('metric/prism_module_rescue_celltype_draw_max', float('nan')):.0f}"
                        + "\n" + "═" * 80,
                        flush=True,
                    )
                else:
                    print(
                        f"[module-rescue epoch {epoch:03d}] "
                        f"loss="
                        f"{train_out.metrics.get('loss/prism_module_rescue', float('nan')):.5f} "
                        f"branch="
                        f"{train_out.metrics.get('loss/prism_module_rescue_branch', float('nan')):.5f} "
                        f"full="
                        f"{train_out.metrics.get('loss/prism_module_rescue_full', float('nan')):.5f} "
                        f"branch_flat_ccc="
                        f"{train_out.metrics.get('metric/prism_module_rescue_branch_flattened_ccc', float('nan')):.3f} "
                        f"full_flat_ccc="
                        f"{train_out.metrics.get('metric/prism_module_rescue_full_flattened_ccc', float('nan')):.3f} "
                        f"branch_r="
                        f"{train_out.metrics.get('metric/prism_module_rescue_branch_pearson', float('nan')):.3f} "
                        f"full_r="
                        f"{train_out.metrics.get('metric/prism_module_rescue_full_pearson', float('nan')):.3f} "
                        f"groups="
                        f"{train_out.metrics.get('metric/prism_module_rescue_groups', float('nan')):.1f} "
                        f"grad="
                        f"{train_out.metrics.get('metric/prism_module_rescue_grad_norm', float('nan')):.3f} "
                        f"ramp="
                        f"{train_out.metrics.get('weight/prism_module_rescue_ramp', float('nan')):.2f} "
                        f"ct_draws="
                        f"{train_out.metrics.get('metric/prism_module_rescue_celltype_draw_min', float('nan')):.0f}-"
                        f"{train_out.metrics.get('metric/prism_module_rescue_celltype_draw_max', float('nan')):.0f}",
                        flush=True,
                    )
            # --- donor x region gradient-SNR audit (DIAGNOSTIC; opt-in via DONOR_SNR_AUDIT=1; fully guarded) ---
            # Default OFF: when the env var is unset this whole block is a no-op, so existing runs are unaffected.
            # rank0 builds (once) + runs the audit on the UNWRAPPED model (no all-reduce); a barrier keeps ranks in sync.
            if os.environ.get("DONOR_SNR_AUDIT", "0") in ("1", "true", "True"):
                _snr_every = max(1, int(os.environ.get("DONOR_SNR_EVERY", "1") or "1"))
                if epoch % _snr_every == 0:
                    try:
                        if runtime.is_main_process and not getattr(trainer, "_donor_snr_off", False):
                            if not hasattr(trainer, "_donor_snr_audit"):
                                import numpy as _np
                                from kmlee_bam.training.donor_snr_audit import DonorSNRAudit as _DSA, make_blocks as _mb
                                from kmlee_bam.preprocessing.data.make_ordinal_data_zarr import read_zarr_dataframe as _rzd
                                _ds = train_loader.dataset
                                _rk = os.environ.get("DONOR_SNR_REGION_KEY", "brain_region")
                                _reg = _rzd(_ds.root, "obs")[_rk].astype(str).values[_np.asarray(_ds.row_idx)]
                                trainer._donor_snr_blocks = _mb(
                                    _ds, _reg,
                                    n_blocks=int(os.environ.get("DONOR_SNR_BLOCKS", "24")),
                                    min_cells=int(os.environ.get("DONOR_SNR_MINCELLS", "64")),
                                    max_cells=int(os.environ.get("DONOR_SNR_MAXCELLS", "256")),
                                )
                                trainer._donor_snr_audit = _DSA()
                                _by_reg = {}
                                for _blk in trainer._donor_snr_blocks:
                                    _by_reg[_blk["region"]] = _by_reg.get(_blk["region"], 0) + 1
                                print(
                                    "[donor-SNR] built "
                                    f"{len(trainer._donor_snr_blocks)} donor x region blocks "
                                    f"{_by_reg}",
                                    flush=True,
                                )
                            _dev = next(trainer.system.parameters()).device
                            trainer._donor_snr_audit.run(
                                trainer.system, trainer._donor_snr_blocks, device=_dev, epoch=epoch
                            )
                    except Exception as _snr_err:  # never let a diagnostic crash training
                        trainer._donor_snr_off = True
                        print(f"[donor-SNR] disabled after error (epoch {epoch}): {_snr_err}", flush=True)
                    finally:
                        try:
                            import torch.distributed as _dist
                            if getattr(runtime, "distributed", False) and _dist.is_available() and _dist.is_initialized():
                                _dist.barrier()
                        except Exception:
                            pass
            generator_search_record: Optional[Dict[str, Any]] = None
            search_due = bool(
                cfg.generator_budget_search.enabled
                and epoch >= int(cfg.generator_budget_search.start_epoch)
                and (
                    epoch - int(cfg.generator_budget_search.start_epoch)
                )
                % int(cfg.generator_budget_search.interval_epochs)
                == 0
            )
            if search_due:
                if runtime.is_main_process:
                    if generator_search is None:
                        raise RuntimeError(
                            "rank 0 generator search driver was not initialized."
                        )
                    print(
                        f"[generator-budget epoch {epoch:03d}] "
                        f"checking active={generator_search.active_count} "
                        "with fixed bypass outputs...",
                        flush=True,
                    )
                    generator_search_record = generator_search.run_check(
                        epoch=epoch
                    )
                    if generator_search_record.get("due", False):
                        decision = generator_search_record["decision"]
                        incumbent = generator_search_record[
                            "incumbent_metrics"
                        ]
                        candidate = generator_search_record[
                            "candidate_metrics"
                        ]
                        print(
                            f"[generator-budget epoch {epoch:03d}] "
                            f"panel={generator_search_record['panel_index']} "
                            f"proposal={len(generator_search_record['proposal'])} "
                            f"eligible={decision['eligible']} "
                            f"applied={decision['applied']} "
                            f"confirm={decision['confirmation_count']}/"
                            f"{decision['required_confirmations']} "
                            f"active={generator_search_record['active_count']} "
                            f"module={incumbent['module_recovery']:.4f}"
                            f"->{candidate['module_recovery']:.4f} "
                            f"isolated_module="
                            f"{incumbent['isolated_module_recovery']:.4f}"
                            f"->{candidate['isolated_module_recovery']:.4f} "
                            f"nll={incumbent['full_nll']:.4f}"
                            f"->{candidate['full_nll']:.4f} "
                            f"isolated_nll="
                            f"{incumbent['isolated_nll']:.4f}"
                            f"->{candidate['isolated_nll']:.4f} "
                            f"reason={','.join(decision['reason_codes'])}",
                            flush=True,
                        )

                if runtime.distributed:
                    decoder = unwrap_system(trainer.system).decoder
                    # Rank 0 has installed the controller mask.  Explicitly
                    # broadcast because DDP was constructed with
                    # broadcast_buffers=False.
                    dist.broadcast(
                        decoder.generator_active_mask,
                        src=0,
                    )
                    barrier(runtime)

            record: Dict[str, Dict[str, float]] = {"train": train_out.metrics}
            if (
                runtime.is_main_process
                and generator_search_record is not None
                and generator_search_record.get("due", False)
            ):
                decision = generator_search_record["decision"]
                record["generator_search"] = {
                    "active_count": float(
                        generator_search_record["active_count"]
                    ),
                    "proposal_count": float(
                        len(generator_search_record["proposal"])
                    ),
                    "eligible": float(bool(decision["eligible"])),
                    "applied": float(bool(decision["applied"])),
                    "module_recovery": float(
                        generator_search_record["incumbent_metrics"][
                            "module_recovery"
                        ]
                    ),
                    "candidate_module_recovery": float(
                        generator_search_record["candidate_metrics"][
                            "module_recovery"
                        ]
                    ),
                }

            metric_source = train_out.metrics
            monitored_split = "train"
            val_time_sec = 0.0
            has_validation = (
                val_ds is not None
                and cfg.train.eval_every is not None
                and int(cfg.train.eval_every) > 0
                and epoch % int(cfg.train.eval_every) == 0
            )
            if (
                bool(cfg.train.require_validation_for_early_stopping)
                and not has_validation
            ):
                raise RuntimeError(
                    "validation-only training control reached an epoch without "
                    "validation; refusing train-metric fallback"
                )

            if has_validation:
                val_eval_start = time.perf_counter()

                if runtime.distributed and bool(cfg.train.use_rank0_eval):
                    # Optional legacy path:
                    # rank 0 evaluates, other ranks wait for broadcast.
                    val_out = evaluate_epoch_rank0_only(
                        trainer,
                        rank0_val_loader,
                        runtime,
                        log_every=(cfg.train.eval_log_every if runtime.is_main_process else None),
                        max_steps=cfg.train.max_eval_steps,
                    )
                    if val_out is None:
                        raise RuntimeError(
                            "Validation dataset exists, but rank0-only validation returned None."
                        )

                else:
                    # Preferred path:
                    # every rank evaluates its own DistributedSampler shard,
                    # then metrics are reduced across ranks.
                    if val_loader is None:
                        raise RuntimeError(
                            "Validation dataset exists, but val_loader is None. "
                            "Check build_loaders()."
                        )

                    if _val_sampler is not None:
                        _val_sampler.set_epoch(epoch)

                    val_out_local = trainer.evaluate_epoch(
                        val_loader,
                        log_every=(cfg.train.eval_log_every if runtime.is_main_process else None),
                        max_steps=cfg.train.max_eval_steps,
                    )
                    val_out = sync_epoch_output(val_out_local, runtime)

                val_time_sec = time.perf_counter() - val_eval_start
                record["val"] = val_out.metrics
                metric_source = val_out.metrics
                monitored_split = "val"
            else:
                val_time_sec = 0.0

            # Reversible hard-mask safety audit.  The gate is allowed to keep
            # a new state only after consecutive train/validation safety
            # evidence; a later violation restores the last confirmed logits.
            if bool(
                cfg.learned_generator_count.enabled
                and str(cfg.learned_generator_count.mode) == "joint"
            ):
                soft_end = int(
                    getattr(cfg.learned_generator_count, "soft_end_epoch", 0)
                    or 0
                )
                hard_start = max(
                    int(cfg.learned_generator_count.start_epoch),
                    soft_end + 1 if soft_end > 0 else int(cfg.learned_generator_count.start_epoch),
                )
                if epoch >= hard_start:
                    full_violation = float(
                        metric_source.get(
                            "metric/generator_violation_full_nll", math.inf
                        )
                    )
                    isolated_violation = float(
                        metric_source.get(
                            "metric/generator_violation_isolated_nll", math.inf
                        )
                    )
                    # Relative BAM/PHU weights are learned from the sealed
                    # training bank.  Validation donors intentionally have no
                    # bank entries, so their default ESS=1 cannot certify a
                    # hard architecture.  Gate safety must use the train ESS.
                    phu_ess = float(
                        train_out.metrics.get(
                            "metric/phu_rel_ess_fraction", 1.0
                        )
                    )
                    rescue_ccc = float(
                        train_out.metrics.get(
                            "metric/prism_module_rescue_full_flattened_ccc",
                            math.nan,
                        )
                    )
                    feasible = bool(
                        math.isfinite(full_violation)
                        and math.isfinite(isolated_violation)
                        and full_violation <= 0.0
                        and isolated_violation <= 0.0
                        and phu_ess >= 0.90
                        and math.isfinite(rescue_ccc)
                        and rescue_ccc > -0.05
                    )
                    gate = unwrap_system(trainer.system).generator_count_gate
                    audit = gate.apply_safety_audit(feasible=feasible)
                    record["generator_safety_audit"] = {
                        key: float(value) if isinstance(value, (bool, int)) else value
                        for key, value in audit.items()
                    }
                    if runtime.is_main_process:
                        print(
                            "[generator-safe-mask] "
                            f"epoch={epoch} feasible={feasible} "
                            f"streak={audit['streak']}/"
                            f"{cfg.learned_generator_count.safe_mask_confirmations} "
                            f"confirmed={audit['has_confirmed_safe_mask']} "
                            f"rollback={audit['rolled_back']}",
                            flush=True,
                        )

            if runtime.is_main_process and rank_gate is not None:
                # One human-readable dashboard plus one machine-readable JSONL
                # row per epoch.  This is intentionally after validation and
                # the generator safety audit, so a tmux tail and downstream
                # monitor see the same finalized evidence.
                capacity_metrics = dict(train_out.metrics)
                capacity_metrics.update(metric_source)
                capacity_system = unwrap_system(trainer.system)
                capacity_generator_gate = getattr(
                    capacity_system, "generator_count_gate", None
                )
                capacity_generator_config = getattr(
                    capacity_system, "generator_count_config", None
                )
                capacity_pathology_scale = float(
                    getattr(
                        getattr(capacity_system, "decoder", None),
                        "pathology_curriculum_scale",
                        torch.tensor(1.0),
                    ).detach().float().cpu()
                )
                capacity_precision_config = getattr(
                    getattr(capacity_system, "precision_head", None),
                    "cfg",
                    None,
                )
                capacity_record = build_architecture_capacity_record(
                    epoch=epoch,
                    total_epochs=cfg.train.epochs,
                    pathology_scale=capacity_pathology_scale,
                    rank_gate=rank_gate,
                    generator_gate=capacity_generator_gate,
                    generator_config=capacity_generator_config,
                    metrics=capacity_metrics,
                    paper_config=capacity_precision_config,
                )
                if bool(
                    capacity_record["coordination_guard"][
                        "simultaneous_cardinality"
                    ]
                ):
                    raise RuntimeError(
                        "pathology-rank and generator cardinality objectives "
                        "became active in the same epoch"
                    )
                print(
                    format_architecture_capacity_record(capacity_record),
                    flush=True,
                )
                with (out_dir / "architecture_capacity_history.jsonl").open(
                    "a", encoding="utf-8"
                ) as handle:
                    handle.write(json.dumps(capacity_record, sort_keys=True))
                    handle.write("\n")
                for capacity_path in (
                    out_dir / "architecture_capacity_latest.json",
                    out_dir / f"architecture_capacity_epoch_{epoch:03d}.json",
                ):
                    with capacity_path.open("w", encoding="utf-8") as handle:
                        json.dump(
                            capacity_record,
                            handle,
                            indent=2,
                            sort_keys=True,
                        )
                        handle.write("\n")

            if runtime.is_main_process:
                epoch_time_sec = time.perf_counter() - epoch_start_time
                elapsed_total_sec = time.perf_counter() - run_start_time
                avg_epoch_sec = elapsed_total_sec / max(epoch, 1)
                eta_sec = avg_epoch_sec * max(0, cfg.train.epochs - epoch)
                compact_experiment_log = (
                    os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower()
                    in {"prism_experiment", "prism_informative"}
                )
                if compact_experiment_log:
                    val_metrics = record.get("val", {})
                    print(
                        f"\n[epoch {epoch:03d}/{cfg.train.epochs:03d}] COMPLETE  "
                        f"train={train_out.metrics.get('loss/total', float('nan')):.4f}  "
                        f"val={val_metrics.get('loss/total', float('nan')):.4f}  "
                        f"PRISM_val={metric_source.get('loss/prism_branch_nll', float('nan')):.4f}  "
                        f"time={_format_seconds(epoch_time_sec)}  "
                        f"eta={_format_seconds(eta_sec)}",
                        flush=True,
                    )
                    print(
                        "  validation  "
                        f"raw={metric_source.get('loss/prism_branch_raw_nll', float('nan')):.4f}  "
                        f"nonzero={metric_source.get('loss/prism_branch_nonzero_nll', float('nan')):.4f}  "
                        "rms(common/personal/response)="
                        f"{metric_source.get('metric/prism_common_module_rms', float('nan')):.4f}/"
                        f"{metric_source.get('metric/prism_personal_module_rms', float('nan')):.4f}/"
                        f"{metric_source.get('metric/prism_response_module_rms', float('nan')):.4f}  "
                        "leak(path/age/state)="
                        f"{metric_source.get('loss/prism_pathology_leak', float('nan')):.4f}/"
                        f"{metric_source.get('loss/prism_age_leak', float('nan')):.4f}/"
                        f"{metric_source.get('loss/prism_state_pathology_leak', float('nan')):.4f}",
                        flush=True,
                    )
                    if "metric/agp_effective_tokens" in metric_source:
                        agp_top1 = metric_source.get(
                            "metric/agp_top1_weight", float("nan")
                        )
                        agp_overlap = metric_source.get(
                            "metric/agp_head_overlap", float("nan")
                        )
                        agp_centered_overlap = metric_source.get(
                            "metric/agp_centered_head_overlap", float("nan")
                        )
                        agp_query_rank = metric_source.get(
                            "metric/agp_query_effective_rank", float("nan")
                        )
                        agp_score_rank = metric_source.get(
                            "metric/agp_score_effective_rank", float("nan")
                        )
                        agp_status = (
                            "WARN uniform+duplicate"
                            if agp_top1 < 0.01 and agp_overlap > 0.90
                            else "WARN duplicate-heads"
                            if (
                                agp_overlap > 0.90
                                and (
                                    not math.isfinite(agp_centered_overlap)
                                    or agp_centered_overlap > 0.75
                                )
                            )
                            else "OK selective"
                        )
                        print(
                            "  AGP validation  "
                            f"effective_tokens={metric_source.get('metric/agp_effective_tokens', float('nan')):.1f}  "
                            f"top1={agp_top1:.3f}  "
                            f"top5={metric_source.get('metric/agp_top5_mass', float('nan')):.3f}  "
                            f"overlap(raw/centered)={agp_overlap:.3f}/{agp_centered_overlap:.3f}  "
                            f"rank(query/score)={agp_query_rank:.2f}/{agp_score_rank:.2f}  "
                            f"temp={metric_source.get('metric/agp_temperature', float('nan')):.3f}  "
                            f"mean_gate={metric_source.get('metric/agp_mean_residual_gate', float('nan')):.3f}  "
                            f"aux={metric_source.get('loss/agp_auxiliary', 0.0):.4f}  "
                            f"status={agp_status}",
                            flush=True,
                        )
                    if float(
                        getattr(
                            cfg.precision_medicine,
                            "lambda_support_query_infonce",
                            0.0,
                        )
                    ) > 0.0:
                        train_metrics = train_out.metrics
                        own = train_metrics.get(
                            "metric/prism_support_query_own_nll", float("nan")
                        )
                        negative = train_metrics.get(
                            "metric/prism_support_query_negative_nll", float("nan")
                        )
                        own_delta = float("nan")
                        negative_delta = float("nan")
                        if history:
                            previous_train = history[-1].get("train", {})
                            own_delta = own - previous_train.get(
                                "metric/prism_support_query_own_nll", float("nan")
                            )
                            negative_delta = negative - previous_train.get(
                                "metric/prism_support_query_negative_nll", float("nan")
                            )
                        if not (
                            math.isfinite(own_delta)
                            and math.isfinite(negative_delta)
                        ):
                            infonce_status = "BASELINE"
                        elif negative - own <= 0.0:
                            infonce_status = "WARN no-separation"
                        elif own_delta >= -1.0e-4 and negative_delta > 1.0e-4:
                            infonce_status = "WARN negative-only"
                        elif own_delta < -1.0e-4:
                            infonce_status = "OK own-improves"
                        else:
                            infonce_status = "MIXED"
                        print(
                            "  InfoNCE train  "
                            f"loss={train_metrics.get('loss/prism_support_query_infonce', float('nan')):.4f}  "
                            f"own={own:.4f} (delta={own_delta:+.4f})  "
                            f"matched_other={negative:.4f} (delta={negative_delta:+.4f})  "
                            f"gap={negative - own:+.4f}  "
                            f"valid={train_metrics.get('metric/prism_support_query_valid_fraction', float('nan')):.2f}  "
                            f"anchor={train_metrics.get('loss/prism_support_query_anchor', float('nan')):.4f}  "
                            f"status={infonce_status}",
                            flush=True,
                        )
                    monitored_now = metric_source.get(
                        cfg.train.metric_name, float("nan")
                    )
                    print(
                        "  selection  "
                        f"{monitored_split}:{cfg.train.metric_name}={monitored_now:.6f}  "
                        f"prior_best={_format_optional_float(best_metric)}  "
                        f"bad_epochs_before={bad_epochs}  "
                        f"test={'SEALED' if cfg.train.skip_final_test else 'configured'}",
                        flush=True,
                    )
                elif "val" in record:
                    print(
                        f"[epoch {epoch:03d}] "
                        f"train_total={train_out.metrics.get('loss/total', float('nan')):.6f}  "
                        f"val_total={record['val'].get('loss/total', float('nan')):.6f}  "
                        f"train_time={_format_seconds(train_out_local.elapsed_sec)}  "
                        f"val_time={_format_seconds(val_time_sec)}  "
                        f"epoch_time={_format_seconds(epoch_time_sec)}  "
                        f"elapsed={_format_seconds(elapsed_total_sec)}  "
                        f"eta={_format_seconds(eta_sec)}"
                    )
                else:
                    print(
                        f"[epoch {epoch:03d}] "
                        f"train_total={train_out.metrics.get('loss/total', float('nan')):.6f}  "
                        f"train_time={_format_seconds(train_out_local.elapsed_sec)}  "
                        f"epoch_time={_format_seconds(epoch_time_sec)}  "
                        f"elapsed={_format_seconds(elapsed_total_sec)}  "
                        f"eta={_format_seconds(eta_sec)}"
                    )
                if (
                    not compact_experiment_log
                    and "loss/prism_branch_nll" in metric_source
                ):
                    print(
                        f"[PRISM epoch {epoch:03d}] split={monitored_split} "
                        f"masked_objective={metric_source.get('loss/prism_branch_nll', float('nan')):.6f} "
                        f"raw_nll={metric_source.get('loss/prism_branch_raw_nll', float('nan')):.6f} "
                        f"nonzero_nll={metric_source.get('loss/prism_branch_nonzero_nll', float('nan')):.6f} "
                        f"normal_region_rms={metric_source.get('metric/prism_normal_region_module_rms', float('nan')):.5f} "
                        f"age_rms={metric_source.get('metric/prism_age_module_rms', float('nan')):.5f} "
                        f"common_rms={metric_source.get('metric/prism_common_module_rms', float('nan')):.5f} "
                        f"personal_rms={metric_source.get('metric/prism_personal_module_rms', float('nan')):.5f} "
                        f"response_rms={metric_source.get('metric/prism_response_module_rms', float('nan')):.5f} "
                        f"code_rms={metric_source.get('metric/prism_code_rms', float('nan')):.4f} "
                        f"support={metric_source.get('metric/prism_support_contexts', float('nan')):.1f} "
                        "ADNC_INPUT=false",
                        flush=True,
                    )

            monitored = float(metric_source.get(cfg.train.metric_name, math.nan))

            # All ranks reach this boundary.  For a committed exact mask this
            # is the once-per-epoch collective agreement audit; for every
            # other training mode it is a no-op.
            run_system_integrity_checks(
                trainer.system,
                stage=f"epoch_{epoch:03d}_complete",
                distributed=runtime.distributed,
            )

            if runtime.is_main_process:
                history.append(record)

                # Always keep a last checkpoint at epoch end.
                save_checkpoint(
                    out_dir / "checkpoint_last.pt",
                    epoch=epoch,
                    global_step=global_step,
                    step_in_epoch=None,
                    cfg=cfg,
                    system=trainer.system,
                    criterion=trainer.criterion,
                    optimizer=trainer.optimizer,
                    scheduler=trainer.scheduler,
                    trainer=trainer,
                    history=history,
                    best_metric=best_metric,
                    best_epoch=best_epoch,
                    stop_reason=None,
                )

                if cfg.train.save_every > 0 and epoch % cfg.train.save_every == 0:
                    save_checkpoint(
                        out_dir / f"checkpoint_epoch_{epoch:03d}.pt",
                        epoch=epoch,
                        global_step=global_step,
                        step_in_epoch=None,
                        cfg=cfg,
                        system=trainer.system,
                        criterion=trainer.criterion,
                        optimizer=trainer.optimizer,
                        scheduler=trainer.scheduler,
                        trainer=trainer,
                        history=history,
                        best_metric=best_metric,
                        best_epoch=best_epoch,
                        stop_reason=None,
                    )

                # Multi-criteria improvement check. For each registered
                # criterion we test against its own rolling best with its own
                # min_delta. The primary determines `save_best` (single
                # canonical checkpoint), while bad_epochs only increments when
                # NO criterion improves (OR mode).
                primary_name = es_criteria[0]["name"]
                joint_selection_feasible, joint_selection_violations = (
                    joint_generator_selection_feasible(cfg, metric_source)
                )
                per_criterion_improved: dict[str, bool] = {}
                per_criterion_values: dict[str, float] = {}
                for crit in es_criteria:
                    name = crit["name"]
                    val = float(metric_source.get(name, math.nan))
                    per_criterion_values[name] = val
                    if math.isnan(val) or not joint_selection_feasible:
                        per_criterion_improved[name] = False
                        continue
                    improved_c = metric_better(
                        val,
                        es_best[name],
                        crit["mode"],
                        min_delta=crit["min_delta"],
                    )
                    per_criterion_improved[name] = improved_c
                    if improved_c:
                        es_best[name] = val

                primary_improved = per_criterion_improved.get(primary_name, False)
                any_improved = any(per_criterion_improved.values())

                # Track primary's best value for reporting/checkpoint metadata.
                if primary_improved:
                    best_metric = per_criterion_values[primary_name]
                    best_epoch = epoch

                # --- `checkpoint_best.pt` decision ------------------------
                # save_best_mode = "composite" uses a weighted z-score sum
                # across criteria, but only after the warm-up window
                # (early_stopping_min_epochs). Before warm-up, or in
                # "primary" mode, fall back to primary-only improvement.
                composite_score = math.nan
                composite_z: Dict[str, float] = {}
                composite_contrib: Dict[str, float] = {}
                composite_saved = False
                # Always-defined so the early-stopping decision below can read
                # these even in warm-up / legacy "primary" mode.
                use_composite = False
                improved_comp = False
                if str(cfg.train.save_best_mode).lower() == "composite":
                    composite_score, composite_z, composite_contrib = compute_composite_score(
                        criteria=es_criteria,
                        values=per_criterion_values,
                        stats=composite_stats,
                        eps=float(cfg.train.composite_std_eps),
                    )
                    use_composite = (
                        epoch >= cfg.train.early_stopping_min_epochs
                        and not math.isnan(composite_score)
                    )
                    if use_composite:
                        improved_comp = (
                            joint_selection_feasible
                            and (
                                best_composite is None
                                or composite_score
                                > best_composite
                                + float(cfg.train.early_stopping_min_delta)
                            )
                        )
                        if improved_comp and cfg.train.save_best:
                            best_composite = composite_score
                            best_composite_epoch = epoch
                            save_checkpoint(
                                out_dir / "checkpoint_best.pt",
                                epoch=epoch,
                                global_step=global_step,
                                step_in_epoch=None,
                                cfg=cfg,
                                system=trainer.system,
                                criterion=trainer.criterion,
                                optimizer=trainer.optimizer,
                                scheduler=trainer.scheduler,
                                trainer=trainer,
                                history=history,
                                best_metric=composite_score,
                                best_epoch=epoch,
                                stop_reason=None,
                            )
                            composite_saved = True
                    else:
                        # Warm-up: keep primary-driven checkpoint_best to
                        # preserve sensible state if learning is cut short.
                        if primary_improved and cfg.train.save_best:
                            save_checkpoint(
                                out_dir / "checkpoint_best.pt",
                                epoch=epoch,
                                global_step=global_step,
                                step_in_epoch=None,
                                cfg=cfg,
                                system=trainer.system,
                                criterion=trainer.criterion,
                                optimizer=trainer.optimizer,
                                scheduler=trainer.scheduler,
                                trainer=trainer,
                                history=history,
                                best_metric=best_metric,
                                best_epoch=best_epoch,
                                stop_reason=None,
                            )
                else:
                    # Legacy "primary" mode.
                    if primary_improved and cfg.train.save_best:
                        save_checkpoint(
                            out_dir / "checkpoint_best.pt",
                            epoch=epoch,
                            global_step=global_step,
                            step_in_epoch=None,
                            cfg=cfg,
                            system=trainer.system,
                            criterion=trainer.criterion,
                            optimizer=trainer.optimizer,
                            scheduler=trainer.scheduler,
                            trainer=trainer,
                            history=history,
                            best_metric=best_metric,
                            best_epoch=best_epoch,
                            stop_reason=None,
                        )

                # --- Per-criterion best checkpoints (independent) --------
                saved_per_crit: list[str] = []
                if cfg.train.save_best_per_criterion and cfg.train.save_best:
                    for crit in es_criteria:
                        name = crit["name"]
                        if not per_criterion_improved.get(name):
                            continue
                        short = crit["short_name"]
                        save_checkpoint(
                            out_dir / f"checkpoint_best_{short}.pt",
                            epoch=epoch,
                            global_step=global_step,
                            step_in_epoch=None,
                            cfg=cfg,
                            system=trainer.system,
                            criterion=trainer.criterion,
                            optimizer=trainer.optimizer,
                            scheduler=trainer.scheduler,
                            trainer=trainer,
                            history=history,
                            best_metric=es_best[name],
                            best_epoch=epoch,
                            stop_reason=None,
                        )
                        saved_per_crit.append(short)

                # --- Update composite running stats AFTER decisions -------
                # This guarantees the z-score at epoch t uses [0..t-1] only.
                if joint_selection_feasible:
                    composite_stats.update(per_criterion_values)

                # Early-stopping criterion must match best-selection: when the
                # composite drives checkpoint_best, count patience on composite
                # improvement (not any_improved) so training stops at the
                # composite peak instead of chasing single-metric noise past it.
                # Warm-up / legacy "primary" mode (or save_best off) falls back.
                stop_on_composite = bool(use_composite and cfg.train.save_best)
                stop_improved = improved_comp if stop_on_composite else any_improved
                if stop_improved:
                    bad_epochs = 0
                else:
                    bad_epochs += 1

                # Backwards-compatible primary-only log line (preserved).
                monitored = per_criterion_values.get(primary_name, math.nan)
                print(
                    f"  monitored[{monitored_split}:{primary_name}]={monitored:.6f}  "
                    f"best={float('nan') if best_metric is None else best_metric:.6f}  "
                    f"bad_epochs={bad_epochs}"
                )
                if joint_selection_violations:
                    print(
                        "  joint_generator_selection_feasible="
                        f"{joint_selection_feasible}  "
                        "full_violation="
                        f"{joint_selection_violations['full_nll']:+.6f}  "
                        "isolated_violation="
                        f"{joint_selection_violations['isolated_nll']:+.6f}"
                    )
                # Per-criterion verdict line. Visible only when extras exist.
                if len(es_criteria) > 1:
                    parts = []
                    for crit in es_criteria:
                        name = crit["name"]
                        val = per_criterion_values[name]
                        bst = es_best[name]
                        flag = "↑" if per_criterion_improved.get(name) else "·"
                        parts.append(
                            f"{name}={val:.6f}({flag},best="
                            f"{float('nan') if bst is None else bst:.6f})"
                        )
                    print(
                        f"  early_stopping[any_improves]={any_improved}  "
                        + "  ".join(parts)
                    )

                # Composite breakdown line.
                if str(cfg.train.save_best_mode).lower() == "composite":
                    use_composite = (
                        epoch >= cfg.train.early_stopping_min_epochs
                        and not math.isnan(composite_score)
                    )
                    if use_composite:
                        breakdown = "  ".join(
                            f"{c['short_name']}={composite_z.get(c['name'], 0.0):+.2f}"
                            f"(×{c.get('weight',1.0):.2f}→{composite_contrib.get(c['name'],0.0):+.2f})"
                            for c in es_criteria
                        )
                        best_str = (
                            f"{best_composite:+.3f}"
                            if best_composite is not None else "nan"
                        )
                        epoch_str = (
                            f"@epoch{best_composite_epoch:03d}"
                            if best_composite_epoch is not None else ""
                        )
                        saved_mark = "  → saved checkpoint_best.pt" if composite_saved else ""
                        print(
                            f"  composite[mode=composite]={composite_score:+.3f}  "
                            f"best={best_str}{epoch_str}  {breakdown}{saved_mark}"
                        )
                    else:
                        # Warm-up: indicate composite is collecting stats.
                        print(
                            f"  composite[mode=composite,warmup]  "
                            f"epoch < min_epochs={cfg.train.early_stopping_min_epochs}; "
                            f"checkpoint_best follows primary"
                        )

                if saved_per_crit:
                    print(
                        f"  per-criterion best saved: "
                        + ", ".join(saved_per_crit)
                    )

                should_stop = bool(
                    cfg.train.early_stopping
                    and epoch >= cfg.train.early_stopping_min_epochs
                    and bad_epochs >= cfg.train.early_stopping_patience
                )
                if should_stop:
                    stop_reason = "early_stopping"
                    if stop_on_composite:
                        es_label = "composite"
                        es_best_epoch = best_composite_epoch
                        es_best_val = best_composite
                    else:
                        es_label = cfg.train.metric_name
                        es_best_epoch = best_epoch
                        es_best_val = best_metric
                    print(
                        f"[early stopping] epoch={epoch:03d}  "
                        f"split={monitored_split}  "
                        f"criterion={es_label}  "
                        f"best_epoch={es_best_epoch}  "
                        f"best={(es_best_val if es_best_val is not None else float('nan')):.6f}  "
                        f"patience={cfg.train.early_stopping_patience}"
                    )
            else:
                should_stop = False

            if runtime.distributed:
                shared_state = [
                    {
                        "best_metric": best_metric,
                        "best_epoch": best_epoch,
                        "bad_epochs": bad_epochs,
                        "stop_reason": stop_reason,
                        "should_stop": should_stop,
                    }
                ]
                broadcast_object_list_control(runtime, shared_state, src=0)
                state = shared_state[0] or {}
                best_metric = state.get("best_metric", best_metric)
                best_epoch = state.get("best_epoch", best_epoch)
                bad_epochs = state.get("bad_epochs", bad_epochs)
                stop_reason = state.get("stop_reason", stop_reason)
                should_stop = state.get("should_stop", should_stop)

            if should_stop:
                break

        control_barrier(runtime)

        if runtime.is_main_process:
            with open(out_dir / "history.json", "w", encoding="utf-8") as f:
                json.dump(history, f, indent=2, ensure_ascii=False)

            train_eval_loader = build_plain_eval_loader(cfg, train_ds, shuffle=False)
            val_eval_loader = (
                build_plain_eval_loader(cfg, val_ds, shuffle=False)
                if val_ds is not None
                else None
            )
            test_eval_loader = None
            if not bool(cfg.train.skip_final_test):
                test_eval_loader = (
                    rank0_test_loader
                    if rank0_test_loader is not None
                    else (
                        build_plain_eval_loader(cfg, test_ds, shuffle=False)
                        if test_ds is not None
                        else None
                    )
                )

            with use_unwrapped_system_for_rank0_only_ops(trainer):
                best_ckpt_path = out_dir / "checkpoint_best.pt"
                use_best_for_eval = bool(
                    cfg.train.reload_best_before_test
                    and cfg.train.save_best
                    and best_ckpt_path.exists()
                )

                if use_best_for_eval:
                    _load_checkpoint_for_analysis(
                        best_ckpt_path,
                        system=trainer.system,
                        criterion=trainer.criterion,
                        device=trainer.device,
                        trainer=trainer,
                    )
                    print("Reloaded best checkpoint for final test evaluation / latent export.")

                if test_eval_loader is not None:
                    test_out = trainer.evaluate_epoch(
                        test_eval_loader,
                        log_every=None,
                        max_steps=cfg.train.max_eval_steps,
                    )
                    print("=" * 72)
                    print("final test metrics")
                    for k, v in sorted(test_out.metrics.items()):
                        print(f"{k:24s}: {v:.6f}")
                    print("=" * 72)

                    with open(out_dir / "test_metrics.json", "w", encoding="utf-8") as f:
                        json.dump(test_out.metrics, f, indent=2, ensure_ascii=False)

                if cfg.train.save_split_latents:
                    prefix = "best" if use_best_for_eval else "last"
                    _save_split_latents(
                        trainer=trainer,
                        out_dir=out_dir,
                        train_loader=train_eval_loader,
                        val_loader=val_eval_loader,
                        test_loader=test_eval_loader,
                        prefix=prefix,
                    )
                    print(f"Saved {prefix}-checkpoint split latents: train / val / test")

            selected_mode = (
                "composite"
                if str(cfg.train.save_best_mode).lower() == "composite"
                and best_composite_epoch is not None
                else "primary"
            )
            selected_epoch = (
                best_composite_epoch if selected_mode == "composite" else best_epoch
            )
            selected_score = (
                best_composite if selected_mode == "composite" else best_metric
            )
            _write_training_summary(
                out_dir=out_dir,
                stopped_epoch=stopped_epoch,
                best_epoch=selected_epoch,
                best_metric=selected_score,
                monitored_split=monitored_split,
                monitored_name=("composite" if selected_mode == "composite" else cfg.train.metric_name),
                stop_reason=stop_reason,
                used_best_checkpoint_for_test=use_best_for_eval,
                runtime=runtime,
                selection_mode=selected_mode,
                primary_best_epoch=best_epoch,
                primary_best_metric=best_metric,
                best_per_criterion=es_best,
                bad_epochs_at_stop=bad_epochs,
                patience=cfg.train.early_stopping_patience,
            )

            print(f"Saved outputs to: {out_dir}")
            print("[OK] KMLEE-BAM integrated training finished.")

        control_barrier(runtime)

    except KeyboardInterrupt:
        if runtime.is_main_process:
            print("[interrupt] KeyboardInterrupt received.")
            if trainer is not None:
                out_dir = Path(cfg.train.out_dir)
                save_checkpoint(
                    out_dir / "checkpoint_interrupt.pt",
                    epoch=stopped_epoch,
                    global_step=global_step,
                    step_in_epoch=None,
                    cfg=cfg,
                    system=trainer.system,
                    criterion=trainer.criterion,
                    optimizer=trainer.optimizer,
                    scheduler=trainer.scheduler,
                    trainer=trainer,
                    history=history,
                    best_metric=best_metric,
                    best_epoch=best_epoch,
                    stop_reason="keyboard_interrupt",
                )
                print(f"[interrupt] saved checkpoint to {out_dir / 'checkpoint_interrupt.pt'}")
        raise

    finally:
        cleanup_runtime(runtime)


# =====================================================================
# CLI
# =====================================================================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run training for the module-token ordinal BAM model with Lie-action decoder (DDP compatible).")
    p.add_argument("--config", type=str, required=True, help="Path to training config JSON.")
    return p


def main() -> None:
    args = build_argparser().parse_args()
    cfg = load_config(args.config)
    run_training(cfg)


if __name__ == "__main__":
    main()
