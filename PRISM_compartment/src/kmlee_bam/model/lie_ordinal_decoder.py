from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


def inverse_softplus(y: float) -> float:
    y = float(y)
    if y <= 0:
        raise ValueError("inverse_softplus requires y > 0.")
    return math.log(math.expm1(y))


# ======================================================================
# Output container
# ======================================================================
@dataclass
class LieActionOrdinalDecoderOutput:
    score: torch.Tensor
    base_score: torch.Tensor
    tech_score: torch.Tensor
    state_score: torch.Tensor
    thresholds: torch.Tensor
    cumprob: torch.Tensor
    probs: torch.Tensor
    nll_per_gene: Optional[torch.Tensor]
    nll_per_cell: Optional[torch.Tensor]
    # Biological-covariate factorization (2026-06): explicit SEX main-effect
    # baseline, additive like tech. ``None`` when the sex term is disabled.
    sex_score: Optional[torch.Tensor] = None
    # PRISM E2E explicit common/personal/pathology score.  Kept separate from
    # ``state_score`` so fingerprints never have to guess which legacy latent
    # coordinates carried disease.  None for legacy runs.
    precision_score: Optional[torch.Tensor] = None


# ======================================================================
# Coefficient head: z_perp -> gamma
# ======================================================================
class LieCoefficientHead(nn.Module):
    """
    Maps z_perp to Lie coefficients gamma_i in R^M.
    """

    def __init__(
        self,
        d_z: int,
        n_generators: int,
        *,
        hidden_dim: Optional[int] = None,
        dropout: float = 0.1,
        init_std: float = 0.02,
        gamma_scale: float = 0.1,
        n_coeff_per_generator: int = 1,
    ) -> None:
        super().__init__()

        if hidden_dim is None:
            hidden_dim = max(64, 2 * d_z)
        if n_coeff_per_generator < 1:
            raise ValueError("n_coeff_per_generator must be >= 1.")

        self.d_z = d_z
        self.n_generators = n_generators
        # Per-rank coefficients (zero-origin / capacity work, 2026-06): when
        # n_coeff_per_generator == 1 this is the original head (one scalar per
        # generator). When == R (the Lie rank), each rank-1 component of every
        # generator gets its OWN per-cell coefficient, so the R components are
        # no longer welded to a single shared scalar — turning rank into real
        # per-cell degrees of freedom. See
        # doc/zero_origin_and_capacity_design_2026-06-01.md.
        self.n_coeff_per_generator = int(n_coeff_per_generator)
        self.hidden_dim = hidden_dim
        self.init_std = init_std
        self.gamma_scale = gamma_scale

        self.fc1 = nn.Linear(d_z, hidden_dim)
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(hidden_dim, n_generators * self.n_coeff_per_generator)

        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.fc1.weight, std=self.init_std)
        nn.init.zeros_(self.fc1.bias)

        # zero-init last layer so decoder starts near baseline-only
        nn.init.normal_(self.fc2.weight, mean=0.0, std=1e-3)
        nn.init.zeros_(self.fc2.bias)

    def forward(self, z_perp: torch.Tensor) -> torch.Tensor:
        h = self.fc1(z_perp)
        h = F.gelu(h)
        h = self.dropout(h)
        gamma = self.fc2(h)
        gamma = self.gamma_scale * gamma
        if self.n_coeff_per_generator > 1:
            # [B, M*R] -> [B, M, R]
            gamma = gamma.reshape(
                gamma.shape[0], self.n_generators, self.n_coeff_per_generator
            )
        return gamma


# ======================================================================
# Score Residual Mixer (v6_min_mixer)
# ======================================================================
class ScoreResidualMixer(nn.Module):
    """
    Replace the fixed sum `score = base + tech + state` with a per-cell
    soft-gated mix:

        gates = gate_clamp * softmax(MLP(z_perp))   # [B, 3], sum == gate_clamp
        score = g_b * base + g_t * tech + g_s * state

    The MLP's last layer is zero-initialised, so at step 0 every gate
    equals gate_clamp * softmax(0,0,0) = gate_clamp/3. With the default
    gate_clamp=3 that is [1, 1, 1] → the mixer is the identity (exact
    match with the unmixed `base + tech + state` decoder). Any deviation
    from identity must be earned by the reconstruction gradient and is
    bounded by regularizers on the loss side (gate_balance + gate_ref_safe
    in V6Trainer). softmax couples the branches — amplifying one forces the
    others down (sum is fixed) — unlike the old independent-sigmoid form.

    Note: only scalar per-branch gates are learned. No gene-axis residual
    is added — that introduces extra capacity that's hard to disentangle
    from the existing branches. The mixer is intentionally minimal.

    See doc/v6_min_mixer_design_2026-05-22.md for full rationale.
    """

    def __init__(
        self,
        d_z: int,
        *,
        hidden_dim: int = 64,
        gate_clamp: float = 3.0,
    ) -> None:
        super().__init__()
        self.d_z = int(d_z)
        self.hidden_dim = int(hidden_dim)
        # gate_clamp is now the TOTAL gate budget (sum across 3 branches).
        # softmax(0,0,0) = [1/3, 1/3, 1/3], multiplied by gate_clamp=3 gives
        # the identity gate [1, 1, 1]. Larger gate_clamp allows wider total
        # amplification while still forcing branches to compete for budget.
        self.gate_clamp = float(gate_clamp)

        self.gate_mlp = nn.Sequential(
            nn.Linear(d_z, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 3),
        )
        self._init_identity()

    def _init_identity(self) -> None:
        # Hidden layer: small std init
        nn.init.normal_(self.gate_mlp[0].weight, std=0.02)
        nn.init.zeros_(self.gate_mlp[0].bias)
        # Output layer: zero init → softmax(0,0,0) = [1/3,1/3,1/3], so with
        # gate_clamp=3 every gate is 1.0 and score == base + tech + state.
        nn.init.zeros_(self.gate_mlp[2].weight)
        nn.init.zeros_(self.gate_mlp[2].bias)

    def forward(
        self,
        z_perp: torch.Tensor,         # [B, d_z]
        base_score: torch.Tensor,      # [B, G]
        tech_score: torch.Tensor,      # [B, G]
        state_score: torch.Tensor,     # [B, G]
    ) -> tuple[torch.Tensor, torch.Tensor]:
        # Softmax-competitive gates (v6_min_mixer review fix 2026-05-22).
        # Original sigmoid form let `base` saturate at clamp=2 because each
        # branch was independent: amplifying base did NOT force tech/state
        # down, and the gate_balance prior was too weak to push back.
        #
        # Softmax form: sum of gates == gate_clamp always. Inflating one
        # branch *requires* deflating others, so degenerate single-branch
        # amplification is structurally impossible. Identity init preserved
        # because softmax(0,0,0) = [1/3, 1/3, 1/3] and gate_clamp=3 brings
        # each back to 1.0.
        gate_logits = self.gate_mlp(z_perp)                       # [B, 3]
        gates = self.gate_clamp * torch.softmax(gate_logits, dim=-1)
        g_base, g_tech, g_state = gates.unbind(-1)                # each [B]

        score = (
            g_base.unsqueeze(-1)  * base_score
            + g_tech.unsqueeze(-1) * tech_score
            + g_state.unsqueeze(-1) * state_score
        )
        return score, gates


# ======================================================================
# Lie-action additive ordinal decoder
# ======================================================================
class LieActionOrdinalDecoder(nn.Module):
    """
    Gene-level ordinal decoder with baseline deformation by Lie generators.

    For each cell i:
        b_i = beta^(c)_{:, c_i}
        t_i = beta^(tech)_{:, t_i}
        gamma_i = g(z_perp_i)

        delta_i = sum_m gamma_{im} * (v_m^T b_i) * u_m
        s_i = b_i + t_i + delta_i

    Here each generator is rank-1:
        A_m = u_m v_m^T
    so
        A_m b_i = u_m (v_m^T b_i)
    """

    def __init__(
        self,
        n_genes: int,
        n_celltypes: int,
        n_tech: int,
        d_z: int,
        n_bins: int,
        *,
        n_sex: Optional[int] = None,
        n_generators: int = 32,
        coeff_hidden_dim: Optional[int] = None,
        coeff_dropout: float = 0.1,
        gamma_scale: float = 0.1,
        generator_init_scale: float = 1e-3,
        translation_init_scale: float = 1e-2,
        use_affine_translation: bool = True,
        use_direct_state_head: bool = False,
        direct_state_init_scale: float = 1e-2,
        module_masks: Optional[torch.Tensor] = None,   # [M, G] or None
        # Full-module generator contract (2026-07).  A singleton module has a
        # one-dimensional support, so rank > 1 is algebraically redundant.
        # ``None`` preserves the historical uniform-rank behaviour.
        singleton_lie_rank: Optional[int] = None,
        # Optional scale stabilization for overlapping module masks:
        #   v input  *= 1/sqrt(module size)
        #   u/a out  *= 1/sqrt(number of modules covering each gene)
        # This keeps initial projection/output variance comparable across
        # small/large and low/high-overlap modules. ``none`` is byte-compatible
        # with all historical checkpoints.
        module_generator_normalization: str = "none",
        center_thresholds: bool = True,
        threshold_init_spacing: float = 1.0,
        prob_eps: float = 1e-8,
        init_std: float = 0.02,
        bin_marginal_freq: Optional[torch.Tensor] = None,
        bin_marginal_cum_eps: float = 1e-3,
        # v6_minimal: rank-R Lie generators. Each generator is a sum of R
        # outer products `u_{m,k} v_{m,k}^T`. R=1 reproduces the v5 decoder
        # exactly. R>1 doubles/triples/etc. the effective generator subspace.
        # See doc/v6_minimal_design_2026-05-21.md §M1.
        lie_max_rank: int = 1,
        # zero-origin / capacity work (2026-06): per-rank coefficients. When
        # True (and lie_max_rank > 1), the coefficient head emits one scalar
        # per (generator, rank) pair instead of one per generator, so the R
        # rank-1 pieces of each generator are independently steerable per cell.
        # Default False reproduces the shared-coefficient decoder exactly.
        # See doc/zero_origin_and_capacity_design_2026-06-01.md.
        per_rank_coeff: bool = False,
        # v6_min_mixer: ScoreResidualMixer. When True, replaces the fixed
        # `score = base + tech + state` sum with a per-cell soft-gated mix
        # initialised to the identity. See doc/v6_min_mixer_design_2026-05-22.md.
        use_score_residual_mixer: bool = False,
        mixer_hidden_dim: int = 64,
        mixer_gate_clamp: float = 3.0,   # softmax identity needs gate_clamp/3 == 1
        # v23: pathology-axis INTERACTION terms. Project z onto k learned pathology
        # coords, add explicit products (pairwise..interaction_max_order) each with a
        # learned gene pattern. Default OFF + zero-init patterns => byte-identical when
        # disabled. Injects CURVATURE (z-dependent Jacobian) + biological synergy.
        use_pathology_interactions: bool = False,
        n_interaction_axes: int = 4,
        interaction_max_order: int = 2,   # 2=pairwise, 3=+triples, 4=+quad
        # v26/v27: pathology-TIED, celltype-conditioned decoder (the *principled* curvature).
        # Unlike v23's FREE interaction_proj, here p is TIED to named pathology axes (soft-tied
        # to donor labels by a loss in the objectives) and the gene patterns are CELLTYPE-
        # CONDITIONED via a low-rank correction: the SAME Braak axis drives DIFFERENT modules in
        # microglia vs neuron. v26 = single-axis only; v27 = + pairwise products p_k·p_l (the
        # synergy = the real disease curvature; this is the term LATE benefits from). Everything
        # zero-init => no-op until trained; default OFF => byte-identical.
        use_pathology_decoder: bool = False,
        n_pathology_axes: int = 4,
        pathology_pairwise: bool = False,   # False = v26 (single-axis), True = v27 (+ pairwise synergy)
        pathology_rank: int = 8,            # low-rank dim of the celltype-specific pattern correction
        # An integrated learned-rank model may allocate more columns while its
        # Phase-I forward must retain the canonical fixed-rank RNG trajectory.
        # Draw this prefix from the main RNG and draw only the extra columns in
        # a forked RNG scope.  Zero preserves every historical constructor.
        pathology_rank_rng_compat_prefix: int = 0,
        # v35: disease-GATED basis. score += (s ⊙ sigmoid(gate_slope·p)) @ disease_base; the sigmoid
        # gate makes ∂score/∂z depend on z => G_disease ≠ G_normal (curvature). s is L1-penalized in
        # the objective so the model self-selects axes. zero-init => no-op; OFF => byte-identical.
        # Requires use_pathology_decoder (provides p). See doc/v35_disease_gate_design.md.
        use_disease_gate: bool = False,
        gate_slope: float = 1.0,
        gate_rank: int = 0,                 # 0 = celltype-SHARED disease_base; >0 + flag adds low-rank corr
        gate_celltype_conditioned: bool = False,
        # learnable Lie rank (2026-07): over-provision lie_max_rank + group-lasso prune the
        # rank-1 pieces so the model self-selects effective rank. OFF => byte-identical.
        learnable_lie_rank: bool = False,
        lambda_lie_rank_sparse: float = 0.0,
        lie_rank_warmup_epochs: int = 2,
    ) -> None:
        super().__init__()

        if n_genes <= 0:
            raise ValueError("n_genes must be positive.")
        if n_celltypes <= 0:
            raise ValueError("n_celltypes must be positive.")
        if n_tech <= 0:
            raise ValueError("n_tech must be positive.")
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if n_bins < 2:
            raise ValueError("n_bins must be >= 2.")
        if n_generators <= 0:
            raise ValueError("n_generators must be positive.")
        if threshold_init_spacing <= 0:
            raise ValueError("threshold_init_spacing must be positive.")
        if prob_eps <= 0:
            raise ValueError("prob_eps must be positive.")
        if lie_max_rank < 1:
            raise ValueError("lie_max_rank must be >= 1.")
        if singleton_lie_rank is not None and not (
            1 <= int(singleton_lie_rank) <= int(lie_max_rank)
        ):
            raise ValueError(
                "singleton_lie_rank must be between 1 and lie_max_rank."
            )
        module_generator_normalization = str(
            module_generator_normalization
        ).lower()
        if module_generator_normalization not in ("none", "size_overlap"):
            raise ValueError(
                "module_generator_normalization must be 'none' or "
                "'size_overlap'."
            )
        if mixer_hidden_dim < 1:
            raise ValueError("mixer_hidden_dim must be >= 1.")
        if mixer_gate_clamp <= 0:
            raise ValueError("mixer_gate_clamp must be > 0.")
        if use_disease_gate and not use_pathology_decoder:
            raise ValueError(
                "use_disease_gate requires use_pathology_decoder=True (needs pathology_head's p)."
            )
        if isinstance(pathology_rank_rng_compat_prefix, bool) or int(
            pathology_rank_rng_compat_prefix
        ) < 0:
            raise ValueError(
                "pathology_rank_rng_compat_prefix must be a non-negative integer"
            )
        if int(pathology_rank_rng_compat_prefix) > int(pathology_rank):
            raise ValueError(
                "pathology_rank_rng_compat_prefix cannot exceed pathology_rank"
            )
        self.n_genes = n_genes
        self.n_celltypes = n_celltypes
        self.n_tech = n_tech
        self.d_z = d_z
        self.n_bins = n_bins
        self.n_generators = n_generators
        self.lie_max_rank = int(lie_max_rank)
        self.singleton_lie_rank = (
            None
            if singleton_lie_rank is None
            else int(singleton_lie_rank)
        )
        self.module_generator_normalization = (
            module_generator_normalization
        )
        # learnable Lie rank (group-lasso on rank-1 pieces). OFF => no stash => byte-identical.
        self.learnable_lie_rank = bool(learnable_lie_rank)
        self.lambda_lie_rank_sparse = float(lambda_lie_rank_sparse)
        self.lie_rank_warmup_epochs = int(lie_rank_warmup_epochs)
        # v23 interaction terms (default OFF). HONEST SCOPE: interaction_proj is a FREELY
        # LEARNED Linear(d_z->k); its k coords are NOT guaranteed to be Braak/Thal/LATE/
        # Lewy — they are just "k decoder-learned interaction coordinates" whose 2/3/4-way
        # products add CURVATURE/expressivity. This gives geometry, NOT named-pathology
        # interpretability (and the curvature may land in non-disease directions). For
        # interpretable + disease-aligned curvature, tie these coords to the pathology_aux
        # head's axis projection (a later refinement). zero-init => no DISRUPTION at step 0,
        # but the patterns DO learn (signal AND noise) once enabled, so order N adds REAL
        # capacity — validate vs a baseline on module_cross / leakage / split-stability,
        # never assume "unused stays 0".
        self.use_pathology_interactions = bool(use_pathology_interactions)
        if self.use_pathology_interactions:
            import itertools as _it
            k = int(n_interaction_axes)
            self.interaction_proj = nn.Linear(d_z, k, bias=False)        # z -> k LEARNED coords (NOT named axes)
            combos = []
            for _order in range(2, int(interaction_max_order) + 1):
                combos += list(_it.combinations(range(k), _order))      # pairwise .. max_order
            self._interaction_combos = combos                            # list[tuple[int,...]]
            self.interaction_patterns = nn.Parameter(torch.zeros(len(combos), n_genes))
        # v26/v27: pathology-tied, celltype-conditioned decoder term (default OFF, all zero-init).
        self.use_pathology_decoder = bool(use_pathology_decoder)
        if self.use_pathology_decoder:
            import itertools as _it2
            na = int(n_pathology_axes); rk = int(pathology_rank)
            self.n_pathology_axes = na
            self.pathology_pairwise = bool(pathology_pairwise)
            self.pathology_head = nn.Linear(d_z, na)                       # z_clean -> p (tied to labels via loss)
            nn.init.normal_(self.pathology_head.weight, std=1e-2); nn.init.zeros_(self.pathology_head.bias)
            # single-axis pattern = shared base [na,G] (0) + celltype low-rank corr U[ct,na,rk] @ V[rk,G]
            self.path_single_base = nn.Parameter(torch.zeros(na, n_genes))
            self.path_single_U = nn.Parameter(torch.zeros(n_celltypes, na, rk))
            compat = int(pathology_rank_rng_compat_prefix)
            if 0 < compat < rk:
                # The first draw is byte-identical to the canonical rank-R
                # constructor.  The fork restores the CPU RNG after the extra
                # capacity is initialised, so every subsequently constructed
                # Phase-I parameter sees the canonical RNG state.
                path_v_prefix = torch.randn(compat, n_genes)
                with torch.random.fork_rng(devices=[]):
                    path_v_extra = torch.randn(rk - compat, n_genes)
                path_v_init = torch.cat((path_v_prefix, path_v_extra), dim=0)
            else:
                path_v_init = torch.randn(rk, n_genes)
            self.path_V = nn.Parameter(path_v_init * 1e-3)                 # shared low-rank gene basis
            # Runtime-only curriculum contribution.  Persistent so a resumed
            # integrated Phase-I/II run cannot silently reactivate this route.
            self.register_buffer(
                "pathology_curriculum_scale",
                torch.tensor(1.0, dtype=torch.float32),
                persistent=True,
            )
            if self.pathology_pairwise:                                    # v27: + pairwise synergy
                self._path_pairs = list(_it2.combinations(range(na), 2))   # 6 pairs for na=4
                self.path_pair_base = nn.Parameter(torch.zeros(len(self._path_pairs), n_genes))
                self.path_pair_U = nn.Parameter(torch.zeros(n_celltypes, len(self._path_pairs), rk))
            # v35: disease-GATED basis (uses pathology_head's p above). Built ONLY when
            # use_disease_gate => OFF leaves zero new params => byte-identical. All zero-init.
            self.use_disease_gate = bool(use_disease_gate)
            self.gate_celltype_conditioned = bool(gate_celltype_conditioned)
            if self.use_disease_gate:
                self.gate_slope = float(gate_slope)
                self.gate_rank = int(gate_rank)
                # per-axis strength s = softplus(raw); init small (s≈0.1) so curvature grows in;
                # L1 (objective) prunes dead axes. math.log(expm1(0.1)) = inverse softplus of 0.1.
                self.gate_strength_raw = nn.Parameter(
                    torch.full((na,), float(math.log(math.expm1(0.1))))
                )
                # per-axis gate threshold: gate = sigmoid(slope·(p - bias)). zero-init => gate=½ at
                # p=bias=0 (current behavior); learns to SHIFT the open/closed boundary so the gate
                # shuts in the normal zone (p≈0) and opens in disease (p high) => curvature lands in
                # the DISEASE region, not the normal one (review finding #3).
                self.gate_bias = nn.Parameter(torch.zeros(na))
                self.disease_base = nn.Parameter(torch.zeros(na, n_genes))               # ZERO-INIT
                if self.gate_celltype_conditioned and self.gate_rank > 0:
                    self.disease_U = nn.Parameter(torch.zeros(n_celltypes, na, self.gate_rank))  # ZERO-INIT
                    self.disease_V = nn.Parameter(torch.randn(self.gate_rank, n_genes) * 1e-3)
        # Per-rank coefficients only have an effect when there is more than one
        # rank component to steer; at rank 1 the flag is a silent no-op so that
        # gamma stays [B, M] and every downstream path is byte-identical.
        self.per_rank_coeff = bool(per_rank_coeff) and self.lie_max_rank > 1
        self._n_coeff_per_generator = self.lie_max_rank if self.per_rank_coeff else 1
        # learnable rank lives in the per_rank branch; fail loudly instead of a silent no-op.
        if self.learnable_lie_rank and not self.per_rank_coeff:
            raise ValueError(
                "learnable_lie_rank=True requires per_rank_coeff=True and lie_max_rank>1 "
                "(the rank group-lasso is computed only in the per-rank-coefficient branch)."
            )
        self.center_thresholds = center_thresholds
        self.prob_eps = prob_eps
        self.init_std = init_std
        self.translation_init_scale = translation_init_scale
        self.use_affine_translation = bool(use_affine_translation)
        self.use_direct_state_head = bool(use_direct_state_head)
        self.direct_state_init_scale = direct_state_init_scale
        if self.use_direct_state_head:
            self.direct_state_head = nn.Linear(d_z, n_genes)
            nn.init.normal_(self.direct_state_head.weight, std=direct_state_init_scale)
            nn.init.zeros_(self.direct_state_head.bias)
        else:
            self.direct_state_head = None
        
        # baseline terms
        self.celltype_baseline = nn.Embedding(n_celltypes, n_genes)
        self.tech_baseline = nn.Embedding(n_tech, n_genes)

        # Biological-covariate factorization (2026-06): an OPTIONAL additive SEX
        # baseline, exactly like ``tech_baseline``. It gives the sex MAIN EFFECT
        # on expression its own home so the encoder need not put it in z, WITHOUT
        # erasing sex (sex×disease interactions remain expressible in z and
        # queryable at readout). Index convention: raw ``sex_id`` in
        # {-1=unknown, 0, 1, ...} is shifted +1 (so -1→row 0 "unknown"), hence
        # ``n_sex`` must be (max_sex_id + 2); e.g. n_sex=3 for {unknown, F, M}.
        # Zero-init ⇒ no-op at step 0. ``n_sex=None`` disables the term and the
        # decoder is byte-identical to base+tech+state. The same pattern extends
        # to age / ancestry / tissue as future additive baselines.
        self.n_sex = int(n_sex) if n_sex else 0
        if self.n_sex > 0:
            self.sex_baseline = nn.Embedding(self.n_sex, n_genes)
        else:
            self.sex_baseline = None

        # gamma(z_perp)
        self.coeff_head = LieCoefficientHead(
            d_z=d_z,
            n_generators=n_generators,
            hidden_dim=coeff_hidden_dim,
            dropout=coeff_dropout,
            init_std=init_std,
            gamma_scale=gamma_scale,
            n_coeff_per_generator=self._n_coeff_per_generator,
        )

        # rank-R generators (v6_minimal): each generator m is a sum of R
        # outer products u_{m,k} v_{m,k}^T. No learnable gate — all R
        # components are active at all times. R=1 reproduces v5 exactly.
        R = self.lie_max_rank
        self.generator_u = nn.Parameter(
            torch.randn(n_generators, R, n_genes) * generator_init_scale
        )
        self.generator_v = nn.Parameter(
            torch.randn(n_generators, R, n_genes) * generator_init_scale
        )

        # `generator_a` is the affine translation, NOT part of the rank-k
        # outer-product decomposition. Keep its [M, G] shape unchanged.
        self.generator_a = nn.Parameter(
            torch.randn(n_generators, n_genes) * translation_init_scale
        )

        # Optional ScoreResidualMixer. Identity-initialised: at step 0,
        # behavior is exactly `score = base + tech + state` (unchanged).
        self.use_score_residual_mixer = bool(use_score_residual_mixer)
        if self.use_score_residual_mixer:
            self.score_mixer = ScoreResidualMixer(
                d_z=d_z,
                hidden_dim=int(mixer_hidden_dim),
                gate_clamp=float(mixer_gate_clamp),
            )
        else:
            self.score_mixer = None

        if module_masks is None:
            self.register_buffer("module_masks", torch.ones(n_generators, n_genes))
        else:
            if module_masks.shape != (n_generators, n_genes):
                raise ValueError(
                    f"module_masks must have shape ({n_generators}, {n_genes}), "
                    f"got {tuple(module_masks.shape)}."
                )
            self.register_buffer("module_masks", module_masks.float())

        # Hard generator-budget switch.  This is deliberately binary and is
        # not a trainable gate: the architecture controller changes it only
        # between optimizer/validation phases.  One row controls the complete
        # generator contract (Lie linear action and its shared-coefficient
        # affine translation), preventing a nominally removed generator from
        # leaking through ``generator_a``.
        self.register_buffer(
            "generator_active_mask",
            torch.ones(
                n_generators,
                dtype=self.module_masks.dtype,
                device=self.module_masks.device,
            ),
        )

        # These buffers are deterministic functions of module_masks + config,
        # so they need not be stored in checkpoints.  Keeping them
        # non-persistent preserves old state_dict contracts while new runs are
        # fully reproducible from their frozen config/registry.
        support = (self.module_masks > 0).to(self.module_masks.dtype)
        module_sizes = support.sum(dim=1)
        overlap_degree = support.sum(dim=0)
        rank_mask = torch.ones(
            n_generators,
            self.lie_max_rank,
            dtype=self.module_masks.dtype,
            device=self.module_masks.device,
        )
        singleton_rows = module_sizes == 1
        if self.singleton_lie_rank is not None and bool(singleton_rows.any()):
            rank_mask[
                singleton_rows, self.singleton_lie_rank :
            ] = 0.0

        if self.module_generator_normalization == "size_overlap":
            input_scale = module_sizes.clamp_min(1.0).rsqrt()
            output_scale = overlap_degree.clamp_min(1.0).rsqrt()
            # A gene outside every selected module must stay outside the Lie
            # and translation paths instead of receiving an arbitrary scale.
            output_scale = torch.where(
                overlap_degree > 0,
                output_scale,
                torch.zeros_like(output_scale),
            )
        else:
            input_scale = torch.ones_like(module_sizes)
            output_scale = torch.ones_like(overlap_degree)

        self.register_buffer(
            "generator_rank_mask", rank_mask, persistent=False
        )
        self.register_buffer(
            "generator_input_scale", input_scale, persistent=False
        )
        self.register_buffer(
            "generator_output_scale", output_scale, persistent=False
        )
        self.register_buffer(
            "module_overlap_degree", overlap_degree, persistent=False
        )
        self.n_singleton_generators = int(singleton_rows.sum().item())
        self.n_module_covered_genes = int(
            (overlap_degree > 0).sum().item()
        )

        # thresholds (default = uniform init; may be overwritten by quantile
        # init below if `bin_marginal_freq` is provided)
        self.raw_threshold_start = nn.Parameter(torch.zeros(n_genes, 1))
        if n_bins > 2:
            raw_delta_init = inverse_softplus(threshold_init_spacing)
            self.raw_threshold_deltas = nn.Parameter(
                torch.full((n_genes, n_bins - 2), raw_delta_init)
            )
        else:
            self.raw_threshold_deltas = None

        self.reset_parameters()

        # Optional data-driven threshold initialisation. See
        # `doc/threshold_quantile_init.md` for motivation and method.
        if bin_marginal_freq is not None:
            if self.center_thresholds:
                # Quantile init places thresholds at the correct absolute
                # positions; subtracting the mean afterwards would destroy
                # that calibration. Auto-disable to avoid silent bugs.
                self.center_thresholds = False
            self._init_thresholds_from_marginal(
                bin_marginal_freq,
                cum_eps=float(bin_marginal_cum_eps),
            )

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
        """Load legacy checkpoints as an all-active generator architecture.

        ``generator_active_mask`` was introduced after the historical CT64 and
        PRISM checkpoints had already been frozen.  Synthesising an all-one
        value only when that exact key is absent preserves strict loading for
        every other tensor while keeping old checkpoints byte-compatible.
        New search checkpoints store the actual hard mask normally.
        """

        key = prefix + "generator_active_mask"
        inserted_legacy_default = key not in state_dict
        if inserted_legacy_default:
            state_dict[key] = torch.ones_like(self.generator_active_mask)
        try:
            super()._load_from_state_dict(
                state_dict,
                prefix,
                local_metadata,
                strict,
                missing_keys,
                unexpected_keys,
                error_msgs,
            )
        finally:
            if inserted_legacy_default:
                state_dict.pop(key, None)

    def reset_parameters(self) -> None:
        nn.init.zeros_(self.celltype_baseline.weight)
        nn.init.zeros_(self.tech_baseline.weight)
        if self.sex_baseline is not None:
            nn.init.zeros_(self.sex_baseline.weight)
        nn.init.zeros_(self.raw_threshold_start)
        # coeff_head and generator factors already initialized

    @torch.no_grad()
    def set_generator_active_mask(self, mask: torch.Tensor) -> None:
        """Install an exact binary generator mask.

        The mask is checkpoint-persistent but never receives gradients.  Use
        this method instead of mutating the buffer directly so malformed or
        soft architecture gates fail closed.
        """

        mask = torch.as_tensor(
            mask,
            device=self.generator_active_mask.device,
        )
        if mask.ndim != 1 or int(mask.numel()) != int(self.n_generators):
            raise ValueError(
                "generator active mask must have shape "
                f"({self.n_generators},), got {tuple(mask.shape)}."
            )
        if mask.dtype == torch.bool:
            exact = mask
        else:
            if not bool(torch.isfinite(mask).all()):
                raise ValueError("generator active mask contains non-finite values.")
            if not bool(((mask == 0) | (mask == 1)).all()):
                raise ValueError(
                    "generator active mask must be exactly binary; "
                    "continuous/soft gates are not permitted."
                )
            exact = mask.bool()
        if not bool(exact.any()):
            raise ValueError("at least one generator must remain active.")
        self.generator_active_mask.copy_(
            exact.to(dtype=self.generator_active_mask.dtype)
        )

    @torch.no_grad()
    def reset_generator_active_mask(self) -> None:
        self.generator_active_mask.fill_(1.0)

    @property
    def n_active_generators(self) -> int:
        return int((self.generator_active_mask > 0).sum().item())

    def _sex_index(self, sex_id: torch.Tensor) -> torch.Tensor:
        """Map raw sex ids {-1=unknown, 0, 1, ...} → embedding rows {0, 1, 2, ...}.

        Shift by +1 so the missing-label sentinel ``-1`` lands on the dedicated
        "unknown" row 0, then clamp into ``[0, n_sex-1]`` so an out-of-range id can
        never index out of bounds. Robust and differentiable-free (index only).
        """
        return (sex_id.long() + 1).clamp_(min=0, max=self.n_sex - 1)

    @torch.no_grad()
    def _init_thresholds_from_marginal(
        self,
        bin_marginal_freq: torch.Tensor,
        *,
        cum_eps: float,
    ) -> None:
        """
        Initialise (raw_threshold_start, raw_threshold_deltas) so that at
        score = 0 the decoder reproduces the per-gene marginal bin
        distribution given by `bin_marginal_freq`.

        Theory
        ------
        With the cumulative-ordinal head,
            P(bin <= k | score) = sigmoid(threshold_k - score).
        At init time, celltype_baseline, tech_baseline are 0 and the Lie
        state contribution is near zero, so score ~= 0. We want
            sigmoid(threshold_k) = cum_P_k,
        hence threshold_k = logit(cum_P_k) for k in {0, ..., K-2}.

        Implementation
        --------------
        `cum_P_k` is clamped into [cum_eps, 1 - cum_eps] to keep logit
        finite. Differences threshold_{k+1} - threshold_k are guaranteed
        positive (since cum_P is non-decreasing) and are stored via
        inverse softplus to match the existing parametrisation
            thresholds = cumsum(softplus(deltas), start = raw_start).
        """
        if bin_marginal_freq.ndim != 2:
            raise ValueError("bin_marginal_freq must be 2-D [G, K].")
        if bin_marginal_freq.shape != (self.n_genes, self.n_bins):
            raise ValueError(
                f"bin_marginal_freq must have shape "
                f"({self.n_genes}, {self.n_bins}), got "
                f"{tuple(bin_marginal_freq.shape)}."
            )

        freq = bin_marginal_freq.to(dtype=torch.float64)
        freq = freq / freq.sum(dim=-1, keepdim=True).clamp_min(1e-12)
        cum = torch.cumsum(freq, dim=-1)                              # [G, K]
        cum_thresh = cum[:, : self.n_bins - 1].clamp(
            min=float(cum_eps),
            max=1.0 - float(cum_eps),
        )                                                              # [G, K-1]
        thresholds = torch.log(cum_thresh / (1.0 - cum_thresh))        # [G, K-1]

        # raw_threshold_start = threshold_0
        start_target = thresholds[:, :1].to(self.raw_threshold_start.dtype)
        self.raw_threshold_start.copy_(start_target)

        if self.n_bins > 2 and self.raw_threshold_deltas is not None:
            deltas = thresholds[:, 1:] - thresholds[:, :-1]            # [G, K-2]
            deltas = deltas.clamp_min(1e-6)
            raw_deltas = torch.log(torch.expm1(deltas))
            self.raw_threshold_deltas.copy_(
                raw_deltas.to(self.raw_threshold_deltas.dtype)
            )

    # ==================================================================
    # Helper: masked generators
    # ==================================================================
    def _masked_u(self) -> torch.Tensor:
        # generator_u: [M, R, G]   module_masks: [M, G] → broadcast over R.
        return (
            self.generator_u
            * self.module_masks.unsqueeze(1)
            * self.generator_rank_mask.unsqueeze(-1)
            * self.generator_output_scale.view(1, 1, -1)
        )

    def _masked_v(self) -> torch.Tensor:
        # generator_v: [M, R, G]   module_masks: [M, G] → broadcast over R.
        return (
            self.generator_v
            * self.module_masks.unsqueeze(1)
            * self.generator_rank_mask.unsqueeze(-1)
            * self.generator_input_scale.view(-1, 1, 1)
        )

    def _masked_a(self) -> torch.Tensor:
        # generator_a stays [M, G] (affine translation is not rank-decomposed).
        return (
            self.generator_a
            * self.module_masks
            * self.generator_output_scale.view(1, -1)
        )

    # ==================================================================
    # Thresholds
    # ==================================================================
    def _compute_thresholds(self) -> torch.Tensor:
        start = self.raw_threshold_start  # [G, 1]

        if self.n_bins == 2:
            thresholds = start
        else:
            deltas = F.softplus(self.raw_threshold_deltas)      # [G, K-2]
            offsets = torch.cumsum(deltas, dim=-1)             # [G, K-2]
            thresholds = torch.cat([start, start + offsets], dim=-1)

        if self.center_thresholds and thresholds.shape[-1] > 1:
            thresholds = thresholds - thresholds.mean(dim=-1, keepdim=True)

        return thresholds

    # ==================================================================
    # Lie state score
    # ==================================================================
    def coeff_from_z(self, z_perp: torch.Tensor) -> torch.Tensor:
        if z_perp.ndim != 2 or z_perp.shape[-1] != self.d_z:
            raise ValueError(
                f"z_perp must have shape [B, {self.d_z}], got {tuple(z_perp.shape)}."
            )
        return self.coeff_head(z_perp)  # [B, M]
    def pathology_state_and_p(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """v26/v27 celltype-conditioned, pathology-tied state term.

        Returns ``(term [B, G], p [B, n_axes])``. ``p`` is the named-pathology projection
        (soft-tied to donor labels by the objective); ``term`` is the celltype-conditioned
        gene shift it drives. zero-init ⇒ ``term == 0`` until trained. The pairwise block
        (v27) adds the ``p_k·p_l`` synergy that makes ∂term/∂z depend on z = real curvature.
        """
        p = self.pathology_head(z_perp)                                       # [B, na] (z_perp = z_clean when v25 on)
        term = p @ self.path_single_base                                      # [B, G] shared per-axis pattern
        Uc = self.path_single_U[celltype_id]                                  # [B, na, rk] celltype-specific
        rank_gate_module = getattr(self, "pathology_rank_gate", None)
        rank_gate = (
            rank_gate_module()
            if rank_gate_module is not None
            else torch.ones(
                self.path_V.shape[0],
                dtype=self.path_V.dtype,
                device=self.path_V.device,
            )
        )
        single_coeff = torch.einsum("bk,bkr->br", p, Uc)
        rank_gate = rank_gate.to(
            device=single_coeff.device,
            dtype=single_coeff.dtype,
        )
        term = term + (single_coeff * rank_gate.unsqueeze(0)) @ self.path_V    # [B, G]
        synergy = z_perp.new_zeros(())                                        # v26: no pairwise term
        if getattr(self, "pathology_pairwise", False):
            pp = torch.stack([p[:, i] * p[:, j] for (i, j) in self._path_pairs], dim=1)  # [B, n_pair]
            pair_term = pp @ self.path_pair_base
            Upc = self.path_pair_U[celltype_id]                               # [B, n_pair, rk]
            pair_coeff = torch.einsum("bk,bkr->br", pp, Upc)
            pair_term = pair_term + (pair_coeff * rank_gate.unsqueeze(0)) @ self.path_V  # [B, G]
            term = term + pair_term
            synergy = pair_term.detach().abs().mean()                         # v27 pairwise = the curvature magnitude
        # v35: disease-GATED basis. gate = sigmoid(slope·p) is NONLINEAR in z, so ∂term/∂z depends
        # on z => G_disease ≠ G_normal (curvature), concentrated where pathology p is high. s = per-
        # axis strength (L1-penalized in the objective => self-selection). zero-init disease_base =>
        # no-op until trained; absent (OFF) => byte-identical.
        if getattr(self, "use_disease_gate", False):
            s = F.softplus(self.gate_strength_raw)                             # [na] ≥0
            # gate = sigmoid(slope·(p - bias)): nonlinear in z (curvature); bias shifts the
            # open/closed boundary so the gate shuts in normal (p≈0) and opens in disease.
            gate = torch.sigmoid(self.gate_slope * (p - self.gate_bias))      # [B, na]
            sg = s.unsqueeze(0) * gate                                        # [B, na]
            dterm = sg @ self.disease_base                                    # [B, G]
            if getattr(self, "gate_celltype_conditioned", False) and self.gate_rank > 0:
                Ud = self.disease_U[celltype_id]                             # [B, na, r]
                dterm = dterm + torch.einsum("bk,bkr->br", sg, Ud) @ self.disease_V   # [B, G]
            term = term + dterm
            if getattr(self, "_stash_diag", True):                          # review #4: skip during the
                self._last_gate_absmean = dterm.detach().abs().mean()       # metric-contrast finite-diff
                self._last_gate_strengths = s.detach()                      # so logs = the REAL forward
        # magnitudes stashed for the progress log (the decoder "turning on" from zero-init)
        if getattr(self, "_stash_diag", True):
            self._last_pathdec_absmean = term.detach().abs().mean()
            self._last_synergy_absmean = synergy
        return term, p

    def state_score_from_z_and_base(
        self,
        z_perp: torch.Tensor,
        base_score: torch.Tensor,
        celltype_id: Optional[torch.Tensor] = None,
        *,
        generator_gate: Optional[torch.Tensor] = None,
        include_direct: bool = True,
        include_pathology: bool = True,
        include_interactions: bool = True,
    ) -> torch.Tensor:
        gamma = self.coeff_from_z(z_perp)       # [B, M] or [B, M, R] (per_rank)

        active_gate = self.generator_active_mask
        if generator_gate is not None:
            if generator_gate.ndim != 1 or int(generator_gate.numel()) != int(
                self.n_generators
            ):
                raise ValueError(
                    "generator_gate must have shape "
                    f"({self.n_generators},), got {tuple(generator_gate.shape)}."
                )
            if not bool(torch.isfinite(generator_gate).all()):
                raise ValueError("generator_gate contains non-finite values.")
            active_gate = active_gate * generator_gate.to(
                device=active_gate.device,
                dtype=active_gate.dtype,
            )

        u = self._masked_u()                    # [M, R, G]
        v = self._masked_v()                    # [M, R, G]

        # Rank-R Lie linear action:
        #   action[b, g] = Σ_m gamma[b, m] · Σ_k (b·v_{m,k}) · u_{m,k,g}
        # Each generator is a sum of R rank-1 outer products. R=1 is the
        # v5 decoder; R>1 simply doubles/triples generator subspace.
        if self.per_rank_coeff:
            # Per-rank coefficients: gamma is [B, M, R]; each rank-1 component
            # is scaled by its OWN per-cell scalar instead of a single shared
            # gamma[b, m]. This "un-welds" the R pieces (the shared-coefficient
            # branch below ties them together) so rank becomes genuine per-cell
            # capacity. See doc/zero_origin_and_capacity_design_2026-06-01.md.
            # Inactive singleton ranks must be removed from both the Lie atom
            # and gamma_gen; otherwise their coefficients could still leak
            # through the affine translation path.
            gamma_effective = (
                gamma
                * self.generator_rank_mask.unsqueeze(0)
                * active_gate.view(1, -1, 1)
            )
            proj = torch.einsum("bg,mkg->bmk", base_score, v)         # [B, M, R]
            weighted = gamma_effective * proj                         # [B, M, R]
            linear_action = torch.einsum("bmk,mkg->bg", weighted, u)  # [B, G]
            # learnable-rank group-lasso: local-batch RMS of each rank-1 piece's output
            # contribution c_{m,r} = (‖u_{m,r}‖/√G)·RMS_b(weighted_{b,m,r}). Un-evadable to the
            # (u<->v),(gamma<->u) rescalings. Train-only (skip eval/probe forwards).
            # _stash_diag guard: skip during metric_contrast's finite-diff re-forwards so the real
            # forward's stash is preserved (matches the pathology-p stash convention).
            if self.learnable_lie_rank and getattr(self, "_stash_diag", True):
                _un = u.norm(dim=2) / (u.shape[-1] ** 0.5)            # [M,R]  ‖u_{m,r}‖/√G
                _ms = weighted.float().pow(2).mean(dim=0)             # [M,R]  local RMS^2 (fp32)
                _c = _un * _ms.clamp_min(self.prob_eps).sqrt()       # [M,R]  effective contribution
                self._last_rank_c = _c.detach()                      # readout: TRAIN and VAL (for history)
                # penalty: TRAIN-only (needs grad). None on val/eval => selection loss stays penalty-free.
                self._lie_rank_pen = _c.mean() if (self.training and torch.is_grad_enabled()) else None
            # Affine translation a is [M, G] (not rank-decomposed); drive it by
            # the total per-generator coefficient mass so it is well-defined.
            gamma_gen = gamma_effective.sum(dim=-1)                   # [B, M]
        elif self.lie_max_rank == 1:
            # Backward-compat fast path: collapse trivial rank axis.
            u2 = u.squeeze(1)                   # [M, G]
            v2 = v.squeeze(1)                   # [M, G]
            proj = base_score @ v2.t()          # [B, M]
            gamma_effective = gamma * active_gate.view(1, -1)
            linear_action = (gamma_effective * proj) @ u2   # [B, G]
            gamma_gen = gamma_effective
        else:
            # General R > 1 path, single shared coefficient per generator.
            proj = torch.einsum("bg,mkg->bmk", base_score, v)        # [B, M, R]
            gamma_effective = gamma * active_gate.view(1, -1)
            weighted = gamma_effective.unsqueeze(-1) * proj           # [B, M, R]
            linear_action = torch.einsum("bmk,mkg->bg", weighted, u)  # [B, G]
            gamma_gen = gamma_effective

        state_score = linear_action

        # Affine translation: a(z)
        # baseline이 0이어도 z가 reconstruction score에 직접 영향을 준다.
        if self.use_affine_translation:
            a = self._masked_a()                # [M, G]
            translation = gamma_gen @ a         # [B, G]
            state_score = state_score + translation

        # Optional direct z -> gene branch
        if include_direct and self.direct_state_head is not None:
            _direct = self.direct_state_head(z_perp)
            state_score = state_score + _direct

        # v23 pathology-axis interactions: explicit products of the k pathology coords.
        # The product terms make ∂state_score/∂z depend on z (curvature); zero-init =>
        # no-op until trained. Default OFF => byte-identical.
        if (
            include_interactions
            and getattr(self, "use_pathology_interactions", False)
        ):
            p = self.interaction_proj(z_perp)                                  # [B, k]
            prods = torch.stack(
                [p[:, list(c)].prod(dim=1) for c in self._interaction_combos], dim=1
            )                                                                  # [B, n_combo]
            state_score = state_score + prods @ self.interaction_patterns      # [B, G]

        # v26/v27 pathology-tied, celltype-conditioned term (needs celltype_id; skipped when
        # called without it, e.g. the pullback probe). zero-init => no-op until trained.
        if (
            include_pathology
            and getattr(self, "use_pathology_decoder", False)
            and celltype_id is not None
        ):
            path_term, p_path = self.pathology_state_and_p(z_perp, celltype_id)
            path_scale = getattr(self, "pathology_curriculum_scale", None)
            if path_scale is not None:
                path_term = path_term * path_scale.to(
                    device=path_term.device,
                    dtype=path_term.dtype,
                )
            state_score = state_score + path_term
            self._last_pathology_p = p_path                                    # stashed for the tie loss

        # learnable-rank VALIDITY: per-channel RMS so we can detect if deformation leaks OUT of the
        # penalized Lie action into the unpenalized channels (direct / pathology / translation).
        # Detached (no grad). Guard flags exactly match each branch above => no NameError.
        # _stash_diag guard: skip during metric_contrast finite-diff re-forwards.
        if getattr(self, "learnable_lie_rank", False) and getattr(self, "_stash_diag", True):
            def _rms(t):
                return float(t.detach().float().pow(2).mean().sqrt())
            _cr = {"linear": _rms(linear_action)}
            _cr["translation"] = _rms(translation) if self.use_affine_translation else 0.0
            _cr["direct"] = (
                _rms(_direct)
                if include_direct and self.direct_state_head is not None
                else 0.0
            )
            _cr["pathology"] = (_rms(path_term)
                                if (
                                    include_pathology
                                    and getattr(self, "use_pathology_decoder", False)
                                    and celltype_id is not None
                                )
                                else 0.0)
            self._chan_rms = _cr

        return state_score

    

    # NOTE:
    # 기존 decoder의 state_score_from_z(z_perp) 인터페이스는
    # Lie action에서는 baseline이 필요하므로 더 이상 본질적으로 맞지 않습니다.
    def state_score_from_z(self, z_perp: torch.Tensor) -> torch.Tensor:
        raise RuntimeError(
            "LieActionOrdinalDecoder needs the cell-type baseline. "
            "Use state_score_from_z_and_base(z_perp, base_score) instead."
        )

    def _compute_score_components(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        tech_id: torch.Tensor,
        *,
        generator_gate: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        self._validate_inputs(
            z_perp=z_perp,
            celltype_id=celltype_id,
            tech_id=tech_id,
            y_ord=None,
        )

        base_score = self.celltype_baseline(celltype_id)    # [B, G]
        tech_score = self.tech_baseline(tech_id)            # [B, G]
        state_score = self.state_score_from_z_and_base(
            z_perp,
            base_score,
            celltype_id=celltype_id,
            generator_gate=generator_gate,
        )

        if self.score_mixer is not None:
            # ScoreResidualMixer: identity-initialised, so at step 0
            # score == base + tech + state (same as the no-mixer path).
            # Gates `[B, 3]` are saved as a buffer so trainers can pull
            # them out for the gate_balance + gate_ref_safe regularizers
            # and for diagnostic logging.
            score, gates = self.score_mixer(
                z_perp=z_perp,
                base_score=base_score,
                tech_score=tech_score,
                state_score=state_score,
            )
            # Stash on a non-Parameter attribute the trainer reads.
            self._last_mixer_gates = gates
        else:
            score = base_score + tech_score + state_score
            self._last_mixer_gates = None

        return score, base_score, tech_score, state_score

    # ==================================================================
    # Probability path
    # ==================================================================
    def _score_to_probs(
        self,
        score: torch.Tensor,
        thresholds: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if score.ndim != 2:
            raise ValueError(f"score must have shape [B, G], got {tuple(score.shape)}.")
        if thresholds.ndim != 2:
            raise ValueError(
                f"thresholds must have shape [G, K-1], got {tuple(thresholds.shape)}."
            )

        score_exp = score.unsqueeze(-1)                         # [B, G, 1]
        thr_exp = thresholds.unsqueeze(0)                       # [1, G, K-1]
        cumprob = torch.sigmoid(thr_exp - score_exp)            # [B, G, K-1]

        probs = []
        probs.append(cumprob[..., 0])

        for k in range(1, self.n_bins - 1):
            probs.append(cumprob[..., k] - cumprob[..., k - 1])

        probs.append(1.0 - cumprob[..., -1])
        probs = torch.stack(probs, dim=-1)                      # [B, G, K]

        probs = probs.clamp(min=self.prob_eps)
        probs = probs / probs.sum(dim=-1, keepdim=True).clamp(min=self.prob_eps)
        return cumprob, probs

    def _nll_from_probs(
        self,
        probs: torch.Tensor,
        y_ord: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        gathered = probs.gather(dim=-1, index=y_ord.unsqueeze(-1)).squeeze(-1)
        nll_per_gene = -torch.log(gathered.clamp(min=self.prob_eps))
        nll_per_cell = nll_per_gene.mean(dim=-1)
        return nll_per_gene, nll_per_cell

    # ==================================================================
    # Main forward
    # ==================================================================
    def forward(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        tech_id: torch.Tensor,
        y_ord: Optional[torch.Tensor] = None,
        sex_id: Optional[torch.Tensor] = None,
        precision_score: Optional[torch.Tensor] = None,
        generator_gate: Optional[torch.Tensor] = None,
    ) -> LieActionOrdinalDecoderOutput:
        score, base_score, tech_score, state_score = self._compute_score_components(
            z_perp=z_perp,
            celltype_id=celltype_id,
            tech_id=tech_id,
            generator_gate=generator_gate,
        )

        # Biological-covariate factorization: add the SEX main-effect baseline as
        # an explicit explained term (additive, like tech / outside the
        # base-tech-state mixer). This lets the encoder keep z free of the sex
        # main effect WITHOUT erasing sex. No-op (byte-identical) when the sex
        # term is disabled (``n_sex`` unset) or no ``sex_id`` is supplied.
        sex_score = None
        if self.sex_baseline is not None and sex_id is not None:
            sex_score = self.sex_baseline(self._sex_index(sex_id))
            score = score + sex_score

        if precision_score is not None:
            if precision_score.shape != score.shape:
                raise ValueError(
                    "precision_score must match decoder score shape, got "
                    f"{tuple(precision_score.shape)} vs {tuple(score.shape)}"
                )
            score = score + precision_score

        thresholds = self._compute_thresholds()
        cumprob, probs = self._score_to_probs(score, thresholds)

        nll_per_gene = None
        nll_per_cell = None
        if y_ord is not None:
            self._validate_inputs(
                z_perp=z_perp,
                celltype_id=celltype_id,
                tech_id=tech_id,
                y_ord=y_ord,
            )
            nll_per_gene, nll_per_cell = self._nll_from_probs(probs, y_ord)

        return LieActionOrdinalDecoderOutput(
            score=score,
            base_score=base_score,
            tech_score=tech_score,
            state_score=state_score,
            thresholds=thresholds,
            cumprob=cumprob,
            probs=probs,
            nll_per_gene=nll_per_gene,
            nll_per_cell=nll_per_cell,
            sex_score=sex_score,
            precision_score=precision_score,
        )

    def nll_from_external_score(
        self,
        score: torch.Tensor,
        y_ord: torch.Tensor,
        *,
        detach_thresholds: bool = False,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Evaluate an audited additive score with the decoder's ordinal head.

        PRISM uses this for its target-latent-free branch likelihood.  It shares
        exactly the same learned thresholds and probability convention as the
        full decoder, but does not run the unrestricted cell-state path.
        """

        if score.ndim != 2 or score.shape != y_ord.shape:
            raise ValueError("external score and y_ord must share shape [B,G]")
        thresholds = self._compute_thresholds()
        if detach_thresholds:
            thresholds = thresholds.detach()
        _, probs = self._score_to_probs(score, thresholds)
        return self._nll_from_probs(probs, y_ord)

    # ==================================================================
    # Validation
    # ==================================================================
    def _validate_inputs(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        tech_id: torch.Tensor,
        y_ord: Optional[torch.Tensor],
    ) -> None:
        if z_perp.ndim != 2 or z_perp.shape[-1] != self.d_z:
            raise ValueError(
                f"z_perp must have shape [B, {self.d_z}], got {tuple(z_perp.shape)}."
            )
        B = z_perp.shape[0]

        if celltype_id.ndim != 1 or celltype_id.shape[0] != B or celltype_id.dtype != torch.long:
            raise ValueError("celltype_id must have shape [B] and dtype torch.long.")
        if tech_id.ndim != 1 or tech_id.shape[0] != B or tech_id.dtype != torch.long:
            raise ValueError("tech_id must have shape [B] and dtype torch.long.")

        if y_ord is not None:
            if y_ord.ndim != 2 or y_ord.shape != (B, self.n_genes):
                raise ValueError(
                    f"y_ord must have shape ({B}, {self.n_genes}), got {tuple(y_ord.shape)}."
                )
            if y_ord.dtype != torch.long:
                raise ValueError("y_ord must be torch.long.")
            if torch.any(y_ord < 0) or torch.any(y_ord >= self.n_bins):
                raise ValueError(f"y_ord contains values outside [0, {self.n_bins - 1}].")

    # ==================================================================
    # Same interface as current loss module expects
    # ==================================================================
    def tech_regularizer(self, tech_weights: Optional[torch.Tensor] = None) -> torch.Tensor:
        beta_tech = self.tech_baseline.weight  # [n_tech, G]

        if tech_weights is None:
            weights = torch.full(
                (self.n_tech,),
                1.0 / self.n_tech,
                dtype=beta_tech.dtype,
                device=beta_tech.device,
            )
        else:
            if tech_weights.ndim != 1 or tech_weights.shape[0] != self.n_tech:
                raise ValueError(
                    f"tech_weights must have shape ({self.n_tech},), got {tuple(tech_weights.shape)}."
                )
            weights = tech_weights.to(device=beta_tech.device, dtype=beta_tech.dtype)
            weights = weights / weights.sum().clamp(min=self.prob_eps)

        weighted_mean = (weights.unsqueeze(-1) * beta_tech).sum(dim=0)
        return weighted_mean.pow(2).sum()

    def gauge_regularizer_db_c(
        self,
        state_score: torch.Tensor,
        celltype_id: torch.Tensor,
        tech_id: torch.Tensor,
    ) -> torch.Tensor:
        if state_score.ndim != 2:
            raise ValueError(
                f"state_score must have shape [B, G], got {tuple(state_score.shape)}."
            )

        penalty = torch.zeros((), device=state_score.device, dtype=state_score.dtype)

        unique_celltypes = torch.unique(celltype_id)
        for c in unique_celltypes:
            mask_c = celltype_id == c
            techs_in_c = torch.unique(tech_id[mask_c])

            donor_means = []
            for t in techs_in_c:
                mask_ct = mask_c & (tech_id == t)
                if torch.any(mask_ct):
                    donor_means.append(state_score[mask_ct].mean(dim=0))

            if len(donor_means) == 0:
                continue

            donor_means = torch.stack(donor_means, dim=0)
            balanced_mean_c = donor_means.mean(dim=0)
            penalty = penalty + balanced_mean_c.pow(2).sum()

        return penalty

    # ==================================================================
    # Lie-specific regularizers
    # ==================================================================
    def generator_l1_penalty(self) -> torch.Tensor:
        # Works for both rank-1 ([M, G]) and rank-R ([M, R, G]) shapes since
        # .abs().sum() reduces over all elements.
        u = self._masked_u()
        v = self._masked_v()
        return u.abs().sum() + v.abs().sum()

    def generator_orth_penalty(self) -> torch.Tensor:
        # Treat each (m, k) rank-1 component as a separate generator for the
        # purpose of orthogonality. For rank-1 this reduces exactly to the
        # original [M, G] case.
        u = self._masked_u()                                # [M, R, G]
        v = self._masked_v()                                # [M, R, G]
        M, R, G = u.shape
        u_flat = u.reshape(M * R, G)
        v_flat = v.reshape(M * R, G)

        u_flat = F.normalize(u_flat, dim=-1)
        v_flat = F.normalize(v_flat, dim=-1)

        gram_u = u_flat @ u_flat.t()
        gram_v = v_flat @ v_flat.t()

        n_components = M * R
        eye = torch.eye(n_components, device=u_flat.device, dtype=u_flat.dtype)
        return (gram_u - eye).pow(2).sum() + (gram_v - eye).pow(2).sum()


    def coeff_zero_mean_penalty(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
    ) -> torch.Tensor:
        gamma = self.coeff_from_z(z_perp)  # [B, M] or [B, M, R] (per_rank)
        if gamma.ndim == 3:
            # Per-rank coefficients: penalise every (generator, rank) channel.
            gamma = gamma.reshape(gamma.shape[0], -1)  # [B, M*R]
        penalty = torch.zeros((), device=gamma.device, dtype=gamma.dtype)

        for c in torch.unique(celltype_id):
            mask = celltype_id == c
            if torch.any(mask):
                penalty = penalty + gamma[mask].mean(dim=0).pow(2).sum()

        return penalty
