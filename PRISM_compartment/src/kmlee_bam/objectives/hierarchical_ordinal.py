"""
Hierarchical ordinal loss for KMLEE-BAM.

Decomposes the 7-bin ordinal reconstruction into three levels:

  L1 (binary)        : zero vs nonzero
  L2 (3-way, cond.)  : low / mid / high *given* nonzero
  L3 (EMD, 7-bin)    : earth-mover's distance on the cumulative distribution

The motivation, philosophy, and mathematical rationale live in
`doc/hierarchical_ordinal_loss.md`. This file is a pure-PyTorch implementation
that can be called from any trainer wrapper.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import torch
import torch.nn.functional as F


# ======================================================================
# Config / output dataclasses
# ======================================================================
@dataclass(frozen=True)
class HierarchicalOrdinalConfig:
    """Configuration for the hierarchical ordinal loss."""

    enabled: bool = False

    # Loss weights (kept conservative — main reconstruction stays dominant).
    lambda_zero: float = 0.08
    lambda_group: float = 0.04
    lambda_emd: float = 0.03

    # L1: weighted BCE between p(nonzero) and (y_true > 0).
    # `zero_pos_weight` upweights nonzero cells in the BCE to counter the
    # 80% zero-mass imbalance. Set to 1.0 for plain BCE.
    zero_pos_weight: float = 2.0

    # L2: group definition (tuple of tuples). Each inner tuple lists the
    # bin indices that belong to that group. Group index = position in
    # the outer tuple. Default = ((1,2), (3,4), (5,6)) i.e. low/mid/high.
    # Bin 0 is NOT in any group (handled by L1).
    group_definition: tuple = ((1, 2), (3, 4), (5, 6))

    # L2 class weighting (within-nonzero imbalance: low ~64%, mid ~22%,
    # high ~14% empirically). Without weighting, L2 collapses again to
    # the low group. Modes:
    #   "uniform"      - no weighting (default for backward compat).
    #   "sqrt_inverse" - weight_g = sqrt(N / (n_groups * n_g)), capped.
    #   "inverse"      - weight_g = N / (n_groups * n_g), capped.
    # If `group_class_weights` is set, it overrides `group_weight_mode`.
    group_weight_mode: str = "sqrt_inverse"
    group_weight_cap: float = 4.0
    group_class_weights: Optional[tuple] = None

    # L2.5: within-group binary CE (bin b_lo vs b_hi inside each group).
    # Only applied to groups of exactly 2 bins. Conditional on being in group.
    # 0.0 disables.
    lambda_within_group: float = 0.0

    # Optional linear warm-up over the first `ramp_epochs`. 0 disables ramp.
    ramp_epochs: int = 0

    # π_tech soft-zero (zero-origin work, 2026-06). When a per-(cell,gene)
    # technical-dropout probability `pi_tech` is supplied to the L1 loss, these
    # control how it relaxes the zero penalty. Both default 0.0 ⇒ no effect,
    # so the loss is byte-identical when pi_tech is absent.
    #   pi_soft_target_zero : blend strength for nudging a suspicious observed
    #                         zero's target toward "nonzero".
    #   pi_downweight_zero  : how much to shrink the BCE weight of a suspicious
    #                         observed zero (never touches true nonzero cells).
    # See doc/zero_origin_and_capacity_design_2026-06-01.md.
    pi_soft_target_zero: float = 0.0
    pi_downweight_zero: float = 0.0

    # Asymmetric under-prediction penalty (nonzero_exact lever, 2026-06). Adds a
    # pinball term on the expected ordinal level for NONZERO cells: predicting
    # BELOW the true tier (under) is penalized `asym_tau`, ABOVE is penalized
    # `1 - asym_tau`. `asym_tau > 0.5` ⇒ under penalized more, fighting the
    # model's downward magnitude bias. `lambda_under_asym = 0` ⇒ byte-identical.
    # See doc/asymmetric_ordinal_loss_design_2026-06-02.md.
    lambda_under_asym: float = 0.0
    asym_tau: float = 0.7

    # Numerical safety.
    eps: float = 1e-8

    def __post_init__(self) -> None:
        # When the asymmetric term is active, asym_tau must keep under-prediction
        # penalized MORE than over (tau>0.5) and stay finite (tau<1.0).
        if float(self.lambda_under_asym) > 0.0 and not (0.5 < float(self.asym_tau) < 1.0):
            raise ValueError(
                f"asym_tau must be in (0.5, 1.0) when lambda_under_asym>0 "
                f"(got asym_tau={self.asym_tau}); tau<=0.5 would not penalize "
                f"under-prediction more and tau>=1.0 is degenerate."
            )


@dataclass
class HierarchicalOrdinalOutput:
    """Returned by `compute_hierarchical_ordinal_loss`."""

    # Unweighted component losses (for logging / inspection).
    loss_zero: torch.Tensor
    loss_group: torch.Tensor
    loss_emd: torch.Tensor

    # Weighted sum that is meant to be added to `total`.
    weighted_total: torch.Tensor

    # Diagnostics. All scalars; safe to log per step.
    loss_within_group: torch.Tensor       # L2.5 unweighted within-group binary CE
    n_nonzero_pairs: int                 # count of (cell, gene) pairs with y_true > 0
    pred_zero_frac: torch.Tensor          # mean of p_zero across (cell, gene)
    pred_group_dist: torch.Tensor         # [n_groups] predicted group mean (conditional)
    true_group_dist: torch.Tensor         # [n_groups] empirical group fraction
    ramp_scale_used: float                # ramp_scale this call applied

    # Asymmetric under-prediction term + its monitoring diagnostics.
    # None when the term is disabled (lambda_under_asym == 0).
    loss_under_asym: Optional[torch.Tensor] = None
    asym_diag: Optional[dict] = None


# ======================================================================
# L1 — binary zero / nonzero head with class weighting
# ======================================================================
def compute_zero_nonzero_loss(
    probs: torch.Tensor,
    y_true: torch.Tensor,
    *,
    zero_pos_weight: float,
    eps: float,
    pi_tech: Optional[torch.Tensor] = None,
    pi_soft_target: float = 0.0,
    pi_downweight: float = 0.0,
) -> torch.Tensor:
    """
    Weighted BCE between predicted P(bin > 0) and the binary label (y_true > 0).

    Parameters
    ----------
    probs : [B, G, K]
        Per-(cell, gene) bin probabilities; bin 0 column is the zero class.
    y_true : [B, G] long
        True bin indices in {0, ..., K-1}.
    zero_pos_weight : float
        Per-sample weight applied to nonzero examples. Counters the heavy
        zero-class imbalance (~80%).
    pi_tech : optional [B, G]
        Per-(cell, gene) probability that an observed zero is a technical
        dropout (only meaningful where ``y_true == 0``; should be 0 elsewhere).
        When ``None`` this function is byte-identical to the original loss.
    pi_soft_target : float
        Blend strength for nudging a suspicious observed zero's target toward
        "nonzero": ``target = y_is_nz + pi_soft_target * pi_tech`` on zeros.
    pi_downweight : float
        Multiplicatively shrinks the BCE weight of suspicious observed zeros by
        ``(1 - pi_downweight * pi_tech)``. Never touches true nonzero cells.

    Returns
    -------
    scalar tensor
        Mean weighted BCE over all (cell, gene) pairs.
    """
    if probs.ndim != 3:
        raise ValueError(f"probs must be [B, G, K], got {tuple(probs.shape)}.")
    if y_true.shape != probs.shape[:2]:
        raise ValueError(
            f"y_true must be [B, G] matching probs prefix; "
            f"got y_true={tuple(y_true.shape)}, probs prefix={tuple(probs.shape[:2])}."
        )

    p_zero = probs[..., 0].float()
    p_nonzero = (1.0 - p_zero).float()

    y_is_nz = (y_true > 0).to(dtype=torch.float32)
    weight = torch.where(
        y_is_nz > 0,
        torch.full_like(y_is_nz, float(zero_pos_weight)),
        torch.full_like(y_is_nz, 1.0),
    )

    # π_tech soft-zero: relax the hard "observed zero == off" supervision on
    # zeros the neighborhood/depth model flags as likely technical dropout.
    # Both knobs default to 0 ⇒ target == y_is_nz and weight unchanged.
    target = y_is_nz
    if pi_tech is not None:
        pi = pi_tech.to(dtype=torch.float32)
        if pi.shape != y_is_nz.shape:
            raise ValueError(
                f"pi_tech must be [B, G] matching y_true; got {tuple(pi.shape)}."
            )
        if float(pi_downweight) > 0.0:
            zero_factor = (1.0 - float(pi_downweight) * pi).clamp_min(0.0)
            # only shrink weight on observed zeros; nonzero cells keep their weight
            weight = torch.where(y_is_nz > 0, weight, weight * zero_factor)
        if float(pi_soft_target) > 0.0:
            target = (y_is_nz + float(pi_soft_target) * pi * (1.0 - y_is_nz)).clamp(0.0, 1.0)

    # Stable binary CE on log-odds P(nonzero) / P(zero). The probability-space
    # form log(p) + log(1-p) is unsafe when p saturates to exactly 0 or 1:
    # eps=1e-8 is below fp32 resolution around 1.0, so 1-eps rounds to 1.0.
    safe_eps = max(float(eps), 1e-8)
    p_nz_safe = torch.nan_to_num(
        p_nonzero, nan=0.0, posinf=1.0, neginf=0.0
    ).clamp_min(safe_eps)
    p_z_safe = torch.nan_to_num(
        p_zero, nan=0.0, posinf=1.0, neginf=0.0
    ).clamp_min(safe_eps)
    logits = (torch.log(p_nz_safe) - torch.log(p_z_safe)).clamp(min=-30.0, max=30.0)
    bce_per = F.binary_cross_entropy_with_logits(logits, target, reduction="none")
    return (weight * bce_per).sum() / weight.sum().clamp_min(safe_eps)


# ======================================================================
# L2 — 3-way conditional CE on nonzero cells
# ======================================================================
def _build_group_index(
    y_true: torch.Tensor,
    group_definition,
) -> torch.Tensor:
    """
    Map y_true (in {0..K-1}) to group index (in {0..n_groups-1}) or -1 if
    not in any group (e.g., y_true == 0).
    """
    g_idx = torch.full_like(y_true, -1)
    for gi, bins_in_group in enumerate(group_definition):
        for b in bins_in_group:
            g_idx = torch.where(y_true == int(b), torch.full_like(y_true, gi), g_idx)
    return g_idx


def _build_group_class_weights(
    g_idx_flat: torch.Tensor,
    n_groups: int,
    mode: str,
    cap: float,
    eps: float,
    explicit: Optional[tuple],
) -> torch.Tensor:
    """
    Build per-group weights for the L2 loss.

    If `explicit` is given (tuple of n_groups floats) it is used directly.
    Otherwise the weights are derived from the empirical group frequency
    within the current nonzero subset:

        sqrt_inverse:  w_g = sqrt( N / (n_groups * n_g) )
        inverse:       w_g = N / (n_groups * n_g)

    Weights are clamped to `cap` and rescaled so the expected weight under
    the observed distribution equals 1 (gradient magnitude preservation).
    """
    device = g_idx_flat.device
    dtype = torch.float32

    if explicit is not None:
        w = torch.tensor(list(explicit), dtype=dtype, device=device)
        if w.numel() != n_groups:
            raise ValueError(
                f"group_class_weights has {w.numel()} entries; expected {n_groups}."
            )
        return w

    if mode == "uniform":
        return torch.ones(n_groups, dtype=dtype, device=device)

    counts = torch.bincount(g_idx_flat, minlength=n_groups).to(dtype=dtype)
    total = counts.sum().clamp_min(float(eps))
    inv = total / (float(n_groups) * counts.clamp_min(float(eps)))
    if mode == "sqrt_inverse":
        w = torch.sqrt(inv)
    elif mode == "inverse":
        w = inv
    else:
        raise ValueError(
            f"group_weight_mode must be 'uniform' | 'sqrt_inverse' | 'inverse', "
            f"got {mode!r}."
        )
    if cap is not None and float(cap) > 0:
        w = w.clamp(max=float(cap))
    # Rescale so expected weight ≈ 1 under the observed distribution.
    expected = (w * (counts / total)).sum()
    w = w / expected.clamp_min(float(eps))
    return w


def compute_group_loss(
    probs: torch.Tensor,
    y_true: torch.Tensor,
    *,
    group_definition,
    eps: float,
    group_weight_mode: str = "uniform",
    group_weight_cap: float = 4.0,
    group_class_weights: Optional[tuple] = None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, int]:
    """
    3-way cross-entropy on (low / mid / high) conditional on nonzero.

    Implementation detail (critical): predicted group probabilities are
    *conditional* on bin > 0:

        p_group_g = (sum_{b in group_g} p_b) / (1 - p_0)

    This factorisation guarantees L2 does NOT recover the zero-shortcut
    gradient through the p_0 mass — L1 owns that.

    Only (cell, gene) pairs with y_true > 0 contribute to the loss
    (filtering both inputs and targets after building conditional probs).

    Returns
    -------
    loss : scalar tensor
    pred_group_dist : [n_groups] mean predicted group probability on
        the nonzero subset
    true_group_dist : [n_groups] empirical group fraction on the nonzero
        subset
    n_nz : number of (cell, gene) pairs used (with y_true > 0)
    """
    n_groups = len(group_definition)
    if n_groups < 2:
        raise ValueError("group_definition must define >= 2 groups.")

    p_zero = probs[..., 0]
    p_nz_total = (1.0 - p_zero).clamp_min(float(eps))                # [B, G]

    group_probs = []
    for bins_in_group in group_definition:
        if len(bins_in_group) == 0:
            raise ValueError("group_definition contains empty group.")
        s = torch.zeros_like(p_zero)
        for b in bins_in_group:
            s = s + probs[..., int(b)]
        group_probs.append(s / p_nz_total)
    p_group = torch.stack(group_probs, dim=-1)                       # [B, G, n_groups]
    p_group = p_group.clamp(min=float(eps), max=1.0)

    g_idx = _build_group_index(y_true, group_definition)             # [B, G]
    nz_mask = g_idx >= 0
    n_nz = int(nz_mask.sum().item())

    if n_nz == 0:
        zero = torch.zeros((), dtype=probs.dtype, device=probs.device)
        dist_zero = torch.zeros(n_groups, dtype=probs.dtype, device=probs.device)
        return zero, dist_zero, dist_zero, 0

    p_group_flat = p_group[nz_mask]                                  # [N_nz, n_groups]
    g_idx_flat = g_idx[nz_mask].long()                               # [N_nz]

    # Build group class weights (uniform / sqrt_inverse / inverse / explicit).
    weights = _build_group_class_weights(
        g_idx_flat=g_idx_flat,
        n_groups=n_groups,
        mode=group_weight_mode,
        cap=group_weight_cap,
        eps=eps,
        explicit=group_class_weights,
    ).to(dtype=p_group_flat.dtype)

    # NLL on log-conditional-prob with per-group class weights.
    log_p = torch.log(p_group_flat)
    loss = F.nll_loss(log_p, g_idx_flat, weight=weights)

    pred_group_dist = p_group_flat.mean(dim=0).detach()
    true_group_dist = torch.zeros(
        n_groups, dtype=probs.dtype, device=probs.device
    )
    for gi in range(n_groups):
        true_group_dist[gi] = (g_idx_flat == gi).to(dtype=probs.dtype).mean()
    return loss, pred_group_dist, true_group_dist.detach(), n_nz


# ======================================================================
# L3 — Earth-mover's distance on the 7-bin cumulative
# ======================================================================
def compute_emd_loss(
    probs: torch.Tensor,
    y_true: torch.Tensor,
) -> torch.Tensor:
    """
    EMD (Wasserstein-1) loss on the cumulative distribution.

    For each (cell, gene),
        cum_pred[k] = sum_{j <= k} probs[..., j]
        cum_true[k] = sum_{j <= k} 1[j == y_true]            (step function)

    EMD = mean_{cells, genes} sum_k (cum_pred[k] - cum_true[k])^2.

    This penalises "wrong by m" with weight ~ m^2 along the ordinal axis,
    so adjacent-bin errors cost much less than far-bin errors. The model
    is rewarded for getting *close* even when not exact.
    """
    n_bins = probs.shape[-1]
    y_oh = F.one_hot(y_true.long(), num_classes=n_bins).to(dtype=probs.dtype)
    cum_pred = torch.cumsum(probs, dim=-1)
    cum_true = torch.cumsum(y_oh, dim=-1)
    diff_sq = (cum_pred - cum_true).pow(2).sum(dim=-1)
    return diff_sq.mean()


# ======================================================================
# L2.5 — within-group binary CE (bin b_lo vs b_hi inside each 2-bin group)
# ======================================================================
def compute_within_group_loss(
    probs: torch.Tensor,
    y_true: torch.Tensor,
    *,
    group_definition,
    eps: float,
) -> torch.Tensor:
    """
    For each group of exactly 2 bins [b_lo, b_hi], compute binary CE between
    p(b_lo | in_group) and the indicator (y_true == b_lo), restricted to
    (cell, gene) pairs where y_true ∈ {b_lo, b_hi}.

    This targets the 1-bin-off error: the model knows the right group (L2)
    but cannot distinguish the lower from upper bin within the group.

    Groups with ≠ 2 bins are skipped. Returns the mean loss over active groups.
    """
    active_losses = []
    for bins_in_group in group_definition:
        if len(bins_in_group) != 2:
            continue
        b_lo, b_hi = int(bins_in_group[0]), int(bins_in_group[1])
        mask = (y_true == b_lo) | (y_true == b_hi)
        n = int(mask.sum().item())
        if n == 0:
            continue
        # Stable conditional binary CE. The previous probability-space
        # formula clamped with max=1-eps; with eps=1e-8 this rounds back
        # to exactly 1.0 in fp32/bfloat16, so log(1-p) becomes log(0)
        # and yields Inf/NaN on saturated predictions. Computing the same
        # conditional CE as BCEWithLogits(log(p_lo) - log(p_hi)) avoids
        # the 0 * -inf path entirely.
        p_lo_raw = probs[mask][..., b_lo].float()
        p_hi_raw = probs[mask][..., b_hi].float()
        safe_eps = max(float(eps), 1e-8)
        p_lo_safe = torch.nan_to_num(
            p_lo_raw, nan=0.0, posinf=1.0, neginf=0.0
        ).clamp_min(safe_eps)
        p_hi_safe = torch.nan_to_num(
            p_hi_raw, nan=0.0, posinf=1.0, neginf=0.0
        ).clamp_min(safe_eps)
        logits = (torch.log(p_lo_safe) - torch.log(p_hi_safe)).clamp(
            min=-30.0, max=30.0
        )
        y_bin = (y_true[mask] == b_lo).to(dtype=logits.dtype)
        bce = F.binary_cross_entropy_with_logits(logits, y_bin, reduction="mean")
        if torch.isfinite(bce):
            active_losses.append(bce.to(dtype=probs.dtype))

    if not active_losses:
        return torch.zeros((), dtype=probs.dtype, device=probs.device)
    return torch.stack(active_losses).mean()


# ======================================================================
# Orchestrator
# ======================================================================
def compute_hierarchical_ordinal_loss(
    probs: torch.Tensor,
    y_true: torch.Tensor,
    config: HierarchicalOrdinalConfig,
    *,
    epoch: Optional[int] = None,
    pi_tech: Optional[torch.Tensor] = None,
) -> HierarchicalOrdinalOutput:
    """
    Compute the full hierarchical loss and return all unweighted components
    plus the weighted sum.

    `epoch` enables a linear warm-up over the first `config.ramp_epochs`
    epochs (1-indexed in spirit; passing None disables ramp).

    `pi_tech` ([B, G], optional) is the per-(cell,gene) technical-dropout
    probability for the π_tech soft-zero relaxation. It is forwarded to the L1
    (zero/nonzero) loss only — L2/L3/L2.5 are unchanged, since π_tech is a
    zero-vs-nonzero concept. ``None`` ⇒ byte-identical to the original loss.
    """
    if not config.enabled:
        zero = torch.zeros((), dtype=probs.dtype, device=probs.device)
        dist = torch.zeros(
            len(config.group_definition), dtype=probs.dtype, device=probs.device
        )
        return HierarchicalOrdinalOutput(
            loss_zero=zero,
            loss_group=zero,
            loss_emd=zero,
            loss_within_group=zero,
            weighted_total=zero,
            n_nonzero_pairs=0,
            pred_zero_frac=zero,
            pred_group_dist=dist,
            true_group_dist=dist,
            ramp_scale_used=0.0,
        )

    if probs.ndim != 3:
        raise ValueError(f"probs must be [B, G, K], got {tuple(probs.shape)}.")
    if y_true.shape != probs.shape[:2]:
        raise ValueError(
            "y_true must have shape matching probs prefix; "
            f"got y_true={tuple(y_true.shape)} probs={tuple(probs.shape)}."
        )

    if config.ramp_epochs > 0 and epoch is not None:
        ramp_scale = min(1.0, float(max(epoch, 0)) / float(config.ramp_epochs))
    else:
        ramp_scale = 1.0

    L1 = compute_zero_nonzero_loss(
        probs=probs,
        y_true=y_true,
        zero_pos_weight=float(config.zero_pos_weight),
        eps=float(config.eps),
        pi_tech=pi_tech,
        pi_soft_target=float(config.pi_soft_target_zero),
        pi_downweight=float(config.pi_downweight_zero),
    )
    L2, pred_g, true_g, n_nz = compute_group_loss(
        probs=probs,
        y_true=y_true,
        group_definition=config.group_definition,
        eps=float(config.eps),
        group_weight_mode=str(config.group_weight_mode),
        group_weight_cap=float(config.group_weight_cap),
        group_class_weights=config.group_class_weights,
    )
    L3 = compute_emd_loss(probs=probs, y_true=y_true)
    # L2.5 within-group conditional binary CE is numerically delicate:
    # `p_lo / (p_lo + p_hi)` can lose precision in bfloat16 when the
    # denominator is very small, producing NaN/Inf that propagates into the
    # total loss and (under DDP) into the L2.5 run's NCCL deadlock. We force
    # fp32 here regardless of the surrounding autocast context.
    if probs.is_cuda:
        with torch.amp.autocast(device_type="cuda", enabled=False):
            L25 = compute_within_group_loss(
                probs=probs.float(),
                y_true=y_true,
                group_definition=config.group_definition,
                eps=float(config.eps),
            )
    else:
        L25 = compute_within_group_loss(
            probs=probs.float() if probs.dtype != torch.float32 else probs,
            y_true=y_true,
            group_definition=config.group_definition,
            eps=float(config.eps),
        )

    weighted_total = (
        ramp_scale * float(config.lambda_zero) * L1
        + ramp_scale * float(config.lambda_group) * L2
        + ramp_scale * float(config.lambda_emd) * L3
        + ramp_scale * float(config.lambda_within_group) * L25
    )

    # ------------------------------------------------------------------
    # Asymmetric under-prediction penalty (nonzero_exact lever).
    # Pinball loss on the expected ordinal level for NONZERO cells; tau>0.5
    # penalizes predicting BELOW the true tier (the model's downward bias) more
    # than above. lambda_under_asym == 0 ⇒ no compute, byte-identical.
    # ------------------------------------------------------------------
    loss_under_asym = None
    asym_diag: Optional[dict] = None
    if float(config.lambda_under_asym) > 0.0:
        K = probs.shape[-1]
        levels = torch.arange(K, device=probs.device, dtype=probs.dtype)
        exp_level = (probs * levels).sum(dim=-1)              # [B, G] predicted level
        nz = (y_true > 0).to(exp_level.dtype)
        diff = y_true.to(probs.dtype) - exp_level             # >0 = under-predicted
        tau = float(config.asym_tau)
        pinball = torch.maximum(tau * diff, (tau - 1.0) * diff)
        # masked mean over nonzero cells (no bool()/.item() → no GPU→CPU sync).
        n_nz_f = nz.sum().clamp_min(1.0)
        loss_under_asym = (pinball * nz).sum() / n_nz_f
        weighted_total = (
            weighted_total
            + ramp_scale * float(config.lambda_under_asym) * loss_under_asym
        )
        # Required monitoring (this term can over-correct into zero under-calling,
        # so we watch the gap, high-tier levels, and zero rates). Built as DETACHED
        # tensors — NO float()/.item()/bool() in the loss hot path (each would force
        # a GPU→CPU sync every step). The trainer converts them in ONE batched sync.
        with torch.no_grad():
            exp_d = exp_level.detach()
            pred_bin = probs.detach().argmax(dim=-1)
            nan = exp_d.new_tensor(float("nan"))

            def _mean_exp(tier: int) -> torch.Tensor:
                m = (y_true == tier).to(exp_d.dtype)
                c = m.sum()
                return torch.where(c > 0, (exp_d * m).sum() / c.clamp_min(1.0), nan)

            pred_zero = (pred_bin == 0).to(exp_d.dtype)
            asym_diag = {
                "metric/asym_nonzero_exp_gap": ((y_true.to(exp_d.dtype) - exp_d) * nz).sum() / n_nz_f,
                "metric/asym_mean_exp_true3": _mean_exp(3),
                "metric/asym_mean_exp_true4": _mean_exp(4),
                "metric/asym_nonzero_to_zero_leak": (pred_zero * nz).sum() / n_nz_f,
                "metric/asym_zero_pred": pred_zero.mean(),
                "metric/asym_zero_true": (y_true == 0).to(exp_d.dtype).mean(),
            }

    pred_zero_frac = probs[..., 0].mean().detach()

    return HierarchicalOrdinalOutput(
        loss_zero=L1.detach(),
        loss_group=L2.detach(),
        loss_emd=L3.detach(),
        loss_within_group=L25.detach(),
        weighted_total=weighted_total,
        n_nonzero_pairs=n_nz,
        pred_zero_frac=pred_zero_frac,
        pred_group_dist=pred_g,
        true_group_dist=true_g,
        ramp_scale_used=float(ramp_scale),
        loss_under_asym=(loss_under_asym.detach() if loss_under_asym is not None else None),
        asym_diag=asym_diag,
    )
