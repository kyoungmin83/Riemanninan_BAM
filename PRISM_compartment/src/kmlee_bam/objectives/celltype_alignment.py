from __future__ import annotations

"""Celltype-ALIGNMENT + WHITENING loss on the disease-state latent ``z_perp``.

Motivation
----------
The residual latent ``z_perp = (z_s - mu_c) / sigma_c`` is meant to be a clean
DISEASE fingerprint, but a linear probe shows it predicts CELLTYPE at +0.90 —
cell-type identity has leaked into the disease-state latent. The previous fix (a
GRL aux-adversary on a dedicated reference batch) is mathematically sound but a
DDP nightmare: it needs a separate forward/backward + a manual all-reduce, which
crashes the DDP reducer and then hangs.

This module replaces that adversary with a PURE LOSS on ``z_perp`` that is
DDP-trivial. Both terms are pure functions of the main forward's ``z_perp``
(already computed every step) and are added to the existing total training loss,
so the EXISTING main backward + DDP reducer handle their gradients exactly like
the reconstruction loss. There is intentionally NO ``torch.distributed`` call,
NO separate forward/backward, NO aux batch, NO manual all-reduce, and NO new
optimizer step anywhere in this file — each rank computes the loss on its local
``z_perp`` and DDP averages the gradients automatically.

The two terms
-------------
(1) ``L_align`` — reference-only celltype-INVARIANCE (removes celltype).
    Computed on the REFERENCE / control cells ``R`` (``is_reference == True``):

        m_c    = mean over (R ∩ c) of z_perp                     # [d_z]
        m      = mean over R of z_perp                           # [d_z]
        L_mean = Σ_c (n_c/|R|) · ‖m_c − m‖²                      # between-celltype
        C_c    = cov over (R ∩ c) of z_perp                      # [d_z, d_z]
        C      = cov over R of z_perp                            # [d_z, d_z]
        L_cov  = Σ_c (n_c/|R|) · ‖C_c − C‖²_F / d_z²   (CORAL)
        L_align = L_mean + L_cov          (L_cov only if use_covariance)

    The per-celltype covariance C_c uses n_c ≥ 1 on the EMA path (it is centred by
    the stable detached EMA mean and divided by n, so a single row is meaningful)
    and n_c ≥ 2 on the legacy per-batch path (which centres by the batch's own
    mean and so needs ≥2 rows).

    This is the between-celltype variance of the per-celltype means (plus a CORAL
    second-moment match): driving it to 0 makes the reference cells' ``z_perp``
    celltype-INVARIANT, i.e. erases cell-type identity from the disease latent on
    the pathology-clean cells while leaving disease cells free to keep their
    legitimate cell-type-specific disease direction.

(2) ``L_white`` — global decorrelation (whitening) of ``z_perp``.
    On all cells (or R if ``whiten_all_cells=False``):

        C       = cov of z_perp over those cells                 # [d_z, d_z]
        L_white = ‖C − I‖²_F / d_z²

    This is the "full-covariance" idea expressed as a LOSS (not a prior): it
    pushes the latent toward an isotropic, decorrelated coordinate system.

EMA-accumulated covariance (the per-batch-sparsity fix)
-------------------------------------------------------
With reference cells sampled ~7 per batch over ~5 celltypes, the per-celltype
covariance ``C_c`` almost never has the n_c ≥ 2 cells it needs, so the per-batch
``L_cov`` above is DEAD (~0.0002, almost always skipped). To activate it without
any cross-rank communication, an :class:`EMAAlignmentState` accumulates DETACHED
per-celltype reference statistics across batches (count, sum, sum-of-outer-
products, plus a global aggregate) with an EMA. From those it derives stable
``mu_c``, ``cov_c``, ``mu_global``, ``cov_global``.

The loss must still be DIFFERENTIABLE w.r.t. ``z_perp`` so the gradient reaches
the encoder; a pure detached EMA would carry no gradient. So the differentiable
term BLENDS the stable detached EMA with the CURRENT batch (which IS connected to
``z_perp``):

        mu_c_used  = decay · mu_c_ema.detach()  + (1-decay) · mu_c_batch
        cov_c_used = decay · cov_c_ema.detach() + (1-decay) · cov_c_batch

where ``mu_c_batch`` / ``cov_c_batch`` come from THIS batch's reference cells of
celltype c (differentiable). On the EMA path BOTH the mean and the covariance use
n_c ≥ 1: a celltype with no current cell (n_c == 0) falls back to the detached
EMA for that step (no gradient, never NaN), but a SINGLE reference cell already
contributes a meaningful, stable, differentiable term. The covariance is centred
by the stable detached EMA mean ``mu_c_ema`` (not the 1-cell batch mean) and
divided by n, which is exactly what makes a 1-row batch covariance meaningful and
safe (the rank-1 ``(z − mu_c_ema)(z − mu_c_ema)ᵀ``). (The legacy per-batch path
keeps its n_c ≥ 2 covariance gate, since it centres by the batch's own mean.)
The same EMA aggregate gives the global ``mu_global_used`` / ``cov_global_used``
and the per-celltype weights ``w_c`` (from the EMA counts). The result is a
covariance contribution that is MEANINGFULLY non-zero and decreasing over steps,
unaffected by per-batch sparsity.

A per-SEX variant of the SAME mechanism is available (align per-sex reference-z
mean + cov) when a per-cell ``sex_id`` is supplied; see
:func:`compute_group_alignment_loss`. (The DLPFC+MTG training batch does NOT
currently carry a sex label — see TotalLoss — so the sex term is gated off by
default and ``TotalLoss`` raises a clear error if asked to use it without one.)

Both use unbiased=False (population) covariances and stay differentiable w.r.t.
``z_perp`` — ``z_perp`` is NOT detached, so the gradient flows to the encoder.
Edge cases return a graph-connected ``z_perp.sum() * 0.0`` (never break autograd,
never NaN).
"""

from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn


__all__ = [
    "compute_celltype_alignment_loss",
    "compute_group_alignment_loss",
    "EMAAlignmentState",
]


def _zero_connected(z_perp: torch.Tensor) -> torch.Tensor:
    """A scalar 0 that stays connected to ``z_perp`` in the autograd graph.

    Using ``z_perp.sum() * 0.0`` (rather than ``torch.zeros(())``) keeps the
    parameter that produced ``z_perp`` in the backward graph, so DDP's
    unused-parameter traversal stays stable and the term never breaks the graph.
    The value is exactly 0 and its gradient is exactly 0.
    """
    return z_perp.sum() * 0.0


def _population_cov(x: torch.Tensor) -> torch.Tensor:
    """Population (unbiased=False) covariance of rows of ``x`` ([n, d] -> [d, d]).

    Centers by the per-feature mean. Differentiable w.r.t. ``x``. Caller is
    responsible for ensuring ``x`` has at least 1 row; with a single row the
    covariance is the zero matrix (centered residual is zero), which is the
    correct population covariance for n=1.
    """
    n = x.shape[0]
    mean = x.mean(dim=0, keepdim=True)
    xc = x - mean
    # unbiased=False ⇒ divide by n (not n-1).
    return (xc.transpose(0, 1) @ xc) / float(n)


def _cov_centered_by(x: torch.Tensor, center: torch.Tensor) -> torch.Tensor:
    """Population second-central-moment of rows of ``x`` about an EXTERNAL ``center``.

    ``(1/n) Σ_i (x_i - center)(x_i - center)^T``. ``center`` is the stable
    DETACHED EMA mean (``[d]``), so even a single-row batch yields a meaningful,
    differentiable contribution (the residual is not forced to zero the way
    centering by the batch's own mean would). Differentiable w.r.t. ``x``.
    """
    n = x.shape[0]
    xc = x - center.unsqueeze(0)
    return (xc.transpose(0, 1) @ xc) / float(n)


# ======================================================================
# EMA accumulator for per-(group) reference-z covariance alignment
# ======================================================================
class EMAAlignmentState(nn.Module):
    """DETACHED EMA accumulator of per-group reference-``z_perp`` statistics.

    Holds, per group g (a celltype id, or a sex id for the optional sex term),
    running EMA estimates of:

        count_g  ≈ EMA of the per-batch reference-cell count of group g
        sum_g    ≈ EMA of the per-batch Σ z_perp over (R ∩ g)        # [d_z]
        outer_g  ≈ EMA of the per-batch Σ z_perp z_perpᵀ over (R∩g)  # [d_z, d_z]

    plus a global aggregate over all reference cells (``count_all`` / ``sum_all``
    / ``outer_all``). From these the stable detached estimates are derived:

        mu_g  = sum_g / count_g
        cov_g = outer_g / count_g − mu_g mu_gᵀ            (second-central moment)

    All buffers are DETACHED running statistics — there is never any autograd
    through the stored state. Buffers are registered NON-PERSISTENT so they do
    NOT enter ``state_dict`` (they are recomputable running stats, and a warm
    start from an older checkpoint that lacks them must not fail); they still
    move with ``.to(device)`` / ``.cuda()`` like the module.

    Lazy sizing: ``d_z`` and the group capacity are inferred on the first
    ``update`` from the data, and the group capacity GROWS automatically if a
    larger group id appears later — so the owning ``TotalLoss`` needs neither
    ``d_z`` nor ``n_groups`` at construction time.

    There are intentionally NO ``torch.distributed`` calls here: each rank keeps
    its own local EMA, exactly like a BatchNorm running mean under DDP.
    """

    def __init__(self, decay: float = 0.99, eps: float = 1e-6) -> None:
        super().__init__()
        if not (0.0 <= decay < 1.0):
            raise ValueError(f"EMA decay must be in [0, 1), got {decay}.")
        self.decay = float(decay)
        self.eps = float(eps)
        self._initialized: bool = False
        self._d_z: int = 0
        # The statistic buffers (count_g, sum_g, outer_g, count_all, sum_all,
        # outer_all) are registered lazily in ``_register_group_buffers`` once
        # ``d_z`` / the group count are known from the first batch. They are NOT
        # pre-declared here: assigning ``self.count_g = None`` would create a
        # plain attribute that then collides with ``register_buffer`` ("attribute
        # already exists"). Until initialised they simply do not exist; every
        # accessor guards on ``self._initialized`` first.

    # -- sizing ---------------------------------------------------------
    def _register_group_buffers(self, n_groups: int, d_z: int, ref: torch.Tensor) -> None:
        dtype = torch.float32  # statistics in fp32 for stability under AMP
        device = ref.device
        self.register_buffer(
            "count_g", torch.zeros(n_groups, dtype=dtype, device=device), persistent=False
        )
        self.register_buffer(
            "sum_g", torch.zeros(n_groups, d_z, dtype=dtype, device=device), persistent=False
        )
        self.register_buffer(
            "outer_g",
            torch.zeros(n_groups, d_z, d_z, dtype=dtype, device=device),
            persistent=False,
        )
        self.register_buffer(
            "count_all", torch.zeros((), dtype=dtype, device=device), persistent=False
        )
        self.register_buffer(
            "sum_all", torch.zeros(d_z, dtype=dtype, device=device), persistent=False
        )
        self.register_buffer(
            "outer_all", torch.zeros(d_z, d_z, dtype=dtype, device=device), persistent=False
        )
        self._d_z = int(d_z)
        self._initialized = True

    def _ensure(self, n_groups_needed: int, d_z: int, ref: torch.Tensor) -> None:
        if not self._initialized:
            self._register_group_buffers(max(n_groups_needed, 1), d_z, ref)
            return
        if d_z != self._d_z:
            raise ValueError(
                f"EMAAlignmentState was initialised with d_z={self._d_z} but got d_z={d_z}."
            )
        # Grow the per-group buffers if a larger group id has appeared.
        cur = int(self.count_g.shape[0])
        if n_groups_needed > cur:
            extra = n_groups_needed - cur
            self.count_g = torch.cat(
                [self.count_g, self.count_g.new_zeros(extra)], dim=0
            )
            self.sum_g = torch.cat(
                [self.sum_g, self.sum_g.new_zeros(extra, self._d_z)], dim=0
            )
            self.outer_g = torch.cat(
                [self.outer_g, self.outer_g.new_zeros(extra, self._d_z, self._d_z)], dim=0
            )

    # -- update ---------------------------------------------------------
    @torch.no_grad()
    def update(self, z_ref: torch.Tensor, group_ref: torch.Tensor) -> None:
        """EMA-update the detached statistics from this batch's reference cells.

        ``z_ref`` : [n_ref, d_z] reference-cell latents (DETACHED internally).
        ``group_ref`` : [n_ref] long group ids for those cells.

        Only groups PRESENT in this batch are decayed/updated; absent groups keep
        their previous EMA untouched (so a group missing from one batch does not
        get pulled toward zero). Robust to n_ref == 0 (no-op after sizing).
        """
        z = z_ref.detach().to(torch.float32)
        g = group_ref.long()
        n_ref = int(z.shape[0])
        d_z = int(z.shape[1])
        n_groups_needed = int(g.max().item()) + 1 if n_ref > 0 else 1
        self._ensure(n_groups_needed, d_z, z)
        if n_ref == 0:
            return

        d = self.decay
        # --- global aggregate (over all reference cells this batch) -----
        b_count_all = float(n_ref)
        b_sum_all = z.sum(dim=0)
        b_outer_all = z.transpose(0, 1) @ z
        self.count_all.mul_(d).add_(b_count_all * (1.0 - d))
        self.sum_all.mul_(d).add_(b_sum_all * (1.0 - d))
        self.outer_all.mul_(d).add_(b_outer_all * (1.0 - d))

        # --- per-group: decay+update only the groups present this batch --
        present = torch.unique(g)
        for c in present.tolist():
            mask = g == c
            zc = z[mask]
            nc = float(zc.shape[0])
            b_sum = zc.sum(dim=0)
            b_outer = zc.transpose(0, 1) @ zc
            self.count_g[c] = d * self.count_g[c] + (1.0 - d) * nc
            self.sum_g[c] = d * self.sum_g[c] + (1.0 - d) * b_sum
            self.outer_g[c] = d * self.outer_g[c] + (1.0 - d) * b_outer

    # -- derived stable (detached) estimates ----------------------------
    def mu_global(self) -> torch.Tensor:
        return self.sum_all / self.count_all.clamp_min(self.eps)

    def cov_global(self) -> torch.Tensor:
        mu = self.mu_global()
        second = self.outer_all / self.count_all.clamp_min(self.eps)
        return second - torch.outer(mu, mu)

    def mu_group(self, c: int) -> torch.Tensor:
        return self.sum_g[c] / self.count_g[c].clamp_min(self.eps)

    def cov_group(self, c: int) -> torch.Tensor:
        mu = self.mu_group(c)
        second = self.outer_g[c] / self.count_g[c].clamp_min(self.eps)
        return second - torch.outer(mu, mu)

    def group_weight(self, c: int) -> float:
        denom = float(self.count_all.clamp_min(self.eps).item())
        return float(self.count_g[c].item()) / denom

    def has_group(self, c: int, min_count: float = 1e-3) -> bool:
        if not self._initialized or c >= int(self.count_g.shape[0]):
            return False
        return float(self.count_g[c].item()) >= min_count


# ======================================================================
# Core: group-alignment (mean [+cov]) on reference cells
# ======================================================================
def compute_group_alignment_loss(
    z_perp: torch.Tensor,
    group_id: torch.Tensor,
    ref_mask: Optional[torch.Tensor],
    *,
    use_covariance: bool,
    ema_state: Optional[EMAAlignmentState],
    ema_decay: float,
    update_ema: bool = True,
) -> Tuple[torch.Tensor, torch.Tensor, int, int]:
    """Reference-only group-INVARIANCE term (mean + optional CORAL covariance).

    Shared engine for the celltype term and the optional sex term. ``group_id``
    is the per-cell group label (celltype id or sex id).

    Returns ``(l_mean, l_cov, n_reference, n_groups_in_ref)`` where ``l_mean`` and
    ``l_cov`` are graph-connected scalars (zero-connected when inactive).

    Two modes:
      * ``ema_state is None`` — the ORIGINAL per-batch behaviour: per-group mean
        and (n_c ≥ 2) per-batch covariance, batch-only. Byte-identical to the
        pre-EMA path. The n_c ≥ 2 covariance gate stays here because this path
        centres by the batch's OWN mean (a single row would be forced to zero).
      * ``ema_state is not None`` — the EMA-blended path. The EMA buffers are
        updated (detached) from this batch, then the differentiable loss blends
        the stable detached EMA with this batch's group statistics so gradient
        flows to ``z_perp`` despite sparse per-batch groups (see module docstring).
        Here the covariance uses n_c ≥ 1 (centred by the detached EMA mean), so a
        SINGLETON-reference group still contributes a differentiable cov gradient.
    """
    l_mean_t = _zero_connected(z_perp)
    l_cov_t = _zero_connected(z_perp)

    if ref_mask is None:
        return l_mean_t, l_cov_t, 0, 0
    # Drop cells with a NEGATIVE group id (the missing-label sentinel, e.g.
    # sex_id == -1 for a cell lacking a sex value): they are IGNORED by this term
    # rather than forming a spurious group or wrapping a buffer index. Celltype
    # ids are always >= 0, so this never changes the celltype path.
    valid_mask = ref_mask.bool() & (group_id >= 0)
    if int(valid_mask.sum().item()) == 0:
        return l_mean_t, l_cov_t, 0, 0
    ref_mask = valid_mask

    d_z = z_perp.shape[1]
    z_ref = z_perp[ref_mask]
    g_ref = group_id[ref_mask]
    n_ref = int(z_ref.shape[0])
    classes = torch.unique(g_ref)
    n_groups_in_ref = int(classes.numel())

    # --------------------------- EMA path -----------------------------
    if ema_state is not None:
        # 1) Update the DETACHED running statistics from this batch (skipped in
        #    eval so validation batches never pollute the running stats — the
        #    same convention as BatchNorm's running mean).
        if update_ema:
            ema_state.update(z_ref, g_ref)
        decay = float(ema_decay)

        # 2) Build the differentiable EMA-blended loss.
        #    Need ≥2 groups WITH accumulated mass to have any between-group
        #    structure to remove; otherwise mu_g == mu_global ⇒ exactly zero.
        active = [c for c in classes.tolist() if ema_state.has_group(int(c))]
        if len(active) >= 2:
            mu_glob = ema_state.mu_global().detach().to(z_ref.dtype)  # [d_z]
            cov_glob_used = None
            if use_covariance:
                cov_glob_used = ema_state.cov_global().detach().to(z_ref.dtype)

            total_w = sum(ema_state.group_weight(int(c)) for c in active)
            total_w = total_w if total_w > 0.0 else 1.0

            mean_terms = []
            cov_terms = []
            for c in active:
                c = int(c)
                w_c = ema_state.group_weight(c) / total_w
                mask = g_ref == c
                n_c = int(mask.sum().item())

                mu_c_ema = ema_state.mu_group(c).detach().to(z_ref.dtype)  # [d_z]
                # --- mean blend (fall back to detached EMA if no current cell)
                if n_c >= 1:
                    mu_c_batch = z_ref[mask].mean(dim=0)  # differentiable
                    # STRAIGHT-THROUGH: value = stable detached EMA mean, gradient = FULL
                    # batch (not scaled by (1-decay)). The convex decay-blend's loss optimum
                    # sits at the batch OVERSHOOTING the target by 1/(1-decay) (=100 at decay
                    # 0.99): the encoder could shrink the align loss by INFLATING z -- a
                    # blow-up-prone incentive baked into the loss landscape. Straight-through
                    # makes the loss VALUE the true EMA offset (logged honestly) while the
                    # gradient nudges this batch by the true offset (1x), converging as the
                    # EMA tracks the genuinely aligned distribution.
                    mu_c_used = mu_c_ema + (mu_c_batch - mu_c_batch.detach())
                else:
                    mu_c_used = mu_c_ema
                mean_terms.append(w_c * (mu_c_used - mu_glob).pow(2).sum())

                # --- covariance blend (centre by the STABLE EMA mean) --------
                # EMA path uses n_c >= 1 (not >= 2): _cov_centered_by centers by
                # the EXTERNAL detached EMA mean and divides by float(n), so a
                # SINGLE reference row gives a meaningful, stable, differentiable
                # rank-1 (z - mu_c_ema)(z - mu_c_ema)^T contribution (weight
                # 1-decay). This is the common case (~7 ref cells over ~5
                # celltypes/batch ⇒ singletons), where the legacy n_c>=2 gate kept
                # the cov term dead.
                if use_covariance:
                    cov_c_ema = ema_state.cov_group(c).detach().to(z_ref.dtype)
                    if n_c >= 1:
                        cov_c_batch = _cov_centered_by(z_ref[mask], mu_c_ema)
                        # STRAIGHT-THROUGH (same overshoot fix as the mean above): value =
                        # detached EMA cov, gradient = full batch cov.
                        cov_c_used = cov_c_ema + (cov_c_batch - cov_c_batch.detach())
                    else:
                        cov_c_used = cov_c_ema
                    cov_terms.append(
                        w_c * (cov_c_used - cov_glob_used).pow(2).sum() / float(d_z * d_z)
                    )

            if len(mean_terms) > 0:
                l_mean_t = torch.stack(mean_terms).sum()
            if use_covariance and len(cov_terms) > 0:
                l_cov_t = torch.stack(cov_terms).sum()

        return l_mean_t, l_cov_t, n_ref, n_groups_in_ref

    # ----------------------- per-batch path (legacy) ------------------
    # EDGE CASE: a single group present ⇒ no between-group structure to remove;
    # m_c == m and C_c == C, so both terms are exactly zero.
    if n_groups_in_ref >= 2:
        m = z_ref.mean(dim=0)  # [d_z], global reference mean
        mean_terms = []
        for c in classes:
            m_c_mask = g_ref == c
            n_c = int(m_c_mask.sum().item())
            if n_c < 1:
                continue
            w_c = float(n_c) / float(n_ref)
            m_c = z_ref[m_c_mask].mean(dim=0)  # [d_z]
            mean_terms.append(w_c * (m_c - m).pow(2).sum())
        if len(mean_terms) > 0:
            l_mean_t = torch.stack(mean_terms).sum()

        if use_covariance:
            C = _population_cov(z_ref)  # [d_z, d_z], global reference cov
            cov_terms = []
            for c in classes:
                m_c_mask = g_ref == c
                n_c = int(m_c_mask.sum().item())
                if n_c < 2:
                    # Population cov undefined-ish for n<2; skip per spec.
                    continue
                w_c = float(n_c) / float(n_ref)
                C_c = _population_cov(z_ref[m_c_mask])  # [d_z, d_z]
                cov_terms.append(w_c * (C_c - C).pow(2).sum() / float(d_z * d_z))
            if len(cov_terms) > 0:
                l_cov_t = torch.stack(cov_terms).sum()

    return l_mean_t, l_cov_t, n_ref, n_groups_in_ref


def compute_celltype_alignment_loss(
    z_perp: torch.Tensor,
    celltype_id: torch.Tensor,
    is_reference: Optional[torch.Tensor],
    *,
    use_covariance: bool = True,
    reference_only: bool = True,
    whiten_all_cells: bool = True,
    ema_state: Optional[EMAAlignmentState] = None,
    ema_decay: float = 0.99,
    update_ema: bool = True,
    sex_id: Optional[torch.Tensor] = None,
    sex_ema_state: Optional[EMAAlignmentState] = None,
    use_sex: bool = False,
    sex_all_cells: bool = False,
    whitening_target: Optional[torch.Tensor] = None,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, Dict[str, float]]:
    """Compute the celltype-alignment and whitening losses on ``z_perp``.

    Parameters
    ----------
    z_perp : [B, d_z] float tensor
        The disease-state residual latent. NOT detached: gradients flow to it.
    celltype_id : [B] long tensor
        Per-cell cell-type index.
    is_reference : [B] bool/0-1 tensor or None
        Mask selecting reference / control (pathology-clean) cells. If None, no
        reference cells are available, so ``L_align`` is a graph-connected zero.
    use_covariance : bool, default True
        Include the CORAL second-moment (covariance) match in ``L_align``.
    reference_only : bool, default True
        Restrict ``L_align`` to reference cells (the intended semantics; see the
        module docstring). When False, ALL cells are treated as "R" for the
        alignment term — this would also erase legitimate disease-induced
        celltype structure, so it is exposed only for ablation.
    whiten_all_cells : bool, default True
        Compute ``L_white`` over all cells. When False, use the reference cells
        only (falls back to all cells when no reference mask is available).
    ema_state : EMAAlignmentState or None, default None
        When provided, the per-celltype covariance/mean alignment uses the
        EMA-accumulated path (activates the covariance term despite sparse
        per-batch reference cells). When None, the original per-batch path is
        used — BYTE-IDENTICAL to the pre-EMA behaviour.
    ema_decay : float, default 0.99
        Blend/decay weight for the EMA path (ignored when ``ema_state is None``).
    sex_id : [B] long tensor or None
        Optional per-cell sex index for the OPTIONAL sex-alignment term. When
        ``use_sex`` is True this must be provided.
    sex_ema_state : EMAAlignmentState or None
        EMA accumulator for the sex term (separate from the celltype one).
    use_sex : bool, default False
        Enable the per-sex alignment term. Requires ``sex_id``. The returned
        ``loss_align_sex`` is a graph-connected zero when disabled.
    sex_all_cells : bool, default False
        Align the sex-mean over ALL cells instead of reference-only. Needed
        because the reference / control cells are heavily sex-imbalanced (≈ a
        single gender), which makes a reference-only sex term a near no-op (no
        sex contrast to align). When True the sex term is MEAN-only over all
        cells (the post-hoc-validated direction); the celltype term is
        unaffected and stays reference-only.
    whitening_target : [d_z, d_z] tensor or None, default None
        Detached covariance target for ``L_white``. ``None`` preserves the
        historical full-space identity target. A hard nuisance projector
        should pass ``P = I - Q Q^T`` so removed directions are not asked to
        have unit variance.

    Returns
    -------
    (loss_align, loss_white, loss_align_sex, stats)
        Three scalar tensors and a ``stats`` dict of detached floats for logging:
        ``n_reference``, ``n_celltypes_in_ref``, ``l_mean``, ``l_cov``,
        ``l_white``, ``l_mean_sex``, ``l_cov_sex``, ``n_sex_in_ref``.
    """
    if z_perp.ndim != 2:
        raise ValueError(f"z_perp must have shape [B, d_z], got {tuple(z_perp.shape)}.")
    if celltype_id.ndim != 1:
        raise ValueError(f"celltype_id must have shape [B], got {tuple(celltype_id.shape)}.")
    if celltype_id.shape[0] != z_perp.shape[0]:
        raise ValueError("Batch size mismatch between z_perp and celltype_id.")

    d_z = z_perp.shape[1]
    celltype_id = celltype_id.long()

    # ------------------------------------------------------------------
    # Resolve the reference-cell mask R for the alignment term.
    # ------------------------------------------------------------------
    if reference_only:
        if is_reference is None:
            ref_mask = None
        else:
            ref_mask = is_reference.bool()
            if ref_mask.ndim != 1:
                raise ValueError(
                    f"is_reference must have shape [B], got {tuple(ref_mask.shape)}."
                )
            if ref_mask.shape[0] != z_perp.shape[0]:
                raise ValueError("Batch size mismatch between z_perp and is_reference.")
    else:
        # Ablation: treat every cell as a reference cell for the alignment term.
        ref_mask = torch.ones(z_perp.shape[0], dtype=torch.bool, device=z_perp.device)

    # ------------------------------------------------------------------
    # (1) L_align — reference-only celltype invariance (mean [+cov]).
    # ------------------------------------------------------------------
    l_mean_t, l_cov_t, n_reference, n_celltypes_in_ref = compute_group_alignment_loss(
        z_perp,
        celltype_id,
        ref_mask,
        use_covariance=use_covariance,
        ema_state=ema_state,
        ema_decay=ema_decay,
        update_ema=update_ema,
    )

    loss_align = l_mean_t + l_cov_t if use_covariance else l_mean_t

    # ------------------------------------------------------------------
    # (1b) Optional L_align_sex — reference-only sex invariance.
    # ------------------------------------------------------------------
    l_mean_sex_t = _zero_connected(z_perp)
    l_cov_sex_t = _zero_connected(z_perp)
    n_sex_in_ref = 0
    if use_sex:
        if sex_id is None:
            raise ValueError(
                "use_sex=True but sex_id is None. The training batch does not carry a "
                "sex label yet (see TotalLoss / dataset). Add a per-cell sex id to the "
                "batch before enabling the sex-alignment term."
            )
        if sex_id.ndim != 1 or sex_id.shape[0] != z_perp.shape[0]:
            raise ValueError(
                f"sex_id must have shape [B]={z_perp.shape[0]}, got {tuple(sex_id.shape)}."
            )
        sex_id = sex_id.long()
        # The reference / control cells are heavily sex-IMBALANCED (≈ a single
        # gender), so a reference-only sex term has essentially no sex contrast to
        # align and is a near no-op. `sex_all_cells=True` aligns the sex-mean over
        # ALL cells instead — validated post-hoc: projecting out the all-cells sex
        # mean-difference direction drops z->sex from +0.41 to ~chance while the
        # within-celltype ADNC disease signal is essentially unchanged (+0.113 ->
        # +0.107). All-cells sex alignment is MEAN-only (the validated direction)
        # and uses the per-batch path (both sexes are present every batch); the
        # covariance / EMA paths are intentionally bypassed there.
        if sex_all_cells:
            sex_mask = torch.ones(
                z_perp.shape[0], dtype=torch.bool, device=z_perp.device
            )
            sex_use_cov = False
            sex_ema = None
        else:
            sex_mask = ref_mask
            sex_use_cov = use_covariance
            sex_ema = sex_ema_state
        l_mean_sex_t, l_cov_sex_t, _, n_sex_in_ref = compute_group_alignment_loss(
            z_perp,
            sex_id,
            sex_mask,
            use_covariance=sex_use_cov,
            ema_state=sex_ema,
            ema_decay=ema_decay,
            update_ema=update_ema,
        )
    loss_align_sex = (
        l_mean_sex_t + l_cov_sex_t if use_covariance else l_mean_sex_t
    )

    # ------------------------------------------------------------------
    # (2) L_white — global decorrelation / whitening.
    # ------------------------------------------------------------------
    if whiten_all_cells or ref_mask is None:
        z_white = z_perp
    else:
        # Reference-only whitening (falls back to all cells if R is empty so the
        # term still does something useful and stays graph-connected).
        if int(ref_mask.sum().item()) > 0:
            z_white = z_perp[ref_mask]
        else:
            z_white = z_perp

    if whitening_target is not None:
        if not isinstance(whitening_target, torch.Tensor):
            raise TypeError("whitening_target must be a torch.Tensor or None.")
        if whitening_target.ndim != 2 or tuple(
            whitening_target.shape
        ) != (d_z, d_z):
            raise ValueError(
                "whitening_target must have shape "
                f"[{d_z}, {d_z}], got {tuple(whitening_target.shape)}."
            )
        if not bool(torch.isfinite(whitening_target).all()):
            raise ValueError("whitening_target contains non-finite values.")
        # CUDA autocast can turn the projector audit into bf16 matmuls and
        # ``eigvalsh`` is not supported for every low-precision backend.  This
        # is a fail-closed API check, so always perform it in fp32.
        with torch.autocast(
            device_type=whitening_target.device.type,
            enabled=False,
        ):
            target_for_rank = whitening_target.detach().float()
            symmetry_error = (
                target_for_rank - target_for_rank.t()
            ).abs().max()
            if float(symmetry_error) > 1e-4:
                raise ValueError(
                    "whitening_target must be symmetric; maximum asymmetry="
                    f"{float(symmetry_error):.3e}."
                )
            target_for_rank = 0.5 * (
                target_for_rank + target_for_rank.t()
            )
            idempotence_error = (
                target_for_rank @ target_for_rank - target_for_rank
            ).abs().max()
            if float(idempotence_error) > 1e-3:
                raise ValueError(
                    "whitening_target must be an orthogonal projector; "
                    f"maximum |P^2-P|={float(idempotence_error):.3e}."
                )
            eigenvalues = torch.linalg.eigvalsh(target_for_rank)
            if (
                float(eigenvalues.min()) < -1e-3
                or float(eigenvalues.max()) > 1.001
            ):
                raise ValueError(
                    "whitening_target projector eigenvalues must lie in [0, 1]."
                )
        white_target_rank = float(round(float(torch.trace(target_for_rank))))
    else:
        white_target_rank = float(d_z)
    white_removed_rank = float(d_z) - white_target_rank
    projection_active = white_removed_rank > 0.5

    if z_white.shape[0] <= 1:
        # Cov is degenerate with ≤1 cell ⇒ graph-connected zero.
        loss_white = _zero_connected(z_perp)
        white_active_error = loss_white.detach()
        white_removed_leakage = loss_white.detach()
        white_cross_leakage = loss_white.detach()
    else:
        # The opt-in projector path is evaluated in fp32 even when the outer
        # model forward uses bf16 autocast.  Otherwise the dense P multiplications
        # introduce a small, hardware-dependent gradient discrepancy relative
        # to the historical identity-target objective.
        if whitening_target is None or not projection_active:
            C_all = _population_cov(z_white)  # [d_z, d_z]
        else:
            with torch.autocast(device_type=z_white.device.type, enabled=False):
                work = (
                    z_white.float()
                    if z_white.dtype in (torch.float16, torch.bfloat16)
                    else z_white
                )
                C_all = _population_cov(work)  # [d_z, d_z]
        # Keep the exact block decomposition in the same precision as C.  If
        # these P@C products fall back into outer bf16 autocast, the logged
        # active+removed+2*cross terms no longer reconstruct the optimized loss.
        with torch.autocast(device_type=C_all.device.type, enabled=False):
            eye = torch.eye(d_z, dtype=C_all.dtype, device=C_all.device)
            if whitening_target is None:
                target = eye
            else:
                target = whitening_target.detach().to(
                    device=C_all.device,
                    dtype=C_all.dtype,
                )
                target = 0.5 * (target + target.t())

            complement = eye - target
            normalizer = float(d_z * d_z)
            loss_white = (C_all - target).pow(2).sum() / normalizer
            white_active_error = (
                target @ C_all @ target - target
            ).pow(2).sum().detach() / normalizer
            white_removed_leakage = (
                complement @ C_all @ complement
            ).pow(2).sum().detach() / normalizer
            white_cross_leakage = (
                target @ C_all @ complement
            ).pow(2).sum().detach() / normalizer
            white_identity_target_loss = (
                C_all - eye
            ).pow(2).sum().detach() / normalizer
            white_target_loss_delta = (
                white_identity_target_loss - loss_white.detach()
            )

    if z_white.shape[0] <= 1:
        white_identity_target_loss = loss_white.detach()
        white_target_loss_delta = loss_white.detach()

    stats: Dict[str, float] = {
        "n_reference": float(n_reference),
        "n_celltypes_in_ref": float(n_celltypes_in_ref),
        "l_mean": float(l_mean_t.detach()),
        "l_cov": float(l_cov_t.detach()),
        "l_white": float(loss_white.detach()),
        "l_mean_sex": float(l_mean_sex_t.detach()),
        "l_cov_sex": float(l_cov_sex_t.detach()),
        "n_sex_in_ref": float(n_sex_in_ref),
        "white_active_error": float(white_active_error),
        "white_removed_leakage": float(white_removed_leakage),
        "white_cross_leakage": float(white_cross_leakage),
        "white_target_rank": float(white_target_rank),
        "white_removed_rank": float(white_removed_rank),
        "white_identity_target_loss": float(white_identity_target_loss),
        "white_target_loss_delta": float(white_target_loss_delta),
    }

    return loss_align, loss_white, loss_align_sex, stats
