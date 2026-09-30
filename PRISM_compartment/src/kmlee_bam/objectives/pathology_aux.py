"""v20 multi-pathology auxiliary head (donor-balanced, aggregate supervision).

Why this exists
---------------
v19 cleaned celltype/sex out of ``z`` but ``z`` collapsed to ~1 effective
dimension (participation_ratio 1.1/32): the reconstruction objective only rewards
the dominant AD-severity axis, so the other 31 latent directions carry no signal.
Whitening alone re-spreads ``z`` but fills the empty axes with *noise*, not biology.

SEA-AD actually carries several partially-independent disease axes (measured at the
donor level: Braak↔Thal r=0.68 [AD-core], but LATE / Lewy / Microinfarct are largely
independent, r<0.4). To make ``z`` rich we must give it an explicit, supervised reason
to represent those axes — that is what this head does.

Pseudo-replication guard (the careful part)
-------------------------------------------
Pathology labels are **donor-level**, not cell-level. A per-cell BCE over ~756k cells
has an effective N of ~89 (the donors), so the model would learn "which donor is this"
(re-leaking donor/batch identity into ``z`` — exactly the nuisance we just removed)
rather than a reproducible cell-state -> pathology relationship. So we **aggregate**:

    z_perp --(EMA bank, per (donor,celltype))--> stable group mean --> pathology head

and the loss is averaged over **groups** (each (donor,celltype) weighted equally — a
10k-cell donor does not dominate a 50-cell donor).

Sparsity / differentiability (mirrors :class:`EMAAlignmentState`)
----------------------------------------------------------------
A 128-cell batch over 89x24 groups is almost all singletons, so a per-batch group mean
is meaningless. An EMA bank accumulates a stable per-group mean over steps. The head
sees the stable detached EMA group mean (a real aggregate, never a single cell), while a
STRAIGHT-THROUGH estimator routes the FULL gradient to the encoder through this batch's
cells:

    value(blended_g) = ema_mean_g.detach()          # head reads the stable aggregate
    grad (blended_g) = d(batch_mean_g) / d(z)        # full-strength signal to z

(An earlier `decay*ema + (1-decay)*batch` blend leaked only a ~(1-decay) fraction of the
gradient to z — far too weak to re-spread the latent; the straight-through form makes
lambda_path_aux the single clean knob for how hard z is shaped. See ``forward``.)

Only groups that already have EMA history are supervised (cold groups just update the
bank), so the head never trains on a raw single cell. The buffers remain per-rank and use
no separate backward, but integrated Phase-I/II runs persist the rank-local bank in the
checkpoint so an epoch-boundary continuation does not discard aggregate history.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F


# Default axes: 4 partially-independent disease dimensions (donor-level measured).
# AD-core is split into Braak(tau)/Thal(amyloid) so the head *can* find the amyloid-tau
# dissociation; LATE(TDP-43) and Lewy are the largely-independent axes that actually add
# new rank. Microinfarct (7/89 donors) is dropped — too few positives to validate.
DEFAULT_PATHOLOGY_AXES: tuple = (
    "axis_braak_high",
    "axis_thal_high",
    "axis_late_positive",
    "axis_lewy_positive",
)


@dataclass  # NOT frozen: runner_base._merge_dataclass populates it via setattr
class PathologyAuxConfig:
    enabled: bool = False
    lambda_path_aux: float = 0.0
    lambda_axis_specificity: float = 0.0  # #2: penalize pred-corr EXCEEDING label-corr
    lambda_resid_axis: float = 0.0  # v23: predict each axis's deviation from overall
                                    # severity (axis - across-axis mean) -> forces
                                    # axis-SPECIFIC signal, blocks 1-D severity collapse
    lambda_interaction_tie: float = 0.0  # v23 (B): cosine-align the decoder's free
                                         # interaction coords to THIS head's supervised
                                         # axis directions -> interaction products become
                                         # named-pathology, curvature lands in disease dirs
    lambda_pathology_tie: float = 0.0    # v26/v27: cosine-align the decoder's pathology_head
                                         # (the named-axis projection p driving the celltype-
                                         # conditioned + pairwise patterns) to THIS head's
                                         # supervised axis directions -> p is pathology-tied,
                                         # not free coords (the v24 failure mode)
    lambda_gate_sparse: float = 0.0      # v35: L1 penalty on the decoder disease-gate per-axis
                                         # strengths s_k (decoder.use_disease_gate) => the model
                                         # SELF-SELECTS which axes curve; dead axes shrink to 0.
                                         # Read by reference_anchored_trainer. See
                                         # doc/v35_disease_gate_design.md.
    axes: Sequence[str] = DEFAULT_PATHOLOGY_AXES
    ema_decay: float = 0.99
    hidden_dim: int = 0  # 0 => linear head
    per_axis_weight: Optional[Sequence[float]] = None
    eps: float = 1e-6


@dataclass
class PathologyAuxOutput:
    loss: torch.Tensor
    per_axis_loss: Dict[str, float] = field(default_factory=dict)
    per_axis_pos_frac: Dict[str, float] = field(default_factory=dict)
    n_groups_total: int = 0
    n_groups_supervised: int = 0
    n_axes_active: int = 0
    specificity_loss: float = 0.0
    residual_frac: float = 0.0
    resid_axis_loss: float = 0.0


class PathologyAuxHead(nn.Module):
    """Donor-balanced multi-pathology head on (donor x celltype)-aggregated z_perp."""

    def __init__(
        self,
        *,
        d_z: int,
        n_donor: int,
        n_celltype: int,
        config: PathologyAuxConfig,
    ) -> None:
        super().__init__()
        self.cfg = config
        self.d_z = int(d_z)
        self.n_celltype = int(n_celltype)
        self.n_group = int(n_donor) * int(n_celltype)
        self.axes: List[str] = list(config.axes)
        self.n_axes = len(self.axes)
        self.decay = float(config.ema_decay)
        if not (0.0 <= self.decay < 1.0):
            raise ValueError(f"ema_decay must be in [0,1), got {self.decay}.")
        if self.n_axes == 0:
            raise ValueError("PathologyAuxHead needs at least one axis.")

        if int(config.hidden_dim) > 0:
            self.head: nn.Module = nn.Sequential(
                nn.Linear(self.d_z, int(config.hidden_dim)),
                nn.GELU(),
                nn.Linear(int(config.hidden_dim), self.n_axes),
            )
        else:
            self.head = nn.Linear(self.d_z, self.n_axes)

        # v23 residualized-axis head (default OFF => no extra params, byte-identical).
        # Predicts each axis's deviation from overall severity (precomputed donor-level
        # targets that are train-fit + standardized + masked); masked MSE on (donor x
        # celltype) aggregates -> forces axis-SPECIFIC signal, blocks 1-D severity collapse.
        self.lambda_resid_axis = float(getattr(config, "lambda_resid_axis", 0.0))
        if self.lambda_resid_axis > 0.0:
            self.resid_head = nn.Linear(self.d_z, self.n_axes)

        w = config.per_axis_weight if config.per_axis_weight is not None else [1.0] * self.n_axes
        if len(w) != self.n_axes:
            raise ValueError(f"per_axis_weight length {len(w)} != n_axes {self.n_axes}.")
        self.register_buffer("axis_weight", torch.tensor(list(w), dtype=torch.float32))

        # EMA bank: first moments per (donor,celltype) group.  Integrated runs
        # persist this state so an epoch-boundary resume is mathematically
        # continuous instead of rewarming donor aggregates from zero.
        self.register_buffer("count_g", torch.zeros(self.n_group), persistent=True)
        self.register_buffer("sum_g", torch.zeros(self.n_group, self.d_z), persistent=True)

    def reset_bank(self) -> None:
        self.count_g.zero_()
        self.sum_g.zero_()

    @torch.no_grad()
    def _update_bank(self, z_det: torch.Tensor, group: torch.Tensor) -> None:
        d = self.decay
        b_count = self.count_g.new_zeros(self.n_group).index_add_(
            0, group, torch.ones_like(group, dtype=self.count_g.dtype)
        )
        b_sum = self.sum_g.new_zeros(self.n_group, self.d_z).index_add_(0, group, z_det)
        present = b_count > 0
        # EMA-update only groups present this batch; absent groups keep their value.
        self.count_g[present] = d * self.count_g[present] + (1.0 - d) * b_count[present]
        self.sum_g[present] = d * self.sum_g[present] + (1.0 - d) * b_sum[present]

    def forward(
        self,
        *,
        z_perp: torch.Tensor,
        donor_id: torch.Tensor,
        celltype_id: torch.Tensor,
        labels: Dict[str, torch.Tensor],
        label_valid: Optional[Dict[str, torch.Tensor]] = None,
        resid_target: Optional[torch.Tensor] = None,   # v23 [B, n_axes] residualized targets
        resid_valid: Optional[torch.Tensor] = None,    # v23 [B, n_axes] validity mask
    ) -> PathologyAuxOutput:
        eps = float(self.cfg.eps)
        dtype = z_perp.dtype
        group = (donor_id.long() * self.n_celltype + celltype_id.long()).clamp_(
            0, self.n_group - 1
        )

        # NB: read the EMA bank BEFORE updating it, so a group's FIRST appearance has
        # no history and is excluded from the loss (never supervise a raw single cell);
        # the bank is updated at the END of this call.

        # 1) per present-group differentiable batch mean + detached EMA mean -> blend
        uniq, inv = torch.unique(group, return_inverse=True)
        G = int(uniq.shape[0])
        ones = z_perp.new_ones(z_perp.shape[0])
        b_sum = z_perp.new_zeros(G, self.d_z).index_add_(0, inv, z_perp)
        b_cnt = z_perp.new_zeros(G).index_add_(0, inv, ones)
        batch_mean = b_sum / b_cnt.clamp_min(1.0).unsqueeze(1)  # [G, d_z], differentiable

        ema_cnt = self.count_g[uniq]
        ema_mean = (self.sum_g[uniq] / ema_cnt.clamp_min(eps).unsqueeze(1)).detach().to(dtype)
        has_hist = (ema_cnt > eps).unsqueeze(1)
        # Straight-through: the head sees the STABLE detached EMA aggregate (value), but
        # z receives the FULL gradient through this batch's cells (value == ema_mean,
        # grad == d(batch_mean)/dz). This decouples "stable target" from "gradient
        # strength" so lambda_path_aux is the single clean knob for how hard z is shaped
        # (the naive d*ema+(1-d)*batch blend leaks only a ~(1-decay) fraction of the
        # gradient to z, far too weak to re-spread the latent). Cold groups (no history)
        # fall back to the plain differentiable batch mean.
        st_blend = ema_mean + (batch_mean - batch_mean.detach())
        blended = torch.where(has_hist, st_blend, batch_mean)

        logits = self.head(blended)  # [G, n_axes]

        # 3) group labels (donor-constant within a group) + masked, group-averaged BCE
        if label_valid is None:
            label_valid = {ax: torch.ones_like(labels[ax], dtype=torch.bool) for ax in self.axes}
        has_hist_g = ema_cnt > eps  # only supervise groups with an aggregated history

        total = z_perp.new_zeros(())
        per_axis_loss: Dict[str, float] = {}
        per_axis_pos: Dict[str, float] = {}
        n_active = 0
        n_sup = 0
        for ai, ax in enumerate(self.axes):
            lab = labels[ax].to(dtype)
            val = label_valid[ax].to(dtype)
            g_lab_sum = z_perp.new_zeros(G).index_add_(0, inv, lab * val)
            g_val = z_perp.new_zeros(G).index_add_(0, inv, val)
            g_use = (g_val > 0) & has_hist_g
            if bool(g_use.any()):
                g_target = (g_lab_sum[g_use] / g_val[g_use].clamp_min(1.0) >= 0.5).to(dtype)
                bce = F.binary_cross_entropy_with_logits(
                    logits[:, ai][g_use], g_target, reduction="mean"
                )
                total = total + float(self.axis_weight[ai]) * bce
                per_axis_loss[ax] = float(bce.detach())
                per_axis_pos[ax] = float(g_target.mean().detach())
                n_active += 1
                n_sup = max(n_sup, int(g_use.sum()))
            else:
                per_axis_loss[ax] = 0.0
                per_axis_pos[ax] = 0.0
        if n_active > 0:
            total = total / n_active

        # v23 RESIDUALIZED-AXIS loss: masked MSE between predicted and precomputed
        # (train-fit, standardized) axis-specific residual targets, on (donor x celltype)
        # AGGREGATES only (never per-cell). Forces deviation-from-severity signal so z
        # cannot collapse to one severity dial. Off (lambda 0) => skipped.
        resid_axis_loss = 0.0
        if self.lambda_resid_axis > 0.0 and resid_target is not None and bool(has_hist_g.any()):
            hh = has_hist_g
            rt = resid_target.to(dtype)                                          # [B, n_axes]
            rv = resid_valid.to(dtype) if resid_valid is not None else torch.ones_like(rt)
            g_rt = z_perp.new_zeros(G, self.n_axes).index_add_(0, inv, rt * rv)
            g_rv = z_perp.new_zeros(G, self.n_axes).index_add_(0, inv, rv)
            g_target_r = (g_rt / g_rv.clamp_min(1.0))[hh].detach()              # [Gh, n_axes]
            g_mask_r = (g_rv[hh] > 0).to(dtype)
            resid_pred = self.resid_head(blended[hh])                           # [Gh, n_axes]
            _se = (resid_pred - g_target_r) ** 2 * g_mask_r
            _rl = _se.sum() / g_mask_r.sum().clamp_min(1.0)
            total = total + self.lambda_resid_axis * _rl
            resid_axis_loss = float(_rl.detach())

        # #2 AXIS-SPECIFICITY: penalize prediction-correlation that EXCEEDS label-
        # correlation. Stops a weak axis (LATE/Lewy) from becoming a SHADOW of the
        # dominant AD-core, WHILE ALLOWING the genuine Braak<->Thal correlation (we
        # penalize only the EXCESS over the label correlation, not correlation itself).
        spec_val = 0.0
        if (
            float(self.cfg.lambda_axis_specificity) > 0.0
            and self.n_axes >= 2
            and int(has_hist_g.sum()) >= 8
        ):
            P = logits[has_hist_g]  # [Gh, K] predicted logits on history groups
            _Lcols = []
            for ax in self.axes:
                _l = labels[ax].to(dtype); _v = label_valid[ax].to(dtype)
                _gl = z_perp.new_zeros(G).index_add_(0, inv, _l * _v)
                _gv = z_perp.new_zeros(G).index_add_(0, inv, _v)
                _Lcols.append((_gl / _gv.clamp_min(1.0))[has_hist_g])
            _Lm = torch.stack(_Lcols, dim=1).detach()  # [Gh, K] group-mean labels

            def _corr(X):
                Xc = X - X.mean(0, keepdim=True)
                Xn = Xc / Xc.std(0, keepdim=True).clamp_min(1e-4)
                return (Xn.t() @ Xn) / max(int(Xn.shape[0]), 1)

            Cp = _corr(P)
            Cl = _corr(_Lm)
            _off = ~torch.eye(self.n_axes, dtype=torch.bool, device=P.device)
            _excess = F.relu(Cp.abs() - Cl.abs())[_off]
            _spec = (_excess * _excess).mean()
            total = total + float(self.cfg.lambda_axis_specificity) * _spec
            spec_val = float(_spec.detach())

        # #1 (light readout): residual fraction — share of the disease aggregate NOT in
        # the span of the K pathology-axis directions. ~0 while z's effective rank <= K;
        # becomes informative only once z carries MORE disease dims than supervised axes
        # (richer data). No params, no loss — a discovery metric (find unexplained disease).
        resid_frac = 0.0
        if isinstance(self.head, nn.Linear) and int(has_hist_g.sum()) >= 4:
            # autocast OFF + float32: under bf16 autocast the matmuls below downcast to
            # bf16, but `_W @ _W.t() + eye(float32)` promotes ONLY A back to float32 →
            # linalg.solve sees A=Float, B=BFloat16 and raises a dtype mismatch (this
            # crashed v22's v7a bootstrap forward). try/except too: #1 is a no-loss
            # discovery metric and must NEVER break the training step.
            try:
                with torch.no_grad(), torch.autocast(device_type=blended.device.type, enabled=False):
                    _Z = blended[has_hist_g].float()
                    _W = self.head.weight.float()  # [K, d_z] axis directions
                    _coef = torch.linalg.solve(
                        _W @ _W.t() + 1e-4 * torch.eye(self.n_axes, device=_W.device, dtype=_W.dtype),
                        _W @ _Z.t(),
                    )
                    _resid = _Z - (_coef.t() @ _W)
                    resid_frac = float(
                        _resid.norm(dim=1).mean() / _Z.norm(dim=1).mean().clamp_min(1e-6)
                    )
            except Exception:
                resid_frac = 0.0

        # update the detached EMA bank LAST (so the loss above used pre-update history)
        self._update_bank(z_perp.detach().to(self.sum_g.dtype), group)

        return PathologyAuxOutput(
            loss=total,
            per_axis_loss=per_axis_loss,
            per_axis_pos_frac=per_axis_pos,
            n_groups_total=G,
            n_groups_supervised=n_sup,
            n_axes_active=n_active,
            specificity_loss=spec_val,
            residual_frac=resid_frac,
            resid_axis_loss=resid_axis_loss,
        )
