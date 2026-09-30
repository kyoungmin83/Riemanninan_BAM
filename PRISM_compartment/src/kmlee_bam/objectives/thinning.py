"""
Binomial-thinning consistency + proven technical-zero supervision.

Why
---
Two jobs that anchor the whole zero-origin machinery:

  1. **Depth-invariance (consistency).** The latent-kNN π_tech estimator
     (`objectives/latent_knn_pi_tech.py`) assumes `z_perp` encodes *biology*,
     not measurement depth. If a shallow cell and a deep cell of the same state
     got different `z_perp` (because the shallow one has more zeros), kNN would
     group cells by depth, not biology, and π_tech would be garbage. Thinning
     fixes this: take a cell, binomially down-sample its counts (a *known*
     depth reduction of the *same* cell), and force the model to produce the
     same `z_perp`. → "measurement depth must not change my read of the biology."

  2. **Proven technical zeros (identifiability anchor).** A gene that was
     nonzero in the full cell but thins to 0 is a *proven* technical zero — we
     removed the molecules ourselves. These positions (`technical_zero_mask`)
     give ground-truth labels for "this 0 is technical", which is exactly what
     a fixed π_tech formula otherwise lacks. An optional supervision term uses
     them to teach the model not to be confident-off at dropped positions.

The thinned *view* (down-sampled counts re-binned with the SAME edges/log1p
code) is produced in the dataset (`OrdinalScDataset(return_thinned=True)`), so
this module only needs the tensors. The thinned forward is used **only** for
consistency + (optional) supervision — NOT as a second reconstruction target
(its label is corrupted by construction).

See doc/zero_origin_and_capacity_design_2026-06-01.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class ThinningConfig:
    enabled: bool = False
    # Must match the dataset's `thin_rho` (kept-molecule probability). Kept here
    # only for logging/sanity; the actual thinning happens in the dataset.
    rho: float = 0.5
    lambda_consistency: float = 0.1          # weight on z_perp depth-invariance
    consistency_metric: str = "smooth_l1"    # "smooth_l1" | "mse" | "cosine"
    detach_full: bool = True                 # full-view z is the (teacher) target
    every_n_steps: int = 1                   # apply on every Nth step (cost control)
    lambda_tech_zero_sup: float = 0.0        # optional proven-tech-zero supervision (0 = off)
    eps: float = 1e-6


@dataclass(frozen=True)
class ZeroReliabilityConfig:
    """Stage D: downweight suspicious *observed zeros* in the MAIN reconstruction
    NLL, using `pi_tech` as the technical-dropout probability, so technical zeros
    do not pull the reconstruction (and the latent) toward a hard biological off.
    Off ⇒ reconstruction is byte-identical. Nonzero labels are never touched."""
    enabled: bool = False
    alpha: float = 0.5     # discount strength: zero weight = 1 - alpha * pi_tech
    floor: float = 0.3     # never discount a zero below this weight (keep a signal)


def zero_reliability_weights(
    y_ord: torch.Tensor,
    pi_tech: Optional[torch.Tensor],
    *,
    alpha: float = 0.5,
    floor: float = 0.3,
) -> Optional[torch.Tensor]:
    """Per-(cell, gene) reconstruction-reliability weights in ``[floor, 1]``.

    Suspicious observed zeros (high `pi_tech`) are downweighted; nonzero labels
    keep weight 1.0 (never downweighted). Returns ``None`` when `pi_tech` is None
    so the caller leaves the reconstruction unweighted (byte-identical)."""
    if pi_tech is None:
        return None
    pi = pi_tech.to(dtype=torch.float32)
    rel = (1.0 - float(alpha) * pi).clamp(min=float(floor), max=1.0)
    is_zero = y_ord == 0
    return torch.where(is_zero, rel, torch.ones_like(rel))


def thinning_consistency_loss(
    z_full: torch.Tensor,
    z_thin: torch.Tensor,
    *,
    metric: str = "smooth_l1",
    detach_full: bool = True,
) -> torch.Tensor:
    """Penalise change in `z_perp` between the full and thinned views of the
    same cells. `detach_full` treats the full view as a fixed teacher so the
    gradient only pulls the thinned (degraded) view toward it."""
    target = z_full.detach() if detach_full else z_full
    if metric == "cosine":
        return (1.0 - F.cosine_similarity(z_thin, target, dim=-1)).mean()
    if metric == "mse":
        return F.mse_loss(z_thin, target)
    return F.smooth_l1_loss(z_thin, target)


def technical_zero_mask(
    y_ord: torch.Tensor,
    y_ord_thin: torch.Tensor,
) -> torch.Tensor:
    """[B, G] bool: positions that were nonzero in the full view but thinned to
    0 — i.e. *proven* technical zeros."""
    return (y_ord > 0) & (y_ord_thin == 0)


def technical_zero_supervision_loss(
    probs_thin: torch.Tensor,
    tech_zero_mask: torch.Tensor,
    *,
    eps: float = 1e-6,
) -> torch.Tensor:
    """At proven technical zeros, the thinned-view input shows 0 but we KNOW the
    gene was expressed; encourage the model's thinned-view P(nonzero) to stay
    high there (i.e. `-log p_nonzero`). Returns 0 when the mask is empty."""
    if not bool(tech_zero_mask.any()):
        return probs_thin.new_zeros(())
    p_zero = probs_thin[..., 0].float()
    p_nz = (1.0 - p_zero).clamp_min(float(eps))
    return (-torch.log(p_nz))[tech_zero_mask].mean()


def build_thinned_batch(
    batch: Dict[str, torch.Tensor],
) -> Optional[Dict[str, torch.Tensor]]:
    """Shallow-copy `batch` with the thinned views swapped into the model-input
    keys (`y_ord`, `x_log1p`, `x_gene_scalar`). Returns None if the dataset did
    not provide thinned views (so callers can no-op gracefully)."""
    if "y_ord_thin" not in batch:
        return None
    thinned = dict(batch)
    thinned["y_ord"] = batch["y_ord_thin"]
    if "x_log1p_thin" in batch:
        thinned["x_log1p"] = batch["x_log1p_thin"]
    if "x_gene_scalar_thin" in batch:
        thinned["x_gene_scalar"] = batch["x_gene_scalar_thin"]
    return thinned
