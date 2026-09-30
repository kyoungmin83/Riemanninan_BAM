"""Tech-invariance training: BAM as a supervised technical-OOD detector + z robustness.

Philosophy (extends thinning [[kmlee-bam-thinning]]): a technical corruption (extra dropout /
size-factor gain / ambient contamination) degrades *how we saw* a cell, not *what the cell is*.
So a good latent and a good reliability signal must satisfy:

  (1) z is INVARIANT to tech corruption                       -> L_z_cons   ||z_clean - z_tech||^2
  (2) the encoder's tech-risk head FLAGS the corruption        -> L_bam_tech BCE(clean=0, tech=1)
  (3) the prediction stays consistent with the CLEAN one       -> L_rec_cons (never trust corrupted target)
  (4) a real DISEASE-direction shift does NOT raise the risk   -> disease-safe negative control (risk=0)

The risk head is a cheap SINGLE-FORWARD amortization of the post-hoc MC-epistemic firewall signal
(already validated tech-specific: tech/disease mc_sens ratio ~3-4x, decoder-complementary on
dropout/gain). We distil that signal into a trained scalar instead of paying N-sample MC per step.

This module is pure helpers (augmenter + loss fns); the trainer owns the wiring + lambdas.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import torch
import torch.nn.functional as F


@dataclass
class TechInvarianceConfig:
    """Config for tech-invariance training. All lambdas default 0 => term off even if enabled."""
    enabled: bool = False
    every_n_steps: int = 1
    lambda_z_cons: float = 0.0       # ||z_clean - z_tech||^2  (z robust to tech corruption)
    lambda_bam_tech: float = 0.0     # BCE(risk_clean=0, risk_tech=1)  (supervised tech detector)
    lambda_rec_cons: float = 0.0     # pred(tech) ~ pred(clean)  (never trust corrupted target)
    lambda_disease_ctrl: float = 0.0 # BCE(risk_disease=0)  (disease negative control; needs disease_dir)
    tech_kinds: Tuple[str, ...] = ("dropout", "gain", "ambient")
    dropout_p: float = 0.35
    gain_lo: float = 0.5
    gain_hi: float = 1.8
    ambient_a: float = 0.4
    disease_gamma: float = 1.0
    rec_metric: str = "etier"        # "etier" (MSE on E[tier]) or "kl"
    detach_clean: bool = True
    warmup_steps: int = 0            # linear ramp of all tech-inv lambdas over the first N steps


def _corrupt_counts(x_log1p: torch.Tensor, kind: str, *, dropout_p: float,
                    gain_lo: float, gain_hi: float, ambient_a: float,
                    gen: Optional[torch.Generator]) -> torch.Tensor:
    """Apply a technical corruption in COUNT space; return corrupted x_log1p. Biology-agnostic."""
    c = torch.expm1(x_log1p).clamp_min(0.0)
    dev = c.device
    if kind == "dropout":                       # extra capture loss (stronger than training thinning)
        r = torch.rand(c.shape, generator=gen, device=dev) if gen is not None else torch.rand_like(c)
        c = c.masked_fill((r < dropout_p) & (c > 0), 0.0)
    elif kind == "gain":                         # size-factor / normalization mismatch (multiplicative)
        u = torch.rand(c.shape[0], 1, generator=gen, device=dev) if gen is not None else torch.rand(c.shape[0], 1, device=dev)
        log_g = torch.log(torch.tensor(gain_lo, device=dev)) + u * (torch.log(torch.tensor(gain_hi, device=dev)) - torch.log(torch.tensor(gain_lo, device=dev)))
        c = c * torch.exp(log_g)
    elif kind == "ambient":                      # ambient-RNA / soup contamination
        c = c + ambient_a * c.mean(dim=0, keepdim=True)
    else:
        raise ValueError(f"unknown tech corruption kind: {kind}")
    return torch.log1p(c)


class TechInvarianceAugmenter:
    """Builds technically-corrupted (and optionally disease-shifted) VIEWS of a batch, keeping
    covariates fixed and recomputing x_gene_scalar consistently (so the shift is 'hidden' from the
    known tech covariates -- the unknown-tech scenario the firewall targets)."""

    def __init__(self, scaler_mean: torch.Tensor, scaler_std_eps: torch.Tensor, clip: Optional[float],
                 *, tech_kinds: List[str], dropout_p: float = 0.35, gain_lo: float = 0.5,
                 gain_hi: float = 1.8, ambient_a: float = 0.4,
                 disease_dir: Optional[torch.Tensor] = None) -> None:
        self.scaler_mean = scaler_mean        # [C, G]
        self.scaler_std_eps = scaler_std_eps  # [C, G]
        self.clip = None if clip is None else float(clip)
        self.tech_kinds = list(tech_kinds)
        self.dropout_p = float(dropout_p)
        self.gain_lo = float(gain_lo)
        self.gain_hi = float(gain_hi)
        self.ambient_a = float(ambient_a)
        self.disease_dir = disease_dir        # [C, G] or None (x_log1p-space hi-AD - lo-AD)

    def to(self, device) -> "TechInvarianceAugmenter":
        self.scaler_mean = self.scaler_mean.to(device)
        self.scaler_std_eps = self.scaler_std_eps.to(device)
        if self.disease_dir is not None:
            self.disease_dir = self.disease_dir.to(device)
        return self

    def _recompute_xgs(self, x_log1p: torch.Tensor, celltype_id: torch.Tensor) -> torch.Tensor:
        mean_b = self.scaler_mean[celltype_id]          # [B, G]
        std_b = self.scaler_std_eps[celltype_id]        # [B, G]
        xgs = (x_log1p - mean_b) / std_b
        if self.clip is not None:
            xgs = xgs.clamp(-self.clip, self.clip)
        return xgs

    def _rebuild(self, batch: Dict[str, torch.Tensor], x_log1p_new: torch.Tensor) -> Dict[str, torch.Tensor]:
        view = dict(batch)
        view["x_log1p"] = x_log1p_new
        if "x_gene_scalar" in batch:
            view["x_gene_scalar"] = self._recompute_xgs(x_log1p_new, batch["celltype_id"])
        return view

    def make_tech_view(self, batch: Dict[str, torch.Tensor], *, kind: Optional[str] = None,
                       gen: Optional[torch.Generator] = None) -> Dict[str, torch.Tensor]:
        if kind is None:  # rotate kinds by a cheap per-call random index (variety across steps)
            i = int(torch.randint(0, len(self.tech_kinds), (1,), generator=gen,
                                   device=batch["x_log1p"].device).item())
            kind = self.tech_kinds[i]
        xl = _corrupt_counts(batch["x_log1p"], kind, dropout_p=self.dropout_p, gain_lo=self.gain_lo,
                             gain_hi=self.gain_hi, ambient_a=self.ambient_a, gen=gen)
        return self._rebuild(batch, xl)

    def make_disease_view(self, batch: Dict[str, torch.Tensor], *, gamma: float = 1.0) -> Optional[Dict[str, torch.Tensor]]:
        if self.disease_dir is None:
            return None
        dd = self.disease_dir[batch["celltype_id"]]     # [B, G]
        xl = batch["x_log1p"] + gamma * dd
        return self._rebuild(batch, xl)


# ======================================================================
# Loss terms
# ======================================================================
def z_consistency_loss(z_clean: torch.Tensor, z_view: torch.Tensor, *, detach_clean: bool = True) -> torch.Tensor:
    """z must not move under tech corruption. detach_clean=True => pull the corrupted view's z
    toward the clean z (treat clean as the anchor), matching the thinning convention."""
    target = z_clean.detach() if detach_clean else z_clean
    return F.mse_loss(z_view, target)


def tech_risk_bce(logit_clean: torch.Tensor, logit_tech: torch.Tensor,
                  logit_disease: Optional[torch.Tensor] = None,
                  *, disease_weight: float = 1.0) -> torch.Tensor:
    """Supervised tech detector: clean -> 0, tech-corrupted -> 1, (real disease shift -> 0).
    Logits are [B]. The disease term is the negative control: real biology must NOT raise risk."""
    loss = (F.binary_cross_entropy_with_logits(logit_clean, torch.zeros_like(logit_clean))
            + F.binary_cross_entropy_with_logits(logit_tech, torch.ones_like(logit_tech)))
    if logit_disease is not None:
        loss = loss + disease_weight * F.binary_cross_entropy_with_logits(
            logit_disease, torch.zeros_like(logit_disease))
    return loss


def rec_consistency_loss(probs_clean: torch.Tensor, probs_view: torch.Tensor,
                         *, n_bins: int, metric: str = "etier") -> torch.Tensor:
    """Prediction should match the CLEAN prediction (never the corrupted input's target).
    'etier' = MSE on E[tier] (cheap, stable); 'kl' = KL(clean_detached || view)."""
    pc = probs_clean.detach()
    if metric == "kl":
        return (pc * (pc.clamp_min(1e-8).log() - probs_view.clamp_min(1e-8).log())).sum(-1).mean()
    lv = torch.arange(n_bins, device=probs_view.device, dtype=probs_view.dtype)
    et_clean = (pc * lv).sum(-1)
    et_view = (probs_view * lv).sum(-1)
    return F.mse_loss(et_view, et_clean)
