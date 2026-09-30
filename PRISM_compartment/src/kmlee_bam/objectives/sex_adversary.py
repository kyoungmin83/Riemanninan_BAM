"""v31 — GRL adversary that ERASES sex from z_clean (the post-projection disease-state latent).

The hard multi-rank sex projection removes a fixed subspace, but the encoder REFILLS sex into the
orthogonal complement over epochs (v30: z_clean->sex +0.07 @ep3 -> +0.21 @ep8). A gradient-reversal
adversary adds a SOFT, direction-agnostic pressure: a small classifier tries to predict sex from
z_clean, and the reversed gradient teaches the encoder to hide sex *anywhere* in z — catching the
refill the linear projection can't. Mirrors ReferenceCelltypeAdversary, but:
  - target = sex (2-class: 0=F, 1=M; -1 missing excluded), applied to ALL valid-sex cells (no ref mask),
  - per-BATCH inverse-frequency class weights ⇒ balanced loss, to resist the majority-class collapse
    that broke the tech adversary (bal-acc pinned at chance).
Default-OFF (None ⇒ byte-identical forward). Use a SMALL lambda + warmup; strong early pressure can
erase pathology signal that correlates with sex.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from kmlee_bam.objectives.tech_adversary import GradientReversal


@dataclass
class SexAdversaryOutput:
    logits: torch.Tensor
    loss: torch.Tensor
    accuracy: torch.Tensor
    balanced_accuracy: torch.Tensor
    n_valid: torch.Tensor


class SexAdversary(nn.Module):
    """Predict SEX from GRL(z_perp); the reversed gradient pushes the encoder to erase sex."""

    def __init__(
        self,
        *,
        d_z: int,
        hidden_dim: int = 64,
        dropout: float = 0.1,
        grl_strength: float = 1.0,
    ) -> None:
        super().__init__()
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        self.grl = GradientReversal(grl_strength)
        self.net = nn.Sequential(
            nn.Linear(d_z, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 2),
        )

    def forward(self, z_perp: torch.Tensor, sex_id: Optional[torch.Tensor]) -> SexAdversaryOutput:
        if z_perp.ndim != 2:
            raise ValueError(f"z_perp must have shape [B, d_z], got {tuple(z_perp.shape)}.")
        logits = self.net(self.grl(z_perp))

        z0 = torch.zeros((), dtype=logits.dtype, device=logits.device)
        zero = logits.sum() * 0.0

        # No-label path: keep the return type stable and the graph DDP-visible.
        if sex_id is None:
            return SexAdversaryOutput(
                logits=logits,
                loss=zero,
                accuracy=z0,
                balanced_accuracy=z0,
                n_valid=z0,
            )

        sex_id = sex_id.long().view(-1)
        if sex_id.shape[0] != z_perp.shape[0]:
            raise ValueError("Batch size mismatch between z_perp and sex_id.")
        valid = (sex_id == 0) | (sex_id == 1)       # -1 = missing; anything else is ignored
        n_valid = valid.sum()

        # <2 labelled cells ⇒ no signal; zero loss but keep the graph attached (DDP stable).
        if int(n_valid.item()) < 2:
            return SexAdversaryOutput(logits=logits, loss=zero, accuracy=z0,
                                      balanced_accuracy=z0, n_valid=n_valid.to(dtype=logits.dtype))

        vl = logits[valid]
        vt = sex_id[valid]
        with torch.no_grad():
            cnt = torch.bincount(vt, minlength=2).to(logits.dtype)
            n_present = int((cnt > 0).sum().item())
        # A one-sex batch cannot teach sex invariance. Skip instead of pushing
        # the encoder against a degenerate single-class classifier.
        if n_present < 2:
            return SexAdversaryOutput(logits=logits, loss=zero, accuracy=z0,
                                      balanced_accuracy=z0, n_valid=n_valid.to(dtype=logits.dtype))

        # per-batch inverse-frequency weights → balanced CE (resists majority-class collapse).
        with torch.no_grad():
            w = torch.where(cnt > 0, cnt.sum() / (2.0 * cnt.clamp_min(1.0)), torch.ones_like(cnt))
        loss = F.cross_entropy(vl, vt, weight=w.to(vl.dtype))

        with torch.no_grad():
            pred = vl.argmax(dim=-1)
            accuracy = (pred == vt).float().mean()
            recs = [((pred[vt == c] == c).float().mean()) for c in (0, 1) if (vt == c).any()]
            balanced_accuracy = torch.stack(recs).mean() if recs else accuracy * 0.0

        return SexAdversaryOutput(
            logits=logits,
            loss=loss,
            accuracy=accuracy,
            balanced_accuracy=balanced_accuracy,
            n_valid=n_valid.to(dtype=logits.dtype),
        )
