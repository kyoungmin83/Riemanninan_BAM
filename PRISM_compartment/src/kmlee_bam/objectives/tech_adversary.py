from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


class _GradientReversalFn(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x: torch.Tensor, strength: float) -> torch.Tensor:
        ctx.strength = float(strength)
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> tuple[torch.Tensor, None]:
        return -ctx.strength * grad_output, None


class GradientReversal(nn.Module):
    """
    Identity in the forward pass, negative scaled gradient in the backward pass.
    """

    def __init__(self, strength: float = 1.0) -> None:
        super().__init__()
        self.strength = float(strength)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return _GradientReversalFn.apply(x, self.strength)


@dataclass
class TechAdversaryOutput:
    logits: torch.Tensor
    loss: torch.Tensor
    accuracy: torch.Tensor
    balanced_accuracy: torch.Tensor
    true_counts: torch.Tensor
    pred_counts: torch.Tensor
    correct_counts: torch.Tensor
    recall_per_class: torch.Tensor


class ConditionalTechAdversary(nn.Module):
    """
    Predicts technology from [GRL(z_dev), celltype_embedding].

    This targets conditional invariance:

        z_dev should not predict tech_id within a cell type.

    The gradient reversal is applied only to z_dev. The cell type embedding is
    a conditioning variable for the classifier, not a representation that should
    be adversarially erased.

    Use a small loss weight and warm it up after the reference origin has started
    to stabilise. Strong early adversarial pressure can erase pathology signal.
    """

    def __init__(
        self,
        *,
        d_z: int,
        n_celltypes: int,
        n_tech: int,
        celltype_embed_dim: int = 16,
        hidden_dim: int = 64,
        dropout: float = 0.1,
        grl_strength: float = 1.0,
        class_weights: Optional[torch.Tensor] = None,
    ) -> None:
        super().__init__()
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if n_celltypes <= 0:
            raise ValueError("n_celltypes must be positive.")
        if n_tech <= 1:
            raise ValueError("n_tech must be greater than one.")

        self.celltype_embedding = nn.Embedding(n_celltypes, celltype_embed_dim)
        self.grl = GradientReversal(grl_strength)
        self.net = nn.Sequential(
            nn.Linear(d_z + celltype_embed_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_tech),
        )
        self.n_tech = int(n_tech)
        if class_weights is None:
            self.register_buffer("class_weights", None, persistent=True)
        else:
            weights = torch.as_tensor(class_weights, dtype=torch.float32)
            if weights.ndim != 1 or weights.numel() != n_tech:
                raise ValueError(
                    "class_weights must have shape [n_tech], got "
                    f"{tuple(weights.shape)} for n_tech={n_tech}."
                )
            if not torch.isfinite(weights).all():
                raise ValueError("class_weights must be finite.")
            if (weights <= 0).any():
                raise ValueError("class_weights must be positive.")
            self.register_buffer("class_weights", weights, persistent=True)

    def forward(
        self,
        z_dev: torch.Tensor,
        celltype_id: torch.Tensor,
        tech_id: Optional[torch.Tensor] = None,
    ) -> TechAdversaryOutput | torch.Tensor:
        if z_dev.ndim != 2:
            raise ValueError(f"z_dev must have shape [B, d_z], got {tuple(z_dev.shape)}.")
        if celltype_id.ndim != 1:
            raise ValueError("celltype_id must have shape [B].")
        if celltype_id.shape[0] != z_dev.shape[0]:
            raise ValueError("Batch size mismatch between z_dev and celltype_id.")

        z = self.grl(z_dev)
        c = self.celltype_embedding(celltype_id.long())
        h = torch.cat([z, c], dim=-1)
        logits = self.net(h)

        if tech_id is None:
            return logits

        tech_id = tech_id.long()
        loss = F.cross_entropy(logits, tech_id, weight=self.class_weights)
        pred = logits.argmax(dim=-1)
        accuracy = (pred == tech_id).float().mean()

        with torch.no_grad():
            true_counts = torch.bincount(tech_id, minlength=self.n_tech).to(logits.device)
            pred_counts = torch.bincount(pred, minlength=self.n_tech).to(logits.device)
            correct_counts = torch.bincount(
                tech_id[pred == tech_id],
                minlength=self.n_tech,
            ).to(logits.device)
            true_counts_f = true_counts.to(dtype=logits.dtype)
            correct_counts_f = correct_counts.to(dtype=logits.dtype)
            present = true_counts_f > 0
            recall_per_class = torch.zeros(self.n_tech, dtype=logits.dtype, device=logits.device)
            recall_per_class[present] = correct_counts_f[present] / true_counts_f[present]
            if present.any():
                balanced_accuracy = recall_per_class[present].mean()
            else:
                balanced_accuracy = torch.zeros((), dtype=logits.dtype, device=logits.device)

        return TechAdversaryOutput(
            logits=logits,
            loss=loss,
            accuracy=accuracy,
            balanced_accuracy=balanced_accuracy,
            true_counts=true_counts,
            pred_counts=pred_counts,
            correct_counts=correct_counts,
            recall_per_class=recall_per_class,
        )
