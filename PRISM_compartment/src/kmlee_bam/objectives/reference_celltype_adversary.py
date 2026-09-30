from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from kmlee_bam.objectives.tech_adversary import GradientReversal
except ImportError:
    from kmlee_bam.objectives.tech_adversary import GradientReversal


@dataclass
class ReferenceCelltypeAdversaryOutput:
    logits: torch.Tensor
    loss: torch.Tensor
    accuracy: torch.Tensor
    balanced_accuracy: torch.Tensor
    n_reference: torch.Tensor
    true_counts: torch.Tensor
    pred_counts: torch.Tensor
    correct_counts: torch.Tensor
    recall_per_class: torch.Tensor


class ReferenceCelltypeAdversary(nn.Module):
    """
    Reference-only adversary that predicts CELLTYPE from GRL(z_perp).

    This targets reference-conditional invariance of the disease-state latent:

        z_perp should NOT encode cell-type identity *on reference/control cells*.

    A probe showed reference cells' z_perp still decodes celltype at
    balanced-acc 0.785 (lift +0.74) — cell-type identity leaks into the
    disease-state latent. The gradient reversal pushes the encoder to erase that
    identity, but the cross-entropy / GRL pressure is applied to REFERENCE cells
    only (mask by ``is_reference``). Disease cells are excluded so they keep their
    legitimate cell-type-specific disease direction.

    Unlike :class:`ConditionalTechAdversary`, celltype is the *target* here, not a
    conditioning variable, so there is no celltype embedding and the classifier
    consumes GRL(z_perp) alone.

    Use a small loss weight and warm it up after the reference origin has started
    to stabilise. Strong early adversarial pressure can erase pathology signal.
    """

    def __init__(
        self,
        *,
        d_z: int,
        n_celltypes: int,
        hidden_dim: int = 64,
        dropout: float = 0.1,
        grl_strength: float = 1.0,
        class_weights: Optional[torch.Tensor] = None,
    ) -> None:
        super().__init__()
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if n_celltypes <= 1:
            raise ValueError("n_celltypes must be greater than one.")

        self.grl = GradientReversal(grl_strength)
        self.net = nn.Sequential(
            nn.Linear(d_z, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_celltypes),
        )
        self.n_celltypes = int(n_celltypes)
        if class_weights is None:
            self.register_buffer("class_weights", None, persistent=True)
        else:
            weights = torch.as_tensor(class_weights, dtype=torch.float32)
            if weights.ndim != 1 or weights.numel() != n_celltypes:
                raise ValueError(
                    "class_weights must have shape [n_celltypes], got "
                    f"{tuple(weights.shape)} for n_celltypes={n_celltypes}."
                )
            if not torch.isfinite(weights).all():
                raise ValueError("class_weights must be finite.")
            if (weights <= 0).any():
                raise ValueError("class_weights must be positive.")
            self.register_buffer("class_weights", weights, persistent=True)

    def forward(
        self,
        z_perp: torch.Tensor,
        celltype_id: torch.Tensor,
        is_reference: Optional[torch.Tensor] = None,
    ) -> ReferenceCelltypeAdversaryOutput | torch.Tensor:
        if z_perp.ndim != 2:
            raise ValueError(f"z_perp must have shape [B, d_z], got {tuple(z_perp.shape)}.")
        if celltype_id.ndim != 1:
            raise ValueError("celltype_id must have shape [B].")
        if celltype_id.shape[0] != z_perp.shape[0]:
            raise ValueError("Batch size mismatch between z_perp and celltype_id.")

        z = self.grl(z_perp)
        logits = self.net(z)

        # Inference / no-supervision path: return raw logits over all cells.
        if is_reference is None:
            return logits

        is_reference = is_reference.bool()
        if is_reference.ndim != 1:
            raise ValueError(
                f"is_reference must have shape [B], got {tuple(is_reference.shape)}."
            )
        if is_reference.shape[0] != z_perp.shape[0]:
            raise ValueError("Batch size mismatch between z_perp and is_reference.")

        celltype_id = celltype_id.long()
        n_reference = is_reference.sum()

        # No reference cells in this batch ⇒ no adversarial signal. Return zero
        # loss / balanced_acc but keep the graph attached to logits (the GRL has
        # already touched z_perp) so DDP unused-parameter traversal stays stable.
        if int(n_reference.item()) == 0:
            zero = logits.sum() * 0.0
            n_zeros = torch.zeros(self.n_celltypes, dtype=torch.long, device=logits.device)
            zero_recall = torch.zeros(self.n_celltypes, dtype=logits.dtype, device=logits.device)
            return ReferenceCelltypeAdversaryOutput(
                logits=logits,
                loss=zero,
                accuracy=torch.zeros((), dtype=logits.dtype, device=logits.device),
                balanced_accuracy=torch.zeros((), dtype=logits.dtype, device=logits.device),
                n_reference=n_reference.to(dtype=logits.dtype),
                true_counts=n_zeros,
                pred_counts=n_zeros,
                correct_counts=n_zeros,
                recall_per_class=zero_recall,
            )

        ref_logits = logits[is_reference]
        ref_targets = celltype_id[is_reference]

        loss = F.cross_entropy(ref_logits, ref_targets, weight=self.class_weights)
        pred = ref_logits.argmax(dim=-1)
        accuracy = (pred == ref_targets).float().mean()

        with torch.no_grad():
            true_counts = torch.bincount(ref_targets, minlength=self.n_celltypes).to(logits.device)
            pred_counts = torch.bincount(pred, minlength=self.n_celltypes).to(logits.device)
            correct_counts = torch.bincount(
                ref_targets[pred == ref_targets],
                minlength=self.n_celltypes,
            ).to(logits.device)
            true_counts_f = true_counts.to(dtype=logits.dtype)
            correct_counts_f = correct_counts.to(dtype=logits.dtype)
            present = true_counts_f > 0
            recall_per_class = torch.zeros(
                self.n_celltypes, dtype=logits.dtype, device=logits.device
            )
            recall_per_class[present] = correct_counts_f[present] / true_counts_f[present]
            if present.any():
                balanced_accuracy = recall_per_class[present].mean()
            else:
                balanced_accuracy = torch.zeros((), dtype=logits.dtype, device=logits.device)

        return ReferenceCelltypeAdversaryOutput(
            logits=logits,
            loss=loss,
            accuracy=accuracy,
            balanced_accuracy=balanced_accuracy,
            n_reference=n_reference.to(dtype=logits.dtype),
            true_counts=true_counts,
            pred_counts=pred_counts,
            correct_counts=correct_counts,
            recall_per_class=recall_per_class,
        )
