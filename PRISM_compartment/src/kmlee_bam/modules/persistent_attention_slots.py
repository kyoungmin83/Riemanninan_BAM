"""Persistent attention-guided pooling over successive encoder stages.

Unlike one-shot AGP, the same learned specialist slots read the tokenizer
output and every Transformer-layer output in sequence.  The module tokens do
not read the slots, so this is a readout path and does not alter the module
self-attention contract.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence

import torch
from torch import nn


def _bounded_logit(value: float, maximum: float) -> float:
    fraction = min(max(float(value) / float(maximum), 1.0e-6), 1.0 - 1.0e-6)
    return math.log(fraction / (1.0 - fraction))


class PersistentAttentionSlots(nn.Module):
    """Read ``[B,T,D]`` states from several stages into persistent slots."""

    def __init__(
        self,
        d_model: int,
        n_slots: int = 8,
        n_attention_heads: int = 4,
        n_stages: int = 4,
        *,
        ffn_dim: Optional[int] = None,
        dropout: float = 0.1,
        init_gate: float = 0.10,
        gate_max: float = 1.0,
        mean_residual_init: float = 0.25,
        mean_residual_max: float = 0.50,
        diversity_weight: float = 0.0,
        diversity_margin: float = 0.50,
        query_orthogonality_weight: float = 0.0,
        entropy_band_weight: float = 0.0,
        min_effective_tokens: float = 1.0,
        max_effective_tokens: float = 1.0e9,
        layer_norm_eps: float = 1.0e-5,
    ) -> None:
        super().__init__()
        if int(d_model) <= 0 or int(n_slots) <= 0 or int(n_stages) <= 0:
            raise ValueError("d_model, n_slots, and n_stages must be positive")
        if int(n_attention_heads) <= 0 or int(d_model) % int(n_attention_heads) != 0:
            raise ValueError("n_attention_heads must be positive and divide d_model")
        if int(d_model) % int(n_slots) != 0:
            raise ValueError("d_model must be divisible by n_slots")
        if not 0.0 <= float(init_gate) <= float(gate_max):
            raise ValueError("init_gate must lie in [0, gate_max]")
        if float(gate_max) <= 0.0:
            raise ValueError("gate_max must be positive")
        if not 0.0 <= float(mean_residual_init) <= float(mean_residual_max):
            raise ValueError("mean_residual_init must lie in [0, mean_residual_max]")
        if float(mean_residual_max) <= 0.0:
            raise ValueError("mean_residual_max must be positive")
        if not -1.0 <= float(diversity_margin) <= 1.0:
            raise ValueError("diversity_margin must lie in [-1, 1]")
        if not 1.0 <= float(min_effective_tokens) <= float(max_effective_tokens):
            raise ValueError("effective-token bounds must satisfy 1 <= min <= max")

        self.d_model = int(d_model)
        self.n_slots = int(n_slots)
        self.n_attention_heads = int(n_attention_heads)
        self.n_stages = int(n_stages)
        self.slot_output_dim = self.d_model // self.n_slots
        self.gate_max = float(gate_max)
        self.mean_residual_max = float(mean_residual_max)
        self.diversity_weight = float(diversity_weight)
        self.diversity_margin = float(diversity_margin)
        self.query_orthogonality_weight = float(query_orthogonality_weight)
        self.entropy_band_weight = float(entropy_band_weight)
        self.min_effective_tokens = float(min_effective_tokens)
        self.max_effective_tokens = float(max_effective_tokens)

        hidden = int(ffn_dim) if ffn_dim is not None else 2 * self.d_model
        self.slot_identity = nn.Parameter(torch.empty(self.n_slots, self.d_model))
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=self.d_model,
            num_heads=self.n_attention_heads,
            dropout=float(dropout),
            batch_first=True,
        )
        self.shared_ffn = nn.Sequential(
            nn.Linear(self.d_model, hidden),
            nn.GELU(),
            nn.Dropout(float(dropout)),
            nn.Linear(hidden, self.d_model),
            nn.Dropout(float(dropout)),
        )
        self.slot_norms = nn.ModuleList(
            nn.LayerNorm(self.d_model, eps=layer_norm_eps) for _ in range(self.n_stages)
        )
        self.memory_norms = nn.ModuleList(
            nn.LayerNorm(self.d_model, eps=layer_norm_eps) for _ in range(self.n_stages)
        )
        self.ffn_norms = nn.ModuleList(
            nn.LayerNorm(self.d_model, eps=layer_norm_eps) for _ in range(self.n_stages)
        )
        gate_logit = _bounded_logit(float(init_gate), self.gate_max)
        self.raw_cross_gates = nn.Parameter(torch.full((self.n_stages,), gate_logit))
        self.raw_ffn_gates = nn.Parameter(torch.full((self.n_stages,), gate_logit))
        self.raw_initial_mean_gate = nn.Parameter(
            torch.tensor(_bounded_logit(mean_residual_init, self.mean_residual_max))
        )
        self.raw_final_mean_gate = nn.Parameter(
            torch.tensor(_bounded_logit(mean_residual_init, self.mean_residual_max))
        )
        self.slot_output = nn.ModuleList(
            nn.Linear(self.d_model, self.slot_output_dim) for _ in range(self.n_slots)
        )
        self.output_norm = nn.LayerNorm(self.d_model, eps=layer_norm_eps)
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.orthogonal_(self.slot_identity)
        for projection in self.slot_output:
            nn.init.xavier_uniform_(projection.weight)
            nn.init.zeros_(projection.bias)

    def cross_gates(self) -> torch.Tensor:
        return self.gate_max * torch.sigmoid(self.raw_cross_gates)

    def ffn_gates(self) -> torch.Tensor:
        return self.gate_max * torch.sigmoid(self.raw_ffn_gates)

    def initial_mean_gate(self) -> torch.Tensor:
        return self.mean_residual_max * torch.sigmoid(self.raw_initial_mean_gate)

    def mean_residual_gate(self) -> torch.Tensor:
        return self.mean_residual_max * torch.sigmoid(self.raw_final_mean_gate)

    @staticmethod
    def _masked_mean(
        values: torch.Tensor, mask: Optional[torch.Tensor]
    ) -> torch.Tensor:
        if mask is None:
            return values.mean(dim=1)
        valid = (~mask).to(values.dtype).unsqueeze(-1)
        return (values * valid).sum(dim=1) / valid.sum(dim=1).clamp_min(1.0)

    def forward(
        self,
        stage_states: Sequence[torch.Tensor],
        *,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if len(stage_states) != self.n_stages:
            raise ValueError(
                f"expected {self.n_stages} stage states, got {len(stage_states)}"
            )
        first = stage_states[0]
        if first.ndim != 3 or int(first.shape[-1]) != self.d_model:
            raise ValueError("each stage state must have shape [B,T,d_model]")
        batch, tokens, _ = first.shape
        if tokens <= 0:
            raise ValueError("persistent AGP received zero tokens")
        if key_padding_mask is not None:
            if key_padding_mask.shape != (batch, tokens):
                raise ValueError("key_padding_mask shape does not match stage states")
            if key_padding_mask.dtype != torch.bool:
                raise TypeError("key_padding_mask must be torch.bool")
            if bool(key_padding_mask.all(dim=1).any()):
                raise ValueError("every sample needs at least one valid token")

        first_mean = self._masked_mean(first, key_padding_mask)
        slots = self.slot_identity.to(first.dtype).unsqueeze(0).expand(batch, -1, -1)
        slots = slots + self.initial_mean_gate().to(first.dtype) * first_mean.unsqueeze(1)
        all_weights = []
        for stage, memory in enumerate(stage_states):
            if memory.shape != first.shape:
                raise ValueError("all stage states must have identical shapes")
            query = self.slot_norms[stage](slots)
            key_value = self.memory_norms[stage](memory)
            update, weights = self.cross_attention(
                query,
                key_value,
                key_value,
                key_padding_mask=key_padding_mask,
                need_weights=True,
                average_attn_weights=True,
            )
            # MultiheadAttention returns post-dropout weights during training,
            # whose rows need not sum exactly to one.  Diagnostics and the AGP
            # auxiliary contract consume probability distributions, so expose
            # an explicitly renormalised copy without changing the update.
            weights = weights.clamp_min(0.0)
            weights = weights / weights.sum(dim=-1, keepdim=True).clamp_min(1.0e-12)
            slots = slots + self.cross_gates()[stage].to(slots.dtype) * update
            slots = slots + self.ffn_gates()[stage].to(slots.dtype) * self.shared_ffn(
                self.ffn_norms[stage](slots)
            )
            all_weights.append(weights)

        slot_blocks = [
            projection(slots[:, index, :])
            for index, projection in enumerate(self.slot_output)
        ]
        slot_readout = torch.cat(slot_blocks, dim=-1)
        final_mean = self._masked_mean(stage_states[-1], key_padding_mask)
        pooled = self.output_norm(
            slot_readout + self.mean_residual_gate().to(slot_readout.dtype) * final_mean
        )
        # Existing diagnostics consume [B,T,H].  H denotes specialist slots.
        final_weights = all_weights[-1].transpose(1, 2).contiguous()
        stage_weights = torch.stack(all_weights, dim=1)
        return pooled, final_weights, stage_weights

    @staticmethod
    def _effective_rank(matrix: torch.Tensor) -> torch.Tensor:
        singular = torch.linalg.svdvals(matrix.to(torch.float32))
        square = singular.square()
        return square.sum().square() / square.square().sum().clamp_min(1.0e-12)

    def auxiliary_terms(
        self,
        weights: torch.Tensor,
        *,
        stage_weights: Optional[torch.Tensor] = None,
    ) -> Dict[str, torch.Tensor]:
        if weights.ndim != 3:
            raise ValueError("weights must have shape [B,T,S]")
        w = weights.float().clamp_min(0.0)
        valid = w.sum(dim=2) > 0.0
        valid_float = valid.to(w.dtype)
        uniform = valid_float / valid_float.sum(dim=1, keepdim=True).clamp_min(1.0)
        by_slot = w.transpose(1, 2)
        centred = by_slot - uniform.unsqueeze(1)
        centred_unit = centred / centred.norm(dim=2, keepdim=True).clamp_min(1.0e-8)
        similarity = torch.bmm(centred_unit, centred_unit.transpose(1, 2))
        if self.n_slots > 1:
            off = ~torch.eye(self.n_slots, dtype=torch.bool, device=w.device)
            centred_overlap = similarity[:, off].mean()
            diversity = torch.relu(similarity[:, off] - self.diversity_margin).square().mean()
        else:
            centred_overlap = w.new_zeros(())
            diversity = w.new_zeros(())

        entropy = -(w.clamp_min(1.0e-12) * w.clamp_min(1.0e-12).log()).sum(dim=1)
        entropy_band = (
            torch.relu(math.log(self.min_effective_tokens) - entropy).square()
            + torch.relu(entropy - math.log(self.max_effective_tokens)).square()
        ).mean()

        query = self.slot_identity.float()
        query_unit = query / query.norm(dim=1, keepdim=True).clamp_min(1.0e-8)
        gram = query_unit @ query_unit.transpose(0, 1)
        if self.n_slots > 1:
            off = ~torch.eye(self.n_slots, dtype=torch.bool, device=w.device)
            query_orthogonality = gram[off].square().mean()
        else:
            query_orthogonality = w.new_zeros(())
        with torch.no_grad(), torch.autocast(device_type=w.device.type, enabled=False):
            query_rank = self._effective_rank(self.slot_identity.detach())
            q_weight = self.cross_attention.in_proj_weight[: self.d_model].detach()
            score_directions = self.slot_identity.detach().to(torch.float32) @ q_weight.to(torch.float32).T
            score_rank = self._effective_rank(score_directions)

        weighted_total = (
            self.diversity_weight * diversity
            + self.query_orthogonality_weight * query_orthogonality
            + self.entropy_band_weight * entropy_band
        )
        result = {
            "weighted_total": weighted_total,
            "diversity": diversity,
            "centred_overlap": centred_overlap,
            "entropy_band": entropy_band,
            "query_orthogonality": query_orthogonality,
            "query_effective_rank": query_rank,
            "score_effective_rank": score_rank,
            "temperature_mean": w.new_tensor(1.0),
            "mean_residual_gate": self.mean_residual_gate().float(),
            "slot_state_effective_rank": score_rank,
        }
        if stage_weights is not None:
            stage_w = stage_weights.float().clamp_min(1.0e-12)
            stage_entropy = -(stage_w * stage_w.log()).sum(dim=-1).exp().mean(dim=(0, 2))
            for stage, value in enumerate(stage_entropy):
                result[f"stage_{stage}_effective_tokens"] = value
        return result
