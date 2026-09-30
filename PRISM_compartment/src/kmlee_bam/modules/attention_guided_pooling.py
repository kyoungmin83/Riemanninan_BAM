"""Multi-head attention-guided sequence pooling.

This is a small readout layer, not another self-attention block.  Each pooling
head assigns one scalar score to every valid token, forms a weighted token
average, and the head summaries are concatenated and projected back to the
model dimension.  Padding and an optional CLS token are excluded before the
softmax.

``variant="legacy"`` preserves the original implementation and state-dict
layout. ``variant="diverse_slots_v2"`` uses learned slot queries together
with head-specific key/value projections, a direct 64-to-64 slot readout, and
a bounded mean residual.  The v2 route also exposes differentiable diversity
and entropy-band penalties; the trainer decides whether to add them.
"""

from __future__ import annotations

import math
from typing import Dict, Optional

import torch
from torch import nn


class AttentionGuidedPooling(nn.Module):
    """Pool ``[B,T,D]`` token states into one ``[B,D]`` representation."""

    def __init__(
        self,
        d_model: int,
        n_heads: int = 8,
        *,
        variant: str = "legacy",
        temperature_init: float = 0.75,
        temperature_min: float = 0.25,
        temperature_max: float = 2.0,
        mean_residual_init: float = 0.25,
        mean_residual_max: float = 0.50,
        diversity_weight: float = 0.0,
        diversity_margin: float = 0.50,
        query_orthogonality_weight: float = 0.0,
        entropy_band_weight: float = 0.0,
        min_effective_tokens: float = 1.0,
        max_effective_tokens: float = 1.0e9,
        init_std: float = 0.02,
    ) -> None:
        super().__init__()
        if isinstance(d_model, bool) or int(d_model) <= 0:
            raise ValueError("d_model must be a positive integer")
        if isinstance(n_heads, bool) or int(n_heads) <= 0:
            raise ValueError("n_heads must be a positive integer")
        if float(init_std) <= 0.0:
            raise ValueError("init_std must be positive")

        variant = str(variant).strip().lower()
        if variant not in {"legacy", "diverse_slots_v2"}:
            raise ValueError(
                "variant must be 'legacy' or 'diverse_slots_v2', got "
                f"{variant!r}"
            )
        if not 0.0 < float(temperature_min) < float(temperature_max):
            raise ValueError("temperature bounds must satisfy 0 < min < max")
        if not float(temperature_min) < float(temperature_init) < float(
            temperature_max
        ):
            raise ValueError("temperature_init must lie strictly inside its bounds")
        if not 0.0 <= float(mean_residual_init) <= float(mean_residual_max):
            raise ValueError(
                "mean_residual_init must be in [0, mean_residual_max]"
            )
        if float(mean_residual_max) <= 0.0:
            raise ValueError("mean_residual_max must be positive")
        for name, value in (
            ("diversity_weight", diversity_weight),
            ("query_orthogonality_weight", query_orthogonality_weight),
            ("entropy_band_weight", entropy_band_weight),
        ):
            if float(value) < 0.0:
                raise ValueError(f"{name} must be non-negative")
        if not -1.0 <= float(diversity_margin) <= 1.0:
            raise ValueError("diversity_margin must be in [-1, 1]")
        if not 1.0 <= float(min_effective_tokens) <= float(max_effective_tokens):
            raise ValueError(
                "effective-token bounds must satisfy 1 <= min <= max"
            )

        self.d_model = int(d_model)
        self.n_heads = int(n_heads)
        self.variant = variant
        self.diversity_weight = float(diversity_weight)
        self.diversity_margin = float(diversity_margin)
        self.query_orthogonality_weight = float(query_orthogonality_weight)
        self.entropy_band_weight = float(entropy_band_weight)
        self.min_effective_tokens = float(min_effective_tokens)
        self.max_effective_tokens = float(max_effective_tokens)

        if self.variant == "legacy":
            # Keep these exact names/shapes so every historical AGP checkpoint
            # remains loadable without a migration.
            self.score = nn.Linear(self.d_model, self.n_heads, bias=False)
            self.combine = nn.Linear(
                self.n_heads * self.d_model,
                self.d_model,
                bias=True,
            )
            self.key_projection = None
            self.value_projection = None
            self.slot_query = None
            self.raw_temperature = None
            self.output_norm = None
            self.raw_mean_residual_gate = None
            self.temperature_min = float(temperature_min)
            self.temperature_max = float(temperature_max)
            self.mean_residual_max = float(mean_residual_max)
        else:
            if self.d_model % self.n_heads != 0:
                raise ValueError(
                    "diverse_slots_v2 requires d_model divisible by n_heads; "
                    f"got d_model={self.d_model}, n_heads={self.n_heads}"
                )
            self.head_dim = self.d_model // self.n_heads
            self.score = None
            self.combine = None
            self.key_projection = nn.Linear(
                self.d_model, self.n_heads * self.head_dim, bias=False
            )
            self.value_projection = nn.Linear(
                self.d_model, self.n_heads * self.head_dim, bias=False
            )
            self.slot_query = nn.Parameter(torch.empty(self.n_heads, self.head_dim))
            self.temperature_min = float(temperature_min)
            self.temperature_max = float(temperature_max)
            temperature_fraction = (
                (float(temperature_init) - self.temperature_min)
                / (self.temperature_max - self.temperature_min)
            )
            temperature_logit = math.log(
                temperature_fraction / (1.0 - temperature_fraction)
            )
            self.raw_temperature = nn.Parameter(
                torch.full((self.n_heads,), temperature_logit)
            )
            self.output_norm = nn.LayerNorm(self.d_model)
            self.mean_residual_max = float(mean_residual_max)
            residual_fraction = min(
                max(float(mean_residual_init) / self.mean_residual_max, 1.0e-6),
                1.0 - 1.0e-6,
            )
            residual_logit = math.log(residual_fraction / (1.0 - residual_fraction))
            self.raw_mean_residual_gate = nn.Parameter(
                torch.tensor(residual_logit)
            )
        self.reset_parameters(float(init_std))

    def reset_parameters(self, init_std: float) -> None:
        if self.variant == "legacy":
            assert self.score is not None and self.combine is not None
            nn.init.normal_(self.score.weight, mean=0.0, std=init_std)
            nn.init.normal_(self.combine.weight, mean=0.0, std=init_std)
            nn.init.zeros_(self.combine.bias)
            return

        assert self.key_projection is not None
        assert self.value_projection is not None
        assert self.slot_query is not None
        assert self.output_norm is not None
        nn.init.xavier_uniform_(self.key_projection.weight)
        # The concatenated slots already have dimension D.  Identity value
        # initialisation keeps those coordinates aligned with the bounded
        # global-mean residual at scratch start; learning can subsequently
        # rotate them without a 512-to-64 bottleneck.
        nn.init.eye_(self.value_projection.weight)
        # Orthogonal rows prevent identical questions at initialisation.  The
        # sqrt(head_dim) norm yields O(1) scaled dot products instead of an
        # almost-uniform 414-token softmax at the first update.
        nn.init.orthogonal_(self.slot_query)
        with torch.no_grad():
            self.slot_query.mul_(math.sqrt(float(self.head_dim)))
        nn.init.ones_(self.output_norm.weight)
        nn.init.zeros_(self.output_norm.bias)

    def temperature(self) -> torch.Tensor:
        """Return bounded, independently learned per-head temperatures."""

        if self.raw_temperature is None:
            return torch.ones(self.n_heads, device=self.score.weight.device)
        return self.temperature_min + (
            self.temperature_max - self.temperature_min
        ) * torch.sigmoid(self.raw_temperature)

    def mean_residual_gate(self) -> torch.Tensor:
        """Return the bounded strength of the information-preserving mean route."""

        if self.raw_mean_residual_gate is None:
            if self.score is None:
                raise RuntimeError("AGP module has no parameters")
            return torch.zeros((), device=self.score.weight.device)
        return self.mean_residual_max * torch.sigmoid(
            self.raw_mean_residual_gate
        )

    def forward(
        self,
        token_states: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor] = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if token_states.ndim != 3:
            raise ValueError(
                "token_states must have shape [B,T,D], got "
                f"{tuple(token_states.shape)}"
            )
        batch, tokens, dim = token_states.shape
        if dim != self.d_model:
            raise ValueError(f"expected D={self.d_model}, got {dim}")
        if tokens == 0:
            raise ValueError("attention-guided pooling received zero tokens")

        if key_padding_mask is not None:
            if key_padding_mask.shape != (batch, tokens):
                raise ValueError(
                    "key_padding_mask must have shape "
                    f"{(batch, tokens)}, got {tuple(key_padding_mask.shape)}"
                )
            if key_padding_mask.dtype != torch.bool:
                raise TypeError("key_padding_mask must be torch.bool")
            if bool(key_padding_mask.all(dim=1).any()):
                raise ValueError("every sample needs at least one valid pooling token")

        if self.variant == "legacy":
            return self._forward_legacy(
                token_states,
                key_padding_mask=key_padding_mask,
            )
        return self._forward_diverse_slots_v2(
            token_states,
            key_padding_mask=key_padding_mask,
        )

    def _forward_legacy(
        self,
        token_states: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        assert self.score is not None and self.combine is not None
        batch = int(token_states.shape[0])
        scores = self.score(token_states).float()
        if key_padding_mask is not None:
            scores = scores.masked_fill(key_padding_mask.unsqueeze(-1), -torch.inf)
        weights = torch.softmax(scores, dim=1)
        head_summary = torch.einsum(
            "bth,btd->bhd",
            weights.to(token_states.dtype),
            token_states,
        )
        pooled = self.combine(head_summary.reshape(batch, -1))
        return pooled, weights

    def _forward_diverse_slots_v2(
        self,
        token_states: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        assert self.key_projection is not None
        assert self.value_projection is not None
        assert self.slot_query is not None
        assert self.output_norm is not None
        batch, tokens, _ = token_states.shape

        keys = self.key_projection(token_states).reshape(
            batch, tokens, self.n_heads, self.head_dim
        )
        values = self.value_projection(token_states).reshape(
            batch, tokens, self.n_heads, self.head_dim
        )
        scores = torch.einsum(
            "bthd,hd->bth",
            keys.float(),
            self.slot_query.float(),
        )
        scores = scores / math.sqrt(float(self.head_dim))
        scores = scores / self.temperature().float().view(1, 1, self.n_heads)
        if key_padding_mask is not None:
            scores = scores.masked_fill(key_padding_mask.unsqueeze(-1), -torch.inf)
        weights = torch.softmax(scores, dim=1)
        slots = torch.einsum(
            "bth,bthd->bhd",
            weights.to(values.dtype),
            values,
        ).reshape(batch, self.d_model)

        if key_padding_mask is None:
            global_mean = token_states.mean(dim=1)
        else:
            valid = (~key_padding_mask).to(token_states.dtype).unsqueeze(-1)
            global_mean = (token_states * valid).sum(dim=1) / valid.sum(
                dim=1
            ).clamp_min(1.0)
        pooled = self.output_norm(
            slots + self.mean_residual_gate().to(slots.dtype) * global_mean
        )
        return pooled, weights

    def auxiliary_terms(self, weights: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Compute v2 anti-collapse penalties and differentiable diagnostics.

        The attention-diversity loss uses *centred* token weights.  This avoids
        declaring two broad heads identical merely because every probability
        distribution shares a uniform background.  Only excessive positive
        similarity is penalised; complementary (negative-correlation) heads
        remain allowed.
        """

        if weights.ndim != 3:
            raise ValueError("weights must have shape [B,T,H]")
        w = weights.float().clamp_min(0.0)
        valid = w.sum(dim=2) > 0.0
        valid_float = valid.to(w.dtype)
        n_valid = valid_float.sum(dim=1, keepdim=True).clamp_min(1.0)
        uniform = valid_float / n_valid
        by_head = w.transpose(1, 2)
        centred = by_head - uniform.unsqueeze(1)
        centred_unit = centred / centred.norm(
            dim=2, keepdim=True
        ).clamp_min(1.0e-8)
        centred_similarity = torch.bmm(
            centred_unit, centred_unit.transpose(1, 2)
        )
        n_heads = int(centred_similarity.shape[1])
        if n_heads > 1:
            off_diagonal = ~torch.eye(
                n_heads,
                dtype=torch.bool,
                device=centred_similarity.device,
            )
            centred_overlap = centred_similarity[:, off_diagonal].mean()
            diversity = torch.relu(
                centred_similarity[:, off_diagonal] - self.diversity_margin
            ).square().mean()
        else:
            centred_overlap = w.new_zeros(())
            diversity = w.new_zeros(())

        entropy = -(w.clamp_min(1.0e-12) * w.clamp_min(1.0e-12).log()).sum(
            dim=1
        )
        log_min = math.log(self.min_effective_tokens)
        log_max = math.log(self.max_effective_tokens)
        entropy_band = (
            torch.relu(log_min - entropy).square()
            + torch.relu(entropy - log_max).square()
        ).mean()

        query_orthogonality = w.new_zeros(())
        query_effective_rank = w.new_tensor(1.0)
        score_effective_rank = w.new_tensor(1.0)
        if self.slot_query is not None:
            query = self.slot_query.float()
            query_unit = query / query.norm(dim=1, keepdim=True).clamp_min(1.0e-8)
            gram = query_unit @ query_unit.transpose(0, 1)
            if self.n_heads > 1:
                off_diagonal = ~torch.eye(
                    self.n_heads,
                    dtype=torch.bool,
                    device=gram.device,
                )
                query_orthogonality = gram[off_diagonal].square().mean()
            # These are diagnostics only.  Explicitly leave AMP because CUDA
            # SVD has no BF16 kernel and the surrounding PHU bootstrap runs
            # under the training autocast context.
            with torch.no_grad(), torch.autocast(
                device_type=w.device.type, enabled=False
            ):
                query_fp32 = self.slot_query.detach().to(torch.float32)
                singular = torch.linalg.svdvals(query_fp32)
                singular_sq = singular.square()
                query_effective_rank = singular_sq.sum().square() / singular_sq.square().sum().clamp_min(
                    1.0e-12
                )
                assert self.key_projection is not None
                key_weight = self.key_projection.weight.detach().to(
                    torch.float32
                ).reshape(self.n_heads, self.head_dim, self.d_model)
                score_directions = torch.einsum(
                    "hd,hdk->hk", query_fp32, key_weight
                )
                score_singular = torch.linalg.svdvals(score_directions)
                score_sq = score_singular.square()
                score_effective_rank = score_sq.sum().square() / score_sq.square().sum().clamp_min(
                    1.0e-12
                )
        elif self.score is not None:
            with torch.no_grad(), torch.autocast(
                device_type=w.device.type, enabled=False
            ):
                score_singular = torch.linalg.svdvals(
                    self.score.weight.detach().to(torch.float32)
                )
                score_sq = score_singular.square()
                score_effective_rank = score_sq.sum().square() / score_sq.square().sum().clamp_min(
                    1.0e-12
                )

        weighted_total = (
            self.diversity_weight * diversity
            + self.query_orthogonality_weight * query_orthogonality
            + self.entropy_band_weight * entropy_band
        )
        return {
            "weighted_total": weighted_total,
            "diversity": diversity,
            "centred_overlap": centred_overlap,
            "entropy_band": entropy_band,
            "query_orthogonality": query_orthogonality,
            "query_effective_rank": query_effective_rank,
            "score_effective_rank": score_effective_rank,
            "temperature_mean": self.temperature().float().mean(),
            "mean_residual_gate": self.mean_residual_gate().float(),
        }
