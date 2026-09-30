"""PRISM-BAM Stage-1 hierarchical personal-response adapter.

The adapter consumes pathology/nuisance-adjusted frozen context summaries.  It
never receives a donor identifier and never sees the target cell type/region in
its source set.  A two-dimensional posterior code is shared by the personal
state and the constrained pathology-response branches.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import torch
from torch import nn


@dataclass(frozen=True)
class PrismAdapterConfig:
    n_contexts: int
    n_celltypes: int
    n_regions: int
    n_modules: int
    latent_dim: int
    n_pathology: int = 5
    personal_rank: int = 2
    hidden_dim: int = 64
    module_token_dim: int = 32
    latent_token_dim: int = 16
    context_embedding_dim: int = 16
    n_heads: int = 4
    n_layers: int = 2
    dropout: float = 0.10
    context_dropout: float = 0.10
    logvar_min: float = -6.0
    logvar_max: float = 2.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _unit_direction(direction: torch.Tensor, eps: float = 1e-8) -> torch.Tensor:
    norm = torch.linalg.vector_norm(direction, dim=-1, keepdim=True)
    return torch.where(norm > eps, direction / norm.clamp_min(eps), torch.zeros_like(direction))


class PrismAdapter(nn.Module):
    """Small masked donor-context adapter with A/B/C/D comparison arms."""

    VALID_ARMS = ("A", "B", "C", "D")

    def __init__(
        self,
        config: PrismAdapterConfig,
        *,
        arm: str,
        context_celltype: torch.Tensor,
        context_region: torch.Tensor,
        axis_direction_module: torch.Tensor,
        axis_direction_latent: torch.Tensor,
    ) -> None:
        super().__init__()
        if arm not in self.VALID_ARMS:
            raise ValueError(f"arm must be one of {self.VALID_ARMS}, got {arm!r}")
        self.config = config
        self.arm = arm
        self.register_buffer("context_celltype", context_celltype.long().clone())
        self.register_buffer("context_region", context_region.long().clone())
        self.register_buffer(
            "axis_direction_module",
            _unit_direction(axis_direction_module.float().clone()),
        )
        self.register_buffer(
            "axis_direction_latent",
            _unit_direction(axis_direction_latent.float().clone()),
        )

        self.module_projection = nn.Sequential(
            nn.Linear(config.n_modules, config.module_token_dim),
            nn.LayerNorm(config.module_token_dim),
            nn.GELU(),
        )
        self.latent_projection = nn.Sequential(
            nn.Linear(config.latent_dim, config.latent_token_dim),
            nn.LayerNorm(config.latent_token_dim),
            nn.GELU(),
        )
        self.celltype_embedding = nn.Embedding(config.n_celltypes, config.context_embedding_dim)
        self.region_embedding = nn.Embedding(config.n_regions, config.context_embedding_dim)
        token_input = (
            config.module_token_dim
            + config.latent_token_dim
            + 2 * config.context_embedding_dim
        )
        self.token_projection = nn.Linear(token_input, config.hidden_dim)
        layer = nn.TransformerEncoderLayer(
            d_model=config.hidden_dim,
            nhead=config.n_heads,
            dim_feedforward=4 * config.hidden_dim,
            dropout=config.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.context_encoder = nn.TransformerEncoder(layer, num_layers=config.n_layers)
        self.posterior = nn.Sequential(
            nn.LayerNorm(config.hidden_dim),
            nn.Linear(config.hidden_dim, config.hidden_dim),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.hidden_dim, 2 * config.personal_rank),
        )

        shape_module = (config.n_contexts, config.personal_rank, config.n_modules)
        shape_latent = (config.n_contexts, config.personal_rank, config.latent_dim)
        self.personal_module_basis = nn.Parameter(torch.empty(shape_module))
        self.personal_latent_basis = nn.Parameter(torch.empty(shape_latent))

        weight_shape = (config.n_contexts, config.n_pathology, config.personal_rank)
        self.response_along_weight = nn.Parameter(torch.empty(weight_shape))
        self.response_orth_weight = nn.Parameter(torch.empty(weight_shape))
        self.response_module_q = nn.Parameter(
            torch.empty(config.n_contexts, config.n_pathology, config.n_modules)
        )
        self.response_latent_q = nn.Parameter(
            torch.empty(config.n_contexts, config.n_pathology, config.latent_dim)
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.personal_module_basis, mean=0.0, std=0.01)
        nn.init.normal_(self.personal_latent_basis, mean=0.0, std=0.01)
        nn.init.normal_(self.response_along_weight, mean=0.0, std=0.01)
        nn.init.normal_(self.response_orth_weight, mean=0.0, std=0.01)
        nn.init.normal_(self.response_module_q, mean=0.0, std=0.01)
        nn.init.normal_(self.response_latent_q, mean=0.0, std=0.01)

    def _drop_source_contexts(self, mask: torch.Tensor) -> torch.Tensor:
        if not self.training or self.config.context_dropout <= 0.0:
            return mask
        keep = torch.rand(mask.shape, device=mask.device) >= self.config.context_dropout
        dropped = mask & keep
        empty = ~dropped.any(dim=1)
        if empty.any():
            for row in torch.nonzero(empty, as_tuple=False).flatten():
                available = torch.nonzero(mask[row], as_tuple=False).flatten()
                if available.numel():
                    dropped[row, available[0]] = True
        return dropped

    def encode(
        self,
        source_module: torch.Tensor,
        source_latent: torch.Tensor,
        source_mask: torch.Tensor,
        source_reliability: torch.Tensor,
        *,
        sample: bool | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Infer the personal posterior without donor ID or pathology labels."""
        mask = self._drop_source_contexts(source_mask.bool())
        reliability = source_reliability.float() * mask.float()
        module_token = self.module_projection(source_module.float())
        latent_token = self.latent_projection(source_latent.float())
        celltype = self.celltype_embedding(self.context_celltype)[None].expand(
            source_module.shape[0], -1, -1
        )
        region = self.region_embedding(self.context_region)[None].expand(
            source_module.shape[0], -1, -1
        )
        token = torch.cat([module_token, latent_token, celltype, region], dim=-1)
        token = self.token_projection(token)
        token = token * reliability.unsqueeze(-1)
        encoded = self.context_encoder(token, src_key_padding_mask=~mask)
        denominator = reliability.sum(dim=1, keepdim=True).clamp_min(1e-6)
        pooled = (encoded * reliability.unsqueeze(-1)).sum(dim=1) / denominator
        posterior = self.posterior(pooled)
        mean, logvar = posterior.chunk(2, dim=-1)
        logvar = logvar.clamp(self.config.logvar_min, self.config.logvar_max)
        if sample is None:
            sample = self.training
        if sample:
            code = mean + torch.exp(0.5 * logvar) * torch.randn_like(mean)
        else:
            code = mean
        return code, mean, logvar

    @staticmethod
    def _orthogonal_q(raw_q: torch.Tensor, common: torch.Tensor) -> torch.Tensor:
        common = _unit_direction(common)
        projected = raw_q - (raw_q * common).sum(dim=-1, keepdim=True) * common
        return _unit_direction(projected)

    def decode(
        self,
        code: torch.Tensor,
        target_context: torch.Tensor,
        pathology: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        batch = code.shape[0]
        zero_module = code.new_zeros((batch, self.config.n_modules))
        zero_latent = code.new_zeros((batch, self.config.latent_dim))
        if self.arm == "A":
            return {
                "module": zero_module,
                "latent": zero_latent,
                "personal_module": zero_module,
                "personal_latent": zero_latent,
                "response_along_module": zero_module,
                "response_orth_module": zero_module,
                "response_along_latent": zero_latent,
                "response_orth_latent": zero_latent,
            }

        personal_module = torch.einsum(
            "br,brm->bm", code, self.personal_module_basis[target_context]
        )
        personal_latent = torch.einsum(
            "br,brz->bz", code, self.personal_latent_basis[target_context]
        )
        along_module = zero_module
        orth_module = zero_module
        along_latent = zero_latent
        orth_latent = zero_latent

        if self.arm in {"C", "D"}:
            along_weight = self.response_along_weight[target_context]
            orth_weight = self.response_orth_weight[target_context]
            if self.arm == "C":
                along_scalar = pathology * torch.einsum("br,bkr->bk", code, along_weight)
                orth_scalar = pathology * torch.einsum("br,bkr->bk", code, orth_weight)
            else:
                # Two fixed nonlinear features per axis.  personal_rank=2 makes
                # C and D exactly parameter matched in the response branch.
                if self.config.personal_rank != 2:
                    raise RuntimeError("matched D control requires personal_rank=2")
                extra = torch.stack([pathology.square(), pathology.pow(3)], dim=-1)
                along_scalar = torch.einsum("bkr,bkr->bk", extra, along_weight)
                orth_scalar = torch.einsum("bkr,bkr->bk", extra, orth_weight)

            common_module = self.axis_direction_module[target_context]
            common_latent = self.axis_direction_latent[target_context]
            q_module = self._orthogonal_q(
                self.response_module_q[target_context], common_module
            )
            q_latent = self._orthogonal_q(
                self.response_latent_q[target_context], common_latent
            )
            along_module = torch.einsum("bk,bkm->bm", along_scalar, common_module)
            orth_module = torch.einsum("bk,bkm->bm", orth_scalar, q_module)
            along_latent = torch.einsum("bk,bkz->bz", along_scalar, common_latent)
            orth_latent = torch.einsum("bk,bkz->bz", orth_scalar, q_latent)

        return {
            "module": personal_module + along_module + orth_module,
            "latent": personal_latent + along_latent + orth_latent,
            "personal_module": personal_module,
            "personal_latent": personal_latent,
            "response_along_module": along_module,
            "response_orth_module": orth_module,
            "response_along_latent": along_latent,
            "response_orth_latent": orth_latent,
        }

    def forward(self, batch: dict[str, torch.Tensor], *, sample: bool | None = None) -> dict[str, torch.Tensor]:
        code, mean, logvar = self.encode(
            batch["source_module"],
            batch["source_latent"],
            batch["source_mask"],
            batch["source_reliability"],
            sample=sample,
        )
        pathology = batch["pathology"].masked_fill(batch["pathology_missing"].bool(), 0.0)
        decoded = self.decode(code, batch["target_context"], pathology)
        decoded.update({"code": code, "posterior_mean": mean, "posterior_logvar": logvar})
        return decoded

    def branch_parameter_count(self) -> dict[str, int]:
        groups = {
            "encoder": [
                self.module_projection,
                self.latent_projection,
                self.celltype_embedding,
                self.region_embedding,
                self.token_projection,
                self.context_encoder,
                self.posterior,
            ],
            "personal": [self.personal_module_basis, self.personal_latent_basis],
            "response": [
                self.response_along_weight,
                self.response_orth_weight,
                self.response_module_q,
                self.response_latent_q,
            ],
        }
        result: dict[str, int] = {}
        for name, objects in groups.items():
            parameters: list[torch.Tensor] = []
            for obj in objects:
                if isinstance(obj, nn.Module):
                    parameters.extend(list(obj.parameters()))
                else:
                    parameters.append(obj)
            result[name] = int(sum(parameter.numel() for parameter in parameters))
        result["total"] = int(sum(parameter.numel() for parameter in self.parameters()))
        return result
