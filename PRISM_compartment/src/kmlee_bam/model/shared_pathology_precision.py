"""PRISM Stage-2 personal-state and hierarchical shared-pathology adapter.

The personal encoder is never given donor identity or pathology tokens.  The
shared branch only sees the five named pathology coordinates and is restricted
to common axis/hinge directions with small cell-type and region amplitude
deviations.  There is deliberately no personal-by-pathology interaction.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import torch
from torch import nn


@dataclass(frozen=True)
class SharedPathologyConfig:
    n_contexts: int
    n_celltypes: int
    n_regions: int
    n_modules: int
    latent_dim: int
    n_pathology: int = 5
    n_pathology_features: int = 3
    personal_rank: int = 2
    hidden_dim: int = 64
    module_token_dim: int = 32
    latent_token_dim: int = 16
    context_embedding_dim: int = 16
    n_heads: int = 4
    n_layers: int = 2
    dropout: float = 0.10
    context_dropout: float = 0.10
    context_gate_scale: float = 0.25
    logvar_min: float = -6.0
    logvar_max: float = 2.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class SharedPathologyAdapter(nn.Module):
    """A/B/R/P/E/S comparison with a frozen-able personal branch."""

    VALID_ARMS = ("A", "B", "R", "P", "E", "S")

    def __init__(
        self,
        config: SharedPathologyConfig,
        *,
        arm: str,
        context_celltype: torch.Tensor,
        context_region: torch.Tensor,
        pathology_feature_center: torch.Tensor,
    ) -> None:
        super().__init__()
        if arm not in self.VALID_ARMS:
            raise ValueError(f"arm must be one of {self.VALID_ARMS}, got {arm!r}")
        if pathology_feature_center.shape != (
            config.n_pathology,
            config.n_pathology_features,
        ):
            raise ValueError("pathology_feature_center has the wrong shape")
        self.config = config
        self.arm = arm
        self.register_buffer("context_celltype", context_celltype.long().clone())
        self.register_buffer("context_region", context_region.long().clone())
        self.register_buffer(
            "pathology_feature_center", pathology_feature_center.float().clone()
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
        self.celltype_embedding = nn.Embedding(
            config.n_celltypes, config.context_embedding_dim
        )
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
        self.personal_module_basis = nn.Parameter(
            torch.empty(config.n_contexts, config.personal_rank, config.n_modules)
        )
        self.personal_latent_basis = nn.Parameter(
            torch.empty(config.n_contexts, config.personal_rank, config.latent_dim)
        )

        # Fixed piecewise-linear features avoid a feature/basis scale ambiguity.
        self.disease_module_basis = nn.Parameter(
            torch.empty(
                config.n_pathology,
                config.n_pathology_features,
                config.n_modules,
            )
        )
        self.disease_latent_basis = nn.Parameter(
            torch.empty(
                config.n_pathology,
                config.n_pathology_features,
                config.latent_dim,
            )
        )
        self.disease_celltype_gate = nn.Parameter(
            torch.empty(
                config.n_celltypes,
                config.n_pathology,
                config.n_pathology_features,
            )
        )
        self.disease_region_gate = nn.Parameter(
            torch.empty(
                config.n_regions,
                config.n_pathology,
                config.n_pathology_features,
            )
        )
        self.reset_parameters()

    @property
    def uses_personal(self) -> bool:
        return self.arm in {"B", "R", "E", "S"}

    @property
    def uses_learned_disease(self) -> bool:
        return self.arm in {"P", "E", "S"}

    def reset_parameters(self) -> None:
        nn.init.normal_(self.personal_module_basis, mean=0.0, std=0.01)
        nn.init.normal_(self.personal_latent_basis, mean=0.0, std=0.01)
        # E/S start exactly at their frozen B prediction; P starts at A.
        nn.init.zeros_(self.disease_module_basis)
        nn.init.zeros_(self.disease_latent_basis)
        nn.init.zeros_(self.disease_celltype_gate)
        nn.init.zeros_(self.disease_region_gate)

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
        token = self.token_projection(
            torch.cat([module_token, latent_token, celltype, region], dim=-1)
        )
        token = token * reliability.unsqueeze(-1)
        encoded = self.context_encoder(token, src_key_padding_mask=~mask)
        denominator = reliability.sum(dim=1, keepdim=True).clamp_min(1e-6)
        pooled = (encoded * reliability.unsqueeze(-1)).sum(dim=1) / denominator
        posterior = self.posterior(pooled)
        mean, logvar = posterior.chunk(2, dim=-1)
        logvar = logvar.clamp(self.config.logvar_min, self.config.logvar_max)
        if sample is None:
            sample = self.training
        code = (
            mean + torch.exp(0.5 * logvar) * torch.randn_like(mean)
            if sample
            else mean
        )
        return code, mean, logvar

    def personal_decode(
        self, code: torch.Tensor, target_context: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        module = torch.einsum(
            "br,brm->bm", code, self.personal_module_basis[target_context]
        )
        latent = torch.einsum(
            "br,brz->bz", code, self.personal_latent_basis[target_context]
        )
        return module, latent

    def _pathology_features(
        self, pathology: torch.Tensor, missing: torch.Tensor
    ) -> torch.Tensor:
        p = pathology.float().clamp(0.0, 1.0)
        features = torch.stack(
            [p, torch.relu(p - 1.0 / 3.0), torch.relu(p - 2.0 / 3.0)],
            dim=-1,
        )
        centered = features - self.pathology_feature_center[None]
        # Missing means "unknown/average", not pathology zero.
        return centered.masked_fill(missing.bool().unsqueeze(-1), 0.0)

    def _context_gate(self, target_context: torch.Tensor) -> torch.Tensor:
        celltype_delta = torch.tanh(self.disease_celltype_gate)
        celltype_delta = celltype_delta - celltype_delta.mean(dim=0, keepdim=True)
        region_delta = torch.tanh(self.disease_region_gate)
        region_delta = region_delta - region_delta.mean(dim=0, keepdim=True)
        celltype = self.context_celltype[target_context]
        region = self.context_region[target_context]
        return 1.0 + self.config.context_gate_scale * (
            celltype_delta[celltype] + region_delta[region]
        )

    def disease_decode(
        self,
        pathology: torch.Tensor,
        missing: torch.Tensor,
        target_context: torch.Tensor,
        *,
        ablate_axis: int | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        features = self._pathology_features(pathology, missing)
        if ablate_axis is not None:
            features = features.clone()
            features[:, int(ablate_axis)] = 0.0
        weighted = features * self._context_gate(target_context)
        axis_module = torch.einsum("bkh,khm->bkm", weighted, self.disease_module_basis)
        axis_latent = torch.einsum("bkh,khz->bkz", weighted, self.disease_latent_basis)
        return (
            axis_module.sum(dim=1),
            axis_latent.sum(dim=1),
            axis_module,
            axis_latent,
        )

    def forward(
        self,
        batch: dict[str, torch.Tensor],
        *,
        sample: bool | None = None,
        ablate_axis: int | None = None,
    ) -> dict[str, torch.Tensor]:
        batch_size = batch["target_context"].shape[0]
        device = batch["target_context"].device
        zero_module = torch.zeros(
            batch_size, self.config.n_modules, device=device, dtype=torch.float32
        )
        zero_latent = torch.zeros(
            batch_size, self.config.latent_dim, device=device, dtype=torch.float32
        )
        zero_code = torch.zeros(
            batch_size, self.config.personal_rank, device=device, dtype=torch.float32
        )

        personal_module = zero_module
        personal_latent = zero_latent
        mean = zero_code
        logvar = zero_code
        if self.uses_personal:
            code, mean, logvar = self.encode(
                batch["source_module"],
                batch["source_latent"],
                batch["source_mask"],
                batch["source_reliability"],
                sample=sample,
            )
            personal_module, personal_latent = self.personal_decode(
                code, batch["target_context"]
            )

        disease_module = zero_module
        disease_latent = zero_latent
        axis_module = zero_module[:, None].expand(-1, self.config.n_pathology, -1)
        axis_latent = zero_latent[:, None].expand(-1, self.config.n_pathology, -1)
        if self.arm == "R":
            disease_module = batch["ridge_disease_module"].float()
            disease_latent = batch["ridge_disease_latent"].float()
        elif self.uses_learned_disease:
            disease_module, disease_latent, axis_module, axis_latent = self.disease_decode(
                batch["pathology"],
                batch["pathology_missing"],
                batch["target_context"],
                ablate_axis=ablate_axis,
            )

        return {
            "module": personal_module + disease_module,
            "latent": personal_latent + disease_latent,
            "personal_module": personal_module,
            "personal_latent": personal_latent,
            "disease_module": disease_module,
            "disease_latent": disease_latent,
            "axis_module": axis_module,
            "axis_latent": axis_latent,
            "posterior_mean": mean,
            "posterior_logvar": logvar,
        }

    def freeze_personal(self) -> None:
        personal_objects: list[nn.Module | nn.Parameter] = [
            self.module_projection,
            self.latent_projection,
            self.celltype_embedding,
            self.region_embedding,
            self.token_projection,
            self.context_encoder,
            self.posterior,
            self.personal_module_basis,
            self.personal_latent_basis,
        ]
        for obj in personal_objects:
            if isinstance(obj, nn.Module):
                for parameter in obj.parameters():
                    parameter.requires_grad_(False)
            else:
                obj.requires_grad_(False)

    def freeze_disease(self) -> None:
        for parameter in (
            self.disease_module_basis,
            self.disease_latent_basis,
            self.disease_celltype_gate,
            self.disease_region_gate,
        ):
            parameter.requires_grad_(False)

    def personal_state_keys(self) -> set[str]:
        prefixes = (
            "module_projection.",
            "latent_projection.",
            "celltype_embedding.",
            "region_embedding.",
            "token_projection.",
            "context_encoder.",
            "posterior.",
            "personal_module_basis",
            "personal_latent_basis",
        )
        return {key for key in self.state_dict() if key.startswith(prefixes)}

    def load_personal_state(self, state: dict[str, torch.Tensor]) -> None:
        expected = self.personal_state_keys()
        available = expected.intersection(state)
        if available != expected:
            raise ValueError(f"personal checkpoint key mismatch: missing={sorted(expected - available)}")
        current = self.state_dict()
        for key in expected:
            if current[key].shape != state[key].shape:
                raise ValueError(f"personal checkpoint shape mismatch for {key}")
            current[key] = state[key]
        self.load_state_dict(current, strict=True)

    def disease_penalty(self) -> torch.Tensor:
        feature_weight = self.disease_module_basis.new_tensor([1.0, 4.0, 4.0])
        module = (self.disease_module_basis.square() * feature_weight[None, :, None]).mean()
        latent = (self.disease_latent_basis.square() * feature_weight[None, :, None]).mean()
        celltype = 4.0 * self.disease_celltype_gate.square().mean()
        region = 8.0 * self.disease_region_gate.square().mean()
        return module + latent + celltype + region

    def branch_parameter_count(self) -> dict[str, int]:
        personal_names = self.personal_state_keys()
        personal = sum(
            parameter.numel()
            for name, parameter in self.named_parameters()
            if name in personal_names
        )
        disease = sum(
            parameter.numel()
            for name, parameter in self.named_parameters()
            if name.startswith("disease_")
        )
        return {
            "personal": int(personal),
            "disease": int(disease),
            "total": int(sum(parameter.numel() for parameter in self.parameters())),
            "trainable": int(
                sum(parameter.numel() for parameter in self.parameters() if parameter.requires_grad)
            ),
        }
