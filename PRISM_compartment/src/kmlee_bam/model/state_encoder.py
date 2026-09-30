from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Literal

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from kmlee_bam.modules.transformer_encoder import TransformerEncoder
    from kmlee_bam.modules.attention_guided_pooling import AttentionGuidedPooling
    from kmlee_bam.modules.persistent_attention_slots import PersistentAttentionSlots
except ImportError:
    from transformer_encoder import TransformerEncoder
    from attention_guided_pooling import AttentionGuidedPooling
    from persistent_attention_slots import PersistentAttentionSlots


# ======================================================================
# Output container
# ======================================================================
@dataclass
class StateEncoderOutput:
    hidden_states: torch.Tensor
    pooled_state: torch.Tensor
    posterior_input: torch.Tensor
    mu_q: torch.Tensor
    logvar_q: torch.Tensor
    z_s: torch.Tensor
    attn_kl: torch.Tensor
    cell_uncertainty: Optional[torch.Tensor]
    all_hidden_states: Optional[List[torch.Tensor]]
    all_attn_weights: Optional[List[torch.Tensor]]
    all_posterior_mean_attn: Optional[List[torch.Tensor]]
    all_prior_mean_attn: Optional[List[Optional[torch.Tensor]]]
    all_attn_uncertainty: Optional[List[Optional[torch.Tensor]]]
    tech_risk_logit: Optional[torch.Tensor] = None
    pooling_attention: Optional[torch.Tensor] = None
    pooling_auxiliary: Optional[Dict[str, torch.Tensor]] = None


# ======================================================================
# Posterior head
# ======================================================================
class PosteriorHead(nn.Module):
    """
    Shared MLP trunk -> separate heads for mu_q and logvar_q.
    """

    def __init__(
        self,
        in_dim: int,
        d_z: int,
        *,
        hidden_dim: Optional[int] = None,
        dropout: float = 0.1,
        init_std: float = 0.02,
    ) -> None:
        super().__init__()

        if in_dim <= 0:
            raise ValueError("in_dim must be positive.")
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if hidden_dim is None:
            hidden_dim = max(in_dim, 2 * d_z)

        self.in_dim = in_dim
        self.d_z = d_z
        self.hidden_dim = hidden_dim
        self.init_std = init_std

        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.dropout = nn.Dropout(dropout)
        self.mu_head = nn.Linear(hidden_dim, d_z)
        self.logvar_head = nn.Linear(hidden_dim, d_z)

        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.fc1.weight, std=self.init_std)
        nn.init.zeros_(self.fc1.bias)

        nn.init.normal_(self.mu_head.weight, std=self.init_std)
        nn.init.zeros_(self.mu_head.bias)

        nn.init.normal_(self.logvar_head.weight, std=self.init_std)
        nn.init.zeros_(self.logvar_head.bias)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.fc1(x)
        h = F.gelu(h)
        h = self.dropout(h)
        mu_q = self.mu_head(h)
        logvar_q = self.logvar_head(h)
        return mu_q, logvar_q


# ======================================================================
# State encoder
# ======================================================================
class StateEncoder(nn.Module):
    """
    BAM-based encoder + pooling + posterior head for

        q_phi(z_s | y_i, c_i)
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        n_layers: int,
        d_z: int,
        *,
        # posterior conditioning
        condition_on_celltype: bool = True,
        n_celltypes: Optional[int] = None,
        celltype_embed_dim: Optional[int] = None,
        posterior_input_norm: Optional[bool] = None,
        # pooling
        pooling: Literal["cls", "mean", "agp", "persistent_agp"] = "cls",
        use_cls_token: bool = True,
        mean_exclude_cls: bool = True,
        agp_heads: int = 8,
        agp_exclude_cls: bool = True,
        agp_variant: str = "legacy",
        agp_temperature_init: float = 0.75,
        agp_temperature_min: float = 0.25,
        agp_temperature_max: float = 2.0,
        agp_mean_residual_init: float = 0.25,
        agp_mean_residual_max: float = 0.50,
        agp_diversity_weight: float = 0.0,
        agp_diversity_margin: float = 0.50,
        agp_query_orthogonality_weight: float = 0.0,
        agp_entropy_band_weight: float = 0.0,
        agp_min_effective_tokens: float = 1.0,
        agp_max_effective_tokens: float = 1.0e9,
        persistent_agp_attention_heads: int = 4,
        persistent_agp_ffn_dim: Optional[int] = None,
        persistent_agp_update_gate_init: float = 0.10,
        persistent_agp_update_gate_max: float = 1.0,
        # uncertainty pooling
        compute_cell_uncertainty: bool = True,
        uncertainty_exclude_cls: bool = True,
        # tech-risk head (supervised technical-OOD detector; off by default)
        tech_risk_head: bool = False,
        # posterior head
        posterior_hidden_dim: Optional[int] = None,
        posterior_dropout: float = 0.1,
        logvar_min: float = -8.0,
        logvar_max: float = 8.0,
        # transformer encoder
        d_ff: Optional[int] = None,
        stochastic_attention: bool = False,
        distribution: Literal["lognormal", "weibull"] = "lognormal",
        sigma_mode: Literal["global", "logit"] = "global",
        sigma: float = 0.3,
        sigma_min: float = 1e-4,
        sigma_max: float = 1.5,
        weibull_k: float = 20.0,
        prior_d_mid: int = 16,
        attn_dropout: float = 0.1,
        proj_dropout: float = 0.1,
        ffn_dropout: float = 0.1,
        qkv_bias: bool = False,
        layer_norm_eps: float = 1e-5,
        init_std: float = 0.02,
        final_norm: bool = True,
    ) -> None:
        super().__init__()

        if d_model <= 0:
            raise ValueError("d_model must be positive.")
        if n_heads <= 0:
            raise ValueError("n_heads must be positive.")
        if n_layers <= 0:
            raise ValueError("n_layers must be positive.")
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if pooling not in {"cls", "mean", "agp", "persistent_agp"}:
            raise ValueError(
                "pooling must be one of 'cls', 'mean', 'agp', or "
                "'persistent_agp'."
            )
        if logvar_max <= logvar_min:
            raise ValueError("logvar_max must be greater than logvar_min.")

        if condition_on_celltype and n_celltypes is None:
            raise ValueError(
                "n_celltypes must be provided when condition_on_celltype=True."
            )

        if celltype_embed_dim is None:
            celltype_embed_dim = max(16, d_model // 2)

        if posterior_input_norm is None:
            posterior_input_norm = condition_on_celltype

        self.d_model = d_model
        self.d_z = d_z
        self.pooling = pooling
        self.use_cls_token = use_cls_token
        self.mean_exclude_cls = mean_exclude_cls
        self.agp_exclude_cls = bool(agp_exclude_cls)
        self.compute_cell_uncertainty = compute_cell_uncertainty
        self.uncertainty_exclude_cls = uncertainty_exclude_cls
        self.condition_on_celltype = condition_on_celltype
        self.logvar_min = logvar_min
        self.logvar_max = logvar_max

        self.encoder = TransformerEncoder(
            d_model=d_model,
            n_heads=n_heads,
            n_layers=n_layers,
            d_ff=d_ff,
            stochastic_attention=stochastic_attention,
            distribution=distribution,
            sigma_mode=sigma_mode,
            sigma=sigma,
            sigma_min=sigma_min,
            sigma_max=sigma_max,
            weibull_k=weibull_k,
            prior_d_mid=prior_d_mid,
            attn_dropout=attn_dropout,
            proj_dropout=proj_dropout,
            ffn_dropout=ffn_dropout,
            qkv_bias=qkv_bias,
            layer_norm_eps=layer_norm_eps,
            init_std=init_std,
            final_norm=final_norm,
        )

        self.attention_pool = (
            AttentionGuidedPooling(
                d_model=d_model,
                n_heads=int(agp_heads),
                variant=agp_variant,
                temperature_init=agp_temperature_init,
                temperature_min=agp_temperature_min,
                temperature_max=agp_temperature_max,
                mean_residual_init=agp_mean_residual_init,
                mean_residual_max=agp_mean_residual_max,
                diversity_weight=agp_diversity_weight,
                diversity_margin=agp_diversity_margin,
                query_orthogonality_weight=agp_query_orthogonality_weight,
                entropy_band_weight=agp_entropy_band_weight,
                min_effective_tokens=agp_min_effective_tokens,
                max_effective_tokens=agp_max_effective_tokens,
                init_std=init_std,
            )
            if pooling == "agp"
            else None
        )
        self.persistent_attention_pool = (
            PersistentAttentionSlots(
                d_model=d_model,
                n_slots=int(agp_heads),
                n_attention_heads=int(persistent_agp_attention_heads),
                n_stages=int(n_layers) + 1,
                ffn_dim=persistent_agp_ffn_dim,
                dropout=ffn_dropout,
                init_gate=persistent_agp_update_gate_init,
                gate_max=persistent_agp_update_gate_max,
                mean_residual_init=agp_mean_residual_init,
                mean_residual_max=agp_mean_residual_max,
                diversity_weight=agp_diversity_weight,
                diversity_margin=agp_diversity_margin,
                query_orthogonality_weight=agp_query_orthogonality_weight,
                entropy_band_weight=agp_entropy_band_weight,
                min_effective_tokens=agp_min_effective_tokens,
                max_effective_tokens=agp_max_effective_tokens,
                layer_norm_eps=layer_norm_eps,
            )
            if pooling == "persistent_agp"
            else None
        )

        if condition_on_celltype:
            self.celltype_embedding = nn.Embedding(n_celltypes, celltype_embed_dim)
            nn.init.normal_(self.celltype_embedding.weight, mean=0.0, std=init_std)
            posterior_in_dim = d_model + celltype_embed_dim
        else:
            self.celltype_embedding = None
            posterior_in_dim = d_model

        self.posterior_input_norm = (
            nn.LayerNorm(posterior_in_dim, eps=layer_norm_eps)
            if posterior_input_norm
            else nn.Identity()
        )

        self.posterior_head = PosteriorHead(
            in_dim=posterior_in_dim,
            d_z=d_z,
            hidden_dim=posterior_hidden_dim,
            dropout=posterior_dropout,
            init_std=init_std,
        )

        # Supervised technical-OOD risk head (single-forward amortization of the MC-epistemic
        # firewall signal). Reads the pooled cell state; trained by tech_risk_bce in the trainer.
        # Off by default => no new params, existing checkpoints load unchanged.
        if tech_risk_head:
            _hid = max(16, d_model // 2)
            self.tech_risk_head = nn.Sequential(
                nn.Linear(d_model, _hid), nn.GELU(), nn.Linear(_hid, 1)
            )
        else:
            self.tech_risk_head = None

    # ==================================================================
    # Main forward
    # ==================================================================
    def forward(
        self,
        x_tokens: torch.Tensor,
        *,
        celltype_id: Optional[torch.Tensor] = None,
        attn_mask: Optional[torch.Tensor] = None,
        key_padding_mask: Optional[torch.Tensor] = None,
        n_samples: int = 1,
        sample_latent: bool = True,
        return_all_hidden_states: bool = False,
        return_attn_diagnostics: bool = False,
    ) -> StateEncoderOutput:
        if x_tokens.ndim != 3:
            raise ValueError(
                f"x_tokens must have shape [B, T, d_model], got {tuple(x_tokens.shape)}."
            )
        if x_tokens.shape[-1] != self.d_model:
            raise ValueError(
                f"Expected last dim = {self.d_model}, got {x_tokens.shape[-1]}."
            )

        B, T, _ = x_tokens.shape
        self._validate_key_padding_mask(key_padding_mask, B=B, T=T)

        if self.condition_on_celltype:
            self._validate_celltype_id(celltype_id, B=B)

        need_full_attn_diagnostics = return_attn_diagnostics
        need_last_attn_uncertainty = self.compute_cell_uncertainty and (not return_attn_diagnostics)

        requested_all_hidden_states = bool(return_all_hidden_states)
        encoder_needs_all_hidden_states = (
            requested_all_hidden_states or self.pooling == "persistent_agp"
        )
        enc_out = self.encoder(
            x_tokens,
            attn_mask=attn_mask,
            key_padding_mask=key_padding_mask,
            n_samples=n_samples,
            return_all_hidden_states=encoder_needs_all_hidden_states,
            return_attn_diagnostics=need_full_attn_diagnostics,
            return_last_attn_uncertainty=need_last_attn_uncertainty,
        )

        pooling_stage_attention = None
        if self.pooling == "persistent_agp":
            if self.persistent_attention_pool is None:
                raise RuntimeError(
                    "pooling='persistent_agp' but persistent_attention_pool is missing"
                )
            if enc_out.all_hidden_states is None:
                raise RuntimeError(
                    "persistent AGP requires all Transformer hidden states"
                )
            transformer_stages = list(enc_out.all_hidden_states)
            # TransformerEncoder records per-layer states before its final
            # LayerNorm.  Use the normalized final output for the last read so
            # the established final-norm parameters remain part of the actual
            # pooled computation (and therefore of DDP reduction).
            transformer_stages[-1] = enc_out.hidden_states
            stage_states = [x_tokens, *transformer_stages]
            stripped_stages = []
            persistent_mask = key_padding_mask
            for stage_state in stage_states:
                stripped, stage_mask = self._maybe_strip_cls(
                    stage_state,
                    key_padding_mask=key_padding_mask,
                    exclude_cls=self.use_cls_token and self.agp_exclude_cls,
                )
                stripped_stages.append(stripped)
                persistent_mask = stage_mask
            pooled_state, pooling_attention, pooling_stage_attention = (
                self.persistent_attention_pool(
                    stripped_stages,
                    key_padding_mask=persistent_mask,
                )
            )
        else:
            pooled_state, pooling_attention = self._pool_hidden_states(
                enc_out.hidden_states,
                key_padding_mask=key_padding_mask,
            )
        pooling_auxiliary = None
        if pooling_attention is not None:
            if self.pooling == "persistent_agp":
                assert self.persistent_attention_pool is not None
                pooling_auxiliary = self.persistent_attention_pool.auxiliary_terms(
                    pooling_attention,
                    stage_weights=pooling_stage_attention,
                )
            else:
                if self.attention_pool is None:
                    raise RuntimeError("pooling attention exists but AGP module is absent")
                pooling_auxiliary = self.attention_pool.auxiliary_terms(
                    pooling_attention
                )

        posterior_input = pooled_state
        if self.condition_on_celltype:
            assert celltype_id is not None
            cell_h = self.celltype_embedding(celltype_id)
            posterior_input = torch.cat([pooled_state, cell_h], dim=-1)

        posterior_input = self.posterior_input_norm(posterior_input)

        mu_q, logvar_q = self.posterior_head(posterior_input)
        logvar_q = logvar_q.clamp(min=self.logvar_min, max=self.logvar_max)

        z_s = self.reparameterize(mu_q, logvar_q, sample=sample_latent)

        cell_uncertainty = None
        if self.compute_cell_uncertainty:
            last_unc = None

            if return_attn_diagnostics:
                if enc_out.all_attn_uncertainty is not None and len(enc_out.all_attn_uncertainty) > 0:
                    last_unc = enc_out.all_attn_uncertainty[-1]
            else:
                last_unc = enc_out.last_attn_uncertainty

            if last_unc is not None:
                cell_uncertainty = self._pool_scalar_map(
                    last_unc,
                    key_padding_mask=key_padding_mask,
                    exclude_cls=self.uncertainty_exclude_cls,
                )

        if not return_attn_diagnostics:
            all_attn_weights = None
            all_posterior_mean_attn = None
            all_prior_mean_attn = None
            all_attn_uncertainty = None
        else:
            all_attn_weights = enc_out.all_attn_weights
            all_posterior_mean_attn = enc_out.all_posterior_mean_attn
            all_prior_mean_attn = enc_out.all_prior_mean_attn
            all_attn_uncertainty = enc_out.all_attn_uncertainty

        tech_risk_logit = None
        if self.tech_risk_head is not None:
            tech_risk_logit = self.tech_risk_head(pooled_state).squeeze(-1)

        return StateEncoderOutput(
            hidden_states=enc_out.hidden_states,
            pooled_state=pooled_state,
            posterior_input=posterior_input,
            mu_q=mu_q,
            logvar_q=logvar_q,
            z_s=z_s,
            attn_kl=enc_out.attn_kl,
            cell_uncertainty=cell_uncertainty,
            tech_risk_logit=tech_risk_logit,
            all_hidden_states=(
                enc_out.all_hidden_states if requested_all_hidden_states else None
            ),
            all_attn_weights=all_attn_weights,
            all_posterior_mean_attn=all_posterior_mean_attn,
            all_prior_mean_attn=all_prior_mean_attn,
            all_attn_uncertainty=all_attn_uncertainty,
            pooling_attention=pooling_attention,
            pooling_auxiliary=pooling_auxiliary,
        )

    # ==================================================================
    # Validation
    # ==================================================================
    @staticmethod
    def _validate_key_padding_mask(
        key_padding_mask: Optional[torch.Tensor],
        *,
        B: int,
        T: int,
    ) -> None:
        if key_padding_mask is None:
            return
        if key_padding_mask.ndim != 2 or key_padding_mask.shape != (B, T):
            raise ValueError(
                f"key_padding_mask must have shape {(B, T)}, got {tuple(key_padding_mask.shape)}."
            )
        if key_padding_mask.dtype != torch.bool:
            raise TypeError(
                f"key_padding_mask must be torch.bool, got {key_padding_mask.dtype}."
            )

    @staticmethod
    def _validate_celltype_id(
        celltype_id: Optional[torch.Tensor],
        *,
        B: int,
    ) -> None:
        if celltype_id is None:
            raise ValueError("celltype_id must be provided when condition_on_celltype=True.")
        if celltype_id.dtype != torch.long:
            raise TypeError("celltype_id must be torch.long.")
        if celltype_id.ndim != 1 or celltype_id.shape[0] != B:
            raise ValueError(
                f"celltype_id must have shape {(B,)}, got {tuple(celltype_id.shape)}."
            )

    # ==================================================================
    # Pooling
    # ==================================================================
    def _pool_hidden_states(
        self,
        hidden_states: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor],
    ) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
        if self.pooling == "cls":
            if not self.use_cls_token:
                raise ValueError(
                    "pooling='cls' but use_cls_token=False. "
                    "Either enable CLS in gene_embedding or switch pooling to 'mean'."
                )
            return hidden_states[:, 0, :], None

        if self.pooling == "agp":
            token_states, token_mask = self._maybe_strip_cls(
                hidden_states,
                key_padding_mask=key_padding_mask,
                exclude_cls=self.use_cls_token and self.agp_exclude_cls,
            )
            if token_states.shape[1] == 0:
                raise ValueError("No tokens remain after CLS stripping for AGP.")
            if self.attention_pool is None:
                raise RuntimeError("pooling='agp' but attention_pool is missing")
            return self.attention_pool(
                token_states,
                key_padding_mask=token_mask,
            )

        token_states, token_mask = self._maybe_strip_cls(
            hidden_states,
            key_padding_mask=key_padding_mask,
            exclude_cls=self.use_cls_token and self.mean_exclude_cls,
        )

        if token_states.shape[1] == 0:
            raise ValueError("No tokens remain after CLS stripping for mean pooling.")

        if token_mask is None:
            return token_states.mean(dim=1), None

        valid = (~token_mask).float().unsqueeze(-1)
        denom = valid.sum(dim=1).clamp_min(1.0)
        return (token_states * valid).sum(dim=1) / denom, None

    def _pool_scalar_map(
        self,
        token_values: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor],
        exclude_cls: bool,
    ) -> torch.Tensor:
        vals, mask = self._maybe_strip_cls(
            token_values,
            key_padding_mask=key_padding_mask,
            exclude_cls=self.use_cls_token and exclude_cls,
        )

        if vals.shape[1] == 0:
            raise ValueError("No tokens remain after CLS stripping for uncertainty pooling.")

        if mask is None:
            return vals.mean(dim=1)

        valid = (~mask).float()
        denom = valid.sum(dim=1).clamp_min(1.0)
        return (vals * valid).sum(dim=1) / denom

    @staticmethod
    def _maybe_strip_cls(
        x: torch.Tensor,
        *,
        key_padding_mask: Optional[torch.Tensor],
        exclude_cls: bool,
    ) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
        if exclude_cls:
            x = x[:, 1:]
            if key_padding_mask is not None:
                key_padding_mask = key_padding_mask[:, 1:]
        return x, key_padding_mask

    # ==================================================================
    # Reparameterization
    # ==================================================================
    @staticmethod
    def reparameterize(
        mu_q: torch.Tensor,
        logvar_q: torch.Tensor,
        *,
        sample: bool = True,
    ) -> torch.Tensor:
        if not sample:
            return mu_q
        std_q = torch.exp(0.5 * logvar_q)
        eps = torch.randn_like(std_q)
        return mu_q + std_q * eps

if __name__ == "__main__":
    torch.manual_seed(42)

    B = 4
    T = 129          # 1 CLS + 128 gene tokens
    d_model = 64
    n_heads = 8
    n_layers = 3
    d_z = 32
    n_celltypes = 24

    x_tokens = torch.randn(B, T, d_model)
    celltype_id = torch.randint(0, n_celltypes, (B,), dtype=torch.long)

    key_padding_mask = torch.zeros(B, T, dtype=torch.bool)
    key_padding_mask[0, -5:] = True
    key_padding_mask[1, -3:] = True

    # --------------------------------------------------------------
    # 1. Deterministic / no full diagnostics / cell_uncertainty ON
    # --------------------------------------------------------------
    enc_det = StateEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_z=d_z,
        n_celltypes=n_celltypes,
        condition_on_celltype=True,
        pooling="cls",
        use_cls_token=True,
        compute_cell_uncertainty=True,
        stochastic_attention=False,
    )

    out_det = enc_det(
        x_tokens,
        celltype_id=celltype_id,
        key_padding_mask=key_padding_mask,
        sample_latent=True,
        return_all_hidden_states=True,
        return_attn_diagnostics=False,
    )

    print("[Deterministic | no full diagnostics]")
    print("hidden_states.shape      =", out_det.hidden_states.shape)
    print("pooled_state.shape       =", out_det.pooled_state.shape)
    print("posterior_input.shape    =", out_det.posterior_input.shape)
    print("mu_q.shape               =", out_det.mu_q.shape)
    print("logvar_q.shape           =", out_det.logvar_q.shape)
    print("z_s.shape                =", out_det.z_s.shape)
    print("attn_kl                  =", float(out_det.attn_kl))
    print("cell_uncertainty         =", out_det.cell_uncertainty)
    print("all_attn_weights         =", out_det.all_attn_weights)
    print("all_posterior_mean_attn  =", out_det.all_posterior_mean_attn)
    print("all_prior_mean_attn      =", out_det.all_prior_mean_attn)
    print("all_attn_uncertainty     =", out_det.all_attn_uncertainty)

    assert out_det.hidden_states.shape == (B, T, d_model)
    assert out_det.pooled_state.shape == (B, d_model)
    assert out_det.mu_q.shape == (B, d_z)
    assert out_det.logvar_q.shape == (B, d_z)
    assert out_det.z_s.shape == (B, d_z)
    assert out_det.attn_kl.ndim == 0
    assert out_det.all_attn_weights is None
    assert out_det.all_posterior_mean_attn is None
    assert out_det.all_prior_mean_attn is None
    assert out_det.all_attn_uncertainty is None

    # deterministic BAM에서는 uncertainty가 None이어야 자연스럽다
    assert out_det.cell_uncertainty is None

    # --------------------------------------------------------------
    # 2. Stochastic / no full diagnostics / lightweight uncertainty path
    # --------------------------------------------------------------
    enc_sto_light = StateEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_z=d_z,
        n_celltypes=n_celltypes,
        condition_on_celltype=True,
        pooling="cls",
        use_cls_token=True,
        compute_cell_uncertainty=True,
        stochastic_attention=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
    )

    out_sto_light = enc_sto_light(
        x_tokens,
        celltype_id=celltype_id,
        key_padding_mask=key_padding_mask,
        sample_latent=True,
        return_all_hidden_states=False,
        return_attn_diagnostics=False,
    )

    print("\n[Stochastic | lightweight uncertainty path]")
    print("hidden_states.shape      =", out_sto_light.hidden_states.shape)
    print("pooled_state.shape       =", out_sto_light.pooled_state.shape)
    print("posterior_input.shape    =", out_sto_light.posterior_input.shape)
    print("mu_q.shape               =", out_sto_light.mu_q.shape)
    print("logvar_q.shape           =", out_sto_light.logvar_q.shape)
    print("z_s.shape                =", out_sto_light.z_s.shape)
    print(
        "attn_kl.shape            =",
        "scalar" if out_sto_light.attn_kl.ndim == 0 else tuple(out_sto_light.attn_kl.shape),
    )
    print(
        "cell_uncertainty.shape   =",
        None if out_sto_light.cell_uncertainty is None else out_sto_light.cell_uncertainty.shape,
    )
    print("all_attn_weights         =", out_sto_light.all_attn_weights)
    print("all_posterior_mean_attn  =", out_sto_light.all_posterior_mean_attn)
    print("all_prior_mean_attn      =", out_sto_light.all_prior_mean_attn)
    print("all_attn_uncertainty     =", out_sto_light.all_attn_uncertainty)

    assert out_sto_light.hidden_states.shape == (B, T, d_model)
    assert out_sto_light.pooled_state.shape == (B, d_model)
    assert out_sto_light.mu_q.shape == (B, d_z)
    assert out_sto_light.logvar_q.shape == (B, d_z)
    assert out_sto_light.z_s.shape == (B, d_z)
    assert out_sto_light.attn_kl.ndim == 0
    assert out_sto_light.cell_uncertainty is not None
    assert out_sto_light.cell_uncertainty.shape == (B,)
    assert out_sto_light.all_attn_weights is None
    assert out_sto_light.all_posterior_mean_attn is None
    assert out_sto_light.all_prior_mean_attn is None
    assert out_sto_light.all_attn_uncertainty is None

    # --------------------------------------------------------------
    # 3. Stochastic / full diagnostics
    # --------------------------------------------------------------
    enc_sto_full = StateEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_z=d_z,
        n_celltypes=n_celltypes,
        condition_on_celltype=True,
        pooling="cls",
        use_cls_token=True,
        compute_cell_uncertainty=True,
        stochastic_attention=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
    )

    out_sto_full = enc_sto_full(
        x_tokens,
        celltype_id=celltype_id,
        key_padding_mask=key_padding_mask,
        sample_latent=True,
        return_all_hidden_states=True,
        return_attn_diagnostics=True,
    )

    print("\n[Stochastic | full diagnostics]")
    print("hidden_states.shape      =", out_sto_full.hidden_states.shape)
    print("pooled_state.shape       =", out_sto_full.pooled_state.shape)
    print("posterior_input.shape    =", out_sto_full.posterior_input.shape)
    print("mu_q.shape               =", out_sto_full.mu_q.shape)
    print("logvar_q.shape           =", out_sto_full.logvar_q.shape)
    print("z_s.shape                =", out_sto_full.z_s.shape)
    print(
        "attn_kl.shape            =",
        "scalar" if out_sto_full.attn_kl.ndim == 0 else tuple(out_sto_full.attn_kl.shape),
    )
    print(
        "cell_uncertainty.shape   =",
        None if out_sto_full.cell_uncertainty is None else out_sto_full.cell_uncertainty.shape,
    )
    print(
        "n_hidden_layers          =",
        None if out_sto_full.all_hidden_states is None else len(out_sto_full.all_hidden_states),
    )
    print(
        "n_attn_layers            =",
        None if out_sto_full.all_attn_weights is None else len(out_sto_full.all_attn_weights),
    )
    print(
        "layer0 attn.shape        =",
        None if out_sto_full.all_attn_weights is None else out_sto_full.all_attn_weights[0].shape,
    )
    print(
        "layer_last unc.shape     =",
        None
        if out_sto_full.all_attn_uncertainty is None or out_sto_full.all_attn_uncertainty[-1] is None
        else out_sto_full.all_attn_uncertainty[-1].shape,
    )

    assert out_sto_full.hidden_states.shape == (B, T, d_model)
    assert out_sto_full.pooled_state.shape == (B, d_model)
    assert out_sto_full.mu_q.shape == (B, d_z)
    assert out_sto_full.logvar_q.shape == (B, d_z)
    assert out_sto_full.z_s.shape == (B, d_z)
    assert out_sto_full.attn_kl.ndim == 0
    assert out_sto_full.cell_uncertainty is not None
    assert out_sto_full.cell_uncertainty.shape == (B,)
    assert out_sto_full.all_hidden_states is not None
    assert len(out_sto_full.all_hidden_states) == n_layers
    assert out_sto_full.all_attn_weights is not None
    assert len(out_sto_full.all_attn_weights) == n_layers
    assert out_sto_full.all_attn_weights[0].shape == (B, n_heads, T, T)
    assert out_sto_full.all_posterior_mean_attn is not None
    assert len(out_sto_full.all_posterior_mean_attn) == n_layers
    assert out_sto_full.all_prior_mean_attn is not None
    assert len(out_sto_full.all_prior_mean_attn) == n_layers
    assert out_sto_full.all_attn_uncertainty is not None
    assert len(out_sto_full.all_attn_uncertainty) == n_layers
    assert out_sto_full.all_attn_uncertainty[-1] is not None
    assert out_sto_full.all_attn_uncertainty[-1].shape == (B, T)

    print("\n[OK] StateEncoder smoke test passed.")
