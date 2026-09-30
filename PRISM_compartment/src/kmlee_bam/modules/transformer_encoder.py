from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Literal

import torch
import torch.nn as nn

try:
    from kmlee_bam.modules.transformer_block import TransformerBlock
except ImportError:
    from transformer_block import TransformerBlock

# ==========================================================================
# Output Container
# ==========================================================================
@dataclass
class TransformerEncoderOutput:
    """
    Attributes
    ----------
    hidden_states : [B, T, d_model]
        Final encoder token representations after all layers.
    attn_kl : scalar tensor
        Sum of attention KL terms across layers.
    all_hidden_states : list[[B, T, d_model]] or None
        Hidden states from each layer after residual update.
    all_attn_weights : list[[B, H, T, T]] or None
        Effective attention probabilities for each layer.
    all_posterior_mean_attn : list[[B, H, T, T]] or None
        Deterministic softmax(Phi) for each layer.
    all_prior_mean_attn : list[[B, H, T, T] or None] or None
        Contextual prior attention for each layer.
    all_attn_uncertainty : list[[B, T] or None] or None
        Token-level attention uncertainty for each layer.
    last_attn_uncertainty : [B, T] or None
        Only the last layer's token-level uncertainty.
        Useful when we want cell-level uncertainty without storing
        full attention diagnostics for every layer.
    """
    hidden_states: torch.Tensor
    attn_kl: torch.Tensor
    all_hidden_states: Optional[List[torch.Tensor]]
    all_attn_weights: Optional[List[torch.Tensor]]
    all_posterior_mean_attn: Optional[List[torch.Tensor]]
    all_prior_mean_attn: Optional[List[Optional[torch.Tensor]]]
    all_attn_uncertainty: Optional[List[Optional[torch.Tensor]]]
    last_attn_uncertainty: Optional[torch.Tensor]


# ======================================================================
# Transformer encoder
# ======================================================================
class TransformerEncoder(nn.Module):
    """
    Stack of BAM-based Transformer blocks for the scRNA ordinal encoder.

    Intended role in your model
    ---------------------------
    1) Input token representations come from GeneExpressionEmbedding
    2) Multiple TransformerBlock layers contextualize gene tokens
    3) The final hidden states are passed to state_encoder.py
       to produce posterior parameters (mu_q, logvar_q)

    Notes
    -----
    - Each block already uses pre-norm internally.
    - We optionally apply a final LayerNorm after the full stack.
    - Attention KL is accumulated across layers as a scalar.
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        n_layers: int,
        *,
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
            raise ValueError(f"d_model must be positive, got {d_model}.")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}.")
        if n_layers <= 0:
            raise ValueError(f"n_layers must be positive, got {n_layers}.")
        if d_model % n_heads != 0:
            raise ValueError(
                f"d_model ({d_model}) must be divisible by n_heads ({n_heads})."
            )

        self.d_model = d_model
        self.n_heads = n_heads
        self.n_layers = n_layers
        self.final_norm_enabled = final_norm

        self.layers = nn.ModuleList(
            [
                TransformerBlock(
                    d_model=d_model,
                    n_heads=n_heads,
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
                )
                for _ in range(n_layers)
            ]
        )

        self.final_norm = (
            nn.LayerNorm(d_model, eps=layer_norm_eps) if final_norm else nn.Identity()
        )

    def forward(
        self,
        x: torch.Tensor,
        *,
        attn_mask: Optional[torch.Tensor] = None,
        key_padding_mask: Optional[torch.Tensor] = None,
        n_samples: int = 1,
        return_all_hidden_states: bool = False,
        return_attn_diagnostics: bool = False,
        return_last_attn_uncertainty: bool = False,
    ) -> TransformerEncoderOutput:
        """
         ----------
        x : [B, T, d_model]
            Input token sequence, usually from GeneExpressionEmbedding.
        attn_mask : optional
            Attention mask forwarded to each block.
        key_padding_mask : optional
            Padding mask forwarded to each block.
        n_samples : int
            Number of MC samples for stochastic BAM at eval time.
        return_all_hidden_states : bool
            Whether to keep hidden states from every layer.
        return_attn_diagnostics : bool
            Whether to keep per-layer attention diagnostics.
        return_last_attn_uncertainty : bool
            Keep only the last layer's token-level uncertainty.
            This is a light-weight alternative to full diagnostics.    
        Returns
        -------
        TransformerEncoderOutput
        """
        if x.ndim != 3:
            raise ValueError(f"x must have shape [B, T, d_model], got {tuple(x.shape)}.")
        if x.shape[-1] != self.d_model:
            raise ValueError(
                f"Expected last dim = {self.d_model}, got {x.shape[-1]}."
            )

        hidden_states = x
        total_attn_kl = torch.zeros((), device=x.device, dtype=x.dtype)

        all_hidden_states: Optional[List[torch.Tensor]] = [] if return_all_hidden_states else None
        all_attn_weights: Optional[List[torch.Tensor]] = [] if return_attn_diagnostics else None
        all_posterior_mean_attn: Optional[List[torch.Tensor]] = [] if return_attn_diagnostics else None
        all_prior_mean_attn: Optional[List[Optional[torch.Tensor]]] = [] if return_attn_diagnostics else None
        all_attn_uncertainty: Optional[List[Optional[torch.Tensor]]] = [] if return_attn_diagnostics else None

        last_attn_uncertainty: Optional[torch.Tensor] = None

        for layer in self.layers:
            layer_out = layer(
                hidden_states,
                attn_mask=attn_mask,
                key_padding_mask=key_padding_mask,
                n_samples=n_samples,
            )

            hidden_states = layer_out.hidden_states
            total_attn_kl = total_attn_kl + layer_out.attn_kl

            if return_all_hidden_states:
                assert all_hidden_states is not None
                all_hidden_states.append(hidden_states)

            if return_attn_diagnostics:
                assert all_attn_weights is not None
                assert all_posterior_mean_attn is not None
                assert all_prior_mean_attn is not None
                assert all_attn_uncertainty is not None

                all_attn_weights.append(layer_out.attn_weights)
                all_posterior_mean_attn.append(layer_out.posterior_mean_attn)
                all_prior_mean_attn.append(layer_out.prior_mean_attn)
                all_attn_uncertainty.append(layer_out.attn_uncertainty)

            if return_last_attn_uncertainty:
                last_attn_uncertainty = layer_out.attn_uncertainty

        hidden_states = self.final_norm(hidden_states)

        return TransformerEncoderOutput(
            hidden_states=hidden_states,
            attn_kl=total_attn_kl,
            all_hidden_states=all_hidden_states,
            all_attn_weights=all_attn_weights,
            all_posterior_mean_attn=all_posterior_mean_attn,
            all_prior_mean_attn=all_prior_mean_attn,
            all_attn_uncertainty=all_attn_uncertainty,
            last_attn_uncertainty=last_attn_uncertainty,
        )
    

if __name__ == "__main__":
    torch.manual_seed(42)

    B, T, d_model = 4, 129, 64   # 1 CLS + 128 gene tokens
    n_heads = 8
    n_layers = 3

    x = torch.randn(B, T, d_model)

    key_padding_mask = torch.zeros(B, T, dtype=torch.bool)
    key_padding_mask[0, -5:] = True
    key_padding_mask[1, -3:] = True

    # --------------------------------------------------------------
    # 1. Deterministic / no diagnostics
    # --------------------------------------------------------------
    enc_det = TransformerEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_ff=256,
        stochastic_attention=False,
        final_norm=True,
    )

    out_det = enc_det(
        x,
        key_padding_mask=key_padding_mask,
        return_all_hidden_states=True,
        return_attn_diagnostics=False,
        return_last_attn_uncertainty=False,
    )

    print("[Deterministic | no diagnostics]")
    print("hidden_states.shape         =", out_det.hidden_states.shape)
    print("attn_kl                     =", float(out_det.attn_kl))
    print(
        "n_hidden_layers             =",
        None if out_det.all_hidden_states is None else len(out_det.all_hidden_states),
    )
    print("all_attn_weights            =", out_det.all_attn_weights)
    print("all_posterior_mean_attn     =", out_det.all_posterior_mean_attn)
    print("all_prior_mean_attn         =", out_det.all_prior_mean_attn)
    print("all_attn_uncertainty        =", out_det.all_attn_uncertainty)
    print("last_attn_uncertainty       =", out_det.last_attn_uncertainty)

    assert out_det.hidden_states.shape == (B, T, d_model)
    assert out_det.attn_kl.ndim == 0
    assert out_det.all_hidden_states is not None
    assert len(out_det.all_hidden_states) == n_layers
    assert out_det.all_attn_weights is None
    assert out_det.all_posterior_mean_attn is None
    assert out_det.all_prior_mean_attn is None
    assert out_det.all_attn_uncertainty is None
    assert out_det.last_attn_uncertainty is None

    # --------------------------------------------------------------
    # 2. Stochastic / uncertainty only
    # --------------------------------------------------------------
    enc_sto_light = TransformerEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_ff=256,
        stochastic_attention=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
        final_norm=True,
    )

    out_sto_light = enc_sto_light(
        x,
        key_padding_mask=key_padding_mask,
        return_all_hidden_states=False,
        return_attn_diagnostics=False,
        return_last_attn_uncertainty=True,
    )

    print("\n[Stochastic | uncertainty only]")
    print("hidden_states.shape         =", out_sto_light.hidden_states.shape)
    print(
        "attn_kl.shape               =",
        "scalar" if out_sto_light.attn_kl.ndim == 0 else tuple(out_sto_light.attn_kl.shape),
    )
    print("all_attn_weights            =", out_sto_light.all_attn_weights)
    print("all_posterior_mean_attn     =", out_sto_light.all_posterior_mean_attn)
    print("all_prior_mean_attn         =", out_sto_light.all_prior_mean_attn)
    print("all_attn_uncertainty        =", out_sto_light.all_attn_uncertainty)
    print(
        "last_attn_uncertainty.shape =",
        None if out_sto_light.last_attn_uncertainty is None
        else out_sto_light.last_attn_uncertainty.shape,
    )

    assert out_sto_light.hidden_states.shape == (B, T, d_model)
    assert out_sto_light.attn_kl.ndim == 0
    assert out_sto_light.all_attn_weights is None
    assert out_sto_light.all_posterior_mean_attn is None
    assert out_sto_light.all_prior_mean_attn is None
    assert out_sto_light.all_attn_uncertainty is None
    assert out_sto_light.last_attn_uncertainty is not None
    assert out_sto_light.last_attn_uncertainty.shape == (B, T)

    # --------------------------------------------------------------
    # 3. Stochastic / full diagnostics
    # --------------------------------------------------------------
    enc_sto_full = TransformerEncoder(
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_ff=256,
        stochastic_attention=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
        final_norm=True,
    )

    out_sto_full = enc_sto_full(
        x,
        key_padding_mask=key_padding_mask,
        return_all_hidden_states=True,
        return_attn_diagnostics=True,
        return_last_attn_uncertainty=False,
    )

    print("\n[Stochastic | full diagnostics]")
    print("hidden_states.shape         =", out_sto_full.hidden_states.shape)
    print(
        "attn_kl.shape               =",
        "scalar" if out_sto_full.attn_kl.ndim == 0 else tuple(out_sto_full.attn_kl.shape),
    )
    print(
        "n_hidden_layers             =",
        None if out_sto_full.all_hidden_states is None else len(out_sto_full.all_hidden_states),
    )
    print(
        "n_attn_layers               =",
        None if out_sto_full.all_attn_weights is None else len(out_sto_full.all_attn_weights),
    )
    print(
        "layer0 attn_weights.shape   =",
        None if out_sto_full.all_attn_weights is None else out_sto_full.all_attn_weights[0].shape,
    )
    print(
        "layer0 prior_mean.shape     =",
        None
        if out_sto_full.all_prior_mean_attn is None or out_sto_full.all_prior_mean_attn[0] is None
        else out_sto_full.all_prior_mean_attn[0].shape,
    )
    print(
        "layer_last unc.shape        =",
        None
        if out_sto_full.all_attn_uncertainty is None or out_sto_full.all_attn_uncertainty[-1] is None
        else out_sto_full.all_attn_uncertainty[-1].shape,
    )
    print("last_attn_uncertainty       =", out_sto_full.last_attn_uncertainty)

    assert out_sto_full.hidden_states.shape == (B, T, d_model)
    assert out_sto_full.attn_kl.ndim == 0
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
    assert out_sto_full.last_attn_uncertainty is None

    print("\n[OK] TransformerEncoder smoke test passed.")