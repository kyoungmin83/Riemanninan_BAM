from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Literal

import torch
import torch.nn as nn
import torch.nn.functional as F

from kmlee_bam.modules.stochastic_attention import BayesianAttention


# ======================================================================
# Output container
# ======================================================================
@dataclass
class TransformerBlockOutput:
    """
    Attributes
    ----------
    hidden_states : [B, T, d_model]
        Output token representations after attention + FFN.
    attn_kl : scalar tensor
        Attention KL from BAM (0 for deterministic attention).
    attn_weights : [B, H, T, T]
        Sampled / effective attention probabilities from BAM, detached.
    posterior_mean_attn : [B, H, T, T]
        Deterministic softmax(Phi), detached.
    prior_mean_attn : [B, H, T, T] or None
        Prior attention probabilities from contextual prior, detached.
    attn_uncertainty : [B, T] or None
        Mean attention entropy over heads from BAM.
    """
    hidden_states: torch.Tensor
    attn_kl: torch.Tensor
    attn_weights: torch.Tensor
    posterior_mean_attn: torch.Tensor
    prior_mean_attn: Optional[torch.Tensor]
    attn_uncertainty: Optional[torch.Tensor]


# ======================================================================
# Feed-forward sublayer
# ======================================================================
class FeedForward(nn.Module):
    """
    Standard transformer FFN:

        x -> Linear(d_model -> d_ff) -> GELU -> Dropout
          -> Linear(d_ff -> d_model) -> Dropout
    """

    def __init__(
        self,
        d_model: int,
        d_ff: int,
        *,
        dropout: float = 0.1,
        init_std: float = 0.02,
    ) -> None:
        super().__init__()

        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}.")
        if d_ff <= 0:
            raise ValueError(f"d_ff must be positive, got {d_ff}.")
        if not (0.0 <= dropout < 1.0):
            raise ValueError(f"dropout must be in [0, 1), got {dropout}.")

        self.fc1 = nn.Linear(d_model, d_ff)
        self.fc2 = nn.Linear(d_ff, d_model)
        self.dropout = nn.Dropout(dropout)
        self.init_std = init_std

        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.fc1.weight, std=self.init_std)
        nn.init.zeros_(self.fc1.bias)
        nn.init.normal_(self.fc2.weight, std=self.init_std)
        nn.init.zeros_(self.fc2.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.fc1(x)
        x = F.gelu(x)
        x = self.dropout(x)
        x = self.fc2(x)
        x = self.dropout(x)
        return x


# ======================================================================
# Transformer block for ordinal scRNA encoder
# ======================================================================
class TransformerBlock(nn.Module):
    """
    Pre-norm Transformer block using BayesianAttention (BAM).

    Structure
    ---------
        x
        -> LN
        -> BAM self-attention
        -> residual add
        -> LN
        -> FFN
        -> residual add

    Notes
    -----
    1) This block is designed for the encoder side of your ordinal scRNA model.
    2) BAM already handles:
         - deterministic / stochastic attention
         - contextual prior
         - attention KL
         - attention uncertainty
    3) The block simply wraps BAM into a stable pre-norm residual block.
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
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
    ) -> None:
        super().__init__()

        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}.")
        if n_heads <= 0:
            raise ValueError(f"n_heads must be positive, got {n_heads}.")
        if d_model % n_heads != 0:
            raise ValueError(
                f"d_model ({d_model}) must be divisible by n_heads ({n_heads})."
            )

        if d_ff is None:
            d_ff = 4 * d_model

        self.d_model = d_model
        self.n_heads = n_heads
        self.d_ff = d_ff

        # Pre-norm before attention
        self.norm_attn = nn.LayerNorm(d_model, eps=layer_norm_eps)

        # BAM attention core
        self.attn = BayesianAttention(
            d_model=d_model,
            n_heads=n_heads,
            stochastic=stochastic_attention,
            distribution=distribution,
            sigma_mode=sigma_mode,
            sigma=sigma,
            sigma_min=sigma_min,
            sigma_max=sigma_max,
            weibull_k=weibull_k,
            prior_d_mid=prior_d_mid,
            attn_dropout=attn_dropout,
            proj_dropout=proj_dropout,
            qkv_bias=qkv_bias,
            init_std=init_std,
        )

        # Pre-norm before FFN
        self.norm_ffn = nn.LayerNorm(d_model, eps=layer_norm_eps)

        # Position-wise FFN
        self.ffn = FeedForward(
            d_model=d_model,
            d_ff=d_ff,
            dropout=ffn_dropout,
            init_std=init_std,
        )

    def forward(
        self,
        x: torch.Tensor,
        *,
        attn_mask: Optional[torch.Tensor] = None,
        key_padding_mask: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> TransformerBlockOutput:
        """
        Parameters
        ----------
        x : [B, T, d_model]
            Token representations.
        attn_mask : optional
            Attention mask passed through to BAM.
        key_padding_mask : optional
            Padding mask passed through to BAM.
        n_samples : int
            MC samples for stochastic attention at eval time.

        Returns
        -------
        TransformerBlockOutput
        """
        if x.ndim != 3:
            raise ValueError(f"x must have shape [B, T, d_model], got {tuple(x.shape)}.")
        if x.shape[-1] != self.d_model:
            raise ValueError(
                f"Expected last dim = {self.d_model}, got {x.shape[-1]}."
            )

        # --------------------------------------------------------------
        # 1) Pre-norm attention sublayer
        # --------------------------------------------------------------
        x_attn_in = self.norm_attn(x)

        attn_out = self.attn(
            x_attn_in,
            attn_mask=attn_mask,
            key_padding_mask=key_padding_mask,
            n_samples=n_samples,
        )

        # Residual add
        x = x + attn_out.output

        # --------------------------------------------------------------
        # 2) Pre-norm FFN sublayer
        # --------------------------------------------------------------
        x_ffn_in = self.norm_ffn(x)
        x = x + self.ffn(x_ffn_in)

        return TransformerBlockOutput(
            hidden_states=x,
            attn_kl=attn_out.attn_kl,
            attn_weights=attn_out.attn_weights,
            posterior_mean_attn=attn_out.posterior_mean_attn,
            prior_mean_attn=attn_out.prior_mean_attn,
            attn_uncertainty=attn_out.attn_uncertainty,
        )


if __name__ == "__main__":
    torch.manual_seed(42)

    B, T, d_model, n_heads = 4, 128, 64, 8

    # deterministic block
    block_det = TransformerBlock(
        d_model=d_model,
        n_heads=n_heads,
        d_ff=256,
        stochastic_attention=False,
    )

    x = torch.randn(B, T, d_model)
    out_det = block_det(x)

    print("[Deterministic]")
    print("hidden_states.shape      =", out_det.hidden_states.shape)
    print("attn_kl                  =", float(out_det.attn_kl))
    print("attn_weights.shape       =", out_det.attn_weights.shape)
    print("posterior_mean_attn.shape=", out_det.posterior_mean_attn.shape)
    print("prior_mean_attn          =", out_det.prior_mean_attn)
    print("attn_uncertainty         =", out_det.attn_uncertainty)

    # stochastic block
    block_sto = TransformerBlock(
        d_model=d_model,
        n_heads=n_heads,
        d_ff=256,
        stochastic_attention=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
    )

    out_sto = block_sto(x)

    print("\n[Stochastic]")
    print("hidden_states.shape      =", out_sto.hidden_states.shape)
    print("attn_kl.shape            =", tuple(out_sto.attn_kl.shape) if out_sto.attn_kl.ndim > 0 else "scalar")
    print("attn_weights.shape       =", out_sto.attn_weights.shape)
    print("posterior_mean_attn.shape=", out_sto.posterior_mean_attn.shape)
    print("prior_mean_attn.shape    =", None if out_sto.prior_mean_attn is None else out_sto.prior_mean_attn.shape)
    print("attn_uncertainty.shape   =", None if out_sto.attn_uncertainty is None else out_sto.attn_uncertainty.shape)

    assert out_det.hidden_states.shape == (B, T, d_model)
    assert out_sto.hidden_states.shape == (B, T, d_model)
    assert out_det.attn_weights.shape == (B, n_heads, T, T)
    assert out_sto.attn_weights.shape == (B, n_heads, T, T)

    print("\n[OK] TransformerBlock smoke test passed.")