"""
Module · BayesianAttention (BAM) — v3
=====================================
Unified deterministic / Bayesian multi-head self-attention.

    1) Prior calibration:
       mu_prior = scale * contextual_prior(K) + bias

    2) Probability vs dropout separation:
       W_prob = true attention distribution
       W_drop = dropout(W_prob) used only for value aggregation

    3) Consistent uncertainty:
       always defined as mean attention entropy over heads,
       computed from W_prob (never from dropout-corrupted weights)

    4) Safe masked softmax:
       supports all-invalid rows safely by returning zero rows

    5) Sigma upper clamp:
       prevents sigma explosion in sigma_mode="logit"

Mathematical Foundation
-----------------------
Deterministic attention:
    W = softmax(Phi),   Phi = Q K^T / sqrt(d_k)

BAM stochastic attention:
    S_{i,j} ~ Distribution(params from Phi_{i,j})
    W_{i,:} = S_{i,:} / sum_j S_{i,j}

Lognormal:
    S ~ LN(mu, sigma^2),   mu = Phi - sigma^2 / 2

Weibull:
    S ~ Weibull(k, lambda),   lambda = exp(Phi) / Gamma(1 + 1/k)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Literal

import torch
import torch.nn as nn
import torch.nn.functional as F


# ======================================================================
# Output container
# ======================================================================
@dataclass
class BAMAttentionOutput:
    """
    Attributes
    ----------
    output : [B, T, d_model]
        Hidden states after attention + output projection.
    attn_kl : scalar tensor
        KL regularisation term (0 when deterministic).
    attn_weights : [B, H, T_q, T_k]
        Proper attention probabilities W_prob (detached).
    posterior_mean_attn : [B, H, T_q, T_k]
        Deterministic softmax(Phi), detached.
    prior_mean_attn : [B, H, T_q, T_k] or None
        softmax(mu_prior), detached.
    attn_uncertainty : [B, T_q] or None
        Mean attention entropy over heads, computed from W_prob.
    """
    output: torch.Tensor
    attn_kl: torch.Tensor
    attn_weights: torch.Tensor
    posterior_mean_attn: torch.Tensor
    prior_mean_attn: Optional[torch.Tensor]
    attn_uncertainty: Optional[torch.Tensor]


# ======================================================================
# Contextual Prior
# ======================================================================
class ContextualPrior(nn.Module):
    """
    Key-dependent prior: K -> per-key log-importance score.

        psi = log softmax(F2(ReLU(F1(K))))

    Output shape:
        [B, H, 1, T_k]
    """

    def __init__(self, d_k: int, d_mid: int = 16) -> None:
        super().__init__()
        self.f1 = nn.Linear(d_k, d_mid)
        self.f2 = nn.Linear(d_mid, 1)
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.xavier_uniform_(self.f1.weight)
        nn.init.zeros_(self.f1.bias)
        nn.init.xavier_uniform_(self.f2.weight)
        nn.init.zeros_(self.f2.bias)

    def forward(self, K: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        K : [B, H, T_k, d_k]

        Returns
        -------
        psi : [B, H, 1, T_k]
            Raw prior log-importance before calibration.
        """
        h = F.relu(self.f1(K))             # [B, H, T_k, d_mid]
        logits = self.f2(h).squeeze(-1)    # [B, H, T_k]
        psi = F.log_softmax(logits, dim=-1)
        return psi.unsqueeze(-2)           # [B, H, 1, T_k]


# ======================================================================
# Main module
# ======================================================================
class BayesianAttention(nn.Module):
    """
    Unified deterministic / Bayesian multi-head self-attention.

    Parameters
    ----------
    d_model : int
        Token embedding dimension.
    n_heads : int
        Number of attention heads.
    stochastic : bool
        False -> standard attention
        True  -> BAM stochastic attention
    distribution : {"lognormal", "weibull"}
        Distribution family for unnormalised scores S.
    sigma_mode : {"global", "logit"}
        "global": one scalar sigma
        "logit" : sigma_ij = sigma_base + softplus(w * Phi_ij + b), then clamped
    sigma : float
        Base sigma.
    sigma_min : float
        Lower bound for sigma.
    sigma_max : float
        Upper bound for sigma.
    weibull_k : float
        Weibull shape parameter.
    prior_d_mid : int
        Hidden dim for contextual prior MLP.
    attn_dropout : float
        Dropout on attention weights (applied only to W_drop).
    proj_dropout : float
        Dropout after output projection.
    qkv_bias : bool
        Whether Q/K/V projections use bias.
    init_std : float
        Std for normal init of projection weights.
    """

    def __init__(
        self,
        d_model: int,
        n_heads: int,
        *,
        stochastic: bool = False,
        distribution: Literal["lognormal", "weibull"] = "lognormal",
        sigma_mode: Literal["global", "logit"] = "global",
        sigma: float = 0.3,
        sigma_min: float = 1e-4,
        sigma_max: float = 1.5,
        weibull_k: float = 20.0,
        prior_d_mid: int = 16,
        attn_dropout: float = 0.1,
        proj_dropout: float = 0.1,
        qkv_bias: bool = False,
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
        if distribution not in ("lognormal", "weibull"):
            raise ValueError(
                f"distribution must be 'lognormal' or 'weibull', got '{distribution}'."
            )
        if sigma_mode not in ("global", "logit"):
            raise ValueError(
                f"sigma_mode must be 'global' or 'logit', got '{sigma_mode}'."
            )
        if sigma <= 0:
            raise ValueError(f"sigma must be positive, got {sigma}.")
        if sigma_min <= 0:
            raise ValueError(f"sigma_min must be positive, got {sigma_min}.")
        if sigma_max <= sigma_min:
            raise ValueError(
                f"sigma_max must be > sigma_min, got sigma_min={sigma_min}, sigma_max={sigma_max}."
            )

        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_model // n_heads
        self.scale = 1.0 / math.sqrt(self.d_k)

        self.stochastic = stochastic
        self.distribution = distribution
        self.sigma_mode = sigma_mode

        self.sigma_base = float(sigma)
        self.sigma_min = float(sigma_min)
        self.sigma_max = float(sigma_max)

        self.weibull_k = weibull_k
        self.init_std = init_std

        # projections
        self.W_q = nn.Linear(d_model, d_model, bias=qkv_bias)
        self.W_k = nn.Linear(d_model, d_model, bias=qkv_bias)
        self.W_v = nn.Linear(d_model, d_model, bias=qkv_bias)
        self.W_o = nn.Linear(d_model, d_model, bias=True)

        self.attn_dropout = nn.Dropout(attn_dropout)
        self.proj_dropout = nn.Dropout(proj_dropout)

        # contextual prior
        self.contextual_prior = ContextualPrior(self.d_k, prior_d_mid)

        # [PATCH-1] prior calibration: positive scale + free bias
        # softplus(0.5413248546) ≈ 1.0
        self.prior_scale_raw = nn.Parameter(torch.tensor(0.5413248546))
        self.prior_bias = nn.Parameter(torch.tensor(0.0))

        # sigma from logits
        if sigma_mode == "logit":
            self.sigma_logit_w = nn.Parameter(torch.tensor(0.0))
            self.sigma_logit_b = nn.Parameter(torch.tensor(-1.0))
        else:
            self.sigma_logit_w = None
            self.sigma_logit_b = None

        # Weibull constant
        self._log_gamma_1pk_inv: Optional[float] = None
        if distribution == "weibull":
            self._log_gamma_1pk_inv = -math.lgamma(1.0 + 1.0 / weibull_k)

        self.reset_parameters()

    def reset_parameters(self) -> None:
        for linear in (self.W_q, self.W_k, self.W_v):
            nn.init.normal_(linear.weight, std=self.init_std)
            if linear.bias is not None:
                nn.init.zeros_(linear.bias)

        nn.init.normal_(self.W_o.weight, std=self.init_std)
        nn.init.zeros_(self.W_o.bias)

    # ==================================================================
    # Forward
    # ==================================================================
    def forward(
        self,
        x: torch.Tensor,
        attn_mask: Optional[torch.Tensor] = None,
        key_padding_mask: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> BAMAttentionOutput:
        """
        Parameters
        ----------
        x : [B, T, d_model]
        attn_mask : [T,T] | [B,T,T] | [B,1,T,T] | [B,H,T,T], optional
            Boolean (True = blocked) or additive (< -1e8 = blocked).
        key_padding_mask : [B, T], optional
            Boolean, True = padding key.
        n_samples : int
            Number of MC samples used only in eval stochastic mode.
        """
        if x.ndim != 3:
            raise ValueError(f"x must have shape [B, T, d_model], got {tuple(x.shape)}.")

        B, T, D = x.shape
        if D != self.d_model:
            raise ValueError(f"Expected d_model={self.d_model}, got {D}.")

        Q = self._to_heads(self.W_q(x))   # [B, H, T, d_k]
        K = self._to_heads(self.W_k(x))   # [B, H, T, d_k]
        V = self._to_heads(self.W_v(x))   # [B, H, T, d_k]

        Phi = torch.matmul(Q, K.transpose(-2, -1)) * self.scale  # [B,H,T,T]

        valid_mask = self._build_valid_mask(
            B=B,
            H=self.n_heads,
            T_q=T,
            T_k=T,
            attn_mask=attn_mask,
            key_padding_mask=key_padding_mask,
            device=x.device,
        )  # [B, H, T, T], True = allowed

        posterior_mean_attn = self._masked_softmax(Phi, valid_mask).detach()

        if not self.stochastic:
            return self._forward_deterministic(
                Phi=Phi,
                V=V,
                valid_mask=valid_mask,
                posterior_mean_attn=posterior_mean_attn,
            )

        if self.training or n_samples <= 1:
            return self._forward_stochastic(
                Phi=Phi,
                K=K,
                V=V,
                valid_mask=valid_mask,
                posterior_mean_attn=posterior_mean_attn,
            )

        return self._forward_mc(
            Phi=Phi,
            K=K,
            V=V,
            valid_mask=valid_mask,
            posterior_mean_attn=posterior_mean_attn,
            n_samples=n_samples,
        )

    # ==================================================================
    # Deterministic path
    # ==================================================================
    def _forward_deterministic(
        self,
        Phi: torch.Tensor,
        V: torch.Tensor,
        valid_mask: torch.Tensor,
        posterior_mean_attn: torch.Tensor,
    ) -> BAMAttentionOutput:
        W_prob = self._masked_softmax(Phi, valid_mask)
        W_drop = self.attn_dropout(W_prob)  # [PATCH-2]

        out = self._project_output(torch.matmul(W_drop, V))

        return BAMAttentionOutput(
            output=out,
            attn_kl=torch.tensor(0.0, device=Phi.device, dtype=Phi.dtype),
            attn_weights=W_prob.detach(),
            posterior_mean_attn=posterior_mean_attn,
            prior_mean_attn=None,
            attn_uncertainty=None,
        )

    # ==================================================================
    # Stochastic path (training / single-sample)
    # ==================================================================
    def _forward_stochastic(
        self,
        Phi: torch.Tensor,
        K: torch.Tensor,
        V: torch.Tensor,
        valid_mask: torch.Tensor,
        posterior_mean_attn: torch.Tensor,
    ) -> BAMAttentionOutput:
        S, mu_q, sigma_q = self._sample(Phi)

        W_prob = self._normalize_positive(S, valid_mask)   # [PATCH-2]
        W_drop = self.attn_dropout(W_prob)

        out = self._project_output(torch.matmul(W_drop, V))

        raw_prior = self.contextual_prior(K)               # [B,H,1,T]
        mu_prior = self._calibrate_prior(raw_prior)        # [PATCH-1]
        sigma_prior = self._base_sigma()

        kl = self._compute_kl(
            mu_q=mu_q,
            sigma_q=sigma_q,
            mu_prior=mu_prior,
            sigma_prior=sigma_prior,
            valid_mask=valid_mask,
        )

        prior_mean_attn = self._masked_softmax(
            mu_prior.expand_as(Phi), valid_mask
        ).detach()

        uncertainty = self._entropy(W_prob)                # [PATCH-3]

        return BAMAttentionOutput(
            output=out,
            attn_kl=kl,
            attn_weights=W_prob.detach(),
            posterior_mean_attn=posterior_mean_attn,
            prior_mean_attn=prior_mean_attn,
            attn_uncertainty=uncertainty,
        )

    # ==================================================================
    # MC path (eval)
    # ==================================================================
    def _forward_mc(
        self,
        Phi: torch.Tensor,
        K: torch.Tensor,
        V: torch.Tensor,
        valid_mask: torch.Tensor,
        posterior_mean_attn: torch.Tensor,
        n_samples: int,
    ) -> BAMAttentionOutput:
        outputs = []
        weight_samples = []

        for _ in range(n_samples):
            S, _, _ = self._sample(Phi)
            W_prob = self._normalize_positive(S, valid_mask)
            outputs.append(torch.matmul(W_prob, V))
            weight_samples.append(W_prob)

        out = self._project_output(torch.stack(outputs, dim=0).mean(dim=0))
        W_mean = torch.stack(weight_samples, dim=0).mean(dim=0)

        raw_prior = self.contextual_prior(K)
        mu_prior = self._calibrate_prior(raw_prior)
        prior_mean_attn = self._masked_softmax(
            mu_prior.expand_as(Phi), valid_mask
        ).detach()

        # [PATCH-3] uncertainty definition unified:
        # always mean attention entropy over heads, computed from a proper probability tensor
        uncertainty = self._entropy(W_mean)

        return BAMAttentionOutput(
            output=out,
            attn_kl=torch.tensor(0.0, device=Phi.device, dtype=Phi.dtype),
            attn_weights=W_mean.detach(),
            posterior_mean_attn=posterior_mean_attn,
            prior_mean_attn=prior_mean_attn,
            attn_uncertainty=uncertainty,
        )

    # ==================================================================
    # Sampling
    # ==================================================================
    def _sample(
        self,
        Phi: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | float]:
        """
        Sample unnormalised scores S.

        Returns
        -------
        S : [B,H,Tq,Tk]
        mu_q : [B,H,Tq,Tk]
        sigma_q : float or [B,H,Tq,Tk]
        """
        sigma_q = self._get_sigma(Phi)

        if self.distribution == "lognormal":
            if isinstance(sigma_q, float):
                mu_q = Phi - 0.5 * (sigma_q ** 2)
            else:
                mu_q = Phi - 0.5 * sigma_q.pow(2)

            eps = torch.randn_like(Phi)
            S = torch.exp(mu_q + sigma_q * eps)
            return S, mu_q, sigma_q

        # Weibull
        k = self.weibull_k
        log_lambda = Phi + self._log_gamma_1pk_inv
        eps = torch.rand_like(Phi).clamp(1e-7, 1.0 - 1e-7)
        S = torch.exp(
            log_lambda + (1.0 / k) * torch.log(-torch.log(1.0 - eps))
        )
        return S, log_lambda, sigma_q

    def _get_sigma(self, Phi: torch.Tensor) -> torch.Tensor | float:
        """
        "global": scalar sigma
        "logit" : sigma_base + softplus(w * Phi + b), then clamped
        """
        if self.sigma_mode == "global":
            return self._base_sigma()

        raw = self.sigma_logit_w * Phi + self.sigma_logit_b
        sigma_q = self.sigma_base + F.softplus(raw)
        sigma_q = sigma_q.clamp(min=self.sigma_min, max=self.sigma_max)  # [PATCH-5]
        return sigma_q

    def _base_sigma(self) -> float:
        return float(min(max(self.sigma_base, self.sigma_min), self.sigma_max))

    # ==================================================================
    # KL divergence
    # ==================================================================
    def _compute_kl(
        self,
        mu_q: torch.Tensor,
        sigma_q,                 # float or Tensor
        mu_prior: torch.Tensor,  # [B,H,1,Tk]
        sigma_prior: float,
        valid_mask: torch.Tensor,
    ) -> torch.Tensor:
        if self.distribution == "lognormal":
            return self._kl_lognormal(
                mu_q=mu_q,
                sigma_q=sigma_q,
                mu_prior=mu_prior,
                sigma_prior=sigma_prior,
                valid_mask=valid_mask,
            )
        return self._kl_weibull_gamma(
            log_lambda_q=mu_q,
            log_alpha_prior=mu_prior,
            valid_mask=valid_mask,
        )

    def _kl_lognormal(
        self,
        mu_q: torch.Tensor,
        sigma_q,                  # float or Tensor
        mu_prior: torch.Tensor,
        sigma_prior: float,
        valid_mask: torch.Tensor,
    ) -> torch.Tensor:
        """
        Analytic KL(LN(mu_q, sigma_q^2) || LN(mu_p, sigma_p^2)).

        General form:
            log(sigma_p / sigma_q)
            + (sigma_q^2 + (mu_q - mu_p)^2) / (2 sigma_p^2)
            - 1/2

        If sigma_q == sigma_p, simplifies to:
            (mu_q - mu_p)^2 / (2 sigma^2)
        """
        diff = mu_q - mu_prior  # broadcast over query dimension
        sigma_p = sigma_prior

        if self.sigma_mode == "global":
            kl_elem = diff.pow(2) / (2.0 * sigma_p ** 2)
        else:
            sigma_q_t = sigma_q.clamp(min=self.sigma_min, max=self.sigma_max)
            sigma_q_sq = sigma_q_t.pow(2)
            sigma_p_sq = sigma_p ** 2

            sigma_p_tensor = torch.tensor(
                sigma_p,
                device=mu_q.device,
                dtype=mu_q.dtype,
            )

            kl_elem = (
                torch.log(sigma_p_tensor / sigma_q_t.clamp(min=1e-8))
                + (sigma_q_sq + diff.pow(2)) / (2.0 * sigma_p_sq)
                - 0.5
            )

        kl_elem = kl_elem * valid_mask.float()
        n_valid = valid_mask.float().sum().clamp(min=1.0)
        return kl_elem.sum() / n_valid

    def _kl_weibull_gamma(
        self,
        log_lambda_q: torch.Tensor,
        log_alpha_prior: torch.Tensor,
        valid_mask: torch.Tensor,
    ) -> torch.Tensor:
        """
        KL(Weibull(k, lambda) || Gamma(alpha, beta=1))
        """
        EULER = 0.5772156649
        k = self.weibull_k
        beta = 1.0

        alpha = torch.exp(log_alpha_prior).clamp(min=1e-6)
        lam = torch.exp(log_lambda_q)
        gamma_1pk = math.exp(-self._log_gamma_1pk_inv)

        kl_elem = (
            EULER * alpha / k
            - alpha * log_lambda_q
            + math.log(k)
            + beta * lam * gamma_1pk
            - EULER
            - 1.0
            - alpha * math.log(beta)
            + torch.lgamma(alpha)
        )

        kl_elem = kl_elem * valid_mask.float()
        n_valid = valid_mask.float().sum().clamp(min=1.0)
        return kl_elem.sum() / n_valid

    # ==================================================================
    # Probability helpers
    # ==================================================================
    def _masked_softmax(
        self,
        logits: torch.Tensor,
        valid_mask: torch.Tensor,
        dim: int = -1,
    ) -> torch.Tensor:
        """
        [PATCH-4] Safe masked softmax.
        - computes in float32 for numerical stability
        - avoids fp16 overflow from large negative fill values
        - invalid positions become 0
        - all-invalid rows become all-zero rows
        - returns exact zero rows when a whole row is invalid
        """
        if logits.shape != valid_mask.shape:
            raise ValueError(
                f"logits and valid_mask must have the same shape, "
                f"got {tuple(logits.shape)} vs {tuble(valid_mask.shape)}."
            )

        #orig_dtype = logits.dtype
        #logits_fp32 = logits.float()
        mask = valid_mask.bool()

        # Using small value in fp32 in order to prevent overflow
        net_large = torch.finfo(logits.dtype).min

        masked_logits = logits.masked_fill(~mask, net_large)
        probs = F.softmax(masked_logits, dim=dim)
        probs = probs * mask.to(probs.dtype)

        denom = probs.sum(dim=dim, keepdim=True)
      

        # if a row is entirely invalid, keep it exactly zero
        probs = torch.where(
            denom > 0,
            probs / denom.clamp_min(1e-8),
            torch.zeros_like(probs),
        )
        return probs

    def _normalize_positive(
        self,
        S: torch.Tensor,
        valid_mask: torch.Tensor,
    ) -> torch.Tensor:
        """
        Safe normalisation for positive stochastic scores.

        Invalid entries are zeroed before normalisation.
        All-invalid rows become all-zero rows.
        """
        S_masked = S * valid_mask.float()
        denom = S_masked.sum(dim=-1, keepdim=True)
        W = S_masked / denom.clamp_min(1e-12)
        W = torch.where(
            denom > 0,
            W,
            torch.zeros_like(W),
        )
        return W

    # ==================================================================
    # Prior helpers
    # ==================================================================
    def _prior_scale(self) -> torch.Tensor:
        # positive scale
        return F.softplus(self.prior_scale_raw) + 1e-4

    def _calibrate_prior(self, raw_prior: torch.Tensor) -> torch.Tensor:
        """
        [PATCH-1] Prior calibration:
            mu_prior = scale * raw_prior + bias
        """
        return self._prior_scale() * raw_prior + self.prior_bias

    # ==================================================================
    # Utility methods
    # ==================================================================
    def _project_output(self, context: torch.Tensor) -> torch.Tensor:
        out = self._from_heads(context)
        out = self.W_o(out)
        out = self.proj_dropout(out)
        return out

    def _entropy(self, W_prob: torch.Tensor) -> torch.Tensor:
        """
        [PATCH-3] Consistent uncertainty:
        mean attention entropy over heads, computed from proper probabilities.
        """
        log_W = torch.log(W_prob.clamp(min=1e-12))
        ent = -(W_prob * log_W).sum(dim=-1)   # [B,H,Tq]
        return ent.mean(dim=1)                # [B,Tq]

    def _to_heads(self, x: torch.Tensor) -> torch.Tensor:
        B, T, _ = x.shape
        return x.view(B, T, self.n_heads, self.d_k).transpose(1, 2)

    def _from_heads(self, x: torch.Tensor) -> torch.Tensor:
        B, H, T, dk = x.shape
        return x.transpose(1, 2).contiguous().view(B, T, H * dk)

    # ==================================================================
    # Mask infrastructure
    # ==================================================================
    def _build_valid_mask(
        self,
        *,
        B: int,
        H: int,
        T_q: int,
        T_k: int,
        attn_mask: Optional[torch.Tensor],
        key_padding_mask: Optional[torch.Tensor],
        device: torch.device,
    ) -> torch.Tensor:
        """
        Build boolean mask: True = attention allowed
        """
        valid = torch.ones(B, H, T_q, T_k, dtype=torch.bool, device=device)

        if attn_mask is not None:
            if attn_mask.dtype == torch.bool:
                blocked = attn_mask
            else:
                blocked = attn_mask < -1e8

            if blocked.ndim == 2:
                valid = valid & ~blocked.unsqueeze(0).unsqueeze(0)
            elif blocked.ndim == 3:
                valid = valid & ~blocked.unsqueeze(1)
            elif blocked.ndim == 4:
                valid = valid & ~blocked
            else:
                raise ValueError(
                    f"attn_mask must be 2D/3D/4D, got ndim={blocked.ndim}."
                )

        if key_padding_mask is not None:
            if key_padding_mask.shape != (B, T_k):
                raise ValueError(
                    f"key_padding_mask shape must be ({B}, {T_k}), got {tuple(key_padding_mask.shape)}."
                )
            if key_padding_mask.dtype != torch.bool:
                raise TypeError("key_padding_mask must be bool.")
            valid = valid & ~key_padding_mask[:, None, None, :]

        return valid

    # ==================================================================
    # Mode switching
    # ==================================================================
    def set_stochastic(self, mode: bool) -> None:
        self.stochastic = mode

    def extra_repr(self) -> str:
        return (
            f"d_model={self.d_model}, n_heads={self.n_heads}, "
            f"d_k={self.d_k}, stochastic={self.stochastic}, "
            f"dist={self.distribution}, sigma_mode={self.sigma_mode}, "
            f"sigma_base={self.sigma_base}, sigma_min={self.sigma_min}, "
            f"sigma_max={self.sigma_max}, weibull_k={self.weibull_k}"
        )


# ======================================================================
# Smoke test
# ======================================================================
if __name__ == "__main__":
    torch.manual_seed(42)

    B, T, D, H = 2, 10, 32, 2
    x = torch.randn(B, T, D)

    kpm = torch.zeros(B, T, dtype=torch.bool)
    kpm[0, -2:] = True

    # 1) Deterministic
    attn = BayesianAttention(d_model=D, n_heads=H, stochastic=False)
    out = attn(x, key_padding_mask=kpm)
    print("=== Deterministic ===")
    print("output      :", out.output.shape)
    print("kl          :", out.attn_kl.item())
    print("uncertainty :", out.attn_uncertainty)

    loss = out.output.sum()
    loss.backward()
    print("grad flows  :", attn.W_q.weight.grad is not None)

    row_sums = out.attn_weights.sum(dim=-1)
    print("simplex err :", (row_sums - 1.0).abs().max().item())

    # 2) Stochastic lognormal, global sigma
    attn_ln = BayesianAttention(
        d_model=D,
        n_heads=H,
        stochastic=True,
        distribution="lognormal",
        sigma_mode="global",
        sigma=0.3,
    )
    out_ln = attn_ln(x, key_padding_mask=kpm)
    print("\n=== Stochastic LN (global sigma) ===")
    print("output      :", out_ln.output.shape)
    print("kl          :", out_ln.attn_kl.item())
    print("uncertainty :", out_ln.attn_uncertainty.shape)
    print("prior_mean  :", out_ln.prior_mean_attn.shape)

    (out_ln.output.sum() + out_ln.attn_kl).backward()
    print("prior grad  :", attn_ln.contextual_prior.f1.weight.grad is not None)

    print("masked attn :", out_ln.attn_weights[0, :, :, -2:].abs().max().item())

    # 3) Stochastic lognormal, logit sigma
    attn_ld = BayesianAttention(
        d_model=D,
        n_heads=H,
        stochastic=True,
        distribution="lognormal",
        sigma_mode="logit",
        sigma=0.1,
        sigma_max=1.2,
    )
    out_ld = attn_ld(x, key_padding_mask=kpm)
    print("\n=== Stochastic LN (logit sigma) ===")
    print("output      :", out_ld.output.shape)
    print("kl          :", out_ld.attn_kl.item())

    (out_ld.output.sum() + out_ld.attn_kl).backward()
    print("sigma_w grad:", attn_ld.sigma_logit_w.grad is not None)

    # 4) Stochastic Weibull
    attn_wb = BayesianAttention(
        d_model=D,
        n_heads=H,
        stochastic=True,
        distribution="weibull",
        weibull_k=20.0,
    )
    out_wb = attn_wb(x)
    print("\n=== Stochastic Weibull ===")
    print("output      :", out_wb.output.shape)
    print("kl          :", out_wb.attn_kl.item())

    (out_wb.output.sum() + out_wb.attn_kl).backward()
    print("grad flows  :", attn_wb.W_q.weight.grad is not None)

    # 5) Mode switch
    attn_ln.set_stochastic(False)
    out_sw = attn_ln(x)
    print("\n=== After set_stochastic(False) ===")
    print("kl          :", out_sw.attn_kl.item())

    # 6) MC eval
    attn_ln.set_stochastic(True)
    attn_ln.eval()
    with torch.no_grad():
        out_mc = attn_ln(x, n_samples=5)
    print("\n=== MC eval (n_samples=5) ===")
    print("output      :", out_mc.output.shape)
    print("uncertainty :", out_mc.attn_uncertainty.shape)

    print("\ndone.")