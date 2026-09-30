from __future__ import annotations

import math
from typing import Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


def inverse_softplus(y: float) -> float:
    """
    Solve softplus(x) = y for x.
    Numerically safer for small positive y.
    """
    y = float(y)
    if y <= 0:
        raise ValueError("inverse_softplus requires y > 0.")
    # softplux(x) = log(1 + exp(x))
    # x = log(exp(y) -1)
    return math.log(math.expm1(y))


class CellTypePrior(nn.Module):
    """
    Cell-type-conditional diagonal Gaussian prior

        p(z_s | c) = N(mu_p(c), diag(sigma_p(c)^2))

    Also provides standardized residualization

        z_perp = (z_s - mu_p(c)) / (sigma_p(c) + eps)
    """

    def __init__(
        self,
        n_celltypes: int,
        d_z: int,
        sigma_init: float = 0.8,
        sigma_min: float = 1e-4,
        residual_eps: float = 1e-5,
        residual_clamp: float = 10.0,
    ) -> None:
        super().__init__()

        if n_celltypes <= 0:
            raise ValueError("n_celltypes must be positive.")
        if d_z <= 0:
            raise ValueError("d_z must be positive.")
        if sigma_init <= 0:
            raise ValueError("sigma_init must be positive.")
        if sigma_min <= 0:
            raise ValueError("sigma_min must be positive.")
        if residual_eps <= 0:
            raise ValueError("residual_eps must be positive.")
        if residual_clamp <= 0:
            raise ValueError("residual_clamp must be positive.")
        
        self.n_celltypes = n_celltypes
        self.d_z = d_z
        self.sigma_min = sigma_min
        self.residual_eps = residual_eps
        self.residual_clamp = residual_clamp

        # Prior mean: one vector per cell type
        self.mu_embedding = nn.Embedding(n_celltypes ,d_z) # mu_p(c_i)
        nn.init.zeros_(self.mu_embedding.weight)

        # Prior raw scale -> softplus -> sigma_p > 0
        raw_init = inverse_softplus(sigma_init)
        # Prior sigma per cell type
        self.raw_sigma_embedding = nn.Embedding(n_celltypes, d_z)
        nn.init.constant_(self.raw_sigma_embedding.weight, raw_init)

    def forward(
        self,
        celltype_id: torch.Tensor,# c_i
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Returns
        -------
        mu_p     : [B, d_z]
        logvar_p : [B, d_z]
        sigma_p  : [B, d_z]
        """
        if celltype_id.dtype != torch.long:
            raise TypeError("celltype_id must be torch.long")

        mu_p = self.mu_embedding(celltype_id) #[B, d_z]

        sigma_p = F.softplus(self.raw_sigma_embedding(celltype_id))
        sigma_p = sigma_p.clamp(min=self.sigma_min)

        logvar_p = 2.0 * torch.log(sigma_p)

        return mu_p, logvar_p, sigma_p

    def residualize(
        self,
        z_s: torch.Tensor,         # [B, d_z]
        celltype_id: torch.Tensor, # [B], long
    ) -> torch.Tensor:
        """
        z_perp = (z_s - mu_p(c)) / (sigma_p(c) + eps)
        """
        mu_p, _, sigma_p = self.forward(celltype_id)
        z_perp = (z_s - mu_p) / (sigma_p + self.residual_eps)
        return z_perp.clamp(-self.residual_clamp, self.residual_clamp)
    
    @staticmethod
    def kl_divergence(
        mu_q: torch.Tensor,     # [B, d_z]
        logvar_q: torch.Tensor, # [B, d_z]
        mu_p: torch.Tensor,     # [B, d_z]
        logvar_p: torch.Tensor, # [B, d_z]
    ) -> torch.Tensor:
        """
       State KL[q || p] for diagonal Gaussians, summed over latent dims.
        Returns [B].
        """
        var_q = torch.exp(logvar_q)
        var_p = torch.exp(logvar_p)

        kl = 0.5 * (
            logvar_p - logvar_q
            + (var_q + (mu_q - mu_p).pow(2)) / var_p
            - 1.0
            )
        return kl.sum(dim=-1)


if __name__ == "__main__":
    torch.manual_seed(42)

    B, n_celltypes, d_z = 8, 24, 32
    prior = CellTypePrior(n_celltypes=n_celltypes, d_z=d_z)
    celltype_id = torch.randint(0, n_celltypes, (B,), dtype=torch.long)

    mu_p, logvar_p, sigma_p = prior(celltype_id)
    print("mu_p.shape     =", mu_p.shape)
    print("logvar_p.shape =", logvar_p.shape)
    print("sigma_p.shape  =", sigma_p.shape)
    print("sigma_p.min()  =", sigma_p.min().item())

    z_s = torch.randn(B, d_z)
    z_perp = prior.residualize(z_s, celltype_id)
    print("z_perp.shape     =", z_perp.shape)

    mu_q = torch.randn(B, d_z)
    logvar_q = torch.randn(B, d_z) * 0.2
    kl = CellTypePrior.kl_divergence(mu_q, logvar_q, mu_p, logvar_p)
    print("kl.shape       =", kl.shape)
    print("kl.min()       =", kl.min().item())

    assert mu_p.shape == (B, d_z)
    assert logvar_p.shape == (B, d_z)
    assert sigma_p.shape == (B, d_z)
    assert z_perp.shape == (B, d_z)
    assert kl.shape == (B,)
    assert torch.all(sigma_p > 0)

    print ("[OK] CellTypePrior smoke test passed.")


