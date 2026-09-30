from __future__ import annotations

from typing import Optional

import torch
import torch.nn as nn


class GeneExpressionEmbedding(nn.Module):
    """
    Gene/token embedding for ordinal single-cell expression data.

    Each token embedding is formed as

        h(g, k, x) = E_gene(g) + E_bin(k) + gate(g) * Proj(x, E_gene(g))

    where
        - E_gene(g): learnable gene identity embedding
        - E_bin(k): learnable ordinal bin embedding
        - Proj(x, E_gene(g)): gene-conditioned projection of continuous value
        - gate(g): learnable gene-wise scalar gate in (0, 1)
    """

    def __init__(
        self,
        n_genes: int,
        d_model: int,
        n_bins: int,
        *,
        use_continuous: bool = True,
        use_cls: bool = True,
        dropout: float = 0.1,
        cont_hidden_dim: Optional[int] = None,
        layer_norm_eps: float = 1e-5,
        init_std: float = 0.02,
    ) -> None:
        super().__init__()

        if n_genes <= 0:
            raise ValueError(f"n_genes must be positive, got {n_genes}.")
        if d_model <= 0:
            raise ValueError(f"d_model must be positive, got {d_model}.")
        if n_bins < 2:
            raise ValueError(f"n_bins must be >= 2, got {n_bins}.")
        if not (0.0 <= dropout < 1.0):
            raise ValueError(f"dropout must be in [0, 1), got {dropout}.")

        self.n_genes = n_genes
        self.d_model = d_model
        self.n_bins = n_bins
        self.use_continuous = use_continuous
        self.use_cls = use_cls
        self.init_std = init_std

        if cont_hidden_dim is None:
            cont_hidden_dim = max(16, d_model // 4)

        # 1) Gene identity embedding
        self.gene_embedding = nn.Embedding(n_genes, d_model)

        # 2) Ordinal bin embedding
        self.bin_embedding = nn.Embedding(n_bins, d_model)

        # 3) Gene-conditioned continuous branch
        if self.use_continuous:
            # input = [gene_embedding, x_log1p] -> output = continuous correction
            self.continuous_proj = nn.Sequential(
                nn.Linear(d_model + 1, cont_hidden_dim),
                nn.GELU(),
                nn.Linear(cont_hidden_dim, d_model),
            )

            # gene-wise gate: one scalar per gene
            self.cont_gate_logit = nn.Parameter(torch.full((n_genes,), -2.0))
        else:
            self.continuous_proj = None
            self.cont_gate_logit = None

        # 4) Optional CLS token
        if self.use_cls:
            self.cls_token = nn.Parameter(torch.zeros(1, 1, d_model))
        else:
            self.cls_token = None

        # 5) Normalization / regularization
        self.layer_norm = nn.LayerNorm(d_model, eps=layer_norm_eps)
        self.dropout = nn.Dropout(dropout)


        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.gene_embedding.weight, mean=0.0, std=self.init_std)
        nn.init.normal_(self.bin_embedding.weight, mean=0.0, std=self.init_std)

        if self.use_continuous:
            for module in self.continuous_proj:
                if isinstance(module, nn.Linear):
                    nn.init.xavier_uniform_(module.weight)
                    nn.init.zeros_(module.bias)

            # start continuous branch gently
            last_linear = self.continuous_proj[-1]
            if isinstance(last_linear, nn.Linear):
                nn.init.zeros_(last_linear.weight)
                nn.init.zeros_(last_linear.bias)

        if self.use_cls:
            nn.init.normal_(self.cls_token, mean=0.0, std=self.init_std)

    def forward(
        self,
        y_ord: torch.Tensor,
        x_log1p: Optional[torch.Tensor] = None,
        gene_ids: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Parameters
        ----------
        y_ord : [B, G] long
            Ordinal bin indices in {0, ..., K-1}

        x_log1p : [B, G] float, optional
            Continuous log1p expression values.
            Required when use_continuous=True.

        gene_ids : None or [G] or [B, G]
            Optional gene identity indices.

        Returns
        -------
        h : [B, G, d_model] or [B, G+1, d_model]
        """
        if y_ord.ndim != 2:
            raise ValueError(f"y_ord must have shape [B, G], got ndim={y_ord.ndim}.")
        if y_ord.dtype != torch.long:
            raise TypeError(f"y_ord must be torch.long, got {y_ord.dtype}.")

        B, G = y_ord.shape

        if G != self.n_genes:
            raise ValueError(
                f"Gene dimension mismatch: expected {self.n_genes}, got {G}."
            )

        if torch.any(y_ord < 0) or torch.any(y_ord >= self.n_bins):
            raise ValueError(
                f"y_ord contains values outside [0, {self.n_bins - 1}]."
            )

        if self.use_continuous:
            if x_log1p is None:
                raise ValueError("x_log1p must be provided when use_continuous=True.")
            if x_log1p.ndim != 2 or x_log1p.shape != (B, G):
                raise ValueError(
                    f"x_log1p must have shape {(B, G)}, got {tuple(x_log1p.shape)}."
                )
            if not torch.is_floating_point(x_log1p):
                raise TypeError(
                    f"x_log1p must be a floating tensor, got {x_log1p.dtype}."
                )

        gene_ids_resolved = self._resolve_gene_ids(gene_ids, B, G, y_ord.device)

        # gene embedding
        if gene_ids_resolved.ndim == 1:
            h_gene = self.gene_embedding(gene_ids_resolved).unsqueeze(0)  # [1, G, d]
        else:
            h_gene = self.gene_embedding(gene_ids_resolved)  # [B, G, d]

        # ordinal bin embedding
        h_bin = self.bin_embedding(y_ord)  # [B, G, d]

        # base token embedding
        h = h_gene + h_bin  # broadcast if h_gene is [1, G, d]

        # gene-conditioned continuous correction
        if self.use_continuous:
            if h_gene.shape[0] == 1:
                h_gene_for_cont = h_gene.expand(B, -1, -1)  # [B, G, d]
            else:
                h_gene_for_cont = h_gene  # [B, G, d]

            x_feat = x_log1p.unsqueeze(-1)  # [B, G, 1]
            cont_input = torch.cat([h_gene_for_cont, x_feat], dim=-1)  # [B, G, d+1]
            h_cont = self.continuous_proj(cont_input)  # [B, G, d]

            cont_gate = torch.sigmoid(self.cont_gate_logit).view(1, G, 1)  # [1, G, 1]
            h = h + cont_gate * h_cont

        h = self.layer_norm(h)
        h = self.dropout(h)

        if self.use_cls:
            cls = self.cls_token.expand(B, -1, -1)  # [B, 1, d]
            h = torch.cat([cls, h], dim=1)  # [B, G+1, d]

        return h

    def _resolve_gene_ids(
        self,
        gene_ids: Optional[torch.Tensor],
        batch_size: int,
        n_genes: int,
        device: torch.device,
    ) -> torch.Tensor:
        if gene_ids is None:
            return torch.arange(self.n_genes, device=device, dtype=torch.long)

        if gene_ids.dtype != torch.long:
            raise TypeError(f"gene_ids must be torch.long, got {gene_ids.dtype}.")

        if gene_ids.ndim == 1:
            if gene_ids.shape[0] != n_genes:
                raise ValueError(
                    f"1D gene_ids must have shape [{n_genes}], got {tuple(gene_ids.shape)}."
                )
            return gene_ids.to(device)

        if gene_ids.ndim == 2:
            if gene_ids.shape != (batch_size, n_genes):
                raise ValueError(
                    f"2D gene_ids must have shape [{batch_size}, {n_genes}], "
                    f"got {tuple(gene_ids.shape)}."
                )
            return gene_ids.to(device)

        raise ValueError(
            f"gene_ids must be None, 1D [G], or 2D [B, G]; got ndim={gene_ids.ndim}."
        )

    @property
    def output_dim(self) -> int:
        return self.d_model

    @property
    def seq_length(self) -> int:
        return self.n_genes + (1 if self.use_cls else 0)

    @property
    def continuous_gate(self) -> Optional[torch.Tensor]:
        if not self.use_continuous:
            return None
        return torch.sigmoid(self.cont_gate_logit).detach().cpu()

    def extra_repr(self) -> str:
        if self.use_continuous:
            gate = self.continuous_gate
            gate_str = (
                f"mean={gate.mean().item():.4f}, "
                f"min={gate.min().item():.4f}, "
                f"max={gate.max().item():.4f}"
            )
        else:
            gate_str = "None"

        return (
            f"n_genes={self.n_genes}, d_model={self.d_model}, "
            f"n_bins={self.n_bins}, use_continuous={self.use_continuous}, "
            f"use_cls={self.use_cls}, continuous_gate={gate_str}"
        )


if __name__ == "__main__":
    torch.manual_seed(42)

    B, G, K, D = 4, 2000, 10, 128

    embed = GeneExpressionEmbedding(
        n_genes=G,
        d_model=D,
        n_bins=K,
        use_continuous=True,
        use_cls=True,
        dropout=0.1,
    )
    print(embed)

    y_ord = torch.randint(0, K, (B, G), dtype=torch.long)
    x_log1p = torch.rand(B, G, dtype=torch.float32) * 6.0

    h = embed(y_ord=y_ord, x_log1p=x_log1p)
    print("h.shape =", h.shape)

    loss = h.sum()
    loss.backward()

    print("gene grad:", embed.gene_embedding.weight.grad is not None)
    print("bin grad:", embed.bin_embedding.weight.grad is not None)
    print("gate grad:", embed.cont_gate_logit.grad is not None)
    print("done")