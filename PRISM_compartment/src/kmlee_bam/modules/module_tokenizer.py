from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


PoolingMode = Literal["mean", "attention"]
ActivityNorm = Literal["l2", "row_sum", "none"]
ActivityDefault = Literal["l2_membership", "mean_membership"]


@dataclass
class GeneModuleTokenizerOutput:
    """
    Output of GeneModuleTokenizer.

    Attributes
    ----------
    tokens:
        [B, M, d] or [B, M+1, d] if preserve_cls=True.
        This tensor is intended for the downstream Transformer/StateEncoder.

    module_tokens:
        [B, M, d], same as tokens but without CLS.

    pooled_gene_state:
        [B, M, d], expression-derived hidden module state before adding
        module identity and scalar activity projection.

    module_activity:
        [B, M] or None. Interpretable scalar activity per cell and module.
        This is the biology-facing module activation score.

    module_attention:
        [B, M, G] or None. Only returned for attention pooling when requested.
        It is detached from the graph because it is intended as a diagnostic.
    """

    tokens: torch.Tensor
    module_tokens: torch.Tensor
    pooled_gene_state: torch.Tensor
    module_activity: Optional[torch.Tensor]
    module_attention: Optional[torch.Tensor]


class GeneModuleTokenizer(nn.Module):
    """
    Convert gene-level tokens into biological module-level tokens.

    Design
    ------
    This module separates three objects that should not be conflated:

        1. Module identity
           E_module[m], learned by nn.Embedding(n_modules, d_model).

        2. Module hidden state token
           A pooled gene-token representation inside each module:

               S[i, m, :] = sum_g membership_weight[m, g] * H_gene[i, g, :]

        3. Module activity scalar
           An interpretable scalar activation score:

               a[i, m] = sum_g activity_weight[m, g] * x_gene_scalar[i, g]

    The output module token is roughly:

        T[i, m, :] = LN(
            pooled_gene_state[i, m, :]
            + module_id_scale * E_module[m]
            + activity_projection(a[i, m])
        )

    Crucial distinction
    -------------------
    membership_weight and activity_weight have different meanings.

    - membership_weight is non-negative and row-sum-normalized. It pools
      gene hidden states into module hidden states.

    - activity_weight can be signed and is preferably L2-normalized. It defines
      the direction along which a cell's gene-expression vector is projected to
      obtain module activity.

    If no activity_weight is supplied, this implementation does NOT use the
    row-normalized mean by default. Instead, it uses L2-normalized positive
    membership weights:

        activity_weight[m, g] = 1 / sqrt(|G_m|), if g in module m.

    This makes activity scores more comparable across module sizes than the
    simple 1/|G_m| row mean. Set activity_default="mean_membership" only when
    you explicitly want the simple row mean.

    x_gene_scalar used in forward
    -----------------------------
    Prefer train-set standardized log1p expression:

        x_gene_scalar[i, g] = (log1p_count[i, g] - mean_control_train[celltype_i, g]) / (std_control_train[cell_type_i, g] + eps)

    For quick smoke tests, normalized ordinal bins are acceptable, but they are
    less suitable for final biological interpretation.

    No source/type embedding is used here. In particular, there is no
    regular/residual/singleton embedding, because that is administrative
    metadata rather than biological module identity.
    """

    def __init__(
        self,
        *,
        membership_weight: torch.Tensor,                  # [M, G]
        d_model: int,
        activity_weight: Optional[torch.Tensor] = None,   # [M, G], signed allowed
        activity_weight_normalization: ActivityNorm = "l2",
        activity_default: ActivityDefault = "l2_membership",
        pooling: PoolingMode = "mean",
        preserve_cls: bool = True,
        use_module_id_embedding: bool = True,
        module_id_init_scale: float = 0.01,
        module_id_scale_min: float = 1e-4,
        module_id_scale_max: float = 10.0,
        use_activity_projection: bool = True,
        dropout: float = 0.1,
        layer_norm_eps: float = 1e-5,
        init_std: float = 0.02,
    ) -> None:
        super().__init__()

        if membership_weight.ndim != 2:
            raise ValueError("membership_weight must have shape [M, G].")
        if d_model <= 0:
            raise ValueError("d_model must be positive.")
        if pooling not in {"mean", "attention"}:
            raise ValueError("pooling must be either 'mean' or 'attention'.")
        if activity_weight_normalization not in {"l2", "row_sum", "none"}:
            raise ValueError("activity_weight_normalization must be 'l2', 'row_sum', or 'none'.")
        if activity_default not in {"l2_membership", "mean_membership"}:
            raise ValueError("activity_default must be 'l2_membership' or 'mean_membership'.")
        if module_id_init_scale <= 0:
            raise ValueError("module_id_init_scale must be positive.")
        if module_id_scale_min <= 0:
            raise ValueError("module_id_scale_min must be positive.")
        if module_id_scale_max <= module_id_scale_min:
            raise ValueError("module_id_scale_max must be greater than module_id_scale_min.")
        if not (module_id_scale_min <= module_id_init_scale <= module_id_scale_max):
            raise ValueError(
                "module_id_init_scale must lie within "
                f"[{module_id_scale_min}, {module_id_scale_max}]."
            )
        if not (0.0 <= dropout < 1.0):
            raise ValueError("dropout must be in [0, 1).")
        if not torch.isfinite(membership_weight).all():
            raise ValueError("membership_weight contains non-finite values.")
        if (membership_weight < 0).any():
            raise ValueError("membership_weight must be non-negative.")

        raw_membership = membership_weight.to(dtype=torch.float32)
        membership_mask = raw_membership > 0
        row_sum = raw_membership.sum(dim=1, keepdim=True)
        if (row_sum <= 0).any():
            bad = torch.where(row_sum.squeeze(1) <= 0)[0].detach().cpu().tolist()[:10]
            raise ValueError(f"membership_weight has empty module rows. Examples: {bad}")

        # Used only for pooling hidden gene tokens. Row-sum normalization is correct here.
        normalized_membership = raw_membership / row_sum.clamp_min(1e-12)

        n_modules, n_genes = normalized_membership.shape
        self.n_modules = int(n_modules)
        self.n_genes = int(n_genes)
        self.d_model = int(d_model)
        self.pooling = pooling
        self.preserve_cls = bool(preserve_cls)
        self.use_module_id_embedding = bool(use_module_id_embedding)
        self.use_activity_projection = bool(use_activity_projection)
        self.module_id_scale_min = float(module_id_scale_min)
        self.module_id_scale_max = float(module_id_scale_max)
        self.activity_weight_normalization = activity_weight_normalization
        self.activity_default = activity_default

        self.register_buffer("membership_weight", normalized_membership)
        self.register_buffer("membership_mask", membership_mask)
        self.register_buffer(
            "module_ids",
            torch.arange(self.n_modules, dtype=torch.long),
            persistent=False,
        )

        # ------------------------------------------------------------------
        # Activity weights: scalar biological module activity.
        # ------------------------------------------------------------------
        if activity_weight is None:
            if activity_default == "l2_membership":
                # Positive L2 gene-set score. For a module with k genes,
                # each member gets 1/sqrt(k), giving ||v_m||_2 = 1.
                act_w = membership_mask.to(dtype=torch.float32)
                denom = act_w.pow(2).sum(dim=1, keepdim=True).sqrt()
                act_w = act_w / denom.clamp_min(1e-12)
            else:
                # Explicit simple mean fallback: 1/k per member.
                act_w = normalized_membership.clone()
        else:
            if activity_weight.shape != normalized_membership.shape:
                raise ValueError(
                    "activity_weight must have the same shape as membership_weight: "
                    f"expected {tuple(normalized_membership.shape)}, got {tuple(activity_weight.shape)}."
                )
            if not torch.isfinite(activity_weight).all():
                raise ValueError("activity_weight contains non-finite values.")

            act_w = activity_weight.to(dtype=torch.float32)

            # Enforce the declared module definition: no activity contribution
            # from genes outside the module membership.
            act_w = act_w * membership_mask.to(dtype=act_w.dtype)

            abs_row_sum = act_w.abs().sum(dim=1, keepdim=True)
            if (abs_row_sum <= 0).any():
                bad = torch.where(abs_row_sum.squeeze(1) <= 0)[0].detach().cpu().tolist()[:10]
                raise ValueError(
                    "activity_weight has empty/zero rows after membership masking. "
                    f"Examples: {bad}"
                )

            if activity_weight_normalization == "l2":
                denom = act_w.pow(2).sum(dim=1, keepdim=True).sqrt()
                act_w = act_w / denom.clamp_min(1e-12)
            elif activity_weight_normalization == "row_sum":
                if (act_w < 0).any():
                    raise ValueError("row_sum normalization requires non-negative activity_weight.")
                denom = act_w.sum(dim=1, keepdim=True)
                if (denom <= 0).any():
                    bad = torch.where(denom.squeeze(1) <= 0)[0].detach().cpu().tolist()[:10]
                    raise ValueError(f"activity_weight has zero row sums after masking. Examples: {bad}")
                act_w = act_w / denom.clamp_min(1e-12)
            # "none" leaves weights as supplied after masking.

        self.register_buffer("activity_weight", act_w)

        # ------------------------------------------------------------------
        # Biological module identity embedding: pathway1/pathway2/hdWGCNA_m...
        # ------------------------------------------------------------------
        if self.use_module_id_embedding:
            self.module_id_embedding = nn.Embedding(self.n_modules, d_model)
            nn.init.normal_(self.module_id_embedding.weight, mean=0.0, std=init_std)

            # Positive scale via exp(log_scale), with a clamp in the property.
            self.module_id_log_scale = nn.Parameter(
                torch.tensor(math.log(float(module_id_init_scale)), dtype=torch.float32)
            )
        else:
            self.module_id_embedding = None
            self.module_id_log_scale = None

        # ------------------------------------------------------------------
        # Optional activity scalar -> d_model projection.
        # This injects interpretable scalar module activity into the token.
        # ------------------------------------------------------------------
        if self.use_activity_projection:
            self.activity_projection = nn.Linear(1, d_model)
            nn.init.normal_(self.activity_projection.weight, mean=0.0, std=init_std)
            nn.init.zeros_(self.activity_projection.bias)
        else:
            self.activity_projection = None

        # ------------------------------------------------------------------
        # Optional attention pooling. Default is mean pooling.
        # ------------------------------------------------------------------
        if self.pooling == "attention":
            self.module_queries = nn.Parameter(torch.empty(self.n_modules, d_model))
            nn.init.normal_(self.module_queries, mean=0.0, std=init_std)
        else:
            self.module_queries = None

        self.layer_norm = nn.LayerNorm(d_model, eps=layer_norm_eps)
        self.dropout = nn.Dropout(dropout)

    @property
    def module_id_scale(self) -> Optional[torch.Tensor]:
        if self.module_id_log_scale is None:
            return None
        log_min = math.log(self.module_id_scale_min)
        log_max = math.log(self.module_id_scale_max)
        return torch.exp(self.module_id_log_scale.clamp(min=log_min, max=log_max))

    def extra_repr(self) -> str:
        scale_str = "None"
        if self.module_id_log_scale is not None:
            scale_str = f"{float(self.module_id_scale.detach().cpu()):.4g}"
        return (
            f"n_modules={self.n_modules}, n_genes={self.n_genes}, "
            f"d_model={self.d_model}, pooling='{self.pooling}', "
            f"preserve_cls={self.preserve_cls}, "
            f"use_module_id_embedding={self.use_module_id_embedding}, "
            f"module_id_scale={scale_str}, "
            f"use_activity_projection={self.use_activity_projection}, "
            f"activity_default='{self.activity_default}', "
            f"activity_weight_normalization='{self.activity_weight_normalization}'"
        )

    def _split_cls(self, gene_tokens: torch.Tensor) -> tuple[Optional[torch.Tensor], torch.Tensor]:
        if self.preserve_cls:
            cls_token = gene_tokens[:, :1, :]   # [B, 1, d]
            h_gene = gene_tokens[:, 1:, :]      # [B, G, d]
        else:
            cls_token = None
            h_gene = gene_tokens                # [B, G, d]
        return cls_token, h_gene

    def _pool_mean(self, h_gene: torch.Tensor) -> torch.Tensor:
        # membership_weight: [M, G], h_gene: [B, G, d] -> [B, M, d]
        weight = self.membership_weight.to(device=h_gene.device, dtype=h_gene.dtype)
        return torch.einsum("mg,bgd->bmd", weight, h_gene)

    def _pool_attention(
        self,
        h_gene: torch.Tensor,
        *,
        return_attention: bool,
    ) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
        if self.module_queries is None:
            raise RuntimeError("module_queries is None; pooling is not attention.")

        queries = self.module_queries.to(device=h_gene.device, dtype=h_gene.dtype)

        # logits: [B, M, G]
        logits = torch.einsum("md,bgd->bmg", queries, h_gene)
        logits = logits / math.sqrt(float(self.d_model))

        mask = self.membership_mask.to(device=h_gene.device).unsqueeze(0)  # [1, M, G]
        logits = logits.masked_fill(~mask, torch.finfo(logits.dtype).min)
        attn = F.softmax(logits, dim=-1)

        # Numerical safety: ensure no leakage outside module membership.
        attn = attn * mask.to(dtype=attn.dtype)
        attn = attn / attn.sum(dim=-1, keepdim=True).clamp_min(1e-12)

        pooled = torch.einsum("bmg,bgd->bmd", attn, h_gene)
        return pooled, attn.detach() if return_attention else None

    def forward(
        self,
        gene_tokens: torch.Tensor,                       # [B, G(+1), d]
        x_gene_scalar: Optional[torch.Tensor] = None,    # [B, G]
        *,
        return_attention: bool = False,
    ) -> GeneModuleTokenizerOutput:
        if gene_tokens.ndim != 3:
            raise ValueError(
                f"gene_tokens must have shape [B, T, d], got {tuple(gene_tokens.shape)}."
            )

        B, T, D = gene_tokens.shape
        if D != self.d_model:
            raise ValueError(f"d_model mismatch: expected {self.d_model}, got {D}.")

        expected_t = self.n_genes + 1 if self.preserve_cls else self.n_genes
        if T != expected_t:
            raise ValueError(
                f"Token length mismatch: expected {expected_t}, got {T}. "
                "This usually means dataset gene order/size does not match the registry."
            )

        cls_token, h_gene = self._split_cls(gene_tokens)

        if h_gene.shape[1] != self.n_genes:
            raise RuntimeError(
                f"Internal gene-token mismatch: expected {self.n_genes}, got {h_gene.shape[1]}."
            )

        # 1) Module hidden state from gene tokens.
        if self.pooling == "mean":
            pooled_gene_state = self._pool_mean(h_gene)
            module_attention = None
        else:
            pooled_gene_state, module_attention = self._pool_attention(
                h_gene,
                return_attention=return_attention,
            )

        module_state = pooled_gene_state

        # 2) Explicit scalar module activity.
        module_activity = None
        if x_gene_scalar is not None:
            if x_gene_scalar.ndim != 2 or x_gene_scalar.shape != (B, self.n_genes):
                raise ValueError(
                    f"x_gene_scalar must have shape [B, G] = {(B, self.n_genes)}, "
                    f"got {tuple(x_gene_scalar.shape)}."
                )

            x_gene_scalar = x_gene_scalar.to(device=gene_tokens.device, dtype=self.activity_weight.dtype)
            act_w = self.activity_weight.to(device=gene_tokens.device, dtype=x_gene_scalar.dtype)
            module_activity = torch.einsum("mg,bg->bm", act_w, x_gene_scalar)

            if self.activity_projection is not None:
                activity_delta = self.activity_projection(module_activity.unsqueeze(-1))
                module_state = module_state + activity_delta.to(dtype=module_state.dtype)

        # 3) Biological module identity.
        if self.module_id_embedding is not None:
            h_id = self.module_id_embedding(self.module_ids)  # [M, d]
            h_id = h_id.to(dtype=module_state.dtype)
            scale = self.module_id_scale.to(dtype=module_state.dtype)
            module_state = module_state + scale * h_id.unsqueeze(0)

        # 4) Final normalization/dropout for downstream Transformer/StateEncoder.
        module_tokens = self.layer_norm(module_state)
        module_tokens = self.dropout(module_tokens)

        if cls_token is not None:
            tokens = torch.cat([cls_token, module_tokens], dim=1)
        else:
            tokens = module_tokens

        return GeneModuleTokenizerOutput(
            tokens=tokens,
            module_tokens=module_tokens,
            pooled_gene_state=pooled_gene_state,
            module_activity=module_activity,
            module_attention=module_attention,
        )


def _sanity_check() -> None:
    torch.manual_seed(42)

    B, G, M, D = 4, 1000, 50, 64

    membership = torch.zeros(M, G)
    for m in range(M):
        n_g = torch.randint(5, 30, (1,)).item()
        idx = torch.randperm(G)[:n_g]
        membership[m, idx] = 1.0

    # Signed eigengene-like activity weights. Zero outside membership.
    activity = torch.randn(M, G) * membership

    tokenizer = GeneModuleTokenizer(
        membership_weight=membership,
        d_model=D,
        activity_weight=activity,
        activity_weight_normalization="l2",
        pooling="mean",
        preserve_cls=True,
        dropout=0.0,
    )
    print(tokenizer)

    gene_tokens = torch.randn(B, G + 1, D, requires_grad=True)
    x_scalar = torch.randn(B, G)

    out = tokenizer(gene_tokens, x_gene_scalar=x_scalar)
    print(f"tokens.shape           = {tuple(out.tokens.shape)}")
    print(f"module_tokens.shape    = {tuple(out.module_tokens.shape)}")
    print(f"pooled_gene_state.shape= {tuple(out.pooled_gene_state.shape)}")
    print(f"module_activity.shape  = {tuple(out.module_activity.shape)}")

    assert out.tokens.shape == (B, M + 1, D)
    assert out.module_tokens.shape == (B, M, D)
    assert out.pooled_gene_state.shape == (B, M, D)
    assert out.module_activity.shape == (B, M)

    loss = out.tokens.pow(2).mean() + 0.01 * out.module_activity.pow(2).mean()
    loss.backward()

    assert gene_tokens.grad is not None
    assert tokenizer.module_id_embedding.weight.grad is not None
    assert tokenizer.module_id_log_scale.grad is not None
    assert tokenizer.activity_projection.weight.grad is not None
    print("OK")


if __name__ == "__main__":
    _sanity_check()