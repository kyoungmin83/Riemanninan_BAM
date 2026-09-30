"""
Per-celltype ANCOVA ridge regression on noise_score.

For each celltype t, we fit:
    bank.mean[d, t] ≈ β_0(t) + Σ_a β_a(t) * pathology_a(d)

The residual (donor-individual deviation after pathology control) drives the
L2 uncertainty component:
    u_donor(c) = residual_mean[d, t] / residual_std[t]

This is *not* a pathology signal extractor. The β coefficients only describe
how a donor's pathology profile shifts their average noise level; the
direction-aware pathology signal is handled in v7b's z_state heads.

See: model_v7a_implementation_plan.md §4.2.5, model_v7_design...md §5.1, §5.2
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import torch

from kmlee_bam.uncertainty.bank import NoiseScoreBank


@dataclass
class ANCOVAConfig:
    ridge_alpha: float = 0.1
    min_donors_per_celltype: int = 5
    pathology_axes: List[str] = field(
        default_factory=lambda: ["Braak_stage", "Thal_phase", "CERAD_score"]
    )
    standardize_axes: bool = True  # z-score donor-pathology across donors


def _ridge_regress(
    X: torch.Tensor,   # [N, P]
    y: torch.Tensor,   # [N]
    *,
    alpha: float,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Closed-form ridge regression. Adds an intercept column internally.

    Returns: (beta [P], intercept scalar). The intercept is NOT penalized.
    """
    X = X.float()
    y = y.float()
    n, p = X.shape

    # Mean-center X and y so the intercept is just y_mean - X_mean·beta.
    x_mean = X.mean(dim=0)
    y_mean = y.mean()
    Xc = X - x_mean.unsqueeze(0)
    yc = y - y_mean

    # Solve (X^T X + αI) β = X^T y
    XtX = Xc.t() @ Xc                                    # [P, P]
    reg = float(alpha) * torch.eye(p, dtype=X.dtype, device=X.device)
    Xty = Xc.t() @ yc                                    # [P]
    try:
        beta = torch.linalg.solve(XtX + reg, Xty)
    except RuntimeError:
        # Numerical fallback: lstsq on the augmented system.
        beta = torch.linalg.lstsq(XtX + reg, Xty.unsqueeze(-1)).solution.squeeze(-1)

    intercept = y_mean - x_mean @ beta
    return beta, intercept


class ANCOVAFit:
    """
    Per-celltype ANCOVA. Refit at epoch boundaries.

    Public tensors:
        beta:           [T, P]   regression coefficients
        intercept:      [T]      per-celltype intercept
        residual_mean:  [D, T]   donor-individual residual (bank.mean - predicted)
        residual_std:   [T]      std of residual_mean across donors
        is_fitted:      [T] bool
    """

    def __init__(
        self,
        config: ANCOVAConfig,
        n_donors: int,
        n_celltypes: int,
        n_axes: int,
        device: Optional[torch.device] = None,
    ):
        self.config = config
        self.n_donors = int(n_donors)
        self.n_celltypes = int(n_celltypes)
        self.n_axes = int(n_axes)
        self.device = device if device is not None else torch.device("cpu")

        self.beta = torch.zeros(n_celltypes, n_axes, dtype=torch.float32, device=self.device)
        self.intercept = torch.zeros(n_celltypes, dtype=torch.float32, device=self.device)
        self.residual_mean = torch.zeros(
            n_donors, n_celltypes, dtype=torch.float32, device=self.device
        )
        self.residual_std = torch.ones(
            n_celltypes, dtype=torch.float32, device=self.device
        )
        self.is_fitted = torch.zeros(n_celltypes, dtype=torch.bool, device=self.device)

    @torch.no_grad()
    def fit(
        self,
        bank: NoiseScoreBank,
        donor_pathology: torch.Tensor,        # [D, P]
        donor_celltype_mask: torch.Tensor,    # [D, T] bool — donor has cells of this celltype
        *,
        donor_pathology_valid: Optional[torch.Tensor] = None,  # [D] bool
    ) -> dict:
        """
        Refit β, intercept, residual_mean, residual_std.

        donor_pathology_valid: cells we trust to use in regression. Donors
            missing pathology labels are excluded from fitting but their
            residual is set to 0 (no donor-level correction).
        """
        device = self.device
        donor_pathology = donor_pathology.to(device).float()       # [D, P]
        donor_celltype_mask = donor_celltype_mask.to(device).bool()  # [D, T]
        if donor_pathology_valid is None:
            donor_pathology_valid = torch.ones(self.n_donors, dtype=torch.bool, device=device)
        else:
            donor_pathology_valid = donor_pathology_valid.to(device).bool()

        bank_mean = bank.get_mean_tensor().to(device)              # [D, T]
        bank_init = bank.get_initialized_mask().to(device)         # [D, T]

        # Optional standardization of pathology axes across the *valid* donors
        # (so β coefficients are comparable across axes).
        if self.config.standardize_axes and bool(donor_pathology_valid.any()):
            X_all = donor_pathology[donor_pathology_valid]         # [Nv, P]
            x_mean_axis = X_all.mean(dim=0, keepdim=True)
            x_std_axis = X_all.std(dim=0, unbiased=False, keepdim=True).clamp_min(1e-6)
            donor_path_z = (donor_pathology - x_mean_axis) / x_std_axis
        else:
            donor_path_z = donor_pathology

        n_fitted = 0
        n_skipped_lowcount = 0
        for t in range(self.n_celltypes):
            valid_mask = (
                donor_celltype_mask[:, t]
                & bank_init[:, t]
                & donor_pathology_valid
            )
            valid_idx = valid_mask.nonzero(as_tuple=False).squeeze(-1)
            if int(valid_idx.numel()) < self.config.min_donors_per_celltype:
                self.is_fitted[t] = False
                self.beta[t].zero_()
                self.intercept[t] = 0.0
                self.residual_mean[:, t] = 0.0
                self.residual_std[t] = 1.0
                n_skipped_lowcount += 1
                continue

            X = donor_path_z[valid_idx]                            # [Nv, P]
            y = bank_mean[valid_idx, t]                            # [Nv]
            beta, intercept = _ridge_regress(X, y, alpha=self.config.ridge_alpha)
            self.beta[t] = beta
            self.intercept[t] = intercept

            predicted_t = intercept + donor_path_z @ beta          # [D]
            # Only donors that were used in the fit get a residual; others
            # get residual = 0 (we cannot judge them).
            self.residual_mean[:, t] = 0.0
            self.residual_mean[valid_idx, t] = (
                bank_mean[valid_idx, t] - predicted_t[valid_idx]
            )

            # Residual std across the fitted donors.
            self.residual_std[t] = self.residual_mean[valid_idx, t].std(
                unbiased=False
            ).clamp_min(1e-6)
            self.is_fitted[t] = True
            n_fitted += 1

        return {
            "n_celltypes_fitted": int(n_fitted),
            "n_celltypes_skipped_lowcount": int(n_skipped_lowcount),
            "n_celltypes_total": int(self.n_celltypes),
            "beta_abs_mean": float(self.beta[self.is_fitted].abs().mean())
            if int(self.is_fitted.sum()) > 0
            else 0.0,
            "residual_std_mean": float(
                self.residual_std[self.is_fitted].mean()
            )
            if int(self.is_fitted.sum()) > 0
            else 0.0,
        }

    # ------------------------------------------------------------------ #
    # State dict
    # ------------------------------------------------------------------ #

    def state_dict(self) -> dict:
        return {
            "beta": self.beta.detach().cpu(),
            "intercept": self.intercept.detach().cpu(),
            "residual_mean": self.residual_mean.detach().cpu(),
            "residual_std": self.residual_std.detach().cpu(),
            "is_fitted": self.is_fitted.detach().cpu(),
        }

    def load_state_dict(self, state: dict) -> None:
        device = self.device
        self.beta = state["beta"].to(device)
        self.intercept = state["intercept"].to(device)
        self.residual_mean = state["residual_mean"].to(device)
        self.residual_std = state["residual_std"].to(device)
        self.is_fitted = state["is_fitted"].to(device)
