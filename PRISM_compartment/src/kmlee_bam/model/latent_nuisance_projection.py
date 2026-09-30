"""v25/v28 — STRUCTURAL nuisance projection-out for z_perp (MANIFEST).

The leak↔spread tension: any spreading pressure (whitening / path_aux) refills the freed z-axes
with the STRONGEST easy signal — in SEA-AD that is SEX (then cell type). A soft penalty only
*pushes* and loses the tug-of-war once z spreads. This module *removes the door*: it estimates the
nuisance subspace in z and HARD-projects it out, so whitening operates on the orthogonal complement.

v28 — MULTI-RANK sex subspace (was 1-rank). The deployed 1-rank u_sex (mean of within-celltype
male-female diffs) left z→sex linearly decodable at ~0.77 and it REFILLED over epochs. Fix: collect
the per-celltype male-female mean-diffs into an EMA SCATTER  Σ_c d_c d_cᵀ  and project out its top
`sex_rank` eigenvectors (the SEX SUBSPACE). The post-hoc sweep: rank-8 → z→sex ~0.52 (≈chance) at no
disease cost. Scatter is sign-invariant (no EMA sign-flip bug). Default-OFF ⇒ identity.

Diagnostics (per forward, for the training log / refill watch): sex eigenvalue spectrum, the rank
ACTUALLY removed, and the fraction of z's energy the projection removes.
"""
from __future__ import annotations
from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.distributed as dist


@dataclass
class NuisanceProjectionConfig:
    enabled: bool = False
    project_sex: bool = True
    project_celltype: bool = False     # remove the next-strongest nuisance too
    sex_rank: int = 8                  # v28: multi-rank sex subspace (1 ≈ old single-direction behaviour)
    celltype_rank: int = 4             # cap on celltype centroid directions removed
    ema_decay: float = 0.99
    eps: float = 1e-6
    projection_aware_whitening: bool = False


class LatentNuisanceProjector(nn.Module):
    """Hard-projects a MULTI-RANK nuisance subspace out of z_perp; EMA-estimated, within-celltype."""

    def __init__(self, *, d_z: int, n_celltype: int, config: NuisanceProjectionConfig) -> None:
        super().__init__()
        if int(d_z) <= 0:
            raise ValueError("d_z must be positive")
        if int(n_celltype) <= 0:
            raise ValueError("n_celltype must be positive")
        if int(config.sex_rank) < 0:
            raise ValueError("nuisance_projection.sex_rank must be non-negative")
        if int(config.celltype_rank) < 0:
            raise ValueError(
                "nuisance_projection.celltype_rank must be non-negative"
            )
        if not 0.0 <= float(config.ema_decay) <= 1.0:
            raise ValueError("nuisance_projection.ema_decay must lie in [0, 1]")
        if float(config.eps) <= 0.0:
            raise ValueError("nuisance_projection.eps must be positive")
        self.cfg = config
        self.d_z = int(d_z)
        self.n_celltype = int(n_celltype)
        self.decay = float(config.ema_decay)
        self.eps = float(config.eps)
        self.sex_rank = int(getattr(config, "sex_rank", 8))
        # EMA scatter of within-celltype male-female mean-diffs; its top-k eigvecs = the sex subspace.
        # _update ALL-REDUCES across DDP ranks so every rank's projection is identical + sees all data.
        self.register_buffer("sex_scatter", torch.zeros(self.d_z, self.d_z))
        self.register_buffer("ct_centroid", torch.zeros(self.n_celltype, self.d_z))
        self._diag: dict = {}          # set each forward; read by the trainer for the log

    @torch.no_grad()
    def _update(self, z: torch.Tensor, sex_id, celltype_id) -> None:
        ddp = dist.is_available() and dist.is_initialized()
        cts = torch.unique(celltype_id)
        if self.cfg.project_sex and sex_id is not None:
            # Global sufficient statistics must be reduced BEFORE forming the
            # nonlinear outer product.  Averaging rank-local d_c d_c^T is not
            # equal to the outer product of the global male-female contrast and
            # spuriously inflates rank when DDP shards have different means.
            male_sum = z.new_zeros(self.n_celltype, self.d_z)
            female_sum = z.new_zeros(self.n_celltype, self.d_z)
            male_count = z.new_zeros(self.n_celltype)
            female_count = z.new_zeros(self.n_celltype)
            valid_ct = (celltype_id >= 0) & (celltype_id < self.n_celltype)
            for sex_value, group_sum, group_count in (
                (1, male_sum, male_count),
                (0, female_sum, female_count),
            ):
                keep = valid_ct & (sex_id == sex_value)
                if bool(keep.any()):
                    index = celltype_id[keep].long()
                    group_sum.index_add_(0, index, z[keep])
                    group_count.index_add_(
                        0,
                        index,
                        torch.ones_like(index, dtype=z.dtype),
                    )
            if ddp:
                for statistic in (
                    male_sum,
                    female_sum,
                    male_count,
                    female_count,
                ):
                    dist.all_reduce(statistic, op=dist.ReduceOp.SUM)
            eligible = (male_count >= 2.0) & (female_count >= 2.0)
            cnt = eligible.to(dtype=z.dtype).sum()
            if float(cnt) > self.eps:
                male_mean = male_sum / male_count.clamp_min(1.0)[:, None]
                female_mean = female_sum / female_count.clamp_min(1.0)[:, None]
                contrast = male_mean[eligible] - female_mean[eligible]
                scatter = contrast.t() @ contrast
                self.sex_scatter.mul_(self.decay).add_(
                    scatter / cnt,
                    alpha=1.0 - self.decay,
                )
        if self.cfg.project_celltype:
            csum = z.new_zeros(self.n_celltype, self.d_z); ccnt = z.new_zeros(self.n_celltype)
            for c in cts:
                ci = int(c)
                if 0 <= ci < self.n_celltype:
                    mc = celltype_id == c
                    csum[ci] = z[mc].sum(0); ccnt[ci] = float(mc.sum())
            if ddp:
                dist.all_reduce(csum); dist.all_reduce(ccnt)
            for ci in range(self.n_celltype):
                if float(ccnt[ci]) > self.eps:
                    self.ct_centroid[ci].mul_(self.decay).add_(csum[ci] / ccnt[ci], alpha=1.0 - self.decay)

    def _basis(self, dtype) -> torch.Tensor | None:
        cols = []
        n_sex = 0
        if self.cfg.project_sex and float(self.sex_scatter.norm()) > self.eps:
            S = self.sex_scatter.to(dtype)
            evals, evecs = torch.linalg.eigh(S)             # ascending eigenvalues, orthonormal eigvecs
            self._diag["sex_eigvals_top"] = evals.flip(0)[:8].tolist()   # descending top-8
            k = min(int(self.sex_rank), self.d_z)
            # In Python ``tensor[-0:]`` means the whole tensor, not an empty
            # slice.  Guard rank zero explicitly so ``sex_rank=0`` is a true
            # no-op rather than accidentally deleting every estimated axis.
            if k > 0:
                ev_top = evals[-k:]
                vec_top = evecs[:, -k:]
                thr = (
                    float(ev_top.max()) * 1e-4
                    if float(ev_top.max()) > 0
                    else 0.0
                )
                for j in range(vec_top.shape[1]):
                    if float(ev_top[j]) > thr:              # drop ~degenerate directions
                        cols.append(vec_top[:, j]); n_sex += 1
        n_ct = 0
        if self.cfg.project_celltype and float(self.ct_centroid.norm()) > self.eps:
            Cc = self.ct_centroid.to(dtype)
            Cc = Cc - Cc.mean(0, keepdim=True)
            try:
                _, S_ct, Vh = torch.linalg.svd(Cc, full_matrices=False)
                k = min(int(self.cfg.celltype_rank), Vh.shape[0])
                thr_ct = float(S_ct[0]) * 1e-4 if S_ct.numel() and float(S_ct[0]) > 0 else 0.0
                for i in range(k):
                    if float(S_ct[i]) > thr_ct:             # skip zero-singular-value PADDING directions
                        cols.append(Vh[i]); n_ct += 1
            except Exception:
                pass
        self._diag["n_sex_dirs"] = n_sex
        self._diag["n_ct_dirs"] = n_ct
        if not cols:
            self._diag["n_eff_dirs"] = 0
            return None
        # The sex (eigh) and celltype (svd) bases are each internally orthonormal, but a sex direction
        # and a celltype direction can OVERLAP (the sex ratio differs by celltype). A bare QR keeps the
        # rank-deficient column as an arbitrary NOISE direction and projects it out of z — wasting
        # capacity and over-counting the rank. SVD is rank-revealing: keep only left-singular vectors
        # with non-negligible singular value ⇒ an orthonormal basis of the TRUE union span, no phantoms.
        M = torch.stack(cols, dim=1)                        # [d_z, n_sex+n_ct] candidate nuisance dirs
        U, S, _ = torch.linalg.svd(M, full_matrices=False)  # S descending
        tol = float(S[0]) * 1e-4 if S.numel() and float(S[0]) > 0 else 0.0
        keep = int((S > tol).sum())
        self._diag["n_eff_dirs"] = keep                     # rank ACTUALLY removed (≤ n_sex+n_ct after dedup)
        if keep == 0:
            return None
        return U[:, :keep]                                  # [d_z, keep] orthonormal nuisance basis

    def forward(self, z_perp, *, sex_id=None, celltype_id=None, donor_id=None, update: bool = True):
        if not self.cfg.enabled:
            return z_perp
        if update and self.training and celltype_id is not None:
            self._update(z_perp.detach().float(), sex_id, celltype_id)
        with torch.autocast(device_type=z_perp.device.type, enabled=False):
            Q = self._basis(torch.float32)
            if Q is None:
                return z_perp
            zf = z_perp.float()
            proj = (zf @ Q) @ Q.t()                         # the nuisance component
            z_clean = zf - proj
            with torch.no_grad():                           # cheap diagnostics for the log / refill watch
                tot = float((zf * zf).sum()) + 1e-12
                self._diag["removed_energy_frac"] = float((proj * proj).sum()) / tot
        return z_clean.to(z_perp.dtype)

    @torch.no_grad()
    def whitening_target(
        self,
        *,
        device: torch.device | str | None = None,
        dtype: torch.dtype = torch.float32,
    ) -> torch.Tensor:
        """Return the projector-valued covariance target used after cleanup.

        The cleaned row-vector latent is ``z_clean = z_perp @ P`` with
        ``P = I - Q Q^T``.  Consequently, whitening must target ``P`` rather
        than the full-space identity.  Recomputing ``Q`` here is deterministic:
        the forward updates the EMA buffers before calling ``_basis`` and no
        projector update occurs between the model forward and the loss call.

        When projection is disabled, or the estimated nuisance basis is empty,
        this returns the identity and therefore preserves the historical loss.
        The returned tensor is detached by construction.
        """

        target_device = (
            self.sex_scatter.device
            if device is None
            else torch.device(device)
        )
        # The trainer calls this method from its AMP autocast region.  Without
        # an explicit guard, CUDA executes Q @ Q.T in bf16 even though Q and
        # eye are fp32 tensors.  The resulting matrix can fail the downstream
        # symmetry/idempotence/eigenvalue audit despite representing the right
        # subspace.  Keep the complete projector construction in fp32, exactly
        # as forward() already protects the projection applied to z_perp.
        with torch.autocast(
            device_type=self.sex_scatter.device.type,
            enabled=False,
        ):
            eye = torch.eye(
                self.d_z,
                dtype=torch.float32,
                device=self.sex_scatter.device,
            )
            Q = self._basis(torch.float32) if bool(self.cfg.enabled) else None
            if Q is None:
                target = eye
                removed_rank = 0
            else:
                target = eye - Q @ Q.t()
                # Numerical symmetrisation prevents tiny eigh/SVD round-off
                # from appearing as an asymmetric whitening target downstream.
                target = 0.5 * (target + target.t())
                removed_rank = int(Q.shape[1])
        self._diag["whitening_removed_rank"] = removed_rank
        self._diag["whitening_effective_dim"] = self.d_z - removed_rank
        return target.to(device=target_device, dtype=dtype)

    @torch.no_grad()
    def diagnostics(self) -> dict:
        """Snapshot for the training log: removed-energy fraction, rank actually removed, sex spectrum."""
        return dict(self._diag)
