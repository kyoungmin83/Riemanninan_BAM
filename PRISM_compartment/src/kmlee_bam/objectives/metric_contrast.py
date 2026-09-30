"""v35 metric-contrast objective — pressure to bend the latent metric where disease is.

The disease-GATED decoder (`decoder.use_disease_gate`) gives the *capacity* to make
``G_disease(z) != G_normal(z)``; this objective gives the *pressure* to make that bending
amyloid-aligned, donor-robust, and anisotropic. Per step, on the DECODER score (which includes the
gate — `_compute_score_components` routes through `pathology_state_and_p`):

  q(z, u) = || score(z + eps·u) - score(z - eps·u) ||^2 / (2 eps)^2      # directional sensitivity
            = an L2 proxy for the Fisher norm ||J(z) u||^2_F. The proper Fisher reweighting is the
            EVAL's job (disease_metric_heldout / r_gf_fingerprint_v2); here we only need the CONTRAST.

We want the AMYLOID-direction sensitivity LARGER in amyloid-high donors than normal donors
(G_disease != G_normal, aligned to amyloid), measured DONOR-BALANCED (per-donor means, equal
weight), and BEYOND a donor-label PERMUTATION baseline so it can't be a donor-memorization artifact.
A random-direction CONTRAST term (not magnitude — so it never fights reconstruction) keeps the
disease anisotropy amyloid-specific rather than isotropic.

  contrast(dir) = mean_{hi donors} q_bar_d(dir)  -  mean_{lo donors} q_bar_d(dir)
  loss = relu(margin - (contrast(u_amyloid) - E_perm[contrast under SHUFFLED hi/lo]))
         + lambda_nuisance * relu(contrast(u_random))

Default disabled (enabled=False / lambda 0) => the head is never called => byte-identical.
NOTE (v1): L2 proxy for Fisher + soft perm-baseline; the real donor-perm gate is the held-out eval.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class MetricContrastConfig:
    enabled: bool = False
    lambda_metric_contrast: float = 0.0
    # amyloid axis: index into the decoder pathology_head rows (DEFAULT_PATHOLOGY_AXES =
    # [Braak, Thal, LATE, Lewy]) and the donor-level label key fed in the batch.
    amyloid_axis_idx: int = 1
    amyloid_label_key: str = "axis_thal_high"
    eps: float = 0.1                 # finite-difference step along the direction
    margin: float = 0.0              # require contrast_real - contrast_perm >= margin
    n_perm: int = 8                  # donor-label shuffles for the permutation baseline
    min_donors: int = 3              # need >= this many hi AND lo donors in the batch, else skip
    lambda_nuisance: float = 1.0     # weight on the random-direction CONTRAST penalty (anisotropy)
    n_nuisance_dirs: int = 1
    detach_direction: bool = True    # move the GATE toward the amyloid dir, not the dir toward easy curvature


@dataclass
class MetricContrastOutput:
    loss: torch.Tensor
    contrast_real: float = 0.0
    contrast_perm: float = 0.0
    contrast_nuisance: float = 0.0
    q_hi: float = 0.0
    q_lo: float = 0.0
    n_hi: int = 0
    n_lo: int = 0
    active: bool = False


class MetricContrastHead(nn.Module):
    """Parameter-free: reads the decoder's (tied) amyloid direction and bends its metric via grads."""

    def __init__(self, config: MetricContrastConfig) -> None:
        super().__init__()
        self.cfg = config

    def _dirnorm(self, decoder, z, ct, tech, direction):
        # review #4: suppress the decoder's diagnostic stashes (_last_*_absmean) during these
        # finite-diff probes so the perturbed values don't overwrite the REAL forward's logs.
        prev = getattr(decoder, "_stash_diag", True)
        decoder._stash_diag = False
        try:
            e = self.cfg.eps * direction                                    # [d_z]
            sp = decoder._compute_score_components(z + e, ct, tech)[0]      # [B, G] (incl. gate)
            sm = decoder._compute_score_components(z - e, ct, tech)[0]
            d = (sp - sm) / (2.0 * self.cfg.eps)                            # [B, G] ≈ J·dir
            # mean (not sum) over genes => scale-stable across n_genes (L2 proxy for Fisher)
            return (d * d).mean(dim=1)                                      # [B] per-gene-mean ‖J·dir‖²
        finally:
            decoder._stash_diag = prev

    @staticmethod
    def _donor_means(q, donor_id):
        uniq = torch.unique(donor_id)
        qs, ds = [], []
        for d in uniq:
            m = donor_id == d
            if int(m.sum()) > 0:
                qs.append(q[m].mean()); ds.append(d)
        return torch.stack(qs), torch.stack(ds)                            # [D], [D]

    def forward(self, *, decoder, z_perp, celltype_id, tech_id, donor_id, amyloid_label):
        cfg = self.cfg
        dev = z_perp.device
        zero = z_perp.new_zeros(())
        ph = getattr(decoder, "pathology_head", None)
        if ph is None or amyloid_label is None:
            return MetricContrastOutput(loss=zero)
        # amyloid direction in z (tied to the supervised Thal axis; detached so the GATE moves)
        u = ph.weight[cfg.amyloid_axis_idx]
        if cfg.detach_direction:
            u = u.detach()
        u = u / (u.norm() + 1e-8)

        don = donor_id.detach()
        q_amy = self._dirnorm(decoder, z_perp, celltype_id, tech_id, u)     # [B]
        qd, dvocab = self._donor_means(q_amy, don)                         # [D]
        # donor-level amyloid hi/lo (label is donor-constant; take per-donor mean then threshold)
        lab = amyloid_label.detach().float()
        labd = self._donor_means(lab, don)[0]                              # [D] ~0/1
        hi = labd >= 0.5
        lo = ~hi
        n_hi, n_lo = int(hi.sum()), int(lo.sum())
        if n_hi < cfg.min_donors or n_lo < cfg.min_donors:
            return MetricContrastOutput(loss=zero, n_hi=n_hi, n_lo=n_lo)

        contrast_real = qd[hi].mean() - qd[lo].mean()
        # permutation baseline: shuffle which donors are "hi" (differentiable through qd)
        D = qd.shape[0]
        perms = []
        for _ in range(cfg.n_perm):
            idx = torch.randperm(D, device=dev)
            hm = torch.zeros(D, dtype=torch.bool, device=dev); hm[idx[:n_hi]] = True
            perms.append(qd[hm].mean() - qd[~hm].mean())
        contrast_perm = torch.stack(perms).mean()
        loss = F.relu(cfg.margin - (contrast_real - contrast_perm))

        # nuisance: a RANDOM direction must NOT show disease contrast (keeps anisotropy amyloid-
        # specific). Penalize the contrast, not the magnitude => never fights reconstruction.
        c_nz_val = 0.0
        for _ in range(max(1, cfg.n_nuisance_dirs)):
            r = torch.randn_like(u); r = r / (r.norm() + 1e-8)
            q_r = self._dirnorm(decoder, z_perp, celltype_id, tech_id, r)
            qrd = self._donor_means(q_r, don)[0]
            c_r = qrd[hi].mean() - qrd[lo].mean()
            # review #2: SYMMETRIC penalty (c_r²) — a random direction must have ~0 disease contrast
            # in EITHER sign, else the model gets free isotropic anisotropy (normal>disease) unblocked.
            loss = loss + (cfg.lambda_nuisance / max(1, cfg.n_nuisance_dirs)) * c_r.pow(2)
            c_nz_val += float(c_r.detach())

        return MetricContrastOutput(
            loss=loss,
            contrast_real=float(contrast_real.detach()),
            contrast_perm=float(contrast_perm.detach()),
            contrast_nuisance=c_nz_val / max(1, cfg.n_nuisance_dirs),
            q_hi=float(qd[hi].mean().detach()), q_lo=float(qd[lo].mean().detach()),
            n_hi=n_hi, n_lo=n_lo, active=True,
        )
