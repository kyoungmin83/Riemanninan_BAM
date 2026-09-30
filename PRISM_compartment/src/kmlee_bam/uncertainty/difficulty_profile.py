"""
Decoder difficulty profile — per-(gene, bin, celltype) median/MAD of rec_err.

This module computes the per-cell decoder difficulty baseline against which we
normalize reconstruction error to derive `noise_score`. The goal is to strip
out structural bias (bin 0 is easy, middle bins are hard, certain genes are
inherently noisy) so that what remains is closer to true measurement noise.

See: model_v7a_implementation_plan.md §4.2.2, model_v7_design...md §5.0
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.distributed as dist


# Fallback level identifiers. The values are *keys* into the profile dict; the
# semantics are documented in §5.0.3 of the v7 design doc and §4.2.2 of the
# implementation plan.
LEVEL_GBT = ("gene", "bin", "celltype")
LEVEL_GB = ("gene", "bin")
LEVEL_BT = ("bin", "celltype")
LEVEL_G = ("gene",)
LEVEL_B = ("bin",)

DEFAULT_FALLBACK_LEVELS: Tuple[Tuple[str, ...], ...] = (
    LEVEL_GBT,
    LEVEL_GB,
    LEVEL_BT,
    LEVEL_G,
    LEVEL_B,
)


@dataclass
class DifficultyProfileConfig:
    min_count: int = 100
    mad_floor: float = 0.01
    winsorize_z: float = 5.0
    mad_to_std: float = 1.4826
    gene_chunk_size: int = 500
    ema_momentum: float = 0.95
    # Cap on the number of (cell, gene) values retained per Level-3 (b, t)
    # bucket or Level-5 (b) bucket. We use a reservoir-style accumulator so
    # peak CPU memory is bounded at ~(n_buckets × max_bucket_samples × 4 B).
    # 100k per bucket gives median CI of ±0.3% (more than enough for fallback).
    max_bucket_samples: int = 100_000
    fallback_levels: Tuple[Tuple[str, ...], ...] = field(
        default_factory=lambda: DEFAULT_FALLBACK_LEVELS
    )


def _median_mad(values: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """Return (median, MAD) along dim 0 of a 1-D or 2-D tensor."""
    if values.numel() == 0:
        return (
            torch.tensor(0.0, dtype=values.dtype, device=values.device),
            torch.tensor(0.0, dtype=values.dtype, device=values.device),
        )
    med = values.median(dim=0).values
    mad = (values - med.unsqueeze(0)).abs().median(dim=0).values
    return med, mad


class _Reservoir:
    """
    Online reservoir for per-bucket value sampling with bounded memory.

    Each `add()` ingests a new chunk of values. If keeping them all would
    exceed `max_size`, we mix old+new and uniformly subsample back down to
    `max_size`. `seen_total` tracks the true total count so callers can
    record the bucket's true cell count for diagnostics.
    """

    __slots__ = ("max_size", "values", "seen_total")

    def __init__(self, max_size: int):
        self.max_size = int(max_size)
        self.values: Optional[torch.Tensor] = None
        self.seen_total: int = 0

    def add(self, vals: torch.Tensor) -> None:
        n = int(vals.numel())
        if n == 0:
            return
        self.seen_total += n
        if self.values is None:
            if n <= self.max_size:
                self.values = vals.clone()
            else:
                idx = torch.randperm(n)[: self.max_size]
                self.values = vals[idx].clone()
            return
        if self.values.numel() + n <= self.max_size:
            self.values = torch.cat([self.values, vals])
            return
        merged = torch.cat([self.values, vals])
        idx = torch.randperm(merged.numel())[: self.max_size]
        self.values = merged[idx]

    def get(self) -> torch.Tensor:
        if self.values is None:
            return torch.empty(0, dtype=torch.float32)
        return self.values


class DifficultyProfile:
    """
    Holds per-level (mu, mad, n) tensors and supports vectorized fallback
    lookup. All tensors live on CPU during compute and are .to(device) when
    used inside normalize_rec_err.

    Shape conventions:
        n_genes      : G
        n_bins       : K (e.g. 7)
        n_celltypes  : T
    """

    def __init__(
        self,
        config: DifficultyProfileConfig,
        n_genes: int,
        n_bins: int,
        n_celltypes: int,
        device: Optional[torch.device] = None,
    ):
        self.config = config
        self.n_genes = int(n_genes)
        self.n_bins = int(n_bins)
        self.n_celltypes = int(n_celltypes)
        self.device = device if device is not None else torch.device("cpu")
        self.fallback_levels: Tuple[Tuple[str, ...], ...] = tuple(
            tuple(level) for level in config.fallback_levels
        )

        # Per-level storage: dict mapping level-key -> dict with 'mu', 'mad', 'n'
        # tensors. Shapes depend on the level:
        #   (gene, bin, celltype): [G, K, T]
        #   (gene, bin)          : [G, K]
        #   (bin, celltype)      : [K, T]
        #   (gene,)              : [G]
        #   (bin,)               : [K]
        self._profiles: dict[Tuple[str, ...], dict[str, torch.Tensor]] = {}
        self._initialized = False

    # ------------------------------------------------------------------ #
    # Computation
    # ------------------------------------------------------------------ #

    @torch.no_grad()
    def compute_from_sample(
        self,
        rec_err: torch.Tensor,        # [N, G]
        true_bins: torch.Tensor,      # [N, G] int64 in [0, K)
        celltypes: torch.Tensor,      # [N] int64 in [0, T)
    ) -> dict:
        """
        Build all fallback-level profiles from a sampled subset of cells.

        Strategy: chunk over genes (bounded peak memory) and use vectorized
        nanquantile per (bin, [celltype]) inside each chunk. This replaces a
        previous implementation that had a Python-level per-gene loop, which
        was O(G * K * T) Python overhead and effectively unusable at the
        50k × 11k scale.

        Returns: stats dict for logging.
        """
        if rec_err.ndim != 2:
            raise ValueError(f"rec_err must be [N, G], got {tuple(rec_err.shape)}")
        if true_bins.shape != rec_err.shape:
            raise ValueError(
                f"true_bins must match rec_err shape, got {tuple(true_bins.shape)} "
                f"vs {tuple(rec_err.shape)}"
            )
        if celltypes.shape[0] != rec_err.shape[0]:
            raise ValueError(
                f"celltypes batch dim must match rec_err, got "
                f"{tuple(celltypes.shape)} vs {tuple(rec_err.shape)}"
            )

        N, G = rec_err.shape
        K = self.n_bins
        T = self.n_celltypes

        # Keep large tensors (rec_err, true_bins) on CPU; move per-gene-chunk
        # subsets to device inside the loop below. This caps peak GPU memory
        # at gene_chunk_size * N * 4 bytes (~100 MB for 50k cells × 500 genes)
        # instead of ~6 GB for the full [50k, 11k] tensors.
        rec_err = rec_err.float()           # CPU
        # int8 storage for bin/celltype tags — equality masks (`tb == b`)
        # work directly with int8, and we save ~8× CPU RAM vs int64.
        # That keeps CPU peak well under 10 GB even at sample_size=50k.
        if true_bins.dtype != torch.int8:
            true_bins = true_bins.to(torch.int8)
        celltypes_cpu = (
            celltypes if celltypes.dtype == torch.int8 else celltypes.to(torch.int8)
        )
        celltypes_dev = celltypes_cpu.to(self.device)

        # Pre-allocated outputs.
        zeros_gkt = lambda: torch.zeros(G, K, T, dtype=torch.float32, device=self.device)
        zeros_gk = lambda: torch.zeros(G, K, dtype=torch.float32, device=self.device)
        zeros_kt = lambda: torch.zeros(K, T, dtype=torch.float32, device=self.device)
        zeros_g = lambda: torch.zeros(G, dtype=torch.float32, device=self.device)
        zeros_k = lambda: torch.zeros(K, dtype=torch.float32, device=self.device)

        prof_gbt = {"mu": zeros_gkt(), "mad": zeros_gkt(), "n": zeros_gkt()}
        prof_gb = {"mu": zeros_gk(), "mad": zeros_gk(), "n": zeros_gk()}
        prof_bt = {"mu": zeros_kt(), "mad": zeros_kt(), "n": zeros_kt()}
        prof_g = {"mu": zeros_g(), "mad": zeros_g(), "n": zeros_g()}
        prof_b = {"mu": zeros_k(), "mad": zeros_k(), "n": zeros_k()}

        chunk = max(1, int(self.config.gene_chunk_size))
        nan_val = float("nan")

        # Accumulators for Level 3 (bin, celltype) and Level 5 (bin). These
        # cross gene-chunk boundaries; use reservoir-style sampling so peak
        # CPU memory is bounded at ~(n_buckets × max_bucket_samples × 4 B).
        # For 168 (b,t) + 7 (b) buckets × 100k samples × 4 B ≈ 70 MB total —
        # huge improvement over the previous "keep everything" approach
        # that hit ~10 GB peak at 50k cells × 11k genes.
        max_bucket = int(self.config.max_bucket_samples)
        _accumulate_bt: "dict[tuple[int, int], _Reservoir]" = {}
        _accumulate_b: "dict[int, _Reservoir]" = {}

        # Iterate over gene chunks to bound memory. We move only this chunk
        # to device so peak memory is bounded.
        for g0 in range(0, G, chunk):
            g1 = min(G, g0 + chunk)
            re_chunk = rec_err[:, g0:g1].to(self.device)    # [N, gc] on device
            tb_chunk = true_bins[:, g0:g1].to(self.device)  # [N, gc] on device
            gc = g1 - g0

            # ---- Level 1: (gene, bin, celltype) ---- #
            # For each (b, t), build a mask over [N, gc] and compute per-gene
            # median/MAD via nanquantile. This vectorizes the previously
            # per-gene Python loop.
            for b in range(K):
                bin_mask_chunk = tb_chunk == b        # [N, gc]
                if not bool(bin_mask_chunk.any()):
                    continue
                for t in range(T):
                    ct_mask = celltypes_dev == t      # [N] on device
                    if not bool(ct_mask.any()):
                        continue
                    mask = bin_mask_chunk & ct_mask.unsqueeze(-1)   # [N, gc]
                    if not bool(mask.any()):
                        continue
                    masked = torch.where(
                        mask, re_chunk, torch.full_like(re_chunk, nan_val)
                    )
                    n_per_gene = mask.sum(dim=0).float()             # [gc]
                    med = torch.nanquantile(masked, 0.5, dim=0)      # [gc]
                    deviation = (masked - med.unsqueeze(0)).abs()
                    mad = torch.nanquantile(deviation, 0.5, dim=0)   # [gc]
                    valid = n_per_gene > 0
                    prof_gbt["mu"][g0:g1, b, t] = torch.where(
                        valid, med, torch.zeros_like(med)
                    )
                    prof_gbt["mad"][g0:g1, b, t] = torch.where(
                        valid, mad, torch.zeros_like(mad)
                    )
                    prof_gbt["n"][g0:g1, b, t] = n_per_gene

            # ---- Level 2: (gene, bin), celltype-marginal ---- #
            for b in range(K):
                bin_mask_chunk = tb_chunk == b
                if not bool(bin_mask_chunk.any()):
                    continue
                masked = torch.where(
                    bin_mask_chunk, re_chunk, torch.full_like(re_chunk, nan_val)
                )
                n_per_gene = bin_mask_chunk.sum(dim=0).float()
                med = torch.nanquantile(masked, 0.5, dim=0)
                deviation = (masked - med.unsqueeze(0)).abs()
                mad = torch.nanquantile(deviation, 0.5, dim=0)
                valid = n_per_gene > 0
                prof_gb["mu"][g0:g1, b] = torch.where(valid, med, torch.zeros_like(med))
                prof_gb["mad"][g0:g1, b] = torch.where(valid, mad, torch.zeros_like(mad))
                prof_gb["n"][g0:g1, b] = n_per_gene

            # ---- Level 4: (gene,) bin+celltype-marginal ---- #
            # No masking; all rows contribute to every gene.
            med = re_chunk.median(dim=0).values
            deviation = (re_chunk - med.unsqueeze(0)).abs()
            mad = deviation.median(dim=0).values
            prof_g["mu"][g0:g1] = med
            prof_g["mad"][g0:g1] = mad
            prof_g["n"][g0:g1] = float(re_chunk.shape[0])

            # ---- Level 3 + 5 reservoir accumulators ---- #
            # For each (b, t) and each (b) bucket, push the chunk's relevant
            # values into a fixed-size reservoir on CPU. Memory per bucket
            # is bounded; per-bin true counts are tracked separately via
            # `Reservoir.seen_total`.
            for b in range(K):
                bm_b = tb_chunk == b
                if not bool(bm_b.any()):
                    continue
                for t in range(T):
                    ct_mask = celltypes_dev == t
                    joint = bm_b & ct_mask.unsqueeze(-1)
                    if not bool(joint.any()):
                        continue
                    vals_chunk = re_chunk[joint].detach().cpu()
                    res = _accumulate_bt.get((b, t))
                    if res is None:
                        res = _Reservoir(max_bucket)
                        _accumulate_bt[(b, t)] = res
                    res.add(vals_chunk)
                # Per-bin reservoir (Level 5)
                vals_b_chunk = re_chunk[bm_b].detach().cpu()
                res_b = _accumulate_b.get(b)
                if res_b is None:
                    res_b = _Reservoir(max_bucket)
                    _accumulate_b[b] = res_b
                res_b.add(vals_b_chunk)

        # Finalize Level 3 + 5 from reservoirs. We report the *true* total
        # count (seen_total) for diagnostic n, but compute median/MAD on the
        # reservoir-subsampled tensor for bounded compute time.
        for (b, t), res in _accumulate_bt.items():
            vals = res.get()
            if vals.numel() == 0:
                continue
            med, mad = _median_mad(vals)
            prof_bt["mu"][b, t] = med.to(self.device)
            prof_bt["mad"][b, t] = mad.to(self.device)
            prof_bt["n"][b, t] = float(res.seen_total)

        for b, res in _accumulate_b.items():
            vals = res.get()
            if vals.numel() == 0:
                continue
            med, mad = _median_mad(vals)
            prof_b["mu"][b] = med.to(self.device)
            prof_b["mad"][b] = mad.to(self.device)
            prof_b["n"][b] = float(res.seen_total)

        self._profiles = {
            LEVEL_GBT: prof_gbt,
            LEVEL_GB: prof_gb,
            LEVEL_BT: prof_bt,
            LEVEL_G: prof_g,
            LEVEL_B: prof_b,
        }
        self._initialized = True

        stats = self._summary_stats(N)
        return stats

    @torch.no_grad()
    def ema_update(self, other: "DifficultyProfile") -> None:
        """
        Exponentially update this profile using the freshly-computed `other`.

        For each level/key:
            new = momentum * old + (1 - momentum) * other
        Cells with `n == 0` in `other` are left untouched.
        """
        if not self._initialized:
            raise RuntimeError("Profile must be initialized before ema_update.")
        if not other._initialized:
            raise RuntimeError("Other profile must be initialized before ema_update.")
        mom = float(self.config.ema_momentum)
        for key, prof in self._profiles.items():
            o = other._profiles[key]
            valid = o["n"] > 0
            if not bool(valid.any()):
                continue
            for stat in ("mu", "mad", "n"):
                target = prof[stat]
                source = o[stat]
                target[valid] = mom * target[valid] + (1.0 - mom) * source[valid]

    # ------------------------------------------------------------------ #
    # Lookup
    # ------------------------------------------------------------------ #

    @torch.no_grad()
    def lookup_vectorized(
        self,
        gene_idx: torch.Tensor,    # [...] long
        bin_idx: torch.Tensor,     # [...] long
        celltype_idx: torch.Tensor,  # [...] long, broadcastable to gene/bin shape
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Vectorized fallback lookup. Returns (mu, mad, level_idx) tensors
        broadcast to gene_idx.shape.

            level_idx == 0: (g, b, t)
            level_idx == 1: (g, b)
            level_idx == 2: (b, t)
            level_idx == 3: (g,)
            level_idx == 4: (b,)

        Fallback rule: a level is "valid" iff its n >= min_count AND mad > floor.
        If no level is valid (extreme corner case), we still use the last
        level's value so caller never gets a NaN/Inf.
        """
        if not self._initialized:
            raise RuntimeError("Profile not initialized; call compute_from_sample first.")

        device = gene_idx.device
        target_shape = torch.broadcast_shapes(
            gene_idx.shape, bin_idx.shape, celltype_idx.shape
        )

        gene_idx = gene_idx.to(device).long().expand(target_shape).reshape(-1)
        bin_idx = bin_idx.to(device).long().expand(target_shape).reshape(-1)
        celltype_idx = celltype_idx.to(device).long().expand(target_shape).reshape(-1)

        min_count = float(self.config.min_count)
        mad_floor = float(self.config.mad_floor)

        # Initialize output with sentinel; we'll fill from the first valid level.
        mu_out = torch.zeros_like(gene_idx, dtype=torch.float32)
        mad_out = torch.full_like(gene_idx, mad_floor, dtype=torch.float32)
        level_out = torch.full_like(
            gene_idx, fill_value=len(self.fallback_levels) - 1, dtype=torch.long
        )
        filled = torch.zeros_like(gene_idx, dtype=torch.bool)

        for level_pos, level_key in enumerate(self.fallback_levels):
            prof = self._profiles[level_key]
            mu_t = prof["mu"].to(device)
            mad_t = prof["mad"].to(device)
            n_t = prof["n"].to(device)

            if level_key == LEVEL_GBT:
                mu = mu_t[gene_idx, bin_idx, celltype_idx]
                mad = mad_t[gene_idx, bin_idx, celltype_idx]
                n = n_t[gene_idx, bin_idx, celltype_idx]
            elif level_key == LEVEL_GB:
                mu = mu_t[gene_idx, bin_idx]
                mad = mad_t[gene_idx, bin_idx]
                n = n_t[gene_idx, bin_idx]
            elif level_key == LEVEL_BT:
                mu = mu_t[bin_idx, celltype_idx]
                mad = mad_t[bin_idx, celltype_idx]
                n = n_t[bin_idx, celltype_idx]
            elif level_key == LEVEL_G:
                mu = mu_t[gene_idx]
                mad = mad_t[gene_idx]
                n = n_t[gene_idx]
            elif level_key == LEVEL_B:
                mu = mu_t[bin_idx]
                mad = mad_t[bin_idx]
                n = n_t[bin_idx]
            else:
                raise ValueError(f"Unknown fallback level key: {level_key}")

            # Fallback validity: only require min_count. A tight (small MAD)
            # but well-populated combo is *more* informative, not less; the
            # mad_floor parameter still guards against division-by-zero
            # downstream when MAD is exactly 0.
            valid = (~filled) & (n >= min_count)
            if bool(valid.any()):
                mu_out[valid] = mu[valid]
                mad_out[valid] = mad[valid].clamp_min(mad_floor)
                level_out[valid] = level_pos
                filled |= valid

            if bool(filled.all()):
                break

        # Anything still unfilled gets the final-level value (which may be
        # below threshold but is at least defined).
        if not bool(filled.all()):
            last_level = self.fallback_levels[-1]
            prof = self._profiles[last_level]
            mu_t = prof["mu"].to(device)
            mad_t = prof["mad"].to(device)
            if last_level == LEVEL_B:
                mu_fallback = mu_t[bin_idx]
                mad_fallback = mad_t[bin_idx]
            else:
                # Defensive; the doc fixes LEVEL_B as the final fallback.
                mu_fallback = mu_t.flatten()[0].expand_as(gene_idx).float()
                mad_fallback = mad_t.flatten()[0].expand_as(gene_idx).float()
            unfilled = ~filled
            mu_out[unfilled] = mu_fallback[unfilled]
            mad_out[unfilled] = mad_fallback[unfilled].clamp_min(mad_floor)
            level_out[unfilled] = len(self.fallback_levels) - 1

        mu_out = mu_out.reshape(target_shape)
        mad_out = mad_out.reshape(target_shape).clamp_min(mad_floor)
        level_out = level_out.reshape(target_shape)
        return mu_out, mad_out, level_out

    # ------------------------------------------------------------------ #
    # State dict (checkpoint persistence)
    # ------------------------------------------------------------------ #

    @torch.no_grad()
    def broadcast_(self, src_rank: int = 0) -> None:
        """
        DDP broadcast: copy the profile state from `src_rank` to all other
        ranks.

        Required because profile init runs on rank 0 only (to avoid each rank
        building a profile from a different DistributedSampler shard, which
        would silently desync the normalization basis across ranks).

        No-op when not running under torch.distributed.
        """
        if not (dist.is_available() and dist.is_initialized()):
            return

        rank = dist.get_rank()

        # Step 1: broadcast "is the source initialized?" flag.
        flag = torch.tensor(
            [1.0 if (rank == src_rank and self._initialized) else 0.0],
            device=self.device,
            dtype=torch.float32,
        )
        dist.broadcast(flag, src=src_rank)
        src_initialized = float(flag.item()) >= 0.5
        if not src_initialized:
            return

        # Step 2: ensure all ranks have correctly-shaped buffers, then broadcast.
        for level_key in self.fallback_levels:
            if level_key not in self._profiles:
                self._profiles[level_key] = self._allocate_empty_level(level_key)
            for stat_key in ("mu", "mad", "n"):
                dist.broadcast(self._profiles[level_key][stat_key], src=src_rank)

        self._initialized = True

    def _allocate_empty_level(self, level_key: Tuple[str, ...]) -> dict:
        """Pre-allocate empty tensors of the correct shape for a given level."""
        G, K, T = self.n_genes, self.n_bins, self.n_celltypes
        if level_key == LEVEL_GBT:
            shape = (G, K, T)
        elif level_key == LEVEL_GB:
            shape = (G, K)
        elif level_key == LEVEL_BT:
            shape = (K, T)
        elif level_key == LEVEL_G:
            shape = (G,)
        elif level_key == LEVEL_B:
            shape = (K,)
        else:
            raise ValueError(f"Unknown level key: {level_key}")
        return {
            "mu": torch.zeros(shape, dtype=torch.float32, device=self.device),
            "mad": torch.zeros(shape, dtype=torch.float32, device=self.device),
            "n": torch.zeros(shape, dtype=torch.float32, device=self.device),
        }

    def state_dict(self) -> dict:
        out = {"initialized": bool(self._initialized)}
        for key, prof in self._profiles.items():
            tag = "_".join(key)
            out[f"{tag}_mu"] = prof["mu"].detach().cpu()
            out[f"{tag}_mad"] = prof["mad"].detach().cpu()
            out[f"{tag}_n"] = prof["n"].detach().cpu()
        return out

    def load_state_dict(self, state: dict) -> None:
        self._initialized = bool(state.get("initialized", False))
        if not self._initialized:
            return
        self._profiles = {}
        for level_key in self.fallback_levels:
            tag = "_".join(level_key)
            mu = state[f"{tag}_mu"].to(self.device)
            mad = state[f"{tag}_mad"].to(self.device)
            n = state[f"{tag}_n"].to(self.device)
            self._profiles[level_key] = {"mu": mu, "mad": mad, "n": n}

    # ------------------------------------------------------------------ #
    # Diagnostics
    # ------------------------------------------------------------------ #

    def _summary_stats(self, n_sample_cells: int) -> dict:
        """Return diagnostic stats about populated bins and fallback hits."""
        stats: dict = {"n_sample_cells": int(n_sample_cells)}
        min_count = float(self.config.min_count)
        mad_floor = float(self.config.mad_floor)
        for level_key, prof in self._profiles.items():
            tag = "_".join(level_key)
            n_t = prof["n"]
            mad_t = prof["mad"]
            n_total = int(n_t.numel())
            n_populated = int((n_t > 0).sum())
            # Validity (diagnostic): mirror the fallback rule — n >= min_count only.
            n_valid = int((n_t >= min_count).sum())
            stats[f"{tag}/cells_total"] = n_total
            stats[f"{tag}/cells_populated"] = n_populated
            stats[f"{tag}/cells_valid"] = n_valid
            if n_populated > 0:
                stats[f"{tag}/mean_n"] = float(n_t[n_t > 0].mean())
                stats[f"{tag}/median_mad"] = float(mad_t[n_t > 0].median())
        return stats

    @property
    def initialized(self) -> bool:
        return self._initialized
