from __future__ import annotations

"""
Donor x Region block gradient SNR audit  — DIAGNOSTIC ONLY.  (v2)

Question
--------
Does each model BRANCH learn a direction that AGREES across DONORS
(-> generalizable population / disease signal) or one that swings donor-to-donor
(-> donor-specific memorization / artifact)?  scRNA cells are NOT i.i.d.; the
*donor* is the exchangeable generalization unit, so agreement is measured across
donors, not across cells in a minibatch (Litman & Guo 2026, with the donor as the
correct exchangeable block).

Why donor x REGION blocks (not bare donor)
------------------------------------------
In SEA-AD the SAME donor contributes BOTH DLPFC and MTG cells, so "region of a
donor" is undefined. We therefore sample (donor, region) BLOCKS — each block is a
single donor in a single region. This lets us cleanly separate two questions:
  * within-region cross-donor SNR  = generalization (do donors agree, region fixed)
  * cross-region cosine            = region leakage (do DLPFC and MTG pull a branch
                                     the same way?)  low / negative => region leaks
                                     into that branch's learning. Relevant because
                                     this DLPFC+MTG config has NO explicit region
                                     term in `base` (celltype-level only).

Metrics  (per branch)
---------------------
For block b, g_b = unit-normalized branch gradient (so a numerous donor does not
dominate).
  within-region consensus_r = || mean_{b in region r} g_b ||^2      in [0,1]
  SNR(branch)               = mean over regions of consensus_r       in [0,1]
  region_cos(branch)        = cos( mean_DLPFC g_b , mean_MTG g_b )   in [-1,1]
1.0 SNR = donors agree (signal); ~0 = cancel (noise).
region_cos ~1 = regions agree (clean);  <~0 = region leaks into the branch.

Safety / correctness  (the v1 fixes)
-----------------------------------
* STATELESS loss: the audit gradient is the plain ordinal NLL computed directly
  from `out.decoder_out.probs` and `y_ord` — it does NOT call the trainer's
  _compute_loss, so it never advances step counters / running EMAs / schedules.
* DDP-UNWRAP: uses `getattr(model, "module", model)` so forward/backward does not
  trigger all-reduce (which would average donors across GPUs and erase the signal).
  Run on a SINGLE process (rank 0) -- the caller guards rank.
* Pure read of `.grad`; grads set to None before and after each block -> never
  pollutes a real optimizer step.
* Default eval() mode (no dropout) -> deterministic per-block gradient direction.
* fp32.  Forward uses sample_latent=False (z = mean) to drop sampling noise.

Caveat: high SNR alone is NOT proof of real pathology -- a SHARED technical
artifact present in all donors also agrees across donors. Read SNR alongside a
held-out *pathology* score (disease_program), not just loss.
"""

from collections import defaultdict
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch

try:
    from torch.utils.data import default_collate
except ImportError:  # older torch
    from torch.utils.data.dataloader import default_collate


# branch prefixes over system.named_parameters()  (see model/system.py)
DEFAULT_BRANCHES: tuple = (
    "state_encoder",
    "decoder.coeff_head",
    "decoder.tech_baseline",
    "prior",
    "gene_embedding",
    "module_tokenizer",
)


BRANCH_LABELS = {
    "state_encoder": "state_encoder",
    "decoder.coeff_head": "state_decoder",
    "decoder.tech_baseline": "tech_branch",
    "prior": "ref_prior",
    "gene_embedding": "gene_embedding",
    "module_tokenizer": "module_tokenizer",
}


STATE_BRANCHES = {"state_encoder", "decoder.coeff_head", "module_tokenizer", "prior"}


def _branch_label(name: str) -> str:
    return BRANCH_LABELS.get(name, name.split(".")[-1])


def _snr_status(v: float) -> str:
    if not np.isfinite(v):
        return "NA"
    if v >= 0.55:
        return "strong donor-consistent"
    if v >= 0.30:
        return "moderate"
    return "weak / donor-noisy"


def _region_cos_status(branch: str, v: float) -> str:
    if not np.isfinite(v):
        return "NA"
    if branch == "decoder.tech_baseline":
        if v < 0.0:
            return "region-specific tech/base signal; acceptable if confined here"
        if v < 0.30:
            return "weakly shared tech signal"
        return "shared tech direction"
    if v < 0.0:
        return "HIGH region leakage risk"
    if v < 0.30:
        return "WARN possible region leakage"
    if v < 0.65:
        return "moderate region agreement"
    return "low region leakage"


def _preferred_region_pair(regions: Sequence[str]) -> tuple[str, str] | None:
    uniq = sorted({str(r) for r in regions})
    if len(uniq) < 2:
        return None
    if "DLPFC" in uniq and "MTG" in uniq:
        return ("DLPFC", "MTG")
    return (uniq[0], uniq[1])


def branch_param_groups(
    model: torch.nn.Module, prefixes: Sequence[str]
) -> Dict[str, List[torch.nn.Parameter]]:
    groups: Dict[str, List[torch.nn.Parameter]] = {p: [] for p in prefixes}
    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue
        for pre in prefixes:
            if name.startswith(pre):
                groups[pre].append(param)
                break
    return {k: v for k, v in groups.items() if v}


def _branch_unit_grad(params: List[torch.nn.Parameter]) -> Optional[torch.Tensor]:
    parts = [p.grad.detach().reshape(-1).float() for p in params if p.grad is not None]
    if not parts:
        return None
    g = torch.cat(parts)
    n = torch.linalg.vector_norm(g)
    if not torch.isfinite(n) or float(n) <= 0.0:
        return None
    return (g / n).cpu()


def ordinal_nll(probs: torch.Tensor, y_ord: torch.Tensor, eps: float = 1e-8) -> torch.Tensor:
    """Stateless ordinal NLL: -mean log P(true bin).  probs [B,G,K], y_ord [B,G]."""
    y = y_ord.long().clamp_min(0).unsqueeze(-1)
    p = probs.gather(-1, y).squeeze(-1).clamp_min(eps)
    return -(p.log().mean())


def _build_batch(dataset, idxs: Sequence[int]) -> Dict[str, torch.Tensor]:
    return default_collate([dataset[int(i)] for i in idxs])


def make_blocks(
    dataset,
    region_arr: Sequence,
    *,
    n_blocks: int = 16,
    min_cells: int = 64,
    max_cells: int = 256,
    seed: int = 0,
) -> List[dict]:
    """
    Build up to n_blocks donor x region blocks, balanced across regions.
    `region_arr` must be aligned to dataset index order (len == len(dataset)).
    donor labels come from the dataset (V3OrdinalScDataset.donor_ids / donor_vocab).
    Returns list of {"donor","region","batch","n"}.
    """
    rng = np.random.default_rng(int(seed))
    donor_arr = np.asarray(dataset.donor_ids)[np.asarray(dataset.row_idx)]  # aligned to dataset idx
    donor_names = np.asarray(dataset.donor_vocab, dtype=object)
    region_arr = np.asarray([str(r) for r in region_arr], dtype=object)
    n = len(dataset)
    if len(donor_arr) != n or len(region_arr) != n:
        raise ValueError(f"alignment: dataset={n} donor={len(donor_arr)} region={len(region_arr)}")

    block_idx: Dict[tuple, List[int]] = defaultdict(list)
    for i in range(n):
        block_idx[(int(donor_arr[i]), str(region_arr[i]))].append(i)

    keys = [k for k, v in block_idx.items() if len(v) >= min_cells]
    regions = sorted({r for _, r in keys})
    per_region = max(1, n_blocks // max(1, len(regions)))
    chosen: List[tuple] = []
    for reg in regions:
        kk = [k for k in keys if k[1] == reg]
        rng.shuffle(kk)
        chosen += kk[:per_region]
    rng.shuffle(chosen)
    chosen = chosen[: n_blocks]

    blocks: List[dict] = []
    for d, reg in chosen:
        idxs = block_idx[(d, reg)]
        if len(idxs) > max_cells:
            idxs = list(rng.choice(np.asarray(idxs), size=max_cells, replace=False))
        # Store INDICES only (+ a dataset reference). The (cells x genes) batch tensors are
        # built lazily per block inside run() and freed immediately -> O(1 block) CPU memory,
        # not O(all blocks). Matters most for the in-training hook, which keeps blocks for the run.
        blocks.append(
            {"donor": str(donor_names[d]), "region": reg, "indices": [int(i) for i in idxs],
             "n": len(idxs), "dataset": dataset}
        )
    return blocks


class DonorSNRAudit:
    def __init__(
        self,
        branch_prefixes: Sequence[str] = DEFAULT_BRANCHES,
        sample_latent: bool = False,
        train_mode: bool = False,
        log_fn: Callable[[str], None] = print,
    ) -> None:
        self.branch_prefixes = list(branch_prefixes)
        self.sample_latent = bool(sample_latent)
        self.train_mode = bool(train_mode)
        self._log = log_fn

    def _zero(self, model: torch.nn.Module) -> None:
        for p in model.parameters():
            p.grad = None

    def _block_grads(self, model, block, groups, device) -> Dict[str, Optional[torch.Tensor]]:
        self._zero(model)
        batch = _build_batch(block["dataset"], block["indices"])  # built here, freed on return
        batch = {k: (v.to(device) if torch.is_tensor(v) else v) for k, v in batch.items()}
        out = model(
            batch,
            sample_latent=self.sample_latent,
            return_all_hidden_states=False,
            return_attn_diagnostics=False,
        )
        loss = ordinal_nll(out.decoder_out.probs, batch["y_ord"])
        loss.backward()
        grads = {b: _branch_unit_grad(ps) for b, ps in groups.items()}
        self._zero(model)
        del batch, out, loss
        return grads

    @staticmethod
    def _consensus(unit_grads: List[torch.Tensor]) -> float:
        m = torch.stack(unit_grads).mean(dim=0)
        return float(torch.dot(m, m))

    def run(self, model, blocks: List[dict], device="cuda", epoch: int = 0) -> Dict[str, dict]:
        model = getattr(model, "module", model)  # DDP unwrap -> no all-reduce
        was_training = model.training
        model.train(self.train_mode)
        groups = branch_param_groups(model, self.branch_prefixes)
        if not groups or len(blocks) < 2:
            self._log(f"[diagnostics | donor-SNR | epoch {epoch}] skipped (branches={len(groups)} blocks={len(blocks)})")
            model.train(was_training)
            return {}

        # per branch: unit grads grouped by region
        by_region: Dict[str, Dict[str, List[torch.Tensor]]] = {b: defaultdict(list) for b in groups}
        used = 0
        cells_used = 0
        for blk in blocks:
            grads = self._block_grads(model, blk, groups, device)
            got = False
            for b, g in grads.items():
                if g is not None:
                    by_region[b][blk["region"]].append(g)
                    got = True
            used += int(got)
            if got:
                cells_used += int(blk.get("n", 0))
        self._zero(model)
        model.train(was_training)

        out: Dict[str, dict] = {}
        region_pair = _preferred_region_pair([blk["region"] for blk in blocks])
        for b, reg_map in by_region.items():
            cons = {r: self._consensus(gs) for r, gs in reg_map.items() if len(gs) >= 2}
            snr = float(np.mean(list(cons.values()))) if cons else float("nan")
            rcos = float("nan")
            if region_pair is not None and all(r in reg_map and len(reg_map[r]) >= 1 for r in region_pair):
                means = {r: torch.stack(gs).mean(dim=0) for r, gs in reg_map.items() if len(gs) >= 1}
                m1, m2 = means[region_pair[0]], means[region_pair[1]]
                den = float(torch.linalg.vector_norm(m1) * torch.linalg.vector_norm(m2))
                if den > 0:
                    rcos = float(torch.dot(m1, m2) / den)
            out[b] = {"snr": snr, "region_cos": rcos, "per_region": cons}

        state_cos = [
            v["region_cos"]
            for b, v in out.items()
            if b in STATE_BRANCHES and np.isfinite(v["region_cos"])
        ]
        state_snr = [
            v["snr"]
            for b, v in out.items()
            if b in STATE_BRANCHES and np.isfinite(v["snr"])
        ]
        min_state_cos = min(state_cos) if state_cos else float("nan")
        mean_state_snr = float(np.mean(state_snr)) if state_snr else float("nan")

        if np.isfinite(min_state_cos) and min_state_cos < 0.0:
            leakage_verdict = "HIGH: state branches disagree across regions"
            action = "add/check region fixed effect in base or celltype_region scaler, then rerun audit"
        elif np.isfinite(min_state_cos) and min_state_cos < 0.30:
            leakage_verdict = "WARN: possible region leakage into state"
            action = "inspect z_perp region separation before changing architecture"
        else:
            leakage_verdict = "LOW: no obvious region leakage into state"
            action = "keep current design; continue disease/module validation"

        if np.isfinite(mean_state_snr) and mean_state_snr < 0.25:
            snr_verdict = "weak donor-consistent state signal"
        elif np.isfinite(mean_state_snr) and mean_state_snr < 0.50:
            snr_verdict = "moderate donor-consistent state signal"
        else:
            snr_verdict = "strong donor-consistent state signal"

        pair_msg = f"{region_pair[0]} vs {region_pair[1]}" if region_pair else "region-pair unavailable"
        lines = [
            f"[diagnostics | donor-region SNR | epoch {epoch}]",
            f"scope: stateless ordinal NLL | mode={'train' if self.train_mode else 'eval'} | "
            f"blocks={used}/{len(blocks)} | cells={cells_used:,} | region_pair={pair_msg}",
            "",
            "QUESTION",
            "  Do donor blocks push each branch in the same direction?",
            "  Do DLPFC and MTG pull z/state in different directions?",
            "",
            "BRANCH SNR",
        ]
        for b, v in sorted(out.items()):
            snr = v["snr"]
            if not np.isfinite(snr):
                continue
            per = " ".join(
                f"{r}={val:.3f}" for r, val in sorted(v.get("per_region", {}).items())
            )
            lines.append(
                f"  {_branch_label(b):<16} {snr:>6.3f}  {_snr_status(snr):<26} {per}"
            )

        lines.extend(["", "REGION LEAKAGE"])
        for b, v in sorted(out.items()):
            rcos = v["region_cos"]
            if not np.isfinite(rcos):
                continue
            lines.append(
                f"  {_branch_label(b):<16} {rcos:>+6.2f}  {_region_cos_status(b, rcos)}"
            )

        lines.extend(
            [
                "",
                "VERDICT",
                f"  z/state region leakage : {leakage_verdict}",
                f"  donor-SNR state signal : {snr_verdict}",
                "",
                "ACTION",
                f"  {action}",
                "  Read with disease_program/module eval; high SNR alone can still be shared artifact.",
            ]
        )
        self._log("\n".join(lines))
        return out
