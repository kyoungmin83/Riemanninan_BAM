#!/usr/bin/env python
"""Read-only audit of the in-batch latent-kNN technical-zero estimator.

The audit deliberately opens only the training split.  It recreates grouped
per-rank batches, measures unrestricted-neighbor purity, evaluates the existing
binomial-thinning positive control, and recomputes two counterfactual scores
without changing model weights:

* same-cell-type neighbors in the current batch; and
* a same-cell-type, donor-balanced train-only latent bank.

The checkpoint, optimizer, scheduler, training process, and data splits are
never mutated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import fields
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset


SEX_LINKED_SENTINELS = (
    "XIST",
    "RPS4Y1",
    "DDX3Y",
    "KDM5D",
    "UTY",
    "EIF1AY",
    "TMSB4Y",
    "ZFY",
    "USP9Y",
    "NLGN4Y",
    "KDM6A",
    "ZFX",
    "RPS4X",
)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if torch.is_tensor(value):
        return value.detach().cpu().tolist()
    if isinstance(value, Path):
        return str(value)
    return value


class HistogramStats:
    def __init__(self, cap: float, bins: int = 1401) -> None:
        self.cap = float(cap)
        self.bins = int(bins)
        self.hist = np.zeros(self.bins, dtype=np.int64)
        self.count = 0
        self.sum = 0.0
        self.sumsq = 0.0
        self.maximum = 0.0
        self.positive = 0
        self.tail = {0.05: 0, 0.10: 0, 0.30: 0, 0.50: 0, 0.69: 0}

    def update(self, values: torch.Tensor) -> None:
        values = values.detach().float().reshape(-1)
        if values.numel() == 0:
            return
        self.count += int(values.numel())
        self.sum += float(values.sum().item())
        self.sumsq += float((values * values).sum().item())
        self.maximum = max(self.maximum, float(values.max().item()))
        self.positive += int((values > 0.0).sum().item())
        for threshold in self.tail:
            self.tail[threshold] += int((values >= threshold).sum().item())
        hist = torch.histc(values, bins=self.bins, min=0.0, max=self.cap)
        self.hist += hist.round().to(torch.int64).cpu().numpy()

    def quantile(self, probability: float) -> float:
        if self.count <= 0:
            return float("nan")
        target = max(1, int(math.ceil(float(probability) * self.count)))
        index = int(np.searchsorted(np.cumsum(self.hist), target, side="left"))
        index = min(max(index, 0), self.bins - 1)
        return float(index) * self.cap / float(self.bins - 1)

    def summary(self) -> dict[str, Any]:
        if self.count <= 0:
            return {"n": 0}
        mean = self.sum / self.count
        variance = max(0.0, self.sumsq / self.count - mean * mean)
        return {
            "n": self.count,
            "mean": mean,
            "sd": math.sqrt(variance),
            "min": 0.0,
            "max": self.maximum,
            "positive_fraction": self.positive / self.count,
            "q50": self.quantile(0.50),
            "q90": self.quantile(0.90),
            "q95": self.quantile(0.95),
            "q99": self.quantile(0.99),
            "tail_fraction": {
                str(threshold): count / self.count
                for threshold, count in self.tail.items()
            },
        }


class BinaryHistogram:
    def __init__(self, cap: float, bins: int = 1401) -> None:
        self.positive = HistogramStats(cap, bins)
        self.negative = HistogramStats(cap, bins)

    def update(
        self,
        scores: torch.Tensor,
        positive_mask: torch.Tensor,
        negative_mask: torch.Tensor,
    ) -> None:
        self.positive.update(scores[positive_mask])
        self.negative.update(scores[negative_mask])

    def summary(self) -> dict[str, Any]:
        pos = self.positive.summary()
        neg = self.negative.summary()
        result: dict[str, Any] = {"proven_dropout": pos, "stable_observed_zero": neg}
        if self.positive.count <= 0 or self.negative.count <= 0:
            return result
        pos_hist = self.positive.hist.astype(np.float64)
        neg_hist = self.negative.hist.astype(np.float64)
        tp = np.cumsum(pos_hist[::-1])
        fp = np.cumsum(neg_hist[::-1])
        recall = tp / max(float(tp[-1]), 1.0)
        fpr = fp / max(float(fp[-1]), 1.0)
        precision = tp / np.maximum(tp + fp, 1.0)
        recall0 = np.concatenate(([0.0], recall))
        fpr0 = np.concatenate(([0.0], fpr))
        tpr0 = recall0
        precision0 = np.concatenate(([1.0], precision))
        auroc = float(np.trapezoid(tpr0, fpr0))
        average_precision = float(np.sum(np.diff(recall0) * precision0[1:]))
        pos_mean = self.positive.sum / self.positive.count
        neg_mean = self.negative.sum / self.negative.count
        result.update(
            {
                "mean_gap": pos_mean - neg_mean,
                "mean_ratio": pos_mean / neg_mean if neg_mean > 0.0 else float("inf"),
                "auroc_histogram_approx": auroc,
                "average_precision_histogram_approx": average_precision,
                "positive_prevalence_in_diagnostic_pool": self.positive.count
                / (self.positive.count + self.negative.count),
            }
        )
        return result


class GroupBinaryAccumulator:
    def __init__(self) -> None:
        self.values: dict[tuple[str, int, str], dict[str, float]] = defaultdict(
            lambda: {"n_proven": 0.0, "sum_proven": 0.0, "n_stable": 0.0, "sum_stable": 0.0}
        )

    def update(
        self,
        stratifier: str,
        group: int,
        variant: str,
        scores: torch.Tensor,
        positive_mask: torch.Tensor,
        negative_mask: torch.Tensor,
    ) -> None:
        bucket = self.values[(str(stratifier), int(group), str(variant))]
        pos = scores[positive_mask]
        neg = scores[negative_mask]
        bucket["n_proven"] += float(pos.numel())
        bucket["sum_proven"] += float(pos.sum().item()) if pos.numel() else 0.0
        bucket["n_stable"] += float(neg.numel())
        bucket["sum_stable"] += float(neg.sum().item()) if neg.numel() else 0.0


class GroupMainAccumulator:
    def __init__(self) -> None:
        self.values: dict[tuple[str, int, str], dict[str, float]] = defaultdict(
            lambda: {"n_zero": 0.0, "sum_pi": 0.0, "n_positive": 0.0, "n_ge_0p3": 0.0}
        )

    def update(
        self,
        stratifier: str,
        group: int,
        variant: str,
        scores: torch.Tensor,
        zero_mask: torch.Tensor,
    ) -> None:
        bucket = self.values[(str(stratifier), int(group), str(variant))]
        values = scores[zero_mask]
        bucket["n_zero"] += float(values.numel())
        bucket["sum_pi"] += float(values.sum().item()) if values.numel() else 0.0
        bucket["n_positive"] += float((values > 0.0).sum().item()) if values.numel() else 0.0
        bucket["n_ge_0p3"] += float((values >= 0.3).sum().item()) if values.numel() else 0.0


def rank_quantile(x: torch.Tensor) -> torch.Tensor:
    n = int(x.shape[0])
    if n <= 1:
        return torch.zeros_like(x)
    order = torch.argsort(x)
    ranks = torch.empty_like(order)
    ranks[order] = torch.arange(n, device=x.device)
    return ranks.to(dtype=x.dtype) / float(n - 1)


def pi_from_neighbor_rate(
    y_ord: torch.Tensor,
    depth: torch.Tensor,
    gene_detect_rate: torch.Tensor,
    neighbor_rate: torch.Tensor,
    config: Any,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    depth_q = rank_quantile(depth.float())
    low_depth = torch.sigmoid(
        (float(config.depth_low_pct) - depth_q)
        / max(float(config.depth_low_temp), 1.0e-6)
    )
    span_g = max(
        float(config.neighbor_on_hi) - float(config.detect_floor),
        float(config.eps),
    )
    detectable = (
        (gene_detect_rate.float() - float(config.detect_floor)) / span_g
    ).clamp(0.0, 1.0)
    measure = low_depth[:, None] * detectable[None, :]
    span_n = max(
        float(config.neighbor_on_hi) - float(config.safety_neighbor_off),
        float(config.eps),
    )
    neighbor = (
        (neighbor_rate.float() - float(config.safety_neighbor_off)) / span_n
    ).clamp(0.0, 1.0)
    raw = float(config.w_measure) * measure + float(config.w_neighbor) * neighbor
    safety_pass = neighbor_rate > float(config.safety_neighbor_off)
    pi = torch.where(safety_pass, raw, torch.zeros_like(raw))
    pi = pi.clamp(0.0, float(config.pi_cap))
    pi = pi * (y_ord == 0).float()
    return pi, {
        "depth_quantile": depth_q,
        "low_depth": low_depth,
        "measure": measure,
        "neighbor": neighbor,
        "neighbor_rate": neighbor_rate,
        "safety_pass": safety_pass,
    }


def unrestricted_neighbors(
    z: torch.Tensor, detect: torch.Tensor, k: int
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    batch = int(z.shape[0])
    use_k = max(1, min(int(k), batch - 1))
    distance = torch.cdist(z.float(), z.float())
    distance = distance + torch.eye(batch, device=z.device) * 1.0e9
    neighbors = torch.topk(distance, use_k, dim=1, largest=False).indices
    rates = detect[neighbors].float().mean(dim=1)
    selected_distance = distance.gather(1, neighbors)
    return rates, neighbors, selected_distance


def same_celltype_inbatch_neighbors(
    z: torch.Tensor,
    detect: torch.Tensor,
    celltype: torch.Tensor,
    k: int,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    batch = int(z.shape[0])
    distance = torch.cdist(z.float(), z.float())
    allowed = celltype[:, None].eq(celltype[None, :])
    allowed.fill_diagonal_(False)
    distance = distance.masked_fill(~allowed, float("inf"))
    use_k = max(1, min(int(k), batch - 1))
    neighbors = torch.topk(distance, use_k, dim=1, largest=False).indices
    selected_distance = distance.gather(1, neighbors)
    valid = torch.isfinite(selected_distance)
    gathered = detect[neighbors].float() * valid[:, :, None].float()
    denominator = valid.sum(dim=1).clamp_min(1).float()
    rates = gathered.sum(dim=1) / denominator[:, None]
    return rates, neighbors, valid


def donor_balanced_bank_neighbors(
    z: torch.Tensor,
    query_row: torch.Tensor,
    query_celltype: torch.Tensor,
    bank_z: torch.Tensor,
    bank_detect: torch.Tensor,
    bank_row: torch.Tensor,
    bank_celltype: torch.Tensor,
    bank_donor: torch.Tensor,
    k: int,
    candidate_multiplier: int = 16,
) -> tuple[torch.Tensor, torch.Tensor]:
    distance = torch.cdist(z.float(), bank_z.float())
    allowed = query_celltype[:, None].eq(bank_celltype[None, :])
    allowed &= query_row[:, None].ne(bank_row[None, :])
    distance = distance.masked_fill(~allowed, float("inf"))
    candidate_k = min(
        int(bank_z.shape[0]), max(int(k), int(k) * int(candidate_multiplier))
    )
    candidates = torch.topk(distance, candidate_k, dim=1, largest=False).indices
    candidates_cpu = candidates.cpu().numpy()
    distances_cpu = distance.gather(1, candidates).cpu().numpy()
    donors_cpu = bank_donor.cpu().numpy()
    selected: list[list[int]] = []
    for row_candidates, row_distances in zip(candidates_cpu, distances_cpu):
        chosen: list[int] = []
        seen_donors: set[int] = set()
        for index, value in zip(row_candidates.tolist(), row_distances.tolist()):
            if not math.isfinite(float(value)):
                break
            donor = int(donors_cpu[index])
            if donor in seen_donors:
                continue
            chosen.append(int(index))
            seen_donors.add(donor)
            if len(chosen) >= int(k):
                break
        if len(chosen) < int(k):
            raise RuntimeError(
                f"same-celltype bank has only {len(chosen)} distinct-donor candidates; need {k}"
            )
        selected.append(chosen)
    neighbor_index = torch.as_tensor(selected, dtype=torch.long, device=z.device)
    rates = bank_detect[neighbor_index].float().mean(dim=1)
    return rates, neighbor_index


def slice_batch(batch: dict[str, Any], start: int, stop: int, total: int, device: torch.device) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in batch.items():
        if torch.is_tensor(value):
            if value.ndim > 0 and int(value.shape[0]) == total:
                value = value[start:stop]
            result[key] = value.to(device, non_blocking=False)
        elif isinstance(value, list) and len(value) == total:
            result[key] = value[start:stop]
        else:
            result[key] = value
    return result


def forward_z(
    system: torch.nn.Module,
    batch: dict[str, Any],
    device: torch.device,
    chunk_size: int,
) -> torch.Tensor:
    total = int(batch["y_ord"].shape[0])
    chunks = []
    for start in range(0, total, int(chunk_size)):
        stop = min(total, start + int(chunk_size))
        device_batch = slice_batch(batch, start, stop, total, device)
        with torch.autocast(
            device_type=device.type,
            dtype=torch.bfloat16,
            enabled=device.type == "cuda",
        ):
            output = system(
                device_batch,
                sample_latent=False,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )
        chunks.append(output.z_perp.detach().float())
        del output, device_batch
    return torch.cat(chunks, dim=0)


def replace_with_thinned(batch: dict[str, Any]) -> dict[str, Any]:
    result = dict(batch)
    result["y_ord"] = batch["y_ord_thin"]
    if "x_log1p_thin" in batch:
        result["x_log1p"] = batch["x_log1p_thin"]
    if "x_gene_scalar_thin" in batch:
        result["x_gene_scalar"] = batch["x_gene_scalar_thin"]
    return result


def select_bank_local_indices(dataset: Any, per_group: int, seed: int) -> np.ndarray:
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for local_index, row in enumerate(np.asarray(dataset.row_idx, dtype=np.int64)):
        groups[(int(dataset.celltype_ids[row]), int(dataset.donor_ids[row]))].append(local_index)
    rng = np.random.default_rng(int(seed))
    selected: list[np.ndarray] = []
    for key in sorted(groups):
        values = np.asarray(groups[key], dtype=np.int64)
        take = min(len(values), int(per_group))
        selected.append(rng.choice(values, size=take, replace=False))
    return np.sort(np.concatenate(selected)).astype(np.int64)


def update_group_diagnostics(
    binary_groups: GroupBinaryAccumulator,
    main_groups: GroupMainAccumulator,
    variants_full: dict[str, torch.Tensor],
    variants_thin: dict[str, torch.Tensor],
    y_full: torch.Tensor,
    y_thin: torch.Tensor,
    labels: dict[str, torch.Tensor],
    gene_detect_rate: torch.Tensor,
) -> None:
    proven = (y_full > 0) & (y_thin == 0)
    stable = (y_full == 0) & (y_thin == 0)
    zero = y_full == 0
    for stratifier, cell_labels in labels.items():
        for group in torch.unique(cell_labels).tolist():
            cell_mask = cell_labels == int(group)
            pair_mask = cell_mask[:, None]
            for variant, scores in variants_thin.items():
                binary_groups.update(
                    stratifier,
                    int(group),
                    variant,
                    scores,
                    proven & pair_mask,
                    stable & pair_mask,
                )
            for variant, scores in variants_full.items():
                main_groups.update(
                    stratifier,
                    int(group),
                    variant,
                    scores,
                    zero & pair_mask,
                )

    for source_bin in (1, 2, 3):
        source_mask = y_full == source_bin
        for variant, scores in variants_thin.items():
            binary_groups.update(
                "source_ordinal_bin",
                source_bin,
                variant,
                scores,
                proven & source_mask,
                stable & source_mask,
            )

    detect_bands = (
        gene_detect_rate < 0.05,
        (gene_detect_rate >= 0.05) & (gene_detect_rate < 0.20),
        (gene_detect_rate >= 0.20) & (gene_detect_rate < 0.50),
        gene_detect_rate >= 0.50,
    )
    for group, gene_mask in enumerate(detect_bands):
        pair_mask = gene_mask[None, :]
        for variant, scores in variants_thin.items():
            binary_groups.update(
                "gene_detect_rate_band",
                group,
                variant,
                scores,
                proven & pair_mask,
                stable & pair_mask,
            )


def group_name(stratifier: str, group: int, dataset: Any) -> str:
    if stratifier == "celltype" and 0 <= group < len(dataset.spec.celltype_vocab):
        return str(dataset.spec.celltype_vocab[group])
    if stratifier == "sex":
        return {-1: "unknown", 0: "female", 1: "male"}.get(group, str(group))
    if stratifier == "tech" and 0 <= group < len(dataset.spec.batch_vocab):
        return str(dataset.spec.batch_vocab[group])
    if stratifier == "source_ordinal_bin":
        return f"ordinal_{group}"
    if stratifier == "gene_detect_rate_band":
        return ("<0.05", "0.05-0.20", "0.20-0.50", ">=0.50")[group]
    return str(group)


def write_group_tables(
    outdir: Path,
    binary: GroupBinaryAccumulator,
    main: GroupMainAccumulator,
    dataset: Any,
) -> None:
    with (outdir / "thinning_stratified.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "stratifier",
                "group_id",
                "group_name",
                "variant",
                "n_proven",
                "n_stable",
                "pi_proven_mean",
                "pi_stable_mean",
                "mean_gap",
                "mean_ratio",
            )
        )
        for (stratifier, group, variant), values in sorted(binary.values.items()):
            np_ = int(values["n_proven"])
            ns = int(values["n_stable"])
            mp = values["sum_proven"] / np_ if np_ else float("nan")
            ms = values["sum_stable"] / ns if ns else float("nan")
            writer.writerow(
                (
                    stratifier,
                    group,
                    group_name(stratifier, group, dataset),
                    variant,
                    np_,
                    ns,
                    mp,
                    ms,
                    mp - ms,
                    mp / ms if ms > 0.0 else float("inf"),
                )
            )

    with (outdir / "main_zero_stratified.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "stratifier",
                "group_id",
                "group_name",
                "variant",
                "n_observed_zero",
                "pi_mean",
                "pi_positive_fraction",
                "pi_ge_0p3_fraction",
            )
        )
        for (stratifier, group, variant), values in sorted(main.values.items()):
            count = int(values["n_zero"])
            writer.writerow(
                (
                    stratifier,
                    group,
                    group_name(stratifier, group, dataset),
                    variant,
                    count,
                    values["sum_pi"] / count if count else float("nan"),
                    values["n_positive"] / count if count else float("nan"),
                    values["n_ge_0p3"] / count if count else float("nan"),
                )
            )


def update_sex_gene_stats(
    accumulator: dict[tuple[str, int, str], dict[str, float]],
    gene_indices: dict[str, int],
    variants: dict[str, torch.Tensor],
    y_ord: torch.Tensor,
    sex: torch.Tensor,
    neighbor_rate: torch.Tensor,
) -> None:
    for gene, index in gene_indices.items():
        for sex_id in torch.unique(sex).tolist():
            cell_mask = sex == int(sex_id)
            zero_mask = cell_mask & (y_ord[:, index] == 0)
            for variant, scores in variants.items():
                values = scores[:, index][zero_mask]
                bucket = accumulator[(gene, int(sex_id), variant)]
                bucket["n_zero"] += float(values.numel())
                bucket["sum_pi"] += float(values.sum().item()) if values.numel() else 0.0
                bucket["n_positive"] += float((values > 0).sum().item()) if values.numel() else 0.0
                if variant == "current_inbatch":
                    nb = neighbor_rate[:, index][zero_mask]
                    bucket["sum_neighbor_rate"] += float(nb.sum().item()) if nb.numel() else 0.0


def write_sex_gene_table(
    outdir: Path,
    accumulator: dict[tuple[str, int, str], dict[str, float]],
) -> None:
    with (outdir / "sex_linked_gene_zero_metrics.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "gene",
                "sex_id",
                "sex_name",
                "variant",
                "n_observed_zero",
                "pi_mean",
                "pi_positive_fraction",
                "current_neighbor_detection_mean",
            )
        )
        for (gene, sex_id, variant), values in sorted(accumulator.items()):
            count = int(values.get("n_zero", 0.0))
            writer.writerow(
                (
                    gene,
                    sex_id,
                    {-1: "unknown", 0: "female", 1: "male"}.get(sex_id, str(sex_id)),
                    variant,
                    count,
                    values.get("sum_pi", 0.0) / count if count else float("nan"),
                    values.get("n_positive", 0.0) / count if count else float("nan"),
                    (
                        values.get("sum_neighbor_rate", 0.0) / count
                        if count and variant == "current_inbatch"
                        else float("nan")
                    ),
                )
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-config", required=True)
    parser.add_argument("--resolved-config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--support-scripts", required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--epoch", type=int, required=True)
    parser.add_argument("--world-size", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=48)
    parser.add_argument("--block-size", type=int, default=8)
    parser.add_argument("--batches-per-rank", type=int, default=8)
    parser.add_argument("--bank-per-donor-celltype", type=int, default=2)
    parser.add_argument("--bank-loader-batch-size", type=int, default=24)
    parser.add_argument("--model-chunk-size", type=int, default=8)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260827)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(Path(args.support_scripts).resolve()))
    from kmlee_bam.data.grouped_batch_sampler import GroupedBatchSampler
    from kmlee_bam.data.ordinal_dataset import read_zarr_dataframe
    from kmlee_bam.objectives.latent_knn_pi_tech import LatentKNNPiTechConfig
    from kmlee_bam.training import run_current as current
    from kmlee_bam.training import runner_base as base
    from kmlee_bam.training.core_trainer import _resolve_tech_id
    from prism_posthoc_system import (
        build_posthoc_system,
        validate_posthoc_checkpoint_compatibility,
    )

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    source = json.load(open(args.source_config, encoding="utf-8"))
    resolved = json.load(open(args.resolved_config, encoding="utf-8"))
    if source.get("train", {}).get("skip_final_test") is not True:
        raise RuntimeError("audit requires skip_final_test=true")
    if source.get("experiment_manifest", {}).get("official_test_donors") != "sealed_not_opened":
        raise RuntimeError("official test donor seal is not declared")
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    checkpoint_epoch = int(payload.get("epoch", -1))
    if checkpoint_epoch != int(args.epoch):
        raise RuntimeError(f"checkpoint epoch {checkpoint_epoch} != requested epoch {args.epoch}")

    current.set_settings(source)
    current.install_hooks()
    cfg = base.load_config(args.resolved_config)
    ds_cfg = current._dataset_section()
    common = dict(
        zarr_path=cfg.data.zarr_path,
        spec_path=cfg.data.spec_path,
        matrix_path=cfg.data.matrix_path,
        counts_layer=cfg.data.counts_layer,
        prefer_raw=cfg.data.prefer_raw,
        split_key=cfg.data.split_key,
        celltype_key=cfg.data.celltype_key,
        batch_key=cfg.data.batch_key,
        return_log_input=cfg.data.return_log_input,
        reference_scaler_path=cfg.data.reference_scaler_path,
        return_x_gene_scalar=cfg.data.return_x_gene_scalar,
        x_gene_scalar_clip=cfg.data.x_gene_scalar_clip,
        x_gene_scalar_eps=cfg.data.x_gene_scalar_eps,
        max_csr_data_span=cfg.data.max_csr_data_span,
        max_span_to_selected_nnz_ratio=cfg.data.max_span_to_selected_nnz_ratio,
        donor_key=ds_cfg.get("donor_key", "donor_id"),
        depth_key=ds_cfg.get("depth_key", "Number of UMIs"),
        detected_genes_key=ds_cfg.get("detected_genes_key", "Genes detected"),
    )
    thin_cfg = source["thinning"]
    train_query = current.V3OrdinalScDataset(
        split="train",
        return_thinned=True,
        thin_rho=float(thin_cfg["rho"]),
        thin_counts_zarr_path=thin_cfg["counts_zarr_path"],
        thin_counts_matrix_path=str(thin_cfg.get("counts_matrix_path", "X")),
        thin_library_size_key=str(thin_cfg.get("library_size_key", "Number of UMIs")),
        **common,
    )
    train_bank = current.V3OrdinalScDataset(split="train", **common)
    obs = read_zarr_dataframe(train_query.root, "obs")

    print(
        f"[pi-audit] train-only cells={len(train_query):,}; test opened=false; checkpoint=e{checkpoint_epoch}",
        flush=True,
    )
    device = torch.device(args.device)
    system = build_posthoc_system(cfg, train_query, obs)
    state = payload.get("system_state_dict", payload)
    missing, unexpected = system.load_state_dict(state, strict=False)
    compatibility = validate_posthoc_checkpoint_compatibility(missing, unexpected)
    print(f"[pi-audit] checkpoint compatibility={compatibility}", flush=True)
    system.to(device).eval()
    for parameter in system.parameters():
        parameter.requires_grad_(False)

    config_keys = {field.name for field in fields(LatentKNNPiTechConfig)}
    pi_config = LatentKNNPiTechConfig(
        **{key: value for key, value in source["pi_tech"].items() if key in config_keys}
    )
    gene_detect_rate = payload["v8_state"]["gene_detect_rate"].detach().float().to(device)
    if int(gene_detect_rate.numel()) != len(train_query.spec.gene_names):
        raise RuntimeError("checkpoint gene_detect_rate dimension mismatch")

    bank_local = select_bank_local_indices(
        train_bank, args.bank_per_donor_celltype, args.seed + 17
    )
    bank_loader = DataLoader(
        Subset(train_bank, bank_local.tolist()),
        batch_size=int(args.bank_loader_batch_size),
        shuffle=False,
        num_workers=int(args.num_workers),
        pin_memory=False,
        persistent_workers=False,
    )
    bank_chunks: dict[str, list[torch.Tensor]] = defaultdict(list)
    with torch.inference_mode():
        for batch_index, batch in enumerate(bank_loader):
            bank_chunks["z"].append(
                forward_z(system, batch, device, args.model_chunk_size).cpu()
            )
            bank_chunks["detect"].append((batch["y_ord"] > 0).to(torch.uint8))
            bank_chunks["row"].append(batch["row_index"].long())
            bank_chunks["celltype"].append(batch["celltype_id"].long())
            bank_chunks["donor"].append(batch["donor_id"].long())
            bank_chunks["sex"].append(batch["sex_id"].long())
            bank_chunks["tech"].append(_resolve_tech_id(batch).long())
            if (batch_index + 1) % 50 == 0 or batch_index + 1 == len(bank_loader):
                print(f"[pi-audit] bank {batch_index + 1}/{len(bank_loader)}", flush=True)
    bank_z = torch.cat(bank_chunks["z"]).to(device)
    bank_detect = torch.cat(bank_chunks["detect"]).to(device)
    bank_row = torch.cat(bank_chunks["row"]).to(device)
    bank_celltype = torch.cat(bank_chunks["celltype"]).to(device)
    bank_donor = torch.cat(bank_chunks["donor"]).to(device)
    bank_sex = torch.cat(bank_chunks["sex"]).to(device)
    bank_tech = torch.cat(bank_chunks["tech"]).to(device)
    print(
        f"[pi-audit] bank cells={len(bank_row):,}; donor-balanced query k={pi_config.k}",
        flush=True,
    )

    variants = ("current_inbatch", "same_celltype_inbatch", "same_celltype_donor_bank")
    main_hist = {name: HistogramStats(pi_config.pi_cap) for name in variants}
    thinning_hist = {name: BinaryHistogram(pi_config.pi_cap) for name in variants}
    component_hist = {
        name: HistogramStats(1.0)
        for name in ("measure", "neighbor", "neighbor_rate")
    }
    diff_hist = {
        "samect_minus_current_main": HistogramStats(1.4),
        "bank_minus_current_main": HistogramStats(1.4),
    }
    binary_groups = GroupBinaryAccumulator()
    main_groups = GroupMainAccumulator()
    neighbor_totals = defaultdict(float)
    neighbor_rank = {
        key: np.zeros(int(pi_config.k), dtype=np.float64)
        for key in ("same_celltype", "same_donor", "same_block", "same_sex", "same_tech")
    }
    neighbor_by_celltype: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    n_celltypes = len(train_query.spec.celltype_vocab)
    cross_celltype = np.zeros((n_celltypes, n_celltypes), dtype=np.int64)
    sex_gene_indices = {
        str(gene): int(index)
        for index, gene in enumerate(np.asarray(train_query.spec.gene_names, dtype=object))
        if str(gene).upper() in SEX_LINKED_SENTINELS
    }
    sex_gene_stats: dict[tuple[str, int, str], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )

    query_cells = 0
    query_batches = 0
    all_pairs = 0
    main_zero_count = 0
    main_weight_sum = 0.0
    main_weight_min = 1.0
    safety_pass_zero = 0
    counterfactual_abs = defaultdict(float)
    counterfactual_signed = defaultdict(float)
    counterfactual_count = 0
    samect_effective_k = np.zeros(int(pi_config.k) + 1, dtype=np.int64)

    with torch.inference_mode():
        for rank in range(int(args.world_size)):
            sampler = GroupedBatchSampler(
                train_query,
                batch_size=int(args.batch_size),
                block_size=int(args.block_size),
                num_replicas=int(args.world_size),
                rank=rank,
                seed=int(source["train"]["seed"]),
                drop_last=False,
            )
            sampler.set_epoch(int(args.epoch))
            loader = DataLoader(
                train_query,
                batch_sampler=sampler,
                num_workers=int(args.num_workers),
                pin_memory=False,
                persistent_workers=False,
            )
            for local_batch_index, cpu_batch in enumerate(loader):
                if local_batch_index >= int(args.batches_per_rank):
                    break
                batch_size = int(cpu_batch["y_ord"].shape[0])
                z_full = forward_z(system, cpu_batch, device, args.model_chunk_size)
                thin_batch = replace_with_thinned(cpu_batch)
                z_thin = forward_z(system, thin_batch, device, args.model_chunk_size)
                y_full = cpu_batch["y_ord"].long().to(device)
                y_thin = cpu_batch["y_ord_thin"].long().to(device)
                detect_full = y_full > 0
                detect_thin = y_thin > 0
                depth_full = detect_full.float().sum(dim=1)
                depth_thin = detect_thin.float().sum(dim=1)
                row = cpu_batch["row_index"].long().to(device)
                celltype = cpu_batch["celltype_id"].long().to(device)
                donor = cpu_batch["donor_id"].long().to(device)
                sex = cpu_batch["sex_id"].long().to(device)
                tech = _resolve_tech_id(cpu_batch).long().to(device)

                nb_full_current, neighbor_index, _ = unrestricted_neighbors(
                    z_full, detect_full, pi_config.k
                )
                nb_thin_current, _, _ = unrestricted_neighbors(
                    z_thin, detect_thin, pi_config.k
                )
                nb_full_samect, _, samect_valid = same_celltype_inbatch_neighbors(
                    z_full, detect_full, celltype, pi_config.k
                )
                nb_thin_samect, _, _ = same_celltype_inbatch_neighbors(
                    z_thin, detect_thin, celltype, pi_config.k
                )
                nb_full_bank, _ = donor_balanced_bank_neighbors(
                    z_full,
                    row,
                    celltype,
                    bank_z,
                    bank_detect,
                    bank_row,
                    bank_celltype,
                    bank_donor,
                    pi_config.k,
                )
                nb_thin_bank, _ = donor_balanced_bank_neighbors(
                    z_thin,
                    row,
                    celltype,
                    bank_z,
                    bank_detect,
                    bank_row,
                    bank_celltype,
                    bank_donor,
                    pi_config.k,
                )

                pi_full_current, current_details = pi_from_neighbor_rate(
                    y_full, depth_full, gene_detect_rate, nb_full_current, pi_config
                )
                pi_full_samect, _ = pi_from_neighbor_rate(
                    y_full, depth_full, gene_detect_rate, nb_full_samect, pi_config
                )
                pi_full_bank, _ = pi_from_neighbor_rate(
                    y_full, depth_full, gene_detect_rate, nb_full_bank, pi_config
                )
                pi_thin_current, _ = pi_from_neighbor_rate(
                    y_thin, depth_thin, gene_detect_rate, nb_thin_current, pi_config
                )
                pi_thin_samect, _ = pi_from_neighbor_rate(
                    y_thin, depth_thin, gene_detect_rate, nb_thin_samect, pi_config
                )
                pi_thin_bank, _ = pi_from_neighbor_rate(
                    y_thin, depth_thin, gene_detect_rate, nb_thin_bank, pi_config
                )
                variants_full = {
                    "current_inbatch": pi_full_current,
                    "same_celltype_inbatch": pi_full_samect,
                    "same_celltype_donor_bank": pi_full_bank,
                }
                variants_thin = {
                    "current_inbatch": pi_thin_current,
                    "same_celltype_inbatch": pi_thin_samect,
                    "same_celltype_donor_bank": pi_thin_bank,
                }
                zero = y_full == 0
                proven = (y_full > 0) & (y_thin == 0)
                stable = (y_full == 0) & (y_thin == 0)
                for name in variants:
                    main_hist[name].update(variants_full[name][zero])
                    thinning_hist[name].update(variants_thin[name], proven, stable)

                component_hist["measure"].update(current_details["measure"][zero])
                component_hist["neighbor"].update(current_details["neighbor"][zero])
                component_hist["neighbor_rate"].update(current_details["neighbor_rate"][zero])
                safety_pass_zero += int(current_details["safety_pass"][zero].sum().item())
                difference_samect = pi_full_samect[zero] - pi_full_current[zero]
                difference_bank = pi_full_bank[zero] - pi_full_current[zero]
                # HistogramStats is non-negative; shift signed differences by the cap.
                diff_hist["samect_minus_current_main"].update(
                    difference_samect + float(pi_config.pi_cap)
                )
                diff_hist["bank_minus_current_main"].update(
                    difference_bank + float(pi_config.pi_cap)
                )
                counterfactual_abs["samect"] += float(difference_samect.abs().sum().item())
                counterfactual_abs["bank"] += float(difference_bank.abs().sum().item())
                counterfactual_signed["samect"] += float(difference_samect.sum().item())
                counterfactual_signed["bank"] += float(difference_bank.sum().item())
                counterfactual_count += int(difference_samect.numel())
                for effective_k in samect_valid.sum(dim=1).cpu().tolist():
                    samect_effective_k[int(effective_k)] += 1

                alpha = float(source["zero_reliability"]["alpha"])
                floor = float(source["zero_reliability"]["floor"])
                weights = (1.0 - alpha * pi_full_current).clamp(floor, 1.0)
                weights = torch.where(zero, weights, torch.ones_like(weights))
                main_weight_sum += float(weights.sum().item())
                main_weight_min = min(main_weight_min, float(weights.min().item()))
                all_pairs += int(weights.numel())
                main_zero_count += int(zero.sum().item())

                neighbor_ct = celltype[neighbor_index]
                neighbor_donor = donor[neighbor_index]
                neighbor_sex = sex[neighbor_index]
                neighbor_tech = tech[neighbor_index]
                same_ct = neighbor_ct.eq(celltype[:, None])
                same_donor = neighbor_donor.eq(donor[:, None])
                same_block = same_ct & same_donor
                same_sex = neighbor_sex.eq(sex[:, None])
                same_tech = neighbor_tech.eq(tech[:, None])
                attribute_masks = {
                    "same_celltype": same_ct,
                    "same_donor": same_donor,
                    "same_block": same_block,
                    "same_sex": same_sex,
                    "same_tech": same_tech,
                }
                for key, mask in attribute_masks.items():
                    neighbor_totals[key] += float(mask.sum().item())
                    neighbor_rank[key] += mask.float().sum(dim=0).cpu().numpy()
                neighbor_totals["edges"] += float(same_ct.numel())
                for query_type in torch.unique(celltype).tolist():
                    query_mask = celltype == int(query_type)
                    bucket = neighbor_by_celltype[int(query_type)]
                    edge_count = int(query_mask.sum().item()) * int(pi_config.k)
                    bucket["edges"] += float(edge_count)
                    bucket["queries"] += float(query_mask.sum().item())
                    for key, mask in attribute_masks.items():
                        bucket[key] += float(mask[query_mask].sum().item())
                for query_type, neighbor_type in zip(
                    celltype[:, None].expand_as(neighbor_ct).reshape(-1).cpu().tolist(),
                    neighbor_ct.reshape(-1).cpu().tolist(),
                ):
                    cross_celltype[int(query_type), int(neighbor_type)] += 1

                eye = torch.eye(batch_size, dtype=torch.bool, device=device)
                candidate = ~eye
                chance_masks = {
                    "chance_same_celltype": celltype[:, None].eq(celltype[None, :]) & candidate,
                    "chance_same_donor": donor[:, None].eq(donor[None, :]) & candidate,
                    "chance_same_block": (
                        celltype[:, None].eq(celltype[None, :])
                        & donor[:, None].eq(donor[None, :])
                        & candidate
                    ),
                    "chance_same_sex": sex[:, None].eq(sex[None, :]) & candidate,
                    "chance_same_tech": tech[:, None].eq(tech[None, :]) & candidate,
                }
                for key, mask in chance_masks.items():
                    neighbor_totals[key] += float(mask.sum().item()) / float(batch_size - 1)
                neighbor_totals["chance_queries"] += float(batch_size)

                labels = {"celltype": celltype, "sex": sex, "tech": tech}
                update_group_diagnostics(
                    binary_groups,
                    main_groups,
                    variants_full,
                    variants_thin,
                    y_full,
                    y_thin,
                    labels,
                    gene_detect_rate,
                )
                update_sex_gene_stats(
                    sex_gene_stats,
                    sex_gene_indices,
                    variants_full,
                    y_full,
                    sex,
                    nb_full_current,
                )

                query_cells += batch_size
                query_batches += 1
                print(
                    f"[pi-audit] query rank={rank} batch={local_batch_index + 1}/{args.batches_per_rank} "
                    f"cells={query_cells:,}",
                    flush=True,
                )
                del z_full, z_thin, y_full, y_thin

    edges = max(neighbor_totals["edges"], 1.0)
    chance_queries = max(neighbor_totals["chance_queries"], 1.0)
    neighbor_summary = {
        key: neighbor_totals[key] / edges
        for key in ("same_celltype", "same_donor", "same_block", "same_sex", "same_tech")
    }
    neighbor_summary["batch_chance"] = {
        key.replace("chance_", ""): neighbor_totals[key] / chance_queries
        for key in (
            "chance_same_celltype",
            "chance_same_donor",
            "chance_same_block",
            "chance_same_sex",
            "chance_same_tech",
        )
    }
    neighbor_summary["enrichment_over_batch_chance"] = {
        key: neighbor_summary[key] / neighbor_summary["batch_chance"][key]
        if neighbor_summary["batch_chance"][key] > 0.0
        else float("nan")
        for key in ("same_celltype", "same_donor", "same_block", "same_sex", "same_tech")
    }
    neighbor_summary["by_neighbor_rank"] = {
        key: (values / max(query_cells, 1)).tolist()
        for key, values in neighbor_rank.items()
    }

    with (outdir / "neighbor_purity_by_celltype.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "celltype_id",
                "celltype",
                "n_queries",
                "same_celltype",
                "same_donor",
                "same_block",
                "same_sex",
                "same_tech",
            )
        )
        for celltype_id, values in sorted(neighbor_by_celltype.items()):
            edge_count = max(values["edges"], 1.0)
            writer.writerow(
                (
                    celltype_id,
                    train_query.spec.celltype_vocab[celltype_id],
                    int(values["queries"]),
                    values["same_celltype"] / edge_count,
                    values["same_donor"] / edge_count,
                    values["same_block"] / edge_count,
                    values["same_sex"] / edge_count,
                    values["same_tech"] / edge_count,
                )
            )

    write_group_tables(outdir, binary_groups, main_groups, train_query)
    write_sex_gene_table(outdir, sex_gene_stats)
    np.savetxt(outdir / "neighbor_celltype_count_matrix.tsv", cross_celltype, fmt="%d", delimiter="\t")

    current_mean_pi_all_pairs = main_hist["current_inbatch"].sum / max(all_pairs, 1)
    summary = {
        "schema": "kmlee_bam.latent_knn_pi_tech_posthoc.v1",
        "scope": {
            "splits_opened": ["train"],
            "official_test_opened": False,
            "model_or_optimizer_mutated": False,
            "checkpoint_epoch": checkpoint_epoch,
            "query_batches": query_batches,
            "query_cells": query_cells,
            "world_size_recreated": int(args.world_size),
            "batches_per_rank": int(args.batches_per_rank),
            "batch_size": int(args.batch_size),
            "block_size": int(args.block_size),
            "bank_cells": int(bank_row.numel()),
            "bank_per_donor_celltype": int(args.bank_per_donor_celltype),
        },
        "pi_config": source["pi_tech"],
        "neighbor_purity": neighbor_summary,
        "main_observed_zero": {name: main_hist[name].summary() for name in variants},
        "thinning_diagnostic": {name: thinning_hist[name].summary() for name in variants},
        "current_components_at_observed_zero": {
            name: stats.summary() for name, stats in component_hist.items()
        },
        "current_safety_pass_fraction_at_observed_zero": safety_pass_zero / max(main_zero_count, 1),
        "current_reconstruction_weight": {
            "mean_over_all_gene_positions": main_weight_sum / max(all_pairs, 1),
            "minimum": main_weight_min,
            "mean_pi_times_zero_indicator_over_all_positions": current_mean_pi_all_pairs,
        },
        "counterfactual": {
            "mean_absolute_pi_change_at_observed_zero": {
                "same_celltype_inbatch_vs_current": counterfactual_abs["samect"]
                / max(counterfactual_count, 1),
                "same_celltype_donor_bank_vs_current": counterfactual_abs["bank"]
                / max(counterfactual_count, 1),
            },
            "mean_signed_pi_change_at_observed_zero": {
                "same_celltype_inbatch_minus_current": counterfactual_signed["samect"]
                / max(counterfactual_count, 1),
                "same_celltype_donor_bank_minus_current": counterfactual_signed["bank"]
                / max(counterfactual_count, 1),
            },
            "same_celltype_inbatch_effective_k_counts": {
                str(index): int(count)
                for index, count in enumerate(samect_effective_k)
                if count > 0
            },
        },
        "sex_linked_genes_present": sorted(sex_gene_indices),
        "checkpoint_compatibility": compatibility,
        "provenance": {
            "source_config": str(Path(args.source_config).resolve()),
            "source_config_sha256": sha256_file(args.source_config),
            "resolved_config": str(Path(args.resolved_config).resolve()),
            "resolved_config_sha256": sha256_file(args.resolved_config),
            "checkpoint": str(Path(args.checkpoint).resolve()),
            "checkpoint_sha256": sha256_file(args.checkpoint),
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(__file__),
        },
    }
    with (outdir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(json_ready(summary), handle, indent=2, sort_keys=True, allow_nan=True)
        handle.write("\n")
    (outdir / "VALIDATION_COMPLETE").touch()
    print(f"[pi-audit] COMPLETE outdir={outdir}", flush=True)


if __name__ == "__main__":
    main()
