#!/usr/bin/env python
"""Paired frozen-checkpoint PRISM path ablation on identical held-out cells.

This evaluator runs the encoder/decoder once, then subtracts one explicit PRISM
score component at a time.  It therefore measures conditional reliance of the
frozen model, not a retrained ablation and not causal importance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import defaultdict

import numpy as np
import torch
from torch.utils.data import DataLoader

from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from disease_program_eval import build_donor_pathology, build_pooled_dataset, move_batch
from prism_posthoc_system import build_posthoc_system, validate_posthoc_checkpoint_compatibility


VARIANTS = (
    "full",
    "common_pathology_off",
    "pathology_interaction_off",
    "all_shared_pathology_off",
    "personal_rank2_baseline_off",
    "personal_pathology_response_off",
    "module_local_off",
    "module_local_nonlinear_off",
    "all_personal_off",
)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--out-npz", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--eval-split", choices=("val", "test"), default="val")
    parser.add_argument("--per-donor-celltype", type=int, default=40)
    parser.add_argument("--sample-seed", type=int, default=20260814)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--min-module-genes", type=int, default=3)
    parser.add_argument("--max-batches", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.per_donor_celltype <= 0:
        raise ValueError("--per-donor-celltype must be positive")
    torch.manual_seed(42)
    np.random.seed(42)

    with open(args.config, encoding="utf-8") as handle:
        raw = json.load(handle)
    current.set_settings(raw)
    current.install_hooks()
    cfg = base.load_config(args.config)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    ds, obs, n_val, n_test = build_pooled_dataset(cfg)
    split_values = obs[cfg.data.split_key].astype(str).to_numpy()
    selected_rows = np.flatnonzero(split_values == args.eval_split).astype(np.int64)
    if selected_rows.size == 0:
        raise SystemExit(f"no rows for split={args.eval_split}")
    sample_rng = np.random.default_rng(args.sample_seed)
    grouped_rows: dict[tuple[int, int], list[int]] = defaultdict(list)
    for row in selected_rows:
        grouped_rows[(int(ds.donor_ids[row]), int(ds.celltype_ids[row]))].append(int(row))
    sampled_rows = []
    for key in sorted(grouped_rows):
        rows = np.asarray(grouped_rows[key], dtype=np.int64)
        if rows.size > args.per_donor_celltype:
            rows = sample_rng.choice(rows, size=args.per_donor_celltype, replace=False)
        sampled_rows.append(rows)
    selected_rows = np.sort(np.concatenate(sampled_rows)).astype(np.int64)
    ds.row_idx = selected_rows

    donor_ids = np.unique(ds.donor_ids[selected_rows]).astype(np.int64)
    donor_pos = np.full(ds.n_donor, -1, dtype=np.int64)
    donor_pos[donor_ids] = np.arange(donor_ids.size)
    donor_names = np.asarray([str(ds.donor_vocab[index]) for index in donor_ids], dtype=str)
    celltype_names = np.asarray(ds.spec.celltype_vocab, dtype=str)

    system = build_posthoc_system(cfg, ds, obs)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = payload.get("system_state_dict", payload) if isinstance(payload, dict) else payload
    missing, unexpected = system.load_state_dict(state, strict=False)
    compatibility = validate_posthoc_checkpoint_compatibility(missing, unexpected)
    system.to(device).eval()
    for parameter in system.parameters():
        parameter.requires_grad_(False)

    registry_path = cfg.module_tokenizer.registry_json_path
    with open(registry_path, encoding="utf-8") as handle:
        registry = json.load(handle)
    registry_genes = np.asarray(registry["gene_names"], dtype=str)
    model_genes = np.asarray(ds.spec.gene_names, dtype=str)
    if not np.array_equal(registry_genes, model_genes):
        raise RuntimeError("module registry gene order does not match dataset")
    membership = np.asarray(registry["membership_binary"], dtype=np.float32)
    keep = membership.sum(axis=1) >= args.min_module_genes
    module_names = np.asarray(registry["module_names"], dtype=str)[keep]
    activity_weight = system.module_tokenizer.activity_weight.detach().to(
        device=device, dtype=torch.float32
    )[torch.as_tensor(keep, device=device)]
    lift = system.module_tokenizer.activity_weight.detach().to(device=device)
    explicit_scale = system.precision_head.explicit_scale.detach().to(device=device)

    variant_names = np.asarray(VARIANTS, dtype=str)
    n_variants = len(VARIANTS)
    n_donors = donor_ids.size
    n_celltypes = len(celltype_names)
    n_modules = len(module_names)
    model_sum = np.zeros((n_variants, n_donors, n_celltypes, n_modules), dtype=np.float64)
    observed_sum = np.zeros((n_donors, n_celltypes, n_modules), dtype=np.float64)
    group_count = np.zeros((n_donors, n_celltypes), dtype=np.int64)
    donor_nll_sum = np.zeros((n_variants, n_donors), dtype=np.float64)
    donor_nll_count = np.zeros(n_donors, dtype=np.int64)
    donor_ct_nll_sum = np.zeros((n_variants, n_donors, n_celltypes), dtype=np.float64)
    donor_ct_nll_count = np.zeros((n_donors, n_celltypes), dtype=np.int64)
    component_coeff_rms_sum = {name: 0.0 for name in VARIANTS[1:]}
    component_gene_rms_sum = {name: 0.0 for name in VARIANTS[1:]}
    component_rms_count = 0

    loader = DataLoader(
        ds,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=device.type == "cuda",
        persistent_workers=args.num_workers > 0,
    )
    cap = len(loader) if args.max_batches <= 0 else min(len(loader), args.max_batches)
    ordinal_values = torch.arange(system.decoder.n_bins, device=device, dtype=torch.float32)
    print(
        f"[path-ablation] split={args.eval_split} cells={len(ds):,} "
        f"donors={n_donors} celltypes={n_celltypes} modules={n_modules} "
        f"batches={cap}/{len(loader)} device={device}",
        flush=True,
    )

    with torch.inference_mode():
        for batch_index, cpu_batch in enumerate(loader):
            if batch_index >= cap:
                break
            batch = move_batch(cpu_batch, device)
            output = system(
                batch,
                sample_latent=False,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )
            precision = output.precision_out
            if precision is None or output.decoder_out is None:
                raise RuntimeError("checkpoint did not produce PRISM/decoder outputs")

            zero = torch.zeros_like(precision.total_module_coeff)
            common = precision.common_axis_coeff.sum(dim=1)
            interaction = precision.interaction_pair_coeff.sum(dim=1)
            rank2 = precision.personal_rank2_coeff
            response = precision.response_axis_coeff.sum(dim=1)
            module_local = precision.module_local_coeff
            nonlinear = precision.module_local_nonlinear_coeff
            removed_coeff = {
                "full": zero,
                "common_pathology_off": common,
                "pathology_interaction_off": interaction,
                "all_shared_pathology_off": common + interaction,
                "personal_rank2_baseline_off": rank2,
                "personal_pathology_response_off": response,
                "module_local_off": module_local,
                "module_local_nonlinear_off": nonlinear,
                "all_personal_off": rank2 + response + module_local,
            }
            removed_score = {
                name: (coeff * explicit_scale.to(dtype=coeff.dtype))
                @ lift.to(dtype=coeff.dtype)
                for name, coeff in removed_coeff.items()
            }
            thresholds = output.decoder_out.thresholds
            expected = []
            cell_nll = []
            target = batch["y_ord"].long()
            for name in VARIANTS:
                score = output.decoder_out.score - removed_score[name]
                probabilities = system.decoder._score_to_probs(score, thresholds)[1]
                expected.append((probabilities * ordinal_values).sum(dim=-1))
                true_probability = probabilities.gather(
                    dim=-1, index=target.unsqueeze(-1)
                ).squeeze(-1)
                cell_nll.append(-torch.log(true_probability.clamp_min(1.0e-8)).mean(dim=1))
            expected_tensor = torch.stack(expected, dim=0).to(dtype=torch.float32)
            model_modules = torch.einsum(
                "vbg,mg->vbm", expected_tensor, activity_weight
            ).cpu().numpy()
            observed_modules = torch.einsum(
                "bg,mg->bm", target.float(), activity_weight
            ).cpu().numpy()
            cell_nll_array = torch.stack(cell_nll, dim=0).cpu().numpy()

            cpu_donor = cpu_batch["donor_id"].numpy().astype(np.int64)
            cpu_celltype = cpu_batch["celltype_id"].numpy().astype(np.int64)
            for local_index, (donor_id, celltype_id) in enumerate(
                zip(cpu_donor, cpu_celltype)
            ):
                dpos = int(donor_pos[donor_id])
                cpos = int(celltype_id)
                model_sum[:, dpos, cpos] += model_modules[:, local_index]
                observed_sum[dpos, cpos] += observed_modules[local_index]
                group_count[dpos, cpos] += 1
                donor_nll_sum[:, dpos] += cell_nll_array[:, local_index]
                donor_nll_count[dpos] += 1
                donor_ct_nll_sum[:, dpos, cpos] += cell_nll_array[:, local_index]
                donor_ct_nll_count[dpos, cpos] += 1

            batch_cells = int(target.shape[0])
            component_rms_count += batch_cells
            for name in VARIANTS[1:]:
                component_coeff_rms_sum[name] += float(
                    removed_coeff[name].float().square().mean(dim=1).sqrt().sum().cpu()
                )
                component_gene_rms_sum[name] += float(
                    removed_score[name].float().square().mean(dim=1).sqrt().sum().cpu()
                )
            if (batch_index + 1) % 100 == 0 or batch_index + 1 == cap:
                print(f"[path-ablation] {batch_index + 1}/{cap}", flush=True)

    valid = group_count > 0
    model_fp = np.full_like(model_sum, np.nan, dtype=np.float32)
    observed_fp = np.full_like(observed_sum, np.nan, dtype=np.float32)
    for dpos, cpos in np.argwhere(valid):
        model_fp[:, dpos, cpos] = (
            model_sum[:, dpos, cpos] / group_count[dpos, cpos]
        ).astype(np.float32)
        observed_fp[dpos, cpos] = (
            observed_sum[dpos, cpos] / group_count[dpos, cpos]
        ).astype(np.float32)
    donor_nll = donor_nll_sum / np.maximum(donor_nll_count[None, :], 1)
    donor_ct_nll = np.full_like(donor_ct_nll_sum, np.nan)
    for dpos, cpos in np.argwhere(donor_ct_nll_count > 0):
        donor_ct_nll[:, dpos, cpos] = (
            donor_ct_nll_sum[:, dpos, cpos] / donor_ct_nll_count[dpos, cpos]
        )

    donor_key = current._dataset_section().get("donor_key", "donor_id")
    donor_pathology, _, _, _ = build_donor_pathology(
        obs, ds.donor_vocab, donor_key, ds.row_idx
    )
    covariate_names = np.asarray(
        ["ADNC", "Braak", "Thal", "CERAD", "LATE", "Lewy", "Age", "Sex_male"],
        dtype=str,
    )
    donor_covariates = np.full((n_donors, len(covariate_names)), np.nan, dtype=np.float64)
    for dpos, donor_name in enumerate(donor_names):
        record = donor_pathology.get(str(donor_name), {})
        donor_covariates[dpos] = [float(record.get(name, np.nan)) for name in covariate_names]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_npz)), exist_ok=True)
    np.savez_compressed(
        args.out_npz,
        variant_names=variant_names,
        donor_names=donor_names,
        celltype_names=celltype_names,
        module_names=module_names,
        model_fp=model_fp,
        observed_fp=observed_fp,
        group_count=group_count,
        donor_nll=donor_nll,
        donor_nll_count=donor_nll_count,
        donor_ct_nll=donor_ct_nll,
        donor_ct_nll_count=donor_ct_nll_count,
        covariate_names=covariate_names,
        donor_covariates=donor_covariates,
    )
    report = {
        "kind": "frozen_checkpoint_conditional_path_reliance",
        "causal_interpretation": False,
        "technical_zero_inference_ablation": "excluded_not_an_inference_path",
        "evaluation_split": args.eval_split,
        "selected_cells": int(len(ds)),
        "total_validation_cells": int(n_val),
        "total_test_cells": int(n_test),
        "per_donor_celltype_cap": int(args.per_donor_celltype),
        "sample_seed": int(args.sample_seed),
        "variants": list(VARIANTS),
        "component_module_coefficient_rms": {
            name: component_coeff_rms_sum[name] / max(component_rms_count, 1)
            for name in VARIANTS[1:]
        },
        "component_gene_score_rms": {
            name: component_gene_rms_sum[name] / max(component_rms_count, 1)
            for name in VARIANTS[1:]
        },
        "checkpoint_compatibility": compatibility,
        "provenance": {
            "config": os.path.abspath(args.config),
            "config_sha256": sha256_file(args.config),
            "checkpoint": os.path.abspath(args.checkpoint),
            "checkpoint_sha256": sha256_file(args.checkpoint),
            "module_registry": os.path.abspath(registry_path),
            "module_registry_sha256": sha256_file(registry_path),
            "out_npz": os.path.abspath(args.out_npz),
        },
    }
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    print(f"[path-ablation] wrote {args.out_npz}", flush=True)
    print(f"[path-ablation] wrote {args.out_json}", flush=True)


if __name__ == "__main__":
    main()
