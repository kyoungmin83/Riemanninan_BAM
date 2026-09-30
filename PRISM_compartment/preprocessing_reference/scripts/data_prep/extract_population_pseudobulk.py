#!/usr/bin/env python3
"""Stream frozen-model outputs into donor x celltype x region sufficient statistics.

The shard command never stores cell-level latent values.  Several GPU shards can
be merged exactly because every output is an additive sum or count.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

try:
    from extract_v31a_allsplits import (
        ADNC,
        BRAAK,
        CERAD,
        LATE,
        THAL,
        _age,
        _codes,
        _lewy,
        _ord,
        build_all_dataset,
    )
    from disease_program_eval import move_batch
except ModuleNotFoundError:
    # Importable package form used by unit tests; direct script execution keeps
    # the historical top-level imports above.
    from scripts.extract_v31a_allsplits import (
        ADNC,
        BRAAK,
        CERAD,
        LATE,
        THAL,
        _age,
        _codes,
        _lewy,
        _ord,
        build_all_dataset,
    )
    from scripts.disease_program_eval import move_batch
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base_mod
from kmlee_bam.training.sealed_donor_allowlist import (
    dataset_train_donor_names,
    load_sealed_donor_allowlist,
    restrict_dataset_to_allowlist,
    validate_train_name_universe,
    validate_w44_source_checkpoint,
)


PATHOLOGY_NAMES = ("thal", "braak", "cerad", "late", "adnc", "lewy")
SPLIT_NAMES = ("train", "val", "test")


def _metadata(obs, region_key: str):
    def col(name):
        return obs[name].values if name in obs.columns else None

    pathology = np.column_stack(
        [
            _ord(col("Thal phase"), THAL),
            _ord(col("Braak stage"), BRAAK),
            _ord(col("CERAD score"), CERAD),
            _ord(col("LATE-NC stage"), LATE),
            _ord(col("ADNC"), ADNC),
            _lewy(col("Lewy body disease pathology")),
        ]
    )
    age_raw = col("Age at death")
    age = np.asarray([_age(value) for value in age_raw], dtype=np.float64)
    sex = _codes(col("sex")).astype(np.float64)
    ref = np.asarray(
        [1.0 if str(value).strip().lower() == "true" else 0.0 for value in col("Neurotypical reference")],
        dtype=np.float64,
    )
    split_text = np.asarray(col("split"), dtype=object).astype(str)
    split = np.asarray(
        [{name: index for index, name in enumerate(SPLIT_NAMES)}.get(value, -1) for value in split_text],
        dtype=np.int64,
    )
    region_text = np.asarray(col(region_key), dtype=object).astype(str)
    region_vocab = np.asarray(sorted(np.unique(region_text).tolist()), dtype=object)
    region_map = {str(value): index for index, value in enumerate(region_vocab)}
    region = np.asarray([region_map[value] for value in region_text], dtype=np.int64)
    return pathology, age, sex, ref, split, region, region_vocab


def _finite_add(target_sum, target_count, group, values):
    finite = np.isfinite(values)
    np.add.at(target_sum, group, np.where(finite, values, 0.0))
    np.add.at(target_count, group, finite.astype(np.int64))


def run_shard(args) -> None:
    torch.manual_seed(0)
    np.random.seed(0)
    raw = json.load(open(args.config, encoding="utf-8"))
    current.set_settings(raw)
    current.install_hooks()
    cfg = base_mod.load_config(args.config)
    ds, obs = build_all_dataset(cfg)
    pathology, age, sex, ref, split, region, region_vocab = _metadata(obs, args.region_key)
    all_rows = np.asarray(ds.row_idx, dtype=np.int64)
    sealed_allowlist = None
    source_checkpoint_sha256 = None
    source_checkpoint_provenance_sha256 = None
    if (
        args.donor_split_manifest is not None
        or args.canonical_fit_partition is not None
    ):
        if not args.donor_split_manifest or not args.canonical_fit_partition:
            raise ValueError(
                "--donor-split-manifest and --canonical-fit-partition must be "
                "supplied together"
            )
        if str(args.canonical_fit_partition).upper() != "W44":
            raise ValueError("canonical context fit partition must be W44")
        sealed_allowlist = load_sealed_donor_allowlist(
            args.donor_split_manifest,
            partition="W44",
        )
        donor_key = str(args.donor_key)
        if donor_key not in obs.columns:
            raise KeyError(f"context obs lacks donor key {donor_key!r}")
        validate_train_name_universe(
            obs.loc[split == 0, donor_key].astype(str).tolist(),
            sealed_allowlist,
        )
        source_checkpoint_sha256, source_checkpoint_provenance = validate_w44_source_checkpoint(
            args.checkpoint, sealed_allowlist
        )
        source_checkpoint_provenance_sha256 = source_checkpoint_provenance[
            "manifest_sha256"
        ]
        expected_checkpoint_sha = args.source_checkpoint_expected_sha256
        if (
            expected_checkpoint_sha is not None
            and str(expected_checkpoint_sha).lower() != source_checkpoint_sha256
        ):
            raise ValueError("source checkpoint expected SHA256 mismatch")
        # R10/C10 expression is permitted only as donor-specific support input;
        # it never fits weights or residualizer statistics.  Official val/test
        # rows are excluded before model inference.
        all_rows = all_rows[split[all_rows] == 0]
        ds.row_idx = all_rows.copy()
        fit_dataset = restrict_dataset_to_allowlist(ds, sealed_allowlist)
        if dataset_train_donor_names(fit_dataset) != sealed_allowlist.donor_names:
            raise RuntimeError("context residualizer fit rows differ from W44")
    else:
        fit_dataset = ds
    shards = np.array_split(all_rows, args.num_shards)
    celltype_vocab = np.asarray(getattr(ds.spec, "celltype_vocab", []), dtype=object)
    donor_vocab = np.asarray(getattr(ds, "donor_vocab", []), dtype=object)
    if len(donor_vocab) == 0:
        donor_vocab = np.asarray(sorted(obs["donor_id"].astype(str).unique().tolist()), dtype=object)
    n_donor = len(donor_vocab)
    n_celltype = len(celltype_vocab)
    n_region = len(region_vocab)
    n_group = n_donor * n_celltype * n_region

    system = current.build_system(cfg, fit_dataset)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = payload.get("system_state_dict", payload) if isinstance(payload, dict) else payload
    system.load_state_dict(state, strict=False)
    # Only after system construction has sampled/inspected the complete fitting
    # view may the inference dataset be narrowed to this additive shard.
    ds.row_idx = np.asarray(shards[args.shard_index], dtype=np.int64)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    system.to(device).eval()
    loader = DataLoader(
        ds,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.workers,
        pin_memory=True,
        persistent_workers=args.workers > 0,
    )

    count = np.zeros(n_group, dtype=np.int64)
    age_sum = np.zeros(n_group, dtype=np.float64)
    age_count = np.zeros(n_group, dtype=np.int64)
    sex_sum = np.zeros(n_group, dtype=np.float64)
    sex_count = np.zeros(n_group, dtype=np.int64)
    ref_sum = np.zeros(n_group, dtype=np.float64)
    pathology_sum = np.zeros((n_group, len(PATHOLOGY_NAMES)), dtype=np.float64)
    pathology_count = np.zeros_like(pathology_sum, dtype=np.int64)
    split_count = np.zeros((n_group, len(SPLIT_NAMES)), dtype=np.int64)
    tech_count = np.zeros((n_group, 8), dtype=np.int64)
    z_sum = None
    module_sum = None

    with torch.no_grad():
        for batch_index, batch in enumerate(loader):
            output = system(
                move_batch(batch, device),
                sample_latent=False,
                return_all_hidden_states=False,
                return_attn_diagnostics=False,
            )
            row = batch["row_index"].numpy().astype(np.int64)
            donor = batch["donor_id"].numpy().astype(np.int64)
            celltype = batch["celltype_id"].numpy().astype(np.int64)
            tech = batch.get("tech_id", batch.get("batch_id")).numpy().astype(np.int64)
            if donor.max(initial=-1) >= n_donor:
                raise ValueError("donor id exceeds donor vocabulary")
            if tech.max(initial=-1) >= tech_count.shape[1]:
                raise ValueError("increase preallocated tech classes")
            group = (donor * n_celltype + celltype) * n_region + region[row]
            latent = output.z_perp.detach().cpu().numpy().astype(np.float64)
            module = getattr(output, "module_activity", None)
            module = None if module is None else module.detach().cpu().numpy().astype(np.float64)
            if z_sum is None:
                z_sum = np.zeros((n_group, latent.shape[1]), dtype=np.float64)
                module_sum = np.zeros((n_group, 0 if module is None else module.shape[1]), dtype=np.float64)
            np.add.at(z_sum, group, latent)
            if module is not None:
                np.add.at(module_sum, group, module)
            np.add.at(count, group, 1)
            _finite_add(age_sum, age_count, group, age[row])
            _finite_add(sex_sum, sex_count, group, sex[row])
            np.add.at(ref_sum, group, ref[row])
            _finite_add(pathology_sum, pathology_count, group, pathology[row])
            for split_id in range(len(SPLIT_NAMES)):
                np.add.at(split_count[:, split_id], group, split[row] == split_id)
            for tech_id in np.unique(tech):
                if tech_id >= 0:
                    np.add.at(tech_count[:, tech_id], group, tech == tech_id)
            if batch_index % 100 == 0:
                print(
                    f"[population-pseudobulk] shard={args.shard_index}/{args.num_shards} "
                    f"cells={min((batch_index + 1) * args.batch_size, len(ds.row_idx))}/{len(ds.row_idx)}",
                    flush=True,
                )

    module_names = []
    try:
        registry = json.load(open(cfg.module_tokenizer.registry_json_path, encoding="utf-8"))
        module_names = registry.get("module_names", [])
    except Exception as exc:
        print(f"[population-pseudobulk] WARNING module names unavailable: {exc}", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        z_sum=z_sum,
        module_sum=module_sum,
        count=count,
        age_sum=age_sum,
        age_count=age_count,
        sex_sum=sex_sum,
        sex_count=sex_count,
        ref_sum=ref_sum,
        pathology_sum=pathology_sum,
        pathology_count=pathology_count,
        split_count=split_count,
        tech_count=tech_count,
        region_vocab=region_vocab,
        celltype_vocab=celltype_vocab,
        donor_vocab=donor_vocab,
        module_names=np.asarray(module_names, dtype=object),
        pathology_names=np.asarray(PATHOLOGY_NAMES, dtype=object),
        n_source_rows=np.asarray(len(ds.row_idx), dtype=np.int64),
        **(
            {
                "allowed_fit_donor_names": np.asarray(
                    sealed_allowlist.donor_names, dtype=str
                ),
                "observed_residualizer_fit_donor_names": np.asarray(
                    sealed_allowlist.donor_names, dtype=str
                ),
                "observed_source_model_fit_donor_names": np.asarray(
                    sealed_allowlist.donor_names, dtype=str
                ),
                "donor_split_manifest_sha256": np.asarray(
                    sealed_allowlist.manifest_sha256
                ),
                "extractor_checkpoint_sha256": np.asarray(
                    source_checkpoint_sha256
                ),
                "extractor_checkpoint_provenance_sha256": np.asarray(
                    source_checkpoint_provenance_sha256
                ),
                "support_donor_scope": np.asarray("train64"),
                "row_level_allowlist_enforced": np.asarray(True),
                "official_validation_used": np.asarray(False),
                "official_test_used": np.asarray(False),
            }
            if sealed_allowlist is not None
            else {}
        ),
    )
    print(f"[population-pseudobulk] wrote {args.out}", flush=True)


def merge_shards(args) -> None:
    files = [np.load(path, allow_pickle=True) for path in args.inputs]
    additive = (
        "z_sum", "module_sum", "count", "age_sum", "age_count", "sex_sum",
        "sex_count", "ref_sum", "pathology_sum", "pathology_count", "split_count", "tech_count",
    )
    merged = {name: sum((item[name] for item in files)) for name in additive}
    count = merged["count"]
    keep = count > 0
    first = files[0]
    n_celltype = len(first["celltype_vocab"])
    n_region = len(first["region_vocab"])
    flat = np.where(keep)[0]
    donor = flat // (n_celltype * n_region)
    remainder = flat % (n_celltype * n_region)
    celltype = remainder // n_region
    region_id = remainder % n_region

    def mean(total, denominator):
        shape = (slice(None),) + (None,) * (total.ndim - 1)
        den = denominator[shape] if denominator.ndim == 1 else denominator
        return np.divide(total, den, out=np.full_like(total, np.nan, dtype=np.float64), where=den > 0)

    group_count = count[keep]
    split_id = np.argmax(merged["split_count"][keep], axis=1)
    pathology = mean(merged["pathology_sum"][keep], merged["pathology_count"][keep])
    tech_count = merged["tech_count"][keep]
    active_tech = np.where(tech_count.sum(axis=0) > 0)[0]
    tech_fraction = tech_count[:, active_tech] / group_count[:, None]
    save = {
        "aggregation_level": np.asarray("donor_celltype_region"),
        "z": merged["z_sum"][keep] / group_count[:, None],
        "module_activity": merged["module_sum"][keep] / group_count[:, None],
        "donor": donor.astype(np.int64),
        "celltype": celltype.astype(np.int64),
        "region": first["region_vocab"][region_id],
        "split": np.asarray(SPLIT_NAMES, dtype=object)[split_id],
        "age": mean(merged["age_sum"][keep], merged["age_count"][keep]),
        "sex": mean(merged["sex_sum"][keep], merged["sex_count"][keep]),
        "isref": merged["ref_sum"][keep] / group_count,
        "cell_count": group_count,
        "tech_fraction": tech_fraction,
        "tech_classes": active_tech.astype(np.int64),
        "celltype_vocab": first["celltype_vocab"],
        "donor_vocab": first["donor_vocab"],
        "module_names": first["module_names"],
        "n_source_cells": np.asarray(sum(int(item["n_source_rows"]) for item in files), dtype=np.int64),
    }
    canonical_fields = (
        "allowed_fit_donor_names",
        "observed_residualizer_fit_donor_names",
        "observed_source_model_fit_donor_names",
        "donor_split_manifest_sha256",
        "extractor_checkpoint_sha256",
        "extractor_checkpoint_provenance_sha256",
        "support_donor_scope",
        "row_level_allowlist_enforced",
        "official_validation_used",
        "official_test_used",
    )
    canonical_present = [
        all(field in item.files for field in canonical_fields) for item in files
    ]
    if any(canonical_present) and not all(canonical_present):
        raise ValueError("cannot merge canonical and legacy context shards")
    if all(canonical_present):
        for field in canonical_fields:
            reference = first[field]
            if any(not np.array_equal(item[field], reference) for item in files[1:]):
                raise ValueError(
                    f"canonical context shard provenance differs for {field}"
                )
            save[field] = reference
    for index, name in enumerate(PATHOLOGY_NAMES):
        save[name] = pathology[:, index]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **save)
    print(
        f"[population-pseudobulk] merged cells={int(save['n_source_cells'])} "
        f"groups={len(group_count)} -> {args.out}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    shard = sub.add_parser("shard")
    shard.add_argument("--config", required=True)
    shard.add_argument("--checkpoint", required=True)
    shard.add_argument("--out", required=True)
    shard.add_argument("--shard-index", type=int, required=True)
    shard.add_argument("--num-shards", type=int, required=True)
    shard.add_argument("--region-key", default="brain_region")
    shard.add_argument("--batch-size", type=int, default=256)
    shard.add_argument("--workers", type=int, default=4)
    shard.add_argument("--donor-key", default="donor_id")
    shard.add_argument("--donor-split-manifest")
    shard.add_argument("--canonical-fit-partition", choices=("W44",))
    shard.add_argument("--source-checkpoint-expected-sha256")
    shard.set_defaults(func=run_shard)
    merge = sub.add_parser("merge")
    merge.add_argument("--inputs", nargs="+", required=True)
    merge.add_argument("--out", required=True)
    merge.set_defaults(func=merge_shards)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
