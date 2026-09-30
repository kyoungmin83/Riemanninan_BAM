#!/usr/bin/env python3
"""Build train-donor-only ordinal module statistics for PRISM v2.

The artifact contains no learned model output. It streams the train split,
aggregates observed ordinal tiers to donor x cell-type pseudobulks, projects
those pseudobulks with the tokenizer's exact activity dictionary, and fits
cell-type-specific robust scales. Validation and test cells are never opened
by this script. Donor centering itself is deliberately *not* stored as a fixed
mean: prediction and observation must be centered separately inside each
cross-donor contrast. The artifact also records a train-only donor×celltype
support-presence mask for the donor-contrast sampler.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

from kmlee_bam.data.ordinal_dataset import read_zarr_dataframe
from kmlee_bam.training import run_current as current
from kmlee_bam.training import runner_base as base
from kmlee_bam.training.sealed_donor_allowlist import (
    bind_allowlist_to_dataset,
    dataset_train_donor_names,
    load_sealed_donor_allowlist,
)


SCHEMA_VERSION = "kmlee_bam.prism_module_rescue_stats.v2"


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_activity_dictionary(cfg, gene_names: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    registry_path = Path(cfg.module_tokenizer.registry_json_path)
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry_genes = np.asarray(registry["gene_names"], dtype=object).astype(str)
    if not np.array_equal(registry_genes, gene_names.astype(str)):
        raise ValueError("registry gene order does not match the ordinal dataset")

    module_names = np.asarray(registry["module_names"], dtype=object).astype(str)
    membership = np.asarray(registry["membership_binary"], dtype=np.float32)
    if membership.shape != (len(module_names), len(gene_names)):
        raise ValueError(
            "registry membership shape mismatch: "
            f"{membership.shape} vs {(len(module_names), len(gene_names))}"
        )

    activity_path = Path(cfg.module_tokenizer.activity_weight_path)
    with np.load(activity_path, allow_pickle=False) as archive:
        key = str(cfg.module_tokenizer.activity_weight_key)
        if key not in archive:
            raise KeyError(f"{activity_path} does not contain activity key {key!r}")
        activity = np.asarray(archive[key], dtype=np.float32)
    if activity.shape != membership.shape:
        raise ValueError(
            f"activity dictionary shape mismatch: {activity.shape} vs {membership.shape}"
        )
    if not np.isfinite(activity).all():
        raise ValueError("activity dictionary contains non-finite values")

    activity *= membership > 0
    normalization = str(cfg.module_tokenizer.activity_weight_normalization)
    if normalization == "l2":
        denominator = np.sqrt(np.square(activity, dtype=np.float64).sum(axis=1))
        if np.any(denominator <= 0):
            raise ValueError("activity dictionary contains an empty row")
        activity /= denominator[:, None].astype(np.float32)
    elif normalization == "row_sum":
        if np.any(activity < 0):
            raise ValueError("row_sum normalization requires non-negative activity weights")
        denominator = activity.sum(axis=1)
        if np.any(denominator <= 0):
            raise ValueError("activity dictionary contains an empty row")
        activity /= denominator[:, None]
    elif normalization != "none":
        raise ValueError(f"unknown activity normalization: {normalization!r}")

    return activity, membership.sum(axis=1).astype(np.int64), module_names


def _robust_sd(values: np.ndarray, *, axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Return median and Gaussian-consistent MAD scale."""

    location = np.median(values, axis=axis)
    mad = np.median(np.abs(values - np.expand_dims(location, axis=axis)), axis=axis)
    return location, 1.4826 * mad


def canonical_stat_fit_donor_mask(
    dataset,
    *,
    donor_split_manifest: str | Path | None,
    canonical_fit_partition: str | None,
) -> tuple[np.ndarray, object | None]:
    """Resolve the donor rows permitted to fit module locations/scales."""

    n_donor = int(len(dataset.donor_vocab))
    fit_mask = np.ones(n_donor, dtype=bool)
    if donor_split_manifest is None and canonical_fit_partition is None:
        return fit_mask, None
    if not donor_split_manifest or not canonical_fit_partition:
        raise ValueError(
            "donor split manifest and canonical fit partition are required together"
        )
    if str(canonical_fit_partition).upper() != "W44":
        raise ValueError("canonical module-stat fit partition must be W44")
    allowlist = load_sealed_donor_allowlist(
        donor_split_manifest,
        partition="W44",
        expected_train_donor_names=dataset_train_donor_names(dataset),
    )
    fit_ids = bind_allowlist_to_dataset(dataset, allowlist)
    fit_mask.fill(False)
    fit_mask[np.asarray(fit_ids, dtype=np.int64)] = True
    return fit_mask, allowlist


def build(args: argparse.Namespace) -> None:
    for name in (
        "batch_size",
        "min_cells_group",
        "min_support_cells",
        "min_donors_celltype",
        "log_every",
    ):
        if int(getattr(args, name)) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if int(args.num_workers) < 0:
        raise ValueError("--num-workers must be non-negative")
    if float(args.scale_floor) <= 0.0:
        raise ValueError("--scale-floor must be positive")

    config_path = Path(args.config).resolve()
    raw = json.loads(config_path.read_text(encoding="utf-8"))
    # This artifact needs only the observed ordinal target.  Production
    # training also asks the dataset to construct x_gene_scalar, x_log1p, and
    # a stochastic thinned companion view; doing that for 2M cells would add
    # substantial I/O/CPU work without changing a single statistic below.
    raw.setdefault("data", {})["return_log_input"] = False
    raw["data"]["return_x_gene_scalar"] = False
    raw.setdefault("thinning", {})["enabled"] = False
    if args.donor_split_manifest is not None:
        # This producer must retain raw R10/C10 targets while fitting all
        # location/scale statistics on W44.  Its own row contract below is
        # therefore authoritative; do not let the source-training dataset hook
        # pre-filter the raw target table to W44.
        raw.setdefault("canonical_w44_source", {})["enabled"] = False
    current.set_settings(raw)
    current.install_hooks()
    cfg = base.load_config(str(config_path))
    cfg.data.return_log_input = False
    cfg.data.return_x_gene_scalar = False
    train_ds, _, _ = current.build_datasets(cfg)

    obs = read_zarr_dataframe(train_ds.root, "obs")
    split_values = np.asarray(obs[cfg.data.split_key].astype(str), dtype=str)
    selected_splits = np.unique(split_values[np.asarray(train_ds.row_idx, dtype=np.int64)])
    if selected_splits.tolist() != ["train"]:
        raise RuntimeError(f"expected train-only dataset, found splits={selected_splits.tolist()}")

    gene_names = np.asarray(train_ds.spec.gene_names, dtype=object).astype(str)
    celltype_names = np.asarray(train_ds.spec.celltype_vocab, dtype=object).astype(str)
    donor_names = np.asarray(train_ds.donor_vocab, dtype=object).astype(str)
    activity, module_size, module_names = _load_activity_dictionary(cfg, gene_names)

    n_donor = len(donor_names)
    n_celltype = len(celltype_names)
    n_gene = len(gene_names)
    n_group = n_donor * n_celltype
    fit_donor_mask, sealed_allowlist = canonical_stat_fit_donor_mask(
        train_ds,
        donor_split_manifest=args.donor_split_manifest,
        canonical_fit_partition=args.canonical_fit_partition,
    )
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false")

    loader = DataLoader(
        train_ds,
        batch_size=int(args.batch_size),
        shuffle=False,
        num_workers=int(args.num_workers),
        pin_memory=device.type == "cuda",
        persistent_workers=int(args.num_workers) > 0,
    )
    gene_sum = torch.zeros((n_group, n_gene), dtype=torch.float32, device=device)
    group_count = torch.zeros(n_group, dtype=torch.float64, device=device)

    with torch.no_grad():
        for batch_index, batch in enumerate(loader, start=1):
            donor = batch["donor_id"].to(device=device, dtype=torch.long, non_blocking=True)
            celltype = batch["celltype_id"].to(
                device=device, dtype=torch.long, non_blocking=True
            )
            ordinal = batch["y_ord"].to(
                device=device, dtype=torch.float32, non_blocking=True
            )
            group = donor * n_celltype + celltype
            gene_sum.index_add_(0, group, ordinal)
            group_count.index_add_(
                0,
                group,
                torch.ones(group.shape[0], dtype=torch.float64, device=device),
            )
            if batch_index % int(args.log_every) == 0:
                seen = min(batch_index * int(args.batch_size), len(train_ds))
                print(
                    f"[prism-module-stats] cells={seen:,}/{len(train_ds):,}",
                    flush=True,
                )

        observed = group_count >= int(args.min_cells_group)
        gene_mean = gene_sum[observed] / group_count[observed].to(torch.float32).unsqueeze(1)
        activity_tensor = torch.from_numpy(activity).to(device=device)
        observed_module = gene_mean @ activity_tensor.transpose(0, 1)

    observed_index = torch.nonzero(observed, as_tuple=False).flatten().cpu().numpy()
    observed_module_np = observed_module.cpu().numpy().astype(np.float32)
    group_count_np = group_count.cpu().numpy().astype(np.int64)
    group_module = np.full((n_group, len(module_names)), np.nan, dtype=np.float32)
    group_module[observed_index] = observed_module_np
    group_module = group_module.reshape(n_donor, n_celltype, len(module_names))
    group_count_2d = group_count_np.reshape(n_donor, n_celltype)
    support_presence = group_count_2d >= int(args.min_support_cells)
    support_context_celltype_ids = np.arange(n_celltype, dtype=np.int64)
    donor_has_non_target_support = np.zeros_like(support_presence)
    for celltype in range(n_celltype):
        donor_has_non_target_support[:, celltype] = np.any(
            support_presence[
                :,
                support_context_celltype_ids != int(celltype),
            ],
            axis=1,
        )
    rescue_eligible_group = (
        (group_count_2d >= int(args.min_cells_group))
        & donor_has_non_target_support
    )

    shrinkage = float(args.scale_shrinkage)
    if not 0.0 <= shrinkage <= 1.0:
        raise ValueError("--scale-shrinkage must lie in [0, 1]")

    # Raw donor×celltype targets remain available for R10/C10 evaluation, but
    # every fitted location/scale is estimated from W44 rows only.
    valid_rows = group_module[fit_donor_mask].reshape(-1, len(module_names))
    valid_rows = valid_rows[np.isfinite(valid_rows).all(axis=1)]
    if len(valid_rows) < int(args.min_donors_celltype):
        raise RuntimeError("too few eligible donor×celltype groups for pooled scale")
    pooled_location, pooled_scale_raw = _robust_sd(valid_rows, axis=0)
    pooled_scale = np.maximum(pooled_scale_raw, float(args.scale_floor)).astype(
        np.float32
    )

    celltype_location = np.zeros(
        (n_celltype, len(module_names)), dtype=np.float32
    )
    celltype_scale_raw = np.zeros_like(celltype_location)
    celltype_scale = np.ones_like(celltype_location)
    n_fit_donor = np.zeros(n_celltype, dtype=np.int64)
    for celltype in range(n_celltype):
        fit = fit_donor_mask & (
            group_count_2d[:, celltype] >= int(args.min_cells_group)
        )
        values = group_module[fit, celltype]
        n_fit_donor[celltype] = int(fit.sum())
        if int(fit.sum()) < int(args.min_donors_celltype):
            raise RuntimeError(
                f"celltype {celltype_names[celltype]} has only {int(fit.sum())} "
                "train donors with sufficient cells"
            )
        location, raw_scale = _robust_sd(values, axis=0)
        shrunk_scale = (1.0 - shrinkage) * raw_scale + shrinkage * pooled_scale
        celltype_location[celltype] = location.astype(np.float32)
        celltype_scale_raw[celltype] = raw_scale.astype(np.float32)
        celltype_scale[celltype] = np.maximum(
            np.where(np.isfinite(shrunk_scale), shrunk_scale, pooled_scale),
            float(args.scale_floor),
        ).astype(np.float32)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "source_split": "train_only",
        "n_train_cells": int(len(train_ds)),
        "n_donors_total_vocab": int(n_donor),
        "n_celltypes": int(n_celltype),
        "n_modules": int(len(module_names)),
        "n_genes": int(n_gene),
        "min_cells_group": int(args.min_cells_group),
        "min_support_cells": int(args.min_support_cells),
        "min_donors_celltype": int(args.min_donors_celltype),
        "scale_floor": float(args.scale_floor),
        "scale_estimator": "1.4826_x_median_absolute_deviation",
        "scale_shrinkage_to_pooled": shrinkage,
        "fixed_center_used_by_loss": False,
        "activity_normalization": str(
            cfg.module_tokenizer.activity_weight_normalization
        ),
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "spec_path": str(Path(cfg.data.spec_path).resolve()),
        "spec_sha256": _sha256(cfg.data.spec_path),
        "registry_path": str(Path(cfg.module_tokenizer.registry_json_path).resolve()),
        "registry_sha256": _sha256(cfg.module_tokenizer.registry_json_path),
        "activity_path": str(Path(cfg.module_tokenizer.activity_weight_path).resolve()),
        "activity_sha256": _sha256(cfg.module_tokenizer.activity_weight_path),
    }
    observed_statistic_donor_names = tuple(
        str(value)
        for value in donor_names[
            fit_donor_mask & np.any(group_count_2d > 0, axis=1)
        ].tolist()
    )
    if sealed_allowlist is not None:
        if observed_statistic_donor_names != sealed_allowlist.donor_names:
            raise RuntimeError(
                "module-stat fitting rows do not cover exactly the sealed W44 "
                "donors"
            )
        metadata.update(
            {
                "canonical_fit_partition": "W44",
                "donor_split_manifest_path": sealed_allowlist.manifest_path,
                "donor_split_manifest_sha256": sealed_allowlist.manifest_sha256,
                "allowed_fit_donor_names": list(sealed_allowlist.donor_names),
                "observed_statistic_donor_names": list(
                    observed_statistic_donor_names
                ),
                "row_level_allowlist_enforced": True,
                "raw_target_support_scope": "train64",
                "official_validation_used": False,
                "official_test_used": False,
            }
        )

    output = Path(args.out).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        schema_version=np.asarray(SCHEMA_VERSION),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
        celltype_location=celltype_location,
        celltype_scale_raw=celltype_scale_raw,
        celltype_scale=celltype_scale,
        pooled_location=pooled_location.astype(np.float32),
        pooled_scale=pooled_scale,
        group_module=group_module,
        group_count=group_count_2d,
        support_presence_by_donor_celltype=support_presence,
        support_context_celltype_ids=support_context_celltype_ids,
        donor_has_non_target_support=donor_has_non_target_support,
        rescue_eligible_group=rescue_eligible_group,
        n_fit_donor=n_fit_donor,
        module_size=module_size,
        module_names=np.asarray(module_names, dtype=str),
        celltype_names=np.asarray(celltype_names, dtype=str),
        donor_names=np.asarray(donor_names, dtype=str),
        gene_names=np.asarray(gene_names, dtype=str),
        **(
            {
                "allowed_fit_donor_names": np.asarray(
                    sealed_allowlist.donor_names, dtype=str
                ),
                "observed_statistic_donor_names": np.asarray(
                    observed_statistic_donor_names, dtype=str
                ),
                "donor_split_manifest_sha256": np.asarray(
                    sealed_allowlist.manifest_sha256
                ),
                "row_level_allowlist_enforced": np.asarray(True),
                "official_validation_used": np.asarray(False),
                "official_test_used": np.asarray(False),
            }
            if sealed_allowlist is not None
            else {}
        ),
    )
    print(
        f"[prism-module-stats] wrote {output} "
        f"sha256={_sha256(output)} train_cells={len(train_ds):,}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--min-cells-group", type=int, default=10)
    parser.add_argument("--min-support-cells", type=int, default=10)
    parser.add_argument("--min-donors-celltype", type=int, default=12)
    parser.add_argument("--scale-floor", type=float, default=0.05)
    parser.add_argument("--scale-shrinkage", type=float, default=0.10)
    parser.add_argument("--log-every", type=int, default=200)
    parser.add_argument("--donor-split-manifest")
    parser.add_argument("--canonical-fit-partition", choices=("W44",))
    build(parser.parse_args())


if __name__ == "__main__":
    main()
