#!/usr/bin/env python3
"""Stream train-only module activities into deterministic split-half summaries.

This is the memory-bounded precursor to
``build_prism_module_local_reliability_20260825.py``.  It never constructs a
validation/test dataset and never stores cell-level expression.  Each shard
reads a disjoint slice of the configured training rows, computes the exact
``GeneModuleTokenizer.module_activity`` input score, and saves additive sums.
The merge command verifies the shard partition and emits the pre-aggregated
half-pseudobulk schema consumed by the sealed reliability builder.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from kmlee_bam.data.dataset_with_covariates import V3OrdinalScDataset
from kmlee_bam.data.donor_context_dataset import DonorContextTable


SCHEMA_SHARD = "kmlee_bam.module_local_split_half_shard.v1"
SCHEMA_MERGED = "kmlee_bam.module_local_split_half_input.v1"
SPLIT_RULE = "blake2b64(seed|cell_id)_least_significant_bit"


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_int64(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values, dtype=np.int64)
    return hashlib.sha256(array.tobytes()).hexdigest()


def _deterministic_half(row_indices: np.ndarray, seed: int) -> np.ndarray:
    result = np.empty(len(row_indices), dtype=np.int8)
    for index, row in enumerate(row_indices):
        payload = f"{int(seed)}|{int(row)}".encode("utf-8")
        digest = hashlib.blake2b(payload, digest_size=8).digest()
        result[index] = int.from_bytes(digest, "little") & 1
    return result


def _dataset_from_config(raw: dict) -> V3OrdinalScDataset:
    data = raw["data"]
    extra = raw.get("dataset") or {}
    return V3OrdinalScDataset(
        zarr_path=data["zarr_path"],
        spec_path=data["spec_path"],
        split="train",
        matrix_path=data.get("matrix_path", "auto"),
        counts_layer=data.get("counts_layer", "counts"),
        prefer_raw=bool(data.get("prefer_raw", True)),
        split_key=data.get("split_key", "split"),
        celltype_key=data.get("celltype_key", "Subclass"),
        batch_key=data.get("batch_key", "donor_id"),
        return_log_input=False,
        reference_scaler_path=data["reference_scaler_path"],
        return_x_gene_scalar=True,
        x_gene_scalar_clip=data.get("x_gene_scalar_clip", 10.0),
        x_gene_scalar_eps=data.get("x_gene_scalar_eps", 1.0e-4),
        max_csr_data_span=data.get("max_csr_data_span", 50_000_000),
        max_span_to_selected_nnz_ratio=data.get(
            "max_span_to_selected_nnz_ratio", 20.0
        ),
        donor_key=extra.get("donor_key", "donor_id"),
        depth_key=extra.get("depth_key", "Number of UMIs"),
        detected_genes_key=extra.get("detected_genes_key", "Genes detected"),
    )


def _activity_weight(raw: dict, dataset: V3OrdinalScDataset) -> tuple[torch.Tensor, np.ndarray]:
    module = raw["module_tokenizer"]
    registry_path = Path(module["registry_json_path"])
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    gene_names = np.asarray(registry["gene_names"], dtype=object).astype(str)
    expected_genes = np.asarray(dataset.spec.gene_names, dtype=object).astype(str)
    if not np.array_equal(gene_names, expected_genes):
        raise ValueError("registry and ordinal-spec gene order differ")
    membership = torch.as_tensor(
        np.asarray(registry["membership_binary"]), dtype=torch.float32
    )
    activity_path = Path(module["activity_weight_path"])
    with np.load(activity_path, allow_pickle=False) as archive:
        activity = torch.as_tensor(
            np.asarray(archive[module.get("activity_weight_key", "activity_weight")]),
            dtype=torch.float32,
        )
    if activity.shape != membership.shape:
        raise ValueError("activity dictionary and module membership shapes differ")
    activity = activity * (membership > 0).to(activity.dtype)
    normalization = module.get("activity_weight_normalization", "l2")
    if normalization == "l2":
        activity = activity / activity.square().sum(1, keepdim=True).sqrt().clamp_min(
            1.0e-12
        )
    elif normalization == "row_sum":
        if bool((activity < 0).any()):
            raise ValueError("row_sum normalization cannot be used with signed weights")
        activity = activity / activity.sum(1, keepdim=True).clamp_min(1.0e-12)
    elif normalization != "none":
        raise ValueError(f"unknown activity normalization: {normalization}")
    if not bool(torch.isfinite(activity).all()) or bool(
        (activity.abs().sum(1) <= 0).any()
    ):
        raise ValueError("normalized activity dictionary is empty or non-finite")
    module_names = np.asarray(registry["module_names"], dtype=object).astype(str)
    if len(module_names) != activity.shape[0]:
        raise ValueError("module name count differs from activity dictionary")
    return activity.contiguous(), module_names


def _validate_layout(
    dataset: V3OrdinalScDataset,
    table: DonorContextTable,
    module_names: np.ndarray,
) -> None:
    if not np.array_equal(
        np.asarray(dataset.donor_vocab).astype(str), table.donor_names.astype(str)
    ):
        raise ValueError("dataset and source-context donor order differ")
    if not np.array_equal(
        np.asarray(dataset.spec.celltype_vocab).astype(str),
        table.celltype_names.astype(str),
    ):
        raise ValueError("dataset and source-context cell-type order differ")
    if not np.array_equal(
        np.asarray(dataset.region_vocab).astype(str), table.region_names.astype(str)
    ):
        raise ValueError("dataset and source-context region order differ")
    if not np.array_equal(module_names.astype(str), table.module_names.astype(str)):
        raise ValueError("activity and source-context module order differ")
    train_rows = np.asarray(dataset.row_idx, dtype=np.int64)
    train_donors = np.unique(np.asarray(dataset.donor_ids)[train_rows])
    expected_train = np.flatnonzero(table.donor_split.astype(str) == "train")
    if not np.array_equal(train_donors, expected_train):
        raise ValueError("training rows do not equal the source-context train donors")


def run_shard(args: argparse.Namespace) -> None:
    if args.num_shards <= 0 or not 0 <= args.shard_index < args.num_shards:
        raise ValueError("invalid shard index/count")
    raw = json.loads(Path(args.config).read_text(encoding="utf-8"))
    if raw.get("train", {}).get("skip_final_test") is not True:
        raise ValueError("split-half extraction requires skip_final_test=true")
    dataset = _dataset_from_config(raw)
    table = DonorContextTable.from_npz(
        args.source_context,
        min_cells=int(raw["precision_medicine"].get("min_context_cells", 10)),
    )
    activity_weight, module_names = _activity_weight(raw, dataset)
    _validate_layout(dataset, table, module_names)

    all_rows = np.asarray(dataset.row_idx, dtype=np.int64)
    rows = np.asarray(
        np.array_split(all_rows, int(args.num_shards))[int(args.shard_index)],
        dtype=np.int64,
    )
    halves = _deterministic_half(rows, int(args.seed))
    n_donor = len(table.donor_names)
    n_context = len(table.context_names)
    n_module = len(module_names)
    sums = np.zeros((2, n_donor * n_context, n_module), dtype=np.float64)
    counts = np.zeros((2, n_donor * n_context), dtype=np.int64)

    device = torch.device(args.device)
    weight_t = activity_weight.to(device=device)
    for start in range(0, len(rows), int(args.batch_size)):
        stop = min(start + int(args.batch_size), len(rows))
        batch_rows = rows[start:stop]
        celltype = np.asarray(dataset.celltype_ids)[batch_rows].astype(np.int64)
        donor = np.asarray(dataset.donor_ids)[batch_rows].astype(np.int64)
        region = np.asarray(dataset.region_ids)[batch_rows].astype(np.int64)
        context = celltype * len(table.region_names) + region
        group = donor * n_context + context
        expression = dataset.X.rows_to_dense(batch_rows).astype(np.float32, copy=False)
        np.log1p(expression, out=expression)
        expression -= dataset._scaler_mean[celltype]
        expression /= dataset._scaler_std_eps[celltype]
        clip = dataset.x_gene_scalar_clip
        if clip is not None:
            np.clip(expression, -float(clip), float(clip), out=expression)
        with torch.no_grad():
            score = torch.matmul(
                torch.from_numpy(expression).to(device=device), weight_t.T
            ).cpu().numpy()
        for side in (0, 1):
            selected = halves[start:stop] == side
            if not bool(selected.any()):
                continue
            np.add.at(sums[side], group[selected], score[selected].astype(np.float64))
            np.add.at(counts[side], group[selected], 1)
        if start == 0 or stop == len(rows) or stop % (100 * int(args.batch_size)) == 0:
            print(
                f"[module-local split-half] shard={args.shard_index}/{args.num_shards} "
                f"cells={stop}/{len(rows)}",
                flush=True,
            )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        schema_version=np.asarray(SCHEMA_SHARD),
        module_sum_half_a=sums[0],
        module_sum_half_b=sums[1],
        cell_count_half_a=counts[0],
        cell_count_half_b=counts[1],
        donor_names=table.donor_names.astype(str),
        donor_split=table.donor_split.astype(str),
        context_names=table.context_names.astype(str),
        module_names=module_names.astype(str),
        split_seed=np.asarray(int(args.seed), dtype=np.int64),
        split_rule=np.asarray(SPLIT_RULE),
        shard_index=np.asarray(int(args.shard_index), dtype=np.int64),
        num_shards=np.asarray(int(args.num_shards), dtype=np.int64),
        shard_row_count=np.asarray(len(rows), dtype=np.int64),
        total_train_row_count=np.asarray(len(all_rows), dtype=np.int64),
        shard_row_index_sha256=np.asarray(_sha256_int64(rows)),
        config_sha256=np.asarray(_sha256_file(args.config)),
        source_context_sha256=np.asarray(_sha256_file(args.source_context)),
        registry_sha256=np.asarray(
            _sha256_file(raw["module_tokenizer"]["registry_json_path"])
        ),
        activity_dictionary_sha256=np.asarray(
            _sha256_file(raw["module_tokenizer"]["activity_weight_path"])
        ),
        validation_donors_used=np.asarray(False),
        test_donors_used=np.asarray(False),
    )
    print(f"wrote {output} sha256={_sha256_file(output)}", flush=True)


def merge_shards(args: argparse.Namespace) -> None:
    if len(args.inputs) == 0:
        raise ValueError("merge requires at least one shard")
    opened = [np.load(path, allow_pickle=False) for path in args.inputs]
    try:
        indices = sorted(int(np.asarray(item["shard_index"]).item()) for item in opened)
        num_shards = int(np.asarray(opened[0]["num_shards"]).item())
        if indices != list(range(num_shards)) or len(opened) != num_shards:
            raise ValueError("input files do not form one complete shard partition")
        reference_fields = (
            "schema_version",
            "donor_names",
            "donor_split",
            "context_names",
            "module_names",
            "split_seed",
            "split_rule",
            "num_shards",
            "total_train_row_count",
            "config_sha256",
            "source_context_sha256",
            "registry_sha256",
            "activity_dictionary_sha256",
            "validation_donors_used",
            "test_donors_used",
        )
        first = opened[0]
        if str(np.asarray(first["schema_version"]).item()) != SCHEMA_SHARD:
            raise ValueError("unsupported split-half shard schema")
        for item in opened:
            for field in reference_fields:
                if not np.array_equal(item[field], first[field]):
                    raise ValueError(f"split-half shards differ at {field}")
            if bool(np.asarray(item["validation_donors_used"]).item()) or bool(
                np.asarray(item["test_donors_used"]).item()
            ):
                raise ValueError("a split-half shard used held-out donors")
        total_rows = sum(int(np.asarray(item["shard_row_count"]).item()) for item in opened)
        expected_rows = int(np.asarray(first["total_train_row_count"]).item())
        if total_rows != expected_rows:
            raise ValueError(f"merged row count {total_rows} != {expected_rows}")
        sum_a = sum((np.asarray(item["module_sum_half_a"]) for item in opened))
        sum_b = sum((np.asarray(item["module_sum_half_b"]) for item in opened))
        count_a = sum((np.asarray(item["cell_count_half_a"]) for item in opened))
        count_b = sum((np.asarray(item["cell_count_half_b"]) for item in opened))
        n_donor = len(first["donor_names"])
        n_context = len(first["context_names"])
        n_module = len(first["module_names"])
        expected_sum_shape = (n_donor * n_context, n_module)
        if sum_a.shape != expected_sum_shape or sum_b.shape != expected_sum_shape:
            raise ValueError("split-half additive array shape is invalid")
        mean_a = np.divide(
            sum_a,
            np.maximum(count_a[:, None], 1),
            out=np.zeros_like(sum_a),
        ).reshape(n_donor, n_context, n_module)
        mean_b = np.divide(
            sum_b,
            np.maximum(count_b[:, None], 1),
            out=np.zeros_like(sum_b),
        ).reshape(n_donor, n_context, n_module)
        observed_a = count_a.reshape(n_donor, n_context) >= int(
            args.minimum_cells_per_half
        )
        observed_b = count_b.reshape(n_donor, n_context) >= int(
            args.minimum_cells_per_half
        )
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output,
            schema_version=np.asarray(SCHEMA_MERGED),
            module_half_a=mean_a.astype(np.float32),
            module_half_b=mean_b.astype(np.float32),
            observed_half_a=observed_a,
            observed_half_b=observed_b,
            donor_names=np.asarray(first["donor_names"]).astype(str),
            donor_split=np.asarray(first["donor_split"]).astype(str),
            context_names=np.asarray(first["context_names"]).astype(str),
            module_names=np.asarray(first["module_names"]).astype(str),
            split_seed=np.asarray(first["split_seed"]),
            split_rule=np.asarray(first["split_rule"]),
            minimum_cells_per_half=np.asarray(
                int(args.minimum_cells_per_half), dtype=np.int64
            ),
            n_source_train_cells=np.asarray(total_rows, dtype=np.int64),
            config_sha256=np.asarray(first["config_sha256"]),
            source_context_sha256=np.asarray(first["source_context_sha256"]),
            registry_sha256=np.asarray(first["registry_sha256"]),
            activity_dictionary_sha256=np.asarray(
                first["activity_dictionary_sha256"]
            ),
            validation_donors_used=np.asarray(False),
            test_donors_used=np.asarray(False),
        )
        print(
            f"wrote {output} sha256={_sha256_file(output)} "
            f"train_cells={total_rows} observed_halves="
            f"{int((observed_a & observed_b).sum())}/{observed_a.size}",
            flush=True,
        )
    finally:
        for item in opened:
            item.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    shard = sub.add_parser("shard")
    shard.add_argument("--config", required=True)
    shard.add_argument("--source-context", required=True)
    shard.add_argument("--output", required=True)
    shard.add_argument("--shard-index", type=int, required=True)
    shard.add_argument("--num-shards", type=int, required=True)
    shard.add_argument("--seed", type=int, default=42)
    shard.add_argument("--batch-size", type=int, default=256)
    shard.add_argument("--device", default="cuda:0")
    merge = sub.add_parser("merge")
    merge.add_argument("--inputs", nargs="+", required=True)
    merge.add_argument("--output", required=True)
    merge.add_argument("--minimum-cells-per-half", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "shard":
        run_shard(args)
    else:
        merge_shards(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
