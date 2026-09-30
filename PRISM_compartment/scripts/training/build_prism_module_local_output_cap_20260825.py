#!/usr/bin/env python3
"""Calibrate the module-local hard cap from a sealed frozen rank-2 model.

Only the checkpoint's target-celltype-excluded personal path is evaluated.
The absolute coefficient quantile is computed across the exact training-donor
allow-list and every target cell type; validation and test donors are never
included.  The resulting small NPZ is fail-closed by the runtime loader.
"""

from __future__ import annotations

import argparse
from dataclasses import fields, replace
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
import torch

from kmlee_bam.data.donor_context_dataset import DonorContextTable
from kmlee_bam.data.module_local_reliability import (
    MODULE_LOCAL_OUTPUT_CAP_SCHEMA,
    sha256_file,
    sha256_string_sequence,
)
from kmlee_bam.model.precision_medicine import (
    PrecisionMedicineConfig,
    PrecisionMedicineHead,
)


COEFFICIENT_ORIGIN = (
    "frozen_rank2_train_donor_personal_coeff_absolute_quantile"
)


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value.lower()
    )


def _precision_config(path: Path) -> PrecisionMedicineConfig:
    raw = json.loads(path.read_text(encoding="utf-8"))
    section = raw.get("precision_medicine")
    if not isinstance(section, dict):
        raise ValueError("source config lacks a precision_medicine object")
    accepted = {field.name for field in fields(PrecisionMedicineConfig)}
    config = PrecisionMedicineConfig(
        **{key: value for key, value in section.items() if key in accepted}
    )
    if not config.context_npz_path:
        raise ValueError("source precision config lacks context_npz_path")
    return replace(
        config,
        enabled=True,
        module_local_enabled=False,
        module_local_reliability_npz_path=None,
        module_local_output_cap_npz_path=None,
        module_local_output_cap_source_checkpoint_sha256=None,
        module_local_output_cap_source_config_sha256=None,
        interaction_enabled=False,
        lambda_support_query_infonce=0.0,
    )


def _precision_state(payload: object) -> dict[str, torch.Tensor]:
    state = (
        payload["system_state_dict"]
        if isinstance(payload, dict) and "system_state_dict" in payload
        else payload
    )
    if not isinstance(state, dict):
        raise ValueError("checkpoint does not contain a state dictionary")
    prefixes = ("precision_head.", "module.precision_head.")
    for prefix in prefixes:
        selected = {
            str(key)[len(prefix) :]: value
            for key, value in state.items()
            if str(key).startswith(prefix) and isinstance(value, torch.Tensor)
        }
        if selected:
            return selected
    raise ValueError("checkpoint contains no precision_head state")


def _build_frozen_head(
    config: PrecisionMedicineConfig,
    state: dict[str, torch.Tensor],
) -> PrecisionMedicineHead:
    required_buffers = (
        "source_module",
        "source_latent",
        "source_observed",
        "source_reliability",
        "context_celltype",
        "context_region",
        "personal_basis",
        "normal_region_delta",
    )
    missing_buffers = [key for key in required_buffers if key not in state]
    if missing_buffers:
        raise ValueError(
            f"checkpoint lacks output-cap calibration state: {missing_buffers}"
        )
    personal_basis = state["personal_basis"]
    if personal_basis.ndim != 3 or int(personal_basis.shape[1]) != 2:
        raise ValueError("output cap must come from a personal-rank2 checkpoint")
    normal_region_delta = state["normal_region_delta"]
    if normal_region_delta.ndim != 3:
        raise ValueError("normal_region_delta has an unexpected shape")
    head = PrecisionMedicineHead(
        config=replace(config, personal_rank=2),
        source_module=state["source_module"],
        source_latent=state["source_latent"],
        source_observed=state["source_observed"],
        source_reliability=state["source_reliability"],
        context_celltype=state["context_celltype"],
        context_region=state["context_region"],
        n_celltypes=int(personal_basis.shape[0]),
        n_regions=int(normal_region_delta.shape[1]),
    )
    target_state = head.state_dict()
    compatible = {
        key: value
        for key, value in state.items()
        if key in target_state and tuple(value.shape) == tuple(target_state[key].shape)
    }
    required_prefixes = (
        "source_",
        "context_celltype",
        "context_region",
        "module_projection.",
        "latent_projection.",
        "context_celltype_embedding.",
        "context_region_embedding.",
        "context_token_projection.",
        "context_encoder.",
        "personal_posterior.",
        "personal_basis",
    )
    required_keys = {
        key
        for key in target_state
        if key.startswith(required_prefixes)
    }
    missing_required = sorted(required_keys.difference(compatible))
    if missing_required:
        raise ValueError(
            "checkpoint cannot reproduce the frozen personal path; missing "
            f"keys={missing_required}"
        )
    head.load_state_dict(compatible, strict=False)
    head.eval()
    for parameter in head.parameters():
        parameter.requires_grad_(False)
    return head


@torch.no_grad()
def _absolute_personal_coefficients(
    head: PrecisionMedicineHead,
    train_donor_ids: np.ndarray,
    *,
    batch_size: int,
) -> np.ndarray:
    donor = torch.from_numpy(
        np.repeat(train_donor_ids.astype(np.int64), head.n_celltypes)
    )
    celltype = torch.arange(head.n_celltypes, dtype=torch.long).repeat(
        len(train_donor_ids)
    )
    values: list[np.ndarray] = []
    for start in range(0, int(donor.numel()), int(batch_size)):
        stop = min(start + int(batch_size), int(donor.numel()))
        donor_batch = donor[start:stop]
        celltype_batch = celltype[start:stop]
        code, _, _ = head.infer_personal_code(donor_batch, celltype_batch)
        basis = head.personal_basis[celltype_batch]
        coefficient = torch.einsum("bp,bpm->bm", code, basis)
        values.append(coefficient.abs().cpu().numpy().reshape(-1))
    result = np.concatenate(values).astype(np.float64, copy=False)
    if result.size == 0 or not np.isfinite(result).all():
        raise ValueError("frozen personal coefficients are empty or non-finite")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--absolute-quantile", type=float, default=0.995)
    parser.add_argument("--expected-train-donors", type=int, default=64)
    parser.add_argument("--batch-size", type=int, default=128)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 0.5 < float(args.absolute_quantile) < 1.0:
        raise ValueError("absolute quantile must lie in (0.5,1)")
    if int(args.expected_train_donors) <= 0 or int(args.batch_size) <= 0:
        raise ValueError("donor count and batch size must be positive")
    expected_checkpoint_sha256 = str(args.expected_checkpoint_sha256).lower()
    if not _is_sha256(expected_checkpoint_sha256):
        raise ValueError("expected checkpoint SHA-256 is invalid")

    config_path = Path(args.config).expanduser().resolve()
    checkpoint_path = Path(args.checkpoint).expanduser().resolve()
    observed_checkpoint_sha256 = sha256_file(checkpoint_path)
    if observed_checkpoint_sha256 != expected_checkpoint_sha256:
        raise RuntimeError(
            "frozen comparator checkpoint SHA-256 mismatch: "
            f"{observed_checkpoint_sha256} != {expected_checkpoint_sha256}"
        )
    config = _precision_config(config_path)
    context_path = Path(config.context_npz_path).expanduser().resolve()
    table = DonorContextTable.from_npz(
        context_path,
        min_cells=int(config.min_context_cells),
    )
    train_donor_ids = np.flatnonzero(
        table.donor_split.astype(str) == "train"
    )
    if int(train_donor_ids.size) != int(args.expected_train_donors):
        raise ValueError(
            "frozen comparator training-donor count mismatch: "
            f"{train_donor_ids.size} != {args.expected_train_donors}"
        )
    train_donor_names = tuple(
        str(value) for value in table.donor_names[train_donor_ids].tolist()
    )

    payload = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )
    state = _precision_state(payload)
    if int(state["source_module"].shape[0]) != int(len(table.donor_names)):
        raise ValueError("checkpoint context donor count differs from source table")
    head = _build_frozen_head(config, state)
    absolute_coefficients = _absolute_personal_coefficients(
        head,
        train_donor_ids,
        batch_size=int(args.batch_size),
    )
    output_cap = float(
        np.quantile(absolute_coefficients, float(args.absolute_quantile))
    )
    if not np.isfinite(output_cap) or output_cap <= 0.0:
        raise RuntimeError("calibrated module-local output cap is invalid")

    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    artifact = {
        "schema_version": np.asarray(MODULE_LOCAL_OUTPUT_CAP_SCHEMA),
        "output_cap": np.asarray(output_cap, dtype=np.float64),
        "absolute_quantile": np.asarray(
            float(args.absolute_quantile), dtype=np.float64
        ),
        "personal_rank": np.asarray(2, dtype=np.int64),
        "personal_coefficient_count": np.asarray(
            int(absolute_coefficients.size), dtype=np.int64
        ),
        "train_donor_allowlist": np.asarray(train_donor_names, dtype=str),
        "train_donor_allowlist_sha256": np.asarray(
            sha256_string_sequence(train_donor_names)
        ),
        "checkpoint_sha256": np.asarray(observed_checkpoint_sha256),
        "source_config_sha256": np.asarray(sha256_file(config_path)),
        "source_context_sha256": np.asarray(sha256_file(context_path)),
        "coefficient_origin": np.asarray(COEFFICIENT_ORIGIN),
        "validation_donors_used": np.asarray(False),
        "test_donors_used": np.asarray(False),
    }
    with tempfile.NamedTemporaryFile(
        prefix=output_path.name + ".",
        suffix=".tmp.npz",
        dir=output_path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
    try:
        np.savez_compressed(temporary_path, **artifact)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()
    print(
        f"wrote {output_path} sha256={sha256_file(output_path)} "
        f"output_cap={output_cap:.8g} quantile={args.absolute_quantile} "
        f"train_donors={len(train_donor_names)} "
        f"coefficients={absolute_coefficients.size}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FATAL: {error}", file=sys.stderr, flush=True)
        raise
