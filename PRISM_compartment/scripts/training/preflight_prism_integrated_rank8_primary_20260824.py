#!/usr/bin/env python3
"""Fail-closed parity, schedule, budget, sealed-test, and resume preflight."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import torch
from torch.utils.data import default_collate

from kmlee_bam.training import run_current
from kmlee_bam.training import runner_base as base
from kmlee_bam.training.integrated_phase_curriculum import (
    IntegratedPhaseCurriculumController,
)


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _build(path: Path, *, integrated: bool, train_ds=None):
    raw = _load(path)
    raw = copy.deepcopy(raw)
    raw.setdefault("train", {})["skip_final_test"] = True
    raw["train"]["reload_best_before_test"] = False
    # System/RNG parity does not depend on the criterion's 8,192-cell class
    # weight estimate.  Use the same deterministic 64-cell estimate in both
    # arms so this preflight tests model construction in seconds, not minutes.
    ordinal = raw.setdefault("v4", {}).setdefault("ordinal_balance", {})
    ordinal["max_cells_for_bin_stats"] = 64
    ordinal["bin_stats_log_every"] = 65
    run_current.set_settings(raw)
    cfg_path = path
    if raw != _load(path):
        cfg_path = Path("/tmp") / f"preflight_{path.name}"
        cfg_path.write_text(json.dumps(raw), encoding="utf-8")
    cfg = base.load_config(str(cfg_path))
    base.set_seed(int(cfg.train.seed))
    val_ds = None
    if train_ds is None:
        train_ds, val_ds, test_ds = base.build_datasets(cfg)
        if test_ds is not None:
            raise RuntimeError("sealed-test audit failed: test dataset was opened")
    system = base.build_system(cfg, train_ds)
    if integrated:
        system.precision_forward_enabled = False
        system.pathology_aux_forward_enabled = True
        if hasattr(system, "generator_count_epoch_state"):
            system.generator_count_epoch_state.fill_(1)
    return raw, cfg, train_ds, val_ds, system, torch.random.get_rng_state().clone()


def _tensor_digest(state):
    digest = hashlib.sha256()
    for key in sorted(state):
        value = state[key].detach().cpu().contiguous()
        digest.update(key.encode("utf-8"))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(bytes(value.numpy().tobytes()))
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase1", required=True, type=Path)
    parser.add_argument("--integrated", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    run_current.install_hooks()

    _, integrated_cfg, integrated_train, _, integrated_system, integrated_rng = _build(
        args.integrated, integrated=True
    )
    rescue_cfg = integrated_cfg.prism_module_rescue
    state_encoder = getattr(integrated_system, "state_encoder", None)
    persistent_pool = (
        getattr(state_encoder, "persistent_attention_pool", None)
        if state_encoder is not None
        else None
    )
    if bool(rescue_cfg.allow_state_readout_updates) and persistent_pool is None:
        raise RuntimeError(
            "CLS/non-persistent pooling cannot enable state-readout rescue updates"
        )
    _, p1_cfg, p1_train, _, p1_system, p1_rng = _build(
        args.phase1, integrated=False, train_ds=integrated_train
    )
    if len(p1_train) != len(integrated_train):
        raise RuntimeError("Phase-I dataset row count differs")
    if not torch.equal(p1_rng, integrated_rng):
        raise RuntimeError("dormant Phase-II construction advanced Phase-I RNG")

    phase1_state = p1_system.state_dict()
    integrated_state = integrated_system.state_dict()
    extra_prefixes = (
        "base_system.precision_head.",
        "base_system.generator_count_gate.",
        "base_system.generator_count_objective.",
        "base_system.generator_count_epoch_state",
        "precision_head.",
        "generator_count_gate.",
        "generator_count_objective.",
        "generator_count_epoch_state",
    )
    shared_integrated = {
        key: value
        for key, value in integrated_state.items()
        if not key.startswith(extra_prefixes)
    }
    if set(phase1_state) != set(shared_integrated):
        missing = sorted(set(phase1_state) - set(shared_integrated))
        unexpected = sorted(set(shared_integrated) - set(phase1_state))
        raise RuntimeError(
            f"shared Phase-I state key mismatch missing={missing[:5]} unexpected={unexpected[:5]}"
        )
    unequal = [
        key
        for key in phase1_state
        if not torch.equal(phase1_state[key], shared_integrated[key])
    ]
    if unequal:
        raise RuntimeError(f"shared Phase-I initialization mismatch: {unequal[:5]}")

    batch = default_collate([p1_train[index] for index in range(2)])
    p1_system.eval()
    integrated_system.eval()
    with torch.no_grad():
        p1_out = p1_system(batch, sample_latent=False)
        integrated_out = integrated_system(batch, sample_latent=False)
    output_diffs = {}
    output_tolerance = {
        "z_perp": 5.0e-4,
        "decoder_score": 2.0e-3,
        "decoder_probs": 5.0e-4,
    }
    for name, left, right in (
        ("z_perp", p1_out.z_perp, integrated_out.z_perp),
        ("decoder_score", p1_out.decoder_out.score, integrated_out.decoder_out.score),
        ("decoder_probs", p1_out.decoder_out.probs, integrated_out.decoder_out.probs),
    ):
        difference = float((left - right).abs().max())
        output_diffs[name] = difference
        if difference > output_tolerance[name]:
            raise RuntimeError(
                f"Phase-I forward exceeded tolerance for {name}: "
                f"{difference} > {output_tolerance[name]}"
            )

    controller = IntegratedPhaseCurriculumController(
        integrated_cfg.integrated_phase_curriculum
    )
    for epoch in range(1, 13):
        state = controller.resolve(epoch)
        expected = {
            "stage": "phase1_exact",
            "pathology_scale": 1.0,
            "learning_rate": 1.0e-4,
            "grad_accum_steps": 1,
            "batch_size": 64,
            "precision_forward_enabled": False,
            "lambda_state_abs": 0.08,
            "lambda_state_fraction": 0.01,
            "projection_aware_whitening": False,
            "phu_alignment_lambda": 0.001,
            "phu_alignment_start_epoch": 1,
            "phu_alignment_ramp_epochs": 5,
            "phu_relative_lambda": 0.0,
        }
        for key, target in expected.items():
            if state[key] != target:
                raise RuntimeError(f"epoch {epoch} Phase-I contract differs at {key}")
    boundary = {epoch: controller.resolve(epoch) for epoch in (12, 13, 16, 17)}
    if boundary[13]["stage"] != "phase1_to_phase2_crossfade":
        raise RuntimeError("epoch-13 boundary did not enter cross-fade")
    if boundary[17]["stage"] != "phase2_continuation":
        raise RuntimeError("epoch-17 boundary did not enter Phase-II continuation")

    world_size = 6
    per_rank = math.ceil(len(integrated_train) / world_size)
    phase1_batches = math.ceil(per_rank / 64)
    phase2_batches = math.ceil(per_rank / 32)
    phase1_updates = math.ceil(phase1_batches / 1)
    phase2_updates = math.ceil(phase2_batches / 2)
    if phase1_updates != phase2_updates:
        raise RuntimeError(
            f"optimizer-update budget differs: Phase-I={phase1_updates} Phase-II={phase2_updates}"
        )
    matched_budget_updates = 12 * phase1_updates + 20 * phase2_updates

    # Checkpointed controller state must restore deterministically at the
    # boundary and continue with the same epoch-13 resolution.
    source_controller = IntegratedPhaseCurriculumController(
        integrated_cfg.integrated_phase_curriculum
    )
    source_controller.epoch = 12
    source_controller.stage = "phase1_exact"
    source_controller.pathology_scale = 1.0
    source_controller.learning_rate = 1.0e-4
    source_controller.grad_accum_steps = 1
    source_controller.batch_size = 64
    restored_controller = IntegratedPhaseCurriculumController(
        integrated_cfg.integrated_phase_curriculum
    )
    restored_controller.load_state_dict(source_controller.state_dict())
    if restored_controller.resolve(13) != source_controller.resolve(13):
        raise RuntimeError("resume boundary controller resolution differs")

    report = {
        "schema_version": "kmlee_bam.prism_integrated_rank8_preflight.v1",
        "status": "PASS",
        "phase1_rows": len(integrated_train),
        "phase1_shared_state_tensor_count": len(phase1_state),
        "phase1_shared_state_sha256": _tensor_digest(phase1_state),
        "phase1_rng_after_construction_bit_identical": True,
        "phase1_forward_max_abs_difference": output_diffs,
        "phase1_forward_absolute_tolerance": output_tolerance,
        "phase1_epochs_contract_checked": [1, 12],
        "boundary_states": boundary,
        "world_size": world_size,
        "optimizer_updates_per_epoch_phase1": phase1_updates,
        "optimizer_updates_per_epoch_phase2": phase2_updates,
        "matched_separated_budget_checkpoint_epoch": 32,
        "matched_separated_budget_optimizer_updates": matched_budget_updates,
        "resume_boundary_state_roundtrip": True,
        "official_test_dataset_opened": False,
        "pathology_rank": int(integrated_cfg.decoder.pathology_rank),
        "personal_rank": int(integrated_cfg.precision_medicine.personal_rank),
        "state_readout_rescue_updates": bool(
            integrated_cfg.prism_module_rescue.allow_state_readout_updates
        ),
        "persistent_agp_pool_present": persistent_pool is not None,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
