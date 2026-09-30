"""Checkpointed Phase-I -> Phase-II controls for one continuous PRISM run.

The controller intentionally changes only runtime contribution scales, learning
rate, accumulation and loader selection.  All model parameters are constructed
before epoch one and remain in one optimizer/checkpoint lineage.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any

import torch


@dataclass
class IntegratedPhaseCurriculumConfig:
    enabled: bool = False
    phase1_end_epoch: int = 12
    phase2_start_epoch: int = 13
    pathology_crossfade_end_epoch: int = 16
    phase1_learning_rate: float = 1.0e-4
    phase2_learning_rate: float = 2.0e-5
    phase1_grad_accum_steps: int = 1
    phase2_grad_accum_steps: int = 2
    phase1_batch_size: int = 64
    phase2_batch_size: int = 32
    phase1_pathology_scale: float = 1.0
    phase2_pathology_scale: float = 0.0
    phase1_lambda_state_abs: float = 0.08
    phase1_lambda_state_fraction: float = 0.01
    phase2_lambda_state_abs: float = 0.0
    phase2_lambda_state_fraction: float = 0.0
    phase1_projection_aware_whitening: bool = False
    phase2_projection_aware_whitening: bool = True
    phase1_phu_alignment_lambda: float = 1.0e-3
    phase1_phu_alignment_start_epoch: int = 1
    phase1_phu_alignment_ramp_epochs: int = 5
    phase2_phu_alignment_lambda: float = 5.0e-3
    phase2_phu_alignment_start_epoch: int = 17
    phase2_phu_alignment_ramp_epochs: int = 4
    phase1_phu_relative_lambda: float = 0.0
    phase2_phu_relative_lambda: float = 0.5
    phase2_phu_relative_start_epoch: int = 25
    phase2_phu_relative_ramp_epochs: int = 5
    skip_dormant_precision_forward: bool = True
    require_exact_phase1_contract: bool = True

    def validate(self) -> None:
        if not bool(self.enabled):
            return
        if int(self.phase1_end_epoch) < 1:
            raise ValueError("phase1_end_epoch must be positive")
        if int(self.phase2_start_epoch) != int(self.phase1_end_epoch) + 1:
            raise ValueError("phase2_start_epoch must immediately follow phase1_end_epoch")
        if int(self.pathology_crossfade_end_epoch) < int(self.phase2_start_epoch):
            raise ValueError("pathology cross-fade cannot end before Phase II starts")
        for name in (
            "phase1_learning_rate",
            "phase2_learning_rate",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for name in (
            "phase1_grad_accum_steps",
            "phase2_grad_accum_steps",
            "phase1_batch_size",
            "phase2_batch_size",
        ):
            if int(getattr(self, name)) <= 0:
                raise ValueError(f"{name} must be positive")
        for name in ("phase1_pathology_scale", "phase2_pathology_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")


class IntegratedPhaseCurriculumController:
    """Resolve and install the deterministic global-epoch curriculum state."""

    schema_version = "kmlee_bam.integrated_phase_curriculum.v2"

    def __init__(self, config: IntegratedPhaseCurriculumConfig) -> None:
        config.validate()
        self.config = config
        self.epoch = 0
        self.stage = "uninitialized"
        self.pathology_scale = float(config.phase1_pathology_scale)
        self.learning_rate = float(config.phase1_learning_rate)
        self.grad_accum_steps = int(config.phase1_grad_accum_steps)
        self.batch_size = int(config.phase1_batch_size)

    @staticmethod
    def _unwrap(system: Any) -> Any:
        return getattr(system, "module", system)

    def pathology_scale_at_epoch(self, epoch: int) -> float:
        cfg = self.config
        epoch = int(epoch)
        if epoch <= int(cfg.phase1_end_epoch):
            return float(cfg.phase1_pathology_scale)
        start = int(cfg.phase2_start_epoch)
        end = int(cfg.pathology_crossfade_end_epoch)
        if epoch > end:
            return float(cfg.phase2_pathology_scale)
        # The first Phase-II epoch takes the first cross-fade step and the end
        # epoch reaches the exact Phase-II scale.
        fraction = (epoch - start + 1) / float(end - start + 1)
        return float(
            cfg.phase1_pathology_scale
            + fraction
            * (cfg.phase2_pathology_scale - cfg.phase1_pathology_scale)
        )

    def resolve(self, epoch: int) -> dict[str, Any]:
        cfg = self.config
        phase1 = int(epoch) <= int(cfg.phase1_end_epoch)
        pathology_scale = self.pathology_scale_at_epoch(int(epoch))
        if phase1:
            stage = "phase1_exact"
        elif int(epoch) <= int(cfg.pathology_crossfade_end_epoch):
            stage = "phase1_to_phase2_crossfade"
        else:
            stage = "phase2_continuation"
        return {
            "epoch": int(epoch),
            "stage": stage,
            "pathology_scale": float(pathology_scale),
            "learning_rate": float(
                cfg.phase1_learning_rate if phase1 else cfg.phase2_learning_rate
            ),
            "grad_accum_steps": int(
                cfg.phase1_grad_accum_steps if phase1 else cfg.phase2_grad_accum_steps
            ),
            "batch_size": int(
                cfg.phase1_batch_size if phase1 else cfg.phase2_batch_size
            ),
            "precision_forward_enabled": bool(not phase1),
            "lambda_state_abs": float(
                cfg.phase1_lambda_state_abs if phase1 else cfg.phase2_lambda_state_abs
            ),
            "lambda_state_fraction": float(
                cfg.phase1_lambda_state_fraction
                if phase1
                else cfg.phase2_lambda_state_fraction
            ),
            "projection_aware_whitening": bool(
                cfg.phase1_projection_aware_whitening
                if phase1
                else cfg.phase2_projection_aware_whitening
            ),
            "phu_alignment_lambda": float(
                cfg.phase1_phu_alignment_lambda
                if phase1
                else cfg.phase2_phu_alignment_lambda
            ),
            "phu_alignment_start_epoch": int(
                cfg.phase1_phu_alignment_start_epoch
                if phase1
                else cfg.phase2_phu_alignment_start_epoch
            ),
            "phu_alignment_ramp_epochs": int(
                cfg.phase1_phu_alignment_ramp_epochs
                if phase1
                else cfg.phase2_phu_alignment_ramp_epochs
            ),
            "phu_relative_lambda": float(
                cfg.phase1_phu_relative_lambda
                if phase1
                else cfg.phase2_phu_relative_lambda
            ),
            "phu_relative_start_epoch": int(
                cfg.phase2_phu_relative_start_epoch
            ),
            "phu_relative_ramp_epochs": int(
                cfg.phase2_phu_relative_ramp_epochs
            ),
        }

    @torch.no_grad()
    def apply_epoch(self, trainer: Any, epoch: int) -> dict[str, Any]:
        state = self.resolve(int(epoch))
        system = self._unwrap(trainer.system)
        decoder = getattr(system, "decoder", None)
        if decoder is not None and hasattr(decoder, "pathology_curriculum_scale"):
            decoder.pathology_curriculum_scale.fill_(state["pathology_scale"])
        system.precision_forward_enabled = bool(state["precision_forward_enabled"])
        system.pathology_aux_forward_enabled = bool(
            float(state["pathology_scale"]) > 0.0
        )
        trainer.pathology_aux_curriculum_scale = float(state["pathology_scale"])
        trainer.grad_accum_steps = int(state["grad_accum_steps"])
        state_usage = getattr(trainer, "v6_state_usage_config", None)
        if state_usage is not None:
            # V5StateUsageConfig is intentionally frozen.  Replace the
            # trainer-owned value atomically instead of mutating the immutable
            # dataclass at the first epoch boundary.
            trainer.v6_state_usage_config = replace(
                state_usage,
                lambda_state_abs=float(state["lambda_state_abs"]),
                lambda_state_fraction=float(state["lambda_state_fraction"]),
            )
        nuisance = getattr(system, "nuisance_projector", None)
        nuisance_config = getattr(nuisance, "cfg", None)
        if nuisance_config is not None:
            nuisance_config.projection_aware_whitening = bool(
                state["projection_aware_whitening"]
            )
        phu = getattr(trainer, "v7a_alignment_config", None)
        if phu is not None:
            phu.lambda_alignment = float(state["phu_alignment_lambda"])
            phu.alignment_start_epoch = int(
                state["phu_alignment_start_epoch"]
            )
            phu.ramp_epochs = int(state["phu_alignment_ramp_epochs"])
            phu.lambda_relative_rec = float(state["phu_relative_lambda"])
            phu.relative_start_epoch = int(
                state["phu_relative_start_epoch"]
            )
            phu.relative_ramp_epochs = int(
                state["phu_relative_ramp_epochs"]
            )

        # Preserve relative LR multipliers (precision and gate groups) while
        # changing the backbone regime at the Phase boundary.
        previous_base = float(self.learning_rate)
        next_base = float(state["learning_rate"])
        ratio = next_base / previous_base if previous_base > 0.0 else 1.0
        if self.epoch == 0:
            configured_initial = float(self.config.phase1_learning_rate)
            for group in trainer.optimizer.param_groups:
                name = str(group.get("name", ""))
                if name in {"joint_generator_gate", "pathology_rank_gate"}:
                    continue
                if name == "precision_new":
                    # The optimizer builder has already applied its declared
                    # multiplier to the Phase-I base LR.
                    continue
                group["lr"] = configured_initial
        elif abs(ratio - 1.0) > 1.0e-15:
            for group in trainer.optimizer.param_groups:
                if str(group.get("name", "")) in {
                    "joint_generator_gate",
                    "pathology_rank_gate",
                }:
                    continue
                group["lr"] = float(group["lr"]) * ratio

        self.epoch = int(state["epoch"])
        self.stage = str(state["stage"])
        self.pathology_scale = float(state["pathology_scale"])
        self.learning_rate = next_base
        self.grad_accum_steps = int(state["grad_accum_steps"])
        self.batch_size = int(state["batch_size"])
        return state

    def state_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "epoch": int(self.epoch),
            "stage": str(self.stage),
            "pathology_scale": float(self.pathology_scale),
            "learning_rate": float(self.learning_rate),
            "grad_accum_steps": int(self.grad_accum_steps),
            "batch_size": int(self.batch_size),
        }

    def load_state_dict(self, state: dict[str, Any]) -> None:
        if not state:
            return
        if state.get("schema_version") != self.schema_version:
            raise ValueError("integrated Phase curriculum state schema mismatch")
        expected = self.resolve(int(state["epoch"]))
        for key in ("stage", "pathology_scale", "learning_rate", "grad_accum_steps", "batch_size"):
            observed = state[key]
            target = expected[key]
            if isinstance(target, float):
                if abs(float(observed) - float(target)) > 1.0e-12:
                    raise ValueError(f"non-deterministic restored curriculum field: {key}")
            elif observed != target:
                raise ValueError(f"non-deterministic restored curriculum field: {key}")
        self.epoch = int(state["epoch"])
        self.stage = str(state["stage"])
        self.pathology_scale = float(state["pathology_scale"])
        self.learning_rate = float(state["learning_rate"])
        self.grad_accum_steps = int(state["grad_accum_steps"])
        self.batch_size = int(state["batch_size"])


def pathology_correction_rank_diagnostics(system: Any) -> dict[str, Any]:
    """Return optimizer-neutral realized-rank diagnostics for each named axis."""

    base = getattr(system, "module", system)
    decoder = getattr(base, "decoder", None)
    if decoder is None or not bool(getattr(decoder, "use_pathology_decoder", False)):
        raise RuntimeError("rank diagnostics require the pathology decoder")
    U = decoder.path_single_U.detach().float().cpu()
    V = decoder.path_V.detach().float().cpu()
    rank_gate = getattr(decoder, "pathology_rank_gate", None)
    gate_report = None
    if rank_gate is not None:
        gate_report = rank_gate.diagnostics()
        gate_mask = rank_gate().detach().float().cpu()
        U = U * gate_mask.view(1, 1, -1)
    rows = []
    for axis in range(int(U.shape[1])):
        correction = U[:, axis, :] @ V
        singular = torch.linalg.svdvals(correction)
        energy = singular.square()
        total = energy.sum()
        if float(total) <= 0.0:
            participation = 0.0
            rank95 = 0
        else:
            participation = float(total.square() / energy.square().sum().clamp_min(1.0e-30))
            rank95 = int(
                torch.searchsorted(
                    torch.cumsum(energy, dim=0) / total,
                    torch.tensor(0.95),
                ).item()
                + 1
            )
        rows.append(
            {
                "axis_index": axis,
                "participation_rank": participation,
                "rank_95_energy": rank95,
                "singular_values": [float(value) for value in singular.tolist()],
            }
        )
    report = {
        "pathology_rank_capacity": int(V.shape[0]),
        "rank_selection": "learned" if rank_gate is not None else "fixed",
        "axes": rows,
    }
    if gate_report is not None:
        report["learned_rank"] = gate_report
    return report
