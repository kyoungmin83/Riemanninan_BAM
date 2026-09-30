from __future__ import annotations

import json
import os
import types
import unittest
from unittest.mock import patch

import torch

from kmlee_bam.training.architecture_capacity_logging import (
    build_architecture_capacity_record,
    format_architecture_capacity_record,
)
from kmlee_bam.training.core_trainer import RunningAverages, Trainer
from kmlee_bam.training.learned_generator_count import (
    HardBinaryConcreteGeneratorGate,
    LearnedGeneratorCountConfig,
)
from kmlee_bam.training.learned_pathology_rank import (
    LearnedPathologyRankConfig,
    LearnedPathologyRankGate,
)


def _rank_config() -> LearnedPathologyRankConfig:
    return LearnedPathologyRankConfig(
        enabled=True,
        initial_keep_probability=0.995,
        warmup_fixed_rank=8,
        warmup_end_epoch=24,
        initial_extra_keep_probability=0.5,
        soft_start_epoch=25,
        sparsity_start_epoch=29,
        sparsity_ramp_epochs=4,
        hard_start_epoch=33,
        freeze_epoch=37,
        require_stable_freeze=True,
    )


def _generator_config(*, start_epoch: int = 37) -> LearnedGeneratorCountConfig:
    return LearnedGeneratorCountConfig(
        enabled=True,
        mode="joint",
        start_epoch=start_epoch,
        shadow_end_epoch=max(start_epoch, 40),
        soft_end_epoch=max(start_epoch, 44),
        base_end_epoch=55,
        minimum_active_generators=None,
    )


class ArchitectureCapacityLoggingTest(unittest.TestCase):
    def test_snapshot_is_json_safe_and_does_not_mutate_gates(self) -> None:
        rank = LearnedPathologyRankGate(12, config=_rank_config())
        rank.set_epoch(30)
        generator_config = _generator_config()
        generator = HardBinaryConcreteGeneratorGate(
            10, config=generator_config
        )
        rank_before = {
            key: value.detach().clone() for key, value in rank.state_dict().items()
        }
        generator_before = {
            key: value.detach().clone()
            for key, value in generator.state_dict().items()
        }

        record = build_architecture_capacity_record(
            epoch=30,
            total_epochs=55,
            pathology_scale=1.0,
            rank_gate=rank,
            generator_gate=generator,
            generator_config=generator_config,
            metrics={
                "metric/v26_pathdec_absmean": 0.0123,
                "loss/total": 1.5,
                "weight/prism_module_local_nonlinear_ramp": 1.0,
                "metric/prism_module_local_module_rms": 0.02,
                "metric/prism_module_local_nonlinear_rms": 0.01,
                "metric/prism_module_local_nonlinear_to_local_ratio": 0.5,
                "metric/prism_module_local_threshold_crossing_fraction": 0.2,
                "metric/prism_module_local_nonlinear_branch_nll_gain": 0.003,
                "metric/prism_module_local_nonlinear_full_nll_gain": 0.002,
            },
            paper_config=types.SimpleNamespace(
                module_local_nonlinear_enabled=True,
                module_local_nonlinear_variant="compartmental_threshold",
            ),
        )
        json.dumps(record)
        self.assertEqual(record["pathology_rank"]["mode"], "soft_rank_learning")
        self.assertEqual(record["generator_count"]["mode"], "warmup_all_on")
        self.assertFalse(record["coordination_guard"]["simultaneous_cardinality"])
        self.assertEqual(
            record["paper_compartmental_branch"]["variant"],
            "compartmental_threshold",
        )
        for key, value in rank.state_dict().items():
            self.assertTrue(torch.equal(value, rank_before[key]))
        for key, value in generator.state_dict().items():
            self.assertTrue(torch.equal(value, generator_before[key]))

        rendered = format_architecture_capacity_record(record)
        self.assertIn("pathology rank", rendered)
        self.assertIn("role separation", rendered)
        self.assertIn("overlap=NO", rendered)
        self.assertIn("paper branch", rendered)
        self.assertIn("paper OFF audit", rendered)

    def test_overlap_guard_is_visible(self) -> None:
        rank = LearnedPathologyRankGate(12, config=_rank_config())
        rank.set_epoch(30)
        generator_config = _generator_config(start_epoch=30)
        generator = HardBinaryConcreteGeneratorGate(
            10, config=generator_config
        )
        record = build_architecture_capacity_record(
            epoch=30,
            total_epochs=55,
            pathology_scale=1.0,
            rank_gate=rank,
            generator_gate=generator,
            generator_config=generator_config,
            metrics={},
        )
        self.assertTrue(record["coordination_guard"]["simultaneous_cardinality"])
        self.assertEqual(record["coordination_guard"]["status"], "violation")

    def test_integrated_step_renderer_contains_live_capacity_block(self) -> None:
        trainer = Trainer.__new__(Trainer)
        trainer.amp_enabled = False
        trainer._prev_diag_metrics = {}
        trainer._verbose_diag = False
        trainer.current_epoch_index = 30
        trainer.integrated_curriculum_schedule = {
            "total_epochs": 55,
            "module_rescue_start": 25,
            "generator_start": 37,
        }
        trainer.system = types.SimpleNamespace(precision_head=None)
        meters = RunningAverages()
        metrics = {
            "loss/total": 1.5,
            "metric/pathology_rank_capacity": 96.0,
            "metric/pathology_rank_expected": 42.0,
            "metric/pathology_rank_hard": 40.0,
            "metric/pathology_rank_search_expected": 42.0,
            "metric/pathology_rank_search_hard": 40.0,
            "metric/pathology_rank_temperature": 0.8,
            "metric/pathology_rank_finalized": 0.0,
            "metric/pathology_rank_freeze_ready": 1.0,
            "metric/pathology_rank_freeze_count_low": 42.0,
            "metric/pathology_rank_freeze_count_mid": 40.0,
            "metric/pathology_rank_freeze_count_high": 39.0,
            "metric/pathology_rank_freeze_count_spread": 3.0,
            "metric/pathology_rank_freeze_near_threshold_uncertain_fraction": 0.03,
            "metric/pathology_rank_mode_soft_rank_learning": 1.0,
            "weight/pathology_rank_sparsity_multiplier": 0.5,
            "weight/pathology_curriculum_scale": 1.0,
            "metric/v26_pathdec_absmean": 0.0123,
            "metric/prism_paper_branch_enabled": 1.0,
            "metric/prism_paper_branch_variant_compartmental_threshold": 1.0,
            "metric/prism_module_local_module_rms": 0.02,
            "metric/prism_module_local_nonlinear_rms": 0.01,
            "metric/prism_module_local_nonlinear_to_local_ratio": 0.5,
            "metric/prism_module_local_threshold_crossing_fraction": 0.2,
            "metric/prism_module_local_nonlinear_branch_nll_gain": 0.003,
            "metric/prism_module_local_nonlinear_full_nll_gain": 0.002,
            "metric/prism_module_local_mix_mean": 0.1,
            "metric/prism_module_local_threshold_mean": 1.2,
            "metric/prism_module_local_slope_mean": 4.0,
            "metric/prism_module_local_gain_mean": 0.4,
            "weight/prism_module_local_nonlinear_ramp": 1.0,
            "metric/generator_warmup_all_on": 1.0,
            "metric/generator_hard_active": 414.0,
            "metric/generator_expected_active": 412.0,
            "metric/generator_candidate_count": 414.0,
        }
        meters.weighted_sums = dict(metrics)
        meters.weight_sums = {key: 1 for key in metrics}
        meters.n_examples = 1
        meters.n_finite_batches = 1
        with patch.dict(os.environ, {"KMLEE_CONSOLE_LOG_STYLE": "prism_integrated"}):
            rendered = trainer._format_step_log(
                prefix="train", step_idx=100, meters=meters
            )
        self.assertIn("구조 용량 자동학습", rendered)
        self.assertIn("Pathology rank live hard/E", rendered)
        self.assertIn("overlap NO", rendered)
        self.assertIn("ref 논문 적용부", rendered)
        self.assertIn("paper_compartmental_threshold", rendered)
        self.assertIn("이 경로를 끄면 생기는 ΔNLL", rendered)


if __name__ == "__main__":
    unittest.main()
