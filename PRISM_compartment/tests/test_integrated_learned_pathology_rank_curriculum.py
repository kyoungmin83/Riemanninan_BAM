from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

import torch

from kmlee_bam.model.lie_ordinal_decoder import LieActionOrdinalDecoder
from kmlee_bam.training.learned_pathology_rank import (
    LearnedPathologyRankConfig,
    LearnedPathologyRankGate,
)
from kmlee_bam.training.runner_base import load_config


ROOT = Path(__file__).resolve().parents[1]


def _integrated_rank_config() -> LearnedPathologyRankConfig:
    return LearnedPathologyRankConfig(
        enabled=True,
        initial_keep_probability=0.995,
        warmup_fixed_rank=8,
        warmup_end_epoch=24,
        initial_extra_keep_probability=0.5,
        require_canonical_phase1_rank8=True,
        require_active_pathology_route_during_search=True,
        temperature_start=2.0,
        temperature_end=0.35,
        soft_start_epoch=25,
        hard_start_epoch=33,
        freeze_epoch=37,
        sparsity_start_epoch=29,
        sparsity_ramp_epochs=4,
        cardinality_loss_fraction=0.0025,
        learning_rate=3.0e-5,
        weight_decay=0.0,
        hard_threshold=0.5,
        require_stable_freeze=True,
        freeze_threshold_low=0.45,
        freeze_threshold_high=0.55,
        maximum_freeze_count_spread=4,
        maximum_freeze_uncertain_fraction=0.15,
    )


class IntegratedLearnedPathologyRankCurriculumTest(unittest.TestCase):
    def test_expanded_rank_preserves_canonical_rng_trajectory(self) -> None:
        common = dict(
            n_genes=17,
            n_celltypes=3,
            n_tech=2,
            d_z=5,
            n_bins=4,
            n_generators=4,
            use_pathology_decoder=True,
            n_pathology_axes=4,
            pathology_pairwise=False,
        )
        torch.manual_seed(420828)
        canonical = LieActionOrdinalDecoder(**common, pathology_rank=8)
        canonical_rng = torch.random.get_rng_state().clone()
        torch.manual_seed(420828)
        expanded = LieActionOrdinalDecoder(
            **common,
            pathology_rank=12,
            pathology_rank_rng_compat_prefix=8,
        )
        expanded_rng = torch.random.get_rng_state().clone()
        self.assertTrue(torch.equal(canonical.path_V, expanded.path_V[:8]))
        self.assertTrue(torch.equal(canonical_rng, expanded_rng))
        for name in (
            "generator_u",
            "generator_v",
            "generator_a",
            "raw_threshold_start",
            "raw_threshold_deltas",
        ):
            self.assertTrue(
                torch.equal(getattr(canonical, name), getattr(expanded, name)),
                msg=name,
            )

    def test_fixed_rank8_warmup_is_exact_and_has_zero_gate_gradient(self) -> None:
        gate = LearnedPathologyRankGate(96, config=_integrated_rank_config())
        gate.set_epoch(12)
        mask = gate()
        self.assertEqual(gate.mode(), "fixed_rank_warmup")
        self.assertTrue(torch.equal(mask[:8], torch.ones(8)))
        self.assertTrue(torch.equal(mask[8:], torch.zeros(88)))
        mask.sum().backward()
        self.assertIsNotNone(gate.log_alpha.grad)
        self.assertTrue(torch.equal(gate.log_alpha.grad, torch.zeros(96)))
        diagnostics = gate.diagnostics()
        self.assertEqual(diagnostics["hard_rank"], 8)
        self.assertEqual(diagnostics["expected_rank"], 8.0)
        self.assertEqual(diagnostics["search_hard_rank"], 96)

    def test_search_opens_after_warmup_and_freezes_checkpointed_mask(self) -> None:
        config = _integrated_rank_config()
        gate = LearnedPathologyRankGate(12, config=config)
        gate.set_epoch(25)
        expected = torch.cat(
            (torch.full((8,), 0.995), torch.full((4,), 0.5))
        )
        self.assertEqual(gate.mode(), "soft_rank_learning")
        self.assertTrue(torch.allclose(gate(), expected, atol=1.0e-6, rtol=0.0))

        with torch.no_grad():
            gate.log_alpha[:5].fill_(2.0)
            gate.log_alpha[5:].fill_(-2.0)
        gate.set_epoch(37)
        self.assertEqual(gate.mode(), "frozen_hard")
        self.assertEqual(gate.diagnostics()["hard_rank"], 5)

        restored = LearnedPathologyRankGate(12, config=config)
        restored.load_state_dict(gate.state_dict())
        self.assertEqual(restored.mode(), "frozen_hard")
        self.assertTrue(torch.equal(restored(), gate()))

    def test_warmup_rank_cannot_exceed_capacity(self) -> None:
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            LearnedPathologyRankGate(7, config=_integrated_rank_config())

    def test_threshold_sensitive_mask_is_not_frozen(self) -> None:
        gate = LearnedPathologyRankGate(12, config=_integrated_rank_config())
        gate.set_epoch(36)
        with self.assertRaisesRegex(RuntimeError, "freeze rejected"):
            gate.set_epoch(37)
        self.assertFalse(gate.diagnostics()["finalized"])

    def _runtime_json(self) -> dict:
        payload = json.loads(
            (
                ROOT
                / "configs/final/train_config_prism_integrated_pathranklearn_s42_sv7_active_source_20260826.json"
            ).read_text(encoding="utf-8")
        )
        payload["learned_pathology_rank"].update(
            {
                "warmup_fixed_rank": 8,
                "warmup_end_epoch": 24,
                "initial_extra_keep_probability": 0.5,
                "require_canonical_phase1_rank8": True,
                "require_active_pathology_route_during_search": True,
                "soft_start_epoch": 25,
                "hard_start_epoch": 33,
                "freeze_epoch": 37,
                "sparsity_start_epoch": 29,
                "sparsity_ramp_epochs": 4,
                "require_stable_freeze": True,
                "freeze_threshold_low": 0.45,
                "freeze_threshold_high": 0.55,
                "maximum_freeze_count_spread": 4,
                "maximum_freeze_uncertain_fraction": 0.15,
            }
        )
        payload["learned_generator_count"].update(
            {
                "start_epoch": 37,
                "shadow_end_epoch": 40,
                "soft_end_epoch": 44,
                "require_frozen_pathology_rank_before_search": True,
            }
        )
        payload["integrated_phase_curriculum"]["phase2_pathology_scale"] = 1.0
        return payload

    def _load_temporary(self, payload: dict):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            return load_config(str(path))

    def test_cross_config_contract_accepts_sequential_rank_then_generator(self) -> None:
        config = self._load_temporary(self._runtime_json())
        self.assertEqual(config.learned_pathology_rank.warmup_fixed_rank, 8)
        self.assertEqual(config.learned_generator_count.start_epoch, 37)
        self.assertEqual(config.integrated_phase_curriculum.phase2_pathology_scale, 1.0)

    def test_cross_config_contract_rejects_dormant_pathology_search(self) -> None:
        payload = self._runtime_json()
        payload["integrated_phase_curriculum"]["phase2_pathology_scale"] = 0.0
        with self.assertRaisesRegex(ValueError, "phase2_pathology_route_active"):
            self._load_temporary(payload)

    def test_cross_config_contract_rejects_generator_overlap(self) -> None:
        payload = self._runtime_json()
        payload["learned_generator_count"]["start_epoch"] = 36
        with self.assertRaisesRegex(ValueError, "cannot start before"):
            self._load_temporary(payload)

    def test_historical_rank_config_remains_loadable(self) -> None:
        payload = self._runtime_json()
        historical = copy.deepcopy(payload)
        for key in (
            "warmup_fixed_rank",
            "warmup_end_epoch",
            "initial_extra_keep_probability",
            "require_canonical_phase1_rank8",
            "require_active_pathology_route_during_search",
        ):
            historical["learned_pathology_rank"].pop(key, None)
        historical["learned_pathology_rank"].update(
            {
                "soft_start_epoch": 1,
                "hard_start_epoch": 13,
                "freeze_epoch": 17,
                "sparsity_start_epoch": 3,
                "sparsity_ramp_epochs": 6,
            }
        )
        historical["learned_generator_count"].pop(
            "require_frozen_pathology_rank_before_search", None
        )
        historical["learned_generator_count"].update(
            {"start_epoch": 13, "shadow_end_epoch": 16, "soft_end_epoch": 24}
        )
        historical["integrated_phase_curriculum"]["phase2_pathology_scale"] = 0.0
        loaded = self._load_temporary(historical)
        self.assertEqual(loaded.learned_pathology_rank.warmup_fixed_rank, 0)


if __name__ == "__main__":
    unittest.main()
