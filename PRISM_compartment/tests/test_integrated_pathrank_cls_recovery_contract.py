from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class IntegratedPathrankClsRecoveryContractTest(unittest.TestCase):
    def _load(self, name: str) -> dict:
        return json.loads(
            (ROOT / "configs/final" / name).read_text(encoding="utf-8")
        )

    def test_only_cls_agp_conflict_changes_outside_provenance(self) -> None:
        source = self._load(
            "train_config_prism_integrated_pathranklearn_s42_sv7_pre_recovery_20260826.json"
        )
        recovery = self._load(
            "train_config_prism_integrated_pathranklearn_s42_sv7_cls_recovery_20260826.json"
        )
        for config in (source, recovery):
            config.pop("experiment_manifest", None)
            config.pop("_launch_guard", None)
        source_out = source["train"].pop("out_dir")
        recovery_out = recovery["train"].pop("out_dir")
        source_resume = source["train"].pop("resume_checkpoint")
        recovery_resume = recovery["train"].pop("resume_checkpoint")
        self.assertNotEqual(source_out, recovery_out)
        self.assertIsNone(source_resume)
        self.assertTrue(recovery_resume.endswith("checkpoint_epoch_016.pt"))
        self.assertTrue(
            source["prism_module_rescue"].pop("allow_state_readout_updates")
        )
        self.assertFalse(
            recovery["prism_module_rescue"].pop("allow_state_readout_updates")
        )
        self.assertEqual(source, recovery)

    def test_integrated_rank_and_sealed_test_contracts_are_preserved(self) -> None:
        recovery = self._load(
            "train_config_prism_integrated_pathranklearn_s42_sv7_cls_recovery_20260826.json"
        )
        self.assertEqual(recovery["encoder"]["pooling"], "cls")
        self.assertTrue(recovery["learned_pathology_rank"]["enabled"])
        self.assertEqual(recovery["decoder"]["pathology_rank"], 0)
        self.assertEqual(recovery["learned_pathology_rank"]["freeze_epoch"], 17)
        self.assertTrue(recovery["train"]["skip_final_test"])
        self.assertTrue(
            recovery["train"]["resume_checkpoint"].endswith(
                "checkpoint_epoch_016.pt"
            )
        )
        self.assertTrue(recovery["_launch_guard"]["launch_allowed"])


if __name__ == "__main__":
    unittest.main()
