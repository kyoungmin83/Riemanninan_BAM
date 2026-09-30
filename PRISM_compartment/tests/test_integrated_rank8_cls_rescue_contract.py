from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class IntegratedRank8ClsRescueContractTest(unittest.TestCase):
    def _load(self, name: str) -> dict:
        path = ROOT / "configs/final" / name
        return json.loads(path.read_text(encoding="utf-8"))

    def test_primary_cls_config_disables_agp_only_state_readout_rescue(self) -> None:
        config = self._load(
            "train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
        )
        self.assertEqual(config["encoder"]["pooling"], "cls")
        self.assertFalse(
            config["prism_module_rescue"]["allow_state_readout_updates"]
        )

    def test_epoch28_resume_is_strict_and_test_sealed(self) -> None:
        config = self._load(
            "train_config_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826.json"
        )
        train = config["train"]
        self.assertTrue(train["resume_checkpoint"].endswith("checkpoint_epoch_028.pt"))
        self.assertTrue(train["resume_optimizer"])
        self.assertTrue(train["resume_history"])
        source = self._load(
            "train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
        )
        self.assertNotEqual(train["out_dir"], source["train"]["out_dir"])
        self.assertTrue(train["skip_final_test"])
        self.assertTrue(config["_launch_guard"]["launch_allowed"])
        self.assertFalse(
            config["_launch_guard"]["explicit_user_reapproval_required"]
        )


if __name__ == "__main__":
    unittest.main()
