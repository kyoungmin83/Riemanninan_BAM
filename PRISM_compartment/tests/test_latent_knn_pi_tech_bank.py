from __future__ import annotations

import unittest
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset

from kmlee_bam.objectives.latent_knn_pi_tech import (
    LatentKNNPiTechConfig,
    TrainOnlyPiTechBank,
    compute_pi_tech_from_bank,
    compute_pi_tech_in_batch,
    select_donor_celltype_bank_indices,
)
from kmlee_bam.training.adaptive_subgroup_trainer import V8Trainer
from kmlee_bam.training.run_current import (
    _resolve_pi_tech_sex_linked_gene_indices,
)


def _two_state_bank() -> TrainOnlyPiTechBank:
    z = []
    detect = []
    donor = []
    sex = []
    row = []
    # Every donor contributes one ON-like and one OFF-like cell.  A nearest
    # donor vote must choose only one of them, preventing within-donor size from
    # becoming pseudoreplication.
    for donor_id in range(6):
        z.extend(([0.01 * donor_id, 0.0], [10.0 + 0.01 * donor_id, 0.0]))
        detect.extend(([1, 1], [0, 0]))
        donor.extend((donor_id, donor_id))
        sex.extend((1, 0))
        row.extend((2 * donor_id, 2 * donor_id + 1))
    # A closer wrong-cell-type cell must never enter the neighbor vote.
    z.append([0.0, 0.0])
    detect.append([0, 0])
    donor.append(6)
    sex.append(0)
    row.append(100)
    return TrainOnlyPiTechBank(
        z=torch.tensor(z, dtype=torch.float32),
        detect=torch.tensor(detect, dtype=torch.uint8),
        row=torch.tensor(row, dtype=torch.long),
        celltype=torch.tensor([0] * 12 + [1], dtype=torch.long),
        region=torch.zeros(13, dtype=torch.long),
        donor=torch.tensor(donor, dtype=torch.long),
        sex=torch.tensor(sex, dtype=torch.long),
        epoch=13,
    )


class TrainOnlyPiTechBankTest(unittest.TestCase):
    def _config(self, **overrides) -> LatentKNNPiTechConfig:
        values = dict(
            enabled=True,
            mode="same_celltype_donor_bank",
            k=4,
            bank_min_distinct_donors=4,
            bank_exclude_query_donor=True,
            pi_cap=0.7,
            sex_linked_gene_indices=(1,),
        )
        values.update(overrides)
        return LatentKNNPiTechConfig(**values)

    def test_same_celltype_distinct_donor_neighbors_separate_proven_and_stable(self) -> None:
        bank = _two_state_bank()
        query_z = torch.tensor([[0.0, 0.0], [10.0, 0.0]])
        y_ord = torch.zeros(2, 2, dtype=torch.long)
        pi, diagnostics = compute_pi_tech_from_bank(
            query_z,
            y_ord,
            torch.tensor([0.0, 1.0]),
            torch.ones(2),
            torch.zeros(2, dtype=torch.long),
            torch.zeros(2, dtype=torch.long),
            torch.tensor([10, 11]),
            torch.tensor([1000, 1001]),
            torch.tensor([1, 0]),
            bank,
            self._config(),
        )
        self.assertTrue(torch.equal(diagnostics["effective_k"], torch.tensor([4, 4])))
        self.assertTrue(bool(diagnostics["valid_query"].all()))
        self.assertGreater(float(pi[0, 0]), 0.5)
        self.assertEqual(float(pi[1, 0]), 0.0)

    def test_separate_bank_k_preserves_canonical_inbatch_k(self) -> None:
        torch.manual_seed(7)
        z = torch.randn(12, 3)
        y = torch.randint(0, 2, (12, 5), dtype=torch.long)
        depth = (y > 0).float().sum(dim=1)
        rate = (y > 0).float().mean(dim=0)
        canonical = LatentKNNPiTechConfig(enabled=True, mode="in_batch", k=8)
        bank_mode = LatentKNNPiTechConfig(
            enabled=True,
            mode="same_celltype_donor_bank",
            k=16,
            in_batch_k=8,
            bank_min_distinct_donors=4,
        )
        self.assertTrue(
            torch.equal(
                compute_pi_tech_in_batch(z, y, depth, rate, canonical),
                compute_pi_tech_in_batch(z, y, depth, rate, bank_mode),
            )
        )

    def test_sex_linked_gene_uses_same_sex_and_unknown_is_fail_closed(self) -> None:
        bank = _two_state_bank()
        query_z = torch.zeros(3, 2)
        y_ord = torch.zeros(3, 2, dtype=torch.long)
        pi, diagnostics = compute_pi_tech_from_bank(
            query_z,
            y_ord,
            torch.tensor([0.0, 1.0, 2.0]),
            torch.ones(2),
            torch.zeros(3, dtype=torch.long),
            torch.zeros(3, dtype=torch.long),
            torch.tensor([10, 11, 12]),
            torch.tensor([1000, 1001, 1002]),
            torch.tensor([1, 0, -1]),
            bank,
            self._config(),
        )
        self.assertGreater(float(pi[0, 1]), 0.5)
        self.assertEqual(float(pi[1, 1]), 0.0)
        self.assertEqual(float(pi[2, 1]), 0.0)
        self.assertFalse(bool(diagnostics["sex_valid_query"][2]))

    def test_wrong_region_candidates_do_not_validate_a_query(self) -> None:
        bank = _two_state_bank()
        pi, diagnostics = compute_pi_tech_from_bank(
            torch.tensor([[0.0, 0.0]]),
            torch.zeros(1, 2, dtype=torch.long),
            torch.tensor([0.0]),
            torch.ones(2),
            torch.tensor([0]),
            torch.tensor([1]),
            torch.tensor([10]),
            torch.tensor([1000]),
            torch.tensor([1]),
            bank,
            self._config(),
        )
        self.assertFalse(bool(diagnostics["valid_query"][0]))
        self.assertTrue(torch.equal(pi, torch.zeros_like(pi)))

    def test_checkpoint_roundtrip_bit_packs_detection_exactly(self) -> None:
        bank = _two_state_bank()
        state = bank.state_dict()
        self.assertLess(
            int(state["detect_packed"].numel()), int(bank.detect.numel())
        )
        restored = TrainOnlyPiTechBank.from_state_dict(state, device="cpu")
        for name in (
            "z",
            "detect",
            "row",
            "celltype",
            "region",
            "donor",
            "sex",
        ):
            self.assertTrue(torch.equal(getattr(bank, name), getattr(restored, name)))
        self.assertEqual(bank.epoch, restored.epoch)

    def test_detection_rate_and_depth_ecdf_give_each_donor_equal_mass(self) -> None:
        # Donor 0 contributes three high-depth ON cells and donor 1 contributes
        # one low-depth OFF cell.  A row mean would be 0.75; donor balancing is
        # exactly (1 + 0) / 2 = 0.5.
        bank = TrainOnlyPiTechBank(
            z=torch.arange(4, dtype=torch.float32).unsqueeze(1),
            detect=torch.tensor(
                [[1, 1], [1, 1], [1, 1], [0, 0]], dtype=torch.uint8
            ),
            row=torch.arange(4, dtype=torch.long),
            celltype=torch.zeros(4, dtype=torch.long),
            region=torch.zeros(4, dtype=torch.long),
            donor=torch.tensor([0, 0, 0, 1], dtype=torch.long),
            sex=torch.zeros(4, dtype=torch.long),
            epoch=13,
        )
        rate = bank.query_gene_detect_rate(
            torch.zeros(2),
            torch.tensor([0]),
            torch.tensor([0]),
            torch.tensor([0]),
            shrinkage_donors=0.0,
            sex_linked_gene_indices=(),
            match_region=True,
        )
        self.assertTrue(torch.equal(rate, torch.full((1, 2), 0.5)))
        depth_q = bank.query_depth_quantile(
            torch.tensor([1.0]),
            torch.tensor([0]),
            torch.tensor([0]),
            match_region=True,
        )
        self.assertTrue(torch.equal(depth_q, torch.tensor([0.5])))

    def test_cached_statistics_preserve_live_ema_and_reference_parity(self) -> None:
        bank = _two_state_bank()
        query_celltype = torch.tensor([0, 0, 1], dtype=torch.long)
        query_region = torch.tensor([0, 0, 0], dtype=torch.long)
        query_sex = torch.tensor([1, 0, 0], dtype=torch.long)
        for global_rate in (
            torch.tensor([0.1, 0.2]),
            torch.tensor([0.8, 0.9]),
        ):
            reference = bank._query_gene_detect_rate_reference(
                global_rate,
                query_celltype,
                query_region,
                query_sex,
                shrinkage_donors=8.0,
                sex_linked_gene_indices=(1,),
                match_region=True,
            )
            cached = bank.query_gene_detect_rate(
                global_rate,
                query_celltype,
                query_region,
                query_sex,
                shrinkage_donors=8.0,
                sex_linked_gene_indices=(1,),
                match_region=True,
            )
            self.assertTrue(
                torch.allclose(cached, reference, atol=1.0e-7, rtol=0.0)
            )

        depth = torch.tensor([0.0, 1.0, 2.0])
        depth_reference = bank._query_depth_quantile_reference(
            depth,
            query_celltype,
            query_region,
            match_region=True,
        )
        depth_cached = bank.query_depth_quantile(
            depth,
            query_celltype,
            query_region,
            match_region=True,
            sex_linked_gene_indices=(1,),
        )
        self.assertTrue(torch.equal(depth_cached, depth_reference))
        self.assertEqual(bank._query_statistics_cache_builds, 1)

    def test_train_index_selection_is_deterministic_and_group_bounded(self) -> None:
        rows = np.asarray([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int64)
        celltype = np.zeros(9, dtype=np.int64)
        donor = np.zeros(9, dtype=np.int64)
        celltype[rows] = np.asarray([0, 0, 0, 0, 1, 1, 1, 1])
        donor[rows] = np.asarray([0, 0, 0, 1, 0, 0, 1, 1])
        first = select_donor_celltype_bank_indices(
            row_indices=rows,
            celltype_ids=celltype,
            donor_ids=donor,
            per_group=2,
            seed=17,
        )
        second = select_donor_celltype_bank_indices(
            row_indices=rows,
            celltype_ids=celltype,
            donor_ids=donor,
            per_group=2,
            seed=17,
        )
        self.assertTrue(np.array_equal(first, second))
        group_counts: dict[tuple[int, int], int] = {}
        for local_index in first.tolist():
            row = rows[local_index]
            key = (int(celltype[row]), int(donor[row]))
            group_counts[key] = group_counts.get(key, 0) + 1
        self.assertTrue(group_counts)
        self.assertLessEqual(max(group_counts.values()), 2)

    def test_two_bank_slots_are_region_stratified_when_both_exist(self) -> None:
        rows = np.arange(6, dtype=np.int64)
        celltype = np.zeros(6, dtype=np.int64)
        donor = np.zeros(6, dtype=np.int64)
        region = np.asarray([0, 0, 0, 1, 1, 1], dtype=np.int64)
        selected = select_donor_celltype_bank_indices(
            row_indices=rows,
            celltype_ids=celltype,
            donor_ids=donor,
            region_ids=region,
            per_group=2,
            seed=19,
        )
        self.assertEqual(set(region[rows[selected]].tolist()), {0, 1})


class _BankDataset(Dataset):
    def __init__(self) -> None:
        self.row_idx = np.arange(12, dtype=np.int64)
        self.celltype_ids = np.zeros(12, dtype=np.int64)
        self.donor_ids = np.repeat(np.arange(6, dtype=np.int64), 2)
        self.region_ids = np.tile(np.asarray([0, 1], dtype=np.int64), 6)
        self.sex_ids = np.tile(np.asarray([0, 1], dtype=np.int64), 6)
        self.return_thinned = True
        self.thin_calls = 0

    def __len__(self) -> int:
        return 12

    def __getitem__(self, index: int):
        if self.return_thinned:
            self.thin_calls += 1
        return {
            "y_ord": torch.tensor([index % 2, (index + 1) % 2]),
            "row_index": torch.tensor(index, dtype=torch.long),
            "celltype_id": torch.tensor(0, dtype=torch.long),
            "region_id": torch.tensor(index % 2, dtype=torch.long),
            "donor_id": torch.tensor(index // 2, dtype=torch.long),
            "sex_id": torch.tensor(index % 2, dtype=torch.long),
        }


class _BankSystem(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.compute_decoder_calls = []

    def forward(
        self, batch, *, sample_latent: bool, compute_decoder: bool = True
    ):
        del sample_latent
        self.compute_decoder_calls.append(bool(compute_decoder))
        row = batch["row_index"].float()
        return SimpleNamespace(z_perp=torch.stack((row, row * 0.0), dim=1))


class PiTechBankTrainerIntegrationTest(unittest.TestCase):
    def test_epoch_diagnostics_use_position_counts_and_log_support(self) -> None:
        trainer = object.__new__(V8Trainer)
        trainer.device = torch.device("cpu")
        trainer.pi_tech_config = LatentKNNPiTechConfig(
            enabled=True,
            mode="same_celltype_donor_bank",
            k=4,
            bank_min_distinct_donors=4,
            diagnostic_min_positions=1,
            diagnostic_min_cells=1,
            diagnostic_min_distinct_donors=1,
        )
        trainer._pi_tech_diag_n_celltypes = 2
        trainer._pi_tech_diag_n_regions = 2
        trainer._pi_tech_diag_n_donors = 3
        trainer._v8_sync_group = None
        trainer._reset_pi_tech_epoch_diagnostics()
        batch = {
            "celltype_id": torch.tensor([0, 0, 1]),
            "region_id": torch.tensor([0, 0, 1]),
            "donor_id": torch.tensor([0, 1, 2]),
        }
        trainer._accumulate_pi_tech_bank_diagnostics(
            batch,
            torch.tensor([True, False, True]),
            torch.tensor([4, 2, 4]),
        )
        pi = torch.tensor([[0.2, 0.9], [0.4, 0.8], [0.6, 0.7]])
        proven = torch.tensor(
            [[True, False], [True, True], [False, False]]
        )
        stable = torch.tensor(
            [[False, True], [False, False], [True, True]]
        )
        trainer._accumulate_pi_tech_thinning_diagnostics(
            batch, pi, proven, stable
        )
        metrics = trainer._finalize_pi_tech_epoch_diagnostics()
        self.assertAlmostEqual(
            metrics["metric/pi_tech_at_proven_dropout"],
            (0.2 + 0.4 + 0.8) / 3.0,
        )
        self.assertAlmostEqual(
            metrics["metric/pi_tech_at_stable_zero"],
            (0.9 + 0.6 + 0.7) / 3.0,
        )
        self.assertEqual(metrics["count/pi_tech_proven_positions"], 3.0)
        self.assertEqual(
            metrics["count/pi_tech_thinning_distinct_donors"], 3.0
        )
        self.assertEqual(
            metrics[
                "metric/pi_tech_bank_valid_query_fraction_ct0_region0"
            ],
            0.5,
        )

    def test_sex_linked_path_preserves_phase1_then_is_fail_closed_at_full_bank(self) -> None:
        cfg = LatentKNNPiTechConfig(
            enabled=True,
            mode="same_celltype_donor_bank",
            k=4,
            bank_start_epoch=13,
            bank_ramp_epochs=4,
            bank_min_distinct_donors=4,
            sex_linked_gene_indices=(1,),
            pi_cap=0.7,
        )
        trainer = object.__new__(V8Trainer)
        trainer.pi_tech_config = cfg
        trainer._pi_tech_epoch = 12
        trainer._pi_tech_bank = None
        batch = {
            "y_ord": torch.tensor([[0, 0], [1, 1]], dtype=torch.long),
        }
        model_out = SimpleNamespace(
            z_perp=torch.tensor([[0.0, 0.0], [0.1, 0.0]])
        )
        pi, _ = trainer._compute_pi_tech(
            batch, model_out, batch["y_ord"], torch.ones(2)
        )
        self.assertGreater(float(pi[0, 0]), 0.0)
        self.assertGreater(float(pi[0, 1]), 0.0)

        trainer._pi_tech_epoch = 16
        trainer._pi_tech_bank = _two_state_bank()
        batch.update(
            {
                "celltype_id": torch.zeros(2, dtype=torch.long),
                # The bank has only region 0, so same-sex bank support is
                # deliberately invalid for these region-1 queries.
                "region_id": torch.ones(2, dtype=torch.long),
                "donor_id": torch.tensor([10, 11], dtype=torch.long),
                "row_index": torch.tensor([1000, 1001], dtype=torch.long),
                "sex_id": torch.tensor([1, 0], dtype=torch.long),
            }
        )
        pi, _ = trainer._compute_pi_tech(
            batch, model_out, batch["y_ord"], torch.ones(2)
        )
        self.assertEqual(float(pi[0, 0]), 0.0)
        self.assertEqual(float(pi[0, 1]), 0.0)

        trainer.pi_tech_config = replace(
            cfg, bank_invalid_fallback="in_batch"
        )
        legacy_pi, _ = trainer._compute_pi_tech(
            batch, model_out, batch["y_ord"], torch.ones(2)
        )
        self.assertGreater(float(legacy_pi[0, 0]), 0.0)

    def test_epoch_refresh_restores_dataset_and_model_state_without_thinning(self) -> None:
        dataset = _BankDataset()
        cfg = LatentKNNPiTechConfig(
            enabled=True,
            mode="same_celltype_donor_bank",
            k=4,
            bank_start_epoch=13,
            bank_ramp_epochs=4,
            bank_per_donor_celltype=2,
            bank_batch_size=3,
            bank_min_distinct_donors=4,
        )
        trainer = object.__new__(V8Trainer)
        trainer.v8_enabled = True
        trainer.pi_tech_config = cfg
        trainer._pi_tech_bank_dataset = dataset
        trainer._pi_tech_bank_local_indices = select_donor_celltype_bank_indices(
            row_indices=dataset.row_idx,
            celltype_ids=dataset.celltype_ids,
            donor_ids=dataset.donor_ids,
            region_ids=dataset.region_ids,
            per_group=2,
            seed=cfg.bank_seed,
        )
        trainer._pi_tech_bank = None
        trainer.system = _BankSystem()
        trainer.system.train()
        trainer.device = torch.device("cpu")
        trainer.amp_enabled = False
        trainer.gene_detect_rate = torch.ones(2)
        trainer._is_rank0 = lambda: False
        trainer._refresh_pi_tech_bank(13)
        self.assertIsNotNone(trainer._pi_tech_bank)
        self.assertEqual(trainer._pi_tech_bank.epoch, 13)
        self.assertTrue(dataset.return_thinned)
        self.assertEqual(dataset.thin_calls, 0)
        self.assertTrue(trainer.system.training)
        self.assertTrue(trainer.system.compute_decoder_calls)
        self.assertFalse(any(trainer.system.compute_decoder_calls))
        restored = TrainOnlyPiTechBank.from_state_dict(
            trainer._pi_tech_bank.state_dict(), device="cpu"
        )
        trainer._validate_pi_tech_bank_against_train_dataset(restored)
        bad_rows = restored.row.clone()
        bad_rows[0] = 999
        tampered = replace(restored, row=bad_rows)
        with self.assertRaisesRegex(RuntimeError, "train-only selection"):
            trainer._validate_pi_tech_bank_against_train_dataset(tampered)

    def test_sex_linked_symbols_resolve_to_exact_ensembl_order(self) -> None:
        dataset = SimpleNamespace(
            spec=SimpleNamespace(
                gene_names=np.asarray(["ENSG_XIST", "ENSG_KDM5D", "ENSG_A"])
            ),
            gene_symbols=np.asarray(["XIST", "KDM5D", "A"]),
        )
        observed = _resolve_pi_tech_sex_linked_gene_indices(
            {
                "sex_linked_gene_names": ["KDM5D", "ENSG_XIST"],
                "sex_linked_expected_resolution": {
                    "KDM5D": {
                        "ensembl_id": "ENSG_KDM5D",
                        "decoder_index": 1,
                    }
                },
            },
            dataset,
        )
        self.assertEqual(observed, (0, 1))
        with self.assertRaisesRegex(ValueError, "resolution contract changed"):
            _resolve_pi_tech_sex_linked_gene_indices(
                {
                    "sex_linked_gene_names": ["KDM5D"],
                    "sex_linked_expected_resolution": {
                        "KDM5D": {
                            "ensembl_id": "ENSG_KDM5D",
                            "decoder_index": 0,
                        }
                    },
                },
                dataset,
            )


if __name__ == "__main__":
    unittest.main()
