from __future__ import annotations

from types import SimpleNamespace
import unittest

import torch
from torch import nn

from kmlee_bam.training.prism_module_rescue_training import (
    _module_rescue_parameter_allowlist,
)
from kmlee_bam.training.core_trainer import Trainer
from kmlee_bam.training.run_current import (
    _precision_warm_start_missing_key_allowed,
    build_optimizer,
)


class _DummyPrecision(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.cfg = SimpleNamespace(
            enabled=True,
            module_local_enabled=True,
            module_local_rank=8,
            module_local_lr_multiplier=10.0,
            module_local_rescue_start_epoch=25,
            lr_multiplier=25.0,
        )
        self.module_local_enabled = True
        self.common_global = nn.Parameter(torch.zeros(1))
        self.personal_basis = nn.Parameter(torch.zeros(1))
        for name in (
            "module_local_read",
            "module_local_write_global",
            "module_local_write_celltype_delta",
            "module_local_diagonal_global",
            "module_local_diagonal_celltype_delta",
            "module_local_output_gate",
        ):
            setattr(self, name, nn.Parameter(torch.ones(1)))
        self.register_buffer("module_local_scale", torch.tensor(1.0))


class _DummyDecoder(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.coeff_head = nn.Linear(1, 1)
        self.direct_state_head = nn.Linear(1, 1)
        self.generator_u = nn.Parameter(torch.ones(1))
        self.generator_v = nn.Parameter(torch.ones(1))
        self.generator_a = nn.Parameter(torch.ones(1))


class _DummySystem(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.backbone = nn.Linear(1, 1)
        self.decoder = _DummyDecoder()
        self.precision_head = _DummyPrecision()


class ModuleLocalTrainingContractTest(unittest.TestCase):
    def test_trainer_regularizers_and_diagnostics_accept_local_output(self) -> None:
        from tests.test_precision_module_local_residual import _forward, _head

        head = _head(local=True)
        with torch.no_grad():
            head.module_local_output_gate.fill_(0.25)
            head.module_local_scale.fill_(1.0)
        precision = _forward(head)
        genes = 5
        batch_size = int(precision.total_module_coeff.shape[0])
        precision.gene_score = torch.zeros(batch_size, genes)
        precision.branch_nll_per_gene = torch.ones(batch_size, genes)
        precision.branch_nll_per_cell = torch.ones(batch_size)
        precision.module_local_off_branch_nll_per_cell = torch.ones(batch_size)
        precision.module_local_off_full_nll_per_cell = torch.ones(batch_size)
        trainer = object.__new__(Trainer)
        trainer.system = SimpleNamespace(precision_head=head, training=True)
        trainer.ordinal_class_weights = None
        model_out = SimpleNamespace(
            precision_out=precision,
            decoder_out=SimpleNamespace(nll_per_cell=torch.ones(batch_size)),
        )
        batch = {
            "donor_balance_weight": torch.ones(batch_size),
            "y_ord": torch.zeros(batch_size, genes, dtype=torch.long),
        }
        loss_out = SimpleNamespace(total=torch.tensor(1.0), details={})
        result = Trainer._maybe_add_precision_medicine(
            trainer, batch, model_out, loss_out
        )
        self.assertTrue(torch.isfinite(result.total))
        for key in (
            "loss/prism_module_local_center",
            "loss/prism_module_local_hierarchy",
            "metric/prism_module_local_module_rms",
            "metric/prism_module_local_code_participation_rank",
            "metric/prism_module_local_off_full_nll",
        ):
            self.assertIn(key, result.details)

    def test_optimizer_has_disjoint_local_group(self) -> None:
        system = _DummySystem()
        cfg = SimpleNamespace(
            precision_medicine=system.precision_head.cfg,
            learned_generator_count=SimpleNamespace(
                enabled=False, mode="gate_only"
            ),
            optim=SimpleNamespace(
                lr=1.0e-4,
                weight_decay=1.0e-2,
                betas=(0.9, 0.999),
                eps=1.0e-8,
            ),
        )
        optimizer = build_optimizer(cfg, system)
        groups = {group.get("name"): group for group in optimizer.param_groups}
        self.assertEqual(
            set(groups),
            {"backbone_unfrozen", "precision_new", "precision_module_local"},
        )
        self.assertEqual(groups["precision_module_local"]["lr"], 1.0e-3)
        ids = [
            id(parameter)
            for group in optimizer.param_groups
            for parameter in group["params"]
        ]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(
            set(ids),
            {id(parameter) for parameter in system.parameters()},
        )

    def test_rescue_local_parameters_are_stage_gated(self) -> None:
        system = _DummySystem()
        before = set(
            _module_rescue_parameter_allowlist(
                system, include_module_local=False
            ).values()
        )
        after = set(
            _module_rescue_parameter_allowlist(
                system, include_module_local=True
            ).values()
        )
        local_names = {
            name
            for name, _ in system.named_parameters()
            if name.startswith("precision_head.module_local_")
        }
        self.assertTrue(local_names)
        self.assertTrue(local_names.isdisjoint(before))
        self.assertTrue(local_names.issubset(after))

    def test_strict_warm_start_only_allows_module_local_namespace(self) -> None:
        protected = {"precision_head.source_module"}
        self.assertTrue(
            _precision_warm_start_missing_key_allowed(
                "precision_head.module_local_read",
                protected_context_keys=protected,
                strict_precision_upgrade=True,
            )
        )
        self.assertTrue(
            _precision_warm_start_missing_key_allowed(
                "precision_head.source_module",
                protected_context_keys=protected,
                strict_precision_upgrade=True,
            )
        )
        self.assertFalse(
            _precision_warm_start_missing_key_allowed(
                "precision_head.unrelated_new_branch",
                protected_context_keys=protected,
                strict_precision_upgrade=True,
            )
        )
        self.assertTrue(
            _precision_warm_start_missing_key_allowed(
                "precision_head.legacy_initialization",
                protected_context_keys=protected,
                strict_precision_upgrade=False,
            )
        )


if __name__ == "__main__":
    unittest.main()
