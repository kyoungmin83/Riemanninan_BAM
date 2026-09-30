from __future__ import annotations

from dataclasses import replace
import inspect
from pathlib import Path
import tempfile
import unittest

import numpy as np
import torch

from kmlee_bam.data.module_local_reliability import (
    MODULE_LOCAL_OUTPUT_CAP_SCHEMA,
    MODULE_LOCAL_RELIABILITY_SCHEMA,
    load_module_local_output_cap_artifact,
    load_module_local_reliability_artifact,
    sha256_string_sequence,
)
from kmlee_bam.model.precision_medicine import (
    PrecisionMedicineConfig,
    PrecisionMedicineHead,
)


def _config(*, local: bool) -> PrecisionMedicineConfig:
    return PrecisionMedicineConfig(
        enabled=True,
        personal_rank=2,
        context_hidden_dim=8,
        module_token_dim=4,
        latent_token_dim=2,
        use_source_latent=False,
        context_embedding_dim=2,
        context_n_heads=2,
        context_n_layers=1,
        context_dropout=0.0,
        context_input_dropout=0.0,
        module_local_enabled=local,
        module_local_rank=2,
        module_local_start_epoch=17,
        module_local_ramp_epochs=8,
        module_local_rescue_start_epoch=25,
    )


def _source() -> dict[str, torch.Tensor | int]:
    donors, celltypes, regions, modules, latent = 3, 3, 2, 10, 4
    contexts = celltypes * regions
    source_module = torch.arange(
        donors * contexts * modules, dtype=torch.float32
    ).reshape(donors, contexts, modules)
    source_module = (source_module - source_module.mean()) / 50.0
    return {
        "source_module": source_module,
        "source_latent": torch.zeros(donors, contexts, latent),
        "source_observed": torch.ones(donors, contexts, dtype=torch.bool),
        "source_reliability": torch.ones(donors, contexts),
        "context_celltype": torch.tensor([0, 0, 1, 1, 2, 2]),
        "context_region": torch.tensor([0, 1, 0, 1, 0, 1]),
        "n_celltypes": celltypes,
        "n_regions": regions,
        "module_local_reliability": torch.ones(contexts, modules),
        "module_local_input_clip": 20.0,
        "module_local_output_cap": 2.0,
    }


def _head(*, local: bool) -> PrecisionMedicineHead:
    torch.manual_seed(1234)
    kwargs = _source()
    if not local:
        kwargs.pop("module_local_reliability")
        kwargs.pop("module_local_input_clip")
        kwargs.pop("module_local_output_cap")
    head = PrecisionMedicineHead(config=_config(local=local), **kwargs)
    head.eval()
    return head


def _forward(head: PrecisionMedicineHead, *, state_value: float = 0.0):
    batch = 4
    return head(
        donor_id=torch.tensor([0, 0, 1, 2]),
        celltype_id=torch.tensor([1, 1, 2, 0]),
        region_id=torch.tensor([0, 1, 0, 1]),
        age_z=torch.zeros(batch),
        age_valid=torch.ones(batch, dtype=torch.bool),
        pathology=torch.zeros(batch, 5),
        pathology_valid=torch.ones(batch, 5, dtype=torch.bool),
        state_latent=torch.full((batch, head.latent_dim), state_value),
        cell_weight=torch.ones(batch),
    )


class ModuleLocalResidualTest(unittest.TestCase):
    def test_default_off_and_exact_zero_checkpoint_identity(self) -> None:
        off = _head(local=False)
        on = _head(local=True)
        off_state = off.state_dict()
        on_state = on.state_dict()
        shared = [key for key in off_state if key in on_state]
        self.assertTrue(shared)
        for key in shared:
            self.assertTrue(torch.equal(off_state[key], on_state[key]), key)
        missing, unexpected = on.load_state_dict(off_state, strict=False)
        self.assertFalse(unexpected)
        self.assertTrue(missing)
        self.assertTrue(
            all(key.startswith("module_local_") for key in missing), missing
        )

        off_output = _forward(off)
        on_output = _forward(on)
        self.assertTrue(
            torch.equal(off_output.total_module_coeff, on_output.total_module_coeff)
        )
        self.assertTrue(torch.equal(on_output.module_local_coeff, torch.zeros_like(on_output.module_local_coeff)))
        self.assertTrue(
            torch.equal(on_output.personal_coeff, on_output.personal_rank2_coeff)
        )

    def test_target_celltype_both_regions_are_an_exact_firewall(self) -> None:
        head = _head(local=True)
        with torch.no_grad():
            head.module_local_output_gate.fill_(0.5)
            head.module_local_scale.fill_(1.0)
        donor = torch.tensor([0])
        target = torch.tensor([1])
        before = head.infer_module_local(donor, target)[0]
        original = head.source_module.clone()
        target_context = head.context_celltype == int(target.item())
        support_context = head.context_celltype != int(target.item())
        self.assertEqual(int(target_context.sum()), head.n_regions)
        with torch.no_grad():
            head.source_module[0, target_context].fill_(1.0e6)
        after_target_change = head.infer_module_local(donor, target)[0]
        self.assertTrue(torch.equal(before, after_target_change))

        with torch.no_grad():
            head.source_module.copy_(original)
            first_support = torch.nonzero(
                support_context, as_tuple=False
            ).flatten()[0]
            head.source_module[0, first_support].add_(5.0)
        after_support_change = head.infer_module_local(donor, target)[0]
        self.assertFalse(torch.equal(before, after_support_change))

    def test_no_target_latent_or_bam_input_and_state_independence(self) -> None:
        parameters = inspect.signature(
            PrecisionMedicineHead.infer_module_local
        ).parameters
        self.assertEqual(tuple(parameters), ("self", "donor_id", "celltype_id"))
        head = _head(local=True)
        with torch.no_grad():
            head.module_local_output_gate.fill_(0.5)
            head.module_local_scale.fill_(1.0)
        first = _forward(head, state_value=-100.0)
        second = _forward(head, state_value=100.0)
        self.assertTrue(torch.equal(first.module_local_coeff, second.module_local_coeff))

    def test_rank2_span_is_ridge_suppressed_without_dense_projection(self) -> None:
        head = _head(local=True)
        with torch.no_grad():
            head.personal_basis.normal_(mean=0.0, std=0.5)
            head.module_local_output_gate.fill_(0.5)
            head.module_local_scale.fill_(1.0)
        donor = torch.tensor([0, 1, 2])
        target = torch.tensor([1, 1, 1])
        residual = head.infer_module_local(donor, target)[0]
        basis = head.personal_basis[target].detach()
        overlap = torch.einsum("bpm,bm->bp", basis, residual)
        denominator = basis.norm(dim=2) * residual.norm(dim=1, keepdim=True)
        normalized_overlap = overlap.abs() / denominator.clamp_min(1.0e-12)
        normalized_overlap = torch.where(
            denominator > 0.0,
            normalized_overlap,
            torch.zeros_like(normalized_overlap),
        )
        self.assertLess(float(normalized_overlap.max()), 5.0e-4)

    def test_gate_gets_first_gradient_then_internal_adapter_opens(self) -> None:
        head = _head(local=True)
        with torch.no_grad():
            head.module_local_scale.fill_(1.0)
        donor = torch.tensor([0, 1, 2])
        target = torch.tensor([1, 1, 1])
        first = head.infer_module_local(donor, target)[0].sum()
        first.backward()
        gate_gradient = head.module_local_output_gate.grad
        self.assertIsNotNone(gate_gradient)
        self.assertTrue(torch.isfinite(gate_gradient).all())
        self.assertGreater(float(gate_gradient.abs().sum()), 0.0)
        for parameter in (
            head.module_local_read,
            head.module_local_write_global,
            head.module_local_diagonal_global,
        ):
            self.assertTrue(
                parameter.grad is None
                or torch.equal(parameter.grad, torch.zeros_like(parameter.grad))
            )

        with torch.no_grad():
            head.module_local_output_gate.add_(-0.05 * gate_gradient)
        head.zero_grad(set_to_none=True)
        second = head.infer_module_local(donor, target)[0].square().mean()
        second.backward()
        for parameter in (
            head.module_local_read,
            head.module_local_write_global,
            head.module_local_diagonal_global,
        ):
            self.assertIsNotNone(parameter.grad)
            self.assertTrue(torch.isfinite(parameter.grad).all())
            self.assertGreater(float(parameter.grad.abs().sum()), 0.0)

    def test_curriculum_scale_is_persistent_and_resumable(self) -> None:
        head = _head(local=True)
        self.assertEqual(head.set_module_local_epoch(16), 0.0)
        self.assertEqual(head.set_module_local_epoch(20), 0.5)
        state = head.state_dict()
        self.assertIn("module_local_scale", state)
        resumed = _head(local=True)
        resumed.load_state_dict(state, strict=True)
        self.assertEqual(float(resumed.module_local_scale), 0.5)
        self.assertEqual(resumed.set_module_local_epoch(24), 1.0)

    def test_primary_adapter_parameter_budget(self) -> None:
        config = replace(_config(local=True), module_local_rank=8)
        source = _source()
        source["source_module"] = torch.zeros(2, 6, 414)
        source["source_latent"] = torch.zeros(2, 6, 4)
        source["source_observed"] = torch.ones(2, 6, dtype=torch.bool)
        source["source_reliability"] = torch.ones(2, 6)
        source["module_local_reliability"] = torch.ones(6, 414)
        source["n_celltypes"] = 24
        head = PrecisionMedicineHead(config=config, **source)
        adapter_names = (
            "module_local_read",
            "module_local_write_global",
            "module_local_write_celltype_delta",
            "module_local_diagonal_global",
            "module_local_diagonal_celltype_delta",
            "module_local_output_gate",
        )
        observed = sum(getattr(head, name).numel() for name in adapter_names)
        self.assertEqual(observed, 106398)
        self.assertFalse(
            any(tuple(parameter.shape) == (414, 414) for parameter in head.parameters())
        )


class ReliabilityArtifactTest(unittest.TestCase):
    def _write(self, path: Path, *, module_names=("m0", "m1")) -> None:
        donors = ("d0", "d1")
        np.savez_compressed(
            path,
            schema_version=np.asarray(MODULE_LOCAL_RELIABILITY_SCHEMA),
            context_module_reliability=np.asarray(
                [[0.9, 0.1], [0.8, 0.7]], dtype=np.float32
            ),
            context_module_reliable_mask=np.asarray(
                [[True, False], [True, True]], dtype=bool
            ),
            context_names=np.asarray(("c0", "c1")),
            module_names=np.asarray(module_names),
            train_donor_allowlist=np.asarray(donors),
            train_donor_allowlist_sha256=np.asarray(
                sha256_string_sequence(donors)
            ),
            split_seed=np.asarray(42),
            split_rule=np.asarray("deterministic-test-rule"),
            source_context_sha256=np.asarray("a" * 64),
            registry_sha256=np.asarray("b" * 64),
            activity_dictionary_sha256=np.asarray("c" * 64),
            validation_donors_used=np.asarray(False),
            test_donors_used=np.asarray(False),
        )

    def test_loader_accepts_exact_contract_and_masks_unreliable_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reliability.npz"
            self._write(path)
            artifact = load_module_local_reliability_artifact(
                path,
                expected_context_names=("c0", "c1"),
                expected_module_names=("m0", "m1"),
                expected_train_donor_names=("d0", "d1"),
                expected_source_context_sha256="a" * 64,
                expected_registry_sha256="b" * 64,
                expected_activity_dictionary_sha256="c" * 64,
            )
            self.assertEqual(float(artifact.effective_reliability[0, 1]), 0.0)

    def test_loader_fails_closed_on_module_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reliability.npz"
            self._write(path, module_names=("m1", "m0"))
            with self.assertRaisesRegex(ValueError, "module names/order"):
                load_module_local_reliability_artifact(
                    path,
                    expected_context_names=("c0", "c1"),
                    expected_module_names=("m0", "m1"),
                    expected_train_donor_names=("d0", "d1"),
                )

    def test_output_cap_requires_frozen_rank2_provenance(self) -> None:
        donors = ("d0", "d1")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "output_cap.npz"
            np.savez_compressed(
                path,
                schema_version=np.asarray(MODULE_LOCAL_OUTPUT_CAP_SCHEMA),
                output_cap=np.asarray(0.37),
                absolute_quantile=np.asarray(0.995),
                personal_rank=np.asarray(2),
                personal_coefficient_count=np.asarray(100),
                train_donor_allowlist=np.asarray(donors),
                train_donor_allowlist_sha256=np.asarray(
                    sha256_string_sequence(donors)
                ),
                checkpoint_sha256=np.asarray("d" * 64),
                source_config_sha256=np.asarray("e" * 64),
                source_context_sha256=np.asarray("f" * 64),
                coefficient_origin=np.asarray(
                    "frozen_rank2_train_donor_personal_coeff_absolute_quantile"
                ),
                validation_donors_used=np.asarray(False),
                test_donors_used=np.asarray(False),
            )
            artifact = load_module_local_output_cap_artifact(
                path,
                expected_train_donor_names=donors,
                expected_checkpoint_sha256="d" * 64,
                expected_source_config_sha256="e" * 64,
                expected_source_context_sha256="f" * 64,
                expected_absolute_quantile=0.995,
                expected_personal_rank=2,
            )
            self.assertEqual(artifact.output_cap, 0.37)
            with self.assertRaisesRegex(ValueError, "checkpoint SHA-256"):
                load_module_local_output_cap_artifact(
                    path,
                    expected_train_donor_names=donors,
                    expected_checkpoint_sha256="a" * 64,
                    expected_source_config_sha256="e" * 64,
                    expected_source_context_sha256="f" * 64,
                    expected_absolute_quantile=0.995,
                )


if __name__ == "__main__":
    unittest.main()
