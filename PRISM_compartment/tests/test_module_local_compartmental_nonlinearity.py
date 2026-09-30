from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np
import torch

from kmlee_bam.data.module_local_reliability import (
    MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA,
    load_module_local_compartment_graph_artifact,
)
from kmlee_bam.model.precision_medicine import PrecisionMedicineHead
from scripts.training.build_prism_module_local_compartment_graph_20260827 import (
    build_adjacency,
)
from tests.test_precision_module_local_residual import _config, _forward, _source


def _adjacency(modules: int, *, ring: bool = True) -> torch.Tensor:
    result = torch.zeros(modules, modules)
    if ring:
        for index in range(modules):
            result[index, (index + 1) % modules] = 1.0
    return result


def _head(
    *,
    ring: bool = True,
    variant: str = "compartmental_threshold",
) -> PrecisionMedicineHead:
    torch.manual_seed(1234)
    source = _source()
    modules = int(source["source_module"].shape[-1])
    config = replace(
        _config(local=True),
        module_local_nonlinear_enabled=True,
        module_local_nonlinear_variant=variant,
        module_local_nonlinear_start_epoch=21,
        module_local_nonlinear_ramp_epochs=4,
    )
    head = PrecisionMedicineHead(
        config=config,
        module_local_compartment_adjacency=_adjacency(modules, ring=ring),
        **source,
    )
    head.eval()
    return head


class CompartmentalNonlinearityTest(unittest.TestCase):
    def test_dormant_full_loss_has_only_finite_module_local_gradients(self) -> None:
        head = _head()
        head.train()
        output = _forward(head)
        loss = (
            output.total_module_coeff.square().mean()
            + output.module_local_center_loss
            + output.module_local_size_loss
            + output.module_local_pathology_leak_loss
            + output.module_local_age_leak_loss
        )

        loss.backward()

        observed = 0
        for name, parameter in head.named_parameters():
            if not name.startswith("module_local_") or parameter.grad is None:
                continue
            observed += 1
            self.assertTrue(torch.isfinite(parameter.grad).all(), name)
        self.assertGreater(observed, 0)

    def test_personal_span_projection_is_bfloat16_safe_and_differentiable(self) -> None:
        head = _head()
        value = torch.randn(
            3,
            head.n_modules,
            dtype=torch.bfloat16,
            requires_grad=True,
        )
        celltype = torch.tensor([0, 1, 2])
        coverage = torch.ones_like(value, dtype=torch.bool)

        residual = head._project_module_local_from_personal_span(
            value,
            celltype,
            coverage,
        )

        self.assertEqual(residual.dtype, torch.bfloat16)
        self.assertTrue(torch.isfinite(residual).all())
        residual.float().square().mean().backward()
        self.assertIsNotNone(value.grad)
        self.assertTrue(torch.isfinite(value.grad).all())

    def test_zero_nonlinear_ramp_is_exact_v1_identity(self) -> None:
        torch.manual_seed(1234)
        v1 = PrecisionMedicineHead(config=_config(local=True), **_source())
        v1.eval()
        v2 = _head()
        missing, unexpected = v2.load_state_dict(v1.state_dict(), strict=False)
        self.assertFalse(unexpected)
        self.assertTrue(missing)
        self.assertTrue(all(name.startswith("module_local_") for name in missing))
        with torch.no_grad():
            v1.module_local_output_gate.fill_(0.25)
            v2.module_local_output_gate.fill_(0.25)
            v1.module_local_scale.fill_(1.0)
            v2.module_local_scale.fill_(1.0)
        donor = torch.tensor([0, 1, 2])
        celltype = torch.tensor([1, 1, 1])
        self.assertTrue(
            torch.equal(
                v1.infer_module_local(donor, celltype)[0],
                v2.infer_module_local(donor, celltype)[0],
            )
        )
        v2.infer_module_local(donor, celltype)[0].sum().backward()
        for family in ("mix", "threshold", "slope", "gain"):
            for suffix in ("global", "celltype_delta"):
                parameter = getattr(
                    v2, f"module_local_nonlinear_{family}_{suffix}"
                )
                self.assertIsNotNone(parameter.grad)
                self.assertTrue(
                    torch.equal(parameter.grad, torch.zeros_like(parameter.grad))
                )

    def test_threshold_innovation_has_zero_value_and_zero_tangent_at_origin(self) -> None:
        head = _head()
        transformed = torch.zeros(2, head.n_modules, requires_grad=True)
        coverage = torch.ones_like(transformed, dtype=torch.bool)
        innovation, _ = head._module_local_nonlinear_transform(
            transformed, torch.tensor([0, 1]), coverage
        )
        self.assertTrue(torch.equal(innovation, torch.zeros_like(innovation)))
        innovation.sum().backward()
        self.assertIsNotNone(transformed.grad)
        self.assertLess(float(transformed.grad.abs().max()), 1.0e-7)

    def test_graph_affects_output_only_after_nonlinear_ramp_opens(self) -> None:
        ring = _head(ring=True)
        isolated = _head(ring=False)
        isolated.load_state_dict(
            {
                key: value
                for key, value in ring.state_dict().items()
                if key != "module_local_compartment_adjacency"
            },
            strict=False,
        )
        donor = torch.tensor([0, 1, 2])
        celltype = torch.tensor([1, 1, 1])
        with torch.no_grad():
            for head in (ring, isolated):
                head.module_local_output_gate.fill_(0.25)
                head.module_local_scale.fill_(1.0)
        ring_off = ring.infer_module_local(donor, celltype)[0]
        isolated_off = isolated.infer_module_local(donor, celltype)[0]
        self.assertTrue(torch.equal(ring_off, isolated_off))
        ring.set_module_local_nonlinear_epoch(24)
        isolated.set_module_local_nonlinear_epoch(24)
        ring_on = ring.infer_module_local(donor, celltype)[0]
        isolated_on = isolated.infer_module_local(donor, celltype)[0]
        self.assertFalse(torch.equal(ring_on, isolated_on))

    def test_all_bounded_nonlinear_families_receive_finite_gradients(self) -> None:
        head = _head()
        with torch.no_grad():
            head.module_local_output_gate.fill_(0.25)
            head.module_local_scale.fill_(1.0)
        head.set_module_local_nonlinear_epoch(24)
        donor = torch.tensor([0, 1, 2])
        celltype = torch.tensor([0, 1, 2])
        residual = head.infer_module_local_with_diagnostics(
            donor, celltype
        )[0]
        residual.square().mean().backward()
        for family in ("mix", "threshold", "slope", "gain"):
            parameter = getattr(head, f"module_local_nonlinear_{family}_global")
            self.assertIsNotNone(parameter.grad)
            self.assertTrue(torch.isfinite(parameter.grad).all())
            self.assertGreater(float(parameter.grad.abs().sum()), 0.0)

    def test_graph_linear_control_freezes_unused_families_and_is_input_linear(self) -> None:
        nonlinear = _head(variant="compartmental_threshold")
        control = _head(variant="graph_linear_control")
        nonlinear_parameters = {
            name: parameter.numel()
            for name, parameter in nonlinear.named_parameters()
            if name.startswith("module_local_nonlinear_")
        }
        control_parameters = {
            name: parameter.numel()
            for name, parameter in control.named_parameters()
            if name.startswith("module_local_nonlinear_")
        }
        self.assertEqual(nonlinear_parameters, control_parameters)
        control_trainable = {
            name
            for name, parameter in control.named_parameters()
            if name.startswith("module_local_nonlinear_")
            and parameter.requires_grad
        }
        self.assertTrue(control_trainable)
        self.assertTrue(
            all(
                "_mix_" in name or "_gain_" in name
                for name in control_trainable
            )
        )
        for family in ("threshold", "slope"):
            for suffix in ("global", "celltype_delta"):
                self.assertFalse(
                    getattr(
                        control,
                        f"module_local_nonlinear_{family}_{suffix}",
                    ).requires_grad
                )

        x = torch.linspace(-2.0, 2.0, nonlinear.n_modules).unsqueeze(0)
        coverage = torch.ones_like(x, dtype=torch.bool)
        celltype = torch.tensor([0])
        control_x, _ = control._module_local_nonlinear_transform(
            x, celltype, coverage
        )
        control_2x, _ = control._module_local_nonlinear_transform(
            2.0 * x, celltype, coverage
        )
        self.assertTrue(torch.allclose(control_2x, 2.0 * control_x))

        nonlinear_x, _ = nonlinear._module_local_nonlinear_transform(
            x, celltype, coverage
        )
        nonlinear_2x, _ = nonlinear._module_local_nonlinear_transform(
            2.0 * x, celltype, coverage
        )
        self.assertFalse(torch.allclose(nonlinear_2x, 2.0 * nonlinear_x))

        diagnostics = control.module_local_nonlinear_parameter_diagnostics()
        self.assertAlmostEqual(float(diagnostics["linear_scale_mean"]), 0.1)

    def test_stacked_input_jacobian_has_two_vs_four_effective_fields(self) -> None:
        def _small_head(variant: str) -> PrecisionMedicineHead:
            source = _source()
            source["source_module"] = source["source_module"][..., :8]
            source["module_local_reliability"] = source[
                "module_local_reliability"
            ][..., :8]
            config = replace(
                _config(local=True),
                module_local_nonlinear_enabled=True,
                module_local_nonlinear_variant=variant,
            )
            return PrecisionMedicineHead(
                config=config,
                module_local_compartment_adjacency=_adjacency(8),
                **source,
            ).double()

        torch.manual_seed(71)
        base = torch.randn(12, 8, dtype=torch.float64)
        scale = torch.tensor(
            [0.25, 0.65, 1.25, 2.5] * 3, dtype=torch.float64
        ).unsqueeze(1)
        transformed = base * scale
        celltype = torch.tensor([0] * 4 + [1] * 4 + [2] * 4)
        coverage = torch.ones_like(transformed, dtype=torch.bool)

        def _rank(variant: str) -> int:
            head = _small_head(variant)
            output = head._module_local_nonlinear_transform(
                transformed, celltype, coverage
            )[0].reshape(-1)
            parameters = [
                parameter
                for name, parameter in head.named_parameters()
                if name.startswith("module_local_nonlinear_")
                and parameter.requires_grad
            ]
            rows = []
            for index in range(int(output.numel())):
                gradients = torch.autograd.grad(
                    output[index],
                    parameters,
                    retain_graph=True,
                    allow_unused=False,
                )
                rows.append(torch.cat([value.reshape(-1) for value in gradients]))
            singular_values = torch.linalg.svdvals(torch.stack(rows))
            tolerance = (
                max(int(output.numel()), sum(p.numel() for p in parameters))
                * torch.finfo(singular_values.dtype).eps
                * singular_values.max()
                * 10.0
            )
            return int((singular_values > tolerance).sum())

        self.assertEqual(_rank("graph_linear_control"), 48)
        self.assertEqual(_rank("compartmental_threshold"), 96)

    def test_two_neighbor_drives_can_cross_threshold_supralinearly(self) -> None:
        head = _head(variant="compartmental_threshold")
        with torch.no_grad():
            head.module_local_compartment_adjacency.zero_()
            head.module_local_compartment_adjacency[0, 1] = 0.5
            head.module_local_compartment_adjacency[0, 2] = 0.5
            for family, raw in {
                "mix": 20.0,
                "threshold": -20.0,
                "slope": 20.0,
                "gain": 20.0,
            }.items():
                getattr(
                    head, f"module_local_nonlinear_{family}_global"
                ).fill_(raw)
                getattr(
                    head, f"module_local_nonlinear_{family}_celltype_delta"
                ).zero_()
        first = torch.zeros(1, head.n_modules)
        second = torch.zeros_like(first)
        first[0, 1] = 1.0
        second[0, 2] = 1.0
        coverage = torch.ones_like(first, dtype=torch.bool)
        celltype = torch.tensor([0])
        first_out, _ = head._module_local_nonlinear_transform(
            first, celltype, coverage
        )
        second_out, _ = head._module_local_nonlinear_transform(
            second, celltype, coverage
        )
        joint_out, _ = head._module_local_nonlinear_transform(
            first + second, celltype, coverage
        )
        self.assertGreater(
            float(joint_out[0, 0]),
            float(first_out[0, 0] + second_out[0, 0]),
        )

    def test_nonlinear_curriculum_is_persistent_and_resumable(self) -> None:
        head = _head()
        self.assertEqual(head.set_module_local_nonlinear_epoch(20), 0.0)
        self.assertEqual(head.set_module_local_nonlinear_epoch(21), 0.25)
        self.assertEqual(head.set_module_local_nonlinear_epoch(24), 1.0)
        resumed = _head()
        resumed.load_state_dict(head.state_dict(), strict=True)
        self.assertEqual(float(resumed.module_local_nonlinear_scale), 1.0)
        self.assertEqual(resumed._module_local_nonlinear_scale_value, 1.0)

    def test_primary_414_by_24_parameter_budget_has_no_trainable_dense_graph(self) -> None:
        config = replace(
            _config(local=True),
            module_local_rank=8,
            module_local_nonlinear_enabled=True,
        )
        source = _source()
        source["source_module"] = torch.zeros(2, 6, 414)
        source["source_latent"] = torch.zeros(2, 6, 4)
        source["source_observed"] = torch.ones(2, 6, dtype=torch.bool)
        source["source_reliability"] = torch.ones(2, 6)
        source["module_local_reliability"] = torch.ones(6, 414)
        source["n_celltypes"] = 24
        head = PrecisionMedicineHead(
            config=config,
            module_local_compartment_adjacency=torch.zeros(414, 414),
            **source,
        )
        nonlinear = [
            parameter
            for name, parameter in head.named_parameters()
            if name.startswith("module_local_nonlinear_")
        ]
        self.assertEqual(sum(parameter.numel() for parameter in nonlinear), 41400)
        graph_control = PrecisionMedicineHead(
            config=replace(
                config,
                module_local_nonlinear_variant="graph_linear_control",
            ),
            module_local_compartment_adjacency=torch.zeros(414, 414),
            **source,
        )
        graph_stored = [
            parameter
            for name, parameter in graph_control.named_parameters()
            if name.startswith("module_local_nonlinear_")
        ]
        graph_trainable = [
            parameter for parameter in graph_stored if parameter.requires_grad
        ]
        self.assertEqual(sum(parameter.numel() for parameter in graph_stored), 41400)
        self.assertEqual(
            sum(parameter.numel() for parameter in graph_trainable), 20700
        )
        self.assertFalse(
            any(tuple(parameter.shape) == (414, 414) for parameter in head.parameters())
        )
        self.assertEqual(
            tuple(head.module_local_compartment_adjacency.shape), (414, 414)
        )


class CompartmentGraphArtifactTest(unittest.TestCase):
    def _write(self, path: Path, *, test_used: bool = False) -> None:
        adjacency = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.float32)
        np.savez_compressed(
            path,
            schema_version=np.asarray(MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA),
            adjacency=adjacency,
            module_names=np.asarray(["m0", "m1"]),
            registry_sha256=np.asarray("a" * 64),
            topk=np.asarray(1),
            minimum_jaccard=np.asarray(0.05),
            validation_donors_used=np.asarray(False),
            test_donors_used=np.asarray(test_used),
        )

    def test_builder_is_sparse_deterministic_and_row_normalized(self) -> None:
        membership = np.asarray(
            [
                [1, 1, 0, 0],
                [1, 1, 1, 0],
                [0, 0, 0, 1],
                [1, 0, 0, 0],
            ],
            dtype=np.float32,
        )
        first = build_adjacency(membership, topk=1, minimum_jaccard=0.1)
        second = build_adjacency(membership, topk=1, minimum_jaccard=0.1)
        self.assertTrue(np.array_equal(first, second))
        self.assertTrue(np.array_equal(np.diag(first), np.zeros(4)))
        self.assertTrue(np.all((first > 0).sum(axis=1) <= 1))
        row_sum = first.sum(axis=1)
        self.assertTrue(np.allclose(row_sum[row_sum > 0], 1.0))
        self.assertEqual(float(row_sum[2]), 0.0)

    def test_loader_checks_order_hash_and_sealed_splits(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            good = Path(directory) / "good.npz"
            self._write(good)
            artifact = load_module_local_compartment_graph_artifact(
                good,
                expected_module_names=("m0", "m1"),
                expected_registry_sha256="a" * 64,
            )
            self.assertEqual(artifact.topk, 1)
            bad = Path(directory) / "bad.npz"
            self._write(bad, test_used=True)
            with self.assertRaisesRegex(ValueError, "test donors"):
                load_module_local_compartment_graph_artifact(
                    bad,
                    expected_module_names=("m0", "m1"),
                    expected_registry_sha256="a" * 64,
                )


if __name__ == "__main__":
    unittest.main()
