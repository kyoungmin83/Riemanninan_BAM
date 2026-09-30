from __future__ import annotations

from types import SimpleNamespace
import unittest

from kmlee_bam.training.core_trainer import RunningAverages, StepOutput


class RunningAveragesSparseMetricTest(unittest.TestCase):
    def test_sparse_metric_uses_only_batches_where_it_was_observed(self) -> None:
        averages = RunningAverages()
        averages.update_from_step(
            StepOutput(
                loss=SimpleNamespace(
                    details={"loss/total": 1.0, "metric/sparse_probe": 0.08}
                ),
                grad_norm=2.0,
                lr=1.0e-4,
                batch_size=2,
            )
        )
        averages.update_from_step(
            StepOutput(
                loss=SimpleNamespace(details={"loss/total": 3.0}),
                grad_norm=None,
                lr=1.0e-4,
                batch_size=2,
            )
        )
        observed = averages.compute()
        self.assertEqual(observed["loss/total"], 2.0)
        self.assertEqual(observed["metric/sparse_probe"], 0.08)
        self.assertEqual(observed["optim/grad_norm"], 2.0)
        self.assertEqual(observed["optim/lr"], 1.0e-4)

    def test_nonfinite_value_does_not_change_that_keys_denominator(self) -> None:
        averages = RunningAverages()
        averages.update_from_step(
            StepOutput(
                loss=SimpleNamespace(
                    details={"loss/total": 1.0, "metric/probe": float("nan")}
                ),
                grad_norm=None,
                lr=1.0e-4,
                batch_size=3,
            )
        )
        averages.update_from_step(
            StepOutput(
                loss=SimpleNamespace(
                    details={"loss/total": 1.0, "metric/probe": 0.5}
                ),
                grad_norm=None,
                lr=1.0e-4,
                batch_size=1,
            )
        )
        observed = averages.compute()
        self.assertEqual(observed["metric/probe"], 0.5)


if __name__ == "__main__":
    unittest.main()
