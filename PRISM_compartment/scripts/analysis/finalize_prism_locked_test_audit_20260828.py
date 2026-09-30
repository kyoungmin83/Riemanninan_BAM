#!/usr/bin/env python3
"""Verify every locked test artifact and mark the recovered audit complete."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path


CANDIDATES = (
    "rank2_e20",
    "sv6_e33_protected",
    "sv6_e42",
    "sv6_e47",
    "sv6_e50",
    "sv6_e55",
)


def load(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--lock-sha256", required=True)
    args = parser.parse_args()
    root = Path(args.root)

    summaries = []
    cell_counts = set()
    for candidate_id in CANDIDATES:
        directory = root / candidate_id
        module = load(directory / "module_disease_test_bal40.json")
        leakage = load(directory / "leakage_test.json")
        by_celltype = load(directory / "leakage_by_celltype_test.json")
        strict = load(directory / "strict_reference_leakage_test.json")

        assert module["evaluation_provenance"]["split"] == "test"
        assert module["evaluation_provenance"]["selected_cells"] == 8478
        assert len(module["per_donor"]) == 9
        assert leakage["test_used"] is True
        assert leakage["cell_counts"].get("val", 0) == 0
        assert leakage["cell_counts"]["test"] > 0
        assert by_celltype["test_used"] is True
        assert by_celltype["selection_locked_before_test"] is True
        assert by_celltype["pretest_lock_sha256"] == args.lock_sha256
        assert by_celltype["target_split"] == "test"
        assert len(by_celltype["per_celltype"]) == 24
        assert strict["test_used"] is True
        assert strict["selection_locked_before_test"] is True
        assert strict["pretest_lock_sha256"] == args.lock_sha256
        assert strict["target_split"] == "test"
        cell_counts.add(
            (leakage["cell_counts"]["train"], leakage["cell_counts"]["test"])
        )
        summaries.append(
            {
                "id": candidate_id,
                "module_test_cells": module["evaluation_provenance"]["selected_cells"],
                "module_test_donors": len(module["per_donor"]),
                "leak_train_cells": leakage["cell_counts"]["train"],
                "leak_test_cells": leakage["cell_counts"]["test"],
                "celltypes": len(by_celltype["per_celltype"]),
                "strict_test_donors": strict["donor_counts"]["test"],
            }
        )
        (directory / "TEST_AUDIT_COMPLETE").touch()

    assert len(cell_counts) == 1, cell_counts
    completed_at = dt.datetime.now(dt.timezone.utc).isoformat()
    receipt = {
        "schema_version": "kmlee_bam.prism_locked_test_completion.v1",
        "completed_at": completed_at,
        "official_test_status": "consumed",
        "selection_locked_before_test": True,
        "test_results_may_change_locked_selection": False,
        "pretest_lock_sha256": args.lock_sha256,
        "candidates": summaries,
        "recovery_notes": [
            "rank2 extraction used explicit_scale=1.0 for the single missing legacy buffer, matching the pre-buffer implementation",
            "cell-type test probes used a dedicated lock-gated test script because the selection-time script is intentionally validation-only",
            "completed GPU extracts were reused after the running worker file changed; no module or latent extract was recomputed for the five SV6 checkpoints",
        ],
    }
    with (root / "audit_completion.json").open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    (root / "TEST_AUDIT_COMPLETE").touch()
    print(json.dumps(receipt, ensure_ascii=False))


if __name__ == "__main__":
    main()
