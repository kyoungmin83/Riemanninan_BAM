#!/usr/bin/env python3
"""Run the frozen extractor with one audited legacy rank2 compatibility allowance."""

from __future__ import annotations

import runpy
import sys

import prism_posthoc_system


_original_validate = prism_posthoc_system.validate_posthoc_checkpoint_compatibility


def _validate_legacy_rank2(missing, unexpected):
    missing = list(missing)
    unexpected = list(unexpected)
    legacy_key = "precision_head.explicit_scale"
    if legacy_key not in missing:
        raise RuntimeError(
            f"legacy rank2 allowance expected missing key {legacy_key!r}; got={missing}"
        )
    remaining_missing = [key for key in missing if key != legacy_key]
    result = _original_validate(remaining_missing, unexpected)
    result["ignored_legacy_missing"] = [legacy_key]
    result["legacy_default_semantics"] = (
        "explicit_scale=1.0, matching the pre-buffer implementation"
    )
    return result


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: run_legacy_rank2_extract_locked_test.py EXTRACT_SCRIPT [ARGS...]")
    extract_script = sys.argv[1]
    prism_posthoc_system.validate_posthoc_checkpoint_compatibility = _validate_legacy_rank2
    sys.argv = sys.argv[1:]
    runpy.run_path(extract_script, run_name="__main__")


if __name__ == "__main__":
    main()
