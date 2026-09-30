#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "src/kmlee_bam/model/latent_nuisance_projection.py": "eaf35a88cccba2d9ec22bb9efda7debae0110ca474f7dda8147b6fb8c65fed49",
    "src/kmlee_bam/model/lie_ordinal_decoder.py": "ba77e3e5cd59703059c093d4ed5807f3c0f47305b567daece6441e511aa59679",
    "src/kmlee_bam/model/system.py": "d6a66b75cb800ff993b87b5f8ae77f12ec0f7fa72daba4ebe8e8cbc430eb3714",
    "src/kmlee_bam/objectives/celltype_alignment.py": "2509e27b3f4a969af9162678a2d5c8261b09ef28ad6aeebfe44432602b7957ec",
    "src/kmlee_bam/training/adversarial_trainer.py": "39f21a13027b20726a4cc833b10ff9f25f7a674850c8ce31ba61f3bfedf4510b",
    "src/kmlee_bam/training/core_trainer.py": "fd90bc3e40aa198794be46f040ed03fb58c85e80cd111ecfa8d364f6753a3219",
    "src/kmlee_bam/training/learned_generator_count.py": "6564af5eaac4ee8c4e8cf0f03f81f23cafcc9833976671a81ad805daa9226053",
    "src/kmlee_bam/training/prism_module_rescue_training.py": "9451867736f4034e3505c33a5e42e0bca036bc847a2e876872bd5aec7c372d54",
    "src/kmlee_bam/training/run_current.py": "8ac9a02dcac57a866ba017d514bdc0b7a70715b270acd334b07303422712774e",
    "src/kmlee_bam/training/runner_base.py": "96ea54dd3a23a549d12b1a5e6d5044087d3bd62ecdee1ac2a3187d1150f16d87",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    failures: list[str] = []
    for relative, expected in EXPECTED.items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        observed = sha256(path)
        if observed != expected:
            failures.append(f"sha256 mismatch: {relative}: {observed} != {expected}")
        else:
            print(f"[OK] {relative}")

    forbidden = [
        ROOT / "src/kmlee_bam/analysis",
        ROOT / "scripts/posthoc",
    ]
    for path in forbidden:
        if path.exists():
            failures.append(f"posthoc/analysis code should be absent: {path.relative_to(ROOT)}")

    sys.path.insert(0, str(ROOT / "src"))
    try:
        module = importlib.import_module("kmlee_bam.training.run_current")
        print(f"[OK] import {module.__name__}")
    except Exception as exc:
        failures.append(f"training import failed: {type(exc).__name__}: {exc}")

    if failures:
        for failure in failures:
            print(f"[FAIL] {failure}", file=sys.stderr)
        raise SystemExit(1)
    print("[PASS] CLS PRISM rank2 epoch-20 training bundle")


if __name__ == "__main__":
    main()
