#!/usr/bin/env python3
"""Verify the published E52 source copy without starting a training run."""
from __future__ import annotations
import argparse
import ast
import hashlib
import importlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
E52_SHA256 = "9fbb7df5119aba87a52d437fe443d656044e3c9c868fe13b1c39b73ac375944d"

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--import-check", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "FILE_MANIFEST.json").read_text())
    failures = []
    for relative, record in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        data = path.read_bytes()
        if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
            failures.append(f"hash/size mismatch: {relative}")
        if path.suffix == ".py":
            try:
                ast.parse(data, filename=relative)
            except SyntaxError as exc:
                failures.append(f"syntax: {relative}: {exc}")
    checkpoint = json.loads((ROOT / "model_card/CHECKPOINT_MANIFEST.json").read_text())
    if checkpoint["epoch"] != 52 or checkpoint["checkpoint_sha256"] != E52_SHA256:
        failures.append("checkpoint identity does not match compartment E52")
    if args.import_check:
        sys.path.insert(0, str(ROOT / "src"))
        for name in ("kmlee_bam.model.precision_medicine", "kmlee_bam.model.system", "kmlee_bam.training.run_current"):
            try:
                importlib.import_module(name)
                print(f"[OK] import {name}")
            except Exception as exc:
                failures.append(f"import {name}: {type(exc).__name__}: {exc}")
    if failures:
        for failure in failures:
            print(f"[FAIL] {failure}", file=sys.stderr)
        raise SystemExit(1)
    print(f"[PASS] compartment E52 source bundle: {len(manifest['files'])} hashed files; no training or data loading")

if __name__ == "__main__":
    main()
