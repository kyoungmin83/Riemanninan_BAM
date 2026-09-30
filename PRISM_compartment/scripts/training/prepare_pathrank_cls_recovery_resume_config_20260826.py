#!/usr/bin/env python3
"""Create the active config from the immutable launch config and latest epoch."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import torch


EPOCH_PATTERN = re.compile(r"^checkpoint_epoch_(\d+)\.pt$")


def latest_epoch_checkpoint(out_dir: Path) -> tuple[int, Path] | None:
    candidates: list[tuple[int, Path]] = []
    if out_dir.is_dir():
        for path in out_dir.iterdir():
            match = EPOCH_PATTERN.match(path.name)
            if match and path.is_file() and path.stat().st_size > 0:
                candidates.append((int(match.group(1)), path.resolve()))
        last = out_dir / "checkpoint_last.pt"
        if last.is_file() and last.stat().st_size > 0:
            candidates.append((10**9, last.resolve()))
    for _, path in sorted(candidates, reverse=True, key=lambda item: item[0]):
        try:
            payload = torch.load(
                path,
                map_location="cpu",
                weights_only=False,
                mmap=True,
            )
            epoch = int(payload.get("epoch", 0))
            if epoch >= 1 and payload.get("step_in_epoch") is None:
                return epoch, path
        except Exception:
            continue
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--active", required=True, type=Path)
    parser.add_argument("--ledger", required=True, type=Path)
    args = parser.parse_args()

    config = json.loads(args.base.read_text(encoding="utf-8"))
    guard = config.get("_launch_guard", {})
    if not bool(guard.get("unattended_same_design_restart_authorized")):
        raise RuntimeError("unattended restart authorization is absent")
    if int(config["decoder"]["pathology_rank"]) != 0:
        raise RuntimeError("fixed pathology rank detected")
    if config["experiment_manifest"].get("pathology_rank_target") is not None:
        raise RuntimeError("pathology rank target detected")

    out_dir = Path(config["train"]["out_dir"])
    latest = latest_epoch_checkpoint(out_dir)
    source_epoch = None
    resume_path = config["train"].get("resume_checkpoint")
    if latest is not None:
        source_epoch, checkpoint = latest
        resume_path = str(checkpoint)
    elif resume_path:
        checkpoint = Path(resume_path)
        payload = torch.load(
            checkpoint,
            map_location="cpu",
            weights_only=False,
            mmap=True,
        )
        source_epoch = int(payload.get("epoch", 0))
        if source_epoch < 1 or payload.get("step_in_epoch") is not None:
            raise RuntimeError("configured bootstrap checkpoint is not an epoch boundary")
    config["train"]["resume_checkpoint"] = resume_path
    config["experiment_manifest"]["active_resume_source_epoch"] = source_epoch
    config["experiment_manifest"]["active_resume_checkpoint"] = resume_path
    args.active.parent.mkdir(parents=True, exist_ok=True)
    args.active.write_text(
        json.dumps(config, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    attempt = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_epoch": source_epoch,
        "resume_checkpoint": resume_path,
        "active_config": str(args.active.resolve()),
    }
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    with args.ledger.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(attempt, sort_keys=True) + "\n")
    print(json.dumps(attempt, sort_keys=True))


if __name__ == "__main__":
    main()
