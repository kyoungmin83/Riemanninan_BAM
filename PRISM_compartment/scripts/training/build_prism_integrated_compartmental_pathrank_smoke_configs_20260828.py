#!/usr/bin/env python3
"""Derive fail-closed one-step DDP/resume smoke configs from the final config."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--first", required=True)
    parser.add_argument("--resume", required=True)
    parser.add_argument("--smoke-run", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = json.loads(Path(args.source).read_text(encoding="utf-8"))
    first = copy.deepcopy(source)
    first["train"].update(
        {
            "epochs": 1,
            "out_dir": args.smoke_run,
            "save_best": False,
            "save_every": 1,
            "eval_every": 1,
            "max_train_steps_per_epoch": 1,
            "max_eval_steps": 1,
            "log_every": 1,
            "progress_log_every": 1,
            "eval_log_every": 1,
            "resume_checkpoint": None,
            "skip_final_test": True,
            "reload_best_before_test": False,
        }
    )
    first["data"].update(
        {"num_workers": 0, "persistent_workers": False, "pin_memory": False}
    )
    # The smoke validates construction, DDP, one optimizer step, checkpoint,
    # and resume.  Re-estimating the production 8,192-cell ordinal weights is
    # unrelated and makes a two-step smoke slower than the steps themselves.
    ordinal = first.setdefault("v4", {}).setdefault("ordinal_balance", {})
    ordinal["max_cells_for_bin_stats"] = 64
    ordinal["bin_stats_log_every"] = 65
    first["experiment_manifest"] = {
        **first["experiment_manifest"],
        "run_name": first["experiment_manifest"]["run_name"] + "_ddp_resume_smoke",
        "smoke_only": True,
        "official_test_dataset_opened_by_training": False,
    }
    first["_launch_guard"]["launch_allowed"] = True
    first_path = Path(args.first)
    first_path.parent.mkdir(parents=True, exist_ok=True)
    first_path.write_text(json.dumps(first, indent=2) + "\n", encoding="utf-8")

    resume = copy.deepcopy(first)
    checkpoint = str(Path(args.smoke_run) / "checkpoint_epoch_001.pt")
    resume["train"].update(
        {
            "epochs": 2,
            "resume_checkpoint": checkpoint,
            "resume_optimizer": True,
            "resume_history": True,
            "allow_in_place_resume": True,
        }
    )
    resume["experiment_manifest"]["resume_smoke_source"] = checkpoint
    resume_path = Path(args.resume)
    resume_path.write_text(json.dumps(resume, indent=2) + "\n", encoding="utf-8")
    print(first_path)
    print(resume_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
