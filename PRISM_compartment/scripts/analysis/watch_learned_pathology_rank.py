#!/usr/bin/env python3
"""Write a compact epoch-by-epoch learned pathology-rank log."""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-jsonl", type=Path, required=True)
    parser.add_argument("--continuation-jsonl", type=Path, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=15.0)
    parser.add_argument("--final-epoch", type=int, default=55)
    return parser.parse_args()


def read_records(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict) and "epoch" in value and "learned_rank" in value:
            records.append(value)
    return records


def compact(record: dict[str, Any]) -> dict[str, Any]:
    learned = record["learned_rank"]
    axes = record.get("axes", [])
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "epoch": int(record["epoch"]),
        "capacity": int(learned["capacity"]),
        "expected_rank": float(learned["expected_rank"]),
        "hard_rank": int(learned["hard_rank"]),
        "minimum_rank": int(learned["minimum_rank"]),
        "finalized": bool(learned["finalized"]),
        "mode": str(learned["mode"]),
        "temperature": float(learned["temperature"]),
        "uncertain_fraction": float(learned["uncertain_fraction"]),
        "axis_participation_ranks": [
            float(axis["participation_rank"]) for axis in axes
        ],
        "axis_rank_95_energy": [int(axis["rank_95_energy"]) for axis in axes],
    }


def main() -> None:
    args = parse_args()
    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if args.output_jsonl.exists() or args.output_tsv.exists():
        raise FileExistsError("refusing to overwrite learned-rank monitor output")

    header = (
        "epoch\tcapacity\texpected_rank\thard_rank\tminimum_rank\tfinalized\tmode"
        "\ttemperature\tuncertain_fraction\taxis_participation_ranks"
        "\taxis_rank_95_energy\ttimestamp_utc\n"
    )
    args.output_tsv.write_text(header, encoding="utf-8")
    seen: set[int] = set()

    while True:
        records_by_epoch: dict[int, dict[str, Any]] = {}
        for path in (args.source_jsonl, args.continuation_jsonl):
            for record in read_records(path):
                records_by_epoch[int(record["epoch"])] = record
        for epoch in sorted(records_by_epoch):
            if epoch in seen:
                continue
            item = compact(records_by_epoch[epoch])
            with args.output_jsonl.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(item, sort_keys=True) + "\n")
            participation = ",".join(
                f"{value:.6f}" for value in item["axis_participation_ranks"]
            )
            energy = ",".join(str(value) for value in item["axis_rank_95_energy"])
            row = (
                f"{item['epoch']}\t{item['capacity']}\t{item['expected_rank']:.6f}"
                f"\t{item['hard_rank']}\t{item['minimum_rank']}"
                f"\t{str(item['finalized']).lower()}\t{item['mode']}"
                f"\t{item['temperature']:.6f}\t{item['uncertain_fraction']:.6f}"
                f"\t{participation}\t{energy}\t{item['timestamp_utc']}\n"
            )
            with args.output_tsv.open("a", encoding="utf-8") as handle:
                handle.write(row)
            print(
                "[pathology-rank] "
                f"epoch={item['epoch']} capacity={item['capacity']} "
                f"expected={item['expected_rank']:.3f} hard={item['hard_rank']} "
                f"finalized={item['finalized']} mode={item['mode']} "
                f"uncertain={item['uncertain_fraction']:.3f}",
                flush=True,
            )
            seen.add(epoch)
        if seen and max(seen) >= int(args.final_epoch):
            return
        time.sleep(float(args.poll_seconds))


if __name__ == "__main__":
    main()
