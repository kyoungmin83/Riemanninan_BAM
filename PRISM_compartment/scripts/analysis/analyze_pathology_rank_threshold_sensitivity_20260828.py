#!/usr/bin/env python3
"""Read-only threshold sensitivity audit for learned pathology-rank gates."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import torch


RANK_LOGIT_KEY = "decoder.pathology_rank_gate.log_alpha"
RANK_FINALIZED_KEY = "decoder.pathology_rank_gate.finalized_state"
THRESHOLDS = (0.45, 0.475, 0.49, 0.50, 0.51, 0.525, 0.55)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--epochs", type=int, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = []
    for epoch in args.epochs:
        checkpoint = args.run / f"checkpoint_epoch_{epoch:03d}.pt"
        payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
        if int(payload.get("epoch", -1)) != epoch:
            raise RuntimeError(f"epoch mismatch for {checkpoint}")
        state = payload["system_state_dict"]
        probability = torch.sigmoid(state[RANK_LOGIT_KEY].detach().float())
        records.append(
            {
                "epoch": epoch,
                "checkpoint": str(checkpoint),
                "capacity": int(probability.numel()),
                "expected_count": float(probability.sum().item()),
                "finalized": bool(state[RANK_FINALIZED_KEY].item()),
                "counts_by_threshold": {
                    str(threshold): int((probability >= threshold).sum().item())
                    for threshold in THRESHOLDS
                },
                "probability_quantiles": {
                    str(quantile): float(torch.quantile(probability, quantile).item())
                    for quantile in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
                },
                "count_between_0.45_and_0.55": int(
                    ((probability >= 0.45) & (probability <= 0.55)).sum().item()
                ),
            }
        )
        del payload, state, probability
        gc.collect()
    result = {
        "schema_version": "kmlee_bam.pathology_rank_threshold_sensitivity.v1",
        "read_only": True,
        "thresholds": list(THRESHOLDS),
        "records": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
