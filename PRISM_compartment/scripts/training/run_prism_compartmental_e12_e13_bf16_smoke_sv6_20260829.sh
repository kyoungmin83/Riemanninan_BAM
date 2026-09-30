#!/usr/bin/env bash
set -Eeuo pipefail

RUNTIME="/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
CONFIG="${PRISM_COMP_RECOVERY_SMOKE_CONFIG:-$RUNTIME/configs/generated/smoke_e12_e13_bf16_recovery_20260829.json}"
RUN="${PRISM_COMP_RECOVERY_SMOKE_RUN:-/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_e12_e13_bf16_smoke_20260829}"
LOG="${PRISM_COMP_RECOVERY_SMOKE_LOG:-$RUNTIME/logs/smoke_e12_e13_bf16_recovery_20260829.log}"
SOURCE_CHECKPOINT="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_20260828/checkpoint_epoch_012.pt"
PORT="${PRISM_COMP_RECOVERY_SMOKE_PORT:-29884}"

if [[ ! -f "$CONFIG" ]]; then
  echo "FATAL: recovery smoke config is missing: $CONFIG" >&2
  exit 2
fi
if [[ -e "$RUN" || -e "$LOG" ]]; then
  echo "FATAL: refusing to mix recovery smoke with existing output" >&2
  exit 2
fi

cd "$RUNTIME"
export PYTHONPATH="$RUNTIME/src:$RUNTIME"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export OMP_NUM_THREADS="8"
export NCCL_P2P_DISABLE="1"
export NCCL_IB_DISABLE="1"
export TORCH_NCCL_ASYNC_ERROR_HANDLING="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_CONSOLE_LOG_STYLE="prism_informative"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_NONFINITE_GRAD_DEBUG="1"

echo "[$(date --iso-8601=seconds)] START epoch12->13 BF16 recovery smoke" | tee "$LOG"
CUDA_VISIBLE_DEVICES="0,1,2,3,4,5" "$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 --rdzv-backend=c10d \
  --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$LOG"

if grep -q "\[skip\] non-finite grad norm" "$LOG"; then
  echo "FATAL: recovery smoke skipped the optimizer update for non-finite gradients" >&2
  exit 3
fi
test -f "$RUN/checkpoint_epoch_013.pt"
"$PYTHON" - "$RUN" "$SOURCE_CHECKPOINT" <<'PY'
import json
import sys
from pathlib import Path

import torch

run = Path(sys.argv[1])
source_checkpoint = Path(sys.argv[2])
checkpoint = torch.load(
    run / "checkpoint_epoch_013.pt", map_location="cpu", weights_only=False
)
source = torch.load(source_checkpoint, map_location="cpu", weights_only=False)
receipt = json.loads((run / "resume_receipt.json").read_text(encoding="utf-8"))
if int(checkpoint.get("epoch", 0)) != 13:
    raise RuntimeError("recovery smoke did not save epoch 13")
if len(checkpoint.get("history") or []) != 13:
    raise RuntimeError("recovery smoke did not preserve and extend history")
if int(receipt.get("source_epoch", 0)) != 12:
    raise RuntimeError("recovery smoke restored the wrong source epoch")
if int(receipt.get("next_epoch", 0)) != 13:
    raise RuntimeError("recovery smoke resumed at the wrong epoch")
if receipt.get("official_test_used") is not False:
    raise RuntimeError("recovery smoke opened the official test split")
if not all(
    bool(receipt.get(key))
    for key in (
        "optimizer_restored",
        "history_restored",
        "v7a_state_restored",
        "v8_state_restored",
        "curriculum_state_restored",
    )
):
    raise RuntimeError("recovery smoke did not restore complete training state")


def optimizer_steps(payload):
    result = []
    for state in payload["optimizer_state_dict"]["state"].values():
        step = state.get("step")
        if step is not None:
            result.append(float(step))
    return result


source_steps = optimizer_steps(source)
smoke_steps = optimizer_steps(checkpoint)
if not source_steps or not smoke_steps or max(smoke_steps) <= max(source_steps):
    raise RuntimeError("recovery smoke did not complete an optimizer update")
print(
    json.dumps(
        {
            "status": "PASS",
            "checkpoint_epoch": int(checkpoint["epoch"]),
            "global_step": int(checkpoint["global_step"]),
            "history_length": len(checkpoint["history"]),
            "source_epoch": int(receipt["source_epoch"]),
            "next_epoch": int(receipt["next_epoch"]),
            "source_optimizer_step_max": max(source_steps),
            "smoke_optimizer_step_max": max(smoke_steps),
            "official_test_used": bool(receipt["official_test_used"]),
        },
        sort_keys=True,
    )
)
PY
echo "[$(date --iso-8601=seconds)] PASS epoch12->13 BF16 recovery smoke" | tee -a "$LOG"
