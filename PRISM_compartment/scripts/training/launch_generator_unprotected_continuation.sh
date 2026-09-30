#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "$#" -ne 10 ]]; then
  echo "usage: $0 RUNTIME PYTHON CONFIG CONFIG_SHA CHECKPOINT CHECKPOINT_SHA RECEIPT LOG NPROC PORT" >&2
  exit 2
fi

RUNTIME="$1"
PYTHON="$2"
CONFIG="$3"
EXPECTED_CONFIG_SHA="$4"
CHECKPOINT="$5"
EXPECTED_CHECKPOINT_SHA="$6"
RECEIPT="$7"
LOG="$8"
NPROC="$9"
PORT="${10}"

if [[ -e "$LOG" ]]; then
  echo "FATAL: launch log already exists: $LOG" >&2
  exit 2
fi

config_sha="$(sha256sum "$CONFIG" | awk '{print $1}')"
checkpoint_sha="$(sha256sum "$CHECKPOINT" | awk '{print $1}')"
if [[ "$config_sha" != "$EXPECTED_CONFIG_SHA" ]]; then
  echo "FATAL: config SHA mismatch: $config_sha" >&2
  exit 2
fi
if [[ "$checkpoint_sha" != "$EXPECTED_CHECKPOINT_SHA" ]]; then
  echo "FATAL: checkpoint SHA mismatch: $checkpoint_sha" >&2
  exit 2
fi

PYTHONPATH="$RUNTIME/src" "$PYTHON" - "$CONFIG" "$RECEIPT" "$CHECKPOINT" <<'PY'
import json
import os
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
receipt = json.load(open(sys.argv[2], encoding="utf-8"))
checkpoint = os.path.realpath(sys.argv[3])
train = config["train"]
gate = config["learned_generator_count"]
guard = config["_launch_guard"]
assert receipt["status"] == "PASS"
assert receipt["output_checkpoint"] == checkpoint
assert receipt["new_protected_count"] == 10
assert receipt["checks"]["log_alpha_bitwise_preserved"] is True
assert receipt["checks"]["optimizer_state_present"] is True
assert receipt["checks"]["scheduler_state_present"] is True
assert receipt["checks"]["curriculum_state_present"] is True
assert train["resume_checkpoint"] == checkpoint
assert train["resume_optimizer"] is True
assert train["resume_history"] is True
assert not train.get("allow_in_place_resume", False)
assert train["skip_final_test"] is True
assert not os.path.exists(train["out_dir"])
assert config["precision_medicine"]["personal_rank"] == 2
assert gate["enabled"] is True
assert gate["protect_unique_gene_coverage"] is False
assert gate["protect_singletons"] is True
assert gate["minimum_active_generators"] is None
assert guard["authorized_design_answer"] == "integrated"
assert guard["official_test_dataset_open_forbidden"] is True
assert guard["expected_generator_protected_count"] == 10
print(
    "[launch-contract] PASS integrated continuation; "
    f"epoch={receipt['source_epoch']}->{receipt['source_epoch'] + 1}; "
    "generator protection=singleton-only(10); official test sealed"
)
PY

if [[ -x /home/kmlee/.local/nvidia-compat-580.159.03/root-complete/usr/bin/nvidia-smi ]]; then
  NVIDIA_COMPAT_ROOT=/home/kmlee/.local/nvidia-compat-580.159.03/root-complete
  export LD_LIBRARY_PATH="$NVIDIA_COMPAT_ROOT/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  NVIDIA_SMI="$NVIDIA_COMPAT_ROOT/usr/bin/nvidia-smi"
else
  NVIDIA_SMI=nvidia-smi
fi
active_gpu_pids="$("$NVIDIA_SMI" --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | xargs || true)"
if [[ -n "$active_gpu_pids" ]]; then
  echo "FATAL: active GPU processes: $active_gpu_pids" >&2
  exit 2
fi

mkdir -p "$(dirname "$LOG")"
cd "$RUNTIME"
export PYTHONPATH="$RUNTIME/src"
if [[ "$NPROC" -eq 6 ]]; then
  export CUDA_VISIBLE_DEVICES="0,1,2,3,4,5"
elif [[ "$NPROC" -eq 4 ]]; then
  export CUDA_VISIBLE_DEVICES="0,1,2,3"
else
  echo "FATAL: unsupported NPROC=$NPROC" >&2
  exit 2
fi
export OMP_NUM_THREADS=8
export NCCL_P2P_DISABLE=1
export NCCL_IB_DISABLE=1
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ASYNC_ERROR_HANDLING=1
export KMLEE_TRAIN_BOOT_DEBUG=1
export KMLEE_NONFINITE_GRAD_DEBUG=1
export KMLEE_CONSOLE_LOG_STYLE=prism_informative
export KMLEE_GROUPED_BLOCK_SIZE=8
export KMLEE_PATH_CONSISTENCY=0
export KMLEE_LAUNCH_SCRIPT="$0"

echo "[$(date --iso-8601=seconds)] START integrated free-generator continuation" | tee -a "$LOG"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node="$NPROC" \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$LOG"
echo "[$(date --iso-8601=seconds)] COMPLETE" | tee -a "$LOG"
