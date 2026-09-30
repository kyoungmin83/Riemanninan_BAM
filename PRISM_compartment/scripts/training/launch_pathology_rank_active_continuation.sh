#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "$#" -ne 11 ]]; then
  echo "usage: $0 RUNTIME PYTHON CONFIG CONFIG_SHA CHECKPOINT CHECKPOINT_SHA RECEIPT RUNNER_BASE_SHA LOG NPROC PORT" >&2
  exit 2
fi

RUNTIME="$1"
PYTHON="$2"
CONFIG="$3"
EXPECTED_CONFIG_SHA="$4"
CHECKPOINT="$5"
EXPECTED_CHECKPOINT_SHA="$6"
RECEIPT="$7"
EXPECTED_RUNNER_BASE_SHA="$8"
LOG="$9"
NPROC="${10}"
PORT="${11}"
RUNNER_BASE="$RUNTIME/src/kmlee_bam/training/runner_base.py"

if [[ -e "$LOG" ]]; then
  echo "FATAL: launch log already exists: $LOG" >&2
  exit 2
fi

config_sha="$(sha256sum "$CONFIG" | awk '{print $1}')"
checkpoint_sha="$(sha256sum "$CHECKPOINT" | awk '{print $1}')"
runner_base_sha="$(sha256sum "$RUNNER_BASE" | awk '{print $1}')"
if [[ "$config_sha" != "$EXPECTED_CONFIG_SHA" ]]; then
  echo "FATAL: config SHA mismatch: $config_sha" >&2
  exit 2
fi
if [[ "$checkpoint_sha" != "$EXPECTED_CHECKPOINT_SHA" ]]; then
  echo "FATAL: checkpoint SHA mismatch: $checkpoint_sha" >&2
  exit 2
fi
if [[ "$runner_base_sha" != "$EXPECTED_RUNNER_BASE_SHA" ]]; then
  echo "FATAL: runner_base SHA mismatch: $runner_base_sha" >&2
  exit 2
fi
if ! grep -Fq '"[pathology-rank] "' "$RUNNER_BASE"; then
  echo "FATAL: runner lacks pathology-rank epoch logging" >&2
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
generator_gate = config["learned_generator_count"]
rank_gate = config["learned_pathology_rank"]
curriculum = config["integrated_phase_curriculum"]
guard = config["_launch_guard"]
checks = receipt["checks"]

assert receipt["status"] == "PASS"
assert receipt["output_checkpoint"] == checkpoint
assert receipt["rank_finalized_after"] is False
assert receipt["generator_protected_count"] == 10
assert checks["rank_log_alpha_bitwise_preserved"] is True
assert checks["rank_gate_reactivated"] is True
assert checks["phase2_pathology_scale_installed"] is True
assert checks["curriculum_scale_migrated"] is True
assert checks["generator_singleton_only_protection_preserved"] is True
assert checks["optimizer_state_present"] is True
assert checks["scheduler_state_present"] is True
assert checks["history_present"] is True
assert checks["v7a_state_present"] is True
assert checks["v8_state_present"] is True
assert train["resume_checkpoint"] == checkpoint
assert train["resume_optimizer"] is True
assert train["resume_history"] is True
assert not train.get("allow_in_place_resume", False)
assert train["skip_final_test"] is True
assert not os.path.exists(train["out_dir"])
assert config["precision_medicine"]["personal_rank"] == 2
assert config["decoder"]["pathology_rank"] == 0
assert rank_gate["enabled"] is True
assert rank_gate["freeze_epoch"] > train["epochs"]
assert curriculum["phase2_pathology_scale"] > 0.0
assert generator_gate["enabled"] is True
assert generator_gate["protect_unique_gene_coverage"] is False
assert generator_gate["protect_singletons"] is True
assert generator_gate["minimum_active_generators"] is None
assert guard["authorized_design_answer"] == "integrated"
assert guard["official_test_dataset_open_forbidden"] is True
assert guard["pathology_rank_learning_required"] is True
print(
    "[launch-contract] PASS integrated continuation; "
    f"epoch={receipt['source_epoch']}->{receipt['source_epoch'] + 1}; "
    f"pathology rank active capacity={receipt['rank_capacity']} "
    f"expected={receipt['rank_expected_at_migration']:.6f} "
    f"hard={receipt['rank_hard_at_migration']}; "
    f"pathology_scale={receipt['phase2_pathology_scale']:.3f}; "
    "generator protection=singleton-only(10); official test sealed"
)
PY

NVIDIA_COMPAT_ROOT=/home/kmlee/.local/nvidia-compat-580.159.03/root-complete
NVIDIA_SMI="$NVIDIA_COMPAT_ROOT/usr/bin/nvidia-smi"
if [[ ! -x "$NVIDIA_SMI" ]]; then
  echo "FATAL: SV7 NVIDIA compatibility binary unavailable" >&2
  exit 2
fi
export LD_LIBRARY_PATH="$NVIDIA_COMPAT_ROOT/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
active_gpu_pids="$("$NVIDIA_SMI" --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | xargs || true)"
if [[ -n "$active_gpu_pids" ]]; then
  echo "FATAL: active GPU processes: $active_gpu_pids" >&2
  exit 2
fi

mkdir -p "$(dirname "$LOG")"
cd "$RUNTIME"
export PYTHONPATH="$RUNTIME/src"
if [[ "$NPROC" -eq 4 ]]; then
  export CUDA_VISIBLE_DEVICES="0,1,2,3"
else
  echo "FATAL: unsupported SV7 NPROC=$NPROC" >&2
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
export KMLEE_GROUPED_BLOCK_SIZE=8
export KMLEE_PATH_CONSISTENCY=0
export KMLEE_LAUNCH_SCRIPT="$0"

echo "[$(date --iso-8601=seconds)] START integrated active-pathology-rank continuation" | tee -a "$LOG"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node="$NPROC" \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$LOG"
echo "[$(date --iso-8601=seconds)] COMPLETE" | tee -a "$LOG"
