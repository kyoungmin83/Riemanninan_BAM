#!/usr/bin/env bash
set -Eeuo pipefail

RUNTIME="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_20260824"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
CONFIG="$RUNTIME/configs/generated/train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
ARTIFACT="$RUNTIME/artifacts/prism_module_rescue_stats_train64_pathaxis_v3_20260824.npz"
PREFLIGHT="$RUNTIME/artifacts/preflight_prism_integrated_rank8_primary_20260824.json"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_retry1_20260824"
LOG="$RUNTIME/logs/prism_integrated_rank8_primary_sv6_retry1_20260824.log"
EXPECTED_CONFIG_SHA="6a1169c96e3182139812d1fc77314f4bf4bf2f0ffd40dca57b40231c6d995a86"
EXPECTED_ARTIFACT_SHA="b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c"
PORT="${PRISM_INTEGRATED_RANK8_PORT:-29984}"

mkdir -p "$RUNTIME/logs"
cd "$RUNTIME"

if [[ -e "$RUN" ]]; then
  echo "FATAL: output path already exists: $RUN" >&2
  exit 2
fi
if [[ -e "$LOG" ]]; then
  echo "FATAL: launch log already exists: $LOG" >&2
  exit 2
fi

config_sha="$(sha256sum "$CONFIG" | awk '{print $1}')"
artifact_sha="$(sha256sum "$ARTIFACT" | awk '{print $1}')"
if [[ "$config_sha" != "$EXPECTED_CONFIG_SHA" ]]; then
  echo "FATAL: config SHA mismatch: $config_sha" >&2
  exit 2
fi
if [[ "$artifact_sha" != "$EXPECTED_ARTIFACT_SHA" ]]; then
  echo "FATAL: train-only pathology artifact SHA mismatch: $artifact_sha" >&2
  exit 2
fi

PYTHONPATH="$RUNTIME/src" "$PYTHON" - "$CONFIG" "$PREFLIGHT" <<'PY'
import json
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
receipt = json.load(open(sys.argv[2], encoding="utf-8"))
gate = config["learned_generator_count"]
align = config["v7a"]["pathology_aware_uncertainty"]["alignment"]
assert receipt["status"] == "PASS" and receipt["official_test_dataset_opened"] is False
assert config["_launch_guard"]["launch_allowed"] is True
assert config["decoder"]["pathology_rank"] == 8
assert config["precision_medicine"]["personal_rank"] == 2
assert (gate["start_epoch"], gate["shadow_end_epoch"], gate["soft_end_epoch"]) == (13, 16, 24)
assert gate["minimum_active_generators"] is None
assert align["normalize_relative_within_celltype"] is True
assert align["minimum_relative_ess_fraction"] == 0.9
assert config["train"]["skip_final_test"] is True
print("[launch-contract] PASS rank=8 generator shadow/soft/hard=13/17/25 sealed-test=true")
PY

active_gpu_pids="$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | xargs || true)"
if [[ -n "$active_gpu_pids" ]]; then
  echo "FATAL: active GPU processes: $active_gpu_pids" >&2
  exit 2
fi

export PYTHONPATH="$RUNTIME/src"
export CUDA_VISIBLE_DEVICES="0,1,2,3,4,5"
export OMP_NUM_THREADS="8"
export NCCL_P2P_DISABLE="1"
export NCCL_IB_DISABLE="1"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export TORCH_NCCL_ASYNC_ERROR_HANDLING="1"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_LAUNCH_SCRIPT="$RUNTIME/scripts/generated/launch_prism_integrated_rank8_primary_sv6_20260824.sh"

echo "[$(date --iso-8601=seconds)] START one-lineage Phase-I -> Phase-II rank8" | tee -a "$LOG"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$LOG"
echo "[$(date --iso-8601=seconds)] COMPLETE" | tee -a "$LOG"
