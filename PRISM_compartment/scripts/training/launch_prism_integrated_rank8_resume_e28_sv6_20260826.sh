#!/usr/bin/env bash
set -Eeuo pipefail

RUNTIME="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
CONFIG="$RUNTIME/configs/generated/train_config_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826.json"
ARTIFACT="$RUNTIME/artifacts/prism_module_rescue_stats_train64_pathaxis_v3_20260824.npz"
PREFLIGHT="$RUNTIME/artifacts/preflight_prism_integrated_rank8_resume_e28_20260826.json"
CHECKPOINT="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_retry1_20260824/checkpoint_epoch_028.pt"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826"
LOG="$RUNTIME/logs/prism_integrated_rank8_primary_resume_e28_retry2_sv6_20260826_attempt2.log"
EXPECTED_CONFIG_SHA="92c066515cff94a33f48e8373e5a8edb96a1c60d3917ed8e795e2a77d2bb6f5e"
EXPECTED_ARTIFACT_SHA="b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c"
EXPECTED_CHECKPOINT_SHA="8ebb2f0644da8bdd985fa11f66dba6719ad7fd20ee6141f54dac71ae064eb108"
PORT="${PRISM_INTEGRATED_RANK8_RESUME_PORT:-29985}"

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
checkpoint_sha="$(sha256sum "$CHECKPOINT" | awk '{print $1}')"
if [[ "$config_sha" != "$EXPECTED_CONFIG_SHA" ]]; then
  echo "FATAL: config SHA mismatch: $config_sha" >&2
  exit 2
fi
if [[ "$artifact_sha" != "$EXPECTED_ARTIFACT_SHA" ]]; then
  echo "FATAL: train-only pathology artifact SHA mismatch: $artifact_sha" >&2
  exit 2
fi
if [[ "$checkpoint_sha" != "$EXPECTED_CHECKPOINT_SHA" ]]; then
  echo "FATAL: epoch-28 checkpoint SHA mismatch: $checkpoint_sha" >&2
  exit 2
fi

PYTHONPATH="$RUNTIME/src" "$PYTHON" - "$CONFIG" "$PREFLIGHT" <<'PY'
import json
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
receipt = json.load(open(sys.argv[2], encoding="utf-8"))
train = config["train"]
guard = config["_launch_guard"]
assert receipt["status"] == "PASS"
assert receipt["source_epoch"] == 28 and receipt["next_epoch"] == 29
assert receipt["official_test_dataset_opened"] is False
assert guard["launch_allowed"] is True
assert guard["explicit_user_reapproval_required"] is False
assert config["encoder"]["pooling"] == "cls"
assert config["prism_module_rescue"]["allow_state_readout_updates"] is False
assert config["decoder"]["pathology_rank"] == 8
assert config["precision_medicine"]["personal_rank"] == 2
assert train["resume_checkpoint"].endswith("checkpoint_epoch_028.pt")
assert train["resume_optimizer"] is True and train["resume_history"] is True
assert train["skip_final_test"] is True
print("[launch-contract] PASS integrated rank8 strict-resume e28->e29 sealed-test=true")
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
export KMLEE_LAUNCH_SCRIPT="$RUNTIME/scripts/generated/launch_prism_integrated_rank8_resume_e28_sv6_20260826.sh"

echo "[$(date --iso-8601=seconds)] START integrated rank8 strict-resume epoch28 -> epoch29" | tee -a "$LOG"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$LOG"
echo "[$(date --iso-8601=seconds)] COMPLETE" | tee -a "$LOG"
