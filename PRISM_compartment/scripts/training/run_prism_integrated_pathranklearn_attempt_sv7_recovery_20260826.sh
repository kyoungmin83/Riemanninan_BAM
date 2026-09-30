#!/usr/bin/env bash
set -Eeuo pipefail

RUNTIME="/home/kmlee/project_local/kmlee_bam_integrated_pathranklearn_20260824"
PYTHON="/home/kmlee/project_local/kmlee_bam/.conda/envs/kmlee/bin/python"
BASE_CONFIG="$RUNTIME/configs/generated/train_config_prism_integrated_pathranklearn_s42_sv7_cls_recovery_20260826.json"
ACTIVE_CONFIG="$RUNTIME/configs/runtime/train_config_cls_recovery_active.json"
ARTIFACT="$RUNTIME/artifacts/prism_module_rescue_stats_train64_pathaxis_v3_20260824.npz"
PREFLIGHT="$RUNTIME/artifacts/preflight_prism_integrated_pathranklearn_cls_recovery_20260826.json"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_pathranklearn_s42_sv7_resume_e16_20260826"
LOG="$RUNTIME/logs/prism_integrated_pathranklearn_cls_recovery_sv7_20260826.log"
LEDGER="$RUNTIME/logs/unattended_cls_recovery_restart_ledger_20260826.jsonl"
COMPLETE="$RUNTIME/artifacts/TRAINING_COMPLETE_CLS_RECOVERY_20260826.json"
EXPECTED_CONFIG_SHA="ed8c248e79580fdadad65b358ac8c8a2f5b9e014eac364c749bd83f7652ec917"
EXPECTED_ARTIFACT_SHA="b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c"
PORT="${PRISM_INTEGRATED_PATHRANK_PORT:-29987}"
NVIDIA_COMPAT_ROOT="/home/kmlee/.local/nvidia-compat-580.159.03/root-complete"

mkdir -p "$RUNTIME/logs" "$RUNTIME/configs/runtime" "$RUNTIME/artifacts"
cd "$RUNTIME"

if [[ -f "$COMPLETE" ]]; then
  echo "[$(date --iso-8601=seconds)] already complete: $COMPLETE"
  exit 0
fi

config_sha="$(sha256sum "$BASE_CONFIG" | awk '{print $1}')"
artifact_sha="$(sha256sum "$ARTIFACT" | awk '{print $1}')"
if [[ "$config_sha" != "$EXPECTED_CONFIG_SHA" ]]; then
  echo "FATAL: immutable config SHA mismatch: $config_sha" >&2
  exit 2
fi
if [[ "$artifact_sha" != "$EXPECTED_ARTIFACT_SHA" ]]; then
  echo "FATAL: train-only pathology artifact SHA mismatch: $artifact_sha" >&2
  exit 2
fi

"$PYTHON" -c 'import json,sys; r=json.load(open(sys.argv[1])); assert r["status"]=="PASS" and r["source_epoch"]==16 and r["next_epoch"]==17 and r["official_test_dataset_opened"] is False and r["agp_only_state_readout_rescue_enabled"] is False and r["learned_rank_state_present"] is True' "$PREFLIGHT"

PYTHONPATH="$RUNTIME/src" "$PYTHON" \
  "$RUNTIME/scripts/generated/prepare_pathrank_cls_recovery_resume_config_20260826.py" \
  --base "$BASE_CONFIG" \
  --active "$ACTIVE_CONFIG" \
  --ledger "$LEDGER"

source_epoch="$("$PYTHON" -c 'import json,sys; x=json.load(open(sys.argv[1])); print(x["experiment_manifest"].get("active_resume_source_epoch") or 0)' "$ACTIVE_CONFIG")"
maximum_epoch="$("$PYTHON" -c 'import json,sys; x=json.load(open(sys.argv[1])); print(x["train"]["epochs"])' "$ACTIVE_CONFIG")"
if (( source_epoch >= maximum_epoch )); then
  "$PYTHON" -c 'import datetime,json,sys; json.dump({"status":"COMPLETE_RECOVERED_FROM_FINAL_EPOCH_CHECKPOINT","epoch":int(sys.argv[2]),"timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()},open(sys.argv[1],"w"),indent=2)' "$COMPLETE" "$source_epoch"
  echo "[$(date --iso-8601=seconds)] recovered completion from epoch $source_epoch" | tee -a "$LOG"
  exit 0
fi

export LD_LIBRARY_PATH="$NVIDIA_COMPAT_ROOT/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
active_gpu_pids="$("$NVIDIA_COMPAT_ROOT/usr/bin/nvidia-smi" --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | xargs || true)"
if [[ -n "$active_gpu_pids" ]]; then
  echo "[$(date --iso-8601=seconds)] GPUs busy with PIDs $active_gpu_pids; systemd will retry" | tee -a "$LOG" >&2
  exit 75
fi

export PYTHONPATH="$RUNTIME/src"
export CUDA_VISIBLE_DEVICES="0,1,2,3"
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
export KMLEE_LAUNCH_SCRIPT="$RUNTIME/scripts/generated/run_prism_integrated_pathranklearn_attempt_sv7_recovery_20260826.sh"

echo "[$(date --iso-8601=seconds)] START integrated learned-pathology-rank source_epoch=$source_epoch" | tee -a "$LOG"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=4 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT" \
  -m kmlee_bam.training.run_current --config "$ACTIVE_CONFIG" \
  2>&1 | tee -a "$LOG"

"$PYTHON" -c 'import datetime,json,sys; json.dump({"status":"COMPLETE","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()},open(sys.argv[1],"w"),indent=2)' "$COMPLETE"
echo "[$(date --iso-8601=seconds)] COMPLETE" | tee -a "$LOG"
