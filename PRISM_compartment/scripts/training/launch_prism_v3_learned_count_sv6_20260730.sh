#!/usr/bin/env bash
set -euo pipefail

REPO="/home/kmlee/project_sv6/kmlee_bam"
CONFIG="$REPO/configs/train_config_kmlee_bam_dlpfc_mtg_bins4_prism_v3_learnedN_modrescue_nointeraction_s42_sv6.json"
SMOKE_CONFIG="$REPO/configs/_smoke_prism_v3_learnedN_modrescue_nointeraction_s42_sv6.json"
RUN="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_v3_learnedN_modrescue_nointeraction_s42_20260730"
SMOKE_RUN="/tmp/prism_v3_learnedN_smoke_nointeraction_sv6"
STATS="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/prism_module_rescue_stats_train64_v2.npz"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
PORT_SMOKE_TRAIN="${PRISM_V3_SV6_SMOKE_TRAIN_PORT:-29720}"
PORT_SMOKE_GATE="${PRISM_V3_SV6_SMOKE_GATE_PORT:-29721}"
PORT_FULL_TRAIN="${PRISM_V3_SV6_FULL_TRAIN_PORT:-29722}"
PORT_FULL_GATE="${PRISM_V3_SV6_FULL_GATE_PORT:-29723}"

cd "$REPO"
mkdir -p "$RUN"
PIPELINE_LOG="$RUN/pipeline.log"

export PYTHONPATH="$REPO/src"
export CUDA_VISIBLE_DEVICES="0,1,2,3,4,5"
export OMP_NUM_THREADS="8"
export NCCL_P2P_DISABLE="1"
export NCCL_IB_DISABLE="1"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_PATH_CONS_MINCELLS="2"

log() {
  echo "[$(date --iso-8601=seconds)] $*" | tee -a "$PIPELINE_LOG"
}

log "PRISM-v3 sv6 pipeline: module-rescue ON, learned generator count ON, pathology interaction OFF"
log "phase0=epochs 1..12 on 52 weight donors; gate-only search=epochs 13..40, conditional max 48 on 12 architecture donors"
log "six GPUs; phase0 effective batch=32x6x2=384; gate blocks/epoch=4x6=24"

for required in "$CONFIG" "$SMOKE_CONFIG"; do
  if [[ ! -f "$required" ]]; then
    log "FATAL missing required file: $required"
    exit 2
  fi
done
if [[ -e "$RUN/checkpoint_last.pt" || -e "$RUN/generator_count_search/learned_generator_state.pt" ]]; then
  log "FATAL prior production state exists; refusing an implicit mixed resume"
  exit 2
fi
if [[ -e "$SMOKE_RUN/checkpoint_last.pt" || -e "$SMOKE_RUN/generator_count_search/learned_generator_state.pt" ]]; then
  log "FATAL prior smoke state exists: $SMOKE_RUN"
  exit 2
fi

log "waiting for train-only module-rescue statistics: $STATS"
for _ in $(seq 1 1440); do
  if [[ -s "$STATS" ]]; then
    break
  fi
  if ! tmux has-session -t prismv3stats 2>/dev/null; then
    log "FATAL statistics builder ended without producing the artifact"
    exit 3
  fi
  sleep 30
done
if [[ ! -s "$STATS" ]]; then
  log "FATAL timed out waiting for module-rescue statistics"
  exit 3
fi

log "statistics ready; recording immutable inputs"
sha256sum "$CONFIG" "$SMOKE_CONFIG" "$STATS" | tee -a "$PIPELINE_LOG"
nvidia-smi --query-gpu=index,name,driver_version,memory.used,utilization.gpu \
  --format=csv,noheader | tee -a "$PIPELINE_LOG"

log "distributed smoke phase0 begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_SMOKE_TRAIN" \
  -m kmlee_bam.training.run_current --config "$SMOKE_CONFIG" \
  2>&1 | tee -a "$PIPELINE_LOG"

log "distributed smoke gate-only phase begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_SMOKE_GATE" \
  -m kmlee_bam.training.learned_generator_search_runner \
  --config "$SMOKE_CONFIG" \
  --checkpoint "$SMOKE_RUN/checkpoint_last.pt" \
  --output-dir "$SMOKE_RUN/generator_count_search" \
  2>&1 | tee -a "$PIPELINE_LOG"

log "SMOKE_OK; full all-on phase0 begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_FULL_TRAIN" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$PIPELINE_LOG"

log "phase0 epoch12 complete; learned generator-count search begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_FULL_GATE" \
  -m kmlee_bam.training.learned_generator_search_runner \
  --config "$CONFIG" \
  --checkpoint "$RUN/checkpoint_last.pt" \
  --output-dir "$RUN/generator_count_search" \
  2>&1 | tee -a "$PIPELINE_LOG"

log "PIPELINE_COMPLETE; selected mask is discovery-only and requires fixed-mask retraining"
