#!/usr/bin/env bash
set -Eeuo pipefail

RUNTIME="/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
BASE_CONFIG="$RUNTIME/configs/final/train_config_prism_integrated_rank8_primary_s42_sv6_retry1_20260824.json"
PHASE1_CONFIG="$RUNTIME/configs/phase1/train_config_kmlee_bam_dlpfc_mtg_bins4_v31a_cons_h8_ct64_gen414_singletonr2_s42.json"
CONTEXT="/home/kmlee/project_sv6/kmlee_bam/outputs/population_fingerprint_20260714/ctemb64_epoch12_full_region_pseudobulk.npz"
REGISTRY="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/final_gene_module_registry_compact.json"
ARTIFACT_DIR="$RUNTIME/artifacts"
LOG_DIR="$RUNTIME/logs"
CAP="$ARTIFACT_DIR/rank2_e20_train64_output_cap_q995_20260825.npz"
RESCUE="$ARTIFACT_DIR/prism_module_rescue_stats_train64_pathaxis_v3_20260824.npz"
GRAPH="$ARTIFACT_DIR/module_local_compartment_graph_registry_only_top4_j005_20260828.npz"
HALVES="$ARTIFACT_DIR/module_local_split_half_train64_s42_20260828.npz"
RELIABILITY="$ARTIFACT_DIR/module_local_reliability_train64_s42_20260828.npz"
FINAL_CONFIG="$RUNTIME/configs/generated/train_config_prism_integrated_compartmental_pathrank_s42_sv6_20260828.json"
PREFLIGHT="$ARTIFACT_DIR/preflight_prism_integrated_compartmental_pathrank_20260828.json"
LAUNCH_RECEIPT="$ARTIFACT_DIR/launch_receipt_prism_integrated_compartmental_pathrank_20260828.json"
SMOKE_FIRST="$RUNTIME/configs/generated/smoke_first_20260828.json"
SMOKE_RESUME="$RUNTIME/configs/generated/smoke_resume_20260828.json"
SMOKE_RUN="$RUNTIME/smoke/ddp_resume"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_20260828"
SUPERVISOR_LOG="$LOG_DIR/prepare_and_launch_20260828.log"
TRAIN_LOG="$LOG_DIR/train_integrated_compartmental_pathrank_20260828.log"
PORT_SMOKE_1="${PRISM_COMP_SMOKE_PORT1:-29881}"
PORT_SMOKE_2="${PRISM_COMP_SMOKE_PORT2:-29882}"
PORT_MAIN="${PRISM_COMP_MAIN_PORT:-29883}"

mkdir -p "$ARTIFACT_DIR" "$LOG_DIR" "$RUNTIME/configs/generated" "$RUNTIME/smoke"
cd "$RUNTIME"
exec > >(tee -a "$SUPERVISOR_LOG") 2>&1

echo "[$(date --iso-8601=seconds)] PREPARE integrated compartmental pathrank on SV6"
if [[ -e "$RUN" ]]; then
  echo "FATAL: main output already exists: $RUN" >&2
  exit 2
fi
if [[ ! -f "$CAP" || ! -f "$RESCUE" ]]; then
  echo "FATAL: pre-staged output-cap/module-rescue artifacts are missing" >&2
  exit 2
fi
if [[ "$(sha256sum "$CAP" | awk '{print $1}')" != "76edbc17c78a29318302468fc0cf045c88cfe9b58b147f2ec57af673bcd8666d" ]]; then
  echo "FATAL: output-cap artifact SHA mismatch" >&2
  exit 2
fi
if [[ "$(sha256sum "$RESCUE" | awk '{print $1}')" != "b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c" ]]; then
  echo "FATAL: module-rescue artifact SHA mismatch" >&2
  exit 2
fi

export PYTHONPATH="$RUNTIME/src:$RUNTIME"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export OMP_NUM_THREADS="8"

if [[ "${PRISM_REUSE_SEALED_ARTIFACTS:-0}" == "1" ]]; then
  echo "[$(date --iso-8601=seconds)] REUSE completed sealed artifacts after preflight-only retry"
  test -f "$GRAPH"
  test -f "$HALVES"
  test -f "$RELIABILITY"
else
  echo "[$(date --iso-8601=seconds)] BUILD registry-only compartment graph"
  "$PYTHON" scripts/training/build_prism_module_local_compartment_graph_20260827.py \
    --registry "$REGISTRY" --output "$GRAPH" --topk 4 --minimum-jaccard 0.05

  echo "[$(date --iso-8601=seconds)] BUILD train-only split-half module summaries on 6 GPUs"
  split_pids=()
  for shard in 0 1 2 3 4 5; do
    CUDA_VISIBLE_DEVICES="$shard" "$PYTHON" \
      scripts/training/build_prism_module_local_split_half_input_20260828.py shard \
      --config "$BASE_CONFIG" --source-context "$CONTEXT" \
      --output "$ARTIFACT_DIR/module_local_split_half_shard_${shard}.npz" \
      --shard-index "$shard" --num-shards 6 --seed 42 --batch-size 256 \
      --device cuda:0 \
      >"$LOG_DIR/module_local_split_half_shard_${shard}.log" 2>&1 &
    split_pids+=("$!")
  done
  for index in 0 1 2 3 4 5; do
    wait "${split_pids[$index]}"
    echo "  split-half shard $index complete"
  done
  "$PYTHON" scripts/training/build_prism_module_local_split_half_input_20260828.py merge \
    --inputs "$ARTIFACT_DIR"/module_local_split_half_shard_*.npz \
    --output "$HALVES" --minimum-cells-per-half 5

  echo "[$(date --iso-8601=seconds)] BUILD sealed reliability artifact"
  "$PYTHON" scripts/training/build_prism_module_local_reliability_20260825.py \
    --input "$HALVES" --output "$RELIABILITY" --source-context "$CONTEXT" \
    --registry "$REGISTRY" \
    --activity-dictionary "/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/activity_weight_kme_or_l2_membership.npz" \
    --seed 42 --minimum-cells-per-half 5 --minimum-donors 12 \
    --minimum-reliability 0.20
fi

RELIABILITY_SHA="$(sha256sum "$RELIABILITY" | awk '{print $1}')"
GRAPH_SHA="$(sha256sum "$GRAPH" | awk '{print $1}')"
CAP_SHA="$(sha256sum "$CAP" | awk '{print $1}')"
RESCUE_SHA="$(sha256sum "$RESCUE" | awk '{print $1}')"

echo "[$(date --iso-8601=seconds)] BUILD exact executable config"
"$PYTHON" scripts/training/build_prism_integrated_compartmental_pathrank_config_20260828.py \
  --output "$FINAL_CONFIG" \
  --reliability "$RELIABILITY" --reliability-sha256 "$RELIABILITY_SHA" \
  --compartment-graph "$GRAPH" --compartment-graph-sha256 "$GRAPH_SHA" \
  --output-cap "$CAP" --output-cap-sha256 "$CAP_SHA" \
  --module-rescue "$RESCUE" --module-rescue-sha256 "$RESCUE_SHA"

if [[ "${PRISM_REUSE_PASSED_PREFLIGHT:-0}" == "1" ]]; then
  echo "[$(date --iso-8601=seconds)] REUSE PASS preflight after smoke-only efficiency retry"
  test -f "$PREFLIGHT"
  grep -q '"status": "PASS"' "$PREFLIGHT"
  grep -q '"official_test_dataset_opened": false' "$PREFLIGHT"
else
  echo "[$(date --iso-8601=seconds)] RUN targeted regression gates"
  for test_pattern in \
    test_integrated_learned_pathology_rank_curriculum.py \
    test_module_local_compartmental_nonlinearity.py \
    test_latent_knn_pi_tech_bank.py \
    test_architecture_capacity_logging.py; do
    "$PYTHON" -m unittest discover -s tests -p "$test_pattern" -q
  done

  echo "[$(date --iso-8601=seconds)] RUN full-data Phase-I parity/schedule/budget/artifact preflight"
  CUDA_VISIBLE_DEVICES="" "$PYTHON" \
    scripts/training/preflight_prism_integrated_compartmental_pathrank_20260828.py \
    --phase1 "$PHASE1_CONFIG" --integrated "$FINAL_CONFIG" --report "$PREFLIGHT"
fi

echo "[$(date --iso-8601=seconds)] RUN real production-world-size 6-GPU checkpoint/resume smoke"
if [[ -e "$SMOKE_RUN" ]]; then
  echo "FATAL: smoke output already exists: $SMOKE_RUN" >&2
  exit 2
fi
"$PYTHON" scripts/training/build_prism_integrated_compartmental_pathrank_smoke_configs_20260828.py \
  --source "$FINAL_CONFIG" --first "$SMOKE_FIRST" --resume "$SMOKE_RESUME" \
  --smoke-run "$SMOKE_RUN"

export NCCL_P2P_DISABLE="1"
export NCCL_IB_DISABLE="1"
export TORCH_NCCL_ASYNC_ERROR_HANDLING="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_CONSOLE_LOG_STYLE="prism_informative"
CUDA_VISIBLE_DEVICES="0,1,2,3,4,5" "$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 --rdzv-backend=c10d \
  --rdzv-endpoint="localhost:$PORT_SMOKE_1" \
  -m kmlee_bam.training.run_current --config "$SMOKE_FIRST" \
  2>&1 | tee "$LOG_DIR/ddp_smoke_first_20260828.log"
test -f "$SMOKE_RUN/checkpoint_epoch_001.pt"
CUDA_VISIBLE_DEVICES="0,1,2,3,4,5" "$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 --rdzv-backend=c10d \
  --rdzv-endpoint="localhost:$PORT_SMOKE_2" \
  -m kmlee_bam.training.run_current --config "$SMOKE_RESUME" \
  2>&1 | tee "$LOG_DIR/ddp_smoke_resume_20260828.log"
test -f "$SMOKE_RUN/checkpoint_epoch_002.pt"

echo "[$(date --iso-8601=seconds)] AUTHORIZE SHA-pinned launch"
"$PYTHON" scripts/training/authorize_prism_integrated_compartmental_pathrank_launch_20260828.py \
  --config "$FINAL_CONFIG" --preflight "$PREFLIGHT" \
  --smoke-checkpoint "$SMOKE_RUN/checkpoint_epoch_002.pt" \
  --output "$LAUNCH_RECEIPT"

active_gpu_pids="$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | xargs || true)"
if [[ -n "$active_gpu_pids" ]]; then
  echo "FATAL: GPU processes remain after preflight: $active_gpu_pids" >&2
  exit 2
fi
if [[ -e "$RUN" || -e "$TRAIN_LOG" ]]; then
  echo "FATAL: main output/log appeared before launch" >&2
  exit 2
fi

echo "[$(date --iso-8601=seconds)] START main 6-GPU integrated training" | tee -a "$TRAIN_LOG"
export CUDA_VISIBLE_DEVICES="0,1,2,3,4,5"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_LAUNCH_SCRIPT="$RUNTIME/scripts/training/prepare_and_launch_prism_integrated_compartmental_pathrank_sv6_20260828.sh"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 --rdzv-backend=c10d \
  --rdzv-endpoint="localhost:$PORT_MAIN" \
  -m kmlee_bam.training.run_current --config "$FINAL_CONFIG" \
  2>&1 | tee -a "$TRAIN_LOG"
echo "[$(date --iso-8601=seconds)] COMPLETE main training" | tee -a "$TRAIN_LOG"
