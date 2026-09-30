#!/usr/bin/env bash
# Validation-only technical-zero audit. Opens the train split only; test stays sealed.
set -Eeuo pipefail

REPO="/home/kmlee/project_local/kmlee_bam_integrated_pathranklearn_20260824"
SUPPORT="/home/kmlee/project_local/kmlee_bam/scripts"
PY="/home/kmlee/project_local/kmlee_bam/.conda/envs/kmlee/bin/python"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_pathranklearn_s42_sv7_freegen_rankactive_resume_e19_20260826"
SOURCE_CONFIG="$RUN/source_config.json"
RESOLVED_CONFIG="$RUN/resolved_config.json"
CHECKPOINT="$RUN/checkpoint_epoch_028.pt"
SCRIPT="$REPO/scripts/generated/analyze_latent_knn_pi_tech_posthoc.py"
OUT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_latent_knn_pi_tech_sv7_e28_posthoc_20260827"
LOG="$OUT/run.log"

for required in "$PY" "$SOURCE_CONFIG" "$RESOLVED_CONFIG" "$CHECKPOINT" "$SCRIPT"; do
  [[ -r "$required" ]] || { echo "FATAL missing=$required" >&2; exit 2; }
done

mkdir -p "$OUT"
cd "$REPO"
export PYTHONPATH="$REPO/src:$SUPPORT"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export MPLCONFIGDIR="/tmp/mpl_prism_pi_tech_e28_20260827"
export XDG_CACHE_HOME="/tmp/cache_prism_pi_tech_e28_20260827"

exec nice -n 15 "$PY" "$SCRIPT" \
  --source-config "$SOURCE_CONFIG" \
  --resolved-config "$RESOLVED_CONFIG" \
  --checkpoint "$CHECKPOINT" \
  --outdir "$OUT" \
  --support-scripts "$SUPPORT" \
  --device cuda:0 \
  --epoch 28 \
  --world-size 4 \
  --batch-size 48 \
  --block-size 8 \
  --batches-per-rank 8 \
  --bank-per-donor-celltype 2 \
  --bank-loader-batch-size 24 \
  --model-chunk-size 8 \
  --num-workers 0 \
  --seed 20260827 >"$LOG" 2>&1
