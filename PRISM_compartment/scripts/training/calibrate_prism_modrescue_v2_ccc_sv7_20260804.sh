#!/usr/bin/env bash
set -euo pipefail

REPO="/home/kmlee/project_local/kmlee_bam"
PYTHON="$REPO/.conda/envs/kmlee/bin/python"
CONFIG="$REPO/configs/_calibration_train_config_kmlee_bam_dlpfc_mtg_bins4_prism_modrescue_v2_ccc_nointeraction_s42_sv7.json"
CALIBRATION_DIR="$REPO/outputs/prism_modrescue_v2_ccc_calibration_20260804"
CALIBRATION_OUT="$CALIBRATION_DIR/calibration_sv7.json"
LOG="$REPO/logs/prism_modrescue_v2_ccc_calibration_sv7_20260804_r3.log"
TEMP_RUN="/tmp/prism_modrescue_v2_ccc_calibration_sv7_r3"
COMPAT_LIB="/home/kmlee/.local/nvidia-compat-580.159.03/root-complete/usr/lib/x86_64-linux-gnu"

for required in "$PYTHON" "$CONFIG"; do
  if [[ ! -f "$required" ]]; then
    echo "FATAL missing required file: $required" >&2
    exit 2
  fi
done
if [[ -e "$CALIBRATION_OUT" || -e "$TEMP_RUN" ]]; then
  echo "FATAL prior calibration state exists; refusing implicit overwrite" >&2
  exit 2
fi

mkdir -p "$CALIBRATION_DIR" "$REPO/logs"
cd "$REPO"

export LD_LIBRARY_PATH="$COMPAT_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$REPO/src"
export CUDA_VISIBLE_DEVICES="0"
export OMP_NUM_THREADS="8"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_PATH_CONS_MINCELLS="2"
export KMLEE_LAUNCH_SCRIPT="$REPO/scripts/calibrate_prism_modrescue_v2_ccc_sv7_20260804.sh"
export KMLEE_PRISM_CCC_CALIBRATION_OUT="$CALIBRATION_OUT"

{
  echo "[$(date --iso-8601=seconds)] sv7 train-only CCC gradient calibration begins"
  sha256sum "$CONFIG"
  "$PYTHON" -m kmlee_bam.training.run_current --config "$CONFIG"
  sha256sum "$CALIBRATION_OUT"
  echo "[$(date --iso-8601=seconds)] CALIBRATION_COMPLETE"
} 2>&1 | tee "$LOG"
