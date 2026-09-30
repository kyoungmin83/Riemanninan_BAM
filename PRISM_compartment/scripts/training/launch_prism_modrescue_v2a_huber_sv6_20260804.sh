#!/usr/bin/env bash
set -euo pipefail

REPO="/home/kmlee/project_sv6/kmlee_bam"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
CONFIG="$REPO/configs/train_config_kmlee_bam_dlpfc_mtg_bins4_prism_modrescue_v2a_huber_nointeraction_s42_sv6.json"
SMOKE_CONFIG="$REPO/configs/_smoke_train_config_kmlee_bam_dlpfc_mtg_bins4_prism_modrescue_v2a_huber_nointeraction_s42_sv6.json"
RUN="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_modrescue_v2a_huber_nointeraction_s42_20260804"
SMOKE_RUN="/tmp/prism_modrescue_v2a_huber_sv6_ddp_smoke_r4"
PIPELINE_LOG="$REPO/logs/prism_modrescue_v2a_huber_sv6_20260804.log"
CHECKPOINT="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_v3_learnedN_modrescue_nointeraction_s42_20260730/checkpoint_epoch_012.pt"
STATS="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/prism_module_rescue_stats_train64_v2.npz"
EXPECTED_CHECKPOINT_SHA="548b0a0fae8beb5c3f2e312384d0ea781a91e13af8927f1375044186d1f9c3e0"
EXPECTED_STATS_SHA="c1d70a9dc88ee78fad047b2a4ebbf723898ebc3cf5dafb639a35b3cbc1d08147"
PORT_SMOKE="${PRISM_V2A_SV6_SMOKE_PORT:-29950}"
PORT_FULL="${PRISM_V2A_SV6_FULL_PORT:-29951}"

mkdir -p "$REPO/logs"
cd "$REPO"
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
export KMLEE_LAUNCH_SCRIPT="$REPO/scripts/launch_prism_modrescue_v2a_huber_sv6_20260804.sh"
unset KMLEE_PRISM_CCC_CALIBRATION_OUT || true

log() {
  echo "[$(date --iso-8601=seconds)] $*" | tee -a "$PIPELINE_LOG"
}

for required in "$PYTHON" "$CONFIG" "$SMOKE_CONFIG" "$CHECKPOINT" "$STATS"; do
  if [[ ! -f "$required" ]]; then
    log "FATAL missing required file: $required"
    exit 2
  fi
done
if [[ -e "$RUN" || -e "$SMOKE_RUN" ]]; then
  log "FATAL prior production or smoke directory exists; refusing overwrite/resume"
  exit 2
fi

checkpoint_sha="$(sha256sum "$CHECKPOINT" | awk '{print $1}')"
stats_sha="$(sha256sum "$STATS" | awk '{print $1}')"
if [[ "$checkpoint_sha" != "$EXPECTED_CHECKPOINT_SHA" ]]; then
  log "FATAL warm-start SHA mismatch: $checkpoint_sha"
  exit 2
fi
if [[ "$stats_sha" != "$EXPECTED_STATS_SHA" ]]; then
  log "FATAL module-rescue stats SHA mismatch: $stats_sha"
  exit 2
fi

read -r ccc_weight calibration_sha interaction_enabled generator_count <<<"$(
  "$PYTHON" - "$CONFIG" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    config = json.load(handle)
print(
    config["prism_module_rescue"]["ccc_weight"],
    config["prism_module_rescue"]["ccc_calibration_sha256"],
    str(config["precision_medicine"]["interaction_enabled"]).lower(),
    config["decoder"]["n_generators"],
)
PY
)"
if [[ "$ccc_weight" != "0" && "$ccc_weight" != "0.0" ]] || \
  [[ "$calibration_sha" != "None" ]]; then
  log "FATAL v2A must be Huber-only: ccc=$ccc_weight calibration=$calibration_sha"
  exit 2
fi
if [[ "$interaction_enabled" != "false" || "$generator_count" != "414" ]]; then
  log "FATAL scientific arm mismatch: interaction=$interaction_enabled generators=$generator_count"
  exit 2
fi

"$PYTHON" - "$CHECKPOINT" <<'PY' | tee -a "$PIPELINE_LOG"
import sys
import torch

payload = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
state = payload.get(
    "system_state_dict",
    payload.get("system", payload.get("model", payload.get("state_dict", payload))),
)
precision = [str(name) for name in state if "precision_head" in str(name)]
if not precision:
    raise SystemExit("FATAL checkpoint has no precision_head parameters")
print(f"[checkpoint-contract] precision_head_tensors={len(precision)}")
PY

log "sv6 v2A preflight: Huber-only interaction=OFF generators=414"
log "checkpoint_sha=$checkpoint_sha stats_sha=$stats_sha"
sha256sum "$CONFIG" "$SMOKE_CONFIG" | tee -a "$PIPELINE_LOG"
nvidia-smi --query-gpu=index,name,driver_version,memory.used,utilization.gpu \
  --format=csv,noheader | tee -a "$PIPELINE_LOG"
"$PYTHON" -c \
  'import torch; print(f"[cuda-contract] torch={torch.__version__} available={torch.cuda.is_available()} devices={torch.cuda.device_count()} name={torch.cuda.get_device_name(0)}")' \
  | tee -a "$PIPELINE_LOG"

log "sv6 v2A 6-GPU production-topology smoke begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_SMOKE" \
  -m kmlee_bam.training.run_current --config "$SMOKE_CONFIG" \
  2>&1 | tee -a "$PIPELINE_LOG"

"$PYTHON" - "$SMOKE_RUN" <<'PY' | tee -a "$PIPELINE_LOG"
import json
import math
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
with (run / "prism_module_rescue_manifest.json").open(encoding="utf-8") as handle:
    manifest = json.load(handle)
schedule = manifest["resolved_schedule"]
expected = {
    "world_size": 6,
    "local_blocks_per_epoch": 4,
    "local_blocks_per_optimizer_update": 2,
    "global_blocks_per_optimizer_update": 12,
    "optimizer_updates_per_epoch": 2,
}
if schedule != expected:
    raise SystemExit(f"FATAL smoke schedule mismatch: {schedule} != {expected}")
with (run / "history.json").open(encoding="utf-8") as handle:
    history = json.load(handle)
train = history[-1]["train"] if isinstance(history, list) else history["history"][-1]["train"]
required = {
    "metric/prism_module_rescue_celltype_draw_min": 1.0,
    "metric/prism_module_rescue_celltype_draw_max": 1.0,
    "metric/prism_module_rescue_updates_applied": 2.0,
    "metric/prism_module_rescue_updates_skipped": 0.0,
    "metric/prism_module_rescue_global_blocks_per_update": 12.0,
}
for key, expected_value in required.items():
    value = float(train[key])
    if not math.isfinite(value) or value != expected_value:
        raise SystemExit(f"FATAL smoke metric {key}={value}, expected={expected_value}")
print(
    "[smoke-contract] PASS "
    f"schedule={schedule} allowlist_sha={manifest['rescue_parameter_names_sha256']}"
)
PY

log "SMOKE_OK; sv6 v2A 12-epoch training begins"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$PORT_FULL" \
  -m kmlee_bam.training.run_current --config "$CONFIG" \
  2>&1 | tee -a "$PIPELINE_LOG"

log "V2A_TRAINING_COMPLETE; compare validation donor-centered recovery with sv7 v2B"
