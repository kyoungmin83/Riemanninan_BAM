#!/usr/bin/env bash
# Validation-only worker for selected completed SV6 rank-8 checkpoints.
# The official nine-donor test split is never requested.
set -Eeuo pipefail

if [[ "$#" -ne 3 ]]; then
  printf 'usage: %s EPOCH MODULE_GPU EXTRACT_GPU\n' "$0" >&2
  exit 2
fi

EPOCH="$1"
MODULE_GPU="$2"
EXTRACT_GPU="$3"
case "$EPOCH" in
  46|50|55) ;;
  *) printf 'unsupported candidate epoch: %s\n' "$EPOCH" >&2; exit 2 ;;
esac
EPOCH3="$(printf '%03d' "$EPOCH")"

REPO="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826"
ANALYSIS_SCRIPTS="/home/kmlee/project_sv6/kmlee_bam/scripts"
GENERATED_SCRIPTS="$REPO/scripts/generated"
PY="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_freegen_resume_e33_retry2_20260826"
CONFIG="$RUN/resolved_config.json"
CONTRACT_CONFIG="$RUN/source_config.json"
CHECKPOINT="$RUN/checkpoint_epoch_${EPOCH3}.pt"
GENERATOR_SNAPSHOT="$RUN/joint_generator_count_epoch_${EPOCH3}.json"
OUT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_integrated_rank8_primary_e${EPOCH}_candidate_posthoc_20260828"
LOG="$OUT/logs"
STATUS="$OUT/status.log"

cd "$REPO"
export PYTHONPATH="$REPO/src:$ANALYSIS_SCRIPTS"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
export KMLEE_GROUPED_BLOCK_SIZE=8
export KMLEE_PATH_CONSISTENCY=0
export KMLEE_PATH_CONS_MINCELLS=2
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export MPLCONFIGDIR="/tmp/mpl_prism_integrated_rank8_e${EPOCH}_20260828"
export XDG_CACHE_HOME="/tmp/cache_prism_integrated_rank8_e${EPOCH}_20260828"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

mkdir -p "$OUT" "$LOG"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }
for required in "$PY" "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" \
  "$RUN/pathology_correction_rank.jsonl" \
  "$GENERATOR_SNAPSHOT" \
  "$ANALYSIS_SCRIPTS/analyze_prism_checkpoint_history.py" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
  "$GENERATED_SCRIPTS/analyze_prism_strict_reference_leakage.py"; do
  [[ -r "$required" ]] || { status "FATAL missing=$required"; exit 2; }
done

"$PY" - "$CONTRACT_CONFIG" "$CONFIG" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    contract = json.load(handle)
with open(sys.argv[2], encoding="utf-8") as handle:
    resolved = json.load(handle)
for config in (contract, resolved):
    assert config.get("train", {}).get("skip_final_test") is True
assert contract.get("_launch_guard", {}).get(
    "official_test_dataset_open_forbidden"
) is True
assert contract.get("experiment_manifest", {}).get(
    "official_test_donors"
) == "sealed_not_opened"
assert resolved["learned_generator_count"]["protect_unique_gene_coverage"] is False
assert resolved["learned_generator_count"]["protect_singletons"] is True
assert resolved["precision_medicine"]["personal_rank"] == 2
assert resolved["decoder"]["pathology_rank"] == 8
PY

run_gpu() {
  local gpu="$1" name="$2" expected="$3"; shift 3
  if [[ -s "$expected" ]]; then
    status "SKIP $name existing=$expected"
    return 0
  fi
  status "START $name gpu=$gpu"
  if CUDA_VISIBLE_DEVICES="$gpu" nice -n 15 "$@" >"$LOG/$name.log" 2>&1; then
    status "DONE $name gpu=$gpu"
  else
    local rc=$?
    status "FAIL $name gpu=$gpu rc=$rc log=$LOG/$name.log"
    return "$rc"
  fi
}

run_cpu() {
  local name="$1" expected="$2"; shift 2
  if [[ -s "$expected" ]]; then
    status "SKIP $name existing=$expected"
    return 0
  fi
  status "START $name cpu"
  if CUDA_VISIBLE_DEVICES='' nice -n 15 "$@" >"$LOG/$name.log" 2>&1; then
    status "DONE $name cpu"
  else
    local rc=$?
    status "FAIL $name cpu rc=$rc log=$LOG/$name.log"
    return "$rc"
  fi
}

status "BEGIN checkpoint=epoch${EPOCH3} validation_only=true locked_test=SEALED"
sha256sum "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" \
  "$RUN/pathology_correction_rank.jsonl" \
  "$GENERATOR_SNAPSHOT" >"$OUT/provenance.sha256"
cp "$RUN/pathology_correction_rank.jsonl" "$OUT/pathology_correction_rank.jsonl"
cp "$GENERATOR_SNAPSHOT" "$OUT/joint_generator_count_epoch_${EPOCH3}.json"

run_cpu history_audit "$OUT/history/interim_summary.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_checkpoint_history.py" \
    --checkpoint "$CHECKPOINT" --outdir "$OUT/history"

run_gpu "$MODULE_GPU" module_recovery "$OUT/module_disease_val_bal40.json" \
  "$PY" "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --eval-split val \
    --per-donor-celltype 40 --sample-seed 20260814 --ceiling-splits 30 \
    --per-donor --adjust-covariates --max-batches 0 \
    --batch-size 8 --num-workers 0 \
    --dump-fingerprint "$OUT/module_disease_val_bal40_fingerprint.npz" \
    --out-json "$OUT/module_disease_val_bal40.json" & p_module=$!

run_gpu "$EXTRACT_GPU" balanced_extract "$OUT/train_val_bal30_extract.npz" \
  "$PY" "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --include-splits train val \
    --per-donor-celltype 30 --seed 20260814 --batch-size 24 --num-workers 0 \
    --out "$OUT/train_val_bal30_extract.npz" & p_extract=$!

wave1_failed=0
for pid in "$p_module" "$p_extract"; do
  wait "$pid" || wave1_failed=1
done
[[ "$wave1_failed" -eq 0 ]] || { status "FATAL wave1_failed"; exit 3; }

status "BEGIN CPU audits fit=train evaluate=validation"
run_cpu leakage "$OUT/leakage_val.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
    --extract "$OUT/train_val_bal30_extract.npz" --config "$CONFIG" \
    --target-splits val --out-json "$OUT/leakage_val.json" & p_leak=$!

run_cpu leakage_by_celltype "$OUT/leakage_by_celltype_val.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
    --extract "$OUT/train_val_bal30_extract.npz" --target-split val \
    --out-json "$OUT/leakage_by_celltype_val.json" & p_ctleak=$!

run_cpu strict_reference_leakage "$OUT/strict_reference_leakage_val.json" \
  "$PY" "$GENERATED_SCRIPTS/analyze_prism_strict_reference_leakage.py" \
    --extract "$OUT/train_val_bal30_extract.npz" --config "$CONFIG" \
    --target-split val --out-json "$OUT/strict_reference_leakage_val.json" & p_strict=$!

wave2_failed=0
for pid in "$p_leak" "$p_ctleak" "$p_strict"; do
  wait "$pid" || wave2_failed=1
done
[[ "$wave2_failed" -eq 0 ]] || { status "FATAL wave2_failed"; exit 4; }

touch "$OUT/VALIDATION_COMPLETE"
status "COMPLETE checkpoint=epoch${EPOCH3} validation_only=true locked_test=SEALED"
