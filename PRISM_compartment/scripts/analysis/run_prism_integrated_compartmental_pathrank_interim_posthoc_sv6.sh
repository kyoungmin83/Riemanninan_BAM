#!/usr/bin/env bash
# Validation-only interim audit for the active SV6 compartmental-threshold run.
set -Eeuo pipefail

REPO="/home/kmlee/project_sv6/kmlee_bam_integrated_compartmental_pathrank_20260828"
ANALYSIS_SCRIPTS="/home/kmlee/project_sv6/kmlee_bam/scripts"
PY="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_compartmental_pathrank_s42_sv6_freegen_recovery_e12_20260829"
EPOCH="${PRISM_POSTHOC_EPOCH:-037}"
EPOCH_TAG="${EPOCH#0}"
POSTHOC_DATE="${PRISM_POSTHOC_DATE:-20260831}"
CONFIG="$RUN/resolved_config.json"
CONTRACT_CONFIG="$RUN/source_config.json"
CHECKPOINT="$RUN/checkpoint_epoch_${EPOCH}.pt"
GENERATOR_SNAPSHOT="$RUN/joint_generator_count_epoch_${EPOCH}.json"
OUT="${PRISM_POSTHOC_OUT:-/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_integrated_compartmental_pathrank_sv6_e${EPOCH_TAG}_interim_posthoc_${POSTHOC_DATE}}"
LOG="$OUT/logs"
STATUS="$OUT/status.log"
STRICT_SCRIPT="$REPO/scripts/analysis/analyze_prism_strict_reference_leakage.py"

cd "$REPO"
export PYTHONPATH="$REPO/src:$ANALYSIS_SCRIPTS"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
export KMLEE_GROUPED_BLOCK_SIZE=8
export KMLEE_PATH_CONSISTENCY=0
export KMLEE_PATH_CONS_MINCELLS=2
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export MPLCONFIGDIR="/tmp/mpl_prism_comp_e${EPOCH_TAG}_${POSTHOC_DATE}"
export XDG_CACHE_HOME="/tmp/cache_prism_comp_e${EPOCH_TAG}_${POSTHOC_DATE}"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

mkdir -p "$OUT" "$LOG" "$OUT/history"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }

for required in "$PY" "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" \
  "$RUN/pathology_correction_rank.jsonl" "$GENERATOR_SNAPSHOT" \
  "$ANALYSIS_SCRIPTS/analyze_prism_checkpoint_history.py" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
  "$STRICT_SCRIPT"; do
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
assert contract.get("_launch_guard", {}).get("official_test_dataset_open_forbidden") is True
assert contract.get("experiment_manifest", {}).get("official_test_donors") == "sealed_not_opened"
assert resolved["learned_generator_count"]["minimum_active_generators"] is None
assert resolved["learned_generator_count"]["protect_unique_gene_coverage"] is False
assert resolved["learned_generator_count"]["protect_singletons"] is True
assert resolved["precision_medicine"]["personal_rank"] == 2
assert resolved["precision_medicine"]["module_local_nonlinear_variant"] == "compartmental_threshold"
assert resolved["decoder"]["pathology_rank"] == 0
assert resolved["learned_pathology_rank"]["enabled"] is True
assert resolved["learned_pathology_rank"]["freeze_epoch"] == 37
assert resolved["integrated_phase_curriculum"]["phase2_pathology_scale"] > 0.0
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

status "BEGIN checkpoint=epoch${EPOCH} validation_only=true paper_arm=compartmental_threshold locked_test=NOT_USED"
cp "$RUN/pathology_correction_rank.jsonl" "$OUT/pathology_correction_rank_snapshot.jsonl"
cp "$GENERATOR_SNAPSHOT" "$OUT/joint_generator_count_epoch_${EPOCH}.json"
sha256sum "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" \
  "$OUT/pathology_correction_rank_snapshot.jsonl" \
  "$OUT/joint_generator_count_epoch_${EPOCH}.json" \
  "$REPO/src/kmlee_bam/model/precision_medicine.py" \
  "$REPO/src/kmlee_bam/training/adaptive_subgroup_trainer.py" \
  "$REPO/src/kmlee_bam/training/runner_base.py" >"$OUT/provenance.sha256"

run_cpu history_audit "$OUT/history/interim_summary.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_checkpoint_history.py" \
    --checkpoint "$CHECKPOINT" --outdir "$OUT/history"

# The training run occupies all GPUs but leaves >30 GiB free per card.  These
# two low-priority validation jobs use separate cards and conservative batches.
run_gpu 0 module_recovery "$OUT/module_disease_val_bal40.json" \
  "$PY" "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --eval-split val \
    --per-donor-celltype 40 --sample-seed 20260814 --ceiling-splits 30 \
    --per-donor --adjust-covariates --max-batches 0 \
    --batch-size 8 --num-workers 0 \
    --dump-fingerprint "$OUT/module_disease_val_bal40_fingerprint.npz" \
    --out-json "$OUT/module_disease_val_bal40.json" & p_module=$!

run_gpu 1 balanced_extract "$OUT/train_val_bal30_extract.npz" \
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
  "$PY" "$STRICT_SCRIPT" \
    --extract "$OUT/train_val_bal30_extract.npz" --config "$CONFIG" \
    --target-split val --out-json "$OUT/strict_reference_leakage_val.json" & p_strict=$!

wave2_failed=0
for pid in "$p_leak" "$p_ctleak" "$p_strict"; do
  wait "$pid" || wave2_failed=1
done
[[ "$wave2_failed" -eq 0 ]] || { status "FATAL wave2_failed"; exit 4; }

touch "$OUT/VALIDATION_COMPLETE"
status "COMPLETE checkpoint=epoch${EPOCH} validation_only=true paper_arm=compartmental_threshold locked_test=NOT_USED"
