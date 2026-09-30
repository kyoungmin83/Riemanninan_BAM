#!/usr/bin/env bash
# Descriptive post-hoc test reuse for one validation-locked SV7 checkpoint.
# The official cohort was already consumed by the prior SV6 audit.
set -Eeuo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "usage: $0 EPOCH MODULE_GPU EXTRACT_GPU" >&2
  exit 64
fi

EPOCH_RAW="$1"
MODULE_GPU="$2"
EXTRACT_GPU="$3"
case "$EPOCH_RAW" in
  33|36|39) ;;
  *) echo "epoch is outside the validation-locked SV7 test scope: $EPOCH_RAW" >&2; exit 64 ;;
esac
[[ "$MODULE_GPU" =~ ^[0-5]$ ]] || { echo "MODULE_GPU must be 0..5" >&2; exit 64; }
[[ "$EXTRACT_GPU" =~ ^[0-5]$ ]] || { echo "EXTRACT_GPU must be 0..5" >&2; exit 64; }
[[ "$MODULE_GPU" != "$EXTRACT_GPU" ]] || { echo "GPU assignments must differ" >&2; exit 64; }
printf -v EPOCH3 '%03d' "$EPOCH_RAW"
CANDIDATE_ID="sv7_pathrank_e${EPOCH_RAW}"

REPO="/home/kmlee/project_sv6/kmlee_bam_integrated_pathranklearn_20260824"
ANALYSIS_SCRIPTS="/home/kmlee/project_sv6/kmlee_bam/scripts"
GENERATED_SCRIPTS="$REPO/scripts/generated"
PY="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_pathranklearn_s42_sv7_freegen_rankactive_resume_e19_20260826"
CONFIG="$RUN/resolved_config.json"
CONTRACT_CONFIG="$RUN/source_config.json"
CHECKPOINT="$RUN/checkpoint_epoch_${EPOCH3}.pt"
SCOPE_LOCK="$GENERATED_SCRIPTS/SV7_CONSUMED_TEST_SCOPE_LOCK.json"
SCOPE_LOCK_SHA256="5105a96025c69f90d7b19128e3a1e52ffc80608d5d546fd574c81f6a87989790"
OUT_ROOT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_sv7_pathrank_consumed_test_posthoc_20260828"
OUT="$OUT_ROOT/$CANDIDATE_ID"
LOG="$OUT/logs"
STATUS="$OUT/status.log"

cd "$REPO"
export PYTHONPATH="$REPO/src:$ANALYSIS_SCRIPTS"
export PRISM_ANALYSIS_SCRIPTS="$ANALYSIS_SCRIPTS"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
export KMLEE_GROUPED_BLOCK_SIZE=8
export KMLEE_PATH_CONSISTENCY=0
export KMLEE_PATH_CONS_MINCELLS=2
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export MPLCONFIGDIR="/tmp/mpl_prism_sv7_consumed_test_e${EPOCH_RAW}_20260828"
export XDG_CACHE_HOME="/tmp/cache_prism_sv7_consumed_test_e${EPOCH_RAW}_20260828"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

mkdir -p "$OUT" "$LOG"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }

for required in \
  "$PY" "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" "$SCOPE_LOCK" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$GENERATED_SCRIPTS/analyze_prism_leakage_by_celltype_consumed_test.py" \
  "$GENERATED_SCRIPTS/analyze_prism_strict_reference_consumed_test.py"; do
  [[ -r "$required" ]] || { status "FATAL missing=$required"; exit 2; }
done

"$PY" - "$SCOPE_LOCK" "$SCOPE_LOCK_SHA256" "$CANDIDATE_ID" "$EPOCH_RAW" \
  "$CHECKPOINT" "$CONFIG" "$CONTRACT_CONFIG" <<'PY'
import hashlib
import json
import sys

scope_path, expected_scope_hash, candidate_id, epoch, checkpoint_path, config_path, contract_path = sys.argv[1:]

def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()

with open(scope_path, "rb") as handle:
    scope_payload = handle.read()
actual_scope_hash = hashlib.sha256(scope_payload).hexdigest()
assert actual_scope_hash == expected_scope_hash, (actual_scope_hash, expected_scope_hash)
scope = json.loads(scope_payload.decode("utf-8"))
assert scope["official_test_status_at_scope_lock"] == "consumed_by_prior_sv6_locked_audit"
assert scope["official_test_independence"] == "not_an_untouched_holdout_for_sv7_posthoc"
assert scope["test_reuse_may_change_validation_roles"] is False
assert scope["test_reuse_may_authorize_baseline_replacement"] is False
record = next(item for item in scope["locked_roles"] if item["id"] == candidate_id)
assert int(record["epoch"]) == int(epoch)
assert digest(checkpoint_path) == record["checkpoint_sha256"]
assert digest(config_path) == scope["runtime"]["resolved_config_sha256"]
with open(config_path, encoding="utf-8") as handle:
    config = json.load(handle)
with open(contract_path, encoding="utf-8") as handle:
    contract = json.load(handle)
for item in (config, contract):
    assert item.get("train", {}).get("skip_final_test") is True
assert config["precision_medicine"]["personal_rank"] == 2
assert config["decoder"]["pathology_rank"] == 0
assert config["learned_pathology_rank"]["enabled"] is True
PY

run_gpu() {
  local gpu="$1" name="$2" expected="$3"; shift 3
  if [[ -s "$expected" ]]; then
    status "SKIP $name existing=$expected"
    return 0
  fi
  status "START $name gpu=$gpu official_test=CONSUMED_REUSE"
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
  status "START $name cpu official_test=CONSUMED_REUSE"
  if CUDA_VISIBLE_DEVICES='' nice -n 15 "$@" >"$LOG/$name.log" 2>&1; then
    status "DONE $name cpu"
  else
    local rc=$?
    status "FAIL $name cpu rc=$rc log=$LOG/$name.log"
    return "$rc"
  fi
}

status "BEGIN candidate=$CANDIDATE_ID module_gpu=$MODULE_GPU extract_gpu=$EXTRACT_GPU official_test=CONSUMED_REUSE independent_holdout=false"
cp "$SCOPE_LOCK" "$OUT/SV7_CONSUMED_TEST_SCOPE_LOCK.json"
sha256sum "$CONFIG" "$CONTRACT_CONFIG" "$CHECKPOINT" "$SCOPE_LOCK" \
  "$REPO/src/kmlee_bam/model/precision_medicine.py" \
  "$REPO/src/kmlee_bam/training/adaptive_subgroup_trainer.py" \
  "$REPO/src/kmlee_bam/training/runner_base.py" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$GENERATED_SCRIPTS/analyze_prism_leakage_by_celltype_consumed_test.py" \
  "$GENERATED_SCRIPTS/analyze_prism_strict_reference_consumed_test.py" \
  >"$OUT/provenance.sha256"

run_gpu "$MODULE_GPU" module_recovery "$OUT/module_disease_test_bal40.json" \
  "$PY" "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --eval-split test \
    --per-donor-celltype 40 --sample-seed 20260814 --ceiling-splits 30 \
    --per-donor --adjust-covariates --max-batches 0 \
    --batch-size 8 --num-workers 0 \
    --dump-fingerprint "$OUT/module_disease_test_bal40_fingerprint.npz" \
    --out-json "$OUT/module_disease_test_bal40.json" & p_module=$!

run_gpu "$EXTRACT_GPU" balanced_extract "$OUT/train_test_bal30_extract.npz" \
  "$PY" "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --include-splits train test \
    --per-donor-celltype 30 --seed 20260814 --batch-size 24 --num-workers 0 \
    --out "$OUT/train_test_bal30_extract.npz" & p_extract=$!

wave1_failed=0
for pid in "$p_module" "$p_extract"; do
  wait "$pid" || wave1_failed=1
done
[[ "$wave1_failed" -eq 0 ]] || { status "FATAL wave1_failed"; exit 3; }

status "BEGIN CPU audits fit=train evaluate=consumed_test"
run_cpu leakage "$OUT/leakage_test.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
    --extract "$OUT/train_test_bal30_extract.npz" --config "$CONFIG" \
    --target-splits test --out-json "$OUT/leakage_test.json" & p_leak=$!

run_cpu leakage_by_celltype "$OUT/leakage_by_celltype_test.json" \
  "$PY" "$GENERATED_SCRIPTS/analyze_prism_leakage_by_celltype_consumed_test.py" \
    --extract "$OUT/train_test_bal30_extract.npz" \
    --scope-lock "$SCOPE_LOCK" --expected-lock-sha256 "$SCOPE_LOCK_SHA256" \
    --out-json "$OUT/leakage_by_celltype_test.json" & p_ctleak=$!

run_cpu strict_reference_leakage "$OUT/strict_reference_leakage_test.json" \
  "$PY" "$GENERATED_SCRIPTS/analyze_prism_strict_reference_consumed_test.py" \
    --extract "$OUT/train_test_bal30_extract.npz" --config "$CONFIG" \
    --scope-lock "$SCOPE_LOCK" --expected-lock-sha256 "$SCOPE_LOCK_SHA256" \
    --out-json "$OUT/strict_reference_leakage_test.json" & p_strict=$!

wave2_failed=0
for pid in "$p_leak" "$p_ctleak" "$p_strict"; do
  wait "$pid" || wave2_failed=1
done
[[ "$wave2_failed" -eq 0 ]] || { status "FATAL wave2_failed"; exit 4; }

"$PY" - "$OUT" "$SCOPE_LOCK_SHA256" "$CANDIDATE_ID" <<'PY'
import json
import os
import sys

out, scope_sha256, candidate_id = sys.argv[1:]
def load(name):
    with open(os.path.join(out, name), encoding="utf-8") as handle:
        return json.load(handle)

module = load("module_disease_test_bal40.json")
leakage = load("leakage_test.json")
by_celltype = load("leakage_by_celltype_test.json")
strict = load("strict_reference_leakage_test.json")
assert module["evaluation_provenance"]["split"] == "test"
assert module["evaluation_provenance"]["selected_cells"] == 8478
assert len(module["per_donor"]) == 9
assert len(module["per_celltype"]) == 24
assert leakage["test_used"] is True
assert leakage["cell_counts"]["test"] == 6411
assert leakage["cell_counts"].get("val", 0) == 0
assert by_celltype["test_used"] is True
assert by_celltype["official_test_status"] == "consumed_reuse"
assert by_celltype["scope_lock_sha256"] == scope_sha256
assert len(by_celltype["per_celltype"]) == 24
assert strict["test_used"] is True
assert strict["official_test_status"] == "consumed_reuse"
assert strict["scope_lock_sha256"] == scope_sha256
receipt = {
    "schema_version": "kmlee_bam.sv7_consumed_test_candidate_completion.v1",
    "candidate_id": candidate_id,
    "official_test_status": "consumed_reuse",
    "independent_holdout": False,
    "test_results_may_change_validation_roles": False,
    "module_test_cells": module["evaluation_provenance"]["selected_cells"],
    "module_test_donors": len(module["per_donor"]),
    "leak_train_cells": leakage["cell_counts"]["train"],
    "leak_test_cells": leakage["cell_counts"]["test"],
    "celltypes": len(by_celltype["per_celltype"]),
    "strict_test_donors": strict["donor_counts"]["test"],
    "scope_lock_sha256": scope_sha256,
}
with open(os.path.join(out, "completion_receipt.json"), "w", encoding="utf-8") as handle:
    json.dump(receipt, handle, indent=2, ensure_ascii=False)
    handle.write("\n")
PY

touch "$OUT/CONSUMED_TEST_POSTHOC_COMPLETE"
status "COMPLETE candidate=$CANDIDATE_ID official_test=CONSUMED_REUSE independent_holdout=false"

exit 0
