#!/usr/bin/env bash
# Post-selection test audit for one pre-locked checkpoint. This consumes test data.
set -Eeuo pipefail

if [[ "$#" -ne 4 ]]; then
  printf 'usage: %s CANDIDATE_ID GPU CONFIG CHECKPOINT\n' "$0" >&2
  exit 2
fi

CANDIDATE_ID="$1"
GPU="$2"
CONFIG="$3"
CHECKPOINT="$4"
case "$CANDIDATE_ID" in
  rank2_e20|sv6_e33_protected|sv6_e42|sv6_e47|sv6_e50|sv6_e55) ;;
  *) printf 'candidate is not in the pre-test lock: %s\n' "$CANDIDATE_ID" >&2; exit 2 ;;
esac

REPO="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826"
ANALYSIS_SCRIPTS="/home/kmlee/project_sv6/kmlee_bam/scripts"
GENERATED_SCRIPTS="$REPO/scripts/generated"
PY="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
LOCK="$GENERATED_SCRIPTS/PRE_TEST_SELECTION_LOCK.json"
LOCK_SHA256="f88f29426e625ceeb2fa3467d000db2e5d6a44a478f71afc21285019354d96ff"
OUT_ROOT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_locked_test_audit_20260828"
OUT="$OUT_ROOT/$CANDIDATE_ID"
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
export MPLCONFIGDIR="/tmp/mpl_prism_locked_test_${CANDIDATE_ID}_20260828"
export XDG_CACHE_HOME="/tmp/cache_prism_locked_test_${CANDIDATE_ID}_20260828"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

mkdir -p "$OUT" "$LOG"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }

for required in \
  "$PY" "$CONFIG" "$CHECKPOINT" "$LOCK" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
  "$GENERATED_SCRIPTS/run_legacy_rank2_extract_locked_test.py" \
  "$GENERATED_SCRIPTS/analyze_prism_strict_reference_locked_test.py"; do
  [[ -r "$required" ]] || { status "FATAL missing=$required"; exit 2; }
done

"$PY" - "$LOCK" "$LOCK_SHA256" "$CANDIDATE_ID" "$CHECKPOINT" <<'PY'
import hashlib
import json
import sys

lock_path, expected_lock_hash, candidate_id, checkpoint_path = sys.argv[1:]
with open(lock_path, "rb") as handle:
    lock_payload = handle.read()
actual_lock_hash = hashlib.sha256(lock_payload).hexdigest()
assert actual_lock_hash == expected_lock_hash, (actual_lock_hash, expected_lock_hash)
lock = json.loads(lock_payload.decode("utf-8"))
assert lock["official_test_previously_used"] is False
assert lock["test_audit_may_change_selection"] is False
record = next(item for item in lock["locked_roles"] if item["id"] == candidate_id)
digest = hashlib.sha256()
with open(checkpoint_path, "rb") as handle:
    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        digest.update(chunk)
assert digest.hexdigest() == record["checkpoint_sha256"], (
    candidate_id,
    digest.hexdigest(),
    record["checkpoint_sha256"],
)
PY

run_gpu() {
  local name="$1" expected="$2"; shift 2
  if [[ -s "$expected" ]]; then
    status "SKIP $name existing=$expected"
    return 0
  fi
  status "START $name gpu=$GPU official_test=AUTHORIZED_CONSUMED"
  if CUDA_VISIBLE_DEVICES="$GPU" nice -n 15 "$@" >"$LOG/$name.log" 2>&1; then
    status "DONE $name gpu=$GPU"
  else
    local rc=$?
    status "FAIL $name gpu=$GPU rc=$rc log=$LOG/$name.log"
    return "$rc"
  fi
}

run_cpu() {
  local name="$1" expected="$2"; shift 2
  if [[ -s "$expected" ]]; then
    status "SKIP $name existing=$expected"
    return 0
  fi
  status "START $name cpu official_test=AUTHORIZED_CONSUMED"
  if CUDA_VISIBLE_DEVICES='' nice -n 15 "$@" >"$LOG/$name.log" 2>&1; then
    status "DONE $name cpu"
  else
    local rc=$?
    status "FAIL $name cpu rc=$rc log=$LOG/$name.log"
    return "$rc"
  fi
}

status "BEGIN candidate=$CANDIDATE_ID gpu=$GPU selection_locked=true official_test=AUTHORIZED_CONSUMED"
cp "$LOCK" "$OUT/PRE_TEST_SELECTION_LOCK.json"
sha256sum "$CONFIG" "$CHECKPOINT" "$LOCK" \
  "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
  "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
  "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
  "$GENERATED_SCRIPTS/run_legacy_rank2_extract_locked_test.py" \
  "$GENERATED_SCRIPTS/analyze_prism_strict_reference_locked_test.py" \
  >"$OUT/provenance.sha256"

run_gpu module_recovery "$OUT/module_disease_test_bal40.json" \
  "$PY" "$ANALYSIS_SCRIPTS/module_disease_eval.py" \
    --config "$CONFIG" --checkpoint "$CHECKPOINT" --eval-split test \
    --per-donor-celltype 40 --sample-seed 20260814 --ceiling-splits 30 \
    --per-donor --adjust-covariates --max-batches 0 \
    --batch-size 8 --num-workers 0 \
    --dump-fingerprint "$OUT/module_disease_test_bal40_fingerprint.npz" \
    --out-json "$OUT/module_disease_test_bal40.json"

if [[ "$CANDIDATE_ID" == "rank2_e20" ]]; then
  run_gpu balanced_extract "$OUT/train_test_bal30_extract.npz" \
    "$PY" "$GENERATED_SCRIPTS/run_legacy_rank2_extract_locked_test.py" \
      "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
      --config "$CONFIG" --checkpoint "$CHECKPOINT" --include-splits train test \
      --per-donor-celltype 30 --seed 20260814 --batch-size 24 --num-workers 0 \
      --out "$OUT/train_test_bal30_extract.npz"
else
  run_gpu balanced_extract "$OUT/train_test_bal30_extract.npz" \
    "$PY" "$ANALYSIS_SCRIPTS/extract_prism_e2e_posthoc.py" \
      --config "$CONFIG" --checkpoint "$CHECKPOINT" --include-splits train test \
      --per-donor-celltype 30 --seed 20260814 --batch-size 24 --num-workers 0 \
      --out "$OUT/train_test_bal30_extract.npz"
fi

status "BEGIN CPU audits fit=train evaluate=locked_test selection_locked=true"
run_cpu leakage "$OUT/leakage_test.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_leakage.py" \
    --extract "$OUT/train_test_bal30_extract.npz" --config "$CONFIG" \
    --target-splits test --out-json "$OUT/leakage_test.json" & p_leak=$!

run_cpu leakage_by_celltype "$OUT/leakage_by_celltype_test.json" \
  "$PY" "$ANALYSIS_SCRIPTS/analyze_prism_leakage_by_celltype.py" \
    --extract "$OUT/train_test_bal30_extract.npz" --target-split test \
    --out-json "$OUT/leakage_by_celltype_test.json" & p_ctleak=$!

run_cpu strict_reference_leakage "$OUT/strict_reference_leakage_test.json" \
  "$PY" "$GENERATED_SCRIPTS/analyze_prism_strict_reference_locked_test.py" \
    --extract "$OUT/train_test_bal30_extract.npz" --config "$CONFIG" \
    --pretest-lock "$LOCK" --expected-lock-sha256 "$LOCK_SHA256" \
    --out-json "$OUT/strict_reference_leakage_test.json" & p_strict=$!

failed=0
for pid in "$p_leak" "$p_ctleak" "$p_strict"; do
  wait "$pid" || failed=1
done
[[ "$failed" -eq 0 ]] || { status "FATAL cpu_audit_failed"; exit 4; }

"$PY" - "$OUT" "$LOCK_SHA256" <<'PY'
import json
import os
import sys

out, lock_sha256 = sys.argv[1:]
with open(os.path.join(out, "module_disease_test_bal40.json"), encoding="utf-8") as handle:
    module = json.load(handle)
with open(os.path.join(out, "leakage_test.json"), encoding="utf-8") as handle:
    leakage = json.load(handle)
with open(os.path.join(out, "leakage_by_celltype_test.json"), encoding="utf-8") as handle:
    by_celltype = json.load(handle)
with open(os.path.join(out, "strict_reference_leakage_test.json"), encoding="utf-8") as handle:
    strict = json.load(handle)
assert module["evaluation_provenance"]["split"] == "test"
assert len(module["per_donor"]) == 9
assert leakage["test_used"] is True
assert leakage["cell_counts"]["test"] > 0
assert leakage["cell_counts"].get("val", 0) == 0
assert by_celltype["test_used"] is True
assert strict["test_used"] is True
assert strict["selection_locked_before_test"] is True
assert strict["target_split"] == "test"
assert strict["pretest_lock_sha256"] == lock_sha256
PY

touch "$OUT/TEST_AUDIT_COMPLETE"
status "COMPLETE candidate=$CANDIDATE_ID selection_locked=true official_test=CONSUMED"
