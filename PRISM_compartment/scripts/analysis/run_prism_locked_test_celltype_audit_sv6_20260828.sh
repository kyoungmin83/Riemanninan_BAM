#!/usr/bin/env bash
# Run the lock-gated cell-type test leakage audit for all six fixed checkpoints.
set -Eeuo pipefail

ROOT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_locked_test_audit_20260828"
GENERATED="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826/scripts/generated"
PY="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
SCRIPT="$GENERATED/analyze_prism_leakage_by_celltype_locked_test.py"
LOCK="$GENERATED/PRE_TEST_SELECTION_LOCK.json"
LOCK_SHA256="f88f29426e625ceeb2fa3467d000db2e5d6a44a478f71afc21285019354d96ff"
STATUS="$ROOT/celltype_audit_status.log"
CANDIDATES=(rank2_e20 sv6_e33_protected sv6_e42 sv6_e47 sv6_e50 sv6_e55)

status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }
for required in "$PY" "$SCRIPT" "$LOCK"; do
  [[ -r "$required" ]] || { status "FATAL missing=$required"; exit 2; }
done

run_one() {
  local id="$1"
  local out="$ROOT/$id"
  local extract="$out/train_test_bal30_extract.npz"
  local result="$out/leakage_by_celltype_test.json"
  local log="$out/logs/leakage_by_celltype_locked_test.log"
  if [[ -s "$result" ]]; then
    status "SKIP candidate=$id existing=$result"
    return 0
  fi
  local waited=0
  while [[ ! -s "$extract" ]]; do
    if [[ "$waited" -ge 2400 ]]; then
      status "FAIL candidate=$id reason=extract_timeout"
      return 3
    fi
    sleep 10
    waited=$((waited + 10))
  done
  status "START candidate=$id waited_seconds=$waited official_test=CONSUMED"
  if CUDA_VISIBLE_DEVICES='' nice -n 15 "$PY" "$SCRIPT" \
      --extract "$extract" \
      --pretest-lock "$LOCK" \
      --expected-lock-sha256 "$LOCK_SHA256" \
      --out-json "$result" >"$log" 2>&1; then
    status "DONE candidate=$id"
  else
    local rc=$?
    status "FAIL candidate=$id rc=$rc log=$log"
    return "$rc"
  fi
}

status "BEGIN candidates=6 selection_locked=true official_test=CONSUMED"
pids=()
for id in "${CANDIDATES[@]}"; do
  run_one "$id" &
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  wait "$pid" || failed=1
done
[[ "$failed" -eq 0 ]] || { status "FATAL one_or_more_candidates_failed"; exit 4; }

"$PY" - "$ROOT" "$LOCK_SHA256" "${CANDIDATES[@]}" <<'PY'
import json
import os
import sys

root, lock_sha256, *candidate_ids = sys.argv[1:]
assert len(candidate_ids) == 6
for candidate_id in candidate_ids:
    path = os.path.join(root, candidate_id, "leakage_by_celltype_test.json")
    with open(path, encoding="utf-8") as handle:
        result = json.load(handle)
    assert result["test_used"] is True
    assert result["selection_locked_before_test"] is True
    assert result["pretest_lock_sha256"] == lock_sha256
    assert result["target_split"] == "test"
    assert len(result["per_celltype"]) == 24
PY

touch "$ROOT/CELLTYPE_TEST_AUDIT_COMPLETE"
status "COMPLETE candidates=6 selection_locked=true official_test=CONSUMED"

