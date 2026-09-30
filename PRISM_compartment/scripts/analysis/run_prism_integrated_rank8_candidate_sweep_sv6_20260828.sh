#!/usr/bin/env bash
# Run the three validation-only SV6 candidate audits concurrently in one tmux job.
set -Eeuo pipefail

WORKER="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826/scripts/generated/run_prism_integrated_rank8_candidate_posthoc_sv6_20260828.sh"
OUT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_integrated_rank8_primary_candidate_sweep_20260828"
STATUS="$OUT/status.log"
mkdir -p "$OUT"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }

[[ -r "$WORKER" ]] || { status "FATAL missing=$WORKER"; exit 2; }
status "BEGIN candidates=46,50,55 validation_only=true locked_test=SEALED"

bash "$WORKER" 46 0 1 >"$OUT/e46_driver.log" 2>&1 & p46=$!
bash "$WORKER" 50 2 3 >"$OUT/e50_driver.log" 2>&1 & p50=$!
bash "$WORKER" 55 4 5 >"$OUT/e55_driver.log" 2>&1 & p55=$!

failed=0
for pid in "$p46" "$p50" "$p55"; do
  wait "$pid" || failed=1
done
[[ "$failed" -eq 0 ]] || { status "FATAL candidate_worker_failed"; exit 3; }

touch "$OUT/VALIDATION_COMPLETE"
status "COMPLETE candidates=46,50,55 validation_only=true locked_test=SEALED"
