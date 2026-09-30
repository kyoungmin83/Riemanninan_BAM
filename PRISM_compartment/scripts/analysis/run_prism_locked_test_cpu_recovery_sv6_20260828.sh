#!/usr/bin/env bash
# Resume the five SV6 test audits after their completed GPU artifacts.
set -Eeuo pipefail

WORKER="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826/scripts/generated/run_prism_locked_test_candidate_sv6_20260828.sh"
OUT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_locked_test_audit_20260828"
STATUS="$OUT/recovery_status.log"
FREE_RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_freegen_resume_e33_retry2_20260826"
PROTECTED_RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826"

status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }
[[ -r "$WORKER" ]] || { status "FATAL missing=$WORKER"; exit 2; }

status "BEGIN recovered_candidates=sv6_e33_protected,sv6_e42,sv6_e47,sv6_e50,sv6_e55 gpu_artifacts=reused official_test=CONSUMED"

bash "$WORKER" sv6_e33_protected 1 "$PROTECTED_RUN/resolved_config.json" "$PROTECTED_RUN/checkpoint_epoch_033.pt" >"$OUT/sv6_e33_protected_recovery.log" 2>&1 & p1=$!
bash "$WORKER" sv6_e42 2 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_042.pt" >"$OUT/sv6_e42_recovery.log" 2>&1 & p2=$!
bash "$WORKER" sv6_e47 3 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_047.pt" >"$OUT/sv6_e47_recovery.log" 2>&1 & p3=$!
bash "$WORKER" sv6_e50 4 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_050.pt" >"$OUT/sv6_e50_recovery.log" 2>&1 & p4=$!
bash "$WORKER" sv6_e55 5 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_055.pt" >"$OUT/sv6_e55_recovery.log" 2>&1 & p5=$!

failed=0
for pid in "$p1" "$p2" "$p3" "$p4" "$p5"; do
  wait "$pid" || failed=1
done
[[ "$failed" -eq 0 ]] || { status "FATAL recovered_candidate_failed"; exit 3; }

touch "$OUT/SV6_RECOVERY_COMPLETE"
status "COMPLETE recovered_candidates=5 gpu_artifacts=reused official_test=CONSUMED"

