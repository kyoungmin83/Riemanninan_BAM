#!/usr/bin/env bash
# Run all six pre-locked checkpoint test audits concurrently, one GPU each.
set -Eeuo pipefail

WORKER="/home/kmlee/project_sv6/kmlee_bam_integrated_rank8_primary_resume_e28_20260826/scripts/generated/run_prism_locked_test_candidate_sv6_20260828.sh"
OUT="/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_locked_test_audit_20260828"
STATUS="$OUT/status.log"
FREE_RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_freegen_resume_e33_retry2_20260826"
PROTECTED_RUN="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_rank8_primary_s42_sv6_resume_e28_retry2_20260826"
RANK2_CONFIG="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_rank2_support_query_infonce_fixed97_s42_20260811/resolved_config.json"
RANK2_CHECKPOINT="/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/prism_rank2_agp_infonce_warmstart_20260811/rank2_e20_fixed97.pt"

mkdir -p "$OUT"
status() { printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "$STATUS"; }
[[ -r "$WORKER" ]] || { status "FATAL missing=$WORKER"; exit 2; }

status "BEGIN candidates=rank2_e20,sv6_e33_protected,sv6_e42,sv6_e47,sv6_e50,sv6_e55 selection_locked=true official_test=AUTHORIZED_CONSUMED"

bash "$WORKER" rank2_e20 0 "$RANK2_CONFIG" "$RANK2_CHECKPOINT" >"$OUT/rank2_e20_driver.log" 2>&1 & p0=$!
bash "$WORKER" sv6_e33_protected 1 "$PROTECTED_RUN/resolved_config.json" "$PROTECTED_RUN/checkpoint_epoch_033.pt" >"$OUT/sv6_e33_protected_driver.log" 2>&1 & p1=$!
bash "$WORKER" sv6_e42 2 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_042.pt" >"$OUT/sv6_e42_driver.log" 2>&1 & p2=$!
bash "$WORKER" sv6_e47 3 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_047.pt" >"$OUT/sv6_e47_driver.log" 2>&1 & p3=$!
bash "$WORKER" sv6_e50 4 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_050.pt" >"$OUT/sv6_e50_driver.log" 2>&1 & p4=$!
bash "$WORKER" sv6_e55 5 "$FREE_RUN/resolved_config.json" "$FREE_RUN/checkpoint_epoch_055.pt" >"$OUT/sv6_e55_driver.log" 2>&1 & p5=$!

failed=0
for pid in "$p0" "$p1" "$p2" "$p3" "$p4" "$p5"; do
  wait "$pid" || failed=1
done
[[ "$failed" -eq 0 ]] || { status "FATAL candidate_worker_failed"; exit 3; }

touch "$OUT/TEST_AUDIT_COMPLETE"
status "COMPLETE candidates=6 selection_locked=true official_test=CONSUMED"

