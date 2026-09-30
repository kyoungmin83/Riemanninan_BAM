#!/usr/bin/env bash
set -u -o pipefail

RUNTIME="/home/kmlee/project_local/kmlee_bam_integrated_pathranklearn_20260824"
ATTEMPT="$RUNTIME/scripts/generated/run_prism_integrated_pathranklearn_attempt_sv7_recovery_20260826.sh"
COMPLETE="$RUNTIME/artifacts/TRAINING_COMPLETE_CLS_RECOVERY_20260826.json"
LOG="$RUNTIME/logs/prism_integrated_pathranklearn_tmux_supervisor_sv7_20260826.log"
LOCK="$RUNTIME/artifacts/prism_integrated_pathranklearn_tmux_supervisor_sv7_20260826.lock"
RETRY_SECONDS="${PRISM_RETRY_SECONDS:-30}"

mkdir -p "$RUNTIME/logs" "$RUNTIME/artifacts"
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "[$(date --iso-8601=seconds)] FATAL another supervisor owns $LOCK" | tee -a "$LOG" >&2
  exit 73
fi

attempt_number=0
while [[ ! -f "$COMPLETE" ]]; do
  attempt_number=$((attempt_number + 1))
  echo "[$(date --iso-8601=seconds)] supervisor attempt=$attempt_number start" | tee -a "$LOG"

  "$ATTEMPT"
  status=$?

  if [[ -f "$COMPLETE" ]]; then
    echo "[$(date --iso-8601=seconds)] supervisor observed completion" | tee -a "$LOG"
    exit 0
  fi

  echo "[$(date --iso-8601=seconds)] supervisor attempt=$attempt_number exit=$status retry_in=${RETRY_SECONDS}s" | tee -a "$LOG" >&2
  sleep "$RETRY_SECONDS"
done

echo "[$(date --iso-8601=seconds)] supervisor found pre-existing completion marker" | tee -a "$LOG"
