#!/usr/bin/env bash
set -Eeuo pipefail

# Fail-closed fresh retry of PRISM v5 on sv6.  The failed run is immutable and
# is used only as provenance.  A dedicated six-GPU rank-0-validation cadence
# regression must complete gated epoch 3 before the new full run is authorized.

REPO="/home/kmlee/project_sv6/kmlee_bam"
PYTHON="/home/kmlee/miniconda3/envs/rie_bam/bin/python"
FULL_REVIEW="$REPO/configs/train_config_prism_v5_sv6_joint_generator_count_s42_retry1_20260806.json"
BOUNDARY_SMOKE_CONFIG="$REPO/configs/_boundary_smoke_train_config_prism_v5_sv6_joint_generator_count_s42_retry1_20260806.json"
DESIGN_DOC="$REPO/doc/prism_v5_online_joint_generator_count_design_20260805.md"

FAILED_RUN="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_v5_sv6_joint_l0_generator_count_s42_20260805"
FAILED_CHECKPOINT="$FAILED_RUN/checkpoint_epoch_001.pt"
FAILED_AUTHORIZED_CONFIG="$FAILED_RUN/authorized_train_config.json"
RETRY_RUN="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_v5_sv6_joint_l0_generator_count_s42_retry1_20260806"
BOUNDARY_SMOKE_RUN="/tmp/prism_v5_joint_generator_count_boundary_smoke_sv6_retry1_20260806"
AUTHORIZED_CONFIG="$RETRY_RUN/authorized_train_config.json"
APPROVAL_RECEIPT="$RETRY_RUN/approval_receipt.json"
BOUNDARY_SMOKE_AUDIT="$REPO/analysis_outputs/prism_v5_retry1_preflight_20260806/sv6_boundary_smoke_audit.json"
LOG="$REPO/logs/prism_v5_joint_generator_count_sv6_retry1_20260806.log"

CT64_ORIGINAL="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_dlpfc_mtg_bins4_v31a_cons_ep22_h8_ctemb64/checkpoint_best.pt"
CT64_BACKUP="$REPO/backups/ct64_canonical_20260805/checkpoint_best.pt"
PRISM_ORIGINAL="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_modrescue_v2a_huber_nointeraction_s42_20260804/checkpoint_epoch_011.pt"
PRISM_BACKUP="$REPO/backups/prism_modrescue_v2a_selected_ep11_20260805/checkpoint_epoch_011.pt"
PRISM_SOURCE_ORIGINAL="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_modrescue_v2a_huber_nointeraction_s42_20260804/source_config.json"
PRISM_SOURCE_BACKUP="$REPO/backups/prism_modrescue_v2a_selected_ep11_20260805/source_config.json"
SELECTION_ORIGINAL="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_modrescue_v2a_huber_nointeraction_s42_20260804/posthoc_v2a_20260805/epoch_selection/epoch_selection_summary.json"
SELECTION_BACKUP="$REPO/backups/prism_modrescue_v2a_selected_ep11_20260805/epoch_selection_summary.json"
PRIOR_SMOKE_AUDIT="$REPO/analysis_outputs/prism_v5_preflight_20260805/sv6_smoke_audit.json"

STATS="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/prism_module_rescue_stats_train64_v2.npz"
SOURCE_CONFIG="$REPO/configs/train_config_kmlee_bam_dlpfc_mtg_bins4_prism_modrescue_v2b_ccc_nointeraction_s42_sv6.json"
CONTEXT_NPZ="$REPO/outputs/population_fingerprint_20260714/ctemb64_epoch12_full_region_pseudobulk.npz"
REGISTRY="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/final_gene_module_registry_compact.json"
ACTIVITY_WEIGHT="/home/kmlee/project/kmlee_bam/outputs/SEA_AD_OUTPUTS/module_mapping_outputs/dlpfc_mtg/combined/final_registry/activity_weight_kme_or_l2_membership.npz"
CCC_CALIBRATION="$REPO/outputs/prism_modrescue_v2_ccc_calibration_20260804/calibration_sv7.json"

EXPECTED_FULL_CONFIG_SHA="625c45399235cdf8457164c0fc2c53fb7813b36f561c513c06bf79422fcf9dd6"
EXPECTED_BOUNDARY_SMOKE_CONFIG_SHA="7ca27affedf2be2874bdeb7631ec5698c99ddefae5bb7cdc0fa80257b8158b99"
EXPECTED_DESIGN_SHA="9faf01974dcee82e327e7ba7cc55c0a8265bbce1e1f9dee6a40ce9280109d32e"
EXPECTED_FAILED_CHECKPOINT_SHA="5cc9d5f2f31edc0fc5975b6df749c677fe42af3e0c24983764f50714f39fe4b4"
EXPECTED_FAILED_AUTHORIZED_CONFIG_SHA="445f5ff01427ef9ede6ccb149726f329b6a382a621711ac3c79828e26a102d06"
EXPECTED_PRIOR_SMOKE_AUDIT_SHA="4f53081d4f8ee4bddc99258076698c6b0814821a9c0827d98134fb5ae0ba1c60"
EXPECTED_CT64_SHA="ef3f6241052b4e691132b9b3fd098fb848f82cfb79a43e9473f30bd230fab44a"
EXPECTED_PRISM_SHA="9bee5aba99e38453499d6cd2333831509b811945a33e77c2b381eb3ea655c555"
EXPECTED_PRISM_SOURCE_SHA="d276439763384e46a9870095115042fe338c991e5266051d31396b6307bf5f9c"
EXPECTED_SELECTION_SHA="6fc9a1e6401e03c69d6dd791b916025a24ab270b76d8f18a87e245774155ec4e"
EXPECTED_STATS_SHA="c1d70a9dc88ee78fad047b2a4ebbf723898ebc3cf5dafb639a35b3cbc1d08147"
EXPECTED_SOURCE_CONFIG_SHA="b8e2aa66ef2c009ab43e9c26b3f8dd1effc1bbe1a79d86fa257b04d7c3381708"
EXPECTED_CONTEXT_SHA="ac39aade72eee687fdaa861dbe06ebd21ce2a5d97db43966397befc6677426f8"
EXPECTED_REGISTRY_SHA="0bc8ab11cec35387b29ee6a6696833a8fb75d7fb0a75614c5488657d73b7d973"
EXPECTED_ACTIVITY_WEIGHT_SHA="2043b5f421a79097f326156d0ed3e23245f29db7ad1a35effc0b90daf750d6df"
EXPECTED_CCC_CALIBRATION_SHA="02bd296bfd2df02de0a54887c5288970d2ca824652e0ee9f7c2cfe6b90ea2f58"
EXPECTED_CADENCE_PATCH_SHA="39f21a13027b20726a4cc833b10ff9f25f7a674850c8ce31ba61f3bfedf4510b"
BOUNDARY_SMOKE_PORT="${PRISM_V5_SV6_RETRY_SMOKE_PORT:-29973}"
FULL_PORT="${PRISM_V5_SV6_RETRY_FULL_PORT:-29974}"

mkdir -p "$REPO/logs" "$(dirname "$BOUNDARY_SMOKE_AUDIT")"
cd "$REPO"

export PYTHONPATH="$REPO/src:$REPO/scripts"
export CUDA_VISIBLE_DEVICES="0,1,2,3,4,5"
export OMP_NUM_THREADS="8"
export NCCL_P2P_DISABLE="1"
export NCCL_IB_DISABLE="1"
export PYTHONUNBUFFERED="1"
export PYTHONFAULTHANDLER="1"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export TORCH_NCCL_ASYNC_ERROR_HANDLING="1"
export TORCH_NCCL_DESYNC_DEBUG="1"
export TORCH_NCCL_TRACE_BUFFER_SIZE="2000"
export TORCH_NCCL_DUMP_ON_TIMEOUT="1"
export KMLEE_TRAIN_BOOT_DEBUG="1"
export KMLEE_GROUPED_BLOCK_SIZE="8"
export KMLEE_PATH_CONSISTENCY="0"
export KMLEE_PATH_CONS_MINCELLS="2"
export KMLEE_LAUNCH_SCRIPT="$REPO/scripts/launch_prism_v5_joint_generator_count_sv6_retry1_20260806.sh"

log() {
  echo "[$(date --iso-8601=seconds)] $*" | tee -a "$LOG"
}

verify_sha() {
  local path="$1" expected="$2" label="$3" observed
  observed="$(sha256sum "$path" | awk '{print $1}')"
  if [[ "$observed" != "$expected" ]]; then
    log "FATAL $label SHA256 mismatch: $observed != $expected"
    exit 2
  fi
}

required_files=(
  "$PYTHON" "$FULL_REVIEW" "$BOUNDARY_SMOKE_CONFIG" "$DESIGN_DOC"
  "$FAILED_CHECKPOINT" "$FAILED_AUTHORIZED_CONFIG" "$PRIOR_SMOKE_AUDIT"
  "$CT64_ORIGINAL" "$CT64_BACKUP" "$PRISM_ORIGINAL" "$PRISM_BACKUP"
  "$PRISM_SOURCE_ORIGINAL" "$PRISM_SOURCE_BACKUP"
  "$SELECTION_ORIGINAL" "$SELECTION_BACKUP" "$STATS" "$SOURCE_CONFIG"
  "$CONTEXT_NPZ" "$REGISTRY" "$ACTIVITY_WEIGHT" "$CCC_CALIBRATION"
)
if [[ -e "$LOG" ]]; then
  echo "FATAL retry log exists; refusing to mix a new audit with stale output: $LOG" >&2
  exit 2
fi
for required in "${required_files[@]}"; do
  if [[ ! -f "$required" ]]; then
    log "FATAL missing required artifact: $required"
    exit 2
  fi
done

if [[ -e "$BOUNDARY_SMOKE_RUN" ]]; then
  log "FATAL boundary-smoke directory exists; refusing overwrite: $BOUNDARY_SMOKE_RUN"
  exit 2
fi
if [[ -e "$BOUNDARY_SMOKE_AUDIT" ]]; then
  log "FATAL boundary-smoke audit exists; refusing overwrite: $BOUNDARY_SMOKE_AUDIT"
  exit 2
fi
if [[ -e "$RETRY_RUN" ]]; then
  log "FATAL retry run directory exists; refusing overwrite/resume: $RETRY_RUN"
  exit 2
fi
ACTIVE_GPU_PIDS="$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits | sort -u | tr '\n' ' ' | xargs || true)"
if [[ -n "$ACTIVE_GPU_PIDS" ]]; then
  log "FATAL GPU compute processes active: $ACTIVE_GPU_PIDS"
  exit 2
fi

verify_sha "$FULL_REVIEW" "$EXPECTED_FULL_CONFIG_SHA" "reviewed retry config"
verify_sha "$BOUNDARY_SMOKE_CONFIG" "$EXPECTED_BOUNDARY_SMOKE_CONFIG_SHA" "boundary-smoke config"
verify_sha "$DESIGN_DOC" "$EXPECTED_DESIGN_SHA" "design document"
verify_sha "$FAILED_CHECKPOINT" "$EXPECTED_FAILED_CHECKPOINT_SHA" "failed-run epoch-1 checkpoint"
verify_sha "$FAILED_AUTHORIZED_CONFIG" "$EXPECTED_FAILED_AUTHORIZED_CONFIG_SHA" "failed-run authorized config"
verify_sha "$PRIOR_SMOKE_AUDIT" "$EXPECTED_PRIOR_SMOKE_AUDIT_SHA" "prior gate smoke audit"
verify_sha "$CT64_ORIGINAL" "$EXPECTED_CT64_SHA" "CT64 canonical"
verify_sha "$CT64_BACKUP" "$EXPECTED_CT64_SHA" "CT64 backup"
verify_sha "$PRISM_ORIGINAL" "$EXPECTED_PRISM_SHA" "PRISM selected warm-start"
verify_sha "$PRISM_BACKUP" "$EXPECTED_PRISM_SHA" "PRISM warm-start backup"
verify_sha "$PRISM_SOURCE_ORIGINAL" "$EXPECTED_PRISM_SOURCE_SHA" "PRISM source config"
verify_sha "$PRISM_SOURCE_BACKUP" "$EXPECTED_PRISM_SOURCE_SHA" "PRISM source-config backup"
verify_sha "$SELECTION_ORIGINAL" "$EXPECTED_SELECTION_SHA" "PRISM selection summary"
verify_sha "$SELECTION_BACKUP" "$EXPECTED_SELECTION_SHA" "PRISM selection-summary backup"
verify_sha "$STATS" "$EXPECTED_STATS_SHA" "module-rescue train64 statistics"
verify_sha "$SOURCE_CONFIG" "$EXPECTED_SOURCE_CONFIG_SHA" "v4 source config"
verify_sha "$CONTEXT_NPZ" "$EXPECTED_CONTEXT_SHA" "PRISM context"
verify_sha "$REGISTRY" "$EXPECTED_REGISTRY_SHA" "module registry"
verify_sha "$ACTIVITY_WEIGHT" "$EXPECTED_ACTIVITY_WEIGHT_SHA" "activity weights"
verify_sha "$CCC_CALIBRATION" "$EXPECTED_CCC_CALIBRATION_SHA" "CCC calibration"
verify_sha "$REPO/src/kmlee_bam/training/adversarial_trainer.py" "$EXPECTED_CADENCE_PATCH_SHA" "DDP cadence patch"

"$PYTHON" - "$FULL_REVIEW" "$BOUNDARY_SMOKE_CONFIG" "$REPO" "$RETRY_RUN" "$BOUNDARY_SMOKE_RUN" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

full_path, smoke_path, repo = map(Path, sys.argv[1:4])
retry_run, smoke_run = sys.argv[4:6]
full = json.loads(full_path.read_text(encoding="utf-8"))
smoke = json.loads(smoke_path.read_text(encoding="utf-8"))

if full["warm_start"]["init_weights_path"] != full["experiment_manifest"]["warm_start"]["path"]:
    raise SystemExit("FATAL retry no longer starts from the original selected PRISM checkpoint")
if full["warm_start"]["expected_sha256"] != "9bee5aba99e38453499d6cd2333831509b811945a33e77c2b381eb3ea655c555":
    raise SystemExit("FATAL original PRISM warm-start SHA changed")
if full["train"]["out_dir"] != retry_run:
    raise SystemExit("FATAL retry output directory mismatch")
if full["train"].get("use_rank0_eval") is not False:
    raise SystemExit("FATAL full retry must use exact six-rank distributed validation")
if full["train"].get("ddp_timeout_minutes") != 30:
    raise SystemExit("FATAL full retry DDP timeout must be 30 minutes")
if full["train"].get("epochs") != 36 or full["train"].get("early_stopping_min_epochs") != 18:
    raise SystemExit("FATAL full retry training horizon changed")
lineage = full["experiment_manifest"].get("retry_lineage", {})
if lineage.get("mode") != "fresh_restart_from_original_prism_warm_start" or lineage.get("failed_checkpoint_loaded") is not False:
    raise SystemExit("FATAL retry lineage permits loading the failed checkpoint")
if full.get("_launch_guard", {}).get("launch_allowed") is not False:
    raise SystemExit("FATAL reviewed retry config must remain launch-locked")

st = smoke["train"]
if (
    st.get("out_dir") != smoke_run
    or st.get("epochs") != 3
    or st.get("use_rank0_eval") is not True
    or st.get("max_train_steps_per_epoch") != 34
    or st.get("max_eval_steps") != 1
):
    raise SystemExit("FATAL dedicated boundary-smoke topology/cadence changed")
if smoke["learned_generator_count"].get("start_epoch") != 3:
    raise SystemExit("FATAL boundary smoke must keep epochs 1-2 all-on and activate the gate in epoch 3")
if smoke.get("_launch_guard", {}).get("launch_allowed") is not True:
    raise SystemExit("FATAL boundary smoke is not launch-authorized")
contract = smoke["experiment_manifest"].get("boundary_smoke_contract", {})
if contract.get("patched_diagnostic_collective_epoch") != 3 or contract.get("patched_diagnostic_collective_train_step") != 32:
    raise SystemExit("FATAL causal cadence regression contract changed")

for config in (full, smoke):
    for relative, expected in config["experiment_manifest"]["runtime_source_sha256"].items():
        observed = hashlib.sha256((repo / relative).read_bytes()).hexdigest()
        if observed != expected:
            raise SystemExit(f"FATAL runtime source mismatch: {relative}")
PY

log "PREFLIGHT_OK immutable failed run + CT64/PRISM backups verified; starting causal six-GPU boundary smoke"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$BOUNDARY_SMOKE_PORT" \
  -m kmlee_bam.training.run_current --config "$BOUNDARY_SMOKE_CONFIG" \
  2>&1 | tee -a "$LOG"

"$PYTHON" - "$BOUNDARY_SMOKE_CONFIG" "$BOUNDARY_SMOKE_RUN" "$BOUNDARY_SMOKE_AUDIT" "$LOG" "$EXPECTED_BOUNDARY_SMOKE_CONFIG_SHA" <<'PY'
import datetime
import hashlib
import json
from pathlib import Path
import sys
import torch

from kmlee_bam.training.learned_generator_count import HardBinaryConcreteGeneratorGate
from kmlee_bam.training.runner_base import load_config

config_path, run_dir, audit_path, log_path = map(Path, sys.argv[1:5])
expected_config_sha = sys.argv[5]
observed_config_sha = hashlib.sha256(config_path.read_bytes()).hexdigest()
if observed_config_sha != expected_config_sha:
    raise SystemExit("FATAL boundary-smoke config changed while running")
checkpoint_path = run_dir / "checkpoint_last.pt"
history_path = run_dir / "history.json"
count_path = run_dir / "joint_generator_count_latest.json"
rescue_manifest_path = run_dir / "prism_module_rescue_manifest.json"
for path in (checkpoint_path, history_path, count_path, rescue_manifest_path):
    if not path.is_file():
        raise SystemExit(f"FATAL boundary-smoke artifact missing: {path}")
checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
history = json.loads(history_path.read_text(encoding="utf-8"))
if int(checkpoint.get("epoch", -1)) != 3 or len(history) != 3:
    raise SystemExit("FATAL boundary smoke did not complete all three epochs")
if int(checkpoint.get("global_step", -1)) != 102:
    raise SystemExit(f"FATAL boundary-smoke global_step is not 102: {checkpoint.get('global_step')}")
log_text = log_path.read_text(encoding="utf-8", errors="replace")
if "[diagnostics | adversary | epoch 3]" not in log_text:
    raise SystemExit("FATAL epoch-3 synchronized diagnostic collective was not observed")
state = checkpoint.get("system_state_dict", {})
if int(state.get("generator_count_epoch_state", torch.tensor(-1)).item()) != 3:
    raise SystemExit("FATAL gated epoch 3 was not checkpointed")
required_state = (
    "generator_count_gate.log_alpha",
    "generator_count_gate.protected_mask",
    "generator_count_objective.dual",
    "generator_count_objective.violation_ema",
    "generator_count_objective.ema_initialized",
)
for key in required_state:
    if key not in state:
        raise SystemExit(f"FATAL persistent joint state missing: {key}")
for key in (
    "generator_count_gate.log_alpha",
    "generator_count_objective.dual",
    "generator_count_objective.violation_ema",
):
    if not bool(torch.isfinite(state[key]).all()):
        raise SystemExit(f"FATAL non-finite joint state: {key}")
if not bool(state["generator_count_objective.ema_initialized"].item()):
    raise SystemExit("FATAL constrained generator-count objective was never activated")

cfg = load_config(str(config_path))
gate = HardBinaryConcreteGeneratorGate(
    int(state["generator_count_gate.log_alpha"].numel()),
    config=cfg.learned_generator_count,
    protected_mask=state["generator_count_gate.protected_mask"],
)
initial = gate.log_alpha.detach().clone()
gate.load_state_dict(
    {
        "log_alpha": state["generator_count_gate.log_alpha"],
        "protected_mask": state["generator_count_gate.protected_mask"],
    }
)
if torch.equal(gate.log_alpha.detach(), initial):
    raise SystemExit("FATAL gated epoch 3 did not update generator logits")
temperature = gate.temperature(1.0)
mask = gate.deterministic_mask(temperature=temperature)
if not bool(((mask == 0) | (mask == 1)).all()):
    raise SystemExit("FATAL deterministic generator architecture is not exactly binary")
k = int(mask.sum().item())
if not 32 <= k <= int(gate.num_generators):
    raise SystemExit(f"FATAL hard generator count violates bounds: {k}")
if not bool(mask[state["generator_count_gate.protected_mask"].bool()].all()):
    raise SystemExit("FATAL a protected singleton generator was disabled")

serialized_history = json.dumps(history)
for metric in (
    "metric/generator_candidate_full_nll",
    "metric/generator_candidate_isolated_nll",
    "metric/generator_expected_active",
    "loss/generator_count",
):
    if metric not in serialized_history:
        raise SystemExit(f"FATAL boundary-smoke history lacks {metric}")
count_payload = json.loads(count_path.read_text(encoding="utf-8"))
if count_payload.get("candidate_k_grid_used") is not False or count_payload.get("official_test_used_for_count") is not False:
    raise SystemExit("FATAL boundary smoke lost its no-grid/sealed-test contract")
rescue_manifest = json.loads(rescue_manifest_path.read_text(encoding="utf-8"))
if "generator_count_gate.log_alpha" not in rescue_manifest.get("rescue_parameter_names", []):
    raise SystemExit("FATAL module-rescue allowlist omitted the gate logits")
audit = {
    "schema_version": "kmlee_bam.prism_v5_sv6_boundary_smoke.v1",
    "status": "passed",
    "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "config_sha256": observed_config_sha,
    "checkpoint_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
    "checkpoint_epoch": 3,
    "global_step": 102,
    "rank0_eval_steps_per_epoch": 1,
    "train_steps_per_epoch": 34,
    "epoch3_synchronized_diagnostic_observed": True,
    "generator_count": k,
    "expected_active": float(gate.expected_active_count(temperature=temperature)),
    "gate_logit_max_abs_change": float((gate.log_alpha.detach() - initial).abs().max()),
    "dual": state["generator_count_objective.dual"].tolist(),
    "test_used": False,
}
audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
print(json.dumps(audit, sort_keys=True))
PY

verify_sha "$FULL_REVIEW" "$EXPECTED_FULL_CONFIG_SHA" "reviewed retry config after smoke"
verify_sha "$BOUNDARY_SMOKE_CONFIG" "$EXPECTED_BOUNDARY_SMOKE_CONFIG_SHA" "boundary-smoke config after smoke"
verify_sha "$FAILED_CHECKPOINT" "$EXPECTED_FAILED_CHECKPOINT_SHA" "failed-run checkpoint after smoke"
verify_sha "$CT64_BACKUP" "$EXPECTED_CT64_SHA" "CT64 backup after smoke"
verify_sha "$PRISM_BACKUP" "$EXPECTED_PRISM_SHA" "PRISM backup after smoke"

mkdir "$RETRY_RUN"
"$PYTHON" - "$FULL_REVIEW" "$AUTHORIZED_CONFIG" "$APPROVAL_RECEIPT" "$BOUNDARY_SMOKE_AUDIT" "$EXPECTED_FULL_CONFIG_SHA" <<'PY'
import datetime
import hashlib
import json
from pathlib import Path
import sys

reviewed_path, authorized_path, receipt_path, smoke_path = map(Path, sys.argv[1:5])
reviewed_sha = sys.argv[5]
config = json.loads(reviewed_path.read_text(encoding="utf-8"))
smoke = json.loads(smoke_path.read_text(encoding="utf-8"))
if smoke.get("status") != "passed" or smoke.get("checkpoint_epoch") != 3 or smoke.get("test_used") is not False:
    raise SystemExit("FATAL full retry authorization requires the passed gated boundary smoke")
smoke_sha = hashlib.sha256(smoke_path.read_bytes()).hexdigest()
config["_launch_guard"] = {
    "launch_allowed": True,
    "approved_parent_config_sha256": reviewed_sha,
    "approval_basis": "user_delegated_retry_after_fail_closed_boundary_smoke_2026-08-06",
    "boundary_smoke_audit_sha256": smoke_sha,
}
manifest = config["experiment_manifest"]
manifest["full_training_authorized"] = True
manifest["approval_received"] = True
manifest["approval_basis"] = "user_delegated_retry_after_fail_closed_boundary_smoke_2026-08-06"
manifest["launch_performed"] = True
authorized_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
authorized_sha = hashlib.sha256(authorized_path.read_bytes()).hexdigest()
receipt = {
    "schema_version": "kmlee_bam.prism_v5_retry1_approval_receipt.v1",
    "reviewed_config_sha256": reviewed_sha,
    "authorized_config_sha256": authorized_sha,
    "boundary_smoke_audit_sha256": smoke_sha,
    "failed_epoch1_checkpoint_sha256": "5cc9d5f2f31edc0fc5975b6df749c677fe42af3e0c24983764f50714f39fe4b4",
    "failed_checkpoint_loaded": False,
    "approval_basis": "user_delegated_retry_after_fail_closed_boundary_smoke_2026-08-06",
    "authorized_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}
receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
PY

log "BOUNDARY_SMOKE_PASS; starting fresh full retry from original PRISM epoch-11 warm-start, distributed validation, max_epochs=36"
"$PYTHON" -m torch.distributed.run \
  --nnodes=1 --nproc_per_node=6 \
  --rdzv-backend=c10d --rdzv-endpoint="localhost:$FULL_PORT" \
  -m kmlee_bam.training.run_current --config "$AUTHORIZED_CONFIG" \
  2>&1 | tee -a "$LOG"

log "TRAINING_COMPLETE; validation-selected checkpoint remains subject to posthoc audit"
