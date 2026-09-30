from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import torch
import torch.distributed as dist
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler

# anndata 0.x emits this dependency deprecation once per DDP process and once
# per data-loader worker.  It does not describe the model or data quality and
# can bury the experiment diagnostics in dozens of identical lines.
warnings.filterwarnings(
    "ignore",
    message=r"Importing read_elem from .* is deprecated.*",
    category=FutureWarning,
    module=r"anndata\.utils",
)

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_SRC_DIR = _PROJECT_ROOT / "src"
for _p in (_PROJECT_ROOT, _SRC_DIR):
    _s = str(_p)
    if _s not in sys.path:
        sys.path.insert(0, _s)

from kmlee_bam.data.dataset_with_covariates import V3OrdinalScDataset
from kmlee_bam.data.reference_anchor_sampler import ReferenceAnchorBatchSampler
from kmlee_bam.data.grouped_batch_sampler import GroupedBatchSampler
from kmlee_bam.model.system import V3OrdinalBAMSystem
from kmlee_bam.objectives.hierarchical_ordinal import HierarchicalOrdinalConfig
from kmlee_bam.objectives.ordinal_uncertainty import (
    GlobalReferenceBankConfig,
    OrdinalBalanceConfig,
    UncertaintyCalibrationConfig,
    make_ordinal_class_weights,
)
from kmlee_bam.objectives.state_usage import V5StateUsageConfig
from kmlee_bam.objectives.tech_adversary import ConditionalTechAdversary
from kmlee_bam.objectives.reference_celltype_adversary import ReferenceCelltypeAdversary
from kmlee_bam.objectives.sex_adversary import SexAdversary
from kmlee_bam.objectives.aux_reference_adversary import AuxReferenceAdversaryConfig
from kmlee_bam.objectives.uncertainty_residual import V5UncertaintyResidualConfig
from kmlee_bam.model.reference_anchored_prior import ReferenceAnchoredPriorConfig
from kmlee_bam.training.reference_anchored_trainer import V6Trainer
from kmlee_bam.training.pathology_aware_uncertainty_trainer import V7aTrainer
from kmlee_bam.training.adaptive_subgroup_trainer import V8Trainer
from kmlee_bam.objectives.latent_knn_pi_tech import LatentKNNPiTechConfig
from kmlee_bam.objectives.thinning import ThinningConfig, ZeroReliabilityConfig
from kmlee_bam.training.tech_invariance import TechInvarianceConfig
from kmlee_bam.training.sealed_donor_allowlist import (
    bind_allowlist_to_dataset,
    dataset_train_donor_names,
    load_sealed_donor_allowlist,
    refit_age_normalization_for_allowlist,
    restrict_dataset_to_allowlist,
    validate_w44_source_checkpoint,
)
from kmlee_bam.training.architecture_donor_split import DonorIndexView
from kmlee_bam.training.generator_exact_training_backend import (
    read_json_object as _read_generator_exact_json,
    validate_fixed_mask_training_job,
)
from kmlee_bam.objectives.subgroup_cvar import SubgroupCVaRConfig
from kmlee_bam.objectives.variance_floor import VarianceFloorConfig
from kmlee_bam.objectives.ranking_loss import HighMarginRankConfig
from kmlee_bam.training import runner_base as base
from kmlee_bam.uncertainty import (
    ANCOVAConfig as V7aANCOVAConfig,
    DifficultyProfileConfig as V7aDifficultyProfileConfig,
)
from kmlee_bam.objectives.uncertainty_alignment import (
    UncertaintyAlignmentConfig as V7aAlignmentConfig,
)

# Keep unpatched base builders. install_hooks() replaces base.* with the
# integrated KMLEE-BAM wrappers below, so wrappers must call these saved
# originals to avoid recursive self-calls.
_BASE_BUILD_LOADERS = base.build_loaders
_BASE_BUILD_SYSTEM = base.build_system
_BASE_BUILD_CRITERION = base.build_criterion


CURRENT_SETTINGS: Dict[str, Any] = {}
ORDINAL_CLASS_WEIGHTS: Optional[torch.Tensor] = None
N_CELLTYPES: Optional[int] = None
N_DONORS: Optional[int] = None
D_Z: Optional[int] = None
N_GENES: Optional[int] = None
N_BINS: Optional[int] = None
N_TECH: Optional[int] = None
# v7a runtime state (populated by _prepare_runtime_state if v7a enabled)
V7A_DONOR_PATHOLOGY: Optional[torch.Tensor] = None        # [D, P]
V7A_DONOR_CELLTYPE_MASK: Optional[torch.Tensor] = None    # [D, T] bool
V7A_DONOR_PATHOLOGY_VALID: Optional[torch.Tensor] = None  # [D] bool
V7A_TRAIN_DATASET = None  # cached for profile init
PI_TECH_TRAIN_DATASET = None  # train split only; used by the detached latent bank
PI_TECH_SEX_LINKED_GENE_INDICES: tuple[int, ...] = ()


_RUNTIME_ENV_KEYS = (
    "KMLEE_GROUPED_BLOCK_SIZE",
    "KMLEE_PATH_CONSISTENCY",
    "KMLEE_PATH_CONS_MINCELLS",
    "WORLD_SIZE",
    "RANK",
    "LOCAL_RANK",
    "CUDA_VISIBLE_DEVICES",
    "OMP_NUM_THREADS",
    "NCCL_DEBUG",
    "NCCL_P2P_DISABLE",
    "NCCL_IB_DISABLE",
    "NCCL_SOCKET_IFNAME",
    "NCCL_ASYNC_ERROR_HANDLING",
    "TORCH_NCCL_ASYNC_ERROR_HANDLING",
    "PYTORCH_CUDA_ALLOC_CONF",
    "LD_LIBRARY_PATH",
    "KMLEE_LAUNCH_SCRIPT",
)

_RUNTIME_SOURCE_FILES = (
    "src/kmlee_bam/model/latent_nuisance_projection.py",
    "src/kmlee_bam/model/lie_ordinal_decoder.py",
    "src/kmlee_bam/model/system.py",
    "src/kmlee_bam/objectives/celltype_alignment.py",
    "src/kmlee_bam/training/core_trainer.py",
    "src/kmlee_bam/training/learned_generator_count.py",
    "src/kmlee_bam/training/learned_pathology_rank.py",
    "src/kmlee_bam/training/prism_module_rescue_training.py",
    "src/kmlee_bam/training/run_current.py",
    "src/kmlee_bam/training/runner_base.py",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validated_sha256(value: Any, *, label: str) -> str:
    """Return a normalized SHA256 digest or reject ambiguous provenance."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be a 64-character hex digest")
    digest = value.lower()
    if len(digest) != 64 or any(
        character not in "0123456789abcdef" for character in digest
    ):
        raise ValueError(f"{label} must be a 64-character hex digest")
    return digest


def _generator_mask_sha256(
    *,
    candidate_count: int,
    selected_local_generator_ids: tuple[int, ...],
) -> str:
    """Hash the exact candidate-order binary mask as raw uint8 bytes."""

    mask = bytearray(int(candidate_count))
    for local_id in selected_local_generator_ids:
        mask[int(local_id)] = 1
    return hashlib.sha256(bytes(mask)).hexdigest()


def _assert_fixed_generator_integrity(
    system,
    *,
    stage: str,
    distributed: bool = False,
) -> dict[str, Any] | None:
    """Fail closed if a committed exact generator mask changes in memory.

    This audit is deliberately dormant for ordinary and discovery/search runs.
    A fixed architecture is recognized only by the exact-confirmation
    provenance installed by :func:`_warm_start_from_checkpoint`.

    When ``distributed`` is true, every rank first reduces local validity and
    mask width, then compares the complete binary mask only if widths agree.
    Local errors are delayed until after the collectives so a single corrupt
    rank cannot strand its peers in NCCL.
    """

    provenance = getattr(system, "fixed_generator_result_provenance", None)
    if provenance is None:
        return None
    if not isinstance(provenance, dict):
        raise RuntimeError("fixed-generator provenance must be an object")
    installed_provenance_sha256 = getattr(
        system,
        "_fixed_generator_integrity_provenance_sha256",
        None,
    )
    observed_provenance_sha256 = _canonical_json_sha256(provenance)
    if (
        installed_provenance_sha256 is not None
        and observed_provenance_sha256 != installed_provenance_sha256
    ):
        provenance_error = (
            "fixed-generator provenance changed after architecture commit"
        )
    else:
        provenance_error = None

    decoder = getattr(system, "decoder", None)
    local_errors: list[str] = []
    if provenance_error is not None:
        local_errors.append(provenance_error)
    if decoder is None or not hasattr(decoder, "generator_active_mask"):
        local_errors.append("decoder generator_active_mask is missing")

    decoder_count_raw = getattr(decoder, "n_generators", None)
    decoder_count_valid = bool(
        not isinstance(decoder_count_raw, bool)
        and isinstance(decoder_count_raw, int)
        and int(decoder_count_raw) > 0
    )
    candidate_count_raw = provenance.get("candidate_count")
    if (
        isinstance(candidate_count_raw, bool)
        or not isinstance(candidate_count_raw, int)
        or int(candidate_count_raw) <= 0
    ):
        local_errors.append("provenance candidate_count is invalid")
        candidate_count = int(decoder_count_raw) if decoder_count_valid else 1
    else:
        candidate_count = int(candidate_count_raw)
    if decoder_count_valid and candidate_count != int(decoder_count_raw):
        local_errors.append(
            "provenance candidate_count differs from decoder n_generators"
        )
        # DDP collectives must use the model's common architecture width, not
        # a possibly corrupted per-rank provenance value.
        candidate_count = int(decoder_count_raw)

    selected_raw = provenance.get("selected_local_generator_ids")
    if not isinstance(selected_raw, list) or any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in (selected_raw if isinstance(selected_raw, list) else [])
    ):
        local_errors.append("provenance selected generator ids are invalid")
        selected_ids: tuple[int, ...] = ()
    else:
        selected_ids = tuple(int(value) for value in selected_raw)
    if (
        not selected_ids
        or len(set(selected_ids)) != len(selected_ids)
        or any(value < 0 or value >= candidate_count for value in selected_ids)
    ):
        local_errors.append("provenance selected generator ids are inconsistent")

    expected_cpu = torch.zeros(candidate_count, dtype=torch.int32)
    if selected_ids and all(
        0 <= value < candidate_count for value in selected_ids
    ):
        expected_cpu[list(selected_ids)] = 1

    mask = (
        decoder.generator_active_mask.detach()
        if decoder is not None and hasattr(decoder, "generator_active_mask")
        else torch.empty(0)
    )
    observed_cpu = torch.full((candidate_count,), -1, dtype=torch.int32)
    if mask.ndim != 1 or int(mask.numel()) != candidate_count:
        local_errors.append(
            "decoder generator mask shape differs from fixed candidate_count"
        )
    else:
        finite = bool(torch.isfinite(mask).all().item())
        binary = finite and bool(((mask == 0) | (mask == 1)).all().item())
        if not binary:
            local_errors.append("decoder generator mask is not exactly binary")
        else:
            observed_cpu.copy_(mask.to(device="cpu", dtype=torch.int32))

    expected_digest = provenance.get("mask_sha256")
    if not isinstance(expected_digest, str):
        local_errors.append("provenance mask_sha256 is missing")
        expected_digest = ""
    observed_digest = hashlib.sha256(
        observed_cpu.clamp_min(0).to(dtype=torch.uint8).numpy().tobytes()
    ).hexdigest()
    if observed_digest != expected_digest:
        local_errors.append(
            "decoder binary mask SHA256 differs from fixed provenance"
        )
    if not torch.equal(observed_cpu, expected_cpu):
        local_errors.append("decoder generator mask differs from selected ids")

    # This is the effective multiplicative route used by both the Lie action
    # and affine translation.  An inactive row may retain learned parameters,
    # but its route coefficient must remain exactly zero.
    inactive_route = observed_cpu * (1 - expected_cpu)
    inactive_route_nonzero = int((inactive_route != 0).sum().item())
    if inactive_route_nonzero:
        local_errors.append("an inactive generator route is nonzero")

    ddp_agreement = True
    distributed_ready = bool(
        dist.is_available() and dist.is_initialized()
    )
    use_collectives = bool(
        distributed and distributed_ready and dist.get_world_size() > 1
    )
    if distributed and not distributed_ready:
        local_errors.append(
            "distributed fixed-mask audit requested without an initialized group"
        )
    if use_collectives:
        if mask.numel() > 0:
            collective_device = mask.device
        elif dist.get_backend() == "nccl":
            collective_device = torch.device("cuda", torch.cuda.current_device())
        else:
            collective_device = torch.device("cpu")
        local_valid = torch.tensor(
            [0 if local_errors else 1],
            dtype=torch.int32,
            device=collective_device,
        )
        dist.all_reduce(local_valid, op=dist.ReduceOp.MIN)

        mask_count = torch.tensor(
            [int(observed_cpu.numel())],
            dtype=torch.int64,
            device=collective_device,
        )
        mask_count_min = mask_count.clone()
        mask_count_max = mask_count.clone()
        dist.all_reduce(mask_count_min, op=dist.ReduceOp.MIN)
        dist.all_reduce(mask_count_max, op=dist.ReduceOp.MAX)
        same_mask_width = bool(torch.equal(mask_count_min, mask_count_max))
        if same_mask_width:
            collective_mask = observed_cpu.to(device=collective_device)
            mask_min = collective_mask.clone()
            mask_max = collective_mask.clone()
            dist.all_reduce(mask_min, op=dist.ReduceOp.MIN)
            dist.all_reduce(mask_max, op=dist.ReduceOp.MAX)
            ddp_agreement = bool(torch.equal(mask_min, mask_max))
        else:
            ddp_agreement = False
        if int(local_valid.item()) != 1:
            local_errors.append(
                "at least one DDP rank failed its local fixed-mask audit"
            )
        if not ddp_agreement:
            local_errors.append("DDP ranks disagree on the fixed generator mask")

    if local_errors:
        raise RuntimeError(
            f"fixed-generator integrity audit failed at {stage}: "
            + "; ".join(dict.fromkeys(local_errors))
        )

    return {
        "stage": str(stage),
        "candidate_count": candidate_count,
        "hard_active": int(expected_cpu.sum().item()),
        "mask_sha256": observed_digest,
        "provenance_sha256": observed_provenance_sha256,
        "inactive_route_nonzero_count": inactive_route_nonzero,
        "ddp_ranks_agree": bool(ddp_agreement),
    }


def _install_fixed_generator_integrity_hook(system) -> None:
    """Attach a non-state-dict audit callback to an exact fixed architecture."""

    provenance = getattr(system, "fixed_generator_result_provenance", None)
    if not isinstance(provenance, dict):
        raise RuntimeError(
            "cannot install fixed-generator integrity hook without provenance"
        )
    system._fixed_generator_integrity_provenance_sha256 = (
        _canonical_json_sha256(provenance)
    )

    def _hook(*, stage: str, distributed: bool = False):
        return _assert_fixed_generator_integrity(
            system,
            stage=stage,
            distributed=distributed,
        )

    system._fixed_generator_integrity_hook = _hook


def _canonical_json_sha256(value: Any) -> str:
    """Hash JSON with the exact serialization used by exact confirmation."""

    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _git_output(*arguments: str) -> bytes | None:
    try:
        result = subprocess.run(
            ["git", *arguments],
            cwd=_PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout if result.returncode == 0 else None


def _build_runtime_manifest(config_path: str | Path) -> dict[str, Any]:
    """Resolve launch-time controls without copying unrelated environment data."""

    config_path = Path(config_path).resolve()
    environment = {key: os.environ.get(key) for key in _RUNTIME_ENV_KEYS}
    launch_script_value = environment["KMLEE_LAUNCH_SCRIPT"]
    launch_script = (
        Path(launch_script_value).resolve()
        if launch_script_value
        else None
    )
    git_commit_raw = _git_output("rev-parse", "HEAD")
    git_status = _git_output("status", "--porcelain=v1", "--untracked-files=normal")
    source_sha256 = {}
    for relative in _RUNTIME_SOURCE_FILES:
        source_path = _PROJECT_ROOT / relative
        if source_path.is_file():
            source_sha256[relative] = _sha256_file(source_path)

    cuda_available = bool(torch.cuda.is_available())
    gpu_names: list[str] = []
    if cuda_available:
        try:
            gpu_names = [
                torch.cuda.get_device_name(index)
                for index in range(torch.cuda.device_count())
            ]
        except Exception:
            gpu_names = []
    try:
        nccl_version = torch.cuda.nccl.version() if cuda_available else None
    except Exception:
        nccl_version = None

    return {
        "schema_version": "kmlee_bam.runtime_manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "launch": {
            "argv": list(sys.argv),
            "config_path": str(config_path),
            "config_sha256": _sha256_file(config_path),
            "script_path": (
                str(launch_script) if launch_script is not None else None
            ),
            "script_sha256": (
                _sha256_file(launch_script)
                if launch_script is not None and launch_script.is_file()
                else None
            ),
        },
        "resolved_execution_controls": {
            "grouped_block_size": int(
                environment["KMLEE_GROUPED_BLOCK_SIZE"] or "0"
            ),
            "path_consistency": float(
                environment["KMLEE_PATH_CONSISTENCY"] or "0"
            ),
            "path_consistency_min_cells": int(
                environment["KMLEE_PATH_CONS_MINCELLS"] or "2"
            ),
            "world_size": int(environment["WORLD_SIZE"] or "1"),
            "rank": int(environment["RANK"] or "0"),
            "local_rank": int(environment["LOCAL_RANK"] or "0"),
        },
        "environment": environment,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "hostname": platform.node(),
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
            "cudnn": torch.backends.cudnn.version(),
            "nccl": nccl_version,
        },
        "accelerator": {
            "cuda_available": cuda_available,
            "visible_device_count": len(gpu_names),
            "gpu_names": gpu_names,
        },
        "source": {
            "git_commit": (
                git_commit_raw.decode("utf-8", errors="replace").strip()
                if git_commit_raw is not None
                else None
            ),
            "working_tree_dirty": (
                bool(git_status) if git_status is not None else None
            ),
            "git_status_entry_count": (
                len(git_status.splitlines()) if git_status is not None else None
            ),
            "git_status_sha256": (
                hashlib.sha256(git_status).hexdigest()
                if git_status is not None
                else None
            ),
            "critical_file_sha256": source_sha256,
        },
    }


def _write_runtime_manifest(out_dir: Path, config_path: str | Path) -> None:
    manifest = _build_runtime_manifest(config_path)
    with (out_dir / "runtime_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _section(name: str) -> Dict[str, Any]:
    return dict(CURRENT_SETTINGS.get(name, {}))


def _dataset_section() -> Dict[str, Any]:
    return dict(_section("v3").get("dataset", {}))


def _load_generator_exact_training_job() -> dict[str, Any] | None:
    section = _section("generator_exact_training")
    if not bool(section.get("enabled", False)):
        return None
    if section.get("fresh_optimizer") is not True or section.get(
        "resume_forbidden"
    ) is not True:
        raise ValueError(
            "generator_exact_training requires fresh_optimizer=true and "
            "resume_forbidden=true"
        )
    path = section.get("job_path")
    if not isinstance(path, str) or not path:
        raise ValueError("generator_exact_training.job_path is required")
    job = validate_fixed_mask_training_job(
        _read_generator_exact_json(path, label="generator exact training job")
    )
    if section.get("expected_job_sha256") != job["job_sha256"]:
        raise ValueError("generator exact training job SHA256 differs from config")
    split_path = Path(str(section.get("donor_split_manifest_path", ""))).resolve()
    if str(split_path) != str(
        Path(job["source_paths"]["donor_split_manifest"]).resolve()
    ):
        raise ValueError("generator exact donor-split path differs from job")
    if _sha256_file(split_path) != section.get("donor_split_manifest_sha256"):
        raise ValueError("generator exact donor-split SHA256 mismatch")
    return job


def _refit_age_for_exact_donors(dataset, donor_ids: tuple[int, ...]) -> dict[str, Any]:
    donor_age = np.asarray(dataset.donor_age_years, dtype=np.float64)
    selected = donor_age[np.asarray(donor_ids, dtype=np.int64)]
    observed = selected[np.isfinite(selected)]
    if observed.size < 2:
        raise ValueError("exact training donors have fewer than two observed ages")
    mean = float(np.mean(observed))
    standard_deviation = float(np.std(observed))
    if not np.isfinite(standard_deviation) or standard_deviation < 1e-6:
        raise ValueError("exact training donor age standard deviation is degenerate")
    dataset.age_train_mean = mean
    dataset.age_train_std = standard_deviation
    dataset.age_z = np.where(
        np.asarray(dataset.age_valid, dtype=bool),
        (np.asarray(dataset.age_years, dtype=np.float64) - mean)
        / standard_deviation,
        0.0,
    ).astype(np.float32)
    return {
        "age_mean": mean,
        "age_standard_deviation": standard_deviation,
        "observed_age_donor_count": int(observed.size),
    }


def _restrict_dataset_for_generator_exact_training(train_ds, job):
    """Apply the W44 or W+R54 union as an actual row-level dataset view."""

    split_path = job["source_paths"]["donor_split_manifest"]
    train64_names = dataset_train_donor_names(train_ds)
    weight = load_sealed_donor_allowlist(
        split_path,
        partition="W44",
        expected_train_donor_names=train64_names,
    )
    donor_ids = list(bind_allowlist_to_dataset(train_ds, weight))
    declared_names = list(weight.donor_names)
    fit_partition = "W44"
    if job["donor_partition"] == "confirmation":
        ranking = load_sealed_donor_allowlist(
            split_path,
            partition="R10",
            expected_train_donor_names=train64_names,
        )
        donor_ids.extend(bind_allowlist_to_dataset(train_ds, ranking))
        declared_names.extend(ranking.donor_names)
        fit_partition = "W_PLUS_R54"
    if declared_names != list(job["training_donor_names"]):
        raise ValueError("sealed donor names differ from exact training job")
    if len(set(donor_ids)) != len(donor_ids):
        raise ValueError("exact training donor partitions overlap")
    age = _refit_age_for_exact_donors(train_ds, tuple(donor_ids))
    view = DonorIndexView(train_ds, donor_ids)
    observed_names = dataset_train_donor_names(view)
    if set(observed_names) != set(declared_names):
        raise RuntimeError("exact optimizer rows differ from sealed donor union")
    evidence: dict[str, Any] = {
        "schema_version": "kmlee_bam.generator_exact_row_allowlist.v1",
        "fit_partition": fit_partition,
        "declared_donor_names": declared_names,
        "declared_donor_names_sha256": job["training_donor_names_sha256"],
        "observed_optimizer_donor_names": list(observed_names),
        "observed_optimizer_donor_count": len(observed_names),
        "exact_training_job_sha256": job["job_sha256"],
        "donor_split_manifest_sha256": weight.manifest_sha256,
        "row_level_donor_allowlist_enforced": True,
        "official_validation_used": False,
        "official_test_used": False,
        **age,
    }
    evidence["evidence_sha256"] = _canonical_json_sha256(evidence)
    view.generator_exact_training_allowlist = evidence
    return view


def _ref_section() -> Dict[str, Any]:
    return dict(_section("v3").get("donor_balanced_ref_center", {}))


def _sampler_section() -> Dict[str, Any]:
    return dict(_section("v3").get("reference_anchor_sampler", {}))


def _tech_adv_section() -> Dict[str, Any]:
    return dict(_section("v3").get("tech_adversary", {}))


def _celltype_adv_section() -> Dict[str, Any]:
    """v3.reference_celltype_adversary — reference-only GRL adversary that erases
    cell-type identity from z_perp on control cells. OFF by default; absent
    section ⇒ celltype adversary is None and the forward/loss path is unchanged.
    See doc memory: reference-only celltype adversary (z_perp celltype leak)."""
    return dict(_section("v3").get("reference_celltype_adversary", {}))


def _sex_adv_section() -> Dict[str, Any]:
    """v3.sex_adversary — v31 GRL adversary that erases SEX from z_clean (soft, complements the
    hard projection's refill). OFF by default; absent section ⇒ sex adversary None, forward unchanged."""
    return dict(_section("v3").get("sex_adversary", {}))


def _ordinal_section() -> Dict[str, Any]:
    return dict(_section("v4").get("ordinal_balance", {}))


def _global_ref_section() -> Dict[str, Any]:
    return dict(_section("v4").get("global_reference_bank", {}))


def _v4_unc_section() -> Dict[str, Any]:
    return dict(_section("v4").get("uncertainty_calibration", {}))


def _reference_prior_section() -> Dict[str, Any]:
    return dict(_section("v6").get("reference_anchored_prior", {}))


def _state_usage_section() -> Dict[str, Any]:
    return dict(_section("v6").get("state_usage", {}))


def _score_mixer_section() -> Dict[str, Any]:
    """v6.score_mixer — gate balance + reference-safe + coeff zero-mean
    regularizers for the ScoreResidualMixer + Lie coefficient head.
    See doc/v6_min_mixer_design_2026-05-22.md."""
    return dict(_section("v6").get("score_mixer", {}))


def _hierarchical_ordinal_section() -> Dict[str, Any]:
    """v8.hierarchical_ordinal — L1 zero/nonzero + L2 group | nonzero + L3 EMD."""
    return dict(_section("v8").get("hierarchical_ordinal", {}))


def _variance_floor_section() -> Dict[str, Any]:
    """v8.variance_floor — within-celltype state-variance floor (Phase 2)."""
    return dict(_section("v8").get("variance_floor", {}))


def _high_margin_rank_section() -> Dict[str, Any]:
    """v8.high_margin_rank — high-tier (tier4) margin ranking loss."""
    return dict(_section("v8").get("high_margin_rank", {}))


def _uncertainty_residual_section() -> Dict[str, Any]:
    return dict(_section("v6").get("uncertainty_residual", {}))


def _v7a_section() -> Dict[str, Any]:
    return dict(_section("v7a").get("pathology_aware_uncertainty", {}))


def _v7a_enabled() -> bool:
    return bool(_v7a_section().get("enabled", False))


def _zero_reliability_section() -> Dict[str, Any]:
    return dict(_section("zero_reliability"))


def _extract_numeric_pathology(values, axis_name: str) -> np.ndarray:
    """
    Convert a pandas-like Series of pathology labels to numeric values.

    SEA-AD typically encodes axes as strings like "Braak III", "Thal 4",
    "Frequent". We use the following rules:
      1. If already numeric, cast to float.
      2. If integer-like string ("3", "0"), parse as float.
      3. If formatted as "<prefix> <integer>" or "<prefix> <Roman>", extract
         the trailing number (Roman → arabic for Braak).
      4. CERAD: {None, Absent, Sparse, Moderate, Frequent} → {0, 0, 1, 2, 3}.
      5. ADNC: {Not, Low, Intermediate, High} → {0, 1, 2, 3}.
      6. Anything else → NaN (donor gets marked invalid).
    """
    import re

    roman_to_int = {
        "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6,
        "0": 0,
    }
    cerad_map = {
        "absent": 0.0, "none": 0.0,
        "sparse": 1.0,
        "moderate": 2.0,
        "frequent": 3.0,
    }
    adnc_map = {
        "not aa": 0.0, "not": 0.0, "none": 0.0,
        "low": 1.0,
        "intermediate": 2.0,
        "high": 3.0,
    }

    raw = np.asarray(values)
    out = np.full(raw.shape[0], np.nan, dtype=np.float64)
    axis_lower = axis_name.lower()

    for i, v in enumerate(raw):
        try:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool):
                out[i] = float(v)
                continue
            s = str(v).strip()
            if s == "" or s.lower() in {"nan", "none", "missing", "n/a"}:
                continue
            # CERAD-specific
            if "cerad" in axis_lower and s.lower() in cerad_map:
                out[i] = cerad_map[s.lower()]
                continue
            # ADNC-specific
            if "adnc" in axis_lower and s.lower() in adnc_map:
                out[i] = adnc_map[s.lower()]
                continue
            # Trailing integer
            m = re.search(r"(-?\d+(?:\.\d+)?)\s*$", s)
            if m:
                out[i] = float(m.group(1))
                continue
            # Trailing Roman numeral (uppercase)
            m_roman = re.search(r"([IVX]+)\s*$", s)
            if m_roman:
                roman = m_roman.group(1)
                if roman in roman_to_int:
                    out[i] = float(roman_to_int[roman])
                    continue
            # Lookup whole string in CERAD/ADNC tables as last resort
            if s.lower() in cerad_map:
                out[i] = cerad_map[s.lower()]
                continue
            if s.lower() in adnc_map:
                out[i] = adnc_map[s.lower()]
                continue
            # Failure: leave as NaN
        except Exception:
            continue
    return out


def _build_v7a_donor_metadata(train_ds) -> Dict[str, torch.Tensor]:
    """
    Build per-donor pathology covariates + (donor, celltype) mask from the
    training dataset's underlying zarr obs DataFrame.
    """
    from kmlee_bam.data.ordinal_dataset import read_zarr_dataframe

    v7a_cfg = _v7a_section()
    ancova_cfg = dict(v7a_cfg.get("ancova", {}))
    axes = list(ancova_cfg.get("pathology_axes", ["Braak_stage", "Thal_phase", "CERAD_score"]))

    # Map config-style names to actual obs column names (SEA-AD uses spaces).
    name_aliases = {
        "Braak_stage": "Braak stage",
        "Thal_phase": "Thal phase",
        "CERAD_score": "CERAD score",
        "ADNC": "Overall AD neuropathological Change",
    }

    obs = read_zarr_dataframe(train_ds.root, "obs")
    donor_vocab = train_ds.donor_vocab
    donor_to_id = {str(v): i for i, v in enumerate(donor_vocab.tolist())}
    n_donors = len(donor_vocab)
    n_axes = len(axes)

    # Donor-level row picker: first row of obs for each donor.
    donor_col_name = train_ds.donor_key
    donor_series = obs[donor_col_name].astype(str).values
    selected_rows = np.asarray(train_ds.row_idx, dtype=np.int64)
    # Build first-occurrence index from optimizer-eligible train rows only.
    # The underlying zarr contains validation/test metadata as well; using its
    # unrestricted first occurrence would silently populate PHU ANCOVA state
    # for sealed donors even when their dataset objects are never opened.
    first_idx = {}
    for i in selected_rows.tolist():
        d = donor_series[i]
        if d not in first_idx:
            first_idx[d] = i

    pathology = np.full((n_donors, n_axes), np.nan, dtype=np.float64)
    valid = np.zeros(n_donors, dtype=bool)

    for d_str, donor_idx in donor_to_id.items():
        row_idx = first_idx.get(d_str, None)
        if row_idx is None:
            continue
        donor_ok = True
        for a_idx, axis_name in enumerate(axes):
            obs_col = name_aliases.get(axis_name, axis_name)
            if obs_col not in obs.columns:
                donor_ok = False
                continue
            val_arr = _extract_numeric_pathology(obs[obs_col].iloc[[row_idx]].values, obs_col)
            v = float(val_arr[0])
            if not np.isfinite(v):
                donor_ok = False
            else:
                pathology[donor_idx, a_idx] = v
        valid[donor_idx] = donor_ok and not np.any(np.isnan(pathology[donor_idx]))

    # Build (donor, celltype) mask. The dataset already exposes per-cell
    # integer arrays (celltype_ids, donor_ids) aligned to its row_idx, so use
    # those directly rather than re-parsing obs columns (the latter approach
    # needs a celltype_key attribute that the dataset does not actually expose).
    n_celltypes = len(train_ds.spec.celltype_vocab)
    ct_series = np.asarray(getattr(train_ds, "celltype_ids", None), dtype=np.int64) \
        if getattr(train_ds, "celltype_ids", None) is not None else None
    dn_series = np.asarray(getattr(train_ds, "donor_ids", None), dtype=np.int64) \
        if getattr(train_ds, "donor_ids", None) is not None else None

    if ct_series is None or dn_series is None:
        # Fallback: derive from obs DataFrame using config-style column names.
        ct_vocab = list(train_ds.spec.celltype_vocab)
        ct_to_id = {str(v): i for i, v in enumerate(ct_vocab)}
        # Common SEA-AD subclass column names.
        candidate_ct_cols = ["Subclass", "subclass", "celltype", "cell_type"]
        ct_col = next((c for c in candidate_ct_cols if c in obs.columns), None)
        if ct_col is None:
            raise KeyError(
                "Could not locate celltype column in obs; expected one of "
                + ", ".join(candidate_ct_cols)
            )
        ct_series = np.array(
            [ct_to_id.get(str(v), -1) for v in obs[ct_col].astype(str).values],
            dtype=np.int64,
        )
        dn_series = np.array(
            [donor_to_id.get(str(v), -1) for v in donor_series],
            dtype=np.int64,
        )

    if ct_series.shape[0] != selected_rows.shape[0]:
        ct_series = ct_series[selected_rows]
    if dn_series.shape[0] != selected_rows.shape[0]:
        dn_series = dn_series[selected_rows]
    donor_celltype_mask = np.zeros((n_donors, n_celltypes), dtype=bool)
    for d_id, t_id in zip(dn_series, ct_series):
        if 0 <= d_id < n_donors and 0 <= t_id < n_celltypes:
            donor_celltype_mask[d_id, t_id] = True

    if int(os.environ.get("RANK", "0")) == 0:
        n_valid = int(valid.sum())
        n_total_pairs = int(donor_celltype_mask.sum())
        print(
            f"[v7a-meta] donors={n_donors} celltypes={n_celltypes} axes={n_axes} "
            f"valid_donors={n_valid}/{n_donors} "
            f"(donor,celltype)_pairs={n_total_pairs}",
            flush=True,
        )

    # Replace NaN with 0 so the tensor is well-defined; the valid mask gates
    # ANCOVA usage so the zeros never enter the regression.
    pathology_safe = np.where(np.isnan(pathology), 0.0, pathology)
    return {
        "donor_pathology": torch.tensor(pathology_safe, dtype=torch.float32),
        "donor_celltype_mask": torch.tensor(donor_celltype_mask, dtype=torch.bool),
        "donor_pathology_valid": torch.tensor(valid, dtype=torch.bool),
    }


def set_settings(raw: Dict[str, Any]) -> None:
    CURRENT_SETTINGS.clear()
    CURRENT_SETTINGS.update(raw)


def _verify_sealed_model_inputs() -> None:
    """Fail before dataset I/O when a config-pinned model input has drifted."""

    expected_inputs = _section("experiment_manifest").get(
        "critical_model_inputs", {}
    )
    if not isinstance(expected_inputs, dict):
        raise ValueError("experiment_manifest.critical_model_inputs must be an object")
    checks = (
        (
            "precision_medicine",
            "context_npz_path",
            "precision_context_npz_sha256",
        ),
        (
            "module_tokenizer",
            "registry_json_path",
            "module_registry_json_sha256",
        ),
        (
            "module_tokenizer",
            "activity_weight_path",
            "activity_weight_npz_sha256",
        ),
    )
    for section_name, path_key, sha_key in checks:
        section = _section(section_name)
        expected_raw = expected_inputs.get(sha_key)
        if expected_raw is None:
            continue
        expected = _validated_sha256(
            expected_raw,
            label=f"experiment_manifest.critical_model_inputs.{sha_key}",
        )
        path_raw = section.get(path_key)
        if not isinstance(path_raw, str) or not path_raw:
            raise ValueError(
                f"experiment_manifest.critical_model_inputs.{sha_key} requires "
                f"{section_name}.{path_key}"
            )
        path = Path(path_raw).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"sealed model input is missing: {path}")
        observed = _sha256_file(path)
        if observed != expected:
            raise RuntimeError(
                f"sealed model input SHA256 mismatch for {section_name}.{path_key}: "
                f"{observed} != {expected}"
            )


def _enforce_launch_guard() -> None:
    """Reject review-only configs before output directories or data are opened."""

    launch_guard = _section("_launch_guard")
    if launch_guard and launch_guard.get("launch_allowed") is not True:
        raise RuntimeError(
            "this generated config is review-only; set "
            "_launch_guard.launch_allowed=true only after explicit user approval"
        )


def build_datasets(cfg):
    _enforce_launch_guard()
    _verify_sealed_model_inputs()
    ds_cfg = _dataset_section()
    common = dict(
        zarr_path=cfg.data.zarr_path,
        spec_path=cfg.data.spec_path,
        matrix_path=cfg.data.matrix_path,
        counts_layer=cfg.data.counts_layer,
        prefer_raw=cfg.data.prefer_raw,
        split_key=cfg.data.split_key,
        celltype_key=cfg.data.celltype_key,
        batch_key=cfg.data.batch_key,
        return_log_input=cfg.data.return_log_input,
        reference_scaler_path=cfg.data.reference_scaler_path,
        return_x_gene_scalar=cfg.data.return_x_gene_scalar,
        x_gene_scalar_clip=cfg.data.x_gene_scalar_clip,
        x_gene_scalar_eps=cfg.data.x_gene_scalar_eps,
        max_csr_data_span=cfg.data.max_csr_data_span,
        max_span_to_selected_nnz_ratio=cfg.data.max_span_to_selected_nnz_ratio,
        donor_key=ds_cfg.get("donor_key", "donor_id"),
        depth_key=ds_cfg.get("depth_key", "Number of UMIs"),
        detected_genes_key=ds_cfg.get("detected_genes_key", "Genes detected"),
    )
    # Thinning (v8 zero-origin): only the train split needs the thinned view —
    # the consistency/supervision losses are training-only. Off ⇒ unchanged.
    thin_cfg = _section("thinning")
    thin_on = bool(thin_cfg.get("enabled", False))
    thin_rho = float(thin_cfg.get("rho", 0.5))
    thin_dataset_kwargs: Dict[str, Any] = {}
    if thin_on:
        thin_dataset_kwargs.update(
            return_thinned=True,
            thin_rho=thin_rho,
            thin_counts_zarr_path=thin_cfg.get("counts_zarr_path"),
            thin_counts_matrix_path=str(thin_cfg.get("counts_matrix_path", "X")),
            thin_library_size_key=str(thin_cfg.get("library_size_key", "Number of UMIs")),
        )
    train_ds = V3OrdinalScDataset(split="train", **thin_dataset_kwargs, **common)
    exact_job = _load_generator_exact_training_job()
    if exact_job is not None:
        if bool(_section("canonical_w44_source").get("enabled", False)):
            raise ValueError(
                "generator_exact_training and canonical_w44_source are mutually exclusive"
            )
        if bool(cfg.train.early_stopping) or bool(
            cfg.train.require_validation_for_early_stopping
        ):
            raise ValueError(
                "generator_exact_training cannot open official validation for "
                "early stopping"
            )
        if bool(cfg.train.save_best) or not bool(cfg.train.skip_final_test):
            raise ValueError(
                "generator_exact_training requires save_best=false and "
                "skip_final_test=true"
            )
        if int(cfg.train.epochs) != int(
            exact_job["training_contract"]["epochs"]
        ) or int(cfg.train.seed) != int(exact_job["training_contract"]["seed"]):
            raise ValueError("exact training epoch/seed schedule differs from job")
        train_ds = _restrict_dataset_for_generator_exact_training(
            train_ds, exact_job
        )
        if int(os.environ.get("RANK", "0")) == 0:
            evidence = train_ds.generator_exact_training_allowlist
            print(
                "[generator exact training] optimizer donors="
                f"{evidence['observed_optimizer_donor_count']} "
                f"partition={evidence['fit_partition']} "
                "official validation/test unopened; fixed mask is committed "
                "before optimizer construction",
                flush=True,
            )
        return train_ds, None, None
    w44_cfg = _section("canonical_w44_source")
    if bool(w44_cfg.get("enabled", False)):
        if bool(cfg.learned_generator_count.enabled) or bool(
            cfg.generator_budget_search.enabled
        ):
            raise ValueError(
                "canonical_w44_source cannot run a generator search while "
                "fitting source weights"
            )
        if bool(cfg.train.early_stopping) or bool(
            cfg.train.require_validation_for_early_stopping
        ):
            raise ValueError(
                "canonical_w44_source uses a fixed epoch budget and cannot "
                "open official validation for early stopping"
            )
        if bool(cfg.train.save_best) or not bool(cfg.train.skip_final_test):
            raise ValueError(
                "canonical_w44_source requires save_best=false and "
                "skip_final_test=true"
            )
        split_path = w44_cfg.get("donor_split_manifest_path")
        age_path = w44_cfg.get("age_normalization_out_path")
        if not split_path or not age_path:
            raise ValueError(
                "canonical_w44_source requires donor_split_manifest_path and "
                "age_normalization_out_path"
            )
        allowlist = load_sealed_donor_allowlist(
            split_path,
            partition="W44",
            expected_train_donor_names=dataset_train_donor_names(train_ds),
        )
        warm_start_path = _section("warm_start").get("init_weights_path")
        warm_start_checkpoint_sha256 = None
        if warm_start_path:
            if not bool(w44_cfg.get("allow_attested_w44_warm_start", False)):
                raise ValueError(
                    "canonical_w44_source warm start requires explicit "
                    "allow_attested_w44_warm_start=true"
                )
            warm_start_checkpoint_sha256, _ = (
                validate_w44_source_checkpoint(warm_start_path, allowlist)
            )
            configured_expected = _section("warm_start").get(
                "expected_sha256"
            )
            if (
                configured_expected is None
                or str(configured_expected).lower()
                != warm_start_checkpoint_sha256
            ):
                raise ValueError(
                    "attested W44 warm-start checkpoint must also match "
                    "warm_start.expected_sha256"
                )
        age_normalization = refit_age_normalization_for_allowlist(
            train_ds, allowlist
        )
        train_ds = restrict_dataset_to_allowlist(train_ds, allowlist)
        observed_names = dataset_train_donor_names(train_ds)
        if observed_names != allowlist.donor_names:
            raise RuntimeError("optimizer rows do not equal the sealed W44 donors")
        source_provenance: dict[str, Any] = {
            "schema_version": "kmlee_bam.canonical_w44_source_checkpoint.v4",
            "fit_partition": "W44",
            "observed_optimizer_donor_names": list(observed_names),
            "observed_optimizer_donor_count": len(observed_names),
            "allowed_fit_donor_names": list(allowlist.donor_names),
            "allowed_fit_donor_sha256": allowlist.donor_names_sha256,
            "donor_split_manifest_path": allowlist.manifest_path,
            "donor_split_manifest_sha256": allowlist.manifest_sha256,
            "age_normalization_manifest_sha256": age_normalization[
                "manifest_sha256"
            ],
            "row_level_allowlist_enforced": True,
            "official_validation_used": False,
            "official_test_used": False,
        }
        if warm_start_checkpoint_sha256 is not None:
            source_provenance["attested_w44_warm_start_checkpoint_sha256"] = (
                warm_start_checkpoint_sha256
            )
        source_provenance["manifest_sha256"] = _canonical_json_sha256(
            source_provenance
        )
        train_ds.canonical_w44_source_provenance = source_provenance
        age_output = Path(age_path).resolve()
        if int(os.environ.get("RANK", "0")) == 0:
            age_output.parent.mkdir(parents=True, exist_ok=True)
            encoded = json.dumps(
                age_normalization, indent=2, sort_keys=True
            ) + "\n"
            if age_output.exists() and age_output.read_text(
                encoding="utf-8"
            ) != encoded:
                raise FileExistsError(
                    "refusing to overwrite a different W44 age artifact: "
                    f"{age_output}"
                )
            age_output.write_text(encoded, encoding="utf-8")
        if int(os.environ.get("RANK", "0")) == 0:
            print(
                "[canonical W44 source] optimizer donors=44; official "
                "validation/test datasets were not opened; "
                f"split_sha256={allowlist.manifest_sha256}",
                flush=True,
            )
        return train_ds, None, None
    val_ds = None
    test_ds = None
    try:
        val_ds = V3OrdinalScDataset(split="val", **common)
    except Exception:
        val_ds = None
    if not bool(cfg.train.skip_final_test):
        try:
            test_ds = V3OrdinalScDataset(split="test", **common)
        except Exception:
            test_ds = None
    elif int(os.environ.get("RANK", "0")) == 0:
        print(
            "[sealed-test] skip_final_test=true; official test dataset was not opened",
            flush=True,
        )
    return train_ds, val_ds, test_ds


def build_loaders(
    cfg,
    train_ds,
    val_ds,
    test_ds,
    *,
    distributed: bool,
):
    # v32: grouped (donor x celltype) block batches for the stochastic-consistency loss.
    # Activated by env KMLEE_GROUPED_BLOCK_SIZE>0 (default off -> existing paths untouched).
    grouped_block = int(os.environ.get("KMLEE_GROUPED_BLOCK_SIZE", "0"))
    if grouped_block > 0:
        return _build_grouped_loaders(
            cfg, train_ds, val_ds, test_ds, distributed=distributed, block_size=grouped_block
        )
    sampler_cfg = _sampler_section()
    use_anchor_sampler = bool(sampler_cfg.get("enabled", False))

    if not use_anchor_sampler:
        return _BASE_BUILD_LOADERS(
            cfg,
            train_ds,
            val_ds,
            test_ds,
            distributed=distributed,
        )

    eval_batch_size = cfg.loader.eval_batch_size or cfg.loader.batch_size
    world_size = int(os.environ.get("WORLD_SIZE", "1")) if distributed else 1
    rank = int(os.environ.get("RANK", "0")) if distributed else 0

    train_sampler = ReferenceAnchorBatchSampler(
        train_ds,
        batch_size=int(cfg.loader.batch_size),
        num_replicas=world_size,
        rank=rank,
        seed=int(cfg.train.seed),
        drop_last=bool(cfg.data.drop_last_train),
        anchor_probability=float(sampler_cfg.get("anchor_probability", 1.0)),
        ref_cells_per_donor=int(sampler_cfg.get("ref_cells_per_donor", 1)),
        ref_donors_per_batch=int(sampler_cfg.get("ref_donors_per_batch", 2)),
        fill_from_reference=bool(sampler_cfg.get("fill_from_reference", False)),
    )
    train_loader = DataLoader(
        train_ds,
        batch_sampler=train_sampler,
        num_workers=cfg.data.num_workers,
        pin_memory=cfg.data.pin_memory,
        persistent_workers=(
            cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
        ),
    )
    if rank == 0:
        print(
            "[kmlee-sampler] ReferenceAnchorBatchSampler active "
            f"batches/epoch/rank={len(train_sampler)} "
            f"batch_size={cfg.loader.batch_size} "
            f"anchor_probability={train_sampler.anchor_probability} "
            f"ref_cells_per_donor={train_sampler.ref_cells_per_donor} "
            f"ref_donors_per_batch={train_sampler.ref_donors_per_batch} "
            f"valid_celltypes={len(train_sampler.valid_celltypes)}",
            flush=True,
        )

    use_rank0_eval = bool(getattr(cfg.train, "use_rank0_eval", False))
    val_sampler = None
    test_sampler = None
    if distributed and use_rank0_eval:
        val_loader = None
        test_loader = None
    else:
        if distributed:
            if val_ds is not None:
                val_sampler = DistributedSampler(
                    val_ds,
                    num_replicas=world_size,
                    rank=rank,
                    shuffle=False,
                    drop_last=False,
                )
            if test_ds is not None:
                test_sampler = DistributedSampler(
                    test_ds,
                    num_replicas=world_size,
                    rank=rank,
                    shuffle=False,
                    drop_last=False,
                )
        val_loader = (
            None
            if val_ds is None
            else DataLoader(
                val_ds,
                batch_size=eval_batch_size,
                shuffle=False,
                sampler=val_sampler,
                num_workers=cfg.data.num_workers,
                pin_memory=cfg.data.pin_memory,
                persistent_workers=(
                    cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
                ),
                drop_last=False,
            )
        )
        test_loader = (
            None
            if test_ds is None
            else DataLoader(
                test_ds,
                batch_size=eval_batch_size,
                shuffle=False,
                sampler=test_sampler,
                num_workers=cfg.data.num_workers,
                pin_memory=cfg.data.pin_memory,
                persistent_workers=(
                    cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
                ),
                drop_last=False,
            )
        )
    return train_loader, val_loader, test_loader, train_sampler, val_sampler, test_sampler


def _build_grouped_loaders(
    cfg,
    train_ds,
    val_ds,
    test_ds,
    *,
    distributed: bool,
    block_size: int,
):
    """v32: grouped (celltype, donor) block batches so the stochastic-consistency loss sees
    DENSE aggregates (n_blocks x block_size cells) instead of ~89% singletons. Train uses
    GroupedBatchSampler; eval loaders are identical to the other paths (DistributedSampler)."""
    eval_batch_size = cfg.loader.eval_batch_size or cfg.loader.batch_size
    world_size = int(os.environ.get("WORLD_SIZE", "1")) if distributed else 1
    rank = int(os.environ.get("RANK", "0")) if distributed else 0

    train_sampler = GroupedBatchSampler(
        train_ds,
        batch_size=int(cfg.loader.batch_size),
        block_size=int(block_size),
        num_replicas=world_size,
        rank=rank,
        seed=int(cfg.train.seed),
        drop_last=bool(cfg.data.drop_last_train),
    )
    train_loader = DataLoader(
        train_ds,
        batch_sampler=train_sampler,
        num_workers=cfg.data.num_workers,
        pin_memory=cfg.data.pin_memory,
        persistent_workers=(
            cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
        ),
    )
    if rank == 0:
        print(
            "[kmlee-sampler] GroupedBatchSampler active "
            f"block_size={train_sampler.block_size} n_blocks={train_sampler.n_blocks} "
            f"effective_batch={train_sampler.effective_batch} "
            f"batches/epoch/rank={len(train_sampler)} "
            f"groups_eligible={train_sampler.n_groups}/{train_sampler.n_groups_total} "
            f"groups_excluded={train_sampler.n_groups_excluded} "
            f"cells_eligible={train_sampler.n_cells_eligible:,}/{len(train_ds):,} "
            f"cells_excluded={train_sampler.n_cells_excluded:,}",
            flush=True,
        )

    use_rank0_eval = bool(getattr(cfg.train, "use_rank0_eval", False))
    val_sampler = None
    test_sampler = None
    if distributed and use_rank0_eval:
        val_loader = None
        test_loader = None
    else:
        if distributed:
            if val_ds is not None:
                val_sampler = DistributedSampler(
                    val_ds,
                    num_replicas=world_size,
                    rank=rank,
                    shuffle=False,
                    drop_last=False,
                )
            if test_ds is not None:
                test_sampler = DistributedSampler(
                    test_ds,
                    num_replicas=world_size,
                    rank=rank,
                    shuffle=False,
                    drop_last=False,
                )
        val_loader = (
            None
            if val_ds is None
            else DataLoader(
                val_ds,
                batch_size=eval_batch_size,
                shuffle=False,
                sampler=val_sampler,
                num_workers=cfg.data.num_workers,
                pin_memory=cfg.data.pin_memory,
                persistent_workers=(
                    cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
                ),
                drop_last=False,
            )
        )
        test_loader = (
            None
            if test_ds is None
            else DataLoader(
                test_ds,
                batch_size=eval_batch_size,
                shuffle=False,
                sampler=test_sampler,
                num_workers=cfg.data.num_workers,
                pin_memory=cfg.data.pin_memory,
                persistent_workers=(
                    cfg.data.persistent_workers if cfg.data.num_workers > 0 else False
                ),
                drop_last=False,
            )
        )
    return train_loader, val_loader, test_loader, train_sampler, val_sampler, test_sampler


def _n_ordinal_bins(train_ds) -> int:
    edges = getattr(train_ds.spec, "edges", None)
    if edges is not None and np.ndim(edges) == 2:
        return int(edges.shape[1]) + 1
    if hasattr(train_ds.spec, "n_bins"):
        return int(train_ds.spec.n_bins)
    raise AttributeError("Could not infer number of ordinal bins from train_ds.spec.")


def _sample_train_bin_counts(train_ds, ordinal_cfg: Dict[str, Any]) -> torch.Tensor:
    explicit = ordinal_cfg.get("bin_counts", None)
    if explicit is not None:
        counts = np.asarray(explicit, dtype=np.float64)
        if counts.ndim != 1:
            raise ValueError("v4.ordinal_balance.bin_counts must be a 1D list.")
        return torch.tensor(counts, dtype=torch.float32)

    n_bins = _n_ordinal_bins(train_ds)
    max_cells = int(ordinal_cfg.get("max_cells_for_bin_stats", 4096))
    max_cells = max(1, min(max_cells, len(train_ds)))
    indices = np.linspace(0, len(train_ds) - 1, num=max_cells, dtype=np.int64)
    bin_stats_log_every = int(ordinal_cfg.get("bin_stats_log_every", 1024))
    if os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower() == "prism_experiment":
        bin_stats_log_every = max(bin_stats_log_every, 2048)

    counts = np.zeros(n_bins, dtype=np.float64)
    for i, idx in enumerate(indices, start=1):
        y_ord = train_ds[int(idx)]["y_ord"].numpy()
        counts += np.bincount(y_ord.astype(np.int64), minlength=n_bins)[:n_bins]
        if (
            int(os.environ.get("RANK", "0")) == 0
            and i % bin_stats_log_every == 0
        ):
            print(f"[kmlee-ordinal] sampled {i}/{max_cells} train cells for bin counts", flush=True)
    return torch.tensor(counts, dtype=torch.float32)


def _resolve_pi_tech_sex_linked_gene_indices(
    pi_tech_cfg: Dict[str, Any], train_ds
) -> tuple[int, ...]:
    """Resolve an explicit sex-linked whitelist against the sealed gene order."""

    gene_names = tuple(str(value) for value in train_ds.spec.gene_names)
    resolved: set[int] = set()
    for raw_index in pi_tech_cfg.get("sex_linked_gene_indices", []):
        if isinstance(raw_index, bool):
            raise ValueError("sex_linked_gene_indices must contain integer indices")
        index = int(raw_index)
        if index < 0 or index >= len(gene_names):
            raise ValueError(
                f"sex-linked gene index {index} is outside [0, {len(gene_names)})"
            )
        resolved.add(index)

    requested_names = tuple(
        str(value) for value in pi_tech_cfg.get("sex_linked_gene_names", [])
    )
    if requested_names:
        name_to_index: Dict[str, int] = {}
        duplicates: set[str] = set()
        symbol_values = getattr(train_ds, "gene_symbols", gene_names)
        gene_symbols = tuple(str(value) for value in symbol_values)
        if len(gene_symbols) != len(gene_names):
            raise ValueError("dataset gene symbols do not match decoder gene order")
        for index, aliases in enumerate(zip(gene_names, gene_symbols)):
            for name in set(aliases):
                if name in name_to_index and name_to_index[name] != index:
                    duplicates.add(name)
                else:
                    name_to_index[name] = index
        ambiguous = sorted(set(requested_names).intersection(duplicates))
        if ambiguous:
            raise ValueError(
                "sex-linked whitelist contains ambiguous gene names: "
                + ", ".join(ambiguous[:10])
            )
        missing = [name for name in requested_names if name not in name_to_index]
        if missing and not bool(
            pi_tech_cfg.get("allow_missing_sex_linked_genes", False)
        ):
            raise ValueError(
                "sex-linked whitelist genes are absent from the decoder gene order: "
                + ", ".join(missing[:20])
            )
        resolved.update(
            name_to_index[name] for name in requested_names if name in name_to_index
        )
        expected_resolution = pi_tech_cfg.get(
            "sex_linked_expected_resolution", {}
        )
        if expected_resolution:
            if not isinstance(expected_resolution, dict):
                raise ValueError("sex_linked_expected_resolution must be an object")
            for alias, contract in expected_resolution.items():
                if not isinstance(contract, dict):
                    raise ValueError(
                        f"sex-linked resolution for {alias!r} must be an object"
                    )
                if str(alias) not in name_to_index:
                    raise ValueError(
                        f"sex-linked resolution alias is absent: {alias!r}"
                    )
                observed_index = int(name_to_index[str(alias)])
                expected_index = int(contract["decoder_index"])
                expected_ensembl = str(contract["ensembl_id"])
                if (
                    observed_index != expected_index
                    or gene_names[observed_index] != expected_ensembl
                ):
                    raise ValueError(
                        "sex-linked gene resolution contract changed for "
                        f"{alias!r}: observed index/id="
                        f"{observed_index}/{gene_names[observed_index]}, "
                        f"expected={expected_index}/{expected_ensembl}"
                    )
    return tuple(sorted(resolved))


def _prepare_runtime_state(cfg, train_ds) -> None:
    global ORDINAL_CLASS_WEIGHTS, N_CELLTYPES, N_DONORS, D_Z, N_GENES, N_BINS, N_TECH
    global V7A_DONOR_PATHOLOGY, V7A_DONOR_CELLTYPE_MASK, V7A_DONOR_PATHOLOGY_VALID
    global V7A_TRAIN_DATASET, PI_TECH_TRAIN_DATASET
    global PI_TECH_SEX_LINKED_GENE_INDICES
    ordinal_cfg = _ordinal_section()
    counts = _sample_train_bin_counts(train_ds, ordinal_cfg)
    ORDINAL_CLASS_WEIGHTS = make_ordinal_class_weights(
        counts,
        mode=str(ordinal_cfg.get("class_weight_mode", "sqrt_inverse")),
        cap=float(ordinal_cfg.get("class_weight_cap", 6.0)),
        zero_bin_weight_floor=float(ordinal_cfg.get("zero_bin_weight_floor", 0.2)),
        eps=float(ordinal_cfg.get("eps", 1e-6)),
    )
    N_CELLTYPES = len(train_ds.spec.celltype_vocab)
    N_TECH = int(len(train_ds.spec.batch_vocab))
    donor_vocab = getattr(train_ds, "donor_vocab", None)
    N_DONORS = int(len(donor_vocab)) if donor_vocab is not None else None
    D_Z = int(cfg.encoder.d_z)

    # Gene + bin counts for v7a difficulty profile.
    n_bins = _n_ordinal_bins(train_ds)
    N_BINS = int(n_bins)
    # n_genes from a single sample probe.
    try:
        probe = train_ds[0]
        if "y_ord" in probe:
            N_GENES = int(probe["y_ord"].shape[-1])
        elif "x" in probe:
            N_GENES = int(probe["x"].shape[-1])
        else:
            N_GENES = None
    except Exception:
        N_GENES = None

    pi_tech_cfg = _section("pi_tech")
    if bool(pi_tech_cfg.get("enabled", False)) and str(
        pi_tech_cfg.get("mode", "in_batch")
    ) == "same_celltype_donor_bank":
        # `train_ds` is constructed with split="train" by the base runner.
        # Keeping this reference separate prevents a validation/test dataset
        # from ever being substituted implicitly during trainer construction.
        PI_TECH_TRAIN_DATASET = train_ds
        PI_TECH_SEX_LINKED_GENE_INDICES = (
            _resolve_pi_tech_sex_linked_gene_indices(pi_tech_cfg, train_ds)
        )
    else:
        PI_TECH_TRAIN_DATASET = None
        PI_TECH_SEX_LINKED_GENE_INDICES = ()

    # v7a donor-level metadata (only if v7a is enabled).
    if _v7a_enabled() and N_DONORS is not None and N_CELLTYPES is not None:
        try:
            meta = _build_v7a_donor_metadata(train_ds)
            V7A_DONOR_PATHOLOGY = meta["donor_pathology"]
            V7A_DONOR_CELLTYPE_MASK = meta["donor_celltype_mask"]
            V7A_DONOR_PATHOLOGY_VALID = meta["donor_pathology_valid"]
            V7A_TRAIN_DATASET = train_ds
        except Exception as e:
            if int(os.environ.get("RANK", "0")) == 0:
                print(f"[v7a-meta] WARNING: failed to build donor metadata: {e}", flush=True)
            V7A_DONOR_PATHOLOGY = None
            V7A_DONOR_CELLTYPE_MASK = None
            V7A_DONOR_PATHOLOGY_VALID = None
            V7A_TRAIN_DATASET = None

    if int(os.environ.get("RANK", "0")) == 0:
        frac = counts / counts.sum().clamp_min(1.0)
        print(
            "[kmlee-ordinal] train bin fractions="
            + ",".join(f"{float(x):.4f}" for x in frac.tolist()),
            flush=True,
        )
        print(
            "[kmlee-ordinal] class weights="
            + ",".join(f"{float(x):.4f}" for x in ORDINAL_CLASS_WEIGHTS.tolist()),
            flush=True,
        )


def _train_tech_counts(train_ds, n_tech: int) -> np.ndarray:
    batch_ids = getattr(train_ds, "batch_ids", None)
    row_idx = getattr(train_ds, "row_idx", None)
    if batch_ids is None or row_idx is None:
        raise AttributeError("train_ds must expose batch_ids and row_idx for tech counts.")
    ids = np.asarray(batch_ids)[np.asarray(row_idx, dtype=np.int64)]
    return np.bincount(ids.astype(np.int64), minlength=int(n_tech)).astype(np.float64)


def _tech_class_weights_from_config(train_ds, n_tech: int, adv_cfg: Dict[str, Any]) -> Optional[torch.Tensor]:
    mode = str(adv_cfg.get("class_weight_mode", "none")).lower()
    if mode in {"none", "off", "false", "0"}:
        return None
    if mode not in {"inverse", "sqrt_inverse"}:
        raise ValueError("tech_adversary.class_weight_mode must be 'none', 'inverse', or 'sqrt_inverse'.")

    counts = _train_tech_counts(train_ds, n_tech)
    eps = float(adv_cfg.get("class_weight_eps", 1.0))
    if eps <= 0:
        raise ValueError("class_weight_eps must be positive.")
    smoothed = np.maximum(counts, eps)
    weights = smoothed.sum() / (float(n_tech) * smoothed)
    if mode == "sqrt_inverse":
        weights = np.sqrt(weights)

    cap = adv_cfg.get("class_weight_cap", None)
    if cap is not None:
        cap_f = float(cap)
        if cap_f <= 0:
            raise ValueError("class_weight_cap must be positive when provided.")
        weights = np.minimum(weights, cap_f)

    if bool(adv_cfg.get("normalize_class_weights", True)):
        mean_weight = float(weights.mean())
        if mean_weight > 0:
            weights = weights / mean_weight

    weights = weights.astype(np.float32)
    if int(os.environ.get("RANK", "0")) == 0:
        names = [str(x) for x in train_ds.spec.batch_vocab.tolist()]
        counts_msg = ",".join(f"{name}:{int(count)}" for name, count in zip(names, counts))
        weights_msg = ",".join(f"{name}:{float(weight):.4f}" for name, weight in zip(names, weights))
        print(
            "[kmlee-adversary] tech class weighting "
            f"mode={mode} counts={counts_msg} weights={weights_msg}",
            flush=True,
        )
    return torch.tensor(weights, dtype=torch.float32)


def _train_celltype_counts(train_ds, n_celltypes: int) -> np.ndarray:
    celltype_ids = getattr(train_ds, "celltype_ids", None)
    row_idx = getattr(train_ds, "row_idx", None)
    if celltype_ids is None or row_idx is None:
        raise AttributeError("train_ds must expose celltype_ids and row_idx for celltype counts.")
    ids = np.asarray(celltype_ids)[np.asarray(row_idx, dtype=np.int64)]
    return np.bincount(ids.astype(np.int64), minlength=int(n_celltypes)).astype(np.float64)


def _celltype_class_weights_from_config(
    train_ds, n_celltypes: int, adv_cfg: Dict[str, Any]
) -> Optional[torch.Tensor]:
    """Class weights for the reference-only celltype adversary, mirroring the
    tech-adversary weighting logic. mode='none' (default) ⇒ unweighted CE.

    Note: counts here are over ALL train cells, not only reference cells. They
    are a rough imbalance prior for the classifier; reference-cell balance is
    additionally reflected by the per-class balanced-accuracy diagnostic."""
    mode = str(adv_cfg.get("class_weight_mode", "none")).lower()
    if mode in {"none", "off", "false", "0"}:
        return None
    if mode not in {"inverse", "sqrt_inverse"}:
        raise ValueError(
            "reference_celltype_adversary.class_weight_mode must be 'none', 'inverse', or 'sqrt_inverse'."
        )

    counts = _train_celltype_counts(train_ds, n_celltypes)
    eps = float(adv_cfg.get("class_weight_eps", 1.0))
    if eps <= 0:
        raise ValueError("class_weight_eps must be positive.")
    smoothed = np.maximum(counts, eps)
    weights = smoothed.sum() / (float(n_celltypes) * smoothed)
    if mode == "sqrt_inverse":
        weights = np.sqrt(weights)

    cap = adv_cfg.get("class_weight_cap", None)
    if cap is not None:
        cap_f = float(cap)
        if cap_f <= 0:
            raise ValueError("class_weight_cap must be positive when provided.")
        weights = np.minimum(weights, cap_f)

    if bool(adv_cfg.get("normalize_class_weights", True)):
        mean_weight = float(weights.mean())
        if mean_weight > 0:
            weights = weights / mean_weight

    weights = weights.astype(np.float32)
    if int(os.environ.get("RANK", "0")) == 0:
        print(
            f"[kmlee-adversary] reference-celltype class weighting mode={mode} "
            f"n_celltypes={n_celltypes}",
            flush=True,
        )
    return torch.tensor(weights, dtype=torch.float32)


def _build_aux_celltype_adv_config(ct_adv_cfg: Dict[str, Any]) -> AuxReferenceAdversaryConfig:
    """v15 auxiliary reference batch config (v3.reference_celltype_adversary.aux_batch).

    OFF by default: returned config has enabled=False unless BOTH the celltype
    adversary is enabled AND its nested ``aux_batch.enabled`` is true (and
    aux_batch_size>0). When disabled the V3Trainer never builds the aux runner,
    so the training loop is byte-identical to the v14 path.
    """
    adv_on = bool(ct_adv_cfg.get("enabled", False))
    aux_cfg = dict(ct_adv_cfg.get("aux_batch", {}))
    enabled = adv_on and bool(aux_cfg.get("enabled", False))
    return AuxReferenceAdversaryConfig(
        enabled=enabled,
        aux_batch_size=int(aux_cfg.get("aux_batch_size", 0)),
        aux_every=int(aux_cfg.get("aux_every", 1)),
        classifier_pretrain_epochs=int(aux_cfg.get("classifier_pretrain_epochs", 0)),
        grl_pretrain_strength=float(aux_cfg.get("grl_pretrain_strength", 0.0)),
        grl_ramp_epochs=int(aux_cfg.get("grl_ramp_epochs", 3)),
        celltype_balanced=bool(aux_cfg.get("celltype_balanced", True)),
        num_workers=int(aux_cfg.get("num_workers", 0)),
        seed=int(aux_cfg.get("seed", 42)),
        min_cells_per_celltype=int(aux_cfg.get("min_cells_per_celltype", 1)),
    )


def _freeze_reference_prior_mu_if_requested(system) -> None:
    ref_cfg = _reference_prior_section()
    if bool(ref_cfg.get("enabled", False)) and bool(ref_cfg.get("freeze_mu", True)):
        system.prior.mu_embedding.weight.requires_grad_(False)
        if int(os.environ.get("RANK", "0")) == 0:
            print(
                "[kmlee-prior] prior.mu_embedding frozen; it will be refreshed from reference cells.",
                flush=True,
            )


_PRECISION_CONTEXT_BUFFER_KEYS = frozenset(
    {
        "precision_head.source_module",
        "precision_head.source_latent",
        "precision_head.source_observed",
        "precision_head.source_reliability",
        "precision_head.context_celltype",
        "precision_head.context_region",
    }
)
_PRECISION_INTERACTION_BUFFER_PREFIX = "precision_head.interaction_basis."
_PRECISION_MODULE_LOCAL_CONTEXT_BUFFER_PREFIXES = (
    "precision_head.module_local_reliability",
    "precision_head.module_local_input_clip",
    "precision_head.module_local_output_cap",
)


def _precision_warm_start_missing_key_allowed(
    key: str,
    *,
    protected_context_keys: set[str],
    strict_precision_upgrade: bool,
) -> bool:
    """Audit an absent key when upgrading an existing PRISM checkpoint.

    Once both source and target contain a precision head, a module-local
    experiment may introduce only its dedicated namespace.  This prevents an
    unrelated accidental architecture change from being silently initialized
    during the same warm start.  Legacy non-upgrade paths retain the historical
    permissive precision-head behavior.
    """

    if key in protected_context_keys:
        return True
    if strict_precision_upgrade:
        return key.startswith("precision_head.module_local_")
    return key.startswith("precision_head.")


def _filter_w44_context_buffers_for_all64_retrain(
    system,
    state,
    *,
    warm_start: dict[str, Any],
    fixed_result_path: Any,
    exact_job: Dict[str, Any] | None = None,
) -> tuple[Any, dict[str, torch.Tensor]]:
    """Keep freshly rebuilt context buffers during a W44 warm start.

    The model is constructed before this function runs.  Therefore its PRISM
    buffers already contain the context residualizer and optional pathology
    interaction basis fitted on the current training dataset.  Loading a W44
    checkpoint naively would silently overwrite those persistent buffers.
    """

    exact_rebuild = bool(
        exact_job is not None
        and exact_job.get("donor_partition") == "confirmation"
    )
    legacy_rebuild = bool(
        warm_start.get("rebuild_precision_context_buffers", False)
    )
    if not exact_rebuild and not legacy_rebuild:
        return state, {}
    if legacy_rebuild and not fixed_result_path:
        raise ValueError(
            "rebuild_precision_context_buffers is restricted to final "
            "fixed-generator retraining"
        )
    if not bool(_section("precision_medicine").get("enabled", False)):
        raise ValueError(
            "rebuild_precision_context_buffers requires precision_medicine.enabled"
        )
    configured_context = _section("precision_medicine").get("context_npz_path")
    if exact_rebuild:
        expected_fit_partition = "W_PLUS_R54"
        expected_fit_count = 54
        declared_context = exact_job["source_paths"]["donor_context"]
        expected_context_sha = exact_job["source_sha256"]["donor_context"]
    else:
        if warm_start.get("rebuilt_precision_context_fit_partition") != "train64":
            raise ValueError(
                "rebuilt_precision_context_fit_partition must be 'train64'"
            )
        expected_fit_partition = "train64"
        expected_fit_count = 64
        declared_context = warm_start.get("rebuilt_precision_context_path")
        expected_context_sha = _validated_sha256(
            warm_start.get("rebuilt_precision_context_expected_sha256"),
            label="warm_start.rebuilt_precision_context_expected_sha256",
        )
    if not configured_context or not declared_context:
        raise ValueError(
            "rebuilt precision context requires configured and warm-start paths"
        )
    if Path(configured_context).resolve() != Path(declared_context).resolve():
        raise ValueError(
            "rebuilt_precision_context_path differs from "
            "precision_medicine.context_npz_path"
        )
    observed_context_sha = _sha256_file(Path(declared_context).resolve())
    if observed_context_sha != expected_context_sha:
        raise RuntimeError(
            "rebuilt precision context SHA256 mismatch: "
            f"{observed_context_sha} != {expected_context_sha}"
        )

    provenance = getattr(system, "precision_context_fit_provenance", None)
    if provenance is None:
        provenance = getattr(
            getattr(system, "precision_head", None),
            "context_fit_provenance",
            None,
        )
    if not isinstance(provenance, dict):
        raise RuntimeError(
            "system lacks runtime-observed precision context fit provenance"
        )
    required_provenance = {
        "fit_partition": expected_fit_partition,
        "observed_fit_donor_count": expected_fit_count,
        "table_train_donor_count": 64,
        "row_level_dataset_view_observed": True,
        "context_npz_path": str(Path(declared_context).resolve()),
    }
    for name, expected in required_provenance.items():
        if provenance.get(name) != expected:
            raise RuntimeError(
                "fresh precision context did not prove the requested fitting set: "
                f"{name}={provenance.get(name)!r}, expected {expected!r}"
            )
    fit_names = provenance.get("observed_fit_donor_names")
    if (
        not isinstance(fit_names, list)
        or len(fit_names) != expected_fit_count
        or len({str(value) for value in fit_names}) != expected_fit_count
    ):
        raise RuntimeError(
            "fresh precision context reports the wrong unique fit donors"
        )

    target_state = system.state_dict()
    missing_target = sorted(_PRECISION_CONTEXT_BUFFER_KEYS - set(target_state))
    if missing_target:
        raise RuntimeError(
            "fresh PRISM system lacks required context buffers: "
            f"{missing_target}"
        )
    checkpoint_keys = set(state)
    missing_checkpoint = sorted(
        _PRECISION_CONTEXT_BUFFER_KEYS - checkpoint_keys
    )
    if missing_checkpoint:
        raise RuntimeError(
            "W44 checkpoint lacks the context buffers that must be replaced: "
            f"{missing_checkpoint}"
        )
    skip = set(_PRECISION_CONTEXT_BUFFER_KEYS)
    skip.update(
        key
        for key in target_state
        if key.startswith(_PRECISION_INTERACTION_BUFFER_PREFIX)
    )
    skip.update(
        key
        for key in target_state
        if key.startswith(_PRECISION_MODULE_LOCAL_CONTEXT_BUFFER_PREFIXES)
    )
    filtered = state.copy()
    skipped_snapshots: dict[str, torch.Tensor] = {}
    for key in sorted(skip):
        if key in filtered:
            skipped_snapshots[key] = target_state[key].detach().clone()
            filtered.pop(key)
    if set(_PRECISION_CONTEXT_BUFFER_KEYS) - set(skipped_snapshots):
        raise RuntimeError("failed to protect every required PRISM context buffer")
    system.rebuilt_precision_context_provenance = {
        **provenance,
        "context_npz_sha256": observed_context_sha,
        "skipped_source_checkpoint_buffer_keys": sorted(skipped_snapshots),
    }
    return filtered, skipped_snapshots


def _commit_generator_exact_training_mask(
    system,
    *,
    job: Dict[str, Any],
    observed_checkpoint_sha256: str,
) -> dict[str, Any]:
    """Commit a sealed evaluation-job mask before optimizer construction."""

    if job["source_sha256"]["checkpoint"] != observed_checkpoint_sha256:
        raise RuntimeError("exact training source checkpoint SHA256 mismatch")
    decoder = getattr(system, "decoder", None)
    if decoder is None or not hasattr(decoder, "set_generator_active_mask"):
        raise RuntimeError("exact training requires a hard-mask decoder")
    candidate_count = int(job["candidate_count"])
    if int(getattr(decoder, "n_generators", 0)) != candidate_count:
        raise ValueError("exact training candidate count differs from decoder")
    selected = tuple(int(value) for value in job["selected_local_generator_ids"])
    expected_mask_sha = _generator_mask_sha256(
        candidate_count=candidate_count,
        selected_local_generator_ids=selected,
    )
    if job["mask_sha256"] != expected_mask_sha:
        raise RuntimeError("exact training mask SHA256 mismatch")
    configured_registry = _section("module_tokenizer").get("registry_json_path")
    if not configured_registry:
        raise ValueError("exact training requires module registry path")
    if Path(configured_registry).resolve() != Path(
        job["source_paths"]["registry"]
    ).resolve() or _sha256_file(Path(configured_registry)) != job["source_sha256"][
        "registry"
    ]:
        raise RuntimeError("exact training registry path/SHA256 mismatch")
    allowlist = getattr(system, "generator_exact_training_allowlist", None)
    if not isinstance(allowlist, dict):
        raise RuntimeError("exact training system lacks row-allowlist evidence")
    if (
        allowlist.get("exact_training_job_sha256") != job["job_sha256"]
        or allowlist.get("row_level_donor_allowlist_enforced") is not True
    ):
        raise RuntimeError("exact training row-allowlist evidence differs from job")

    active_mask = torch.zeros(candidate_count, dtype=torch.float32)
    active_mask[list(selected)] = 1.0
    decoder.set_generator_active_mask(active_mask)
    if int(getattr(decoder, "n_active_generators", -1)) != len(selected):
        raise RuntimeError("exact training decoder active count mismatch")
    provenance: dict[str, Any] = {
        "schema_version": "kmlee_bam.generator_exact_fixed_mask_checkpoint.v1",
        "exact_training_job_sha256": job["job_sha256"],
        "donor_partition": job["donor_partition"],
        "candidate_count": candidate_count,
        "hard_active": len(selected),
        "selected_local_generator_ids": list(selected),
        "mask_sha256": job["mask_sha256"],
        "source_checkpoint_sha256": observed_checkpoint_sha256,
        "registry_sha256": job["source_sha256"]["registry"],
        "training_donor_names_sha256": job["training_donor_names_sha256"],
        "row_level_donor_allowlist_enforced": True,
        "fresh_optimizer": True,
        "fixed_mask_for_entire_run": True,
        "mask_trainable": False,
        "same_seed_schedule_across_k": True,
        "official_validation_used": False,
        "official_test_used": False,
        "allowlist_evidence_sha256": allowlist["evidence_sha256"],
    }
    provenance["provenance_sha256"] = _canonical_json_sha256(provenance)
    system.fixed_generator_result_provenance = provenance
    decoder.fixed_generator_result_provenance = provenance
    _install_fixed_generator_integrity_hook(system)
    _assert_fixed_generator_integrity(
        system, stage="generator_exact_warm_start_commit", distributed=False
    )
    return provenance


def _expand_precision_personal_rank_for_warm_start(
    system,
    state: dict[str, torch.Tensor],
    *,
    warm_start: dict[str, Any],
) -> dict[str, torch.Tensor]:
    """Expand an audited PRISM personal rank while preserving prior axes.

    ``load_state_dict(strict=False)`` still rejects tensors whose keys match but
    whose shapes differ.  A personal-rank experiment intentionally changes six
    such tensors.  When explicitly requested in ``warm_start``, copy every
    source-rank slice bit-for-bit into the freshly initialized target tensor and
    leave only the appended dimensions at their deterministic target
    initialization.  Any undeclared shape or rank mismatch fails closed.
    """

    contract = warm_start.get("personal_rank_expansion")
    if contract is None:
        return state
    if not isinstance(contract, dict):
        raise ValueError("warm_start.personal_rank_expansion must be an object")
    allowed_fields = {"source_rank", "target_rank"}
    unexpected_fields = sorted(set(contract) - allowed_fields)
    if unexpected_fields:
        raise ValueError(
            "warm_start.personal_rank_expansion has unsupported fields: "
            f"{unexpected_fields}"
        )
    source_rank = int(contract.get("source_rank", 0))
    target_rank = int(contract.get("target_rank", 0))
    if source_rank <= 0 or target_rank <= source_rank:
        raise ValueError(
            "warm_start.personal_rank_expansion requires "
            "0 < source_rank < target_rank"
        )
    precision_head = getattr(system, "precision_head", None)
    observed_target_rank = int(getattr(precision_head, "personal_rank", 0))
    if observed_target_rank != target_rank:
        raise RuntimeError(
            "personal-rank expansion target does not match the constructed "
            f"precision head: {target_rank} != {observed_target_rank}"
        )

    target_state = system.state_dict()
    expansion_axes = {
        "precision_head.personal_basis": 1,
        "precision_head.response_basis": 2,
        "precision_head.personal_posterior.4.weight": 0,
        "precision_head.personal_posterior.4.bias": 0,
        "precision_head.pathology_adversary.weight": 1,
        "precision_head.age_adversary.weight": 1,
    }
    expanded = dict(state)
    audit: dict[str, dict[str, Any]] = {}
    for key, axis in expansion_axes.items():
        if key not in state or key not in target_state:
            raise RuntimeError(
                f"personal-rank expansion requires checkpoint and target key {key}"
            )
        source = state[key]
        target = target_state[key]
        if source.ndim != target.ndim:
            raise RuntimeError(
                f"personal-rank expansion ndim mismatch for {key}: "
                f"{source.ndim} != {target.ndim}"
            )
        expected_source_shape = list(target.shape)
        expected_source_shape[axis] = source_rank
        if list(source.shape) != expected_source_shape:
            raise RuntimeError(
                f"personal-rank expansion source shape mismatch for {key}: "
                f"{tuple(source.shape)} != {tuple(expected_source_shape)}"
            )
        if int(target.shape[axis]) != target_rank:
            raise RuntimeError(
                f"personal-rank expansion target shape mismatch for {key}: "
                f"axis {axis} has {target.shape[axis]} != {target_rank}"
            )
        migrated = target.detach().clone()
        source_slice = [slice(None)] * source.ndim
        source_slice[axis] = slice(0, source_rank)
        migrated[tuple(source_slice)].copy_(source.to(migrated.dtype))
        if not torch.equal(
            migrated[tuple(source_slice)].detach().cpu(), source.detach().cpu()
        ):
            raise RuntimeError(
                f"personal-rank expansion failed to preserve source slice for {key}"
            )
        expanded[key] = migrated
        audit[key] = {
            "axis": int(axis),
            "source_shape": list(source.shape),
            "target_shape": list(target.shape),
        }

    remaining_shape_mismatches = sorted(
        key
        for key, value in expanded.items()
        if key in target_state and tuple(value.shape) != tuple(target_state[key].shape)
    )
    if remaining_shape_mismatches:
        raise RuntimeError(
            "personal-rank expansion left undeclared checkpoint shape mismatches: "
            f"{remaining_shape_mismatches}"
        )
    system.personal_rank_expansion_provenance = {
        "schema_version": "kmlee_bam.personal_rank_expansion.v1",
        "source_rank": source_rank,
        "target_rank": target_rank,
        "migrated_tensors": audit,
    }
    if int(os.environ.get("RANK", "0")) == 0:
        print(
            "[warm-start audit] expanded PRISM personal rank "
            f"{source_rank}->{target_rank}; preserved prior axes in "
            f"{len(expansion_axes)} tensors",
            flush=True,
        )
    return expanded


def _warm_start_from_checkpoint(system) -> None:
    """Warm-start: initialise model weights from a prior run's checkpoint
    (top-level config ``warm_start.init_weights_path``) before training. Loads
    ``system_state_dict`` ONLY (not optimizer/epoch/v8 state), so a fresh run with
    a new loss (e.g. the Phase-2 variance floor) continues from a trained base on a
    clean schedule. Architecture must match. No-op when the key is absent.
    A top-level section (like ``_doc_pointer``) is ignored by ``load_config`` but
    still readable via ``_section`` from the raw settings."""
    warm_start = _section("warm_start")
    path = warm_start.get("init_weights_path")
    fixed_result_path = warm_start.get("fixed_generator_result_path")
    exact_job = _load_generator_exact_training_job()
    if not path:
        if fixed_result_path:
            raise ValueError(
                "warm_start.fixed_generator_result_path requires "
                "warm_start.init_weights_path"
            )
        return
    if fixed_result_path:
        incompatible_searches = [
            section
            for section in ("learned_generator_count", "generator_budget_search")
            if bool(_section(section).get("enabled", False))
        ]
        if incompatible_searches:
            raise ValueError(
                "warm_start.fixed_generator_result_path cannot be combined with "
                "an enabled generator search: "
                f"{incompatible_searches}"
            )
    if exact_job is not None and fixed_result_path:
        raise ValueError(
            "generator_exact_training cannot consume a canonical fixed result"
        )
    expected_sha256 = warm_start.get("expected_sha256")
    observed_sha256 = None
    if expected_sha256 is not None or fixed_result_path:
        observed_sha256 = _sha256_file(Path(path))
    if expected_sha256 is not None:
        expected_sha256 = _validated_sha256(
            expected_sha256,
            label="warm_start.expected_sha256",
        )
        if observed_sha256 != expected_sha256:
            raise RuntimeError(
                "warm-start checkpoint SHA256 mismatch: "
                f"{observed_sha256} != {expected_sha256}"
            )
    rank0 = int(os.environ.get("RANK", "0")) == 0
    payload = torch.load(path, map_location="cpu", weights_only=False)
    state = (
        payload["system_state_dict"]
        if isinstance(payload, dict) and "system_state_dict" in payload
        else payload
    )
    checkpoint_has_precision = any(
        str(key).startswith("precision_head.") for key in state
    )
    if bool(warm_start.get("require_precision_head", False)) and not any(
        str(key).startswith("precision_head.") for key in state
    ):
        raise RuntimeError(
            "warm-start checkpoint is required to contain precision_head keys"
        )
    state, protected_context_buffers = (
        _filter_w44_context_buffers_for_all64_retrain(
            system,
            state,
            warm_start=warm_start,
            fixed_result_path=fixed_result_path,
            exact_job=exact_job,
        )
    )
    state = _expand_precision_personal_rank_for_warm_start(
        system,
        state,
        warm_start=warm_start,
    )
    missing, unexpected = system.load_state_dict(state, strict=False)
    if protected_context_buffers:
        current_state = system.state_dict()
        changed = [
            key
            for key, expected in protected_context_buffers.items()
            if key not in current_state
            or not torch.equal(current_state[key].detach().cpu(), expected.cpu())
        ]
        if changed:
            raise RuntimeError(
                "fresh train64 precision context buffers changed during W44 "
                f"warm start: {changed}"
            )
    prism_enabled = bool(_section("precision_medicine").get("enabled", False))
    if prism_enabled:
        strict_precision_upgrade = bool(
            checkpoint_has_precision
            and _section("precision_medicine").get(
                "module_local_enabled", False
            )
        )
        allowed_missing_exact = {
            # Hard, non-trainable architecture state introduced after the
            # canonical CT64/all414 checkpoints.  It is intentionally
            # initialised to all-active before the search controller starts.
            "decoder.generator_active_mask",
        }
        agp_warm_start = str(_section("encoder").get("pooling", "cls")) == "agp"
        bad_missing = [
            key
            for key in missing
            if not _precision_warm_start_missing_key_allowed(
                key,
                protected_context_keys=set(protected_context_buffers),
                strict_precision_upgrade=strict_precision_upgrade,
            )
            and key not in allowed_missing_exact
            and not key.startswith("generator_count_gate.")
            and not key.startswith("generator_count_objective.")
            and key != "generator_count_epoch_state"
            and not (
                agp_warm_start
                and key.startswith("state_encoder.attention_pool.")
            )
        ]
        allowed_unexpected_prefixes = (
            "pathology_aux_head.",
            "decoder.path_single_",
            "decoder.path_V",
            "decoder.pathology_head.",
        )
        bad_unexpected = [
            key
            for key in unexpected
            if not key.startswith(allowed_unexpected_prefixes)
        ]
        if bad_missing or bad_unexpected:
            raise RuntimeError(
                "PRISM warm-start mismatch outside the audited new/removed paths: "
                f"bad_missing={bad_missing} bad_unexpected={bad_unexpected}"
            )

    # A learned-generator search is discovery-only until its exact binary mask
    # is committed to the warm-started decoder.  Apply it only after loading the
    # checkpoint, because the checkpoint's persistent all-on mask would
    # otherwise overwrite the selected architecture.  Every selection contract
    # is checked before mutating the decoder.
    fixed_result_provenance = None
    if exact_job is not None:
        if observed_sha256 is None:
            observed_sha256 = _sha256_file(Path(path))
        fixed_result_provenance = _commit_generator_exact_training_mask(
            system,
            job=exact_job,
            observed_checkpoint_sha256=observed_sha256,
        )
    if fixed_result_path:
        if observed_sha256 is None:
            raise RuntimeError("internal error: checkpoint SHA256 was not computed")
        result_path = Path(fixed_result_path).resolve()
        result_sha256 = _sha256_file(result_path)
        expected_result_sha256 = warm_start.get(
            "fixed_generator_result_expected_sha256"
        )
        if expected_result_sha256 is not None:
            expected_result_sha256 = _validated_sha256(
                expected_result_sha256,
                label="warm_start.fixed_generator_result_expected_sha256",
            )
            if result_sha256 != expected_result_sha256:
                raise RuntimeError(
                    "fixed-generator result SHA256 mismatch: "
                    f"{result_sha256} != {expected_result_sha256}"
                )
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"failed to read fixed-generator result: {result_path}"
            ) from exc
        if not isinstance(result, dict):
            raise ValueError("fixed-generator result must be a JSON object")
        if result.get("schema_version") != "kmlee_bam.learned_generator_result.v1":
            raise ValueError("fixed-generator result schema_version mismatch")
        canonical_payload_sha256 = _validated_sha256(
            result.get("canonical_payload_sha256"),
            label="fixed-generator canonical_payload_sha256",
        )
        observed_canonical_payload_sha256 = _canonical_json_sha256(
            {
                key: value
                for key, value in result.items()
                if key != "canonical_payload_sha256"
            }
        )
        if canonical_payload_sha256 != observed_canonical_payload_sha256:
            raise RuntimeError(
                "fixed-generator canonical payload SHA256 mismatch: "
                f"{canonical_payload_sha256} != "
                f"{observed_canonical_payload_sha256}"
            )
        if result.get("selection_stage") != "exact_confirmation":
            raise ValueError(
                "fixed-generator result requires "
                "selection_stage='exact_confirmation'"
            )
        for field in (
            "exact_confirmation_passed",
            "stable",
            "canonical_mask_allowed",
            "every_confirmation_constraint_passed",
        ):
            if result.get(field) is not True:
                raise ValueError(f"fixed-generator result requires {field}=true")
        for field in ("official_validation_used", "official_test_used"):
            if result.get(field) is not False:
                raise ValueError(f"fixed-generator result requires {field}=false")

        decoder = getattr(system, "decoder", None)
        if decoder is None or not hasattr(decoder, "set_generator_active_mask"):
            raise RuntimeError(
                "fixed-generator result requires a decoder with "
                "set_generator_active_mask"
            )
        n_generators = int(getattr(decoder, "n_generators", 0))
        if n_generators <= 0:
            raise RuntimeError(
                "fixed-generator result requires a positive decoder generator count"
            )
        candidate_count = result.get("candidate_count")
        if isinstance(candidate_count, bool) or not isinstance(candidate_count, int):
            raise ValueError("fixed-generator candidate_count must be an integer")
        if int(candidate_count) != n_generators:
            raise ValueError(
                "fixed-generator candidate_count does not match the decoder "
                f"generator count ({candidate_count} != {n_generators})"
            )
        selected_raw = result.get("selected_local_generator_ids")
        if not isinstance(selected_raw, list) or not selected_raw:
            raise ValueError(
                "fixed-generator selected_local_generator_ids must be a "
                "non-empty list"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in selected_raw
        ):
            raise ValueError(
                "fixed-generator selected_local_generator_ids must contain integers"
            )
        selected_ids = tuple(int(value) for value in selected_raw)
        if len(set(selected_ids)) != len(selected_ids):
            raise ValueError(
                "fixed-generator selected_local_generator_ids must be unique"
            )
        if any(value < 0 or value >= n_generators for value in selected_ids):
            raise ValueError(
                "fixed-generator selected_local_generator_ids contain an "
                "out-of-range id"
            )
        hard_active = result.get("hard_active")
        if isinstance(hard_active, bool) or not isinstance(hard_active, int):
            raise ValueError("fixed-generator hard_active must be an integer")
        if int(hard_active) != len(selected_ids):
            raise ValueError(
                "fixed-generator hard_active does not match the selected id count"
            )

        source_checkpoint_sha256 = _validated_sha256(
            result.get("source_checkpoint_sha256"),
            label="fixed-generator source_checkpoint_sha256",
        )
        if source_checkpoint_sha256 != observed_sha256:
            raise RuntimeError(
                "fixed-generator source checkpoint SHA256 mismatch: "
                f"{source_checkpoint_sha256} != {observed_sha256}"
            )
        registry_sha256 = _validated_sha256(
            result.get("registry_sha256"),
            label="fixed-generator registry_sha256",
        )
        ranking_result_sha256 = _validated_sha256(
            result.get("ranking_result_sha256"),
            label="fixed-generator ranking_result_sha256",
        )
        source_hashes = result.get("source_hashes")
        if source_hashes is not None:
            if not isinstance(source_hashes, dict):
                raise ValueError("fixed-generator source_hashes must be an object")
            top_level_hashes = {
                "source_checkpoint_sha256": source_checkpoint_sha256,
                "registry_sha256": registry_sha256,
                "ranking_result_sha256": ranking_result_sha256,
            }
            for field, top_level_digest in top_level_hashes.items():
                nested_digest = _validated_sha256(
                    source_hashes.get(field),
                    label=f"fixed-generator source_hashes.{field}",
                )
                if nested_digest != top_level_digest:
                    raise RuntimeError(
                        "fixed-generator source_hashes disagrees with the "
                        f"top-level {field}"
                    )
        mask_sha256 = _validated_sha256(
            result.get("mask_sha256"),
            label="fixed-generator mask_sha256",
        )
        expected_mask_sha256 = _generator_mask_sha256(
            candidate_count=n_generators,
            selected_local_generator_ids=selected_ids,
        )
        if mask_sha256 != expected_mask_sha256:
            raise RuntimeError(
                "fixed-generator binary mask SHA256 mismatch: "
                f"{mask_sha256} != {expected_mask_sha256}"
            )

        metadata = getattr(decoder, "decoder_architecture_metadata", None)
        metadata_registry_path = (
            metadata.get("registry_path")
            if isinstance(metadata, dict)
            else None
        )
        configured_registry_path = _section("module_tokenizer").get(
            "registry_json_path"
        )
        if metadata_registry_path and configured_registry_path:
            if Path(metadata_registry_path).resolve() != Path(
                configured_registry_path
            ).resolve():
                raise RuntimeError(
                    "decoder and config disagree on the current registry path"
                )
        current_registry_path = metadata_registry_path or configured_registry_path
        if not current_registry_path:
            raise RuntimeError(
                "fixed-generator result requires the current registry path in "
                "decoder architecture metadata or module_tokenizer.registry_json_path"
            )
        try:
            observed_registry_sha256 = _sha256_file(Path(current_registry_path))
        except OSError as exc:
            raise RuntimeError(
                "failed to hash the current generator registry file: "
                f"{current_registry_path}"
            ) from exc
        if registry_sha256 != observed_registry_sha256:
            raise RuntimeError(
                "fixed-generator registry SHA256 mismatch: "
                f"{registry_sha256} != {observed_registry_sha256}"
            )

        selected_registry_indices_raw = result.get("selected_registry_indices")
        if not isinstance(selected_registry_indices_raw, list):
            raise ValueError(
                "fixed-generator selected_registry_indices must be a list"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in selected_registry_indices_raw
        ):
            raise ValueError(
                "fixed-generator selected_registry_indices must contain integers"
            )
        selected_registry_indices = tuple(
            int(value) for value in selected_registry_indices_raw
        )
        if len(selected_registry_indices) != len(selected_ids):
            raise ValueError(
                "fixed-generator selected_registry_indices count does not match "
                "hard_active"
            )
        if len(set(selected_registry_indices)) != len(selected_registry_indices):
            raise ValueError(
                "fixed-generator selected_registry_indices must be unique"
            )
        if any(value < 0 for value in selected_registry_indices):
            raise ValueError(
                "fixed-generator selected_registry_indices must be non-negative"
            )
        selected_registry_names_raw = result.get("selected_registry_module_names")
        if not isinstance(selected_registry_names_raw, list) or any(
            not isinstance(value, str) for value in selected_registry_names_raw
        ):
            raise ValueError(
                "fixed-generator selected_registry_module_names must be a list "
                "of strings"
            )
        selected_registry_names = tuple(selected_registry_names_raw)
        if len(selected_registry_names) != len(selected_ids):
            raise ValueError(
                "fixed-generator selected_registry_module_names count does not "
                "match hard_active"
            )

        decoder_registry_indices = getattr(
            decoder,
            "selected_registry_module_indices",
            None,
        )
        if decoder_registry_indices is not None and len(decoder_registry_indices) > 0:
            if len(decoder_registry_indices) != n_generators:
                raise RuntimeError(
                    "decoder selected_registry_module_indices length does not "
                    "match candidate_count"
                )
            expected_registry_indices = tuple(
                int(decoder_registry_indices[local_id]) for local_id in selected_ids
            )
            if selected_registry_indices != expected_registry_indices:
                raise RuntimeError(
                    "fixed-generator selected_registry_indices do not match the "
                    "decoder candidate mapping"
                )
        decoder_registry_names = getattr(
            decoder,
            "selected_registry_module_names",
            None,
        )
        if decoder_registry_names is not None and len(decoder_registry_names) > 0:
            if len(decoder_registry_names) != n_generators:
                raise RuntimeError(
                    "decoder selected_registry_module_names length does not match "
                    "candidate_count"
                )
            expected_registry_names = tuple(
                str(decoder_registry_names[local_id]) for local_id in selected_ids
            )
            if selected_registry_names != expected_registry_names:
                raise RuntimeError(
                    "fixed-generator selected_registry_module_names do not match "
                    "the decoder candidate mapping"
                )

        # This is the first mutation of the learned fixed architecture.  Every
        # exact-confirmation and provenance check above has already passed.
        active_mask = torch.zeros(
            n_generators,
            dtype=torch.float32,
        )
        active_mask[list(selected_ids)] = 1.0
        decoder.set_generator_active_mask(active_mask)
        if int(getattr(decoder, "n_active_generators", -1)) != int(hard_active):
            raise RuntimeError(
                "decoder active generator count differs after fixed-mask commit"
            )
        fixed_result_provenance = {
            "schema_version": result["schema_version"],
            "path": str(result_path),
            "sha256": result_sha256,
            "hard_active": int(hard_active),
            "selected_local_generator_ids": list(selected_ids),
            "selected_registry_indices": list(selected_registry_indices),
            "selected_registry_module_names": list(selected_registry_names),
            "candidate_count": int(candidate_count),
            "source_checkpoint_sha256": source_checkpoint_sha256,
            "registry_sha256": registry_sha256,
            "ranking_result_sha256": ranking_result_sha256,
            "mask_sha256": mask_sha256,
            "canonical_payload_sha256": canonical_payload_sha256,
        }
        # Plain metadata deliberately does not enter state_dict.  The committed
        # binary decoder buffer is checkpoint-persistent; this sidecar metadata
        # lets callers/manifests record the exact selection artifact as well.
        system.fixed_generator_result_provenance = fixed_result_provenance
        decoder.fixed_generator_result_provenance = fixed_result_provenance
        _install_fixed_generator_integrity_hook(system)
        _assert_fixed_generator_integrity(
            system,
            stage="warm_start_commit",
            distributed=False,
        )
    if rank0:
        print(
            f"[warm-start] loaded system_state_dict from {path} "
            f"(missing={len(missing)}, unexpected={len(unexpected)})",
            flush=True,
        )
        if expected_sha256 is not None:
            print(
                f"[warm-start audit] SHA256={expected_sha256}",
                flush=True,
            )
        if fixed_result_provenance is not None:
            print(
                "[warm-start fixed-generator] "
                f"active={fixed_result_provenance['hard_active']}/"
                f"{system.decoder.n_generators} "
                f"SHA256={fixed_result_provenance['sha256']} "
                f"result={fixed_result_provenance['path']}",
                flush=True,
            )
        if protected_context_buffers:
            print(
                "[warm-start context rebuild] retained fresh train64 PRISM "
                f"buffers={len(protected_context_buffers)} instead of W44 "
                "checkpoint buffers",
                flush=True,
            )
        if missing or unexpected:
            if prism_enabled:
                print(
                    "[warm-start audit] PASS: every missing key belongs to the new "
                    "precision_head or the audited hard generator mask; every "
                    "unexpected key belongs to the intentionally removed legacy "
                    "pathology route",
                    flush=True,
                )
                print(f"[warm-start audit] allowed_missing={list(missing)}", flush=True)
                print(f"[warm-start audit] allowed_unexpected={list(unexpected)}", flush=True)
            else:
                print(
                    f"[warm-start] WARNING missing[:6]={list(missing)[:6]} "
                    f"unexpected[:6]={list(unexpected)[:6]}",
                    flush=True,
                )


def _warm_start_trainer_state_from_checkpoint(trainer) -> None:
    """Restore warm-start trainer sidecars that are safe for a fresh schedule.

    ``_warm_start_from_checkpoint`` intentionally loads only model weights because
    optimizer/history should restart cleanly.  V7a's difficulty profile is
    different: without restoring it, a warm-started PHU run must rebuild the
    profile and DDP-broadcast it before step 1.  That collective is expensive and
    has proven fragile on long production runs.  Restoring the saved V7a profile
    keeps the fresh optimizer schedule while avoiding redundant bootstrap.
    """
    path = _section("warm_start").get("init_weights_path")
    if not path or not hasattr(trainer, "v7a_load_state_dict"):
        return

    rank0 = int(os.environ.get("RANK", "0")) == 0
    try:
        payload = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:  # noqa: BLE001
        if rank0:
            print(
                f"[warm-start] WARNING failed to read trainer state from {path}: {exc}",
                flush=True,
            )
        return

    v7a_state = payload.get("v7a_state") if isinstance(payload, dict) else None
    if not v7a_state:
        if rank0:
            print("[warm-start] no v7a_state found in checkpoint; PHU profile will bootstrap.")
        return

    try:
        trainer.v7a_load_state_dict(v7a_state)
    except Exception as exc:  # noqa: BLE001
        if rank0:
            print(
                f"[warm-start] WARNING failed to restore v7a_state from {path}: {exc}",
                flush=True,
            )
        return

    if rank0:
        profile_loaded = bool(v7a_state.get("profile"))
        bank_loaded = bool(v7a_state.get("bank"))
        ancova_loaded = bool(v7a_state.get("ancova"))
        print(
            "[warm-start] restored v7a_state "
            f"(profile={profile_loaded}, bank={bank_loaded}, ancova={ancova_loaded}) "
            f"from {path}",
            flush=True,
        )


def build_system(cfg, train_ds):
    _prepare_runtime_state(cfg, train_ds)
    base_system = _BASE_BUILD_SYSTEM(cfg, train_ds)
    adv_cfg = _tech_adv_section()
    tech_adversary = None
    if bool(adv_cfg.get("enabled", False)):
        n_tech = len(train_ds.spec.batch_vocab)
        tech_adversary = ConditionalTechAdversary(
            d_z=cfg.encoder.d_z,
            n_celltypes=len(train_ds.spec.celltype_vocab),
            n_tech=n_tech,
            celltype_embed_dim=int(adv_cfg.get("celltype_embed_dim", 16)),
            hidden_dim=int(adv_cfg.get("hidden_dim", 64)),
            dropout=float(adv_cfg.get("dropout", 0.1)),
            grl_strength=float(adv_cfg.get("grl_strength", 1.0)),
            class_weights=_tech_class_weights_from_config(train_ds, n_tech, adv_cfg),
        )
    ct_adv_cfg = _celltype_adv_section()
    celltype_adversary = None
    if bool(ct_adv_cfg.get("enabled", False)):
        n_celltypes = len(train_ds.spec.celltype_vocab)
        celltype_adversary = ReferenceCelltypeAdversary(
            d_z=cfg.encoder.d_z,
            n_celltypes=n_celltypes,
            hidden_dim=int(ct_adv_cfg.get("hidden_dim", 64)),
            dropout=float(ct_adv_cfg.get("dropout", 0.1)),
            grl_strength=float(ct_adv_cfg.get("grl_strength", 1.0)),
            class_weights=_celltype_class_weights_from_config(train_ds, n_celltypes, ct_adv_cfg),
        )
    sex_adv_cfg = _sex_adv_section()
    sex_adversary = None
    if bool(sex_adv_cfg.get("enabled", False)):
        sex_adversary = SexAdversary(
            d_z=cfg.encoder.d_z,
            hidden_dim=int(sex_adv_cfg.get("hidden_dim", 64)),
            dropout=float(sex_adv_cfg.get("dropout", 0.1)),
            grl_strength=float(sex_adv_cfg.get("grl_strength", 1.0)),
        )
    system = V3OrdinalBAMSystem(
        base_system=base_system,
        tech_adversary=tech_adversary,
        celltype_adversary=celltype_adversary,
        sex_adversary=sex_adversary,
    )
    source_provenance = getattr(
        train_ds, "canonical_w44_source_provenance", None
    )
    if source_provenance is not None:
        system.canonical_w44_source_provenance = dict(source_provenance)
    exact_allowlist = getattr(
        train_ds, "generator_exact_training_allowlist", None
    )
    if exact_allowlist is not None:
        if not isinstance(exact_allowlist, dict):
            raise RuntimeError("generator exact allowlist evidence must be an object")
        system.generator_exact_training_allowlist = dict(exact_allowlist)
    _freeze_reference_prior_mu_if_requested(system)
    _warm_start_from_checkpoint(system)
    return system


def build_criterion(cfg):
    criterion = _BASE_BUILD_CRITERION(cfg)
    ref_cfg = _ref_section()
    if bool(ref_cfg.get("enabled", False)):
        criterion.lambda_ref_center = 0.0
    return criterion


def build_optimizer(cfg, system):
    integrity_hook = getattr(system, "_fixed_generator_integrity_hook", None)
    if integrity_hook is not None:
        integrity_hook(stage="optimizer_construction", distributed=False)
    named = [(name, p) for name, p in system.named_parameters() if p.requires_grad]
    if not named:
        raise RuntimeError("No trainable parameters found for optimizer.")
    joint_count_enabled = bool(
        getattr(cfg.learned_generator_count, "enabled", False)
        and str(getattr(cfg.learned_generator_count, "mode", "gate_only"))
        == "joint"
    )
    gate_params = [
        p for name, p in named if name.startswith("generator_count_gate.")
    ]
    pathology_rank_cfg = getattr(cfg, "learned_pathology_rank", None)
    pathology_rank_enabled = bool(
        getattr(pathology_rank_cfg, "enabled", False)
    )
    pathology_rank_params = [
        p
        for name, p in named
        if name.startswith("decoder.pathology_rank_gate.")
    ]
    if joint_count_enabled and len(gate_params) != 1:
        raise RuntimeError(
            "joint generator-count optimizer requires exactly one gate-logit "
            f"parameter, found {len(gate_params)}"
        )
    if not joint_count_enabled and gate_params:
        raise RuntimeError("generator-count gate exists while joint mode is disabled")
    if pathology_rank_enabled and len(pathology_rank_params) != 1:
        raise RuntimeError(
            "learned pathology-rank optimizer requires exactly one gate-logit "
            f"parameter, found {len(pathology_rank_params)}"
        )
    if not pathology_rank_enabled and pathology_rank_params:
        raise RuntimeError(
            "pathology-rank gate exists while learned rank is disabled"
        )
    pm_cfg = getattr(cfg, "precision_medicine", None)
    if pm_cfg is not None and bool(getattr(pm_cfg, "enabled", False)):
        module_local_params = [
            p
            for name, p in named
            if name.startswith("precision_head.module_local_")
        ]
        precision_params = [
            p
            for name, p in named
            if name.startswith("precision_head.")
            and not name.startswith("precision_head.module_local_")
        ]
        backbone_params = [
            p
            for name, p in named
            if not name.startswith("precision_head.")
            and not name.startswith("generator_count_gate.")
            and not name.startswith("decoder.pathology_rank_gate.")
        ]
        if not precision_params or not backbone_params:
            raise RuntimeError("PRISM optimiser split failed: missing precision or backbone parameters")
        module_local_enabled = bool(
            getattr(pm_cfg, "module_local_enabled", False)
        )
        if module_local_enabled and not module_local_params:
            raise RuntimeError(
                "module-local PRISM is enabled but its optimizer group is empty"
            )
        if not module_local_enabled and module_local_params:
            raise RuntimeError(
                "module-local PRISM parameters exist while the feature is disabled"
            )
        multiplier = float(getattr(pm_cfg, "lr_multiplier", 1.0))
        param_groups = [
            {"params": backbone_params, "lr": float(cfg.optim.lr), "name": "backbone_unfrozen"},
            {
                "params": precision_params,
                "lr": float(cfg.optim.lr) * multiplier,
                "name": "precision_new",
            },
        ]
        if module_local_enabled:
            module_local_multiplier = float(
                getattr(pm_cfg, "module_local_lr_multiplier", 10.0)
            )
            param_groups.append(
                {
                    "params": module_local_params,
                    "lr": float(cfg.optim.lr) * module_local_multiplier,
                    "name": "precision_module_local",
                }
            )
        if joint_count_enabled:
            param_groups.append(
                {
                    "params": gate_params,
                    "lr": float(cfg.learned_generator_count.learning_rate),
                    "weight_decay": float(
                        cfg.learned_generator_count.weight_decay
                    ),
                    "name": "joint_generator_gate",
                }
            )
        if pathology_rank_enabled:
            param_groups.append(
                {
                    "params": pathology_rank_params,
                    "lr": float(pathology_rank_cfg.learning_rate),
                    "weight_decay": float(
                        pathology_rank_cfg.weight_decay
                    ),
                    "name": "pathology_rank_gate",
                }
            )
        if int(os.environ.get("RANK", "0")) == 0:
            component_prefixes = (
                "gene_embedding",
                "module_tokenizer",
                "state_encoder",
                "prior",
                "decoder",
                "precision_head",
                "generator_count_gate",
            )
            component_counts = {
                prefix: sum(
                    parameter.numel()
                    for name, parameter in system.named_parameters()
                    if name.startswith(prefix + ".") and parameter.requires_grad
                )
                for prefix in component_prefixes
            }
            frozen = [
                (name, int(parameter.numel()))
                for name, parameter in system.named_parameters()
                if not parameter.requires_grad
            ]
            print(
                "[PRISM-E2E optimizer] "
                f"backbone_unfrozen={sum(p.numel() for p in backbone_params):,} "
                f"lr={cfg.optim.lr:.3g}; precision_new={sum(p.numel() for p in precision_params):,} "
                f"lr={float(cfg.optim.lr) * multiplier:.3g}",
                flush=True,
            )
            if module_local_enabled:
                print(
                    "[PRISM module-local optimizer] "
                    f"params={sum(p.numel() for p in module_local_params):,} "
                    f"lr={float(cfg.optim.lr) * module_local_multiplier:.3g} "
                    f"rank={int(pm_cfg.module_local_rank)}",
                    flush=True,
                )
            if joint_count_enabled:
                print(
                    "[joint-generator optimizer] "
                    f"logits={sum(p.numel() for p in gate_params):,} "
                    f"lr={float(cfg.learned_generator_count.learning_rate):.3g} "
                    f"weight_decay={float(cfg.learned_generator_count.weight_decay):.3g} "
                    "selection=online_hard_concrete candidate_K_grid=false",
                    flush=True,
                )
            if pathology_rank_enabled:
                print(
                    "[pathology-rank optimizer] "
                    f"logits={sum(p.numel() for p in pathology_rank_params):,} "
                    f"lr={float(pathology_rank_cfg.learning_rate):.3g} "
                    f"weight_decay={float(pathology_rank_cfg.weight_decay):.3g} "
                    "fixed_rank=false target_rank=none minimum_rank=0",
                    flush=True,
                )
            print(
                f"[PRISM-E2E trainable-by-component] {component_counts}",
                flush=True,
            )
            print(
                f"[PRISM-E2E intentionally-frozen] {frozen}",
                flush=True,
            )
    else:
        if joint_count_enabled or pathology_rank_enabled:
            backbone_params = [
                p
                for name, p in named
                if not name.startswith("generator_count_gate.")
                and not name.startswith("decoder.pathology_rank_gate.")
            ]
            param_groups = [
                {
                    "params": backbone_params,
                    "lr": float(cfg.optim.lr),
                    "name": "backbone_unfrozen",
                }
            ]
            if joint_count_enabled:
                param_groups.append(
                    {
                        "params": gate_params,
                        "lr": float(cfg.learned_generator_count.learning_rate),
                        "weight_decay": float(
                            cfg.learned_generator_count.weight_decay
                        ),
                        "name": "joint_generator_gate",
                    }
                )
            if pathology_rank_enabled:
                param_groups.append(
                    {
                        "params": pathology_rank_params,
                        "lr": float(pathology_rank_cfg.learning_rate),
                        "weight_decay": float(
                            pathology_rank_cfg.weight_decay
                        ),
                        "name": "pathology_rank_gate",
                    }
                )
        else:
            param_groups = [p for _, p in named]
    if param_groups and isinstance(param_groups[0], dict):
        grouped_ids = [
            id(parameter)
            for group in param_groups
            for parameter in group["params"]
        ]
        expected_ids = {id(parameter) for _, parameter in named}
        if len(grouped_ids) != len(set(grouped_ids)):
            raise RuntimeError("optimizer parameter groups overlap")
        if set(grouped_ids) != expected_ids:
            raise RuntimeError(
                "optimizer parameter groups do not cover every trainable parameter"
            )
    return torch.optim.AdamW(
        param_groups,
        lr=cfg.optim.lr,
        weight_decay=cfg.optim.weight_decay,
        betas=cfg.optim.betas,
        eps=cfg.optim.eps,
        foreach=False,
        fused=False,
    )


def make_trainer(*args, **kwargs):
    ref_cfg = _ref_section()
    adv_cfg = _tech_adv_section()
    ct_adv_cfg = _celltype_adv_section()
    sex_adv_cfg = _sex_adv_section()
    ord_cfg = _ordinal_section()
    global_ref_cfg = _global_ref_section()
    v4_unc_cfg = _v4_unc_section()
    ref_prior_cfg = _reference_prior_section()
    state_cfg = _state_usage_section()
    hier_cfg = _hierarchical_ordinal_section()
    unc_cfg = _uncertainty_residual_section()
    mixer_cfg = _score_mixer_section()
    console_style = os.environ.get("KMLEE_CONSOLE_LOG_STYLE", "").lower()
    # The informative PRISM step block already includes aggregated adversary
    # balance/collapse metrics.  Reprinting the raw six-line v3 block every 100
    # micro-batches made DDP logs look stalled and buried the experiment arm.
    # An explicit v3.diag_print_every still overrides this default for a deep
    # adversary investigation.
    v3_diag_default = (
        0
        if console_style in {"prism_experiment", "prism_informative"}
        else 100
    )

    v7a_cfg = _v7a_section()
    v7a_enabled = bool(v7a_cfg.get("enabled", False))

    if v7a_enabled and unc_cfg.get("enabled", False) and not unc_cfg.get("allow_with_v7a", False):
        # v7a replaces v6 uncertainty residual; do not run both at once UNLESS the config explicitly
        # opts in via v6.uncertainty_residual.allow_with_v7a=true. v28 opts in: the rec_alignment
        # FIREWALL (depth-controlled) is wanted as a guardrail ALONGSIDE PHU — deliberate, minimal
        # (only lambda_rec_alignment; others 0) double-supervision of the BAM head. Default stays safe.
        if int(os.environ.get("RANK", "0")) == 0:
            print(
                "[v7a] WARNING: v7a is enabled while v6.uncertainty_residual is also enabled. "
                "Disabling v6.uncertainty_residual to avoid double-supervision of BAM head.",
                flush=True,
            )
        unc_cfg = dict(unc_cfg)
        unc_cfg["enabled"] = False
    elif v7a_enabled and unc_cfg.get("enabled", False) and unc_cfg.get("allow_with_v7a", False):
        if int(os.environ.get("RANK", "0")) == 0:
            print(
                "[v7a] v6.uncertainty_residual.allow_with_v7a=true → running the rec_alignment "
                "firewall ALONGSIDE PHU (intentional). Watch the BAM uncertainty metrics for conflict.",
                flush=True,
            )

    # v8 (zero-origin / capacity) sits on top of v7a. It is enabled whenever any
    # of its three loss-side mechanisms is on. (Per-rank coefficients are a
    # decoder flag — decoder.per_rank_coeff — and need no trainer involvement.)
    pi_tech_cfg = _section("pi_tech")
    thinning_cfg = _section("thinning")
    cvar_cfg = _section("subgroup_cvar")
    zero_rel_cfg = _zero_reliability_section()
    vf_cfg = _variance_floor_section()
    hmr_cfg = _high_margin_rank_section()
    tech_inv_cfg = _section("tech_invariance")
    v8_enabled = (
        bool(pi_tech_cfg.get("enabled", False))
        or bool(thinning_cfg.get("enabled", False))
        or bool(cvar_cfg.get("enabled", False))
        or bool(zero_rel_cfg.get("enabled", False))
        or bool(vf_cfg.get("enabled", False))
        or bool(hmr_cfg.get("enabled", False))
        or bool(tech_inv_cfg.get("enabled", False))
    )

    if v8_enabled:
        trainer_cls = V8Trainer
    elif v7a_enabled:
        trainer_cls = V7aTrainer
    else:
        trainer_cls = V6Trainer

    # Build v7a-specific kwargs only when active.
    v7a_kwargs: Dict[str, Any] = {}
    if v7a_enabled:
        profile_cfg = dict(v7a_cfg.get("difficulty_profile", {}))
        ancova_cfg = dict(v7a_cfg.get("ancova", {}))
        align_cfg = dict(v7a_cfg.get("alignment", {}))
        bank_cfg = dict(v7a_cfg.get("bank", {}))

        v7a_kwargs.update(
            v7a_enabled=True,
            v7a_profile_config=V7aDifficultyProfileConfig(
                min_count=int(profile_cfg.get("min_count", 100)),
                mad_floor=float(profile_cfg.get("mad_floor", 0.01)),
                winsorize_z=float(profile_cfg.get("winsorize_z", 5.0)),
                mad_to_std=float(profile_cfg.get("mad_to_std", 1.4826)),
                gene_chunk_size=int(profile_cfg.get("gene_chunk_size", 500)),
                ema_momentum=float(profile_cfg.get("ema_momentum", 0.95)),
                max_bucket_samples=int(
                    profile_cfg.get("max_bucket_samples", 100_000)
                ),
            ),
            v7a_bank_momentum=float(bank_cfg.get("ema_momentum", 0.95)),
            v7a_ancova_config=V7aANCOVAConfig(
                ridge_alpha=float(ancova_cfg.get("ridge_alpha", 0.1)),
                min_donors_per_celltype=int(
                    ancova_cfg.get("min_donors_per_celltype", 5)
                ),
                pathology_axes=list(
                    ancova_cfg.get(
                        "pathology_axes",
                        ["Braak_stage", "Thal_phase", "CERAD_score"],
                    )
                ),
                standardize_axes=bool(ancova_cfg.get("standardize_axes", True)),
            ),
            v7a_alignment_config=V7aAlignmentConfig(
                enabled=bool(align_cfg.get("enabled", True)),
                # Production-friendly defaults derived from the L2.5 smoke
                # feedback (alignment→encoder feedback loop). Smoke config can
                # still override these. See model_v7a_implementation_plan.md
                # §9 (open decisions) for tuning rationale.
                lambda_alignment=float(align_cfg.get("lambda_alignment", 0.001)),
                detach_target=bool(align_cfg.get("detach_target", True)),
                loss_form=str(align_cfg.get("loss_form", "huber")),
                huber_delta=float(align_cfg.get("huber_delta", 1.0)),
                alignment_start_epoch=int(
                    align_cfg.get("alignment_start_epoch", 1)
                ),
                ramp_epochs=int(align_cfg.get("ramp_epochs", 5)),
                u_total_cap=(
                    float(align_cfg["u_total_cap"])
                    if align_cfg.get("u_total_cap") is not None
                    else 5.0
                ),
                lambda_relative_rec=float(
                    align_cfg.get("lambda_relative_rec", 0.0)
                ),
                tau_relative=float(align_cfg.get("tau_relative", 0.25)),
                relative_start_epoch=int(
                    align_cfg.get("relative_start_epoch", 6)
                ),
                relative_ramp_epochs=int(
                    align_cfg.get("relative_ramp_epochs", 5)
                ),
                relative_weight_min=float(
                    align_cfg.get("relative_weight_min", 0.67)
                ),
                relative_weight_max=float(
                    align_cfg.get("relative_weight_max", 1.50)
                ),
                relative_weight_source=str(
                    align_cfg.get("relative_weight_source", "target")
                ),
                detach_relative_source=bool(
                    align_cfg.get("detach_relative_source", True)
                ),
                normalize_relative_within_celltype=bool(
                    align_cfg.get("normalize_relative_within_celltype", False)
                ),
                minimum_relative_ess_fraction=float(
                    align_cfg.get("minimum_relative_ess_fraction", 0.0)
                ),
                eps=float(align_cfg.get("eps", 1e-6)),
            ),
            v7a_donor_pathology=V7A_DONOR_PATHOLOGY,
            v7a_donor_celltype_mask=V7A_DONOR_CELLTYPE_MASK,
            v7a_donor_pathology_valid=V7A_DONOR_PATHOLOGY_VALID,
            v7a_sample_size=int(profile_cfg.get("sample_size", 50_000)),
            v7a_profile_init_seed=int(profile_cfg.get("init_seed", 42)),
            v7a_diag_print_every=int(v7a_cfg.get("diag_print_every", 100)),
        )

    # n_genes / n_bins are needed by V7a (difficulty profile) and V8 (gene
    # detection-rate buffer). Pass them exactly once whenever either layer is
    # active so V6-only runs are unaffected.
    dims_kwargs: Dict[str, Any] = {}
    if v7a_enabled or v8_enabled:
        dims_kwargs.update(n_genes=N_GENES, n_bins=N_BINS)

    v8_kwargs: Dict[str, Any] = {}
    if v8_enabled:
        v8_kwargs.update(
            v8_enabled=True,
            n_tech_max=N_TECH,
            pi_tech_config=LatentKNNPiTechConfig(
                enabled=bool(pi_tech_cfg.get("enabled", False)),
                mode=str(pi_tech_cfg.get("mode", "in_batch")),
                k=int(pi_tech_cfg.get("k", 32)),
                in_batch_k=(
                    int(pi_tech_cfg["in_batch_k"])
                    if pi_tech_cfg.get("in_batch_k") is not None
                    else None
                ),
                bank_start_epoch=int(pi_tech_cfg.get("bank_start_epoch", 13)),
                bank_ramp_epochs=int(pi_tech_cfg.get("bank_ramp_epochs", 4)),
                bank_refresh_epochs=int(
                    pi_tech_cfg.get("bank_refresh_epochs", 1)
                ),
                bank_per_donor_celltype=int(
                    pi_tech_cfg.get("bank_per_donor_celltype", 2)
                ),
                bank_batch_size=int(pi_tech_cfg.get("bank_batch_size", 48)),
                bank_num_workers=int(pi_tech_cfg.get("bank_num_workers", 0)),
                bank_min_distinct_donors=int(
                    pi_tech_cfg.get("bank_min_distinct_donors", 4)
                ),
                bank_detect_shrinkage_donors=float(
                    pi_tech_cfg.get("bank_detect_shrinkage_donors", 8.0)
                ),
                bank_seed=int(pi_tech_cfg.get("bank_seed", 20260827)),
                bank_exclude_query_donor=bool(
                    pi_tech_cfg.get("bank_exclude_query_donor", True)
                ),
                bank_invalid_fallback=str(
                    pi_tech_cfg.get("bank_invalid_fallback", "zero")
                ),
                bank_match_region=bool(
                    pi_tech_cfg.get("bank_match_region", True)
                ),
                sex_linked_gene_indices=PI_TECH_SEX_LINKED_GENE_INDICES,
                diagnostic_min_positions=int(
                    pi_tech_cfg.get("diagnostic_min_positions", 256)
                ),
                diagnostic_min_cells=int(
                    pi_tech_cfg.get("diagnostic_min_cells", 8)
                ),
                diagnostic_min_distinct_donors=int(
                    pi_tech_cfg.get(
                        "diagnostic_min_distinct_donors", 2
                    )
                ),
                pi_cap=float(pi_tech_cfg.get("pi_cap", 0.6)),
                depth_proxy=str(pi_tech_cfg.get("depth_proxy", "n_detected")),
                depth_low_pct=float(pi_tech_cfg.get("depth_low_pct", 0.25)),
                depth_low_temp=float(pi_tech_cfg.get("depth_low_temp", 0.05)),
                detect_floor=float(pi_tech_cfg.get("detect_floor", 0.05)),
                neighbor_on_hi=float(pi_tech_cfg.get("neighbor_on_hi", 0.5)),
                safety_neighbor_off=float(pi_tech_cfg.get("safety_neighbor_off", 0.2)),
                w_measure=float(pi_tech_cfg.get("w_measure", 0.5)),
                w_neighbor=float(pi_tech_cfg.get("w_neighbor", 0.5)),
                eps=float(pi_tech_cfg.get("eps", 1e-6)),
            ),
            pi_tech_bank_dataset=PI_TECH_TRAIN_DATASET,
            thinning_config=ThinningConfig(
                enabled=bool(thinning_cfg.get("enabled", False)),
                rho=float(thinning_cfg.get("rho", 0.5)),
                lambda_consistency=float(thinning_cfg.get("lambda_consistency", 0.1)),
                consistency_metric=str(thinning_cfg.get("consistency_metric", "smooth_l1")),
                detach_full=bool(thinning_cfg.get("detach_full", True)),
                every_n_steps=int(thinning_cfg.get("every_n_steps", 1)),
                lambda_tech_zero_sup=float(thinning_cfg.get("lambda_tech_zero_sup", 0.0)),
                eps=float(thinning_cfg.get("eps", 1e-6)),
            ),
            subgroup_cvar_config=SubgroupCVaRConfig(
                enabled=bool(cvar_cfg.get("enabled", False)),
                group_by=str(cvar_cfg.get("group_by", "celltype_tech_depth")),
                n_depth_bins=int(cvar_cfg.get("n_depth_bins", 4)),
                metric=str(cvar_cfg.get("metric", "nll")),
                ema_momentum=float(cvar_cfg.get("ema_momentum", 0.9)),
                min_group_cells=int(cvar_cfg.get("min_group_cells", 50)),
                warmup_epochs=int(cvar_cfg.get("warmup_epochs", 3)),
                mad_k=float(cvar_cfg.get("mad_k", 3.0)),
                worst_alpha=float(cvar_cfg.get("worst_alpha", 0.10)),
                persist_epochs=int(cvar_cfg.get("persist_epochs", 2)),
                release_epochs=int(cvar_cfg.get("release_epochs", 3)),
                mode=str(cvar_cfg.get("mode", "reweight")),
                w_max=float(cvar_cfg.get("w_max", 4.0)),
                ramp_epochs=int(cvar_cfg.get("ramp_epochs", 3)),
                lambda_cvar=float(cvar_cfg.get("lambda_cvar", 1.0)),
                eps=float(cvar_cfg.get("eps", 1e-6)),
            ),
            zero_reliability_config=ZeroReliabilityConfig(
                enabled=bool(zero_rel_cfg.get("enabled", False)),
                alpha=float(zero_rel_cfg.get("alpha", 0.5)),
                floor=float(zero_rel_cfg.get("floor", 0.3)),
            ),
            variance_floor_config=VarianceFloorConfig(
                enabled=bool(vf_cfg.get("enabled", False)),
                beta=float(vf_cfg.get("beta", 0.5)),
                lambda_var=float(vf_cfg.get("lambda_var", 0.0)),
                lambda_depth_decorr=float(vf_cfg.get("lambda_depth_decorr", 0.0)),
                min_cells_per_group=int(vf_cfg.get("min_cells_per_group", 16)),
                min_nonzero_frac=float(vf_cfg.get("min_nonzero_frac", 0.1)),
                nonzero_only=bool(vf_cfg.get("nonzero_only", False)),
                target_mode=str(vf_cfg.get("target_mode", "in_batch")),
                min_target_count=int(vf_cfg.get("min_target_count", 50)),
                min_nonzero_cells=int(vf_cfg.get("min_nonzero_cells", 4)),
                n_bins=int(vf_cfg.get("n_bins", N_BINS)),
                eps=float(vf_cfg.get("eps", 1e-6)),
            ),
            high_margin_config=HighMarginRankConfig(
                enabled=bool(hmr_cfg.get("enabled", False)),
                high_tier=int(hmr_cfg.get("high_tier", 4)),
                low_tiers=tuple(hmr_cfg.get("low_tiers", [1, 2])),
                margin=float(hmr_cfg.get("margin", 1.0)),
                lambda_high_margin=float(hmr_cfg.get("lambda_high_margin", 0.0)),
                min_high_cells=int(hmr_cfg.get("min_high_cells", 2)),
                min_low_ref_count=int(hmr_cfg.get("min_low_ref_count", 8)),
                ema_momentum=float(hmr_cfg.get("ema_momentum", 0.9)),
                t4_floor_target=float(hmr_cfg.get("t4_floor_target", 0.0)),
                lambda_t4_floor=float(hmr_cfg.get("lambda_t4_floor", 0.0)),
                lambda_tier4_focal=float(hmr_cfg.get("lambda_tier4_focal", 0.0)),
                lambda_zero_focal=float(hmr_cfg.get("lambda_zero_focal", 0.0)),
                focal_gamma=float(hmr_cfg.get("focal_gamma", 2.0)),
                tier4_focal_pos_weight=float(hmr_cfg.get("tier4_focal_pos_weight", 1.0)),
                tier4_focal_neg_weight=float(hmr_cfg.get("tier4_focal_neg_weight", 1.0)),
                eps=float(hmr_cfg.get("eps", 1e-6)),
            ),
            tech_invariance_config=TechInvarianceConfig(
                enabled=bool(tech_inv_cfg.get("enabled", False)),
                every_n_steps=int(tech_inv_cfg.get("every_n_steps", 1)),
                lambda_z_cons=float(tech_inv_cfg.get("lambda_z_cons", 0.0)),
                lambda_bam_tech=float(tech_inv_cfg.get("lambda_bam_tech", 0.0)),
                lambda_rec_cons=float(tech_inv_cfg.get("lambda_rec_cons", 0.0)),
                lambda_disease_ctrl=float(tech_inv_cfg.get("lambda_disease_ctrl", 0.0)),
                tech_kinds=tuple(tech_inv_cfg.get("tech_kinds", ["dropout", "gain", "ambient"])),
                dropout_p=float(tech_inv_cfg.get("dropout_p", 0.35)),
                gain_lo=float(tech_inv_cfg.get("gain_lo", 0.5)),
                gain_hi=float(tech_inv_cfg.get("gain_hi", 1.8)),
                ambient_a=float(tech_inv_cfg.get("ambient_a", 0.4)),
                disease_gamma=float(tech_inv_cfg.get("disease_gamma", 1.0)),
                rec_metric=str(tech_inv_cfg.get("rec_metric", "etier")),
                detach_clean=bool(tech_inv_cfg.get("detach_clean", True)),
                warmup_steps=int(tech_inv_cfg.get("warmup_steps", 0)),
            ),
        )

    trainer = trainer_cls(
        *args,
        lambda_ref_center_db=(
            float(ref_cfg.get("lambda", 0.0)) if bool(ref_cfg.get("enabled", False)) else 0.0
        ),
        ref_center_min_cells_per_donor=int(ref_cfg.get("min_cells_per_donor", 2)),
        ref_center_min_donors_per_celltype=int(ref_cfg.get("min_donors_per_celltype", 2)),
        lambda_tech_adv=(
            float(adv_cfg.get("lambda", 0.0)) if bool(adv_cfg.get("enabled", False)) else 0.0
        ),
        tech_adv_warmup_epochs=int(adv_cfg.get("warmup_epochs", 5)),
        tech_adv_ramp_epochs=int(adv_cfg.get("ramp_epochs", 0)),
        tech_adv_max_lambda=float(adv_cfg.get("max_lambda", adv_cfg.get("lambda", 0.0))),
        lambda_celltype_adv=(
            float(ct_adv_cfg.get("lambda", 0.0)) if bool(ct_adv_cfg.get("enabled", False)) else 0.0
        ),
        celltype_adv_warmup_epochs=int(ct_adv_cfg.get("warmup_epochs", 5)),
        celltype_adv_ramp_epochs=int(ct_adv_cfg.get("ramp_epochs", 0)),
        celltype_adv_max_lambda=float(ct_adv_cfg.get("max_lambda", ct_adv_cfg.get("lambda", 0.0))),
        # v31 sex adversary (v3.sex_adversary) — soft sex erasure on z_clean.
        lambda_sex_adv=(
            float(sex_adv_cfg.get("lambda", 0.0)) if bool(sex_adv_cfg.get("enabled", False)) else 0.0
        ),
        sex_adv_warmup_epochs=int(sex_adv_cfg.get("warmup_epochs", 5)),
        sex_adv_ramp_epochs=int(sex_adv_cfg.get("ramp_epochs", 0)),
        sex_adv_max_lambda=float(sex_adv_cfg.get("max_lambda", sex_adv_cfg.get("lambda", 0.0))),
        v3_diag_print_every=int(
            _section("v3").get("diag_print_every", v3_diag_default)
        ),
        # v15 auxiliary reference batch (v3.reference_celltype_adversary.aux_batch).
        # OFF unless the adversary is enabled AND aux_batch.enabled is true; an
        # absent section ⇒ AuxReferenceAdversaryConfig() (enabled=False) ⇒ the
        # trainer never builds the aux runner (byte-identical v14 path).
        aux_celltype_adv_config=_build_aux_celltype_adv_config(ct_adv_cfg),
        ordinal_class_weights=ORDINAL_CLASS_WEIGHTS,
        ordinal_balance_config=OrdinalBalanceConfig(
            enabled=bool(ord_cfg.get("enabled", False)),
            lambda_balanced=float(ord_cfg.get("lambda_balanced", 0.0)),
            lambda_nonzero=float(ord_cfg.get("lambda_nonzero", 0.0)),
            class_weight_mode=str(ord_cfg.get("class_weight_mode", "sqrt_inverse")),
            class_weight_cap=float(ord_cfg.get("class_weight_cap", 6.0)),
            zero_bin_weight_floor=float(ord_cfg.get("zero_bin_weight_floor", 0.2)),
            eps=float(ord_cfg.get("eps", 1e-6)),
        ),
        hierarchical_ordinal_config=HierarchicalOrdinalConfig(
            enabled=bool(hier_cfg.get("enabled", False)),
            lambda_zero=float(hier_cfg.get("lambda_zero", 0.0)),
            lambda_group=float(hier_cfg.get("lambda_group", 0.0)),
            lambda_emd=float(hier_cfg.get("lambda_emd", 0.0)),
            zero_pos_weight=float(hier_cfg.get("zero_pos_weight", 2.0)),
            group_definition=tuple(
                tuple(int(b) for b in g)
                for g in hier_cfg.get("group_definition", [[1, 2], [3, 4], [5, 6]])
            ),
            group_weight_mode=str(hier_cfg.get("group_weight_mode", "sqrt_inverse")),
            group_weight_cap=float(hier_cfg.get("group_weight_cap", 4.0)),
            group_class_weights=(
                tuple(float(w) for w in hier_cfg["group_class_weights"])
                if hier_cfg.get("group_class_weights") is not None
                else None
            ),
            lambda_within_group=float(hier_cfg.get("lambda_within_group", 0.0)),
            ramp_epochs=int(hier_cfg.get("ramp_epochs", 0)),
            # π_tech soft-zero strengths (read from the top-level pi_tech section
            # so all π_tech knobs live in one place). 0.0 ⇒ no effect.
            pi_soft_target_zero=float(pi_tech_cfg.get("pi_soft_target_zero", 0.0)),
            pi_downweight_zero=float(pi_tech_cfg.get("pi_downweight_zero", 0.0)),
            # Asymmetric under-prediction penalty (nonzero_exact lever). 0 ⇒ off.
            lambda_under_asym=float(hier_cfg.get("lambda_under_asym", 0.0)),
            asym_tau=float(hier_cfg.get("asym_tau", 0.7)),
            eps=float(hier_cfg.get("eps", 1e-8)),
        ),
        global_ref_config=GlobalReferenceBankConfig(
            enabled=bool(global_ref_cfg.get("enabled", False)),
            lambda_global=float(global_ref_cfg.get("lambda", 0.0)),
            ema_momentum=float(global_ref_cfg.get("ema_momentum", 0.95)),
            shrinkage_k=float(global_ref_cfg.get("shrinkage_k", 8.0)),
            min_reference_cells=int(global_ref_cfg.get("min_reference_cells", 1)),
            update_during_eval=bool(global_ref_cfg.get("update_during_eval", False)),
        ),
        uncertainty_calibration_config=UncertaintyCalibrationConfig(
            enabled=bool(v4_unc_cfg.get("enabled", False)),
            lambda_rec_alignment=float(v4_unc_cfg.get("lambda_rec_alignment", 0.0)),
            lambda_relative_rec=float(v4_unc_cfg.get("lambda_relative_rec", 0.0)),
            tau_relative=float(v4_unc_cfg.get("tau_relative", 0.1)),
            group_by=str(v4_unc_cfg.get("group_by", "celltype_tech_depth")),
            n_depth_bins=int(v4_unc_cfg.get("n_depth_bins", 4)),
            min_group_size=int(v4_unc_cfg.get("min_group_size", 4)),
            clamp_z=float(v4_unc_cfg.get("clamp_z", 3.0)),
            eps=float(v4_unc_cfg.get("eps", 1e-6)),
        ),
        v6_reference_prior_config=ReferenceAnchoredPriorConfig(
            enabled=bool(ref_prior_cfg.get("enabled", False)),
            freeze_mu=bool(ref_prior_cfg.get("freeze_mu", True)),
            refresh_before_epoch=bool(ref_prior_cfg.get("refresh_before_epoch", True)),
            refresh_every_epochs=int(ref_prior_cfg.get("refresh_every_epochs", 1)),
            refresh_batch_size=int(ref_prior_cfg.get("refresh_batch_size", 64)),
            refresh_num_workers=int(ref_prior_cfg.get("refresh_num_workers", 0)),
            max_refresh_cells=int(ref_prior_cfg.get("max_refresh_cells", 131072)),
            ema_momentum=float(ref_prior_cfg.get("ema_momentum", 0.8)),
            min_cells_per_donor=int(ref_prior_cfg.get("min_cells_per_donor", 1)),
            min_donors_per_celltype=int(ref_prior_cfg.get("min_donors_per_celltype", 2)),
            min_reference_cells=int(ref_prior_cfg.get("min_reference_cells", 4)),
            sync_ddp=bool(ref_prior_cfg.get("sync_ddp", True)),
            eps=float(ref_prior_cfg.get("eps", 1e-6)),
        ),
        v6_state_usage_config=V5StateUsageConfig(
            enabled=bool(state_cfg.get("enabled", False)),
            lambda_state_abs=float(state_cfg.get("lambda_state_abs", 0.0)),
            target_state_abs=float(state_cfg.get("target_state_abs", 0.0)),
            lambda_state_fraction=float(state_cfg.get("lambda_state_fraction", 0.0)),
            target_state_fraction=float(state_cfg.get("target_state_fraction", 0.0)),
            lambda_tech_to_state=float(state_cfg.get("lambda_tech_to_state", 0.0)),
            max_tech_to_state=float(state_cfg.get("max_tech_to_state", 0.0)),
            eps=float(state_cfg.get("eps", 1e-6)),
        ),
        v6_uncertainty_config=V5UncertaintyResidualConfig(
            enabled=bool(unc_cfg.get("enabled", False)),
            lambda_relative_rec=float(unc_cfg.get("lambda_relative_rec", 0.0)),
            lambda_rec_alignment=float(unc_cfg.get("lambda_rec_alignment", 0.0)),
            lambda_depth_corr=float(unc_cfg.get("lambda_depth_corr", 0.0)),
            lambda_saturation=float(unc_cfg.get("lambda_saturation", 0.0)),
            tau_relative=float(unc_cfg.get("tau_relative", 0.1)),
            group_by=str(unc_cfg.get("group_by", "celltype_tech_depth")),
            n_depth_bins=int(unc_cfg.get("n_depth_bins", 4)),
            min_group_size=int(unc_cfg.get("min_group_size", 8)),
            clamp_z=float(unc_cfg.get("clamp_z", 2.5)),
            raw_std_floor=float(unc_cfg.get("raw_std_floor", 0.0)),
            detach_weight_uncertainty=bool(unc_cfg.get("detach_weight_uncertainty", True)),
            eps=float(unc_cfg.get("eps", 1e-6)),
        ),
        # v6_min_mixer regularizers. Only fire when the corresponding
        # decoder.use_score_residual_mixer / coeff_zero_mean_penalty paths
        # are active; otherwise lambda=0 means a no-op.
        v6_lambda_gate_balance=float(mixer_cfg.get("lambda_gate_balance", 0.0)),
        v6_lambda_gate_ref_safe=float(mixer_cfg.get("lambda_gate_ref_safe", 0.0)),
        v6_lambda_coeff_zero_mean=float(mixer_cfg.get("lambda_coeff_zero_mean", 0.0)),
        n_celltypes=N_CELLTYPES,
        n_donors=N_DONORS,
        d_z=D_Z,
        v4_diag_print_every=int(_section("v4").get("diag_print_every", 100)),
        v6_diag_print_every=int(_section("v6").get("diag_print_every", 100)),
        **dims_kwargs,
        **v7a_kwargs,
        **v8_kwargs,
        **kwargs,
    )
    _warm_start_trainer_state_from_checkpoint(trainer)
    return trainer


def install_hooks() -> None:
    base.OrdinalScDataset = V3OrdinalScDataset
    base.build_datasets = build_datasets
    base.build_loaders = build_loaders
    base.build_system = build_system
    base.build_criterion = build_criterion
    base.build_optimizer = build_optimizer
    base.Trainer = make_trainer


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run integrated KMLEE-BAM training.")
    p.add_argument("--config", type=str, required=True, help="Path to KMLEE-BAM training config JSON.")
    return p


def main() -> None:
    args = build_argparser().parse_args()
    with open(args.config, "r", encoding="utf-8") as f:
        raw = json.load(f)
    set_settings(raw)
    _enforce_launch_guard()
    install_hooks()
    cfg = base.load_config(args.config)

    out_dir = Path(cfg.train.out_dir)
    rank = int(os.environ.get("RANK", "0"))
    if rank == 0:
        out_dir.mkdir(parents=True, exist_ok=True)
        # Preserve the complete launch input, including top-level provenance
        # sections that the legacy RunConfig dataclass does not model.
        with open(out_dir / "source_config.json", "w", encoding="utf-8") as f:
            json.dump(raw, f, indent=2, ensure_ascii=False)
            f.write("\n")
        _write_runtime_manifest(out_dir, args.config)
        # Settings sidecar — include v7a + v8 alongside the legacy v3/v4/v6
        # blocks so reruns / latent extraction tools can recover the full
        # training config from disk.
        with open(out_dir / "kmlee_bam_settings.json", "w", encoding="utf-8") as f:
            json.dump(
                    {
                        k: raw.get(k, {})
                        for k in (
                            "v3",
                            "v4",
                            "v6",
                            "v7a",
                            "v8",
                            "pi_tech",
                            "thinning",
                            "subgroup_cvar",
                            "zero_reliability",
                            "tech_invariance",
                            "generator_budget_search",
                            "learned_generator_count",
                        )
                    },
                f,
                indent=2,
                ensure_ascii=False,
            )
            f.write("\n")
        # Compatibility sidecars for old latent extraction tooling.
        for name in (
            "v3",
            "v4",
            "v6",
            "v7a",
            "v8",
            "pi_tech",
            "thinning",
            "subgroup_cvar",
            "zero_reliability",
            "generator_budget_search",
            "learned_generator_count",
        ):
            with open(out_dir / f"{name}_settings.json", "w", encoding="utf-8") as f:
                json.dump(raw.get(name, {}), f, indent=2, ensure_ascii=False)
                f.write("\n")

    base.run_training(cfg)


if __name__ == "__main__":
    main()
