"""Production grouped updates for the PRISM module-rescue objective.

Ordinary random minibatches are invalid for donor-centred module rescue:
the donor centre would change with whichever donors happened to enter the
minibatch.  This module builds one complete cell-type contrast at a time,
sampling the same number of cells from every eligible weight-training donor.

The update is deliberately auxiliary to the ordinary end-to-end epoch.  It
uses the same optimiser, but no scheduler step, and therefore adds a small,
auditable number of donor-balanced module gradients per epoch.
"""

from __future__ import annotations

from contextlib import nullcontext
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import torch
import torch.distributed as dist
import torch.nn.functional as F
from torch.utils.data import default_collate

from kmlee_bam.objectives.prism_module_rescue import (
    GroupedModuleActivity,
    PrismModuleRescueConfig,
    PrismModuleRescueOutput,
    dual_view_module_rescue_loss,
    dual_view_grouped_module_rescue_loss,
    grouped_expected_module_activity,
)


PRISM_NON_NEURONAL_CELLTYPE_NAMES = (
    "Astrocyte",
    "Endothelial",
    "Microglia-PVM",
    "OPC",
    "Oligodendrocyte",
    "VLMC",
)
PRISM_NEURONAL_CELLTYPE_NAMES = (
    "Chandelier",
    "L2/3 IT",
    "L4 IT",
    "L5 ET",
    "L5 IT",
    "L5/6 NP",
    "L6 CT",
    "L6 IT",
    "L6 IT Car3",
    "L6b",
    "Lamp5",
    "Lamp5 Lhx6",
    "Pax6",
    "Pvalb",
    "Sncg",
    "Sst",
    "Sst Chodl",
    "Vip",
)
LINEAGE_BALANCE_UNIFORM = "uniform"
LINEAGE_BALANCE_PRISM_18_6 = "prism_18_neuronal_6_non_neuronal_equal_mass"
NON_NEURONAL_PSEUDOBULK_LEGACY = "legacy"
NON_NEURONAL_PSEUDOBULK_SEPARATE_CONTROL = "separate_control"
NON_NEURONAL_PSEUDOBULK_AGGREGATE = "aggregate"


@dataclass
class PrismModuleRescueTrainingConfig:
    """Runtime configuration; absent/disabled remains a strict no-op."""

    enabled: bool = False
    stats_path: str | None = None
    expected_stats_sha256: str | None = None
    lambda_module: float = 0.10
    branch_fraction: float = 2.0 / 3.0
    huber_delta: float = 1.0
    ccc_weight: float = 0.0
    steps_per_epoch: int = 4
    draws_per_celltype_per_epoch: int = 0
    cells_per_donor_celltype: int = 2
    minimum_donors_per_celltype: int = 8
    minimum_target_cells: int = 10
    scale_floor: float = 0.05
    scale_shrinkage: float = 0.10
    start_epoch: int = 1
    warmup_epochs: int = 2
    seed: int = 420731
    nonoverlap_sampling: bool = False
    deterministic_forward: bool = False
    restrict_parameter_updates: bool = False
    # Opt-in v3 route: after the configured epoch, allow the rescue-only full
    # view to teach the persistent pooling readout and posterior mapping.  The
    # module Transformer, tokenizer, cell-type embedding, and nuisances remain
    # excluded.
    allow_state_readout_updates: bool = False
    state_readout_start_epoch: int = 4
    contract_version: str = "legacy"
    global_blocks_per_optimizer_step: int = 0
    equal_donor_rescue_gradient: bool = False
    ccc_target_variance_floor: float = 1e-4
    ccc_calibration_sha256: str | None = None
    # Optional, fail-closed PRISM lineage contract.  The named mode uses the
    # six canonical glial/vascular cell types above and requires all 24 cell
    # types to be rescue-eligible.  It gives the 18 neuronal types raw weight
    # 1 and the six non-neuronal types raw weight 3, then normalizes the mean
    # block weight back to one.
    lineage_balance_mode: str = LINEAGE_BALANCE_UNIFORM
    # v6A is opt-in.  Both new modes expose the same four deterministic K=2
    # draws (eight nuclei per donor×celltype) in the same packed schedule.
    # ``separate_control`` scores four K=2 losses; ``aggregate`` scores one K=8
    # pseudobulk loss and replays its exact group-mean VJP.
    non_neuronal_pseudobulk_mode: str = NON_NEURONAL_PSEUDOBULK_LEGACY
    non_neuronal_aggregate_draws: int = 4
    # Train-only named-pathology coefficient rescue and weakest-celltype tail
    # protection.  Both are opt-in and start only after the decoder/PRISM path
    # has matured.
    pathology_axis_enabled: bool = False
    pathology_axis_start_epoch: int = 25
    lambda_pathology_axis: float = 0.05
    pathology_axis_magnitude_weight: float = 0.10
    celltype_tail_enabled: bool = False
    celltype_tail_start_epoch: int = 25
    celltype_deficit_ema_decay: float = 0.80
    celltype_weight_min: float = 0.50
    celltype_weight_max: float = 2.00
    celltype_tail_fraction: float = 0.25
    ad_recovery_baseline_by_celltype: dict[str, float] = field(default_factory=dict)

    def validate(self) -> None:
        if bool(self.enabled) and not self.stats_path:
            raise ValueError(
                "prism_module_rescue.stats_path is required when enabled"
            )
        for name in (
            "steps_per_epoch",
            "cells_per_donor_celltype",
            "minimum_donors_per_celltype",
            "minimum_target_cells",
        ):
            if int(getattr(self, name)) <= 0:
                raise ValueError(f"prism_module_rescue.{name} must be positive")
        if int(self.start_epoch) <= 0:
            raise ValueError(
                "prism_module_rescue.start_epoch must be positive"
            )
        if int(self.warmup_epochs) < 0:
            raise ValueError(
                "prism_module_rescue.warmup_epochs must be non-negative"
            )
        if int(self.state_readout_start_epoch) <= 0:
            raise ValueError(
                "prism_module_rescue.state_readout_start_epoch must be positive"
            )
        if int(self.pathology_axis_start_epoch) <= 0:
            raise ValueError("pathology_axis_start_epoch must be positive")
        if int(self.celltype_tail_start_epoch) <= 0:
            raise ValueError("celltype_tail_start_epoch must be positive")
        if not math.isfinite(float(self.lambda_pathology_axis)) or float(
            self.lambda_pathology_axis
        ) < 0.0:
            raise ValueError("lambda_pathology_axis must be finite/non-negative")
        if not 0.0 <= float(self.pathology_axis_magnitude_weight) <= 1.0:
            raise ValueError("pathology_axis_magnitude_weight must lie in [0, 1]")
        if not 0.0 <= float(self.celltype_deficit_ema_decay) < 1.0:
            raise ValueError("celltype_deficit_ema_decay must lie in [0, 1)")
        if not (
            0.0 < float(self.celltype_weight_min)
            <= 1.0
            <= float(self.celltype_weight_max)
        ):
            raise ValueError("celltype weight bounds must bracket one")
        if not 0.0 < float(self.celltype_tail_fraction) <= 1.0:
            raise ValueError("celltype_tail_fraction must lie in (0, 1]")
        if int(self.draws_per_celltype_per_epoch) < 0:
            raise ValueError(
                "prism_module_rescue.draws_per_celltype_per_epoch must be "
                "non-negative"
            )
        if not math.isfinite(float(self.lambda_module)) or float(
            self.lambda_module
        ) < 0.0:
            raise ValueError(
                "prism_module_rescue.lambda_module must be finite/non-negative"
            )
        if not 0.0 <= float(self.branch_fraction) <= 1.0:
            raise ValueError(
                "prism_module_rescue.branch_fraction must lie in [0, 1]"
            )

        if float(self.huber_delta) <= 0.0:
            raise ValueError(
                "prism_module_rescue.huber_delta must be positive"
            )
        if not math.isfinite(float(self.ccc_weight)) or float(
            self.ccc_weight
        ) < 0.0:
            raise ValueError(
                "prism_module_rescue.ccc_weight must be finite/non-negative"
            )
        if float(self.scale_floor) <= 0.0:
            raise ValueError(
                "prism_module_rescue.scale_floor must be positive"
            )
        if not 0.0 <= float(self.scale_shrinkage) <= 1.0:
            raise ValueError(
                "prism_module_rescue.scale_shrinkage must lie in [0, 1]"
            )
        if int(self.global_blocks_per_optimizer_step) < 0:
            raise ValueError(
                "prism_module_rescue.global_blocks_per_optimizer_step must "
                "be non-negative"
            )
        if not math.isfinite(float(self.ccc_target_variance_floor)) or float(
            self.ccc_target_variance_floor
        ) <= 0.0:
            raise ValueError(
                "prism_module_rescue.ccc_target_variance_floor must be "
                "finite and positive"
            )
        contract = str(self.contract_version)
        if contract not in {"legacy", "v2"}:
            raise ValueError(
                "prism_module_rescue.contract_version must be 'legacy' or 'v2'"
            )
        lineage_mode = str(self.lineage_balance_mode)
        if lineage_mode not in {
            LINEAGE_BALANCE_UNIFORM,
            LINEAGE_BALANCE_PRISM_18_6,
        }:
            raise ValueError(
                "prism_module_rescue.lineage_balance_mode must be "
                f"'{LINEAGE_BALANCE_UNIFORM}' or "
                f"'{LINEAGE_BALANCE_PRISM_18_6}'"
            )
        pseudobulk_mode = str(self.non_neuronal_pseudobulk_mode)
        allowed_pseudobulk_modes = {
            NON_NEURONAL_PSEUDOBULK_LEGACY,
            NON_NEURONAL_PSEUDOBULK_SEPARATE_CONTROL,
            NON_NEURONAL_PSEUDOBULK_AGGREGATE,
        }
        if pseudobulk_mode not in allowed_pseudobulk_modes:
            raise ValueError(
                "prism_module_rescue.non_neuronal_pseudobulk_mode must be "
                f"one of {sorted(allowed_pseudobulk_modes)}"
            )
        if int(self.non_neuronal_aggregate_draws) <= 0:
            raise ValueError(
                "prism_module_rescue.non_neuronal_aggregate_draws must be "
                "positive"
            )
        if pseudobulk_mode != NON_NEURONAL_PSEUDOBULK_LEGACY:
            required_v6a = {
                "enabled": bool(self.enabled),
                "contract_version_v2": contract == "v2",
                "lineage_balance_18_6": (
                    lineage_mode == LINEAGE_BALANCE_PRISM_18_6
                ),
                "deterministic_forward": bool(self.deterministic_forward),
                "restrict_parameter_updates": bool(
                    self.restrict_parameter_updates
                ),
                "nonoverlap_sampling": bool(self.nonoverlap_sampling),
            }
            failed_v6a = [
                name for name, enabled in required_v6a.items() if not enabled
            ]
            if failed_v6a:
                raise ValueError(
                    "v6A non-neuronal pseudobulk requires: "
                    + ", ".join(failed_v6a)
                )
            if int(self.cells_per_donor_celltype) != 2:
                raise ValueError(
                    "v6A fixes cells_per_donor_celltype=2 so control and "
                    "aggregate use the same nuclei"
                )
            if int(self.non_neuronal_aggregate_draws) != 4:
                raise ValueError(
                    "the production v6A K=8 contract fixes "
                    "non_neuronal_aggregate_draws=4"
                )
            if int(self.draws_per_celltype_per_epoch) != int(
                self.non_neuronal_aggregate_draws
            ):
                raise ValueError(
                    "v6A requires draws_per_celltype_per_epoch to equal "
                    "non_neuronal_aggregate_draws"
                )
            if int(self.global_blocks_per_optimizer_step) != 12:
                raise ValueError(
                    "v6A fixes global_blocks_per_optimizer_step=12"
                )
        if lineage_mode != LINEAGE_BALANCE_UNIFORM:
            if not bool(self.enabled):
                raise ValueError(
                    "lineage-balanced module rescue requires enabled=True"
                )
            if contract != "v2":
                raise ValueError(
                    "lineage-balanced module rescue requires contract_version='v2'"
                )
        if contract == "v2" and bool(self.enabled):
            if self.expected_stats_sha256 is not None:
                digest = str(self.expected_stats_sha256).lower()
                if len(digest) != 64 or any(
                    character not in "0123456789abcdef"
                    for character in digest
                ):
                    raise ValueError(
                        "prism_module_rescue.expected_stats_sha256 must be "
                        "a 64-character hex digest"
                    )
            required_flags = {
                "nonoverlap_sampling": bool(self.nonoverlap_sampling),
                "deterministic_forward": bool(self.deterministic_forward),
                "restrict_parameter_updates": bool(
                    self.restrict_parameter_updates
                ),
                "equal_donor_rescue_gradient": bool(
                    self.equal_donor_rescue_gradient
                ),
            }
            disabled = [name for name, value in required_flags.items() if not value]
            if disabled:
                raise ValueError(
                    "PRISM module-rescue v2 requires enabled flags: "
                    + ", ".join(disabled)
                )
            if int(self.draws_per_celltype_per_epoch) <= 0:
                raise ValueError(
                    "PRISM module-rescue v2 requires a positive "
                    "draws_per_celltype_per_epoch"
                )
            if int(self.global_blocks_per_optimizer_step) <= 0:
                raise ValueError(
                    "PRISM module-rescue v2 requires a positive "
                    "global_blocks_per_optimizer_step"
                )
            if abs(float(self.branch_fraction) - 0.25) > 1e-12:
                raise ValueError(
                    "PRISM module-rescue v2 fixes branch_fraction=0.25"
                )
            if float(self.ccc_weight) > 0.0:
                digest = str(self.ccc_calibration_sha256 or "").lower()
                if len(digest) != 64 or any(
                    character not in "0123456789abcdef"
                    for character in digest
                ):
                    raise ValueError(
                        "CCC-enabled v2 requires a 64-character hex "
                        "ccc_calibration_sha256 from train-only gradient "
                        "calibration"
                    )

    def ramp_at_epoch(self, epoch: int) -> float:
        """Return zero before start, then a 1-based linear rescue ramp."""

        if int(epoch) < int(self.start_epoch):
            return 0.0
        if int(self.warmup_epochs) <= 0:
            return 1.0
        return float(
            min(
                1.0,
                (int(epoch) - int(self.start_epoch) + 1)
                / float(self.warmup_epochs),
            )
        )


@dataclass(frozen=True)
class LineageBalancePlan:
    """Resolved weights and GPU-count-independent global block order."""

    mode: str
    celltype_raw_weights: tuple[float, ...]
    celltype_normalized_weights: tuple[float, ...]
    neuronal_celltype_ids: tuple[int, ...]
    non_neuronal_celltype_ids: tuple[int, ...]
    global_schedule: tuple[tuple[int, int], ...]
    neuronal_blocks_per_update: int
    non_neuronal_blocks_per_update: int
    neuronal_blocks_by_update: tuple[int, ...] = ()
    non_neuronal_blocks_by_update: tuple[int, ...] = ()
    non_neuronal_pseudobulk_mode: str = NON_NEURONAL_PSEUDOBULK_LEGACY
    aggregate_draws: int = 1

    @property
    def enabled(self) -> bool:
        return self.mode != LINEAGE_BALANCE_UNIFORM


def _interleave_lineage_update(
    neuronal: list[tuple[int, int]],
    non_neuronal: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Spread non-neuronal blocks through one already-balanced update."""

    total = len(neuronal) + len(non_neuronal)
    if total <= 0:
        return []
    output: list[tuple[int, int]] = []
    neuronal_cursor = 0
    non_neuronal_cursor = 0
    n_non_neuronal = len(non_neuronal)
    for position in range(total):
        previous_target = (position * n_non_neuronal) // total
        next_target = ((position + 1) * n_non_neuronal) // total
        if next_target > previous_target:
            output.append(non_neuronal[non_neuronal_cursor])
            non_neuronal_cursor += 1
        else:
            output.append(neuronal[neuronal_cursor])
            neuronal_cursor += 1
    if (
        neuronal_cursor != len(neuronal)
        or non_neuronal_cursor != len(non_neuronal)
    ):
        raise RuntimeError("lineage interleaving did not consume every block")
    return output


def build_lineage_balance_plan(
    *,
    config: PrismModuleRescueTrainingConfig,
    celltype_names: Iterable[str],
    eligible_celltypes: Iterable[int],
) -> LineageBalancePlan:
    """Resolve and audit the optional 18-neuron/6-non-neuron contract.

    The schedule is expressed in global-block order, before DDP shards blocks
    by ``global_block = local_slot * world_size + rank``.  Consequently every
    optimizer update sees the same lineage composition on 4 or 6 GPUs.
    """

    names = tuple(str(value) for value in celltype_names)
    eligible = tuple(int(value) for value in eligible_celltypes)
    if len(set(names)) != len(names):
        raise ValueError("module-rescue celltype names must be unique")
    if len(set(eligible)) != len(eligible):
        raise ValueError("module-rescue eligible celltype ids must be unique")
    invalid_ids = [value for value in eligible if not 0 <= value < len(names)]
    if invalid_ids:
        raise ValueError(
            "module-rescue eligible celltype ids are out of range: "
            f"{invalid_ids}"
        )

    mode = str(config.lineage_balance_mode)
    if mode == LINEAGE_BALANCE_UNIFORM:
        return LineageBalancePlan(
            mode=mode,
            celltype_raw_weights=tuple(1.0 for _ in names),
            celltype_normalized_weights=tuple(1.0 for _ in names),
            neuronal_celltype_ids=(),
            non_neuronal_celltype_ids=(),
            global_schedule=(),
            neuronal_blocks_per_update=0,
            non_neuronal_blocks_per_update=0,
            non_neuronal_pseudobulk_mode=str(
                config.non_neuronal_pseudobulk_mode
            ),
            aggregate_draws=int(config.non_neuronal_aggregate_draws),
        )
    if mode != LINEAGE_BALANCE_PRISM_18_6:
        # ``validate`` normally catches this first, but this resolver is also
        # used directly by focused preflight tests and must fail closed itself.
        raise ValueError(f"unsupported lineage_balance_mode: {mode}")

    canonical_neuronal = set(PRISM_NEURONAL_CELLTYPE_NAMES)
    canonical_non_neuronal = set(PRISM_NON_NEURONAL_CELLTYPE_NAMES)
    canonical_all = canonical_neuronal | canonical_non_neuronal
    available = set(names)
    missing_names = sorted(canonical_all - available)
    unknown_names = sorted(available - canonical_all)
    if missing_names or unknown_names:
        raise ValueError(
            "lineage-balanced module rescue celltype vocabulary differs "
            "from the canonical 18+6 names: "
            f"missing={missing_names}, unknown={unknown_names}"
        )
    eligible_set = set(eligible)
    missing_eligible = [
        names[index]
        for index in range(len(names))
        if index not in eligible_set
    ]
    if missing_eligible:
        raise ValueError(
            "lineage-balanced module rescue requires every named cell type "
            f"to be eligible; missing={missing_eligible}"
        )
    if len(names) != 24:
        raise ValueError(
            "the PRISM 18+6 lineage contract requires exactly 24 named "
            f"cell types, found {len(names)}"
        )

    non_neuronal_ids = tuple(
        index for index, name in enumerate(names)
        if name in canonical_non_neuronal
    )
    neuronal_ids = tuple(
        index for index, name in enumerate(names)
        if name in canonical_neuronal
    )
    if len(neuronal_ids) != 18 or len(non_neuronal_ids) != 6:
        raise ValueError(
            "the PRISM lineage contract requires 18 neuronal and six named "
            "non-neuronal cell types; found "
            f"{len(neuronal_ids)} and {len(non_neuronal_ids)}"
        )

    draws = int(config.draws_per_celltype_per_epoch)
    blocks_per_update = int(config.global_blocks_per_optimizer_step)
    if draws <= 0 or blocks_per_update <= 0:
        raise ValueError(
            "lineage-balanced scheduling requires positive draws and "
            "global_blocks_per_optimizer_step"
        )
    total_blocks = len(names) * draws
    if total_blocks % blocks_per_update != 0:
        raise ValueError(
            "lineage-balanced global blocks must divide into complete "
            f"optimizer updates: {total_blocks} / {blocks_per_update}"
        )
    updates = total_blocks // blocks_per_update
    neuronal_blocks = len(neuronal_ids) * draws
    non_neuronal_blocks = len(non_neuronal_ids) * draws
    if neuronal_blocks % updates != 0 or non_neuronal_blocks % updates != 0:
        raise ValueError(
            "each optimizer update cannot receive an identical lineage mix: "
            f"neuronal={neuronal_blocks}, non_neuronal={non_neuronal_blocks}, "
            f"updates={updates}"
        )
    neuronal_per_update = neuronal_blocks // updates
    non_neuronal_per_update = non_neuronal_blocks // updates
    if neuronal_per_update + non_neuronal_per_update != blocks_per_update:
        raise RuntimeError("resolved lineage blocks do not fill an update")

    neuronal_stream = [
        (celltype, draw)
        for draw in range(draws)
        for celltype in neuronal_ids
    ]
    non_neuronal_stream = [
        (celltype, draw)
        for draw in range(draws)
        for celltype in non_neuronal_ids
    ]
    pseudobulk_mode = str(config.non_neuronal_pseudobulk_mode)
    schedule: list[tuple[int, int]] = []
    neuronal_by_update: list[int] = []
    non_neuronal_by_update: list[int] = []
    if pseudobulk_mode == NON_NEURONAL_PSEUDOBULK_LEGACY:
        for update in range(updates):
            n_start = update * neuronal_per_update
            g_start = update * non_neuronal_per_update
            update_blocks = _interleave_lineage_update(
                neuronal_stream[n_start : n_start + neuronal_per_update],
                non_neuronal_stream[
                    g_start : g_start + non_neuronal_per_update
                ],
            )
            if len(update_blocks) != blocks_per_update:
                raise RuntimeError("lineage schedule produced a short update")
            schedule.extend(update_blocks)
            neuronal_by_update.append(neuronal_per_update)
            non_neuronal_by_update.append(non_neuronal_per_update)
    else:
        aggregate_draws = int(config.non_neuronal_aggregate_draws)
        if updates < len(non_neuronal_ids):
            raise ValueError(
                "v6A schedule needs at least one optimizer update per "
                "non-neuronal cell type"
            )
        neuronal_cursor = 0
        # Interleave the two neuronal-only updates so Adam does not see all
        # six non-neuronal units before the final two updates.  This fixed
        # pattern is shared by control and aggregate arms.
        neuronal_only_updates = {2, 5}
        mixed_updates = [
            update for update in range(updates)
            if update not in neuronal_only_updates
        ]
        if len(mixed_updates) != len(non_neuronal_ids):
            raise RuntimeError("v6A mixed-update pattern must contain six slots")
        mixed_position = {
            update: position for position, update in enumerate(mixed_updates)
        }
        for update in range(updates):
            if update in mixed_position:
                aggregate_index = mixed_position[update]
                celltype = non_neuronal_ids[aggregate_index]
                aggregate = [
                    (celltype, draw) for draw in range(aggregate_draws)
                ]
            else:
                aggregate = []
            neuronal_count = blocks_per_update - len(aggregate)
            update_neuronal = neuronal_stream[
                neuronal_cursor : neuronal_cursor + neuronal_count
            ]
            if len(update_neuronal) != neuronal_count:
                raise RuntimeError("v6A schedule exhausted neuronal blocks")
            neuronal_cursor += neuronal_count
            update_schedule: list[tuple[int, int] | None] = [None] * (
                blocks_per_update
            )
            if aggregate:
                # Rotate a contiguous four-slot window across the six mixed
                # updates.  This assigns exactly four non-neuronal draws to
                # every rank over an epoch at world sizes three and six.
                aggregate_start = mixed_position[update]
                for offset, item in enumerate(aggregate):
                    update_schedule[aggregate_start + offset] = item
            neuronal_iterator = iter(update_neuronal)
            for position, item in enumerate(update_schedule):
                if item is None:
                    update_schedule[position] = next(neuronal_iterator)
            try:
                next(neuronal_iterator)
            except StopIteration:
                pass
            else:
                raise RuntimeError("v6A update left neuronal blocks unused")
            schedule.extend(
                item for item in update_schedule if item is not None
            )
            neuronal_by_update.append(neuronal_count)
            non_neuronal_by_update.append(len(aggregate))
        if neuronal_cursor != len(neuronal_stream):
            raise RuntimeError("v6A schedule left neuronal blocks unused")

    occurrence: dict[int, list[int]] = {index: [] for index in range(len(names))}
    for celltype, draw in schedule:
        occurrence[celltype].append(draw)
    expected_draws = list(range(draws))
    bad_occurrence = {
        names[celltype]: values
        for celltype, values in occurrence.items()
        if values != expected_draws
    }
    if bad_occurrence:
        raise RuntimeError(
            "lineage schedule did not assign every draw exactly once: "
            f"{bad_occurrence}"
        )

    raw_non_neuronal_weight = len(neuronal_ids) / len(non_neuronal_ids)
    raw_weights = tuple(
        raw_non_neuronal_weight
        if index in set(non_neuronal_ids)
        else 1.0
        for index in range(len(names))
    )
    mean_raw_weight = sum(raw_weights) / len(raw_weights)
    normalized_weights = tuple(
        value / mean_raw_weight for value in raw_weights
    )
    epoch_weight_sum = (
        neuronal_blocks * normalized_weights[neuronal_ids[0]]
        + non_neuronal_blocks * normalized_weights[non_neuronal_ids[0]]
    )
    if not math.isclose(
        epoch_weight_sum,
        float(total_blocks),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("normalized lineage weights do not preserve epoch scale")
    if pseudobulk_mode == NON_NEURONAL_PSEUDOBULK_LEGACY:
        update_weight_sum = (
            neuronal_per_update * normalized_weights[neuronal_ids[0]]
            + non_neuronal_per_update
            * normalized_weights[non_neuronal_ids[0]]
        )
        if not math.isclose(
            update_weight_sum,
            float(blocks_per_update),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise RuntimeError(
                "normalized lineage weights do not preserve optimizer scale"
            )

    return LineageBalancePlan(
        mode=mode,
        celltype_raw_weights=raw_weights,
        celltype_normalized_weights=normalized_weights,
        neuronal_celltype_ids=neuronal_ids,
        non_neuronal_celltype_ids=non_neuronal_ids,
        global_schedule=tuple(schedule),
        neuronal_blocks_per_update=neuronal_per_update,
        non_neuronal_blocks_per_update=non_neuronal_per_update,
        neuronal_blocks_by_update=tuple(neuronal_by_update),
        non_neuronal_blocks_by_update=tuple(non_neuronal_by_update),
        non_neuronal_pseudobulk_mode=pseudobulk_mode,
        aggregate_draws=int(config.non_neuronal_aggregate_draws),
    )


MODULE_RESCUE_ARTIFACT_SCHEMA = "kmlee_bam.prism_module_rescue_stats.v2"


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _robust_scale(values: np.ndarray, axis: int = 0) -> np.ndarray:
    location = np.median(values, axis=axis)
    mad = np.median(
        np.abs(values - np.expand_dims(location, axis=axis)),
        axis=axis,
    )
    return 1.4826 * mad


@dataclass(frozen=True)
class ModuleRescueArtifact:
    group_module: np.ndarray
    group_count: np.ndarray
    celltype_scale: np.ndarray
    module_names: tuple[str, ...]
    celltype_names: tuple[str, ...]
    donor_names: tuple[str, ...]
    ccc_valid_mask: np.ndarray | None = None
    pathology_residual: np.ndarray | None = None
    pathology_residual_valid: np.ndarray | None = None
    pathology_reliable_mask: np.ndarray | None = None
    pathology_axis_names: tuple[str, ...] = ()
    schema_version: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    artifact_sha256: str = ""

    @classmethod
    def load_for_donors(
        cls,
        path: str | Path,
        *,
        donor_ids: Iterable[int],
        expected_donor_names: Iterable[str],
        expected_celltype_names: Iterable[str],
        expected_module_names: Iterable[str],
        expected_gene_names: Iterable[str] | None = None,
        expected_spec_path: str | Path | None = None,
        expected_registry_path: str | Path | None = None,
        expected_activity_path: str | Path | None = None,
        expected_activity_normalization: str | None = None,
        expected_artifact_sha256: str | None = None,
        minimum_target_cells: int,
        minimum_donors_per_celltype: int,
        scale_floor: float,
        scale_shrinkage: float,
        ccc_target_variance_floor: float,
        require_provenance: bool = False,
    ) -> "ModuleRescueArtifact":
        source = Path(path)
        artifact_sha256 = _sha256_file(source)
        if expected_artifact_sha256 is not None:
            expected_digest = str(expected_artifact_sha256).lower()
            if artifact_sha256 != expected_digest:
                raise ValueError(
                    "module-rescue artifact SHA256 mismatch: "
                    f"{artifact_sha256} != {expected_digest}"
                )
        with np.load(source, allow_pickle=False) as archive:
            required = {
                "group_module",
                "group_count",
                "module_names",
                "celltype_names",
                "donor_names",
            }
            missing = sorted(required.difference(archive.files))
            if missing:
                raise ValueError(
                    f"module-rescue artifact is missing arrays: {missing}"
                )
            group_module = np.asarray(
                archive["group_module"], dtype=np.float32
            )
            group_count = np.asarray(archive["group_count"], dtype=np.int64)
            module_names = tuple(
                np.asarray(archive["module_names"]).astype(str).tolist()
            )
            celltype_names = tuple(
                np.asarray(archive["celltype_names"]).astype(str).tolist()
            )
            donor_names = tuple(
                np.asarray(archive["donor_names"]).astype(str).tolist()
            )
            gene_names = (
                tuple(np.asarray(archive["gene_names"]).astype(str).tolist())
                if "gene_names" in archive.files
                else None
            )
            schema_version = (
                str(np.asarray(archive["schema_version"]).item())
                if "schema_version" in archive.files
                else ""
            )
            if "metadata_json" in archive.files:
                try:
                    metadata = json.loads(
                        str(np.asarray(archive["metadata_json"]).item())
                    )
                except Exception as exc:
                    raise ValueError(
                        "module-rescue metadata_json is invalid JSON"
                    ) from exc
            else:
                metadata = {}
            pathology_residual = (
                np.asarray(archive["pathology_residual"], dtype=np.float32)
                if "pathology_residual" in archive.files
                else None
            )
            pathology_residual_valid = (
                np.asarray(archive["pathology_residual_valid"], dtype=bool)
                if "pathology_residual_valid" in archive.files
                else None
            )
            pathology_reliable_mask = (
                np.asarray(archive["pathology_reliable_mask"], dtype=bool)
                if "pathology_reliable_mask" in archive.files
                else None
            )
            pathology_axis_names = (
                tuple(np.asarray(archive["pathology_axis_names"]).astype(str).tolist())
                if "pathology_axis_names" in archive.files
                else ()
            )

        expected_donor_names = tuple(str(x) for x in expected_donor_names)
        expected_celltype_names = tuple(
            str(x) for x in expected_celltype_names
        )
        expected_module_names = tuple(str(x) for x in expected_module_names)
        if donor_names != expected_donor_names:
            raise ValueError(
                "module-rescue donor order differs from the cell dataset"
            )
        if celltype_names != expected_celltype_names:
            raise ValueError(
                "module-rescue celltype order differs from the ordinal spec"
            )
        if module_names != expected_module_names:
            raise ValueError(
                "module-rescue module order differs from the registry"
            )
        if expected_gene_names is not None:
            expected_gene_names = tuple(str(x) for x in expected_gene_names)
            if gene_names != expected_gene_names:
                raise ValueError(
                    "module-rescue gene order differs from the ordinal spec"
                )
        if group_module.shape[:2] != group_count.shape:
            raise ValueError(
                "module-rescue group_module/group_count dimensions differ"
            )
        if pathology_residual is not None:
            expected_path_shape = (len(donor_names), len(pathology_axis_names))
            if pathology_residual.shape != expected_path_shape:
                raise ValueError("pathology_residual shape mismatch")
            if (
                pathology_residual_valid is None
                or pathology_residual_valid.shape != expected_path_shape
            ):
                raise ValueError("pathology_residual_valid shape mismatch")
            if (
                pathology_reliable_mask is None
                or pathology_reliable_mask.shape
                != (len(celltype_names), len(pathology_axis_names), len(module_names))
            ):
                raise ValueError("pathology_reliable_mask shape mismatch")
        if bool(require_provenance):
            if schema_version != MODULE_RESCUE_ARTIFACT_SCHEMA:
                raise ValueError(
                    "module-rescue artifact schema mismatch: "
                    f"{schema_version!r} != {MODULE_RESCUE_ARTIFACT_SCHEMA!r}"
                )
            if metadata.get("schema_version") != MODULE_RESCUE_ARTIFACT_SCHEMA:
                raise ValueError("module-rescue metadata schema mismatch")
            if metadata.get("source_split") != "train_only":
                raise ValueError(
                    "module-rescue artifact must declare source_split=train_only"
                )
            if expected_activity_normalization is None:
                raise ValueError(
                    "v2 provenance validation requires expected activity "
                    "normalization"
                )
            if metadata.get("activity_normalization") != str(
                expected_activity_normalization
            ):
                raise ValueError(
                    "module-rescue activity normalization mismatch: "
                    f"{metadata.get('activity_normalization')!r} != "
                    f"{str(expected_activity_normalization)!r}"
                )
            for label, expected_path, metadata_key in (
                ("spec", expected_spec_path, "spec_sha256"),
                ("registry", expected_registry_path, "registry_sha256"),
                ("activity", expected_activity_path, "activity_sha256"),
            ):
                if expected_path is None:
                    raise ValueError(
                        f"v2 provenance validation requires expected {label} path"
                    )
                expected_hash = _sha256_file(expected_path)
                if metadata.get(metadata_key) != expected_hash:
                    raise ValueError(
                        f"module-rescue {label} SHA256 mismatch"
                    )

        selected = np.asarray(sorted({int(x) for x in donor_ids}), dtype=np.int64)
        if selected.size == 0:
            raise ValueError("module-rescue requires weight-training donors")
        eligible = group_count[selected] >= int(minimum_target_cells)
        flat = group_module[selected][eligible]
        flat = flat[np.isfinite(flat).all(axis=1)]
        if flat.shape[0] < 2:
            raise ValueError("too few donor-context rows to fit module scale")
        pooled = np.maximum(
            _robust_scale(flat, axis=0),
            float(scale_floor),
        )

        n_celltype = group_module.shape[1]
        n_module = group_module.shape[2]
        scale = np.empty((n_celltype, n_module), dtype=np.float32)
        for celltype in range(n_celltype):
            keep = (
                group_count[selected, celltype] >= int(minimum_target_cells)
            )
            values = group_module[selected[keep], celltype]
            values = values[np.isfinite(values).all(axis=1)]
            if values.shape[0] < 2:
                scale[celltype] = pooled
                continue
            raw = _robust_scale(values, axis=0)
            shrunk = (
                (1.0 - float(scale_shrinkage)) * raw
                + float(scale_shrinkage) * pooled
            )
            scale[celltype] = np.maximum(
                np.where(np.isfinite(shrunk), shrunk, pooled),
                float(scale_floor),
            ).astype(np.float32)

        ccc_valid = np.zeros((n_celltype, n_module), dtype=bool)
        variance_floor = float(ccc_target_variance_floor)
        for celltype in range(n_celltype):
            keep = group_count[selected, celltype] >= int(minimum_target_cells)
            values = group_module[selected[keep], celltype]
            values = values[np.isfinite(values).all(axis=1)]
            if values.shape[0] < int(minimum_donors_per_celltype):
                continue
            centered = values.astype(np.float64) - values.astype(
                np.float64
            ).mean(axis=0, keepdims=True)
            standardized = centered / scale[celltype].astype(np.float64)
            variance = np.mean(np.square(standardized), axis=0)
            ccc_valid[celltype] = (
                np.isfinite(variance) & (variance >= variance_floor)
            )

        return cls(
            group_module=group_module,
            group_count=group_count,
            celltype_scale=scale,
            module_names=module_names,
            celltype_names=celltype_names,
            donor_names=donor_names,
            ccc_valid_mask=ccc_valid,
            pathology_residual=pathology_residual,
            pathology_residual_valid=pathology_residual_valid,
            pathology_reliable_mask=pathology_reliable_mask,
            pathology_axis_names=pathology_axis_names,
            schema_version=schema_version,
            metadata=metadata,
            artifact_sha256=artifact_sha256,
        )


class DonorCelltypeContrastSampler:
    """Deterministic complete donor blocks for one cell type at a time."""

    def __init__(
        self,
        dataset: Any,
        *,
        artifact: ModuleRescueArtifact,
        donor_ids: Iterable[int],
        cells_per_group: int,
        minimum_donors: int,
        minimum_target_cells: int,
        seed: int,
        nonoverlap_sampling: bool = False,
    ) -> None:
        self.dataset = dataset
        self.artifact = artifact
        self.donor_ids = tuple(sorted({int(x) for x in donor_ids}))
        self.cells_per_group = int(cells_per_group)
        self.minimum_donors = int(minimum_donors)
        self.minimum_target_cells = int(minimum_target_cells)
        self.seed = int(seed)
        self.nonoverlap_sampling = bool(nonoverlap_sampling)

        row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
        donor = np.asarray(dataset.donor_ids, dtype=np.int64)[row_idx]
        celltype = np.asarray(dataset.celltype_ids, dtype=np.int64)[row_idx]
        selected_set = set(self.donor_ids)
        self.groups: dict[tuple[int, int], np.ndarray] = {}
        for ct in range(len(dataset.spec.celltype_vocab)):
            for dn in self.donor_ids:
                if dn not in selected_set:
                    continue
                local = np.flatnonzero((donor == dn) & (celltype == ct))
                if local.size < self.cells_per_group:
                    continue
                if (
                    int(artifact.group_count[dn, ct])
                    < self.minimum_target_cells
                ):
                    continue
                self.groups[(ct, dn)] = local.astype(np.int64)
        self.eligible_celltypes = tuple(
            ct
            for ct in range(len(dataset.spec.celltype_vocab))
            if sum((ct, dn) in self.groups for dn in self.donor_ids)
            >= self.minimum_donors
        )
        if not self.eligible_celltypes:
            raise ValueError(
                "no cell type has enough weight-training donor groups for "
                "module rescue"
            )

    def sample(
        self,
        *,
        epoch: int,
        slot: int,
        rank: int,
        world_size: int,
    ) -> dict[str, Any]:
        global_block = int(slot) * int(world_size) + int(rank)
        if self.nonoverlap_sampling:
            # The global sequence is independent of how blocks are sharded
            # over 4 or 6 GPUs.  ``epoch`` changes the cell permutation seed,
            # not the schedule order or draw number.
            schedule_start = 0
            schedule_index = global_block
        else:
            # Historical schedule retained for legacy configurations.
            schedule_start = (int(epoch) - 1) * int(world_size)
            schedule_index = schedule_start + global_block
        celltype_position = schedule_index % len(self.eligible_celltypes)
        celltype = self.eligible_celltypes[celltype_position]
        first_block = (
            int(celltype_position) - schedule_start
        ) % len(self.eligible_celltypes)
        draw_index = (
            global_block - first_block
        ) // len(self.eligible_celltypes)
        return self.sample_celltype(
            celltype=celltype,
            epoch=epoch,
            slot=slot,
            rank=rank,
            draw_index=draw_index,
        )

    def _sample_candidate_indices(
        self,
        candidates: np.ndarray,
        *,
        epoch: int,
        celltype: int,
        donor: int,
        draw_index: int,
        fallback_seed: int,
    ) -> np.ndarray:
        """Choose one draw, exhausting a deterministic permutation before reuse."""

        if not self.nonoverlap_sampling:
            rng = np.random.default_rng(int(fallback_seed))
            return rng.choice(
                candidates,
                size=self.cells_per_group,
                replace=False,
            )

        count = int(candidates.size)
        start = int(draw_index) * self.cells_per_group
        stop = start + self.cells_per_group
        selected: list[int] = []
        cursor = start
        while cursor < stop:
            cycle = cursor // count
            offset = cursor % count
            cycle_seed = (
                self.seed
                + 1_000_003 * int(epoch)
                + 97 * int(celltype)
                + 53 * int(donor)
                + 15_485_863 * int(cycle)
            )
            permutation = np.random.default_rng(cycle_seed).permutation(
                candidates
            )
            take = min(stop - cursor, count - offset)
            selected.extend(
                int(value)
                for value in permutation[offset : offset + take].tolist()
            )
            cursor += take
        return np.asarray(selected, dtype=np.int64)

    def sample_celltype(
        self,
        *,
        celltype: int,
        epoch: int,
        slot: int,
        rank: int,
        draw_index: int | None = None,
    ) -> dict[str, Any]:
        """Sample one named eligible cell type with complete donor support."""

        celltype = int(celltype)
        if celltype not in self.eligible_celltypes:
            raise ValueError(f"celltype {celltype} is not rescue-eligible")
        selected_local: list[int] = []
        resolved_draw_index = int(slot) if draw_index is None else int(draw_index)
        if resolved_draw_index < 0:
            raise ValueError("draw_index must be non-negative")
        for donor in self.donor_ids:
            candidates = self.groups.get((celltype, donor))
            if candidates is None:
                continue
            group_seed = (
                self.seed
                + 1_000_003 * int(epoch)
                + 10_007 * int(slot)
                + 1_009 * int(rank)
                + 97 * int(celltype)
                + 53 * int(donor)
            )
            chosen = self._sample_candidate_indices(
                candidates,
                epoch=int(epoch),
                celltype=int(celltype),
                donor=int(donor),
                draw_index=resolved_draw_index,
                fallback_seed=group_seed,
            )
            selected_local.extend(int(value) for value in chosen.tolist())
        if len(selected_local) < self.minimum_donors * self.cells_per_group:
            raise RuntimeError(
                f"celltype {celltype} produced too few sampled cells"
            )
        batch = default_collate(
            [self.dataset[index] for index in selected_local]
        )
        batch["_module_rescue_celltype"] = int(celltype)
        batch["_module_rescue_draw_index"] = resolved_draw_index
        batch["_module_rescue_local_indices"] = torch.as_tensor(
            selected_local,
            dtype=torch.long,
        )
        return batch


def _unwrap_system(system: Any) -> Any:
    if hasattr(system, "module"):
        return system.module
    if hasattr(system, "_ddp_model") and hasattr(system._ddp_model, "module"):
        return system._ddp_model.module
    return system


def _manual_average_gradients(
    module: torch.nn.Module,
    *,
    world_size: int,
) -> None:
    """Average auxiliary-step gradients without arming DDP's reducer.

    The full system forward returns many diagnostic tensors, while the
    module-rescue loss intentionally consumes only two decoder views.  DDP's
    unused-output discovery consequently expects gradients for outputs that
    are irrelevant to this auxiliary objective and fails on the next forward.
    The auxiliary forward/backward therefore runs under ``DDP.no_sync()`` and
    this function performs the same cross-rank mean explicitly.

    A parameter may be used on only a subset of ranks because ranks sample
    different cell types.  The global-used bitmap makes those ranks contribute
    zeros, matching ordinary DDP semantics; parameters unused on every rank
    retain ``grad=None`` and are not spuriously weight-decayed.
    """

    if int(world_size) <= 1:
        return
    if not dist.is_available() or not dist.is_initialized():
        raise RuntimeError(
            "manual module-rescue gradient averaging requires an initialized "
            "distributed process group"
        )

    parameters = [
        parameter
        for parameter in module.parameters()
        if parameter.requires_grad
    ]
    if not parameters:
        return
    device = parameters[0].device
    local_used = torch.tensor(
        [parameter.grad is not None for parameter in parameters],
        dtype=torch.uint8,
        device=device,
    )
    global_used = local_used.clone()
    dist.all_reduce(global_used, op=dist.ReduceOp.MAX)

    buckets: dict[
        tuple[torch.device, torch.dtype],
        list[tuple[torch.nn.Parameter, torch.Tensor]],
    ] = {}
    for parameter, used in zip(parameters, global_used.tolist()):
        if not bool(used):
            continue
        gradient = parameter.grad
        if gradient is None:
            gradient = torch.zeros_like(
                parameter, memory_format=torch.preserve_format
            )
        elif gradient.is_sparse:
            gradient = gradient.to_dense()
        buckets.setdefault(
            (gradient.device, gradient.dtype), []
        ).append((parameter, gradient))

    for entries in buckets.values():
        flat = torch.cat(
            [gradient.detach().reshape(-1) for _, gradient in entries],
            dim=0,
        )
        dist.all_reduce(flat, op=dist.ReduceOp.SUM)
        flat.div_(float(world_size))
        offset = 0
        for parameter, gradient in entries:
            count = int(gradient.numel())
            averaged = flat[offset : offset + count].view_as(gradient)
            if parameter.grad is None or parameter.grad.is_sparse:
                parameter.grad = averaged.clone()
            else:
                parameter.grad.copy_(averaged)
            offset += count


def _module_rescue_parameter_allowlist(
    system: torch.nn.Module,
    *,
    include_state_readout: bool = False,
    include_module_local: bool = False,
) -> dict[int, str]:
    """Return parameters that an auxiliary rescue step may update.

    Ordinary end-to-end steps still update the complete model.  This allowlist
    applies only to the extra rescue optimizer step, preventing a donor-module
    auxiliary target from moving the encoder, nuisance projector, baselines,
    ordinal thresholds, or the shared base/tech/state mixer.
    """

    base = _unwrap_system(system)
    decoder = getattr(base, "decoder", None)
    precision = getattr(base, "precision_head", None)
    if decoder is None:
        raise RuntimeError("module-rescue parameter restriction requires decoder")
    if precision is None:
        raise RuntimeError(
            "module-rescue parameter restriction requires precision_head"
        )

    allowed: set[int] = set()

    def add_module(module: torch.nn.Module | None) -> None:
        if module is None:
            return
        for parameter in module.parameters():
            if parameter.requires_grad:
                allowed.add(id(parameter))

    def add_parameter(parameter: torch.nn.Parameter | None) -> None:
        if parameter is not None and parameter.requires_grad:
            allowed.add(id(parameter))

    add_module(getattr(decoder, "coeff_head", None))
    add_parameter(getattr(decoder, "generator_u", None))
    add_parameter(getattr(decoder, "generator_v", None))
    add_parameter(getattr(decoder, "generator_a", None))
    add_module(getattr(decoder, "direct_state_head", None))

    # In online joint-count mode, module rescue must also tell the architecture
    # gate which generator routes are needed for cell-type/module recovery.
    # The rescue forward uses a deterministic exact-hard mask with a
    # straight-through probability gradient, so allowing only ``log_alpha``
    # here preserves the v2 isolation contract while retaining that signal.
    # Legacy/gate-only runs do not attach this parameter to the allowlist.
    count_cfg = getattr(base, "generator_count_config", None)
    count_gate = getattr(base, "generator_count_gate", None)
    if (
        count_cfg is not None
        and bool(getattr(count_cfg, "enabled", False))
        and str(getattr(count_cfg, "mode", "gate_only")) == "joint"
    ):
        if count_gate is None:
            raise RuntimeError(
                "joint generator-count rescue requires generator_count_gate"
            )
        add_parameter(getattr(count_gate, "log_alpha", None))

    # Keep rescue updates on the explicit module-to-gene dictionaries.  Adding
    # the complete PrecisionMedicineHead here would also update its donor
    # context Transformer, personal posterior, regional/age baselines and
    # adversaries, contradicting the v2 isolation contract.  Interaction bases
    # are explicit output dictionaries too, so include them only when present.
    precision_output_parameters = (
        "common_global",
        "common_celltype_delta",
        "common_region_gate",
        "personal_basis",
        "response_basis",
        "interaction_global",
        "interaction_module_basis",
        "interaction_celltype_loading",
    )
    for name in precision_output_parameters:
        add_parameter(getattr(precision, name, None))

    if bool(include_module_local):
        if not bool(getattr(precision, "module_local_enabled", False)):
            raise RuntimeError(
                "module-local rescue updates require an enabled local path"
            )
        local_scale = getattr(precision, "module_local_scale", None)
        if local_scale is None or float(local_scale.detach().cpu()) < 1.0 - 1.0e-8:
            raise RuntimeError(
                "module-local rescue updates require the completed local ramp"
            )
        for name in (
            "module_local_read",
            "module_local_write_global",
            "module_local_write_celltype_delta",
            "module_local_diagonal_global",
            "module_local_diagonal_celltype_delta",
            "module_local_output_gate",
        ):
            parameter = getattr(precision, name, None)
            if parameter is None:
                raise RuntimeError(
                    f"module-local rescue parameter is missing: {name}"
                )
            add_parameter(parameter)
        if bool(getattr(precision, "module_local_nonlinear_enabled", False)):
            nonlinear_scale = getattr(
                precision, "module_local_nonlinear_scale", None
            )
            if (
                nonlinear_scale is None
                or float(nonlinear_scale.detach().cpu()) < 1.0 - 1.0e-8
            ):
                raise RuntimeError(
                    "nonlinear module-local rescue requires its completed ramp"
                )
            trainable_families = getattr(
                precision,
                "module_local_nonlinear_trainable_families",
                ("mix", "threshold", "slope", "gain"),
            )
            for family in trainable_families:
                for suffix in ("global", "celltype_delta"):
                    name = f"module_local_nonlinear_{family}_{suffix}"
                    parameter = getattr(precision, name, None)
                    if parameter is None:
                        raise RuntimeError(
                            "nonlinear module-local rescue parameter is missing: "
                            f"{name}"
                        )
                    add_parameter(parameter)

    if bool(include_state_readout):
        state_encoder = getattr(base, "state_encoder", None)
        if state_encoder is None:
            raise RuntimeError(
                "state-readout rescue updates require state_encoder"
            )
        persistent_pool = getattr(
            state_encoder, "persistent_attention_pool", None
        )
        if persistent_pool is None:
            raise RuntimeError(
                "state-readout rescue updates require persistent AGP pooling"
            )
        add_module(persistent_pool)
        add_module(getattr(state_encoder, "posterior_input_norm", None))
        add_module(getattr(state_encoder, "posterior_head", None))

    # ScoreResidualMixer has one softmax that jointly changes base, technical,
    # and state gates.  There is no state-only parameter subset, so it is
    # intentionally frozen in rescue-only steps.
    named = {
        id(parameter): name
        for name, parameter in base.named_parameters()
        if parameter.requires_grad
    }
    missing = sorted(parameter_id for parameter_id in allowed if parameter_id not in named)
    if missing:
        raise RuntimeError("module-rescue allowlist contains unregistered parameters")
    selected = {
        parameter_id: named[parameter_id]
        for parameter_id in allowed
    }
    if not selected:
        raise RuntimeError("module-rescue parameter allowlist is empty")
    return selected


def _apply_module_rescue_gradient_allowlist(
    system: torch.nn.Module,
    allowed: dict[int, str],
) -> tuple[int, int]:
    """Clear disallowed gradients and return used-allowed/cleared counts."""

    used_allowed = 0
    cleared = 0
    for parameter in _unwrap_system(system).parameters():
        if parameter.grad is None:
            continue
        if id(parameter) in allowed:
            used_allowed += 1
        else:
            parameter.grad = None
            cleared += 1
    if used_allowed <= 0:
        raise RuntimeError(
            "module-rescue loss produced no gradient for any allowed parameter"
        )
    return used_allowed, cleared


@contextmanager
def _temporarily_freeze_outside_allowlist(
    system: torch.nn.Module,
    allowed: dict[int, str] | None,
):
    """Avoid building rescue graphs for parameters cleared before the step."""

    if allowed is None:
        yield
        return
    snapshot = [
        (parameter, bool(parameter.requires_grad))
        for parameter in _unwrap_system(system).parameters()
    ]
    try:
        for parameter, required in snapshot:
            if required and id(parameter) not in allowed:
                parameter.requires_grad_(False)
        yield
    finally:
        for parameter, required in snapshot:
            parameter.requires_grad_(required)


@dataclass(frozen=True)
class ExactGroupMeanVJP:
    """Small first-pass loss graph and its derivatives with respect to means."""

    output: PrismModuleRescueOutput
    branch_gradient: torch.Tensor
    full_gradient: torch.Tensor


def _exact_group_mean_leaf_vjp(
    *,
    branch_group_mean: torch.Tensor,
    full_group_mean: torch.Tensor,
    loss_builder: Callable[
        [torch.Tensor, torch.Tensor], PrismModuleRescueOutput
    ],
    loss_scale: float,
) -> ExactGroupMeanVJP:
    """Differentiate only the nonlinear loss over detached group means."""

    branch_leaf = branch_group_mean.detach().requires_grad_(True)
    full_leaf = full_group_mean.detach().requires_grad_(True)
    output = loss_builder(branch_leaf, full_leaf)
    scaled_loss = output.loss * float(loss_scale)
    branch_gradient, full_gradient = torch.autograd.grad(
        scaled_loss,
        (branch_leaf, full_leaf),
        retain_graph=False,
        create_graph=False,
    )
    return ExactGroupMeanVJP(
        output=output,
        branch_gradient=branch_gradient.detach(),
        full_gradient=full_gradient.detach(),
    )


def _grouped_mean_replay_proxy(
    *,
    branch_chunk: GroupedModuleActivity,
    full_chunk: GroupedModuleActivity,
    aggregate_group_id: torch.Tensor,
    aggregate_cell_count: torch.Tensor,
    branch_gradient: torch.Tensor,
    full_gradient: torch.Tensor,
) -> torch.Tensor:
    """Return the exact microbatch VJP proxy for one aggregate mean."""

    if not bool(torch.equal(branch_chunk.group_id, aggregate_group_id)):
        raise RuntimeError("branch replay group order differs from first pass")
    if not bool(torch.equal(full_chunk.group_id, aggregate_group_id)):
        raise RuntimeError("full replay group order differs from first pass")
    if not bool(torch.equal(branch_chunk.cell_count, full_chunk.cell_count)):
        raise RuntimeError("branch/full replay group counts differ")
    denominator = aggregate_cell_count.to(
        device=branch_chunk.mean.device,
        dtype=branch_chunk.mean.dtype,
    ).clamp_min(1.0)
    fraction = branch_chunk.cell_count.to(
        device=branch_chunk.mean.device,
        dtype=branch_chunk.mean.dtype,
    ) / denominator
    return (
        branch_chunk.mean * branch_gradient.to(branch_chunk.mean.dtype)
        + full_chunk.mean * full_gradient.to(full_chunk.mean.dtype)
    ).mul(fraction[:, None]).sum()


class PrismModuleRescueUpdater:
    """Run a fixed number of grouped module-rescue optimiser updates."""

    def __init__(
        self,
        *,
        config: PrismModuleRescueTrainingConfig,
        dataset: Any,
        donor_ids: Iterable[int],
        module_names: Iterable[str],
        spec_path: str | Path | None = None,
        registry_path: str | Path | None = None,
        activity_path: str | Path | None = None,
        activity_normalization: str | None = None,
    ) -> None:
        config.validate()
        if not bool(config.enabled):
            raise ValueError("cannot construct a disabled module-rescue updater")
        self.config = config
        self.donor_ids = tuple(sorted({int(x) for x in donor_ids}))
        self.artifact = ModuleRescueArtifact.load_for_donors(
            str(config.stats_path),
            donor_ids=self.donor_ids,
            expected_donor_names=np.asarray(dataset.donor_vocab).astype(str),
            expected_celltype_names=np.asarray(
                dataset.spec.celltype_vocab
            ).astype(str),
            expected_module_names=module_names,
            expected_gene_names=np.asarray(dataset.spec.gene_names).astype(str),
            expected_spec_path=spec_path,
            expected_registry_path=registry_path,
            expected_activity_path=activity_path,
            expected_activity_normalization=activity_normalization,
            expected_artifact_sha256=config.expected_stats_sha256,
            minimum_target_cells=int(config.minimum_target_cells),
            minimum_donors_per_celltype=int(
                config.minimum_donors_per_celltype
            ),
            scale_floor=float(config.scale_floor),
            scale_shrinkage=float(config.scale_shrinkage),
            ccc_target_variance_floor=float(
                config.ccc_target_variance_floor
            ),
            require_provenance=str(config.contract_version) == "v2",
        )
        if bool(config.pathology_axis_enabled):
            if (
                self.artifact.pathology_residual is None
                or self.artifact.pathology_residual_valid is None
                or self.artifact.pathology_reliable_mask is None
                or not self.artifact.pathology_axis_names
            ):
                raise ValueError(
                    "pathology-axis rescue requires train-only residual and reliability arrays"
                )
        n_celltypes = len(self.artifact.celltype_names)
        self.celltype_deficit_ema = np.zeros(n_celltypes, dtype=np.float64)
        self.celltype_dynamic_weight = np.ones(n_celltypes, dtype=np.float64)
        self._last_celltype_recovery = np.full(n_celltypes, np.nan, dtype=np.float64)
        if bool(config.celltype_tail_enabled):
            missing_baselines = [
                name
                for name in self.artifact.celltype_names
                if name not in config.ad_recovery_baseline_by_celltype
            ]
            if missing_baselines:
                raise ValueError(
                    "celltype-tail protection lacks frozen baseline values: "
                    + ", ".join(missing_baselines)
                )
        # A declared train-only artifact is not enough: compare every
        # donor×celltype cell count against the local training dataset.  This
        # fail-closed check catches same-shape artifacts contaminated by
        # validation/test rows or built from a different train split.
        row_idx = np.asarray(dataset.row_idx, dtype=np.int64)
        donor = np.asarray(dataset.donor_ids, dtype=np.int64)[row_idx]
        celltype = np.asarray(dataset.celltype_ids, dtype=np.int64)[row_idx]
        local_count = np.zeros_like(self.artifact.group_count, dtype=np.int64)
        np.add.at(local_count, (donor, celltype), 1)
        if not np.array_equal(local_count, self.artifact.group_count):
            mismatch = int(np.count_nonzero(local_count != self.artifact.group_count))
            raise ValueError(
                "module-rescue artifact group_count differs from the local "
                f"training split in {mismatch} donor×celltype rows"
            )
        selected = np.asarray(self.donor_ids, dtype=np.int64)
        eligible = (
            self.artifact.group_count[selected]
            >= int(config.minimum_target_cells)
        )
        target = self.artifact.group_module[selected]
        if not np.isfinite(target[eligible]).all():
            raise ValueError(
                "module-rescue artifact has non-finite eligible train targets"
            )
        self.sampler = DonorCelltypeContrastSampler(
            dataset,
            artifact=self.artifact,
            donor_ids=self.donor_ids,
            cells_per_group=int(config.cells_per_donor_celltype),
            minimum_donors=int(config.minimum_donors_per_celltype),
            minimum_target_cells=int(config.minimum_target_cells),
            seed=int(config.seed),
            nonoverlap_sampling=bool(config.nonoverlap_sampling),
        )
        candidate_sizes = np.asarray(
            [
                len(indices)
                for (celltype, _), indices in self.sampler.groups.items()
                if celltype in self.sampler.eligible_celltypes
            ],
            dtype=np.int64,
        )
        self.candidate_size_min = int(candidate_sizes.min())
        self.candidate_size_max = int(candidate_sizes.max())
        required_unique = (
            int(config.draws_per_celltype_per_epoch)
            * int(config.cells_per_donor_celltype)
        )
        if (
            str(config.contract_version) == "v2"
            and self.candidate_size_min < required_unique
        ):
            raise ValueError(
                "module-rescue v2 non-overlap contract requires at least "
                f"{required_unique} cells per eligible donor×celltype; "
                f"minimum={self.candidate_size_min}"
            )
        self.lineage_balance_plan = build_lineage_balance_plan(
            config=config,
            celltype_names=self.artifact.celltype_names,
            eligible_celltypes=self.sampler.eligible_celltypes,
        )
        self._parameter_allowlist: dict[int, str] | None = None
        self._base_parameter_allowlist: dict[int, str] | None = None
        self._state_parameter_allowlist: dict[int, str] | None = None
        self._base_module_local_parameter_allowlist: dict[int, str] | None = None
        self._state_module_local_parameter_allowlist: dict[int, str] | None = None

    def state_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "kmlee_bam.module_rescue_ad_tail.v1",
            "celltype_deficit_ema": self.celltype_deficit_ema.tolist(),
            "celltype_dynamic_weight": self.celltype_dynamic_weight.tolist(),
            "last_celltype_recovery": self._last_celltype_recovery.tolist(),
        }

    def load_state_dict(self, state: dict[str, Any]) -> None:
        if not state:
            return
        if state.get("schema_version") != "kmlee_bam.module_rescue_ad_tail.v1":
            raise ValueError("module-rescue tail state schema mismatch")
        n = len(self.artifact.celltype_names)
        for key, target in (
            ("celltype_deficit_ema", self.celltype_deficit_ema),
            ("celltype_dynamic_weight", self.celltype_dynamic_weight),
            ("last_celltype_recovery", self._last_celltype_recovery),
        ):
            value = np.asarray(state[key], dtype=np.float64)
            if value.shape != (n,):
                raise ValueError(f"module-rescue restored {key} shape mismatch")
            target[:] = value

    def _combined_celltype_weight(self, celltype: int) -> float:
        base = float(
            self.lineage_balance_plan.celltype_normalized_weights[int(celltype)]
        )
        return base * float(self.celltype_dynamic_weight[int(celltype)])

    def _update_tail_weights(self, recovery: np.ndarray, *, epoch: int) -> None:
        if not bool(self.config.celltype_tail_enabled) or int(epoch) < int(
            self.config.celltype_tail_start_epoch
        ):
            return
        finite = np.isfinite(recovery)
        if not bool(finite.any()):
            return
        baseline = np.asarray(
            [
                float(self.config.ad_recovery_baseline_by_celltype[name])
                for name in self.artifact.celltype_names
            ],
            dtype=np.float64,
        )
        deficit = np.maximum(baseline - recovery, 0.0)
        deficit[~finite] = self.celltype_deficit_ema[~finite]
        decay = float(self.config.celltype_deficit_ema_decay)
        self.celltype_deficit_ema = (
            decay * self.celltype_deficit_ema + (1.0 - decay) * deficit
        )
        raw = np.ones_like(self.celltype_deficit_ema)
        tail_n = max(
            1,
            int(math.ceil(len(raw) * float(self.config.celltype_tail_fraction))),
        )
        tail = np.argsort(self.celltype_deficit_ema)[-tail_n:]
        positive = self.celltype_deficit_ema[self.celltype_deficit_ema > 0]
        scale = float(np.median(positive)) if positive.size else 1.0
        raw[tail] += self.celltype_deficit_ema[tail] / max(scale, 1.0e-8)
        raw = np.clip(
            raw,
            float(self.config.celltype_weight_min),
            float(self.config.celltype_weight_max),
        )
        # Preserve exactly half the total objective mass in each lineage.
        for ids in (
            self.lineage_balance_plan.neuronal_celltype_ids,
            self.lineage_balance_plan.non_neuronal_celltype_ids,
        ):
            index = np.asarray(ids, dtype=np.int64)
            raw[index] /= max(float(raw[index].mean()), 1.0e-8)
        self.celltype_dynamic_weight[:] = raw
        self._last_celltype_recovery[:] = recovery

    def _add_pathology_axis_loss(
        self,
        output: Any,
        *,
        group_id: torch.Tensor,
        celltype: int,
        epoch: int,
    ) -> Any:
        """Add train-only named-axis coefficient recovery to one full donor block."""

        if (
            not bool(self.config.pathology_axis_enabled)
            or int(epoch) < int(self.config.pathology_axis_start_epoch)
            or output.branch_activity is None
            or output.full_activity is None
            or output.target_activity is None
        ):
            return output
        donor = torch.div(
            group_id.long(),
            len(self.artifact.celltype_names),
            rounding_mode="floor",
        )
        residual = torch.as_tensor(
            self.artifact.pathology_residual,
            device=group_id.device,
            dtype=torch.float32,
        ).index_select(0, donor)
        residual_valid = torch.as_tensor(
            self.artifact.pathology_residual_valid,
            device=group_id.device,
            dtype=torch.bool,
        ).index_select(0, donor)
        reliable = torch.as_tensor(
            self.artifact.pathology_reliable_mask[int(celltype)],
            device=group_id.device,
            dtype=torch.bool,
        )
        branch_losses: list[torch.Tensor] = []
        full_losses: list[torch.Tensor] = []
        correlations: list[torch.Tensor] = []
        magnitude_weight = float(self.config.pathology_axis_magnitude_weight)
        for axis in range(int(residual.shape[1])):
            rows = residual_valid[:, axis] & torch.isfinite(residual[:, axis])
            modules = reliable[axis]
            if int(rows.sum()) < int(self.config.minimum_donors_per_celltype):
                continue
            if int(modules.sum()) < 8:
                continue
            x = residual[rows, axis]
            x = x - x.mean()
            x = x / x.square().mean().sqrt().clamp_min(1.0e-6)
            denominator = x.square().sum().clamp_min(1.0e-6)

            def coefficient(values: torch.Tensor) -> torch.Tensor:
                return (x[:, None] * values[rows]).sum(dim=0) / denominator

            target = coefficient(output.target_activity).detach()[modules]
            target_center = target - target.mean()
            target_scale = target_center.square().mean().sqrt().clamp_min(1.0e-6)

            def one_view(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
                pred = coefficient(values)[modules]
                pred_center = pred - pred.mean()
                corr = (
                    (pred_center * target_center).sum()
                    / (
                        pred_center.square().sum().sqrt()
                        * target_center.square().sum().sqrt()
                    ).clamp_min(1.0e-6)
                ).clamp(-1.0, 1.0)
                magnitude = F.smooth_l1_loss(
                    pred / target_scale,
                    target / target_scale,
                )
                loss = (1.0 - corr) + magnitude_weight * magnitude
                return loss, corr

            branch_loss, branch_corr = one_view(output.branch_activity)
            full_loss, full_corr = one_view(output.full_activity)
            branch_losses.append(branch_loss)
            full_losses.append(full_loss)
            correlations.extend((branch_corr, full_corr))
        if not branch_losses:
            return output
        branch_loss = torch.stack(branch_losses).mean()
        full_loss = torch.stack(full_losses).mean()
        rho = float(self.config.branch_fraction)
        unweighted = rho * branch_loss + (1.0 - rho) * full_loss
        ramp = min(
            1.0,
            (int(epoch) - int(self.config.pathology_axis_start_epoch) + 1) / 4.0,
        )
        weighted = (
            float(self.config.lambda_pathology_axis) * float(ramp) * unweighted
        )
        return replace(
            output,
            loss=output.loss + weighted,
            pathology_axis_loss=weighted,
            pathology_axis_branch_loss=branch_loss,
            pathology_axis_full_loss=full_loss,
            pathology_axis_mean_correlation=torch.stack(correlations).mean().detach(),
        )

    def _sample_rescue_block(
        self,
        *,
        epoch: int,
        slot: int,
        rank: int,
        world_size: int,
    ) -> dict[str, Any]:
        """Sample the legacy order or the audited global lineage schedule."""

        plan = self.lineage_balance_plan
        if not plan.enabled:
            return self.sampler.sample(
                epoch=int(epoch),
                slot=int(slot),
                rank=int(rank),
                world_size=int(world_size),
            )
        global_block = int(slot) * int(world_size) + int(rank)
        if not 0 <= global_block < len(plan.global_schedule):
            raise RuntimeError(
                "lineage-balanced rescue global block is out of range: "
                f"{global_block} / {len(plan.global_schedule)}"
            )
        celltype, draw_index = plan.global_schedule[global_block]
        return self.sampler.sample_celltype(
            celltype=int(celltype),
            epoch=int(epoch),
            slot=int(slot),
            rank=int(rank),
            draw_index=int(draw_index),
        )

    def lineage_balance_manifest(self) -> dict[str, Any]:
        """Return an auditable, JSON-safe description of loss and schedule."""

        plan = self.lineage_balance_plan
        names = self.artifact.celltype_names
        schedule_payload = json.dumps(
            [list(value) for value in plan.global_schedule],
            separators=(",", ":"),
        ).encode("utf-8")
        return {
            "mode": plan.mode,
            "enabled": bool(plan.enabled),
            "canonical_non_neuronal_celltype_names": list(
                PRISM_NON_NEURONAL_CELLTYPE_NAMES
            ),
            "canonical_neuronal_celltype_names": list(
                PRISM_NEURONAL_CELLTYPE_NAMES
            ),
            "neuronal_celltype_names": [
                names[index] for index in plan.neuronal_celltype_ids
            ],
            "non_neuronal_celltype_names": [
                names[index] for index in plan.non_neuronal_celltype_ids
            ],
            "raw_celltype_weights": {
                name: float(plan.celltype_raw_weights[index])
                for index, name in enumerate(names)
            },
            "normalized_celltype_weights": {
                name: float(plan.celltype_normalized_weights[index])
                for index, name in enumerate(names)
            },
            "neuronal_blocks_per_optimizer_update": int(
                plan.neuronal_blocks_per_update
            ),
            "non_neuronal_blocks_per_optimizer_update": int(
                plan.non_neuronal_blocks_per_update
            ),
            "neuronal_block_equivalents_by_optimizer_update": list(
                plan.neuronal_blocks_by_update
            ),
            "non_neuronal_block_equivalents_by_optimizer_update": list(
                plan.non_neuronal_blocks_by_update
            ),
            "non_neuronal_pseudobulk_mode": (
                plan.non_neuronal_pseudobulk_mode
            ),
            "non_neuronal_aggregate_draws": int(plan.aggregate_draws),
            "non_neuronal_prediction_cells": int(
                self.config.cells_per_donor_celltype
                * (
                    plan.aggregate_draws
                    if plan.non_neuronal_pseudobulk_mode
                    != NON_NEURONAL_PSEUDOBULK_LEGACY
                    else 1
                )
            ),
            "neuronal_prediction_cells": int(
                self.config.cells_per_donor_celltype
            ),
            "two_pass_exact_group_mean_vjp": bool(
                plan.non_neuronal_pseudobulk_mode
                != NON_NEURONAL_PSEUDOBULK_LEGACY
            ),
            "global_schedule_block_count": len(plan.global_schedule),
            "global_schedule_sha256": hashlib.sha256(
                schedule_payload
            ).hexdigest(),
        }

    def steps_per_rank(self, *, world_size: int) -> int:
        """Resolve a GPU-count-independent global celltype exposure budget."""

        world_size = int(world_size)
        if world_size <= 0:
            raise ValueError("world_size must be positive")
        draws = int(self.config.draws_per_celltype_per_epoch)
        if draws <= 0:
            return int(self.config.steps_per_epoch)
        global_blocks = len(self.sampler.eligible_celltypes) * draws
        if global_blocks % world_size != 0:
            raise ValueError(
                "eligible_celltypes * draws_per_celltype_per_epoch must be "
                "divisible by world_size; got "
                f"{len(self.sampler.eligible_celltypes)} * {draws} / {world_size}."
            )
        resolved = global_blocks // world_size
        configured = int(self.config.steps_per_epoch)
        if configured != resolved:
            raise ValueError(
                "steps_per_epoch conflicts with the GPU-independent draw "
                f"contract: configured={configured}, required={resolved}."
            )
        return resolved

    def local_blocks_per_optimizer_update(self, *, world_size: int) -> int:
        """Number of local contrast blocks accumulated before one update."""

        world_size = int(world_size)
        global_per_update = int(
            self.config.global_blocks_per_optimizer_step
        )
        if global_per_update <= 0:
            return 1
        if global_per_update % world_size != 0:
            raise ValueError(
                "global_blocks_per_optimizer_step must be divisible by "
                f"world_size; got {global_per_update} / {world_size}."
            )
        local = global_per_update // world_size
        local_steps = self.steps_per_rank(world_size=world_size)
        if local_steps % local != 0:
            raise ValueError(
                "local rescue blocks per epoch must be divisible by local "
                f"accumulation: {local_steps} / {local}."
            )
        return local

    def optimizer_updates_per_epoch(self, *, world_size: int) -> int:
        local_steps = self.steps_per_rank(world_size=int(world_size))
        local_per_update = self.local_blocks_per_optimizer_update(
            world_size=int(world_size)
        )
        return local_steps // local_per_update

    def resolve_parameter_allowlist(
        self,
        system: torch.nn.Module,
        *,
        epoch: int | None = None,
    ) -> dict[int, str] | None:
        if not bool(self.config.restrict_parameter_updates):
            self._parameter_allowlist = None
            return None
        include_state = (
            bool(self.config.allow_state_readout_updates)
            and epoch is not None
            and int(epoch) >= int(self.config.state_readout_start_epoch)
        )
        base = _unwrap_system(system)
        precision = getattr(base, "precision_head", None)
        precision_cfg = getattr(precision, "cfg", None)
        include_module_local = bool(
            precision is not None
            and bool(getattr(precision_cfg, "module_local_enabled", False))
            and epoch is not None
            and int(epoch)
            >= int(getattr(precision_cfg, "module_local_rescue_start_epoch", 25))
        )
        if include_state and include_module_local:
            if self._state_module_local_parameter_allowlist is None:
                self._state_module_local_parameter_allowlist = (
                    _module_rescue_parameter_allowlist(
                        base,
                        include_state_readout=True,
                        include_module_local=True,
                    )
                )
            self._parameter_allowlist = (
                self._state_module_local_parameter_allowlist
            )
        elif include_state:
            if self._state_parameter_allowlist is None:
                self._state_parameter_allowlist = (
                    _module_rescue_parameter_allowlist(
                        base,
                        include_state_readout=True,
                        include_module_local=False,
                    )
                )
            self._parameter_allowlist = self._state_parameter_allowlist
        elif include_module_local:
            if self._base_module_local_parameter_allowlist is None:
                self._base_module_local_parameter_allowlist = (
                    _module_rescue_parameter_allowlist(
                        base,
                        include_state_readout=False,
                        include_module_local=True,
                    )
                )
            self._parameter_allowlist = (
                self._base_module_local_parameter_allowlist
            )
        else:
            if self._base_parameter_allowlist is None:
                self._base_parameter_allowlist = (
                    _module_rescue_parameter_allowlist(
                        base,
                        include_state_readout=False,
                        include_module_local=False,
                    )
                )
            self._parameter_allowlist = self._base_parameter_allowlist
        return self._parameter_allowlist

    def _branch_probabilities(
        self,
        trainer: Any,
        batch: dict[str, torch.Tensor],
        model_out: Any,
    ) -> torch.Tensor:
        system = _unwrap_system(trainer.system)
        decoder = system.decoder
        precision = getattr(model_out, "precision_out", None)
        if precision is None or precision.gene_score is None:
            raise RuntimeError(
                "module-rescue branch view requires enabled PRISM precision head"
            )
        tech_id = batch.get("tech_id", batch.get("batch_id"))
        if tech_id is None:
            raise KeyError("module-rescue batch lacks tech/batch id")
        score = (
            decoder.celltype_baseline(batch["celltype_id"]).detach()
            + decoder.tech_baseline(tech_id).detach()
            + precision.gene_score
        )
        if decoder.sex_baseline is not None and "sex_id" in batch:
            score = score + decoder.sex_baseline(
                decoder._sex_index(batch["sex_id"])
            ).detach()
        _, probabilities = decoder._score_to_probs(
            score, decoder._compute_thresholds()
        )
        return probabilities

    def calibrate_ccc_gradient_weight(
        self,
        trainer: Any,
        *,
        calibration_epoch: int = 1,
        calibration_blocks: int | None = None,
        clip_min: float = 0.1,
        clip_max: float = 1.0,
    ) -> dict[str, Any]:
        """Calibrate the CCC coefficient on deterministic train-only blocks.

        The calibration never takes an optimiser step.  It measures the L2
        gradient norm of the unweighted Huber and flattened-CCC rescue terms
        over one deterministic draw of every eligible cell type, then chooses
        ``lambda_ccc = clip(huber_grad_rms / ccc_grad_rms, 0.1, 1.0)``.
        Validation and test rows are unavailable to this updater by contract.
        """

        if int(calibration_epoch) <= 0:
            raise ValueError("calibration_epoch must be positive")
        if not 0.0 < float(clip_min) <= float(clip_max):
            raise ValueError("CCC calibration clip bounds are invalid")
        eligible = tuple(int(x) for x in self.sampler.eligible_celltypes)
        n_blocks = len(eligible) if calibration_blocks is None else int(
            calibration_blocks
        )
        if n_blocks <= 0 or n_blocks > len(eligible):
            raise ValueError(
                "calibration_blocks must lie in [1, number of eligible cell types]"
            )

        base_system = _unwrap_system(trainer.system)
        allowed = self.resolve_parameter_allowlist(base_system)
        if not allowed:
            raise RuntimeError(
                "CCC calibration requires the restricted module-rescue allowlist"
            )
        named_parameters = {
            id(parameter): (name, parameter)
            for name, parameter in base_system.named_parameters()
            if parameter.requires_grad
        }
        selected = [
            named_parameters[parameter_id]
            for parameter_id in sorted(allowed, key=lambda value: allowed[value])
        ]
        if not selected:
            raise RuntimeError("CCC calibration parameter list is empty")
        parameter_names = [name for name, _ in selected]
        parameters = [parameter for _, parameter in selected]

        activity_weight = base_system.module_tokenizer.activity_weight
        scale = torch.from_numpy(self.artifact.celltype_scale).to(
            device=trainer.device,
            dtype=torch.float32,
        )
        ccc_valid_mask = torch.from_numpy(
            self.artifact.ccc_valid_mask
        ).to(device=trainer.device, dtype=torch.bool)
        objective_config = PrismModuleRescueConfig(
            enabled=True,
            lambda_module=1.0,
            branch_fraction=float(self.config.branch_fraction),
            huber_delta=float(self.config.huber_delta),
            ccc_weight=0.0,
        )

        def gradient_norm(
            loss: torch.Tensor,
            *,
            retain_graph: bool,
        ) -> float:
            gradients = torch.autograd.grad(
                loss,
                parameters,
                retain_graph=retain_graph,
                allow_unused=True,
            )
            squared = torch.zeros((), device=trainer.device, dtype=torch.float64)
            used = 0
            for gradient in gradients:
                if gradient is None:
                    continue
                used += 1
                squared = squared + gradient.detach().double().square().sum()
            if used <= 0:
                raise RuntimeError("CCC calibration loss produced no allowed gradient")
            value = float(torch.sqrt(squared).cpu())
            if not math.isfinite(value):
                raise RuntimeError("CCC calibration produced a non-finite gradient norm")
            return value

        mode_snapshot = [
            (module, bool(module.training))
            for module in trainer.system.modules()
        ]
        rows: list[dict[str, Any]] = []
        trainer.system.eval()
        trainer.optimizer.zero_grad(set_to_none=True)
        try:
            for slot in range(n_blocks):
                raw_batch = self._sample_rescue_block(
                    epoch=int(calibration_epoch),
                    slot=int(slot),
                    rank=0,
                    world_size=1,
                )
                celltype = int(raw_batch.pop("_module_rescue_celltype"))
                draw_index = int(raw_batch.pop("_module_rescue_draw_index"))
                sampled_local_indices = raw_batch.pop(
                    "_module_rescue_local_indices"
                )
                batch = trainer._move_batch_to_device(raw_batch)
                batch["donor_balance_weight"] = torch.ones(
                    batch["donor_id"].shape,
                    dtype=torch.float32,
                    device=trainer.device,
                )
                with trainer._autocast_context():
                    model_out = trainer.system(
                        batch,
                        sample_latent=False,
                        return_all_hidden_states=False,
                        return_attn_diagnostics=False,
                    )
                    if model_out.decoder_out is None:
                        raise RuntimeError(
                            "CCC calibration requires decoder output"
                        )
                    branch_probabilities = self._branch_probabilities(
                        trainer, batch, model_out
                    )
                    donor = batch["donor_id"].long()
                    group_id = (
                        donor * len(self.artifact.celltype_names)
                        + batch["celltype_id"].long()
                    )
                    unique_group = torch.unique(group_id, sorted=True)
                    target_module = []
                    target_count = []
                    for group in unique_group.detach().cpu().tolist():
                        donor_id = int(group) // len(
                            self.artifact.celltype_names
                        )
                        ct_id = int(group) % len(
                            self.artifact.celltype_names
                        )
                        target_module.append(
                            self.artifact.group_module[donor_id, ct_id]
                        )
                        target_count.append(
                            self.artifact.group_count[donor_id, ct_id]
                        )
                    rescue = dual_view_module_rescue_loss(
                        config=objective_config,
                        branch_probabilities=branch_probabilities,
                        full_probabilities=model_out.decoder_out.probs,
                        target_group_activity=torch.as_tensor(
                            np.asarray(target_module, dtype=np.float32),
                            device=trainer.device,
                        ),
                        target_group_id=unique_group,
                        target_group_cell_count=torch.as_tensor(
                            np.asarray(target_count, dtype=np.int64),
                            device=trainer.device,
                        ),
                        activity_weight=activity_weight,
                        group_id=group_id,
                        celltype_id=batch["celltype_id"],
                        train_scale=scale,
                        ccc_valid_mask=ccc_valid_mask,
                        min_cells_per_group=1,
                        min_donor_groups_per_celltype=int(
                            self.config.minimum_donors_per_celltype
                        ),
                    )
                celltype_weight = float(
                    self.lineage_balance_plan.celltype_normalized_weights[
                        celltype
                    ]
                )
                weighted_huber_loss = rescue.huber_rescue_loss
                weighted_ccc_loss = rescue.ccc_rescue_loss
                if self.lineage_balance_plan.enabled:
                    weighted_huber_loss = (
                        weighted_huber_loss * celltype_weight
                    )
                    weighted_ccc_loss = weighted_ccc_loss * celltype_weight
                huber_norm = gradient_norm(
                    weighted_huber_loss,
                    retain_graph=True,
                )
                ccc_norm = gradient_norm(
                    weighted_ccc_loss,
                    retain_graph=False,
                )
                rows.append(
                    {
                        "slot": int(slot),
                        "celltype_id": celltype,
                        "celltype": self.artifact.celltype_names[celltype],
                        "draw_index": draw_index,
                        "normalized_celltype_weight": celltype_weight,
                        "sampled_cells": int(sampled_local_indices.numel()),
                        "huber_loss": float(
                            rescue.huber_rescue_loss.detach().float().cpu()
                        ),
                        "ccc_loss": float(
                            rescue.ccc_rescue_loss.detach().float().cpu()
                        ),
                        "huber_gradient_l2": huber_norm,
                        "ccc_gradient_l2": ccc_norm,
                    }
                )
        finally:
            trainer.optimizer.zero_grad(set_to_none=True)
            for module, training in mode_snapshot:
                module.training = training

        huber_rms = math.sqrt(
            sum(row["huber_gradient_l2"] ** 2 for row in rows) / len(rows)
        )
        ccc_rms = math.sqrt(
            sum(row["ccc_gradient_l2"] ** 2 for row in rows) / len(rows)
        )
        if not math.isfinite(ccc_rms) or ccc_rms <= 0.0:
            raise RuntimeError("CCC calibration gradient RMS is zero or non-finite")
        raw_weight = huber_rms / ccc_rms
        calibrated_weight = min(max(raw_weight, float(clip_min)), float(clip_max))
        allowlist_sha256 = hashlib.sha256(
            "\n".join(parameter_names).encode("utf-8")
        ).hexdigest()
        return {
            "schema_version": "kmlee_bam.prism_module_rescue_ccc_calibration.v1",
            "source_split": "train_only",
            "validation_or_test_used": False,
            "calibration_epoch": int(calibration_epoch),
            "calibration_blocks": int(n_blocks),
            "calibration_celltypes": [row["celltype"] for row in rows],
            "lineage_balance_mode": self.lineage_balance_plan.mode,
            "branch_fraction": float(self.config.branch_fraction),
            "formula": (
                "ccc_weight = clip(huber_gradient_rms / "
                "ccc_gradient_rms, 0.1, 1.0)"
            ),
            "huber_gradient_rms": float(huber_rms),
            "ccc_gradient_rms": float(ccc_rms),
            "raw_ccc_weight": float(raw_weight),
            "calibrated_ccc_weight": float(calibrated_weight),
            "clip_min": float(clip_min),
            "clip_max": float(clip_max),
            "target_artifact_sha256": self.artifact.artifact_sha256,
            "target_artifact_schema": self.artifact.schema_version,
            "ccc_valid_celltype_module_count": int(
                self.artifact.ccc_valid_mask.sum()
            ),
            "allowed_parameter_count": len(parameter_names),
            "allowed_parameter_names_sha256": allowlist_sha256,
            "blocks": rows,
        }

    def _run_epoch_v6a(
        self,
        trainer: Any,
        *,
        epoch: int,
        rank: int,
        world_size: int,
    ) -> dict[str, float]:
        """Run the schedule-matched v6A control or exact K=8 aggregate."""

        if epoch <= 0:
            raise ValueError("epoch must be positive")
        if world_size <= 0:
            raise ValueError("world_size must be positive")
        if world_size > 1 and (
            not dist.is_available() or not dist.is_initialized()
        ):
            raise RuntimeError("v6A DDP requires an initialized process group")
        mode = str(self.config.non_neuronal_pseudobulk_mode)
        if mode not in {
            NON_NEURONAL_PSEUDOBULK_SEPARATE_CONTROL,
            NON_NEURONAL_PSEUDOBULK_AGGREGATE,
        }:
            raise RuntimeError(f"invalid v6A mode: {mode}")

        base_system = _unwrap_system(trainer.system)
        activity_weight = base_system.module_tokenizer.activity_weight
        scale = torch.from_numpy(self.artifact.celltype_scale).to(
            device=trainer.device,
            dtype=torch.float32,
        )
        ccc_valid_mask = torch.from_numpy(
            self.artifact.ccc_valid_mask
        ).to(device=trainer.device, dtype=torch.bool)
        warmup = self.config.ramp_at_epoch(int(epoch))
        objective_config = PrismModuleRescueConfig(
            enabled=True,
            lambda_module=float(self.config.lambda_module) * warmup,
            branch_fraction=float(self.config.branch_fraction),
            huber_delta=float(self.config.huber_delta),
            ccc_weight=float(self.config.ccc_weight),
        )
        local_steps = self.steps_per_rank(world_size=world_size)
        local_per_update = self.local_blocks_per_optimizer_update(
            world_size=world_size
        )
        planned_updates = self.optimizer_updates_per_epoch(
            world_size=world_size
        )
        if local_per_update * world_size != int(
            self.config.global_blocks_per_optimizer_step
        ):
            raise RuntimeError("v6A objective normalization is not global /12")
        plan = self.lineage_balance_plan
        if len(plan.global_schedule) != local_steps * world_size:
            raise RuntimeError("v6A global schedule length mismatch")
        if len(plan.neuronal_blocks_by_update) != planned_updates:
            raise RuntimeError("v6A lineage update schedule is incomplete")

        allowed = self.resolve_parameter_allowlist(base_system, epoch=epoch)
        if allowed is None:
            raise RuntimeError("v6A requires the restricted rescue allowlist")
        non_neuronal = set(plan.non_neuronal_celltype_ids)
        metric_keys = (
            "loss", "unweighted_loss", "rescue", "huber", "ccc", "branch",
            "full", "branch_huber", "full_huber", "branch_ccc", "full_ccc",
            "branch_mean_ccc", "full_mean_ccc", "branch_corr", "full_corr",
            "branch_module_ccc", "full_module_ccc", "donor_groups",
            "pathology_axis", "pathology_axis_branch", "pathology_axis_full",
            "pathology_axis_correlation",
            "sampled_cells", "grad_norm", "allowed_grad_tensors",
            "cleared_grad_tensors",
        )
        sums = {key: 0.0 for key in metric_keys}
        block_metric_keys = tuple(
            key for key in metric_keys
            if key not in {
                "grad_norm", "allowed_grad_tensors", "cleared_grad_tensors"
            }
        )
        completed_blocks = 0
        attempted_updates = 0
        applied_updates = 0
        skipped_updates = 0
        padded_first_pass_forwards = 0
        celltype_draw_count = torch.zeros(
            len(self.sampler.eligible_celltypes),
            dtype=torch.long,
            device=trainer.device,
        )
        eligible_position = {
            int(celltype): position
            for position, celltype in enumerate(self.sampler.eligible_celltypes)
        }
        lineage_update_count = torch.zeros(
            (planned_updates, 2),
            dtype=torch.long,
            device=trainer.device,
        )
        celltype_recovery_sum = torch.zeros(
            len(self.artifact.celltype_names),
            dtype=torch.float64,
            device=trainer.device,
        )
        celltype_recovery_count = torch.zeros_like(celltype_recovery_sum)

        def sync_context():
            ddp = getattr(trainer.system, "_ddp_model", None)
            return (
                ddp.no_sync()
                if ddp is not None and world_size > 1
                else nullcontext()
            )

        def move_batch(raw_batch: dict[str, Any]) -> dict[str, Any]:
            batch = trainer._move_batch_to_device(raw_batch)
            batch["donor_balance_weight"] = torch.ones(
                batch["donor_id"].shape,
                dtype=torch.float32,
                device=trainer.device,
            )
            return batch

        def grouped_forward(
            raw_batch: dict[str, Any],
        ) -> tuple[GroupedModuleActivity, GroupedModuleActivity]:
            batch = move_batch(raw_batch)
            with trainer._autocast_context():
                model_out = trainer.system(
                    batch,
                    sample_latent=False,
                    return_all_hidden_states=False,
                    return_attn_diagnostics=False,
                )
                if model_out.decoder_out is None:
                    raise RuntimeError("v6A rescue requires decoder output")
                branch_probability = self._branch_probabilities(
                    trainer, batch, model_out
                )
                group_id = (
                    batch["donor_id"].long()
                    * len(self.artifact.celltype_names)
                    + batch["celltype_id"].long()
                )
                branch_group = grouped_expected_module_activity(
                    branch_probability,
                    activity_weight,
                    group_id,
                    batch["celltype_id"],
                    validate_probabilities=True,
                    probability_atol=5e-3,
                )
                full_group = grouped_expected_module_activity(
                    model_out.decoder_out.probs,
                    activity_weight,
                    group_id,
                    batch["celltype_id"],
                    validate_probabilities=True,
                    probability_atol=5e-3,
                )
            return branch_group, full_group

        def group_contract(celltype: int):
            donors = [
                donor for donor in self.donor_ids
                if (int(celltype), int(donor)) in self.sampler.groups
            ]
            group_id = torch.as_tensor(
                [
                    donor * len(self.artifact.celltype_names) + int(celltype)
                    for donor in donors
                ],
                dtype=torch.long,
                device=trainer.device,
            )
            group_celltype = torch.full_like(group_id, int(celltype))
            target = torch.as_tensor(
                self.artifact.group_module[donors, int(celltype)],
                dtype=torch.float32,
                device=trainer.device,
            )
            target_count = torch.as_tensor(
                self.artifact.group_count[donors, int(celltype)],
                dtype=torch.long,
                device=trainer.device,
            )
            return group_id, group_celltype, target, target_count

        def loss_builder(
            *,
            group_id: torch.Tensor,
            group_celltype: torch.Tensor,
            prediction_count: torch.Tensor,
            target: torch.Tensor,
            target_count: torch.Tensor,
        ):
            celltype = int(group_celltype[0].item())

            def build(branch_mean, full_mean):
                output = dual_view_grouped_module_rescue_loss(
                    config=objective_config,
                    branch_group_activity=branch_mean,
                    full_group_activity=full_mean,
                    prediction_group_id=group_id,
                    prediction_group_celltype_id=group_celltype,
                    prediction_group_cell_count=prediction_count,
                    target_group_activity=target,
                    target_group_id=group_id,
                    target_group_cell_count=target_count,
                    train_scale=scale,
                    ccc_valid_mask=ccc_valid_mask,
                    min_cells_per_group=1,
                    min_donor_groups_per_celltype=int(
                        self.config.minimum_donors_per_celltype
                    ),
                )
                return self._add_pathology_axis_loss(
                    output,
                    group_id=group_id,
                    celltype=celltype,
                    epoch=epoch,
                )

            return build

        def record_output(
            output: PrismModuleRescueOutput,
            *,
            celltype_weight: float,
            multiplier: float,
        ) -> None:
            sums["loss"] += multiplier * celltype_weight * float(
                output.loss.detach().float().cpu()
            )
            scalar_fields = {
                "unweighted_loss": output.loss,
                "rescue": output.rescue_loss,
                "huber": output.huber_rescue_loss,
                "ccc": output.ccc_rescue_loss,
                "branch": output.branch_loss,
                "full": output.full_loss,
                "branch_huber": output.branch_huber_loss,
                "full_huber": output.full_huber_loss,
                "branch_ccc": output.branch_ccc_loss,
                "full_ccc": output.full_ccc_loss,
                "branch_mean_ccc": output.branch_mean_ccc,
                "full_mean_ccc": output.full_mean_ccc,
            }
            optional_fields = {
                "pathology_axis": output.pathology_axis_loss,
                "pathology_axis_branch": output.pathology_axis_branch_loss,
                "pathology_axis_full": output.pathology_axis_full_loss,
                "pathology_axis_correlation": (
                    output.pathology_axis_mean_correlation
                ),
            }
            scalar_fields.update(
                {key: value for key, value in optional_fields.items() if value is not None}
            )
            for key, value in scalar_fields.items():
                sums[key] += multiplier * float(value.detach().float().cpu())
            if output.branch_view is not None:
                valid = output.branch_view.metrics.valid
                if bool(valid.any()):
                    sums["branch_corr"] += multiplier * float(
                        output.branch_view.metrics.correlation[valid]
                        .mean().detach().cpu()
                    )
                    sums["branch_module_ccc"] += multiplier * float(
                        output.branch_view.metrics.ccc[valid]
                        .mean().detach().cpu()
                    )
                sums["donor_groups"] += multiplier * float(
                    output.branch_view.valid_group.sum().detach().cpu()
                )
            if output.full_view is not None:
                valid = output.full_view.metrics.valid
                if bool(valid.any()):
                    sums["full_corr"] += multiplier * float(
                        output.full_view.metrics.correlation[valid]
                        .mean().detach().cpu()
                    )
                    sums["full_module_ccc"] += multiplier * float(
                        output.full_view.metrics.ccc[valid]
                        .mean().detach().cpu()
                    )
            celltype = int(output.full_view.group_celltype_id[0].item())
            recovery = output.pathology_axis_mean_correlation
            if recovery is None:
                valid = output.full_view.metrics.valid
                if bool(valid.any()):
                    recovery = output.full_view.metrics.correlation[valid].mean()
            if recovery is not None and bool(torch.isfinite(recovery)):
                celltype_recovery_sum[celltype] += multiplier * recovery.detach().double()
                celltype_recovery_count[celltype] += multiplier

        def sample_position(global_position: int):
            local_slot = global_position // world_size
            raw = self._sample_rescue_block(
                epoch=epoch,
                slot=local_slot,
                rank=rank,
                world_size=world_size,
            )
            celltype = int(raw.pop("_module_rescue_celltype"))
            draw = int(raw.pop("_module_rescue_draw_index"))
            indices = raw.pop("_module_rescue_local_indices")
            return {
                "raw": raw,
                "celltype": celltype,
                "draw": draw,
                "indices": indices,
                "local_slot": local_slot,
            }

        ddp_model = getattr(trainer.system, "_ddp_model", None)
        ddp_training_snapshot = (
            bool(ddp_model.training) if ddp_model is not None else None
        )
        mode_snapshot = [
            (module, bool(module.training))
            for module in trainer.system.modules()
        ]
        trainer.system.eval()
        try:
            with _temporarily_freeze_outside_allowlist(base_system, allowed):
                for update_index in range(planned_updates):
                    trainer.optimizer.zero_grad(set_to_none=True)
                    start = update_index * int(
                        self.config.global_blocks_per_optimizer_step
                    )
                    stop = start + int(
                        self.config.global_blocks_per_optimizer_step
                    )
                    local_items = [
                        sample_position(position)
                        for position in range(start, stop)
                        if position % world_size == rank
                    ]
                    if len(local_items) != local_per_update:
                        raise RuntimeError("v6A local schedule is imbalanced")
                    for item in local_items:
                        celltype_draw_count[
                            eligible_position[item["celltype"]]
                        ] += 1
                        lineage = 1 if item["celltype"] in non_neuronal else 0
                        lineage_update_count[update_index, lineage] += 1
                        sums["sampled_cells"] += float(
                            item["indices"].numel()
                        )
                    non_neuronal_items = [
                        item for item in local_items
                        if item["celltype"] in non_neuronal
                    ]
                    neuronal_items = [
                        item for item in local_items
                        if item["celltype"] not in non_neuronal
                    ]
                    scheduled_non_neuronal = int(
                        plan.non_neuronal_blocks_by_update[update_index]
                    )

                    if scheduled_non_neuronal:
                        aggregate_celltype = next(
                            int(celltype)
                            for celltype, _ in plan.global_schedule[start:stop]
                            if int(celltype) in non_neuronal
                        )
                        if any(
                            item["celltype"] != aggregate_celltype
                            for item in non_neuronal_items
                        ):
                            raise RuntimeError(
                                "v6A update mixed non-neuronal cell types"
                            )
                        group_id, group_celltype, target, target_count = (
                            group_contract(aggregate_celltype)
                        )
                        first_pass = []
                        for item in non_neuronal_items:
                            with sync_context(), torch.no_grad():
                                branch_group, full_group = grouped_forward(
                                    item["raw"]
                                )
                            if not bool(torch.equal(branch_group.group_id, group_id)):
                                raise RuntimeError("v6A sampled donor groups changed")
                            first_pass.append((item, branch_group, full_group))
                        local_g_counts = []
                        for candidate_rank in range(world_size):
                            local_g_counts.append(sum(
                                1
                                for position in range(start, stop)
                                if position % world_size == candidate_rank
                                and plan.global_schedule[position][0]
                                in non_neuronal
                            ))
                        max_local_g = max(local_g_counts)
                        padding = max_local_g - len(first_pass)
                        if padding < 0:
                            raise RuntimeError("v6A first-pass cadence underflow")
                        if padding:
                            # DDP may broadcast buffers in pre-forward even
                            # inside no_sync.  Match the number of first-pass
                            # forwards on every rank before aggregate all-reduce.
                            dummy_raw = local_items[0]["raw"]
                            for _ in range(padding):
                                with sync_context(), torch.no_grad():
                                    grouped_forward(dummy_raw)
                            padded_first_pass_forwards += padding

                        celltype_weight = self._combined_celltype_weight(
                            aggregate_celltype
                        )
                        if mode == NON_NEURONAL_PSEUDOBULK_AGGREGATE:
                            modules = int(activity_weight.shape[0])
                            branch_sum = torch.zeros(
                                (int(group_id.numel()), modules),
                                dtype=torch.float32,
                                device=trainer.device,
                            )
                            full_sum = torch.zeros_like(branch_sum)
                            prediction_count = torch.zeros(
                                int(group_id.numel()),
                                dtype=torch.long,
                                device=trainer.device,
                            )
                            for _, branch_group, full_group in first_pass:
                                count = branch_group.cell_count
                                branch_sum.add_(
                                    branch_group.mean.detach()
                                    * count.to(branch_group.mean.dtype)[:, None]
                                )
                                full_sum.add_(
                                    full_group.mean.detach()
                                    * count.to(full_group.mean.dtype)[:, None]
                                )
                                prediction_count.add_(count)
                            if world_size > 1:
                                dist.all_reduce(branch_sum, op=dist.ReduceOp.SUM)
                                dist.all_reduce(full_sum, op=dist.ReduceOp.SUM)
                                dist.all_reduce(
                                    prediction_count, op=dist.ReduceOp.SUM
                                )
                            denominator = prediction_count.clamp_min(1).to(
                                branch_sum.dtype
                            )[:, None]
                            leaf_vjp = _exact_group_mean_leaf_vjp(
                                branch_group_mean=branch_sum / denominator,
                                full_group_mean=full_sum / denominator,
                                loss_builder=loss_builder(
                                    group_id=group_id,
                                    group_celltype=group_celltype,
                                    prediction_count=prediction_count,
                                    target=target,
                                    target_count=target_count,
                                ),
                                loss_scale=(
                                    celltype_weight
                                    * float(scheduled_non_neuronal)
                                    / float(local_per_update)
                                ),
                            )
                            for item, _, _ in first_pass:
                                with sync_context():
                                    branch_group, full_group = grouped_forward(
                                        item["raw"]
                                    )
                                    proxy = _grouped_mean_replay_proxy(
                                        branch_chunk=branch_group,
                                        full_chunk=full_group,
                                        aggregate_group_id=group_id,
                                        aggregate_cell_count=prediction_count,
                                        branch_gradient=leaf_vjp.branch_gradient,
                                        full_gradient=leaf_vjp.full_gradient,
                                    )
                                    if not trainer._backward(
                                        proxy,
                                        step_idx=int(item["local_slot"]) + 1,
                                    ):
                                        raise RuntimeError(
                                            "v6A aggregate replay was non-finite"
                                        )
                            # Every rank computes the same group loss.  Record
                            # it once before block-metric all-reduce.
                            if rank == 0:
                                record_output(
                                    leaf_vjp.output,
                                    celltype_weight=celltype_weight,
                                    multiplier=float(scheduled_non_neuronal),
                                )
                        else:
                            for item, branch_group, full_group in first_pass:
                                prediction_count = branch_group.cell_count
                                leaf_vjp = _exact_group_mean_leaf_vjp(
                                    branch_group_mean=branch_group.mean,
                                    full_group_mean=full_group.mean,
                                    loss_builder=loss_builder(
                                        group_id=group_id,
                                        group_celltype=group_celltype,
                                        prediction_count=prediction_count,
                                        target=target,
                                        target_count=target_count,
                                    ),
                                    loss_scale=(
                                        celltype_weight
                                        / float(local_per_update)
                                    ),
                                )
                                with sync_context():
                                    replay_branch, replay_full = grouped_forward(
                                        item["raw"]
                                    )
                                    proxy = _grouped_mean_replay_proxy(
                                        branch_chunk=replay_branch,
                                        full_chunk=replay_full,
                                        aggregate_group_id=group_id,
                                        aggregate_cell_count=prediction_count,
                                        branch_gradient=leaf_vjp.branch_gradient,
                                        full_gradient=leaf_vjp.full_gradient,
                                    )
                                    if not trainer._backward(
                                        proxy,
                                        step_idx=int(item["local_slot"]) + 1,
                                    ):
                                        raise RuntimeError(
                                            "v6A control replay was non-finite"
                                        )
                                record_output(
                                    leaf_vjp.output,
                                    celltype_weight=celltype_weight,
                                    multiplier=1.0,
                                )

                    for item in neuronal_items:
                        celltype = int(item["celltype"])
                        with sync_context():
                            branch_group, full_group = grouped_forward(item["raw"])
                            group_id, group_celltype, target, target_count = (
                                group_contract(celltype)
                            )
                            if not bool(torch.equal(branch_group.group_id, group_id)):
                                raise RuntimeError("neuronal donor groups changed")
                            output = dual_view_grouped_module_rescue_loss(
                                config=objective_config,
                                branch_group_activity=branch_group.mean,
                                full_group_activity=full_group.mean,
                                prediction_group_id=group_id,
                                prediction_group_celltype_id=group_celltype,
                                prediction_group_cell_count=(
                                    branch_group.cell_count
                                ),
                                target_group_activity=target,
                                target_group_id=group_id,
                                target_group_cell_count=target_count,
                                train_scale=scale,
                                ccc_valid_mask=ccc_valid_mask,
                                min_cells_per_group=1,
                                min_donor_groups_per_celltype=int(
                                    self.config.minimum_donors_per_celltype
                                ),
                            )
                            output = self._add_pathology_axis_loss(
                                output,
                                group_id=group_id,
                                celltype=celltype,
                                epoch=epoch,
                            )
                            celltype_weight = self._combined_celltype_weight(celltype)
                            if not trainer._backward(
                                output.loss
                                * celltype_weight
                                / float(local_per_update),
                                step_idx=int(item["local_slot"]) + 1,
                            ):
                                raise RuntimeError(
                                    "v6A neuronal rescue was non-finite"
                                )
                        record_output(
                            output,
                            celltype_weight=celltype_weight,
                            multiplier=1.0,
                        )

                    completed_blocks += len(local_items)
                    allowed_used, cleared = _apply_module_rescue_gradient_allowlist(
                        base_system, allowed
                    )
                    _manual_average_gradients(
                        base_system,
                        world_size=world_size,
                    )
                    attempted_updates += 1
                    original_weight_decay = [
                        group.get("weight_decay", 0.0)
                        for group in trainer.optimizer.param_groups
                    ]
                    try:
                        for group in trainer.optimizer.param_groups:
                            group["weight_decay"] = 0.0
                        grad_norm = trainer._finish_optimizer_step(
                            advance_scheduler=False
                        )
                    finally:
                        for group, value in zip(
                            trainer.optimizer.param_groups,
                            original_weight_decay,
                        ):
                            group["weight_decay"] = value
                    if bool(getattr(trainer, "_last_step_was_skipped", False)):
                        skipped_updates += 1
                        raise RuntimeError("v6A optimizer update was skipped")
                    applied_updates += 1
                    if rank == 0:
                        sums["grad_norm"] += float(grad_norm or 0.0)
                        sums["allowed_grad_tensors"] += float(allowed_used)
                        sums["cleared_grad_tensors"] += float(cleared)
        finally:
            trainer.optimizer.zero_grad(set_to_none=True)
            for module, training in mode_snapshot:
                module.training = training
            if ddp_model is not None and ddp_training_snapshot is not None:
                ddp_model.training = ddp_training_snapshot

        if completed_blocks != local_steps:
            raise RuntimeError(
                f"v6A completed {completed_blocks} != {local_steps} local blocks"
            )
        if world_size > 1:
            dist.all_reduce(celltype_draw_count, op=dist.ReduceOp.SUM)
            dist.all_reduce(lineage_update_count, op=dist.ReduceOp.SUM)
            dist.all_reduce(celltype_recovery_sum, op=dist.ReduceOp.SUM)
            dist.all_reduce(celltype_recovery_count, op=dist.ReduceOp.SUM)
            block_values = torch.as_tensor(
                [sums[key] for key in block_metric_keys],
                dtype=torch.float64,
                device=trainer.device,
            )
            dist.all_reduce(block_values, op=dist.ReduceOp.SUM)
            for key, value in zip(block_metric_keys, block_values.tolist()):
                sums[key] = float(value)
            update_values = torch.as_tensor(
                [
                    sums["grad_norm"],
                    sums["allowed_grad_tensors"],
                    sums["cleared_grad_tensors"],
                ],
                dtype=torch.float64,
                device=trainer.device,
            )
            dist.all_reduce(update_values, op=dist.ReduceOp.SUM)
            for key, value in zip(
                ("grad_norm", "allowed_grad_tensors", "cleared_grad_tensors"),
                update_values.tolist(),
            ):
                sums[key] = float(value)
            padding_tensor = torch.tensor(
                padded_first_pass_forwards,
                dtype=torch.long,
                device=trainer.device,
            )
            dist.all_reduce(padding_tensor, op=dist.ReduceOp.SUM)
            padded_first_pass_forwards = int(padding_tensor.item())
        expected_draws = int(self.config.draws_per_celltype_per_epoch)
        if not bool((celltype_draw_count == expected_draws).all()):
            raise RuntimeError(
                "v6A celltype draw audit failed: "
                f"{celltype_draw_count.detach().cpu().tolist()}"
            )
        expected_lineage = torch.as_tensor(
            list(zip(
                plan.neuronal_blocks_by_update,
                plan.non_neuronal_blocks_by_update,
            )),
            dtype=torch.long,
            device=trainer.device,
        )
        if not bool(torch.equal(lineage_update_count, expected_lineage)):
            raise RuntimeError(
                "v6A lineage schedule audit failed: expected "
                f"{expected_lineage.detach().cpu().tolist()} observed "
                f"{lineage_update_count.detach().cpu().tolist()}"
            )
        if applied_updates != planned_updates or skipped_updates:
            raise RuntimeError("v6A optimizer update count contract failed")

        recovery = torch.full_like(celltype_recovery_sum, float("nan"))
        recovery_valid = celltype_recovery_count > 0
        recovery[recovery_valid] = (
            celltype_recovery_sum[recovery_valid]
            / celltype_recovery_count[recovery_valid]
        )
        self._update_tail_weights(
            recovery.detach().cpu().numpy(),
            epoch=epoch,
        )

        global_blocks = completed_blocks * world_size
        block_denominator = max(global_blocks, 1)
        update_denominator = max(applied_updates, 1)
        neuronal_blocks = float(lineage_update_count[:, 0].sum().cpu())
        non_neuronal_blocks = float(lineage_update_count[:, 1].sum().cpu())
        neuronal_weight = float(
            plan.celltype_normalized_weights[plan.neuronal_celltype_ids[0]]
        )
        non_neuronal_weight = float(
            plan.celltype_normalized_weights[
                plan.non_neuronal_celltype_ids[0]
            ]
        )
        neuronal_mass = neuronal_blocks * neuronal_weight
        non_neuronal_mass = non_neuronal_blocks * non_neuronal_weight
        total_mass = neuronal_mass + non_neuronal_mass
        metrics = {
            "loss/prism_module_rescue": sums["loss"] / block_denominator,
            "loss/prism_module_rescue_objective_unweighted": sums["unweighted_loss"] / block_denominator,
            "metric/prism_module_rescue_unweighted": sums["rescue"] / block_denominator,
            "loss/prism_module_rescue_huber_unweighted": sums["huber"] / block_denominator,
            "loss/prism_module_rescue_ccc_unweighted": sums["ccc"] / block_denominator,
            "loss/prism_module_rescue_pathology_axis": sums["pathology_axis"] / block_denominator,
            "loss/prism_module_rescue_pathology_axis_branch": sums["pathology_axis_branch"] / block_denominator,
            "loss/prism_module_rescue_pathology_axis_full": sums["pathology_axis_full"] / block_denominator,
            "metric/prism_module_rescue_pathology_axis_correlation": sums["pathology_axis_correlation"] / block_denominator,
            "loss/prism_module_rescue_branch": sums["branch"] / block_denominator,
            "loss/prism_module_rescue_full": sums["full"] / block_denominator,
            "loss/prism_module_rescue_branch_huber": sums["branch_huber"] / block_denominator,
            "loss/prism_module_rescue_full_huber": sums["full_huber"] / block_denominator,
            "loss/prism_module_rescue_branch_ccc": sums["branch_ccc"] / block_denominator,
            "loss/prism_module_rescue_full_ccc": sums["full_ccc"] / block_denominator,
            "metric/prism_module_rescue_branch_ccc": sums["branch_module_ccc"] / block_denominator,
            "metric/prism_module_rescue_full_ccc": sums["full_module_ccc"] / block_denominator,
            "metric/prism_module_rescue_branch_flattened_ccc": sums["branch_mean_ccc"] / block_denominator,
            "metric/prism_module_rescue_full_flattened_ccc": sums["full_mean_ccc"] / block_denominator,
            "metric/prism_module_rescue_branch_pearson": sums["branch_corr"] / block_denominator,
            "metric/prism_module_rescue_full_pearson": sums["full_corr"] / block_denominator,
            "metric/prism_module_rescue_groups": sums["donor_groups"] / block_denominator,
            "metric/prism_module_rescue_grad_norm": sums["grad_norm"] / update_denominator,
            "metric/prism_module_rescue_allowed_grad_tensors": sums["allowed_grad_tensors"] / update_denominator,
            "metric/prism_module_rescue_cleared_grad_tensors": sums["cleared_grad_tensors"] / update_denominator,
            "metric/prism_module_rescue_sampled_cells_per_block": sums["sampled_cells"] / block_denominator,
            "metric/prism_module_rescue_celltype_draw_min": float(celltype_draw_count.min().cpu()),
            "metric/prism_module_rescue_celltype_draw_max": float(celltype_draw_count.max().cpu()),
            "metric/prism_module_rescue_blocks": float(completed_blocks),
            "metric/prism_module_rescue_global_block_equivalents": float(global_blocks),
            "metric/prism_module_rescue_updates_attempted": float(attempted_updates),
            "metric/prism_module_rescue_updates_applied": float(applied_updates),
            "metric/prism_module_rescue_updates_skipped": float(skipped_updates),
            "metric/prism_module_rescue_global_blocks_per_update": float(local_per_update * world_size),
            "metric/prism_module_rescue_lineage_balance_enabled": 1.0,
            "metric/prism_module_rescue_neuronal_blocks": neuronal_blocks,
            "metric/prism_module_rescue_non_neuronal_blocks": non_neuronal_blocks,
            "metric/prism_module_rescue_neuronal_blocks_per_update_min": float(lineage_update_count[:, 0].min().cpu()),
            "metric/prism_module_rescue_neuronal_blocks_per_update_max": float(lineage_update_count[:, 0].max().cpu()),
            "metric/prism_module_rescue_non_neuronal_blocks_per_update_min": float(lineage_update_count[:, 1].min().cpu()),
            "metric/prism_module_rescue_non_neuronal_blocks_per_update_max": float(lineage_update_count[:, 1].max().cpu()),
            "metric/prism_module_rescue_neuronal_objective_mass_fraction": neuronal_mass / total_mass,
            "metric/prism_module_rescue_non_neuronal_objective_mass_fraction": non_neuronal_mass / total_mass,
            "metric/prism_module_rescue_v6a_aggregate_enabled": float(mode == NON_NEURONAL_PSEUDOBULK_AGGREGATE),
            "metric/prism_module_rescue_v6a_aggregate_k": float(int(self.config.cells_per_donor_celltype) * int(self.config.non_neuronal_aggregate_draws)),
            "metric/prism_module_rescue_v6a_padded_first_pass_forwards": float(padded_first_pass_forwards),
            "weight/prism_module_rescue_neuronal_raw": float(plan.celltype_raw_weights[plan.neuronal_celltype_ids[0]]),
            "weight/prism_module_rescue_non_neuronal_raw": float(plan.celltype_raw_weights[plan.non_neuronal_celltype_ids[0]]),
            "weight/prism_module_rescue_neuronal_normalized": neuronal_weight,
            "weight/prism_module_rescue_non_neuronal_normalized": non_neuronal_weight,
            "weight/prism_module_rescue_ramp": float(warmup),
            "weight/prism_module_rescue_ccc": float(self.config.ccc_weight),
            "metric/prism_module_rescue_steps": float(applied_updates),
        }
        for celltype, name in enumerate(self.artifact.celltype_names):
            metrics[f"weight/prism_module_rescue_celltype/{name}"] = float(
                self._combined_celltype_weight(celltype)
            )
            metrics[f"metric/prism_module_rescue_recovery/{name}"] = float(
                recovery[celltype].cpu()
            )
        return metrics

    def run_epoch(
        self,
        trainer: Any,
        *,
        epoch: int,
        rank: int,
        world_size: int,
    ) -> dict[str, float]:
        if int(epoch) <= 0:
            raise ValueError("epoch must be positive")
        if int(epoch) < int(self.config.start_epoch):
            return {
                "loss/prism_module_rescue": 0.0,
                "loss/prism_module_rescue_branch": 0.0,
                "loss/prism_module_rescue_full": 0.0,
                "metric/prism_module_rescue_branch_flattened_ccc": 0.0,
                "metric/prism_module_rescue_full_flattened_ccc": 0.0,
                "metric/prism_module_rescue_branch_pearson": 0.0,
                "metric/prism_module_rescue_full_pearson": 0.0,
                "metric/prism_module_rescue_groups": 0.0,
                "metric/prism_module_rescue_grad_norm": 0.0,
                "metric/prism_module_rescue_celltype_draw_min": 0.0,
                "metric/prism_module_rescue_celltype_draw_max": 0.0,
                "metric/prism_module_rescue_steps": 0.0,
                "weight/prism_module_rescue_ramp": 0.0,
            }
        if (
            str(self.config.non_neuronal_pseudobulk_mode)
            != NON_NEURONAL_PSEUDOBULK_LEGACY
        ):
            return self._run_epoch_v6a(
                trainer,
                epoch=int(epoch),
                rank=int(rank),
                world_size=int(world_size),
            )
        world_size = int(world_size)
        base_system = _unwrap_system(trainer.system)
        activity_weight = base_system.module_tokenizer.activity_weight
        scale = torch.from_numpy(self.artifact.celltype_scale).to(
            device=trainer.device,
            dtype=torch.float32,
        )
        ccc_valid_mask = torch.from_numpy(
            self.artifact.ccc_valid_mask
        ).to(device=trainer.device, dtype=torch.bool)
        warmup = self.config.ramp_at_epoch(int(epoch))
        objective_config = PrismModuleRescueConfig(
            enabled=True,
            lambda_module=float(self.config.lambda_module) * warmup,
            branch_fraction=float(self.config.branch_fraction),
            huber_delta=float(self.config.huber_delta),
            ccc_weight=float(self.config.ccc_weight),
        )
        metric_keys = (
            "loss", "unweighted_loss", "rescue", "huber", "ccc", "branch", "full",
            "branch_huber", "full_huber", "branch_ccc", "full_ccc",
            "branch_mean_ccc", "full_mean_ccc", "branch_corr",
            "full_corr", "branch_module_ccc", "full_module_ccc",
            "grad_norm", "donor_groups",
            "allowed_grad_tensors", "cleared_grad_tensors", "sampled_cells",
        )
        sums = {key: 0.0 for key in metric_keys}
        completed_blocks = 0
        attempted_updates = 0
        applied_updates = 0
        skipped_updates = 0
        local_steps = self.steps_per_rank(world_size=world_size)
        local_per_update = self.local_blocks_per_optimizer_update(
            world_size=world_size
        )
        planned_updates = self.optimizer_updates_per_epoch(
            world_size=world_size
        )
        self.resolve_parameter_allowlist(base_system, epoch=epoch)

        eligible_position = {
            int(celltype): position
            for position, celltype in enumerate(self.sampler.eligible_celltypes)
        }
        celltype_draw_count = torch.zeros(
            len(self.sampler.eligible_celltypes),
            dtype=torch.long,
            device=trainer.device,
        )
        lineage_update_count = torch.zeros(
            (planned_updates, 2),
            dtype=torch.long,
            device=trainer.device,
        )
        neuronal_celltypes = set(
            self.lineage_balance_plan.neuronal_celltype_ids
        )
        non_neuronal_celltypes = set(
            self.lineage_balance_plan.non_neuronal_celltype_ids
        )
        # Preserve heterogeneous train/eval flags exactly, even if an
        # exception aborts an auxiliary block halfway through the epoch.
        ddp_model = getattr(trainer.system, "_ddp_model", None)
        ddp_training_snapshot = (
            bool(ddp_model.training) if ddp_model is not None else None
        )
        mode_snapshot = [
            (module, bool(module.training))
            for module in trainer.system.modules()
        ]
        if bool(self.config.deterministic_forward):
            trainer.system.eval()
        try:
            for update_index in range(planned_updates):
                trainer.optimizer.zero_grad(set_to_none=True)
                for local_offset in range(local_per_update):
                    slot = update_index * local_per_update + local_offset
                    raw_batch = self._sample_rescue_block(
                        epoch=int(epoch),
                        slot=slot,
                        rank=int(rank),
                        world_size=world_size,
                    )
                    celltype = int(raw_batch.pop("_module_rescue_celltype"))
                    raw_batch.pop("_module_rescue_draw_index")
                    sampled_local_indices = raw_batch.pop(
                        "_module_rescue_local_indices"
                    )
                    sums["sampled_cells"] += float(
                        sampled_local_indices.numel()
                    )
                    celltype_draw_count[eligible_position[celltype]] += 1
                    if self.lineage_balance_plan.enabled:
                        if celltype in neuronal_celltypes:
                            lineage_update_count[update_index, 0] += 1
                        elif celltype in non_neuronal_celltypes:
                            lineage_update_count[update_index, 1] += 1
                        else:
                            raise RuntimeError(
                                "lineage-balanced rescue sampled an "
                                f"unclassified cell type: {celltype}"
                            )
                    batch = trainer._move_batch_to_device(raw_batch)
                    if bool(self.config.equal_donor_rescue_gradient):
                        # The grouped sampler already gives every donor equal
                        # cell count.  Override the ordinary inverse-frequency
                        # wrapper so the full view is not donor-weighted twice.
                        batch["donor_balance_weight"] = torch.ones(
                            batch["donor_id"].shape,
                            dtype=torch.float32,
                            device=trainer.device,
                        )
                        if not bool(
                            (batch["donor_balance_weight"] == 1).all()
                        ):
                            raise RuntimeError(
                                "rescue donor-balance override failed"
                            )
                    sync_context = (
                        ddp_model.no_sync()
                        if ddp_model is not None and world_size > 1
                        else nullcontext()
                    )
                    # The DDP reducer stays disarmed for every accumulated
                    # block; one explicit mean is performed at update boundary.
                    with sync_context:
                        with trainer._autocast_context():
                            model_out = trainer.system(
                                batch,
                                sample_latent=not bool(
                                    self.config.deterministic_forward
                                ),
                                return_all_hidden_states=False,
                                return_attn_diagnostics=False,
                            )
                            if model_out.decoder_out is None:
                                raise RuntimeError(
                                    "module-rescue update requires decoder output"
                                )
                            branch_probabilities = self._branch_probabilities(
                                trainer, batch, model_out
                            )
                            donor = batch["donor_id"].long()
                            group_id = (
                                donor * len(self.artifact.celltype_names)
                                + batch["celltype_id"].long()
                            )
                            unique_group = torch.unique(group_id, sorted=True)
                            target_module = []
                            target_count = []
                            for group in unique_group.detach().cpu().tolist():
                                donor_id = int(group) // len(
                                    self.artifact.celltype_names
                                )
                                ct_id = int(group) % len(
                                    self.artifact.celltype_names
                                )
                                target_module.append(
                                    self.artifact.group_module[donor_id, ct_id]
                                )
                                target_count.append(
                                    self.artifact.group_count[donor_id, ct_id]
                                )
                            target_module_tensor = torch.as_tensor(
                                np.asarray(target_module, dtype=np.float32),
                                device=trainer.device,
                            )
                            target_count_tensor = torch.as_tensor(
                                np.asarray(target_count, dtype=np.int64),
                                device=trainer.device,
                            )
                            rescue = dual_view_module_rescue_loss(
                                config=objective_config,
                                branch_probabilities=branch_probabilities,
                                full_probabilities=model_out.decoder_out.probs,
                                target_group_activity=target_module_tensor,
                                target_group_id=unique_group,
                                target_group_cell_count=target_count_tensor,
                                activity_weight=activity_weight,
                                group_id=group_id,
                                celltype_id=batch["celltype_id"],
                                train_scale=scale,
                                ccc_valid_mask=ccc_valid_mask,
                                min_cells_per_group=1,
                                min_donor_groups_per_celltype=int(
                                    self.config.minimum_donors_per_celltype
                                ),
                            )

                        celltype_weight = float(
                            self.lineage_balance_plan
                            .celltype_normalized_weights[celltype]
                        )
                        weighted_loss = rescue.loss
                        if self.lineage_balance_plan.enabled:
                            weighted_loss = weighted_loss * celltype_weight
                        backward_ok = trainer._backward(
                            weighted_loss / float(local_per_update),
                            step_idx=slot + 1,
                        )
                    if not backward_ok:
                        trainer.optimizer.zero_grad(set_to_none=True)
                        if str(self.config.contract_version) == "v2":
                            raise RuntimeError(
                                "module-rescue v2 encountered a non-finite "
                                f"loss at local block {slot}"
                            )
                        break

                    completed_blocks += 1
                    sums["loss"] += float(
                        weighted_loss.detach().float().cpu()
                    )
                    sums["unweighted_loss"] += float(
                        rescue.loss.detach().float().cpu()
                    )
                    sums["rescue"] += float(
                        rescue.rescue_loss.detach().float().cpu()
                    )
                    sums["huber"] += float(
                        rescue.huber_rescue_loss.detach().float().cpu()
                    )
                    sums["ccc"] += float(
                        rescue.ccc_rescue_loss.detach().float().cpu()
                    )
                    sums["branch"] += float(
                        rescue.branch_loss.detach().float().cpu()
                    )
                    sums["full"] += float(
                        rescue.full_loss.detach().float().cpu()
                    )
                    sums["branch_huber"] += float(
                        rescue.branch_huber_loss.detach().float().cpu()
                    )
                    sums["full_huber"] += float(
                        rescue.full_huber_loss.detach().float().cpu()
                    )
                    sums["branch_ccc"] += float(
                        rescue.branch_ccc_loss.detach().float().cpu()
                    )
                    sums["full_ccc"] += float(
                        rescue.full_ccc_loss.detach().float().cpu()
                    )
                    sums["branch_mean_ccc"] += float(
                        rescue.branch_mean_ccc.detach().float().cpu()
                    )
                    sums["full_mean_ccc"] += float(
                        rescue.full_mean_ccc.detach().float().cpu()
                    )
                    if rescue.branch_view is not None:
                        valid = rescue.branch_view.metrics.valid
                        if bool(valid.any()):
                            sums["branch_corr"] += float(
                                rescue.branch_view.metrics.correlation[valid]
                                .mean().detach().cpu()
                            )
                            sums["branch_module_ccc"] += float(
                                rescue.branch_view.metrics.ccc[valid]
                                .mean().detach().cpu()
                            )
                        sums["donor_groups"] += float(
                            rescue.branch_view.valid_group.sum().detach().cpu()
                        )
                    if rescue.full_view is not None:
                        valid = rescue.full_view.metrics.valid
                        if bool(valid.any()):
                            sums["full_corr"] += float(
                                rescue.full_view.metrics.correlation[valid]
                                .mean().detach().cpu()
                            )
                            sums["full_module_ccc"] += float(
                                rescue.full_view.metrics.ccc[valid]
                                .mean().detach().cpu()
                            )

                if self._parameter_allowlist is not None:
                    allowed_used, cleared = (
                        _apply_module_rescue_gradient_allowlist(
                            base_system,
                            self._parameter_allowlist,
                        )
                    )
                    sums["allowed_grad_tensors"] += float(allowed_used)
                    sums["cleared_grad_tensors"] += float(cleared)
                _manual_average_gradients(
                    base_system,
                    world_size=world_size,
                )
                attempted_updates += 1
                # Auxiliary loss must not silently apply extra AdamW decay or
                # advance a step-wise scheduler.  Adam moments for the allowed
                # parameters still update, as they should for a real loss step.
                original_weight_decay = [
                    group.get("weight_decay", 0.0)
                    for group in trainer.optimizer.param_groups
                ]
                try:
                    for group in trainer.optimizer.param_groups:
                        group["weight_decay"] = 0.0
                    grad_norm = trainer._finish_optimizer_step(
                        advance_scheduler=False
                    )
                finally:
                    for group, value in zip(
                        trainer.optimizer.param_groups,
                        original_weight_decay,
                    ):
                        group["weight_decay"] = value
                if bool(getattr(trainer, "_last_step_was_skipped", False)):
                    skipped_updates += 1
                    if str(self.config.contract_version) == "v2":
                        raise RuntimeError(
                            "module-rescue v2 optimizer update was skipped"
                        )
                else:
                    applied_updates += 1
                    sums["grad_norm"] += float(grad_norm or 0.0)
        finally:
            trainer.optimizer.zero_grad(set_to_none=True)
            for module, training in mode_snapshot:
                module.training = training
            if ddp_model is not None and ddp_training_snapshot is not None:
                ddp_model.training = ddp_training_snapshot

        if completed_blocks != local_steps:
            raise RuntimeError(
                "module-rescue did not complete the planned local blocks: "
                f"{completed_blocks} != {local_steps}"
            )
        if int(world_size) > 1:
            if not dist.is_available() or not dist.is_initialized():
                raise RuntimeError(
                    "module-rescue draw audit requires initialized DDP"
                )
            dist.all_reduce(celltype_draw_count, op=dist.ReduceOp.SUM)
            if self.lineage_balance_plan.enabled:
                dist.all_reduce(lineage_update_count, op=dist.ReduceOp.SUM)
        expected_draws = int(self.config.draws_per_celltype_per_epoch)
        if expected_draws > 0 and not bool(
            (celltype_draw_count == expected_draws).all()
        ):
            raise RuntimeError(
                "module-rescue celltype exposure contract failed: expected "
                f"{expected_draws} draws, observed "
                f"{celltype_draw_count.detach().cpu().tolist()}."
            )
        if self.lineage_balance_plan.enabled:
            expected_lineage_count = torch.tensor(
                [
                    self.lineage_balance_plan.neuronal_blocks_per_update,
                    self.lineage_balance_plan.non_neuronal_blocks_per_update,
                ],
                dtype=torch.long,
                device=trainer.device,
            )
            if not bool(
                (lineage_update_count == expected_lineage_count).all()
            ):
                raise RuntimeError(
                    "module-rescue lineage exposure contract failed: "
                    f"expected={expected_lineage_count.detach().cpu().tolist()}, "
                    f"observed={lineage_update_count.detach().cpu().tolist()}"
                )
        if str(self.config.contract_version) == "v2" and (
            applied_updates != planned_updates or skipped_updates != 0
        ):
            raise RuntimeError(
                "module-rescue v2 update contract failed: "
                f"applied={applied_updates}, planned={planned_updates}, "
                f"skipped={skipped_updates}"
            )

        block_denominator = max(completed_blocks, 1)
        update_denominator = max(applied_updates, 1)
        if self.lineage_balance_plan.enabled:
            neuronal_blocks = float(lineage_update_count[:, 0].sum().cpu())
            non_neuronal_blocks = float(
                lineage_update_count[:, 1].sum().cpu()
            )
            raw_neuronal_weight = float(
                self.lineage_balance_plan.celltype_raw_weights[
                    self.lineage_balance_plan.neuronal_celltype_ids[0]
                ]
            )
            raw_non_neuronal_weight = float(
                self.lineage_balance_plan.celltype_raw_weights[
                    self.lineage_balance_plan.non_neuronal_celltype_ids[0]
                ]
            )
            normalized_neuronal_weight = float(
                self.lineage_balance_plan.celltype_normalized_weights[
                    self.lineage_balance_plan.neuronal_celltype_ids[0]
                ]
            )
            normalized_non_neuronal_weight = float(
                self.lineage_balance_plan.celltype_normalized_weights[
                    self.lineage_balance_plan.non_neuronal_celltype_ids[0]
                ]
            )
            neuronal_mass = neuronal_blocks * normalized_neuronal_weight
            non_neuronal_mass = (
                non_neuronal_blocks * normalized_non_neuronal_weight
            )
            total_lineage_mass = neuronal_mass + non_neuronal_mass
            neuronal_mass_fraction = neuronal_mass / total_lineage_mass
            non_neuronal_mass_fraction = (
                non_neuronal_mass / total_lineage_mass
            )
            neuronal_per_update_min = float(
                lineage_update_count[:, 0].min().cpu()
            )
            neuronal_per_update_max = float(
                lineage_update_count[:, 0].max().cpu()
            )
            non_neuronal_per_update_min = float(
                lineage_update_count[:, 1].min().cpu()
            )
            non_neuronal_per_update_max = float(
                lineage_update_count[:, 1].max().cpu()
            )
        else:
            neuronal_blocks = 0.0
            non_neuronal_blocks = 0.0
            raw_neuronal_weight = 1.0
            raw_non_neuronal_weight = 1.0
            normalized_neuronal_weight = 1.0
            normalized_non_neuronal_weight = 1.0
            neuronal_mass_fraction = 0.0
            non_neuronal_mass_fraction = 0.0
            neuronal_per_update_min = 0.0
            neuronal_per_update_max = 0.0
            non_neuronal_per_update_min = 0.0
            non_neuronal_per_update_max = 0.0
        return {
            "loss/prism_module_rescue": sums["loss"] / block_denominator,
            "loss/prism_module_rescue_objective_unweighted": sums["unweighted_loss"] / block_denominator,
            "metric/prism_module_rescue_unweighted": sums["rescue"] / block_denominator,
            "loss/prism_module_rescue_huber_unweighted": sums["huber"] / block_denominator,
            "loss/prism_module_rescue_ccc_unweighted": sums["ccc"] / block_denominator,
            "loss/prism_module_rescue_branch": sums["branch"] / block_denominator,
            "loss/prism_module_rescue_full": sums["full"] / block_denominator,
            "loss/prism_module_rescue_branch_huber": sums["branch_huber"] / block_denominator,
            "loss/prism_module_rescue_full_huber": sums["full_huber"] / block_denominator,
            "loss/prism_module_rescue_branch_ccc": sums["branch_ccc"] / block_denominator,
            "loss/prism_module_rescue_full_ccc": sums["full_ccc"] / block_denominator,
            # Historical keys retain their original module-wise standardized
            # CCC meaning so old and new logs remain comparable.
            "metric/prism_module_rescue_branch_ccc": sums["branch_module_ccc"] / block_denominator,
            "metric/prism_module_rescue_full_ccc": sums["full_module_ccc"] / block_denominator,
            "metric/prism_module_rescue_branch_flattened_ccc": sums["branch_mean_ccc"] / block_denominator,
            "metric/prism_module_rescue_full_flattened_ccc": sums["full_mean_ccc"] / block_denominator,
            "metric/prism_module_rescue_branch_pearson": sums["branch_corr"] / block_denominator,
            "metric/prism_module_rescue_full_pearson": sums["full_corr"] / block_denominator,
            "metric/prism_module_rescue_groups": sums["donor_groups"] / block_denominator,
            "metric/prism_module_rescue_grad_norm": sums["grad_norm"] / update_denominator,
            "metric/prism_module_rescue_allowed_grad_tensors": sums["allowed_grad_tensors"] / update_denominator,
            "metric/prism_module_rescue_cleared_grad_tensors": sums["cleared_grad_tensors"] / update_denominator,
            "metric/prism_module_rescue_sampled_cells_per_block": sums["sampled_cells"] / block_denominator,
            "metric/prism_module_rescue_celltype_draw_min": float(celltype_draw_count.min().detach().cpu()),
            "metric/prism_module_rescue_celltype_draw_max": float(celltype_draw_count.max().detach().cpu()),
            "metric/prism_module_rescue_blocks": float(completed_blocks),
            "metric/prism_module_rescue_updates_attempted": float(attempted_updates),
            "metric/prism_module_rescue_updates_applied": float(applied_updates),
            "metric/prism_module_rescue_updates_skipped": float(skipped_updates),
            "metric/prism_module_rescue_global_blocks_per_update": float(local_per_update * world_size),
            "metric/prism_module_rescue_lineage_balance_enabled": float(
                self.lineage_balance_plan.enabled
            ),
            "metric/prism_module_rescue_neuronal_blocks": neuronal_blocks,
            "metric/prism_module_rescue_non_neuronal_blocks": non_neuronal_blocks,
            "metric/prism_module_rescue_neuronal_blocks_per_update_min": neuronal_per_update_min,
            "metric/prism_module_rescue_neuronal_blocks_per_update_max": neuronal_per_update_max,
            "metric/prism_module_rescue_non_neuronal_blocks_per_update_min": non_neuronal_per_update_min,
            "metric/prism_module_rescue_non_neuronal_blocks_per_update_max": non_neuronal_per_update_max,
            "metric/prism_module_rescue_neuronal_objective_mass_fraction": neuronal_mass_fraction,
            "metric/prism_module_rescue_non_neuronal_objective_mass_fraction": non_neuronal_mass_fraction,
            "weight/prism_module_rescue_neuronal_raw": raw_neuronal_weight,
            "weight/prism_module_rescue_non_neuronal_raw": raw_non_neuronal_weight,
            "weight/prism_module_rescue_neuronal_normalized": normalized_neuronal_weight,
            "weight/prism_module_rescue_non_neuronal_normalized": normalized_non_neuronal_weight,
            "weight/prism_module_rescue_ramp": float(warmup),
            "weight/prism_module_rescue_ccc": float(self.config.ccc_weight),
            # Historical alias now means applied optimiser updates, not blocks.
            "metric/prism_module_rescue_steps": float(applied_updates),
        }


def write_module_rescue_manifest(
    path: str | Path,
    *,
    config: PrismModuleRescueTrainingConfig,
    updater: PrismModuleRescueUpdater,
    world_size: int | None = None,
) -> None:
    allowed_names = (
        sorted(updater._parameter_allowlist.values())
        if updater._parameter_allowlist is not None
        else []
    )
    allowlist_sha256 = hashlib.sha256(
        "\n".join(allowed_names).encode("utf-8")
    ).hexdigest()
    payload = {
        "schema_version": "kmlee_bam.prism_module_rescue_training.v2",
        "config": {
            key: value
            for key, value in vars(config).items()
        },
        "weight_donor_ids": list(updater.donor_ids),
        "eligible_celltypes": list(updater.sampler.eligible_celltypes),
        "lineage_balance": updater.lineage_balance_manifest(),
        "module_count": len(updater.artifact.module_names),
        "target_source": str(config.stats_path),
        "target_artifact_sha256": updater.artifact.artifact_sha256,
        "target_artifact_schema": updater.artifact.schema_version,
        "target_source_split": updater.artifact.metadata.get("source_split"),
        "target_activity_normalization": updater.artifact.metadata.get(
            "activity_normalization"
        ),
        "ccc_valid_celltype_module_count": int(
            updater.artifact.ccc_valid_mask.sum()
        ),
        "candidate_cells_min": updater.candidate_size_min,
        "candidate_cells_max": updater.candidate_size_max,
        "nonoverlap_reuse_fraction": 0.0 if bool(
            config.nonoverlap_sampling
        ) else None,
        "ordinary_random_minibatch_centering_used": False,
        "projection_statistics_updated_by_rescue": not bool(
            config.deterministic_forward
        ),
        "rescue_parameter_scope": (
            (
                "decoder_state_precision_and_joint_gate_only"
                if "generator_count_gate.log_alpha" in allowed_names
                else "decoder_state_and_precision_only"
            )
            if bool(config.restrict_parameter_updates)
            else "all_gradient_reachable_parameters"
        ),
        "rescue_parameter_names": allowed_names,
        "rescue_parameter_names_sha256": allowlist_sha256,
    }
    if world_size is not None:
        payload["resolved_schedule"] = {
            "world_size": int(world_size),
            "local_blocks_per_epoch": updater.steps_per_rank(
                world_size=int(world_size)
            ),
            "local_blocks_per_optimizer_update": (
                updater.local_blocks_per_optimizer_update(
                    world_size=int(world_size)
                )
            ),
            "global_blocks_per_optimizer_update": (
                updater.local_blocks_per_optimizer_update(
                    world_size=int(world_size)
                )
                * int(world_size)
            ),
            "optimizer_updates_per_epoch": (
                updater.optimizer_updates_per_epoch(
                    world_size=int(world_size)
                )
            ),
        }
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        try:
            existing = json.loads(destination.read_text(encoding="utf-8"))
        except Exception as exc:
            raise ValueError(
                f"existing module-rescue manifest is unreadable: {destination}"
            ) from exc
        if existing != payload:
            raise ValueError(
                "refusing to overwrite an incompatible module-rescue manifest: "
                f"{destination}"
            )
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
