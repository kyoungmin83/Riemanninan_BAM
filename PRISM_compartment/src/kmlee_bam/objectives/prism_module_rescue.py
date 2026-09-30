"""Pure mathematical utilities for the PRISM v2 module-rescue objective.

The module-rescue target is deliberately *not* the PRISM module coefficient.
Predicted ordinal probabilities and observed ordinal tiers are first mapped to
the same gene-tier space and projected with the same fixed activity dictionary:

    probabilities -> expected ordinal tier -> donor x cell-type gene mean
                  -> fixed W_act projection
                  -> separately donor-centered prediction and observation
                  -> cell-type train-scale standardization -> weighted Huber

This module contains no model, trainer, DDP, artifact-loading, or checkpoint
integration. It estimates only the two separate donor centers from one
complete, sampler-defined contrast block. Callers must supply the frozen
train-only scale and any precomputed donor/context weights.

The high-level objective requires a grouped block containing the complete set
of donor groups used to define each cell-type center.  Centering independently
inside ordinary minibatches or DDP ranks is mathematically wrong.  A future
trainer integration must either construct complete grouped blocks or globally
synchronize group sums before calling the centering step.

See ``doc/prism_v2_module_rescue_generator_budget_design_20260729.md`` §6.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Optional

import torch
import torch.nn.functional as F


@dataclass(frozen=True)
class PrismModuleRescueConfig:
    """Configuration for the dual-view PRISM module-rescue objective.

    ``enabled=False`` is the compatibility default.  The high-level loss
    function returns scalar zeros without inspecting any optional rescue
    inputs in that mode.
    """

    enabled: bool = False
    lambda_module: float = 0.0
    branch_fraction: float = 2.0 / 3.0
    huber_delta: float = 1.0
    ccc_weight: float = 0.0
    eps: float = 1e-8
    validate_probabilities: bool = True
    probability_atol: float = 5e-3


@dataclass(frozen=True)
class WeightedHuberOutput:
    """A weighted mean together with DDP-safe sufficient statistics."""

    loss: torch.Tensor
    numerator: torch.Tensor
    denominator: torch.Tensor


@dataclass(frozen=True)
class GroupedModuleActivity:
    """Differentiable donor-context means with deterministic sorted group ids."""

    mean: torch.Tensor
    group_id: torch.Tensor
    celltype_id: torch.Tensor
    cell_count: torch.Tensor


@dataclass(frozen=True)
class CorrelationCCCSufficientStats:
    """Additive, module-wise weighted statistics for correlation and CCC."""

    weight_sum: torch.Tensor
    observation_count: torch.Tensor
    prediction_sum: torch.Tensor
    target_sum: torch.Tensor
    prediction_square_sum: torch.Tensor
    target_square_sum: torch.Tensor
    cross_sum: torch.Tensor


@dataclass(frozen=True)
class CorrelationCCCOutput:
    """Finite module-wise metrics plus an authoritative validity mask."""

    correlation: torch.Tensor
    ccc: torch.Tensor
    valid: torch.Tensor


@dataclass(frozen=True)
class DifferentiableCCCOutput:
    """Differentiable celltype-level concordance loss and diagnostics."""

    loss: torch.Tensor
    mean_ccc: torch.Tensor
    valid_count: torch.Tensor


@dataclass(frozen=True)
class CenteredModuleViewOutput:
    """One prediction view evaluated on donor×celltype pseudobulk rows."""

    huber: WeightedHuberOutput
    prediction: torch.Tensor
    target: torch.Tensor
    group_id: torch.Tensor
    group_celltype_id: torch.Tensor
    group_cell_count: torch.Tensor
    group_weight: torch.Tensor
    valid_group: torch.Tensor
    celltype_group_count: torch.Tensor
    correlation_stats: CorrelationCCCSufficientStats
    metrics: CorrelationCCCOutput
    celltype_correlation_stats: CorrelationCCCSufficientStats
    celltype_metrics: CorrelationCCCOutput
    concordance: DifferentiableCCCOutput


@dataclass(frozen=True)
class PrismModuleRescueOutput:
    """Dual-view loss values and the standardized activities used to form them."""

    loss: torch.Tensor
    rescue_loss: torch.Tensor
    branch_loss: torch.Tensor
    full_loss: torch.Tensor
    branch_numerator: torch.Tensor
    branch_denominator: torch.Tensor
    full_numerator: torch.Tensor
    full_denominator: torch.Tensor
    huber_rescue_loss: torch.Tensor
    ccc_rescue_loss: torch.Tensor
    branch_huber_loss: torch.Tensor
    full_huber_loss: torch.Tensor
    branch_ccc_loss: torch.Tensor
    full_ccc_loss: torch.Tensor
    branch_mean_ccc: torch.Tensor
    full_mean_ccc: torch.Tensor
    branch_activity: Optional[torch.Tensor] = None
    full_activity: Optional[torch.Tensor] = None
    target_activity: Optional[torch.Tensor] = None
    branch_view: Optional[CenteredModuleViewOutput] = None
    full_view: Optional[CenteredModuleViewOutput] = None
    pathology_axis_loss: Optional[torch.Tensor] = None
    pathology_axis_branch_loss: Optional[torch.Tensor] = None
    pathology_axis_full_loss: Optional[torch.Tensor] = None
    pathology_axis_mean_correlation: Optional[torch.Tensor] = None


def _shape(tensor: torch.Tensor) -> tuple[int, ...]:
    return tuple(int(x) for x in tensor.shape)


def _require_tensor(name: str, tensor: object) -> torch.Tensor:
    if not isinstance(tensor, torch.Tensor):
        raise TypeError(f"{name} must be a torch.Tensor, got {type(tensor).__name__}.")
    return tensor


def _require_finite(name: str, tensor: torch.Tensor) -> None:
    if not bool(torch.isfinite(tensor).all()):
        raise ValueError(f"{name} contains non-finite values.")


def _stable_float(tensor: torch.Tensor) -> torch.Tensor:
    """Use float32 for integer/half inputs while preserving float64 tests."""

    if not (tensor.is_floating_point() or tensor.is_complex()):
        return tensor.to(dtype=torch.float32)
    if tensor.dtype in (torch.float16, torch.bfloat16):
        return tensor.to(dtype=torch.float32)
    return tensor


def _zero(reference: Optional[torch.Tensor]) -> torch.Tensor:
    if isinstance(reference, torch.Tensor):
        dtype = reference.dtype if reference.is_floating_point() else torch.float32
        return torch.zeros((), dtype=dtype, device=reference.device)
    return torch.zeros((), dtype=torch.float32)


def expected_ordinal_tier(
    probabilities: torch.Tensor,
    *,
    validate: bool = True,
    atol: float = 5e-3,
) -> torch.Tensor:
    """Return :math:`E[c]` from ordinal probabilities of shape ``[B, G, C]``.

    The bin coordinates are exactly ``0, ..., C-1``.  Half-precision inputs
    are accumulated in float32 so summing many genes in the later activity
    projection does not begin with avoidable tier-rounding error.
    """

    probabilities = _require_tensor("probabilities", probabilities)
    if probabilities.ndim != 3:
        raise ValueError(
            "probabilities must have shape [B, G, C], "
            f"got {_shape(probabilities)}."
        )
    if int(probabilities.shape[0]) <= 0 or int(probabilities.shape[1]) <= 0:
        raise ValueError("probabilities must contain at least one cell and one gene.")
    if int(probabilities.shape[2]) < 2:
        raise ValueError("probabilities must contain at least two ordinal bins.")
    if not probabilities.is_floating_point():
        raise TypeError(f"probabilities must be floating point, got {probabilities.dtype}.")
    _require_finite("probabilities", probabilities)
    if float(atol) < 0.0:
        raise ValueError(f"atol must be non-negative, got {atol}.")

    probs = _stable_float(probabilities)
    if validate:
        tolerance = float(atol)
        if bool((probs < -tolerance).any()) or bool((probs > 1.0 + tolerance).any()):
            raise ValueError("probabilities must lie in [0, 1] within atol.")
        total = probs.sum(dim=-1)
        if not bool(
            torch.allclose(
                total,
                torch.ones_like(total),
                rtol=0.0,
                atol=tolerance,
            )
        ):
            max_error = float((total - 1.0).abs().max().detach().cpu())
            raise ValueError(
                "probabilities must sum to one across ordinal bins; "
                f"maximum absolute error is {max_error:.6g}."
            )

    tiers = torch.arange(
        int(probs.shape[-1]),
        dtype=probs.dtype,
        device=probs.device,
    )
    expected = (probs * tiers).sum(dim=-1)
    _require_finite("expected ordinal tier", expected)
    return expected


def project_fixed_module_activity(
    gene_tier: torch.Tensor,
    activity_weight: torch.Tensor,
) -> torch.Tensor:
    """Project ``[B, G]`` gene tiers with a fixed ``W_act[M, G]`` dictionary.

    ``activity_weight`` is detached inside this function by design.  Gradients
    flow to predicted gene tiers, never into the biological module registry.
    """

    gene_tier = _require_tensor("gene_tier", gene_tier)
    activity_weight = _require_tensor("activity_weight", activity_weight)
    if gene_tier.ndim != 2:
        raise ValueError(
            f"gene_tier must have shape [B, G], got {_shape(gene_tier)}."
        )
    if activity_weight.ndim != 2:
        raise ValueError(
            "activity_weight must have shape [M, G], "
            f"got {_shape(activity_weight)}."
        )
    if int(gene_tier.shape[0]) <= 0 or int(gene_tier.shape[1]) <= 0:
        raise ValueError("gene_tier must contain at least one cell and one gene.")
    if int(activity_weight.shape[0]) <= 0:
        raise ValueError("activity_weight must contain at least one module.")
    if int(gene_tier.shape[1]) != int(activity_weight.shape[1]):
        raise ValueError(
            "gene dimension mismatch: "
            f"gene_tier has G={int(gene_tier.shape[1])}, "
            f"activity_weight has G={int(activity_weight.shape[1])}."
        )
    _require_finite("gene_tier", gene_tier)
    _require_finite("activity_weight", activity_weight)

    values = _stable_float(gene_tier)
    weight = activity_weight.detach().to(device=values.device, dtype=values.dtype)
    activity = F.linear(values, weight)
    _require_finite("projected module activity", activity)
    return activity


def standardize_module_activity(
    activity: torch.Tensor,
    celltype_id: torch.Tensor,
    train_center: torch.Tensor,
    train_scale: torch.Tensor,
) -> torch.Tensor:
    """Apply frozen cell-type train statistics to activity ``[B, M]``.

    ``train_center`` and ``train_scale`` must both have shape ``[T, M]`` and
    must have been fitted using training donors only.  This function cannot
    infer or certify the artifact split, but it detaches both tensors so the
    loss cannot alter them.
    """

    activity = _require_tensor("activity", activity)
    celltype_id = _require_tensor("celltype_id", celltype_id)
    train_center = _require_tensor("train_center", train_center)
    train_scale = _require_tensor("train_scale", train_scale)
    if activity.ndim != 2:
        raise ValueError(f"activity must have shape [B, M], got {_shape(activity)}.")
    if celltype_id.ndim != 1:
        raise ValueError(
            f"celltype_id must have shape [B], got {_shape(celltype_id)}."
        )
    if train_center.ndim != 2 or train_scale.ndim != 2:
        raise ValueError(
            "train_center and train_scale must both have shape [T, M], "
            f"got {_shape(train_center)} and {_shape(train_scale)}."
        )
    if _shape(train_center) != _shape(train_scale):
        raise ValueError(
            "train_center and train_scale shapes differ: "
            f"{_shape(train_center)} vs {_shape(train_scale)}."
        )
    batch, modules = int(activity.shape[0]), int(activity.shape[1])
    if batch <= 0 or modules <= 0:
        raise ValueError("activity must contain at least one cell and one module.")
    if int(celltype_id.numel()) != batch:
        raise ValueError(
            f"celltype_id length {int(celltype_id.numel())} != batch size {batch}."
        )
    if int(train_center.shape[0]) <= 0 or int(train_center.shape[1]) != modules:
        raise ValueError(
            "train statistic module dimension mismatch: "
            f"activity has M={modules}, train_center has {_shape(train_center)}."
        )
    if celltype_id.dtype not in (
        torch.int8,
        torch.int16,
        torch.int32,
        torch.int64,
        torch.uint8,
    ):
        raise TypeError(f"celltype_id must have integer dtype, got {celltype_id.dtype}.")
    _require_finite("activity", activity)
    _require_finite("train_center", train_center)
    _require_finite("train_scale", train_scale)
    if not bool((train_scale > 0).all()):
        raise ValueError("train_scale must be strictly positive for every celltype/module.")

    ids = celltype_id.to(device=activity.device, dtype=torch.long)
    n_celltypes = int(train_center.shape[0])
    if bool((ids < 0).any()) or bool((ids >= n_celltypes).any()):
        low = int(ids.min().detach().cpu())
        high = int(ids.max().detach().cpu())
        raise ValueError(
            f"celltype_id values must be in [0, {n_celltypes}), got range [{low}, {high}]."
        )

    values = _stable_float(activity)
    center_table = train_center.detach().to(device=values.device, dtype=values.dtype)
    scale_table = train_scale.detach().to(device=values.device, dtype=values.dtype)
    center = center_table.index_select(0, ids)
    scale = scale_table.index_select(0, ids)
    standardized = (values - center) / scale
    _require_finite("standardized module activity", standardized)
    return standardized


def _integer_id_vector(
    name: str,
    value: torch.Tensor,
    *,
    batch: int,
    device: torch.device,
) -> torch.Tensor:
    value = _require_tensor(name, value)
    if value.ndim != 1 or int(value.numel()) != int(batch):
        raise ValueError(f"{name} must have shape [{batch}], got {_shape(value)}.")
    if value.dtype not in (
        torch.int8,
        torch.int16,
        torch.int32,
        torch.int64,
        torch.uint8,
    ):
        raise TypeError(f"{name} must have integer dtype, got {value.dtype}.")
    result = value.to(device=device, dtype=torch.long)
    if bool((result < 0).any()):
        raise ValueError(f"{name} must be non-negative.")
    return result


def aggregate_donor_context_mean(
    activity: torch.Tensor,
    group_id: torch.Tensor,
    celltype_id: torch.Tensor,
    *,
    cell_valid_mask: Optional[torch.Tensor] = None,
) -> GroupedModuleActivity:
    """Differentiably aggregate cells to sorted donor-context group means.

    A group id must identify exactly one cell type.  For the final PRISM v2
    endpoint it should encode one ``(donor, celltype)`` pair, after any
    pre-specified region pooling.  ``index_add`` keeps gradients from every
    group mean back to its contributing predicted cells.
    """

    activity = _require_tensor("activity", activity)
    if activity.ndim != 2:
        raise ValueError(f"activity must have shape [B, M], got {_shape(activity)}.")
    batch, modules = int(activity.shape[0]), int(activity.shape[1])
    if batch <= 0 or modules <= 0:
        raise ValueError("activity must contain at least one cell and one module.")
    _require_finite("activity", activity)
    values = _stable_float(activity)
    groups = _integer_id_vector(
        "group_id",
        group_id,
        batch=batch,
        device=values.device,
    )
    celltypes = _integer_id_vector(
        "celltype_id",
        celltype_id,
        batch=batch,
        device=values.device,
    )

    if cell_valid_mask is None:
        cell_valid = torch.ones(batch, dtype=torch.bool, device=values.device)
    else:
        cell_valid_mask = _require_tensor("cell_valid_mask", cell_valid_mask)
        if cell_valid_mask.ndim != 1 or int(cell_valid_mask.numel()) != batch:
            raise ValueError(
                f"cell_valid_mask must have shape [{batch}], "
                f"got {_shape(cell_valid_mask)}."
            )
        if cell_valid_mask.dtype != torch.bool:
            raise TypeError(
                f"cell_valid_mask must be bool, got {cell_valid_mask.dtype}."
            )
        cell_valid = cell_valid_mask.to(device=values.device)

    unique_group, inverse = torch.unique(
        groups,
        sorted=True,
        return_inverse=True,
    )
    n_groups = int(unique_group.numel())
    group_celltype = torch.zeros(
        n_groups,
        dtype=torch.long,
        device=values.device,
    )
    group_celltype.scatter_(0, inverse, celltypes)
    if not bool(
        torch.equal(
            group_celltype.index_select(0, inverse),
            celltypes,
        )
    ):
        raise ValueError("each group_id must map to exactly one celltype_id.")

    valid_float = cell_valid.to(dtype=values.dtype)
    cell_count = torch.zeros(
        n_groups,
        dtype=torch.long,
        device=values.device,
    ).index_add_(0, inverse, cell_valid.to(dtype=torch.long))
    activity_sum = values.new_zeros(n_groups, modules).index_add_(
        0,
        inverse,
        values * valid_float[:, None],
    )
    mean = activity_sum / cell_count.clamp_min(1).to(values.dtype)[:, None]
    _require_finite("donor-context activity mean", mean)
    return GroupedModuleActivity(
        mean=mean,
        group_id=unique_group,
        celltype_id=group_celltype,
        cell_count=cell_count,
    )


def grouped_expected_module_activity(
    probabilities: torch.Tensor,
    activity_weight: torch.Tensor,
    group_id: torch.Tensor,
    celltype_id: torch.Tensor,
    *,
    cell_valid_mask: Optional[torch.Tensor] = None,
    validate_probabilities: bool = True,
    probability_atol: float = 5e-3,
) -> GroupedModuleActivity:
    """Map probabilities to donor-context module means efficiently.

    Expected gene tiers are averaged within groups *before* applying
    ``W_act``.  The projection is linear, so this is exactly equivalent to
    projecting every cell and averaging, while avoiding a potentially huge
    ``[n_cells, M]`` matrix multiplication.
    """

    expected_tier = expected_ordinal_tier(
        probabilities,
        validate=bool(validate_probabilities),
        atol=float(probability_atol),
    )
    grouped_tier = aggregate_donor_context_mean(
        expected_tier,
        group_id,
        celltype_id,
        cell_valid_mask=cell_valid_mask,
    )
    grouped_activity = project_fixed_module_activity(
        grouped_tier.mean,
        activity_weight,
    )
    return GroupedModuleActivity(
        mean=grouped_activity,
        group_id=grouped_tier.group_id,
        celltype_id=grouped_tier.celltype_id,
        cell_count=grouped_tier.cell_count,
    )


def _expand_group_rows_to_cells(
    grouped: GroupedModuleActivity,
    cell_group_id: torch.Tensor,
) -> torch.Tensor:
    """Differentiably repeat sorted group rows for their source cells."""

    if not isinstance(grouped, GroupedModuleActivity):
        raise TypeError("grouped must be GroupedModuleActivity.")
    cells = int(cell_group_id.numel())
    ids = _integer_id_vector(
        "group_id",
        cell_group_id,
        batch=cells,
        device=grouped.mean.device,
    )
    position = torch.searchsorted(grouped.group_id, ids)
    safe_position = position.clamp(max=max(0, int(grouped.group_id.numel()) - 1))
    found = (position < int(grouped.group_id.numel())) & (
        grouped.group_id.index_select(0, safe_position) == ids
    )
    if not bool(found.all()):
        raise RuntimeError("grouped activity does not cover every source cell group.")
    return grouped.mean.index_select(0, position)


def donor_context_weighted_huber(
    prediction: torch.Tensor,
    target: torch.Tensor,
    *,
    donor_weight: Optional[torch.Tensor] = None,
    context_weight: Optional[torch.Tensor] = None,
    valid_mask: Optional[torch.Tensor] = None,
    delta: float = 1.0,
    eps: float = 1e-8,
) -> WeightedHuberOutput:
    """Compute elementwise Huber averaged with frozen donor/context weights.

    Inputs have shape ``[N, M]``.  In the final module-rescue endpoint, ``N``
    must index donor×celltype aggregates, never raw cells.  ``donor_weight`` and
    ``context_weight`` are optional ``[N]`` vectors and are multiplied.  They
    must be computed from the complete training split (or supplied by a
    donor-context grouped loader), not estimated from the minibatch.

    The returned numerator and denominator are additive sufficient statistics:
    sum them across batches and DDP ranks before division for an unbiased
    epoch-level value.
    """

    prediction = _require_tensor("prediction", prediction)
    target = _require_tensor("target", target)
    if prediction.ndim != 2 or target.ndim != 2:
        raise ValueError(
            "prediction and target must both have shape [B, M], "
            f"got {_shape(prediction)} and {_shape(target)}."
        )
    if _shape(prediction) != _shape(target):
        raise ValueError(
            f"prediction and target shapes differ: {_shape(prediction)} vs {_shape(target)}."
        )
    if int(prediction.shape[0]) <= 0 or int(prediction.shape[1]) <= 0:
        raise ValueError("prediction must contain at least one row and one module.")
    if not prediction.is_floating_point():
        raise TypeError(f"prediction must be floating point, got {prediction.dtype}.")
    if float(delta) <= 0.0:
        raise ValueError(f"delta must be positive, got {delta}.")
    if float(eps) <= 0.0:
        raise ValueError(f"eps must be positive, got {eps}.")
    _require_finite("prediction", prediction)
    _require_finite("target", target)

    pred = _stable_float(prediction)
    fixed_target = target.detach().to(device=pred.device, dtype=pred.dtype)
    element_loss = F.huber_loss(
        pred,
        fixed_target,
        reduction="none",
        delta=float(delta),
    )
    batch, modules = int(pred.shape[0]), int(pred.shape[1])
    row_weight = pred.new_ones(batch)

    for name, supplied in (
        ("donor_weight", donor_weight),
        ("context_weight", context_weight),
    ):
        if supplied is None:
            continue
        supplied = _require_tensor(name, supplied)
        if supplied.ndim != 1 or int(supplied.numel()) != batch:
            raise ValueError(
                f"{name} must have shape [{batch}], got {_shape(supplied)}."
            )
        _require_finite(name, supplied)
        fixed_weight = supplied.detach().to(device=pred.device, dtype=pred.dtype)
        if bool((fixed_weight < 0).any()):
            raise ValueError(f"{name} must be non-negative.")
        row_weight = row_weight * fixed_weight

    element_weight = row_weight[:, None].expand(batch, modules)
    if valid_mask is not None:
        valid_mask = _require_tensor("valid_mask", valid_mask)
        if valid_mask.ndim == 1 and int(valid_mask.numel()) == batch:
            mask = valid_mask[:, None].expand(batch, modules)
        elif _shape(valid_mask) == (batch, modules):
            mask = valid_mask
        else:
            raise ValueError(
                "valid_mask must have shape [B] or [B, M], "
                f"got {_shape(valid_mask)} for B={batch}, M={modules}."
            )
        if mask.dtype != torch.bool:
            raise TypeError(f"valid_mask must be bool, got {mask.dtype}.")
        element_weight = element_weight * mask.to(
            device=pred.device,
            dtype=pred.dtype,
        )

    denominator = element_weight.sum()
    if not bool(denominator > 0):
        raise ValueError("weighted Huber denominator is zero; no valid positive-weight data.")
    numerator = (element_loss * element_weight).sum()
    loss = numerator / denominator.clamp_min(float(eps))
    _require_finite("weighted Huber numerator", numerator)
    _require_finite("weighted Huber denominator", denominator)
    _require_finite("weighted Huber loss", loss)
    return WeightedHuberOutput(
        loss=loss,
        numerator=numerator,
        denominator=denominator,
    )


def correlation_ccc_sufficient_statistics(
    prediction: torch.Tensor,
    target: torch.Tensor,
    *,
    row_weight: Optional[torch.Tensor] = None,
    valid_mask: Optional[torch.Tensor] = None,
) -> CorrelationCCCSufficientStats:
    """Return additive weighted statistics for module-wise correlation/CCC.

    Both inputs must already be donor-centered and train-scale standardized.
    All returned tensors are detached because these are validation/reporting
    statistics, not an additional differentiable objective.
    """

    prediction = _require_tensor("prediction", prediction)
    target = _require_tensor("target", target)
    if prediction.ndim != 2 or target.ndim != 2:
        raise ValueError(
            "prediction and target must both have shape [N, M], "
            f"got {_shape(prediction)} and {_shape(target)}."
        )
    if _shape(prediction) != _shape(target):
        raise ValueError(
            f"prediction and target shapes differ: {_shape(prediction)} vs {_shape(target)}."
        )
    rows, modules = int(prediction.shape[0]), int(prediction.shape[1])
    if rows <= 0 or modules <= 0:
        raise ValueError("prediction must contain at least one row and one module.")
    _require_finite("prediction", prediction)
    _require_finite("target", target)
    pred = _stable_float(prediction).detach()
    fixed_target = target.detach().to(device=pred.device, dtype=pred.dtype)

    if row_weight is None:
        weight = pred.new_ones(rows)
    else:
        row_weight = _require_tensor("row_weight", row_weight)
        if row_weight.ndim != 1 or int(row_weight.numel()) != rows:
            raise ValueError(
                f"row_weight must have shape [{rows}], got {_shape(row_weight)}."
            )
        _require_finite("row_weight", row_weight)
        weight = row_weight.detach().to(device=pred.device, dtype=pred.dtype)
        if bool((weight < 0).any()):
            raise ValueError("row_weight must be non-negative.")

    if valid_mask is None:
        mask = torch.ones(rows, modules, dtype=torch.bool, device=pred.device)
    else:
        valid_mask = _require_tensor("valid_mask", valid_mask)
        if valid_mask.ndim == 1 and int(valid_mask.numel()) == rows:
            mask = valid_mask[:, None].expand(rows, modules)
        elif _shape(valid_mask) == (rows, modules):
            mask = valid_mask
        else:
            raise ValueError(
                "valid_mask must have shape [N] or [N, M], "
                f"got {_shape(valid_mask)} for N={rows}, M={modules}."
            )
        if mask.dtype != torch.bool:
            raise TypeError(f"valid_mask must be bool, got {mask.dtype}.")
        mask = mask.to(device=pred.device)

    positive = weight > 0
    effective_mask = mask & positive[:, None]
    element_weight = weight[:, None] * effective_mask.to(dtype=pred.dtype)
    stats = CorrelationCCCSufficientStats(
        weight_sum=element_weight.sum(dim=0),
        observation_count=effective_mask.sum(dim=0).to(dtype=pred.dtype),
        prediction_sum=(element_weight * pred).sum(dim=0),
        target_sum=(element_weight * fixed_target).sum(dim=0),
        prediction_square_sum=(element_weight * pred.square()).sum(dim=0),
        target_square_sum=(element_weight * fixed_target.square()).sum(dim=0),
        cross_sum=(element_weight * pred * fixed_target).sum(dim=0),
    )
    for name, value in (
        ("weight_sum", stats.weight_sum),
        ("observation_count", stats.observation_count),
        ("prediction_sum", stats.prediction_sum),
        ("target_sum", stats.target_sum),
        ("prediction_square_sum", stats.prediction_square_sum),
        ("target_square_sum", stats.target_square_sum),
        ("cross_sum", stats.cross_sum),
    ):
        _require_finite(f"correlation statistic {name}", value)
    return stats


def celltype_correlation_ccc_sufficient_statistics(
    prediction: torch.Tensor,
    target: torch.Tensor,
    celltype_id: torch.Tensor,
    *,
    n_celltypes: int,
    row_weight: Optional[torch.Tensor] = None,
    valid_mask: Optional[torch.Tensor] = None,
) -> CorrelationCCCSufficientStats:
    """Return per-celltype stats after flattening donor×module elements.

    This is the additive Pearson/CCC analogue of the canonical post-hoc
    ``blind_centered`` endpoint, which computes a Spearman correlation after
    flattening the centered donor×module matrix separately for each cell type.
    Spearman ranks are not additive sufficient statistics, so the exact
    canonical Spearman remains a post-hoc metric.
    """

    prediction = _require_tensor("prediction", prediction)
    target = _require_tensor("target", target)
    if prediction.ndim != 2 or _shape(prediction) != _shape(target):
        raise ValueError(
            "prediction and target must have the same shape [N, M], "
            f"got {_shape(prediction)} and {_shape(target)}."
        )
    rows, modules = int(prediction.shape[0]), int(prediction.shape[1])
    if rows <= 0 or modules <= 0:
        raise ValueError("prediction must contain at least one row and one module.")
    if int(n_celltypes) <= 0:
        raise ValueError("n_celltypes must be positive.")
    _require_finite("prediction", prediction)
    _require_finite("target", target)
    pred = _stable_float(prediction).detach()
    fixed_target = target.detach().to(device=pred.device, dtype=pred.dtype)
    ct_id = _integer_id_vector(
        "celltype_id",
        celltype_id,
        batch=rows,
        device=pred.device,
    )
    if bool((ct_id >= int(n_celltypes)).any()):
        raise ValueError(
            f"celltype_id values must be in [0, {int(n_celltypes)})."
        )

    if row_weight is None:
        weight = pred.new_ones(rows)
    else:
        row_weight = _require_tensor("row_weight", row_weight)
        if row_weight.ndim != 1 or int(row_weight.numel()) != rows:
            raise ValueError(
                f"row_weight must have shape [{rows}], got {_shape(row_weight)}."
            )
        _require_finite("row_weight", row_weight)
        weight = row_weight.detach().to(device=pred.device, dtype=pred.dtype)
        if bool((weight < 0).any()):
            raise ValueError("row_weight must be non-negative.")

    if valid_mask is None:
        mask = torch.ones(rows, modules, dtype=torch.bool, device=pred.device)
    else:
        valid_mask = _require_tensor("valid_mask", valid_mask)
        if valid_mask.ndim == 1 and int(valid_mask.numel()) == rows:
            mask = valid_mask[:, None].expand(rows, modules)
        elif _shape(valid_mask) == (rows, modules):
            mask = valid_mask
        else:
            raise ValueError(
                "valid_mask must have shape [N] or [N, M], "
                f"got {_shape(valid_mask)} for N={rows}, M={modules}."
            )
        if mask.dtype != torch.bool:
            raise TypeError(f"valid_mask must be bool, got {mask.dtype}.")
        mask = mask.to(device=pred.device)

    positive = weight > 0
    effective_mask = mask & positive[:, None]
    element_weight = weight[:, None] * effective_mask.to(dtype=pred.dtype)
    flat_celltype = ct_id[:, None].expand(rows, modules).reshape(-1)
    flat_weight = element_weight.reshape(-1)
    flat_pred = pred.reshape(-1)
    flat_target = fixed_target.reshape(-1)
    flat_valid = effective_mask.reshape(-1)

    def reduce_sum(value: torch.Tensor) -> torch.Tensor:
        return pred.new_zeros(int(n_celltypes)).index_add_(
            0,
            flat_celltype,
            value,
        )

    stats = CorrelationCCCSufficientStats(
        weight_sum=reduce_sum(flat_weight),
        observation_count=reduce_sum(flat_valid.to(dtype=pred.dtype)),
        prediction_sum=reduce_sum(flat_weight * flat_pred),
        target_sum=reduce_sum(flat_weight * flat_target),
        prediction_square_sum=reduce_sum(flat_weight * flat_pred.square()),
        target_square_sum=reduce_sum(flat_weight * flat_target.square()),
        cross_sum=reduce_sum(flat_weight * flat_pred * flat_target),
    )
    for value in (
        stats.weight_sum,
        stats.observation_count,
        stats.prediction_sum,
        stats.target_sum,
        stats.prediction_square_sum,
        stats.target_square_sum,
        stats.cross_sum,
    ):
        _require_finite("celltype correlation sufficient statistics", value)
    return stats


def correlation_ccc_from_sufficient_statistics(
    stats: CorrelationCCCSufficientStats,
    *,
    min_observations: int = 2,
    eps: float = 1e-8,
) -> CorrelationCCCOutput:
    """Compute finite correlation/CCC values after optional DDP summation."""

    if not isinstance(stats, CorrelationCCCSufficientStats):
        raise TypeError(
            "stats must be CorrelationCCCSufficientStats, "
            f"got {type(stats).__name__}."
        )
    if int(min_observations) < 2:
        raise ValueError("min_observations must be at least 2.")
    if float(eps) <= 0.0:
        raise ValueError("eps must be positive.")
    tensors = (
        stats.weight_sum,
        stats.observation_count,
        stats.prediction_sum,
        stats.target_sum,
        stats.prediction_square_sum,
        stats.target_square_sum,
        stats.cross_sum,
    )
    reference_shape = _shape(tensors[0])
    if len(reference_shape) != 1:
        raise ValueError(
            f"correlation sufficient statistics must have shape [M], got {reference_shape}."
        )
    for value in tensors:
        if _shape(value) != reference_shape:
            raise ValueError("correlation sufficient statistic shapes do not match.")
        _require_finite("correlation sufficient statistics", value)

    weight = stats.weight_sum
    safe_weight = weight.clamp_min(float(eps))
    pred_mean = stats.prediction_sum / safe_weight
    target_mean = stats.target_sum / safe_weight
    pred_ss = (
        stats.prediction_square_sum
        - stats.prediction_sum.square() / safe_weight
    ).clamp_min(0.0)
    target_ss = (
        stats.target_square_sum
        - stats.target_sum.square() / safe_weight
    ).clamp_min(0.0)
    cross = (
        stats.cross_sum
        - stats.prediction_sum * stats.target_sum / safe_weight
    )
    corr_denominator = torch.sqrt(pred_ss * target_ss)
    ccc_denominator = (
        pred_ss
        + target_ss
        + weight * (pred_mean - target_mean).square()
    )
    valid = (
        (stats.observation_count >= int(min_observations))
        & (weight > float(eps))
        & (corr_denominator > float(eps))
        & (ccc_denominator > float(eps))
    )
    zero = torch.zeros_like(weight)
    correlation = torch.where(
        valid,
        cross / corr_denominator.clamp_min(float(eps)),
        zero,
    ).clamp(min=-1.0, max=1.0)
    ccc = torch.where(
        valid,
        2.0 * cross / ccc_denominator.clamp_min(float(eps)),
        zero,
    ).clamp(min=-1.0, max=1.0)
    _require_finite("correlation", correlation)
    _require_finite("CCC", ccc)
    return CorrelationCCCOutput(
        correlation=correlation,
        ccc=ccc,
        valid=valid,
    )


def differentiable_celltype_module_ccc_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    celltype_id: torch.Tensor,
    *,
    row_weight: Optional[torch.Tensor] = None,
    valid_mask: Optional[torch.Tensor] = None,
    min_observations: int = 3,
    eps: float = 1e-8,
) -> DifferentiableCCCOutput:
    """Return a differentiable donor-axis CCC loss per celltype and module.

    ``prediction`` and ``target`` must already be separately donor-centered
    and train-scale standardized.  The exact held-out endpoint remains
    Spearman; this CCC term is its stable, differentiable training surrogate.
    Constant or insufficiently observed targets are excluded rather than
    contributing unstable zero-variance gradients.
    """

    prediction = _require_tensor("prediction", prediction)
    target = _require_tensor("target", target)
    if prediction.ndim != 2 or _shape(prediction) != _shape(target):
        raise ValueError(
            "prediction and target must have the same shape [N, M], "
            f"got {_shape(prediction)} and {_shape(target)}."
        )
    if int(min_observations) < 2:
        raise ValueError("min_observations must be at least 2.")
    if float(eps) <= 0.0:
        raise ValueError("eps must be positive.")
    rows, modules = int(prediction.shape[0]), int(prediction.shape[1])
    if rows <= 0 or modules <= 0:
        raise ValueError("prediction must contain at least one row and one module.")
    _require_finite("prediction", prediction)
    _require_finite("target", target)
    pred = _stable_float(prediction)
    fixed_target = target.detach().to(device=pred.device, dtype=pred.dtype)
    ct_id = _integer_id_vector(
        "celltype_id",
        celltype_id,
        batch=rows,
        device=pred.device,
    )

    if row_weight is None:
        weight = pred.new_ones(rows)
    else:
        row_weight = _require_tensor("row_weight", row_weight)
        if row_weight.ndim != 1 or int(row_weight.numel()) != rows:
            raise ValueError(
                f"row_weight must have shape [{rows}], got {_shape(row_weight)}."
            )
        _require_finite("row_weight", row_weight)
        weight = row_weight.detach().to(device=pred.device, dtype=pred.dtype)
        if bool((weight < 0).any()):
            raise ValueError("row_weight must be non-negative.")

    if valid_mask is None:
        mask = torch.ones(rows, modules, dtype=torch.bool, device=pred.device)
    else:
        valid_mask = _require_tensor("valid_mask", valid_mask)
        if valid_mask.ndim == 1 and int(valid_mask.numel()) == rows:
            mask = valid_mask[:, None].expand(rows, modules)
        elif _shape(valid_mask) == (rows, modules):
            mask = valid_mask
        else:
            raise ValueError(
                "valid_mask must have shape [N] or [N, M], "
                f"got {_shape(valid_mask)} for N={rows}, M={modules}."
            )
        if mask.dtype != torch.bool:
            raise TypeError(f"valid_mask must be bool, got {mask.dtype}.")
        mask = mask.to(device=pred.device)

    valid_ccc: list[torch.Tensor] = []
    for celltype in torch.unique(ct_id, sorted=True):
        row = ct_id == celltype
        local_pred = pred[row]
        local_target = fixed_target[row]
        local_row_weight = weight[row][:, None]
        local_mask = mask[row] & (local_row_weight > 0)
        local_weight = local_row_weight * local_mask.to(dtype=pred.dtype)
        weight_sum = local_weight.sum(dim=0)
        observation_count = local_mask.sum(dim=0)
        safe_weight = weight_sum.clamp_min(float(eps))
        pred_mean = (local_weight * local_pred).sum(dim=0) / safe_weight
        target_mean = (local_weight * local_target).sum(dim=0) / safe_weight
        pred_delta = local_pred - pred_mean
        target_delta = local_target - target_mean
        covariance = (
            local_weight * pred_delta * target_delta
        ).sum(dim=0) / safe_weight
        pred_variance = (
            local_weight * pred_delta.square()
        ).sum(dim=0) / safe_weight
        target_variance = (
            local_weight * target_delta.square()
        ).sum(dim=0) / safe_weight
        denominator = (
            pred_variance
            + target_variance
            + (pred_mean - target_mean).square()
        )
        # Validity is determined by the frozen observation, never by the
        # prediction.  In particular, a constant prediction against a varying
        # target must receive CCC=0, loss=1, and a gradient that can un-collapse
        # it.  Gating on ``pred_variance`` made collapse a zero-cost solution.
        valid = (
            (observation_count >= int(min_observations))
            & (weight_sum > float(eps))
            & (target_variance > float(eps))
            & (denominator > float(eps))
        )
        if bool(valid.any()):
            ccc = 2.0 * covariance / denominator.clamp_min(float(eps))
            valid_ccc.append(ccc[valid].clamp(min=-1.0, max=1.0))

    if valid_ccc:
        values = torch.cat(valid_ccc)
        mean_ccc = values.mean()
        loss = 1.0 - mean_ccc
        valid_count = pred.new_tensor(float(values.numel()))
    else:
        # Keep the zero graph-connected.  This lets Huber continue training a
        # rare/constant block without manufacturing a CCC gradient.
        loss = pred.sum() * 0.0
        mean_ccc = pred.sum().detach() * 0.0
        valid_count = pred.new_zeros(())
    _require_finite("differentiable CCC loss", loss)
    return DifferentiableCCCOutput(
        loss=loss,
        mean_ccc=mean_ccc.detach(),
        valid_count=valid_count.detach(),
    )


def differentiable_celltype_flattened_ccc_loss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    celltype_id: torch.Tensor,
    *,
    row_weight: Optional[torch.Tensor] = None,
    valid_mask: Optional[torch.Tensor] = None,
    min_observations: int = 3,
    eps: float = 1e-8,
) -> DifferentiableCCCOutput:
    """CCC surrogate aligned to the held-out ``blind_centered`` layout.

    For each cell type, the donor-by-module table is flattened and one CCC is
    computed over all retained entries.  This mirrors the layout of the exact
    endpoint, which computes a Spearman correlation on the same flattened
    table.  CCC remains an amplitude-sensitive differentiable surrogate; the
    exact raw-centered Spearman is still authoritative for checkpoint choice.

    ``valid_mask`` must be fixed from training targets and shared by every
    prediction view.  Prediction variance is deliberately *not* a validity
    condition, so a collapsed prediction receives loss rather than disappearing
    from the objective.
    """

    prediction = _require_tensor("prediction", prediction)
    target = _require_tensor("target", target)
    if prediction.ndim != 2 or _shape(prediction) != _shape(target):
        raise ValueError(
            "prediction and target must have the same shape [N, M], "
            f"got {_shape(prediction)} and {_shape(target)}."
        )
    if int(min_observations) < 2:
        raise ValueError("min_observations must be at least 2.")
    if float(eps) <= 0.0:
        raise ValueError("eps must be positive.")
    rows, modules = int(prediction.shape[0]), int(prediction.shape[1])
    if rows <= 0 or modules <= 0:
        raise ValueError("prediction must contain at least one row and one module.")
    _require_finite("prediction", prediction)
    _require_finite("target", target)
    pred = _stable_float(prediction)
    fixed_target = target.detach().to(device=pred.device, dtype=pred.dtype)
    ct_id = _integer_id_vector(
        "celltype_id", celltype_id, batch=rows, device=pred.device
    )

    if row_weight is None:
        weight = pred.new_ones(rows)
    else:
        row_weight = _require_tensor("row_weight", row_weight)
        if row_weight.ndim != 1 or int(row_weight.numel()) != rows:
            raise ValueError(
                f"row_weight must have shape [{rows}], got {_shape(row_weight)}."
            )
        _require_finite("row_weight", row_weight)
        weight = row_weight.detach().to(device=pred.device, dtype=pred.dtype)
        if bool((weight < 0).any()):
            raise ValueError("row_weight must be non-negative.")

    if valid_mask is None:
        mask = torch.ones(rows, modules, dtype=torch.bool, device=pred.device)
    else:
        valid_mask = _require_tensor("valid_mask", valid_mask)
        if _shape(valid_mask) != (rows, modules):
            raise ValueError(
                "valid_mask must have shape [N, M], "
                f"got {_shape(valid_mask)} for N={rows}, M={modules}."
            )
        if valid_mask.dtype != torch.bool:
            raise TypeError(f"valid_mask must be bool, got {valid_mask.dtype}.")
        mask = valid_mask.to(device=pred.device)

    ccc_values: list[torch.Tensor] = []
    for celltype in torch.unique(ct_id, sorted=True):
        row = ct_id == celltype
        local_mask = mask[row]
        local_weight = weight[row][:, None].expand(-1, modules)
        keep = local_mask & (local_weight > 0)
        observation_count = keep.sum()
        if int(observation_count) < int(min_observations):
            continue
        flat_weight = local_weight[keep]
        flat_pred = pred[row][keep]
        flat_target = fixed_target[row][keep]
        weight_sum = flat_weight.sum()
        if float(weight_sum.detach()) <= float(eps):
            continue
        safe_weight = weight_sum.clamp_min(float(eps))
        pred_mean = (flat_weight * flat_pred).sum() / safe_weight
        target_mean = (flat_weight * flat_target).sum() / safe_weight
        pred_delta = flat_pred - pred_mean
        target_delta = flat_target - target_mean
        covariance = (
            flat_weight * pred_delta * target_delta
        ).sum() / safe_weight
        pred_variance = (
            flat_weight * pred_delta.square()
        ).sum() / safe_weight
        target_variance = (
            flat_weight * target_delta.square()
        ).sum() / safe_weight
        # Target-only validity means branch and full views use the same set.
        if float(target_variance.detach()) <= float(eps):
            continue
        denominator = (
            pred_variance
            + target_variance
            + (pred_mean - target_mean).square()
        )
        ccc = 2.0 * covariance / denominator.clamp_min(float(eps))
        ccc_values.append(ccc.clamp(min=-1.0, max=1.0))

    if ccc_values:
        values = torch.stack(ccc_values)
        mean_ccc = values.mean()
        loss = 1.0 - mean_ccc
        valid_count = pred.new_tensor(float(values.numel()))
    else:
        loss = pred.sum() * 0.0
        mean_ccc = pred.sum().detach() * 0.0
        valid_count = pred.new_zeros(())
    _require_finite("differentiable flattened CCC loss", loss)
    return DifferentiableCCCOutput(
        loss=loss,
        mean_ccc=mean_ccc.detach(),
        valid_count=valid_count.detach(),
    )


def _align_group_table(
    *,
    table_value: torch.Tensor,
    table_group_id: torch.Tensor,
    requested_group_id: torch.Tensor,
    name: str,
    allow_nonfinite: bool = False,
) -> torch.Tensor:
    """Lookup a unique group-keyed table and return requested sorted rows."""

    table_value = _require_tensor(name, table_value)
    if table_value.ndim not in (1, 2):
        raise ValueError(
            f"{name} must have shape [N] or [N, D], got {_shape(table_value)}."
        )
    rows = int(table_value.shape[0])
    if rows <= 0:
        raise ValueError(f"{name} must contain at least one group.")
    ids = _integer_id_vector(
        f"{name}_group_id",
        table_group_id,
        batch=rows,
        device=requested_group_id.device,
    )
    sorted_id, order = torch.sort(ids)
    if int(sorted_id.numel()) > 1 and bool((sorted_id[1:] == sorted_id[:-1]).any()):
        raise ValueError(f"{name}_group_id must contain unique group keys.")
    requested = requested_group_id.to(device=sorted_id.device, dtype=torch.long)
    position = torch.searchsorted(sorted_id, requested)
    safe_position = position.clamp(max=max(0, int(sorted_id.numel()) - 1))
    found = (position < int(sorted_id.numel())) & (
        sorted_id.index_select(0, safe_position) == requested
    )
    if not bool(found.all()):
        missing = requested[~found].detach().cpu().tolist()
        raise ValueError(
            f"{name} is missing {len(missing)} requested group ids; "
            f"examples={missing[:8]}."
        )
    sorted_value = table_value.detach().to(
        device=requested_group_id.device,
    ).index_select(0, order)
    result = sorted_value.index_select(0, position)
    if not bool(allow_nonfinite):
        _require_finite(name, result)
    return result


def donor_celltype_centered_huber(
    prediction_activity: torch.Tensor,
    target_activity: Optional[torch.Tensor] = None,
    *,
    group_id: torch.Tensor,
    celltype_id: torch.Tensor,
    train_scale: torch.Tensor,
    target_group_activity: Optional[torch.Tensor] = None,
    target_group_id: Optional[torch.Tensor] = None,
    target_group_cell_count: Optional[torch.Tensor] = None,
    donor_weight: Optional[torch.Tensor] = None,
    context_weight: Optional[torch.Tensor] = None,
    cell_valid_mask: Optional[torch.Tensor] = None,
    module_valid_mask: Optional[torch.Tensor] = None,
    ccc_valid_mask: Optional[torch.Tensor] = None,
    min_cells_per_group: int = 1,
    min_donor_groups_per_celltype: int = 2,
    min_metric_observations: int = 2,
    delta: float = 1.0,
    eps: float = 1e-8,
) -> CenteredModuleViewOutput:
    """Compare separately donor-centered prediction/observation pseudobulks.

    The order of operations is a scientific contract:

    1. aggregate predicted cells to donor×celltype means;
    2. exclude groups below ``min_cells_per_group``;
    3. within each cell type, center prediction and observation *separately*
       across eligible donor groups;
    4. divide both deviations by frozen train-only scale;
    5. apply donor/context-weighted Huber and emit correlation/CCC statistics.

    The preferred observed target is ``target_group_activity`` from the frozen
    full donor×celltype artifact, keyed by ``target_group_id``.  The alternative
    ``target_activity`` path aggregates sampled target cells and is retained
    for mathematical tests and controlled fallback only.  Exactly one target
    source is required.

    The prediction input must cover complete cell-type donor blocks.  This
    function must not be called independently on arbitrary minibatches or DDP
    rank shards.
    """

    if int(min_cells_per_group) <= 0:
        raise ValueError("min_cells_per_group must be positive.")
    if int(min_donor_groups_per_celltype) < 2:
        raise ValueError("min_donor_groups_per_celltype must be at least 2.")
    if float(eps) <= 0.0:
        raise ValueError("eps must be positive.")
    prediction_activity = _require_tensor(
        "prediction_activity",
        prediction_activity,
    )
    if prediction_activity.ndim != 2:
        raise ValueError("prediction_activity must have shape [B, M].")
    batch, modules = (
        int(prediction_activity.shape[0]),
        int(prediction_activity.shape[1]),
    )
    if batch <= 0 or modules <= 0:
        raise ValueError("activity must contain at least one cell and one module.")

    prediction_group = aggregate_donor_context_mean(
        prediction_activity,
        group_id,
        celltype_id,
        cell_valid_mask=cell_valid_mask,
    )
    has_cell_target = target_activity is not None
    has_group_target = target_group_activity is not None or target_group_id is not None
    if has_cell_target == has_group_target:
        raise ValueError(
            "provide exactly one observed target source: target_activity or "
            "(target_group_activity, target_group_id)."
        )
    if has_cell_target:
        target_activity = _require_tensor("target_activity", target_activity)
        if _shape(prediction_activity) != _shape(target_activity):
            raise ValueError(
                "prediction_activity and target_activity shapes differ: "
                f"{_shape(prediction_activity)} vs {_shape(target_activity)}."
            )
        target_group = aggregate_donor_context_mean(
            target_activity.detach(),
            group_id,
            celltype_id,
            cell_valid_mask=cell_valid_mask,
        )
        target_support_count = target_group.cell_count
    else:
        if target_group_activity is None or target_group_id is None:
            raise ValueError(
                "target_group_activity and target_group_id must be supplied together."
            )
        target_group_activity = _require_tensor(
            "target_group_activity",
            target_group_activity,
        )
        if target_group_activity.ndim != 2 or int(target_group_activity.shape[1]) != modules:
            raise ValueError(
                f"target_group_activity must have shape [N, {modules}], "
                f"got {_shape(target_group_activity)}."
            )
        aligned_target = _align_group_table(
            table_value=target_group_activity,
            table_group_id=target_group_id,
            requested_group_id=prediction_group.group_id,
            name="target_group_activity",
            allow_nonfinite=target_group_cell_count is not None,
        ).to(dtype=prediction_group.mean.dtype)
        if target_group_cell_count is None:
            target_support_count = prediction_group.cell_count
        else:
            aligned_count = _align_group_table(
                table_value=target_group_cell_count,
                table_group_id=target_group_id,
                requested_group_id=prediction_group.group_id,
                name="target_group_cell_count",
            )
            if aligned_count.ndim != 1:
                raise ValueError("target_group_cell_count must have shape [N].")
            if bool((aligned_count < 0).any()):
                raise ValueError("target_group_cell_count must be non-negative.")
            if aligned_count.is_floating_point() and not bool(
                torch.equal(aligned_count, aligned_count.round())
            ):
                raise ValueError("target_group_cell_count must be integer-valued.")
            target_support_count = aligned_count.to(dtype=torch.long)
        target_group = GroupedModuleActivity(
            mean=aligned_target.detach(),
            group_id=prediction_group.group_id,
            celltype_id=prediction_group.celltype_id,
            cell_count=target_support_count,
        )

    train_scale = _require_tensor("train_scale", train_scale)
    if train_scale.ndim != 2 or int(train_scale.shape[1]) != modules:
        raise ValueError(
            f"train_scale must have shape [T, {modules}], got {_shape(train_scale)}."
        )
    if int(train_scale.shape[0]) <= 0:
        raise ValueError("train_scale must contain at least one cell type.")
    _require_finite("train_scale", train_scale)
    if not bool((train_scale > 0).all()):
        raise ValueError("train_scale must be strictly positive for every celltype/module.")

    group_celltype = prediction_group.celltype_id
    n_celltypes = int(train_scale.shape[0])
    if bool((group_celltype >= n_celltypes).any()):
        raise ValueError(
            f"group celltype ids exceed train_scale cell types T={n_celltypes}."
        )

    # Group-level weights are supplied repeated on cells and averaged only to
    # recover their group value.  They are detached data definitions.
    group_weight = prediction_group.mean.new_ones(
        int(prediction_group.group_id.numel())
    )
    for name, supplied in (
        ("donor_weight", donor_weight),
        ("context_weight", context_weight),
    ):
        if supplied is None:
            continue
        supplied = _require_tensor(name, supplied)
        if supplied.ndim != 1 or int(supplied.numel()) != batch:
            raise ValueError(f"{name} must have shape [{batch}], got {_shape(supplied)}.")
        _require_finite(name, supplied)
        fixed = supplied.detach()
        if bool((fixed < 0).any()):
            raise ValueError(f"{name} must be non-negative.")
        grouped = aggregate_donor_context_mean(
            fixed[:, None],
            group_id,
            celltype_id,
            cell_valid_mask=cell_valid_mask,
        )
        group_weight = group_weight * grouped.mean[:, 0].to(
            device=group_weight.device,
            dtype=group_weight.dtype,
        )

    enough_cells = (
        (prediction_group.cell_count >= int(min_cells_per_group))
        & (target_support_count >= int(min_cells_per_group))
    )
    if bool(enough_cells.any()):
        _require_finite(
            "eligible target_group_activity",
            target_group.mean[enough_cells],
        )
    target_group = GroupedModuleActivity(
        mean=torch.where(
            enough_cells[:, None],
            target_group.mean,
            torch.zeros_like(target_group.mean),
        ),
        group_id=target_group.group_id,
        celltype_id=target_group.celltype_id,
        cell_count=target_group.cell_count,
    )
    positive_weight = group_weight > 0
    center_eligible = enough_cells & positive_weight
    celltype_group_count = torch.zeros(
        n_celltypes,
        dtype=torch.long,
        device=group_celltype.device,
    ).index_add_(0, group_celltype, center_eligible.to(dtype=torch.long))
    eligible_celltype = (
        celltype_group_count >= int(min_donor_groups_per_celltype)
    )
    valid_group = center_eligible & eligible_celltype.index_select(
        0,
        group_celltype,
    )
    if not bool(valid_group.any()):
        raise ValueError(
            "no donor×celltype group passes min_cells_per_group and "
            "min_donor_groups_per_celltype."
        )

    n_groups = int(prediction_group.group_id.numel())
    if module_valid_mask is None:
        module_valid = torch.ones(
            n_groups,
            modules,
            dtype=torch.bool,
            device=prediction_group.mean.device,
        )
    else:
        module_valid_mask = _require_tensor(
            "module_valid_mask",
            module_valid_mask,
        )
        if module_valid_mask.ndim == 1 and int(module_valid_mask.numel()) == modules:
            module_valid = module_valid_mask[None, :].expand(n_groups, modules)
        elif _shape(module_valid_mask) == (n_groups, modules):
            module_valid = module_valid_mask
        else:
            raise ValueError(
                "module_valid_mask must have shape [M] or [N_groups, M], "
                f"got {_shape(module_valid_mask)}."
            )
        if module_valid.dtype != torch.bool:
            raise TypeError(
                f"module_valid_mask must be bool, got {module_valid.dtype}."
            )
        module_valid = module_valid.to(device=prediction_group.mean.device)

    # Centers are uniform over eligible donor groups within each cell type.
    # This is intentionally separate for prediction and target; subtracting
    # one shared train mean from both would algebraically cancel from the error.
    # A group×module validity mask participates in the center itself; applying
    # it only after centering would let excluded donor-module entries contaminate
    # every retained residual in that cell type.
    center_mask = center_eligible[:, None] & module_valid
    center_weight = center_mask.to(dtype=prediction_group.mean.dtype)
    center_denominator = prediction_group.mean.new_zeros(
        n_celltypes,
        modules,
    ).index_add_(
        0,
        group_celltype,
        center_weight,
    )
    prediction_center = prediction_group.mean.new_zeros(
        n_celltypes,
        modules,
    ).index_add_(
        0,
        group_celltype,
        prediction_group.mean * center_weight,
    )
    prediction_center = prediction_center / center_denominator.clamp_min(1.0)
    target_center = target_group.mean.new_zeros(
        n_celltypes,
        modules,
    ).index_add_(
        0,
        group_celltype,
        target_group.mean * center_weight,
    )
    target_center = target_center / center_denominator.clamp_min(1.0)

    scale = train_scale.detach().to(
        device=prediction_group.mean.device,
        dtype=prediction_group.mean.dtype,
    ).index_select(0, group_celltype)
    prediction_standardized = (
        prediction_group.mean
        - prediction_center.index_select(0, group_celltype)
    ) / scale
    target_standardized = (
        target_group.mean
        - target_center.index_select(0, group_celltype)
    ) / scale
    target_standardized = target_standardized.detach()
    _require_finite(
        "centered standardized prediction activity",
        prediction_standardized,
    )
    _require_finite(
        "centered standardized target activity",
        target_standardized,
    )

    eligible_celltype_module = center_denominator >= float(
        min_donor_groups_per_celltype
    )
    loss_mask = (
        module_valid
        & valid_group[:, None]
        & eligible_celltype_module.index_select(0, group_celltype)
    )

    huber = donor_context_weighted_huber(
        prediction_standardized,
        target_standardized,
        donor_weight=group_weight,
        valid_mask=loss_mask,
        delta=float(delta),
        eps=float(eps),
    )
    stats = correlation_ccc_sufficient_statistics(
        prediction_standardized,
        target_standardized,
        row_weight=group_weight,
        valid_mask=loss_mask,
    )
    metrics = correlation_ccc_from_sufficient_statistics(
        stats,
        min_observations=int(min_metric_observations),
        eps=float(eps),
    )
    celltype_stats = celltype_correlation_ccc_sufficient_statistics(
        prediction_standardized,
        target_standardized,
        group_celltype,
        n_celltypes=n_celltypes,
        row_weight=group_weight,
        valid_mask=loss_mask,
    )
    celltype_metrics = correlation_ccc_from_sufficient_statistics(
        celltype_stats,
        min_observations=int(min_metric_observations),
        eps=float(eps),
    )
    if ccc_valid_mask is None:
        ccc_module_valid = torch.ones(
            n_celltypes,
            modules,
            dtype=torch.bool,
            device=prediction_standardized.device,
        )
    else:
        ccc_valid_mask = _require_tensor("ccc_valid_mask", ccc_valid_mask)
        if _shape(ccc_valid_mask) != (n_celltypes, modules):
            raise ValueError(
                "ccc_valid_mask must have shape [T, M], "
                f"got {_shape(ccc_valid_mask)} for T={n_celltypes}, "
                f"M={modules}."
            )
        if ccc_valid_mask.dtype != torch.bool:
            raise TypeError(
                f"ccc_valid_mask must be bool, got {ccc_valid_mask.dtype}."
            )
        ccc_module_valid = ccc_valid_mask.detach().to(
            device=prediction_standardized.device
        )
    ccc_loss_mask = loss_mask & ccc_module_valid.index_select(
        0, group_celltype
    )
    prediction_centered = prediction_group.mean - prediction_center.index_select(
        0, group_celltype
    )
    target_centered = (
        target_group.mean - target_center.index_select(0, group_celltype)
    ).detach()
    concordance = differentiable_celltype_flattened_ccc_loss(
        prediction_centered,
        target_centered,
        group_celltype,
        row_weight=group_weight,
        valid_mask=ccc_loss_mask,
        min_observations=max(3, int(min_metric_observations)),
        eps=float(eps),
    )
    return CenteredModuleViewOutput(
        huber=huber,
        prediction=prediction_standardized,
        target=target_standardized,
        group_id=prediction_group.group_id,
        group_celltype_id=group_celltype,
        group_cell_count=prediction_group.cell_count,
        group_weight=group_weight,
        valid_group=valid_group,
        celltype_group_count=celltype_group_count,
        correlation_stats=stats,
        metrics=metrics,
        celltype_correlation_stats=celltype_stats,
        celltype_metrics=celltype_metrics,
        concordance=concordance,
    )


def dual_view_module_rescue_loss(
    *,
    config: PrismModuleRescueConfig,
    branch_probabilities: Optional[torch.Tensor] = None,
    full_probabilities: Optional[torch.Tensor] = None,
    observed_tier: Optional[torch.Tensor] = None,
    target_group_activity: Optional[torch.Tensor] = None,
    target_group_id: Optional[torch.Tensor] = None,
    target_group_cell_count: Optional[torch.Tensor] = None,
    activity_weight: Optional[torch.Tensor] = None,
    group_id: Optional[torch.Tensor] = None,
    celltype_id: Optional[torch.Tensor] = None,
    train_scale: Optional[torch.Tensor] = None,
    donor_weight: Optional[torch.Tensor] = None,
    context_weight: Optional[torch.Tensor] = None,
    cell_valid_mask: Optional[torch.Tensor] = None,
    module_valid_mask: Optional[torch.Tensor] = None,
    ccc_valid_mask: Optional[torch.Tensor] = None,
    min_cells_per_group: int = 1,
    min_donor_groups_per_celltype: int = 2,
    min_metric_observations: int = 2,
) -> PrismModuleRescueOutput:
    """Compute dual-view rescue on donor-centered donor×celltype activities.

    ``branch_probabilities`` are the target-latent-free cross-context PRISM
    view. ``full_probabilities`` are the ordinary unmasked reconstruction view.
    Both are evaluated in the same fixed activity space, but they retain their
    different scientific meanings.

    The preferred target is the full donor×celltype artifact supplied through
    ``target_group_activity`` and ``target_group_id``.  ``observed_tier`` is a
    sampled-cell fallback; exactly one target source is allowed.

    When disabled, all rescue inputs may be ``None``. ``loss`` includes
    ``lambda_module``; ``rescue_loss`` is the unweighted branch/full mixture.
    """

    if not isinstance(config, PrismModuleRescueConfig):
        raise TypeError(
            "config must be PrismModuleRescueConfig, "
            f"got {type(config).__name__}."
        )
    reference = (
        branch_probabilities
        if isinstance(branch_probabilities, torch.Tensor)
        else full_probabilities
        if isinstance(full_probabilities, torch.Tensor)
        else observed_tier
        if isinstance(observed_tier, torch.Tensor)
        else target_group_activity
        if isinstance(target_group_activity, torch.Tensor)
        else None
    )
    if not bool(config.enabled):
        zero = _zero(reference)
        return PrismModuleRescueOutput(
            loss=zero,
            rescue_loss=zero,
            branch_loss=zero,
            full_loss=zero,
            branch_numerator=zero,
            branch_denominator=zero,
            full_numerator=zero,
            full_denominator=zero,
            huber_rescue_loss=zero,
            ccc_rescue_loss=zero,
            branch_huber_loss=zero,
            full_huber_loss=zero,
            branch_ccc_loss=zero,
            full_ccc_loss=zero,
            branch_mean_ccc=zero,
            full_mean_ccc=zero,
        )
    if float(config.ccc_weight) > 0.0 and ccc_valid_mask is None:
        raise ValueError(
            "CCC-enabled module rescue requires a frozen train-only "
            "ccc_valid_mask with shape [T, M]."
        )

    required = {
        "branch_probabilities": branch_probabilities,
        "full_probabilities": full_probabilities,
        "activity_weight": activity_weight,
        "group_id": group_id,
        "celltype_id": celltype_id,
        "train_scale": train_scale,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError(
            "enabled module rescue requires: " + ", ".join(sorted(missing)) + "."
        )
    has_cell_target = observed_tier is not None
    has_group_target = target_group_activity is not None or target_group_id is not None
    if has_cell_target == has_group_target:
        raise ValueError(
            "enabled module rescue requires exactly one observed target source: "
            "observed_tier or (target_group_activity, target_group_id)."
        )
    if has_group_target and (
        target_group_activity is None or target_group_id is None
    ):
        raise ValueError(
            "target_group_activity and target_group_id must be supplied together."
        )
    if float(config.lambda_module) < 0.0:
        raise ValueError("lambda_module must be non-negative.")
    rho = float(config.branch_fraction)
    if not (0.0 <= rho <= 1.0):
        raise ValueError("branch_fraction must be in [0, 1].")
    if float(config.huber_delta) <= 0.0:
        raise ValueError("huber_delta must be positive.")
    if not torch.isfinite(torch.tensor(float(config.ccc_weight))):
        raise ValueError("ccc_weight must be finite.")
    if float(config.ccc_weight) < 0.0:
        raise ValueError("ccc_weight must be non-negative.")
    if float(config.eps) <= 0.0:
        raise ValueError("eps must be positive.")

    # Narrow the Optional types after the explicit contract check above.
    branch_probs = _require_tensor("branch_probabilities", branch_probabilities)
    full_probs = _require_tensor("full_probabilities", full_probabilities)
    fixed_weight = _require_tensor("activity_weight", activity_weight)
    fixed_group_id = _require_tensor("group_id", group_id)
    ct_id = _require_tensor("celltype_id", celltype_id)
    fixed_scale = _require_tensor("train_scale", train_scale)

    if _shape(branch_probs) != _shape(full_probs):
        raise ValueError(
            "branch_probabilities and full_probabilities shapes differ: "
            f"{_shape(branch_probs)} vs {_shape(full_probs)}."
        )
    if branch_probs.ndim != 3:
        raise ValueError(
            "branch_probabilities and full_probabilities must have shape [B, G, C]."
        )
    batch, genes, bins = (
        int(branch_probs.shape[0]),
        int(branch_probs.shape[1]),
        int(branch_probs.shape[2]),
    )
    branch_group_activity = grouped_expected_module_activity(
        branch_probs,
        fixed_weight,
        fixed_group_id,
        ct_id,
        cell_valid_mask=cell_valid_mask,
        validate_probabilities=bool(config.validate_probabilities),
        probability_atol=float(config.probability_atol),
    )
    full_group_activity = grouped_expected_module_activity(
        full_probs,
        fixed_weight,
        fixed_group_id,
        ct_id,
        cell_valid_mask=cell_valid_mask,
        validate_probabilities=bool(config.validate_probabilities),
        probability_atol=float(config.probability_atol),
    )
    if not bool(
        torch.equal(
            branch_group_activity.group_id,
            full_group_activity.group_id,
        )
    ):
        raise RuntimeError("branch and full group order unexpectedly differs.")
    branch_activity = _expand_group_rows_to_cells(
        branch_group_activity,
        fixed_group_id,
    )
    full_activity = _expand_group_rows_to_cells(
        full_group_activity,
        fixed_group_id,
    )

    target_activity: Optional[torch.Tensor]
    if observed_tier is not None:
        target_tier = _require_tensor("observed_tier", observed_tier)
        if target_tier.ndim != 2 or _shape(target_tier) != (batch, genes):
            raise ValueError(
                f"observed_tier must have shape [{batch}, {genes}], "
                f"got {_shape(target_tier)}."
            )
        _require_finite("observed_tier", target_tier)
        target_float = _stable_float(target_tier)
        if bool((target_float < 0).any()) or bool(
            (target_float > bins - 1).any()
        ):
            raise ValueError(f"observed_tier values must lie in [0, {bins - 1}].")
        if target_tier.is_floating_point() and not bool(
            torch.equal(target_float, target_float.round())
        ):
            raise ValueError("observed_tier must contain integer-valued ordinal bins.")
        target_group_tier = aggregate_donor_context_mean(
            target_float,
            fixed_group_id,
            ct_id,
            cell_valid_mask=cell_valid_mask,
        )
        target_group_projected = GroupedModuleActivity(
            mean=project_fixed_module_activity(
                target_group_tier.mean,
                fixed_weight,
            ),
            group_id=target_group_tier.group_id,
            celltype_id=target_group_tier.celltype_id,
            cell_count=target_group_tier.cell_count,
        )
        target_activity = _expand_group_rows_to_cells(
            target_group_projected,
            fixed_group_id,
        )
    else:
        target_activity = None

    view_arguments = dict(
        target_activity=target_activity,
        group_id=fixed_group_id,
        celltype_id=ct_id,
        train_scale=fixed_scale,
        target_group_activity=target_group_activity,
        target_group_id=target_group_id,
        target_group_cell_count=target_group_cell_count,
        donor_weight=donor_weight,
        context_weight=context_weight,
        cell_valid_mask=cell_valid_mask,
        module_valid_mask=module_valid_mask,
        ccc_valid_mask=ccc_valid_mask,
        min_cells_per_group=int(min_cells_per_group),
        min_donor_groups_per_celltype=int(min_donor_groups_per_celltype),
        min_metric_observations=int(min_metric_observations),
        delta=float(config.huber_delta),
        eps=float(config.eps),
    )
    branch = donor_celltype_centered_huber(
        branch_activity,
        **view_arguments,
    )
    full = donor_celltype_centered_huber(
        full_activity,
        **view_arguments,
    )
    if not bool(torch.equal(branch.group_id, full.group_id)):
        raise RuntimeError("branch and full grouped outputs unexpectedly differ.")
    huber_rescue = rho * branch.huber.loss + (1.0 - rho) * full.huber.loss
    ccc_rescue = (
        rho * branch.concordance.loss
        + (1.0 - rho) * full.concordance.loss
    )
    branch_loss = (
        branch.huber.loss
        + float(config.ccc_weight) * branch.concordance.loss
    )
    full_loss = (
        full.huber.loss
        + float(config.ccc_weight) * full.concordance.loss
    )
    rescue = huber_rescue + float(config.ccc_weight) * ccc_rescue
    loss = float(config.lambda_module) * rescue
    _require_finite("dual-view module rescue loss", loss)
    return PrismModuleRescueOutput(
        loss=loss,
        rescue_loss=rescue,
        branch_loss=branch_loss,
        full_loss=full_loss,
        branch_numerator=branch.huber.numerator,
        branch_denominator=branch.huber.denominator,
        full_numerator=full.huber.numerator,
        full_denominator=full.huber.denominator,
        huber_rescue_loss=huber_rescue,
        ccc_rescue_loss=ccc_rescue,
        branch_huber_loss=branch.huber.loss,
        full_huber_loss=full.huber.loss,
        branch_ccc_loss=branch.concordance.loss,
        full_ccc_loss=full.concordance.loss,
        branch_mean_ccc=branch.concordance.mean_ccc,
        full_mean_ccc=full.concordance.mean_ccc,
        branch_activity=branch.prediction,
        full_activity=full.prediction,
        target_activity=branch.target,
        branch_view=branch,
        full_view=full,
    )


def dual_view_grouped_module_rescue_loss(
    *,
    config: PrismModuleRescueConfig,
    branch_group_activity: torch.Tensor,
    full_group_activity: torch.Tensor,
    prediction_group_id: torch.Tensor,
    prediction_group_celltype_id: torch.Tensor,
    prediction_group_cell_count: torch.Tensor,
    target_group_activity: torch.Tensor,
    target_group_id: torch.Tensor,
    target_group_cell_count: Optional[torch.Tensor],
    train_scale: torch.Tensor,
    donor_group_weight: Optional[torch.Tensor] = None,
    context_group_weight: Optional[torch.Tensor] = None,
    module_valid_mask: Optional[torch.Tensor] = None,
    ccc_valid_mask: Optional[torch.Tensor] = None,
    min_cells_per_group: int = 1,
    min_donor_groups_per_celltype: int = 2,
    min_metric_observations: int = 2,
) -> PrismModuleRescueOutput:
    """Evaluate rescue from pre-aggregated prediction group means.

    This is the low-memory endpoint used by deterministic two-pass replay.
    ``branch_group_activity`` and ``full_group_activity`` contain one row per
    unique donor×celltype group.  The nonlinear donor-centering, Huber and CCC
    terms are therefore built only on the small group table; a caller may then
    replay cell microbatches with the returned mean VJP.

    The mathematical result is identical to :func:`dual_view_module_rescue_loss`
    after its probability-to-group-activity projection.  Actual prediction
    group counts are used for eligibility and restored on the diagnostic view;
    the internal one-row-per-group representation is never interpreted as one
    sampled nucleus.
    """

    if not isinstance(config, PrismModuleRescueConfig):
        raise TypeError(
            "config must be PrismModuleRescueConfig, "
            f"got {type(config).__name__}."
        )
    branch = _require_tensor("branch_group_activity", branch_group_activity)
    if not bool(config.enabled):
        zero = _zero(branch)
        return PrismModuleRescueOutput(
            loss=zero,
            rescue_loss=zero,
            branch_loss=zero,
            full_loss=zero,
            branch_numerator=zero,
            branch_denominator=zero,
            full_numerator=zero,
            full_denominator=zero,
            huber_rescue_loss=zero,
            ccc_rescue_loss=zero,
            branch_huber_loss=zero,
            full_huber_loss=zero,
            branch_ccc_loss=zero,
            full_ccc_loss=zero,
            branch_mean_ccc=zero,
            full_mean_ccc=zero,
        )
    full = _require_tensor("full_group_activity", full_group_activity)
    group = _require_tensor("prediction_group_id", prediction_group_id)
    celltype = _require_tensor(
        "prediction_group_celltype_id",
        prediction_group_celltype_id,
    )
    count = _require_tensor(
        "prediction_group_cell_count",
        prediction_group_cell_count,
    )
    if branch.ndim != 2 or _shape(full) != _shape(branch):
        raise ValueError(
            "branch_group_activity and full_group_activity must share shape "
            f"[N, M], got {_shape(branch)} and {_shape(full)}."
        )
    groups, modules = int(branch.shape[0]), int(branch.shape[1])
    if groups <= 0 or modules <= 0:
        raise ValueError("group activity must contain at least one row and module.")
    if _shape(group) != (groups,) or _shape(celltype) != (groups,):
        raise ValueError(
            "prediction group and celltype ids must have shape "
            f"[{groups}]."
        )
    if _shape(count) != (groups,):
        raise ValueError(
            f"prediction_group_cell_count must have shape [{groups}]."
        )
    if count.is_floating_point() and not bool(torch.equal(count, count.round())):
        raise ValueError("prediction_group_cell_count must be integer-valued.")
    count = count.to(device=branch.device, dtype=torch.long)
    if bool((count < 0).any()):
        raise ValueError("prediction_group_cell_count must be non-negative.")
    group = _integer_id_vector(
        "prediction_group_id",
        group,
        batch=groups,
        device=branch.device,
    )
    celltype = _integer_id_vector(
        "prediction_group_celltype_id",
        celltype,
        batch=groups,
        device=branch.device,
    )
    if int(torch.unique(group).numel()) != groups:
        raise ValueError("prediction_group_id must contain one unique row per group.")
    if int(min_cells_per_group) <= 0:
        raise ValueError("min_cells_per_group must be positive.")
    if float(config.lambda_module) < 0.0:
        raise ValueError("lambda_module must be non-negative.")
    rho = float(config.branch_fraction)
    if not 0.0 <= rho <= 1.0:
        raise ValueError("branch_fraction must be in [0, 1].")
    if float(config.huber_delta) <= 0.0:
        raise ValueError("huber_delta must be positive.")
    if not math.isfinite(float(config.ccc_weight)) or float(config.ccc_weight) < 0.0:
        raise ValueError("ccc_weight must be finite and non-negative.")
    if float(config.ccc_weight) > 0.0 and ccc_valid_mask is None:
        raise ValueError(
            "CCC-enabled module rescue requires a frozen train-only "
            "ccc_valid_mask with shape [T, M]."
        )
    _require_finite("branch_group_activity", branch)
    _require_finite("full_group_activity", full)

    group_valid = count >= int(min_cells_per_group)
    view_arguments = dict(
        target_group_activity=target_group_activity,
        target_group_id=target_group_id,
        target_group_cell_count=target_group_cell_count,
        group_id=group,
        celltype_id=celltype,
        train_scale=train_scale,
        donor_weight=donor_group_weight,
        context_weight=context_group_weight,
        cell_valid_mask=group_valid,
        module_valid_mask=module_valid_mask,
        ccc_valid_mask=ccc_valid_mask,
        # Each valid prediction row is already the mean of ``count`` cells.
        min_cells_per_group=1,
        min_donor_groups_per_celltype=int(min_donor_groups_per_celltype),
        min_metric_observations=int(min_metric_observations),
        delta=float(config.huber_delta),
        eps=float(config.eps),
    )
    branch_view = donor_celltype_centered_huber(branch, **view_arguments)
    full_view = donor_celltype_centered_huber(full, **view_arguments)
    if not bool(torch.equal(branch_view.group_id, group)):
        raise RuntimeError("grouped branch rows were unexpectedly reordered.")
    if not bool(torch.equal(full_view.group_id, group)):
        raise RuntimeError("grouped full rows were unexpectedly reordered.")
    branch_view = replace(branch_view, group_cell_count=count)
    full_view = replace(full_view, group_cell_count=count)

    huber_rescue = (
        rho * branch_view.huber.loss
        + (1.0 - rho) * full_view.huber.loss
    )
    ccc_rescue = (
        rho * branch_view.concordance.loss
        + (1.0 - rho) * full_view.concordance.loss
    )
    branch_loss = (
        branch_view.huber.loss
        + float(config.ccc_weight) * branch_view.concordance.loss
    )
    full_loss = (
        full_view.huber.loss
        + float(config.ccc_weight) * full_view.concordance.loss
    )
    rescue = huber_rescue + float(config.ccc_weight) * ccc_rescue
    loss = float(config.lambda_module) * rescue
    _require_finite("dual-view grouped module rescue loss", loss)
    return PrismModuleRescueOutput(
        loss=loss,
        rescue_loss=rescue,
        branch_loss=branch_loss,
        full_loss=full_loss,
        branch_numerator=branch_view.huber.numerator,
        branch_denominator=branch_view.huber.denominator,
        full_numerator=full_view.huber.numerator,
        full_denominator=full_view.huber.denominator,
        huber_rescue_loss=huber_rescue,
        ccc_rescue_loss=ccc_rescue,
        branch_huber_loss=branch_view.huber.loss,
        full_huber_loss=full_view.huber.loss,
        branch_ccc_loss=branch_view.concordance.loss,
        full_ccc_loss=full_view.concordance.loss,
        branch_mean_ccc=branch_view.concordance.mean_ccc,
        full_mean_ccc=full_view.concordance.mean_ccc,
        branch_activity=branch_view.prediction,
        full_activity=full_view.prediction,
        target_activity=branch_view.target,
        branch_view=branch_view,
        full_view=full_view,
    )


__all__ = [
    "CenteredModuleViewOutput",
    "CorrelationCCCOutput",
    "CorrelationCCCSufficientStats",
    "DifferentiableCCCOutput",
    "GroupedModuleActivity",
    "PrismModuleRescueConfig",
    "PrismModuleRescueOutput",
    "WeightedHuberOutput",
    "aggregate_donor_context_mean",
    "celltype_correlation_ccc_sufficient_statistics",
    "correlation_ccc_from_sufficient_statistics",
    "correlation_ccc_sufficient_statistics",
    "donor_celltype_centered_huber",
    "donor_context_weighted_huber",
    "differentiable_celltype_module_ccc_loss",
    "differentiable_celltype_flattened_ccc_loss",
    "dual_view_module_rescue_loss",
    "dual_view_grouped_module_rescue_loss",
    "expected_ordinal_tier",
    "grouped_expected_module_activity",
    "project_fixed_module_activity",
    "standardize_module_activity",
]
