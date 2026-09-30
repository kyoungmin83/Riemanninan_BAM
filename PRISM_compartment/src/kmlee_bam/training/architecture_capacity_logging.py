"""Optimizer-neutral live logging for learned PRISM architecture capacity."""

from __future__ import annotations

import math
from typing import Any, Mapping

import torch


def _finite_metric(metrics: Mapping[str, float], key: str) -> float | None:
    try:
        value = float(metrics[key])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _generator_mode(config: Any, epoch: int) -> str:
    start = int(config.start_epoch)
    shadow_end = int(getattr(config, "shadow_end_epoch", 0) or 0)
    soft_end = int(getattr(config, "soft_end_epoch", 0) or 0)
    if int(epoch) < start:
        return "warmup_all_on"
    if shadow_end >= start and int(epoch) <= shadow_end:
        return "shadow_all_on"
    if soft_end > 0 and int(epoch) <= soft_end:
        return "soft_adaptation"
    return "protected_hard"


@torch.no_grad()
def build_architecture_capacity_record(
    *,
    epoch: int,
    total_epochs: int,
    pathology_scale: float,
    rank_gate: Any | None,
    generator_gate: Any | None,
    generator_config: Any | None,
    metrics: Mapping[str, float],
    paper_config: Any | None = None,
) -> dict[str, Any]:
    """Build a JSON-safe snapshot without changing parameters or RNG state."""

    record: dict[str, Any] = {
        "schema_version": "kmlee_bam.architecture_capacity_live.v2",
        "epoch": int(epoch),
        "total_epochs": int(total_epochs),
        "pathology_route": {
            "scale": float(pathology_scale),
            "decoder_absmean": _finite_metric(metrics, "metric/v26_pathdec_absmean"),
            "auxiliary_loss": _finite_metric(metrics, "loss/v20_path_aux"),
        },
    }

    paper_enabled = bool(
        paper_config is not None
        and getattr(paper_config, "module_local_nonlinear_enabled", False)
    )
    if paper_enabled:
        record["paper_compartmental_branch"] = {
            "enabled": True,
            "source_pdf": (
                "ref/Dendritic morphology and synaptic nonlinearities "
                "enhancefunctional complexity in human cortical neurons.pdf"
            ),
            "interpretation": (
                "computational analogy; not an observed dendrite or NMDA measurement"
            ),
            "variant": str(
                getattr(
                    paper_config,
                    "module_local_nonlinear_variant",
                    "unknown",
                )
            ),
            "ramp": _finite_metric(
                metrics, "weight/prism_module_local_nonlinear_ramp"
            ),
            "local_rms": _finite_metric(
                metrics, "metric/prism_module_local_module_rms"
            ),
            "nonlinear_rms": _finite_metric(
                metrics, "metric/prism_module_local_nonlinear_rms"
            ),
            "nonlinear_to_local_ratio": _finite_metric(
                metrics, "metric/prism_module_local_nonlinear_to_local_ratio"
            ),
            "threshold_crossing_fraction": _finite_metric(
                metrics,
                "metric/prism_module_local_threshold_crossing_fraction",
            ),
            "branch_nll_gain_when_enabled": _finite_metric(
                metrics,
                "metric/prism_module_local_nonlinear_branch_nll_gain",
            ),
            "full_nll_gain_when_enabled": _finite_metric(
                metrics,
                "metric/prism_module_local_nonlinear_full_nll_gain",
            ),
            "mix_mean": _finite_metric(
                metrics, "metric/prism_module_local_mix_mean"
            ),
            "threshold_mean": _finite_metric(
                metrics, "metric/prism_module_local_threshold_mean"
            ),
            "slope_mean": _finite_metric(
                metrics, "metric/prism_module_local_slope_mean"
            ),
            "gain_mean": _finite_metric(
                metrics, "metric/prism_module_local_gain_mean"
            ),
            "linear_scale_mean": _finite_metric(
                metrics, "metric/prism_module_local_linear_scale_mean"
            ),
        }
    else:
        record["paper_compartmental_branch"] = {"enabled": False}

    rank_penalty_active = False
    if rank_gate is None:
        record["pathology_rank"] = {"enabled": False}
    else:
        diagnostics = dict(rank_gate.diagnostics())
        multiplier = float(rank_gate.sparsity_multiplier())
        rank_penalty_active = multiplier > 0.0
        diagnostics.update(
            {
                "enabled": True,
                "sparsity_multiplier": multiplier,
            }
        )
        record["pathology_rank"] = diagnostics

    generator_objective_active = False
    if generator_gate is None or generator_config is None or not bool(
        getattr(generator_config, "enabled", False)
    ):
        record["generator_count"] = {"enabled": False}
    else:
        mode = _generator_mode(generator_config, int(epoch))
        start = int(generator_config.start_epoch)
        end = max(int(generator_config.base_end_epoch), start)
        progress = (
            0.0
            if int(epoch) <= start
            else min(1.0, (int(epoch) - start) / max(end - start, 1))
        )
        temperature = float(generator_gate.temperature(progress))
        probability = generator_gate.keep_probability(
            temperature=temperature
        ).detach().float()
        candidate = generator_gate.deterministic_mask(
            temperature=temperature
        ).detach().float()
        generator_objective_active = int(epoch) >= start
        live_hard_count: int | None
        if mode in {"warmup_all_on", "shadow_all_on"}:
            live_hard_count = int(generator_gate.num_generators)
        elif mode == "protected_hard":
            live_hard_count = int(candidate.sum().item())
        else:
            live_hard_count = None
        record["generator_count"] = {
            "enabled": True,
            "mode": mode,
            "candidate_count": int(generator_gate.num_generators),
            "live_hard_count": live_hard_count,
            "deterministic_candidate_hard_count": int(candidate.sum().item()),
            "expected_active": float(probability.sum().cpu()),
            "uncertain_fraction": float(
                ((probability > 0.1) & (probability < 0.9))
                .float()
                .mean()
                .cpu()
            ),
            "temperature": temperature,
            "full_nll_violation": _finite_metric(
                metrics, "metric/generator_violation_full_nll"
            ),
            "isolated_nll_violation": _finite_metric(
                metrics, "metric/generator_violation_isolated_nll"
            ),
        }

    overlap = bool(rank_penalty_active and generator_objective_active)
    record["coordination_guard"] = {
        "rank_cardinality_active": rank_penalty_active,
        "generator_objective_active": generator_objective_active,
        "simultaneous_cardinality": overlap,
        "status": "violation" if overlap else "separated",
    }
    record["validation_context"] = {
        "total_loss": _finite_metric(metrics, "loss/total"),
        "prism_branch_nll": _finite_metric(metrics, "loss/prism_branch_nll"),
        "module_rescue_full_ccc": _finite_metric(
            metrics, "metric/prism_module_rescue_full_flattened_ccc"
        ),
        "phu_relative_ess": _finite_metric(
            metrics, "metric/phu_rel_ess_fraction"
        ),
    }
    return record


def format_architecture_capacity_record(record: Mapping[str, Any]) -> str:
    """Render a compact Korean dashboard suitable for a live tmux tail."""

    epoch = int(record["epoch"])
    total = int(record["total_epochs"])
    route = record["pathology_route"]
    rank = record["pathology_rank"]
    generator = record["generator_count"]
    guard = record["coordination_guard"]
    paper = record.get("paper_compartmental_branch", {"enabled": False})

    def value_text(value: Any, digits: int = 3) -> str:
        if value is None:
            return "-"
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)
        return f"{number:.{digits}f}" if math.isfinite(number) else "-"

    lines = [
        f"[architecture-capacity epoch {epoch:03d}/{total:03d}]",
        "  pathology route  "
        f"scale={value_text(route.get('scale'))} "
        f"absmean={value_text(route.get('decoder_absmean'), 5)} "
        f"aux={value_text(route.get('auxiliary_loss'), 4)}",
    ]
    if bool(rank.get("enabled", False)):
        lines.append(
            "  pathology rank   "
            f"live hard/expected={int(rank['hard_rank'])}/"
            f"{value_text(rank['expected_rank'], 2)} of {int(rank['capacity'])} "
            f"search hard/expected={int(rank['search_hard_rank'])}/"
            f"{value_text(rank['search_expected_rank'], 2)} "
            f"mode={rank['mode']} tau={value_text(rank['temperature'])} "
            f"sparse={value_text(rank['sparsity_multiplier'], 2)}"
        )
        lines.append(
            "  rank commit      "
            f"low/mid/high={int(rank['freeze_count_at_low_threshold'])}/"
            f"{int(rank['freeze_count_at_hard_threshold'])}/"
            f"{int(rank['freeze_count_at_high_threshold'])} "
            f"spread={int(rank['freeze_count_spread'])} "
            f"near={100.0 * float(rank['freeze_near_threshold_uncertain_fraction']):.1f}% "
            f"ready={bool(rank['freeze_ready'])} finalized={bool(rank['finalized'])}"
        )
    else:
        lines.append("  pathology rank   disabled")

    if bool(paper.get("enabled", False)):
        lines.append(
            "  paper branch     "
            f"variant={paper['variant']} "
            f"ramp={value_text(paper.get('ramp'), 2)} "
            f"nonlinear/local={value_text(paper.get('nonlinear_rms'), 5)}/"
            f"{value_text(paper.get('local_rms'), 5)} "
            f"ratio={value_text(paper.get('nonlinear_to_local_ratio'))} "
            f"crossing={value_text(paper.get('threshold_crossing_fraction'))}"
        )
        lines.append(
            "  paper OFF audit  "
            f"branch/full dNLL="
            f"{value_text(paper.get('branch_nll_gain_when_enabled'), 5)}/"
            f"{value_text(paper.get('full_nll_gain_when_enabled'), 5)} "
            "positive=helpful"
        )
    else:
        lines.append("  paper branch     disabled")

    if bool(generator.get("enabled", False)):
        live = generator.get("live_hard_count")
        live_text = "soft" if live is None else str(int(live))
        lines.append(
            "  generators      "
            f"live={live_text}/{int(generator['candidate_count'])} "
            f"candidate={int(generator['deterministic_candidate_hard_count'])} "
            f"expected={value_text(generator['expected_active'], 2)} "
            f"mode={generator['mode']} tau={value_text(generator['temperature'])} "
            f"uncertain={100.0 * float(generator['uncertain_fraction']):.1f}%"
        )
    else:
        lines.append("  generators       disabled")
    lines.append(
        "  role separation  "
        f"rank_penalty={'ON' if guard['rank_cardinality_active'] else 'OFF'} "
        f"generator_objective={'ON' if guard['generator_objective_active'] else 'OFF'} "
        f"overlap={'YES' if guard['simultaneous_cardinality'] else 'NO'} "
        f"status={guard['status']}"
    )
    return "\n".join(lines)


__all__ = [
    "build_architecture_capacity_record",
    "format_architecture_capacity_record",
]
