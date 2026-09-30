"""Deterministic hard-mask controller for generator birth/death search.

The controller intentionally does *not* learn continuous gates and does not
apply L1/L0 penalties.  A generator is either present or absent.  Candidate
models are evaluated externally with exact counterfactual forward passes, and
this module decides whether an atomic block birth/death is admissible.

The decision policy is:

1. every candidate must satisfy absolute metric limits;
2. every candidate metric must be non-inferior to the incumbent;
3. a death is useful by reducing cardinality once (1) and (2) hold;
4. a birth must additionally give a lexicographic metric improvement;
5. the same admissible proposal must pass on consecutive checks before the
   hard active set changes (hysteresis);
6. protected generators can never be removed.

Metric values, not parameter magnitudes, drive the decision.  The caller is
responsible for producing leakage-safe train/inner-validation measurements
and for freezing bypass routes while measuring generator counterfactuals.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence, cast


Direction = Literal["maximize", "minimize"]
ProposalKind = Literal["birth", "death"]

_COMPONENT_NAME = "hard_generator_budget_controller"
_SCHEMA_VERSION = 1


def _integer(name: str, value: object, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, got {type(value).__name__}.")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {result}.")
    return result


def _finite_nonnegative(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be numeric, got {type(value).__name__}.")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative, got {result}.")
    return result


def _normalise_ids(
    name: str,
    values: Sequence[int],
    *,
    num_generators: int,
    allow_empty: bool,
) -> tuple[int, ...]:
    raw = tuple(values)
    if not allow_empty and not raw:
        raise ValueError(f"{name} must contain at least one generator ID.")
    checked = tuple(_integer(f"{name} item", item) for item in raw)
    if len(set(checked)) != len(checked):
        raise ValueError(f"{name} contains duplicate generator IDs.")
    out_of_range = tuple(
        item for item in checked if item < 0 or item >= num_generators
    )
    if out_of_range:
        raise ValueError(
            f"{name} contains IDs outside [0, {num_generators}): "
            f"{out_of_range}."
        )
    return tuple(sorted(checked))


@dataclass(frozen=True)
class MetricRule:
    """One ordered safety/non-inferiority rule.

    Parameters
    ----------
    name:
        Key expected in incumbent and candidate metric dictionaries.
    direction:
        Whether larger or smaller values are preferable.
    hard_limit:
        Optional absolute floor (``maximize``) or ceiling (``minimize``).
    noninferiority_margin:
        Maximum allowed worsening relative to the incumbent.
    improvement_margin:
        Change that must be exceeded to be lexicographically decisive for a
        birth.  Values inside ``+/- improvement_margin`` are treated as tied.
    atol:
        Numerical comparison tolerance.

    Notes
    -----
    Rules are evaluated in the order provided to the controller.  This order
    is the lexicographic priority order; it is therefore part of the manifest.
    """

    name: str
    direction: Direction
    hard_limit: float | None = None
    noninferiority_margin: float = 0.0
    improvement_margin: float = 0.0
    atol: float = 1e-12

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("MetricRule.name must be a non-empty string.")
        if self.direction not in ("maximize", "minimize"):
            raise ValueError(
                "MetricRule.direction must be 'maximize' or 'minimize', "
                f"got {self.direction!r}."
            )
        if self.hard_limit is not None:
            if isinstance(self.hard_limit, bool) or not isinstance(
                self.hard_limit, (int, float)
            ):
                raise TypeError("MetricRule.hard_limit must be numeric or None.")
            if not math.isfinite(float(self.hard_limit)):
                raise ValueError("MetricRule.hard_limit must be finite.")
        object.__setattr__(
            self,
            "noninferiority_margin",
            _finite_nonnegative(
                "MetricRule.noninferiority_margin",
                self.noninferiority_margin,
            ),
        )
        object.__setattr__(
            self,
            "improvement_margin",
            _finite_nonnegative(
                "MetricRule.improvement_margin", self.improvement_margin
            ),
        )
        object.__setattr__(
            self,
            "atol",
            _finite_nonnegative("MetricRule.atol", self.atol),
        )
        if self.hard_limit is not None:
            object.__setattr__(self, "hard_limit", float(self.hard_limit))

    def oriented_gain(self, incumbent: float, candidate: float) -> float:
        """Return positive values when the candidate is better."""

        if self.direction == "maximize":
            return candidate - incumbent
        return incumbent - candidate

    def satisfies_hard_limit(self, value: float) -> bool:
        if self.hard_limit is None:
            return True
        if self.direction == "maximize":
            return value + self.atol >= self.hard_limit
        return value - self.atol <= self.hard_limit

    def incumbent_violates_hard_limit(self, value: float) -> bool:
        return not self.satisfies_hard_limit(value)


@dataclass(frozen=True)
class GeneratorBudgetDecision:
    """Auditable result of one proposal check."""

    step: int
    kind: ProposalKind
    generator_ids: tuple[int, ...]
    eligible: bool
    applied: bool
    confirmation_count: int
    required_confirmations: int
    active_count_before: int
    active_count_after: int
    reason_codes: tuple[str, ...]
    first_decisive_metric: str | None
    oriented_metric_gains: tuple[tuple[str, float], ...]

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["generator_ids"] = list(self.generator_ids)
        result["reason_codes"] = list(self.reason_codes)
        result["oriented_metric_gains"] = [
            [name, value] for name, value in self.oriented_metric_gains
        ]
        return result

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "GeneratorBudgetDecision":
        return cls(
            step=int(raw["step"]),
            kind=str(raw["kind"]),  # type: ignore[arg-type]
            generator_ids=tuple(int(x) for x in raw["generator_ids"]),
            eligible=bool(raw["eligible"]),
            applied=bool(raw["applied"]),
            confirmation_count=int(raw["confirmation_count"]),
            required_confirmations=int(raw["required_confirmations"]),
            active_count_before=int(raw["active_count_before"]),
            active_count_after=int(raw["active_count_after"]),
            reason_codes=tuple(str(x) for x in raw["reason_codes"]),
            first_decisive_metric=(
                None
                if raw["first_decisive_metric"] is None
                else str(raw["first_decisive_metric"])
            ),
            oriented_metric_gains=tuple(
                (str(name), float(value))
                for name, value in raw["oriented_metric_gains"]
            ),
        )


class HardGeneratorBudgetController:
    """Maintain and update a hard boolean generator active set.

    This class does not rank generators and never reads model weights.  The
    trainer or analysis driver proposes an exact block and supplies paired
    incumbent/candidate metrics.  Calling :meth:`consider` is the sole state
    transition operation.
    """

    def __init__(
        self,
        *,
        num_generators: int,
        metric_rules: Sequence[MetricRule],
        active_ids: Sequence[int] | None = None,
        protected_ids: Sequence[int] = (),
        minimum_active_generators: int = 1,
        maximum_active_generators: int | None = None,
        required_confirmations: int = 2,
        cooldown_checks: int = 1,
    ) -> None:
        self.num_generators = _integer(
            "num_generators", num_generators, minimum=1
        )
        if not metric_rules:
            raise ValueError("metric_rules must contain at least one rule.")
        self.metric_rules = tuple(metric_rules)
        if not all(isinstance(rule, MetricRule) for rule in self.metric_rules):
            raise TypeError("Every metric_rules item must be a MetricRule.")
        names = tuple(rule.name for rule in self.metric_rules)
        if len(set(names)) != len(names):
            raise ValueError("metric_rules names must be unique.")

        self.protected_ids = _normalise_ids(
            "protected_ids",
            protected_ids,
            num_generators=self.num_generators,
            allow_empty=True,
        )
        if active_ids is None:
            normalised_active = tuple(range(self.num_generators))
        else:
            normalised_active = _normalise_ids(
                "active_ids",
                active_ids,
                num_generators=self.num_generators,
                allow_empty=True,
            )
        missing_protected = tuple(
            item
            for item in self.protected_ids
            if item not in set(normalised_active)
        )
        if missing_protected:
            raise ValueError(
                "Every protected generator must initially be active; missing "
                f"{missing_protected}."
            )

        self.minimum_active_generators = _integer(
            "minimum_active_generators",
            minimum_active_generators,
            minimum=1,
        )
        if maximum_active_generators is None:
            self.maximum_active_generators = self.num_generators
        else:
            self.maximum_active_generators = _integer(
                "maximum_active_generators",
                maximum_active_generators,
                minimum=1,
            )
        if self.minimum_active_generators > self.maximum_active_generators:
            raise ValueError(
                "minimum_active_generators cannot exceed "
                "maximum_active_generators."
            )
        if self.maximum_active_generators > self.num_generators:
            raise ValueError(
                "maximum_active_generators cannot exceed num_generators."
            )
        if not (
            self.minimum_active_generators
            <= len(normalised_active)
            <= self.maximum_active_generators
        ):
            raise ValueError(
                "Initial active count must lie within the configured bounds."
            )
        if len(self.protected_ids) > self.maximum_active_generators:
            raise ValueError(
                "The protected set exceeds maximum_active_generators."
            )

        self.required_confirmations = _integer(
            "required_confirmations", required_confirmations, minimum=1
        )
        self.cooldown_checks = _integer(
            "cooldown_checks", cooldown_checks, minimum=0
        )

        self._active = [False] * self.num_generators
        for item in normalised_active:
            self._active[item] = True
        self._pending_signature: tuple[ProposalKind, tuple[int, ...]] | None = (
            None
        )
        self._pending_confirmations = 0
        self._last_evaluation_step: int | None = None
        self._last_change_step: int | None = None
        self._history: list[GeneratorBudgetDecision] = []

    @property
    def active_mask(self) -> tuple[bool, ...]:
        """Immutable *hard* boolean mask; no fractional values are possible."""

        return tuple(self._active)

    @property
    def active_ids(self) -> tuple[int, ...]:
        return tuple(i for i, active in enumerate(self._active) if active)

    @property
    def inactive_ids(self) -> tuple[int, ...]:
        return tuple(i for i, active in enumerate(self._active) if not active)

    @property
    def history(self) -> tuple[GeneratorBudgetDecision, ...]:
        return tuple(self._history)

    def _metrics(
        self, name: str, values: Mapping[str, float]
    ) -> dict[str, float]:
        missing = tuple(
            rule.name for rule in self.metric_rules if rule.name not in values
        )
        if missing:
            raise KeyError(f"{name} is missing required metrics: {missing}.")
        result: dict[str, float] = {}
        for rule in self.metric_rules:
            value = values[rule.name]
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(
                    f"{name}[{rule.name!r}] must be numeric, got "
                    f"{type(value).__name__}."
                )
            result[rule.name] = float(value)
        return result

    def _validate_proposal(
        self, kind: ProposalKind, generator_ids: Sequence[int]
    ) -> tuple[int, ...]:
        if kind not in ("birth", "death"):
            raise ValueError(
                f"kind must be 'birth' or 'death', got {kind!r}."
            )
        block = _normalise_ids(
            "generator_ids",
            generator_ids,
            num_generators=self.num_generators,
            allow_empty=False,
        )
        if kind == "death":
            protected = tuple(
                item for item in block if item in set(self.protected_ids)
            )
            if protected:
                raise ValueError(
                    f"A death proposal contains protected IDs: {protected}."
                )
            inactive = tuple(item for item in block if not self._active[item])
            if inactive:
                raise ValueError(
                    f"A death proposal contains inactive IDs: {inactive}."
                )
            after = len(self.active_ids) - len(block)
            if after < self.minimum_active_generators:
                raise ValueError(
                    "Death proposal would violate minimum_active_generators: "
                    f"{after} < {self.minimum_active_generators}."
                )
        else:
            active = tuple(item for item in block if self._active[item])
            if active:
                raise ValueError(
                    f"A birth proposal contains active IDs: {active}."
                )
            after = len(self.active_ids) + len(block)
            if after > self.maximum_active_generators:
                raise ValueError(
                    "Birth proposal would violate maximum_active_generators: "
                    f"{after} > {self.maximum_active_generators}."
                )
        return block

    def _cooldown_active(self, step: int) -> bool:
        if self._last_change_step is None:
            return False
        return step - self._last_change_step <= self.cooldown_checks

    def _record(
        self,
        *,
        step: int,
        kind: ProposalKind,
        block: tuple[int, ...],
        eligible: bool,
        applied: bool,
        confirmations: int,
        active_count_before: int,
        reasons: tuple[str, ...],
        first_decisive_metric: str | None,
        gains: tuple[tuple[str, float], ...],
    ) -> GeneratorBudgetDecision:
        # Non-finite proposals are fail-closed above.  Keep their explicit
        # ``nonfinite:<metric>`` reason code, but store a finite placeholder
        # gain so the strict-JSON resumable audit trail itself cannot crash.
        json_safe_gains = tuple(
            (name, float(gain) if math.isfinite(float(gain)) else 0.0)
            for name, gain in gains
        )
        decision = GeneratorBudgetDecision(
            step=step,
            kind=kind,
            generator_ids=block,
            eligible=eligible,
            applied=applied,
            confirmation_count=confirmations,
            required_confirmations=self.required_confirmations,
            active_count_before=active_count_before,
            active_count_after=len(self.active_ids),
            reason_codes=reasons,
            first_decisive_metric=first_decisive_metric,
            oriented_metric_gains=json_safe_gains,
        )
        self._history.append(decision)
        self._last_evaluation_step = step
        return decision

    def consider(
        self,
        *,
        step: int,
        kind: ProposalKind,
        generator_ids: Sequence[int],
        incumbent_metrics: Mapping[str, float],
        candidate_metrics: Mapping[str, float],
        external_reason_codes: Sequence[str] = (),
    ) -> GeneratorBudgetDecision:
        """Evaluate one atomic birth/death proposal and possibly apply it.

        ``step`` must increase strictly between calls.  This prevents a caller
        from satisfying hysteresis by submitting the same measurements twice.
        A different proposal, any failed safety check, or cooldown resets the
        pending confirmation streak.
        """

        checked_step = _integer("step", step)
        if (
            self._last_evaluation_step is not None
            and checked_step <= self._last_evaluation_step
        ):
            raise ValueError(
                "step must increase strictly between consider() calls; "
                f"last={self._last_evaluation_step}, got {checked_step}."
            )
        block = self._validate_proposal(kind, generator_ids)
        incumbent = self._metrics("incumbent_metrics", incumbent_metrics)
        candidate = self._metrics("candidate_metrics", candidate_metrics)
        active_before = len(self.active_ids)

        gains = tuple(
            (
                rule.name,
                rule.oriented_gain(
                    incumbent[rule.name], candidate[rule.name]
                ),
            )
            for rule in self.metric_rules
        )

        if self._cooldown_active(checked_step):
            self._pending_signature = None
            self._pending_confirmations = 0
            return self._record(
                step=checked_step,
                kind=kind,
                block=block,
                eligible=False,
                applied=False,
                confirmations=0,
                active_count_before=active_before,
                reasons=("cooldown",),
                first_decisive_metric=None,
                gains=gains,
            )

        reasons: list[str] = []
        for reason in external_reason_codes:
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError(
                    "external_reason_codes must contain non-empty strings."
                )
            reasons.append(reason.strip())
        for rule, (_, gain) in zip(self.metric_rules, gains):
            incumbent_value = incumbent[rule.name]
            candidate_value = candidate[rule.name]
            if not (
                math.isfinite(incumbent_value)
                and math.isfinite(candidate_value)
                and math.isfinite(gain)
            ):
                reasons.append(f"nonfinite:{rule.name}")
                continue
            if not rule.satisfies_hard_limit(candidate_value):
                reasons.append(f"hard_limit:{rule.name}")
            if gain + rule.atol < -rule.noninferiority_margin:
                reasons.append(f"noninferiority:{rule.name}")

        first_decisive_metric: str | None = None
        lexicographic_result = 0
        repaired_hard_limit = False
        if not reasons:
            for rule, (_, gain) in zip(self.metric_rules, gains):
                if (
                    rule.incumbent_violates_hard_limit(incumbent[rule.name])
                    and rule.satisfies_hard_limit(candidate[rule.name])
                ):
                    first_decisive_metric = rule.name
                    lexicographic_result = 1
                    repaired_hard_limit = True
                    break
                if gain > rule.improvement_margin + rule.atol:
                    first_decisive_metric = rule.name
                    lexicographic_result = 1
                    break
                if gain < -rule.improvement_margin - rule.atol:
                    first_decisive_metric = rule.name
                    lexicographic_result = -1
                    break

            if kind == "birth" and lexicographic_result <= 0:
                if lexicographic_result < 0:
                    reasons.append(
                        "birth_lexicographic_harm:"
                        f"{first_decisive_metric}"
                    )
                else:
                    reasons.append("birth_no_lexicographic_improvement")

        if reasons:
            self._pending_signature = None
            self._pending_confirmations = 0
            return self._record(
                step=checked_step,
                kind=kind,
                block=block,
                eligible=False,
                applied=False,
                confirmations=0,
                active_count_before=active_before,
                reasons=tuple(reasons),
                first_decisive_metric=first_decisive_metric,
                gains=gains,
            )

        signature = (kind, block)
        if self._pending_signature == signature:
            self._pending_confirmations += 1
        else:
            self._pending_signature = signature
            self._pending_confirmations = 1

        confirmations = self._pending_confirmations
        if confirmations < self.required_confirmations:
            return self._record(
                step=checked_step,
                kind=kind,
                block=block,
                eligible=True,
                applied=False,
                confirmations=confirmations,
                active_count_before=active_before,
                reasons=("awaiting_hysteresis",),
                first_decisive_metric=first_decisive_metric,
                gains=gains,
            )

        target_value = kind == "birth"
        for item in block:
            self._active[item] = target_value
        self._last_change_step = checked_step
        self._pending_signature = None
        self._pending_confirmations = 0

        reason = (
            "applied:repaired_hard_limit"
            if repaired_hard_limit
            else "applied"
        )
        return self._record(
            step=checked_step,
            kind=kind,
            block=block,
            eligible=True,
            applied=True,
            confirmations=confirmations,
            active_count_before=active_before,
            reasons=(reason,),
            first_decisive_metric=first_decisive_metric,
            gains=gains,
        )

    def _config_dict(self) -> dict[str, Any]:
        return {
            "num_generators": self.num_generators,
            "metric_rules": [asdict(rule) for rule in self.metric_rules],
            "protected_ids": list(self.protected_ids),
            "minimum_active_generators": self.minimum_active_generators,
            "maximum_active_generators": self.maximum_active_generators,
            "required_confirmations": self.required_confirmations,
            "cooldown_checks": self.cooldown_checks,
        }

    def config_fingerprint(self) -> str:
        canonical = json.dumps(
            self._config_dict(),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def state_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable, versioned controller checkpoint."""

        return {
            "component": _COMPONENT_NAME,
            "schema_version": _SCHEMA_VERSION,
            "config": self._config_dict(),
            "config_sha256": self.config_fingerprint(),
            "state": {
                "active_ids": list(self.active_ids),
                "pending": (
                    None
                    if self._pending_signature is None
                    else {
                        "kind": self._pending_signature[0],
                        "generator_ids": list(self._pending_signature[1]),
                        "confirmations": self._pending_confirmations,
                    }
                ),
                "last_evaluation_step": self._last_evaluation_step,
                "last_change_step": self._last_change_step,
                "history": [item.to_dict() for item in self._history],
            },
        }

    @classmethod
    def from_state_dict(
        cls, raw: Mapping[str, Any]
    ) -> "HardGeneratorBudgetController":
        """Restore a controller and validate its configuration fingerprint."""

        if raw.get("component") != _COMPONENT_NAME:
            raise ValueError("State component name does not match controller.")
        if raw.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError(
                "Unsupported controller state schema version: "
                f"{raw.get('schema_version')!r}."
            )
        config = dict(raw["config"])
        rules = tuple(
            MetricRule(**dict(rule)) for rule in config.pop("metric_rules")
        )
        state = dict(raw["state"])
        controller = cls(
            metric_rules=rules,
            active_ids=state["active_ids"],
            **config,
        )
        expected_fingerprint = str(raw["config_sha256"])
        if controller.config_fingerprint() != expected_fingerprint:
            raise ValueError(
                "Controller configuration fingerprint mismatch; the state "
                "may have been edited or corrupted."
            )

        pending = state.get("pending")
        if pending is not None:
            pending_kind = str(pending["kind"])
            if pending_kind not in ("birth", "death"):
                raise ValueError("Invalid pending proposal kind in state.")
            pending_ids = _normalise_ids(
                "pending.generator_ids",
                pending["generator_ids"],
                num_generators=controller.num_generators,
                allow_empty=False,
            )
            controller._pending_signature = (  # noqa: SLF001
                cast(ProposalKind, pending_kind),
                pending_ids,
            )
            controller._pending_confirmations = _integer(  # noqa: SLF001
                "pending.confirmations",
                pending["confirmations"],
                minimum=1,
            )
        controller._last_evaluation_step = state.get(  # noqa: SLF001
            "last_evaluation_step"
        )
        controller._last_change_step = state.get(  # noqa: SLF001
            "last_change_step"
        )
        controller._history = [  # noqa: SLF001
            GeneratorBudgetDecision.from_dict(item)
            for item in state.get("history", ())
        ]
        return controller

    def manifest(self) -> dict[str, Any]:
        """Return a compact JSON-safe audit manifest for a run artifact."""

        return {
            "component": _COMPONENT_NAME,
            "schema_version": _SCHEMA_VERSION,
            "selection_type": "hard_boolean_birth_death",
            "uses_soft_gates": False,
            "uses_l1_or_l0_penalty": False,
            "config_sha256": self.config_fingerprint(),
            "num_generators": self.num_generators,
            "active_count": len(self.active_ids),
            "active_ids": list(self.active_ids),
            "protected_ids": list(self.protected_ids),
            "metric_priority": [rule.name for rule in self.metric_rules],
            "required_confirmations": self.required_confirmations,
            "cooldown_checks": self.cooldown_checks,
            "last_evaluation_step": self._last_evaluation_step,
            "last_change_step": self._last_change_step,
            "num_decisions": len(self._history),
        }

    def write_state_json(self, path: str | Path) -> None:
        """Atomically write the complete resumable state as strict JSON."""

        destination = Path(path)
        temporary = destination.with_name(destination.name + ".tmp")
        payload = json.dumps(
            self.state_dict(),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        temporary.write_text(payload + "\n", encoding="utf-8")
        temporary.replace(destination)

    @classmethod
    def read_state_json(
        cls, path: str | Path
    ) -> "HardGeneratorBudgetController":
        return cls.from_state_dict(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )
