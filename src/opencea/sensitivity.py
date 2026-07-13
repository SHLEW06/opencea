"""One-way deterministic sensitivity analysis (DSA).

For each parameter in turn, the DSA sets that parameter to its low value
then its high value (all others held at base case), rebuilds the model
through :func:`opencea.builders.build_darth_sick_sicker`, runs the
validated deterministic engine, and records the outcome at each end.

The outcome of interest is the **incremental net monetary benefit** of a
chosen comparator versus the baseline (SoC by default) at a fixed
willingness-to-pay threshold. NMB is chosen over the ICER because ICERs
flip sign and go undefined when sweeps push incremental QALYs through
zero, making them unusable as a tornado outcome.

Parameter low / high values are pulled from the percentiles of each
parameter's existing PSA marginal in :data:`opencea.psa.PSA_PARAM_SPECS`,
so the DSA ranges descend from the same DARTH-cited distributions as the
PSA — no fabricated bounds. The default percentiles are 2.5 / 97.5; the
range is extended to include the deterministic base case when the base
lies outside the percentile interval (relevant for ``u_H``, whose
deterministic value of 1.0 sits above the Beta(200, 3) 97.5th percentile).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ._types import DSAOutcome, ParameterInput
from .builders import build_darth_sick_sicker
from .cea import nmb as _nmb
from .engine import run_model
from .psa import PSA_PARAM_SPECS, DistSpec, sample_psa_params

EvaluatorFn = Callable[[Mapping[str, Any]], Mapping[str, DSAOutcome]]


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParameterSweep:
    """Result of sweeping one parameter from its low to its high value.

    ``low_outcome`` and ``high_outcome`` are the outcome when the parameter
    is held at ``low_value`` / ``high_value`` respectively (not necessarily
    the minimum / maximum outcome — the outcome function may decrease in
    the parameter).
    """

    parameter: str
    base_value: float
    low_value: float
    high_value: float
    base_outcome: float
    low_outcome: float
    high_outcome: float

    @property
    def swing(self) -> float:
        """Absolute width of the tornado bar for this parameter."""
        return float(abs(self.high_outcome - self.low_outcome))


@dataclass
class DSAResult:
    """Tidy result of a one-way DSA."""

    base_outcome: float
    wtp: float
    comparator: str
    baseline: str
    sweeps: List[ParameterSweep]
    base_params: Mapping[str, Any]
    percentiles: Tuple[float, float]

    def to_dataframe(self) -> pd.DataFrame:
        """Long DataFrame with one row per swept parameter, sorted by swing."""
        return pd.DataFrame(
            [
                {
                    "parameter": s.parameter,
                    "base_value": s.base_value,
                    "low_value": s.low_value,
                    "high_value": s.high_value,
                    "low_outcome": s.low_outcome,
                    "high_outcome": s.high_outcome,
                    "base_outcome": s.base_outcome,
                    "swing": s.swing,
                }
                for s in self.sweeps
            ]
        )


# ---------------------------------------------------------------------------
# Parameter ranges from PSA marginals
# ---------------------------------------------------------------------------


def compute_parameter_ranges(
    base_params: Mapping[str, Any],
    specs: Optional[Mapping[str, DistSpec]] = None,
    percentiles: Tuple[float, float] = (2.5, 97.5),
    n_sample: int = 200_000,
    seed: int = 20260625,
    extend_to_include_base: bool = True,
) -> Dict[str, Tuple[float, float]]:
    """Compute one-way DSA low/high values from PSA marginal percentiles.

    The percentiles are taken empirically from a single seeded draw of the
    existing :data:`opencea.psa.PSA_PARAM_SPECS` distributions (200k samples
    by default, so the percentile is stable to ~4 sig figs). When
    ``extend_to_include_base`` is True, each range is clamped to include
    the deterministic base value — this keeps the sweep bracketing the
    base case even for parameters like ``u_H`` whose deterministic value
    sits outside the PSA distribution's support tail.
    """
    specs = specs if specs is not None else PSA_PARAM_SPECS
    draws = sample_psa_params(n_sim=n_sample, seed=seed, specs=specs)

    p_low, p_high = percentiles
    out: Dict[str, Tuple[float, float]] = {}
    for name in specs.keys():
        lo, hi = np.percentile(draws[name].to_numpy(), [p_low, p_high])
        lo, hi = float(lo), float(hi)
        if extend_to_include_base and name in base_params:
            base = float(base_params[name])
            lo = min(lo, base)
            hi = max(hi, base)
        out[name] = (lo, hi)
    return out


# ---------------------------------------------------------------------------
# Outcome helpers
# ---------------------------------------------------------------------------


def _eval_strategies(params: Mapping[str, Any]) -> dict[str, DSAOutcome]:
    """Build + run the DARTH model; return ``{name: {cost, qaly}}``."""
    model = build_darth_sick_sicker(params)
    return {
        r["name"]: {"cost": float(r["total_cost"]), "qaly": float(r["total_qaly"])}
        for r in run_model(model)
    }


def _inc_nmb(
    results: Mapping[str, DSAOutcome],
    comparator: str,
    baseline: str,
    wtp: float,
) -> float:
    base = results[baseline]
    comp = results[comparator]
    return _nmb(comp["cost"], comp["qaly"], wtp) - _nmb(base["cost"], base["qaly"], wtp)


def _optimal_comparator(
    results: Mapping[str, DSAOutcome],
    baseline: str,
    wtp: float,
) -> str:
    """Strategy with highest NMB at ``wtp``, excluding ``baseline``."""
    candidates = [n for n in results if n != baseline]
    return max(
        candidates, key=lambda n: _nmb(results[n]["cost"], results[n]["qaly"], wtp)
    )


# ---------------------------------------------------------------------------
# Public DSA driver
# ---------------------------------------------------------------------------


def _load_base_params(base_params: ParameterInput) -> Dict[str, Any]:
    if isinstance(base_params, (str, Path)):
        import yaml

        with open(base_params, "r") as f:
            return dict(yaml.safe_load(f))
    return dict(base_params)


def run_dsa(
    base_params: ParameterInput,
    wtp: float = 100_000.0,
    comparator: Optional[str] = None,
    baseline: str = "Standard of care",
    percentiles: Tuple[float, float] = (2.5, 97.5),
    sweep_params: Optional[Iterable[str]] = None,
    n_sample: int = 200_000,
    seed: int = 20260625,
    *,
    evaluator: Optional[EvaluatorFn] = None,
    param_specs: Optional[Mapping[str, DistSpec]] = None,
    param_ranges: Optional[Mapping[str, Tuple[float, float]]] = None,
) -> DSAResult:
    """Run a one-way DSA on a cohort cost-effectiveness model.

    Defaults to the DARTH Sick-Sicker setup; reusable on any model that
    can produce a ``{strategy_name: {"cost": float, "qaly": float}}``
    dictionary from a parameter mapping.

    Parameters
    ----------
    base_params
        Mapping or YAML path with the deterministic base parameter set.
    wtp
        Willingness-to-pay threshold for the NMB outcome (default $100k).
    comparator
        Strategy whose incremental NMB vs ``baseline`` is the outcome.
        Defaults to the strategy with the highest base-case NMB at
        ``wtp`` (excluding ``baseline``) — the expected-NMB-optimal
        strategy at base case.
    baseline
        Comparator's reference strategy (default "Standard of care").
    percentiles
        ``(low, high)`` percentiles of the PSA marginal used as the
        parameter sweep bounds (default 2.5 / 97.5).
    sweep_params
        Iterable of parameter names to sweep. Defaults to every key in
        ``param_specs``.
    n_sample, seed
        Sample size and seed for the empirical percentile estimate.
    evaluator
        Optional callable that takes a parameter mapping and returns
        ``{name: {"cost": float, "qaly": float}}``. Defaults to building
        the DARTH model + running the engine. Overriding it lets the
        same DSA driver power any cohort model that exposes the same
        outcome shape (e.g., the empagliflozin case study).
    param_specs
        Optional override of :data:`opencea.psa.PSA_PARAM_SPECS` —
        used to derive percentile-based ranges for parameters not
        explicitly given in ``param_ranges``.
    param_ranges
        Optional per-parameter ``(low, high)`` overrides. Takes
        precedence over any percentile-derived range from
        ``param_specs``. Parameters with explicit literature bounds
        (e.g., trial 95% CIs) are typically passed here.
    """
    bp = _load_base_params(base_params)
    evaluator = evaluator if evaluator is not None else _eval_strategies
    param_specs = param_specs if param_specs is not None else PSA_PARAM_SPECS

    base_results = evaluator(bp)
    if comparator is None:
        comparator = _optimal_comparator(base_results, baseline, wtp)
    if comparator not in base_results:
        raise KeyError(f"comparator {comparator!r} not among {list(base_results)}")
    if baseline not in base_results:
        raise KeyError(f"baseline {baseline!r} not among {list(base_results)}")

    base_outcome = _inc_nmb(base_results, comparator, baseline, wtp)

    ranges = compute_parameter_ranges(
        bp,
        specs=param_specs,
        percentiles=percentiles,
        n_sample=n_sample,
        seed=seed,
    )
    if param_ranges is not None:
        # Per-parameter overrides; extend to include the base if needed.
        for name, (lo, hi) in param_ranges.items():
            if name in bp:
                base_v = float(bp[name])
                lo = min(float(lo), base_v)
                hi = max(float(hi), base_v)
            ranges[name] = (float(lo), float(hi))

    names = list(sweep_params) if sweep_params is not None else list(param_specs.keys())

    sweeps: List[ParameterSweep] = []
    for name in names:
        if name not in ranges:
            raise KeyError(f"no range available for {name!r}")
        if name not in bp:
            raise KeyError(f"{name!r} not in base parameters")

        lo_v, hi_v = ranges[name]
        base_v = float(bp[name])

        params_low = dict(bp)
        params_low[name] = lo_v
        out_low = _inc_nmb(evaluator(params_low), comparator, baseline, wtp)

        params_high = dict(bp)
        params_high[name] = hi_v
        out_high = _inc_nmb(evaluator(params_high), comparator, baseline, wtp)

        sweeps.append(
            ParameterSweep(
                parameter=name,
                base_value=base_v,
                low_value=lo_v,
                high_value=hi_v,
                base_outcome=base_outcome,
                low_outcome=out_low,
                high_outcome=out_high,
            )
        )

    sweeps.sort(key=lambda s: s.swing, reverse=True)

    return DSAResult(
        base_outcome=base_outcome,
        wtp=wtp,
        comparator=comparator,
        baseline=baseline,
        sweeps=sweeps,
        base_params=bp,
        percentiles=percentiles,
    )
