"""Probabilistic sensitivity analysis (PSA) for the DARTH Sick-Sicker model.

The distribution table is transcribed verbatim from ``generate_psa_params()``
in the DARTH ``cohort-modeling-tutorial-intro`` repo
(``R/Functions_cSTM_time_indep.R``, lines 278-309). R's ``rgamma`` defaults
to a *rate* parameterization; numpy's :func:`numpy.random.Generator.gamma`
uses *scale*. Every rate-parameterized draw is converted to ``scale = 1/rate``
before being handed to numpy.

The Monte Carlo runner builds all four strategies' transition matrices and
reward vectors per draw using the same competing-risks construction as
:mod:`opencea.builders`, then propagates the cohort trace in batched form:
``trace[t+1] = trace[t] @ P`` becomes a single ``einsum`` over an
``(n_sim, n_states, n_states)`` tensor, so the inner loop runs ``T``
matmuls of width ``n_sim`` rather than ``n_sim * T`` scalar matmuls. The
discount weights and within-cycle correction are produced by the same
:func:`opencea.engine.gen_wcc` used by the deterministic engine.

The result object exposes the sampled parameter draws plus per-strategy
``(n_sim,)`` cost / QALY arrays — the raw material for the CEAC, the
CE plane, and any downstream expected-NMB analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Literal, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ._types import ParameterInput
from .engine import WCCMethod, gen_wcc

# ---------------------------------------------------------------------------
# Distribution specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DistSpec:
    """A single parameter's PSA distribution.

    ``family`` is one of ``gamma``, ``lognormal``, or ``beta``. ``params``
    holds the numpy parameterization (``scale`` for gamma, ``mean``/``sigma``
    for lognormal, ``a``/``b`` for beta). ``target_mean`` is the analytic
    mean (for gamma / beta) or median (for lognormal) used by the
    self-check test that the sampler is wired up correctly.
    """

    family: Literal["gamma", "lognormal", "beta"]
    params: Dict[str, float]
    target_mean: float
    is_median: bool = False  # lognormal target is the median, not the mean


# DARTH generate_psa_params() lines 278-309, converted to numpy parameterization.
#   R rgamma(shape, rate)        -> numpy gamma(shape, scale=1/rate)
#   R rgamma(shape, scale=c)     -> numpy gamma(shape, scale=c)
#   R rlnorm(meanlog, sdlog)     -> numpy lognormal(mean=meanlog, sigma=sdlog)
#   R rbeta(shape1, shape2)      -> numpy beta(a=shape1, b=shape2)
# Fixed (not sampled) per the DARTH function: c_D = 0, u_D = 0.
PSA_PARAM_SPECS: Dict[str, DistSpec] = {
    # transition rates and hazard ratios
    "r_HD": DistSpec("gamma", {"shape": 20.0, "scale": 1.0 / 10000.0}, 20.0 / 10000.0),
    "r_HS1": DistSpec("gamma", {"shape": 30.0, "scale": 1.0 / 200.0}, 30.0 / 200.0),
    "r_S1H": DistSpec("gamma", {"shape": 60.0, "scale": 1.0 / 120.0}, 60.0 / 120.0),
    "r_S1S2": DistSpec("gamma", {"shape": 84.0, "scale": 1.0 / 800.0}, 84.0 / 800.0),
    "hr_S1": DistSpec(
        "lognormal", {"mean": np.log(3.0), "sigma": 0.01}, 3.0, is_median=True
    ),
    "hr_S2": DistSpec(
        "lognormal", {"mean": np.log(10.0), "sigma": 0.02}, 10.0, is_median=True
    ),
    "hr_S1S2_trtB": DistSpec(
        "lognormal", {"mean": np.log(0.6), "sigma": 0.02}, 0.6, is_median=True
    ),
    # costs
    "c_H": DistSpec("gamma", {"shape": 100.0, "scale": 20.0}, 100.0 * 20.0),
    "c_S1": DistSpec("gamma", {"shape": 177.8, "scale": 22.5}, 177.8 * 22.5),
    "c_S2": DistSpec("gamma", {"shape": 225.0, "scale": 66.7}, 225.0 * 66.7),
    "c_trtA": DistSpec("gamma", {"shape": 73.5, "scale": 163.3}, 73.5 * 163.3),
    "c_trtB": DistSpec("gamma", {"shape": 86.2, "scale": 150.8}, 86.2 * 150.8),
    # utilities
    "u_H": DistSpec("beta", {"a": 200.0, "b": 3.0}, 200.0 / 203.0),
    "u_S1": DistSpec("beta", {"a": 130.0, "b": 45.0}, 130.0 / 175.0),
    "u_S2": DistSpec("beta", {"a": 230.0, "b": 230.0}, 0.5),
    "u_trtA": DistSpec("beta", {"a": 300.0, "b": 15.0}, 300.0 / 315.0),
}

# Parameters that are fixed (not sampled), per the DARTH reference.
PSA_FIXED_PARAMS: Dict[str, float] = {"c_D": 0.0, "u_D": 0.0}


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


def sample_psa_params(
    n_sim: int = 10_000,
    seed: int = 20260625,
    specs: Optional[Mapping[str, DistSpec]] = None,
) -> pd.DataFrame:
    """Draw ``n_sim`` PSA parameter sets.

    Uses a single :class:`numpy.random.Generator` seeded with ``seed``; the
    parameters are drawn in dict-insertion order, so identical ``(n_sim,
    seed)`` yields bit-identical results.

    Returns a DataFrame with one row per draw and one column per parameter,
    plus the fixed ``c_D`` / ``u_D`` columns for downstream convenience.
    """
    if n_sim <= 0:
        raise ValueError("n_sim must be positive")
    specs = specs if specs is not None else PSA_PARAM_SPECS

    rng = np.random.default_rng(seed)
    cols: Dict[str, np.ndarray] = {}
    for name, spec in specs.items():
        if spec.family == "gamma":
            cols[name] = rng.gamma(
                shape=spec.params["shape"], scale=spec.params["scale"], size=n_sim
            )
        elif spec.family == "lognormal":
            cols[name] = rng.lognormal(
                mean=spec.params["mean"], sigma=spec.params["sigma"], size=n_sim
            )
        elif spec.family == "beta":
            cols[name] = rng.beta(a=spec.params["a"], b=spec.params["b"], size=n_sim)
        else:
            raise ValueError(f"unknown distribution family: {spec.family!r}")

    for name, value in PSA_FIXED_PARAMS.items():
        cols[name] = np.full(n_sim, value, dtype=float)

    return pd.DataFrame(cols)


# ---------------------------------------------------------------------------
# Vectorized model construction
# ---------------------------------------------------------------------------


# Order matters: indices used for the (n_sim, 4, 4) transition tensor.
_STATES = ("H", "S1", "S2", "D")
_STRATEGY_NAMES = ("Standard of care", "Strategy A", "Strategy B", "Strategy AB")


def _rate_to_prob_vec(r: np.ndarray, t: float = 1.0) -> np.ndarray:
    """Vectorized DARTH rate -> probability: ``1 - exp(-r * t)``."""
    return 1.0 - np.exp(-r * t)


def _build_transition_tensor(
    p_HD: np.ndarray,
    p_HS1: np.ndarray,
    p_S1H: np.ndarray,
    p_S1S2: np.ndarray,
    p_S1D: np.ndarray,
    p_S2D: np.ndarray,
) -> np.ndarray:
    """Construct an ``(n_sim, 4, 4)`` batch of Sick-Sicker transition matrices.

    Mirrors the per-draw construction in :mod:`opencea.builders` but writes
    every entry into a 3-D tensor so that the engine's per-cycle matmul can
    be batched over draws.
    """
    n_sim = p_HD.size
    P = np.zeros((n_sim, 4, 4), dtype=float)
    # From H
    P[:, 0, 0] = (1.0 - p_HD) * (1.0 - p_HS1)
    P[:, 0, 1] = (1.0 - p_HD) * p_HS1
    P[:, 0, 3] = p_HD
    # From S1
    P[:, 1, 0] = (1.0 - p_S1D) * p_S1H
    P[:, 1, 1] = (1.0 - p_S1D) * (1.0 - (p_S1H + p_S1S2))
    P[:, 1, 2] = (1.0 - p_S1D) * p_S1S2
    P[:, 1, 3] = p_S1D
    # From S2
    P[:, 2, 2] = 1.0 - p_S2D
    P[:, 2, 3] = p_S2D
    # From D (absorbing)
    P[:, 3, 3] = 1.0
    return P


def _build_batched_strategies(
    draws: pd.DataFrame, cycle_length: float
) -> Dict[str, Dict[str, np.ndarray]]:
    """Assemble per-strategy ``(P, costs, utils)`` tensors for all draws."""
    r_HD = draws["r_HD"].to_numpy()
    r_HS1 = draws["r_HS1"].to_numpy()
    r_S1H = draws["r_S1H"].to_numpy()
    r_S1S2 = draws["r_S1S2"].to_numpy()
    hr_S1 = draws["hr_S1"].to_numpy()
    hr_S2 = draws["hr_S2"].to_numpy()
    hr_S1S2_trtB = draws["hr_S1S2_trtB"].to_numpy()

    # Mortality rates in Sick / Sicker = baseline rate * hazard ratio.
    r_S1D = r_HD * hr_S1
    r_S2D = r_HD * hr_S2

    p_HS1 = _rate_to_prob_vec(r_HS1, cycle_length)
    p_S1H = _rate_to_prob_vec(r_S1H, cycle_length)
    p_S1S2 = _rate_to_prob_vec(r_S1S2, cycle_length)
    p_HD = _rate_to_prob_vec(r_HD, cycle_length)
    p_S1D = _rate_to_prob_vec(r_S1D, cycle_length)
    p_S2D = _rate_to_prob_vec(r_S2D, cycle_length)

    # Treatment B reduces the *rate* of S1 -> S2 before the rate -> prob conversion.
    r_S1S2_trtB = r_S1S2 * hr_S1S2_trtB
    p_S1S2_trtB = _rate_to_prob_vec(r_S1S2_trtB, cycle_length)

    P_SoC = _build_transition_tensor(p_HD, p_HS1, p_S1H, p_S1S2, p_S1D, p_S2D)
    P_A = P_SoC  # Strategy A leaves transitions unchanged
    P_B = _build_transition_tensor(p_HD, p_HS1, p_S1H, p_S1S2_trtB, p_S1D, p_S2D)
    P_AB = P_B

    c_H = draws["c_H"].to_numpy()
    c_S1 = draws["c_S1"].to_numpy()
    c_S2 = draws["c_S2"].to_numpy()
    c_D = draws["c_D"].to_numpy()
    c_trtA = draws["c_trtA"].to_numpy()
    c_trtB = draws["c_trtB"].to_numpy()

    u_H = draws["u_H"].to_numpy()
    u_S1 = draws["u_S1"].to_numpy()
    u_S2 = draws["u_S2"].to_numpy()
    u_D = draws["u_D"].to_numpy()
    u_trtA = draws["u_trtA"].to_numpy()

    def stack(*cols: np.ndarray) -> np.ndarray:
        return np.column_stack(cols)

    return {
        "Standard of care": {
            "P": P_SoC,
            "costs": stack(c_H, c_S1, c_S2, c_D),
            "utils": stack(u_H, u_S1, u_S2, u_D),
        },
        "Strategy A": {
            "P": P_A,
            "costs": stack(c_H, c_S1 + c_trtA, c_S2 + c_trtA, c_D),
            "utils": stack(u_H, u_trtA, u_S2, u_D),
        },
        "Strategy B": {
            "P": P_B,
            "costs": stack(c_H, c_S1 + c_trtB, c_S2 + c_trtB, c_D),
            "utils": stack(u_H, u_S1, u_S2, u_D),
        },
        "Strategy AB": {
            "P": P_AB,
            "costs": stack(c_H, c_S1 + c_trtA + c_trtB, c_S2 + c_trtA + c_trtB, c_D),
            "utils": stack(u_H, u_trtA, u_S2, u_D),
        },
    }


# ---------------------------------------------------------------------------
# Batched simulation
# ---------------------------------------------------------------------------


def _discount_weights(rate: float, n_cycles: int, cycle_length: float) -> np.ndarray:
    """DARTH-style discount weights ``1 / (1 + d * cycle_length) ** t``.

    Identical to ``opencea.engine._discount_weights``; reproduced here so
    the PSA module stays self-contained for the formula it accumulates
    against the batched trace.
    """
    t = np.arange(n_cycles + 1, dtype=float)
    return 1.0 / (1.0 + rate * cycle_length) ** t


def _simulate_strategy_batched(
    P: np.ndarray,
    costs: np.ndarray,
    utils: np.ndarray,
    init: np.ndarray,
    n_cycles: int,
    cycle_length: float,
    d_c: float,
    d_e: float,
    wcc_method: WCCMethod,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run one strategy across all draws; return ``(total_cost, total_qaly)``.

    The per-cycle step ``trace[t+1] = trace[t] @ P`` becomes a single
    ``einsum('is,isk->ik', trace[t], P)`` over the ``(n_sim, n_states,
    n_states)`` transition tensor.
    """
    n_sim, n_states, _ = P.shape
    trace = np.empty((n_cycles + 1, n_sim, n_states), dtype=float)
    trace[0] = init  # broadcasts (n_states,) -> (n_sim, n_states)
    for t in range(n_cycles):
        trace[t + 1] = np.einsum("is,isk->ik", trace[t], P)

    # Per-cycle reward over the cohort, scaled by cycle_length (DARTH convention).
    cycle_costs = np.einsum("tis,is->ti", trace, costs) * cycle_length
    cycle_qalys = np.einsum("tis,is->ti", trace, utils) * cycle_length

    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    dw_e = _discount_weights(d_e, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=wcc_method)

    total_cost = (cycle_costs * (dw_c * wcc)[:, None]).sum(axis=0)
    total_qaly = (cycle_qalys * (dw_e * wcc)[:, None]).sum(axis=0)
    return total_cost, total_qaly


# ---------------------------------------------------------------------------
# Public PSA result and driver
# ---------------------------------------------------------------------------


@dataclass
class PSAResult:
    """Tidy PSA output.

    Attributes
    ----------
    params
        ``(n_sim, n_params)`` DataFrame of the sampled draws.
    costs, qalys
        ``(n_sim, n_strategies)`` DataFrames keyed by strategy name.
    strategy_names
        Order of strategies as columns in ``costs`` / ``qalys``.
    n_sim
        Number of draws.
    """

    params: pd.DataFrame
    costs: pd.DataFrame
    qalys: pd.DataFrame
    strategy_names: Tuple[str, ...]
    n_sim: int


def _load_base_params(
    base_params: ParameterInput,
) -> Mapping[str, Any]:
    if isinstance(base_params, (str, Path)):
        import yaml

        with open(base_params, "r") as f:
            return yaml.safe_load(f)
    return base_params


def run_psa(
    base_params: ParameterInput,
    n_sim: int = 10_000,
    seed: int = 20260625,
    wcc_method: WCCMethod = "simpson_1_3",
) -> PSAResult:
    """Run the full DARTH Sick-Sicker PSA in vectorized form.

    Parameters
    ----------
    base_params
        Either a mapping or a path to a YAML file holding the model's
        non-sampled fixed parameters (``cycle_length``, ``n_age_init``,
        ``n_age_max``, ``d_c``, ``d_e``). The rate / cost / utility
        parameters in the file are ignored — they are replaced by the
        PSA draws.
    n_sim
        Number of parameter sets to draw (default 10000).
    seed
        Seed for :class:`numpy.random.default_rng`.
    wcc_method
        Within-cycle-correction method passed to
        :func:`opencea.engine.gen_wcc`. Defaults to Simpson's 1/3 to
        match the deterministic engine.
    """
    bp = _load_base_params(base_params)

    cycle_length = float(bp["cycle_length"])
    n_cycles = int((float(bp["n_age_max"]) - float(bp["n_age_init"])) / cycle_length)
    if n_cycles <= 0:
        raise ValueError("n_age_max must exceed n_age_init")
    d_c = float(bp["d_c"])
    d_e = float(bp["d_e"])

    draws = sample_psa_params(n_sim=n_sim, seed=seed)
    batched = _build_batched_strategies(draws, cycle_length=cycle_length)

    init = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    cost_cols: Dict[str, np.ndarray] = {}
    qaly_cols: Dict[str, np.ndarray] = {}
    for name in _STRATEGY_NAMES:
        bundle = batched[name]
        total_cost, total_qaly = _simulate_strategy_batched(
            P=bundle["P"],
            costs=bundle["costs"],
            utils=bundle["utils"],
            init=init,
            n_cycles=n_cycles,
            cycle_length=cycle_length,
            d_c=d_c,
            d_e=d_e,
            wcc_method=wcc_method,
        )
        cost_cols[name] = total_cost
        qaly_cols[name] = total_qaly

    return PSAResult(
        params=draws,
        costs=pd.DataFrame(cost_cols),
        qalys=pd.DataFrame(qaly_cols),
        strategy_names=_STRATEGY_NAMES,
        n_sim=n_sim,
    )


# ---------------------------------------------------------------------------
# PSA cost-effectiveness outputs
# ---------------------------------------------------------------------------


def default_wtp_grid(
    low: float = 0.0,
    high: float = 200_000.0,
    step: float = 1_000.0,
) -> np.ndarray:
    """Default willingness-to-pay grid: 0 to 200000 in 1000-unit steps."""
    n = int(round((high - low) / step)) + 1
    return np.linspace(low, high, n)


def compute_nmb(result: PSAResult, wtp_grid: np.ndarray) -> np.ndarray:
    """Net monetary benefit per draw per strategy per WTP.

    Returns an array of shape ``(n_sim, n_strategies, n_wtp)`` ordered to
    match ``result.strategy_names`` and ``wtp_grid``.
    """
    wtp = np.asarray(wtp_grid, dtype=float)
    costs = result.costs[list(result.strategy_names)].to_numpy()  # (n_sim, n_strat)
    qalys = result.qalys[list(result.strategy_names)].to_numpy()
    # NMB[i, s, w] = wtp[w] * qalys[i, s] - costs[i, s]
    return qalys[:, :, None] * wtp[None, None, :] - costs[:, :, None]


def compute_ceac(result: PSAResult, wtp_grid: np.ndarray) -> pd.DataFrame:
    """Cost-effectiveness acceptability curve.

    Returns a DataFrame indexed by WTP with one column per strategy whose
    value is the share of draws in which that strategy attains the maximum
    NMB at the given WTP. Ties are broken in favour of the lowest-indexed
    strategy (negligible at typical scales).
    """
    wtp = np.asarray(wtp_grid, dtype=float)
    nmb = compute_nmb(result, wtp)  # (n_sim, n_strat, n_wtp)
    winners = nmb.argmax(axis=1)  # (n_sim, n_wtp)
    n_strat = len(result.strategy_names)
    n_wtp = wtp.size
    counts = np.zeros((n_wtp, n_strat), dtype=float)
    for s in range(n_strat):
        counts[:, s] = (winners == s).sum(axis=0)
    probs = counts / result.n_sim
    return pd.DataFrame(
        probs, index=pd.Index(wtp, name="wtp"), columns=list(result.strategy_names)
    )


def expected_nmb_frontier(result: PSAResult, wtp_grid: np.ndarray) -> pd.DataFrame:
    """Expected-NMB-optimal strategy at each WTP (cost-effectiveness frontier).

    Returns a DataFrame indexed by WTP with columns ``best_strategy`` (name
    of the strategy maximizing expected NMB) and ``expected_nmb`` (its
    expected NMB at that WTP). Also includes a column per strategy with its
    expected NMB.
    """
    wtp = np.asarray(wtp_grid, dtype=float)
    nmb = compute_nmb(result, wtp)  # (n_sim, n_strat, n_wtp)
    expected = nmb.mean(axis=0)  # (n_strat, n_wtp)
    best_idx = expected.argmax(axis=0)  # (n_wtp,)
    names = list(result.strategy_names)
    best_name = np.array(names)[best_idx]
    out = pd.DataFrame(expected.T, index=pd.Index(wtp, name="wtp"), columns=names)
    out["best_strategy"] = best_name
    out["expected_nmb"] = expected[best_idx, np.arange(wtp.size)]
    return out


def incremental_vs_baseline(
    result: PSAResult, baseline: str = "Standard of care"
) -> pd.DataFrame:
    """Per-draw incremental cost and QALY of each non-baseline strategy.

    Returns a long-form DataFrame with columns ``draw``, ``strategy``,
    ``inc_cost``, ``inc_qaly`` — the raw material for a cost-effectiveness
    plane scatter.
    """
    if baseline not in result.strategy_names:
        raise KeyError(
            f"baseline {baseline!r} not in strategies {result.strategy_names}"
        )
    base_cost = result.costs[baseline].to_numpy()
    base_qaly = result.qalys[baseline].to_numpy()
    rows: List[pd.DataFrame] = []
    n_sim = result.n_sim
    for name in result.strategy_names:
        if name == baseline:
            continue
        rows.append(
            pd.DataFrame(
                {
                    "draw": np.arange(n_sim),
                    "strategy": name,
                    "inc_cost": result.costs[name].to_numpy() - base_cost,
                    "inc_qaly": result.qalys[name].to_numpy() - base_qaly,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)
