"""Applied case study: Empagliflozin vs SoC in T2D + established CVD.

Builds a 3-state cohort model (EF = event-free T2D+CVD, PE = post major CV
event, D = dead) with two strategies, reusing the validated engine,
discount weights, within-cycle correction, PSA sampler, and DSA driver.
Nothing in the simulation logic is reimplemented here.

The one-time acute-event cost on entering PE is implemented as a
**transition cost**: the trace is unchanged, and the engine's per-cycle
discount + WCC weights are reused to discount the cost at the cycle when
new PE entrants appear (cycle ``t + 1`` for transitions that occur during
cycle ``t``).

References for every parameter are inline in
``examples/empagliflozin_t2d.yaml``; the CHEERS-structured writeup is in
``examples/empagliflozin_case_study.md``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Tuple, Union

import numpy as np
import pandas as pd

from .builders import rate_to_prob
from .engine import evaluate_strategy, gen_wcc, run_model
from .model import CohortModel, Strategy
from .psa import (
    DistSpec,
    PSAResult,
    _discount_weights,  # private helpers reused; the alternative would be
    _simulate_strategy_batched,  # to copy the formulas, which the user
    sample_psa_params,             # explicitly forbade
)


# Strategy names — used by builders, DSA, PSA, and the test suite.
SOC = "Standard of care"
EMPA = "Empagliflozin"
STRATEGY_NAMES = (SOC, EMPA)
STATES = ("EF", "PE", "D")
EF_IDX, PE_IDX, D_IDX = 0, 1, 2

EMPA_REQUIRED_KEYS = {
    "cycle_length", "n_age_init", "n_age_max",
    "d_c", "d_e",
    "r_EF_PE", "r_EF_D", "r_PE_D",
    "hr_event", "hr_death",
    "c_drug", "c_EF", "c_PE", "c_acute_PE", "c_D",
    "u_EF", "u_PE", "u_D",
}


def _check_params(params: Mapping[str, Any]) -> None:
    missing = EMPA_REQUIRED_KEYS - set(params)
    if missing:
        raise KeyError(
            "Empagliflozin parameter set is missing required keys: "
            f"{sorted(missing)}"
        )


def _build_transition_matrix(
    p_EF_PE_event: float, p_EF_D: float, p_PE_D: float
) -> list[list[float]]:
    """3x3 competing-risks transition matrix (EF, PE, D).

    Same construction as ``opencea.builders._build_transition_matrix`` —
    non-death transitions are scaled by ``(1 - p_kD)`` so that the row
    sums to 1 by construction.
    """
    return [
        [
            (1.0 - p_EF_D) * (1.0 - p_EF_PE_event),  # EF -> EF
            (1.0 - p_EF_D) * p_EF_PE_event,          # EF -> PE
            p_EF_D,                                  # EF -> D
        ],
        [
            0.0,                                     # PE -> EF (no recovery)
            1.0 - p_PE_D,                            # PE -> PE
            p_PE_D,                                  # PE -> D
        ],
        [0.0, 0.0, 1.0],                             # D absorbing
    ]


def _load_params(p: Union[Mapping[str, Any], str, Path]) -> Dict[str, Any]:
    if isinstance(p, (str, Path)):
        import yaml

        with open(p, "r") as f:
            return dict(yaml.safe_load(f))
    return dict(p)


def build_empagliflozin_t2d(
    params: Union[Mapping[str, Any], str, Path],
) -> Tuple[CohortModel, float]:
    """Construct the empagliflozin case-study ``CohortModel``.

    Returns a tuple ``(model, c_acute_PE)``. The acute-event cost is
    returned alongside the model because it is a **transition** cost —
    not a per-cycle state cost — and the ``CohortModel`` schema only
    holds per-cycle state costs. The caller (typically
    :func:`evaluate_empagliflozin_case`) folds the acute cost in
    post-hoc using the same discount / WCC weights the engine uses.
    """
    p = _load_params(params)
    _check_params(p)

    cycle_length = float(p["cycle_length"])

    # ---- transition probabilities ---------------------------------------
    r_EF_PE = float(p["r_EF_PE"])
    r_EF_D = float(p["r_EF_D"])
    r_PE_D = float(p["r_PE_D"])
    hr_event = float(p["hr_event"])
    hr_death = float(p["hr_death"])

    # SoC
    p_EF_PE_soc = rate_to_prob(r_EF_PE, cycle_length)
    p_EF_D_soc = rate_to_prob(r_EF_D, cycle_length)
    p_PE_D_soc = rate_to_prob(r_PE_D, cycle_length)
    P_SoC = _build_transition_matrix(p_EF_PE_soc, p_EF_D_soc, p_PE_D_soc)

    # Empagliflozin: HRs scale the rates BEFORE the rate->prob conversion.
    p_EF_PE_empa = rate_to_prob(r_EF_PE * hr_event, cycle_length)
    p_EF_D_empa = rate_to_prob(r_EF_D * hr_death, cycle_length)
    p_PE_D_empa = rate_to_prob(r_PE_D * hr_death, cycle_length)
    P_Empa = _build_transition_matrix(p_EF_PE_empa, p_EF_D_empa, p_PE_D_empa)

    # ---- state costs and utilities --------------------------------------
    c_drug = float(p["c_drug"])
    c_EF = float(p["c_EF"])
    c_PE = float(p["c_PE"])
    c_D = float(p["c_D"])
    c_acute_PE = float(p["c_acute_PE"])

    u_EF = float(p["u_EF"])
    u_PE = float(p["u_PE"])
    u_D = float(p["u_D"])

    costs_SoC = [c_EF, c_PE, c_D]
    costs_Empa = [c_EF + c_drug, c_PE + c_drug, c_D]
    util_SoC = [u_EF, u_PE, u_D]
    util_Empa = [u_EF, u_PE, u_D]

    n_cycles = int(
        (float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length
    )
    if n_cycles <= 0:
        raise ValueError("n_age_max must exceed n_age_init")

    model = CohortModel(
        states=list(STATES),
        initial_distribution=[1.0, 0.0, 0.0],
        cycle_length=cycle_length,
        time_horizon=n_cycles,
        discount_rate_costs=float(p["d_c"]),
        discount_rate_qalys=float(p["d_e"]),
        wcc_method="simpson_1_3",
        strategies=[
            Strategy(
                name=SOC,
                transition_matrix=P_SoC,
                state_costs=costs_SoC,
                state_utilities=util_SoC,
            ),
            Strategy(
                name=EMPA,
                transition_matrix=P_Empa,
                state_costs=costs_Empa,
                state_utilities=util_Empa,
            ),
        ],
    )
    return model, c_acute_PE


# ---------------------------------------------------------------------------
# Deterministic evaluator (engine + acute-cost transition contribution)
# ---------------------------------------------------------------------------


def _acute_cost_contribution(
    trace: np.ndarray,
    P: np.ndarray,
    c_acute_PE: float,
    dw_c: np.ndarray,
    wcc: np.ndarray,
) -> float:
    """Discounted total of one-time acute-event costs on EF -> PE transitions.

    New PE entrants during cycle ``t`` (the flow from ``t`` to ``t + 1``)
    equal ``trace[t, EF] * P[EF, PE]``. They appear in PE at cycle ``t + 1``
    and the acute cost is charged then, weighted by ``dw_c[t+1] * wcc[t+1]``
    so the discount / WCC treatment matches the engine's per-cycle
    convention exactly.
    """
    T = trace.shape[0] - 1
    entries = trace[:T, EF_IDX] * P[EF_IDX, PE_IDX]
    weights = dw_c[1:T + 1] * wcc[1:T + 1]
    return float(np.sum(entries * weights) * c_acute_PE)


def evaluate_empagliflozin_case(
    params: Union[Mapping[str, Any], str, Path],
) -> Dict[str, Dict[str, float]]:
    """Run the case study deterministically; return per-strategy totals.

    For each strategy the returned dict carries:

      - ``total_cost`` — engine state-cost total **plus** the discounted
        acute-event transition cost.
      - ``total_qaly`` — engine QALY total (acute event has no utility
        decrement of its own here).
      - ``state_cost`` — the engine's state-cost-only total, for audit.
      - ``acute_event_cost`` — the discounted transition-cost contribution.
      - ``trace`` — the cohort trace (T+1 x 3).
    """
    p = _load_params(params)
    model, c_acute_PE = build_empagliflozin_t2d(p)

    dw_c = _discount_weights(
        model.discount_rate_costs, model.time_horizon, model.cycle_length
    )
    wcc = gen_wcc(model.time_horizon, method=model.wcc_method)

    out: Dict[str, Dict[str, float]] = {}
    for strat in model.strategies:
        r = evaluate_strategy(strat, model)
        P = np.asarray(strat.transition_matrix, dtype=float)
        acute = _acute_cost_contribution(r["trace"], P, c_acute_PE, dw_c, wcc)
        out[r["name"]] = {
            "total_cost": float(r["total_cost"]) + acute,
            "total_qaly": float(r["total_qaly"]),
            "state_cost": float(r["total_cost"]),
            "acute_event_cost": acute,
            "trace": r["trace"],
        }
    return out


def case_results_for_cea(
    params: Union[Mapping[str, Any], str, Path],
) -> list[Dict[str, object]]:
    """Adapter from :func:`evaluate_empagliflozin_case` to the shape
    consumed by :func:`opencea.cea.cea_table`.
    """
    out = evaluate_empagliflozin_case(params)
    return [
        {"name": k, "total_cost": v["total_cost"], "total_qaly": v["total_qaly"]}
        for k, v in out.items()
    ]


# ---------------------------------------------------------------------------
# DSA / PSA glue
# ---------------------------------------------------------------------------


def dsa_evaluator(params: Mapping[str, Any]) -> Dict[str, Dict[str, float]]:
    """Evaluator compatible with :func:`opencea.sensitivity.run_dsa`."""
    out = evaluate_empagliflozin_case(params)
    return {
        name: {"cost": v["total_cost"], "qaly": v["total_qaly"]}
        for name, v in out.items()
    }


# Distributions for PSA. HRs are lognormal, parameterized so that the
# distribution's median equals the trial point estimate and the 95% range
# matches the published 95% CI:
#     median = exp(meanlog),  CI = median * exp(+/- 1.96 * sdlog)
#  => sdlog = (ln(UCL) - ln(LCL)) / (2 * 1.96)
def _ln_sigma_from_ci(ucl: float, lcl: float) -> float:
    return (np.log(ucl) - np.log(lcl)) / (2.0 * 1.96)


_HR_EVENT_SIGMA = _ln_sigma_from_ci(0.99, 0.74)  # ~ 0.0743
_HR_DEATH_SIGMA = _ln_sigma_from_ci(0.82, 0.57)  # ~ 0.0928
_R_PE_D_SIGMA = _ln_sigma_from_ci(0.10, 0.04)    # ~ 0.234, the wide soft band


# Gamma with mean m and CV 0.25 => shape = 1/CV^2 = 16, scale = m/16.
def _gamma_mean_cv(mean: float, cv: float = 0.25) -> Dict[str, float]:
    shape = 1.0 / (cv ** 2)
    return {"shape": shape, "scale": mean / shape}


EMPA_PSA_SPECS: Dict[str, DistSpec] = {
    # Treatment effect — trial-derived lognormals.
    "hr_event": DistSpec(
        "lognormal",
        {"mean": float(np.log(0.86)), "sigma": _HR_EVENT_SIGMA},
        target_mean=0.86,
        is_median=True,
    ),
    "hr_death": DistSpec(
        "lognormal",
        {"mean": float(np.log(0.68)), "sigma": _HR_DEATH_SIGMA},
        target_mean=0.68,
        is_median=True,
    ),
    # Costs — gamma with 25% CV (the spec the user asked for).
    "c_drug":     DistSpec("gamma", _gamma_mean_cv(6264.0),  target_mean=6264.0),
    "c_EF":       DistSpec("gamma", _gamma_mean_cv(16000.0), target_mean=16000.0),
    "c_PE":       DistSpec("gamma", _gamma_mean_cv(20000.0), target_mean=20000.0),
    "c_acute_PE": DistSpec("gamma", _gamma_mean_cv(11650.0), target_mean=11650.0),
    # Utilities — beta with a + b = 100 (modest precision around the mean).
    "u_EF": DistSpec("beta", {"a": 75.0, "b": 25.0}, target_mean=0.75),
    "u_PE": DistSpec("beta", {"a": 68.0, "b": 32.0}, target_mean=0.68),
    # Soft post-event mortality — wide lognormal spanning 0.04 - 0.10 / yr.
    "r_PE_D": DistSpec(
        "lognormal",
        {"mean": float(np.log(0.06)), "sigma": _R_PE_D_SIGMA},
        target_mean=0.06,
        is_median=True,
    ),
}


# DSA range overrides. Where the user spec gives explicit low/high values
# we use them verbatim; the remaining parameters (PE cost, acute PE cost,
# utilities) fall back to the PSA percentile machinery in run_dsa.
EMPA_DSA_RANGE_OVERRIDES: Dict[str, Tuple[float, float]] = {
    "hr_event": (0.74, 0.99),   # trial 95% CI on the composite event
    "hr_death": (0.57, 0.82),   # trial 95% CI on all-cause mortality
    "c_drug":   (4872.0, 7596.0),
    "c_EF":     (13000.0, 19700.0),
    "r_PE_D":   (0.04, 0.10),   # wide soft band (user spec)
}


# ---------------------------------------------------------------------------
# Vectorized PSA — produces a generic PSAResult
# ---------------------------------------------------------------------------


def _empa_transition_tensor(
    p_EF_PE_event: np.ndarray, p_EF_D: np.ndarray, p_PE_D: np.ndarray
) -> np.ndarray:
    """Batched 3x3 competing-risks transition matrices."""
    n_sim = p_EF_PE_event.size
    P = np.zeros((n_sim, 3, 3), dtype=float)
    P[:, 0, 0] = (1.0 - p_EF_D) * (1.0 - p_EF_PE_event)
    P[:, 0, 1] = (1.0 - p_EF_D) * p_EF_PE_event
    P[:, 0, 2] = p_EF_D
    P[:, 1, 1] = 1.0 - p_PE_D
    P[:, 1, 2] = p_PE_D
    P[:, 2, 2] = 1.0
    return P


def _batched_acute_cost(
    P: np.ndarray,
    c_acute_PE: np.ndarray,
    init: np.ndarray,
    n_cycles: int,
    cycle_length: float,
    d_c: float,
    wcc_method: str,
) -> np.ndarray:
    """Discounted acute-event cost per draw — vectorized over ``n_sim``.

    Reuses the engine's discount + WCC weights (no reimplementation), and
    propagates the trace via the same batched einsum step that the PSA
    module uses for the DARTH model.
    """
    n_sim = P.shape[0]
    trace = np.empty((n_cycles + 1, n_sim, 3), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = np.einsum("is,isk->ik", trace[t], P)
    # Flow from EF into PE per cycle: (n_cycles, n_sim)
    entries = trace[:n_cycles, :, EF_IDX] * P[:, EF_IDX, PE_IDX][None, :]
    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=wcc_method)
    weights = dw_c[1:n_cycles + 1] * wcc[1:n_cycles + 1]
    return (entries * weights[:, None]).sum(axis=0) * c_acute_PE


def run_empa_psa(
    base_params: Union[Mapping[str, Any], str, Path],
    n_sim: int = 10_000,
    seed: int = 20260626,
    wcc_method: str = "simpson_1_3",
) -> PSAResult:
    """Vectorized PSA for the empagliflozin case study.

    The transition tensors, the cycle-trace propagation, and the rewards
    accumulation are all done with the same helpers used by the DARTH
    PSA (``_simulate_strategy_batched``, ``_discount_weights``,
    :func:`opencea.engine.gen_wcc`). The only case-study-specific piece
    is the acute-event transition-cost contribution, which is added per
    draw on top of the state-cost total.
    """
    p = _load_params(base_params)
    _check_params(p)

    cycle_length = float(p["cycle_length"])
    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
    d_c = float(p["d_c"])
    d_e = float(p["d_e"])

    # Baseline event rates are held fixed at their point estimates; trial
    # uncertainty is carried through the HRs and through r_PE_D (the
    # explicitly soft input).
    r_EF_PE = float(p["r_EF_PE"])
    r_EF_D = float(p["r_EF_D"])

    draws = sample_psa_params(n_sim=n_sim, seed=seed, specs=EMPA_PSA_SPECS)
    hr_event = draws["hr_event"].to_numpy()
    hr_death = draws["hr_death"].to_numpy()
    r_PE_D = draws["r_PE_D"].to_numpy()
    c_drug = draws["c_drug"].to_numpy()
    c_EF = draws["c_EF"].to_numpy()
    c_PE = draws["c_PE"].to_numpy()
    c_acute_PE = draws["c_acute_PE"].to_numpy()
    u_EF = draws["u_EF"].to_numpy()
    u_PE = draws["u_PE"].to_numpy()

    zero = np.zeros(n_sim)

    # --- SoC ---
    p_EF_PE_soc = np.full(n_sim, 1.0 - np.exp(-r_EF_PE * cycle_length))
    p_EF_D_soc = np.full(n_sim, 1.0 - np.exp(-r_EF_D * cycle_length))
    p_PE_D_soc = 1.0 - np.exp(-r_PE_D * cycle_length)
    P_SoC = _empa_transition_tensor(p_EF_PE_soc, p_EF_D_soc, p_PE_D_soc)
    costs_SoC = np.column_stack([c_EF, c_PE, zero])
    util_SoC = np.column_stack([u_EF, u_PE, zero])

    # --- Empagliflozin ---
    p_EF_PE_empa = 1.0 - np.exp(-r_EF_PE * hr_event * cycle_length)
    p_EF_D_empa = 1.0 - np.exp(-r_EF_D * hr_death * cycle_length)
    p_PE_D_empa = 1.0 - np.exp(-r_PE_D * hr_death * cycle_length)
    P_Empa = _empa_transition_tensor(p_EF_PE_empa, p_EF_D_empa, p_PE_D_empa)
    costs_Empa = np.column_stack([c_EF + c_drug, c_PE + c_drug, zero])
    util_Empa = np.column_stack([u_EF, u_PE, zero])

    init = np.array([1.0, 0.0, 0.0], dtype=float)

    cost_soc_state, qaly_soc = _simulate_strategy_batched(
        P_SoC, costs_SoC, util_SoC, init,
        n_cycles, cycle_length, d_c, d_e, wcc_method,
    )
    cost_empa_state, qaly_empa = _simulate_strategy_batched(
        P_Empa, costs_Empa, util_Empa, init,
        n_cycles, cycle_length, d_c, d_e, wcc_method,
    )

    acute_soc = _batched_acute_cost(
        P_SoC, c_acute_PE, init, n_cycles, cycle_length, d_c, wcc_method
    )
    acute_empa = _batched_acute_cost(
        P_Empa, c_acute_PE, init, n_cycles, cycle_length, d_c, wcc_method
    )

    cost_soc = cost_soc_state + acute_soc
    cost_empa = cost_empa_state + acute_empa

    return PSAResult(
        params=draws,
        costs=pd.DataFrame({SOC: cost_soc, EMPA: cost_empa}),
        qalys=pd.DataFrame({SOC: qaly_soc, EMPA: qaly_empa}),
        strategy_names=STRATEGY_NAMES,
        n_sim=n_sim,
    )
