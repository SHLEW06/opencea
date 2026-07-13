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

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ._types import CaseStrategyResult, CEAStrategyResult, DSAOutcome, ParameterInput
from .builders import rate_to_prob
from .engine import (
    WCCMethod,
    evaluate_sequence,
    evaluate_strategy,
    gen_wcc,
)
from .model import CohortModel, Strategy
from .psa import (
    DistSpec,
    PSAResult,
    _discount_weights,  # private helpers reused; the alternative would be
    _simulate_strategy_batched,  # to copy the formulas, which the user
    sample_psa_params,  # explicitly forbade
)

# Strategy names — used by builders, DSA, PSA, and the test suite.
SOC = "Standard of care"
EMPA = "Empagliflozin"
STRATEGY_NAMES = (SOC, EMPA)
STATES = ("EF", "PE", "D")
EF_IDX, PE_IDX, D_IDX = 0, 1, 2

EMPA_REQUIRED_KEYS = {
    "cycle_length",
    "n_age_init",
    "n_age_max",
    "d_c",
    "d_e",
    "r_EF_PE",
    "r_EF_D",
    "r_PE_D",
    "hr_event",
    "hr_death",
    "c_drug",
    "c_EF",
    "c_PE",
    "c_acute_PE",
    "c_D",
    "u_EF",
    "u_PE",
    "u_D",
}


def _check_params(params: Mapping[str, Any]) -> None:
    missing = EMPA_REQUIRED_KEYS - set(params)
    if missing:
        raise KeyError(
            f"Empagliflozin parameter set is missing required keys: {sorted(missing)}"
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
            (1.0 - p_EF_D) * p_EF_PE_event,  # EF -> PE
            p_EF_D,  # EF -> D
        ],
        [
            0.0,  # PE -> EF (no recovery)
            1.0 - p_PE_D,  # PE -> PE
            p_PE_D,  # PE -> D
        ],
        [0.0, 0.0, 1.0],  # D absorbing
    ]


def _load_params(p: ParameterInput) -> Dict[str, Any]:
    if isinstance(p, (str, Path)):
        import yaml

        with open(p, "r") as f:
            return dict(yaml.safe_load(f))
    return dict(p)


def build_empagliflozin_t2d(params: ParameterInput) -> Tuple[CohortModel, float]:
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

    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
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
    weights = dw_c[1 : T + 1] * wcc[1 : T + 1]
    return float(np.sum(entries * weights) * c_acute_PE)


def evaluate_empagliflozin_case(
    params: ParameterInput,
) -> dict[str, CaseStrategyResult]:
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

    out: dict[str, CaseStrategyResult] = {}
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


def case_results_for_cea(params: ParameterInput) -> list[CEAStrategyResult]:
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


def dsa_evaluator(params: Mapping[str, Any]) -> dict[str, DSAOutcome]:
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
_R_PE_D_SIGMA = _ln_sigma_from_ci(0.10, 0.04)  # ~ 0.234, the wide soft band


# Gamma with mean m and CV 0.25 => shape = 1/CV^2 = 16, scale = m/16.
def _gamma_mean_cv(mean: float, cv: float = 0.25) -> Dict[str, float]:
    shape = 1.0 / (cv**2)
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
    "c_drug": DistSpec("gamma", _gamma_mean_cv(6264.0), target_mean=6264.0),
    "c_EF": DistSpec("gamma", _gamma_mean_cv(16000.0), target_mean=16000.0),
    "c_PE": DistSpec("gamma", _gamma_mean_cv(20000.0), target_mean=20000.0),
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
    "hr_event": (0.74, 0.99),  # trial 95% CI on the composite event
    "hr_death": (0.57, 0.82),  # trial 95% CI on all-cause mortality
    "c_drug": (4872.0, 7596.0),
    "c_EF": (13000.0, 19700.0),
    "r_PE_D": (0.04, 0.10),  # wide soft band (user spec)
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
    wcc_method: WCCMethod,
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
    weights = dw_c[1 : n_cycles + 1] * wcc[1 : n_cycles + 1]
    return (entries * weights[:, None]).sum(axis=0) * c_acute_PE


def run_empa_psa(
    base_params: ParameterInput,
    n_sim: int = 10_000,
    seed: int = 20260626,
    wcc_method: WCCMethod = "simpson_1_3",
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
        P_SoC,
        costs_SoC,
        util_SoC,
        init,
        n_cycles,
        cycle_length,
        d_c,
        d_e,
        wcc_method,
    )
    cost_empa_state, qaly_empa = _simulate_strategy_batched(
        P_Empa,
        costs_Empa,
        util_Empa,
        init,
        n_cycles,
        cycle_length,
        d_c,
        d_e,
        wcc_method,
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


# ---------------------------------------------------------------------------
# Scenario analysis: net-price and treatment-effect waning
# ---------------------------------------------------------------------------


def _run_empa_psa_fixed_drug_price(
    p: Mapping[str, Any],
    n_sim: int,
    seed: int,
    wcc_method: WCCMethod,
    fixed_drug_price: float,
) -> PSAResult:
    """Sustained-effect PSA with ``c_drug`` pinned to ``fixed_drug_price``.

    Identical to :func:`run_empa_psa` except the sampled ``c_drug``
    column is replaced by a constant. Other parameters keep their PSA
    distributions, so cost-effectiveness uncertainty still propagates
    through HRs / state costs / utilities / r_PE_D.
    """
    result = run_empa_psa(p, n_sim=n_sim, seed=seed, wcc_method=wcc_method)
    # Recompute Empa costs with the fixed drug price. Empa cost differs
    # from SoC cost only through (i) the extra annual drug cost in EF
    # and PE, and (ii) the different transition probabilities (hr_event,
    # hr_death). The transition piece is already baked into result;
    # we only need to swap in the fixed drug price by computing the
    # drug-cost differential explicitly.
    cycle_length = float(p["cycle_length"])
    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
    d_c = float(p["d_c"])

    draws = result.params
    # Restore the trace under Empa exactly as run_empa_psa built it so
    # we can recompute drug spend as (occupied person-time in EF + PE)
    # * fixed_drug_price * discount weights.
    r_EF_PE = float(p["r_EF_PE"])
    r_EF_D = float(p["r_EF_D"])
    hr_event = draws["hr_event"].to_numpy()
    hr_death = draws["hr_death"].to_numpy()
    r_PE_D = draws["r_PE_D"].to_numpy()
    p_EF_PE_empa = 1.0 - np.exp(-r_EF_PE * hr_event * cycle_length)
    p_EF_D_empa = 1.0 - np.exp(-r_EF_D * hr_death * cycle_length)
    p_PE_D_empa = 1.0 - np.exp(-r_PE_D * hr_death * cycle_length)
    P_empa = _empa_transition_tensor(p_EF_PE_empa, p_EF_D_empa, p_PE_D_empa)
    init = np.array([1.0, 0.0, 0.0])
    n_sim = result.n_sim
    trace = np.empty((n_cycles + 1, n_sim, 3), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = np.einsum("is,isk->ik", trace[t], P_empa)
    alive = trace[:, :, EF_IDX] + trace[:, :, PE_IDX]  # (T+1, n_sim)
    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=wcc_method)
    weights = dw_c * wcc  # (T+1,)

    drug_spend_sampled = (
        alive * weights[:, None] * draws["c_drug"].to_numpy()[None, :]
    ).sum(axis=0) * cycle_length
    drug_spend_fixed = (alive * weights[:, None] * fixed_drug_price).sum(
        axis=0
    ) * cycle_length

    empa_costs_fixed = (
        result.costs[EMPA].to_numpy() - drug_spend_sampled + drug_spend_fixed
    )
    new_costs = result.costs.copy()
    new_costs[EMPA] = empa_costs_fixed
    new_params = result.params.copy()
    new_params["c_drug"] = fixed_drug_price
    return PSAResult(
        params=new_params,
        costs=new_costs,
        qalys=result.qalys,
        strategy_names=result.strategy_names,
        n_sim=result.n_sim,
    )


@dataclass(frozen=True)
class WaningSpec:
    """Linear treatment-effect waning.

    ``HR(t)`` is held fully at the base hazard ratio (``hr_event`` or
    ``hr_death``) for the first ``start_year`` years, interpolated
    linearly to 1.0 by ``end_year``, then 1.0 thereafter — i.e., the
    treatment effect tapers off after the trial follow-up horizon and
    eventually disappears.

    Defaults: trial follow-up of 3 years (EMPA-REG OUTCOME median 3.1
    yr) and full disappearance by year 10.
    """

    start_year: float = 3.0
    end_year: float = 10.0

    def __post_init__(self) -> None:
        if self.end_year <= self.start_year:
            raise ValueError("end_year must exceed start_year")

    def effect_fraction(self, t: np.ndarray) -> np.ndarray:
        """Fraction of the HR-vs-1.0 effect retained at time ``t``."""
        t = np.asarray(t, dtype=float)
        frac = (self.end_year - t) / (self.end_year - self.start_year)
        return np.clip(frac, 0.0, 1.0)

    def hr_path(self, hr_base: float, t: np.ndarray) -> np.ndarray:
        """Time-varying HR(t) = 1 - effect_fraction(t) * (1 - hr_base)."""
        return 1.0 - self.effect_fraction(t) * (1.0 - hr_base)


def _build_empa_transition_sequence(
    p: Mapping[str, Any], waning: WaningSpec
) -> np.ndarray:
    """Build the (T, 3, 3) per-cycle transition matrices for Empagliflozin
    under treatment-effect waning. SoC is unchanged across cycles.
    """
    cycle_length = float(p["cycle_length"])
    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
    r_EF_PE = float(p["r_EF_PE"])
    r_EF_D = float(p["r_EF_D"])
    r_PE_D = float(p["r_PE_D"])
    hr_event = float(p["hr_event"])
    hr_death = float(p["hr_death"])

    # HR(t) sampled at the start of each cycle (t = 0, 1, ..., T-1).
    t = np.arange(n_cycles, dtype=float) * cycle_length
    hr_event_t = waning.hr_path(hr_event, t)
    hr_death_t = waning.hr_path(hr_death, t)

    p_EF_PE = 1.0 - np.exp(-r_EF_PE * hr_event_t * cycle_length)
    p_EF_D = 1.0 - np.exp(-r_EF_D * hr_death_t * cycle_length)
    p_PE_D = 1.0 - np.exp(-r_PE_D * hr_death_t * cycle_length)

    P = np.zeros((n_cycles, 3, 3), dtype=float)
    P[:, 0, 0] = (1.0 - p_EF_D) * (1.0 - p_EF_PE)
    P[:, 0, 1] = (1.0 - p_EF_D) * p_EF_PE
    P[:, 0, 2] = p_EF_D
    P[:, 1, 1] = 1.0 - p_PE_D
    P[:, 1, 2] = p_PE_D
    P[:, 2, 2] = 1.0
    return P


def _acute_cost_contribution_sequence(
    trace: np.ndarray,
    P_seq: np.ndarray,
    c_acute_PE: float,
    dw_c: np.ndarray,
    wcc: np.ndarray,
) -> float:
    """Time-varying analogue of :func:`_acute_cost_contribution`.

    Uses ``P_seq[t][EF, PE]`` for the flow into PE during cycle ``t``.
    """
    T = trace.shape[0] - 1
    entries = trace[:T, EF_IDX] * P_seq[:, EF_IDX, PE_IDX]
    weights = dw_c[1 : T + 1] * wcc[1 : T + 1]
    return float(np.sum(entries * weights) * c_acute_PE)


def evaluate_scenario(
    params: ParameterInput,
    drug_price: Optional[float] = None,
    waning: Optional[WaningSpec] = None,
) -> dict[str, CaseStrategyResult]:
    """Run the case study under a scenario.

    Parameters
    ----------
    drug_price
        Annual empagliflozin acquisition cost (USD). If ``None``, uses
        the YAML ``c_drug`` (WAC). A typical illustrative net price is
        around $4,500/yr (~28% rebate off the WAC); real net prices are
        confidential and vary by payer.
    waning
        Optional :class:`WaningSpec`. If ``None``, the trial-derived
        ``hr_event`` and ``hr_death`` are applied for the full lifetime
        horizon (sustained-effect base case).

    Returns the same shape as :func:`evaluate_empagliflozin_case`.
    """
    p = _load_params(params)
    if drug_price is not None:
        p["c_drug"] = float(drug_price)

    if waning is None:
        return evaluate_empagliflozin_case(p)

    # --- waning case ----------------------------------------------------
    cycle_length = float(p["cycle_length"])
    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
    d_c = float(p["d_c"])
    d_e = float(p["d_e"])

    # SoC: unchanged across cycles -> build it via the existing builder
    # and run through the time-homogeneous engine. Empa under waning gets
    # the per-cycle sequence.
    model, c_acute_PE = build_empagliflozin_t2d(p)
    soc_strat = next(s for s in model.strategies if s.name == SOC)
    soc_r = evaluate_strategy(soc_strat, model)
    P_soc = np.asarray(soc_strat.transition_matrix, dtype=float)

    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=model.wcc_method)

    soc_acute = _acute_cost_contribution(soc_r["trace"], P_soc, c_acute_PE, dw_c, wcc)

    P_empa_seq = _build_empa_transition_sequence(p, waning)
    c_drug = float(p["c_drug"])
    c_EF = float(p["c_EF"])
    c_PE = float(p["c_PE"])
    c_D = float(p["c_D"])
    u_EF = float(p["u_EF"])
    u_PE = float(p["u_PE"])
    u_D = float(p["u_D"])

    empa_r = evaluate_sequence(
        name=EMPA,
        transition_sequence=P_empa_seq,
        state_costs=np.array([c_EF + c_drug, c_PE + c_drug, c_D]),
        state_utilities=np.array([u_EF, u_PE, u_D]),
        initial_distribution=np.array([1.0, 0.0, 0.0]),
        cycle_length=cycle_length,
        discount_rate_costs=d_c,
        discount_rate_qalys=d_e,
        wcc_method=model.wcc_method,
    )
    empa_acute = _acute_cost_contribution_sequence(
        empa_r["trace"], P_empa_seq, c_acute_PE, dw_c, wcc
    )

    return {
        SOC: {
            "total_cost": float(soc_r["total_cost"]) + soc_acute,
            "total_qaly": float(soc_r["total_qaly"]),
            "state_cost": float(soc_r["total_cost"]),
            "acute_event_cost": soc_acute,
            "trace": soc_r["trace"],
        },
        EMPA: {
            "total_cost": float(empa_r["total_cost"]) + empa_acute,
            "total_qaly": float(empa_r["total_qaly"]),
            "state_cost": float(empa_r["total_cost"]),
            "acute_event_cost": empa_acute,
            "trace": empa_r["trace"],
        },
    }


def scenario_icer(
    params: ParameterInput,
    drug_price: Optional[float] = None,
    waning: Optional[WaningSpec] = None,
) -> float:
    """Deterministic Empa-vs-SoC ICER under a scenario."""
    out = evaluate_scenario(params, drug_price=drug_price, waning=waning)
    dc = out[EMPA]["total_cost"] - out[SOC]["total_cost"]
    dq = out[EMPA]["total_qaly"] - out[SOC]["total_qaly"]
    if dq == 0:
        return float("nan")
    return float(dc / dq)


def breakeven_drug_price(
    params: ParameterInput,
    target_icer: float = 100_000.0,
    waning: Optional[WaningSpec] = None,
    bracket: Tuple[float, float] = (0.0, 50_000.0),
    tol: float = 1.0,
    max_iter: int = 80,
) -> float:
    """Annual drug price at which the Empa-vs-SoC ICER equals ``target_icer``.

    The incremental QALYs are independent of drug price (the cost only
    enters the cost ledger), and incremental cost is monotone increasing
    in drug price. So ``ICER(price)`` is a strictly increasing affine
    function — bisection converges fast on any bracketing interval where
    ``ICER(low) < target < ICER(high)``.
    """
    lo, hi = float(bracket[0]), float(bracket[1])
    f_lo = scenario_icer(params, drug_price=lo, waning=waning) - target_icer
    f_hi = scenario_icer(params, drug_price=hi, waning=waning) - target_icer
    if f_lo * f_hi > 0:
        raise ValueError(
            "target_icer not bracketed by [bracket[0], bracket[1]]: "
            f"ICER(low)={f_lo + target_icer:.2f}, ICER(high)={f_hi + target_icer:.2f}"
        )
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = scenario_icer(params, drug_price=mid, waning=waning) - target_icer
        if abs(f_mid) < tol or (hi - lo) < 1e-6:
            return float(mid)
        if f_mid * f_lo < 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return float(0.5 * (lo + hi))


# ---------------------------------------------------------------------------
# Scenario-aware PSA
# ---------------------------------------------------------------------------


def _batched_acute_cost_sequence(
    P_seq: np.ndarray,
    c_acute_PE: np.ndarray,
    init: np.ndarray,
    cycle_length: float,
    d_c: float,
    wcc_method: WCCMethod,
) -> np.ndarray:
    """Per-draw acute cost when the transition tensor is time-varying.

    ``P_seq`` has shape ``(n_cycles, n_sim, 3, 3)``. Trace propagation
    becomes ``trace[t+1] = einsum('is,isk->ik', trace[t], P_seq[t])``.
    """
    n_cycles, n_sim, _, _ = P_seq.shape
    trace = np.empty((n_cycles + 1, n_sim, 3), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = np.einsum("is,isk->ik", trace[t], P_seq[t])
    entries = trace[:n_cycles, :, EF_IDX] * P_seq[:, :, EF_IDX, PE_IDX]
    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=wcc_method)
    weights = dw_c[1 : n_cycles + 1] * wcc[1 : n_cycles + 1]
    return (entries * weights[:, None]).sum(axis=0) * c_acute_PE


def _simulate_sequence_batched(
    P_seq: np.ndarray,
    costs: np.ndarray,
    utils: np.ndarray,
    init: np.ndarray,
    cycle_length: float,
    d_c: float,
    d_e: float,
    wcc_method: WCCMethod,
) -> Tuple[np.ndarray, np.ndarray]:
    """Batched analogue of ``_simulate_strategy_batched`` for time-varying P.

    ``P_seq`` shape ``(n_cycles, n_sim, n_states, n_states)``.
    Returns ``(total_cost, total_qaly)`` of shape ``(n_sim,)``.
    """
    n_cycles, n_sim, n_states, _ = P_seq.shape
    trace = np.empty((n_cycles + 1, n_sim, n_states), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = np.einsum("is,isk->ik", trace[t], P_seq[t])
    cycle_costs = np.einsum("tis,is->ti", trace, costs) * cycle_length
    cycle_qalys = np.einsum("tis,is->ti", trace, utils) * cycle_length
    dw_c = _discount_weights(d_c, n_cycles, cycle_length)
    dw_e = _discount_weights(d_e, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, method=wcc_method)
    total_cost = (cycle_costs * (dw_c * wcc)[:, None]).sum(axis=0)
    total_qaly = (cycle_qalys * (dw_e * wcc)[:, None]).sum(axis=0)
    return total_cost, total_qaly


def run_empa_psa_scenario(
    base_params: ParameterInput,
    n_sim: int = 10_000,
    seed: int = 20260626,
    wcc_method: WCCMethod = "simpson_1_3",
    drug_price: Optional[float] = None,
    waning: Optional[WaningSpec] = None,
) -> PSAResult:
    """Vectorized PSA for the empagliflozin case under a scenario.

    Reuses the same sampler + Beta / Gamma / lognormal specs as
    :func:`run_empa_psa`. The sustained-effect path (no ``waning``)
    delegates straight to :func:`run_empa_psa` with the drug price
    optionally overridden; the waning path builds a time-varying
    transition tensor per draw and reuses the batched einsum
    propagation.
    """
    p = _load_params(base_params)
    _check_params(p)
    if drug_price is not None:
        p = dict(p)
        p["c_drug"] = float(drug_price)

    if waning is None and drug_price is None:
        # No time-variation and no drug-price override — straight reuse.
        return run_empa_psa(p, n_sim=n_sim, seed=seed, wcc_method=wcc_method)

    if waning is None:
        # Drug-price override without waning: sustained-effect PSA but with
        # ``c_drug`` held fixed at the scenario price across all draws,
        # rather than sampled around the WAC mean. The rest of the PSA
        # plumbing comes from run_empa_psa via reproduced code below.
        # (We can't simply delegate because run_empa_psa samples c_drug
        # from EMPA_PSA_SPECS, which would re-center it at the WAC mean.)
        assert drug_price is not None
        return _run_empa_psa_fixed_drug_price(
            p,
            n_sim=n_sim,
            seed=seed,
            wcc_method=wcc_method,
            fixed_drug_price=float(drug_price),
        )

    cycle_length = float(p["cycle_length"])
    n_cycles = int((float(p["n_age_max"]) - float(p["n_age_init"])) / cycle_length)
    d_c = float(p["d_c"])
    d_e = float(p["d_e"])

    r_EF_PE = float(p["r_EF_PE"])
    r_EF_D = float(p["r_EF_D"])

    draws = sample_psa_params(n_sim=n_sim, seed=seed, specs=EMPA_PSA_SPECS)
    hr_event = draws["hr_event"].to_numpy()
    hr_death = draws["hr_death"].to_numpy()
    r_PE_D = draws["r_PE_D"].to_numpy()
    c_drug = (
        draws["c_drug"].to_numpy()
        if drug_price is None
        else np.full(n_sim, float(drug_price))
    )
    c_EF = draws["c_EF"].to_numpy()
    c_PE = draws["c_PE"].to_numpy()
    c_acute_PE = draws["c_acute_PE"].to_numpy()
    u_EF = draws["u_EF"].to_numpy()
    u_PE = draws["u_PE"].to_numpy()
    zero = np.zeros(n_sim)
    init = np.array([1.0, 0.0, 0.0])

    # --- SoC: time-homogeneous; reuse the sustained-effect machinery ---
    p_EF_PE_soc = np.full(n_sim, 1.0 - np.exp(-r_EF_PE * cycle_length))
    p_EF_D_soc = np.full(n_sim, 1.0 - np.exp(-r_EF_D * cycle_length))
    p_PE_D_soc = 1.0 - np.exp(-r_PE_D * cycle_length)
    P_SoC = _empa_transition_tensor(p_EF_PE_soc, p_EF_D_soc, p_PE_D_soc)
    costs_SoC = np.column_stack([c_EF, c_PE, zero])
    util_SoC = np.column_stack([u_EF, u_PE, zero])
    cost_soc_state, qaly_soc = _simulate_strategy_batched(
        P_SoC,
        costs_SoC,
        util_SoC,
        init,
        n_cycles,
        cycle_length,
        d_c,
        d_e,
        wcc_method,
    )
    acute_soc = _batched_acute_cost(
        P_SoC,
        c_acute_PE,
        init,
        n_cycles,
        cycle_length,
        d_c,
        wcc_method,
    )

    # --- Empa: time-varying HRs ---
    t_grid = np.arange(n_cycles, dtype=float) * cycle_length
    frac = waning.effect_fraction(t_grid)  # (n_cycles,)
    hr_event_t = 1.0 - frac[:, None] * (1.0 - hr_event[None, :])  # (n_cycles, n_sim)
    hr_death_t = 1.0 - frac[:, None] * (1.0 - hr_death[None, :])

    p_EF_PE_empa = 1.0 - np.exp(
        -r_EF_PE * hr_event_t * cycle_length
    )  # (n_cycles, n_sim)
    p_EF_D_empa = 1.0 - np.exp(-r_EF_D * hr_death_t * cycle_length)
    p_PE_D_empa = 1.0 - np.exp(-r_PE_D[None, :] * hr_death_t * cycle_length)

    P_empa_seq = np.zeros((n_cycles, n_sim, 3, 3), dtype=float)
    P_empa_seq[:, :, 0, 0] = (1.0 - p_EF_D_empa) * (1.0 - p_EF_PE_empa)
    P_empa_seq[:, :, 0, 1] = (1.0 - p_EF_D_empa) * p_EF_PE_empa
    P_empa_seq[:, :, 0, 2] = p_EF_D_empa
    P_empa_seq[:, :, 1, 1] = 1.0 - p_PE_D_empa
    P_empa_seq[:, :, 1, 2] = p_PE_D_empa
    P_empa_seq[:, :, 2, 2] = 1.0

    costs_Empa = np.column_stack([c_EF + c_drug, c_PE + c_drug, zero])
    util_Empa = np.column_stack([u_EF, u_PE, zero])
    cost_empa_state, qaly_empa = _simulate_sequence_batched(
        P_empa_seq,
        costs_Empa,
        util_Empa,
        init,
        cycle_length,
        d_c,
        d_e,
        wcc_method,
    )
    acute_empa = _batched_acute_cost_sequence(
        P_empa_seq,
        c_acute_PE,
        init,
        cycle_length,
        d_c,
        wcc_method,
    )

    return PSAResult(
        params=draws,
        costs=pd.DataFrame(
            {SOC: cost_soc_state + acute_soc, EMPA: cost_empa_state + acute_empa}
        ),
        qalys=pd.DataFrame({SOC: qaly_soc, EMPA: qaly_empa}),
        strategy_names=STRATEGY_NAMES,
        n_sim=n_sim,
    )
