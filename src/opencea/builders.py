"""Model builders for published reference models.

This module turns published rate-based parameter specifications into the
generic ``CohortModel`` form consumed by the engine. The first (and
currently only) builder reproduces the DARTH Sick-Sicker time-independent
tutorial — keeping the rate-to-probability and competing-risks
construction in one place so that ``examples/sick_sicker.yaml`` can
encode the *published* rates rather than hand-derived probabilities.

Source of truth for every step below:

  Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. "An Introductory Tutorial on Cohort
  State-Transition Models in R Using a Cost-Effectiveness Analysis
  Example." Medical Decision Making, 2023;43(1):3-20.

  Repo: https://github.com/DARTH-git/cohort-modeling-tutorial-intro
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

from ._types import ParameterInput
from .model import CohortModel, Strategy


def rate_to_prob(r: float, t: float = 1.0) -> float:
    """Convert an instantaneous rate to a discrete-time probability.

    Faithful to ``rate_to_prob`` in DARTH ``R/Functions.R`` (lines 366-372)
    and to the equation in the manuscript section "Rates versus
    probabilities" (eq. rate-to-prob-ann):

        p = 1 - exp(-r * t)
    """
    if r < 0:
        raise ValueError("rate must be non-negative")
    return float(1.0 - np.exp(-r * t))


# Required keys for a DARTH Sick-Sicker parameter set. Matches the
# ``l_params_all`` list in ``analysis/cSTM_time_indep.R`` (lines 351-383).
_DARTH_REQUIRED_KEYS = {
    # transition rates & hazard ratios
    "r_HD",
    "r_HS1",
    "r_S1H",
    "r_S1S2",
    "hr_S1",
    "hr_S2",
    "hr_S1S2_trtB",
    # state costs (per cycle, before treatment add-ons)
    "c_H",
    "c_S1",
    "c_S2",
    "c_D",
    "c_trtA",
    "c_trtB",
    # state utilities
    "u_H",
    "u_S1",
    "u_S2",
    "u_D",
    "u_trtA",
    # cohort horizon and discounting
    "n_age_init",
    "n_age_max",
    "cycle_length",
    "d_c",
    "d_e",
}


def _check_params(params: Mapping[str, Any]) -> None:
    missing = _DARTH_REQUIRED_KEYS - set(params)
    if missing:
        raise KeyError(
            "DARTH Sick-Sicker parameter set is missing required keys: "
            f"{sorted(missing)}"
        )


def _build_transition_matrix(
    p_HD: float,
    p_HS1: float,
    p_S1H: float,
    p_S1S2: float,
    p_S1D: float,
    p_S2D: float,
) -> list[list[float]]:
    """Construct the 4x4 Sick-Sicker transition matrix (H, S1, S2, D).

    Reproduces the DARTH construction in ``analysis/cSTM_time_indep.R``
    lines 179-191 (and equivalently ``Functions_cSTM_time_indep.R``
    lines 49-66): transitions to non-death states are conditional on
    surviving the cycle, so they pick up a ``(1 - p_kD)`` factor.
    """
    # H, S1, S2, D
    P = [
        # From H
        [
            (1 - p_HD) * (1 - p_HS1),  # H -> H
            (1 - p_HD) * p_HS1,  # H -> S1
            0.0,  # H -> S2
            p_HD,  # H -> D
        ],
        # From S1
        [
            (1 - p_S1D) * p_S1H,  # S1 -> H
            (1 - p_S1D) * (1 - (p_S1H + p_S1S2)),  # S1 -> S1
            (1 - p_S1D) * p_S1S2,  # S1 -> S2
            p_S1D,  # S1 -> D
        ],
        # From S2
        [
            0.0,
            0.0,
            1 - p_S2D,
            p_S2D,
        ],
        # From D (absorbing)
        [0.0, 0.0, 0.0, 1.0],
    ]
    return P


def build_darth_sick_sicker(params: ParameterInput) -> CohortModel:
    """Construct a ``CohortModel`` for the DARTH Sick-Sicker tutorial.

    Parameters
    ----------
    params
        Either a mapping of parameter name -> value (matching the keys
        defined by the DARTH ``l_params_all`` list) or a path to a YAML
        file containing those parameters.

    Returns
    -------
    CohortModel
        A fully validated cohort model with four strategies (Standard
        of care, Strategy A, Strategy B, Strategy AB), ready to pass
        through ``opencea.engine.run_model``.
    """
    if isinstance(params, (str, Path)):
        import yaml

        with open(params, "r") as f:
            params = yaml.safe_load(f)
    assert isinstance(params, Mapping)
    _check_params(params)

    # --- rate -> probability conversion (DARTH lines 127-138) ----------
    cycle_length = float(params["cycle_length"])
    r_HD = float(params["r_HD"])
    r_HS1 = float(params["r_HS1"])
    r_S1H = float(params["r_S1H"])
    r_S1S2 = float(params["r_S1S2"])
    hr_S1 = float(params["hr_S1"])
    hr_S2 = float(params["hr_S2"])
    hr_S1S2_trtB = float(params["hr_S1S2_trtB"])

    # Mortality rates in Sick / Sicker states = hazard ratio * baseline.
    r_S1D = r_HD * hr_S1
    r_S2D = r_HD * hr_S2

    p_HS1 = rate_to_prob(r_HS1, cycle_length)
    p_S1H = rate_to_prob(r_S1H, cycle_length)
    p_S1S2 = rate_to_prob(r_S1S2, cycle_length)
    p_HD = rate_to_prob(r_HD, cycle_length)
    p_S1D = rate_to_prob(r_S1D, cycle_length)
    p_S2D = rate_to_prob(r_S2D, cycle_length)

    # Treatment B reduces the *rate* of S1 -> S2 (DARTH lines 140-148).
    r_S1S2_trtB = r_S1S2 * hr_S1S2_trtB
    p_S1S2_trtB = rate_to_prob(r_S1S2_trtB, cycle_length)

    # --- transition matrices (DARTH lines 178-203) ---------------------
    P_SoC = _build_transition_matrix(p_HD, p_HS1, p_S1H, p_S1S2, p_S1D, p_S2D)
    P_A = _build_transition_matrix(
        p_HD, p_HS1, p_S1H, p_S1S2, p_S1D, p_S2D
    )  # A leaves transitions unchanged
    P_B = _build_transition_matrix(p_HD, p_HS1, p_S1H, p_S1S2_trtB, p_S1D, p_S2D)
    P_AB = _build_transition_matrix(p_HD, p_HS1, p_S1H, p_S1S2_trtB, p_S1D, p_S2D)

    # --- state rewards (DARTH lines 244-284, manuscript lines 568-604) -
    # State cost / utility vectors here are PER CYCLE at cycle_length = 1
    # year (i.e. the annual values shown in Table 1 of the manuscript).
    # The engine multiplies by ``model.cycle_length`` internally, so the
    # YAML stays in the natural "annual" units.
    c_H = float(params["c_H"])
    c_S1 = float(params["c_S1"])
    c_S2 = float(params["c_S2"])
    c_D = float(params["c_D"])
    c_trtA = float(params["c_trtA"])
    c_trtB = float(params["c_trtB"])
    u_H = float(params["u_H"])
    u_S1 = float(params["u_S1"])
    u_S2 = float(params["u_S2"])
    u_D = float(params["u_D"])
    u_trtA = float(params["u_trtA"])

    # Costs: treatments add to both Sick and Sicker (lines 261-264, 271-283).
    costs_SoC = [c_H, c_S1, c_S2, c_D]
    costs_A = [c_H, c_S1 + c_trtA, c_S2 + c_trtA, c_D]
    costs_B = [c_H, c_S1 + c_trtB, c_S2 + c_trtB, c_D]
    costs_AB = [c_H, c_S1 + c_trtA + c_trtB, c_S2 + c_trtA + c_trtB, c_D]

    # Utilities: A and AB raise the Sick utility to u_trtA;
    # Sicker utility is unaffected because Strategy A treats only Sick.
    util_SoC = [u_H, u_S1, u_S2, u_D]
    util_A = [u_H, u_trtA, u_S2, u_D]
    util_B = [u_H, u_S1, u_S2, u_D]
    util_AB = [u_H, u_trtA, u_S2, u_D]

    n_cycles = int(
        (float(params["n_age_max"]) - float(params["n_age_init"])) / cycle_length
    )
    if n_cycles <= 0:
        raise ValueError("n_age_max must exceed n_age_init")

    return CohortModel(
        states=["H", "S1", "S2", "D"],
        initial_distribution=[1.0, 0.0, 0.0, 0.0],
        cycle_length=cycle_length,
        time_horizon=n_cycles,
        discount_rate_costs=float(params["d_c"]),
        discount_rate_qalys=float(params["d_e"]),
        wcc_method="simpson_1_3",
        strategies=[
            Strategy(
                name="Standard of care",
                transition_matrix=P_SoC,
                state_costs=costs_SoC,
                state_utilities=util_SoC,
            ),
            Strategy(
                name="Strategy A",
                transition_matrix=P_A,
                state_costs=costs_A,
                state_utilities=util_A,
            ),
            Strategy(
                name="Strategy B",
                transition_matrix=P_B,
                state_costs=costs_B,
                state_utilities=util_B,
            ),
            Strategy(
                name="Strategy AB",
                transition_matrix=P_AB,
                state_costs=costs_AB,
                state_utilities=util_AB,
            ),
        ],
    )
