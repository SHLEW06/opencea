"""Cohort trace simulation, discounting, and within-cycle correction (WCC).

The engine deliberately mirrors the DARTH time-independent tutorial's
numerical conventions so that golden tests can validate against the
published reference outputs:

  - Per-state cost / utility vectors are multiplied by ``cycle_length``
    when accumulating per-cycle rewards (cf. ``cSTM_time_indep.R`` lines
    246-284 in the DARTH intro repo).
  - Cycle-specific discount weights are computed from the DARTH form
    ``1 / (1 + d * cycle_length) ^ t`` (cf. ``cSTM_time_indep.R`` lines
    123-124) so that the engine and a rerun of the DARTH code agree
    exactly for any cycle length.
  - Three within-cycle-correction methods are supported, matching
    ``darthtools::gen_wcc``: Simpson's 1/3 rule (the DARTH default),
    half-cycle correction, and none.

Conventions
-----------
* The cohort trace has shape ``(time_horizon + 1, n_states)``. Row 0 is
  the initial distribution; row ``t`` is the distribution at the
  *start* of cycle ``t``. The propagation rule is
  ``trace[t+1] = trace[t] @ P``.
* Per-cycle rewards (costs / utilities) are interpreted as accruing to
  occupants of each state during the cycle.
* Total reward = sum over t of (per-cycle reward) * (discount weight) *
  (within-cycle-correction weight).
"""

from __future__ import annotations

from typing import Dict, List, Literal

import numpy as np

from .model import CohortModel, Strategy

WCCMethod = Literal["simpson_1_3", "half_cycle", "none"]


def simulate_trace(strategy: Strategy, model: CohortModel) -> np.ndarray:
    """Simulate the cohort trace for one strategy.

    Returns an array of shape ``(time_horizon + 1, n_states)`` whose
    rows are the cohort distributions at cycles ``0, 1, ..., T``.
    """
    P = np.asarray(strategy.transition_matrix, dtype=float)
    init = np.asarray(model.initial_distribution, dtype=float)
    n_cycles = model.time_horizon
    n_states = init.size

    trace = np.zeros((n_cycles + 1, n_states), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = trace[t] @ P
    return trace


def _discount_weights(rate: float, n_cycles: int, cycle_length: float) -> np.ndarray:
    """DARTH-style discount weights: 1 / (1 + d * cycle_length) ** t.

    Mirrors ``cSTM_time_indep.R`` lines 123-124. Reduces to the standard
    ``1 / (1 + d) ** t`` when ``cycle_length == 1``.
    """
    t = np.arange(n_cycles + 1, dtype=float)
    return 1.0 / (1.0 + rate * cycle_length) ** t


def gen_wcc(n_cycles: int, method: WCCMethod = "simpson_1_3") -> np.ndarray:
    """Within-cycle-correction weights.

    Implementation faithful to ``darthtools::gen_wcc`` (DARTH intro repo,
    ``R/Functions.R`` lines 321-342):

    * ``"simpson_1_3"`` — Simpson's 1/3 rule. R uses 1-based indices, so
      after translating to 0-based Python indices ``k``:
      - ``k`` even (R index ``k+1`` odd) -> 4/3
      - ``k`` odd (R index ``k+1`` even) -> 2/3
      - first (``k = 0``) and last (``k = n_cycles``) overridden to 1/3.
    * ``"half_cycle"`` — trapezoidal rule, endpoints 0.5, interior 1.
    * ``"none"`` — uniform weights of 1.
    """
    if n_cycles <= 0:
        raise ValueError("n_cycles must be positive")

    n = n_cycles + 1

    if method == "simpson_1_3":
        # k even -> 4/3, k odd -> 2/3, then override endpoints to 1/3.
        k = np.arange(n)
        w = np.where(k % 2 == 0, 4.0 / 3.0, 2.0 / 3.0)
        w[0] = 1.0 / 3.0
        w[-1] = 1.0 / 3.0
        return w

    if method == "half_cycle":
        w = np.ones(n)
        w[0] = 0.5
        w[-1] = 0.5
        return w

    if method == "none":
        return np.ones(n)

    raise ValueError(f"unknown WCC method: {method!r}")


def evaluate_strategy(strategy: Strategy, model: CohortModel) -> Dict[str, object]:
    """Compute the cohort trace and total discounted costs / QALYs for one strategy.

    Returns a dict with keys:

    - ``name``: strategy name
    - ``trace``: ``np.ndarray`` of shape ``(T+1, n_states)``
    - ``cycle_costs``: per-cycle expected costs over the cohort, ``(T+1,)``
    - ``cycle_qalys``: per-cycle expected QALYs over the cohort, ``(T+1,)``
    - ``total_cost``: total discounted cost (scalar)
    - ``total_qaly``: total discounted QALYs (scalar)
    """
    trace = simulate_trace(strategy, model)

    costs = np.asarray(strategy.state_costs, dtype=float)
    utilities = np.asarray(strategy.state_utilities, dtype=float)

    # DARTH convention: per-state rewards are scaled by cycle_length before
    # the cohort dot product, so the resulting per-cycle reward is the
    # expected reward accumulated over one cycle of length cycle_length.
    cycle_costs = trace @ (costs * model.cycle_length)
    cycle_qalys = trace @ (utilities * model.cycle_length)

    dw_c = _discount_weights(
        model.discount_rate_costs, model.time_horizon, model.cycle_length
    )
    dw_e = _discount_weights(
        model.discount_rate_qalys, model.time_horizon, model.cycle_length
    )
    wcc = gen_wcc(model.time_horizon, model.wcc_method)

    total_cost = float(np.sum(cycle_costs * dw_c * wcc))
    total_qaly = float(np.sum(cycle_qalys * dw_e * wcc))

    return {
        "name": strategy.name,
        "trace": trace,
        "cycle_costs": cycle_costs,
        "cycle_qalys": cycle_qalys,
        "total_cost": total_cost,
        "total_qaly": total_qaly,
    }


def run_model(model: CohortModel) -> List[Dict[str, object]]:
    """Evaluate every strategy in the model and return their results."""
    return [evaluate_strategy(s, model) for s in model.strategies]


# ---------------------------------------------------------------------------
# Time-varying transitions (additive helpers)
#
# The Strategy / CohortModel schema is intentionally time-homogeneous (one
# transition matrix per strategy). These helpers add the minimum needed to
# evaluate a strategy whose transition matrix varies per cycle — used by
# the empagliflozin case study to model treatment-effect waning. The
# validated time-homogeneous path above is untouched.
# ---------------------------------------------------------------------------


def simulate_trace_sequence(
    transition_sequence: np.ndarray,
    initial_distribution: np.ndarray,
) -> np.ndarray:
    """Cohort trace under a per-cycle sequence of transition matrices.

    Parameters
    ----------
    transition_sequence
        Array of shape ``(n_cycles, n_states, n_states)``. ``P[t]`` is
        the transition matrix applied during cycle ``t`` (i.e., the one
        that takes ``trace[t]`` to ``trace[t + 1]``).
    initial_distribution
        Shape ``(n_states,)``.

    Returns
    -------
    np.ndarray
        Shape ``(n_cycles + 1, n_states)``, same convention as
        :func:`simulate_trace`.
    """
    Pseq = np.asarray(transition_sequence, dtype=float)
    if Pseq.ndim != 3 or Pseq.shape[1] != Pseq.shape[2]:
        raise ValueError(
            f"transition_sequence must have shape (T, n, n); got {Pseq.shape}"
        )
    init = np.asarray(initial_distribution, dtype=float)
    n_cycles, n_states, _ = Pseq.shape
    if init.size != n_states:
        raise ValueError(
            f"initial_distribution has {init.size} entries; expected {n_states}"
        )

    trace = np.empty((n_cycles + 1, n_states), dtype=float)
    trace[0] = init
    for t in range(n_cycles):
        trace[t + 1] = trace[t] @ Pseq[t]
    return trace


def evaluate_sequence(
    name: str,
    transition_sequence: np.ndarray,
    state_costs: np.ndarray,
    state_utilities: np.ndarray,
    initial_distribution: np.ndarray,
    cycle_length: float,
    discount_rate_costs: float,
    discount_rate_qalys: float,
    wcc_method: WCCMethod = "simpson_1_3",
) -> Dict[str, object]:
    """Time-varying analogue of :func:`evaluate_strategy`.

    Uses the same per-cycle reward / discount / WCC accumulation as the
    time-homogeneous evaluator so totals are directly comparable when
    ``transition_sequence`` happens to be a stack of identical matrices.
    """
    trace = simulate_trace_sequence(transition_sequence, initial_distribution)
    n_cycles = trace.shape[0] - 1

    costs = np.asarray(state_costs, dtype=float)
    utilities = np.asarray(state_utilities, dtype=float)

    cycle_costs = trace @ (costs * cycle_length)
    cycle_qalys = trace @ (utilities * cycle_length)

    dw_c = _discount_weights(discount_rate_costs, n_cycles, cycle_length)
    dw_e = _discount_weights(discount_rate_qalys, n_cycles, cycle_length)
    wcc = gen_wcc(n_cycles, wcc_method)

    total_cost = float(np.sum(cycle_costs * dw_c * wcc))
    total_qaly = float(np.sum(cycle_qalys * dw_e * wcc))

    return {
        "name": name,
        "trace": trace,
        "cycle_costs": cycle_costs,
        "cycle_qalys": cycle_qalys,
        "total_cost": total_cost,
        "total_qaly": total_qaly,
    }
