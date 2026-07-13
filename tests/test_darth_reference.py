"""Golden-reference tests against the DARTH Sick-Sicker time-independent tutorial.

Source of truth (single authoritative reference for both parameters and
expected outputs):

  Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. "An Introductory Tutorial on Cohort
  State-Transition Models in R Using a Cost-Effectiveness Analysis
  Example." Medical Decision Making, 2023; 43(1):3-20.
  Repo: https://github.com/DARTH-git/cohort-modeling-tutorial-intro

The expected totals are taken verbatim from the typeset manuscript
(``manuscript/cSTM_Tutorial_Intro.tex`` in that repo, lines 746-749 of
the version archived on Zenodo / Medical Decision Making):

    Standard of care   $151,580    20.711
    Strategy A         $284,805    21.499   (dominated)
    Strategy B         $259,100    22.184
    Strategy AB        $378,875    23.137

    ICER B vs SoC = $72,988 / QALY
    ICER AB vs B  = $125,764 / QALY

These are the round-to-the-dollar totals printed in the manuscript's
Table 5 (Expected-outcomes-table) and Table 6 (CEA-table). The Python
implementation under test reproduces them to the cent when fed the
DARTH base-case parameters in ``examples/sick_sicker.yaml`` through
``opencea.builders.build_darth_sick_sicker``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from opencea.builders import build_darth_sick_sicker
from opencea.cea import cea_table, icer
from opencea.engine import gen_wcc, run_model, simulate_trace
from opencea.model import CohortModel

EXAMPLE_PATH = Path(__file__).resolve().parents[1] / "examples" / "sick_sicker.yaml"


# Expected outputs from the manuscript's Table 5 / Table 6.
# Costs are integer dollars; QALYs are reported to three decimals.
EXPECTED = {
    "Standard of care": {"cost": 151_580.0, "qaly": 20.711},
    "Strategy A": {"cost": 284_805.0, "qaly": 21.499},
    "Strategy B": {"cost": 259_100.0, "qaly": 22.184},
    "Strategy AB": {"cost": 378_875.0, "qaly": 23.137},
}

# Manuscript rounds costs to the nearest dollar and QALYs to 3 decimals.
# Allow a tight absolute tolerance reflecting the rounding precision plus
# fp noise: $1 on costs (≈ 7e-6 relative), 5e-4 on QALYs.
COST_ATOL = 1.0
QALY_ATOL = 5e-4


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def model() -> CohortModel:
    return build_darth_sick_sicker(EXAMPLE_PATH)


@pytest.fixture(scope="module")
def results(model: CohortModel):
    return run_model(model)


def _by_name(results, name):
    return next(r for r in results if r["name"] == name)


# ---------------------------------------------------------------------------
# 1. Strict structural properties (must hold to machine precision)
# ---------------------------------------------------------------------------


def test_transition_matrices_square_and_row_stochastic(model: CohortModel):
    n = len(model.states)
    for strat in model.strategies:
        P = np.asarray(strat.transition_matrix, dtype=float)
        assert P.shape == (n, n), f"{strat.name}: shape {P.shape} != ({n}, {n})"
        assert (P >= 0).all() and (P <= 1).all(), (
            f"{strat.name}: probabilities out of [0, 1]"
        )
        np.testing.assert_allclose(
            P.sum(axis=1),
            1.0,
            atol=1e-12,
            err_msg=f"{strat.name}: rows do not sum to 1",
        )


def test_simpson_1_3_wcc_matches_darth(model: CohortModel):
    """gen_wcc must replicate darthtools::gen_wcc for the base-case horizon.

    For n_cycles = 75 the DARTH function yields [1/3, 2/3, 4/3, 2/3, ...,
    4/3, 1/3] of length 76. (Endpoints 1/3; interior alternates 2/3 / 4/3
    starting at 2/3 at index 1 = R index 2 = "even".)
    """
    w = gen_wcc(model.time_horizon, method="simpson_1_3")
    assert w.shape == (model.time_horizon + 1,)
    assert w[0] == pytest.approx(1 / 3)
    assert w[-1] == pytest.approx(1 / 3)
    assert w[1] == pytest.approx(2 / 3)  # R index 2, even -> 2/3
    assert w[2] == pytest.approx(4 / 3)  # R index 3, odd  -> 4/3
    # Interior alternation
    interior = w[1:-1]
    np.testing.assert_allclose(interior[::2], 2 / 3, atol=1e-15)
    np.testing.assert_allclose(interior[1::2], 4 / 3, atol=1e-15)


def test_cohort_trace_shape_and_mass_conservation(model: CohortModel):
    soc = next(s for s in model.strategies if s.name == "Standard of care")
    trace = simulate_trace(soc, model)
    assert trace.shape == (model.time_horizon + 1, len(model.states))
    np.testing.assert_allclose(trace.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_array_equal(trace[0], np.array(model.initial_distribution))
    # Dead is absorbing -> monotone non-decreasing in the D column.
    d = model.states.index("D")
    assert np.all(np.diff(trace[:, d]) >= -1e-15)


# ---------------------------------------------------------------------------
# 2. Manuscript-pinned totals (Table 5 of cSTM_Tutorial_Intro)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", list(EXPECTED.keys()))
def test_total_discounted_cost_matches_manuscript(results, name):
    r = _by_name(results, name)
    expected_cost = EXPECTED[name]["cost"]
    assert r["total_cost"] == pytest.approx(expected_cost, abs=COST_ATOL), (
        f"{name}: computed cost {r['total_cost']:.4f} vs manuscript {expected_cost:.4f}"
    )


@pytest.mark.parametrize("name", list(EXPECTED.keys()))
def test_total_discounted_qaly_matches_manuscript(results, name):
    r = _by_name(results, name)
    expected_qaly = EXPECTED[name]["qaly"]
    assert r["total_qaly"] == pytest.approx(expected_qaly, abs=QALY_ATOL), (
        f"{name}: computed QALY {r['total_qaly']:.6f} vs manuscript {expected_qaly:.6f}"
    )


# ---------------------------------------------------------------------------
# 3. CEA: ICERs, dominance, frontier (Table 6 of cSTM_Tutorial_Intro)
# ---------------------------------------------------------------------------


def test_icer_b_vs_soc_matches_manuscript(results):
    soc = _by_name(results, "Standard of care")
    b = _by_name(results, "Strategy B")
    icer_val = icer(
        soc["total_cost"], soc["total_qaly"], b["total_cost"], b["total_qaly"]
    )
    # Manuscript prints $72,988 rounded to the dollar.
    assert icer_val == pytest.approx(72988.0, abs=1.0)


def test_icer_ab_vs_b_matches_manuscript(results):
    b = _by_name(results, "Strategy B")
    ab = _by_name(results, "Strategy AB")
    icer_val = icer(
        b["total_cost"], b["total_qaly"], ab["total_cost"], ab["total_qaly"]
    )
    assert icer_val == pytest.approx(125764.0, abs=1.0)


def test_strategy_a_is_strictly_dominated(results):
    """Manuscript Table 6: Strategy A is dominated; SoC / B / AB form the frontier."""
    table = cea_table(results)
    dominated = set(table.loc[table["dominated"], "strategy"])
    non_dominated = set(table.loc[~table["dominated"], "strategy"])
    assert dominated == {"Strategy A"}
    assert non_dominated == {"Standard of care", "Strategy B", "Strategy AB"}


def test_cea_table_orders_by_cost_and_exposes_required_columns(results):
    table = cea_table(results, wtp=100_000.0)
    required = {
        "strategy",
        "total_cost",
        "total_qaly",
        "dominated",
        "inc_cost",
        "inc_qaly",
        "icer",
        "nmb",
    }
    assert required.issubset(table.columns)
    assert list(table["total_cost"]) == sorted(table["total_cost"])
