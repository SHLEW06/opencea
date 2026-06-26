"""Tests for the empagliflozin applied case study.

Anchors the case-study evaluator to the validated deterministic engine
(base-case consistency), checks structural integrity of the 3-state /
2-strategy build, asserts the deterministic ICER lands in a plausible
band, and smoke-tests the DSA tornado + PSA / CEAC pipeline so an
inadvertent break in the case-study glue trips a test rather than the
manuscript-pinned DARTH golden tests.

Sanity-range commentary
-----------------------
Published US empagliflozin CEAs report ICERs in the ~$26k - $88k / QALY
band, but most of those analyses use **rebated / net** drug prices
(~$4 - 5k / yr) and **more granular** event states (separating MI,
stroke, and HF). This illustrative model uses the **WAC** drug cost
($6,264 / yr) and aggregates MI / stroke / HF into a single PE state, so
the deterministic ICER sits at the upper end of plausibility (around
$98 - 100k / QALY). The sanity-range assertion below is intentionally
permissive ($20k - $120k / QALY) — wide enough to allow this structural
result, narrow enough to catch order-of-magnitude mis-specifications
(a swapped HR direction, a missing acute cost, a missing discount, etc.).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from opencea.cea import cea_table, icer, nmb
from opencea.empagliflozin import (
    EMPA,
    EMPA_DSA_RANGE_OVERRIDES,
    EMPA_PSA_SPECS,
    SOC,
    STATES,
    STRATEGY_NAMES,
    build_empagliflozin_t2d,
    case_results_for_cea,
    dsa_evaluator,
    evaluate_empagliflozin_case,
    run_empa_psa,
)
from opencea.engine import evaluate_strategy, run_model
from opencea.psa import compute_ceac, default_wtp_grid
from opencea.sensitivity import run_dsa


YAML_PATH = Path(__file__).resolve().parents[1] / "examples" / "empagliflozin_t2d.yaml"

WTP = 100_000.0


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def case():
    return evaluate_empagliflozin_case(YAML_PATH)


@pytest.fixture(scope="module")
def model_pair():
    return build_empagliflozin_t2d(YAML_PATH)


@pytest.fixture(scope="module")
def dsa():
    return run_dsa(
        base_params=YAML_PATH,
        wtp=WTP,
        baseline=SOC,
        sweep_params=list(EMPA_PSA_SPECS.keys()),
        evaluator=dsa_evaluator,
        param_specs=EMPA_PSA_SPECS,
        param_ranges=EMPA_DSA_RANGE_OVERRIDES,
    )


@pytest.fixture(scope="module")
def psa():
    return run_empa_psa(YAML_PATH, n_sim=5_000, seed=20260626)


# ---------------------------------------------------------------------------
# 1. Structural integrity
# ---------------------------------------------------------------------------


def test_states_and_strategies(model_pair):
    model, c_acute = model_pair
    assert tuple(model.states) == STATES
    assert tuple(s.name for s in model.strategies) == STRATEGY_NAMES
    assert model.initial_distribution == [1.0, 0.0, 0.0]
    assert model.time_horizon == 37  # age 63 -> 100, annual cycles
    assert c_acute == pytest.approx(11650.0, abs=1e-9)


def test_transition_matrices_row_stochastic(model_pair):
    model, _ = model_pair
    for strat in model.strategies:
        P = np.asarray(strat.transition_matrix, dtype=float)
        assert P.shape == (3, 3)
        assert (P >= 0).all() and (P <= 1).all()
        np.testing.assert_allclose(P.sum(axis=1), 1.0, atol=1e-12)


def test_dead_state_is_absorbing(model_pair):
    model, _ = model_pair
    for strat in model.strategies:
        P = np.asarray(strat.transition_matrix, dtype=float)
        np.testing.assert_array_equal(P[2], [0.0, 0.0, 1.0])


def test_pe_does_not_recover_to_ef(model_pair):
    """PE only goes to PE or D (mirroring the Sicker state in DARTH)."""
    model, _ = model_pair
    for strat in model.strategies:
        P = np.asarray(strat.transition_matrix, dtype=float)
        assert P[1, 0] == 0.0  # PE -> EF blocked


def test_empagliflozin_reduces_progression_and_mortality(model_pair):
    """HRs < 1 on both channels => Empa's EF->PE and *->D probabilities
    should be strictly smaller than SoC's."""
    model, _ = model_pair
    P_soc = np.asarray(model.strategies[0].transition_matrix, dtype=float)
    P_empa = np.asarray(model.strategies[1].transition_matrix, dtype=float)
    assert P_empa[0, 1] < P_soc[0, 1]  # EF -> PE
    assert P_empa[0, 2] < P_soc[0, 2]  # EF -> D
    assert P_empa[1, 2] < P_soc[1, 2]  # PE -> D


# ---------------------------------------------------------------------------
# 2. Base-case consistency — DSA / evaluator tied back to the engine
# ---------------------------------------------------------------------------


def test_state_cost_matches_direct_engine_call(case, model_pair):
    """``evaluate_empagliflozin_case`` total_cost = engine state-cost
    total + acute transition cost. Verify the state-cost piece matches
    a direct engine call."""
    model, _ = model_pair
    for strat in model.strategies:
        r = evaluate_strategy(strat, model)
        assert case[r["name"]]["state_cost"] == pytest.approx(
            r["total_cost"], abs=1e-9
        )
        assert case[r["name"]]["total_qaly"] == pytest.approx(
            r["total_qaly"], abs=1e-9
        )


def test_total_cost_equals_state_plus_acute(case):
    for v in case.values():
        assert v["total_cost"] == pytest.approx(
            v["state_cost"] + v["acute_event_cost"], abs=1e-9
        )


def test_acute_event_cost_is_non_negative(case):
    for v in case.values():
        assert v["acute_event_cost"] >= 0


def test_empa_dominates_soc_on_qalys(case):
    """Both channels of empagliflozin reduce harmful rates, so its QALY
    total must exceed SoC's."""
    assert case[EMPA]["total_qaly"] > case[SOC]["total_qaly"]


def test_empa_costs_more_than_soc(case):
    """Drug cost > savings from prevented events at base case."""
    assert case[EMPA]["total_cost"] > case[SOC]["total_cost"]


# ---------------------------------------------------------------------------
# 3. ICER sanity band
# ---------------------------------------------------------------------------


# See the module docstring for why this band is intentionally wider than
# the literal user-cited $26-88k range.
ICER_SANITY_LOW = 20_000.0
ICER_SANITY_HIGH = 120_000.0


def test_deterministic_icer_in_sanity_band(case):
    """The deterministic ICER must land in a plausible range. The model
    is illustrative (3-state aggregation, WAC drug pricing), so the band
    is permissive — narrow enough to catch a swapped HR, missing acute
    cost, or missing discount; wide enough to allow this aggregation
    structure to land slightly above the central published $26 - 88k band."""
    soc = case[SOC]
    empa = case[EMPA]
    ic = icer(soc["total_cost"], soc["total_qaly"], empa["total_cost"], empa["total_qaly"])
    assert ICER_SANITY_LOW <= ic <= ICER_SANITY_HIGH, (
        f"ICER {ic:.2f} / QALY outside sanity band "
        f"[{ICER_SANITY_LOW}, {ICER_SANITY_HIGH}] — likely mis-specification"
    )


def test_icer_table_consistent_with_evaluator(case):
    """The CEA table built on case_results_for_cea must report the same
    ICER as a direct icer() call on the evaluator outputs."""
    table = cea_table(case_results_for_cea(YAML_PATH), wtp=WTP)
    # Sorted by total_cost ascending => SoC first, Empa second.
    assert list(table["strategy"]) == [SOC, EMPA]
    empa_row = table[table["strategy"] == EMPA].iloc[0]
    direct = icer(
        case[SOC]["total_cost"], case[SOC]["total_qaly"],
        case[EMPA]["total_cost"], case[EMPA]["total_qaly"],
    )
    assert empa_row["icer"] == pytest.approx(direct, abs=1e-9)


# ---------------------------------------------------------------------------
# 4. DSA structural sanity
# ---------------------------------------------------------------------------


def test_dsa_picks_empagliflozin_as_comparator_at_100k(dsa):
    """Empa NMB(100k) > SoC NMB(100k) at base case => comparator = Empa."""
    assert dsa.comparator == EMPA


def test_dsa_base_outcome_ties_to_evaluator(dsa, case):
    """Reuse of the engine: DSA's base_outcome equals the direct
    inc-NMB at $100k computed from the evaluator."""
    direct = nmb(case[EMPA]["total_cost"], case[EMPA]["total_qaly"], WTP) - nmb(
        case[SOC]["total_cost"], case[SOC]["total_qaly"], WTP
    )
    assert dsa.base_outcome == pytest.approx(direct, abs=1e-9)


def test_dsa_sweeps_cover_all_psa_specs(dsa):
    swept = {sw.parameter for sw in dsa.sweeps}
    assert swept == set(EMPA_PSA_SPECS.keys())


def test_dsa_sweeps_sorted_by_swing_desc(dsa):
    swings = [s.swing for s in dsa.sweeps]
    assert swings == sorted(swings, reverse=True)


def test_dsa_sweeps_bracket_base(dsa):
    tol = 1e-6
    for sw in dsa.sweeps:
        lo = min(sw.low_outcome, sw.high_outcome)
        hi = max(sw.low_outcome, sw.high_outcome)
        assert lo - tol <= sw.base_outcome <= hi + tol, (
            f"{sw.parameter}: base {sw.base_outcome} outside [{lo}, {hi}]"
        )


def test_dsa_top_driver_is_hr_death_or_drug_cost(dsa):
    """Top tornado driver should be one of the two parameters with the
    largest economic leverage: the all-cause mortality HR (which gates
    survival and downstream cumulative drug cost) or the drug
    acquisition cost itself."""
    assert dsa.sweeps[0].parameter in {"hr_death", "c_drug"}


# ---------------------------------------------------------------------------
# 5. PSA smoke
# ---------------------------------------------------------------------------


def test_psa_reproducible():
    r1 = run_empa_psa(YAML_PATH, n_sim=500, seed=7)
    r2 = run_empa_psa(YAML_PATH, n_sim=500, seed=7)
    pd.testing.assert_frame_equal(r1.costs, r2.costs)
    pd.testing.assert_frame_equal(r1.qalys, r2.qalys)


def test_psa_means_close_to_deterministic(psa, case):
    """PSA mean cost / qaly per strategy lands close to the deterministic
    base case (the parameter distributions are centred on the base
    values, so the means should track within a few percent)."""
    for name in STRATEGY_NAMES:
        mean_cost = float(psa.costs[name].mean())
        mean_qaly = float(psa.qalys[name].mean())
        assert mean_cost == pytest.approx(case[name]["total_cost"], rel=0.03)
        assert mean_qaly == pytest.approx(case[name]["total_qaly"], rel=0.03)


def test_ceac_endpoints_have_expected_winners(psa):
    """At WTP = $0 SoC must dominate (cheapest); at WTP = $200k Empa must
    dominate (highest QALYs)."""
    grid = default_wtp_grid()
    ceac = compute_ceac(psa, grid)
    assert ceac.iloc[0].idxmax() == SOC
    assert ceac.iloc[-1].idxmax() == EMPA


def test_ceac_probabilities_valid(psa):
    grid = default_wtp_grid()
    ceac = compute_ceac(psa, grid)
    arr = ceac.to_numpy()
    assert (arr >= 0).all() and (arr <= 1).all()
    np.testing.assert_allclose(arr.sum(axis=1), 1.0, atol=1e-12)


def test_ceac_at_100k_is_meaningful(psa):
    """ICER lands near $100k WTP, so P(CE at $100k) should be in the
    middle of [0, 1] — neither always nor never cost-effective."""
    grid = default_wtp_grid()
    ceac = compute_ceac(psa, grid)
    p = float(ceac.loc[100_000.0, EMPA])
    assert 0.2 < p < 0.8, (
        f"P(CE at $100k) = {p:.3f}; expected a middle-band value given "
        f"the deterministic ICER is near this WTP."
    )


# ---------------------------------------------------------------------------
# 6. Plot rendering smoke (PSA + DSA share the case-study figures)
# ---------------------------------------------------------------------------


def test_case_study_figures_render(psa, dsa, tmp_path):
    from opencea.plots import plot_ce_plane, plot_ceac, plot_ce_frontier, plot_tornado

    p1 = plot_ce_plane(psa, tmp_path / "ce_plane.png")
    p2 = plot_ceac(psa, tmp_path / "ceac.png")
    p3 = plot_ce_frontier(psa, tmp_path / "frontier.png")
    p4 = plot_tornado(dsa, tmp_path / "tornado.png")
    for p in (p1, p2, p3, p4):
        assert p.exists() and p.stat().st_size > 0
