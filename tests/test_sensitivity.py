"""Tests for one-way deterministic sensitivity analysis and the tornado plot.

These tests anchor the DSA to the same validated deterministic engine the
manuscript-pinned tests already cover: the base-case outcome must equal
the value you would compute by running the engine directly on the YAML
parameters. Beyond that, the tests assert structural properties of the
tornado (bracketing, non-negative swings, descending sort) and a sanity
ordering — a parameter with real uncertainty must outrank a tightly
specified one (hr_S1 has sdlog = 0.01 by construction).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from opencea.builders import build_darth_sick_sicker
from opencea.cea import nmb
from opencea.engine import run_model
from opencea.psa import PSA_PARAM_SPECS
from opencea.sensitivity import compute_parameter_ranges, run_dsa

EXAMPLE_PATH = Path(__file__).resolve().parents[1] / "examples" / "sick_sicker.yaml"

WTP = 100_000.0


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dsa():
    """Default DSA: optimal-vs-SoC at WTP=$100k across every PSA parameter."""
    return run_dsa(EXAMPLE_PATH, wtp=WTP)


@pytest.fixture(scope="module")
def deterministic_results():
    model = build_darth_sick_sicker(EXAMPLE_PATH)
    return {r["name"]: r for r in run_model(model)}


# ---------------------------------------------------------------------------
# Base-case consistency: ties DSA back to the validated engine
# ---------------------------------------------------------------------------


def test_optimal_comparator_at_base_is_strategy_b(dsa):
    """At WTP=$100k the deterministic NMB-optimal strategy vs SoC is B
    (ICER B vs SoC = $72,988 < $100k; ICER AB vs B = $125,764 > $100k)."""
    assert dsa.comparator == "Strategy B"


def test_dsa_base_outcome_matches_engine_directly(dsa, deterministic_results):
    """The DSA's base_outcome must equal the incremental NMB computed by
    running the validated engine on the YAML parameters."""
    soc = deterministic_results["Standard of care"]
    b = deterministic_results["Strategy B"]
    expected = nmb(b["total_cost"], b["total_qaly"], WTP) - nmb(
        soc["total_cost"], soc["total_qaly"], WTP
    )
    assert dsa.base_outcome == pytest.approx(expected, abs=1e-9)


def test_dsa_endpoint_eval_uses_same_builder_engine_path(dsa, deterministic_results):
    """Evaluating one endpoint manually — modifying a single parameter in the
    YAML dict, then running build_darth_sick_sicker + run_model — must
    reproduce the DSA's recorded endpoint outcome exactly. This is the
    sharpest tie between the DSA driver and the validated engine: the
    DSA goes through the same code path the deterministic golden tests
    exercise."""
    import yaml

    bp = yaml.safe_load(EXAMPLE_PATH.read_text())
    target_param = "c_trtB"  # arbitrary high-swing parameter
    sw = next(s for s in dsa.sweeps if s.parameter == target_param)

    bp_low = dict(bp)
    bp_low[target_param] = sw.low_value
    res = {r["name"]: r for r in run_model(build_darth_sick_sicker(bp_low))}
    expected_low = nmb(
        res["Strategy B"]["total_cost"], res["Strategy B"]["total_qaly"], WTP
    ) - nmb(
        res["Standard of care"]["total_cost"],
        res["Standard of care"]["total_qaly"],
        WTP,
    )
    assert sw.low_outcome == pytest.approx(expected_low, abs=1e-9)


# ---------------------------------------------------------------------------
# Structural properties of the tornado
# ---------------------------------------------------------------------------


def test_sweeps_cover_all_psa_parameters(dsa):
    """Every PSA-sampled parameter should be in the tornado (the fixed
    c_D / u_D are excluded by construction)."""
    swept = {sw.parameter for sw in dsa.sweeps}
    assert swept == set(PSA_PARAM_SPECS.keys())


def test_sweeps_sorted_by_swing_descending(dsa):
    swings = [sw.swing for sw in dsa.sweeps]
    assert swings == sorted(swings, reverse=True)


def test_all_swings_non_negative(dsa):
    for sw in dsa.sweeps:
        assert sw.swing >= 0


def test_each_sweep_brackets_base_outcome(dsa):
    """For each parameter, the base outcome lies within [min, max] of the
    low / high outcomes. Holds because the DSA ranges are extended to
    include the deterministic base when it falls outside [P2.5, P97.5]
    (relevant for u_H)."""
    tol = 1e-6
    for sw in dsa.sweeps:
        lo = min(sw.low_outcome, sw.high_outcome)
        hi = max(sw.low_outcome, sw.high_outcome)
        assert lo - tol <= sw.base_outcome <= hi + tol, (
            f"{sw.parameter}: base {sw.base_outcome} outside [{lo}, {hi}]"
        )


def test_each_sweep_brackets_base_parameter_value(dsa):
    """Each parameter's base value must lie between its low and high
    sweep values (this is what makes 'extend to include base' work)."""
    for sw in dsa.sweeps:
        assert sw.low_value <= sw.base_value <= sw.high_value, (
            f"{sw.parameter}: base value {sw.base_value} not in "
            f"[{sw.low_value}, {sw.high_value}]"
        )


# ---------------------------------------------------------------------------
# Sensitivity ordering sanity
# ---------------------------------------------------------------------------


def test_high_uncertainty_param_outranks_tight_hazard_ratio(dsa):
    """A genuinely uncertain parameter (cost of treatment B, the dominant
    driver of incremental NMB vs SoC) must produce a larger swing than
    hr_S1, whose lognormal sdlog of 0.01 makes it essentially fixed."""
    by_name = {sw.parameter: sw for sw in dsa.sweeps}
    assert by_name["c_trtB"].swing > by_name["hr_S1"].swing
    # Comfortable margin: c_trtB swings tens of thousands of dollars while
    # hr_S1 barely moves the outcome.
    assert by_name["c_trtB"].swing > 100 * by_name["hr_S1"].swing


def test_c_trtA_and_u_trtA_have_zero_swing_in_b_vs_soc(dsa):
    """Strategy-A-only parameters cannot affect Strategy B vs SoC.

    c_trtA and u_trtA only enter Strategy A and Strategy AB. With Strategy
    B as the comparator (chosen at base because B is NMB-optimal at $100k)
    these parameters must have exactly zero swing — a strong check that
    the DSA is wiring parameters into the right strategies."""
    by_name = {sw.parameter: sw for sw in dsa.sweeps}
    assert by_name["c_trtA"].swing == pytest.approx(0.0, abs=1e-6)
    assert by_name["u_trtA"].swing == pytest.approx(0.0, abs=1e-6)


def test_c_trtB_is_top_driver(dsa):
    """At WTP=$100k the cost of treatment B is the largest single driver
    of incremental NMB B vs SoC (it scales linearly with person-cycles in
    S1/S2 under B; no offsetting term exists in SoC)."""
    assert dsa.sweeps[0].parameter == "c_trtB"


# ---------------------------------------------------------------------------
# Parameter-range plumbing
# ---------------------------------------------------------------------------


def test_compute_parameter_ranges_reproducible():
    """Same seed -> identical ranges."""
    bp = {**{k: v.target_mean for k, v in PSA_PARAM_SPECS.items()}}
    r1 = compute_parameter_ranges(bp, seed=42, n_sample=10_000)
    r2 = compute_parameter_ranges(bp, seed=42, n_sample=10_000)
    assert r1 == r2


def test_compute_parameter_ranges_extends_for_u_H():
    """u_H has deterministic base 1.0 > 97.5%ile of Beta(200, 3) ~ 0.998.
    With extension enabled the range high must equal the base (1.0)."""
    import yaml

    bp = yaml.safe_load(EXAMPLE_PATH.read_text())
    ranges = compute_parameter_ranges(bp, extend_to_include_base=True)
    lo, hi = ranges["u_H"]
    assert hi == pytest.approx(1.0, abs=1e-9)
    assert lo < 1.0


def test_dsa_to_dataframe_columns(dsa):
    df = dsa.to_dataframe()
    required = {
        "parameter",
        "base_value",
        "low_value",
        "high_value",
        "low_outcome",
        "high_outcome",
        "base_outcome",
        "swing",
    }
    assert required.issubset(df.columns)
    # Sorted by swing descending in the row order
    assert list(df["swing"]) == sorted(df["swing"], reverse=True)


# ---------------------------------------------------------------------------
# Tornado plot smoke test
# ---------------------------------------------------------------------------


def test_tornado_plot_renders(dsa, tmp_path):
    from opencea.plots import plot_tornado

    out = plot_tornado(dsa, tmp_path / "tornado.png")
    assert out.exists() and out.stat().st_size > 0


def test_tornado_plot_top_n(dsa, tmp_path):
    """top_n cap should still produce a valid file with the largest swings."""
    from opencea.plots import plot_tornado

    out = plot_tornado(dsa, tmp_path / "tornado_top5.png", top_n=5)
    assert out.exists() and out.stat().st_size > 0
