"""Structural and statistical validation for the PSA layer.

Exact draw-for-draw matching with the DARTH R reference is impossible
across RNGs (numpy and R use different streams), so these tests validate
*properties* instead of values:

  - the sampler is reproducible given a seed,
  - the per-parameter Monte Carlo means recover the published targets
    in ``generate_psa_params()`` (DARTH ``R/Functions_cSTM_time_indep.R``
    lines 278-309) within ~3x MC error,
  - the PSA mean cost / QALY per strategy lands within ~1-2% of the
    deterministic Table 5 totals (QALYs come in slightly low because
    u_H samples around 0.985 rather than 1.0),
  - the CEAC is a valid probability matrix (entries in [0, 1], rows sum
    to 1) with the expected ordering at the extremes of the WTP grid
    (SoC dominates at WTP = 0; AB dominates at very high WTP),
  - Strategy A is dominated *in expectation* (matches the deterministic
    Table 6 result).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from opencea.psa import (
    PSA_PARAM_SPECS,
    compute_ceac,
    compute_nmb,
    default_wtp_grid,
    expected_nmb_frontier,
    incremental_vs_baseline,
    run_psa,
    sample_psa_params,
)

EXAMPLE_PATH = Path(__file__).resolve().parents[1] / "examples" / "sick_sicker.yaml"

# Used by the "PSA mean vs deterministic" tests. Numbers from
# tests/test_darth_reference.py — the manuscript's Table 5.
DETERMINISTIC = {
    "Standard of care": {"cost": 151_580.0, "qaly": 20.711},
    "Strategy A": {"cost": 284_805.0, "qaly": 21.499},
    "Strategy B": {"cost": 259_100.0, "qaly": 22.184},
    "Strategy AB": {"cost": 378_875.0, "qaly": 23.137},
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


PSA_SEED = 20260625
PSA_N_SIM = 10_000  # default used in the production run


@pytest.fixture(scope="module")
def psa_result():
    return run_psa(EXAMPLE_PATH, n_sim=PSA_N_SIM, seed=PSA_SEED)


@pytest.fixture(scope="module")
def wtp_grid():
    return default_wtp_grid()


# ---------------------------------------------------------------------------
# 1. Reproducibility
# ---------------------------------------------------------------------------


def test_sampler_is_reproducible_under_seed():
    """Same seed -> bit-identical draws across calls."""
    a = sample_psa_params(n_sim=1000, seed=12345)
    b = sample_psa_params(n_sim=1000, seed=12345)
    pd.testing.assert_frame_equal(a, b)


def test_sampler_differs_with_different_seeds():
    """Different seeds -> different draws (sanity)."""
    a = sample_psa_params(n_sim=1000, seed=1)
    b = sample_psa_params(n_sim=1000, seed=2)
    assert not a.equals(b)


def test_full_psa_run_is_reproducible_under_seed():
    """Identical seed -> identical PSA cost / QALY arrays."""
    r1 = run_psa(EXAMPLE_PATH, n_sim=500, seed=999)
    r2 = run_psa(EXAMPLE_PATH, n_sim=500, seed=999)
    pd.testing.assert_frame_equal(r1.costs, r2.costs)
    pd.testing.assert_frame_equal(r1.qalys, r2.qalys)


# ---------------------------------------------------------------------------
# 2. Sampler statistical correctness — wires the right R distribution to
#    the right numpy parameterization
# ---------------------------------------------------------------------------


# Use a larger draw count for this test so the central-limit error band
# is tight enough to detect a swapped scale / rate parameterization.
_MEAN_TEST_N = 50_000


@pytest.fixture(scope="module")
def big_draws():
    return sample_psa_params(n_sim=_MEAN_TEST_N, seed=20260625)


def _analytic_variance(name: str) -> float:
    """Closed-form variance of each PSA distribution for an MC error band."""
    spec = PSA_PARAM_SPECS[name]
    if spec.family == "gamma":
        shape, scale = spec.params["shape"], spec.params["scale"]
        return shape * scale**2
    if spec.family == "beta":
        a, b = spec.params["a"], spec.params["b"]
        return (a * b) / ((a + b) ** 2 * (a + b + 1.0))
    if spec.family == "lognormal":
        # Var of lognormal(mu, sigma) = (exp(sigma^2) - 1) * exp(2 mu + sigma^2)
        mu, sigma = spec.params["mean"], spec.params["sigma"]
        return (np.exp(sigma**2) - 1.0) * np.exp(2.0 * mu + sigma**2)
    raise ValueError(spec.family)


@pytest.mark.parametrize("name", list(PSA_PARAM_SPECS.keys()))
def test_sampled_mean_recovers_target(name, big_draws):
    """Sample mean must land within ~4 standard errors of the analytic mean.

    For lognormal parameters the manuscript / DARTH spec the *median*
    rather than the mean; we still compare the sample mean to the
    distribution's analytic mean (= exp(mu + sigma^2/2)) — that's the
    quantity the SE band describes.
    """
    spec = PSA_PARAM_SPECS[name]
    if spec.family == "lognormal":
        analytic_mean = float(
            np.exp(spec.params["mean"] + 0.5 * spec.params["sigma"] ** 2)
        )
    else:
        analytic_mean = spec.target_mean

    var = _analytic_variance(name)
    se = (var / _MEAN_TEST_N) ** 0.5
    sample_mean = float(big_draws[name].mean())

    # 4 SE one-sided -> P(|Z| > 4) ~ 6e-5, comfortable margin to avoid
    # flaky failures while still catching a swapped parameterization
    # (those typically shift the mean by orders of magnitude).
    assert abs(sample_mean - analytic_mean) < 4.0 * se, (
        f"{name}: sample mean {sample_mean:.6g} vs analytic {analytic_mean:.6g} "
        f"(4 SE band = {4 * se:.6g})"
    )


def test_fixed_parameters_are_constant_zero():
    """c_D and u_D are not sampled; they should be exactly 0 in every draw."""
    draws = sample_psa_params(n_sim=200, seed=1)
    assert (draws["c_D"] == 0).all()
    assert (draws["u_D"] == 0).all()


def test_lognormal_medians_match_specification():
    """For lognormal parameters the *median* equals exp(mu); verify the
    sample median tracks it."""
    draws = sample_psa_params(n_sim=_MEAN_TEST_N, seed=2026)
    for name in ("hr_S1", "hr_S2", "hr_S1S2_trtB"):
        spec = PSA_PARAM_SPECS[name]
        target_median = float(np.exp(spec.params["mean"]))
        sample_median = float(np.median(draws[name]))
        # Lognormal medians are very tight at these sigmas; 0.5% rel is loose.
        assert sample_median == pytest.approx(target_median, rel=0.005), (
            f"{name}: sample median {sample_median:.5f} vs target {target_median:.5f}"
        )


# ---------------------------------------------------------------------------
# 3. PSA mean cost / QALY vs deterministic Table 5
# ---------------------------------------------------------------------------


# PSA means won't match Table 5 exactly:
#   - u_H samples around 0.9852 (vs deterministic 1.0), so all QALY totals
#     come in roughly 1% low.
#   - Other parameters are symmetric / centred on the deterministic values,
#     so costs sit within ~0.5%.
# Tolerances cover MC error + the structural u_H gap.
COST_REL_TOL = 0.02  # 2% relative for costs
QALY_REL_TOL = 0.025  # 2.5% relative for QALYs (u_H drift = ~1.5% on its own)


@pytest.mark.parametrize("name", list(DETERMINISTIC.keys()))
def test_psa_mean_cost_close_to_table5(psa_result, name):
    mean_cost = float(psa_result.costs[name].mean())
    expected = DETERMINISTIC[name]["cost"]
    assert mean_cost == pytest.approx(expected, rel=COST_REL_TOL), (
        f"{name}: PSA mean cost {mean_cost:.2f} vs deterministic {expected:.2f}"
    )


@pytest.mark.parametrize("name", list(DETERMINISTIC.keys()))
def test_psa_mean_qaly_close_to_table5(psa_result, name):
    mean_qaly = float(psa_result.qalys[name].mean())
    expected = DETERMINISTIC[name]["qaly"]
    assert mean_qaly == pytest.approx(expected, rel=QALY_REL_TOL), (
        f"{name}: PSA mean QALY {mean_qaly:.4f} vs deterministic {expected:.4f}"
    )


def test_psa_qalys_come_in_slightly_below_deterministic(psa_result):
    """Documented structural drift from u_H sampling around 0.985."""
    for name in DETERMINISTIC:
        mean_qaly = float(psa_result.qalys[name].mean())
        assert mean_qaly < DETERMINISTIC[name]["qaly"], (
            f"{name}: expected PSA mean QALY < deterministic, got {mean_qaly:.4f}"
        )


# ---------------------------------------------------------------------------
# 4. CEAC sanity
# ---------------------------------------------------------------------------


def test_ceac_probabilities_valid(psa_result, wtp_grid):
    """At every WTP each strategy probability is in [0, 1] and the row
    sums to 1 (each draw has exactly one NMB-maximizing strategy)."""
    ceac = compute_ceac(psa_result, wtp_grid)
    arr = ceac.to_numpy()
    assert (arr >= 0).all() and (arr <= 1).all()
    np.testing.assert_allclose(arr.sum(axis=1), 1.0, atol=1e-12)


def test_ceac_endpoints_have_expected_winners(psa_result, wtp_grid):
    """SoC has the highest probability at very low WTP (cost matters most);
    AB has the highest probability at very high WTP (QALYs matter most)."""
    ceac = compute_ceac(psa_result, wtp_grid)
    row_low = ceac.iloc[0]
    row_high = ceac.iloc[-1]
    assert row_low.idxmax() == "Standard of care", (
        f"At WTP={wtp_grid[0]}, expected SoC highest; got {row_low.to_dict()}"
    )
    assert row_high.idxmax() == "Strategy AB", (
        f"At WTP={wtp_grid[-1]}, expected AB highest; got {row_high.to_dict()}"
    )


def test_compute_nmb_shape_and_definition(psa_result):
    """NMB tensor has shape (n_sim, n_strat, n_wtp) and matches
    ``wtp*q - c`` for a hand-picked WTP."""
    wtp = np.array([0.0, 100_000.0])
    nmb = compute_nmb(psa_result, wtp)
    assert nmb.shape == (psa_result.n_sim, len(psa_result.strategy_names), 2)
    # At wtp = 0, NMB == -cost
    for j, name in enumerate(psa_result.strategy_names):
        np.testing.assert_allclose(
            nmb[:, j, 0],
            -psa_result.costs[name].to_numpy(),
            atol=1e-9,
        )


# ---------------------------------------------------------------------------
# 5. Frontier / dominance
# ---------------------------------------------------------------------------


def test_strategy_a_is_dominated_in_expectation(psa_result):
    """Mean cost(A) > mean cost(B) AND mean QALY(A) < mean QALY(B): A is
    strictly dominated by B in expectation, matching the deterministic
    Table 6 result."""
    mean_cost = psa_result.costs.mean()
    mean_qaly = psa_result.qalys.mean()
    assert mean_cost["Strategy A"] > mean_cost["Strategy B"]
    assert mean_qaly["Strategy A"] < mean_qaly["Strategy B"]


def test_expected_nmb_frontier_never_picks_strategy_a(psa_result, wtp_grid):
    """Because A is dominated in expectation by B (lower mean cost, higher
    mean QALYs), no WTP should pick A on the expected-NMB frontier."""
    frontier = expected_nmb_frontier(psa_result, wtp_grid)
    assert "Strategy A" not in set(frontier["best_strategy"].unique())


def test_incremental_vs_baseline_shape(psa_result):
    inc = incremental_vs_baseline(psa_result, baseline="Standard of care")
    assert set(inc["strategy"].unique()) == {"Strategy A", "Strategy B", "Strategy AB"}
    assert len(inc) == 3 * psa_result.n_sim
    assert {"draw", "strategy", "inc_cost", "inc_qaly"}.issubset(inc.columns)


# ---------------------------------------------------------------------------
# 6. Plotting smoke test — files get written, axes don't blow up
# ---------------------------------------------------------------------------


def test_plots_render_to_disk(psa_result, wtp_grid, tmp_path):
    """The plotting module produces non-empty PNG files for all three
    figures. Content correctness is exercised by the CEAC / frontier
    statistical tests above; this test just guards the matplotlib path."""
    from opencea.plots import plot_ce_frontier, plot_ce_plane, plot_ceac

    p1 = plot_ce_plane(psa_result, tmp_path / "ce_plane.png")
    p2 = plot_ceac(psa_result, tmp_path / "ceac.png", wtp_grid=wtp_grid)
    p3 = plot_ce_frontier(psa_result, tmp_path / "frontier.png", wtp_grid=wtp_grid)
    for p in (p1, p2, p3):
        assert p.exists() and p.stat().st_size > 0
