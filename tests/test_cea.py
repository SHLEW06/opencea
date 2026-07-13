"""Hand-checked tests for strict and extended dominance in :mod:`opencea.cea`.

Every cost / QALY pair below is constructed so the frontier is computable
by hand; the tests assert both the dominance classification and the
exact ICER reported on the frontier.

Notation: the standard CEA frontier algorithm classifies each strategy as

  * ``"strict"`` — strictly dominated (more costly *and* fewer QALYs than
    some other strategy).
  * ``"extended"`` — weakly dominated; not on the convex lower-right hull
    of the cost / QALY plane after strict dominants are removed.
  * ``"frontier"`` — on the efficient frontier; its incremental cost,
    incremental QALY, and ICER columns are computed against the previous
    *on-frontier* strategy (not the raw cost-sorted neighbor).
"""

from __future__ import annotations

import numpy as np
import pytest

from opencea.cea import (
    STATUS_EXTENDED,
    STATUS_FRONTIER,
    STATUS_STRICT,
    cea_table,
    icer,
)


def _make_results(rows):
    """Helper: turn a list of (name, cost, qaly) tuples into the engine's
    result shape consumed by ``cea_table``."""
    return [
        {"name": n, "total_cost": float(c), "total_qaly": float(q)} for n, c, q in rows
    ]


def _row(df, name):
    return df.loc[df["strategy"] == name].iloc[0]


# ---------------------------------------------------------------------------
# 1. Hand-checked extended-dominance case
# ---------------------------------------------------------------------------


def test_extended_dominance_three_strategies_hand_checked():
    """Canonical pattern:

        A   $10,000   1.0 QALY
        B   $30,000   1.5 QALY   <- extended-dominated by A & C
        C   $50,000   2.5 QALY

    ICER B vs A = $40,000 / QALY; ICER C vs B = $20,000 / QALY.
    Because ICER(C vs B) < ICER(B vs A), B is extended-dominated:
    a 50/50 mix of A and C would deliver the same QALYs as B at lower
    total cost. The efficient frontier is {A, C} with
    ICER(C vs A) = $40,000 / $1.5 = **$26,666.67 / QALY** (exact).

    B is NOT strictly dominated — C costs more than B, A delivers fewer
    QALYs than B — so this isolates extended dominance from strict.
    """
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 30_000, 1.5),
                ("C", 50_000, 2.5),
            ]
        )
    )

    # Classification
    assert _row(df, "A")["status"] == STATUS_FRONTIER
    assert _row(df, "B")["status"] == STATUS_EXTENDED
    assert _row(df, "C")["status"] == STATUS_FRONTIER
    assert bool(_row(df, "B")["dominated"]) is True

    # B is NOT strictly dominated by either A or C (would-be dominators
    # would need both lower-or-equal cost AND higher-or-equal QALY).
    assert not (10_000 <= 30_000 and 1.0 >= 1.5)  # A doesn't dominate B
    assert not (50_000 <= 30_000 and 2.5 >= 1.5)  # C doesn't dominate B

    # Frontier ICERs: B has NaN; C is computed vs A (the previous on-frontier
    # point), not vs B.
    assert np.isnan(_row(df, "B")["icer"])
    assert _row(df, "C")["icer"] == pytest.approx(40_000.0 / 1.5)
    assert _row(df, "C")["inc_cost"] == pytest.approx(40_000.0)
    assert _row(df, "C")["inc_qaly"] == pytest.approx(1.5)


def test_iterative_pruning_with_four_strategies():
    """One sweep can expose new extended-dominance violations; the
    algorithm must iterate to a fixed point.

        A   $10,000   1.0
        B   $20,000   1.2   <- ICER vs A = $50k / QALY
        C   $35,000   1.5   <- ICER vs B = $50k / QALY; vs A = $50k / QALY
        D   $50,000   3.0   <- ICER vs C = $10k / QALY; vs A = $20k / QALY

    D's much higher cost-effectiveness (ICER $20k vs A directly) makes
    both B and C extended-dominated. The frontier is {A, D} with
    ICER(D vs A) = $40,000 / 2.0 = **$20,000 / QALY** exactly.
    """
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 20_000, 1.2),
                ("C", 35_000, 1.5),
                ("D", 50_000, 3.0),
            ]
        )
    )
    assert _row(df, "A")["status"] == STATUS_FRONTIER
    assert _row(df, "B")["status"] == STATUS_EXTENDED
    assert _row(df, "C")["status"] == STATUS_EXTENDED
    assert _row(df, "D")["status"] == STATUS_FRONTIER
    assert _row(df, "D")["icer"] == pytest.approx(20_000.0)
    assert _row(df, "D")["inc_qaly"] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# 2. Monotonic-ICER case — every strategy on the frontier
# ---------------------------------------------------------------------------


def test_monotonic_no_dominance_all_on_frontier():
    """ICERs strictly increasing => no dominance, every strategy on frontier.

    A   $10,000   1.0
    B   $25,000   1.5   ICER vs A = $30,000 / QALY
    C   $50,000   2.0   ICER vs B = $50,000 / QALY
    """
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 25_000, 1.5),
                ("C", 50_000, 2.0),
            ]
        )
    )
    assert (df["status"] == STATUS_FRONTIER).all()
    assert not df["dominated"].any()
    assert _row(df, "B")["icer"] == pytest.approx(30_000.0)
    assert _row(df, "C")["icer"] == pytest.approx(50_000.0)


# ---------------------------------------------------------------------------
# 3. Strict-dominance-only regression
# ---------------------------------------------------------------------------


def test_strict_dominance_only_regression():
    """No extended dominance, but B is strictly dominated by A
    (more costly AND fewer QALYs).

        A   $10,000   1.0
        B   $20,000   0.8   <- strictly dominated by A
        C   $30,000   2.0

    Frontier {A, C}; ICER(C vs A) = $20,000 / 1.0 = **$20,000 / QALY**.
    """
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 20_000, 0.8),
                ("C", 30_000, 2.0),
            ]
        )
    )
    assert _row(df, "A")["status"] == STATUS_FRONTIER
    assert _row(df, "B")["status"] == STATUS_STRICT
    assert _row(df, "C")["status"] == STATUS_FRONTIER
    assert _row(df, "C")["icer"] == pytest.approx(20_000.0)


def test_strict_and_extended_mixed():
    """Mix: B strictly dominated by A, C extended-dominated, D on frontier.

    A   $10,000   1.0
    B   $15,000   0.9   <- strictly dominated by A
    C   $30,000   1.5   <- extended-dominated (ICER vs A 40k > ICER D vs C 30k)
    D   $60,000   2.5   <- frontier; ICER vs A = $50k / $1.5 = 33,333.33
    """
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 15_000, 0.9),
                ("C", 30_000, 1.5),
                ("D", 60_000, 2.5),
            ]
        )
    )
    assert _row(df, "B")["status"] == STATUS_STRICT
    assert _row(df, "C")["status"] == STATUS_EXTENDED
    assert _row(df, "A")["status"] == STATUS_FRONTIER
    assert _row(df, "D")["status"] == STATUS_FRONTIER
    assert _row(df, "D")["icer"] == pytest.approx(50_000.0 / 1.5)


# ---------------------------------------------------------------------------
# 4. Post-algorithm invariant: frontier ICERs strictly increase
# ---------------------------------------------------------------------------


def test_frontier_icers_strictly_increase_after_pruning():
    """The defining property of an efficient frontier: after pruning both
    strict and extended dominance, ICERs along the frontier (in cost
    order) are strictly increasing."""
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 20_000, 1.2),  # extended later
                ("C", 25_000, 1.3),  # extended later
                ("D", 30_000, 1.4),  # extended later
                ("E", 50_000, 3.0),  # frontier
                ("F", 80_000, 4.0),  # frontier
                ("G", 90_000, 3.9),  # strictly dominated by F
            ]
        )
    )
    frontier_rows = df[df["status"] == STATUS_FRONTIER].sort_values("total_cost")
    icers = frontier_rows["icer"].dropna().to_numpy()
    assert icers.size >= 1
    diffs = np.diff(icers)
    assert (diffs > 0).all(), (
        f"Frontier ICERs must be strictly increasing; got {icers.tolist()}"
    )


def test_two_strategies_only():
    """Trivial 2-strategy case — nothing can be extended-dominated."""
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 20_000, 1.5),
            ]
        )
    )
    assert list(df["status"]) == [STATUS_FRONTIER, STATUS_FRONTIER]
    assert _row(df, "B")["icer"] == pytest.approx(icer(10_000, 1.0, 20_000, 1.5))


def test_single_strategy():
    """No comparator — no incrementals."""
    df = cea_table(_make_results([("A", 10_000, 1.0)]))
    assert df.iloc[0]["status"] == STATUS_FRONTIER
    assert np.isnan(df.iloc[0]["icer"])


# ---------------------------------------------------------------------------
# 5. Back-compat: ``dominated`` boolean = (status != frontier)
# ---------------------------------------------------------------------------


def test_dominated_column_back_compat():
    df = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 30_000, 1.5),
                ("C", 50_000, 2.5),
            ]
        )
    )
    assert bool(_row(df, "A")["dominated"]) is False
    assert bool(_row(df, "B")["dominated"]) is True  # extended
    assert bool(_row(df, "C")["dominated"]) is False
    # The boolean still flips on strict dominance too.
    df2 = cea_table(
        _make_results(
            [
                ("A", 10_000, 1.0),
                ("B", 20_000, 0.8),
                ("C", 30_000, 2.0),
            ]
        )
    )
    assert bool(_row(df2, "B")["dominated"]) is True  # strict


# ---------------------------------------------------------------------------
# 6. DARTH classification consistency
# ---------------------------------------------------------------------------


def test_darth_strategy_a_is_strictly_dominated_not_extended():
    """In the DARTH Sick-Sicker manuscript Strategy A costs MORE than B
    AND yields FEWER QALYs than B — that is strict dominance, not
    extended. Confirms the new classifier preserves the manuscript's
    intent."""
    df = cea_table(
        _make_results(
            [
                ("Standard of care", 151_580.0, 20.711),
                ("Strategy A", 284_805.0, 21.499),
                ("Strategy B", 259_100.0, 22.184),
                ("Strategy AB", 378_875.0, 23.137),
            ]
        )
    )
    assert _row(df, "Strategy A")["status"] == STATUS_STRICT
    # Manuscript ICERs (to the dollar) — must come out of the frontier walk.
    assert _row(df, "Strategy B")["icer"] == pytest.approx(72_988.0, abs=10.0)
    assert _row(df, "Strategy AB")["icer"] == pytest.approx(125_764.0, abs=200.0)
