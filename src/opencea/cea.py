"""Basic cost-effectiveness outputs: totals, incrementals, ICER, NMB, dominance.

Dominance handling matches the standard CEA frontier algorithm:

* **Strict dominance** — strategy ``i`` is strictly dominated iff some
  other strategy ``j`` achieves ``cost_j <= cost_i`` and ``qaly_j >=
  qaly_i`` with at least one strict inequality.
* **Extended (weak) dominance** — after removing strictly dominated
  strategies, sort the survivors by cost and iteratively drop any
  strategy whose ICER vs the previous on-frontier point exceeds the
  ICER of the *next* on-frontier point (equivalently, any survivor not
  on the convex lower-right hull of the cost / QALY plane). Repeat until
  the per-step ICERs are monotonically increasing along the frontier.

ICERs reported in :func:`cea_table` are **on-frontier** ICERs — each
computed against the previous strategy on the efficient frontier, not
against the raw cost-sorted neighbor — so strictly and extended-dominated
rows show NaN for ``inc_cost`` / ``inc_qaly`` / ``icer``.

CEAC and PSA-derived outputs live in :mod:`opencea.psa`.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


# Status labels used by the cea_table ``status`` column.
STATUS_FRONTIER = "frontier"
STATUS_STRICT = "strict"
STATUS_EXTENDED = "extended"


def icer(cost0: float, qaly0: float, cost1: float, qaly1: float) -> float:
    """Incremental cost-effectiveness ratio of strategy 1 vs strategy 0.

    Returns ``nan`` if the incremental QALYs are zero (the ratio is
    undefined and the comparison should be made via NMB instead).
    """
    dq = qaly1 - qaly0
    if dq == 0:
        return math.nan
    return (cost1 - cost0) / dq


def nmb(cost: float, qaly: float, wtp: float) -> float:
    """Net monetary benefit: ``wtp * qaly - cost``."""
    return wtp * qaly - cost


def _flag_strictly_dominated(costs: np.ndarray, qalys: np.ndarray) -> np.ndarray:
    """Return a boolean mask marking strictly-dominated strategies.

    Strategy ``i`` is strictly dominated iff some other strategy ``j``
    achieves ``cost_j <= cost_i`` and ``qaly_j >= qaly_i`` with at least
    one strict inequality.
    """
    n = costs.size
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (
                costs[j] <= costs[i]
                and qalys[j] >= qalys[i]
                and (costs[j] < costs[i] or qalys[j] > qalys[i])
            ):
                dominated[i] = True
                break
    return dominated


def _prune_extended_dominated(
    candidate_idx: List[int], costs: np.ndarray, qalys: np.ndarray
) -> tuple[List[int], List[int]]:
    """Iteratively prune extended-dominated strategies from a cost-sorted list.

    Given ``candidate_idx`` (already sorted by cost and free of strictly
    dominated points), walk left-to-right and drop any interior point
    whose ICER vs the previous on-frontier point exceeds the ICER of the
    next on-frontier point. One sweep can expose new violations, so the
    procedure is repeated to a fixed point.

    Returns ``(frontier_idx, extended_dominated_idx)``.
    """
    extended: List[int] = []
    work = list(candidate_idx)
    while True:
        # Build sequential ICERs along the current candidate frontier.
        # A point j is extended-dominated when ICER(j vs prev) > ICER(next vs j).
        changed = False
        i = 1
        while i < len(work) - 1:
            prev = work[i - 1]
            cur = work[i]
            nxt = work[i + 1]
            dq_cur = qalys[cur] - qalys[prev]
            dq_nxt = qalys[nxt] - qalys[cur]
            # On a cost-sorted, non-strictly-dominated frontier QALYs are
            # strictly increasing; guard against equal-QALY ties anyway so
            # we never divide by zero.
            if dq_cur <= 0 or dq_nxt <= 0:
                i += 1
                continue
            icer_cur = (costs[cur] - costs[prev]) / dq_cur
            icer_nxt = (costs[nxt] - costs[cur]) / dq_nxt
            if icer_cur > icer_nxt:
                extended.append(cur)
                del work[i]
                changed = True
                # Restart the scan from the affected segment.
                i = max(1, i - 1)
            else:
                i += 1
        if not changed:
            break
    return work, extended


def _classify_dominance(costs: np.ndarray, qalys: np.ndarray) -> np.ndarray:
    """Classify each strategy as ``frontier``, ``strict``, or ``extended``.

    Operates on the order of ``costs`` / ``qalys`` as given — call sites
    can sort first if they want frontier-order output.
    """
    n = costs.size
    status = np.array([STATUS_FRONTIER] * n, dtype=object)

    strict_mask = _flag_strictly_dominated(costs, qalys)
    for i in np.where(strict_mask)[0]:
        status[i] = STATUS_STRICT

    survivors = [i for i in range(n) if not strict_mask[i]]
    # Cost ascending; QALY ascending on ties so the "next more effective"
    # neighbor is well-defined when costs coincide.
    survivors.sort(key=lambda i: (costs[i], qalys[i]))

    _, extended_idx = _prune_extended_dominated(survivors, costs, qalys)
    for i in extended_idx:
        status[i] = STATUS_EXTENDED
    return status


def cea_table(
    results: Iterable[Dict[str, object]],
    wtp: Optional[float] = None,
) -> pd.DataFrame:
    """Build a cost-effectiveness table from a list of strategy results.

    Each result must expose ``name``, ``total_cost``, and ``total_qaly``
    (as produced by :func:`opencea.engine.evaluate_strategy`).

    The returned DataFrame is sorted by total cost ascending and always
    contains:

    - ``strategy``, ``total_cost``, ``total_qaly``
    - ``status`` — one of ``"frontier"``, ``"strict"``, ``"extended"``.
    - ``dominated`` — boolean back-compat shortcut, ``status != "frontier"``.
    - ``inc_cost``, ``inc_qaly``, ``icer`` — computed by walking the
      **efficient frontier** (after pruning both strict and extended
      dominance); dominated rows show NaN.

    If ``wtp`` is provided, a ``nmb`` column is added (computed for all
    strategies, dominated or not).
    """
    rows: List[Dict[str, object]] = []
    for r in results:
        rows.append(
            {
                "strategy": r["name"],
                "total_cost": float(r["total_cost"]),
                "total_qaly": float(r["total_qaly"]),
            }
        )

    df = pd.DataFrame(rows).sort_values("total_cost").reset_index(drop=True)

    costs = df["total_cost"].to_numpy()
    qalys = df["total_qaly"].to_numpy()

    status = _classify_dominance(costs, qalys)
    df["status"] = status
    df["dominated"] = status != STATUS_FRONTIER

    # Incremental analysis along the efficient frontier (in cost order).
    inc_cost = np.full(len(df), np.nan)
    inc_qaly = np.full(len(df), np.nan)
    icer_col = np.full(len(df), np.nan)

    frontier_idx = [i for i in range(len(df)) if status[i] == STATUS_FRONTIER]
    for k, i in enumerate(frontier_idx):
        if k == 0:
            continue
        prev = frontier_idx[k - 1]
        dc = float(df.at[i, "total_cost"] - df.at[prev, "total_cost"])
        dq = float(df.at[i, "total_qaly"] - df.at[prev, "total_qaly"])
        inc_cost[i] = dc
        inc_qaly[i] = dq
        icer_col[i] = dc / dq if dq != 0 else np.nan

    df["inc_cost"] = inc_cost
    df["inc_qaly"] = inc_qaly
    df["icer"] = icer_col

    if wtp is not None:
        df["nmb"] = wtp * df["total_qaly"] - df["total_cost"]

    return df
