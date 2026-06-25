"""Basic cost-effectiveness outputs: totals, incrementals, ICER, NMB, dominance.

Outputs are intentionally minimal in stage 1: a tabular CEA summary with
strict dominance flagging plus a couple of standalone helpers (``icer``,
``nmb``). Extended dominance, frontier construction, and PSA-derived
CEACs belong in later stages.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


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


def cea_table(
    results: Iterable[Dict[str, object]],
    wtp: Optional[float] = None,
) -> pd.DataFrame:
    """Build a cost-effectiveness table from a list of strategy results.

    Each result must expose ``name``, ``total_cost``, and ``total_qaly``
    (as produced by :func:`opencea.engine.evaluate_strategy`).

    The returned DataFrame is sorted by total cost ascending. It always
    contains:

    - ``strategy``, ``total_cost``, ``total_qaly``
    - ``dominated`` (strict dominance flag)
    - ``inc_cost``, ``inc_qaly``, ``icer`` — computed by walking the
      cost-sorted non-dominated frontier, so dominated rows show NaN.

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
    df["dominated"] = _flag_strictly_dominated(costs, qalys)

    # Incremental analysis on the non-dominated, cost-sorted frontier.
    inc_cost = np.full(len(df), np.nan)
    inc_qaly = np.full(len(df), np.nan)
    icer_col = np.full(len(df), np.nan)

    frontier_idx = [i for i in range(len(df)) if not df.at[i, "dominated"]]
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
