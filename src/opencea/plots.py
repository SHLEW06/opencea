"""Matplotlib figures for PSA outputs.

All figures use the headless ``Agg`` backend and are saved to a file path
chosen by the caller. The CEAC is the project's hero figure — one
probability curve per strategy across the WTP grid — and is intentionally
the cleanest of the three.

The functions take a :class:`opencea.psa.PSAResult` and a WTP grid,
internally calling the analytic helpers in :mod:`opencea.psa` so plots
and tests share one source of truth for NMB / CEAC values.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import matplotlib

matplotlib.use("Agg")  # headless: must be set before importing pyplot

import matplotlib.pyplot as plt
import numpy as np

from .psa import (
    PSAResult,
    compute_ceac,
    default_wtp_grid,
    expected_nmb_frontier,
    incremental_vs_baseline,
)


PathLike = Union[str, Path]


def _format_money(x: float, _pos: int = 0) -> str:
    return f"${x/1000:,.0f}k" if abs(x) >= 1000 else f"${x:,.0f}"


def plot_ce_plane(
    result: PSAResult,
    out_path: PathLike,
    baseline: str = "Standard of care",
    dpi: int = 150,
) -> Path:
    """Cost-effectiveness plane: incremental cost vs QALY scatter per draw.

    Each non-baseline strategy gets its own colour. A single semi-transparent
    point per draw keeps the cloud readable at n_sim = 10000.
    """
    inc = incremental_vs_baseline(result, baseline=baseline)

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = plt.get_cmap("tab10")
    for i, name in enumerate(sorted(inc["strategy"].unique())):
        sub = inc[inc["strategy"] == name]
        ax.scatter(
            sub["inc_qaly"],
            sub["inc_cost"],
            s=6,
            alpha=0.25,
            color=colors(i),
            label=name,
        )

    ax.axhline(0, color="grey", lw=0.6)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel(f"Incremental QALYs vs {baseline}")
    ax.set_ylabel(f"Incremental cost vs {baseline} (USD)")
    ax.set_title("Cost-effectiveness plane (PSA)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(_format_money))
    leg = ax.legend(loc="best", frameon=False, markerscale=2)
    for handle in leg.legend_handles:
        handle.set_alpha(1.0)
    fig.tight_layout()

    out = Path(out_path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_ceac(
    result: PSAResult,
    out_path: PathLike,
    wtp_grid: Optional[np.ndarray] = None,
    dpi: int = 150,
) -> Path:
    """Cost-effectiveness acceptability curve — the project hero figure.

    One probability curve per strategy across the WTP grid. Clean axes,
    labelled lines, no extraneous chrome.
    """
    grid = default_wtp_grid() if wtp_grid is None else np.asarray(wtp_grid, dtype=float)
    ceac = compute_ceac(result, grid)

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.get_cmap("tab10")
    for i, name in enumerate(result.strategy_names):
        ax.plot(grid, ceac[name].to_numpy(), lw=2.0, color=colors(i), label=name)

    ax.set_xlabel("Willingness to pay (USD / QALY)")
    ax.set_ylabel("Probability cost-effective")
    ax.set_title("Cost-effectiveness acceptability curve")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(grid.min(), grid.max())
    ax.xaxis.set_major_formatter(plt.FuncFormatter(_format_money))
    ax.grid(True, alpha=0.25)
    ax.legend(loc="center right", frameon=False)
    fig.tight_layout()

    out = Path(out_path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_ce_frontier(
    result: PSAResult,
    out_path: PathLike,
    wtp_grid: Optional[np.ndarray] = None,
    dpi: int = 150,
) -> Path:
    """Cost-effectiveness frontier from expected NMB across the WTP grid.

    Plots expected NMB curves and marks where the optimal strategy
    switches as WTP rises.
    """
    grid = default_wtp_grid() if wtp_grid is None else np.asarray(wtp_grid, dtype=float)
    frontier = expected_nmb_frontier(result, grid)

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.get_cmap("tab10")
    for i, name in enumerate(result.strategy_names):
        ax.plot(grid, frontier[name].to_numpy(), lw=2.0, color=colors(i), label=name)

    # Mark the switch points where the optimal strategy changes.
    best = frontier["best_strategy"].to_numpy()
    switches = np.where(best[:-1] != best[1:])[0]
    for k in switches:
        wtp_switch = grid[k + 1]
        ax.axvline(wtp_switch, color="black", lw=0.7, ls="--")
        ax.annotate(
            f"{best[k]} -> {best[k + 1]}",
            xy=(wtp_switch, ax.get_ylim()[1]),
            xytext=(5, -10),
            textcoords="offset points",
            fontsize=8,
            rotation=90,
            va="top",
        )

    ax.set_xlabel("Willingness to pay (USD / QALY)")
    ax.set_ylabel("Expected net monetary benefit (USD)")
    ax.set_title("Cost-effectiveness frontier (expected NMB)")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(_format_money))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(_format_money))
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()

    out = Path(out_path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out
