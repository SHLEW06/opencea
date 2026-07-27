"""Matplotlib figures for PSA and DSA outputs.

All figures use the headless ``Agg`` backend and are saved to a file path
chosen by the caller. The CEAC is the project's main PSA figure, with one
probability curve per strategy across the WTP grid. The tornado is the
main DSA figure.

The PSA helpers take a :class:`opencea.psa.PSAResult` and call the
analytic helpers in :mod:`opencea.psa` so plots and tests share one
source of truth for NMB / CEAC values. The tornado takes a
:class:`opencea.sensitivity.DSAResult` directly.
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
from .sensitivity import DSAResult

PathLike = Union[str, Path]


def _format_money(x: float, _pos: int = 0) -> str:
    return f"${x / 1000:,.0f}k" if abs(x) >= 1000 else f"${x:,.0f}"


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
        if handle is not None:
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
    """Cost-effectiveness acceptability curve.

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


# ---------------------------------------------------------------------------
# DSA tornado
# ---------------------------------------------------------------------------


def _fmt_param_value(v: float) -> str:
    """Compact parameter-value label used at the bar ends."""
    av = abs(v)
    if av >= 100:
        return f"{v:,.0f}"
    if av >= 1:
        return f"{v:.3g}"
    return f"{v:.4g}"


def plot_tornado(
    dsa_result: DSAResult,
    out_path: PathLike,
    top_n: Optional[int] = None,
    dpi: int = 150,
) -> Path:
    """Classic CEA tornado: horizontal bars sorted by swing, centered on base.

    Each parameter gets one bar spanning its low-parameter outcome to its
    high-parameter outcome. The bar is split at the base-case outcome and
    coloured by which side of the base each end falls on. Bars are
    annotated with the parameter's low / high input values.

    Parameters
    ----------
    dsa_result
        Output of :func:`opencea.sensitivity.run_dsa`.
    out_path
        Destination path for the PNG.
    top_n
        Optional cap on the number of parameters shown (largest swings).
    """
    sweeps = list(dsa_result.sweeps)
    if top_n is not None:
        sweeps = sweeps[:top_n]
    # matplotlib horizontal bars draw from the bottom up; reverse so the
    # largest-swing parameter sits at the top of the plot.
    sweeps_plot = list(reversed(sweeps))

    base = dsa_result.base_outcome
    n = len(sweeps_plot)

    fig, ax = plt.subplots(figsize=(11, max(3.5, 0.45 * n + 2.5)))

    color_low = "#4C78A8"  # outcome at the parameter's LOW value
    color_high = "#E15759"  # outcome at the parameter's HIGH value
    bar_h = 0.7

    for i, sw in enumerate(sweeps_plot):
        # Draw both bars from base out to each endpoint. When both ends sit
        # on the same side of base (rare, only when base coincides with a
        # bound) the shorter bar is fully contained inside the longer one;
        # a thin black edge keeps both visible.
        ax.barh(
            i,
            sw.low_outcome - base,
            left=base,
            height=bar_h,
            color=color_low,
            edgecolor="black",
            linewidth=0.4,
            label="parameter at low value" if i == n - 1 else None,
        )
        ax.barh(
            i,
            sw.high_outcome - base,
            left=base,
            height=bar_h,
            color=color_high,
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
            label="parameter at high value" if i == n - 1 else None,
        )

    # Annotate each bar end with the corresponding parameter value, outside
    # the bar to keep the chart readable.
    xmin = min(min(s.low_outcome, s.high_outcome) for s in sweeps_plot)
    xmax = max(max(s.low_outcome, s.high_outcome) for s in sweeps_plot)
    xmin = min(xmin, base)
    xmax = max(xmax, base)
    pad = 0.04 * (xmax - xmin if xmax > xmin else max(abs(base), 1.0))

    for i, sw in enumerate(sweeps_plot):
        left_end = min(sw.low_outcome, sw.high_outcome)
        right_end = max(sw.low_outcome, sw.high_outcome)
        left_val = sw.low_value if sw.low_outcome <= sw.high_outcome else sw.high_value
        right_val = sw.high_value if sw.high_outcome >= sw.low_outcome else sw.low_value
        ax.text(
            left_end - pad * 0.3,
            i,
            _fmt_param_value(left_val),
            ha="right",
            va="center",
            fontsize=8,
            color="#333",
        )
        ax.text(
            right_end + pad * 0.3,
            i,
            _fmt_param_value(right_val),
            ha="left",
            va="center",
            fontsize=8,
            color="#333",
        )

    ax.axvline(base, color="black", lw=1.2)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels([sw.parameter for sw in sweeps_plot])
    ax.set_xlim(xmin - 2.5 * pad, xmax + 2.5 * pad)

    ax.set_xlabel(
        f"Incremental NMB ({dsa_result.comparator} vs {dsa_result.baseline}) "
        f"at WTP = ${dsa_result.wtp:,.0f} / QALY (USD)"
    )
    ax.set_title("One-way deterministic sensitivity analysis")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(_format_money))
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()

    out = Path(out_path)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out
