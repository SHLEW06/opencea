"""Regenerate the quantitative assets used in the OpenCEA README."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from opencea.empagliflozin import (
    EMPA,
    EMPA_DSA_RANGE_OVERRIDES,
    EMPA_PSA_SPECS,
    SOC,
    WaningSpec,
    breakeven_drug_price,
    dsa_evaluator,
    run_empa_psa,
    run_empa_psa_scenario,
    scenario_icer,
)
from opencea.plots import plot_ceac, plot_tornado
from opencea.psa import PSAResult, compute_ceac, default_wtp_grid
from opencea.sensitivity import run_dsa

ROOT = Path(__file__).resolve().parents[1]
PARAMS = ROOT / "examples" / "empagliflozin_t2d.yaml"
ASSET_DIR = ROOT / "docs" / "assets"

N_SIM = 10_000
SEED = 20260626
WTP = 100_000.0
NET_PRICE = 4_500.0
WANING = WaningSpec(start_year=3.0, end_year=10.0)
TORNADO_PARAMETER_LABELS = {
    "hr_death": "Mortality hazard ratio",
    "c_drug": "Annual drug cost",
    "hr_event": "Major CV event hazard ratio",
    "u_EF": "Event-free utility",
    "r_PE_D": "Post-event mortality rate",
    "c_EF": "Event-free annual cost",
    "c_PE": "Post-event annual cost",
    "u_PE": "Post-event utility",
    "c_acute_PE": "Acute event cost",
}


@dataclass(frozen=True)
class ScenarioResult:
    """One deterministic and probabilistic case-study result."""

    name: str
    label: str
    annual_drug_price: Optional[float]
    waning_start_year: Optional[float]
    waning_end_year: Optional[float]
    icer_per_qaly: float
    probability_cost_effective: float


def _scenario_psa(name: str) -> PSAResult:
    if name == "wac_sustained":
        return run_empa_psa(PARAMS, n_sim=N_SIM, seed=SEED)
    if name == "net_sustained":
        return run_empa_psa_scenario(
            PARAMS,
            n_sim=N_SIM,
            seed=SEED,
            drug_price=NET_PRICE,
        )
    if name == "wac_waning":
        return run_empa_psa_scenario(
            PARAMS,
            n_sim=N_SIM,
            seed=SEED,
            waning=WANING,
        )
    if name == "net_waning":
        return run_empa_psa_scenario(
            PARAMS,
            n_sim=N_SIM,
            seed=SEED,
            drug_price=NET_PRICE,
            waning=WANING,
        )
    raise ValueError(f"unknown scenario: {name}")


def _probability_cost_effective(psa: PSAResult) -> float:
    ceac = compute_ceac(psa, np.array([WTP]))
    return float(ceac.loc[WTP, EMPA])


def compute_results() -> tuple[list[ScenarioResult], dict[str, float], PSAResult]:
    """Run the four public scenarios with the README's fixed seed."""

    definitions = [
        ("wac_sustained", "WAC, sustained", None, None),
        ("net_sustained", "$4,500, sustained", NET_PRICE, None),
        ("wac_waning", "WAC, waning", None, WANING),
        ("net_waning", "$4,500, waning", NET_PRICE, WANING),
    ]

    scenarios: list[ScenarioResult] = []
    base_psa: Optional[PSAResult] = None
    for name, label, price, waning in definitions:
        psa = _scenario_psa(name)
        if name == "wac_sustained":
            base_psa = psa
        scenarios.append(
            ScenarioResult(
                name=name,
                label=label,
                annual_drug_price=price,
                waning_start_year=None if waning is None else waning.start_year,
                waning_end_year=None if waning is None else waning.end_year,
                icer_per_qaly=scenario_icer(
                    PARAMS,
                    drug_price=price,
                    waning=waning,
                ),
                probability_cost_effective=_probability_cost_effective(psa),
            )
        )

    if base_psa is None:  # pragma: no cover
        raise RuntimeError("base-case PSA was not created")

    break_even = {
        "sustained_annual_drug_price": breakeven_drug_price(
            PARAMS,
            target_icer=WTP,
        ),
        "waning_annual_drug_price": breakeven_drug_price(
            PARAMS,
            target_icer=WTP,
            waning=WANING,
        ),
    }
    return scenarios, break_even, base_psa


def plot_scenario_icers(
    scenarios: list[ScenarioResult],
    out_path: Path,
) -> Path:
    """Plot deterministic ICERs against the stated decision threshold."""

    values = np.array([scenario.icer_per_qaly for scenario in scenarios])
    labels = [scenario.label.replace(", ", ",\n") for scenario in scenarios]
    colors = ["#2F6690" if value <= WTP else "#C65D21" for value in values]

    fig, ax = plt.subplots(figsize=(9, 5.2))
    bars = ax.bar(labels, values, color=colors, width=0.66)
    ax.axhline(
        WTP,
        color="#222222",
        linewidth=1.2,
        linestyle="--",
        label="$100k/QALY threshold",
    )
    ax.bar_label(
        bars,
        labels=[f"${value / 1000:,.1f}k" for value in values],
        padding=4,
        fontsize=10,
    )
    ax.set_ylabel("Incremental cost-effectiveness ratio (USD/QALY)")
    ax.set_title("Empagliflozin scenario ICERs")
    ax.set_ylim(0, max(values) * 1.17)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda value, _position: f"${value / 1000:,.0f}k")
    )
    ax.grid(axis="y", alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_result_summary(
    scenarios: list[ScenarioResult],
    break_even: dict[str, float],
    out_path: Path,
) -> Path:
    """Write the exact inputs and outputs behind the README table."""

    payload = {
        "model": "examples/empagliflozin_t2d.yaml",
        "seed": SEED,
        "psa_draws": N_SIM,
        "willingness_to_pay_per_qaly": WTP,
        "scenarios": [asdict(scenario) for scenario in scenarios],
        "break_even": break_even,
    }
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    return out_path


def main() -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    scenarios, break_even, base_psa = compute_results()

    plot_scenario_icers(
        scenarios,
        ASSET_DIR / "empagliflozin-scenario-icers.png",
    )
    plot_ceac(
        base_psa,
        ASSET_DIR / "empagliflozin-ceac.png",
        wtp_grid=default_wtp_grid(),
        dpi=180,
    )

    dsa = run_dsa(
        base_params=PARAMS,
        wtp=WTP,
        baseline=SOC,
        sweep_params=list(EMPA_PSA_SPECS),
        evaluator=dsa_evaluator,
        param_specs=EMPA_PSA_SPECS,
        param_ranges=EMPA_DSA_RANGE_OVERRIDES,
    )
    plot_tornado(
        dsa,
        ASSET_DIR / "empagliflozin-tornado.png",
        dpi=180,
        parameter_labels=TORNADO_PARAMETER_LABELS,
    )
    write_result_summary(
        scenarios,
        break_even,
        ASSET_DIR / "empagliflozin-results.json",
    )

    for scenario in scenarios:
        print(
            f"{scenario.name}: ICER=${scenario.icer_per_qaly:,.0f}/QALY, "
            f"P(CE)={scenario.probability_cost_effective:.4f}"
        )
    print(
        "break_even: "
        f"sustained=${break_even['sustained_annual_drug_price']:,.0f}/year, "
        f"waning=${break_even['waning_annual_drug_price']:,.0f}/year"
    )


if __name__ == "__main__":
    main()
