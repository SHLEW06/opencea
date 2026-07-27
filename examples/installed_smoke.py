"""Run a self-contained calculation with an installed OpenCEA package."""

from __future__ import annotations

import opencea
from opencea import CohortModel, Strategy, run_model


def main() -> None:
    strategy = Strategy(
        name="care",
        transition_matrix=[[0.8, 0.2], [0.0, 1.0]],
        state_costs=[100.0, 0.0],
        state_utilities=[1.0, 0.0],
    )
    model = CohortModel(
        states=["well", "dead"],
        strategies=[strategy],
        initial_distribution=[1.0, 0.0],
        time_horizon=2,
        discount_rate_costs=0.0,
        discount_rate_qalys=0.0,
        wcc_method="none",
    )
    result = run_model(model)[0]

    assert abs(result["total_cost"] - 244.0) < 1e-9
    assert abs(result["total_qaly"] - 2.44) < 1e-12
    print(
        f"OpenCEA {opencea.__version__}: "
        f"total_cost={result['total_cost']:.2f}, "
        f"total_qaly={result['total_qaly']:.2f}"
    )


if __name__ == "__main__":
    main()
