"""OpenCEA — open, reproducible health-economic decision modeling.

Surface: Pydantic model spec, cohort simulation engine, basic CEA,
probabilistic sensitivity analysis (PSA) with CEAC + frontier, and a
matplotlib plotting layer for the CE plane, CEAC, and frontier figures.
"""
from .model import CohortModel, Strategy
from .engine import simulate_trace, evaluate_strategy, run_model, gen_wcc
from .cea import cea_table, icer, nmb
from .builders import build_darth_sick_sicker, rate_to_prob
from .psa import (
    DistSpec,
    PSAResult,
    PSA_PARAM_SPECS,
    compute_ceac,
    compute_nmb,
    default_wtp_grid,
    expected_nmb_frontier,
    incremental_vs_baseline,
    run_psa,
    sample_psa_params,
)

__version__ = "0.2.0"

__all__ = [
    "CohortModel",
    "Strategy",
    "simulate_trace",
    "evaluate_strategy",
    "run_model",
    "gen_wcc",
    "cea_table",
    "icer",
    "nmb",
    "build_darth_sick_sicker",
    "rate_to_prob",
    "DistSpec",
    "PSAResult",
    "PSA_PARAM_SPECS",
    "compute_ceac",
    "compute_nmb",
    "default_wtp_grid",
    "expected_nmb_frontier",
    "incremental_vs_baseline",
    "run_psa",
    "sample_psa_params",
]
