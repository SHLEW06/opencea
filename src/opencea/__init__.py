"""OpenCEA — open, reproducible health-economic decision modeling.

Stage 1 surface: Pydantic model spec, cohort simulation engine, and basic
cost-effectiveness analysis. PSA, DSA, plots, app, and reporting come later.
"""
from .model import CohortModel, Strategy
from .engine import simulate_trace, evaluate_strategy, run_model, gen_wcc
from .cea import cea_table, icer, nmb
from .builders import build_darth_sick_sicker, rate_to_prob

__version__ = "0.1.0"

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
]
