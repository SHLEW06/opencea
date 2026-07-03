"""OpenCEA — open, reproducible health-economic decision modeling.

Surface: Pydantic model spec, cohort simulation engine, basic CEA,
probabilistic sensitivity analysis (PSA) with CEAC + frontier, one-way
deterministic sensitivity analysis (DSA) with tornado, and a matplotlib
plotting layer for all figures.
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
from .sensitivity import (
    DSAResult,
    ParameterSweep,
    compute_parameter_ranges,
    run_dsa,
)
from .empagliflozin import (
    EMPA,
    EMPA_DSA_RANGE_OVERRIDES,
    EMPA_PSA_SPECS,
    SOC,
    WaningSpec,
    breakeven_drug_price,
    build_empagliflozin_t2d,
    case_results_for_cea,
    dsa_evaluator,
    evaluate_empagliflozin_case,
    evaluate_scenario,
    run_empa_psa,
    run_empa_psa_scenario,
    scenario_icer,
)
from .engine import evaluate_sequence, simulate_trace_sequence

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
    "DSAResult",
    "ParameterSweep",
    "compute_parameter_ranges",
    "run_dsa",
    "SOC",
    "EMPA",
    "EMPA_PSA_SPECS",
    "EMPA_DSA_RANGE_OVERRIDES",
    "WaningSpec",
    "breakeven_drug_price",
    "build_empagliflozin_t2d",
    "case_results_for_cea",
    "dsa_evaluator",
    "evaluate_empagliflozin_case",
    "evaluate_scenario",
    "run_empa_psa",
    "run_empa_psa_scenario",
    "scenario_icer",
    "evaluate_sequence",
    "simulate_trace_sequence",
]
