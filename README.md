# OpenCEA

[![CI](https://github.com/SHLEW06/OpenCEA/actions/workflows/ci.yml/badge.svg)](https://github.com/SHLEW06/OpenCEA/actions/workflows/ci.yml)

> Open, reproducible health-economic decision modeling in Python.

OpenCEA is an open-source Python toolkit for **cohort state-transition cost-effectiveness modeling**. It targets the gap left by R-only stacks (heemod, hesim, dampack) and proprietary tools (TreeAge) by providing a transparent, validated, vectorized engine in Python.

## Validation & tests

107 tests run on every push to `main` and on every pull request, across Python 3.10, 3.11, and 3.12:

- **Deterministic golden tests** pinned to the published DARTH Sick-Sicker manuscript (Table 5 totals to the cent; Table 6 ICERs to the dollar).
- **PSA structural tests** for sampler reproducibility, per-parameter Monte Carlo mean recovery, PSA-mean vs deterministic Table 5 within ~1-2%, and CEAC sanity at the WTP extremes.
- **DSA structural tests** for base-case consistency with the validated engine, bracketing, swing ordering, and the tornado plot.
- **Empagliflozin case-study tests** for structural integrity, evaluator-vs-engine consistency, ICER sanity, DSA tornado sanity, PSA reproducibility, and figure rendering.

## Status

### Stage 1 — validated deterministic core

All four DARTH Sick-Sicker strategies (SoC, A, B, AB) reproduce the published manuscript totals to the cent:

- `opencea.model` — Pydantic specification for a cohort state-transition model with input validation (square transition matrices, rows sum to 1, probabilities/utilities in [0, 1], non-negative costs).
- `opencea.engine` — NumPy cohort trace, discounting, and within-cycle correction (Simpson's 1/3 rule, half-cycle, or none); total discounted costs and QALYs.
- `opencea.cea` — Total cost / QALY, incremental cost / QALY, ICER, net monetary benefit, strict dominance flagging.
- `opencea.builders.build_darth_sick_sicker` — Constructs a `CohortModel` for the DARTH Sick-Sicker tutorial from a rate-based parameter spec using DARTH's rate-to-probability and competing-risks transition construction.
- `examples/sick_sicker.yaml` — Rate-based parameter file for the DARTH intro-tutorial model.
- `tests/test_darth_reference.py` — Golden-reference tests pinned to manuscript Table 5 / Table 6 (SoC $151,580 / 20.711, A $284,805 / 21.499, B $259,100 / 22.184, AB $378,875 / 23.137; ICER B vs SoC $72,988/QALY, ICER AB vs B $125,764/QALY).

### Stage 2 — probabilistic sensitivity analysis (PSA) and CEAC

- `opencea.psa` — PSA distribution spec, seeded sampler, and a vectorized Monte Carlo runner. Distributions are transcribed verbatim from the DARTH `generate_psa_params()` function (`R/Functions_cSTM_time_indep.R`, lines 278-309), with R's rate-parameterized `rgamma` converted to numpy's scale parameterization. The per-cycle trace step `trace[t+1] = trace[t] @ P` is batched via `einsum` over an `(n_sim, n_states, n_states)` transition tensor so the inner loop runs `T` matmuls of width `n_sim` instead of `n_sim * T` scalar matmuls.
- `compute_nmb`, `compute_ceac`, `expected_nmb_frontier`, `incremental_vs_baseline` — net monetary benefit per draw per strategy across a WTP grid (default 0 to 200,000), the cost-effectiveness acceptability curve, the expected-NMB-optimal strategy at each WTP, and the per-draw incremental cost / QALY vs SoC for the CE plane.
- `opencea.plots` — Matplotlib (headless `Agg` backend) figures: cost-effectiveness plane scatter, CEAC, and cost-effectiveness frontier. Output is saved to file paths chosen by the caller.
- `tests/test_psa.py` — Reproducibility under seed, per-parameter sample-mean recovery within Monte Carlo error, PSA mean cost / QALY within ~1-2% of the deterministic Table 5 totals (QALYs sit slightly low because `u_H` samples around 0.985 vs the deterministic 1.0), CEAC validity (row sums to 1, SoC dominates at WTP = 0, AB dominates at WTP = 200,000), and Strategy A dominated in expectation.

### Stage 3 — one-way deterministic sensitivity analysis (DSA) and tornado

- `opencea.sensitivity` — One-way DSA runner. For each parameter, the runner sets the value to its low then its high (others held at base), rebuilds the model through the validated `build_darth_sick_sicker` builder, runs the deterministic engine, and records the outcome at each end. Returns a tidy `DSAResult` sorted by swing descending.
- **Outcome of interest**: incremental net monetary benefit of a chosen comparator vs the baseline (SoC by default) at a fixed WTP (default $100k). NMB rather than ICER because ICERs flip sign and go undefined under parameter sweeps, making them unusable in a tornado. The comparator defaults to the strategy with the highest base-case NMB at the WTP (Strategy B at $100k for the DARTH model).
- **Parameter ranges** are derived from the existing PSA distributions in `PSA_PARAM_SPECS` (default 2.5 / 97.5 percentiles, configurable). The range is extended to include the deterministic base value when it falls outside the percentile interval — relevant for `u_H`, whose deterministic value of 1.0 sits above the 97.5th percentile of `Beta(200, 3)`. The tightly specified lognormal HRs (sdlog 0.01-0.02) produce negligible swings, by construction; the tornado is correctly dominated by parameters with real uncertainty (`c_trtB`, `hr_S1S2_trtB`, `r_HD`, disease-progression rates, utilities).
- `opencea.plots.plot_tornado` — Classic CEA tornado. Horizontal bars sorted by swing, centered on the base outcome, two-tone coloring distinguishes outcomes at the low vs high parameter value. Annotated with each parameter's low / high input value.
- `tests/test_sensitivity.py` — Base-case consistency (DSA's `base_outcome` matches a direct engine run on the YAML), per-parameter bracketing of the base, descending-swing ordering, sensitivity ordering sanity (`c_trtB` swing > `hr_S1` swing by >100x), structural checks (Strategy-A-only parameters `c_trtA` and `u_trtA` have exactly zero swing in the B vs SoC comparison), and a tornado plot smoke test.

### Stage 4 — applied case study: empagliflozin vs SoC in T2D + CVD

- `opencea.empagliflozin` — Illustrative 3-state (EF / PE / D), 2-strategy (SoC, Empagliflozin) cohort cost-effectiveness model anchored on EMPA-REG OUTCOME (Zinman 2015). Reuses `rate_to_prob`, the competing-risks construction, the engine, the discount + WCC weights, the PSA sampler, the DSA driver, and the plotting layer — nothing in the simulation logic is reimplemented. The one-time acute-event cost on transitioning EF → PE is implemented as a transition cost discounted at the engine's `dw_c * wcc` weight.
- `examples/empagliflozin_t2d.yaml` — Every parameter cited inline (EMPA-REG OUTCOME rates / HRs, ADA / Red Book / MEPS / Nicholson 2016 costs, UKPDS 62 / Janssen 2022 utilities).
- **Scenario analysis** — `evaluate_scenario`, `scenario_icer`, `breakeven_drug_price`, and `run_empa_psa_scenario` for net / rebated drug pricing and treatment-effect waning. Waning runs via a per-cycle transition-matrix sequence; the additive `simulate_trace_sequence` / `evaluate_sequence` helpers in `opencea.engine` are parallel to the validated time-homogeneous path and ship with a regression test that confirms they reduce to `evaluate_strategy` when fed a stack of identical matrices.
- `examples/empagliflozin_case_study.md` — CHEERS-structured writeup with a scenario grid (WAC sustained / net price sustained / waning WAC / waning net), breakeven price under each effect-duration assumption, and a scenario-driven conclusion: the dominant uncertainty is **the durability of the EMPA-REG all-cause mortality benefit, not drug pricing**. Sustained effect: cost-effective at $100k/QALY at WAC, decisively so at $4,500/yr net. Waning effect: ICER > $200k/QALY at WAC and > $150k/QALY even at $4,500/yr.
- `examples/figures/empa_*.png` — generated tornado, CEAC, CE plane, and frontier figures (base case; scenarios are reproducible from the writeup's snippet).
- `tests/test_empagliflozin_case.py` — 39 tests: structural integrity (states, row-stochastic matrices, PE non-recovery), base-case consistency tying the case evaluator back to a direct engine call, ICER sanity-band check (\$20k - \$120k / QALY), DSA + PSA structural sanity, plot smoke, and scenario tests (engine sequence reduces to time-homogeneous; net-price lowers ICER, waning raises ICER, combined lies between; breakeven recovers $100k target; waning breakeven < sustained breakeven; CEAC ordering net > sustained > waning).

Not yet implemented (later stages): FastAPI backend, Streamlit / Next.js front end, LLM "assumption critic".

## Install

```bash
pip install -e ".[dev]"
```

Requires Python 3.10+.

## Quickstart

Run the test suite:

```bash
pytest
```

Use the engine interactively:

```python
from opencea import build_darth_sick_sicker, run_model
from opencea.cea import cea_table

model = build_darth_sick_sicker("examples/sick_sicker.yaml")
results = run_model(model)
print(cea_table(results, wtp=100_000))
```

Run a PSA and render the CEAC:

```python
from opencea import run_psa, compute_ceac, default_wtp_grid
from opencea.plots import plot_ceac, plot_ce_plane, plot_ce_frontier

psa = run_psa("examples/sick_sicker.yaml", n_sim=10_000, seed=20260625)
grid = default_wtp_grid()                     # 0 to 200,000 in 1,000 steps
ceac = compute_ceac(psa, grid)                # DataFrame of P(cost-effective)
plot_ceac(psa, "ceac.png", wtp_grid=grid)     # hero figure
plot_ce_plane(psa, "ce_plane.png")
plot_ce_frontier(psa, "frontier.png", wtp_grid=grid)
```

Run a one-way DSA and render the tornado:

```python
from opencea import run_dsa
from opencea.plots import plot_tornado

dsa = run_dsa("examples/sick_sicker.yaml", wtp=100_000)
# dsa.comparator -> "Strategy B" (optimal vs SoC at $100k);
# dsa.base_outcome -> ~$39,793 incremental NMB.
print(dsa.to_dataframe().head())
plot_tornado(dsa, "tornado.png")
```

## Reference model

`examples/sick_sicker.yaml` reproduces the DARTH "Sick-Sicker" cohort tutorial:

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H.
> *An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example.*
> Medical Decision Making, 2023; 43(1):3–20.
> Code: <https://github.com/DARTH-git/cohort-modeling-tutorial-intro>.

All four strategies (SoC, A, B, AB) are constructed from a single rate-based parameter file using DARTH's rate-to-probability and competing-risks transition construction. The engine reproduces the manuscript's published total discounted costs and QALYs to the cent, and the published ICERs to the dollar.

## License

MIT. Copyright (c) 2026 Shunji Lewandowski. See [LICENSE](LICENSE).
