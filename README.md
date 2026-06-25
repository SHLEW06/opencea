# OpenCEA

> Open, reproducible health-economic decision modeling in Python.

OpenCEA is an open-source Python toolkit for **cohort state-transition cost-effectiveness modeling**. It targets the gap left by R-only stacks (heemod, hesim, dampack) and proprietary tools (TreeAge) by providing a transparent, validated, vectorized engine in Python.

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

Not yet implemented (later stages): one-way DSA / tornado, FastAPI backend, Streamlit / Next.js front end, CHEERS reporting, LLM "assumption critic".

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

## Reference model

`examples/sick_sicker.yaml` reproduces the DARTH "Sick-Sicker" cohort tutorial:

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H.
> *An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example.*
> Medical Decision Making, 2023; 43(1):3–20.
> Code: <https://github.com/DARTH-git/cohort-modeling-tutorial-intro>.

All four strategies (SoC, A, B, AB) are constructed from a single rate-based parameter file using DARTH's rate-to-probability and competing-risks transition construction. The engine reproduces the manuscript's published total discounted costs and QALYs to the cent, and the published ICERs to the dollar.

## License

MIT.
