# OpenCEA

> Open, reproducible health-economic decision modeling in Python.

OpenCEA is an open-source Python toolkit for **cohort state-transition cost-effectiveness modeling**. It targets the gap left by R-only stacks (heemod, hesim, dampack) and proprietary tools (TreeAge) by providing a transparent, validated, vectorized engine in Python.

## Status — Stage 1: validated modeling core

This first stage delivers the correctness anchor for the rest of the project. It is intentionally small:

- `opencea.model` — Pydantic specification for a cohort state-transition model with input validation (square transition matrices, rows sum to 1, probabilities/utilities in [0, 1], non-negative costs).
- `opencea.engine` — NumPy cohort trace, discounting, and within-cycle correction (Simpson's 1/3 rule, half-cycle, or none); total discounted costs and QALYs.
- `opencea.cea` — Total cost / QALY, incremental cost / QALY, ICER, net monetary benefit, strict dominance flagging.
- `opencea.builders.build_darth_sick_sicker` — Constructs a `CohortModel` for the DARTH Sick-Sicker tutorial from a rate-based parameter spec (matching the DARTH variable names verbatim) using DARTH's rate-to-probability and competing-risks transition construction.
- `examples/sick_sicker.yaml` — Rate-based parameter file for the DARTH intro-tutorial model (SoC, A, B, AB).
- `tests/test_darth_reference.py` — Golden-reference tests pinned to the published Table 5 / Table 6 totals in the DARTH manuscript (SoC $151,580 / 20.711, A $284,805 / 21.499, B $259,100 / 22.184, AB $378,875 / 23.137; ICER B vs SoC $72,988/QALY, ICER AB vs B $125,764/QALY).

Not yet implemented (later stages): one-way DSA / tornado, PSA, plotting, FastAPI backend, Streamlit / Next.js front end, CHEERS reporting, LLM "assumption critic".

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

## Reference model

`examples/sick_sicker.yaml` reproduces the DARTH "Sick-Sicker" cohort tutorial:

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H.
> *An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example.*
> Medical Decision Making, 2023; 43(1):3–20.
> Code: <https://github.com/DARTH-git/cohort-modeling-tutorial-intro>.

All four strategies (SoC, A, B, AB) are constructed from a single rate-based parameter file using DARTH's rate-to-probability and competing-risks transition construction. The engine reproduces the manuscript's published total discounted costs and QALYs to the cent, and the published ICERs to the dollar.

## License

MIT.
