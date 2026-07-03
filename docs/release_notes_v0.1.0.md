# OpenCEA 0.1.0 — first public release

## What it is

**OpenCEA** is an open, reproducible, Python-native library for
cohort state-transition cost-effectiveness analysis. It targets the same
work that has historically been done in R (`heemod`, `hesim`, `dampack`)
or proprietary tools (TreeAge), but in the language most data and ML
teams already use.

Version 0.1.0 packages the four internal build phases — validated
deterministic core, PSA/CEAC, DSA/tornado, and an applied case study —
into a single installable library.

## Validation

The deterministic engine reproduces the peer-reviewed **DARTH Sick-Sicker
reference model** (Alarid-Escudero et al., *Medical Decision Making*
2023) to the cent on manuscript Table 5 totals and to the dollar on
manuscript Table 6 ICERs, for all four strategies (SoC, A, B, AB):

| Strategy | Total cost | Total QALYs |
|---|---:|---:|
| SoC | $151,580 | 20.711 |
| A | $284,805 | 21.499 |
| B | $259,100 | 22.184 |
| AB | $378,875 | 23.137 |

| Comparison | ICER ($/QALY) |
|---|---:|
| B vs SoC | 72,988 |
| AB vs B | 125,764 |

These numbers are pinned as golden tests that run on every push and on
every pull request across Python 3.10, 3.11, and 3.12. **117 tests, 95%
line coverage.**

## What's in the box

- **Deterministic engine** — cohort trace, discounting, within-cycle
  correction (Simpson's 1/3, half-cycle, or none); optional time-varying
  transitions for scenarios like treatment-effect waning.
- **Probabilistic sensitivity analysis (PSA)** — vectorized Monte Carlo,
  seeded reproducible sampling, CEAC, cost-effectiveness plane and
  frontier.
- **Deterministic sensitivity analysis (DSA)** — one-way sweeps on
  incremental NMB at a chosen WTP with a clean tornado plot.
- **Applied case study** — empagliflozin vs standard of care for T2D +
  established CVD, anchored on EMPA-REG OUTCOME (Zinman 2015), with
  scenario analysis (net price, treatment-effect waning) and a breakeven
  price solver. Conclusion: the dominant uncertainty is **the durability
  of the mortality benefit, not the drug price**.
- **Strict input validation** via Pydantic (row-stochastic transition
  matrices, probabilities and utilities in [0, 1], non-negative costs,
  dimension consistency).
- **Headless matplotlib plotting** for CE plane, CEAC, frontier, and
  tornado.

## Install

```bash
pip install opencea
```

## Links

- **Repository:** <https://github.com/SHLEW06/opencea>
- **Changelog:** <https://github.com/SHLEW06/opencea/blob/main/CHANGELOG.md>
- **Citation:** [`CITATION.cff`](https://github.com/SHLEW06/opencea/blob/main/CITATION.cff)
- **License:** MIT

## Reference

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
> Pechlivanoglou P, Jalal H. *An Introductory Tutorial on Cohort
> State-Transition Models in R Using a Cost-Effectiveness Analysis
> Example.* Medical Decision Making, 2023; 43(1):3–20.
> DOI: [10.1177/0272989X221103163](https://doi.org/10.1177/0272989X221103163).
