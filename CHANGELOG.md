# Project history

OpenCEA was built incrementally, with each stage anchored on a validation criterion. The list below is the build sequence; see the [README](README.md) for the current capability summary and [docs/release_notes_v0.1.1.md](docs/release_notes_v0.1.1.md) for the current release notes.

## 0.1.1 - release candidate (unreleased)

This candidate contains the commits made after the local `v0.1.0` tag
and repairs the release path without changing model calculations.

- The public API includes type annotations and the wheel includes the
  `py.typed` marker.
- Ruff, MyPy, pytest, coverage, build, Twine, pre-commit, and pip use
  controlled versions. Ruff uses an explicit rule set in local checks,
  pre-commit, and CI.
- Package metadata now uses an SPDX license expression, declares Python
  3.10 through 3.12, and bounds runtime dependencies by major version.
- The sdist includes the tests, YAML fixtures, examples, release
  documents, and distribution validator needed to test it independently.
- CI and the publishing workflow inspect both archives, test the sdist,
  install the wheel in a clean environment, and run a self-contained
  public API calculation before publication.
- Tornado plots accept optional reader-friendly parameter labels, and the
  public case-study figure uses them.
- The case-study writeup distinguishes cited inputs from modeling assumptions
  and describes uncertainty without overstating the result.

## 0.1.0 — first tagged release (2026-07-03)

First public release. Consolidates the four internal build phases (v0.1–v0.4 below) into a single installable library.

- **Packaging.** `pyproject.toml` with full metadata (Intended Audience :: Science/Research, Medical Science Apps, Python 3.10–3.12 classifiers, keywords, project URLs); dev extra now includes `pytest-cov`, `build`, and `twine`. `src/opencea/__init__.py` `__version__` synced to `0.1.0`. Fresh-venv wheel install + smoke evaluation verified locally.
- **Coverage.** `pytest --cov=opencea` reports **95%** line coverage across 117 tests (uncovered lines are Pydantic validator error branches and a few edge cases in the engine and model).
- **Citation.** `CITATION.cff` (CFF 1.2) with author, version, MIT, repo URL, and a `references` entry for Alarid-Escudero et al. 2023 (the DARTH tutorial).
- **PyPI publishing.** `.github/workflows/publish.yml` uses PyPI Trusted Publishing (OIDC, no stored token); one-time setup documented in `docs/PUBLISHING.md`.
- **Housekeeping.** Repo-root planning doc moved to `docs/planning/`; README picks up a Roadmap section (FastAPI, web UI, LLM assumption critic, microsimulation, planned `opencea-evals` companion) and normalized `SHLEW06/opencea` link casing.
- **No behavior changes.** The DARTH golden tests (manuscript Tables 5 and 6, to the cent and to the dollar) remain the release gate and are untouched.

## v0.4 — applied case study: empagliflozin vs SoC in T2D + CVD

- `opencea.empagliflozin` — illustrative 3-state (EF / PE / D), 2-strategy cohort cost-effectiveness model anchored on EMPA-REG OUTCOME (Zinman 2015). Reuses `rate_to_prob`, the competing-risks construction, the engine, the discount + WCC weights, the PSA sampler, the DSA driver, and the plotting layer — nothing in the simulation logic is reimplemented. The one-time acute-event cost on transitioning EF → PE is implemented as a transition cost discounted at the engine's `dw_c * wcc` weight.
- `examples/empagliflozin_t2d.yaml` — every parameter cited inline (EMPA-REG OUTCOME rates / HRs, ADA / Red Book / MEPS / Nicholson 2016 costs, UKPDS 62 / Janssen 2022 utilities).
- **Scenario analysis** — `evaluate_scenario`, `scenario_icer`, `breakeven_drug_price`, and `run_empa_psa_scenario` for net / rebated drug pricing and treatment-effect waning. Waning runs via a per-cycle transition-matrix sequence; the additive `simulate_trace_sequence` / `evaluate_sequence` helpers in `opencea.engine` are parallel to the validated time-homogeneous path and ship with a regression test that confirms they reduce to `evaluate_strategy` when fed a stack of identical matrices.
- `examples/empagliflozin_case_study.md` — CHEERS-structured writeup with a scenario grid (WAC sustained / net price sustained / waning WAC / waning net), breakeven price under each effect-duration assumption, and a scenario-driven conclusion: the dominant uncertainty is **the durability of the EMPA-REG all-cause mortality benefit, not drug pricing**.
- `examples/figures/empa_*.png` — generated tornado, CEAC, CE plane, and frontier figures (base case; scenarios are reproducible from the writeup's snippet).
- `tests/test_empagliflozin_case.py` — 39 tests: structural integrity, base-case consistency, ICER sanity-band check, DSA + PSA structural sanity, plot smoke, and scenario tests (engine sequence reduces to time-homogeneous; net-price lowers ICER, waning raises ICER, combined lies between; breakeven recovers $100k target; waning breakeven < sustained breakeven; CEAC ordering net > sustained > waning).

## v0.3 — one-way deterministic sensitivity analysis (DSA) and tornado

- `opencea.sensitivity` — DSA runner. For each parameter, the runner sets the value to its low then its high (others held at base), rebuilds the model through the validated builder, runs the deterministic engine, and records the outcome at each end. Returns a tidy `DSAResult` sorted by swing descending.
- **Outcome of interest**: incremental net monetary benefit of a chosen comparator vs the baseline (SoC by default) at a fixed WTP (default $100k). NMB rather than ICER because ICERs flip sign and go undefined under parameter sweeps, making them unusable in a tornado. The comparator defaults to the strategy with the highest base-case NMB at the WTP (Strategy B at $100k for the DARTH model).
- **Parameter ranges** are derived from the existing PSA distributions in `PSA_PARAM_SPECS` (default 2.5 / 97.5 percentiles, configurable). The range is extended to include the deterministic base value when it falls outside the percentile interval — relevant for `u_H`, whose deterministic value of 1.0 sits above the 97.5th percentile of `Beta(200, 3)`.
- `opencea.plots.plot_tornado` — classic CEA tornado. Horizontal bars sorted by swing, centered on the base outcome, two-tone coloring distinguishes outcomes at the low vs high parameter value. Annotated with each parameter's low / high input value.
- `tests/test_sensitivity.py` — base-case consistency (DSA's `base_outcome` matches a direct engine run on the YAML), per-parameter bracketing of the base, descending-swing ordering, sensitivity ordering sanity (`c_trtB` swing > `hr_S1` swing by >100x), structural checks (Strategy-A-only parameters `c_trtA` and `u_trtA` have exactly zero swing in the B vs SoC comparison), and a tornado plot smoke test.

## v0.2 — probabilistic sensitivity analysis (PSA) and CEAC

- `opencea.psa` — PSA distribution spec, seeded sampler, and a vectorized Monte Carlo runner. Distributions are transcribed verbatim from the DARTH `generate_psa_params()` function (`R/Functions_cSTM_time_indep.R`, lines 278-309), with R's rate-parameterized `rgamma` converted to numpy's scale parameterization. The per-cycle trace step `trace[t+1] = trace[t] @ P` is batched via `einsum` over an `(n_sim, n_states, n_states)` transition tensor so the inner loop runs `T` matmuls of width `n_sim` instead of `n_sim * T` scalar matmuls.
- `compute_nmb`, `compute_ceac`, `expected_nmb_frontier`, `incremental_vs_baseline` — net monetary benefit per draw per strategy across a WTP grid (default 0 to 200,000), the cost-effectiveness acceptability curve, the expected-NMB-optimal strategy at each WTP, and the per-draw incremental cost / QALY vs SoC for the CE plane.
- `opencea.plots` — matplotlib (headless `Agg` backend) figures: cost-effectiveness plane scatter, CEAC, and cost-effectiveness frontier. Output is saved to file paths chosen by the caller.
- `tests/test_psa.py` — reproducibility under seed, per-parameter sample-mean recovery within Monte Carlo error, PSA mean cost / QALY within ~1-2% of the deterministic Table 5 totals (QALYs sit slightly low because `u_H` samples around 0.985 vs the deterministic 1.0), CEAC validity (row sums to 1, SoC dominates at WTP = 0, AB dominates at WTP = 200,000), and Strategy A dominated in expectation.

## v0.1 — validated deterministic core

All four DARTH Sick-Sicker strategies (SoC, A, B, AB) reproduce the published manuscript totals to the cent:

- `opencea.model` — Pydantic specification for a cohort state-transition model with input validation (square transition matrices, rows sum to 1, probabilities/utilities in [0, 1], non-negative costs).
- `opencea.engine` — NumPy cohort trace, discounting, and within-cycle correction (Simpson's 1/3 rule, half-cycle, or none); total discounted costs and QALYs.
- `opencea.cea` — total cost / QALY, incremental cost / QALY, ICER, net monetary benefit, strict dominance flagging.
- `opencea.builders.build_darth_sick_sicker` — constructs a `CohortModel` for the DARTH Sick-Sicker tutorial from a rate-based parameter spec using DARTH's rate-to-probability and competing-risks transition construction.
- `examples/sick_sicker.yaml` — rate-based parameter file for the DARTH intro-tutorial model.
- `tests/test_darth_reference.py` — golden-reference tests pinned to manuscript Table 5 / Table 6 (SoC $151,580 / 20.711, A $284,805 / 21.499, B $259,100 / 22.184, AB $378,875 / 23.137; ICER B vs SoC $72,988/QALY, ICER AB vs B $125,764/QALY).
