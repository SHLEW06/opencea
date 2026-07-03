# Flagship Build Plan: OpenCEA — Open, Reproducible Health-Economic Decision Modeling

> **Use this file to seed a new chat or Claude Code session.** It supersedes the prior "HealthIntel + Second Flagship" handoff. It records the decision, the rationale, and a staged, feasible build plan.

---

## 0. Decision summary

After evaluating the candidate list, the chosen flagship is a **Python-native, open-source health-economic decision-modeling toolkit**: Markov / cohort state-transition cost-effectiveness models with probabilistic sensitivity analysis (PSA), cost-effectiveness acceptability curves (CEACs), and CHEERS-compliant reporting.

**Working name:** OpenCEA (alternatives: CEAflow, DecisionRx, HEORkit).

**Why this one (and why the others were cut):**
- The LLM-evaluation and RAG space is fully commoditized (DeepEval, Langfuse, Phoenix, Promptfoo, RAGAS, Braintrust). A generic eval platform or document-RAG assistant reads as "didn't know what already exists."
- Recommenders, anomaly detectors, search engines, and forecasters are the standard bootcamp portfolio — they signal conformity, not originality.
- Serious cost-effectiveness modeling lives almost entirely in **R** (heemod, hesim, BCEA, dampack) or expensive commercial software (**TreeAge**). The field's documented pain points are reproducibility and transparency. **Python is conspicuously under-served.** This is a real, current gap that sits exactly on the Economics + Human Health background.
- It is the single candidate that serves **consulting/HEOR and tech/data/AI-product equally**: to a HEOR reader it *is* the discipline; to an engineer it's a vectorized simulation engine + Monte Carlo layer + API + tests + CI, validated against a peer-reviewed model.
- **Moat:** a strong coder cannot easily build this without knowing what a QALY, ICER, discount rate, or value-of-information analysis is.

**Time budget:** ~12–18 focused hrs/week for ~12 weeks (some heavier weeks). Staged so a strong MVP ships before expansion.

---

## 1. Positioning (dual-audience)

**HEOR / consulting framing:**
> An open, transparent, reproducible cost-effectiveness modeling toolkit — the Python answer to TreeAge — with full probabilistic sensitivity analysis and CHEERS-compliant reporting.

**Tech / data / AI-product framing:**
> A vectorized numerical simulation engine with a Monte Carlo uncertainty layer, exposed via a FastAPI backend and web app, validated against a peer-reviewed reference model, with a full pytest suite and CI.

Same repo. Lead with whichever framing matches the reader.

---

## 2. The credibility anchor: validate against a published reference model

The most important strategic move in the whole project.

Reproduce the **DARTH "Sick-Sicker" cohort state-transition model** in Python and verify the outputs (cohort trace, QALYs, costs, ICERs, PSA results) **numerically match the published reference**.

- Reference: Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. *An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example.* Medical Decision Making, 2023;43(1):3–20.
- Code (R) + parameters: `github.com/DARTH-git/cohort-modeling-tutorial-intro` (archived on Zenodo). A time-dependent companion exists at `.../cohort-modeling-tutorial-timedep`.
- Model: 4 states — Healthy (H), Sick (S1), Sicker (S2), Dead (D); strategies SoC, A, B, AB; 4×4 transition matrix; cohort trace × state-reward vectors → costs and QALYs.

**Why this matters:** it proves the engine is correct (guarding against the off-by-one cycle, half-cycle-correction, and discount-timing bugs that quietly break these models), and it shows the actual academic literature was read. Build this validation in the MVP, not at the end.

---

## 3. Architecture & stack (with rationale)

- **Core library — Python + NumPy.** Vectorized matrix algebra for the cohort trace (`trace[t+1] = trace[t] @ P`) and the PSA loop. This is where the real engineering lives. SciPy for distributions. pandas for I/O only.
- **Model spec — Pydantic dataclasses (or YAML/JSON validated by Pydantic).** Define states, transition matrix (constant or parameterized), state rewards (cost + utility), strategies, discount rate, cycle length, horizon. Pydantic validation enforces integrity (transition rows sum to 1, probabilities and utilities in valid ranges) — a clean "data integrity" signal that also connects to prior reconciliation/QA experience.
- **Plotting — matplotlib.** Cohort trace, tornado diagram, cost-effectiveness plane, CEAC, frontier.
- **API — FastAPI.** Submit a model spec, get back CEA results and PSA outputs. The right, modern choice.
- **Frontend — Streamlit for the MVP demo; optional Next.js/React for the advanced version.** Streamlit ships fast. The Next.js upgrade leverages existing front-end skill and raises the tech-lane signal — but only after the core is solid. Do not start here.
- **Reproducibility — fixed seeds, config-driven runs, `pyproject.toml`, a CLI (`python -m opencea run config.yaml`), optional Docker.**
- **Testing — pytest.** Golden-output tests vs. DARTH numbers; property tests (rows sum to 1; total QALYs fall as discount rate rises; dominated strategies flagged correctly).
- **CI — GitHub Actions** running tests on every push. Cheap, high-signal "real engineer" move.

---

## 4. Repo structure

```
opencea/
  pyproject.toml
  README.md
  LICENSE                 # MIT
  .github/workflows/ci.yml
  src/opencea/
    __init__.py
    model.py              # Pydantic model spec + validation
    engine.py             # cohort trace, half-cycle correction, discounting
    cea.py                # ICER, NMB, dominance / extended dominance
    sensitivity.py        # one-way DSA, tornado, two-way
    psa.py                # distributions, Monte Carlo, vectorized PSA
    plots.py              # trace, tornado, CE plane, CEAC, frontier
    report.py             # CHEERS-structured markdown/PDF report
    cli.py
  tests/
    test_engine.py
    test_cea.py
    test_psa.py
    test_darth_reference.py   # golden tests vs. published Sick-Sicker outputs
  examples/
    sick_sicker.yaml          # the validation model
    <applied_case>.yaml       # the real case study
    walkthrough.ipynb         # end-to-end narrated notebook
  app/                        # Streamlit (MVP) / Next.js (advanced)
  docs/
    methodology.md            # assumptions, equations, limitations
    architecture.md           # diagram + design decisions
```

---

## 5. Feature stages

### MVP — weeks 1–5 (the strong, shippable core)
1. **Cohort engine** (time-homogeneous): cohort trace, half-cycle correction, discounting of costs and QALYs.
2. **CEA outputs:** total discounted cost & QALYs per strategy; incremental analysis; **ICER = ΔCost / ΔQALY**; **Net Monetary Benefit = WTP × QALY − Cost**; strict and extended dominance detection.
3. **DARTH validation:** reproduce Sick-Sicker (SoC, A, B, AB); outputs match the published reference. *(Correctness story.)*
4. **One-way deterministic sensitivity analysis** + tornado diagram.
5. **Narrated notebook** walking one model end to end.
6. **Basic Streamlit app:** load a config, view trace + CEA table + tornado.
7. **Tests + CI + a real README.**

### Advanced — weeks 6–9
8. **Probabilistic sensitivity analysis:** parameter distributions — Beta for probabilities/utilities (bounded 0–1), Gamma or LogNormal for costs (positive, right-skewed), Dirichlet for transition probabilities out of a state, LogNormal for hazard/relative risks. Vectorized Monte Carlo.
9. **Cost-effectiveness plane** (scatter of PSA draws) + **CEAC** (probability each strategy is cost-effective across WTP thresholds) + cost-effectiveness frontier.
10. **Two-way sensitivity analysis / scenario comparison.**
11. **CHEERS-structured auto-report** (markdown → PDF): assumptions, parameters, results, the iconic figures. *(Domain-credibility signal.)*
12. **FastAPI backend**; optionally a **Next.js** front end for a polished demo.
13. **One real applied case study** using public, cited parameters (see §6). *This is what makes it memorable.*

### Stretch — weeks 10–12 (pick 1–2)
14. **Expected Value of Perfect Information (EVPI)** — value-of-information analysis. Graduate-level signal.
15. **Time-dependent transitions** (state-residence / tunnel states), matching the DARTH time-dependent tutorial.
16. **LLM "assumption critic"** — given a model spec, generate a plain-language summary of assumptions and flag the parameters driving the most uncertainty (tie to DSA/PSA output). A *garnish*, not a chatbot. Cap at ~2 days.

---

## 6. The applied case study (do not skip)

A sterile library is forgettable; a tool used to answer a real question is not. In the advanced tier, parameterize one real intervention from published literature, e.g.:
- a vaccination or screening program (well-documented parameters),
- a chronic-disease treatment vs. standard of care.

**Rules:** use published, cited values for transition probabilities, costs, and utilities. State explicitly that values are illustrative and drawn from the literature. **Never fabricate parameters** — fabrication is the single most common way these projects lose credibility. Frame the writeup around a decision question ("Is intervention X cost-effective at a $100k/QALY threshold, and how robust is that conclusion to uncertainty?").

---

## 7. Honest risks (be critical)

- **Parameter sourcing** is the real time sink. Mitigate by using the given Sick-Sicker parameters for the core and a single well-documented intervention for the case study.
- **Scope creep.** This domain is bottomless (microsimulation, individual-level models, survival modeling, network meta-analysis). The cohort model + PSA + CEAC + one real case is a *complete* flagship. Everything else is stretch.
- **Correctness bugs** (cycle indexing, half-cycle correction, discount timing) are easy and embarrassing. The DARTH validation is the guardrail — build it early.
- **The AI feature** will tempt over-investment. It is optional and small.
- **Reality check on SWE:** this does not substitute for data-structures-and-algorithms prep if chasing pure software-engineering interviews. It makes for a standout applied/quant/HEOR/AI-product candidate with one genuinely deep technical artifact, and it makes the profile *interesting*.

---

## 8. README outline (recruiters scan in ~30 seconds — front-load impact)

1. **One-line pitch** + a single hero figure (the CEAC).
2. **What it does**, in three sentences, with the validation claim up front ("validated against the peer-reviewed DARTH reference model — outputs match").
3. **Why it exists** (the R/TreeAge gap; reproducibility and transparency).
4. **Quickstart** (install, run the example, see results).
5. **Features** (engine, PSA, CEAC, CHEERS report, API).
6. **The applied case study** (the decision question + headline result).
7. **Validation & tests** (how correctness is established; CI badge).
8. **Methodology & limitations** (link to `docs/methodology.md`).
9. **Architecture** (diagram + key design decisions).

---

## 9. Resume bullets (non-inflated, dual-framed)

- Built **OpenCEA**, an open-source Python toolkit for health-economic cost-effectiveness modeling (Markov cohort state-transition models, probabilistic sensitivity analysis, cost-effectiveness acceptability curves), **validated against the peer-reviewed DARTH reference model with numerically matching outputs**.
- Implemented a **vectorized Monte Carlo PSA engine** (Beta/Gamma/Dirichlet parameter distributions) and automated **CHEERS-compliant reporting**, served via a **FastAPI** backend and interactive web app.
- Engineered a **reproducible, test-covered pipeline** (pytest + GitHub Actions CI) and applied it to a **[intervention] cost-effectiveness case study** using published clinical and cost parameters, quantifying decision uncertainty across willingness-to-pay thresholds.

---

## 10. How to pitch it in interviews

- **Consulting / HEOR:** lead with the decision question and the uncertainty story ("the ICER was $X/QALY, but the CEAC showed only a 60% probability of cost-effectiveness at $100k — so the recommendation depends on the threshold and on parameter Y"). Emphasize transparency, reproducibility, and that it matches a published model.
- **Tech / data / AI-product:** lead with the engineering ("vectorized simulation, Monte Carlo over parameter distributions, validated against a reference, tested, CI'd, exposed via an API"). Mention the design decision you're proudest of and one hard bug you caught via the validation.

---

## 11. First task to hand Claude Code

> "Scaffold the `opencea` repo per the structure in this plan. Implement `model.py` (Pydantic spec with validation: transition rows sum to 1, utilities in [0,1]) and `engine.py` (cohort trace with half-cycle correction and discounting). Then implement the DARTH Sick-Sicker model as `examples/sick_sicker.yaml` and write `tests/test_darth_reference.py` asserting the cohort trace, total discounted costs, and total discounted QALYs for strategies SoC and A. Use the published parameters from the DARTH cohort-modeling-tutorial-intro repo."

---

## 12. The rest of the portfolio (context)

This flagship is one of a small, deliberately uncluttered set. The supporting pieces are **finish-and-frame jobs, not new builds**:
- **Chinese Adaptive Reader** (Next.js / Firebase / Gemini): add an evaluation layer (does difficulty adaptation measurably work?) + a strong README → a credible "AI product beyond a demo chatbot."
- **YOLO defect detection** (manufacturing CV): document the image-processing pipeline and report precision/recall → a real applied-ML artifact.
- **HealthIntel:** if kept at all, amputate to a single analytical spine (e.g., a county preventive-care priority score with honest sensitivity analysis and one good map). Do not build seven modules.

Quality over quantity. 2–4 sharp, well-presented projects beat a cluttered profile every time.
