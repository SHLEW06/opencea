# Empagliflozin vs Standard of Care in T2D with established CVD

**Illustrative cost-effectiveness analysis** built on the validated OpenCEA engine. Every parameter is drawn from cited literature (see [`examples/empagliflozin_t2d.yaml`](empagliflozin_t2d.yaml)); no values are invented.

This writeup follows the CHEERS 2022 reporting structure adapted for a brief working document.

---

## Decision problem

Is empagliflozin added to standard of care cost-effective versus standard of care alone for adults with type 2 diabetes (T2D) and established cardiovascular disease (CVD) at a willingness-to-pay (WTP) of **$100,000 / QALY**, and how robust is that conclusion to plausible parameter uncertainty?

| Element | Choice |
|---|---|
| Comparators | Standard of care (SoC); Empagliflozin + SoC |
| Population | Adults with T2D + established CVD, mean age 63 (EMPA-REG OUTCOME placebo arm) |
| Perspective | US healthcare-sector (drug + medical care; non-medical costs excluded) |
| Time horizon | Lifetime — start age 63, terminate at age 100 (37 annual cycles) |
| Discount rate | 3% per year for costs and QALYs |
| Cycle length | 1 year, Simpson's 1/3 within-cycle correction |
| Outcome | Incremental cost / QALY (ICER); net monetary benefit at $100k / QALY; CEAC over 0-$200k / QALY |

---

## Model structure

Three health states:

```
        +-- HR_event ---+
        v               |
   [ EF ]----acute-->[ PE ]
     |                  |
     +--HR_death---+    +--HR_death---+
                   v                  v
                 [ D ] <----------- [ D ]
```

- **EF** — event-free: T2D with CVD but no major CV event yet.
- **PE** — post major CV event (composite of MI / stroke / HF hospitalization).
- **D** — dead, absorbing.

PE cannot recover to EF (mirroring the "Sicker" structure in the DARTH tutorial). Transition probabilities are built from annual rates via `p = 1 - exp(-r * t)` with the same competing-risks construction the engine already validates — non-death transitions are scaled by `(1 - p_death)` so each row sums to 1.

The acute event cost on entering PE is implemented as a **transition cost**: discounted at the engine's `dw_c[t+1] * wcc[t+1]` weight, matching how state costs are discounted at the cycle when occupants are present.

---

## Two-channel treatment effect

Empagliflozin acts on the model through two independent channels, each anchored on a different EMPA-REG OUTCOME endpoint:

1. **Progression channel (`HR_event = 0.86`)** — composite hazard ratio for first major CV event (`EF -> PE`). Constructed as the EMPA-REG composition-weighted blend `0.60 * 1.00 (MI/stroke, no significant trial effect) + 0.40 * 0.65 (HF hospitalization HR 0.65 [0.50-0.85])`. 95% CI for the composite ~ `[0.74, 0.99]`.
2. **Mortality channel (`HR_death = 0.68`)** — all-cause mortality HR (95% CI `[0.57, 0.82]`, Zinman 2015). Applied to **both** `EF -> D` and `PE -> D` rates, on the assumption (per the case-study spec) that the trial-reported all-cause mortality benefit persists across both health states. *This is the assumption with the largest economic leverage* (see tornado below); a more conservative assumption that applies the HR only in EF would shift the ICER downward.

No utility effect — empagliflozin's QALY benefit comes through event avoidance and survival, not through state-utility change.

---

## Parameters

All values, with citations, are inline in [`empagliflozin_t2d.yaml`](empagliflozin_t2d.yaml). Summary:

| Parameter | Base | Sensitivity range | Source |
|---|---:|---:|---|
| `r_EF_PE` (first event / yr) | 0.0351 | trial-derived | EMPA-REG OUTCOME placebo (Zinman 2015, Tables 2, S5) |
| `r_EF_D` (event-free CV death / yr) | 0.0196 | trial-derived | EMPA-REG OUTCOME placebo CV death |
| `r_PE_D` (post-event death / yr) | 0.060 | 0.04 - 0.10 | Brinkert 2017 (post-MI/stroke) + US HF registries; soft input |
| `HR_event` | 0.86 | 0.74 - 0.99 | EMPA-REG composite + HF HR (Zinman 2015) |
| `HR_death` | 0.68 | 0.57 - 0.82 | EMPA-REG all-cause mortality (Zinman 2015) |
| `c_drug` (USD / yr) | 6,264 | 4,872 - 7,596 | $522/mo WAC, Red Book / JAHA 2024; FSS $406 - AWP $633 |
| `c_EF` (USD / yr) | 16,000 | 13,000 - 19,700 | ADA Economic Costs of Diabetes 2022 |
| `c_PE` (USD / yr ongoing) | 20,000 | gamma 25% CV | 60/40 blend post-MI/stroke + post-HF (MEPS / medRxiv 2025) |
| `c_acute_PE` (USD one-time) | 11,650 | gamma 25% CV | Nicholson 2016 acute-care costs (MI / stroke / HF avg) |
| `u_EF` | 0.75 | beta(75, 25) | Janssen 2022 pooled T2D EQ-5D; UKPDS 62 baseline |
| `u_PE` | 0.68 | beta(68, 32) | UKPDS 62 event decrements weighted 60/40 |

PSA distributions: HRs lognormal parameterized so that `median = point estimate` and `95% interval = trial CI`; costs gamma with mean = base value and 25% CV; utilities beta with `a + b = 100` precision; `r_PE_D` wide lognormal spanning the 0.04 - 0.10 / yr band. PSA holds the baseline event and EF-mortality rates fixed at their trial-derived point estimates.

---

## Deterministic results (base case)

| Strategy | Total cost (USD) | Total QALY | Inc cost | Inc QALY | ICER ($/QALY) | NMB @ $100k |
|---|---:|---:|---:|---:|---:|---:|
| Standard of care | 259,540 | 10.944 | — | — | — | 834,855 |
| Empagliflozin | 397,423 | 12.338 | 137,883 | 1.394 | **98,900** | 836,388 |

Incremental NMB of Empa vs SoC at $100k / QALY: **$1,534** (very thin margin).

The acute event cost contributes $4,449 to SoC and $4,325 to Empa (Empa's lower EF→PE transition rate slightly reduces acute event spend, but the difference is tiny relative to drug cost).

---

## One-way DSA — tornado on incremental NMB at $100k / QALY

Top 5 drivers, sorted by swing (full table in `tornado.png`):

| Rank | Parameter | Low → High value | Inc NMB at low | Inc NMB at high | Swing (USD) |
|---:|---|---|---:|---:|---:|
| 1 | `hr_death` | 0.57 → 0.82 | +35,815 | −33,080 | **68,895** |
| 2 | `c_drug` | 4,872 → 7,596 | +24,501 | −21,434 | 45,935 |
| 3 | `hr_event` | 0.74 → 0.99 | +14,891 | −12,301 | 27,192 |
| 4 | `u_EF` | 0.66 → 0.83 | −10,632 | +13,329 | 23,961 |
| 5 | `r_PE_D` | 0.04 → 0.10 | +12,041 | −7,943 | 19,984 |

Reading: `hr_death` flips the decision (the bar straddles zero). Empa goes from clearly cost-effective at the optimistic CI bound to clearly not cost-effective at the pessimistic CI bound. The drug cost is the second-largest single lever — bringing the price down to $4,872 (FSS) makes Empa cost-effective at $100k WTP without any other change. `c_acute_PE` has near-zero swing because the difference in EF→PE flow between strategies is small.

`hr_event` ranks third even though the composite HR is closer to 1 than `hr_death`. That is consistent with empagliflozin's CV benefit in EMPA-REG being driven more by mortality than by non-fatal events.

---

## PSA — CEAC and CE plane (n_sim = 10,000)

| WTP ($/QALY) | P(Empa cost-effective) |
|---:|---:|
| 50,000 | 0.002 |
| **100,000** | **0.517** |
| 150,000 | 0.945 |

The deterministic ICER ($98,900) sits a hair below $100k WTP, so the CEAC crosses 50% almost exactly at the standard US threshold — confirming the result is genuinely on the cost-effectiveness boundary rather than safely on one side. PSA mean ICER (cost diff / qaly diff in expectation) tracks the deterministic value within ~3%.

The CE plane scatter (`empa_ce_plane.png`) lies almost entirely in the upper-right quadrant — Empa is dominantly more costly and more effective than SoC; the question is only whether the QALY gain is worth the spend.

---

## Scenario analysis

The base case sits *right at* the $100k/QALY threshold, so the conclusion turns on two assumptions flagged in the limitations: WAC drug pricing and lifetime extrapolation of the EMPA-REG mortality benefit. Both deserve explicit scenarios.

### Scenarios

- **Net price** — annual drug cost lowered to **$4,500/yr**, illustrative of a typical ~28% rebate off the $6,264 WAC. *Real net prices are confidential and vary by payer*; this is a single round number, not a quotation.
- **Treatment-effect waning** — both `hr_event` and `hr_death` are held at their trial point estimates for the first **3 years** (EMPA-REG follow-up), then linearly interpolated to **HR = 1.0 by year 10**, then no effect for the remaining 27 cycles. Implemented via a per-cycle transition-matrix sequence; the time-homogeneous engine path is untouched (a regression test confirms the sequence helper reduces to the validated `evaluate_strategy` when fed a stack of identical matrices).

### Scenario grid (n_sim = 10,000)

| Scenario | ICER ($/QALY) | P(CE at $100k) |
|---|---:|---:|
| Base case (WAC, sustained effect) | **98,900** | **0.517** |
| Net price ($4,500/yr, sustained) | 77,564 | 0.920 |
| Waning effect (WAC) | 206,778 | 0.005 |
| Waning + net price | 154,665 | 0.001 |

### Breakeven price

The annual empagliflozin price at which the deterministic ICER equals exactly $100k/QALY:

| Effect assumption | Breakeven price (USD/yr) |
|---|---:|
| Sustained effect | **6,355** |
| Waning effect | **2,650** |

Read: under the sustained-effect base case, empagliflozin would need to come in below $6,355/yr (a ~1% rebate off WAC) to clear $100k/QALY. Under waning, the price would need to fall by **~58%** off WAC to clear the same threshold.

### Interpretation

The two assumptions move in opposite directions and the size of the move is asymmetric. The net-price scenario reduces incremental cost by roughly the discounted lifetime drug-cost differential (~$37k) and pulls the ICER 21% below the base case. Waning, in contrast, eliminates most of the long-run mortality benefit — the incremental QALY gain drops from 1.39 to about 0.49 (a ~65% reduction). Cost savings from fewer high-cost PE years partially offset this, but not enough; the waning ICER more than doubles.

The combined scenario shows that **rebates alone do not rescue a waning-effect interpretation** of the trial. If the EMPA-REG mortality benefit truly persists for life, empagliflozin is cost-effective at $100k/QALY at WAC (with healthy PSA support — 52% at WAC, 92% at $4,500 net). If the mortality benefit fades within ~10 years, the drug is not cost-effective at $100k/QALY even after a ~28% rebate — and would need to drop below ~$2,650/yr to clear the threshold.

The cost-effectiveness conclusion is therefore **conditional on the durability of the EMPA-REG all-cause mortality benefit, not on drug pricing**. Trial-effect durability — not the price — is the dominant uncertainty.

---

## Limitations

The user spec explicitly framed this as an **illustrative** model. Limitations material to the result:

1. **Aggregated event state.** MI, stroke, and HF are collapsed into a single "PE" state with a blended ongoing cost, acute cost, and utility. Published US empagliflozin CEAs in the **$26k - $88k / QALY** range typically separate these (often into 4-6 states), which yields a tighter cost ledger and different mortality trajectories for each event type.
2. **Cohort, not microsimulation.** No individual heterogeneity, no time since first event, no second-event modeling. The 60/40 MI-stroke vs HF blend is fixed at the trial composition forever — in reality the mix would evolve with age and prior events.
3. **Lifetime extrapolation of trial-derived HRs.** EMPA-REG OUTCOME measured outcomes over a median 3.1 years. The base case applies both HRs (`hr_event` = 0.86, `hr_death` = 0.68) for all 37 cycles. The **waning scenario** above quantifies the impact of that assumption — the deterministic ICER more than doubles, from $98,900 to $206,778/QALY, when the effect is allowed to taper to zero by year 10.
4. **All-cause mortality HR applied to PE→D.** The trial measured HR 0.68 in the trial population (which started event-free). Applying that HR to post-event mortality is a strong assumption — the DSA shows it's the single largest driver of incremental NMB, with a tornado bar that straddles zero across the trial's 95% CI.
5. **WAC drug price.** $6,264/yr is the listed wholesale price. Most US payers see meaningfully lower **net** prices after rebates; the net-price scenario above shows that a $4,500/yr price (~28% rebate) lowers the ICER to $77,564/QALY under sustained effect. Real net prices are confidential.
6. **No adverse-event modeling.** Genital infections (the most notable EMPA-REG side effect) and DKA risk are not separately costed.
7. **Linear time-homogeneity in baseline rates.** SoC rates do not change with age, time since diagnosis, or cumulative drug exposure. (The Empa transition matrix can vary per cycle via the additive sequence engine — used here for waning.)

## Bottom line

The base-case ICER ($98,900/QALY) lands at the upper end of plausibility, just below a $100k WTP. The scenario grid makes the driver of that conclusion explicit:

- **If the EMPA-REG mortality benefit is sustained**, empagliflozin is cost-effective at $100k/QALY at WAC pricing (with PSA support of ~52% rising to ~92% at a $4,500 net price).
- **If the effect wanes** (linearly to zero between years 3 and 10), the ICER rises above $200k/QALY at WAC and stays above $150k/QALY even at a $4,500 net price.

The dominant uncertainty is therefore the **durability of the all-cause mortality benefit**, not drug pricing. Two structural choices distinguish this analysis from the central published $26-88k/QALY range: WAC drug pricing (relaxed in the net-price scenario) and the 3-state aggregation (left in place as an illustrative simplification).

---

## Figures

Generated by the OpenCEA plotting module from the same `PSAResult` / `DSAResult` objects the tests verify. Saved to `examples/figures/`:

- `empa_tornado.png` — one-way DSA tornado on incremental NMB at $100k.
- `empa_ceac.png` — CEAC for SoC vs Empa across $0 - $200k WTP.
- `empa_ce_plane.png` — incremental cost / QALY scatter, 10,000 draws.
- `empa_ce_frontier.png` — expected NMB vs WTP, with optimal-strategy switch marked.

To regenerate:

```python
from opencea.empagliflozin import (
    evaluate_empagliflozin_case, run_empa_psa, run_empa_psa_scenario,
    scenario_icer, breakeven_drug_price, WaningSpec,
    dsa_evaluator, EMPA_PSA_SPECS, EMPA_DSA_RANGE_OVERRIDES, SOC,
)
from opencea.sensitivity import run_dsa
from opencea.plots import plot_tornado, plot_ceac, plot_ce_plane, plot_ce_frontier
from opencea.psa import default_wtp_grid

YAML = "examples/empagliflozin_t2d.yaml"

# Base-case figures
psa = run_empa_psa(YAML, n_sim=10_000, seed=20260626)
dsa = run_dsa(
    base_params=YAML, wtp=100_000, baseline=SOC,
    sweep_params=list(EMPA_PSA_SPECS.keys()),
    evaluator=dsa_evaluator, param_specs=EMPA_PSA_SPECS,
    param_ranges=EMPA_DSA_RANGE_OVERRIDES,
)
plot_tornado(dsa, "examples/figures/empa_tornado.png")
plot_ceac(psa, "examples/figures/empa_ceac.png", wtp_grid=default_wtp_grid())
plot_ce_plane(psa, "examples/figures/empa_ce_plane.png")
plot_ce_frontier(psa, "examples/figures/empa_ce_frontier.png")

# Scenarios
waning = WaningSpec(start_year=3.0, end_year=10.0)
for label, dp, wan in [
    ("net_price",    4500.0, None),
    ("waning",        None,  waning),
    ("waning_net",  4500.0, waning),
]:
    icer = scenario_icer(YAML, drug_price=dp, waning=wan)
    print(f"{label}: ICER ${icer:,.0f}/QALY")
print(f"Breakeven (sustained): ${breakeven_drug_price(YAML, 100_000):,.2f}/yr")
print(f"Breakeven (waning):    ${breakeven_drug_price(YAML, 100_000, waning=waning):,.2f}/yr")
```

---

## References

- Zinman B et al. *Empagliflozin, Cardiovascular Outcomes, and Mortality in Type 2 Diabetes.* NEJM 2015; 373:2117-2128.
- American Diabetes Association. *Economic Costs of Diabetes in the U.S. in 2022.* Diabetes Care 2024.
- Clarke P et al. *Estimating the cost-effectiveness of UKPDS interventions.* (UKPDS 62) Med Decis Making 2002; 22:340-9.
- Janssen LMM et al. *Health-related quality of life in type 2 diabetes mellitus: a pooled meta-analysis.* (utility values used here) 2022.
- Nicholson G et al. *Direct cost burden of major adverse cardiovascular events in commercially-insured US patients.* 2016.
- Brinkert M et al. *Long-term survival after MI in patients with diabetes.* Am J Cardiol 2017.
- Red Book / JAHA 2024 SGLT2i cost-effectiveness review (drug acquisition cost).
- MEPS / medRxiv 2025 cost-of-illness review (post-event ongoing cost).
