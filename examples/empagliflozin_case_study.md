# Empagliflozin vs Standard of Care in T2D with established CVD

This illustrative cost-effectiveness analysis runs on the validated OpenCEA
engine. The model inputs are documented in
[`examples/empagliflozin_t2d.yaml`](empagliflozin_t2d.yaml) with a source or
modeling rationale. The $4,500 net-price scenario is an explicit modeling
assumption, not an observed price.

This writeup follows the CHEERS 2022 reporting structure adapted for a brief working document.

---

## Decision problem

Is empagliflozin added to standard of care cost-effective versus standard of care alone for adults with type 2 diabetes (T2D) and established cardiovascular disease (CVD) at a willingness-to-pay (WTP) of **$100,000 / QALY**, and how robust is that conclusion to plausible parameter uncertainty?

| Element | Choice |
|---|---|
| Comparators | Standard of care (SoC); Empagliflozin + SoC |
| Population | Adults with T2D + established CVD, mean age 63 (EMPA-REG OUTCOME placebo arm) |
| Perspective | US healthcare-sector (drug + medical care; non-medical costs excluded) |
| Time horizon | Lifetime: start age 63, terminate at age 100 (37 annual cycles) |
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

- **EF:** event-free, with T2D and CVD but no major CV event yet.
- **PE:** post major CV event (a composite of MI, stroke, and HF hospitalization).
- **D:** dead, absorbing.

PE cannot recover to EF, matching the "Sicker" structure in the DARTH
tutorial. Transition probabilities are built from annual rates with
`p = 1 - exp(-r * t)`. The engine scales non-death transitions by
`(1 - p_death)` so each row sums to 1.

The acute event cost on entering PE is implemented as a **transition cost**: discounted at the engine's `dw_c[t+1] * wcc[t+1]` weight, matching how state costs are discounted at the cycle when occupants are present.

---

## Two-channel treatment effect

Empagliflozin acts on the model through two independent channels, each anchored on a different EMPA-REG OUTCOME endpoint:

1. **Progression channel (`HR_event = 0.86`):** composite hazard ratio for first major CV event (`EF -> PE`). Constructed as the EMPA-REG composition-weighted blend `0.60 * 1.00 (MI/stroke, no significant trial effect) + 0.40 * 0.65 (HF hospitalization HR 0.65 [0.50-0.85])`. 95% CI for the composite is approximately `[0.74, 0.99]`.
2. **Mortality channel (`HR_death = 0.68`):** all-cause mortality HR (95% CI `[0.57, 0.82]`, Zinman 2015). The model applies it to both `EF -> D` and `PE -> D`, assuming that the trial-reported all-cause mortality benefit persists in both health states. This assumption has the largest economic leverage in the tornado analysis. Applying the HR only in EF would materially change the model and requires a separate scenario.

The model applies no separate utility effect. Empagliflozin's QALY benefit
comes through event avoidance and survival.

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
| Standard of care | 259,540 | 10.944 | Not applicable | Not applicable | Not applicable | 834,855 |
| Empagliflozin | 397,423 | 12.338 | 137,883 | 1.394 | **98,900** | 836,388 |

Incremental NMB of Empa vs SoC at $100k / QALY: **$1,534** (very thin margin).

The acute event cost contributes $4,449 to SoC and $4,325 to Empa (Empa's lower EF→PE transition rate slightly reduces acute event spend, but the difference is tiny relative to drug cost).

---

## One-way DSA at $100k / QALY

Top 5 drivers, sorted by swing (full table in `tornado.png`):

| Rank | Parameter | Low → High value | Inc NMB at low | Inc NMB at high | Swing (USD) |
|---:|---|---|---:|---:|---:|
| 1 | `hr_death` | 0.57 → 0.82 | +35,815 | −33,080 | **68,895** |
| 2 | `c_drug` | 4,872 → 7,596 | +24,501 | −21,434 | 45,935 |
| 3 | `hr_event` | 0.74 → 0.99 | +14,891 | −12,301 | 27,192 |
| 4 | `u_EF` | 0.66 → 0.83 | −10,632 | +13,329 | 23,961 |
| 5 | `r_PE_D` | 0.04 → 0.10 | +12,041 | −7,943 | 19,984 |

`hr_death` flips the decision because its bar straddles zero. Empa goes from
cost-effective at the optimistic CI bound to not cost-effective at the
pessimistic CI bound. Drug cost is the second-largest single lever. Lowering
the price to $4,872 (FSS) makes Empa cost-effective at $100k WTP without any
other change. `c_acute_PE` has near-zero swing because the difference in
EF→PE flow between strategies is small.

`hr_event` ranks third even though the composite HR is closer to 1 than `hr_death`. That is consistent with empagliflozin's CV benefit in EMPA-REG being driven more by mortality than by non-fatal events.

---

## PSA, CEAC, and CE plane (n_sim = 10,000)

| WTP ($/QALY) | P(Empa cost-effective) |
|---:|---:|
| 50,000 | 0.002 |
| **100,000** | **0.517** |
| 150,000 | 0.945 |

The deterministic ICER ($98,900) is just below $100k WTP, and the CEAC
crosses 50% near that threshold. The result sits on the cost-effectiveness
boundary rather than safely on one side. The PSA mean ICER (cost difference
divided by QALY difference in expectation) is within about 3% of the
deterministic value.

The CE plane scatter (`empa_ce_plane.png`) lies almost entirely in the
upper-right quadrant. Empa is usually more costly and more effective than SoC;
the decision depends on whether the QALY gain justifies the added cost.

---

## Scenario analysis

The base case sits *right at* the $100k/QALY threshold, so the conclusion turns on two assumptions flagged in the limitations: WAC drug pricing and lifetime extrapolation of the EMPA-REG mortality benefit. Both deserve explicit scenarios.

### Scenarios

- **Net price:** annual drug cost lowered to **$4,500/yr**, an illustrative
  ~28% rebate from the $6,264 WAC. Real net prices are confidential and vary
  by payer; this is a round modeling assumption, not a quotation.
- **Treatment-effect waning:** both `hr_event` and `hr_death` remain at their
  trial point estimates for the first **3 years** (EMPA-REG follow-up), then
  move linearly to **HR = 1.0 by year 10**. The model applies no treatment
  effect for the remaining 27 cycles. A per-cycle transition-matrix sequence
  implements the scenario, and a regression test confirms that identical
  matrices reduce to the validated `evaluate_strategy` result.

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

The two assumptions move in opposite directions, but not by the same amount.
The net-price scenario reduces incremental cost by roughly the discounted
lifetime drug-cost difference (~$37k) and puts the ICER 21% below the base
case. Under waning, the incremental QALY gain falls from 1.39 to about 0.49,
a reduction of roughly 65%. Cost savings from fewer high-cost PE years offset
part of that loss, but the waning ICER still more than doubles.

Under the sustained-effect assumption, the deterministic result is
cost-effective at $100k/QALY at WAC, but only 51.7% of PSA draws favor
empagliflozin. The $4,500 net-price scenario raises that probability to 92.0%.
If the mortality benefit fades within about 10 years, the drug is not
cost-effective at $100k/QALY even after the modeled ~28% rebate. Its annual
price would need to fall below roughly $2,650 to clear the threshold.

The modeled conclusion is most sensitive to the durability of the EMPA-REG
all-cause mortality benefit. Drug price still matters and is the second-largest
driver in the one-way sensitivity analysis.

---

## Limitations

This is an **illustrative** model. The following limitations are material to
the result:

1. **Aggregated event state.** MI, stroke, and HF are collapsed into a single "PE" state with a blended ongoing cost, acute cost, and utility. Published US empagliflozin CEAs in the **$26k - $88k / QALY** range typically separate these (often into 4-6 states), which yields a tighter cost ledger and different mortality trajectories for each event type.
2. **Cohort, not microsimulation.** No individual heterogeneity, no time since first event, no second-event modeling. The 60/40 MI-stroke vs HF blend stays fixed at the trial composition, although the mix could change with age and prior events.
3. **Lifetime extrapolation of trial-derived HRs.** EMPA-REG OUTCOME measured outcomes over a median 3.1 years. The base case applies both HRs (`hr_event` = 0.86, `hr_death` = 0.68) for all 37 cycles. The **waning scenario** quantifies the impact of that assumption: the deterministic ICER more than doubles, from $98,900 to $206,778/QALY, when the effect tapers to zero by year 10.
4. **All-cause mortality HR applied to PE→D.** The trial measured HR 0.68 in the trial population, which started event-free. Applying that HR to post-event mortality is a strong assumption. The DSA shows it is the largest driver of incremental NMB, with a tornado bar that straddles zero across the trial's 95% CI.
5. **WAC drug price.** $6,264/yr is the listed wholesale price. Most US payers see meaningfully lower **net** prices after rebates; the net-price scenario above shows that a $4,500/yr price (~28% rebate) lowers the ICER to $77,564/QALY under sustained effect. Real net prices are confidential.
6. **No adverse-event modeling.** Genital infections (the most notable EMPA-REG side effect) and DKA risk are not separately costed.
7. **Linear time-homogeneity in baseline rates.** SoC rates do not change with age, time since diagnosis, or cumulative drug exposure. The Empa transition matrix can vary by cycle through the additive sequence engine used for waning.

## Bottom line

The base-case ICER ($98,900/QALY) lands at the upper end of plausibility, just below a $100k WTP. The scenario grid makes the driver of that conclusion explicit:

- **If the EMPA-REG mortality benefit is sustained**, the deterministic ICER
  is $98,900/QALY at WAC pricing. The probability of cost-effectiveness is
  51.7% at WAC and 92.0% at the modeled $4,500 net price.
- **If the effect wanes** (linearly to zero between years 3 and 10), the ICER rises above $200k/QALY at WAC and stays above $150k/QALY even at a $4,500 net price.

Within this scenario grid, the largest uncertainty is the durability of the
all-cause mortality benefit. Drug pricing remains important. Two structural
choices distinguish this analysis from the central published $26-88k/QALY
range: WAC drug pricing, which the net-price scenario relaxes, and the
three-state aggregation, which remains an illustrative simplification.

---

## Figures

Generated by the OpenCEA plotting module from the same `PSAResult` / `DSAResult` objects the tests verify. Saved to `examples/figures/`:

- `empa_tornado.png`: one-way DSA tornado on incremental NMB at $100k.
- `empa_ceac.png`: CEAC for SoC vs Empa across $0 - $200k WTP.
- `empa_ce_plane.png`: incremental cost / QALY scatter, 10,000 draws.
- `empa_ce_frontier.png`: expected NMB vs WTP, with the optimal-strategy switch marked.

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
