# Bayesian Information Borrowing in Basket Trials

A from-scratch **Bayesian trial-design and simulation study** of information
borrowing across tumor subtypes in a basket trial, motivated by VE-BASKET
(vemurafenib in BRAF V600–mutant cancers). The project quantifies, by
simulation, the central trade-off of borrowing: it sharpens estimates and
buys power when baskets are truly similar, but drags truly-inactive baskets
upward and inflates their false-positive rate when baskets differ.

> **Thesis.** Small per-basket samples force borrowing; borrowing is a *risky
> knob*; homogeneous truth → borrowing pays (power ↑), heterogeneous truth →
> blind borrowing inflates per-basket type I error. We measure the trade-off
> with operating characteristics (OC) and show how calibration controls it.

This is a methods/portfolio project, not a regulatory deliverable. Everything
is reproducible from the scripts below.

---

## Trial design (locked specification)

- **Endpoint:** objective response rate (ORR), binary. `y_k ~ Binomial(n_k, p_k)`.
- **Hypotheses (per basket):** `H0: p_k ≤ 0.10` vs `H1: p_k > 0.10`. Clinically
  meaningful rate `p1 = 0.30`.
- **Baskets:** `K = 4`, uniform `n_k = 20` (unequal-n left as future work).
  True anchors taken from VE-BASKET (Hyman et al. 2015, Table 2):

  | Basket | n | Real ORR anchor | Role |
  |---|---|---|---|
  | NSCLC | 20 | 42% | truly high-activity |
  | Anaplastic thyroid | 20 | 29% | mid-activity |
  | Biliary | 20 | 12% | low / gray |
  | Colorectal (CRC), monotherapy | 20 | 0% | truly inactive |

- **Simulation scenarios (true ORRs):**

  | Scenario | NSCLC | Thyroid | Biliary | CRC | Purpose |
  |---|---|---|---|---|---|
  | 1. global null | 0.10 | 0.10 | 0.10 | 0.10 | calibrate `c`; per-basket type I error |
  | 2. all active | 0.30 | 0.30 | 0.30 | 0.30 | power; homogeneous (borrowing should pay) |
  | 3. mixed | 0.30 | 0.30 | 0.10 | 0.10 | borrowing across a mixed structure |
  | 4. heterogeneous | 0.40 | 0.30 | 0.15 | 0.05 | flagship; inactive CRC next to high-activity NSCLC |

  Scenario-4 biliary (0.15) sits in the gray zone (`p0 < p < p1`) and is
  reported separately, not as power or type I error.

- **Methods compared (columns of the OC table):**
  1. **No borrowing** — independent Beta-Binomial (lower bracket).
  2. **BHM** — Bayesian hierarchical model (main method, partial borrowing).
  3. **Complete pooling** — single shared rate (upper bracket / strawman).
  4. **EXNEX** — robust mixture (planned; column 4).

- **Decision rule:** `Pr(p_k > p0 | data) > c  ⇒  Go`, else No-Go.
  `c` is **not fixed at 0.95**; it is calibrated.

- **Calibration:** each method is calibrated **separately** so that the
  **per-basket type I error = 0.10** under scenario 1 (loose, single-sided,
  appropriate for a phase-II go/no-go screen — not the 0.025 of a confirmatory
  trial). Power is then read off scenarios 2–4 at this matched type I error.

- **Reporting rule:** type I error is reported **per basket**, never
  family-wise. The same `Pr(Go)` is labeled *power* when the basket is truly
  active and *type I error* when it is truly null.

### BHM specification

`logit(p_k) = θ_k`, `θ_k ~ Normal(μ, τ²)`, `μ ~ Normal(m0, s0²)`,
`τ ~ Half-Normal(0, s)`. Implemented **non-centered** (`θ_k = μ + τ·z_k`,
`z_k ~ N(0,1)`) to avoid the funnel/divergences that bite HMC with few baskets
and small τ. Likelihood uses `binomial_logit` for numerical stability.

---

## Repository structure

```
bayesian-basket-borrowing/
├── bhm.stan                      # BHM: logit-normal, half-normal τ, non-centered
├── no_borrow.stan                # independent Beta-Binomial per basket
├── pool.stan                     # complete pooling: one shared p
├── run_bhm_step1.R               # demo: single-scenario BHM fit + decision (step 1)
├── calibrate_c_fast.R            # FAST pass: locate BHM c* (M=300) — optional
├── calibrate_c_full.R            # BHM c* calibration under scenario 1 (M=2000)
├── step3a_calibrate_brackets.R   # calibrate c* for no-borrow & complete-pool
├── step3b_oc_sweep.R             # OC sweep over scenarios 2/3/4, all methods
├── basket_notes.tex              # companion study/derivation notes (Chinese; xelatex)
└── README.md
```

Generated artifacts (git-ignored): `pihat_scenario1_full.rds`, `c_star.rds`,
`c_star_all.rds`, `oc_results.rds`, and CmdStan-compiled model executables.

---

## Reproducing the results

**Dependencies**

- R (tested on 4.4.x) with a C++ toolchain (Rtools on Windows).
- [`cmdstanr`](https://mc-stan.org/cmdstanr/) + CmdStan
  (`cmdstanr::install_cmdstan()`).
- `future` + `furrr` for parallelism (optional; scripts fall back to serial).

**Run order** (each step consumes the previous step's `.rds`):

```r
# 1. BHM threshold calibration under the global null
source("calibrate_c_full.R")        # -> pihat_scenario1_full.rds, c_star.rds

# 2. Calibrate the two bracket methods on the SAME null datasets
source("step3a_calibrate_brackets.R")   # -> c_star_all.rds

# 3. Operating-characteristics sweep (scenarios 2/3/4, all 3 methods)
source("step3b_oc_sweep.R")         # -> oc_results.rds
```

`calibrate_c_fast.R` and `run_bhm_step1.R` are optional warm-ups and are not
required by the pipeline. Set the working directory to the repo root so the
relative paths to `*.stan` resolve.

---

## Results so far

Calibrated thresholds (per-basket type I error = 0.10 under scenario 1,
`M = 2000`, per-basket MC SE ≈ 0.0067):

| Method | calibrated `c*` |
|---|---|
| No borrowing | 0.920 |
| BHM | 0.820 |
| Complete pooling | 0.914 |

BHM's threshold is materially lower because borrowing shrinks the per-basket
posteriors toward a common mean, compressing the null distribution of
`Pr(p_k > p0 | data)`. This is exactly why each method must be calibrated
separately: reusing one method's `c` for another would mismatch the type I
error and confound calibration with method. (At the conventional `c = 0.95`,
BHM's per-basket type I error would be ≈ 1.3%, wasting ~7/8 of the 0.10 budget
as lost power.)

HMC diagnostics under the (funnel-prone) global null are clean: ≈ 0.3% of fits
have any divergence, ≈ 0.0001% of transitions overall (non-centered + high
`adapt_delta`).

---

## Methodological notes & open items

These are deliberate, flagged choices — good targets for review:

- **`μ` prior (`mu_mean`, `mu_sd`)** is a placeholder weakly-informative choice
  (`m0 = logit(0.10)`, `s0 = 2`), **not** yet locked to a basket-literature
  convention.
- **`τ` prior scale `s = 1`** (logit scale) is the main-analysis setting; the
  prior-sensitivity sweep (`s`, half-Cauchy, inverse-gamma) is planned and is
  itself a selling point given τ is weakly identified at `K = 4`.
- **Half-Normal attribution.** `τ ~ Half-Normal(0, s)` is attributed to the
  basket-borrowing literature / Neuenschwander et al. (2016). It is **not**
  attributed to Gelman (2006), who recommends half-Cauchy / folded-t; Gelman
  (2006) is cited only to justify *avoiding* inverse-gamma priors with few
  groups.
- **Berry (2013) τ² prior form** is referenced as the canonical BHM but its
  exact hyperprior is to be verified against the source.
- **No-borrow prior** is independent Jeffreys `Beta(0.5, 0.5)`. The prior forms
  differ across methods, but each is calibrated to the same type I error, so the
  decision-level comparison is fair regardless of prior shape.
- **Variance reduction.** Within each scenario all three methods are fit on the
  **same** simulated datasets (paired); calibration reuses the exact scenario-1
  datasets across methods.
- **Hand-rolled Stan vs packages.** `bhmbasket` (JAGS) covers berry/exnex/
  pooled/stratified and is intended only as an independent cross-check; the main
  pipeline is hand-rolled Stan because the τ-prior *family* sensitivity (the
  project's selling point) requires controlling the prior, which the packaged
  models hard-code.

---

## Notes for reviewers

Specific things worth checking:

1. **Stan models.** Is the non-centered BHM equivalent to the centered spec? Is
   the half-normal prior correctly induced by `<lower=0>` + `normal(0, s)`? Is
   `binomial_logit` used correctly? Does `pool.stan` correctly impose a single
   shared `p`, and is its single decision correctly broadcast to all baskets in R?
2. **Calibration logic** (`calibrate_c_*`, `step3a`). Is the per-basket
   false-positive curve monotone in `c`, and is the bracketed interpolation for
   `c*` correct at the boundary? Is targeting **per-basket** (not family-wise)
   type I error the intended and correctly-implemented criterion?
3. **Decision quantity.** `Pr(p_k > p0 | data)` is computed as the posterior
   mean of an indicator (`exceeds_p0`). Is the power/type-I-error labeling by
   true ORR (`role_of`) correct, including the gray-zone handling?
4. **Pairing / seeding.** Are datasets correctly shared across methods within a
   scenario, and is reproducibility independent of parallel worker scheduling?
5. **Estimation metrics.** Are bias / MSE / posterior SD computed against the
   correct per-basket truth (note pooling's single estimate is compared to each
   basket's own truth on purpose)?
6. **Statistical critique welcome.** Is `n = 20`, `p0 = 0.10`, `p1 = 0.30`,
   target type I error 0.10 a coherent phase-II screening setup? Are there
   borrowing pathologies or calibration edge cases not yet covered?

---

## Status / roadmap

- [x] Step 1 — BHM in Stan, single-scenario posterior + decision rule
- [x] Step 2 — calibrate BHM `c*` under scenario 1 (per-basket type I error 0.10)
- [x] Step 3a — calibrate `c*` for no-borrow and complete-pool
- [ ] Step 3b — OC sweep (scenarios 2/3/4): power + per-basket type I error +
      bias/MSE/posterior SD (ready to run)
- [ ] Step 4 — add EXNEX (column 4); focus on scenario-4 CRC false-positive
      protection vs BHM
- [ ] Step 5 — prior-sensitivity analysis (`s`, half-Cauchy, inverse-gamma; π for EXNEX)
- [ ] Step 6 — single interim futility analysis (note: borrowing makes interim
      decisions across baskets non-independent)
- [ ] Optional — SAS (PROC MCMC / BGLIMM) cross-check of core posteriors

---

## References

- Hyman DM, et al. *Vemurafenib in multiple nonmelanoma cancers with BRAF V600
  mutations.* N Engl J Med. 2015;373:726–736. (VE-BASKET)
- Berry SM, Broglio KR, Groshen S, Berry DA. *Bayesian hierarchical modeling of
  patient subpopulations: efficient designs of phase II oncology clinical
  trials.* Clin Trials. 2013;10(5):720–734.
- Neuenschwander B, Wandel S, Roychoudhury S, Bailey S. *Robust exchangeability
  designs for early phase clinical trials with multiple strata.* Pharm Stat.
  2016;15(2):123–134. (EXNEX)
- Gelman A. *Prior distributions for variance parameters in hierarchical models.*
  Bayesian Anal. 2006;1(3):515–534.
