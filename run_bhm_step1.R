## ============================================================================
## Step 1 (per handoff §7): single-scenario BHM fit + decision rule.
## Engine: cmdstanr (hand-rolled Stan spine). bhmbasket/JAGS is the later
## independent cross-check, not this.
##
## Setup once:
##   install.packages("cmdstanr",
##     repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
##   cmdstanr::install_cmdstan()
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- locked design constants ----
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets)
n  <- rep(20L, K)
p0 <- 0.10                              # null / historical ORR

## ---- example data: ONE simulated trial under scenario 4 (0.40/0.30/0.15/0.05) ----
## Hardcoded to match the reference table. To draw a fresh trial instead:
##   set.seed(20240601); y <- rbinom(K, n, c(0.40, 0.30, 0.15, 0.05))
y <- c(10L, 6L, 3L, 3L)

## ---- hyperpriors (NOT yet locked -- see checklist items 4 & 5) ----
## tau ~ Half-Normal(0, s), s on logit scale. s = 1 already encodes
## substantial heterogeneity; this is THE knob the sensitivity analysis sweeps.
tau_scale <- 1.0
## mu ~ N(mu_mean, mu_sd) on logit scale. Placeholder weakly-informative prior
## centered near logit(p0). Neutral alternative: mu_mean = 0.
mu_mean <- logit(p0)                    # = -2.197
mu_sd   <- 2.0

stan_data <- list(K = K, n = n, y = y,
                  mu_mean = mu_mean, mu_sd = mu_sd,
                  tau_scale = tau_scale, p0 = p0)

## ---- fit ----
mod <- cmdstan_model("bhm.stan")
fit <- mod$sample(
  data = stan_data, seed = 1,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000,
  adapt_delta = 0.99, max_treedepth = 12,   # tight: small K + small tau => funnel risk
  refresh = 0
)

## ---- diagnostics FIRST (handoff §5) ----
fit$cmdstan_diagnose()                       # divergences, E-BFMI, treedepth
print(fit$summary(c("mu", "tau")), n = Inf)  # want Rhat ~ 1.00, healthy ESS

## ---- per-basket posterior + decision ----
draws_p <- fit$draws("p",          format = "draws_matrix")  # iters x K, cols p[1..K]
pr_go   <- colMeans(fit$draws("exceeds_p0", format = "draws_matrix"))  # Pr(p_k > p0)

post <- data.frame(
  basket    = baskets, y = y, n = n,
  post_mean = apply(draws_p, 2, mean),
  q025      = apply(draws_p, 2, quantile, 0.025),
  q975      = apply(draws_p, 2, quantile, 0.975),
  Pr_gt_p0  = as.numeric(pr_go)
)

## Decision: Pr(p_k > p0 | data) > c  ->  Go.
## c is CALIBRATED in step 2 to hit per-basket type I error = 0.10 under
## scenario 1. Do NOT fix it at 0.95. The value below is a single-run
## illustration only.
c_threshold   <- 0.90
post$decision <- ifelse(post$Pr_gt_p0 > c_threshold, "Go", "No-Go")

print(post, digits = 3, row.names = FALSE)

## ---- what to expect on this dataset ----
## Borrowing pulls the four no-borrow means (0.50 / 0.31 / 0.167 / 0.167)
## toward the common mean; how hard depends on the tau posterior. Sanity check
## that NSCLC shrinks down a little and Biliary/CRC shrink up. The per-basket
## false-positive question (does CRC's Pr(Go) get inflated?) is answered by the
## OC simulation in steps 2-4 averaging over many trials with CRC ~ Binom(20, 0.05),
## NOT by this single realized draw (which happens to have CRC at 3/20).



