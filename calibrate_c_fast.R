## ============================================================================
## Step 2 (FAST PASS): calibrate the Go threshold c under scenario 1 (all null).
##
## Goal of the fast pass (NOT report-grade numbers):
##   (1) confirm the pipeline runs end to end without crashing,
##   (2) check the all-null datasets don't blow up HMC (divergences),
##   (3) confirm per-basket false-positive rate decreases monotonically in c
##       and crosses the 0.10 target -> c* exists, and locate it roughly,
##       so the FULL run can use a tight c-grid and a large M.
##
## Key efficiency point: pi_hat = Pr(p_k > p0 | y) depends only on the data, NOT
## on c. So we fit M times ONCE, cache the M x K matrix of pi_hat, and the search
## over c is just thresholding the cache (free). Cost = M fits only.
##
## Needs bhm.stan in the same directory. Engine: cmdstanr.
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- locked design ----
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets)
n  <- rep(20L, K)
p0 <- 0.10

## scenario 1: all baskets truly null
p_true_scen1 <- rep(0.10, K)

## ---- hyperpriors (fixed during calibration; only c is calibrated) ----
tau_scale <- 1.0
mu_mean   <- logit(p0)   # -2.197
mu_sd     <- 2.0

## ---- FAST-PASS knobs (deliberately cheap; bump these for the full run) ----
M            <- 300L     # number of simulated trials  (full: 2000-5000)
n_chains     <- 2L       # full: 4
warmup       <- 500L     # keep >= 500: all-null data is near the boundary
sampling     <- 500L     # per chain; full: 2000+
adapt_delta  <- 0.95     # full: 0.99
master_seed  <- 2026L

mod <- cmdstan_model("bhm.stan")

## ---- simulate M datasets and fit once each; CACHE pi_hat ----
set.seed(master_seed)
pi_hat   <- matrix(NA_real_, nrow = M, ncol = K,
                   dimnames = list(NULL, baskets))   # Pr(p_k > p0 | y) per sim
n_diverg <- integer(M)                                # divergences per fit
y_store  <- matrix(NA_integer_, nrow = M, ncol = K)   # keep the data too

t0 <- Sys.time()
for (m in seq_len(M)) {
  y_m <- rbinom(K, n, p_true_scen1)
  y_store[m, ] <- y_m
  
  fit <- mod$sample(
    data = list(K = K, n = n, y = y_m,
                mu_mean = mu_mean, mu_sd = mu_sd,
                tau_scale = tau_scale, p0 = p0),
    seed = master_seed + m,
    chains = n_chains, parallel_chains = n_chains,
    iter_warmup = warmup, iter_sampling = sampling,
    adapt_delta = adapt_delta, max_treedepth = 12,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  
  ## pi_hat_k = posterior mean of the indicator {p_k > p0}
  pi_hat[m, ] <- colMeans(fit$draws("exceeds_p0", format = "draws_matrix"))
  n_diverg[m] <- tryCatch(
    sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]),
    error = function(e) NA_integer_)
  
  if (m %% 50 == 0) cat(sprintf("  fit %d / %d\n", m, M))
}
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## cache to disk so the c-search can be re-run without refitting
saveRDS(list(pi_hat = pi_hat, y = y_store, n_diverg = n_diverg,
             hyper = list(tau_scale = tau_scale, mu_mean = mu_mean, mu_sd = mu_sd)),
        file = "pihat_scenario1_fast.rds")

## ---- sanity: did the all-null data behave for HMC? ----
cat(sprintf("\nfits with >=1 divergence: %d / %d (%.1f%%); max divergences in a fit: %d\n",
            sum(n_diverg > 0, na.rm = TRUE), M,
            100 * mean(n_diverg > 0, na.rm = TRUE), max(n_diverg, na.rm = TRUE)))
## If this is more than a few %, bump adapt_delta (and warmup) before trusting numbers.

## ---- search c over a grid (FREE: just thresholding the cache) ----
c_grid <- seq(0.80, 0.99, by = 0.005)

## per-basket false-positive rate at each c: alpha_k(c) = mean_m 1{pi_hat[m,k] > c}
alpha_by_c <- sapply(c_grid, function(cc) colMeans(pi_hat > cc))  # K x length(c_grid)
rownames(alpha_by_c) <- baskets
mean_alpha <- colMeans(alpha_by_c)   # averaged over the 4 (exchangeable) baskets

tab <- data.frame(c = c_grid, mean_alpha = mean_alpha, t(alpha_by_c))
cat("\nper-basket false-positive rate vs c (scenario 1):\n")
print(tab, digits = 3, row.names = FALSE)

## monotonic? (mean_alpha should be non-increasing in c)
cat(sprintf("\nmean_alpha monotone non-increasing in c: %s\n",
            all(diff(mean_alpha) <= 1e-9)))

## ---- locate c* where mean per-basket FP rate crosses 0.10 ----
target <- 0.10
if (min(mean_alpha) > target || max(mean_alpha) < target) {
  cat(sprintf("\n0.10 NOT bracketed by c in [%.2f, %.2f] (range of mean_alpha: %.3f-%.3f).\n",
              min(c_grid), max(c_grid), min(mean_alpha), max(mean_alpha)))
  cat("=> widen c_grid before trusting c*.\n")
} else {
  ## mean_alpha is decreasing in c; invert by interpolation
  c_star <- approx(x = rev(mean_alpha), y = rev(c_grid), xout = target)$y
  ## MC error of a per-basket rate at alpha=0.10 with M sims, for context:
  se_alpha <- sqrt(target * (1 - target) / M)
  cat(sprintf("\nc* (mean per-basket FP = 0.10)  ~  %.3f\n", c_star))
  cat(sprintf("per-basket alpha at nearest grid c: \n"))
  k_near <- which.min(abs(c_grid - c_star))
  print(round(alpha_by_c[, k_near], 3))
  cat(sprintf("\n[fast-pass MC note] SE of a per-basket rate ~ %.3f at M=%d; the 4\n",
              se_alpha, M))
  cat("per-basket alphas should agree within a couple SEs (they are exchangeable\n")
  cat("under scenario 1). This c* is a LOCATION estimate only -- re-run FULL\n")
  cat("(large M, 4 chains, adapt_delta=0.99, tight c-grid around it) for the value\n")
  cat("you carry into the power runs.\n")
}