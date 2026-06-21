## ============================================================================
## Step 3a-EXNEX: calibrate the Go threshold c for EXNEX (column 4) under
##                scenario 1 (all null), then APPEND it to c_star_all.rds.
##
## Same logic as step3a_calibrate_brackets.R: calibrate EXNEX SEPARATELY to
## per-basket type I error = 0.10 so the later power comparison is at matched
## type I error. Reuses the EXACT scenario-1 datasets cached in
## pihat_scenario1_full.rds (the same Y that calibrated BHM / no-borrow / pool),
## so all FOUR methods are calibrated on identical null data (paired).
##
## Needs: exnex.stan, pihat_scenario1_full.rds, c_star_all.rds  (in this dir).
## Writes: c_star_all.rds  (adds $exnex and $hyper_exnex; keeps the other three).
## Engine: cmdstanr.
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- locked design ----
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets)
n  <- rep(20L, K)
p0 <- 0.10
target <- 0.10

## ---- EX hyperprior: reuse EXACTLY what BHM was calibrated under ----
cs <- readRDS("c_star_all.rds")
mu_mean   <- cs$hyper_bhm$mu_mean
mu_sd     <- cs$hyper_bhm$mu_sd
tau_scale <- cs$hyper_bhm$tau_scale

## ---- EXNEX design constants (main analysis) ----
m_nex  <- 0.0                       # NEX mean, logit scale  (ORR 0.5 center)
s_nex  <- 2.0                       # NEX sd,   logit scale  (+/-2SD ~ ORR 0.02-0.98)
pi_val <- 0.5                       # fixed EX prior weight, all baskets (main analysis)
pi_ex  <- rep(pi_val, K)

## ---- reuse the SAME scenario-1 datasets used to calibrate the other 3 ----
cache1 <- readRDS("pihat_scenario1_full.rds")
Y <- cache1$Y                       # M x K matrix of null datasets
M <- nrow(Y)
cat(sprintf("reusing %d scenario-1 datasets (paired with BHM/no-borrow/pool)\n", M))

## ---- fit knobs (EXNEX keeps the EX-component funnel -> BHM-grade adapt_delta) ----
n_chains    <- 4L
warmup      <- 1000L
sampling    <- 1000L
adapt_delta <- 0.99
max_td      <- 12L
master_seed <- 2026L
seeds <- master_seed + seq_len(M)

PARALLEL  <- TRUE
n_workers <- max(1L, parallel::detectCores() - 1L)
cmdstan_dir <- cmdstanr::cmdstan_path()

## compile once in main; workers just load the exe
stan_file <- normalizePath("exnex.stan")
cmdstan_model(stan_file)

## ---- one fit -> per-basket pi_hat vector (length K) ----
fit_pi <- function(m) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = list(K = K, n = n, y = Y[m, ], p0 = p0,
                mu_mean = mu_mean, mu_sd = mu_sd, tau_scale = tau_scale,
                m_nex = m_nex, s_nex = s_nex, pi_ex = pi_ex),
    seed = seeds[m], chains = n_chains, parallel_chains = 1,
    iter_warmup = warmup, iter_sampling = sampling,
    adapt_delta = adapt_delta, max_treedepth = max_td,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  list(pi   = colMeans(fit$draws("exceeds_p0", format = "draws_matrix")),
       ndiv = sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]))
}

## ---- run over all M datasets (parallel or serial fallback) ----
cat(sprintf("\n[EXNEX] fitting %d datasets ...\n", M))
t0 <- Sys.time()
if (PARALLEL && requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)) {
  library(furrr); library(future)
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)
  res <- furrr::future_map(seq_len(M), fit_pi,
           .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
  plan(sequential)
} else {
  if (PARALLEL) message("furrr/future not found -> SERIAL.")
  res <- lapply(seq_len(M), function(m) {
    if (m %% 200 == 0) cat(sprintf("  %d / %d\n", m, M)); fit_pi(m) })
}
cat(sprintf("[EXNEX] done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

pi_hat   <- t(vapply(res, `[[`, numeric(K), "pi")); colnames(pi_hat) <- baskets
n_diverg <- vapply(res, `[[`, numeric(1), "ndiv")
total_draws <- n_chains * sampling
cat(sprintf("divergences: fits with >=1 = %d/%d (%.1f%%); overall rate = %.4f%% of draws\n",
            sum(n_diverg > 0), M, 100 * mean(n_diverg > 0),
            100 * sum(n_diverg) / (M * total_draws)))

## ---- robust c*: explicit bracketed interpolation on mean per-basket FP ----
## (identical procedure to step3a_calibrate_brackets.R)
calibrate <- function(pi_hat, c_grid) {
  alpha_by_c <- sapply(c_grid, function(cc) colMeans(pi_hat > cc))
  rownames(alpha_by_c) <- baskets
  mean_alpha <- colMeans(alpha_by_c)
  monotone <- all(diff(mean_alpha) <= 1e-9)
  if (min(mean_alpha) > target || max(mean_alpha) < target)
    return(list(c_star = NA_real_, bracketed = FALSE, boundary_hit = NA,
                alpha_by_c = alpha_by_c, mean_alpha = mean_alpha, monotone = monotone))
  i_hi <- max(which(mean_alpha >= target))
  if (i_hi == length(c_grid)) {
    c_star <- c_grid[i_hi]; boundary_hit <- TRUE
  } else {
    a1 <- mean_alpha[i_hi]; a2 <- mean_alpha[i_hi + 1]
    c_star <- c_grid[i_hi] + (target - a1) * (c_grid[i_hi + 1] - c_grid[i_hi]) / (a2 - a1)
    boundary_hit <- FALSE
  }
  list(c_star = c_star, bracketed = TRUE, boundary_hit = boundary_hit,
       alpha_by_c = alpha_by_c, mean_alpha = mean_alpha, monotone = monotone)
}

c_grid <- seq(0.40, 0.99, by = 0.0025)   # wide: EXNEX c* expected between BHM and no-borrow
cal <- calibrate(pi_hat, c_grid)

## ---- report ----
se_alpha <- sqrt(target * (1 - target) / M)
cat(sprintf("\n=== EXNEX (pi=%.2f, NEX N(%.1f,%.1f^2)) ===\n", pi_val, m_nex, s_nex))
cat(sprintf("  monotone non-increasing: %s\n", cal$monotone))
cat(sprintf("  c* (per-basket FP = %.2f) = %.4f\n", target, cal$c_star))
if (!is.na(cal$c_star)) {
  if (isTRUE(cal$boundary_hit))
    cat("  *** WARNING: c* at grid boundary -- extend c_grid and re-run.\n")
  k <- which.min(abs(c_grid - cal$c_star))
  cat("  per-basket alpha at c*: "); print(round(cal$alpha_by_c[, k], 4))
} else cat("  -> 0.10 NOT bracketed; widen c_grid.\n")
cat(sprintf("  [per-basket MC SE ~ %.4f at M=%d]\n", se_alpha, M))
cat(sprintf("\nfor reference: no-borrow %.4f | BHM %.4f | pool %.4f\n",
            cs$no_borrow, cs$bhm, cs$pool))

## ---- APPEND exnex to c_star_all.rds (keep the existing three untouched) ----
cs$exnex <- cal$c_star
cs$boundary_hit$exnex <- cal$boundary_hit
cs$hyper_exnex <- list(mu_mean = mu_mean, mu_sd = mu_sd, tau_scale = tau_scale,
                       m_nex = m_nex, s_nex = s_nex, pi_val = pi_val)
saveRDS(cs, file = "c_star_all.rds")
cat("\nsaved c_star_all.rds with $exnex added (bhm / no_borrow / pool / exnex).\n")
