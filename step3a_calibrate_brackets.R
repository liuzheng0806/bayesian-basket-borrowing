## ============================================================================
## Step 3a: calibrate the Go threshold c for the two BRACKET methods
##          (no-borrow, complete-pool) under scenario 1 (all null).
##
## Each method is calibrated SEPARATELY to per-basket type I error = 0.10, so
## the later power comparison is at matched type I error (apples to apples).
## BHM's c* was already done (loaded from c_star.rds).
##
## Reuses the EXACT scenario-1 datasets cached in pihat_scenario1_full.rds, so
## all three methods are calibrated on identical null data (paired comparison).
##
## Needs no_borrow.stan, pool.stan, c_star.rds, pihat_scenario1_full.rds in dir.
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- locked design ----
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets)
n  <- rep(20L, K)
p0 <- 0.10
target <- 0.10

## ---- Beta prior for the bracket models (Jeffreys) ----
a0 <- 0.5
b0 <- 0.5

## ---- reuse the SAME scenario-1 datasets used to calibrate BHM ----
cache1 <- readRDS("pihat_scenario1_full.rds")
Y <- cache1$Y                      # M x K matrix of null datasets
M <- nrow(Y)
cat(sprintf("reusing %d scenario-1 datasets from BHM calibration\n", M))

## ---- fit knobs (these models are tiny/fast; no funnel) ----
n_chains    <- 4L
warmup      <- 1000L
sampling    <- 1000L
adapt_delta <- 0.95
master_seed <- 2026L
seeds <- master_seed + seq_len(M)

PARALLEL  <- TRUE
n_workers <- max(1L, parallel::detectCores() - 1L)
cmdstan_dir <- cmdstanr::cmdstan_path()

## compile both models once in main
nb_file   <- normalizePath("no_borrow.stan"); cmdstan_model(nb_file)
pool_file <- normalizePath("pool.stan");      cmdstan_model(pool_file)

## ---- one fit -> per-basket pi_hat vector (length K) ----
## shared = TRUE for pool (one decision replicated to all K baskets)
fit_pi <- function(m, stan_file, shared) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = list(K = K, n = n, y = Y[m, ], a0 = a0, b0 = b0, p0 = p0),
    seed = seeds[m], chains = n_chains, parallel_chains = 1,
    iter_warmup = warmup, iter_sampling = sampling, adapt_delta = adapt_delta,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  ex <- fit$draws("exceeds_p0", format = "draws_matrix")
  if (shared) rep(mean(ex[, 1]), K) else colMeans(ex)
}

## ---- run one method over all M datasets (parallel or serial) ----
run_method <- function(stan_file, shared, label) {
  cat(sprintf("\n[%s] fitting %d datasets ...\n", label, M))
  t0 <- Sys.time()
  if (PARALLEL && requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    library(furrr); library(future)
    plan(multisession, workers = n_workers)
    on.exit(plan(sequential), add = TRUE)
    pis <- furrr::future_map(seq_len(M), ~ fit_pi(.x, stan_file, shared),
                             .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
    plan(sequential)
  } else {
    if (PARALLEL) message("furrr/future not found -> SERIAL.")
    pis <- lapply(seq_len(M), function(m) {
      if (m %% 200 == 0) cat(sprintf("  %d / %d\n", m, M)); fit_pi(m, stan_file, shared) })
  }
  cat(sprintf("[%s] done in %.1f min\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  pi_hat <- do.call(rbind, pis); colnames(pi_hat) <- baskets
  pi_hat
}

## ---- robust c*: explicit bracketed interpolation on mean per-basket FP ----
calibrate <- function(pi_hat, c_grid) {
  alpha_by_c <- sapply(c_grid, function(cc) colMeans(pi_hat > cc))
  rownames(alpha_by_c) <- baskets
  mean_alpha <- colMeans(alpha_by_c)
  if (min(mean_alpha) > target || max(mean_alpha) < target)
    return(list(c_star = NA_real_, alpha_by_c = alpha_by_c, mean_alpha = mean_alpha,
                monotone = all(diff(mean_alpha) <= 1e-9)))
  i_hi <- max(which(mean_alpha >= target))
  if (i_hi == length(c_grid)) { c_star <- c_grid[i_hi]
  } else { a1 <- mean_alpha[i_hi]; a2 <- mean_alpha[i_hi+1]
  c_star <- c_grid[i_hi] + (target-a1)*(c_grid[i_hi+1]-c_grid[i_hi])/(a2-a1) }
  list(c_star = c_star, alpha_by_c = alpha_by_c, mean_alpha = mean_alpha,
       monotone = all(diff(mean_alpha) <= 1e-9))
}

c_grid <- seq(0.40, 0.99, by = 0.0025)   # wide: bracket c* might sit lower than BHM's

## ---- calibrate the two brackets ----
pi_nb   <- run_method(nb_file,   shared = FALSE, label = "no-borrow")
pi_pool <- run_method(pool_file, shared = TRUE,  label = "complete-pool")

cal_nb   <- calibrate(pi_nb,   c_grid)
cal_pool <- calibrate(pi_pool, c_grid)

## ---- load BHM c* (already calibrated) ----
c_bhm <- readRDS("c_star.rds")$c_star

## ---- report ----
se_alpha <- sqrt(target * (1 - target) / M)
report <- function(label, cal) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  monotone non-increasing: %s\n", cal$monotone))
  cat(sprintf("  c* (per-basket FP = %.2f) = %.4f\n", target, cal$c_star))
  if (!is.na(cal$c_star)) {
    k <- which.min(abs(c_grid - cal$c_star))
    cat("  per-basket alpha at c*: "); print(round(cal$alpha_by_c[, k], 4))
  } else cat("  -> 0.10 NOT bracketed; widen c_grid.\n")
}
report("no-borrow",     cal_nb)
report("complete-pool", cal_pool)
cat(sprintf("\n=== BHM (from c_star.rds) ===\n  c* = %.4f\n", c_bhm))
cat(sprintf("\n[per-basket MC SE ~ %.4f at M=%d]\n", se_alpha, M))

## ---- save all three c* for step 3b ----
saveRDS(list(bhm = c_bhm, no_borrow = cal_nb$c_star, pool = cal_pool$c_star,
             target = target, a0 = a0, b0 = b0,
             hyper_bhm = readRDS("c_star.rds")$hyper),
        file = "c_star_all.rds")
cat("\nsaved c_star_all.rds (bhm / no_borrow / pool) for step 3b.\n")