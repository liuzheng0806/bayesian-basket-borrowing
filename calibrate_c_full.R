## ============================================================================
## Step 2 (FULL RUN): calibrate the Go threshold c under scenario 1 (all null).
##
## Produces the report-grade c* you carry into the power runs (step 3).
## Changes vs the fast pass:
##   - M = 2000           (per-basket FP SE ~ 0.0067 at alpha=0.10)
##   - 4 chains, warmup 1000, sampling 1000, adapt_delta = 0.99
##   - c-grid extended DOWN to 0.70 (fast pass put c* at the 0.80 edge) and finer
##   - robust c* via explicit bracketed interpolation (no approx ties warning)
##   - PARALLEL across the M datasets (furrr/multisession; serial fallback)
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
p_true_scen1 <- rep(0.10, K)

## ---- hyperpriors (fixed during calibration) ----
tau_scale <- 1.0
mu_mean   <- logit(p0)
mu_sd     <- 2.0

## ---- FULL-RUN knobs ----
M            <- 2000L
n_chains     <- 4L
warmup       <- 1000L
sampling     <- 1000L      # 4 chains x 1000 = 4000 draws -> pi_hat very precise
adapt_delta  <- 0.99
master_seed  <- 2026L
target       <- 0.10       # per-basket type I error target

## ---- parallel settings ----
PARALLEL  <- TRUE
n_workers <- max(1L, parallel::detectCores() - 1L)

## compile once in the main session so workers just load the existing exe
stan_file   <- normalizePath("bhm.stan")
mod_main    <- cmdstan_model(stan_file)
cmdstan_dir <- cmdstanr::cmdstan_path()   # pass to workers so they find CmdStan

## ---- pre-generate ALL datasets (reproducible, independent of worker order) ----
set.seed(master_seed)
Y <- t(replicate(M, rbinom(K, n, p_true_scen1)))   # M x K
seeds <- master_seed + seq_len(M)

## ---- one fit -> (pi_hat row, divergence count) ----
fit_one <- function(m) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(stan_file)         # loads up-to-date exe, no recompile
  fit <- mod$sample(
    data = list(K = K, n = n, y = Y[m, ],
                mu_mean = mu_mean, mu_sd = mu_sd,
                tau_scale = tau_scale, p0 = p0),
    seed = seeds[m],
    chains = n_chains, parallel_chains = 1,          # chains sequential INSIDE a worker
    iter_warmup = warmup, iter_sampling = sampling,
    adapt_delta = adapt_delta, max_treedepth = 12,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE
  )
  list(pi   = colMeans(fit$draws("exceeds_p0", format = "draws_matrix")),
       ndiv = sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]))
}

## ---- run (parallel across datasets, or serial fallback) ----
t0 <- Sys.time()
if (PARALLEL && requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)) {
  library(furrr); library(future)
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)
  res <- future_map(seq_len(M), fit_one,
                    .options = furrr_options(seed = TRUE,
                                             globals = c("Y","seeds","K","n","p0",
                                                         "mu_mean","mu_sd","tau_scale",
                                                         "n_chains","warmup","sampling",
                                                         "adapt_delta","stan_file","cmdstan_dir"),
                                             packages = "cmdstanr"))
  plan(sequential)
} else {
  if (PARALLEL) message("furrr/future not installed -> running SERIAL. ",
                        "install.packages(c('future','furrr')) to parallelize.")
  res <- lapply(seq_len(M), function(m) {
    if (m %% 100 == 0) cat(sprintf("  fit %d / %d\n", m, M))
    fit_one(m)
  })
}
cat(sprintf("done in %.1f min (%s, workers=%d)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            if (PARALLEL) "parallel" else "serial",
            if (PARALLEL) n_workers else 1L))

## ---- assemble cache ----
pi_hat   <- t(vapply(res, `[[`, numeric(K), "pi")); colnames(pi_hat) <- baskets
n_diverg <- vapply(res, `[[`, numeric(1), "ndiv")
saveRDS(list(pi_hat = pi_hat, Y = Y, n_diverg = n_diverg,
             hyper = list(tau_scale = tau_scale, mu_mean = mu_mean, mu_sd = mu_sd),
             knobs = list(M = M, n_chains = n_chains, warmup = warmup,
                          sampling = sampling, adapt_delta = adapt_delta)),
        file = "pihat_scenario1_full.rds")

## ---- divergence sanity ----
total_draws <- n_chains * sampling
cat(sprintf("\ndivergences: fits with >=1 = %d/%d (%.1f%%); overall rate = %.4f%% of draws\n",
            sum(n_diverg > 0), M, 100 * mean(n_diverg > 0),
            100 * sum(n_diverg) / (M * total_draws)))

## ---- c search (free: threshold the cache) ----
c_grid     <- seq(0.70, 0.99, by = 0.0025)
alpha_by_c <- sapply(c_grid, function(cc) colMeans(pi_hat > cc))
rownames(alpha_by_c) <- baskets
mean_alpha <- colMeans(alpha_by_c)

cat(sprintf("\nmean_alpha monotone non-increasing: %s\n", all(diff(mean_alpha) <= 1e-9)))

## ---- robust c*: explicit bracketed linear interpolation (mean_alpha decreasing) ----
locate_c <- function(cc, a, tgt) {
  if (min(a) > tgt || max(a) < tgt)
    return(c(c_star = NA_real_, note = NA_real_))   # not bracketed -> widen grid
  i_hi <- max(which(a >= tgt))                        # largest c with FP still >= target
  if (i_hi == length(cc)) return(c(c_star = cc[i_hi], note = 1))  # at edge
  a1 <- a[i_hi]; a2 <- a[i_hi + 1]; c1 <- cc[i_hi]; c2 <- cc[i_hi + 1]
  c(c_star = c1 + (tgt - a1) * (c2 - c1) / (a2 - a1), note = 0)
}
loc      <- locate_c(c_grid, mean_alpha, target)
c_star   <- unname(loc["c_star"])
se_alpha <- sqrt(target * (1 - target) / M)

cat(sprintf("\nc* (mean per-basket FP = %.2f)  =  %.4f\n", target, c_star))
if (is.na(c_star)) {
  cat(sprintf("  -> 0.10 NOT bracketed by c in [%.2f, %.2f]; widen the grid.\n",
              min(c_grid), max(c_grid)))
} else {
  k_near <- which.min(abs(c_grid - c_star))
  cat("per-basket alpha at nearest grid c (should agree within ~2 SE):\n")
  print(round(alpha_by_c[, k_near], 4))
  cat(sprintf("  [per-basket MC SE ~ %.4f at M=%d]\n", se_alpha, M))
}

## ---- carry forward: save c* for the power runs (step 3) ----
saveRDS(list(c_star = c_star, target = target,
             hyper = list(tau_scale = tau_scale, mu_mean = mu_mean, mu_sd = mu_sd)),
        file = "c_star.rds")
cat("\nsaved c_star.rds for step 3 (power across scenarios).\n")

