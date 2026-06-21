## ============================================================================
## Step 3b (4-METHOD): operating-characteristics sweep, now including EXNEX.
##
## Drop-in successor to step3b_oc_sweep.R. Scenarios 2/3/4 x methods
##   {no-borrow, BHM, complete-pool, EXNEX},  each using its OWN calibrated c*.
## Within a scenario all FOUR methods are fit on the SAME simulated datasets
## (paired). The data-generation seeds are UNCHANGED from the 3-method script, so
## the no-borrow / BHM / pool columns reproduce the previously reported numbers
## exactly, and EXNEX is evaluated on identical data -- this is why we re-run all
## four together rather than bolting EXNEX on separately (avoids data mismatch).
##
## Outputs per (scenario, method, basket): Go frequency (= power if active,
## type I error if null, Go-rate if gray), bias / MSE / mean posterior SD, and --
## for EXNEX only -- the mean posterior EX weight (the borrow-vs-self diagnostic).
##
## Needs: no_borrow.stan, bhm.stan, pool.stan, exnex.stan, c_star_all.rds  (with
##        $exnex present -- run step3a_exnex_calibrate.R first).
## Writes: oc_results_4m.rds   (kept separate from the 3-method oc_results.rds).
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- locked design ----
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets)
n  <- rep(20L, K)
p0 <- 0.10
p1 <- 0.30                       # clinically meaningful rate (defines the gray zone)

## ---- load calibrated thresholds + the priors they were calibrated under ----
cs <- readRDS("c_star_all.rds")
stopifnot(!is.null(cs$exnex))    # run step3a_exnex_calibrate.R first
c_by_method <- list(no_borrow = cs$no_borrow, bhm = cs$bhm,
                    pool = cs$pool, exnex = cs$exnex)
a0 <- cs$a0; b0 <- cs$b0                                   # Beta prior (brackets)
tau_scale <- cs$hyper_bhm$tau_scale                        # BHM / EX hyperprior
mu_mean   <- cs$hyper_bhm$mu_mean
mu_sd     <- cs$hyper_bhm$mu_sd
m_nex  <- cs$hyper_exnex$m_nex                             # NEX prior
s_nex  <- cs$hyper_exnex$s_nex
pi_ex  <- rep(cs$hyper_exnex$pi_val, K)                    # fixed EX weight (main analysis)
cat(sprintf("thresholds:  no-borrow %.4f | BHM %.4f | pool %.4f | EXNEX %.4f\n",
            cs$no_borrow, cs$bhm, cs$pool, cs$exnex))
cat(sprintf("EXNEX: pi=%.2f, NEX N(%.1f, %.1f^2)\n", pi_ex[1], m_nex, s_nex))

## ---- knobs (UNCHANGED from the 3-method script: reproduces the first 3 cols) ----
SCEN_M      <- 2000L
n_chains    <- 4L
warmup      <- 1000L
sampling    <- 1000L
master_seed <- 2026L
PARALLEL    <- TRUE
n_workers   <- max(1L, parallel::detectCores() - 1L)
cmdstan_dir <- cmdstanr::cmdstan_path()

## ---- compile all four models once ----
nb_file    <- normalizePath("no_borrow.stan"); cmdstan_model(nb_file)
bhm_file   <- normalizePath("bhm.stan");       cmdstan_model(bhm_file)
pool_file  <- normalizePath("pool.stan");      cmdstan_model(pool_file)
exnex_file <- normalizePath("exnex.stan");     cmdstan_model(exnex_file)

methods <- list(
  no_borrow = list(file = nb_file,    shared = FALSE, type = "beta",  adapt_delta = 0.95),
  bhm       = list(file = bhm_file,   shared = FALSE, type = "bhm",   adapt_delta = 0.99),
  pool      = list(file = pool_file,  shared = TRUE,  type = "beta",  adapt_delta = 0.95),
  exnex     = list(file = exnex_file, shared = FALSE, type = "exnex", adapt_delta = 0.99)
)

scenarios <- list(
  s2 = c(0.30, 0.30, 0.30, 0.30),
  s3 = c(0.30, 0.30, 0.10, 0.10),
  s4 = c(0.40, 0.30, 0.15, 0.05)
)
scen_label <- c(s2 = "Scenario 2 (all 0.30, homogeneous active)",
                s3 = "Scenario 3 (0.30/0.30/0.10/0.10, mixed)",
                s4 = "Scenario 4 (0.40/0.30/0.15/0.05, heterogeneous)")

role_of <- function(pt) ifelse(pt <= p0 + 1e-9, "null",
                        ifelse(pt >= p1 - 1e-9, "active", "gray"))

## ---- one fit -> per-basket (pi_hat, post mean, post sd, EX weight) + divergences ----
fit_one <- function(m, Ymat, spec, seed_base) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(spec$file)
  y_m <- Ymat[m, ]
  dat <- if (spec$type == "bhm")
    list(K = K, n = n, y = y_m, mu_mean = mu_mean, mu_sd = mu_sd,
         tau_scale = tau_scale, p0 = p0)
  else if (spec$type == "exnex")
    list(K = K, n = n, y = y_m, p0 = p0, mu_mean = mu_mean, mu_sd = mu_sd,
         tau_scale = tau_scale, m_nex = m_nex, s_nex = s_nex, pi_ex = pi_ex)
  else
    list(K = K, n = n, y = y_m, a0 = a0, b0 = b0, p0 = p0)
  fit <- mod$sample(
    data = dat, seed = seed_base + m, chains = n_chains, parallel_chains = 1,
    iter_warmup = warmup, iter_sampling = sampling, adapt_delta = spec$adapt_delta,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  pe <- fit$draws("p",          format = "draws_matrix")
  ex <- fit$draws("exceeds_p0", format = "draws_matrix")
  if (spec$shared) {                     # complete pooling: one value -> all K baskets
    pii <- rep(mean(ex[, 1]), K); pmn <- rep(mean(pe[, 1]), K); psd <- rep(sd(pe[, 1]), K)
  } else {
    pii <- colMeans(ex); pmn <- colMeans(pe); psd <- apply(pe, 2, sd)
  }
  wex <- if (spec$type == "exnex")       # EXNEX-only diagnostic; NA elsewhere
    colMeans(fit$draws("w_ex", format = "draws_matrix")) else rep(NA_real_, K)
  list(pi = pii, pmn = pmn, psd = psd, wex = wex,
       ndiv = sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]))
}

run_combo <- function(Ymat, spec, seed_base) {
  if (PARALLEL && requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    library(furrr); library(future)
    plan(multisession, workers = n_workers); on.exit(plan(sequential), add = TRUE)
    res <- furrr::future_map(seq_len(SCEN_M), ~ fit_one(.x, Ymat, spec, seed_base),
            .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
    plan(sequential); res
  } else {
    if (PARALLEL) message("furrr/future not found -> SERIAL.")
    lapply(seq_len(SCEN_M), function(m) fit_one(m, Ymat, spec, seed_base))
  }
}

## ---- sweep ----
rows <- list(); raw_PI <- list()
for (si in names(scenarios)) {
  pt <- scenarios[[si]]
  set.seed(master_seed + 100 * match(si, names(scenarios)))   # UNCHANGED seed
  Y_s       <- t(replicate(SCEN_M, rbinom(K, n, pt)))         # same data for all 4 methods
  seed_base <- master_seed + 1000 * match(si, names(scenarios))
  cat(sprintf("\n##### %s #####\n", scen_label[si]))
  for (mn in names(methods)) {
    t0  <- Sys.time()
    res <- run_combo(Y_s, methods[[mn]], seed_base)
    PI  <- do.call(rbind, lapply(res, `[[`, "pi"))
    PMN <- do.call(rbind, lapply(res, `[[`, "pmn"))
    PSD <- do.call(rbind, lapply(res, `[[`, "psd"))
    WEX <- do.call(rbind, lapply(res, `[[`, "wex"))
    ndv <- vapply(res, `[[`, numeric(1), "ndiv")
    cm  <- c_by_method[[mn]]
    Go  <- colMeans(PI > cm)
    bias <- colMeans(sweep(PMN, 2, pt)); mse <- colMeans(sweep(PMN, 2, pt)^2)
    psdm <- colMeans(PSD); wexm <- colMeans(WEX)
    cat(sprintf("  [%s] %.1f min, divergent fits %.1f%%\n", mn,
                as.numeric(difftime(Sys.time(), t0, units = "mins")), 100 * mean(ndv > 0)))
    for (k in seq_len(K))
      rows[[length(rows) + 1]] <- data.frame(
        scenario = si, method = mn, basket = baskets[k], p_true = pt[k],
        role = role_of(pt[k]), Go = Go[k], Go_se = sqrt(Go[k] * (1 - Go[k]) / SCEN_M),
        bias = bias[k], MSE = mse[k], post_sd = psdm[k], ex_weight = wexm[k])
    raw_PI[[paste(si, mn, sep = "_")]] <- PI
  }
}
OC <- do.call(rbind, rows)
saveRDS(list(OC = OC, raw_PI = raw_PI, M = SCEN_M,
             c_by_method = c_by_method, scenarios = scenarios),
        file = "oc_results_4m.rds")
cat("\nsaved oc_results_4m.rds (4-method table + raw posterior probs for figures).\n")

## ---- pretty per-scenario tables ----
as_mat <- function(sub, col) {
  m <- matrix(NA_real_, length(methods), K, dimnames = list(names(methods), baskets))
  for (i in seq_len(nrow(sub))) m[sub$method[i], sub$basket[i]] <- sub[[col]][i]
  m
}
for (si in names(scenarios)) {
  sub <- OC[OC$scenario == si, ]
  cat(sprintf("\n================ %s ================\n", scen_label[si]))
  cat("true ORR: ", paste(sprintf("%s=%.2f", baskets, scenarios[[si]]), collapse = "  "), "\n")
  cat("role:     ", paste(sprintf("%s=%s", baskets, sapply(scenarios[[si]], role_of)),
                          collapse = "  "), "\n")
  cat("\nGo frequency  (active->power | null->type I error | gray->Go-rate):\n")
  print(round(as_mat(sub, "Go"), 3))
  cat("\nbias (post mean - truth):\n"); print(round(as_mat(sub, "bias"), 3))
  cat("\nMSE:\n");                       print(round(as_mat(sub, "MSE"), 4))
  cat("\nmean posterior SD:\n");         print(round(as_mat(sub, "post_sd"), 3))
  cat("\nEXNEX posterior EX weight (borrow<->self; NA for non-EXNEX):\n")
  print(round(as_mat(sub, "ex_weight")["exnex", , drop = FALSE], 3))
}

## ---- headline numbers for the interview story ----
cat("\n================ HEADLINE ================\n")
s2 <- OC[OC$scenario == "s2", ]
pw <- tapply(s2$Go, s2$method, mean)[names(methods)]
cat("Scenario 2 power (avg over baskets, all truly active):\n"); print(round(pw, 3))
cat("  -> borrowing benefit = each borrowing method's power minus no-borrow power.\n")
crc <- OC[OC$scenario == "s4" & OC$basket == "CRC", ]
ti  <- setNames(crc$Go, crc$method)[names(methods)]
cat("\nScenario 4 CRC type I error (truly null, p=0.05; target was 0.10):\n"); print(round(ti, 3))
cat("  -> KEY EXNEX TEST: does EXNEX pull CRC's false positive below BHM's?\n")
exw <- setNames(crc$ex_weight, crc$method)["exnex"]
cat(sprintf("\nScenario 4 CRC posterior EX weight (EXNEX) = %.3f\n", exw))
cat("  -> the lower this is, the more CRC has defected to NEX (self-estimates,\n",
    "     escaping the upward drag) -- the mechanism behind any FP reduction.\n", sep = "")
