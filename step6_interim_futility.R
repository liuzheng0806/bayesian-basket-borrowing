## ============================================================================
## Step 6: SINGLE INTERIM FUTILITY analysis.
##
## Design (one interim look):
##   - Enroll to n=10, take an INTERIM look. A basket is STOPPED FOR FUTILITY if
##     Pr(p_k > p0 | interim data) < c_fut. Stopped baskets get final No-Go.
##   - Surviving baskets enroll to n=20; the FINAL look decides Go/No-Go by
##     Pr(p_k > p0 | final data) > c_final.
##   - Stopped baskets are FROZEN at their n=10 interim data inside the final model
##     (they still inform the shared pool, but cannot Go). So K stays 4, the per-basket
##     n vector is mixed (10 for stopped, 20 for survivors), which the locked Stan
##     models already accept as data -- no new Stan file needed.
##
## Two design knobs:
##   - c_fut (futility threshold): a DESIGN CHOICE, fixed (main: 0.10). Sweepable.
##   - c_final (efficacy threshold): RE-CALIBRATED under the two-stage procedure so
##     per-basket type I error returns to 0.10. Binding futility removes some null
##     baskets -> conservatism -> c_final calibrates LOWER than the no-interim c*
##     (you "spend back" the budget). This is the point of two-stage calibration.
##
## Why only borrowing methods matter here: with BHM/EXNEX the interim decision for
## basket k uses Pr(p_k>p0 | ALL baskets' interim data) through the shared (mu,tau),
## so interim stops are CORRELATED across baskets and the final pool depends on who
## survived -- the non-independence the README flags. no_borrow is included as the
## INDEPENDENT baseline that makes that contrast visible.
##
## Reuses locked bhm.stan / exnex.stan / no_borrow.stan and c_star_all.rds (priors +
## no-interim c* for reference). Writes interim_futility_results.rds.
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- run controls --------------------------------------------------------
DRY_RUN <- FALSE                  # TRUE -> M=50, short chains: local path check.
CFUT    <- 0.10                  # interim futility threshold (design choice).
                                 # To sweep, wrap the whole run in a loop over CFUT_GRID.
METHODS <- c("bhm", "exnex", "no_borrow")   # borrowing pair + independent baseline

M           <- if (DRY_RUN) 50L  else 2000L
n_chains    <- if (DRY_RUN) 2L   else 4L
warmup      <- if (DRY_RUN) 300L else 1000L
sampling    <- if (DRY_RUN) 300L else 1000L
master_seed <- 2026L
target      <- 0.10
interim_n   <- 10L
final_n     <- 20L

PARALLEL  <- TRUE
n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "55"))  # match allocation
cmdstan_dir <- cmdstanr::cmdstan_path()

## ---- locked design + priors ----------------------------------------------
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets); n_full <- rep(final_n, K); p0 <- 0.10; p1 <- 0.30
cs <- readRDS("c_star_all.rds")
mu_mean <- cs$hyper_bhm$mu_mean; mu_sd <- cs$hyper_bhm$mu_sd
tau_scale <- cs$hyper_bhm$tau_scale
a0 <- cs$a0; b0 <- cs$b0
m_nex <- cs$hyper_exnex$m_nex; s_nex <- cs$hyper_exnex$s_nex; pi_val <- cs$hyper_exnex$pi_val
cstar_no_interim <- list(bhm = cs$bhm, exnex = cs$exnex, no_borrow = cs$no_borrow)

scenarios <- list(s1 = c(0.10,0.10,0.10,0.10),   # null: calibration + type I error
                  s2 = c(0.30,0.30,0.30,0.30),
                  s3 = c(0.30,0.30,0.10,0.10),
                  s4 = c(0.40,0.30,0.15,0.05))
role_of <- function(pt) ifelse(pt <= p0+1e-9, "null", ifelse(pt >= p1-1e-9, "active", "gray"))

## ---- compile locked models ----------------------------------------------
spec_of <- list(
  bhm       = list(file = normalizePath("bhm.stan"),       type = "bhm",   shared = FALSE, adapt_delta = 0.99),
  exnex     = list(file = normalizePath("exnex.stan"),     type = "exnex", shared = FALSE, adapt_delta = 0.99),
  no_borrow = list(file = normalizePath("no_borrow.stan"), type = "beta",  shared = FALSE, adapt_delta = 0.95)
)
for (mn in METHODS) cmdstan_model(spec_of[[mn]]$file)

## ---- data list for a given (n vector, y vector, method) ------------------
make_data <- function(nv, yv, type) {
  if (type == "bhm")   return(list(K=K, n=nv, y=yv, mu_mean=mu_mean, mu_sd=mu_sd,
                                   tau_scale=tau_scale, p0=p0))
  if (type == "exnex") return(list(K=K, n=nv, y=yv, p0=p0, mu_mean=mu_mean, mu_sd=mu_sd,
                                   tau_scale=tau_scale, m_nex=m_nex, s_nex=s_nex,
                                   pi_ex=rep(pi_val, K)))
  list(K=K, n=nv, y=yv, a0=a0, b0=b0, p0=p0)                       # beta (no_borrow)
}

## ---- one model fit -> per-basket Pr(p_k>p0|data) -------------------------
fit_pi <- function(nv, yv, spec, seed) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(spec$file)
  fit <- mod$sample(data = make_data(nv, yv, spec$type), seed = seed,
                    chains = n_chains, parallel_chains = 1,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = spec$adapt_delta, max_treedepth = 12,
                    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  ex <- fit$draws("exceeds_p0", format = "draws_matrix")
  if (spec$shared) rep(mean(ex[,1]), K) else colMeans(ex)
}

## ---- TWO-STAGE procedure for one simulated trial -------------------------
## interim fit (all baskets, n=10) -> stops; final fit (mixed n) -> survivor Go probs.
fit_trial <- function(y_int, y_add, spec, c_fut, seed_int, seed_fin) {
  pi_int  <- fit_pi(rep(interim_n, K), y_int, spec, seed_int)      # Pr(p>p0 | interim)
  stopped <- pi_int < c_fut                                       # futility stop per basket
  n_fin <- ifelse(stopped, interim_n, final_n)                    # frozen at 10, else 20
  y_fin <- ifelse(stopped, y_int, y_int + y_add)
  pi_fin <- fit_pi(n_fin, y_fin, spec, seed_fin)                  # Pr(p>p0 | final), all K
  ## survivors decide on pi_fin; stopped baskets are forced No-Go (Go prob set to -1)
  go_prob <- ifelse(stopped, -1, pi_fin)
  list(pi_int = pi_int, stopped = stopped, go_prob = go_prob)
}

run_scenario <- function(pt, spec, c_fut, scen_idx) {
  set.seed(master_seed + 100 * scen_idx)
  Y_int <- t(replicate(M, rbinom(K, interim_n,          pt)))     # interim half
  Y_add <- t(replicate(M, rbinom(K, final_n - interim_n, pt)))    # second half
  sb_i <- master_seed + 1000 * scen_idx
  sb_f <- master_seed + 1000 * scen_idx + 500000L
  one <- function(m) fit_trial(Y_int[m,], Y_add[m,], spec, c_fut, sb_i + m, sb_f + m)
  if (PARALLEL && requireNamespace("furrr", quietly=TRUE) && requireNamespace("future", quietly=TRUE)) {
    library(furrr); library(future)
    plan(multisession, workers = n_workers); on.exit(plan(sequential), add = TRUE)
    res <- furrr::future_map(seq_len(M), one,
            .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
    plan(sequential)
  } else res <- lapply(seq_len(M), one)
  list(GO  = do.call(rbind, lapply(res, `[[`, "go_prob")),        # M x K (-1 = stopped)
       STP = do.call(rbind, lapply(res, `[[`, "stopped")),        # M x K logical
       PINT= do.call(rbind, lapply(res, `[[`, "pi_int")))
}

## ---- per-basket type I error at threshold c (stopped -> never Go) --------
alpha_at <- function(GO, c) colMeans(GO > c)        # GO=-1 for stopped never exceeds c>0

calibrate_cfinal <- function(GO_null, c_grid = seq(0.30, 0.99, by = 0.0025)) {
  mean_alpha <- sapply(c_grid, function(cc) mean(alpha_at(GO_null, cc)))
  if (min(mean_alpha) > target || max(mean_alpha) < target)
    return(list(c_final = NA_real_, bracketed = FALSE, ceiling = max(mean_alpha)))
  i_hi <- max(which(mean_alpha >= target))
  if (i_hi == length(c_grid)) return(list(c_final = c_grid[i_hi], bracketed = TRUE, boundary = TRUE))
  a1 <- mean_alpha[i_hi]; a2 <- mean_alpha[i_hi+1]
  list(c_final = c_grid[i_hi] + (target-a1)*(c_grid[i_hi+1]-c_grid[i_hi])/(a2-a1),
       bracketed = TRUE, boundary = FALSE)
}

## ==========================================================================
## RUN: per method -> calibrate c_final on s1, then OC on s2/s3/s4 (+ reuse s1)
## ==========================================================================
cat(sprintf("MODE: %s | M=%d | c_fut=%.2f | methods: %s | workers=%d\n",
            if (DRY_RUN) "DRY RUN" else "FULL", M, CFUT, paste(METHODS, collapse=", "), n_workers))

rows <- list(); cfinal <- list()
for (mn in METHODS) {
  spec <- spec_of[[mn]]
  cat(sprintf("\n========== [%s] ==========\n", mn))

  ## --- calibration on scenario 1 (two-stage), cache its OC too ---
  t0 <- Sys.time()
  s1 <- run_scenario(scenarios$s1, spec, CFUT, scen_idx = match("s1", names(scenarios)))
  cal <- calibrate_cfinal(s1$GO)
  cat(sprintf("  calib %.1f min | c_final=%s%s | (no-interim c*=%.3f)\n",
              as.numeric(difftime(Sys.time(), t0, units="mins")),
              ifelse(is.na(cal$c_final), "NA", sprintf("%.4f", cal$c_final)),
              if (!isTRUE(cal$bracketed)) sprintf(" [0.10 NOT reachable; ceiling=%.3f -> loosen c_fut]", cal$ceiling) else "",
              cstar_no_interim[[mn]]))
  if (is.na(cal$c_final)) { cat("  -> skipping sweep for this method.\n"); next }
  cfinal[[mn]] <- cal$c_final

  ## helper to turn a scenario run into OC rows
  emit <- function(si, pt, run) {
    Go  <- alpha_at(run$GO, cal$c_final)        # final Go freq (power / type I / gray)
    stp <- colMeans(run$STP)                    # interim futility-stop freq
    EN  <- interim_n + (final_n - interim_n) * (1 - stp)   # expected N per basket
    for (k in seq_len(K))
      rows[[length(rows)+1]] <<- data.frame(
        method = mn, scenario = si, basket = baskets[k], p_true = pt[k],
        role = role_of(pt[k]), c_final = cal$c_final,
        Go = Go[k], stop_interim = stp[k], EN = EN[k])
  }
  emit("s1", scenarios$s1, s1)                  # reuse calibration run

  ## --- sweep the remaining scenarios ---
  for (si in c("s2","s3","s4")) {
    t1 <- Sys.time()
    run <- run_scenario(scenarios[[si]], spec, CFUT, scen_idx = match(si, names(scenarios)))
    emit(si, scenarios[[si]], run)
    cat(sprintf("  %s %.1f min\n", si, as.numeric(difftime(Sys.time(), t1, units="mins"))))
  }
}
OC6 <- do.call(rbind, rows)
saveRDS(list(OC6 = OC6, c_final = cfinal, c_fut = CFUT, M = M,
             cstar_no_interim = cstar_no_interim, dry_run = DRY_RUN),
        file = "interim_futility_results.rds")
cat("\nsaved interim_futility_results.rds\n")

## ==========================================================================
## SUMMARY
## ==========================================================================
as_mat <- function(sub, col) {
  m <- matrix(NA_real_, length(METHODS), K, dimnames = list(METHODS, baskets))
  for (i in seq_len(nrow(sub))) m[sub$method[i], sub$basket[i]] <- sub[[col]][i]
  m
}
for (si in names(scenarios)) {
  sub <- OC6[OC6$scenario == si, ]
  if (!nrow(sub)) next
  cat(sprintf("\n================ %s  (true ORR: %s) ================\n", si,
              paste(sprintf("%.2f", scenarios[[si]]), collapse="/")))
  cat("final Go freq (active->power | null->type I error | gray->Go-rate):\n"); print(round(as_mat(sub,"Go"),3))
  cat("\ninterim futility-stop freq:\n"); print(round(as_mat(sub,"stop_interim"),3))
  cat("\nexpected N per basket (max 20):\n"); print(round(as_mat(sub,"EN"),1))
}

cat("\n================ HEADLINE ================\n")
cat("Recalibrated final thresholds vs no-interim c*:\n")
for (mn in names(cfinal))
  cat(sprintf("  %-9s c_final=%.4f  (no-interim %.4f, delta %+.4f)\n",
              mn, cfinal[[mn]], cstar_no_interim[[mn]], cfinal[[mn]] - cstar_no_interim[[mn]]))

cat("\nExpected total N (sum over baskets, max 80) by scenario:\n")
for (si in names(scenarios)) {
  sub <- OC6[OC6$scenario == si, ]; if (!nrow(sub)) next
  tot <- tapply(sub$EN, sub$method, sum)[names(cfinal)]
  cat(sprintf("  %s: %s\n", si, paste(sprintf("%s=%.1f", names(tot), tot), collapse="  ")))
}

## non-independence demo: CRC truth = 0.10 in BOTH s1 and s3; its interim-stop rate
## should be ~equal under no_borrow but DIFFER under BHM/EXNEX (active neighbours in s3
## lift the shared pool -> CRC looks less futile at interim).
cat("\nNON-INDEPENDENCE (CRC stop rate, same truth 0.10 in s1 & s3):\n")
for (mn in names(cfinal)) {
  a <- OC6$stop_interim[OC6$method==mn & OC6$scenario=="s1" & OC6$basket=="CRC"]
  b <- OC6$stop_interim[OC6$method==mn & OC6$scenario=="s3" & OC6$basket=="CRC"]
  cat(sprintf("  %-9s s1=%.3f  s3=%.3f  (diff %+.3f)%s\n", mn, a, b, b-a,
              if (mn=="no_borrow") "  <- baseline: should be ~0" else "  <- borrowing -> shifts"))
}
## the core trade-off: futility costs power but saves patients (vs single-stage)
cat("\nFUTILITY TRADE-OFF vs single-stage (active baskets only; baseline = oc_results_4m.rds):\n")
base_OC <- tryCatch(readRDS("oc_results_4m.rds")$OC, error = function(e) NULL)
if (!is.null(base_OC)) {
  for (si in c("s2","s3","s4")) {
    cat(sprintf("  %s:\n", si))
    for (mn in names(cfinal)) {
      act <- OC6$role=="active" & OC6$scenario==si & OC6$method==mn
      pw2 <- mean(OC6$Go[act])                                   # two-stage power
      bb  <- base_OC$role=="active" & base_OC$scenario==si & base_OC$method==mn
      pw1 <- if (any(bb)) mean(base_OC$Go[bb]) else NA_real_     # single-stage power
      en  <- sum(OC6$EN[OC6$scenario==si & OC6$method==mn])      # expected total N (max 80)
      cat(sprintf("    %-9s power %0.3f vs %0.3f (cost %+.3f) | expected total N %.1f / 80\n",
                  mn, pw2, pw1, pw2 - pw1, en))
    }
  }
  cat("  -> read as: how much power futility gives up, against how many patients it saves.\n")
} else cat("  (oc_results_4m.rds not found; skipping single-stage comparison.)\n")

cat(if (DRY_RUN) "\n*** DRY RUN: M=50, numbers are noise; checking the pipeline only. ***\n" else "")
