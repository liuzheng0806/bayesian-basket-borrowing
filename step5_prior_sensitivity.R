## ============================================================================
## Step 5: PRIOR SENSITIVITY (config-table driven).
##
## Each prior change moves the null distribution of Pr(p_k>p0|y), so every config
## must RECALIBRATE its own c* (scenario 1) and THEN run the OC sweep with that c*.
## This script automates "calibrate -> sweep" per config.
##
## Scope note: no-borrow and complete-pool have NO hierarchical/mixture prior, so
## they are prior-INVARIANT -- their c* and OC do not move. They are therefore NOT
## in the config table; their locked values live in oc_results_4m.rds for reference.
## Only the borrowing methods (BHM, EXNEX) are swept here.
##
## Pairing: calibration reuses the cached scenario-1 null datasets
## (pihat_scenario1_full.rds$Y); the sweep regenerates each scenario's data with the
## SAME seeds as step3b. So at full M each config differs from the main analysis ONLY
## in the prior, on identical datasets -> the deltas are clean.
##
## Needs: bhm_psens.stan, exnex_psens.stan, exnex_pi_random.stan,
##        pihat_scenario1_full.rds, c_star_all.rds, oc_results_4m.rds.
## Writes: prior_sensitivity_results.rds
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- run controls --------------------------------------------------------
DRY_RUN  <- TRUE        # TRUE -> M=50, 2 chains: validate the code path locally, fast.
                        # set FALSE for the cluster run (M=2000, paired with main analysis).
RUN_ONLY <- NULL        # NULL = all configs; or a character vector of labels to subset,
                        # e.g. c("bhm_tau_hcauchy","exnex_tau_hcauchy","exnex_pi_beta11",
                        #        "exnex_nex_logit_p0") for the high-value first pass.

M_calib <- if (DRY_RUN)  50L else 2000L
M_sweep <- if (DRY_RUN)  50L else 2000L
n_chains    <- if (DRY_RUN) 2L   else 4L
warmup      <- if (DRY_RUN) 300L else 1000L
sampling    <- if (DRY_RUN) 300L else 1000L
adapt_delta <- 0.99     # all swept methods borrow -> funnel-grade
max_td      <- 12L
master_seed <- 2026L    # MUST match step3b for paired sweep data
target      <- 0.10

PARALLEL  <- TRUE
n_workers <- max(1L, parallel::detectCores() - 1L)
cmdstan_dir <- cmdstanr::cmdstan_path()

## ---- locked design (same as step3b) -------------------------------------
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets); n <- rep(20L, K); p0 <- 0.10; p1 <- 0.30
scenarios <- list(s2 = c(0.30,0.30,0.30,0.30),
                  s3 = c(0.30,0.30,0.10,0.10),
                  s4 = c(0.40,0.30,0.15,0.05))
role_of <- function(pt) ifelse(pt <= p0+1e-9, "null", ifelse(pt >= p1-1e-9, "active", "gray"))

## ---- fixed EX hyperprior (mu held at main-analysis values; only tau prior varies) ----
cs <- readRDS("c_star_all.rds")
mu_mean <- cs$hyper_bhm$mu_mean
mu_sd   <- cs$hyper_bhm$mu_sd

## ---- cached scenario-1 null data (paired with main calibration) ----------
Y1_full <- readRDS("pihat_scenario1_full.rds")$Y
stopifnot(M_calib <= nrow(Y1_full))
Y1 <- Y1_full[seq_len(M_calib), , drop = FALSE]

## ---- compile the three psens models once ---------------------------------
bhm_ps   <- normalizePath("bhm_psens.stan");      cmdstan_model(bhm_ps)
exnex_ps <- normalizePath("exnex_psens.stan");    cmdstan_model(exnex_ps)
exnex_pr <- normalizePath("exnex_pi_random.stan");cmdstan_model(exnex_pr)

## ==========================================================================
## CONFIG TABLE  (add a row = add a sensitivity config; one line each)
##   type    : "bhm" | "exnex" | "exnex_pi_random"   (selects data layout + model)
##   ref     : "bhm" | "exnex"  (which main-analysis column to delta against)
##   tau_family: 1 HN | 2 HC | 3 IG(on tau^2);  tau_scale used by HN/HC;  ig_a/ig_b by IG
##   pi_val  : fixed EX weight (exnex);  nex_center: NEX mean on logit scale (0 or logit(p0))
##   a_pi/b_pi: Beta hyperprior (exnex_pi_random only)
## Main analysis (HN s=1, pi=0.5, NEX center 0) is the REFERENCE, already in
## oc_results_4m.rds -- not re-run here.
## ==========================================================================
mk <- function(label, type, ref, file, tau_family=1L, tau_scale=1.0,
               ig_a=0.001, ig_b=0.001, pi_val=0.5, nex_center=0.0, s_nex=2.0,
               a_pi=1.0, b_pi=1.0, adapt_delta=0.99, max_tree=12L)
  list(label=label, type=type, ref=ref, file=file, tau_family=as.integer(tau_family),
       tau_scale=tau_scale, ig_a=ig_a, ig_b=ig_b, pi_val=pi_val,
       nex_center=nex_center, s_nex=s_nex, a_pi=a_pi, b_pi=b_pi,
       adapt_delta=adapt_delta, max_tree=max_tree)

## NOTE on sampler settings (set by model geometry, not after-the-fact patching):
##   BHM is NON-CENTERED (theta=mu+tau*z) -> adapt_delta 0.99 suffices (clean).
##   EXNEX is a CENTERED mixture (log_mix; non-centering does not apply) -> funnel-
##   prone BY DESIGN, so all EXNEX configs use adapt_delta=0.999, max_tree=15 as the
##   standard for that model class. (IG-EXNEX still diverges heavily even at 0.999 ->
##   that is PRIOR pathology, not a sampler issue: the cleanest reverse-teaching point.)
configs <- list(
  ## --- tau scale (borrowing strength) : BHM + EXNEX EX-component ---
  mk("bhm_tau_hn_s0.5",   "bhm",   "bhm",   bhm_ps,   tau_family=1, tau_scale=0.5),
  mk("bhm_tau_hn_s2",     "bhm",   "bhm",   bhm_ps,   tau_family=1, tau_scale=2.0),
  mk("exnex_tau_hn_s0.5", "exnex", "exnex", exnex_ps, tau_family=1, tau_scale=0.5, adapt_delta=0.999, max_tree=15),
  mk("exnex_tau_hn_s2",   "exnex", "exnex", exnex_ps, tau_family=1, tau_scale=2.0, adapt_delta=0.999, max_tree=15),
  ## --- tau prior FAMILY (the selling point) ---
  mk("bhm_tau_hcauchy",   "bhm",   "bhm",   bhm_ps,   tau_family=2, tau_scale=1.0),
  mk("bhm_tau_invgamma",  "bhm",   "bhm",   bhm_ps,   tau_family=3, ig_a=0.001, ig_b=0.001),
  mk("exnex_tau_hcauchy", "exnex", "exnex", exnex_ps, tau_family=2, tau_scale=1.0, adapt_delta=0.999, max_tree=15),
  mk("exnex_tau_invgamma","exnex", "exnex", exnex_ps, tau_family=3, ig_a=0.001, ig_b=0.001, adapt_delta=0.999, max_tree=15),
  ## --- EXNEX mixing weight pi ---
  mk("exnex_pi_0.3",      "exnex", "exnex", exnex_ps, pi_val=0.3, adapt_delta=0.999, max_tree=15),
  mk("exnex_pi_0.7",      "exnex", "exnex", exnex_ps, pi_val=0.7, adapt_delta=0.999, max_tree=15),
  mk("exnex_pi_beta11",   "exnex_pi_random", "exnex", exnex_pr, a_pi=1.0, b_pi=1.0, adapt_delta=0.999, max_tree=15),
  ## --- EXNEX NEX center (the diagnostic-driven axis) ---
  mk("exnex_nex_logit_p0","exnex", "exnex", exnex_ps, nex_center=logit(p0), adapt_delta=0.999, max_tree=15)
)
names(configs) <- vapply(configs, `[[`, "", "label")
if (!is.null(RUN_ONLY)) configs <- configs[intersect(RUN_ONLY, names(configs))]
cat(sprintf("configs to run (%d): %s\n", length(configs), paste(names(configs), collapse=", ")))
cat(sprintf("MODE: %s | M_calib=%d M_sweep=%d chains=%d\n",
            if (DRY_RUN) "DRY RUN (path check only)" else "FULL", M_calib, M_sweep, n_chains))

## ---- build the data list for a given config + dataset --------------------
make_data <- function(cfg, y) {
  base <- list(K=K, n=n, y=y, p0=p0, mu_mean=mu_mean, mu_sd=mu_sd,
               tau_prior_family=cfg$tau_family, tau_scale=cfg$tau_scale,
               ig_a=cfg$ig_a, ig_b=cfg$ig_b)
  if (cfg$type == "bhm") return(base)
  if (cfg$type == "exnex")
    return(c(base, list(m_nex=cfg$nex_center, s_nex=cfg$s_nex, pi_ex=rep(cfg$pi_val, K))))
  ## exnex_pi_random:
  c(base, list(m_nex=cfg$nex_center, s_nex=cfg$s_nex, a_pi=cfg$a_pi, b_pi=cfg$b_pi))
}

## ---- one fit -> per-basket summaries -------------------------------------
fit_one <- function(y, cfg, seed) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(cfg$file)
  fit <- mod$sample(data = make_data(cfg, y), seed = seed,
                    chains = n_chains, parallel_chains = 1,
                    iter_warmup = warmup, iter_sampling = sampling,
                    adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_tree,
                    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  pe <- fit$draws("p",          format = "draws_matrix")
  ex <- fit$draws("exceeds_p0", format = "draws_matrix")
  wex <- if (cfg$type != "bhm")
    colMeans(fit$draws("w_ex", format = "draws_matrix")) else rep(NA_real_, K)
  list(pi = colMeans(ex), pmn = colMeans(pe), psd = apply(pe, 2, sd), wex = wex,
       ndiv = sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]))
}

run_over <- function(Ymat, cfg, seed_base) {
  M <- nrow(Ymat)
  if (PARALLEL && requireNamespace("furrr", quietly=TRUE) && requireNamespace("future", quietly=TRUE)) {
    library(furrr); library(future)
    plan(multisession, workers = n_workers); on.exit(plan(sequential), add = TRUE)
    res <- furrr::future_map(seq_len(M), ~ fit_one(Ymat[.x, ], cfg, seed_base + .x),
            .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
    plan(sequential); res
  } else {
    if (PARALLEL) message("furrr/future not found -> SERIAL.")
    lapply(seq_len(M), function(m) fit_one(Ymat[m, ], cfg, seed_base + m))
  }
}

## ---- bracketed c* (same procedure as step3a) -----------------------------
calibrate_cstar <- function(pi_hat, c_grid = seq(0.40, 0.99, by = 0.0025)) {
  alpha_by_c <- sapply(c_grid, function(cc) colMeans(pi_hat > cc))
  mean_alpha <- colMeans(alpha_by_c)
  monotone <- all(diff(mean_alpha) <= 1e-9)
  if (min(mean_alpha) > target || max(mean_alpha) < target)
    return(list(c_star = NA_real_, boundary = NA, monotone = monotone))
  i_hi <- max(which(mean_alpha >= target))
  if (i_hi == length(c_grid)) return(list(c_star = c_grid[i_hi], boundary = TRUE, monotone = monotone))
  a1 <- mean_alpha[i_hi]; a2 <- mean_alpha[i_hi+1]
  list(c_star = c_grid[i_hi] + (target-a1)*(c_grid[i_hi+1]-c_grid[i_hi])/(a2-a1),
       boundary = FALSE, monotone = monotone)
}

## ---- per-config: calibrate then sweep ------------------------------------
run_config <- function(cfg) {
  cat(sprintf("\n========== [%s]  type=%s  tau_family=%d ==========\n",
              cfg$label, cfg$type, cfg$tau_family))
  ## 1) calibrate c* on cached scenario-1 nulls
  t0 <- Sys.time()
  cal_res <- run_over(Y1, cfg, seed_base = master_seed)
  pi_hat  <- do.call(rbind, lapply(cal_res, `[[`, "pi")); colnames(pi_hat) <- baskets
  cal     <- calibrate_cstar(pi_hat)
  ndiv_c  <- vapply(cal_res, `[[`, numeric(1), "ndiv")
  draws_per_fit <- n_chains * sampling
  div_rows <- list()
  add_div <- function(stage, ndv) data.frame(config=cfg$label, type=cfg$type, stage=stage,
               pct_fits_div = 100*mean(ndv > 0),
               draw_div_pct = 100*sum(ndv)/(length(ndv)*draws_per_fit))
  div_rows[["calib"]] <- add_div("calib", ndiv_c)
  cat(sprintf("  calib %.1f min | c*=%.4f%s | monotone=%s | divergent fits %.1f%% (draw rate %.3f%%)\n",
              as.numeric(difftime(Sys.time(), t0, units="mins")), cal$c_star,
              if (isTRUE(cal$boundary)) " (BOUNDARY-extend grid!)" else "",
              cal$monotone, 100*mean(ndiv_c > 0), 100*sum(ndiv_c)/(length(ndiv_c)*draws_per_fit)))
  if (is.na(cal$c_star)) { cat("  -> 0.10 not bracketed; skipping sweep.\n")
    return(structure(NULL, div = do.call(rbind, div_rows))) }

  ## 2) sweep scenarios 2/3/4 with this c* (paired data via step3b seeds)
  rows <- list()
  for (si in names(scenarios)) {
    pt <- scenarios[[si]]
    set.seed(master_seed + 100 * match(si, names(scenarios)))    # SAME as step3b
    Y_s <- t(replicate(M_sweep, rbinom(K, n, pt)))
    seed_base <- master_seed + 1000 * match(si, names(scenarios))
    t1 <- Sys.time()
    res <- run_over(Y_s, cfg, seed_base)
    PI  <- do.call(rbind, lapply(res, `[[`, "pi"))
    PMN <- do.call(rbind, lapply(res, `[[`, "pmn"))
    PSD <- do.call(rbind, lapply(res, `[[`, "psd"))
    WEX <- do.call(rbind, lapply(res, `[[`, "wex"))
    ndv <- vapply(res, `[[`, numeric(1), "ndiv")
    div_rows[[si]] <- add_div(si, ndv)
    Go  <- colMeans(PI > cal$c_star)
    bias <- colMeans(sweep(PMN, 2, pt)); mse <- colMeans(sweep(PMN, 2, pt)^2)
    cat(sprintf("  sweep %s %.1f min | divergent fits %.1f%% (draw rate %.3f%%)\n", si,
                as.numeric(difftime(Sys.time(), t1, units="mins")),
                100*mean(ndv > 0), 100*sum(ndv)/(length(ndv)*draws_per_fit)))
    for (k in seq_len(K))
      rows[[length(rows)+1]] <- data.frame(
        config=cfg$label, type=cfg$type, ref=cfg$ref, c_star=cal$c_star,
        scenario=si, basket=baskets[k], p_true=pt[k], role=role_of(pt[k]),
        Go=Go[k], bias=bias[k], MSE=mse[k], post_sd=colMeans(PSD)[k],
        ex_weight=colMeans(WEX)[k])
  }
  out <- do.call(rbind, rows)
  structure(out, div = do.call(rbind, div_rows))
}

## ---- run all configs -----------------------------------------------------
all_rows <- list(); all_div <- list()
for (lab in names(configs)) {
  r <- run_config(configs[[lab]])
  all_div[[lab]] <- attr(r, "div")
  if (!is.null(r)) all_rows[[lab]] <- r
}
PS  <- do.call(rbind, all_rows)
DIV <- do.call(rbind, all_div); rownames(DIV) <- NULL
saveRDS(list(PS = PS, DIV = DIV, configs = configs, M_calib = M_calib, M_sweep = M_sweep,
             dry_run = DRY_RUN),
        file = "prior_sensitivity_results.rds")
cat("\nsaved prior_sensitivity_results.rds (PS = OC table, DIV = per-config divergences)\n")

## ---- divergence table (also stored as $DIV) ------------------------------
cat("\n================ DIVERGENCES (per config x stage) ================\n")
cat("pct_fits_div = % of fits with >=1 divergence; draw_div_pct = % of ALL draws divergent\n")
DIVp <- DIV; DIVp$pct_fits_div <- round(DIVp$pct_fits_div,1); DIVp$draw_div_pct <- round(DIVp$draw_div_pct,3)
print(DIVp, row.names = FALSE)
divmax <- tapply(DIV$pct_fits_div, DIV$config, max)   # worst stage per config, for summary

## ==========================================================================
## SUMMARY: each config's headline vs the MAIN-ANALYSIS baseline
##   benefit  = scenario-2 power (avg over baskets)
##   cost     = scenario-4 CRC type I error (truly null, p=0.05)
##   (exnex)  = scenario-4 CRC posterior EX weight
## Baseline pulled from oc_results_4m.rds so nothing is hard-coded.
## ==========================================================================
base_OC <- readRDS("oc_results_4m.rds")$OC
base_cstar <- list(bhm = cs$bhm, exnex = cs$exnex)
base_metrics <- function(meth) {
  s2 <- base_OC[base_OC$scenario=="s2" & base_OC$method==meth, ]
  crc <- base_OC[base_OC$scenario=="s4" & base_OC$basket=="CRC" & base_OC$method==meth, ]
  list(pow = mean(s2$Go), crc_ti = crc$Go, crc_exw = crc$ex_weight)
}
cfg_metrics <- function(lab) {
  sub <- PS[PS$config==lab, ]
  s2 <- sub[sub$scenario=="s2", ]
  crc <- sub[sub$scenario=="s4" & sub$basket=="CRC", ]
  list(cstar = sub$c_star[1], ref = sub$ref[1],
       pow = mean(s2$Go), crc_ti = crc$Go, crc_exw = crc$ex_weight)
}

cat("\n================ PRIOR-SENSITIVITY SUMMARY (delta vs main analysis) ================\n")
if (DRY_RUN) cat("*** DRY RUN: M=50, NOT paired with the M=2000 baseline -- deltas are noise. ***\n")
cat(sprintf("%-22s %6s %8s | %6s %7s | %7s %7s | %7s | %7s\n",
            "config","c*","dc*","s2pow","dPow","s4CRCti","dTI","CRCexw","maxdiv%"))
for (lab in names(all_rows)) {
  m <- cfg_metrics(lab); b <- base_metrics(m$ref); bc <- base_cstar[[m$ref]]
  cat(sprintf("%-22s %6.3f %+8.3f | %6.3f %+7.3f | %7.3f %+7.3f | %7s | %6.1f%%\n",
              lab, m$cstar, m$cstar - bc, m$pow, m$pow - b$pow,
              m$crc_ti, m$crc_ti - b$crc_ti,
              if (is.na(m$crc_exw)) "  --" else sprintf("%.3f", m$crc_exw),
              if (lab %in% names(divmax)) as.numeric(divmax[[lab]]) else NA_real_))
}
cat("\nReading: dTI > 0 means the prior INFLATES CRC false positive vs the locked prior;\n")
cat("dPow is the power change. A robust prior keeps dTI small without crushing dPow.\n")
cat(sprintf("Baseline (main analysis): BHM c*=%.3f pow=%.3f CRCti=%.3f | EXNEX c*=%.3f pow=%.3f CRCti=%.3f CRCexw=%.3f\n",
            cs$bhm, base_metrics("bhm")$pow, base_metrics("bhm")$crc_ti,
            cs$exnex, base_metrics("exnex")$pow, base_metrics("exnex")$crc_ti,
            base_metrics("exnex")$crc_exw))
