## ============================================================================
## Step 6b: PREDICTIVE-PROBABILITY (BOP2-style) interim futility, WITH borrowing
##          inside the prediction (option A), full nested MCMC (option c).
##
## WHY this exists: step6 used a POSTERIOR-probability futility rule (same family as
## the final Go rule: Pr(p_k>p0|interim) < c_fut). At n=10 the posterior is too wide
## for that to trigger, so futility was a near-no-op. BOP2 instead asks the forward
## question -- "if we enroll the rest, how likely is this basket to PASS at the end?"
## -- via a PREDICTIVE probability, which is far more sensitive at small n.
##
## ---------------------------------------------------------------------------
## DESIGN (every choice is a decision; the subtle ones are flagged):
##
## Predictive probability for basket k:
##   PPP_k = P( basket k passes the final Go gate | interim data, all baskets
##              complete to n=20 ),
##   estimated by nested Monte Carlo:
##     1. fit the JOINT borrowing model on interim data (n1=10)            [1 MCMC]
##     2. for r=1..R: draw p^(r) from the interim posterior, simulate the future
##        y2^(r) ~ Binom(n2, p^(r)) for ALL baskets, form the complete data
##        (y1+y2^(r), n=20), and FIT THE JOINT borrowing model again        [R MCMC]
##        -> one joint fit yields Pr(p_k>p0 | final^(r)) for every basket at once.
##     3. PPP_k = mean_r 1[ Pr(p_k>p0 | final^(r)) > c_ref ].
##   Futility: stop basket k at interim if PPP_k < theta_fut (= 0.05, fixed).
##
##   ** Borrowing IS inside the prediction (option A): the predictive fits are the
##      joint BHM/EXNEX model, so an inactive basket's PPP is lifted by active
##      neighbours through the shared (mu,tau). This is exactly what lets us MEASURE
##      whether borrowing makes futility fail in s4 (the project's main thread). **
##
## THE CIRCULARITY, and how it is broken (important):
##   PPP needs a definition of "final success", which is the gate Pr(p>p0)>c_final;
##   but c_final is itself calibrated using the stops that PPP produces -> circular.
##   Resolution: inside PPP, "success" is judged against a FIXED reference
##   c_ref = the single-stage c* (BHM 0.820 / EXNEX 0.852, from c_star_all.rds).
##   The ACTUAL final decision then uses a SEPARATELY calibrated c_final. So:
##     - PPP / futility question = "will it clear the usual efficacy bar (c_ref)?"
##     - final decision threshold c_final = calibrated for two-stage type-I = 0.10.
##   This decouples the loop cleanly and is documented as a design simplification of
##   full joint BOP2 (theta_fut, c_final) calibration.
##
## Efficiency: the R predictive fits feed an R-AVERAGED indicator, so per-fit MCMC
##   error averages out -> they use LIGHTER sampler settings (pred_chains/pred_iter).
##   The decision-relevant fits (interim, actual final) use full settings.
##
## Methods: BHM, EXNEX (the borrowing pair; no_borrow/pool omitted for cost -- add
##   later for the independent baseline). Scenarios: s1 (null: calibrates c_final +
##   shows futility CAN trigger) and s4 (heterogeneous: shows whether borrowing keeps
##   the dud alive). M=500, R=30.
##
## Divergences are captured into a table and the rds FROM THE START (interim + final
##   decision fits). Reuses locked bhm.stan / exnex.stan + c_star_all.rds.
## Writes: predprob_futility_results.rds
## ============================================================================

library(cmdstanr)
logit <- qlogis

## ---- run controls --------------------------------------------------------
DRY_RUN <- FALSE                  # TRUE -> tiny M/R, short chains: path check.
METHODS <- c("bhm", "exnex")
SCEN    <- c("s1", "s4")
theta_fut <- 0.05                # futility threshold on the predictive probability
n1 <- 10L; n2 <- 10L             # interim size, stage-2 increment (total 20)

M <- if (DRY_RUN)  20L else 500L     # simulated trials per scenario
R <- if (DRY_RUN)   5L else  30L     # posterior-predictive replicates per trial

## decision-relevant fits (interim, actual final): full settings
n_chains <- if (DRY_RUN) 2L else 4L
warmup   <- if (DRY_RUN) 300L else 1000L
sampling <- if (DRY_RUN) 300L else 1000L
## predictive fits (feed an R-averaged indicator): lighter
pred_chains <- if (DRY_RUN) 1L else 2L
pred_warmup <- if (DRY_RUN) 200L else 500L
pred_iter   <- if (DRY_RUN) 200L else 500L
max_td <- 12L

master_seed <- 8126L
target <- 0.10
PARALLEL  <- TRUE
n_workers <- max(1L, parallel::detectCores() - 1L)
cmdstan_dir <- cmdstanr::cmdstan_path()

## ---- locked design + priors ----------------------------------------------
baskets <- c("NSCLC", "Thyroid", "Biliary", "CRC")
K  <- length(baskets); p0 <- 0.10; p1 <- 0.30
cs <- readRDS("c_star_all.rds")
mu_mean <- cs$hyper_bhm$mu_mean; mu_sd <- cs$hyper_bhm$mu_sd; tau_scale <- cs$hyper_bhm$tau_scale
m_nex <- 0; s_nex <- 2; pi_val <- 0.5            # EXNEX locked design constants
c_ref <- list(bhm = cs$bhm, exnex = cs$exnex)    # fixed success bar INSIDE PPP

scenarios <- list(s1 = c(0.10,0.10,0.10,0.10),
                  s4 = c(0.40,0.30,0.15,0.05))
role_of <- function(pt) ifelse(pt <= p0+1e-9, "null", ifelse(pt >= p1-1e-9, "active", "gray"))

spec_of <- list(
  bhm   = list(file = normalizePath("bhm.stan"),   type = "bhm",   adapt_delta = 0.99),
  exnex = list(file = normalizePath("exnex.stan"), type = "exnex", adapt_delta = 0.99)
)
for (mn in METHODS) cmdstan_model(spec_of[[mn]]$file)

make_data <- function(type, nv, yv) {
  if (type == "bhm")   return(list(K=K, n=nv, y=yv, mu_mean=mu_mean, mu_sd=mu_sd,
                                   tau_scale=tau_scale, p0=p0))
  list(K=K, n=nv, y=yv, p0=p0, mu_mean=mu_mean, mu_sd=mu_sd, tau_scale=tau_scale,
       m_nex=m_nex, s_nex=s_nex, pi_ex=rep(pi_val, K))                        # exnex
}

## ---- one fit: returns Pr(p>p0) per basket, p draws, divergences -----------
fit_joint <- function(nv, yv, spec, seed, light = FALSE) {
  cmdstanr::set_cmdstan_path(cmdstan_dir)
  mod <- cmdstanr::cmdstan_model(spec$file)
  ch <- if (light) pred_chains else n_chains
  wu <- if (light) pred_warmup else warmup
  it <- if (light) pred_iter   else sampling
  fit <- mod$sample(data = make_data(spec$type, nv, yv), seed = seed,
                    chains = ch, parallel_chains = 1, iter_warmup = wu, iter_sampling = it,
                    adapt_delta = spec$adapt_delta, max_treedepth = max_td,
                    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  list(pi = colMeans(fit$draws("exceeds_p0", format = "draws_matrix")),
       pdraws = fit$draws("p", format = "draws_matrix"),
       ndiv = sum(fit$sampler_diagnostics(format = "draws_matrix")[, "divergent__"]))
}

## ---- predictive probability of final success for all baskets -------------
## nested MC: R joint future datasets drawn from the interim posterior predictive.
ppp_all <- function(y1, interim_pdraws, spec, cref, seed) {
  D <- nrow(interim_pdraws)
  idx <- sample.int(D, R, replace = (D < R))      # R posterior draws of p
  PrR <- matrix(NA_real_, R, K)
  ndv <- 0L
  for (r in seq_len(R)) {
    p_r  <- interim_pdraws[idx[r], ]
    y2   <- rbinom(K, n2, p_r)                    # future data, ALL baskets
    fin  <- fit_joint(rep(n1 + n2, K), y1 + y2, spec, seed + r, light = TRUE)
    PrR[r, ] <- fin$pi
    ndv <- ndv + fin$ndiv
  }
  list(ppp = colMeans(PrR > cref), ndiv_pred = ndv)   # PPP_k per basket
}

## ---- one full two-stage trial (predictive-prob futility) -----------------
one_trial <- function(y1, y2_real, spec, cref, seed) {
  intr <- fit_joint(rep(n1, K), y1, spec, seed)                       # interim [full]
  pp   <- ppp_all(y1, intr$pdraws, spec, cref, seed + 1000L)          # R predictive [light]
  stop <- pp$ppp < theta_fut                                          # futility
  yf <- y1 + ifelse(stop, 0L, y2_real)                               # frozen if stopped
  nf <- ifelse(stop, n1, n1 + n2)
  fin <- fit_joint(nf, yf, spec, seed + 2000L)                        # actual final [full]
  list(stop = stop, ppp = pp$ppp, pi_fin = fin$pi,
       ndiv_int = intr$ndiv, ndiv_fin = fin$ndiv, ndiv_pred = pp$ndiv_pred)
}

run_scen <- function(pt, spec, cref, seed_base) {
  set.seed(seed_base)
  Y1 <- t(replicate(M, rbinom(K, n1, pt)))        # interim data
  Y2 <- t(replicate(M, rbinom(K, n2, pt)))        # actual stage-2 data
  one <- function(m) one_trial(Y1[m,], Y2[m,], spec, cref, seed_base + 17L*m)
  if (PARALLEL && requireNamespace("furrr", quietly=TRUE) && requireNamespace("future", quietly=TRUE)) {
    library(furrr); library(future)
    plan(multisession, workers = n_workers); on.exit(plan(sequential), add = TRUE)
    res <- furrr::future_map(seq_len(M), one,
            .options = furrr::furrr_options(seed = TRUE, packages = "cmdstanr"))
    plan(sequential)
  } else res <- lapply(seq_len(M), one)
  list(STOP = do.call(rbind, lapply(res, `[[`, "stop")),
       PIF  = do.call(rbind, lapply(res, `[[`, "pi_fin")),
       PPP  = do.call(rbind, lapply(res, `[[`, "ppp")),
       ndiv_int = vapply(res, `[[`, numeric(1), "ndiv_int"),
       ndiv_fin = vapply(res, `[[`, numeric(1), "ndiv_fin"))
}

## ---- bracketed c_final calibration (stops are fixed; c_final is free) -----
calibrate_cfinal <- function(STOP, PIF, c_grid = seq(0.30, 0.99, by = 0.0025)) {
  ma <- sapply(c_grid, function(cc) mean(colMeans((!STOP) & (PIF > cc))))
  if (min(ma) > target || max(ma) < target) return(list(c_final = NA_real_, ceiling = max(ma)))
  i <- max(which(ma >= target))
  if (i == length(c_grid)) return(list(c_final = c_grid[i], boundary = TRUE))
  a1 <- ma[i]; a2 <- ma[i+1]
  list(c_final = c_grid[i] + (target-a1)*(c_grid[i+1]-c_grid[i])/(a2-a1), boundary = FALSE)
}

## ==========================================================================
## RUN
## ==========================================================================
draws_per_fit <- n_chains * sampling
cat(sprintf("MODE: %s | M=%d R=%d | theta_fut=%.2f | methods: %s | scen: %s | workers=%d\n",
            if (DRY_RUN) "DRY RUN" else "FULL", M, R, theta_fut,
            paste(METHODS, collapse=","), paste(SCEN, collapse=","), n_workers))

rows <- list(); div_rows <- list(); cfin <- list()
for (mn in METHODS) {
  spec <- spec_of[[mn]]; cref <- c_ref[[mn]]
  cat(sprintf("\n========== [%s]  c_ref=%.4f ==========\n", mn, cref))
  scen_runs <- list()
  for (si in SCEN) {
    t0 <- Sys.time()
    run <- run_scen(scenarios[[si]], spec, cref, master_seed + 100*match(si, names(scenarios)))
    scen_runs[[si]] <- run
    div_rows[[length(div_rows)+1]] <- data.frame(method=mn, scenario=si, stage="interim",
        pct_fits_div=100*mean(run$ndiv_int>0), draw_div_pct=100*sum(run$ndiv_int)/(M*draws_per_fit))
    div_rows[[length(div_rows)+1]] <- data.frame(method=mn, scenario=si, stage="final",
        pct_fits_div=100*mean(run$ndiv_fin>0), draw_div_pct=100*sum(run$ndiv_fin)/(M*draws_per_fit))
    cat(sprintf("  %s %.1f min | interim div %.1f%% | final div %.1f%%\n", si,
                as.numeric(difftime(Sys.time(), t0, units="mins")),
                100*mean(run$ndiv_int>0), 100*mean(run$ndiv_fin>0)))
  }
  ## calibrate c_final on s1 (must be in SCEN)
  if (!"s1" %in% SCEN) stop("s1 must be included to calibrate c_final.")
  cal <- calibrate_cfinal(scen_runs[["s1"]]$STOP, scen_runs[["s1"]]$PIF)
  cfin[[mn]] <- cal$c_final
  cat(sprintf("  c_final=%s (c_ref %.4f)%s\n",
              ifelse(is.na(cal$c_final),"NA",sprintf("%.4f",cal$c_final)), cref,
              if (is.na(cal$c_final)) " [0.10 not bracketed -> loosen theta_fut]" else ""))
  if (is.na(cal$c_final)) next
  for (si in SCEN) {
    run <- scen_runs[[si]]; pt <- scenarios[[si]]
    Go <- colMeans((!run$STOP) & (run$PIF > cal$c_final))
    stp <- colMeans(run$STOP); ppp <- colMeans(run$PPP)
    ess <- n1 + (n2)*(1 - stp)
    for (k in seq_len(K))
      rows[[length(rows)+1]] <- data.frame(
        method=mn, scenario=si, basket=baskets[k], p_true=pt[k], role=role_of(pt[k]),
        c_final=cal$c_final, Go=Go[k], stop_interim=stp[k], mean_PPP=ppp[k], ESS=ess[k])
  }
}
OC <- do.call(rbind, rows); DIV <- do.call(rbind, div_rows); rownames(DIV) <- NULL
saveRDS(list(OC=OC, DIV=DIV, c_final=cfin, c_ref=c_ref, theta_fut=theta_fut,
             M=M, R=R, dry_run=DRY_RUN),
        file = "predprob_futility_results.rds")
cat("\nsaved predprob_futility_results.rds (OC, DIV, c_final, c_ref)\n")

## ==========================================================================
## SUMMARY
## ==========================================================================
as_mat <- function(sub, col) {
  m <- matrix(NA_real_, length(METHODS), K, dimnames=list(METHODS, baskets))
  for (i in seq_len(nrow(sub))) m[sub$method[i], sub$basket[i]] <- sub[[col]][i]; m
}
for (si in SCEN) {
  sub <- OC[OC$scenario==si, ]; if (!nrow(sub)) next
  cat(sprintf("\n================ %s (true ORR: %s) ================\n", si,
              paste(sprintf("%.2f", scenarios[[si]]), collapse="/")))
  cat("final Go freq (active->power | null->type I | gray->Go-rate):\n"); print(round(as_mat(sub,"Go"),3))
  cat("\ninterim futility-stop freq (PREDICTIVE-prob rule):\n");           print(round(as_mat(sub,"stop_interim"),3))
  cat("\nmean predictive prob (PPP) per basket:\n");                       print(round(as_mat(sub,"mean_PPP"),3))
  cat("\nexpected N per basket (max 20):\n");                              print(round(as_mat(sub,"ESS"),1))
}

cat("\n================ DIVERGENCES ================\n")
print(within(DIV, {pct_fits_div<-round(pct_fits_div,1); draw_div_pct<-round(draw_div_pct,3)}), row.names=FALSE)

cat("\n================ HEADLINE ================\n")
cat("Recalibrated c_final vs c_ref (single-stage c*):\n")
for (mn in names(cfin)) if(!is.na(cfin[[mn]]))
  cat(sprintf("  %-7s c_final=%.4f  (c_ref=%.4f, delta %+.4f)\n", mn, cfin[[mn]], c_ref[[mn]], cfin[[mn]]-c_ref[[mn]]))

cat("\nDoes predictive-prob futility now TRIGGER? (s1, all null -> want non-trivial stops)\n")
for (mn in names(cfin)) {
  s1 <- OC[OC$scenario=="s1" & OC$method==mn, ]
  if (nrow(s1)) cat(sprintf("  %-7s s1 interim-stop (avg) = %.3f | avg ESS = %.1f / 20\n",
                            mn, mean(s1$stop_interim), mean(s1$ESS)))
}
cat("\nDoes BORROWING keep the dud alive? (s4 CRC, truly null 0.05)\n")
for (mn in names(cfin)) {
  crc <- OC[OC$scenario=="s4" & OC$basket=="CRC" & OC$method==mn, ]
  if (nrow(crc)) cat(sprintf("  %-7s CRC: interim-stop=%.3f  PPP=%.3f  final type-I=%.3f  ESS=%.1f\n",
                            mn, crc$stop_interim, crc$mean_PPP, crc$Go, crc$ESS))
}
cat("  -> compare s1 stop (no neighbours) vs s4 CRC stop (active neighbours lift PPP via borrowing).\n")
cat(if (DRY_RUN) "\n*** DRY RUN: tiny M/R, numbers are noise; checking the pipeline only. ***\n" else "")
