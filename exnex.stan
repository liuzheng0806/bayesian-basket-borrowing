// Basket-trial EXNEX (column 4): robust EXchangeable / Non-EXchangeable mixture.
// Neuenschwander, Wandel, Roychoudhury & Bailey (2016), Pharm Stat 15(2):123-134.
//
// Each basket's log-odds theta_k is a 2-component mixture:
//   with prior weight   pi_k :  EX  -- theta_k ~ Normal(mu, tau^2)      (shared, = BHM)
//   with prior weight 1-pi_k :  NEX -- theta_k ~ Normal(m_nex, s_nex^2) (independent, vague)
// The EX component reuses the BHM hierarchy verbatim (same mu, tau and their priors);
// the NEX component is a wide stratum-specific prior that lets a basket "opt out" of
// borrowing when its data conflict with the rest. That opt-out is the whole point:
// it is meant to protect a truly-inactive basket (scenario-4 CRC) from being dragged
// up by the active baskets the way blind BHM borrowing does.
//
// The discrete EX/NEX indicator is MARGINALIZED OUT (HMC cannot sample discrete
// parameters): the prior contribution of each theta_k is the log of the 2-component
// mixture density, added with log_mix. Because of this marginalization theta_k is a
// DIRECTLY-sampled (centered) parameter -- the non-centered reparameterization used in
// bhm.stan does not carry over to a mixture prior, so a mild EX-component funnel can
// remain. Run with high adapt_delta (0.99) and monitor divergences.
//
// Likelihood is binomial_logit (numerically stable), identical to bhm.stan, so the
// decision quantity p / exceeds_p0 is produced exactly as in the other three models
// and EXNEX drops straight into step3b's fit_one() as a non-shared method.
//
// generated quantities also returns the per-draw EX responsibility w_ex[k] =
// Pr(basket k came from EX | theta_k, mu, tau). Its posterior MEAN is the per-basket
// "posterior EX weight": the core EXNEX diagnostic. It slides toward 0 for a basket
// whose data conflict with the others (that basket self-estimates via NEX) and stays
// near pi_k for a basket consistent with the pool.

data {
  int<lower=1> K;                          // number of baskets
  array[K] int<lower=0> n;                 // patients per basket
  array[K] int<lower=0> y;                 // responders per basket
  real<lower=0, upper=1> p0;               // null ORR for the decision quantity

  // EX hyperprior (IDENTICAL to bhm.stan; kept as data so nothing is hard-coded):
  real mu_mean;                            // prior mean for mu (logit scale)
  real<lower=0> mu_sd;                     // prior sd   for mu (logit scale)
  real<lower=0> tau_scale;                 // s in tau ~ Half-Normal(0, s)

  // NEX prior (wide, stratum-specific; design spec: Normal(0, 2) on the logit scale,
  // i.e. +/-2SD ~ ORR 0.02-0.98 -- a deliberately vague "own estimate" prior):
  real m_nex;                              // NEX mean (logit scale), e.g. 0
  real<lower=0> s_nex;                     // NEX sd   (logit scale), e.g. 2

  // mixing weights, per basket so sensitivity can vary them (main analysis = 0.5):
  vector<lower=0, upper=1>[K] pi_ex;       // prior Pr(basket k is EX)
}

parameters {
  real mu;                                 // EX population mean log-odds  (= BHM mu)
  real<lower=0> tau;                       // EX between-basket SD, logit  (= BHM tau)
  vector[K] theta;                         // basket log-odds (sampled directly; see header)
}

model {
  // --- EX hyperpriors (same as BHM) ---
  mu  ~ normal(mu_mean, mu_sd);
  tau ~ normal(0, tau_scale);              // half-normal via the <lower=0> on tau

  // --- mixture prior on each theta_k (discrete EX/NEX indicator marginalized) ---
  // log_mix(w, a, b) = log( w*exp(a) + (1-w)*exp(b) ); first arg weights the EX piece.
  for (k in 1:K)
    target += log_mix(pi_ex[k],
                      normal_lpdf(theta[k] | mu,    tau),      // EX  component
                      normal_lpdf(theta[k] | m_nex, s_nex));   // NEX component

  // --- likelihood (identical to bhm.stan) ---
  y ~ binomial_logit(n, theta);            // numerically stable
}

generated quantities {
  vector<lower=0, upper=1>[K] p = inv_logit(theta);   // basket ORRs
  array[K] int exceeds_p0;                            // 1 if p_k > p0 this draw
  vector<lower=0, upper=1>[K] w_ex;                   // EX responsibility this draw
  for (k in 1:K) {
    exceeds_p0[k] = (p[k] > p0);
    // posterior Pr(EX | theta_k, mu, tau) for this draw, computed in log space:
    real lex  = log(pi_ex[k])   + normal_lpdf(theta[k] | mu,    tau);
    real lnex = log1m(pi_ex[k])  + normal_lpdf(theta[k] | m_nex, s_nex);
    w_ex[k]   = exp(lex - log_sum_exp(lex, lnex));
  }
  // posterior mean of exceeds_p0[k] == Pr(p_k > p0 | data)  (decision quantity, as BHM)
  // posterior mean of w_ex[k]       == per-basket posterior EX weight (EXNEX diagnostic)
}
