// EXNEX for PRIOR SENSITIVITY (step 5): identical to exnex.stan EXCEPT the EX-component
// tau prior family is switchable via data (same switch as bhm_psens.stan). This one
// file covers THREE sensitivity axes with no further code changes:
//   - tau prior family / scale : via tau_prior_family, tau_scale, ig_a, ig_b
//   - fixed mixing weight pi    : via the data vector pi_ex (e.g. 0.3 / 0.5 / 0.7)
//   - NEX center / scale        : via m_nex, s_nex (e.g. m_nex = 0 vs logit(0.10))
// The Beta(1,1) RANDOM-pi variant is a different model (pi becomes a parameter) and
// lives in exnex_pi_random.stan.
//
//   tau_prior_family = 1 Half-Normal | 2 Half-Cauchy | 3 Inverse-Gamma(on tau^2)
// (See bhm_psens.stan header for the family math and the IG Jacobian.)
// theta is sampled directly (centered) as in exnex.stan; the discrete EX/NEX
// indicator is marginalized via log_mix. Run with adapt_delta 0.99.

data {
  int<lower=1> K;
  array[K] int<lower=0> n;
  array[K] int<lower=0> y;
  real<lower=0, upper=1> p0;

  real mu_mean;
  real<lower=0> mu_sd;
  int<lower=1, upper=3> tau_prior_family;  // EX-component tau prior
  real<lower=0> tau_scale;
  real<lower=0> ig_a;
  real<lower=0> ig_b;

  real m_nex;                              // NEX mean  (logit scale): 0 or logit(0.10)
  real<lower=0> s_nex;                     // NEX sd    (logit scale)
  vector<lower=0, upper=1>[K] pi_ex;       // fixed EX prior weight per basket
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[K] theta;
}

model {
  mu ~ normal(mu_mean, mu_sd);

  if (tau_prior_family == 1)
    tau ~ normal(0, tau_scale);
  else if (tau_prior_family == 2)
    tau ~ cauchy(0, tau_scale);
  else
    target += inv_gamma_lpdf(square(tau) | ig_a, ig_b) + log(2) + log(tau);

  for (k in 1:K)
    target += log_mix(pi_ex[k],
                      normal_lpdf(theta[k] | mu,    tau),
                      normal_lpdf(theta[k] | m_nex, s_nex));

  y ~ binomial_logit(n, theta);
}

generated quantities {
  vector<lower=0, upper=1>[K] p = inv_logit(theta);
  array[K] int exceeds_p0;
  vector<lower=0, upper=1>[K] w_ex;
  for (k in 1:K) {
    exceeds_p0[k] = (p[k] > p0);
    real lex  = log(pi_ex[k])  + normal_lpdf(theta[k] | mu,    tau);
    real lnex = log1m(pi_ex[k]) + normal_lpdf(theta[k] | m_nex, s_nex);
    w_ex[k]   = exp(lex - log_sum_exp(lex, lnex));
  }
}
