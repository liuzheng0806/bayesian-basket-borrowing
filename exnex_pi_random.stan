// EXNEX with RANDOM mixing weight (step-5 sensitivity): pi_k is ESTIMATED, not fixed.
//   pi_k ~ Beta(a_pi, b_pi)    (a_pi = b_pi = 1  ->  Beta(1,1) = Uniform(0,1), the spec'd variant)
// Everything else matches exnex_psens.stan (same EX/NEX mixture via log_mix, same
// switchable EX-component tau prior). The question this answers: does letting the
// data estimate the EX/NEX weight change the OC vs fixing pi = 0.5? With n = 20 per
// basket pi_k is weakly identified, so a near-flat pi posterior is itself a finding.
//
// w_ex now uses the SAMPLED pi_ex[k] per draw, so its posterior mean is the posterior
// EX weight under a learned (rather than fixed) prior weight.

data {
  int<lower=1> K;
  array[K] int<lower=0> n;
  array[K] int<lower=0> y;
  real<lower=0, upper=1> p0;

  real mu_mean;
  real<lower=0> mu_sd;
  int<lower=1, upper=3> tau_prior_family;
  real<lower=0> tau_scale;
  real<lower=0> ig_a;
  real<lower=0> ig_b;

  real m_nex;
  real<lower=0> s_nex;
  real<lower=0> a_pi;                       // Beta hyperprior shape (1 for Beta(1,1))
  real<lower=0> b_pi;
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[K] theta;
  vector<lower=0, upper=1>[K] pi_ex;        // EX weight, now a PARAMETER
}

model {
  mu ~ normal(mu_mean, mu_sd);

  if (tau_prior_family == 1)
    tau ~ normal(0, tau_scale);
  else if (tau_prior_family == 2)
    tau ~ cauchy(0, tau_scale);
  else
    target += inv_gamma_lpdf(square(tau) | ig_a, ig_b) + log(2) + log(tau);

  pi_ex ~ beta(a_pi, b_pi);                 // estimate the mixing weights

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
