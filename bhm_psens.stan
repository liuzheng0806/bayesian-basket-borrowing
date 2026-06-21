// BHM for PRIOR SENSITIVITY (step 5): identical to bhm.stan EXCEPT the tau prior
// family is switchable via data, so the tau-prior-FAMILY sweep (the project's
// selling point) runs off one compiled model. Non-centered theta = mu + tau*z is
// kept for all families (the funnel mitigation concerns theta|mu,tau, independent
// of tau's marginal prior).
//
//   tau_prior_family = 1 : tau ~ Half-Normal(0, tau_scale)     (== bhm.stan; main analysis)
//                      2 : tau ~ Half-Cauchy(0, tau_scale)     (Gelman 2006 recommendation)
//                      3 : tau^2 ~ Inverse-Gamma(ig_a, ig_b)   (cautionary; Gelman 2006 warns)
//
// For family 3 the prior is placed on tau^2, so the Jacobian for tau^2 -> tau is
// added explicitly: log|d(tau^2)/dtau| = log(2*tau). (Verified: the resulting tau
// density integrates to 1.) Half-Normal / Half-Cauchy normalizing constants are
// dropped (constant in the posterior; same as the locked bhm.stan).
//
// All hyperprior inputs stay DATA. ig_a/ig_b/tau_scale are always supplied (unused
// ones are dummies) because Stan requires every declared data variable.

data {
  int<lower=1> K;
  array[K] int<lower=0> n;
  array[K] int<lower=0> y;
  real mu_mean;
  real<lower=0> mu_sd;
  real<lower=0, upper=1> p0;

  int<lower=1, upper=3> tau_prior_family;  // 1 HN | 2 HC | 3 IG(on tau^2)
  real<lower=0> tau_scale;                 // scale for HN (fam 1) or HC (fam 2)
  real<lower=0> ig_a;                      // IG shape (fam 3)
  real<lower=0> ig_b;                      // IG scale (fam 3)
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[K] z;
}

transformed parameters {
  vector[K] theta = mu + tau * z;          // non-centered (as bhm.stan)
}

model {
  mu ~ normal(mu_mean, mu_sd);
  z  ~ std_normal();

  if (tau_prior_family == 1)
    tau ~ normal(0, tau_scale);            // Half-Normal (<lower=0> truncates)
  else if (tau_prior_family == 2)
    tau ~ cauchy(0, tau_scale);            // Half-Cauchy (<lower=0> folds)
  else
    target += inv_gamma_lpdf(square(tau) | ig_a, ig_b)
              + log(2) + log(tau);         // tau^2 ~ IG, with Jacobian for tau

  y ~ binomial_logit(n, theta);
}

generated quantities {
  vector<lower=0, upper=1>[K] p = inv_logit(theta);
  array[K] int exceeds_p0;
  for (k in 1:K) exceeds_p0[k] = (p[k] > p0);
}
