// Basket-trial BHM, step 1: two-level logit-normal with partial borrowing.
//   y_k ~ Binomial(n_k, p_k)
//   logit(p_k) = theta_k,  theta_k ~ Normal(mu, tau^2)
//   tau ~ Half-Normal(0, s)         (between-basket SD, logit scale)
//   mu  ~ Normal(mu_mean, mu_sd)    (population mean, logit scale)
//
// Non-centered parameterization (theta = mu + tau*z) to avoid the
// Neal's-funnel / divergence pathology that bites HMC when K is small
// and tau is near 0 -- exactly this regime. This is why HMC+reparam is
// preferable to Gibbs here.

data {
  int<lower=1> K;                       // number of baskets
  array[K] int<lower=0> n;              // patients per basket
  array[K] int<lower=0> y;              // responders per basket
  // hyperprior inputs kept as DATA so nothing is silently hard-coded:
  real mu_mean;                         // prior mean for mu  (logit scale)
  real<lower=0> mu_sd;                  // prior sd   for mu  (logit scale)
  real<lower=0> tau_scale;              // s in tau ~ Half-Normal(0, s)
  real<lower=0, upper=1> p0;            // null ORR for the decision quantity
}

parameters {
  real mu;                              // population mean log-odds
  real<lower=0> tau;                    // between-basket SD (logit scale)
  vector[K] z;                          // standardized basket offsets
}

transformed parameters {
  vector[K] theta = mu + tau * z;       // basket log-odds
}

model {
  mu  ~ normal(mu_mean, mu_sd);
  tau ~ normal(0, tau_scale);           // half-normal: <lower=0> truncates
  z   ~ std_normal();                   // => theta_k ~ N(mu, tau^2)
  y   ~ binomial_logit(n, theta);       // numerically stable likelihood
}

generated quantities {
  vector<lower=0, upper=1>[K] p = inv_logit(theta);  // basket ORRs
  array[K] int exceeds_p0;                           // 1 if p_k > p0 this draw
  for (k in 1:K) exceeds_p0[k] = (p[k] > p0);
  // posterior mean of exceeds_p0[k] over draws == Pr(p_k > p0 | data)
}
