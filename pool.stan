// Complete pooling (upper bracket / strawman): ONE shared response rate p
// across all baskets. The decision is identical for every basket.
data {
  int<lower=1> K;
  array[K] int<lower=0> n;
  array[K] int<lower=0> y;
  real<lower=0> a0;                 // Beta prior shape (Jeffreys: a0 = b0 = 0.5)
  real<lower=0> b0;
  real<lower=0, upper=1> p0;
}
parameters {
  real<lower=0, upper=1> p;         // single shared rate
}
model {
  p ~ beta(a0, b0);
  for (k in 1:K) y[k] ~ binomial(n[k], p);
}
generated quantities {
  int exceeds_p0 = (p > p0);        // ONE decision; applied to every basket in R
  // posterior mean of exceeds_p0 == Pr(p > p0 | all data)
}
