// No borrowing (lower bracket): independent Beta-Binomial per basket.
// Each p_k estimated alone; no information shared across baskets.
data {
  int<lower=1> K;
  array[K] int<lower=0> n;
  array[K] int<lower=0> y;
  real<lower=0> a0;                 // Beta prior shape (Jeffreys: a0 = b0 = 0.5)
  real<lower=0> b0;
  real<lower=0, upper=1> p0;        // null ORR for the decision quantity
}
parameters {
  vector<lower=0, upper=1>[K] p;    // one independent rate per basket
}
model {
  p ~ beta(a0, b0);
  y ~ binomial(n, p);
}
generated quantities {
  array[K] int exceeds_p0;
  for (k in 1:K) exceeds_p0[k] = (p[k] > p0);
  // posterior mean of exceeds_p0[k] == Pr(p_k > p0 | y_k)  (per basket, no borrowing)
}
