functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}

data {
  real t0;    // Initial time in days;
  real C0[1]; // Initial concentration at t0 in mg/L

  real D;   // Total dosage in mg

  int<lower=1> N_t;

  real times[N_t];   // Measurement times in days

  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_t];
}

parameters {
  real<lower=0> k_a; // Dosing rate in 1/day
  real<lower=0> k;   // Elimination rate in 1/day
  real<lower=0> sigma;
  real<lower=0> V; // Volume of the absorptive compartment
}

model {
  // Priors
  k_a ~ cauchy(0, .1);
  k ~ cauchy(0, .1);
  sigma ~ cauchy(0, 1);
  V ~ cauchy(1, 1); // 0.1 Litres, or 100 millileters

  // Likelihood
  for (n in 1:N_t)
    C_hat[n] ~ lognormal(log(C_anal(times[n], D, V, k_a, k)), sigma);
}

generated quantities {
  real C_ppc[N_t];
  real C_pred[N_t];
  vector[N_t] log_lik;
  for (n in 1:N_t) {
    C_pred[n] = C_anal(times[n], D, V, k_a, k);
    C_ppc[n] = lognormal_rng(log(C_anal(times[n], D, V, k_a, k)), sigma);
    log_lik[n] = lognormal_lpdf(C_hat[n] | log(C_anal(times[n], D, V, k_a, k)), sigma);
  }
}
