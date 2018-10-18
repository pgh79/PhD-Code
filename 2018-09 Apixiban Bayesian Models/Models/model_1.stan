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
  int<lower=1> N_ut;
  real times[N_t];   // Measurement times in days
  real utimes[N_ut];

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
  k_a ~ cauchy(0, 100);
  k ~ cauchy(0, 100);
  sigma ~ cauchy(0, 100);
  V ~ normal(0.1, 0.075); // 0.1 Litres, or 100 millileters

  // Likelihood
  for (n in 1:N_t)
    // C_hat[n] ~ normal(C_anal(times[n], D, V, k_a, k), sigma);
    C_hat[n] ~ lognormal(log(C_anal(times[n], D, V, k_a, k)) - sigma/2, sigma);
}

generated quantities {
  real C[N_ut];
  real C_ppc[N_ut];
  for (n in 1:N_ut) {
    C[n] = C_anal(utimes[n], D, V, k_a, k);
    C_ppc[n] = lognormal_rng(log(C[n]) - sigma/2, sigma);
  }
}
