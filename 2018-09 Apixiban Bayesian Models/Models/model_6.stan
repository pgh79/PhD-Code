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
  int<lower=1> N_patients;
  real C_hat[N_patients, N_t];
}
parameters {
  real<lower=0> k[N_patients];
  real<lower=0> k_a[N_patients];
  real<lower=0> V[N_patients];
  real<lower=0>sigma;
}

transformed parameters {

  real C[N_patients, N_t];
  for (n in 1:N_patients) {
    for (t in 1:N_t)
      C[n,t] = C_anal(times[t], D, V[n], k_a[n], k[n]);
  }
}

model {
  // Priors

  k ~ cauchy(0,1);
  k_a ~ cauchy(0,1);
  V ~ cauchy(0,1);
  sigma ~ cauchy(0, 1);

  // Likelihood
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_hat[n, t] ~ lognormal(log(C[n, t]) - sigma/2, sigma);
}

generated quantities {
  real C_ppc[N_patients, N_t];
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_ppc[n, t] = lognormal_rng(log(C[n, t]) - sigma/2, sigma);
}
