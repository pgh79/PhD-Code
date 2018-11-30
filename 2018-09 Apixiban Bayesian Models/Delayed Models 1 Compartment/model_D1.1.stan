// Same as model 1.7 except with a normal likelihood

functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data {
  real D;   // Total dosage in mg

  int<lower=1> N_t;
  real times[N_t];   // Measurement times in days

  
  real<lower=0> A;
  
  int<lower=1>p;


  // Measured concentrations in effect compartment in mg/L
  int<lower=1> N_patients;
  real C_hat[N_patients, N_t];
  matrix[N_patients,p] X;
}
parameters {
  
  //Volume
  vector[p] BETA_V;
  real<lower=0> SIGMA_V;
  real z_V[N_patients];
  
  //Elimination Rate
  vector[p] BETA_K;
  real<lower=0> SIGMA_K;
  real z_k[N_patients];
  
  //Absorption
  vector[p] BETA_KA;
  real<lower=0> SIGMA_KA;
  real z_ka[N_patients];
  
  real<lower=0>sigma;
  real<lower=0, upper=1> alpha;
}

transformed parameters {

  real C[N_patients, N_t];
  real<lower=0> k_a[N_patients];
  real<lower=0> k[N_patients];
  real<lower=0> V[N_patients];
  
  vector[N_patients] MU_KA;
  vector[N_patients] MU_K;
  vector[N_patients] MU_V;

  MU_KA = X*BETA_KA;
  MU_K = X*BETA_K;
  MU_V = X*BETA_V;
  
  
  for (n in 1:N_patients) {
    k_a[n] = exp(MU_KA[n] + z_ka[n]*SIGMA_KA);
    k[n] = exp(MU_K[n] + z_k[n]*SIGMA_K);
    V[n] = exp(MU_V[n] + z_V[n]*SIGMA_V);
    
    for (t in 1:N_t)
      C[n,t] = C_anal(times[t], D, V[n], k_a[n], k[n]);

  }
}

model {
  // Priors

  BETA_V ~ normal(0,1);
  SIGMA_V ~ cauchy(0,1);
  z_V ~ normal(0,1);
  
  BETA_K ~ normal(0,1);
  SIGMA_K ~ cauchy(0,1);
  z_k ~ normal(0,1);

  BETA_KA ~ normal(0,1);
  SIGMA_KA ~ cauchy(0,1);
  z_ka ~ normal(0,1);

  sigma ~ cauchy(0,1);
  alpha ~ uniform(0,1);

  // Likelihood
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_hat[n, t] ~ normal(C[n, t], pow(C[n,t]/A,2*alpha)*sigma);
      // C_hat[n, t] ~ normal(C[n, t], sigma);
}

generated quantities {
  real C_ppc[N_patients, N_t];
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_ppc[n, t] = normal_rng(C[n, t], pow(C[n,t]/A,2*alpha)*sigma);
      // C_ppc[n, t] = normal_rng(C[n, t], sigma);
}
