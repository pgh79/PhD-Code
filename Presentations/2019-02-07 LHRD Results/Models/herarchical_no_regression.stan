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
  int<lower=1> N_patients;

  real times[N_t];   // Measurement times in days

  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_patients,N_t];
}

parameters {
  real<lower=0> MU_k_a; 
  real<lower=0> SIGMA_k_a;
  real<lower=0> z_k_a[N_patients]; // Dosing rate in 1/day
  
  
  real<lower=0> MU_k; 
  real<lower=0> SIGMA_k;
  real<lower=0> z_k[N_patients];   // Elimination rate in 1/day

  real<lower=0> MU_V; 
  real<lower=0> SIGMA_V;
  real<lower=0> z_V[N_patients]; // Volume of the absorptive compartment
  
  real<lower=0> sigma;
}

transformed parameters{
  
  real C[N_patients, N_t];
  real k[N_patients];
  real k_a[N_patients];
  real V[N_patients];
  
  for (n in 1:N_patients){
    
    k[n] = exp(MU_k + SIGMA_k*z_k[n]);
    k_a[n] = exp(MU_k_a + SIGMA_k_a*z_k_a[n]);
    V[n] = exp(MU_V + SIGMA_V*z_V[n]);
    
    for (t in 1:N_t){
      C[n,t] = C_anal(times[t], D, V[n], k_a[n], k[n]);
    }
  }
  
}

model {
  // Priors
  sigma ~ cauchy(0, 1);
  
  MU_k_a ~ normal(0, 1);
  SIGMA_k_a ~ cauchy(0, 1);
  z_k_a ~ normal(0,1);
  
  MU_k ~ normal(0, 1);
  SIGMA_k ~ cauchy(0, 1);
  z_k ~ normal(0,1);
  
  MU_V ~ normal(0, 1);
  SIGMA_V ~ cauchy(0, 1);
  z_V ~ normal(0,1);

  // Likelihood
    for (n in 1:N_patients){
      for (t in 1:N_t){
        C_hat[n,t] ~ lognormal(log(C_anal(times[t], D, V[n], k_a[n], k[n])), sigma );
    }
  }
}

generated quantities {
  real C_ppc[N_patients, N_t];
  vector[N_t*N_patients] log_lik;
  for (n in 1:N_patients){
    for (t in 1:N_t){
      C_ppc[n, t] = lognormal_rng(log(C[n, t]), sigma );
      log_lik[N_t*(n-1) + t] = lognormal_lpdf( C_hat[n, t] | log(C[n,t]), sigma);
    }
  }
}
