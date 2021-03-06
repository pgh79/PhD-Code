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
  int sex[N_t];
  
  int<lower=1> N_ut;
  real utimes[N_ut];
  
  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_t];
}


parameters {
  real<lower=0> k_a; // Dosing rate in 1/day
  real<lower=0> k;

  real beta_V[2];
  

  real<lower=0> sigma;
}

transformed parameters {

  real C[N_t];
  real V[N_t];
  real C_pred[2,N_ut];
  real V_ppc;
  for (n in 1:N_t) 
    V[n] = exp(beta_V[1] + sex[n]*beta_V[2]);
  for (n in 1:N_t)
    C[n] = C_anal(times[n], D, V[n], k_a, k);
    
  for (i in 1:2){
    V_ppc = exp(beta_V[1]+ (i-1)*beta_V[2]);
    for (n in 1:N_ut)    
      C_pred[i,n] = C_anal(utimes[n], D, V_ppc, k_a, k);
  }
}

model {
  // Priors
  k_a ~ cauchy(0, 1);
  k ~ cauchy(0, 1);
  beta_V ~ normal(0,1);
  sigma ~ cauchy(0, 1);

  // Likelihood
  for (n in 1:N_t)
    C_hat[n] ~ lognormal(log(C[n]) - sigma/2, sigma);
    // C_hat[n] ~ normal(C[n], sigma);
}

generated quantities {
  real C_ppc[2,N_ut];
  for (i in 1:2){
    for (n in 1:N_ut)
      C_ppc[i,n] = lognormal_rng(log(C_pred[i,n])- sigma/2, sigma);
    // C_ppc[i,n] = normal_rng(C_pred[i,n], sigma);
  }
}
