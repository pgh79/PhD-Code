functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data{
  real t0;    // Initial time in days;
  real C0[1]; // Initial concentration at t0 in mg/L
  real D;   // Total dosage in mg.  2.5 mg of Apixaban

  int<lower=1> N_t;
  int<lower=1> N_patients;
  real times[N_t];   // Measurement times in days
  matrix[N_patients,4] X;
  
  real<lower=0> A;
  
  
  //Use these for generating samples at time unique times
  
  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_patients,N_t];
}
parameters {

  vector[4] beta_v;
  vector[4] beta_k;
  vector[4] beta_ka;
  real<lower=0> sigma;
  real<lower=0, upper=1> alpha;
}
transformed parameters {

  // Store analytic concentration
  real<lower=0> C[N_patients,N_t];
  vector[N_patients] V;
  vector[N_patients] k;
  vector[N_patients] k_a;
  
  V = exp(X*beta_v);
  k = exp(X*beta_k);
  k_a = exp(X*beta_ka);
  // REMOVE THE ZEROS FOR INTERACTION
  for (p in 1:N_patients){
    for (n in 1:N_t)
      C[p,n] = C_anal(times[n], D, V[p], k_a[p], k[p]);
  }
}

model {
  
  
  beta_v ~ normal(0,2);
  beta_k ~ normal(0,2);
  beta_ka~ normal(0,2);
  alpha ~ uniform(0,1);
  
  
  sigma ~ cauchy(0,1);
  // Likelihood
  for (p in 1:N_patients){
   for (n in 1:N_t)
    C_hat[p,n] ~ normal(C[p,n], pow(C[p,n]/A,2*alpha)*sigma); 
  }
}

generated quantities {
  real C_ppc[N_patients,N_t];
  for (p in 1:N_patients){
   for (n in 1:N_t)
    C_ppc[p,n] = normal_rng(C[p,n], pow(C[p,n]/A,2*alpha)*sigma); 
  }

}


