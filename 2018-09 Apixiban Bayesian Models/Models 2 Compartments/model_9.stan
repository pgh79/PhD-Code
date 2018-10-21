functions {
  real[] two_comp(real t,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
    real dydt[2];
    real k_a = theta[1];  
    real k = theta[2];  
    real V = theta[3];  // Maximum elimination rate in 1/day
    real k_12 = theta[4]; // Compartment exchange rate in 1/day
    real k_21 = theta[5]; // Compartment exchange rate in 1/day
    real D = x_r[1];
    real dose = 0;
    real elim =  -(k_12 +k)*y[1] +k_21*y[2];

    if (t > 0)
      dose = exp(- k_a * t) * D * k_a / V;

    dydt[1] = dose + elim;
    dydt[2] = k_12*y[1] - k_21*y[2];

    return dydt;
  }
}

data {
  real t0;    // Initial time in days;
  real C0[2]; // Initial concentration at t0 in mg/L

  real D;   // Total dosage in mg

  int<lower=1> N_t;
  int<lower=1> N_patients;
  real times[N_t];   // Measurement times in days


  matrix[N_patients, 4] X;
  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_patients,N_t];
}

transformed data {
  real x_r[1] = {D};
  int x_i[0];
}

parameters {
  vector[4] beta_v;
  vector[4] beta_k;
  vector[4] beta_ka;
  vector[4] beta_k12; 
  vector[4] beta_k21; 
  real<lower=0> sigma;
}

transformed parameters {
  real C_sol[N_t,2];
  real C[N_patients, N_t];
  vector[N_patients] V;
  vector[N_patients] k;
  vector[N_patients] k_a;
  vector[N_patients] k_12;
  vector[N_patients] k_21;
  
  V = exp(X*beta_v);
  k = exp(X*beta_k);
  k_a = exp(X*beta_ka);
  k_12 = exp(X*beta_k12);
  k_21 = exp(X*beta_k21);
  {
    for (p in 1:N_patients){
    real theta[5] = {k_a[p], k[p], V[p], k_12[p], k_21[p]};
    C_sol= integrate_ode_bdf(two_comp, C0, t0, times, theta, x_r, x_i);
    C[p,1:N_t] = C_sol[1:N_t,1];
  }
}
}

model {
  // Priors
  
  beta_ka ~ normal(0, 1);
  beta_k ~ normal(0, 1);
  beta_v ~ normal(0,1);
  beta_k12 ~ normal(0,1);
  beta_k21 ~ normal(0,1);
  sigma ~ cauchy(0, 1);

  // Likelihood
  for (p in 1:N_patients){
  for (n in 1:N_t)
    C_hat[p,n] ~ lognormal(log(C[p,n]),sigma);
    }
}

generated quantities {
  real C_ppc[N_patients,N_t];
  for (p in 1:N_patients){
   for (n in 1:N_t)
    C_ppc[p,n] = lognormal_rng(log(C[p,n]), sigma);
  }
}
