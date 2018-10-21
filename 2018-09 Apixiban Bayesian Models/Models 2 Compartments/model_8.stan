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
  real times[N_t];   // Measurement times in days
  
  int<lower=1> N_it;
  int itimes[N_it];

  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_it];
}

transformed data {
  real x_r[1] = {D};
  int x_i[0];
}

parameters {
  real<lower=0> k_a;  
  real<lower=0> k;  
  real<lower=0> V;  
  real<lower=0> k_12; 
  real<lower=0> k_21; 
  real<lower=0> sigma;
}

transformed parameters {
  real C[N_t, 2];
  {
    real theta[5] = {k_a, k, V, k_12, k_21};
    C = integrate_ode_bdf(two_comp, C0, t0, times, theta, x_r, x_i);
  }
}

model {
  // Priors
  
  k_a ~ cauchy(0, 1);
  k ~ cauchy(0, 1);
  V ~ normal(0.250, 0.075);
  k_12 ~ cauchy(0, 1);
  k_21 ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);

  // Likelihood
  for (n in 1:N_it){
    C_hat[n] ~ lognormal(log(C[itimes[n],1]),sigma);
    }
}

generated quantities {
  real C_pred[N_t];
  real C_ppc[N_t];
  for (n in 1:N_t){
    C_pred[n] = C[n,1];
    C_ppc[n] = lognormal_rng(log(C[n, 1]), sigma) ;
  }
}
