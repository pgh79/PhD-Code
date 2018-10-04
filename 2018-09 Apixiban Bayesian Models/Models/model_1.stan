functions {
  //This describes the dynamics of the ODE
  real[] one_comp_lin_elim_abs(real t,
                               real[] y,
                               real[] theta,
                               real[] x_r,
                               int[] x_i) {
    real dydt[1];
    
    //PKPD params
    real k_a = theta[1]; // Dosing rate in 1/day
    real k = theta[2];// Elimination rate in 1/day
    
    //Put a prior on the volume of the absorptive compartment volume (i.e the volume of the gut)
    real V = theta[3];
    
    //Dose is known
    real D = x_r[1];
    real dose = 0;
    real elim = k * y[1];

    if (t > 0)
      dose = exp(- k_a * t) * D * k_a / V;

    dydt[1] = dose - elim;

    return dydt;
  }
}

data {
  real t0;    // Initial time in Hours post dose;
  real C0[1]; // Initial concentration at t0 in mg/L
  real D;   // Total dosage in mg

  //This is the number of unique sampling points (here, N_tt = 8)
  int<lower=1> N_tt;
  real times[N_tt];   // Measurement times in days

  //Number of observations (equal to N_tt x number of patients)
  int N_t;
  //An index for accessing times.  times[ix] are the sampling times.  So we observe C_hat at times[ix]
  int ix[N_t];
  
  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_t];

  
}

transformed data {
  real x_r[1] = {D};
  int x_i[0];
}

parameters {
  real<lower=0> k_a; // Dosing rate in 1/day
  real<lower=0> k;   // Elimination rate in 1/day
  real<lower=0> sigma;
  real<lower=0> V;
}

transformed parameters {

  real C[N_tt, 1];
  {
    real theta[3] = {k_a, k,V};
    
    //C is now a vector of solutions to the ODE with params listed.
    //C(t) = < C(t = times[1]), C(t = times[2]), ... >
    C = integrate_ode_bdf(one_comp_lin_elim_abs, C0, t0, times, theta, x_r, x_i);
  }
}

model {
  // Priors
  k_a ~ cauchy(0, 1);
  k ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);
  
  //Following based off of some googling
  V ~ normal(.1, 0.075);

  // Likelihood
  // for (n in 1:N_t)
  //   C_hat[n] ~ lognormal(log(C[ix[n], 1]), sigma);
  
  C_hat ~ lognormal(log(C[ix, 1]), sigma);
}

generated quantities {
  real C_ppc[N_tt];
  for (n in 1:N_tt)
    C_ppc[n] = lognormal_rng(log(C[n, 1]), sigma);

  
}

