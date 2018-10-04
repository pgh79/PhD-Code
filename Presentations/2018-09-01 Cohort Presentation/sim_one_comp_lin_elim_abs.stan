functions {
  real[] one_comp_lin_elim_abs(real t,
                               real[] y,
                               real[] theta,
                               real[] x_r,
                               int[] x_i) {
    real dydt[1];
    real k_a = theta[1]; // Dosing rate in 1/day
    real k = theta[2];   // Elimination rate in 1/day
    real D = x_r[1];
    real V = x_r[2];
    real dose = 0;
    real elim = k * y[1];

    if (t > 0)
      dose = exp(- k_a * t) * D * k_a / V;

    dydt[1] = dose - elim;

    return dydt;
  }
}

transformed data {
  int N_t = 20;
  real times[N_t];
  real t0 = 0;
  real C0[1] = {0.0};
  real theta[2] = {0.75, 0.5}; // Dosing rate in 1/day, Elimination rate in 1/day
  real sigma = 0.1;
  real x_r[2] = {30.0, 2.0}; // Total dosage in mg,  Comparment volume in L
  int x_i[0];

  for (n in 1:N_t)
    times[n] =  0.5*n;
}

model {}

generated quantities {
  real t_init = t0;
  real C_init[1] = {C0[1]};

  real D = x_r[1];
  real V = x_r[2];
  real ts[N_t];

  real C[N_t, 1];
  real C_hat[N_t];

  for (n in 1:N_t)
    ts[n] = times[n];

  C = integrate_ode_bdf(one_comp_lin_elim_abs, C0, t0, times, theta, x_r, x_i);

  for (n in 1:N_t)
    C_hat[n] = lognormal_rng(log(C[n, 1]), sigma);
}
