functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
transformed data {
  int N_t = 100;
  real delta_t = (12.5-0.5)/100;
  real times[N_t];
  real t0 = 0;
  real C0[1] = {0.0};
  real D = 2.5;
  
  for (i in 1:N_t)
    times[i] = 0.5 + i*delta_t;
}
model {}
generated quantities {
  real t_init = t0;
  real C_init[1] = {C0[1]};
  
  
  real<lower=0> V = sqrt(cauchy_rng(0,1)^2) ;
  real<lower=0> k_a = sqrt(cauchy_rng(0,1)^2) ;
  real<lower=0> k = sqrt(cauchy_rng(0,1)^2) ;

  real C[N_t];
  
  for (t in 1:N_t)
    C[t] = C_anal(times[t],D,V,k_a,k);

}
