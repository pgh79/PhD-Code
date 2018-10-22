functions {
  real C_anal(real t, real D, real V, real ka, real k, real k12, real k21) {

      real t1;
      real t4;
      real t6;
      real t9;
      real t11;
      real t12;
      real t15;
      real t17;
      real t19;
      real t20;
      real t25;
      real t26;
      real t28;
      real t31;
      real t34;
      real t40; 
      real t44;
      real t46;
      real t51;
      real t66;

      t1 = k*k;
      t4 = k21*k;
      t6 = k12*k12;
      t9 = k21*k21;
      t11 = sqrt(2.0*k*k12+2.0*k12*k21+t1-2.0*t4+t6+t9);
      t12 = -k-k12-k21+t11;
      t15 = exp(t12*t/2.0);
      t17 = 1/t11;
      t19 = k12*ka;
      t20 = 2.0*ka;
      t25 = ka*ka;
      t26 = -k*ka-k21*ka-t19+t25+t4;
      t28 = t19/(t20-k-k12-k21+t11)*t26;
      t31 = -k-k12-k21-t11;
      t34 = exp(t31*t/2.0);
      t40 = t19/(-t11-k-k12-k21+t20)*t26;
      t44 = exp(-ka*t);
      t46 = D*k12;
      t51 = D*t17;
      t66 = (k21*(2.0*t15*D*t17*t28-2.0*t34*D*t17*t40+t44*ka*t46)+t12*t15*t51*t28-t31*t34*t51*t40-t25*t44*t46)/t26/V/k12;

      return(t66);
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

parameters {
  vector[4] beta_v;
  real<lower=0> sigma_v;
  real z_v[N_patients];
  
  vector[4] beta_k;
  real<lower=0> sigma_k;
  real z_k[N_patients];
  
  vector[4] beta_ka;
  real<lower=0> sigma_ka;
  real z_ka[N_patients];
  
  vector[4] beta_k12;
  real<lower=0> sigma_k12;
  real z_k12[N_patients];
  
  vector[4] beta_k21;
  real<lower=0> sigma_k21;
  real z_k21[N_patients];
  
  real<lower=0> sigma;
  real<lower=0, upper=1> alpha;
}

transformed parameters {
  real C[N_patients, N_t];
  vector[N_patients] mu_V;
  vector[N_patients] mu_k;
  vector[N_patients] mu_ka;
  vector[N_patients] mu_k12;
  vector[N_patients] mu_k21;
  
  real<lower=0> V;
  real<lower=0> k;
  real<lower=0> ka;
  real<lower=0> k12;
  real<lower=0> k21;
  
  mu_V = X*beta_v;
  mu_k = X*beta_k;
  mu_ka = X*beta_ka;
  mu_k12 = X*beta_k12;
  mu_k21 = X*beta_k21;
  {
    for (p in 1:N_patients){
      V = exp(mu_V[p] + z_v[p]*sigma_v);
      k = exp(mu_k[p] + z_k[p]*sigma_k);
      ka = exp(mu_ka[p] + z_ka[p]*sigma_ka);
      k12 = exp(mu_k12[p] + z_k12[p]*sigma_k12);
      k21 = exp(mu_k12[p] + z_k21[p]*sigma_k21);
    for (n in 1:N_t)
    C[p,n] = C_anal(times[n],D, V, ka, k, k12, k21);
  }
}
}

model {
  // Priors
  
  beta_ka ~ normal(0, 1);
  z_ka ~ normal(0,1);
  sigma_ka ~ cauchy(0,1);

  beta_k ~ normal(0, 1);
  z_k ~ normal(0,1);
  sigma_k ~ cauchy(0,1);
  
  beta_v ~ normal(0,1);
  z_v ~ normal(0,1);
  sigma_v ~ cauchy(0,1);

  beta_k12 ~ normal(0,1);
  z_k12 ~ normal(0,1);
  sigma_k12 ~ cauchy(0,1);
  
  beta_k21 ~ normal(0,1);
  z_k21 ~ normal(0,1);
  sigma_k21 ~ normal(0,1);
  
  sigma ~ cauchy(0, 1);
  alpha ~ uniform(0, 1);

  // Likelihood
  for (p in 1:N_patients){
  for (n in 1:N_t)
    C_hat[p,n] ~ normal(C[p,n], sigma*pow(C[p,n]/0.4,2*alpha));
    }
}

generated quantities {
  real C_ppc[N_patients,N_t];
  for (p in 1:N_patients){
   for (n in 1:N_t)
    C_ppc[p,n] = normal_rng(C[p,n], sigma*pow(C[p,n]/0.4,2*alpha));
  }
}
