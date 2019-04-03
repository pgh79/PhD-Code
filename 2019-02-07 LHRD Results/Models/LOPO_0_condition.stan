// Model with delay in absoprtion
// Time delay suggested by Rommel
// Log-normal likelihood
// Author: Demetri Pananos

//Model should fit on N-1 patients
//Then predict on held out patient
//without conditioning on any held out patient data

functions {
  
  // PK Function.  Solution to differential equation
  // y' = ka*(D/V)exp(-ka*t) - k*y, y(0) = 0
  real PK_profile(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data{
  
  //Total Dose in mg
  real D;
  
  //Num obs
  int<lower=1> N_train;
  
  //Num distinct patients
  int<lower=1> N_patients_train;
  
  //Patient IDs
  int<lower=1> patient_ID_train[N_train];
  
  //Num covariates
  int<lower=1> p;
  
  //Design Matrix
  matrix[N_patients_train,p] X_train;
  
  //Sample times
  real<lower=0> times_train[N_train];
  
  //Scaled concentration in mg/L
  real<lower=0> C_hat_train[N_train];
  
  
  //Held out data
  int<lower=0> N_test;
  vector[N_test] C_hat_test;
  real<lower=0> times_test[N_test];
  row_vector[p] X_test;
  
}
parameters{
  // parameter: V Volume
  vector[3] BETA_V;
  real<lower=0> SIGMA_V;
  vector[N_patients_train] z_V;
  
  // parameter: k Elimination
  vector[p] BETA_k;
  real<lower=0> SIGMA_k;
  vector[N_patients_train] z_k;
  
  // parameter: ka Absorption
  vector[4] BETA_ka;
  real<lower=0> SIGMA_ka;
  vector[N_patients_train] z_ka;
  
  //parameter: noise in likelihood
  real<lower=0> sigma;
  
  //Delay effects
  //See "Modeling of delays in PKPD: classical approaches and 
  //a tutorial for delay differential equations" -- Koch G., et. al
  //From a ODE model, we can cook up a delay.
  //Delay, in this context, is an estimate of the mean transit time.
  //Hence, times are interpreted as times after absorption, 
  //rather than times after administration
  
  real<lower=0,upper=1> phi; //Use to measure population delay.
  real<lower=10> lambda;
  real<lower=0,upper=1> delay_raw[N_patients_train]; //Each patient has their own delay
}
transformed parameters{
  //Predicted concentrations
  real C[N_train];
  
  //Parameters for each patient
  vector[N_patients_train] k_a;
  vector[N_patients_train] k;
  vector[N_patients_train] V;
  
  //Parameter means
  vector[N_patients_train] MU_KA;
  vector[N_patients_train] MU_K;
  vector[N_patients_train] MU_V;
  
  MU_KA = X_train[,{1,3,4,5}]*BETA_ka;
  MU_K = X_train*BETA_k;
  MU_V = X_train[,{1,3,5}]*BETA_V;  //Only Baseline, Ismale, and Weight

  k_a = exp(MU_KA + z_ka*SIGMA_ka);
  k = exp(MU_K + z_k*SIGMA_k);
  V = exp(MU_V + z_V*SIGMA_V);
  
  for (i in 1:N_train){
    C[i] = PK_profile(times_train[i] - 0.5*delay_raw[patient_ID_train[i]],
                      D,
                      V[patient_ID_train[i]],
                      k_a[patient_ID_train[i]],
                      k[patient_ID_train[i]]
                      );  
  }
  
}
model{
  
  // Priors
  //Coefficients
  BETA_V  ~ normal(0,0.5);
  BETA_V[1] ~ normal(log(5.5),0.1);
  BETA_k  ~ normal(0,0.5);
  BETA_ka ~ normal(0,0.5);
  
  //Noise
  SIGMA_V ~ cauchy(0,2);
  SIGMA_k ~ cauchy(0,2);
  SIGMA_ka ~cauchy(0,2);
  sigma ~ cauchy(0,2);
  
  //Random Effects
  z_V ~ normal(0,1);
  z_k ~ normal(0,1);
  z_ka ~ normal(0,1);

  
  //ADDED: Delay
  //These priors make more sense to me and
  //result in posteriors which are not much different.
  phi ~ beta(2.5, 2.5);
  lambda ~ pareto(10, 1.5);
  delay_raw ~ beta(lambda * phi, lambda * (1 - phi));

  // Likelihood
  C_hat_train ~ lognormal(log(C), sigma);
  
}
generated quantities{
  real C_ppc[N_train];
  vector[N_test] C_pred;
  
  
  C_ppc = lognormal_rng(log(C), sigma);
  
  for (i in 1:N_test)
    C_pred[i] = PK_profile(times_test[i] - phi,
                    D,
                    exp(X_test[{1,3,5}]*BETA_V),
                    exp(X_test[{1,3,4,5}]*BETA_ka),
                    exp(X_test*BETA_k)
                    ); 
}
