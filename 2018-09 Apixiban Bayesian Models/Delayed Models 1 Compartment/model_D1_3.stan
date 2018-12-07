// Model with delay in absoprtion
// Time delay suggested by Rommel
// Log-normal likelihood
// Prior on the Dose
// Author: Demetri Pananos

functions {
  
  // PK Function.  Solution to differential equation
  // y' = ka*(D/V)exp(-ka*t) - k*y, y(0) = 0
  real PK_profile(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data {
  
  real D;   // Total dosage in mg

  int<lower=1> N_t; //Number of time points sampeld per patient
  real times[N_t];   //Measurement times. 
                     //Same for all patients.  Measured in hours post dose.

  real<lower=0> A; // Median observed concentration.
                   // To be used in the variance of likelihood
                   // See Gelman BDA3 chapter 19 or Jon Wakefield 1996
  
  
  int<lower=1>N_covariates;
  int<lower=1> N_patients;
  real C_hat[N_patients, N_t]; //Observed concentration.  
                              //  Originally in ng/ml.  
                              //Convert to mg/L bc Dose is in mg.
  
  
  
  matrix[N_patients,N_covariates] X; //Design matrix for patients
                                     // Covariates include 
                                     //Age, 
                                     //Weight, 
                                     //Creatinine, 
                                     //Sex (binary), 
                                     //Disease (binary)
}
parameters {
  
  //Volume
  vector[N_covariates] BETA_V;
  real<lower=0> SIGMA_V;
  real z_V[N_patients];
  
  //Elimination Rate
  vector[N_covariates] BETA_K;
  real<lower=0> SIGMA_K;
  real z_k[N_patients];
  
  //Absorption
  vector[N_covariates] BETA_KA;
  real<lower=0> SIGMA_KA;
  real z_ka[N_patients];
  
  real<lower=0>sigma;
  real<lower=0, upper=1> alpha;
  
  //ADDED: Delay effects
  //See "Modeling of delays in PKPD: classical approaches and 
  //a tutorial for delay differential equations" -- Koch G., et. al
  //From a ODE model, we can cook up a delay.
  //Delay, in this context, is an estimate of the mean transit time.
  //Hence, times are interpreted as times after absorption, 
  //rather than times after administration
  
  real<lower=0,upper=1> phi; //Use to measure population delay.
  real<lower=10> lambda;
  real<lower=0, upper=1> delay_raw[N_patients]; //Patients have their own delay
  
  //Bioavailability;
  real<lower=0,upper=1> phi_F; //Use to measure population delay.
  real<lower=10> lambda_F;
  real<lower=0, upper=1> F[N_patients];
}

transformed parameters {

  real C[N_patients, N_t];
  real<lower=0> k_a[N_patients];
  real<lower=0> k[N_patients];
  real<lower=0> V[N_patients];
  
  vector[N_patients] MU_KA;
  vector[N_patients] MU_K;
  vector[N_patients] MU_V;
  
  real delay[N_patients];
  real dose[N_patients];
  
  real alpha_F;
  real beta_F;
  
  alpha_F = lambda_F * phi_F;
  beta_F = lambda_F * (1 - phi_F);
  
  MU_KA = X*BETA_KA;
  MU_K = X*BETA_K;
  MU_V = X*BETA_V;
  
  for (n in 1:N_patients) {
    
    //PK parameters have lognormal distribution
    k_a[n] = exp(MU_KA[n] + z_ka[n]*SIGMA_KA);
    k[n] = exp(MU_K[n] + z_k[n]*SIGMA_K);
    V[n] = exp(MU_V[n] + z_V[n]*SIGMA_V);
    delay[n] = 0.5*delay_raw[n];
    dose[n] = D*F[n];
    
    for (t in 1:N_t)
      //Our pharmacologist says that there can be a delay 
      //in the absorption of the drug
      //So although the ODE model assumes the drug is instantaneously absorbed, 
      //this is def not the case.
      //We posit that the drus is absorbed a little later than we think, 
      //which translates to an error in our time measurements
      //Alternatively, that the times are times after absorption, 
      //not times after administration.
      
      C[n,t] = PK_profile(times[t] - delay[n], dose[n], V[n], k_a[n], k[n]);

  }
}

model {
  
  // Priors
  //Coefficients
  BETA_V  ~ normal(0,1);
  BETA_K  ~ normal(0,1);
  BETA_KA ~ normal(0,1);
  
  //Noise
  SIGMA_V ~ normal(0,1);
  SIGMA_K ~ normal(0,1);
  SIGMA_KA ~normal(0,1);
  
  //Random Effects
  z_V ~ normal(0,1);
  z_k ~ normal(0,1);
  z_ka ~ normal(0,1);

  sigma ~ normal(0,1);
  alpha ~ uniform(0,1); //See Gelman BDA, chapter 19.
  
  //These priors make more sense to me and
  //result in posteriors which are not much different.
  phi ~ beta(2.5, 2.5);
  lambda ~ pareto(10, 1.5);
  delay_raw ~ beta(lambda * phi, lambda * (1 - phi));
  
  //ADDED: Bioavailability
  phi_F ~ beta(2.5, 2.5);
  lambda_F ~ pareto(10, 1.5);
  F ~ beta(lambda_F * phi_F, lambda_F * (1 - phi_F));

  // Likelihood
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_hat[n, t] ~ lognormal(log(C[n, t]), sigma);
}

generated quantities {
  real C_ppc[N_patients, N_t];
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_ppc[n, t] = lognormal_rng(log(C[n, t]), sigma );
}
