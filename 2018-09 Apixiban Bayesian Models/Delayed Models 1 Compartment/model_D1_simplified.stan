// Model with delay in absoprtion
// Time delay suggested by Rommel
// Log-normal likelihood
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
                     
  
  
  int<lower=1>N_covariates;
  int<lower=1> N_patients;
  int<lower=1> N;
  int<lower=1> N_v_c;
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
                                     //BMI
                                     
  matrix[N_patients,N_v_c] X_v;
}
parameters {
  
  //Volume
  vector[N_v_c] BETA_V;
  real<lower=0> SIGMA_V;
  real z_V[N_patients];
  
  //Elimination Rate
  real k;
  
  //Absorption
  vector[N_covariates] BETA_KA;
  real<lower=0> SIGMA_KA;
  real z_ka[N_patients];
  
  real<lower=0>sigma;
  
  //ADDED: Delay effects
  //See "Modeling of delays in PKPD: classical approaches and 
  //a tutorial for delay differential equations" -- Koch G., et. al
  //From a ODE model, we can cook up a delay.
  //Delay, in this context, is an estimate of the mean transit time.
  //Hence, times are interpreted as times after absorption, 
  //rather than times after administration
  
  real<lower=0,upper=1> phi; //Use to measure population delay.
  real<lower=10> lambda;
  real<lower=0, upper=1> delay_raw[N_patients]; //Each patient has their own delay
}

transformed parameters {

  real C[N_patients, N_t];
  real<lower=0> k_a[N_patients];
  real<lower=0> V[N_patients];
  
  vector[N_patients] MU_KA;
  vector[N_patients] MU_V;
  
  real delay[N_patients];
  
  MU_KA = X*BETA_KA;
  
  MU_V = X_v*BETA_V;
  
  for (n in 1:N_patients) {
    
    //PK parameters have lognormal distribution
    k_a[n] = exp(MU_KA[n] + z_ka[n]*SIGMA_KA);
    V[n] = exp(MU_V[n] + z_V[n]*SIGMA_V);
    delay[n] = 0.5*delay_raw[n];
    
    for (t in 1:N_t)
      //Our pharmacologist says that there can be a delay 
      //in the absorption of the drug
      //So although the ODE model assumes the drug is instantaneously absorbed, 
      //this is def not the case.
      //We posit that the drus is absorbed a little later than we think, 
      //which translates to an error in our time measurements
      //Alternatively, that the times are times after absorption, 
      //not times after administration.
      
      C[n,t] = PK_profile(times[t] - delay[n], D, V[n], k_a[n], k);

  }
}

model {
  
  // Priors
  //Coefficients
  BETA_V  ~ normal(0,1);
  BETA_KA ~ normal(0,1);
  k ~ normal(0,1);
  //Noise
  SIGMA_V ~ normal(0,3);
  SIGMA_KA ~normal(0,3);
  sigma ~ normal(0,3);  
  //Random Effects
  z_V ~ normal(0,1);
  z_ka ~ normal(0,1);

  
  //ADDED: Delay
  //These priors make more sense to me and
  //result in posteriors which are not much different.
  phi ~ beta(2.5, 2.5);
  lambda ~ pareto(10, 1.5);
  delay_raw ~ beta(lambda * phi, lambda * (1 - phi));

  // Likelihood
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_hat[n, t] ~ lognormal(log(C[n, t]), sigma);
}

generated quantities {
  real C_ppc[N_patients, N_t];
  vector[N] log_lik;
  for (n in 1:N_patients){
    for (t in 1:N_t){
      C_ppc[n, t] = lognormal_rng(log(C[n, t]), sigma );
      log_lik[N_t*(n-1) + t] = lognormal_lpdf( C_hat[n, t] | log(C[n,t]), sigma);
    }
  }
      
  
          
  
}
