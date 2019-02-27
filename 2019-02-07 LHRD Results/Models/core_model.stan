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
  int<lower=1> N;
  
  //Num distinct patients
  int<lower=1> N_patients;
  
  //Patient IDs
  int<lower=1> patient_ID[N];
  
  //Num covariates
  int<lower=1> p;
  
  //Design Matrix
  matrix[N_patients,p] X;
  
  //Sample times
  real<lower=0> times[N];
  
  //Scaled concentration in mg/L
  real<lower=0> C_hat[N];
}
parameters{
  // parameter: V Volume
  vector[3] BETA_V;
  real<lower=0> SIGMA_V;
  vector[N_patients] z_V;
  
  // parameter: k Elimination
  vector[p] BETA_k;
  real<lower=0> SIGMA_k;
  vector[N_patients] z_k;
  
  // parameter: ka Absorption
  vector[5] BETA_ka;
  real<lower=0> SIGMA_ka;
  vector[N_patients] z_ka;
  
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
  real<lower=0,upper=1> delay_raw[N_patients]; //Each patient has their own delay
}
transformed parameters{
  //Predicted concentrations
  real C[N];
  
  //Parameters for each patient
  vector[N_patients] k_a;
  vector[N_patients] k;
  vector[N_patients] V;
  
  //Parameter means
  vector[N_patients] MU_KA;
  vector[N_patients] MU_K;
  vector[N_patients] MU_V;
  
  MU_KA = X[,{1,3,4,5,6}]*BETA_ka; //Only Baseline, IsMale, Age,Weight, Creatinine
  MU_K = X*BETA_k;
  MU_V = X[,{1,3,5}]*BETA_V;  //Only Baseline, Ismale, and Weight

  k_a = exp(MU_KA + z_ka*SIGMA_ka);
  k = exp(MU_K + z_k*SIGMA_k);
  V = exp(MU_V + z_V*SIGMA_V);
  
  for (i in 1:N){
      //Our pharmacologist says that there can be a delay 
      //in the absorption of the drug
      //So although the ODE model assumes the drug is instantaneously absorbed, 
      //this is def not the case.
      //We posit that the drus is absorbed a little later than we think, 
      //which translates to an error in our time measurements
      //Alternatively, that the times are times after absorption, 
      //not times after administration.
    C[i] = PK_profile(times[i] - 0.5*delay_raw[patient_ID[i]],
                      D,
                      V[patient_ID[i]],
                      k_a[patient_ID[i]],
                      k[patient_ID[i]]
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

  
  //Delay
  //These priors make more sense to me and
  //result in posteriors which are not much different.
  
  //Phi is the mean of the betra
  //lambda sum of successes and losses
  //See stan manual See section 20.2 : Reparameteriztions
  phi ~ beta(2.5, 2.5);
  lambda ~ pareto(10, 1.5);
  delay_raw ~ beta(lambda * phi, lambda * (1 - phi));

  // Likelihood
  C_hat ~ lognormal(log(C), sigma);
  
}
generated quantities{
  real C_ppc[N];
  vector[N] log_lik;
  
  C_ppc = lognormal_rng(log(C), sigma);
  
  for (i in 1:N)
    log_lik[i] = lognormal_lpdf(C_hat[i] | log(C[i]), sigma);
}
