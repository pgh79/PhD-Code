// Model with delay in absoprtion
// Time delay suggested by Rommel
// Log-normal likelihood
// Author: Demetri Pananos


//------------------------Introduction------------------------
//The goal of the analysis is as follows:  I have 36 patients, each have 8 
// time points at which their blood is sampled and the concentration of a drug
// is measured.  I'm intereted in predictive accuracy as measured by MAE or RMSE

//In this model, I am leaving some data out to do a cross validation.  This
//model will condition on some of the data from a patient and will test on the 
//remaining data from that patient.
//So let's say I fit my model on 35 patients + 1 sample from the patient I've
//held out, and then the model rpedicts on the remianing 7 points from 
//the left out patient.


functions {
  
  // PK Function.  Solution to differential equation
  // y' = ka*(D/V)exp(-ka*t) - k*y, y(0) = 0
  //This function is what we think the drug concentrtions should look like
  //over time.  The model estimates the parameters V,k_a, and k.
  real PK_profile(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k)) * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data{
  //The following are what I provide the program
  //Total Dose in mg
  real D;
  
  //Num obs for training data
  int<lower=1> N_train;
  
  //Num distinct patients
  int<lower=1> N_patients_train;
  
  //Patient IDs
  //used to index parameters in the model.
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
  
  //Number of test points
  int<lower=0> N_test;
  //test times
  real<lower=0> times_test[N_test];
  //test covariates
  row_vector[p] X_test;
  //the ID of the patient left out.
  int<lower=0> ID;
  
}
parameters{
  
  //These are model parameters we will estimate
  // parameter: V Volume
  vector[3] BETA_V; //Fixed effets
  real<lower=0> SIGMA_V;
  vector[N_patients_train] z_V; //Random Effects
  
  // parameter: k Drug Elimination
  vector[p] BETA_k; //Fixed effects
  real<lower=0> SIGMA_k;
  vector[N_patients_train] z_k; //Random effects
  
  // parameter: ka Drug Absorption
  vector[4] BETA_ka; //Fixed effects
  real<lower=0> SIGMA_ka;
  vector[N_patients_train] z_ka; //random effects
  
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
  
  //Here is where we compute the concentrations for each patient
  
  //Training concentrations
  real C[N_train];
  
  //Parameters for each patient.  Each patient has their own V,k,k_
  vector[N_patients_train] k_a;
  vector[N_patients_train] k;
  vector[N_patients_train] V;
  
  //Parameter means
  vector[N_patients_train] MU_KA;
  vector[N_patients_train] MU_K;
  vector[N_patients_train] MU_V;
  
  MU_KA = X_train[,{1,3,4,5}]*BETA_ka;  //Mean of the ka distribution
  MU_K = X_train*BETA_k; //mean of the k distribution
  MU_V = X_train[,{1,3,5}]*BETA_V;  //Mean of the V distribution.
                                    //Only uses Baseline, Ismale, and Weight

  //Posit that the parameters are lognormal with mean MU
  k_a = exp(MU_KA + z_ka*SIGMA_ka); 
  k = exp(MU_K + z_k*SIGMA_k);
  V = exp(MU_V + z_V*SIGMA_V);
  
  
  //compute the concentrations at each time.
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
  
  //Here is where we do the testing
  vector[N_test] C_pred;
  //Below, I predict the concentration for the held out data
  //From the data I used to train on, I have an estimate of:
  // -Patient's delay
  // -Patient's V
  // -Patient's k
  // -Patient's k_a
  // Use these estimates in the prediction
  // Note that X_test*BETA_i is the mean of the log normal for somoene with
  // covariates X_test, and z_i[ID]SIGMA_i is the random effect for that
  // particular person
  
  real delay_test;
  real V_test;
  real ka_test;
  real k_test;
  
  delay_test = 0.5*delay_raw[ID];
  V_test = exp(X_test[{1,3,5}]*BETA_V + z_V[ID]*SIGMA_V);
  ka_test = exp(X_test[{1,3,4,5}]*BETA_ka + z_ka[ID]*SIGMA_ka);
  k_test = exp(X_test*BETA_k + z_k[ID]*SIGMA_k);
  
  for (i in 1:N_test)
    C_pred[i] = PK_profile(times_test[i] -delay_test,
                    D,
                    V_test,
                    ka_test,
                    k_test
                    ); 
                    
}
