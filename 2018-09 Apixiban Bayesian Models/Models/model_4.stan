functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data{
  real t0;    // Initial time in days;
  real C0[1]; // Initial concentration at t0 in mg/L

  real D;   // Total dosage in mg.  2.5 mg of Apixaban

  int<lower=1> N_t;
  real times[N_t];   // Measurement times in days
  int sex[N_t];
  int group[N_t];
  
  //Use these for generating samples at time unique times
  int<lower=1> N_ut;
  real utimes[N_ut];
  
  // Measured concentrations in effect compartment in mg/L
  real C_hat[N_t];
}
parameters {

  //V effects
  real b_V_int;
  real b_V_sex;
  real b_V_group;
  real b_V_intxn;
  
  //k effects
  real b_k_int;
  real b_k_sex;
  real b_k_group;
  real b_k_intxn;
  
  //k_a effects
  real b_k_a_int;
  real b_k_a_sex;
  real b_k_a_group;
  real b_k_a_intxn;
  
  real<lower=0> sigma;
}
transformed parameters {

  // Store analytic concentration
  real C[N_t];
  // PKPD params. 
  real V[N_t];
  real k_a[N_t];
  real k[N_t];
  
  //Will use in PPC
  real C_pred[2,N_ut,2]; //[sex,time,group]
  real V_ppc;
  real k_ppc;
  real k_a_ppc;
  

  // For the likelihood
  
  // REMOVE THE ZEROS FOR INTERACTION
  for (n in 1:N_t){
    V[n] = exp(b_V_int + 
                sex[n]*b_V_sex + 
                group[n]*b_V_group + 
                sex[n]*group[n]*b_V_intxn );
                
    k_a[n] = exp(b_k_a_int + 
                 sex[n]*b_k_a_sex + 
                 group[n]*b_k_a_group + 
                 sex[n]*group[n]*b_k_a_intxn);
                 
    k[n] = exp(b_k_int + 
                sex[n]*b_k_sex + 
                group[n]*b_k_group + 
                sex[n]*group[n]*b_k_intxn);
  }
  for (n in 1:N_t)
    C[n] = C_anal(times[n], D, V[n], k_a[n], k[n]);
    
  //For plots
  for (i in 1:2){
    for (j in 1:2){
    V_ppc = exp(b_V_int + (i-1)*b_V_sex + (j-1)*b_V_group + b_V_intxn*(i-1)*(j-1) );
    k_a_ppc = exp(b_k_a_int + (i-1)*b_k_a_sex + (j-1)*b_k_a_group + b_k_a_intxn*(i-1)*(j-1));
    k_ppc = exp(b_k_int + (i-1)*b_k_sex + (j-1)*b_k_group + b_k_intxn*(i-1)*(j-1));
      for (n in 1:N_ut)    
        C_pred[i,n,j] = C_anal(utimes[n], D, V_ppc, k_a_ppc, k_ppc);
    }
  }
}

model {
  b_V_int ~ normal(0,1);
  b_V_sex ~ normal(0,1);
  b_V_group ~ normal(0,1);
  b_V_intxn ~ normal(0,1);
  
  //k effects
  b_k_int ~ normal(0,1);
  b_k_sex ~ normal(0,1);
  b_k_group ~ normal(0,1);
  b_k_intxn ~ normal(0,1);
  
  //k_a effects
  b_k_a_int ~ normal(0,1);
  b_k_a_sex ~ normal(0,1);
  b_k_a_group ~ normal(0,1);
  b_k_a_intxn ~ normal(0,1);
  
  sigma ~ cauchy(0,1);
  // Likelihood
  for (n in 1:N_t)
    C_hat[n] ~ lognormal(log(C[n]) - sigma/2, sigma);
    // C_hat[n] ~ normal(C[n], sigma);
}

generated quantities {
  real C_ppc[2,N_ut,2];
  for (i in 1:2){
    for (j in 1:2){
      for (n in 1:N_ut)
        C_ppc[i,n,j] = lognormal_rng(log(C_pred[i,n,j]) - sigma/2, sigma);
        // C_ppc[i,n,j] = normal_rng(C_pred[i,n,j], sigma);
    }
  }
}

