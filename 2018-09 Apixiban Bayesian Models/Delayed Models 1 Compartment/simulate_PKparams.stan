
data {
  int<lower=1> N_covariates;
  int<lower=1> N_patients;
  matrix[N_patients,N_covariates] X;
}
model{
  
}
generated quantities{
  
  vector[N_covariates] BETA;
  real<lower=0> SIGMA;
  vector[N_patients] z;
  vector[N_patients] PKP;
  
  
  for (i in 1:N_covariates) BETA[i] = normal_rng(0,1);
  for (i in 1:N_patients) z[i] = normal_rng(0,1);
  

  SIGMA = fabs(normal_rng(0,0.1));

  
  PKP = exp(X*BETA + z*SIGMA);

}
