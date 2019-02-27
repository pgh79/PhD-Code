functions {
  
  // PK Function.  Solution to differential equation
  // y' = ka*(D/V)exp(-ka*t) - k*y, y(0) = 0
  real PK_profile(real t, real alpha) {
    return ((1 / (alpha - 1) * exp(t * (alpha - 1)) - 1 / (alpha - 1)) * exp(-alpha * t) );
  }
}
data{
  int<lower=0> train_N;
  int<lower=0> train_ID[train_N];
  real<lower=0> train_t[train_N];
  real<lower=0> train_y[train_N];
  int num_patients;
  
  int<lower=0> test_N;
  int<lower=0> test_ID[test_N];
  real<lower=0> test_t[test_N];
  real<lower=0> test_y[test_N];
}
parameters{
  real<lower=0, upper = 1> alpha[num_patients];
  real<lower=0> sigma;
}
transformed parameters{
  real C[train_N];
  
  for(i in 1:train_N){
    C[i] = PK_profile(train_t[i], 0.5*alpha[train_ID[i]]);
  }
}
model{
  
  alpha ~ beta(2,2);
  sigma ~ cauchy(0,2);
  train_y ~ lognormal(log(C), sigma);
  
  }
generated quantities{
  real C_pred[test_N];
  real err[test_N];
  
  for(i in 1:test_N){
    C_pred[i] = PK_profile(test_t[i], alpha[test_ID[i]]);
    err[i] = C_pred[i] - test_y[i];
  }
  
}
