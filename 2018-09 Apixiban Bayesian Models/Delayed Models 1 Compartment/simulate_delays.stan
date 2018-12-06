
data {
  real A;
  real B;
  real ymin;
  real Z;
  
}
model{
 
}
generated quantities{
  real<lower=0, upper=1> phi;
  real<lower=0> lambda;
  
  real<lower=0> alpha;
  real<lower=0> beta;
  
  real delay;
  
  phi = beta_rng(A,B);
  lambda = pareto_rng(ymin,Z);
  alpha = lambda * phi;
  beta = lambda * (1 - phi);
  
  delay = beta_rng(alpha,beta);  
  
  
}
