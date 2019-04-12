library(drake)
library(rstan)
source('01_simulate_covariates.R')
source('02_simulate_pk_params.R')
rstan_options(auto_write= T)

plan<- drake_plan(
  n_patients = 5000,
  p_continuous = 2,
  p_binary = 2,
  theta = 0.5,
  X = simulate_covariates(n_patients, p_continuous, p_binary, theta),
  simulated_params = simulate_pk_params(n_patients,p_continuous , p_binary , X)
)

clean()
make(plan, lock_envir = F)
##Recode sampling for pk params so that excretion is sampled via the non dimensionalization
p = readd(simulated_params)
r = p$pk_params[1,,2]/p$pk_params[1,,3]
hist(r) 




