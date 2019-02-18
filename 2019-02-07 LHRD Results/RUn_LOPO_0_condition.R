library(rstan)
library(glue)
source('stan_utilities.R')
source('LOPO_tools.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

results = list()
input_data = rstan::read_rdump('Model Data/core_model')
N = input_data$N_patients
for (i in 1:N){
  
  phrase = glue('Working on Leave Patient {i} Out')
  print(phrase)
  iter = list()
  lopo_data = LeaveOnePatientOut(input_data = input_data, i = i)
  
  
  fit = stan('Models/LOPO_0_condition.stan',
             data = lopo_data,
             chains = 12,
             seed = 10090908,
             refresh = 0,
             control = list(max_treedepth = 13,adapt_delta = 0.8))
  
  
  check_all_diagnostics(fit)
  
  p = rstan::extract(fit)
  
  iter$RMSE = mean(p$RMSE)
  iter$MAE = mean(p$MAE)
  iter$relative_error = apply(p$relative_error,2, mean)
  
  results[[i]] = iter
}


saveRDS(results, 'LOPO Results/results_0_condition')