library(rstan)
library(glue)
source('stan_utilities.R')
source('LOPO.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

results = list()
input_data = rstan::read_rdump('Model Data/core_model_raw_data')
N = input_data$N_patients

divs = 0
t = c(0.5,1,2,4,6,10)
  for (i in 1:N){
    
    phrase = glue('Working on Leave Patient {i} Out')
    print(phrase)
    iter = list()
    lopo_data = LeaveOnePatientOut(input_data,i, t)
    lopo_data$ID = i
    
    
    fit = stan('Models/LOPO_1_condition.stan',
               data = lopo_data,
               chains = 12,
               seed = 10090908,
               refresh = 0,
               control = list(max_treedepth = 20,adapt_delta = 0.999))
    
    
    check_all_diagnostics(fit)
    p = rstan::extract(fit)
    
    iter$RMSE = mean(p$RMSE)
    iter$MAE = mean(p$MAE)
    iter$relative_error = apply(p$relative_error,2, mean)
    
    results[[i]] = iter
  }
  
  
  saveRDS(results, glue('LOPO Results/results_4_condition_test') )


