library(tidyverse)
library(rstan)
source('stan_utilities.R')
source('LOPO.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

t = c(0.5,1,2,4,6,8,10,12)
indata = rstan::read_rdump('Model Data/core_model_raw_data')

fit.model<- function(f){
  
  
  
  fit = rstan::stan('Models/LOPO_0_condition.stan', 
                    data = f,
                    chains = 16,
                    control = list(max_treedepth = 13,adapt_delta = 0.99),)
  
  # check_all_diagnostics(fit)
  
  p = rstan::extract(fit)
  params = lapply(p, mean)
  ypred = apply(p$C_pred,2,mean)
  
  return(list(ypred = ypred, params = params))
  
}

tibble(patients = 1:36, 
       times = list(t), 
)  %>% 
  mutate(LOPOData = map2(patients,times,~LeaveOnePatientOut(indata,.x,.y)),
         ytest = map(LOPOData,'C_hat_test'),
         results = map(LOPOData,fit.model),
         ypred = map(results, 'ypred'),
         params = map(results,'params'),
         err = map2(ytest,ypred,~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_none.RDS')

