# Dan thinks we should measure error a different way.
# This script does that

library(tidyverse)
library(rstan)
source('stan_utilities.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

LOPO<- function(data,test_patient,train_times){

  env = environment()

  list2env(indata, envir = env)
  
  condition_on = ( (patient_ID==test_patient)&(times %in% train_times) )|(patient_ID!=test_patient)
  for_testing = patient_ID==test_patient
  
  
  D = 2.5
  patient_ID_train = patient_ID[condition_on]
  times_train = times[condition_on]
  C_hat_train = C_hat[condition_on]
  N_train = sum(condition_on)
  N_patients_train = length(unique(patient_ID_train))
  p<-dim(X)[2]
  X_train = X
  
  ID = test_patient
  times_test = times[for_testing]
  C_hat_test = C_hat[for_testing]
  X_test = X[ID,]
  N_test = sum(for_testing)
  
  training_data = list(D = D, 
                       patient_ID_train = patient_ID_train,
                       times_train = times_train,
                       C_hat_train = C_hat_train,
                       N_train = N_train,
                       N_patients_train = N_patients_train,
                       p = p,
                       ID = ID,
                       times_test = times_test,
                       X_test = X_test,
                       X_train = X_train,
                       N_test = N_test
                       )
  
  testing_data = C_hat_test
  
  return(list(training_data = training_data, testing_data = testing_data))
}



fit.model<- function(f){
      
      
      
      fit = rstan::stan('Models/LOPO_1_condition_new.stan', 
                        data = f,
                        chains = 4,
                        control = list(max_treedepth = 13,adapt_delta = 0.99),)
      
      # check_all_diagnostics(fit)
      
      p = rstan::extract(fit)
      ypred = apply(p$C_pred,2,mean)
      
      return(ypred)
      
    }


indata = rstan::read_rdump('Model Data/core_model_scaled_data')


#Condition on first----
tibble(paitents = 1:36, 
       times = list(c(0.5)), 
       )  %>% 
  mutate(LOPOData = map2(paitents,times,~LOPO(indata,.x,.y)),
         training_data = map(LOPOData,'training_data'),
         ypred = map(training_data, fit.model),
         ytest = map(LOPOData,'testing_data')
         ) %>% 
  mutate(err = map2(ytest, ypred, ~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_first.RDS')

  

#Condition on first two----
tibble(paitents = 1:36, 
       times = list(c(0.5,1)), 
)  %>% 
  mutate(LOPOData = map2(paitents,times,~LOPO(indata,.x,.y)),
         training_data = map(LOPOData,'training_data'),
         ypred = map(training_data, fit.model),
         ytest = map(LOPOData,'testing_data')
  ) %>% 
  mutate(err = map2(ytest, ypred, ~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_first_two.RDS')


#Condition on last----
tibble(paitents = 1:36, 
       times = list(c(12)), 
)  %>% 
  mutate(LOPOData = map2(paitents,times,~LOPO(indata,.x,.y)),
         training_data = map(LOPOData,'training_data'),
         ypred = map(training_data, fit.model),
         ytest = map(LOPOData,'testing_data')
  ) %>% 
  mutate(err = map2(ytest, ypred, ~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_last.RDS')



#Condition on last two----
tibble(paitents = 1:36, 
       times = list(c(10,12)), 
)  %>% 
  mutate(LOPOData = map2(paitents,times,~LOPO(indata,.x,.y)),
         training_data = map(LOPOData,'training_data'),
         ypred = map(training_data, fit.model),
         ytest = map(LOPOData,'testing_data')
  ) %>% 
  mutate(err = map2(ytest, ypred, ~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_last_two.RDS')



#Condition on first and last----
tibble(paitents = 1:36, 
       times = list(c(0.5,12)), 
)  %>% 
  mutate(LOPOData = map2(paitents,times,~LOPO(indata,.x,.y)),
         training_data = map(LOPOData,'training_data'),
         ypred = map(training_data, fit.model),
         ytest = map(LOPOData,'testing_data')
  ) %>% 
  mutate(err = map2(ytest, ypred, ~.x-.y)) %>% 
  saveRDS('LOPO Data/condition_on_last_two.RDS')

