## Dan wants me to check if the phenomenon we observe 
#(that sampling earlier is worse than sampling later)
#is from model mispecification or a real thing.
#This script fits my model,
#creates a new dataset where the target is generated from the posterior
#If the model is mispecified, this phenom will vanish
#If the pheom is real, we will observe it in generated data


library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('LOPO.R')
source('stan_utilities.R')

#Get the scaled data
input_data = rstan::read_rdump('Model Data/core_model_scaled_data')

#Fit a model
fit = stan('Models/core_model.stan',
           data = input_data,
           chains = 8,
           seed = 10090908,
           control = list(max_treedepth = 13,adapt_delta = 0.8),
           verbose = F)

#Extract the samples
p = rstan::extract(fit)


#REPLACE THE ORIGINAL CONCENTRATIONS WITH THE SAMPLED ONES
#NOW, DO LOPO.  DO WE STILL SEE THE SAME PHENOM?
raw_data = read_rdump('Model Data/core_model_raw_data')
raw_data$C_hat = p$C_ppc[1,]

#----Leave One Patient Out: No Samples----

t = c(0.5,1,2,4,6,8,10,12)
results = list()
N = raw_data$N_patients

for (i in 1:N){
  
  iter = list()
  lopo_data = LeaveOnePatientOut(raw_data,i, t)
  
  
  fit = stan('Models/LOPO_0_condition.stan',
             data = lopo_data,
             chains = 8,
             seed = 10090908,
             # refresh = 0,
             control = list(max_treedepth = 13,adapt_delta = 0.8))
  
  
  check_all_diagnostics(fit)
  
  p = rstan::extract(fit)
  
  iter$RMSE = mean(p$RMSE)
  iter$MAE = mean(p$MAE)
  iter$relative_error = apply(p$relative_error,2, mean)
  
  results[[i]] = iter
}


saveRDS(results, 'LOPO Results/results_0_condition_test_on_posterior_C')

print('Finished LOPO NO SAMPLE')

#----Leave One Patient Out: With Samples----

sample_times = list(
                    # c(1,2,4,6,8,10,12), 
                    # c(2,4,6,8,10,12),
                    c(0.5,1,2,4,6,8,10),
                    c(0.5,1,2,4,6,8)
                    )
file_names = list(
  # 'LOPO Results/results_posterior_C_conditioned_on_fist_point',
  # 'LOPO Results/results_posterior_C_conditioned_on_fist_two_points',
  'LOPO Results/results_posterior_C_conditioned_on_last_point',
  'LOPO Results/results_posterior_C_conditioned_on_last_two_points'
  
)

for( j in 1:length(sample_times)){
  results = list()
  for(i in 1:N){
    iter = list()
    lopo_data = LeaveOnePatientOut(raw_data,i, sample_times[[j]])
    lopo_data$ID = i
    
    
    fit = stan('Models/LOPO_1_condition.stan',
               data = lopo_data, 
               chains = 8,
               seed = 10090908,
               control = list(max_treedepth = 20,adapt_delta = 0.999))
    
    
    check_all_diagnostics(fit)
    p = rstan::extract(fit)
    
    iter$RMSE = mean(p$RMSE)
    iter$MAE = mean(p$MAE)
    iter$relative_error = apply(p$relative_error,2, mean)
    
    results[[i]] = iter
  }
  saveRDS(results, file_names[[j]] )
  
}

print('FINISHED SCRIPT')
