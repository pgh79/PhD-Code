

simulate_pk_params<-function(n_patients,p_continuous,p_binary,X){
  
  data = list(
    X = X,
    n_patients = n_patients,
    p_continuous = p_continuous,
    p_binary = p_binary
  )
  
  simulation_results = stan('stan files/2_simulate_pk_params.stan',
                             data= data,
                             chains = 1, 
                             iter = 1,
                             seed = 19920908,
                             cores = 1,
                             algorithm = 'Fixed_param')
  
  simulated_params = rstan::extract(simulation_results)
  
  return(simulated_params)
  
}