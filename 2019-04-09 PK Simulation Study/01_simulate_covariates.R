
simulate_covariates<-function(n_patients, p_continuous, p_binary, theta){
#This script will simulate covariates for the PK study from a stan model
  

  
  simulation_parameters = list(n_patients = n_patients,
                                p_continuous = p_continuous,
                                p_binary = p_binary,
                                theta = theta)

  simulation_results = stan('stan files/1_simulate_data.stan',
                            data =simulation_parameters,
                            iter = 1,
                            chains = 1,
                            seed = 19920908,
                            cores = 1,
                            algorithm = 'Fixed_param',
                            )
  
  simulation_params = rstan::extract(simulation_results)
  
  
  #Check matrix makes sense
  design_matrix = simulation_params$X[1,,]  
  
  return(design_matrix)
}
