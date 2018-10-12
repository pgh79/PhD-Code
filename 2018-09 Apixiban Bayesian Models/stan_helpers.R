stan_fitter <- function(model_name, input_data){
  
  #model fit
  model_dir = '2018-09 Apixiban Bayesian Models/Models/'
  model = paste(model_dir, model_name,'.stan', sep = '')
  fit = rstan::stan(file = model,
                    data = input_data,
                    chains = 4)
  
  #Save fit
  data_dir = '2018-09 Apixiban Bayesian Models/Data/'
  data_name = paste(data_dir,model_name,'_fit.rds', sep = '')
  saveRDS(fit, file = data_name)
  
  #return params
  params = rstan::extract(fit)
  
  return(list(params = params, fit = fit))
  
}