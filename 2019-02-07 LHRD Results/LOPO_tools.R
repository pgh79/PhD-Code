
x = rstan::read_rdump('Model Data/core_model')

LeaveOnePatientOut <- function(input_data,i,t = 0){
  
  env = environment()
  list2env(input_data, envir = env)
  
  #Make LOPO data
  LOPOix = which(patient_ID==i)
  times_test = times[LOPOix]
  X_test = X[i,]
  C_hat_test = C_hat[LOPOix]
  p_test = length(X_test)
  N_test = length(LOPOix)
  
  #Training Data
  
  N = N - length(LOPOix)
  N_patients = N_patients - 1
  times = times[-LOPOix]
  C_hat = C_hat[-LOPOix]
  X = X[-i,]
  
  #Here is the hard part
  
  patient_ID = factor(patient_ID[-LOPOix])
  patient_ID = as.numeric(patient_ID)

  train_names = c('N','p','times','N_patients','X', 'D', 'patient_ID','C_hat')
  test_names = c('times_test','X_test','C_hat_test','p_test','N_test')
  names = c(train_names, test_names)
  model_data = mget(names, env)
return(model_data)

}

LeaveOnePatientOut(x,1)
