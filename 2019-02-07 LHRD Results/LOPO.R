input_data = rstan::read_rdump('Model Data/core_model_raw_data')

LeaveOnePatientOut <- function(data, test_patient, test_time){
    
  #make a place to put this crap
  env = environment()
  list2env(input_data, envir = env)
  
  #test critera for the LOPO
  no_condition_on =  (patient_ID==test_patient)
  
  #Make the test_data
  C_hat_test = C_hat[no_condition_on]
  times_test = times[no_condition_on]
  X_test = X[test_patient,]
  N_test = sum(no_condition_on)
  p = length(X_test)
  
  #Make the train data
  C_hat_train = C_hat[!no_condition_on]
  times_train = times[!no_condition_on]
  N_train = N - N_test
  
  patient_ID_train = as.numeric(as.factor(patient_ID[!no_condition_on]) ) 
  X_train = X[unique(patient_ID[!no_condition_on]),]
  N_patients_train = dim(X_train)[1]
  #Now scale the data
  X_mu = apply(X_train, 2, mean)
  X_sd = apply(X_train, 2, sd)
  X_train[,4:6] =  t( (t(X_train[,4:6]) - X_mu[4:6] )/X_sd[4:6])
  apply(X_train[,4:6], 2, sd)
  
  X_test[4:6] = t( (t(X_test[4:6]) - X_mu[4:6] )/X_sd[4:6])
  
  
  
  test_names = c('C_hat_test','times_test', 'X_test','N_test', 'p')
  train_names = c('C_hat_train','times_train', 'X_train','N_train','D', 'patient_ID_train','N_patients_train')
  names = c(test_names, train_names)
  model_data = mget(names,envir = env)
  
  
  return(model_data)
}

