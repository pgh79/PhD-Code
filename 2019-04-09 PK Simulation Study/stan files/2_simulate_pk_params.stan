data{
  int<lower=1> n_patients; //Number of patients to simulate.  Rows of design matrix
  int<lower=1> p_continuous; //How many continuous variables?
  int<lower=1> p_binary; //How many binary variables?
  matrix[n_patients,1+p_continuous+p_binary] X; //design matrix
}
transformed data{
  int<lower=1> p_covariates = 1+p_continuous+p_binary;
  
}
model{}
generated quantities{
  // Step 1: Generate reg_coefs, coefficients of the regression
  // Step 2: Generate Linear predictor Xreg_coefs
  // Step 3: Generate a Cholesky factor for the covariance between PK parameters
  // Step 4: Generate some standard normals
  // Step 5: Generate PK params exp(eta+AZ)
    
  matrix[p_covariates,3] reg_coefs = rep_matrix(0, p_covariates,3);
  matrix[n_patients,3]  linear_predictors;
  cov_matrix[3] Sigma; //Covariance
  matrix[3,3] sigma_chol_factor;
  matrix[n_patients,3] pk_params;
  
  Sigma = rep_matrix(0,3,3);
  Sigma[1,1] = 0.0025;
  Sigma[2,2] = 0.1;
  Sigma[3,3] = 0.1;
  
  
  //Effects for V.  Remainder will be sparse
  reg_coefs[1,1] = log(5);
  //Effects for remaining parms
  for(i in 2:3){
    for(j in 1:p_covariates){
      reg_coefs[j,i] = normal_rng(0, 0.2);
    }
  }
  
  linear_predictors = X*reg_coefs;
  
  sigma_chol_factor = cholesky_decompose(Sigma);

  for (i in 1:n_patients){
    pk_params[i] = exp(multi_normal_cholesky_rng(linear_predictors[i],sigma_chol_factor))';
  }
}
