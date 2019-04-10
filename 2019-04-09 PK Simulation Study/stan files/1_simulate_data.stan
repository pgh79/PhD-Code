//This stan file is used to generate the covariates
//for simulated patients

//Output:
// X - matrix of columns n_patients x (p_covariates+1)
// Some will be continuous
// Some will be binary

data{
  int<lower=1> n_patients; //Number of patients to simulate.  Rows of design matrix
  int<lower=1> p_continuous; //How many continuous variables?
  int<lower=1> p_binary; //How many binary variables?
  real<lower=0,upper=1> theta ; //What is probability of 1 for binary variables?
}
transformed data{
  //Variable for how many columns in total.  Add 1 at beggining for a column of ones (i.e. intercept)
  int<lower=1> p_covariates = 1+p_continuous + p_binary;
  //Mean vector for continuous variables.
  //Need mean vector because simulating from MVN.
  //Since we can always center variables, choose 0.
  vector[p_continuous]  mu= rep_vector(0,p_continuous);
}
model{}
generated quantities {

  cov_matrix[p_continuous] Sigma_continuous; //cov matrix
  matrix[n_patients,p_covariates] X; //Design matrix
  vector[p_continuous] x_continuous; //somewhere to place draws from MVN
  
  
  //Random covariance matrix.  Actally just a correlation matrix
  Sigma_continuous = lkj_corr_rng(p_continuous,0.5);
  
  for (i in 1:n_patients){
    
    //Draw the continuous variables
    x_continuous =multi_normal_rng(mu,Sigma_continuous);
    //Fill in the column of ones
    X[i,1] =1;
    //Fill in the continuous variables
    X[i,2:(1+p_continuous)] = x_continuous';
    //In the remaining columns, fill with int.
    //Do this because Stan doesn't like ints in vectors for some reason.
    for(j in (1+p_continuous+1):p_covariates){
      X[i,j] = bernoulli_rng(theta);
    }
  }

}
