functions {
  real C_anal(real t, real D, real V, real k_a, real k) {
    return (D / V) * (k_a / (k_a - k))
            * ( exp(- k * t) - exp(-k_a * t) );
  }
}
data {
  real t0;    // Initial time in days;
  real C0[1]; // Initial concentration at t0 in mg/L

  real D;   // Total dosage in mg

  int<lower=1> N_t;
  real times[N_t];   // Measurement times in days
  int<lower=1> N_ut;
  real utimes[N_ut];
  
  real<lower=0> A;


  // Measured concentrations in effect compartment in mg/L
  int<lower=1> N_patients;
  real C_hat[N_patients, N_t];
  matrix[N_patients,4] X;
}
parameters {
 
 // phi = log(k_a);
<<<<<<< HEAD
 real<lower=0> ka;
=======
 real phi;
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
 
 //psi_1 = log(V);
 vector[4] BETA_psi_1;
 real<lower=0> SIGMA_psi_1;
 real psi_1[N_patients];
 
 //psi_2 = log(CL);
 //CL = V*k
 vector[4] BETA_psi_2;
 real<lower=0> SIGMA_psi_2;
<<<<<<< HEAD
 real<lower=psi_1> psi_2[N_patients];
=======
 real psi_2[N_patients];
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
 
 real<lower=0, upper=1> alpha;
 real<lower=0> SIGMA_y;
  

}

transformed parameters {
  
  real C[N_patients, N_t];
  
  vector[N_patients] MU_psi_1;
  vector[N_patients] MU_psi_2;
  
  real<lower=0> k[N_patients];
  real<lower=0> V[N_patients];
<<<<<<< HEAD
=======
  real<lower=0> ka;
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
  
  MU_psi_1 = X*BETA_psi_1;
  MU_psi_2 = X*BETA_psi_2;
  
<<<<<<< HEAD
  
=======
  ka = phi;
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
  for (n in 1:N_patients) {
    V[n] = exp(psi_1[n]);
    k[n] = exp(psi_2[n] - psi_1[n]);
    for (t in 1:N_t)
      C[n,t] = C_anal(times[t], D, V[n] , ka, k[n]);
  }
}

model {
<<<<<<< HEAD
  // 
  // BETA_phi ~ normal(0,2);
  // SIGMA_phi ~ cauchy(0,1);
  // phi ~ normal(MU_phi, SIGMA_phi);
  
  ka ~ cauchy(0,1);
  
  BETA_psi_1 ~ normal(0,2);
  SIGMA_psi_1 ~ cauchy(0,1);
  psi_1 ~student_t(10,MU_psi_1, SIGMA_psi_1);
=======

  
  phi ~ lognormal(0,1);
  
  BETA_psi_1 ~ normal(0,2);
  SIGMA_psi_1 ~ cauchy(0,1);
  psi_1 ~student_t(7,MU_psi_1, SIGMA_psi_1);
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
  
  
  BETA_psi_2 ~ normal(0,2);
  SIGMA_psi_2 ~ cauchy(0,1);
<<<<<<< HEAD
  psi_2 ~student_t(10,MU_psi_2, SIGMA_psi_2);
=======
  psi_2 ~student_t(7,MU_psi_2, SIGMA_psi_2);
>>>>>>> f21a23f9e02e9a38369c560bd798488a68cee9f3
  
  SIGMA_y ~ cauchy(0,1);
  alpha ~ uniform(0,1);
  
  
  for (n in 1:N_patients)
  for (t in 1:N_t)
    C_hat[n, t] ~ normal(C[n, t], pow(C[n,t]/A,2*alpha)*SIGMA_y);
}

generated quantities {
  
  real C_ppc[N_patients, N_t];
  for (n in 1:N_patients)
    for (t in 1:N_t)
      C_ppc[n, t] = normal_rng(C[n, t], pow(C[n,t]/A,2*alpha)*SIGMA_y);

}
