library(rstan)
rstan_options(auto_write= T)

#This script will simulate covariates for the PK study from a stan model

n_patients = 50
p_continuous = 2
p_binary = 2
theta = 0.5

data = list(n_patients = n_patients,
            p_continuous = p_continuous,
            p_binary = p_binary,
            theta = theta)

sims = stan('stan files/1_simulate_data.stan',
           data =data,
           iter = 1,
           chains = 1,
           algorithm = 'Fixed_param')

params = rstan::extract(sims)


#Check matrix makes sense
params$X[1,,]
