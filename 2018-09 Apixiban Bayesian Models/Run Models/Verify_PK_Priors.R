library(rstan)
library(tidyverse)
library(bayesplot)
options(mc.cores = parallel::detectCores())

apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv')

#Covariate information like Weight and Age.  Maybe useful in the regression
covariates = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentCovariates.csv')

#Join the covariates with a left join.
apixaban.data = apixaban.data %>% 
  mutate(Concentration_scaled = 0.001*Concentration) %>%  #conver to mg/L for stability
  left_join(covariates) %>% 
  arrange(Subject, Time) %>% 
  filter(Time>0) %>% 
  replace_na(list(Creatinine = mean(covariates$Creatinine,na.rm = T)))


X= apixaban.data %>% 
  distinct(Subject,Sex,Group, Age, Weight,Creatinine) %>% 
  mutate(Age= (Age - mean(Age))/sd(Age), 
         Weight = (Weight - mean(Weight))/sd(Weight),
         Creatinine = (Creatinine - mean(Creatinine))/sd(Creatinine)) %>% 
  model.matrix(~Sex*Group + Weight + Age+ Creatinine, data = .)

N_patients = dim(X)[1]
N_covariates = dim(X)[2]

input_data = list(N_covariates = N_covariates, 
                  N_patients = N_patients,
                  X = X)


stan_program = '2018-09 Apixiban Bayesian Models/Delayed Models 1 Compartment/simulate_PKparams.stan'
sims = stan(stan_program,
            data = input_data,
            init = 19920908,
            algorithm = 'Fixed_param',
            chains =1,
            iter=5000
)


params = rstan::extract(sims)

hist(params$PKP[params$PKP[,1]<1,1], breaks = 50)
quantile(params$PKP[,1], probs = seq(0,1,0.1)) %>% log
