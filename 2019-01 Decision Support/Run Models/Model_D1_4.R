library(glue)
library(tidyverse)
library(bayesplot)
library(rstan)
options(mc.cores = parallel::detectCores())
source('2018-09 Apixiban Bayesian Models/stan_utilities.R')


#Use model_D1.2 Lognormal likelihood. Seems to be better in terms of diagnostiscs
###########################
which.model <- 'model_D1_4'#
###########################

apixaban.data = read_csv('2019-01 Decision Support/Data/ApixibanExperimentData.csv')

#Covariate information like Weight and Age.  Maybe useful in the regression
covariates = read_csv('2019-01 Decision Support/Data/ApixibanExperimentCovariates.csv')

#Join the covariates with a left join.
apixaban.data = apixaban.data %>% 
  mutate(Concentration_scaled = 0.001*Concentration) %>%  #conver to mg/L for stability
  left_join(covariates) %>% 
  arrange(Subject, Time) %>% 
  filter(Time>0) %>% 
  replace_na(list(Creatinine = mean(covariates$Creatinine,na.rm = T)))

#---- Prepare Data ----

D = 2.5 #mg of drug

times = sort(unique(apixaban.data$Time))
subject_names = apixaban.data$Subject %>% unique

subjects = apixaban.data$Subject %>%
  factor %>%
  as.numeric %>%
  unique

C_hat = apixaban.data %>%
  select(Time, Subject, Concentration_scaled) %>%
  spread(Time, Concentration_scaled) %>%
  arrange(Subject) %>%
  select(-Subject) %>%
  as.matrix

N_t = length(times)

X= apixaban.data %>% 
  distinct(Subject,Sex,Group, Age, Weight,Creatinine) %>% 
  mutate(Age= (Age - mean(Age))/sd(Age), 
         Weight = (Weight - mean(Weight))/sd(Weight),
         Creatinine = (Creatinine - mean(Creatinine))/sd(Creatinine)) %>% 
  model.matrix(~Sex*Group + Weight + Age+ Creatinine, data = .)


X_v = X[,c(1,2,4)]
N_patients = dim(X)[1]
N_covariates = dim(X)[2]


file.name = "2019-01 Decision Support/Data/model_data.R"
rstan::stan_rdump(c(
  "D",
  "times",
  "N_t",
  "N_patients",
  "C_hat",
  'N_covariates',
  'X',
  'X_v'
),file = file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)


#---- Fit Model ----
model = glue('2019-01 Decision Support/Delayed Models 1 Compartment/{which.model}.stan')
fit = rstan::stan(file = model,
                  data = input_data,
                  chains = 10,
                  seed = 10090908,
                  control = list(max_treedepth = 13,adapt_delta = 0.8)
)

check_rhat(fit)

params = rstan::extract(fit)


