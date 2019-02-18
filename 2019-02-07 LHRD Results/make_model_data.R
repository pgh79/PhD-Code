library(tidyverse)
library(rstan)
library(loo)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source('stan_utilities.R')


##---- Scaled Data ----
df = read_csv('Data/apixiban_regression_data.csv') %>% 
  mutate(SubjectID = as.numeric(factor(Subject)) ) %>% 
  select(Subject,SubjectID, everything()) %>% 
  mutate_at(vars(Age,Weight,Creatinine), scale)

D = 2.5
N_patients = df$Subject %>% unique %>% length
patient_ID = df$SubjectID

X = df %>% 
  distinct(SubjectID, .keep_all = T) %>% 
  model.matrix(~Group*Sex + Age + Weight + Creatinine, data =.)
N = dim(df)[1]
p = dim(X)[2]
times = df$Time
C_hat = df$Concentration_scaled

file.name = "Model Data/core_model_scaled_data"
rstan::stan_rdump(c(
  "D",
  "N_patients",
  "patient_ID",
  "X",
  "N",
  'p',
  'times',
  'C_hat'
),file = file.name)


#---- Raw Data ----

df = read_csv('Data/apixiban_regression_data.csv') %>% 
  mutate(SubjectID = as.numeric(factor(Subject)) ) %>% 
  select(Subject,SubjectID, everything())

D = 2.5
N_patients = df$Subject %>% unique %>% length
patient_ID = df$SubjectID

X = df %>% 
  distinct(SubjectID, .keep_all = T) %>% 
  model.matrix(~Group*Sex + Age + Weight + Creatinine, data =.)
N = dim(df)[1]
p = dim(X)[2]
times = df$Time
C_hat = df$Concentration_scaled

file.name = "Model Data/core_model_raw_data"
rstan::stan_rdump(c(
  "D",
  "N_patients",
  "patient_ID",
  "X",
  "N",
  'p',
  'times',
  'C_hat'
),file = file.name)

