library(tidyverse)
library(rstan)
library(loo)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source('stan_utilities.R')


df = read_csv('Data/apixiban_regression_data.csv') %>% 
  mutate(SubjectID = as.numeric(factor(Subject)) ) %>% 
  select(Subject,SubjectID, everything()) %>% 
  mutate_at(vars(Age,Weight,Creatinine,BMI),scale)


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

file.name = "Model Data/LOPO_no_conditioning"
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
input_data <- read_rdump(file.name)

fit = stan('Models/LOPO_no_conditioning.stan',
           data = input_data,
           chains = 12,
           seed = 10090908,
           control = list(max_treedepth = 13,adapt_delta = 0.8),
           verbose = F)

check_all_diagnostics(fit)

p = rstan::extract(fit)

C_pred = p$C %>% apply(., 2 , mean)
C_q = p$C_ppc %>% 
  apply(., 2 , function(x) quantile(x, c(0.025, 0.975))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(ymin = `2.5%`, ymax = `97.5%`)

dfp = p$C_ppc[sample(1:dim(p$C_ppc)[1],3),] %>% 
      t() %>% 
      as.data.frame() %>% 
      bind_cols(df %>% select(Time, SubjectID)) %>% 
      gather(round,val,-SubjectID,-Time)

df$pred = C_pred

df %>% 
  cbind(C_q) %>% 
  ggplot()+
  geom_ribbon(aes(Time, ymin = ymin, ymax = ymax), fill = 'red', alpha = 0.5)+
  geom_line(aes(Time,pred), color = 'red') +
  geom_line(data = dfp, aes(Time,val,group = round), color = 'blue')+
  geom_line(aes(Time,Concentration_scaled))+
  facet_wrap(~SubjectID, scale = 'free_y')

df %>% 
  cbind(C_q) %>% 
  ggplot()+
  geom_pointrange(aes(Concentration_scaled, pred, ymin = ymin, ymax = ymax), alpha = 0.25)+
  geom_abline()+
  geom_smooth(aes(Concentration_scaled, pred))

