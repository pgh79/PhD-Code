library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
source('2018-09 Apixiban Bayesian Models/stan_utilities.R')

options(mc.cores = parallel::detectCores())

theme_set(theme_minimal())

# Need to rescale the concentration from ng/ML to mg/L since dose is in mg.
apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>% 
  arrange(Subject,Time) %>% 
  filter(Time>0) %>% 
  mutate(Concentration_orig = Concentration,
         Concentration = 10e-3*Concentration_orig) #resalced from ng/ml to mg/L

apixaban.plot = apixaban.data %>% 
  ggplot(aes(Time,Concentration))+
  geom_point()+
  facet_wrap(~Subject)

#---- Prepare Data ----#


t0 = 0
C0 = array(c(0), dim = 1)

D = 2.5
times = sort(unique(apixaban.data$Time))

subject_names = apixaban.data$Subject %>% unique
subjects = apixaban.data$Subject %>% 
  factor %>% 
  as.numeric %>% 
  unique

C_hat = apixaban.data %>% 
  select(Time,Subject,Concentration) %>% 
  spread(Time,Concentration) %>% 
  arrange(Subject) %>% 
  select(-Subject) %>% 
  as.matrix

D = 2.5
N_t = length(times)
N_patients = length(subjects)

N_ut = 100;
utimes = seq(min(times), max(times), length.out = N_ut)


rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "N_patients", "C_hat", 'utimes','N_ut'), 
                  file="2018-09 Apixiban Bayesian Models/Data/model_6_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_6_data.data.R")

file.remove("2018-09 Apixiban Bayesian Models/Data/model_6_data.data.R")

fit = stan('2018-09 Apixiban Bayesian Models/Models 1 Compartment/model_1.6.stan',
           data = input_data,
           chains = 4)

params = rstan::extract(fit)


check_rhat(fit)

#---- plots ---- 
x = params$C

conc = apply(x,c(2,3),median) %>%  
  as.table %>% 
  `dimnames<-`(list(Subject = subject_names[subjects],Time = utimes)) %>% 
  as.data.frame.table(stringsAsFactors = F, responseName = 'Concentration') %>% 
  mutate(Time = as.numeric(Time))


lim = apply(x,c(2,3),function(x) quantile(x,c(0.025, 0.975))) %>%  
  as.table %>% 
  `dimnames<-`( list(Lims = c('ymin','ymax'), Subject = subject_names ,Time = utimes) ) %>% 
  as.data.frame.table(stringsAsFactors = F, responseName = 'Concentration') %>% 
  mutate(Time = as.numeric(Time)) %>% 
  spread(Lims,Concentration)


mcmc = conc %>% 
  left_join(lim)


apixaban.plot+
  geom_line(data = mcmc, color = 'red', aes(group = Subject))+
  geom_ribbon(data = mcmc, aes(ymin = ymin, ymax = ymax,group = Subject), fill = 'red', alpha = 0.5)




mcmc_areas(as.matrix(fit),
           regex_pars = c('k'),
           prob = 0.95, # 80% intervals
           prob_outer = 0.99, # 99%
            )



