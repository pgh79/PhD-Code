library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
source('2018-09 Apixiban Bayesian Models/stan_helpers.R')

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
  geom_line()+
  facet_wrap(~Subject, nrow = 6)

apixaban.plot
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


rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "N_patients", "C_hat"), 
                  file="2018-09 Apixiban Bayesian Models/Data/model_5_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_5_data.data.R")

L = stan_fitter('model_5', input_data = input_data)

params = L$params
fit = L$fit


x = params$C

conc = apply(x,c(2,3),median) %>%  
        as.table %>% 
        `dimnames<-`(list(Subject = subject_names[subjects],Time = times)) %>% 
        as.data.frame.table(stringsAsFactors = F, responseName = 'Concentration') %>% 
        mutate(Time = as.numeric(Time))


lim = apply(x,c(2,3),function(x) quantile(x,c(0.025, 0.975))) %>%  
        as.table %>% 
       `dimnames<-`( list(Lims = c('ymin','ymax'), Subject = subject_names ,Time = times) ) %>% 
        as.data.frame.table(stringsAsFactors = F, responseName = 'Concentration') %>% 
         mutate(Time = as.numeric(Time)) %>% 
        spread(Lims,Concentration)
    

mcmc = conc %>% 
       left_join(lim)


results = apixaban.plot+
  geom_line(data = mcmc, color = 'red')+
  geom_ribbon(data = mcmc, aes(ymin = ymin, ymax = ymax), fill = 'red', alpha = 0.5)


ggsave('no_pooling.pdf', results, path = '~/Desktop/', height = 8, width = 13)



preds = apply(x,c(2,3),mean)

(preds - C_hat) %>% hist
