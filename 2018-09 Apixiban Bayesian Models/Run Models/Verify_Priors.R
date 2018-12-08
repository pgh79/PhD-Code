library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())

theme_set(theme_minimal())

# Need to rescale the concentration from ng/ML to mg/L since dose is in mg.
apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>% 
  arrange(Subject,Time) %>% 
  filter(Time>0) %>% 
  mutate(Concentration_orig = Concentration,
         Concentration = 10e-3*Concentration_orig) #resalced from ng/ml to mg/L

fit = rstan::stan(file = '2018-09 Apixiban Bayesian Models//Models/veryify_priors.stan',iter=100,chains=1, seed=19920908, algorithm="Fixed_param")

params = rstan::extract(fit)

utimes = seq(0.5,12.5, length.out = 100)

sims= data_frame(time = utimes) %>% 
  bind_cols(t(params$C) %>% as.data.frame()) %>% 
  gather(round,val,-time)

sims %>% 
  ggplot(aes(time,val))+
  geom_line(aes( group= round), size = 0.5)+
  geom_point(data = apixaban.data, aes(Time, Concentration), color = 'red')+
  scale_y_continuous(limits = c(0,5))
