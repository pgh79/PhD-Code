library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
source('2018-09 Apixiban Bayesian Models/stan_helpers.R')
options(mc.cores = parallel::detectCores())

theme_set(theme_minimal())

## This model uses a normal likelihood.  Seems to fit better.

# Need to rescale the concentration from ng/ML to mg/L since dose is in mg.
apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>% 
  arrange(Subject,Time) %>% 
  filter(Time>0) %>% 
  mutate(Concentration_orig = Concentration,
         Concentration = 10e-3*Concentration_orig) #resalced from ng/ml to mg/L


apixaban.plot =  apixaban.data %>% 
  ggplot(aes(Time,Concentration))+
  stat_summary(fun.data = function(x) mean_se(x,1.96), geom = 'errorbar')+
  stat_summary(fun.y = mean, geom = 'point')+
  stat_summary(fun.y = mean, geom = 'line')

#Prepare the data to be passed into stan
C_hat = apixaban.data$Concentration
C0 = array(c(0), dim = 1)

t0 = 0
times = apixaban.data$Time
utimes = sort(unique(times))
N_t = length(times)
N_ut = length(utimes)

D = 2.5 #mg/L


rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "C_hat",'utimes'), 
                  file="2018-09 Apixiban Bayesian Models/Data/model_1_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_1_data.data.R")



#Fit Model ----
L = stan_fitter('model_1', input_data = input_data)
fit = L$fit
params = L$params


color_scheme_set("darkgray")

mcmc_areas(as.matrix(fit), pars = c('k','k_a','V','sigma'))

mcmc_pairs(as.matrix(fit), pars = c('k_a','k','V'))

C = apply(params$C,2, median)
C_interval = apply(params$C,2, function(x) quantile(x,c(0.025,0.975))) %>% t()

preds = data_frame(Concentration = C, 
                   C_L = C_interval[,1], 
                   C_U = C_interval[,2],
                   Time = utimes)

preds %>% 
  ggplot()+
  geom_line(aes(Time, Concentration,color = 'Bayesian'))+
  geom_ribbon(aes(x = Time, ymin = C_L, ymax = C_U), alpha = 0.25, fill = 'red' )+
  stat_summary(data = apixaban.data, aes(Time, Concentration, color = 'Mean' ), geom = 'line', fun.data = function(x) mean_se(x,1.96))+
  stat_summary(data = apixaban.data, aes(Time, Concentration), geom = 'errorbar', fun.data = function(x) mean_se(x,1.96))+
  labs(title = 'Credible Interval')+
  scale_color_manual('', 
                     breaks = c('Bayesian','Mean'),
                     values = c('red','black')
  )


C = apply(params$C,2, median)
C_interval = apply(params$C,2, function(x) quantile(x,c(0.025,0.975))) %>% t()



preds = data_frame(Concentration = C, 
                   C_L = C_interval[,1], 
                   C_U = C_interval[,2],
                   Time = utimes)


apixaban.plot+
  stat_summary(fun.y = median, geom = 'point', shape = 21, size = 5)+
  geom_line(data = preds, color = 'red')+
  geom_ribbon(data= preds,aes(ymin = C_L, ymax = C_U), fill = 'red', alpha= 0.5)



