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

apixiban.plot = apixaban.data %>% 
  ggplot(aes(Time,Concentration))+
  stat_summary(fun.data = function(x) mean_se(x,1.96), geom = 'errorbar')+
  stat_summary(fun.y = mean, geom = 'line')+
  stat_summary(fun.y = mean, geom = 'point')+
  scale_color_brewer(palette = 'Set1')+
  facet_wrap(~Sex)


#Prepare the data to be passed into stan
t0 = 0
C0 = array(c(0), dim = 1)

D = 2.5

times = apixaban.data$Time
utimes = unique(times) %>% sort
sex = apixaban.data$Sex %>% factor %>% as.numeric - 1
N_t = length(times)
N_ut = length(utimes)

C_hat = apixaban.data$Concentration

rstan::stan_rdump(c("t0", "C0", "D", "times", "utimes", "sex", "N_t", "N_ut", "C_hat"), 
                  file="2018-09 Apixiban Bayesian Models/Data/model_3_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_3_data.data.R")

L = stan_fitter('model_3', input_data = input_data)

params = L$params
fit = L$fit



#Plots----

#Female is first column
pred_c = apply(params$C_pred,c(2,3), median) %>% 
  t() %>% 
  as.data.frame() %>% 
  cbind(utimes) %>% 
  rename(Female = V1, Male = V2, Time = utimes) %>% 
  gather(Sex,Concentration,-Time)

c_low = apply(params$C_pred,c(2,3), function(x) quantile(x, probs = 0.025)) %>%
  t() %>% 
  as.data.frame() %>% 
  cbind(utimes) %>% 
  rename(Female = V1, Male = V2, Time = utimes) %>% 
  gather(Sex,ymin,-Time)


c_high = apply(params$C_pred,c(2,3), function(x) quantile(x, probs = 0.975)) %>%
  t() %>% 
  as.data.frame() %>% 
  cbind(utimes) %>% 
  rename(Female = V1, Male = V2, Time = utimes) %>% 
  gather(Sex,ymax,-Time)

pred = pred_c %>% 
  left_join(c_low, by = c('Time','Sex')) %>% 
  left_join(c_high, by = c('Time','Sex'))


pix = apixiban.plot+
  geom_line(data = pred, aes(color= Sex))+
  geom_ribbon(data = pred,aes(color = Sex, fill = Sex, ymin = ymin, ymax = ymax),alpha = 0.5)+
  theme(aspect.ratio = 1)

pix$layers = pix$layers[c(4,5,1,2,3)]

pix

mcmc_areas(as.matrix(fit),  
           pars = c('beta_k[1]','beta_k[2]','beta_k_a[1]','beta_k_a[2]','beta_V[1]','beta_V[2]'),
           prob = 0.95, # 80% intervals
           prob_outer = 0.95, # 99%
           point_est = "mean")

mcmc_pairs(
  as.matrix(fit),
  pars = c('beta_k[1]','beta_k[2]','beta_k_a[1]','beta_k_a[2]','beta_V[1]','beta_V[2]'), 
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)
