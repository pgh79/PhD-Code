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
  ggplot(aes(Time,Concentration, color = interaction(Sex,Group), fill = interaction(Sex,Group)))+
  stat_summary(fun.data = function(x) mean_se(x,1.96), geom = 'errorbar')+
  stat_summary(fun.y = mean, geom = 'line')+
  stat_summary(fun.y = mean, geom = 'point')+
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')+
  facet_grid(Sex~Group)+
  theme(aspect.ratio = 1)

#---- Prepare Data ----#


t0 = 0
C0 = array(c(0), dim = 1)

D = 2.5

times = apixaban.data$Time
utimes = unique(times) %>% sort
utimes = seq(min(utimes), max(utimes), length.out = 100)
sex = apixaban.data$Sex %>% factor %>% as.numeric - 1
group = apixaban.data$Group %>% factor %>% as.numeric -1
N_t = length(times)
N_ut = length(utimes)

C_hat = apixaban.data$Concentration

rstan::stan_rdump(c("t0", "C0", "D", "times", "utimes", "sex", "N_t", "N_ut", "C_hat"), 
                  file="2018-09 Apixiban Bayesian Models/Data/model_4_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_4_data.data.R")

L = stan_fitter('model_4', input_data = input_data)

params = L$params
fit = L$fit



#---- Plots ----

mcmc_areas(as.matrix(fit),  
           regex_pars = 'b',
           prob = 0.95, # 80% intervals
           prob_outer = 0.99, # 99%
           point_est = "mean")

bpl = mcmc_pairs(
  as.matrix(fit),
  regex_pars = 'b',
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)


x = params$C_pred
y = apply(x,c(2,3,4),mean)
dimnames(y) <- list(Sex = c('Female','Male'),Time = utimes,Group =  c('Control','NAFLD'))

conc = as.data.frame.table(y, responseName = 'Concentration', stringsAsFactors = F) %>% 
  mutate(Time = as.numeric(Time))

z = apply(x,c(2,3,4),function(x) quantile(x,probs = c(0.025,0.975)))
dimnames(z) <- list(quant = c('ymin','ymax'),
                    Sex = c('Female','Male'),
                    Time = utimes,
                    Group =  c('Control','NAFLD'))

lims = as.data.frame.table(z, responseName = 'Concentration', stringsAsFactors = F) %>% 
       spread(quant,Concentration) %>% 
       mutate(Time = as.numeric(Time))

mcmc = conc %>% left_join(lims)

apixaban.plot+
  geom_line(data= mcmc)+
  geom_ribbon(data = mcmc, aes(ymin = ymin, ymax = ymax), alpha = 0.25)
