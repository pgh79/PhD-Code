library(tidyverse)
library(rstan)
source('LOPO.R')
options(mc.cores = parallel::detectCores())

p = 1:36

input_data = rstan::read_rdump('Model Data/core_model_raw_data')
data_sets = map(p, ~LeaveOnePatientOut(input_data,.x,1))

for (i in p){
  data_sets[[i]]$ID = i
}

seed = 19920908
s1 = stan_model(file = 'Models/LOPO_1_condition.stan')

m<-map2(data_sets,p, ~sampling(object = s1, 
             data = .x,
             chains = 1,
             seed = seed,
             chain_id = .y) )

fit = sflist2stanfit(m)


p = rstan::extract(fit)

sampler_params = get_sampler_params(fit, inc_warmup = F)
sapply(m, function(x) mean(as.matrix(x)[,'MAE']) ) %>% hist
