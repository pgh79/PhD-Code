library(tidyverse)
library(rstan)
library(loo)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source('2018-09 Apixiban Bayesian Models/stan_utilities.R')

#TODO 
#Write non-herarchical regression

#TODO
#Write completely heirarhcical model.