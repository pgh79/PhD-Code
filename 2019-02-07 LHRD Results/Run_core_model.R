library(tidyverse)
library(rstan)
library(loo)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source('stan_utilities.R')


df = read_csv('Data/apixiban_regression_data.csv') %>% 
  mutate(SubjectID = as.numeric(factor(Subject)) ) %>% 
  select(Subject,SubjectID, everything()) %>% 
  mutate_at(vars(Age,Weight,Creatinine), scale)


input_data <- read_rdump('Model Data/core_model_scaled_data')

fit = stan('Models/core_model.stan',
           data = input_data,
           chains = 12,
           seed = 10090908,
           control = list(max_treedepth = 13,adapt_delta = 0.8),
           verbose = F)

check_all_diagnostics(fit)

p = rstan::extract(fit)




# Vis posterior ----

bayesplot::mcmc_areas(as.matrix(fit), regex_pars = 'SIGMA')

b = colnames(as.matrix(fit))
b = b[map_lgl(b,~grepl('BETA',.x))]
as.matrix(fit)[,b] %>% cov() %>% corrplot::corrplot()

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
  geom_smooth(aes(Concentration_scaled, pred), method= 'lm')+
  theme(aspect.ratio = 1)


library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(fit, merge_chains = FALSE)

# as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- relative_eff(exp(log_lik_1)) 

loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
