library(rstan)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('LOPO.R')
#Lopo 0
patient = 2
model_data = rstan::read_rdump('Model Data/core_model_raw_data')
pdata = LeaveOnePatientOut(model_data,patient,c(0.5,1,2,4,6,8,10,12))

fit_1 = stan('Models/LOPO_0_condition.stan',
     data = pdata,
     chains = 4,
     refresh = -1)

p_1 = rstan::extract(fit_1)


MAE_1 = mean(p_1$MAE)

#---- Condition on 1 ----
pdata = LeaveOnePatientOut(model_data,patient,c(4,6,8,10,12))
pdata$ID = patient

fit_2 = stan('Models/LOPO_1_condition.stan',
           data = pdata,
           chains = 4,
           refresh = -1)
p_2 = rstan::extract(fit_2)

MAE_2 = mean(p_2$MAE)

print(MAE_1)
print(MAE_2)

plot(pdata$times_test, pdata$C_hat_test)
lines(pdata$times_test, apply(p_2$C_pred,2,mean))
