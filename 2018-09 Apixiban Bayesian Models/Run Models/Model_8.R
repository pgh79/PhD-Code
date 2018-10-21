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
         Concentration = 10e-3*Concentration_orig)  #resalced from ng/ml to mg/L


apixaban.plot =  apixaban.data %>% 
  ggplot(aes(Time,Concentration))+
  stat_summary(fun.data = function(x) mean_se(x,1.96), geom = 'errorbar')+
  stat_summary(fun.y = mean, geom = 'point')+
  stat_summary(fun.y = mean, geom = 'line')

#---- Data Prep ----

C_hat = apixaban.data$Concentration 
C0 = array(c(0,0), dim = 2)
t0 = 0

times = unique(apixaban.data$Time)
N_t = length(times)
itimes = match(apixaban.data$Time, times)
N_it = length(itimes)

D = 2.5 #mg/L

file.name = "2018-09 Apixiban Bayesian Models/Data/model_data.data.R"
rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "C_hat",'itimes','N_it'), 
                  file=file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)

#---- Fit Model ----
model = '2018-09 Apixiban Bayesian Models/Models 2 Compartments/model_8.stan'
fit = rstan::stan(file = model,
                  data = input_data,
                  chains=2,
                  control=list(adapt_delta=0.9))

params = rstan::extract(fit)


x = params$C_pred
y = apply(x, 2, mean)
yl = apply(x, 2, function(x) quantile(x,0.025))
yu = apply(x, 2, function(x) quantile(x,0.975))
conc =data.frame(Time = times,
                 Concentration = y,
                 ymin = yl,
                 ymax = yu)


apixaban.plot+
  geom_line(data = conc, color = 'red')+
  geom_ribbon(data = conc, aes(Time, ymin = ymin, ymax = ymax), fill = 'red', alpha = 0.5)

