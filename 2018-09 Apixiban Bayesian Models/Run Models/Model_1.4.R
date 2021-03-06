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
  theme(aspect.ratio = 1)+
  guides(color=F, fill= F)

#---- Prepare Data ----#


t0 = 0
C0 = array(c(0), dim = 1)
D = 2.5

cdata = apixaban.data %>%
  select(Subject, Time, Concentration) %>%
  spread(Time, Concentration) %>%
  arrange(Subject)

subjects = cdata$Subject
N_patients = length(subjects)
C_hat = cdata %>% select(-Subject) %>% as.matrix
A = median(C_hat)
times = apixaban.data$Time %>% unique %>% sort
N_t = length(times)

Xs = apixaban.data %>%
  select(Subject, Sex, Group) %>%
  arrange(Subject) %>%
  distinct() 
  

X = Xs %>% 
  mutate(
    Intercept = 1,
    Sex = Sex %>% factor %>% as.numeric - 1,
    Group = Group %>% factor %>% as.numeric - 1,
    INTXN = Sex * Group) %>% 
  select(-Subject) %>%
  select(Intercept, Sex, Group, INTXN) %>% 
  as.matrix


file.name = "2018-09 Apixiban Bayesian Models/Data/model_data.data.R"
rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "X", "C_hat", 'N_patients','A'), 
                  file=file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)
model = '2018-09 Apixiban Bayesian Models/Models 1 Compartment/model_4.stan'
fit = rstan::stan(file = model,
                  data = input_data,
                  chains=2,
                  control=list(adapt_delta=0.9))

params = rstan::extract(fit)



#----Bayesian Credible Intervals----

y = params$C_ppc
N =  dim(y)[1]

simulations = y %>%
  as.table() %>%
  `dimnames<-`(list(
    Round = 1:N ,
    Subject = subjects,
    Time = times
  )) %>%
  as.data.frame.table(responseName = 'Concentration', stringsAsFactors = F) %>%
  mutate(
    Time = as.numeric(Time),
    Subject = as.numeric(Subject),
    Round = as.numeric(Round)
  ) %>%
  bind_rows(apixaban.data %>% select(Time, Concentration, Subject)) %>%
  mutate(
    Round = replace_na(Round, N + 1),
    Kind = if_else(Round > N, 'Data', 'Simulation')
  ) %>% 
  left_join(Xs, by = 'Subject')


mcmc = simulations %>% 
  group_by(Time,Sex,Group) %>% 
  summarise(C = mean(Concentration),
            ymin = quantile(Concentration, 0.025),
            ymax = quantile(Concentration, 0.975) ) %>% 
  rename("Concentration" = "C")



apixaban.plot+
  geom_line(data = mcmc)+
  geom_ribbon(data = mcmc, aes(ymin = ymin, ymax = ymax), alpha = 0.5)


# mcmc_areas(as.matrix(fit), regex_pars = c('beta'))
