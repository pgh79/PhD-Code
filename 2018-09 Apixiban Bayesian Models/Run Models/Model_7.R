library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
source('2018-09 Apixiban Bayesian Models/stan_helpers.R')

options(mc.cores = parallel::detectCores())

theme_set(theme_minimal())

# Need to rescale the concentration from ng/ML to mg/L since dose is in mg.
apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>%
  arrange(Subject, Time) %>%
  filter(Time > 0) %>%
  mutate(Concentration_orig = Concentration,
         Concentration = 10e-3 * Concentration_orig) #resalced from ng/ml to mg/L

apixaban.plot = apixaban.data %>%
  ggplot(aes(Time, Concentration)) +
  geom_point(aes(color = interaction(Sex, Group))) +
  facet_wrap( ~ Subject) +
  scale_color_brewer(palette = 'Set1') +
  theme(legend.position = 'bottom')


#---- Prepare Data ----#


t0 = 0
C0 = array(c(0), dim = 1)

D = 2.5
times = sort(unique(apixaban.data$Time))

subject_names = apixaban.data$Subject %>% unique
subjects = apixaban.data$Subject %>%
  factor %>%
  as.numeric %>%
  unique

C_hat = apixaban.data %>%
  select(Time, Subject, Concentration) %>%
  spread(Time, Concentration) %>%
  arrange(Subject) %>%
  select(-Subject) %>%
  as.matrix

D = 2.5
N_t = length(times)
N_patients = length(subjects)

N_ut = 50

utimes = seq(min(times), max(times), length.out = N_ut)

X = apixaban.data %>%
  mutate(Intercept = 1) %>%
  distinct(Subject, Intercept, Sex, Group) %>%
  mutate(
    Sex = Sex %>% factor %>% as.numeric - 1,
    Group = Group %>% factor %>% as.numeric - 1,
    ITXN = Sex * Group
  ) %>%
  select(-Subject) %>%
  as.matrix

X = X[, c('Intercept', 'Sex', 'Group', 'ITXN')]



rstan::stan_rdump(c(
  "t0",
  "C0",
  "D",
  "times",
  "N_t",
  "N_patients",
  "C_hat",
  'utimes',
  'N_ut',
  'X'
),
file = "2018-09 Apixiban Bayesian Models/Data/model_7_data.data.R")

input_data <- read_rdump("2018-09 Apixiban Bayesian Models/Data/model_7_data.data.R")

L = stan_fitter('model_7', input_data = input_data)
params = L$params
fit = L$fit

file.remove("2018-09 Apixiban Bayesian Models/Data/model_7_data.data.R")
#---- plots ----

#----Bayesian Credible Intervals ----


y = params$C[ , , ]
N =  dim(y)[1]

simulations = y %>%
  as.table() %>%
  `dimnames<-`(list(
    Round = 1:N ,
    Subject = subject_names,
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
  )

simulations %>%
  group_by(Time,Subject,Kind) %>% 
  summarise(C = median(Concentration), 
            ymin = quantile(Concentration, 0.025),
            ymax = quantile(Concentration, 0.975)) %>% 
  ggplot(aes(
    Time,
    C,
    color = Kind,
    fill = Kind,
    ymin = ymin,
    ymax = ymax
  )) +
  geom_line()+
  geom_ribbon(alpha = 0.5)+
  facet_wrap(~ Subject, scale = 'free_y') +
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')+
  labs(title = 'Bayesian Credible Interval', color = 'Data Source', fill = 'Data Source')


#----Posterior Predictive Interval----


y = params$C_ppc[ , , ]
N =  dim(y)[1]

simulations = y %>%
  as.table() %>%
  `dimnames<-`(list(
    Round = 1:N ,
    Subject = subject_names,
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
  )

simulations %>%
  group_by(Time,Subject,Kind) %>% 
  summarise(C = median(Concentration), 
            ymin = quantile(Concentration, 0.025),
            ymax = quantile(Concentration, 0.975)) %>% 
  ggplot(aes(
    Time,
    C,
    color = Kind,
    fill = Kind,
    ymin = ymin,
    ymax = ymax
  )) +
  geom_line()+
  geom_ribbon(alpha = 0.5)+
  facet_wrap(~ Subject, scale = 'free_y') +
  scale_color_brewer(palette = 'Set1')+
  scale_fill_brewer(palette = 'Set1')+
  labs(title = 'Posterior Predictive Interval', color = 'Data Source', fill = 'Data Source')

#----Draws From Posterior----

y = params$C_ppc[ , , ]
N =  dim(y)[1]

simulations = y %>%
  as.table() %>%
  `dimnames<-`(list(
    Round = 1:N ,
    Subject = subject_names,
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
  )

simulations %>%
  filter(Round %in% sample(1:N,size = 10)|Round==N+1) %>% 
  ggplot(aes(
    Time,
    Concentration,
    color = Kind,
    group = Round
  )) +
  geom_line()+
  facet_wrap(~ Subject, scale = 'free_y') +
  scale_color_brewer(palette = 'Set1')+
  labs(title = 'Posterior Draws', color = 'Data Source', fill = 'Data Source')


