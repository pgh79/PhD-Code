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
         Concentration = 10e-3*Concentration_orig) #resalced from ng/ml to mg/L

apixaban.plot = apixaban.data %>%
  ggplot(aes(Time, Concentration)) +
  geom_point(aes(color = interaction(Sex, Group))) +
  facet_wrap( ~ Subject) +
  scale_color_brewer(palette = 'Set1') +
  theme(legend.position = 'bottom')

#---- Prepare Data ----#


t0 = 0
C0 = array(c(0,0), dim = 2)
D = 2.5

cdata = apixaban.data %>%
  select(Subject, Time, Concentration) %>%
  spread(Time, Concentration) %>%
  arrange(Subject)

subjects = cdata$Subject
N_patients = length(subjects)
C_hat = cdata %>% select(-Subject) %>% as.matrix
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
rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "X", "C_hat", 'N_patients'), 
                  file=file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)
#---- Fit Model -----

model = '2018-09 Apixiban Bayesian Models/Models 2 Compartments/model_10_a.stan'
fit = rstan::stan(file = model,
                  data = input_data,
                  chains=2,
                  control=list(adapt_delta=0.9))

params = rstan::extract(fit)
#---- plots ----

#----Bayesian Credible Intervals ----


y = params$C[ , , ]
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
  )

simulations %>%
  group_by(Time,Subject,Kind) %>% 
  summarise(C = mean(Concentration), 
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
  facet_wrap(~ Subject) +
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
  )

simulations %>%
  group_by(Time,Subject,Kind) %>% 
  summarise(C = mean(Concentration), 
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
  )

simulations %>%
  filter(Round %in% sample(1:N,size = 3)|Round==N+1) %>% 
  ggplot(aes(
    Time,
    Concentration,
    color = Kind,
    group = Round
  )) +
  geom_line()+
  facet_wrap(~ Subject) +
  scale_color_brewer(palette = 'Set1')+
  labs(title = 'Posterior Draws', color = 'Data Source', fill = 'Data Source')


# #---- Stuff ----
# 
# y = params$C_ppc[ , , ]
# N =  dim(y)[1]
# 
# 
# simulations = y %>%
#   as.table() %>%
#   `dimnames<-`(list(
#     Round = 1:N ,
#     Subject = subjects,
#     Time = times
#   )) %>%
#   as.data.frame.table(responseName = 'SimulatedConcentration', stringsAsFactors = F) %>%
#   mutate(
#     Time = as.numeric(Time),
#     Subject = as.numeric(Subject),
#     Round = as.numeric(Round)
#   ) %>%
#   left_join(apixaban.data %>% select(Time, Concentration, Subject) %>% rename('RealConcentration' = 'Concentration'), by = c('Subject','Time'))
# 
# 
# simulations %>%
#   ggplot()+
#   geom_histogram(aes(SimulatedConcentration))+
#   geom_rug(aes(RealConcentration), color = 'red', size = 2)+
#   facet_grid(Time~Subject)
