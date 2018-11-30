library(glue)
library(tidyverse)
library(rstan)
library(measurements)
options(mc.cores = parallel::detectCores())

###########################
which.model <- 'model_D1.1'#
###########################

apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv')

#Covariate information like Weight and Age.  Maybe useful in the regression
covariates = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentCovariates.csv')

#Join the covariates with a left join.
apixaban.data = apixaban.data %>% 
                mutate(Concentration_scaled = 0.001*Concentration) %>%  #conver to mg/L for stability
                left_join(covariates) %>% 
                arrange(Subject, Time) %>% 
                filter(Time>0) %>% 
                replace_na(list(Creatinine = mean(covariates$Creatinine,na.rm = T)))

#---- Prepare Data ----



D = 2.5

times = sort(unique(apixaban.data$Time))
subject_names = apixaban.data$Subject %>% unique

subjects = apixaban.data$Subject %>%
  factor %>%
  as.numeric %>%
  unique

C_hat = apixaban.data %>%
  select(Time, Subject, Concentration_scaled) %>%
  spread(Time, Concentration_scaled) %>%
  arrange(Subject) %>%
  select(-Subject) %>%
  as.matrix

A = median(C_hat)
N_t = length(times)

X= apixaban.data %>% 
  distinct(Subject,Sex,Group, Age, Weight,Creatinine) %>% 
  mutate(Age= Age - median(Age), 
         Weight = (Weight - mean(Weight))/sd(Weight),
         Creatinine = (Creatinine - mean(Creatinine))/sd(Creatinine)) %>% 
  model.matrix(~Sex*Group + Weight + Age+ Creatinine, data = .)

N_patients = dim(X)[1]
p = dim(X)[2]


file.name = "2018-09 Apixiban Bayesian Models/Data/model_data.data.R"
rstan::stan_rdump(c(
  "D",
  "times",
  "N_t",
  "N_patients",
  "C_hat",
  'p',
  'X',
  'A'
),file = file.name)

rstan::stan_rdump(c(
  "D",
  "times",
  "N_t",
  "N_patients",
  "C_hat",
  'p',
  'X',
  'A'
), file=file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)


#---- Fit Model ----
model = glue('2018-09 Apixiban Bayesian Models/Delayed Models 1 Compartment/{which.model}.stan')
fit = rstan::stan(file = model,
                  data = input_data,
                  chains = 2,
                  # control = list(adapt_delta = 0.99),
                  seed=19920908)

params = rstan::extract(fit)


#---- plots ----

#----Bayesian Credible Intervals ----


y = params$C[ , , ]*1000
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

bci<-simulations %>%
  left_join(apixaban.data %>% select(Subject,Sex,Group)) %>% 
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
  labs(title = 'Bayesian Credible Interval', color = 'Data Source', fill = 'Data Source')+
  theme(legend.position = 'None')

bci

#---- err ----

y = params$C[, ,]*1000
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
  mutate(Round = replace_na(Round, N + 1),
         Kind = if_else(Round > N, 'Data', 'Simulation'))

err<-simulations %>%
  group_by(Time, Subject, Kind) %>%
  summarise(Concentration = mean(Concentration)) %>% 
  spread(Kind, Concentration) %>%
  ungroup %>% 
  left_join(apixaban.data) %>% 
  mutate(err =(Data - Simulation)) %>%
  ggplot(aes(factor(Time), err)) +
  # geom_line(aes(group = Subject)) +
  geom_hline(aes(yintercept = 0))+
  geom_jitter(width = 0.2)+
  stat_summary(geom = 'pointrange', 
               fun.data = function(x) mean_se(x,1.96), 
               color = 'red')+
  facet_grid(Group~Sex)+
  theme_bw()+
  labs(x = 'Time', y = 'Actual - Predicted')


err+theme(aspect.ratio = 1)

#---
simulations %>%
  group_by(Time, Subject, Kind) %>%
  summarise(Concentration = mean(Concentration)) %>%
  spread(Kind, Concentration) %>%
  ungroup %>% 
  left_join(apixaban.data) %>% 
  mutate(res = Concentration - Simulation) %>% 
  # summarise(sqrt(mean(res^2)))
  ggplot(aes(sample = res))+geom_qq_line() + stat_qq()


