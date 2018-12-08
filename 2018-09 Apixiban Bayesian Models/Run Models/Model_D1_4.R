library(glue)
library(tidyverse)
library(bayesplot)
library(rstan)
options(mc.cores = parallel::detectCores())
source('2018-09 Apixiban Bayesian Models/stan_utilities.R')


#Use model_D1.2 Lognormal likelihood. Seems to be better in terms of diagnostiscs
###########################
which.model <- 'model_D1_4'#
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

D = 2.5 #mg of drug

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

N_t = length(times)

X= apixaban.data %>% 
  distinct(Subject,Sex,Group, Age, Weight,Creatinine) %>% 
  mutate(Age= (Age - mean(Age))/sd(Age), 
         Weight = (Weight - mean(Weight))/sd(Weight),
         Creatinine = (Creatinine - mean(Creatinine))/sd(Creatinine)) %>% 
  model.matrix(~Sex*Group + Weight + Age+ Creatinine, data = .)


X_v = X[,1:2]
N_patients = dim(X)[1]
N_covariates = dim(X)[2]


file.name = "2018-09 Apixiban Bayesian Models/Data/model_data.R"
rstan::stan_rdump(c(
  "D",
  "times",
  "N_t",
  "N_patients",
  "C_hat",
  'N_covariates',
  'X',
  'X_v'
),file = file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)


#---- Fit Model ----
model = glue('2018-09 Apixiban Bayesian Models/Delayed Models 1 Compartment/{which.model}.stan')
fit = rstan::stan(file = model,
                  data = input_data,
                  chains = 4,
                  control = list(max_treedepth = 13,adapt_delta = 0.8)
)

check_rhat(fit)

params = rstan::extract(fit)



#---- plots ----
#----Bayesian Credible Intervals ----


y = params$C_ppc[ , , ]*1000 #Converting back to ng/ml
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
  facet_wrap(~ Subject) +
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

errs<-simulations %>%
  group_by(Time, Subject, Kind) %>%
  summarise(Concentration = mean(Concentration)) %>% 
  spread(Kind, Concentration) %>%
  ungroup %>% 
  left_join(apixaban.data) %>% 
  mutate(err =(Data - Simulation))

err = errs %>% 
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
  ggplot(aes(sample = res))+
  geom_qq_line(distribution = stats::qlnorm ) +
  stat_qq(distribution = stats::qlnorm)


#--- POSTERIOR DRAWS ---
y = params$C_ppc[ , , ]*1000
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

dfp<-simulations %>%
  mutate(alpha = if_else(Kind=='Data',1,0.5)) %>% 
  filter(Round %in% sample(1:N,size = 5)|Round==N+1) %>% 
  ggplot(aes(
    Time,
    Concentration,
    color = Kind,
    group = Round
  )) +
  geom_line()+
  facet_wrap(~Subject, scale = 'free_y') +
  scale_color_brewer(palette = 'Set1')+
  theme(legend.position = 'bottom')+
  scale_alpha_identity()+
  labs(title = 'Posterior Draws', color = 'Data Source', fill = 'Data Source')+
  theme(legend.position = 'None')

dfp

