library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
library(glue)
library(corrplot)
# source('2018-09 Apixiban Bayesian Models/stan_helpers.R')





#----WhichModel----
which.model = 'model_7_a'
#------------------

options(mc.cores = parallel::detectCores())

theme_set(theme_minimal())

# Need to rescale the concentration from ng/ML to mg/L since dose is in mg.
apixaban.data = read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>%
  arrange(Subject, Time) %>%
  filter(Time > 0) %>%
  mutate(Concentration_orig = Concentration,
         Concentration = 0.001* Concentration_orig) #resalced from ng/ml to mg/L

apixaban.plot = apixaban.data %>%
  ggplot(aes(Time, Concentration)) +
  geom_point(aes(color = interaction(Sex, Group))) +
  facet_wrap( ~ Subject) +
  scale_color_brewer(palette = 'Set1') +
  theme(legend.position = 'bottom')



f.plot = apixaban.data %>% 
  ggplot(aes(Time,Concentration, color = Group))+
  stat_summary(geom = 'pointrange',fun.data = function(x) mean_se(x,1.96))+
  stat_summary(geom = 'line',fun.y = 'mean')+
  facet_wrap(~Sex)+
  scale_color_brewer(palette = 'Set1')+
  theme(aspect.ratio = 1)

  ggsave('apixaban_viz.png', path = '2018-09 Apixiban Bayesian Models/Figures/')


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

A = median(C_hat)

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


file.name = "2018-09 Apixiban Bayesian Models/Data/model_data.data.R"
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
  'X',
  'A'
),file = file.name)

rstan::stan_rdump(c("t0", "C0", "D", "times", "N_t", "X", "C_hat", 'N_patients'), file=file.name)

input_data <- read_rdump(file.name)

file.remove(file.name)


#---- Fit Model ----
model = glue('2018-09 Apixiban Bayesian Models/Models 1 Compartment/{which.model}.stan')
fit = rstan::stan(file = model,
                  data = input_data,
                  chains=2)

params = rstan::extract(fit)


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

bci<-simulations %>%
  left_join(apixaban.data %>% select(Subject,Sex,Group)) %>% 
  filter(Group == 'NAFLD', Sex == "Female") %>% 
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

ggsave(glue('{which.model}_BCI.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi =400,
        )

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

ppi<-simulations %>%
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
  labs(title = 'Posterior Predictive Interval', color = 'Data Source', fill = 'Data Source')+
  theme(legend.position = 'None')

ggsave(glue('{which.model}_PPI.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)

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

ggsave(glue('{which.model}_DFP.png'),path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)

# ---- Pred vs Data ----


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




pvd<-simulations %>% 
  group_by(Time,Subject,Kind) %>%  
  summarise(C = mean(Concentration)) %>% 
  spread(Kind,C) %>% 
  ggplot(aes(Simulation, Data - Simulation))+
  geom_point()+
  geom_smooth()+
  #stat_function(fun = function(x) qnorm(0.975,0,(x/0.3833)^(0.66)*0.11), color = 'red')+
  #stat_function(fun = function(x) qnorm(0.025,0,(x/0.3833)^(0.66)*0.11), color = 'red')+
  geom_hline(aes(yintercept = 0), color = 'black')
  theme(aspect.ratio = 1)

ggsave(glue('{which.model}_PVD.png'),path = '2018-09 Apixiban Bayesian Models/Figures/')
  
#---- err ----

y = params$C[, ,]
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
  mutate(err = Data - Simulation) %>%
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

err

ggsave(glue('{which.model}_err.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)




#----Histograms----

params$V%>% 
  apply(2,mean) %>% 
  as.table() %>% 
  `dimnames<-`(list(Subject = subject_names)) %>% 
  as.data.frame.table(responseName = 'V', stringsAsFactors = F) %>% 
  mutate(Subject = as.numeric(Subject)) %>% 
  left_join(apixaban.data %>% select(Subject,Group,Sex) %>% distinct) %>% 
  ggplot(aes(V))+
  geom_histogram(color = 'white')+
  facet_grid(Group~Sex)+
  labs(x = 'V (Litres)')+
  theme_bw()


ggsave(glue('{which.model}_V.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)


params$k%>% 
  apply(2,mean) %>% 
  as.table() %>% 
  `dimnames<-`(list(Subject = subject_names)) %>% 
  as.data.frame.table(responseName = 'V', stringsAsFactors = F) %>% 
  mutate(Subject = as.numeric(Subject)) %>% 
  left_join(apixaban.data %>% select(Subject,Group,Sex) %>% distinct) %>% 
  ggplot(aes(V))+
  geom_histogram(color = 'white')+
  facet_grid(Group~Sex)+
  labs(x = expression(paste(k,' (mg/L/Hour)')))+
  theme_bw()


ggsave(glue('{which.model}_k.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)


params$k_a%>% 
  apply(2,mean) %>% 
  as.table() %>% 
  `dimnames<-`(list(Subject = subject_names)) %>% 
  as.data.frame.table(responseName = 'V', stringsAsFactors = F) %>% 
  mutate(Subject = as.numeric(Subject)) %>% 
  left_join(apixaban.data %>% select(Subject,Group,Sex) %>% distinct) %>% 
  ggplot(aes(V))+
  geom_histogram(color = 'white')+
  facet_grid(Group~Sex)+
  labs(x = expression(paste(k[a],' (mg/L/Hour)')))+
  theme_bw()


ggsave(glue('{which.model}_k_a.png'),
       path = '2018-09 Apixiban Bayesian Models/Figures/',
       dpi = 400)


#---
simulations %>%
group_by(Time, Subject, Kind) %>%
  summarise(Concentration = mean(Concentration)) %>%
  spread(Kind, Concentration) %>%
  ungroup %>% 
  left_join(apixaban.data) %>% 
  mutate(res = Concentration - Simulation) %>% 
  summarise(sqrt(mean(res^2)))
