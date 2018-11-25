library(nlmeODE)
library(nlme)
library(tidyverse)

data(Theoph)

apixiban <- read_csv('2018-09 Apixiban Bayesian Models/Data/ApixibanExperimentData.csv') %>% 
            mutate(conc = 0.001*Concentration,
                   Dose = 2.5) %>% 
            as.data.frame() %>% 
            arrange(Subject,Time) %>% 
            mutate(Sex_i = if_else(Sex=='Female',1,0),
                   Group_i = if_else(Group=='Control',1,0))

TheophODE <- apixiban %>% groupedData(conc ~ Time|Subject, data = .)
# TheophODE <- Theoph
TheophODE$Dose[TheophODE$Time != 0] <- 0
TheophODE$Cmt <- rep(1, dim(TheophODE)[1])

OneComp <- list(
  DiffEq = list(dy1dt = ~ -ka * y1 ,
                dy2dt = ~ ka * y1 - ke * y2),
  ObsEq = list(c1 = ~ 0,
               c2 = ~ y2 / CL * ke),
  Parms = c("ka", "ke", "CL"),
  States = c("y1", "y2"),
  Init = list(0, 0)
)

TheophModel <- nlmeODE(OneComp, TheophODE)



Theoph.nlme <- nlme(
  conc ~ TheophModel(ka, ke, CL, Time, Subject),
  data = TheophODE,
  fixed = list(ka+ke~1, CL~Sex_i+Group_i),
  random = pdDiag(ka + CL ~ 1),
  start = c(ka = 0.5, ke = -2.5, CL = -3.2, Sex_i = -1, Group_i = -1),
  control = list(returnObject = TRUE, msVerbose = TRUE),
  verbose = TRUE
)


plot(augPred(Theoph.nlme, level = 0:1))


apixiban$pred = predict(Theoph.nlme)
apixiban %>%
  filter(Time>0) %>%
  mutate(res = conc - pred) %>%
  ggplot(aes(factor(Time),res))+
  geom_jitter(width = 0.25)+
  stat_summary(geom = 'pointrange', fun.data = function(x) mean_se(x,1.96), color = 'orange')+
  facet_grid(Group~Sex)+
  geom_hline(aes(yintercept = 0))+
  theme_bw()

ggsave('NLME_ODE.png',path = '2018-09 Apixiban Bayesian Models/Figures')
