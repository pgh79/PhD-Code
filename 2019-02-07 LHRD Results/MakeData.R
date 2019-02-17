library(tidyverse)

apixiban.data = read_csv('2019-02-07 LHRD Results/Data/ApixibanExperimentData.csv')

#Covariate information like Weight and Age.  Maybe useful in the regression
covariates = read_csv('2019-02-07 LHRD Results/Data/ApixibanExperimentCovariates.csv')

#Join the covariates with a left join.
apixiban.data = apixiban.data %>% 
  mutate(Concentration_scaled = 0.001*Concentration) %>%  #conver to mg/L for stability
  left_join(covariates) %>% 
  arrange(Subject, Time) %>% 
  filter(Time>0) %>% 
  replace_na(list(Creatinine = mean(covariates$Creatinine,na.rm = T)))


#Create a plot to verify everything went OK

apixiban.data %>% 
  ggplot(aes(Time,Concentration, color = Group))+
  stat_summary(geom = 'point', fun.y = mean)+
  stat_summary(geom = 'line', fun.y = mean)+ 
  stat_summary(geom = 'pointrange', fun.data = function(x) mean_se(x, 1.96))+
  facet_wrap(~Sex)+
  theme_classic()+
  theme(aspect.ratio = 1/1.61)+
  scale_color_brewer(palette = 'Set1', direction = -1)+
  labs(title = 'Mean Apixiban PK-Curves')

apixiban.data %>% 
  mutate(Subject = as.factor(Subject)) %>% 
  write_csv('2019-02-07 LHRD Results/Data/apixiban_regression_data.csv')
