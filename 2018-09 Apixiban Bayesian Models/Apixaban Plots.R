library(tidyverse)


conc = read_csv('ApixibanConc.csv') %>% 
       na.omit() %>% 
       gather(Subject,Concentration,-Time) %>% 
       mutate(Subject = as.factor(Subject))

groups = read_csv('ApixabanGroups.csv') %>% 
          mutate(Subject = as.factor(Subject))



data = conc %>% 
      inner_join(groups) %>% 
      mutate(Sex = if_else(Sex=='m','Male','Female'))

data %>% 
  ggplot(aes(Time,Concentration,color = interaction(Group,Sex)))+
  stat_summary(fun.data = function(x) mean_se(x,1.96), size = 1, geom = 'errorbar')+
  stat_summary(fun.y = mean, geom = 'line')+
  stat_summary(fun.y = mean, geom = 'point', size = 2)+
  scale_color_brewer(palette = 'Set1')+
  theme_minimal()

data %>% 
  ggplot(aes(Time,Concentration, color= Group))+
  geom_point(size = 2, alpha = 0.25)+
  scale_color_brewer(palette = 'Set1')+
  facet_wrap(~Sex, nrow = 2)+
  theme_minimal()

data %>% 
  write_csv('ApixibanExperimentData.csv')
