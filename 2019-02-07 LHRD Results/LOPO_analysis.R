library(tidyverse)
library(patchwork)

files = list.files('LOPO Results', full.names = T)

datasets = map(files, readRDS)


map_dfc(datasets,~map(.x,'MAE') %>% unlist) %>% 
  `colnames<-`(c('Conditioned on 0','Conditioned on 1','Conditioned on 2')) %>% 
  mutate(patient = 1:36) %>%  
  gather(time_point,val,-patient)  %>% 
  ggplot(aes(time_point,val))+
  geom_line(aes(group = patient),alpha = 0.25)+
  stat_summary(geom = 'point', fun.y = mean, color = 'red')+
  labs(x = '', y = 'MAE')+
  theme(aspect.ratio = 1/1.61)



datasets %>% 
  map(~map(.x,'relative_error') %>% data_frame, .id = 'pateint')
