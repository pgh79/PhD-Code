library(tidyverse)

files = list.files('LOPO Results', full.names = T)

datasets = map(files, readRDS)


map_dfc(datasets,~map(.x,'MAE') %>% unlist) %>% 
  `colnames<-`(c('Conditioned on none',
                 'Conditioned on first',
                 'Conditioned on first 2',
                 'Conditioned on first and last',
                 'Conditioned on last 2')) %>% 
  mutate(patient = 1:36) %>%  
  gather(time_point,val,-patient) %>% 
  group_by(time_point) %>% 
  summarise(err = mean(val)) %>% 
  arrange(err)


  # ggplot(aes(time_point,val))+
  # geom_line(aes(group = patient),alpha = 0.25)+
  # stat_summary(geom = 'point', fun.y = mean, color = 'red')+
  # labs(x = '', y = 'MAE')+
  # theme(aspect.ratio = 1/1.61)
  # 
  # 
