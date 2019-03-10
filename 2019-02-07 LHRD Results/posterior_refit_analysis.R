library(tidyverse)

files = list.files('LOPO Results', full.names = T)

is.posterior.refit = map_lgl(files, ~grepl('posterior',.x) )

files = files[is.posterior.refit]

datasets = map(files, readRDS)


map_dfc(datasets,~map(.x,'MAE') %>% unlist) %>% 
  `colnames<-`(c('Conditioned on none',
                 'Conditioned on first',
                 'Conditioned on first 2',
                 'Conditioned on last',
                 'Conditioned on last 2')) %>% 
  mutate(patient = 1:36) %>%  
  gather(time_point,val,-patient) %>% 
  group_by(time_point) %>% 
  summarise(err = mean(val)) %>% 
  arrange(err)
