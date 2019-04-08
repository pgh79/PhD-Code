library(tidyverse)
library(Metrics)



experiment_names = tibble(name = list.files('LOPO Data', full.names = F)) %>% 
                   mutate(experiment = as.character(seq_along(name)) )

results = list.files('LOPO Data', full.names = T) %>% 
  map(readRDS) %>% 
  map_df(., ~.x %>% 
          select(patients,err) %>% 
          unnest() %>% 
          group_by(patients) %>% 
          summarise(rmse = mean(abs(err)))
                    , .id = 'experiment') 



results %>% 
  inner_join(experiment_names) %>% 
  mutate(name = stringr::str_split(name,'.RDS') %>% 
                  map_chr(., 1)
         ) %>%
  group_by(name) %>% 
  summarise(error.metric = mean(rmse)) %>% 
  arrange(desc(error.metric))
  

results = list.files('LOPO Data', full.names = T)


readRDS('LOPO Data/condition_on_last.RDS') %>% 
  select(patients,ypred,ytest) %>% 
  mutate(t = list(c(0.5,1,2,4,6,8,10,12))) %>% 
  unnest() %>% 
  gather(type,val,-patients,-t) %>% 
  ggplot(aes(t,val, color = type))+
  geom_line()+
  facet_wrap(~patients, scale = 'free_y')+
  labs(title = 'condition on last')

readRDS('LOPO Data/condition_on_first.RDS') %>% 
  select(patients,ypred,ytest) %>% 
  mutate(t = list(c(0.5,1,2,4,6,8,10,12))) %>% 
  unnest() %>% 
  gather(type,val,-patients,-t) %>% 
  ggplot(aes(t,val, color = type))+
  geom_line()+
  facet_wrap(~patients, scale = 'free_y')+
  labs(title = 'condition on first')

