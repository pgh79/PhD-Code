library(tidyverse)

experiment_names = tibble(name = list.files('LOPO Data', full.names = F)) %>% 
                   mutate(experiment = as.character(seq_along(name)) )

results = list.files('LOPO Data', full.names = T) %>% 
  map(readRDS) %>% 
  map_df(., ~.x %>% 
          select(1,7) %>% 
          unnest() %>% 
          group_by(patients) %>% 
          summarise(mae = mean(abs(err)))
                    , .id = 'experiment') 



results %>% 
  inner_join(experiment_names) %>% 
  mutate(name = stringr::str_split(name,'.RDS') %>% 
                  map_chr(., 1)
         ) %>%
  group_by(name) %>% 
  summarise(m = mean(mae)) %>% 
  arrange(m)
  
