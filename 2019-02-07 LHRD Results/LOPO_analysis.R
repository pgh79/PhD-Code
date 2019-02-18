library(tidyverse)


##No Conditioning
loadRDS('LOPO Results/results_0_condition')

mae = purrr::map(results,'MAE') %>% unlist
hist(mae)
           
rmse = purrr::map(results,'RMSE') %>% unlist
hist(rmse)

relerr = purrr::map(results, 'relative_error') %>% 
         unlist() %>% 
          matrix(ncol = 36)

matplot(relerr, type = 'l')
