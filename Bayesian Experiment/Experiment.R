library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Function for non-dimensionalized PK
PK_profile <- function(t,alpha){
  y = (1 / (alpha - 1) * exp(t * (alpha - 1)) - 1 / (alpha - 1)) * exp(-alpha * t)
  
  return(y)
}
#Vectorize for ease of use
PK_profile = Vectorize(PK_profile)

#Time points to sample, as in real data
t = c(0.5,1,2,4,6,8,10,12)
#Ratio of excretive and absorptive rates.  Higher means you get rid of it faster
alpha = seq(0.1, 0.5, 0.01)
#Make the dataset with some noise
set.seed(19920908)
data = tibble(
  ID = 1:length(alpha),
  alpha = alpha,
  t = list(t),
  y_true = map2(t,alpha, ~PK_profile(.x,.y)),
  y = map(y_true, ~rlnorm(length(.x), log(.x), 0.15))
) %>% 
unnest()

data %>% 
  ggplot(aes(t,y_true, color = factor(alpha), group = ID))+
  geom_line()+
  guides(color = F)+
  scale_color_brewer(palette = 'RdBu')


#---- Stan First ----
data = data %>% mutate(train = t %in% c(0.5,1))
train = data %>% filter(train)
train_N = dim(train)[1]
train_ID = train$ID
train_t = train$t
train_y = train$y

test = data %>% filter(!train)
test_N = dim(test)[1]
test_ID = test$ID
test_t = test$t
test_y = test$y

input_data = list(
  num_patients = length(alpha),
  train_N = train_N,
  test_N = test_N,
  train_ID = train_ID,
  train_t = train_t,
  train_y = train_y,
  test_ID = test_ID,
  test_t = test_t,
  test_y = test_y
)
fit = stan('model.stan', 
           data = input_data,
           control = list(adapt_delta = 0.99))

p = rstan::extract(fit)

err = p$err %>% apply(.,2,mean)

test$err = err
test %>% write_csv('test_sample_first.csv')

#---- Stan Last ---- 

data = data %>% mutate(train = t %in% c(10,12))
train = data %>% filter(train)
train_N = dim(train)[1]
train_ID = train$ID
train_t = train$t
train_y = train$y

test = data %>% filter(!train)
test_N = dim(test)[1]
test_ID = test$ID
test_t = test$t
test_y = test$y

input_data = list(
  num_patients = length(alpha),
  train_N = train_N,
  test_N = test_N,
  train_ID = train_ID,
  train_t = train_t,
  train_y = train_y,
  test_ID = test_ID,
  test_t = test_t,
  test_y = test_y
)
fit = stan('model.stan', 
           data = input_data,
           control = list(adapt_delta = 0.99))

p = rstan::extract(fit)

err = p$err %>% apply(.,2,mean)

test$err = err
test %>% write_csv('test_sample_last.csv')

#---- Results ----

test1 = read_csv('test_sample_first.csv')
test2 = read_csv('test_sample_last.csv')


x = test1 %>% 
  group_by(ID,alpha) %>% 
  summarise(MAE_first = mean(abs(err)))
y = test2 %>% 
  group_by(ID,alpha) %>% 
  summarise(MAE_last = mean(abs(err)))


x %>% 
  inner_join(y) %>% 
  mutate(diff = MAE_first - MAE_last) %>% 
  ggplot(aes(alpha,diff))+
  geom_line()+
  geom_hline(yintercept = 0)
