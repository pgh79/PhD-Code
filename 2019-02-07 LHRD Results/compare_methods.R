library(tidyverse)

PKcurve<-function(t,V,ka,k){
  prop_to = (2.5/V)*(ka/(ka-k))
  
  return(prop_to*(exp(-k*t) - exp(-ka*t)))
}

find_t_crit<-function(model){
  
  k = coef(model)['k']
  ka = coef(model)['ka']
  
  names(k)=NULL
  names(ka)=NULL
  
  top = log(k) - log(ka)
  bot = k - ka
  
  return(top/bot)
}


do_model<-function(t,y){
  model = nls(y~PKcurve(t,V,ka,k),
              start = list(V = 5, ka = 0.4, k = 0.5),
              lower = list(V = 0,ka = 0,k = 0),
              control = list(printEval = T,maxiter = 1000, warnOnly = T),
              algorithm = 'port')
  
  return(model)
}

t<-c(0.5,1,2,4,6,8,10,12)
d = readRDS(file = 'LOPO Data/condition_on_first_and_last.RDS') 