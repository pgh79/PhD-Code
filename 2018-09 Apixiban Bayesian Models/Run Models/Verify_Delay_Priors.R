library(rstan)
library(tidyverse)

stan_program = '2018-09 Apixiban Bayesian Models/Delayed Models 1 Compartment/simulate_delays.stan'

sims = stan(stan_program,
            data = list(A=2.5, B = 2.5, ymin = 10, Z =1.5),
            init = 19920908,
            algorithm = 'Fixed_param',
            chains =1,
            iter=5000
            )


params = extract(sims)


hist(params$delay/2)
