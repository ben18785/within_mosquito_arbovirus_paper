library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores=4)
source("src/r/helper.R")

fit_optimise <- function(stan_data, filename_stan_model,
                         n_opts=3) {
  
  stan_model <- stan_model(filename_stan_model)
  fit <- optimise_repeat(n_opts, stan_data, stan_model)
  
  fit
}