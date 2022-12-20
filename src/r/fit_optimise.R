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


fit_optimise_kappa <- function(stan_data, filename_stan_model,
                         n_opts=3) {
  
  dilutions <- rev(c(1, 5, 12, 20, 25)) # need reverse since need smallest kappa for highest dilution
  stan_data$dilutions_lookup <- match(stan_data$dilutions, dilutions)
  stan_data$dilutions_binary_lookup <- match(stan_data$dilutions_binary, dilutions)
  
  stan_model <- stan_model(filename_stan_model)
  fit <- optimise_repeat(n_opts, stan_data, stan_model)
  
  fit
}
