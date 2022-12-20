
fit_mcmc <- function(opt, stan_model_file, data_in, n_iterations,
                     n_chains) {
  
  rstan_options(auto_write = TRUE)
  options(mc.cores=4)
  
  init_fun <- function(...){
    list(l0 = opt$par$l0,
         gamma = opt$par$gamma,
         k_lm = opt$par$k_lm,
         a = opt$par$a,
         alpha_m = opt$par$alpha_m,
         k_m = opt$par$k_m,
         k_mh = opt$par$k_mh,
         alpha_h = opt$par$alpha_h,
         k_h = opt$par$k_h,
         sigma = opt$par$sigma,
         phi_m = opt$par$phi_m,
         phi_d = opt$par$phi_d,
         zeta = opt$par$zeta,
         b1 = opt$par$b1,
         b2 = opt$par$b2,
         b3 = opt$par$b3,
         b4 = opt$par$b4,
         x_star = opt$par$x_star,
         eta = opt$par$eta,
         chp_vals = opt$par$chp_vals,
         chp_sigma = opt$par$chp_sigma)}
  
  fit <- stan(stan_model_file,
              data=data_in,
              iter=n_iterations,
              chains=n_chains,
              init=init_fun)
  
  fit
}


fit_mcmc_kappa <- function(opt, stan_model_file, data_in, n_iterations,
                     n_chains) {
  
  rstan_options(auto_write = TRUE)
  options(mc.cores=4)
  
  init_fun <- function(...){
    list(l0 = opt$par$l0,
         gamma = opt$par$gamma,
         k_lm = opt$par$k_lm,
         a = opt$par$a,
         alpha_m = opt$par$alpha_m,
         k_m = opt$par$k_m,
         k_mh = opt$par$k_mh,
         alpha_h = opt$par$alpha_h,
         k_h = opt$par$k_h,
         sigma = opt$par$sigma,
         phi_m = opt$par$phi_m,
         phi_d = opt$par$phi_d,
         zeta = opt$par$zeta,
         b1 = opt$par$b1,
         b2 = opt$par$b2,
         b3 = opt$par$b3,
         b4 = opt$par$b4,
         x_star = opt$par$x_star,
         eta = opt$par$eta,
         chp_vals = opt$par$chp_vals,
         chp_sigma = opt$par$chp_sigma)}
  
  dilutions <- rev(c(1, 5, 12, 20, 25)) # need reverse since need smallest kappa for highest dilution
  data_in$dilutions_lookup <- match(data_in$dilutions, dilutions)
  data_in$dilutions_binary_lookup <- match(data_in$dilutions_binary, dilutions)
  
  fit <- stan(stan_model_file,
              data=data_in,
              iter=n_iterations,
              chains=n_chains,
              init=init_fun)
  
  fit
}


