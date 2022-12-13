
create_sampling_diagnostics <- function(fit) {
  
  # extract only estimated parameters, not derived ones
  params_to_extract <- c(
    "l0", "gamma", "k_lm", "a",
    "alpha_m", "k_m", "k_mh",
    "alpha_h", "k_h", "sigma",
    "phi_d", "zeta", "b1", "b2",
    "b3", "b4", "x_star", "eta",
    "chp_vals", "chp_sigma")
  df <- rstan::extract(fit, params_to_extract) %>% 
    as.data.frame()
  
  sum_df <- posterior::summarise_draws(df)
  
  rhats <- sum(sum_df$rhat>1.01, na.rm = TRUE)
  ess_bulks <- sum(sum_df$ess_bulk<400, na.rm = TRUE)
  ess_tails <- sum(sum_df$ess_tail<400, na.rm = TRUE)
  is_converged <- rhats == 0 & ess_bulks == 0 & ess_tails == 0
  
  list(is_converged=is_converged,
       rhat=rhats,
       ess_bulk=ess_bulks,
       ess_tail=ess_tails)
}