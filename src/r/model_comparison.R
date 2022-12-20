get_full_log_likelihood <- function(stanfit) {
  log_likelihood <- loo::extract_log_lik(stanfit, "log_likelihood")
  log_likelihood_binary <- loo::extract_log_lik(stanfit, "log_likelihood_binary")
  log_likelihood_chp <- loo::extract_log_lik(stanfit, "log_likelihood_chp")
  cbind(log_likelihood, log_likelihood_binary, log_likelihood_chp)
}

model_comparison <- function(stanfit_base, stanfit_kappa) {
  
  log_likelihood_base <- get_full_log_likelihood(stanfit_kappa_base)
  log_likelihood_kappa <- get_full_log_likelihood(stanfit_kappa)
  
  loo(log_likelihood_base, log_likelihood_kappa)
}