logistic_curve <- function(dilution, a, k, b, q) {
  log_concentration = log(1.0 / dilution);
  a + (k - a)/(1 + q * exp(-b * log_concentration))
}

optimise_repeat <- function(n_opts, data_in, stan_model) {
  current_best <- -Inf
  for(i in 1:n_opts) {
    print(paste0("optimisation #: ", i))
    fit <- optimizing(object = stan_model,
                      data = data_in,
                      as_vector=FALSE
    )
    val <- fit$value
    if(val > current_best) {
      current_best <- val
      best_fit <- fit
    }
  }
  best_fit
}