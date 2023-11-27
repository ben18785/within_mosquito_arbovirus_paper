
logistic_curve_conc_simple <- function(concentration, a, k, b, q) {
  
  log_concentration = log(concentration)
  
  a + (k - a)/(1 + q * exp(-b * log_concentration))
}

rnorm_truncated_single <- function(mu, sigma, lower=0, upper=Inf) {
  x <- rnorm(1, mu, sigma)
  if(x < lower | x > upper)
    x <- rnorm_truncated_single(mu, sigma, lower, upper)
  x
}

rnorm_truncated <- function(n, mu, sigma, lower=0, upper=Inf) {
  x <- vector(length = n)
  for(i in seq_along(x))
    x[i] <- rnorm_truncated_single(mu, sigma, lower, upper)
  x
}

prior_predictive_phi <- function(fit) {
  concentration <- seq(0, sqrt(2), length.out=200)^2
  n <- 1000
  b <- rnorm_truncated(n, 1, 1)
  q <- rnorm_truncated(n, 0, 0.5)
  
  nreps <- 1000
  for(i in 1:nreps) {
    probs <- logistic_curve_conc_simple(concentration, 0, 1, b[i], q[i])
    tmp <- tibble(x=concentration, prob=probs) %>% 
      mutate(iteration=i)
    if(i == 1)
      big_df <- tmp
    else
      big_df <- big_df %>% bind_rows(tmp)
  }
  
  df_summary <- big_df %>% 
    group_by(x) %>% 
    summarise(
      lower=quantile(prob, 0.025),
      middle=quantile(prob, 0.5),
      upper=quantile(prob, 0.975),
      )
  
  # get estimated parameters
  b <- mean(rstan::extract(fit, "b3")[[1]])
  q <- mean(rstan::extract(fit, "b4")[[1]])
  probs <- logistic_curve(concentration, 0, 1, b, q)
  df_est <- tibble(x=concentration, middle=probs)
  
  df_summary %>% 
    ggplot(aes(x=x, y=middle)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,
                fill="blue") +
    geom_line() +
    geom_line(data=df_est, colour="orange") +
    scale_x_sqrt() +
    scale_y_continuous(labels = scales::percent) +
    xlab("Concentration") +
    ylab("Positive")
}

prior_predictive_logistic_growth <- function(alpha_est, kappa_est, mu_alpha, sigma_alpha,
                                             mu_kappa, sigma_kappa) {

  logistic_growth <- function(t, y0, alpha, kappa) {
    kappa * y0 / (y0 + (kappa - y0) * exp(-alpha * t))
  }
  
  times <- seq(0, 20, 0.1)
  n <- 1000
  alpha <- rnorm_truncated(n, mu_alpha, sigma_alpha)
  kappa <- rnorm_truncated(n, mu_kappa, sigma_kappa)

  nreps <- 1000
  y0 <- 0.0001
  for(i in 1:nreps) {
    y <- logistic_growth(times, y0, alpha[i], kappa[i])
    tmp <- tibble(time=times, prob=y) %>% 
      mutate(iteration=i)
    if(i == 1)
      big_df <- tmp
    else
      big_df <- big_df %>% bind_rows(tmp)
  }
  
  df_summary <- big_df %>% 
    group_by(time) %>% 
    summarise(
      lower=quantile(prob, 0.025),
      middle=quantile(prob, 0.5),
      upper=quantile(prob, 0.975),
    )
  
  # get estimated parameters
  y <- logistic_growth(times, y0, alpha_est, kappa_est)
  df_est <- tibble(time=times, middle=y)
  
  df_summary %>% 
    ggplot(aes(x=time, y=middle)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,
                fill="blue") +
    geom_line() +
    geom_line(data=df_est, colour="orange") +
    xlab("Time, DPI") +
    ylab("Normalised DENV titer")
}
