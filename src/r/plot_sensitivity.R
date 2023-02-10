
rhs_complex <- function(t, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    x_t = x_star + exp(-eta * t) * (x_0 - x_star)
    if(t < t_refeed)
      x_t <- x_star + exp(-eta*t) * (x_0 - x_star)
    else
      x_t <- x_star + exp(-eta*(t - t_refeed)) * (x_0 - x_star)
    dl_dt = -gamma * l - x_t * (k_lm * l^2) / (a^2 + l^2)
    dm_dt = x_t * (k_lm * l^2) / (a^2 + l^2) + alpha_m * m * (1 - m / k_m) - k_mh * (x_t - x_star) * m
    dh_dt = k_mh * (x_t - x_star) * m + alpha_h * h * (1 - h/k_h)
    list(c(dl_dt, dm_dt, dh_dt))
  })
  
}


simulate_model <- function(ode_parameters, inits, times=seq(0, 15, 0.01)) {
  library(deSolve)
  sol <- ode(y=inits, func=rhs_complex, parms=ode_parameters, times=times) %>% 
    as.data.frame()
  
  sol
}


get_summary_parameters <- function(fit, x0, t_refeed, summary_func=mean) {
  
  gamma <- summary_func(rstan::extract(fit, "gamma")[[1]])
  k_lm <- summary_func(rstan::extract(fit, "k_lm")[[1]])
  a <- summary_func(rstan::extract(fit, "a")[[1]])
  eta <- summary_func(rstan::extract(fit, "eta")[[1]])
  x_star <- summary_func(rstan::extract(fit, "x_star")[[1]])
  alpha_m <- summary_func(rstan::extract(fit, "alpha_m")[[1]])
  k_m <- summary_func(rstan::extract(fit, "k_m")[[1]])
  k_mh <- summary_func(rstan::extract(fit, "k_mh")[[1]])
  alpha_h <- summary_func(rstan::extract(fit, "alpha_h")[[1]])
  k_h <- summary_func(rstan::extract(fit, "k_h")[[1]])
  l_0 <- summary_func(rstan::extract(fit, "l0")[[1]])
  
  parameters <- c(
    gamma=gamma,
    k_lm=k_lm,
    a=a,
    eta=eta,
    x_star=x_star,
    alpha_m=alpha_m,
    k_m=k_m,
    k_mh=k_mh,
    alpha_h=alpha_h,
    k_h=k_h,
    l0=l_0,
    x_0=x0,
    t_refeed=t_refeed
  )
  
  parameters
}

calculate_sensitivity <- function(parameter_name, multipliers, fit, stan_data) {
  
  base_parameters <- get_summary_parameters(fit, stan_data$x_0, 100)
  inits <- c(l=as.numeric(base_parameters["l0"]), m=0, h=0)
  parameter_values <- base_parameters[parameter_name] * multipliers
  times <- seq(0.01, 10, 0.01)
  base_sol <- simulate_model(base_parameters, inits, times=times) %>% 
    mutate(multiplier=1)
  
  for(i in seq_along(parameter_values)) {
    parameters <- base_parameters
    parameters[parameter_name] <- parameter_values[i]
    sol <- simulate_model(parameters, inits, times=times) %>% 
      mutate(multiplier=multipliers[i])
    if(i == 1)
      big_df <- sol
    else
      big_df <- big_df %>% bind_rows(sol)
  }
  big_df <- big_df %>% 
    bind_rows(base_sol)
  
  big_df %>% 
    mutate(parameter=parameter_name)
}

plot_sensitivities_midgut_invasion <- function(fit, list_stan_datasets) {
  
  stan_data <- list_stan_datasets$stan_data
  multipliers <- seq(0.1, 11, 1)
  df_gamma <- calculate_sensitivity("gamma", multipliers, fit, stan_data)
  df_klm <- calculate_sensitivity("k_lm", multipliers, fit, stan_data) %>% 
    mutate(parameter="kappa[lm]")
  df_a <- calculate_sensitivity("a", multipliers, fit, stan_data)
  
  df_all <- df_gamma %>% 
    bind_rows(
      df_klm,
      df_a
    ) %>% 
    dplyr::select(time, m, multiplier, parameter)
  
  titer_multiplier <- list_stan_datasets$dataset_denv_sum %>% 
    filter(tissue=="midgut") %>% 
    pull(overall_denv_titer)
  df_all <- df_all %>% 
    mutate(m = m * titer_multiplier)
  
  df_all %>% 
    ggplot(aes(x=time, y=m, colour=multiplier, group=multiplier)) +
    geom_line() +
    geom_line(data=df_all %>% filter(multiplier==1),
              colour="black", linetype=2) +
    theme(
      legend.position = "none"
    ) +
    xlab("DPI") +
    ylab("Midgut DENV titer") +
    theme_bw() +
    scale_color_viridis_c("Multiplier",
                          trans = "log",
                          breaks=c(0.1, 1, 10)) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    facet_wrap(~parameter, labeller = label_parsed)
}


