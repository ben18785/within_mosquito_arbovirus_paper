
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

calculate_sensitivity <- function(parameter_name, multipliers, fit, stan_data,
                                  times, t_refeed=100) {
  
  base_parameters <- get_summary_parameters(fit, stan_data$x_0, t_refeed)
  inits <- c(l=as.numeric(base_parameters["l0"]), m=0, h=0)
  parameter_values <- base_parameters[parameter_name] * multipliers
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
  lower <- log10(0.1)
  upper <- log10(11)
  multipliers <- 10^seq(lower, upper, length.out=10)
  times <- seq(0.01, 15, 0.01)
  df_gamma <- calculate_sensitivity("gamma", multipliers, fit, stan_data, times=times)
  df_klm <- calculate_sensitivity("k_lm", multipliers, fit, stan_data, times=times) %>% 
    mutate(parameter="k[lm]")
  df_a <- calculate_sensitivity("a", multipliers, fit, stan_data, times=times)
  
  df_all <- df_gamma %>% 
    bind_rows(
      df_klm,
      df_a
    ) %>% 
    dplyr::select(time, m, h, multiplier, parameter) %>% 
    pivot_longer(c(m, h)) %>% 
    rename(tissue=name) %>% 
    mutate(
      tissue=if_else(tissue=="m", "midgut", "legs")
    )
  
  titer_multipliers <- list_stan_datasets$dataset_denv_sum
  
  df_all <- df_all %>% 
    left_join(titer_multipliers) %>% 
    mutate(value = value * overall_denv_titer) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
    mutate(parameter=as.factor(parameter)) %>% 
    mutate(parameter=fct_relevel(parameter, "gamma", "a", "k[lm]"))
  
  df_all %>% 
    ggplot(aes(x=time, y=value, colour=multiplier, group=multiplier)) +
    geom_line() +
    geom_line(data=df_all %>% filter(multiplier==1),
              colour="black", linetype=2) +
    theme(
      legend.position = "none"
    ) +
    xlab("DPI") +
    ylab("DENV titer") +
    theme_bw() +
    scale_color_viridis_c("Multiplier",
                          trans = "log",
                          breaks=c(0.1, 1, 10)) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14)) +
    facet_grid(vars(tissue), vars(parameter), labeller = label_parsed,
               scales = "free")
}

plot_sensitivities_single_double_feed <- function(fit, list_stan_datasets) {
  stan_data <- list_stan_datasets$stan_data
  multipliers <- c(0.1, 1, 10)
  times <- seq(0.01, 15, 0.01)
  df_single_eta <- calculate_sensitivity("eta", multipliers, fit, stan_data, times=times, t_refeed=100) %>% 
    mutate(type="single")
  df_double_eta <- calculate_sensitivity("eta", multipliers, fit, stan_data, times=times, t_refeed=3) %>% 
    mutate(type="double")
  df_single_alpha <- calculate_sensitivity("alpha_m", multipliers, fit, stan_data, times=times, t_refeed=100) %>% 
    mutate(type="single") %>% 
    mutate(parameter="alpha[m]")
  df_double_alpha <- calculate_sensitivity("alpha_m", multipliers, fit, stan_data, times=times, t_refeed=3) %>% 
    mutate(type="double") %>% 
    mutate(parameter="alpha[m]")
  df_single_kappa <- calculate_sensitivity("k_mh", multipliers, fit, stan_data, times=times, t_refeed=100) %>% 
    mutate(type="single") %>% 
    mutate(parameter="k[mh]")
  df_double_kappa <- calculate_sensitivity("k_mh", multipliers, fit, stan_data, times=times, t_refeed=3) %>% 
    mutate(type="double") %>% 
    mutate(parameter="k[mh]")
  df_both <- df_single_eta %>% 
    bind_rows(
      df_single_alpha,
      df_single_kappa,
      df_double_eta,
      df_double_alpha,
      df_double_kappa
      )
  
  df_all <- df_both %>% 
    dplyr::select(time, m, h, multiplier, parameter, type) %>% 
    pivot_longer(c(m, h)) %>% 
    rename(tissue=name) %>% 
    mutate(
      tissue=if_else(tissue=="m", "midgut", "legs")
    )
  titer_multipliers <- list_stan_datasets$dataset_denv_sum
  
  df_all <- df_all %>% 
    left_join(titer_multipliers) %>% 
    mutate(value = value * overall_denv_titer) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs"))
  
  df_all %>% 
    filter(tissue == "legs") %>% 
    mutate(type=as.factor(type)) %>% 
    mutate(type=fct_relevel(type, "single", "double")) %>% 
    ggplot(aes(x=time, y=value, colour=type, group=type)) +
    geom_line() +
    theme(
      legend.position = "none"
    ) +
    xlab("DPI") +
    ylab("DENV titer") +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14)) +
    facet_grid(vars(parameter), vars(multiplier), labeller = label_parsed,
               scales = "free")
}

calculate_sensitivity_2d <- function(
    parameter_names, multipliers1, multipliers2, fit, stan_data,
    times, t_refeed=100) {
  
  base_parameters <- get_summary_parameters(fit, stan_data$x_0, t_refeed)
  inits <- c(l=as.numeric(base_parameters["l0"]), m=0, h=0)
  
  parameter1_values <- base_parameters[parameter_names[1]] * multipliers1
  parameter2_values <- base_parameters[parameter_names[2]] * multipliers2
  base_sol <- simulate_model(base_parameters, inits, times=times) %>% 
    mutate(multiplier=1)
  k <- 1
  for(i in seq_along(parameter1_values)) {
    for(j in seq_along(parameter2_values)) {
      parameters <- base_parameters
      parameters[parameter_names[1]] <- parameter1_values[i]
      parameters[parameter_names[2]] <- parameter1_values[j]
      sol <- simulate_model(parameters, inits, times=times) %>% 
        mutate(
          multiplier1=multipliers1[i],
          multiplier2=multipliers2[j],
          )
      if(k == 1)
        big_df <- sol
      else
        big_df <- big_df %>% bind_rows(sol)
      
      k <- k + 1
    }
  }
  big_df %>% 
    mutate(
      parameter1=parameter_names[1],
      parameter2=parameter_names[2])
}

single_double_sensitivity_2d <- function(parameter_names, multipliers1, multipliers2, fit, list_stan_datasets) {
  
  stan_data <- list_stan_datasets$stan_data
  times <- seq(0.01, 15, 0.01)
  df_single <- calculate_sensitivity_2d(
    parameter_names, multipliers1, multipliers2, fit, stan_data,
    times, t_refeed=100) %>% 
    mutate(type="single") %>% 
    select(-c(l, m))
  df_double <- calculate_sensitivity_2d(
    parameter_names, multipliers1, multipliers2, fit, stan_data,
    times, t_refeed=3) %>% 
    mutate(type="double") %>% 
    select(-c(l, m))
  df_both <- df_single %>% 
    bind_rows(df_double)
  
  m_double_feed_effect <- matrix(
    nrow = length(multipliers1),
    ncol = length(multipliers2),
    )

  for(i in seq_along(multipliers1)) {
    for(j in seq_along(multipliers2)) {
      df_tmp_s <- df_single %>% 
        filter(
          multiplier1==multipliers1[i],
          multiplier2==multipliers2[j],
          )
      h_max <- last(df_tmp_s$h)
      f <- approxfun(df_tmp_s$h, df_tmp_s$time)
      s_time <- f(h_max / 2)
      df_tmp_d <- df_double %>% 
        filter(
          multiplier1==multipliers1[i],
          multiplier2==multipliers2[j],
        )
      h_max <- last(df_tmp_d$h)
      f <- approxfun(df_tmp_d$h, df_tmp_d$time)
      d_time <- f(h_max / 2)
      m_double_feed_effect[i, j] <- s_time - d_time
    }
  }
  colnames(m_double_feed_effect) <- multipliers2
  m_double_feed_effect <- m_double_feed_effect %>% 
    as.data.frame() %>% 
    mutate(multiplier1=multipliers1)
  df_long <- m_double_feed_effect %>% 
    pivot_longer(as.character(multipliers2)) %>% 
    rename(multiplier2=name) %>% 
    mutate(multiplier2=as.numeric(multiplier2))
  
  colnames(df_long)[1:2] <- parameter_names
  
  df_long
}

plot_single_double_sensitivity_2d <- function(df_2d_sensitivity, parameter_names) {
  colnames(df_2d_sensitivity)[1:2] <- c("V1", "V2")
  max_val <- max(df_2d_sensitivity$value)
  df_2d_sensitivity <- df_2d_sensitivity %>% 
    mutate(value=if_else(value < 0, 0, value))
  
  breaks_lower <- seq(0, ceiling(max_val) - 1, length.out=10)
  breaks_upper <- seq(1, ceiling(max_val), length.out=10)
  breaks_mid <- 0.5 * (breaks_lower + breaks_upper)
  
  ggplot(df_2d_sensitivity, aes(x=V1, y=V2)) +
    geom_contour_filled(aes(z=value)) +
    xlab(TeX(paste0("$", "\\", parameter_names[1], "$"))) +
    ylab(TeX(paste0("$", "\\", parameter_names[2], "$"))) +
    geom_point(data=tibble(V1=1, V2=1), colour="red", size=3) +
    theme_bw() +
    scale_fill_viridis_d("Double vs\nsingle feeding\neffect, days",
                         guide = guide_bins(title.position = "right", reverse=TRUE))
}
