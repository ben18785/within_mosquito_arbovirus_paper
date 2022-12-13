
exp_decline <- function(t, chp_0, chp_inf, eta) {
  chp_inf + (chp_0 - chp_inf) * exp(- eta * t)
}

plot_chp_damage_fit <- function(fit, df_chp_damage) {
  eta <- fit$par$eta
  chp_vals <- fit$par$chp_vals
  chp_inf <- chp_vals[1]
  chp_0 <- chp_vals[2]
  
  times <- seq(0, 12.5, 0.1)
  chp <- map_dbl(times, ~exp_decline(., chp_0, chp_inf, eta))
  df_sim <- tibble(
    time=times,
    chp=chp,
    type="simulated"
  )
  df_chp_damage <- df_chp_damage %>% 
    mutate(type="real") %>% 
    mutate(assumption=if_else(name=="Unfed", 1, 0))
  
  ggplot(df_chp_damage, aes(x=time, y=chp)) +
    geom_jitter(width = 0.2, aes(shape=as.factor(assumption))) +
    geom_line(data=df_sim, aes(colour=type)) +
    xlab("DPI") +
    ylab("CHP") +
    scale_x_continuous(limits=c(0, 12.5)) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none")
}


plot_chp_damage_fit_mcmc <- function(fit, df_chp_damage) {
  
  eta <- rstan::extract(fit, "eta")[[1]]
  chp_sigma <- rstan::extract(fit, "chp_sigma")[[1]]
  chp_vals <- rstan::extract(fit, "chp_vals")[[1]]
  chp_inf <- chp_vals[, 1]
  chp_0 <- chp_vals[, 2]
  
  total_iterations <- length(eta)
  n_iterations <- 100
  n_rnorm_iterations <- 100
  times <- seq(0, 12.5, 0.1)
  k <- 1
  for(i in 1:n_iterations) {
    print(i)
    idx1 <- sample(total_iterations, 1)
    
    for(j in 1:n_rnorm_iterations) {
      
      chp <- map_dbl(times, ~exp_decline(
        ., chp_0[idx1], chp_inf[idx1], eta[idx1]))
      
      idx <- sample(n_iterations, 1)
      chp_sigma_tmp <- chp_sigma[idx]
      chp <- rlnorm(length(chp), log(chp), chp_sigma_tmp)
      tmp_df <- tibble(
        time=times,
        chp=chp,
        iteration=k
      )
      
      if(k == 1)
        big_df <- tmp_df
      else
        big_df <- big_df %>% bind_rows(tmp_df)
      
      k <- k + 1
    }
    print(nrow(big_df))

  }
  
  df_sim <- big_df %>% 
    group_by(time) %>% 
    summarise(middle=median(chp),
              lower=quantile(chp, 0.025),
              upper=quantile(chp, 0.975)) %>% 
    mutate(type="simulated")

  df_chp_damage <- df_chp_damage %>% 
    rename(middle=chp) %>% 
    mutate(type="real") %>% 
    mutate(assumption=if_else(name=="Unfed", 1, 0))
  
  ggplot(df_chp_damage, aes(x=time, y=middle)) +
    geom_jitter(width = 0.2,
                aes(shape=as.factor(assumption))) +
    geom_ribbon(data=df_sim,
                aes(ymin=lower, ymax=upper, fill=type),
                alpha=0.3) +
    geom_line(data=df_sim, aes(colour=type)) +
    xlab("DPI") +
    ylab("CHP") +
    scale_x_continuous(limits=c(0, 12.5)) +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "none")
}
