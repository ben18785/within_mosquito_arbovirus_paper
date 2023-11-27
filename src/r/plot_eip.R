

rhs <- function(t, state, parameters) {

  with(as.list(c(state, parameters)), {
    x_t = x_star + exp(-eta * t) * (x_0 - x_star)
    dl_dt = -gamma * l - x_t * (k_lm * l) / (a + l)
    dm_dt = x_t * (k_lm * l) / (a + l) + alpha_m * m * (1 - m / k_m) - k_mh * (x_t - x_star) * m
    dh_dt = k_mh * (x_t - x_star) * m + alpha_h * h * (1 - h/k_h)
    list(c(dl_dt, dm_dt, dh_dt))
  })
  
}

eip_quantile_calculator <- function(
    ode_parameters, sigma, inits, stan_data,
    times=seq(0, 15, 0.01), q=0.5) {
  
  sol <- ode(y=inits, func=rhs, parms=ode_parameters, times=times) %>% 
    as.data.frame()
  
  threshold <- stan_data$titer_lower_bound[2]
  prev <- (1 - plnorm(threshold, log(sol$h), sigma))
  prev_to_time <- approxfun(prev, times)
  
  prev_to_time(q)
}

eip_quantile_samples <- function(dose_multiplier, q, fit, stan_data, nsamples=200) {
  
  gamma <- rstan::extract(fit, "gamma")[[1]]
  k_lm <- rstan::extract(fit, "k_lm")[[1]]
  a <- rstan::extract(fit, "a")[[1]]
  eta <- rstan::extract(fit, "eta")[[1]]
  x_star <- rstan::extract(fit, "x_star")[[1]]
  alpha_m <- rstan::extract(fit, "alpha_m")[[1]]
  k_m <- rstan::extract(fit, "k_m")[[1]]
  k_mh <- rstan::extract(fit, "k_mh")[[1]]
  alpha_h <- rstan::extract(fit, "alpha_h")[[1]]
  k_h <- rstan::extract(fit, "k_h")[[1]]
  l_0 <- rstan::extract(fit, "l0")[[1]]
  
  eips <- vector(length = nsamples)
  for(j in seq_along(eips)) {
    
    i <- sample(length(gamma), size=1)
    
    parameters <- c(
      gamma=gamma[i],
      k_lm=k_lm[i],
      a=a[i],
      eta=eta[i],
      x_star=x_star[i],
      alpha_m=alpha_m[i],
      k_m=k_m[i],
      k_mh=k_mh[i],
      alpha_h=alpha_h[i],
      k_h=k_h[i],
      x_0=stan_data$x_0
    )
    
    sigmas <- rstan::extract(fit, "sigma")[[1]]
    sigmas <- sigmas[i, ]
    
    inits <- c(
      l=l_0[i] * dose_multiplier,
      m=0,
      h=0
    )

    eips[j] <- eip_quantile_calculator(
      parameters, sigmas[2], inits, stan_data,
      q=q)
  }
  eips
}


eip_dose_response <- function(fit, stan_data, nsamples_val) {
  
  doses <- c(seq(0.05, 2, 0.1), 2)
  quantiles <- c(0.25, 0.5, 0.75)
  k <- 1
  for(j in seq_along(quantiles)) {
    print(paste0("quantile = ", quantiles[j]))
    for(i in seq_along(doses)) {
      print(paste0("dose = ", doses[i]))
      eips <- eip_quantile_samples(doses[i], q=quantiles[j], fit, stan_data,
                                   nsamples=nsamples_val)
      tmp <- tibble(eip=eips, dose=doses[i], q=quantiles[j]) %>% 
        mutate(iteration=seq_along(eip))
      if(k == 1)
        big_df <- tmp
      else
        big_df <- big_df %>% bind_rows(tmp)
      k <- k + 1
    }
  }
  big_df
}

plot_eip_dose_response <- function(df_eip_dose_respose) {
  
  df_eip_dose_respose %>% 
    mutate(q = q * 100) %>% 
    mutate(q = as.character(q)) %>% 
    mutate(q = paste0(q, "%")) %>% 
    group_by(dose, q) %>% 
    summarise(lower=quantile(eip, 0.25, na.rm=T),
              middle=median(eip, na.rm=T),
              upper=quantile(eip, 0.75, na.rm=T)) %>% 
    ggplot(aes(x=dose, y=middle)) +
    geom_line(aes(colour=as.factor(q))) +
    geom_ribbon(aes(ymin=lower, ymax=upper,
                    fill=as.factor(q)),
                alpha=0.3) +
    scale_colour_brewer("EIP quant.",
                      palette = "Dark2") +
    scale_fill_brewer("EIP quant.",
                       palette = "Dark2") +
    xlab("Concentration") +
    ylab("Time for dissemination, DPI") +
    theme_bw() +
    geom_vline(xintercept = 0.08, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    theme(
      legend.position = c(0.75, 0.7)
    )
}

