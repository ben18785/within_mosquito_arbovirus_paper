bl_permeability <- function(t, eta, x_0, x_star, t_refeed) {
  if(t < t_refeed)
    x_star + exp(-eta*t) * (x_0 - x_star)
  else
    x_star + exp(-eta*(t - t_refeed)) * (x_0 - x_star)
}

plot_bl_permeability <- function(fit, list_stan_datasets) {
  
  eta <- fit$par$eta
  x_star <- fit$par$x_star
  data_in <- list_stan_datasets$stan_data
  x_0 <- data_in$x_0
  
  times <- seq(0, 12.5, 0.1)
  single <- map_dbl(times, ~bl_permeability(., eta, x_0, x_star, 100))
  double <- map_dbl(times, ~bl_permeability(., eta, x_0, x_star, 3))
  
  df <- tibble(time=times, value=single, type="single") %>% 
    bind_rows(tibble(time=times, value=double, type="double") %>% 
                filter(time > 2.9)) %>% 
    mutate(type=as.factor(type)) %>% 
    mutate(type=fct_rev(type))
  
  ggplot(df, aes(x=time, y=value, colour=type)) +
    geom_line() +
    scale_color_brewer("Feed type", palette = "Dark2") +
    xlab("DPI") +
    ylab("Basal lamina permeability") +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          axis.ticks.y=element_blank())
}

plot_bl_permeability_mcmc <- function(fit, list_stan_datasets,
                                      include_uncertainty=FALSE) {
  
  eta <- rstan::extract(fit, "eta")[[1]]
  x_star <- rstan::extract(fit, "x_star")[[1]]
  
  data_in <- list_stan_datasets$stan_data
  x_0 <- data_in$x_0
  times <- seq(0, 12.5, 0.1)
  
  if(include_uncertainty) {
    n_iterations <- length(eta)
    
    for(i in 1:n_iterations) {
      single <- map_dbl(times, ~bl_permeability(., eta[i], x_0, x_star[i], 100))
      double <- map_dbl(times, ~bl_permeability(., eta[i], x_0, x_star[i], 3))
      tmp_single <- tibble(time=times, value=single, type="single")
      tmp_double <- tibble(time=times, value=double, type="double") %>% 
        filter(time > 2.9)
      tmp_both <- tmp_single %>%
        bind_rows(tmp_double) %>% 
        mutate(iteration=i)
      if(i == 1)
        big_df <- tmp_both
      else
        big_df <- big_df %>% bind_rows(tmp_both)
    }
    df <- big_df %>% 
      mutate(type=as.factor(type)) %>% 
      mutate(type=fct_rev(type)) %>% 
      group_by(time, type) %>% 
      summarise(middle=quantile(value, 0.5),
                lower=quantile(value, 0.025),
                upper=quantile(value, 0.975))
    
    g  <- ggplot(df, aes(x=time, y=middle)) +
      geom_ribbon(aes(ymin=lower, ymax=upper,
                      fill=type), alpha=0.1) +
      geom_line(aes(colour=type)) +
      scale_color_brewer("Feed type", palette = "Dark2") +
      scale_fill_brewer("Feed type", palette = "Dark2") +
      xlab("DPI") +
      ylab("Basal lamina permeability") +
      theme(axis.text.y = element_blank(),
            legend.position = "none",
            axis.ticks.y=element_blank())
  } else {
    times <- seq(0, 12.5, 0.1)
    eta <- median(eta)
    x_star <- median(x_star)
    single <- map_dbl(times, ~bl_permeability(., eta, x_0, x_star, 100))
    double <- map_dbl(times, ~bl_permeability(., eta, x_0, x_star, 3))
    
    df <- tibble(time=times, value=single, type="single") %>% 
      bind_rows(tibble(time=times, value=double, type="double") %>% 
                  filter(time > 2.9)) %>% 
      mutate(type=as.factor(type)) %>% 
      mutate(type=fct_rev(type))
    
    ggplot(df, aes(x=time, y=value, colour=type)) +
      geom_line() +
      scale_color_brewer("Feed type", palette = "Dark2") +
      xlab("DPI") +
      ylab("Basal lamina permeability") +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            legend.position = "none",
            axis.ticks.y=element_blank())
    }
}