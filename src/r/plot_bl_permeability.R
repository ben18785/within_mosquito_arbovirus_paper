
plot_bl_permeability <- function(fit, list_stan_datasets) {
  
  eta <- fit$par$eta
  x_star <- fit$par$x_star
  data_in <- list_stan_datasets$stan_data
  x_0 <- data_in$x_0
  
  bl_permeability <- function(t, eta, x_0, x_star, t_refeed) {
    if(t < t_refeed)
      x_star + exp(-eta*t) * (x_0 - x_star)
    else
      x_star + exp(-eta*(t - t_refeed)) * (x_0 - x_star)
  }
  
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