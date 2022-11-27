
plot_chp_damage_fit <- function(fit, df_chp_damage) {
  eta <- fit$par$eta
  chp_vals <- fit$par$chp_vals
  chp_inf <- chp_vals[1]
  chp_0 <- chp_vals[2]
  
  exp_decline <- function(t, chp_0, chp_inf, eta) {
    chp_inf + (chp_0 - chp_inf) * exp(- eta * t)
  }
  
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