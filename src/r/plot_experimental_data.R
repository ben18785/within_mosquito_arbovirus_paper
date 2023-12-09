
plot_experimental_data <- function(df_midgut_legs) {
  
  df_midgut_legs %>% 
    ggplot(aes(x=day, y=denv_titer, colour=sample)) +
    geom_jitter(width=0.3) +
    scale_y_log10() +
    geom_smooth(se=FALSE) +
    ylab("DENV titer") +
    scale_color_brewer(palette = "Spectral") +
    facet_grid(~tissue) +
    theme(text=element_text(size=14))
}

plot_experimental_data_midgut <- function(df_midgut_legs) {
  
  df_midgut_legs %>% 
    mutate(Concentration = 1 / dilutions) %>% 
    mutate(Concentration=format(round(Concentration, 2), nsmall = 2)) %>%
    mutate(Concentration=as.factor(Concentration)) %>% 
    mutate(Concentration=fct_rev(Concentration)) %>% 
    filter(tissue=="midgut") %>% 
    ggplot(aes(x=day, y=denv_titer, colour=Concentration)) +
    geom_jitter(width=0.2, aes(shape=as.factor(Concentration))) +
    scale_y_log10() +
    geom_smooth(se=FALSE) +
    scale_shape("Concentration") +
    ylab("DENV titer, GE/tissue") +
    xlab("Days post infection") +
    scale_color_viridis_d("Concentration") +
    theme_bw() +
    theme(text=element_text(size=14))
}