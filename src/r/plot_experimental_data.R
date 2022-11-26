
plot_experimental_data <- function(df_midgut_legs) {
  
  df_midgut_legs %>% 
    ggplot(aes(x=day, y=denv_titer, colour=sample)) +
    geom_jitter(width=0.3) +
    scale_y_log10() +
    geom_smooth(se=FALSE) +
    ylab("DENV titer") +
    scale_color_brewer(palette = "Spectral") +
    facet_grid(~tissue)
}