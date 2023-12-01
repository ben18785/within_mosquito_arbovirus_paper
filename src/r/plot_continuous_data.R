
plot_continuous_data <- function(sampling_fit, list_stan_datasets) {
  
  df <- list_stan_datasets$dataset_denv_expanded
  data_in <- list_stan_datasets$stan_data
  
  pred_y_hat <- rstan::extract(sampling_fit, "pred_y_hat")[[1]]
  pred_y_hat <- pred_y_hat[, 1, 2, , , ]
  lowers <- apply(pred_y_hat, c(2, 3, 4), function(x) quantile(x, 0.025))
  middles <- apply(pred_y_hat, c(2, 3, 4), function(x) quantile(x, 0.5))
  uppers <- apply(pred_y_hat, c(2, 3, 4), function(x) quantile(x, 0.975))
  
  midgut_1 <- tibble(
    middle=middles[, 2, 1],
    lower=lowers[, 2, 1],
    upper=uppers[, 2, 1],
    dilutions=1,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_2 <- tibble(
    middle=middles[, 2, 2],
    lower=lowers[, 2, 2],
    upper=uppers[, 2, 2],
    dilutions=5,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_3 <- tibble(
    middle=middles[, 2, 3],
    lower=lowers[, 2, 3],
    upper=uppers[, 2, 3],
    dilutions=12,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_df <- midgut_1 %>% 
    bind_rows(midgut_2) %>% 
    bind_rows(midgut_3)
  legs_1 <- tibble(
    middle=middles[, 3, 1],
    lower=lowers[, 3, 1],
    upper=uppers[, 3, 1],
    dilutions=1,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_2 <- tibble(
    middle=middles[, 3, 2],
    lower=lowers[, 3, 2],
    upper=uppers[, 3, 2],
    dilutions=5,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_3 <- tibble(
    middle=middles[, 3, 3],
    lower=lowers[, 3, 3],
    upper=uppers[, 3, 3],
    dilutions=12,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_df <- legs_1 %>% 
    bind_rows(legs_2) %>% 
    bind_rows(legs_3)
  
  sim_df <- midgut_df %>% 
    bind_rows(legs_df)
  
  # calculate titers on raw scale
  df_scales <- list_stan_datasets$dataset_denv_sum
  sim_df <- sim_df %>% 
    left_join(df_scales)
  sim_df <- sim_df %>% 
    mutate(
      lower=lower * overall_denv_titer,
      middle=middle * overall_denv_titer,
      upper=upper * overall_denv_titer,
    )
  sim_df <- sim_df %>% 
    filter(day <= 10)
  
  df <- df %>% 
    rename(middle=denv_titer) %>% 
    mutate(middle=if_else(middle < 0, 0, middle)) %>% 
    mutate(middle=middle * overall_denv_titer)
  
  sim_df <- sim_df %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
    mutate(concentration=1/dilutions) %>% 
    mutate(concentration=format(round(concentration, 2), nsmall = 2)) 
  df <- df %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
    mutate(concentration=1/dilutions) %>% 
    mutate(concentration=format(round(concentration, 2), nsmall = 2))
  
  df_sum <- df %>% 
    group_by(day, concentration, tissue) %>% 
    summarise(
      lower=quantile(middle, 0.1),
      upper=quantile(middle, 0.9),
      middle=quantile(middle, 0.5),
    )
  
  ggplot(df, aes(x=day, y=middle)) +
    geom_line(data=sim_df, colour="blue") +
    geom_ribbon(data=sim_df, aes(ymin=lower, ymax=upper),
                alpha=0.4, fill="blue") +
    geom_pointrange(data=df_sum, aes(ymin=lower, ymax=upper)) +
    facet_grid(vars(tissue), vars(concentration),
               scales="free") +
    ylab("DENV titer, GE/uL") +
    xlab("Days post infection") +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) + 
    scale_y_sqrt(labels = function(x) format(x, scientific = TRUE)) +
    theme_bw() +
    theme(text=element_text(size=14))
}
