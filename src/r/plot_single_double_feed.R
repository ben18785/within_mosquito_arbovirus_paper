
plot_single_double_feed <- function(fit, list_stan_datasets) {
    
  pred_y_hat <- fit$par$pred_y_hat[2,,,,]
  both <- list_stan_datasets$dataset_denv
  df_multiple <- list_stan_datasets$dataset_binary_denv
  data_in <- list_stan_datasets$stan_data
  sum_df <- list_stan_datasets$dataset_denv_sum
  
  midgut_1 <- tibble(
      denv_titer=pred_y_hat[1,,2,2],
      dilutions=5,
      feed_type="double",
      tissue="midgut",
      day=data_in$g_t,
      type="simulated")
    midgut_2 <- tibble(
      denv_titer=pred_y_hat[2,,2,2],
      dilutions=5,
      feed_type="single",
      tissue="midgut",
      day=data_in$g_t,
      type="simulated")
    legs_1 <- tibble(
      denv_titer=pred_y_hat[1,,3,2],
      dilutions=5,
      feed_type="double",
      tissue="legs",
      day=data_in$g_t,
      type="simulated")
    legs_2 <- tibble(
      denv_titer=pred_y_hat[2,,3,2],
      dilutions=5,
      feed_type="single",
      tissue="legs",
      day=data_in$g_t,
      type="simulated")
    thresholds <- data_in$titer_lower_bound
    sigmas <- fit$par$sigma
    prob_infect_midgut <- logistic_curve(5 * (1 / fit$par$zeta), fit$par$b1, fit$par$b2, fit$par$b3, fit$par$b4) 
    prob_escape_midgut <- prob_infect_midgut * fit$par$phi_d
    
    lookup <- tribble(
      ~tissue,
      "midgut",
      "legs"
    ) %>% 
      mutate(threshold=thresholds,
             sigma=sigmas) %>% 
      mutate(prob_escape=c(prob_infect_midgut, prob_escape_midgut))
    
    sim_df <- midgut_1 %>% 
      bind_rows(midgut_2) %>% 
      bind_rows(legs_1) %>% 
      bind_rows(legs_2) %>% 
      left_join(lookup) %>% 
      mutate(denv_titer=if_else(denv_titer < 0, 1e-10, denv_titer)) %>% 
      mutate(middle=prob_escape * (1-plnorm(threshold, log(denv_titer), sigma))) %>% 
      rename(dilution=dilutions)
    temp_1 <- df_multiple %>% 
      rename(n_positive=n_infected) %>% 
      filter(dilutions==5) %>%  
      mutate(dilution=dilutions) %>% 
      select(type, tissue, dilution, day, n_positive, n_dissected) %>% 
      rename(feed_type=type) %>% 
      mutate(type="actual") %>% 
      mutate(middle=qbeta(0.5, 1 + n_positive, 1 + n_dissected - n_positive),
             lower=qbeta(0.025, 1 + n_positive, 1 + n_dissected - n_positive),
             upper=qbeta(0.975, 1 + n_positive, 1 + n_dissected - n_positive))
    all_df <- sim_df %>% 
      bind_rows(temp_1) %>% 
      mutate(tissue=as.factor(tissue)) %>% 
      mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
      filter(tissue=="legs") %>% 
      mutate(feed_type=as.factor(feed_type)) %>% 
      mutate(feed_type=fct_rev(feed_type))
      
    
    ggplot(all_df,
                aes(x=day, y=middle, colour=feed_type)) +
      geom_vline(xintercept = 3, linetype=2) +
      geom_line(data=all_df %>% filter(type=="simulated")) +
      geom_pointrange(data=all_df %>% filter(type!="simulated"),
                      aes(ymin=lower, ymax=upper),
                      position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2)) +
      scale_color_brewer("Feed type", palette = "Dark2") +
      xlab("DPI") +
      ylab("Positive") +
      scale_y_continuous(labels = scales::percent) +
      theme(legend.position = c(0.8, 0.4)) +
      scale_x_continuous(limits=c(0, 12.5))
}

plot_noninfectious_then_infectious_double <- function(fit, list_stan_datasets) {
  
}
