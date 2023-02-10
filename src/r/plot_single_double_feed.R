
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
  
  data_in <- list_stan_datasets$stan_data
  preds <- fit$par$pred_y_hat_other_order
  
  # get simulations where non-infectious first
  # remove last term in series since this is repeated in the first term of next
  g_t_before <- data_in$g_t_before[-length(data_in$g_t_before)]
  preds <- preds[-length(g_t_before), ]
  g_t_after <- data_in$g_t_after + max(g_t_before)
  
  times <- c(g_t_before, g_t_after)
  df_double_wrong <- tibble(time=times, midgut=preds[, 2], legs=preds[, 3]) %>% 
    pivot_longer(-time) %>% 
    mutate(type="double:\nnoninfectious first")
  
  # get double where infectious first
  preds <- fit$par$pred_y_hat[2, 1, , , 2]
  df_double <- tibble(time=times, midgut=preds[, 2], legs=preds[, 3]) %>% 
    pivot_longer(-time) %>% 
    mutate(type="double:\ninfectious first")
  
  # get single feeds
  preds <- fit$par$pred_y_hat[2, 2, , , 2]
  df_single <- tibble(time=times, midgut=preds[, 2], legs=preds[, 3]) %>% 
    pivot_longer(-time) %>% 
    mutate(type="single")
  
  df <- df_double_wrong %>% 
    bind_rows(df_double,
              df_single) %>% 
    mutate(type=as.factor(type)) %>% 
    mutate(type=fct_relevel(
      type,
      "single",
      "double:\ninfectious first",
      "double:\nnoninfectious first")) %>% 
    rename(tissue=name)
  
  thresholds <- data_in$titer_lower_bound
  prob_infect_midgut <- logistic_curve(5 * (1 / fit$par$zeta), fit$par$b1, fit$par$b2, fit$par$b3, fit$par$b4) 
  prob_escape_midgut <- prob_infect_midgut * fit$par$phi_d
  sigmas <- fit$par$sigma

  lookup <- tribble(
    ~tissue,
    "midgut",
    "legs"
  ) %>% 
    mutate(threshold=thresholds,
           sigma=sigmas) %>% 
    mutate(prob_escape=c(prob_infect_midgut, prob_escape_midgut))
  
  df_all <- df %>% 
    left_join(lookup) %>% 
    mutate(middle=prob_escape * (1-plnorm(threshold, log(value), sigma))) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs"))
  
  df_all %>% 
    ggplot(aes(x=time, y=middle, colour=type)) +
    geom_vline(xintercept = 3, linetype=2) +
    geom_line() +
    facet_wrap(~tissue) +
    xlab("DPI") +
    scale_color_brewer("Feed type", palette = "Dark2") +
    ylab("Positive") +
    scale_y_continuous(labels = scales::percent,
                       limits=c(0, 1)) +
    theme_bw() +
    scale_x_continuous(limits=c(0, 15))
}

plot_single_double_feed_mcmc <- function(
    fit, list_stan_datasets,
    quantiles=c(0.025, 0.975)) {
  
  pred_y_hat <- rstan::extract(fit, "pred_y_hat")[[1]]
  pred_y_hat <- pred_y_hat[, 2, , , , ]
  middles <- apply(pred_y_hat, c(2, 3, 4, 5), function(x) quantile(x, 0.5))
  lowers <- apply(pred_y_hat, c(2, 3, 4, 5), function(x) quantile(x, quantiles[1]))
  uppers <- apply(pred_y_hat, c(2, 3, 4, 5), function(x) quantile(x, quantiles[2]))
  
  both <- list_stan_datasets$dataset_denv
  df_multiple <- list_stan_datasets$dataset_binary_denv
  data_in <- list_stan_datasets$stan_data
  sum_df <- list_stan_datasets$dataset_denv_sum
  
  midgut_1 <- tibble(
    middle=middles[1,,2,2],
    upper=uppers[1,,2,2],
    lower=lowers[1,,2,2],
    dilutions=5,
    feed_type="double",
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_2 <- tibble(
    middle=middles[2,,2,2],
    upper=uppers[2,,2,2],
    lower=lowers[2,,2,2],
    dilutions=5,
    feed_type="single",
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  legs_1 <- tibble(
    middle=middles[1,,3,2],
    upper=uppers[1,,3,2],
    lower=lowers[1,,3,2],
    dilutions=5,
    feed_type="double",
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_2 <- tibble(
    middle=middles[2,,3,2],
    upper=uppers[2,,3,2],
    lower=lowers[2,,3,2],
    dilutions=5,
    feed_type="single",
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  thresholds <- data_in$titer_lower_bound
  
  # use median values for sigmas, b1-b4, zetam and phi_d
  sigmas <- apply(rstan::extract(fit, "sigma")[[1]], 2, median)
  b1 <- median(rstan::extract(fit, "b1")[[1]])
  b2 <- median(rstan::extract(fit, "b2")[[1]])
  b3 <- median(rstan::extract(fit, "b3")[[1]])
  b4 <- median(rstan::extract(fit, "b4")[[1]])
  phi_d <- median(rstan::extract(fit, "phi_d")[[1]])
  zeta <- median(rstan::extract(fit, "zeta")[[1]])
  prob_infect_midgut <- logistic_curve(5 * (1 / zeta), b1, b2, b3, b4) 
  prob_escape_midgut <- prob_infect_midgut * phi_d
  
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
    mutate(middle=if_else(middle < 0, 1e-10, middle),
           lower=if_else(lower < 0, 1e-10, lower),
           upper=if_else(upper < 0, 1e-10, upper)) %>%
    mutate(
      middle=prob_escape * (1-plnorm(threshold, log(middle), sigma)),
      lower=prob_escape * (1-plnorm(threshold, log(lower), sigma)),
      upper=prob_escape * (1-plnorm(threshold, log(upper), sigma))
      ) %>% 
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
         aes(x=day, y=middle)) +
    geom_vline(xintercept = 3, linetype=2) +
    geom_ribbon(data=all_df %>% filter(type=="simulated"),
                aes(ymin=lower, ymax=upper, fill=feed_type),
                alpha=0.3) +
    geom_line(data=all_df %>% filter(type=="simulated"),
              aes(colour=feed_type)) +
    geom_pointrange(data=all_df %>% filter(type!="simulated"),
                    aes(ymin=lower, ymax=upper, colour=feed_type),
                    position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2)) +
    scale_color_brewer("Feed type", palette = "Dark2") +
    scale_fill_brewer("Feed type", palette = "Dark2") +
    xlab("DPI") +
    ylab("Positive") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.4)) +
    scale_x_continuous(limits=c(0, 12.5))
}
