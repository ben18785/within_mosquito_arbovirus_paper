
prepare_simulated_data <- function(first_index, pred_y_hat, data_in) {
  midgut_1 <- tibble(
    DENV=pred_y_hat[first_index, , 2, 1],
    dilutions=1,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_2 <- tibble(
    DENV=pred_y_hat[first_index, , 2, 2],
    dilutions=5,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_3 <- tibble(
    DENV=pred_y_hat[first_index, , 2, 3],
    dilutions=12,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_4 <- tibble(
    DENV=pred_y_hat[first_index, , 2, 4],
    dilutions=20,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_5 <- tibble(
    DENV=pred_y_hat[first_index, , 2, 5],
    dilutions=25,
    tissue="midgut",
    day=data_in$g_t,
    type="simulated")
  midgut_df <- midgut_1 %>% 
    bind_rows(midgut_2) %>% 
    bind_rows(midgut_3) %>% 
    bind_rows(midgut_4) %>% 
    bind_rows(midgut_5)
  
  legs_1 <- tibble(
    DENV=pred_y_hat[first_index, , 3, 1],
    dilutions=1,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_2 <- tibble(
    DENV=pred_y_hat[first_index, , 3, 2],
    dilutions=5,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_3 <- tibble(
    DENV=pred_y_hat[first_index, , 3, 3],
    dilutions=12,
    tissue="legs",
    day=data_in$g_t,
    type="simulated")
  legs_df <- legs_1 %>% 
    bind_rows(legs_2) %>% 
    bind_rows(legs_3)
  
  sim_df <- midgut_df %>% 
    bind_rows(legs_df) %>% 
    pivot_longer(c(DENV)) %>% 
    rename(data_type=name) %>% 
    rename(denv_titer=value)
  
  sim_df
}

dose_response_curve <- function(dilutions, fit) {
  phi_m <- vector(length = length(dilutions))
  for(i in seq_along(dilutions)) {
    phi_m[i] <- logistic_curve(
      dilutions[i], fit$par$b1, fit$par$b2, fit$par$b3, fit$par$b4)
  }
  phi_m
}

prepare_data_prefit <- function(fit, list_stan_datasets) {
  
  pred_y_hat <- fit$par$pred_y_hat[,2,,,]
  both <- list_stan_datasets$dataset_denv
  df_multiple <- list_stan_datasets$dataset_binary_denv
  data_in <- list_stan_datasets$stan_data
  sum_df <- list_stan_datasets$dataset_denv_sum
  
  sim_df_c <- prepare_simulated_data(1, pred_y_hat, data_in) %>% 
    mutate(category="continuous")
  sim_df_d <- prepare_simulated_data(2, pred_y_hat, data_in) %>% 
    mutate(category="dichtonomous")
  sim_df <- sim_df_c %>% 
    bind_rows(sim_df_d)
  
  real_df <- both %>% 
    mutate(type="real")
  combined_df <- sim_df %>% 
    bind_rows(real_df) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
    mutate(denv_above_zero=if_else(denv_titer > 0, "non-censored", "below detection\nthreshold"))
  
  temp <- combined_df %>% 
    filter(type != "simulated") %>% 
    group_by(tissue, dilutions, day) %>% 
    summarise(n_positive=sum(denv_above_zero=="non-censored", na.rm = TRUE),
              n_dissected=n()) %>% 
    mutate(type="single") %>% 
    rename(dilution=dilutions) %>% 
    mutate(data_type="DENV") %>% 
    mutate(category="continuous")
  
  # only plot single feeds here
  temp_1 <- df_multiple %>% 
    rename(n_positive=n_infected) %>% 
    mutate(dilution=dilutions) %>% 
    select(type, tissue, dilution, day, n_positive, n_dissected, virus) %>% 
    rename(data_type=virus) %>% 
    filter(type=="single") %>% 
    mutate(category="dichtonomous")
  
  count_combined <- temp_1 %>% 
    bind_rows(temp) %>% 
    mutate(middle=qbeta(0.5, 1 + n_positive, 1 + n_dissected - n_positive),
           lower=qbeta(0.025, 1 + n_positive, 1 + n_dissected - n_positive),
           upper=qbeta(0.975, 1 + n_positive, 1 + n_dissected - n_positive))
  sum_df <- sum_df %>% 
    mutate(threshold=data_in$titer_lower_bound)
  
  # use mean sigma and phis
  sigma <- tibble(sigma=fit$par$sigma,
                  tissue=c("midgut", "legs"))
  dilutions <- c(1, 5, 12, 20, 25)
  
  phi_m <- dose_response_curve(dilutions * (1 / fit$par$zeta), fit)
  phi_1 <- tribble(~phi, ~tissue, ~dilutions,
                   phi_m[1], "midgut", 1,
                   phi_m[2], "midgut", 5,
                   phi_m[3], "midgut", 12,
                   phi_m[4], "midgut", 20,
                   phi_m[5], "midgut", 25)
  phi_legs <- tribble(~phi, ~tissue, ~dilutions,
                      fit$par$phi_d, "legs", 1,
                      fit$par$phi_d, "legs", 5,
                      fit$par$phi_d, "legs", 12)
  phi_1 <- phi_1 %>%
    bind_rows(phi_legs) %>% 
    mutate(category="dichtonomous")
  
  phi_m <- dose_response_curve(dilutions, fit)
  phi_2 <- tribble(~phi, ~tissue, ~dilutions,
                   phi_m[1], "midgut", 1,
                   phi_m[2], "midgut", 5,
                   phi_m[3], "midgut", 12,
                   phi_m[4], "midgut", 20,
                   phi_m[5], "midgut", 25)
  phi_legs <- tribble(~phi, ~tissue, ~dilutions,
                      fit$par$phi_d, "legs", 1,
                      fit$par$phi_d, "legs", 5,
                      fit$par$phi_d, "legs", 12)
  phi_2 <- phi_2 %>%
    bind_rows(phi_legs) %>% 
    mutate(category="continuous")
  phi <- phi_1 %>% 
    bind_rows(phi_2)
  simulated_df <- combined_df %>% 
    filter(type=="simulated") %>% 
    left_join(sum_df, by="tissue") %>% 
    left_join(sigma, by="tissue") %>%
    left_join(phi, by=c("tissue", "dilutions", "category")) %>%
    mutate(denv_titer=if_else(denv_titer < 0, 1e-10, denv_titer)) %>% 
    mutate(middle= phi * (1-plnorm(threshold, log(denv_titer), sigma))) %>% 
    rename(dilution=dilutions)
  
  all_df <- count_combined %>% 
    bind_rows(simulated_df) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs"))
  
  all_df
}

plot_fit_prevalence_midgut_legs <- function(fit, list_stan_datasets) {
  
  all_df <- prepare_data_prefit(fit, list_stan_datasets) %>% 
    mutate(concentration=1/dilution) %>% 
    mutate(concentration=format(round(concentration, 2), nsmall = 2)) 
  
  g <- ggplot(all_df,
              aes(x=day, y=middle)) +
    geom_line(data=all_df %>% filter(type=="simulated"),
              aes(colour=category)) +
    geom_pointrange(data=all_df %>% filter(type!="simulated"),
                    aes(ymin=lower, ymax=upper, colour=category)) +
    xlab("DPI") +
    ylab("Positive") +
    scale_y_continuous(labels = scales::percent) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(vars(tissue),
               vars(concentration))
  g
}


plot_fit_prevalence_dose_response <- function(fit, list_stan_datasets) {
  
  data_in <- list_stan_datasets$stan_data
  all_df <- prepare_data_prefit(fit, list_stan_datasets) %>% 
    filter(tissue=="midgut") %>% 
    filter(day == 5) %>% 
    filter(type != "simulated") %>% 
    mutate(concentration=1/dilution)
  
  # calculate dose response curve
  phi_1 <- dose_response_curve(data_in$dilutions_sim_fine, fit)
  phi_2 <- dose_response_curve(data_in$dilutions_sim_fine * (1 / fit$par$zeta), fit)
  sigma <- fit$par$sigma[1]
  threshold <- data_in$titer_lower_bound[1]
  prob <- (1 - plnorm(threshold, log(fit$par$dose_response_midgut), sigma))
  prob[1, ] = prob[1, ] * phi_1
  prob[2, ] = prob[2, ] * phi_2
  df_probs_c <- tibble(
    dilution=data_in$dilutions_sim_fine,
    middle=prob[1, ],
    category="continuous")
  df_probs_d <- tibble(
    dilution=data_in$dilutions_sim_fine,
    middle=prob[2, ],
    category="dichtonomous")
  df_probs <- df_probs_c %>% 
    bind_rows(df_probs_d) %>% 
    mutate(concentration=1/dilution)
  
  g <- ggplot(df_probs %>% filter(category=="dichtonomous"),
              aes(x=concentration, y=middle, group=category)) +
    geom_line(colour="grey") +
    geom_pointrange(data=all_df%>% filter(category=="dichtonomous"),
                    aes(ymin=lower, ymax=upper),
                    colour="black") +
    xlab("Concentration") +
    ylab("Positive") +
    scale_y_continuous(labels = scales::percent,
                       limits=c(0, 1)) +
    scale_x_log10(limits=c(0.02, 2))
  g
}

plot_fit_prevalence_midgut_only <- function(fit, list_stan_datasets) {
  
  all_df <- prepare_data_prefit(fit, list_stan_datasets) %>% 
    filter(tissue == "midgut") %>% 
    mutate(concentration=1/dilution) %>% 
    mutate(concentration=format(round(concentration, 2), nsmall = 2)) %>% 
    filter(category == "dichtonomous") %>% 
    mutate(concentration=fct_rev(concentration))
  
  g <- ggplot(all_df,
              aes(x=day, y=middle, colour=concentration,
                  group=as.factor(concentration))) +
    geom_line(data=all_df %>% filter(type=="simulated")) +
    geom_pointrange(data=all_df %>% filter(type!="simulated"),
                    aes(ymin=lower, ymax=upper),
                    position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2)) +
    xlab("DPI") +
    ylab("Positive") +
    scale_y_continuous(labels = scales::percent,
                       limits=c(0, 1)) +
    scale_color_brewer("Concentration",
                       palette = "Spectral") +
    scale_x_continuous(limits=c(0, 13)) +
    theme(legend.position = c(0.8, 0.4))
  g
}

plot_midgut_dose_response_combined <- function(fit, list_stan_datasets) {
  g1 <- plot_fit_prevalence_midgut_only(fit, list_stan_datasets)
  g2 <- plot_fit_prevalence_dose_response(fit, list_stan_datasets)
  g <- plot_grid(g1, g2, nrow = 1,
                 labels = c("A.", "B."),
                 label_x = -0.01)
  g
}

plot_fit_prevalence_legs <- function(fit, list_stan_datasets) {
  
  all_df <- prepare_data_prefit(fit, list_stan_datasets) %>% 
    mutate(concentration=1/dilution) %>% 
    mutate(concentration=format(round(concentration, 2), nsmall = 2)) %>% 
    filter(tissue=="legs") %>% 
    mutate(concentration=as.factor(concentration)) %>% 
    mutate(concentration=fct_rev(concentration))
  
  g <- ggplot(all_df,
              aes(x=day, y=middle, colour=concentration)) +
    geom_line(data=all_df %>% filter(type=="simulated", category=="continuous")) +
    geom_pointrange(data=all_df %>% filter(type!="simulated"),
                    aes(ymin=lower, ymax=upper, shape=category),
                    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)) +
    xlab("DPI") +
    ylab("Positive") +
    scale_shape(guide="none") +
    scale_y_continuous(labels = scales::percent) +
    scale_color_brewer("Concentration", palette = "Spectral") +
    scale_x_continuous(limits = c(0, 15))
  g
}
