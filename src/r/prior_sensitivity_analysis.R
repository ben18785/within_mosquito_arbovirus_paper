
prior_sensitivity <- function(fit) {
  
  sensitivity_analysis <- function(old_prior_formula, fit, multiplier) {
    forms <- Reduce(paste, deparse(old_prior_formula))
    varbl_end <- unname(str_locate(forms, "~")[1, 1]) - 2
    varbl <<- substr(forms, 1, varbl_end) # for some reason varbl is only found it its a global
    new_form <- paste0(substr(forms, 1, nchar(forms) - 1), "* multiplier)")
    df <- make_spec(as.formula(new_form), multiplier=multiplier) %>% 
      adjust_weights(fit, keep_bad = TRUE) %>% # keep bad since otherwise summary functions throw error if pareto k > 0.7
      summarize(
        diff_l0=wasserstein(l0), mu_l0=mean(l0), sigma_l0=sd(l0),
        diff_zeta=wasserstein(zeta), mu_zeta=mean(zeta), sigma_zeta=sd(zeta),
        diff_gamma=wasserstein(gamma), mu_gamma=mean(gamma), sigma_gamma=sd(gamma),
        diff_a=wasserstein(a), mu_a=mean(a), sigma_a=sd(a),
        diff_alpha_m=wasserstein(alpha_m), mu_alpha_m=mean(alpha_m), sigma_alpha_m=sd(alpha_m),
        diff_k_m=wasserstein(k_m), mu_k_m=mean(k_m), sigma_k_m=sd(k_m),
        diff_k_lm=wasserstein(k_lm), mu_k_lm=mean(k_lm), sigma_k_lm=sd(k_lm),
        diff_k_mh=wasserstein(k_mh), mu_k_mh=mean(k_mh), sigma_k_mh=sd(k_mh),
        diff_alpha_h=wasserstein(alpha_h), mu_alpha_h=mean(alpha_h), sigma_alpha_h=sd(alpha_h),
        diff_k_h=wasserstein(k_h), mu_k_h=mean(k_h), sigma_k_h=sd(k_h),
        diff_eta=wasserstein(eta), mu_eta=mean(eta), sigma_eta=sd(eta),
        diff_x_star=wasserstein(x_star), mu_x_star=mean(x_star), sigma_x_star=sd(x_star),
        diff_b1=wasserstein(b1), mu_b1=mean(b1), sigma_b1=sd(b1),
        diff_b2=wasserstein(b2), mu_b2=mean(b2), sigma_b2=sd(b2),
        diff_b3=wasserstein(b3), mu_b3=mean(b3), sigma_b3=sd(b3),
        diff_b4=wasserstein(b4), mu_b4=mean(b4), sigma_b4=sd(b4),
        diff_chp_sigma=wasserstein(chp_sigma), mu_chp_sigma=mean(chp_sigma), sigma_chp_sigma=sd(chp_sigma),
      )
    
    df %>% 
      mutate(variable=varbl)
  }
  
  a <- extract_samp_stmts(fit)
  multiplier <- 5
  a_normal_only <- a[str_detect(a, "normal")]
  
  for(i in seq_along(a_normal_only)) {
    print(i)
    tmp <- sensitivity_analysis(a_normal_only[[i]], fit, multiplier)
    if(i == 1)
      big_df <- tmp
    else
      big_df <- big_df %>% bind_rows(tmp)
  }
  
  # get mean parameter estimates
  big_df1 <- big_df %>% 
    mutate(multiplier=replace_na(multiplier, 1)) %>% 
    dplyr::select(variable, multiplier, .pareto_k, contains("mu")) %>% 
    pivot_longer(-c(variable, multiplier, .pareto_k)) %>% 
    mutate(name=substr(name, 4, nchar(name)))
  
  original_estimates <- big_df1 %>% 
    filter(multiplier==1) %>% 
    select(-multiplier) %>% 
    pivot_wider(id_cols = c(variable, .pareto_k)) %>% 
    slice(1:1) %>% 
    mutate(variable="original")
  
  new_estimates <- big_df1 %>% 
    filter(multiplier==5) %>% 
    select(-multiplier) %>% 
    pivot_wider(id_cols = c(variable, .pareto_k)) %>% 
    arrange(desc(.pareto_k))
  
  complete_estimates <- original_estimates %>% 
    bind_rows(new_estimates)
  mean_estimates <- complete_estimates
  
  # get wassterstein distances
  big_df1 <- big_df %>% 
    mutate(multiplier=replace_na(multiplier, 1)) %>% 
    dplyr::select(variable, multiplier, .pareto_k, contains("diff")) %>% 
    pivot_longer(-c(variable, multiplier, .pareto_k)) %>% 
    mutate(name=substr(name, 6, nchar(name)))
  
  new_estimates <- big_df1 %>% 
    filter(multiplier==5) %>% 
    select(-multiplier) %>% 
    pivot_wider(id_cols = c(variable, .pareto_k)) %>% 
    arrange(desc(.pareto_k))
  
  wasserstein_estimates <- new_estimates
  
  list(mean=mean_estimates, wasserstein=wasserstein_estimates)
}

wasserstein_heatmap <- function(wasserstein_estimates) {
  
  var_order <- rev(wasserstein_estimates$variable)
  greek_labels <- c(
    'l0'=expression(l[0]),
    'gamma'=expression(gamma),
    'k_h'=expression(kappa[h]),
    'k_lm'=expression(k[lm]),
    'k_mh'=expression(k[mh]),
    'k_m'=expression(kappa[m]),
    'x_star'='x*',
    'eta'=expression(eta),
    'zeta'=expression(zeta),
    'alpha_h'=expression(alpha[h]),
    'alpha_m'=expression(alpha[m]),
    'b1'='c',
    'b2'='k',
    'b3'='b',
    'b4'='q',
    'chp_sigma'=expression(sigma[c])
  )
  df <- wasserstein_estimates %>% 
    dplyr::select(-.pareto_k) %>% 
    pivot_longer(-variable)
  df %>% 
    mutate(
      name=as.factor(name),
      variable=as.factor(variable)
      ) %>% 
    mutate(
      name=fct_rev(fct_relevel(name, var_order)),
      variable=fct_relevel(variable, var_order)
      ) %>% 
    ggplot(aes(x=name, y=variable)) +
    geom_tile(aes(fill=value)) +
    scale_x_discrete(
      labels=greek_labels
      ) +
    scale_y_discrete(
      labels=greek_labels
    ) +
    xlab("") +
    ylab("") +
    scale_fill_viridis_c("Wasserstein\ndist.")
}


