pars_to_keep <- c(
  "l0", "zeta", "gamma", "a",
  "alpha_m", "k_m", "k_mh",
  "alpha_h", "k_h", "eta",
  "x_star", "chp_vals",
  "b1", "b2", "b3", "b4",
  "phi_d", "chp_sigma",
  "sigma")

posterior_summary <- function(fit) {
  df_pars <- rstan::extract(fit, pars_to_keep) %>% 
    as.data.frame()
  
  summarize_draws(df_pars)
}

posteriors_correlation <- function(fit) {
  
  df_pars <- rstan::extract(fit, pars_to_keep) %>% 
    as.data.frame()
  cc <- cor(df_pars) %>% 
    as.matrix()
  
  # find high correlations
  w <- which(abs(cc)>0.57 & row(cc)<col(cc), arr.ind=TRUE)
  high_cor <- matrix(colnames(cc)[w], ncol=2) %>% 
    as.data.frame() %>% 
    mutate(name_long=paste0(V1, "_", V2))
  all_vars <- unique(c(high_cor$V1, high_cor$V2))
  all_vars[all_vars=="chp_vals.2"] <- "chp_vals[2]"
  
  df <- df_pars %>% 
    dplyr::select(all_of(all_vars)) %>% 
    rename(
      c=b1,
      b=b3,
      q=b4,
      "x*"=x_star
      )
  
  GGally::ggpairs(df,
                  lower = list(continuous = GGally::wrap("points", alpha = 0.3, size=0.5))
  )
}
