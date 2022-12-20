
logistic_regression_comparison <- function(list_stan_datasets) {
  
  df_continuous <- list_stan_datasets$dataset_denv %>% 
    mutate(denv_above_zero=if_else(denv_titer > 0,
                                   "non-censored", "below detection\nthreshold")) %>% 
    group_by(tissue, dilutions, day) %>% 
    summarise(n_positive=sum(denv_above_zero=="non-censored", na.rm = TRUE),
              n_dissected=n()) %>% 
    mutate(type="single") %>% 
    rename(dilution=dilutions) %>% 
    mutate(data_type="DENV") %>% 
    mutate(category="continuous")
  
  df_dichtonomous <- list_stan_datasets$dataset_binary_denv %>% 
    rename(n_positive=n_infected) %>% 
    mutate(dilution=dilutions) %>% 
    select(type, tissue, dilution, day, n_positive, n_dissected, virus) %>% 
    rename(data_type=virus) %>% 
    filter(type=="single") %>% 
    mutate(category="dichtonomous")
  
  df_both <- df_dichtonomous %>% 
    bind_rows(df_continuous) %>% 
    filter(tissue=="legs")
  
  model_0 <- glm(
    cbind(n_positive, n_dissected-n_positive)~day + category + dilution,
    data=df_both, family = binomial)
  model_1 <- glm(
    cbind(n_positive, n_dissected-n_positive)~day + category,
    data=df_both, family = binomial)
  
  tmp <- lrtest(model_0, model_1)
  list(model_0=model_0,
       model_1=model_1,
       lr_test=tmp)
}

