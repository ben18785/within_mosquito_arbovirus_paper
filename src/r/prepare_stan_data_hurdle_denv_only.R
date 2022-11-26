
prepare_stan_data_hurdle_denv_only <- function(
    df_midgut_legs,
    df_denv_dilutions_infected,
    df_disseminated_infection_time_course) {
  
  partial_refeed_amount <- 0.56
  both <- df_midgut_legs %>%
    mutate(t_ind = as.numeric(as.factor(day)),
           difeq_ind = ifelse(tissue == "midgut", 2, 3),
           dilution_ind = as.numeric(as.factor(dilutions)))
  minimal_denv_detectable <- 35
  
  # set observations (only one) < 35 to 35.001
  both <- both %>% 
    mutate(denv_titer=if_else(denv_titer < minimal_denv_detectable, minimal_denv_detectable + 0.001, denv_titer))
  
  # look for individuals with -ve midgut but positive legs
  mismatch <- both %>% 
    group_by(day, sample, individual) %>%
    arrange(day, sample, individual) %>% 
    summarise(midgut_na=is.na(denv_titer[1]),
              legs_na=is.na(denv_titer[2])) %>% 
    filter(midgut_na) %>% 
    filter(!legs_na)
  stopifnot(nrow(mismatch)==0)
  # both NAs
  indivs_nas <- both %>% 
    group_by(day, sample, individual) %>% 
    summarise(n_na=sum(is.na(denv_titer))) %>% 
    filter(n_na==2)
  
  # since those specimens with midgut == 0 also have zero disseminated infection, include these
  a <- both %>%
    select(individual, day, dilutions, tissue, denv_titer) %>%
    pivot_wider(id_cols = c(individual, dilutions, day),
                names_from = tissue,
                values_from=denv_titer) %>% 
    pivot_longer(c(midgut, legs)) %>% 
    rename(tissue=name,
           denv_titer=value) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs"))
  
  # transform denv titer data to be nearer unit scale
  sum_df <- a %>% 
    group_by(tissue) %>% 
    summarise(overall_denv_titer=mean(denv_titer, na.rm=T))
  titer_lower_bound <- minimal_denv_detectable / sum_df$overall_denv_titer
  threshold <- -0.1
  a <- a %>% 
    left_join(sum_df) %>% 
    mutate(denv_titer=denv_titer / overall_denv_titer) %>% 
    mutate(denv_titer=if_else(is.na(denv_titer), threshold, denv_titer)) %>% 
    mutate(t_ind = as.numeric(as.factor(day)),
           difeq_ind = ifelse(tissue == "midgut", 2, 3),
           dilution_ind = as.numeric(as.factor(dilutions)))
  
  # process binary titer data: positive or negative
  df_body <- df_denv_dilutions_infected %>% 
    mutate(tissue="midgut") # assume that if body infected so is midgut
  df_multiple <- df_disseminated_infection_time_course %>%
    bind_rows(df_body) %>% 
    mutate(refeed_time=if_else(is.na(refeed_time), 100, refeed_time)) %>% 
    mutate(index=seq_along(refeed_time))
  
  n_unq_t_binary = n_distinct(df_multiple$day)
  ts_binary <- sort(unique(df_multiple$day))
  unique_experiments <- df_multiple %>% 
    select(day, type, refeed_time, dilutions) %>% 
    unique() %>% 
    mutate(refeed_amount=if_else(type=="partial", partial_refeed_amount, 1.0)) %>%
    mutate(index_new=seq_along(refeed_amount))
  n_experiment_types_binary <- nrow(unique_experiments)
  dilutions_binary <- unique_experiments$dilutions
  t_refeed_binary <- unique_experiments$refeed_time
  refeed_amount_binary <- unique_experiments$refeed_amount
  n_obs_binary <- nrow(df_multiple)
  df_multiple <- df_multiple %>% 
    left_join(unique_experiments)
  t_ind_binary <- match(df_multiple$day, ts_binary)
  n_infected_binary <- df_multiple$n_infected
  n_dissected_binary <- df_multiple$n_dissected
  ind_binary <- df_multiple$index_new
  difeq_ind_binary <- if_else(df_multiple$tissue=="midgut", 2, 3)
  data_in <- list(n_unq_t = as.integer(length(unique(a$day))),
                  n_theta = as.integer(12),
                  n_difeq = as.integer(3),
                  n_dilutions = as.integer(length(unique(a$dilutions))),
                  t0 = as.double(0),
                  ts = unique(a$day),
                  d_r_in = as.vector(max(unique(a$day))+10), # simulating without the refeeding
                  dilutions = as.double(unique(a$dilutions)),
                  n_obs = as.integer(nrow(a)),
                  y = as.double(a$denv_titer),
                  t_ind = as.integer(a$t_ind),
                  difeq_ind = as.integer(a$difeq_ind),
                  dilution_ind = as.integer(a$dilution_ind),
                  x_star = as.double(0.1),
                  eta = as.double(2),
                  g_t = seq(0.5, max(unique(a$day)) + 10, 0.5),
                  n_g_t = length(seq(0.5, max(unique(a$day))+10, 0.5)),
                  x_0=1,
                  x_r=1,
                  titer_lower_bound=titer_lower_bound,
                  n_unq_t_binary=n_unq_t_binary,
                  ts_binary=ts_binary,
                  n_experiment_types_binary=n_experiment_types_binary,
                  dilutions_binary=dilutions_binary,
                  t_refeed_binary=t_refeed_binary,
                  refeed_amount_binary=refeed_amount_binary,
                  n_obs_binary=n_obs_binary,
                  t_ind_binary=t_ind_binary,
                  n_infected_binary=n_infected_binary,
                  n_dissected_binary=n_dissected_binary,
                  ind_binary=ind_binary,
                  difeq_ind_binary=difeq_ind_binary,
                  include_continuous=1,
                  include_continuous_censored=1,
                  include_binary=1,
                  include_hurdle=1,
                  n_refeed_types=2,
                  t_refeeds_unique=c(3, 100),
                  refeed_amounts_unique=c(1.0, 1.0),
                  n_dilutions_sim=5,
                  dilutions_sim=c(1, 5, 12, 20, 25))
  
  list(dataset_binary_denv=df_multiple,
       dataset_denv=both,
       dataset_denv_expanded=a,
       dataset_denv_sum=sum_df,
       stan_data=data_in)
}