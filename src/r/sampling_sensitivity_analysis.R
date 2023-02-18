library(adjustr)

sensitivity_analysis <- function(old_prior_formula, fit, multiplier) {
  forms <- Reduce(paste, deparse(old_prior_formula))
  varbl_end <- unname(str_locate(forms, "~")[1, 1]) - 2
  varbl <<- substr(forms, 1, varbl_end) # for some reason varbl is only found it its a global
  new_form <- paste0(substr(forms, 1, nchar(forms) - 1), "* multiplier)")
  df <- make_spec(as.formula(new_form), multiplier=multiplier) %>% 
    adjust_weights(fit, keep_bad = TRUE) %>% # keep bad since otherwise summary functions throw error if pareto k > 0.7
    summarize(diff=wasserstein(eval(as.symbol(varbl))), mu=mean(eval(as.symbol(varbl))), sigma=sd(eval(as.symbol(varbl))))
  
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

big_df1 <- big_df %>% 
  mutate(multiplier=replace_na(multiplier, 1)) %>% 
  dplyr::select(variable, multiplier, .pareto_k, diff, mu, sigma)


make_spec(eta ~ normal(1, 2)) %>%
  adjust_weights(fit) %>% 
  summarize(wasserstein(eta), mean(eta), sd(eta))
