pars_to_keep <- c(
  "l0", "zeta", "gamma", "a", "k_lm",
  "alpha_m", "k_m", "k_mh",
  "alpha_h", "k_h", "eta",
  "x_star", "chp_vals",
  "b3", "b4",
  "phi_d", "chp_sigma",
  "sigma")

f_create_kde <- function(x, method="bkde", bandwidth=NULL, ...){
  bandwidth1 <- bandwidth
  if(method=="bkde") {
    if(is.null(bandwidth))
      a <- bkde(x, ...)
    else
      a <- bkde(x, bandwidth = bandwidth1, ...)
  }else if(method=="ash1") {
    if(is.null(bandwidth))
      a <- ash2(bin1(x), ...)
    else
      a <- ash2(bin1(x, nbin = round(bandwidth1 * 50)), ...)
  }else{
    stop("Method must either be bkde or ash1.")
  }
  f <- approxfun(a$x, a$y, method="constant")
  kde <- function(x) {
    f_val <- f(x)
    return(if_else(is.na(f_val), 0, f_val))
  }
  kde
}


# to handle nas
f_compare_estimates <- function(prob1, prob2) {
  k <- 0
  n <- length(prob1)
  for(i in seq_along(prob1))
    if((!is.na(prob1[i])) && (!is.na(prob2[i])))
      k <- k + if_else(prob1[i] > prob2[i], 1, 0)
  else if((!is.na(prob1[i])))
    k <- k + 1
  else if((is.na(prob1[i])) && (is.na(prob2[i])))
    n <- n - 1
  k / n
}

# calculate accuracy across training and testing sets
dstar <- function(x1, x2, f1=NULL, f2=NULL,
                  method="bkde",
                  training_percent=0.7,
                  bandwidth_1=NULL,
                  bandwidth_2=NULL,
                  ...) {
  # check if f1 or f2 are defined (i.e. a probability density)
  # if so check they must be functions
  if(!is.null(f1))
    isfunc1 <- is.function(f1)
  else
    isfunc1 <- NULL
  if(!is.null(f2))
    isfunc2 <- is.function(f2)
  else
    isfunc2 <- NULL
  
  if(!is.null(isfunc1))
    if(!isfunc1)
      stop("f1 must be a function or NULL type.")
  if(!is.null(isfunc2))
    if(!isfunc2)
      stop("f2 must be a function or NULL type.")
  if(is.null(isfunc1))
    isfunc1 <- FALSE
  if(is.null(isfunc2))
    isfunc2 <- FALSE
  
  if(isfunc1 && isfunc2) {
    # use whole sets given
    prob11 <- map_dbl(x1, f1)
    prob12 <- map_dbl(x1, f2)
    prob21 <- map_dbl(x2, f1)
    prob22 <- map_dbl(x2, f2)
    
  } else {
    # training and testing sets
    idxs <- seq_along(x1)
    training_ids <- sample(idxs,
                           size = round(training_percent * length(x1)))
    testing_ids <- setdiff(idxs, training_ids)
    x1_train <- x1[training_ids]
    x2_train <- x2[training_ids]
    x1_test <- x1[testing_ids]
    x2_test <- x2[testing_ids]
    if(isfunc1) {
      f2 <- f_create_kde(x2_train, 
                         method,
                         bandwidth=bandwidth_2,
                         ...)
    } else if(isfunc2) {
      f1 <- f_create_kde(x1_train,
                         method,
                         bandwidth=bandwidth_1,
                         ...)
    } else {
      f1 <- f_create_kde(x1_train,
                         method,
                         bandwidth=bandwidth_1,
                         ...)
      f2 <- f_create_kde(x2_train,
                         method,
                         bandwidth=bandwidth_2,
                         ...)
    }
    prob11 <- map_dbl(x1_test, f1)
    prob12 <- map_dbl(x1_test, f2)
    prob21 <- map_dbl(x2_test, f1)
    prob22 <- map_dbl(x2_test, f2)
  }
  # calculate accuracy
  acc1 <- f_compare_estimates(prob11, prob12)
  acc2 <- f_compare_estimates(prob22, prob21)
  acc_overall <- 0.5 * (acc1 + acc2)
  
  # calculate score
  score <- 2 * (acc_overall - 0.5)
  score <- if_else(score < 0, 0, score)
  score
}

sample_prior <- function(n, prior_formula) {
  str_prior <- Reduce(paste, deparse(prior_formula))
  loc_normal <- str_locate(str_prior, "normal")
  str_normal <- substr(str_prior, loc_normal[1], nchar(str_prior))
  str_sample <- str_replace(str_normal, "normal\\(", "rnorm\\(n, ")
  eval(parse(text=str_sample))
}

find_variable_name <- function(prior_formula) {
  str_prior <- Reduce(paste, deparse(prior_formula))
  varbl_end <- unname(str_locate(str_prior, "~")[1, 1]) - 2
  print(varbl_end)
  substr(str_prior, 1, varbl_end)
}

prior_to_posterior_summary <- function(fit) {
  
  library(KernSmooth)
  
  # find priors
  a <- extract_samp_stmts(fit)
  a_normal_only <- a[str_detect(a, "normal")]
  
  # for each parameter calculate d*
  dstars <- vector(length = length(a_normal_only))
  param_names <- vector(length = length(a_normal_only))
  n <- 1000
  for(i in seq_along(dstars)) {
    print(a_normal_only[[i]])
    a_var <- find_variable_name(a_normal_only[[i]])
    param_names[i] <- a_var
    x_prior <- sample_prior(n, a_normal_only[[i]])
    x_post <- rstan::extract(fit, a_var)[[1]]
    dstars[i] <- dstar(x_prior, x_post)
  }
  df_normal <- tibble(variable=param_names, dstar=dstars) %>% 
    mutate(variable=case_when(
      variable=="b3"~"b",
      variable=="b4"~"q",
      TRUE~variable
    ))
  print("hiya")
  # bespoke ones where prior extraction doesn't work
  ## phi_d
  x_prior <- runif(n)
  x_posterior <- rstan::extract(fit, "phi_d")[[1]]
  d_val <- dstar(x_prior, x_posterior)
  df_phi <- tibble(variable="phi_d", dstar=d_val)
  
  ## sigmas
  x_prior <- rcauchy(2 * n, 0, 5)
  x_prior <- x_prior[x_prior > 0]
  x_posterior <- rstan::extract(fit, "sigma[1]")[[1]]
  d_val1 <- dstar(x_prior, x_posterior)
  x_posterior <- rstan::extract(fit, "sigma[2]")[[1]]
  d_val2 <- dstar(x_prior, x_posterior)
  df_sigmas <- tibble(variable=c("sigma.1", "sigma.2"),
                      dstar=c(d_val1, d_val2))
  
  ## CHP vals
  x_prior <- rnorm(n, 300, 100)
  x_posterior <- rstan::extract(fit, "chp_vals[1]")[[1]]
  chp1 <- dstar(x_prior, x_posterior)
  x_prior <- rnorm(n, 3000, 500)
  x_posterior <- rstan::extract(fit, "chp_vals[2]")[[1]]
  chp2 <- dstar(x_prior, x_posterior)
  df_chps <- tibble(variable=c("chp_vals.1", "chp_vals.2"),
                    dstar=c(chp1, chp2))
  
  df_dstar <- df_normal %>% 
    bind_rows(
      df_phi,
      df_sigmas,
      df_chps
    )
  
  df_dstar
}

posterior_summary <- function(fit) {
  
  df_pars <- rstan::extract(fit, pars_to_keep) %>% 
    as.data.frame()
  
  df_sum <- summarize_draws(df_pars) %>% 
    mutate(variable=case_when(
      variable=="b3"~"b",
      variable=="b4"~"q",
      TRUE~variable
    ))
  
  df_dstar <- prior_to_posterior_summary(fit)
  
  df_sum <- df_sum %>% 
    left_join(df_dstar) %>% 
    mutate(variable=case_when(
      variable=="chp_vals.1"~"c*",
      variable=="chp_vals.2"~"c_0",
      variable=="sigma.1"~"sigma_m",
      variable=="sigma.2"~"sigma_h",
      TRUE~variable
    )) %>% 
    mutate(across(!variable,
                  function(x) format(round(x, 3), nsmall = 3)))
  
  df_sum
}

posteriors_correlation <- function(fit) {
  
  df_pars <- rstan::extract(fit, pars_to_keep) %>% 
    as.data.frame()
  cc <- cor(df_pars) %>% 
    as.matrix()
  rowms <- rowMeans(cc)
  vars <- tibble(
    variable=row.names(cc),
    value=rowms
  ) %>% 
    arrange(desc(value))
  
  # # find high correlations
  # w <- which(abs(cc)>0.57 & row(cc)<col(cc), arr.ind=TRUE)
  # high_cor <- matrix(colnames(cc)[w], ncol=2) %>% 
  #   as.data.frame() %>% 
  #   mutate(name_long=paste0(V1, "_", V2))
  # all_vars <- unique(c(high_cor$V1, high_cor$V2))
  # all_vars[all_vars=="chp_vals.2"] <- "chp_vals[2]"
  
  
  df <- df_pars %>% 
    dplyr::select(all_of(vars$variable[1:7])) %>% 
    rename(
      sigma_h=sigma.2
      )
  
  GGally::ggpairs(df,
                  lower = list(continuous = GGally::wrap("points", alpha = 0.3, size=0.5)),
                  upper = list(continuous = GGally::wrap("cor", size = 4))
  ) + theme(strip.text = element_text(size = 7))
}
