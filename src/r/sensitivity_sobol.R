
sensitivity_sobol <- function(fit) {
  
  N <- 2 ^ 8 # number of points at which ODE is solved
  params <- c(
    "gamma",
    "k_lm",
    "a",
    "eta",
    "x_star",
    "alpha_m",
    "k_m",
    "k_mh",
    "alpha_h",
    "k_h")
  order <- "first"
  R <- 10 ^ 3
  type <- "norm"
  conf <- 0.95
  times <- seq(0, 15, 1)
  timeOutput <- seq(1, 15, 1)
  
  rhs_complex <- function(t, state, parameters) {
    
    x_0 <- 1
    t_refeed <- 100
    
    with(as.list(c(state, parameters)), {
      x_t = x_star + exp(-eta * t) * (x_0 - x_star)
      if(t < t_refeed)
        x_t <- x_star + exp(-eta*t) * (x_0 - x_star)
      else
        x_t <- x_star + exp(-eta*(t - t_refeed)) * (x_0 - x_star)
      dl_dt = -gamma * l - x_t * (k_lm * l^2) / (a^2 + l^2)
      dm_dt = x_t * (k_lm * l^2) / (a^2 + l^2) + alpha_m * m * (1 - m / k_m) - k_mh * (x_t - x_star) * m
      dh_dt = k_mh * (x_t - x_star) * m + alpha_h * h * (1 - h/k_h)
      list(c(dl_dt, dm_dt, dh_dt))
    })
    
  }
  
  summary_func <- median
  gamma <- summary_func(rstan::extract(fit, "gamma")[[1]])
  k_lm <- summary_func(rstan::extract(fit, "k_lm")[[1]])
  a <- summary_func(rstan::extract(fit, "a")[[1]])
  eta <- summary_func(rstan::extract(fit, "eta")[[1]])
  x_star <- summary_func(rstan::extract(fit, "x_star")[[1]])
  alpha_m <- summary_func(rstan::extract(fit, "alpha_m")[[1]])
  k_m <- summary_func(rstan::extract(fit, "k_m")[[1]])
  k_mh <- summary_func(rstan::extract(fit, "k_mh")[[1]])
  alpha_h <- summary_func(rstan::extract(fit, "alpha_h")[[1]])
  k_h <- summary_func(rstan::extract(fit, "k_h")[[1]])
  l0 <- summary_func(rstan::extract(fit, "l0")[[1]])
  
  mat <- sobol_matrices(N = N, params = params, order = order)
  lower <- 2 / 3
  upper <- 3 / 2
  mat[, "gamma"] <- qunif(mat[, "gamma"], gamma * lower, gamma * upper)
  mat[, "k_lm"] <- qunif(mat[, "k_lm"], k_lm * lower, k_lm * upper)
  mat[, "a"] <- qunif(mat[, "a"], a * lower, a * upper)
  mat[, "eta"] <- qunif(mat[, "eta"], eta * lower, eta * upper)
  mat[, "alpha_m"] <- qunif(mat[, "alpha_m"], alpha_m * lower, alpha_m * upper)
  mat[, "k_m"] <- qunif(mat[, "k_m"], k_m * lower, k_m * upper)
  mat[, "k_mh"] <- qunif(mat[, "k_mh"], k_mh * lower, k_mh * upper)
  mat[, "alpha_h"] <- qunif(mat[, "alpha_h"], alpha_h * lower, alpha_h * upper)
  mat[, "x_star"] <- qunif(mat[, "x_star"], x_star * lower, x_star * upper)
  
  
  for(i in 1:nrow(mat)) {
    print(i)
      temp_df <- sobol_ode(d = mat[i, ], times = times, timeOutput = timeOutput,
                state = c(l = l0, m=0, h=0), func = rhs_complex) %>% 
        as.data.frame()
      if(i == 1)
        big_df <- temp_df
      else
        big_df <- big_df %>% bind_rows(temp_df)
  }
  
  long_df <- big_df %>% 
    pivot_longer(c(l, m, h))
  
  indices <- long_df %>% 
    group_by(name, time) %>% 
    summarise(
      sobol_indices(Y = value, N = N, params = params,
                          order = order, boot = TRUE, first = "jansen", R = R,
                          parallel = "multicore", ncpus = 4)$results
    )
  
  indices
}

plot_sensitivity_sobol <- function(indices) {
  
  indices1 <- indices %>% 
    ungroup() %>% 
    mutate(
      name=as.factor(name),
      name=fct_relevel(name, "l", "m", "h")
    ) %>% 
    mutate(parameters=case_when(
      parameters == "k_lm"~"k[lm]",
      parameters == "alpha_m"~"alpha[m]",
      parameters == "k_m"~"kappa[m]",
      parameters == "k_mh"~"k[mh]",
      parameters == "alpha_h"~"alpha[h]",
      parameters == "k_h"~"kappa[h]",
      parameters == "x_star"~"xstar",
      parameters == "eta"~"eta",
      TRUE~parameters
    )) %>% 
    mutate(
      parameters=as.factor(parameters),
      parameters=fct_relevel(parameters,
                             "gamma", "a", "k[lm]",
                             "alpha[m]", "kappa[m]", "k[mh]",
                             "alpha[h]", "kappa[h]",
                             "eta", "xstar")
    )
  
  
  
  indices_filt <- indices1 %>% 
    filter(!(name == "l" & time > 3))
  
  
  ggplot(indices_filt,
         aes(x=time, y=original,
         fill=sensitivity,
         group=sensitivity)) +
    geom_ribbon(aes(ymin=low.ci, ymax=high.ci),
                alpha=0.3) +
    geom_line(aes(colour=sensitivity)) +
    facet_grid(parameters ~ name, labeller = label_parsed) +
    scale_color_brewer("Sensitivity", palette = "Dark2") +
    scale_fill_brewer("Sensitivity", palette = "Dark2") +
    theme_bw() +
    xlab("DPI") +
    ylab("Sobol' indices")
}
