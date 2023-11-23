# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan",
               "cowplot", "posterior",
               "loo", "deSolve", "lmtest",
               "stargazer", "latex2exp",
               "adjustr", "bayesplot",
               "sensobol"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
source("src/r/process_experimental_data.R")
source("src/r/plot_experimental_data.R")
source("src/r/process_denv_dilutions_infected.R")
source("src/r/process_disseminated_infection_time_course.R")
source("src/r/prepare_stan_data_hurdle_denv_only.R")
source("src/r/fit_optimise.R")
source("src/r/plot_fit_prevalence.R")
source("src/r/process_damage_data.R")
source("src/r/plot_chp_damage_fit.R")
source("src/r/plot_single_double_feed.R")
source("src/r/plot_bl_permeability.R")
source("src/r/helper.R")
source("src/r/fit_sampling.R")
source("src/r/sampling_diagnostics.R")
source("src/r/model_comparison.R")
source("src/r/plot_eip.R")
source("src/r/logistic_regression.R")
source("src/r/plot_continuous_data.R")
source("src/r/plot_sensitivity.R")
source("src/r/prior_sensitivity_analysis.R")
source("src/r/posterior_summary_statistics.R")
source("src/r/prior_predictive.R")
source("src/r/sensitivity_sobol.R")

list(
  
  # overall list of figures and outputs that go into the paper
  tar_target(overall_list, {
    list(file_graph_midgut_dose_response_combined_mcmc,
         file_graph_experimental_data_midgut,
         file_graph_continuous,
         file_graph_fit_prevalence_legs_eip,
         file_graph_single_double_chp_mcmc,
         file_plot_sd_sensitivity_2d_alpha_kappa,
         file_plot_sd_sensitivity_2d_alpha_l0,
         file_graph_prior_predictives,
         file_graph_posterior_correlations,
         file_graph_fit_prevalence_midgut_only_mcmc_both_types,
         file_graph_midgut_invasion_sensitivities,
         file_graph_noninfectious_then_infectious_double,
         file_graph_sensitivities_single_double_feed,
         file_graph_sobol_indices,
         file_graph_wasserstein_estimate,
         sampling_fit_diagnostics,
         file_posteriors_summary
    )
  }),
  
  # raw experimental data processing
  tar_target(filename_midgut, "data/raw/Compiled midgut data.xlsx",
             format = "file"),
  tar_target(filename_legs, "data/raw/Compiled leg data.xlsx",
             format = "file"),
  tar_target(df_midgut_legs,
             process_experimental_data(
               filename_midgut, filename_legs)),
  tar_target(filename_denv_dilutions, "data/raw/DENV test of dilutions 8-26-22.xlsx",
             format="file"),
  tar_target(df_denv_dilutions_infected,
             process_denv_dilutions_infected(filename_denv_dilutions)),
  tar_target(filename_disseminated_time_course, "data/raw/DENV SF DF dissemination time course 8-26-22.xlsx",
             format="file"),
  tar_target(df_disseminated_infection_time_course,
             process_disseminated_infection_time_course(
               filename_disseminated_time_course)),
  
  # Plot raw experimental data
  tar_target(graph_experimental_data, plot_experimental_data(df_midgut_legs)),
  tar_target(file_graph_experimental_data, {
    ggsave("figures/experimental_titers.pdf", graph_experimental_data, width=10, height=4);
    "figures/experimental_titers.pdf"}, format="file"),
  tar_target(graph_experimental_data_midgut,
             plot_experimental_data_midgut(df_midgut_legs)),
  tar_target(file_graph_experimental_data_midgut, {
    ggsave("figures/experimental_titers_midgut.pdf", graph_experimental_data_midgut, width=6, height=4);
    "figures/experimental_titers_midgut.pdf"}, format="file"),
  tar_target(anova_tests_midgut, { # not appropriate since data aren't normal
    df <- df_midgut_legs %>% 
      filter(tissue=="midgut",
             day==3) %>% # represents negative cases
      mutate(denv_titer=if_else(is.na(denv_titer), 35, denv_titer))
    fit <- aov(denv_titer~sample, data=df)
    summary(fit)
  }),
  tar_target(t_tests_midgut, { # not appropriate since data aren't normal
    df <- df_midgut_legs %>% 
      filter(tissue=="midgut",
             day==3) %>% 
      mutate(denv_titer=if_else(is.na(denv_titer), 35, denv_titer)) # represents negative cases
    df_1 <- df %>% 
      filter(sample %in% c("1:1", "1:5"))
    fit_1 <- t.test(denv_titer~sample, data=df_1)
    df_2 <- df %>% 
      filter(sample %in% c("1:1", "1:12"))
    fit_2 <- t.test(denv_titer~sample, data=df_2)
    df_3 <- df %>% 
      filter(sample %in% c("1:5", "1:12"))
    fit_3 <- t.test(denv_titer~sample, data=df_3)
    list(one_five=fit_1, one_12=fit_2, five_12=fit_3)
  }),
  tar_target(kruskal_tests_midgut, { # most appropriate since the observations are not normally distributed
    df <- df_midgut_legs %>% 
      filter(tissue=="midgut",
             day==3) %>% 
      mutate(denv_titer=if_else(is.na(denv_titer), 35, denv_titer)) # represents negative cases
    df_1 <- df %>% 
      filter(sample %in% c("1:1", "1:5"))
    fit_1 <- kruskal.test(denv_titer~sample, data=df_1)
    df_2 <- df %>% 
      filter(sample %in% c("1:1", "1:12"))
    fit_2 <- kruskal.test(denv_titer~sample, data=df_2)
    df_3 <- df %>% 
      filter(sample %in% c("1:5", "1:12"))
    fit_3 <- kruskal.test(denv_titer~sample, data=df_3)
    list(one_five=fit_1, one_12=fit_2, five_12=fit_3)
  }),
  
  # process CHP damage data
  tar_target(filename_chp_damage, "data/raw/midgut damage over time.xlsx",
             format = "file"),
  tar_target(df_chp_damage,
             process_damage_data(filename_chp_damage)),
  
  # stan data inputs processing
  tar_target(list_stan_datasets,
             prepare_stan_data_hurdle_denv_only(
               df_midgut_legs,
               df_denv_dilutions_infected,
               df_disseminated_infection_time_course,
               df_chp_damage)),
  
  # fit stan model via optimisation
  tar_target(stan_model, "src/stan/model_full.stan",
             format="file"),
  tar_target(opt_fit,
             fit_optimise(
               list_stan_datasets$stan_data,
               stan_model, 5)),
  
  # plots based on optimisation
  tar_target(graph_fit_prevalence_midgut_legs,
             plot_fit_prevalence_midgut_legs(
               opt_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_fit_prevalence_midgut_legs, {
    ggsave("figures/prevalence_midgut_legs.pdf", graph_fit_prevalence_midgut_legs, width=10, height=6);
    "figures/prevalence_midgut_legs.pdf"}, format="file"),
  tar_target(graph_fit_prevalence_midgut_only,
             plot_fit_prevalence_midgut_only(
               opt_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_fit_prevalence_midgut_only, {
    ggsave("figures/prevalence_midgut_only.pdf", graph_fit_prevalence_midgut_only, width=6, height=3);
    "figures/prevalence_midgut_only.pdf"}, format="file"),
  tar_target(file_graph_fit_prevalence_midgut_only_png, {
    ggsave("figures/prevalence_midgut_only.png", graph_fit_prevalence_midgut_only, width=4, height=3);
    "figures/prevalence_midgut_only.png"}, format="file"),
  tar_target(graph_fit_prevalence_dose_response,
             plot_fit_prevalence_dose_response(
               opt_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_fit_prevalence_dose_response, {
    ggsave("figures/prevalence_dose_response.pdf", graph_fit_prevalence_dose_response, width=10, height=6);
    "figures/prevalence_dose_response.pdf"}, format="file"),
  tar_target(graph_midgut_dose_response_combined,
             plot_midgut_dose_response_combined(
               opt_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_midgut_dose_response_combined, {
    ggsave("figures/prevalence_midgut_dose_response.pdf", graph_midgut_dose_response_combined, width=10, height=4);
    "figures/prevalence_midgut_dose_response.pdf"}, format="file"),
  tar_target(graph_chp_damage, plot_chp_damage_fit(opt_fit, df_chp_damage)),
  tar_target(file_graph_chp_damage, {
    ggsave("figures/chp_damage.pdf", graph_chp_damage, width=10, height=6);
    "figures/chp_damage.pdf"}, format="file"),
  tar_target(graph_single_double, plot_single_double_feed(opt_fit, list_stan_datasets)),
  tar_target(file_graph_single_double, {
    ggsave("figures/single_double.png", graph_single_double, width=4, height=3);
    "figures/single_double.png"}, format="file"),
  tar_target(graph_bl_permeability,
             plot_bl_permeability(
               opt_fit,
               list_stan_datasets)),
  tar_target(graph_single_double_chp, 
             {
               plot_grid(graph_chp_damage, graph_bl_permeability, graph_single_double, nrow = 1,
                         labels = c("A.", "B.", "C."))
             }),
  tar_target(file_graph_single_double_chp, {
    ggsave("figures/graph_single_double_chp.pdf", graph_single_double_chp, width=12, height=4);
    "figures/graph_single_double_chp.pdf"}, format="file"),
  tar_target(graph_noninfectious_then_infectious_double,
             plot_noninfectious_then_infectious_double(
               opt_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_noninfectious_then_infectious_double, {
    ggsave("figures/double_other_order.pdf", graph_noninfectious_then_infectious_double, width=8, height=4);
    "figures/double_other_order.pdf"}, format="file"),
  tar_target(graph_fit_prevalence_legs,
             plot_fit_prevalence_legs(
               opt_fit,
               list_stan_datasets)),
  tar_target(file_graph_fit_prevalence_legs, {
    ggsave("figures/prevalence_legs.pdf", graph_fit_prevalence_legs, width=6, height=4);
    "figures/prevalence_legs.pdf"}, format="file"),
  
  
  # fit stan model via sampling
  tar_target(sampling_fit,
             fit_mcmc(opt_fit,
                      stan_model,
                      list_stan_datasets$stan_data,
                      n_iterations=600,
                      n_chains=4)
             ),
  
  # MCMC diagnostics
  tar_target(sampling_fit_diagnostics,
             create_sampling_diagnostics(sampling_fit)),
  
  # Plots based on MCMC fitting
  
  # Single-double feed
  tar_target(graph_single_double_mcmc,
             plot_single_double_feed_mcmc(
               sampling_fit, list_stan_datasets)),
  tar_target(graph_chp_damage_mcmc, plot_chp_damage_fit_mcmc(sampling_fit, df_chp_damage)),
  tar_target(graph_bl_permeability_mcmc,
             plot_bl_permeability_mcmc(
               sampling_fit,
               list_stan_datasets)),
  tar_target(graph_single_double_chp_mcmc, 
             {
               plot_grid(graph_chp_damage_mcmc, graph_bl_permeability_mcmc, graph_single_double_mcmc, nrow = 1,
                         labels = c("A.", "B.", "C."))
             }),
  tar_target(file_graph_single_double_chp_mcmc, {
    ggsave("figures/graph_single_double_chp_mcmc.pdf", graph_single_double_chp_mcmc, width=12, height=4);
    "figures/graph_single_double_chp_mcmc.pdf"}, format="file"),
  
  ## Concentration plots
  tar_target(graph_fit_prevalence_midgut_only_mcmc,
             plot_fit_prevalence_midgut_only_mcmc(
               sampling_fit,
               list_stan_datasets
             )),
  tar_target(graph_fit_prevalence_dose_response_mcmc,
             plot_fit_prevalence_dose_response_mcmc(
               sampling_fit,
               list_stan_datasets
             )),
  tar_target(graph_midgut_dose_response_combined_mcmc, 
             {
               plot_grid(graph_fit_prevalence_midgut_only_mcmc,
                         graph_fit_prevalence_dose_response_mcmc,
                         nrow = 1,
                         labels = c("A.", "B."),
                         label_x = -0.01)
             }),
  tar_target(file_graph_midgut_dose_response_combined_mcmc, {
    ggsave("figures/prevalence_midgut_dose_response_mcmc.pdf",
           graph_midgut_dose_response_combined_mcmc, width=10, height=4);
    "figures/prevalence_midgut_dose_response_mcmc.pdf"}, format="file"),
  tar_target(graph_fit_prevalence_midgut_only_mcmc_both_types,
             plot_fit_prevalence_midgut_only_mcmc_both_types(
               sampling_fit,
               list_stan_datasets
             )),
  tar_target(file_graph_fit_prevalence_midgut_only_mcmc_both_types, {
    ggsave("figures/prevalence_midgut_both_types.pdf",
           graph_fit_prevalence_midgut_only_mcmc_both_types,
           width=10, height=4)
  }),
  
  # eip dose-response
  tar_target(df_eip_dose_response,
             eip_dose_response(
               sampling_fit, list_stan_datasets$stan_data,
               500)),
  tar_target(graph_eip_dose_response,
             plot_eip_dose_response(df_eip_dose_response)),
  
  # legs prevalence mcmc
  tar_target(graph_fit_prevalence_legs_mcmc,
                        plot_fit_prevalence_legs_mcmc(
                          sampling_fit,
                          list_stan_datasets)),
  
  # combine eip plot with dissemination plot
  tar_target(graph_fit_prevalence_legs_eip,
             {
               plot_grid(graph_fit_prevalence_legs_mcmc,
                         graph_eip_dose_response,
                         nrow = 1,
                         labels = c("A.", "B."),
                         label_x = -0.01)
             }),
  tar_target(file_graph_fit_prevalence_legs_eip, {
    ggsave("figures/prevalence_legs_eip_mcmc.pdf",
           graph_fit_prevalence_legs_eip, width=10, height=4);
    "figures/prevalence_legs_eip_mcmc.pdf"}, format="file"),
  
  # fitting for model comparison
  ## fit model without different kappa_m values per dose
  tar_target(stan_model_bare, "src/stan/model_log_likelihood.stan",
             format="file"),
  tar_target(sampling_fit_log_likelihood,
             fit_mcmc(opt_fit,
                      stan_model_bare,
                      list_stan_datasets$stan_data,
                      n_iterations=400,
                      n_chains=4)),
  ## fit model with different kappa_m values per dose
  tar_target(stan_model_kappa, "src/stan/model_log_likelihood_different_kappa.stan",
             format="file"),
  tar_target(opt_fit_kappa,
             fit_optimise_kappa(
               list_stan_datasets$stan_data,
               stan_model_kappa, 5)),
  tar_target(sampling_fit_log_likelihood_kappa,
             fit_mcmc_kappa(opt_fit_kappa,
                      stan_model_kappa,
                      list_stan_datasets$stan_data,
                      n_iterations=400,
                      n_chains=4)),
  tar_target(model_comparison_kappa,
             model_comparison(sampling_fit_log_likelihood,
                              sampling_fit_log_likelihood_kappa)),
  
  # fit logistic model to determine if dose effect
  tar_target(logistic_results, logistic_regression_comparison(list_stan_datasets)),
  tar_target(logistic_results_tex, {stargazer::stargazer(logistic_results$model_1, logistic_results$model_0, out="figures/logistic_regression.tex")}),
  
  # plot continuous data and fit
  tar_target(graph_continuous, plot_continuous_data(sampling_fit, list_stan_datasets)),
  tar_target(file_graph_continuous, {
    ggsave("figures/fit_vs_continuous.pdf", graph_continuous, width=8, height=5);
    "figures/fit_vs_continuous.pdf"}, format="file"),
  
  # prior predictive
  tar_target(graph_prior_phi, prior_predictive_phi(sampling_fit)),
  tar_target(graph_prior_logistic_growth_midgut,
             prior_predictive_logistic_growth(
               alpha_est=mean(rstan::extract(sampling_fit, "alpha_m")[[1]]),
               kappa_est=mean(rstan::extract(sampling_fit, "k_m")[[1]]),
               mu_alpha=3, sigma_alpha=10,
               mu_kappa=1, sigma_kappa=10
             )),
  tar_target(graph_prior_logistic_growth_legs,
             prior_predictive_logistic_growth(
               alpha_est=mean(rstan::extract(sampling_fit, "alpha_h")[[1]]),
               kappa_est=mean(rstan::extract(sampling_fit, "k_h")[[1]]),
               mu_alpha=1.5, sigma_alpha=10,
               mu_kappa=1, sigma_kappa=10
             )),
  tar_target(graph_prior_predictives, {
    plot_grid(graph_prior_phi, graph_prior_logistic_growth_midgut, graph_prior_logistic_growth_legs,
              labels=c("A.", "B.", "C."),
              nrow = 1,
              label_x = -0.01)
  }),
  tar_target(file_graph_prior_predictives, {
    ggsave("figures/prior_predictive.pdf", graph_prior_predictives, width = 10, height = 4);
    "figures/prior_predictive.pdf"
  }),
  
  # sensitivity plots
  tar_target(graph_midgut_invasion_sensitivities,
             plot_sensitivities_midgut_invasion(sampling_fit, list_stan_datasets)),
  tar_target(file_graph_midgut_invasion_sensitivities, {
    ggsave("figures/sensitivities_midgut_invasion.pdf", graph_midgut_invasion_sensitivities, width=8, height=5);
    "figures/sensitivities_midgut_invasion.pdf"}, format="file"),
  tar_target(graph_sensitivities_single_double_feed,
             plot_sensitivities_single_double_feed(sampling_fit, list_stan_datasets)),
  tar_target(file_graph_sensitivities_single_double_feed, {
    ggsave("figures/sensitivities_single_double_feed.pdf", graph_sensitivities_single_double_feed, width=8, height=5);
    "figures/sensitivities_single_double_feed.pdf"}, format="file"),
  tar_target(sd_sensitivity_2d_alpha_kappa, single_double_sensitivity_2d(c("alpha_m", "k_mh"),
                                                                         seq(0.1, 5, length.out=25),
                                                                         seq(0.1, 10, length.out=25),
                                                                         sampling_fit, list_stan_datasets)),
  tar_target(plot_sd_sensitivity_2d_alpha_kappa, plot_single_double_sensitivity_2d(sd_sensitivity_2d_alpha_kappa, c("alpha_m", "k_{mh}"))),
  tar_target(file_plot_sd_sensitivity_2d_alpha_kappa, {
    ggsave("figures/sensitivity_2d_single_double_feed.pdf", plot_sd_sensitivity_2d_alpha_kappa, width=8, height=5);
    "figures/sensitivity_2d_single_double_feed.pdf"}, format="file"),
  tar_target(sd_sensitivity_2d_alpha_l0, single_double_sensitivity_2d(c("alpha_m", "l0"),
                                                                         seq(0.1, 5, length.out=25),
                                                                         seq(0.1, 10, length.out=25),
                                                                         sampling_fit,
                                                                         list_stan_datasets,
                                                                         non_l0_2nd=FALSE)),
  tar_target(plot_sd_sensitivity_2d_alpha_l0, plot_single_double_sensitivity_2d(sd_sensitivity_2d_alpha_l0, c("alpha_m", "l_{0}"))),
  tar_target(file_plot_sd_sensitivity_2d_alpha_l0, {
    filename <- "figures/sensitivity_2d_single_double_feed_l0.pdf"
    ggsave(filename, plot_sd_sensitivity_2d_alpha_l0, width=8, height=5);
    filename}, format="file"),
  ## combine both 2d sensitivities
  tar_target(sd_sensitivity_both, {
    df_1 <- sd_sensitivity_2d_alpha_kappa
    colnames(df_1)[1:2] <- c("V1", "V2")
    df_2 <- sd_sensitivity_2d_alpha_l0
    colnames(df_2)[1:2] <- c("V1", "V2")
    df_1 <- df_1 %>% 
      mutate(type="a")
    df_2 <- df_2 %>% 
      mutate(type="b")
    df_1 %>% 
      bind_rows(df_2)
  }),
  tar_target(plot_sd_sensitivity_both, {
    
    max_val <- max(sd_sensitivity_both$value)
    sd_sensitivity_both1 <- sd_sensitivity_both %>%
      mutate(value=if_else(value < 0, 0, value))

    breaks_lower <- seq(0, ceiling(max_val) - 1, length.out=10)
    breaks_upper <- seq(1, ceiling(max_val), length.out=10)
    breaks_mid <- 0.5 * (breaks_lower + breaks_upper)

    g <- ggplot(sd_sensitivity_both1, aes(x=V1, y=V2)) +
      geom_contour_filled(aes(z=value), breaks=seq(0, 2.2, 0.2)) +
      geom_point(data=tibble(V1=1, V2=1), colour="red", size=3) +
      theme_bw() +
      scale_fill_viridis_d("Double vs\nsingle feeding\neffect, days",
                           guide = guide_bins(title.position = "right", reverse=TRUE)) +
      facet_wrap(~type, 
                 strip.position = "left", 
                 labeller = as_labeller(c(a = "k[mh]", b = "l[0]"),
                                        default = label_parsed)) +
      ylab(NULL) +
      xlab(TeX("$\\alpha$")) +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(size=14),
            axis.title.x = element_text(size=14))
    g
  }),
  tar_target(file_plot_sd_sensitivity_both, {
    filename <- "figures/single_double_feed_sensitivity_both.pdf"
    ggsave(filename, plot_sd_sensitivity_both, width=10, height=4);
    filename}, format="file"),
  
  # prior sensitivity analysis for fit
  tar_target(prior_sensitivities, prior_sensitivity(sampling_fit)),
  tar_target(files_prior_sensitivities, {
    write.csv(prior_sensitivities$mean, "data/processed/prior_sensitivities_mean.csv",
              row.names = FALSE)
    write.csv(prior_sensitivities$wasserstein, "data/processed/prior_sensitivities_wasserstein.csv",
              row.names = FALSE)
    "data/processed/prior_sensitivities_mean.csv"
  }),
  tar_target(graph_wasserstein_heatmap, wasserstein_heatmap(prior_sensitivities$wasserstein)),
  tar_target(file_graph_wasserstein_heatmap, {
    ggsave("figures/wasserstein_heatmap.pdf", graph_wasserstein_heatmap, width = 5, height = 3);
  }),
  tar_target(graph_estimates_heatmap, estimates_heatmap(prior_sensitivities$mean)),
  tar_target(graph_wasserstein_estimates, {
    plot_grid(graph_wasserstein_heatmap, graph_estimates_heatmap, nrow = 2,
              labels = c("A.", "B."))
  }),
  tar_target(file_graph_wasserstein_estimate, {
    ggsave("figures/wasserstein_estimates_heatmaps.pdf", graph_wasserstein_estimates)
    "figures/wasserstein_estimates_heatmaps.pdf"
  }),
  
  tar_target(posteriors_summary, posterior_summary(sampling_fit)),
  tar_target(file_posteriors_summary, {
    write.csv(posteriors_summary, "data/processed/posterior_summary.csv",
              row.names = FALSE)
  }),
  tar_target(graph_posterior_correlations, posteriors_correlation(sampling_fit)),
  tar_target(file_graph_posterior_correlations, {
    ggsave("figures/posterior_correlations.pdf", graph_posterior_correlations,
           width = 8, height = 4);
    "figures/posterior_correlations.pdf"
  }),
  
  # global sensitivity analysis
  tar_target(sobol_indices, sensitivity_sobol(sampling_fit)),
  tar_target(graph_sobol_indices, plot_sensitivity_sobol(sobol_indices)),
  tar_target(file_graph_sobol_indices, {
    ggsave("figures/sobol_indices.pdf", graph_sobol_indices,
           width = 8, height = 6);
    "figures/sobol_indices.pdf"
  }),
  
  # outputted parameter values for Alex
  tar_target(mean_parameter_values, get_summary_parameters(sampling_fit, list_stan_datasets$stan_data$x_0, 100)),
  tar_target(file_mean_parameter_values, {
    write.csv(mean_parameter_values, "data/processed/mean_parameter_values.csv");
    "data/processed/mean_parameter_values.csv"
  }
  )
)
