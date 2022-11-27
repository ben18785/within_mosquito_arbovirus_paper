# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan", "cowplot"), # packages that your targets need to run
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


list(
  
  # raw experimental data processing
  tar_target(filename_midgut, "data/raw/Compiled midgut data.xlsx",
             format = "file"),
  tar_target(filename_legs, "data/raw/Compiled leg data.xlsx",
             format = "file"),
  tar_target(df_midgut_legs,
             process_experimental_data(
               filename_midgut, filename_legs)),
  tar_target(graph_experimental_data, plot_experimental_data(df_midgut_legs)),
  tar_target(file_graph_experimental_data, {
    ggsave("figures/experimental_titers.pdf", graph_experimental_data, width=10, height=4);
    "figures/experimental_titers.pdf"}, format="file"),
  tar_target(filename_denv_dilutions, "data/raw/DENV test of dilutions 8-26-22.xlsx",
             format="file"),
  tar_target(df_denv_dilutions_infected,
             process_denv_dilutions_infected(filename_denv_dilutions)),
  tar_target(filename_disseminated_time_course, "data/raw/DENV SF DF dissemination time course 8-26-22.xlsx",
             format="file"),
  tar_target(df_disseminated_infection_time_course,
             process_disseminated_infection_time_course(
               filename_disseminated_time_course)),
  
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
  tar_target(stan_model, "src/stan/model_hurdle_binary_richard.stan",
             format="file"),
  tar_target(opt_fit,
             fit_optimise(
               list_stan_datasets$stan_data,
               stan_model, 2)),
  
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
    "figures/graph_single_double_chp.pdf"}, format="file")
)
