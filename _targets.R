# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse", "rstan"), # packages that your targets need to run
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
  
  # stan data inputs processing
  tar_target(list_stan_datasets,
             prepare_stan_data_hurdle_denv_only(
               df_midgut_legs,
               df_denv_dilutions_infected,
               df_disseminated_infection_time_course)),
  
  # fit stan model via optimisation
  tar_target(stan_model, "src/stan/model_hurdle_binary_richard.stan",
             format="file"),
  tar_target(opt_fit,
             fit_optimise(
               list_stan_datasets$stan_data,
               stan_model, 5)),
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
    ggsave("figures/prevalence_midgut_only.pdf", graph_fit_prevalence_midgut_only, width=10, height=6);
    "figures/prevalence_midgut_only.pdf"}, format="file")
)
