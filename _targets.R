# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse"), # packages that your targets need to run
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

# Replace the target list below with your own:
list(
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
    "figures/experimental_titers.pdf"}, format="file")
)
