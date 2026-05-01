# GeoMix Paper - Analysis Reproduction Script
# ============================================
#
# This script reproduces all results in the paper:
#   "GeoMix: A Bayesian Hierarchical Model for Joint Inversion of
#    Geotechnical and Geophysical Data"
#
# PREREQUISITES:
#   Data files must be placed in data/application/ and data/simulation/ before running.
#   See data/README.md for the full list of required files.
#
# WORKFLOW:
#   Step 1: Process raw CPT and geophysics data
#           Outputs: data/processed/data3D.RData, data/processed/cell_grid.rds, data/processed/locs.rds
#
#   Step 2: Fit models and generate all paper results
#           Includes simulation study and IJV wind farm application
#           Internally runs:
#             Step 3 - simulation study figures (03_simulation_results.R)
#             Step 4 - IJV application model fitting
#             Step 5 - study area map figures (05_map_figures.R)
#                      saves data/processed/line_df.rds and data/processed/cell_grid_map.rds
#             Step 6 - IJV application results figures (04_application_results.R)
#           Outputs: results/simulation/, results/application/, results/figures/
#
# RUNTIME:
#   Running the full MCMC chains is computationally intensive and was originally
#   run on a multi-core server (~50 cores). The confirm_run() prompts in
#   02_fit_models.R allow you to skip MCMC sampling and load pre-saved results.
#
# OPTIONAL SCRIPTS (run independently for additional diagnostics):
#   scripts/exploratory/cross_section_plots.R  - Additional cross-section diagnostics
#   scripts/exploratory/misspec.R              - Misspecification illustration figure

# --- Install geomix if needed ---
if (!requireNamespace("geomix", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("bradleywakefield/geomix")
}
library(geomix)
library(tidyverse)

# --- Step 1: Data processing ---
source("scripts/01_process_data.R")

# --- Step 2: Model fitting and results ---
# This sources 03_simulation_results.R and 04_application_results.R internally
source("scripts/02_fit_models.R")
