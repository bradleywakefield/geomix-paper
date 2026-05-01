# GeoMix Paper - Analysis Reproduction Script
# ============================================
#
# This script reproduces all results in the paper:
#   "GeoMix: A Bayesian Hierarchical Model for Joint Inversion of
#    Geotechnical and Geophysical Data"
#
# PREREQUISITES:
#   Data files must be placed in data/ and data/geophys/ before running.
#   See data/README.md for the full list of required files.
#
# WORKFLOW:
#   Step 1: Process raw CPT and geophysics data
#           Outputs: data/processed/data3D.RData, data/cell_grid.rds, data/locs.rds
#
#   Step 2: Fit models and generate all paper results
#           Includes simulation study (Section X) and IJV wind farm application (Section Y)
#           Outputs: results/simulation/, results/application/, results/figures/
#
# RUNTIME:
#   Running the full MCMC chains is computationally intensive and was originally
#   run on a multi-core server (~50 cores). The confirm_run() prompts in
#   02_fit_models.R allow you to skip MCMC sampling and load pre-saved results.
#
# OPTIONAL SCRIPTS (run independently after Step 2):
#   scripts/exploratory/eda.R           - Map figures (requires geomix_setup in memory)
#   scripts/exploratory/misspec.R       - Misspecification illustration figure
#   scripts/exploratory/model_plot.R    - Additional cross-section diagnostics

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
