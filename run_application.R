# GeoMix Paper — IJV Wind Farm Application
# ==========================================
#
# Reproduces all application results in the paper.
# Requires the proprietary IJV dataset (see data/README.md).
#
# PREREQUISITES:
#   Data files must be placed in data/application/ before running.
#   See data/README.md for the full list of required files.
#
# WORKFLOW:
#   Step 1: Process raw CPT and geophysics data.
#           Outputs: data/processed/data3D.RData
#                    data/processed/cell_grid.rds
#                    data/processed/locs.rds
#                    data/processed/line_df.rds
#
#   Step 2: Fit GeoMix, LGFM, and competing models on the IJV dataset;
#           generate all application figures, tables, and map outputs.
#           Outputs: results/application/, results/figures/ijv_*
#
# RUNTIME:
#   Running the full MCMC chains is computationally intensive and was originally
#   run on a multi-core server (~50 cores). The ask_run_mode() / ask_pred_mode()
#   prompts in Step 2 allow you to skip MCMC sampling and load pre-saved results.

# --- Install geomix if needed ---
if (!requireNamespace("geomix", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("bradleywakefield/geomix")
}

# --- Install geowarp if needed ---
if (!requireNamespace("geowarp", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("mbertolacci/geowarp")
}

library(geomix)
library(tidyverse)

# --- Step 1: Data processing ---
source("scripts/application/01_process_data.R")

# --- Step 2: Produce map and lines ---
source('scripts/application/02_map_figures.R')
rm(list = ls())

# --- Step 3: Fit models and run predictions ---
source("scripts/application/03_fit_models.R")
rm(list = ls())

# --- Step 4: Compute results ---
source('scripts/application/04_results.R')
