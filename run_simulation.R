# GeoMix Paper — Simulation Study
# =================================
#
# Reproduces all simulation study results in the paper.
# Does NOT require the IJV application data.
#
# PREREQUISITES:
#   Install the geomix package (done automatically below if needed).
#
# WORKFLOW:
#   Step 1: Generate synthetic datasets
#           Outputs: data/simulation/synthetic_data.rds
#                    data/simulation/offshore.rds
#
#   Step 2: Fit GeoMix, LGFM, and competing models on synthetic data;
#           generate all simulation figures and tables.
#           Outputs: results/simulation/,  results/figures/syn_*
#
# RUNTIME:
#   Running the full MCMC chains is computationally intensive.
#   The ask_run_mode() / ask_pred_mode() prompts in Step 2 let you
#   skip sampling and load pre-saved results instead.

# --- Install geomix if needed ---
if (!requireNamespace("geomix", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("bradleywakefield/geomix")
}
library(geomix)
library(tidyverse)

# --- Step 1: Generate synthetic data ---
source("scripts/simulation/01_generate_data.R")

# --- Step 2: Fit models and produce results ---
source("scripts/simulation/02_fit_models.R")
