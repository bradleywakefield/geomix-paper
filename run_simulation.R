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

# --- Install geowarp if needed ---
if (!requireNamespace("geowarp", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("mbertolacci/geowarp")
}

library(geomix)
library(tidyverse)

# --- Step 1: Generate synthetic data ---
## Note although this script generates the synthetic data, the full MCMC sampling
## in Step 2 is computationally intensive. To save time, you can skip Step 1 and 
## load the pre-saved synthetic data instead (see Step 2 prompts).
# source("scripts/simulation/01_generate_data.R")  

# --- Step 2: Fit models and compute predictions ---
source("scripts/simulation/02_fit_models.R")

# --- Step 3: Produce results ---
source('scripts/simulation/03_results.R')
