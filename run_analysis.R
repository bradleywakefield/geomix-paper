# GeoMix Paper — Full Analysis Reproduction Script
# ==================================================
#
# Reproduces ALL results in the paper:
#   "GeoMix: A Bayesian Hierarchical Model for Joint Inversion of
#    Geotechnical and Geophysical Data"
#
# To run only the simulation study (no proprietary data required):
#   source("run_simulation.R")
#
# To run only the IJV application (requires proprietary data):
#   source("run_application.R")
#
# PREREQUISITES:
#   Application data files must be placed in data/application/ before running.
#   See data/README.md for the full list of required files.
#
# OPTIONAL SCRIPTS (run independently for additional diagnostics/figures):
#   scripts/exploratory/cross_section_plots.R  - Additional cross-section diagnostics
#   scripts/exploratory/misspec.R              - Misspecification illustration figure
#   scripts/exploratory/eda.R                  - Study area map and EDA figures

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

# --- Simulation study ---
source("run_simulation.R")

# --- IJV wind farm application ---
source("run_application.R")
