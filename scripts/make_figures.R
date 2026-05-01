# Convenience script to regenerate all paper figures.
#
# Run all analysis scripts in sequence:
#
#   Rscript scripts/make_figures.R

scripts <- c(
  "scripts/simulation_study.R",
  "scripts/ijmuiden_ver.R"
)

for (script in scripts) {
  message("Running ", script, " ...")
  result <- tryCatch(
    source(script),
    error = function(e) {
      message("  WARNING: ", script, " failed with: ", conditionMessage(e))
      NULL
    }
  )
}

message("Done. Figures written to figures/.")
