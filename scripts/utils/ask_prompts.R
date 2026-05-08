#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helpers for 02_fit_models.R ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## ask_run_mode: prompt for MCMC chain run mode ----
## Returns "s" (scratch), "c" (continue), or "l" (load existing samples)
ask_run_mode <- function(label) {
  cat(paste0("\n--- ", label, " ---\n"))
  cat("  [s] Run chains from scratch\n")
  cat("  [c] Continue from previous saved state\n")
  cat("  [l] Skip — load existing samples only\n")
  ans <- tolower(trimws(readline("Choice [s/c/l]: ")))
  while (!ans %in% c("s", "c", "l")) {
    message("  Please enter s, c, or l.")
    ans <- tolower(trimws(readline("Choice [s/c/l]: ")))
  }
  ans
}

## ask_pred_mode: prompt for prediction / estimation run mode ----
## Returns "r" (run) or "l" (load existing output)
ask_pred_mode <- function(label) {
  cat(paste0("\n--- ", label, " ---\n"))
  cat("  [r] Run\n")
  cat("  [l] Load existing output\n")
  ans <- tolower(trimws(readline("Choice [r/l]: ")))
  while (!ans %in% c("r", "l")) {
    message("  Please enter r or l.")
    ans <- tolower(trimws(readline("Choice [r/l]: ")))
  }
  ans
}

## save_diag_tables: save MCMC diagnostic tables to {path}/tables/ ----
save_diag_tables <- function(diag, path, prefix) {
  dir.create(file.path(path, "tables"), showWarnings = FALSE)
  write.csv(diag$tables$overall_diagnostics,
            file.path(path, "tables", paste0(prefix, "_overall.csv")),    row.names = FALSE)
  write.csv(diag$tables$class_diagnostics,
            file.path(path, "tables", paste0(prefix, "_class.csv")),      row.names = FALSE)
  write.csv(diag$tables$diagnostics,
            file.path(path, "tables", paste0(prefix, "_params.csv")),     row.names = FALSE)
  write.csv(diag$tables$worst_rhat,
            file.path(path, "tables", paste0(prefix, "_worst_rhat.csv")), row.names = FALSE)
  write.csv(diag$tables$worst_ess,
            file.path(path, "tables", paste0(prefix, "_worst_ess.csv")),  row.names = FALSE)
}
