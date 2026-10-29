#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sensitivity of the IJV application to the Vecchia conditioning set size m
#
# In the grouped Vecchia approximation of the Z2 likelihood each group (CPT
# profile) conditions on the m nearest previously ordered observations. The
# paper uses m = 160. The grid concentrates on smaller m, where the
# approximation is coarser, with one larger value (m = 320).
#
# Refits GeoMix to the IJV data over a grid of m values, with everything else as
# in scripts/application/03_fit_models.R (beta fixed at its pseudo-likelihood
# estimate in results/application/beta.rds, m = 0.9, lateral length-scales
# fixed at their MAP estimates, cold-started chains, seed 16). The setup changes
# only in m and in the neighbour sets built from it (groupNeighbours); the
# predictions use the same m as the fit.
# The first 2,000 iterations of the existing chains in results/application/
# (m = 160) are used as the baseline comparator. Mirrors
# scripts/sensitivity/sensitivity_m.R.
#
# Usage in RStudio (project root as working directory, no interactive prompts):
#   Set source_task in section 1.0.2 and click Source.
#   "launch" starts one background R process per grid value and returns straight
#   away; the fits keep running if the session is closed. Progress is written to
#   results/sensitivity/m/logs/, and "list" shows the batches completed so far.
#
# Usage from a shell (project root):
#   Rscript scripts/sensitivity/sensitivity_m.R list        # show the grid and progress
#   Rscript scripts/sensitivity/sensitivity_m.R run <i>     # fit grid value i (n_chains cores)
#   Rscript scripts/sensitivity/sensitivity_m.R launch      # fit all grid values in parallel and wait
#                                                              #   (length(m_grid) * n_chains cores)
#   Rscript scripts/sensitivity/sensitivity_m.R summarise   # diagnostics, predictions, tables, figures
#
# "run" saves batches of 125 iterations; if it is interrupted, calling it again
# continues each chain from its last saved batch via run_chains(load_previous_state = TRUE).
# Continuing needs a geomix with the reload fix (bradleywakefield/geomix#7); "run"
# stops if the installed version refills alpha and gammaMat in the wrong order.
#
# Reads : data/processed/data3D.RData, results/application/MAP_covariance.rds,
#         results/application/beta.rds, results/application/GeoMix_*/ (baseline, read-only)
# Writes: results/sensitivity/m/
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 1.0.1 Settings ----
m_grid      <- c(10, 20, 40, 80, 120, 320)   # baseline 160 is the existing fit
n_chains    <- 4
n_iter      <- 2000
n_batches   <- 16        # 125 iterations per batch, as in the baseline run
burn_batches <- 4        # 500 burn-in iterations, as in the paper
thin        <- 10        # as in the paper
seed        <- 16        # as in the baseline run
lambda      <- 0.6       # Box-Cox parameter used in the application

base_path   <- "results/application"       # baseline chains and fixed estimates (read-only)
out_path    <- "results/sensitivity/m"

config_dir  <- function(m) file.path(out_path, sprintf("m_%03d", m))
keep_index  <- (burn_batches + 1):n_batches

### 1.0.2 What to do ----
# Used when the script is sourced (or run with no arguments):
#   "launch"    start a background fit for each grid value in source_index
#   "list"      show the grid and the batches completed so far
#   "summarise" diagnostics, predictions, tables, figures
source_task  <- "launch"
source_index <- seq_along(m_grid)   # grid values to launch, e.g. 1:4 or c(2, 5)

args <- commandArgs(trailingOnly = TRUE)
sourced <- length(args) == 0
mode <- if (sourced) source_task else args[1]
if (!mode %in% c("list", "run", "launch", "summarise") || (sourced && mode == "run")) {
  stop("Unknown mode '", mode, "'. Use one of: list, run <i> (shell only), launch, summarise.")
}

# Ends the script here: returns to the console when sourced, closes R under Rscript.
end_script <- function() if (interactive()) invokeRestart("abort") else quit(save = "no")

# A fit writes a RUNNING file into its folder while it is going.
marker_file <- function(m) file.path(config_dir(m), "RUNNING")
grid_status <- function() {
  data.frame(
    i = seq_along(m_grid), m = m_grid, path = sapply(m_grid, config_dir),
    batches_done = sapply(m_grid, function(b) min(sapply(seq_len(n_chains), function(ch) {
      length(list.files(file.path(config_dir(b), paste0("GeoMix_", ch)), pattern = "^batch_\\d+\\.rds$"))
    }))),
    running = file.exists(sapply(m_grid, marker_file)))
}

if (mode == "list") {
  print(grid_status(), row.names = FALSE)
  end_script()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Launch all grid values ------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Starts one "run <i>" process per grid value. From a shell it waits for all of
# them; when sourced it leaves them running in the background and returns.
if (mode == "launch") {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) == 0) script <- "scripts/sensitivity/sensitivity_m.R"
  if (!file.exists(script)) stop("Cannot find ", script, "; set the working directory to the project root.")
  dir.create(file.path(out_path, "logs"), recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), "Rscript")
  log_of <- function(i) file.path(out_path, "logs", sprintf("m_%03d.log", m_grid[i]))

  if (sourced) {
    status <- grid_status()
    status$action <- ifelse(!status$i %in% source_index, "not selected",
                     ifelse(status$running, "skipped: already running",
                     ifelse(status$batches_done >= n_batches, "skipped: complete", "started")))
    detach_cmd <- if (nzchar(Sys.which("setsid"))) c("setsid", "nohup") else "nohup"
    for (i in status$i[status$action == "started"]) {
      system2(detach_cmd[1], c(detach_cmd[-1], shQuote(rscript), shQuote(script), "run", i),
              stdout = log_of(i), stderr = log_of(i), wait = FALSE)
    }
    print(status[, c("i", "m", "batches_done", "action")], row.names = FALSE)
    message("Logs: ", file.path(out_path, "logs"), ". Set source_task to \"list\" to check progress.")
  } else {
    status <- parallel::mclapply(seq_along(m_grid), function(i) {
      system2(rscript, c(script, "run", i), stdout = log_of(i), stderr = log_of(i))
    }, mc.cores = length(m_grid), mc.preschedule = FALSE)
    print(data.frame(m = m_grid, exit_status = unlist(status)))
  }
  end_script()
}

library(geomix)
library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Fit one grid value ----------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (mode == "run") {
  i <- suppressWarnings(as.integer(args[2]))
  if (is.na(i) || !i %in% seq_along(m_grid)) {
    stop("Supply a grid index between 1 and ", length(m_grid), ", e.g. 'run 1'.")
  }
  m <- m_grid[i]
  path <- config_dir(m)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  message("m = ", m, " -> ", path)

  marker <- marker_file(m)
  if (file.exists(marker)) {
    stop("m = ", m, " looks to be running already. If that fit has died, delete ", marker, ".")
  }
  writeLines(c(paste("pid:", Sys.getpid()), paste("host:", Sys.info()[["nodename"]]),
               paste("started:", format(Sys.time()))), marker)
  reg.finalizer(globalenv(), function(e) unlink(marker), onexit = TRUE)

  ## 3.1 Model setup (as in scripts/application/03_fit_models.R, with m replaced) ----
  load("data/processed/data3D.RData")
  cov_par <- readRDS(file.path(base_path, "MAP_covariance.rds"))

  geomix_setup <- setupGeoMixModel(
    data,
    K = K,
    dims = dims,
    beta = readRDS(file.path(base_path, "beta.rds")),
    m = m,
    kappa = 0.9,
    aformula = ~1+d,
    variables = list(
      loc = "loc_id",
      xID = "xid",
      yID = "yid",
      dID = "dID",
      groups = "groups"),
    fix_lateral = T,
    inits = list(lL = cov_par$lL)
  )

  ## 3.2 Check the setup matches the baseline in everything except m ----
  baseline_file <- file.path(base_path, "GeoMix_1/geomix_setup.rds")
  if (file.exists(baseline_file)) {
    baseline_setup <- readRDS(baseline_file)
    same <- c(
      # groupNeighbours holds the conditioning sets, so it changes with m
      data_list = isTRUE(all.equal(geomix_setup$data_list[names(geomix_setup$data_list) != "groupNeighbours"],
                                   baseline_setup$data_list[names(baseline_setup$data_list) != "groupNeighbours"])),
      constants = isTRUE(all.equal(geomix_setup$constants[names(geomix_setup$constants) != "m"],
                                    baseline_setup$constants[names(baseline_setup$constants) != "m"])),
      controlHMC = isTRUE(all.equal(geomix_setup$controlHMC[names(geomix_setup$controlHMC) != "m"],
                                    baseline_setup$controlHMC[names(baseline_setup$controlHMC) != "m"])),
      # groupDeps/groupDepL list the groups whose conditioning sets reach each group, so they change with m
      controlGibbs = isTRUE(all.equal(
        geomix_setup$controlGibbs[!names(geomix_setup$controlGibbs) %in% c("m", "groupDeps", "groupDepL")],
        baseline_setup$controlGibbs[!names(baseline_setup$controlGibbs) %in% c("m", "groupDeps", "groupDepL")])),
      lL = isTRUE(all.equal(geomix_setup$inits$lL, baseline_setup$inits$lL))
    )
    if (!all(same)) {
      stop("Setup differs from the baseline in: ", paste(names(same)[!same], collapse = ", "))
    }
    message("Setup matches the baseline run in everything except m.")
    rm(baseline_setup)
  } else {
    warning("Baseline setup not found at ", baseline_file, "; equivalence with the baseline was not checked.")
  }

  ## 3.3 Decide whether to start or continue ----
  n_done <- sapply(seq_len(n_chains), function(ch) {
    length(list.files(file.path(path, paste0("GeoMix_", ch)), pattern = "^batch_\\d+\\.rds$"))
  })
  resume <- all(n_done > 0)
  batches_left <- if (resume) n_batches - min(n_done) else n_batches
  if (batches_left <= 0) {
    message("All ", n_batches, " batches already present for every chain; nothing to run.")
    end_script()
  }
  if (resume) {
    message("Continuing from saved batches (completed per chain: ", paste(n_done, collapse = ", "), ").")
    if (any(grepl("byrow = TRUE", deparse(run_chains), fixed = TRUE))) {
      stop("The installed geomix reloads alpha and gammaMat in the wrong order; install a version with bradleywakefield/geomix#7.")
    }
  }

  writeLines(c(
    paste("m:", m), paste("started:", format(Sys.time())),
    paste("resumed:", resume), paste("geomix:", as.character(packageVersion("geomix"))),
    paste("nimble:", as.character(packageVersion("nimble"))), R.version.string,
    paste("host:", Sys.info()[["nodename"]])
  ), file.path(path, paste0("run_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")))

  ## 3.4 Run chains (parallel, cold start unless continuing) ----
  run_chains(geomix_setup,
             nchains = n_chains,
             path = path,
             controlMCMC = list(niter = batches_left * (n_iter / n_batches), thin = 1,
                                nbatches = batches_left, save_batches = T,
                                retain_draws = F),
             LGFM = F,
             run_parallel = TRUE,
             load_previous_state = resume,
             mc.cores = n_chains,
             seed = seed)
  message("Finished m = ", m, " at ", format(Sys.time()))
  end_script()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Summarise -------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(patchwork)
library(posterior)
library(scoringRules)

n_cores <- max(1, parallel::detectCores() - 2)
sum_path <- file.path(out_path, "summary")
dir.create(file.path(sum_path, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sum_path, "tables"), showWarnings = FALSE)
dir.create(file.path(sum_path, "predictions"), showWarnings = FALSE)

## 4.1 Scoring (as in scripts/utils/prediction_scoring.R) ----
box_cox <- function(y, lambda = 0.6) (y^lambda - 1) / lambda

dss_p_score <- function(observed, posterior_samples, p) {
  sapply(seq_len(floor(length(observed) / p)), function(i) {
    indices <- p * (i - 1) + seq_len(p)
    samples <- posterior_samples[indices, ]
    mu <- rowMeans(samples)
    Sigma <- cov(t(samples)) + diag(length(mu)) * 1e-6
    chol_Sigma <- chol(Sigma)
    sum(log(diag(chol_Sigma))) +
      0.5 * sum(backsolve(chol_Sigma, observed[indices] - mu, transpose = TRUE)^2)
  })
}

# MSPE is the quantity reported in the "RMSPE" column of the paper's tables
# (prediction_scoring.R returns the mean squared error); RMSPE is its square root.
score_predictions <- function(obs, samples, loc_id) {
  samples[is.nan(samples)] <- 0
  pred_mean <- rowMeans(samples)
  lci <- apply(samples, 1, quantile, 0.025)
  uci <- apply(samples, 1, quantile, 0.975)
  dss2 <- mean(unlist(lapply(unique(loc_id), function(i) {
    indices <- which(loc_id == i)
    dss_p_score(obs[indices], samples[indices, , drop = FALSE], p = 2)
  })))
  data.frame(
    MSPE     = mean((pred_mean - obs)^2),
    RMSPE    = sqrt(mean((pred_mean - obs)^2)),
    coverage = mean(obs > lci & obs < uci),
    R2       = 1 - sum((obs - pred_mean)^2) / sum((obs - mean(obs))^2),
    bias     = mean(pred_mean - obs),
    IS       = mean(ints_sample(y = obs, dat = samples, 0.95)),
    CRPS     = mean(crps_sample(y = obs, dat = samples)),
    DSS      = mean(dss_sample(y = obs, dat = samples)),
    DSS2     = dss2,
    logS     = mean(logs_sample(y = obs, dat = samples))
  )
}

## 4.2 Configurations to summarise ----
configs <- tibble(m = m_grid, path = sapply(m_grid, config_dir), baseline = FALSE)
if (dir.exists(file.path(base_path, "GeoMix_1"))) {
  configs <- bind_rows(configs, tibble(m = readRDS(file.path(base_path, "GeoMix_1/geomix_setup.rds"))$constants$m,
                                       path = base_path, baseline = TRUE))
} else {
  warning("Baseline chains not found in ", base_path, "; summarising the grid without a baseline.")
}
configs <- arrange(configs, m) %>% mutate(label = sprintf("%d", as.integer(m)))

# Reuse the baseline predictive draws from the kappa study, so the baseline row
# matches the kappa (and first beta) tables exactly
base_pred <- file.path(sum_path, "predictions", "m_160.rds")
kappa_base_pred <- "results/sensitivity/kappa/summary/predictions/kappa_0.900.rds"
if (!file.exists(base_pred) && file.exists(kappa_base_pred)) file.copy(kappa_base_pred, base_pred)

batch_files <- function(path, ch) {
  f <- list.files(file.path(path, paste0("GeoMix_", ch)), pattern = "^batch_\\d+\\.rds$", full.names = TRUE)
  f[order(as.numeric(str_extract(basename(f), "\\d+")))]
}
configs$complete <- sapply(configs$path, function(p) {
  all(sapply(seq_len(n_chains), function(ch) length(batch_files(p, ch)) >= n_batches))
})
if (any(!configs$complete)) {
  warning("Skipping incomplete configurations (fewer than ", n_batches, " batches in some chain): m = ",
          paste(configs$label[!configs$complete], collapse = ", "))
}
configs <- filter(configs, complete)
if (nrow(configs) == 0) stop("No complete configurations to summarise.")

## 4.3 Loop over configurations ----
trace_list <- time_list <- diag_list <- overall_list <- metric_list <- list()
prob_list <- list()

for (r in seq_len(nrow(configs))) {
  cfg <- configs[r, ]
  message("\n=== m = ", cfg$label, if (cfg$baseline) " (baseline)" else "", " ===")

  ### 4.3.1 Log-likelihood traces and run time (all iterations, including burn-in) ----
  for (ch in seq_len(n_chains)) {
    files <- batch_files(cfg$path, ch)[seq_len(n_batches)]
    lp <- do.call(rbind, lapply(files, function(f) {
      readRDS(f)[, c("logProb_Z1[1]", "logProb_Z2[1]"), drop = FALSE]
    }))
    trace_list[[length(trace_list) + 1]] <- tibble(
      m = cfg$m, label = cfg$label, chain = ch, iteration = seq_len(nrow(lp)),
      logProb_Z1 = lp[, 1], logProb_Z2 = lp[, 2])
    time_list[[length(time_list) + 1]] <- tibble(
      m = cfg$m, label = cfg$label, chain = ch,
      hours_per_batch = median(as.numeric(diff(file.mtime(files)), units = "hours")))
  }

  ### 4.3.2 Load samples and extract parameters ----
  geomix_setup <- readRDS(file.path(cfg$path, "GeoMix_1/geomix_setup.rds"))
  samples_post <- load_mcmc_samples(cfg$path, index = keep_index, thin = thin)
  params_post <- extract_parameters(samples_post)

  ### 4.3.3 Diagnostics ----
  diagnostics <- run_mcmc_diagnostics(
    exclude = c("sigma2_L", paste0("lL[", 1:8, "]")),
    params_post, Y1index = geomix_setup$controlGibbs$Z2_ind)
  diag_list[[r]] <- cbind(m = cfg$m, label = cfg$label, diagnostics$tables$diagnostics)
  overall_list[[r]] <- cbind(m = cfg$m, label = cfg$label, diagnostics$tables$overall_diagnostics)

  ### 4.3.4 Latent strata: posterior class probabilities ----
  prob_list[[cfg$label]] <- Reduce(`+`, lapply(params_post, function(x) x$params$Y1prob)) / length(params_post)

  ### 4.3.5 Test-set predictions ----
  predict_index <- which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
  pred_file <- file.path(sum_path, "predictions", paste0("m_", cfg$label, ".rds"))
  if (file.exists(pred_file)) {
    test_pred <- readRDS(pred_file)
  } else {
    set.seed(seed)
    test_pred <- produce_prediction(
      samples_post,
      geomix_setup,
      nugget = T,
      include_samples = T,
      run_parallel = TRUE,
      mc.cores = n_cores,
      predict_index = predict_index
    )
    saveRDS(test_pred, pred_file)
  }
  test_obs <- box_cox(geomix_setup$df$qc[predict_index], lambda)
  test_loc <- geomix_setup$df$loc_id[predict_index]
  chain_of_col <- rep(seq_along(samples_post), sapply(samples_post, nrow))
  metric_list[[r]] <- bind_rows(
    cbind(m = cfg$m, label = cfg$label, chain = "all",
          score_predictions(test_obs, test_pred$samples, test_loc)),
    map_dfr(seq_along(samples_post), function(ch) {
      cbind(m = cfg$m, label = cfg$label, chain = as.character(ch),
            score_predictions(test_obs, test_pred$samples[, chain_of_col == ch, drop = FALSE], test_loc))
    })
  )

  Z1 <- geomix_setup$data_list$Z1
  Z2_ind <- geomix_setup$constants$Z2_ind
  rm(samples_post, params_post, diagnostics, test_pred); gc()
}

## 4.4 Latent strata summaries ----
N1 <- length(Z1)
is_cpt <- seq_len(N1) %in% Z2_ind
site_sets <- list("All lattice sites" = rep(TRUE, N1), "Training CPT sites" = is_cpt, "Other sites" = !is_cpt)
base_label <- configs$label[configs$baseline]
strata_summary <- map_dfr(configs$label, function(lab) {
  prob <- prob_list[[lab]]
  map_class <- max.col(prob, ties.method = "first")
  p_not_Z1 <- 1 - prob[cbind(seq_len(N1), Z1)]
  map_dfr(names(site_sets), function(s) {
    idx <- site_sets[[s]]
    out <- tibble(label = lab, sites = s,
                  p_differs_from_Z1 = mean(p_not_Z1[idx]),
                  share_map_differs_from_Z1 = mean((map_class != Z1)[idx]),
                  mean_max_prob = mean(apply(prob[idx, , drop = FALSE], 1, max)))
    if (length(base_label) == 1) {
      base_prob <- prob_list[[base_label]]
      out$share_map_differs_from_baseline <- mean((map_class != max.col(base_prob, ties.method = "first"))[idx])
      out$mean_tv_distance_from_baseline <- mean(0.5 * rowSums(abs(prob - base_prob))[idx])
    }
    out
  })
}) %>%
  left_join(select(configs, label, m, baseline), by = "label") %>%
  relocate(m, label, baseline)

## 4.5 Save tables ----
trace_df   <- bind_rows(trace_list)
time_df    <- bind_rows(time_list)
diag_df    <- bind_rows(diag_list)
overall_df <- bind_rows(overall_list)
metric_df  <- bind_rows(metric_list)

write.csv(metric_df,      file.path(sum_path, "tables", "test_metrics.csv"), row.names = FALSE)
write.csv(strata_summary, file.path(sum_path, "tables", "strata_summary.csv"), row.names = FALSE)
write.csv(diag_df,        file.path(sum_path, "tables", "parameter_diagnostics.csv"), row.names = FALSE)
write.csv(overall_df,     file.path(sum_path, "tables", "overall_diagnostics.csv"), row.names = FALSE)
write.csv(time_df,        file.path(sum_path, "tables", "run_time.csv"), row.names = FALSE)
saveRDS(trace_df,         file.path(sum_path, "tables", "logProb_trace.rds"))
saveRDS(prob_list,        file.path(sum_path, "predictions", "Y1_probabilities.rds"))

## 4.6 Figures ----
plot_theme <- theme_bw() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        axis.text = element_text(size = 8), axis.title = element_text(size = 10),
        strip.text = element_text(size = 8))
burn_iter <- burn_batches * n_iter / n_batches
m_hat <- configs$m[configs$baseline]

### 4.6.1 Log-likelihood traces ----
trace_plot <- function(var, ylab, drop_first = 0) {
  trace_df %>%
    filter(iteration > drop_first) %>%
    ggplot() +
    geom_line(aes(x = iteration, y = .data[[var]], col = factor(chain)), linewidth = 0.2, alpha = 0.7) +
    geom_vline(xintercept = burn_iter, linetype = "dashed") +
    facet_wrap(vars(paste0("m = ", label)), scales = "free_y", ncol = 3) +
    labs(x = "Iteration", y = ylab, col = "Chain") + plot_theme
}
ggsave(file.path(sum_path, "figures", "trace_logProb_Z2.pdf"), trace_plot("logProb_Z2", "Z2 log-likelihood"),
       width = 8, height = 7)
ggsave(file.path(sum_path, "figures", "trace_logProb_Z2_after100.pdf"),
       trace_plot("logProb_Z2", "Z2 log-likelihood", drop_first = 100), width = 8, height = 7)
ggsave(file.path(sum_path, "figures", "trace_logProb_Z1.pdf"), trace_plot("logProb_Z1", "Z1 log-likelihood"),
       width = 8, height = 7)
ggsave(file.path(sum_path, "figures", "trace_logProb_Z1_after100.pdf"),
       trace_plot("logProb_Z1", "Z1 log-likelihood", drop_first = 100), width = 8, height = 7)

### 4.6.2 Test-set metrics against m ----
metric_long <- metric_df %>%
  pivot_longer(-c(m, label, chain), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = unique(metric)))
g_metrics <- ggplot() +
  geom_vline(xintercept = m_hat, linetype = "dashed") +
  geom_point(data = filter(metric_long, chain != "all"), aes(x = m, y = value), size = 0.6, alpha = 0.5) +
  geom_line(data = filter(metric_long, chain == "all"), aes(x = m, y = value)) +
  geom_point(data = filter(metric_long, chain == "all"), aes(x = m, y = value), col = "red", size = 1) +
  facet_wrap(vars(metric), scales = "free_y", ncol = 5) +
  scale_x_continuous(trans = "log2", breaks = c(m_grid, m_hat)) +
  labs(x = expression(m), y = "") + plot_theme
ggsave(file.path(sum_path, "figures", "metrics_vs_m.pdf"), g_metrics, width = 10, height = 4.5)

### 4.6.3 Parameter posteriors against m ----
g_params <- diag_df %>%
  mutate(type = str_remove(param, "\\[.*\\]")) %>%
  filter(type %in% c("sigma2", "lD", "tau2", "h")) %>%
  ggplot(aes(x = m, y = q50, ymin = q025, ymax = q975)) +
  geom_vline(xintercept = m_hat, linetype = "dashed") +
  geom_pointrange(size = 0.15) +
  facet_wrap(vars(param), scales = "free_y", ncol = 6) +
  scale_x_continuous(trans = "log2", breaks = c(m_grid, m_hat)) +
  labs(x = expression(m), y = "Posterior median and 95% interval") + plot_theme
ggsave(file.path(sum_path, "figures", "parameters_vs_m.pdf"), g_params, width = 11, height = 6)

### 4.6.4 Latent strata against m ----
g_strata <- strata_summary %>%
  pivot_longer(-c(m, label, baseline, sites), names_to = "summary", values_to = "value") %>%
  ggplot(aes(x = m, y = value, col = sites)) +
  geom_vline(xintercept = m_hat, linetype = "dashed") +
  geom_line() + geom_point(size = 0.8) +
  facet_wrap(vars(summary), scales = "free_y", ncol = 3) +
  scale_x_continuous(trans = "log2", breaks = c(m_grid, m_hat)) +
  labs(x = expression(m), y = "", col = "") + plot_theme
ggsave(file.path(sum_path, "figures", "strata_vs_m.pdf"), g_strata, width = 9, height = 5.5)

## 4.7 Report ----
cat("\n=== Kappa sensitivity ===\n")
cat("\nTest-set metrics (all chains):\n")
print(filter(metric_df, chain == "all") %>% select(-chain), digits = 4)
cat("\nBetween-chain sd of each metric:\n")
print(metric_df %>% filter(chain != "all") %>% group_by(label) %>%
        summarise(across(MSPE:logS, sd)) %>% as.data.frame(), digits = 3)
cat("\nConvergence:\n")
print(overall_df, digits = 4)
cat("\nLatent strata:\n")
print(as.data.frame(strata_summary), digits = 3)
cat("\nMedian hours per batch of ", n_iter / n_batches, " iterations:\n", sep = "")
print(time_df %>% group_by(label) %>% summarise(hours_per_batch = median(hours_per_batch)) %>% as.data.frame())
