#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exact-likelihood check of the grouped Vecchia approximation (IJV application)
#
# Uses the existing GeoMix chains in results/application/ (no new MCMC).
# For every retained posterior draw this script
#   (i)   recomputes the grouped Vecchia Z2 log-likelihood with the package
#         function and checks it against the value saved by the sampler,
#   (ii)  computes the exact Z2 log-likelihood (dense per-class covariances),
#   (iii) forms importance weights exact / Vecchia and, only if these pass the
#         reliability check in Section 5, reweights parameters, latent strata
#         and test-set predictions to the exact-likelihood posterior.
#
# Usage (from the project root):
#   Rscript scripts/sensitivity/exact_likelihood_check.R
#
# Reads : results/application/GeoMix_*/   (read-only)
# Writes: results/sensitivity/exact_likelihood/
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 1.0.1 Libraries ----
library(geomix)
library(tidyverse)
library(patchwork)
library(posterior)
library(scoringRules)

### 1.0.2 Settings ----
base_path   <- "results/application"                     # existing chains (read-only)
out_path    <- "results/sensitivity/exact_likelihood"
batch_index <- 5:32      # batches retained in the paper (500 burn-in iterations discarded)
thin        <- 10        # as in the paper
n_chains    <- 4
n_cores     <- max(1, parallel::detectCores() - 2)
n_resample  <- 5000      # size of the importance resample used for reweighted summaries
n_full_check <- 2        # draws on which the dense likelihood is checked against the package
seed        <- 16
lambda      <- 0.6       # Box-Cox parameter used in the application

dir.create(file.path(out_path, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_path, "tables"), showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Helper functions ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Low-memory batch reader ----
# Returns the same rows as load_mcmc_samples(path, index = index, thin = thin)
# for one chain, but thins batch by batch so the full unthinned matrix is never
# held in memory.
# The MCMC iteration of each retained row is attached as attribute "iteration"
# (assuming equal batch sizes before the first batch read).
read_thinned_chain <- function(dir, index, thin) {
  offset <- 0
  out <- vector("list", length(index))
  iteration <- vector("list", length(index))
  for (i in seq_along(index)) {
    batch <- readRDS(file.path(dir, paste0("batch_", index[i], ".rds")))
    keep <- which((offset + seq_len(nrow(batch))) %% thin == 0)
    out[[i]] <- batch[keep, , drop = FALSE]
    iteration[[i]] <- (index[1] - 1) * nrow(batch) + offset + keep
    offset <- offset + nrow(batch)
  }
  out <- do.call(rbind, out)
  attr(out, "iteration") <- unlist(iteration)
  out
}

## 2.2 Exact Z2 log-likelihood ----
# Given the latent strata, Z2 is independent across classes and Gaussian within
# class with a Matern-3/2 covariance plus nugget. Returns one term per class.
exact_loglik <- function(Y1_valid, alpha, sigma2, tau2, lL, lD, st) {
  out <- numeric(st$K)
  for (k in seq_len(st$K)) {
    idx <- which(Y1_valid == k)
    n <- length(idx)
    if (n == 0) next
    r2 <- st$distD[st$dIDv[idx], st$dIDv[idx], drop = FALSE] / lD[k]^2
    if (st$LFlag[k] != 0) {
      r2 <- r2 + st$distL[st$locIDv[idx], st$locIDv[idx], drop = FALSE] / lL[k]^2
    }
    r <- sqrt(3 * r2)
    C <- sigma2[k] * (1 + r) * exp(-r)
    diag(C) <- diag(C) + tau2
    U <- chol(C)
    res <- st$Z2[idx] - drop(st$X[idx, , drop = FALSE] %*% alpha[k, ])
    v <- backsolve(U, res, transpose = TRUE)
    out[k] <- -sum(log(diag(U))) - 0.5 * sum(v^2) - 0.5 * n * log(2 * pi)
  }
  out
}

## 2.3 Systematic resampling ----
systematic_resample <- function(w, n) {
  u <- (runif(1) + 0:(n - 1)) / n
  pmin(findInterval(u, cumsum(w)) + 1, length(w))
}

## 2.4 Scoring (as in scripts/utils/prediction_scoring.R) ----
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Load setup and samples ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geomix_setup <- readRDS(file.path(base_path, "GeoMix_1/geomix_setup.rds"))
dl <- geomix_setup$data_list
cs <- geomix_setup$constants

st <- list(
  K = cs$K, Z2 = dl$Z2, X = dl$X, Z2_ind = cs$Z2_ind,
  dID = dl$dID, locID = dl$locID,
  dIDv = dl$dID[cs$Z2_ind], locIDv = dl$locID[cs$Z2_ind],
  distD = dl$distD, distL = dl$distL, LFlag = dl$LFlag
)

message("Loading thinned samples...")
samples_post <- lapply(seq_len(n_chains), function(ch) {
  read_thinned_chain(file.path(base_path, paste0("GeoMix_", ch)), batch_index, thin)
})
n_draws <- sapply(samples_post, nrow)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Vecchia and exact log-likelihood per draw ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loglik_file <- file.path(out_path, "loglik_draws.rds")
if (file.exists(loglik_file)) {
  message("Loading cached log-likelihood evaluations from ", loglik_file)
  loglik <- readRDS(loglik_file)
} else {
  CevaluateGroupGPvec <- nimble::compileNimble(evaluateGroupGPvec)

  vecchia_loglik <- function(Y1, alpha, sigma2, tau2, lL, lD, vecchia) {
    colSums(CevaluateGroupGPvec(
      Z2 = st$Z2, alpha = alpha, sigma2 = sigma2, tau2 = tau2, lL = lL, lD = lD,
      LFlag = st$LFlag, Y1 = Y1, X = st$X, K = st$K, Z2_ind = st$Z2_ind,
      dID = st$dID, locID = st$locID, distD = st$distD, distL = st$distL,
      m = vecchia$m, groupLookup = vecchia$groupLookup,
      groupNum = vecchia$groupNum, groupNeighbours = vecchia$groupNeighbours
    ))
  }
  vecchia_fit <- list(m = cs$m, groupLookup = dl$groupLookup,
                      groupNum = dl$groupNum, groupNeighbours = dl$groupNeighbours)

  ## 4.1 Evaluate both likelihoods at every retained draw ----
  loglik <- list()
  for (ch in seq_len(n_chains)) {
    message("Chain ", ch, ": evaluating ", n_draws[ch], " draws on ", n_cores, " cores...")
    par <- extract_parameters(samples_post[[ch]])$samples
    res <- parallel::mclapply(seq_len(n_draws[ch]), function(j) {
      vec <- vecchia_loglik(par$Y1[j, ], par$alpha[j, , ], par$sigma2[j, ], par$tau2[j],
                            par$lL[j, ], par$lD[j, ], vecchia_fit)
      exa <- exact_loglik(par$Y1[j, st$Z2_ind], par$alpha[j, , ], par$sigma2[j, ], par$tau2[j],
                          par$lL[j, ], par$lD[j, ], st)
      list(vec = vec, exa = exa)
    }, mc.cores = n_cores)
    failed <- !sapply(res, is.list)
    if (any(failed)) stop("Log-likelihood evaluation failed for ", sum(failed), " draws in chain ", ch)
    loglik[[ch]] <- list(
      chain = ch,
      iteration = attr(samples_post[[ch]], "iteration"),
      saved   = par$logProbZ2,
      vecchia = t(sapply(res, `[[`, "vec")),
      exact   = t(sapply(res, `[[`, "exa"))
    )
    rm(par, res); gc()
  }

  ## 4.2 Check the dense likelihood against the package's own exact form ----
  # A grouped Vecchia approximation in which every group conditions on all
  # earlier groups, with no cap on the conditioning set, is the exact likelihood.
  message("Checking dense likelihood against the full-dependency grouped form...")
  vecchia_full <- setupVecchiaGeoMix(
    geomix_setup$lattice_coords[cs$Z2_ind, ], m = cs$N2,
    groups = geomix_setup$vecchia$groups, depStructure = "full"
  )
  vecchia_full <- c(list(m = cs$N2), vecchia_full[c("groupLookup", "groupNum", "groupNeighbours")])
  par <- extract_parameters(samples_post[[1]][seq_len(n_full_check), , drop = FALSE])$samples
  full_check <- do.call(rbind, lapply(seq_len(n_full_check), function(j) {
    full <- sum(vecchia_loglik(par$Y1[j, ], par$alpha[j, , ], par$sigma2[j, ], par$tau2[j],
                               par$lL[j, ], par$lD[j, ], vecchia_full))
    data.frame(draw = j, dense = sum(loglik[[1]]$exact[j, ]), full_dependency = full)
  }))
  full_check$abs_diff <- abs(full_check$dense - full_check$full_dependency)
  attr(loglik, "full_check") <- full_check
  rm(par)

  saveRDS(loglik, loglik_file)
}

## 4.3 Tidy and check ----
loglik_df <- map_dfr(loglik, function(x) {
  tibble(chain = x$chain, iteration = x$iteration, saved = x$saved,
         vecchia = rowSums(x$vecchia), exact = rowSums(x$exact))
}) %>%
  mutate(draw = row_number(), gap = exact - vecchia)

saved_check <- max(abs(loglik_df$saved - loglik_df$vecchia))
full_check  <- attr(loglik, "full_check")
message(sprintf("Max |saved - recomputed| Vecchia log-likelihood: %.2e", saved_check))
message(sprintf("Max |dense - full-dependency| exact log-likelihood: %.2e", max(full_check$abs_diff)))
if (saved_check > 1e-4) warning("Recomputed Vecchia log-likelihood does not match the saved values.")
if (max(full_check$abs_diff) > 1e-4) warning("Dense exact log-likelihood does not match the full-dependency form.")

class_gap <- map_dfr(loglik, function(x) {
  as_tibble(x$exact - x$vecchia, .name_repair = ~ paste0("class_", seq_len(st$K)))
}) %>%
  pivot_longer(everything(), names_to = "class", values_to = "gap") %>%
  group_by(class) %>%
  summarise(mean_gap = mean(gap), sd_gap = sd(gap), min_gap = min(gap), max_gap = max(gap))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Importance weights ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_w <- loglik_df$gap
w <- exp(log_w - max(log_w))
w <- w / sum(w)
loglik_df$weight <- w

khat <- tryCatch(
  as.numeric(unlist(posterior::pareto_khat(log_w, are_log_weights = TRUE))[1]),
  error = function(e) NA_real_
)
weight_diag <- tibble(
  n_draws    = length(w),
  mean_gap   = mean(log_w),
  sd_gap     = sd(log_w),
  cor_exact_vecchia = cor(loglik_df$exact, loglik_df$vecchia),
  ess        = 1 / sum(w^2),
  max_weight = max(w),
  pareto_k   = khat
) %>%
  mutate(reliable = ess >= 100 & (is.na(pareto_k) | pareto_k <= 0.7))

chain_diag <- loglik_df %>%
  group_by(chain) %>%
  summarise(mean_gap = mean(gap), sd_gap = sd(gap),
            ess_within = {ww <- exp(gap - max(gap)); sum(ww)^2 / sum(ww^2)},
            share_of_total_weight = sum(weight))

set.seed(seed)
resample_idx <- systematic_resample(w, n_resample)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Reweighted posterior summaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 Parameters ----
param_draws <- do.call(rbind, lapply(samples_post, function(x) {
  x[, !grepl("^Y1\\[|^logProb|^lL\\[|^sigma2_L$", colnames(x)), drop = FALSE]
}))
summarise_draws_mat <- function(x, label) {
  tibble(parameter = colnames(x), posterior = label,
         mean = colMeans(x), sd = apply(x, 2, sd),
         q025 = apply(x, 2, quantile, 0.025), q50 = apply(x, 2, quantile, 0.5),
         q975 = apply(x, 2, quantile, 0.975))
}
param_summary <- summarise_draws_mat(param_draws, "Vecchia")
if (weight_diag$reliable) {
  param_summary <- bind_rows(
    param_summary,
    summarise_draws_mat(param_draws[resample_idx, , drop = FALSE], "Exact (reweighted)")
  )
}

## 6.2 Latent strata ----
# Posterior class probabilities at every lattice site, unweighted and weighted.
chain_of_draw <- rep(seq_len(n_chains), n_draws)
Y1_cols <- grep("^Y1\\[", colnames(samples_post[[1]]))
Y1_site <- as.integer(sub("^Y1\\[(\\d+)\\]$", "\\1", colnames(samples_post[[1]])[Y1_cols]))
prob_vecchia <- prob_exact <- matrix(0, nrow = cs$N1, ncol = st$K)
for (ch in seq_len(n_chains)) {
  Y1_ch <- samples_post[[ch]][, Y1_cols[order(Y1_site)], drop = FALSE]
  w_ch <- w[chain_of_draw == ch]
  for (k in seq_len(st$K)) {
    ind <- Y1_ch == k
    prob_vecchia[, k] <- prob_vecchia[, k] + colSums(ind) / length(w)
    prob_exact[, k]   <- prob_exact[, k] + drop(crossprod(w_ch, ind))
  }
  rm(Y1_ch, ind)
}
map_vecchia <- max.col(prob_vecchia, ties.method = "first")
map_exact   <- max.col(prob_exact, ties.method = "first")
tv_dist     <- 0.5 * rowSums(abs(prob_vecchia - prob_exact))
is_cpt      <- seq_len(cs$N1) %in% cs$Z2_ind
Z1          <- dl$Z1
strata_summary <- tibble(
  sites = c("All lattice sites", "Training CPT sites", "Other sites"),
  n = c(cs$N1, sum(is_cpt), sum(!is_cpt)),
  share_map_differs = c(mean(map_vecchia != map_exact),
                        mean((map_vecchia != map_exact)[is_cpt]),
                        mean((map_vecchia != map_exact)[!is_cpt])),
  mean_tv_distance = c(mean(tv_dist), mean(tv_dist[is_cpt]), mean(tv_dist[!is_cpt])),
  p_differs_from_Z1_vecchia = c(mean(1 - prob_vecchia[cbind(seq_len(cs$N1), Z1)]),
                                mean(1 - prob_vecchia[cbind(seq_len(cs$N1), Z1)][is_cpt]),
                                mean(1 - prob_vecchia[cbind(seq_len(cs$N1), Z1)][!is_cpt])),
  p_differs_from_Z1_exact = c(mean(1 - prob_exact[cbind(seq_len(cs$N1), Z1)]),
                              mean(1 - prob_exact[cbind(seq_len(cs$N1), Z1)][is_cpt]),
                              mean(1 - prob_exact[cbind(seq_len(cs$N1), Z1)][!is_cpt]))
)
if (!weight_diag$reliable) {
  strata_summary <- select(strata_summary, sites, n, p_differs_from_Z1_vecchia)
}

## 6.3 Test-set predictions ----
predict_index <- which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
pred_file <- file.path(out_path, "test_predictions.rds")
if (file.exists(pred_file)) {
  message("Loading cached test predictions from ", pred_file)
  test_pred <- readRDS(pred_file)
} else {
  set.seed(seed)
  test_pred <- produce_prediction(
    samples_post,
    geomix_setup,
    nugget = TRUE,
    include_samples = TRUE,
    run_parallel = TRUE,
    mc.cores = n_cores,
    predict_index = predict_index
  )
  saveRDS(test_pred, pred_file)
}
stopifnot(ncol(test_pred$samples) == length(w))

test_obs <- box_cox(geomix_setup$df$qc[predict_index], lambda)
test_loc <- geomix_setup$df$loc_id[predict_index]
metrics <- cbind(posterior = "Vecchia", score_predictions(test_obs, test_pred$samples, test_loc))
if (weight_diag$reliable) {
  metrics <- bind_rows(
    metrics,
    cbind(posterior = "Exact (reweighted)",
          score_predictions(test_obs, test_pred$samples[, resample_idx, drop = FALSE], test_loc))
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 Save tables and figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write.csv(loglik_df,      file.path(out_path, "tables", "loglik_draws.csv"), row.names = FALSE)
write.csv(full_check,     file.path(out_path, "tables", "exact_form_check.csv"), row.names = FALSE)
write.csv(class_gap,      file.path(out_path, "tables", "gap_by_class.csv"), row.names = FALSE)
write.csv(weight_diag,    file.path(out_path, "tables", "weight_diagnostics.csv"), row.names = FALSE)
write.csv(chain_diag,     file.path(out_path, "tables", "weight_diagnostics_by_chain.csv"), row.names = FALSE)
write.csv(param_summary,  file.path(out_path, "tables", "parameter_summary.csv"), row.names = FALSE)
write.csv(strata_summary, file.path(out_path, "tables", "strata_summary.csv"), row.names = FALSE)
write.csv(metrics,        file.path(out_path, "tables", "test_metrics.csv"), row.names = FALSE)

plot_theme <- theme_bw() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        axis.text = element_text(size = 8), axis.title = element_text(size = 10))

g_trace <- loglik_df %>%
  pivot_longer(c(vecchia, exact), names_to = "likelihood", values_to = "value") %>%
  ggplot() +
  geom_line(aes(x = iteration, y = value, col = likelihood), linewidth = 0.3) +
  facet_wrap(vars(chain), ncol = 2, labeller = label_both) +
  labs(x = "Iteration", y = "Z2 log-likelihood", col = "") + plot_theme

g_gap <- ggplot(loglik_df) +
  geom_line(aes(x = iteration, y = gap, col = factor(chain)), linewidth = 0.3, alpha = 0.7) +
  labs(x = "Iteration", y = "Exact - Vecchia log-likelihood", col = "Chain") + plot_theme

g_scatter <- ggplot(loglik_df) +
  geom_point(aes(x = vecchia, y = exact, col = factor(chain)), size = 0.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = mean(loglik_df$gap), linetype = "dashed") +
  labs(x = "Vecchia log-likelihood", y = "Exact log-likelihood") + guides(col = "none") + plot_theme

g_weights <- ggplot(loglik_df) +
  geom_col(aes(x = draw, y = weight, fill = factor(chain)), width = 1) +
  labs(x = "Draw", y = "Normalised importance weight", fill = "Chain") + plot_theme

ggsave(file.path(out_path, "figures", "loglik_trace.pdf"), g_trace, width = 7, height = 5)
ggsave(file.path(out_path, "figures", "loglik_gap.pdf"), (g_gap | g_scatter) + plot_layout(guides = "collect") &
         theme(legend.position = "bottom"), width = 8, height = 4)
ggsave(file.path(out_path, "figures", "importance_weights.pdf"), g_weights, width = 7, height = 3.5)

if (weight_diag$reliable) {
g_params <- param_summary %>%
  mutate(type = str_remove(parameter, "\\[.*\\]")) %>%
  filter(type %in% c("sigma2", "lD", "tau2", "h")) %>%
  ggplot(aes(x = parameter, y = q50, ymin = q025, ymax = q975, col = posterior)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.2) +
  facet_wrap(vars(type), scales = "free", ncol = 2) +
  labs(x = "", y = "Posterior median and 95% interval", col = "") + plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_path, "figures", "parameter_comparison.pdf"), g_params, width = 8, height = 6)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Report ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Exact-likelihood check ===\n")
cat("\nImportance-weight diagnostics:\n");  print(as.data.frame(weight_diag))
cat("\nBy chain:\n");                       print(as.data.frame(chain_diag))
cat("\nGap (exact - Vecchia) by class:\n"); print(as.data.frame(class_gap))
if (!weight_diag$reliable) {
  cat("\nNOTE: the importance weights are degenerate (see ESS and Pareto k above).\n",
      "Reweighted (exact-likelihood) summaries were therefore not produced; only Vecchia results follow.\n", sep = "")
}
cat("\nTest-set metrics:\n");               print(metrics)
cat("\nLatent strata:\n");                  print(as.data.frame(strata_summary))
