#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Profile comparison of the exact and grouped Vecchia Z2 likelihoods (IJV)
#
# Uses the existing GeoMix chains in results/application/ (no new MCMC).
# The latent strata are fixed at their most probable class and all other
# parameters at their posterior means. The exact and the grouped Vecchia Z2
# log-likelihoods are then compared as functions of the GP parameters:
#   (i)   one-dimensional curves in each class's vertical length-scale lD[k]
#         and process variance sigma2[k], and in the nugget tau2;
#   (ii)  the maximiser and curvature of each curve;
#   (iii) the joint maximiser over (sigma2[k], lD[k]) for each class.
#
# Usage (from the project root):
#   Rscript scripts/sensitivity/vecchia_profile_check.R
#
# Reads : results/application/GeoMix_*/   (read-only)
# Writes: results/sensitivity/vecchia_profile/
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 1.0.1 Libraries ----
library(geomix)
library(tidyverse)
library(patchwork)

### 1.0.2 Settings ----
base_path   <- "results/application"                     # existing chains (read-only)
out_path    <- "results/sensitivity/vecchia_profile"
batch_index <- 5:32      # batches retained in the paper (500 burn-in iterations discarded)
thin        <- 10        # as in the paper
n_chains    <- 4
n_cores     <- max(1, parallel::detectCores() - 2)
grid_range  <- c(0.5, 2) # curves span this multiple of the posterior mean
n_grid      <- 41

dir.create(file.path(out_path, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_path, "tables"), showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Load setup and posterior ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geomix_setup <- readRDS(file.path(base_path, "GeoMix_1/geomix_setup.rds"))
dl <- geomix_setup$data_list
cs <- geomix_setup$constants
K  <- cs$K

st <- list(
  K = K, Z2 = dl$Z2, X = dl$X, Z2_ind = cs$Z2_ind,
  dID = dl$dID, locID = dl$locID,
  dIDv = dl$dID[cs$Z2_ind], locIDv = dl$locID[cs$Z2_ind],
  distD = dl$distD, distL = dl$distL, LFlag = dl$LFlag
)

## 2.1 Retained draws: parameters and the latent strata at the CPT sites ----
# Same rows as load_mcmc_samples(path, index = batch_index, thin = thin), but
# only the columns needed here are kept.
message("Loading thinned samples...")
Y1_names <- paste0("Y1[", cs$Z2_ind, "]")
draws <- do.call(rbind, lapply(seq_len(n_chains), function(ch) {
  offset <- 0
  out <- vector("list", length(batch_index))
  for (i in seq_along(batch_index)) {
    batch <- readRDS(file.path(base_path, paste0("GeoMix_", ch), paste0("batch_", batch_index[i], ".rds")))
    keep <- which((offset + seq_len(nrow(batch))) %% thin == 0)
    cols <- c(grep("^Y1\\[", colnames(batch), invert = TRUE, value = TRUE), Y1_names)
    out[[i]] <- batch[keep, cols, drop = FALSE]
    offset <- offset + nrow(batch)
  }
  do.call(rbind, out)
}))

## 2.2 Reference point: most probable strata and posterior means ----
Y1_map_cpt <- apply(draws[, Y1_names], 2, function(x) which.max(tabulate(x, nbins = K)))
Y1_ref <- rep(1, cs$N1)          # only the entries at Z2_ind enter the Z2 likelihood
Y1_ref[cs$Z2_ind] <- Y1_map_cpt

post <- extract_parameters(draws[, !colnames(draws) %in% Y1_names])
ref <- list(
  alpha  = apply(post$samples$alpha, c(2, 3), mean),
  sigma2 = unname(post$params$sigma2),
  lD     = unname(post$params$lD),
  lL     = unname(post$params$lL),
  tau2   = post$params$tau2
)
post_sd <- list(sigma2 = apply(post$samples$sigma2, 2, sd), lD = apply(post$samples$lD, 2, sd),
                tau2 = sd(post$samples$tau2))
n_class <- tabulate(Y1_map_cpt, nbins = K)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Likelihood functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Exact log-likelihood of class k (dense Matern-3/2 covariance plus nugget) ----
exact_loglik_k <- function(k, par) {
  idx <- which(Y1_map_cpt == k)
  n <- length(idx)
  if (n == 0) return(0)
  r2 <- st$distD[st$dIDv[idx], st$dIDv[idx], drop = FALSE] / par$lD[k]^2
  if (st$LFlag[k] != 0) {
    r2 <- r2 + st$distL[st$locIDv[idx], st$locIDv[idx], drop = FALSE] / par$lL[k]^2
  }
  r <- sqrt(3 * r2)
  C <- par$sigma2[k] * (1 + r) * exp(-r)
  diag(C) <- diag(C) + par$tau2
  U <- chol(C)
  res <- st$Z2[idx] - drop(st$X[idx, , drop = FALSE] %*% par$alpha[k, ])
  v <- backsolve(U, res, transpose = TRUE)
  -sum(log(diag(U))) - 0.5 * sum(v^2) - 0.5 * n * log(2 * pi)
}

## 3.2 Grouped Vecchia log-likelihood by class (package function, as used in the MCMC) ----
CevaluateGroupGPvec <- nimble::compileNimble(evaluateGroupGPvec)
vecchia_loglik <- function(par) {
  colSums(CevaluateGroupGPvec(
    Z2 = st$Z2, alpha = par$alpha, sigma2 = par$sigma2, tau2 = par$tau2, lL = par$lL, lD = par$lD,
    LFlag = st$LFlag, Y1 = Y1_ref, X = st$X, K = K, Z2_ind = st$Z2_ind,
    dID = st$dID, locID = st$locID, distD = st$distD, distL = st$distL,
    m = cs$m, groupLookup = dl$groupLookup, groupNum = dl$groupNum,
    groupNeighbours = dl$groupNeighbours
  ))
}

## 3.3 Log-likelihood as a function of one parameter ----
# For lD[k] and sigma2[k] only class k changes, so the class-k term is returned;
# for tau2 (shared by all classes) the total is returned.
set_par <- function(name, k, value) {
  par <- ref
  if (name == "tau2") par$tau2 <- value else par[[name]][k] <- value
  par
}
loglik_at <- function(likelihood, name, k, value) {
  par <- set_par(name, k, value)
  if (likelihood == "Vecchia") {
    ll <- vecchia_loglik(par)
    if (name == "tau2") sum(ll) else ll[k]
  } else {
    if (name == "tau2") sum(sapply(seq_len(K), exact_loglik_k, par = par)) else exact_loglik_k(k, par)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 One-dimensional curves ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

targets <- bind_rows(
  tibble(name = "lD", k = seq_len(K), ref_value = ref$lD, post_sd = post_sd$lD),
  tibble(name = "sigma2", k = seq_len(K), ref_value = ref$sigma2, post_sd = post_sd$sigma2),
  tibble(name = "tau2", k = NA_integer_, ref_value = ref$tau2, post_sd = post_sd$tau2)
) %>%
  mutate(parameter = if_else(name == "tau2", "tau2", paste0(name, "[", k, "]")))

curve_file <- file.path(out_path, "curves.rds")
if (file.exists(curve_file)) {
  message("Loading cached curves from ", curve_file)
  curves <- readRDS(curve_file)
} else {
  tasks <- targets %>%
    crossing(likelihood = c("Vecchia", "Exact"),
             multiple = exp(seq(log(grid_range[1]), log(grid_range[2]), length.out = n_grid))) %>%
    mutate(value = ref_value * multiple)
  message("Evaluating ", nrow(tasks), " curve points on ", n_cores, " cores...")
  tasks$loglik <- unlist(parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    loglik_at(tasks$likelihood[i], tasks$name[i], tasks$k[i], tasks$value[i])
  }, mc.cores = n_cores, mc.preschedule = FALSE))
  curves <- tasks
  saveRDS(curves, curve_file)
}
curves <- curves %>%
  group_by(parameter, likelihood) %>%
  mutate(rel_loglik = loglik - max(loglik)) %>%
  ungroup()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Maximisers and curvature ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Each curve is maximised on the log scale; the curvature there gives the
# standard deviation of the corresponding normal approximation.
max_file <- file.path(out_path, "maximisers.rds")
if (file.exists(max_file)) {
  message("Loading cached maximisers from ", max_file)
  maximisers <- readRDS(max_file)
} else {
  tasks <- crossing(targets, likelihood = c("Vecchia", "Exact"))
  message("Maximising ", nrow(tasks), " curves...")
  res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    f <- function(lv) loglik_at(tasks$likelihood[i], tasks$name[i], tasks$k[i], exp(lv))
    opt <- optimize(f, log(tasks$ref_value[i] * grid_range), maximum = TRUE, tol = 1e-4)
    step <- 0.02
    d2 <- (f(opt$maximum + step) - 2 * opt$objective + f(opt$maximum - step)) / step^2
    c(maximiser = exp(opt$maximum), sd_log = if (d2 < 0) 1 / sqrt(-d2) else NA_real_)
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  maximisers <- bind_cols(tasks, as_tibble(do.call(rbind, res)))
  saveRDS(maximisers, max_file)
}

max_table <- maximisers %>%
  mutate(sd = maximiser * sd_log) %>%
  select(parameter, name, k, posterior_mean = ref_value, posterior_sd = post_sd, likelihood, maximiser, sd) %>%
  pivot_wider(names_from = likelihood, values_from = c(maximiser, sd)) %>%
  mutate(
    n_obs = if_else(is.na(k), sum(n_class), n_class[k]),
    ratio_exact_to_vecchia = maximiser_Exact / maximiser_Vecchia,
    shift_in_posterior_sd  = (maximiser_Exact - maximiser_Vecchia) / posterior_sd,
    shift_in_curve_sd      = (maximiser_Exact - maximiser_Vecchia) / sd_Vecchia,
    curve_sd_ratio         = sd_Exact / sd_Vecchia,
    at_grid_edge = pmin(abs(log(maximiser_Exact / posterior_mean / grid_range[1])),
                        abs(log(maximiser_Exact / posterior_mean / grid_range[2])),
                        abs(log(maximiser_Vecchia / posterior_mean / grid_range[1])),
                        abs(log(maximiser_Vecchia / posterior_mean / grid_range[2]))) < 1e-3
  ) %>%
  arrange(name, k)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Joint maximiser over (sigma2[k], lD[k]) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

joint_file <- file.path(out_path, "joint_maximisers.rds")
if (file.exists(joint_file)) {
  message("Loading cached joint maximisers from ", joint_file)
  joint <- readRDS(joint_file)
} else {
  tasks <- crossing(k = seq_len(K), likelihood = c("Vecchia", "Exact"))
  message("Jointly maximising over (sigma2, lD) for each class...")
  res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    k <- tasks$k[i]
    f <- function(lv) {
      par <- ref
      par$sigma2[k] <- exp(lv[1]); par$lD[k] <- exp(lv[2])
      if (tasks$likelihood[i] == "Vecchia") -vecchia_loglik(par)[k] else -exact_loglik_k(k, par)
    }
    opt <- optim(log(c(ref$sigma2[k], ref$lD[k])), f, method = "Nelder-Mead",
                 control = list(reltol = 1e-9, maxit = 500))
    c(sigma2 = exp(opt$par[1]), lD = exp(opt$par[2]), loglik = -opt$value, convergence = opt$convergence)
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  joint <- bind_cols(tasks, as_tibble(do.call(rbind, res)))
  saveRDS(joint, joint_file)
}

joint_table <- joint %>%
  select(-loglik) %>%
  pivot_wider(names_from = likelihood, values_from = c(sigma2, lD, convergence)) %>%
  mutate(
    n_obs = n_class[k],
    posterior_mean_sigma2 = ref$sigma2[k], posterior_sd_sigma2 = post_sd$sigma2[k],
    posterior_mean_lD = ref$lD[k], posterior_sd_lD = post_sd$lD[k],
    sigma2_shift_in_posterior_sd = (sigma2_Exact - sigma2_Vecchia) / posterior_sd_sigma2,
    lD_shift_in_posterior_sd = (lD_Exact - lD_Vecchia) / posterior_sd_lD
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 Save tables and figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write.csv(curves,      file.path(out_path, "tables", "curves.csv"), row.names = FALSE)
write.csv(max_table,   file.path(out_path, "tables", "maximisers.csv"), row.names = FALSE)
write.csv(joint_table, file.path(out_path, "tables", "joint_maximisers.csv"), row.names = FALSE)

plot_theme <- theme_bw() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        axis.text = element_text(size = 8), axis.title = element_text(size = 10))

# Curves are shown relative to their own maximum and cut at 10 log units below it;
# the shaded band is the MCMC posterior mean +/- 2 posterior sd.
curve_plot <- function(par_name, xlab) {
  band <- filter(targets, name == par_name) %>%
    mutate(lo = ref_value - 2 * post_sd, hi = ref_value + 2 * post_sd)
  curves %>%
    filter(name == par_name, rel_loglik > -10) %>%
    ggplot() +
    geom_rect(data = band, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "grey85") +
    geom_line(aes(x = value, y = rel_loglik, col = likelihood)) +
    geom_point(aes(x = value, y = rel_loglik, col = likelihood), size = 0.6) +
    facet_wrap(vars(parameter), scales = "free_x", ncol = 4) +
    labs(x = xlab, y = "Log-likelihood relative to maximum", col = "") + plot_theme
}
ggsave(file.path(out_path, "figures", "profile_lD.pdf"), curve_plot("lD", "Vertical length-scale"),
       width = 9, height = 5)
ggsave(file.path(out_path, "figures", "profile_sigma2.pdf"), curve_plot("sigma2", "Process variance"),
       width = 9, height = 5)
ggsave(file.path(out_path, "figures", "profile_tau2.pdf"), curve_plot("tau2", "Nugget variance"),
       width = 4, height = 3.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Report ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("\n=== Exact versus Vecchia likelihood: profile comparison ===\n")
cat("\nOne-dimensional maximisers (other parameters at posterior means, strata at MAP):\n")
print(as.data.frame(select(max_table, parameter, n_obs, posterior_mean, posterior_sd, maximiser_Vecchia,
                           maximiser_Exact, ratio_exact_to_vecchia, shift_in_posterior_sd,
                           curve_sd_ratio, at_grid_edge)), digits = 3)
cat("\nJoint maximisers over (sigma2[k], lD[k]):\n")
print(as.data.frame(select(joint_table, k, n_obs, sigma2_Vecchia, sigma2_Exact, sigma2_shift_in_posterior_sd,
                           lD_Vecchia, lD_Exact, lD_shift_in_posterior_sd,
                           convergence_Vecchia, convergence_Exact)), digits = 3)
