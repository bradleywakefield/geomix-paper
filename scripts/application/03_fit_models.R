#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 1.0.1 Libraries and environment ----
library(geomix)
library(tidyverse)
library(patchwork)
library(posterior)

### 1.0.2 Helpers ----
source('scripts/utils/ask_prompts.R')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Setup paths and data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path <- "results/application/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path, "predictions"), showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

load("data/processed/data3D.RData")

geomix_setup <- setupGeoMixModel(
  data,
  K = K,
  dims = dims,
  beta = 1,
  m = 160,
  kappa = 0.9,
  aformula = ~1+d,
  variables = list(
    loc = "loc_id",
    xID = "xid",
    yID = "yid",
    dID = "dID",
    groups = "groups"
  ),
  mcmc_control = list(niter = 4000, thin = 1,
                      nbatches = 32, save_batches =T,
                      retain_draws = F)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Estimate spatial parameters ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ijv_beta_mode <- ask_pred_mode("IJV wind farm — estimate beta")
if (ijv_beta_mode == "l") {
  beta <- readRDS(file.path(path, 'beta.rds'))
} else {
  beta <- estimate_beta(geomix_setup)$estimate # - Obtained 1.237988
  saveRDS(beta, file.path(path, 'beta.rds'))
}

ijv_cov_mode <- ask_pred_mode("IJV wind farm — estimate MAP covariance")
if (ijv_cov_mode == "l") {
  cov_par <- readRDS(file.path(path, 'MAP_covariance.rds'))
} else {
  cov_par <- estimate_MAP_covariance(geomix_setup)
  saveRDS(cov_par, file.path(path, 'MAP_covariance.rds'))
}

geomix_setup <- setupGeoMixModel(
  data,
  K = K,
  dims = dims,
  beta = beta,
  m = 160,
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 GeoMix MCMC ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Run chains (parallel) ----
ijv_mode <- ask_run_mode("IJV wind farm — GeoMix chains")
if (ijv_mode %in% c("s", "c")) {
  run_chains(geomix_setup,
             nchains = 4,
             path = path,
             LGFM = F,
             run_parallel = TRUE,
             load_previous_state = ijv_mode == "c",
             mc.cores = 4,
             seed = 16)
}

## 4.2 Process samples ----
geomix_setup <- readRDS(file.path(path,"GeoMix_1/geomix_setup.rds"))

samples_post <- load_mcmc_samples(path, index = 5:32, thin = 10)

params_post <- extract_parameters(samples_post)
merged_params <- combine_chains(params_post)

## 4.3 Diagnostics ----
diagnostics <- run_mcmc_diagnostics(
  exclude = c("sigma2_L",paste0("lL[",1:8,"]")),
  params_post, Y1index = geomix_setup$controlGibbs$Z2_ind)
save_diag_tables(diagnostics, path, "GeoMix")

### 4.3.1 Trace plots ----
default_theme <-   theme(legend.position = "bottom",
                         legend.key.size = unit(0.5,"cm"),
                         axis.text = element_text(size=8),
                         axis.title = element_text(size=10),
                         legend.title = element_text(size=10),
                         plot.subtitle = element_text(size=8),
                         legend.title.align = 0.5,
                         legend.spacing.x  = unit(0.02, "cm"),
                         legend.margin = margin(0,0,0,0),
                         plot.margin = margin(0, 0, 0, 0),
                         strip.background = element_blank(),
                         strip.placement = "inside",
                         strip.switch.pad.wrap = unit(0.0, "cm"),
                         strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                         legend.text = element_text(size=9))

#### 4.3.1.1 Data setup ----
GPc_params <- c("alpha0","alpha1","sigma2","lD")
GPnc_params <- c("tau2","sigma2_D","h")
GP_params <- c(GPc_params,GPnc_params)

diag_df <- as_draws_df(diagnostics$draws$draws_array) %>%
  group_by(.chain) %>%
  mutate(.iteration = row_number()) %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "parameter",
    values_to = "value"
  )

diagGP_df <- diag_df %>%
  filter(str_detect(parameter,paste0(GP_params,collapse = "|"))) %>%
  mutate(class = as.integer(str_match(parameter, "\\[(\\d+)\\]")[,2]),
         type = str_remove(parameter,"\\[\\d+\\]"),
         type = factor(type,levels = GP_params),
         parameter = factor(parameter,levels = c(paste0(rep(GPc_params,8),"[",rep(1:8,each = 4),"]"),GPnc_params)))

diagCM_df <- diag_df %>%
  filter(str_detect(parameter,"gamma")) %>%
  mutate(
    idx = str_match(parameter, "\\[(\\d+), (\\d+)\\]"),
    k = as.integer(idx[,2]),
    l = as.integer(idx[,3]),
    parameter = str_remove(parameter,"Mat")
  ) %>%
  select(-idx)

thin <- 10

#### 4.3.1.2 Plot 1 - alpha + other params ----
plot1 <- diagGP_df %>%
  filter(type %in% c("alpha0","alpha1")) %>%
ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 2)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(), plot.margin = margin(0,4,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

#### 4.3.1.3 Plot 2 - covariance params ----
g1 <- diagGP_df %>%
  filter(type %in% c("sigma2","lD")) %>%
  mutate(parameter = factor(parameter,levels = paste0(rep(c("sigma2","lD"),8),"[",rep(1:8,each = 2),"]"))) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 2, drop = FALSE)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(),plot.margin = margin(0,5,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))
g2 <- diagGP_df %>%
  filter(type %in% GPnc_params) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2, show.legend = F)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 1)+
  theme_bw()+default_theme+labs(x="",y="")+theme(axis.title = element_blank())+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'   '),
                                           format(x*thin, big.mark = ",") ))
right_col <- (g2 / plot_spacer()) +
  plot_layout(heights = c(3, 5))

plot2 <- (g1 | right_col) +
  plot_layout(
    widths = c(2, 1),
    guides = "collect"
  ) &
  theme(legend.position = "bottom")

#### 4.3.1.4 Plot 3 - gamma page 1 ----
plot3 <- diagCM_df %>%
  filter(k <= 4) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 3)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(legend.position = c(0.7,0.01), legend.direction = "horizontal",
        axis.title.x = element_text(margin = margin(t = -45)),
        axis.title.y = element_blank(),
        plot.margin = margin(0,5,40,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

#### 4.3.1.5 Plot 4 - gamma page 2 ----
plot4 <- diagCM_df %>%
  filter(k > 4) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 3)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(legend.position = c(0.7,0.01), legend.direction = "horizontal",
        axis.title.x = element_text(margin = margin(t = -45)),
        axis.title.y = element_blank(),
        plot.margin = margin(0,5,40,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

plot5 <-  diagnostics$plots$Y1 +
  guides(col = guide_legend(ncol=2))+
  theme(legend.position = c(0.85,0.075),
        legend.title.position = "top",
        plot.margin = margin(0,5,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==350,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

#### 4.3.1.6 Save plots ----
ggsave("results/figures/ijv_trace1.png", plot1, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace2.png", plot2, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace3.png", plot3, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace4.png", plot4, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_traceY1.png", diagnostics$plots$Y1, width = 4.6, height = 3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 GeoMix predictions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ijv_pred_mode <- ask_pred_mode("IJV wind farm — GeoMix test predictions")
if (ijv_pred_mode == "l") {
  test_pred <- readRDS(file.path(path, "predictions", "GeoMix_predictions.rds"))
} else {
  test_pred <- produce_prediction(
    samples_post,
    geomix_setup,
    nugget = T,
    include_samples = T,
    mc.cores =  parallel::detectCores(),
    predict_index = which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
  )
  saveRDS(test_pred, file.path(path, "predictions", "GeoMix_predictions.rds"))
}

ijv_data <- drop_na(geomix_setup$df) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .))

ijv_pred <- geomix_setup$df %>%
  filter(!is.na(qc) & is.na(Z2)) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .)) %>%
  select(!Z2)

### 5.1 Save outputs ----
if (ijv_pred_mode == "r") {
  saveRDS(full_df, file.path(path, 'predictions/full_df.rds'))
  saveRDS(ijv_data, paste0(path, 'predictions/data.rds'))
  saveRDS(ijv_pred, paste0(path, 'predictions/pred_df.rds'))
  saveRDS(merged_params, file.path(path, 'predictions/GeoMix_params.rds'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Grid predictions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 Load grid ----
library(sf)

points_grid <- select(readRDS("data/processed/points_grid.rds"),xid,yid,grid_id) %>%
  st_set_geometry(value=NULL) %>%
  distinct(xid,yid,.keep_all = T) %>%
  right_join(geomix_setup$lattice_coords)

grid_cells <- unique(points_grid$grid_id)

geomix_setup$df <- geomix_setup$df %>%
  left_join(points_grid) %>%
  relocate(grid_id,.after = yid)

## 6.2 Run predictions ----
ijv_grid_mode <- ask_pred_mode("IJV wind farm — GeoMix grid predictions")
if (ijv_grid_mode == "l") {
  full_grid <- readRDS(paste0(path, 'gridpredictions/GeoMix_predictions.rds'))
} else {
  dir.create(file.path(path, "gridpredictions"), showWarnings = FALSE)
  grid_pred <- list()
  for(i in seq_along(grid_cells)){
    cat(paste0(' |========================================| \n |        Grid Prediction ',i,' of 77         |\n |========================================| '))

    pred_index <- which(geomix_setup$df$grid_id == grid_cells[i])

    grid_pred[[i]] <- produce_prediction(
      samples_post,
      geomix_setup,
      nugget = T,
      include_samples = T,
      mc.cores =  parallel::detectCores(),
      predict_index = pred_index
    )

    grid_pred[[i]]$pred_index <- pred_index

    saveRDS(grid_pred[[i]],
            paste0(path,'gridpredictions/GeoMix_predictions_',i,'.rds'))
  }

  ## 6.3 Combine and save ----
  full_grid <- list()
  full_grid$mean    <- unlist(map(grid_pred, ~.x$mean))
  full_grid$sd      <- unlist(map(grid_pred, ~.x$sd))
  full_grid$samples <- do.call(rbind, map(grid_pred, ~.x$samples))
  full_grid$order   <- unlist(map(grid_pred, ~.x$pred_index))

  full_grid$mean    <- full_grid$mean[order(full_grid$order)]
  full_grid$sd      <- full_grid$sd[order(full_grid$order)]
  full_grid$samples <- full_grid$samples[order(full_grid$order), ]

  saveRDS(full_grid, paste0(path, 'gridpredictions/GeoMix_predictions.rds'))
  saveRDS(full_grid$mean, paste0(path, 'gridpredictions/GeoMix_predictions_mean.rds'))
  saveRDS(full_grid$sd,   paste0(path, 'gridpredictions/GeoMix_predictions_sd.rds'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 LGFM ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 7.1 Run chains ----
ijv_lgfm_mode <- ask_run_mode("IJV wind farm — LGFM chains")
if (ijv_lgfm_mode %in% c("s", "c")) {
  run_chains(geomix_setup,
             nchains = 4,
             path = path,
             LGFM = T,
             run_parallel = TRUE,
             load_previous_state = ijv_lgfm_mode == "c",
             mc.cores = parallel::detectCores(),
             seed = 16)
}

## 7.2 Load and extract ----
samples_LGFM <- load_mcmc_samples(path, name = "LGFM", index=5:32, thin = 10)
params_LGFM <- extract_parameters(samples_LGFM)
merged_params_LGFM <- combine_chains(params_LGFM)

## 7.3 Diagnostics ----
diagnostics_LGFM <- run_mcmc_diagnostics(params_LGFM, name = "LGFM",
                                         exclude = c("sigma2_L",paste0("lL[",1:8,"]")),
                                         Y1index = geomix_setup$controlGibbs$Z2_ind)
save_diag_tables(diagnostics_LGFM, path, "LGFM")

## 7.4 Predictions ----
ijv_lgfm_pred_mode <- ask_pred_mode("IJV wind farm — LGFM predictions")
if (ijv_lgfm_pred_mode == "l") {
  pred_LGFM <- readRDS(file.path(path, "predictions", "LGFM_predictions.rds"))
} else {
  pred_LGFM <- produce_prediction(
    samples_LGFM,
    geomix_setup,
    nugget = T,
    include_samples = T,
    mc.cores = parallel::detectCores(),
    predict_index = which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
  )
  saveRDS(pred_LGFM, file.path(path, "predictions", "LGFM_predictions.rds"))
}

if (ijv_lgfm_pred_mode == "r") saveRDS(merged_params_LGFM, file.path(path, "predictions", "LGFM_params.rds"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Competing models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source('scripts/utils/fit_competing_models.R')
ijv_comp_mode <- ask_pred_mode("IJV wind farm — competing models")
if (ijv_comp_mode == "r") {
  run_comparisons(path, depth_interval = c(24,50))
}


