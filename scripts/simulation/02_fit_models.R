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

path <- "results/simulation/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

synthetic_list <- readRDS("data/simulation/synthetic_data.rds")

geomix_setup <- setupGeoMixModel(
  data = synthetic_list$data,
  K = synthetic_list$K,
  dims = synthetic_list$dims,
  beta = synthetic_list$beta,
  m = 160,
  kappa = 0.9,
  aformula = ~ 1 + d,
  variables = list(loc = "locID", groups = "group", xID = "x", yID = "y", dID = "d"),
  penalty = synthetic_list$penalty,
  mcmc_control = list(niter = 4000, thin = 1,
                      nbatches = 32, save_batches =T,
                      retain_draws = F)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 GeoMix MCMC ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Run chains (parallel) ----
sim_mode <- ask_run_mode("Simulation study — GeoMix chains")
if (sim_mode %in% c("s", "c")) {
  run_chains(geomix_setup,
             nchains = 4,
             path = path,
             LGFM = F,
             run_parallel = TRUE,
             load_previous_state = sim_mode == "c",
             mc.cores = 4,
             seed = 16)
}

## 3.2 Post-processing samples ----
samples <- load_mcmc_samples(path,index = 5:32, thin = 10)

params <- extract_parameters(samples)
merged_params <- combine_chains(params)

## 3.3 Diagnostics ----
diagnostics <- run_mcmc_diagnostics(params)
save_diag_tables(diagnostics, path, "GeoMix")

### 3.3.1 Trace plots ----
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

#### 3.3.1.1 Data setup ----
GPc_params <- c("alpha0","alpha1","sigma2","lL","lD")
GPnc_params <- c("tau2","sigma2_L","sigma2_D","h")
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
         parameter = factor(parameter,levels = c(paste0(rep(GPc_params,8),"[",rep(1:8,each = 5),"]"),GPnc_params)))

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

#### 3.3.1.2 Plot 1 - alpha + other params ----
g1 <- diagGP_df %>%
  filter(type %in% c("alpha0","alpha1")) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 2)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(), plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )

g2 <- diagGP_df %>%
  filter(is.na(class)) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2, show.legend = F)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 1)+
  theme_bw()+default_theme+labs(x="",y="")+theme(axis.title = element_blank())+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )
plot1 <- (g1 | (g2 / plot_spacer())) +
  plot_layout(heights = c(1, 1),widths = c(2.3,1))

#### 3.3.1.3 Plot 2 - covariance params ----
plot2 <- diagGP_df %>%
  filter(type %in% c("sigma2","lL","lD")) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 3)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(),plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 3.3.1.4 Plot 3 - gamma page 1 ----
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
        plot.margin = margin(0,7,40,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 3.3.1.5 Plot 4 - gamma page 2 ----
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
        plot.margin = margin(0,7,40,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 3.3.1.6 Plot 5 - Y1 count ----
plot5 <-  diagnostics$plots$Y1 +
  guides(col = guide_legend(ncol=2))+
  theme(legend.position = c(0.85,0.075),
        legend.title.position = "top",
        plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 3.3.1.7 Save plots ----
ggsave("results/figures/syn_trace1.png", plot1, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace2.png", plot2, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace3.png", plot3, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace4.png", plot4, width = 4.6, height = 6.7)
ggsave("results/figures/syn_traceY1.png", plot5, width = 4.6, height = 3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 GeoMix predictions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir.create(file.path(path, "predictions"), showWarnings = FALSE)
sim_pred_mode <- ask_pred_mode("Simulation study — GeoMix predictions")
if (sim_pred_mode == "l") {
  pred <- readRDS(file.path(path, "predictions", "GeoMix_predictions.rds"))
} else {
  pred <- produce_prediction(
    samples,
    geomix_setup,
    nugget = F,
    include_samples = T,
    mc.cores =  parallel::detectCores()
  )
  saveRDS(pred, file.path(path, "predictions", "GeoMix_predictions.rds"))
}

### 4.1 Save outputs ----
syn_data <- drop_na(synthetic_list$data) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .))

syn_pred<- synthetic_list$data %>%
  filter(is.na(Z2)) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .)) %>%
  select(!Z2)

if (sim_pred_mode == "r"){
  saveRDS(syn_data, paste0(path, 'predictions/data.rds'))
  saveRDS(syn_pred, paste0(path, 'predictions/pred_df.rds'))
  saveRDS(merged_params, file.path(path, "predictions", "GeoMix_params.rds"))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 LGFM ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 5.1 Run chains ----
sim_lgfm_mode <- ask_run_mode("Simulation study — LGFM chains")
if (sim_lgfm_mode %in% c("s", "c")) {
  run_chains(geomix_setup,
             nchains = 4,
             path = path,
             LGFM = T,
             run_parallel = TRUE,
             load_previous_state = sim_lgfm_mode == "c",
             mc.cores = 4,
             seed = 16)
}

## 5.2 Load and extract ----
samples_LGFM <- load_mcmc_samples(path, name = "LGFM", index = 5:32, thin = 10)
params_LGFM <- extract_parameters(samples_LGFM)
merged_params_LGFM <- combine_chains(params_LGFM)

## 5.3 Diagnostics ----
diagnostics_LGFM <- run_mcmc_diagnostics(params_LGFM, name = "LGFM")
save_diag_tables(diagnostics_LGFM, path, "LGFM")

## 5.4 Predictions ----
sim_lgfm_pred_mode <- ask_pred_mode("Simulation study — LGFM predictions")
if (sim_lgfm_pred_mode == "l") {
  pred_LGFM <- readRDS(file.path(path, "predictions", "LGFM_predictions.rds"))
} else {
  pred_LGFM <- produce_prediction(
    samples_LGFM,
    geomix_setup,
    nugget = F,
    include_samples = T,
    mc.cores =  parallel::detectCores()
  )
  saveRDS(pred_LGFM, file.path(path, "predictions", "LGFM_predictions.rds"))
}

if (sim_lgfm_pred_mode == "r") saveRDS(merged_params_LGFM, file.path(path, "predictions", "LGFM_params.rds"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Competing models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source('scripts/utils/fit_competing_models.R')
sim_comp_mode <- ask_pred_mode("Simulation study — competing models")
if (sim_comp_mode == "r") {
  run_comparisons(path, depth_interval = c(0,21))
}


