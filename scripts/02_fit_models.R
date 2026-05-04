#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 1.0.1 Libraries and environment ----
library(geomix)
library(tidyverse)
library(patchwork)
library(posterior)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Simulation Study -----------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Setup paths and data ---------

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
  mcmc_control = list(niter = 2750, thin = 1,
                      nbatches = 22, save_batches =T,
                      retain_draws = F)
)

## 2.2 Run MCMC chains (parallel) ----
run_sim <- confirm_run("the simulation study example...")
if(run_sim){
run_chains(geomix_setup,
 nchains = 4,
 path = path,
 LGFM = F,
 run_parallel = TRUE,
 load_previous_state = T,
 mc.cores = 50,
 seed = 16)
}
## 2.3 Post-processing samples ----
nburnin <- 250

samples <- load_mcmc_samples(path) %>%
  map( ~.x[-(1:nburnin),]) %>%
  map(~.x[seq(10,nrow(.x),10),])

params <- extract_parameters(samples)
merged_params <- combine_chains(params)

## 2.4 Diagnostics ----
diagnostics <- run_mcmc_diagnostics(params)

### 2.4.1 Trace plots ----
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
                         strip.placement = "inside",   # optional
                         strip.switch.pad.wrap = unit(0.0, "cm"),
                         strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                         legend.text = element_text(size=9))

#### 2.4.1.1 Data setup ----
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

thin <- 1
#### 2.4.1.2 Plot 1 - alpha + other params ----
g1 <- diagGP_df %>%
  filter(type %in% c("alpha0","alpha1")) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 2)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(), plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )

g2 <- diagGP_df %>%
  filter(is.na(class)) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2, show.legend = F)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 1)+
  theme_bw()+default_theme+labs(x="",y="")+theme(axis.title = element_blank())+
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )
plot1 <- (g1 | (g2 / plot_spacer())) +
  plot_layout(heights = c(1, 1),widths = c(2.3,1))

#### 2.4.1.3 Plot 2 - covariance params ----
plot2 <- diagGP_df %>%
  filter(type %in% c("sigma2","lL","lD")) %>%
  ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 3)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(),plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 2.4.1.4 Plot 3 - gamma page 1 ----
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
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 2.4.1.5 Plot 4 - gamma page 2 ----
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
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 2.4.1.5 Plot 5 - Y1 count ----

plot5 <-  diagnostics$plots$Y1 +
  guides(col = guide_legend(ncol=2))+
  theme(legend.position = c(0.85,0.075),
        legend.title.position = "top",
        plot.margin = margin(0,7,0,0))+
  scale_x_continuous(breaks = c(0,100,200),
                     labels = \(x) format(x*thin, big.mark = ",") )

#### 2.4.1.7 Save plots ----
ggsave("results/figures/syn_trace1.png", plot1, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace2.png", plot2, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace3.png", plot3, width = 4.6, height = 6.7)
ggsave("results/figures/syn_trace4.png", plot4, width = 4.6, height = 6.7)
ggsave("results/figures/syn_traceY1.png", plot5, width = 4.6, height = 3)

## 2.5 Prediction ----
pred <- produce_prediction(
  samples,
  geomix_setup,
  nugget = F,
  include_samples = T,
  mc.cores = 50
)

### 2.5.1 Save outputs ----
dir.create(file.path(path, "predictions"), showWarnings = FALSE)

syn_data <- drop_na(synthetic_list$data) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .))

syn_pred<- synthetic_list$data %>%
  filter(is.na(Z2)) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .)) %>%
  select(!Z2)

saveRDS(syn_data,paste0(path,'predictions/data.rds'))
saveRDS(syn_pred,paste0(path,'predictions/pred_df.rds'))
saveRDS(pred, file.path(path,"predictions","GeoMix_predictions.rds"))
saveRDS(merged_params, file.path(path,"predictions","GeoMix_params.rds"))

## 2.6 Run LGFM ----
if(run_sim){
run_chains(geomix_setup,
           nchains = 4,
           path = path,
           LGFM = T,
           run_parallel = TRUE,
           load_previous_state = T,
           mc.cores = 50,
           seed = 16)
}
### 2.6.1 Load samples ----
samples_LGFM <- load_mcmc_samples(path, name = "LGFM") %>%
  map( ~.x[-(1:nburnin),]) %>%
  map(~.x[seq(10,nrow(.x),10),])

### 2.6.2 Extract parameters ----
params_LGFM <- extract_parameters(samples_LGFM)
merged_params_LGFM <- combine_chains(params_LGFM)

### 2.6.3 Diagnostics ----
diagnostics_LGFM <- run_mcmc_diagnostics(params_LGFM)

### 2.6.4 Compute predictions ----
pred_LGFM <- produce_prediction(
  samples_LGFM,
  geomix_setup,
  nugget = F,
  include_samples = T,
  mc.cores = 50
)

### 2.6.5 Save outputs ----
saveRDS(pred_LGFM, file.path(path,"predictions","LGFM_predictions.rds"))
saveRDS(merged_params_LGFM, file.path(path,"predictions","LGFM_params.rds"))

## 2.7 Run competing models ----
source('scripts/utils/fit_competing_models.R')
run_comparisons(path, depth_interval = c(0,21))

## 2.8 Compute results ----
source('scripts/03_simulation_results.R')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Application (IJV Wind Farm) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Setup paths and data ----
path <- "results/application/"

dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path, "predictions"), showWarnings = FALSE)

load("data/processed/data3D.RData")
data$dID <- as.numeric(factor(data$d))

saveRDS(full_df,file.path(path,'predictions/full_df.rds'))

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
  mcmc_control = list(niter = 5000, thin = 1,
                      nbatches = 40, save_batches =T,
                      retain_draws = F)
)

beta <- estimate_beta(geomix_setup) # - Obtained 1.237988
cov_par <- estimate_MAP_covariance(geomix_setup)
saveRDS(beta,file.path(path,'beta.rds'))
saveRDS(cov_par,file.path(path,'MAP_covariance.rds'))
cov_par <- readRDS(file.path(path,'MAP_covariance.rds'))

geomix_setup <- setupGeoMixModel(
  data,
  K = K,
  dims = dims,
  beta = 1.237988,
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
run_ijv <- confirm_run("the IJV wind farm zone example...")
## 3.2 Run MCMC chains (parallel) ----
if(run_ijv){
run_chains(geomix_setup,
           nchains = 4,
           path = path,
           LGFM = F,
           run_parallel = TRUE,
           load_previous_state = F,
           mc.cores = 50,
           seed = 16)
}
## 3.3 Process samples ----
geomix_setup <- readRDS(file.path(path,"GeoMix_1/geomix_setup.rds"))

samples_post <- load_mcmc_samples(path, index = 5:32) %>%
  map(~.x[seq(10,nrow(.x),10),]) # Thin

params_post <- extract_parameters(samples_post)
merged_params <- combine_chains(params_post)

## 3.4 Diagnostics ----
diagnostics <- run_mcmc_diagnostics(
  exclude = c("sigma2_L",paste0("lL[",1:8,"]")),
  params_post, Y1index = geomix_setup$controlGibbs$Z2_ind)

### 3.4.1 Trace plots ----
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
                         strip.placement = "inside",   # optional
                         strip.switch.pad.wrap = unit(0.0, "cm"),
                         strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                         legend.text = element_text(size=9))

#### 3.4.1.1 Data setup ----
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

#### 3.4.1.2 Plot 1 - alpha + other params ----
thin <- 10
plot1 <- diagGP_df %>%
  filter(type %in% c("alpha0","alpha1")) %>%
ggplot()+
  geom_line(aes(x=.iteration, y = value, col = factor(.chain)),
            alpha = 0.6, linewidth = 0.2)+
  facet_wrap(vars(parameter), scales = "free_y", ncol = 2)+
  theme_bw()+default_theme+labs(x="Iteration", y="",col = "Chain")+
  theme(axis.title.y = element_blank(), plot.margin = margin(0,4,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==300,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

#### 3.4.1.3 Plot 2 - covariance params ----
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
                     labels = \(x) if_else(x==300,
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
                     labels = \(x) if_else(x==300,
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

#### 3.4.1.4 Plot 3 - gamma page 1 ----
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
                     labels = \(x) if_else(x==300,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

#### 3.4.1.5 Plot 4 - gamma page 2 ----
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
                     labels = \(x) if_else(x==300,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))

plot5 <-  diagnostics$plots$Y1 +
  guides(col = guide_legend(ncol=2))+
  theme(legend.position = c(0.85,0.075),
        legend.title.position = "top",
        plot.margin = margin(0,5,0,0))+
  scale_x_continuous(breaks = c(0,100,200,300),
                     labels = \(x) if_else(x==300,
                                           paste0(format(x*thin, big.mark = ","),'  '),
                                           format(x*thin, big.mark = ",") ))
#### 3.4.1.6 Save plots ----
ggsave("results/figures/ijv_trace1.png", plot1, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace2.png", plot2, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace3.png", plot3, width = 4.6, height = 6.7)
ggsave("results/figures/ijv_trace4.png", plot4, width = 4.6, height = 6.7)

ggsave("results/figures/ijv_traceY1.png", diagnostics$plots$Y1, width = 4.6, height = 3)

## 3.5 Prediction ----
test_pred <- produce_prediction(
  samples_post,
  geomix_setup,
  nugget = T,
  include_samples = T,
  mc.cores = 50,
  predict_index = which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
)

ijv_data <- drop_na(geomix_setup$df) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .))

ijv_pred <- geomix_setup$df %>%
  filter(!is.na(qc) & is.na(Z2)) %>%
  mutate(SU = factor(Z1)) %>%
  cbind(.,model.matrix(~-1+SU,data = .)) %>%
  select(!Z2)

### 3.5.1 Save outputs ----
saveRDS(ijv_data,paste0(path,'predictions/data.rds'))
saveRDS(ijv_pred,paste0(path,'predictions/pred_df.rds'))
saveRDS(test_pred, file.path(path,"predictions","GeoMix_predictions.rds"))
saveRDS(merged_params,file.path(path,'predictions/GeoMix_params.rds'))

## 3.6 Grid predictions ----
### 3.6.1 Load Grid  ----
library(sf)

points_grid <- select(readRDS("data/processed/points_grid.rds"),xid,yid,grid_id) %>%
  st_set_geometry(value=NULL) %>%
  distinct(xid,yid,.keep_all = T) %>%
  right_join(geomix_setup$lattice_coords)

grid_cells <- unique(points_grid$grid_id)

geomix_setup$df <- geomix_setup$df %>%
  left_join(points_grid) %>%
  relocate(grid_id,.after = yid)

dir.create(file.path(path, "gridpredictions"), showWarnings = FALSE)

### 3.6.2 Run Predictions ----
grid_pred <- list()
for(i in seq_along(grid_cells)){
  cat(paste0(' |========================================| \n |        Grid Prediction ',i,' of 77         |\n |========================================| '))

  pred_index <- which(geomix_setup$df$grid_id == grid_cells[i])

  grid_pred[[i]] <- produce_prediction(
    samples_post,
    geomix_setup,
    nugget = T,
    include_samples = T,
    mc.cores = 50,
    predict_index = pred_index
  )

  grid_pred[[i]]$pred_index <- pred_index

  saveRDS(grid_pred[[i]],
          paste0(path,'gridpredictions/GeoMix_predictions_',i,'.rds'))
}

### 3.6.3 Load Grid Predictions ----
grid_pred <- lapply(
  paste0(path,'gridpredictions/GeoMix_predictions_',seq_along(grid_cells),'.rds'),
  readRDS
)

full_grid <- list()
full_grid$mean <- unlist(map(grid_pred,~.x$mean))
full_grid$sd <- unlist(map(grid_pred,~.x$sd))
full_grid$samples <- do.call(rbind,map(grid_pred,~.x$samples))
full_grid$order <- unlist(map(grid_pred,~.x$pred_index))

full_grid$mean <- full_grid$mean[order(full_grid$order)]
full_grid$sd <- full_grid$sd[order(full_grid$order)]
full_grid$samples <- full_grid$samples[order(full_grid$order),]

### 3.6.4 Save Grid Predictions ----
saveRDS(full_grid, paste0(path,'gridpredictions/GeoMix_predictions.rds'))
saveRDS(full_grid$mean, paste0(path,'gridpredictions/GeoMix_predictions_mean.rds'))
saveRDS(full_grid$sd, paste0(path,'gridpredictions/GeoMix_predictions_sd.rds'))

## 3.7 Run LGFM ----
if(run_ijv){
  run_chains(geomix_setup,
             nchains = 4,
             path = path,
             LGFM = T,
             run_parallel = TRUE,
             load_previous_state = F,
             mc.cores = 50,
             seed = 16)
}

### 3.7.1 Load samples ----
samples_LGFM <- load_mcmc_samples(path, name = "LGFM", index=1:18) %>%
  map( ~.x[-(1:nburnin),]) %>%
  map(~.x[seq(10,nrow(.x),10),])

### 3.7.2 Extract parameters ----
params_LGFM <- extract_parameters(samples_LGFM)
merged_params_LGFM <- combine_chains(params_LGFM)

### 3.7.3 Diagnostics ----
diagnostics_LGFM <- run_mcmc_diagnostics(params_LGFM)

### 3.7.4 Compute predictions ----
pred_LGFM <- produce_prediction(
  samples_LGFM,
  geomix_setup,
  nugget = T,
  include_samples = T,
  mc.cores = 50,
  predict_index = which(!is.na(geomix_setup$df$qc) & is.na(geomix_setup$df$Z2))
)

### 3.7.5 Save outputs ----
saveRDS(pred_LGFM, file.path(path,"predictions","LGFM_predictions.rds"))
saveRDS(merged_params_LGFM, file.path(path,"predictions","LGFM_params.rds"))

## 3.8 Run competing models ----
source('scripts/utils/fit_competing_models.R')
run_comparisons(path, depth_interval = c(24,50))

## 3.9 Map figures (saves data/processed/line_df.rds, data/processed/cell_grid_map.rds) ----
source('scripts/05_map_figures.R')

## 3.10 Compute results ----
source('scripts/04_application_results.R')
