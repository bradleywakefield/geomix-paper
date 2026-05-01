#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Libraries and helpers ----
library(sf)
library(tidyverse)
box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)
colours <-c("gold","darkorange3", "forestgreen","darkblue","maroon", "red4", "dodgerblue4","gray30")

# Requires model results in memory — run after 02_fit_models.R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Load data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Processed 3D dataset and model output ----
load("data/processed/data3D.RData")
path = "results/application/"

geomix_setup <- readRDS(paste0(path,"GeoMix_1/geomix_setup.rds"))

## 2.2 Shared plot theme ----
default_theme <- theme_bw()+ theme(
  legend.position = "bottom",
  legend.key.size = unit(0.3,"cm"),
  axis.text = element_text(size=8),
  axis.title = element_text(size=10),
  legend.title = element_text(size=10),
  legend.title.position = "top",
  legend.title.align = 0.5,
  legend.spacing.x  = unit(0.01, "cm"),
  plot.margin = margin(0, 0, 0, 0),
  strip.background = element_blank(),
  strip.placement = "inside",
  strip.switch.pad.wrap = unit(0.0, "cm"),
  strip.text = element_text(size = 8, margin = margin(0,0,0,0)),
  legend.text = element_text(size=8, margin = margin(0,0,0,1))
)

## 2.3 Load predictions and attach to setup ----
GeoMix_predictions_mean <- readRDS(paste0(path,"gridpredictions/GeoMix_predictions_mean.rds"))
GeoMix_predictions_sd <- readRDS(paste0(path,"gridpredictions/GeoMix_predictions_sd.rds"))
params <- readRDS(paste0(path,"predictions/GeoMix_params.rds"))

geomix_setup$df$Z2_mean <- GeoMix_predictions_mean
geomix_setup$df$Z2_sd <- GeoMix_predictions_sd
geomix_setup$df$Y1 <- params$params$Y1
geomix_setup$df$Y1p <- sapply(1:length(params$params$Y1),function(j) params$params$Y1prob[j,params$params$Y1[j]])
geomix_setup$df <- left_join(
  geomix_setup$df,
  distinct(full_df,xid,yid,easting=x,northing=y)
)

## 2.4 Spatial grid and cross-section lines ----
cell_grid <- readRDS('data/processed/cell_grid.rds')
line_df <- readRDS("data/processed/line_df.rds")

## 2.5 Centred coordinates for plotting ----
centred_locations <- distinct(geomix_setup$df,loc_id,x,y) %>%
  right_join(
    data.frame(
      loc_id = cell_grid$loc_id,
      st_coordinates(cell_grid$centre)
    ) %>%
      mutate(
        X = (X-min(full_df$x))/100000,
        Y = (Y-min(full_df$y))/100000
      )
  )

## 2.6 Build base cross-section data frame ----
base_df <- line_df %>%
  left_join(geomix_setup$df) %>%
  left_join(centred_locations) %>%
  group_by(line) %>%
  mutate(r = sqrt(2) * (Y - min(Y))) %>%
  ungroup() %>%
  dplyr::filter(d <= 50) %>%
  mutate(Y1 = factor(
    Y1,
    levels = 1:8,
  ))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Cross-section plots ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Select example line ----
l0 = 2

## 3.2 CPT profile along line ----
g1 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_line(aes(x=d, y=(box_cox(qc)/10)+r-1, group = loc_id))+
  scale_x_reverse()+
  coord_flip()+
  default_theme+ theme(plot.margin = margin(3,0,20,0))+
  labs(y = "Distance along line (km)", x = "Depth (m)")

## 3.3 Ground model stratigraphy (Z1) ----
g2 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_tile(aes(x=r, y=d, fill = factor(Z1)))+
  guides(fill = guide_legend(nrow = 1))+
  scale_y_reverse()+
  scale_fill_manual(values = colours)+
  default_theme+theme(legend.position = c(0.5,-0.275), plot.margin = margin(3,0,20,0),
                      legend.title = element_blank())+
  labs(x = "Distance along line (km)", y = "Depth (m)", fill = "")

## 3.4 MAP posterior stratigraphy (Y1) ----
g3 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_tile(aes(x=r, y=d, fill = Y1))+
  guides(fill = guide_legend(nrow = 1))+
  scale_y_reverse()+
  scale_fill_manual(values = colours)+
  default_theme+theme(legend.position = c(0.5,-0.275), plot.margin = margin(2,0,20,0),
                      legend.title = element_blank())+
  labs(x = "Distance along line (km)", y = "Depth (m)", fill = "")

## 3.5 Posterior mean cone tip resistance (Z2) ----
g4 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_tile(aes(x=r, y=d, fill = Z2_mean))+
  scale_y_reverse()+
  scale_fill_distiller(palette = "Spectral",
                       label = ~round(inv_box_cox(.x)),
                       breaks = box_cox(c(0, 5, 15, 30, 50, 75, 100)))+
  default_theme+
  theme(legend.position = c(0.5,-0.3),
        legend.direction = "horizontal",
        plot.margin = margin(3,0,20),
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank())+
  labs(x = "Distance along line (km)", y = "Depth (m)", fill = "")

## 3.6 Posterior standard deviation ----
g5 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_tile(aes(x = r, y = d, fill = Z2_sd),show.legend = T)+
  scale_fill_distiller(palette = "BrBG")+
  scale_y_reverse()+ default_theme+
  theme(legend.position = c(0.5,-0.3),
        legend.direction = "horizontal",
        plot.margin = margin(3,0,20),
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank())+
  labs(x="Distance along line (km)",y="",
       fill="")

## 3.7 MAP stratigraphy probability ----
g6 <- base_df %>%
  filter(line == l0) %>%
  ggplot()+
  geom_tile(aes(x = r, y = d, fill = Y1p),show.legend = T)+
  scale_fill_distiller(palette = "Blues")+
  scale_y_reverse()+ default_theme+
  theme(legend.position = c(0.5,-0.3),
        legend.direction = "horizontal",
        plot.margin = margin(3,0,20),
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank())+
  labs(x="Distance along line (km)",y="",
       fill="")

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Save figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave('results/figures/model_plot_p1.png',g1, width = 3, height = 2.5)
ggsave('results/figures/model_plot_p2.png',g2, width = 3, height = 2.5)
ggsave('results/figures/model_plot_p3.png',g3, width = 3, height = 2.5)
ggsave('results/figures/model_plot_p4.png',g4, width = 3, height = 2.5)
ggsave('results/figures/model_plot_p5.png',g5, width = 3, height = 2.5)
ggsave('results/figures/model_plot_p6.png',g6, width = 3, height = 2.5)
