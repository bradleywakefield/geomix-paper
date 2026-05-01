#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(sf)
library(patchwork)
library(tidyverse)

source('scripts/utils/prediction_scoring.R')
load("data/processed/data3D.RData")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   2. Data setup: Load data, run evaluation functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Global parameters ----

### 2.1.1 Randomly select CPTs for visualisation ----
set.seed(16)
select_locID <- test_locs[sample(1:length(test_locs),16,replace = F)]

### 2.1.2 Depth resolution ----
Ddelta <- max(diff(unique(data$d)))

### 2.1.3 Load GeoMix setup object ----
geomix_setup <- readRDS(paste0(path,"GeoMix_1/geomix_setup.rds"))

### 2.1.4 Extract training locations ----
train_locs <- distinct(drop_na(geomix_setup$df),x,y)

## 2.2 Evaluate test performance ----
### 2.2.1 Model names to evaluate ----
names <- "GeoMix"
### 2.2.2  Evaluate predictions for each model ----
prediction_list <- lapply(names,eval_prediction_Z2,path = path)

### 2.2.3  Combine outputs into tidy data frames ----
stats_df <- bind_rows(map(prediction_list,~.x$stats))
pred_df <- bind_rows(map(prediction_list,~.x$pred_df))

## 2.4 Evaluate classification performance ----
prob_list <- list(
  eval_probability_Z1("GeoMix",path),
  eval_probability_Z1("LGFM",path)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Test Plots ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 3.0 Shared ggplot theme across all figures ----
default_theme <-   theme(legend.position = "bottom",
                         legend.key.size = unit(0.4,"cm"),
                         axis.text = element_text(size=8),
                         axis.title = element_text(size=10),
                         legend.title = element_text(size=10),
                         plot.subtitle = element_text(size=8),
                         legend.title.align = 0.5,
                         legend.spacing.x  = unit(0.01, "cm"),
                         legend.margin = margin(0,0,0,0),
                         plot.margin = margin(0, 0, 0, 0),
                         strip.background = element_blank(),
                         strip.placement = "inside",   # optional
                         strip.switch.pad.wrap = unit(0.0, "cm"),
                         strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                         legend.text = element_text(size=9))

## 3.1 Test CPT profiles qc ----
gtest <- pred_df %>%
  filter(loc_id %in% select_locID) %>%
  filter(model %in% c("GeoMix", "LM")) %>%
  group_by(loc_id) %>%
  mutate(
    dist = min(sqrt((x - train_locs$x)^2 + (y - train_locs$y)^2)),
    label = paste0(
      "CPT",
      str_remove(str_remove(name, "-TNO"), "IJV"),
      " - ",
      sprintf("%1.1f", dist),
      " km*"
    )
  ) %>%
  ungroup() %>%
  ggplot(aes(x = d)) +
  geom_line(
    aes(y = Z2, colour = "actual"),
    linewidth = 0.4
  ) +
  geom_ribbon(
    data = \(df) dplyr::filter(df, model == "GeoMix"),
    aes(
      ymin = Z2_lci,
      ymax = Z2_uci,
      fill = "interval"
    ),
    alpha = 0.2,
    colour = "red",
    linetype = 3
  ) +
  geom_line(
    data = \(df) dplyr::filter(df, model == "GeoMix"),
    aes(y = Z2_mean, colour = "geomix"),
    linewidth = 0.4
  ) +
  facet_wrap(vars(label), nrow = 4) +
  scale_x_reverse() +
  scale_colour_manual(
    values = c(
      actual = "black",
      geomix = "red"
    ),
    labels = c(
      actual = expression(atop("Actual cone tip","resistance " ~ q[c])),
      geomix = expression(atop("GeoMix predicted","cone tip resistance"))
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(interval = "red"),
    labels = c(
      interval = expression(atop("GeoMix 95%","prediction interval"))
    ),
    name = NULL
  ) +
  coord_flip() +
  theme_bw() +
  default_theme +
  theme(legend.key.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.7,"cm"),
        legend.spacing.x  = unit(0.3, "cm"),)+
  labs(
    x = "Depth (m)",
    y = expression("Cone tip resistance " ~ q[c]),
    subtitle = "*Distance to nearest training CPT"
  )

## 3.2 Small test CPT profiles qc ----
gtest_small <- pred_df %>%
  filter(loc_id %in% select_locID[c(1:4,13:16)]) %>%
  filter(model %in% c("GeoMix", "LM")) %>%
  group_by(loc_id) %>%
  mutate(
    dist = min(sqrt((x - train_locs$x)^2 + (y - train_locs$y)^2)),
    label = paste0(
      "CPT",
      str_remove(str_remove(name, "-TNO"), "IJV"),
      " - ",
      sprintf("%1.1f", dist),
      " km*"
    )
  ) %>%
  ungroup() %>%
  ggplot(aes(x = d)) +
  geom_line(
    aes(y = Z2, colour = "actual"),
    linewidth = 0.4
  ) +
  geom_ribbon(
    data = \(df) dplyr::filter(df, model == "GeoMix"),
    aes(
      ymin = Z2_lci,
      ymax = Z2_uci,
      fill = "interval"
    ),
    alpha = 0.2,
    colour = "red",
    linetype = 3
  ) +
  geom_line(
    data = \(df) dplyr::filter(df, model == "GeoMix"),
    aes(y = Z2_mean, colour = "geomix"),
    linewidth = 0.4
  ) +
  facet_wrap(vars(label), nrow = 2) +
  scale_x_reverse() +
  scale_colour_manual(
    values = c(
      actual = "black",
      geomix = "red"
    ),
    labels = c(
      actual = expression(atop("Actual cone tip","resistance " ~ q[c])),
      geomix = expression(atop("GeoMix predicted","cone tip resistance"))
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(interval = "red"),
    labels = c(
      interval = expression(atop("GeoMix 95%","prediction interval"))
    ),
    name = NULL
  ) +
  coord_flip() +
  theme_bw() +
  default_theme +
  theme(legend.key.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.7,"cm"),
        legend.spacing.x  = unit(0.3, "cm"),)+
  labs(
    x = "Depth (m)",
    y = expression("Cone tip resistance " ~ q[c]),
    subtitle = "*Distance to nearest training CPT"
  )

## 3.3 Save plots ----
ggsave("results/figures/ijv-test.pdf",gtest,width = 4.78, height = 4.8)
ggsave("results/figures/ijv-test-small.pdf",gtest_small,width = 4.78, height = 2.8)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Grid Plots ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Load predictions and parameters ----
GeoMix_predictions_mean <- readRDS(paste0(path,"gridpredictions/GeoMix_predictions_mean.rds"))
GeoMix_predictions_sd <- readRDS(paste0(path,"gridpredictions/GeoMix_predictions_sd.rds"))
params <- readRDS(paste0(path,"predictions/GeoMix_params.rds"))

# Attach predictions to GeoMix data frame
geomix_setup$df$Z2_mean <- GeoMix_predictions_mean
geomix_setup$df$Z2_sd <- GeoMix_predictions_sd
geomix_setup$df$Y1 <- params$params$Y1

# Join spatial coordinates
geomix_setup$df <- left_join(
  geomix_setup$df,
  distinct(full_df,xid,yid,easting=x,northing=y)
)

## 4.2 Setup spatial gridding ----
cell_grid <- readRDS('data/processed/cell_grid.rds')
line_df <- readRDS("data/processed/line_df.rds")

# Compute centred coordinates for plotting
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

## 4.3 Predicted qc (Z2) plot ----
g0 <- geomix_setup$df %>%
  dplyr::filter(d > 27, d < 48) %>%
  filter(d == floor(d)) %>%
  left_join(centred_locations) %>%
  ggplot() +
  geom_tile(aes(x = X, y = Y, fill = Z2_mean), width = 0.36, height = 0.36) +
  facet_wrap(vars(factor(d, levels = seq(27.5, 48, 0.5),
                         labels = paste0("", seq(27.5, 48, 0.5), "m BSL")))) +
  scale_fill_distiller(
    palette = "Spectral",
    label = ~round(inv_box_cox(.x)),
    breaks = box_cox(c(0, 5, 15, 30, 50, 75, 100)),
    guide = guide_colorbar(
      barwidth = unit(8, "cm"),
      barheight = unit(4, "mm")
    )
  ) +
  theme_bw() +default_theme+
  coord_fixed()+
  labs(x = "Easting (km)", y = "Northing (km)",
       fill = expression(atop("Net Cone Tip","Resistance"~q[c])),
       title = "")

## 4.4 Ground model (Z1) plot ----
gap <- 0.35
g1 <- geomix_setup$df %>%
  dplyr::filter(d > 27, d < 48) %>%
  left_join(centred_locations) %>%
  mutate(Y1 = factor(Z1,levels = 1:8, labels = c("GT1","GT2c","GT2","GT3","GT4","GT5","GT5*","GT6"))) %>%
  filter(floor(d) == d) %>%
  ggplot() +
  geom_tile(aes(x = X, y = Y, fill = Y1), width = 0.36, height = 0.36) +
  facet_wrap(vars(factor(d, levels = seq(27.5, 48, 0.5),
                         labels = paste0("", seq(27.5, 48, 0.5), "m BSL")))) +
  theme_bw() + default_theme + theme(plot.subtitle = element_text(hjust = 1))+
  scale_fill_manual(values = colours)+
  coord_fixed()+
  labs(x = "Easting (km)", y = "Northing (km)",
       subtitle = "BSL = Below Seal Level",
       fill = expression(atop("Ground Model","Geological Strata")))

## 4.5 MAP stratigraphy (Y1) plot ----
g2 <- geomix_setup$df %>%
  dplyr::filter(d > 27, d < 48) %>%
  left_join(centred_locations) %>%
  mutate(Y1 = factor(Y1,levels = 1:8, labels = c("GT1","GT2c","GT2","GT3","GT4","GT5","GT5*","GT6"))) %>%
  filter(floor(d) == d) %>%
  ggplot() +
  geom_tile(aes(x = X, y = Y, fill = Y1), width = 0.36, height = 0.36) +
  facet_wrap(vars(factor(d, levels = seq(27.5, 48, 0.5),
                         labels = paste0("", seq(27.5, 48, 0.5), "m BSL")))) +
  theme_bw() + default_theme +
  coord_fixed()+
  scale_fill_manual(values = colours)+
  labs(x = "Easting (km)", y = "Northing (km)",
       fill = expression(atop("Maximum A Posteriori","Geological Strata")))

## 4.6 Save  plots ----
ggsave("results/figures/ijv-z1.pdf",g1,width = 4.78, height = 5)
ggsave("results/figures/ijv-z2.pdf",g0,width = 4.78, height = 5)
ggsave("results/figures/ijv-y1.pdf",g2,width = 4.78, height = 5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. Cross-section plots  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_line_stack <- function(line_vals, show_legend = TRUE){
  default_theme <- theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25,"cm"),
    axis.text = element_text(size=6),
    axis.title = element_text(size=8),
    legend.title = element_text(size=8),
    legend.title.position = "top",
    legend.title.align = 0.5,
    legend.spacing.x  = unit(0.01, "cm"),
    plot.margin = margin(0, 0, 0, 0),
    strip.background = element_blank(),
    strip.placement = "inside",
    strip.switch.pad.wrap = unit(0.0, "cm"),
    strip.text = element_text(size = 6, margin = margin(0,0,0,0)),
    legend.text = element_text(size=6, margin = margin(0,0,0,1))
  )

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
      labels = c("GT1","GT2c","GT2","GT3","GT4","GT5","GT5*","GT6")
    )) %>%
    filter(line %in% line_vals)

  # PANEL 1: POSTERIOR MEAN qc
  g6p1 <- base_df %>%
    ggplot() +
    geom_tile(
      aes(x = d, y = r, fill = Z2_mean),
      width = 0.5,
      height = sqrt(2) * 0.35,
      show.legend = show_legend
    ) +
    geom_tile(
      aes(x = d, y = r, fill = Z2),
      data = . %>% drop_na(),
      width = 0.5,
      height = sqrt(2) * 0.35,
      show.legend = show_legend
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      data = . %>%
        filter(!is.na(Z2)) %>%
        group_by(line, xid, yid) %>%
        summarise(
          xmin = min(d) - sqrt(2)*0.35/2,
          xmax = max(d) + sqrt(2)*0.35/2,
          ymin = r[1] - 0.25,
          ymax = r[1] + 0.25
        ),
      colour = 'black',
      size = 0.1,
      fill = NA,
      show.legend = show_legend
    ) +
    facet_wrap(
      vars(factor(line, levels = 1:16, labels = paste0("Line ",1:16/2))),
      ncol = 1
    ) +
    scale_x_reverse() +
    coord_flip() +
    theme_bw() +
    default_theme +
    scale_fill_distiller(
      palette = "Spectral",
      label = ~round(inv_box_cox(.x)),
      breaks = box_cox(c(0, 5, 15, 30, 50, 75, 100)),
      guide = guide_colorbar(
        barwidth = unit(3.7, "cm"),
        barheight = unit(4, "mm")
      )
    ) +
    labs(
      x = "Depth (m) BSL",
      y = if_else(show_legend, "Distance along line (km)", ""),
      fill = expression(atop("Posterior Mean Estimate of",
                             "Cone Tip Resistance " ~ q[c]))
    )

  # PANEL 2: POSTERIOR STANDARD DEVIATION qc
  g6p2 <- base_df %>%
    ggplot() +
    geom_tile(
      aes(x = d, y = r, fill = Z2_sd),
      width = 0.5,
      height = sqrt(2) * 0.35,
      show.legend = show_legend
    ) +
    geom_tile(
      aes(x = d, y = r, fill = Z2_sd),
      data = . %>% drop_na(),
      width = 0.5,
      height = sqrt(2) * 0.35,
      show.legend = show_legend
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      data = . %>%
        filter(!is.na(Z2)) %>%
        group_by(line, xid, yid) %>%
        summarise(
          xmin = min(d) - sqrt(2)*0.35/2,
          xmax = max(d) + sqrt(2)*0.35/2,
          ymin = r[1] - 0.25,
          ymax = r[1] + 0.25
        ),
      colour = 'black',
      size = 0.1,
      fill = NA,
      show.legend = show_legend
    ) +
    facet_wrap(
      vars(factor(line, levels = 1:16, labels = paste0("Line ",1:16/2))),
      ncol = 1
    ) +
    scale_x_reverse() +
    coord_flip() +
    theme_bw() +
    default_theme +
    scale_fill_distiller(
      palette = "BrBG",
      guide = guide_colorbar(
        barwidth = unit(3.7, "cm"),
        barheight = unit(4, "mm")
      )
    ) +
    labs(
      x = "",
      y = if_else(show_legend, "Distance along line (km)", ""),
      fill = expression(atop("Posterior Standard Deviation of",
                             paste("Transformed ", q[c])))
    )

  # PANEL 3: MAP STRATIGRAPHY (Y1)
  g6p3 <- base_df %>%
    ggplot() +
    geom_tile(
      aes(x = d, y = r, fill = Y1),
      width = 0.5,
      height = sqrt(2) * 0.35,
      show.legend = show_legend
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      data = . %>%
        filter(!is.na(Z2)) %>%
        group_by(line, xid, yid) %>%
        summarise(
          xmin = min(d) - sqrt(2)*0.35/2,
          xmax = max(d) + sqrt(2)*0.35/2,
          ymin = r[1] - 0.25,
          ymax = r[1] + 0.25
        ),
      colour = 'black',
      size = 0.1,
      fill = NA,
      show.legend = show_legend
    ) +
    facet_wrap(
      vars(factor(line, levels = 1:16, labels = paste0("Line ",1:16/2))),
      ncol = 1
    ) +
    scale_x_reverse() +
    coord_flip() +
    theme_bw() +
    default_theme +
    scale_fill_manual(values = colours) +
    labs(
      x = "",
      y = if_else(show_legend, "Distance along line (km)", ""),
      fill = expression(atop("Maximum A Posteriori",
                             "Geological Strata"))
    )

  g6 <- g6p1 + g6p2 + g6p3
  return(g6)
}

## 5.1 Cross section plot Lines 1,2 ----
g6 <- make_line_stack(c(2,4),show_legend = T)
## 5.2 Cross section plot Lines 3-8 ----
g7 <- make_line_stack(c(3:8)*2,show_legend = T)

## 5.3 Cross section plot Lines 1-8 ----
g6nl <- make_line_stack(c(2,4),show_legend = F)
g8 <- (g6nl / g7) + plot_layout(heights = c(1, 3))

## 5.4 Save plots ----
ggsave("results/figures/ijv-cross.pdf",g8,width = 4.78, height = 6.8)
ggsave("results/figures/ijv-crossA.pdf",g6,width = 4.78, height = 3.2)
ggsave("results/figures/ijv-crossB.pdf",g7,width = 4.78, height = 6)
