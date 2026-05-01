#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Libraries and helpers ----
library(patchwork)
library(slider)
library(ggplot2)
library(tidyverse)
select <- dplyr::select
box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)

colours <-c("gold","darkorange3","gray30", "forestgreen","darkblue","maroon", "red4", "dodgerblue4")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Simulate misspecification data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Load inputs ----
synthetic_list <- readRDS("data/simulation/synthetic_data.rds")

## 2.2 Matern-3/2 simulation function ----
sample_matern32_vec <- function(depths, mu, sigma2, ell, nugget = 0) {
  D <- as.matrix(dist(depths))
  sqrt3 <- sqrt(3)

  K <- sigma2 * (1 + sqrt3 * D / ell) * exp(-sqrt3 * D / ell)
  diag(K) <- diag(K) + nugget

  drop(chol(K) %*% rnorm(NROW(depths))) + mu
}

## 2.3 Generate CPT profile with misspecified stratigraphy ----
set.seed(1234)
data_base <- read_csv("data/application/misspec_cpt.csv") %>%
  left_join(data.frame(Y1 = 1:8, a0 = synthetic_list$alpha[,1], a1 = synthetic_list$alpha[,2],
            sigma2 = synthetic_list$sigma2, l = synthetic_list$lD)) %>%
  group_by(Y1) %>%
  mutate(Y2 = sample_matern32_vec(depths=d,mu = a0[1] + a1[1]*d,sigma2[1], ell = l[1], nugget = 0),
         Y2 = c(mean(Y2),Y2[-1]) + error/sd(Y2)) %>%
  ungroup() %>%
  mutate(Y2_2 = slide_dbl(Y2,mean,.before = 2,.after  = 2,.complete = TRUE),
         Y2 = if_else(is.na(Y2_2),Y2,Y2_2)) %>%
  mutate(Z2 = Y2 + rnorm(n(),0,sd=synthetic_list$tau)) %>%
  select(d,Y1,Z1,Y2,Z2) %>%
  mutate(Y1 = factor(Y1,levels = 1:8),Z1 = factor(Z1,levels = 1:8))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Boundaries and annotations ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Stratigraphy boundary locations ----
bnd_z1 <- data_base %>%
  filter(Z1!=1) %>%
  group_by(Z1) %>%
  summarise(ld = min(d), ud = max(d), d = ld - 0.05) %>%
  mutate(id = row_number())
bnd_y1 <-  data_base %>%
  filter(Y1!=1) %>%
  group_by(Y1) %>%
  summarise(ld = min(d), ud = max(d), d = ld - 0.05) %>%
  mutate(id = row_number())

delta <- 0.1/2
bar_locs <- c(90,115,120,145)

## 3.2 Text annotations for misspecification examples ----
ann1 <- filter(left_join(select(bnd_z1,Z1,d,ld,ud,id),select(bnd_y1,Y1,d)),is.na(Y1)) %>%
  mutate(label = c("→ prematurely  \ndetected boundary","→ failed to detect\nboundary until later")) %>%
  rename(class = Z1) %>%
  bind_rows(
    filter(bnd_y1,!(Y1 %in% bnd_z1$Z1)) %>%
      rename(class = Y1) %>% mutate(label = c("→ omitted a stratum","→ omitted another\nstratum") ,
                                    d = (ld+ud)/2 + c(0,1))
  ) %>% select(!Y1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Plot 1: Full misspec figure ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mis <- ggplot(data_base) +
  geom_rect(xmin = 0.5, xmax = 0.6,ymin = 0.1, ymax = 0.2,
            aes(fill = factor(3,levels =1:8))) +
  geom_rect(xmin = -Inf, xmax = Inf,ymin = -Inf, ymax = box_cox(14.5),
            fill = 'white') +
  geom_rect(xmin = -Inf, xmax = Inf,ymin = box_cox(75.5), ymax = Inf,
            fill = 'white') +
  geom_rect(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = 'white') +
  geom_rect(aes(xmin = d - delta, xmax = d + delta,
                ymin = box_cox(bar_locs[1]), ymax = box_cox(bar_locs[2]), fill = factor(Y1,levels = 1:8)),
            alpha = 0.8) +
  geom_rect(aes(xmin = d - delta, xmax = d + delta,
                ymin = box_cox(bar_locs[3]), ymax = box_cox(bar_locs[4]), fill = factor(Z1,levels = 1:8)),
            alpha = 0.8) +
  geom_segment(
    aes(x = d, xend = d, y = box_cox(15), yend = box_cox(145)),
    col= 'gray70', linewidth = 0.8, linetype = 3,
    data = bnd_z1
  ) +
  geom_segment(
    aes(x = d,xend=d,y=box_cox(15), yend=box_cox(115)),
    col= 'gray70', linewidth = 1, linetype = 1,
    data = bnd_y1
  ) +
  geom_segment(
    aes(x = d, xend = d, y = box_cox(15), yend = box_cox(145)),
    col= 'gray70', linewidth = 1, linetype = 1,
    data = drop_na(left_join(select(bnd_z1,Z1,d),select(bnd_y1,Y1,d)))
  ) +
  geom_label(aes(x = d, label = label, vjust = if_else(str_detect(label,"premature"),0.8,0.5)),
             y = box_cox(14), hjust = 1, data = ann1, size = 3,border.colour = 'NA')+
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='white',size=2) +
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='black',size=0.8) +
  annotate("text",x = -2, y= -Inf, label = " Problem:", size=4, hjust =0, fontface = "bold")+
  annotate("label",x = 0, y= -Inf, label = " The ground model...", size=3, hjust =0,border.colour = NA)+
  annotate("text",x = -2, y= box_cox(45), label = "CPT", size=5, fontface = 'bold')+
  annotate("text",x = -2, y= box_cox(mean(bar_locs[1:2])), label = "True", size=5, fontface = 'bold')+
  annotate("text",x = -2, y= box_cox(mean(bar_locs[3:4])), label = "GM ", size=5, fontface = 'bold')+
  coord_flip() +
  scale_x_reverse(breaks = seq(0,25,5)) +
  scale_fill_manual(values = colours, drop = FALSE) +
  scale_colour_manual(values = colours,drop = FALSE) +
  scale_y_continuous(label = ~round(inv_box_cox(.x)),breaks = box_cox(c(15,30,50,75)),limits = c(-4,box_cox(150)))+
  guides(fill = guide_legend(nrow = 1))+
  theme_minimal() +
  theme(
    legend.title.position = "top",
    legend.position = c(0.5,-0.3),
    legend.box = "vertical",
    axis.line = element_line(size=0.5),
    axis.ticks = element_line(size=0.5),
    legend.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 5, r = 0, b = 40, l = 5),
  ) +
  labs(
    x = expression("Depth " * d * " (m)"),
    y = expression("Cone tip resistance " * q[c]),
    fill = "Geological Strata",
    col = "",
    title = ""
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Plot 2: Simplified Y1/Z1 panels ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 5.1 Y1 (true stratigraphy) panel ----
mis2 <- ggplot(data_base) +
  geom_rect(xmin = 0.5, xmax = 0.6,ymin = 0.1, ymax = 0.2,
            aes(fill = factor(3,levels =1:8))) +
  geom_rect(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = 'white') +
  geom_rect(aes(xmin = d - delta, xmax = d + delta,
                ymin = box_cox(15), ymax = box_cox(75), fill = factor(Y1,levels = 1:8)),
            alpha = 0.8) +
  geom_segment(
    aes(x = d, xend = d, y = box_cox(9), yend =Inf),
    col= 'gray70', linewidth = 0.8, linetype = 3,
    data = bnd_z1
  ) +
  geom_segment(
    aes(x = d,xend=d,y=box_cox(9), yend=box_cox(75)),
    col= 'gray70', linewidth = 1, linetype = 1,
    data = bnd_y1
  ) +
  geom_segment(
    aes(x = d,xend=d,y=box_cox(9), yend=Inf),
    col= 'gray70', linewidth = 1, linetype = 1,
    data = drop_na(left_join(select(bnd_z1,Z1,d),select(bnd_y1,Y1,d)))
  ) +
  geom_label(aes(x = d, label = label), y = box_cox(9), hjust = 1, data = ann1, size = 3.5,border.colour = 'NA')+
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='white',size=2) +
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='black',size=0.8) +
  annotate("text",x = -2, y= box_cox(45), label = "Y1", size=5, fontface = 'bold')+
  coord_flip() +
  scale_x_reverse(breaks = seq(0,25,5)) +
  scale_fill_manual(values = colours, drop = FALSE) +
  scale_colour_manual(values = colours,drop = FALSE) +
  scale_y_continuous(label = ~round(inv_box_cox(.x)),breaks = box_cox(c(15,30,50,75)),limits = c(-10,box_cox(90)))+
  guides(fill = guide_legend(nrow = 1))+
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.75,-0.25),
    legend.title.position = "top",
    legend.box = "vertical",
    axis.line = element_line(size=0.5),
    axis.title.x = element_text(hjust = 1.8),
    axis.ticks = element_line(size=0.5),
    legend.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 5, r = 0, b = 30, l = 5)
  ) +
  labs(
    x = "Depth (m)",
    y = "Cone tip resistance (qc)",
    fill = "Geological Strata",
    col = "",
    title = "")

## 5.2 Z1 (ground model stratigraphy) panel ----
mis3 <- ggplot(data_base) +
  geom_rect(xmin = -Inf, xmax = Inf,ymin = -Inf, ymax = box_cox(14.5),
            fill = 'white') +
  geom_rect(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = 'white') +
  geom_rect(aes(xmin = d - delta, xmax = d + delta,
                ymin = box_cox(15), ymax = box_cox(75), fill = factor(Z1,levels = 1:8)),
            alpha = 0.8, show.legend = F) +
  geom_segment(
    aes(x = d, xend = d, y = -Inf, yend =box_cox(75)),
    col= 'gray70', linewidth = 0.8, linetype = 3,
    data = bnd_z1
  ) +
  geom_segment(
    aes(x = d, xend = d, y = -Inf, yend = box_cox(75)),
    col= 'gray70', linewidth = 1, linetype = 1,
    data = drop_na(left_join(select(bnd_z1,Z1,d),select(bnd_y1,Y1,d)))
  ) +
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='white',size=2) +
  geom_line(aes(x = d, y = Y2, col = "Actual qc"),col='black',size=0.8) +
  annotate("text",x = -2, y= box_cox(45), label = "Z1", size=5, fontface = 'bold')+
  coord_flip() +
  scale_x_reverse(breaks = seq(0,25,5)) +
  scale_fill_manual(values = colours, drop = FALSE) +
  scale_colour_manual(values = colours,drop = FALSE) +
  scale_y_continuous(label = ~round(inv_box_cox(.x)),breaks = box_cox(c(15,30,50,75)),limits = c(box_cox(10),box_cox(75)))+
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.ticks.x = element_line(size=0.5),
    legend.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 30, l = 0)
  ) +
  labs(y = "")

## 5.3 Combine panels ----
mis4 <-  mis2+mis3+
  plot_layout(widths = c(180, 80))

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Save figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

ggsave("results/figures/misspec.pdf",mis,width = 4.6,height = 3.6, device = cairo_pdf)
ggsave("results/figures/misspec2.pdf",mis4,width = 4.6,height = 4, device = cairo_pdf)
