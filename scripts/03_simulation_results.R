#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Setup and input data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Preliminaries  ----
source('scripts/utils/prediction_scoring.R')

## 1.2 Load data and parmeters ----
names <- c("GP","GW","LM","GPwC","GWwC","LGFM","GeoMix")
GeoMix_params <- readRDS(paste0(path,"predictions/GeoMix_params.rds"))
GeoMix_predictions <- readRDS(paste0(path,"predictions/GeoMix_predictions.rds"))
pred_df <- readRDS(paste0(path,"predictions/pred_df.rds")) %>%
  mutate(id = 'test')
train_df <- readRDS(paste0(path,"predictions/data.rds")) %>%
  mutate(id = 'train')
df <- arrange(bind_rows(train_df,pred_df),y,x,d)

## 1.3 Compute test evaluations ----
prediction_list <- lapply(names,eval_prediction_Y2,path=path)
prob_list <- c(list(compute_Z1(path)),lapply(c("LGFM","GeoMix"),eval_probability_Y1,path=path))

## 1.4 Assemble test statistics ----
stats_df <- bind_rows(map(prediction_list,~.x$stats))
pred_df <- bind_rows(map(prediction_list,~.x$pred_df))
prob_stats_df <- bind_rows(map(prob_list,~.x$stats))

## 1.5 Quick check Z1 Y1hat agreement ----
filter(prob_list[[3]]$prob_df,id == "test") %>%
  summarise(mean(Y1hat == Z1))

#%%%%%%%%%%%%%%%%%%
# 2 Plot theme ----
#%%%%%%%%%%%%%%%%%%

## 2.1 General theme ----
default_theme <-   theme(legend.position = "bottom",
                         legend.title.position = "top",
                         axis.text = element_text(size=8),
                         axis.title = element_text(size=9),
                         legend.title = element_text(size=9),
                         plot.title = element_text(size=9),
                         plot.subtitle = element_text(size=8),
                         legend.title.align = 0.5,
                         legend.spacing.x  = unit(0.01, "cm"),
                         legend.margin = margin(0,0,0,-20),
                         plot.margin = margin(0, 0, 0, 0),
                         strip.background = element_blank(),
                         strip.placement = "inside",   # optional
                         strip.switch.pad.wrap = unit(0.0, "cm"),
                         strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                         legend.text = element_text(size=9))

## 2.2 Choose an example x-y profile ----
xex <- 8; yex <- 7

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Setup Predictions Dataframe  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Construct  prediction data frame ----
Y2pred <- cbind(pred_df,as.data.frame(GeoMix_predictions$samples)) %>%
  pivot_longer(starts_with("V"),names_prefix = "V",
               names_to = "Iter",values_to = "Y2hat") %>%
  group_by(x,y,d,Y2) %>%
  summarise(Y2hat_m = mean(Y2hat), Y2hat_sd = sd(Y2hat)) %>%
  bind_rows(mutate(select(train_df,x,y,d,Y2,Z2),Y2hat_m = Z2,
                   Y2hat_sd = sqrt(GeoMix_params$params$tau2))) %>%
  rename(Y2hat = Y2hat_m) %>%
  filter(y == yex)

## 3.2 Extract Y1 probabilities ----
Y1probs <- prob_list[[3]]$prob_df %>%
  mutate(Y1hat_prob = do.call(pmax, c(select(., starts_with("X",ignore.case = F)), na.rm = TRUE))) %>%
  select(ID=sample_index,x,y,d,Y1,Z1,flag=id,Y1hat,Y1hat_prob) %>%
  mutate(ID = as.integer(ID),Y1 = as.numeric(Y1), Z1 = as.numeric(Z1), Y1hat = as.numeric(Y1hat))%>%
  filter(y == yex)

## 3.3 Merge data sets ----
full_pred_df <- left_join(Y1probs,Y2pred) %>%
  left_join(Y1probs)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Plot 1 - Cross-section panels ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Plot latent, observed, MAP strata ----
g1 <-
  full_pred_df %>%
  select(x,y,d,Y1_val = Y1,Z1_val = Z1,Y1hat_val = Y1hat,Y1hat_prob) %>%
  mutate(Y1_prob = 1, Z1_prob = 1) %>%
  pivot_longer(!c(x,y,d),names_to = c("facet","type"), names_sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(facet = factor(facet,levels = c("Y1","Z1","Y1hat"),
                        labels = c("Latent Geological Strata",
                                   "Ground Model Strata",
                                   "MAP Geological Strata"))) %>%
  ggplot()+
  geom_tile(aes(x = x, y = d, fill = factor(val)),show.legend = T)+
  scale_fill_manual(values = colours)+
  scale_y_reverse(limits = c(-0.5,20))+theme_bw()+default_theme+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4,"cm"),
        legend.key.spacing = unit(0.05,"cm"),
        plot.margin = margin(0, 0, 0, 0),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = guide_legend(nrow=2), alpha = "none")+
  facet_wrap(vars(facet),ncol = 1)+
  labs(x="",y="Depth (m)",
       fill="Geological\nStrata")

## 4.2 Plot latent, observed, and predicted geotechnical ----
g2 <-  full_pred_df %>%
  select(x,y,d,Y2,Z2,Y2hat) %>%
  pivot_longer(!c(x,y,d),names_to = c("facet")) %>%
  mutate(facet = factor(facet,levels = c("Y2","Z2","Y2hat"),
                        labels = c("Latent Geotechnical Property",
                                   "Geotechnical Measurements",
                                   "Predicted Geotechnical Property"))) %>%
  drop_na() %>%
  ggplot()+
  geom_tile(aes(x = x, y = d, fill = value),show.legend = T)+
  scale_fill_distiller(palette = "Spectral")+
  scale_y_reverse(limits = c(-0.5,20))+theme_bw()+ default_theme+
  theme(legend.position = "bottom",legend.key.size = unit(0.5,"cm"),
        legend.key.spacing = unit(0.05,"cm"),
        plot.margin = margin(0, 0, 0, 0),axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_wrap(vars(facet),ncol = 1)+
  labs(x="",y="",
       fill="Geotechnics")

## 4.3 Probability of MAP ----
g3 <-
  full_pred_df %>%
  ggplot()+
  geom_tile(aes(x = x, y = d, fill = Y1hat_prob),show.legend = T)+
  scale_fill_distiller(palette = "Blues", breaks = c(0.6,0.8,1))+
  scale_y_reverse(limits = c(-0.5,20))+theme_bw()+default_theme+
  facet_wrap(vars("Probability of the MAP Strata"),ncol = 1)+
  theme(legend.position = "bottom", legend.key.size = unit(0.5,"cm"),
        legend.key.spacing = unit(0.05,"cm"))+
  labs(x="Distance along x axis",y="",
       fill="MAP\nprobability")

## 4.4 Posterior standard deviation ----
g4 <-  full_pred_df %>%
  ggplot()+
  geom_tile(aes(x = x, y = d, fill = Y2hat_sd),show.legend = T)+
  scale_fill_distiller(palette = "BrBG")+
  scale_y_reverse(limits = c(-0.5,20))+theme_bw()+ default_theme+
  theme(legend.position = "bottom",legend.key.size = unit(0.5,"cm"),
        legend.key.spacing = unit(0.05,"cm"))+
  facet_wrap(vars("Std Deviation of Predictions"),ncol = 1)+
  labs(x="Distance along x axis",y="",
       fill="Standard\ndeviation")

## 4.5 Extract legends (vertical, for right-column layout) ----
leg_theme <- theme(
  legend.position       = "right",
  legend.box            = "vertical",
  legend.title          = element_text(hjust = 0.5, size = 8),
  legend.title.position = "top",
  legend.key.width      = unit(0.4, "cm")
)
get_leg_vert <- function(p) cowplot::get_legend(p + leg_theme)
leg_strata <- get_leg_vert(g1 + guides(fill = guide_legend(ncol = 2, keywidth = unit(0.4, "cm"))))
leg_prop   <- get_leg_vert(g2 + guides(fill = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(1.3, "cm"))))
leg_prob   <- get_leg_vert(g3 + guides(fill = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(1.2, "cm"))) +
                theme(legend.title = element_text(hjust = 0.5, size = 8, margin = margin(b = 4))))
leg_sd     <- get_leg_vert(g4 + guides(fill = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(1.4, "cm"))))

## 4.6 Remove legends ----
g1n <- g1 + theme(legend.position = "none")
g2n <- g2 + theme(legend.position = "none")
g3n <- g3 + theme(legend.position = "none")
g4n <- g4 + theme(legend.position = "none")

## 4.7  Format plots ----
p_leg_strata <- wrap_elements(leg_strata)
p_leg_prop   <- wrap_elements(leg_prop)
p_leg_prob   <- wrap_elements(leg_prob)
p_leg_sd     <- wrap_elements(leg_sd)

## 4.8  Make columns and assemble plots ----
left_col  <- ((g1n / g3n) + plot_layout(heights = c(4, 1)))
right_col <- ((g2n / g4n) + plot_layout(heights = c(4, 1)))
leg_col   <- (p_leg_strata / plot_spacer() / p_leg_prop / plot_spacer() / p_leg_prob / plot_spacer() / p_leg_sd) +
  plot_layout(heights = c(0.95, 0.05, 0.95, 0.05, 0.95, 0.1, 0.95))

## 4.9  Final plot ----
syn_cross <-
  ((left_col | right_col | plot_spacer() | leg_col) +
     plot_layout(widths = c(5, 5, 0, 2)))

## 4.10  Save plot ----
ggsave("results/figures/syn-cross.pdf",syn_cross,width=4.78,height=3.5,device=cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Plot 2 - Posterior samples.       ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 5.1 Process posterior data sets ----
GM_sampY2 <- cbind(pred_df,as.data.frame(GeoMix_predictions$samples)) %>%
  filter(x== xex,y ==yex) %>%
  pivot_longer(starts_with("V"),names_prefix = "V",names_to = "Iter",values_to = "Y2hat")
GM_sampY1 <- cbind(df,as.data.frame(t(GeoMix_params$samples$Y1))) %>%
  filter(x== xex,y ==yex,id == "test") %>%
  pivot_longer(starts_with("V"),names_prefix = "V",names_to = "Iter",values_to = "Y1hat")
GM_samp <- left_join(GM_sampY1,GM_sampY2)

## 5.2 Posterior sample plot ----
syn_plot <- GM_samp %>%
  filter(d %in% seq(1,20,2)) %>%
  ggplot()+
  geom_histogram(aes(x=Y2hat,fill=factor(Y1hat,levels = 1:8), y =after_stat(count)),alpha = 0.75, position = "stack")+
  geom_vline(aes(xintercept = Y2,col = factor(Y1,levels = 1:8)), data = filter(pred_df,x== xex,y ==yex,d %in% seq(1,20,2)),size=1.5)+
  facet_wrap(vars(factor(d,levels = 1:20, label = paste0("",1:20,"m Depth"))),ncol=5)+
  theme_bw()+
  guides(fill = guide_legend(nrow=1), col = guide_legend(nrow=1))+
  labs(x="Predicted Latent Cone Tip Resistance ",y="Stacked Frequency",
       fill=expression("Geological Strata "*Y[1](bold(s))*" and "*Y[2](bold(s))*" (line)"),
       col=expression("Geological Strata "*Y[1](bold(s))*" and "*Y[2](bold(s))*" (line)")) +
  default_theme+theme(legend.key.size = unit(0.4,"cm"), legend.position = c(0.5,-0.39),
                      legend.title = element_text(hjust = 0.5, margin = margin(b = 0)),
                      plot.margin = margin(0, 0, 30, 0))

## 5.3 Save plot ----
ggsave("results/figures/syn-post.pdf",syn_plot,width=4.78,height=2.3,device=cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Plot 3 - Parameter recovery ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 Collate posterior samples ----
param_samples <- data.frame(GeoMix_params$samples$a0) %>%
  mutate(iter = 1:n(), param = "alpha1") %>%
  pivot_longer(starts_with("X"),names_to = "k", values_to = "estimate", names_prefix = "X") %>%
  bind_rows(
    data.frame(GeoMix_params$samples$a1) %>%
      mutate(iter = 1:n(), param = "alpha2") %>%
      pivot_longer(starts_with("X"),names_to = "k", values_to = "estimate", names_prefix = "X")
  ) %>%
  bind_rows(
    data.frame(GeoMix_params$samples$sigma2) %>%
      mutate(iter = 1:n(), param = "sigma2") %>%
      pivot_longer(starts_with("X"),names_to = "k", values_to = "estimate", names_prefix = "X")
  ) %>%
  bind_rows(
    data.frame(GeoMix_params$samples$lD) %>%
      mutate(iter = 1:n(), param = "lD") %>%
      pivot_longer(starts_with("X"),names_to = "k", values_to = "estimate", names_prefix = "X")
  )%>%
  bind_rows(
    data.frame(GeoMix_params$samples$lL) %>%
      mutate(iter = 1:n(), param = "lL") %>%
      pivot_longer(starts_with("X"),names_to = "k", values_to = "estimate", names_prefix = "X")
  )%>%
  bind_rows(
    data.frame(estimate = GeoMix_params$samples$tau2) %>%
      mutate(iter = 1:n(), param = "tau2", k="1")
  )

## 6.2 True parameters ----
true_params <- data.frame(param = c(rep(c("alpha1","alpha2"),each = 8), rep("sigma2",8), rep("lL",8),rep("lD",8),"tau2"),
                          k = c(rep(1:8,5),NA),
                          value = c(c(10,17,14,19,7,11,6,12,0.2,0.05,0.2,0.05,0.1,0.15,0.2,0.15),c(5,4,6,5,10,4,3,7),
                                    c(0.5,0.5,0.6,0.4,0.6,0.5,0.7,0.3),c(1,0.8,0.6,0.7,0.6,1,0.7,0.5)+1,0.25))

## 6.3 Plot samples ----
syn_param <- param_samples %>%
  filter(param!="tau2") %>%
  ggplot()+
  geom_histogram(aes(x=estimate,fill = k),alpha = 0.75,show.legend = F)+
  geom_vline(aes(xintercept = value), colour = 'black',show.legend = F,
             data = filter(true_params,param!="tau2"))+
  facet_grid(factor(k,levels=1:8,labels = paste0("Class ",1:8))
             ~factor(param,levels=c("alpha1","alpha2","sigma2","lL","lD"),
                     labels = c("α₁","α₂","σ²","l₁","l₂")),scales = "free_x")+
  scale_x_continuous(breaks = scales::breaks_extended(n = 4))+
  labs(y = "Frequency",x="Estimate Value",title = "Posterior Parameter Samples",subtitle = "Actual value is indicated on the plot as a vertical line.")+
  theme_bw()+theme(panel.grid.minor = element_blank(),strip.text = element_text(size=10),
                   plot.title = element_text(size = 14))

## 6.4 Save plot ----
ggsave("results/figures/syn-params.pdf",syn_param,width=4.78,height=6.5,device=cairo_pdf)
