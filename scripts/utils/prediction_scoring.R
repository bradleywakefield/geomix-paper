#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Libraries ----
library(patchwork)
library(cowplot)
library(scoringRules)
library(tidyverse)
library(reshape2)
library(ggh4x)
library(scales)
library(grid)

select <- dplyr::select

## 1.2 Shared helpers and colour palette ----
box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)
get_mode <- function(x) names(sort(table(x), decreasing = TRUE))[1]
get_mode_prob <- function(x) max(table(x))/length(x)
colours <-c("gold","darkorange3", "forestgreen","darkblue","maroon", "red4", "dodgerblue4","gray30")
get_leg <- function(p) {
  legs <- cowplot::get_plot_component(
    p + theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    ),
    "guide-box",
    return_all = TRUE
  )
  
  legs[[which.max(vapply(legs, function(x) length(x$grobs), integer(1)))]]
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Scoring functions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Block diagonal scoring rule (DSS at lag p) ----
dss_p_score <- function(observed, posterior_samples, p) {
  sapply(seq_len(floor(length(observed) / p)), function(i) {
    indices <- p * (i - 1) + seq_len(p)
    samples <- posterior_samples[indices, ]
    mu <- rowMeans(samples)
    Sigma <- cov(t(samples))+diag(length(mu))*1e-6
    chol_Sigma <- chol(Sigma)
    (
      sum(log(diag(chol_Sigma)))
      + 0.5 * sum(
        backsolve(chol_Sigma, observed[indices] - mu, transpose = TRUE) ^ 2
      )
    )
  })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 IJV application evaluation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Evaluate cone tip resistance predictions (Z2) ----
eval_prediction_Z2 <- function(name, path = "results/application/"){
  require(scoringRules)
  pred_df <- readRDS(paste0(path,"predictions/pred_df.rds"))
  predictions <- readRDS(paste0(path,"predictions/",name,"_predictions.rds"))
  pred_df$Z2 <-  box_cox(pred_df$qc)
  na_vals <- which(!is.na(pred_df$qc))
  pred_df <- pred_df[na_vals,]
  predictions$samples <- (predictions$samples[na_vals,])
  predictions$samples[is.nan(predictions$samples)] <- 0
  pred_df$Z2_mean <- apply(predictions$samples,1,mean)
  pred_df$Z2_med <- apply(predictions$samples,1,median)
  pred_df$Z2_sd <- apply(predictions$samples,1,sd)
  pred_df$Z2_lci <- apply(predictions$samples,1,quantile,1-0.975)
  pred_df$Z2_uci <- apply(predictions$samples,1,quantile,0.975)

  pred_df$Z2_crps <- crps_sample(y =pred_df$Z2, dat = predictions$samples)
  pred_df$Z2_logS <- logs_sample(y = pred_df$Z2, dat = predictions$samples)
  pred_df$Z2_dss <- dss_sample(y = pred_df$Z2, dat = predictions$samples)
  pred_df$Z2_int <- ints_sample(y = pred_df$Z2, dat = predictions$samples, 0.95)
  pred_df$cov <- (pred_df$Z2 > pred_df$Z2_lci) & (pred_df$Z2 < pred_df$Z2_uci)
  pred_df$bias <- pred_df$Z2_mean - pred_df$Z2
  pred_df$sqdev <- (pred_df$Z2_mean - pred_df$Z2)^2

  dss2_value <- mean(unlist(lapply(unique(pred_df$loc_id),
                                   function(i){
                                     indices <- which(pred_df$loc_id==i)
                                     score <- dss_p_score(pred_df$Z2[indices],predictions$samples[indices,],p=2)
                                   })))

  ss_tot <- sum((pred_df$Z2 - mean(pred_df$Z2))^2)
  ss_res <- sum((pred_df$Z2 - pred_df$Z2_mean)^2)
  R2 <- 1 - ss_res/ss_tot

  stats <- data.frame(model = name,
                      RMSE = mean(pred_df$sqdev), coverage = mean(pred_df$cov),
                      R2 = R2, bias = mean(pred_df$bias), int_score = mean(pred_df$Z2_int),
                      CRPS = mean(pred_df$Z2_crps), DSS = mean(pred_df$Z2_dss), DSS2 = dss2_value,
                      logS = mean(pred_df$Z2_logS))

  return(list(stats = stats,pred_df = data.frame(model=name,pred_df)))
}

## 3.2 Evaluate stratigraphy classification probabilities (Z1) ----
eval_probability_Z1 <- function(name, path = "results/application/"){
  require(yardstick)
  full_df <- readRDS(paste0(path,"predictions/full_df.rds"))
  pred_df <- readRDS(paste0(path,"predictions/pred_df.rds"))
  params <- readRDS(paste0(path,"predictions/",name,"_params.rds"))
  K <- length(params$params$sigma2)
  Y1_prob <- data.frame(select(full_df,x,y,d,Z1,Z2,loc_id),
                        Y1_id = rownames(params$params$Y1prob),
                        params$params$Y1prob,
                        Y1hat = factor(params$params$Y1,levels = 1:K), model = name) %>%
    mutate(sample_index = str_extract(str_remove(Y1_id,"Y1"),"\\d+"),.before = everything()) %>%
    select(!Y1_id) %>% mutate(Z1 = factor(Z1,levels = 1:K)) %>%
    mutate(id = if_else(loc_id %in% unique(pred_df$loc_id) | is.na(Z2),"test","train"))

  stats_all <- data.frame(Subset = "All",
                          class_accu = accuracy(Y1_prob,Z1,Y1hat)$.estimate,
                          class_sens = sensitivity(Y1_prob,Z1,Y1hat)$.estimate,
                          class_spec = specificity(Y1_prob,Z1,Y1hat)$.estimate,
                          brier = brier_class(Y1_prob, Z1,paste0("X",1:K))$.estimate,
                          loss = mn_log_loss(Y1_prob, Z1,paste0("X",1:K))$.estimate,
                          auc = roc_auc(Y1_prob, Z1,paste0("X",1:K))$.estimate)

  stats_train <- data.frame(Subset = "Train",
                            class_accu = accuracy(filter(Y1_prob,id=="train"),Z1,Y1hat)$.estimate,
                            class_sens = sensitivity(filter(Y1_prob,id=="train"),Z1,Y1hat)$.estimate,
                            class_spec = specificity(filter(Y1_prob,id=="train"),Z1,Y1hat)$.estimate,
                            brier = brier_class(filter(Y1_prob,id=="train"), Z1,paste0("X",1:K))$.estimate,
                            loss = mn_log_loss(filter(Y1_prob,id=="train"), Z1,paste0("X",1:K))$.estimate,
                            auc = roc_auc(filter(Y1_prob,id=="train"), Z1,paste0("X",1:K))$.estimate)

  stats_test <- data.frame(Subset = "Test",
                           class_accu = accuracy(filter(Y1_prob,id=="test"),Z1,Y1hat)$.estimate,
                           class_sens = sensitivity(filter(Y1_prob,id=="test"),Z1,Y1hat)$.estimate,
                           class_spec = specificity(filter(Y1_prob,id=="test"),Z1,Y1hat)$.estimate,
                           brier = brier_class(filter(Y1_prob,id=="test"), Z1,paste0("X",1:K))$.estimate,
                           loss = mn_log_loss(filter(Y1_prob,id=="test"), Z1,paste0("X",1:K))$.estimate,
                           auc = roc_auc(filter(Y1_prob,id=="test"), Z1,paste0("X",1:K))$.estimate)

  stats <- mutate(bind_rows(stats_train,stats_test,stats_all),Model = name,.before =everything())

  boundaries <- Y1_prob %>%
    select(loc_id,d,Y1hat,Z1) %>%
    pivot_longer(!c(loc_id,d), names_to = "bound", values_to = "class") %>%
    group_by(bound,loc_id,class) %>%
    summarise(ldepth = min(d), udepth = max(d)) %>%
    left_join(distinct(full_df,loc_id,x,y))

  return(list(stats = stats,prob_df = Y1_prob, boundaries = boundaries))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Simulation study evaluation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Evaluate latent geotechnical predictions (Y2) ----
eval_prediction_Y2 <- function(name, path = "results/simulation/"){
  require(scoringRules)
  pred_df <- readRDS(paste0(path,"predictions/pred_df.rds"))
  predictions <- readRDS(paste0(path,"predictions/",name,"_predictions.rds"))
  pred_df$Z2_mean <- predictions$mean
  pred_df$Z2_med <- apply(predictions$samples,1,median)
  pred_df$Z2_sd <- apply(predictions$samples,1,sd)
  pred_df$Z2_lci <- apply(predictions$samples,1,quantile,1-0.975)
  pred_df$Z2_uci <- apply(predictions$samples,1,quantile,0.975)
  pred_df$Z2_crps <- crps_sample(y = pred_df$Y2, dat = predictions$samples)
  pred_df$Z2_logS <- logs_sample(y = pred_df$Y2, dat = predictions$samples)
  pred_df$Z2_dss <- dss_sample(y = pred_df$Y2, dat = predictions$samples)
  pred_df$Z2_int <- ints_sample(y = pred_df$Y2, dat = predictions$samples, 0.95)
  pred_df$cov <- (pred_df$Y2 > pred_df$Z2_lci) & (pred_df$Y2 < pred_df$Z2_uci)
  pred_df$bias <- pred_df$Z2_mean - pred_df$Y2
  pred_df$sqdev <- (pred_df$Z2_mean - pred_df$Y2)^2

  ss_tot <- sum((pred_df$Y2 - mean(pred_df$Y2))^2)
  ss_res <- sum((pred_df$Y2 - pred_df$Z2_mean)^2)
  R2 <- 1 - ss_res/ss_tot

  dss2_value <- mean(unlist(lapply(1:400,
                                   function(i){
                                     indices <- which(pred_df$locID==i)
                                     score <- dss_p_score(pred_df$Y2[indices],predictions$samples[indices,],p=2)
                                   })))

  stats <- data.frame(model = name,
                      RMSE = mean(pred_df$sqdev), coverage = mean(pred_df$cov),
                      R2 = R2, bias = mean(pred_df$bias), int_score = mean(pred_df$Z2_int),
                      CRPS = mean(pred_df$Z2_crps), DSS = mean(pred_df$Z2_dss), DSS2 = dss2_value,
                      logS = mean(pred_df$Z2_logS))

  return(list(stats = stats,pred_df = data.frame(model=name,pred_df)))
}

## 4.2 Evaluate latent stratigraphy classification (Y1) ----
eval_probability_Y1 <- function(name, path = "results/simulation/"){
  require(yardstick)
  train_df <- readRDS(paste0(path,"predictions/data.rds")) %>%
    mutate(id = 'train')
  pred_df <- readRDS(paste0(path,"predictions/pred_df.rds")) %>%
    mutate(id = 'test')
  df <- arrange(bind_rows(train_df,pred_df),y,x,d)
  params <- readRDS(paste0(path,"predictions/",name,"_params.rds"))
  K <- length(params$params$sigma2)
  Y1_prob <- data.frame(select(df,x,y,d), Y1_id = rownames(params$params$Y1prob),
                        params$params$Y1prob,Y1 = factor(df$Y1), Z1 = factor(df$Z1), Y1hat = factor(params$params$Y1,levels = 1:K), id = df$id, model = name) %>%
    mutate(sample_index = str_extract(str_remove(Y1_id,"Y1"),"\\d+"),.before = everything()) %>%
    select(!Y1_id)

  stats_all <- data.frame(Subset = "All",
                          class_accu = accuracy(Y1_prob,Y1,Y1hat)$.estimate,
                          class_sens = sensitivity(Y1_prob,Y1,Y1hat)$.estimate,
                          class_spec = specificity(Y1_prob,Y1,Y1hat)$.estimate,
                          brier = brier_class(Y1_prob, Y1,paste0("X",1:K))$.estimate,
                          loss = mn_log_loss(Y1_prob, Y1,paste0("X",1:K))$.estimate,
                          auc = roc_auc(Y1_prob, Y1,paste0("X",1:K))$.estimate)

  stats_train <- data.frame(Subset = "Train",
                            class_accu = accuracy(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            class_sens = sensitivity(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            class_spec = specificity(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            brier = brier_class(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate,
                            loss = mn_log_loss(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate,
                            auc = roc_auc(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate)

  stats_test <- data.frame(Subset = "Test",
                           class_accu = accuracy(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           class_sens = sensitivity(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           class_spec = specificity(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           brier = brier_class(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate,
                           loss = mn_log_loss(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate,
                           auc = roc_auc(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate)

  boundaries <- Y1_prob %>%
    select(x,y,d,Y1hat,Y1,Z1) %>%
    pivot_longer(!c(x,y,d), names_to = "bound", values_to = "class") %>%
    group_by(bound,x,y,class) %>%
    summarise(ldepth = min(d), udepth = max(d)) %>%
    left_join(distinct(df,x,y))

  stats <- mutate(bind_rows(stats_train,stats_test,stats_all),Model = name,.before =everything())
  return(list(stats = stats,prob_df = Y1_prob,boundaries=boundaries))
}

## 4.3 Ground model (Z1) as baseline classifier ----
compute_Z1 <- function(path = "results/simulation/"){
  require(yardstick)
  train_df <- readRDS(paste0(path,"predictions/data.rds")) %>%
    mutate(id = 'train')
  pred_df <- readRDS(paste0(path,"predictions/pred_df.rds")) %>%
    mutate(id = 'test')
  df <- arrange(bind_rows(train_df,pred_df),y,x,d)
  K <- max(train_df$Y1)
  Z1_prob <- data.frame(model.matrix(~-1+Z1,data=mutate(df,Z1 = factor(Z1,levels = 1:K))))
  names(Z1_prob) <- paste0("X",1:K)
  Y1_prob <- data.frame(select(df,x,y,d),Z1_prob,Y1 = factor(df$Y1), Y1hat = factor(df$Z1,levels = 1:K), id = df$id, model = "Z1") %>%
    mutate(sample_index = 1:n(),.before = everything())

  stats_all <- data.frame(Subset = "All",
                          class_accu = accuracy(Y1_prob,Y1,Y1hat)$.estimate,
                          class_sens = sensitivity(Y1_prob,Y1,Y1hat)$.estimate,
                          class_spec = specificity(Y1_prob,Y1,Y1hat)$.estimate,
                          brier = brier_class(Y1_prob, Y1,paste0("X",1:K))$.estimate,
                          loss = mn_log_loss(Y1_prob, Y1,paste0("X",1:K))$.estimate,
                          auc = roc_auc(Y1_prob, Y1,paste0("X",1:K))$.estimate)

  stats_train <- data.frame(Subset = "Train",
                            class_accu = accuracy(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            class_sens = sensitivity(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            class_spec = specificity(filter(Y1_prob,id=="train"),Y1,Y1hat)$.estimate,
                            brier = brier_class(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate,
                            loss = mn_log_loss(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate,
                            auc = roc_auc(filter(Y1_prob,id=="train"), Y1,paste0("X",1:K))$.estimate)

  stats_test <- data.frame(Subset = "Test",
                           class_accu = accuracy(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           class_sens = sensitivity(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           class_spec = specificity(filter(Y1_prob,id=="test"),Y1,Y1hat)$.estimate,
                           brier = brier_class(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate,
                           loss = mn_log_loss(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate,
                           auc = roc_auc(filter(Y1_prob,id=="test"), Y1,paste0("X",1:K))$.estimate)

  stats <- mutate(bind_rows(stats_train,stats_test,stats_all),Model = "Z1",.before =everything())
  return(list(stats = stats,prob_df = Y1_prob))
}
