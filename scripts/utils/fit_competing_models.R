#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

library(geowarp)
library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 run_comparisons function ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Fit GeoWarp / GP / LM competing models and save predictions ----
#
# Fits five competing models (GWwC, GW, GPwC, GP, LM) using the geowarp package
# with Vecchia approximation. Reads train/pred data from path, saves predictions.
#
# Arguments:
#   path           - results directory (e.g. "results/application/")
#   depth_interval - c(min, max) depth range for model fitting
run_comparisons <- function(path, depth_interval){
  source("scripts/utils/geowarp_builder.R")
  GWdata <- readRDS(paste0(path,'predictions/data.rds'))
  GWpred <- readRDS(paste0(path,'predictions/pred_df.rds'))

  ## 2.2 Loop over models ----
  for(model in c("GWwC","GW","GPwC","GP","LM")){
    cat("\n Model: Running Predictions for Model ",model,". \n --------------------------------\n")

    ### 2.2.1 Set mean formula and covariate flag ----
    if(str_detect(model,"wC") | str_detect(model,"LM")){
      mformula <- as.formula(paste0("~d+",paste0("SU",1:8,collapse = "+"),collapse = ""))
      wcov <- T
    }else{
      mformula <- as.formula("~d")
      wcov <- F
    }

    ### 2.2.2 Build model specification ----
    if(str_detect(model,"GW")){
      GWmodel <- model_factory(df = GWdata,
                               variable_name = "Z2",
                               m_formula = mformula,
                               min_depth = depth_interval[1], max_depth = depth_interval[2])
      cat(paste0("\n Fitting GeoWarp Model ",if_else(wcov,"with Z1 as covariate","without Z1 as covariate","\n")))
    }else if(str_detect(model,"GP")){
      GWmodel <- model_factory(df = GWdata,
                               variable_name = "Z2",
                               m_formula = mformula,
                               vert_warp = F,
                               geom_warp = F,
                               mbasis = F,
                               vbasis = F,
                               min_depth = depth_interval[1], max_depth = depth_interval[2])
      cat(paste0("\n Fitting Gaussian Process Model ",if_else(wcov,"with Z1 as covariate","without Z1 as covariate","\n")))
    }else if(str_detect(model,"LM")){
      GWmodel <- model_factory(df = GWdata,
                               variable_name = "Z2",
                               m_formula = mformula,
                               vert_warp = F,
                               geom_warp = F,
                               mbasis = F,
                               vbasis = F,
                               dev_model = "white",
                               min_depth = depth_interval[1], max_depth = depth_interval[2])
      cat(paste0("\n Fitting Linear Model ",if_else(wcov,"with Z1 as covariate","without Z1 as covariate","\n")))
    }

    ### 2.2.3 Optimise model ----
    fit <- geowarp_optimise(
      GWdata,
      GWmodel,
      n_parents = 150,
      best_of = 1,
      max_attempts = 10,
      grain_size = 32,
      threads = -1,
      trace = 1,
      control = list(
        maxit = 1000,
        prec = 0.001
      )
    )

    ### 2.2.4 Predict with Vecchia approximation ----
    if(str_detect(model,"G")){
      parent_structure <- vecchia_parent_structure(
        fit$observed_df,
        fit$model,
        observed_options = list(n_parents = 150 + 250),
        prediction_df = GWpred,
        prediction_type = 'regular',
        prediction_within_options = list(
          n_parents = 150
        ),
        prediction_between_options = list(
          n_parents = 250
        )
      )

      predictions <- predict(
        fit,
        GWpred,
        include_mean = TRUE,
        include_samples = TRUE,
        n_samples = 500,
        nugget = FALSE,
        parent_structure = parent_structure
      )
    }else{
      predictions <- predict(
        fit,
        GWpred,
        include_mean = TRUE,
        include_samples = TRUE,
        n_samples = 500,
        nugget = FALSE
      )
    }

    ### 2.2.5 Save predictions ----
    saveRDS(predictions,paste0(path,'predictions/',model,'_predictions.rds'))
  }
}
