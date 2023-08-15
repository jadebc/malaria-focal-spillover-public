################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence
# Data adaptively choose if Q is 
# individual or cohort level model 
################################################ 
rm(list=ls())
library(sl3)
library(tmle3)
library(origami)
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(GGally)

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
data_hm_list <- readRDS(namibia_hm_process_path)

# load covariate list ---------------------------------------
covarname <- readRDS(namibia_inc_covar_path)
covarname_g = c("pre_spray_cover", "pre_incidence", "pre_rainfall",
                "pop_size_ea")
covarname_Q = covarname[!covarname %in% c(covarname_g, "pre_evi", "ea_elevation","surface_temp")]

# SuperLearner set up ---------------------------------------
source(paste0(here::here(), '/2-analysis/0-namibia-incidence-learners.R'))

# Run hTMLE - Total effect ---------------------------------------


## Human & mosquito ---------------------------------------

res_total_hm <- list()
set.seed(123)
res_total_hm <- lapply(data_hm_list, function(x)
  choose_Qlevel(df = x,
          effecttype = 'total',
          dependency = TRUE, # TRUE if we want CIs adjusted for dependency
          unadj_est = FALSE, # TRUE if we want a unadjusted estimate
          y_l = 0,
          y_u = 1,
          covarname_Q = covarname_Q,    # set of names for adjustment variables
          SL_lib_Q = SL_lib_Q,
          SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
          SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
          verbose = FALSE,
          reservoir_name = "Human & mosquito",
          Vfolds = 30)
)

names(res_total_hm) = c("Long","Short")
saveRDS(res_total_hm, file = paste0(results_path, "namibia_htmle_inc_total_hm_adaptQ.RDS"))



