################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence
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
data_human_list <- readRDS(namibia_human_process_path)

# load covariate list ---------------------------------------
covarname <- readRDS(namibia_inc_covar_path)
covarname_g = c("pre_spray_cover", "pre_incidence", "pre_rainfall",
                "pop_size_ea","pre_evi", "ea_elevation","surface_temp")
covarname_Q = covarname[!covarname %in% c(covarname_g)]

# SuperLearner set up ---------------------------------------
source(paste0(here::here(), '/2-analysis/0-namibia-incidence-learners.R'))

# Run hTMLE - Total effect ---------------------------------------

res_total_human <- list()

## Human ---------------------------------------
tic()
set.seed(123)
res_total_human <- lapply(data_human_list, function(x)
  run_hTMLE(df = x,
            effecttype = 'total',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_Q = SL_lib_Q,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
)
toc()

names(res_total_human) = c("Long","Short")
saveRDS(res_total_human, file = paste0(results_path, "namibia_htmle_inc_total_human_cohort.RDS"))


# Run hTMLE - Spillover effect ---------------------------------------
res_spillover_human <- list()

## Human ---------------------------------------
tic()
set.seed(123)
res_spillover_human <- lapply(data_human_list, function(x)
  run_hTMLE(df = x,
            effecttype = 'spillover',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_Q = SL_lib_Q,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
)
toc()

names(res_spillover_human) = c("Long","Short")
saveRDS(res_spillover_human, file = paste0(results_path, "namibia_htmle_inc_spillover_human_cohort.RDS"))


# Run hTMLE - Direct effect ---------------------------------------
res_direct_human <- list()

## Human ---------------------------------------
tic()
set.seed(123)
res_direct_human <- lapply(data_human_list, function(x)
  run_hTMLE(df = x,
            effecttype = 'direct',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_Q = SL_lib_Q,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
)
toc()

names(res_direct_human) = c("Long","Short")
saveRDS(res_direct_human, file = paste0(results_path, "namibia_htmle_inc_direct_human_cohort.RDS"))





