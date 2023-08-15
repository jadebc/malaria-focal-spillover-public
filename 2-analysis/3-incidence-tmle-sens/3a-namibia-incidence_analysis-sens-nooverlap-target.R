################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence

# Sensitivity analysis
# Exclude cohorts with target areas
# that overlap with target areas

# Setting Q level to cohort to match 
# level selected for primary analysis 
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
data_human_list <- readRDS(namibia_human_process_sens_nooverlap_target_path)
data_mosq_list <- readRDS(namibia_mosq_process_sens_nooverlap_target_path)
data_hm_list <- readRDS(namibia_hm_process_sens_nooverlap_target_path)

# load covariate list ---------------------------------------
covarname <- readRDS(namibia_inc_covar_sens_nooverlap_target_path)
covarname_g = c("pre_spray_cover", "pre_incidence", "pre_rainfall",
                "pop_size_ea","pre_evi", "ea_elevation","surface_temp")
covarname_Q = covarname[!covarname %in% c(covarname_g)]

# SuperLearner set up ---------------------------------------
source(paste0(here::here(), '/2-analysis/0-namibia-incidence-learners.R'))

# Run hTMLE - Total effect ---------------------------------------

## Human ---------------------------------------
tic()
set.seed(123)
res_total_human <- run_hTMLE(
            df = data_human_list$short,
            effecttype = 'total',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
toc()

saveRDS(res_total_human, file = paste0(results_path, "namibia_htmle_inc_total_human_sens_nooverlap_target.RDS"))

## Mosquito ---------------------------------------
set.seed(123)
res_total_mosq <- 
  run_hTMLE(df = data_mosq_list$long,
            effecttype = 'total',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Mosquito")
toc()

saveRDS(res_total_mosq, file = paste0(results_path, "namibia_htmle_inc_total_mosq_sens_nooverlap_target.RDS"))

## Human & mosquito ---------------------------------------
set.seed(123)
res_total_hm_long <- run_hTMLE(df = data_hm_list$long,
            effecttype = 'total',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = Lrnr_mean$new(), # setting to mean to avoid unstable estimates
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human & mosquito",
            Vfolds = 30)

saveRDS(res_total_hm_long, file = paste0(results_path, "namibia_htmle_inc_total_hm_sens_nooverlap_target.RDS"))


# Run hTMLE - Spillover effect ---------------------------------------

## Human ---------------------------------------
tic()
set.seed(123)
res_spillover_human <- 
  run_hTMLE(df = data_human_list$short,
            effecttype = 'spillover',
            Qlevel = "individual", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
toc()

saveRDS(res_spillover_human, file = paste0(results_path, "namibia_htmle_inc_spillover_human_sens_nooverlap_target.RDS"))

## Mosquito ---------------------------------------
tic()
set.seed(123)
res_spillover_mosq <- 
  run_hTMLE(df = data_mosq_list$long,
            effecttype = 'spillover',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Mosquito")
toc()

saveRDS(res_spillover_mosq, file = paste0(results_path, "namibia_htmle_inc_spillover_mosq_sens_nooverlap_target.RDS"))

## Human & mosquito---------------------------------------
tic()
set.seed(123)
res_spillover_hm_long <- run_hTMLE(df = data_hm_list$long,
            effecttype = 'spillover',
            Qlevel = "individual", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g, # setting to mean to avoid unstable estimates
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human & mosquito",
            Vfolds = 30)

saveRDS(res_spillover_hm_long, file = paste0(results_path, "namibia_htmle_inc_spillover_hm_sens_nooverlap_target.RDS"))



# Run hTMLE - Direct effect ---------------------------------------

## Human ---------------------------------------
tic()
set.seed(123)
res_direct_human <- 
  run_hTMLE(df = data_human_list$short,
            effecttype = 'direct',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Human")
toc()

saveRDS(res_direct_human, file = paste0(results_path, "namibia_htmle_inc_direct_human_sens_nooverlap_target.RDS"))


## Mosquito ---------------------------------------
tic()
set.seed(123)
res_direct_mosq <- 
  run_hTMLE(df = data_mosq_list$long,
            effecttype = 'direct',
            Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
            glevel = "cohort", # always cohort for incidence, individual for prev
            dependency = TRUE, # TRUE if we want CIs adjusted for dependency
            unadj_est = FALSE, # TRUE if we want a unadjusted estimate
            y_l = 0,
            y_u = 1,
            covarname_g = covarname_g,    # set of names for adjustment variables
            covarname_Q = covarname_Q,
            SL_lib_Q = SL_lib_Q,
            SL_lib_g = SL_lib_g,
            SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
            SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
            verbose = FALSE,
            reservoir_name = "Mosquito")
toc()

saveRDS(res_direct_mosq, file = paste0(results_path, "namibia_htmle_inc_direct_mosq_sens_nooverlap_target.RDS"))


## Human & mosquito ---------------------------------------
# 
# ## DATA TOO SPARSE TO FIT MODELS
set.seed(123)
res_direct_hm_long <- run_hTMLE(df = data_hm_list$long,
                                   effecttype = 'direct',
                                   Qlevel = "cohort", # cohort/individual/adaptive for incidence, individual for prev
                                   glevel = "cohort", # always cohort for incidence, individual for prev
                                   dependency = TRUE, # TRUE if we want CIs adjusted for dependency
                                   unadj_est = FALSE, # TRUE if we want a unadjusted estimate
                                   y_l = 0,
                                   y_u = 1,
                                   covarname_g = covarname_g,    # set of names for adjustment variables
                                   covarname_Q = covarname_Q,
                                   SL_lib_Q = SL_lib_Q,
                                   SL_lib_g = SL_lib_g, # setting to mean to avoid unstable estimates
                                   SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
                                   SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
                                   verbose = FALSE,
                                   reservoir_name = "Human & mosquito",
                                   Vfolds = 30)

saveRDS(res_direct_hm_long, file = paste0(results_path, "namibia_htmle_inc_direct_hm_sens_nooverlap_target.RDS"))


