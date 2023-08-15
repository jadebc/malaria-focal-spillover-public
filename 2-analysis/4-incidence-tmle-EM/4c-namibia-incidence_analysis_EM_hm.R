################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence

# All unadjusted due to data sparsity

# All cohort-level Q model to be consistent with
# primary analysis 
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

data_hm_list$long$sex = ifelse(data_hm_list$long$sex == "Male", 1, 0)

# define modifier list  ---------------------------------------
modifier_list <- c("pre_incidence_abovemed", "tx_cov_cohort_abovemed", "pre_spray_cover_abovemed",
  "pre_rainfall_abovemed", "pre_evi_abovemed", "ea_elevation_abovemed", 
  "surface_temp_abovemed", "sex", "min_dist_hf_abovemed")

# make list of modifier subset data  ---------------------------------------
# long
data_hm_long_llist_em <- map(modifier_list, function(x) 
  filter_data(data = data_hm_list$long,
              modifier_name = x))
data_hm_long_list_em <- unlist(data_hm_long_llist_em, recursive = F)
names(data_hm_long_list_em) <-  paste0(rep(modifier_list, each =2), "_", 
                                         names(data_hm_long_list_em), "_",
                                         "long")

# load covariate list ---------------------------------------
covarname <- readRDS(namibia_inc_covar_path)
covarname_g = c("pre_spray_cover", "pre_incidence", "pre_rainfall",
                "pop_size_ea","pre_evi", "ea_elevation","surface_temp")
covarname_Q = covarname[!covarname %in% c(covarname_g)]

# SuperLearner set up ---------------------------------------
source(paste0(here::here(), '/2-analysis/0-namibia-incidence-learners.R'))


# Run hTMLE - Spillover effect ---------------------------------------
res_spillover_hm <- list()

for(i in 1:length(data_hm_long_list_em)){
  print("=====================================================================")
  print(paste0("Human & mosquito - spillover effect #", i, ": ", names(data_hm_long_list_em)[[i]]))
  print("=====================================================================")

  set.seed(123)
  fit <- run_hTMLE(df = data_hm_long_list_em[[i]],
                   effecttype = 'spillover',
                   Qlevel = "individual", # cohort/individual/adaptive for incidence, individual for prev
                   glevel = "cohort", # always cohort for incidence, individual for prev
                   dependency = TRUE, # TRUE if we want CIs adjusted for dependency
                   unadj_est = TRUE, # TRUE if we want a unadjusted estimate
                   y_l = 0,
                   y_u = 1,
                   covarname_g = covarname_g,    # set of names for adjustment variables
                   covarname_Q = covarname_Q,
                   SL_lib_Q = SL_lib_Q,
                   SL_lib_g = SL_lib_g,
                   SL_lib_adaQ = SL_lib_simple, # need to specify if for Qlevel = "adaptive"
                   SL_lib_depen = SL_lib_simple, # need to specify if dependency = T
                   verbose = FALSE,
                   reservoir_name = "Human & mosquito",
                   EM = T,
                   Vfolds = 30,
                   EMname = names(data_human_short_list_em)[i])

  saveRDS(fit, file = paste0(results_path,
                             paste0("namibia_htmle_inc_spillover_hm_",
                                    names(data_hm_long_list_em)[[i]],
                                    ".RDS")))
  rm(fit)
}
