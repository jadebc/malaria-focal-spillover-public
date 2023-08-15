################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Prevalence analysis
################################################  
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))
library(sl3)
library(tmle3)
library(origami)
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(GGally)


# load and process analysis dataset ---------------------------------------------------
xs_clean <- readRDS(namibia_analysis_prev)

# specify the variables used for identification, but not for estimation
dropme <- c("eaid", "hhid", "xsid", "iid", "latitude", "longitude",
            "date", "spillover_zone")

xs <- xs_clean[,!colnames(xs_clean) %in% dropme] 

# clean RDT results
xs <- xs %>% mutate(
  rdtyn = case_when(
    rdtresult == "Negative" ~ 0,
    rdtresult == "Pf only" ~ 1,
    rdtresult == "Pf or Mixed" ~ 1, 
    rdtresult == "Invalid" ~ NA_real_, 
    rdtresult == "Not conducted" ~ NA_real_
  ), 
  
  hsrdtyn = case_when(
    hsrdtresult == "Negative" ~ 0,
    hsrdtresult == "Positive" ~ 1,
    hsrdtresult == "Invalid" ~ NA_real_, 
    hsrdtresult == "Not conducted" ~ NA_real_
  )
)



# run tmle function ---------------------------------------------------
outcome_list <- list("qPCRposneg", "rdtyn", "hsrdtyn")

## human reservoir ---- 
set.seed(123)
res_human_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Direct effect", reservoir = "human", yname= x))
res_human_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Spillover effect", reservoir = "human", yname= x)) 
res_human_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Total effect", reservoir = "human", yname= x)) 

names(res_human_de) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_human_sp) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_human_te) = c("qPCRposneg", "rdtyn", "hsrdtyn")

## mosquito reservoir ---- 
set.seed(123)
res_mosq_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Direct effect", reservoir = "mosquito", yname= x))
res_mosq_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Spillover effect", reservoir = "mosquito", yname= x)) 
res_mosq_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Total effect", reservoir = "mosquito", yname= x))

names(res_mosq_de) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_mosq_sp) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_mosq_te) = c("qPCRposneg", "rdtyn", "hsrdtyn")

## human and mosquito reservoir ---- 
set.seed(123)
res_both_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Direct effect", reservoir = "human & mosquito", yname= x)) 
res_both_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Spillover effect", reservoir = "human & mosquito", yname= x)) 
res_both_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = xs, parameter = "Total effect", reservoir = "human & mosquito", yname= x))

names(res_both_de) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_both_sp) = c("qPCRposneg", "rdtyn", "hsrdtyn")
names(res_both_te) = c("qPCRposneg", "rdtyn", "hsrdtyn")

# save results ---------------------------------------------------
prev_tmle_results <- list(
  res_human_de = res_human_de, 
  res_human_sp = res_human_sp,
  res_human_te = res_human_te,
  
  res_mosq_de = res_mosq_de,
  res_mosq_sp = res_mosq_sp,
  res_mosq_te = res_mosq_te,
  
  res_both_de = res_both_de,
  res_both_sp = res_both_sp,
  res_both_te = res_both_te
  
)

saveRDS(prev_tmle_results, file = paste0(results_path, "prevalence-tmle-results.RDS"))
