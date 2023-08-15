################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Prevalence analysis
# Serological outcomes
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
library(readxl)

# load and process analysis dataset ---------------------------------------------------
xs_clean <- readRDS(namibia_analysis_prev) %>% arrange(eaid, hhid, iid)
raw_data <- read_excel(namibia_raw_sero_path)

sero_cols <- raw_data %>% dplyr::select(eaid, hhid, sample, Etramp5.Ag1, Etramp5.Ag1_pos) %>% 
  rename(iid = sample) %>% 
  arrange(eaid, hhid, iid) %>% 
  distinct()

sero <- full_join(xs_clean, sero_cols, by = c("eaid", "hhid", "iid"))

# drop if no spillover zone information
sero <- sero %>% filter(!is.na(spall))

# drop if no serology
sero <- sero %>% filter(!is.na(Etramp5.Ag1))


# specify the variables used for identification, but not for estimation
dropme <- c("eaid", "hhid", "xsid", "iid", "latitude", "longitude",
            "date", "spillover_zone")

sero <- sero[,!colnames(sero) %in% dropme] 

# run tmle function ---------------------------------------------------
reservoir_list = c("human", "mosquito", "human & mosquito")

## human reservoir ---- 
set.seed(123)
res_de <- lapply(reservoir_list, function(x) fit_tmle_prev(
  data = sero, parameter = "Direct effect", reservoir = x, yname= "Etramp5.Ag1_pos")) 
res_sp <- lapply(reservoir_list, function(x) fit_tmle_prev(
  data = sero, parameter = "Spillover effect", reservoir = x, yname= "Etramp5.Ag1_pos")) 
res_te <- lapply(reservoir_list, function(x) fit_tmle_prev(
  data = sero, parameter = "Total effect", reservoir = x, yname= "Etramp5.Ag1_pos")) 

names(res_de) = c("human", "mosquito", "human & mosquito")
names(res_sp) = c("human", "mosquito", "human & mosquito")
names(res_te) = c("human", "mosquito", "human & mosquito")

# save results ---------------------------------------------------
prev_tmle_results <- list(
  res_de = res_de, 
  res_sp = res_sp,
  res_te = res_te
)

saveRDS(prev_tmle_results, file = paste0(results_path, "prevalence-tmle-results-sero.RDS"))

