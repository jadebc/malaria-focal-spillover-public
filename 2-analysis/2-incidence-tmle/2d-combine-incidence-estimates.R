################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Combine incidence HTMLE estimates
# for primary and sensitivity analyses
################################################ 
library(purrr)

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# define file root names --------------------------------------------------------
file_list <- c(
  "namibia_htmle_inc_direct_hm_indiv",
               "namibia_htmle_inc_direct_mosq_indiv", 
               "namibia_htmle_inc_direct_human_indiv", "namibia_htmle_inc_spillover_hm_indiv", 
               "namibia_htmle_inc_spillover_mosq_indiv", "namibia_htmle_inc_spillover_human_indiv", 
               "namibia_htmle_inc_total_hm_indiv", "namibia_htmle_inc_total_mosq_indiv", 
               "namibia_htmle_inc_total_human_indiv",
               "namibia_htmle_inc_direct_hm_cohort",
               "namibia_htmle_inc_direct_mosq_cohort", 
               "namibia_htmle_inc_direct_human_cohort", "namibia_htmle_inc_spillover_hm_cohort", 
               "namibia_htmle_inc_spillover_mosq_cohort", "namibia_htmle_inc_spillover_human_cohort", 
               "namibia_htmle_inc_total_hm_cohort", "namibia_htmle_inc_total_mosq_cohort", 
               "namibia_htmle_inc_total_human_cohort")

# primary analysis --------------------------------------------------------
# manually processing combined intervention results
x = readRDS(paste0(results_path, "namibia_htmle_inc_direct_hm_long_indiv.RDS"))
estimates = data.frame(
  Risk1 = NA, Risk0 = NA, Psi_hat = NA, CI_l = NA, CI_u = NA,
  CI_l_unadj = NA, CI_u_unadj = NA, Psi_type = NA,
  Dependency = NA, Qlevel = NA, glevel = NA, reservoir = "Human & mosquito",
  parameter = "direct"
)
res_estX = list(estimates = estimates, SL_coefs = NA)
short = list(res_screen = NA, res_adap_Q = NA, res_est = res_estX)
y = list(Long = x, Short = short)
saveRDS(y, paste0(results_path, "namibia_htmle_inc_direct_hm_long_list_indiv.RDS"))

x = readRDS(paste0(results_path, "namibia_htmle_inc_direct_hm_long_cohort.RDS"))
estimates = data.frame(
  Risk1 = NA, Risk0 = NA, Psi_hat = NA, CI_l = NA, CI_u = NA,
  CI_l_unadj = NA, CI_u_unadj = NA, Psi_type = NA,
  Dependency = NA, Qlevel = NA, glevel = NA, reservoir = "Human & mosquito",
  parameter = "direct"
)
res_estX = list(estimates = estimates, SL_coefs = NA)
short = list(res_screen = NA, res_adap_Q = NA, res_est = res_estX)
y = list(Long = x, Short = short)
saveRDS(y, paste0(results_path, "namibia_htmle_inc_direct_hm_long_list_cohort.RDS"))


files_primary <- paste0(file_list, ".RDS")

files_primary = files_primary[-which(files_primary == "namibia_htmle_inc_direct_hm_indiv.RDS")]
files_primary = files_primary[-which(files_primary == "namibia_htmle_inc_direct_hm_cohort.RDS")]
files_primary = c("namibia_htmle_inc_direct_hm_long_list_indiv.RDS",
                  "namibia_htmle_inc_direct_hm_long_list_cohort.RDS", files_primary)

estimates_primary_list <- lapply(files_primary, read_file_results)
estimates_primary <- bind_rows(estimates_primary_list)

saveRDS(estimates_primary, file = paste0(results_path, "namibia_htmle_inc_all.RDS"))




