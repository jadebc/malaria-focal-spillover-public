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
file_list <- c("namibia_htmle_inc_direct_hm", "namibia_htmle_inc_direct_mosq", 
               "namibia_htmle_inc_direct_human", "namibia_htmle_inc_spillover_hm", 
               "namibia_htmle_inc_spillover_mosq", "namibia_htmle_inc_spillover_human", 
               "namibia_htmle_inc_total_hm", "namibia_htmle_inc_total_mosq", 
               "namibia_htmle_inc_total_human")

# sensitivity analysis: no overlap of target --------------------------------------------------------
files_sens_nooverlap_target <- paste0(file_list, "_sens_nooverlap_target.RDS")
estimates_sens_nooverlap_target_list <- lapply(files_sens_nooverlap_target, read_file_results, both_obs=F)
estimates_sens_nooverlap_target <- bind_rows(estimates_sens_nooverlap_target_list)
saveRDS(estimates_sens_nooverlap_target, file = paste0(results_path, "namibia_htmle_inc_all_sens_nooverlap_target.RDS"))

# sensitivity analysis: no overlap of spill --------------------------------------------------------
# human direct: data too sparse to fit model
file_list_no_spill <- file_list[-which("namibia_htmle_inc_direct_human" == file_list)]
files_sens_nooverlap_spill <- paste0(file_list_no_spill, "_sens_nooverlap_spill.RDS")
estimates_sens_nooverlap_spill_list <- lapply(files_sens_nooverlap_spill, read_file_results, both_obs=F)
estimates_sens_nooverlap_spill <- bind_rows(estimates_sens_nooverlap_spill_list)
saveRDS(estimates_sens_nooverlap_spill, file = paste0(results_path, "namibia_htmle_inc_all_sens_nooverlap_spill.RDS"))

# sensitivity analysis: 2km spillover zone --------------------------------------------------------
files_sens_spzone_2km <- paste0(file_list, "_sens_spzone_2km.RDS")
# drop direct effects
files_sens_spzone_2km <- files_sens_spzone_2km[-grep("direct", files_sens_spzone_2km)]
estimates_sens_spzone_2km_list <- lapply(files_sens_spzone_2km, read_file_results, both_obs=F)
estimates_sens_spzone_2km <- bind_rows(estimates_sens_spzone_2km_list)
saveRDS(estimates_sens_spzone_2km, file = paste0(results_path, "namibia_htmle_inc_all_sens_spzone_2km.RDS"))

# sensitivity analysis: 3km spillover zone --------------------------------------------------------
files_sens_spzone_3km <- paste0(file_list, "_sens_spzone_3km.RDS")
# drop direct effects
files_sens_spzone_3km <- files_sens_spzone_3km[-grep("direct", files_sens_spzone_3km)]
estimates_sens_spzone_3km_list <- lapply(files_sens_spzone_3km, read_file_results, both_obs=F)
estimates_sens_spzone_3km <- bind_rows(estimates_sens_spzone_3km_list)
saveRDS(estimates_sens_spzone_3km, file = paste0(results_path, "namibia_htmle_inc_all_sens_spzone_3km.RDS"))

# sensitivity analysis: different observation period --------------------------------------------------------
# mosq and hm direct: data too sparse to fit model
file_list_obs <- file_list[-which("namibia_htmle_inc_direct_human" == file_list)]
file_list_obs <- file_list_obs[-which("namibia_htmle_inc_direct_mosq" == file_list_obs)]
file_list_obs <- file_list_obs[-which("namibia_htmle_inc_direct_hm" == file_list_obs)]
files_sens_obs <- paste0(file_list_obs, "_sens_obs.RDS")
estimates_sens_obs_list <- lapply(files_sens_obs, read_file_results, both_obs=F)
estimates_sens_obs <- bind_rows(estimates_sens_obs_list)
saveRDS(estimates_sens_obs, file = paste0(results_path, "namibia_htmle_inc_all_sens_obs.RDS"))



