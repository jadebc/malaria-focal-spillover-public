################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Configure libraries and paths
# Source base functions
################################################

# Load libraries ----------------------------------------------------------
library(ggplot2)
library(readstata13)
library(rgdal)
library(assertthat)
library(tidyr)
library(geosphere)
library(RColorBrewer)
library(lubridate)
library(reshape2)
library(purrr)
library(tictoc)
library(data.table)
library(stringr)
library(dplyr)
library(GGally)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# for survival curves
library(ggfortify)
library(survival)

# Source base functions ----------------------------------------------------------
file_sources = list.files(path = paste0(here::here(), "/0-base-functions"), pattern="*.R")
file_names = as.list(paste0(here::here(),"/0-base-functions/",file_sources))
lapply(file_names, source, .GlobalEnv)

# Detect system, identify root directory  ----------------------------------------------------------
if(Sys.getenv("LOGNAME")=="jadebc"){
  box_path = "/Users/jadebc/Library/CloudStorage/Box-Box/Jade Benjamin-Chung's Externally Shareable Files/"
}
# else{
#   box_path = "~/Box Sync/"
# }

if(Sys.getenv("LOGNAME")=="jadebc"){
  box_shared_path = "/Users/jadebc/Library/CloudStorage/Box-Box/NAMEP SMT Shared Folder/"
}
# # else{
#   box_shared_path = "~/Box Sync/"
# }

if(Sys.getenv("LOGNAME")=="gabriellabh"){
  box_shared_path = "/Users/gabriellabh/Library/CloudStorage/Box-Box/"
  box_path = "/Users/gabriellabh/Library/CloudStorage/Box-Box/"
  namibia_raw_path = paste0(box_path, "K01-data/")
}
# else{
#   box_shared_path = "~/Box Sync/"
# }

if(Sys.getenv("LMOD_SYSHOST")=="sherlock"){
  sherlock_path = paste0(Sys.getenv("HOME"), "/K01/")
}

# Define raw data paths  ----------------------------------------------------------
# Namibia paths 
namibia_raw_path = paste0(box_shared_path, "Datasets/Final trial data for collaborators/")

namibia_raw_gr_path = paste0(namibia_raw_path, "GR/GR_2017_updates_EA_joined.csv")
namibia_raw_xs_path = paste0(box_shared_path, "Datasets/2017 x-sectional/XS_v2.dta")
namibia_raw_sero_path = paste0(box_shared_path, "Datasets/2017 x-sectional/namibia_sero_forhotspot.xlsx")
namibia_raw_index_path = paste0(namibia_raw_path, "Index/SDSS_index_deidentified.dta")
namibia_raw_intervention_path = paste0(namibia_raw_path, "Individual/Individual_deidentified.dta")
namibia_raw_ea_path = paste0(namibia_raw_path, "EA level/EA_level_dataset.dta")
namibia_raw_ravc_path = paste0(namibia_raw_path, "RAVC/RAVC_data_deidentified.dta")

namibia_shapefile_dsn = paste0(namibia_raw_path,"Shapefiles/Study Area/")
namibia_shapefile_layer = "EA_study_area_Y2"



# Eswatini paths
eswatini_raw_path = paste0(box_shared_path, "Datasets/Final data for collaborators/")

eswatini_raw_gr_path = paste0(eswatini_raw_path, "Household/allhouse_eajoined.dta")
eswatini_raw_index_path = paste0(eswatini_raw_path, "Index/swaziland_indexes.dta")
eswatini_raw_intervention_path = paste0(eswatini_raw_path, "Individual/individual_clean.dta")

eswatini_new_shapefile_dsn = paste0(eswatini_raw_path, 'Shapefiles/EA level- "reprojected"/')
eswatini_new_shapefile_layer = "reprojected shape file "

eswatini_old_shapefile_dsn = paste0(eswatini_raw_path, 'Shapefiles/')
eswatini_old_shapefile_layer = "LocalitiesSWA_new2017_verified_arm"

eswatini_raw_ea_level_path = paste0(eswatini_raw_path, "Incidence & IRRs/locality_incidence_cases.csv")

eswatini_ea_shp_path = "~/Box/Swaziland TPE Study/Shape Files/LocalitiesProjectComplete-LocCode (1)/LocalitiesProjectComplete-LocCode.shp"

eswatini_raw_ea_path = "~/Box/Swaziland TPE Study/Data sharing/Final data for collaborators/Incidence & IRRs/locality_incidence_cases.dta"

# Define clean data paths  ----------------------------------------------------------
## Namibia paths   ----------------------------------------------------------
namibia_clean_path = paste0(box_path, "K01-Data/Data/Namibia/clean-data/")
namibia_matched_path = paste0(box_path, "K01-Data/Data/Namibia/matched-data/")
namibia_cov_path = paste0(box_path, "K01-Data/Data/Namibia/spatial-covariates/")

namibia_clean_gr_path = paste0(namibia_clean_path, "namibia_gr_clean.RDS")

namibia_clean_index_intervention_path = paste0(namibia_clean_path, "namibia_index_intervention_clean.RDS")
namibia_clean_index_nonintervention_path = paste0(namibia_clean_path, "namibia_index_nonintervention_clean.RDS")
namibia_clean_xs_path = paste0(namibia_clean_path, "namibia_xs_clean.RDS")

namibia_partial_clean_intervention_path = paste0(namibia_clean_path, "namibia_intervention_partial_clean.RDS")

namibia_int_ind_data_path = paste0(namibia_clean_path, "namibia_merged_intervention_index.RDS")
namibia_intervention_case_counts_path = paste0(namibia_clean_path, "namibia_intervention_case_counts.RDS")
namibia_analysis_data_path = paste0(namibia_clean_path, "namibia_analysis_dataset.RDS")

namibia_ea_data_path = paste0(namibia_clean_path, "namibia_ea_data.RDS")
namibia_ea_day_data_path = paste0(namibia_clean_path, "namibia_ea_day_data.RDS")



### sherlock paths   ----------------------------------------------------------
if(Sys.getenv("LMOD_SYSHOST")=="sherlock"){
  sherlock_clean_path = paste0(sherlock_path, "clean-data/")
  sherlock_matched_path = paste0(sherlock_path, "matched-data/")
  sherlock_cov_path = paste0(sherlock_path, "spatial-covariates/")
  
  namibia_raw_ea_level_path = paste0(sherlock_path, "data/EA_level_dataset.csv")
  
  namibia_indiv_imputed_path = paste0(sherlock_clean_path, "namibia_indiv_imputed.RDS")
  namibia_clean_index_path = paste0(sherlock_clean_path, "namibia_index_clean.RDS")
  namibia_clean_intervention_path = paste0(sherlock_clean_path, "namibia_intervention_clean.RDS")
  
  namibia_index_matched_path = paste0(sherlock_matched_path, "index_matched_df.RDS")
  namibia_pt_path = paste0(sherlock_matched_path, "person_time.RDS")
  namibia_temp_int_match_path = paste0(sherlock_matched_path, "intervention_matched_df.RDS")
  namibia_temp_int_unmatch_path = paste0(sherlock_matched_path, "intervention_unmatched_df.RDS")
  namibia_int_match_path = paste0(sherlock_matched_path, "intervention_matched.RDS")
  namibia_temp_nontrig_match_path = paste0(sherlock_matched_path, "intervention_nonint_matched_df.RDS")
  namibia_temp_nontrig_unmatch_path = paste0(sherlock_matched_path, "intervention_nonint_unmatched_df.RDS")
  
  # spatial covariates
  count_int_pre_recei_path = paste0(sherlock_cov_path, "df_count_int_pre_recei.RDS")
  count_int_pre_deliv_path = paste0(sherlock_cov_path, "df_count_int_pre_deliv.RDS")
  count_int_concu_deliv_path = paste0(sherlock_cov_path, "df_count_int_concu_deliv.RDS")
  count_pre_index_path = paste0(sherlock_cov_path, "df_count_pre_index.RDS")
  count_pop_path = paste0(sherlock_cov_path, "df_count_pop.RDS")
  count_int_concu_deliv_s_path = paste0(sherlock_cov_path, "df_count_int_concu_deliv_s.RDS")
  count_int_concu_deliv_l_path = paste0(sherlock_cov_path, "df_count_int_concu_deliv_l.RDS")
  
  # primary analysis 
  namibia_cohort_path = paste0(sherlock_clean_path, "namibia_cohorts.RDS")
  namibia_sub_coho_index_int_path = paste0(sherlock_clean_path, "namibia_sub_coho_index_int.RDS")
  namibia_sub_coho_index_nonint_path = paste0(sherlock_clean_path, "namibia_sub_coho_index_nonint.RDS")
  namibia_sub_coho_int_recip_path = paste0(sherlock_clean_path, "namibia_sub_coho_int_recip.RDS")
  namibia_coho_person_list_path = paste0(sherlock_clean_path, "namibia_coho_person_list.RDS")
  
  namibia_df_base_path = paste0(sherlock_clean_path, "namibia_df_base.RDS")
  namibia_df_base_path_sens_spzone_2km = paste0(sherlock_clean_path, "namibia_df_base_sens_spzone_2km.RDS")
  namibia_df_base_path_sens_spzone_3km = paste0(sherlock_clean_path, "namibia_df_base_sens_spzone_3km.RDS")
  
  namibia_df_all_path = paste0(sherlock_clean_path, "namibia_df_all.RDS")
  namibia_df_long_path = paste0(sherlock_clean_path, "namibia_df_long.RDS")
  namibia_df_short_path = paste0(sherlock_clean_path, "namibia_df_short.RDS")
  namibia_df_short_permute_path = paste0(sherlock_clean_path, "namibia_df_short_permute.RDS")
  namibia_df_long_permute_path = paste0(sherlock_clean_path, "namibia_df_long_permute.RDS")
  namibia_df_short_deid_path = paste0(sherlock_clean_path, "namibia_df_short_deid.RDS")
  namibia_df_long_deid_path = paste0(sherlock_clean_path, "namibia_df_long_deid.RDS")
  
  # paths for sensitivity analyses
  namibia_df_short_sens_obs_path = paste0(sherlock_clean_path, "namibia_df_sens_obs_short.RDS")
  namibia_df_long_sens_obs_path = paste0(sherlock_clean_path, "namibia_df_sens_obs_long.RDS")
  namibia_df_short_sens_spzone_2km_path = paste0(sherlock_clean_path, "namibia_df_sens_spzone_2km_short.RDS")
  namibia_df_short_sens_spzone_3km_path = paste0(sherlock_clean_path, "namibia_df_sens_spzone_3km_short.RDS")
  namibia_df_long_sens_spzone_2km_path = paste0(sherlock_clean_path, "namibia_df_sens_spzone_2km_long.RDS")
  namibia_df_long_sens_spzone_3km_path = paste0(sherlock_clean_path, "namibia_df_sens_spzone_3km_long.RDS")
  
  namibia_df_short_sens_nooverlap_target_path = paste0(sherlock_clean_path, "namibia_df_sens_nooverlap_target_short.RDS")
  namibia_df_short_sens_nooverlap_spill_path = paste0(sherlock_clean_path, "namibia_df_sens_nooverlap_spill_short.RDS")
  namibia_df_long_sens_nooverlap_target_path = paste0(sherlock_clean_path, "namibia_df_sens_nooverlap_target_long.RDS")
  namibia_df_long_sens_nooverlap_spill_path = paste0(sherlock_clean_path, "namibia_df_sens_nooverlap_spill_long.RDS")
  
  
  ### local paths   ----------------------------------------------------------
}else{
  
  namibia_raw_ea_level_path = paste0(namibia_raw_path, "EA level/EA_level_dataset.csv")
  
  namibia_indiv_imputed_path = paste0(namibia_clean_path, "namibia_indiv_imputed.RDS")
  namibia_clean_index_path = paste0(namibia_clean_path, "namibia_index_clean.RDS")
  namibia_clean_intervention_path = paste0(namibia_clean_path, "namibia_intervention_clean.RDS")
  
  namibia_index_matched_path = paste0(namibia_matched_path, "index_matched_df.RDS")
  namibia_pt_path = paste0(namibia_matched_path, "person_time.RDS")
  namibia_temp_int_match_path = paste0(namibia_matched_path, "intervention_matched_df.RDS")
  namibia_temp_int_unmatch_path = paste0(namibia_matched_path, "intervention_unmatched_df.RDS")
  namibia_int_match_path = paste0(namibia_matched_path, "intervention_matched.RDS")
  namibia_temp_nontrig_match_path = paste0(namibia_matched_path, "intervention_nonint_matched_df.RDS")
  namibia_temp_nontrig_unmatch_path = paste0(namibia_matched_path, "intervention_nonint_unmatched_df.RDS")
  
  # spatial covariates
  count_int_pre_recei_path = paste0(namibia_cov_path, "df_count_int_pre_recei.RDS")
  count_int_pre_deliv_path = paste0(namibia_cov_path, "df_count_int_pre_deliv.RDS")
  count_int_concu_deliv_path = paste0(namibia_cov_path, "df_count_int_concu_deliv.RDS")
  count_int_concu_deliv_s_path = paste0(namibia_cov_path, "df_count_int_concu_deliv_s.RDS")
  count_int_concu_deliv_l_path = paste0(namibia_cov_path, "df_count_int_concu_deliv_l.RDS")
  count_pre_index_path = paste0(namibia_cov_path, "df_count_pre_index.RDS")
  count_pop_path = paste0(namibia_cov_path, "df_count_pop.RDS")
  
  # primary analysis 
  namibia_cohort_path = paste0(namibia_clean_path, "namibia_cohorts.RDS")
  namibia_sub_coho_index_int_path = paste0(namibia_clean_path, "namibia_sub_coho_index_int.RDS")
  namibia_sub_coho_index_nonint_path = paste0(namibia_clean_path, "namibia_sub_coho_index_nonint.RDS")
  namibia_sub_coho_int_recip_path = paste0(namibia_clean_path, "namibia_sub_coho_int_recip.RDS")
  namibia_coho_person_list_path = paste0(namibia_clean_path, "namibia_coho_person_list.RDS")
  
  namibia_df_base_path = paste0(namibia_clean_path, "namibia_df_base.RDS")
  namibia_df_base_path_sens_spzone_2km = paste0(namibia_clean_path, "namibia_df_base_sens_spzone_2km.RDS")
  namibia_df_base_path_sens_spzone_3km = paste0(namibia_clean_path, "namibia_df_base_sens_spzone_3km.RDS")
  
  namibia_df_all_path = paste0(namibia_clean_path, "namibia_df_all.RDS")
  namibia_df_short_path = paste0(namibia_clean_path, "namibia_df_short.RDS")
  namibia_df_long_path = paste0(namibia_clean_path, "namibia_df_long.RDS")
  namibia_df_short_permute_path = paste0(namibia_clean_path, "namibia_df_short_permute.RDS")
  namibia_df_long_permute_path = paste0(namibia_clean_path, "namibia_df_long_permute.RDS")
  namibia_df_short_deid_path = paste0(namibia_clean_path, "namibia_df_short_deid.RDS")
  namibia_df_long_deid_path = paste0(namibia_clean_path, "namibia_df_long_deid.RDS")
  
  # paths for sensitivity analyses
  namibia_df_short_sens_obs_path = paste0(namibia_clean_path, "namibia_df_sens_obs_short.RDS")
  namibia_df_long_sens_obs_path = paste0(namibia_clean_path, "namibia_df_sens_obs_long.RDS")
  namibia_df_short_sens_spzone_2km_path = paste0(namibia_clean_path, "namibia_df_sens_spzone_2km_short.RDS")
  namibia_df_short_sens_spzone_3km_path = paste0(namibia_clean_path, "namibia_df_sens_spzone_3km_short.RDS")
  namibia_df_long_sens_spzone_2km_path = paste0(namibia_clean_path, "namibia_df_sens_spzone_2km_long.RDS")
  namibia_df_long_sens_spzone_3km_path = paste0(namibia_clean_path, "namibia_df_sens_spzone_3km_long.RDS")
  
  namibia_df_short_sens_nooverlap_target_path = paste0(namibia_clean_path, "namibia_df_sens_nooverlap_target_short.RDS")
  namibia_df_short_sens_nooverlap_spill_path = paste0(namibia_clean_path, "namibia_df_sens_nooverlap_spill_short.RDS")
  namibia_df_long_sens_nooverlap_target_path = paste0(namibia_clean_path, "namibia_df_sens_nooverlap_target_long.RDS")
  namibia_df_long_sens_nooverlap_spill_path = paste0(namibia_clean_path, "namibia_df_sens_nooverlap_spill_long.RDS")
  
}


## Eswatini paths   ----------------------------------------------------------
eswatini_clean_path = paste0(box_path, "K01-Data/Data/Eswatini/")

eswatini_clean_gr_path = paste0(eswatini_clean_path, "eswatini_gr_clean.RDS")
eswatini_indiv_imputed_path = paste0(eswatini_clean_path, "eswatini_indiv_imputed.RDS")

eswatini_partial_clean_intervention_path = paste0(eswatini_clean_path, "eswatini_intervention_partial_clean.RDS")
eswatini_clean_intervention_path = paste0(eswatini_clean_path, "eswatini_intervention_clean.RDS")

eswatini_clean_index_path = paste0(eswatini_clean_path, "eswatini_index_clean.RDS")

eswatini_cohort_path = paste0(eswatini_clean_path, "eswatini_cohorts.RDS")


# Define processed data paths  ----------------------------------------------------------

### sherlock paths   ----------------------------------------------------------
if(Sys.getenv("LMOD_SYSHOST")=="sherlock"){
  sherlock_process_path = paste0(sherlock_path, "processed-data/")
  
  namibia_human_process_path = paste0(sherlock_process_path, "namibia_human_process.RDS")
  namibia_mosq_process_path = paste0(sherlock_process_path, "namibia_mosq_process.RDS")
  namibia_hm_process_path = paste0(sherlock_process_path, "namibia_hm_process.RDS")
  
  namibia_human_process_sens_obs_path = paste0(sherlock_process_path, "namibia_human_process_sens_obs.RDS")
  namibia_mosq_process_sens_obs_path = paste0(sherlock_process_path, "namibia_mosq_process_sens_obs.RDS")
  namibia_hm_process_sens_obs_path = paste0(sherlock_process_path, "namibia_hm_process_sens_obs.RDS")
  
  namibia_human_process_sens_nooverlap_target_path = paste0(sherlock_process_path, "namibia_human_process_sens_nooverlap_target.RDS")
  namibia_mosq_process_sens_nooverlap_target_path = paste0(sherlock_process_path, "namibia_mosq_process_sens_nooverlap_target.RDS")
  namibia_hm_process_sens_nooverlap_target_path = paste0(sherlock_process_path, "namibia_hm_process_sens_nooverlap_target.RDS")
  
  namibia_human_process_sens_nooverlap_spill_path = paste0(sherlock_process_path, "namibia_human_process_sens_nooverlap_spill.RDS")
  namibia_mosq_process_sens_nooverlap_spill_path = paste0(sherlock_process_path, "namibia_mosq_process_sens_nooverlap_spill.RDS")
  namibia_hm_process_sens_nooverlap_spill_path = paste0(sherlock_process_path, "namibia_hm_process_sens_nooverlap_spill.RDS")
  
  namibia_human_process_sens_spzone_2km_path = paste0(sherlock_process_path, "namibia_human_process_sens_spzone_2km.RDS")
  namibia_mosq_process_sens_spzone_2km_path = paste0(sherlock_process_path, "namibia_mosq_process_sens_spzone_2km.RDS")
  namibia_hm_process_sens_spzone_2km_path = paste0(sherlock_process_path, "namibia_hm_process_sens_spzone_2km.RDS")
  
  namibia_human_process_sens_spzone_3km_path = paste0(sherlock_process_path, "namibia_human_process_sens_spzone_3km.RDS")
  namibia_mosq_process_sens_spzone_3km_path = paste0(sherlock_process_path, "namibia_mosq_process_sens_spzone_3km.RDS")
  namibia_hm_process_sens_spzone_3km_path = paste0(sherlock_process_path, "namibia_hm_process_sens_spzone_3km.RDS")
  
  # covariate name paths
  sherlock_covar_path = paste0(sherlock_path, "covariates/")
  
  namibia_inc_covar_path = paste0(sherlock_covar_path, "namibia_inc_covar.RDS")
  namibia_inc_covar_sens_obs_path = paste0(sherlock_covar_path, "namibia_inc_covar_sens_obs.RDS")
  namibia_inc_covar_sens_nooverlap_target_path = paste0(sherlock_covar_path, "namibia_inc_covar_sens_nooverlap_target.RDS")
  namibia_inc_covar_sens_nooverlap_spill_path = paste0(sherlock_covar_path, "namibia_inc_covar_sens_nooverlap_spill.RDS")
  namibia_inc_covar_sens_spzone_2km_path = paste0(sherlock_covar_path, "namibia_inc_covar_sens_spzone_2km.RDS")
  namibia_inc_covar_sens_spzone_3km_path = paste0(sherlock_covar_path, "namibia_inc_covar_sens_spzone_3km.RDS")
  
  
}else{
  ### local paths   ----------------------------------------------------------
  namibia_process_path = paste0(box_path, "K01-Data/Data/Namibia/processed-data/")
  
  namibia_human_process_path = paste0(namibia_process_path, "namibia_human_process.RDS")
  namibia_mosq_process_path = paste0(namibia_process_path, "namibia_mosq_process.RDS")
  namibia_hm_process_path = paste0(namibia_process_path, "namibia_hm_process.RDS")
  
  namibia_human_process_sens_obs_path = paste0(namibia_process_path, "namibia_human_process_sens_obs.RDS")
  namibia_mosq_process_sens_obs_path = paste0(namibia_process_path, "namibia_mosq_process_sens_obs.RDS")
  namibia_hm_process_sens_obs_path = paste0(namibia_process_path, "namibia_hm_process_sens_obs.RDS")
  
  namibia_human_process_sens_nooverlap_target_path = paste0(namibia_process_path, "namibia_human_process_sens_nooverlap_target.RDS")
  namibia_mosq_process_sens_nooverlap_target_path = paste0(namibia_process_path, "namibia_mosq_process_sens_nooverlap_target.RDS")
  namibia_hm_process_sens_nooverlap_target_path = paste0(namibia_process_path, "namibia_hm_process_sens_nooverlap_target.RDS")
  
  namibia_human_process_sens_nooverlap_spill_path = paste0(namibia_process_path, "namibia_human_process_sens_nooverlap_spill.RDS")
  namibia_mosq_process_sens_nooverlap_spill_path = paste0(namibia_process_path, "namibia_mosq_process_sens_nooverlap_spill.RDS")
  namibia_hm_process_sens_nooverlap_spill_path = paste0(namibia_process_path, "namibia_hm_process_sens_nooverlap_spill.RDS")
  
  namibia_human_process_sens_spzone_2km_path = paste0(namibia_process_path, "namibia_human_process_sens_spzone_2km.RDS")
  namibia_mosq_process_sens_spzone_2km_path = paste0(namibia_process_path, "namibia_mosq_process_sens_spzone_2km.RDS")
  namibia_hm_process_sens_spzone_2km_path = paste0(namibia_process_path, "namibia_hm_process_sens_spzone_2km.RDS")
  
  namibia_human_process_sens_spzone_3km_path = paste0(namibia_process_path, "namibia_human_process_sens_spzone_3km.RDS")
  namibia_mosq_process_sens_spzone_3km_path = paste0(namibia_process_path, "namibia_mosq_process_sens_spzone_3km.RDS")
  namibia_hm_process_sens_spzone_3km_path = paste0(namibia_process_path, "namibia_hm_process_sens_spzone_3km.RDS")
  
  namibia_analysis_prev = paste0(namibia_process_path, "namibia_prev_analysis.RDS")
  
  # covariate name paths
  namibia_covar_path = paste0(box_path, "K01-Data/Data/Namibia/covariates/")
  
  namibia_inc_covar_path = paste0(namibia_covar_path, "namibia_inc_covar.RDS")
  namibia_inc_covar_sens_obs_path = paste0(namibia_covar_path, "namibia_inc_covar_sens_obs.RDS")
  namibia_inc_covar_sens_nooverlap_target_path = paste0(namibia_covar_path, "namibia_inc_covar_sens_nooverlap_target.RDS")
  namibia_inc_covar_sens_nooverlap_spill_path = paste0(namibia_covar_path, "namibia_inc_covar_sens_nooverlap_spill.RDS")
  namibia_inc_covar_sens_spzone_2km_path = paste0(namibia_covar_path, "namibia_inc_covar_sens_spzone_2km.RDS")
  namibia_inc_covar_sens_spzone_3km_path = paste0(namibia_covar_path, "namibia_inc_covar_sens_spzone_3km.RDS")
}



# Define results paths   ----------------------------------------------------------
# sherlock paths
if(Sys.getenv("LMOD_SYSHOST")=="sherlock"){
  results_path = paste0(sherlock_path, "results/")
  figure_path = paste0(sherlock_path, "/figures/")
  table_path = paste0(sherlock_path, "/tables/")
}else{
  results_path = paste0(box_path, "K01-Data/Results/Namibia/")
  figure_path = paste0(here::here(), "/6-figures/")
  table_path = paste0(here::here(), "/7-tables/")
}




