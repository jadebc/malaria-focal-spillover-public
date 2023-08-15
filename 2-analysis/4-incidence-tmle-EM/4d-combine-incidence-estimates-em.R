################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Combine incidence HTMLE estimates
# for effect modification analyses
################################################ 
library(purrr)

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# define file root names --------------------------------------------------------
all_files <- list.files(path = results_path)
modifier_files <- all_files[grep('modifier', all_files)]

# function to read and process files --------------------------------------------------------
read_file_results <- function(file_path){
  results <- readRDS(paste0(results_path, file_path))
  
  # if model wasn't fit, return data frame with missings
  if(length(results$res_est)==1){
    print(paste0("File '", file_path, "' contains no results."))
    
    if(length(grep("_human_", file_path)==1)) reservoir = "Human"
    if(length(grep("_mosq_", file_path)==1)) reservoir = "Mosquito"
    if(length(grep("_hm_", file_path)==1)) reservoir = "Human & mosquito"

    if(length(grep("spillover", file_path)==1)) parameter = "spillover"
    if(length(grep("total", file_path)==1)) parameter = "total"
    if(length(grep("direct", file_path)==1)) parameter = "direct"
    
    if(length(grep("short", file_path)==1)) period = "5 weeks"
    if(length(grep("long", file_path)==1)) period = "6 months"
    
    if(reservoir=="Human") start =  str_locate(file_path, "human_")[2] + 1
    if(reservoir=="Mosquito") start =  str_locate(file_path, "mosq_")[2] + 1
    if(reservoir=="Human & mosquito") start = str_locate(file_path, "hm_")[2] + 1
    end = str_locate(file_path, "_abovemed")[1] - 1
    modifier = substr(file_path, start, end)
      
    if(length(grep("modifier0", file_path)==1)) modifier_level = 0
    if(length(grep("modifier1", file_path)==1)) modifier_level = 1
    
    out = data.frame(
      Risk1 = NA, 
      Risk0 = NA,
      Psi_hat = NA,
      CI_l = NA,
      CI_u = NA,
      CI_l_unadj = NA,
      CI_u_unadj = NA, 
      Psi_type = NA,
      Dependency = NA,
      Qlevel = NA,
      glevel = NA, 
      reservoir = reservoir,
      parameter = parameter,
      period = period,
      modifier = modifier,
      modifier_level = modifier_level 
    )
    
    return(out)
  } 
  
  # if model was fit, format results
  if(length(results$res_est)>1){
    estimates <- results$res_est$estimates
    
    if(!is.null(estimates)){
      # identify observation period
      if(length(grep("long", file_path))==1) estimates = estimates %>% mutate(period = "6 months")
      if(length(grep("short", file_path))==1) estimates = estimates %>% mutate(period = "5 weeks")
      
      # define modifier list  ---------------------------------------
      modifier_list <- list("pre_incidence_abovemed", "tx_cov_cohort_abovemed", "pre_spray_cover_abovemed",
                            "pre_rainfall_abovemed", "pre_evi_abovemed", "ea_elevation_abovemed", 
                            "surface_temp_abovemed" ,"sex", "min_dist_hf_abovemed")
      
      # get modifier name
      is_modifier <- lapply(modifier_list, function(x) grep(x, file_path))
      estimates <- estimates %>% mutate(modifier= modifier_list[which(is_modifier==1)] %>% unlist())
      
      # get modifier level
      if(length(grep("modifier1", file_path))==1) estimates = estimates %>% mutate(modifier_level = 1)
      if(length(grep("modifier0", file_path))==1) estimates = estimates %>% mutate(modifier_level = 0)
      
      return(estimates)
    }
  } 

}


# sensitivity analysis: different outcome bounds --------------------------------------------------------
estimates_modifier_list <- lapply(modifier_files, read_file_results)
estimates_modifier <- bind_rows(estimates_modifier_list)
saveRDS(estimates_modifier, file = paste0(results_path, "namibia_htmle_inc_all_em.RDS"))


