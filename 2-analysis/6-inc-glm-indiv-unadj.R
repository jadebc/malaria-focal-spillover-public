################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence
# GLM model - individual-level analysis
# Unadjusted
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
data_human_list <- readRDS(namibia_human_process_path)
data_mosq_list <- readRDS(namibia_mosq_process_path)
data_hm_list <- readRDS(namibia_hm_process_path)




# human reservoir ---------------------------------------
res_total_human <- lapply(data_human_list, function(x)
  fit_indiv_glm_unadj(data = x,
                      effecttype = "total"))

res_spillover_human <- lapply(data_human_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "spillover"))

res_direct_human <- lapply(data_human_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "direct"))

# save results
saveRDS(res_total_human, file = paste0(results_path, "namibia_glm_inc_total_human_unadj.RDS"))
saveRDS(res_spillover_human, file = paste0(results_path, "namibia_glm_inc_spillover_human_unadj.RDS"))
saveRDS(res_direct_human, file = paste0(results_path, "namibia_glm_inc_direct_human_unadj.RDS"))


# mosquito reservoir ---------------------------------------
res_total_mosq <- lapply(data_mosq_list, function(x)
  fit_indiv_glm_unadj(data = x,
                      effecttype = "total"))

res_spillover_mosq <- lapply(data_mosq_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "spillover"))

res_direct_mosq <- lapply(data_mosq_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "direct"))

# save results
saveRDS(res_total_mosq, file = paste0(results_path, "namibia_glm_inc_total_mosq_unadj.RDS"))
saveRDS(res_spillover_mosq, file = paste0(results_path, "namibia_glm_inc_spillover_mosq_unadj.RDS"))
saveRDS(res_direct_mosq, file = paste0(results_path, "namibia_glm_inc_direct_mosq_unadj.RDS"))


# human & mosquito reservoir ---------------------------------------
res_total_hm <- lapply(data_hm_list, function(x)
  fit_indiv_glm_unadj(data = x,
                      effecttype = "total"))

res_spillover_hm <- lapply(data_hm_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "spillover"))

res_direct_hm <- lapply(data_hm_list, function(x)
  fit_indiv_glm_unadj(data = x,
                effecttype = "direct"))

# save results
saveRDS(res_total_hm, file = paste0(results_path, "namibia_glm_inc_total_hm_unadj.RDS"))
saveRDS(res_spillover_hm, file = paste0(results_path, "namibia_glm_inc_spillover_hm_unadj.RDS"))
saveRDS(res_direct_hm, file = paste0(results_path, "namibia_glm_inc_direct_hm_unadj.RDS"))

