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

x = readRDS(paste0(results_path,  "namibia_htmle_inc_direct_hm_long_adaptQ.RDS"))
y = list(
  Long = x,
  Short = NA
)
saveRDS(y, paste0(results_path, "namibia_htmle_inc_direct_hm_adaptQ.RDS"))

file_list <- c("namibia_htmle_inc_direct_mosq_adaptQ.RDS",
               "namibia_htmle_inc_direct_human_adaptQ.RDS",
               "namibia_htmle_inc_direct_hm_adaptQ.RDS",
               "namibia_htmle_inc_spillover_hm_adaptQ.RDS",
               "namibia_htmle_inc_spillover_mosq_adaptQ.RDS",
               "namibia_htmle_inc_spillover_human_adaptQ.RDS",
               "namibia_htmle_inc_total_mosq_adaptQ.RDS",
               "namibia_htmle_inc_total_human_adaptQ.RDS",
               "namibia_htmle_inc_total_hm_adaptQ.RDS")



read_adaptQ_results <- function(file_path){
  
  results <- readRDS(paste0(results_path, file_path)) %>% unlist()
  df = as.data.frame(t(results))
  
  if(length(grep("_human_", file_path)==1)) reservoir = "Human"
  if(length(grep("_mosq_", file_path)==1)) reservoir = "Mosquito"
  if(length(grep("_hm_", file_path)==1)) reservoir = "Human & mosquito"
  
  if(length(grep("spillover", file_path)==1)) parameter = "spillover"
  if(length(grep("total", file_path)==1)) parameter = "total"
  if(length(grep("direct", file_path)==1)) parameter = "direct"
  
  df$reservoir = reservoir
  df$parameter = parameter
  df = df %>% dplyr::select(reservoir, parameter, everything())
  
  return(df)
}

# primary analysis --------------------------------------------------------
estimates_list <- lapply(file_list, read_adaptQ_results)
estimates <- bind_rows(estimates_list) %>% arrange(reservoir, parameter)

saveRDS(estimates, file = paste0(results_path, "namibia_htmle_inc_adaptQ.RDS"))




