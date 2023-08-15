################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Primary analysis: incidence
################################################ 
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
data_human_list <- readRDS(namibia_human_process_path)
data_mosq_list <- readRDS(namibia_mosq_process_path)
data_hm_list <- readRDS(namibia_hm_process_path)

data_human_sens_list <- readRDS(namibia_human_process_sens_obs_path)
data_mosq_sens_list <- readRDS(namibia_mosq_process_sens_obs_path)
data_hm_sens_list <- readRDS(namibia_hm_process_sens_obs_path)



get_percent_overlap <- function(data){
  overlap_output <- check_overlap_all(df=data, r=1000) 
  
  overlap_summary <- overlap_output %>% 
    mutate(overlap_spill=ifelse(overlap_type=="spill", 1, 0),
           overlap_target=ifelse(overlap_type=="target", 1, 0)) %>% 
    group_by(coho_i) %>% 
    summarise(n_spill_overlap = sum(overlap_spill),
              n_target_overlap = sum(overlap_target)) %>% 
    rename(id = coho_i)
  
  overlap <- left_join(data, overlap_summary, by = "id") %>% 
    dplyr::select(id, n_spill_overlap, n_target_overlap) %>% 
    distinct() %>% 
    mutate(any_overlap_target = ifelse(n_target_overlap>1, 1, 0),
           any_overlap_spill = ifelse(n_spill_overlap>1, 1, 0))
  
  target <- sprintf("%0.01f", mean(overlap$any_overlap_target, na.rm=T)*100)
  spill <- sprintf("%0.01f", mean(overlap$any_overlap_spill, na.rm=T)*100)
  
  return(data.frame(target=target, spill=spill))
  
}

human <- get_percent_overlap(data_human_list$short)
mosq <- get_percent_overlap(data_mosq_list$long)
hm <- get_percent_overlap(data_hm_list$long)

human_sens <- get_percent_overlap(data_human_sens_list$short)
mosq_sens <- get_percent_overlap(data_mosq_sens_list$long)
hm_sens <- get_percent_overlap(data_hm_sens_list$long)

human_row <- cbind(human, human_sens)
mosq_row <- cbind(mosq, mosq_sens)
hm_row <- cbind(hm, hm_sens)

colnames(human_row) <- c("target", "spill", "target-sens", "spill-sens")
colnames(mosq_row) <- c("target", "spill", "target-sens", "spill-sens")
colnames(hm_row) <- c("target", "spill", "target-sens", "spill-sens")

table <- bind_rows(human_row, mosq_row, hm_row) %>% 
  mutate(comparison = c("Human", "Mosquito", "Human & mosquito")) %>% 
  dplyr::select(comparison, everything())

write.csv(table, file = paste0(table_path, "table-cohort-overlap.csv"))

