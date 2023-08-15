################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# table of number of interventions delivered per cluster and per arm 
################################################ 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ------
df_all <- readRDS(namibia_df_all_path)

N_int_clus_arm <- df_all %>% dplyr::select(intarm, ea, cohort_id) %>% 
  distinct() %>% 
  group_by(intarm, ea) %>% 
  summarise(N_int = n()) %>% 
  group_by(intarm) %>% 
  summarise(min_int = min(N_int),
            mean_int = mean(N_int), 
            max_int = max(N_int)) %>% 
  mutate(mean_int = sprintf("%0.01f", mean_int))

write.csv(N_int_clus_arm, file = paste0(table_path, "table-interventions-per-arm.csv"))


# reported in text
# number of interventions delivered per cluster and per arm 
N_int_clus_arm <- df_all %>% dplyr::select(intarm, ea, cohort_id) %>% 
  distinct() %>% 
  group_by(intarm, ea) %>% 
  summarise(N_int = n()) %>% 
  group_by(intarm) %>% 
  summarise(min_int = min(N_int),
            mean_int = mean(N_int), 
            max_int = max(N_int))

N_int_clus <- df_all %>% dplyr::select(ea, cohort_id) %>% 
  distinct() %>% 
  group_by(ea) %>% 
  summarise(N_int = n()) %>% 
  ungroup() %>% 
  summarise(min_int = min(N_int),
            mean_int = mean(N_int), 
            max_int = max(N_int))




