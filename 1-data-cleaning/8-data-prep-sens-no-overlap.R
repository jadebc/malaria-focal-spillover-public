# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create analysis data structure
# Sensitivity analysis 
# of cohort with no overlap
# ............................................


rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))
library(fastDummies)

# load data --------------------------------------------
short <- readRDS(namibia_df_short_path) %>% 
  rename(A = intarm,
         id = cohort_id)
long <- readRDS(namibia_df_long_path) %>% 
  rename(A = intarm,
         id = cohort_id)


# identify overlapping cohorts --------------------------------------------
overlap_short <- check_overlap_all(df = short, r = 1000)
overlap_long <- check_overlap_all(df = long, r = 1000)


# list overlapping cohorts --------------------------------------------
overlap_short_target_df <- overlap_short %>% 
  filter(overlap_type == "target") %>% 
  dplyr::select(coho_i, coho_j)

overlap_short_target <- c(
  overlap_short_target_df[,1],
  overlap_short_target_df[,2]
) %>% unique() %>% as.character()

overlap_long_target_df <- overlap_long %>% 
  filter(overlap_type == "target") %>% 
  dplyr::select(coho_i, coho_j)

overlap_long_target <- c(
  overlap_long_target_df[,1],
  overlap_long_target_df[,2]
) %>% unique() %>% as.character()

overlap_short_spill_df <- overlap_short %>% 
  filter(overlap_type == "spill") %>% 
  dplyr::select(coho_i, coho_j)

overlap_short_spill <- c(
  overlap_short_spill_df[,1],
  overlap_short_spill_df[,2]
) %>% unique() %>% as.character()

overlap_long_spill_df <- overlap_long %>% 
  filter(overlap_type == "spill") %>% 
  dplyr::select(coho_i, coho_j)

overlap_long_spill <- c(
  overlap_long_spill_df[,1],
  overlap_long_spill_df[,2]
) %>% unique() %>% as.character()

# create datasets with no overlapping cohorts --------------------------------------------

nonoverlap_short_target <- short %>% 
  filter(id %in% overlap_short_target) %>% 
  rename(intarm = A,
         cohort_id = id)

nonoverlap_long_target <- long %>% 
  filter(id %in% overlap_long_target) %>% 
  rename(intarm = A,
         cohort_id = id)

nonoverlap_short_spill <- short %>% 
  filter(id %in% overlap_short_spill) %>% 
  rename(intarm = A,
         cohort_id = id)

nonoverlap_long_spill <- long %>% 
  filter(id %in% overlap_long_spill) %>% 
  rename(intarm = A,
         cohort_id = id)

# save datasets --------------------------------------------

saveRDS(nonoverlap_short_target, file = namibia_df_short_sens_nooverlap_target_path)
saveRDS(nonoverlap_short_spill, file = namibia_df_short_sens_nooverlap_spill_path)
saveRDS(nonoverlap_long_target, file = namibia_df_long_sens_nooverlap_target_path)
saveRDS(nonoverlap_long_spill, file = namibia_df_long_sens_nooverlap_spill_path)
