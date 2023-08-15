################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Inputs for flowchart
################################################ 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ------
df_all <- readRDS(namibia_df_all_path)
index_df <- readRDS(namibia_clean_index_path)


# to ensure that N cohorts is same for direct and spillover ---------------------------------------
# drop the 11 cohorts that either 1) never had int_recip==1 or 
# 2) had int_recip==1 but they were never in a target area
direct_ids = df_all %>% filter(target_area == 1, int_recip == 1) %>% 
  dplyr::select(cohort_id) %>% distinct() %>% pull(cohort_id)

spill_ids = df_all %>%  filter(target_area == 0 | (target_area == 1 & int_recip == 0)) %>% 
  dplyr::select(cohort_id) %>% distinct()%>% pull(cohort_id)

# ids that are not in direct but are in spill 
drop_ids <- spill_ids[!spill_ids %in% direct_ids]

df_sub <- df_all %>% filter(!cohort_id %in% drop_ids)


# get counts for flowchart ---------------------------------------


# number of EAs with interventions 
N_EA_per_arm <- index_df %>% filter(first_case_yn == 1) %>% 
  dplyr::select(id, intarm) %>% 
  group_by(intarm) %>% 
  summarise(N = n())

# number of clusters with index cases that triggered interventions per arm 
Nclus_withindex_per_arm <- index_df %>% filter(first_case_yn == 1) %>% 
  dplyr::select(ea_index, intarm) %>% 
  group_by(intarm) %>% 
  summarise(N = length(unique(ea_index)))

# number of index cases that triggered interventions per arm 
Nindex_per_arm <- index_df %>% filter(first_case_yn == 1) %>% 
  dplyr::select(ea_index, intarm) %>% 
  group_by(intarm) %>% 
  summarise(N = n())

# number of cohorts per arm after subsetting 
Ncohort_sub_per_arm <- df_sub %>% dplyr::select(cohort_id, intarm) %>% 
  group_by(intarm) %>% 
  distinct() %>% 
  summarise(N = n())

# number of cohorts per arm dropped because no intervention recip nearby
cohort_drop = left_join(Nindex_per_arm, 
                      Ncohort_sub_per_arm %>% rename(Nsub=N), 
                      by = "intarm") %>% 
  mutate(diff = N-Nsub) 

# number of clusters included in analysis in each arm
Nclus_per_arm <- df_all %>% dplyr::select(ea, intarm) %>% 
  group_by(intarm) %>% 
  distinct() %>% 
  summarise(N = n())

# N population per arm
pop_per_arm <- df_sub %>% dplyr::select(indiv_id, intarm) %>% 
  group_by(intarm) %>% 
  distinct() %>% 
  summarise(Npop = n())

# N non triggering index cases per arm
N_index_nontrig_per_arm <- df_sub %>% 
  filter(indexcase==1) %>% 
  dplyr::select(indiv_id, intarm) %>% 
  group_by(intarm) %>% 
  distinct() %>% 
  summarise(Npop = n())



