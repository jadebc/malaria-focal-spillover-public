# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Classify cohorts using different
# radii for spillover zones
# ............................................

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))

install.packages("fastDummies",repos = "http://cran.us.r-project.org")

library(fastDummies)

# load data ------------------------------------------
cohorts <- readRDS(namibia_cohort_path)
list_people <- unique(cohorts[,c("indiv_id", "longitude", "latitude")])

index_int <- readRDS(namibia_sub_coho_index_int_path)
int_recip <- readRDS(namibia_sub_coho_int_recip_path)
list_people <- readRDS(namibia_coho_person_list_path)

# classify cohort and set spillover zone radius ------------------------------------------
df_base <- coho_classify(r = 1000, index_int, list_people)
saveRDS(df_base, namibia_df_base_path)

df_base_2km <- coho_classify(r = 2000, index_int, list_people)
saveRDS(df_base_2km, namibia_df_base_path_sens_spzone_2km)

df_base_3km <- coho_classify(r = 3000, index_int, list_people)
saveRDS(df_base_3km, namibia_df_base_path_sens_spzone_3km)


