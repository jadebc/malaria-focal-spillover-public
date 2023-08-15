# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create matched cohort subsets
# ............................................

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))

# ............................................
# Load data
# ............................................
cohorts <- readRDS(namibia_cohort_path)

# only need matched people to show on the map
sub_coho <- cohorts %>% filter(evermatchday == 1)

# ............................................
# Data process
# ............................................
# replace NA with 0 for some variables
sub_coho$int_trig[which(is.na(sub_coho$int_trig))] <- 0
sub_coho$clinic_day[which(is.na(sub_coho$clinic_day))] <- 0
sub_coho$int_day[which(is.na(sub_coho$int_day))] <- 0

sub_coho <- sub_coho %>%
  mutate(person_type = case_when(clinic_day == 1 & int_trig == 1 ~ "index_int_clinic",
                                 int_day == 1 & int_trig == 1 ~ "index_int_visit",
                                 indexcase == 1 & int_trig == 0 ~ "index_nonint",
                                 TRUE ~ "int_recip"))

sub_coho_index_int <- sub_coho %>% filter(int_day == 1 & int_trig == 1)
sub_coho_int_recip <- sub_coho %>% filter(person_type == "int_recip")
sub_coho_index_nonint <- sub_coho %>% filter(person_type == "index_nonint")

# generate a "int_date" (min(visit_day, first_int_day)) to indicate the intervention day of each cohort
list_iid_combined <- sub_coho_index_int$iid_combined
int_date <- rep(NA, length(list_iid_combined))
for (i in 1:length(list_iid_combined)){
  int_recip <- sub_coho_int_recip %>% filter(iid_combined == list_iid_combined[i])
  if (nrow(int_recip) == 0){
    int_date[i] = sub_coho_index_int$date[i]
  }else{
    int_date[i] = min(int_recip$date)
  }
}
int_date <- as.Date(int_date, origin = "1970-01-01")

sub_coho_index_int$int_date <- int_date
list_people <- unique(cohorts[,c("indiv_id", "longitude", "latitude")])


saveRDS(sub_coho_index_int,
        namibia_sub_coho_index_int_path)
saveRDS(list_people,
        namibia_coho_person_list_path)
saveRDS(sub_coho_int_recip,
        namibia_sub_coho_int_recip_path)
saveRDS(sub_coho_index_nonint,
        namibia_sub_coho_index_nonint_path)

