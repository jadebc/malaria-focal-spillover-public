# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create analysis data structure
# ............................................


rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))
library(fastDummies)


# ............................................
# Load data
# ............................................
cohorts <- readRDS(namibia_cohort_path)


# ............................................
# Generate individual-level covariates (W)
# ............................................
## 1. Generate individual-level covariates (W) ------------------------------------------
index_int <- readRDS(namibia_sub_coho_index_int_path)
int_recip <- readRDS(namibia_sub_coho_int_recip_path)
list_people <- readRDS(namibia_coho_person_list_path)

## read in cohort classification for spillover zone radii of interest ------------------------------------------
df_base <- readRDS(namibia_df_base_path)
list_people_coho <- unique(df_base[,c("indiv_id", "longitude", "latitude")])


### 1.1 count_int_pre_recei ------------------------------------------
df_count_int_pre_recei <- count_int_pre_recei(cohorts = cohorts,
                                              df_base = df_base)
saveRDS(df_count_int_pre_recei, file = count_int_pre_recei_path)
df_count_int_pre_recei = readRDS(count_int_pre_recei_path)

### 1.2 count_int_pre_deliv ------------------------------------------
df_count_int_pre_deliv <- count_int_pre_deliv(index_int = index_int,
                                              int_recip = int_recip,
                                              df_base = df_base,
                                              list_people_coho = list_people_coho)

dist_cols_pre_deliv <- grep("pre_int", colnames(df_count_int_pre_deliv))
df_count_int_pre_deliv <- df_count_int_pre_deliv %>%
  mutate_at(dist_cols_pre_deliv, ~replace(., is.na(.), 0))

saveRDS(df_count_int_pre_deliv,  file = count_int_pre_deliv_path)
df_count_int_pre_deliv = readRDS(count_int_pre_deliv_path)

### 1.3 count_int_concu_deliv ------------------------------------------
df_count_int_concu_deliv <- count_int_concu_deliv (r = c(500, 1000, 2000, 3000),
                                                   index_int = index_int,
                                                   int_recip = int_recip,
                                                   df_base = df_base,
                                                   list_people_coho = list_people_coho,
                                                   t1_l = 0,
                                                   t1_u =35,
                                                   t2_l = 0,
                                                   t2_u =180)

# temp1 <- df_count_int_concu_deliv %>% mutate_at(c(10:25), ~replace(., is.na(.), 0))
dist_cols_concu_deliv <- grep("pre_int", colnames(df_count_int_concu_deliv))
df_count_int_concu_deliv <- df_count_int_concu_deliv %>%
  mutate_at(dist_cols_concu_deliv, ~replace(., is.na(.), 0))

saveRDS(df_count_int_concu_deliv, file = count_int_concu_deliv_path)
df_count_int_concu_deliv = readRDS(count_int_concu_deliv_path)

### 1.4 count_pre_index ------------------------------------------
df_count_pre_index <- count_pre_index(df_base = df_base,
                                      cohorts = cohorts)

saveRDS(df_count_pre_index, count_pre_index_path)
df_count_pre_index = readRDS(count_pre_index_path)

### 1.5 count_pop ------------------------------------------
df_count_pop <- count_pop(cohorts = cohorts)

saveRDS(df_count_pop,  file = count_pop_path)
df_count_pop = readRDS(count_pop_path)


### 1.6 age and gender ------------------------------------------
age_gender <- cohorts %>% select(c("indiv_id","age","sex")) %>% unique()



# ............................................
# Skim individual-level covariates (W)
# ............................................
summary(df_count_int_pre_recei)

get_summary <- function(df){
  df_summary <- data.frame('var_name' = names(df)[-1],
                           'min' = NA, 
                           'q05' = NA,
                           'median' = NA, 
                           'q95' = NA,
                           'max' = NA)
  
  for (i in 1:nrow(df_summary)){
    df_summary[i,2] = min(df[,i + 1])
    df_summary[i,3] = quantile(df[,i + 1], probs = c(0.05))
    df_summary[i,4] = median(df[,i + 1])
    df_summary[i,5] =  quantile(df[,i + 1], probs = c(0.95))
    df_summary[i,6] = max(df[,i + 1])
  }
  return(df_summary)
}

# get_summary(df_count_int_pre_recei)
# get_summary(df_count_int_pre_deliv)
# get_summary(df_count_int_concu_deliv)
# get_summary(df_prop_pre_int)
# get_summary(df_count_pre_index)
# get_summary(df_count_pop)


# ............................................
# Genarate analysis data
# ............................................
# merge individual-level covariates (W)
df_base <- readRDS(namibia_df_base_path)
df_count_int_pre_recei = readRDS(count_int_pre_recei_path)
df_count_int_pre_deliv <- readRDS(count_int_pre_deliv_path) 
df_count_int_concu_deliv <- readRDS(count_int_concu_deliv_path)
df_count_pre_index <- readRDS(count_pre_index_path)
df_count_pop <- readRDS(count_pop_path)

df_all <- left_join(df_count_int_pre_deliv, 
                    df_count_int_concu_deliv,
                    by = c("indiv_id","cohort_id",
                           "date", "intarm", 
                           "longitude", "latitude",
                           "dist_to_index", "target_area", "int_recip"))

df_all <- left_join(df_all, df_count_int_pre_recei, by = c("indiv_id", "cohort_id", "date"))
df_all <- left_join(df_all, df_count_pop, by = c("indiv_id"))
df_all <- left_join(df_all, df_count_pre_index, by = c("indiv_id", "cohort_id", "date"))

sub_coho_index_nonint <- readRDS(namibia_sub_coho_index_nonint_path)
sub_coho_index_nonint <- sub_coho_index_nonint %>% 
  dplyr::select(c("indiv_id", "indexcase", "date")) %>% 
  unique() %>% 
  mutate(infected_date = date) %>% 
  dplyr::select(c("indiv_id", "indexcase", "infected_date"))

df_all <- left_join(df_all, sub_coho_index_nonint, by = c("indiv_id"))
df_all <- df_all %>% mutate(indexcase = ifelse(is.na(indexcase), 0, indexcase))

df_all <- df_all %>% 
  mutate(t0 = case_when(target_area == 1 ~ date,
                        target_area != 1 & (intarm == "RO" | intarm == "TO") ~ (date + 21),
                        target_area != 1 & (intarm == "RV" | intarm == "TV") ~ (date + 17)),
         
         t1 = case_when(target_area == 1 & (intarm == "RO" | intarm == "TO") ~ (date + 35),
                        target_area == 1 & (intarm == "RV" | intarm == "TV") ~ (date + 180),
                        target_area != 1 & (intarm == "RO" | intarm == "TO") ~ (date + 56),
                        target_area != 1 & (intarm == "RV" | intarm == "TV") ~ (date + 180))
  )

df_all_index <- df_all %>% filter(indexcase == 1)
df_all_nonindex <- df_all %>% filter(indexcase != 1)

df_all_index <- df_all_index %>% mutate(indexcase = ifelse(infected_date >= t0 & 
                                                             infected_date <= t1,
                                                           indexcase, 0))

df_all_index <- df_all_index %>%  group_by(cohort_id, indiv_id) %>% mutate(occurence = n()) 

df_all_index <- df_all_index %>% 
  group_by(indiv_id, cohort_id) %>% 
  arrange(desc(indexcase), infected_date) %>% 
  filter(row_number() == 1) %>% 
  select(-c("occurence"))

df_all <- bind_rows(df_all_index, df_all_nonindex)

df_all$infected_date[df_all$indexcase == 0] <- NA

# add time-to-event variable

df_all <- df_all %>% mutate(tte = ifelse(indexcase == 1, 
                                         infected_date - t0, 
                                         t1 - t0))

age_gender <- cohorts %>% select(c("indiv_id","age","sex")) %>% unique()
age_gender <- age_gender[complete.cases(age_gender),]
age_gender <- age_gender %>% group_by(indiv_id) %>% filter(row_number()==1)
df_all <- left_join(df_all, age_gender, by = c("indiv_id"))

# ............................................
# Import cohort-level covariates (E)
# ............................................
## 2. Generate cohort-level covariates (E) ------------------------------------------
### 2.1 EA id ------------------------------------------
ea_id_merge <- cohorts %>% 
  mutate(cohort_id = iid_combined) %>%
  filter(!is.na(cohort_id) & cohort_id!= "") %>% 
  dplyr::select(c("cohort_id", "ea")) %>% unique()

ea_id_merge <- ea_id_merge %>% 
  mutate(ea = ifelse(cohort_id == "RO27-6401", "RO27", ea)) %>% 
  unique()

df_all <- left_join(df_all, ea_id_merge, by = c("cohort_id"))


### 2.2 - 2.10 ------------------------------------------
ea_level_raw <- read_df(read.csv, namibia_raw_ea_level_path)
df_ea <- ea_level_raw %>% 
  filter(ea_actual != "") %>% 
  mutate(ea = ea_actual,
         ### 2.2 malaria incidence in the season prior to the trial ------------------------------------------
         pre_incidence = ea_incidence_2016,
         ### 2.3 pre-season indoor residual spray coverage ------------------------------------------
         pre_spray_cover = mohss_spray_coverage_ea,
         ### 2.4 population size of EA ------------------------------------------
         pop_size_ea = ea_pop_update,
         ### 2.5 distance to nearest neighboring household  ------------------------------------------
         dist_hh = mean_nearest_hh_distance_meters,
         ### 2.6 distance to nearest healthcare facility  ------------------------------------------
         dist_hf = mean_hh_distance_to_hf,
         ### 2.7 median monthly rainfall   ------------------------------------------
         pre_rainfall = meanprecmedianealag2,
         ### 2.8 median enhanced vegetative index   ------------------------------------------
         pre_evi = medianevimedianealag2,
         ### 2.9 median elevation   ------------------------------------------
         ea_elevation = elevmedianea,
         ### 2.10 median daytime land surface temperature   ------------------------------------------
         surface_temp = medianlstmedianealag2_fixed) %>% 
  dplyr::select(c("ea",
                  "pre_incidence",
                  "pre_spray_cover",
                  "pop_size_ea",
                  "dist_hh",
                  "dist_hf",
                  "pre_rainfall",
                  "pre_evi",
                  "ea_elevation",
                  "surface_temp"))

df_all <- left_join(df_all, df_ea, by = c("ea"))

### 2.11 time lag from index case detection to intervention in the cohort  ------------------------------------------

# only need matched people to show on the map
sub_coho <- cohorts %>% filter(evermatchday == 1)

# replace NA with 0 for some variables
sub_coho$int_trig[which(is.na(sub_coho$int_trig))] <- 0
sub_coho$clinic_day[which(is.na(sub_coho$clinic_day))] <- 0
sub_coho$int_day[which(is.na(sub_coho$int_day))] <- 0

sub_coho <- sub_coho %>%
  mutate(person_type = case_when(clinic_day == 1 & int_trig == 1 ~ "index_int_clinic",
                                 int_day == 1 & int_trig == 1 ~ "index_int_visit",
                                 indexcase == 1 & int_trig == 0 ~ "index_nonint",
                                 TRUE ~ "int_recip"))

sub_coho_index_clinic <- sub_coho %>% 
  filter(person_type == "index_int_clinic") %>% 
  mutate(clinic_date = date) %>% 
  dplyr::select(c("indiv_id", "clinic_date"))
index_int_visit <- readRDS(namibia_sub_coho_index_int_path)

index_int_meantime_merge <- left_join(index_int_visit, sub_coho_index_clinic, by ="indiv_id") %>% 
  mutate(cohort_id = iid_combined) %>% 
  dplyr::select(c("cohort_id", "clinic_date", "int_date"))
index_int_meantime_merge <- index_int_meantime_merge %>% 
  mutate(time_lag = as.numeric(int_date - clinic_date)) %>% 
  dplyr::select(c("cohort_id", "time_lag"))
index_int_meantime_merge <- index_int_meantime_merge %>% 
  mutate(time_lag = ifelse(time_lag < 0, 0, time_lag))

df_all <- left_join(df_all, index_int_meantime_merge, by = c("cohort_id"))

### 2.12 population size of cohorts ------------------------------------------

pop_cohort_merge <- df_base %>% 
  group_by(cohort_id) %>%   
  summarise(pop_size_cohort = n())

df_all <- left_join(df_all, pop_cohort_merge, by = c("cohort_id"))


# add an indicator for index cases that triggered the intervention
temp <- cohorts %>% filter(int_trig == 1) %>% 
  mutate("cohort_id"= iid_combined) %>% 
  select(c("cohort_id","indiv_id","int_trig")) %>% 
  unique()

df_all <- left_join(df_all, temp, by = c("indiv_id", "cohort_id"))
df_all$int_trig[is.na(df_all$int_trig)] <- 0


# distance to closest health facility data  ------------------
hf <- read.csv(paste0(box_shared_path, "Datasets/Final trial data for collaborators/11_health_facilities_location.csv"))
hf_ll <- hf %>% dplyr::select(Longitude, Latitude)
hf_ids <-  hf %>% pull(FacID)

dist_to_hf <- function(data, id_name, id_value){
  
  ll = data %>% 
    ungroup() %>% 
    dplyr::filter(!!rlang::sym(id_name) == id_value) %>% 
    dplyr::select(longitude, latitude) %>% 
    distinct()
  
  dists <- apply(hf_ll, 1, function(x) distm(x, as.matrix(ll)))
  
  # save in km
  min_dist <- dists[which.min(dists)]/1000
  
  out <- data.frame(min_dist_hf = min_dist) %>% 
    mutate({{id_name}} := id_value) 
  
  return(out)
}

list_inc_ids <- as.list(unique(df_all$indiv_id))

result_inc <- lapply(list_inc_ids, function(x) 
  dist_to_hf(data = df_all, id_name = "indiv_id", id_value = x)) %>% bind_rows()

df_all <- left_join(df_all, result_inc, by = "indiv_id")

median_min_dist_hf <- median((df_all$min_dist_hf),na.rm=T)

df_all <- df_all %>% mutate(
  min_dist_hf_abovemed = ifelse(min_dist_hf > median_min_dist_hf, 1, 0)
)


df_all <- df_all %>% ungroup

## define variable class
df_all$cohort_id <- as.factor(df_all$cohort_id)
df_all$intarm <- as.factor(df_all$intarm)

df_all$sex <- as.character(df_all$sex)
df_all$sex[is.na(df_all$sex)] <- "unknown"
df_all$sex <- as.factor(df_all$sex)

saveRDS(df_all, namibia_df_all_path)

