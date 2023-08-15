# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create analysis data structure
# Short observation period for 
# interventions on human reservoir

# Sensitivity analysis removing overlap
# between target areas and spillover zones
# ............................................


rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))
library(fastDummies)


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

# intersect(sub_coho_int_recip$indiv_id, sub_coho_index_nonint$indiv_id)

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



## 09-08-2021 for df_short, re-calc this variable 
### 1.3 count_int_concu_deliv ------------------------------------------
df_count_int_concu_deliv <- count_int_concu_deliv (r = c(500, 1000, 2000, 3000),
                                                   index_int = index_int,
                                                   int_recip = int_recip,
                                                   df_base = df_base,
                                                   list_people_coho = list_people_coho,
                                                   t1_l = 0,
                                                   t1_u =35,
                                                   t2_l = 0,
                                                   t2_u =35)

dist_cols_concu_deliv <- grep("pre_int", colnames(df_count_int_concu_deliv))
df_count_int_concu_deliv <- df_count_int_concu_deliv %>%
  mutate_at(dist_cols_concu_deliv, ~replace(., is.na(.), 0))
saveRDS(df_count_int_concu_deliv, file = count_int_concu_deliv_s_path)

# ............................................
# Genarate analysis data
# ............................................
# merge individual-level covariates (W)
df_base <- readRDS(namibia_df_base_path)
df_count_int_pre_recei = readRDS(count_int_pre_recei_path)
df_count_int_pre_deliv <- readRDS(count_int_pre_deliv_path) 

# load df_short version "df_count_int_concu_deliv"
df_count_int_concu_deliv <- readRDS(count_int_concu_deliv_s_path)

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

# define observation period -----------------------------------------------
df_all <- df_all %>% 
                mutate(t0 = case_when(target_area == 1 & int_recip == 1 ~ date,
                                      target_area == 1 & int_recip == 0 ~ (date + 21),
                                      target_area != 1 ~ (date + 21)),
                                      
                       t1 = case_when(target_area == 1 & int_recip == 1 ~ (date + 35),
                                      target_area == 1 & int_recip == 0 ~ (date + 56),
                                      target_area != 1 ~ (date + 56))
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

# add age and gender
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


### 2.2 - 2.10
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

df_all <- df_all %>% ungroup

## 3. define variable class------------------------------------------
df_all$cohort_id <- as.factor(df_all$cohort_id)
df_all$intarm <- as.factor(df_all$intarm)

df_all$sex <- as.character(df_all$sex)
df_all$sex[is.na(df_all$sex)] <- "unknown"
df_all$sex <- as.factor(df_all$sex)

## 4. generate variables for effect modification assessment ------------------------------------------
### calculate percentage of cohort treated -------------
cohort_cov <- df_all %>% group_by(cohort_id) %>% 
  summarise(Ntx = sum(int_recip),
            Npop = mean(pop_size_cohort)) %>% 
  mutate(tx_cov_cohort = Ntx/Npop)

df_all <- left_join(df_all, cohort_cov, by = "cohort_id")

median_tx_cov_cohort <- median((df_all$tx_cov_cohort),na.rm=T)

df_all <- df_all %>% mutate(
  tx_cov_cohort_abovemed = ifelse(tx_cov_cohort > median_tx_cov_cohort, 1, 0)
)

### calculate response time -------------
# merge in index case date
index_data <- readRDS(namibia_clean_index_path) %>% 
  filter(first_case_yn == 1) %>% 
  dplyr::select(iid_combined, clinic_date_index) %>% 
  distinct() %>% 
  rename(cohort_id = iid_combined)

df_all <- left_join(df_all, index_data, by = "cohort_id") %>% 
  mutate(response_time = date - clinic_date_index) %>% 
  mutate(response_time = as.numeric(response_time)) %>% 
  mutate(response_time = ifelse(response_time<0, NA, response_time)) 

median_response_time <- median((df_all$response_time),na.rm=T)

df_all <- df_all %>% 
  mutate(response_time_abovemed = ifelse(
    response_time >median_response_time,1,0 ))

### distance to index case -------------

median_dist_to_index <- median((df_all$dist_to_index),na.rm=T)
df_all <- df_all %>% 
  mutate(dist_index_abovemed = ifelse(
    dist_to_index >median_dist_to_index,1,0 )) %>% 
  mutate(dist_index_quant = cut(
    dist_to_index, breaks = quantile(df_all$dist_to_index, probs = c(0,.25,.5,.75,1))
  ))



## create above/below median indicators for other EM variables -------------
df_all <- df_all %>% 
  mutate(pre_incidence_abovemed = ifelse(
    pre_incidence > median(df_all$pre_incidence, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_spray_cover_abovemed = ifelse(
    pre_spray_cover > median(df_all$pre_spray_cover, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_rainfall_abovemed = ifelse(
    pre_rainfall > median(df_all$pre_rainfall, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_evi_abovemed = ifelse(
    pre_evi > median(df_all$pre_evi, na.rm=T), 1, 0
  )) %>% 
  mutate(surface_temp_abovemed = ifelse(
    surface_temp > median(df_all$surface_temp, na.rm=T), 1, 0
  )) %>% 
  mutate(ea_elevation_abovemed = ifelse(
    ea_elevation > median(df_all$ea_elevation, na.rm=T), 1, 0
  ))

## high transmission season -------------
df_all <- df_all %>% mutate(high_season = ifelse(month(date) >=2 & month(date) <=5 , 1, 0))

## Proportion of population previously treated  -------------
prop_prev_tx_df <- df_all %>% 
  mutate(pre_recei_RV = ifelse(n_int_pre_recei_RV>0, 1, 0),
         pre_recei_RO = ifelse(n_int_pre_recei_RO>0, 1, 0),
         pre_recei_TV = ifelse(n_int_pre_recei_TV>0, 1, 0),
         pre_recei_TO = ifelse(n_int_pre_recei_TV>0, 1, 0)) %>% 
  mutate(n_int_pre_received = 
           pre_recei_RV+ pre_recei_RO +  
           pre_recei_TV+ pre_recei_TO) %>% 
  group_by(cohort_id) %>% 
  summarise(
    n_prev_tx = sum(n_int_pre_received),
    NPop = mean(pop_size_cohort)
  ) %>% 
  mutate(prop_prev_tx = n_prev_tx / NPop) %>% 
  mutate(prop_prev_tx = ifelse(prop_prev_tx>1, 1, prop_prev_tx)) %>% 
  dplyr::select(-NPop)
  
df_all <- df_all %>% left_join(prop_prev_tx_df, by = "cohort_id")


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

saveRDS(df_all, namibia_df_short_path)
  
## 5. permute intervention at ea level------------------------------------------
set.seed(123)
df_all <- readRDS(namibia_df_short_path)
temp_ea <- df_all %>% select(ea, intarm) %>% unique()
temp_ea$intarm <- sample(temp_ea$intarm, size = length(temp_ea$intarm))
df_all_permute <- df_all %>% select(-intarm)
df_all_permute <- left_join(df_all_permute, temp_ea, by = "ea") 

saveRDS(df_all_permute, namibia_df_short_permute_path)

## make de-identified version ------------------------------------------
df_all <- readRDS(namibia_df_short_path) 
df_deid <- df_all %>% dplyr::select(-c(longitude, latitude, age, sex))

saveRDS(df_deid, namibia_df_short_deid_path)

