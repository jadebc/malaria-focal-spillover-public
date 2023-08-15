################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Preprocess data for 
# Primary analysis
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ------
df_long <- readRDS(namibia_df_long_path)
df_short <- readRDS(namibia_df_short_path)

# add month variable -----
df_long <- df_long %>% mutate(month = month(date))
df_short <- df_short %>% mutate(month = month(date))

# to ensure that N cohorts is same for direct and spillover ---------------------------------------
# drop the 11 cohorts that either 1) never had int_recip==1 or 
# 2) had int_recip==1 but they were never in a target area
long_direct_ids = df_long %>% filter(target_area == 1, int_recip == 1) %>% 
  dplyr::select(cohort_id) %>% distinct() %>% pull(cohort_id)

long_spill_ids = df_long %>%  filter(target_area == 0 | (target_area == 1 & int_recip == 0)) %>% 
  dplyr::select(cohort_id) %>% distinct()%>% pull(cohort_id)

short_direct_ids = df_short %>% filter(target_area == 1, int_recip == 1) %>% 
  dplyr::select(cohort_id) %>% distinct() %>% pull(cohort_id)

short_spill_ids = df_short %>%  filter(target_area == 0 | (target_area == 1 & int_recip == 0)) %>% 
  dplyr::select(cohort_id) %>% distinct()%>% pull(cohort_id)

# ids that are not in direct but are in spill 
long_drop_ids <- long_spill_ids[!long_spill_ids %in% long_direct_ids]
short_drop_ids <- short_spill_ids[!short_spill_ids %in% short_direct_ids]

df_long_sub <- df_long %>% filter(!cohort_id %in% long_drop_ids)
df_short_sub <- df_short %>% filter(!cohort_id %in% short_drop_ids)

# preprocess data ------
df_long_process <- preprocess_data(df_long_sub)
df_short_process <- preprocess_data(df_short_sub)

# define covariates ---------------------------------------

## names of covariates to drop due to collinearity ------------------
dropnames = c(
  "delta_concur_int_1000",
  "delta_concur_int_500",             
  "delta_concur_int_propdiffarm_500",
  "delta_concur_int_propdiffarm_1000",
  'concur_int_propdiffarm_500',
  'concur_int_propdiffarm_1000',
  'concur_int_propdiffarm_2000',
  'concur_int_propdiffarm_3000',
  'concur_int_3000',
  'concur_int_2000',
  'concur_int_1000',
  'concur_int_500',
  'Npop', 'NPop',
  'time_lag', 
  'n_prev_tx',
  'pre_index_500m',
  'pre_index_1km',
  'pre_index_2km',
  'pre_index_3km',
  'min_dist_hf'
)

covarname_long <- process_covariates(df_long_process, dropnames)
covarname_short <- process_covariates(df_short_process, dropnames)

assert_that(setequal(covarname_long, covarname_short), 
            msg = "covariates differ between long and short data")

# check collinearity after dropping variables
check_cor(data = df_long_process, covarname_long)

plot_cor <- ggcorr(df_long_process[,covarname_long],
                   method = c("everything", "pearson"))


# dropping delta_response_time, delta_dist_index_quant because
# <5% prevalence in the individual data 
prop.table(table(df_long_process$delta_response_time))
prop.table(table(df_long_process$delta_dist_index_quant))

covarname_long = covarname_long[-which(covarname_long %in% c("delta_response_time", 
                                                             "delta_dist_index_quant",
                                                             "min_dist_hf"))]

# save covariates
saveRDS(covarname_long, namibia_inc_covar_path) 

# format data for analysis ---------------------------------------

data_human_long <- prep_htmle_data(data = df_long_process, reservoir = "human")
data_mosq_long <- prep_htmle_data(data = df_long_process, reservoir = "mosquito")
data_hm_long <- prep_htmle_data(data = df_long_process, reservoir = "human & mosquito")

data_human_short <- prep_htmle_data(data = df_short_process, reservoir = "human")
data_mosq_short <- prep_htmle_data(data = df_short_process, reservoir = "mosquito")
data_hm_short <- prep_htmle_data(data = df_short_process, reservoir = "human & mosquito")

data_human_list <- list(data_human_long, data_human_short)
data_mosq_list <- list(data_mosq_long, data_mosq_short)
data_hm_list <- list(data_hm_long, data_hm_short)

names(data_human_list) <- c("long", "short")
names(data_mosq_list) <- c("long", "short")
names(data_hm_list) <- c("long", "short")

# save data for analysis ---------------------------------------
saveRDS(data_human_list, namibia_human_process_path)
saveRDS(data_mosq_list, namibia_mosq_process_path)
saveRDS(data_hm_list, namibia_hm_process_path)


