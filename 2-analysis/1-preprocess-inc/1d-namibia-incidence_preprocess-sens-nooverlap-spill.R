################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Preprocess data for 
# Sensitivity analysis: no overlap between spillover zones
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ------
df_long <- readRDS(namibia_df_long_sens_nooverlap_spill_path)
df_short <- readRDS(namibia_df_short_sens_nooverlap_spill_path)

# preprocess data ------
df_long_process <- preprocess_data(df_long)
df_short_process <- preprocess_data(df_short)

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

# plot_cor <- ggcorr(df_long_process[,covarname_long], 
#                    method = c("everything", "pearson")) 

saveRDS(covarname_long, namibia_inc_covar_sens_nooverlap_spill_path)

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
saveRDS(data_human_list, namibia_human_process_sens_nooverlap_spill_path)
saveRDS(data_mosq_list, namibia_mosq_process_sens_nooverlap_spill_path)
saveRDS(data_hm_list, namibia_hm_process_sens_nooverlap_spill_path)


