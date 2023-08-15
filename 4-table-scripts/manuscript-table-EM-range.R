################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Table of ranges for effect modifier variables
# Above/below median 
################################################ 
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

ea_level_raw <- read_df(read.csv, namibia_raw_ea_level_path)

# drop ea that had no index cases ------------------------------------------- 
df_long <- readRDS(namibia_df_long_path)
included_eas <- unique(df_long$ea)

ea_level_raw <- ea_level_raw %>% filter(ea_actual != "")
ea_level <- ea_level_raw %>% filter(ea_actual %in% included_eas)

df_all <- readRDS(namibia_df_long_path)

df_all$tx_cov_cohort = df_all$tx_cov_cohort*100

# process data ------------------------------------------- 
df_ea <- ea_level %>% 
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

# create above/below median indicators for other EM variables -------------
df_ea <- df_ea %>% 
  mutate(pre_incidence_abovemed = ifelse(
    pre_incidence > median(df_ea$pre_incidence, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_spray_cover_abovemed = ifelse(
    pre_spray_cover > median(df_ea$pre_spray_cover, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_rainfall_abovemed = ifelse(
    pre_rainfall > median(df_ea$pre_rainfall, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_evi_abovemed = ifelse(
    pre_evi > median(df_ea$pre_evi, na.rm=T), 1, 0
  )) %>% 
  mutate(surface_temp_abovemed = ifelse(
    surface_temp > median(df_ea$surface_temp, na.rm=T), 1, 0
  )) %>% 
  mutate(ea_elevation_abovemed = ifelse(
    ea_elevation > median(df_ea$ea_elevation, na.rm=T), 1, 0
  )) %>% 
  mutate(pre_spray_cover = pre_spray_cover*100) %>% 
  mutate(pre_rainfall = pre_rainfall*1000)


# make table -------------
make_row <- function(EMvar){
  
  if(EMvar == "pre_incidence") digits = 1
  if(EMvar == "pre_spray_cover") digits = 1
  if(EMvar == "pre_rainfall") digits = 1
  if(EMvar == "pre_evi") digits = 2
  if(EMvar == "surface_temp") digits = 1
  if(EMvar == "ea_elevation") digits = 0
  if(EMvar == "tx_cov_cohort") digits = 1
  
  EMvar_abovemed = sym(paste0(EMvar, "_abovemed"))
  
  if(EMvar!="tx_cov_cohort") data = df_ea 
  if(EMvar=="tx_cov_cohort") data = df_all
  row <- data %>% group_by(!!EMvar_abovemed) %>% 
    summarise(min=min(!!sym(EMvar)), max=max(!!sym(EMvar))) %>%  
    mutate(min = sprintf(paste0("%0.", digits, "f"), min),
           max = sprintf(paste0("%0.", digits, "f"), max)) %>% 
    pivot_wider(names_from = !!EMvar_abovemed, names_prefix = "test",
                values_from = c(min, max)) %>% 
    mutate(modifier = EMvar) %>% 
    dplyr::select(modifier, min_test0, max_test0, min_test1, max_test1) 
  
  return(row)
}

EMvar_list = list("pre_incidence", "pre_spray_cover", "surface_temp",
                  "pre_rainfall",   "pre_evi",  "ea_elevation",
                  "tx_cov_cohort")


table <- lapply(EMvar_list, make_row) %>% bind_rows()

write.csv(table, file = paste0(table_path, "table-EM-var-ranges.csv"),
          row.names = F)





