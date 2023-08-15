################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# percentage of population 
#in target areas and spillover zones
################################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
df_long <- readRDS(namibia_df_long_path)

# mean percentage across cohorts
# of individuals in target area who were treated
df_long %>% filter(target_area==1) %>% 
  group_by(cohort_id) %>% 
  summarise(percent = mean(int_recip)) %>% 
  ungroup() %>% 
  summarise(mean_percent = mean(percent))

# mean percentage across cohorts
# of individuals within 1km who were not treated 
df_long %>% 
  group_by(cohort_id) %>% 
  summarise(percent = mean(int_recip)) %>% 
  ungroup() %>% 
  summarise(mean_percent = mean(percent)) %>% 
  mutate(mean_percent_untx = 1-mean_percent)



# load prevalence data ------
prev <- readRDS(namibia_analysis_prev) 

mean(prev$target_area)
mean(prev$spall)

nrow(prev[prev$target_area==1,])
nrow(prev[prev$spall==1,])


