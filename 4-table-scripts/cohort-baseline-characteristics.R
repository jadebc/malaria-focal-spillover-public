################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Table - Cohort characteristics at baseline
################################################ 
rm(list=ls())
library(tidyverse)

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
df_long <- readRDS(namibia_df_long_path)
df_short <- readRDS(namibia_df_short_path)

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

long <- df_long %>% filter(!cohort_id %in% long_drop_ids)
short <- df_short %>% filter(!cohort_id %in% short_drop_ids)

# pre-process covariates
long$tx_cov_cohort = long$tx_cov_cohort*100
short$tx_cov_cohort = short$tx_cov_cohort*100

long$prop_prev_tx = long$prop_prev_tx*100
short$prop_prev_tx = short$prop_prev_tx*100

long$pre_rainfall = long$pre_rainfall*1000
short$pre_rainfall = short$pre_rainfall*1000

long$pre_spray_cover = long$pre_spray_cover*100
short$pre_spray_cover = short$pre_spray_cover*100

long$dist_hf = long$dist_hf/1000
short$dist_hf = short$dist_hf/1000

long_direct <- long %>% filter(target_area == 1 & int_recip == 1)
long_direct <- as.data.frame(long_direct)
long_spillover <- long %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
long_spillover <- as.data.frame(long_spillover)

short_direct <- short %>% filter(target_area == 1 & int_recip == 1)
short_direct <- as.data.frame(short_direct)
short_spillover <- short %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
short_spillover <- as.data.frame(short_spillover)

##############################################
##############################################
# Documentation: calc_outcome
# Usage: calc_outcome(data, resevoir, outcome, zone)
# Description: Make a vector of length 2 including the control and treatment 
# arm for the desired reservoir and outcome

# Args/Options:
# data: a dataframe with study characteristics variables
# resevoir: a string for the resevoir name (human, mosquito, both)
# outcome: a string specifying which outcome will be generated
# zone: a number (1 or 0) describing a spillover zone (0) or target area (1)

# Returns: a vector of length 2 with the calculated outcome for the treatment and control
# Output: prints the vector

calc_outcome <- function(reservoir, outcome, zone) {
  
  if(zone == 0) {
    data_short = short_spillover
    data_long = long_spillover
  }

  if(zone == 1) {
    data_short = short_direct
    data_long = long_direct
  }
  
  # set treatment & control based on reservoir
  if(reservoir == "human") {
    cell_data <- data_short  %>% 
      mutate(treatment_ind = ifelse(intarm == "TO"| intarm=="TV", 1, 0))
  }
  
  if(reservoir == "mosquito") {
    cell_data <- data_long %>% 
      mutate(treatment_ind = ifelse(intarm == "RV" | intarm == "TV", 1, 0))
  }
  
  if(reservoir == "both") {
    cell_data_human <- data_short %>% 
      filter(intarm == "RO") %>% 
      mutate(treatment_ind = 0)
    cell_data_mosquito <- data_long %>% 
      filter(intarm == "TV") %>% 
      mutate(treatment_ind = 1)
    cell_data = rbind(cell_data_human, cell_data_mosquito)
  }
  
  # spillover_filtered <- cell_data %>% filter(target_area == zone)
  
  # N cohorts 
  if (outcome == "n_cohorts") {
    cohorts <- cell_data %>% 
      group_by(treatment_ind) %>% 
      summarize(n_cohorts = length(unique(cohort_id)))
    
    return(cohorts$n_cohorts)
  }
  
  # Cohort size (not filtering by spillover zone)
  if (outcome == "cohort_size") {
    size_table <- cell_data %>% group_by(treatment_ind, cohort_id) %>% 
      summarize(each_cohort_size = length(indiv_id)) %>%
      summarize(mean_size = round(mean(each_cohort_size, na.rm = T)),
                se_size = round((sd(each_cohort_size, na.rm = T))/sqrt(length(each_cohort_size))))

    format_size_se <- paste("(", size_table$se_size, ")", sep = "")
    
    return(paste(size_table$mean_size, format_size_se, sep = " "))
  }
  
  # number of index cases 
  if (outcome == "indexcase") {
    indexcases <- cell_data %>% 
      group_by(treatment_ind) %>% 
      summarize(n_indexcases = sum(indexcase))
    
    return(indexcases$n_indexcases)
  }

  # calculate & format mean aand se for other parameters
  parameter_table <- data.frame(cell_data$treatment_ind, cell_data[[outcome]])  
  colnames(parameter_table) <- c("treatment_ind", "parameter")
    
  if(outcome %in% c("n_cohorts", "cohort_size",
                    "pop_size_ea", "pre_incidence","pre_spray_cover","dist_hf")){
    tab <- parameter_table %>%
      group_by(treatment_ind) %>%
      summarize(mean = sprintf("%2.1f", mean(parameter, na.rm = T)),
                se = (sd(parameter, na.rm = T))/sqrt(length(parameter)))
    formatted_se <- paste("(", sprintf("%2.2f", tab$se), ")", sep = "")
    return(paste(tab$mean, formatted_se, sep = " "))
  }else{
    if(outcome %in% c("pre_rainfall", "surface_temp")){
      tab <- parameter_table %>%
        group_by(treatment_ind) %>%
        summarize(med = sprintf("%2.1f", median(parameter, na.rm = T)),
                  min = sprintf("%2.1f", min(parameter, na.rm = T)),
                  max = sprintf("%2.1f", max(parameter, na.rm = T))) %>% 
        mutate(med_range = paste0(med, " (", min, ", ", max, ")")) %>% 
        dplyr::select(-c(min, max))
      return(pull(tab, med_range))
    }
    if(outcome == "pre_evi"){
      tab <- parameter_table %>%
        group_by(treatment_ind) %>%
        summarize(med = sprintf("%2.2f", median(parameter, na.rm = T)),
                  min = sprintf("%2.2f", min(parameter, na.rm = T)),
                  max = sprintf("%2.2f", max(parameter, na.rm = T))) %>% 
        mutate(med_range = paste0(med, " (", min, ", ", max, ")")) %>% 
        dplyr::select(-c(min, max))
      return(pull(tab, med_range))
    }
    if(outcome == "ea_elevation"){
      tab <- parameter_table %>%
        group_by(treatment_ind) %>%
        summarize(med = sprintf("%2.0f", median(parameter, na.rm = T)),
                  min = sprintf("%2.0f", min(parameter, na.rm = T)),
                  max = sprintf("%2.0f", max(parameter, na.rm = T))) %>% 
        mutate(med_range = paste0(med, " (", min, ", ", max, ")")) %>% 
        dplyr::select(-c(min, max))
      return(pull(tab, med_range))
    }

  }



}

# set up lists of parameters to perform next function on 
characteristics <- c("n_cohorts", "cohort_size",
                     "pop_size_ea", "pre_incidence", "pre_spray_cover", 
                     "dist_hf", "pre_rainfall", "pre_evi", "ea_elevation","surface_temp")


##############################################
##############################################
# Documentation: make_table
# Usage: make_table(data, zone)
# Description: Create a six-columned table of treatment and control outcomes
# for the three treatment reservoirs in the trial. Requires calc_outcome above and 
# a set of parameters from the dataframe passed into data to perform calc_outcome on. 

# Args/Options:
# data: a dataframe with study characteristics variables; the specific parameters
# to be included in the table must be set up as vectors that calc_outcome will call within
# the make_table function
# zone: a number (1 or 0) describing a spillover zone (0) or target area (1)

# Returns: a dataframe with 6 variable columns(treatment & control for each reservoir) + a names column. 
# The number of rows = the length of the vector passed in as the parameters + 3. 
# Output: Prints the dataframe.

make_table <- function(zone) {
  row <- list()
  
  human_table <- c("RACD", "rfMDA")
  for (i in 1:length(characteristics)) {
    row[[i]] <- calc_outcome("human", characteristics[i], zone)
    human_table <- rbind(human_table, row[[i]])
  }
  
  mosquito_table <- c("No RAVC", "RAVC")
  for (i in 1:length(characteristics)) {
    row[[i]] <- calc_outcome("mosquito", characteristics[i], zone)
    mosquito_table <- rbind(mosquito_table, row[[i]])
  }
  
  both_table <- c("RACD", "rfMDA + RAVC")
  for (i in 1:length(characteristics)) {
    row[[i]] <- calc_outcome("both", characteristics[i], zone)
    both_table <- rbind(both_table, row[[i]])
  }
  
  # Bind reservoir tables together into complete table
  complete_table <- cbind(human_table, mosquito_table, both_table)
  
  rownames(complete_table) <-   c("",characteristics)

  rownames(complete_table) <- c(" ", "Number of cohorts", 
                             "Mean cohort population size (SE)",
                             "Mean cluster population size (SE)",
                             "Malaria incidence per 1,000 in 2016",
                             "Pre-season indoor residual spray coverage 2016",
                             "Distance to nearest healthcare facility (km)",
                             "Median monthly rainfall November 2016-April 2017 (mm)",
                             "Median enhanced vegetative index January 2017-July 2017",
                             "Median elevation (m)",
                             "Median daytime land surface temperature (C)"
)
  as.data.frame(complete_table)
  complete_table
}

# Make spillover & target zone tables

direct_table <- write.csv(make_table(1), file = paste0(table_path, "cohort_characteristics_direct.csv"))

spillover_table <- write.csv(make_table(0), file = paste0(table_path, "cohort_characteristics_spillover.csv"))
