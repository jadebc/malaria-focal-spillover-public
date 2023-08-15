################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Cost-effectiveness analysis 

# N's and costs from Ntuku et al., 2022

# rfMDA, RACD: Short observation period
# RAVC: Long observation period 
# rfMDA+RAVC vs. RACD: Long observation period
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# get population proportion: prevalence --------------------------------------
## target area by arm   --------------------------------------
xs_clean <- readRDS(namibia_analysis_prev) %>% 
  mutate(human = ifelse(arm=="TO" |arm=="TV", "tx", "control"),
         mosq = ifelse(arm=="TV" | arm=="RV", "tx", "control"))

# human intervention 
pop_prop_human_de <- xs_clean %>% 
  group_by(human) %>% 
  summarise(pop_prop_de = mean(target_area)) %>% 
  select(human, pop_prop_de) %>% mutate("res" = "human") %>% rename("arm" = "human")

# mosquito intervention 
pop_prop_mosq_de <- xs_clean %>% 
  group_by(mosq) %>% 
  summarise(pop_prop_de = mean(target_area)) %>% 
  select(mosq, pop_prop_de) %>% mutate("res" = "mosq") %>% rename("arm" = "mosq")

# combined intervention 
pop_prop_both_de <- xs_clean %>% filter(arm %in% c("RO", "TV")) %>% 
  group_by(arm) %>%
  summarise(pop_prop_de = mean(target_area)) %>% 
  mutate(int_arm = ifelse(arm=="TV", "tx", "control")) %>% 
  select(int_arm, pop_prop_de) %>% mutate("res" = "both") %>% rename("arm" = "int_arm")

pop_prop_de <- rbind(pop_prop_human_de, pop_prop_mosq_de, pop_prop_both_de)

## spillover zone by arm  --------------------------------------

# human intervention 
pop_prop_human_sp <- xs_clean %>% 
  group_by(human) %>% 
  summarise(pop_prop_sp = mean(spall)) %>% 
  select(human, pop_prop_sp) %>% mutate("res" = "human") %>% rename("arm" = "human")

# mosquito intervention 
pop_prop_mosq_sp <- xs_clean %>% 
  group_by(mosq) %>% 
  summarise(pop_prop_sp = mean(spall)) %>% 
  select(mosq, pop_prop_sp) %>% mutate("res" = "mosq") %>% rename("arm" = "mosq")

# combined intervention 
pop_prop_both_sp <- xs_clean %>% filter(arm %in% c("RO", "TV")) %>% 
  group_by(arm) %>% summarise(pop_prop_sp = mean(spall)) %>% 
  mutate(int_arm = ifelse(arm=="TV", "tx", "control")) %>%
  select(int_arm, pop_prop_sp) %>% mutate("res" = "both") %>% rename("arm" = "int_arm")

pop_prop_sp <- rbind(pop_prop_human_sp, pop_prop_mosq_sp, pop_prop_both_sp)

pop_prop_df <- left_join(pop_prop_de, pop_prop_sp)

# read in spillover results for prevalence --------------------------------------
prev <- readRDS(paste0(results_path, "prevalence-tmle-results.RDS"))

# Make table, calculate cost effectiveness --------------------------------------------------------------
## TMLE results prevalence table, by arm and treatment -------------------------
write_cost_effectiveness_table <- function(tmle_bound) {
  if(tmle_bound == "upper") {
    prev_est = "upper_transformed"
  }
  
  if(tmle_bound == "main") {
    prev_est = "tmle_est"
  }
  
  if(tmle_bound == "lower") {
    prev_est = "lower_transformed"
  }
  
  ### human intervention arms
  prev_human_de <- prev$res_human_de$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_de" = prev_est)
  prev_human_sp <- prev$res_human_sp$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_sp" = prev_est)
  prev_human <- left_join(prev_human_de, prev_human_sp, by = c("tx" = "tx")) %>% mutate(res = "human")
  
  ### mosqutio intervention arms 
  prev_mosq_de <- prev$res_mosq_de$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_de" = prev_est)
  prev_mosq_sp <- prev$res_mosq_sp$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_sp" = prev_est)
  prev_mosq <- left_join(prev_mosq_de, prev_mosq_sp, by = c("tx" = "tx")) %>% mutate("res" = "mosq")
  
  # combined intervention arms
  prev_both_de <- prev$res_both_de$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_de" = prev_est)
  prev_both_sp <- prev$res_both_sp$qPCRposneg$res_est %>% filter(type=="TSM" & model=="Unadjusted") %>% 
    dplyr::select(param, !!sym(prev_est)) %>% rename("tx" = "param", "pop_prev_sp" = prev_est)
  prev_both <- left_join(prev_both_de, prev_both_sp, by = c("tx" = "tx")) %>% mutate("res" = "both")
  
  all_prev <- rbind(prev_human, prev_mosq, prev_both) %>% 
    mutate(arm = ifelse(tx=="E[Y_{A=0}]", "control", "tx")) %>% select(res, arm, pop_prev_de, pop_prev_sp)
  
  ## perform calculations
  table <- left_join(all_prev, pop_prop_df) %>% 
    mutate("N" = c(9406.1757719715, 9389.83050847458, 9396.75174013921, 9405.59440559441, 4533.333333, 4512.195122),
           "total_cost" = c(354750, 368321, 261409, 461661, 127312, 234223)) # hard-coded costs from Ntuku et al. 2022
  
  ### number of individuals in each zone (target & spillover)
  table <- table %>% mutate("N_de" = N*pop_prop_de, "N_sp" = N*pop_prop_sp)
  
  ### number of cases in each zone (target & spillover)
  table <- table %>% mutate("prev_cases_de" = N_de*pop_prev_de, "prev_cases_sp" = N_sp*pop_prev_sp)
  
  ### number of cases averted in each zone (target & spillover)
  cases_averted_de <- table %>% group_by(res) %>% select(arm, prev_cases_de) %>% 
    pivot_wider(names_from = c(arm), values_from = prev_cases_de) %>% 
    mutate("cases_averted_de" = (tx - control)*-1) %>% select(res, cases_averted_de)
  
  cases_averted_sp <- table %>% group_by(res) %>% select(arm, prev_cases_sp) %>% 
    pivot_wider(names_from = c(arm), values_from = prev_cases_sp) %>% 
    mutate("cases_averted_sp" = (tx - control)*-1) %>% select(res, cases_averted_sp)
  
  ### cost difference compared to Ntuku et al. 2022
  cost_difference <- table %>% group_by(res) %>% select(arm, total_cost) %>% 
    pivot_wider(names_from = c(arm), values_from = total_cost) %>% 
    mutate("cost_difference" = (tx - control)) %>% select(res, cost_difference) %>% mutate("arm" = "tx")
  
  #### hard-coded costs of interventions (human, mosquito, combined -- respectively)
  cost_difference <- cost_difference %>% ungroup() %>% mutate("icer_ntuku" = c(162, 2670, 1812))
  
  cases_averted = left_join(cases_averted_de, cases_averted_sp) %>% 
    mutate(arm = "tx") %>% group_by(res) %>% mutate("total_cases_averted" = cases_averted_de + cases_averted_sp)
  
  ### merge together prevalence & cost calculations with cases averted & calculate costs
  table <- left_join(table, cases_averted)
  
  table <- left_join(table, cost_difference) %>% group_by(res) %>% 
    mutate("icer_prev_de_sp" = cost_difference/total_cases_averted) %>% ungroup() %>% 
    mutate("icer_perc_diff_from_ntuku" = ((icer_prev_de_sp - icer_ntuku)/icer_ntuku)*100)
  
  ### add in complete intervention arm names
  table_human <- table %>% filter(res=="human") %>% 
    mutate(res = "Human", int_arm = ifelse(arm=="tx", "rfMDA", "RACD"))
  
  table_mosq <- table %>% filter(res=="mosq") %>% 
    mutate(res = "Mosquito", int_arm = ifelse(arm=="tx", "RAVC", "No RAVC"))
  
  table_both <- table %>% filter(res=="both") %>% 
    mutate(res = "Human + Mosquito", int_arm = ifelse(arm=="tx", "rfMDA + RAVC", "RACD"))
  
  table <- rbind(table_human, table_mosq, table_both)
  
  saveRDS(table, file = paste0("cost_effectiveness_analysis-", tmle_bound, ".RDS"))
}

tmle_bounds <- c("upper", "main", "lower")
lapply(tmle_bounds, write_cost_effectiveness_table)

# Formal table formatting & cleaning -------------------------------------------------------------------
upper_table <- readRDS(paste0("cost_effectiveness_analysis-", "upper", ".RDS"))
main_table <- readRDS(paste0("cost_effectiveness_analysis-", "main", ".RDS"))
lower_table <- readRDS(paste0("cost_effectiveness_analysis-", "lower", ".RDS"))

round_var_whole <- function(var_to_round, df) {
  df[var_to_round] = round(df[var_to_round], digits = 0)
  return(df[var_to_round])
}

round_var_3sf <- function(var_to_round, df) {
  df[var_to_round] = round(df[var_to_round], digits = 3)
  return(df[var_to_round])
}

whole_number_var_list <- c("N", "N_de", "N_sp", "prev_cases_de", "cases_averted_de", "prev_cases_sp", 
                           "cases_averted_sp", "total_cases_averted", "icer_perc_diff_from_ntuku")

decimal_var_list <- c("pop_prop_de", "pop_prop_sp", "pop_prev_de", "pop_prev_sp")

clean_table <- function(df) {
  df_round <- map_dfc(.x = whole_number_var_list, df = df, .f = round_var_whole)
  
  df_3sf <- map_dfc(.x = decimal_var_list, df = df, .f = round_var_3sf)
  
  df_total_cost <- df %>% mutate(total_cost = paste0("$", formatC(total_cost, format="f", digits=0, big.mark=","))) %>% 
    select(res, int_arm, total_cost)
  
  df_other_costs <- df %>% filter(!is.na(cost_difference)) %>% 
    mutate(cost_difference = paste0("$", formatC(cost_difference, format="f", digits=0, big.mark=",")),
           icer_ntuku = paste0("$", formatC(icer_ntuku, format="f", digits=0, big.mark=",")),
           icer_prev_de_sp = paste0("$", formatC(icer_prev_de_sp, format="f", digits=0, big.mark=","))) %>% 
    select(res, int_arm, cost_difference, icer_ntuku, icer_prev_de_sp)
  
  df_costs <- left_join(df_total_cost, df_other_costs)
  
  df <- cbind(df_round, df_3sf) %>% cbind(df_costs)
  
  df <- df %>% select(res, int_arm, total_cost, cost_difference, icer_ntuku, N, pop_prop_de, pop_prop_sp, N_de, N_sp,
                      pop_prev_de, prev_cases_de, cases_averted_de, pop_prev_sp, prev_cases_sp, cases_averted_sp,
                      total_cases_averted, icer_prev_de_sp, icer_perc_diff_from_ntuku)
  
  return(df)
}

main_table <- clean_table(main_table)
upper_table <- clean_table(upper_table) %>% filter(!is.na(total_cases_averted)) # only need cost CIs
lower_table <- clean_table(lower_table) %>% filter(!is.na(total_cases_averted)) # only need cost CIs

## add CIs to main table
main_table_ci <- main_table %>% filter(!is.na(total_cases_averted)) %>% 
  mutate("total_cases_averted_ci" = paste0(total_cases_averted, " (", lower_table$total_cases_averted, ", ",
                                           upper_table$total_cases_averted, ")"),
         "icer_prev_de_sp_ci" = paste0(icer_prev_de_sp, " (", lower_table$icer_prev_de_sp, ", ",
                                       upper_table$icer_prev_de_sp, ")")) %>% 
  select(res, int_arm, total_cases_averted_ci, icer_prev_de_sp_ci)

main_table <- left_join(main_table, main_table_ci)

## add in intervention rows, rm unnecessary/repeat columns
main_table <- main_table %>% select(int_arm, total_cost, N_de, N_sp, pop_prev_de, pop_prev_sp, prev_cases_de, 
                                    prev_cases_sp, total_cases_averted_ci, icer_prev_de_sp_ci, icer_perc_diff_from_ntuku)

### create dummy rows for aesthetics/labelling purposes
human_row <- as.data.frame(t(c("Human intervention", rep(" ", 10))))
mosq_row <- as.data.frame(t(c("Mosquito intervention", rep(" ", 10))))
both_row <- as.data.frame(t(c("Combined intervention", rep(" ", 10))))

colnames(human_row) = colnames(main_table)
colnames(mosq_row) = colnames(main_table)
colnames(both_row) = colnames(main_table)

main_table <- rbind(human_row, main_table[1:2,], mosq_row, main_table[3:4,], both_row, main_table[5:6,])
main_table[is.na(main_table)] <- "(ref)"

main_table <- main_table %>% add_row(.before = 1, int_arm = " ", total_cost = " ", N_de = "N Individuals", N_sp = " ", 
                                     pop_prev_de = "Prevalence", pop_prev_sp = " ", prev_cases_de = "Prevalent Cases", 
                                     prev_cases_sp = " ", total_cases_averted_ci = "Total prevalent cases averted (95% CI)",
                                     icer_prev_de_sp_ci = "Incremental cost-effectiveness ratio (95% CI)", 
                                     icer_perc_diff_from_ntuku = "% change from original estimate")

main_table <- main_table %>% add_row(.before = 2, int_arm = " ", total_cost = "Intervention Cost", N_de = "Target area",
                                     N_sp = "Spillover zone", pop_prev_de = "Target area", pop_prev_sp = "Spillover zone",
                                     prev_cases_de = "Target area", prev_cases_sp = "Spillover zone", 
                                     total_cases_averted_ci = " ", icer_prev_de_sp_ci = " ", 
                                     icer_perc_diff_from_ntuku = " ")

write.csv(main_table, file = paste0(table_path, "cost-effectiveness-table.csv"))
