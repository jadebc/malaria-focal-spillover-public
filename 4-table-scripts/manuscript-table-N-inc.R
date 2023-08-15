################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Table - incidence primary analysis results

# rfMDA, RACD: Short observation period
# RAVC: Long observation period 
# rfMDA+RAVC vs. RACD: Long observation period
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
data_human_list <- readRDS(namibia_human_process_path)
data_mosq_list <- readRDS(namibia_mosq_process_path)
data_hm_list <- readRDS(namibia_hm_process_path)

# load GLM unadjusted output - individual level ---------------------------------------
## direct effects ---------------------------------------
unadj_direct_human_indiv = readRDS(paste0(results_path, "namibia_glm_inc_direct_human_unadj.RDS"))
unadj_direct_mosq_indiv = readRDS(paste0(results_path, "namibia_glm_inc_direct_mosq_unadj.RDS"))
unadj_direct_hm_indiv = readRDS(paste0(results_path, "namibia_glm_inc_direct_hm_unadj.RDS"))

## spillover effects ---------------------------------------
unadj_spillover_human_indiv = readRDS(paste0(results_path, "namibia_glm_inc_spillover_human_unadj.RDS"))
unadj_spillover_mosq_indiv = readRDS(paste0(results_path, "namibia_glm_inc_spillover_mosq_unadj.RDS"))
unadj_spillover_hm_indiv = readRDS(paste0(results_path, "namibia_glm_inc_spillover_hm_unadj.RDS"))

## total effects ---------------------------------------
unadj_total_human_indiv = readRDS(paste0(results_path, "namibia_glm_inc_total_human_unadj.RDS"))
unadj_total_mosq_indiv = readRDS(paste0(results_path, "namibia_glm_inc_total_mosq_unadj.RDS"))
unadj_total_hm_indiv = readRDS(paste0(results_path, "namibia_glm_inc_total_hm_unadj.RDS"))


# load TMLE output ---------------------------------------
## direct effects ---------------------------------------
res_direct_human_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_human_cohort.RDS"))
res_direct_mosq_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_mosq_cohort.RDS"))
res_direct_hm_long_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_hm_long_cohort.RDS"))
res_direct_hm_cohort <- list(Long = res_direct_hm_long_cohort,
                            Short = NA)

## spillover effects ---------------------------------------
res_spillover_human_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_spillover_human_cohort.RDS"))
res_spillover_mosq_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_spillover_mosq_cohort.RDS"))
res_spillover_hm_indiv <- readRDS(paste0(results_path, "namibia_htmle_inc_spillover_hm_indiv.RDS"))

## total effects ---------------------------------------
res_total_human_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_total_human_cohort.RDS"))
res_total_mosq_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_total_mosq_cohort.RDS"))
res_total_hm_cohort <- readRDS(paste0(results_path, "namibia_htmle_inc_total_hm_cohort.RDS"))


# make table ---------------------------------------
make_sub_table <- function(data_short, data_long, reservoir_name, parameter_name, obs, data_level){
  print(paste0(reservoir_name, " - ", parameter_name))
  # individual-level data
  if(parameter_name == "Direct effect" & reservoir_name == "Human" & data_level == "indiv"){
    unadj_model = unadj_direct_human_indiv$short
    model = res_direct_human_indiv$Short
  } 
  if(parameter_name == "Direct effect" & reservoir_name == "Mosquito" & data_level == "indiv"){
    unadj_model = unadj_direct_mosq_indiv$long
    model = res_direct_mosq_indiv$Long
  } 
  if(parameter_name == "Direct effect" & reservoir_name == "Human & mosquito" & data_level == "indiv"){
    unadj_model = unadj_direct_hm_indiv$long
    model = res_direct_hm_indiv$Long
  } 

  if(parameter_name == "Spillover effect" & reservoir_name == "Human" & data_level == "indiv"){
    unadj_model = unadj_spillover_human_indiv$short
    model = res_spillover_human_indiv$Short
  } 
  if(parameter_name == "Spillover effect" & reservoir_name == "Mosquito" & data_level == "indiv"){
    unadj_model = unadj_spillover_mosq_indiv$long
    model = res_spillover_mosq_indiv$Long
  } 
  if(parameter_name == "Spillover effect" & reservoir_name == "Human & mosquito" & data_level == "indiv"){
    unadj_model = unadj_spillover_hm_indiv$long
    model = res_spillover_hm_indiv$Long
  } 
  
  if(parameter_name == "Total effect" & reservoir_name == "Human" & data_level == "indiv"){
    unadj_model = unadj_total_human_indiv$short
    model = res_total_human_indiv$Short
  } 
  if(parameter_name == "Total effect" & reservoir_name == "Mosquito" & data_level == "indiv"){
    unadj_model = unadj_total_mosq_indiv$long
    model = res_total_mosq_indiv$Long
  } 
  if(parameter_name == "Total effect" & reservoir_name == "Human & mosquito" & data_level == "indiv"){
    unadj_model = unadj_total_hm_indiv$long
    model = res_total_hm_indiv$Long
  } 
  
  # cohort-level data
  if(parameter_name == "Direct effect" & reservoir_name == "Human" & data_level == "cohort"){
    unadj_model = unadj_direct_human_indiv$short
    model = res_direct_human_cohort$Short
  } 
  if(parameter_name == "Direct effect" & reservoir_name == "Mosquito" & data_level == "cohort"){
    unadj_model = unadj_direct_mosq_indiv$long
    model = res_direct_mosq_cohort$Long
  } 
  if(parameter_name == "Direct effect" & reservoir_name == "Human & mosquito" & data_level == "cohort"){
    unadj_model = unadj_direct_hm_indiv$long
    model = res_direct_hm_cohort$Long
  } 
  
  if(parameter_name == "Spillover effect" & reservoir_name == "Human" & data_level == "cohort"){
    unadj_model = unadj_spillover_human_indiv$short
    model = res_spillover_human_cohort$Short
  } 
  if(parameter_name == "Spillover effect" & reservoir_name == "Mosquito" & data_level == "cohort"){
    unadj_model = unadj_spillover_mosq_indiv$long
    model = res_spillover_mosq_cohort$Long
  } 
  if(parameter_name == "Spillover effect" & reservoir_name == "Human & mosquito" & data_level == "cohort"){
    unadj_model = unadj_spillover_hm_indiv$long
    model = res_spillover_hm_cohort$Long
  } 
  
  if(parameter_name == "Total effect" & reservoir_name == "Human" & data_level == "cohort"){
    unadj_model = unadj_total_human_indiv$short
    model = res_total_human_cohort$Short
  } 
  if(parameter_name == "Total effect" & reservoir_name == "Mosquito" & data_level == "cohort"){
    unadj_model = unadj_total_mosq_indiv$long
    model = res_total_mosq_cohort$Long
  } 
  if(parameter_name == "Total effect" & reservoir_name == "Human & mosquito" & data_level == "cohort"){
    unadj_model = unadj_total_hm_indiv$long
    model = res_total_hm_cohort$Long
  } 
  
  if(length(model)>1){
      model_fit = model$res_est$estimates %>% filter(Psi_type == "RR") %>% 
        dplyr::select(Risk1, Risk0, Psi_hat, CI_l_unadj, CI_u_unadj, CI_l, CI_u)
  } 
  
  if(parameter_name == "Total effect"){
    data_short_group <- data_short %>% group_by(A) 
    data_long_group <- data_long %>% group_by(A) 
  }
  if(parameter_name == "Direct effect"){
    data_short_group <- data_short %>% filter(target_area == 1, int_recip == 1) %>% 
      group_by(A) 
    data_long_group <- data_long %>% filter(target_area == 1, int_recip == 1) %>% 
      group_by(A) 
  }
  if(parameter_name == "Spillover effect"){
    data_short_group <- data_short %>% 
      filter(target_area == 0 | (target_area == 1 & int_recip == 0)) %>% 
      group_by(A) 
    data_long_group <- data_long %>% 
      filter(target_area == 0 | (target_area == 1 & int_recip == 0)) %>% 
      group_by(A) 
  }
  
  if(reservoir_name=="Human") data = data_short_group
  if(reservoir_name!="Human") data = data_long_group
  
  # calculate individual N -------------------------
  n_table <- data %>% 
    summarise(n_cases = sum(Y),
              N = n()) %>% 
    mutate(inc = n_cases/N*1000) 
  
  n_table_wide <- cbind(
    n_table %>% filter(A==1) %>% dplyr::select(n_cases, N, inc),
    n_table %>% filter(A==0) %>% dplyr::select(n_cases, N, inc)
  )
  
  colnames(n_table_wide) = c(
    "n_case_tx", "N_tx", "inc_tx",
    "n_case_cont", "N_cont", "inc_cont"
  ) 
  
  n_table_wide <- n_table_wide %>% 
    mutate(reservoir = reservoir_name,
           parameter = parameter_name) %>% 
    mutate(IRR_crude = inc_tx/ inc_cont) 

  
  # calculate cohort N -------------------------
  cohort_table <- data %>% 
    summarise(n_cohorts = length(unique(id)))
  
  cohort_data <- data %>%
    group_by(A) %>%
    summarise(ncase = sum(Y),
              N = n()) %>%
    mutate(mean_inc = (ncase/N)  *1000)
  
  n_table_wide_cohort <- as.data.frame(cbind(
    cohorts_tx = cohort_table %>% filter(A==1) %>% pull(n_cohorts),
    cohorts_cont = cohort_table %>% filter(A==0) %>% pull(n_cohorts),
    inc_tx = cohort_data %>% filter(A==1) %>% pull(mean_inc),
    inc_cont = cohort_data %>% filter(A==0) %>% pull(mean_inc))) %>% 
    mutate(reservoir = reservoir_name,
           parameter = parameter_name) 
  
  # individual-level incidence ratio -------------------------
  if(data_level=="indiv"){
    if(reservoir_name == "Human") data <- data_short_group
    if(reservoir_name != "Human") data <- data_long_group

   # add model estimates
    if(!is.null(model)){
      n_table_out <- cbind(n_table_wide, model_fit, 
                           n_table_wide_cohort %>% dplyr::select(cohorts_tx, cohorts_cont)) %>% 
        rename(IRR_adj = Psi_hat) %>% 
        dplyr::select(reservoir, parameter, 
                      n_case_tx, N_tx, 
                      cohorts_tx, cohorts_cont,
                      inc_tx, 
                      n_case_cont, N_cont,  inc_cont, 
                      IRR_crude, IRR_adj, CI_l_unadj, CI_u_unadj, 
                      CI_l, CI_u) %>% 
        mutate(
          N_cohort = cohorts_tx + cohorts_cont,
               inc_tx = sprintf("%0.01f", inc_tx),
               inc_cont = sprintf("%0.01f", inc_cont),
               IRR_crude = sprintf("%0.02f", IRR_crude),
               IRR_adj = sprintf("%0.02f", IRR_adj),
               CI_l_unadj = sprintf("%0.02f", CI_l_unadj),
               CI_u_unadj = sprintf("%0.02f", CI_u_unadj),
               CI_l = sprintf("%0.02f", CI_l),
               CI_u = sprintf("%0.02f", CI_u),
               N= formatC(N_tx + N_cont, big.mark = ",")) %>% 
        mutate(
          result_IRR = paste0(IRR_adj, " (", CI_l_unadj, ", ",CI_u_unadj, ")"),
          result_IRR_adj = paste0(IRR_adj, " (", CI_l, ", ",CI_u, ")")
        )   %>% mutate(
        lb_crude = sprintf("%0.02f", unadj_model$lb),
        ub_crude = sprintf("%0.02f", unadj_model$ub)) %>% 
        mutate(result_IRR_crude = paste0(IRR_crude, " (", lb_crude, ", ", ub_crude, ")")) %>% 
        dplyr::select(
          reservoir, parameter, N_cohort, N, inc_tx, inc_cont, result_IRR_crude, result_IRR, result_IRR_adj)
    }else{
      n_table_out <- n_table_wide %>% mutate(
        inc_tx = sprintf("%0.01f", inc_tx),
        inc_cont = sprintf("%0.01f", inc_cont),
        IRR_crude = sprintf("%0.02f", IRR_crude),
        IRR_adj = NA,
        CI_l_unadj = NA,
        CI_u_unadj = NA,
        CI_l = NA,
        CI_u = NA,
        N= formatC(N_tx + N_cont, big.mark=",")
      ) %>%  mutate(
          lb_crude = sprintf("%0.02f", unadj_model$lb),
          ub_crude = sprintf("%0.02f", unadj_model$ub)) %>% 
        mutate(result_IRR_crude = paste0(IRR_crude, " (", lb_crude, ", ", ub_crude, ")")) %>% 
        mutate(
         result_IRR = NA,
         result_IRR_adj = NA) %>%  
        dplyr::select(
          reservoir, parameter, N, inc_tx, inc_cont, result_IRR_crude, 
          result_IRR, result_IRR_adj)
      
    }
  }
  
  
  
  # cohort-specific incidence ratio  -------------------------
  if(data_level=="cohort"){
    # N per cohort
    if(reservoir_name == "Human") data <- data_short_group
    if(reservoir_name != "Human") data <- data_long_group

    # add model estimates 
    if(!is.null(model)){
      n_table_out <- cbind(n_table_wide_cohort, model_fit, 
                           n_table_wide %>% dplyr::select(N_tx, N_cont)) %>% 
        rename(IRR_adj = Psi_hat) %>% 
        dplyr::select(reservoir, parameter, 
                      cohorts_tx,N_tx,N_cont,
                      inc_tx,
                      cohorts_cont,inc_cont,
                      IRR_adj, CI_l_unadj, CI_u_unadj, 
                      CI_l, CI_u) %>% 
        mutate(IRR_adj = sprintf("%0.02f", IRR_adj),
               CI_l_unadj = sprintf("%0.02f", CI_l_unadj),
               CI_u_unadj = sprintf("%0.02f", CI_u_unadj),
               CI_l = sprintf("%0.02f", CI_l),
               CI_u = sprintf("%0.02f", CI_u)) %>% 
        mutate(
          result_IRR = paste0(IRR_adj, " (", CI_l_unadj, ", ",CI_u_unadj, ")"),
          result_IRR_adj = paste0(IRR_adj, " (", CI_l, ", ",CI_u, ")"))  %>% 
        mutate(
          N_cohort = cohorts_tx + cohorts_cont,
          N= formatC(N_tx + N_cont, big.mark=","),
          inc_tx = sprintf("%0.01f", inc_tx),
          inc_cont = sprintf("%0.01f", inc_cont),
          # cohort_IRR_crude = unadj_model$RR_format$ptest,
          cohort_IRR_crude = sprintf("%0.02f",unadj_model$ptest),
          # lb_crude = unadj_model$RR_format$lb,
          lb_crude = sprintf("%0.02f",unadj_model$lb),
          # ub_crude = unadj_model$RR_format$ub,
          ub_crude = sprintf("%0.02f",unadj_model$ub)) %>% 
        mutate(
          result_IRR_crude = paste0(cohort_IRR_crude, " (", lb_crude, ", ", ub_crude, ")")) %>% 
        dplyr::select(
          reservoir, parameter, N_cohort, N, inc_tx, inc_cont, 
          result_IRR_crude, result_IRR, result_IRR_adj)
    }else{
      n_table_out <- n_table_wide %>% mutate(
        N_cohort = cohorts_tx + cohorts_cont,
        inc_tx = sprintf("%0.01f", inc_tx),
        inc_cont = sprintf("%0.01f", inc_cont),
        # cohort_IRR_crude = unadj_model$RR_format$ptest,
        cohort_IRR_crude = sprintf("%0.02f",unadj_model$ptest),
        # lb_crude = unadj_model$RR_format$lb,
        lb_crude = sprintf("%0.02f",unadj_model$lb),
        # ub_crude = unadj_model$RR_format$ub,
        ub_crude = sprintf("%0.02f",unadj_model$ub),
        result_IRR_crude = NA,
        result_IRR = NA,
        result_IRR_adj = NA
      ) %>%  mutate(
        result_IRR_crude = paste0(cohort_IRR_crude, " (", lb_crude, ", ", ub_crude, ")")) %>% 
        dplyr::select(
          reservoir, parameter, N_cohort, inc_tx, inc_cont, 
          result_IRR_crude, result_IRR, result_IRR_adj)
    }
   
      
  }
  
  return(n_table_out)
}

parameter_list <- list("Direct effect", "Spillover effect", "Total effect")


## make table (individual data) ---------------------------------------
hm_indiv_table <-  make_sub_table(
  data_short = data_hm_list$short,
  data_long = data_hm_list$long,
  reservoir_name = "Human & mosquito",
  parameter_name = "Spillover effect",
  data_level = "indiv"
)

## make table (cohort data) ---------------------------------------
human_cohort_table <- map_dfr(parameter_list,
                                    function(x)
                                      make_sub_table(
                                        data_short = data_human_list$short,
                                        data_long = data_human_list$long,
                                        reservoir_name = "Human",
                                        parameter_name = x,
                                        data_level = "cohort"
                                      ))
mosq_cohort_table <- map_dfr(parameter_list,
                                   function(x)
                                     make_sub_table(
                                       data_short = data_mosq_list$short,
                                       data_long = data_mosq_list$long,
                                       reservoir_name = "Mosquito",
                                       parameter_name = x,
                                       data_level = "cohort"
                                     ))
hm_cohort_table <- map_dfr(list("Direct effect", "Total effect"),
                                 function(x)
                                   make_sub_table(
                                     data_short = data_hm_list$short,
                                     data_long = data_hm_list$long,
                                     reservoir_name = "Human & mosquito",
                                     parameter_name = x,
                                     data_level = "cohort"
                                   ))

table <- bind_rows(
  human_cohort_table, 
  mosq_cohort_table, 
  hm_cohort_table,
  hm_indiv_table
) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito"))) %>% 
  arrange(reservoir, parameter)

# save tables ---------------------------------------
write.csv(table, paste0(table_path, "table-inc-primary.csv"))




