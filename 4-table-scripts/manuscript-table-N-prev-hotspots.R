################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Table - prevalence hot spot analysis results

# rfMDA, RACD: Short observation period
# RAVC: Long observation period 
# rfMDA+RAVC vs. RACD: Long observation period
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load data ---------------------------------------
# xs_clean <- readRDS(namibia_analysis_prev)
xs <- readRDS(paste0(namibia_process_path, "namibia_prev_hotspot.RDS"))

# load TMLE output ---------------------------------------
prev <- readRDS(paste0(results_path, "prevalence-hotspot-tmle-results.RDS"))


make_prev_table <- function(data, reservoir_name, parameter_name, yname){
  
  print(paste0(reservoir_name, " - ", parameter_name , " - ", yname))
  
  # Ns for prevalence analysis --------------------
  
  # Drop if outcome is missing
  data <- data %>% filter(!is.na(!!sym(yname)))
  
  ## Drop observations outside spillover zone ----------------------
  if(parameter_name=="Direct effect"){
    data_sub <- data %>% filter(target_area==1)
  }
  if(parameter_name=="Spillover effect"){
    data_sub <- data %>% filter(spall==1)
  }
  if(parameter_name=="Total effect"){
    data_sub <- data %>% filter(target_area==1 |spall==1)
  }
  
  # binarize treatment variable
  # human reservoir
  if (reservoir_name == "Human") {
    data_sub <- data_sub %>%
      mutate(             
        ## add a column for constructing interaction term
        inter = as.numeric(arm %in% c("RV", "TV")),
        arm_yn = case_when(arm %in% c("TO", "TV") ~ 1,
                           arm %in% c("RO", "RV") ~ 0)
      )
  }
  
  # mosquito reservoir
  if(reservoir_name == "Mosquito"){
    data_sub <- data_sub %>%
      mutate(
        ## add a column for constructing interaction term
        inter = as.numeric(arm %in% c("TO", "TV")),
        arm_yn = case_when(arm %in% c("RV", "TV") ~ 1,
                           arm %in% c("RO", "TO") ~ 0)
      )
  }
  
  # human & mosquito reservoir
  if(reservoir_name == "Human & mosquito"){
    data_sub <- data_sub %>% filter(arm %in% c("TV", "RO"))
    data_sub <- data_sub %>% mutate(arm_yn = case_when(arm %in% c("TV") ~ 1,
                                                       arm %in% c("RO") ~ 0))
    
  }
  
  N <- nrow(data_sub)
  N_int <- nrow(data_sub %>% filter(arm_yn==1))
  N_cont <- nrow(data_sub %>% filter(arm_yn==0))
  
  prev_arm <- data_sub %>% group_by(arm_yn) %>% 
    summarise(prev = mean(!!sym(yname), na.rm = T))
  if(reservoir_name == "Human"){
    prev_control = prev_arm %>% filter(arm_yn==0) %>% pull(prev)
    prev_int = prev_arm %>% filter(arm_yn==1) %>% pull(prev)
  } 
  if(reservoir_name == "Mosquito"){
    prev_control = prev_arm %>% filter(arm_yn==0) %>% pull(prev)
    prev_int = prev_arm %>% filter(arm_yn==1) %>% pull(prev)
  } 
  if(reservoir_name == "Human & mosquito"){
    prev_control = prev_arm %>% filter(arm_yn==0) %>% pull(prev)
    prev_int = prev_arm %>% filter(arm_yn==1) %>% pull(prev)
  } 
  
  # process model outputs -----------------------------------------
  if(reservoir_name=="Human"){
    if(parameter_name=="Direct effect"){
       if(yname == "cases_qpcr")  model <- prev$res_human_de$cases_qpcr$res_est
       if(yname == "cases_rdt")  model <- prev$res_human_de$cases_rdt$res_est
       if(yname == "cases_hsrdt")  model <- prev$res_human_de$cases_hsrdt$res_est
    }
    if(parameter_name=="Spillover effect"){
      if(yname == "cases_qpcr")  model <- prev$res_human_sp$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_human_sp$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_human_sp$cases_hsrdt$res_est
    }
    if(parameter_name=="Total effect"){
      if(yname == "cases_qpcr")  model <- prev$res_human_te$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_human_te$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_human_te$cases_hsrdt$res_est
    }
  }
  
  if(reservoir_name=="Mosquito"){
    if(parameter_name=="Direct effect"){
      if(yname == "cases_qpcr")  model <- prev$res_mosq_de$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_mosq_de$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_mosq_de$cases_hsrdt$res_est
    }
    if(parameter_name=="Spillover effect"){
      if(yname == "cases_qpcr")  model <- prev$res_mosq_sp$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_mosq_sp$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_mosq_sp$cases_hsrdt$res_est
    }
    if(parameter_name=="Total effect"){
      if(yname == "cases_qpcr")  model <- prev$res_mosq_te$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_mosq_te$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_mosq_te$cases_hsrdt$res_est
    }
  }
  
  
  if(reservoir_name=="Human & mosquito"){
    if(parameter_name=="Direct effect"){
      if(yname == "cases_qpcr")  model <- prev$res_both_de$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_both_de$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_both_de$cases_hsrdt$res_est
    }
    if(parameter_name=="Spillover effect"){
      if(yname == "cases_qpcr")  model <- prev$res_both_sp$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_both_sp$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_both_sp$cases_hsrdt$res_est
    }
    if(parameter_name=="Total effect"){
      if(yname == "cases_qpcr")  model <- prev$res_both_te$cases_qpcr$res_est
      if(yname == "cases_rdt")  model <- prev$res_both_te$cases_rdt$res_est
      if(yname == "cases_hsrdt")  model <- prev$res_both_te$cases_hsrdt$res_est
    }
  }
  
  
  if(all(!is.na(model))){
    estimate_unadj <- model %>% filter(type=="RR" & model=="Unadjusted") %>% 
      dplyr::select(psi_transformed, lower_transformed, upper_transformed) %>% 
      rename(PR_unadj = psi_transformed, 
             lb_unadj = lower_transformed,
             ub_unadj = upper_transformed)
    estimate_adj <- model %>% filter(type=="RR" & model=="Adjusted") %>% 
      dplyr::select(psi_transformed, lower_transformed, upper_transformed) %>% 
      rename(PR_adj = psi_transformed, 
             lb_adj = lower_transformed,
             ub_adj = upper_transformed)
    
    if(nrow(estimate_unadj)>0){
      if(estimate_unadj$PR_unadj>100){
        estimate_unadj$PR_unadj=NA
        estimate_unadj$lb_unadj=NA
        estimate_unadj$ub_unadj=NA
      }
    }
    
    if(nrow(estimate_adj)>0){
      if(estimate_adj$PR_adj>100){
        estimate_adj$PR_adj=NA
        estimate_adj$lb_adj=NA
        estimate_adj$ub_adj=NA
      }
    }
    
    out = cbind(N, prev_int, prev_control, estimate_unadj, estimate_adj)
    
 
    out = out %>% 
      mutate(prev_int = sprintf("%0.03f", prev_int),
             prev_control = sprintf("%0.03f", prev_control),
             PR_unadj_f = sprintf("%0.02f", PR_unadj),
             lb_unadj_f = sprintf("%0.02f", lb_unadj),
             ub_unadj_f = sprintf("%0.02f", ub_unadj),
             PR_adj_f = sprintf("%0.02f", PR_adj),
             lb_adj_f = sprintf("%0.02f", lb_adj),
             ub_adj_f = sprintf("%0.02f", ub_adj)) %>% 
      
      mutate(
        result_PR = paste0(PR_unadj_f, " (", lb_unadj_f, ", ",ub_unadj_f, ")"),
        result_PR_adj = paste0(PR_adj_f, " (", lb_adj_f, ", ",ub_adj_f, ")")
      ) %>% 
      mutate(reservoir = reservoir_name, 
             parameter = parameter_name,
             N_int = N_int,
             N_cont = N_cont) %>% 
      dplyr::select(reservoir, parameter, N_int, N_cont, prev_int, prev_control, 
                    PR_unadj, lb_unadj, ub_unadj,
                    PR_adj, lb_adj, ub_adj,
                    result_PR, result_PR_adj)
    
  }else{
    out <- data.table(
      reservoir = reservoir_name,
      parameter = parameter_name, 
      N_int = N_int,
      N_cont = N_cont,
      prev_int = sprintf("%0.03f", prev_int), 
      prev_control = sprintf("%0.03f", prev_control), 
      result_PR = NA, 
      result_PR_adj = NA
    )
  }
 
  return(out)
}



parameter_list <- list("Direct effect", "Spillover effect","Total effect")
reservoir_list <- list("Human", "Mosquito", "Human & mosquito")

## make table  ---------------------------------------
qpcr_list = list()
for(i in 1:length(parameter_list)){
  qpcr_list[[i]] = map_dfr(reservoir_list,
          function(x)
            make_prev_table(
              data = xs,
              reservoir_name = x,
              parameter_name = parameter_list[[i]],
              yname = "cases_qpcr"
            ))
}

qpcr_output = bind_rows(qpcr_list) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito"))) %>% 
  arrange(reservoir, parameter)

qpcr_table <- qpcr_output %>% dplyr::select(-c(PR_unadj, lb_unadj, ub_unadj,
                                               PR_adj, lb_adj, ub_adj))

rdt_list = list()
for(i in 1:length(parameter_list)){
  rdt_list[[i]] = map_dfr(reservoir_list,
                           function(x)
                             make_prev_table(
                               data = xs,
                               reservoir_name = x,
                               parameter_name = parameter_list[[i]],
                               yname = "cases_rdt"
                             ))
}

rdt_output = bind_rows(rdt_list) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito"))) %>% 
  arrange(reservoir, parameter)

rdt_table <- rdt_output %>% dplyr::select(-c(PR_unadj, lb_unadj, ub_unadj,
                                               PR_adj, lb_adj, ub_adj))


hsrdt_list = list()
for(i in 1:length(parameter_list)){
  hsrdt_list[[i]] = map_dfr(reservoir_list,
                           function(x)
                             make_prev_table(
                               data = xs,
                               reservoir_name = x,
                               parameter_name = parameter_list[[i]],
                               yname = "cases_hsrdt"
                             ))
}

hsrdt_output = bind_rows(hsrdt_list) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito"))) %>% 
  arrange(reservoir, parameter)

hsrdt_table <- hsrdt_output %>% dplyr::select(-c(PR_unadj, lb_unadj, ub_unadj,
                                             PR_adj, lb_adj, ub_adj))

# save data for plot
saveRDS(qpcr_output, paste0(namibia_process_path, "plot-data-prev-primary-qpcr-hotspot.RDS"))


# save tables ---------------------------------------
write.csv(qpcr_table, paste0(table_path, "table-prev-primary-qpcr-hotspot.csv"))
write.csv(rdt_table, paste0(table_path, "table-prev-primary-rdt-hotspot.csv"))
write.csv(hsrdt_table, paste0(table_path, "table-prev-primary-hsrdt-hotspot.csv"))




