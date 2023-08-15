################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Table - prevalence serology results

# rfMDA, RACD: Short observation period
# RAVC: Long observation period 
# rfMDA+RAVC vs. RACD: Long observation period
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(readxl)

# load and process analysis dataset ---------------------------------------------------
xs_clean <- readRDS(namibia_analysis_prev) %>% arrange(eaid, hhid, iid)
raw_data <- read_excel(namibia_raw_sero_path)

sero_cols <- raw_data %>% dplyr::select(eaid, hhid, sample, Etramp5.Ag1, Etramp5.Ag1_pos) %>% 
  rename(iid = sample) %>% 
  arrange(eaid, hhid, iid) %>% 
  distinct()

sero <- full_join(xs_clean, sero_cols, by = c("eaid", "hhid", "iid"))

# drop if no spillover zone information
sero <- sero %>% filter(!is.na(spall))

# drop if no serology
sero <- sero %>% filter(!is.na(Etramp5.Ag1))


# specify the variables used for identification, but not for estimation
dropme <- c("eaid", "ea", "xsid", "iid", "latitude", "longitude",
            "date", "spillover_zone")

sero <- sero[,!colnames(sero) %in% dropme] 

# load TMLE output ---------------------------------------
sero_models <- readRDS(paste0(results_path, "prevalence-tmle-results-sero.RDS"))


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
    if(parameter_name=="Direct effect")    model <- sero_models$res_de$human$res_est
    if(parameter_name=="Spillover effect") model <- sero_models$res_sp$human$res_est
    if(parameter_name=="Total effect")     model <- sero_models$res_te$human$res_est
  }
  
  if(reservoir_name=="Mosquito"){
    if(parameter_name=="Direct effect")   model <- sero_models$res_de$mosquito$res_est
    if(parameter_name=="Spillover effect") model <- sero_models$res_sp$mosquito$res_est
    if(parameter_name=="Total effect")     model <- sero_models$res_te$mosquito$res_est
  }
  
  if(reservoir_name=="Human & mosquito"){
    if(parameter_name=="Direct effect")   model <- sero_models$res_de$`human & mosquito`$res_est
    if(parameter_name=="Spillover effect") model <- sero_models$res_te$`human & mosquito`$res_est
    if(parameter_name=="Total effect")     model <- sero_models$res_sp$`human & mosquito`$res_est
  }
  
  if(!is.null(model)){
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
      dplyr::select(reservoir, parameter, N_int, N_cont, prev_int, prev_control, PR_unadj, lb_unadj, ub_unadj,
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
}



parameter_list <- list("Direct effect", "Spillover effect","Total effect")
reservoir_list <- list("Human", "Mosquito", "Human & mosquito")

## make table  ---------------------------------------
sero_list = list()
for(i in 1:length(parameter_list)){
  sero_list[[i]] = map_dfr(reservoir_list,
          function(x)
            make_prev_table(
              data = sero,
              reservoir_name = x,
              parameter_name = parameter_list[[i]],
              yname = "Etramp5.Ag1_pos"
            ))
}

sero_output = bind_rows(sero_list) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito"))) %>% 
  arrange(reservoir, parameter)

saveRDS(sero_output, paste0(namibia_process_path, "plot-data-prev-primary-sero.RDS"))

sero_table <- sero_output %>% dplyr::select(-c(PR_unadj, lb_unadj, ub_unadj))


# save tables ---------------------------------------
write.csv(sero_table, paste0(table_path, "table-prev-primary-sero.csv"))


