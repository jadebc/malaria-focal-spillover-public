################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Prevalence analysis
################################################  
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions/1-helper-functions-covar.R"))
library(sl3)
library(tmle3)
library(origami)
library(SuperLearner)
library(glmnet)
library(gam)
library(xgboost)
library(GGally)


# load and process analysis dataset ---------------------------------------------------
xs_clean <- readRDS(namibia_analysis_prev)

# specify the variables used for identification, but not for estimation
dropme <- c("eaid", "xsid", "iid", "latitude", "longitude",
            "date", "spillover_zone")

xs <- xs_clean[,!colnames(xs_clean) %in% dropme] 

# clean RDT results
xs <- xs %>% mutate(
  rdtyn = case_when(
    rdtresult == "Negative" ~ 0,
    rdtresult == "Pf only" ~ 1,
    rdtresult == "Pf or Mixed" ~ 1, 
    rdtresult == "Invalid" ~ NA_real_, 
    rdtresult == "Not conducted" ~ NA_real_
  ), 
  
  hsrdtyn = case_when(
    hsrdtresult == "Negative" ~ 0,
    hsrdtresult == "Positive" ~ 1,
    hsrdtresult == "Invalid" ~ NA_real_, 
    hsrdtresult == "Not conducted" ~ NA_real_
  )
)

# identify hot spots
hotspot1 <- xs %>% 
  mutate(travelyn = ifelse(traveled_recently=="Yes", 1, 0),
         occup_highrisk = ifelse(occupation=="High risk", 1, 0),
         sleep_net_yn = ifelse(sleep_net == "Always" | sleep_net == "Sometimes", 1, 0),
         slept_net_yest_yn = ifelse(slept_net_yest == "Yes", 1, 0)) %>% 
  group_by(arm, ea, hhid, target_area, spall) %>% 
  summarise(
    cases_qpcr = sum(qPCRposneg, na.rm=T), 
    cases_rdt = sum(rdtyn, na.rm=T),
    cases_hsrdt = sum(hsrdtyn, na.rm=T),
    total_travel = sum(travelyn, na.rm=T),
    total_occup_highrisk = sum(occup_highrisk, na.rm=T),
    total_sleep_net = sum(sleep_net_yn, na.rm=T),
    total_slept_net_yest = sum(slept_net_yest_yn, na.rm=T)
  )  %>% 
  mutate(cases_qpcr = ifelse(cases_qpcr>1, 1, 0),
         cases_rdt = ifelse(cases_rdt>1, 1, 0),
         cases_hsrdt = ifelse(cases_hsrdt>1, 1, 0))
  
hotspot2 <- xs %>% 
  group_by(arm, ea, hhid, target_area, spall) %>% 
  summarise_at(c("dist_popN_target",      
                 "dist_arm_target",        "dist_Ntx_target",       
                 "dist_Ntx_sametx_target", "dist_Ntx_difftx_target",
                 "dist_popN_1km",          "dist_arm_1km",          
                 "dist_Ntx_1km",           "dist_Ntx_sametx_1km",   
                 "dist_Ntx_difftx_1km",    "dist_popN_2km",         
                 "dist_arm_2km",           "dist_Ntx_2km",          
                 "dist_Ntx_sametx_2km",    "dist_Ntx_difftx_2km" ,  
                 "dist_popN_3km",          "dist_arm_3km" ,         
                 "dist_Ntx_3km",           "dist_Ntx_sametx_3km" ,  
                 "dist_Ntx_difftx_3km",       
                 "pre_incidence",          "pre_spray_cover",       
                 "pop_size_ea",            "dist_hh",               
                 "dist_hf",                "pre_rainfall",          
                 "pre_evi",                "ea_elevation",          
                 "surface_temp"  ), mean, na.rm=T) 

hotspot <- full_join(hotspot1, hotspot2, 
                     by = c("arm", "ea", "hhid", "target_area", "spall"))

saveRDS(hotspot, file = paste0(namibia_process_path, "namibia_prev_hotspot.RDS"))

# make tmle function ---------------------------------------------------
fit_tmle_prev <- function(data, parameter, reservoir, yname){
  
  print(paste0(parameter, "- ", reservoir, "- ", yname, "---------"))
  
  # Drop if outcome is missing
  data <- data %>% filter(!is.na(!!sym(yname)))
  
  # Drop observations outside spillover zone ----------------------
  if(parameter=="Direct effect"){
    data_sub <- data %>% filter(target_area==1)
  }
  if(parameter=="Spillover effect"){
    data_sub <- data %>% filter(spall==1)
  }
  if(parameter=="Total effect"){
    data_sub <- data %>% filter(target_area==1 | spall==1)
  }
  
  
  # binarize treatment variable
  # human reservoir
  if (reservoir == "human") {
    data_sub <- data_sub %>%
      mutate(             
        ## add a column for constructing interaction term
        inter = as.numeric(arm %in% c("RV", "TV")),
        arm_yn = case_when(arm %in% c("TO", "TV") ~ 1,
                           arm %in% c("RO", "RV") ~ 0)
      )
  }
  
  # mosquito reservoir
  if(reservoir == "mosquito"){
    data_sub <- data_sub %>%
      mutate(
        ## add a column for constructing interaction term
        inter = as.numeric(arm %in% c("TO", "TV")),
        arm_yn = case_when(arm %in% c("RV", "TV") ~ 1,
                           arm %in% c("RO", "TO") ~ 0)
      )
  }
  
  # human & mosquito reservoir
  if(reservoir == "human & mosquito"){
    data_sub <- data_sub %>% filter(arm %in% c("TV", "RO"))
    data_sub <- data_sub %>% mutate(arm_yn = case_when(arm %in% c("TV") ~ 1,
                                           arm %in% c("RO") ~ 0))
    
  }
  
  data_sub <- data_sub %>% ungroup() %>% dplyr::select(-c(arm, 
                                            spall, target_area, hhid, ea))
  
  # impute missing values  -------------------------------------
  nodes <- list(
    W = colnames(data_sub)[-which(names(data_sub) %in% c("arm_yn",
                                    "cases_qpcr", "cases_rdt", "cases_hsrdt"))],
    A = "arm_yn",
    Y = yname)
  
  data_nomiss <- process_missing(data = data_sub,
                                  node_list = nodes)$data
  
  data_nomiss <- as.data.frame(data_nomiss)
  
  # sparsity check  -------------------------------------
  checkN = data_nomiss %>% group_by(!!sym(yname), arm_yn) %>% summarise(n=n())
  
  # TMLE specs --------------------------------------
  tsm_spec <- tmle_TSM_all()
  
  rr_spec <- tmle_RR(
    baseline_level = 0,
    contrast_level = 1
  )
  
  rd_spec <- tmle_ATE(
    control_level = 0,
    treatment_level = 1
  )
  
  # SuperLearner set up ---------------------------------------
  lrn_mean <- Lrnr_mean$new()
  lrn_glm <- Lrnr_glm$new()
  lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
  lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
  lrn_elnet <- Lrnr_glmnet$new(alpha = 0.5)
  lrn_xgboost <- Lrnr_xgboost$new(eval_metric = "logloss")
  
  interactions <- list(c("inter", "arm_yn"))
  # main terms as well as the interactions above will be included
  lrn_interaction <- make_learner(Lrnr_define_interactions, interactions)
  # use a pipeline to combine the interaction learn with other learners 
  lrn_glm_interaction <- make_learner(Pipeline, lrn_interaction, lrn_glm)
  lrn_lasso_interaction <- make_learner(Pipeline, lrn_interaction, lrn_lasso)
  lrn_elnet_interaction <- make_learner(Pipeline, lrn_interaction, lrn_elnet)
  lrn_xgboost_interaction <- make_learner(Pipeline, lrn_interaction, lrn_xgboost)
  
  SL_lib_Q <- make_learner(Stack,
                           lrn_mean,
                           lrn_glm,
                           lrn_lasso,
                           lrn_elnet,
                           lrn_xgboost,
                           lrn_glm_interaction,
                           lrn_lasso_interaction,
                           lrn_elnet_interaction,
                           lrn_xgboost_interaction)
  
  SL_lib_g <- make_learner(Stack,
                           lrn_mean,
                           lrn_glm,
                           lrn_lasso,
                           lrn_elnet)
  
  
  ls_metalearner <- make_learner(Lrnr_nnls)
  lb_metalearner <- make_learner(Lrnr_solnp,
                                 learner_function = metalearner_logistic_binomial,
                                 loss_function = loss_loglik_binomial)
  sl_Y <- Lrnr_sl$new(
    learners = SL_lib_Q,
    metalearner = ls_metalearner
  )
  
  sl_A <- Lrnr_sl$new(
    learners = SL_lib_g,
    metalearner = lb_metalearner
  )
  
  learner_list <- list(A = sl_A, Y = sl_Y)
  
  # proceed with estimation as long as there are some positive cases
  # within this analysis dataset
  if(any(data_nomiss %>% pull(!!sym(yname)) != 0)){
    
    # covariate adjustment if at least 30 obs per Y, A strata -----------------
    if(min(checkN$n)>30){
  
      # list covariates  -------------------------------------
      covarname <- setdiff(names(data_nomiss), c(dropme,yname, 'arm','arm_yn', 
                                                  'target_area_sp1km','eaid','ea','inter'))
      
      # covariate screening  -------------------------------------
      df_Y = data_nomiss %>% pull(!!sym(yname)) 
      df_Ws = data_nomiss[,c(covarname)]
      
  
      
        res_screen <- cov_prescreen(df = data_nomiss,
                                      Y = df_Y,
                                      Ws = df_Ws,
                                      yname = yname,
                                      covarname = covarname,
                                      family="binomial",
                                      pval=0.2,
                                      print=TRUE)
       
        
        # update adjustment sets based on screening results ----------
        node_list_tmle3 <- list(
          W = c(res_screen, "inter"),
          A = 'arm_yn',
          Y = yname
        )
        
         # fit tmle --------------------------------------
        tmle_task <- tsm_spec$make_tmle_task(data_nomiss, node_list_tmle3, id = "ea")
        initial_likelihood <- tsm_spec$make_initial_likelihood(tmle_task, learner_list)
        targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = list(cvtmle = TRUE, 
                                                                                          convergence_type = "sample_size"))
        
        tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
        
        rd_param <- define_param(
          Param_delta, targeted_likelihood,
          delta_param_ATE,
          list(tsm_params[[1]], tsm_params[[2]])
        )
        
        rr_param <- define_param(
          Param_delta, targeted_likelihood,
          delta_param_RR,
          list(tsm_params[[1]], tsm_params[[2]])
        )
        
        all_params <- c(tsm_params, rd_param, rr_param)
        tmle_fit_all <- fit_tmle3(tmle_task, targeted_likelihood, all_params, targeted_likelihood$updater)
        
        tmle_adj_df <- tmle_fit_all$summary %>% mutate(model = "Adjusted")
        

    }else{
    
    # unadjusted models  ----------
    node_list_tmle3 <- list(
      W = NULL,
      A = 'arm_yn',
      Y = yname
    )

    lrnr_sl_a <- make_learner(Lrnr_sl,
                              learners =make_learner(Lrnr_mean))
    lrnr_sl_y <- make_learner(Lrnr_sl,
                              learners = make_learner(Lrnr_glm))
    
    learner_list_unadj <- list(Y = lrnr_sl_y, A = lrnr_sl_a)
    
    tmle_task <- tsm_spec$make_tmle_task(data_nomiss, node_list_tmle3, id = "ea")
    initial_likelihood <- tsm_spec$make_initial_likelihood(tmle_task, learner_list_unadj)
    targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = list(cvtmle = FALSE, 
                                                                                      convergence_type = "sample_size"))
    
    tsm_params <- tsm_spec$make_params(tmle_task, targeted_likelihood)
    
    rd_param <- define_param(
      Param_delta, targeted_likelihood,
      delta_param_ATE,
      list(tsm_params[[1]], tsm_params[[2]])
    )
    
    rr_param <- define_param(
      Param_delta, targeted_likelihood,
      delta_param_RR,
      list(tsm_params[[1]], tsm_params[[2]])
    )
    
    all_params <- c(tsm_params, rd_param, rr_param)
    tmle_fit_all <- fit_tmle3(tmle_task, targeted_likelihood, all_params, targeted_likelihood$updater)
    
    ## fit tmle --------------------------------------
    tmle_unadj_df <- tmle_fit_all$summary %>% mutate(model = "Unadjusted")
    res_screen <- NA
    }
    
    ## output results --------------------------------------
    if(min(checkN$n)>30){
      fits <- bind_rows(tmle_adj_df, tmle_unadj_df) %>% 
        mutate(param_type = case_when(
          param == "E[Y_{A=0}]" ~ "E[Y_{A=0}]",
          param == "E[Y_{A=1}]" ~ "E[Y_{A=1}]",
          type == "RR" ~ "RR",
          type == "ATE" ~ "RD"
        ), 
        outcome = yname)
    }else{
      fits <- tmle_unadj_df %>%
        mutate(
          param_type = case_when(
            param == "E[Y_{A=0}]" ~ "E[Y_{A=0}]",
            param == "E[Y_{A=1}]" ~ "E[Y_{A=1}]",
            type == "RR" ~ "RR",
            type == "ATE" ~ "RD"
          ), 
          outcome = yname
        )
    }
    
  }else{
    res_screen = NA
    fits = NA
  }
 
  
  
  return(list('res_screen' = res_screen,
              'res_est' = fits))
}


# run tmle function ---------------------------------------------------
outcome_list <- list("cases_qpcr", "cases_rdt", "cases_hsrdt")

## human reservoir ---- 
set.seed(123)
res_human_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Direct effect", reservoir = "human", yname= x)) 
res_human_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Spillover effect", reservoir = "human", yname= x))
res_human_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Total effect", reservoir = "human", yname= x))

names(res_human_de) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_human_sp) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_human_te) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")

## mosquito reservoir ---- 
set.seed(123)
res_mosq_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Direct effect", reservoir = "mosquito", yname= x)) 
res_mosq_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Spillover effect", reservoir = "mosquito", yname= x)) 
res_mosq_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Total effect", reservoir = "mosquito", yname= x)) 

names(res_mosq_de) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_mosq_sp) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_mosq_te) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")

## human and mosquito reservoir ---- 
set.seed(123)
res_both_de <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Direct effect", reservoir = "human & mosquito", yname= x)) 
res_both_sp <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Spillover effect", reservoir = "human & mosquito", yname= x)) 
res_both_te <- lapply(outcome_list, function(x) fit_tmle_prev(
  data = hotspot, parameter = "Total effect", reservoir = "human & mosquito", yname= x)) 

names(res_both_de) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_both_sp) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")
names(res_both_te) = c("cases_qpcr", "cases_rdt", "cases_hsrdt")

# save results ---------------------------------------------------
prev_tmle_results <- list(
  res_human_de = res_human_de, 
  res_human_sp = res_human_sp,
  res_human_te = res_human_te,
  
  res_mosq_de = res_mosq_de,
  res_mosq_sp = res_mosq_sp,
  res_mosq_te = res_mosq_te,
  
  res_both_de = res_both_de,
  res_both_sp = res_both_sp,
  res_both_te = res_both_te
  
)

saveRDS(prev_tmle_results, file = paste0(results_path, "prevalence-hotspot-tmle-results.RDS"))
