#======================================================================================
#  A wrapper function that first do screening, then estimation and inference.
#  component scripts: pre_screen, hTMLE, inference
# 
#======================================================================================

run_hTMLE <- function(df,
                      effecttype,
                      Qlevel = 'individual', # cohort/individual/adaptive for incidence, individual for prev
                      glevel = 'individual', # always cohort for incidence, individual for prev
                      unadj_est = FALSE, # TRUE if we want a unadjusted estimate
                      dependency = FALSE, # TRUE if we want CIs adjusted for dependency
                      y_l = 0,
                      y_u = 1, ## TEMP think more
                      covarname_g, # set of names for adjustment variables for g
                      covarname_Q, # set of names for adjustment variables for Q
                      intername = NULL, # set of names for interaction variables
                      SL_lib_Q,
                      SL_lib_g,
                      SL_lib_adaQ = NULL, # need to specify if for Qlevel = 'adaptive'
                      SL_lib_depen = NULL, # need to specify if dependency = T
                      verbose = FALSE,
                      reservoir_name,
                      Vfolds = 10,
                      EM = F, # indicator for whether it is an effect modification analysis,
                      EMname = NULL # name of effect modifier list element 
){
  print("Fitting htmle model =======================================")
  
  # initial checks ---------  =========================================
  ## check the validity of input parameters ---------
  assert_that(effecttype %in% c('total', 'direct', 'direct_prev', 'spillover'),
              msg = paste0('effecttype must be one of: total, direct, direct_prev, or spillover'))
  
  assert_that(Qlevel %in% c('cohort', 'individual', 'adaptive'),
              msg = paste0('Qlevel must be one of: cohort, individual, adaptive'))
  
  if (dependency){
    assert_that(!is.null(SL_lib_depen),
                msg = paste0('SL_lib_depen must be specified if dependency = T'))
  }
  
  if (Qlevel == "adaptive"){
    assert_that(!is.null(SL_lib_adaQ),
                msg = paste0("SL_lib_adaQ must be specified if Qlevel == 'adaptive'"))
  }
  
  ## check if data contains all necessary variables ---------
  assert_that('Y' %in% names(df),
              msg = paste0("did not find 'Y' variable in df"))
  
  assert_that('A' %in% names(df),
              msg = paste0("did not find 'A' variable in df"))
  
  assert_that('id' %in% names(df),
              msg = paste0("did not find 'id' variable in df"))
  
  assert_that('alpha' %in% names(df),
              msg = paste0("did not find 'alpha' variable in df"))
  
  ## check treatment coding
  assert_that(length(unique(df$A)) == 2,
              msg = paste0("A should be binary"))
  
  res_screen <- NULL
  res_adap_Q <- NULL
  res_est <- NULL
  
  ## select data based on 'dependency', 'effecttype'  ---------
  if (!dependency){
    df_overlap = NULL
  }else{
    df_overlap = check_overlap_all(df = df, r = 1000)
    if (effecttype == 'direct' | effecttype == 'direct_prev'){
      df_overlap = df_overlap %>% filter(overlap_type =='target')
    }
  }
  
  if (effecttype == 'direct'){
    df <- df %>% filter(target_area == 1, int_recip == 1)
  }else if (effecttype == 'spillover'){
    df <- df %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
  }else if (effecttype == 'direct_prev'){
    df <- df %>% filter(target_area == 1)
  }
  
  ## check at least 10 events per variable  ---------
  check_EPV = df %>% group_by(Y, A) %>% summarise(n=n())
  if(min(check_EPV$n)<10) print("Too few observations to fit model")
  
  # determine whether a model should be fit  --------- =========================================
  fit_model = min(check_EPV$n)>=10
  
  if(fit_model){
    
    # check sufficient data for adjusted model
    if(min(check_EPV$n)<30){
      print("Too few observations to fit adjusted model")
      unadj_est = TRUE
      intername = NULL
    } 
    
    # check for data sparsity with interaction term
    if(!is.null(intername)) check_inter = df %>% group_by(Y, A, !!sym(intername)) %>% summarise(n=n())
    if(!is.null(intername)) assert_that(min(check_inter$n)>=10, msg = "Too few observations to do interaction")
    
    
    ## adjusted model workflow (individual-level)  ---------
    if(!unadj_est){
      
      ## drop effect modifiers from covariate list  ---------
      if(EM){
        if(length(grep("sex",EMname))>0){
          if(length(which(covarname_g == "sex"))>0) covarname_g = covarname_g[-which(covarname_g == "sex")]
          if(length(which(covarname_Q == "sex"))>0) covarname_Q = covarname_Q[-which(covarname_Q == "sex")]
          
        }else{
          EM_name = substr(unique(df$modifier_var), 1, 
                           str_locate(unique(df$modifier_var), "_abovemed")[1] - 1)
          if(length(which(covarname_g == EM_name))>0) covarname_g = covarname_g[-which(covarname_g == EM_name)]
          if(length(which(covarname_Q == EM_name))>0) covarname_Q = covarname_Q[-which(covarname_Q == EM_name)]
          
        }
      }
      
      ## check for positivity violations  ---------
      # pos_data = positivity_screening(data = df, varlist = c(covarname_g, covarname_Q))
      # nopos_cov = pos_data$covariates
      # nopos_df = pos_data$data
      nopos_cov = c(covarname_g, covarname_Q)
      nopos_df = df
      
      covar_selected_g <- apply(as.data.frame(covarname_g), 1, 
                                function(x) nopos_cov[grep(x, nopos_cov)]) %>% unlist()
      covar_nopos_Q <- apply(as.data.frame(covarname_Q), 1, 
                             function(x) nopos_cov[grep(x, nopos_cov)]) %>% unlist()
      
      
      ## covariate screening  ---------
      print("Covariate screening ---------------------------------")
      
      
      # covariate screening for g model 
      df_Ws = nopos_df %>% dplyr::select(all_of(covar_selected_g))
      covar_selected_g <- cov_prescreen(df = nopos_df,
                                               Y = nopos_df$A,
                                               Ws = df_Ws,
                                               yname = 'A',
                                               covarname = covar_selected_g,
                                               family="binomial",
                                               pval=0.2,
                                               print=TRUE)
        
      
      
      # covariate screening for Q model 
      tic()
      df_Y = nopos_df$Y
      df_Ws = nopos_df[,c(covar_nopos_Q)]
      covar_selected_Q <- cov_prescreen(df = nopos_df,
                                        Y = df_Y,
                                        Ws = df_Ws,
                                        yname = 'Y',
                                        covarname = covar_nopos_Q,
                                        family="binomial", 
                                        pval=0.2, 
                                        print=TRUE)
      toc()
      
      
      # check that there are no covariates of character class
      Ws = nopos_df[,c(covar_selected_g, covar_selected_Q)]
      class_list =  map(Ws, class) %>% t()
      classes = names(table(unlist(class_list)))
      assert_that(!"character" %in% classes, msg = "W includes characters")
      
      # check that A is binary
      assert_that(setequal(unique(nopos_df$A), c(0,1)),
                  msg = "A is not binary")
    }
    
    ## adaptively choose individual vs. cohort level data  -----------------
    print("Choose Qlevel ---------------------------------")
    tic()
    if (Qlevel == "adaptive"){
      
      # data-adaptively choose Q model
      if(!unadj_est){
        res_adap_Q <- do_adap_Q(train = nopos_df,
                                # QAdj = covarname,
                                QAdj = covar_selected_Q,
                                SL_lib = SL_lib_adaQ, 
                                V = Vfolds,
                                verbose = F)
      }else{
        res_adap_Q <- do_adap_Q(train = df,
                                QAdj = NULL,
                                SL_lib = SL_lib_adaQ, 
                                V = Vfolds,
                                verbose = F)
      }
      
      
      mse_coho <- mean((res_adap_Q$res_cv$Y - res_adap_Q$res_cv$pred_coho)^2)
      mse_indiv <- mean((res_adap_Q$res_cv$Y - res_adap_Q$res_cv$pred_indiv)^2)
      
      Qlevel = ifelse(mse_coho <= mse_indiv, "cohort", "individual")
      
      print(paste0('Data-adaptive Q model selection result: ',
                   Qlevel, " level"))
    }
    toc()
    
    
    
    
    ## make cohort-level data for g and/or Q ------------------------------
    df_coho <- aggregate(x = df,
                         by = list(id=df$id),
                         hp_aggregate)[,2:(ncol(df)+1)]
    
    # transform Y on cohort level since it is continuous
    df_coho$Y = (df_coho$Y - y_l) / (y_u - y_l)
    
    ## adjusted model workflow (cohort-level) ------------------------------
    if(!unadj_est){
      
      ### check for positivity violations (cohort-level) ---------------------------
      # pos_coho_data = positivity_screening(data = df_coho, varlist = c(covarname_g, covarname_Q))
      # nopos_cov_coho = pos_coho_data$covariates
      # nopos_coho_df = pos_coho_data$data
      nopos_cov_coho = c(covarname_g, covarname_Q)
      nopos_coho_df = df_coho
      
      covar_selected_g_coho <- apply(as.data.frame(covarname_g), 1, 
                                     function(x) nopos_cov_coho[grep(x, nopos_cov_coho)]) %>% unlist()
      covar_nopos_Q_coho <- apply(as.data.frame(covarname_Q), 1, 
                                  function(x) nopos_cov_coho[grep(x, nopos_cov_coho)]) %>% unlist()
      
      
      ## covariate screening (cohort-level) ------------------------------
      print("Cohort-level covariate screening ---------------------------------")
      df_Y = nopos_coho_df$Y
      
      # covariate screening for g model 
      if(length(covar_selected_g_coho)>0){
        df_Ws = nopos_coho_df %>% dplyr::select(all_of(covar_selected_g_coho))
        covar_selected_g_coho <- cov_prescreen(df = nopos_coho_df,
                                               Y = nopos_coho_df$A,
                                               Ws = df_Ws,
                                               yname = 'A',
                                               covarname = covar_selected_g_coho,
                                               family="binomial",
                                               pval=0.2,
                                               print=TRUE)
        
      }
      
      # covariate screening for Q model 
      if(length(covar_nopos_Q_coho)>0){
        df_Ws = nopos_coho_df %>% dplyr::select(all_of(covar_nopos_Q_coho))
        covar_selected_Q_coho <- cov_prescreen(df = nopos_coho_df,
                                               Y = df_Y,
                                               Ws = df_Ws,
                                               yname = 'Y',
                                               covarname = covar_nopos_Q_coho,
                                               family="gaussian",
                                               pval=0.2,
                                               print=TRUE)
        
        df_coho <- nopos_coho_df
      }
    }
    
    
    ## create small regression learners for g model -----------------
    
    # if number of observations per covariate <20, 
    # use smaller regression models in SL library 
    # to prevent positivity violations
    
    ### check if small learners are needed --------
    #### g model  --------
    if(!unadj_est){
      check_n_per_cov_g <- nrow(nopos_df) / length(covar_selected_g_coho) > 20
    }else{
      if(unadj_est) check_n_per_cov_g = T
    } 
    
    if(!check_n_per_cov_g) print("<20 obs per covariate for g model; creating small regression learners")
    
    #### Q model  --------
    if(Qlevel != "cohort" & !unadj_est) check_n_per_cov_Q <-
      nrow(nopos_df) / length(covar_selected_Q) > 20
    
    if(Qlevel == "cohort" & !unadj_est)  check_n_per_cov_Q <-
      length(df_coho) / length(covar_selected_Q_coho) > 20
    
    if(unadj_est) check_n_per_cov_Q = T
    
    if(!check_n_per_cov_Q) print("<20 obs per covariate for g model; creating small regression learners")
    
    #### create small learners (g model) --------
    if(!unadj_est & !check_n_per_cov_g){
      
      for(i in 1:length(covar_selected_g_coho)){
        lrner_g_name <- paste0("lrn_glm_g_small_", i)
        assign(lrner_g_name, Lrnr_glm_fast$new(covariates = covar_selected_g_coho[i]))
      }
      
      original <- list(Stack)
      
      g_small_regressions_chr <- as.list(objects()[grep("lrn_glm_g_small_",
                                                        objects())])
      g_small_regressions = lapply(g_small_regressions_chr, function(x)
        eval(parse(text=x)))
      
      g_params <- c(original, g_small_regressions)
      SL_lib_g_small <- do.call(make_learner, g_params)
      
      SL_lib_g <- SL_lib_g_small
    }
    
    
    ## create small regressions learners for Q model (individual-level)  -----------------
    # for individual level Q model 
    if(!unadj_est & Qlevel != "cohort" & !check_n_per_cov_Q){
      
      for(i in 1:length(covar_selected_Q)){
        lrner_Q_name <- paste0("lrn_glm_Q_small_", i)
        assign(lrner_Q_name, Lrnr_glm_fast$new(covariates = c("A", covar_selected_Q[i])))
        
        original <- list(Stack)
        
        Q_small_regressions_chr <- as.list(objects()[grep("lrn_glm_Q_small_",
                                                          objects())])
        Q_small_regressions = lapply(Q_small_regressions_chr, function(x)
          eval(parse(text=x)))
        
        Q_params <- c(original, Q_small_regressions)
        SL_lib_Q_small <- do.call(make_learner, Q_params)
        
        SL_lib_Q <- SL_lib_Q_small
        SL_lib_adaQ <- SL_lib_Q_small
      }
    }
    
    ## create small regressions learners for Q model (cohort-level)  -----------------
    
    # create learners for small regressions
    # for cohort level Q model
    if(!unadj_est & Qlevel == "cohort" & !check_n_per_cov_Q){
      
      for(i in 1:length(covar_selected_Q_coho)){
        lrner_Q_name <- paste0("lrn_glm_Q_small_", i)
        assign(lrner_Q_name, Lrnr_glm_fast$new(covariates = c("A", covar_selected_Q_coho[i])))
      }
      
      original <- list(Stack)
      
      Q_small_regressions_chr <- as.list(objects()[grep("lrn_glm_Q_small_",
                                                        objects())])
      Q_small_regressions = lapply(Q_small_regressions_chr, function(x)
        eval(parse(text=x)))
      
      Q_params <- c(original, Q_small_regressions)
      SL_lib_Q_small <- do.call(make_learner, Q_params)
      
      SL_lib_Q <- SL_lib_Q_small
    }
    
    
    ##  estimation and inference  ----------------- =========================================
    print("hTMLE estimation ---------------------------------")
    tic()
    
    ###  adjusted models  ----------------- 
    
    if(!unadj_est){
      if(Qlevel=="cohort"){
        res_est <- try(do_hTMLE(df = nopos_df,
                                df_coho = df_coho,
                                Qlevel = Qlevel, # cohort/individual for incidence, individual for prev
                                glevel = glevel, # always cohort for incidence, individual for prev
                                QAdj = unique(c(covar_selected_Q_coho,intername)), 
                                gAdj = covar_selected_g_coho, # usually = QAdj. use any covar if unadj_est = T
                                SL_lib_Q = SL_lib_Q, 
                                SL_lib_g = SL_lib_g, 
                                dependency = dependency,
                                SL_lib_depen = SL_lib_depen, # need if dependency = T
                                df_overlap = df_overlap, # need if dependency = T
                                unadj_est = unadj_est, 
                                y_l = y_l,
                                y_u = y_u, 
                                Vfolds = Vfolds,
                                EM = EM))
        
      }else{
        res_est <- try(do_hTMLE(df = nopos_df,
                                df_coho = df_coho,
                                Qlevel = Qlevel, # cohort/individual for incidence, individual for prev
                                glevel = glevel, # always cohort for incidence, individual for prev
                                QAdj = unique(c(covar_selected_Q,intername)), 
                                gAdj = covar_selected_g_coho, # usually = QAdj. use any covar if unadj_est = T
                                SL_lib_Q = SL_lib_Q, 
                                SL_lib_g = SL_lib_g, 
                                dependency = dependency,
                                SL_lib_depen = SL_lib_depen, # need if dependency = T
                                df_overlap = df_overlap, # need if dependency = T
                                unadj_est = unadj_est, 
                                y_l = y_l,
                                y_u = y_u, 
                                Vfolds = Vfolds,
                                EM = EM))
        
      }
      
    }else{
      ###  unadjusted models  ----------------- 
      
      res_est <- try(do_hTMLE(df = df,
                              df_coho = df_coho,
                              Qlevel = Qlevel, # cohort/individual for incidence, individual for prev
                              glevel = glevel, # always cohort for incidence, individual for prev
                              QAdj = NULL,
                              gAdj = NULL, 
                              SL_lib_Q = SL_lib_Q, 
                              SL_lib_g = SL_lib_g, 
                              dependency = dependency,
                              SL_lib_depen = SL_lib_depen, # need if dependency = T
                              df_overlap = df_overlap, # need if dependency = T
                              unadj_est = T, 
                              y_l = y_l,
                              y_u = y_u, 
                              Vfolds = Vfolds,
                              EM = EM
      ))
    }
    
    toc()
    
    ### Add labels to results data frame  -----------------
    
    if(all(!is.na(res_est))){
      res_est$estimates = res_est$estimates %>% mutate(
        reservoir = reservoir_name,
        parameter = effecttype
      )
    }
    
    if(!unadj_est){
      if(Qlevel =="cohort"){
        return(list('res_screen_g' = covar_selected_g_coho,
                    'res_screen_Q' = covar_selected_Q_coho,
                    'res_adap_Q' = res_adap_Q,
                    'res_est' = res_est))  
      }else{
        return(list('res_screen_g' = covar_selected_g_coho,
                    'res_screen_Q' = covar_selected_Q,
                    'res_adap_Q' = res_adap_Q,
                    'res_est' = res_est))  
      }
      
    }else{
      return(list('res_screen' = NA,
                  'res_adap_Q' = res_adap_Q,
                  'res_est' = res_est))   
    }
    
    
    print("Estimation complete ---------------------------------")
  }
}







