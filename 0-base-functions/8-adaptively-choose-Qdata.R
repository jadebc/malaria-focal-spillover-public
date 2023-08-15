#======================================================================================
#  Data-adaptively choose Q model level (individual or cohort)
#======================================================================================

choose_Qlevel <- function(df,
                      effecttype,
                      unadj_est = FALSE, # TRUE if we want a unadjusted estimate
                      dependency = FALSE, # TRUE if we want CIs adjusted for dependency
                      y_l = 0,
                      y_u = 1, ## TEMP think more
                      covarname_Q, # set of names for adjustment variables for Q
                      intername = NULL, # set of names for interaction variables
                      SL_lib_Q,
                      SL_lib_adaQ = NULL, # need to specify if for Qlevel = 'adaptive'
                      SL_lib_depen = NULL, # need to specify if dependency = T
                      verbose = FALSE,
                      reservoir_name,
                      EM = F, # indicator for whether it is an effect modification analysis
                      Vfolds = 10
){

  # initial checks ---------  =========================================
  ## check the validity of input parameters ---------
  assert_that(effecttype %in% c('total', 'direct', 'direct_prev', 'spillover'),
              msg = paste0('effecttype must be one of: total, direct, direct_prev, or spillover'))

  if (dependency){
    assert_that(!is.null(SL_lib_depen),
                msg = paste0('SL_lib_depen must be specified if dependency = T'))
  }
  

  assert_that(!is.null(SL_lib_adaQ),
                msg = paste0("SL_lib_adaQ must be specified if Qlevel == 'adaptive'"))
  
  
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
        EM_name = substr(unique(df$modifier_var), 1, 
                         str_locate(unique(df$modifier_var), "_abovemed")[1] - 1)
        if(length(which(covarname_g == EM_name))>0) covarname_g = covarname_g[-which(covarname_g == EM_name)]
        if(length(which(covarname_Q == EM_name))>0) covarname_Q = covarname_Q[-which(covarname_Q == EM_name)]
      }
      
      ## check for positivity violations  ---------
      pos_data = positivity_screening(data = df, varlist = c(covarname_g, covarname_Q))
      nopos_cov = pos_data$covariates
      nopos_df = pos_data$data
      
      covar_selected_g <- apply(as.data.frame(covarname_g), 1, 
                                function(x) pos_data$covariates[grep(x, pos_data$covariates)]) %>% unlist()
      covar_nopos_Q <- apply(as.data.frame(covarname_Q), 1, 
                             function(x) pos_data$covariates[grep(x, pos_data$covariates)]) %>% unlist()
      
      
      ## covariate screening  ---------
      print("Covariate screening ---------------------------------")
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
    # tic()

    # data-adaptively choose Q model
    if(!unadj_est){
      tic()
      res_adap_Q <- do_adap_Q(train = nopos_df,
                              QAdj = covar_selected_Q,
                              SL_lib = SL_lib_adaQ, 
                              V = Vfolds,
                              verbose = F)
      toc()
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
    
    # toc()
    
    return(Qlevel)
    
    
  }
}







