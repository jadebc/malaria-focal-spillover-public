
fit_indiv_glm = function(data, effecttype, covarname, unadj_est=F){
  if (effecttype == 'direct'){
    df <- data %>% filter(target_area == 1, int_recip == 1)
  }else if (effecttype == 'spillover'){
    df <- data %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
  }else if (effecttype == 'direct_prev'){
    df <- data %>% filter(target_area == 1)
  }
  if (effecttype == 'total') df = data
  
  # check at least 10 events per variable
  check_EPV = df %>% group_by(Y, A) %>% summarise(n=n())
  if(min(check_EPV$n)<10) print("Too few observations to fit model")
  
  fit_model = min(check_EPV$n)>=10
  
  if(fit_model){
    
    # check sufficient data for adjusted model
    if(!unadj_est){
      if(min(check_EPV$n)<= 30){
        print("Too few observations to fit adjusted model")
        unadj_est = TRUE
      }else{
        unadj_est = FALSE
      } 
    }

    
    # check that A is binary
    assert_that(setequal(unique(df$A), c(0,1)),
                msg = "A is not binary")
    
    # Step 1. screening
    if(!unadj_est){
      
      dmodel = df %>% droplevels()
      
      print("Covariate screening ---------------------------------")
      tic()
      df_Y = dmodel$Y
      df_Ws = dmodel[,c(covarname)]
      
      res_screen <- cov_prescreen(
        df = dmodel,
                                  Y = df_Y,
                                  Ws = df_Ws,
                                  yname = 'Y',
                                  covarname = covarname,
                                  family="binomial", 
                                  pval=0.2, 
                                  print=TRUE)
      toc()
      
      # check that there are no covariates of character class
      Ws = dmodel[,res_screen]
      class_list =  map(Ws, class) %>% t()
      classes = names(table(unlist(class_list)))
      assert_that(!"character" %in% classes, msg = "W includes characters")
   
      # fit individual model
      indiv_formula = paste0("Y ~ A +", paste0(res_screen, collapse = "+")) %>% as.formula()
      indiv_glm_fit = glm(indiv_formula, data = dmodel, family = "poisson"(link=log))
      vcovCL <- sandwichSE(fm = indiv_glm_fit, cluster = dmodel$id)
      rfit <- coeftest(indiv_glm_fit, vcovCL)
    }else{
      dmodel = df %>% dplyr::select(Y, A, id) %>% droplevels()
      # fit individual model
      indiv_glm_fit = glm("Y~A", data = dmodel, family = "poisson"(link=log))
      
      vcovCL <- sandwichSE(fm = indiv_glm_fit, cluster = dmodel$id)
      rfit <- coeftest(indiv_glm_fit, vcovCL)
      
    }

    
    RR <- data.frame(N = nrow(dmodel),
                     py1 = mean(dmodel$Y[dmodel$A==1]),
                     py0 = mean(dmodel$Y[dmodel$A==0]),
                     ptest = round(exp(rfit[2, 1]), 8), 
                     lb = round(exp(rfit[2, 1] - 1.96 * rfit[2, 2]), 8), 
                     ub = round(exp(rfit[2, 1] + 1.96 * rfit[2, 2]), 8))
    
    RR_format <- RR %>% mutate(
      py1 = sprintf("%0.03f", py1),
      py0 = sprintf("%0.03f", py0),
      ptest = sprintf("%0.02f", ptest),
      lb = sprintf("%0.02f", lb),
      ub = sprintf("%0.02f", ub)
    )
    
    return(list(RR = RR, RR_format = RR_format  ))
    
  }else{
    print("Data too sparse to fit model.")
  }
}



fit_indiv_glm_unadj = function(data, effecttype){
  if (effecttype == 'direct'){
    df <- data %>% filter(target_area == 1, int_recip == 1)
  }else if (effecttype == 'spillover'){
    df <- data %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
  }else if (effecttype == 'direct_prev'){
    df <- data %>% filter(target_area == 1)
  }
  if (effecttype == 'total') df = data
  
  dmodel = df %>% dplyr::select(Y, A, id) %>% droplevels()
  
  # fit individual model
  indiv_glm_fit = glm("Y~A", data = dmodel, family = "poisson"(link=log))
  res=summary(indiv_glm_fit)$coeff
  vcovCL <- sandwichSE(fm = indiv_glm_fit, cluster = dmodel$id)
  rfit <- coeftest(indiv_glm_fit, vcovCL)
  
  RR <- data.frame(N = nrow(dmodel),
                   py1 = mean(dmodel$Y[dmodel$A==1]),
                   py0 = mean(dmodel$Y[dmodel$A==0]),
                   ptest = round(exp(rfit[2, 1]), 8), 
                   lb = round(exp(rfit[2, 1] - 1.96 * rfit[2, 2]), 8), 
                   ub = round(exp(rfit[2, 1] + 1.96 * rfit[2, 2]), 8))
  
  return(RR)
  
}
