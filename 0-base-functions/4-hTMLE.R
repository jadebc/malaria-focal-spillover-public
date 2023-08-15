# 0909 Let us try use sl2 and tmle3 
# use Laura's methodology, but not her previous implementation

# ......................................................
# Adaptively fit the Q model using the idea of discrete SL
# ......................................................
# major inputs: 
#       train, individual level data set (need "id", "Y", "A, "alpha")
#       QAdj, covar set for outcome
do_adap_Q <- function(train, QAdj, SL_lib, V = 10, verbose=F){
  # aggregate to cohort level
  train_coho <- aggregate(x=train, by=list(id=train$id), hp_aggregate)[,2:(ncol(train)+1)]
  # weights for aggregated data should be 1
  train_coho$alpha<- 1
  
  folds = make_folds(train_coho, fold_fun = folds_vfold, V = V)
  
  res_cv <- data.frame('id' = train_coho$id,
                       'Y' = train_coho$Y,
                       'pred_coho' = NA,
                       'pred_indiv' = NA)
  
  for (i in 1:length(folds)){
    # cohort level
    df_t_coho <- train_coho[folds[[i]]$training_set,]
    df_v_coho <- train_coho[folds[[i]]$validation_set,]
    
    df_t_coho_covar <- df_t_coho[, c(QAdj, 'A')] 
    df_v_coho_covar <- df_v_coho[, c(QAdj, 'A')] 
    
    # use sl3
    # define training task
    sl3_task <- make_sl3_Task(
      data = df_t_coho,
      covariates = c(QAdj,'A'),
      outcome = 'Y',
      weights = 'alpha',
      id = 'id'
    )
    
    # train
    sl <- make_learner(Lrnr_sl, learners = SL_lib)
    sl_fit <- sl$train(sl3_task)
    
    # define predicting task
    sl3_task <- make_sl3_Task(
      data = df_v_coho,
      covariates = c(QAdj,'A'),
      outcome = 'Y',
      weights = 'alpha',
      id = 'id'
    )
    
    pred_Q_coho <- sl_fit$predict(sl3_task)
    
    df_pred_Q_coho <- data.frame('id' = df_v_coho$id, 
                                 'pred_Q' = pred_Q_coho)
    
    res_cv$pred_coho[res_cv$id %in% df_pred_Q_coho$id] <- df_pred_Q_coho$pred_Q
    
    # individual level
    df_t <- train %>% filter(id %in% train_coho$id[folds[[i]]$training_set]) 
    df_v <- train %>% filter(id %in% train_coho$id[folds[[i]]$validation_set]) 
    
    df_t_covar <- df_t[, c(QAdj, 'A')] 
    df_v_covar <- df_v[, c(QAdj, 'A')] 
    
    # use sl3
    # define training task
    sl3_task <- make_sl3_Task(
      data = df_t,
      covariates = c(QAdj,'A'),
      outcome = 'Y',
      weights = 'alpha',
      id = 'id'
    )
    
    # train
    sl <- make_learner(Lrnr_sl, learners = SL_lib)
    sl_fit <- sl$train(sl3_task)
    
    # define predicting task
    sl3_task <- make_sl3_Task(
      data = df_v,
      covariates = c(QAdj,'A'),
      outcome = 'Y',
      weights = 'alpha',
      id = 'id'
    )
    
    pred_Q <- sl_fit$predict(sl3_task)
    
    df_pred_Q <- data.frame('id' = df_v$id, 
                            'pred_Q' = pred_Q)
    df_pred_Q <- aggregate(df_pred_Q$pred_Q, by=list(id=df_pred_Q$id), mean)		
    
    res_cv$pred_indiv[res_cv$id %in% df_pred_Q$id] <- df_pred_Q[,2]
  }
  
  mse_coho <- mean((res_cv$Y - res_cv$pred_coho)^2)
  mse_indiv <- mean((res_cv$Y - res_cv$pred_indiv)^2)
  
  res <- list('res_cv' = res_cv,
              'mse_coho' = mse_coho,
              'mse_indiv' = mse_indiv)
  
  return(res)
}	

# ......................................................
# Helper functions
# ......................................................
hp_aggregate <- function(x){
  if (class(x) == "numeric"){
    mean(x)
  }else{
    count_df <- aggregate(count ~ x, data = data.frame(count = 1, x = x), sum)
    mode <- count_df$x[which.max(count_df$count)]
    mode
  }
}

# ......................................................
# The main function for estimation 
# ......................................................
do_hTMLE <- function(df, 
                     df_coho,
                     Qlevel, # cohort/individual for incidence, individual for prev
                     glevel, # always cohort for incidence, individual for prev
                     QAdj, # set of covariate, not include A
                     gAdj, # usually = QAdj
                     SL_lib_Q,
                     SL_lib_g, 
                     dependency = F,
                     SL_lib_depen = NULL, # need if dependency = T
                     df_overlap = NULL, # need if dependency = T
                     unadj_est = F, 
                     y_l = 0,
                     y_u = 1,
                     Vfolds = 10,
                     EM = F){
  
  # Part 1. Estimation .................................
  
  # the number of cohorts
  n_coho <- length(unique(df$id))

  # # cohort level data df_coho
  # if (Qlevel=="cohort"){
  #   if (length(unique(df$Y)) > 2){
  #     # transform Y if Y is continuous
  #     df$Y = (df$Y - y_l) / (y_u - y_l)
  #   }else{
  #     # identity transform (no transform) if Y is binary
  #     y_u = 1
  #     y_l = 0
  #   }
  #   df_coho = df
  # }else{
  #   # if Q level is individual, make cohort level data for g estimation
  #   df_coho <- aggregate(x = df,
  #                        by = list(id=df$id),
  #                        hp_aggregate)[,2:(ncol(df)+1)]
  #   
  #   # transform Y on cohort level since it is continuous
  #   df_coho$Y = (df_coho$Y - y_l) / (y_u - y_l)
  # }


  if (unadj_est){
    SL_lib_Q = make_learner(Lrnr_glm)
    SL_lib_g = make_learner(Lrnr_mean)
    QAdj = NULL
  }
  
  # 1.1. fit Q .................................
  if (Qlevel == "cohort"){
    # fit Q model on cohort level
    res_fit_Q <- fit_Q(df = df_coho, SL_lib_Q = SL_lib_Q, 
                       QAdj = QAdj, unadj_est = unadj_est,
                       Vfolds = Vfolds)
    df_Q <- res_fit_Q$df_Q
    
  }else{
    # fit Q model on individual level
    res_fit_Q <- fit_Q(df = df, SL_lib_Q = SL_lib_Q, 
                       QAdj = QAdj, unadj_est = unadj_est,
                       Vfolds = Vfolds)
    df_Q <- res_fit_Q$df_Q
    # aggregate data to cohort level after fitting (if there is any cohort)
    if (n_coho != nrow(df)){
      df_Q <- aggregate(x=df_Q, by=list(id=df$id), mean)
      if (length(unique(df$Y)) == 2){
        # transform
        df_Q$QbarAW <- (df_Q$QbarAW - y_l) / (y_u - y_l)
        df_Q$Qbar1W <- (df_Q$Qbar1W - y_l) / (y_u - y_l)
        df_Q$Qbar0W <- (df_Q$Qbar0W - y_l) / (y_u - y_l)
      }
    }
  }
  
  QbarAW <- df_Q$QbarAW
  Qbar1W <- df_Q$Qbar1W
  Qbar0W <- df_Q$Qbar0W
  
  # print(paste0("initial Qbar1W, Qbar0W: ", mean(Qbar1W),", ", mean(Qbar0W)))
  
  # bound pred above 0 ## TEMP!!!
  QbarAW<- pmax(QbarAW, 10^(-20))
  Qbar1W<- pmax(Qbar1W, 10^(-20))
  Qbar0W<- pmax(Qbar0W, 10^(-20))
  
  # print(paste0("initial bounded Qbar1W, Qbar0W: ", mean(Qbar1W),", ", mean(Qbar0W)))
  
  # 1.2. fit g .................................
  if (glevel == "cohort"){
    # fit Q model on cohort level
    res_fit_g <- fit_g(df = df_coho, 
                       SL_lib_g = SL_lib_g, 
                       gAdj = gAdj, 
                       unadj_est = unadj_est,
                       Vfolds = Vfolds)
    g_hat <- res_fit_g$g_hat
  }else{
    # fit g model on individual level 
    res_fit_g <- fit_g(df = df, 
                       SL_lib_g = SL_lib_g, 
                       gAdj = gAdj, 
                       unadj_est = unadj_est,
                       Vfolds = Vfolds)
    g_hat <- res_fit_g$g_hat
    # aggregate data to cohort level after fitting (if there is any cohort)
    if (n_coho != nrow(df)){
      g_hat <- aggregate(x=g_hat, by=list(id=df$id), mean)[,2]
    }
  }
  
  # 1.3. update Q .................................
  H_AW <- df_coho$A/g_hat - (1-df_coho$A)/(1-g_hat)
  suppressWarnings({
    logitUpdate = glm(df_coho$Y ~ -1 + H_AW, 
                      offset = qlogis(QbarAW),
                      family = binomial())
    })
  epsilon <- logitUpdate$coef
  QbarAW_star = plogis(qlogis(QbarAW) + epsilon*H_AW)
  Qbar1W_star = plogis(qlogis(Qbar1W) + epsilon/g_hat)
  Qbar0W_star = plogis(qlogis(Qbar0W) - epsilon/(1-g_hat))
  
  TSM1 = mean(Qbar1W_star) 
  TSM0 = mean(Qbar0W_star)
  
  TSM1_init = mean(Qbar1W)
  TSM0_init = mean(Qbar0W)
  
  print(paste0("initial Qbar1W, Qbar0W: ", round(TSM1_init, 3), ", ", round(TSM0_init, 3)))
  
  print(paste0("updated Qbar1W, Qbar0W: ", round(TSM1, 3),", ", round(TSM0, 3)))
  
  RD_tmle = TSM1 - TSM0
  RR_tmle = TSM1 / TSM0
  logRR_tmle = log(RR_tmle)

  # Part 2. Inference .................................
  IC1 <- (df_coho$A/g_hat)*(df_coho$Y - QbarAW_star) + Qbar1W_star - TSM1
  IC0 <- ((1-df_coho$A)/(1-g_hat))*(df_coho$Y - QbarAW_star) + Qbar0W_star - TSM0
  
  # IC_RD <- H_AW*(df_coho$Y - QbarAW_star) + Qbar1W_star - Qbar0W_star - RD_tmle
  IC_RD <- IC1 - IC0
  IC_logRR <- IC1/TSM1 - IC0/TSM0
  
  # adjust for potential dependency 
  if (dependency){
    
    df_ic_RD <- data.frame('id' = df_coho$id,
                           'ic' = IC_RD)
    
    res_add_inf_RD <- adj_depen(df_ic = df_ic_RD, 
                                df_overlap = df_overlap, 
                                SL_lib = SL_lib_depen)
    
    df_ic_logRR <- data.frame('id' = df_coho$id,
                              'ic' = IC_logRR)
    
    res_add_inf_logRR <- adj_depen(df_ic = df_ic_logRR, 
                                df_overlap = df_overlap, 
                                SL_lib = SL_lib_depen)

  }else{
    res_add_inf_RD = 0
    res_add_inf_logRR = 0
  }
  
  var_IC_RD <- 1/(n_coho^2) * (n_coho*var(IC_RD) + res_add_inf_RD)
  var_IC_RD_unadj <- 1/(n_coho^2) * (n_coho*var(IC_RD))
  
  var_IC_logRR <- 1/(n_coho^2) * (n_coho*var(IC_logRR) + res_add_inf_logRR)
  var_IC_logRR_unadj <- 1/(n_coho^2) * (n_coho*var(IC_logRR))
  
  res_CI_RD <- get_CI(psi_hat = RD_tmle, var = var_IC_RD)
  res_CI_logRR <- get_CI(psi_hat = logRR_tmle, var = var_IC_logRR)
  res_CI_RR <- exp(res_CI_logRR)
  
  res_CI_RD_uadj <- get_CI(psi_hat = RD_tmle, var = var_IC_RD_unadj)
  res_CI_logRR_uadj <- get_CI(psi_hat = logRR_tmle, var = var_IC_logRR_unadj)
  res_CI_RR_uadj <- exp(res_CI_logRR_uadj)
  
  # transform RD back ## TEMP
  # print(TSM1)
  TSM1 = TSM1*(y_u - y_l) + y_l
  # print(TSM1)
  TSM0 = TSM0*(y_u - y_l) + y_l
  RD_tmle = TSM1 - TSM0
  res_CI_RD$CI_l <- (y_u - y_l)*res_CI_RD$CI_l
  res_CI_RD$CI_u <- (y_u - y_l)*res_CI_RD$CI_u
  res_CI_RD_uadj$CI_l = (y_u - y_l)*res_CI_RD_uadj$CI_l
  res_CI_RD_uadj$CI_u = (y_u - y_l)*res_CI_RD_uadj$CI_u

  # Part 3. Output .................................
  # save sl coefs
  SL_Q_coef <- res_fit_Q$Q_fit$coefficients
  SL_g_coef <- res_fit_g$g_fit$coefficients
  list_SL_coefs <- list("SL_Q_coef" = SL_Q_coef,
                        "SL_g_coef" = SL_g_coef)
  # save est
  df_res <- data.frame('Risk1' = TSM1, 
                       'Risk0' = TSM0, 
                       'Psi_hat' = c(RD_tmle, RR_tmle), 
                       'CI_l' = c(res_CI_RD$CI_l, res_CI_RR$CI_l), 
                       'CI_u' = c(res_CI_RD$CI_u, res_CI_RR$CI_u), 
                       'CI_l_unadj' = c(res_CI_RD_uadj$CI_l, res_CI_RR_uadj$CI_l), 
                       'CI_u_unadj' = c(res_CI_RD_uadj$CI_u, res_CI_RR_uadj$CI_u),
                       'Psi_type' = c("RD", "RR"), 
                       'Dependency' = dependency, 
                       'Qlevel' = Qlevel, 
                       'glevel' = glevel)
  
  res <- list("estimates" = df_res,
              "SL_coefs" = list_SL_coefs,
              "g_hat" = g_hat)
  
  if(any(prop.table(table(df$Y, df$A)) < 0.01) & 
     (mean(res$g_hat<0.2)>0.33 | mean(res$g_hat>0.8)>0.33)){
    print("WARNING: rare outcome; check for positivity violation.")
  }

  return(res)
}


fit_Q <- function(df, SL_lib_Q, QAdj, unadj_est = F, Vfolds){
  # now do initial estimation at the cluster level
  df1 <- df %>% mutate(A = 1)
  df0 <- df %>% mutate(A = 0)
  
  # use sl3
  sl3_task <- make_sl3_Task(
    data = df,
    covariates = c(QAdj,'A'),
    outcome = 'Y',
    weights = 'alpha',
    id = 'id',
    folds = origami::make_folds(df, fold_fun = folds_vfold, V = Vfolds)
  )
  
  # define leaner library
  sl <- make_learner(Lrnr_sl, 
                     learners = SL_lib_Q)

  sl_fit <- sl$train(sl3_task)
  
  sl3_task1 <- make_sl3_Task(
    data = df1,
    covariates = c(QAdj,'A'),
    outcome = 'Y',
    weights = 'alpha',
    id = 'id'
  )
  
  sl3_task0 <- make_sl3_Task(
    data = df0,
    covariates = c(QAdj,'A'),
    outcome = 'Y',
    weights = 'alpha',
    id = 'id'
  )
  
  if (unadj_est){
    # do not use cv predictions for unadjusted est
    QbarAW <- sl_fit$predict()
    Qbar1W <- sl_fit$predict(sl3_task1)
    Qbar0W <- sl_fit$predict(sl3_task0)
  }else{
    QbarAW <- sl_fit$predict_fold(sl3_task,"validation")
    Qbar1W <- sl_fit$predict_fold(sl3_task1,"validation")
    Qbar0W <- sl_fit$predict_fold(sl3_task0,"validation")
  }
  
  df_Q <- data.frame("QbarAW" = QbarAW,
                     "Qbar1W" = Qbar1W,
                     "Qbar0W" = Qbar0W)
  
  return(list("Q_fit" = sl_fit,
              "df_Q" = df_Q))
}

fit_g <- function(df, SL_lib_g, gAdj, Vfolds, unadj_est = F, EM = F){
  # use sl3
  sl3_task <- make_sl3_Task(
    data = df,
    covariates = gAdj,
    outcome = 'A',
    weights = 'alpha',
    id = 'id',
    outcome_type = 'binomial',
    folds = origami::make_folds(df, fold_fun = folds_vfold, V = Vfolds)
  )
  
  # define leaner library
  if(!EM) lb_metalearner <- make_learner(Lrnr_solnp,
                                 learner_function = metalearner_logistic_binomial,
                                 loss_function = loss_loglik_binomial)
  
  if(EM) lb_metalearner = Lrnr_cv_selector$new()
  
  
  sl <- make_learner(Lrnr_sl, 
                     learners = SL_lib_g, 
                     outcome_type = 'binomial',
                     metalearner = lb_metalearner)
  
  sl_fit <- sl$train(sl3_task)
  
  if (unadj_est){
    g_hat <- sl_fit$predict()
  }else{
    g_hat <- sl_fit$predict_fold(sl3_task,"validation")
  }
  
  # bound g
  g_hat[g_hat< 0.025] <- 0.025
  g_hat[g_hat> 0.975] <- 0.975
  
  return(list("g_fit" = sl_fit,
              "g_hat" = g_hat))
}

get_CI <- function(psi_hat, var){
  # standard error (square root of the variance)
  se <- sqrt(var)
  # 95% confidence interval
  CI_l <- psi_hat - 1.96*se
  CI_u <-  psi_hat + 1.96*se
  
  res <- data.frame(psi_hat, CI_l, CI_u)
  return(res) 
}

# ......................................................
# adj_depen: function to estimate correlation term with sl3
# ......................................................
# input: 
#     df_ic, a data frame with: id, ic_product

adj_depen <- function(df_ic, 
                      df_overlap, 
                      SL_lib){
  
  df_ic$id <- as.character(df_ic$id)
  df_overlap$coho_i <- as.character(df_overlap$coho_i)
  df_overlap$coho_j <- as.character(df_overlap$coho_j)
  
  # convert date to numeric and centralize it
  df_overlap$date_i <- as.numeric(df_overlap$date_i) - 
    mean(as.numeric(df_overlap$date_i))
  df_overlap$date_j <- as.numeric(df_overlap$date_j) - 
    mean(as.numeric(df_overlap$date_j))
  
  # 
  df_overlap <- df_overlap[which(df_overlap$coho_i %in% df_ic$id &
                                   df_overlap$coho_j %in% df_ic$id),]
  
  for (k in 1:nrow(df_overlap)){
    df_overlap$ic_product[k] <- df_ic$ic[df_ic$id == df_overlap$coho_i[k]] * 
      df_ic$ic[df_ic$id == df_overlap$coho_j[k]]
  }
  
  overlap_covarname <- c("dist_ij", "date_i", "date_j", "intarm_i", "intarm_j")
  
  # use sl3
  sl3_task <- make_sl3_Task(
    data = df_overlap,
    covariates = overlap_covarname,
    outcome = 'ic_product'
  )
  
  # define leaner library
  sl <- make_learner(Lrnr_sl, learners = SL_lib)
  sl_fit <- sl$train(sl3_task)
  ic_pred <- sl_fit$predict()
  
  return(2*sum(ic_pred))
}


# check the overlaps among cohorts
check_overlap_all <- function(df, r= 1000){
  
  df_coho <- df %>% filter(int_trig == 1)
  
  res_overlap <- data.frame('coho_i' = NA,
                            'coho_j' = NA,
                            'intarm_i' = NA,
                            'intarm_j' = NA,
                            'date_i' = NA,
                            'date_j' = NA,
                            'dist_ij' = NA,
                            'overlap_type' = NA)
  
  for (i in 1:(nrow(df_coho)-1)){
    this_coho <- df_coho[i,]
    this_ll <- this_coho %>% dplyr::select(latitude, longitude)
    this_arm <- this_coho$A
    this_date <- this_coho$date
    
    other_coho <- df_coho[-c(1:i),]
    other_ll <- other_coho %>% dplyr::select(latitude, longitude)
    
    id_product <- expand.grid(this_coho$id, other_coho$id)
    
    this_dist <- distGeo(this_ll, other_ll)
    
    df_overlap_loop <- data.frame('coho_i' = id_product[,1],
                                  'coho_j' = id_product[,2],
                                  'intarm_i' = this_arm,
                                  'intarm_j' = other_coho$A,
                                  'date_i' = this_date,
                                  'date_j' = other_coho$date,
                                  'dist_ij' = this_dist,
                                  'overlap_type' = NA,
                                  't0_i' = this_coho$t0,
                                  't1_i' = this_coho$t1,
                                  't0_j' = other_coho$t0,
                                  't1_j' = other_coho$t1)
    
    # filter by space
    df_overlap_loop <- df_overlap_loop %>% filter(dist_ij < 2*r)
    # filter by time
    df_overlap_loop <- df_overlap_loop %>% filter(!(t1_i < t0_j | 
                                                      t1_j < t0_i))
    
    if (nrow(df_overlap_loop) == 0){
      next
    }
    
    df_overlap_loop <- df_overlap_loop %>% 
      mutate(overlap_type = ifelse(dist_ij <= 2*500, 
                                   'target',
                                   'spill'))
    
    df_overlap_loop <- df_overlap_loop[,c(1:8)]
    res_overlap <- bind_rows(res_overlap, df_overlap_loop)
  }
  # remove the NA row
  res_overlap <- res_overlap[-1,]
  return(res_overlap)
}
