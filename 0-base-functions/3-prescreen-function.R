# pre-screening with cov_prescreen
# source code (https://github.com/ben-arnold/washb/blob/master/R/washb_prescreen.R)

cov_prescreen <- function(df, 
                            Y, 
                            Ws,
                            yname,
                            covarname,
                            family="binomial", 
                            pval=0.2, 
                            print=TRUE) {
  # Y   : outcome variable of interest
  # Ws  : data frame of candidate covariates to screen
  # family : exponential model family (gaussian for continuous outcomes, binomial for binary outcomes, poisson for counts, and neg.binom for negative binomial models)
  require(lmtest)
  if(family[[1]]=="neg.binom"){
    require(MASS)
  }
  
  #Check pvalue
  if(pval>1|pval<0){
    stop("P-value threshold not set between 0 and 1.")
  }
  
  # ensure Ws are a data frame
  Ws <- as.data.frame(Ws)
  
  dat <- data.frame(Ws,Y)
  dat <- dat[complete.cases(dat),]
  nW <- ncol(Ws)
  LRp <- matrix(rep(NA,nW),nrow=nW,ncol=1)
  rownames(LRp) <- names(Ws)
  colnames(LRp) <- "P-value"
  if(family[[1]]!="neg.binom"){
    for(i in 1:nW) {
      dat$W <- dat[,i]
      if(class(dat$W)=="factor") dat$W = droplevels(dat$W)
      if(class(dat$W)=="factor" & dim(table(dat$W))==1){
        #skip factors with a single level to avoid error
        fit1 <- fit0 <- glm(Y~1,data=dat,family=family)
      }else{
        fit1 <- glm(Y~W,data=dat,family=family)
        fit0 <- glm(Y~1,data=dat,family=family)
      }
      LRp[i] <- lrtest(fit1,fit0)[2,5]
    }
  }else{
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Pkg needed for this function to work. Please install it.",
           call. = FALSE)
    }else{
      for(i in 1:nW) {
        dat$W <- dat[,i]
        if(class(dat$W)=="factor" & dim(table(dat$W))==1){
          #skip factors with a single level to avoid error
          fit1 <- fit0 <- glm(Y~1,data=dat,family=family)
        }else{
          fit1 <- glm.nb(Y~W,data=dat,family=family)
          fit0 <- glm.nb(Y~1,data=dat,family=family)
        }
        LRp[i] <- lrtest(fit1,fit0)[2,5]
      }
    }
  }
  p20 <- ifelse(LRp<pval,1,0)
  
  if(print==TRUE){
    cat("\nLikelihood Ratio Test P-values:\n")
    print(round(LRp,5))
    if(sum(p20)>0) {
      LRps <- matrix(LRp[p20==1,],ncol=1)
      rownames(LRps) <- names(Ws)[p20==1]
      colnames(LRps) <- "P-value"
      cat(paste("\n\nCovariates selected (P<",pval,"):\n",sep=""))
      print(LRps)
    } else{
      cat(paste("\nNo covariates were associated with the outcome at P<",pval))
    }
  }
  
  var_screen <-names(Ws)[p20==1]
  
  # print the names of variable excluded from lmtest
  print("===== variables excluded after lmtest ===== ")
  print(setdiff(covarname, var_screen))
  
  # check variable sparsity
  df_covar <- df %>% dplyr::select(all_of(covarname))
  temp <- sapply(df_covar,
                 function(x) length(x[!is.na(x) & x!=0])/length(x))
  
  var_nonsparse <- names(df_covar)[temp >= 0.05]
  
  # print the names of variable excluded from sparsity check
  print("===== variables excluded after sparsity check ===== ")
  print(setdiff(covarname, var_nonsparse))
  
  # check constant variable 
  temp <- sapply(df_covar,
                 function(x) length(unique(x))!=1)
  
  var_nonconstant <- names(df_covar)[temp]
  
  # print the names of constant variable
  print("===== excluded constant variables ===== ")
  print(setdiff(covarname, var_nonconstant))
  
  # selected variables
  var_selected <- intersect(intersect(var_nonsparse, var_screen),var_nonconstant)
  
  # print the names of variable left
  print("===== variables selected in the end ===== ")
  print(var_selected)
  
  if(family=="binomial"){
    # check EPV (number of events per variable)
    if(length(var_selected)>0){
      p = length(var_selected)
      df_Y <- df %>% dplyr::select(all_of(yname))
      epv = length(nrow(df_Y)[df_Y != 0])/p
      
      if(epv<10) print("Warning: number of events per variable < 10")
      
    }

  }


  
  return(var_selected)
}

library(sl3)
library(tmle3)
# Y   : outcome variable of interest
# Ws  : data frame of candidate covariates to screen
lasso_prescreen <- function(df, yname,covarname, family="gaussian") {

  
  # define lasso learner
  lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
  screen_lasso <- Lrnr_screener_coefs$new(learner = lrn_lasso, threshold = 0)

  
  # create the sl3 task
  screen_task <- make_sl3_Task(
    data = df,
    covariates = covarname,
    outcome = yname
  )
  
  res_screen <- screen_lasso$train(screen_task)
  
  var_screen <-  res_screen$fit_object$selected
  
  
  # print the names of variable excluded from lasso
  print("===== variables excluded from lasso ===== ")
  print(setdiff(covarname, var_screen))
  
  # check variable sparsity
  df_covar <- df %>% dplyr::select(all_of(covarname))
  temp <- sapply(df_covar,
                 function(x) length(x[!is.na(x) & x!=0])/length(x))
  
  var_nonsparse <- names(df_covar)[temp >= 0.05]
  
  # print the names of variable excluded from sparsity check
  print("===== variables excluded from sparsity check ===== ")
  print(setdiff(covarname, var_nonsparse))
  
  # selected variables
  var_selected <- intersect(var_nonsparse, var_screen)
  
  # print the names of variable left
  print("===== variables selected in the end ===== ")
  print(var_selected)
  
  # check EPV (number of events per variable)
  p = length(var_selected)
  df_Y <- df %>% dplyr::select(all_of(yname))
  epv = length(nrow(df_Y)[df_Y != 0])/p
  
  if (epv < 10){
    print("Warning: number of events per variable is less than 10")
  }
  
  return(var_selected)
  
}


