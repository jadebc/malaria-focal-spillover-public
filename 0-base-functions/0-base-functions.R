
##############################################
##############################################

# Documentation: format_glm
# Usage: format_glm(data, rfit, outcome, treatment, year, adj, did, family)
# Description: nicely format glm model output
# Args/Options: 
# data:                    the data used to fit GLM
# cluster_coefficients:    coefficients and SEs adjusted for clustering
# family:                  option for the type of GLM family being used
# 
# Returns: an appropriately formatted tibble that includes: sample size, a point estimate, 
# robust standard errors, the lower and upper bounds of a 95% CI, the outcome, the treatment, 
# and the year for which the glm was fitted

# Output: none

format_glm = function(data,
                      cluster_coefficients,
                      family = "gaussian") {
  
  row_index = ifelse(test = did,
                     yes = which(names(cluster_coefficients[, 1]) == "did"),
                     no = 2)
  
  n      = nrow(data)
  pt.est = cluster_coefficients[row_index, 1]
  se     = cluster_coefficients[row_index, 2]
  
  lb     = pt.est - (qnorm(0.975) * se)
  ub     = pt.est + (qnorm(0.975) * se)
  
  if (family == "binomial") {
    formatted_glm_results = tibble(
      "Sample Size" = n,
      "Point Estimate" = pt.est,
      "Standard Error" = se,
      "95% CI Lower Bound" = lb,
      "95% CI Upper Bound" = ub
    ) %>% 
      exp %>%
      select(Outcome, Year, `Sample Size`, `Point Estimate`, `Robust Standard Error`,
             `95% CI Lower Bound`, `95% CI Upper Bound`)
    
  } else {
    formatted_glm_results = tibble(
      "Outcome" = outcome,
      "Year" = year,
      "Sample Size" = n,
      "Point Estimate" = pt.est,
      "Standard Error" = se,
      "95% CI Lower Bound" = lb,
      "95% CI Upper Bound" = ub
    ) 
    
  }
  
  return(formatted_glm_results)
}





sandwichSE=function (fm, cluster) 
{
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj <- apply(estfun(fm), 2, function(x) tapply(x, cluster, 
                                                sum))
  vcovCL <- dfc * sandwich(fm, meat = crossprod(uj)/N)
  return(vcovCL)
}

# vcovCL <- sandwichSE(fm = fit, cluster = mudat$id)
# rfit <- coeftest(fit, vcovCL)
# lb <- rfit[1, 1] - 1.96 * rfit[1, 2]
# ub <- rfit[1, 1] + 1.96 * rfit[1, 2]



get_ci = function(fit, delta_formula, pr, clustervec){
  
  vcovCL <- sandwichSE(fm = fit, cluster = clustervec)
  
  se = msm::deltamethod(g = delta_formula, coef(fit), vcovCL)
  lb = exp(log(pr) - se*qnorm(0.975))
  ub = exp(log(pr) + se*qnorm(0.975))
  
  out = data.frame(
    pr = round(pr, 2),
    se = round(se, 2),
    lb = round(lb, 2),
    ub = round(ub, 2)
  )
  
  print(out)
  
  return(out)
}

g_comp = function(data1, data0, fit){
  
  
}




mean_se = function(Y, id, print = TRUE){
  mudat <- data.frame(id = id, Y = Y)
  n.orig <- dim(mudat)[1]
  mudat <- mudat[complete.cases(mudat), ]
  n.sub <- dim(mudat)[1]
  if (n.orig > n.sub) 
    cat("\n-----------------------------------------\nDropping", 
        n.orig - n.sub, "observations\ndue to missing values in the outcome\n", 
        "Final sample size:", n.sub, "\n-----------------------------------------\n")
  fit <- glm(Y ~ 1, family = binomial, data = mudat)
  vcovCL <- washb::sandwichSE(dat = mudat, fm = fit, cluster = mudat$id)
  rfit <- coeftest(fit, vcovCL)
  lb <- rfit[1, 1] - 1.96 * rfit[1, 2]
  ub <- rfit[1, 1] + 1.96 * rfit[1, 2]
  mean_ci <- matrix(c(n.sub, exp(rfit[1, 1]),  exp(lb), exp(ub)), nrow = 1, ncol = 4)
  colnames(mean_ci) <- c("N", "Prev", "lb", "ub")
  
  return(mean_ci)
}



#-----------------------------------------------
# Helper functions
#-----------------------------------------------


# Helper function 1
#--------------------
# A helper function, try the alternative directory 
# if the default fails when importing data.
# Input: read_func: a importing function such as read.csv
#        path: the default path
# Output: df: a data set
read_df <- function (read_func, path) {
  df <- try(suppressWarnings(read_func(path)))
  
  if (is.null(df) | class(df) == "try-error") {
    cat("Defaul path not found, trying the alternative")
    path = str_remove(path, "/Datasets")
    df <- read_func(path)
  }
  
  return(df)
}

read_shp <- function(dsn, layer){
  df <- try(suppressWarnings(rgdal::readOGR(dsn = dsn,
                                            layer = layer)))
  if (is.null(df) | class(df) == "try-error") {
    cat("Defaul dsn not found, trying the alternative")
    dsn = str_remove(dsn, "/Datasets")
    df <- rgdal::readOGR(dsn = dsn, layer = layer)
  }
  return(df)
}


# Helper function 2
#--------------------
# helper functions for wash-out check

get_last_day <- function(date_list, date_x){
  last_day <- NULL
  last_position <- 0
  last_diff <- 0
  for (i in 1:length(date_list)){
    if (date_list[i] - date_x > 0){
      break
    }
    last_day <- date_list[i]
    last_position <- last_position + 1
  }
  if (!is.null(last_day)){
    last_diff <- as.numeric(date_x - last_day)
  }
  return(list("last_day" = last_day,
              "last_position" = last_position,
              "last_diff" = last_diff))
}

get_next_day <- function(date_list, date_x){
  next_day <- NULL
  next_position <- 0
  next_diff <- 0
  for (i in 1:length(date_list)){
    next_position <- next_position + 1
    if (date_list[i] - date_x > 0){
      next_day <- date_list[i]
      break
    }
  }
  if (!is.null(next_day)){
    next_diff <- as.numeric(next_day - date_x)
  }
  return(list("next_day" = next_day,
              "next_position" = next_position,
              "next_diff" = next_diff))
  
}


check_washout_arm <- function(prev_match_list, this_date, mindays0 = 35, mindays1 = 56){
  
  check_prev = TRUE
  check_next = TRUE
  
  prev_days <- prev_match_list$date
  last_day_res <- get_last_day(prev_days, this_date)
  next_day_res <- get_next_day(prev_days, this_date)
  
  if (!is.null(last_day_res$last_day)){
    last_type <- prev_match_list[last_day_res$last_position, ]$indexcase + 
      
      prev_match_list[last_day_res$last_position, ]$int_trig
    
    # last_type can be 2 (index_int), 1 (index_nonint), 0 (int_recip)
    last_type <- ifelse(is.na(last_type),0,last_type)
    diff_days <- last_day_res$last_diff
    # if last match is within the wash-out period
    if ((this_arm == 1 & diff_days <= mindays1) |
        (this_arm == 0 & diff_days <= mindays0)){
      check_prev = FALSE
    }
    # if last match is index_nonint, ignore check
    if (last_type == 1){
      check_prev = TRUE
    }
  }
  
  if (!is.null(next_day_res$next_day)){
    next_type <- prev_match_list[next_day_res$next_position, ]$indexcase + 
      
      prev_match_list[next_day_res$next_position, ]$int_trig
    
    # next_type can be 2 (index_int), 1 (index_nonint), 0 (int_recip)
    next_type <- ifelse(is.na(next_type),0,next_type)
    diff_days <- next_day_res$next_diff
    # if last match is within the wash-out period
    if ((this_arm == 1 & diff_days <= mindays1) |
        (this_arm == 0 & diff_days <= mindays0)){
      check_next = FALSE
    }
    # if last match is index_nonint, ignore check
    if (next_type == 1){
      check_next = TRUE
    }
  }
  washout_check = all(check_prev, check_next)
  return(washout_check)
}


check_washout <- function(prev_match_list, this_date, mindays = 21){
  
  check_prev = TRUE
  check_next = TRUE
  
  prev_days <- prev_match_list$date
  last_day_res <- get_last_day(prev_days, this_date)
  next_day_res <- get_next_day(prev_days, this_date)
  
  if (!is.null(last_day_res$last_day)){
    last_type <- prev_match_list[last_day_res$last_position, ]$int_day
    last_type <- ifelse(is.na(last_type),0,last_type)
    diff_days <- last_day_res$last_diff
    # if last match is index match and within the 3 weeks wash-out period
    if (last_type != 1 & diff_days <= mindays){
      check_prev = FALSE
    }
  }
  
  if (!is.null(next_day_res$next_day)){
    next_type <- prev_match_list[next_day_res$next_position, ]$int_day
    next_type <- ifelse(is.na(next_type),0,next_type)
    diff_days <- next_day_res$next_diff
    # if last match is index match and within the 3 weeks wash-out period
    if (next_type != 1 & diff_days <= mindays){
      check_next = FALSE
    }
  }
  
  washout_check = all(check_prev, check_next)
  return(washout_check)
}








# Note that "target" overlap means purely 
# "target zone overlaps with target zone", 
# all other kind of overlaps was classified as "spill") 

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



positivity_screening = function(data, varlist){
  
  Wdata = data %>% ungroup() %>% dplyr::select(all_of(varlist))
  
  # identify categorical variables ------------------
  is.categorical = function(data, varname){
    x = data %>% dplyr::select(!!sym(varname)) %>% pull()
    factor = is.factor(x)
    twocats = length(unique(x))<=2
    return(ifelse(factor | twocats, T, F))
  }
  categorical_vars = varlist[map_lgl(as.list(varlist), function(x) is.categorical(data = data, varname = x))]
  
  # identify which categorical variables have positivity violations
  cat_violations = lapply(as.list(categorical_vars), function(x)
    check_positivity_categorical(data = data,
                                 varname = x)) %>% unlist()
  
  # drop covariates with violation from list
  if(all(!cat_violations)){
    keep_cat = categorical_vars
  }else{
    keep_cat = categorical_vars[-which(categorical_vars%in% categorical_vars[cat_violations])]
  }
  drop_cat = categorical_vars[which(categorical_vars%in% categorical_vars[cat_violations])]
  assert_that(length(categorical_vars)==length(c(keep_cat, drop_cat)))
  
  # identify continuous variables ------------------
  continuous_vars = varlist[!varlist %in% categorical_vars]
  
  # identify which continuous variables have positivity violations
  cont_violations = lapply(as.list(continuous_vars), function(x)
    check_positivity_continuous(data = data,
                                varname = x)) %>% unlist()
  cont_violation_vars = continuous_vars[which(continuous_vars%in% continuous_vars[cont_violations])]
  cont_noviol_vars = continuous_vars[-which(continuous_vars%in% continuous_vars[cont_violations])]
  assert_that(length(continuous_vars)==length(c(cont_violation_vars, cont_noviol_vars)))
  
  # recode continuous variables with positivity violations
  newdata_list = list()
  data = as.data.frame(data)
  for(i in 1:length(cont_violation_vars)){
    newdata_list[[i]] = cont_to_cat(data[,cont_violation_vars[i]])
  }
  names(newdata_list) = paste0(cont_violation_vars, "_cat")
  newdf = as.data.frame(bind_cols(newdata_list))
  
  newdata = bind_cols(
    data %>% dplyr::select(-all_of(varlist)),
    data %>% dplyr::select(all_of(keep_cat)),
    data %>% dplyr::select(all_of(cont_noviol_vars)),
    newdf
  )
  
  # recheck positivity after continuous recoded to categorical 
  # identify categorical variables ------------------
  varlist2 = colnames(newdf)
  categorical_vars2 = varlist2[map_lgl(as.list(varlist2), function(x) is.categorical(data = newdata, varname = x))]
  
  cat_violations2 = lapply(as.list(categorical_vars2), function(x)
    check_positivity_categorical(data = newdata,
                                 varname = x)) %>% unlist()
  
  # drop covariates with violation from list
  if(all(!cat_violations2)){
    keep_cat2 = categorical_vars2
  }else{
    keep_cat2 = categorical_vars2[-which(categorical_vars2%in% categorical_vars2[cat_violations2])]
  }
  drop_cat2 = categorical_vars2[which(categorical_vars2%in% categorical_vars2[cat_violations2])]
  assert_that(length(categorical_vars2)==length(c(keep_cat2, drop_cat2)))
  
  finaldata = newdata %>% dplyr::select(-all_of(drop_cat2))
  finalcovars = c(keep_cat, cont_noviol_vars, keep_cat2)
  
  return(list(data = finaldata,
              covariates = finalcovars))
}





check_positivity_categorical = function(data, varname){
  x = data %>% group_by(A, !!sym(varname)) %>% 
    summarise(n=n())
  if(min(x$n)<20){
    return(T)
  }else{
    return(F)
  } 
}



check_positivity_continuous = function(data, varname){
  dvar = data %>% dplyr::select(!!sym(varname)) %>% pull()
  if(length(unique(dvar))>=20){
    mybreaks = seq(min(data %>% pull(!!sym(varname))), max(data %>% pull(!!sym(varname))), length=20)
    d = data %>% mutate(varcut = cut(!!sym(varname), 
                                     breaks = mybreaks))
  }else{
    d = data %>% mutate(varcut = !!sym(varname)) %>% 
      mutate(as.factor(varcut))
  }
  
  
  x = as.data.frame(table(d$varcut, d$A))
  
  if(min(x$Freq)==0){
    return(T)
  }else{
    return(F)
  }
}

# convert continuous to categorical
cont_to_cat <- function(x){
  quartiles = quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  tertiles = quantile(x, probs = c(0, 0.33, 0.67, 1))
  median = quantile(x, probs = c(0, 0.5, 1), digits = 30)
  
  if(length(unique(quartiles))==5){
    out = cut(x, quartiles, include.lowest = T)
    return(out)
  }  
  if(length(unique(quartiles))!=5 & length(unique(tertiles))==4){
    out = cut(x, tertiles, include.lowest = T)
    return(out)
  } 
  if(length(unique(tertiles))!=4 & length(unique(median))==3){
    out = cut(x, median, include.lowest = T)
    return(out)
  }
  if(length(unique(median))!=3){
    if(median[1]==0 & median[2]==0){
      out = as.factor(ifelse(x ==0, "zero", "non-zero"))
      return(out)
    }else{
      if(median[2]==median[3]){
        out = as.factor(ifelse(x==median[2], "median","non-median"))
        return(out)
      }else{
        out = rep(NA, length = length(x))
        return(out)
      }
    }
  }
  
}



# function to read and process files --------------------------------------------------------
# both_obs: T if the results file is a list with both the short and long observation period 
#           F if only one type was saved (the one that matches reservoir)
read_file_results <- function(file_path, both_obs=T){
  
  results <- readRDS(paste0(results_path, file_path))
  
  if(length(grep("_human_", file_path)==1)) reservoir = "Human"
  if(length(grep("_mosq_", file_path)==1)) reservoir = "Mosquito"
  if(length(grep("_hm_", file_path)==1)) reservoir = "Human & mosquito"
  
  if(length(grep("spillover", file_path)==1)) parameter = "spillover"
  if(length(grep("total", file_path)==1)) parameter = "total"
  if(length(grep("direct", file_path)==1)) parameter = "direct"
  
  
  if(both_obs){
    
    if(length(grep("short", file_path)==1)) period = "5 weeks"
    if(length(grep("long", file_path)==1)) period = "6 months"
    
    if(all(!is.na(results$Long$res_est))){
      estimates_long <- results$Long$res_est$estimates
    }else{
      print(paste0("File '", file_path, "' contains no results."))
      
      estimates_long <- data.frame(
        Risk1 = NA, 
        Risk0 = NA,
        Psi_hat = NA,
        CI_l = NA,
        CI_u = NA,
        CI_l_unadj = NA,
        CI_u_unadj = NA, 
        Psi_type = NA,
        Dependency = NA,
        Qlevel = NA,
        glevel = NA, 
        reservoir = reservoir,
        parameter = parameter
      )
    }
    
    estimates_short <- results$Short$res_est$estimates
    
    estimates <- bind_rows(
      estimates_long %>% mutate(period = "6 months"),
      estimates_short %>% mutate(period = "5 weeks")
    )
    
  }
  if(!both_obs){
    estimates <- results$res_est$estimates 
    if(reservoir=="Human") estimates <- estimates %>% mutate(period = "5 weeks")
    if(reservoir!="Human") estimates <- estimates %>% mutate(period = "6 months")
  }

  return(estimates)
}


filter_data <- function(data, modifier_name){
  data <- data %>% rename(modifier =  !!sym(modifier_name))
  modifier1 = data %>% filter(modifier==1) %>% 
    mutate(modifier_var = modifier_name,
           modifier = 1)
  modifier0 = data %>% filter(modifier==0) %>% 
    mutate(modifier_var = modifier_name,
           modifier = 0)
  return(list(modifier1 = modifier1, modifier0 = modifier0))
}
