# combine distance variables to avoid colinearity and
# impute missings

preprocess_data <- function(data){
  # combine dist variable to avoid colinearity
  df_temp <- comb_distvar(df_all = data)
  
  # drop infected_date, which we do not need for estimation
  df_temp <- df_temp %>% dplyr::select(-infected_date)
  
  # prepare data for tmle ---------------------------------------
  # impute missingness (age)
  nodes <- list(W = colnames(df_temp)[-which(names(df_temp) %in% c("indexcase", "intarm"))],
                A = "intarm",
                Y = "indexcase")
  
  df_nomiss <- process_missing(data = df_temp,
                               node_list = nodes)$data
  
  df_nomiss <- as.data.frame(df_nomiss)
  
  df_nomiss$cohort_id = as.factor(df_nomiss$cohort_id)
  
  return(df_nomiss)
}



process_covariates <- function(data, dropnames){
  # specify the variables used for identification, but not for estimation
  drop_EM <- colnames(data)[grep("abovemed", colnames(data))]
  
  dropme <- c("indiv_id", "date", "clinic_date_index", "longitude", "latitude",
              "target_area", "int_recip", "infected_date", "t0", "t1", 
              "tte", "ea", "int_trig", "indexcase", "intarm", "cohort_id",
              "dist_index_quant",
              drop_EM)
  
  covarname_initial <- setdiff(names(data), dropme)
  
  # drop covariates to avoid collinearity
  covarname <- setdiff(covarname_initial, dropnames)
  return(covarname)
}




## format data for reservoir ---------------------------------------
prep_htmle_data <- function(data, reservoir){
  
  assert_that(reservoir %in% c("human", "mosquito", "human & mosquito"))
  
  # binarize treatment variable
  # human reservoir
  if (reservoir == "human") {
    df <- data %>%
      mutate(             
        ## add a column for constructing interaction term
        inter = as.numeric(intarm %in% c("RV", "TV")),
        intarm = case_when(intarm %in% c("TO", "TV") ~ 1,
                           intarm %in% c("RO", "RV") ~ 0)
      )
  }
  
  # mosquito reservoir
  if(reservoir == "mosquito"){
    df <- data %>%
      mutate(
        ## add a column for constructing interaction term
        inter = as.numeric(intarm %in% c("TO", "TV")),
        intarm = case_when(intarm %in% c("RV", "TV") ~ 1,
                           intarm %in% c("RO", "TO") ~ 0)
      )
  }
  
  # human & mosquito reservoir
  if(reservoir == "human & mosquito"){
    df <- data %>% filter(intarm %in% c("TV", "RO"))
    df <- df %>% mutate(intarm = case_when(intarm %in% c("TV") ~ 1,
                                           intarm %in% c("RO") ~ 0))
    
  }
  
  # rename variables for TMLE
  df <- df %>% rename(A = intarm, Y = indexcase, id = cohort_id)
  
  # add weights alpha = 1/N_j where N_j is the size of j-th cohort
  df_weights <- df %>% group_by(id) %>% summarise(alpha = 1/n())
  df <- left_join(df, df_weights, by = 'id')
  
  # recode sex to avoid sparsity issues
  df = df %>% mutate(sex = recode_factor(sex, `unknown`= "Male"))
  
  return(df)
}


# get upper bound of incidence
# (empirical maximum of cohort-level incidence
# x 150% [default])
get_y_u <- function(df, effecttype, scale=1.5){
  
  # subset data by parameter
  if (effecttype == 'direct'){
    df <- df %>% filter(target_area == 1, int_recip == 1)
  }else if (effecttype == 'spillover'){
    df <- df %>% filter(target_area == 0 | (target_area == 1 & int_recip == 0))
  }else if (effecttype == 'direct_prev'){
    df <- df %>% filter(target_area == 1)
  }
  
  # calculate cluster incidence
  inc <- df %>%
    group_by(id) %>%
    summarise(ncases = sum(Y),
              Npop = n()) %>%
    mutate(inc = ncases / Npop)
  
  # empirical max of cluster incidence * 150%
  y_u <- round(max(inc$inc)*scale,2)
  
  return(y_u)
  
}

# check that correlation between each covariate is less
# than a certain threshold
check_cor <- function(data, covarname, threshold = 0.9){
  x = data %>% dplyr::select(all_of(covarname_long)) %>% 
    dplyr::select(-sex)
  cormat <- cor(x)
  diag(cormat)=0
  assert_that(max(abs(cormat))<threshold, msg = paste0("Some covariates have >. ",
                                                       threshold*100, "% correlation!"))
}

