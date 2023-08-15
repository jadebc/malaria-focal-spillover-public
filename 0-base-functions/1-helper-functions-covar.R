

##############################################
##############################################

# Documentation: count_int_pre_recei
# Usage: count_int_pre_recei(cohorts)
# Description: Number of interventions* an individual previously received 
# Args/Options: 
# data: the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable
count_int_pre_recei <- function(cohorts, df_base){
  
  sub_df_base <- df_base %>% 
                 arrange(indiv_id, cohort_id, date) %>% 
                 dplyr::select(c("indiv_id", "cohort_id", "date"))
  
  person_list <- unique(sub_df_base$indiv_id)
  res <- sub_df_base %>% 
         mutate(n_int_pre_recei_RO = 0,
                n_int_pre_recei_RV = 0,
                n_int_pre_recei_TO = 0,
                n_int_pre_recei_TV = 0)
  
  # subset cohorts data to people with int match (n = 8612)
  subdf <- cohorts %>% filter(int_recip == 1)
  subdf <- dummy_cols(subdf, select_columns = c("intarm"))
  
  all_res <- map_dfr(person_list, function(x) 
    hp_pre_recei(res, subdf, x))
  
  return(all_res)
}


hp_pre_recei <- function(res, subdf, this_person){
  this_res <- res %>% filter(indiv_id == this_person)
  for (i in 1:nrow(this_res)){
    this_subdf <- subdf %>% filter(indiv_id == this_person & date <= this_res$date[i])
    if (nrow(this_subdf) != 0) {
      this_subdf <- this_subdf %>% 
        group_by(indiv_id) %>% 
        summarise("n_int_pre_recei_RO" = sum(intarm_RO),
                  "n_int_pre_recei_RV" = sum(intarm_RV),
                  "n_int_pre_recei_TO" = sum(intarm_TO),
                  "n_int_pre_recei_TV" = sum(intarm_TV)) %>% 
        dplyr::select(c("indiv_id", "n_int_pre_recei_RO",
                        "n_int_pre_recei_RV", "n_int_pre_recei_TO",
                        "n_int_pre_recei_TV")) 
      
      this_res$n_int_pre_recei_RO[i] = this_subdf$n_int_pre_recei_RO
      this_res$n_int_pre_recei_RV[i] = this_subdf$n_int_pre_recei_RV
      this_res$n_int_pre_recei_TO[i] = this_subdf$n_int_pre_recei_TO
      this_res$n_int_pre_recei_TV[i] = this_subdf$n_int_pre_recei_TV
    }
  }
  return(this_res)
}


##############################################
##############################################
# Documentation: count_int_pre_deliv
# Usage: count_int_pre_deliv(cohorts, r)
# Description: Number of interventions* previously delivered near each individual
# Args/Options: cohorts data and a sequence of radius
# data:  the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable

count_int_pre_deliv <- function(r = c(500, 1000, 2000, 3000), 
                                index_int, 
                                int_recip,
                                df_base,
                                list_people_coho){
  
  res_pre_int <- data.frame()
  for (i in 1:nrow(list_people_coho)){
    
    print(paste0("---------", round(i/nrow(list_people_coho), 5)*100,"%","---------"))
    
    this_person <- list_people_coho[i,]
    this_ll <- list_people_coho[i,c(2,3)]
    
    cohort_ll <- index_int %>% dplyr::select(c("longitude", "latitude"))
    dist <-  distGeo(this_ll, cohort_ll)
    
    df_base_this <- df_base %>% filter(indiv_id == this_person$indiv_id)
    
    for (k in 1:length(r)){
      
      cohort_nearby <- index_int[which(dist <= r[k]),]
      
      if (nrow(cohort_nearby) == 0){
        df_base_this$pre_int_RO <- 0
        df_base_this$pre_int_RV <- 0
        df_base_this$pre_int_TO <- 0
        df_base_this$pre_int_TV <- 0
        
        names(df_base_this)[(ncol(df_base_this)-3):ncol(df_base_this)] <- 
          sapply(names(df_base_this)[(ncol(df_base_this)-3):ncol(df_base_this)],
                 function(x) paste0(x, "_", r[k]))
      }else{
        add_coho_df <- data.frame("indiv_id" = rep(this_person$indiv_id, 
                                                   nrow(cohort_nearby)), 
                                  "cohort_id" = cohort_nearby$iid_combined, 
                                  "date" = cohort_nearby$int_date,
                                  "intarm" = cohort_nearby$intarm,
                                  "pre_int_RO" = 0,
                                  "pre_int_RV" = 0,
                                  "pre_int_TO" = 0,
                                  "pre_int_TV" = 0)
        
        
        names(add_coho_df)[5:8] <- sapply(names(add_coho_df)[5:8],
                                           function(x) paste0(x, "_", r[k]))
        
        # check prev int
        if (nrow(add_coho_df) == 1){
          add_coho_df[1,c(5:8)] = 0
        }else{
          for (j in 1:nrow(add_coho_df)){
            this_coho_int_day <- add_coho_df$date[j]
            other_coho <- add_coho_df[-j,]
            prev_coho <- other_coho[which(other_coho$date <= this_coho_int_day), ]
            if (nrow(prev_coho) == 0){
              add_coho_df[j,c(5:8)] = 0
            }else{
              # count individual
              pre_RO <- vector()
              pre_RV <- vector()
              pre_TO <- vector()
              pre_TV <- vector()
              for (m in 1:nrow(prev_coho)){
                pre_int_recip <- int_recip %>% filter(iid_combined == prev_coho$cohort_id[m])
                if (nrow(pre_int_recip) == 0){
                  pre_RO[m] <- 0
                  pre_RV[m] <- 0
                  pre_TO[m] <- 0
                  pre_TV[m] <- 0
                }else{
                  pre_int_recip_ll <- pre_int_recip %>% dplyr::select(c("longitude", "latitude"))
                  dist_individual <-  distGeo(this_ll, pre_int_recip_ll)
                  near_int_recip <- pre_int_recip[which(dist_individual <= r[k]),]
                  pre_RO[m] <- sum(near_int_recip$intarm == "RO")
                  pre_RV[m] <- sum(near_int_recip$intarm == "RV")
                  pre_TO[m] <- sum(near_int_recip$intarm == "TO")
                  pre_TV[m] <- sum(near_int_recip$intarm == "TV")
                }
              }
              add_coho_df[j,5] = sum(pre_RO)
              add_coho_df[j,6] = sum(pre_RV)
              add_coho_df[j,7] = sum(pre_TO)
              add_coho_df[j,8] = sum(pre_TV)
            }
            # print(overlap_coho$cohort_id)
          }
        }
        add_coho_df <- add_coho_df[, c(1,2,5:8)]
        df_base_this <- left_join(df_base_this, add_coho_df, by = c("indiv_id", "cohort_id"))
      }
    }
    res_pre_int <- bind_rows(res_pre_int, df_base_this)
  }
  return(res_pre_int)
}




##############################################
##############################################
# Documentation: count_int_concu_deliv
# Usage: count_int_concu_deliv(cohorts, r, t)
# Description: Number of interventions* an individual's neighbors 
#              "concurrently"(during the obs period) received 
# Args/Options: cohorts data, a sequence of radius r and a scalar value t for time range
# data:  the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable

count_int_concu_deliv <- function(r = c(500, 1000, 2000, 3000), 
                                  index_int, 
                                  int_recip,
                                  df_base,
                                  list_people_coho,
                                  t1_l = 0, 
                                  t1_u =35,
                                  t2_l = 0, 
                                  t2_u =180){
  
  res_concur_int <- data.frame()
  for (i in 1:nrow(list_people_coho)){
    
    print(paste0("---------", round(i/nrow(list_people_coho), 5)*100,"%","---------"))
    
    this_person <- list_people_coho[i,]
    this_ll <- list_people_coho[i,c(2,3)]
    
    cohort_ll <- index_int %>% dplyr::select(c("longitude", "latitude"))
    dist <-  distGeo(this_ll, cohort_ll)
    
    df_base_this <- df_base %>% filter(indiv_id == this_person$indiv_id)
    
    for (k in 1:length(r)){
      
      cohort_nearby <- index_int[which(dist <= r[k]),]
      
      if (nrow(cohort_nearby) == 0){
        df_base_this$concur_int_RO <- 0
        df_base_this$concur_int_RV <- 0
        df_base_this$concur_int_TO <- 0
        df_base_this$concur_int_TV <- 0
        
        names(df_base_this)[(ncol(df_base_this)-3):ncol(df_base_this)] <- 
              sapply(names(df_base_this)[(ncol(df_base_this)-3):ncol(df_base_this)],
                     function(x) paste0(x, "_", r[k]))
      }else{
        add_coho_df <- data.frame("indiv_id" = rep(this_person$indiv_id, 
                                                   nrow(cohort_nearby)), 
                                  "cohort_id" = cohort_nearby$iid_combined, 
                                  "date" = cohort_nearby$int_date,
                                  "intarm" = cohort_nearby$intarm,
                                  "t_l" = 0,
                                  "t_u" = 0,
                                  "concur_int_RO" = 0,
                                  "concur_int_RV" = 0,
                                  "concur_int_TO" = 0,
                                  "concur_int_TV" = 0)
        
        add_coho_df <- add_coho_df %>% 
          mutate(t_l = ifelse(intarm %in% c("TV", "RV"), date + t2_l, date + t1_l),
                 t_u = ifelse(intarm %in% c("TV", "RV"), date + t2_u, date + t1_u))
        
        names(add_coho_df)[7:10] <- sapply(names(add_coho_df)[7:10],
                                          function(x) paste0(x, "_", r[k]))
        
        # check overlapping in time
        if (nrow(add_coho_df) == 1){
          add_coho_df[1,c(7:10)] = 0
        }else{
          for (j in 1:nrow(add_coho_df)){
            this_coho_int_day <- add_coho_df$date[j]
            other_coho <- add_coho_df[-j,]
            overlap_coho <- other_coho[which(other_coho$t_l <= this_coho_int_day &
                                               other_coho$t_u >= this_coho_int_day), ]
            if (nrow(overlap_coho) == 0){
              add_coho_df[j,c(7:10)] = 0
            }else{
              # count individual
              overlap_RO <- vector()
              overlap_RV <- vector()
              overlap_TO <- vector()
              overlap_TV <- vector()
              for (m in 1:nrow(overlap_coho)){
                over_int_recip <- int_recip %>% filter(iid_combined == overlap_coho$cohort_id[m])
                if (nrow(over_int_recip) == 0){
                  overlap_RO[m] <- 0
                  overlap_RV[m] <- 0
                  overlap_TO[m] <- 0
                  overlap_TV[m] <- 0
                }else{
                  over_int_recip_ll <- over_int_recip %>% dplyr::select(c("longitude", "latitude"))
                  dist_individual <-  distGeo(this_ll, over_int_recip_ll)
                  near_int_recip <- over_int_recip[which(dist_individual <= r[k]),]
                  overlap_RO[m] <- sum(near_int_recip$intarm == "RO")
                  overlap_RV[m] <- sum(near_int_recip$intarm == "RV")
                  overlap_TO[m] <- sum(near_int_recip$intarm == "TO")
                  overlap_TV[m] <- sum(near_int_recip$intarm == "TV")
                }
              }
              add_coho_df[j,7] = sum(overlap_RO)
              add_coho_df[j,8] = sum(overlap_RV)
              add_coho_df[j,9] = sum(overlap_TO)
              add_coho_df[j,10] = sum(overlap_TV)
            }
            # print(overlap_coho$cohort_id)
          }
        }
        add_coho_df <- add_coho_df[, c(1,2,7:10)]
        df_base_this <- left_join(df_base_this, add_coho_df, by = c("indiv_id", "cohort_id"))
      }
    }
    res_concur_int <- bind_rows(res_concur_int, df_base_this)
  }
  return(res_concur_int)
}



##############################################
##############################################
# Documentation: prop_pre_int
# Usage: prop_pre_int(cohorts, r)
# Description: The proportion of individuals who previously participated in an intervention 
#              within 0.5, 1, 2, and 3km of the individual’s residence
# Args/Options: cohorts data and a sequence of radius
# data:  the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable

prop_pre_int <- function(cohorts, r = c(500, 1000, 2000, 3000)){
  # a list of unique id for each people (n = 19127)
  person_list <- unique(cohorts$indiv_id)
  person_df <- unique(cohorts[,c("indiv_id", "longitude", "latitude")])
  res <- data.frame("indiv_id" = person_list,
                    "prop_pre_int_500m_RO" = 0,
                    "prop_pre_int_1km_RO" = 0,
                    "prop_pre_int_2km_RO" = 0,
                    "prop_pre_int_3km_RO" = 0,
                    "prop_pre_int_500m_RV" = 0,
                    "prop_pre_int_1km_RV" = 0,
                    "prop_pre_int_2km_RV" = 0,
                    "prop_pre_int_3km_RV" = 0,
                    "prop_pre_int_500m_TO" = 0,
                    "prop_pre_int_1km_TO" = 0,
                    "prop_pre_int_2km_TO" = 0,
                    "prop_pre_int_3km_TO" = 0,
                    "prop_pre_int_500m_TV" = 0,
                    "prop_pre_int_1km_TV" = 0,
                    "prop_pre_int_2km_TV" = 0,
                    "prop_pre_int_3km_TV" = 0)
  
  # subset cohorts data to people with int match (n = 7468)
  subdf <- cohorts %>% filter(int_recip == 1)
  
  for (i in 1:length(person_list)){
    print(paste0("---------", round(i/length(person_list), 5)*100,"%","---------"))
    
    this_ll <- cohorts %>% filter(indiv_id == person_list[i]) %>% 
      dplyr::select(c("longitude", "latitude")) %>% 
      unique()
    
    other_df <- subdf %>% filter(indiv_id != person_list[i])
    other_ll <- other_df %>% dplyr::select(c("longitude", "latitude")) 
    
    all_other_df <- person_df %>% filter(indiv_id != person_list[i])
    all_other_ll <- all_other_df %>% dplyr::select(c("longitude", "latitude")) 
    
    for (j in 1:length(r)){
      nearby <- which(distGeo(this_ll, other_ll) <= r[j])
      all_nearby <- which(distGeo(this_ll, all_other_ll) <= r[j])
      
      if (length(nearby) <= 0){
        # should they be discerned ?
        next
      }else{
        nearby_df <- other_df[nearby,] %>% dplyr::select("indiv_id", "intarm") %>% unique()
        all_nearby_df <- all_other_df[all_nearby,] %>% dplyr::select("indiv_id") %>% unique()
        res[i,j+1] = sum(nearby_df$intarm=="RO")/nrow(all_nearby_df)
        res[i,j+5] = sum(nearby_df$intarm=="RV")/nrow(all_nearby_df)
        res[i,j+9] = sum(nearby_df$intarm=="TO")/nrow(all_nearby_df)
        res[i,j+13] = sum(nearby_df$intarm=="TV")/nrow(all_nearby_df)
      }
    }
  }
  return(res)
}


##############################################
##############################################
# Documentation: count_pre_index
# Usage: count_pre_index(cohorts, r)
# Description: The total number of previous index cases within 
#              0.5, 1, 2, and 3km of the individual’s residence
# Args/Options: cohorts data and a sequence of radius
# data:  the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable

count_pre_index <- function(df_base, cohorts, r=c(500,1000,2000,3000)){
  sub_df_base <- df_base %>% 
    arrange(indiv_id, cohort_id, date) %>% 
    dplyr::select(c("indiv_id", "cohort_id", "date"))
  
  person_list <- unique(sub_df_base$indiv_id)
  res <- sub_df_base %>% 
    mutate(pre_index_500m = 0,
           pre_index_1km = 0,
           pre_index_2km = 0,
           pre_index_3km = 0)
  
  # subset cohorts data to people with index match 
  subdf <- cohorts %>% filter(indexcase == 1)
  
  all_res <- map_dfr(person_list, function(x) 
    hp_pre_index(res=res,subdf=subdf, this_person=x))
  
  return(all_res)
}


hp_pre_index <- function(res, subdf, this_person,r=c(500,1000,2000,3000)){
  this_res <- res %>% filter(indiv_id == this_person)
  
  this_ll <- df_base %>% filter(indiv_id == this_person) %>% 
    dplyr::select(c("longitude", "latitude")) %>% unique()
  
  other_df <- subdf %>% filter(indiv_id != this_person)
  other_ll <- other_df %>% dplyr::select(c("longitude", "latitude")) 
  other_df$dist <- distGeo(this_ll, other_ll)
  
  for (i in 1:nrow(this_res)){
    this_res$pre_index_500m[i] = nrow(other_df %>% filter(dist <= r[1] & date <= this_res$date[i]))
    this_res$pre_index_1km[i] = nrow(other_df %>% filter(dist <= r[2] & date <= this_res$date[i]))
    this_res$pre_index_2km[i] = nrow(other_df %>% filter(dist <= r[3] & date <= this_res$date[i]))
    this_res$pre_index_3km[i] = nrow(other_df %>% filter(dist <= r[4] & date <= this_res$date[i]))
  }
  return(this_res)
}


##############################################
##############################################
# Documentation: count_pop
# Usage: count_pop(cohorts, r)
# Description: Population size within 0.5, 1, 2, and 3km of the individual’s residence
# Args/Options: cohorts data and a sequence of radius
# data:  the cohorts data obtained from matching 
# 
# Returns: a data frame with indiv_id and the desired variable

count_pop <- function(cohorts, r=c(500,1000,2000,3000)){
  # a list of unique id for each people (n = 19127)
  person_list <- unique(cohorts$indiv_id)
  person_df <- unique(cohorts[,c("indiv_id", "longitude", "latitude")])
  res <- data.frame("indiv_id" = person_list,
                    "pop_500m" = 0,
                    "pop_1km" = 0,
                    "pop_2km" = 0,
                    "pop_3km" = 0)
  
  
  for (i in 1:length(person_list)){
    print(paste0("---------", round(i/length(person_list), 5)*100,"%","---------"))
    
    this_ll <- cohorts %>% filter(indiv_id == person_list[i]) %>% 
      dplyr::select(c("longitude", "latitude")) %>% 
      unique()
    
    all_other_df <- person_df %>% filter(indiv_id != person_list[i])
    all_other_ll <- all_other_df %>% dplyr::select(c("longitude", "latitude")) 
    
    for (j in 1:length(r)){
      all_nearby <- which(distGeo(this_ll, all_other_ll) <= r[j])
      
      if (length(all_nearby) <= 0){
        # should they be discerned ?
        next
      }else{
        all_nearby_df <- all_other_df[all_nearby,] %>% dplyr::select("indiv_id") %>% unique()
        res[i,j+1] = nrow(all_nearby_df)
      }
    }
  }
  return(res)
}


##############################################
##############################################
# Documentation: coho_classify
# Usage: count_pop(r = 1000, index_int, list_people)
# Description: classify every unique person in the cohort data into different cohorts
# Args/Options: radius of the cohort, 
#               a list of index cases triggering interventions, 
#               a list of unique people in the cohort data
# 
# Returns: a data frame 
#          which is called "df_base", the initial data set for analysis
#          more covariates will be added into this data set

coho_classify <- function(r = 1000, index_int, list_people){
  
  coho_df <- data.frame()
  for (i in 1:nrow(list_people)){
    
    print(paste0("---------", round(i/nrow(list_people), 5)*100,"%","---------"))
    
    this_person <- list_people[i,]
    this_ll <- list_people[i,c(2,3)]
    
    cohort_ll <- index_int %>% dplyr::select(c("longitude", "latitude"))
    dist <-  distGeo(this_ll, cohort_ll)
    cohort_nearby <- index_int[which(dist <= r),]
    
    if (nrow(cohort_nearby) == 0){
      next
    }else{
      
      dist_nearby <- dist[which(dist <= r)]
      
      iid_candidates <- int_recip %>% 
                        filter(indiv_id == this_person$indiv_id) %>% 
                        dplyr::select("iid_combined") %>% pull()
      
      add_coho_df <- data.frame("indiv_id" = rep(this_person$indiv_id, 
                                                 nrow(cohort_nearby)), 
                                "cohort_id" = cohort_nearby$iid_combined, 
                                "date" = cohort_nearby$int_date,
                                "intarm" = cohort_nearby$intarm,
                                "longitude" = this_ll$longitude,
                                "latitude" = this_ll$latitude,
                                "dist_to_index" = dist_nearby,
                                "target_area" = ifelse(dist_nearby <= 500, 1,0),
                                "int_recip" = as.numeric(cohort_nearby$iid_combined %in% iid_candidates))

      coho_df <- bind_rows(coho_df, add_coho_df)
    }
  }
  return(coho_df)
}







