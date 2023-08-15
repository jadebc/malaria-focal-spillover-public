# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create cohort data structure

# Match index cases that did not trigger interventions
# to individual imputed dataset 
# ............................................

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# ............................................
# Load data
# ............................................
# cleaned index data with imputed GPS for index cases
# missing this info
index = readRDS(namibia_clean_index_path)

# data generated in 4b with intervention individuals
# and indexes that triggered intervention matched 
# to individual person-day dataset
all_index_int_imputed = readRDS(namibia_int_match_path)

person_time = readRDS(namibia_pt_path)

# ............................................
# impute individuals who received interventions
# who were not in the original georeconaissance
# survey
# ............................................

intervention =readRDS(namibia_clean_intervention_path)
intervention_unmatched_df = readRDS(namibia_temp_int_unmatch_path)

temp_int <- intervention %>% filter(id %in% intervention_unmatched_df$interventionid)

generate_index <- function(prefix, id){
  n = length(id)
  res <- rep(NA, n)
  for(i in 1:n){
    res[i] = paste0(prefix, "_", id[i])
  }
  return(res)
}

temp_int <- temp_int %>%
  mutate(indiv_id = generate_index("000", temp_int$id),
         longitude = longitude_int,
         latitude = latitude_int,
         interventionid = id) %>%
  dplyr::select(all_of(c("indiv_id", "longitude", "latitude",
                         "ea", "intarm", "iid_combined",
                         "age", "sex", "start_date_int", "sha_id", "interventionid")))

seed <- all_index_int_imputed[1:319, ]
seed <- seed %>%
  dplyr::select(c("indiv_id", setdiff(names(seed), names(temp_int))))

seed[, -c(2,7)] <- as.numeric(NA)

df_rep <- function(df, n){
  df0 = df
  while( n > 0){
    df = bind_rows(df, df0)
    n = n -1
  }
  return(df)
}

leaf <- df_rep(seed, nrow(temp_int) - 1)
leaf$indiv_id <- rep(temp_int$indiv_id, each = 319)

temp_add <- left_join(leaf, temp_int, by = "indiv_id")

temp_add[which(temp_add$date == temp_add$start_date_int), c(8,10,14)] <- 1
temp_add[which(temp_add$date == temp_add$start_date_int), c(11)] <- 0


temp_add <- temp_add %>% dplyr::select(all_of(names(all_index_int_imputed)))

#
all_index_int_imputed <- bind_rows(all_index_int_imputed, temp_add)

# sapply(all_index_int_imputed, class)


# ............................................
# clean data 
# ............................................
# clean index data
index = index %>% 
  dplyr::select(
    iid_combined,
    longitude_index_imp, latitude_index_imp,
    ea_index, intarm, 
    clinic_date_index, visit_date_index, 
    first_case_yn,
    age_index, sex_index
  ) %>%
  dplyr::rename(ea = ea_index)

# create unique identifier for index cases because
# iid is used for more than one index case if they 
# occurred in the same place in the same period 
index$indexid = seq(1, nrow(index))

# manually drop index that has ea is not present 
# in the intervention ea list
index = index[-which(index$ea=="TO41"),]

# drop if missing lat/long
index = index %>% filter(!is.na(longitude_index_imp))

# ............................................
# match non intervention triggering index cases
# to the person-time dataset based on date and
# location
# ............................................
index_nonint_matched_list = list()
index_nonint_unmatched_list = list()

# create master list data frame with ids and age and sex for 
# each individual 
master_ppl3 = all_index_int_imputed %>% 
  dplyr::select(indiv_id, longitude, latitude, ea,
                age, sex, evermatch) %>% 
  distinct()

master_pdays3 = all_index_int_imputed 

# subset to index cases that did not trigger intervention
index_nonint =  index %>% filter(first_case_yn!=1)

# check that all index cases have clinic date
assert_that(all(!is.na(index_nonint$clinic_date)))

# i - 685 - something wrong with how person_match is generated
# the clinic date match check failed
# population available = 320

# ............................................
# match non triggering index cases to unmatched
# person-time observations based on who is
# closest and the date that matches the clinic date
# ............................................

# Matching ------------------------------------------
tic()
for(i in 1:length(index_nonint$indexid)){
  
  pct_complete = round(i / length(index_nonint$indexid) *100, 1)
  print(paste0("Index case #: ", i, "; Completed: ", pct_complete, "% --------------------"))
  
  # index id
  this_id_df = index_nonint %>% filter(indexid == unique(index_nonint$indexid)[i])
  
  # identify index EA
  ea = this_id_df$ea
  
  # identify iid
  iid_comb = this_id_df$iid_combined
  
  # identify arm
  this_arm = this_id_df$intarm
  
  # identify dates for this index case
  clinic_date_index = this_id_df$clinic_date
  visit_date_index = this_id_df$visit_date
  this_date = clinic_date_index
  
  # grab lat and long for intervention person and each imputed individual
  nonint_index_ll = data.frame(long = this_id_df$longitude_index_imp,
                               lat = this_id_df$latitude_index_imp)
  
  #people_ea = master_ppl3[master_ppl3$ea == ea,]
  people_ea = master_ppl3
  print(paste("nrow(people_ea)", nrow(people_ea)))
  
  reason_msg = paste0("")
  repeat{
    
    print(paste("Number of people available to match:", nrow(people_ea)))
    
    if(nrow(people_ea)==0){
      
      print("No one available in this ea. No match.")
      index_nonint_unmatched_list[[i]] = data.frame(indexid = this_id_df$indexid, 
                                                    reason = paste0(reason_msg,
                                                                    "No age/sex matches in this EA.")
      )
      
      person_match=NULL
      break
    }
    
    # create data frame with the lat/long for unmatched individuals
    # in the same EA on the same date as the intervention
    people_ea_ll = data.frame(long = people_ea$longitude,
                              lat = people_ea$latitude)
    
    # calculate meters between index case and each imputed individual 
    dist <- distGeo(nonint_index_ll[, c("long", "lat")], people_ea_ll[, c("long", "lat")])
    
    ## 1. Check distance ------------------------------------------
    if(dist[which.min(dist)] > 100){
      ## 1.1 If no one within 100m, failed ------------------------------------------
      
      print("No match available within 100m in this EA on this date")
      
      if (str_count(reason_msg, pattern = "A/S not match") > 0){
        temp_msg = 
          str_replace_all(reason_msg, 
                          "A/S not match.+A/S not match", 
                          paste0("A/S not match *", str_count(reason_msg, pattern = "A/S not match")))
        
      }else{
        temp_msg = reason_msg
      }
      
      index_nonint_unmatched_list[[i]] = data.frame(indexid = this_id_df$indexid, 
                                                    reason = paste0(temp_msg,
                                                                    "No match within 100m.")
      )
      
      person_match=NULL
      break
      
      ## 1.2 Else,  match closest individual ------------------------------------------
    }else{
      # match closest individual on that clinic date in that EA
      matched_indiv_id = people_ea[which.min(dist),] %>%  
        dplyr::select(indiv_id) %>% pull() 
      
      # 2. Check previous match ------------------------------------------
      # check if there was a previous match for this index id
      pre_match = people_ea %>% filter(indiv_id==matched_indiv_id) %>% pull(evermatch)
      if(!is.na(pre_match) & pre_match == 1){
        ## 2.1 If yes, check age, sex match and date not taken ------------------------------------------
        # check if sex age match, and date has not been taken
        pre_age = people_ea %>% 
          filter(indiv_id == matched_indiv_id) %>%
          pull(age)
        
        if (length(pre_age) > 1 & any(!is.na(pre_age))){
          pre_age <- max(pre_age[!is.na(pre_age)])
        }
        
        pre_sex = people_ea %>% 
          filter(indiv_id == matched_indiv_id) %>%
          pull(sex)
        
        if (length(pre_sex) > 1 & any(!is.na(pre_sex))){
          pre_sex <- unique(pre_sex[!is.na(pre_sex)])
          assert_that(length(pre_sex) == 1,
                      msg = "Different genders for the same person")
        }
        
        this_age <- this_id_df$age_index
        
        prev_match_list <- master_pdays3 %>%  
          filter(indiv_id==matched_indiv_id &
                   evermatchday == 1)
        
        last_date <- prev_match_list$date[nrow(prev_match_list)]
        age_diff <- ifelse(last_date <= this_date, this_age - pre_age, pre_age - this_age)
        
        # 06/21/2021 age tolerance
        agematch = is.na(pre_age) | (0 <= age_diff & age_diff <= 1)
        sexmatch = is.na(pre_sex) | pre_sex == this_id_df$sex_index
        
        print(paste0("pre_age: ", people_ea %>% 
                       filter(indiv_id == matched_indiv_id) %>%
                       pull(age)))
        print(paste0("this_age: ", this_id_df$age_index))
        
        print(paste0("pre_sex: ", people_ea %>% 
                       filter(indiv_id == matched_indiv_id) %>%
                       pull(sex)))
        print(paste0("this_sex: ", this_id_df$sex_index))
        
        print(this_id_df$indexid)
        
        prev_day_match = identical(master_pdays3 %>% 
                          filter(indiv_id == matched_indiv_id & 
                                   (date == clinic_date_index)) %>% 
                          pull(evermatchday), 1)

        
        if (agematch & sexmatch & !is.na(agematch) & !is.na(sexmatch) & !prev_day_match){
          ### 2.1.1 If yes, check wash-out  ------------------------------------------
          # if age and sex match and the date is available, check wash-out 
          print("Checking wash-out period")
          
          prev_match_list <-  master_pdays3 %>% 
            filter(indiv_id==matched_indiv_id,
                   evermatchday ==1) 
          
          washout_check <- check_washout(prev_match_list = prev_match_list, 
                                         this_date = clinic_date_index)
          
          if (!washout_check){
            ### 2.1.1,1 If no, next closest person  ------------------------------------------
            temp_msg = "Index appeared in the disease wash-out period; proceed to next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "wash-out fail.")
            people_ea = people_ea[people_ea$indiv_id != matched_indiv_id,]
            person_match = NULL
          }else{
            ### 2.1.1.2 If yes, success  ------------------------------------------
            print("Successful match")
            person_match = master_pdays3 %>% 
              filter(indiv_id==matched_indiv_id & (date == clinic_date_index))
            
            break
          }
          
        }else{
          ### 2.1.2 If no, next closest person  ------------------------------------------
          # if age and sex do not match or date has been taken, continue search 
          if (!is.na(prev_day_match) & prev_day_match){
            temp_msg = "Date has been taken. Try next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "Date taken.")
          }else{
            temp_msg = "Age and/or sex do not match with previous assigned person. Try next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "A/S not match.", paste())
          }
          
          people_ea = people_ea[people_ea$indiv_id != matched_indiv_id,]
          person_match = NULL
          
        }
        
        # ............................................
        # if no previous match, check whether person-days are 
        # available for this clinic and visit day
      }else{
        ## 2.2 If no, success ------------------------------------------
        
        
        print("Successful match")
        person_match = master_pdays3 %>% 
          filter(indiv_id==matched_indiv_id & (date == clinic_date_index))
        
        break
        
        View(master_pdays3 %>% 
               filter(indiv_id==matched_indiv_id))
        # clinic_day_match = master_pdays3 %>% 
        #   filter(indiv_id==matched_indiv_id & date == clinic_date_index) %>% nrow() ==1
        # 
        # 
        # # if clinic and visit days are unassigned for this person, proceed with match
        # if(clinic_day_match){
        #   
        #   print("Successful match")
        #   person_match = master_pdays3 %>% 
        #     filter(indiv_id==matched_indiv_id & (date == clinic_date_index))
        #   
        #   break
        #   
        #   # if clinic or visit days not available, find next closest person 
        # }else{
        #   temp_msg ="These dates not available for this person; find next nearest person"
        #   print(temp_msg)
        #   reason_msg = paste0(reason_msg, temp_msg)
        #   people_ea = people_ea[people_ea$indiv_id != matched_indiv_id,]
        #   person_match = NULL
        #   
        # }
        
        # end of check for previous match 
      }
      
      # end of distance check  
    }
    # end of repeat    
  }
  
  
  if(!is.null(person_match)){
    # TEMPORARY for debugging
    # print("this_id_df----")
    # print(this_id_df)
    # print("master_person----")
    # print( people_ea %>% filter(indiv_id==matched_indiv_id))
    # print("master_pdays3----")
    # print( master_pdays3 %>% 
    #          filter(indiv_id==matched_indiv_id & (date == clinic_date_index)))
    # 
    # print("person_match----")
    # print(person_match)
    
    
    person_match = person_match %>% 
      mutate(
        indexid = this_id_df$indexid,
        iid_combined = iid_comb,
        first_case_yn = this_id_df$first_case_yn,
        age = this_id_df$age_index,
        sex = this_id_df$sex_index,
        dist_to_match = dist[which.min(dist)],
        evermatch = 1,
        evermatchday = 1,
        clinic_day = case_when(
          date == clinic_date_index ~ 1,
          date != clinic_date_index ~ 0,
          is.na(clinic_date_index) ~ NA_real_
        ),
        int_day = case_when(
          date == visit_date_index ~ 1,
          date != visit_date_index ~ 0,
          is.na(visit_date_index) ~ NA_real_
        ),
        indexcase = 1,
        int_trig = 0
      )
    
    # check that only 1 person-day was selected
    assert_that(length(unique(person_match$indiv_id))==1,
                msg = "More than 1 individual id matched")
    
    index_nonint_matched_list[[i]] = person_match
    
    # ............................................
    # assign age and sex to this indiv_id so that 
    # it is factored into subsequent matches 
    # ............................................
    if (!is.null(person_match)){
      master_ppl3 = master_ppl3 %>% mutate(
        age = dplyr::if_else(indiv_id == unique(person_match$indiv_id), unique(person_match$age), age),
        sex = dplyr::if_else(indiv_id == unique(person_match$indiv_id), unique(person_match$sex), sex),
        evermatch = dplyr::if_else(indiv_id == unique(person_match$indiv_id), 1, evermatch))
      
      
      master_pdays3 = master_pdays3 %>% mutate(
        indexid = ifelse(indiv_id == matched_indiv_id, unique(person_match$indexid), indexid),
        evermatch = dplyr::if_else(indiv_id == matched_indiv_id, 1, evermatch),
        evermatchday = dplyr::if_else(indiv_id == matched_indiv_id & (date == clinic_date_index), 
                              1, evermatchday)
      )
    }
    
    # end of check if match occurred 
    rm(person_match)
    
  }
  
}
toc()


index_nonint_matched_df = bind_rows(index_nonint_matched_list)
index_nonint_unmatched_df = bind_rows(index_nonint_unmatched_list)

saveRDS(index_nonint_matched_df,
        file = namibia_temp_nontrig_match_path)
saveRDS(index_nonint_unmatched_df,
        file = namibia_temp_nontrig_unmatch_path)

index_nonint_matched_df = readRDS(namibia_temp_nontrig_match_path)
index_nonint_unmatched_df = readRDS(namibia_temp_nontrig_unmatch_path)


# check that number of ids matched + unmatched = total
assert_that(
  index_nonint_matched_df %>% dplyr::select(indexid) %>% distinct() %>% nrow() +
    nrow(index_nonint_unmatched_df) == 
    length(index_nonint$indexid),
  msg = "Some ids unaccounted for following matching."
)

# ensure that the same person was not selected for more than
# one index case
assert_that(any(duplicated(unique(index_nonint_matched_df$indiv_id)))==F)

# how many non-triggering index cases also received interventions
index_nonint_matched_df %>% filter(indexcase==1 & int_recip==1) %>% nrow()
index_nonint_matched_df %>% filter(indexcase==1 & int_recip==0) %>% nrow()

# ............................................
# combine non-intervention triggering index cases 
# with matched intervention and index data 
# ............................................
# create person-day dataset that only includes
# individuals that weren't matched to intervention recipients
matched_nonint_ids = index_nonint_matched_df %>%
  dplyr::select(indiv_id, date) %>%
  mutate(drop=1)

all_index_int_sub = full_join(all_index_int_imputed, 
                              matched_nonint_ids, 
                              by = c("indiv_id", "date")) %>%
  filter(is.na(drop)) %>%
  dplyr::select(-drop)

# combine intervention matches into data frame with
# index matches and non-matches
all = bind_rows(all_index_int_sub, index_nonint_matched_df) %>%
  arrange(indiv_id, date, ea)

# check for duplicates in individual ids
assert_that(any(duplicated(unique(all$indiv_id)))==F)

# count number of people matched

# check that number of index cases that triggered intervention is correct
assert_that(
  all %>% 
    filter(int_trig==1) %>% 
    dplyr::select(indiv_id, indexcase, indexid) %>%
    distinct() %>%  nrow() ==

  all_index_int_imputed %>% 
    filter(int_trig==1) %>% 
    dplyr::select(indiv_id, indexcase, indexid) %>%
    distinct() %>%  nrow()
)

# ............................................
# impute individuals who were non triggering index
# cases who were not in the original georeconaissance
# survey
# ............................................

temp_index <- index_nonint %>% filter(indexid %in% index_nonint_unmatched_df$indexid)

generate_index <- function(prefix, id){
  n = length(id)
  res <- rep(NA, n)
  for(i in 1:n){
    res[i] = paste0(prefix, "_", id[i])
  }
  return(res)
}

temp_index <- temp_index %>%
  mutate(indiv_id = generate_index("0000", temp_index$indexid),
         longitude = longitude_index_imp,
         latitude = latitude_index_imp,
         indexid = indexid,
         age = age_index,
         sex = sex_index) %>%
  dplyr::select(all_of(c("indiv_id", "longitude", "latitude",
                         "ea", "intarm", "iid_combined",
                         "age", "sex", "indexid", "clinic_date_index")))

seed <- all[1:319, ]
seed <- seed %>%
  dplyr::select(c("indiv_id", setdiff(names(seed), names(temp_index))))
seed[, -c(2,6)] <- as.numeric(NA)

df_rep <- function(df, n){
  df0 = df
  while( n > 0){
    df = bind_rows(df, df0)
    n = n -1
  }
  return(df)
}

leaf <- df_rep(seed, nrow(temp_index) - 1)
leaf$indiv_id <- rep(temp_index$indiv_id, each = 319)

temp_add <- left_join(leaf, temp_index, by = "indiv_id")

temp_add[which(temp_add$date == temp_add$clinic_date_index), c(7,8,10)] <- 1
temp_add[which(temp_add$date == temp_add$clinic_date_index), c(11)] <- 0


temp_add <- temp_add %>% dplyr::select(all_of(names(all)))


all <- bind_rows(all, temp_add)


saveRDS(all, namibia_cohort_path)
