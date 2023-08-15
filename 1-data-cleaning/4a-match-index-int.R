# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create cohort data structure

# Match index cases that triggered interventions
# to individual imputed dataset 
# ............................................

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# ............................................
# load data
# ............................................

# baseline data to be used to create a list of IDs for
# all individuals in the study area
indiv_clean = readRDS(namibia_indiv_imputed_path) 

# cleaned index data with imputed GPS for index cases
# missing this info
index = readRDS(namibia_clean_index_path)

# cleaned intervention data 
intervention =readRDS(namibia_clean_intervention_path)

# ............................................
# drop if missing lat/long 
# ............................................

# how many iids have missing lat/long
# number of index case ids that triggered intervention missing lat/long
index %>% filter(is.na(longitude_index_imp) & first_case_yn==1) %>% nrow()

# number of other index case ids missing lat/long
index %>% filter(is.na(longitude_index_imp) & first_case_yn!=1) %>% nrow()

# number of intervention recipient ids missing lat/long
intervention %>% filter(is.na(longitude_int)) %>% 
  dplyr::select(iid_combined) %>% distinct() %>%  nrow()

# drop if missing lat/long
intervention = intervention %>% filter(!is.na(longitude_int))
index = index %>% filter(!is.na(longitude_index_imp))

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
  rename(ea = ea_index)

#....................................................
# check that all the index cases have matching intervention people
#....................................................

index_iids = index %>% filter(first_case_yn == 1) %>%
  dplyr::select(iid_combined) %>% distinct() 

intervention_iids_N = intervention %>%
  group_by(iid_combined) %>% 
  summarise(N = n())

index_int_N = full_join(index_iids, intervention_iids_N, by = "iid_combined")

index_iids_list = index_iids$iid_combined
intervention_iids_list = intervention_iids_N$iid_combined

# create unique identifier for index cases because
# iid is used for more than one index case if they 
# occurred in the same place in the same period 
index$indexid = seq(1, nrow(index))

# subset to index cases that triggered intervention
index_int = index %>% filter(first_case_yn==1)

# impute missing visit_date by the earliest int_day
temp_i <- which(is.na(index_int$visit_date_index))

temp_iid_combined <- index_int[temp_i, ]$iid_combined

temp <- intervention %>% filter(iid_combined %in% temp_iid_combined)

for (i in 1:length(temp_iid_combined)){
  index_int[temp_i, ]$visit_date_index[i] = 
    min(temp[which(temp$iid_combined == temp_iid_combined[i]), ]$start_date_int)
}


index_int <- index_int %>% mutate(clinic_date_index = 
                                    if_else(clinic_date_index > visit_date_index,
                                           visit_date_index, clinic_date_index))


# temp <- intervention %>% filter(iid_combined == "TV34-5270")

# ............................................
# create dataset with each follow up day and each person 
# ............................................

# list of dates and ids for entire study
# study_dates = seq.Date(from = min(index$clinic_date_index, na.rm=T),
#                        to = max(index$visit_date_index, na.rm=T)+60,
#                        by = "day")

study_dates = seq.Date(from = min(index$clinic_date_index, na.rm=T),
                       to = max(index$visit_date_index, na.rm=T)+90,
                       by = "day")

person_ids = unique(indiv_clean$indiv_id)

# create dataset with rows for each person and date
id_date = expand_grid(person_ids, study_dates) %>%
  rename(indiv_id = person_ids,
         date = study_dates)

# merge lat/long onto dataset
person_time = full_join(id_date, indiv_clean, by = "indiv_id")
assert_that(nrow(id_date) == nrow(person_time))

# ............................................
# match index cases to imputed person-time dataset 
# based on location and clinic or visit date
# ............................................

index_unmatched_list = list()
index_matched_list = list()
master_ppl1 = indiv_clean %>% 
  mutate(indexid = as.integer(NA),
         evermatch = 0,
         age = as.integer(NA),
         sex = as.factor(NA)) %>% ungroup()
master_pdays1 = person_time %>% 
  mutate(indexid = as.integer(NA),
         evermatch = 0,
         age = as.integer(NA),
         sex = as.factor(NA),
         evermatchday = 0) %>% ungroup()



# matched is an indicator for whether this individual id was
# ever previously matched to a person

# Matching ------------------------------------------
tic()
for(i in 1:length(index_int$indexid)){
  
  print(paste(i,"----------------------------"))
  
  # index id
  this_id_df = index_int %>% filter(indexid == unique(index_int$indexid)[i])
  
  # identify dates for this index case
  clinic_date_index = this_id_df$clinic_date
  visit_date_index = this_id_df$visit_date
  this_date <- clinic_date_index
  
  # grab lat and long for index case and each imputed individual
  index_ll = this_id_df %>%
    dplyr::select(longitude_index_imp,
                  latitude_index_imp) %>%
    rename(long = longitude_index_imp,
           lat = latitude_index_imp)
  
  # set master list for this iteration only 
  master_ppl1_i = master_ppl1
  master_pdays1_i = master_pdays1
  
  # 02-17-2021 calculate distance
  indiv_ll <- data.frame(long = master_ppl1_i$longitude,
                         lat = master_ppl1_i$latitude)
  
  dist <- distGeo(index_ll[, c("long", "lat")], indiv_ll[, c("long", "lat")])
  
  reason_msg = paste0("")
  
  repeat{
    print(paste0("Available people inside repeat: ", nrow(master_ppl1_i)))
    
    # ............................................
    # check that minimum distance is <= 100m
    # if not, do not match
    # ............................................
    
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
      
      index_unmatched_list[[i]] = data.frame(indexid = this_id_df$indexid, 
                                             reason = paste0(temp_msg,
                                                             "No match within 100m."))
      
      person_match_index=NULL
      break
      
    }else{
      ## 1.2 Else,  match closest individual ------------------------------------------
      matched_indiv_id = master_ppl1_i[which.min(dist),] %>%  
        dplyr::select(indiv_id) %>% pull() 
      
      # 2. Check if there was a previous match ------------------------------------------
      if(master_ppl1_i %>% filter(indiv_id==matched_indiv_id) %>% pull(evermatch) == 1){
        
        ## 2.1 If yes, check age, sex match and date not taken ------------------------------------------
        pre_age = master_ppl1_i %>% 
          filter(indiv_id == matched_indiv_id) %>%
          pull(age)
        
        pre_sex = master_ppl1_i %>% 
          filter(indiv_id == matched_indiv_id) %>%
          pull(sex)
        
        this_age <- this_id_df$age_index
        
        prev_match_list <- master_pdays1 %>%  
                           filter(indiv_id==matched_indiv_id &
                                  evermatchday == 1)
        
        last_date <- prev_match_list$date[nrow(prev_match_list)]
        
        age_diff <- ifelse(last_date <= this_date, this_age - pre_age, pre_age - this_age)
        
        # 06/21/2021 age tolerance
        agematch = is.na(pre_age) | (0 <= age_diff & age_diff <= 1)
        sexmatch = is.na(pre_sex) | pre_sex == this_id_df$sex_index
        
        prev_day_match = master_pdays1 %>% 
          filter(indiv_id == matched_indiv_id & 
                   (date == clinic_date_index)) %>% 
          pull(evermatchday) == 1
        
        print(paste0("pre_age: ", pre_age))
        print(paste0("this_age: ", this_id_df$age_index))
        
        print(paste0("pre_sex: ", pre_sex))
        print(paste0("this_sex: ", this_id_df$sex_index))
        
        # if age and sex match and the date is available, check wash-out 
        if (agematch & sexmatch & !is.na(agematch) & !is.na(sexmatch) & is.na(prev_day_match)){
          
          ### 2.1.1 If yes, check wash-out ------------------------------------------
          print("Checking wash-out period")
          
          prev_match_list <-  master_pdays1 %>% 
            filter(indiv_id==matched_indiv_id,
                   evermatchday ==1) 
          
          washout_check <- check_washout(prev_match_list = prev_match_list, 
                                         this_date = clinic_date_index,
                                         mindays = 35)
          
          if (!washout_check){
            temp_msg = "Index appeared in the disease wash-out period; proceed to next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "wash-out fail.")
            
            index <- master_ppl1_i$indiv_id != matched_indiv_id
            master_ppl1_i = master_ppl1_i[index,]
            indiv_ll <- indiv_ll[index,]
            dist <- dist[index]
            
            person_match_index = NULL
          }else{
            print("Successful match")
            person_match_index = master_pdays1 %>% 
              filter(indiv_id==matched_indiv_id & (date == clinic_date_index |
                                                     date == visit_date_index))
            
            break
          }
          
          # ............................................
          
          ## 2.1.2 If no, next closest person ------------------------------------------
          # if age and sex do not match or date has been taken, continue search 
        }else{
          
          if (!is.na(prev_day_match) & prev_day_match){
            temp_msg = "Date has been taken. Try next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "Date taken.")
          }else{
            temp_msg = "Age and/or sex do not match with previous assigned person. Try next nearest person."
            print(temp_msg)
            reason_msg = paste0(reason_msg, "A/S not match.", paste())
          }
          
          index <- master_ppl1_i$indiv_id != matched_indiv_id
          master_ppl1_i = master_ppl1_i[index,]
          indiv_ll <- indiv_ll[index,]
          dist <- dist[index]
          
          person_match_index = NULL
          
        }
        # ............................................
        
        ## 2.2 If no, success -----------------------------------------
      }else{
        
        print("Successful match")
        person_match_index = master_pdays1 %>% 
          filter(indiv_id==matched_indiv_id & (date == clinic_date_index |
                                                 date == visit_date_index))
        break
        
        # end of check for previous match 
      }
      # end of distance check  
    }
    # end of repeat    
  }
  # ............................................
  
  if(!is.null(person_match_index)){
    person_match_index = person_match_index %>%
      mutate(event = case_when(
        date == index_int$clinic_date_index[i] ~ "Clinic visit",
        date == index_int$visit_date_index[i] ~ "Intervention",
        TRUE ~ ""
      )) %>%
      mutate(indexid = index_int$indexid[i],
             iid_combined = index_int$iid_combined[i],
             first_case_yn = index_int$first_case_yn[i],
             age = index_int$age_index[i],
             sex = index_int$sex_index[i],
             clinic_date_index = index_int$clinic_date_index[i],
             visit_date_index = index_int$visit_date_index[i],
             dist_to_match = dist[which.min(dist)],
             evermatch = 1,
             evermatchday = 1
      )
    
    index_matched_list[[i]] = person_match_index
    
    # ............................................
    # assign age and sex to this indiv_id so that 
    # it is factored into subsequent matches 
    # ............................................
    
    master_pdays1 = master_pdays1 %>% mutate(
      indexid = if_else(indiv_id == matched_indiv_id, unique(person_match_index$indexid), indexid),
      evermatch = if_else(indiv_id == matched_indiv_id, 1, evermatch),
      evermatchday = if_else(indiv_id == matched_indiv_id & (date == clinic_date_index | 
                                                              date == visit_date_index), 
                            1, evermatchday)
    )
    
    master_ppl1 = master_ppl1 %>% mutate(
      indexid = if_else(indiv_id == matched_indiv_id, unique(person_match_index$indexid), indexid),
      age = if_else(indiv_id == matched_indiv_id, unique(person_match_index$age), age),
      sex = if_else(indiv_id == matched_indiv_id, unique(person_match_index$sex), sex),
      evermatch = if_else(indiv_id == matched_indiv_id, 1, evermatch)
    )
  }
}
toc()


indiv_matched_df = bind_rows(index_matched_list)
indiv_unmatched_df = bind_rows(index_unmatched_list)

# check that master list and newly matched objects have 
# same number of matches
assert_that(nrow(master_ppl1[master_ppl1$evermatch==1,]) ==
              length(unique(indiv_matched_df$indexid)))

# check that number of rows equals number of index
# cases with unique dates for clinic and visit
assert_that(
  length(unique(indiv_matched_df$indexid))+ 
    length(unique(indiv_unmatched_df$indexid)) ==
    index_int %>% nrow()
)

# ensure that the same person was not selected for more than
# one index case
assert_that(any(duplicated(unique(indiv_matched_df$indiv_id)))==F)

# ............................................
# combine index and individual data into a single data frame
# ............................................

# identify individuals not identified as index cases
indiv_matched_df = indiv_matched_df %>%
  dplyr::select(new_hhid, personid, indiv_id, indexid, iid_combined,
                date, clinic_date_index, visit_date_index,
                first_case_yn, age, sex, evermatch, evermatchday) %>%
  mutate(clinic_day = ifelse(date == clinic_date_index, 1, 0),
         int_day = ifelse(date == visit_date_index, 1, 0)) %>%
  mutate(indexcase = 1,
         int_trig = 1,) %>%
  dplyr::select(-c(clinic_date_index, visit_date_index))

all_index = full_join(person_time, indiv_matched_df, by = c("new_hhid", "personid", "indiv_id", "date"))

# check that merge did not add rows
assert_that(nrow(person_time) == nrow(all_index))

# check for duplicates in individual ids
assert_that(any(duplicated(unique(all_index$indiv_id)))==F)

# for a given indiv_id, if match occurred, 
# impute index id, age, sex for all other dates 
all_index_matched = all_index %>% 
  group_by(indiv_id) %>% 
  fill(indexid, .direction = "downup") %>% 
  fill(age, .direction = "downup") %>% 
  fill(sex, .direction = "downup") %>% 
  fill(evermatch, .direction = "downup") %>% 
  ungroup()


saveRDS(person_time, file = namibia_pt_path)
saveRDS(all_index_matched, file = namibia_index_matched_path)

