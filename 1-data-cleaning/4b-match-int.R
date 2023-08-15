# ............................................
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Create cohort data structure

# Match intervention individuals 
# to individual imputed dataset 
# ............................................

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

all_index_matched = readRDS(namibia_index_matched_path)
all_index_matched$sha_id <- 0

# cleaned intervention data 
intervention =readRDS(namibia_clean_intervention_path)

# drop if missing lat/long
intervention = intervention %>% filter(!is.na(longitude_int))

person_time = readRDS(namibia_pt_path)

# ............................................
# identify intervention individuals in each ea
# and match them to imputed dataset 
# ............................................
int_matched_list = list()
int_unmatched_list = list()
master_pdays2 = all_index_matched %>% 
  mutate(interventionid = as.integer(NA))
master_ppl2 = all_index_matched %>% 
  group_by(indiv_id, ea, longitude, latitude, age, sex) %>% 
  summarise(evermatch = max(evermatch, na.rm=T)) %>% 
  mutate(evermatch = ifelse(evermatch!=1, 0, evermatch)) %>% 
  ungroup() %>% 
  mutate(interventionid = as.integer(NA),
         sha_id = 0) 


# check that all intervention EAs are included
# in the individual dataset
assert_that(all(unique(intervention$ea) %in% all_index_matched$ea),
            msg = "Some EAs in intervention dataset are not in the individual dataset.")

# loop over each unique intervention id and
# match each intervention individual to the
# individual in the imputed dataset with
# the closest distance on a given date

# Matching ------------------------------------------
tic()
for(i in 1:length(unique(intervention$id))){
    
  ptm <- proc.time()[3]
  print(paste0("i: ", i, "; Complete: ", 
               round(i/length(unique(intervention$id))*100, 2), "%--------------------------------"))
  
  # intervention id
  this_id_df = intervention %>% filter(id == unique(intervention$id)[i])
  
  # identify intervention EA
  this_ea = this_id_df$ea
  
  # identify iid
  iid_comb = this_id_df$iid_combined
  
  # identify sha_id
  this_sha_id = this_id_df$sha_id
  
  # identify arm
  this_arm = ifelse(this_id_df$intarm %in% c("RO", "TO"), 0, 1)
  
  # identify intervention start date for this intervention individual
  int_date = this_id_df$start_date_int
  this_date <- int_date
  
  # data frame with unmatched individuals in the same ea 
  master_ppl2_loop = master_ppl2 %>% filter(ea == this_ea)
  master_pdays2_loop = master_pdays2 %>% filter(ea == this_ea)
  
  # grab lat and long for intervention person and each imputed individual
  intervention_ll = data.frame(long = this_id_df$longitude_int,
                               lat = this_id_df$latitude_int)
  
  # create data frame with the lat/long for unmatched individuals
  # in the same EA on the same date as the intervention
  indiv_ll = data.frame(long = master_ppl2_loop$longitude,
                        lat = master_ppl2_loop$latitude)
  
  # calculate meters between intervention person and each imputed individual
  dist <- distGeo(intervention_ll[, c("long", "lat")],
                  indiv_ll[, c("long", "lat")])
  
  reason_msg = paste0("")
  repeat{
    print(paste("Number of people available to match:", nrow(master_ppl2_loop)))
    
    assert_that(nrow(indiv_ll) > 0,
                msg = "No individuals available to match in this EA on this date.")
    
    
    ## 1. Check distance ------------------------------------------
    # check that minimum distance is <= 100m
    # if not, do not match
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
      
      int_unmatched_list[[i]] = data.frame(interventionid = this_id_df$id, 
                                           reason = paste0(temp_msg, "No match within 100m."))
      
      person_match_int=NULL
      
      break
      
      ## 1.2 Else,  match closest individual ------------------------------------------
    }else{
      # identify individual with closest GPS to intervention recipient
      matched_indiv_id = master_ppl2_loop[which.min(dist),] %>%
        dplyr::select(indiv_id) %>% pull()
      
      # matched_indiv_id = master_ppl2_loop[26,] %>%
      #   dplyr::select(indiv_id) %>% pull()

      
      # 2. Check date availability ------------------------------------------
      if(master_pdays2_loop %>% 
         filter(indiv_id == matched_indiv_id & date==int_date) %>% 
         nrow() ==0 ){
        
        ## 2.1 If no, next closest person ------------------------------------------
        # if there are no available observations for this person-day, 
        # proceed to next closest match
        print("No person-days available for this person and day, proceed to next closest person.")
        reason_msg = paste0(reason_msg, "No available date.")
        
        index <- master_ppl2_loop$indiv_id != matched_indiv_id
        master_ppl2_loop = 
          master_ppl2_loop[index,]
        # 02-17-2021
        indiv_ll <- indiv_ll[index,]
        dist <- dist[index]
        
      }else{
        ## 2.2 If yes, check previous match ------------------------------------------
        # if this person was already matched to intervention recipient, 
        # check that age and sex match
        # if not, go to next closest person
        prev_match = master_ppl2_loop %>% 
          filter(indiv_id == matched_indiv_id) %>% 
          pull(evermatch) == 1
        
        if(prev_match){
          is_index_match <- nrow(master_pdays2 %>% 
                                 filter(indiv_id == matched_indiv_id & date == int_date &
                                         indexcase==1)) != 0
          
          if(is_index_match){
            print("Person-day already matched to index case")
            reason_msg = paste0(reason_msg, "Matched index case")
            
            index <- master_ppl2_loop$indiv_id != matched_indiv_id
            master_ppl2_loop = 
              master_ppl2_loop[index,]
            indiv_ll <- indiv_ll[index,]
            dist <- dist[index]
            person_match_int = NULL
          }else{
            
          
          ### HERE
          is_int_match <- nrow(master_pdays2 %>% 
                                 filter(indiv_id == matched_indiv_id & !is.na(interventionid))) != 0
          
          candid_sha_id = master_ppl2 %>%
            filter(indiv_id == matched_indiv_id & !is.na(interventionid)) %>%
            pull(sha_id)
          
          candid_sha_id <- ifelse(length(candid_sha_id) == 0, 0, candid_sha_id[1])
          
          # if a previous int match exists, check sha_id
          if (is_int_match & this_sha_id != candid_sha_id){
            print("Different sha_id.")
            # print(paste0("this_sha_id", this_sha_id))
            # print(paste0("candid_sha_id", candid_sha_id))
            reason_msg = paste0(reason_msg, "Different sha_id.")
            
            index <- master_ppl2_loop$indiv_id != matched_indiv_id
            master_ppl2_loop = 
              master_ppl2_loop[index,]
            # 02-17-2021
            indiv_ll <- indiv_ll[index,]
            dist <- dist[index]
            person_match_int = NULL
            
          }else if (is_int_match & this_sha_id != 0){
            person_match_int = master_pdays2 %>%
              filter(indiv_id==matched_indiv_id & date == int_date) 
            print("Successful match via sha_id")
            break
            
          }else {
            ### 2.2.1 If yes, check whether age, sex match ------------------------------------------
            pre_age = master_ppl2_loop %>% 
              filter(indiv_id == matched_indiv_id) %>% pull(age)
            
            if (length(pre_age) > 1 & any(!is.na(pre_age))){
              pre_age <- max(pre_age[!is.na(pre_age)])
            }
            
            pre_sex = master_ppl2_loop %>% 
              filter(indiv_id == matched_indiv_id) %>% pull(sex)
            
            if (length(pre_sex) > 1 & any(!is.na(pre_sex))){
              pre_sex <- unique(pre_sex[!is.na(pre_sex)])
              assert_that(length(pre_sex) == 1,
                          msg = "Different genders for the same person")
            }
            
            this_age <- this_id_df$age
            
            prev_match_list <- master_pdays2 %>%  
              filter(indiv_id==matched_indiv_id &
                       evermatchday == 1)
            
            last_date <- prev_match_list$date[nrow(prev_match_list)]
            age_diff <- ifelse(last_date <= this_date, this_age - pre_age, pre_age - this_age)
            
            # 06/21/2021 age tolerance
            age_match = is.na(pre_age) | (0 <= age_diff & age_diff <= 1)
            sex_match = is.na(pre_sex) | this_id_df$sex == pre_sex
            
            print(paste0("pre_age: ", pre_age))
            print(paste0("this_age: ", this_id_df$age))
            
            print(paste0("pre_sex: ", pre_sex))
            print(paste0("this_sex: ", this_id_df$sex))
            
            print(paste0("last_date: ", last_date))
            print(paste0("this_date: ", this_date))
            
            if(!age_match | !sex_match | is.na(age_match) | is.na(sex_match)){
              #### 2.2.1.1 If no, next closest person  ------------------------------------------
              # if age and sex do not match, find next closest person 
              print("Age and sex do not match for this id, find next closest person.")
              reason_msg = paste0(reason_msg, "A/S not match.", paste())
              
              index <- master_ppl2_loop$indiv_id != matched_indiv_id
              master_ppl2_loop = 
                master_ppl2_loop[index,]
              # 02-17-2021
              indiv_ll <- indiv_ll[index,]
              dist <- dist[index]
              
              # if age and sex DO match, then keep this match
            }else{
              #### 2.2.1.2 If yes, check whether date is taken------------------------------------------
              # check whether date is already taken for this person
              prev_day_match = master_pdays2 %>% 
                filter(indiv_id == matched_indiv_id & 
                         (date == int_date)) %>% 
                pull(evermatchday) == 1
              
              if(!is.na(prev_day_match) & prev_day_match){
                #### 2.2.1.2.1 If yes, next closest person------------------------------------------
                print("Age and sex match but date already taken. Proceed to next nearest person.")
                reason_msg = paste0(reason_msg, "Date taken.")
                
                index <- master_ppl2_loop$indiv_id != matched_indiv_id
                master_ppl2_loop = 
                  master_ppl2_loop[index,]
                # 02-17-2021
                indiv_ll <- indiv_ll[index,]
                dist <- dist[index]
                
                person_match_int = NULL
                
                # if date is available, check last int at  
                # least 5/8 weeks ago (arm specific arm 0 (RO, TO): 35 days, arm 1 (RV, TV):  35 days)
              }else{
                #### 2.2.1.2.2 If no, check wash-out------------------------------------------
                print("Checking wash-out period")
                
                prev_match_list <-  master_pdays2 %>% 
                  filter(indiv_id== matched_indiv_id,
                         evermatchday ==1)
                
                washout_check <- check_washout_arm(prev_match_list = prev_match_list, 
                                                   this_date = int_date,
                                                   mindays0 = 35, 
                                                   mindays1 = 35)
                
                if (!washout_check){
                  ##### 2.2.1.2.2.1 If not pass, next closest person------------------------------------------
                  print("Two matches are two close in time. Proceed to next nearest person.")
                  reason_msg = paste0(reason_msg, "wash-out fail.")
                  
                  index <- master_ppl2_loop$indiv_id != matched_indiv_id
                  master_ppl2_loop = 
                    master_ppl2_loop[index,]
                  indiv_ll <- indiv_ll[index,]
                  dist <- dist[index]
                  person_match_int = NULL
                }else {
                  ##### 2.2.1.2.2.2 If pass, success------------------------------------------
                  person_match_int = master_pdays2 %>%
                    filter(indiv_id==matched_indiv_id & date == int_date) 
                  print("Successful match via wash")
                  break
                }
              } # end checking wash-out period
            } # end checking date taken
           } # NEW
          }

          
        }else{
          ##### 2.2.2 If no, success ------------------------------------------
          # if no previous match, proceed with this match
          person_match_int = master_pdays2 %>%
            filter(indiv_id==matched_indiv_id & date == int_date) 
          print("Successful match via no")
          break
        }
        
        # end check of available person-days   
      }
      # end of distance check loop
    }
    # end of repeat  
  }
  
  if(!is.null(person_match_int)){
    # save variables for matched person 
    person_match_int = person_match_int %>%
      mutate(interventionid = this_id_df$id,
             iid_combined = iid_comb,
             age = this_id_df$age,
             sex = this_id_df$sex,
             dist_to_match = dist[which.min(dist)],
             evermatch = 1,
             evermatchday = 1,
             sha_id = this_id_df$sha_id
      )
    
    # check that only 1 person-day was selected
    assert_that(nrow(person_match_int)==1)
    
    int_matched_list[[i]] = person_match_int
    
    # ............................................
    # assign age and sex to this indiv_id so that 
    # it is factored into subsequent matches 
    # ............................................
    if (!is.null(person_match_int)){
      master_pdays2 = master_pdays2 %>% mutate(
        interventionid = if_else(indiv_id == matched_indiv_id, 
                                unique(person_match_int$interventionid), interventionid),
        evermatch = if_else(indiv_id == matched_indiv_id, 1, evermatch),
        evermatchday = if_else(indiv_id == matched_indiv_id & date == int_date, 1, evermatchday),
        sha_id = if_else(indiv_id == matched_indiv_id & date == int_date, this_sha_id, sha_id)
      )
      
      master_ppl2 = master_ppl2 %>% mutate(
        interventionid = if_else(indiv_id == matched_indiv_id, 
                                unique(person_match_int$interventionid), interventionid),
        age = dplyr::if_else(indiv_id == matched_indiv_id, unique(person_match_int$age), age),
        sex = if_else(indiv_id == matched_indiv_id, unique(person_match_int$sex), sex),
        evermatch = if_else(indiv_id == matched_indiv_id, 1, evermatch),
        sha_id = if_else(indiv_id == matched_indiv_id, this_sha_id, sha_id)
      )
      
    }
    
    
    # end of check for prior match loop
  }
}
toc()


intervention_matched_df = bind_rows(int_matched_list)
intervention_unmatched_df = bind_rows(int_unmatched_list)

# check that there are no duplicated individuals and dates
assert_that(
  intervention_matched_df %>% dplyr::select(indiv_id, date) %>% nrow() ==
    intervention_matched_df %>% dplyr::select(indiv_id, date) %>% distinct() %>% nrow()
)

# ............................................
# Save initial matching output
# ............................................

saveRDS(intervention_matched_df,
        file = namibia_temp_int_match_path)
saveRDS(intervention_unmatched_df,
        file = namibia_temp_int_unmatch_path)

intervention_matched_df = readRDS(namibia_temp_int_match_path)
intervention_unmatched_df = readRDS(namibia_temp_int_unmatch_path)

# check that number of ids matched + unmatched = total
assert_that(
  intervention_matched_df %>% dplyr::select(interventionid) %>% distinct() %>% nrow() +
    nrow(intervention_unmatched_df) == 
    length(intervention$id),
  msg = "Some ids unaccounted for following matching."
)

# ............................................
# ensure that the same person was not selected for more than
# one intervention iid_combined
# ............................................
assert_that(intervention_matched_df %>% 
              group_by(iid_combined, interventionid, date) %>%
              summarise(n=n()) %>% filter(n!=1) %>% nrow() == 0)


# ............................................
# combine intervention and individual data into a single data frame
# each row contains intervention individuals matched to 
# individuals from the GR survey 
# ............................................
intervention_matched_df = intervention_matched_df %>%
  mutate(
    clinic_day = NA,
    int_day = 1,
    indexcase = 0,
    int_trig = NA,
    int_recip = 1
  )

# ............................................
# merge dataset with index to individual matches and
# dataset with intervention to individual matches. 
# ............................................
# check that there isn't any overlap. 
check_overlap_index = all_index_matched %>% 
  filter(indexcase==1) %>% 
  dplyr::select(indiv_id, date) %>% 
  mutate(merged_index = 1)

check_overlap_int = intervention_matched_df %>% 
  dplyr::select(indiv_id, date) %>% 
  mutate(merged_int = 1)

check_overlap = full_join(check_overlap_index, check_overlap_int, 
                          by = c("indiv_id", "date")) %>% 
  mutate(check = merged_index + merged_int)

assert_that(all(is.na(check_overlap$check)))

# create person-day dataset that only includes
# individuals that weren't matched to intervention recipients
matched_int_ids = intervention_matched_df %>%
  dplyr::select(indiv_id, date) %>%
  mutate(drop=1)

all_index_sub = full_join(all_index_matched, matched_int_ids, by = c("indiv_id", "date")) %>%
  filter(is.na(drop)) %>%
  dplyr::select(-drop)

# combine intervention matches into data frame with
# index matches and non-matches
all_index_int = bind_rows(all_index_sub, intervention_matched_df) %>%
                arrange(indiv_id, date, ea)

all_index_int = all_index_int %>% 
                mutate(evermatch = ifelse(evermatch!=1, 0, evermatch))

# check that the number of person-days matches original data frame
assert_that(nrow(person_time) == nrow(all_index_int))

# check for duplicates in individual ids
assert_that(any(duplicated(unique(all_index_int$indiv_id)))==F)

# for a given indiv_id, if match occurred, 
# impute index id, age, sex for all other dates 
all_index_int_matched = all_index_int %>% 
  group_by(indiv_id) %>% 
  mutate(evermatch = max(evermatch, na.rm = T)) %>% 
  filter(evermatch==1) %>% 
  group_by(indiv_id) %>% 
  fill(age, .direction = "downup") %>% 
  fill(sex, .direction = "downup") %>% 
  ungroup()

all_index_int_unmatched = all_index_int %>% 
  group_by(indiv_id) %>% 
  mutate(evermatch = max(evermatch, na.rm = T)) %>% 
  filter(evermatch!=1)

all_index_int_imputed = bind_rows(
  all_index_int_matched, all_index_int_unmatched) %>% 
  arrange(indiv_id, date)

#....................................................
# check that the number of index cases triggering
# intervention = the number of matched index cases
#....................................................
assert_that(
  all_index_int_imputed %>% filter(indexcase==1 & int_trig==1)  %>%
    dplyr::select(indiv_id) %>% distinct() %>% nrow() ==
    all_index_matched %>% filter(indexcase==1) %>%  dplyr::select(indiv_id)  %>%
    distinct() %>% nrow(),
  msg = "# index cases triggering intervention not equal to # matched index cases"
)


saveRDS(all_index_int_imputed, file = namibia_int_match_path)
