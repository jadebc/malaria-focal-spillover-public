################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Clean intervention data
################################################
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(raster)

#--------------------------------------------
# load raw intervention data
#--------------------------------------------
intervention_raw = read_df(read_func = read.dta13,
                           path = namibia_raw_intervention_path)

# 06/24/2021 drop obs with indexcase = Yes
intervention_raw <- intervention_raw %>% filter(indexcase == "No")

#--------------------------------------------
# clean intervention data
#--------------------------------------------
intervention = intervention_raw %>%
  
  # keep required variables only
  dplyr::select(c(gpslongitude, gpslatitude, start_date, 
                  individual_combined, gr_combined,
                  iid_combined, 
                  intarm, age, sex,
                  c_questionnaire, c_bloodrapidtesting, c_antimalarialmed ,
                  rdt_result, lamp_result)) %>%
  
  # rename variables for convenience
  rename(longitude_int = gpslongitude,
         latitude_int = gpslatitude,
         start_date_int = start_date) %>% 
  
  # create ea variable from iid per Brooke's advice 
  mutate(ea = substr(iid_combined, 1, 4)) %>%
  
  # sort by date and gr code
  arrange(start_date_int, gr_combined)

# create unique id for intervention
intervention$id = seq(1, nrow(intervention))

#--------------------------------------------
# manual id corrections 
#--------------------------------------------
intervention = intervention %>%
  mutate(intarm = ifelse(iid_combined == "RV01-4396", "RV", intarm)) %>%
  mutate(intarm = ifelse(iid_combined == "RV46-6502", "RV", intarm)) %>%
  mutate(intarm = ifelse(iid_combined == "RV48-5782", "RV", intarm)) %>%
  mutate(intarm = ifelse(iid_combined == "TO38-4732", "TO", intarm)) %>% 
  mutate(ea = ifelse(ea == "TO32", "RV32", ea))

################################################
# check EA against shape file
################################################
# read in shape file 
ea_shape <- read_shp(dsn = namibia_shapefile_dsn,
                     layer = namibia_shapefile_layer)

# drop if missing longitude
intervention_sub = intervention %>% filter(!is.na(longitude_int))

# create a SpatalPointsDataFrame for intervention data
intervention_SPDF <- SpatialPointsDataFrame(
  coords = intervention_sub[, c("longitude_int", "latitude_int")],
  data = intervention_sub[,c("intarm","ea", "id")],
  proj4string = CRS("+init=epsg:4326"))

# confirm that point data and polygon data use same CRS
crs(intervention_SPDF)
crs(ea_shape)

# obtain polygon data for each point 
# each row is a household in the GR survey
point_poly = over(intervention_SPDF, ea_shape)
assert_that(nrow(point_poly) == nrow(intervention_SPDF))
point_poly = point_poly %>% 
  mutate(ea = intervention_sub$ea,
         intarm = intervention_sub$intarm,
         iid_combined = intervention_sub$iid_combined,
         id = intervention_sub$id) %>%
  rename(shapefile_ID = ID2,
         shapefile_ea = Eas_StudyA,
         shapefile_intarm = Eas_Stud_1)

# check for discrepancies
point_poly = point_poly %>% mutate(
  diff_ea = ifelse(shapefile_ea != ea, "Yes", "No")
)


# create list of discrepancies
discrep_ea = point_poly %>%
  filter(diff_ea == "Yes" | is.na(diff_ea)) %>%
  group_by(iid_combined, diff_ea) %>%
  summarise(n = n())

# save discrepant EAs
write.csv(discrep_ea, file = paste0(here::here(),
       "/1-data-cleaning/1-namibia/temp/intervention_discrepant_eas.csv"), row.names = F)


#--------------------------------------------
# merge discrepant eas back onto intervention data
# and drop GPS values for those individuals
#--------------------------------------------
# <1% of observations missing GPS 
prop.table(table(is.na(intervention$longitude_int)))

point_poly_m = point_poly %>% 
  dplyr::select(iid_combined, id, diff_ea)

intervention_clean = full_join(intervention, point_poly_m, 
              by = c("iid_combined", "id")) %>%
  mutate(longitude_int = ifelse(diff_ea=="Yes" | is.na(diff_ea), NA, longitude_int),
         latitude_int = ifelse(diff_ea=="Yes" | is.na(diff_ea), NA, latitude_int)) 

# % of individuals who received interventions
# whose GPS are likely incorrect
prop.table(table(intervention_clean$diff_ea))

intervention_clean = intervention_clean %>%
  dplyr::select(-diff_ea) %>%
  dplyr::select(ea, id, everything())

# drop obs with missing lat/long
intervention_clean = intervention_clean %>% filter(!is.na(longitude_int) & 
                                                   !is.na(latitude_int))

# identify same person if: age, sex are same, iid_comb different, distance <= 20m, outside wash-out

intervention_clean$sha_id <- 0
count_link <- 0

# View(intervention_clean %>%  filter(id %in% c(6011, 8223, 8224, 8225)))

for (i in 1:nrow(intervention_clean)){
  cat("loop", i, "\n")
  
  this_df <- intervention_clean[i, ]
  
  # if already linked, jump to next person
  if (this_df$sha_id != 0){
    next
  }
  
  this_ll <- data.frame(long = this_df$longitude_int,
                        lat = this_df$latitude_int)
  
  this_iid_comb <- this_df$iid_combined
  this_date <- this_df$start_date_int
  
  # candidates are those with same age, sex, 
  # but different iid_combined and outside the wash-out period
  # other_df <- intervention_clean %>% 
  #                 filter((age >= this_df$age - 1) & 
  #                 (age <= this_df$age + 1) &
  #                 sex == this_df$sex &
  #                 iid_combined != this_iid_comb &
  #                 (start_date_int >= this_date + 28 |
  #                 start_date_int <= this_date - 28))
  
  other_df <- intervention_clean %>% 
    filter((((start_date_int < this_date) & 
              (age <= this_df$age) & 
              (age >= this_df$age - 1)) |
             ((start_date_int > this_date) & 
                (age >= this_df$age) & 
                (age <= this_df$age + 1))) &
           sex == this_df$sex &
           iid_combined != this_iid_comb &
           (start_date_int >= this_date + 28 |
              start_date_int <= this_date - 28))
  
  
  # if no avalable people, next
  if (nrow(other_df) == 0){
    next
  }
  
  other_ll <- data.frame(long = other_df$longitude_int,
                         lat = other_df$latitude_int)
  
  dist <- distGeo(this_ll[, c("long", "lat")],
                  other_ll[, c("long", "lat")])
  
  p <- which(dist <= 20)
  # if no people within 20m circle, next
  if (length(p) == 0){
    next
  }

  sub_other_df <- other_df[p,]
  sub_other_df$dist <- dist[p]
  sub_other_df <- arrange(sub_other_df, dist)
  
  for (j in 1:nrow(sub_other_df)){
    candidate <- sub_other_df[j, ]
    # check previous link
    if (candidate$sha_id != 0){
      candid_group <- intervention_clean %>% 
                              filter(sha_id == candidate$sha_id)
      candid_group_iid_comb <- candid_group$iid_combined
      
      candid_group_ll <- data.frame(long = candid_group$longitude_int,
                                    lat = candid_group$latitude_int)
      
      candid_group_dist <- distGeo(this_ll[, c("long", "lat")],
                                   candid_group_ll[, c("long", "lat")])
      
      if (this_iid_comb %in% candid_group_iid_comb | max(candid_group_dist) > 20 ){
        next
      }else{
        # link success
        intervention_clean <- intervention_clean %>% 
          mutate(sha_id = ifelse(id == this_df$id, candidate$sha_id, sha_id))
        print("link success")
        count_link = count_link + 1
        break
      }
    }else{
      # link success
      link_id <- c(this_df$id, candidate$id)
      intervention_clean <- intervention_clean %>% 
                              mutate(sha_id = ifelse(id %in% link_id, this_df$id, sha_id))
      print("link success")
      count_link = count_link + 1
      break
    }
  }
}

temp <- intervention_clean %>% filter(sha_id != 0)
# 
temp2 <- temp %>%
  group_by(sha_id) %>%
  mutate(occurrences = n()) %>%
  ungroup()

#--------------------------------------------
# save clean intervention data
#--------------------------------------------
saveRDS(intervention_clean, namibia_partial_clean_intervention_path)


