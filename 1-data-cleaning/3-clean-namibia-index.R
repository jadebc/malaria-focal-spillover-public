################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Clean index data
################################################
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(raster)
library(geosphere)
library(sp)
library(leaflet)

## NOTES on variables
# Even if an index case has an iid_combined, 
# that doesn’t mean that it was involved in or 
# covered by an actual or "real" intervention.

# If village_covered==1, assume they were covered by an intervention
# based on the stipulations set forth in the protocol. 
# This is the variable used for index coverage for the paper.

#--------------------------------------------
# load raw index data
#--------------------------------------------
index_raw = read_df(read_func = read.dta13,
                    path = namibia_raw_index_path)

intervention = readRDS(namibia_partial_clean_intervention_path)

################################################
# clean index data
################################################
index = index_raw %>%
  
  # keep required variables only
  dplyr::select(-c(nationality, spray, ownnet, sleptnet,
                   preg, travhist, travcountry1, imported,
                   imported_from, recid,
                   region_name,
                   village, villageID,
                   constituency,
                   week, has_gr,
                   ea_type, ea_subtype, sq_m,
                   sq_km, pop_size, occ_clean,
                   comments)) %>%
  
  # rename variables for convenience
  rename(longitude_index = longitude_combined,
         latitude_index = latitude_combined,
         visit_date_index = visit_date,
         iid_combined = iid_combined_c,
         gr_combined_index = gr_combined,
         intarm = intcode,
         ea_index = ea,
         eano_index = ea_no,
         age_index = age,
         sex_index = sex,
         clinic_date_index = date) %>%
 
  # sort by date and gr code 
  arrange(visit_date_index, iid_combined) 

#---------------------------------------------------
# If first_case_yn==1 and outside_studyarea==1
# manually correcting outside_study area coding to 0
# after checking lat long for iid for index and intervention 
#---------------------------------------------------
# list index iids with first_case_yn==1 and outside_studyarea==1
index = index %>%
  mutate(outside_studyarea = 
           ifelse(first_case_yn == 1 & outside_studyarea == 1, 0, outside_studyarea))

assert_that(names(table(index$outside_studyarea[index$first_case_yn==1]))==0)

#---------------------------------------------------
# According to Brooke, first_case_yn==2 means they were in the study area 
# but were not assigned to a "real" intervention ID. 
# They could have been covered by a different intervention 
# ID than the one to which they were assigned, so imputing
# NA for their iid_combined. 
#---------------------------------------------------
index = index %>% mutate(
  iid_combined = ifelse(first_case_yn==2, "", iid_combined)
)

# create unique id for index
index$id = seq(1, nrow(index))


#--------------------------------------------
# manual correction of visit dates that occurred 
# after clinic date based on intervention dataset
# with matching iid
#--------------------------------------------
index$visit_date_index[index$iid_combined=="RO39-4254" & index$first_case_yn==1] =
  unique(intervention$start_date_int[intervention$iid_combined=="RO39-4254"])

index$visit_date_index[index$iid_combined=="TV41-4404" & index$first_case_yn==1] =
  min(intervention$start_date_int[intervention$iid_combined=="TV41-4404"])

index$visit_date_index[index$iid_combined=="RV07-4042" & index$first_case_yn==1] =
  min(intervention$start_date_int[intervention$iid_combined=="RV07-4042"])

index$clinic_date_index[index$iid_combined=="TO10-3975" & index$first_case_yn==1] =
  as.Date("2017-02-02")

check_dates = index %>% 
  filter(first_case_yn==1) %>%
  mutate(date_diff = visit_date_index - clinic_date_index) %>% 
  filter(date_diff < 0) %>%
  dplyr::select(iid_combined, visit_date_index, clinic_date_index)

assert_that(check_dates %>% nrow() == 0)

#--------------------------------------------
# manual id corrections
#--------------------------------------------
index = index %>%
  mutate(intarm = ifelse(iid_combined == "RV32-3944", "RV", intarm)) %>%
  mutate(ea_index = ifelse(iid_combined == "TV16-3912", "TV16", ea_index)) %>%
  mutate(intarm = ifelse(iid_combined == "TV16-3912", "TV", intarm)) %>% 
  mutate(ea_index = ifelse(ea_index == "TO32" & intarm=="RV", "RV32", ea_index)) %>% 
  mutate(iid_combined = ifelse(iid_combined == "RO24-6721", "RV41-6721", iid_combined)) %>% 
  mutate(iid_combined = ifelse(iid_combined == "TO37-4923", "RO52-4923", iid_combined)) %>% 
  mutate(iid_combined = ifelse(iid_combined == "RV41-6721", "TV41-6721", iid_combined)) %>% 
  mutate(ea_index = ifelse(iid_combined=="RO15-4242" & ea_index == "", "RO15", ea_index))

intervention = intervention  %>%
  # mutate(iid_combined = ifelse(iid_combined == "RV32-3944", "TO32-3944", iid_combined)) %>% 
  mutate(ea = ifelse(iid_combined == "RV32-4225", "RV32", ea))


# check that first two letters match in intarm and iid and ea
index_check = index %>% 
  dplyr::select(ea_index, iid_combined, intarm) %>%
  mutate(ea_prefix = substr(ea_index, 1, 2),
         iid_prefix = substr(iid_combined, 1,2)) 

assert_that(index_check %>% 
  filter(ea_prefix!="") %>%
  filter(intarm!=ea_prefix) %>%  nrow() ==0)

assert_that(index_check %>% 
  filter(iid_prefix!="") %>%
  filter(intarm!=iid_prefix) %>%  nrow() ==0)


################################################
# check EA against shape file
################################################
# read in shape file 
ea_shape <- read_shp(dsn = namibia_shapefile_dsn,
                     layer = namibia_shapefile_layer)

# drop if missing longitude
index_sub = index %>% filter(!is.na(longitude_index))

# drop if no ea
index_sub = index_sub %>% filter(ea_index != "")

# create a SpatialPointsDataFrame for intervention data
index_SPDF <- SpatialPointsDataFrame(
  coords = index_sub[, c("longitude_index", "latitude_index")],
  data = index_sub[,c("intarm","ea_index", "id")],
  proj4string = CRS("+init=epsg:4326"))

# confirm that point data and polygon data use same CRS
crs(index_SPDF) = crs(ea_shape)

# obtain polygon data for each point 
# each row is a household in the GR survey
point_poly = over(index_SPDF, ea_shape)
assert_that(nrow(point_poly) == nrow(index_SPDF))
point_poly = point_poly %>% 
  mutate(ea_index = index_sub$ea_index,
         intarm = index_sub$intarm,
         iid_combined = index_sub$iid_combined,
         id = index_sub$id) %>%
  rename(shapefile_ID = ID2,
         shapefile_ea = Eas_StudyA,
         shapefile_intarm = Eas_Stud_1)

# check for discrepancies
point_poly = point_poly %>% mutate(
  diff_ea = ifelse(shapefile_ea != ea_index, "Yes", "No")
)

# create list of discrepancies
discrep_ea = point_poly %>%
  filter(diff_ea == "Yes" | is.na(diff_ea)) %>%
  group_by(ea_index, intarm, iid_combined, diff_ea) %>%
  rename(ea = ea_index) %>% 
  summarise(n = n())

# save discrepant EAs
# write.csv(discrep_ea, file = paste0(here::here(), 
#                                    "/1-data-cleaning/temp/index_discrepant_eas.csv"), row.names = F)

#--------------------------------------------
# merge discrepant eas back onto index data
# and drop GPS values for those individuals
#--------------------------------------------
point_poly_m = point_poly %>% 
  dplyr::select(iid_combined, id, diff_ea)

index_clean = full_join(index, point_poly_m, 
                               by = c("iid_combined", "id")) %>%
  mutate(longitude_index = ifelse(diff_ea=="Yes" | is.na(diff_ea), NA, longitude_index),
         latitude_index = ifelse(diff_ea=="Yes" | is.na(diff_ea), NA, latitude_index)) 

# % of index cases who received interventions
# whose GPS are likely incorrect
prop.table(table(index_clean$diff_ea))
summary(as.factor(index_clean$diff_ea))

index_clean = index_clean %>%
  dplyr::select(-diff_ea) %>%
  dplyr::select(ea_index, id, everything())

#--------------------------------------------
# merge discrepant eas back onto intervention data
# to create map
#--------------------------------------------
mapdata = index_sub %>% 
  dplyr::select(
    ea_index, intarm, iid_combined,
    longitude_index, latitude_index) %>%
  rename(ea = ea_index)

mapdata = left_join(discrep_ea, mapdata, by = c("ea", "intarm", "iid_combined"))

ea_sub = ea_shape[ea_shape$Eas_StudyA %in% c(mapdata$ea),]
ea_sub$Eas_StudyA = as.factor(ea_sub$Eas_StudyA)
ea_sub$Eas_StudyA = droplevels(ea_sub$Eas_StudyA)

# color palette for EA intervention assignment
rainbow_pal = rainbow(length(unique(ea_sub$Eas_StudyA)))
int_pal <- colorFactor(rainbow_pal, ea_sub$Eas_StudyA,
                       na.color = "#E8EAEA")

mapdata_SPDF <- SpatialPointsDataFrame(
  coords = mapdata[,c("longitude_index", "latitude_index")],
  data = mapdata[,c("intarm", "ea")],
  proj4string = CRS("+init=epsg:4326"))

# map of discrepant EAs
eamap =  leaflet() %>% addProviderTiles("CartoDB.Positron") %>% 
  # set zoom center point and level
  setView(lng = 23.8, lat = -17.9, zoom = 08) %>%
  
  # add ea polygons
  addPolygons(
    data=ea_shape,
    weight = 1,
    fillOpacity = 1,
    fillColor = ~int_pal(ea_sub$Eas_StudyA),
    col = "black"
  )  %>%
  
  # legend for intervention arm
  addLegend(pal = int_pal,
            title = "EA - shape file",
            values = ea_sub$Eas_StudyA,
            opacity = 1) %>%
  
  # add markers for individuals with discrepant EAs
  addCircleMarkers(data=mapdata_SPDF,
                   color="black",
                   fillColor = ~int_pal(mapdata_SPDF$ea),
                   fillOpacity = 1,
                   weight = 0.2,
                   radius = 2) 

# eamap

# check ea missing index_nonint and other index_nonint

studyarea_palette <- colorFactor(c("#8d99ae", rep("#fb8500", 4)),
                                 ea_shape$Arm2)

index_nonint_miss <- index %>% filter(first_case_yn != 1 &
                                      ea_index == "")
index_nonint_ea <- index %>% filter(first_case_yn != 1 &
                                    ea_index != "")


index_nonint_miss_SPDF <- SpatialPointsDataFrame(
  coords = index_nonint_miss[,c("longitude_index", "latitude_index")],
  data = index_nonint_miss,
  proj4string = CRS("+init=epsg:4326"))

index_nonint_ea_SPDF <- SpatialPointsDataFrame(
  coords = index_nonint_ea[,c("longitude_index", "latitude_index")],
  data = index_nonint_ea,
  proj4string = CRS("+init=epsg:4326"))


leaflet() %>% 
  addProviderTiles("CartoDB.Positron")  %>% 
  # set zoom center point and level
  setView(lng = 23.6, lat = -17.9, zoom = 9) %>% 
  # add ea polygons
  addPolygons(data = ea_shape,
              weight = 1.5,
              color = "black",
              fillColor = ~studyarea_palette(Arm2)
  ) %>% 
  addCircleMarkers(data = index_nonint_miss_SPDF,
                   color = "red",
                   fillOpacity = 1,
                   radius = 2) %>% 
  addCircleMarkers(data = index_nonint_ea_SPDF,
                   color = "blue",
                   fillOpacity = 1,
                   radius = 2)




################################################
# for index cases with missing GPS

# 14% of triggering index cases missing GPS (n=49)
# 32% of first_case_yn==0 missing GPS (n=242)
# 73% of first_case_yn==2 missing GPS (n=904) -- FOR NOW EXCLUDING THEM 
################################################

#----------------------------------------------
# list of intervention iid_combineds to exclude from 
# imputation because of outlier GPS or there are 
# fewer than 3 people for an intervention iid
#----------------------------------------------
dist_list = list()
iid_list = unique(intervention$iid_combined[!is.na(intervention$iid_combined)])

for(i in 1:length(iid_list)) {
  dm = dist(intervention[intervention$iid_combined == iid_list[i],
                         c("longitude_int", "latitude_int")])
  dist_list[[i]] = max(dm, na.rm = T)
}

max_dist_df = unlist(dist_list)


excludes = c(
  # outlier
  iid_list[which(max_dist_df > 0.12)],
  
  # fewer than 3 people for intervention
  intervention %>% group_by(iid_combined) %>% 
    summarise(n=n()) %>% filter(n<3) %>% 
    dplyr::select(iid_combined) %>% pull()
)

#----------------------------------------------
# for index cases that triggered intervention, 
# impute coordinates at the center of the intervention GPS points 
# for the intervention GPS points with the same iid
#----------------------------------------------
index_missing_GPS_trig_sub = index_clean %>% 
  filter(is.na(longitude_index)) %>%
  filter(first_case_yn==1)

# number of intervention triggering index cases
# with missing lat long: 
nrow(index_missing_GPS_trig_sub)

# calculate center of intervention area for each iid_combined
intervention_sub = intervention[!intervention$iid_combined %in% excludes,]
# 06/22/2021 try median
centers = intervention_sub %>% 
  group_by(iid_combined) %>%
  dplyr::select(iid_combined, longitude_int, latitude_int) %>% 
  summarise(
    longitude_center = median(longitude_int, na.rm=T),
    latitude_center = median(latitude_int, na.rm=T)
  ) %>%
  filter(!is.na(longitude_center))

index_missing_GPS_trig = left_join(
  index_missing_GPS_trig_sub, centers, by = "iid_combined") %>%
  mutate(
    longitude_index_imp = ifelse(is.na(longitude_index), longitude_center, longitude_index),
    latitude_index_imp = ifelse(is.na(latitude_index), latitude_center, latitude_index)
  )

assert_that(nrow(index_missing_GPS_trig) == nrow(index_missing_GPS_trig_sub))

# percentage of index cases with newly imputed lat/long
sum(!is.na(index_missing_GPS_trig$longitude_index_imp))/nrow(index_missing_GPS_trig)

# number of index cases with newly imputed lat/long
sum(!is.na(index_missing_GPS_trig$longitude_index_imp))

# number of index cases without lat/long
sum(is.na(index_missing_GPS_trig$longitude_index_imp))

# make sure that for all index cases without lat/long,
# intervention has no lat/long too
check = index_missing_GPS_trig$iid_combined[is.na(index_missing_GPS_trig$longitude_index_imp)]

if(intervention_sub %>% filter(iid_combined %in% check) %>% nrow() >0){
  assert_that(is.na(
    intervention_sub %>% filter(iid_combined %in%
                                  check) %>%
      dplyr::select(longitude_int) %>% unique() %>% pull()
  ))
}

#----------------------------------------------
# index cases that did not trigger intervention 
# and that have an iid_combined,
# impute coordinates at random location inside
# of polygon of GPS points for an intervention 
# N = 66 index cases
#----------------------------------------------
index_missing_GPS_nontrig_sub = index_clean %>% 
  filter(is.na(longitude_index)) %>%
  filter(first_case_yn==0) %>%
  filter(iid_combined!="")

# drop observations with missing lat/long
random_loc = intervention %>% 
  group_by(iid_combined) %>%
  dplyr::select(iid_combined, longitude_int, latitude_int)  %>%
  filter(!is.na(longitude_int)) 

random_loc = random_loc[!random_loc$iid_combined %in% excludes,]

# create polygon for each intervention 
polygon_list = list() 

for (i in 1:length(unique(random_loc$iid_combined))) {
  block_coords = random_loc[random_loc$iid_combined == unique(random_loc$iid_combined)[i], 
                            c("longitude_int", "latitude_int")]
  ch <- chull(block_coords)
  coords <- block_coords[c(ch, ch[1]), ]
  p = Polygon(coords)
  polygon_list[[i]] = p
}

# spatial_polys = SpatialPolygons(polygon_list) #- this ran before; not sure why not working now
# iid_df = data.frame(iid_combined = unique(random_loc$iid_combined))
# intervention_polygons = SpatialPolygonsDataFrame(spatial_polys, iid_df)
# 
# # plot polygons
# spplot(intervention_polygons, col=1:326)

# randomly sample from each polygon 
rand_list = list()
for(i in 1:length(polygon_list)){
  if(polygon_list[[i]]@area>0){
    rand = spsample(x = polygon_list[[i]], n=1, type = "random") 
    rand_list[[i]] = rand %>% coordinates() %>% 
      as.data.frame() %>%
      mutate(iid_combined = unique(random_loc$iid_combined)[i])
  } 
}

rand_df = bind_rows(rand_list) %>%
  rename(longitude_rand = x,
         latitude_rand = y)

# impute 
index_missing_GPS_nontrig = left_join(index_missing_GPS_nontrig_sub, rand_df, by = "iid_combined") %>%
  mutate(longitude_index_imp = ifelse(is.na(longitude_index), longitude_rand, longitude_index),
         latitude_index_imp = ifelse(is.na(latitude_index), latitude_rand, latitude_index))
assert_that(nrow(index_missing_GPS_nontrig) == nrow(index_missing_GPS_nontrig_sub))

# percentage of index cases with newly imputed lat/long
sum(!is.na(index_missing_GPS_nontrig$longitude_index_imp))/nrow(index_missing_GPS_nontrig)

# number of index cases with newly imputed lat/long
sum(!is.na(index_missing_GPS_nontrig$longitude_index_imp))

#----------------------------------------------
# combine datasets after imputation
#----------------------------------------------
# imputed ids
imputed_ids = c(index_missing_GPS_trig$id,
                index_missing_GPS_nontrig$id)

# grab ids not in the imputed datasets
non_imputed_data = index_clean[!index_clean$id %in% imputed_ids,] %>%
  mutate(longitude_index_imp = longitude_index,
         latitude_index_imp = latitude_index)

assert_that(nrow(index_clean) == nrow(non_imputed_data)+nrow(index_missing_GPS_nontrig)+
  nrow(index_missing_GPS_trig))

index_imputed = bind_rows(
  index_missing_GPS_nontrig,
  index_missing_GPS_trig,
  non_imputed_data
) %>%
  mutate(longitude_index = longitude_index_imp, 
         latitude_index = latitude_index_imp) %>%
  dplyr::select(-c(longitude_center, latitude_center, 
                   longitude_rand, latitude_rand)) %>%
  dplyr::select(id, ea_index, iid_combined,
                longitude_index, latitude_index, everything()) %>%
  arrange(id)

#----------------------------------------------
# manual id correction
#----------------------------------------------

#--------------------------------------------
# save clean index data
#--------------------------------------------
saveRDS(index_imputed, namibia_clean_index_path)
saveRDS(intervention, namibia_clean_intervention_path)


