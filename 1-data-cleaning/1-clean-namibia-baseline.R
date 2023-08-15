################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Clean baseline geographical reconaissance data
################################################

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(leaflet)
library(sp)

#--------------------------------------------
# Read in datasets
#--------------------------------------------

gr_raw = read_df(read_func = read.csv,
                 path = namibia_raw_gr_path)

index_raw = read_df(read_func = read.dta13,
                    path = namibia_raw_index_path)

# read in shape file 
ea_shape <- read_shp(dsn = namibia_shapefile_dsn,
                     layer = namibia_shapefile_layer)

#--------------------------------------------
# subset to relevant variables
#--------------------------------------------
gr = gr_raw %>% 
  
  # ADD: Filter on ea_actual to ensure they are in the study area!!!!
  # waiting for brooke to clarify where this variable is 
  
  dplyr::select(
    latitude, longitude, 
    hh_id, hh_pop,
    Eas_Stud29, Eas_Stud30
  ) %>%
  rename(gr_ea = Eas_Stud29,
         gr_intarm = Eas_Stud30) 

# assign new unique id
gr$new_hhid = seq(1, nrow(gr))


#--------------------------------------------
# if HH==0, impute population as mean hh size
# for
#--------------------------------------------
mean_hh_size = mean(gr$hh_pop)

gr = gr %>% mutate(hh_pop = ifelse(hh_pop==0, mean_hh_size, hh_pop))

#--------------------------------------------
# create a SpatialPointsDataFrame for GR data
#--------------------------------------------
namibia_GR_SPDF <- SpatialPointsDataFrame(
  coords = gr[,c("longitude", "latitude")],
  data = gr[,c("hh_id","hh_pop", "gr_intarm")],
  proj4string = CRS("+init=epsg:4326"))


#--------------------------------------------
# reassign intervention arm and EA based on 
# location within shape file
# https://cran.r-project.org/web/packages/sp/vignettes/over.pdf
# http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS2_MergingSpatialData_part1_Joins.html
# https://stackoverflow.com/questions/3647744/intersecting-points-and-polygons-in-r
#--------------------------------------------
# confirm that point data and polygon data use same CRS
raster::crs(namibia_GR_SPDF)
raster::crs(ea_shape)

# set same CRS
raster::crs(ea_shape) = raster::crs(namibia_GR_SPDF)

# obtain polygon data for each point 
# each row is a household in the GR survey
point_poly = over(namibia_GR_SPDF, ea_shape)
assert_that(nrow(point_poly) == length(gr$new_hhid))
point_poly = point_poly %>% 
  mutate(new_hhid = gr$new_hhid) %>%
  rename(shapefile_ID = ID2,
         shapefile_ea = Eas_StudyA,
         shapefile_intarm = Eas_Stud_1)

# merge point poly data back onto main GR data frame
gr_shapefile = full_join(gr, point_poly, by = c("new_hhid"))
assert_that(nrow(gr_shapefile) == nrow(gr))

# manual corrections to eas in the GR but not the shape file
gr_shapefile = gr_shapefile %>%
  mutate(gr_ea = as.character(gr_ea),
         gr_intarm = as.character(gr_intarm)) %>%
  mutate(
    gr_ea = ifelse(new_hhid == 2139 & shapefile_ea == "TV16" & shapefile_ID == "211", "TV16", gr_ea),
    gr_intarm = ifelse(new_hhid == 2139 & shapefile_ea == "TV16" & shapefile_ID == "211", "TV", gr_intarm)
  )

# manual corrections to eas in the GR but not the shape file
gr_shapefile = gr_shapefile %>%
  mutate(
    gr_ea = ifelse(new_hhid == 112 & gr_ea == "RO15" & shapefile_ID == "201", "", gr_ea),
    gr_intarm = ifelse(new_hhid == 112 & gr_ea == "RO15" & shapefile_ID == "201", "", gr_intarm),

    gr_ea = ifelse(new_hhid == 2833 & gr_ea == "TO02" & shapefile_ID == "238", "", gr_ea),
    gr_intarm = ifelse(new_hhid == 2833 & gr_ea == "TO02" & shapefile_ID == "238", "", gr_intarm),

    gr_ea = ifelse(new_hhid == 3840 & gr_ea == "TO38" & shapefile_ID == "426", "", gr_ea),
    gr_intarm = ifelse(new_hhid == 3840 & gr_ea == "TO38" & shapefile_ID == "426", "", gr_intarm),
    
    gr_ea = ifelse(new_hhid == 2206 & shapefile_ea == "TV16" & shapefile_ID == "211", "TV16", gr_ea),
    gr_intarm = ifelse(new_hhid == 2206 & shapefile_ea == "TV16" & shapefile_ID == "211", "TV", gr_intarm)
  )

# check that all hhs with a missing EA / intarm from the shape file
# also are missing EA / intarm from the GR survey
assert_that(all(names(table(gr_shapefile$gr_ea[is.na(gr_shapefile$shapefile_ea)]))==""))


# check that no eas are discrepant 
assert_that(gr_shapefile %>% 
  dplyr::select(new_hhid, shapefile_ID, gr_ea, gr_intarm, shapefile_ea, shapefile_intarm) %>% 
  mutate(gr_ea = as.character(gr_ea),
         shapefile_ea = as.character(shapefile_ea)) %>%
  filter(gr_ea!=shapefile_ea) %>%
  nrow() == 0)

# drop extra ea / intarm variables
gr_clean = gr_shapefile %>%
  mutate(ea = gr_ea, 
         intarm = gr_intarm) %>%
  mutate(ea = ifelse(ea =="", NA, ea),
         intarm = ifelse(intarm=="", NA, intarm)) %>%
  mutate(studyarea = ifelse(ea=="", 0, 1)) %>%
  dplyr::select(-c(gr_ea, gr_intarm, shapefile_ea, shapefile_intarm))


#--------------------------------------------
# create map of study area
#--------------------------------------------

# create a SpatalPointsDataFrame for GR data
namibia_GR_clean_SPDF <- SpatialPointsDataFrame(
  coords = gr_clean[,c("longitude", "latitude")],
  data = gr_clean[,c("hh_id","hh_pop", "intarm")],
  proj4string = CRS("+init=epsg:4326"))

# color palette for EA intervention assignment
int_pal <- colorFactor("BrBG", ea_shape$Eas_Stud_1, 
                       na.color = "#E8EAEA")

int_pal2 <- colorFactor("BrBG", gr_clean$intarm, 
                        na.color = "red")

eamap =  leaflet() %>% addProviderTiles("CartoDB.Positron") %>% 
  # set zoom center point and level
  setView(lng = 23.8, lat = -17.9, zoom = 08) %>%
  
  # add ea polygons
  addPolygons(
    data=ea_shape, 
    weight = 1,
    fillOpacity = 1,
    fillColor = ~int_pal(ea_shape$Eas_Stud_1),
    col = "black"
  )  %>%
  
  # legend for intervention arm 
  addLegend(pal = int_pal,
            title = "Arm",
            values = ea_shape$Eas_Stud_1,
            opacity = 1) %>%
  
  # add markers for each household
  addCircleMarkers(data=namibia_GR_clean_SPDF,
                   color="black",
                   weight = 0.2,
                   fillOpacity = 1,
                   fillColor = ~int_pal2(gr_clean$intarm),
                   radius = 2)   


eamap

#--------------------------------------------
# create individual level dataset 
#--------------------------------------------
gr_clean = gr_clean %>% 
  dplyr::select(new_hhid, ea, intarm, studyarea, longitude, latitude, hh_pop)

gr_study = gr_clean  %>% filter(studyarea==1)

# drop households outside study area
gr_study = gr_study %>% mutate(hh_pop = round(hh_pop, 0))
imputed = uncount(data = gr_study, weights = hh_pop)

# check that number of households is correct
assert_that(sum(gr_study$hh_pop) == nrow(imputed))

# create individual ids
rownames(imputed) = NULL
imputed = imputed %>% 
  group_by(new_hhid) %>%
  mutate(personid = seq_along(new_hhid)) %>%
  mutate(indiv_id = paste0(new_hhid, "_", personid)) %>%
  dplyr::select(new_hhid, personid, indiv_id, everything()) %>%
  dplyr::select(-studyarea)
  

#--------------------------------------------
# save
#--------------------------------------------
saveRDS(gr_clean, file = namibia_clean_gr_path)
saveRDS(imputed, file = namibia_indiv_imputed_path)




