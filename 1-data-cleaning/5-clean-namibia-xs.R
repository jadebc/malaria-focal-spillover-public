################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Clean cross-sectional data
################################################
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))

# Load data ---------------------------------------------------------------
xs_raw <- read.dta13(namibia_raw_xs_path)
xs_gps <- read.dta13(paste0(box_shared_path,"Datasets/2017 x-sectional/XS_HH_cleaned_with_samplestage_for_misuk.dta"))
cohorts <- readRDS(namibia_df_all_path)
ea_level_raw <- read_df(read.csv, namibia_raw_ea_level_path)

# Process XS data ---------------------------------------------------------------
xs_sub = xs_raw %>% 
  dplyr::select(hhid, eaid, iid, arm, date, age, gender, rdtresult, hsrdtresult,
                qPCRposneg, traveled_recently, occupation, sleep_net,slept_outdoors,
                slept_net_yest)

xs_gps_ll = xs_gps %>% 
  dplyr::select(hhid, arm, latitude, longitude)

xs_clean = left_join(xs_sub, xs_gps_ll, by = c("hhid","arm")) %>% 
  arrange(hhid) 

assert_that(nrow(xs_clean) == nrow(xs_sub)) 

xs_clean = xs_clean %>% filter(!is.na(longitude)) 
xs_clean$xsid = seq(1, nrow(xs_clean), 1)

# get EA for xs dataset ---------------------------------------------------------------
# read in shape file 
ea_shape <- rgdal::readOGR(dsn = namibia_shapefile_dsn,
                     layer = namibia_shapefile_layer)

xs_hh <- xs_clean %>% dplyr::select(hhid, eaid, arm, longitude, latitude) %>% 
  distinct()

# create a SpatialPointsDataFrame for intervention data
xs_SPDF <- SpatialPointsDataFrame(
  coords = xs_hh[, c("longitude", "latitude")],
  data = xs_hh[,c("hhid","eaid","arm")],
  proj4string = CRS("+init=epsg:4326"))

# confirm that point data and polygon data use same CRS
raster::crs(xs_SPDF) = raster::crs(ea_shape)

# obtain polygon data for each point 
# each row is a household in the GR survey
point_poly = over(xs_SPDF, ea_shape)
assert_that(nrow(point_poly) == nrow(xs_SPDF))
point_poly = point_poly %>% 
  mutate(hhid = xs_hh$hhid) %>%
  dplyr::rename(ea = Eas_StudyA,
         shapefile_intarm = Eas_Stud_1) %>% 
  dplyr::select(hhid, ea, shapefile_intarm)

# merge onto main xs dataset
xs_hh_checked <- full_join(xs_hh, point_poly, by = "hhid") %>% 
  mutate(discrepant = as.factor(case_when(shapefile_intarm!=arm ~ "Yes",
                                          shapefile_intarm==arm ~ "No",
                                          is.na(shapefile_intarm) ~ "Yes"))) 
  
xs_correct <- full_join(xs_clean, xs_hh_checked %>% 
                          dplyr::select(hhid, ea, discrepant), by = "hhid") %>% 
# drop if arm is discrepant, lat/long might be wrong
  filter(discrepant == "No") %>% 
  dplyr::select(-discrepant)


# create distance-based covariates ---------------------------------------------------------------
make_dist_covs <- function(hh_id, min_dist, max_dist){
  print(hh_id)
  
  # subset to single cross-sectionally sampled household
  xs_ll <- xs_correct %>% filter(hhid==hh_id) %>% ungroup() %>%
    dplyr::select(latitude, longitude) %>% 
    distinct() 
  
  hh_arm <- xs_correct %>% filter(hhid==hh_id) %>% 
    dplyr::select(arm) %>% distinct() %>% pull()
    
  # population within distance d
  dist_pop_vector <- distGeo(xs_ll[, c("longitude", "latitude")], 
                             cohorts[, c("longitude", "latitude")]) 
  
  dist_popN <- length(dist_pop_vector[dist_pop_vector>min_dist & 
                                        dist_pop_vector<=max_dist])
  
  # indicator for residing within distance d from a treated person
  intervention_recips = cohorts %>% filter(int_recip==1)
  
  # filter to dates prior to xs survey
  xs_date <- xs_correct %>% filter(hhid==hh_id) %>% dplyr::select(date) %>% 
    distinct() %>%  pull(date) %>% as.Date()
  intervention_recips = intervention_recips %>% filter(date < xs_date)
  
  dist_tx <- distGeo(xs_ll[, c("longitude", "latitude")], 
                     intervention_recips[, c("longitude", "latitude")])
  
  dist_tx_d <- intervention_recips %>% 
    mutate(dist_to_xs = dist_tx) %>% 
    filter(dist_to_xs > min_dist & dist_to_xs <= max_dist)
  

  if(nrow(dist_tx_d)>0){
    match_arm <- dist_tx_d$intarm[dist_tx_d$intarm==hh_arm]
    
    # proportion of treated individuals in distance d with the same treatment 
    if(length(match_arm)>0){
      prop_arm <- prop.table(table(droplevels(match_arm))) %>% as.numeric()
      assert_that(length(prop_arm)==1)
    }else{
      prop_arm <- 0    
    }
    
    # number of treated individuals within distance d 
    dist_Ntx <- nrow(dist_tx_d)
    
    # number of treated individuals with same treatment assignment within distance d 
    dist_Ntx_sametx <- dist_tx_d %>% filter(intarm == hh_arm) %>% nrow()
    
    # number of treated individuals with different treatment assignment within distance d 
    dist_Ntx_difftx <- dist_tx_d %>% filter(intarm != hh_arm) %>% nrow()
    
  if(dist_popN < dist_Ntx) dist_Ntx = dist_popN
    
    out <- data.frame(
      hhid = hh_id, 
      dist_popN = dist_popN,
      dist_arm = prop_arm,
      dist_Ntx = dist_Ntx,
      dist_Ntx_sametx = dist_Ntx_sametx,
      dist_Ntx_difftx = dist_Ntx_difftx
    )
    
  }else{
    out <- data.frame(
      hhid = hh_id, 
      dist_popN = 0,
      dist_arm = NA,
      dist_Ntx = 0,
      dist_Ntx_sametx = 0,
      dist_Ntx_difftx = 0
    )
  }

  return(out)
}

dist_tx_target <- bind_rows(lapply(as.list(unique(xs_correct$hhid)), function(x)
  make_dist_covs(hh_id = x, min_dist = 0, max_dist = 500)))

dist_tx_spillover_1km <- bind_rows(lapply(as.list(unique(xs_correct$hhid)), function(x)
  make_dist_covs(hh_id = x, min_dist = 500, max_dist = 1000)))

dist_tx_spillover_2km <- bind_rows(lapply(as.list(unique(xs_correct$hhid)), function(x)
  make_dist_covs(hh_id = x, min_dist = 1000, max_dist = 2000)))

dist_tx_spillover_3km <- bind_rows(lapply(as.list(unique(xs_correct$hhid)), function(x)
  make_dist_covs(hh_id = x, min_dist = 2000, max_dist = 3000)))

nrow(xs_correct)
xs_cov <- xs_correct %>% 
  left_join(dist_tx_target, by = c("hhid")) %>% 
  dplyr::rename_at(vars(dist_popN:dist_Ntx_difftx), function(x) paste0(x, "_target"))

nrow(xs_cov)

xs_cov = xs_cov %>% 
  full_join(dist_tx_spillover_1km, by = c("hhid")) %>% 
  dplyr::rename_at(vars(dist_popN:dist_Ntx_difftx), function(x) paste0(x, "_1km"))

nrow(xs_cov)

xs_cov = xs_cov %>% 
  full_join(dist_tx_spillover_2km, by = c("hhid"))  %>% 
  dplyr::rename_at(vars(dist_popN:dist_Ntx_difftx), function(x) paste0(x, "_2km"))
  
nrow(xs_cov)

xs_cov = xs_cov %>% 
  full_join(dist_tx_spillover_3km, by = c("hhid")) %>% 
  dplyr::rename_at(vars(dist_popN:dist_Ntx_difftx), function(x) paste0(x, "_3km"))
  
nrow(xs_cov)

assert_that(nrow(xs_correct) == nrow(xs_cov))
  
# indicator for whether an individual was ever within distance d of an intervention 
xs_cov <- xs_cov %>%  mutate(
  # target zone 
  target_area = ifelse(dist_Ntx_target > 0, 1, 0),

  # set to missing if no treatment within 500m-1km
  sp1km = ifelse(dist_Ntx_target == 0 & dist_Ntx_1km > 0, 1, 0),
  
  # set to missing if no treatment within 1km-2km
  sp2km = ifelse(dist_Ntx_target == 0 & 
                   dist_Ntx_2km > 0, 1, 0),
  
  # set to missing if no treatment within 2km-3km
  sp3km = ifelse(dist_Ntx_target == 0 & dist_Ntx_3km > 0, 1, 0),

  # set to missing if no treatment within 500m-3km
  spall = ifelse(dist_Ntx_target == 0 & (dist_Ntx_1km > 0 |
                                           dist_Ntx_2km > 0 |
                                             dist_Ntx_3km > 0), 1, 0)
) 

# E covariates ---------------------------------------------------------------
df_ea <- ea_level_raw %>% 
  filter(ea_actual != "") %>% 
  mutate(ea = ea_actual,
         ### 2.2 malaria incidence in the season prior to the trial ------------------------------------------
         pre_incidence = ea_incidence_2016,
         ### 2.3 pre-season indoor residual spray coverage ------------------------------------------
         pre_spray_cover = mohss_spray_coverage_ea,
         ### 2.4 population size of EA ------------------------------------------
         pop_size_ea = ea_pop_update,
         ### 2.5 distance to nearest neighboring household  ------------------------------------------
         dist_hh = mean_nearest_hh_distance_meters,
         ### 2.6 distance to nearest healthcare facility  ------------------------------------------
         dist_hf = mean_hh_distance_to_hf,
         ### 2.7 median monthly rainfall   ------------------------------------------
         pre_rainfall = meanprecmedianealag2,
         ### 2.8 median enhanced vegetative index   ------------------------------------------
         pre_evi = medianevimedianealag2,
         ### 2.9 median elevation   ------------------------------------------
         ea_elevation = elevmedianea,
         ### 2.10 median daytime land surface temperature   ------------------------------------------
         surface_temp = medianlstmedianealag2_fixed) %>% 
  dplyr::select(c("ea",
                  "pre_incidence",
                  "pre_spray_cover",
                  "pop_size_ea",
                  "dist_hh",
                  "dist_hf",
                  "pre_rainfall",
                  "pre_evi",
                  "ea_elevation",
                  "surface_temp"))

xs_cov <- left_join(xs_cov, df_ea, by = c("ea"))

# reorder columns
xs_cov <- xs_cov %>% dplyr::select(
  eaid, ea, arm, hhid, xsid, 
  sp1km, sp2km, sp3km, spall,
  latitude, longitude, everything()
)

# recode covariates with missing and sparse categories
xs_recode <- xs_cov %>% mutate(
  traveled_recently = ifelse(traveled_recently=="Yes" & !is.na(traveled_recently), "Yes", "No/NA"),
  occupation = case_when(
    occupation == "Cattle herder" ~ "High risk",
    occupation == "Child, not student" ~ "Low risk/NA",
    occupation == "Construction worker" ~ "Low risk/NA",
    occupation == "Domestic worker" ~ "Low risk/NA",
    occupation == "Farm worker" ~ "High risk",
    occupation == "Nurse/Teacher/Professional" ~ "Low risk/NA",
    occupation == "Retired" ~ "Low risk/NA",
    occupation == "Piece work" ~ "Low risk/NA",
    occupation == "Reed/wood cutter" ~ "High risk",
    occupation == "Security guard" ~ "High risk",
    occupation == "Self-employed" ~ "Low risk/NA",
    occupation == "Small business" ~ "Low risk/NA",
    occupation == "Student" ~ "Low risk/NA",
    occupation == "Unemployed" ~ "Low risk/NA",
    occupation == "Other" ~ "Low risk/NA",
    occupation == "No response" ~ "Low risk/NA",
    is.na(occupation) ~ "Low risk/NA"
  )) %>% 
  mutate(sleep_net = as.character(sleep_net),
         slept_net_yest = as.character(slept_net_yest)) %>% 
  mutate(sleep_net = as.factor(case_when(
      is.na(sleep_net) ~ "Never/NA",
      sleep_net == "No response" ~ "Never/NA",
      sleep_net == "Never" ~ "Never/NA",
      sleep_net == "Always" ~ "Always",
      sleep_net == "Sometimes" ~ "Sometimes"
    ))
) %>% mutate(slept_net_yest = as.factor(case_when(
  slept_net_yest == "No" ~ "No/NA",
  slept_net_yest == "No response" ~ "No/NA",
  slept_net_yest == "Yes" ~ "Yes",
  is.na(slept_net_yest)  ~ "No/NA"
))) 

saveRDS(xs_recode, namibia_analysis_prev)





