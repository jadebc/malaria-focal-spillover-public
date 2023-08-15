################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Assess whether incidence and prevalence are 
# associated with distance to health facility 
################################################

rm(list=ls())
library(geosphere)

source(paste0(here::here(), "/0-config.R"))

# load data ------------------

# endline prevalence survey
xs_clean <- readRDS(namibia_analysis_prev)

# incidence datasets
df_long <- readRDS(namibia_df_long_path)

# health facility data 
hf <- read.csv(paste0(box_shared_path, "Datasets/Final trial data for collaborators/11_health_facilities_location.csv"))
hf_ll <- hf %>% dplyr::select(Longitude, Latitude)
hf_ids <-  hf %>% pull(FacID)

# calculate distance to closest health care facility ------------------
hh_id = xs_clean$hhid[1]

dist_to_hf <- function(data, id_name, id_value){
  
  ll = data %>% 
    dplyr::filter(!!rlang::sym(id_name) == id_value) %>% 
    dplyr::select(longitude, latitude) %>% 
    distinct()
  
  dists <- apply(hf_ll, 1, function(x) distm(x, as.matrix(ll)))
  
  # save in km
  min_dist <- dists[which.min(dists)]/1000
  
  out <- data.frame(min_dist_hf = min_dist) %>% 
    mutate({{id_name}} := id_value) 
  
  return(out)
}


# association between prevalence and distance to closest health care facility ------------------

## obtain distances for prevalence data ---
list_hhids <- as.list(unique(xs_clean$hhid))

result_prev <- lapply(list_hhids, function(x) 
  dist_to_hf(data = xs_clean, id_name = "hhid", id_value = x)) %>% bind_rows()

xs_dist <- left_join(xs_clean, result_prev, by = "hhid")
saveRDS(xs_dist, file = paste0(namibia_clean_path, "xs_dist_to_hf.RDS"))

## fit model ---
fit <- glm(qPCRposneg ~ min_dist_hf, data = xs_dist, family = "binomial")
summary(fit)

ggplot(xs_dist, aes(y=qPCRposneg, x = min_dist_hf))+
  geom_smooth()

# association between incidence and distance to closest health care facility ------------------
## obtain distances for incidence data ---
## fit model ---
fit <- glm(indexcase ~ min_dist_hf, data = df_long, family = "binomial")
summary(fit)

ggplot(df_long, aes(y=indexcase, x = min_dist_hf))+
  stat_smooth(method = 'gam',
              formula = y~s(x, bs= "cs", k=2)) 



