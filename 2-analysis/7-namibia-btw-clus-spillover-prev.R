################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Test for between-cluster spillovers 
################################################

rm(list=ls())
library(rgeos)
library(readstata13)
library(MASS)
library(lmtest)

source(paste0(here::here(), "/0-config.R"))

# load locality shape file 
shp <- rgdal::readOGR(dsn =namibia_shapefile_dsn,
                      layer = namibia_shapefile_layer)

prev <- readRDS(namibia_analysis_prev)

prev_ea <- prev %>% group_by(eaid) %>% 
  summarise(prev = mean(qPCRposneg, na.rm=T))

# identify neighboring locality to each locality ------------------------------------
neighbors = gTouches(shp, byid=TRUE)

# get list of all EAs
eas = unique(shp$EA_No)

# get list of EAs in the trial 
trial_eas = eas = unique(shp$EA[shp$Arm2!=0])

# subset matrix of neighboring localities to those in 
# the trial 
ind_trial_eas = apply(as.matrix(trial_eas), 1, function(x) 
  which(eas == x))
neighbors_trial = neighbors[ind_trial_eas, ind_trial_eas]
assert_that(setequal(dim(neighbors_trial), c(56,56)))

# identify neighbors of each trial EA 
# create list where each element is an EA and
# contains a vector of that EA's neighbors
ea_neighbor = list()

for(i in 1:length(trial_eas)){
  x = neighbors_trial[which(trial_eas == trial_eas[i]),]
  ea_neighbor[[i]] = trial_eas[which(x)]
}
names(ea_neighbor) = trial_eas

 
# transform data to prep for test of contamination ------------------------------------
untx_ea = shp$EA[is.na(shp$Eas_Stud_1)]
RO_ea = shp$EA[shp$Eas_Stud_1=="RO" & !is.na(shp$Eas_Stud_1)]
RV_ea = shp$EA[shp$Eas_Stud_1=="RV" & !is.na(shp$Eas_Stud_1)]
TO_ea = shp$EA[shp$Eas_Stud_1=="TO" & !is.na(shp$Eas_Stud_1)]
TV_ea = shp$EA[shp$Eas_Stud_1=="TV" & !is.na(shp$Eas_Stud_1)]

RO_list = ea_neighbor[names(ea_neighbor) %in% RO_ea]
RV_list = ea_neighbor[names(ea_neighbor) %in% RV_ea]
TO_list = ea_neighbor[names(ea_neighbor) %in% TO_ea]
TV_list = ea_neighbor[names(ea_neighbor) %in% TV_ea]

# classify each neighbor as untreated, same, or different tx
classify_number <- function(d, tx){
  x = matrix(NA, length(d), 1)
  y = matrix(NA, length(d), 1)
  
  for(i in 1:length(d)){
    
    x[i,1] = case_when(
      d[i] %in% RO_ea ~ "RO",
      d[i] %in% RV_ea ~ "RV",
      d[i] %in% TO_ea ~ "TO",
      d[i] %in% TV_ea ~ "TV"
    )
    
    y[i,1] = case_when(
      d[i] %in% RO_ea & tx == "RO" ~ "Concordant",
      d[i] %in% RO_ea & tx != "RO" ~ "Discordant",
      d[i] %in% RV_ea & tx == "RV" ~ "Concordant",
      d[i] %in% RV_ea & tx != "RV" ~ "Discordant",
      d[i] %in% TO_ea & tx == "TO" ~ "Concordant",
      d[i] %in% TO_ea & tx != "TO" ~ "Discordant",
      d[i] %in% TV_ea & tx == "TV" ~ "Concordant",
      d[i] %in% TV_ea & tx != "TV" ~ "Discordant"
    )
  }
  
  neighb_prev <- prev %>% filter(eaid %in% as.numeric(d)) %>% 
    summarise(mean_prev = mean(qPCRposneg, na.rm=T))

  # summarize neighbor characteristics
  out = data.frame(
    n_neighbors = length(d),
    n_disc_neighbors = length(y[y=="Discordant"]),
    n_disc_neighb_RO = sum(x[y=="Discordant"]=="RO"),
    n_disc_neighb_RV = sum(x[y=="Discordant"]=="RV"),
    n_disc_neighb_TO = sum(x[y=="Discordant"]=="TO"),
    n_disc_neighb_TV = sum(x[y=="Discordant"]=="TV"),
    neighb_prev = neighb_prev$mean_prev
  )

  return(out)
}

RO_neighbors = lapply(RO_list, function(x) classify_number(x, tx = "RO")) %>% 
  bind_rows() %>%
  mutate(EA = names(RO_list))
RV_neighbors = lapply(RV_list, function(x) classify_number(x, tx = "RV")) %>% 
  bind_rows() %>% 
  mutate(EA = names(RV_list))
TO_neighbors = lapply(TO_list, function(x) classify_number(x, tx = "TO")) %>% 
  bind_rows() %>% 
  mutate(EA = names(TO_list))
TV_neighbors = lapply(TV_list, function(x) classify_number(x, tx = "TV")) %>% 
  bind_rows() %>% 
  mutate(EA = names(TV_list))

# one row for each EA with summary covariates
# about neighboring EAs
all_neighb_data = bind_rows(
  RO_neighbors, RV_neighbors, TO_neighbors, TV_neighbors)

# merge in other EA-level variables 
neighb_data = merge(all_neighb_data, 
          shp@data %>% dplyr::select(EA, Eas_StudyA, Eas_Stud_1), 
          by = "EA") %>% 
  rename(arm = Eas_Stud_1) %>% 
  mutate(EA = as.numeric(EA))


# merge in incidence data ------------------------------------
m <- merge(neighb_data, prev_ea,
           by.x = "EA", by.y = "eaid", 
           all.x = TRUE, all.y = FALSE)
assert_that(nrow(m) == nrow(neighb_data))

m <- m %>% 
  # filter(ea_actual!="RO49") %>% 
  mutate(disc_neigh_cat = case_when(
    n_disc_neighbors == 1 ~ "1", 
    n_disc_neighbors == 2 ~ "2", 
    n_disc_neighbors == 3 ~ "3", 
    n_disc_neighbors >= 4 ~ "4+"
  )) %>% 
  mutate(disc_neigh_cat = as.factor(disc_neigh_cat))


# plot incidence by neighbor incidence ------------------------------------
ggplot(m , aes(x = neighb_prev, y= prev)) + 
  geom_point() + 
  geom_smooth() 

ggplot(m , aes(x = neighb_prev, y= prev)) + 
  geom_point(aes(col=arm)) + 
  facet_wrap(~arm) +
  geom_smooth() 

m = m %>% mutate(armlabel = case_when(
  arm == "RO" ~ "RACD only",
  arm == "TO" ~ "rfMDA only",
  arm == "TV" ~ "rfMDA + RAVC",
  arm == "RV" ~ "RACD + RAVC"
)) %>% mutate(armlabel = factor(armlabel,
                                levels = c("RACD only", "rfMDA only", "RACD + RAVC", "rfMDA + RAVC")))

mypalette <- c("#1D6996","#0F8554","#EDAD08","#94346E")

plot <- ggplot(m , aes(x = neighb_prev, y= prev)) + 
  geom_point(aes(col=armlabel)) + 
  facet_grid(~armlabel) +
  xlab("Prevalence in adjacent clusters") + 
  ylab("Prevalence") +
  scale_color_manual(values= mypalette) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(plot, filename = paste0(figure_path, "plot-between-clus-spillover-prev.png"),
       width = 8, height =3)



# model prevalence by neighbor prevalence ------------------------------------
m = m %>% mutate(
  tpe = ifelse(arm=="TV" | arm=="TO", 1, 0), 
  ravc = ifelse(arm=="RV" | arm=="TV", 1, 0)
)

# model without neighboring EA prevalence
fit = glm(prev ~  tpe + ravc + tpe*ravc, data = m)
summary(fit)

# model with neighboring EA prevalence (continuous)
fit_neighb = glm(prev ~  tpe + ravc + tpe*ravc + neighb_prev,
                data = m)
summary(fit_neighb)

# model with neighboring EA prevalence (categorical)
fit_disc_neighb = glm(prev ~  tpe + ravc + tpe*ravc + disc_neigh_cat ,
                       data = m)
summary(fit_disc_neighb)

# likelihood ratio test comparing models
lrtest(fit, fit_neighb)
lrtest(fit, fit_disc_neighb)

# model prevalence by neighbor prevalence ------------------------------------
# drop outlier 
msub <- m %>% filter(EA!="10599007" & EA!="10599024")

# model without neighboring EA prevalence
fit_sub = glm(prev ~  tpe + ravc + tpe*ravc, data = msub)
summary(fit_sub)

# model with neighboring EA prevalence (continuous)
fit_neighb_sub = glm(prev ~  tpe + ravc + tpe*ravc + neighb_prev,
                    data = msub)
summary(fit_neighb_sub)

# model with neighboring EA prevalence (categorical)
fit_disc_neighb_sub = glm(prev ~  tpe + ravc + tpe*ravc + disc_neigh_cat ,
                      data = msub)
summary(fit_disc_neighb_sub)

# likelihood ratio test comparing models
lrtest(fit_sub, fit_neighb_sub)
lrtest(fit_sub, fit_disc_neighb_sub)

