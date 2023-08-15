################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial map of cohorts
################################################
rm(list=ls())
library(sp)
library(raster)
library(rgdal)
library(leaflet)
library(sf) # for working with shape files
library(htmlwidgets) # for saving html
library(webshot2) # for saving portion of html

source(paste0(here::here(), "/0-config.R"))

#------------------------------------------------
# load polygons
#------------------------------------------------
ea_shape <- rgdal::readOGR(dsn = namibia_shapefile_dsn,
                           layer = namibia_shapefile_layer)

namibia_shape <- rgdal::readOGR(dsn = paste0(box_path, "K01-data/Shapefiles/nam_adm_nsa_ocha_20200109_shp/nam_admbnda_adm2_nsa_ocha_20200109.shp"),
                                layer = "nam_admbnda_adm2_nsa_ocha_20200109")


#------------------------------------------------
# load and process data
#------------------------------------------------
data = readRDS(paste0(namibia_clean_path, "namibia_index_clean.RDS")) %>% 
  filter(!is.na(longitude_index_imp)) %>% 
  filter(first_case_yn>0) %>% 
  mutate(month = month(visit_date_index))  %>% 
  mutate(
    tx_human = case_when(
      intarm=="TV" | intarm=="TO" ~ "rfMDA",
      intarm=="RV" | intarm=="RO" ~ "RACD"),
    tx_mosq = case_when(
      intarm=="TV" | intarm=="RV" ~ "RAVC",
      intarm=="RO" | intarm=="TO" ~ "No RAVC"),
    tx_both = case_when(
      intarm=="TV" ~ "rfMDA + RAVC",
      intarm=="RO" ~ "RACD only",
      TRUE ~ "")
  ) %>% 
  mutate(intarm = case_when(
    intarm == "RO" ~ "RACD",
    intarm == "RV" ~ "RACD+RAVC",
    intarm == "TO" ~ "rfMDA",
    intarm == "TV" ~ "rfMDA+RAVC"
  ))

# identify months with most index cases (1-6)
table(data$month)

data_6m = data %>% filter(month<=6) 
data_5w = data %>% filter(visit_date_index >= as.Date("2017-04-25") & 
                            visit_date_index <= as.Date("2017-05-30"))

indexcase_spdf <- SpatialPointsDataFrame(
  coords = data[,c("longitude_index_imp", "latitude_index_imp")],
  data = data[,c("intarm","tx_human", "tx_mosq","tx_both","id")],
  proj4string = CRS("+init=epsg:4326"))

indexcase_6m_spdf <- SpatialPointsDataFrame(
  coords = data_6m[,c("longitude_index_imp", "latitude_index_imp")],
  data = data_6m[,c("intarm","tx_human", "tx_mosq","tx_both","id")],
  proj4string = CRS("+init=epsg:4326"))

indexcase_5w_spdf <- SpatialPointsDataFrame(
  coords = data_5w[,c("longitude_index_imp", "latitude_index_imp")],
  data = data_5w[,c("intarm","tx_human", "tx_mosq","tx_both","id")],
  proj4string = CRS("+init=epsg:4326"))

indexcase_6m_spdf_both <- SpatialPointsDataFrame(
  coords = data_6m[data_6m$intarm %in% c("rfMDA+RAVC", "RACD"),
                   c("longitude_index_imp", "latitude_index_imp")],
  data = data_6m[data_6m$intarm %in% c("rfMDA+RAVC", "RACD")
                 ,c("intarm","tx_human", "tx_mosq","tx_both","id")],
  proj4string = CRS("+init=epsg:4326"))

#------------------------------------------------
# map style
#------------------------------------------------
ea_shape@data$intrial = ifelse(ea_shape@data$Arm2>0, "Yes", "No")
studyarea_palette <- colorFactor(c("#8d99ae","#FFFFFF"),
                                 ea_shape$intrial)

# arm_palette <- colorFactor(c("#f6e8c3","#d8b365","#c7eae5","#5ab4ac"),
#                            indexcase_spdf$intarm)

arm_palette_human <- colorFactor(c("#b5c8b8","#2887a1"),
                                 indexcase_5w_spdf$tx_human)

arm_palette_mosq <- colorFactor(c("#d6bd8d","#A16928"),
                           indexcase_spdf$tx_mosq)

arm_palette_both <- colorFactor(c("#edeac2","#798234"),
                                indexcase_6m_spdf_both$tx_both)



#------------------------------------------------
# map 
#------------------------------------------------
map_human <-
  leaflet(data = ea_shape, options = leafletOptions(zoomControl = FALSE,
          attributionControl = FALSE)) %>% 
  # set zoom center point and level
  setView(lng = 23.37, lat = -17.86, zoom = 11) %>% 
  
  # set background white rectangle
  addRectangles(lng1 = 22.5, lng2 = 24.5, 
                lat1 = -17.4, lat2 = -18.8,
                opacity = 0,
                fillOpacity = 1,
                fillColor = "white") %>% 
  
  # add ea polygons
  addPolygons(data = ea_shape,
              weight = 1.5,
              color = "black",
              fillOpacity = 0
  ) %>%
  
  # circle for spillover zones - 5 weeks
  addCircles(data = indexcase_5w_spdf,
             weight = 1.5,
             radius = 1000,
             color = "black",
             fillOpacity = 0.8,
             fillColor = ~arm_palette_human(indexcase_5w_spdf$tx_human)) %>%
  
  # circle for target areas - 5 weeks
  addCircles(data = indexcase_5w_spdf,
             weight = 1.5,
             radius = 500,
             fillColor = "#8a8781",
             color = "black") 

saveWidget(map_human, file=paste0(figure_path, "map_human.html"))
webshot2::webshot(paste0(figure_path, "map_human.html"), 
                  file=paste0(figure_path, "map_human.png"), 
                  cliprect="viewport",
                  vwidth = 550,
                  vheight=700)

map_mosq <-
  leaflet(ea_shape, options = leafletOptions(zoomControl = FALSE,
                                             attributionControl = FALSE)) %>% 
  # set zoom center point and level
  setView(lng = 23.37, lat = -17.86, zoom = 11) %>% 
  
  # set background white rectangle
  addRectangles(lng1 = 22.5, lng2 = 24.5, 
                lat1 = -17.4, lat2 = -18.8,
                opacity = 0,
                fillOpacity = 1,
                fillColor = "white") %>% 
  
  # add ea polygons
  addPolygons(data = ea_shape,
              weight = 1.5,
              color = "black",
              fillOpacity = 0
  ) %>%

  # circle for spillover zones - 6 months
  addCircles(data = indexcase_6m_spdf,
           weight = 1.5,
           radius = 1000,
           color = "black",
           fillOpacity = 0.8,
  fillColor = ~arm_palette_mosq(indexcase_6m_spdf$tx_mosq)) %>%
  
  # circle for target areas - 6 months
  addCircles(data = indexcase_6m_spdf,
             weight = 1.5,
             radius = 500,
             fillColor = "#8a8781",
             color = "black")

saveWidget(map_mosq, file=paste0(figure_path, "map_mosq.html"))
webshot2::webshot(paste0(figure_path, "map_mosq.html"), 
                  file=paste0(figure_path, "map_mosq.png"), 
                  cliprect="viewport",
                  vwidth = 550,
                  vheight=700)


map_both <-
  leaflet(ea_shape, options = leafletOptions(zoomControl = FALSE,
                                             attributionControl = FALSE)) %>% 
  # set zoom center point and level
  setView(lng = 23.37, lat = -17.86, zoom = 11) %>% 
  
  # set background white rectangle
  addRectangles(lng1 = 22.5, lng2 = 24.5, 
                lat1 = -17.4, lat2 = -18.8,
                opacity = 0,
                fillOpacity = 1,
                fillColor = "white") %>% 
  
  # add ea polygons
  addPolygons(data = ea_shape,
              weight = 1.5,
              color = "black",
              fillOpacity = 0
  ) %>%
  
  # circle for spillover zones - 6 months
  addCircles(data = indexcase_6m_spdf_both,
             weight = 1.5,
             radius = 1000,
             color = "black",
             fillOpacity = 0.8,
             fillColor = ~arm_palette_both(indexcase_6m_spdf_both$tx_both)) %>%
  
  # circle for target areas - 6 months
  addCircles(data = indexcase_6m_spdf_both,
             weight = 1.5,
             radius = 500,
             fillColor = "#8a8781",
             color = "black")

saveWidget(map_both, file=paste0(figure_path, "map_both.html"))
webshot2::webshot(paste0(figure_path, "map_both.html"), 
                  file=paste0(figure_path, "map_both.png"), 
                  cliprect="viewport",
                  vwidth = 550,
                  vheight=700)

map_all <- leaflet(ea_shape, options = leafletOptions(zoomControl = FALSE,
                                                      attributionControl = FALSE)) %>% 
  # set zoom center point and level
  setView(lng = 23.83, lat = -17.9, zoom = 9.5) %>% 
  
  # set background white rectangle
  addRectangles(lng1 = 22.5, lng2 = 24.8, 
                lat1 = -17.4, lat2 = -18.8,
                opacity = 0,
                fillOpacity = 1,
                fillColor = "white") %>% 
  
  # add namibia polygons
  addPolygons(data = namibia_shape,
              weight = 1.5, 
              color = "black",
              fillOpacity = 0
              ) %>% 
  
  # add ea polygons
  addPolygons(data = ea_shape,
              weight = 1.5,
              color = "black",
              fillOpacity = 0
  ) %>%

  # circle for spillover zones
  addCircles(data = indexcase_spdf,
             weight = 1.5,
             radius = 1000,
             color = "black",
             fillOpacity = 0.1) %>%

  # circle for target areas
  addCircles(data = indexcase_spdf,
              weight = 1.5,
              radius = 500,
             fillColor = "#8a8781",
             fillOpacity = 0.8,
             color = "black")

saveWidget(map_all, file=paste0(figure_path, "map_all.html"))
webshot2::webshot(paste0(figure_path, "map_all.html"), 
                  file=paste0(figure_path, "map_all.png"), 
                  cliprect="viewport",
                  vwidth = 1000,
                  vheight=700)
