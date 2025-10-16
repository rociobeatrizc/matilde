library(sf)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(terra)
library(geodata)
library(osmdata)
library(osmextract)


# upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# plot region
plot(aoi_abruzzo)

# bounding box 
abruzzo_bb <- st_bbox(aoi_abruzzo)

# from OSM select type of roads: primary, secondary, tertiary (paths)
ht_secondary <- "secondary"

# download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_secondary, ]
