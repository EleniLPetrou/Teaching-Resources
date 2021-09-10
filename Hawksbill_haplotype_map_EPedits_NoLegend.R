# Load libraries
library(maps)
library(mapdata)
library(tidyverse)
library(readxl)
library(maptools)
library(RColorBrewer)
library(classInt)
library(scatterpie)
library(cowplot)
library(ggplot2)
library(rgdal)
library(sf)



# Specify working directory where input file is found
#setwd("C:/Users/lhaus/Dropbox (MERLAB)/Projects/Turtles/ZA Stranding/Ei")
setwd("C:/Users/Eleni/Downloads")

####################################################################
# Read in data and haplotype metadata
EiMap <- read_excel("ArantesSI2_EP.xls", sheet = "Eimap")
haplo_metadata <- read_excel("ArantesSI2_EP.xls", sheet = "hawkbill_hap_metadata")

# save some metadata information as vectors (for plotting)
Ei_cols <- haplo_metadata$Ei_cols
Ei_haps <- haplo_metadata$Ei_haps

####################################################################
# Read in the shape file containing management units (save both .shp and .shx files in the data directory)

Ei_RMU <- st_read("Ei_RMU_20101013.shp")

####################################################################
# Create a database of coastline coordinates for the world, for plotting
world <- map_data('worldHires')
head(world)


# plot the management units


test_map <- ggplot() +
  geom_map(data = world, map = world, aes(map_id = region), fill = "gray", color = "gray") +
  geom_sf(data = Ei_RMU, colour = "blue", fill = "lightblue1") +
  theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-120, 180), ylim = c(-40, 40)) #specify map limits

test_map



# add haplotype pie charts to the test map
final_plot <- test_map + 
  geom_scatterpie(data = EiMap,
                  aes(x = long, y = lat, group = region, r = lnsum),
                  legend_name = "haplotype",
                  sorted_by_radius = FALSE, 
                  cols = Ei_haps,
                  color = "black", 
                  alpha = 0.8) +
  scale_fill_manual(values = Ei_cols) + #custom color scale
  theme(legend.position = "none") #remove legend
  

final_plot
