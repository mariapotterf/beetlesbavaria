
#######################################
#
#     Analyze the beetle data         #
#
#######################################

#-	Do beetle data correlate with Cornelius data? By hexagons??
#    - 
#    - 	Is there less beetles if there is more deciduous trees?
# 
#-	Effects of drought: 
#    -	Is there more beetles from 2018? Eg.g. drought will trigger beetle population
#    - compare the counts between years   
# 
# -	Effect of sampling design: do I get the same results given the random selection of two traps?


##########################################################
#
# 1. Do beetle counts data correlate with RS data??         #
#
##########################################################

# Process:

# get input data
# get a grid
# get info of mortality over years from RS per cell
# aggregate traps into cells, standardize by number of traps per year/cell
# merge databases by cell_id
# does it correlate over Bavaria? over years?


# Input data:

# bavaria shp
# disturbance raster
# xy traps
# beetle data
# raster species composition: deciduous/coniferous


rm(list = ls())


# Read libs  --------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
library(rgdal)
library(tidyverse)
library(lubridate)
library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)



# Read my paths -----------------------------------------------------------
source('myPaths.R')


# Get sf data: bavaria, grid and XY traps
bav.shp <- st_read(paste(myPath, inFolder, "outline_bavaria.gpkg", sep = "/"), 
                   layer = 'outline_bavaria') # read watershed
grid    <- st_read(paste(myPath, inFolder, 'grid_12.shp', sep = '/'))
xy      <- st_read(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                   layer = 'xy_3035') # read watershed

# Get rasters 
disturbance <- raster(paste(myPath, inFolder,  "disturbance_map_bavaria.tif", sep = "/"))
forest      <- raster(paste(myPath, inFolder,  "forestcover_bav.tif", sep = "/"))
forest_type <- raster(paste(myPath, outFolder, "bav_fortype_ext30_int2u_LZW.tif", sep = "/"))


# Get beetle counts
dat      <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))


# Maybe run analyses on teh adjacency of the traps??
# try it out
# create the buffers~ 500m around the trap location;
# extract the disturbance data from the disturbance maps
# correlate teh values


# Process vector (sf) data ------------------------------------------------------

# Simplify bavaria geometry, keep only geometry
bav.simple <- st_simplify(bav.shp, preserveTopology = FALSE, dTolerance = 1000)

# remove unnecessary colums:
bav.simple <- bav.simple %>% 
  dplyr::select('geom')


# Project data to the same system
#st_crs(ext) <- st_crs(grid)
xy         <- st_transform(xy, st_crs(grid)) 
bav.simple <- st_transform(bav.simple, st_crs(grid))  


#plot(grid, add = T)
#plot(bav.simple, add = T)
#plot(xy['OBJECTID'], add = T)

# Create buffers around XY point: r = 500 m to link the estimate RS damage
xy_500 = xy %>% 
  #dplyr::filter(grepl("_1",falsto_name)) %>%   # Keep only the first trap from the pair (to limit the sample size)
  st_buffer(500)   # 500 m

plot(xy_500['OBJECTID'])



# Process raster data -----------------------------------------------------

# Mask disturbances only to the coniferous forest
# limit the forest damage to the conifeorous forests, to corresponds to beetle data

# convert first to terra format:
forest_terra  <- rast(forest_type)
disturb_terra <- rast(disturbance)

# # simply substitute the values by NA, keep only ceniferous = 2
forest_mask <-subst(forest_terra, 1, NA) # change deciduous to NA
forest_mask <-subst(forest_terra, 0, NA)# change background to NA


# Keep only disturbances > 2014, replace other values by NA
disturb_terra14 <- disturb_terra
values(disturb_terra14)[values(disturb_terra14) < 2013 ] <- NA



# Make sure that the extend of both rasters fit:
# Resample forest raster to match/snap/align the disturbance raster
# Buffer raster is already resampled
forest_mask_resample <- terra::resample(forest_mask,  # raster to be resampled 
                                     disturb_terra14,      # Master raster
                                     method = 'near')


# Extract by mask; creates a list of datasets
#forest_ex <- terra::extract(forest_mask_resample, buff, list = F)
#dist_ex   <- terra::extract(dist,     buff, list = F)
#buff_ex   <- terra::extract(buff_ras, buff, list = F)



# Mask the disturbance data by the coniferous forest extent
disturb13 <- mask(disturb_terra14, forest_mask_resample)


windows()
plot(disturb13)


