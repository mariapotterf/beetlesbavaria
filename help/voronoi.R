# Get voronoi tesselation

# for each trap, for each location (two traps)
# get the sums of beetle counts per polygon/location (sum beetle counts on location level)

# Read my paths -----------------------------------------------------------
source('myPaths.R')

library(terra)
library(rnaturalearth)
library(dplyr)
library(data.table)

# load cleaned data
load("outData/ips.Rdata")
load("outData/spatial.Rdata")


xy_year <-   vect(paste(myPath, "outSpatial/xy_fin_years_3035.gpkg", sep = '/'))
# Download the shapefile for Bavaria
bavaria <- vect('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/outSpatial/bavaria.gpkg')
bav_proj <- project(bavaria, crs(xy_year))

# Compute Voronoi polygons
voronoi <- voronoi(xy_year)

crs(voronoi) == crs(bav_proj)


# Clip Voronoi polygons by Bavaria
clipped_voronoi <- crop(voronoi, bav_proj)

# Plot clipped Voronoi polygons
plot(clipped_voronoi, main="Clipped Voronoi Tessellation")
plot(bav_proj, add=TRUE)
plot(xy_year, add=TRUE)


# Plot Voronoi polygons
plot(voronoi, main="Voronoi Tessellation")
points(terra_data, col="red", pch=19)


# Read data
disturb_year <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/germany114/disturbance_year_1986-2020_germany.tif')
disturb_type <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/fire_wind_barkbeetle_germany.tif')
forest_cover <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/germany114/forestcover_germany.tif')

# tree species classification: from 2017/2018 - some disturbed areas previously can be already removed???
spruce_cover <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/tree_species_map/spruce_bavaria_resampl30.tif')


# ID,species
# 2,Birch
# 3,Beech
# 4,Douglas fir
# 5,Oak
# 6,Alder
# 8,Spruce
# 9,Pine
# 10,Larch
# 14,Fir
# 16,ODH
# 17,ODL
crs(disturb_type) <- crs(disturb_year)
crs(spruce_cover) <- crs(disturb_year)


# crop rasters on Bavaria
disturb_year_bav <- crop(disturb_year, bav_proj)
disturb_type_bav <- crop(disturb_type, bav_proj)
forest_cover_bav <- crop(forest_cover, bav_proj)




#  Process  ---------------------------------------------------------------

xy <- xy_year[xy_year$year == 2015, ]
voronoi <- voronoi(xy)

# Clip Voronoi polygons by Bavaria
clipped_voronoi <- crop(voronoi, bav_proj)
clipped_voronoi$ID = 1:nrow(clipped_voronoi)


# split by years
# make vronoi
# rasterize to raster to fit disturbance data
# get unique ID
# calculate the damage, forest and agent per polygon
# use hexagon approach, as teh area is fully covered

# 
# Create a grid index for each pixel
grid_sel     <- terra::intersect(clipped_voronoi, bav_proj)
grid_sel_ras <- rasterize(clipped_voronoi, disturb_year_bav, field = "ID")  # Name the new raster as polygon ids


dim(grid_sel_ras)
dim(disturb_year_bav)
dim(disturb_type_bav)
dim(forest_cover_bav)

plot(grid_sel_ras)
plot(xy, add = T)


# Get in the same way the disturbance data by year aggregated by grid
dist_df <-
  data.frame(gridindex =   as.vector(values(grid_sel_ras)),
             dist_year = as.vector(values(disturb_year_bav))) %>%
  na.omit(.) %>%
  group_by(gridindex, dist_year) %>%   # 'dist' is the year of disturbanec here
  dplyr::summarize(disturbance_ha = n() * 0.09) %>%
  ungroup(.) %>% 
  rename(year = dist_year)


