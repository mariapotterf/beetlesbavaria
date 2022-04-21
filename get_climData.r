# 
#
# Get climatic data
#
#

rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')


# Read libs  --------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(raster)
library(tidyverse)
library(lubridate)
library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
#library(easyclimate) # does not work on R 4.1
# WorldClim - ends in 2018!!


# Get spatial data for each trap
xy      <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                layer = 'xy_3035') # read watershed
xy_latlng <- terra::project(xy, 
                          "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# Starting from netCFD: 
# https://cran.r-project.org/web/packages/futureheatwaves/vignettes/starting_from_netcdf.html
# Packages for climate variables:
library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

# Load ERA5 soil moisture data --------------------------------------------

# Data was downloaded as netCDF from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview

# Read soil water content, 4 levels:
vars <- c('u10', # ,  10 m component of wind
          'src', 'tp',  "swvl1", "swvl2", "swvl3", "swvl4")  # , 

out <- vector("list", length(vars))

# Conevrt vect from terra to sf object
xy_latlng_sf <- sf::st_as_sf(xy_latlng)


# Crop dataset for the study region, for each of teh rasters:
for (i in 1:length(vars)) {
  
  print(vars[i])
  
  
  # Try with Raster
  dat_ras <- raster::brick(paste(myPath, inFolder,  "ERA_data_soilMoisture.nc", sep = "/"), varname = vars[i])
  dat_ras <- raster::stack(dat_ras)
 
  dat_ras <- extract(dat_ras, xy_latlng_sf)
  
  out[[k]] <- dat_ras %>%
      as.matrix()
  
}



# Tr with terra:

for (i in 1:length(vars)) {
  
dat_ras <- rast(paste(myPath, inFolder,  "ERA_data_soilMoisture.nc", sep = "/"))

#  str(dat_ras)

dat_ras2 <- extract(dat_ras, xy_latlng)
#  dat_ras <- raster::mask(dat_ras, studyregion_latlng)

# out[[k]] <- dat_ras %>%
#    as.matrix()

}


# terra works well, but require correct naming of variables:
# has ID = OBJECTID from teh vect

# naming: need to be make manually:
# years: 2015-2021: 7
# months: april-October: 7
# hour: at 12
# variables: 7

# Check time from terra
time(dat_ras)

# Total number of names: 7 (years)*7(months)*7(variables)
# 343

# How to link the variables names with the dates???





r<-raster::brick(paste(myPath, inFolder,  "ERA_data_soilMoisture.nc", sep = "/"), var = vars[1])




# Nc-open = just opens the coccention to the file, does not read the file
# needs to be closed after it is not needed by 'nc_close()'
climate_output <- nc_open(paste(myPath, inFolder,  "ERA_data_soilMoisture.nc", sep = "/"))


# Get the index of the long-lat values:
lon <- ncvar_get(climate_output, varid = "longitude")
lat <- ncvar_get(climate_output, varid = "latitude")


# What is the time interval?
climate_output$dim$time$units
# "hours since 1900-01-01 00:00:00.0"


# Check what is my calendar:
climate_output$dim$time$calendar
# "gregorian"










r <- getData("worldclim", var="bio", res=0.5, lon=-72, lat=44)


# -----------------------------------
# Dumy exmaple

set.seed(802)
long <- runif(10, -72.85, -72.78)
lat <- runif(10, 44.5, 44.6)
# a vector of ID numbers for these coordinates
ID <- 1:10
# bind the long and lat into a dataframe
coords <- data.frame(long, lat)

# Check wehere aare the points located?
library(leaflet)
library(maps)
# visualizing the ten coordinates we generated:
leaflet(data=coords) %>% 
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>% 
  addCircleMarkers(~long, ~lat, label=as.character(ID))


# downloading the bioclimatic variables from worldclim at a resolution of 30 seconds (.5 minutes)
r_bio  <- raster::getData("worldclim", var="bio", res=0.5, lon=11.142, lat=48.89)
r_tmin <- raster::getData("worldclim", var="tmin", res=0.5, lon=11.142, lat=48.89)
r_tmax <- raster::getData("worldclim", var="tmax", res=0.5, lon=11.142, lat=48.89)
r_prec <- raster::getData("worldclim", var="prec", res=0.5, lon=11.142, lat=48.89)
# lets also get the elevational data associated with the climate data
alt <- raster::getData("worldclim", var="alt", res=.5, lon=11.142, lat=48.89) # middle bavaria


# how many layers?
nlayers(r) 

windows()
plot(r)

# Download the file for bavaria: no need!
# First get the coordinates:
# lat-long
# NW: 50.43 9.14
# NE: 50.36 13.68
# SE  47.47 14.09
# SW: 47.57 7.56

# seems that data end in 2018? !

# Check the observation data from Meteorological stations:
# https://cran.r-project.org/web/packages/climate/vignettes/getstarted.html


# Try climate package -------------------------------------------------------------



library(climate)
DE = stations_ogimet(country = "Germany", add_map = TRUE)

# Download the the annual summary of air temperature:
df = meteo_imgw(interval = "monthly", rank = "synop", year = 2014:2021, station = 10007) 



# Read data from ERA



