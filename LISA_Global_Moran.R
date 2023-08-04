

# Spatial patters analysis: 

# LISA - local indicator of spatial association
# Global Moran's I
# semivariance = analyze the spatial dependence or variability of a dataset;  quantifies the average dissimilarity between pairs of observations at different distances or lags. 
# variogram


# from ChatGPT:

# Split data by years, have XY trap locations
# 1. Local Variation Analysis: investigate  concept of spatial autocorrelation. 
#    Compute (LISA) using the localmoran() function from the spdep package.
#    - will find the hotspots (High-High, Low-Low)/coldspots, clusters  
#    - can visualizeLISA using map 

# 2. Regional Variation Analysis: overall spatial autocorrelation in  data. 
#  - Calculate the global Moran's I statistic using the spdep::moran.test() 
#  - This test will determine whether the beetle counts are spatially autocorrelated at a broader regional level.

# 3. Scale-dependent Analysis: 
#  - compute semivariograms or variograms to assess the spatial dependence of beetle counts at different lag distances. 
#  - The gstat::variogram() 
# Plot the semivariograms to observe how the spatial dependence changes across different scales.


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

library(sf)
library(plyr)
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
library(ggplot2)
library(ggpubr)

# Spatstat
library(spdep)

# load cleaned data
load("outData/ips.Rdata")
load("outData/spatial.Rdata")

# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.year           # sum beetles/year, avg/betles/year per trap

# Spatial data: 
# - xy_sf  # XY as sf data
xy_terra <- vect(xy_sf)

# rename the coordinates into xy
xy_df <- data.frame(x = geom(xy_terra)[,'x'],
                        y = geom(xy_terra)[,'y'],
                        globalid = xy_terra$globalid)


# Get sums of IPS beetle per year/trap: April 31 to October 30
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name, globalid) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)) %>% 
  left_join(xy_df, by = c("globalid")) # df_xy3035


nrow(ips_sum)

# Run LISA for each year separately
years <- 2015:2021

get_lisa <- function(i, ...) {
 
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  #
  #  # Local variation analysis (LISA) per year:
  #  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  #
  #  # Create queen contiguity neighbors object
  nb <- knn2nb(knearneigh(coordinates(ips_sum_sub),
                          k = 1),  # nearest neighbor
               sym = TRUE)
  
  #  # Calculate LISA
  lisa_res<-localmoran(ips_sum_sub$sum_beetle, nb2listw(nb))
  #
  # add LISA to points:
  ips_sum_sub$Morans_I <-lisa_res[,1] # get the first column: Ii - local moran  stats
  #
  # convert to sf for plotting
  ips_sum_sub_sf <- st_as_sf(ips_sum_sub)
  
  return(ips_sum_sub_sf)
}

# Run LISA on all years separately
lisa_out <- lapply(years, get_lisa )

# Merge all in one sf
lisa_merged <- dplyr::bind_rows(lisa_out)

# Create map fr each one of them and save as saparate object


ggplot() +
  geom_sf(data = lisa_merged, 
          aes(#x = x, 
            #y = y, 
            color = Morans_I)) +
  scale_color_gradient(low = "white", 
                       high = "red", 
                       name = "Moran's I") +
  facet_wrap(year~.) +
  theme_void() +
  ggtitle('LISA: Moran I')












ips_sum_15 <-ips_sum %>% 
  filter(year == 2015)


ips_sum_20 <-ips_sum %>% 
  filter(year == 2020)

# LISA -----------------------------------------------------------------
# Local variation analysis (LISA):
# can be done per year


# Get LISA: 

# Create a spatial points data frame
coordinates(ips_sum_15) <- ~ x + y
coordinates(ips_sum_20) <- ~ x + y

# Create queen contiguity neighbors object
nb_15 <- knn2nb(knearneigh(coordinates(ips_sum_15),
                           k = 1),  # nearest neighbor
                sym = TRUE)
nb_20 <- knn2nb(knearneigh(coordinates(ips_sum_20),
                           k = 1),  # nearest neighbor
                sym = TRUE)



lisa_res_15<-localmoran(ips_sum_15$sum_beetle, nb2listw(nb_15))
lisa_res_20<-localmoran(ips_sum_20$sum_beetle, nb2listw(nb_20))


# add Lisa to points:
ips_sum_15$Morans_I <-lisa_res_15[,1] # get the first column: Ii - local moran  stats
ips_sum_20$Morans_I <-lisa_res_20[,1] # get the first column: Ii - local moran  stats



# convert to sf object:
ips_sum_15_sf <- st_as_sf(ips_sum_15)
ips_sum_20_sf <- st_as_sf(ips_sum_20)





ggplot() +
  geom_sf(data = ips_sum_15_sf, 
          aes(#x = x, 
            #y = y, 
            color = Morans_I)) +
  scale_color_gradient(low = "white", 
                       high = "red", 
                       name = "Moran's I") +
  theme_void() +
  ggtitle('LISA: Moran I')



# Global Moran =================================================================

# get more distant neighbors:
nb_g_15 <- knn2nb(knearneigh(coordinates(ips_sum_15),
                             k = 20),  # nearest neighbor
                  sym = TRUE)
nb_g_20 <- knn2nb(knearneigh(coordinates(ips_sum_20),
                             k = 20),  # nearest neighbor
                  sym = TRUE)


# Calculate the global Moran's I
global_moran_15 <- moran.test(ips_sum_15$sum_beetle, nb2listw(nb_g_15) )
global_moran_20 <- moran.test(ips_sum_20$sum_beetle, nb2listw(nb_g_20) )

# !!!
