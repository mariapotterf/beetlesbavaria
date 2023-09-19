

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
# - ips.year.sum       # sum beetles/year per trap

# Vars
n_neighbors = 10      # number of nearest neighbors

# Spatial data: 
# should have 158 regions (trap pairs)
sort(unique(xy_sf_fin$falsto_name))


# get coordinates from sf object
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_fin)[,"X"],
                    y = sf::st_coordinates(xy_sf_fin)[,"Y"],
                    falsto_name = xy_sf_fin$falsto_name)


# Get sums of IPS beetle per year/trap: April 31 to October 30
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)) %>% 
  left_join(xy_df, by = c("falsto_name")) # df_xy3035


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
                          k = n_neighbors),  # nearest neighbor
               sym = TRUE)
  
  #  # Calculate LISA
  lisa_res<-localmoran(ips_sum_sub$sum_beetle, nb2listw(nb))
  #
  # add LISA to points:
  ips_sum_sub$Morans_I <-lisa_res[,1] # get the first column: Ii - local moran  stats
  ips_sum_sub$clust <-attributes(lisa_res)$quadr$mean  # get classified data
  #
  # convert to sf for plotting
  ips_sum_sub_sf <- st_as_sf(ips_sum_sub)
  
  return(ips_sum_sub_sf)
}

# Run LISA on all years separately
lisa_out <- lapply(years, get_lisa )

# Merge all in one sf
lisa_merged <- dplyr::bind_rows(lisa_out)


#st_crs(lisa_merged) <- crs(bav_sf)
#lisa_merged_proj <- st_crs(crs(bav_sf))
#  st_transform(lisa_merged, projection(bav_sf))
# Create map fr each one of them and save as saparate object


p_lisa_sub <- ggplot() +
  geom_sf(data = filter(lisa_merged, 
                        clust %in% c("Low-Low", "High-High")),
          aes(color = clust)) +
  scale_color_manual(breaks = c("Low-Low", "High-Low", "Low-High", "High-High"),
                     values=c("blue", "grey90", "grey90", 'red')) +
  facet_wrap(year~.) +
  theme_void() +
  ggtitle('LISA: Moran I')



p_lisa_all <- ggplot() +
  geom_sf(data = lisa_merged, #filter(lisa_merged, 
                        #clust %in% c("Low-Low", "High-High")),
          aes(color = clust)) +
  scale_color_manual(breaks = c("Low-Low", "High-Low", "Low-High", "High-High"),
                     values=c("blue", "grey90", "grey95", 'red')) +
  facet_wrap(year~.) +
  theme_void() +
  ggtitle('LISA: Moran I')


# convert to df 
lisa_merged_df <- as.data.frame(lisa_merged)


p_lisa_freq <- 
  lisa_merged_df %>%
  group_by(clust, year) %>% 
    dplyr::summarize(freq = n()) %>%
    ggplot(aes(x = year,
               y = freq,
               fill = clust)) +
    geom_col()
    






  
# LISA single year -----------------------------------------------------------------
# Local variation analysis (LISA):
  # can be done per year
  

ips_sum_15 <-ips_sum %>% 
  filter(year == 2015)


ips_sum_20 <-ips_sum %>% 
  filter(year == 2020)


# Create a spatial points data frame
coordinates(ips_sum_15) <- ~ x + y
coordinates(ips_sum_20) <- ~ x + y

# Create queen contiguity neighbors object
nb_15 <- knn2nb(knearneigh(coordinates(ips_sum_15),
                           k = n_neighbors),  # nearest neighbor
                sym = TRUE)
nb_20 <- knn2nb(knearneigh(coordinates(ips_sum_20),
                           k = n_neighbors),  # nearest neighbor
                sym = TRUE)



lisa_res_15<-localmoran(ips_sum_15$sum_beetle, nb2listw(nb_15))
lisa_res_20<-localmoran(ips_sum_20$sum_beetle, nb2listw(nb_20))


# add Lisa to points:
ips_sum_15$Morans_I <-lisa_res_15[,1] # get the first column: Ii - local moran  stats
ips_sum_20$Morans_I <-lisa_res_20[,1] # get the first column: Ii - local moran  stats



# convert to sf object:
ips_sum_15_sf <- st_as_sf(ips_sum_15)
ips_sum_20_sf <- st_as_sf(ips_sum_20)



# get classifiued values
ips_sum_15$cl_mean <-attributes(lisa_res_15)$quadr$mean
ips_sum_15_sf <- st_as_sf(ips_sum_15)














# Global Moran =================================================================

# get more distant neighbors:
nb_g_15 <- knn2nb(knearneigh(coordinates(ips_sum_15),
                             k = n_neighbors),  # nearest neighbor
                  sym = TRUE)
nb_g_20 <- knn2nb(knearneigh(coordinates(ips_sum_20),
                             k = n_neighbors),  # nearest neighbor
                  sym = TRUE)


# Calculate the global Moran's I
global_moran_15 <- moran.test(ips_sum_15$sum_beetle, nb2listw(nb_g_15) )
global_moran_20 <- moran.test(ips_sum_20$sum_beetle, nb2listw(nb_g_20) )

# !!!
# for all years:


get_global_moran <- function(i, ...) {
  
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  #
  #  # Global Moran's per year:
  #  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  #
  #  # Create queen contiguity neighbors object
  nb <- knn2nb(knearneigh(coordinates(ips_sum_sub),
                          k = n_neighbors),  # nearest neighbor
               sym = TRUE)
  
  #  # Calculate Global Moran's I
  global_moran_res <- moran.test(ips_sum_sub$sum_beetle, nb2listw(nb) )
  #
  # export as df
  df <- data.frame(year = i,
                   stat = global_moran_res$statistic,
                   p_val = global_moran_res$p.value )
  
  # export table with results
  return(df)
}

# Run Global mora on all years separately
global_moran_out <- lapply(years, get_global_moran )

# Merge all in one sf
glob_merged <- dplyr::bind_rows(global_moran_out)




# Get variograms: -----------------------------------------------------------------

library(gstat)


# how to make a function, with output raw data and model?? for each year??
# get only raw date, the models need to be done visually: adjust sill, nugget and range
get_variogram <- function(i, ...) {
  
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  #
  #  # Global Moran's per year:
  #  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
 
  # export as df
  variogram_raw <- variogram(log10(sum_beetle)  ~ 1, 
                             cutoff=cutoff,
                             data = ips_sum_sub)  # all directions are assumed equal
  # add year
  variogram_raw$year = i
  
   # export table with results
  return(variogram_raw)
}


# Run on all years separately
variogram_out <- lapply(years, get_variogram )

# Merge all in one sf
var_merged <- dplyr::bind_rows(variogram_out)


# plot
ggplot(var_merged,aes(x = dist/1000,
                      y = gamma,
                      color = year)) +
  geom_point() #+
 # geom_line() 







# Load dataset sliced by year 

# 2015
ips_sum_15 <-ips_sum %>% 
  filter(year == 2015)


# Create a spatial points data frame
coordinates(ips_sum_15) <- ~ x + y


#data(ips_sum_15)

# Create a "SpatialPointsDataFrame" object
#coordinates(your_data) <- c("X", "Y")

# Define the variogram model
windows()

cutoff = 300000

variogram_raw_directions <- variogram(log10(sum_beetle)  ~ 1, 
                             data = ips_sum_15,  
                             alpha = c(0, 45, + 90, 135)) # check different directions (angles)

variogram_raw <- variogram(log10(sum_beetle)  ~ 1, 
                             cutoff=cutoff,
                             data = ips_sum_15)  # all directions are assumed equal


# Plot the semivariogram
plot(variogram_raw)

# plot variogram cloud:
variogram_cloud <- variogram(log10(sum_beetle)~1, data=ips_sum_15, 
                             cutoff=cutoff,width = 5000, cloud=TRUE)

plot(variogram_cloud)

# Fit a variogram model (spherical model in this example)
# show different functions: 
show.vgms(par.strip.text=list(cex=0.75))

m1 <-  fit.variogram(variogram_raw, vgm(1, 'Sph', 1), fit.kappa = TRUE)

m2 <-  fit.variogram(variogram_raw, 
                     vgm(psill = 0.2, model="Cir", range = 200000, nugget = 0.11))


# Plot the fitted semivariogram
plot(variogram_raw, # raw data 
     m2,            # fitted model 
     main = "Fitted Semivariogram", cutoff = cutoff)

# Print the model parameters
print(m2)

# 2020 

ips_sum_20 <-ips_sum %>% 
  filter(year == 2020)


# Create a spatial points data frame
coordinates(ips_sum_20) <- ~ x + y


#data(ips_sum_15)

# Create a "SpatialPointsDataFrame" object
#coordinates(your_data) <- c("X", "Y")

# Define the variogram model
windows()

cutoff = 300000

variogram_raw_directions <- variogram(log10(sum_beetle)  ~ 1, 
                                      data = ips_sum_20,  
                                      alpha = c(0, 45, + 90, 135)) # check different directions (angles)

variogram_raw <- variogram(log10(sum_beetle)  ~ 1, 
                           cutoff=cutoff,
                           data = ips_sum_20) # check different directions (angles)


# Plot the semivariogram
plot(variogram_raw)

# plot variogram cloud:
variogram_cloud <- variogram(log10(sum_beetle)~1, data=ips_sum_15, 
                             cutoff=cutoff,width = 5000, cloud=TRUE)

plot(variogram_cloud)

# Fit a variogram model (spherical model in this example)
# show different functions: 

m2 <-  fit.variogram(variogram_raw, 
                     vgm(psill = 0.2, model="Lin", range = 200000, nugget = 0.11))


# Plot the fitted semivariogram
windows()
plot(variogram_raw, m2, main = "Fitted Semivariogram", cutoff = cutoff)


# Humblt: fitting variograms: https://gsp.humboldt.edu/olm/R/04_01_Variograms.html


# Example ----------------------
# eye-ball modell fitting



library(geoR)
library(sf)

library(sp)        # for meuse dataset

library("gstat")   # geostatistics
library("mapview") # map plot
library("sf")      # spatial vector data
library("stars")   # spatial-temporal data
library("terra")   # raster data handling 
library("ggplot2") # plotting
mapviewOptions(fgb = FALSE)

meuse <- ips_sum_15 #read_sf('data/meuse.gpkg')
v.eye <- eyefit(variog(as.geodata(meuse["sum_beetle"]), max.dist = 10000))
ve.fit <- as.vgm.variomodel(v.eye[[1]])


# from WSL

library(gstat)
# https://gsp.humboldt.edu/olm/R/04_01_Variograms.html

library(sp)
data(meuse)
# no trend:
coordinates(meuse) = ~x+y
plot(variogram(log(zinc)~1, meuse))



# Save outputs ------------------------------------------------------------



save(p_lisa_sub, 
     p_lisa_all,
     glob_merged,
     file = "outData/lisa.Rdata")




