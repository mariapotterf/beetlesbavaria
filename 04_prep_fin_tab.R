

#  Explore dependent vs predictors variables
# on raw beetle counts (per month) vs summarized values over year (sum, DOY of max increase, DOY of aggregation)


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
library(rgdal)
library(fasterize)
library(terra)
library(ggplot2)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# Stats
library('here')
library('mgcv')
library('gratia')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
library('gganimate')
library('cowplot')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('itsadug')
library(DHARMa)
library(MASS)
library(car)     # for VIF

# colors
library(RColorBrewer)


# for anomalies: 
reference_period = 1986:2010



# load cleaned data
load("outData/ips.Rdata")
load("outData/spatial.Rdata")

# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.aggreg         # get DOY of max increase
# - ips.year.avg       # sum beetles/year, 
# - df.daily           # avg/betles/year per trap
 
# Spatial data: 
# - xy_sf_fin  # XY as sf data, filtered to have only one globalid per trap

# Get SPEI and clim data: they are for whole year! check only veg season?
#df_spei       <- fread(paste(myPath, outTable, 'xy_spei.csv', sep = '/'))
df_prec       <- fread(paste(myPath, outTable, 'xy_precip.csv', sep = '/'))  # from DWD, 1 km res
df_temp       <- fread(paste(myPath, outTable, 'xy_temp.csv', sep = '/'))    # from DWD, 1 km res
#df_clim_ERA   <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = '/')) # need to filter soils volume content index!
df_clim_ERA_ID<- fread(paste(myPath, outTable, 'xy_clim_IDs.csv', sep = '/'), encoding = 'Latin-1') # get the ID to join falsto_name
df_anom       <- fread(paste(myPath, outTable, 'xy_anom.csv', sep = '/')) # from ERA, 10 km res, temp,. prec, soils water, spei,

# set correct encoding for german characters
Encoding(df_anom$falsto_name)        <-  "UTF-8"
Encoding(df_clim_ERA_ID$falsto_name) <-  "UTF-8"

# Get coniferous cover data: coniferous == 2! 
# spruce: share of spruce
df_tree   <- fread(paste(myPath, outTable, 'xy_treeComp.csv', sep = '/'))
#df_spruce <- fread(paste(myPath, outTable, 'xy_spruce.csv', sep = '/'))  # is this spruce of coniferous???

# Get geo data: elevation, slope, roughness...
xy_topo      <- vect(paste(myPath, outFolder, "xy_3035_topo.gpkg", sep = "/"), 
                  layer = 'xy_3035_topo') # read trap location


# Get climate data for traps: --------------------------------------------------
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_fin)[,"X"],
                    y = sf::st_coordinates(xy_sf_fin)[,"Y"],
                    falsto_name = xy_sf_fin$falsto_name,
                    globalid =    xy_sf_fin$globalid)

# get the final subset of globalids
trap_globid <- unique(xy_df$globalid)  # 158, selected


# convert to DF
df_topo <- as.data.frame(xy_topo)
df_topo <- df_topo %>% 
  filter(globalid %in% trap_globid) %>% 
  dplyr::select(c('globalid', 'elev','slope', 'aspect')) %>% 
  full_join(xy_df %>% dplyr::select(falsto_name, globalid))
#dplyr::select(globalid:roughness)


# Process predictors  data ------------------------------------

# select only coniferous: == 2, 0 is no forest (eg. covers the whole 500 m buffer)
df_conif <- df_tree %>% 
  filter(species == 2 ) %>% # 2 = conifers, can be spruce and pine
  filter(globalid %in% trap_globid) %>% 
  dplyr::select(-c(species_n, species)) %>%
  dplyr::rename(conif_prop = freq) %>% 
  full_join(xy_df %>% dplyr::select(falsto_name, globalid))

  
# whole year
df_predictors <- 
  df_anom %>% 
    filter(year %in% 2015:2021) %>%
    left_join(xy_df, by = c('falsto_name')) %>%
    left_join(df_topo, by = c("falsto_name", 'globalid')) %>%
    left_join(df_conif, by = c("falsto_name", 'globalid')) %>%
    dplyr::select(-c(OBJECTID, globalid, ID))
 

# Clean up dependent variables------------------------------

# add also my dependent variables: 
# [1] sum beetles 

# [2] Max Diff DOY
max.diff.doy <- max.diff.doy %>% 
  dplyr::select(c(falsto_name, year, doy, diff)) %>% 
  dplyr::rename(peak_doy = doy,
                peak_diff = diff)

# [3] Agg DOY !!!!: only 1080 rows!!! some traps dis not reached the values in time; need to increase threshold or somehow account for this!
ips.aggreg <- ips.aggreg %>% 
  dplyr::select(c(falsto_name, year, doy)) %>% 
  dplyr::rename(agg_doy = doy )
  
# [4] avg number of beetles per trap/catch
ips.year.avg <- ips.year.avg %>% 
  ungroup(.) %>% 
  dplyr::select(c(year, falsto_name, avg_beetles_trap))


# make a final table with all Ys, and all Xs to see a correlation matrix. reduce predictors to keep meaningful ones.----------
dat_fin <- df_predictors %>% 
  right_join(ips.year.sum) %>%
  #right_join(ips.year.avg) %>%
  right_join(max.diff.doy) %>%
  left_join(ips.aggreg) %>%
  dplyr::select(c(falsto_name,
                  year,
                 # prec,temp, 
                 # swvl, 
                 conif_prop,
                 elev,
                  x, y,
                  sm, vpd, 
                  sm_z, vpd_z,
                  tmp, tmp_z,
                  spei, spei_z,
                 # avg_beetles_trap,
                  sum_ips, 
                  peak_doy,
                  peak_diff, #, 
                  agg_doy
                 )) %>% 
  mutate(location = gsub('.{2}$', '', falsto_name),
         tmp      = tmp - 273.15) %>% # convert Kelvin to Celsius
  dplyr::rename(trap_name = falsto_name) %>% 
  mutate(location = factor(location),
         trap_name = factor(trap_name))



theme_set(theme_classic())

theme_update(aspect.ratio = 1, 
                panel.background = element_rect(fill = "white", colour = "black"))









p1<- dat_fin %>% 
  ggplot(aes(x = tmp,   # temp is in Kelvin ->convert to C by substrating 273.15, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview
             y = sum_ips,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = tmp,   # temp is in Kelvin ->convert to C by substrating 273.15, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview
                  y = sum_ips,
                  color= NULL),
              color = 'black')


p2<- 
  dat_fin %>% 
  ggplot(aes(x = tmp,
             y = peak_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = tmp,
                  y = peak_doy,
                  color = NULL),
              color = 'black')

p3<- dat_fin %>% 
  ggplot(aes(x = tmp,
             y = agg_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = tmp,
                  y = agg_doy,
                  color= NULL),
              color = 'black')


p4<- dat_fin %>% 
  ggplot(aes(x = tmp_z,
             y = sum_ips,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = tmp_z,
                  y = sum_ips,
                  color= NULL),
              color = 'black') 

p5<- dat_fin %>% 
  ggplot(aes(x = tmp_z,
             y = peak_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = tmp_z,
                  y = peak_doy,
                  color= NULL),
              color = 'black')
#facet_wrap(.~year)


p6<- dat_fin %>% 
  ggplot(aes(x = tmp_z,
             y = agg_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = tmp_z,
                  y = agg_doy,
                  color= NULL),
              color = 'black')




ggarrange(p1, p2, p3, 
          p4, p5, p6, common.legend = TRUE)





# plots Spei_z ------------------------------------------------------------

p1<- dat_fin %>% 
  ggplot(aes(x = spei,   # temp is in Kelvin ->convert to C by substrating 273.15, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview
             y = sum_ips,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = spei,   # temp is in Kelvin ->convert to C by substrating 273.15, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview
                  y = sum_ips,
                  color= NULL),
              color = 'black')


p2<- 
  dat_fin %>% 
  ggplot(aes(x = spei,
             y = peak_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = spei,
                  y = peak_doy,
                  color = NULL),
              color = 'black')

p3<- dat_fin %>% 
  ggplot(aes(x = spei,
             y = agg_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = spei,
                  y = agg_doy,
                  color= NULL),
              color = 'black')


p4<- dat_fin %>% 
  ggplot(aes(x = spei_z,
             y = sum_ips,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = spei_z,
                  y = sum_ips,
                  color= NULL),
              color = 'black') 

p5<- dat_fin %>% 
  ggplot(aes(x = spei_z,
             y = peak_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = spei_z,
                  y = peak_doy,
                  color= NULL),
              color = 'black')
#facet_wrap(.~year)


p6<- dat_fin %>% 
  ggplot(aes(x = spei_z,
             y = agg_doy,
             color= factor(year))) +
  geom_point(alpha = 0.5)  +
  geom_smooth(aes(x = spei_z,
                  y = agg_doy,
                  color= NULL),
              color = 'black')




ggarrange(p1, p2, p3, 
          p4, p5, p6, common.legend = TRUE)


# save final table --------------------------------------------------------


save(dat_fin, file="outData/final.Rdata")





# Save data ---------------------------------------------------------------

# save(ips_sum_preds,           # final df: predictors aggregated by year, sum beetle by year
#      df_predictors_year,      # df predictors per year, can be merged with different dependent variables
#      p_temp,                  # plot temp
#      p_prec,                  # plot prec
#      p_spei12,
#      p_spei12_bav, 
#      p_spei_freq, 
#      vif_values_fin,          # VIF values for final set of predictors
#      file="outData/predict.Rdata") 













