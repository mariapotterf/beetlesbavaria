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
xy        <- vect(paste(myPath, outFolder, "xy_fin_3035.gpkg", sep = "/"), 
                layer = 'xy_fin_3035') # read trap location
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


# Read .nc data as a raster in terra - way faster!  
dat_ras <- terra::rast(paste(myPath, inFolder,  "ERA_Bav.nc", sep = "/"))

# extract all values to the xy coordinates:
dat_ext_df <- terra::extract(dat_ras, xy_latlng)

# get which ID equals which falsto_name
xy_latlng$ID = 1:nrow(xy_latlng)

xy_names <- data.frame(ID = xy_latlng$ID, 
                       falsto_name = xy_latlng$falsto_name)


# Check the variables of interest:
var_names <- varnames(dat_ras)    # 9 variables
ras_time  <- time(dat_ras)        # 441
n_layer   <- nlyr(dat_ras)        # 441 


# terra works well, but require correct naming of variables:

sort(unique(xy_latlng$OBJECTID)) == sort(unique(dat_ext_df$ID))

length(unique(xy_latlng$falsto_name))
# naming: need to be make manually:
# years:  2015-2021: 7
# months: april-October: 7
# hour:   at 12:00
# variables: 9

# Create a datatable for each site (ID), variable, and time
df_melt <- data.table::melt(dat_ext_df, 
                            id.vars = c('ID'))

# Add time  to df and split in months:
df <- df_melt %>% 
  arrange(ID, variable) %>%
  mutate(time = rep(ras_time, nrow(dat_ext_df))) %>% 
  separate(variable, 
           c("var", "time_num"), "_") %>% 
  dplyr::mutate(year  = lubridate::year(time), 
                month = lubridate::month(time), 
                day   = lubridate::day(time),
                doy   =  lubridate::yday(time) + 1)  # as POXIT data has January 1st at 0


# combine soil moisture data: sum the 4 layers:
soil_vars <- c("swvl1", "swvl2", "swvl3", "swvl4")


# Split df in two tales: having all variables  besides soil water content,
# and only soil water content
df_vars <- df %>% 
  filter(!(var %in% soil_vars))
  

# For soil water content: need to get sums as Cornelius?
df_soil <- df %>% 
  filter(var %in% soil_vars ) %>%
  group_by(ID, year, month, day) %>% 
  mutate(#sum_swv = sum(value),
         value = sum(value)) 


# remove the variable and create a new one:
df_soil <- df_soil %>% 
  #dplyr::select(-c(var)) %>% 
  mutate(var = c('swv')) %>% 
  dplyr::select("ID",   # correctly order the columns
                "var",
                "time_num",
                "value",
                "time", 
                "year",
                "month" ,
                "day",
                "doy"  )
  


# Merge the df_soil and df_vars onto single file
df_out <- rbind(df_soil, df_vars)



# Export the final table:

data.table::fwrite(df_out, 
                   paste(myPath, outTable, 'xy_clim.csv', sep = "/"))



data.table::fwrite(xy_names, 
                   paste(myPath, outTable, 'xy_clim_IDs.csv', sep = "/"))





# Check data attribution: -----------------------------------------------------










# check temperatures: if seems correct
df %>% 
  filter(var %in% soil_vars)  %>%  # temperature
  ggplot(aes(x = time,
             y = value,
             group = as.factor(year)
             )) +
  geom_point() +
  geom_smooth() +
  facet_grid(var ~ .)


# Add horizontal ine for means:
X.mean <- df %>%
  filter(var %in% soil_vars)  %>%  
#  group_by(grp) %>%
  summarize(y = mean(value))

df %>% 
  filter(var %in% soil_vars)  %>%  # temperature
  ggplot(aes(x = as.factor(year),
             y = value,
             fill = var)) +
  geom_boxplot(outlier.colour = NA ) +
  geom_hline(data = X.mean, aes(yintercept = y), col = 'red', linetype = 'dashed') 


  

# Check per individual months:
df %>% 
  filter(var %in% soil_vars)  %>%  # temperature
  ggplot(aes(x = as.factor(month),
             y = value,
             fill = var)) +
 # stat_summary(fun.data = "mean_cl_boot",  size = 0.5) +
  geom_boxplot(outlier.colour = NA ) +
  #geom_hline(data = X.mean, aes(yintercept = y), col = 'red', linetype = 'dashed')  +
  facet_grid(.~factor(year))



df %>% 
  filter(var %in% soil_vars)  %>%  # temperature
  ggplot(aes(x = as.factor(month),
             y = value,
             fill = var)) +
  stat_summary(fun.data = "mean_cl_boot",  size = 0.5, aes(group = var,col = var)) +
 # geom_boxplot(outlier.colour = NA ) +
  #geom_hline(data = X.mean, aes(yintercept = y), col = 'red', linetype = 'dashed')  +
  facet_grid(.~factor(year)) +
  stat_summary(fun = median,
               geom = "line",
               aes(group = var,col = var)) + 
  theme(legend.position = 'bottom')
