
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
source('process_input_data.R')





# Aggregate raster disturbances to grid -------------------------------------------------

# Create grid index for each pixel and export it as a long vector of values
grid_sel     <- st_intersection(st_as_sf(grid), st_as_sf(ext))
grid_sel_ras <- fasterize(grid_sel, disturbance14, field = "FID") # name the new rasters as polygon ids
grid_values  <- values(grid_sel_ras)     # resolution is still 30 m, grid value for each pixel; this is a vector of values




# Check the length of each raster:  144738748;
# all values have the same length
length(values(forest))
length(grid_values)
length(disturbance14)
length(values(forest_type))


# Create dataframes containing grid id and disturbaance vs forest data
forest_df <-
  data.frame(gridindex = grid_values,
             forest    = values(forest)) %>%
  group_by(gridindex) %>%
  summarize(forest_ha = sum(forest == 1, na.rm = TRUE) * 0.09,
            land_ha = n() * 0.09) %>%
  ungroup(.)


# Get in the same way the disturbance data aggregated by grid
dist_df <-
  data.frame(gridindex = grid_values,
             dist = values(disturbance14)) %>%
  na.omit(.) %>%
  # filter(!is.na(disturbance)) %>%
  group_by(gridindex) %>%
  summarize(disturbance_ha = n() * 0.09) %>%
  ungroup(.)

# merge forest and disturbance data
dist_df <- forest_df %>% 
  left_join(dist_df, by = "gridindex") %>% 
  rename(OBJECTID = gridindex)


# Export files:
#write_csv(dist_df, paste(myPath, outTable, "dist_grid_15.csv", sep = "/"))





#####################################################
#              Beetle data                          #
#####################################################

# Spatial join between XY data, beetle counts by hexagon, by year 

# get input data
# overlay with the grid values to get the new attribute - interssect?
# Calculate the beetle counts by number of active traps
# 

# Add grid FID to the XY data 
xy <- st_join(xy, grid)

# get rid of geometry
xy_df <- xy %>%  st_drop_geometry()


# Join beetle count data with XY and grid attributes
dat_df <- dat %>% 
  full_join(xy_df, by = 'globalid')


# How many traps do I have per hexagon?
# they vary by year;
# n = number of traps by hexagon, year and species 
trap_n <- 
  dat_df %>% 
  filter(year > 2014 & representativ == 'Ja' ) %>% # & art == 'Buchdrucker' 
  group_by( #monsto_name,  # trap location name
    OBJECTID,     # hexagon ID
    year,
    art) %>% 
  tally() %>%
  rename(trap_n = n)
#print(n = 40)


# Let's see if the counts are right: 
# filter just specific values:  Aldersbach          224  2017    23
dat_df %>% 
  filter(monsto_name == 'Aldersbach' & year == 2017 & art == 'Buchdrucker' & OBJECTID == 218 & representativ == 'Ja')


dat_df %>% 
  distinct(einheit)  # Stuck = number/count, NA?


# get number of beetles per hexagon by years
beetle_n <- dat_df %>% 
  group_by(art, year, OBJECTID) %>% 
  summarize(beetle_sum = sum(fangmenge, na.rm = T))

# Add the beetle count data to the number of traps & standardize the beetle counts by a trap number
out_df <- beetle_n %>% 
  left_join(trap_n, by = c('OBJECTID', 'art', 'year')) %>% 
  mutate(beetle_by_trap = beetle_sum/trap_n) %>% 
  left_join(dist_df) #%>% 
#mutate(lag_dist_sum = lag(dist_sum))

# Add lagged disturbance data by one year:
# need to be correctly arranged;
out_df <- out_df %>% 
  filter(year > 2014) %>% 
  group_by(art, OBJECTID) %>% # , OBJECTID 
  arrange(year, .by_group = TRUE) %>% 
  mutate(lag_dist_sum = lag(disturbance_ha))




# Check the correlation between the beetle numbers and damage from RS:
# just using stat_smooth
out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = disturbance_ha,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_grid(art ~ year, scales = 'free')


out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Kupferstecher') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = disturbance_ha,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_grid(art ~ year, scales = 'free')



# Check out the effects of the lagged values
windows()
out_df %>% 
  filter((art == 'Kupferstecher' & beetle_by_trap < 30000) | 
           (art == 'Buchdrucker' & beetle_by_trap < 6000 )) %>%  #filter(year > 2014 & beetle_by_trap< 30000) %>% 
  #filter(art == 'Buchdrucker'   & beetle_by_trap < 6000) #%>% 
  # filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = lag_dist_sum,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_wrap(art ~ ., scales = 'free')


windows()
out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Kupferstecher') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = lag_dist_sum,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") #+
#facet_grid(art ~ year, scales = 'free')






# Check the overall counts and damages over years? also, the high rates of beetles can predict mortality in 
# following years?
year_sum_df <- 
  out_df %>%
  filter(year > 2014) %>% 
  group_by(art, year) %>% 
  summarize(beetle_sum = sum(beetle_by_trap , na.rm = T),
            dist_sum   = sum(disturbance_ha, na.rm = T)) %>% 
  mutate(lag_dist_sum = lag(dist_sum))


# Plot the yearly values nad lagged values
windows()

# Plots yearly values
p1<- year_sum_df %>% 
  ggplot(aes(x = beetle_sum,
             y = dist_sum,
             color = art)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(.~art, scales = 'free')


# Plot lagged values
p2 <- year_sum_df %>% 
  ggplot(aes(x = beetle_sum,
             y = lag_dist_sum,
             color = art)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(.~art, scales = 'free')


windows()
ggarrange(p1, p2, ncol = 1)
