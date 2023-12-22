

##########################################################
#
# 1. Do beetle counts data correlate with RS data??         #
#
##########################################################

# Process:

# get input disturbance data: drought vs windthrow
# get a XY of the final trap set: one trap can have several locations! 
#            needs to make a new location file, with 
#            locations of the shifting traps: maybe then rerun LISA et GLobal Moran's?
# get one XY file per year!

# make a buffers around XY data
# get info of mortality over years from RS per cell
# aggregate traps into cells, standardize by number of traps per year/cell
# merge databases by cell_id
# does it correlate over Bavaria? over years?


# Input data:

# bavaria shp
# disturbance raster
# xy traps
# beetle count data
# raster species composition: deciduous/coniferous

rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')


# Read libs  --------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
library(rgdal)
library(lubridate)
library(fasterize)
library(ggpubr)
library(terra)


# load cleaned data
load("outData/ips_counts.Rdata")
load("outData/spatial.Rdata")

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


#rast(paste(myPath, outFolder,  "disturbance14_conif.tif", sep = "/"))
crs(disturb_year)
crs(disturb_type)
# disturbance BIOTIC: 
# 1 = wind & beetles
# 2 = fire
# 3 = harvest

crs(disturb_type) <- crs(disturb_year)
crs(spruce_cover) <- crs(disturb_year)

# get trap data
xy_sf_all  # all XY traps
xy_sf_fin  # selected traps with consistent monitoring , one trap have only one XY per all years (I neglected thet they moved over years)


# Clean & expand trap cooordinates ------------------------------------------
# to do: get the shifting XY per trap to work with the disturbance map data!
trap_names = unique(xy_sf_fin$falsto_name)


# get study period
all_years    = 2006:2021
study_period = 2015:2021

# Expand coordinates table to have a geometry for each trap and year
xy_sf_expand <- 
  xy_sf_all %>% 
  filter(falsto_name %in% trap_names ) %>% 
    dplyr::mutate(x = sf::st_coordinates(.)[,1],
                  y = sf::st_coordinates(.)[,2]) %>% 
  dplyr::select(von_dat, bis_dat, falsto_name, globalid, x, y) %>% #, geometry
  mutate(from_year = lubridate::year(dmy(von_dat)))  %>% # get the year from the starting date
  group_by(falsto_name) %>% 
  complete(from_year = all_years)  %>% 
  arrange(falsto_name, from_year) %>% 
 # df %>%
   tidyr::fill(globalid, .direction = 'down') %>%
    tidyr::fill(x, .direction = 'down') %>%
    tidyr::fill(y, .direction = 'down')%>%
   # dplyr::rename(from_year = year) #%>% 
  filter(from_year %in% study_period) %>% # filter only years 2015:2021
    st_as_sf(coords = c("x", "y"), 
             crs = 3035) %>% 
  mutate(year = from_year) %>% 
  dplyr::select(-c(from_year, von_dat, bis_dat))


xy_sf_expand$id <- 1:nrow(xy_sf_expand)

# export XZ file with shifting traps
xy_sf_expand %>%  
  st_write(paste(myPath, "outSpatial/xy_fin_years_3035.gpkg", sep = '/'), append=FALSE)


# Alternatively, save only specific objects to the .Rdata file
save(list=c("xy_sf_expand"), file="outData/spatial.Rdata")




#  ---------------------------------------------------------------------
# How much mortality happened in each buffer? 
# ----------------------------------------------------------------------

# Check the projection systems: 
# Get the correct projection: need to be a projected system! 
# better to project the vector than the raster as the vector maintains the geometry and precision for each vertex: not possible for each pixel


# Project the vector (XY) into raster data
# Check if the CRS is the same
#crs(disturbance) == crs(xy)

# convert to terra object; keep XY naminf for conveniece
xy <- vect(xy_sf_expand)

xy$id <- 1:nrow(xy)

# change the projection of the XY data
xy <- terra::project(xy, crs(disturb_year))


# ---------------------------------------------------------------
#                       Process
# ---------------------------------------------------------------

# Make buffers
# Intersect buffers with forest and disturbance maps:
# forest values: calculate prop of forest type from all forest
#                         prop of forest type from whole buffer 
# export as values from rasters to have a dataframe
# for each trap: 

# create a list to look over
xy_ls <- terra::split(xy, c('id')) #,# "OBJECTID" # 

buff_width = 1000  # 1 km

# make a function to export df from each buffer: keep years, disturbance type and forest:
# filter by spruce cover (30 m)
extract_rst_val <- function(xy, ...) {
  
  # xy<- xy_ls[[951]]
  # xy<- xy_ls[[952]]
   id <- xy$id
   #print(id)
   
   # get buffer
   buff <- buffer(xy, buff_width)
   #print('buff')
   
   # disturbance year
   # crop data and then mask them to have a circle
   dist_year_crop <- crop(disturb_year, buff)
   dist_year_mask <- mask(dist_year_crop, buff)
   #print('2')
   
   
   # disturbance type: harvest, fire, beetle
   dist_type_crop <- crop(disturb_type, buff)
   dist_type_mask <- mask(dist_type_crop, buff)
   #print('3')
   
   # tree class raster
   # crop data and then mask them to have a circle
   spruce_cover_crop <- crop(spruce_cover, buff)
   spruce_cover_mask <- mask(spruce_cover_crop, buff)
   
   # Get vector of values
   val_year <- as.vector(values(dist_year_mask))
   val_type <- as.vector(values(dist_type_mask))
   val_spruce <- as.vector(values(spruce_cover_mask))
   
   # merge data and drop NA values
   df <- data.frame(year = val_year,
                    type = val_type,
                    spruce = val_spruce,
                    id = id) 
   df <- df %>% 
     drop_na()

   
   return(df)
  
  }


# make a function to export df from each buffer: keep years, disturbance type and forest;
# this forest is from 1986, and include all species
extract_forest <- function(xy, ...) {
 # xy1<- xy_ls[[8]]
  id <- xy$id
  
  # get buffer
  buff <- buffer(xy, buff_width)
  
  # disturbance year
  # crop data and then mask them to have a circle
  forest_crop <- terra::crop(spruce_cover, buff)
  forest_mask <- terra::mask(forest_crop, buff)
  
  # Get vector of values
  val <- as.vector(values(forest_mask))
  
  # count number of cells:
  df <- as.data.frame(table(val))
  
  # add name
  df$id <- id
  
  return(df)
  
}

options(terra.pardegree = 1)

# run over all locations
# extract disturbances
out_dist_ls <- lapply(xy_ls, extract_rst_val) #extract_rst_val(xy_ls[[10]])

# merge partial tables into one df
dist_df      <-do.call("rbind", out_dist_ls)


# pre-process table, calculate individual pixels into categories
dist_df_out <- 
  dist_df %>% 
  #mutate(year = as.numeric(as.character(year))) %>% # convert factor to numeric
    group_by(year, type, id) %>%
      dplyr::select(-spruce) %>% 
      dplyr::summarize(Freq = n()) %>% # get count of rows
     spread(type, Freq) 
  
# rplelace NA by 0
dist_df_out[is.na(dist_df_out)] <- 0
  
# output needed: for every year, I need how many cells was disturbed by what agent
# coding
# 1 = wind & beetles
# 2 = fire
# 3 = harvest

head(dist_df_out)

names(dist_df_out) <- c("year" , "id",'wind_beetle', "harvest")



# extract spruce cover per buffer  ------------------------------------------------
out_forest_ls <- lapply(xy_ls, extract_forest ) #extract_rst_val(xy_ls[[10]])

# merge partial tables into one df
forest_df      <-do.call("rbind", out_forest_ls)

spruce_df <- forest_df %>% 
  dplyr::rename(spruce_1986 =Freq) %>% 
  dplyr::select(-val)


forest_df %>% 
  filter(id == 8)
  



#  Merge forest and yearly disturbances -----------------------------------
# calculate remaining spruce forest, first need to expand to all years
dist_df_exp <- 
  dist_df_out %>% 
  group_by(id) %>% 
  complete(year = 1986:2020,  # expand to all years
           fill = list(id = id))  #%>%

 

# Merge disturance data with spruce per buffer, to know how much spruce I have left
df_merge <- 
  dist_df_exp %>% 
  mutate(wind_beetle  = replace_na(wind_beetle , 0),
         harvest  = replace_na(harvest , 0)) %>% 
  full_join(spruce_df, by = 'id') %>%
  mutate(all_dist_sum        = wind_beetle + harvest,
           cum_removed       = cumsum(all_dist_sum),
           remained_spruce   = spruce_1986 - cum_removed,
           beetle_rate       = wind_beetle/spruce_1986,
           harvest_rate      = harvest/spruce_1986) #%>%

summary(df_merge)

# check how anomalies look like?
ggplot(df_merge, aes(x = year,
                    y = cum_removed)) +
  # geom_smooth() +
  geom_point()


# calculate anomalies ----------------------------------------------------------
reference_period <- 1986:2014

df_anom <- 
  df_merge %>% 
    group_by(id) %>%
  mutate(ref_all_dist       = mean(all_dist_sum[year %in% reference_period], na.rm = T),
         anom_all_dist     = all_dist_sum / mean(all_dist_sum[year %in% reference_period], na.rm = TRUE) - 1,
         ref_wind_beetle    = mean(wind_beetle [year %in% reference_period], na.rm = T),
         anom_wind_beetle  = wind_beetle  / mean(wind_beetle [year %in% reference_period], na.rm = TRUE) - 1,
         ref_harvest        = mean(harvest[year %in% reference_period], na.rm = T),
         anom_harvest      = harvest / mean(harvest[year %in% reference_period], na.rm = TRUE) - 1) %>% 
  dplyr::select(c(id, year, ref_all_dist,anom_all_dist,
                  ref_wind_beetle,anom_wind_beetle,
                  ref_harvest, anom_harvest))
      

# how to handle if I have the same trap (same position) recorded over multiple buffers? if I have a set of buffer 
# per year? I think it is ok, as I will just merge it with beetle counts per year

# check how anomalies look like?
ggplot(df_anom, aes(x = year,
                    y = wind_beetle_z)) +
 # geom_smooth() +
  geom_point()
  
 
  #@filter(remained_forest< 0)
  
df_merge %>% 
 filter(id == 8) 





# Merge ips counts per year with buffer info  ---------------------------------------
# use left_right join to merge individual shifting traps to the buffers! 


# 
df_RS_out <- 
  xy %>% 
  as.data.frame() %>%
  left_join(df_merge, c('id', 'year')) %>% # filter disturbance years that corresponds to trap data/year
  left_join(df_anom, c('id', 'year')) %>% # filter disturbance years that corresponds to trap data/year
  right_join(ips.year.sum, by = c('falsto_name', 'year')) # %>%
 # group_by(falsto_name) %>%
 # mutate(lag1_dist_sum  = dplyr::lag(all_dist_sum, n=1),
 #         lag1_harvest   = dplyr::lag(harvest,      n=1),
 #         lag2_harvest   = dplyr::lag(harvest,      n=2),
  #        lag1_beetles   = dplyr::lag(wind_beetle,      n=1),
  #        lag2_beetles   = dplyr::lag(wind_beetle,      n=2)) 


head(df_RS_out)


save(spruce_df,         # spruce cover by id, from 2017
     dist_df_out,       # disturbances by year
     df_RS_out,            # final table with coordinates, counts and Rs disturbances and present spruce per year'
     df_anom,           # get harvest and wind anomalies per buffer                   
     xy_sf_expand,      # all xy trap coordinates (shifted over years)
     file =  "outData/buffers.Rdata")


