

##########################################################
#
# 1. Do beetle counts data correlate with RS data??         #
#
##########################################################

# Process:

# get input data
# get a XY data
# make a buffers around XY data
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
library(tidyverse)
library(lubridate)
library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)



# Read data
disturbance <- rast(paste(myPath, outFolder,  "disturbance14_conif.tif", sep = "/"))


# Get beetle counts
dat      <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))

# Get spatial data for each trap
xy      <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                   layer = 'xy_3035') # read watershed

# Get climate data for traps:
xy_clim <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = "/"))



# Get forest type: 
forest_type <- rast(paste(myPath, outFolder, "bav_fortype_ext30_int2u_LZW.tif", sep = "/"))
# code: 
# 2 - coniferous
# 1 - deciduous
# 0 - background

# check the projection of teh raster!
crs(forest_type) == crs(disturbance)

forest_type <- terra::project(forest_type, crs(disturbance))


#  ---------------------------------------------------------------------
# How much mortality happened in each buffer? 
# ----------------------------------------------------------------------

# Check the projection systems: 
# Get the correct projection: need to be a projected system! 
# better to project the vector than the raster as the vector maintains the geometry and precision for each vertex: not possible for each pixel


# Project the vector (XY) into raster data
# Check if the CRS is the same
#crs(disturbance) == crs(xy)
#
# change the projection of the XY data
xy <- terra::project(xy, crs(disturbance))

# Create buffers around XY point: r = 500 m to link the estimate RS damage
buff_500 = xy %>% 
  buffer(500)   # 500 m

#windows()
#plot(xy_500['OBJECTID'])
#plot(xy['OBJECTID'], add = T, col = 'red')



# ---------------------------------------------------------------
#                       Process
# ---------------------------------------------------------------

# Make buffers
# Intersect buffers with forest and disturbance maps
# export as values from rasters to have a dataframe
# for each trap: 

buff_ls <- terra::split(buff_500, "OBJECTID")


mortality_by_buff <- function(spatVect, ...) {

  # Crop the disturbance raster by the buffer
  r <-crop(disturbance, 
           spatVect)

  # Get raster values
  rst_values <- values(r)
  
  # Get a dataframe and count the disturbance cells
  r_df <- data.frame(year = rst_values)
  
  r_df <- r_df %>%
   filter(!is.na(r_df)) %>%
    rename(year = disturbance_map_bavaria) %>% 
   group_by(year) %>%
   tally()
  
  return(r_df)
}


species_comp_by_buff <- function(spatVect, ...) {
  
  # Crop the disturbance raster by the buffer
  r <-crop(forest_type, 
           spatVect)
  
  # Get raster values
  rst_values <- values(r)
  
  # Get a dataframe and count the disturbance cells
  r_df <- data.frame(species = rst_values)
  
  r_df <- r_df %>%
    filter(!is.na(r_df)) %>%
    rename(species = bav_fortype_ext30_int2u_LZW) %>% 
    group_by(species) %>%
    summarise(species_n = n()) %>% 
    mutate(freq = species_n / sum(species_n))
  
  return(r_df)
}



# Run the function across all buffers: for RS disturbances, for tree species composition
dist_ls <- 
  lapply(buff_ls, mortality_by_buff)

tree_species_ls <- 
  lapply(buff_ls, species_comp_by_buff)





# Add OBJECTID to link data to XY beetle data
dist_ls2 <- map2(dist_ls, 
                 buff_500$OBJECTID, ~.x %>% mutate(OBJECTID = .y))

dist_ls2 <- map2(dist_ls2, 
                 buff_500$globalid, ~.x %>% mutate(globalid = .y))



# Add OBJECTID to link data to XY beetle data
tree_species_ls2 <- map2(tree_species_ls, 
                 buff_500$OBJECTID, ~.x %>% mutate(OBJECTID = .y))

tree_species_ls2 <- map2(tree_species_ls2, 
                         buff_500$globalid, ~.x %>% mutate(globalid = .y))




# Merge all dataframes in a single one:
dist_rs_df      <-do.call("rbind", dist_ls2)
tree_species_df <-do.call("rbind", tree_species_ls2)


# Remove the background: 0, keep only deciduous = 1
tree_species_df <- tree_species_df %>% 
  filter(species == 1) 


# Merge RS data to XY trap data:
# get rid of geometry
xy_df <- xy %>%  
  sf::st_as_sf() %>%  # convert first to sf object 
  st_drop_geometry() 


# Link disturbance data to geometry data:
dist_rs_df2 <- dist_rs_df %>% 
  left_join(xy_df, by = c("OBJECTID", 'globalid')) %>% 
 # left_join(tree_species_df, by = "OBJECTID") %>% 
  dplyr::select(-c(bearbeiter,
                   bearbeit_dat, 
                   von_dat,    
                   bis_dat,
                   falsto_name,
                   fk_monsto,
                   bemerk,
                   aelf,
                   pk_globalid))



# ---------------------------------------------------------
# Standardize beetle counts per trap: ---------------------
# ---------------------------------------------------------


# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year  = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day   = lubridate::day(kont_dat),
                doy   = lubridate::yday(kont_dat) + 1)  # as POXIT data has January 1st at 0



# Add drought classification:
dat <- dat %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) 


# Correctly order levels
dat <- dat %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))



# Calculate the mean number of beetles per trap over the season
# corrected for revisit frequency
dat_avg <- dat %>% 
  group_by(art, year, globalid) %>% 
  summarize(sum_beetles = sum(fangmenge, na.rm = T),
            sum_trap    = length(unique(globalid)),
            freq_visit  = length(unique(kont_dat)),
            avg_beetles_trap = sum_beetles/freq_visit) #%>%


# Classify the drought season
dat_avg <- dat_avg %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))



# -----------------------------------------------------------
# Link RS disturbance data to Beetle counts data 
# -----------------------------------------------------------

out_df <- 
  dat_avg %>% 
  left_join(tree_species_df, by = "globalid") %>% 
  left_join(dist_rs_df2, by = c("globalid", 'year', 'OBJECTID')) #%>% 
 

# Link to beetle data, by globalid
head(out_df)


# Change pixel counts to ha:
out_df <- 
  out_df %>% 
  mutate(n = n*0.09,
         species_n = species_n*0.09) %>% # convert to hectares
  rename(RS_area_ha = n,
         deciduous_area_ha = species_n)  %>% 
  mutate(RS_area_ha = replace_na(RS_area_ha, 0))  # replace NA values by 0, as no damage have been detected
  


# Add lagged RS mortality --------------------------------------------
out_df <- 
  out_df %>% 
  #filter(year > 2014) %>% 
  group_by(art, OBJECTID) %>% # , OBJECTID 
  arrange(year, .by_group = TRUE) %>% 
  mutate(lag_dist_sum = lag(RS_area_ha)) #%>% 




# Does the beetle trap data predict tree mortality? 
# for current year
out_df %>% 
  filter(year > 2014 & year < 2021) %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = avg_beetles_trap ,
             y = RS_area_ha,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "gam") +
  facet_wrap(art ~ drought_period, scales = 'free')


# For lagged year
out_df %>% 
  filter(year > 2014 ) %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = avg_beetles_trap ,
             y = lag_dist_sum,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "gam") +
  facet_wrap(art ~ drought_period, scales = 'free')



# Share of deciduous vs betle counts?


out_df %>% 
  filter(year > 2014 ) %>% 
 # filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = freq ,
             y = avg_beetles_trap )) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  facet_grid(.~art) +
  xlab('Deciduous share [%]')
  






##################################################################
#             Statistically evaluate data
##################################################################

# What is the variability between traps?
#                         between traps pairs?
#                         between years?
#                         within traps (repeated measures): longitudinal data
# Dependencies: for one trap over time
#
# split in two data: IT, PC
dat_IT <- dat %>% 
  filter(art == 'Buchdrucker')

#
# split in two data: IT, PC
dat_PC <- dat %>% 
  filter(art == 'Kupferstecher')


# -------------------------------------------

library(nlme)

# Dies the counts change over years?

# lok at counts per years
str(dat_IT)

# Variabilita betweeen years
plot(y = dat_IT$fangmenge, x = as.factor(dat_IT$year ))

# Variability betweeen sites
dat_IT %>% 
ggplot(aes(x = as.factor(doy), 
           y = fangmenge)) +
  geom_boxplot() +
  facet_grid(year~.)


# Check the variability between counts per one trap over time
dat_IT %>% 
  filter(monsto_name == 'Traunstein') %>% 
  ggplot(aes(x = as.factor(doy), 
             y = fangmenge,
             color = as.factor(year))) +
  geom_point()

# Check variability between locations?
dat_IT %>% 
  #filter(monsto_name == 'Traunstein') %>% 
  ggplot(aes(x = as.factor(monsto_name), 
             y = fangmenge)) +
  geom_boxplot()





m1 <- lme(fangmenge ~ 1, dat_IT, random = ~1|year) # Create a linear model: alpha = mean 

coef(m1)

anova(m1)
ggqqplot(residuals(m1))

summary(m1)




# ----------------------------------------------------------------
# LInk beetle counts data to climate XY
# ----------------------------------------------------------------

head(dat)
head(xy_clim)


xy_clim2 <- xy_clim %>%
  dplyr::rename(objectid = ID ) %>% 
  dplyr::select(-c(day, doy))


# convert long to wide format



# Merg data
dat_clim <- dat %>% 
  right_join(xy_clim2, by = c("objectid","year", "month" ))


# relatioship between counts and drought?
dat_clim %>% 
  filter(var == 'swv' & month == 7) %>% 
 ggplot(aes(y = fangmenge,
           x = value)) +
  geom_point() +
  facet_grid(art~year)









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
