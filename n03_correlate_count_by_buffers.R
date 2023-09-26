

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
#library(tidyverse)
library(lubridate)
#library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)


# load cleaned data
load("outData/ips.Rdata")
load("outData/spatial.Rdata")

# Read data
disturb_year <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/germany114/disturbance_year_1986-2020_germany.tif')
disturb_type <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/storm_fire_harvest/storm_fire_other_classification_germany.tif')
forest_cover <- rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/germany114/forestcover_germany.tif')


#rast(paste(myPath, outFolder,  "disturbance14_conif.tif", sep = "/"))
crs(disturb_year)
crs(disturb_type)
# disturbance coding: 
# 0 = no disturbance
# 1 = other disturbances (not storm- or fire-related)
# 2 = storm disturbances
# 3 = fire disturbances


# Get beetle counts
dat      <- ips.year.sum #fread(paste(myPath, outFolder, "dat.csv", sep = "/"))

# Get spatial data for each trap: each trap has to have only one location per one year! 
# 158 traps in total (79*2), 7
#xy      <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
#                   layer = 'xy_3035') # read watershed
xy_sf_all  # all XY traps
xy_sf_fin  # selected traps with consistent monitoring , one trap have only one XY per all years (I neglected thet they moved over years)


# Get climate data for traps:
#xy_clim <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = "/"))

# Clean trap cooordinates ------------------------------------------
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

















#  ---------------------------------------------------------------------
# How much mortality happened in each buffer? 
# ----------------------------------------------------------------------

# Check the projection systems: 
# Get the correct projection: need to be a projected system! 
# better to project the vector than the raster as the vector maintains the geometry and precision for each vertex: not possible for each pixel


# Project the vector (XY) into raster data
# Check if the CRS is the same
#crs(disturbance) == crs(xy)

# convert to terra object
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

# make a function to export df from each buffer: keep years, disturbance type and forest
extract_rst_val <- function(xy, ...) {
   
   #xy<- xy_ls[[13]]
   id <- xy$id
   
   # get buffer
   buff <- buffer(xy, 500)
   
   # disturbance year
   # crop data and then mask them to have a circle
   dist_year_crop <- crop(disturb_year, buff)
   dist_year_mask <- mask(dist_year_crop, buff)
   
   
   # disturbance type
   dist_type_crop <- crop(disturb_type, buff)
   dist_type_mask <- mask(dist_type_crop, buff)
   
   # Get vector of values
   val_year <- values(dist_year_mask)
   val_type <- values(dist_type_mask)
   
   # count number of cells:
   df <- as.data.frame(table(val_type, val_year))
   
   # add name
   df$id <- id
   
   return(df)
  
  }


# make a function to export df from each buffer: keep years, disturbance type and forest
extract_forest <- function(xy, ...) {
  #xy<- xy_ls[[10]]
  id <- xy$id
  
  # get buffer
  buff <- buffer(xy, 500)
  
  # disturbance year
  # crop data and then mask them to have a circle
  forest_crop <- crop(forest_cover, buff)
  forest_mask <- mask(forest_crop, buff)
  
  # Get vector of values
  val <- values(forest_mask)
  
  # count number of cells:
  df <- as.data.frame(table(val))
  
  # add name
  df$id <- id
  
  return(df)
  
}



# run over all locations
# extract disturbances
out_dist_ls <- lapply(xy_ls, extract_rst_val) #extract_rst_val(xy_ls[[10]])

# merge partial tables into one df
dist_df      <-do.call("rbind", out_dist_ls)

# check data
dist_df[order(dist_df$id),]

#dist_df_wide <- 
  #dist_df %>% 
  #spread(val_type, Freq) %>% 
  
    
table(dist_df$id, dist_df$val_year)

# remove years before study period 
dist_df_out <- 
  dist_df %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate(val_year = as.numeric(as.character(val_year))) %>% # convert factor to numeric
    group_by(val_year, id) %>%
    filter(val_type != 0) %>% # remove 0 = No disturbance!!!!  
    mutate(all_dist_sum = sum(Freq)) %>% 
    spread(val_type, Freq) %>% 
  #filter(val_year  > 2014 ) %>% 
#  filter(val_type == '1') %>%  # 1 = other disturbances (no wind or fire), 2 wind, 3-fire
  dplyr::rename(year = val_year)
# output needed: for every year, I need how many cells was disturbed by what agent
# forest cover = 1, it copntains all disturbances

names(dist_df_out) <- c("year" , "id","all_dist_sum","harvest","wind", "fire")

# investigate teh disturbances per years and trap buffer
ggplot(dist_df_out, aes(x = year,
                        y = all_dist_sum,
                        color = id)) +
  geom_point()


# extract forest ------------------------------------------------
out_forest_ls <- lapply(xy_ls, extract_forest ) #extract_rst_val(xy_ls[[10]])

# merge partial tables into one df
forest_df      <-do.call("rbind", out_forest_ls)

forest_df <- forest_df %>% 
  filter(val == '1') %>%   # 1 = forest from 1986!!!!
  dplyr::select(-val) %>% 
  dplyr::rename(forest_freq =Freq)


# Merge ips counts per year withbuffer info  ---------------------------------------

df_fin <- 
  xy %>% 
  as.data.frame() %>% 
  left_join(forest_df, c('id')) %>% 
  left_join(dist_df_out, c('id', 'year')) %>%
  left_join(ips.year.sum, by = c('falsto_name', 'year'))  %>% 
  group_by(falsto_name) %>% 
  arrange(falsto_name, year) %>% 
   mutate(lag1_dist_sum  = dplyr::lag(all_dist_sum, n=1),
          lag1_harvest   = dplyr::lag(harvest,      n=1),
          lag2_harvest   = dplyr::lag(harvest,      n=2)) 



df_fin %>% 
  ggplot(aes(harvest*0.09))+ 
  geom_histogram(binwidth = 0.1) + 
  facet_wrap(year~.)

df_fin %>% 
  ggplot(aes(sum_ips))+ 
  geom_histogram() + 
  facet_wrap(year~.)


library(ggpmisc)
df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = all_dist_sum )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~. ) +
  ggtitle('All disturbances')
 


df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = harvest*0.09  )) +
 # stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.,) +
  ggtitle('Harvest') + 
  theme_bw()



# need to lag values by one-two year! eg high number of beetles in year n
# triggers high mortality in next one to two years
# export all disturbances = sum, maybe that will work?


df_fin %>%
  as.data.frame() %>% 
  ggplot(aes(x = sum_ips/10000,
             y = lag1_dist_sum  )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  # geom_smooth('lm') +
  facet_wrap(year ~., scales = 'free') +
  ggtitle('Lag1 all disturb')


# buffer area = 780 ha
df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = lag1_harvest*0.09  )) +
  #stat_poly_line() +
 # stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.) + # , scales = 'free'
  ggtitle('Harvest +1year') + 
  theme_bw()



df_fin %>% 
  ggplot(aes(x = sum_ips,
             y = lag2_harvest  )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  # geom_smooth('lm') +
  facet_wrap(year ~., scales = 'free') +
  ggtitle('Lag2 Harvest')















# Get forest type: 
forest_type <- rast(paste(myPath, outFolder, "bav_fortype_ext30_int2u_LZW.tif", sep = "/"))
# code: 
# 2 - coniferous
# 1 - deciduous
# 0 - background

spruce <- rast(paste(myPath, outFolder, "bav_spruce.tif", sep = "/"))
# keep as spruce only the pixel values 88-100 (the % probability of spruce occurence: visually tested with teh coniferous data extent)

# check the projection of teh raster!
crs(forest_type) == crs(disturbance)
crs(spruce) == crs(disturbance)

forest_type <- terra::project(forest_type, crs(disturbance))
spruce <- terra::project(spruce, crs(disturbance))

# reclassify spruce: keep only values 88-100 as 1
unique(values(spruce))
m <- c(1, 90, NA,
       90, 100, 1,
       100, 256, NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
spruce1 <- classify(spruce, rclmat, include.lowest=TRUE)














# get functions to calculate sum mortality and forest composition:
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

# buffer size:
buff_size <- pi*500^2/30^2 # number of pixels by buffer

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
    mutate(freq = species_n / sum(species_n)) # contains 0 as well, so shows the whole buffer
    #mutate(sp_prop = species_n / buff_size)
  
  return(r_df)
}

spruce_comp_by_buff <- function(spatVect, ...) {
  
  #spatVect = buff_ls[[5]]
  # Crop the disturbance raster by the buffer
  r <-crop(spruce1, 
           spatVect)
  
  # Get raster values
  rst_values <- values(r)
  
  # Get a dataframe and count the disturbance cells
  r_df <- data.frame(spruce = rst_values)
  
  r_df <- 
    r_df %>%
    filter(!is.na(r_df)) %>%
    rename(spruce = bav_spruce)  %>% 
    group_by(spruce) %>%
    summarise(spruce_n = n()) %>% 
    #mutate(freq = species_n / sum(species_n)) # contains 0 as well, so shows the whole buffer
  mutate(sp_prop = spruce_n / buff_size)
  
  return(r_df)
}

# Run the function across all buffers: for RS disturbances, for tree species composition, for spruce share
dist_ls <- 
  lapply(buff_ls, mortality_by_buff)

tree_species_ls <- 
  lapply(buff_ls, species_comp_by_buff)


spruce_ls <- 
  lapply(buff_ls, spruce_comp_by_buff)


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

# Add OBJECTID to link data to XY beetle data
spruce_ls2 <- map2(spruce_ls, 
                         buff_500$OBJECTID, ~.x %>% mutate(OBJECTID = .y))

spruce_ls2 <- map2(spruce_ls, 
                         buff_500$globalid, ~.x %>% mutate(globalid = .y))





# Merge all dataframes in a single one:
dist_rs_df      <-do.call("rbind", dist_ls2)
tree_species_df <-do.call("rbind", tree_species_ls2)
spruce_df       <-do.call("rbind", spruce_ls2)

# export data: speies composition
fwrite(tree_species_df, paste(myPath, outTable, 'xy_treeComp.csv', sep = '/'))
fwrite(spruce_df,       paste(myPath, outTable, 'xy_spruce.csv', sep = '/'))


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
