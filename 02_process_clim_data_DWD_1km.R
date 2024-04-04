# Get SPEI data

# at daily resolution
# for each trap

# SPEI data monitor the drought worldwide:

# Calculation of the Standardised Precipitation-Evapotranspiration Index
# Global SPEI database


# The Global SPEI database, SPEIbase, offers long-time, 
# robust information about drought conditions at the global scale, 
# with a 0.5 degrees spatial resolution and a monthly time resolution. 
# It has a multi-scale character, providing SPEI time-scales between 1 and 48 months. 
# Currently it covers the period between January 1901 and December 2020.

# Beguería S. (2017) SPEIbase: R code used in generating the SPEI global database, doi:10.5281/zenodo.834462.

gc()
rm(list = ls())

#spring.months         = 3:5
#veg.months            = 4:9 #3:10
#study.period.extended = 2012:2021




# Read my paths -----------------------------------------------------------
source('myPaths.R')


# Read libs  --------------------------------------------------------------

library(SPEI)
library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(raster)
#library(tidyverse)
library(lubridate)
#library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
library(R.utils)
library(stringr)

# SPEI scales:
my_spei_scales = c(1,3, 6, 12, 24) #c(3,12) # 3,6,12
# filter through years: >1980 ------------------------
# 
pattern_years = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[0-1]).*\\.gz$"



# Get spatial data for each trap
#xy        <- vect(paste(myPath, outFolder, "/xy_all_3035.gpkg", sep = "/"), 
#                  layer = 'xy_all_3035') # read trap location
# Convert data to DWD coordinate system:
#xy2 <- terra::project(xy, "EPSG:31467")  # coordinate system from the DWD data: Germany
 

# Get spatial data for each trap
xy        <- vect(paste(myPath, outFolder, "xy_fin_3035.gpkg", sep = "/")) # read trap location
xy2       <- terra::project(xy, "EPSG:31467")  # coordinate system from the DWD data: Germany
xy_latlng <- terra::project(xy2, 
                            "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
unique(xy2$falsto_name)

# remove all special characters
xy2$falsto_name <- gsub("[^A-Za-z0-9_]", "", xy2$falsto_name) # Remove all non-alphanumeric characters except underscore
xy_latlng$falsto_name <- gsub("[^A-Za-z0-9_]", "", xy_latlng$falsto_name) # Remove all non-alphanumeric characters except underscore

length(unique(xy_latlng$falsto_name))

#Encoding(xy2$falsto_name)        <-  "UTF-8"

# get lat long values
crds(xy_latlng)

# keep only usnique localtions - 158 traps
xy2 <- xy2[, !names(xy2) %in% c('year', 'id')] # remove columns
xy2 <- unique(xy2)

# get xy as wgs for teh SPEI
xy2$x_wgs <- crds(xy_latlng)[,'x']
xy2$y_wgs <- crds(xy_latlng)[,'y']

# select just one location to inspect it for SPEI

# Find rows where falsto_name is "Offenhausen_2"
# Find rows where falsto_name is either "Offenhausen_2" or "Blaichach_1"
#rows_to_keep <- which(xy2$falsto_name == "Offenhausen_2" | xy2$falsto_name == "Blaichach_1")

# Subset those rows and keep only 'year' and 'id' columns
#xy2 <- xy2[rows_to_keep, ]









# stack rasters firsst
process_raster_file <- function(file, path_prefix, var) {
  ras_path <- paste(path_prefix, var, file, sep = "/")
  unzipped_path <- gsub(".gz$", "", ras_path)
  
  # Check if unzipped file already exists
  if (!file.exists(unzipped_path)) {
    R.utils::gunzip(ras_path, remove = FALSE, overwrite = TRUE)
  }
  
  z <- rast(unzipped_path)
  crs(z) <- "EPSG:31467"
  return(z)
}


result_list <- list()
vars <- c("temp", "precip")  # Assuming these are your variables

for (i in vars){
  file_ls <- list.files(paste(myPath, "rawData/DeutschWetter", i, sep = "/"),
                        pattern = pattern_years,#"^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[01]).*\\.gz$",
                        recursive=TRUE)
  
  #(file_ls)
  # additional filtering to keep only files with '.gz' at the end
  file_ls_gz_only <- file_ls[grepl("\\.gz$", file_ls)]
  #(file_ls_gz_only)
  
  # Assuming file_ls contains the paths of all raster files
   ras_ls <- lapply(file_ls_gz_only, process_raster_file, path_prefix = paste(myPath, "rawData/DeutschWetter", sep = "/"), var = i)

  ras_stack <- rast(ras_ls)

  #ras_stack <- process_raster_file(file_ls)  # Stack all rasters
  extracted_values <- terra::extract(ras_stack, xy2)

  extracted_values$falsto_name <- xy2$falsto_name
  extracted_values$globalid    <- xy2$globalid

  names(extracted_values) <- gsub('asc', '', names(extracted_values))

  long.df <- extracted_values %>%
    pivot_longer(!c(ID, globalid, falsto_name), names_to = "time", values_to = 'vals') %>%
    mutate(month = as.numeric(str_sub(time, -2, -1)),
           year = as.numeric(str_sub(time, 1, 4))) %>%
    dplyr::select(-c(time))

  result_list[[i]] <- long.df

}

# merge data together to calculate SPEI -----------------------------
df_prec <- result_list$precip
df_temp <- result_list$temp

df_prec <- df_prec %>% 
  rename(prcp = vals)

df_temp <- df_temp %>% 
  rename(tmp = vals) %>% 
  mutate(tmp = tmp/10)  # correct for DWD recording: values are in 1/10 C format


hist(df_temp$tmp)
hist(df_prec$prcp)


# merge data
df_clim <- df_prec %>% 
  left_join(df_temp, by = join_by(globalid, falsto_name, month, year, ID)) %>% 
  right_join(as.data.frame(xy2  ), by = c('falsto_name',  'globalid'                )) %>% 
  dplyr::select(c(falsto_name, prcp, month,  year,   tmp, y_wgs)) %>% 
  distinct()

summary(df_clim)
(df_clim)

unique(df_clim$falsto_name)


# get SPEI -----------------------------------------------------------
# SPEI: Standardized potential evapotranspiration index --------------------------------

# Function to calculate PET using Thornthwaite method
calculate_pet_thornthwaite <- function(temp, latitude) {
  # Using the built-in Thornthwaite function from the SPEI package
    PET <- thornthwaite(temp, latitude)
  return(PET)
}

#Get values for SPEI
df1 <-
  df_clim %>% 
  na.omit() %>% # remove duplicated values
  group_by(falsto_name) %>%
  mutate(PET = calculate_pet_thornthwaite(tmp, unique(y_wgs)[1]),
           BAL = prcp - PET) 

# calculate SPEI for each falsto location:
df_ls <- df1 %>%
  #group_by(falsto_name, .add = TRUE) %>% 
  group_split(falsto_name)

#df.ts$SPEI <- dd




# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  #df <- test
  # get XY name
  id = unique(df$falsto_name)
  
  # convert df to time series
  df.ts <- df %>% 
    arrange(year, month) %>% 
    ts(df, start = c(1980, 01), end=c(2021,12), frequency=12) 
    #ts(df, start = c(2008, 01), end=c(2021,12), frequency=12) 
  
  # Calculate spei or different time intervals:
  spei_ls <- lapply(my_spei_scales, function(s) {
    
    # extract just values from SPEI object:
    dd = spei(df.ts[,'BAL'], scale = s)$fitted
    
    # covert to dataframe, convert date to format
    df.out <- data.frame(spei=as.matrix(dd), 
                         date=zoo::as.Date(time(dd)))
    
    # add scale indication
    df.out <-df.out %>% 
      mutate(scale = rep(s, nrow(df.out)))
    return(df.out)
  })
  # merge scaes tables
  out_scales = do.call('rbind', spei_ls)
  
  # add location indication
  out_scales <-out_scales %>% 
    mutate(falsto_name = rep(id, nrow(out_scales)))
 # (out_scales)
  return(out_scales)
  
}

# apply over the list of df (locations):
df_ls2<- lapply(df_ls, get_SPEI)

# merge into one file:
df_spei_ID <- do.call('rbind', df_ls2)

df_spei_ID2 <- df_spei_ID

df_spei_ID2 %>% 
 # filter(falsto_name == 'Bleichach_1') %>% 
  filter(scale == 12) %>%
  ggplot(aes(x = date, y = spei)) +
  #geom_line() +
  stat_summary(fun = mean, geom = "point", color = "black", size = 1) +
  scale_x_date(limits = as.Date(c("1980-01-01", "2021-12-01"))) # Replace with actual dates


hist(df_prec$prcp)


df_spei_ID2 %>% 
  filter(falsto_name == 'Blaichach_1') %>% 
  filter(scale == 12) %>%
  ggplot(aes(x = date, y = spei)) +
  geom_line() 
  

df_prec %>% 
  group_by(year, falsto_name) %>% 
  summarize(sum_prcp = sum(prcp)) %>%
  View()
  ggplot(aes(x = year,
             y = sum_prcp)) +
  geom_line(alpha = 0.5)

# summarize spei per year - one SPEI value per year and ID! 
df_spei_ID <-  df_spei_ID %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(date))


# plot if correct??
df_spei_ID %>%
  filter(year %in% 2015:2021) %>% 
  ggplot(aes(x = month ,
             y = spei,
             group = factor(year),
             color = factor(year))) +

  #geom_point(alpha = 0.5)  +
  geom_smooth()+
  #geom_point(alpha = 0.5)  +
  facet_grid(year~scale) +
  geom_hline(yintercept = 0, col = 'red', lty = 'dashed') 





# export file:
#fwrite(out.df, paste(myPath, outTable, 'xy_spei.csv', sep = "/"))

# spei does only up to 2021
df_spei_veg_season <- 
  df_spei_ID %>% 
  dplyr::filter(month %in% veg.months & year %in% study.period.extended) %>% 
  ungroup(.) %>% 
  group_by(falsto_name, year, scale) %>% 
  summarise(spei = median(spei)) 




# Export data -------------------------------------------------------------


data.table::fwrite(df_clim, 
                   'outTable/xy_clim_DWD.csv')
data.table::fwrite(df_spei_veg_season, 
                   'outTable/xy_spei_veg_season_DWD.csv')

data.table::fwrite(df_spei_ID, 
                   'outTable/xy_spei_all_DWD.csv')








  
# convert to time series data ----------------------------------------
df.ts <- ts(test, start = c(1980, 01), end=c(2021,12), frequency=12) 

head(df.ts)
plot(df.ts)






#  Run examples: --------------------------------------


# calculate drought indices:
# Given a time series of the climatic water balance (precipitation minus potential evapotranspiration), 
# gives a time series of the Standardized Precipitation-Evapotranspiration Index (SPEI).
# SPEI - input: 


# Variables needed:
# MONTH monthlyprecipitation totals, in mm. 
# PRCP monthly precipitation totals, in mm. 
# TMED monthly mean temperature, in ºC. 




# Example:

# Load data 
data(wichita) 
tail(wichita)
str(wichita)
plot(wichita)

# Compute potential evapotranspiration (PET) and climatic water balance (BAL) 
wichita$PET <- thornthwaite(wichita$TMED, lat = 37.6475) # 48.777500
wichita$BAL <- wichita$PRCP-wichita$PET 

# Convert to a ts (time series) object for convenience 
wichita.ts <- ts(wichita[,-c(1,2)], end=c(2011,10), frequency=12) 
str(wichita.ts)
head(wichita.ts)
plot(wichita.ts) 

# One and tvelwe-months SPEI 
spei1 <- spei(wichita[,'BAL'], 3) 
spei12 <- spei(wichita[,'BAL'], 12) 
class(spei1) 

# Extract information from spei object: 
# summary, call function, fitted values, and coefficients 
summary(spei1) 
names(spei1) 
spei1$call 
spei1$fitted 
spei1$coefficients 

# Plot spei object 
par(mfrow=c(2,1)) 
plot(spei1, main='Wichita, SPEI-1') 
plot(spei12, main='Wichita, SPEI-12') 

# One and tvelwe-months SPI 
spi_1 <- spi(wichita[,'PRCP'], 1) 
spi_12 <- spi(wichita[,'PRCP'], 12) 

par(mfrow=c(2,1)) 
plot(spi_1, 'Wichita, SPI-1') 
plot(spi_12, 'Wichita, SPI-12') 

# Time series not starting in January 
par(mfrow=c(1,1)) 
plot(spei(ts(wichita[,'BAL'], freq=12, start=c(1980,6)), 12)) 

# Using a particular reference period (1980-2000) for computing the parameters 
plot(spei(ts(wichita[,'BAL'], 
             freq=12, 
             start=c(1980,6)), 12, 
          ref.start=c(1980,1), 
          ref.end=c(2000,1))) 


# Using different kernels 
spei24 <- spei(wichita[,'BAL'],24) 
spei24_gau <- spei(wichita[,'BAL'], 24, 
                   kernel=list(type='gaussian', shift=0)) 


par(mfrow=c(2,1)) 
plot(spei24, main='SPEI-24 with rectangular kernel') 
plot(spei24_gau, main='SPEI-24 with gaussian kernel')

