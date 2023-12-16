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


rm(list = ls())


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
library(tidyverse)
library(lubridate)
#library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
library(R.utils)
library(stringr)

# Get spatial data for each trap
xy        <- vect(paste(myPath, outFolder, "/xy_all_3035.gpkg", sep = "/"), 
                  layer = 'xy_all_3035') # read trap location
# Convert data to DWD coordinate system:
xy2 <- terra::project(xy, "EPSG:31467")  # coordinate system from the DWD data: Germany


# filter through years: >1970 ------------------------
# 
pattern_years = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[0-1]).*\\.gz$"
# List files: from 2000 onwards:
i = 'temp'
file_ls <- list.files(paste(myPath, "rawData/DeutschWetter", i, sep = "/"),
                      #pattern = "^20.*\\.gz$",
                     # pattern = "^20(1[3-9]|2[01]).*\\.gz$", # match 2013-2019, or 2020-2021
                      pattern = pattern_years,#"^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[01]).*\\.gz$",
                      #pattern = "*\\.gz$",
                      recursive=TRUE)

(file_ls)
# additional filtering to keep only files with '.gz' at the end
file_ls_gz_only <- file_ls[grepl("\\.gz$", file_ls)]
(file_ls_gz_only)

s <- c("jan/193501asc.gz", "feb/188209asc.gz", "mar/197501asc.gz", "apr/202107asc.gz")

f_1970 <- function(x, y = 1970) {
  first4 <- substr(x, 8, 11)
  print(first4)
  #year <- as.numeric(first4)
 # year >= y
}

file_ls[f_1970(file_ls)]

# List all files > than 2000
# C:\Users\ge45lep\Documents\2022_BarkBeetles_Bavaria\rawData\DeutschWetter

# Get vector of folders by climate variable
vars <- c('temp', 'precip')

result_list <- list()

for (i in vars){
  #print(i)
 # vars = 'temp'
  
  # List files: from 2000 onwards:
  file_ls <- list.files(paste(myPath, "rawData/DeutschWetter", i, sep = "/"),
                            #pattern = "^20.*\\.gz$",
                            #pattern = "^20(1[3-9]|2[01]).*\\.gz$",
                        #pattern = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[01]).*\\.gz$",
                        pattern = pattern_years,
                            recursive=TRUE)
  
  file_ls_gz_only <- file_ls[grepl("\\.gz$", file_ls)]
 # (file_ls_gz_only)
  
  
  # read in rasters 
  ras_ls <- lapply(file_ls_gz_only, function(file, ...) {
    # set raster file
    print(file)
    #file = '11_Nov'
    ras_path = paste(myPath, 'rawData/DeutschWetter',  i, file , sep = "/")
    
   # print(file)
    # unzip file
    R.utils::gunzip(ras_path, remove = FALSE, overwrite = T)
    
    # read file
    ff <- gsub(".gz$", "", ras_path)
    # read raster
    z <- rast(ff)
    # set the reference system:
    crs(z) <- "EPSG:31467"
    return(z)
    
  })
  
  #plot(ras_ls[[50]])
  
  # extract raster values to a vector:
  ext_ls <- lapply(ras_ls, function(ras) {
    df <- terra::extract(ras, xy2)
    return(df)
  })
  
  
  # merge data by columns:
  df <- do.call(cbind, ext_ls)
  
  # add the globalID XY data 
  out.df <- cbind(df, 
                  globalid           = xy$globalid,
                  falsto_name        = xy$falsto_name       )
  #names(out.df)
  
  # remove IDs: every raster column now has each ID column: duplicated
  out.df <- out.df %>% 
    dplyr::select(-c(ID))
  
  
  # organize the data in time: 
  # codng for the names XXXXYY - XXXX - year, YY - month
  # data represent the monthly mean values at the grid of 1km2
  # 
  # Mean of the monthly averaged mean daily air temperature in 2 m height above ground, 
  # given in 1/10 °C.
  names(out.df) <- gsub('asc', '', names(out.df))
  names(out.df)
  
  
  # Convert to long format
  long.df <- 
    out.df %>% 
    pivot_longer(!c(globalid, falsto_name), names_to = "time", values_to = 'vals') %>% 
    # split year and months
    mutate(month = as.numeric(str_sub(time, -2, -1)),  # extract last two characters (month)
           year = as.numeric(str_sub(time, 1,4))) %>%      # extract first 4 characters (year)  
    dplyr::select(-c(time))
  
  
  # Export file
  outName = paste0('xy_dwd_', i, '.csv')
  #print(paste(myPath, outTable, outName, sep = '/'))
  fwrite(long.df, paste(myPath, outTable, outName, sep = '/'))
  
  # Store the dataframe in the result list
  result_list[[i]] <- long.df
}



# # For single file: ---------------------------
# Get rasters in .asc format:
# ras_path = paste(myPath, 'rawData/DeutschWetter/01_Jan',  "200501asc.gz", sep = "/")
# R.utils::gunzip(ras_path, remove=FALSE)
# ff <- gsub(".gz$", "", ras_path)
# z <- rast(ff)
# # set the reference system:
# crs(z) <- "EPSG:31467"
# 
# plot(z)
# plot(xy2, add = T)

# merge data together to calculate SPEI -----------------------------
df_prec <- result_list$precip
df_temp <- result_list$temp

df_prec <- df_prec %>% 
  rename(prcp = vals)

df_temp <- df_temp %>% 
  rename(tmp = vals) %>% 
  mutate(tmp = tmp/10)  # correct for DWD recording: values are in 1/10 C format


nrow(df_temp)
nrow(df_prec)

# merge data
df_clim <- df_prec %>% 
  left_join(df_temp, by = join_by(globalid, falsto_name, month, year))


(df_clim)
# get SPEI -----------------------------------------------------------
# SPEI: Standardized water precipitation index --------------------------------
#SPEI_vars <- c("t2m", "tp")

#temp_convert = 273.15  # convert temperature from Kelvin to Celsius


#Get values for SPEI
df1 <-
  df_clim %>% 
 # dplyr::filter(var %in% SPEI_vars ) %>% # filter only soil water content
  na.omit() %>% # remove duplicated values
  #dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  #spread(var, value) %>%
  #mutate(t2m = t2m - temp_convert) %>%
  mutate(PET = thornthwaite(tmp, 48.7775), 
         BAL = prcp - PET) #%>%   # Längengrad Mitte Bayern
  #mutate(ID = 1:nrow(df))

# calculate SPEI for each falsto location:
df_ls <- df1 %>%
  group_split(falsto_name)

test <- df_ls[[2]]

# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  #df <- test
  # get XY name
  id = unique(df$falsto_name)
  
  # convert df to time series
  df.ts <- df %>% 
    ts(df, start = c(1980, 01), end=c(2021,12), frequency=12) 
  
  # Calculate spei or different time intervals:
  my_scales = c(3,12) #c(3,12) # 3,6,12
  spei_ls <- lapply(my_scales, function(s) {
    
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

# summarize spei per year - one SPEI value per year and ID! 
df_spei_ID <-  df_spei_ID %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(scale, date))

# export file:
#fwrite(out.df, paste(myPath, outTable, 'xy_spei.csv', sep = "/"))

# spei does only up to 2021
df_spei_year <- 
  df_spei_ID %>% 
  ungroup(.) %>% 
  group_by(ID, year) %>% 
  summarise(spei = mean(spei)) 









  
# convert to time series data ----------------------------------------
df.ts <- ts(df, start = c(1980, 01), end=c(2021,12), frequency=12) 

head(df.ts)
plot(df.ts)






#  Run examples: --------------------------------------


# calculate drought indices:
# Given a time series of the climatic water balance (precipitation minus potential evapotranspiration), 
# gives a time series of the Standardized Precipitation-Evapotranspiration Index (SPEI).
# SPEI - input: 


# Variables needed:
# YEAR monthly precipitation totals, in mm. 
# MONTH monthlyprecipitation totals, in mm. 
# PRCP monthly precipitation totals, in mm. 
# TMAX monthlymean daily maximum temperature, in ºC. 
# TMIN monthly mean daily minimum temperature, in ºC. 
# TMED monthly mean temperature, in ºC. 
# AWND monthly mean wind speed, in km h-1 
# ACSH monthly mean sun hours, in h. 
# ACSH monthly mean cloud cover, in %.








# Example:

# Load data 
data(wichita) 
head(wichita)
str(wichita)

# Compute potential evapotranspiration (PET) and climatic water balance (BAL) 
wichita$PET <- thornthwaite(wichita$TMED, lat = 37.6475) # 48.777500
wichita$BAL <- wichita$PRCP-wichita$PET 

# Convert to a ts (time series) object for convenience 
wichita <- ts(wichita[,-c(1,2)], end=c(2011,10), frequency=12) 
plot(wichita) 

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

