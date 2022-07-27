# Calculate SPEI data

# at monthly resolution
# for each trap

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
library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
library(R.utils)
library(zoo)

# Get spatial data for each trap
xy        <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                  layer = 'xy_3035') # read trap location
# Convert data to DWD coordinate system:
xy2 <- terra::project(xy, "EPSG:31467")  # coordinate system from the DWD data: Germany


# Get precip and temp data:
df_prec   <- fread(paste(myPath, outTable, 'xy_precip.csv', sep = '/'))
df_temp   <- fread(paste(myPath, outTable, 'xy_temp.csv', sep = '/'))

head(df_prec)
head(df_temp)

# rename column names to join the datasets:
df_prec <- df_prec %>% 
  rename(PRCP = vals)

df_temp <- df_temp %>% 
  rename(TMED = vals) %>% 
  mutate(TMED = TMED/10)  # as in the description

# join data
df <- df_prec %>% 
  full_join(df_temp, by = c("globalid", "month", "year"))


# data check
par(mfrow=c(1,1)) 
plot(df$PRCP, df$TMED)

plot(df$month, df$TMED, pch=".")
plot(df$month, df$PRCP, pch=".")
plot(df$year, df$TMED, pch=".")
abline(lm(df$TMED~df$year), col=2) #  



# Calculate SPEI: -----------------------------------------------------------
# Given a time series of the climatic water balance (precipitation minus potential evapotranspiration), 
# gives a time series of the Standardized Precipitation-Evapotranspiration Index (SPEI).
# SPEI - input: 


# Variables needed:
# YEAR monthly precipitation totals, in mm. 
# MONTH monthlyprecipitation totals, in mm. 
# PRCP monthly precipitation totals, in mm. 
# TMED monthly mean temperature, in ºC. 
# other data are missing: use Thornthwaite transformation for middle Bavaria

# Compute potential evapotranspiration (PET) and climatic water balance (BAL) 
df$PET <- thornthwaite(df$TMED, lat = 48.777500 ) # 48.777500 # Ingolstadt
df$BAL <- df$PRCP-df$PET 

spei11 <- spei(df[,'BAL'], scale = 1) # calculate SPEI for current month

# Convert to a ts (time series) object for convenience
# convert to time series data ----------------------------------------

# calculate SPEI for each location:
df_ls <- df %>%
  group_split(globalid)


# Calculate the SPEI for each location:

get_SPEI <- function(df, ...){
  # get XY name
  id = unique(df$globalid)
  
  # convert df to time series
  df.ts <- df %>% 
    ts(df, start = c(2000, 01), end=c(2021,12), frequency=12) 
  
    # One and tvelwe-months SPEI 
  my_scales = c(1,3,6,12)
  for (scale in my_scales) {
    
  }
  spei1 <- spei(df.ts[,'BAL'], scale = 1) # calculate SPEI for current month
  #spei3 <- spei(df.ts[,'BAL'], scale = 3) # calculate SPEI for current + 2 previous  months
  #spei6 <- spei(df.ts[,'BAL'], scale = 6) # calculate SPEI for current + 6 previous  months
  #spei12 <- spei(df.ts[,'BAL'], scale = 12) # calculates SPEI for current and 11 previous months
  #class(spei1) 
  # extract just values from SPEI object:
  dd <- spei1$fitted
  
  # covert to dataframe
  df.out <- data.frame(spei=as.matrix(dd), 
                       date=zoo::as.Date(time(dd)))
  # add location indication
  df.out <-df.out %>% 
    mutate(globalid = rep(id, nrow(df.out)))
  
  return(df.out)
  
}

# apply over the list of df (locations):

df_ls2<- lapply(df_ls, get_SPEI)

# merge into one file:
out.df <- do.call('rbind', df_ls2)

# export file:
fwrite(out.df, paste(myPath, outTable, 'xy_spei.csv'))











# !!! to complete later!!! Get SPEI for several scales: 1,3,6,12 -------------------------------
spei_ls <- c()
my_scales = c(1,3,6,12)
for (s in my_scales) {
  s = 1
  print(s)
  spei_s <- spei(df.ts[,'BAL'], scale = s) # calculate SPEI for current month
    # extract just values from SPEI object:
  dd <- spei_s$fitted
  
  # covert to dataframe
  df.out <- data.frame(spei=as.matrix(dd), 
                       date=zoo::as.Date(time(dd)))
  # add location indication
  df.out <-df.out %>% 
    mutate(globalid = rep(id, nrow(df.out)))
  
  return(df.out)
  
}



spei1 <- spei(df.ts[,'BAL'], scale = 1) # calculate SPEI for current month
#spei3 <- spei(df.ts[,'BAL'], scale = 3) # calculate SPEI for current + 2 previous  months
#spei6 <- spei(df.ts[,'BAL'], scale = 6) # calculate SPEI for current + 6 previous  months
#spei12 <- spei(df.ts[,'BAL'], scale = 12) # calculates SPEI for current and 11 previous months
#class(spei1) 
# extract just values from SPEI object:
dd <- spei1$fitted




# -----------------------------------------------------------------------
df.ts <- df %>% 
  ts(df, start = c(2000, 01), end=c(2021,12), frequency=12) 

head(df.ts)
plot(df.ts)


# One and tvelwe-months SPEI 
spei1 <- spei(df.ts[,'BAL'], scale = 1) # calculate SPEI for current month
#spei3 <- spei(df.ts[,'BAL'], scale = 3) # calculate SPEI for current + 2 previous  months
#spei6 <- spei(df.ts[,'BAL'], scale = 6) # calculate SPEI for current + 6 previous  months
#spei12 <- spei(df.ts[,'BAL'], scale = 12) # calculates SPEI for current and 11 previous months
class(spei1$fitted) 

# export df:
dd <- spei1$fitted

  data.frame(spei=as.matrix(dd), date=zoo::as.Date(time(dd)))



df.ts=df.ts[order(df.ts$globalid, df.ts$year, df.ts$month,decreasing=F),]
klima_zeit=data.frame()
for (i in unique(klima$Tnr)) {
  klima_zeit_i=ts(klima[klima$Tnr==i,c(6,7)], start = c(1971, 1), end=c(2020,12), frequency=12)
  klima_zeit_i=as.data.frame(klima_zeit_i)
  klima_zeit_i$spei1=spei(klima_zeit_i[,'BAL'], 1)$fitted # one-months SPEI
  klima_zeit_i$spei3=spei(klima_zeit_i[,'BAL'], 3)$fitted
  klima_zeit_i$spei6=spei(klima_zeit_i[,'BAL'], 6)$fitted
  klima_zeit_i$spei12=spei(klima_zeit_i[,'BAL'], 12)$fitted # tvelwe-months SPEI
  klima_zeit_i$spei24=spei(klima_zeit_i[,'BAL'], 24)$fitted
  klima_zeit_i$Tnr=i # Traktnummer aus i
  klima_zeit_i$Jahr=klima$Jahr[klima$Tnr==i]
  klima_zeit_i$Monat=klima$Monat[klima$Tnr==i]
  klima_zeit=rbind(klima_zeit, klima_zeit_i)
}



# SPEI values can be categorized into extremely dry (SPEI ≤ −2), 
# severely dry (−2 < SPEI ≤ −1.5), 
# moderately dry (−1.5 < SPEI ≤ −1), and 
# near normal conditions (−1 < SPEI < + 1) (Slette et al., 2019).

# Extract information from spei object: 
# summary, call function, fitted values, and coefficients 
summary(spei1) 
names(spei1) 
spei1$call 
spei1$fitted 
spei1$coefficients 

# Plot spei object 
par(mfrow=c(2,1)) 
plot(spei1, main='Bavaria, SPEI-1') 
plot(spei12, main='Bavaria, SPEI-12') 


par(mfrow=c(2,1)) 
plot(spi_1, 'Bavaria, SPI-1') 
plot(spi_12, 'Bavaria, SPI-12') 

# Time series not starting in January 
par(mfrow=c(1,1)) 
plot(spei(ts(wichita[,'BAL'], freq=12, start=c(1980,6)), 12)) 

# Using a particular reference period (1980-2000) for computing the parameters 
plot(spei(ts(wichita[,'BAL'], 
             freq=12, 
             start=c(1980,6)), 12, 
          ref.start=c(1980,1), 
          ref.end=c(2000,1))) 










# Example Wichita: -----------------------------------------------------------------------

# Load data 
data(wichita) 

# Compute potential evapotranspiration (PET) and climatic water balance (BAL) 
wichita$PET <- thornthwaite(wichita$TMED, lat = 37.6475) # 48.777500
wichita$BAL <- wichita$PRCP-wichita$PET 

# Convert to a ts (time series) object for convenience 
wichita <- ts(wichita[,-c(1,2)], end=c(2011,10), frequency=12) 
plot(wichita) 

# One and tvelwe-months SPEI 
spei1 <- spei(wichita[,'BAL'], 1) 
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
