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

library(SPEI)



#  Run examples:


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

# Compute potential evapotranspiration (PET) and climatic water balance (BAL) 
wichita$PET <- thornthwaite(wichita$TMED, lat = 37.6475) 
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
