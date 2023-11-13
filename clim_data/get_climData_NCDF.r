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
#library(tidyverse)
library(lubridate)
library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
library(SPEI)

library(lubridate)
#library(easyclimate) # does not work on R 4.1
# WorldClim - ends in 2018!!


# Get spatial data for each trap
xy        <- vect(paste(myPath, outFolder, "xy_fin_years_3035.gpkg", sep = "/")) # read trap location
xy_latlng <- terra::project(xy, 
                          "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# Starting from netCFD: 
# https://cran.r-project.org/web/packages/futureheatwaves/vignettes/starting_from_netcdf.html
# Packages for climate variables:
library(PCICt)
library(zoo)                # for as.Date() specification

# Load ERA5 soil moisture data --------------------------------------------

# Data was downloaded as netCDF from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview

# Read .nc data as a raster in terra - way faster!  
#dat_ras <- terra::rast(paste(myPath, inFolder,  "ERA_Bav_12months.nc", sep = "/"))
dat_ras <- terra::rast('C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/ERA_NET/all_ind/data.nc')

# extract all values to the xy coordinates:
dat_ext_df <- terra::extract(dat_ras, xy_latlng)

# get which ID equals which falsto_name
xy_latlng$ID = 1:nrow(xy_latlng)

xy_names <- data.frame(ID = xy_latlng$ID, 
                       falsto_name = xy_latlng$falsto_name,
                       globalid = xy_latlng$globalid,
                       year = xy_latlng$year)


# Check the variables of interest:
(var_names <- varnames(dat_ras))    # 6 variables
(ras_time  <- time(dat_ras))        # 1584
(ras_source<- sources(dat_ras))        # 1584
(n_layer   <- nlyr(dat_ras))        # 1584 


# terra works well, but require correct naming of variables:

#sort(unique(xy_latlng$OBJECTID)) == sort(unique(dat_ext_df$ID))

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
df <- 
  df_melt %>% 
  arrange(ID, variable) %>%
  mutate(time = rep(ras_time, nrow(dat_ext_df))) %>% 
  separate(variable, 
           c("var", "time_num", 'xx'), "_") %>% 
  dplyr::mutate(year  = lubridate::year(time), 
                month = lubridate::month(time)#, 
                #day   = lubridate::day(time),
                #doy   =  lubridate::yday(time) + 1
                ) # %>%  # as POXIT data has January 1st at 0
  #dplyr::select(-day)




# SPEI: Standardized water precipitation index --------------------------------
SPEI_vars <- c("t2m", "tp")

temp_convert = 273.15  # convert temperature from Kelvin to Celsius


#Get values for SPEI
df1 <-
  df %>% 
  dplyr::filter(var %in% SPEI_vars ) %>% # filter only soil water content
  na.omit() %>% # remove duplicated values
  dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  spread(var, value) %>%
    mutate(t2m = t2m - temp_convert) %>%
    mutate(PET = thornthwaite(t2m, 48.7775),
           BAL = tp - PET) #%>%   # Längengrad Mitte Bayern
 
# test SPEI if improved total precipitation - multiply by 30.41
number_days <- 365/12
df2 <-
  df %>% 
  dplyr::filter(var %in% SPEI_vars ) %>% # filter only soil water content
  na.omit() %>% # remove duplicated values
  dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  spread(var, value) %>%
  mutate(t2m = t2m - temp_convert,
         tp  = tp*number_days*1000) %>% # convert to monthly sums and to mm from meters
  mutate(PET = thornthwaite(t2m, 48.7775),
         BAL = tp - PET) #%>%   # Längengrad Mitte Bayern

  
# calculate SPEI for each XY location:
df_ls <- df2 %>%
    group_split(ID)


# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  
  # get XY name
  id = unique(df$ID)
  
  # convert df to time series
  df.ts <- df %>% 
    ts(df, start = c(1980, 01), end=c(2021,12), frequency=12) 
  
  # Calculate spei or different time intervals:
  my_scales = c(12) # 3,6,12
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
    mutate(ID = rep(id, nrow(out_scales)))
  
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
  

  
# Soil moisture ----------------------------------------------------------


# combine soil moisture data: sum the 4 layers:
soil_vars <- c("swvl1", "swvl2", "swvl3", "swvl4")


# Split df in two tales: having all variables  besides soil water content,
# and only soil water content
df_vars <- df %>% 
  filter(!(var %in% soil_vars))
  

# For soil water content: need to get sums of swvl1:swvl4
df_soil <-
  df %>% 
  dplyr::filter(var %in% soil_vars ) %>% # filter only soil water content
    na.omit() %>% 
    dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  group_by(ID, year, month) %>% 
  summarize(#sum_swv = sum(value),
         sm = sum(value)) 




# Vapour pressure deficit data --------------------------------------------------


# Formula from https://journals.ametsoc.org/view/journals/apme/54/6/jamc-d-14-0321.1.xml?tab_body=pdf
# Based on Allen, R. G., et al. (1998), Crop evapotranspiration. Guidelines for computing crop water requirements. Irrigation and Drainage Paper 56, FAO,Rome, Italy, 27 pp.

vpd_calc <- function(t, td, c1 = 0.611, c2 = 17.67, c3 = 243.5) {
  c1 * exp((c2 * (t - 273.15)) / (t - 273.15 + c3)) - c1 * exp((c2 * (td - 273.15)) / (td - 273.15 + c3))
}


vdp_vars <- c("t2m", "d2m")


# Vapour pressure deficit ------------------------------------------------------
df_vpd <-
  df %>% 
  dplyr::filter(var %in% vdp_vars ) %>% # filter only soil water content
  na.omit() %>% 
  dplyr::select(-c(time_num,time, xx)) %>% # remove unnecessary cols
  spread(var, value) %>% 
  mutate(vpd = vpd_calc(t = t2m, td = d2m))


# Temp -------------------------------------------------------------------------
df_temp <-
  df %>% 
  dplyr::filter(var %in% c('t2m') ) %>% # filter only soil water content
  na.omit() %>% 
  dplyr::select(-c(time_num,time, xx)) %>% # remove unnecessary cols
  spread(var, value) %>% 
  dplyr::rename(tmp = t2m ) %>% 
  mutate(tmp = tmp - temp_convert)  # convert to celsius



# Prec -------------------------------------------------------------------------
df_prcp <-
  df %>% 
  dplyr::filter(var %in% c('tp') ) %>% # total precipitation - here, in one day of a months!
  na.omit() %>% 
  dplyr::select(-c(time_num,time, xx)) %>% # remove unnecessary cols
  spread(var, value) %>% 
  dplyr::rename(prcp = tp )  %>% 
  mutate(prcp = prcp*number_days*1000)  # convert to monthly sums in mm
#mutate(tmp = vpd_calc(t = t2m, td = d2m))




# calculate anomalies ----------------------------------------------------------
reference_period <- 1980:2010
veg.period <- 4:6

df_anom_year <- 
  df_soil %>% 
  left_join(df_vpd) %>% 
  left_join(df_temp) %>% 
  left_join(df_prcp) %>% 
  left_join(df_spei_ID) %>%
 # filter(month %in% veg.period) %>% # filter chunk of veg period
  ungroup(.) %>% 
  group_by(year, ID) %>%
  summarize(sm = mean(sm),
            vpd = mean(vpd),
            tmp = mean(tmp),
            prcp = sum(prcp),  # sum precipitation per year
            spei = mean(spei)) %>%
  ungroup()
  

table(df_anom_year$year, df_anom_year$ID)

df_anom <- 
  df_anom_year %>% 
  group_by(ID) %>%
  mutate(sm_z   = (sm - mean(sm[year %in% reference_period])) / sd(sm[year %in% reference_period]),
         vpd_z  = (vpd - mean(vpd[year %in% reference_period])) / sd(vpd[year %in% reference_period]),
         tmp_z  = (tmp - mean(tmp[year %in% reference_period])) / sd(tmp[year %in% reference_period]),
         prcp_z = (prcp - mean(prcp[year %in% reference_period])) / sd(prcp[year %in% reference_period]),
         spei_z = (spei - mean(spei[year %in% reference_period], na.rm = T)) / sd(spei[year %in% reference_period], na.rm = T)) %>%
  ungroup(.) %>% 
  left_join(xy_names, by = join_by(ID, year)) %>% 
  filter(year %in% c(2015:2021))

length(unique(df_anom_year$ID))  # 1106
length(unique(df_anom_year$year))  # 7 - 1980-2023
range(unique(df_anom$year))  # 2015:2021
length(unique(df_anom$falsto_name))  # 158
length(unique(df_anom$prcp))  # 158
length(unique(df_anom$spei))
length(unique(df_anom$spei_z))

# check some plots---------------------------------------------------------------
p0 <-ggplot(df_anom, aes(x = year,
                         y = spei)) +
  geom_point()+
  geom_smooth() 


p1 <-ggplot(df_anom, aes(x = year,
                    y = spei_z)) +
  geom_point()+
  geom_smooth() 

p2<-ggplot(df_anom, aes(x = year,
                    y = vpd_z)) +
  geom_point()+
  geom_smooth() 


p3<-ggplot(df_anom, aes(x = year,
                    y = tmp_z)) +
  geom_point()+
  geom_smooth() 


p4<-ggplot(df_anom, aes(x = year,
                        y = tmp)) +
  geom_point() +
  geom_smooth()

p5<-ggplot(df_anom, aes(x = year,
                        y = prcp)) +
  geom_point() +
  geom_smooth() 

ggarrange(p0, p1, p2, p3, p4, p5, nrow = 3, ncol = 2)




# add naming to ID and filter data to have only one row per ID ----------------
# set correct encoding for german characters
Encoding(df_anom$falsto_name)        <-  "UTF-8"
Encoding(xy_names$falsto_name) <-  "UTF-8"


# Merge the df_soil and df_vars onto single file
#df_out <- rbind(df_soil, df_vars)
df_anom_named <- df_anom %>% 
  right_join(xy_names, by = join_by(ID, falsto_name, globalid, year)) %>% 
  dplyr::select(-ID)# %>% 
 # unique()

View(df_anom_named)
# Export the final table:

# data.table::fwrite(df_out, 
#                    paste(myPath, outTable, 'xy_clim.csv', sep = "/"))
data.table::fwrite(xy_names, 
                   paste(myPath, outTable, 'xy_clim_IDs.csv', sep = "/"))
data.table::fwrite(df_anom_named, 
                   paste(myPath, outTable, 'xy_anom.csv', sep = "/"))







