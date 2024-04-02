

#  Explore dependent vs predictors variables
# on raw beetle counts (per month) vs summarized values over year (sum, DOY of max increase, DOY of aggregation)


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')

spring.months         = 3:5
veg.months            = 4:9  # change climate period to fit the data
study.period.extended = 2012:2021


# get libs ----------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
library(rgdal)
library(fasterize)
library(terra)
library(ggplot2)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# colors
library(RColorBrewer)


# load cleaned data
load("outData/ips_counts.Rdata")
load("outData/spatial.Rdata")

# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.aggreg         # get DOY of max increase
# - ips.year.avg       # sum beetles/year, 
# - df.daily           # avg/betles/year per trap
 
# Spatial data: 
# - xy_sf_fin  # XY as sf data, filtered to have only one globalid per trap
# - xy_sf_expand  # XY as sf data, one point for trap per every year, 1106 rows

# Get SPEI and clim data: they are for whole year! check only veg season? now 3-10 (includes march)
df_clim <- fread('outTable/xy_clim_DWD.csv')
df_spei <- fread('outTable/xy_spei_veg_season_DWD.csv')
df_spei_months <- fread('outTable/xy_spei_all_DWD.csv')

anyNA(df_spei_months)

#df_clim$falsto_name <- iconv(df_clim$falsto_name, from = "UTF-8", to = "")

#View(df_spei)

# make a plot with speis
df_spei_months <- df_spei_months %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))


# Calculate median SPEI for each scale
median_spei <- 
  df_spei_months %>%
  group_by(scale, date) %>%
  summarise(median_spei = median(spei), .groups = 'drop') %>% 
  mutate(sign = case_when(median_spei >=0 ~ 'pos',
                          median_spei < 0 ~ 'neg')) %>% 
  pivot_wider(values_from = median_spei, names_from = sign) %>% 
    replace(is.na(.), 0)


# plot SPEI over different scales ----------------------------------------------
ggplot2::ggplot(median_spei) + 
  geom_area(aes(x = date, y = pos), fill = "blue", col = "grey8") +
  geom_area(aes(x = date, y = neg), fill = "red",  col = "grey8") +
   ylab("SPEI") + ggtitle("SPEI") +
  facet_grid(scale~.) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))



# plot barplot of spei - veg season vs no veg season, temperature and precipitation
df_spei_months <- df_spei_months %>% 
  mutate(season = case_when(month %in% veg.months ~ 'veg_season',
                            month %in% c(1,2, 3,10, 11, 12) ~ 'winter'))# %>% 


  
mean.prcp = mean(df_clim[,year %in% 1980:2015]$prcp)
mean.tmp = mean(df_clim$tmp)
  
  p.prcp <- df_clim %>%
    filter(year %in% study.period.extended) %>% 
    ggplot(aes(x = as.factor(year), 
               y = prcp)) +
    geom_boxplot() +
    xlab('') +
    ylab('Monthly precipitation [mm]') + 
    theme(axis.text=element_text(angle = 90, vjust = 0, hjust = 0.5)) +
    geom_hline(yintercept = mean.prcp, lty = 'dashed', col = "grey10" ) 
 
  
  p.tmp <- df_clim %>%
    filter(year %in% study.period.extended) %>% 
    ggplot(aes(x = as.factor(year), 
               y = tmp)) +
    geom_boxplot() +
    xlab('') +
    ylab('Monthly temperature [C]') + 
    theme(axis.text=element_text(angle = 90, vjust = 0, hjust = 0.5)) +
    geom_hline(yintercept = mean.tmp, lty = 'dashed', col = "grey10" ) 
  
 ggarrange(p.prcp, p.tmp) 
  
  


# Get geo data: elevation, slope, roughness...
xy_topo      <- vect(paste(myPath, outFolder, "xy_3035_topo.gpkg", sep = "/"), 
                  layer = 'xy_3035_topo') # read trap pairID


# Get climate data for traps: --------------------------------------------------
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    falsto_name = xy_sf_expand$falsto_name,
                    globalid =    xy_sf_expand$globalid,
                    year =    xy_sf_expand$year)

# remove special characters
xy_df$falsto_name <- gsub("[^A-Za-z0-9_]", "", xy_df$falsto_name) # Remove all non-alphanumeric characters except underscore
xy_topo$falsto_name <- gsub("[^A-Za-z0-9_]", "", xy_topo$falsto_name) # Remove all non-alphanumeric characters except underscore


# get the final subset of globalids
trap_globid <- unique(xy_df$globalid)  # 313, selected
(trap_globid)

# convert to DF
df_topo <- as.data.frame(xy_topo)

dim(df_topo)
df_topo <- df_topo %>% 
  filter(globalid %in% trap_globid) %>% 
  dplyr::select(c('globalid', 'elev'
                  )) %>% 
  full_join(xy_df %>% dplyr::select(falsto_name, globalid))
#dplyr::select(globalid:roughness)


# get average spring temp
# get temp over veg season
# order spei data: spei3, spei12
# spring temperature
df_temp_spring = df_clim %>% 
  filter(year %in% study.period.extended & month %in% spring.months) %>%
  group_by(year, falsto_name) %>% 
  summarise(spring_tmp = mean(tmp))
  
# temperature over vegetation season
df_temp_veg_season = df_clim %>% 
  filter(year %in% study.period.extended & month %in% veg.months) %>%
  group_by(year, falsto_name) %>% 
  summarise(veg_tmp = mean(tmp))

# precipitation
df_prcp_veg_season = df_clim %>% 
  filter(year %in% study.period.extended & month %in% veg.months) %>%
  group_by(year, falsto_name) %>% 
  summarise(veg_prcp = mean(prcp))


# get median spei for all months, for the spearmans correlations
df_spei_year <- df_spei_months %>% 
  group_by(falsto_name, year, scale) %>% 
  summarise(spei = median(spei)) %>% 
  pivot_wider(names_from = scale, values_from = spei) %>% 
  rename(
    annual_spei1 = `1`, 
    annual_spei3 = `3`, 
    annual_spei6 = `6`, 
    annual_spei12 = `12`, 
    annual_spei24 = `24`
  )

head(df_spei)#SPEI is already filtered over years and to vegetation season, just needs to be converted to wide format to havendle two types of SPEI1, 12
# convert to wide format
df_spei_wide <- 
  df_spei %>% 
  pivot_wider(names_from = scale, values_from = spei)# %>% 
  #filter(year %in% study.period.extended) #%>% already set up for veg season

# rename columns
df_spei_wide <- df_spei_wide %>%
  rename(
    spei1 = `1`, 
    spei3 = `3`, 
    spei6 = `6`, 
    spei12 = `12`, 
    spei24 = `24`
  )

# Process predictors  data ------------------------------------


# remove duplicates for proper merging
xy_df          <- unique(xy_df)
df_topo        <- unique(df_topo)
df_temp_spring <- unique(df_temp_spring)

dim(xy_df)
dim(df_topo)


# whole year
df_predictors <- 
  df_temp_spring %>%  
  right_join(df_temp_veg_season, by = c("falsto_name", 'year')) %>%
  right_join(df_prcp_veg_season, by = c("falsto_name",  'year')) %>%
  right_join(df_spei_wide, by = c("falsto_name", 'year')) %>%  # spei vegetation season
 # right_join(df_spei_year, by = c("falsto_name", 'year')) %>%
   # View()
    left_join(xy_df, by = c('falsto_name', 'year')) %>% # add XY coordinates for each trap  and filter the data
    left_join(df_topo, by = c("falsto_name", 'globalid')) %>%
  ungroup(.) %>% 
    dplyr::select(-c(globalid))
 

unique(df_predictors$falsto_name)
# Clean up dependent variables------------------------------

# add also my dependent variables: 
# [1] sum beetles 
ips.year.sum <- ips.year.sum %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) 

# [2] Max Diff DOY
max.diff.doy <- max.diff.doy %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>%  # Remove all non-alphanumeric characters except underscore
  dplyr::select(c(falsto_name, year, doy, diff)) %>% 
  dplyr::rename(peak_doy = doy,
                peak_diff = diff)

# [3] Agg DOY !!!!: only 1080 rows!!! some traps dis not reached the values in time; need to increase threshold or somehow account for this!
ips.aggreg <- ips.aggreg %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>% 
  dplyr::select(c(falsto_name, year, doy)) %>% 
  dplyr::rename(agg_doy = doy )
  
# [4] avg number of beetles per trap/catch
ips.year.avg <- ips.year.avg %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>% 
  ungroup(.) %>% 
  dplyr::select(c(year, falsto_name, avg_beetles_trap))


# make a final table with all Ys, and all Xs to see a correlation matrix. reduce predictors to keep meaningful ones.----------
dat_fin <- 
  df_predictors %>% 
  left_join(ips.year.sum) %>%
  left_join(max.diff.doy) %>%
  left_join(ips.aggreg) %>%
  mutate(pairID = gsub('.{2}$', '', falsto_name)) %>% 
  dplyr::rename(trapID = falsto_name) %>% 
  mutate(pairID = factor(pairID),
         trapID = factor(trapID))

unique(dat_fin$trapID)

# Plots -------------------------------------------------------------------


theme_set(theme_classic())

theme_update(#aspect.ratio = 1, 
                panel.background = element_rect(fill = "white", colour = "black"))

# Function to create a plot
create_plot <- function(data, x_var, y_var, x_label = x_var, y_label = y_var) {
  ggplot(data, aes_string(x = x_var, y = y_var, color = "factor(year)")) +
    geom_point(alpha = 0.5) +
    geom_smooth(aes_string(x = x_var, y = y_var, color = NULL), color = 'black') +
    labs(x = x_label, y = y_label, color = "Year")
}

# Using the function to create plots
p1 <- create_plot(dat_fin, "veg_tmp", "sum_ips")
p2 <- create_plot(dat_fin, "veg_tmp", "peak_diff")
p3 <- create_plot(dat_fin, "veg_tmp", "peak_doy")
p4 <- create_plot(dat_fin, "veg_tmp", "agg_doy")
p5 <- create_plot(dat_fin, "spring_tmp", "sum_ips")
p6 <- create_plot(dat_fin, "spring_tmp", "peak_diff")
p7 <- create_plot(dat_fin, "spring_tmp", "peak_doy")
p8 <- create_plot(dat_fin, "spring_tmp", "agg_doy")

windows()
ggarrange(p1, p2, p3, p4, 
          p5, p6, p7, p8, ncol = 4, nrow = 2, common.legend = TRUE)





# plots Spei ------------------------------------------------------------




# Function to create a plot for given SPEI and y-axis variable
create_spei_plot <- function(data, spei_col, y_var, y_label) {
  spei_col_sym <- sym(spei_col)
  y_var_sym <- sym(y_var)
  
  ggplot(data, aes(!!spei_col_sym, !!y_var_sym, color = factor(year))) +
    geom_point(alpha = 0.5) +
    xlim(-2.5, 1) +
    geom_smooth(aes(!!spei_col_sym, !!y_var_sym, color = NULL), color = 'black') +
    labs(x = spei_col, y = y_label)
}

# Using the function to create plots
p1 <- create_spei_plot(dat_fin, "spei1", "sum_ips", "Sum IPS")
p2 <- create_spei_plot(dat_fin, "spei1", "peak_diff", "Peak diff")
p3 <- create_spei_plot(dat_fin, "spei1", "peak_doy", "Peak DOY")
p4 <- create_spei_plot(dat_fin, "spei1", "agg_doy", "Aggregate DOY")
p5 <- create_spei_plot(dat_fin, "spei3", "sum_ips", "Sum IPS")
p6 <- create_spei_plot(dat_fin, "spei3", "peak_diff", "Peak diff")
p7 <- create_spei_plot(dat_fin, "spei3", "peak_doy", "Peak DOY")
p8 <- create_spei_plot(dat_fin, "spei3", "agg_doy", "Aggregate DOY")
p9 <- create_spei_plot(dat_fin, "spei12", "sum_ips", "Sum IPS")
p10 <- create_spei_plot(dat_fin, "spei12", "peak_diff", "Peak diff")
p11 <- create_spei_plot(dat_fin, "spei12", "peak_doy", "Peak DOY")
p12 <- create_spei_plot(dat_fin, "spei12", "agg_doy", "Aggregate DOY")
p13 <- create_spei_plot(dat_fin, "spei24", "sum_ips", "Sum IPS")
p14 <- create_spei_plot(dat_fin, "spei24", "peak_diff", "Peak diff")
p15 <- create_spei_plot(dat_fin, "spei24", "peak_doy", "Peak DOY")
p16 <- create_spei_plot(dat_fin, "spei24", "agg_doy", "Aggregate DOY")





ggarrange(p1, p2, p3,p4,  
          p5, p6, p7, p8, 
          p9, p10, p11, p12,
          p13, p14, p15, p16,
          ncol = 4, nrow = 4,
          common.legend = TRUE)



# get plots for climatic predictors: -------------------------------
df_clim %>% 
  ggplot(aes(x = year,
             y = tmp,
             group = year)) +
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = df_clim,
    mapping = aes(x = year, 
                  y = tmp  ),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) #+



# plots spei over years
df_clim %>% 
  filter(month %in% veg.months) %>% 
  ggplot(aes(x = year,
             y = tmp,
             group = year)) +
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = filter(df_clim, month %in% veg.months), 
    mapping = aes(x = year, 
                  y = tmp  ),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) #+



# SPEI1
df_spei %>% 
  filter(scale == 24) %>% 
  ggplot(aes(x = year,
             y = spei,
             group = falsto_name)) +
  geom_line(alpha = 0.4)



# SPEI1
df_spei %>% 
 # filter(scale == 12) %>% 
  ggplot(aes(x = year,
             y = spei,
             group = scale)) +
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = df_spei,
    mapping = aes(x = year, 
                  y = spei  ),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
 facet_grid(.~scale)






# save final table --------------------------------------------------------


save(dat_fin, file="outData/final_table.Rdata")





#










