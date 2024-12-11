

#  Explore dependent vs predictors variables
# on raw beetle counts (per month) vs 
#summarized values over year (sum, DOY of max increase, DOY of aggregation)


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')


# get libs ----------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
#library(rgdal)
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
df_clim         <- fread('outTable/xy_clim_DWD.csv')
df_spei_season  <- fread('outTable/xy_spei_veg_season_DWD.csv')
df_spei_year    <- fread('outTable/xy_spei_all_DWD.csv')

anyNA(df_spei_year)

df_spei_year %>% 
  ungroup(.) %>% 
  group_by(scale) %>% 
  summarise(mean = mean(spei, na.rm = T),
         sd = sd(spei, na.rm = T),
         median = median(spei, na.rm = T))

#df_clim$falsto_name <- iconv(df_clim$falsto_name, from = "UTF-8", to = "")

#View(df_spei_season)

# make a plot with speis take data fr whole year ----------------------
df_spei_year <- df_spei_year %>% 
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))


# Calculate median SPEI for each scale
median_spei <- 
  df_spei_year %>%
  group_by(scale, date) %>%
  summarise(median_spei = median(spei), .groups = 'drop') %>% 
  mutate(sign = case_when(median_spei >=0 ~ 'pos',
                          median_spei < 0 ~ 'neg')) %>% 
  pivot_wider(values_from = median_spei, names_from = sign) %>% 
    replace(is.na(.), 0)


# Get geo data: elevation, slope, roughness...
xy_topo      <- vect(paste(myPath, outFolder, "xy_3035_topo.gpkg", sep = "/"), 
                  layer = 'xy_3035_topo') # read trap pairID


# Get climate data for traps: --------------------------------------------------
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    falsto_name = xy_sf_expand$falsto_name,
                    globalid =    xy_sf_expand$globalid,
                    year =    xy_sf_expand$year)



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
  summarise(veg_prcp = sum(prcp))


# get median spei for all months, for the spearmans correlations
df_spei_year_median <- df_spei_year %>% 
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

head(df_spei_season)#SPEI is already filtered over years and to vegetation season, just needs to be converted to wide format to havendle two types of SPEI1, 12
# convert to wide format
df_spei_wide <- 
  df_spei_season %>% 
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
unique(df_spei_wide$falsto_name)

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
    dplyr::select(-c(globalid, x, y, elev))
 


# calculate anomalies: for temp and spei ----------------------------------------


df_anom_tmp_prcp_veg_season <- 
  df_clim %>% 
  dplyr::filter(month %in% veg.months) %>% # filter chunk of veg period
  ungroup(.) %>% 
  group_by(year, falsto_name  ) %>%
  summarize(prcp = sum(prcp),
            tmp = mean(tmp)) %>%
  ungroup()

df_tmp_anom <- 
  df_anom_tmp_prcp_veg_season %>% 
  group_by(falsto_name) %>%
  mutate(mean_tmp  = mean(tmp[year %in% reference_period]),
         sd_tmp = sd(tmp[year %in% reference_period]),
         tmp_z  = (tmp - mean(tmp[year %in% reference_period])) / sd(tmp[year %in% reference_period]),
         prcp_z = (prcp - mean(prcp[year %in% reference_period])) / sd(prcp[year %in% reference_period])) %>%
  ungroup(.) 

##### anomalies for spei --------------------------
df_anom_spei_veg_season <- 
  df_spei_year %>% 
  dplyr::filter(month %in% veg.months) %>% # filter chunk of veg period
  dplyr::filter(scale == 3) %>%  # select only spei3
  ungroup(.) %>% 
  group_by(year, falsto_name, scale  ) %>%
  dplyr::summarize(spei = mean(spei)) %>%
  ungroup()

df_spei_anom <- 
  df_anom_spei_veg_season %>% 
  group_by(falsto_name, scale) %>%
  mutate(spei_z  = (spei - mean(spei[year %in% reference_period])) / sd(spei[year %in% reference_period])) %>%
  ungroup(.) 

plot(df_spei_anom$year, df_spei_anom$spei_z)
plot(df_spei_anom$year, df_spei_anom$spei)


#### START

# time series plot: ----------------------------------------
df_anom_all <- df_tmp_anom %>% 
  left_join(df_spei_anom,by = join_by(year, falsto_name)) %>% 
  mutate(cat = case_when(year %in% study_period ~ 'study_period',
                         .default = "ref")) 


# get reference values for tmp and spei
tmp_ref  = df_anom_all %>% 
  ungroup() %>% 
  filter(year %in% reference_period) %>% 
  summarise(tmp = mean(tmp)) %>% 
  pull()

spei_ref  = df_anom_all %>% 
  ungroup() %>% 
  filter(year %in% reference_period) %>% 
  summarise(spei = mean(spei)) %>% 
  pull()


# undersstand what does the z scored mean:
# for TMP, for SPEI, for PRCP
# calculated per each location: eg mean and sd are locatin specific! 
# but i can get average values of Bavaria



# get a table for reference period:
df_anom_all %>% 
  ungroup(.) %>% 
  mutate(type = ifelse( year %in% reference_period, 'ref', 'current')) %>%
  group_by(type) %>% 
  summarise(
    # mean_tmp = mean(tmp, na.rm = TRUE),
    # sd_tmp = sd(tmp, na.rm = TRUE),
    # min_tmp = min(tmp, na.rm = TRUE),
    # max_tmp = max(tmp, na.rm = TRUE),
    # mean_prcp = mean(prcp, na.rm = TRUE),
    # sd_prcp = sd(prcp, na.rm = TRUE),
    # min_prcp = min(prcp, na.rm = TRUE),
    # max_prcp = max(prcp, na.rm = TRUE),
    mean_spei = mean(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE),
    min_spei = min(spei, na.rm = TRUE),
    max_spei = max(spei, na.rm = TRUE)
    
  )



df_anom_all %>% 
  ungroup() %>% 
  filter(!year %in% reference_period) %>% 
  summarise(mean_spei = mean(spei),
            sd_spei = sd(spei),
            mean_tmp = mean(tmp),
            sd_tmp = sd(tmp),
            mean_prcp = mean(prcp),
            sd_prcp = sd(prcp))




p.time.series.spei <- df_anom_all %>% 
  ggplot(aes(x = year,
             y = spei,
             col = cat)) + 
  scale_color_manual(values = c('black', 'red')) +
  geom_rect(
    aes(xmin = 2018, xmax = 2020, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.5, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = spei_ref, lty = 'dashed', col = 'grey50' ) +
  stat_summary(fun.data = "mean_sdl") +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(y =  "SPEI [dim.]", 
       x = '') +
  scale_x_continuous(n.breaks = 5)
p.time.series.spei


p.time.series.tmp <- df_anom_all %>% 
  ggplot(aes(x = year,
             y = tmp,
             col = cat)) + 
  geom_rect(
    aes(xmin = 2018, xmax = 2020, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.5, inherit.aes = FALSE
  ) +
  stat_summary(fun.data = "mean_sdl") +
  geom_hline(yintercept = tmp_ref, lty = 'dashed', col = 'grey50' ) +
  scale_color_manual(values = c('black', 'red')) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(y =  bquote('Temperature [' * degree * 'C]'), 
       x = '') +
  scale_x_continuous(n.breaks = 5)


p.time <- ggarrange(p.time.series.tmp,
                    p.time.series.spei, nrow = 2)

windows(4, 4)
p.time

ggsave(filename = 'outFigs/time_series.png', 
       plot = p.time, width = 4, height = 4, dpi = 300)





# add anoalies to predictors:
df_predictors <- df_predictors %>% 
  left_join(df_spei_anom, by = c('falsto_name', 'year')) %>% 
  left_join(df_tmp_anom,  by = c('falsto_name', 'year'))

plot(df_predictors$year, df_predictors$tmp_z)
plot(df_predictors$year, df_predictors$spei_z)
plot(df_predictors$year, df_predictors$prcp_z)

unique(df_predictors$falsto_name)
# Clean up dependent variables------------------------------

# add also my dependent variables: 
# [1] sum beetles 
#ips.year.sum <- ips.year.sum %>% 
 # mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) 

# [2] Max Diff DOY
max.diff.doy <- max.diff.doy %>% 
  dplyr::select(c(falsto_name, year, doy, diff)) %>% 
  dplyr::rename(peak_doy = doy,
                peak_diff = diff)

# [3] Agg DOY !!!!: only 1080 rows!!! some traps dis not reached the values in time; need to increase threshold or somehow account for this!
ips.aggreg <- ips.aggreg %>% 
  dplyr::select(c(falsto_name, year, doy)) %>% 
  dplyr::rename(agg_doy = doy )
  
# [4] avg number of beetles per trap/catch
ips.year.avg <- ips.year.avg %>% 
  ungroup(.) %>% 
  dplyr::select(c(year, falsto_name, avg_beetles_trap))


# make a final table with all Ys, and all Xs to see a correlation matrix. reduce predictors to keep meaningful ones.----------
dat_fin <- 
  df_predictors %>% 
  left_join(ips.year.sum, by = join_by(year, falsto_name)) %>%
  #View()
  left_join(max.diff.doy,  by = join_by(year, falsto_name)) %>%
  left_join(ips.aggreg, by = join_by(year, falsto_name)) %>%
 # View()
  mutate(pairID = gsub('.{2}$', '', falsto_name)) %>% 
  dplyr::rename(trapID = falsto_name) %>% 
  mutate(pairID = factor(pairID),
         trapID = factor(trapID))

View(dat_fin)

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








# save final table --------------------------------------------------------
save(dat_fin, file="outData/final_table.Rdata")










