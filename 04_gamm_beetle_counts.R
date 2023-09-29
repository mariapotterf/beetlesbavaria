

#  Explore dependent vs predictors variables
# on raw beetle counts (per month) vs summarized values over year (sum, DOY of max increase, DOY of aggregation)


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(raster)
library(rgdal)
library(tidyverse)
library(lubridate)
#library(patchwork)
library(fasterize)
library(terra)
library(ggplot2)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# Stats
library('here')
library('mgcv')
library('gratia')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
library('gganimate')
library('cowplot')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('itsadug')
library(DHARMa)
library(MASS)
library(car)     # for VIF

# colors
library(RColorBrewer)



# load cleaned data
load("outData/ips.Rdata")
load("outData/spatial.Rdata")

# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.aggreg         # get DOY of max increase
# - ips.year           # sum beetles/year, avg/betles/year per trap
 
# Spatial data: 
# - xy_sf_fin  # XY as sf data, filtered to have only one globalid per trap

#!!! need to filter data for my 158 fin traps!!!!!!!!

# Get SPEI and clim data: they are for whole year! check only veg season?
df_spei       <- fread(paste(myPath, outTable, 'xy_spei.csv', sep = '/'))
df_prec       <- fread(paste(myPath, outTable, 'xy_precip.csv', sep = '/'))
df_temp       <- fread(paste(myPath, outTable, 'xy_temp.csv', sep = '/'))
df_clim_ERA   <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = '/')) # need to filter soils volume content index!
df_clim_ERA_ID<- fread(paste(myPath, outTable, 'xy_clim_IDs.csv', sep = '/')) # get the ID to join falsto_name
df_anom       <- fread(paste(myPath, outTable, 'xy_anom.csv', sep = '/')) #


# Get coniferous cover data: coniferous == 2! 
# spruce: share of spruce
df_tree   <- fread(paste(myPath, outTable, 'xy_treeComp.csv', sep = '/'))
df_spruce <- fread(paste(myPath, outTable, 'xy_spruce.csv', sep = '/'))

# Get geo data: elevation, slope, roughness...
xy_topo      <- vect(paste(myPath, outFolder, "xy_3035_topo.gpkg", sep = "/"), 
                  layer = 'xy_3035_topo') # read trap location


# Get climate data for traps: --------------------------------------------------
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_fin)[,"X"],
                    y = sf::st_coordinates(xy_sf_fin)[,"Y"],
                    falsto_name = xy_sf_fin$falsto_name,
                    globalid =    xy_sf_fin$globalid)

# get the final subset of globalids
trap_globid <- unique(xy_df$globalid)  # 158, selected


# convert to DF
df_topo <- as.data.frame(xy_topo)
df_topo <- df_topo %>% 
  filter(globalid %in% trap_globid) %>% 
  dplyr::select(c('globalid', 'elev','slope', 'aspect', 'tpi', 'tri', 'roughness')) %>% 
  full_join(xy_df %>% dplyr::select(falsto_name, globalid))
#dplyr::select(globalid:roughness)


# create ID to match the extracted soil data (NCDF)
#xy_sf_fin$ID = 1:nrow(xy_sf_fin)



# Filter soil drought : soil water content
df_swvl <- df_clim_ERA %>% 
  filter(var == "swv") %>% 
  full_join(df_clim_ERA_ID, by = 'ID')

#View(df_swvl)

length(unique(xy_df$globalid))

# Do neighboring traps correlate in beetle numbers? --------------------------
head(dat.ips.clean)


# Analyze data ------------------------------------

# dplyr::rename column names to join the datasets:
df_prec <- df_prec %>% 
  dplyr::rename(PRCP = vals)

df_temp <- df_temp %>% 
  dplyr::rename(TMED = vals) %>% 
  mutate(TMED = TMED/10)  # as in the description

# select only coniferous: == 2, 0 is no forest (eg. covers the whole 500 m buffer)
df_conif <- df_tree %>% 
  filter(species == 2 ) %>% 
  filter(globalid %in% trap_globid) %>% 
  full_join(xy_df %>% dplyr::select(falsto_name, globalid))

# species = 2 = conif
# species_n = how many cells are coniferous
# freq = what is share of spruce

length(unique(df_conif$globalid))

# plot TEMP and SPEI --------------------------------------------------------
df_spei <- df_spei %>%
  mutate(date = format(as.Date(date, "%d/%m/%Y"), "%m.%Y")) %>%
  separate(date, c('month', 'year')) %>%
  mutate(month = as.numeric(month),
         year = as.numeric(year)) #%>%
  #filter(year > 2014)

reference_period = 1986:2010

df_spei_summer_anom <- df_spei %>% 
  filter(scale == 1) %>%
  filter(month %in% c(6,7,8)) %>%
  group_by(year, globalid) %>%
  summarize(spei = mean(spei)) %>%
  group_by(globalid) %>%
  mutate(spei_z = (spei - mean(spei[year %in% reference_period])) / sd(spei[year %in% reference_period]))# %>%
  

ggplot(df_spei_summer_anom, aes(x = year,
                           y = spei_z)) +
  geom_point()

# Get df SPEI for summer months: scale = 1

# check scales -----------------------------------------------------------------
unique(df_spei$scale)
spei_scale = 12  # visually tested all '1,3,6,12', and the most contrasting is spei12:
# very low 2016, and 2014-2015m 2017-2018

p_spei12<- df_spei %>% 
  filter(scale == spei_scale) %>% 
  ggplot(aes(x = factor(year),
             y = spei)) +
  geom_boxplot(#outlier.shape = NA, #,
               fill = 'grey80' ) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1, lty = 'dashed') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab(paste('SPEI', spei_scale))


# SPEI SEEMS weird!!!!!
# I have the lowest SPEI in 2016, and 2018-2019 seem quite wet!! - remove SPEI for now
df_spei %>% 
  filter(year == 2016 & month == 7)


# SPEI plot on map: ----------------------------------------------
df_spei_avg <- df_spei %>% 
  filter(scale == spei_scale) %>%
  group_by(globalid, year ) %>% 
  dplyr::summarize(spei_med = median(spei))
           

# merge SPEI data
df_spei_avg_sf <- xy_sf_fin %>% 
  right_join(df_spei_avg, 'globalid') 


# map: SPEI:lower number = dryer conditions 
p_spei12_bav <- ggplot(bav_sf) +
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = df_spei_avg_sf,
          aes(color = spei_med)) + # , size = Age, size = 0.8size by factor itself!
   scale_color_distiller(type = 'div', 
                        palette = "RdBu", 
                        direction = 1, 
                        "SPEI") + 
  facet_wrap(.~year) +  # , ncol = 4 
  theme_void() +
  labs(subtitle = '')





# Classify the SPEI data:
# SPEI values can be categorized into extremely dry (SPEI ≤ −2), 
# severely dry (−2 < SPEI ≤ −1.5), 
# moderately dry (−1.5 < SPEI ≤ −1), and 
# near normal conditions (−1 < SPEI < + 1) (Slette et al., 2019).
range(df_spei$spei, na.rm = T)
#  -3.076544  2.980012


# Understand different SPEI:

df_spei %>% 
  filter(globalid == "84617C77-1ECE-4119-8C04-452E41126A05") %>% 
  ggplot(aes(x = month,
             y = spei,
             color = scale)) + 
  geom_point() +
  geom_line() +
  facet_grid(year~scale)

# SPEI 12 seems quite even; find locations with etreme drought to see that?
df_spei %>%  
  filter(spei < -2) 
  


df_spei <- df_spei %>% 
  mutate(spei_cat = case_when(spei <= -2 ~ 'ext_dry',
                              spei > -2 & spei <= -1.5 ~ 'sev_dry',
                              spei > -1.5 & spei <= -1 ~ 'mod_dry',
                              spei > -1 & spei <= +1 ~ 'normal',
                              spei > +1 & spei <= +1.5 ~ 'mod_wet',
                              spei > +1.5 & spei <= +2 ~ 'sev_wet',
                              spei > 2 ~ 'ext_wet'
                              ))

# order spei classification
df_spei <- df_spei %>%
  mutate(spei_cat = factor(
    spei_cat,
    levels = c(
      'ext_dry',
      'sev_dry',
      'mod_dry',
      'normal',
      'mod_wet',
      'sev_wet',
      'ext_wet'
    )
  ))


# convert spei long to wide format, to keep scale values as columns 
df_spei_all <- df_spei %>% 
  pivot_wider(id_cols = !spei_cat, 
              names_from =  scale, 
              values_from = spei, 
              names_prefix = "spei")



# Get frequency of extremely dry locations: -----------------------------------------

# This SPEI has all 12 months!!!
p_spei_freq <- df_spei %>% 
  group_by(year, spei_cat) %>% 
  tally() %>% 
  #filter(spei_cat == 'ext_dry'|spei_cat == 'sev_dry') %>% 
  filter(!is.na(spei_cat     )) %>% 
  ggplot(aes(x = year,
             y = n,
             fill = spei_cat)) +
  geom_col(position = "stack") + 
  scale_fill_brewer(type = 'div', palette = "RdBu", "SPEI") + #"RdBu") #BrBG
 theme_bw()




# join PREC, TEMP, SPEI data --------------------------------------------------
df_swvl_sub <- df_swvl %>% 
  dplyr::select(value, falsto_name, month, year) %>% 
  dplyr::rename(swvl = value)

# add SWVL - now only from 2015, not from 2000 as the others
df_clim_month <- df_prec %>% 
  full_join(df_temp,  by = c("globalid", "month", "year")) %>% 
  full_join(df_spei_all, by = c("globalid", "month", "year")) %>%
  full_join(xy_df, by = c("globalid")) %>% 
  full_join(df_swvl_sub, by = c("falsto_name", "month", "year")) %>%
  dplyr::filter(globalid %in% trap_globid)

length(unique(df_clim_month$globalid))

# Boxplots: TEMP and PRCP -------------------------------------------------
# calculate mean summer season values over year
temp_mean = 
  df_temp %>%
  filter(month %in% 4:10) %>%
  as.data.frame() %>% 
  dplyr::summarize(avg = mean(TMED)) %>%
  pull()

prec_mean = 
  df_prec %>%
  filter(month %in% 4:10) %>%
  as.data.frame() %>% 
  dplyr::summarize(avg = mean(PRCP)) %>%
  pull()


# make plots per year
p_temp <- df_temp %>% 
  filter(month %in% 4:10 & year %in% 2015:2021) %>% 
  ggplot(aes(x = factor(year),
             y = TMED)) +
  geom_boxplot(outlier.shape = 1, 
               fill = 'lightgrey'  ) +
  geom_hline(yintercept = temp_mean, col="red") +
  #geom_hline(yintercept = -1, lty = 'dashed') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))

p_prec <- df_prec %>% 
  filter(month %in% 4:10 & year %in% 2015:2021) %>% 
  ggplot(aes(x = factor(year),
             y = PRCP)) +
  geom_boxplot(outlier.shape = 1, 
               fill = 'lightgrey', 
               size = 0.5 ) +
  geom_hline(yintercept = prec_mean, col="red") +
  #geom_hline(yintercept = -1, lty = 'dashed') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


ggarrange(p_temp, p_prec, p_spei12, ncol = 3)


# DOY aggregation  --------------------------------------------------------

df_swvl_year <- df_swvl %>% 
  ungroup(.) %>% 
  group_by(year, falsto_name) %>% 
  dplyr::summarize(swvl = max(value, na.rm = T)) # , .groups = falsto_name


# merge soils drought with beetle counts
df_swvl_year %>% 
  full_join(dplyr::select(ips.aggreg, 
                          c(doy, year, falsto_name))) %>% 
  ggplot(aes(x = swvl,
             y= doy)) +
  geom_point(alpha = 0.5) +
  facet_wrap(year~.)
  


# Check if there is any trends in beetles numbers? -----------------------------------------
# effect: location, effect within year, between years, between locations:




# Summarize climate per veg season ----------------------------------------------

# avg temp, spei, PRCP :april 31 - oct 30 
df_clim_veg <- df_clim_month %>% 
  dplyr::filter(year > 2014 & year < 2022  ) %>% 
  dplyr::filter(month %in% 4:10) %>% # vegetation period: April 1 - Oct 31
  group_by(falsto_name, year) %>% 
  dplyr::summarise(prec = median(PRCP),
            temp = median(TMED),
            swvl = median(swvl),
            spei1 = median(spei1, na.rm = T),
            spei3 = median(spei3, na.rm = T),
            spei6 = median(spei6, na.rm = T),
            spei12 = median(spei12, na.rm = T)
            )


#### add site conditions -------------------------------------
# have a table with all predictors: as the dependent vriable will be changing (likely)

df_predictors_year <- 
  df_clim_veg %>% 
  left_join(df_conif , by = c("falsto_name")) %>% 
  left_join(df_topo, by = c("falsto_name", 'globalid')) %>% 
  left_join(xy_df, by = c('falsto_name',"globalid")) %>% 
  left_join(df_anom, by = c('falsto_name',"year")) %>% 
  left_join(df_spei_summer_anom, by = c('globalid',"year")) %>% 
  ungroup(.) %>%
  dplyr::select(-c('globalid', 'OBJECTID')) # 


# make lm predictors ----------------------------------------------------------------------
# temp vs elevation - nice trend
df_predictors_year %>% 
  ggplot(aes(x = temp,
             y = elev,
             group = year)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(.~year) + #, scales = 'free'
  theme_bw()


# precip vs elevation - quite ok trend
df_predictors_year %>% 
  ggplot(aes(x = prec,
             y = elev,
             group = year)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(.~year) + #, scales = 'free'
  theme_bw()

# spei12 vs soil water content - weird trend with many locations showing not variation in spei??
df_predictors_year %>% 
  ggplot(aes(x = spei12,
             y = swvl,
             group = year)) +
  stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(.~year) + #, scales = 'free'
  theme_bw()

# precipitation vs spruce share
df_predictors_year %>% 
  ggplot(aes(x = prec,
             y = freq,  # share of spruce
             group = year)) +
  stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(.~year) + #, scales = 'free'
  theme_bw()

# temp vs spruce share
df_predictors_year %>% 
  ggplot(aes(x = temp,
             y = freq,  # share of spruce
             group = year)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(.~year) + #, scales = 'free'
  theme_bw()

# temp vs spruce share - together
# temp vs spruce share
df_predictors_year %>% 
  ggplot(aes(x = temp,
             y = freq #,  # share of spruce
             #group = year
             )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  #facet_wrap(.~year) + #, scales = 'free'
  theme_bw()



# Merge dependent and predictors df ------------------------------------------
# [1] sum beetles year
ips_sum_preds <- 
  ips.year.sum %>% 
  left_join(df_predictors_year, by = join_by(year, falsto_name)) #, by = c("falsto_name", 'globalid', 'year'))


# [2] Max Diff DOY
df_max_diff_doy <-    # merge max difference by doy with the yearly predictors
  max.diff.doy %>% 
  left_join(df_predictors_year, by = join_by(year, falsto_name))
  

# [3] Agg DOY
df_agg_doy <-    # merge max difference by doy with the yearly predictors
  ips.aggreg %>% 
  left_join(df_predictors_year, by = join_by(year, falsto_name))

# LM anomalies ------------------------------------------------------------

ggplot(ips_sum_preds, aes(x = sm_z,
                          y = sum_ips,
                          color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic()# + facet_wrap(.~year)


ggplot(ips_sum_preds, aes(x = vpd_z,
                          y = sum_ips,
                          color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + 
  facet_wrap(.~year)



# earlier peak culmination 
ggplot(df_max_diff_doy, aes(x = sm_z,
                          y = doy,
                          color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + facet_wrap(.~year)


ggplot(df_max_diff_doy, aes(x = sm_z,
                          y = doy,
                          color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + 
  facet_wrap(.~year)


# earlier agregation
ggplot(df_agg_doy , aes(x = sm_z,
                            y = doy,
                            color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + facet_wrap(.~year)


# earlier agregation
ggplot(df_agg_doy , aes(x = spei_z,
                        y = doy,
                        color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + facet_wrap(.~year)


ggplot(df_agg_doy , aes(x = sm_z,
                        y = doy,
                        color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  theme_classic() + facet_wrap(.~year)



# line plot - IPS sum count over years per trap 
p1 <- ips_sum_preds %>% 
  group_by(year, falsto_name) %>% 
  summarize(sum_ips = sum(sum_ips)) %>% 
  ggplot(aes(x = year,
             y = sum_ips/1000,
             color = falsto_name)) +
  geom_line(alpha = 0.2, show.legend = FALSE) +
  stat_summary(fun = median, geom = "path",
               mapping = aes(group = -1), show.legend = FALSE) +
  theme_classic() +
  ggtitle('Sum beetles/year') +
  theme(aspect.ratio = 1, 
        panel.background = element_rect(fill = "white", colour = "black"))

p2<-df_max_diff_doy %>% 
  group_by(year, falsto_name) %>% 
  ggplot(aes(x = year,
             y = doy,
             color = falsto_name)) +
  geom_line(alpha = 0.2, show.legend = FALSE) +
  stat_summary(fun = median, geom = "path",
               mapping = aes(group = -1), show.legend = FALSE) +
  theme_classic() +
  ggtitle('Beetle peak population') +
  theme(aspect.ratio = 1, 
        panel.background = element_rect(fill = "white", colour = "black"))

p3<-df_agg_doy %>% 
  group_by(year, falsto_name) %>% 
  ggplot(aes(x = year,
             y = doy,
             color = falsto_name)) +
  geom_line(alpha = 0.2, show.legend = FALSE) +
  stat_summary(fun = median, geom = "path",
               mapping = aes(group = -1), show.legend = FALSE) +
  theme_classic() +
  ggtitle('Beetle aggregation') +
  theme(aspect.ratio = 1, 
        panel.background = element_rect(fill = "white", colour = "black"))


ggarrange(p1, p2, p3, ncol = 3)




# # variance inflation factor (VIF) ---------------------------------------------
# https://www.statology.org/variance-inflation-factor-r/
v_predictors <- c('year', # if analyse predictors, 'year' needs to be removed! 
                      'temp', 'prec', 'swvl',
                      'spei1', 
                      'spei3', 'spei6', 'spei12', 
                      'elev', 'slope', 'aspect', 'tpi', 'tri', 'roughness')
  

# Compute the VIF values: remove multicollinearity between predictors
# https://www.statology.org/variance-inflation-factor-r/
vif_model_all <- lm(sum_ips ~ temp + elev + prec + swvl+ 
                  spei12 + 
                  spei1 + spei3 + spei6 +  
                  # species_n + 
                  freq+ 
                  slope+  aspect+
                  tpi, # + 
                 tri +
                 roughness, 
                data = ips_sum_preds)


vif_model_fin <- lm(sum_ips ~ temp + elev + prec + swvl+ 
                  spei12 + 
                  #spei1 + spei3 + spei6 +  
                  # species_n + 
                  freq+ 
                  slope+  aspect+
                  tpi, # + 
                #tri +
                #roughness, 
                data = ips_sum_preds)


# get a vector of vif values
vif_values_all <- car::vif(vif_model_all)
vif_values_fin <- car::vif(vif_model_fin)



#create horizontal bar chart to display each VIF value
p <- barplot(vif_values_all, main = "VIF Values", horiz = TRUE, col = "steelblue")

#add vertical line at 5
p_vif_all <- p + abline(v = 5, lwd = 3, lty = 2)

#create horizontal bar chart to display each VIF value
barplot(vif_values_fin, main = "VIF Values", horiz = TRUE, col = "steelblue")

#add vertical line at 5
abline(v = 5, lwd = 3, lty = 2)

# Get overal plots of all predictors by site:
# density plot + points()

summary(df_predictors_doy)



# get correlation coefficient across all predictors:------------------------------
# (to run this, need to remove the y-vars from above: doy, diff, cumsum!)
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
ips2_predictors <- ips_sum_preds %>% 
  ungroup(.) %>%  
  dplyr::select(-c('year', 'falsto_name', 'x', 'y', 'species'))

# Check out which predictors are highly correlated to remove them from teh VIF
preds.res <- cor(ips2_predictors, use = "complete.obs")
round(preds.res, 2)  # simplify visualization

# p-values are missing: get them from different package:
library("Hmisc")
res2 <- Hmisc::rcorr(as.matrix(ips2_predictors))
res2


library(corrplot)

# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")



# Check for scatter plots between precistors:
#install.packages("PerformanceAnalytics")

library("PerformanceAnalytics")
my_data <- mtcars[, c(1,3,4,5,6,7)]
chart.Correlation(ips2_predictors, histogram=TRUE, pch=19)


# export final table ------------------------------------------------------------
write.csv(ips_sum_preds, 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/outTable/ips_sum.csv')


summary(ips_sum_preds)

# Save data ---------------------------------------------------------------

save(ips_sum_preds,           # final df: predictors aggregated by year, sum beetle by year
     df_predictors_year,      # df predictors per year, can be merged with different dependent variables
     p_temp,                  # plot temp
     p_prec,                  # plot prec
     p_spei12,
     p_spei12_bav, 
     p_spei_freq, 
     vif_values_fin,          # VIF values for final set of predictors
     file="outData/predict.Rdata") 






# Do drivers importance change over time? --------------
# ChatGPT: 
# To test if the importance of drivers explaining variability changed over time: 
# perform temporal interaction analysis - 
# allows to assess whether the relationships between the drivers 
# and the response variable differ across different time periods. 

# check hist of conts
windows()
hist(ips_sum_preds$sum_ips)
median(ips_sum_preds$sum_ips)
mean(ips_sum_preds$sum_ips)
sd(ips_sum_preds$sum_ips)


ggplot(ips_sum_preds, aes(x = sum_ips)) + 
  geom_histogram(colour = 4, fill = "white", 
                 bins = 2000)

# check for 0
ips_sum_preds %>% 
  filter(sum_ips == 0) %>% # 0 no!! 
  distinct(falsto_name)
# no 0

fitdistr(ips_sum_preds$sum_ips, 'Poisson')


# poisson: mean and variance are equal
# negative-binomial - allows for overdispersion (variance excees the mean) - not my case

glm1 <- glm(sum_ips ~ temp + prec + elev + temp:as.numeric(year) + spei12,
   data = ips_sum_preds,
   family = "poisson",
   na.action = "na.fail")

# Explore teh model summary:
summary(glm1)
anova(glm1)
MuMIn::dredge(glm1)

# explore residuals
windows()
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm1, 
                                      plot = T)
# Yay! the QQ plot shows that we have some issues, residuals are not independent

# test if the residuals are the same as random?
testDispersion(simulationOutput)




# # try zero inflated model:
# library(pscl)
# 
# install.packages("countreg", repos = "http://R-Forge.R-project.org")
# library("countreg")
# 
# #https://stackoverflow.com/questions/43075911/examining-residuals-and-visualizing-zero-inflated-poission-r
# 
# # Fit a zero-inflated negative binomial model
# m_zero <- zeroinfl(sum_ips ~ temp*year + m_spei3, 
#                    dist = "negbin", data = ips_sum_preds)
# 
# # View the model summary
# summary(m_zero)
# 
# 
# # simply plot effects: 
# plot(allEffects(m_zero))


# Compute the correlation matrix
cor_matrix <- cor(ips_sum_preds[, c(#"year",
                           "temp", 
                           'prec',
                           #"m_spei1",
                           "spei12"#,
                           #"m_spei6",
                           #"m_spei12"
                           )])

# View the correlation matrix
print(cor_matrix)




#AIC(m_zero,glm1)



ggplot(ips_sum_preds, aes(x = year, 
                 y = sum_ips)) + 
  geom_point() +
  stat_smooth(method = "glm")


# Fit a linear regression model for each time period
model <- glm(sum_ips ~ temp:year + elev + prec:year + swvl:year + 
              spei12:year + 
                 freq+ 
                 slope+  
                 aspect+
                 tpi, 
               data = ips_sum_preds,
             "poisson")

# Perform an analysis of variance (ANOVA) to assess the significance of the interaction terms
anova(model)

# Assess the significance of the interaction terms
summary(model)

windows()
plot(model, 2)








unique(dat.ips.clean$falsto_name)












# try the prediction using the raw data --------------------------------

M <- list(c(1, 0.5), NA)
m_counts <- bam(fangmenge ~
            s(spei1, k = 80)+ # + # drought
            s(TMED, k = 30) +  # temperature
            #s(PRCP, k = 20) +         # precip 
            s(freq, k = 50), #+         # spruce %
            #s(month, k =7) +  # months
            #s(year, k = 8) +  # year
            #s(x, y, k = 10, bs = 'ds', m = c(1, 0.5)) + # 2D smooth
            #ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
            #ti(x, y, TMED, d = c(2,1), bs = c('ds','tp'), 
            #   m = M, k = c(25, 10)) +
            #ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
             #  m = M, k = c(25, 15)) +
            #ti(x, y, spei, d = c(2,1), bs = c('ds','tp'),
             #  m = M, k = c(25, 15)),
          data = ips_sum2, 
          method = 'fREML',
          #select = TRUE,
          #family = Tweedie(p=1.1, link = power(0)),
          tw(link = "log"),
         # knots = knots,
          nthreads = 4, 
          discrete = TRUE)

# check for k: edf and k should be balanced ~ 1
m_counts
summary(m_counts)
k.check(m_counts)
plot(m_counts, plot.all = T)


ggplot(data = ips_sum_preds, aes(y = fangmenge, 
                        x = spei1)) +
  geom_point() +
  geom_line(aes(y = fitted(m_counts)),
            colour = "red", size = 1.2) +
  theme_bw()




windows()
appraise(m_counts, method = 'simulate')
plot(m5, page = 1, shade = T)

gam.check(m5)
k.check(m5)
summary(m5)





# Test XY relationship pattern: observed vs modelled data ------------------
# do not run
# linear_model <- gam(Sources ~ SampleDepth, data = isit2)
# summary(linear_model)
# 
# data_plot <- ggplot(data = isit2, aes(y = Sources, x = SampleDepth)) + 
#   geom_point() +
#   geom_line(aes(y = fitted(linear_model)),
#             colour = "red", size = 1.2) + 
#   theme_bw()
# data_plot
# 
# gam_model <- gam(Sources ~ s(SampleDepth), data = isit2)
# 



# get together montly ips counts and spei: ------------------------------------------
ips_sum_month <- dat.ips.clean %>%
  group_by(year, month, falsto_name ) %>% 
  summarise(ips_sum = sum(fangmenge, na.rm = T))


#  Add SPEI, temp, precip, spruce share data --------------------------------------------
# need to check for proper log()!!
ips_sum2 <- 
  ips_sum_month %>%
  left_join(xy_df, by = c('falsto_name')) %>% 
  left_join(filter(df_clim_month, year > 2014), by = c("falsto_name", "year", 'month','globalid','x', 'y'))  %>% # , "month" 
  left_join(df_conif , by = c("falsto_name", 'globalid')) %>% 
  left_join(df_topo, by = c("falsto_name", 'globalid')) %>% #mutate(year = as.factor(as.character(year)))
  dplyr::select(-c('globalid', 'OBJECTID'))

  #mutate(log_sum = log(ips_sum +1 ))  # get log of the monthly values



# Convert to factors to use in GAM: ---------------------------------------
ips_sum2 <- ips_sum2 %>% 
  mutate(#globalid = factor(globalid),
         falsto_name = factor(falsto_name)) #%>% 

# add temporal autocorrelation  
ips_sum2$AR.START <-ips_sum2$year==2015

# use month-year as new variable: 
ips_sum2 <- ips_sum2 %>% 
  mutate(year_month = paste(year, month, sep = '_')) 


ips_sum2 %>% 
  filter(month != 10) %>% 
  ggplot(aes(x = year_month,
             y = ips_sum)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90))


# Investigate the Y distribution -----------------------------------
# check for zeros, for NA
range(ips_sum2$ips_sum)

ips_sum2 %>% 
  filter(ips_sum == 0) %>% 
  tally() %>% 
  distinct(falsto_name) %>% 
  tally()



# What distribution? ---------------------------------------------
  # https://towardsdatascience.com/insurance-risk-pricing-tweedie-approach-1d71207268fc
  # poisson - only for counts
  # gamma - does not take zero values
  # tweedie - can handle zeros values
  
  # Select optimal 'p' for tweetie - varies between 1 to2 ---------------------------------------
m1 <- gam(ips_sum ~ 1,ips_sum2, family = tw() )  
m2 <- gam(ips_sum ~ 1,ips_sum2, family = tw(link = 'log') )  # those area teh same
m3 <- gam(ips_sum ~ 1,ips_sum2, family = tw(link = 'log', a=1.85,b=1.93) )  # those area teh same, low and uppre limit for p optimization
m4 <- gam(ips_sum ~ 1,ips_sum2, family = tw(link = 'log', theta = 1.85, a=1.82,b=1.99) )# the best
m5 <- gam(ips_sum ~ 1,ips_sum2, family = Tweedie(p = 1.9, link = power(0.1)) ) 

appraise(m1)
appraise(m2)
appraise(m3)
appraise(m4)
appraise(m5) # no

# Families in bam: 
# ocat for ordered categorical data.
# tw for Tweedie distributed data, when the power parameter relating the variance to the mean is to be estimated.
# nb for negative binomial data when the theta parameter is to be estimated.
# betar for proportions data on (0,1) when the binomial is not appropriate.
# scat scaled t for heavy tailed data that would otherwise be modelled as Gaussian.
# ziP for zero inflated Poisson data, when the zero inflation rate depends simply on the Poisson mean.
m <- gam(ips_sum ~ 1,ips_sum2, family = nb )  
appraise(m)

m <- gam(ips_sum ~ 1,ips_sum2, family = nb(link = 'log' ))  
appraise(m)

m <- gam(ips_sum ~ 1,ips_sum2, family = nb(link = 'log', theta = 0.5 ))  
appraise(m)

m <- gam(ips_sum ~ 1,ips_sum2, family = nb(link = 'sqrt', theta = 0.5 ))  
appraise(m)

# wrong ones:
#m <- gam(ips_sum ~ 1,ips_sum2, family = scat )  
#appraise(m)

#m <- gam(ips_sum ~ 1,ips_sum2, family = ziP )  
#appraise(m)

# Families to test: counts or aboundances:
# poisson
# negative binomial
# quasi-poisson

# 
#Test negative binomial
m <- gam(ips_sum ~ 1,ips_sum2, family = nb ) # negative binomial ,theta to be estimated
m <- gam(log_sum ~ 1, ips_sum2 ) # gaussian
appraise(m)
  
# Check parameters for Tweedie
  

hist(log(ips_sum2$ips_sum+1))
hist(ips_sum2$freq, breaks = seq(0, 1, 0.1) )
hist(ips_sum2$ips_sum, breaks = seq(0, 55000, 500) )


# Check zero-inflated regression data:
fm_zip



# Get model with gam: ------------------------------------------------

ips_sum2 <- ips_sum2 %>% 
  mutate(monsto_name = gsub('.{2}$', '', falsto_name)) %>% # remove last two characters (indicating trap pair)
  mutate(monsto_name = factor(monsto_name))
# keep only the not correlated predictors:
vif_values

m <- gam(ips_sum ~ 
           s(TMED, k = 5) +
           s(elev, k = 5)+
         s(PRCP, k = 5)+
         s(spei12, k = 5),
         data = ips_sum2,
         family = tw()
         )
appraise(m)

knots_month <- list(doy = c(4.2, 10.5))

M <- list(c(1, 0.5), NA)
m7 <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei12, k = 5) +
    s(TMED, k = 5) +  # temperature
    s(elev, k = 5)  +         # elevation
    s(PRCP, k = 5) +         # precipitation
    s(swvl) +                # soil water volumetric content
    s(freq, k = 5) +         # % coniferous
    s(slope, k = 5) +        # slope
    s(aspect, k = 5) +        # aspect
    s(tpi, k = 5) +             # terrain data
    s(month, k = 7, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(falsto_name, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      year,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 5)
    ),#  +
  # ti(
  #   x,
  #   y,
  #   spei1,
  #   d = c(2, 1),
  #   bs = c('ds', 'tp'),
  #   m = M,
  #   k = c(20, 15)
  # )    ,
  data = ips_sum2,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)





















# Check zero inflated count regession data? -------------------------------------
library(pscl)

# NOT RUN {
## data
data("bioChemists", package = "pscl")

## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois  <- glm(art ~ ., data = bioChemists, family = poisson)
fm_qpois <- glm(art ~ ., data = bioChemists, family = quasipoisson)
fm_nb    <- MASS::glm.nb(art ~ ., data = bioChemists)

## with simple inflation (no regressors for zero component)
fm_zip  <- zeroinfl(art ~ . | 1, data = bioChemists)
fm_zinb <- zeroinfl(art ~ . | 1, data = bioChemists, dist = "negbin")

## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2  <- zeroinfl(art ~ . | ., data = bioChemists)
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
# }






plot(ips_sum2$spei1, ips_sum2$month )


# Get GAM with all predictors: temp, spei, conif, XY ------------------------------------
# check the dstribution:
hist(ips_sum2$ips_sum, breaks = seq(0, 51000, 50))

plot(ips_sum2$spei6, ips_sum2$ips_sum)

# standardize teh data: Z score: mean = 0, sd = 1 ---------------
# will ths lead to better fit? no!
ips_standard <- ips_sum2 %>% 
  mutate(ips_sum = as.vector(scale(ips_sum)),
         PRCP = as.vector(scale(PRCP)),
         TMED = as.vector(scale(TMED)),
         spei1 = as.vector(scale(spei1)),
         spei3 = as.vector(scale(spei3)),
         spei6 = as.vector(scale(spei6)),
         spei12 = as.vector(scale(spei12)),
         freq = as.vector(scale(freq) )) # %>% 
 # filter(!is.na) %>% 
  #complete.cases(.) %>% 
  as.data.frame()



# Test gams with standardized data: ----------------------
m1 <- gam(ips_sum ~ s(spei3),
          data = ips_standard,
          family = tw(link = 'log') )

### GAM standardized y value: -------------
m1 <- gam(ips_sum_stan ~ s(spei3), 
          data = ips_sum2)

length(fitted(m1))
length(ips_standard$spei3)

appraise(m1)

ggplot(ips_sum2, aes(x = spei3, y = ips_sum_stan)) +
  geom_point(alpha=0.7)+
  geom_line(colour = "red", size = 1.2,
            aes(y = fitted(m1))) #+
  geom_line(colour = "blue", size = 1.2,
            aes(y = fitted(m1))) +
  theme_bw()


# GAM: raw sums: monthly sums -------------------------------------------
knots_month <- list(doy = c(4.2, 10.5))

M <- list(c(1, 0.5), NA)
m7 <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 5) +  # temperature
    s(PRCP, k = 5)  +         # precip
     s(freq, k = 5) +         # conif %
   s(sp_prop) +             # spruce proportion
    s(month, k = 7, bs = 'cc') +  # months
     s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    # s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
     ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
  ti(
    x,
    y,
    TMED,
    d = c(2, 1),
    bs = c('ds', 'tp'),
    m = M,
    k = c(20, 10)
  ) +
    ti(
      x,
      y,
      year,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 5)
    ),#  +
    # ti(
    #   x,
    #   y,
    #   spei1,
    #   d = c(2, 1),
    #   bs = c('ds', 'tp'),
    #   m = M,
    #   k = c(20, 15)
    # )    ,
  data = ips_sum2,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)


# !! does not work !!! Test Tweedie and hurdle/zero inflated :from glmmTMB package
library(glmmTMB)
library(MuMIn)

m8. <- glmmTMB(
  ips_sum ~ TMED#,
    ,
  data = ips_sum2)

summary(m8.)

dredge(m8.)








ips_sum2_sept <- filter(ips_sum2, month < 10)

M <- list(c(1, 0.5), NA)
m_rem10 <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 20) +  # temperature
    s(PRCP, k = 20)  +         # precip
    s(month, k = 6, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
     s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      year,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 5)
    ),#  +
  # ti(
  #   x,
  #   y,
  #   spei1,
  #   d = c(2, 1),
  #   bs = c('ds', 'tp'),
  #   m = M,
  #   k = c(20, 15)
  # )    ,
  data = ips_sum2_sept,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)


plot(m_rem10, page = 1, scale = 0, shade = T)
summary(m_rem10)
appraise(m_rem10)

draw(m_rem10, select = 13)  # xy



# m7.noAR --------------------------------------------------------
m7.noAR <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 5) +  # temperature
    s(PRCP, k = 5)  +         # precip
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(month, k = 7, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      PRCP,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 15)
    )  +
    ti(
      x,
      y,
      spei1,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 15)
    ),
  data = ips_sum2,
  #AR.start=AR.START, 
 # rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)

# AIC --------------------------------------------
AIC(m7, m7.noAR )


# m7.sub <- bam(
#   ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
#     #s(spei1, k = 5)  + # drought
#     s(spei3, k = 5) +
#     s(spei6, k = 5) +
#     s(spei12, k = 5) +
#     s(TMED, k = 5) +  # temperature
#     #s(PRCP, k = 5)  +         # precip
#     s(freq, k = 5) +         # conif %
#    # s(sp_prop) +             # spruce proportion
#     s(month, k = 7, bs = 'cc') +  # months
#     s(year, k = 5)+  # year
#     s(globalid, k = 5, bs = 're') + #random effect of trap
#     s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
#     s(
#       x,
#       y,
#       k = 5,
#       bs = 'ds',
#       m = c(1, 0.5)
#     ) + # 2D smooth
#     ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
#     ti(
#       x,
#       y,
#       TMED,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 10)
#     ) +
#     ti(
#       x,
#       y,
#       PRCP,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     )  +
#     ti(
#       x,
#       y,
#       spei1,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     ) +
#     ti(
#       x,
#       y,
#       ,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     )
#   data = ips_sum2,
#   AR.start=AR.START, 
#   rho=0.15,
#   method = 'fREML',
#   family = tw(),
#   knots = knots_month ,
#   nthreads = 4,
#   discrete = TRUE
# )
# 
# 
# 



# compare more complex with less complex model: which one is better?
AIC(m7, m7.sub)


windows()
summary(m7)
plot(m7, page = 1, shade = T, scale = 0)
appraise(m7)
gratia::draw(m7, select = 13)

windows()
gratia::draw(m7, select = 16)


# Make some spatial predictions: !!! does not work!!!
ips_sum3 <- ungroup(ips_sum2)
exp_data <- with(ungroup(ips_sum3),
              expand.grid(month = 6,
                          #DoY = 180,
                          year = seq(min(year), max(year), by = 1),
                          x  = seq(min(x), max(x), length = 100),
                          y  = seq(min(y), max(y), length = 100)))
head(exp_data)
fit <- predict(m7, exp_data)

hist(predict(m7$fitted.values))
save.image("models.RData")




# Test small model:
m8 <- gam(ips_sum ~ 
            s(month, k = 6) +
            s(year, k = 7) + 
            s(TMED) + 
            s(x,y),
          data = ips_sum2,
          family = tw()
            )

summary(m8)
appraise(m8)


m.logy <- bam(
  log_sum  ~
    s(month, k = 6, bs = 'cc') +
    #s(year, k = 7) +
    s(TMED, k = 80) +
    s(PRCP, k = 10) +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(x, y), # +
    data = ips_sum2,
  AR.start = AR.START,
  rho = 0.15,
  family = scat,
  method = 'fREML',
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)

appraise(m.logy)
summary(m.logy)
plot(m.logy, page = 1, )

AIC(m7, m.logy)



# 
exp_data <- with(ips_sum2,
                 expand.grid(month = 6,
                             #DoY = 180,
                             year = seq(min(year), max(year), by = 1),
                             x  = seq(min(x), max(x), length = 100),
                             y  = seq(min(y), max(y), length = 100)))

exp_data<- exp_data %>% 
  left_join(ips_sum2)


exp_data %>% 
  filter(!is.na(globalid))

head(exp_data)
fit <- predict(m8, exp_data)



# with function -----------------------------------------------------
str(mtcars)

a2 <- with(mtcars, mpg[cyl == 8  &  disp > 350])




# Make maps: --------------------------------------------------------------
library(tidyverse)
# ips_sum2 <- ips_sum2 %>% 
#   mutate(x = as.numeric(round(x/100000, 2)),
#          y = as.numeric(round(y/100000, 2)))#

ips_sum2 <- ips_sum2  %>% 
  mutate(across(c(x, y), round, digits = 4))


ggplot(ips_sum2, 
       aes(x = x, y = y)) +
  geom_raster(aes(fill = ips_sum/10000))  + 
  facet_wrap(~ year_month, ncol = 3) # +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma',
                     na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'top', legend.key.width = unit(2, 'cm'))


# Galveston
ggplot(pred.g, 
       aes(x = LONGITUDE, 
           y = LATITUDE)) +
  geom_raster(aes(fill = Fitted)) + facet_wrap(~ YEAR, ncol = 12) +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma',
                     na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'top', legend.key.width = unit(2, 'cm'))

# test for Tweedie family ---------------------------- --------------------

library(mgcv)
set.seed(3)
n<-400
## Simulate data...
dat <- gamSim(1,n=n,dist="poisson",scale=.2)
dat$y <- rTweedie(exp(dat$f),p=1.3,phi=.5) ## Tweedie response

## Fit a fixed p Tweedie, with wrong link ...
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
         data=dat)
plot(b,pages=1)
print(b)
gam.check(b)
gratia::appraise(b, method = 'simulate')

## Same by approximate REML...
b1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
          data=dat,method="REML")
plot(b1,pages=1)
print(b1)
gam.check(b1)
gratia::appraise(b1, method = 'simulate')


## estimate p as part of fitting

b2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=tw(),
          data=dat,method="REML")
plot(b2,pages=1)
print(b2)

gam.check(b2)
gratia::appraise(b2, method = 'simulate')


rm(dat)


# check the k-values: bases functions:
gam.check(m6)
summary(m6)
# plot to have all 4 plots on one page

plot(m5, pages = 1, scheme = 2, shade = TRUE) # , scale = 0

#summary(m4)
gratia::appraise(m5, method = 'simulate')

# check the k value:

gratia::draw(m5, scales = 'free')

plot(m5, pages = 1, scheme = 2, shade = TRUE)




anova(m1, m2, m3)
AIC(m1, m2, m3)




# ---------------------------------------------------------
# test gam models on full catch dataset 
# --------------------------------------------------------

m1 <- gam(fangmenge ~ s(doy, k = 100, bs = 'cc') +
            s(year, k = 8), #+
            #s(globalid, k = 300),
          data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          methods = 'REML')
windows()
plot(m1, shade = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m1, method = 'simulate')

# check the k value:
gam.check(m1)


knots <- list(doy = c(0.5, 365.5))
m2 <- bam(fangmenge ~ s(doy, k = 110, bs = 'cc') +
            s(year, k = 8), # +
          data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)



m3 <- bam(fangmenge ~ s(doy, k = 110, bs = 'cc') +
            s(year, k = 8) +
            s(globalid, bs = 're'),  # 're' = random effect
            data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m3, method = 'simulate')

# check the k value:
gam.check(m3)
plot(m3)
gratia::draw(m3, scales = 'free')

plot(m3, pages = 1, scheme = 2, shade = TRUE)

anova(m1, m2, m3)
AIC(m1, m2, m3)
           
  


# GAM with medians counts per month ---------------------------------------

# still have a lot of wiggliness over the ips catch data:
# maybe I can use the monthly medians?
# will fit with the temperature data as well

# 
ips_med <- ips %>% 
  group_by(year, month, globalid) %>% 
  summarize(fang_med = median(fangmenge, na.rm = T)) %>% 
  filter(year> 2014) %>% 
  as.data.frame()


str(ips_med)

# 
knots_month <- list(doy = c(3.5, 10.5))
m_med1 <- gam(fang_med ~ s(month, k = 7, bs = 'cc'),#+
           # s(year, k = 7), #+
          #s(globalid, k = 300),
          knots = knots_month,    # start and end of the cyclic term 'cc'
          data = ips_med,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          methods = 'REML')

gam.check(m_med1)




windows()
plot(m_med1, shade = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m_med1, method = 'simulate')

# check the k value:
gam.check(m_med1)



m_med2 <- bam(fang_med ~ s(month, k = 7, bs = 'cc')+
            s(year, k = 7), # +
          data = ips_med,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots_month,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)



m_med3 <- bam(fang_med ~ s(month, k = 7, bs = 'cc')+
            s(year, k = 7) +
            s(globalid, k = 500, bs = 're'),  # 're' = random effect
          data = ips_med,
          family = nb, #,nb,#tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots_month,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m_med3, method = 'simulate')
gam.check(m_med3, pages = 1)

summary(m_med3)

# check the k value:
gratia::draw(m_med3, scales = 'free')

plot(m3, pages = 1, scheme = 2, shade = TRUE)

anova(m1, m2, m3)
AIC(m1, m2, m3)
















