

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LInk trap data with damage volume and RS observation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Goal: do beetle counts/dynamics preds (+ lags) predict tree damage?  

# from LWF: the damage volume+ geometries
#  - read tables
#  - merge geometries
#  explore damage per year
# RS disturbance and forests: 
#  - read RS data from Cornelius: 
# - do amage data correlate with RS beetle estimation? YES!
# beetle counts data:
# - how to link them with beetle trap data?
# - select the respective region per year
# - select in terra oes not work well: convert location to raster, and extract forest_reviewr1 for trap data
# - merge by years with beetle counts and damage
# - do beetle counts fit over forest damage??

# to do:
# read data
# expland teh geometries for teh damage data, to plot the damage with ggplot
# how to extrapolated beetle counts from the point-base data? - create a map and overlay this with the reporte damage data
# use a spruce layer for it
# 
  

rm(list=ls()) 

# Libs --------------------------------------------------------------------------
library(dplyr)

library(sf)
library(raster)
library(terra)


library(data.table)
library(ggplot2)
library(ggpubr)
#library(ggspatvector)

library(lubridate)
library(fasterize)
library(tidyr)


# Stats
#library('here')
#library('gamair')
library('purrr')
#library('mvnfast')
library("tibble")
library('cowplot')
library('tidyr')
library("knitr")

library('readr')

library(ggeffects)

library(MASS)
library(car)     # for VIF
library(glmmTMB) # many families



library(lme4) #pseudo-R2 and and mixed model functionality
library(MuMIn) #dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc) #pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa) #model diagnostics
library(effects) #what do my marginal effects look like?
library(performance) #binomial model diagnostics
library(emmeans) #post hoc for categorical predictors

# test for autocorrelation
library(lmtest)

# colors
library(RColorBrewer)

library(mgcv)
library(mgcViz)
library(lme4)
library(gratia)
#library(knitr)   # for table outputs

library(forcats)  # reorder by factor

library(stringr) 



# Data --------------------------------------------------------------------------
### read beetle data 
# load cleaned data
load("outData/ips_counts.Rdata")
load("outData/spatial.Rdata")
dat_dynamics <- fread( 'outTable/beetle_dynamic_indicators.csv')


source("my_functions.R")
# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.aggreg         # get DOY of max increase
# - ips.year.avg       # sum beetles/year, 
# - df.daily           # avg/betles/year per trap

# Spatial data: 
# - xy_sf_fin  # XY as sf data, filtered to have only one globalid per trap
# - xy_sf_expand  # XY as sf data, one point for trap per every year, 1106 rows



### read tree damage volume & districts shp 
sf_districts <-  st_read('rawData/damage_volume/AELF_Revier_20240502/AELF_Revier_bis20240424_LWF_20240502.shp')
df_damage    <-  fread('rawData/damage_volume/Bavaria_TreeDamageData_IpsTypographus_2015-2021.csv', dec = ',')

### disturbance rasters 
disturb_year <- rast('rawData/germany114/disturbance_year_1986-2020_germany.tif')
disturb_type <- rast('rawData/fire_wind_barkbeetle_germany.tif')
forest       <- rast('rawData/germany114/forestcover_germany.tif')

crs(disturb_year) <- "EPSG:3035"
# check data: how many regions I have?
sf_ID   <- unique(sf_districts$forstrev_1)  # 336 districts numbers
df_ID   <- unique(df_damage$AELF_district_ID)  # 337 district numbers

#length(sf_ID)
#length(df_ID)



#### Process beetle sums/year/trap --------------------------------------------
# add beetle counts to trap locations per year 
sum_ips <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_ips = sum(fangmenge, na.rm = T)) 

# convert to sf object, crs 3035
sum_ips_sf <-  xy_sf_expand %>%  left_join(sum_ips, by = c("falsto_name", 'year')) 
#st_is_valid(sum_ips_sf_trans)

# change projection & convert to terra
my_crs = st_crs(disturb_year)

# change projection & convert to terra
# districts = polygons
sf_districts_trans <- st_transform(sf_districts, crs = my_crs) # st_crs(disturb_year))
v_districts        <-  terra::vect(sf_districts_trans) # convert to terra

sum_ips_sf_trans  <- st_transform(sum_ips_sf, crs = my_crs)
v_sum_ips         <- vect(sum_ips_sf_trans)


# works well for terra! 
plot(forest)
plot(disturb_type)
plot(v_sum_ips, add = T)



# simplify terra object:
v_districts_simpl <- simplifyGeom(v_districts, tolerance=20, preserveTopology=TRUE, makeValid=TRUE)
v_sum_ips_simpl   <- simplifyGeom(v_sum_ips, tolerance=20, preserveTopology=TRUE, makeValid=TRUE)


writeVector(v_sum_ips_simpl, 'outSpatial/XY_beetle_sum_years.gpkg', overwrite=TRUE) 

#  does simplification work? yes!!
plot(v_districts_simpl)
plot(v_sum_ips_simpl, add = T, col = 'red')

identical(crs(v_districts_simpl), crs(v_sum_ips_simpl))


# Read rasters ------------------------------------------------------------
# there is an aditiona row in the damage data with district ID 0 - I can remove that


###### crop  rasters to Bavaria --------------------------
dist_year_crop <- terra::crop(disturb_year, v_districts_simpl)
dist_type_crop <- terra::crop(disturb_type, v_districts_simpl)
forest_crop    <- terra::crop(forest,       v_districts_simpl)

# process spatial data: check geometries, change projection to Corneliuses data, convert to raster the AELF districts
crs(disturb_type) <- crs(dist_year_crop)
crs(forest_crop)  <- crs(dist_year_crop)


### Process department shp and extract revier ID  -------------------------------
sf_simpled <- st_as_sf(v_districts_simpl) #st_simplify(sf_districts_trans, dTolerance = 20, preserveTopology = TRUE)

identical(crs(sf_simpled), crs(v_districts_simpl))

# get extend and set correct projectipn
ext <- as(extent(sf_simpled), 'SpatialPolygons')  # treat raster as from raster, not terra package
ext <- st_as_sf(ext)

identical(crs(sum_ips_sf_trans), crs(ext))

# efine teh crs first - it was NA
st_crs(ext) <- my_crs# st_crs(raster(disturb_year))

#ext  <- st_transform(ext, crs = my_crs)

st_crs(sf_simpled) == st_crs(ext) 
st_crs(disturb_year) == st_crs(ext) 
st_crs(sum_ips_sf_trans) == st_crs(ext) 

# create a base raster to rasteroize to it
base_raster = raster(dist_year_crop)

# Create grid index for each pixel
grid_sel     <- st_intersection(sf_simpled, ext)
grid_sel_ras <- fasterize::fasterize(grid_sel,
                                     base_raster,
                                     #raster(dist_year_crop), 
                                     field = "forstrev_1") # name the new rasters as polygon ids  #, field = "forstrev_1"

# assigne the crs manually! differet assigenement for vector and raster
crs(grid_sel_ras)              <- crs(sf_simpled)
st_crs(sum_ips_sf_trans)       <- crs(sf_simpled)

# Extract CRS again after transformation
crs_raster <- crs(grid_sel_ras)
crs_vector <- st_crs(sum_ips_sf_trans)

# Confirm they are now the same
identical(as.character(crs_raster), crs_vector$wkt)  # Should return TRUE


# check if the crs fits! YES
st_crs(grid_sel_ras)       == st_crs(sf_simpled)
st_crs(grid_sel_ras)       == st_crs(sum_ips_sf_trans)


# extract the 'reviers ID' (from the districts) to the beetle counts data 
district_ID <- terra::extract(rast(grid_sel_ras),
                              v_sum_ips_simpl,
                            #  vect(sum_ips_sf_trans), 
                              df = TRUE, bind = TRUE )

#### Link beetle counts vs damage volume per department -------------------------------------------------------
ips_damage_district <- 
  district_ID %>% 
  as.data.frame() %>% 
  dplyr::rename(., forstrev_1 = layer)


# join trap counts with volume damage data and beetle dynamics indicators;
# average per trap pair to not have duplicate records
ips_damage_merge <- 
  ips_damage_district %>% 
  left_join(df_damage, by = c('forstrev_1' = 'AELF_district_ID',
                              'year' = 'Year')) %>% 
  dplyr::select(-c(id, AELF_district_name,AELF_name, AELF_ID,
                 pest_species, tree_species, unit,
                 globalid )) %>% 
  dplyr::rename(trapID = falsto_name) %>% 
   # add beetle dynamics variables 
  left_join(dat_dynamics) %>%
  group_by(trapID)  %>% 
  dplyr::arrange(year,.by_group = TRUE) %>% 
   mutate(sum_ips_lag1  = dplyr::lag(sum_ips, n = 1, default = NA),  # 
         sum_ips_lag2  = dplyr::lag(sum_ips, n = 2, default = NA),  # , order_by = year
         agg_doy_lag1  = lag(agg_doy, n = 1, default = NA),
         agg_doy_lag2  = lag(agg_doy, n = 2, default = NA),
         peak_doy_lag1  = lag(peak_doy, n = 1, default = NA),
         peak_doy_lag2  = lag(peak_doy, n = 2, default = NA),
         peak_diff_lag1  = lag(peak_diff, n = 1, default = NA),
         peak_diff_lag2  = lag(peak_diff, n = 2, default = NA),
         lag1_damaged_volume_total  = lag(damaged_volume_total_m3, n = 1, default = NA)
        )  %>%
  
  ungroup(.) %>% 
  mutate(pairID = as.factor(pairID),
         trapID = as.factor(trapID)) #%>% 
    
  

# Check if lags are correct!!!! YES!!!

ips_damage_merge  %>% 
   dplyr::filter(trapID == "Anzinger_Forst_1") %>% 
     dplyr::select(trapID, year, sum_ips, sum_ips_lag1, sum_ips_lag2,
                   agg_doy, agg_doy_lag2, damaged_volume_total_m3, lag1_damaged_volume_total )
  


# remove NA values
ips_damage_clean <- na.omit(ips_damage_merge)



# Extract RS data   -------------------------------------------------------------
# merge disturbance rasters together; check n of cells 
ncell(grid_sel_ras)
ncell(dist_year_crop)
ncell(dist_type_crop)
ncell(forest_crop)


# extract values:
ID     = terra::values(grid_sel_ras)
year   = terra::values(dist_year_crop, mat = F)
type   = terra::values(dist_type_crop, mat = F)
forest = terra::values(forest_crop, mat = F)

# merge into one df
df_disturb <- data.frame(ID = ID,
                         year = year,
                         type = type,
                         forest = forest
                         ) %>% 
  dplyr::filter(!is.na(ID))

# What to get: per district: how many cells, forest, disturbances, harvest& bark beetles
df_forest <- df_disturb %>% 
  group_by(ID) %>%
  dplyr::summarize(forest86 = sum(forest == 1, na.rm = TRUE),
                   area_cells = n())

# RS get disturbances and type 2015-2020
df_disturb_sum <- 
  df_disturb %>% 
  dplyr::filter(year %in% 2015:2021) %>% 
  dplyr::select(!forest) %>% 
  dplyr::filter(!is.nan(type)) %>%        # Exclude rows where 'type' is NaN
  group_by(ID, year, type) %>% 
    dplyr::summarize(count = n(), .groups = "drop") %>% 
  pivot_wider(
    names_from = type,          # Create new column names from the 'type' values
    values_from = count,        # Fill these new columns with values from 'count'
    names_prefix = "type_",     # Optional: prefix for new column names for clarity
    values_fill = list(count = 0)  # Fill missing values with 0
  )


# merge forest and disturbances from RS data per ASFL districts
df_RS <- df_disturb_sum %>% 
  full_join(df_forest, by = join_by(ID)) %>% 
  dplyr::rename(RS_wind_beetle = type_1,
                RS_fire    = type_2,
                RS_harvest = type_3) %>% 
  mutate(RS_sum = RS_wind_beetle + RS_fire +RS_harvest) %>% 
  dplyr::select(c(ID,  year,  RS_wind_beetle)) %>% 
  group_by(ID)
  

# type_1 = wind & beetles
# type_2 = fire
# type_3 = harvest



#### Merge RS with damage data -----------------------------------------------
df_all <- df_RS %>% 
  full_join(df_damage, by = c("ID"   = "AELF_district_ID",
                              "year" = "Year")) %>% 
 # dplyr::filter(year %in% 2015:2020) %>% 
  #dplyr::filter(RS_wind_beetle <2000) %>% 
  dplyr::select(!c( AELF_ID, tree_species, unit, pest_species ))# %>%  # keep: AELF_district_name, AELF_name,



# Get a summary table for year & per district --------------------------------------

# Represent results using quantiles, as they are skewed?
#qntils = c(0, 0.25, 0.5, 0.75, 1)
df_all %>% 
  ungroup() %>% 
  summarize(sum_volume = sum(damaged_volume_total_m3),
            sum_RS = sum(RS_wind_beetle, na.rm = T)*0.09)

# data per years ----------------------------------------------------------------
observation_mortality_year <- 
  df_all %>% 
  ungroup(.) %>% 
  filter(year %in% 2015:2021) %>%
    dplyr::select(ID, year, damaged_volume_total_m3, RS_wind_beetle) %>% 
    dplyr::mutate(RS_wind_beetle = RS_wind_beetle*0.09) %>% # convert pixels to hectares
    group_by(year) %>% 
    dplyr::summarize(mean_dam   = mean(damaged_volume_total_m3  , na.rm = T),
                     sd_dam     = sd(damaged_volume_total_m3, na.rm = T),
                     mean_RS   = mean(RS_wind_beetle  , na.rm = T),
                     sd_RS     = sd(RS_wind_beetle, na.rm = T)) %>% 
  mutate(Field_volume       = stringr::str_glue("{round(mean_dam,1)}±{round(sd_dam,1)}"),
         RS_area            = stringr::str_glue("{round(mean_RS,1)}±{round(sd_RS,1)}")) %>% 
  dplyr::select(year, Field_volume, RS_area) 


(observation_mortality_year)


# Export as a nice table in word:
sjPlot::tab_df(observation_mortality_year,
               #               col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/summary_out_observation_damage_year.doc",
               digits = 0) 






#### MAP: correlations per XY: beetle counts vs damage volume --------------------
# df_cor_counts_damage <- 
#   ips_damage_clean %>% 
#   #dplyr::filter(!ID %in% c(0, 50104, 60613,62705)) %>% # exlude if there is missing data
#   group_by(trapID,forstrev_1        ) %>% 
#   dplyr::summarize(spearm_cor_lag0 = cor(damaged_volume_total_m3, sum_ips, 
#                                            method = "spearman", use = "complete.obs"),
#                    spearm_cor_lag1 = cor(damaged_volume_total_m3, sum_ips_lag1 , 
#                                             method = "spearman", use = "complete.obs"),
#                    spearm_cor_lag2 = cor(damaged_volume_total_m3, sum_ips_lag2 , 
#                                         method = "spearman", use = "complete.obs"))


#### MAP: correlations per district: beetle counts vs RS damage --------------------------
# merge based on the district ID, subset the data as the traps do not cover the whole 
# bavaria's districts

# Final table: for RS and for damaged volume ----------------------------------

### Tab for spagetti plot -----------------

df_full_years_to_spagetti_plot <-  
  ips_damage_merge %>%  
  left_join(df_all,
            by = c("forstrev_1" =   "ID",
                   "year"       = "year",
                   "damaged_volume_total_m3" = "damaged_volume_total_m3") ) %>% 
  mutate(trapID = as.factor(trapID)) %>% 
  group_by(pairID, year) %>% 
  summarise(sum_ips    = mean(sum_ips, na.rm = T),
            agg_doy    = mean(agg_doy, na.rm = T),
            peak_doy   = mean(peak_doy, na.rm = T),
            peak_diff  = mean(peak_diff, na.rm = T),
            damage_vol = mean(damaged_volume_total_m3, na.rm = T),
            RS_wind_beetle = mean(RS_wind_beetle, na.rm = T)#,
  ) %>%
  ungroup(.) %>% 
  mutate(#log_sum_ips = log(sum_ips),
         RS_wind_beetle_ha = RS_wind_beetle*0.09,  # update values for plotting
         damage_vol_k = damage_vol/1000) %>% 
#  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% c('Peiting', "Eschenbach_idOPf"))

# Replace or filter out zero values before plotting
df_full_years_to_spagetti_plot$damage_vol_k <- ifelse(df_full_years_to_spagetti_plot$damage_vol_k <= 0.01, 
                                                    0.01, 
                                                    df_full_years_to_spagetti_plot$damage_vol_k)




# keep datasets separated, as RS has lower number of years
# run analysis for trap vs RS damage
# trap vs damaged volume
df_all_lag1 <- df_all %>% 
  group_by(ID) %>% 
  arrange(year, .by_group = TRUE) %>%
  mutate(lag1_RS_wind_beetle  = lag(RS_wind_beetle , n = 1, default = NA)) #%>% 
  

# check if correct:
df_all_lag1 %>% 
  dplyr::filter(ID == 10504) %>% 
  dplyr::select(ID, year, RS_wind_beetle , lag1_RS_wind_beetle) #%>% # sum_ips, sum_ips_lag1


df_traps_RS_damage <-  
  ips_damage_clean %>%  
  left_join(df_all_lag1,
             by = c("forstrev_1" =   "ID",
                    "year"       = "year",
                    "damaged_volume_total_m3" = "damaged_volume_total_m3") ) %>% 
  mutate(trapID = as.factor(trapID)) %>% 
  group_by(pairID, year) %>% 
 summarise(sum_ips    = mean(sum_ips, na.rm = T),
           agg_doy    = mean(agg_doy, na.rm = T),
           peak_doy   = mean(peak_doy, na.rm = T),
           peak_diff  = mean(peak_diff, na.rm = T),
           damage_vol = mean(damaged_volume_total_m3, na.rm = T),
           lag1_damage_vol = mean(lag1_damaged_volume_total, na.rm = T),
           x              = mean(x, na.rm = T),
           y              = mean(y, na.rm = T),
           RS_wind_beetle = mean(RS_wind_beetle, na.rm = T),
           lag1_RS_wind_beetle = mean(lag1_RS_wind_beetle, na.rm = T),
            sum_ips_lag1  = mean(sum_ips_lag1, na.rm = T),
           sum_ips_lag2   = mean(sum_ips_lag2, na.rm = T),
           log_sum_ips_lag1 = mean(log(sum_ips_lag1), na.rm = T),
           log_sum_ips_lag2 = mean(log(sum_ips_lag2), na.rm = T),
           agg_doy_lag1     = mean(agg_doy_lag1, na.rm = T),
           agg_doy_lag2     = mean(agg_doy_lag2, na.rm = T),
           peak_doy_lag1    = mean(peak_doy_lag1, na.rm = T),
           peak_doy_lag2    = mean(peak_doy_lag2, na.rm = T),
           peak_diff_lag1   = mean(peak_diff_lag1, na.rm = T),
           peak_diff_lag2   = mean(peak_diff_lag2, na.rm = T)
           ) %>%
 ungroup(.) %>% 
 mutate(log_sum_ips = log(sum_ips)) %>% 
  mutate(f_year = as.factor(year))
  


summary(df_traps_RS_damage)

###### create final tables for RS data and volume data 

fin_dat_damage <- df_traps_RS_damage %>% 
  dplyr::select(-RS_wind_beetle) %>% 
  na.omit() 

dim(fin_dat_damage)


fin_dat_damage %>% 
  ggplot(aes(x = log(sum_ips),
             y = log(damage_vol))) +
  geom_point() +
  geom_smooth()

# do not remove outliers here


fin_dat_RS <- df_traps_RS_damage %>% 
  na.omit()

dim(fin_dat_RS)





# Find lags and predictors -------------------------------------------------------
##### Find lags: damage volume  -------------
dependent_damage <-  c("damage_vol")

# Initialize a data frame to store AIC values and deviance explained
model_metrics_damage <- data.frame(Predictor = character(), 
                                  Dependent = character(), 
                                  AIC = numeric(), 
                                  R_squared = numeric())

# list predictors to test
selected_predictors <- c('sum_ips', 'sum_ips_lag1','sum_ips_lag2',
                         #'log_sum_ips', 'log_sum_ips_lag1','sum_ips_lag2',
                         'agg_doy'   , 'agg_doy_lag1', 'agg_doy_lag2' ,
                         'peak_doy', 'peak_doy_lag1','peak_doy_lag2',
                         'peak_diff', 'peak_diff_lag1',  'peak_diff_lag2'
) 



# Loop over each dependent variable
for (dep in dependent_damage) {
  #print(dep)
  # Loop over each predictor
  for (pred in selected_predictors) {
    print(pred)
    # Fit the model
   # formula <- as.formula(paste(dep, "~ s(", pred, ", by = f_year, k = 8) + s(year, k=5) +s(x, y, bs = 'gp', k = 20)"))
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 8)"))
    
    #print(formula)
    model <- gamm(formula,  family = tw,
                  method = 'REML',  
                  data = fin_dat_damage,
                  correlation = corAR1(form = ~ year | pairID))
    
    # Extract model summary
    model_summary <- summary(model$gam)
    
    # Store the AIC value and deviance explained
    model_metrics_damage <- rbind(model_metrics_damage, data.frame(Predictor = pred, 
                                                                 Dependent = dep, 
                                                                 AIC = AIC(model$lme), 
                                                                 R_squared = round(model_summary$r.sq*100,1)
    ))
  }
}

# View the AIC values and deviance explained

(model_metrics_damage)



####### Find predicors  lags: RS   -------------
dependent_RS <-  c("RS_wind_beetle")

# Initialize a data frame to store AIC values and deviance explained
model_metrics_RS <- data.frame(Predictor = character(), 
                                   Dependent = character(), 
                                   AIC = numeric(), 
                                   R_squared = numeric())


# Loop over each dependent variable
for (dep in dependent_RS) {
  #print(dep)
  # Loop over each predictor
  for (pred in selected_predictors) {
    print(pred)
    # Fit the model
    #formula <- as.formula(paste(dep, "~ s(", pred, ", by = f_year, k = 8) + s(year, k=4) +s(x, y, bs = 'gp', k = 20)"))
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 8)"))
    #print(formula)
    model <- gamm(formula,  family = tw,
                  method = 'REML',  
                  data = fin_dat_RS,
                  correlation = corAR1(form = ~ year | pairID))
    
    # Extract model summary
    model_summary <- summary(model$gam)
    
    # Store the AIC value and deviance explained
    model_metrics_RS <- rbind(model_metrics_RS, data.frame(Predictor = pred, 
                                                                   Dependent = dep, 
                                                                   AIC = AIC(model$lme), 
                                                                   R_squared = round(model_summary$r.sq*100,1)
    ))
  }
}

# View the AIC values and deviance explained

(model_metrics_RS)

# merge two tables into one
model_lag_counts_damage_RS <- rbind(model_metrics_RS,
                           model_metrics_damage) %>% 
  distinct()


####  Export as a nice table in word: ------------------------------------------
sjPlot::tab_df(model_lag_counts_damage_RS,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_lag_predict_damage_by_counts.doc",
               digits = 1) 


# plot the best predictors and lags: ----------------------------------------

# first order the data:

average_R_squared <- model_lag_counts_damage_RS %>%
  group_by(Predictor) %>%
  summarize(avg_R_squared = mean(R_squared)) %>%
  arrange(desc(avg_R_squared))


# Merge the average R_squared back to the original data
model_lag_counts_damage_RS_avg <- model_lag_counts_damage_RS %>%
  left_join(average_R_squared, by = "Predictor") %>%
  mutate(Predictor = fct_reorder(Predictor, R_squared, .desc = TRUE))



model_lag_counts_damage_RS_avg %>% 
  ggplot(aes(x = Predictor, y = R_squared, col = Dependent)) +
  geom_segment(aes(x = Predictor, 
                   xend = Predictor, 
                   y = 0, 
                   yend = R_squared), color = "grey") +
  geom_point(size = 4) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



library(stringr)

p_lags_observed_damage <- 
  model_lag_counts_damage_RS %>% 
  dplyr::filter(!str_detect(Predictor, "log")) %>%
  mutate(Predictor = ifelse(str_detect(Predictor, "lag"), 
                            Predictor, 
                            paste0(Predictor, "_lag0"))) %>% 
  dplyr::select(-AIC) %>% 
  pivot_wider(names_from = Dependent, 
              values_from = c(R_squared)) %>% 
  arrange(-damage_vol ) %>% 
 # mutate(Predictor = case_when())
  mutate(Predictor = gsub("sum_ips", "Pop. lev.", Predictor)) %>% 
  mutate(Predictor = gsub("agg_doy", "Agg. timing", Predictor)) %>% 
  mutate(Predictor = gsub("peak_doy", "Peak s. timing", Predictor)) %>% 
  mutate(Predictor = gsub("peak_diff", "Peak s. intens.", Predictor)) %>% 
  mutate(Predictor = as.factor(Predictor)) %>%
  mutate(Predictor = fct_reorder(Predictor, -damage_vol)) %>% 
  ggplot() +
    geom_segment( aes(x=Predictor, xend=Predictor, 
                      y=damage_vol, yend=RS_wind_beetle  ), color="grey",lwd = 1.2) +
    geom_point( aes(x=Predictor, y=RS_wind_beetle, color = "RS_wind_beetle"), size=3 ) +
    geom_point( aes(x=Predictor, y=damage_vol,color = "damage_vol" ),  size=3 )  +
    scale_color_manual(values = c("RS_wind_beetle" = rgb(0.2, 0.7, 0.1, 0.5), 
                                  "damage_vol" = rgb(0.7, 0.2, 0.1, 0.5)),
                       labels = c("RS_wind_beetle" = "Remote sensing\n[# pixels]",
                                  "damage_vol" =  expression("Field observation [m"^3*"]"))) +
    labs(y = expression("Adjusted R"^2*"[%]"),
         x = "",
         color = 'Tree mortality\nsource')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = 'right')
  
windows(6,3.5)
p_lags_observed_damage

ggsave(filename = 'outFigs/observation_vs_predictors_lags.png', 
       plot = p_lags_observed_damage, width = 6, 
       height = 3.5, dpi = 300, bg = 'white')


#### Spearman correlation matrix between predictors -------------------------------------------------------

# Extract the relevant columns from the data
data_selected <- df_traps_RS_damage[selected_predictors]

cor_matrix <- cor(data_selected, use = "pairwise.complete.obs", method = "spearman")
windows()
corrplot::corrplot(cor_matrix)

# inspect relationships:
pairs(damage_vol ~  sum_ips + sum_ips_lag1 + sum_ips_lag2 +
        agg_doy    + agg_doy_lag1 +  agg_doy_lag2  +
        peak_doy +  peak_doy_lag1 + peak_doy_lag2 +
        peak_diff +  peak_diff_lag1 +   peak_diff_lag2, data = df_traps_RS_damage)
# fin_dat_damage

# check for outliers? -------------------------------------------------------------

# filter outliers from dependent and the predictors variables! 
ips_damage_clean_no_outliers <- df_traps_RS_damage %>% 
  mutate(log_damage_vol = log(damage_vol + 1)) %>% 
  dplyr::filter(log_damage_vol>5,
                log_sum_ips > 8.4) #%>% 
  






# DAMAGE VOLUME model: play with model diagnostics ---------------------------------------
# only sum ips
m1 <- gamm(damage_vol ~ #s(year,k = 4) + 
             s(log_sum_ips,k = 5)+ 
             #s(pairID, bs = 're') +
             s(x, y, bs = 'gp', k = 30),
           family = tw, #Gamma(link = "log"), # nb,
           method = 'REML',  
           data = fin_dat_damage,
           correlation = corAR1(form = ~ year | pairID))

# only years
m1.1 <- gamm(damage_vol ~ s(year,k = 4) + 
            # s(log_sum_ips,k = 5)+ 
             #s(pairID, bs = 're') +
             s(x, y, bs = 'gp', k = 30),
           family = tw, #Gamma(link = "log"), # nb,
           method = 'REML',  
           data = fin_dat_damage,
           correlation = corAR1(form = ~ year | pairID))


AIC(m1$lme, m1.1$lme)

# both years and years sums
m2 <- gamm(damage_vol ~ s(year,k = 4) + 
             s(log_sum_ips,k = 5)+ 
             #s(pairID, bs = 're') +
             s(x, y, bs = 'gp', k = 30),
           family = tw, #Gamma(link = "log"), # nb,
           method = 'REML',  
           data = fin_dat_damage,
           correlation = corAR1(form = ~ year | pairID))


# add year as factor to allow relatioship to change over years
m3 <- gamm(damage_vol ~ s(year,k = 5) + 
             s(log_sum_ips, by = f_year, k = 8)+ 
             #s(pairID, bs = 're') +
             s(x, y, bs = 'gp', k = 20),
           family = tw, #Gamma(link = "log"), # nb,
           method = 'REML',  
           data = fin_dat_damage,
           correlation = corAR1(form = ~ year | pairID))


# including both is teh best
AIC(m1$lme, m1.1$lme, m2$lme, m3$lme)

#fin.m.damage <- m3$gam



AIC(m1$lme, m2$lme, m3$lme)
appraise(m3$gam)
summary(m3$gam)
plot(m3$gam, page = 1)
gam.check(m3$gam)
k.check(m3$gam)


### PREVIOUS DAMAGE ---------------------------------------------------------
# increase the number of iterations to converge the model
control <- list(niterPQL = 50)

fin_dat_damage_no_extremes <- fin_dat_damage %>% 
  dplyr::filter(damage_vol<150000) %>% 
  dplyr::filter(lag1_damage_vol<150000) 


#### test effect of previuos year damage ----------------------------------------------
m.damage.previous <- gamm(damage_vol ~ 
                            s(sum_ips, k = 3) +
                            s(lag1_damage_vol, k = 15), # + # Added term for previous beetle counts
                           # s(pairID, bs = 're') +
                            #s(f_year, bs = 're'),
                          data = fin_dat_damage_no_extremes, 
                          family = tw,
                          control = control)

m.damage.previous2 <- gamm(damage_vol ~ 
                            s(sum_ips, k = 3) +
                            s(lag1_damage_vol, k = 15) + # Added term for previous beetle counts
                           s(pairID, bs = 're'),# +
                          #s(f_year, bs = 're'),
                          data = fin_dat_damage_no_extremes, 
                          family = tw,
                          control = control)

m.damage.previous3 <- gamm(damage_vol ~ 
                             s(sum_ips, k = 3) +
                             s(lag1_damage_vol, k = 4) + # Added term for previous beetle counts
                             s(pairID, bs = 're') +
                             #s(x, y, bs = 'tp', k = 5),
                           s(f_year, bs = 're'),
                           data = fin_dat_damage_no_extremes, 
                           family = tw,
                           control = control)


AIC(m.damage.previous2, m.damage.previous, m.damage.previous3)

fin.m.damage <- m.damage.previous3$gam

appraise(m.damage.previous3$gam)
summary(m.damage.previous$gam)
plot(m.damage.previous3$gam, page = 1)
gam.check(m.damage.previous3$gam)
k.check(m.damage.previous3$gam)
draw(m.damage.previous3$gam)

# check for autocorrelation
# Extract residuals from the model
residuals <- resid(m.damage.previous3$lme, type = "normalized")

# Plot ACF of the residuals
acf(residuals, main="ACF of Model Residuals")
pacf(residuals, main="ACF of Model Residuals")





# RS model: play with model diagnostics ---------------------------------------
# use teh most important predictor: sum_ips_lag1
boxplot(log(fin_dat_RS$RS_wind_beetle + 0.01))
boxplot(fin_dat_RS$RS_wind_beetle)
boxplot(log(fin_dat_damage$damage_vol))
hist(fin_dat_RS$RS_wind_beetle)
# if using log values, the data do not seems as outliers

plot(x = sum_ips_lag1, y = RS_wind_beetle, data = fin_dat_RS)

fin_dat_RS %>% 
  ggplot(aes(x = sum_ips_lag1,
             y = RS_wind_beetle)) + 
  geom_point()+
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10()


fin_dat_damage %>% 
  ggplot(aes(x = sum_ips,
             y = damage_vol)) + 
  geom_point()+
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10()


fin_dat_RS_clean <- fin_dat_RS %>% 
  dplyr::filter(RS_wind_beetle <1700) %>% 
  dplyr::filter(lag1_RS_wind_beetle <1700)


# PREVIOUS RS --------------------------------------

mRS.previous <- gamm(RS_wind_beetle ~ #s(year,k = 4) + 
             s(lag1_RS_wind_beetle, k = 9) + 
             s(sum_ips_lag1 , k=10), #+ 
             # #s(pairID, bs = 're') +
             #s(x, y, bs = 'tp', k = 50),
           family = tw, 
           method = 'REML',  
           data = fin_dat_RS_clean)


mRS.previous1 <- gamm(RS_wind_beetle ~ #s(year,k = 4) + 
                       s(lag1_RS_wind_beetle, k = 6) + 
                       s(sum_ips_lag1 , k=5)+ 
                       s(pairID, bs = 're'),# +
                     #s(x, y, bs = 'tp', k = 50),
                     family = tw, 
                     method = 'REML',  
                     data = fin_dat_RS_clean)


mRS.previous2 <- gamm(RS_wind_beetle ~ #s(year,k = 4) + 
                        s(lag1_RS_wind_beetle, k = 6) + 
                        s(sum_ips_lag1 , k=5)+ 
                        s(pairID, bs = 're') +
                       s(x, y, bs = 'tp', k = 5),
                      family = tw, 
                      method = 'REML',  
                      data = fin_dat_RS_clean)


AIC(mRS.previous1, mRS.previous, mRS.previous2)

appraise(mRS.previous1$gam)
summary(mRS.previous1$gam)
plot(mRS.previous1$gam, page = 1)
gam.check(mRS.previous1$gam)
k.check(mRS.previous1$gam)
draw(mRS.previous1$gam)

# check for autocorrelation
# Extract residuals from the model
residuals <- resid(mRS.previous1$lme, type = "normalized")

# Plot ACF of the residuals
acf(residuals, main="ACF of Model Residuals")
pacf(residuals, main="ACF of Model Residuals")









# the best model!
fin.m.RS <- mRS.previous1$gam

sjPlot::tab_model(fin.m.RS,     file = "outTable/model_trap_RS.doc")
sjPlot::tab_model(fin.m.damage, file = "outTable/model_trap_damage.doc")



# prepare two plots as output: 
# merge together the predicted data by year, log sum ips for RS, NFI
# export in a single plot
# use log scales 



# Effect plots ----------------------------------------------------
summary(fin.m.damage)
summary(fin.m.RS)


# Spagetti plots for years -----------------------------------------------------

library(scales)
# update data to scatter plot 
fin_dat_damage_plot <- fin_dat_damage %>% 
  dplyr::filter(lag1_damage_vol < 150000) %>% 
  dplyr::filter(damage_vol < 150000) %>% 
  mutate(damage_vol       = damage_vol/1000,
         sum_ips          = sum_ips/1000,
         lag1_damage_vol  = lag1_damage_vol/1000)


fin_dat_RS_clean_plot <- fin_dat_RS_clean %>% 
  mutate(RS_wind_beetle_ha      = RS_wind_beetle*0.09,
         lag1_RS_wind_beetle_ha = lag1_RS_wind_beetle*0.09,
         sum_ips_lag1           = sum_ips_lag1/1000)




p_spagett_damage <- plot_data_with_average(df_full_years_to_spagetti_plot, "damage_vol_k", 'my_lab',   
                                              my_title = expression("Ground survey [*1000 m"^3*"]"))+ 
  #scale_y_log10() +
  scale_y_log10(labels = label_number()) +  # Adjust y-axis labels to avoid scientific notation
  lims(x = c(2015,2021))

(p_spagett_damage)
range(df_full_years_to_spagetti_plot$damage_vol_k)
hist(df_full_years_to_spagetti_plot$damage_vol_k)


p_spagett_RS<- plot_data_with_average(df_full_years_to_spagetti_plot, "RS_wind_beetle_ha", 'my_lab',   
                                       my_title = paste('Remote sensing', '[ha]')) +  
                                         scale_y_log10() + 
  lims(x = c(2015,2021))

(p_spagett_RS)

# PLOT: field based damage --------------------------

# Assuming 'model' is your glm.nb model
summary(fin.m.damage)
p0 <- ggpredict(fin.m.damage, terms = "lag1_damage_vol [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.damage, terms = "sum_ips [all]", allow.new.levels = TRUE)

p0$x <- p0$x/1000 # convert to original values
p0$predicted <- p0$predicted/1000
p0$conf.low <- p0$conf.low /1000
p0$conf.high <- p0$conf.high /1000

p1$x <- p1$x/1000 # convert to original values
p1$predicted <- p1$predicted/1000
p1$conf.low <- p1$conf.low /1000
p1$conf.high <- p1$conf.high /1000


p0.damage <- 
  create_effect_plot(data = p0, 
                     avg_data = fin_dat_damage_plot, 
                     x_col = "lag1_damage_vol", 
                     y_col = "damage_vol", 
                     line_color = "darkgreen", 
                     x_title = expression('Tree mortality lag1' ~ "[" * m^3 * " * 1000]"), 
                     #y_title = paste(lab_popul_level, '*1000'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 50, 
                     lab_annotate = "***")

(p0.damage)




p1.damage <- 
  create_effect_plot(data = p1, 
                     avg_data = fin_dat_damage_plot, 
                     x_col = "sum_ips", 
                     y_col = "damage_vol", 
                     line_color = "grey30", 
                     x_title = 'Population level [#*1000]' , 
                     y_title = paste(lab_popul_level, '*1000'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 50, 
                     lab_annotate = "***")

(p1.damage)

##### PLOT RS based damage -----------------------------------------

# Assuming 'model' is your glm.nb model
summary(fin.m.RS)
p0 <- ggpredict(fin.m.RS, terms = "lag1_RS_wind_beetle [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.RS, terms = "sum_ips_lag1 [all]", allow.new.levels = TRUE)

p0$x <- p0$x*0.09 # convert piel counts to hectares
p0$predicted <- p0$predicted*0.09 #/100
p0$conf.low <- p0$conf.low*0.09 #/100
p0$conf.high <- p0$conf.high*0.09# /100

p1$x <- p1$x/1000 # convert to original values
p1$predicted <- p1$predicted*0.09
p1$conf.low <- p1$conf.low *0.09
p1$conf.high <- p1$conf.high *0.09


p0.RS <- 
  create_effect_plot(data = p0, 
                     avg_data = fin_dat_RS_clean_plot, 
                     x_col = "lag1_RS_wind_beetle_ha", 
                     y_col = "RS_wind_beetle_ha", 
                     line_color = "darkgreen", 
                     x_title = 'RS Tree mortality lag1 [ha]' , 
                     #y_title = paste(lab_popul_level, '*1000'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 75, 
                     lab_annotate = "***")

(p0.RS)




p1.RS <- 
  create_effect_plot(data = p1, 
                     avg_data = fin_dat_RS_clean_plot, 
                     x_col = "sum_ips_lag1", 
                     y_col = "RS_wind_beetle_ha", 
                     line_color = "grey30", 
                     x_title = 'Population level lag1 [#*1000]' , 
                     #y_title = paste(lab_popul_level, '*1000'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 40, 
                     lab_annotate = "0.08")

(p1.RS)

p.comp2 <- ggarrange(
  p_spagett_damage, p0.damage, p1.damage, 
  p_spagett_RS, p0.RS, p1.RS, 
  ncol = 3, nrow = 2, align = 'hv',
  labels = c('[a]', '[b]', '[c]', '[d]', '[e]', '[f]'),
  font.label = list(size = 10, face = "plain"),  # Set font to normal (plain) and specify size
  label.x = 0.15,  # Adjust these values to position the labels within the plot
  label.y = 0.85
)
(p.comp2)


ggsave(filename = 'outFigs/observation_vs_traps.png', 
       plot = p.comp2, width = 6, height = 4, dpi = 300, bg = 'white')


##### quick plotting -----------------------------------------------

# to exclude  RS wind beetle in Passau 2018??? no!

# create area plots: sum damage per year
df_sum <- df_all %>% 
  group_by(year) %>%
  dplyr::summarise(RS_harvest = sum(RS_harvest, na.rm = T),
                   RS_wind_beetle  = sum(RS_wind_beetle , na.rm = T),
                   damaged_volume_total_m3  = sum(damaged_volume_total_m3,na.rm = T ))
  

p1 <- df_sum %>% 
  ggplot(aes(x = year,
             y = RS_wind_beetle*0.09)) +
  geom_line(col = 'green', lty = 'dashed', lwd = 1.5) + 
  labs(y = 'RS tree mortality [ha]') + 
  theme_classic()


p2<- df_sum %>% 
  ggplot(aes(x = year,
             y = damaged_volume_total_m3/10000)) +
  geom_line(col = 'red', lwd = 1.5) +
  labs(y = 'Reported tree damage [*10000 m^3]') + 
  theme_classic()

# df_sum
windows(7,3)
ggarrange(p1, p2, nrow = 1, ncol = 2, align = 'hv', labels = c('[a]', '[b]'))





##### merge damage data to geometry data ----------------------------------

# Correlate between RS and volume data m3 ---------------------------------------------------
df_cor <- df_all %>% 
  dplyr::filter(!ID %in% c(0, 50104, 60613,62705)) %>% # exlude if there is missing data
  group_by(ID) %>% 
  dplyr::summarize(spearm_cor_beetle = cor(damaged_volume_total_m3, RS_wind_beetle, 
                                           method = "spearman", use = "complete.obs"))


summary(df_cor)

df_all %>% 
  dplyr::filter(RS_wind_beetle < 10000) %>% 
  ggplot(aes(x = RS_wind_beetle, y = damaged_volume_total_m3)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_cor(method = "pearson", label.x = 3000, label.y = max(df_all$damaged_volume_total_m3, na.rm = TRUE))


cor_shp <- sf_simpled %>% 
  left_join(df_cor, by = c("forstrev_1" = "ID"))

all_shp <- sf_simpled %>% 
  left_join(df_all, by = c("forstrev_1" = "ID"))

df_damage_sum <- df_all %>% 
  group_by(ID) %>% 
  summarise(sum_dmg  =sum(damaged_volume_total_m3, na.rm = T),
            sum_RS  =sum(RS_wind_beetle, na.rm = T)) %>% 
  dplyr::filter(ID!=0) %>% 
  ungroup(.)

# add data to the full geometry to run teh intersection with the trap data
#all_shp_complex <- sf_districts_trans %>% 
#  left_join(df_all, by = c("forstrev_1" = "ID"))

dmg_shp <- 
  sf_simpled %>% 
    dplyr::select(c(forstrev_1)) %>% 
  full_join(df_damage_sum, by = c("forstrev_1" = "ID")) %>% 
  dplyr::select(c(forstrev_1, sum_dmg,sum_RS)) %>% 
  mutate(log_sum_dmg )
  
  

# cap values at 2e05

hist(dmg_shp$sum_dmg)
hist(dmg_shp$sum_RS)

 

 # Calculate the XX percentile of the sum_dmg column
 perc95_dmg <- quantile(dmg_shp$sum_dmg, 0.95)
 perc95_RS <- quantile(dmg_shp$sum_RS, 0.95)
 


# cap values at 95% interval
dmg_shp <- dmg_shp %>% 
  mutate(dmg_cap = case_when(sum_dmg > perc95_dmg~perc95_dmg, TRUE ~ sum_dmg),
         RS_cap = case_when(sum_RS > perc95_RS ~ perc95_RS, TRUE ~ sum_RS))

 #sum coundt of beetles by locations 

df_sum_ips_years <- sum_ips %>% 
   mutate(pairID = str_sub(falsto_name, 1, -3)) %>% 
  group_by(pairID) %>% 
  summarise(sum = sum(sum_ips, na.rm = T)/2)  # divide by two to have average sum per pair

 # merge with simpler trap pairs (sites)
 df_sum_ips_years_sf <- xy_sf_fin %>% 
   mutate(pairID = str_sub(falsto_name, 1, -3)) %>% 
   group_by(pairID) %>% 
   slice(1) %>%  # Keep only the on trap per pair
   ungroup(.) %>% 
   full_join(df_sum_ips_years,by = join_by(pairID)) %>% 
   dplyr::select(pairID, sum) #%>% 
   
 
 # Calculate the XX percentile of the sum_dmg column
 perc95_sum_ips <- quantile(df_sum_ips_years_sf$sum, 0.95)
 
 df_sum_ips_years_sf <- df_sum_ips_years_sf %>% 
   mutate(sum_cap = case_when(sum > perc95_sum_ips~perc95_sum_ips, TRUE ~ sum))
 
 
 hist(df_sum_ips_years_sf$sum)
 
# MAPS -----------------------------------------------------------------------------------
 bav_shp <- st_read('rawData/outline_bavaria.gpkg')
 bav_3035 <- st_transform(bav_shp, crs(dmg_shp))
 
 #4298150 
 
 
# ad damage data to district geometry
# MAP: plot correlation with beetle counts & tree damage ----------------------------------
# Extract colors from RColorBrewer's YlOrRd palette
colors <- rev(brewer.pal(9,  "RdYlBu" )) # "YlOrRd"  "RdYlBu"


map.damage <- 
  ggplot(dmg_shp) +
   geom_sf(data = dmg_shp,
          aes(fill = dmg_cap/1000)) +
    geom_sf(data = bav_3035,
            color = 'black', fill = NA, lwd = 0.8)  +
    geom_sf(data = df_sum_ips_years_sf,
            aes(color = sum_cap/1000), color = 'grey20', size = 2.5) +
    geom_sf(data = df_sum_ips_years_sf,
            aes(color = sum_cap/1000), size = 2) +
  scale_color_gradientn(
    name = "Beetle catch [#*1000]",
    colors = colors,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )  +
    scale_fill_gradientn(
    name = bquote(Tree~mortality~('1000'~m^3)), #expression(Tree~mortality~('1000'm^3)) , #"Tree mortality [m3]",
    colors = colors,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )  +
  labs(title = "[a]") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(angle = 0, size = 8),  # Smaller title
    legend.text = element_text(size =8),              # Smaller text
    legend.key.width = unit(0.5, "cm"),                  # Adjust width
    legend.key.height = unit(0.2, "cm"),               # Adjust height
    plot.title = element_text(hjust = 0, size = 8 )             # Center title
  ) #+
  #annotation_scale(location = "bl", width_hint = 0.5) +  # Add scale bar at bottom left
  #annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering)  # Add north arrow at top left
map.damage



map.RS <- 
  ggplot(dmg_shp) +
   geom_sf(data = dmg_shp,
          aes(fill = RS_cap*0.09)) +
  geom_sf(data = bav_3035,
          color = 'black', fill = NA, lwd = 0.8)  +
  geom_sf(data = df_sum_ips_years_sf,
          aes(color = sum_cap/1000), color = 'grey20', size = 2.5) +
  geom_sf(data = df_sum_ips_years_sf,
          aes(color = sum_cap/1000), size = 2) +
  scale_color_gradientn(
    name = "Beetle catch [#*1000]",
    colors = colors,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )  +
  scale_fill_gradientn(
    name = "Tree mortality [ha]",
    colors = colors,
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )  +
  labs(title = "[b]") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(angle = 0, size = 8),  # Smaller title
    legend.text = element_text(size = 8),              # Smaller text
    legend.key.width = unit(0.5, "cm"),                  # Adjust width
    legend.key.height = unit(0.2, "cm") ,               # Adjust height
    plot.title = element_text(hjust = 0, size = 8)             # Center title
  )

windows(7,5)
damage_maps <- ggarrange(map.damage, map.RS, align = 'hv')
damage_maps

ggsave(filename = 'outFigs/damage_maps.png', 
       plot = damage_maps, width = 7, 
       height = 5, dpi = 300, bg = 'white')



# MAP: plot correlation with beetle counts & tree damage ----------------------------------
ggplot(cor_shp) +
   geom_sf(color = 'grey93', 
            fill  = 'grey93') + 
  geom_sf(data = cor_shp_counts_dmg,
          aes(fill = value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,  # Set the midpoint at 0 for white color
                       limit = c(min(cor_shp_counts_dmg$value, na.rm = TRUE), 
                                 max(cor_shp_counts_dmg$value, na.rm = TRUE)),
                       name = "Spearman\nCorrelation") +
  labs(title = "Spearman Correlation Coefficient",
       subtitle = "between spruce damage [m3] and beetle sum/year by region") +
  theme_void()  + 
  facet_wrap(.~name)#





# Plots & maps ------------------------------------------------------------


### RS plots -----------------------------------------------------------------


df_damage %>%
  ggplot(aes(x = Year, 
             y = damaged_volume_total_m3/10000, 
             group = Year,
             fill = factor(Year))) + 
  # geom_violin() +
  geom_boxplot()


df_damage %>%
  ggplot(aes(x = Year, 
             y = damaged_volume_total_m3)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))




#### show spearman correlation coefficient  -------
# between volume damage and RS observation (wind, beetle)

# wind_beetle vs damage volume 

p1 <- ggplot(cor_shp) +
  # geom_sf(color = 'black', 
  #          fill  = 'grey93') #+ 
  geom_sf(data = cor_shp,
          aes(fill = spearm_cor_beetle)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,  # Set the midpoint at 0 for white color
                       limit = c(min(cor_shp$spearm_cor_beetle, na.rm = TRUE), 
                                 max(cor_shp$spearm_cor_beetle, na.rm = TRUE)),
                       name = "Spearman\nCorrelation") +
  labs(title = "Spearman Correlation Coefficient",
       subtitle = "Correlation between ground surevy data and RS based estimation [ha] by forest district") +
  theme_minimal()  #


summary(cor_shp$spearm_cor_beetle)
mean(cor_shp$spearm_cor_beetle, na.rm = T)
sd(cor_shp$spearm_cor_beetle, na.rm = T)
median(cor_shp$spearm_cor_beetle, na.rm = T)


# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.9856  0.1518  0.4928  0.4412  0.8117  1.0000       3 

# RS harvest vs damage volume 

p2 <- ggplot(cor_shp) +
  # geom_sf(color = 'black', 
  #          fill  = 'grey93') #+ 
  geom_sf(data = cor_shp,
          aes(fill = spearm_cor_harvest)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,  # Set the midpoint at 0 for white color
                       limit = c(min(cor_shp$spearm_cor_harvest, na.rm = TRUE), 
                                 max(cor_shp$spearm_cor_harvest, na.rm = TRUE)),
                       name = "Spearman\nCorrelation") +
  labs(title = "Spearman Correlation Coefficient",
       subtitle = "Correlation between spruce damage [m3] and RS harvest [pixel counts] by region") +
  theme_minimal()  #



p3 <- ggplot(cor_shp) +
  # geom_sf(color = 'black', 
  #          fill  = 'grey93') #+ 
  geom_sf(data = cor_shp,
          aes(fill = spearm_cor_sum )) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,  # Set the midpoint at 0 for white color
                       limit = c(min(cor_shp$spearm_cor_sum, na.rm = TRUE), 
                                 max(cor_shp$spearm_cor_sum, na.rm = TRUE)),
                       name = "Spearman\nCorrelation") +
  labs(title = "Spearman Correlation Coefficient",
       subtitle = "Correlation between spruce damage [m3] and RS sum damage [pixel counts] by region") +
  theme_minimal()  #


ggarrange(p1, p2, p3)




### Map damage [m3] per year ------------------------------------------------

summary(all_shp$damaged_volume_total_m3)

# Define the 75th percentile
damage_75 <- quantile(all_shp$damaged_volume_total_m3, probs = 0.75)
RS_wind_beetle_75 <- quantile(all_shp$RS_wind_beetle , probs = 0.75, na.rm = T)

# cap tthe values for better visualization of data
all_shp <- all_shp %>% 
  mutate(damaged_volume_cap = case_when(damaged_volume_total_m3 >= damage_75 ~ damage_75,
                                        TRUE ~ damaged_volume_total_m3),
         RS_wind_beetle_cap = case_when(RS_wind_beetle >= RS_wind_beetle_75 ~ RS_wind_beetle_75,
                                        TRUE ~ RS_wind_beetle))


# Plot the map
p1<-ggplot(all_shp) +
  geom_sf(aes(fill = damaged_volume_cap)) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      name = expression("Tree damage [m"^3*"]"),  
                      na.value = 'transparent') +
  facet_wrap(. ~ year) +
  theme(legend.position = 'bottom') +
  theme_void() 


p2<-ggplot(all_shp) +
  geom_sf(aes(fill = RS_wind_beetle_cap*0.09)) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      name = expression("RS damage [ha]"),  
                      na.value = 'transparent') +
  facet_wrap(. ~ year) +
  theme(legend.position = 'bottom') +
  theme_void() 


ggarrange(p1, p2, labels = c('[a] Damaged volume', '[b] RS Damaged area'))

# plot gam between damage and RS beetle ------------------------------------------
df_all %>% 
  ggplot(aes(x = log(RS_wind_beetle +1) ,
             y = log(damaged_volume_total_m3) )) +
  geom_point() +
  geom_smooth() + 
  facet_wrap(.~year, scales = 'free')




