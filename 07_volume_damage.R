## -----------------------------------------------------
#
# Process volume data : tree damage from LWF
#
# -----------------------------------------------------

# Goal: do beetle counts (+ lags) predict tree damage?  

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

#library(knitr)   # for table outputs



# Data --------------------------------------------------------------------------
### read beetle data 
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
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)) 

# convert to sf object, crs 3035
ips_sum_sf <-  xy_sf_expand %>%  left_join(ips_sum, by = c("falsto_name", 'year')) 
#st_is_valid(ips_sum_sf_trans)

# change projection & convert to terra
my_crs = st_crs(disturb_year)

# change projection & convert to terra
# districts = polygons
sf_districts_trans <- st_transform(sf_districts, crs = my_crs) # st_crs(disturb_year))
v_districts        <-  terra::vect(sf_districts_trans) # convert to terra

ips_sum_sf_trans  <- st_transform(ips_sum_sf, crs = my_crs)
v_ips_sum         <- vect(ips_sum_sf_trans)


# works well for terra! 
plot(forest)
plot(disturb_type)
plot(v_ips_sum, add = T)



# simplify terra object:
v_districts_simpl <- simplifyGeom(v_districts, tolerance=20, preserveTopology=TRUE, makeValid=TRUE)
v_ips_sum_simpl   <- simplifyGeom(v_ips_sum, tolerance=20, preserveTopology=TRUE, makeValid=TRUE)


writeVector(v_ips_sum_simpl, 'outSpatial/XY_beetle_sum_years.gpkg', overwrite=TRUE) 

#  does simplification work? yes!!
plot(v_districts_simpl)
plot(v_ips_sum_simpl, add = T, col = 'red')

identical(crs(v_districts_simpl), crs(v_ips_sum_simpl))


# Process rasters ------------------------------------------------------------
# there is an aditiona row in the damage data with district ID 0 - I can remove that


###### crop the rasters: they are for whole germany, crop to Bavaria --------------------------
dist_year_crop <- terra::crop(disturb_year, v_districts_simpl)
dist_type_crop <- terra::crop(disturb_type, v_districts_simpl)
forest_crop    <- terra::crop(forest,       v_districts_simpl)

# process spatial data: check geometries, change projection to Corneliuses data, convert to raster the AELF districts
crs(disturb_type) <- crs(dist_year_crop)
crs(forest_crop)  <- crs(dist_year_crop)


# Simplify the df districts geometry -------------------------------
sf_simpled <- st_as_sf(v_districts_simpl) #st_simplify(sf_districts_trans, dTolerance = 20, preserveTopology = TRUE)

identical(crs(sf_simpled), crs(v_districts_simpl))

###### convert districts poly to raster ------------------------------------

# get extend and set correct projectipn
ext <- as(extent(sf_simpled), 'SpatialPolygons')  # treat raster as from raster, not terra package
ext <- st_as_sf(ext)

identical(crs(ips_sum_sf_trans), crs(ext))

# efine teh crs first - it was NA
st_crs(ext) <- my_crs# st_crs(raster(disturb_year))

#ext  <- st_transform(ext, crs = my_crs)

st_crs(sf_simpled) == st_crs(ext) 
st_crs(disturb_year) == st_crs(ext) 
st_crs(ips_sum_sf_trans) == st_crs(ext) 

# create a base raster to rasteroize to it
base_raster = raster(dist_year_crop)

# Create grid index for each pixel
grid_sel     <- st_intersection(sf_simpled, ext)
grid_sel_ras <- fasterize::fasterize(grid_sel,
                                     base_raster,
                                     #raster(dist_year_crop), 
                                     field = "forstrev_1") # name the new rasters as polygon ids  #, field = "forstrev_1"


#ips_sum_sf_trans2  <- st_transform(ips_sum_sf_trans, crs = my_crs)
# try of geometry fits - does not fit!!!
plot(grid_sel_ras)
plot(ips_sum_sf_trans ,add = T)

plot(rast(grid_sel_ras))
plot(vect(ips_sum_sf_trans) ,add = T)

# assigne the crs manually! differet assigenement for vector and raster
crs(grid_sel_ras)              <- crs(sf_simpled)
st_crs(ips_sum_sf_trans)       <- crs(sf_simpled)


#ips_sum_sf_trans <- st_transform(ips_sum_sf_trans, crs_raster@projargs)


# Extract CRS again after transformation
crs_raster <- crs(grid_sel_ras)
crs_vector <- st_crs(ips_sum_sf_trans)

# Confirm they are now the same
identical(as.character(crs_raster), crs_vector$wkt)  # Should return TRUE


# check if the crs fits! YES
st_crs(grid_sel_ras)       == st_crs(sf_simpled)
st_crs(grid_sel_ras)       == st_crs(ips_sum_sf_trans)
#st_crs(rast(grid_sel_ras)) == st_crs(vect(ips_sum_sf_trans))
#crs(rast(grid_sel_ras))    == crs(vect(ips_sum_sf_trans))


# extract the 'reviers ID' (from the districts) to the beetle counts data 
district_ID <- terra::extract(rast(grid_sel_ras),
                              v_ips_sum_simpl,
                            #  vect(ips_sum_sf_trans), 
                              df = TRUE, bind = TRUE )

#### Link beetle counts vs damage volume per department (sub district)  -------------------------------------------------------
ips_damage_district <- district_ID %>% 
  as.data.frame() %>% 
  rename(forstrev_1 = layer)


# join trap counts with damage data
ips_damage_merge <- ips_damage_district %>% 
  left_join(df_damage, by = c('forstrev_1' = 'AELF_district_ID',
                              'year' = 'Year')) %>% 
  dplyr::select(-c(id, AELF_district_name,AELF_name, AELF_ID,
                 pest_species, tree_species, unit,
                 globalid )) %>% 
  # lag tree damage: does beetle counts (or their lags) predict tree damage?:
  group_by(falsto_name) %>% 
  arrange(year, .by_group = TRUE) %>%
  mutate(sum_beetle_lag1  = lag(sum_beetle, n = 1, default = NA),
         sum_beetle_lag2  = lag(sum_beetle, n = 2, default = NA, order_by = year))

# remove NAs
ips_damage_clean <- na.omit(ips_damage_merge)

# get log values - makes it nice for damaged volume
ips_damage_clean$log_sum_beetle              <- log(ips_damage_clean$sum_beetle + 1)
ips_damage_clean$log_sum_beetle_lag1         <- log(ips_damage_clean$sum_beetle_lag1 + 1)
ips_damage_clean$log_sum_beetle_lag2         <- log(ips_damage_clean$sum_beetle_lag2 + 1)
ips_damage_clean$log_damaged_volume_total_m3 <- log(ips_damage_clean$damaged_volume_total_m3 + 1)


# try correlation
cor(ips_damage_clean$sum_beetle, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')
cor(ips_damage_clean$sum_beetle_lag1, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')
cor(ips_damage_clean$sum_beetle_lag2, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')



# try correlation: of the log values: the same values!
cor(ips_damage_clean$log_sum_beetle, ips_damage_clean$log_damaged_volume_total_m3)
cor(ips_damage_clean$log_sum_beetle_lag1, ips_damage_clean$log_damaged_volume_total_m3)
cor(ips_damage_clean$log_sum_beetle_lag2, ips_damage_clean$log_damaged_volume_total_m3)

#library(glm)
m_lag0 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle, data = ips_damage_clean, family = gaussian())
m_lag1 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle_lag1, data = ips_damage_clean, family = gaussian())
m_lag2 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle_lag2, data = ips_damage_clean, family = gaussian())

m_poly <- glm(formula = log_damaged_volume_total_m3 ~ log_sum_beetle + I(log_sum_beetle^2), family = gaussian(), data = ips_damage_clean)
# the lag0 is teh best from AIC, also teh variance explained: 6% for lag0
AIC(m_lag0,m_lag1,m_lag2,m_poly)
BIC(m_lag0,m_lag1,m_lag2)
#r2(m_lag0,m_lag1,m_lag2)
summary(m_poly)
simulateResiduals(m_poly, plot = T)
plot(allEffects(m_poly))

##### gam ---------------
library(mgcv)
gam_model_lag0 <- gam(log_damaged_volume_total_m3 ~ s(log_sum_beetle), family = gaussian(), data = ips_damage_clean)
gam_model_lag1 <- gam(log_damaged_volume_total_m3 ~ s(log_sum_beetle_lag1), family = gaussian(), data = ips_damage_clean)
gam_model_lag2 <- gam(log_damaged_volume_total_m3 ~ s(log_sum_beetle_lag2), family = gaussian(), data = ips_damage_clean)

plot(gam_model_lag0)

AIC(gam_model_lag0,gam_model_lag1, gam_model_lag2, m_lag0,m_lag1,m_lag2,m_poly, gam_model_lag0_6) # ,gam_model_0_6

# test k size
gam_model_lag0_3 <- gam(log_damaged_volume_total_m3 ~ s(log_sum_beetle, k =3), family = gaussian(), data = ips_damage_clean)
gam_model_lag0_6 <- gam(log_damaged_volume_total_m3 ~ s(log_sum_beetle, k =4), family = gaussian(), data = ips_damage_clean)

# try without logs - wrong! better to stick with 
gam_model_0_6 <- gam(damaged_volume_total_m3 ~ s(sum_beetle, k =6), family = tw(), data = ips_damage_clean)

plot(gam_model_lag0_3)
plot(gam_model_lag0_6)
plot(gam_model_0_6)
summary(gam_model_0_6)

AIC(gam_model_0_6,gam_model_lag0_6 )
summary(gam_model_lag0_6) # the best!!!

# include random effects
gam_tw1 <- gam(damaged_volume_total_m3 ~ s(sum_beetle, k = 3) + s(year,
                                                                               k = 3) , family = tw(),method = "REML", data = ips_damage_clean)
AIC(gam_model_0_6, gam_tw1)
plot(gamm_fit, page = 1)


# !!! try the account for autocorrelation -------------
library(nlme)

# Fit a GAMM model with AR(1) autocorrelation structure
gamm_fit <- gam(damaged_volume_total_m3 ~ s(sum_beetle, k = 3) + s(year, k = 5),
                 data = ips_damage_clean,
                # random = list(forstrev_1  = ~1),
                # correlation = corAR1(form = ~ year | forstrev_1 ),  # AR(1) within each location
                 method = "REML",
                 family = tw())

# Summary of the fit
summary(gamm_fit$gam)
summary(gamm_fit$lme)  # Check the output of the mixed effects component, especially the correlation structure

# plot using ggeffects:
p1 <- ggpredict(gam_model_lag0_6, terms = "log_sum_beetle [all]", allow.new.levels = TRUE)

ggplot(p1, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_line(color = 'blue') +
  geom_ribbon(alpha = 0.1, fill = "blue") +
  labs(y = "Effect on Log Damaged Volume (m^3)", x = "Log Sum Beetle") +
  theme_classic()
#pred_effects <- ggpredict(gam_model_lag0_6, terms = "log_sum_beetle")

#p <- plot(pred_effects) + 
#  labs(y = "Effect on Log Damaged Volume (m^3)", x = "Log Sum Beetle") +
#  ggtitle("Predicted Effects of Log Sum Beetle on Tree Damage Volume")
#print(p)

# Function to convert log values to original scale and interpret
interpret_effect <- function(log_beetle) {
  beetle_count <- exp(log_beetle)
  damaged_volume <- exp(predict(gam_model_lag0_6, newdata = data.frame(log_sum_beetle = log_beetle)))
  cat("With a beetle count of", round(beetle_count), "we predict a tree damage volume of approximately", round(damaged_volume), "cubic meters.\n")
}

# Example usage
interpret_effect(8)


library(lme4)
mixed_model <- lmer(log_damaged_volume_total_m3 ~ log_sum_beetle + (1 | falsto_name  ) + (1 | year), data = ips_damage_clean)
mixed_model_poly <- lmer(log_damaged_volume_total_m3 ~ poly(log_sum_beetle,2) + (1 | falsto_name  ) + (1 | year), data = ips_damage_clean)

#mixed_model_poly <- lmer(log_damaged_volume_total_m3 ~ log_sum_beetle*falsto_name + (1 | falsto_name  ), data = ips_damage_clean)


summary(mixed_model_poly)
r2(mixed_model_poly)
plot(mixed_model)
plot(allEffects(mixed_model_poly))
AIC(mixed_model,gam_model_lag0_6,mixed_model_poly)

# try glm again with teh factor
m_level_lag0 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle +  pop_level, data = ips_damage_clean, family = gaussian())
#m_level_poly <- glm(log_damaged_volume_total_m3 ~ poly(log_sum_beetle,2) +  pop_level, data = ips_damage_clean, family = gaussian())
m_level_poly2 <- glm(log_damaged_volume_total_m3 ~ poly(log_sum_beetle,2), data = ips_damage_clean, family = gaussian())


summary(m_level_lag0)
plot(allEffects(m_level_poly2))

AIC(m_level_poly,m_level_lag0,m_level_poly2)

boxplot(ips_damage_clean$log_sum_beetle)
# check relationship by years
ips_damage_clean %>% 
  ggplot(aes(x = log_sum_beetle,
             y = log_damaged_volume_total_m3)) +
  geom_point() +
  geom_smooth() #+
  #facet_wrap(.~year)


# use a gamma to predictthe tree damage (continuous,positive)
# Fitting a Gamma model
gamma_m_lag0 <- glm(damaged_volume_total_m3 ~ sum_beetle, family = Gamma(link = "log"), data = ips_damage_clean)
gamma_m_lag1 <- glm(damaged_volume_total_m3 ~ sum_beetle_lag1, family = Gamma(link = "log"), data = ips_damage_clean)
gamma_m_lag2 <- glm(damaged_volume_total_m3 ~ sum_beetle_lag2, family = Gamma(link = "log"), data = ips_damage_clean)

AIC(gamma_m_lag0, gamma_m_lag1, gamma_m_lag2)
simulateResiduals(gamma_m_lag0, plot = T)

# r2(gamma_m_lag0, gamma_m_lag1, gamma_m_lag2)

# try tweedie???
ips_damage_clean$sum_beetle_scaled <- scale(ips_damage_clean$sum_beetle)
ips_damage_clean$sum_beetle_lag1_scaled <- scale(ips_damage_clean$sum_beetle_lag1)
ips_damage_clean$sum_beetle_lag2_scaled <- scale(ips_damage_clean$sum_beetle_lag2)

# run without scaling: running into convergence issues; too long runnning time
#' tw_lag0 <- glmmTMB(damaged_volume_total_m3 ~ sum_beetle_scaled, family = tweedie(link = "log"), 
#'                    data = ips_damage_clean, 
#'                       na.action = na.exclude)
#' tw_lag1 <- glmmTMB(damaged_volume_total_m3 ~ sum_beetle_lag1_scaled, family = tweedie(link = "log"), 
#'                    data = ips_damage_clean, 
#'                    na.action = na.exclude)
#' tw_lag2 <- glmmTMB(damaged_volume_total_m3 ~ sum_beetle_lag2_scaled, family = tweedie(link = "log"), 
#'                    data = ips_damage_clean, 
#'                    na.action = na.exclude)

#AIC(tw_lag0, tw_lag1,tw_lag2)
#simulateResiduals(gamma_m_lag0, plot = T)
#plot(allEffects(m3))


# get variance explained: for Gamma distribution, this is a pseudo r2




# Kriging -----------------------------------------------------------------------

# rasterzie
# xv <- rasterize(v_districts_simpl, dist_year_crop, fun=sum)
# plot(xv)
# 
# # convert spatRast to raster format
# back_rst <- raster(xv)
# 
# plot(v_ips_sum_simpl, add = T )
# plot(v_ips_sum_2015, add = T , col = 'red')
# 
# 
# plot(back_rst )
# plot(trap_data, add = T , col = 'red')
# 
# 
# 
# # test kriging: chat gpt
# 
# # Load necessary libraries
# library(terra)
# library(gstat)
# 
# # Create some sample spatial data
# set.seed(123)
# coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
# values <- sin(coords[,1]) + cos(coords[,2]) + rnorm(50, 0, 0.1)
# 
# # Create a SpatVector object with the sample data
# sp_data <- vect(data.frame(x = coords[,1], y = coords[,2], z = values), geom = c("x", "y"))
# 
# # Convert to SpatialPointsDataFrame
# sp_data_spdf <- as(sp_data, "Spatial")
# 
# # Create a raster template
# rast_template <- rast(xmin=0, xmax=10, ymin=0, ymax=10, resolution=0.1)
# 
# # Convert raster template to SpatialPixelsDataFrame
# rast_template_spdf <- as(rast_template, "Spatial")
# 
# # Create a gstat object for ordinary kriging
# g <- gstat(id="z", formula=z~1, data=sp_data_spdf)
# 
# # Compute and fit variogram
# vgm_model <- variogram(z~1, sp_data_spdf)
# fit_vgm <- fit.variogram(vgm_model, model=vgm("Sph"))
# 
# # Perform kriging
# krig_result <- predict(g, newdata=rast_template_spdf, model=fit_vgm)
# 
# # Convert the result to a SpatRaster
# krig_raster <- rast(krig_result)
# 
# # Plot the result
# plot(krig_raster, main="Kriging Interpolation")
# points(sp_data, col="red", pch=16)
# 
# 
# 
# 
# 
# # interpolate the beetle trap ata overBavaria/ year
# # link better with the reported damage data on teh district levels
# library(gstat)
# 
# v_ips_sum_2015 <- v_ips_sum[v_ips_sum$year == 2015, ]
# 
# #convert to sf object
# trap_data <- st_as_sf(v_ips_sum_2015)
# 
# 
# 
# interpolate(p, tps)
# 
# 
# # Create a vast_as_sf()# Create a variogram
# vgm_model <- variogram(sum_beetle ~ 1, trap_data)
# 
# 
# # Fit a variogram model
# fit_model1 <- fit.variogram(vgm_model, model = vgm(psill=1, model="Exp", range=1, nugget=0.1))
# 
# # Fit a variogram model with the estimated parameters from the semivariogram plot
# fit_model2 <- fit.variogram(vgm_model, model = vgm(psill = 8e7, model = "Exp", range = 100000, nugget = 2e7))
# 
# #grid_sel_ras_sf <- raster::raster(grid_sel_ras)
# #grid_spdf <- as(grid_sel_ras, "SpatialPixelsDataFrame")
# #crs(grid_spdf) <- crs(raster(dist_year_crop))
# #identical(crs(grid_spdf), crs(raster(dist_year_crop)))
# 
# # Step 2: Perform Kriging
# kriging_result <- krige(sum_beetle                 ~ 1, trap_data, newdata = back_rst, model = fit_model2)
# 
# # Optional: Convert kriging result back to RasterLayer if needed
# kriging_raster <- raster(kriging_result)
# 
# 
# 
# # Perform Kriging
# kriging_result <- krige(count ~ 1, trap_data, newdata =grid_sel_ras, model = fit_model2)
# grid_sel_ras
# 
# # Convert the result to a raster for visualization
# r_interpolated <- rast(kriging_result)
# 
# # Plot the interpolated raster
# plot(r_interpolated, main = "Interpolated Bark Beetle Counts")
# 
# 
# # test krigin on Meuse ---------------------
# library(sp)
# data(meuse)
# coordinates(meuse) = ~x+y
# data(meuse.grid)
# gridded(meuse.grid) = ~x+y
# m <- vgm(.59, "Sph", 874, .04)
# # ordinary kriging:
# x <- krige(log(zinc)~1, meuse, meuse.grid, model = m)
# spplot(x["var1.pred"], main = "ordinary kriging predictions")
# spplot(x["var1.var"],  main = "ordinary kriging variance")
# # simple kriging:
# x <- krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)
# # residual variogram:
# m <- vgm(.4, "Sph", 954, .06)
# # universal block kriging:
# x <- krige(log(zinc)~x+y, meuse, meuse.grid, model = m, block = c(40,40))
# spplot(x["var1.pred"], main = "universal kriging predictions")
# 
# 
# 
# 
# 

# merge disturbance rasters together; check n of cells---------------------------------------
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
  mutate(RS_sum = RS_wind_beetle + RS_fire +RS_harvest)

# type_1 = wind & beetles
# type_2 = fire
# type_3 = harvest



#### Merge RS with damage data -----------------------------------------------
df_all <- df_RS %>% 
  full_join(df_damage, by = c("ID" = "AELF_district_ID",
                              "year" = "Year")) %>% 
 # dplyr::filter(year %in% 2015:2020) %>% 
  #dplyr::filter(RS_wind_beetle <2000) %>% 
  dplyr::select(!c( AELF_ID, tree_species, unit, pest_species ))# %>%  # keep: AELF_district_name, AELF_name,

# get correlations
df_cor <- df_all %>% 
  dplyr::filter(!ID %in% c(0, 50104, 60613,62705)) %>% # exlude if there is missing data
  group_by(ID) %>% 
  dplyr::summarize(spearm_cor_beetle = cor(damaged_volume_total_m3, RS_wind_beetle, 
                                    method = "spearman", use = "complete.obs"),
                   spearm_cor_harvest = cor(damaged_volume_total_m3, RS_harvest, 
                                           method = "spearman", use = "complete.obs"),
                   spearm_cor_sum = cor(damaged_volume_total_m3, RS_sum, 
                                            method = "spearman", use = "complete.obs"))


#### correlations per location: beetle counts vs damage volume --------------------
df_cor_counts_damage <- 
  ips_damage_clean %>% 
  #dplyr::filter(!ID %in% c(0, 50104, 60613,62705)) %>% # exlude if there is missing data
  group_by(falsto_name,forstrev_1        ) %>% 
  dplyr::summarize(spearm_cor_lag0 = cor(damaged_volume_total_m3, sum_beetle, 
                                           method = "spearman", use = "complete.obs"),
                   spearm_cor_lag1 = cor(damaged_volume_total_m3, sum_beetle_lag1 , 
                                            method = "spearman", use = "complete.obs"),
                   spearm_cor_lag2 = cor(damaged_volume_total_m3, sum_beetle_lag2 , 
                                        method = "spearman", use = "complete.obs"))


#### correlations per disctrict: beetle counts vs RS damage --------------------------
# merge based on the district ID, subset the data as the traps do not cover the whole 
# bavaria's districts
df_traps_RS_damage <-  
ips_damage_clean %>%  
  left_join(df_all,
             by = c("forstrev_1" =   "ID",
                              "year" = "year",
                    "damaged_volume_total_m3" = "damaged_volume_total_m3") ) %>% 
  mutate(falsto_name = as.factor(falsto_name))


df_traps_RS_damage_mean <- df_traps_RS_damage %>% 
  mutate(pairID = as.factor(str_sub(falsto_name, 1, -3)),
         year_fact = as.factor(year)) %>% 
  group_by(pairID, year_fact) %>% 
  summarise(log_sum_beetle = mean(log_sum_beetle, na.rm =T),
            RS_wind_beetle = mean(RS_wind_beetle, na.rm =T)) %>% 
  na.omit() %>% 
  distinct()


# try with avrages
m0<- gam(RS_wind_beetle ~ s(log_sum_beetle) + s(pairID, bs = 're'), 
         family = nb(),data =df_traps_RS_damage_mean) 

m_fixed<- gam(RS_wind_beetle ~ s(log_sum_beetle, k = 15) + year_fact, 
         family = nb(),data =df_traps_RS_damage_mean) 

m_fixed<- gam(RS_wind_beetle ~ s(log_sum_beetle, k = 15) + year_fact, 
              family = nb(),data =df_traps_RS_damage_mean) 


check_concurvity(m_fixed)
  summary(m_fixed)
  plot(m_fixed, page =1)
  
# try simple gams: how do beetle counts represent the RS damage?

m0<- gam(RS_wind_beetle ~ s(sum_beetle) + s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 
m1<- gam(RS_wind_beetle ~ s(sum_beetle_lag1 ) + s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 
m2<- gam(RS_wind_beetle ~ s(sum_beetle_lag1)+ s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 


m0_log<- gam(RS_wind_beetle ~ s(log_sum_beetle) + s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 
m1_log<- gam(RS_wind_beetle ~ s(log_sum_beetle_lag1 ) + s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 
m2_log<- gam(RS_wind_beetle ~ s(log_sum_beetle_lag1)+s(falsto_name, bs = 're'), family = nb(),data =df_traps_RS_damage) 


summary(m)
AIC(m0, m1, m2,m0_log,m1_log,m2_log)
plot(m0, page = 1, shade = T)
plot(m1, page = 1, shade = T)
plot(m2_log, page = 1, shade = T)

summary(m2_log)

check_concurvity(m2_log)


df_cor_counts_damage_long <- df_cor_counts_damage %>% 
  pivot_longer(-c(falsto_name, forstrev_1))

df_cor_counts_damage_long %>% 
  ggplot(aes(y = value, group =factor(name), fill = name)) +
  geom_boxplot()


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
cor_shp <- sf_simpled %>% 
  left_join(df_cor, by = c("forstrev_1" = "ID"))

all_shp <- sf_simpled %>% 
  left_join(df_all, by = c("forstrev_1" = "ID"))

# add data to the full geometry to run teh intersection with the trap data
#all_shp_complex <- sf_districts_trans %>% 
#  left_join(df_all, by = c("forstrev_1" = "ID"))

cor_shp_counts_dmg <- sf_simpled %>% 
  right_join(df_cor_counts_damage_long, by = c("forstrev_1"))




### --------------------------

# Check correlation: is there any correlation, overall pattern? oevrall sum over years for damage/RS data? 
df_all %>% 
  ggplot(aes(x = RS_harvest,
             y= damaged_volume_total_m3)) +
  geom_point() +
  geom_smooth()


df_all %>% 
  ggplot(aes(x = RS_wind_beetle,
             y= damaged_volume_total_m3)) +
  geom_point() +
  geom_smooth()  # extreme is RS_wind_beetle

df_all %>% 
  ggplot(aes(x = RS_sum,
             y= damaged_volume_total_m3)) +
  geom_point() +
  geom_smooth()



# plot correlation with beetle counts & tree damage ----------------------------------
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


# ad damage data to district geometry



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
       subtitle = "Correlation between spruce damage [m3] and RS wind-beetle [pixel counts] by region") +
  theme_minimal()  #




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

#damaged_volume_total_m3

View(df_all)

# log values
df_all$log_RS_wind_beetle          <- log(df_all$RS_wind_beetle + 1)
df_all$log_damaged_volume_total_m3<- log(df_all$damaged_volume_total_m3 + 1)

# try different functions and seettings - not finished!
gam_RS_1 <- gam(log_RS_wind_beetle ~ s(log_damaged_volume_total_m3, k = 3) + s(year,
                                                                               k = 3) + s(ID, bs = 'fs', k = 30) +, family = gaussian(),method = "REML", data = df_all)



# GAM Trap vs RS damage data ----------------------------------------------



summary(gam_RS_1)
plot(gam_RS_1, page = 1)
library(gratia)
gam.check(gam_RS_1)

