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
  


# Libs --------------------------------------------------------------------------
library(dplyr)
library(raster)
library(terra)

library(sf)
library(data.table)
library(ggplot2)
library(ggpubr)
#library(ggspatvector)

library(lubridate)
library(fasterize)
library(tidyr)


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

length(sf_ID)
length(df_ID)


# Process rasters ------------------------------------------------------------
# there is an aditiona row in the damage data with district ID 0 - I can remove that
hist(df_damage$damaged_volume_total_m3)

my_crs = st_crs(disturb_year)

# change projection & convert to terra
sf_districts_trans <- st_transform(sf_districts, crs = my_crs) # st_crs(disturb_year))
v_districts  <-  terra::vect(sf_districts_trans) # convert to terra


# add beetle counts to trap locations per year 
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)) 

# convert to sf object, crs 3035
ips_sum_sf <-  xy_sf_expand %>%  left_join(ips_sum, by = c("falsto_name", 'year')) 
#st_is_valid(ips_sum_sf_trans)

# change projection 
ips_sum_sf_trans  <- st_transform(ips_sum_sf, crs = st_crs(disturb_year))
v_ips_sum         <- vect(ips_sum_sf_trans)


# works well for terra! 
plot(forest)
plot(disturb_type)
plot(v_ips_sum, add = T)


plot(v_districts)
plot(v_ips_sum, add = T)


# simplify terra object:
v_districts_simpl <- simplifyGeom(v_districts, tolerance=20, preserveTopology=TRUE, makeValid=TRUE)


# crop the rasters: they are for whole germany, crop to Bavaria
dist_year_crop <- terra::crop(disturb_year, v_districts_simpl)
dist_type_crop <- terra::crop(disturb_type, v_districts_simpl)
forest_crop    <- terra::crop(forest,       v_districts_simpl)

# process spatial data: check geometries, change projection to Corneliuses data, convert to raster the AELF districts
crs(disturb_type) <- crs(dist_year_crop)
crs(forest_crop)  <- crs(dist_year_crop)


# Simplify the df districts geometry -------------------------------
sf_simpled <- st_simplify(sf_districts_trans, dTolerance = 20, preserveTopology = TRUE)

# check validity
st_is_valid(sf_simpled)

sf_simpled = st_make_valid(sf_simpled)
st_is_valid(sf_simpled)

#plot(sf_simpled, max.plot = 1)


# Convert damage data to raster --------------------------------------------

# Assuming dist_year_crop is a SpatRaster and sf_districts_trans is an sf object
# Convert sf object to SpatVector
# add data to the full geometry to run teh intersection with the trap data
v_damage <- sf_districts_trans %>% 
  left_join(df_damage, by = c("forstrev_1" = "AELF_district_ID")) %>% 
  vect() 


# Assuming 'forstrev_1' in SpatVector and 'AELF_district_ID' in DataFrame are the matching join keys
#v_damage2 <- merge(v_districts_simpl, df_damage, by.x="forstrev_1", by.y="AELF_district_ID")

#plot(v_damage2)

# Check if identical --------------
identical(crs(v_damage), crs(v_damage2))  # yes
identical(ext(v_damage), ext(v_damage2))

# !!!! convert damage data to raster???? -----------------------------------------
# maybe skip? 
v_districts_trans <- vect(sf_districts_trans)

plot(v_districts_trans)
# Set the CRS of the vector to match the raster's CRS
crs(v_districts_trans) <- crs(dist_year_crop)

# Get the extent of the vector
ext <- ext(v_districts_trans)

# Convert extent to a SpatExtent object (if needed, usually ext() returns SpatExtent)
ext <- as(ext, "SpatExtent")


# Intersect the SpatVector with the extent (converted to a SpatRaster for compatibility)
# First, create a SpatRaster from the extent using the original raster's resolution
ext_ras <- rast(ext, nrow=nrow(dist_year_crop), 
                ncol=ncol(dist_year_crop), crs=crs(dist_year_crop))

# Rasterize the intersected vector back onto the original raster grid
# Assuming 'forstrev_1' is a field in your vector data
grid_damage <- rasterize(v_districts_trans, 
                          dist_year_crop, 
                          field="damaged_volume_total_m3")

# -*----------------------------------------------------


###### convert districts poly to raster ------------------------------------

# get extend and set correct projectipn
ext <- as(extent(sf_simpled), 'SpatialPolygons')  # treat raster as from raster, not terra package
ext <- st_as_sf(ext)

# efine teh crs first - it was NA
st_crs(ext) <- st_crs(raster(disturb_year))

ext  <- st_transform(ext, crs = st_crs(disturb_year))

st_crs(sf_simpled) == st_crs(ext) 
st_crs(disturb_year) == st_crs(ext) 

# Create grid index for each pixel
grid_sel     <- st_intersection(st_as_sf(sf_simpled), st_as_sf(ext))
grid_sel_ras <- fasterize::fasterize(grid_sel, raster(dist_year_crop), field = "forstrev_1") # name the new rasters as polygon ids  #, field = "forstrev_1"

# extract beetle sums positions
ips_sum_sf_trans <- st_transform(ips_sum_sf, crs = st_crs(disturb_year))

# try of geometry fits
plot(grid_sel_ras)
plot(ips_sum_sf_trans ,add = T)


# extract the reviers ID (from the districts) to the beetle counts data 
district_ID <- terra::extract(rast(grid_sel_ras), vect(ips_sum_sf_trans), df = TRUE )

# merge district number data to the beetle counts  
ips_damage_district<- cbind(district_ID, ips_sum_sf_trans)

ips_damage_district <- ips_damage_district %>% 
  rename(forstrev_1 = layer)


# join with damage data
ips_damage_merge <- ips_damage_district %>% 
  left_join(df_damage, by = c('forstrev_1' = 'AELF_district_ID',
                              'year' = 'Year')) %>% 
  dplyr::select(-c(ID, id, AELF_district_name,AELF_name, AELF_ID,
                 pest_species, tree_species, unit,
                 globalid,geometry )) %>% 
  # lag tree damage: does beetle counts (or their lags) predict tree damage?:
  group_by(falsto_name) %>% 
  arrange(year, .by_group = TRUE) %>%
  mutate(sum_beetle_lag1  = lag(sum_beetle, n = 1, default = NA),
         sum_beetle_lag2  = lag(sum_beetle, n = 2, default = NA, order_by = year))

# remove NAs
ips_damage_clean <- na.omit(ips_damage_merge)

# get log values - makes it nice for damaged volume
ips_damage_clean$log_sum_beetle                   <- log(ips_damage_clean$sum_beetle + 1)
ips_damage_clean$log_sum_beetle_lag1              <- log(ips_damage_clean$sum_beetle_lag1 + 1)
ips_damage_clean$log_sum_beetle_lag2              <- log(ips_damage_clean$sum_beetle_lag2 + 1)
ips_damage_clean$log_damaged_volume_total_m3 <- log(ips_damage_clean$damaged_volume_total_m3 + 1)


# try correlation
cor(ips_damage_clean$sum_beetle, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')
cor(ips_damage_clean$sum_beetle_lag1, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')
cor(ips_damage_clean$sum_beetle_lag2, ips_damage_clean$damaged_volume_total_m3, method = 'spearman')



# try correlation: pearson
cor(ips_damage_clean$log_sum_beetle, ips_damage_clean$log_damaged_volume_total_m3)
cor(ips_damage_clean$log_sum_beetle_lag1, ips_damage_clean$log_damaged_volume_total_m3)
cor(ips_damage_clean$log_sum_beetle_lag2, ips_damage_clean$log_damaged_volume_total_m3)

#library(glm)
m_lag0 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle, data = ips_damage_clean, family = gaussian())
m_lag1 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle_lag1, data = ips_damage_clean, family = gaussian())
m_lag2 <- glm(log_damaged_volume_total_m3 ~ log_sum_beetle_lag2, data = ips_damage_clean, family = gaussian())

# 
AIC(m_lag0,m_lag1,m_lag2)


# use a gamma to predictthe tree damage (continuous,positive)
# Fitting a Gamma model
gamma_m_lag0 <- glm(damaged_volume_total_m3 ~ sum_beetle, family = Gamma(link = "log"), data = ips_damage_clean)
gamma_m_lag1 <- glm(damaged_volume_total_m3 ~ sum_beetle_lag1, family = Gamma(link = "log"), data = ips_damage_clean)
gamma_m_lag2 <- glm(damaged_volume_total_m3 ~ sum_beetle_lag2, family = Gamma(link = "log"), data = ips_damage_clean)

AIC(gamma_m_lag0, gamma_m_lag1, gamma_m_lag2)
# r2(gamma_m_lag0, gamma_m_lag1, gamma_m_lag2)

# try tweedie???

# get variance explained: for Gamma distribution, this is a pseudo r2


# does beetle counts and their lags predict tree damage???
hist(log(ips_damage_clean$sum_beetle_lag1))
hist(log(ips_damage_clean$damaged_volume_total_m3))
  

pairs(damaged_volume_total_m3~ sum_beetle + sum_beetle_lag1 + sum_beetle_lag2, data =  ips_damage_merge)
#cor(damaged_volume_total_m3~ sum_beetle + sum_beetle_lag1 + sum_beetle_lag2, data =  ips_damage_merge)


# 
ips_damage_merge %>% 
  ggplot(aes(x = sum_beetle_lag1,
             y = damaged_volume_total_m3)) + 
  geom_point() + 
  geom_smooth(method = 'lm')








# merge all raster data together; check n of cells
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

# What to get: per region: how many cells, forest, disturbances, harvest& bark beetles

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
all_shp_complex <- sf_districts_trans %>% 
  left_join(df_all, by = c("forstrev_1" = "ID"))

# split by years
all_shp_complex_list <- split(all_shp_complex, all_shp_complex$year)



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




# merge damage data with the beetle counts (IPS) ---------------------------------------------
all_shp_complex  #

all_shp_complex_list <- split(all_shp_complex, all_shp_complex$year)

ips_sum_sf_trans_list <- split(ips_sum_sf_trans, ips_sum_sf_trans$year)


plot(all_shp_complex_list[[5]]['damaged_volume_total_m3'])
plot(ips_sum_sf_trans_list[[5]]['sum_beetle'], col = "red", add = T)

# Assuming the lists poly_list and point_list correspond to each year from 2015 to 2021
results <- lapply(seq_along(all_shp_complex_list), function(i) {
  # Perform a spatial join where points are matched within polygons
  joined_sf <- st_join(ips_sum_sf_trans_list[[i]], all_shp_complex_list[[i]], join = st_within)
  
  # Convert to a data frame for safer manipulation with dplyr, if necessary
  # and immediately use dplyr to handle attributes
  joined_df <- as.data.frame(joined_sf) %>%
    dplyr::select(polygon_data = "damaged_volume_total_m3", 
                  point_data = "sum_beetle", 
                  everything())
  
  # Optionally, filter out non-matches if there are any
  joined_df <- filter(joined_df, !is.na(polygon_data))
  
  return(joined_df)
})

# Combine all yearly data into one data frame
final_df <- bind_rows(results, .id = "year")

# !!! continue from here
#extracted_data <- terra::extract(x = v_districts_trans, y = v_ips_sum) # works, but first I need to add there a damage data!









# convert to terra
v_ips_sum  <- vect(ips_sum_sf_trans)


ggplot(cor_shp) +
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


# ad damage data to district geometry




# select polygons IPS_sum -----------------------------------------------------------------

#selected_dist <- terra::intersect(v_ips_sum, v_districts_trans )

plot(v_districts_trans)
plot(v_ips_sum, add = T, col = "red")









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
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3))




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





