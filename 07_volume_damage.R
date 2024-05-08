# Proces volume data

# from LWF: the damage volume+ geometries
# read tables
# merge geometries
# explore damage per year
#  read ata from Cornelius: do they correlate somehow? YES!
# how to link them with beetle trap data?
# select the respective region per year
# do beetle counts fit over forest damage??



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

library(lubridate)
library(fasterize)
library(tidyr)


# read beetle data -------------------------------------------------------------
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



# read districts shp ------------------------------------------------------------
sf_districts <-  st_read('rawData/damage_volume/AELF_Revier_20240502/AELF_Revier_bis20240424_LWF_20240502.shp')
df_damage    <-  fread('rawData/damage_volume/Bavaria_TreeDamageData_IpsTypographus_2015-2021.csv', dec = ',')
v_districts  <-  terra::vect(sf_districts) # convert to terra

# Read data
disturb_year <- rast('rawData/germany114/disturbance_year_1986-2020_germany.tif')
disturb_type <- rast('rawData/fire_wind_barkbeetle_germany.tif')
forest       <- rast('rawData/germany114/forestcover_germany.tif')

# check data: how many regions I have?
sf_ID   <- unique(sf_districts$forstrev_1)  # 336 districts numbers
df_ID   <- unique(df_damage$AELF_district_ID)  # 337 district numbers

length(sf_ID)
length(df_ID)


# there is an aditiona row in the damage data with district ID 0 - I can remove that
hist(df_damage$damaged_volume_total_m3)

# change projection & convert to terra
sf_districts_trans <- st_transform(sf_districts, crs = st_crs(disturb_year))
ips_sum_sf_trans <- st_transform(ips_sum_sf, crs = st_crs(disturb_year))

v_districts_trans  <- vect(sf_districts_trans)
v_ips_sum  <- vect(ips_sum_sf_trans)


# crop the rasters: they are for whole germany, crop to Bavaria
dist_year_crop <- terra::crop(disturb_year, v_districts_trans)
dist_type_crop <- terra::crop(disturb_type, v_districts_trans)
forest_crop    <- terra::crop(forest, v_districts_trans)

# process spatial data: check geometries, change projection to Corneliuses data, convert to raster the AELF districts
#crs(v_districts) <- crs(dist_year_crop) 
#crs(v_districts)
crs(disturb_type) <- crs(dist_year_crop)
crs(forest_crop) <- crs(dist_year_crop)


# select polygons IPS_sum -----------------------------------------------------------------
# split beetle by years: create a list by years
xy_ls <- terra::split(v_ips_sum, c('year')) #,# "OBJECTID" # 

# extract values from polygons
extracted_data <- terra::extract(x = v_districts_trans, y = v_ips_sum) # works, but first I need to add there a damage data!

#selected_dist <- terra::intersect(v_ips_sum, v_districts_trans )

plot(v_districts_trans)
plot(v_ips_sum, add = T, col = "red")


# convert to districts poly to raster ------------------------------------

# get extend and set correct projectipn
ext <- as(extent(sf_districts_trans), 'SpatialPolygons')  # treat raster as from raster, not terra package
ext <- st_as_sf(ext)
st_crs(ext) <- st_crs(raster(dist_year_crop))

# Create grid index for each pixel
grid_sel     <- st_intersection(st_as_sf(sf_districts_trans), st_as_sf(ext))
grid_sel_ras <- fasterize::fasterize(grid_sel, raster(dist_year_crop), field = "forstrev_1") # name the new rasters as polygon ids  #, field = "forstrev_1"

plot(grid_sel_ras)
# 

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



# Merge RS with damage data -----------------------------------------------
df_all <- df_RS %>% 
  full_join(df_damage, by = c("ID" = "AELF_district_ID",
                              "year" = "Year")) %>% 
 # dplyr::filter(year %in% 2015:2020) %>% 
  #dplyr::filter(RS_wind_beetle <2000) %>% 
  dplyr::select(!c( AELF_ID, tree_species, unit, pest_species ))# %>%  # keep: AELF_district_name, AELF_name,


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




# is there any correlation, overall pattern? oevrall sum over years for damage/RS data? ----
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

# merge counts to spatial data
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)#,
                   #log_sum_beetle = log(sum_beetle)
  ) #%>% 

# convert to sf object
ips_sum_sf <-  xy_sf_expand %>%  left_join(ips_sum, by = c("falsto_name", 'year')) #%>% # df_xy3035
# mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) #%>% 

# 3035 crs








# RS plots -----------------------------------------------------------------


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



# merge damage data to geometry data ----------------------------------
# Simplify the geometry of your sf object
sf_simpled10 <- st_simplify(sf_districts, dTolerance = 20, preserveTopology = TRUE)

plot(sf_simpled10, max.plot = 1)


cor_shp <- sf_simpled10 %>% 
  left_join(df_cor, by = c("forstrev_1" = "ID"))

all_shp <- sf_simpled10 %>% 
  left_join(df_all, by = c("forstrev_1" = "ID"))



# show spearman correlation coefficient between volume damage and RS observation (wind, beetle)

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




# Map damage [m3] per year ------------------------------------------------

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





