# Make a map

library(ggplot2)
library(rnaturalearth)
library(sf)
library(ggspatial)
library(cowplot)
library(terra)
library(raster)
library(lubridate)
library(tidyterra)
library(fasterize)

# Load world data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Filter for Germany
germany <- world[world$admin == "Germany",]


# Get spatial data for each trap
xy        <- vect("C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/outSpatial/xy_fin_3035.gpkg") # read trap location
xy2       <- terra::project(xy, crs(germany))  # coordinate system from the DWD data: Germany
xy2 <- st_as_sf(xy2)


#forest   <- rast("C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/outSpatial/bav_fortype_ext30_int2u_LZW.tif") # read trap location
#forest2  <- terra::project(forest, crs(germany))  # coordinate system from the DWD data: Germany
#forest2  <- raster::raster(forest)  # coordinate system from the DWD data: Germany
size(forest)


# Get state boundaries within Germany, if available in rnaturalearth; otherwise, focus on Germany
tryCatch({
  states <- ne_states(country = "germany", returnclass = "sf")
  bavaria <- states[states$name == "Bayern",]
}, error = function(e) {
  bavaria <- germany # Fallback to Germany if state-level data isn't available
})


# make a raster from bavaria shp to properly fit the forest data
bav_3035 <- terra::project(vect(bavaria), 'EPSG:3035')
r <- raster(st_as_sf(bav_3035), res = 30)
r <- fasterize(st_as_sf(bav_3035), r)

# convert sf to terra object
r<-rast(r)

# Resample bav_forest to match r's resolution
r_resampled <- terra::resample(r, forest)

# maks bavaria forest to only Bavaria shape
bav_forest <- crop(forest, r_resampled)
bav_forest_mask <- mask(bav_forest, r_resampled)

# Create a new raster where only the cells with a value of '2' are kept; all others are set to NA
reclass_matrix <- matrix(c(0, 1, NA,
                           2, 2, 2), byrow = TRUE, ncol = 3)
bav_forest_recl <- classify(bav_forest_mask, rcl = reclass_matrix, right=NA, include.lowest=TRUE)


# increase the pixel size to improve plottin
bav_forest_mask_simpl <- terra::aggregate(bav_forest_recl, fact = 5, fun="max" ) 
hist(bav_forest_mask_simpl)
plot(bav_forest_mask_simpl)
dim(bav_forest_mask_simpl)

#bav_forest_recl <- classify(bav_forest_mask_simpl, rcl = reclass_matrix)
ggplot() +
  geom_spatraster(data = bav_forest_mask_simpl) +
  geom_sf(bavaria)


names(bav_forest_mask_simpl)

bav_forest_mask_simpl <- terra::resample(bav_forest_mask_simpl, method = 'near')

sum(values(bav_forest_recl) == 2, na.rm = TRUE)

terra::dType(bav_forest_recl)

# Convert NaN to NA
values(bav_forest_mask_simpl)[is.nan(values(bav_forest_mask_simpl))] <- NA


# Convert raster to data frame
bav_forest_df <- as.data.frame(bav_forest_recl, xy = TRUE) %>%
  mutate(value = factor(bav_fortype_ext30_int2u_LZW))


# Plotting with geom_tile
ggplot(bav_forest_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_d(name = "Forest Type", drop=FALSE) +
  coord_equal() +
  theme_minimal() +
  labs(fill = "Value")

# Now, plot the filtered raster with ggplot2
ggplot() +
  #geom_raster(aes(data = bav_forest_recl, fill= bav_fortype_ext30_int2u_LZW))
  geom_spatraster(data = bav_forest_mask_simpl, fill= bav_forest_mask_simpl)# +
  scale_fill_manual(values = c(`2` = "blue", `0` = NA), na.value = NA) #+
  theme_minimal() +
  labs(fill = "Forest Type") +
  theme(legend.position = "none")  # Optional: Remove the legend if not needed



# Define Central Europe bounds manually for the inset map
central_europe_bounds <- st_bbox(c(xmin = 5, xmax = 20, ymin = 44, ymax = 55), crs = st_crs(world))

# Main Map of Bavaria
main_map <- ggplot() +
  #geom_sf(data = world, fill = "lightgrey") +
  geom_sf(data = bavaria, fill = "grey") +
  #geom_sf(data = forest2, fill = 'darkgreen') +
  #geom_spatraster(data = bav_forest_mask_simpl, fill = 'darkgreen') +
  #geom_spatraster(data = bav_forest_recl) +
  #scale_fill_manual(values = c("2" = "darkgrey"), na.value = NA) +
  geom_sf(data = xy2) +
  theme_void() +
  #theme_minimal() +
  ggtitle("") +
  # Add a north arrow
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_minimal) +
  # Add a scale bar
  annotation_scale(location = "br", width_hint = 0.5, style = 'ticks')

(main_map)
unique(values(bav_forest_mask_simpl))
# Inset Map of Central Europe
inset_map <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = 'white') +
  geom_sf(data = germany, fill = "darkgrey", color = 'darkgrey') +
  geom_sf(data = bavaria, fill = "grey50", color = 'darkgrey') +
  coord_sf(xlim = c(-2,20), 
           ylim = c(46,59)) +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))  # Adding a frame

inset_map
plot_grid(inset_map, main_map, nrow = 1, rel_widths = c(1, 2.3))

# Combine both maps
main_map +
  annotation_custom(ggplotGrob(inset_map), xmin = 14, xmax = 18, ymin = 48, ymax = 49.5) #+
  #theme(plot.margin = margin(1, 1, 1, 1, "cm"))
