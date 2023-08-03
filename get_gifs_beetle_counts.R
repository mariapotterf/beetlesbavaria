

# Explore GIFS


# Read data
load(file = "outData/ips.Rdata")


# R animate: ------------------------------------------------

# Get spatial data libs -----------------------------------------------------------------
library(gapminder)

# Charge libraries:
library(ggplot2)
library(gganimate)
library(transformr)
library(gifski)  # needed for correct rendering of the animation

library(rnaturalearth) # for map data
library(ggspatial)


# Get colors ------------------------------------
library(RColorBrewer)
library(scales)
library(viridis)


# Prepare data --------------------------------------------------

# sum beetle counts by months& year
ips.year.sum <-   dat.ips.clean %>%     
  group_by(year, month, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
  dplyr::select(bav_sum, globalid, year) #%>% 


# Merge spatial data with sum data by globalid
ips.year.sum_sf <- xy_sf %>% 
  right_join(ips.year.sum, 'globalid')




# get spatial data ------------------------------------------------------------------ 
de_sf <- ne_states(country = "germany", returnclass = "sf")

# Get only bavaria
bav_sf <- de_sf %>% 
  dplyr::filter(name_en == "Bavaria")



# make map ----------------------------------------------------------------

display.brewer.pal(7, "BrBG")

buch_df_f <- ips.year.sum_sf %>% 
  filter(bav_sum > 3000)


# GIF: months and years -------------------------------------------------------- 
p<-ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df_f,
          aes(color = bav_sum,
              size = bav_sum)) + # , size = Age, size = 0.8size by factor itself!
  # viridis::scale_color_viridis(discrete = FALSE, option = "viridis",
  #                             na.value = "red") +
  scale_fill_gradientn(colours = colorspace::heat_hcl(7)) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  # gganimate specific bits: ------------------
labs(title = 'IT beetle population counts {current_frame}',
     color  = "IT Beetle [count]") +
  transition_manual(time) +
  #transition_time(year) +
  ease_aes('linear')


# animate with the gifski renderer
animate(p, renderer = gifski_renderer())




