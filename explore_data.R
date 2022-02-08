

# -----------------------------------
# Explore bark beetle data
# -----------------------------------
rm(list=ls()) 


# Libs -----------------------------------------------------------
library(lubridate)
library(dplyr)
library(ggplot2)
#library(rgdal)   # for spatial data
library(sf)

# get data
path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/Fwd__Borkenkäferforschung__Datentransfer'

# Load RData table
load(paste(path, 'BoMo_2015_2021_Rohdaten.RData', sep = "/"))

# Read XY coordinates -----------------------------------------------------------
# need to save the csv fil as a new file, otherwise it did not work
# or use teh  version 3: now it worked with the original one! 
xy <- read.table(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit.csv", sep = '/'), 
                 dec=',', sep = ';', 
                 header = T, quote="\"", 
                 skip = 1, #header = T,
                 skipNul = TRUE) # ,fileEncoding="UTF-16LE"

#readLines(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit.csv", sep = '/'))
# Investigate beetle data: Bavaria
head(xy)


# Convert coordinates to sf object
xy_sf <- st_as_sf(xy, coords = c('x_coord', 'y_coord' ), crs = 32623 )  # What is a coordinate system??? for forestry in bavaria?

plot(xy_sf["OBJECTID"])


# Get dynamics data -------------------------------------------------------------- 

# Change name of the table
# for dat
dat <- Daten_B01


# "Buchdrucker"   "Kupferstecher"

# what year?
unique(dat$kont_dat)

str(dat)  # datum is in POSIXct format, ut does not contain time: can convert just to date

# convert to date
dat$kont_dat <- as.Date(dat$kont_dat)

str(dat)

# Species:
# -----------------------------------
# Buchdrucker  : Ips typographus 
# Kupferstecher: Pityogenes chalcographus      


# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year  = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day   = lubridate::day(kont_dat))


# get basic stats: count by beetles over year, counts by records
nrow(dat)

# Some basic plotting
dat %>% 
  group_by(art, year) %>% 
  summarise(my_mean = mean(fangmenge, na.rm = T)) %>% 
  ggplot(aes(x = year,
             y = my_mean,
             fill = art,
             group = year)) +
  geom_point() +
 # ylim(0, 250000) +
  facet_grid(.~art)


# 
summary(dat)

# Join data with XY coordinates by globalid   ---------------------------------------------------

dat <- dat %>% 
  mutate(globalid =  gsub("\\{|\\}", "", globalid))  


# for coordinates table

xy_sf <-
  xy_sf %>% 
  mutate(globalid = toupper(globalid))


# Merge data by globalid
# 
merged_df <- xy_sf %>% 
  right_join(dat, 'globalid')






# Plots spatial data ------------------------------------

# 
windows()
merged_df %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot() + 
  geom_sf(aes(color = fangmenge ))  + # , size = 0.5 , size = AREA
  scale_color_continuous(low = "lightgreen", 
                         high = "darkgreen",
                         space = "Lab", 
                         na.value = "red", guide = "colourbar")#+




# R animate: ------------------------------------------------

# Sums by months and years

buch_sum_ysr <- dat %>%
  filter(art == 'Buchdrucker') %>% 
  group_by(globalid, year) %>% 
  summarize(sum = sum(fangmenge, na.rm = T))

# Merge summarized data with XY 
buch_df <- xy_sf %>% 
  right_join(buch_sum_ysr, 'globalid')




# Get data:
library(gapminder)

# Charge libraries:
library(ggplot2)
library(gganimate)
library(transformr)


# need to recheck!! 
ggplot(buch_df) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df,
          aes(color = sum)) + # , size = Age, size = 0.8size by factor itself!
  viridis::scale_color_viridis(discrete = FALSE, option = "viridis",
                               na.value = "red") +
 # annotation_scale(location = "bl", width_hint = 0.4) +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  # gganimate specific bits: ------------------
labs(title = 'Changes in Tree Height over years {current_frame}',
     color  = "Tree height [m]") +
  transition_manual(year) +
  #transition_time(year) +
  ease_aes('linear')





