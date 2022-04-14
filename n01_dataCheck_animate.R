
# here are my git changes in Home dir


# -----------------------------------
# Explore bark beetle data
# -----------------------------------
rm(list=ls()) 


source('myPaths.R')

# Libs -----------------------------------------------------------
library(lubridate)
library(dplyr)
library(sf)
library(tidyverse)
library(ggplot2)
library(zoo)
library(data.table)

# get data
path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/Fwd__Borkenkäferforschung__Datentransfer'
out_path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria'

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
xy_sf <- st_as_sf(xy, coords = c('x_coord', 'y_coord' ), 
                  crs = 32632 )  # What is a coordinate system??? for forestry in bavaria?

plot(xy_sf["OBJECTID"])

# Change projection to common european one: 3035
xy_sf <- st_transform(xy_sf, crs = 3035)


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

# The counts of Kupferstecher needs to be divided by 10!!


# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year  = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day   = lubridate::day(kont_dat),
                doy   =  lubridate::yday(kont_dat) + 1)  # as POXIT data has January 1st at 0


# get basic stats: count by beetles over year, counts by records/traps
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


# Split data between years: , until 2018, > 2018 - drought years, 
# need to standardize data between years: get mean yearly counts?
# maybe also standarziye by number of plots??? does number od plot vary over years?
dat %>% 
  mutate(time_cl = case_when(
    year <  2018 ~  '2014-2017',
    year >= 2018 ~  '2018-2020')) %>% 
  group_by(art, time_cl) %>% 
  summarize(my_mean = mean(fangmenge, na.rm = T),
            my_med = median(fangmenge, na.rm = T))# %>% 
  

# Get number of traps per year to calculate the mean beetle count per trap
#trap_n <- 
  dat %>% 
    group_by(art, year, globalid) %>% 
    count()
  
# How may traps in total?  
length(unique(dat$globalid)) # 490 traps

# how many traps per year?
dat %>% 
  group_by(art, year) %>% 
  distinct(globalid) %>% 
  count()


# --------------------------------------------------------------------
#              Standardize beetle counts/trap/year
# --------------------------------------------------------------------


# how to get the average number of beetles per trap per year?

# Calculate sum of beetles/year
# sum of traps/year
# divide by frequency of recording: revisit time? : 
# how many times the traps have been revisited? 'kont_dat': lenght(unique)..


# Try on filtered data: globalid  == {0718EEEE-A45E-40C5-9BE6-83EDE98B7035}
dat %>% 
  filter(globalid == '{0718EEEE-A45E-40C5-9BE6-83EDE98B7035}') %>% 
  group_by(art, year) %>% 
  summarize(sum_beetle = sum(fangmenge, na.rm = T),
            freq_visit  = length(unique(kont_dat))) #%>% 
  #distinct(monsto_name, year)

#     art            year sum_beetle freq_visit
#   <chr>         <dbl>      <int>      <int>
# 1 Buchdrucker    2015      12628         26
# 2 Buchdrucker    2016      14907         23
# 3 Kupferstecher  2015      19632         26
# 4 Kupferstecher  2016      14506         23


# Check if calculation corect of revisit times:
dat %>% 
  filter(globalid == '{0718EEEE-A45E-40C5-9BE6-83EDE98B7035}' & art == "Buchdrucker" & year == 2015) %>% 
  #group_by(art, year) %>% 
  distinct(kont_dat) %>% 
  pull()

# some sites have been reisited every 7 days!
# calculation works great!

# Calculate the mean number of beetles per trap over the season
# corrected for revisit frequency
dat_avg <- dat %>% 
  group_by(art, year, globalid) %>% 
  summarize(sum_beetles = sum(fangmenge, na.rm = T),
            sum_trap    = length(unique(globalid)),
            freq_visit  = length(unique(kont_dat)),
            avg_beetles_trap = sum_beetles/freq_visit)

# does the average number of beetles per trap differ between locations?
# between years?
dat_avg %>% 
ggplot(aes(x = factor(year),
           y = avg_beetles_trap,
           fill = art)) +
  geom_boxplot() +
  facet_wrap(.~art , scales = 'free_y') +
  theme(legend.position = 'bottom')


dat <- dat %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) 

# Correctly order levels
dat <- dat %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))






# 
# Make a boxplot of mean values between populations and years intervals:
dat_avg <- dat_avg %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) %>% #,
                         
    #) %>% 
  group_by(art, time_cl) # %>%

# Correctly order levels
dat_avg <- dat_avg %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))
                                 
     


# geom smooth plot
dat_avg %>% 
  ggplot() +
  geom_boxplot(aes(x = drought_period,
                   y = avg_beetles_trap,
                   fill = art)) +
#  scale_y_continuous(limits = c(0,25000)) +
  facet_wrap(.~art, scales = 'free') +
  theme(legend.position = 'bottom')
  # summarize(my_mean = mean(fangmenge, na.rm = T),
  #           my_med = median(fangmenge, na.rm = T))# %>% 



  
# Check development of counts over one year
windows()
dat %>% 
  filter(month > 4 ) %>%             #  &  art == "Buchdrucker"
  ggplot(aes(x = doy,                # DOY = Day of the Year
             y = fangmenge,
             color = factor(year))) +
    #geom_point() +
    geom_smooth(method = 'loess') +  #'gam' 
  scale_y_log10()+
  scale_color_viridis_d() + 
  facet_grid(.~ art)
  #coord_trans(ytrans="pow10")


# for overall counts:
windows()
dat %>% 
  filter(month > 4 ) %>%             #  &  art == "Buchdrucker"
  ggplot(aes(x = doy,                # DOY = Day of the Year
             y = fangmenge,
             color = drought_period )) +
  #geom_point() +
  geom_smooth(method = 'loess') +  #'gam' 
  scale_y_log10()+
  scale_color_viridis_d() + 
  facet_grid(.~ art)
#coord_trans(ytrans="pow10")









# Calculate the sums per day & year, then plot the data
dat %>% 
  filter(month > 4 & year == 2017) %>%             #  &  art == "Buchdrucker"
  group_by(doy, art, year)  %>% 
  summarize(my_sum = sum(fangmenge, na.rm = T)) #  %>% 
  ggplot(aes(x = doy,                # DOY = Day of the Year
             y = my_sum,
             color = factor(year))) +
  geom_point() +
  geom_line() + 
  #scale_y_log10()+
  scale_color_viridis_d() + 
  facet_wrap(.~ art, scales = 'free')
#coord_trans(ytrans="pow10")



  
# May data: time series data:
  
# Plot only time series of IT if there can be seens some seasonality: ----------------
  
dat %>% 
 # filter(month > 4 & year == 2017) %>%             #  &  art == "Buchdrucker"
  #filter(art == 'Buchdrucker') %>% 
  group_by(kont_dat, art)  %>%
    summarize(bav_sum = sum(fangmenge, na.rm = T))   %>%
    ggplot(aes(x = kont_dat,
               y = bav_sum/100000)) +
    #geom_line() +
    geom_point(alpha = 0.5) +
    facet_wrap(art~., scales = 'free')
  
  
  
# ----------------------------------------------
#              ROLLING AVERAGES
# ----------------------------------------------
  
   
  
# Get the rolling averages & and plot population counts by DOY over years --------------
  
  # April 1st  = DOY 91
 dat %>%  
  dplyr::filter(doy > 91) %>% 
  dplyr::arrange(doy) %>% 
    dplyr::group_by(year, art) %>% 
    dplyr::mutate(mean_25     = zoo::rollmean(fangmenge , k = 25, fill = NA)) %>% 
    dplyr::ungroup() %>% 
   ggplot(aes(x = doy, 
              y = mean_25,
              group = factor(year),
              color = factor(year))) + 
   geom_line(alpha = 0.2, color = 'grey20') + 
   geom_smooth() +
   facet_grid(art~., scales = 'free') +
   scale_color_viridis_d()


  

# Some data have counts in Jan-Feb??
dat %>% 
  filter(month < 4) %>% 
  ungroup() %>% 
  distinct(kont_dat) %>%
  pull() 
  print(n = 40)

  
  
# 
summary(dat)

# Join data with XY coordinates by globalid   ---------------------------------------------------

dat <- dat %>% 
  mutate(globalid =  gsub("\\{|\\}", "", globalid))  


# for coordinates table --------------------------------------------------------------------------

# Get the letters to upper cases
xy_sf <-
  xy_sf %>% 
  mutate(globalid = toupper(globalid))



# Export XY and DAT table -------------------------------------------------


# Make a use of gpkg in R: https://inbo.github.io/tutorials/tutorials/spatial_standards_vector/
# Write the output XY files in 3035 coordinate
xy_sf %>% 
  st_write(paste(out_path, "outSpatial/xy_3035.gpkg", sep = '/'), append=FALSE)


data.table::fwrite(dat, paste(out_path, 'outSpatial/dat.csv', sep = "/"))



# Get total sum of buchdrucker over year --------------------------------------------------------
#dat_sum <- 
  dat %>%  
  dplyr::filter(doy > 91  & year != 2014) %>%  # & art == 'Buchdrucker'
  group_by(year, art) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>% 
  ggplot(aes(y = bav_sum,
             x = year,
             color = art)) + 
  geom_line()
  

# get the summary data to merge with XY coordinates!  
dat.sum.IT <-   dat %>%  
  dplyr::filter(doy > 91  & year != 2014 & art == 'Buchdrucker') %>%  # & art == 'Buchdrucker'
  group_by(year, month, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
  mutate(time = paste(year, month, sep = '_')) %>% 
  select(bav_sum, globalid, time) #%>% 
 



# Merge spatial data with sum data by globalid
buch_df <- xy_sf %>% 
  right_join(dat.sum.IT, 'globalid')






# Plot spatial data ------------------------------------

# 
windows()
buch_df %>% 
  filter(year == 2021) %>% 
  ggplot() + 
  geom_sf(aes(color = bav_sum ))  + # , size = 0.5 , size = AREA
  scale_color_continuous(low = "lightgreen", 
                         high = "darkgreen",
                         space = "Lab", 
                         na.value = "red", guide = "colourbar")#+




# R animate: ------------------------------------------------


# Get data:
library(gapminder)

# Charge libraries:
library(ggplot2)
library(gganimate)
library(transformr)
library(gifski)  # needed for correct rendering of the animation

library(rnaturalearth) # for map data
library(ggspatial)




# get spatial data: 
# get and bind country data
de_sf <- ne_states(country = "germany", returnclass = "sf")

# Get only bavaria
bav_sf <- de_sf %>% 
  dplyr::filter(name_en == "Bavaria")


# 
# need to recheck!! gg anime makes a png, but not a video 
p<-ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df,
          aes(color = bav_sum,
              size = bav_sum)) + # , size = Age, size = 0.8size by factor itself!
  viridis::scale_color_viridis(discrete = FALSE, option = "viridis",
                               na.value = "red") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  # gganimate specific bits: ------------------
  labs(title = 'beetle population counts {current_frame}',
       color  = "Beetle [count]") +
  transition_manual(time) +
  #transition_time(year) +
  ease_aes('linear')


# animate with the gifski renderer
animate(p, renderer = gifski_renderer())




# Try if the gganimate actually works: now works widh gifski renderer!!!

library(ggplot2)
library(gganimate)
theme_set(theme_bw())

library(gapminder)
head(gapminder)

# make a plot
p <- ggplot(
  gapminder, 
  aes(x = gdpPercap, y=lifeExp, size = pop, colour = country)
) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  scale_color_viridis_d() +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  labs(x = "GDP per capita", y = "Life expectancy")
p

p<-  p + transition_time(year) +
  labs(title = "Year: {frame_time}")

install.packages("gifski")
library(gifski)

animate(p, renderer = gifski_renderer())
