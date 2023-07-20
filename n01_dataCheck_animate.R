
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
                doy   = lubridate::yday(kont_dat) + 1)  # as POXIT data has January 1st at 0


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
dat <- dat %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) 

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

# usefull period: 2015-2021 (2014 does not have counts per trap, only sum by year)
#   art            year     n
#   <chr>         <dbl> <int>
# 1 Buchdrucker    2014     1
# 2 Buchdrucker    2015   203
# 3 Buchdrucker    2016   234
# 4 Buchdrucker    2017   224
# 5 Buchdrucker    2018   247
# 6 Buchdrucker    2019   244
# 7 Buchdrucker    2020   251
# 8 Buchdrucker    2021   246
# 9  Kupferstecher  2014     1
# 10 Kupferstecher  2015   203
# 11 Kupferstecher  2016   234
# 12 Kupferstecher  2017   224
# 13 Kupferstecher  2018   247
# 14 Kupferstecher  2019   244
# 15 Kupferstecher  2020   251
# 16 Kupferstecher  2021   246


# --------------------------------------------------------------------
#              Standardize beetle counts/trap/year
# --------------------------------------------------------------------


# how to get the average number of beetles per trap per year?

# Calculate sum of beetles/year
# sum of traps/year
# divide by frequency of recording: revisit time? : 
# how many times the traps have been revisited? 'kont_dat': lenght(unique)..



#     art            year sum_beetle freq_visit
#   <chr>         <dbl>      <int>      <int>
# 1 Buchdrucker    2015      12628         26
# 2 Buchdrucker    2016      14907         23
# 3 Kupferstecher  2015      19632         26
# 4 Kupferstecher  2016      14506         23


# Calculate the mean number of beetles per trap over the season
# corrected for revisit frequency
dat_avg <- dat %>% 
  group_by(art, year, globalid) %>% 
  summarize(sum_beetles = sum(fangmenge, na.rm = T),
            sum_trap    = length(unique(globalid)),
            freq_visit  = length(unique(kont_dat)),
            avg_beetles_trap = sum_beetles/freq_visit) #%>%


# check revisit times:
dat_avg %>% 
  ungroup() %>% 
  distinct(freq_visit) %>% 
  pull() %>% 
  hist()

# 6 recrding months: need at least 12 recordings!! or 10..
# 


# about 93 records have 0 sum beetles per year (11 ips, 83 pityogenes), ad checked only once per year! 
# check:
# Buchdrucker 2015 {4826784B-6D9B-4CB6-B431-14132184BCB7}
dat %>% 
  filter(globalid == "{4826784B-6D9B-4CB6-B431-14132184BCB7}")

# or have 0 sum beetles, but checked 22 times!! per season! (Pityogenes)
#Kupferstecher 2016 {9FD90E97-9C4E-478F-9A89-724C74C95FDC}
dat %>% 
  filter(globalid == "{9FD90E97-9C4E-478F-9A89-724C74C95FDC}" & year == 2016 
         #& art == "Kupferstecher"
         ) %>% 
  arrange(kont_dat)


# get unique site numbers: one globalid can have two traps: one for Ips, one for Pityogenes;
# seems that traps with zeros are consistent...
zero_ips_traps <- dat_avg %>% 
  filter(sum_beetles == 0 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(globalid) %>% 
  pull()
  
# 12 locations of zero beetles for IPS:
dat_avg %>% 
  filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263") %>% 
  # dplyr::filter(globalid %in% zero_ips_traps) #%>% 
  filter(art == 'Buchdrucker')

dat %>% 
  filter(art == 'Buchdrucker') %>% 
  filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263" & year == 2020)# %>% 
  

# check distribution and zeros:
dat_avg %>% 
  ggplot(aes(x = avg_beetles_trap)) +
  geom_histogram(bins = 500)


# how many zeros I have for mean beetles per trap??
dat_avg %>% 
  filter(avg_beetles_trap == 0) %>% 
  View()

# what is average revisit time?
# number of days: April to October: 212 days
range(212/dat_avg$freq_visit)

hist(212/dat_avg$freq_visit)

# some sites have been revisited every 7 days!
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
  facet_wrap(. ~ art , scales = 'free_y') +
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
# Correctly order levels
dat_avg <- dat_avg %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) 



dat_avg <- dat_avg %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))
                                 
     


# geom boxplot plot
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
  
# Plot only time series of IT if there can be seems some seasonality: ----------------
  
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
  
  
  
  
  
  
# -------------------------------------------------------------
# variations between traps???
# ---------------------------------------------------------
  
# spliut in two groups: 1&2
head(dat)
  
# 
sort(unique(dat$falsto_name2))

# First, replace all spaces by '_'
dat <- dat %>% 
  mutate(falsto_name2 = gsub(' ', '_', falsto_name)) %>% 
  mutate(pair_grp = as.numeric(gsub("\\D", "", falsto_name2)))


# compare counts by groups:
dat %>% 
  filter(pair_grp != 3 & year > 2014) %>% 
  ggplot(aes(y = fangmenge/10000,
             x = factor(year),
             fill = factor(pair_grp))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  facet_grid(.~art)
  

# convert data for counts on 1 and counts on 2:
# in yx format:

dat_pairs <- dat %>% 
  filter(pair_grp != 3 & year > 2014) %>% 
  dplyr::select('year', 'fangmenge', 'pair_grp', 'art', 'monsto_name')


# split in two grousp
dat_pairs1 <- dat_pairs %>% 
  filter(pair_grp == 1) %>% 
  rename(fangmenge1 = fangmenge)

dat_pairs2 <- dat_pairs %>% 
  filter(pair_grp == 2) %>% 
  rename(fangmenge2 = fangmenge)

# join data tables
dat_pairs_out <- dat_pairs1 %>% 
  left_join(dat_pairs2, by=c('year', 'art', 'monsto_name'))
  
  
# variation between catched numbers:

dat %>% 
  filter(pair_grp != 3 & year > 2014) %>% 
  ggplot(aes(y = fangmenge/10000,
             x = factor(year),
             fill = factor(pair_grp))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  facet_grid(.~art)


# Check if individual captures correlate between paired locations?
windows()
dat_pairs_out %>% 
  filter(art == 'Buchdrucker') %>%  
  ggplot(aes(y = fangmenge1,
             x = fangmenge2,
             color = art)) +
  geom_point() #+
  scale_y_continuous(trans='log10') +
  facet_grid(.~art)


# seems weird: just check on single locations?
  dat %>% 
    filter(monsto_name == 'Hemau' & year == 2021)
  
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
dat_sum <- 
  dat %>%  
  dplyr::filter(doy > 91  & year != 2014) %>%  # & art == 'Buchdrucker'
  group_by(year, art) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T))# %>% 
  #ggplot(aes(y = bav_sum,
  #           x = year,
  #           color = art)) + 
  #geom_line()
  

# get the summary data to merge with XY coordinates!  
dat.sum.IT <-   dat %>%  
  dplyr::filter(doy > 91  & year != 2014 & art == 'Buchdrucker') %>%  # & art == 'Buchdrucker'
  group_by(year, month, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
  mutate(time = paste(year, month, sep = '_')) %>% 
  dplyr::select(bav_sum, globalid, time) #%>% 
 



# Merge spatial data with sum data by globalid
buch_df <- xy_sf %>% 
  right_join(dat.sum.IT, 'globalid')






# Plot spatial data ------------------------------------
library(RColorBrewer )
display.brewer.pal(7, "BrBG")

buch_df_f <- buch_df %>% 
  filter(bav_sum > 3000)
# 
windows()
#buch_df %>% 
  #filter(year == 2021) %>% 
ggplot(bav_sf) +
  geom_sf(color = 'black', 
        fill  = 'grey93') + 
  geom_sf(data = buch_df_f,
          aes(color = bav_sum,
              size = bav_sum)) + # , size = Age, size = 0.8size by factor itself!
  #geom_sf(buch_df, aes(color = bav_sum ))  + # , size = 0.5 , size = AREA
  #scale_color_continuous(low = "#FFFF00",
                       #  high = "#990000",
                       #  space = "Lab",
                       #  na.value = "transparent",
                       # guide = "colourbar") +
  scale_color_viridis(name = 'Beetle count/trap', 
                      alpha = 0.8,
                     option = 'magma',
                     direction = -1,
                     na.value = 'transparent') +
  facet_wrap(.~year) + 
  theme_bw()




# try to plot variogram:  ----------------------------------
library(gstat)

#buch_df2 <- data.frame(buch_df) %>% 
  

?variogram
vgm1<- variogram(log(bav_sum)~1,buch_df)
plot(vgm1, type='b', main='Co-variogram')





# Test fr Morans'I: ---------------------------
Moran.I(ozone$Av8top, ozone.dists.inv) 




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





# ===============================================================================
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
