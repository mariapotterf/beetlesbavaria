
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


# Vars
doy.start  =  91
doy.end    = 304
veg.period = doy.start:doy.end

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

# Get the letters to upper cases - to merge with beetle counts
xy_sf <-
  xy_sf %>% 
  mutate(globalid = toupper(globalid))



# Get beetle count data -------------------------------------------------------------- 

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
nrow(dat) # >73.000 rows


# Some basic plotting
dat %>% 
  group_by(art, year) %>% 
  summarise(my_mean = mean(fangmenge, na.rm = T)) %>% 
  ggplot(aes(x = year,
             y = my_mean,
             fill = art,
             group = year)) +
  geom_point() +
  geom_line() +
 # ylim(0, 250000) +
  facet_grid(.~art)


# Split data between years: , until 2018, > 2018 - drought years, 
# need to standardize data between years: get mean yearly counts?
# maybe also standarziye by number of plots??? does number od plot vary over years?
dat <- dat %>% 
  mutate(drought_period = case_when(
    year < 2018 ~  'before',
    year >= 2018 ~ 'after')) 

# Correctly order levels
dat <- dat %>% 
  mutate(drought_period = factor(drought_period,
                                 levels=c('before', 'after')))




# How may traps in total?  
length(unique(dat$globalid)) # 490 traps

# how many traps per year?
dat %>% 
  group_by(art, year) %>% 
  distinct(globalid) %>% 
  count()

# First, replace all spaces by '_'
dat <- dat %>% 
  mutate(falsto_name2 = gsub(' ', '_', falsto_name)) %>% 
  mutate(trap_pair = as.numeric(str_extract(falsto_name, "[0-9]+"))) # get the trap pair number

# remove the parantheses
dat <- dat %>% 
  mutate(globalid =  gsub("\\{|\\}", "", globalid))  

 
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



# Standardize beetle counts/trap/year --------------------------------------------------------------------


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


# Calculate the mean number of beetles per trap over the season ----------------
dat_avg <- dat %>% 
  group_by(art, year, globalid,drought_period ) %>% 
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
# check: ---------------------------------------------------
# Buchdrucker 2015 {4826784B-6D9B-4CB6-B431-14132184BCB7}
dat %>% 
  filter(globalid == "4826784B-6D9B-4CB6-B431-14132184BCB7")

# or have 0 sum beetles, but checked 22 times!! per season! (Pityogenes)
#Kupferstecher 2016 {9FD90E97-9C4E-478F-9A89-724C74C95FDC}
dat %>% 
  filter(globalid == "9FD90E97-9C4E-478F-9A89-724C74C95FDC" & year == 2016 
         #& art == "Kupferstecher"
         ) %>% 
  arrange(kont_dat)


# get unique site numbers: one globalid can have two traps: one for Ips, one for Pityogenes;
# seems that traps with zeros are consistent...
zero_catch_id <- dat_avg %>% 
  filter(sum_beetles == 0 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(globalid) %>% 
  pull()

# low revisit time
low_visit_id <- dat_avg %>% 
  filter(freq_visit < 10 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(globalid) %>% 
  pull()
  
# 12 locations of zero beetles for IPS:
dat_avg %>% 
  filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263") %>% 
  # dplyr::filter(globalid %in% zero_catch_id) #%>% 
  filter(art == 'Buchdrucker')

dat %>% 
  filter(art == 'Buchdrucker') %>% 
  filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263" & year == 2020)# %>% 

# buchdrucker has here 23 revisits ad still 0 beetles? - check location
# Buchdrucker 2016 1557D362-F231-4A47-9976-7063820397C1 0 1 23
dat %>% 
  filter(art == 'Buchdrucker') %>% # has zero, check the Brunnthal 2
  filter(globalid == "1557D362-F231-4A47-9976-7063820397C1" & year == 2016)# %>% 


dat %>% 
  filter(art == 'Buchdrucker') %>% # has zero, check the Brunnthal 2: has normal values
  filter(falsto_name          == "Brunnthal_2" )# %>% 


# check distribution and zeros:
dat_avg %>% 
  ggplot(aes(x = avg_beetles_trap)) +
  geom_histogram(bins = 1000)



# Inspect data ------------------------------------------------------------

# check for zeros presence - present! 
# how many zeros I have for sum beetles per trap?? - if 0 beetles over whole year, it is suspicious..

# filter data:
#  - set dates: April 1st (DOY 91, DOY 92 lap year) to Oct 31  (304) - from Phenips
#  - check revisit times: exclude if traps not collecte reularly! (eg. 10 times over the season: April 1st to Oct 30)


# what is average revisit time?
# number of days: April to October: 212 days
range(212/dat_avg$freq_visit)

hist(212/dat_avg$freq_visit)

# some sites have been revisited every 7 days!
# calculation works great!



# does the average number of beetles per trap differ between locations?
# between years?
dat_avg %>%
  ggplot(aes(x = factor(year),
             y = avg_beetles_trap,
             fill = art)) +
  geom_boxplot() +
  facet_wrap(. ~ art , scales = 'free_y') +
  theme(legend.position = 'bottom')







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
  filter(doy %in% veg.period ) %>%             #  &  art == "Buchdrucker"
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
  filter(doy %in% veg.period ) %>%           #  &  art == "Buchdrucker"
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
  filter(doy %in% veg.period) %>%             #  &  art == "Buchdrucker"
  filter(year == 2017 ) %>% 
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
  group_by(kont_dat, art)  %>%
    summarize(bav_sum = sum(fangmenge, na.rm = T))   %>%
    ggplot(aes(x = kont_dat,
               y = bav_sum/100000)) +
    #geom_line() +
    geom_point(alpha = 0.5) +
    facet_wrap(art~., scales = 'free')
  
  
  
  


# Clean data: IPS ---------------------------------------------------------
# - keep IPS
# - vegetation period: doy
# - years
# - frequent recording: 10 revisit times at least
# - remove the '3' groups???
dat.ips.clean <- dat %>% 
  filter(doy %in% veg.period) %>% 
  filter(year > 2014) %>% 
  filter(art == 'Buchdrucker') %>%
  #mutate(globalid = gsub(paste(c("[{]", "[}]"), collapse = "|"), "", globalid) ) %>% # remove parantheses 
  ungroup(.) %>% 
  dplyr::filter(!globalid %in% low_visit_id) %>%  # remove traps with low visit time (nee to check this for year?)
  dplyr::filter(!globalid %in% zero_catch_id) %>% # exclude if zero beetles caught per whole year
  dplyr::filter(trap_pair  %in% c(1,2))  %>%      # exclude if not trap pair (1-2)
  filter(representativ == 'Ja') 
  
  
# Get daily counts
df.daily <-  
  dat.ips.clean %>% 
    group_by(globalid,art, year) %>% 
    arrange(kont_dat) %>% # order the data
    mutate(
      period = diff(c(doy.start, doy))+1,  # get difference between two consecutive numbers
      daily  = fangmenge/period,           # average beetle counts per day
      sum = daily*period,                  # to check counts
      cumsum = cumsum(sum),
      diff   = daily - lag(daily) ) #%>%   # difference between average daily values
   # filter(globalid == '0010C7C3-C8BC-44D3-8F6E-4445CB8B1DC9' & art == 'Buchdrucker') %>%
   # select(kont_dat, fangmenge, falsto_name, cumsum,daily,
  #         diff,
   #        doy) %>%
  #  arrange(year) %>%
  #  View()

  
  

 
##### example -----------------------------------------------------

# https://stackoverflow.com/questions/76806372/count-number-caught-per-day-from-longer-intervals


  df <- data.frame(gr = rep(c('a','b'), each = 3),
                   count = c(10,15,23, 
                             5, 10, 100),
                   doy = c(6, 11,15,
                           7,12,19))

# example base R
  first_day <- 5
  df$period <- diff(c(first_day, df$doy))+1
  df$daily <- df$count/df$period  

# Example tidyr with groups
  df %>% 
    group_by(gr) %>% 
    mutate(
           period = diff(c(first_day, doy))+1,
           daily  = count/period,
           sum = daily*period,
           cumsum = cumsum(sum))
    
    
 
  
  
# check plot of differences: ---------------------------------------
# compare counts by groups:
  df.daily %>% 
  ggplot(aes(y = diff,
             x = doy,
             color = year)) +
  #geom_point() + 
    coord_cartesian(ylim = c(0,500)) + # only '+' values (increase)
  stat_summary(fun = mean, geom="line") + 
  facet_grid(year~.)



length(unique(df.daily$globalid))

windows()
hist(df.daily$daily)

filter(df.daily)

# convert data for counts on 1 and counts on 2: ---------------------
# are 1 and 3 consistently indicating different groups? or is it 2-3, or 1-3?
# quick overview: subset only data that have coorrect 1/2:
df1 <-dat.ips.clean %>% 
  filter(trap_pair == 1 ) %>% 
  dplyr::select(monsto_name, year, month, doy, fangmenge, art)

df2 <-dat.ips.clean %>% 
  filter(trap_pair == 2 ) %>% 
  dplyr::select(monsto_name, year, month, doy, fangmenge, art)

# merge df1 & df2
dd <- df1 %>% 
  left_join(df2, by = c("monsto_name", "year", "month", "doy", "art"))


# aggeraget by months instead of teh DOY:
dd_sum <- 
  dd %>% 
  group_by(year, month, monsto_name) %>% 
  summarise(trap1 = sum(fangmenge.x),
            trap2 = sum(fangmenge.y)) %>%
  ungroup() %>%
  filter(complete.cases(.))

# od Laca: ----------
# tau autokoralacia - korelacia v case
# pozriet differences v casoch
# elevation?


ggplot(dd_sum, aes(x = trap1, 
                       y = trap2)) +
  geom_point()+
  geom_smooth(method = 'gam')


  




# Export XY and DAT table -------------------------------------------------


# Make a use of gpkg in R: https://inbo.github.io/tutorials/tutorials/spatial_standards_vector/
# Write the output XY files in 3035 coordinate
xy_sf %>% 
  st_write(paste(out_path, "outSpatial/xy_3035.gpkg", sep = '/'), append=FALSE)


data.table::fwrite(dat, paste(out_path, 'outSpatial/dat.csv', sep = "/"))



# Get total sum of buchdrucker over year --------------------------------------------------------
dat_sum <- 
  dat %>%  
  dplyr::filter(doy >= doy.start &  doy <= doy.end & # April 1st to Oct 30
            year != 2014) %>%  # & art == 'Buchdrucker'
  group_by(year, art) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T))# %>% 
  #ggplot(aes(y = bav_sum,
  #           x = year,
  #           color = art)) + 
  #geom_line()
  

# get the summary data to merge with XY coordinates!  
dat.sum.IT <-   dat %>%  
  dplyr::filter(doy >= doy.start &  doy <= doy.end  & year != 2014 & art == 'Buchdrucker') %>%  # & art == 'Buchdrucker'
  group_by(year, month, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
  mutate(time = paste(year, month, sep = '_')) %>% 
  dplyr::select(bav_sum, globalid, time) #%>% 
 


dat.sum.IT.year <-   dat %>%  
  dplyr::filter(doy >= doy.start &  doy <= doy.end & year != 2014 & art == 'Buchdrucker') %>%  # & art == 'Buchdrucker'
  group_by(year, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
 # mutate(time = paste(year, month, sep = '_')) %>% 
  dplyr::select(bav_sum, globalid, year) #%>% 



# Merge spatial data with sum data by globalid
buch_df <- xy_sf %>% 
  right_join(dat.sum.IT, 'globalid')

buch_df_year <- xy_sf %>% 
  right_join(dat.sum.IT.year, 'globalid')






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






# Max increase in counts per year & gap ----------------------------------------------------

# find min DOY of the max diff per year and location: plot on map ----------------
max.diff.doy <- df.daily %>% 
  group_by(globalid, year, falsto_name2) %>% 
  slice(which.max(diff))

max.diff.doy.sf <- xy_sf %>% 
  right_join(max.diff.doy, 'globalid')


# check range DOY:
max.diff.doy.sf %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarise(min = min(doy),
            max = max(doy),
            mean = mean(doy),
            sd = sd(doy),
            median = median(doy))


# Check Max increase population
max.diff.doy.sf %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarise(min = min(diff),
            max = max(diff),
            mean = mean(diff),
            sd = sd(diff),
            median = median(diff))


# the main diference in counts happends earlier in the season
avg.doy <- mean(df.daily$doy)

windows()
df.daily %>% 
  ggplot(aes(y = doy,
             x = factor(year))) +
  geom_violin() +
  stat_summary() +
  geom_hline(yintercept = avg.doy, lty = 'dashed', col = "red")
  #stat_summary(fun = mean, geom="line")



# facet map plots -------------------------------------------------------------------
## overall: DOY f max increase
ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey33') + 
  geom_sf(data = max.diff.doy.sf,
          aes(color = doy#,
              #size = diff
          )
  ) + # , size = Age, size = 0.8size by factor itself!
  colorspace::scale_color_continuous_sequential(palette = "Heat", alpha = 0.8,
                                                rev = FALSE) + # reverse order  = sooner = darker color
  # annotation_scale(location = "bl", 
  #                   width_hint = 0.4) +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = 'DOY of max increase/trap/day',
       color  = "DOY",
       size = "Ips [sum]")



## filter  DOY 150 = July 31 ----------------------------------------------------
ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey33') + 
  geom_sf(data = filter(max.diff.doy.sf, doy <= 150),
          aes(color = doy#,
              #size = diff
              )
  ) + # , size = Age, size = 0.8size by factor itself!
  colorspace::scale_color_continuous_sequential(palette = "Heat", alpha = 0.8,
                                                rev = FALSE) + # reverse order  = sooner = darker color
  # annotation_scale(location = "bl", 
  #                   width_hint = 0.4) +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = 'DOY of max increase/trap/day',
       color  = "DOY",
       size = "Ips [sum]")


# 2.population increase -------------------------------------------------------------
ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey33') + 
  geom_sf(data = max.diff.doy.sf,
          aes(color = diff#,
              #size = diff
          )
  ) + # , size = Age, size = 0.8size by factor itself!
  colorspace::scale_color_continuous_sequential(palette = "Heat", alpha = 0.8,
                                                rev = TRUE) + # reverse order  = sooner = darker color
  # annotation_scale(location = "bl", 
  #                   width_hint = 0.4) +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = 'Max population size/trap/day',
       #color  = "DOY",
       size = "Ips increase [day/trap]")


# two aspects: the max increase in season
#              the size of population
# driver: previous year population??





# ANOVA : # does the DOY of max increase differes between years?  ----------------
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r#check-your-data
# 
# one way anova - finds that there are differenes, we don't know where
res.aov <- aov(doy ~ factor(year), max.diff.doy.sf)

summary(res.aov)

# multiple way anova = pairwise comparison between groups
TukeyHSD(res.aov) # challengin to read differences, if many groups

pairwise.t.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                p.adjust.method = "BH")


### check ANOVA assumptions:
# 1. Homogeneity of variances
# 2. Normality


# 1. Homogeneity of variances
plot(res.aov, 1)

# test for homonenity f variance: 
library(car)
# high p-value = variances are the same (H0) 
leveneTest(doy ~ factor(year), data = max.diff.doy.sf)
# low p-values, so variancs are different ->

# use ANOVA with no assumption of equal variances
res.aov_no_equal <- oneway.test(doy ~ factor(year), data = max.diff.doy.sf)



### t-test for not equal variances : ---------------------------
pairwise.t.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                p.adjust.method = "BH", pool.sd = FALSE)


# 2. Normality
plot(res.aov, 2)

# check by Shapiro-Wilk test for residuals
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test # does not work, have more values then 5000
shapiro.test(x = aov_residuals )  # not normaly distributed


# run rather Kruskal-Wallis test - not equal variances
kruskal.test(doy ~ factor(year), data = max.diff.doy.sf)


# post-hoc test: pairwise-wilcox
wilk <-pairwise.wilcox.test(max.diff.doy.sf$doy, factor(max.diff.doy.sf$year),
                     p.adjust.method = "BH")
# get matrix of p-values for further ploting
wilk$p.value

p.vals <- replace(wilk$p.value,wilk$p.value>0.05,1)


# show p-values - first transformation is needed, 
# as this is the opposite of the corr() plot - larger is better,
# while in p-values smaller is better
# https://stackoverflow.com/questions/61861049/how-plot-results-from-a-pairwise-wilcox-test
my_transformed_pvals=-log10(p.vals)  # negative of 'log' works, as all p-vals are <1

library ( corrplot)
corrplot(as.matrix(my_transformed_pvals),is.corr=F)
corrplot(as.matrix(p.vals),is.corr=F)


# !!! my groups are not independent! they are colected on teh same traps, 
# they are autocorrelated!
# use anova with repeated measures:

aov_rep <- aov(doy~factor(year)+Error(factor(globalid)), data = max.diff.doy.sf)
summary(aov_rep)



# Coefficient of variation -------------------------------------------------------
# EXAMPLE : Sample beetle count data (mean counts per year)
years <- c(2017, 2018, 2019, 2020, 2021)
mean_counts <- c(50, 60, 55, 75, 70)

# Sample standard deviation data for each year
std_dev <- c(10, 8, 12, 15, 10)

# Calculate coefficient of variation (CV) for each year
cv_values <- (std_dev / mean_counts) * 100

# Create a plot to illustrate temporal variability using CV
windows()
plot(years, cv_values, type='o', col='blue', pch=16, lty=1, xlab='Year', ylab='Coefficient of Variation (CV %)',
     main='Temporal Variability of Beetle Counts per Year')
grid()

# CV for daily counts data ---------------------------------------------------------------
# for whole data, not just the max increase
df.daily %>%
  group_by(year, globalid, art) %>% 
  summarize(mean = mean(daily),
            sd  = sd(daily),
            cv  = sd/mean) %>%  # can be turned to % by '*100'
  ggplot(aes(x = factor(year),
             y = cv)) +
  geom_violin() +
  stat_summary() +
  ggtitle('CV of daily beetle counts')


# CV of DOY of the max beetle increase ------------------
max.diff.doy %>%
  ungroup(.) %>% 
  group_by(year) %>% 
  summarize(mean = mean(doy),
            sd   = sd(doy, na.rm=TRUE) ,
            cv   = sd/mean) %>%  # can be turned to % by '*100'
  ggplot(aes(x = year,
             y = cv)) +
  geom_line() +
 # stat_summary() +
  ggtitle('CV of DOY of max beetle increase')

  

# Interpret CV: smaller CV is better
# 1/10 (sd/mean)  = 0.1
# 5/10 (sd/mean)  = 0.5







# GIf by years --------------------------------------------------------

# 
# need to recheck!! gg anime makes a png, but not a video 
p<-ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df_year,
          aes(color = bav_sum,
              size = bav_sum)) + # , size = Age, size = 0.8size by factor itself!
  # viridis::scale_color_viridis(discrete = FALSE, option = "viridis",
  #                             na.value = "red") +
  #scale_color_gradientn(colours = #colorspace::hcl_palettes(palette = "Dark 2")
   #                      colorspace::heat_hcl(7)
    #                   ) +
  colorspace::scale_color_continuous_sequential(palette = "Heat") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  # gganimate specific bits: ------------------
labs(title = 'IT beetle population counts {current_frame}',
     color  = "IT Beetle [count]") +
  transition_manual(year) +
  
  ease_aes('linear')


# animate with the gifski renderer
animate(p, renderer = gifski_renderer())


hcl_palettes("sequential (single-hue)", n = 7, plot = TRUE)















# GIF: months and years -------------------------------------------------------- 
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
