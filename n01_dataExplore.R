

# -----------------------------------
# Explore bark beetle data
# -----------------------------------

# Notes from data exploratory:

# Individual trap = indicated by falsto_name = region_name + trap 2
# - some traps have A,. B, and C! - do tehy have pairs?
# - globalid - trap position, one trap can have several globalids if moved
# trap_n = 3 exlude! 




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


# Colors libs 
library(RColorBrewer)
library(scales)
library(viridis)


# get data
path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/Fwd__Borkenkäferforschung__Datentransfer'
out_path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria'

# Load RData table
load(paste(path, 'BoMo_2015_2021_Rohdaten.RData', sep = "/"))



# Vars
doy.start  =  91 # April 1st
doy.end    = 304 # Oct 30
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

#str(dat)  # datum is in POSIXct format, ut does not contain time: can convert just to date

# convert to date
dat$kont_dat <- as.Date(dat$kont_dat)

#str(dat)

# Species:
# -----------------------------------
# Buchdrucker  : Ips typographus 
# Kupferstecher: Pityogenes chalcographus   


# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year  = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day   = lubridate::day(kont_dat),
                doy   = lubridate::yday(kont_dat) + 1)  # as POXIT data has January 1st at 0


# get basic stats: count by beetles over year, counts by records/traps
nrow(dat) # >73.000 rows


# # Some basic plotting
# dat %>% 
#   group_by(art, year) %>% 
#   summarise(my_mean = mean(fangmenge, na.rm = T)) %>% 
#   ggplot(aes(x = year,
#              y = my_mean,
#              fill = art,
#              group = year)) +
#   geom_point() +
#   geom_line() +
#  # ylim(0, 250000) +
#   facet_grid(.~art)


# How may traps in total? falsto_name = unique trap name  
length(unique(dat$falsto_name)) # 302 traps

# how many traps per year?
dat %>% 
  filter(art == 'Buchdrucker') %>% 
  group_by(trap_n, year) %>% 
  distinct(falsto_name) %>% 
  count()

# REname falsto_name trap to merge with sf data ------------------------------------------ 
# replace all spaces by '_', get trap_pair number 
dat <- dat %>% 
  mutate(falsto_name2 = gsub(' ', '_', falsto_name)) %>% 
  mutate(trap_pair = as.numeric(str_extract(falsto_name, "[0-9]+"))) # get the trap pair number

# remove the parantheses from globalid to merge with coordinates
dat <- dat %>% 
  mutate(globalid =  gsub("\\{|\\}", "", globalid))  

# how many traps per year, per pair trap?
dat %>% 
  group_by(art, year, trap_pair) %>% 
  distinct(falsto_name, globalid) %>% 
  count() %>% 
  View()

# how to handle if one trap  has several globalids? - just select the first in the row to calculate the LISA
  
# trap_pair == 3 I an exclude, there is only one trap added per year


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


# Get average count of beetles per trap: 

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


ips.year.avg <- dat %>% 
  group_by(art, year, globalid) %>% 
  summarize(sum_beetles = sum(fangmenge, na.rm = T),
            sum_trap    = length(unique(globalid)),
            freq_visit  = length(unique(kont_dat)),
            avg_beetles_trap = sum_beetles/freq_visit) #%>%


# check revisit times:
ips.year.avg %>% 
  ungroup() %>% 
  distinct(freq_visit) %>% 
  pull() %>% 
  hist()

# 6 recrding months: need at least 12 recordings!! or 10..
# 


# about 93 records have 0 sum beetles per year (11 ips, 83 pityogenes), ad checked only once per year! 
# check: ---------------------------------------------------

# Bayerisch Eisenstein C_1 - in 2017
# dat %>% 
#   filter(globalid == "3B99DB06-1520-49E0-9707-7C1916CFFFC3" & art == "Buchdrucker") %>% 
#   View()




# Buchdrucker 2015 {4826784B-6D9B-4CB6-B431-14132184BCB7}
# dat %>% 
#   filter(globalid == "4826784B-6D9B-4CB6-B431-14132184BCB7")
# 
# # or have 0 sum beetles, but checked 22 times!! per season! (Pityogenes)
# #Kupferstecher 2016 {9FD90E97-9C4E-478F-9A89-724C74C95FDC}
# dat %>% 
#   filter(globalid == "9FD90E97-9C4E-478F-9A89-724C74C95FDC" & year == 2016 
#          #& art == "Kupferstecher"
#          ) %>% 
#   arrange(kont_dat)
# 

# get unique site numbers: one globalid can have two traps: one for Ips, one for Pityogenes;

  
# 12 locations of zero beetles for IPS: -----------------------------------------
# ips.year.avg %>% 
#   filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263") %>% 
#   filter(art == 'Buchdrucker')
# 
# dat %>% 
#   filter(art == 'Buchdrucker') %>% 
#   filter(globalid == "9A355BAA-43CC-4004-BFF8-8A3BEF2C1263" & year == 2020)# %>% 
# 
# # buchdrucker has here 23 revisits ad still 0 beetles? - check location
# # Buchdrucker 2016 1557D362-F231-4A47-9976-7063820397C1 0 1 23
# dat %>% 
#   filter(art == 'Buchdrucker') %>% # has zero, check the Brunnthal 2
#   filter(globalid == "1557D362-F231-4A47-9976-7063820397C1" & year == 2016)# %>% 
# 
# 
# dat %>% 
#   filter(art == 'Buchdrucker') %>% # has zero, check the Brunnthal 2: has normal values
#   filter(falsto_name          == "Brunnthal_2" )# %>% 
# 
# 
# # check distribution and zeros:
# ips.year.avg %>% 
#   ggplot(aes(x = avg_beetles_trap)) +
#   geom_histogram(bins = 1000)
# 



# Inspect data ------------------------------------------------------------

# check for zeros presence - present! 
# how many zeros I have for sum beetles per trap?? - if 0 beetles over whole year, it is suspicious..

# filter data:
#  - set dates: April 1st (DOY 91, DOY 92 lap year) to Oct 31  (304) - from Phenips
#  - check revisit times: exclude if traps not collecte reularly! (eg. 10 times over the season: April 1st to Oct 30)

# Filter traps: -------------------------------------------------------------- 
# 1) by 0 sum beetles per year
# 2) by low revisit frequency - once per year

# zero count for whole year, with single recording times... seems fishy
# seems that traps with zeros are consistent...
zero_catch_id <- ips.year.avg %>% 
  filter(sum_beetles == 0 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(globalid) %>% 
  pull()

# low revisit time
low_visit_id <- ips.year.avg %>% 
  filter(freq_visit < 10 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(globalid) %>% 
  pull()



# does the average number of beetles per trap differ between locations?
# between years?
dev.off()

ips.year.avg %>%
  ggplot(aes(x = factor(year),
             y = avg_beetles_trap,
             fill = art)) +
  geom_boxplot() +
  facet_wrap(. ~ art , scales = 'free_y') +
  theme(legend.position = 'bottom')






  
# My data: time series data: -----------------------------------------------------
  
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
# - keep only IPS
# - vegetation period: doy
# - years: 2015-2021
# - frequent recording: 10 revisit times at least
# - remove the '3rd' trap
# - keep only records that have all years and all traps??




dat.ips.clean <- dat %>% 
  filter(doy %in% veg.period) %>% 
  filter(year > 2014) %>% 
  filter(art == 'Buchdrucker') %>%
  ungroup(.) %>% 
  dplyr::filter(!globalid %in% low_visit_id) %>%  # remove traps with low visit time (nee to check this for year?)
  dplyr::filter(!globalid %in% zero_catch_id) %>% # exclude if zero beetles caught per whole year
  dplyr::filter(trap_pair  %in% c(1,2))  #%>%      # exclude if not trap pair (1-2)
  #filter(representativ == 'Ja') 


# filter only traps that have records over all years
dat2 <- dat.ips.clean %>% 
  # group_by(years, loc) %>% 
  group_by(falsto_name) %>%
  dplyr::filter(n_distinct(year) == n_distinct(dat.ips.clean$year))

dat2 %>% 
  group_by(year, trap_pair) %>% 
  distinct(falsto_name) %>% 
  count() #%>% 

# if keeping only traps that are cpontinously monitored over time = only 67*2 traps


  
# Get daily counts per trap (averaged by the revisit time)
df.daily <-  
  dat.ips.clean %>% 
    group_by(globalid,art, year) %>% 
    arrange(kont_dat) %>% # order the data
    mutate(
      period = diff(c(doy.start, doy))+1,  # get difference between two consecutive numbers
      avg_day_count  = fangmenge/period,           # average beetle counts per day
      sum = avg_day_count*period,                  # to check counts
      cumsum = cumsum(sum),
      diff   = avg_day_count - lag(avg_day_count) ) #%>%   # difference between average avg_day_count values
   # filter(globalid == '0010C7C3-C8BC-44D3-8F6E-4445CB8B1DC9' & art == 'Buchdrucker') %>%
   # select(kont_dat, fangmenge, falsto_name, cumsum,avg_day_count,
  #         diff,
   #        doy) %>%
  #  arrange(year) %>%
  #  View()

  
  

 
# check plot of differences: ---------------------------------------
# compare counts by groups:
p_count_diff <- df.daily %>% 
  ggplot(aes(y = diff,
             x = doy,
             color = factor(year))) +
  #geom_point() + 
 #   coord_cartesian(ylim = c(0,500)) + # only '+' values (increase)
  stat_summary(fun = mean, geom="line") + 
  facet_grid(year~.) +
  theme_bw() +
  ylab('Difference in IPS counts')




windows()
hist(df.daily$avg_day_count)

# convert data for counts on 1 and counts on 2: ---------------------
# are 1 and 3 consistently indicating different groups? or is it 2-3, or 1-3?
# 3 is removed (since 08/01/2023)
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


## visualize trap1 vs trap2 ----------------------------------------
ggplot(dd_sum, aes(x = trap1, 
                       y = trap2)) +
  geom_point()+
  geom_smooth(method = 'lm')


m_lm <- lm(trap1~trap2, dd_sum)
summary(m_lm)


m_glm <- glm(trap1~trap2, dd_sum, family = poisson)
summary(m_glm)
plot(m_gam, 1)


# Export XY and DAT table -------------------------------------------------


# Make a use of gpkg in R: https://inbo.github.io/tutorials/tutorials/spatial_standards_vector/
# Write the output XY files in 3035 coordinate
xy_sf %>% 
 st_write(paste(out_path, "outSpatial/xy_3035.gpkg", sep = '/'), append=FALSE)

# export RAW data
data.table::fwrite(dat, paste(out_path, 'outSpatial/dat.csv', sep = "/"))



# get the summary IPS data to merge with XY coordinates!  
ips.year.sum <-   dat.ips.clean %>%     
  group_by(year, art, globalid) %>% 
  mutate(bav_sum = sum(fangmenge, na.rm =T)) %>%  
  dplyr::select(bav_sum, globalid, year) #%>% 


# Merge spatial data with sum data by globalid
ips.year.sum_sf <- xy_sf %>% 
  right_join(ips.year.sum, 'globalid')




# Get spatial data libs -----------------------------------------------------------------
library(rnaturalearth) # for map data
library(ggspatial)


# get spatial data ------------------------------------------------------------------ 
de_sf <- ne_states(country = "germany", returnclass = "sf")

# Get only bavaria
bav_sf <- de_sf %>% 
  dplyr::filter(name_en == "Bavaria")

buch_df_f <- ips.year.sum_sf %>% 
  filter(bav_sum > 3000)

# for visibility, only beetles sums > 3000/year/trap are shown!
ggplot(bav_sf) +
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df_f,
          aes(color = bav_sum,
              size = bav_sum)) + # , size = Age, size = 0.8size by factor itself!
  scale_color_viridis(name = 'Beetle sum/year/trap', 
                      alpha = 0.8,
                      option = 'magma',
                      direction = -1,
                      na.value = 'transparent') +
  facet_wrap(.~year) + 
  theme_void()

# Max increase in counts per year & trap ----------------------------------------------------
# find min DOY of the max diff per year and location: plot on map ----------------
max.diff.doy <- df.daily %>% 
  group_by(globalid, year, falsto_name2) %>% 
  slice(which.max(diff))

max.diff.doy.sf <- xy_sf %>%   # to keep the sf structure, need to add df to sf
  right_join(max.diff.doy, 'globalid')


# check range DOY:
max.diff.doy %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarise(min = min(doy),
            max = max(doy),
            mean = mean(doy),
            sd = sd(doy),
            cv = sd/mean,
            median = median(doy))


# Check Max increase population
max.diff.doy %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  summarise(min = min(diff),
            max = max(diff),
            mean = mean(diff),
            sd = sd(diff),
            cv = sd/mean,
            median = median(diff)
            )


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
p_doy_max_increase <- ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey80') + 
  geom_sf(data = max.diff.doy.sf,
          aes(color = doy,
              size = diff
          )
  ) + 
  scale_color_viridis(name = 'DOY of max increase', 
                      alpha = 0.8,
                      option = 'magma',
                      direction = 1,
                      na.value = 'transparent') +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = 'DOY of max increase/trap/day',
       color  = "DOY",
       size = "Ips [sum]")



## filter  DOY 150 = May 30 ----------------------------------------------------
windows()
p_doy_max_increase150 <- ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey85') + 
  geom_sf(data = filter(max.diff.doy.sf, doy <= 150),
          aes(color = doy,
              size = diff
              )
  ) + 
  scale_color_viridis(name = 'DOY of max increase', 
                      alpha = 0.8,
                      option = 'magma',
                      direction = 1,
                      na.value = 'transparent') +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = 'DOY of max increase/trap/day',
       subtitle = 'DOY <= 150',
       color  = "DOY",
       size = "Ips [sum]")


# does earlier DOY occurs with the higher size?
windows()
p_diff_doy <- max.diff.doy %>% 
  #filter(doy < 151) %>% 
  ggplot(aes(x = doy,
             y = diff,
             color = factor(year)))+
  geom_point() +
  facet_grid(year~.) +
  ylab('Max increase of beetle counts\n[diff]') +
  theme_bw()
  


#DOY
# 100 April 10
# 150 May 30
# 200 July 19
# 211 July 30
# 250 SEpt 7


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
  labs(title = 'Max increase/trap/day',
       #color  = "DOY",
       size = "Ips increase [day/trap]")


# two aspects: the max increase in season
#              the size of population
# driver: previous year population??






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
  summarize(mean = mean(avg_day_count),
            sd  = sd(avg_day_count),
            cv  = sd/mean) %>%  # can be turned to % by '*100'
  ggplot(aes(x = factor(year),
             y = cv)) +
  geom_violin() +
  stat_summary() +
  ggtitle('CV of daily beetle counts')


# CV of DOY of the max beetle increase ------------------
max.diff.doy %>%
#  ungroup(.) %>% 
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







# Save selected dfs in R object: ------------------------------------------------------------
save(
     df.daily,                   # avg daily beetle counts (adjusted by revisit times), over DOY 
     ips.year.avg,                   # sum & avg number per beetles/trap per vegetation season (), over whole year
     ips.year.sum,
     dat.ips.clean,              # filtered raw IPS counts
     max.diff.doy,               # max increase in beetle counts per DOY
     max.diff.doy.sf,            # sf: max increase in beetle counts per DOY
     p_count_diff,               # plot: difference in beetle daily counts  
     p_doy_max_increase,         # map: all locations
     p_doy_max_increase150,      # map: filter early locations
     p_diff_doy,                 # scatter plot per year
     file="outData/ips.Rdata")


# Spatial data for maps:
save(xy_sf,                      # trap location (one Globalid is ised for both IPS & Chalcographus!)
     de_sf,                      # germany shp
     bav_sf,                     # bavaria shp
     file="outData/spatial.Rdata") 
