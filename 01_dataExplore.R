


# Explore bark beetle data


# Notes from data exploratory:

# Individual trap = indicated by falsto_name = region_name + trap 2
# - some traps have A,. B, and C! - do tehy have pairs?
# - globalid - trap position, one trap can have several globalids if moved
# trap_n = 3 exlude! 

# check in how many traps beetle start to swarm early?


# Input ------------------------------------

rm(list=ls()) 

### Libs -----------------------------------------------------------
library(lubridate)
library(dplyr)
library(tidyr)
library(sf)
library(stringr)
library(ggplot2)
library(zoo)
library(data.table)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# Colors libs 
library(RColorBrewer)
library(scales)
library(viridis)


### get data -----------------------------------
# Set the root path dynamically (relative to the R project folder)
root_path <- here::here()

# Define relative paths for data and output
raw_data_path <- file.path(root_path, "rawData", "Fwd__Borkenkäferforschung__Datentransfer")
output_path <- file.path(root_path)

source('myVars.R')

#path = 'rawData/Fwd__Borkenkäferforschung__Datentransfer'
#out_path ='outSpatial' #'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria'

# Load RData table
load(file.path(raw_data_path, "BoMo_2015_2021_Rohdaten.RData"))


## Read XY coordinates -----------------------------------------------------------
# need to save the csv fil as a new file, otherwise it did not work
# or use teh  version 3: now it worked with the original one! 

xy <- read.csv2(paste(raw_data_path, "BoMo_Fallenstandorte_Lage_Laufzeit.csv", sep = '/'), 
                dec = ',', header = TRUE, quote = "\"", 
                skip = 1, fileEncoding = "UTF-16LE")


#readLines(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit.csv", sep = '/'))
# Investigate beetle data: Bavaria
head(xy)


# Convert coordinates to sf object
xy_sf <- st_as_sf(xy, coords = c('x_coord', 'y_coord' ), 
                  crs = 32632 )  # What is a coordinate system??? for forestry in bavaria?

plot(xy_sf["OBJECTID"])
#st_distance(xy_sf, by_element = TRUE)

# Change projection to common european one: 3035
xy_sf <- st_transform(xy_sf, crs = 3035)


# Change global ID and falsto names to merge them with beetle counts 
# Get the GLOBALID letters to upper cases - to merge with beetle counts
xy_sf <-
  xy_sf %>% 
  mutate(globalid = toupper(globalid),
         falsto_name = gsub(' ', '_', falsto_name)) %>%
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) #%>% #remove special characters from falsto_name 
  

sort(unique(xy_sf$falsto_name))

# export the shp of teh trap files
xy_sf %>%  
  st_write(paste(output_path, "outSpatial/all_traps_3035.gpkg", sep = '/'), append=FALSE)


# Beetle data processing -------------------------------------------------------------- 

# Change name of the table
# for dat
dat <- Daten_B01

#length(unique(dat$falsto_name)) # 302 regions

# "Buchdrucker"   "Kupferstecher"

# convert to date
dat$kont_dat <- as.Date(dat$kont_dat)

#str(dat)

# Species: 
# Buchdrucker  : Ips typographus 
# Kupferstecher: Pityogenes chalcographus   


# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year  = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day   = lubridate::day(kont_dat),
                doy   = lubridate::yday(kont_dat) + 1)  #%>% # as POXIT data has January 1st at 0
   
#filter(dat, year == 2014)  # I do not have trap names neither globalid for 2014 year
# get basic stats: count by beetles over year, counts by records/traps
nrow(dat) # >73.000 rows

# How may traps in total? falsto_name = unique trap name  
length(unique(dat$falsto_name)) # 302 traps


# REname falsto_name trap to merge with sf data 
# replace all spaces by '_', get trap_n number 
dat <- dat %>% 
   mutate(falsto_name = gsub(' ', '_', falsto_name)) %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>%  #remove special characters from falsto_name 
  mutate(trap_n = as.numeric(str_extract(falsto_name, "[0-9]+"))) # get the trap pair number

# remove the parantheses from globalid to merge with coordinates
dat <- dat %>% 
  mutate(globalid =  gsub("\\{|\\}", "", globalid))  


dat %>% 
  filter(art == "Buchdrucker") %>% 
  group_by(falsto_name) %>% 
  distinct(globalid) %>%
  dplyr::summarize(globalid_count = n()) %>%
  filter(globalid_count> 1)
  #arrange(-globalid_count) 

# 138 traps from 302 have changed location 2-4 times over 2015-2021
dat %>% 
  filter(falsto_name == "Falkenstein_1") %>% #   'Aldersbach_1') %>% Falkenstein_1 has trap 4x changed
  filter(art == "Buchdrucker") %>% 
  group_by(art, falsto_name, year) %>% 
  distinct(globalid) %>%
  arrange(year)
  
# trap_n == 3 I an exclude, there is only one trap added per year


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


# get median dates between two recordings
# Calculate the number of days between recordings per 'falsto_name'
median_days_between <- dat %>%
  dplyr::filter(art == 'Buchdrucker') %>% 
  arrange(falsto_name, kont_dat) %>%  # Sort data by 'falsto_name' and 'kont_dat'
  group_by(falsto_name, year) %>%           # Group by 'falsto_name'
  mutate(days_diff = as.numeric(difftime(kont_dat, lag(kont_dat), units = "days"))) %>%  # Calculate days difference
  summarise(median_days = median(days_diff, na.rm = TRUE),
            mean_days = mean(days_diff, na.rm = TRUE),
            sd_days = sd(days_diff, na.rm = TRUE))  # Calculate median days difference, excluding NAs

(median_days_between)
hist(median_days_between$median_days)
range(median_days_between$median_days, na.rm  =T)
range(median_days_between$mean_days, na.rm  =T)

# how many records have median of <7?
# Calculate the number of median_days values that are less than 7
med_less_than_7 <- sum(median_days_between$median_days <= 7, na.rm = TRUE)

# Print the result
print(med_less_than_7)


# get raw data with DOY counts -

dat_DOY_count <- dat %>%
  dplyr::filter(art == 'Buchdrucker') %>%
  dplyr::filter(year !=2014) %>% 
  dplyr::filter(month %in% 4:9) %>%  # filter months during vegetation season: some has record in januay, march ..
  dplyr::select(
    year,
    doy,
    #month,
    fangmenge,
    art,
    falsto_name,
    monsto_name,
    trap_n
  ) %>%
  dplyr::rename(
    count = fangmenge,
    species = art,
    trapID = falsto_name,
    pairID = monsto_name
  )



## Standardize beetle counts/trap/year --------------------------------------------------------------------


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


ips.year.avg <- 
  dat %>% 
  dplyr::filter(year > 2014) %>% 
  group_by(art, year, falsto_name ) %>% 
  dplyr::summarize(sum_beetles = sum(fangmenge, na.rm = T),
            freq_visit  = length(unique(kont_dat)),
            avg_beetles_trap = sum_beetles/freq_visit) #%>%


# coount by weeks:
ips.week.agg <- dat %>% 
  dplyr::filter(year > 2014) %>% 
  dplyr::filter(art == "Buchdrucker") %>% 
  mutate(week_number = week(kont_dat)) %>% 
    group_by(art, year,week_number, falsto_name ) %>% 
  dplyr::summarize(week_sum_beetles = sum(fangmenge, na.rm = T)) #%>%

# gte the maximal beetle population groth per week: - can i compare with my data???
# no necessary: I have used teh daily values to compensate for different catching interval between the data, as they are not collected regularly
View(ips.week.agg)


ips.week.agg %>% 
  dplyr::filter(week_number %in% 10:40) %>% 
  ggplot(aes(x = week_number,
             y = week_sum_beetles,
             group = falsto_name)) + 
  geom_line(alpha = 0.1) + 
  geom_hline(yintercept = 3000) + 
  facet_wrap(~year, scales = 'free')


# check revisit times:
ips.year.avg %>% 
  ungroup() %>% 
  distinct(freq_visit) %>%
  pull() %>% 
  hist()




# 6 recrding months: need at least 12 recordings!! or 10..
# 

  # Check for 0 beetles:
  
  #  # need to check! "Brunnthal_1"   "Weißenstadt_1" "Weißenstadt_2" "Holzgünz_2"  
  ips.year.avg %>% 
    filter(falsto_name == "Holzgünz_2") 

  dat %>% 
    filter(falsto_name == "Brunnthal_1") %>% 
    filter(year == 2016:2017) 
  
# Filter beetle data 
  # 12 locations of zero beetles for IPS: 
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



## Inspect data ------------------------------------------------------------

# check for zeros presence - present! 
# how many zeros I have for sum beetles per trap?? - if 0 beetles over whole year, it is suspicious..

# filter data:
#  - set dates: April 1st (DOY 91, DOY 92 lap year) to Oct 31  (304) - from Phenips
#  - check revisit times: exclude if traps not collecte reularly! (eg. 10 times over the season: April 1st to Oct 30)

# Filter traps:
# 1) by 0 sum beetles per year
# 2) by low revisit frequency - once per year

# zero count for whole year, with single recording times... seems fishy
# seems that traps with zeros are consistent...
  # 08/05/2023 -> exclude zeros, but adrees this later
zero_catch_id <- ips.year.avg %>% 
  filter(sum_beetles == 0 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(falsto_name) %>% 
  pull()

  # 4 locations: [1] "Brunnthal_1"   "Weißenstadt_1" "Weißenstadt_2" "Holzgünz_2" 
  
 

# low revisit time
low_visit_id <- 
  ips.year.avg %>% 
  filter(freq_visit < 5 & art == 'Buchdrucker') %>% 
  ungroup(.) %>% 
  dplyr::distinct(falsto_name) %>% 
  pull()

(low_visit_id)

# does the average number of beetles per trap differ between locations?
# between years?
#dev.off()

ips.year.avg %>%
  ggplot(aes(x = factor(year),
             y = avg_beetles_trap,
             fill = art)) +
  geom_boxplot() +
  facet_wrap(. ~ art , scales = 'free_y') +
  theme(legend.position = 'bottom')




## Clean data: IPS ---------------------------------------------------------
# - keep only IPS
# - vegetation period: doy
# - years: 2015-2021
# - frequent recording: 10 revisit times at least
# - remove the '3rd' trap
# - keep only records that have all years and all traps

# check how many traps have records in March ---------------------------------
# get first a clean dataset for the whole year
# check how many records i have in March 

# get clean data for each month in a year
  dat.ips.clean.year <- dat %>% 
   # filter(doy %in% veg.period) %>% 
    filter(year > 2014) %>% 
    filter(art == 'Buchdrucker') %>%
    ungroup(.) %>% 
    dplyr::filter(!falsto_name %in% low_visit_id) %>%  # remove traps with low visit time (nee to check this for year?)
    dplyr::filter(!falsto_name %in% zero_catch_id) #%>% # exclude if zero beetles caught per whole year
  
#dat.counts.month <- 
dat.ips.clean.year %>% 
  #filter(doy %in% veg.period) %>% 
   ungroup(.) %>%
  group_by(month, year) %>% 
  summarize(n = n()) %>%
  pivot_wider(names_from = year, values_from = n, values_fill = list(n = 0))


dat.ips.clean.year %>% 
  #filter(doy %in% veg.period) %>% 
  ungroup(.) %>%
  group_by(year, representativ) %>% 
  summarize(n = n()) %>%
  ungroup(.) %>% 
  group_by(year) %>% 
  mutate(sum = sum(n),
         prop = n/sum*100)
 # pivot_wider(names_from = year, values_from = n, values_fill = list(n = 0))

#table(month, year)
#View()
#dplyr::filter(!falsto_name %in% low_visit_id) %>%  # remove traps with low visit time (nee to check this for year?)
#dplyr::filter(!falsto_name %in% zero_catch_id) #%>% # exclude if zero beetles caught per whole year



dat.ips.clean <- dat %>% 
  filter(doy %in% veg.period) %>% 
  filter(year > 2014) %>% 
  filter(art == 'Buchdrucker') %>%
  ungroup(.) %>% 
  dplyr::filter(!falsto_name %in% low_visit_id) %>%  # remove traps with low visit time (nee to check this for year?)
  dplyr::filter(!falsto_name %in% zero_catch_id) #%>% # exclude if zero beetles caught per whole year
 
# filter only for traps that have both 1&2:
# https://stackoverflow.com/questions/61060750/r-better-way-to-filter-by-group-and-vector-of-values
dat.ips.clean <- dat.ips.clean %>% 
  group_by(monsto_name) %>%      # filter by region
  slice(which(all(c(1,2) %in% trap_n) & trap_n %in% c(1,2)))   # every trap needs to have both 1&2 traps 


# filter only traps that have records over all years
dat.ips.clean <- dat.ips.clean %>% 
  group_by(falsto_name) %>%
  dplyr::filter(n_distinct(year) == n_distinct(dat.ips.clean$year))

dat.ips.clean %>% 
  group_by(year, trap_n) %>% 
  distinct(falsto_name) %>% 
  count() #%>% 

# if keeping only traps that are continously monitored over time 
# ->  79*2 traps/year
consistent_trapIDs <- unique(dat.ips.clean$falsto_name)

# filter the data counts on doy level to have also just consisten trap recordings ---------------------

dat_DOY_count_consist <- dat_DOY_count %>% 
  dplyr::filter(trapID %in% consistent_trapIDs)


data.table::fwrite(dat_DOY_count_consist, 
                   paste(output_path, 'outTable/dat_DOY_counts_IPS.csv', sep = "/"))

# Get daily counts per trap (averaged by the revisit time) ---------------------------------------------
df.daily <-  
  dat.ips.clean %>% 
    group_by(falsto_name, year) %>% 
    arrange(kont_dat) %>% # order the data
    mutate(
      period = diff(c(doy.start, doy))+1,  # get difference between two consecutive numbers
      avg_day_count  = fangmenge/period,           # average beetle counts per day
      sum = avg_day_count*period,                  # to check counts
      cumsum = cumsum(sum),
      diff   = avg_day_count - lag(avg_day_count) ) %>%   # difference between average avg_day_count values
  dplyr::select(-c("objectid",
                   "aelf",  
                   "art",
                    "einheit",
                   "koederwechsel",
                   "monsto_name",
                   "representativ")) 


length(unique(df.daily$falsto_name))

# find DOY when the cumsum overpasses values: XX beetles per season
# Warning levels from fovgis: 
# per week: 1000 beetles - warning: expected spread
#           3000 beetles - danger: high risk of standing infestation
#           5000 beetles
beetle_threshold = 1000

ips.aggreg <- df.daily %>% 
  group_by(year, falsto_name) %>% 
  arrange(doy) %>% 
  dplyr::filter(cumsum > beetle_threshold) %>%
  dplyr::filter(row_number()==1)     # filter first record per group > 2000
  
ips.aggreg %>% 
  ggplot(aes(y=doy,
             x = year,
             group = year)) +
  geom_boxplot()


# Sensitivity analysis  ---------------------------------------------------


# Loop over threshold values
beetle_thresholds <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3:10)*1000

# Assuming df.daily is your dataframe and it contains columns: year, falsto_name, doy, cumsum
# Replace 'cumsum' with the actual column name in your dataframe that should exceed the threshold

results <- lapply(beetle_thresholds, function(threshold) {
  df.daily %>%
    group_by(year, falsto_name) %>%
    arrange(doy) %>%
    filter(cumsum > threshold) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(beetle_threshold = threshold)  # Add threshold value to the results for identification
})

# Combine all results into a single dataframe
final_results <- bind_rows(results)

tab_trap_counts <-as.data.frame(table(final_results$beetle_threshold))
tab_trap_counts$Var1 <- as.numeric(as.character(tab_trap_counts$Var1))

# Define a custom function to calculate mean and sd
mean_sd <- function(x) {
  data.frame(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
}



p.supp.agg <- 
  final_results %>% 
  ggplot(aes(x = beetle_threshold,  # Ensure beetle_threshold is treated as a categorical variable
             y = doy,
             color = ifelse(beetle_threshold == 1000, 'thresh 1000', 'Other'))) + 
  stat_summary(fun.data = mean_sd) +
  stat_summary(fun = mean, geom = "point", size = 3) +
    scale_color_manual(values = c('thresh 1000' = 'red', 'Other' = 'black')) +
  #geom_vline(xintercept = 1000, col = 'red', lty = 'dashed') +
  theme_classic() + 
  xlab('Beetle colonization threshold') + 
  ylab('DOY') + 
  theme(aspect.ratio = 1,
        text = element_text(size = 14), # Set general text size
        axis.title = element_text(size = 14), # Set axis titles text size
        axis.text = element_text(size = 14), # Set axis text (ticks) size
        plot.title = element_text(size = 14)) +
    guides(color = FALSE) # This removes the legend for color

p.supp.agg


p.supp.trap_counts <- 
  tab_trap_counts %>% 
  ggplot(aes(x = Var1,
             y = Freq,
             color = ifelse(Var1 == 1000, 'thresh 1000', 'Other'))) + 
  stat_summary() +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(values = c('thresh 1000' = 'red', 'Other' = 'black')) +
  theme_classic() + 
  #geom_vline(xintercept = 1000, col = 'red', lty = 'dashed') +
  xlab('Beetle colonization threshold') + 
  ylab('Trap count [#]') +
  theme(aspect.ratio = 1,
        text = element_text(size = 14), # Set general text size
        axis.title = element_text(size = 14), # Set axis titles text size
        axis.text = element_text(size = 14), # Set axis text (ticks) size
        plot.title = element_text(size = 14)) + 
    guides(color = FALSE)

ggarrange(p.supp.agg, p.supp.trap_counts, ncol = 2, labels = 'auto')



# Plots -------------------------------------------------------------------

 
## Daily increases: ---------------------------------------
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




#windows()
hist(df.daily$avg_day_count)



# Export XY and DAT table -------------------------------------------------

# only for the selected traps - 78*2 per year
# and only one globalid per trap
trap_names    <- unique(dat.ips.clean$falsto_name)
trap_names_sf <- unique(xy_sf$falsto_name)

intersect(trap_names,trap_names_sf)

# Make a use of gpkg in R: https://inbo.github.io/tutorials/tutorials/spatial_standards_vector/
# Write the output XY files in 3035 coordinate

# filter XY for unique falso_name  - one falsto_name per globalid, keep just one row
xy_sf_fin <- xy_sf %>% 
  filter(falsto_name %in% trap_names) %>% 
  group_by(falsto_name) %>% 
  slice(1)   # keep only 1 globalid per group

unique(xy_sf_fin$falsto_name)

# save as a new object
xy_sf_all <- xy_sf

#rm(xy_sf)  # remove the old one

xy_sf_all %>%  
 st_write(paste(output_path, "outSpatial/xy_all_3035.gpkg", sep = '/'), append=FALSE)

xy_sf_fin %>%  
  st_write(paste(output_path, "outSpatial/xy_fin_3035.gpkg", sep = '/'), append=FALSE)

# export RAW data
# data.table::fwrite(dat, paste(output_path, 'outSpatial/dat.csv', sep = "/"))



# get study period
all_years    = 2006:2021
study_period = 2015:2021

# Expand coordinates table to have a geometry for each trap and year
xy_sf_expand <- 
  xy_sf_all %>% 
  filter(falsto_name %in% trap_names ) %>% 
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2]) %>% 
  dplyr::select(von_dat, bis_dat, falsto_name, globalid, x, y) %>% #, geometry
  mutate(from_year = lubridate::year(dmy(von_dat)))  %>% # get the year from the starting date
  group_by(falsto_name) %>% 
  complete(from_year = all_years)  %>% 
  arrange(falsto_name, from_year) %>% 
  # df %>%
  tidyr::fill(globalid, .direction = 'down') %>%
  tidyr::fill(x, .direction = 'down') %>%
  tidyr::fill(y, .direction = 'down')%>%
  # dplyr::rename(from_year = year) #%>% 
  filter(from_year %in% study_period) %>% # filter only years 2015:2021
  st_as_sf(coords = c("x", "y"), 
           crs = 3035) %>% 
  mutate(year = from_year) %>% 
  dplyr::select(-c(from_year, von_dat, bis_dat))


xy_sf_expand$id <- 1:nrow(xy_sf_expand)

# export XZ file with shifting traps
xy_sf_expand %>%  
  st_write(paste(myPath, "outSpatial/xy_fin_years_3035.gpkg", sep = '/'), append=FALSE)


# Alternatively, save only specific objects to the .Rdata file
#save(list=c("xy_sf_expand"), file="outData/spatial.Rdata")







# Prepare for map -------------------------------------------------------


# get the summary IPS data to merge with XY coordinates!  
ips.year.sum <-   dat.ips.clean %>%     
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_ips = sum(fangmenge, na.rm =T)) %>%  
  dplyr::select(sum_ips, year, falsto_name) #%>% 


# Merge spatial data with sum data by trap name (falsto_name)
ips.year.sum_sf <- xy_sf_fin %>% 
  right_join(ips.year.sum,  'falsto_name')

unique(ips.year.sum_sf$falsto_name)


# Get spatial data libs -----------------------------------------------------------------
library(rnaturalearth) # for map data
library(ggspatial)


# get spatial data ------------------------------------------------------------------ 
de_sf <- ne_states(country = "germany", returnclass = "sf")

# Get only bavaria
bav_sf <- de_sf %>% 
  dplyr::filter(name_en == "Bavaria")

bav_sf %>%  
  st_write(paste(output_path, "outSpatial/bavaria.gpkg", sep = '/'), append=FALSE)


buch_df_f <- ips.year.sum_sf %>% 
  filter(sum_ips > 3000)

# for visibility, only beetles sums > 3000/year/trap are shown!
ggplot(bav_sf) +
  geom_sf(color = 'black', 
          fill  = 'grey93') + 
  geom_sf(data = buch_df_f,
          aes(color = sum_ips,
              size = sum_ips)) + # , size = Age, size = 0.8size by factor itself!
  scale_color_viridis(name = 'Beetle sum/year/trap', 
                      alpha = 0.8,
                      option = 'magma',
                      direction = -1,
                      na.value = 'transparent') +
  facet_wrap(.~year) + 
  theme_void()

# Max increase in counts per year & trap ----------------------------------------------------
# find min DOY of the max diff per year and location: plot on map 
max.diff.doy <- df.daily %>% 
  group_by(year, falsto_name) %>% 
  slice(which.max(diff))

max.diff.doy.sf <- xy_sf_fin %>%   # to keep the sf structure, need to add df to sf
  right_join(max.diff.doy, 'falsto_name')


# Get spatial data: cumulative DOY -----------------------------
ips.aggreg.sf <- xy_sf_fin %>%   # to keep the sf structure, need to add df to sf
  dplyr::select(falsto_name, geometry) %>% 
  right_join(ips.aggreg, 'falsto_name')


head(ips.aggreg)


### Table: Max increase check range DOY: -------------------------------------
max.diff.doy %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  dplyr::summarise(min    = min(doy),
            max    = max(doy),
            mean   = mean(doy),
            sd     = sd(doy),
            cv     = sd/mean,
            median = median(doy))


# Check Max increase population
max.diff.doy %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  dplyr::summarise(min = min(diff),
            max = max(diff),
            mean = mean(diff),
            sd = sd(diff),
            cv = sd/mean,
            median = median(diff)
            )


# the main diference in counts happends earlier in the season
avg.doy <- mean(df.daily$doy)




ips.aggreg %>% 
  as.data.frame() %>% 
  group_by(year) %>% 
  dplyr::summarise(min    = min(doy),
                   max    = max(doy),
                   mean   = mean(doy),
                   sd     = sd(doy),
                   cv     = sd/mean,
                   median = median(doy))


# Maps facet  -------------------------------------------------------------------
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
#windows()
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
#windows()
p_diff_doy <- max.diff.doy %>% 
  #filter(doy < 151) %>% 
  ggplot(aes(x = doy,
             y = diff,
             color = factor(year)))+
  geom_point() +
  facet_wrap(year~.) +
  ylab('Max increase of beetle counts\n[diff]') +
  theme_bw()
  


#DOY
# 91-92 April 1st (92 in lap year, etc)
# 100 April 10
# 150 May 30
# 200 July 19
# 211 July 30
# 250 SEpt 7

# 1. beetle sume per year
ips.year.sum


p_ips.year.sum <- ggplot(data = ips.year.sum,
                             aes(x = year, 
                                 y = sum_ips  )) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = ips.year.sum,
    mapping = aes(x = year, 
                  y = sum_ips  ),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ggtitle("Beetle sum/year") + 
  theme_bw()

p_ips.year.sum

# 2.population increase -------------------------------------------------------------


# three aspects: the max increase in season
#                the size of population
#                beetle aggregations: DOY of getting beetle_threshold beetles per trap
# driver: previous year population??




## 3. Aggeragte beetle numbers (beetle_threshold per trap) IPS -------------------------------
lab_doy_cumulative = paste('DOY of cumulative \n', beetle_threshold)

p_aggreg <- 
  ggplot(bav_sf) +   # base map
  geom_sf(color = 'black', 
          fill  = 'grey80') + 
  geom_sf(data = filter(ips.aggreg.sf, doy < 130),
          aes(color = doy,
              size  = doy
          )
  ) + 
  scale_color_viridis(name = lab_doy_cumulative, #paste('DOY of cumulative ', beetle_threshold),
                      alpha = 0.8,
                      option = 'magma',
                      direction = 1,
                      na.value = 'transparent') +
  scale_size(trans = 'reverse') +
  theme_void() +
  facet_wrap(~year) +
  xlab("Longitude") + 
  ylab("Latitude") +
  labs(title = paste('Beetle aggregation ', beetle_threshold),
       subtitle = 'doy <130',
       color  = "DOY",
       size = "DOY")

(p_aggreg)


ips.aggreg


# barplot with medians + IQR
p_ips.agg <- ggplot(data = ips.aggreg,
       aes(x = year, 
           y = doy)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = ips.aggreg,
    mapping = aes(x = year, 
                  y = doy),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ggtitle(paste("Aggregation: DOY ", beetle_threshold,  "beetles")) + 
  theme_bw()
  #coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)

# Coefficient of variation:  -------------------------------------------------------

# for whole data, not just the max increase
df.daily %>%
  group_by(year, falsto_name) %>% 
  dplyr::summarize(mean = mean(avg_day_count),
            sd  = sd(avg_day_count),
            cv  = sd/mean) %>%  # can be turned to % by '*100'
 # as.data.frame() %>% 
  ggplot(aes(x = factor(year),
             y = cv)) +
  geom_violin() +
  stat_summary() +
  ggtitle('CV of daily beetle counts')


# CV of DOY of the max beetle increase ------------------
max.diff.doy %>%
#  ungroup(.) %>% 
  group_by(year) %>% 
  dplyr::summarize(mean = mean(doy),
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





# DOY max increase: barplot + medians

p_ips.max.diff.doy <- ggplot(data = max.diff.doy,
                    aes(x = year, 
                        y = doy)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = max.diff.doy,
    mapping = aes(x = year, 
                  y = doy),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ggtitle("Emergence: DOY of max increase") + 
  theme_bw()



# DOY max increase: amount of beetles
p_ips.max.diff <- ggplot(data = max.diff.doy,
                             aes(x = year, 
                                 y = diff)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = max.diff.doy,
    mapping = aes(x = year, 
                  y = diff),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ggtitle("Emergence: Max increase (diff in #beetles)") + 
  theme_bw()

# get correlation between trap1 and trap 2
# get function to plot the lm equation and R squared:
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)


# Linear regression between trap 1&2 --------------------------------------------
df_compare_traps <- ips.year.sum %>% 
  mutate(trap_n =  substr(falsto_name, nchar(falsto_name)-1+1, nchar(falsto_name)),
         loc =  substr(falsto_name, 1, nchar(falsto_name)-2))%>%
  dplyr::select(-falsto_name) %>% 
  tidyr::spread(trap_n, sum_ips) %>% 
  dplyr::rename(x = '1',         
                y = '2') #%>% 
  

# get lm between trap 1 & 2
p_lm_traps <- 
  ggplot(df_compare_traps, aes(x = x/10^3,
                               y = y/10^3,
                              # color = year,
                               group = year)) +
  stat_poly_line() +
    stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
    geom_point(alpha = 0.5) +
    facet_wrap(.~year) + #, scales = 'free'
    theme_bw()
  
  



# Save  ------------------------------------------------------------
save(trap_names,                 # vector of final trap names (79*2)
     df.daily,                   # avg daily beetle counts (adjusted by revisit times), over DOY 
     ips.year.avg,               # sum & avg number per beetles/trap per vegetation season (), over whole year
     ips.year.sum,               # simple beetle sum per trap
     ips.year.sum_sf,            # simple beetl sum per trap sf
     dat.ips.clean,              # cleaned IPS counts (only consisten trap selected, checked for revisit time,...)
     max.diff.doy,               # max increase in beetle counts per DOY
     max.diff.doy.sf,            # sf: max increase in beetle counts per DOY
    ips.aggreg,                  # beetle agregation per trap: DOY of overpassing XX beetles 
     file="outData/ips_counts.Rdata")


# Spatial data for maps:
save(xy_sf_all,                  # all locations (one trap has several locations!!)
     xy_sf_fin,                  # filtered trap location (one Globalid is ised for both IPS & Chalcographus!)
     xy_sf_expand,               # all XY data for shifting traps over years
    # de_sf,                      # germany shp
    # bav_sf,                     # bavaria shp
     file="outData/spatial.Rdata") 
