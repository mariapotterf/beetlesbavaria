

# -----------------------------------
# Explore bark beetle data
# -----------------------------------

# Libs -----------------------------------------------------------
library(lubridate)
library(dplyr)
library(ggplot2)
library(rgdal)   # for spatial data


# get data
path = 'C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria/rawData/Fwd__Borkenkäferforschung__Datentransfer'

# Load RData table
load(paste(path, 'BoMo_2015_2021_Rohdaten.RData', sep = "/"))

# Read XY coordinates
xy <- read.csv2(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit.csv", sep = '/'))

head(xy)

xy <- read.delim(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit2.txt", sep = '/'), 
               sep = ';', dec = ".")

head(xy)


xy <- read.csv(paste(path, "BoMo_Fallenstandorte_Lage_Laufzeit2.csv", sep = '/'), 
                 sep = ',', dec = ",",skip = 0, encoding = "UTF-16LE")

head(xy)

str(xy)




# Test decimals reading
xy <- read.delim(paste(path, "test_decimals2.txt", sep = '/'), 
               sep = '\t', dec = ",",skip = 0)
head(xy)


# Investigate beetle data: Bavaria

# Change name of the table
dat <- Daten_B01


# Explore the table
head(dat)

# what beetles species?
unique(dat$art)

# "Buchdrucker"   "Kupferstecher"

# what year?
unique(dat$kont_dat)

str(dat)  # datum is in POSIXct format, ut does not contain time: can convert just to date

# convert to date
dat$kont_dat <- as.Date(dat$kont_dat)

str(dat)

# Species
# Buchdrucker  : Ips typographus 
# Kupferstecher: Pityogenes chalcographus      
# Variation within a year?
# yes, need to split in year, months, dates
dat <- dat %>% 
  dplyr::mutate(year = lubridate::year(kont_dat), 
                month = lubridate::month(kont_dat), 
                day = lubridate::day(kont_dat))


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

# 



