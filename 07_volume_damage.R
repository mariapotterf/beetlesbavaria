# Proces volume data

# from LWF: the damage volume+ geometries
# read tables
# merge geometries
# explore damage per year
#  read ata from Cornelius: do they correlate somehow?
# how to link them with beetle trap data?


# Libs
library(dplyr)
library(terra)
library(data.table)
library(ggplot2)

# read districts shp
districts <-  vect('rawData/damage_volume/AELF_Revier_20240502/AELF_Revier_bis20240424_LWF_20240502.shp')
damage    <- fread('rawData/damage_volume/Bavaria_TreeDamageData_IpsTypographus_2015-2021.csv', dec = ',')

shp_ID   <- unique(districts$forstrev_1)  # 336 districts numbers
district_ID <- unique(damage$AELF_district_ID)  # 337 district numbers

length(shp_ID)
length(district_ID)

setdiff( c(1,2), c(1) )

setdiff(district_ID, shp_ID)
summary(districts)
summary(damage)

# there is an aditiona row in the damage data with district ID 0 - I can remove that

hist(damage$damaged_volume_total_m3)


damage %>%
  ggplot(aes(x = Year, 
             y = damaged_volume_total_m3, 
             group = Year,
             fill = factor(Year))) + 
  geom_violin() #+
  #geom_boxplot()

