

##########################################################
#
# 1. Do beetle counts data correlate with RS data??         #
#
##########################################################

# Process:

# get input data
# get a XY data
# make a buffers around XY data
# get info of mortality over years from RS per cell
# aggregate traps into cells, standardize by number of traps per year/cell
# merge databases by cell_id
# does it correlate over Bavaria? over years?


# Input data:

# bavaria shp
# disturbance raster
# xy traps
# beetle data
# raster species composition: deciduous/coniferous


# Read my paths -----------------------------------------------------------
source('myPaths.R')
source('process_input_data.R')


# Process:

# Make buffers
# Intersect buffers with forest and disturbance maps
# export as values from rasters to have a dataframe
# for each trap: 
# 

# Convert objects into terra formats to speed up extraction
buff     <- vect(xy_500)
for_type <- rast(forest_type)
dist     <- rast(disturbance14)

# Convert buffers to raster
buff_ras <- fasterize(xy_500, disturbance14, field = "OBJECTID") # name the new rasters as polygon ids
buff_ras <- rast(buff_ras)   # convert to terra format
#grid_values  <- values(grid_sel_ras)     # resolution is still 30 m, grid value for each pixel; this is a vector of values


# Export raster to check it on the QGIS
# writeRaster(buff_ras, paste(myPath, outFolder, 'buff_ras.tif', sep = '/'), overwrite=TRUE)


# Resample forest raster to match/snap/align the disturbance raster
# Buffer raster is already resampled
for_type_resample <- terra::resample(for_type,  # raster to be resampled 
                                     dist,      # Master raster
                                     method = 'near')


# Extract by mask; creates a list of datasets
forest_ex <- terra::extract(for_type_resample, buff, list = F)
dist_ex   <- terra::extract(dist,     buff, list = F)
buff_ex   <- terra::extract(buff_ras, buff, list = F)



# Does the extend fit?? YES
nrow(forest_ex)
nrow(dist_ex)
nrow(buff_ex)


# Get sums of pixels by buffers
forest_df <- 
  forest_ex %>% 
  #filter(bav_fortype_ext30_int2u_LZW != 0) %>% 
  group_by(ID) %>%
  summarize(decid_ha = sum(bav_fortype_ext30_int2u_LZW == 1, na.rm = TRUE) * 0.09,
            conif_ha = sum(bav_fortype_ext30_int2u_LZW == 2, na.rm = TRUE) * 0.09 #,
            #land_ha = n() * 0.09  # not sure if calculating land is meaningful here?
            ) %>%
  ungroup(.)


# Get in the same way the disturbance data aggregated by grid
dist_df <- 
  dist_ex %>%
  na.omit(.) %>%
  group_by(ID, disturbance_map_bavaria) %>%
  summarize(disturbance_ha = n() * 0.09) %>%
  rename(year = disturbance_map_bavaria ) %>% 
  ungroup(.)


# Compare two vectors to fond out which ones are missing? ----------------------------------
#a <- 1:497


setdiff(a, unique(buff_ex$layer))
#  [1]  61  62 162 166 167 257 291 296 297 298 300 301 306 308 314 324 325 327 331 332 337 338 350 363    # in total, 24 traps (5% of all traps)
# for now I can just proceed with the buffer_id/ or double the infos with the other databases?

# Buffer:
# - layer is a raster number
# - id - the number in a list

# Chck it for the buffer: 
buff_ex %>% 
  #filter(layer == 61)  # no rows to return
  filter(layer == 370)  

# Need to merge the data by the ID !
# my ID is already in teh forest_df and dist_df, maybe I do not need to have also buffer indication?
# just directly link with the XY value 
  


# Seems that some data have the same coordinates?
# Eg. 61 and 370?
# Check teh coordinates?
xy %>% filter(OBJECTID == 61 | OBJECTID == 370) #%>% # They have exactly teh same geometry, but different globalid
# Subset the data by the 
duplic_globid_v <- xy %>% 
  filter(OBJECTID == 61 | OBJECTID == 370) %>% 
  distinct(globalid) %>% 
  pull()
  

# Look at the beetle database, why the data can be duplicated
dat %>% 
 # filter(globalid == "A94A9D61-E305-495D-83E1-9DBCBD395977") %>% 
  filter(globalid %in% duplic_globid_v) %>% 
  filter(representativ != 'Ja') %>% 
  print(n = 40)



# Need to aggregate the buff values as well, as now they are for each pixel
buff_df <- 
  buff_ex %>% 
    distinct()    # one ID can have several buffer numbers
  
  
 # group_by(layer, ID) %>% 
#  summarize(buff_id = mean(layer))
  
length(unique(buff_ex$ID))    # General ID, order,              497
length(unique(buff_ex$layer)) # indicator of the XY OBJECTID    473

# Merge forest and disturbance data by years, by the buff_id or ID - should be the same ------------------
dist_for_df <- 
  forest_df %>% 
  right_join(dist_df, by = "ID") %>%
  #right_join(buff_ex, by = "ID") %>% 
  rename(OBJECTID = ID)
    #print(n=40)
  
# CHeck the location oif idd 446 and 1: 25 m, close to each other; therefore the values are recorded two times, for value of record
#rename(OBJECTID = gridindex)

unique(buff_ex$layer)
sort(unique(dist_for_df$buff_id))
unique(xy$OBJECTID)


# The ID and buff_id are thow different things
#dist_for_df$ID == dist_for_df$buff_id

length(unique(dist_for_df$ID))      # 497
length(unique(dist_for_df$buff_id)) # 473



# get rid of geometry
xy_df <- xy %>%  
  st_drop_geometry() #%>% 
#rename(ID = OBJECTID)


# Link disturbance data to geometry data:
out_df <- dist_for_df %>% 
  left_join(xy_df, by = "OBJECTID")


# Link to beetle data, by globalid
head(dat)




# Standardize the beetle counts -------------------------------------------




# Join between XY data, RS mortality (by buffer ID) and beetle counts (globalid) ---------------
out_df <- 
  dist_df %>% 
  left_join(xy_df, by = "OBJECTID") %>%
  rename(year = dist_year )
#  tail()




# Get number of traps by buffer - always one here, unless they were very close to each other

# How many traps do I have per hexagon?
# they vary by year;
# n = number of traps by hexagon, year and species 
trap_n <- 
  dat_df %>% 
  filter(year > 2014 & representativ == 'Ja' ) %>% # & art == 'Buchdrucker' 
  group_by( #monsto_name,  # trap location name
    OBJECTID,     # hexagon ID
    year,
    art) %>% 
  tally() %>%
  rename(trap_n = n)
#print(n = 40)


# # Let's see if the counts are right: 
# # filter just specific values:  Aldersbach          224  2017    23
# dat_df %>% 
#   filter(monsto_name == 'Aldersbach' & year == 2017 & art == 'Buchdrucker' & OBJECTID == 218 & representativ == 'Ja')
# 

# dat_df %>% 
#   distinct(einheit)  # Stuck = number/count, NA?


# get number of beetles per hexagon by years
beetle_n <- dat_df %>% 
  group_by(art, year, OBJECTID) %>% 
  summarize(beetle_sum = sum(fangmenge, na.rm = T))

# Add the beetle count data to the number of traps & standardize the beetle counts by a trap number
out_df <- beetle_n %>% 
  left_join(trap_n, by = c('OBJECTID', 'art', 'year')) %>% 
  mutate(beetle_by_trap = beetle_sum/trap_n) %>% 
  left_join(dist_df) #%>% 
#mutate(lag_dist_sum = lag(dist_sum))


# Need to check out this part!!! how to merge teh data! keep the id as buff_id!
out_df2 <- 
  out_df %>% 
  right_join(dat, by = c('globalid', 'year', 'falsto_name', 'aelf'))  

# the beetle data have 490 globalids
# the XY has 497 ids

# which ones are missing?
setdiff(unique(xy$globalid), unique(dat$globalid))

# [1] "CE935FFE-0B77-41E8-8B8A-A48338954BD3" "BF168CB3-750D-4135-866A-15C8477FEB28" "CAE0A588-C73D-4389-AF43-EA22E31B0117" "16622072-13C8-4CCD-9A11-501A19B7A0A8" "EC1B6963-339B-40AC-B87F-FC8DF8CD3566"
# [6] "6CECBC3C-F69B-4A5C-9E41-B053F92DA5B4" "37E3B5F6-76E4-4CDA-91B4-FADA03A033DD"

# Are those really missing from beetle data?
dat %>% 
  filter(globalid == "CE935FFE-0B77-41E8-8B8A-A48338954BD3")  # thisis explained the notes, that several sites have been moved, but have the same name, althought unique globalid


# FINALIZE THE DATABASE!!

# Questions: do counts corresponds to total yearly mortality?
# -----------------------------------------------------------------------------

# get beetle sums per year: IT, PC individually, and together
# lag the RS data by one year
# get beetle sums until/before july - need to calculate the daily counts from the dates?

# Get yearly sums of beetle counts:
out_df2 %>% 
  group_by(OBJECTID,  # this is a buffer indicator
           year,
           art
  )



# Add lagged disturbance data by one year:
# need to be correctly arranged;
out_df <- out_df %>% 
  filter(year > 2014) %>% 
  group_by(art, OBJECTID) %>% # , OBJECTID 
  arrange(year, .by_group = TRUE) %>% 
  mutate(lag_dist_sum = lag(disturbance_ha))




# Check the correlation between the beetle numbers and damage from RS:
# just using stat_smooth
out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = disturbance_ha,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_grid(art ~ year, scales = 'free')


out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Kupferstecher') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = disturbance_ha,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_grid(art ~ year, scales = 'free')



# Check out the effects of the lagged values
windows()
out_df %>% 
  filter((art == 'Kupferstecher' & beetle_by_trap < 30000) | 
           (art == 'Buchdrucker' & beetle_by_trap < 6000 )) %>%  #filter(year > 2014 & beetle_by_trap< 30000) %>% 
  #filter(art == 'Buchdrucker'   & beetle_by_trap < 6000) #%>% 
  # filter(art == 'Buchdrucker') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = lag_dist_sum,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") +
  facet_wrap(art ~ ., scales = 'free')


windows()
out_df %>% 
  filter(year > 2014) %>% 
  filter(art == 'Kupferstecher') %>% 
  ggplot(aes(x = beetle_by_trap,
             y = lag_dist_sum,
             color = art)) +
  geom_point() + 
  geom_smooth(method= "loess") #+
#facet_grid(art ~ year, scales = 'free')






# Check the overall counts and damages over years? also, the high rates of beetles can predict mortality in 
# following years?
year_sum_df <- 
  out_df %>%
  filter(year > 2014) %>% 
  group_by(art, year) %>% 
  summarize(beetle_sum = sum(beetle_by_trap , na.rm = T),
            dist_sum   = sum(disturbance_ha, na.rm = T)) %>% 
  mutate(lag_dist_sum = lag(dist_sum))


# Plot the yearly values nad lagged values
windows()

# Plots yearly values
p1<- year_sum_df %>% 
  ggplot(aes(x = beetle_sum,
             y = dist_sum,
             color = art)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(.~art, scales = 'free')


# Plot lagged values
p2 <- year_sum_df %>% 
  ggplot(aes(x = beetle_sum,
             y = lag_dist_sum,
             color = art)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(.~art, scales = 'free')


windows()
ggarrange(p1, p2, ncol = 1)
