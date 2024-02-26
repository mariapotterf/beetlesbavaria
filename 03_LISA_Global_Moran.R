

# Spatial patters analysis: 

# LISA - local indicator of spatial association
# Global Moran's I
# semivariance = analyze the spatial dependence or variability of a dataset;  quantifies the average dissimilarity between pairs of observations at different distances or lags. 
# variogram


# from ChatGPT:

# Split data by years, have XY trap locations
# 1. Local Variation Analysis: investigate  concept of spatial autocorrelation. 
#    Compute (LISA) using the localmoran() function from the spdep package.
#    - will find the hotspots (High-High, Low-Low)/coldspots, clusters  
#    - can visualizeLISA using map 

# 2. Regional Variation Analysis: overall spatial autocorrelation in  data. 
#  - Calculate the global Moran's I statistic using the spdep::moran.test() 
#  - This test will determine whether the beetle counts are spatially autocorrelated at a broader regional level.

# 3. Scale-dependent Analysis: 
#  - compute semivariograms or variograms to assess the spatial dependence of beetle counts at different lag distances. 
#  - The gstat::variogram() 
# Plot the semivariograms to observe how the spatial dependence changes across different scales.


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

library(sf)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)
#library(raster)
library(rgdal) # for sp data
#library(tidyverse)
#library(lubridate)
#library(patchwork)
library(fasterize)
library(terra)
library(ggplot2)
library(ggpubr)
library('cowplot') # arrange plots


# Spatstat
library(spdep)

# load cleaned data
load("outData/ips_counts.Rdata")
load("outData/spatial.Rdata")

# Get beetle counts - corrected, instead of previous 'dat'
# - dat.ips.clean      # dat <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))
# - max.diff.doy       # get DOY of max increase
# - ips.year.sum       # sum beetles/year per trap

# Vars
n_neighbors = 5      # number of nearest neighbors

# Spatial data: 
# should have 158 regions (trap pairs)
sort(unique(xy_sf_expand$falsto_name))


# get coordinates from sf object
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    year = xy_sf_expand$year,
                    falsto_name = xy_sf_expand$falsto_name)

xy_df <- distinct(xy_df)


# Get sums of IPS beetle per year/trap: March 31 to October 30
ips_sum <- 
  dat.ips.clean %>% 
  group_by(year,falsto_name) %>% 
  dplyr::summarize(sum_beetle = sum(fangmenge, na.rm = T)) %>% 
  left_join(xy_df, by = c("falsto_name", 'year')) # df_xy3035


nrow(ips_sum)

# Run LISA for each year separately
years <- 2015:2021

get_lisa <- function(i, ...) {
 
  #i = 2015
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  #
  #  # Local variation analysis (LISA) per year:
  #  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  #
  #  # Create queen contiguity neighbors object
  nb <- knn2nb(knearneigh(coordinates(ips_sum_sub),
                          k = n_neighbors),  # nearest neighbor
               sym = TRUE)
  
  #  # Calculate LISA
  lisa_res<-localmoran(ips_sum_sub$sum_beetle, nb2listw(nb,style="W"))
  #
  # add LISA to points:
  ips_sum_sub$Morans_I <-lisa_res[,1] # get the first column: Ii - local moran  stats
  ips_sum_sub$clust <-attributes(lisa_res)$quadr$mean  # get classified data
  #
  # convert to sf for plotting
  ips_sum_sub_sf <- st_as_sf(ips_sum_sub)
  
  return(ips_sum_sub_sf)
}

# Run LISA on all years separately
lisa_out <- lapply(years, get_lisa )

# Merge all in one sf
lisa_merged <- dplyr::bind_rows(lisa_out)




# START sensitivity to number of neighbors ------------------------------------------

get_lisa_neighbours <- function(year, n_neighbors, ips_sum) {
  # Slice specific year
  ips_sum_sub <- ips_sum %>% filter(year == year)
  
  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  
  # Create knn neighbors object
  nb <- knn2nb(knearneigh(coordinates(ips_sum_sub), k = n_neighbors), sym = TRUE)
  
  # Calculate LISA
  lisa_res <- localmoran(ips_sum_sub$sum_beetle, nb2listw(nb, style="W"))
  
  # Add LISA to points
  ips_sum_sub$Morans_I <- lisa_res[,1] # Local Moran's I statistic
  ips_sum_sub$clust <- attributes(lisa_res)$quadr$mean  # Quadrant mean for LISA clusters
  ips_sum_sub$n_neighbors <- n_neighbors
  
  # Convert to sf for plotting
  ips_sum_sub_sf <- st_as_sf(ips_sum_sub)
  
  return(ips_sum_sub_sf)
}

# Vector of neighbors to test
neighbor_counts <- c(1, 3, 5, 7, 9, 11, 15, 21)

# Years to loop over
years <- 2015:2021

# Perform sensitivity analysis over years and number of neighbors
sensitivity_results <- lapply(years, function(yr) {
  lapply(neighbor_counts, function(n_nb) {
    get_lisa_neighbours(yr, n_nb, ips_sum)
  })
})



# Merge all in one sf
sensitivity_merged <- dplyr::bind_rows(sensitivity_results)
sensitivity_merged_df <- sensitivity_merged %>% 
  as.data.frame()


p.supp.neighb <- 
sensitivity_merged_df %>% 
  ggplot(aes(x = n_neighbors                ,  # Ensure beetle_threshold is treated as a categorical variable
             y = Morans_I ,   #,
             color = ifelse(n_neighbors == 5, 'thresh 5', 'Other')
             )) + 
  stat_summary(fun.data = mean_sd) +
  stat_summary(fun = mean, geom = "point", size = 3) +

  scale_color_manual(values = c('thresh 5' = 'red', 'Other' = 'black')) +
  #geom_vline(xintercept = 1000, col = 'red', lty = 'dashed') +
  theme_classic() + 
  xlab('Nearest neighbours number threshold') + 
  ylab('Local Morans I' ) + 
  theme(aspect.ratio = 1,
        text = element_text(size = 14), # Set general text size
        axis.title = element_text(size = 14), # Set axis titles text size
        axis.text = element_text(size = 14), # Set axis text (ticks) size
        plot.title = element_text(size = 14)) +
  guides(color = FALSE) # This removes the legend for color


p.supp.neighb


 ##############END sensitivity analysis 


# can MoransI explain beetle counts??? - not a meaningful predictors for beetles counts
lisa_merged %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Morans_I ,
             y = sum_beetle)) +
  geom_point()

m1<-glm.nb(sum_beetle~Morans_I, lisa_merged)
plot(m1)
#st_crs(lisa_merged) <- crs(bav_sf)
#lisa_merged_proj <- st_crs(crs(bav_sf))
#  st_transform(lisa_merged, projection(bav_sf))
# Create map fr each one of them and save as saparate object


p_lisa_sub <- ggplot() +
  geom_sf(data = filter(lisa_merged, 
                        clust %in% c("Low-Low", "High-High")),
          aes(color = clust)) +
  scale_color_manual(breaks = c("Low-Low", "High-Low", "Low-High", "High-High"),
                     values=c("blue", "grey90", "grey90", 'red')) +
  facet_wrap(year~.) +
  theme_void() +
  ggtitle('LISA: Moran I')



p_lisa_all <- ggplot() +
  geom_sf(data = lisa_merged, #filter(lisa_merged, 
                        #clust %in% c("Low-Low", "High-High")),
          aes(color = clust)) +
  scale_color_manual(breaks = c("Low-Low", "High-Low", "Low-High", "High-High"),
                     values=c("blue", "grey90", "grey95", 'red')) +
  facet_wrap(year~.) +
  theme_void() +
  ggtitle('LISA: Moran I')


# improve plotting: cap values 

# Cap values greater than 5 to 5

hist(lisa_merged$Morans_I_capped)
range(hist(lisa_merged$Morans_I))
lisa_merged$Morans_I_capped <- pmin(lisa_merged$Morans_I, 1)
lisa_merged$Morans_I_capped <- pmax(lisa_merged$Morans_I_capped, -1)

# Define a color palette for the divergence scale
color_palette <- c("blue", "white", "red")

# Calculate the minimum and maximum Moran's I values
min_moran <- min(lisa_merged$Morans_I_capped)
max_moran <- max(lisa_merged$Morans_I_capped)


lisa_merged %>% 
  ggplot(aes(x = year,
             y = Morans_I)) +
  #geom_bar() +
  geom_bar(fun = "mean", stat = "summary")

# how can I see from this plot the my spatial synchonization is increasing ???
p_lisa_Morans <- ggplot() +
  geom_sf(data = lisa_merged,
          aes(color = Morans_I)) +
  # scale_color_gradientn(colors = color_palette,
  #                       values = scales::rescale(c(min_moran, 0, max_moran), c(0, 0.5, 1)),
  #                       breaks = c(min_moran, 0, max_moran),
  #                       labels = c("Low-Low", "No Spatial Autocorrelation", "High-High"),
  #                       limits = c(min_moran, max_moran)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, 
                        breaks = c(min_moran,0, max_moran),
                        labels = c("Low-Low", "No Spatial Autocorrelation", "High-High"),
                        limits = c(min_moran, max_moran)) +
  
    facet_wrap(year ~ .) +
  theme_void() +
  ggtitle('LISA: Moran I')

windows()
(p_lisa_Morans)

# convert to df 
lisa_merged_df <- as.data.frame(lisa_merged)


# get LISA distribution
lisa_merged %>% 
  ggplot(aes(Morans_I,
             group = year)) + 
  geom_density() 
  


lisa_merged %>% 
  ggplot(aes(y = sum_beetle ,
             x = Morans_I  ,
             color = year)) + 
  geom_point() + 
  geom_smooth()



p_lisa_freq <- 
  lisa_merged_df %>%
  group_by(clust, year) %>% 
    dplyr::summarize(freq = n()) %>%
    ggplot(aes(x = year,
               y = freq,
               fill = clust)) +
    geom_col()
    






  






# Global Moran =================================================================

# get more distant neighbors:
nb_g_15 <- knn2nb(knearneigh(coordinates(ips_sum_15),
                             k = n_neighbors),  # nearest neighbor
                  sym = TRUE)
nb_g_20 <- knn2nb(knearneigh(coordinates(ips_sum_20),
                             k = n_neighbors),  # nearest neighbor
                  sym = TRUE)


# Calculate the global Moran's I
global_moran_15 <- moran.test(ips_sum_15$sum_beetle, nb2listw(nb_g_15) )
global_moran_20 <- moran.test(ips_sum_20$sum_beetle, nb2listw(nb_g_20) )

# for all years:
get_global_moran <- function(i, ...) {
  
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  
  #  # Global Moran's per year: Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  #
  #  # Create queen contiguity neighbors object
  nb <- knn2nb(knearneigh(coordinates(ips_sum_sub),
                          k = n_neighbors),  # nearest neighbor
               sym = TRUE)
  
  #  # Calculate Global Moran's I
  global_moran_res <- moran.test(ips_sum_sub$sum_beetle, nb2listw(nb) )
  #
  # export as df
  df <- data.frame(year = i,
                   stat = global_moran_res$statistic,
                   p_val = global_moran_res$p.value )
  
  # export table with results
  return(df)
}

# Run Global mora on all years separately
global_moran_out <- lapply(years, get_global_moran )

# Merge all in one sf
glob_merged <- dplyr::bind_rows(global_moran_out)

# add significance
glob_merged <- glob_merged %>% 
  mutate(sign = case_when(p_val > 0.05 ~ 'n.s',
                          p_val < 0.05 & p_val > 0.01 ~ '*',
                          p_val < 0.01 & p_val > 0.001 ~ '**',
                          p_val < 0.001  ~ '***'
                          )) 

# Create a barplot of the values
p_glob_moran_bar <-ggplot(glob_merged, aes(x = year,
                        y = stat,
                        label = sign)) +
  geom_col(fill = 'lightgrey', col = 'black') +
  xlim(2014.3,2021.6) +
  ylab('Moran I statistic\nstandard deviate')+
  ylim(0,13) +
  geom_text(aes(x = year,
            y = stat + 1 )) + #rep(12.2, 7))) +
  theme_bw() 



# Get variograms: -----------------------------------------------------------------

library(gstat)

cutoff = 300000


# how to make a function, with output raw data and model?? for each year??
# get only raw date, the models need to be done visually: adjust sill, nugget and range
get_variogram <- function(i, ...) {
  
  #  # Slice specific year
  ips_sum_sub <-ips_sum %>%
    filter(year == i)
  #
  #  # Global Moran's per year:
  #  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
 
  # export as df
  variogram_raw <- variogram(log10(sum_beetle)  ~ 1, 
                             cutoff=cutoff,
                             data = ips_sum_sub)  # all directions are assumed equal
  # add year
  variogram_raw$year = i
  
   # export table with results
  return(variogram_raw)
}


# Run on all years separately
variogram_out <- lapply(years, get_variogram )

# Merge all in one df
var_merged <- dplyr::bind_rows(variogram_out)

# can't plot them in the one plot, as I can't get fitted models;
# make one by one instead


# plot the variograms and their fitted models one by one
cutoff = 300000

m15 <-  fit.variogram(variogram_out[[1]], 
                     vgm(psill = 0.2, model="Sph", range = 200000, nugget = 0.11))

#  # Slice specific year
ips_sum_16 <-ips_sum %>% filter(year == 2016)
coordinates(ips_sum_16) <- ~ x + y

variogram_out[[2]] <- variogram(log10(sum_beetle)  ~ 1, 
                           cutoff=cutoff,
                           data = ips_sum_16)  # all directions are assumed equal

m16 <-  fit.variogram(variogram_out[[2]], 
                      vgm(psill = 0.17, model="Sph", range = 40000, nugget = 0.12))

m17 <-  fit.variogram(variogram_out[[3]], 
                      vgm(psill = 0.2, model="Lin", range = 200000, nugget = 0.11))

m18 <-  fit.variogram(variogram_out[[4]], 
                      vgm(psill = 0.2, model="Sph", range = 200000, nugget = 0.11))

m19 <-  fit.variogram(variogram_out[[5]], 
                      vgm(psill = 0.2, model="Lin", range = 200000, nugget = 0.11))

#  # Slice specific year
ips_sum_20 <-ips_sum %>% filter(year == 2020)
coordinates(ips_sum_20) <- ~ x + y
variogram_out[[6]] <- variogram(log10(sum_beetle)  ~ 1, 
                                cutoff=400000,
                                data = ips_sum_20)  # all directions are assumed equal  # ,
#alpha = c(0, 45, + 90, 135)

variogram_cloud_20 <- variogram(log10(sum_beetle)  ~ 1, 
                                #cutoff=400000,
                                data = ips_sum_20,
                                cloud = T,
                                alpha = c(0, 45, + 90, 135))

# here cutoff needs to be 300000
m20 <-  fit.variogram(variogram_out[[6]], 
                      vgm(psill = 0.25, model="Gau", range = 400000, nugget = 0.11))

plot(variogram_out[[6]], m20, main = "2020", cutoff = cutoff)

m21 <-  fit.variogram(variogram_out[[7]], 
                      vgm(psill = 0.2, model="Sph", range = 200000, nugget = 0.11))



# Plot the fitted semivariogram
p15<-plot(variogram_out[[1]], m15, main = "2015", cutoff = cutoff)
p16<-plot(variogram_out[[2]], m16, main = "2016", cutoff = cutoff)
p17<-plot(variogram_out[[3]], m17, main = "2017", cutoff = cutoff)
p18<-plot(variogram_out[[4]], m18, main = "2018", cutoff = cutoff)
p19<-plot(variogram_out[[5]], m19, main = "2019", cutoff = cutoff)
p20<-plot(variogram_out[[6]], m20, main = "2020", cutoff = cutoff)
p21<-plot(variogram_out[[7]], m21, main = "2021", cutoff = cutoff)

windows()
semi_out <- plot_grid(p15, p16, p17, p18, p19, p20, p21, nrow = 3)# , label_size = 8


# year 2020 has many outliers:
hist(ips_sum_20$sum_beetle)

ips_sum_20 %>% 
  data.frame() %>% 
  dplyr::filter(sum_beetle > 40000)
  dplyr::filter(falsto_name %in% c('Auerbach_2', 'Auerbach_1',
                                   'Heinrichsthaler_Forst_2', 
                                   'Heinrichsthaler_Forst_1'))

# values over 100.000 beetles per year:
#     year             falsto_name sum_beetle       x       y optional
# 4   2020              Auerbach_2     121600 4548500 2858706     TRUE
# 52  2020 Heinrichsthaler_Forst_2     113920 4276581 2992330     TRUE
# 125 2020             Sugenheim_1     106408 4354050 2942596     TRUE
# 137 2020             Wegscheid_1      84930 4599349 2838491     TRUE


# Humblt: fitting variograms: https://gsp.humboldt.edu/olm/R/04_01_Variograms.html


# Example ----------------------
# eye-ball modell fitting


# 
# library(geoR)
# library(sf)
# 
# library(sp)        # for meuse dataset
# 
# library("gstat")   # geostatistics
# library("mapview") # map plot
# library("sf")      # spatial vector data
# library("stars")   # spatial-temporal data
# library("terra")   # raster data handling 
# library("ggplot2") # plotting
# mapviewOptions(fgb = FALSE)
# 
# meuse <- ips_sum_15 #read_sf('data/meuse.gpkg')
# v.eye <- eyefit(variog(as.geodata(meuse["sum_beetle"]), max.dist = 10000))
# ve.fit <- as.vgm.variomodel(v.eye[[1]])


# from WSL

#library(gstat)
# https://gsp.humboldt.edu/olm/R/04_01_Variograms.html

#library(sp)
data(meuse)
# no trend:
coordinates(meuse) = ~x+y
plot(variogram(log(zinc)~1, meuse))



# Save outputs ------------------------------------------------------------



save(lisa_merged_df,
  p_lisa_sub, 
     p_lisa_all,
     p_lisa_freq,
     glob_merged,
     semi_out,
     p_glob_moran_bar,
     file = "outData/lisa.Rdata")




