
# get variograms for all indicators of beetle dynamics
# adress better spatial synchonization


# get libs ----------------------------------------------------------------

library(sf)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)
library(rgdal) # for sp data
library(sp)
library(fasterize)
library(terra)
library(ggplot2)
library(ggpubr)
library('cowplot') # arrange plots


# Spatstat
library(spdep)
library(gstat)  # 

# load cleaned data
#load("outData/ips_counts.Rdata")
load("outData/spatial.Rdata")
dat_dynamics <- fread( 'outTable/beetle_dynamic_indicators.csv')


# Spatial data: 
sort(unique(xy_sf_expand$falsto_name))


# get coordinates from sf object
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    year = xy_sf_expand$year,
                    trapID = xy_sf_expand$falsto_name)

xy_df <- distinct(xy_df)
my_crs <- crs(xy_sf_expand)


# Merge beetle dynamics inicators with XYs
dat_dynamics_xy <- 
  dat_dynamics %>% 
    left_join(xy_df, by = c("trapID", 'year')) %>% 
  # remove years woith NAs:
  dplyr::filter(year %in% 2015:2021) %>% 
  mutate(log_sum_ips = log(sum_ips),
         log_peak_diff = log(peak_diff)) %>% 
  na.omit() # as agg_doy have some NAs




# try variogram ---------------------

cutoff = 300000


# for all variables:
get_variogram <- function(i, var_name, cutoff = 300000) {
  # Slice specific year
  ips_sum_sub <- dat_dynamics_xy %>%
    filter(year == i)
  
  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  
  # Dynamically create the formula
  formula <- as.formula(paste(var_name, "~ 1"))
  
  # Calculate the empirical variogram
  variogram_raw <- variogram(formula, cutoff = cutoff, data = ips_sum_sub)  # all directions are assumed equal
  
  # Add year and variable name
  variogram_raw$year <- i
  variogram_raw$variable <- var_name
  
  # Return the variogram data frame
  return(variogram_raw)
}


# List of years
years <- 2015:2021

# List of dependent variables
dependent_vars <- c("log_sum_ips",  "tr_agg_doy", "tr_peak_doy", "log_peak_diff")

# Initialize an empty list to store results
variogram_results <- list()

# Loop over years and dependent variables
for (year in years) {
  for (var in dependent_vars) {
    variogram_result <- get_variogram(year, var, cutoff = 300000)
    variogram_results[[paste(year, var, sep = "_")]] <- variogram_result
  }
}

# Combine all results into a single data frame
combined_variogram_results <- do.call(rbind, variogram_results)



ggplot(combined_variogram_results, aes(x = dist, y = gamma, 
                                       group = factor(year), 
                                       color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'loess') + 
  facet_wrap(.~variable, scales = 'free')








# fit variogram ------------------

variogram_data <- variogram_out[[5]]

# Calculate the Nugget (gamma at the smallest distance)
nugget <- mean(variogram_data$gamma[1:5])

# Estimate the Sill (average gamma at the highest distances, e.g., last 3-5 points)
sill <- mean(variogram_data$gamma[13:15])

# Estimate the Range (distance where the variogram levels off, e.g., first row where gamma is near the sill)
range <- variogram_data$dist[which.max(variogram_data$gamma >= sill * 0.95)]

# Display the estimates
nugget
sill
range

# Specify the initial variogram model
initial_model <- vgm(psill = sill, model = "Lin", range = range, nugget = nugget)

# Fit the variogram model to the empirical data
fitted_variogram <- fit.variogram(variogram_data, model = initial_model)

print(fitted_variogram)

# Plot the empirical and fitted variogram
plot(variogram_data, model = fitted_variogram)

