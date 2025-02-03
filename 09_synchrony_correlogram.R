

# try synchrony analysis: correlogram space and time between beetle populations?

# correlogram analysis: space and time synchrony:
# spatial covariance functions 

library(readr)
library(dplyr)
library(ncf)
library(ggplot2)

dat <- read_csv("outTable/beetle_dynamic_indicators.csv")
#View(dat)

str(dat)

# Remove NA values for sum_ips
dat_clean <- dat %>%
  dplyr::filter(!is.na(sum_ips)) 



# Check unique years
unique_years <- unique(dat_clean$year)
print(unique_years)  # Ensure you have 7 years

# Extract required variables
x_coords <- dat_clean$x
y_coords <- dat_clean$y
time_var <- dat_clean$year
sum_ips_values <- dat_clean$sum_ips

# Compute spatial correlogram (Moran's I)
correlog_result <- correlog(x=x_coords, y=y_coords, z=sum_ips_values, increment=5000, resamp=100)

# Plot the synchrony matrix
plot(correlog_result)


# run separate correlograms per year --------------

# Get unique years
years <- unique(dat_clean$year)

# Create an empty list to store correlogram results
correlog_list <- list()

# Loop through each year and compute correlograms
for (yr in years) {
  subset_data <- dat_clean[dat_clean$year == yr, ]
  
  # Run correlogram analysis only if enough data points exist
  if (nrow(subset_data) > 10) {
    correlog_res <- correlog(x=subset_data$x, 
                             y=subset_data$y, 
                             z=subset_data$sum_ips, 
                             increment=5000, 
                            # max.distance = 170000, # half distance
                             resamp=100)
    
    # Store results in a data frame
    correlog_df <- data.frame(
      distance = correlog_res$mean.of.class,
      correlation = correlog_res$correlation,
      year = as.factor(yr)  # Convert year to factor for ggplot facets
    )
    
    # Append to list
    correlog_list[[as.character(yr)]] <- correlog_df
  }
}

# Combine all years into one data frame
correlog_data <- bind_rows(correlog_list)

# Plot using ggplot with facets
correlog_data %>% 
  dplyr::filter(distance < 170000) %>% 
  ggplot( aes(x = distance, y = correlation, color = year)) +
  geom_line() + 
  geom_point() +
 #facet_wrap(~ year, scales = "free_x") +  # Separate plots per year
  labs(title = "Correlogram by Year", x = "Distance", y = "Correlation") +
  theme_minimal()



library(vegan)  # For Kendall's W calculation

# Reshape data to wide format: Rows = Years, Columns = Locations
sum_ips_wide <- dat_clean %>%
  dplyr::select(year, trapID, sum_ips) %>%
  tidyr::pivot_wider(names_from = trapID, values_from = sum_ips)

# Compute Kendall's W to assess synchrony
sync_w <- kendall.global(sum_ips_wide[,-1])  # Remove year column

print(sync_w)

# Convert data to long format for ggplot
sum_ips_long <- tidyr::pivot_longer(sum_ips_wide, cols = -year, names_to = "trapID", values_to = "sum_ips")

# Plot synchrony trends across traps
ggplot(sum_ips_long, aes(x = year, y = sum_ips, color = trapID, group = trapID)) +
  geom_line(alpha = 0.5) +
  labs(title = "Beetle Outbreak Trends (sum_ips) Across Traps",
       x = "Year", y = "sum_ips") +
  theme_minimal()





# test corrplot -----------------

# correlog: estimate 



# test cross-correlation structure ---------------------------------------------

# Create a dummy dataset with spatial and temporal structure
set.seed(123)
dat <- expand.grid(trapID = 1:10, year = 2015:2021) %>%
  mutate(
    x = runif(n(), 1000, 5000),  # Random spatial X coordinates
    y = runif(n(), 1000, 5000),  # Random spatial Y coordinates
    spei = rnorm(n(), mean = 0, sd = 1),  # Simulated drought index
    tmp = rnorm(n(), mean = 15, sd = 2),  # Temperature
    ips_sum = rpois(n(), lambda = 500) * (1 + spei * 0.2)  # Simulated beetle population response to drought
  )

# Check structure
str(dat)


