

# try synchrony analysis: correlogram space and time between beetle populations?

# correlogram analysis: space and time synchrony:
# spatial covariance functions 

library(readr)
library(dplyr)
library(ncf)
library(ggplot2)

dat <- read_csv("outTable/beetle_dynamic_indicators.csv") # yearly level
dat_doy <- read_csv("outTable/dat_DOY_counts_IPS.csv")

# get uinique coordinates: if teh trap shifted, use average of teh coordinate per trap
dat_xy <- dat %>% 
  dplyr::select(trapID, x, y) %>% 
  distinct() %>% 
  group_by(trapID) %>% 
  summarise(x = mean(x, na.rm =T),
            y = mean(y, na.rm =T))

# add coordinates to DOY data
dat_doy <- dat_doy %>% 
  left_join(dat_xy, by = join_by(trapID))
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



# Investigate cross-correlation structure: dummy example ---------------------------------------------
# ips counts and climate predictors ------------------------------------
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
head(dat)


library(ncf)  # For spatial cross-correlation analysis

# Define spatial cross-correlation function
spatial_ccf <- function(df, variable1, variable2, distance_increment = 500) {
  
  # Remove NA values
  df <- df %>% filter(!is.na(!!sym(variable1)), !is.na(!!sym(variable2)))
  
  # Run cross-correlation at increasing spatial distances
  correlog_res <- correlog(
    x = df$x,
    y = df$y,
    z = df[[variable1]],
    w = df[[variable2]],  # Second variable for cross-correlation
    increment = distance_increment,
    resamp = 100
  )
  
  return(data.frame(distance = correlog_res$mean.of.class, correlation = correlog_res$correlation))
}

# Compute spatial cross-correlation for beetle population (ips_sum) vs. drought (spei)
correlog_data <- spatial_ccf(dat, "ips_sum", "spei")

# Plot the cross-correlation function
ggplot(correlog_data, aes(x = distance, y = correlation)) +
  geom_line() + 
  geom_point() +
  labs(title = "Spatial Cross-Correlation (Beetles vs. Drought)", x = "Distance", y = "Spearman Correlation") +
  theme_minimal()


# Define drought years
drought_years <- c(2018, 2019, 2020)

# Compute CCF separately for drought vs. non-drought years
correlog_drought <- spatial_ccf(filter(dat, year %in% drought_years), "ips_sum", "spei")
correlog_nondrought <- spatial_ccf(filter(dat, !year %in% drought_years), "ips_sum", "spei")

# Add labels for comparison
correlog_drought$group <- "Drought Years (2018-2020)"
correlog_nondrought$group <- "Non-Drought Years (2015-2017, 2021)"

# Combine results
correlog_comparison <- bind_rows(correlog_drought, correlog_nondrought)

# Plot comparison
ggplot(correlog_comparison, aes(x = distance, y = correlation, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Spatial Cross-Correlation: Drought vs. Non-Drought Years",
       x = "Distance", y = "Spearman Correlation") +
  theme_minimal()




# pairwise correlation beetle-beetle -------------------------------------------
library(geosphere)  # For distance calculation
library(dplyr)
library(tidyr)
library(ncf)


dat <- expand.grid(trapID = 1:10, year = 2015:2021, doy = sample(100:300, 20, replace = TRUE)) %>%
  arrange(trapID, year, doy) %>%
  mutate(
    x = runif(n(), 1000, 5000),  # Random spatial X coordinates
    y = runif(n(), 1000, 5000),  # Random spatial Y coordinates
    spei = rnorm(n(), mean = 0, sd = 1),  # Simulated drought index
    tmp = rnorm(n(), mean = 15, sd = 2),  # Temperature
    ips_sum = rpois(n(), lambda = 500) * (1 + spei * 0.2) + 50 * sin(doy / 30)  # Seasonal beetle activity
  )



# Compute pairwise distances between traps
trap_locations <- dat %>%
 # dplyr::filter(is.na(ips_sum)) %>% 
  dplyr::select(trapID, x, y) %>%
   distinct()

distance_matrix <- dist(trap_locations[, c("x", "y")])  # Euclidean distance
distance_df <- as.data.frame(as.table(as.matrix(distance_matrix)))  # Convert to dataframe
colnames(distance_df) <- c("trap1", "trap2", "distance")

# Remove self-distances (trap compared to itself)
distance_df <- dplyr::filter(distance_df, trap1 != trap2)

# Merge beetle population data for each trap
dat_wide <- dat %>%
  group_by(year, doy, trapID) %>%
  summarise(ips_sum = mean(ips_sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = trapID, values_from = ips_sum)


# Select two traps for cross-correlation
trap1 <- dat_wide$`1`
trap2 <- dat_wide$`2`

# Compute Cross-Correlation Function (CCF)
ccf_result <- ccf(trap1, trap2, lag.max = 10, plot = TRUE)

# Print correlation values at different time lags
print(ccf_result)


ggplot(dat, aes(x = doy, y = ips_sum, color = as.factor(trapID))) +
  geom_line(alpha = 0.5) +
  facet_wrap(~year) +
  labs(title = "Beetle Population Trends Over Time",
       x = "Day of Year (DOY)", y = "Beetle Count (ips_sum)") +
  theme_minimal()



# old

# Compute pairwise correlations between traps for each year
cor_results <- list()
for (yr in unique(dat$year)) {
  year_data <- dplyr::filter(dat_wide, year == yr)[, -1]  # Remove year column
  
  # Compute Spearman correlation between traps
  corr_matrix <- cor(year_data, use = "pairwise.complete.obs", method = "spearman")
  corr_df <- as.data.frame(as.table(corr_matrix))
  colnames(corr_df) <- c("trap1", "trap2", "correlation")
  
  # Merge with distance data
  merged_df <- merge(distance_df, corr_df, by = c("trap1", "trap2"))
  merged_df$year <- yr  # Add year column
  
  cor_results[[as.character(yr)]] <- merged_df
}

# Combine all years into one dataframe
correlation_data <- bind_rows(cor_results)

# Check structure
head(correlation_data)

