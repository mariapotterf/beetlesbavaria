

# try synchrony analysis: correlogram space and time between beetle populations?

# correlogram analysis: space and time synchrony:
# spatial covariance functions 

library(readr)
library(dplyr)
library(ncf)
library(ggplot2)

# Vars ----------
drought_years <- c(2018, 2019, 2020)  # Add actual drought years


# read two input data type: -------------------------------------------------
# on yearly level (7 observations/ trap)
# on 
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

str(dat_doy)
head(dat_doy)


# Remove NA values for sum_ips
dat_clean <- dat %>%
  dplyr::filter(!is.na(sum_ips)) 


# get values on pair level, not on a trap level ------------
# Compute mean values at the trap pair level
dat_pairID <- dat_clean %>% 
  group_by(pairID, year) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


# what is teh average distance between two traps? 

# Load necessary package
library(stats)  # dist() function is from base R


# Compute pairwise Euclidean distances
distance_matrix <- as.matrix(dist(cbind(x_coords, y_coords), method = "euclidean"))

# Extract only lower triangle values (excluding diagonal)
pairwise_distances <- distance_matrix[lower.tri(distance_matrix, diag = FALSE)]

# Compute the average Euclidean distance
average_distance <- mean(pairwise_distances)

# 152 km between all trap pairs

###  get average nearest neighbour distance: -----------------

# Replace diagonal values (self-distances) with a large number so they aren't selected
diag(distance_matrix) <- Inf  

# Find the nearest neighbor for each trap
nearest_distances <- apply(distance_matrix, 1, min)

# Compute the average nearest neighbor distance
average_nn_distance <- mean(nearest_distances)

average_nn_distance



# Correlogram ---------------------------------------------
# run on pair level
dat_slice = dat_pairID %>% 
  dplyr::filter(year == 2018)
# Check unique years ----------------------------------------
unique_years <- unique(dat_slice$year)
print(unique_years)  # Ensure you have 7 years

# Extract required variables
x_coords <- dat_slice$x
y_coords <- dat_slice$y
time_var <- dat_slice$year
sum_ips_values <- dat_slice$sum_ips



# Compute spatial correlogram (Moran's I)
correlog_result <- correlog(x=x_coords, 
                            y=y_coords, 
                            z=sum_ips_values, 
                            increment=30000, resamp=1000)

# Plot the synchrony matrix
plot(correlog_result)


# run separate correlograms per year --------------

# Get unique years
years <- unique(dat_pairID$year)

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
                             increment=20000, 
                             resamp=500)
    
    # Store results in a data frame
    correlog_df <- data.frame(
      year = yr,
     # n = correlog_res$mean.of.class,
      distance = correlog_res$mean.of.class,
      correlation = correlog_res$correlation,
      p_val = correlog_res$p
    ) %>% 
      mutate(row_number = 1:n(),
             significant = ifelse(p_val < 0.05, "Yes", "No") ) 
    
    # Append to list
    correlog_list[[as.character(yr)]] <- correlog_df
  }
}

# Combine all years into one data frame
correlog_data <- bind_rows(correlog_list)

# calculate teh average across bins
corr_data_drought_avg <- correlog_data %>%
  dplyr::filter(distance < 160000) %>% 
  mutate(drought_status = ifelse(year %in% drought_years, "Drought", "Non-Drought")) %>% 
  group_by(drought_status, row_number  ) %>%
  summarise(
    year= mean(year, na.rm = T),
    distance = mean(distance, na.rm = TRUE),
    mean_corr = mean(correlation, na.rm = TRUE),
    sd_corr = sd(correlation, na.rm = TRUE),  # Standard deviation
    lower_ci = quantile(correlation, 0.05, na.rm = TRUE),  # 5th percentile
    upper_ci = quantile(correlation, 0.95, na.rm = TRUE),  # 95th percentile
    p_val  = mean(p_val , na.rm = TRUE)
  ) %>% 
  mutate(significant = ifelse(p_val < 0.05, "Yes", "No") ) 
  #dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE))

  

# Plot using ggplot with facets
corr_data_drought_avg %>% 
  ggplot( aes(x = distance, y = mean_corr,  fill = drought_status)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  geom_point( aes(size = significant,  color = drought_status)) +
  geom_line(aes(color = drought_status)) + 
 # geom_point() +
 #facet_wrap(~ year, scales = "free_x") +  # Separate plots per year
  labs(title = "Correlogram by Year", x = "Distance", y = "Correlation") +
  theme_classic2()


# run correlogram for all variables ---------------------------------------



# pairwise correlation beetle-beetle: DOY -------------------------------------------

# temporal corss-correlation: averaging withn years and DOY ------------------------

# Ensure correct column types
dat_doy <- dat_doy %>%
  mutate(year = as.numeric(year),
         doy = as.numeric(doy),
         count = as.numeric(count),
         trapID = as.factor(trapID),
         species = as.factor(species))


# make on teh trap pair level: 
dat_doy_pairID <- dat_doy %>% 
  group_by(pairID, year, doy) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE))


## analyse all data across years ---------------------
# Convert to time series
ts_data <- ts(dat_doy_pairID$count, start = min(dat_doy_pairID$year), frequency = 365)

# Compute Cross-Correlation Function (CCF):
# check how beetle data are correlated over time: eg if counts today qaffect counts in later day (lag)
ccf_result <- ccf(ts_data, ts_data, lag.max = 30, na.action = na.omit)

# Plot the CCF
plot(ccf_result, main="Cross-Correlation Function of Beetle Counts", ylab="CCF", xlab="Lag (days)")     


## differentiate betweend rought and no drought years ------------

# Create separate datasets for drought and non-drought years
drought_years <- c(2018, 2019, 2020)
non_drought_years <- c(2015, 2016, 2017, 2021)

# Filter for drought years
dat_drought <- dat_doy_pairID %>%
  dplyr::filter(year %in% drought_years)

# Filter for non-drought years
dat_no_drought <- dat_doy_pairID %>%
  dplyr::filter(year %in% non_drought_years)

# Convert both datasets to time series
ts_drought <- ts(dat_drought$count, start = min(dat_drought$year), frequency = 365)
ts_no_drought <- ts(dat_no_drought$count, start = min(dat_no_drought$year), frequency = 365)

# Compute and plot CCF for drought years
ccf_drought <- ccf(ts_drought, ts_drought, lag.max = 30, na.action = na.omit)
plot(ccf_drought, main="CCF of Beetle Counts (Drought Years)", ylab="CCF", xlab="Lag (days)")

# Compute and plot CCF for non-drought years
ccf_no_drought <- ccf(ts_no_drought, ts_no_drought, lag.max = 30, na.action = na.omit)
plot(ccf_no_drought, main="CCF of Beetle Counts (Non-Drought Years)", ylab="CCF", xlab="Lag (days)")

#### compare temporal synchronization: drought vs no drought years: ----------

# Function to extract cross-correlation values
extract_ccf_values <- function(ccf_result, drought_label) {
  data.frame(
    lag = ccf_result$lag,
    correlation = abs(ccf_result$acf),  # Use absolute values to capture synchronization strength
    condition = drought_label  # Label for drought vs. non-drought
  )
}

# Extract values from drought CCF
drought_ccf_data <- extract_ccf_values(ccf_drought, "Drought")

# Extract values from non-drought CCF
non_drought_ccf_data <- extract_ccf_values(ccf_no_drought, "Non-Drought")

# Combine both datasets
ccf_summary_data <- rbind(drought_ccf_data, non_drought_ccf_data)

# select only specific lags:
#ccf_summary_data <- ccf_summary_data %>% 
#  dplyr::filter(lag %in% c(0, 5, 10, -5, -10)) #10,, -10

# Create boxplot comparing synchronization across conditions
ggplot(ccf_summary_data, aes(x = condition, y = correlation, fill = condition)) +
  geom_boxplot(notch = T) +
  labs(title = "Comparison of Temporal Synchronization",
       x = "Condition",
       y = "Absolute Cross-Correlation (CCF)") +
  theme_minimal()



### CCF temporal synchronization per year -------------
library(purrr)

# Extract unique years
years <- unique(dat_doy_pairID$year)

# Function to compute CCF for each year
compute_ccf_for_year <- function(year) {
  # Filter data for the specific year
  dat_year <- dat_doy_pairID %>% filter(year == !!year)
  
  # Convert to time series
  ts_year <- ts(dat_year$count, start = min(dat_year$year), frequency = 365)
  
  # Compute CCF
  ccf_year <- ccf(ts_year, ts_year, lag.max = 30, na.action = na.omit, plot = FALSE)
  
  # Extract relevant CCF values
  data.frame(
    lag = ccf_year$lag,
    correlation = abs(ccf_year$acf),  # Use absolute values to measure synchronization strength
    year = as.factor(year)  # Convert year to a categorical variable for plotting
  )
}

# Compute CCF for each year and combine results
ccf_results_per_year <- map_dfr(years, compute_ccf_for_year)

# Ensure 'year' is treated as a factor
ccf_summary_data <- ccf_results_per_year %>%
  mutate(year = as.factor(year))  # Convert year to categorical factor

ccf_summary_data <- ccf_summary_data %>%
  dplyr::filter(lag > -0.02 & lag < 0.02)

nrow(ccf_summary_data)
# Create boxplot comparing synchronization across years
ggplot(ccf_results_per_year , aes(x = year, y = correlation, fill = year)) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Comparison of Temporal Synchronization Across Years",
       x = "Year",
       y = "Absolute Cross-Correlation (CCF)")  +  # Remove legend since x-axis represents years
 # stat_compare_means(method = "kruskal.test", 
     #                label.y = max(ccf_summary_data$correlation) * 1.05) + # Overall test
#  stat_compare_means(comparisons = list(c("2015", "2016"), c("2018", "2019"), c("2020", "2021")), 
      #               method = "wilcox.test", label = "p.signif")+  # Pairwise comparisons
 ylim(0,0.5) + 
   theme_classic2() +
  theme(legend.position = "none")  # Remove legend since x-axis represents years





# spatial cross correlation analysis:drought year level ---------------------------------------------

# Define drought and non-drought years
drought_years <- c(2018) # , 2019, 2020
non_drought_years <- c(2016) #2015,   2017, 2021


# Function to compute spatial correlations for a given period (drought or non-drought)
compute_spatial_corr_for_period <- function(years_group, label) {
 # year = 2018
  dat_subset <- 
    dat_clean %>%
    #dplyr::filter(year %in% drought_years)  %>% 
    dplyr::filter(year %in% years_group) %>% #drought_years
   # dplyr::select(year, pairID, sum_ips) %>%
    dplyr::select(year, pairID, sum_ips) %>%
    pivot_wider(names_from = pairID, values_from = sum_ips) %>%
    unnest(cols = -year) 
   # mutate(across(-year, ~ map(.x, ~ .x %>% unlist() %>% as.numeric()), .names = "{.col}")) #%>%
    #pivot_wider(names_from = pairID, values_from = sum_ips) #%>%
  #  drop_na()  # Remove rows with missing data
  
  pair_names <- colnames(dat_subset)[-1]  # Exclude 'year' column
  
  # Function to compute correlation for each trap pair
  compute_corr <- function(trap1, trap2) {
    cor(dat_subset[[trap1]], 
        dat_subset[[trap2]], use = "complete.obs", method = "spearman")
  }
  
  # Compute correlation for all trap pairs
  spatial_corr_results <- expand.grid(pair1 = pair_names, pair2 = pair_names) %>%
    dplyr::filter(pair1 != pair2) %>%
    mutate(corr_value = map2_dbl(pair1, pair2, compute_corr),
           period = label)  # Add drought vs. non-drought label
  
  return(spatial_corr_results)
}

# Compute correlations for drought and non-drought years
spatial_corr_drought <- compute_spatial_corr_for_period(drought_years, "Drought")
spatial_corr_non_drought <- compute_spatial_corr_for_period(non_drought_years, "Non-Drought")

# Combine results
spatial_corr_results <- bind_rows(spatial_corr_drought, spatial_corr_non_drought)

# Create boxplot comparing drought vs. non-drought spatial synchronization
ggplot(spatial_corr_results, aes(x = period, y = corr_value, fill = period)) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Comparison of Spatial Synchronization: Drought vs. Non-Drought",
       x = "Period",
       y = "Spearman Correlation (Spatial Synchronization)") +
  theme_minimal() +
  theme(legend.position = "none")




# include spatial correlation with distance? -------------------

# Define drought and non-drought years
drought_years <- c(2018, 2019, 2020)
non_drought_years <- c(2015, 2016, 2017, 2021)

# Compute mean values at the trap pair level
dat_pairID <- dat_clean %>% 
  group_by(pairID, year) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Function to compute spatial correlations for a given period (drought or non-drought)
compute_spatial_corr_for_period <- function(years_group, label) {
  dat_subset <- dat_pairID %>%
    filter(year %in% years_group) %>%
    select(year, pairID, sum_ips) %>%
    pivot_wider(names_from = pairID, values_from = sum_ips) %>%
    drop_na()  # Remove rows with missing data
  
  pair_names <- colnames(dat_subset)[-1]  # Exclude 'year' column
  
  # Function to compute correlation for each trap pair
  compute_corr <- function(trap1, trap2) {
    cor(dat_subset[[trap1]], dat_subset[[trap2]], use = "complete.obs", method = "spearman")
  }
  
  # Compute correlation for all trap pairs
  spatial_corr_results <- expand.grid(pair1 = pair_names, pair2 = pair_names) %>%
    filter(pair1 != pair2) %>%
    mutate(corr_value = map2_dbl(pair1, pair2, compute_corr),
           period = label)  # Add drought vs. non-drought label
  
  return(spatial_corr_results)
}

# Compute correlations for drought and non-drought years
spatial_corr_drought <- compute_spatial_corr_for_period(drought_years, "Drought")
spatial_corr_non_drought <- compute_spatial_corr_for_period(non_drought_years, "Non-Drought")

# Combine results
spatial_corr_results <- bind_rows(spatial_corr_drought, spatial_corr_non_drought)

# Create boxplot comparing drought vs. non-drought spatial synchronization
ggplot(spatial_corr_results, aes(x = period, y = corr_value, fill = period)) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Comparison of Spatial Synchronization: Drought vs. Non-Drought",
       x = "Period",
       y = "Spearman Correlation (Spatial Synchronization)") +
  theme_minimal() +
  theme(legend.position = "none")
