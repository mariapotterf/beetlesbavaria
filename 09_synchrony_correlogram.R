

# try synchrony analysis: correlogram space and time between beetle populations?

# correlogram analysis: space and time synchrony:
# spatial covariance functions 

library(readr)
library(dplyr)
library(ncf)
library(ggplot2)

source("my_functions.R")

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
  drop_na()
  #dplyr::filter(!is.na(sum_ips)) 


# get values on pair level, not on a trap level ------------
# Compute mean values at the trap pair level
dat_pairID <- dat_clean %>% 
  group_by(pairID, year) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))




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



# what is teh average distance between two traps? 
# Compute pairwise Euclidean distances
distance_matrix <- as.matrix(dist(cbind(x_coords, y_coords), method = "euclidean"))

# Extract only lower triangle values (excluding diagonal)
pairwise_distances <- distance_matrix[lower.tri(distance_matrix, diag = FALSE)]

# Compute the average Euclidean distance
average_distance <- mean(pairwise_distances)

# 152 km for pairwise distance calculation between all trap pairs

###  get average nearest neighbour distance: -----------------

# Replace diagonal values (self-distances) with a large number so they aren't selected
diag(distance_matrix) <- Inf  

# Find the nearest neighbor for each trap
nearest_distances <- apply(distance_matrix, 1, min)

# Compute the average nearest neighbor distance
average_nn_distance <- mean(nearest_distances)

average_nn_distance
# 23 km



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

# Define a function to compute correlograms
compute_correlogram <- function(data, variable, years, increment = 20000, resamp = 500) {
  
  # Initialize an empty list to store results
  correlog_list <- list()
  
  # Loop through each year and compute correlograms
  for (yr in years) {
    subset_data <- data %>% dplyr::filter(year == yr)
    
    # Run correlogram analysis only if enough data points exist
    if (nrow(subset_data) > 10) {
      correlog_res <- correlog(
        x = subset_data$x, 
        y = subset_data$y, 
        z = subset_data[[variable]],  # Use variable dynamically
        increment = increment, 
        resamp = resamp
      )
      
      # Store results in a data frame
      correlog_df <- data.frame(
        variable = variable,
        year = yr,
        distance = correlog_res$mean.of.class,
        correlation = correlog_res$correlation,
        p_val = correlog_res$p
      ) %>% 
        mutate(row_number = 1:n(),
               significant = ifelse(p_val < 0.05, "Yes", "No"))
      
      # Append to list
      correlog_list[[paste(variable, yr, sep = "_")]] <- correlog_df
    }
  }
  
  # Combine results into a single dataframe
  correlog_data <- bind_rows(correlog_list)
  
  return(correlog_data)
}

# Define variables to analyze
variables <- c("sum_ips", "agg_doy", "peak_doy", "peak_diff") #,

# Get unique years
years <- unique(dat_clean$year)

# Run the function for each variable and combine results
correlog_results <- map_dfr(variables, ~ compute_correlogram(dat_clean, .x, years))


# avegare values across bins: 

corr_data_drought_avg <- correlog_results %>%
  dplyr::filter(distance < 160000) %>% 
  mutate(drought_status = ifelse(year %in% drought_years, "Hotter drought", "Other")) %>% 
  group_by(variable, drought_status, row_number) %>%
  summarise(
    year = mean(year, na.rm = TRUE),
    distance = mean(distance, na.rm = TRUE),
    mean_corr = mean(correlation, na.rm = TRUE),
    sd_corr = sd(correlation, na.rm = TRUE),
    lower_ci = quantile(correlation, 0.05, na.rm = TRUE),
    upper_ci = quantile(correlation, 0.95, na.rm = TRUE),
    p_val = mean(p_val, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(significant = ifelse(p_val < 0.05, "Yes", "No"))


# Define factor levels for variable ordering
corr_data_drought_avg <- corr_data_drought_avg %>%
  mutate(variable = factor(variable, levels = c("sum_ips", "agg_doy", 
                                                "peak_doy", "peak_diff"),
                           labels = c(lab_popul_level, 
                                      lab_colonization_time, 
                                      lab_peak_time, 
                                      lab_peak_growth)))


# Define the color scheme based on extracted colors
color_palette <- c("Hotter drought" = "#8b0000", "Other" = "#b0b0b0")  # Red & Gray
fill_palette <- c("Hotter drought" = "#8b0000", "Other" = "#b0b0b0")   # Same for ribbon shading



# Create a vector with custom facet labels
facet_labels <- c(
  "Population level [#]" = "[a] Population level [#]", # lab_popul_level = paste("[a] ", lab_popul_level),
  "Aggregation timing [DOY]" = "[b] Aggregation timing [DOY]",
  "Peak sw. timing [DOY]" = "[c] Peak sw. timing [DOY]",
  "Peak sw. intensity [#]" = "[d] Peak sw. intensity [#]"
)
p_correlograms_lines <- corr_data_drought_avg %>% 
  ggplot(aes(x = distance/1000, y = mean_corr, fill = drought_status)) +
  geom_line(aes(color = drought_status), lwd = 0.8) +  # Line plot
  geom_point(color = 'white', size = 2) +  # Points with size for significance
  
   geom_point(aes(shape = significant, color = drought_status),  size = 1) +  # Points with size for significance
  scale_color_manual(values = color_palette) +  # Apply extracted colors to lines and points
  scale_fill_manual(values = fill_palette) +  # Apply extracted colors to shaded CIs
  scale_shape_manual(values = c("Yes" = 16, "No" = 1)) +  # 16 = closed circle, 1 = open circle 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey') +
  facet_wrap(~ variable, labeller = labeller(variable = facet_labels)) + 
  # facet_wrap(~ variable) +  # Facet by variable , scales = "free_y"
  labs(title = "",
       x = "Distance [km]",
       y = "Spatial Correlation",
       fill = "Drought Status",
       color = "Drought Status",
       shape = "Significance") +
 # theme_bw() +
  #theme(legend.position = "bottom")  + # Move legend to the right
  theme_classic(base_size = 8) +
  theme(
    aspect.ratio = 1, 
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "plain"),  # Adjusts facet title text styling
    #panel.background = element_rect(fill = NA, colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    #axis.title.y = element_blank(),
    plot.title = element_text(size = 8) # Set the plot title size to 10
  ) 
p_correlograms_lines
#p_correlograms_ribbons
ggsave(filename = 'outFigs/p_correlograms_lines.png', plot = p_correlograms_lines, 
       width = 5, height = 5, dpi = 300, 
       bg = 'white')

# ggsave(filename = 'outFigs/p_correlograms_ribbons.png', 
#        plot = p_correlograms_ribbons, 
#        width = 5, height = 5, dpi = 300, 
#        bg = 'white')

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









