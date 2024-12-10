
# get variograms for all indicators of beetle dynamics
# adress better spatial synchonization

rm(list=ls()) 
gc()

# get libs ----------------------------------------------------------------

library(sf)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)
#library(rgdal) # for sp data
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


# Merge beetle dynamics inicators with XYs
dat_dynamics_xy <- 
  dat_dynamics %>% 
   # remove years woith NAs:
  dplyr::filter(year %in% 2015:2021) %>% 
  mutate(log_sum_ips = log(sum_ips),
         log_peak_diff = log(peak_diff)) %>% 
  na.omit() # as agg_doy have some NAs




# try variogram ---------------------

cutoff = 300000


# Generalized function to calculate variograms for any input data
get_variogram <- function(data, group_value, group_column, var_name, cutoff = 300000) {
  # Filter the input data based on the specified group column and value
  ips_sum_sub <- data %>%
    dplyr::filter(!!rlang::sym(group_column) == group_value)  # Filter dynamically
  
  # Create a spatial points data frame
  coordinates(ips_sum_sub) <- ~ x + y
  
  # Dynamically create the formula
  formula <- as.formula(paste(var_name, "~ 1"))
  
  # Calculate the empirical variogram
  variogram_raw <- variogram(formula, cutoff = cutoff, data = ips_sum_sub)  # all directions are assumed equal
  
  # Add group value and variable name for reference
  variogram_raw[[group_column]] <- group_value
  variogram_raw$variable <- var_name
  
  # Return the variogram data frame
  return(variogram_raw)
}


# List of years
years <- 2015:2021

# List of dependent variables
dependent_vars <- c("log_sum_ips",  "tr_agg_doy", "tr_peak_doy", "log_peak_diff", "spei12", "veg_tmp")

# Create a data frame of all combinations of years and variables;
# to run lapply afterwards
combinations <- expand.grid(year = years, 
                            variable = dependent_vars, 
                            stringsAsFactors = FALSE)

# Update get_variogram call to include the data and group column
variogram_results_ls <- lapply(1:nrow(combinations), function(i) {
  get_variogram(
    data = dat_dynamics_xy,                     # Pass your data table
    group_value = combinations$year[i],         # Use the year from combinations
    group_column = "year",                      # Specify the group column
    var_name = combinations$variable[i],        # Use the variable from combinations
    cutoff = 300000                             # Set the cutoff distance
  )
})

# Combine all results into a single data frame
combined_variogram_results <- do.call(rbind, variogram_results_ls)

# Preview the results
head(combined_variogram_results2)



# Combine all results into a single data frame
combined_variogram_results <- do.call(rbind, variogram_results_ls)

combined_variogram_results$drought_status <- ifelse(combined_variogram_results$year %in% 2018:2020, "Drought", "No Drought")

# plot years
ggplot(combined_variogram_results, aes(x = dist, y = gamma, 
                                       group = factor(year), 
                                       color = factor(year))) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'loess') + 
  facet_wrap(.~variable, scales = 'free')

# test plot over drought period
ggplot(combined_variogram_results, aes(x = dist, y = gamma, 
                                       group = factor(drought_status), 
                                       color = factor(drought_status))) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'loess') + 
  facet_wrap(.~variable, scales = 'free')







# fit variogram per years ------------------


# Define the fit_variogram function
fit_variogram <- function(variogram_data) {
  #
 # variogram_data <- variogram_results_ls[[25]]
  year = unique(variogram_data$year)
  variable = unique(variogram_data$variable)
  print(variable)
  print(year)
  
  # Calculate the Nugget (gamma at the smallest distance)
  nugget <- mean(variogram_data$gamma[1:2])
  
  # Estimate the Sill (average gamma at the highest distances, e.g., last 3-5 points)
  sill <- mean(variogram_data$gamma[(nrow(variogram_data)-2):nrow(variogram_data)])
  
  # Estimate the Range (distance where the variogram levels off, e.g., 
  # first row where gamma is near the sill)
  range <- variogram_data$dist[which.max(variogram_data$gamma >= sill * 0.95)]
  
  # Specify the initial variogram model
  initial_model <- vgm(psill = sill - nugget, model = "Mat", range = range, nugget = nugget)
  
  # Fit the variogram model to the empirical data
  fitted_variogram <- fit.variogram(variogram_data, model = initial_model) #   initial_model_peak_diff
  plot(variogram_data,fitted_variogram, cutoff = cutoff)
  
  # extract the parameters
  fitted_nugget <- fitted_variogram$psill[1]
  fitted_sill <- sum(fitted_variogram$psill)
  fitted_range <- fitted_variogram$range[2]
  
  
  # Generate a sequence of distances for prediction
  dist_seq <- seq(nugget, max(variogram_data$dist), length.out = 100)
  
  # Predict the semivariance values using the fitted model
  predicted_gamma <- variogramLine(fitted_variogram, dist_vector = dist_seq)$gamma
  
  # Create a data frame for plotting
  df <- data.frame(
    year = year,
    variable = variable,
    dist = c(variogram_data$dist, dist_seq),
    gamma = c(variogram_data$gamma, predicted_gamma),
    type = rep(c("Empirical", "Fitted"), c(nrow(variogram_data), length(dist_seq))),
    fitted_nugget = fitted_nugget,
    fitted_sill = fitted_sill,
    fitted_range = fitted_range
  )
  
  return(df)
}




# Example usage
# Assuming variogram_results_ls is your list of variogram data
# Apply the function to each variogram dataset in the list
results_years <- lapply(variogram_results_ls, fit_variogram)

# Combine the results into a single data frame
results_years <- do.call(rbind, results_years)

# make a tab;e for silt, nugget and range 
out_tab_variograms_years <- results_years %>% 
  dplyr::select(c(year,    variable , fitted_nugget, fitted_sill, fitted_range)) %>% 
  distinct() %>% 
  mutate(fitted_range = case_when(fitted_range > 300000 ~ 300000,
                                      TRUE~fitted_range)) %>% 
  dplyr::rename(nugget = fitted_nugget,
                sill = fitted_sill,
                range = fitted_range)

sjPlot::tab_df(out_tab_variograms,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/variogram_year.doc",
               digits = 3) 







# Plot years ------------------
combined_results_sub <- combined_results %>% 
 # mutate(year_color = ifelse(year %in% 2018:2020, "red", "grey")) %>% 
  dplyr::filter(dist <300000) %>% 
  mutate(variable = factor(variable, 
                           levels = c('log_sum_ips', 'tr_agg_doy', 'veg_tmp', 
                                      'tr_peak_doy', 'log_peak_diff', 'spei12'),
                           labels = c('Population level\n[#]', 'Aggregation timing\n[DOY]', expression("Temperature [°C]"),
                                      'Peak swarming\ntiming [DOY]', "Peak swarming\nintensity [#]",
                                       'SPEI [dim.]'))) %>%
  mutate(drought = factor(ifelse(year %in% 2018:2020, "yes", "no"))) #%>%

# Reverse the color palette and map to the species in the desired order

my_colors <- c(
  "grey90",  # Light gray (2015)
  "grey75",  # Medium gray (2016)
  "grey60",  # Dark gray (2017)
  "#A50026",  # Bright yellow-green (2018, emphasized)
  "#FDBE6E",  # Bright green (2019, emphasized)
  "#3366CC",  # Turquoise-green (2020, emphasized)
  "grey50"   # Charcoal gray (2021)
)
p_vario_color <- 
  combined_results_sub %>% 
    ggplot(aes(x = dist/1000, 
             y = gamma, 
             color = factor(year),
             size = drought
            # linewidth = drought
            # group = factor(year)
             )) +
  geom_point(data = subset(combined_results_sub, type == "Empirical"), 
             alpha = 0.9) +
  geom_line(data = subset(combined_results_sub, type == "Fitted"),
            aes(linewidth = drought)) +
  labs(#title = "Empirical and Fitted Variogram", 
       x = "Distance [km]", 
       y = "Semivariance",
       color = "Year"
       ) +
  scale_color_manual(values = my_colors) +
  scale_linewidth_manual(values = c("no" = 0.4, "yes" = 0.8
                                    ),
                         name = "Drought",                     # Legend title
                         labels = c("no" = "No Drought", "yes" = "Drought")) +  # Custom labels) +
  scale_size_manual(values = c("no" = 0.2, "yes" = 0.6),
                    name = "Drought",                     # Legend title
                    labels = c("no" = "No Drought", "yes" = "Drought")) +  # Custom labels) +
  #scale_color_brewer(palette = "Set1", direction = -1) +
  theme_minimal(base_size = 10) +
   theme(aspect.ratio = 1, 
         legend.position = 'right',
         legend.title =element_blank() #,
         #panel.border = element_rect(color = "black")
         ) +
  facet_wrap(variable~., scales = 'free') +
   # theme_classic()+
  theme_minimal(base_size = 10) +
  theme(
    aspect.ratio = 1, 
    legend.position = 'right',
    legend.title = element_blank(),
    panel.background = element_rect(fill = "grey99", color = 'black'), # Light grey background
    strip.background = element_blank(), # Remove boxes around facet labels
    strip.text = element_text(face = "plain"), # Optional: make facet labels bold
    panel.grid.major = element_blank(), #element_line(color = "white", linetype = "solid"), # Show major grid lines
    panel.grid.minor = element_blank(),                                    # Hide minor grid lines
    
  ) 
p_vario_color
ggsave(filename = 'outFigs/p_vario_color.png', plot = p_vario_color, 
       width = 7, height = 5, dpi = 300, 
       bg = 'white')



fwrite(combined_results, 'outTable/variogram.csv')


## test variogram years: correlation climate vs beetles  -------------------------------------------

# Combine all variogram results into one data frame
variogram_results <- bind_rows(variogram_results_ls, .id = "list_id")


# Filter relevant years and variables
drought_2018 <- variogram_results %>%
  dplyr::filter(year == 2018, variable == "spei12") %>% 
  mutate(n = 1:15)

beetles_2019 <- variogram_results %>%
  dplyr::filter(year == 2019, variable %in% c("log_sum_ips", "tr_agg_doy", "tr_peak_doy", "log_peak_diff")) %>% 
  mutate(n = rep(1:15, 4))


# Filter relevant years and variables
# Filter relevant columns and rename
drought_2018_filtered <- drought_2018 %>%
  dplyr::select(n, gamma) %>%
  dplyr::rename(gamma_spei12 = gamma)

beetles_2019_filtered <- beetles_2019 %>%
  dplyr::select(n,variable,  gamma) %>%
  dplyr::rename(gamma_beetle = gamma)

# Merge drought (2018) and beetles (2019) by distance
# Merge on 'dist' (common distance column)
combined_data <- beetles_2019_filtered %>% 
  left_join(drought_2018_filtered, by = join_by(n)) 

print(combined_data)

# Correlation analysis grouped by beetle variable
cor_results <- combined_data %>%
  dplyr::filter(variable == 'tr_peak_doy') %>% 
  group_by(variable) %>%
  summarise(
    gamma_corr = cor(gamma_spei12, gamma_beetle, use = "complete.obs", method = 'spearman'),
    .groups = "drop"
  )


# Initialize a data frame to store results
cor_results_table <- data.frame(
  variable = character(),
  pearson_corr = numeric(),
  spearman_corr = numeric(),
  stringsAsFactors = FALSE
)

# Unique variables
unique_variables <- unique(combined_data$variable)

# Loop over variables
for (var in unique_variables) {
  # Filter data for the current variable
  filtered_data <- combined_data %>%
    dplyr::filter(variable == var)
  
  # Compute Pearson and Spearman correlations
  pearson_corr <- cor(filtered_data$gamma_spei12, filtered_data$gamma_beetle, 
                      use = "complete.obs", method = "pearson")
  spearman_corr <- cor(filtered_data$gamma_spei12, filtered_data$gamma_beetle, 
                       use = "complete.obs", method = "spearman")
  
  # Append results to the table
  cor_results_table <- rbind(
    cor_results_table,
    data.frame(
      variable = var,
      pearson_corr = pearson_corr,
      spearman_corr = spearman_corr,
      stringsAsFactors = FALSE
    )
  )
}

# Print the results
print(cor_results_table)




print(cor_results)

# Scatterplot for each beetle variable
ggplot(combined_data, aes(x = gamma_spei12, y = gamma_beetle)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~variable, scales = "free") +
  labs(
    title = "Correlation Between SPEI12 (2018) and Beetle Synchronization (2019)",
    x = "Empirical Gamma (SPEI12 2018)",
    y = "Empirical Gamma (Beetle Variables 2019)"
  )



### Range plotting: -------------------------------------------------------------------
windows(6,6)
p_range <- combined_results_sub %>% 
  mutate(cap_fitted_range = case_when(fitted_range > 200000 ~ 200000,
                                      TRUE~fitted_range)) %>% 
 # filter(variable == 'log_sum_ips') %>% 
  ggplot(aes(x = year,
             y = cap_fitted_range/1000)) +
  geom_segment( aes(x=year, xend=year, y=0, yend=cap_fitted_range/1000),
                color="grey") +
  geom_point( size=2, 
              color="black", 
              #fill=alpha("orange", 3), 
              alpha=0.7, shape=16, stroke=0.5) +
  facet_wrap(.~variable, scales = 'free') +
  theme_minimal() +
  theme(aspect.ratio = 1, 
        legend.position = 'right',
        legend.title =element_blank()  ,
        panel.border = element_rect(color = "black", fill = NA)) +
  labs(y = 'Semivariogram range [km]',
       x = '')
  

ggsave(filename = 'outFigs/variogram_range.png', plot = p_range, width = 6, height = 6, dpi = 300, bg = 'white')




# Test variogram for groupped years: drought vs no drought -----------------------------------------
# Group data by drought period and aggregate
dat_dynamics_xy_grouped <- dat_dynamics_xy %>%
  mutate(drought_status = ifelse(year %in% 2018:2020, "Drought", "No Drought")) %>%
  group_by(x, y, drought_status) %>%
  summarise(across(c(log_sum_ips, tr_agg_doy, tr_peak_doy, log_peak_diff, spei12, veg_tmp), mean, na.rm = TRUE)) %>%
  ungroup()


cutoff = 300000


# for all variables:
get_variogram <- function(i, var_name, cutoff = 300000) {
  #i = 2015
  #var_name = 'log_sum_ips'
  # Slice specific year
  ips_sum_sub <- dat_dynamics_xy_grouped %>%
    dplyr::filter(drought_status == i)
  
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

# List of drought groups
drought_groups <- c("Drought", "No Drought")

# Create a data frame of all combinations of years and variables;
# to run lapply afterwards
combinations <- expand.grid(drought_group = drought_groups, 
                            variable = dependent_vars, 
                            stringsAsFactors = FALSE)


# Use lapply to apply the function to each combination
variogram_results_ls <- lapply(1:nrow(combinations), function(i) {
  get_variogram(combinations$drought_status[i], combinations$variable[i], cutoff = 300000)
})




# Apply function for each drought group and variable
grouped_variogram_results_ls <- lapply(drought_groups, function(group) {
  lapply(dependent_vars, function(var) {
    get_variogram(group, var)
  })
})

# Combine results into a single data frame
grouped_variogram_results <- do.call(rbind, do.call(c, grouped_variogram_results_ls))



# Define the fit_variogram function
#fit_variogram_group <- function(variogram_data) {
  #
  variogram_data <- grouped_variogram_results_ls[[1]]
  
  year = unique(variogram_data$year)
  variable = unique(variogram_data$variable)
  print(variable)
  print(year)
  
  
  # Calculate the Nugget (gamma at the smallest distance)
  nugget <- mean(variogram_data$gamma[1:2])
  
  # Estimate the Sill (average gamma at the highest distances, e.g., last 3-5 points)
  sill <- mean(variogram_data$gamma[(nrow(variogram_data)-2):nrow(variogram_data)])
  
  # Estimate the Range (distance where the variogram levels off, e.g., 
  # first row where gamma is near the sill)
  range <- variogram_data$dist[which.max(variogram_data$gamma >= sill * 0.95)]
  
  # Specify the initial variogram model
  initial_model <- vgm(psill = sill - nugget, model = "Mat", range = range, nugget = nugget)
  
  # Fit the variogram model to the empirical data
  fitted_variogram <- fit.variogram(variogram_data, model = initial_model) #   initial_model_peak_diff
  plot(variogram_data,fitted_variogram, cutoff = cutoff)
  
  # extract the parameters
  fitted_nugget <- fitted_variogram$psill[1]
  fitted_sill <- sum(fitted_variogram$psill)
  fitted_range <- fitted_variogram$range[2]
  
  
  # Generate a sequence of distances for prediction
  dist_seq <- seq(nugget, max(variogram_data$dist), length.out = 100)
  
  # Predict the semivariance values using the fitted model
  predicted_gamma <- variogramLine(fitted_variogram, dist_vector = dist_seq)$gamma
  
  # Create a data frame for plotting
  df <- data.frame(
    drought_status = drought_status,
    variable = variable,
    dist = c(variogram_data$dist, dist_seq),
    gamma = c(variogram_data$gamma, predicted_gamma),
    type = rep(c("Empirical", "Fitted"), c(nrow(variogram_data), length(dist_seq))),
    fitted_nugget = fitted_nugget,
    fitted_sill = fitted_sill,
    fitted_range = fitted_range
  )
  
  return(df)
}




# Example usage
# Assuming variogram_results_ls is your list of variogram data
# Apply the function to each variogram dataset in the list
results_groups <- lapply(grouped_variogram_results_ls, fit_variogram)

# Combine the results into a single data frame
results_groups<- do.call(rbind, results_groups)

# make a tab;e for silt, nugget and range 
out_tab_variograms_years <- results_years %>% 
  dplyr::select(c(year,    variable , fitted_nugget, fitted_sill, fitted_range)) %>% 
  distinct() %>% 
  mutate(fitted_range = case_when(fitted_range > 300000 ~ 300000,
                                  TRUE~fitted_range)) %>% 
  dplyr::rename(nugget = fitted_nugget,
                sill = fitted_sill,
                range = fitted_range)


## PPLOT variogram by groups -----------------------------------------------------

# Adjust the variable names for plotting (optional)
grouped_variogram_results$variable <- factor(grouped_variogram_results$variable, 
                                             levels = c('log_sum_ips', 'tr_agg_doy', 'tr_peak_doy', 
                                                        'log_peak_diff', 'spei12', 'veg_tmp'),
                                             labels = c('Population level\n[#]', 'Aggregation timing\n[DOY]', 
                                                        'Peak swarming\ntiming [DOY]', 
                                                        'Peak swarming\nintensity [#]', 
                                                        'SPEI [dim.]', 'Temperature [°C]'))

# Adjusted grouped data for plotting
grouped_variogram_results <- grouped_variogram_results %>%
  dplyr::filter(dist < 300000) %>%
  mutate(variable = factor(variable, 
                           levels = c('log_sum_ips', 'tr_agg_doy', 'veg_tmp', 
                                      'tr_peak_doy', 'log_peak_diff', 'spei12'),
                           labels = c('Population level\n[#]', 'Aggregation timing\n[DOY]', 
                                      expression("Temperature [°C]"), 'Peak swarming\ntiming [DOY]', 
                                      "Peak swarming\nintensity [#]", 'SPEI [dim.]')),
         drought_status = factor(drought_status, levels = c("No Drought", "Drought")))

# Plot with points and lines
#p_drought_variogram <- 
grouped_variogram_results %>%
  ggplot(aes(x = dist/1000, 
             y = gamma, 
             color = drought_status,
             size = drought_status
             # linewidth = drought
             # group = factor(year)
  )) +
  geom_point(data = subset(grouped_variogram_results, type == "Empirical"), 
             alpha = 0.9) +
  geom_line(data = subset(grouped_variogram_results, type == "Fitted"),
            aes(linewidth = drought_status))# +

geom_point(data = subset(grouped_variogram_results, type == "Empirical"), 
           alpha = 0.8, size = 1.5) +
  geom_line(data = subset(grouped_variogram_results, type == "Fitted"), 
            aes(linewidth = drought_status), alpha = 0.9) +
  scale_color_manual(values = c("No Drought" = "#377EB8", "Drought" = "#E41A1C")) +
  scale_linewidth_manual(values = c("No Drought" = 0.6, "Drought" = 1.0),
                         name = "Drought Status",
                         labels = c("No Drought", "Drought")) +
  labs(x = "Distance [km]",
       y = "Semivariance",
       color = "Drought Status",
       title = "Empirical and Fitted Variograms by Drought Status") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    panel.background = element_rect(fill = "grey95", color = "black")
  )

# Display the plot
print(p_drought_variogram)
# Save the plot
ggsave(filename = "outFigs/drought_variograms_points_lines.png", 
       plot = p_drought_variogram, 
       width = 8, height = 6, dpi = 300, bg = "white")


# Display the plot
print(p_drought_variogram)



# END



