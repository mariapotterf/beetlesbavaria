# R HElp


# filetr only locations that occurs oevr all years:

# Example: keep only records that are present over all categories
df <- data.frame(loc = c('a', "b", "c",
                         'a', "b",
                         "b", "c",
                         "b", "c"
),
year = c(1,1,1,
         2,2,
         3,3,
         4,4))

# keep only rows that have the same amount of observations
loc.unique <- unique(df$loc)
year.unique <- unique(df$year)

df %>% 
  # group_by(years, loc) %>% 
  group_by(loc) %>%
  filter(n_distinct(year) == n_distinct(df$year))





##################################################################################

##################################################################################


# Coefficient of Variation

# EXAMPLE : Sample beetle count data (mean counts per year)
years <- c(2017, 2018, 2019, 2020, 2021)
mean_counts <- c(50, 60, 55, 75, 70)

# Sample standard deviation data for each year
std_dev <- c(10, 8, 12, 15, 10)

# Calculate coefficient of variation (CV) for each year
cv_values <- (std_dev / mean_counts) * 100

# Create a plot to illustrate temporal variability using CV
windows()
plot(years, cv_values, type='o', col='blue', pch=16, lty=1, xlab='Year', ylab='Coefficient of Variation (CV %)',
     main='Temporal Variability of Beetle Counts per Year')
grid()
