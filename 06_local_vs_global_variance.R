# Compare local vs global variance

# is it changing over changing temperatures??


library(lme4)  # for mixed-effects models
library(lmerTest)  # for getting p-values in lmer models
library(ggplot2)  # for plotting


#data <- read.csv("your_data_file.csv")
#head(data)




# Create a dummay dataset
# Load necessary library
library(dplyr)

# Set parameters
num_years <- 7          # e.g., 2015-2021
num_traps <- 158        # Total number of traps
num_pairs <- 79         # Number of trap pairs
samples_per_season <- 15  # Number of samples per trap per season

# Create a sequence of years
years <- seq(2015, 2015 + num_years - 1)

# Create trap IDs and pair them
trap_ids <- sprintf("%03d", 1:num_traps)
pair_ids <- rep(sprintf("%02d", 1:num_pairs), length.out = num_traps)

# Ensure even distribution of traps in pairs (if possible)
if (length(unique(pair_ids)) < length(trap_ids)) {
  pair_ids <- c(pair_ids, rep(NA, length(trap_ids) - length(pair_ids)))
}


# Create a data frame
data <- expand.grid(Year = years, TrapID = trap_ids, Sample = 1:samples_per_season) %>%
  mutate(
    PairID = pair_ids[TrapID],
    BeetleCount = rpois(n(), lambda = 20),  # Random counts, change lambda as needed
    Temperature = rnorm(n(), mean = 15, sd = 5),  # Random temperature (°C)
    SPEI = rnorm(n(), mean = 0, sd = 1)  # Random SPEI values
  )

# View the first few rows of the dataset
head(data)



# Basic model:
# mixed effects model: year, temp, spei as fixed effects, trap id as random 
model <- lmer(BeetleCount ~ Year + Temperature + SPEI + (1 | TrapID), data = data)
summary(model)

# compare local vs global variability --------------------------------

# calculate standsrd deviation (variance) per trap
# Local variability
local_var <- aggregate(BeetleCount ~ PairID + Year, data, var)

# Global variability
global_var <- aggregate(BeetleCount ~ Year, data, var)

# Merge the two for comparison
combined_var <- merge(local_var, global_var, by = "Year")
names(combined_var) <- c("Year", "PairID", "LocalVariance", "GlobalVariance")

head(combined_var)


# stat test local vs global -----------------------------------------------

# You might need to reshape your data for this test
anova_test <- aov(LocalVariance ~ GlobalVariance, data = combined_var)
summary(anova_test)

# plot changing variance
ggplot(combined_var, aes(x = Year, y = LocalVariance)) +
  geom_line(aes(color = "Local Variance")) +
  geom_line(aes(y = GlobalVariance, color = "Global Variance")) +
  labs(title = "Local vs Global Variance over Years")
