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

# read climate table
# Get SPEI and clim data: they are for whole year! check only veg season? now 3-10 (includes march)
df_clim <- fread('outTable/xy_clim_DWD.csv')
df_spei_months <- fread('outTable/xy_spei_all_DWD.csv')


# Create a data frame
data <- 
  dat.ips.clean %>% 
  filter(month %in% 3:10) %>% 
  dplyr::select(fangmenge, falsto_name, monsto_name, year, month) %>% 
  group_by(falsto_name, monsto_name, year, month) %>% 
  summarise(sum_ips = sum(fangmenge, na.rm  = T)) %>% # sum up the numbers per months first
  #expand.grid(Year = years, TrapID = trap_ids, Sample = 1:samples_per_season) %>%
   
  left_join(df_clim, by = join_by(falsto_name, year, month)) %>% 
  left_join(dplyr::filter(df_spei_months, scale == 3), 
            by = join_by(falsto_name, year, month)) %>% 
  dplyr::rename(
    pairID = monsto_name,
    trapID = falsto_name#,
    #sum_ips = fangmenge
    ) %>%
  mutate(pairID = factor(pairID),
         trapID = factor(trapID)) %>% 
  # add spei lags
  ungroup(.) %>% 
  arrange(trapID, month, year) %>%
  group_by(trapID, month) %>%
  mutate(spei1 = lag(spei, order_by = year),
         spei2 = lag(spei, order_by = year, n = 2))




# View the first few rows of the dataset
head(data)



# Basic model:
#



library(ggeffects)

library(MASS)
library(car)     # for VIF
library(glmmTMB) # many families
#library(plyr)


library(lme4) #pseudo-R2 and and mixed model functionality
library(MuMIn) #dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc) #pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa) #model diagnostics
library(effects) #what do my marginal effects look like?
library(performance) #binomial model diagnostics
library(emmeans) #post hoc for categorical predictors

# test for autocorrelation
library(lmtest)

# colors
library(RColorBrewer)

library(knitr)   # for table outputs




# mixed effects model: year, temp, spei as fixed effects, trap id as random 
m.var1 <- glmmTMB(sum_ips ~ tmp + spei + (1 | pairID),
               family = nbinom2,
               data = data)

plot(m.var1)

output<- simulateResiduals(m.var2, plot = T)

m.var2 <- glmmTMB(sum_ips ~ tmp + spei + (1 | pairID) + (1 | year),
                  family = nbinom2,
                  data = data)

m.var3 <- glmmTMB(sum_ips ~ tmp + spei1 + (1 | pairID) + (1 | year),
                  family = nbinom2,
                  data = data)

m.var4 <- glmmTMB(sum_ips ~ tmp + spei2 + (1 | pairID) + (1 | year),
                  family = nbinom2,
                  data = data)


# include interaction in time
m.var5 <- glmmTMB(sum_ips ~ tmp * year + spei2 * year + (1 | pairID) + (1 | year),
                  family = nbinom2, data = data)

AIC(m.var4, m.var3, m.var2, m.var1)

summary(m.var4)
plot(allEffects(m.var4))

# compare local vs global variability --------------------------------

# calculate standsrd deviation (variance) per trap
# Local variability
local_var <- aggregate(sum_ips ~ pairID + year, data, var)

# Global variability
global_var <- aggregate(sum_ips ~ year, data, var)

# Merge the two for comparison
combined_var <- merge(local_var, global_var, by = "year")
names(combined_var) <- c("Year", "PairID", "LocalVariance", "GlobalVariance")

head(combined_var)


# stat test local vs global -----------------------------------------------
hist(combined_var$LocalVariance)
hist(combined_var$GlobalVariance)

# You might need to reshape your data for this test
kruskal_test <- kruskal.test(LocalVariance ~ GlobalVariance, data = combined_var)
summary(kruskal_test)

# plot changing variance
ggplot(combined_var, aes(x = Year, y = LocalVariance)) +
  # This will add error bars for mean +/- SD
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color = "blue") +
  # This will add a point for the mean of Local Variance
  stat_summary(fun = mean, geom = "point", color = "blue") +
  # This will add a line for the mean of Local Variance
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue") +
  # This will add a line for Global Variance
  geom_line(aes(y = GlobalVariance, color = "Global Variance"), color = "red") +
  # Optional: Color legend adjustments
  scale_color_manual(values = c("Global Variance" = "red", "Local Variance" = "blue")) +
  labs(title = "Local vs Global Variance over Years", color = "Variance Type") +
  theme_classic()
summary(combined_var)


# summarize values to explain variance by climatic factors
local_var <- data %>%
  group_by(pairID, year) %>%
  summarize(
    LocalVariance = var(sum_ips, na.rm = TRUE),  # Calculate variance of beetle counts
    MeanTemperature = mean(tmp, na.rm = TRUE),  # Calculate mean temperature
    MeanSPEI = mean(spei, na.rm = TRUE),  # Calculate mean SPEI
    .groups = 'drop'  # This option removes the grouping structure after summarizing
  ) %>% 
  arrange(pairID, year) %>%
  group_by(pairID) %>%
  mutate(spei1 = lag(MeanSPEI, order_by = year),
         spei2 = lag(MeanSPEI, order_by = year, n = 2))


pairs(LocalVariance ~ MeanTemperature + MeanSPEI + spei1 +spei2, data = local_var)

# predict local valriance based on climatic predcitors
m1 <- glmmTMB(LocalVariance ~ MeanTemperature +spei2,# + #(1 | trapID), 
              data = local_var, 
              family = tweedie)

summary(m1)
plot(allEffects(m1))
simulateResiduals(m1, plot = T)
