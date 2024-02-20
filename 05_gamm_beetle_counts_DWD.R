


rm(list=ls()) 
gc()

# Intro -----------------------------------------------------------------------
### Read my paths -----------------------------------------------------------
source('myPaths.R')



### get libs ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(ggplot2)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# Stats
library('here')
#library('mgcv')
#library('gratia')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
#library('gganimate')
library('cowplot')
library('tidyr')
library("knitr")
#library("viridis")
library('readr')
#library('itsadug')

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





##### read data ----------------------------------------------------------------

load(file=   "outData/final_table.Rdata")
load(file =  "outData/buffers.Rdata")  # df_RS_out
load(file =  "outData/lisa.Rdata")     # read LISA Moran's I stats



# prepare data RS --------------------------------------------
df_RS_out <- df_RS_out %>% 
  dplyr::select(-c(globalid, id)) %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>% 
  dplyr::rename(trapID = falsto_name)

# change column names
lisa_merged_df <- lisa_merged_df %>% 
  mutate(falsto_name = gsub("[^A-Za-z0-9_]", "", falsto_name)) %>% 
  dplyr::rename(trapID = falsto_name,
                sum_ips = sum_beetle) %>% 
  dplyr::select(c(year, trapID, sum_ips, Morans_I))

# add lag: previous year counts, previous year temperature ----------------------
dat_lag <-   dat_fin %>%
    ungroup(.) %>% 
   group_by(trapID) %>%
    dplyr::mutate(previous_sum_ips    =  dplyr::lag(sum_ips,    order_by = year),
                  previous_sum_ips2   =  dplyr::lag(sum_ips,    order_by = year, n = 2),
                  previous_agg1       =  dplyr::lag(agg_doy ,   order_by = year, n = 1), 
                  previous_agg2       =  dplyr::lag(agg_doy ,   order_by = year, n = 2), 
                  previous_peak_diff1 =  dplyr::lag(peak_diff , order_by = year, n = 1), 
                  previous_peak_diff2 =  dplyr::lag(peak_diff , order_by = year, n = 2),
                  previous_veg_tmp    =  dplyr::lag(veg_tmp,    order_by = year),
                  previous_spring_tmp =  dplyr::lag(spring_tmp, order_by = year),
                  previous_veg_prcp   =  dplyr::lag(veg_prcp,   order_by = year),
                  previous_spei1      =  dplyr::lag(spei1,      order_by = year),
                  previous_spei3      =  dplyr::lag(spei3,      order_by = year),
                  previous_spei6      =  dplyr::lag(spei6,     order_by = year),
                  previous_spei12     =  dplyr::lag(spei12,     order_by = year),
                  previous_spei1_2    =  dplyr::lag(spei1,      order_by = year, n = 2),
                  previous_spei3_2    =  dplyr::lag(spei3,      order_by = year, n = 2),
                  previous_spei6_2    =  dplyr::lag(spei6,     order_by = year, n = 2),
                  previous_spei12_2   =  dplyr::lag(spei12,     order_by = year, n = 2)) %>% 
  dplyr::mutate(population_growth     = (sum_ips - previous_sum_ips) / previous_sum_ips * 100,
                population_growth2    = dplyr::lag(population_growth, order_by = year)) %>%  # lag population growth by one more year
  left_join(df_RS_out, 
            by = join_by(trapID, year, sum_ips)) %>%
  left_join(lisa_merged_df, 
            by = join_by(trapID, year, sum_ips)) %>%
  dplyr::mutate(previous_Moran       =  dplyr::lag(Morans_I , order_by = year),
                previous_Moran2      =  dplyr::lag(Morans_I , order_by = year, n = 2),
                ) %>% 
  dplyr::mutate(trapID = as.factor(trapID)) %>%
  dplyr::filter(year %in% 2015:2021)



### Spearman - select SPEI --------------------------------------------------
keep_speis <- c(#'sum_ips', 
                "spei1"  ,            
               "spei3", 'spei6',  "spei12",  "spei24", 
               # "annual_spei1",
               # "annual_spei3",
               # 'annual_spei6',
               #  "annual_spei12" ,  
               #  "annual_spei24" ,
                "previous_spei1",
               "previous_spei1_2",
               "previous_spei3",
               "previous_spei3_2",
               "previous_spei6",
               "previous_spei6_2",
               "previous_spei12",
               "previous_spei12_2")
  


# Spearman between SPEIs nd their lags
df_spearman_spei <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c(keep_speis)))

spearman_correlation_matrix <- cor(df_spearman_spei, 
                                   method = "spearman", 
                                   use = "complete.obs")

# Print the correlation matrix
#print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables

# Convert the correlation matrix to a long format
long_cor <- spearman_correlation_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable1") %>%
  pivot_longer(cols = -Variable1, names_to = "Variable2", values_to = "Correlation")

# Filter out redundant pairs (upper triangle of the matrix)
long_cor <- long_cor %>%
  filter(!duplicated(paste(pmin(Variable1, Variable2), pmax(Variable1, Variable2))))

# Sort by the absolute values of correlations
sorted_cor <- long_cor %>%
  arrange(abs(Correlation))

# Select the lowest correlations
# You can adjust the number of rows you want to view
lowest_correlations <- head(sorted_cor, n = 10)

# View the result
lowest_correlations





##### IPS spearmans -----------  
df_spearman_spei <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('sum_ips', keep_speis)))


spearman_correlation_matrix <- cor(df_spearman_spei, 
                                   method = "spearman", 
                                   use = "complete.obs")


# Extract and sort correlations of sum_ips with SPEI variables
cor_with_sum_ips <- spearman_correlation_matrix["sum_ips", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_sum_ips, decreasing = TRUE)

# Print the sorted correlations
print(sorted_correlations)



# Assuming spearman_correlation_matrix is your correlation matrix
#formatted_table <- kable(spearman_correlation_matrix, format = "pipe", digits = 2)

# Print the formatted table
#cat(formatted_table)

# th3e best is SPEI3, veg season!! for beetle sums!



##### DOY aggregation Spearman -----------  
df_spearman_spei <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('agg_doy', keep_speis)))


spearman_correlation_matrix <- cor(df_spearman_spei, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables
cor_with_y <- spearman_correlation_matrix["agg_doy", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_y, decreasing = TRUE)

# Print the sorted correlations
print(sorted_correlations)

# th3e best is annula SPEI1 for beetle aggregation date



##### DOY peak Spearmans -----------  
df_spearman_spei <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('peak_doy', keep_speis)))


spearman_correlation_matrix <- cor(df_spearman_spei, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables
cor_with_y <- spearman_correlation_matrix["peak_doy", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_y, decreasing = TRUE)

# Print the sorted correlations
print(sorted_correlations)





##### DOY peak difference Spearmans -----------  
df_spearman_spei <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('peak_diff', keep_speis)))


spearman_correlation_matrix <- cor(df_spearman_spei, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables
cor_with_y <- spearman_correlation_matrix["peak_diff", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_y, decreasing = TRUE)

# Print the sorted correlations
print(sorted_correlations)



### Spearman - select spring vs veg temperature ------------------------------------------------------

keep_temps <- c( 
  "spring_tmp", "veg_tmp",
  "previous_veg_tmp",    "previous_spring_tmp",
  "previous_veg_prcp")


##### IPS separmans -----------  
df_spearman_temps <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('sum_ips', keep_temps)))


spearman_correlation_matrix <- cor(df_spearman_temps, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables
cor_with_sum_ips <- spearman_correlation_matrix["sum_ips", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_sum_ips, decreasing = TRUE)

# Print the sorted correlations
print(sorted_correlations)


# previous_veg_tmp             veg_tmp          spring_tmp previous_spring_tmp 
# 0.08339555          0.07393616          0.07088885          0.05713606 
# previous_veg_prcp 
# -0.04463953 




##### DOY aggregation Spearman -----------  
df_spearman_temps <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('agg_doy', keep_temps)))


spearman_correlation_matrix <- cor(df_spearman_temps, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations of sum_ips with SPEI variables
cor_with_sum_ips <- spearman_correlation_matrix["agg_doy", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_sum_ips, decreasing = F)

# Print the sorted correlations
print(sorted_correlations)

# best: spring temp, veg_temp



##### DOY peak Spearmans -----------  
df_spearman_temps <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('peak_doy', keep_temps)))


spearman_correlation_matrix <- cor(df_spearman_temps, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations
cor_with_sum_ips <- spearman_correlation_matrix["peak_doy", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_sum_ips, decreasing = F)

# Print the sorted correlations
print(sorted_correlations)

print(sorted_correlations)
# spring_tmp             veg_tmp previous_spring_tmp    previous_veg_tmp 
# -0.21839001         -0.20928259         -0.15303943         -0.13813963 
# previous_veg_prcp 
# -0.06085891 

##### DOY peak difference Spearmans -----------  
df_spearman_temps <- dat_lag %>% 
  ungroup(.) %>% 
  dplyr::select(all_of(c('peak_diff', keep_temps)))


spearman_correlation_matrix <- cor(df_spearman_temps, method = "spearman", use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation_matrix)

# Extract and sort correlations
cor_with_sum_ips <- spearman_correlation_matrix["peak_diff", -1]  # Exclude the first column (self-correlation)
sorted_correlations <- sort(cor_with_sum_ips, decreasing = F)

# Print the sorted correlations
print(sorted_correlations)

# previous_spring_tmp    previous_veg_tmp   previous_veg_prcp             veg_tmp 
# -0.01454686         -0.01821267         -0.06324645         -0.09972453 
# spring_tmp 
# -0.11110834 






# Summary table  ----------------------------------------------------------

# Represent results using quantiles, as they are skewed?
qntils = c(0, 0.25, 0.5, 0.75, 1)

qs_dat_fin <- 
  dat_fin %>% 
    ungroup(.) %>% 
    filter(year %in% 2015:2021) %>% 
  dplyr::reframe(sum_ips   = quantile(sum_ips, qntils, na.rm = T ),
            peak_doy  = quantile(peak_doy, qntils, na.rm = T ),
            agg_doy   = quantile(agg_doy, qntils, na.rm = T ),
            peak_diff = quantile(peak_diff, qntils, na.rm = T )) %>% 
  t() %>%
  round(1) %>% 
  as.data.frame()

(qs_dat_fin)  


means_dat_fin <- 
  dat_fin %>% 
  ungroup(.) %>% 
  filter(year %in% 2015:2021) %>% 
  dplyr::reframe(sum_ips   = mean(sum_ips, na.rm = T),
                 peak_doy  = mean(peak_doy, na.rm = T),
                 agg_doy   = mean(agg_doy, na.rm = T),
                 peak_diff = mean(peak_diff, na.rm = T)) %>% 
  t() %>%
  round(1) %>% 
  as.data.frame() 


str(means_dat_fin)

# merge qs and mean tables
summary_out <- cbind(qs_dat_fin, means_dat_fin)# %>%

# Export as a nice table in word:
sjPlot::tab_df(summary_out,
               col.header = c(as.character(qntils), 'mean'),
               show.rownames = TRUE,
               file="outTable/summary_out.doc",
               digits = 1) 



### summary table per year: -------------------------------------------------------


means_dat_fin_year <- 
  dat_fin %>% 
  ungroup(.) %>% 
  filter(year %in% 2015:2021) %>% 
  group_by(year) %>% 
  dplyr::reframe(mean_sum_ips   = mean(sum_ips, na.rm = T),
                 mean_agg_doy   = mean(agg_doy, na.rm = T),
                  mean_peak_doy  = mean(peak_doy, na.rm = T),
                 mean_peak_diff = mean(peak_diff, na.rm = T),
                 sd_sum_ips   = sd(sum_ips, na.rm = T),
                 sd_agg_doy   = sd(agg_doy, na.rm = T),
                  sd_peak_doy  = sd(peak_doy, na.rm = T),
                 sd_peak_diff = sd(peak_diff, na.rm = T)) %>% 
  mutate(Population_level    = stringr::str_glue("{round(mean_sum_ips,1)}±{round(sd_sum_ips,1)}"),
         Colonization_DOY    = stringr::str_glue("{round(mean_agg_doy,1)}±{round(sd_agg_doy,1)}"),
         Peak_growth_DOY    = stringr::str_glue("{round(mean_peak_doy,1)}±{round(sd_peak_doy,1)}"),
         Peak_growth    = stringr::str_glue("{round(mean_peak_diff,1)}±{round(sd_peak_diff,1)}")) %>% 
  dplyr::select(year, Population_level, Colonization_DOY, Peak_growth_DOY,Peak_growth) 


(means_dat_fin_year)


# Export as a nice table in word:
sjPlot::tab_df(means_dat_fin_year,
#               col.header = c(as.character(qntils), 'mean'),
               show.rownames = TRUE,
               file="outTable/summary_out_year.doc",
               digits = 0) 


# Spagetti plots: development over time ----------------------------------------



plot_data_with_average <- function(data, y_var, y_label, my_title) {
  # Convert the y_var argument to a symbol to use in aes()
  y_var_sym <- rlang::sym(y_var)
  
  data %>%
    ungroup() %>%
    filter(year %in% 2015:2021) %>%
    ggplot(aes(x = year, y = !!y_var_sym, group = trapID)) +
    labs(x = 'Year', y = y_label, title = my_title) +
    geom_line(alpha = 0.1) +  
    stat_summary(
      aes(x = year, y = !!y_var_sym, group = 1), 
      fun = mean,  # Calculate the mean for each year
      geom = "line", 
      color = "red",  # Ensure the average line is also red
      size = 1  # Make the average line slightly thicker than individual lines
    ) +
    # stat_summary(
    #   aes(x = year, y = !!y_var_sym, group = 1),
    #   fun.data = mean_sdl,  # Calculate mean and standard deviation
    #   fun.args = list(mult = 1),  # 'mult = 1' for 1 SD, 'mult = 2' for 2 SDs, etc.
    #   geom = "errorbar",  # Use error bars to represent the SD
    #   width = 0.1,
    #   col = 'red'  # SD indicators in red
    # ) +
    
    theme_minimal(base_size = 10) +
    theme(
      aspect.ratio = 1, 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "white", colour = "black")
    )
}

# Example usage:
p_spagett_ips       <- plot_data_with_average(dat_fin, "sum_ips", "Beetle counts [#]", '[a]')
p_spagett_agg_doy   <- plot_data_with_average(dat_fin, "agg_doy", "Colonization [DOY]", '[b]')
p_spagett_peak_doy  <- plot_data_with_average(dat_fin, "peak_doy", "Peak population Growth [DOY]", '[c]')
p_spagett_peak_diff <- plot_data_with_average(dat_fin, "peak_diff", "Peak population Growth [#]", '[d]')

ggarrange(p_spagett_ips, 
          p_spagett_agg_doy, 
          p_spagett_peak_doy, 
          p_spagett_peak_diff, ncol = 2,nrow = 2, align = 'hv')


## Scale predictors ================================================================================

# skip columns if not for scaling
#spei_cols <- grepl("spei", names(dat_lag))

additional_skip_cols <- c('trapID', 'pairID', "x", "y", 
                          'year', 'agg_doy', 'peak_doy', 'peak_diff',
                          'sum_ips', 'previous_sum_ips', 'previous_sum_ips2', 
                          'wind_beetle', 'Morans_I')

# Combine spei columns with additional columns to skip
skip_col <- additional_skip_cols #c(names(dat_lag)[spei_cols], additional_skip_cols)

# scale also SPEIs

# export new df
dat_lag_scaled <-
  dat_lag %>%
  ungroup(.) %>% 
  mutate(sc_sum_ips = sum_ips,  # add beetle counts as scaled values
         sc_previous_sum_ips = previous_sum_ips,
         sc_previous_sum_ips2 = previous_sum_ips2) %>% 
  dplyr::select(-all_of(skip_col )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_lag %>% dplyr::select(all_of(skip_col)), .)

#simplify, keep all columns as vectors 
dat_lag_scaled <- dat_lag_scaled %>% 
  ungroup(.) %>% 
  dplyr::mutate(across(where(~is.numeric(.) && is.matrix(.)), as.vector))


# remove additionsl NAs
dat_lag_scaled_complete <- dat_lag_scaled %>%
  dplyr::select(-anom_wind_beetle) %>% 
  na.omit()
#  filter(across(everything(), ~ !is.na(.)))



#### get a table for RS data --------------------

skip_cols_RS <- c('trapID', 'pairID', 'year', 'wind_beetle')



dat_lag_RS <- dat_lag %>% 
  dplyr::select('trapID', 'pairID', 'year', 'wind_beetle',
    'agg_doy', 'peak_doy', 'peak_diff', 
                'sum_ips', 'previous_sum_ips', 'previous_sum_ips2',
                "previous_agg1" ,      "previous_agg2",
    "previous_peak_diff1", "previous_peak_diff2", 
     "veg_tmp", "previous_veg_tmp",  
                "spei3",  "previous_spei3", "previous_spei3_2",
                "population_growth",   "population_growth2",
    "Morans_I") %>% 
  mutate(beetle_sum2yrs = previous_sum_ips + previous_sum_ips2)  # new variable having cumulative values ob beetle population per two years
  

dd_RS_scaled <- dat_lag_RS %>% 
  ungroup(.) %>% 
  dplyr::select(-all_of(skip_cols_RS )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_lag_RS %>% dplyr::select(all_of(skip_cols_RS)), .)




# Analyses =====================================================================

# beetle population sum vs climate drivers
# DOY aggregation vs climate drivers
# DOY peak vs climate drivers
# peak difference vs climate drivers




#### spatial synchronization with veg_tmp? ---------------------------------------

hist(dat_lag_scaled$Morans_I)

pairs(Morans_I ~  previous_veg_tmp + previous_sum_ips + spei1, dat_lag_scaled)

p.temp <- dat_lag %>% 
  ggplot(aes(x = veg_tmp,
             y = Morans_I)) +
  geom_smooth() +
  theme_bw()


p.previous.temp <- dat_lag %>% 
  ggplot(aes(x = previous_veg_tmp,
             y = Morans_I)) +
  geom_smooth()+
  theme_bw()


dat_lag %>% 
  ggplot(aes(x = veg_prcp,
             y = Morans_I)) +
  geom_point()


p.ips_sum <- dat_lag %>% 
  ggplot(aes(x = sum_ips/10000,
             y = Morans_I)) +
  geom_smooth()+
  theme_bw()

dat_lag %>% 
  ggplot(aes(x = year,
             y = Morans_I)) +
  geom_point()

p.spei <- dat_lag %>% 
  ggplot(aes(x = spei3,
             y = Morans_I)) +
  geom_smooth() +
  theme_bw()


ggarrange(p.temp, p.previous.temp, p.spei, p.ips_sum)

# Fit a linear regression model to explain Moran's I ------------------------------------
m_moran1 <- lm(Morans_I ~ poly(veg_tmp,2) + previous_sum_ips + poly(spei3,2), data = dat_lag_scaled_complete)

plot(allEffects(m_moran1))

summary(m_moran1)
r2(m_moran1)

# lmer: account for the random effect of teh trap: the effect of the trap can vary between years
m_moran2 <- lmer(Morans_I ~ veg_tmp + previous_spei3_2 + (1|pairID), data = dat_lag_scaled)

# remove prec as it is correlated with veg_tmp
m_moran3 <- lmer(Morans_I ~ veg_tmp + (1|pairID), data = dat_lag_scaled)

#add pairID as random
m_moran4 <- lmer(Morans_I ~ veg_tmp + spei3 + (1|pairID), data = dat_lag_scaled)

# if temp and prec are correlated, maybe there is an interaction effect?
m_moran5 <- lmer(Morans_I ~ veg_tmp + previous_spei3+ (1 | pairID), data = dat_lag_scaled)


# try different family: tweedie, gaussian and linear model does not work well

m_moran6 <- glmmTMB(Morans_I ~ veg_tmp + previous_spei3_2,# + #(1 | trapID), 
                 data = dat_lag_scaled, 
                 family = tweedie)

# add randm effect
m_moran7 <- glmmTMB(Morans_I ~ previous_spei3_2 + (1 | pairID), 
                    data = dat_lag_scaled, 
                    family = tweedie)

m_moran8 <- glmmTMB(Morans_I ~ veg_tmp + previous_spei3_2 + (1 | pairID), 
                    data = dat_lag_scaled, 
                    family = tweedie)

m_moran9 <- glmmTMB(Morans_I ~ veg_tmp*previous_spei3_2 + (1 | pairID), 
                     data = dat_lag_scaled, 
                     family = tweedie)

m_moran10 <- glmmTMB(Morans_I ~ previous_spei3_2 + (1 | pairID), 
                    data = dat_lag_scaled, 
                    family = tweedie)
# use only gaussian family, as tweedie is suitable for 0s and positive values, not negative
m_moran11 <- glmmTMB(Morans_I ~ previous_spei3_2 + (1 | year), 
                     data = dat_lag_scaled, 
                     family = tweedie)

m_moran12 <- lm(Morans_I ~ previous_spei3_2 + veg_tmp, 
                     data = dat_lag_scaled)

m_moran13 <- lmer(Morans_I ~ previous_spei3_2 + veg_tmp + (1 | pairID), 
                data = dat_lag_scaled)

m_moran14 <- lmer(Morans_I ~ previous_spei3_2*veg_tmp + (1 | pairID), 
                  data = dat_lag_scaled)

m_moran15 <- lmer(Morans_I ~ spei3  + veg_tmp + sc_previous_sum_ips + (1 | pairID), 
                  data = dat_lag_scaled)

m_moran16 <- lmer(Morans_I ~ spei3  + veg_tmp + sc_previous_sum_ips + (1 | pairID), 
                  data = dat_lag_scaled)

AIC(m_moran2, m_moran3,m_moran4,m_moran5,m_moran6,m_moran7,m_moran8, 
    m_moran9, m_moran10, m_moran11, m_moran12,m_moran13,m_moran14)

output <- simulateResiduals(m_moran15, plot = T)

acf(residuals(m_moran15))
dw_result <- dwtest(residuals(m_moran15) ~ fitted(m_moran15))
(dw_result)

# check for autocorrelation


# choose teh right family:
summary(dat_lag_scaled$Morans_I)
hist(dat_lag_scaled$Morans_I)


library(MASS)

# Assuming the minimum value is negative, shift all data to be positive
shifted_value = abs(min(dat_lag_scaled$Morans_I, na.rm = TRUE)) + 1
dat_lag_scaled$Morans_I_shifted = dat_lag_scaled$Morans_I + shifted_value

# Apply Box-Cox Transformation
bc_trans = boxcox(dat_lag_scaled$Morans_I_shifted ~ 1, lambda = seq(-2, 2, by = 0.1))

# Find the lambda that maximizes the log-likelihood
lambda_opt = bc_trans$x[which.max(bc_trans$y)]
dat_lag_scaled$Morans_I_transformed = (dat_lag_scaled$Morans_I_shifted^lambda_opt - 1) / lambda_opt

model = lmer(Morans_I_transformed ~ veg_tmp + previous_spei3_2 + (1 | pairID), data = dat_lag_scaled)
summary(model)

par(mfrow = c(2, 2))
plot(model)

# add previous year
m_moran9 <- glmmTMB(Morans_I ~ previous_veg_tmp + previous_spei3_2 + (1 | year), 
                    data = dat_lag_scaled, 
                    family = tweedie)


# remove random effects 
m_moran10 <- glmmTMB(Morans_I ~ previous_veg_tmp + veg_prcp, #+ (1 | year), 
                    data = dat_lag_scaled, 
                    family = tweedie)


# remove random effects 
m_moran11 <- glmmTMB(Morans_I ~ previous_veg_tmp, # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled, 
                     family = tweedie)


# gaussian with poly link
m_moran12 <- glmmTMB(Morans_I ~ poly(previous_veg_tmp,2), # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled_complete)


# gaussian with poly link
m_moran13 <- glmmTMB(Morans_I ~ poly(previous_veg_tmp,3), # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled_complete)


# remove poly 
m_moran14 <- glmmTMB(Morans_I ~ previous_veg_tmp, # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled_complete)


# null model 
m_moran15 <- glmmTMB(Morans_I ~ 1, # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled_complete)

# add current year beetle counts 
m_moran16 <- glmmTMB(Morans_I ~ previous_veg_tmp + sum_ips, # + veg_prcp, #+ (1 | year), 
                     data = dat_lag_scaled_complete)


# add previous year beetle counts 
m_moran17 <- glmmTMB(Morans_I ~ previous_veg_tmp + sum_ips + 
                       previous_sum_ips + (1 | year), 
                     data = dat_lag_scaled_complete,
                     family = tweedie)


# add lagged Morans vals 
m_moran18 <- glmmTMB(Morans_I ~ previous_veg_tmp + sum_ips + 
                       previous_sum_ips + previous_Moran, # + (1 | year), 
                     data = dat_lag_scaled_complete,
                     family = tweedie)


# add lagged Morans vals 
m_moran19 <- glmmTMB(Morans_I ~ previous_veg_tmp + sum_ips_scaled + 
                       previous_sum_ips + previous_Moran + previous_Moran2, # + (1 | year), 
                     data = dat_lag_scaled_complete,
                     family = tweedie)


# add lagged Morans vals 
m_moran20 <- glmmTMB(Morans_I ~ previous_veg_tmp + #sum_ips_scaled + 
                       #previous_sum_ips + 
                       spei3, #+
                       #previous_Moran + previous_Moran2, # + (1 | year), 
                     data = dat_lag_scaled_complete,
                     family = tweedie)


# add lagged Morans vals 
m_moran21 <- glmmTMB(Morans_I ~ previous_veg_tmp + #sum_ips_scaled + 
                       #previous_sum_ips + 
                       spei3, #+
                     #previous_Moran + previous_Moran2, # + (1 | year), 
                     data = dat_lag_scaled_complete)


# add lagged Morans vals 
m_moran22 <- glm(Morans_I ~ previous_veg_tmp + #sum_ips_scaled + 
                       #previous_sum_ips + 
                       spei3, #+
                     #previous_Moran + previous_Moran2, # + (1 | year), 
                     data = dat_lag_scaled_complete)


m_moran22 <- glm(Morans_I ~ spei3, 
                 #previous_Moran + previous_Moran2, # + (1 | year), 
                 data = dat_lag_scaled_complete)





plot(allEffects(m_moran22))

r2(m_moran22)
cor(dat_lag_scaled_complete$previous_veg_tmp,dat_lag_scaled_complete$spei3 )
cor(dat_lag_scaled_complete$previous_veg_prcp,dat_lag_scaled_complete$spei3 )


write.csv(dat_lag_scaled_complete, 'test.csv')



# Linear Regression Model
model <- lm(Morans_I ~ veg_tmp + veg_prcp + spei3, data = dat_lag_scaled_complete)
summary(model)

# Diagnostic Plots
par(mfrow = c(2, 2))
plot(m_moran20)



AIC(m_moran13, m_moran12,m_moran14,m_moran15, m_moran16 , m_moran18, m_moran19) # for tweedie
AICc(m_moran6, m_moran7, m_moran8, m_moran9, m_moran11) # for gaussian
summary(m_moran21)
plot(allEffects(m_moran21))
r2(m_moran21)

resid_outpu <- simulateResiduals(m_moran21, plot = T)

### RS:  Link beetle data with observed RS damage ------------------------------------------------------------------------


###### get glm RS vs ips sums: DREDGE  ----------------------------------------------------------------------
# investigate all variables & dredge

# keep only complete cases - can alter this to keep more rows? can remove unnecessary columns
# remove unnecessare predictors, also lag values only by 1 year
#col_remove <- c( "x",  "y", "sm" ,"vpd", "sm_z","vpd_z", "previous_sum_ips2","previous_agg2", "population_growth2" )

# remove extra columns !!!! - not necessary to remove those columns, as those are complete!
#dat_lag_scaled <- dat_lag_scaled %>%
#  dplyr::select(-all_of(col_remove)) #%>% 
 # filter_all(all_vars(!is.na(.))) %>% 
#  filter(across(everything(), ~ !is.na(.)))

###### get RS with beetle population predictors: ----------------------------------
# add predictors one by one


dat_lag_RS_compl <- dd_RS_scaled  %>% 
  na.omit()

hist(dat_lag_RS_compl$wind_beetle)

simple_model <- glm.nb(wind_beetle ~ #sum_ips + 
                         beetle_sum2yrs +
                         previous_veg_tmp +
                         #previous_sum_ips +
                         #previous_sum_ips2 +
                         #veg_tmp + 
                         #spring_tmp + 
                         
                         #spei3 + 
                         #previous_spei3 +
                         previous_spei3_2 +
                         #veg_prcp + 
  
                         #agg_doy + # correlated with sum_ips
                         #peak_doy +  ## correlated with sum_ips
                         # peak_diff + # correlated with sum_ips and its lags
                         
                         #previous_agg1 + 
                         #previous_agg2 +
                         #previous_peak_diff1 + 
                         #previous_peak_diff2 +
                         #previous_Moran +
                         #Morans_I +    # correlated with sum_ips and its lags
                         population_growth2,
                       na.action = 'na.fail',
                       data = dat_lag_RS_compl )


#  identify the most influential predictors

# Subsetting the dataframe to include only relevant predictors
predictors_df <- dat_lag_RS_compl[, c("wind_beetle", "agg_doy", "peak_doy", "peak_diff", 
                                      "sum_ips", "previous_sum_ips", "previous_sum_ips2", 
                                      "previous_agg1", "previous_agg2", "previous_peak_diff1", 
                                      "previous_peak_diff2", "veg_tmp", "previous_veg_tmp", 
                                      "spei3", "previous_spei3", "previous_spei3_2", 
                                      "population_growth", "population_growth2", "Morans_I")]

# Calculate the correlation matrix
correlation_matrix <- cor(predictors_df, use = "complete.obs") # 'complete.obs' handles missing values by case-wise deletion

# View the correlation matrix
print(correlation_matrix)

# identifiued teh most influentials predictors, while minimising collinearity between predictors



# Define the global model
simple_model <- glm.nb(wind_beetle ~ previous_sum_ips + 
                         beetle_sum2yrs +
                         previous_veg_tmp +
                         previous_spei3_2 +
                         population_growth2,
                       na.action = 'na.fail',
                       data = dat_lag_RS_compl)

# Run dredge to get all possible model combinations
model_candidates <- dredge(simple_model)

# Extract the top 10 models
top_models <- head(model_candidates, 10)

# Print the top 10 models
print(top_models)


mRS1<-  glm.nb(wind_beetle ~ #previous_sum_ips + 
                 beetle_sum2yrs,# +
                 #previous_veg_tmp +
                 #previous_spei3_2 +
                 #population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)
# mRS1 has perfect residuals!!

mRS2<-  glm.nb(wind_beetle ~ previous_sum_ips,# + 
                # beetle_sum2yrs,# +
               #previous_veg_tmp +
               #previous_spei3_2 +
               #population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)



mRS3<-  glm.nb(wind_beetle ~ #previous_sum_ips# + 
               # beetle_sum2yrs,# +
               #previous_veg_tmp +
               #previous_spei3_2 +
               population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)

mRS4<-  glm.nb(wind_beetle ~ #previous_sum_ips# + 
                  beetle_sum2yrs +
                 previous_veg_tmp +
                 previous_spei3_2 +
                 population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)

mRS5<-  glm.nb(wind_beetle ~ previous_sum_ips + 
                # beetle_sum2yrs +
                 previous_veg_tmp +
                 previous_spei3_2 +
                 population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)

fin.m.RS <-mRS5 
r2(mRS5)
summary(mRS5)
simulateResiduals(mRS5, plot = T)

AIC(mRS2, mRS1,mRS3, mRS4, mRS5)


anyNA(dat_lag_RS )

# run approach similar to DREDGE 

library(pscl)
library(MuMIn)

# Define a set of predictors
predictors <- c("sum_ips", 
                "previous_sum_ips",
               # "previous_sum_ips2",
                "previous_veg_tmp", 
               # "veg_tmp", 
                "spei3",
               "previous_spei3",   
               "previous_spei3_2", 
               
              # "agg_doy",
               "previous_agg1",
               "peak_doy", 
                "population_growth2")

# Create a dataframe to store model results

model_results <- data.frame(model = character(), AIC = numeric(), stringsAsFactors = FALSE)

# Loop over all possible predictor combinations
for (i in 1:length(predictors)) {
  combos <- combn(predictors, i, simplify = FALSE)
  for (combo in combos) {
    formula <- as.formula(paste("wind_beetle ~", paste(combo, collapse = "+")))
    #model <- hurdle(formula, data = dat_lag_RS_compl, dist = "negbin", zero.dist = "binomial")
    model <- glm.nb(formula, data = dat_lag_RS_compl,na.action = 'na.fail')
    model_results <- rbind(model_results, data.frame(model = deparse(formula), AIC = AIC(model)))
  }
}
# Sort the models by AIC and select the top 10
top_models <- model_results[order(model_results$AIC), ][1:10, ]
print(top_models)



# get a vector of vif values
(vif_model_ips <- car::vif(simple_model))


# Run dredge on the global model
model_set <- dredge(simple_model)

# Sort by AIC
model_set_sorted <- model_set[order(model_set$AIC),]

# View the top models with the lowest AIC
model_set_sorted[1:30]

# the most important: 


# Model selection table 
# (Int)       ppl_gr2   prv_ag1 prv_Mrn prv_pek_df1 prv_sum_ips prv_sum_ip2 df   logLik   AICc delta weight
# 17 1.855                                               0.3648              3 -757.416 1520.9  0.00  0.124
# 18 1.901 -0.0009242                                    0.4377              4 -756.940 1522.0  1.08  0.072
# 1  1.940                                                                   2 -759.002 1522.0  1.15  0.070
# 25 1.839                                  -0.3625      0.6019              4 -757.027 1522.1  1.25  0.066
# 49 2.140                                               0.4307  -1.516e-05  4 -757.155 1522.4  1.51  0.058
# 50 2.374 -0.0009695                                    0.5638  -2.558e-05  5 -756.273 1522.7  1.79  0.051
# 19 1.850            -0.050400                          0.3435              4 -757.401 1522.9  2.00  0.045

###### Manually create models -------------------------------------------------------

m1 <- glm.nb(wind_beetle ~ previous_sum_ips,
                       na.action = 'na.fail',
                       data = dat_lag_scaled_complete)

m2 <- glm.nb(wind_beetle ~ previous_sum_ips +
               population_growth2,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)

m3 <- glm.nb(wind_beetle ~ 1,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)


m4 <- glm.nb(wind_beetle ~ previous_sum_ips +
               #population_growth2 +
               previous_peak_diff1,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)


m5 <- glm.nb(wind_beetle ~ previous_sum_ips +
               previous_sum_ips2 ,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)

m6 <- glm.nb(wind_beetle ~ previous_sum_ips +
               population_growth2 +
               previous_sum_ips2,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)

m7 <- glm.nb(wind_beetle ~ previous_sum_ips +
               previous_agg1 ,
             na.action = 'na.fail',
             data = dat_lag_scaled_complete)

simulationOutput <- simulateResiduals(fittedModel = m0, plot = T) # residuals are bnetter for m_beetle6
AIC(m1, m2, m3, m4, m5, m6, m7)

r2(m1)
r2(m2)
r2(m3)
r2(m4)
r2(m5)
r2(m6)
r2(m7)

plot(allEffects(m6))






######  RS harvest =================================================================
# add predictors one by one
simple_model <- glm.nb(harvest ~ previous_sum_ips +
                         previous_sum_ips2 +
                         # agg_doy +
                         previous_agg1 + 
                         #previous_agg2 +
                         previous_peak_diff1 + 
                         # previous_peak_diff2 +
                         previous_Moran +
                         # Morans_I +
                         population_growth2,
                       na.action = 'na.fail',
                       data = dat_lag_scaled_complete)

# get a vector of vif values
(vif_model_ips <- car::vif(simple_model))


# Run dredge on the global model
model_set <- dredge(simple_model)

# Sort by AIC
model_set_sorted <- model_set[order(model_set$AIC),]

# View the top models with the lowest AIC
model_set_sorted[1:30]


simulationOutput <- simulateResiduals(fittedModel = simple_model  , plot = T) # m3.7.1.aut2

r2(simple_model)

m1 <-  glm.nb(harvest ~ #previous_sum_ips +
                previous_sum_ips2 +
                previous_Moran +
                # Morans_I +
                population_growth2,
              na.action = 'na.fail',
              data = dat_lag_scaled_complete)

summary(m1)
plot(allEffects(m1))



######  RS explore both climatic and beetles  variables -----------------------------------------------------



RS_global_model1 <- glm.nb(wind_beetle ~ sum_ips_scaled     + 
                           previous_sum_ips + 
                           previous_sum_ips2 + 
                             agg_doy +
                             previous_agg1 + 
                             previous_agg2 +
                               peak_doy +
                             peak_diff +
                             previous_peak_diff1 +
                             previous_peak_diff2 +
                            Morans_I+
                             previous_Moran +
                             previous_Moran2 +
                             population_growth +
                           population_growth2, #+ (1 | pairID/trapID), 
                      data = dat_lag_scaled_complete,
                      na.action = 'na.fail',
                      start = start_values#,
                     # family = negative.binomial(2.8588)
                      #, na.action = "na.omit"  "na.fail"
)

# update model by removing unnecessary variables
RS_global_model2 <- glm.nb(wind_beetle ~ #sum_ips_scaled     + 
                            previous_sum_ips + 
                            previous_sum_ips2 + 
                            #peak_diff +
                            previous_agg1 + 
                            #previous_agg2 +
                            veg_prcp +
                            veg_tmp +
                            spei3 +
                            #previous_veg_tmp +
                            previous_spei3 +
                            #previous_veg_prcp +
                         
                            #remained_spruce +
                            population_growth2, #+ (1 | pairID/trapID), 
                          data = dat_lag_scaled_complete,
                          na.action = 'na.fail'#,
                          # family = negative.binomial(2.8588)
                          #, na.action = "na.omit"  "na.fail"
)

# further narrow model down to only beetle pop and environm variables
RS_global_model3 <- glm.nb(wind_beetle ~ #sum_ips_scaled     + 
                             previous_sum_ips + 
                             #previous_sum_ips2 + 
                             #peak_diff +
                            # previous_agg1 + 
                             #previous_agg2 +
                             veg_prcp +
                             veg_tmp +
                             spei3 +
                             #previous_veg_tmp +
                             previous_spei3,# +
                             #previous_veg_prcp +
                             
                             #remained_spruce +
                           #  population_growth2
                           #, #+ (1 | pairID/trapID), 
                           data = dat_lag_scaled_complete,
                           na.action = 'na.fail'#,
                           # family = negative.binomial(2.8588)
                           #, na.action = "na.omit"  "na.fail"
)




###### ## examine best models from dredge -----------------------------------------------

# Assuming you want to examine model 207 from your model set
i = 2
selected_model <- get.models(model_set, subset = i)  # subset = 1 for the best model; adjust if needed

# Get a summary of the selected model
summary(selected_model[[1]])  # selected_model is a list, [[1]] accesses the first element

plot(allEffects(selected_model[[1]]))

simulationOutput <- simulateResiduals(fittedModel = selected_model[[1]]  , plot = T) # m3.7.1.aut2

testOutliers(selected_model , type = 'bootstrap') # i have outliers dectected, need to simplify the model! 
# Calculate pseudo-R^2 for the selected model
r.squared <- r.squaredGLMM(selected_model[[1]])
print(r.squared)  # prints both marginal and conditional R^2; for GLM, marginal R^2 is more relevant


## Manual models: models by hand: RS_beetle, RS_harvest, RS_dist ---------------------------------------
# always inclkude lagged variables
# dependent variable: sum of pixels: integer, I can use neg bin distribution
m_beetle1 <- glm.nb(wind_beetle ~ sum_ips, dat_lag_scaled)

# select: veg_prcp, temp, previous_sum beetles, spei3

# add previous beetle counts 
m_beetle2 <- glm.nb(wind_beetle ~ sum_ips + previous_sum_ips, dat_lag_scaled)
m_beetle3 <- glm.nb(wind_beetle ~  previous_sum_ips, dat_lag_scaled)

# add precip
m_beetle4 <- glm.nb(wind_beetle ~  previous_sum_ips + veg_prcp, dat_lag_scaled)

# add poly term
#m_beetle5 <- glm.nb(wind_beetle ~  poly(previous_sum_ips,2) + veg_prcp, dat_lag_scaled) # works only woiith dat_lag_scaled_completed

# poly(prev_ips) is not better, remove add temp
m_beetle6 <- glm.nb(wind_beetle ~  previous_sum_ips + veg_prcp + veg_tmp, dat_lag_scaled)  # the vbest!!!

# add sum_ips = (current year)
m_beetle6.2 <- glm.nb(wind_beetle ~  sum_ips_scaled + previous_sum_ips + veg_prcp + veg_tmp, dat_lag_scaled) 

# add previous veg_tmp
m_beetle6.3 <- glm.nb(wind_beetle ~  sum_ips_scaled + previous_sum_ips + veg_prcp + veg_tmp +
                        previous_veg_tmp, dat_lag_scaled)  

# remove surrent year sum beetles
m_beetle6.4 <- glm.nb(wind_beetle ~   previous_sum_ips + veg_prcp + veg_tmp +
                        previous_veg_tmp +spei3, dat_lag_scaled)  




# add poly(temp,2)
m_beetle7 <- glm.nb(wind_beetle ~  previous_sum_ips + veg_prcp + poly(veg_tmp), dat_lag_scaled)

# not better with poly(veg_tmp), add spei3
m_beetle8 <- glm.nb(wind_beetle ~  previous_sum_ips + veg_prcp + poly(veg_tmp,2) + poly(spei3,2), dat_lag_scaled_complete)


# add interaction
m_beetle9 <- glm.nb(wind_beetle ~  previous_sum_ips + veg_prcp + poly(veg_tmp,2) + poly(spei3,2) +
                      previous_sum_ips:veg_tmp, dat_lag_scaled_complete)


AIC(m_beetle8, m_beetle9)


# !!!!!
summary(m_beetle9)
r2(m_beetle9) 
# # R2 for Generalized Linear Regression - pseudo-R2 - not directly the % of variance explained but close
# Nagelkerke's R2: 0.16


# test interaction: beetles vs temp
m_beetle_interaction1 <- glm.nb(wind_beetle ~ previous_sum_ips + veg_prcp + veg_tmp + previous_sum_ips:veg_tmp, 
                               data = dat_lag_scaled)

# great looking residuals!! in m_beetle_interaction1.1
m_beetle_interaction1.1 <- glm.nb(wind_beetle ~ previous_sum_ips + 
                                    veg_prcp + veg_tmp + spei3 + 
                                    previous_sum_ips:veg_tmp, 
                                data = dat_lag_scaled_complete)


m_beetle_interaction2 <- glm.nb(wind_beetle ~ previous_sum_ips + veg_prcp + veg_tmp + previous_sum_ips:veg_tmp +
                                  previous_sum_ips:veg_prcp, 
                                data = dat_lag_scaled_complete)

m_beetle_interaction3 <- glm.nb(wind_beetle ~ previous_sum_ips + veg_prcp + poly(veg_tmp,2) + previous_sum_ips:veg_tmp, 
                                data = dat_lag_scaled)

m_beetle_interaction4 <- glm.nb(wind_beetle ~ previous_sum_ips + veg_prcp + poly(veg_tmp,2) + previous_sum_ips:poly(veg_tmp,2), 
                                data = dat_lag_scaled)


# !!!add interaction with previous temperatures
m_beetle_interaction5 <- glm.nb(wind_beetle ~ previous_sum_ips + previous_veg_tmp + previous_sum_ips:previous_veg_tmp, 
                                data = dat_lag_scaled)

# use 
m_beetle_interaction6 <- glm.nb(wind_beetle ~ previous_sum_ips + poly(previous_veg_tmp,2) + previous_sum_ips:previous_veg_tmp, 
                                data = dat_lag_scaled_complete)


# evaluate models

AIC(m_beetle1, m_beetle2, m_beetle3, m_beetle4,m_beetle6 , m_beetle6.2, m_beetle7, m_beetle8, m_beetle_interaction1, m_beetle_interaction2)
AICc(m_beetle1, m_beetle2, m_beetle3, m_beetle4,m_beetle6 ,m_beetle6.2, m_beetle7, m_beetle8, m_beetle_interaction1, m_beetle_interaction2)
BIC(m_beetle1, m_beetle2, m_beetle3, m_beetle4,m_beetle6 ,m_beetle6.2, m_beetle7, m_beetle8, m_beetle_interaction1, m_beetle_interaction2)

simulationOutput <- simulateResiduals(fittedModel = m_beetle_interaction1.1, plot = T) # residuals are bnetter for m_beetle6


AIC(m_beetle_interaction1, m_beetle_interaction1.1, 
    m_beetle_interaction3,m_beetle_interaction4, m_beetle_interaction5)

windows()
layout(matrix(1:3, nrow = 1))
plot(allEffects(m_beetle_interaction1.1))



summary(m_beetle_interaction1.1)
r2(m_beetle_interaction1.1)




# Predict beetle sums/year by environmental predictors -------------------------
#

windows()

pairs(sum_ips ~ veg_tmp  + annual_spei6, dat_lag, panel = panel.smooth)
pairs(sum_ips ~ veg_tmp + veg_prcp + previous_spei6_2 + previous_spei6  + spei6, dat_lag, panel = panel.smooth)



# make glm - beetle sums by temp and spei
# test with the raw and scaled predictors;
# glm, - define the polynomial relationship
### explore plots ------------------------------------------------------------------------

## test GLM (no random effects) on raw data ---------------------------------------------------------------------

##apparoces: GLM - no random effects (glm.nb)
# random effects : (glmer.nb) - not working really well
# final: use glmmTMB -allow for specifying random structure and dependencies in data collection

# select only columns of interests, to not remove more data as necessary 
# eg lagged values, that were NA
dd <- dat_lag_scaled %>% 
  dplyr::select(c(sum_ips, 
                  agg_doy, 
                  peak_doy,
                  peak_diff,
                  veg_tmp,
                  spring_tmp,
                  previous_veg_tmp,
                  spei1,
                  previous_spei1,
                  previous_spei1_2,
                  spei3,
                  previous_spei3,
                  previous_spei3_2,
                  spei12,
                  spei6,
                  previous_spei6,
                  previous_spei6_2,
                  spei12,
                  previous_spei12,
                  previous_spei12_2,
                  pairID, trapID, year))




### simplify analysis: AVG get average per pairID --------------------------
dd_simpl <- dd %>% 
  ungroup(.) %>% 
  group_by(year, pairID) %>% 
  summarise(sum_ips = round(mean(sum_ips)), 
            veg_tmp =mean(veg_tmp),
            previous_veg_tmp =mean(previous_veg_tmp),
            spei1 = mean(spei1),
            previous_spei1 = mean(previous_spei1),
            previous_spei1_2 = mean(previous_spei1_2),
            spei3 = mean(spei3),
            previous_spei3 = mean(previous_spei3),
            previous_spei3_2 = mean(previous_spei3_2),
            spei6 = mean(spei6),
            previous_spei6 = mean(previous_spei6),
            previous_spei6_2 = mean(previous_spei6_2),
            
            spei12 = mean(spei12),
            previous_spei12 = mean(previous_spei12),
            previous_spei12_2 = mean(previous_spei12_2))




###### Predict Ips sum/year FINAL ONE! start simple with glmmTMB, then add random effects ============================
m1 <- glmmTMB(sum_ips ~ veg_tmp,
              family = nbinom2,
              data = dd)

m2 <- glmmTMB(sum_ips ~ veg_tmp +previous_spei3_2 ,
              family = nbinom2,
              data = dd)

m3 <- glmmTMB(sum_ips ~ veg_tmp +previous_spei12_2 ,
              family = nbinom2,
              data = dd)
# add random year
m4 <- glmmTMB(sum_ips ~ veg_tmp +previous_spei12_2 + (1|year),
              family = nbinom2,
              data = dd)


# add random pairs
m5 <- glmmTMB(sum_ips ~ veg_tmp +previous_spei12_2 + (1|year) + (1 | pairID),
              family = nbinom2,
              data = dd)

m6 <- glmmTMB(sum_ips ~ veg_tmp +previous_spei12_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)

# add interaction
m7 <- glmmTMB(sum_ips ~ veg_tmp*previous_spei12_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)


# use poly fr temprature
m8 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei12_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)

# FINAl !!!use spei3 instead of spei 12 - m9 is the best !!!!!
m9 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei3_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)

# remove the I term to siimplify plotting - check how much it is better?
m9.2 <- glmmTMB(sum_ips ~ veg_tmp  + previous_spei3_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)

# save the final model
fin.m.counts <- m9

AIC(m9, m9.2)



# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.counts, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)

p1.counts <- create_effect_plot(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Beetle counts", show_y_axis = TRUE)
p2.counts <- create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "", show_y_axis = FALSE)

# Arrange the plots side by side with the same size
p.effect.counts <- ggarrange(p1.counts, p2.counts, ncol = 2, align = "v")
print(p.effect.counts)



# use spei1 instead of spei3
m9.1 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei1_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)



# use both SPEIs spei3 instead of spei 12
m10 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei3_2 + previous_spei12_2 + (1|year) + (1 | pairID/trapID),
              family = nbinom2,
              data = dd)


# remove spei 12, use poly term for spei
m11 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei3_2 + I(previous_spei3_2^2) + (1|year) + (1 | pairID/trapID),
               family = nbinom2,
               data = dd)



m12 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + spei3 + (1|year) + (1 | pairID/trapID),
               family = nbinom2,
               data = dd)

m13 <- glmmTMB(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei3 + (1|year) + (1 | pairID/trapID),
               family = nbinom2,
               data = dd)


arrange(AIC(m1, m2, m3, m4, m5, m6, m7,m8, m9, m9.1, m10, m11, m12, m13), AIC)




m9.tmb <- m9

summary(m9)
summary(m9.2)

r2(m9)
r2(m9.2)
m9.tmb <- m9

simRs <- simulateResiduals(m9, plot = T)
simRs <- simulateResiduals(m9.2, plot = T)

plot(allEffects(m9))
testOutliers(m9)
acf(residuals(m9.tmb))

dw_result <- dwtest(residuals(m9.tmb) ~ fitted(m9.tmb))
(dw_result)





# > summary(m9.tmb)
# Family: nbinom2  ( log )
# Formula:          sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei3_2 + (1 | year) +      (1 | pairID/trapID)
# Data: dd
# 
# AIC      BIC   logLik deviance df.resid 
# 23796.9  23837.0 -11890.5  23780.9     1098 
# 
# Random effects:
#   
#   Conditional model:
#   Groups        Name        Variance Std.Dev.
# year          (Intercept) 0.06265  0.2503  
# trapID:pairID (Intercept) 0.04338  0.2083  
# pairID        (Intercept) 0.25363  0.5036  
# Number of obs: 1106, groups:  year, 7; trapID:pairID, 158; pairID, 79
# 
# Dispersion parameter for nbinom2 family (): 2.89 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       9.80490    0.11453   85.61  < 2e-16 ***
#   veg_tmp           0.18368    0.06350    2.89 0.003819 ** 
#   I(veg_tmp^2)      0.06230    0.01897    3.28 0.001021 ** 
#   previous_spei3_2 -0.14683    0.04239   -3.46 0.000532 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#








# remove additionsl NAs
dd_complete <- dd %>%
  na.omit()





##### include 2 SPEIs and lagged values: ---------------------------------------
pairs(sum_ips ~ veg_tmp + previous_spei3_2  + previous_spei12_2, dd_complete, panel = panel.smooth) 

m1 <-  glm.nb(sum_ips ~ veg_tmp, dd_complete) 

# inspect short term drought
m2 <-  glm.nb(sum_ips ~ veg_tmp + spei1, dd_complete) # spei3 is not importantt 
m3 <-  glm.nb(sum_ips ~ veg_tmp + spei1 + previous_spei1, dd_complete)   # one year lag is not iportant
m4 <-  glm.nb(sum_ips ~ veg_tmp + spei1 + previous_spei1_2, dd_complete)

m5 <-  glm.nb(sum_ips ~ veg_tmp + spei1 + previous_spei1 + previous_spei1_2 + previous_spei12_2, dd_complete)


m6 <-  glm.nb(sum_ips ~ veg_tmp + spei1 + previous_spei12_2, dd_complete)


AIC(m1, m2, m3, m4, m5)
summary(m6)
simulateOutput <- simulateResiduals(m6_1, plot = T)

# addd long term grought
m6 <-  glm.nb(sum_ips ~ veg_tmp + previous_spei6_2 + spei12, dd) # spei12 is not important

m6_1 <- glmer.nb(sum_ips ~ veg_tmp + spei3 + previous_spei3_2 + (1 | pairID), data = dd)

AIC(m6, m6_1)
m7 <-  glm.nb(sum_ips ~ veg_tmp + previous_spei1_2 + previous_spei12, dd_complete)
m8 <-  glm.nb(sum_ips ~ veg_tmp + previous_spei1_2 + previous_spei12 + previous_spei12_2, dd_complete)
m9 <-  glm.nb(sum_ips ~ veg_tmp + previous_spei1_2  + previous_spei12_2, dd)  # remove non signoficant lag1

# add poly or quadratic tersm
m9_1 <-  glm.nb(sum_ips ~ veg_tmp + I(veg_tmp^2) + previous_spei1_2  + previous_spei12_2, dd_complete)  # remove non signoficant lag1

m9_interaction <- glm.nb(formula = sum_ips ~ veg_tmp * previous_spei3_2 + veg_tmp * previous_spei12_2, 
                         data = dd_complete, init.theta = 1.773405852, link = log)


m9_interaction2 <- glm.nb(formula = sum_ips ~ veg_tmp+ previous_spei3_2*previous_spei12_2, 
                         data = dd_complete, init.theta = 1.773405852, link = log)


m9_poly <-  glm.nb(sum_ips ~poly(veg_tmp,2) + poly(previous_spei3_2,2)  + poly(previous_spei12_2,2), dd_complete)  # remove non signoficant lag1

m9_poly2 <-  glm.nb(sum_ips ~I(veg_tmp^2) + poly(previous_spei3_2,2)  + poly(previous_spei12_2,2), dd_complete)  # remove non signoficant lag1



# remove tmp
m10 <-  glm.nb(sum_ips ~ previous_spei3_2  + previous_spei12_2, dd_complete)  # remove non signoficant lag1


# simplify, keep only spei12
m11 <-  glm.nb(sum_ips ~ veg_tmp +  previous_spei12_2, dd_complete)  # does not improve the model, I should keep both SPEIs


AIC(m11,m9, m9_poly2)



# add random effects
mm8 <-  glmer.nb(sum_ips ~ veg_tmp + previous_spei3_2 + previous_spei12_2 +  (1 | pairID), dd_complete)  # too complex
mm9 <-  glmer.nb(sum_ips ~ veg_tmp + previous_spei3_2 + previous_spei12 + previous_spei12_2 +  (1 | year), dd)    #
mm10 <-  glmer.nb(sum_ips ~ veg_tmp + previous_spei12_2 +  (1 | pairID), dd_complete) # too complex

vif(m9_poly)

AIC(m9.tmb, m9, mm9)

AICc(m5, m6, m7, m8, m9, m10, m9_1, m9_interaction, m9_poly, m9_poly2, m9_interaction2, m11)
summary(m9_poly2)
r2(m9)
plot(allEffects(m9_interaction2))
simulateResiduals(fittedModel = m9  , plot = T)

mm1 <- glm.nb(sum_ips ~ veg_tmp + spei3 + previous_spei3+ previous_spei3_2 + spei12+ previous_spei12+ previous_spei12_2, dd_complete) 


glmm1 <- glmer.nb(sum_ips ~ veg_tmp + previous_spei3_2 + previous_spei12_2 + 
                    #    (1 | pairID) + 
                    (1 | year), 
                      data = dd_complete,
                      na.action = 'na.fail'#,
                      #family = negative.binomial(2.8588)
                      #, na.action = "na.omit"  "na.fail"
)
# remove year as random
glmm2 <- glmer(sum_ips ~ veg_tmp + previous_spei12_2 + #+ spei12+ previous_spei12
                 (1 | pairID) , 
               data = dd_complete,
               na.action = 'na.fail',
               family = negative.binomial(2.4025)
               #, na.action = "na.omit"  "na.fail"
)

glmm3 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + poly(previous_spei12_2,2) + 
                    (1 | pairID), 
                  data = dd_complete, 
                  na.action = 'na.fail')

summary(glmm3)
AICc(mm1, glmm1,glmm2, glmm3, m9_interaction2,glmm3  )

simulationOutput <- simulateResiduals(fittedModel = glmm3, plot = T)
testOutliers(glmm2, type = 'bootstrap') 
plot(allEffects(glmm2))


# Variance Inflation Factor for fixed effects
vif_values <- vif(glmm1)

(vif_values)

# Autocorrelation check
acf(residuals(m9.tmb))
dw_result <- dwtest(residuals(m9.tmb) ~ fitted(m9.tmb))
(dw_result)

summary(m3.7.1.aut6)



# predict DOY aggregation ---------------------------------------------------


### simplify analysis: AVG get average per pairID --------------------------
dd_simpl <- dd %>% 
  ungroup(.) %>% 
  group_by(year, pairID) %>% 
  summarise(agg_doy = round(mean(agg_doy)),
            peak_doy = round(mean(peak_doy)),
            spring_tmp =mean(spring_tmp),
            veg_tmp = mean(veg_tmp),
            spei1 = mean(spei1),
            previous_spei1 = mean(previous_spei1),
            previous_spei1_2 = mean(previous_spei1_2),
            spei3 = mean(spei3),
            previous_spei3 = mean(previous_spei3),
            previous_spei3_2 = mean(previous_spei3_2),
            spei6 = mean(spei6),
            previous_spei6 = mean(previous_spei6),
            previous_spei6_2 = mean(previous_spei6_2),
            
            spei12 = mean(spei12),
            previous_spei12 = mean(previous_spei12),
            previous_spei12_2 = mean(previous_spei12_2))






# try modelsL # if bounded values, skewed: potentially beta regression with tranformation
dd$tr_agg_doy  <- (dd$agg_doy - 60) / (304 - 60)
dd$tr_peak_doy <- (dd$peak_doy - 60) / (304 - 60)

# Fit a beta regression model

library(glmmTMB)

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dd$tr_agg_doy  <- pmin(pmax(dd$tr_agg_doy, 1e-4),  1 - 1e-4)
dd$tr_peak_doy <- pmin(pmax(dd$tr_peak_doy, 1e-4), 1 - 1e-4)




# Transform your response variable to fit within 0 to 1
#dd_simpl$tr_agg_doy  <- (dd_simpl$agg_doy - 60) / (300 - 60)
#dd_simpl$tr_peak_doy <- (dd_simpl$peak_doy - 60) / (300 - 60)

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
#dd_simpl$tr_agg_doy  <- pmin(pmax(dd_simpl$tr_agg_doy, 1e-4),  1 - 1e-4)
#dd_simpl$tr_peak_doy <- pmin(pmax(dd_simpl$tr_peak_doy, 1e-4), 1 - 1e-4)


# try GLM aggregation  -----------------

hist(dd$agg_doy)
hist(dd$peak_doy)


m.agg0 <- glmmTMB(tr_agg_doy ~ veg_tmp + spei1  + (1 | pairID),
                   family = beta_family(link = "logit"),
                   dd)

m.agg01 <- glmmTMB(tr_agg_doy ~ veg_tmp + previous_spei3_2  + (1 | pairID),
                  family = beta_family(link = "logit"),
                  dd)

m.agg02 <- glmmTMB(tr_agg_doy ~ spring_tmp + previous_spei3_2  + (1 | pairID),
                   family = beta_family(link = "logit"),
                   dd)


m.agg1 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 ,
                  family = beta_family(link = "logit"),
                  data = dd)


m.agg2 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 + previous_spei3_2 ,
                  family = beta_family(link = "logit"),
                  data = dd)

m.agg3 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 + previous_spei3 + previous_spei3_2 ,
                  family = beta_family(link = "logit"),
                  data = dd)

# add random effects
m.agg4 <- glmmTMB(tr_agg_doy ~ spring_tmp + previous_spei3_2 + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)

m.agg5 <- glmmTMB(tr_agg_doy ~ spring_tmp + previous_spei3_2 + previous_spei12_2 + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)

# add nested trap within pairID - does not improve model
m.agg6 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 + previous_spei3_2 + previous_spei12_2 + (1 | pairID/trapID),
                  family = beta_family(link = "logit"),
                  data = dd)

# add spei
m.agg7 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 +previous_spei3 + previous_spei3_2 + previous_spei12_2 +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)


m.agg8 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 +previous_spei3 + previous_spei3_2 + previous_spei12_2 +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)


m.agg9 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 +previous_spei1 + previous_spei1_2 + previous_spei12_2 +(1 | pairID/trapID),
                  family = beta_family(link = "logit"),
                  data = dd)

m.agg10 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 +previous_spei1 + previous_spei1_2 + previous_spei12_2 + (1 | pairID/trapID) +
                    (1|year),
                  family = beta_family(link = "logit"),
                  data = dd)

# remove non signoificant vars - ned to keep random structure!  m.agg11 is teh best
m.agg11 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 + (1 | pairID), #+
                   #  (1|year)
                   family = beta_family(link = "logit"),
                   data = dd)

m.agg11.sq <- glmmTMB(tr_agg_doy ~ spring_tmp + I(spring_tmp^2) + spei1 + (1 | pairID)   ,
                   family = beta_family(link = "logit"),
                   data = dd)

m.agg11.poly <- glmmTMB(tr_agg_doy ~ poly(spring_tmp,2) + spei1 + (1 | pairID),
                      family = beta_family(link = "logit"),
                      data = dd)

# use SPEI3
m.agg11_1 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 + (1 | pairID),
                   family = beta_family(link = "logit"),
                   data = dd)
# exclude random eff
m.agg12 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1,
                   family = beta_family(link = "logit"),
                   data = dd)

# exclude random eff, add lagged vals
m.agg12.2 <- glmmTMB(tr_agg_doy ~ veg_tmp + previous_spei3_2 + (1| pairID),
                   family = beta_family(link = "logit"),
                   data = dd)

r2(m.agg12.2)
cor(dd$spring_tmp,dd$previous_spei3_2, method = "spearman")

m.agg12.3 <- glmmTMB(tr_agg_doy ~ poly(veg_tmp,2) + poly(previous_spei3_2,2) + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)


m.agg12.3.12 <- glmmTMB(tr_agg_doy ~ poly(spring_tmp,2) + spei1 + poly(previous_spei12_2,2) + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)

m.agg12.4 <- glmmTMB(tr_agg_doy ~ poly(spring_tmp,2) + previous_spei3_2 + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)

m.agg12.5 <- glmmTMB(tr_agg_doy ~ spring_tmp + poly(previous_spei3_2,2) + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)

m.agg12.6 <- glmmTMB(tr_agg_doy ~ poly(spring_tmp,2) + poly(previous_spei3_2,2) + as.factor(year) +(1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)



m.agg12.spei1 <- glmmTMB(tr_agg_doy ~ poly(spring_tmp,2) + poly(spei1,2) +(1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)


# try interaction - has worse AIC, go for the m.agg12.3
m.agg12.3.int <- glmmTMB(tr_agg_doy ~ spring_tmp*previous_spei3_2 + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)




# simplify random effects
m.agg13 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei1 + (1|pairID),
                   family = beta_family(link = "logit"),
                   data = dd)

# use SPEI3
m.agg14 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3 + (1|pairID),
                   family = beta_family(link = "logit"),
                   data = dd)

# use previous speis
m.agg15 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3+ previous_spei3 + (1|pairID),
                   family = beta_family(link = "logit"),
                   data = dd)


# use previous speis
m.agg16 <- glmmTMB(tr_agg_doy ~ spring_tmp + spei3+ previous_spei3 + previous_spei3_2 + (1|pairID),
                   family = beta_family(link = "logit"),
                   data = dd)




summary(m.agg12.3)
AIC(m.agg1, m.agg2, m.agg3, m.agg4, m.agg5, m.agg6, m.agg7, m.agg8, m.agg9, m.agg10, m.agg11,m.agg11_1, m.agg12)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.agg12.3, #m.agg10, #,#m.agg12.3,  #m.agg01 
                                              plot = T)
AICc( m.agg01,m.agg13, m.agg14, m.agg15, m.agg16, m.agg12, m.agg12.2, m.agg12.3, m.agg12.4,m.agg12.5, m.agg12.3.int, m.agg12.3.12, m.agg0)


AIC(m.agg12.2, m.agg12.3)
windows()
plot(allEffects(m.agg12.3 ))
r2(m.agg12.3)
summary(m.agg12.3)

acf(residuals(m.agg12.3))

(dw_result <- dwtest(residuals(m.agg12.3) ~ fitted(m.agg12.3)))


summary(m.agg12.3)  # is teh winner!! m.agg12.3



# > summary(m.agg12.3)  # is teh winner!! m.agg12.3
# Family: beta  ( logit )
# Formula:          tr_agg_doy ~ poly(veg_tmp, 2) + poly(previous_spei3_2, 2) + (1 |      pairID)
# Data: dd
# 
# AIC      BIC   logLik deviance df.resid 
# -2272.2  -2237.3   1143.1  -2286.2     1073 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# pairID (Intercept) 0.04945  0.2224  
# Number of obs: 1080, groups:  pairID, 79
# 
# Dispersion parameter for beta family (): 32.8 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                -0.75822    0.02748 -27.587  < 2e-16 ***
#   poly(veg_tmp, 2)1          -6.81460    0.61985 -10.994  < 2e-16 ***
#   poly(veg_tmp, 2)2          -3.18407    0.48070  -6.624 3.50e-11 ***
#   poly(previous_spei3_2, 2)1  5.38866    0.41051  13.127  < 2e-16 ***
#   poly(previous_spei3_2, 2)2  1.73480    0.39732   4.366 1.26e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# best model ------------------------------------------------------------

fin.m.agg <- m.agg12.3



### test fr peak population ----------------------------------------------------


# add random effects
m.peak4 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_spei6  + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)

m.peak5 <- glmmTMB(tr_peak_doy ~ poly(veg_tmp,2) + previous_spei6 + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)


# add nested trap within pairID - does not improve model
m.peak6 <- glmmTMB(tr_peak_doy ~ veg_tmp + poly(previous_spei6,2) + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)

# add spei
m.peak7 <- glmmTMB(tr_peak_doy ~ veg_tmp  + spei6 + previous_spei6 + previous_spei6_2 + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dd)


# remove non significant
m.peak8 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_spei6 + previous_spei6_2 + (1 | pairID),
                   family = beta_family(link = "logit"),
                   dd)
                   

# remove less important spei:
m.peak9 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_spei6  + (1 | pairID),
                   family = beta_family(link = "logit"),
                   dd)



AIC(m.peak4, m.peak5, m.peak6, m.peak7, m.peak8, m.peak9)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.peak9, 
                                              plot = T)

fin.m.peak <- m.peak9

windows()
plot(allEffects(m.peak9 ))
r2(m.peak9)

summary(m.peak9)

acf(residuals(m.peak9))

(dw_result <- dwtest(residuals(m.peak9) ~ fitted(m.peak9)))


summary(m.peak9)  # is teh winner!! 

# > summary(m.peak9)
# Family: beta  ( logit )
# Formula:          tr_peak_doy ~ veg_tmp + previous_spei6 + (1 | pairID)
# Data: dd
# 
# AIC      BIC   logLik deviance df.resid 
# -1672.1  -1647.1    841.1  -1682.1     1101 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# pairID (Intercept) 0.006043 0.07773 
# Number of obs: 1106, groups:  pairID, 79
# 
# Dispersion parameter for beta family (): 18.6 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.19810    0.01624 -12.200  < 2e-16 ***
#   veg_tmp        -0.09209    0.01518  -6.066 1.31e-09 ***
#   previous_spei6 -0.10219    0.01405  -7.275 3.45e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



cor(dd$previous_spei6,dd$veg_tmp )
#cor(dd$previous_spei3_2,dd$spring_tmp )

# GLM Peak difference -----------------------------------------------------
#pairs(peak_diff ~ spri)

# convert first to integer values: difference betweeen two consecutive dates

dd$peak_diff <- as.integer(dd$peak_diff)
cor(dd$spei6, dd$previous_spei3_2)

hist(dd$peak_diff)

m.peak.diff1 <- glmmTMB(peak_diff ~ spring_tmp,
              family = nbinom2,
              data = dd)

m.peak.diff2 <- glmmTMB(peak_diff ~ spring_tmp + spei6 ,
              family = nbinom2,
              data = dd)


m.peak.diff2.1 <- glmmTMB(peak_diff ~ spei6 ,
                        family = nbinom2,
                        data = dd)


m.peak.diff2.2 <- glmmTMB(peak_diff ~ veg_tmp + spei6 +(1 | pairID),
                          family = nbinom2,
                          data = dd)

m.peak.diff2.3 <- glmmTMB(peak_diff ~ spring_tmp + previous_spei3_2  +(1 | pairID),
                          family = nbinom2,
                          data = dd)

m.peak.diff2.4 <- glmmTMB(peak_diff ~ veg_tmp + spei3 + previous_spei3 + previous_spei3_2  +(1 | pairID),
                          family = nbinom2,
                          data = dd)

m.peak.diff2.5 <- glmmTMB(peak_diff ~ veg_tmp + previous_spei3 + previous_spei3_2  +(1 | pairID),
                          family = nbinom2,
                          data = dd)

m.peak.diff2.6 <- glmmTMB(peak_diff ~ veg_tmp  + previous_spei3_2  +(1 | pairID),
                          family = nbinom2,
                          data = dd)


AICc(m.peak.diff2.3, m.peak.diff2.4, m.peak.diff2.5,m.peak.diff2.6)
summary(m.peak.diff2.6)  # teh best!!!!
fin.m.peak.diff <- m.peak.diff2.6

cor(dd$spring_tmp, dd$previous_spei3_2)

AIC(m.peak.diff2.1, m.peak.diff2,m.peak.diff1,m.peak.diff2.2,m.peak.diff2.3)

# add random trap effect
m.peak.diff3 <- glmmTMB(peak_diff ~ spring_tmp + spei6 +  (1|year)+(1 | pairID),
                        family = nbinom2,
                        data = dd)

# add random trap effect
m.peak.diff3.1 <- glmmTMB(peak_diff ~ spring_tmp + previous_spei3_2 +(1 | pairID),
                        family = nbinom2,
                        data = dd)


# add random year
m.peak.diff4 <- glmmTMB(peak_diff ~ spring_tmp + spei6 + (1|year),
              family = nbinom2,
              data = dd)

# add random trap effect
m.peak.diff5 <- glmmTMB(peak_diff ~ spring_tmp + spei6 +  (1 | pairID),
              family = nbinom2,
              data = dd)

# add random trap effect
m.peak.diff5.1 <- glmmTMB(peak_diff ~ spring_tmp + previous_spei6_2 +  (1 | pairID),
                        family = nbinom2,
                        data = dd)


m.peak.diff5.2 <- glmmTMB(peak_diff ~ spring_tmp + previous_spei3_2 +  (1 | pairID),
                          family = nbinom2,
                          data = dd)


m.peak.diff5.3 <- glmmTMB(peak_diff ~ spring_tmp + I(spring_tmp^2) + previous_spei3_2 +  (1 | pairID),
                          family = nbinom2,
                          data = dd)
AIC(m.peak.diff5, m.peak.diff5.1,m.peak.diff5.2,m.peak.diff5.3)

# add interaction
m.peak.diff6 <- glmmTMB(peak_diff ~ spring_tmp*spei6  + (1 | pairID),
              family = nbinom2,
              data = dd)


# use poly fr temprature
m.peak.diff7 <- glmmTMB(peak_diff ~ spring_tmp + I(spring_tmp^2) + spei6 + (1 | pairID),
              family = nbinom2,
              data = dd)



# use annual temp and previous spei3 (as in teh count model)
m.peak.diff8 <- glmmTMB(peak_diff ~ veg_tmp + I(veg_tmp^2) + previous_spei3_2  + (1 | pairID),
                        family = nbinom2,
                        data = dd)

# ass SPEI 6 
m.peak.diff9 <- glmmTMB(peak_diff ~ veg_tmp + I(veg_tmp^2) + spei6 + previous_spei3_2  + (1 | pairID),
                        family = nbinom2,
                        data = dd)


m.peak.diff10 <- glmmTMB(peak_diff ~ spring_tmp + I(spring_tmp^2) + spei3 + previous_spei3_2  + (1 | pairID),
                        family = nbinom2,
                        data = dd)


m.peak.diff11 <- glmmTMB(peak_diff ~ spring_tmp + I(spring_tmp^2) + spei6 + previous_spei6_2  + (1 | pairID),
                         family = nbinom2,
                         data = dd)

m.peak.diff12 <- glmmTMB(peak_diff ~ spring_tmp + I(spring_tmp^2) + previous_spei6_2  + (1 | pairID),
                         family = nbinom2,
                         data = dd)

m.peak.diff13 <- glmmTMB(peak_diff ~ spring_tmp  + previous_spei6_2  + (1 | pairID),
                         family = nbinom2,
                         data = dd)


arrange(AIC(m.peak.diff1,m.peak.diff2,m.peak.diff3,m.peak.diff4,m.peak.diff5,m.peak.diff6,m.peak.diff7,m.peak.diff8,m.peak.diff9, m.peak.diff10, m.peak.diff11,m.peak.diff12,m.peak.diff13,
            
            m.peak.diff5.1,m.peak.diff5.2,m.peak.diff5.3,
            m.peak.diff3.1, 
            m.peak.diff2.1), AIC)



summary(m.peak.diff2.3)
r2(m.peak.diff2.6)
simRs <- simulateResiduals(m.peak.diff2.3, plot = T)
plot(allEffects(m.peak.diff2.3))
testOutliers(m.peak.diff8)
acf(residuals(m.peak.diff2.6))

dw_result <- dwtest(residuals(m.peak.diff2.6) ~ fitted(m.peak.diff2.6))
(dw_result)




# save models -------------------------------------------------------------

# 
# 
# save( dd,
#       fin.m.agg,
#       fin.m.counts,
#       fin.m.peak,
#       fin.m.peak.diff,
#       file="outData/fin_models.Rdata")
# 
# 
# 


# Effect plots ------------------------------------------------------------


# define plotiing functions


# Create effect plot function with an additional argument to control y-axis labels
create_effect_plot <- function(data, line_color = "blue", x_title = "X-axis", y_title = "Y-axis", y_lim = c(100,80000), show_y_axis = TRUE,  my_title = '') {
  p <- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title,
         title = my_title) +
    ylim(y_lim) +
  #  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.7, 2.7)) + # Set x-axis breaks and limits
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black"))
  
  if (show_y_axis) {
    p <- p + labs(y = y_title)
  } else {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  }
  
  return(p)
}




# change y-axis into days (from 0-1)
create_effect_Y_trans <- function(data, line_color = "blue", x_title = "X-axis", y_title = "Y-axis", y_lim = c(0, 1), show_y_axis = TRUE, my_title = '') {
  data$predicted <- (data$predicted * (304 - 60)) + 60  # Reverse transformation for y-axis
  data$conf.low  <- (data$conf.low * (304 - 60)) + 60
  data$conf.high <- (data$conf.high * (304 - 60)) + 60
  
  p<- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title, y = y_title, title = my_title) +
    ylim(y_lim) +
  #  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.7, 2.7)) + # Set x-axis breaks and limits
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black"))
  
  if (show_y_axis) {
    p <- p + labs(y = y_title)
  } else {
    p <- p + theme(axis.title.y = element_blank(), 
                   axis.text.y = element_blank())
  }
}



# Beetle counts -----------------------------------------------------------


# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.counts, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)

# change the values to allow y lables to fit 
p1$predicted <- p1$predicted/100
p1$conf.low <- p1$conf.low/100
p1$conf.high <- p1$conf.high/100

p2$predicted <- p2$predicted/100
p2$conf.low <- p2$conf.low/100
p2$conf.high <- p2$conf.high/100

p1.counts <- create_effect_plot(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Beetle population level (*100)", show_y_axis = TRUE, y_lim = c(80,800))
p2.counts <- create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Beetle population level (*100)", show_y_axis = TRUE, y_lim = c(80,800))

# Arrange the plots side by side with the same size
p.effect.counts <- ggarrange(p1.counts, p2.counts, ncol = 1, align = "hv")
(p.effect.counts)





# DOY aggregation ---------------------------------------------------------

p1 <- ggpredict(fin.m.agg, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)


p1.agg <- create_effect_Y_trans(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Colonization DOY", y_lim = c(80,250),show_y_axis = TRUE  )
p2.agg <- create_effect_Y_trans(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Colonization DOY" , y_lim = c(80,250), show_y_axis = TRUE)



p.effect.agg <- ggarrange(p1.agg,p2.agg,ncol = 1, align = "hv")

(p.effect.agg)

###### PEak Effect plots ------------------------------------------------------------

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak, terms = "previous_spei6 [all]", allow.new.levels = TRUE)


p1.peak <- create_effect_Y_trans(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Peak population growth DOY", y_lim = c(130,200), show_y_axis = TRUE )
p2.peak <-create_effect_Y_trans(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Peak population growth DOY", y_lim = c(130,200), show_y_axis = TRUE)

# effect plots  ----------------------------------------------------
p.effect.peak <- ggarrange(p1.peak,p2.peak, ncol = 1, align = "hv")

p.effect.peak



# effect plot peak dif ----------------------------------------------------
# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak.diff, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)


p1.peak.diff <- create_effect_plot(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Peak difference DOY", y_lim = c(220,650), show_y_axis = TRUE)
p2.peak.diff <-create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Peak difference DOY", y_lim = c(220,650), show_y_axis = TRUE)

# effect plots Sum IPS ----------------------------------------------------
p.effect.peak.diff <- ggarrange(p1.peak.diff,p2.peak.diff, ncol = 1, align = "hv")

p.effect.peak.diff


# all Effect plots --------------------------------------------------------
windows(9,5)
ggarrange(p.effect.counts, p.effect.agg, p.effect.peak, p.effect.peak.diff, ncol=4, nrow = 1 , align = 'hv', font.label = list(size = 10, color = "black", face = "bold", family = NULL),
          labels = c( "[a]","[b]","[c]","[d]"))


# get characteristics
r2(fin.m.counts )
r2(fin.m.agg)
r2(fin.m.peak)
r2(fin.m.peak.diff)



summary(fin.m.counts )
summary(fin.m.agg)
summary(fin.m.peak)
summary(fin.m.peak.diff)



# print models outputs:



# export them all in one document:
# Temporarily save tables as HTML
sjPlot::tab_model(fin.m.counts, file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg, file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak, file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff, file = "outTable/model_peak_diff.doc")





# Effect plots RS ----------------------------------------------------------



# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.RS, terms = "previous_sum_ips [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.RS, terms = "previous_veg_tmp [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.RS, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m.RS, terms = "population_growth2 [all]", allow.new.levels = TRUE)


p1.RS <- create_effect_plot(p1, line_color = "red", x_title = "Beetle population level (lag1) [z-score]", y_title = "Tree mortality [# pixels]", y_lim = c(0,1000), show_y_axis = TRUE)
p2.RS <- create_effect_plot(p2, line_color = "blue", x_title = "Temperature (lag1)  [z-score]", y_title = "Tree mortality [# pixels]", y_lim = c(0,320), show_y_axis = TRUE)
p3.RS <- create_effect_plot(p3, line_color = "yellow", x_title = "SPEI (lag2) [z-score]", y_title = "Tree mortality [# pixels]", y_lim = c(0,20), show_y_axis = TRUE)
p4.RS <- create_effect_plot(p4, line_color = "grey", x_title = "Popul. growth (lag2) [z-score]", y_title = "Tree mortality [# pixels]", y_lim = c(0,50), show_y_axis = TRUE)



# effect plots Sum IPS ----------------------------------------------------
p.effect.RS <- ggarrange(p1.RS,p2.RS,p3.RS,p4.RS, ncol = 2, nrow = 2, align = 'hv')

p.effect.RS


# calcualte DOYs for teh temp-z acroe: aggergation day

# Coefficients from the model
coef1 = -6.76424
coef2 = -3.15288
intercept = -0.78305

# Temperature z-scores
z_scores <- c(-2, -1, 0, 1, 2)

# Function to calculate predicted DOY
calculate_DOY <- function(z_score) {
  # Calculate the polynomial terms
  poly1 = z_score  # First polynomial term (linear)
  poly2 = z_score^2  # Second polynomial term (squared)
  
  # Calculate the linear predictor
  linear_predictor = intercept + coef1 * poly1 + coef2 * poly2
  
  # Applying logistic transformation
  predicted_proportion = exp(linear_predictor) / (1 + exp(linear_predictor))
  
  # Converting the proportion back to the DOY scale
  predicted_DOY = predicted_proportion * (304 - 60) + 60
  return(predicted_DOY)
}

# Calculate DOY for each z-score
predicted_DOYs <- sapply(z_scores, calculate_DOY)

# Create a data frame for the table
doys_table <- data.frame(Temperature_Z_Score = z_scores, Estimated_DOY = predicted_DOYs)
doys_table


