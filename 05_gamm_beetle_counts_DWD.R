


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
  dplyr::select(c(year, trapID, sum_ips, 
                  Morans_I, # calculated from beetles sums 
                  Morans_I_log)) # calculated from log(beetles sum)

# lisa_merged_agg_doy_df <- 
#   lisa_merged_agg_doy_df %>% 
#   as.data.frame() %>% 
#     dplyr::rename(Morans_I_agg = Morans_I) %>% 
#   dplyr::select(year, trapID, agg_doy, Morans_I_agg)
  


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
                  previous_peak_doy1 =  dplyr::lag(peak_doy , order_by = year, n = 1), 
                  previous_peak_doy2 =  dplyr::lag(peak_doy , order_by = year, n = 2),
                  previous_veg_tmp    =  dplyr::lag(veg_tmp,    order_by = year),
                  previous_veg_tmp2    =  dplyr::lag(veg_tmp,    order_by = year, n = 2),
                  #previous_spring_tmp =  dplyr::lag(spring_tmp, order_by = year),
                  previous_veg_prcp   =  dplyr::lag(veg_prcp,   order_by = year),
                  #previous_spei1      =  dplyr::lag(spei1,      order_by = year),
                  previous_spei3      =  dplyr::lag(spei3,      order_by = year),
                  #previous_spei6      =  dplyr::lag(spei6,     order_by = year),
                  #previous_spei12     =  dplyr::lag(spei12,     order_by = year),
                  #previous_spei1_2    =  dplyr::lag(spei1,      order_by = year, n = 2),
                  previous_spei3_2    =  dplyr::lag(spei3,      order_by = year, n = 2),
                  #previous_spei6_2    =  dplyr::lag(spei6,     order_by = year, n = 2),
                  #previous_spei12_2   =  dplyr::lag(spei12,     order_by = year, n = 2)
                  ) %>% 
  dplyr::mutate(population_growth     = (sum_ips - previous_sum_ips) / previous_sum_ips * 100,
                population_growth2    = dplyr::lag(population_growth, order_by = year)) %>%  # lag population growth by one more year
  left_join(df_RS_out, 
            by = join_by(trapID, year, sum_ips)) %>%
  left_join(lisa_merged_df, 
            by = join_by(trapID, year, sum_ips)) %>%
#  left_join(lisa_merged_agg_doy_df, 
#            by = join_by(trapID, year, agg_doy)) %>%
    dplyr::mutate(previous_Moran       =  dplyr::lag(Morans_I , order_by = year),
                previous_Moran2      =  dplyr::lag(Morans_I , order_by = year, n = 2),
                previous_Moran_log1       =  dplyr::lag(Morans_I_log , order_by = year),
                previous_Moran_log2      =  dplyr::lag(Morans_I_log , order_by = year, n = 2),
                ) %>% 
  dplyr::mutate(trapID = as.factor(trapID)) %>%
  dplyr::filter(year %in% 2015:2021)

# check morans from ips sum and agg doy:
#pairs(Morans_I ~ Morans_I_agg, dat_lag )


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


# Scale predictors ================================================================================

# skip columns if not for scaling
#spei_cols <- grepl("spei", names(dat_lag))

additional_skip_cols <- c('trapID', 'pairID', "x", "y", 
                          'year', 'agg_doy', 'peak_doy', 'peak_diff',
                          'sum_ips', 'previous_sum_ips', 'previous_sum_ips2', 
                          'wind_beetle', 'Morans_I', 'Morans_I_log')

# Combine spei columns with additional columns to skip
skip_col <- additional_skip_cols #c(names(dat_lag)[spei_cols], additional_skip_cols)

# scale also SPEIs

# export new df
dat_lag_scaled <-
  dat_lag %>%
  ungroup(.) %>% 
 # na.omit() %>% 
  mutate(sc_sum_ips = sum_ips,  # add beetle counts as scaled values
         sc_previous_sum_ips = previous_sum_ips,
         sc_previous_sum_ips2 = previous_sum_ips2) %>% 
  dplyr::select(-all_of(skip_col )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_lag %>% dplyr::select(all_of(skip_col)), .)

# #simplify, keep all columns as vectors 
# dat_lag_scaled <- dat_lag_scaled %>% 
#   ungroup(.) %>% 
#   dplyr::mutate(across(where(~is.numeric(.) && is.matrix(.)), as.vector))


# remove additionsl NAs
dat_lag_scaled_complete <- 
  dat_lag_scaled %>%
  dplyr::select(-anom_wind_beetle, -anom_all_dist, - anom_harvest) %>%
  na.omit()
#  filter(across(everything(), ~ !is.na(.)))



#### get a table for RS data --------------------

# insspect beetle rate: get cumulative rates!
skip_cols_RS <- c('trapID', 'pairID', 'year', 'wind_beetle', 'cumSum_beetlerate')


dat_lag_RS <- dat_lag %>% 
  dplyr::select('trapID', 'pairID', 'year', 'wind_beetle',
    'agg_doy', 'peak_doy', 'peak_diff', 
                'sum_ips', 'previous_sum_ips', 'previous_sum_ips2',
                "previous_agg1" ,      "previous_agg2",
                "previous_peak_diff1", "previous_peak_diff2", 
                "veg_tmp", "previous_veg_tmp",  "previous_veg_tmp2",
                "spei3",  "previous_spei3", "previous_spei3_2",
                "population_growth",   "population_growth2",
    "previous_Moran", "previous_Moran2", 
    "Morans_I", 
    "previous_Moran_log1", "previous_Moran_log2", 
    "Morans_I_log",
    "beetle_rate") %>% 
  mutate(beetle_sum2yrs = previous_sum_ips + previous_sum_ips2) %>%   # new variable having cumulative values ob beetle population per two years
  group_by(trapID) %>% 
  arrange(year) %>% 
  mutate(cumSum_beetlerate = base::cumsum(beetle_rate)) %>% 
  ungroup(.)

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

pairs(Morans_I_agg ~  previous_veg_tmp + previous_sum_ips + spei1, dat_lag_scaled_complete)

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

# perform EDA: MOrans_log -----------------------------
library(GGally)
ggpairs(dat_lag_scaled[, c("Morans_I_log", "veg_tmp", "spei3", "previous_spei3_2")])


### explain Moran's I: LM ------------------------------------
m_moran1 <- lm(Morans_I ~ poly(veg_tmp,2) + previous_sum_ips + poly(spei3,2), data = dat_lag_scaled)

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

# filter only the most important predictors
pairs(Morans_I~ previous_spei3_2 + 
        previous_veg_tmp + 
        previous_agg2 +
        #previous_peak_doy2 +
        previous_peak_diff2 +
        population_growth2
        , dat_lag_scaled_complete)


pairs(Morans_I_agg~ previous_spei3_2 + 
        previous_veg_tmp + 
        previous_agg2 +
        #previous_peak_doy2 +
        previous_peak_diff2 +
        population_growth2
      , dat_lag_scaled_complete)


m_moran17 <- glmmTMB(Morans_I ~ 
                       #spei3 +
                       #previous_spei3 + 
                       previous_spei3_2   + 
                       #veg_tmp + 
                       previous_veg_tmp + 
                       #sc_previous_sum_ips + 
                       #previous_agg1 + 
                       previous_agg2 +
                       previous_peak_doy2 +
                       #previous_peak_diff1  + 
                       previous_peak_diff2 +
                       population_growth2  + 
                    (1 | pairID),  data = dat_lag_scaled_complete,
                    na.action = "na.fail")


m_moran18 <- glmmTMB(Morans_I ~ 
                       #spei3 +
                       #previous_spei3 + 
                       previous_spei3_2   + 
                       #veg_tmp + 
                       poly(previous_veg_tmp,2) + 
                       #sc_previous_sum_ips + 
                       #previous_agg1 + 
                       previous_agg2 +
                       previous_peak_doy2 +
                       #previous_peak_diff1  + 
                       previous_peak_diff2 +
                       population_growth2  + 
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")



m_moran19 <- glmmTMB(Morans_I ~ 
                       #spei3 +
                       #previous_spei3 + 
                       previous_spei3_2   + 
                       #veg_tmp + 
                       previous_veg_tmp + 
                       #sc_previous_sum_ips + 
                       #previous_agg1 + 
                       previous_agg2 +
                       previous_peak_doy2 +
                       #previous_peak_diff1  + 
                       previous_peak_diff2 +
                       poly(population_growth2,2)  + 
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")

m_moran20 <- glmmTMB(Morans_I ~ 
                       #spei3 +
                       #previous_spei3 + 
                       previous_spei3_2   + 
                       #veg_tmp + 
                       previous_veg_tmp + 
                       #sc_previous_sum_ips + 
                       #previous_agg1 + 
                       previous_agg2 +
                       previous_peak_doy2 +
                       #previous_peak_diff1  + 
                       previous_peak_diff2 +
                       population_growth2  + 
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")

# remove no significant ones
m_moran21 <- glmmTMB(Morans_I ~ 
                       #spei3 +
                       #previous_spei3 + 
                       previous_spei3_2   + 
                       #veg_tmp + 
                       previous_veg_tmp + 
                       #sc_previous_sum_ips + 
                       #previous_agg1 + 
                       previous_agg2 +
                       #previous_peak_doy2 +
                       #previous_peak_diff1  + 
                       #previous_peak_diff2 +
                       population_growth2  + 
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")
# remove spei
m_moran22 <- glmmTMB(Morans_I ~ 
                          previous_veg_tmp + 
                       previous_agg2 +
                       population_growth2  + 
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")

# addd year as random
m_moran23 <- glmmTMB(Morans_I ~ 
                       previous_veg_tmp + 
                       previous_agg2 +
                       population_growth2  + (1 | year) +
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     na.action = "na.fail")

# addd year as random
m_moran24 <- glmmTMB(Morans_I ~ 
                       previous_veg_tmp + 
                       previous_agg2 +
                       population_growth2  + #(1 | year) +
                       (1 | pairID),  data = dat_lag_scaled_complete,
                     family = tweedie,
                     na.action = "na.fail")



summary(m_moran24)
AIC(m_moran17,m_moran18,m_moran19,m_moran20,m_moran21,m_moran22,m_moran23)
plot(allEffects(m_moran24))
simulateResiduals(m_moran24, plot = T)


##### test only on positve MOra's I values: ---------------------------------
morans_sub <- dat_lag_scaled_complete %>% 
  dplyr::filter(Morans_I > 0)







m_moran_sub1 <- glmmTMB(Morans_I ~ 
                       previous_veg_tmp + 
                       previous_agg2 +
                       population_growth2  + #(1 | year) +
                       (1 | pairID),  data = morans_sub,
                     family = tweedie,
                     na.action = "na.fail")

# get only weather condistions - the best!!!
m_moran_sub2 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2 +
                          previous_veg_tmp + 
                          previous_agg2 +
                          population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")

# add poly - does not improve model, keep simple
m_moran_sub3 <- glmmTMB(Morans_I ~ 
                          poly(previous_spei3_2,2) +
                          poly(previous_veg_tmp,2) + 
                          previous_agg2 +
                          population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")

# add interaction between spei and tmp
m_moran_sub4 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2*previous_veg_tmp + 
                          previous_agg2 +
                          population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")

m_moran_sub5 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2 +
                          previous_veg_tmp + 
                          previous_agg2 +
                          population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")

m_moran_sub6 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2 +
                          veg_tmp + 
                          previous_agg2 +
                          population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")


m_moran_sub7 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2 +
                          veg_tmp + 
                          previous_agg2 +
                         # population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")


m_moran_sub8 <- glmmTMB(Morans_I ~ 
                          previous_spei3_2*veg_tmp + 
                          previous_agg2 +
                          # population_growth2  + #(1 | year) +
                          (1 | pairID),  data = morans_sub,
                        family = tweedie,
                        na.action = "na.fail")



summary(m_moran_sub8)
summary(m_moran_sub5)

r2(m_moran_sub8)
r2(m_moran_sub5)
fin.m.moran <- m_moran_sub8
AIC(m_moran_sub2, m_moran_sub1, m_moran_sub3, m_moran_sub4, m_moran_sub5, m_moran_sub6, m_moran_sub7, m_moran_sub8)
plot(allEffects(m_moran_sub8))
simulateResiduals(m_moran_sub8, plot = T)



### RS:  Link beetle data with observed RS damage ------------------------------------------------------------------------


dat_lag_RS_compl <- dd_RS_scaled  %>% 
  dplyr::filter(Morans_I> 0) %>%  # keep only postive MOran's I values
  na.omit() #%>% 
 
# perform EDA: MOrans_log -----------------------------
#library(GGally)
ggpairs(dat_lag_RS_compl[, c("wind_beetle", 
                             "Morans_I", "previous_Moran", "previous_Moran2",
                             "Morans_I_log", "previous_Moran_log1", "previous_Moran_log2")])
# choose Morans_I  and previous_Moran2


ggpairs(dat_lag_RS_compl[, c("wind_beetle", 
                             "population_growth", "population_growth2")])
# does not matter: choose population_growth2, +,-

ggpairs(dat_lag_RS_compl[, c("wind_beetle", 
                             "sum_ips", "previous_sum_ips", "previous_sum_ips2")])
# sum_ips seems the best

ggpairs(dat_lag_RS_compl[, c("wind_beetle", 
                             "agg_doy", "previous_agg1", "previous_agg2")])

# does not matter, signs is changing: -, +, +

ggpairs(dat_lag_RS_compl[, c("wind_beetle", 
                             "peak_diff", "previous_peak_diff1", "previous_peak_diff2")])
# peak_diff, +,+,+

hist(dat_lag_RS_compl$wind_beetle)
mean(dat_lag_RS_compl$wind_beetle)
sd(dat_lag_RS_compl$wind_beetle)
var(dat_lag_RS_compl$wind_beetle)
density(dat_lag_RS_compl$cumSum_beetlerate )


dat_lag_RS_compl$cumSum_beetlerate_adj <- pmin(pmax(dat_lag_RS_compl$cumSum_beetlerate, 1e-4), 1 - 1e-4)

summary(dat_lag_RS_compl$cumSum_beetlerate_adj)

# which distributio to use? beta, poisson or negbin?
# use beta distribution - residuals are completely wrong! continues with negbin  
mRS01.tmb <-  glmmTMB(cumSum_beetlerate_adj ~ #previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran + (1|pairID/trapID),
                     #na.action = 'na.fail',
                     family = beta_family(link = "logit"),
                     data = dat_lag_RS_compl)

summary(mRS01.tmb)
simulateResiduals(mRS01.tmb, plot = T)


###### get RS with beetle population predictors: ----------------------------------
# add predictors one by one

hist(dat_lag_RS_compl$wind_beetle)

# decide between MOran from log or not:
mm1 <- glm.nb(wind_beetle ~ Morans_I, dat_lag_RS_compl )
mm1.log <- glm.nb(wind_beetle ~ Morans_I_log , dat_lag_RS_compl )

mm2 <- glm.nb(wind_beetle ~ previous_Moran2, dat_lag_RS_compl )
mm2.log <- glm.nb(wind_beetle ~ previous_Moran_log2 , dat_lag_RS_compl )

AIC(mm1, mm1.log, mm2, mm2.log)

simulateResiduals(mm2.log, plot = T)

mm3 <- glm.nb(wind_beetle ~ previous_Moran2 + previous_sum_ips, dat_lag_RS_compl)

AIC(mm3, mm2, mm4, mm5)

mm4 <- glm.nb(wind_beetle ~ previous_Moran2 + previous_sum_ips2 +
                previous_agg1, dat_lag_RS_compl)

mm5 <- glm.nb(wind_beetle ~ previous_Moran2 + previous_sum_ips2 +
                previous_agg2, dat_lag_RS_compl)

# logs are worse! ise Morans from teh simple counts. Continue wit lag2, as it is more meaningful

simple_model <- glm.nb(wind_beetle ~ sum_ips + 
                         #beetle_sum2yrs +
                         previous_veg_tmp +
                         previous_sum_ips +
                         previous_sum_ips2 +
                         veg_tmp + 
                         spei3 + 
                         previous_spei3 +
                         previous_spei3_2 +
                         
                         agg_doy + # correlated with sum_ips
                         peak_doy +  ## correlated with sum_ips
                          peak_diff + # correlated with sum_ips and its lags
                         
                         previous_agg1 + 
                         previous_agg2 +
                         previous_peak_diff1 + 
                         previous_peak_diff2 +
                         Morans_I + 
                         previous_Moran +
                         previous_Moran + 
                         # +    # correlated with sum_ips and its lags
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





mRS.nb1<-  glm.nb(wind_beetle ~ #previous_sum_ips + 
                 beetle_sum2yrs,# +
                 #previous_veg_tmp +
                 #previous_spei3_2 +
                 #population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)
# mRS1 has perfect residuals!!

mRS.nb2<-  glm.nb(wind_beetle ~ previous_sum_ips,# + 
                # beetle_sum2yrs,# +
               #previous_veg_tmp +
               #previous_spei3_2 +
               #population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)



mRS.nb3<-  glm.nb(wind_beetle ~ #previous_sum_ips# + 
               # beetle_sum2yrs,# +
               #previous_veg_tmp +
               #previous_spei3_2 +
               population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)


mRS.nb4<-  glm.nb(wind_beetle ~ previous_sum_ips + 
                # beetle_sum2yrs +
                 #previous_veg_tmp +
                 #previous_spei3_2 +
                 population_growth2,
               na.action = 'na.fail',
               data = dat_lag_RS_compl)


mRS.nb5<-  glm.nb(wind_beetle ~ previous_sum_ips + 
                    population_growth2 +
                    previous_agg1 + 
                    previous_Moran2,
                  na.action = 'na.fail',
                  data = dat_lag_RS_compl)

AIC(mRS.nb5,mRS.tmb5)


# use mixed effects, as I have paired traps ----------------------------------

# start with one variable at time
dat_lag_RS_compl %>% 
  ggplot(aes(x = previous_agg2,
             y = wind_beetle)) + 
  geom_smooth()


mRS.tmb0<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                  +Morans_I, #+ (1|pairID),
                na.action = 'na.fail',
                family = nbinom1,
                data = dat_lag_RS_compl)
allEffects(mRS.tmb0)
AIC(mRS.tmb1, mRS.nb1, mRS.tmb2,mRS.nb4,mRS.nb3, mRS.nb2)
# add random effects
mRS.tmb1<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      +previous_Moran2+ (1|pairID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

# adding random effects inproved the model!

mRS.tmb2<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                 +previous_Moran2+ (1|pairID),
                na.action = 'na.fail',
                family = nbinom1,
                data = dat_lag_RS_compl)

mRS.tmb3<-  glmmTMB(wind_beetle ~ previous_sum_ips2 + 
                      +previous_Moran2+ (1|pairID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

# stay with last year beetle counts, add colonization date
mRS.tmb4<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      previous_agg1 +
                      +previous_Moran2+ (1|pairID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)


# add year as random
mRS.tmb5<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      previous_agg1 +
                      +previous_Moran2+ (1|pairID) +(1|year) ,
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

# keep adding population characteristics
mRS.tmb6<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      population_growth2 +
                      previous_agg2 +
                      previous_Moran2+ (1|pairID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

# specify nested structure
mRS.tmb7<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      population_growth2 +
                      previous_agg1 +
                      previous_Moran2+ (1|pairID/trapID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)


mRS.tmb8<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      population_growth2 +
                      previous_agg1 +
                      +previous_Moran2+  (1|pairID/trapID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

mRS.tmb9<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      population_growth2 +
                      previous_agg2 +
                      previous_Moran2 + (1|pairID/trapID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)


mRS.tmb10<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                      population_growth2 +
                      previous_agg2 +
                      previous_Moran2 + (1|pairID/trapID),
                    na.action = 'na.fail',
                    family = nbinom1,
                    data = dat_lag_RS_compl)

# remove population growth
mRS.tmb11<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                       #population_growth2 +
                       previous_agg2 +
                       previous_Moran + (1|pairID/trapID),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)

# keep population growth, remove beetle sum
mRS.tmb12<-  glmmTMB(wind_beetle ~ #previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran + (1|pairID/trapID),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)

# keep population growth, remove beetle sum
mRS.tmb13<-  glmmTMB(wind_beetle ~ #previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + (1|pairID/trapID),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)


# keep random only years
mRS.tmb15<-  glmmTMB(wind_beetle ~ #previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + 
                       previous_Moran + (1|year),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)

mRS.tmb16<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + 
                       previous_Moran + (1|year),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)


mRS.tmb17<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + 
                       #previous_Moran + 
                       (1|year),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)

mRS.tmb17.int<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                         previous_sum_ips*previous_agg2 +
                       previous_Moran2 + 
                       #previous_Moran + 
                       (1|year),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)
AIC(mRS.tmb17.int, mRS.tmb17)

mRS.tmb17.2 <-  glmmTMB(wind_beetle ~ #previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + 
                       #previous_Moran + 
                       (1|year),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)

AIC(mRS.tmb17.2, mRS.tmb17)
mRS.tmb18<-  glmmTMB(wind_beetle ~ previous_sum_ips + 
                       population_growth2 +
                       previous_agg2 +
                       previous_Moran2 + 
                       #previous_Moran + 
                       (1|year) +
                       (1|pairID/trapID),
                     na.action = 'na.fail',
                     family = nbinom1,
                     data = dat_lag_RS_compl)





AIC(mRS.tmb1, mRS.tmb2,mRS.tmb3,mRS.tmb4,mRS.tmb5, 
    mRS.tmb6, mRS.tmb7, mRS.tmb8, mRS.tmb9,mRS.tmb10,mRS.tmb11, mRS.tmb12, mRS.tmb13,
    mRS.tmb15,mRS.tmb16,mRS.tmb17,mRS.tmb18)

AICc(mRS1, mRS2, mRS6)
simulateResiduals(mRS.tmb17, plot = T)
summary(mRS.tmb17) # final one!! just include year as random effect! not trap pairs
# as here, we d nt account explicitely for weather condistions, that can explain the variation between years
fin.m.RS <- mRS.tmb17

(vif_values <- performance::check_collinearity(mRS.tmb17))
r2(mRS.tmb17)
summary(mRS.tmb17)
simulateResiduals(mRS.tmb4, plot = T)








######  RS explore both climatic and beetles  variables -----------------------------------------------------

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

m9.int <- glmmTMB(sum_ips ~ veg_tmp + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 + (1 | pairID),
              family = nbinom2,
              data = dd)

# remove year as random
m9.int2 <- glmmTMB(sum_ips ~ veg_tmp + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 +  (1 | pairID/trapID),
                  family = nbinom2,
                  data = dd)

# simplify random effect on trap level
m9.int3 <- glmmTMB(sum_ips ~ veg_tmp + previous_spei3_2 + 
                     veg_tmp*previous_spei3_2 +  (1 | pairID),
                   family = nbinom2,
                   data = dd)

# simplify random effect on trap level
m9.tw <- glmmTMB(sum_ips ~ veg_tmp + previous_spei3_2 + 
                     veg_tmp*previous_spei3_2 +  (1 | pairID),
                   family =  tweedie() ,
                   data = dd)

# remove interaction, as it is not significant
m9.tw1 <- glmmTMB(sum_ips ~ veg_tmp + previous_spei3_2 + 
                   (1 | pairID),
                 family =  tweedie() ,
                 data = dd)

# adding interaction improves model: keepit!
m9.tw2 <- glmmTMB(sum_ips ~ poly(veg_tmp,2) + previous_spei3_2 + 
                   veg_tmp*previous_spei3_2 +  (1 | pairID),
                 family =  tweedie() ,
                 data = dd)

# improve random structure
m9.tw3 <- glmmTMB(sum_ips ~ poly(veg_tmp,2) + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 +  (1 | pairID/trapID),
                  family =  tweedie() ,
                  data = dd)
# convergence
m9.tw4 <- glmmTMB(sum_ips ~ poly(veg_tmp,2) + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 +  (1 | pairID/trapID) + (1|year),
                  family =  tweedie() ,
                  data = dd)
# simplify random structure - convergence
m9.tw5 <- glmmTMB(sum_ips ~ poly(veg_tmp,2) + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 +  (1 | pairID) + (1|year),
                  family =  tweedie() ,
                  data = dd)

# simplify random structure: keep only year
m9.tw6 <- glmmTMB(sum_ips ~ poly(veg_tmp,2) + previous_spei3_2 + 
                    veg_tmp*previous_spei3_2 +(1|year),
                  family =  tweedie() ,
                  data = dd)




plot(allEffects(m9.int))
summary(m9.int)
AICc(m9.int, m9,m9.int2,m9.int3, m9.tw, m9.tw1, m9.tw2, m9.tw3, m9.tw4,m9.tw5, m9.tw6)
simulateResiduals(m9.int, plot = T)
simulateResiduals(m9.tw4, plot = T)

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


m.agg12.3.int <- glmmTMB(tr_agg_doy ~ veg_tmp + previous_spei3_2 + 
                           veg_tmp*previous_spei3_2 + (1| pairID),
                     family = beta_family(link = "logit"),
                     data = dd)

AIC(m.agg12.3.int, m.agg12.3)
summary(m.agg12.3)
summary(m.agg12.3.int)
plot(allEffects(m.agg12.3.int))
plot(allEffects(m.agg12.3))


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

m.peak9.spei3 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_spei3_2  + (1 | pairID),
                   family = beta_family(link = "logit"),
                   dd)

m.peak9.spei3_lag1 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_spei3  + (1 | pairID),
                         family = beta_family(link = "logit"),
                         dd)


AIC(m.peak9, m.peak9.spei3, m.peak9.spei3_lag1)

AIC(m.peak4, m.peak5, m.peak6, m.peak7, m.peak8, m.peak9)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.peak9.spei3, 
                                              plot = T)

fin.m.peak <- m.peak9.spei3

windows()
plot(allEffects(m.peak9.spei3 ))
r2(m.peak9.spei3)

summary(m.peak9.spei3)

acf(residuals(m.peak9.spei3))
acf(residuals(m.peak9))

(dw_result <- dwtest(residuals(m.peak9.spei3) ~ fitted(m.peak9.spei3)))


# Chek for multicollinearity
library(performance)

# Assuming your model is an 'lme4' or 'glmmTMB' object
vif_values <- performance::check_collinearity(m.peak9.spei3)
vif_values


# 
# summary(m.peak9.spei3)  # is teh winner!! 
# 
# > summary(m.peak9.spei3)  # is teh winner!! 
# Family: beta  ( logit )
# Formula:          tr_peak_doy ~ veg_tmp + previous_spei3_2 + (1 | pairID)
# Data: dd
# 
# AIC      BIC   logLik deviance df.resid 
# -1694.4  -1669.3    852.2  -1704.4     1101 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# pairID (Intercept) 0.005386 0.07339 
# Number of obs: 1106, groups:  pairID, 79
# 
# Dispersion parameter for beta family (): 18.9 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.22796    0.01590 -14.336  < 2e-16 ***
#   veg_tmp          -0.10535    0.01531  -6.879 6.01e-12 ***
#   previous_spei3_2  0.08414    0.01410   5.967 2.42e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


cor(dd$previous_spei6,dd$veg_tmp )
plot(dd$previous_spei3_2,dd$veg_tmp )
cor(dd$previous_spei3_2,dd$veg_tmp , method = 'spearman')
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

m.peak.diff2.6.poly <- glmmTMB(peak_diff ~ poly(veg_tmp,2)  + poly(previous_spei3_2,2)  +(1 | pairID),
                          family = nbinom2,
                          data = dd)


m.peak.diff2.6.int <- glmmTMB(peak_diff ~ veg_tmp  + previous_spei3_2  +
                                veg_tmp*previous_spei3_2  + (1 | pairID),
                          family = nbinom2,
                          data = dd)





AIC(m.peak.diff2.6, m.peak.diff2.6.int, m.peak.diff2.6.poly2)
summary(m.peak.diff2.6.int)
summary(m.peak.diff2.6)

plot(allEffects(m.peak.diff2.6.int))
allEffects(m.peak.diff2.6)

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



### effect plot peak dif ----------------------------------------------------
# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak.diff, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)


p1.peak.diff <- create_effect_plot(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Peak difference DOY", y_lim = c(220,650), show_y_axis = TRUE)
p2.peak.diff <-create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Peak difference DOY", y_lim = c(220,650), show_y_axis = TRUE)

### effect plots Sum IPS ----------------------------------------------------
p.effect.peak.diff <- ggarrange(p1.peak.diff,p2.peak.diff, ncol = 1, align = "hv")

p.effect.peak.diff


### Moran's I --------------------------------------------------------------------
# Generate predictions
pred <- ggpredict(fin.m.moran, terms = c("previous_spei3_2", "veg_tmp"))

ggplot(pred, aes(x = x, y = predicted, group = group)) +
  
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1.2) +
  #scale_color_manual(values = c('low' = "#00BFC4", 'mean' = "#F8766D", 'high' = "#A2C8EC")) +
  #scale_fill_manual(values = c('low' = "#00BFC4", 'mean' = "#F8766D", 'high' = "#A2C8EC")) +
  #scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  
  labs(x = "Previous SPEI-3", y = "Predicted Moran's I") +
  # facet_wrap(~group) + #, scales = "free"
  theme_classic()




### all Effect plots --------------------------------------------------------
windows(9,5)
ggarrange(p.effect.counts, p.effect.agg, p.effect.peak, p.effect.peak.diff, ncol=4, nrow = 1 , align = 'hv', font.label = list(size = 10, color = "black", face = "bold", family = NULL),
          labels = c( "[a]","[b]","[c]","[d]"))


# get characteristics
r2(fin.m.counts )
r2(fin.m.agg)
r2(fin.m.peak)
r2(fin.m.peak.diff)
r2(fin.m.moran)
r2(fin.m.RS)



summary(fin.m.counts )
summary(fin.m.agg)
summary(fin.m.peak)
summary(fin.m.peak.diff)
summary(fin.m.moran)
summary(fin.m.RS)



# print models outputs:



# export all models
sjPlot::tab_model(fin.m.counts,    file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg,       file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak,      file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff, file = "outTable/model_peak_diff.doc")
sjPlot::tab_model(fin.m.moran,     file = "outTable/model_moran.doc")
sjPlot::tab_model(fin.m.RS,        file = "outTable/model_RS.doc")





# Effect plots RS ----------------------------------------------------------

fin.m.RS <- mRS.tmb4

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.RS, terms = "previous_sum_ips [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.RS, terms = "previous_agg2   [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.RS, terms = "previous_Moran2   [all]", allow.new.levels = TRUE)
#p4 <- ggpredict(fin.m.RS, terms = "population_growth2 [all]", allow.new.levels = TRUE)


p1.RS <- create_effect_plot(p1, line_color = "red", 
                            x_title = "Beetle population level (lag1) [z-score]", 
                            y_title = "Tree mortality\n[# pixels]", y_lim = c(0,65), 
                            show_y_axis = TRUE)
p2.RS <- create_effect_plot(p2, line_color = "blue", 
                            x_title = "Colonization DOY (lag2) [z-score]", 
                            y_title = "Tree mortality\n[# pixels]", y_lim = c(0,65), 
                            show_y_axis = TRUE)
p3.RS <- create_effect_plot(p3, line_color = "yellow", 
                            x_title = "Local Moran's I (lag2) [z-score]", 
                            y_title = "Tree mortality\n[# pixels]", y_lim = c(0,80), show_y_axis = TRUE)
p4.RS <- create_effect_plot(p4, line_color = "grey", 
                            x_title = "Popul. growth (lag2) [z-score]", 
                            y_title = "Tree mortality\n[# pixels]", 
                            y_lim = c(0,80), show_y_axis = TRUE)



# effect plots Sum IPS ----------------------------------------------------
p.effect.RS <- ggarrange(p1.RS,p2.RS,p3.RS,
                         #p4.RS, 
                         ncol = 2, nrow = 2, align = 'hv')

windows(7,7)
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


# effect plot MOran's I -------------------------------


# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.moran , terms = "previous_veg_tmp  [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.moran, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.moran, terms = "previous_agg2  [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m.moran, terms = "population_growth2   [all]", allow.new.levels = TRUE) 


p1.counts <- create_effect_plot(p1, line_color = "red", x_title = "Temperature (lag2) [z-score]", y_title = "Local Moran's I", show_y_axis = TRUE, y_lim = c(0,1))
p2.counts <- create_effect_plot(p2, line_color = "blue", x_title = "SPEI (lag2) [z-score]", y_title = "Local Moran's I", show_y_axis = TRUE, y_lim = c(0,1))
p3.counts <- create_effect_plot(p3, line_color = "orange", x_title = "Colonization  [DOY, z-score]", y_title = "Local Moran's I", show_y_axis = TRUE, y_lim = c(0,1))
p4.counts <- create_effect_plot(p2, line_color = "grey", x_title = "Population growth [z-score]", y_title = "Local Moran's I", show_y_axis = TRUE, y_lim = c(0,1))



# Arrange the plots side by side with the same size
p.effect.moran <- ggarrange(p1.counts, p2.counts,
                            p3.counts, p4.counts,
                            ncol = 2, nrow=2, align = "hv")
(p.effect.moran)




