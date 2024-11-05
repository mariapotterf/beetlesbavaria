


rm(list=ls()) 
gc()

# Intro -----------------------------------------------------------------------
### Read my paths -----------------------------------------------------------
source('myPaths.R')
source('my_functions.R')



### get libs ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(tidyr)
#library(rgdal)
library(ggplot2)
library(ggpubr)
#library(ggpmisc)  # add equation to plots smooth 


# Stats
library('here')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
#library('cowplot')
library('tidyr')
library("knitr")

library('readr')

library(ggeffects)
library(stringr)


library(performance)  # Chek for multicollinearity

library(MASS)
library(car)     # for VIF
library(glmmTMB) # many families

library(lme4) #pseudo-R2 and and mixed model functionality
library(MuMIn) #dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc) #pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa) #model diagnostics
library(effects) #what do my marginal effects look like?
library(performance) #binomial model diagnostics
library(emmeans) #post hoc for categorical predictors


library(mgcv)  # for gams
library(gratia)

# test for autocorrelation
library(lmtest)

# colors
library(RColorBrewer)

library(knitr)   # for table outputs


# plot labels: 
lab_popul_level       = "Population level [#]"
lab_colonization_time = "Aggregation timing [DOY]"
lab_peak_time         = "Peak swarming timing [DOY]"
lab_peak_growth       = "Peak swarming intensity [#]"




##### read data ----------------------------------------------------------------
load(file=   "outData/final_table.Rdata")
load(file =  "outData/buffers.Rdata")  # df_RS_out
load(file =  "outData/lisa.Rdata")     # read LISA Moran's I stats
load(file =  "outData/spatial.Rdata")  # read xy coordinates




# Spatial data: 
sort(unique(xy_sf_expand$falsto_name))


# get coordinates from sf object
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    year = xy_sf_expand$year,
                    trapID = xy_sf_expand$falsto_name)

xy_df <- distinct(xy_df)


# check data: length 1106 or more!
nrow(df_RS_out)
nrow(lisa_merged_df)
nrow(dat_fin)
nrow(xy_df)

# prepare data RS --------------------------------------------
df_RS_out <- df_RS_out %>% 
  dplyr::select(-c(globalid, id,  "all_dist_sum",     
                   "cum_removed",      "ref_all_dist",     "anom_all_dist" ,   "ref_wind_beetle", 
                   "anom_wind_beetle", "ref_harvest",      "anom_harvest")) %>% 
  dplyr::rename(trapID = falsto_name)

# change column names
lisa_merged_df <- lisa_merged_df %>% 
  dplyr::select(year,      falsto_name,          sum_beetle,     log_sum_beetle,      Morans_I,              clust,      Morans_I_log) %>% 
  dplyr::rename(trapID = falsto_name,
                sum_ips = sum_beetle) %>% 
  dplyr::select(c(year, trapID, sum_ips, 
                  Morans_I, # calculated from beetles sums 
                  Morans_I_log)) # calculated from log(beetles sum)

# add MOran's I values, transform the agg doy and peak doy between 0-1
dat_fin <-   dat_fin %>% 
  left_join(lisa_merged_df, 
            by = join_by(trapID, year, sum_ips)) %>%
  mutate(tr_agg_doy   = (agg_doy - doy.start) / (doy.end - doy.start),
         tr_peak_doy  = (peak_doy - doy.start) / (doy.end - doy.start))

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dat_fin$tr_agg_doy  <- pmin(pmax(dat_fin$tr_agg_doy, 1e-4),  1 - 1e-4)
dat_fin$tr_peak_doy <- pmin(pmax(dat_fin$tr_peak_doy, 1e-4), 1 - 1e-4)


dat_fin <-   dat_fin %>% 
  left_join(df_RS_out, 
            by = join_by(trapID, year, sum_ips)) 
  
# add previous year lag
dat_fin <- dat_fin %>% 
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(peak_diff = as.integer(round(peak_diff))) %>% 
  mutate(sum_ips_lag1          = lag(sum_ips, n = 1, default = NA)) %>%  # lag population growth by one more year
 ungroup(.) 

# Merge beetle dynamics inicators with XYs
dat_fin <- 
  dat_fin %>% 
  left_join(xy_df, by = c("trapID", 'year')) %>% 
  # remove years woith NAs:
  dplyr::filter(year %in% 2012:2021) 

nrow(dat_fin)





# Specify predictors and lags ----------------------------
# dopes not remove the number of larg = annd number of observations constant
# List of predictors
predictors_spei       <- c("spei1", "spei3","spei6", "spei12", "spei24")
predictors_beetle_dyn <- c("spring_tmp", "veg_tmp",  "spei3", "spei12")
predictors_Morans     <- c("spring_tmp", "veg_tmp",  "spei3", "spei12", "spei24", "sum_ips", "peak_doy", "peak_diff", "agg_doy")
lags <- 0:3
  

# dat_fin: calculate climate lags --------------
dat_spei_lags <-  dat_fin %>% 
  dplyr::select(c(year, pairID, trapID, tmp, 
                  spei,  # spei3
                  spei12,
                  tmp_z, 
                  sum_ips, 
                  #sum_ips_lag1, 
                  tr_agg_doy, tr_peak_doy,peak_diff,
                  agg_doy, peak_doy, x, y, Morans_I_log )) %>% 
  group_by(trapID) %>%
  mutate(trapID = as.factor(trapID),
         pairID = as.factor(pairID),
         year_fact = as.factor(year)) %>% 
  arrange(year, .by_group = TRUE) %>%
  mutate(#sum_ips_lag1 = lag(sum_ips_log, n = 1, default = NA),
         tmp_lag1 = lag(tmp, n = 1, default = NA),
         tmp_lag2 = lag(tmp, n = 2, default = NA),
         tmp_lag3 = lag(tmp, n = 3, default = NA),
         tmp_z_lag1 = lag(tmp_z, n = 1, default = NA),
         tmp_z_lag2 = lag(tmp_z, n = 2, default = NA),
         tmp_z_lag3 = lag(tmp_z, n = 3, default = NA),
         spei_lag1 = lag(spei, n = 1, default = NA),
         spei_lag2 = lag(spei, n = 2, default = NA),
         spei_lag3 = lag(spei, n = 3, default = NA),
         spei12_lag1 = lag(spei12, n = 1, default = NA),
         spei12_lag2 = lag(spei12, n = 2, default = NA),
         spei12_lag3 = lag(spei12, n = 3, default = NA),
          ) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup(.) %>% 
  mutate(log_sum_ips = log(sum_ips+1),
         log_peak_diff = log(peak_diff + 1))

nrow(dat_spei_lags)


# list predictors to test
selected_predictors <- c('sum_ips', 'sum_ips_lag1','sum_ips_lag2',
                         #'log_sum_ips', 'log_sum_ips_lag1','sum_ips_lag2',
                         'tr_agg_doy'   , 'tr_agg_doy_lag1', 'tr_agg_doy_lag2' ,
                         'tr_peak_doy', 'tr_peak_doy_lag1','tr_peak_doy_lag2',
                         'peak_diff', 'peak_diff_lag1',  'peak_diff_lag2',
                         'spei',  'spei_lag1','spei_lag2',
                         'tmp_z',  'tmp_z_lag1', 'tmp_z_lag2'  
) 



# check correlation between tmp and spei ---------------------------------------
# Calculate pairwise correlations for numeric columns
numeric_cols <- dat_spei_lags %>%
  dplyr::select_if(is.numeric)  %>% # Select only numeric columns
  dplyr::select(tmp, tmp_lag1,tmp_lag2,tmp_lag3, spei,spei_lag1, spei_lag2, spei_lag3)


correlations <- cor(numeric_cols, use = "complete.obs")  # Calculate correlations, excluding missing values

# Display the correlation matrix
print(correlations)


dat_spei_lags %>% 
  dplyr::filter(trapID == "Anzinger_Forst_1") %>% 
  dplyr::select(trapID, year, tmp_z, tmp_z_lag1,tmp_z_lag2 ) # sum_ips, sum_ips_lag1

# check lags if correct. YES!

# Outliers handling: is it erroneours?? ------------------------------------------------------------------------------

# outliers can be cause by probelsm in the trap (low pheromones,,..) or by external effect - local windthrow
# investigate outliers for 4 dependent variables
# current observations: igf I have outliers present, i need higher k to converge the model
# also, maybe i can olny remove erroneous trap, not the whole trap pair

# Define the columns to check for outliers
outlier_columns <- c("log_sum_ips", "tr_agg_doy", "tr_peak_doy", "log_peak_diff")

# Function to detect outliers based on IQR
detect_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Filter the rows where the column value is an outlier
  outlier_data <- data %>%
    filter(!is.na(data[[column]]) & (data[[column]] < lower_bound | data[[column]] > upper_bound))
  
  if (nrow(outlier_data) == 0) {
    return(tibble(year = integer(), pairID = factor(), trapID = factor(), 
                  outlier_column = character(), outlier_value = numeric()))
  }
  
  outlier_data %>%
    mutate(outlier_column = column, outlier_value = .data[[column]]) %>%
    dplyr::select(year, pairID, trapID, outlier_column, outlier_value)
}

# Apply the outlier detection function to each column and combine the results
df_outliers <- lapply(outlier_columns, function(col) detect_outliers(dat_spei_lags, col)) %>%
  bind_rows()

# View the outliers
print(df_outliers)

# Count unique trapIDs and keep them in long format
outlier_details_long <- df_outliers %>%
  group_by(outlier_column, trapID, pairID) %>%
  summarise(outlier_count = n(), .groups = "drop") %>%
  group_by(outlier_column) %>%
  summarise(trap_count = n(), trapID = list(trapID), pairID = list(pairID)) %>%
  unnest_longer(c(trapID, pairID))

# View the detailed counts and individual trapIDs
print(outlier_details_long)

# Extract distinct pairIDs for each outlier column
pairID_outliers <- outliers %>%
  group_by(outlier_column, pairID) %>%
  summarise(.groups = "drop") %>%
  distinct()

# View the result
print(pairID_outliers)


# remove too low values: peiting and Eschenbach_idOPf: identify loutlier from log_sum_ips
df_outliers <- dat_spei_lags %>% 
  group_by(year) %>% 
  mutate(sum_ips_log = log(sum_ips)) %>% 
  mutate(is_outlier = ifelse(sum_ips_log > quantile(sum_ips_log, 0.75,na.rm=T) + 1.5 * IQR(sum_ips_log,na.rm=T) |
                               sum_ips_log < quantile(sum_ips_log, 0.25,na.rm=T) - 1.5 * IQR(sum_ips_log,na.rm=T), TRUE, FALSE)) # %>% 

# count how many outliers i have?
df_outliers %>% 
  dplyr::filter(is_outlier == TRUE) %>% 
  ungroup() %>% 
  summarise(n = n()) #%>%

# 73 are outliers, from 1106, which is 6.6%

# check which ones are outliers and how often?
df_outliers %>% 
  dplyr::filter(is_outlier == TRUE) %>% 
  ungroup() %>% 
  distinct(pairID)


pair_outliers <- df_outliers %>%
  dplyr::filter(is_outlier) %>%
  group_by(trapID, pairID) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  dplyr::filter(n > 3) %>% 
  pull(pairID)

unique(pair_outliers)

# if outlier > 4 time, remove the whole pair; keep the correct traps
# otherwise i need to remove 35 pairs, which is too much


# export table to merge it with the spatial data: for variograms:
dat_dynamics <- dat_fin %>% 
  dplyr::select(c(year, trapID, pairID, spring_tmp, veg_tmp, veg_prcp,   
                  #spei3, 
                  spei, # spei for veg season
                  tmp_z, prcp_z,# spei_z, # z score fr veg season, sd from teh mean reference conditions
                sum_ips,  peak_doy, peak_diff, tr_agg_doy, tr_peak_doy, agg_doy, peak_doy,
                x,y
                ))
fwrite(dat_dynamics, 'outTable/beetle_dynamic_indicators.csv')




hist(dat_dynamics$veg_tmp)
hist(dat_dynamics$spei)


# get summary statistics for mean and sd for weather condistions: --------------------------
desc_stat_climate_summary <- dat_dynamics %>%
  mutate(class = ifelse(year %in% c(2018), "2018", "Other Years")) %>%
  group_by(class) %>%
  summarise(
    mean_veg_tmp = mean(veg_tmp, na.rm = TRUE),
    med_veg_tmp = median(veg_tmp, na.rm = TRUE),
    sd_veg_tmp = sd(veg_tmp, na.rm = TRUE),
    mean_spei = mean(spei, na.rm = TRUE),
    med_spei_tmp = median(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE),
    mean_veg_prcp = mean(veg_prcp, na.rm = TRUE),
    med_veg_prcp = median(veg_prcp, na.rm = TRUE),
    sd_veg_prcp = sd(veg_prcp, na.rm = TRUE)
  )

(desc_stat_climate_summary)

desc_stat_climate_summary_year <- dat_dynamics %>%
  group_by(year) %>%
  summarise(
    mean_veg_tmp = mean(veg_tmp, na.rm = TRUE),
    med_veg_tmp = median(veg_tmp, na.rm = TRUE),
    sd_veg_tmp = sd(veg_tmp, na.rm = TRUE),
    mean_spei = mean(spei, na.rm = TRUE),
    med_spei_tmp = median(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE),
    mean_veg_prcp = mean(veg_prcp, na.rm = TRUE),
    med_veg_prcp = median(veg_prcp, na.rm = TRUE),
    sd_veg_prcp = sd(veg_prcp, na.rm = TRUE)
  )

(desc_stat_climate_summary_year)


# make a plot for each variable per year 
# Melt the data to long format for easier plotting with facets
# Rename columns for clarity
dat_dynamics2 <- dat_dynamics %>%
  dplyr::rename(temperature = veg_tmp, precipitation = veg_prcp, 
                SPEI = spei)

# Prepare the data for plotting by pivoting it longer
plot_data <- dat_dynamics2 %>%
  dplyr::select(year, temperature, precipitation, SPEI) %>%
  pivot_longer(cols = c(temperature, precipitation, SPEI), names_to = "variable", values_to = "value")

# Set the order of facets
plot_data$variable <- factor(plot_data$variable, levels = c("temperature", "precipitation", "SPEI"))

# Calculate the mean value for each variable for horizontal lines
mean_values <- plot_data %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

# Create the plot with violin plots, jittered points, mean with error bars, and horizontal mean lines
weather_summary_plot <- ggplot(plot_data, aes(x = as.factor(year), 
                                              y = value#, 
                                              #color = variable, 
                                              #fill = variable
                                              )) +
  geom_violin(alpha = 0.5, fill = 'grey60', color = 'grey60') +  # Violin plot to show distribution
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.5, color = 'grey20') +  # Jitter points for individual data
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", size = 0.3,
               position = position_dodge(width = 0.9), color = "red") +  # Mean and standard error bars
  facet_wrap(variable ~ ., scales = "free_y") +  # Facet by variable
  labs(x = "Year", y = "Value", title = "") +  # Labels and title
  theme_classic() +  # Use a classic theme for better visualization
  theme(
    legend.position = "none",  # Remove legend for cleaner appearance
    text = element_text(size = 8),  # Set all font sizes to 8 points
    axis.title = element_text(size = 8),  # Axis title font size
    axis.text = element_text(size = 8),  # Axis text font size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Axis text font size
    strip.text = element_text(size = 8)  # Facet strip text font size
  ) +
  # Add horizontal dashed line for the mean of each facet
  geom_hline(data = mean_values, aes(yintercept = mean_value), linetype = "dashed", color = "grey30")

# Display the plot
print(weather_summary_plot)


ggsave(filename = 'outFigs/vegetation_weather_summary3.png', 
       plot = weather_summary_plot, width = 7, height = 3, dpi = 300, bg = 'white')




# get final tables fo the model -------------------------------------------------


# table for ips counts, peak diff
dat_fin_counts_m <- 
  dat_spei_lags %>%  
  dplyr::select( -c(agg_doy, tr_agg_doy)) %>% 
  na.omit()
 
nrow(dat_fin_counts_m)  # 1106

# create new table with only agg values
dat_fin_agg_m <- 
  dat_spei_lags %>% 
   ungroup() %>%
  na.omit()
nrow(dat_fin_agg_m)  # 1080


#### Calculate the correlation matrix for all predictors and their lags ---------------------------

predictor_vars <- dat_fin_counts_m[, selected_predictors  ] # unique(best_predictors$Predictor)
cor_matrix <- cor(predictor_vars, method = "spearman")
(cor_matrix)
corrplot::corrplot(cor_matrix)

# mask highly correlated values
# Mask high correlations (set them to NA)
cor_matrix_masked <- cor_matrix
cor_matrix_masked[cor_matrix > 0.4 | cor_matrix < -0.4] <- NA
corrplot::corrplot(cor_matrix_masked)







# Calculate average counts for each trap pair ------------------------------------
avg_data <- dat_spei_lags %>%
  group_by(pairID, year) %>%
  summarise(Morans_I_log  = mean(Morans_I_log , na.rm = TRUE),
            
            # dependents
            sum_ips     = mean(sum_ips, na.rm = TRUE),
            peak_diff   = mean(peak_diff, na.rm = TRUE),
            tr_agg_doy  = mean(tr_agg_doy, na.rm = TRUE),
            tr_peak_doy = mean(tr_peak_doy, na.rm = TRUE), 
            agg_doy     = mean(agg_doy, na.rm = TRUE),
            peak_doy    = mean(peak_doy, na.rm = TRUE),
            # SPEIS - spei3
            # spei_lag0        = mean(spei, na.rm = TRUE),
            # spei_lag1   = mean(spei_lag1, na.rm = TRUE), # for agg_doy
            # spei_lag2   = mean(spei_lag2, na.rm = TRUE),
            # # SPEI 12
            spei12_lag0        = mean(spei12, na.rm = TRUE),
            spei12_lag1   = mean(spei12_lag1, na.rm = TRUE), # for agg_doy
            spei12_lag2   = mean(spei12_lag2, na.rm = TRUE),
            # TMP
            tmp_lag0         = mean(tmp, na.rm = T),
            tmp_lag1    = mean(tmp_lag1, na.rm = TRUE),
            tmp_lag2    = mean(tmp_lag2, na.rm = TRUE),  # for agg_doy
            # TMP_z
            # tmp_z_lag0       = mean(tmp_z, na.rm = TRUE),
            # tmp_z_lag1  = mean(tmp_z_lag1, na.rm = TRUE),
            # tmp_z_lag2  = mean(tmp_z_lag2, na.rm = TRUE),  # for agg_doy
            # 
            x           = mean(x, na.rm = TRUE),
            y           = mean(y, na.rm = TRUE)
            ) %>%
  ungroup(.) %>% 
  na.omit() %>% 
  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% c('Eschenbach_idOPf', 'Peiting') ) # remove two traps with always having too low trap counts

fwrite(avg_data, 'outTable/fin_tab_avg.csv')




# gamme is better then quasi and betar!
# filter extreme values: occuring too late in a year ()
avg_data_agg <- avg_data %>% 
  dplyr::filter(tr_agg_doy < 0.54)

# remove outliers: Piding
avg_data_peak_diff <- avg_data %>% 
  dplyr::filter(!pairID %in% c('Piding', 'Gangkofen'))


###### SUM IPS TW  ------------------------------------

###  Use previous year ------------------------------------------------

avg_data_filt <- avg_data %>% 
  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% pair_outliers )


avg_data_lagged <- avg_data %>% 
  group_by(pairID) %>%
  arrange(year, .by_group = TRUE) %>%
  # mutate(peak_diff = as.integer(round(peak_diff))) %>% 
  mutate(sum_ips_lag1          = lag(sum_ips, n = 1, default = NA),  # laf previous year values to adress autocorrelation
         peak_diff_lag1        = lag(peak_diff , n = 1, default = NA),
         tr_agg_doy_lag1       = lag(tr_agg_doy , n = 1, default = NA),
         tr_peak_doy_lag1      = lag(tr_peak_doy , n = 1, default = NA)) #%>%  
  

nrow(avg_data_filt) # 518
nrow(avg_data) # 539
nrow(avg_data_lagged) # 539

# create extra table for sum_ips_lag1

avg_data_filt_lagged <- avg_data_filt %>% #  avg_data %>% # 
  group_by(pairID) %>%
  arrange(year, .by_group = TRUE) %>%
  # mutate(peak_diff = as.integer(round(peak_diff))) %>% 
  mutate(sum_ips_lag1          = lag(sum_ips, n = 1, default = NA),  # laf previous year values to adress autocorrelation
         peak_diff_lag1        = lag(peak_diff , n = 1, default = NA),
         tr_agg_doy_lag1       = lag(tr_agg_doy , n = 1, default = NA),
         tr_peak_doy_lag1      = lag(tr_peak_doy , n = 1, default = NA)) %>%  
 na.omit()


# Remove empty levels
avg_data_filt_lagged$f_year <- droplevels(avg_data_filt_lagged$f_year)


#avg_data_filt_lagged %>% 
  avg_data_lagged %>%
  dplyr::filter(pairID == 'Anzinger_Forst') %>% 
  dplyr::select(year, sum_ips, sum_ips_lag1, 
                tr_peak_doy, tr_peak_doy_lag1 )



# TEST : discuss with Rupert and Dominic about best predictors ------------------
# 
# increase the number of iterations to converge the model
control <- list(niterPQL = 50)

cor(avg_data$spei12_lag1, avg_data$tmp_lag0)

# makde models with diffr autocorrelation structur (direct vs values from previous years) and compare their performacce 
#### IPS_SUM ---------------------------------------------------------------------- 


##### use previous years values: ----------------------------------------------------

##### tmp ---------------------

m.previous.tmp0_spei12_1_te <- gamm(sum_ips ~ 
                                     # tmp_lag0 +
                                     s(tmp_lag0, k = 3) +
                                     s(spei12_lag1, k = 3) + 
                                     te(tmp_lag0, spei12_lag1, k = 5) + 
                                     s(sum_ips_lag1, k = 5) + # Added term for previous beetle counts
                                     s(x,y),
                                   data = avg_data_lagged, #avg_data_filt_lagged, 
                                   family = tw,
                                   control = control)
vis.gam(m.previous.tmp0_spei12_1_te$gam)
AIC(m.previous.tmp0_spei12_1_te, m.previous.tmp0_spei12_1_ti)
gam.check(m.previous.tmp0_spei12_1_te$gam)
k.check(m.previous.tmp0_spei12_1_te$gam)
summary(m.previous.tmp0_spei12_1_te$gam)
summary(m.previous.tmp0_spei12_1_ti$gam)

predict_data <- ggpredict(m.previous.tmp0_spei12_1_te$gam, terms = c('tmp_lag0'))
predict_data <- ggpredict(m.previous.tmp0_spei12_1_te$gam, terms = c('spei12_lag1'))
predict_data <- ggpredict(m.previous.tmp0_spei12_1_te$gam, terms = c('tmp_lag0', 'spei12_lag1'))
plot(predict_data)

#### previous peak_diff  ------------------------------------------------------

m.peak.diff.previous.tmp_0_spei12_1 <- gamm(peak_diff ~ 
                                                s(tmp_lag0, k = 3) +
                                                s(spei12_lag1, k = 3) + 
                                                te(tmp_lag0, spei12_lag1, k = 10) + 
                                                s(peak_diff_lag1 , k =5) + # Added term for previous beetle counts
                                                s(x,y) + #,
                                              s(pairID, bs = 're') +
                                              s(f_year, bs = 're'),
                                              data = avg_data_lagged, #avg_data_filt_lagged, 
                                              family = tw,
                                              control = control)

m<-m.peak.diff.previous.tmp_0_spei12_1$gam

appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('peak_diff_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1'))
plot(predict_data)


#### previuous AGG DOY ----------------------------------------------------------


###### agg dy previous - test which model alternative is more stable - the specific one ------
avg_data_agg_no_out <- avg_data_filt_lagged %>% 
  dplyr::filter(tr_agg_doy < 0.54)


avg_data_lagged_noNA <- avg_data_lagged[complete.cases(avg_data_lagged), ]

anyNA(avg_data_lagged_noNA)

m.agg.previous_tmp0_spei12_1 <- gamm(tr_agg_doy ~
                                               s(tmp_lag0, k = 3) +
                                               s(spei12_lag1, k = 3) +
                                               te(tmp_lag0, spei12_lag1, k = 5) +
                                               s(tr_agg_doy_lag1 , k = 5) +# Added term for previous year numbers to account ofr autocorrelation
                                               s(x,y) +  
                                               s(pairID, bs = 're') +
                                               s(f_year, bs = 're')
                                             ,
                                             data = avg_data_lagged_noNA, #avg_data_agg_no_out, #avg_data_filt_lagged, #, #,
                                             family = Gamma(link = "log"),
                                             control = control)





summary(m.agg.previous_tmp0_spei12_1$gam)
summary(m.agg.previous_tmp0_spei12_1_year_re$gam)
AIC(m.agg.previous_tmp0_spei12_1_year_re, 
    m.agg.previous_tmp0_spei12_1, 
    m.agg.previous_tmp_z0_spei12_1)

m<-m.agg.previous_tmp0_spei12_1$gam
k.check(m)

appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tr_agg_doy_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1'))
plot(predict_data)






###### Previous PEak DoY ----------------------------------------------------


avg_data_peak_no_out <- avg_data_filt_lagged %>% 
  dplyr::filter(tr_peak_doy < 0.72) %>% 
  dplyr::filter(tr_peak_doy_lag1 < 0.72)

# Plot the variables to look for outliers
plot(avg_data_peak_no_out$tmp_lag0, avg_data_peak_no_out$tr_peak_doy)
plot(avg_data_peak_no_out$spei12_lag1, avg_data_peak_no_out$tr_peak_doy)


m.peak.doy.previous_tmp0_spei12_1 <- gamm(tr_peak_doy ~ 
                                              s(tmp_lag0, k = 5) +
                                              s(spei12_lag1, k = 5) +
                                              te(tmp_lag0, spei12_lag1, k = 3) +
                                              s(tr_peak_doy_lag1 , k = 5) + # Added term for previous year numbers to account ofr autocorrelation
                                              s(x,y, bs = 'gp')   +
                                              s(pairID, bs = 're') +
                                            s(f_year, bs = 're')
                                            ,
                                            data = avg_data_lagged_noNA, # avg_data_peak_no_out, 
                                            family = Gamma(link = "log"),
                                            control = control)


m<-m.peak.doy.previous_tmp0_spei12_1$gam
k.check(m)

appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tr_peak_doy_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1'))
plot(predict_data)


# save final models : current year TMP and previous year SPEI12 ---------------------------------------
fin.m.counts.previous.tw    <- m.previous.tmp0_spei12_1_te$gam
fin.m.peak.diff.previous.tw <- m.peak.diff.previous.tmp_0_spei12_1$gam
fin.m.agg.doy.gamma         <- m.agg.previous_tmp0_spei12_1$gam
fin.m.peak.doy.gamma        <- m.peak.doy.previous_tmp0_spei12_1$gam





# check for spatial and tmp autocorrelaton ----------------------------------
# Extract residuals
# Extract residuals from the model
residuals <- resid(m.peak.diff.previous.tmp_z_0_spei12_1$lme, type = "normalized")


# Check for temporal autocorrelation
par(mfrow = c(1, 2))  # Arrange plots in one row, two columns
acf(residuals, main = "ACF of Residuals")  # Autocorrelation function
pacf(residuals, main = "PACF of Residuals")  # Partial autocorrelation function


# Create a spatial weights matrix using the coordinates
library(spdep)
# Test for spatial autocorrelation
model_residuals <- residuals(m.previous.tmp_z0_spei1$lme, type = "pearson")

avg_data_filt_2019 <- avg_data_filt %>% 
  dplyr::filter(year == 2019)

# Create coordinates matrix
coords <- cbind(avg_data_filt$x, avg_data_filt$y)
# Create a spatial weights matrix using a distance threshold
dists <- as.matrix(dist(coords))  # Calculate distances between coordinates
nb <- dnearneigh(coords, 0, 10000)  # You can adjust the distance threshold as needed

# If you prefer a more dynamic approach, consider using a distance threshold
# to select neighbors, or use a combination of methods.
listw <- nb2listw(nb, style = "W")
# Perform Moran's I test
moran_test <- moran.test(model_residuals, listw)
print(moran_test)

m <- m.previous.tmp_z0_spei1$gam
# Generate predicted values for temperature and SPEI
pred_temp <- ggpredict(m, terms = "tmp_z")
pred_spei <- ggpredict(m, terms = "spei_lag1")
pred_interaction <- ggpredict(m, terms = c("tmp_z", "spei_lag1[-1.5,-0.5,0]"))


# Plot for temperature predictor
p1 <- plot(pred_temp)
# Plot for SPEI predictor
p2 <- plot(pred_spei)

# Plot for interaction between temperature and SPEI
p3 <- plot(pred_interaction)

# Print the plots
print(p1)
print(p2)
print(p3)




##### TW PREV COUNTS --------------------------------------------------------------


# the best!
m.counts.tw <- gamm(sum_ips ~ #s(year, k = 6) +
                         s(tmp_z_lag1, k = 5) +
                         s(spei_lag2, k = 8) + 
                         te(tmp_z_lag1, spei_lag2, k = 7) + 
                         #s(x, y, bs = 'gp', k = 70) + 
                         s(pairID, bs = 're') +
                      s(f_year, bs = 're'),
                       data = avg_data_filt, 
                       family = tw,
                       correlation = corAR1(form = ~ year | pairID),
                    control = control)

gam.check(m.counts.tw$gam)
summary(m.counts.tw$gam)

# test different k 
m.counts.tw1 <- gamm(sum_ips ~ #s(year, k = 6) +
                      s(tmp_z_lag1, k = 5) +
                      s(spei_lag2, k = 8) + 
                      te(tmp_z_lag1, spei_lag2, k = 7) + 
                      #s(x, y, bs = 'gp', k = 70) + 
                      s(pairID, bs = 're') +
                      s(f_year, bs = 're'),
                    data = avg_data_filt, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID),
                    control = control)

# test different k 
m.counts.tw2 <- gamm(sum_ips ~ #s(year, k = 6) +
                       s(tmp_z_lag1, k = 5) +
                       s(spei_lag2, k = 8) + 
                       te(tmp_z_lag1, spei_lag2, k = 7) + 
                       #s(x, y, bs = 'gp', k = 70) + 
                       s(pairID, bs = 're') +
                       s(f_year, bs = 're', k = 5),
                     data = avg_data_filt, 
                     family = tw,
                     correlation = corAR1(form = ~ year | pairID),
                     control = control)

AIC(m.counts.tw2, m.counts.tw, m.counts.tw1)



summary(m.counts.tw$gam)
AIC(m.counts.tw, m.counts.tw.year)

m.counts.tw.year <- gamm(sum_ips ~ s(year, k = 6) +
                      s(tmp_z_lag1, k = 5) +
                      s(spei_lag2, k = 8) + 
                      te(tmp_z_lag1, spei_lag2, k = 20) + 
                     # s(x, y, bs = 'gp', k = 70) + 
                      s(pairID, bs = 're'),
                    data = avg_data_filt, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID),
                    control = control)



# m.counts.tw.test.f <- gamm(sum_ips ~ f_year +
#                            s(tmp_lag1, k = 5) +
#                            s(spei_lag2, k = 8) + 
#                            te(tmp_lag1, spei_lag2, k = 20) + 
#                            # s(x, y, bs = 'gp', k = 70) + 
#                            s(pairID, bs = 're'),
#                          data = avg_data_filt, 
#                          family = tw,
#                          correlation = corAR1(form = ~ year | pairID))

# add year as random
# m.counts.tw.test.f.rndm <- gamm(sum_ips ~ s(f_year,bs = 're')  +
#                              s(tmp_z_lag1, k = 5) +
#                              s(spei_lag2, k = 8) + 
#                              te(tmp_z_lag1, spei_lag2, k = 10) + 
#                              # s(x, y, bs = 'gp', k = 70) + 
#                              s(pairID, bs = 're'),
#                            data = avg_data_filt, 
#                            family = tw,
#                            correlation = corAR1(form = ~ year | pairID))



# check fr multicollinearity
library(car)
vif_model <- lm(sum_ips ~ tmp_z_lag1 + spei_lag2 + year + pairID, data = avg_data_filt)
vif(vif_model)






appraise(m.counts.tw$gam)
summary(m.counts.tw$gam)
k.check(m.counts.tw$gam)
#gam.check(m.counts.tw$gam)
plot(m.counts.tw$gam, page = 1)



# for testing:
appraise(m.counts.tw.test$gam)
summary(m.counts.tw.test$gam)
k.check(m.counts.tw.test$gam)
#gam.check(m.counts.tw$gam)
plot(m.counts.tw.test$gam, page = 1)


# for testing:
appraise(m.counts.tw.test.f$gam)
summary(m.counts.tw.test.f$gam)
k.check(m.counts.tw.test.f$gam)
#gam.check(m.counts.tw$gam)
plot(m.counts.tw.test.f$gam, page = 1)


# with random effects of years
appraise(m.counts.tw.test.f.rndm$gam)
summary(m.counts.tw.test.f.rndm$gam)
k.check(m.counts.tw.test.f.rndm$gam)
#gam.check(m.counts.tw$gam)
plot(m.counts.tw.test.f.rndm$gam, page = 1)
#anova(m.counts.tw.test.f.rndm$gam)




### TEST models with different temp autocorrelation and proper SPEI --------------------------------------------------------------

# increase the number of iterations to converge the model
control <- list(niterPQL = 50)


#### test effect of previuos year cunts ----------------------------------------------
# the best!
m.counts.previous <- gamm(sum_ips ~ 
                            s(tmp_z_lag1, k = 5) +
                            s(spei_lag2, k = 8) + 
                            te(tmp_z_lag1, spei_lag2, k = 7) + 
                            s(sum_ips_lag1, k = 5) + # Added term for previous beetle counts
                            s(pairID, bs = 're') +
                            s(f_year, bs = 're'),
                          data = avg_data_filt_lagged, 
                          family = tw,
                          control = control)


##### TW PREV PEAK_DIFF --------------------------------------------------------------
# for peak_diff
m.peak.diff.previous <- gamm(peak_diff ~ 
                               s(tmp_z_lag1, k = 5) +
                               s(spei_lag2, k = 8) + 
                               te(tmp_z_lag1, spei_lag2, k = 7) + 
                               s(peak_diff_lag1 , k = 5) + # Added term for previous beetle counts
                               s(pairID, bs = 're') +
                               s(f_year, bs = 're'),
                             data = avg_data_filt_lagged, 
                             family = tw,
                             control = control)

avg_data_agg_no_out <- avg_data_filt_lagged %>% 
  dplyr::filter(tr_agg_doy < 0.54)
 
###### agg dy previous - test which model alternative is more stable - the specific one ------
m.agg.previous <- gamm(tr_agg_doy ~
                         s(tmp_z_lag2, k = 4,bs = "tp") +
                         s(spei_lag1, k = 3) +
                         te(tmp_z_lag2, spei_lag1, k = 3,bs = "cs") +
                         s(tr_agg_doy_lag1 , k = 3)# + # Added term for previous year numbers to account ofr autocorrelation
                       #s(x,y, bs = 'gp')  
                       #s(pairID, bs = 're')+
                         #s(f_year, bs = 're')
                         ,
                       data = avg_data_agg_no_out, #avg_data_filt_lagged,
                       family = Gamma(link = "log"),
                       control = control)
# 
# with random effects of years
appraise(m.agg.previous$gam)
summary(m.agg.previous$gam)
k.check(m.agg.previous$gam)
gam.check(m.agg.previous$gam)
plot(m.agg.previous$gam, page = 1)
#anova(m.counts.tw.test.f.rndm$gam)





# keep he one with previous year values
# better! 
# m.agg.gamma <- gamm(tr_agg_doy ~ #s(year, k = 6) + 
#                       s(tmp_z_lag2, k = 3, bs = "cs") + # 
#                       s(spei_lag1, k = 3,bs = "cs") + 
#                       te(tmp_z_lag2, spei_lag1, k = 4, bs = "cs") + #
#                       #s(x, y, bs = 'gp', k = 30) + 
#                       s(pairID, bs = 're') +
#                       s(f_year, bs = 're'),
#                     data =  avg_data_agg , 
#                     family = Gamma(link = "log"),
#                     correlation = corAR1(form = ~ year | pairID))



# test for peak day ----------------------
# previous vs autocorrelation specified
# agg dy previous - test which model alternative is more stable - the specific one

# Calculate Q1, Q3, and IQR
Q1 <- quantile(avg_data_filt$tr_peak_doy, 0.25)
Q3 <- quantile(avg_data_filt$tr_peak_doy, 0.75)
IQR <- Q3 - Q1

# Define outlier thresholds
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Filter out the outliers
avg_data_filt_peak_doy <- avg_data_filt[avg_data_filt$tr_peak_doy >= lower_bound & avg_data_filt$tr_peak_doy <= upper_bound, ]

nrow(avg_data_filt_peak_doy) # 512
nrow(avg_data_filt)




# previous PEAK_DOY -----------------
avg_data_peak_no_out <- avg_data_filt_lagged %>% 
  dplyr::filter(tr_peak_doy < 0.72) %>% 
  dplyr::filter(tr_peak_doy_lag1 < 0.72)

# agg dy previous - test which model alternative is more stable - the specific one
m.peak.previous.gamm <- gamm(tr_peak_doy ~ 
                         s(tmp_z_lag2, k = 4,bs = "tp") +
                         s(spei_lag1, k = 4) +
                         te(tmp_z_lag2, spei_lag1, k = 7,bs = "tp") +
                         s(tr_peak_doy_lag1 , k = 3) + # Added term for previous year numbers to account ofr autocorrelation
                       #s(x,y, bs = 'gp')  
                       s(pairID, bs = 're')# +
                       #s(f_year, bs = 're')
                       ,
                       data = avg_data_peak_no_out, 
                       family = Gamma(link = "log"),
                       control = control)

summary(m.peak.previous.gamm$gam)
appraise(m.peak.previous.gamm$gam)
k.check(m.peak.previous.gamm$gam)
plot(m.peak.previous.gamm$gam, page = 1, shade = TRUE)


# final model: m.peak.previous.gamm

AIC(m.peak.previous.year, m.peak.previous)








####### Find predicors &  lags: MOrans  -------------

# Load the saved Rdata file
load("outData/lisa_avg.Rdata")         # LISA and tables from averaged values


lisa_merged_df_avg <- lisa_merged_df_avg %>% 
  mutate(pairID = factor(pairID)) %>% 
  group_by(pairID) %>% 
  # get lags of beetle indicators
  arrange(year, .by_group = TRUE) %>%
  mutate(Morans_I_log_lag1 = lag(Morans_I_log , n = 1, default = NA),
    sum_ips_lag1 = lag(sum_ips , n = 1, default = NA),
         sum_ips_lag2 = lag(sum_ips , n = 2, default = NA),
         peak_diff_lag1 = lag(peak_diff , n = 1, default = NA),
         peak_diff_lag2 = lag(peak_diff , n = 2, default = NA),
         tr_agg_doy_lag1 = lag(tr_agg_doy , n = 1, default = NA),
         tr_agg_doy_lag2 = lag(tr_agg_doy , n = 2, default = NA),
         tr_peak_doy_lag1 = lag(tr_peak_doy , n = 1, default = NA),
         tr_peak_doy_lag2 = lag(tr_peak_doy , n = 2, default = NA)) %>% 
  na.omit() %>% 
  # convert to log values for facter calculation
  mutate(log_sum_ips = log(sum_ips),
         log_sum_ips_lag1 = log(sum_ips_lag1),
         log_sum_ips_lag2 = log(sum_ips_lag2),
         log_peak_diff    = log(peak_diff),
         log_peak_diff_lag1 = log(peak_diff_lag1),
         log_peak_diff_lag2 = log(peak_diff_lag2)) %>% 
  mutate(f_year = factor(year)) %>% 
  dplyr::select(-geometry) %>% 
  as.data.frame()



dependent_moran <-  c("Morans_I_log")



plot(lisa_merged_df_avg$sum_ips, lisa_merged_df_avg$Morans_I_log )

lisa_merged_df_avg %>% 
  ggplot(aes(x = sum_ips, 
             y = Morans_I_log,
             color = year)) +
  geom_point() + 
  facet_wrap(.~year)

# test best families
m.gauss <- gam(Morans_I_log ~ s(sum_ips, k = 5),  
             #family = gaussian(),
             method = 'REML',  
             data = lisa_merged_df_avg,
             correlation = corAR1(form = ~ year | pairID))

plot(m.gauss)
summary(m.gauss)

# does not work with tw!
m.tw <- gam(Morans_I_log ~ s(sum_ips),  
                family = tw,
                method = 'REML',  
                data = lisa_merged_df_avg,
                correlation = corAR1(form = ~ year | pairID))

# ----------------------------------

# Initialize a data frame to store AIC values and deviance explained
model_metrics_moran <- data.frame(Predictor = character(), 
                                  Dependent = character(), 
                                  AIC = numeric(), 
                                  R_squared = numeric())



# Loop over each dependent variable
for (dep in dependent_moran) {
  #print(dep)
  # Loop over each predictor
  for (pred in selected_predictors) {
    print(pred)
    # Fit the model
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 10)"))
    print(formula)

    model <- gam(formula,  
                  family = gaussian(),
                  method = 'REML',  
                  data = lisa_merged_df_avg)
    
    # Extract model summary
    model_summary <- summary(model)  # model$gam
    
    # Store the AIC value and deviance explained
    model_metrics_moran <- rbind(model_metrics_moran, data.frame(Predictor = pred, 
                                                           Dependent = dep, 
                                                           AIC = AIC(model), 
                                                           R_squared = round(model_summary$r.sq*100,1)
    ))
  }
}

# View the AIC values and deviance explained

(model_metrics_moran)


####  Export as a nice table in word: ------------------------------------------
sjPlot::tab_df(model_metrics_moran,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_lag_predictors_moran.doc",
               digits = 1) 


#fin.m.moran

m.morans.previous <- gam(Morans_I_log ~ 
                            s(tmp_z, k = 8) +
                            s(spei_lag1, k = 6) + 
                            te(tmp_z, spei_lag1, k = 5) + 
                           s(sum_ips, k =5, bs = 'tp') + # Added term for  beetle counts
                          # s(sum_ips_lag1, k = 5) + # Added term for previous beetle counts
                            s(Morans_I_log_lag1, k = 5),# + # Added term for previous MOrans I 
                           # s(pairID, bs = 're') +
                           # s(f_year, bs = 're'),
                          data = lisa_merged_df_avg, 
                          family = gaussian())


summary(m.morans.previous)
k.check(m.morans.previous)
gam.check(m.morans.previous)
appraise(m.morans.previous)
draw(m.morans.previous)

fin.m.moran <- m.morans.previous 




# check boxplots ----------------------------------------------------------------
boxplot(avg_data_moran$sum_ips)
boxplot(avg_data_moran_sub$Morans_I_log)


# Assuming avg_data_moran is your data frame
# Select relevant columns excluding Morans_I_log
predictor_columns <- avg_data_moran_sub %>%
  ungroup(.) %>% 
  dplyr::select(-Morans_I_log,-pairID  ,-year) %>%
  dplyr::select(all_of(selected_predictors))


# Compute the correlation matrix
correlation_matrix <- cor(predictor_columns, use = "pairwise.complete.obs", method = "spearman")

# Print the correlation matrix
print(correlation_matrix)

corrplot::corrplot(correlation_matrix)






##### quick plotting -----------------------------------------------
fin.m <- m.agg.previous$gam # m2.agg.gamma$gam #m.counts.tw$gam # m.counts.tw.bam.year#m.counts.tw.test.f.rndm#  fin.m.agg#$gam

# check for tempoeal autocorrelation

# Extract residuals from the model
residuals <- resid(m.agg.previous$lme, type = "normalized")

# Plot ACF of the residuals
acf(residuals, main="ACF of Model Residuals")
pacf(residuals, main="ACF of Model Residuals")


# plot in easy way how tdoes the k value affect interaction
#p1 <- ggpredict(fin.m, terms = "year [all]", allow.new.levels = TRUE)
#p1 <- ggpredict(fin.m, terms = "peak_diff_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "tmp_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = "spei_lag1 [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m, terms = c("tmp_z_lag2", "spei_lag1 [-1, 0, 1]"), allow.new.levels = TRUE)

#p_df <- as.data.frame(p3)
# test simple plot:
plot1<-ggplot(p1, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot2<-ggplot(p2, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot3<-ggplot(p3, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot4<-ggplot(p4, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

ggarrange(#plot1,
          plot2,plot3,plot4, ncol = 2, nrow = 2)




# update& save the best models: --------------------------------------------------------

# # explore results for individual dependent variables: select the onse with the lowest AIC/with the all random effects
# fin.m.counts     <- m.counts.tw$gam #result$sum_ips$models$random_effect$gam  # lowest AIC
# fin.m.agg        <- m2.agg.gamma$gam #m1.gamma$gam#result$tr_agg_doy$models$random_effect$gam  # lowest AIC
# fin.m.peak       <- m.peak.gamma$gam #m.peak.diff.tw$gam #result$tr_peak_doy$models$random_effect$gam  # lowest AIC
# fin.m.peak.diff  <- m.peak.diff.tw$gam #result_peak_diff$peak_diff$models$random_effect$gam  # lowest AIC
# 
# 
# save(fin.m.counts.previous.tw,
#      fin.m.peak.diff.previous.tw,
#      fin.m.agg.doy.gamma,
#      fin.m.peak.doy.gamma,
#          avg_data,
#      file = "outData/fin_models.RData")









# read spatial ata to check where are the significnat locations
load("outData/spatial.Rdata")
dat_dynamics <- fread( 'outTable/beetle_dynamic_indicators.csv')

# Spatial data: 
sort(unique(xy_sf_expand$falsto_name))


# get coordinates from sf object
xy_df <- data.frame(x = sf::st_coordinates(xy_sf_expand)[,"X"],
                    y = sf::st_coordinates(xy_sf_expand)[,"Y"],
                    year = xy_sf_expand$year,
                    trapID = xy_sf_expand$falsto_name)

xy_df <- distinct(xy_df)
my_crs <- crs(xy_sf_expand)


# Merge beetle dynamics inicators with XYs
dat_dynamics_sf <- 
  dat_dynamics %>% 
  left_join(xy_df, by = c("trapID", 'year')) %>% 
  # remove years woith NAs:
  dplyr::filter(year %in% 2015:2021) %>% 
  mutate(log_sum_ips = log(sum_ips),
         log_peak_diff = log(peak_diff)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 3035)
#na.omit() # as agg_doy have some NAs

# select only signif ones
dat_dynamics_sf_sign <- dat_dynamics_sf %>% 
  dplyr::filter(pairID %in% significant_traps$pairID)

# Ensure CRS match
dat_dynamics_sf_sign <- st_transform(dat_dynamics_sf_sign, st_crs(dat_dynamics_sf))

v_all <- vect(dat_dynamics_sf)
v_sub <- vect(dat_dynamics_sf_sign)

# Get the extents of both layers
extent_all <- ext(v_all)
extent_sub <- ext(v_sub)

# Combine extents to get the overall extent
combined_extent <- union(extent_all, extent_sub)

# Define a yellow-to-red color palette
yellow_to_red_palette <- colorRampPalette(c("blue", "red"))


# Plot the layers with the combined extent
plot(v_all, ext = combined_extent, col = 'black', main = "Sign difference plots")
plot(v_sub['sum_ips'],  ext = combined_extent, add = TRUE, col = 'red')
# Plot the subset layer with 'log_sum_ips' values
plot(v_sub, 'log_sum_ips', ext = combined_extent, 
     col = colorRampPalette(c("blue", "red"))(10), add = T)


# plotting with sf values --------------------

plot(dat_dynamics_sf['log_sum_ips'], col = 'black')
plot(dat_dynamics_sf_sign['log_sum_ips'])

plot(dat_dynamics_sf_sign['log_sum_ips'])
# Add the points from the second dataset
plot(st_geometry(dat_dynamics_sf_sign), col = 'red', add = TRUE)

# Initial plot
plot(dat_dynamics_sf['log_sum_ips'], col = 'black')

# Add additional points
plot(st_geometry(dat_dynamics_sf_sign), col = 'red', add = TRUE)











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
                 sd_sum_ips     = sd(sum_ips, na.rm = T),
                 sd_agg_doy     = sd(agg_doy, na.rm = T),
                 sd_peak_doy    = sd(peak_doy, na.rm = T),
                 sd_peak_diff   = sd(peak_diff, na.rm = T)) %>% 
  mutate(Population_level       = stringr::str_glue("{round(mean_sum_ips,1)}±{round(sd_sum_ips,1)}"),
         Aggregation_timing       = stringr::str_glue("{round(mean_agg_doy,1)}±{round(sd_agg_doy,1)}"),
         Peak_swarming_timing     = stringr::str_glue("{round(mean_peak_doy,1)}±{round(sd_peak_doy,1)}"),
         Peak_swarming_intensity  = stringr::str_glue("{round(mean_peak_diff,1)}±{round(sd_peak_diff,1)}")) %>% 
  dplyr::select(year, Population_level, Aggregation_timing, Peak_swarming_timing,Peak_swarming_intensity) 


(means_dat_fin_year)


# Export as a nice table in word:
sjPlot::tab_df(means_dat_fin_year,
#               col.header = c(as.character(qntils), 'mean'),
               show.rownames = TRUE,
               file="outTable/summary_out_year.doc",
               digits = 0) 





#cor(avg_data_sum_ips$tmp_z, avg_data_sum_ips$sum_ips)
# Make table for response and predictor variables for models:





# Effect plots ------------------------------------------------------------


# define plotiing functions


plot_effect_interactions(p3, avg_data = avg_data_filt_lagged_plotting)


##### Spagetting plot: development over time -----------------------------------------


# get update table for the line an pint plotting
avg_data_sum_ips <- avg_data %>%
  mutate(sum_ips = sum_ips/1000,
         peak_diff = peak_diff/10)

# point table for potting 
avg_data_filt_lagged_plotting <- avg_data_filt_lagged %>% 
  mutate(sum_ips        = sum_ips/1000,
         sum_ips_lag1   = sum_ips_lag1/1000,
         peak_diff      = peak_diff/10,
         peak_diff_lag1 = peak_diff_lag1/10
  )





# Example usage:
p_spagett_ips       <- plot_data_with_average(avg_data_filt_lagged_plotting, "sum_ips", lab_popul_level,   
                                              my_title = paste("[a]", 'Population level', '\n[#*1000]'))
p_spagett_agg_doy   <- plot_data_with_average(avg_data_filt_lagged_plotting, "agg_doy", lab_colonization_time, 
                                              my_title = paste("[b]", 'Aggregation timing', '\n[DOY]'))
p_spagett_peak_doy  <- plot_data_with_average(avg_data_filt_lagged_plotting, "peak_doy", lab_peak_time, 
                                              my_title = paste("[c]", 'Peak swarming timing', '\n[DOY]'))
p_spagett_peak_diff <- plot_data_with_average(avg_data_filt_lagged_plotting, "peak_diff", lab_peak_growth, 
                                              my_title = paste("[d]", 'Peak swarming intensity', '\n[#*10]'))




#windows(7,6)
p_spagetti <- ggarrange(p_spagett_ips, 
                        p_spagett_agg_doy, 
                        p_spagett_peak_doy, 
                        p_spagett_peak_diff, ncol = 4,nrow = 1, align = 'hv')
p_spagetti
#ggsave(filename = 'outFigs/Fig1.png', plot = p_spagetti, width = 7, height = 6, dpi = 300, bg = 'white')






##### PLOT: IPS SUM counts -----------------------------------------------------------

temp_label <- "Temp. [z-score]" #expression(paste("Temperature [", degree, "C]", sep=""))
spei_label <- 'SPEI [dim.]'


#my_lab <-  expression(paste("Temp. [", "lag1 ", "z-score]"))
# Define the transformation function
adjust_predictions_counts <- function(df, divisor) {
  df <- as.data.frame(df)
  #df$x <- df$x  / divisor
  df$predicted <- df$predicted  / divisor
  df$conf.low <- df$conf.low / divisor
  df$conf.high <- df$conf.high / divisor
  return(df)
}


transform_predictions_DOY <- function(predictions, doy.start, doy.end) {
  scale_factor <- doy.end - doy.start
  predictions$predicted <- (predictions$predicted * scale_factor) + doy.start
  predictions$conf.low  <- (predictions$conf.low * scale_factor) + doy.start
  predictions$conf.high <- (predictions$conf.high * scale_factor) + doy.start
  return(predictions)
}


transform_single_prediction_DOY <- function(predicted_value, doy.start, doy.end) {
  scale_factor <- doy.end - doy.start
  transformed_value <- (predicted_value * scale_factor) + doy.start
  return(transformed_value)
}


# Transform a single value
transform_single_prediction_DOY(0.43, doy.start, doy.end)


# Assuming 'model' is your glm.nb model
summary(fin.m.counts.previous.tw)
p0 <- ggpredict(fin.m.counts.previous.tw, terms = "sum_ips_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.counts.previous.tw, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts.previous.tw, terms = "spei_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.counts.previous.tw, terms = c("tmp_z_lag1", "spei_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

divisor_population <- 1000
p0 <- adjust_predictions_counts(p0, divisor_population)
p0$x <- p0$x/1000
p1 <- adjust_predictions_counts(p1, divisor_population)
p2 <- adjust_predictions_counts(p2, divisor_population)
p3 <- adjust_predictions_counts(p3, divisor_population)



p0.count <- create_effect_previous_year(data = p0, 
                              avg_data = avg_data_filt_lagged_plotting,
                              x_col = "sum_ips_lag1",
                              y_col = "sum_ips",
                              line_color = "darkgreen", 
                              x_title = "Pop. level lag1 [*1000]", 
                              y_title = paste(lab_popul_level, '*1000'),# "Population level\n[# beetle*100]",
                              #my_title = paste("[a]", 'Population level', '\n[#*100]'), 
                              x_annotate = 50, lab_annotate = "***")

(p0.count)
p1.count <- 
  create_effect_plot(data = p1, 
                     avg_data = avg_data_filt_lagged_plotting, 
                     x_col = "tmp_z_lag1", 
                     y_col = "sum_ips", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag1 [z-score]")),
                     y_title = paste(lab_popul_level, '*1000'), 
                     x_annotate = 2, 
                     lab_annotate = "**")

(p1.count)
p2.count <- create_effect_plot(data = p2, avg_data = avg_data_filt_lagged_plotting, 
                               x_col = "spei_lag2", y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI. lag2 [dim.]")),#spei_label , 
                               y_title = paste(lab_popul_level, '*1000'), 
                               x_annotate = -0.5, 
                               lab_annotate = "***")
(p2.count)
p3.count <- plot_effect_interactions(data = p3, 
                                     avg_data = avg_data_filt_lagged_plotting, 
                                     x_col = "tmp_z_lag1", 
                                     y_col = "sum_ips", 
                                     temp_label = expression(paste("Temp. lag1 [z-score]")), 
                                     y_title = paste(lab_popul_level, '*1000'),
                                     x_annotate = 2,
                                     lab_annotate = "**") 

(p3.count)
ggarrange(p0.count,p1.count,p2.count, p3.count)




##### PLOT DOY aggregation ---------------------------------------------------------
summary(fin.m.agg.doy.gamma)
p0 <- ggpredict(fin.m.agg.doy.gamma, terms = "tr_agg_doy_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.agg.doy.gamma, terms = "tmp_z_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg.doy.gamma, terms = "spei_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.agg.doy.gamma, terms = c("tmp_z_lag2", "spei_lag1 [-1, 0, 1]"), allow.new.levels = TRUE)

# Apply transformation to each prediction data frame
p0 <- transform_predictions_DOY(p0, doy.start, doy.end)
p1 <- transform_predictions_DOY(p1, doy.start, doy.end)
p2 <- transform_predictions_DOY(p2, doy.start, doy.end)
p3 <- transform_predictions_DOY(p3, doy.start, doy.end)


# Calculate the scale factor
scale_factor <- doy.end - doy.start

# Apply the transformation directly to the column
avg_data_filt_lagged_plotting$tr_agg_doy_lag1_plot <- (avg_data_filt_lagged_plotting$tr_agg_doy_lag1 * scale_factor) + doy.start
avg_data_filt_lagged_plotting$tr_peak_doy_lag1_plot <- (avg_data_filt_lagged_plotting$tr_peak_doy_lag1 * scale_factor) + doy.start

p0$x <- (p0$x * scale_factor) + doy.start


p0.agg <- create_effect_previous_year(data = p0,
                               avg_data = avg_data_filt_lagged_plotting,
                               x_col = "tr_agg_doy_lag1_plot",
                               y_col = "agg_doy",
                               line_color = "darkgreen",
                               x_title = "Agg. timing lag1 [DOY]",
                               y_title = lab_colonization_time,
                               x_annotate = 175, lab_annotate = "***")

(p0.agg)
p1.agg <- 
  create_effect_plot(data = p1,  avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),  
                     x_col = "tmp_z_lag2", y_col = "agg_doy", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag2 [z-score]")) , 
                     y_title = lab_colonization_time, 
                      x_annotate = 2, 
                     lab_annotate = "***")

(p1.agg)
p2.agg <- create_effect_plot(data = p2, avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),    
                               x_col = "spei_lag1", y_col = "agg_doy", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI lag1 [dim.]")) , 
                               y_title = lab_colonization_time, 
                               x_annotate = -0.5, 
                               lab_annotate = "0.49")
(p2.agg)
p3.agg <- plot_effect_interactions(p3,
                                   avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),
                                   x_col = "tmp_z_lag2", 
                                   y_col = "agg_doy", 
                                     temp_label = expression(paste("Temp. lag2 [z-score]")), 
                                     y_title = lab_colonization_time,
                                     x_annotate = 2,
                                     lab_annotate = "***") 

(p3.agg)


#### PLOT: Peak DOY  ------------------------------------------------------------
summary(fin.m.peak.doy.gamma)

p0 <- ggpredict(fin.m.peak.doy.gamma, terms = "tr_peak_doy_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak.doy.gamma, terms = "tmp_z_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.doy.gamma, terms = "spei_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.doy.gamma, terms = c("tmp_z_lag2" ,"spei_lag1 [-1, 0, 1]"), allow.new.levels = TRUE)


# transform values back to DOY (from 0-1)
p0 <- transform_predictions_DOY(p0, doy.start, doy.end)
p1 <- transform_predictions_DOY(p1, doy.start, doy.end)
p2 <- transform_predictions_DOY(p2, doy.start, doy.end)
p3 <- transform_predictions_DOY(p3, doy.start, doy.end)

p0$x <- (p0$x * scale_factor) + doy.start

p0.peak <- create_effect_previous_year(data = p0,
                                       avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                                      x_col = "tr_peak_doy_lag1_plot",
                                      y_col = "peak_doy",
                                      line_color = "darkgreen",
                                      x_title = "Peak sw. timing lag1 [DOY]",
                                      y_title = lab_colonization_time,
                                      x_annotate = 160, lab_annotate = "***")

(p0.peak)
p1.peak <- 
  create_effect_plot(data = p1,
                     avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                     x_col = "tmp_z_lag2", y_col = "peak_doy", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag2 [z-score]")) , 
                     y_title = lab_peak_time, 
                       x_annotate = 2, 
                     lab_annotate = "0.71")

(p1.peak)
p2.peak <- create_effect_plot(data = p2,  avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220), 
                             x_col = "spei_lag1", y_col = "peak_doy", 
                             line_color = "blue", 
                             x_title = expression(paste("SPEI lag1 [dim.]")) , 
                             y_title = lab_peak_time, 
                             x_annotate = -0.5, 
                             lab_annotate = "**")
(p2.peak)
p3.peak <- plot_effect_interactions(p3, 
                                    avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                                    x_col = "tmp_z_lag2", 
                                    y_col = "peak_doy", 
                                   temp_label = expression(paste("Temp. lag2 [z-score]")), 
                                   y_title = lab_peak_time,
                                   x_annotate = 2,
                                   lab_annotate = "0.09") 

(p3.peak)



##### PLOT: Peak diff  ----------------------------------------------------
fin.m.peak.diff.previous.tw #<- fin.m.peak.diff$gam
summary(fin.m.peak.diff.previous.tw)
p0 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "peak_diff_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "spei_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.diff.previous.tw, terms = c("tmp_z_lag1", "spei_lag2 [-1, 0, 1]"), allow.new.levels = T)


divisor_diff <- 10
p0 <- adjust_predictions_counts(p0, divisor_diff)
p0$x <- p0$x/10
p1 <- adjust_predictions_counts(p1, divisor_diff)
p2 <- adjust_predictions_counts(p2, divisor_diff)
p3 <- adjust_predictions_counts(p3, divisor_diff)





p0.peak.diff <- create_effect_previous_year(data = p0, 
                                            avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4),
                                           x_col = "peak_diff_lag1",
                               y_col = "peak_diff",
                               line_color = "darkgreen", 
                               x_title = "Peak sw. intensity lag1 [*10]", 
                               y_title = lab_peak_growth,
                               x_annotate = 90, lab_annotate = "***")

(p0.peak.diff)
p1.peak.diff <- 
  create_effect_plot(data = p1, 
                     avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4), 
                     x_col = "tmp_z_lag1", y_col = "peak_diff", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag1 [z-score]")) , 
                     y_title = lab_peak_growth, 
                     x_annotate = 2, 
                     lab_annotate = "0.14")

(p1.peak.diff)
p2.peak.diff <- create_effect_plot(data = p2,
                                  avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4),
                                  x_col = "spei_lag2", y_col = "peak_diff", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI lag2 [dim.]")) , 
                               y_title = lab_peak_growth, 
                               x_annotate = -0.5, 
                               lab_annotate = "0.14")
(p2.peak.diff)
p3.peak.diff <- plot_effect_interactions(p3,
                                         avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4),
                                         x_col = "tmp_z_lag1", 
                                         y_col = "peak_diff", 
                                     temp_label = expression(paste("Temp. lag1 [z-score]")), 
                                     y_title = lab_peak_growth,
                                     x_annotate = 2,
                                     lab_annotate = "*") 

(p3.peak.diff)
ggarrange(p0.peak.diff,p1.peak.diff,p2.peak.diff,p3.peak.diff)


# put in the same figure: tmp, spei, interaction, previous years on teh bottom

p.previous <- ggarrange(p0.count, p0.agg, p0.peak, p0.peak.diff, 
                         ncol=4, nrow = 1 , #align = 'hv', 
                         font.label = list(size = 8, color = "black", face = "plain", family = NULL) )

(p.previous)

p.temp <-  ggarrange(p1.count, p1.agg, p1.peak, p1.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv', 
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.spei <-  ggarrange(p2.count, p2.agg, p2.peak, p2.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv', 
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.int <- ggarrange(p3.count,  p3.agg, p3.peak, p3.peak.diff,
                   ncol=4, nrow = 1 , align = 'hv',common.legend = TRUE, legend = 'bottom',
                   font.label = list(size = 8, color = "black", face = "plain", family = NULL))
windows(7,10)
full_preds <- ggarrange(p_spagetti,
                        p.temp, 
                        p.spei,  
                        p.previous,
                        p.int,
                        ncol = 1, nrow= 5, 
                        align = 'hv', heights = c(1.1, 1, 1, 1, 1.1),
                        widths = c(1, 1, 1, 1, 1))
#"#009E73" 
(full_preds)
ggsave(filename = 'outFigs/Fig_full_eff3.png', plot = full_preds, 
       width = 7.5, height = 10.5, dpi = 300, bg = 'white')








#### PLOT Moran's I -------------------------------
summary(fin.m.moran)

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.moran , terms = "tmp_z  [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.moran, terms = "spei_lag1  [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.moran, terms = "sum_ips [all]", allow.new.levels = TRUE)
p5 <- ggpredict(fin.m.moran, terms = "Morans_I_log_lag1 [all]", allow.new.levels = TRUE)

p4 <- ggpredict(fin.m.moran, terms = c("tmp_z", "spei_lag1 [-1,0,1]"), allow.new.levels = TRUE) 


p1.moran <- create_effect_plot(p1,
                               avg_data =  filter(lisa_merged_df_avg, Morans_I_log  > -1 & Morans_I_log  < 2), 
                               x_col = "tmp_z", 
                               y_col = "Morans_I_log", 
                               line_color = "red", 
                               x_title = "Temp. lag0 [z-score]", 
                               y_title = "Local Moran's I", 
                               x_annotate = 2,
                               lab_annotate = "**"
                               #y_lim = c(0,1)
                               ) +
  theme(axis.title.y = element_text("Local Moran's I", angle = 90))

p1.moran
p2.moran <- create_effect_plot(p2,
                               avg_data =  filter(lisa_merged_df_avg, Morans_I_log  > -1 & Morans_I_log  < 2), 
                               x_col = "spei_lag1", 
                               y_col = "Morans_I_log", 
                               line_color = "blue", 
                               x_title = "SPEI lag1 [dim.]", 
                               y_title = "Local Moran's I",  
                               #y_lim = c(0,1)
                               x_annotate = -0.5,
                               lab_annotate = "*"
                                                         ) +
  theme(axis.title.y = element_text("Local Moran's I", angle = 90))

p2.moran

lisa_merged_df_avg$sum_ips_scaled <- lisa_merged_df_avg$sum_ips / 1000
p3$x <- p3$x/1000 
p3.moran <- create_effect_plot(p3, 
                               avg_data =  filter(lisa_merged_df_avg, Morans_I_log  > -1 & Morans_I_log  < 2), 
                               x_col = "sum_ips_scaled", 
                               y_col = "Morans_I_log", 
                               line_color = "grey50", 
                               x_title = "Pop. level [#*1000]", 
                               y_title = "Local Moran's I", # y_lim = c(0,1),
                               x_annotate = 50,
                               lab_annotate = "***") +
  theme(axis.title.y = element_text("Local Moran's I", angle = 90))
p3.moran
p5.prev.Morans <- create_effect_previous_year(data = p5, 
                                            avg_data = filter(lisa_merged_df_avg, Morans_I_log  > -1 & Morans_I_log  < 2),
                                            x_col = "Morans_I_log_lag1",
                                            y_col = "Morans_I_log",
                                            line_color = "darkgreen", 
                                            x_title = "Local Moran's I lag1 [dim.]", 
                                            y_title =  "Local Moran's I [dim.]",
                                            x_annotate = 0, lab_annotate = "***") +
  theme(axis.title.y = element_text("Local Moran's I", angle = 90))

(p5.prev.Morans)


p4.moran <- plot_effect_interactions(p4,
                                     avg_data = filter(lisa_merged_df_avg, Morans_I_log  > -1 & Morans_I_log  < 2),
                                     
                                     x_col = "tmp_z", 
                                     y_col = "Morans_I_log", 
                                     temp_label = "Temp. lag0 [z-score]", 
                                     y_title = "Local Moran's I",
                                     x_annotate = 0.3,
                                     lab_annotate = "***") +
  theme(legend.position = c(0.6, 0.8),
        legend.title = element_text(hjust = 0.5),       # Align the legend title with the legend items (centered)
        #legend.title = "SPEI levels" 
        legend.background = element_rect(fill = "white", color = NA)  # White background with alpha
        ) +
  guides(color = guide_legend(ncol = 1)) +
  theme(axis.title.y = element_text("Local Moran's I", angle = 90))

 
p4.moran
  #geom_point(data = avg_data_moran_sub, aes(x = tmp_z_lag1, y = Morans_I_log))
p4.moran.no.leg <- p4.moran + theme(legend.position = 'none')

p4.moran.no.leg
# Step 2: Extract the legend as its own plot
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

my_legend <- get_legend(p4.moran)


##### test statistical differences and add to plot --------------- 
kruskal.test(Morans_I_log ~ year,dat_fin )

library(dunn.test)
posthoc_results <- dunn.test(dat_fin$Morans_I_log, dat_fin$year, method = "bonferroni")  # method can be adjusted for p-value correction
print(posthoc_results)


# Create a dataframe for annotation
annotation_data <- data.frame(year = 2015:2021, label = c('a', 'b',
                                                          'a', 'a', 'ac',
                                                          'c', 'a'))

p_morans_log_summary <- 
  ggplot(dat_fin, aes(x = year,
                      group = year,
                      y = Morans_I_log )) +
  stat_summary() +
  geom_hline(yintercept = 0, lty = 'dashed') +
  labs(y = "Local Moran's I") +
  coord_cartesian(ylim = c(-0.1, 0.5)) +
  # Ensure annotate uses explicit values for x and label
  geom_text(data = annotation_data, aes(x = year, y = Inf, label = label), 
            hjust = 0.5, vjust = 1.5, position = position_nudge(y = 0.05)) +
  theme_bw() + #geom_boxplot()
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"))  


p_morans_log_summary




# Arrange the plots side by side with the same size
p.effect.moran <- ggarrange( 
                            p1.moran, p2.moran,p4.moran, #p4.moran.no.leg, my_legend,
                            p3.moran,  
                            p5.prev.Morans,  #, 
                            p_morans_log_summary,
                            align = 'h',#p5.moran,
                            ncol = 3, nrow=2,
                            widths=c(1, 1, 1),
                            heights = c(1,1),
                            font.label = list(size = 8, color = "black", face = "plain", family = NULL),
                            labels = c( "[a]",
                                        "[b]",
                                        "[c]",
                                        "[d]",
                                        "[e]",
                                        "[f]"))


windows(7, 5)
(p.effect.moran)
ggsave(filename = 'outFigs/Fig5.png', plot = p.effect.moran, 
       width = 7, height = 5, dpi = 300, bg = 'white')




###### print models outputs: ---------------------------------------------------------

# export all models
sjPlot::tab_model(fin.m.counts.previous.tw,    file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg.doy.gamma,         file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak.doy.gamma,        file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff.previous.tw, file = "outTable/model_peak_diff.doc")
sjPlot::tab_model(fin.m.moran,                 file = "outTable/model_moran.doc")
#sjPlot::tab_model(fin.m.RS,        file = "outTable/model_RS_3.doc")



##### add significance ----------------------------------------------------
glob_merged <- glob_merged %>% 
  mutate(sign = case_when(p_val > 0.05 ~ 'n.s',
                          p_val < 0.05 & p_val > 0.01 ~ '*',
                          p_val < 0.01 & p_val > 0.001 ~ '**',
                          p_val < 0.001  ~ '***'
  )) 

# Create a barplot of the values
p_glob_moran_bar <-ggplot(glob_merged, aes(x = year,
                                           y = stat,
                                           label = sign)) +
  geom_col(fill = 'lightgrey', col = 'black') +
  xlim(2014.3,2021.6) +
  ylab('Moran I statistic\nstandard deviate')+
  ylim(0,13) +
  geom_text(aes(x = year,
                y = stat + 1 )) + #rep(12.2, 7))) +
  theme_bw() 



# END Global morans




