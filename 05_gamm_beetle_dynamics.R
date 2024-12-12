


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
#load(file =  "outData/lisa.Rdata")     # read LISA Moran's I stats
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

#nrow(lisa_merged_df)
nrow(dat_fin)
nrow(xy_df)

# add MOran's I values, transform the agg doy and peak doy between 0-1
dat_fin <-   dat_fin %>% 
  mutate(tr_agg_doy   = (agg_doy - doy.start) / (doy.end - doy.start),
         tr_peak_doy  = (peak_doy - doy.start) / (doy.end - doy.start))

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dat_fin$tr_agg_doy  <- pmin(pmax(dat_fin$tr_agg_doy, 1e-4),  1 - 1e-4)
dat_fin$tr_peak_doy <- pmin(pmax(dat_fin$tr_peak_doy, 1e-4), 1 - 1e-4)



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





## Pre-processs table: get lags ----------------------------

dat_spei_lags <-  dat_fin %>% 
  dplyr::select(c(year, pairID, trapID, tmp, 
                  spei,  # spei3
                  spei12,
                  spei24,
                  #tmp_z, 
                  sum_ips, 
                  #sum_ips_lag1, 
                  tr_agg_doy, tr_peak_doy,peak_diff,
                  agg_doy, peak_doy, 
                  x, y #,
                  #Morans_I_log 
                  )) %>% 
  group_by(trapID) %>%
  mutate(trapID = as.factor(trapID),
         pairID = as.factor(pairID),
         f_year = as.factor(year)) %>% 
  arrange(year, .by_group = TRUE) %>%
  # add lags
  mutate(sum_ips_lag1      = lag(sum_ips, n = 1, default = NA),
         tr_agg_doy_lag1   = lag(tr_agg_doy, n = 1, default = NA),
         tr_peak_doy_lag1  = lag(tr_peak_doy, n = 1, default = NA),
         peak_diff_lag1    = lag(peak_diff, n = 1, default = NA),
         tmp_lag1          = lag(tmp, n = 1, default = NA),
         spei12_lag1       = lag(spei12, n = 1, default = NA),
          ) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup(.) %>% 
  mutate(log_sum_ips = log(sum_ips+1),
         log_peak_diff = log(peak_diff + 1)) %>% 
  rename(tmp_lag0 = tmp,
         spei12_lag0 = spei12)

nrow(dat_spei_lags)
head(dat_spei_lags)


### remove the NAs for analysis 
dat_spei_lags_clean <- dat_spei_lags %>% 
  na.omit(dat_spei_lags)



# list predictors to test
selected_predictors_beetle <- c('sum_ips', 'sum_ips_lag1','sum_ips_lag2',
                         'tr_agg_doy'   , 'tr_agg_doy_lag1', 'tr_agg_doy_lag2' ,
                         'tr_peak_doy', 'tr_peak_doy_lag1','tr_peak_doy_lag2',
                         'peak_diff', 'peak_diff_lag1',  'peak_diff_lag2'
                         
) 
selected_predictors_climate <- c("spei12_lag0","spei12_lag1" , #"spei12_lag2",
                                 "tmp_lag0" , "tmp_lag1"#, "tmp_lag2" 
                                 )   


### check correlation between tmp and spei ---------------------------------------
# Calculate pairwise correlations for numeric columns
numeric_cols <- dat_spei_lags %>%
  dplyr::select_if(is.numeric)  %>% # Select only numeric columns
  dplyr::select(selected_predictors_climate #tmp_lag0,# tmp_lag1,tmp_lag2, 
                #spei12_lag0,  spei12_lag1 #, spei12_lag2
                )

# get correlation matrix between clim preictors
cor_matrix <- cor(numeric_cols, use = "pairwise.complete.obs", method = "spearman")
windows()
corrplot::corrplot(cor_matrix)


# export table to merge it with the spatial data: for variograms:
dat_dynamics <- dat_fin %>% 
  dplyr::select(c(year, trapID, pairID, spring_tmp, veg_tmp, veg_prcp,   
                  spei1,
                  spei3, 
                  spei12, # spei for veg season
                  #tmp, prcp,
                sum_ips,  
                agg_doy,
                peak_doy, 
                peak_diff, 
                tr_agg_doy, 
                tr_peak_doy, 
                                x,y
                ))
fwrite(dat_dynamics, 'outTable/beetle_dynamic_indicators.csv')


## Calculate average counts for each trap pair ------------------------------------
avg_data <- dat_spei_lags %>%
  group_by(pairID, year) %>%
  summarise(#Morans_I_log  = mean(Morans_I_log , na.rm = TRUE),
            
            # dependents
            sum_ips     = mean(sum_ips, na.rm = TRUE),
            peak_diff   = mean(peak_diff, na.rm = TRUE),
            tr_agg_doy  = mean(tr_agg_doy, na.rm = TRUE),
            tr_peak_doy = mean(tr_peak_doy, na.rm = TRUE), 
            agg_doy     = mean(agg_doy, na.rm = TRUE),
            peak_doy    = mean(peak_doy, na.rm = TRUE),
            # SPEIS - spei3
            spei_lag0   = mean(spei, na.rm = TRUE),
           # spei_lag1   = mean(spei_lag1, na.rm = TRUE), # for agg_doy
           # spei_lag2   = mean(spei_lag2, na.rm = TRUE),
            # # SPEI 12
            #spei24_lag0        = mean(spei24, na.rm = TRUE),
            spei12_lag0        = mean(spei12_lag0, na.rm = TRUE),
            spei12_lag1   = mean(spei12_lag1, na.rm = TRUE), # for agg_doy
            #spei12_lag2   = mean(spei12_lag2, na.rm = TRUE),
            # TMP
            tmp_lag0         = mean(tmp_lag0, na.rm = T),
            #tmp_lag1    = mean(tmp_lag1, na.rm = TRUE),
            #tmp_lag2    = mean(tmp_lag2, na.rm = TRUE),  # for agg_doy
            
            x           = mean(x, na.rm = TRUE),
            y           = mean(y, na.rm = TRUE)
  ) %>%
  ungroup(.) %>% 
  na.omit() %>% 
  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% c('Eschenbach_idOPf', 'Peiting') ) # remove two traps with always having too low trap counts

fwrite(avg_data, 'outTable/fin_tab_avg.csv')




hist(dat_dynamics$veg_tmp)
hist(dat_dynamics$spei12)


### make a plot for each variable per year ------------------------------------------
# Melt the data to long format for easier plotting with facets
# Rename columns for clarity
dat_dynamics2 <- dat_dynamics %>%
  dplyr::rename(temperature = veg_tmp, 
                precipitation = veg_prcp, 
                SPEI = spei3)

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
weather_summary_plot <- ggplot(plot_data, aes(x = year, 
                                              y = value#, 
                                              #color = variable, 
                                              #fill = variable
)) +
  geom_rect(
    aes(xmin = 2017.9, xmax = 2020.1, ymin = -Inf, ymax = Inf),
    fill = "grey90", alpha = 0.5, inherit.aes = FALSE
  ) +  # Add grey rectangle for 2018-2020
  geom_violin(aes(group = as.factor(year)), alpha = 0.5, fill = 'grey60', color = 'grey60') +  # Violin plot to show distribution
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





# Summary stats ---------------------------------------------------
## weather condistions: --------------------------
desc_stat_climate_summary <- dat_dynamics %>%
  mutate(class = ifelse(year %in% c(2018:2020), "2018-2020", "Other Years")) %>%
  group_by(class) %>%
  summarise(
    mean_veg_tmp = mean(veg_tmp, na.rm = TRUE),
    med_veg_tmp = median(veg_tmp, na.rm = TRUE),
    sd_veg_tmp = sd(veg_tmp, na.rm = TRUE),
    mean_spei1 = mean(spei1, na.rm = TRUE),
    med_spei1 = median(spei1, na.rm = TRUE),
    sd_spei1 = sd(spei1, na.rm = TRUE),
    mean_spei12 = mean(spei12, na.rm = TRUE),
    med_spei12 = median(spei12, na.rm = TRUE),
    sd_spei12 = sd(spei12, na.rm = TRUE),
    
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
   
     mean_spei1 = mean(spei1, na.rm = TRUE),
    med_spei1 = median(spei1, na.rm = TRUE),
    sd_spei1 = sd(spei1, na.rm = TRUE),
    
    mean_spei12 = mean(spei12, na.rm = TRUE),
    med_spei12 = median(spei12, na.rm = TRUE),
    sd_spei12 = sd(spei12, na.rm = TRUE),
    
    mean_veg_prcp = mean(veg_prcp, na.rm = TRUE),
    med_veg_prcp = median(veg_prcp, na.rm = TRUE),
    sd_veg_prcp = sd(veg_prcp, na.rm = TRUE)
  )

(desc_stat_climate_summary_year)



## for beetle population indicators ------------------------------------

### quantils  ----------------------------------------------------------

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



### Beetle indcators per year: -------------------------------------------------------


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
  mutate(Population_level       = stringr::str_glue("{round(mean_sum_ips,0)}±{round(sd_sum_ips,0)}"),
         Aggregation_timing       = stringr::str_glue("{round(mean_agg_doy,0)}±{round(sd_agg_doy,0)}"),
         Peak_swarming_timing     = stringr::str_glue("{round(mean_peak_doy,0)}±{round(sd_peak_doy,0)}"),
         Peak_swarming_intensity  = stringr::str_glue("{round(mean_peak_diff,0)}±{round(sd_peak_diff,0)}")) %>% 
  dplyr::select(year, Population_level, Aggregation_timing, Peak_swarming_timing,Peak_swarming_intensity) 


(means_dat_fin_year)


# Export as a nice table in word:
sjPlot::tab_df(means_dat_fin_year,
               #               col.header = c(as.character(qntils), 'mean'),
               show.rownames = TRUE,
               file="outTable/summary_out_year.doc",
               digits = 0) 

### Beetle indicators per drought/no drought ----------------------------------------

means_dat_fin_drought <- 
  dat_fin %>% 
  mutate(drought_status = ifelse(year %in% 2018:2020, "Hotter drought", "No drought")) %>% 
  ungroup(.) %>% 
  filter(year %in% 2015:2021) %>% 
  group_by(drought_status) %>% 
  dplyr::reframe(mean_sum_ips   = mean(sum_ips, na.rm = T),
                 mean_agg_doy   = mean(agg_doy, na.rm = T),
                 mean_peak_doy  = mean(peak_doy, na.rm = T),
                 mean_peak_diff = mean(peak_diff, na.rm = T),
                 sd_sum_ips     = sd(sum_ips, na.rm = T),
                 sd_agg_doy     = sd(agg_doy, na.rm = T),
                 sd_peak_doy    = sd(peak_doy, na.rm = T),
                 sd_peak_diff   = sd(peak_diff, na.rm = T)) %>% 
  mutate(Population_level       = stringr::str_glue("{round(mean_sum_ips,0)}±{round(sd_sum_ips,0)}"),
         Aggregation_timing       = stringr::str_glue("{round(mean_agg_doy,0)}±{round(sd_agg_doy,0)}"),
         Peak_swarming_timing     = stringr::str_glue("{round(mean_peak_doy,0)}±{round(sd_peak_doy,0)}"),
         Peak_swarming_intensity  = stringr::str_glue("{round(mean_peak_diff,0)}±{round(sd_peak_diff,0)}")) %>% 
  dplyr::select(drought_status, Population_level, Aggregation_timing, Peak_swarming_timing,Peak_swarming_intensity) 


(means_dat_fin_drought)


# Export as a nice table in word:
sjPlot::tab_df(means_dat_fin_drought,
               show.rownames = TRUE,
               file="outTable/summary_out_drought.doc",
               digits = 0) 


# get medians and IQR
medians_dat_fin_drought <- 
  dat_fin %>% 
  mutate(drought_status = ifelse(year %in% 2018:2020, "Hotter drought", "No drought")) %>% 
  ungroup() %>% 
  dplyr::filter(year %in% 2015:2021) %>% 
  group_by(drought_status) %>% 
  dplyr::reframe(
    median_sum_ips   = median(sum_ips, na.rm = TRUE),
    iqr_sum_ips      = IQR(sum_ips, na.rm = TRUE),
    median_agg_doy   = median(agg_doy, na.rm = TRUE),
    iqr_agg_doy      = IQR(agg_doy, na.rm = TRUE),
    median_peak_doy  = median(peak_doy, na.rm = TRUE),
    iqr_peak_doy     = IQR(peak_doy, na.rm = TRUE),
    median_peak_diff = median(peak_diff, na.rm = TRUE),
    iqr_peak_diff    = IQR(peak_diff, na.rm = TRUE)
  ) %>% 
  mutate(
    Population_level       = stringr::str_glue("{round(median_sum_ips, 0)} ({round(iqr_sum_ips, 0)})"),
    Aggregation_timing      = stringr::str_glue("{round(median_agg_doy, 0)} ({round(iqr_agg_doy, 0)})"),
    Peak_swarming_timing    = stringr::str_glue("{round(median_peak_doy, 0)} ({round(iqr_peak_doy, 0)})"),
    Peak_swarming_intensity = stringr::str_glue("{round(median_peak_diff, 0)} ({round(iqr_peak_diff, 0)})")
  ) %>% 
  dplyr::select(drought_status, Population_level, Aggregation_timing, Peak_swarming_timing, Peak_swarming_intensity)

# Preview the table
print(medians_dat_fin_drought)

# Export as a nice table in Word
sjPlot::tab_df(
  medians_dat_fin_drought,
  show.rownames = TRUE,
  file = "outTable/summary_out_drought_medians.doc",
  digits = 0
)






#Filter data for plotting --------------------------------------------------------
# filter extreme values: occuring too late in a year ()
avg_data_agg <- avg_data %>% 
  dplyr::filter(tr_agg_doy < 0.54)

# remove outliers: Piding
avg_data_peak_diff <- avg_data %>% 
  dplyr::filter(!pairID %in% c('Piding', 'Gangkofen'))

avg_data_filt <- avg_data #%>% 
  #mutate(f_year = as.factor(year)) %>% 
  #dplyr::filter(!pairID %in% pair_outliers )


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




# makde models with diffr autocorrelation structur (direct vs values from previous years) and compare their performacce 

# Drivers -----------------------------------------------------------------------
# Final predictors: current year tmp_lag0 & spei12_lag1 best predictors
# 
cor(avg_data$spei12_lag1, avg_data$tmp_lag0, method = 'spearman')
#cor(avg_data$spei24_lag0, avg_data$tmp_lag0)

## IPS_SUM ---------------------------------------------------------------------- 

m2_full_ti <- gam(sum_ips ~ 
                 tmp_lag0 +
                 #s(tmp_lag0, k = 3) +
                 s(spei12_lag1, k = 3) + 
                 ti(tmp_lag0, spei12_lag1, k = 3) + 
                 s(sum_ips_lag1, k = 5) + # Added term for previous beetle counts
                 s(x,y) +
                 s(pairID, bs = 're') +
                 s(f_year, bs = 're')
               ,
               data = dat_spei_lags_clean, #avg_data_filt_lagged, 
               family = tw,
               method = "REML",
               select = TRUE)

AIC(m2_full_ti, m2_full)
summary(m2_full_ti)

fin.m.counts.previous.tw <- m2_full_ti

## previous peak_diff  ------------------------------------------------------

m.peak.diff.previous.tmp_0_spei12_1 <- gam(peak_diff ~ 
                                             tmp_lag0 +
                                              #  s(tmp_lag0, k = 3) +
                                                s(spei12_lag1, k = 3) + 
                                                ti(tmp_lag0, spei12_lag1, k = 3) + 
                                                s(peak_diff_lag1 , k =5) + # Added term for previous beetle counts
                                                s(x,y) + #,
                                              s(pairID, bs = 're') +
                                              s(f_year, bs = 're'),
                                              data = dat_spei_lags_clean, #avg_data_filt_lagged, 
                                              family = tw,
                                           select = TRUE)

m<-m.peak.diff.previous.tmp_0_spei12_1

#appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('peak_diff_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1[-1.5,0]'))
plot(predict_data)


## AGG DOY ----------------------------------------------------------


m.agg.previous_tmp0_spei12_1 <- gam(tr_agg_doy ~
                                               s(tmp_lag0, k = 3) +
                                               s(spei12_lag1, k = 3) +
                                               ti(tmp_lag0, spei12_lag1, k = 3) +
                                               s(tr_agg_doy_lag1 , k = 5) +# Added term for previous year numbers to account ofr autocorrelation
                                               s(x,y) +  
                                               s(pairID, bs = 're') +
                                               s(f_year, bs = 're')
                                             ,
                                             data = dat_spei_lags_clean, 
                                             family = Gamma(link = "log"),
                                    select = TRUE)


m<-m.agg.previous_tmp0_spei12_1
k.check(m)

appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tr_agg_doy_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1[-1.5,0]'))
plot(predict_data)






## Previous PEak DoY ----------------------------------------------------

m.peak.doy.previous_tmp0_spei12_1 <- gam(tr_peak_doy ~ 
                                              s(tmp_lag0, k = 5) +
                                              s(spei12_lag1, k = 5) +
                                              ti(tmp_lag0, spei12_lag1, k = 3) +
                                              s(tr_peak_doy_lag1 , k = 5) + # Added term for previous year numbers to account ofr autocorrelation
                                              s(x,y, bs = 'gp')   +
                                              s(pairID, bs = 're') +
                                            s(f_year, bs = 're')
                                            ,
                                            data = dat_spei_lags_clean, # avg_data_peak_no_out, 
                                            family = Gamma(link = "log"))


m<-m.peak.doy.previous_tmp0_spei12_1
k.check(m)

appraise(m)
summary(m)

predict_data <- ggpredict(m, terms = c('tmp_lag0'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('spei12_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tr_peak_doy_lag1'))
plot(predict_data)

predict_data <- ggpredict(m, terms = c('tmp_lag0', 'spei12_lag1[-1.5,0]'))
plot(predict_data)


# save final models : current year TMP and previous year SPEI12 ---------------------------------------
fin.m.counts.previous.tw    <- m2_full_ti
fin.m.peak.diff.previous.tw <- m.peak.diff.previous.tmp_0_spei12_1
fin.m.agg.doy.gamma         <- m.agg.previous_tmp0_spei12_1
fin.m.peak.doy.gamma        <- m.peak.doy.previous_tmp0_spei12_1

summary(fin.m.counts.previous.tw)   
summary(fin.m.peak.diff.previous.tw)
summary(fin.m.agg.doy.gamma)         
summary(fin.m.peak.doy.gamma)        


# check for spatial and tmp autocorrelaton ----------------------------------
# Extract residuals
# Extract residuals from the model
residuals <- resid(m.peak.diff.previous.tmp_z_0_spei12_1, type = "normalized")


# Check for temporal autocorrelation
par(mfrow = c(1, 2))  # Arrange plots in one row, two columns
acf(residuals, main = "ACF of Residuals")  # Autocorrelation function
pacf(residuals, main = "PACF of Residuals")  # Partial autocorrelation function



####### Find predicors &  lags: MOrans  -------------

# update& save the best models: --------------------------------------------------------

save(fin.m.counts.previous.tw,
     fin.m.peak.diff.previous.tw,
     fin.m.agg.doy.gamma,
     fin.m.peak.doy.gamma,
         avg_data,
     dat_spei_lags_clean,
     file = "outData/fin_models.RData")




# Effect plots ------------------------------------------------------------


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


# get predictions per one dergee celsius/per change in SPEI  -------------

### IPS COUNT -----------------------
predictions_tmp_counts <- ggpredict(fin.m.counts.previous.tw, terms = "tmp_lag0[11:18]")
predictions_spei_counts <- ggpredict(fin.m.counts.previous.tw, terms = "spei12_lag1[-1,0]")

# Calculate the increase per 1°C for TMP
increase_per_degree_tmp <- (predictions_tmp_counts$predicted[8] - predictions_tmp_counts$predicted[1]) / 7
total_increase_tmp <- predictions_tmp_counts$predicted[8] - predictions_tmp_counts$predicted[1]
increase_per_degree_tmp
total_increase_tmp

# Predictions for SPEI
predictions_spei_counts <- ggpredict(fin.m.counts.previous.tw, terms = "spei12_lag1[-1,0]")

# Calculate the increase per unit change in SPEI
increase_per_unit_spei <- (predictions_spei_counts$predicted[2] - predictions_spei_counts$predicted[1]) / 1
total_increase_spei <- predictions_spei_counts$predicted[2] - predictions_spei_counts$predicted[1]
increase_per_unit_spei
total_increase_spei

# int-------------
# Generate predictions for the interaction effect of tmp_lag0 and spei12_lag1 at SPEI = -1
predictions_interaction <- ggpredict(fin.m.counts.previous.tw, terms = c("tmp_lag0[11:18]", "spei12_lag1[-1]"))

# Calculate the increase per 1°C for TMP when SPEI = -1
increase_per_degree_tmp_spei <- (predictions_interaction$predicted[8] - predictions_interaction$predicted[1]) / 7
total_increase_tmp_spei <- predictions_interaction$predicted[8] - predictions_interaction$predicted[1]

# Print the results
cat("Increase per 1°C (TMP | SPEI = -1):", increase_per_degree_tmp_spei, "\n")
cat("Total increase (TMP 11°C to 18°C | SPEI = -1):", total_increase_tmp_spei, "\n")

### AGG DOY --------------


# Generate predictions for temperature range (TMP)
predictions_tmp_agg_doy <- ggpredict(fin.m.agg.doy.gamma, terms = "tmp_lag0[11:18]")

# Transform predictions to DOY scale
predictions_tmp_agg_doy_transformed <- transform_predictions_DOY(
  predictions_tmp_agg_doy,
  doy.start = 91,
  doy.end = 273
)

# Calculate total shift in aggregation timing for TMP
total_shift_tmp <- predictions_tmp_agg_doy_transformed$predicted[8] - predictions_tmp_agg_doy_transformed$predicted[1]

# Calculate per-degree shift for TMP
shift_per_degree_tmp <- total_shift_tmp / (18 - 11)

# Print TMP results
cat("AGG DOY TMP Total shift (11°C to 18°C):", total_shift_tmp, "\n")
cat("AGG DOY TMP Shift per degree:", shift_per_degree_tmp, "\n")

# Generate predictions for SPEI range (-1 to 0)
predictions_spei_agg_doy <- ggpredict(fin.m.agg.doy.gamma, terms = "spei12_lag1[-1,0]")

# Transform predictions to DOY scale
predictions_spei_agg_doy_transformed <- transform_predictions_DOY(
  predictions_spei_agg_doy,
  doy.start = 91,
  doy.end = 273
)

# Calculate total shift in aggregation timing for SPEI
total_shift_spei <- predictions_spei_agg_doy_transformed$predicted[2] - predictions_spei_agg_doy_transformed$predicted[1]

# Calculate per-unit shift for SPEI
shift_per_unit_spei <- total_shift_spei / 1  # Change is between -1 and 0, so divide by 1

# Print SPEI results
cat("AGG DOY SPEI Total shift (-1 to 0):", total_shift_spei, "\n")
cat("AGG DOY SPEI Shift per unit SPEI:", shift_per_unit_spei, "\n")


# for peak doy ---------------------
# Generate predictions for temperature range (TMP)
predictions_tmp_peak_doy <- ggpredict(fin.m.peak.doy.gamma, terms = "tmp_lag0[11:18]")

# Transform predictions to DOY scale
predictions_tmp_peak_doy_transformed <- transform_predictions_DOY(
  predictions_tmp_peak_doy,
  doy.start = 91,
  doy.end = 273
)

# Calculate total shift in peak timing for TMP
total_shift_tmp_peak <- predictions_tmp_peak_doy_transformed$predicted[8] - predictions_tmp_peak_doy_transformed$predicted[1]

# Calculate per-degree shift for TMP
shift_per_degree_tmp_peak <- total_shift_tmp_peak / (18 - 11)

# Print TMP results
cat("PEAK DOY TMP Total shift (11°C to 18°C):", total_shift_tmp_peak, "\n")
cat("PEAK DOY TMP Shift per degree:", shift_per_degree_tmp_peak, "\n")

# Generate predictions for SPEI range (-1 to 0)
predictions_spei_peak_doy <- ggpredict(fin.m.peak.doy.gamma, terms = "spei12_lag1[-1,0]")

# Transform predictions to DOY scale
predictions_spei_peak_doy_transformed <- transform_predictions_DOY(
  predictions_spei_peak_doy,
  doy.start = 91,
  doy.end = 273
)

# Calculate total shift in peak timing for SPEI
total_shift_spei_peak <- predictions_spei_peak_doy_transformed$predicted[2] - predictions_spei_peak_doy_transformed$predicted[1]

# Calculate per-unit shift for SPEI
shift_per_unit_spei_peak <- total_shift_spei_peak / 1  # Change is between -1 and 0, so divide by 1

# Print SPEI results
cat("PEAK DOY SPEI Total shift (-1 to 0):", total_shift_spei_peak, "\n")
cat("PEAK DOY SPEI Shift per unit SPEI:", shift_per_unit_spei_peak, "\n")


# ### PEAK DIFFERENCE -----------------------

# Generate predictions for TMP (temperature) range
predictions_tmp_diff <- ggpredict(fin.m.peak.diff.previous.tw, terms = "tmp_lag0[11:18]")

# Calculate the increase in difference per 1°C for TMP
increase_per_degree_tmp_diff <- (predictions_tmp_diff$predicted[8] - predictions_tmp_diff$predicted[1]) / 7
total_increase_tmp_diff <- predictions_tmp_diff$predicted[8] - predictions_tmp_diff$predicted[1]

# Print TMP results
cat("PEAK DIFF TMP Increase per degree (11°C to 18°C):", increase_per_degree_tmp_diff, "\n")
cat("PEAK DIFF TMP Total increase (11°C to 18°C):", total_increase_tmp_diff, "\n")

# Generate predictions for SPEI range (-1 to 0)
predictions_spei_diff <- ggpredict(fin.m.peak.diff.previous.tw, terms = "spei12_lag1[-1,0]")

# Calculate the increase in difference per unit change in SPEI
increase_per_unit_spei_diff <- (predictions_spei_diff$predicted[2] - predictions_spei_diff$predicted[1]) / 1
total_increase_spei_diff <- predictions_spei_diff$predicted[2] - predictions_spei_diff$predicted[1]

# Print SPEI results
cat("PEAK DIFF SPEI Increase per unit (-1 to 0):", increase_per_unit_spei_diff, "\n")
cat("PEAK DIFF SPEI Total increase (-1 to 0):", total_increase_spei_diff, "\n")






## PLOT: IPS SUM counts -----------------------------------------------------------


# Assuming 'model' is your glm.nb model
summary(fin.m.counts.previous.tw)



#### plot --------------------------------

m<-fin.m.counts.previous.tw

summary(m)

p_val_tmp           <- format_p_value_label(summary(m), "tmp_lag0")
p_val_spei          <- format_p_value_label(summary(m), "s(spei12_lag1)")
p_val_int           <- format_p_value_label(summary(m), "ti(tmp_lag0,spei12_lag1)")
p_val_lag_ips_count     <- format_p_value_label(summary(m), "s(sum_ips_lag1)")


p0 <- ggpredict(fin.m.counts.previous.tw, terms = "sum_ips_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.counts.previous.tw, terms = "tmp_lag0 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts.previous.tw, terms = "spei12_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.counts.previous.tw, terms = c("tmp_lag0", "spei12_lag1 [-1.5, 0]"), allow.new.levels = TRUE)

divisor_population <- 1000
p0 <- adjust_predictions_counts(p0, divisor_population)
p0$x <- p0$x/1000
p1 <- adjust_predictions_counts(p1, divisor_population)
p2 <- adjust_predictions_counts(p2, divisor_population)
p3 <- adjust_predictions_counts(p3, divisor_population)





my_colors_interaction <- c("#A50026", "#3366CC" ) # "#FDAE61"
#my_colors_interaction <- c("#A50026", "#FDAE61" )


p0.count <- create_effect_previous_year(data = p0, 
                              avg_data = avg_data_filt_lagged_plotting,
                              x_col = "sum_ips_lag1",
                              y_col = "sum_ips",
                              line_color =  "darkgrey",#"#229C52", 
                              x_title = "Pop. level lag1 [*1000]", 
                              y_title = paste(lab_popul_level, '*1000'),# "Population level\n[# beetle*100]",
                              #my_title = paste("[a]", 'Population level', '\n[#*100]'), 
                              x_annotate = 75, lab_annotate = p_val_lag_ips_count)

(p0.count)
p1.count <- 
  create_effect_plot(data = p1, 
                     avg_data = avg_data_filt_lagged_plotting, 
                     x_col = "tmp_lag0", 
                     y_col = "sum_ips", 
                     line_color = "#A50026", 
                     x_title = temp_label,
                     y_title = paste(lab_popul_level, '*1000'), 
                     x_annotate = 15, 
                     lab_annotate = p_val_tmp)

(p1.count)
p2.count <- create_effect_plot(data = p2, avg_data = avg_data_filt_lagged_plotting, 
                               x_col = "spei12_lag1", y_col = "sum_ips", 
                               line_color = "#FDAE61", 
                               x_title = spei_label,
                               y_title = paste(lab_popul_level, '*1000'), 
                               x_annotate = -0.5, 
                               lab_annotate = p_val_spei)
(p2.count)
p3.count <- plot_effect_interactions(data = p3, 
                                     avg_data = avg_data_filt_lagged_plotting, 
                                     x_col = "tmp_lag0", 
                                     y_col = "sum_ips", 
                                     temp_label = temp_label, 
                                     y_title = paste(lab_popul_level, '*1000'),
                                     x_annotate = 15,
                                     lab_annotate = p_val_int)  +
  scale_color_manual(values = my_colors_interaction, name = "SPEI") +
  scale_fill_manual(values = my_colors_interaction, name = "SPEI") #+

(p3.count)
ggarrange(p0.count,p1.count,p2.count, p3.count,
          labels = c("[a]", "[b]", "[c]", "[d]"), # Custom labels
          label.x = 0.1, # Adjust x-position (0 is left, 1 is right)
          label.y = 0.9, # Adjust y-position (1 is top, 0 is bottom)
          font.label = list(size = 8, face = "plain"))




##### PLOT DOY aggregation ---------------------------------------------------------
summary(fin.m.agg.doy.gamma)




m<-fin.m.agg.doy.gamma

summary(m)

# get p values 
p_val_tmp           <- format_p_value_label(summary(m), "s(tmp_lag0)")
p_val_spei          <- format_p_value_label(summary(m), "s(spei12_lag1)")
p_val_int           <- format_p_value_label(summary(m), "ti(tmp_lag0,spei12_lag1)")
p_val_agg_doy     <- format_p_value_label(summary(m), "s(tr_agg_doy_lag1)")
(p_val_agg_doy)


# predict
p0 <- ggpredict(fin.m.agg.doy.gamma, terms = "tr_agg_doy_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.agg.doy.gamma, terms = "tmp_lag0 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg.doy.gamma, terms = "spei12_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.agg.doy.gamma, terms = c("tmp_lag0", "spei12_lag1 [-1.5, 0]"), allow.new.levels = TRUE)

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
                               line_color = "darkgrey",
                               x_title = "Agg. timing lag1 [DOY]",
                               y_title = lab_colonization_time,
                               x_annotate = 175, lab_annotate = p_val_agg_doy )

(p0.agg)
p1.agg <- 
  create_effect_plot(data = p1,  avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),  
                     x_col = "tmp_lag0", y_col = "agg_doy", 
                     line_color = "#A50026", 
                     x_title = temp_label,
                     y_title = lab_colonization_time, 
                      x_annotate = 15, 
                     lab_annotate = p_val_tmp)

(p1.agg)
p2.agg <- create_effect_plot(data = p2, avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),    
                               x_col = "spei12_lag1", y_col = "agg_doy", 
                               line_color = "#FDAE61", 
                               x_title = spei_label,
                               y_title = lab_colonization_time, 
                               x_annotate = -0.5, 
                               lab_annotate = p_val_spei)
(p2.agg)
p3.agg <- plot_effect_interactions(p3,
                                   avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),
                                   x_col = "tmp_lag0", 
                                   y_col = "agg_doy", 
                                     temp_label = temp_label, #expression(paste("Temp. lag2 [z-score]")), 
                                     y_title = lab_colonization_time,
                                     x_annotate = 15,
                                     lab_annotate = p_val_int) +
  scale_color_manual(values = my_colors_interaction, name = "SPEI") +
  scale_fill_manual(values = my_colors_interaction, name = "SPEI") #+

(p3.agg)


#### PLOT: Peak DOY  ------------------------------------------------------------
summary(fin.m.peak.doy.gamma)


m<-fin.m.peak.doy.gamma

summary(m)

# get p values 
p_val_tmp           <- format_p_value_label(summary(m), "s(tmp_lag0)")
p_val_spei          <- format_p_value_label(summary(m), "s(spei12_lag1)")
p_val_int           <- format_p_value_label(summary(m), "ti(tmp_lag0,spei12_lag1)")
p_val_peak_doy     <- format_p_value_label(summary(m), "s(tr_peak_doy_lag1)")




p0 <- ggpredict(fin.m.peak.doy.gamma, terms = "tr_peak_doy_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak.doy.gamma, terms = "tmp_lag0 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.doy.gamma, terms = "spei12_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.doy.gamma, terms = c("tmp_lag0" ,"spei12_lag1 [-1.5, 0]"), allow.new.levels = TRUE)


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
                                      line_color = "darkgrey",
                                      x_title = "Peak sw. timing lag1 [DOY]",
                                      y_title = lab_colonization_time,
                                      x_annotate = 160, lab_annotate = p_val_peak_doy)

(p0.peak)
p1.peak <- 
  create_effect_plot(data = p1,
                     avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                     x_col = "tmp_lag0", y_col = "peak_doy", 
                     line_color = "#A50026", 
                     x_title = temp_label, #expression(paste("Temp. lag2 [z-score]")) , 
                     y_title = lab_peak_time, 
                       x_annotate = 15, 
                     lab_annotate = p_val_tmp)

(p1.peak)
p2.peak <- create_effect_plot(data = p2,  avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220), 
                             x_col = "spei12_lag1", y_col = "peak_doy", 
                             line_color = "#FDAE61", 
                             x_title = spei_label, 
                             y_title = lab_peak_time, 
                             x_annotate = -0.5, 
                             lab_annotate = p_val_spei)
(p2.peak)
p3.peak <- plot_effect_interactions(p3, 
                                    avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                                    x_col = "tmp_lag0", 
                                    y_col = "peak_doy", 
                                   temp_label = temp_label, #expression(paste("Temp. lag2 [z-score]")), 
                                   y_title = lab_peak_time,
                                   x_annotate = 15,
                                   lab_annotate = p_val_int) +
  scale_color_manual(values = my_colors_interaction, name = "SPEI") +
  scale_fill_manual(values = my_colors_interaction, name = "SPEI") #+

(p3.peak)



##### PLOT: Peak diff  ----------------------------------------------------

summary(fin.m.peak.diff.previous.tw)


m<-fin.m.peak.diff.previous.tw

summary(m)

# get p values 
p_val_tmp           <- format_p_value_label(summary(m), "tmp_lag0")
p_val_spei          <- format_p_value_label(summary(m), "s(spei12_lag1)")
p_val_int           <- format_p_value_label(summary(m), "ti(tmp_lag0,spei12_lag1)")
p_val_peak_diff     <- format_p_value_label(summary(m), "s(peak_diff_lag1)")



p0 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "peak_diff_lag1 [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "tmp_lag0 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff.previous.tw, terms = "spei12_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.diff.previous.tw, terms = c("tmp_lag0", "spei12_lag1 [-1.5, 0]"), allow.new.levels = T)


divisor_diff <- 10
p0 <- adjust_predictions_counts(p0, divisor_diff)
p0$x <- p0$x/10
p1 <- adjust_predictions_counts(p1, divisor_diff)
p2 <- adjust_predictions_counts(p2, divisor_diff)
p3 <- adjust_predictions_counts(p3, divisor_diff)





p0.peak.diff <- create_effect_previous_year(data = p0, 
                                            avg_data = avg_data_filt_lagged_plotting, #dplyr::filter(avg_data_filt_lagged_plotting, tmp_lag0  < 18),
                                           x_col = "peak_diff_lag1",
                               y_col = "peak_diff",
                               line_color = "darkgrey", 
                               x_title = "Peak sw. intensity lag1 [*10]", 
                               y_title = lab_peak_growth,
                               x_annotate = 150, lab_annotate = p_val_peak_diff)

(p0.peak.diff)
p1.peak.diff <- 
  create_effect_plot(data = p1, 
                     avg_data = avg_data_filt_lagged_plotting, #filter(avg_data_filt_lagged_plotting, tmp_lag0  <18), 
                     x_col = "tmp_lag0", y_col = "peak_diff", 
                     line_color = "#A50026", 
                     x_title = temp_label, 
                     y_title = lab_peak_growth, 
                     x_annotate = 15, 
                     lab_annotate = p_val_tmp)

(p1.peak.diff)
p2.peak.diff <- create_effect_plot(data = p2,
                                  avg_data = avg_data_filt_lagged_plotting, #filter(avg_data_filt_lagged_plotting, tmp_lag0  < 4),
                                  x_col = "spei12_lag1", y_col = "peak_diff", 
                               line_color = "#FDAE61", 
                               x_title = spei_label,  
                               y_title = lab_peak_growth, 
                               x_annotate = -0.5, 
                               lab_annotate = p_val_spei)
(p2.peak.diff)
p3.peak.diff <- plot_effect_interactions(p3,
                                         avg_data = avg_data_filt_lagged_plotting, #filter(avg_data_filt_lagged_plotting, tmp_lag0  < 4),
                                         x_col = "tmp_lag0", 
                                         y_col = "peak_diff", 
                                     temp_label = temp_label, 
                                     y_title = lab_peak_growth,
                                     x_annotate = 15,
                                     lab_annotate = p_val_int)  +
  scale_color_manual(values = my_colors_interaction, name = "SPEI") +
  scale_fill_manual(values = my_colors_interaction, name = "SPEI") #+

(p3.peak.diff)
ggarrange(p0.peak.diff,p1.peak.diff,p2.peak.diff,p3.peak.diff)


### Aggregeate plots  -----
dat_fin_plot <- dat_fin 

dat_fin_plot$sum_ips <- dat_fin$sum_ips/1000  
dat_fin_plot$peak_diff <- dat_fin$peak_diff/10 

p_spagett_ips       <- plot_data_with_average(dat_fin_plot, "sum_ips", lab_popul_level,   
                                              my_title = paste( 'Pop. level', '[#*1000]'))
p_spagett_ips
p_spagett_agg_doy   <- plot_data_with_average(dat_fin_plot, "agg_doy", lab_colonization_time, 
                                              my_title = paste('Aggregation timing', '[DOY]'))
p_spagett_peak_doy  <- plot_data_with_average(dat_fin_plot, "peak_doy", lab_peak_time, 
                                              my_title = paste('Peak sw. timing', '[DOY]'))
p_spagett_peak_diff <- plot_data_with_average(dat_fin_plot, "peak_diff", lab_peak_growth, 
                                              my_title = paste('Peak sw. intensity', '[#*10]'))




#windows(7,6)
p_spagetti <- ggarrange(p_spagett_ips, 
                        p_spagett_agg_doy, 
                        p_spagett_peak_doy, 
                        p_spagett_peak_diff, ncol = 4,nrow = 1, align = 'hv',
                        labels = c("[a]", "[b]", "[c]", "[d]"), # Custom labels
                        label.x = 0.15, # Adjust x-position (0 is left, 1 is right)
                        label.y = 0.85, # Adjust y-position (1 is top, 0 is bottom)
                        font.label = list(size = 8, face = "plain"))
p_spagetti
#ggsave(filename = 'outFigs/Fig1.png', plot = p_spagetti, width = 7, height = 6, dpi = 300, bg = 'white')



p.previous <- ggarrange(p0.count, p0.agg, p0.peak, p0.peak.diff, 
                        labels = c("[e]", "[f]", "[g]", "[h]"), # Custom labels
                        ncol=4, nrow = 1 , #align = 'hv', 
                        label.x = 0.15, # Adjust x-position (0 is left, 1 is right)
                        label.y = 0.85, # Adjust y-position (1 is top, 0 is bottom)
                        font.label = list(size = 8, color = "black", face = "plain", family = NULL) )

(p.previous)

p.temp <-  ggarrange(p1.count, p1.agg, p1.peak, p1.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv',
                     label.x = 0.15, # Adjust x-position (0 is left, 1 is right)
                     label.y = 0.85, # Adjust y-position (1 is top, 0 is bottom)
                     labels = c("[e]", "[f]", "[g]", "[h]"), # Custom labels
                     
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.spei <-  ggarrange(p2.count, p2.agg, p2.peak, p2.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv', 
                     label.x = 0.15, # Adjust x-position (0 is left, 1 is right)
                     label.y = 0.85, # Adjust y-position (1 is top, 0 is bottom)
                     labels = c("[i]", "[j]", "[k]", "[l]"),
                   
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.int <- ggarrange(p3.count,  p3.agg, p3.peak, p3.peak.diff,
                   labels = c("[m]", "[n]", "[o]", "[p]"),
                   #labels = c("[q]", "[r]", "[s]", "[t]"),
                   label.x = 0.15, # Adjust x-position (0 is left, 1 is right)
                   label.y = 0.9, # Adjust y-position (1 is top, 0 is bottom)
                   ncol=4, nrow = 1 , align = 'hv',common.legend = TRUE, legend = 'bottom',
                   font.label = list(size = 8, color = "black", face = "plain", family = NULL))
full_preds <- ggarrange(p_spagetti,
                        #p.previous,
                        p.temp, 
                        p.spei,  
                        p.int,
                        ncol = 1, nrow= 4, 
                        align = 'hv', heights = c(1, 1, 1,  1.1),
                        widths = c(1, 1, 1, 1))
#"#009E73" 
windows(7,9)

(full_preds)
ggsave(filename = 'outFigs/Fig_full_main_effects.png', plot = full_preds, 
       width = 7.5, height = 9, dpi = 300, bg = 'white')










##### test Morans I (from logged beetles numbers) for differences between years --------------- 
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




