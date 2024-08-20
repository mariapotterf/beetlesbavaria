


rm(list=ls()) 
gc()

# Intro -----------------------------------------------------------------------
### Read my paths -----------------------------------------------------------
source('myPaths.R')



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
#predictors_RS         <- c("spring_tmp", "veg_tmp", "spei3", "spei12", "sum_ips", "peak_doy", "peak_diff", "agg_doy", "Morans_I", "Morans_I_log")
#dependent_var <- 'sum_ips'
lags <- 0:3
#data <- dat_fin # Assuming dat_fin is your dataset

  

# dat_fin: calculate climate lags --------------
dat_spei_lags <-  dat_fin %>% 
  dplyr::select(c(year, pairID, trapID, tmp, spei, 
                  tmp_z, #spei_z, 
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
   #     spei_z_lag1 = lag(spei_z, n = 1, default = NA),
  #       spei_z_lag2 = lag(spei_z, n = 2, default = NA),
  #       spei_z_lag3 = lag(spei_z, n = 3, default = NA),
    spei_lag1 = lag(spei, n = 1, default = NA),
    spei_lag2 = lag(spei, n = 2, default = NA),
    spei_lag3 = lag(spei, n = 3, default = NA),
          ) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup(.) 

nrow(dat_spei_lags)


dat_spei_lags %>% 
  dplyr::filter(trapID == "Anzinger_Forst_1") %>% 
  dplyr::select(trapID, year, tmp_z, tmp_z_lag1,tmp_z_lag2 ) # sum_ips, sum_ips_lag1

# check lags if correct. YES!


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

# scale IPS_count data and AGG_doy tables --------
# Specify the columns to skip
columns_to_skip_counts <- c("year", "trapID", "sum_ips", "tr_peak_doy", "peak_diff","peak_doy",  "pairID", "x", "y", "year_fact")
columns_to_skip_agg <- c("year", "trapID", "tr_agg_doy", "tr_peak_doy","agg_doy",  "pairID", "x", "y", "year_fact")


# Scale and center the remaining numeric columns
dat_fin_counts_m_scaled <- dat_fin_counts_m %>%
  mutate(across(
    .cols = where(is.numeric) & !all_of(columns_to_skip_counts),
    .fns = ~ scale(.) %>% as.vector()
  ))


# Scale and center the remaining numeric columns
dat_fin_agg_m_scaled <- dat_fin_agg_m %>%
  mutate(across(
    .cols = where(is.numeric) & !all_of(columns_to_skip_agg),
    .fns = ~ scale(.) %>% as.vector()
  ))



####GAM: Find the best predictor lags  (simple) --------------------------------------------------------------- 
dependent_vars_counts <-  c("sum_ips", "peak_diff")
dependent_vars_doy    <-  c("tr_agg_doy", "tr_peak_doy")


# Initialize a data frame to store AIC values and deviance explained
model_metrics_count <- data.frame(Predictor = character(), 
                            Dependent = character(), 
                            AIC = numeric(), DevianceExplained = numeric())

# List of dependent variables and predictors:
# keep only spei3 and spei 12 - for the short term vs long term effect
selected_predictors <- c(
  #"tmp", "tmp_lag1", "tmp_lag2","tmp_lag3" ,
  "tmp_z", "tmp_z_lag1", "tmp_z_lag2","tmp_z_lag3" ,
  "spei", "spei_lag1", "spei_lag2","spei_lag3" )
  #"spei_z", "spei_z_lag1", "spei_z_lag2","spei_z_lag3") 



# Loop over each dependent variable
for (dep in dependent_vars_counts) {
  #print(dep)
  # Loop over each predictor
  for (pred in selected_predictors) {
    #print(pred)
    # Fit the model
    #formula <- as.formula(paste(dep, "~ s(", pred, ", k = 4) + s(pairID, bs = 're')"))
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 5)"))
    #print(formula)
    model <- gam(formula, family = nb(), method = 'REML', data = dat_fin_counts_m_scaled)
    
    # Extract model summary
    model_summary <- summary(model)
    
    # Store the AIC value and deviance explained
    model_metrics_count <- rbind(model_metrics_count, data.frame(Predictor = pred, Dependent = dep, AIC = AIC(model), DevianceExplained = round(model_summary$dev.expl*100,1)))
  }
}

# View the AIC values and deviance explained
print(model_metrics_count)

# Select the best predictor for each dependent variable based on the lowest AIC
best_predictors_counts <- model_metrics_count %>% 
  mutate(category = case_when(
    grepl("tmp", Predictor  ) ~ "tmp",
    grepl("spei", Predictor  ) ~ "spei",
    TRUE ~ "other"
  )) %>% 
  group_by(Dependent, category) %>% 
  slice_min(AIC, n = 1)   # Select the best 3 based on AIC
 # slice(which.max(DevianceExplained))

best_predictors_counts
####### for DOY values ---------------------------------------------------------------

# Initialize a data frame to store AIC values and deviance explained
model_metrics_doy <- data.frame(Predictor = character(), 
                                  Dependent = character(), 
                                  AIC = numeric(), DevianceExplained = numeric())

# Loop over each dependent variable
for (dep in dependent_vars_doy ) {
  #print(dep)
  # Loop over each predictor
  for (pred in selected_predictors) {
    #print(pred)
    # Fit the model
    #formula <- as.formula(paste(dep, "~ s(", pred, ", k = 4) + s(pairID, bs = 're')"))
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 5)"))
    #print(formula)
    model <- gam(formula, family = betar(link = "logit"), method = 'REML', data = dat_fin_agg_m_scaled)
    
    # Extract model summary
    model_summary <- summary(model)
    
    # Store the AIC value and deviance explained
    model_metrics_doy <- rbind(model_metrics_doy, data.frame(Predictor = pred, Dependent = dep, AIC = AIC(model), DevianceExplained = round(model_summary$dev.expl*100,1)))
  }
}

# View the AIC values and deviance explained
print(model_metrics_doy)

# Select the best predictor for each dependent variable based on the lowest AIC
best_predictors_doy <- model_metrics_doy %>% 
  mutate(category = case_when(
    grepl("tmp", Predictor  ) ~ "tmp",
    grepl("spei", Predictor  ) ~ "spei",
    TRUE ~ "other"
  )) %>% 
  group_by(Dependent, category) %>% 
 
  slice_min(AIC, n = 1) #%>%  # Select the best 3 based on AIC
  #slice(which.min(AIC))

best_predictors_doy


# merge both tables:
full_lag_models <- rbind(model_metrics_count,model_metrics_doy) #%>% 
full_lag_models <- full_lag_models %>% 
  mutate(category = case_when(
    grepl("tmp", Predictor  ) ~ "tmp",
    grepl("spei", Predictor  ) ~ "spei",
    TRUE ~ "other"
  ))  %>% 
  dplyr::filter(category != 'other') 


  
View(full_lag_models)
best_predictors <- rbind(best_predictors_counts, best_predictors_doy )

best_predictors <- best_predictors %>% 
  dplyr::filter(category != 'other') 

View(best_predictors)


#fwrite(full_lag_models, 'outTable/fin_lag_predictors.csv')
sjPlot::tab_df(full_lag_models,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_lag_predictors.doc",
               digits = 1) 

# Calculate the correlation matrix for all predictors and their lags ---------------------------

predictor_vars <- dat_fin_counts_m[, selected_predictors  ] # unique(best_predictors$Predictor)
cor_matrix <- cor(predictor_vars, method = "spearman")
(cor_matrix)
corrplot::corrplot(cor_matrix)

# mask highly correlated values
# Mask high correlations (set them to NA)
cor_matrix_masked <- cor_matrix
cor_matrix_masked[cor_matrix > 0.4 | cor_matrix < -0.4] <- NA
corrplot::corrplot(cor_matrix_masked)


#### adreass lags in different way: provide plot with 0,1,2,3 lag ----------------------
# plot for SPEIs (1,3,6,12)and for veg_tmp
# Define the function
run_models_and_collect_results <- function(data, predictors, response, k = 4, family = nb(), method = 'REML') {
  # Define a function to fit the model and get predictions
  fit_model_and_predict <- function(predictor) {
    # Construct the formula
    formula <- as.formula(paste(response, "~ s(", predictor, ", k = ", k, ")", sep = ""))
    
    # Fit the model
    model <- gam(formula, family = family, method = method, data = data)
    
    # Get predictions using ggpredict
    predictions <- ggpredict(model, terms = paste(predictor, "[all]", sep = ""), allow.new.levels = TRUE)
    
    # Add the predictor name to the results
    predictions$predictor <- predictor
    
    return(predictions)
  }
  
  # Use lapply to apply the function to each predictor
  results_list <- lapply(predictors, fit_model_and_predict)
  
  # Combine all results into a single data frame
  combined_results <- bind_rows(results_list)
  
  return(combined_results)
}

response <- "sum_ips"

results_df <- run_models_and_collect_results(dat_fin_counts_m_scaled, 
                                             selected_predictors  , response)

# Modify the predictors
results_df2 <- 
  results_df %>% 
    as.data.frame() %>%  
  mutate(predictor = str_replace_all(predictor, "veg_tmp", "vegtmp"),
         predictor = str_replace_all(predictor, "spring_tmp", "springtmp"),
         predictor = str_replace_all(predictor, "tmp_z", "tmpz"),
         predictor = str_replace_all(predictor, "spei", "spei")) %>% 
    mutate(predictor_full = case_when(
      str_detect(predictor, "_lag") ~ predictor,
      TRUE ~ paste0(predictor, "_lag0")
    )) %>%
    # Split the strings into their components
  mutate(predictor_split = stringr::str_split(predictor_full, "_")) %>%
  # Unnest the list column into separate columns
  mutate(type = map_chr(predictor_split, 1),
         lag = map_chr(predictor_split, 2)) %>%
  # Drop the temporary columns
  dplyr::select(-predictor_split) 
  


# plot: effect of lags on beetle population level
ggplot(results_df2, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = lag  ), alpha = 0.3) +
  geom_line(aes(color = lag  , linetype = lag  ), linewidth = 1)  +
  facet_grid(lag~type) +
  theme_bw()



# identify the 'hot and dry conditions'
plot(dat_fin_counts_m_scaled$veg_tmp,dat_fin_counts_m_scaled$spei)
plot(dat_fin_counts_m_scaled$tmp_z,dat_fin_counts_m_scaled$spei)

pairs(sum_ips ~tmp_z_lag2 + spei_z_lag2,   data  = dat_fin_counts_m_scaled, panel=panel.smooth)

# identify howt and ry conditions visually
# Load your data
dat <- dat_fin_counts_m_scaled_no_outliers

# Identify the points in the lower right cluster
cluster_points <- dat[dat$tmp_z > 1.5 & dat$spei < 0, ]

# Print or inspect the identified points
summary(cluster_points)

# they are all in 2018: split data between before and after drought?
dat_before <- dat %>% dplyr::filter(year%in% 2015:2017)
dat_after <- dat %>% dplyr::filter(year%in% 2018:2021)

mm <- gam(sum_ips ~s(tmp_z, k = 4), 
          family = nb(), method = 'REML', data = dat[,year == 2018])

plot(mm, page = 1)


# check multicollinearity between SPEI and tmp: -------------------------------
# Function to fit model and calculate VIF
# Function to fit model and check multicollinearity
fit_model_and_check_vif <- function(tmp_var, spei_var, data) {
  formula <- as.formula(paste("sum_ips ~ s(", tmp_var, ", k = 5) + s(", spei_var, ", k = 5)", sep = ""))
  model <- gam(formula, family = nb(), method = 'REML', data = data)
  collinearity <- check_collinearity(model)
  return(collinearity)
}

# Find combinations and check for multicollinearity
results <- data.frame(tmp_var = character(), spei_var = character(), vif_tmp = numeric(), vif_spei = numeric(), stringsAsFactors = FALSE)

tmp_vars <- selected_predictors[str_detect(selected_predictors, "veg_tmp|spring_tmp")]
spei_vars <- selected_predictors[str_detect(selected_predictors, "spei3|spei12")]

for (tmp in tmp_vars) {
  for (spei in spei_vars) {
    collinearity <- fit_model_and_check_vif(tmp, spei, dat_fin_counts_m_scaled)
    vif_values <- collinearity$VIF
    results <- rbind(results, data.frame(tmp_var = tmp, spei_var = spei, vif_tmp = vif_values[1], vif_spei = vif_values[2]))
  }
}


# Filter combinations with acceptable VIF values (e.g., VIF < 5)
acceptable_results <- results %>% filter(vif_tmp < 5 & vif_spei < 5)

(acceptable_results)

# based on this results I have removes the spei 12 from teh study - 
# to avoid multicllinearity with tmp
# spei 1, 6 were removed due to strong non-linear relaionship before





###### check correlation between traps -------------------------- 
# check correlation between trap counts: does it explains anything?
# need to conevert to wide format for correlation
df_wide <- dat_fin_counts_m %>%
  dplyr::select(year, trapID, pairID, sum_ips) %>%
  mutate(trap_number = str_sub(trapID, -1, -1)) %>%
  dplyr::select(-trapID) %>% 
  pivot_wider(
    names_from = trap_number,
    values_from = sum_ips,
    names_prefix = "trap_"
  ) %>% 
  # calculkate difference between two traps
  mutate(abs_trap_diff = trap_1 - trap_2,
         ratio = trap_1/trap_2) # larger number shows higher difference


df_corr <- df_wide %>%
  group_by(pairID) %>%
  mutate(
    cor = cor(trap_1, trap_2, method = 'spearman', use = "complete.obs")  ) 



df_wide %>% 
  na.omit() %>% 
  ggplot(aes(x = trap_1, y = trap_2)) + 
  geom_point() + 
  geom_smooth()


plot(df_wide$trap_1, df_wide$trap_2 )


df_corr %>% 
  ggplot(aes(x = trap_1, y = trap_2)) + 
  geom_point() + 
  geom_smooth(method = 'lm') #+ 
  #facet_wrap(.~pairID, scale = 'free')


df_corr %>% 
  ggplot(aes(x = trap_1, y = trap_2)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.8, lty = 'dashed') +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs", k = 4)) + 
  facet_wrap(.~year, scale = 'free')



#data_filtered <- data %>% 
 # na.omit() # This removes rows with any NAs across all columns

# check outliers in beetle counts
boxplot(dat_fin_counts_m$sum_ips)
boxplot(log(dat_fin_counts_m$sum_ips))




# Calculate average counts for each trap pair ------------------------------------
avg_data <- dat_spei_lags %>%
  group_by(pairID, year) %>%
  summarise(Morans_I_log  = mean(Morans_I_log , na.rm = TRUE),
            sum_ips     = mean(sum_ips, na.rm = TRUE),
            peak_diff   = mean(peak_diff, na.rm = TRUE),
            tr_agg_doy  = mean(tr_agg_doy, na.rm = TRUE),
            tr_peak_doy = mean(tr_peak_doy, na.rm = TRUE), 
            agg_doy     = mean(agg_doy, na.rm = TRUE),
            peak_doy    = mean(peak_doy, na.rm = TRUE),
            spei        = mean(spei, na.rm = TRUE),
            tmp_z       = mean(tmp_z, na.rm = TRUE),
            spei_lag2 = mean(spei_lag2, na.rm = TRUE),
            tmp_z_lag1  = mean(tmp_z_lag1, na.rm = TRUE),
           # tmp_lag1  = mean(tmp_lag1, na.rm = TRUE),
            spei_lag1 = mean(spei_lag1, na.rm = TRUE), # for agg_doy
            tmp_z_lag2  = mean(tmp_z_lag2, na.rm = TRUE),  # for agg_doy
          #  tmp_lag2  = mean(tmp_lag2, na.rm = TRUE),  # for agg_doy
            x           = mean(x, na.rm = TRUE),
            y           = mean(y, na.rm = TRUE),
            ) %>%
  ungroup(.) %>% 
  na.omit() %>% 
  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% c('Eschenbach_idOPf', 'Peiting') ) # remove two traps with always having too low trap counts

fwrite(avg_data, 'outTable/fin_tab_avg.csv')


### automate models over all dependent variables and over increasing model complexity  -----------------------------------------
# START 

# Define a function to build and compare models with descriptive names for multiple dependent variables
compare_models <- function(data, dependent_vars) {
  results <- list()
  
  for (dep_var in dependent_vars) {
    models <- list()
    AIC_values <- data.frame(Model = character(), AIC = numeric(), stringsAsFactors = FALSE)
    
    # Determine the appropriate family for the dependent variable
    if (dep_var %in% c("sum_ips", "peak_diff")) {
      family <- nb()
    } else if (dep_var %in% c("tr_agg_doy", "tr_peak_doy")) {
      family <- betar()
    } else {
      stop("Unknown dependent variable family.")
    }
    
    # Function to safely fit a model and catch errors
    safe_gamm <- function(formula, data, family, correlation) {
      tryCatch({
        gamm(formula, data = data, family = family, correlation = correlation)
      }, error = function(e) {
        cat("Error in model fitting:", conditionMessage(e), "\n")
        return(NULL)
      })
    }
    
    # Step 1: Initial Model (Base Model with Year)
    models$base_model <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6)")),
                                   data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$base_model)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Base Model (Year)", AIC = AIC(models$base_model$lme)))
    }
    
    # Step 2: Add s(tmp_z_lag1, k = 8)
    models$tmp_lag1 <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6) + s(tmp_z_lag1, k = 8)")),
                                 data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$tmp_lag1)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Add tmp_z_lag1", AIC = AIC(models$tmp_lag1$lme)))
    }
    
    # Step 3: Add s(spei_z_lag2, k = 4)
    models$spei_lag2 <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4)")),
                                  data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$spei_lag2)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Add spei_z_lag2", AIC = AIC(models$spei_lag2$lme)))
    }
    
    # Step 4: Add te(tmp_z_lag1, spei_z_lag2, k = 15)
    models$tmp_spei_interaction <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) + te(tmp_z_lag1, spei_z_lag2, k = 15)")),
                                             data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$tmp_spei_interaction)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Add tmp_spei_interaction", AIC = AIC(models$tmp_spei_interaction$lme)))
    }
    
    # Step 5: Add s(x, y, bs = "gp")
    models$spatial <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) + te(tmp_z_lag1, spei_z_lag2, k = 15) + s(x, y, bs = 'gp')")),
                                data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$spatial)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Add spatial", AIC = AIC(models$spatial$lme)))
    }
    
    # Step 6: Add s(pairID, bs = "re")
    models$random_effect <- safe_gamm(as.formula(paste(dep_var, "~ s(year, k = 6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) + te(tmp_z_lag1, spei_z_lag2, k = 15) + s(x, y, bs = 'gp') + s(pairID, bs = 're')")),
                                      data, family, corAR1(form = ~ year | pairID))
    
    if (!is.null(models$random_effect)) {
      AIC_values <- rbind(AIC_values, data.frame(Model = "Add random effect", AIC = AIC(models$random_effect$lme)))
    }
    
    # Store results for each dependent variable
    results[[dep_var]] <- list(models = models, AIC_values = AIC_values)
    
    # Print AIC values for each dependent variable
    cat("\nAIC values for", dep_var, ":\n")
    print(AIC_values)
  }
  
  # Return the results for further inspection if needed
  return(results)
}

# Define the list of dependent variables
dependent_vars <- c("sum_ips", "peak_diff", "tr_agg_doy", "tr_peak_doy")
dependent_vars  <- c("peak_diff")
# Run the function on your dataset
result <- compare_models(avg_data, dependent_vars)
result_peak_diff <- compare_models(avg_data, dependent_vars =  c("peak_diff"))


### Specify temp autocorrelation ------------------------------------------------

###### AGG DOY ------------------------------------

m.gamma <- gamm(tr_agg_doy ~ s(year, k = 7) + 
             s(tmp_z_lag2, k = 10)  + #
             s(spei_z_lag1, k = 3),# + 
           data = avg_data_filt, 
           family = Gamma(link = "log"),
           correlation = corAR1(form = ~ year | pairID))

appraise(m.gamma$gam)


# gamme is better then quasi and betar!
# filter extreme values: occuring too late in a year ()
avg_data_agg <- avg_data %>% 
  dplyr::filter(tr_agg_doy < 0.54)

# m1.agg.gamma one converges well!! - but seems very wiggly in teh iteraction
m1.agg.gamma <- gamm(tr_agg_doy ~ s(year, k = 6) + 
                       s(tmp_z_lag2, k = 7) + 
                       s(spei_z_lag1, k = 5) + 
             te(tmp_z_lag2, spei_z_lag1, k = 20) + 
             #s(x, y, bs = 'gp', k = 15) + 
             s(pairID, bs = 're'),
             data = avg_data_agg , 
           family = Gamma(link = "log"),
           correlation = corAR1(form = ~ year | pairID))

# test simpler model, change spline type: works well without xy!
m2.agg.gamma <- gamm(tr_agg_doy ~ s(year, k = 6) + 
                       s(tmp_z_lag2, k = 3,bs = "cs") + 
                       s(spei_z_lag1, k = 3) + 
                       te(tmp_z_lag2, spei_z_lag1, k = 8, bs = "cs") + 
                       #s(x, y, bs = 'gp', k = 30) + 
                       s(pairID, bs = 're'),
                     data = avg_data_agg , 
                     family = Gamma(link = "log"),
                     correlation = corAR1(form = ~ year | pairID))

AIC(m2.agg.gamma$lme, m1.agg.gamma$lme)
#boxplot(avg_data_agg$tr_agg_doy)
summary(m2.agg.gamma$gam)
appraise(m2.agg.gamma$gam)
k.check(m2.agg.gamma$gam)
gam.check(m2.agg.gamma$gam)
plot(m2.agg.gamma$gam, page = 1)


###### PEAK DOY ------------------------------------
m.peak.gamma <- gamm(tr_peak_doy ~ s(year, k = 6) + s(tmp_z_lag1, k = 5) + 
                       s(spei_z_lag2, k = 5) + 
                   te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                   s(x, y, bs = 'gp', k = 50) + 
                   s(pairID, bs = 're'),
                   data = avg_data_filt, 
                 family = Gamma(link = "log"),
                 correlation = corAR1(form = ~ year | pairID))

appraise(m.peak.gamma$gam)
plot(m.peak.gamma$gam, page = 1)
k.check(m.peak.gamma$gam)
gam.check(m.peak.gamma$gam)
summary(m.peak.gamma$gam)

###### PEAK DIFF TW  ------------------------------------

m.peak.diff <- gamm(peak_diff ~ s(year, k = 6) + s(tmp_z_lag1, k = 5) + s(spei_z_lag2, k = 5) + 
                   te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                   s(x, y, bs = 'gp', k = 15) + 
                   s(pairID, bs = 're'),data = avg_data, 
                 family = nb, #Gamma(link = "log"),
                 #family = betar, 
                 correlation = corAR1(form = ~ year | pairID))


# remove outliers: Piding
avg_data_peak_diff <- avg_data %>% 
  dplyr::filter(!pairID %in% c('Piding', 'Gangkofen'))

m.peak.diff.tw <- gamm(peak_diff ~ s(year, k = 6) +
                         s(tmp_z_lag1, k = 5) +
                         s(spei_z_lag2, k = 5) + 
                      te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      s(x, y, bs = 'gp', k = 50) + 
                      s(pairID, bs = 're'),
                        data = avg_data_peak_diff, 
                      family = tw,
                      correlation = corAR1(form = ~ year | pairID))

appraise(m.peak.diff.tw$gam)
summary(m.peak.diff.tw$gam)
k.check(m.peak.diff.tw$gam)
gam.check(m.peak.diff.tw$gam)
plot(m.peak.diff.tw$gam, page = 1)




###### SUM IPS TW  ------------------------------------

###  Use previous year ------------------------------------------------

avg_data_filt <- avg_data %>% 
  mutate(f_year = as.factor(year)) %>% 
  dplyr::filter(!pairID %in% pair_outliers )


nrow(avg_data_filt) # 518
nrow(avg_data) # 539

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


avg_data_filt_lagged %>% 
  dplyr::filter(pairID == 'Anzinger_Forst') %>% 
  dplyr::select(year, sum_ips, sum_ips_lag1, 
                tr_peak_doy, tr_peak_doy_lag1 )


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





vif_data <- avg_data_peak_diff[, c("year", "tmp_z_lag1", "spei_z_lag2", "x", "y", "peak_diff")]
vif_result <- vif(lm(peak_diff ~ ., data = vif_data))


appraise(fin.m.agg)

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


# keep simple one for plottig purposes
# m.peak.gamma.simple <- gamm(tr_peak_doy ~ 
#                               s(tmp_z_lag2, k = 3)  + # 
#                               s(spei_lag1, k = 3) + 
#                               te(tmp_z_lag2, spei_lag1, k = 5, bs = 'tp') + #
#                               #s(x, y, bs = 'gp', k = 10) + 
#                               s(pairID, bs = 're')# +
#                             # s(f_year, bs = 're')
#                             ,
#                             data =  avg_data_filt_peak_doy,#avg_data_filt , 
#                             family = Gamma(link = "log"),
#                             correlation = corAR1(form = ~ year | pairID),
#                             control = control)
# 
# summary(m.peak.gamma.simple$gam)
# appraise(m.peak.gamma.simple$gam)
# k.check(m.peak.gamma.simple$gam)
# plot(m.peak.gamma.simple$gam, page = 1, shade = TRUE)
# 


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
                       data = avg_data_peak_no_out, ,
                       family = Gamma(link = "log"),
                       control = control)

summary(m.peak.previous.gamm$gam)
appraise(m.peak.previous.gamm$gam)
k.check(m.peak.previous.gamm$gam)
plot(m.peak.previous.gamm$gam, page = 1, shade = TRUE)


# final model: m.peak.previous.gamm

AIC(m.peak.previous.year, m.peak.previous)



###### PREVIOUS get final models ----------------
fin.m.counts.previous.tw    <- m.counts.previous$gam
fin.m.peak.diff.previous.tw <- m.peak.diff.previous$gam
fin.m.agg.doy.gamma         <- m.agg.previous$gam
fin.m.peak.doy.gamma        <- m.peak.previous.gamm$gam






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

# list predictors to test
selected_predictors <- c('sum_ips', 'sum_ips_lag1','sum_ips_lag2',
                         #'log_sum_ips', 'log_sum_ips_lag1','sum_ips_lag2',
                         'tr_agg_doy'   , 'tr_agg_doy_lag1', 'tr_agg_doy_lag2' ,
                         'tr_peak_doy', 'tr_peak_doy_lag1','tr_peak_doy_lag2',
                         'peak_diff', 'peak_diff_lag1',  'peak_diff_lag2',
                        'spei',  'spei_lag1','spei_lag2',
                        'tmp_z',  'tmp_z_lag1', 'tmp_z_lag2'  
) 


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
    #print(formula)
   # model <- gamm(formula,  
  #                family = tw,
  #                method = 'REML',  
  #                data = lisa_merged_df_avg,#avg_data_moran_sub,
  #                correlation = corAR1(form = ~ year | pairID))
    
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

sjPlot::tab_model(fin.m.moran,
              #col.header = c(as.character(qntils), 'mean'),
             # show.rownames = FALSE,
              file="outTable/model_moran.doc") 



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
fin.m <- fin.m.moran #m.agg.previous$gam # m2.agg.gamma$gam #m.counts.tw$gam # m.counts.tw.bam.year#m.counts.tw.test.f.rndm#  fin.m.agg#$gam

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


save(fin.m.counts.previous.tw,
     fin.m.peak.diff.previous.tw,
     fin.m.agg.doy.gamma,
     fin.m.peak.doy.gamma,
         avg_data,
     file = "outData/fin_models.RData")









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

my_theme_square <- function() {
  theme_minimal(base_size = 8) +
    theme(
      #aspect.ratio = 1, 
      axis.ticks.y = element_line(),
      axis.ticks.x = element_line(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.position = "bottom",
      axis.title.y = element_blank()
    ) 
}






# Create effect plot function with an additional argument to control y-axis labels
create_effect_previous_year <- function(data, avg_data, line_color = "grey40", 
                               x_col = "sum_ips_lag1", 
                               y_col = "sum_ips", 
                               x_title = "X-axis", 
                               y_title = "Y-axis", 
                               my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() + 
    geom_point(data = avg_data, aes(x = !!x_col, y = !!y_col, group = pairID), col = "gray60", alpha = 0.4) +
    
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color) +
    
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}


# Create effect plot function with additional arguments to select columns from avg_data
create_effect_plot <- function(data, 
                               avg_data, 
                               x_col = "tmp_z_lag1", 
                               y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = "X-axis", 
                               y_title = "Y-axis", my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() +
    geom_point(data = avg_data, aes(x = !!x_col, y = !!y_col), col = "gray60", alpha = 0.3) +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color) +
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}



# plot interactions
plot_effect_interactions <- function(data, 
                                     avg_data, 
                                     x_col = "tmp_z_lag1", 
                                     y_col = "sum_ips", 
                                     temp_label = 'Temp',
                                     y_title = 'y_title',
                                     x_annotate = 0, 
                                     lab_annotate = "lab ann") {
  #library(ggplot2)
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
 p<- ggplot() +
     geom_point(data = avg_data, aes(x = !!x_col, y = !!y_col), col = "gray60", alpha = 0.4) +
    geom_ribbon(data=data, aes(x = x,
                               ymin = conf.low, 
                               ymax = conf.high, 
                               fill = group), alpha = 0.25) +
   geom_line(data = data, 
             aes(x = x, y = predicted,
                 color = group, linetype = group), linewidth = 1) +
    labs(x = temp_label,
         y = y_title,
         fill = "SPEI levels",
         color = "SPEI levels",
         linetype = "SPEI levels") +  # Fixed "y_title" to "y" for correct y-axis label argument
  
   guides(color = guide_legend(nrow = 1), 
           fill = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1)) +
    annotate("text", x = x_annotate, y = Inf, 
             label = lab_annotate, hjust = 0.5, vjust = 1.5) +
   my_theme_square() +
   theme(legend.position = "bottom")
  
  return(p)
}

plot_effect_interactions(p3, avg_data = avg_data_filt_lagged_plotting)


##### Spagetting plot: development over time -----------------------------------------




plot_data_with_average <- function(data, y_var, y_label, my_title) {
  # Convert the y_var argument to a symbol to use in aes()
  y_var_sym <- rlang::sym(y_var)
  
  data %>%
    ungroup() %>%
    filter(year %in% 2015:2021) %>%
    ggplot(aes(x = year, y = !!y_var_sym, group = pairID)) +
    labs(x = 'Year', 
         y = y_label, 
         title = my_title) +
    geom_line(alpha = 0.1) +  
    stat_summary(
      aes(x = year, y = !!y_var_sym, group = 1), 
      fun = mean,  # Calculate the mean for each year
      geom = "line", 
      color = "#E69F00",  # Ensure the average line is also red
      linewidth = 1  # Make the average line slightly thicker than individual lines
    ) + my_theme_square()
   
}

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


p0.count <- create_effect_previous_year(data = p0, 
                              avg_data = avg_data_filt_lagged_plotting,
                              x_col = "sum_ips_lag1",
                              y_col = "sum_ips",
                              line_color = "darkgreen", 
                              x_title = "Pop. level [*1000, lag1]", 
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
                     x_title = expression(paste("Temp. lag1 [", "z-score]")),
                    # x_title = temp_label , 
                     y_title = paste(lab_popul_level, '*1000'), 
                    # my_title = paste("[a]", 'Population level', '\n[#*1000]'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "**")

(p1.count)
p2.count <- create_effect_plot(data = p2, avg_data = avg_data_filt_lagged_plotting, 
                               x_col = "spei_lag2", y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI. lag2 [", "dim.]")),#spei_label , 
                               y_title = paste(lab_popul_level, '*1000'), 
                               #my_title = "Effect of Year on Sum IPS", 
                               x_annotate = -0.5, 
                               lab_annotate = "***")
(p2.count)
p3.count <- plot_effect_interactions(data = p3, 
                                     avg_data = avg_data_filt_lagged_plotting, 
                                     x_col = "tmp_z_lag1", 
                                     y_col = "sum_ips", 
                                     temp_label = expression(paste("Temp. lag1 [", "z-score]")), 
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
                               x_title = "Aggregation timing[DOY, lag1]",
                               y_title = lab_colonization_time,
                              # my_title = paste("[b]", 'Aggregation timing', '\n[DOY]'),
                               x_annotate = 175, lab_annotate = "***")

(p0.agg)
p1.agg <- 
  create_effect_plot(data = p1,  avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),  
                     x_col = "tmp_z_lag1", y_col = "agg_doy", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag1 [", "z-score]")) , 
                     y_title = lab_colonization_time, 
                     #my_title = paste("[b]", 'Aggregation timing', '\n[DOY]'), 
                     # my_title = paste("[b]", 'Aggregation timing', '\n[DOY]'),
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "***")

(p1.agg)
p2.agg <- create_effect_plot(data = p2, avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),    
                               x_col = "spei_lag2", y_col = "agg_doy", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI. lag2 [", "dim.]")) , 
                               y_title = lab_colonization_time, 
                               x_annotate = -0.5, 
                               lab_annotate = "0.49")
(p2.agg)
p3.agg <- plot_effect_interactions(p3,
                                   avg_data = filter(avg_data_filt_lagged_plotting, agg_doy < 200),
                                   x_col = "tmp_z_lag1", 
                                   y_col = "agg_doy", 
                                     temp_label = expression(paste("Temp. lag1 [", "z-score]")), 
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
                                      #avg_data = avg_data_filt_lagged_plotting,
                                      x_col = "tr_peak_doy_lag1_plot",
                                      y_col = "peak_doy",
                                      line_color = "darkgreen",
                                      x_title = "Peak swarming timing [DOY, lag1]",
                                      y_title = lab_colonization_time,
                                      #my_title = paste("[b]", 'Peak swarming timing', '\n[DOY]'),
                                      x_annotate = 160, lab_annotate = "***")

(p0.peak)
p1.peak <- 
  create_effect_plot(data = p1,
                     avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                     x_col = "tmp_z_lag2", y_col = "peak_doy", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag2 [", "z-score]")) , 
                     y_title = lab_peak_time, 
                       x_annotate = 2, 
                     lab_annotate = "0.71")

(p1.peak)
p2.peak <- create_effect_plot(data = p2,  avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220), 
                             x_col = "spei_lag1", y_col = "peak_doy", 
                             line_color = "blue", 
                             x_title = expression(paste("SPEI. lag1 [", "dim.]")) , 
                             y_title = lab_peak_time, 
                             x_annotate = -0.5, 
                             lab_annotate = "**")
(p2.peak)
p3.peak <- plot_effect_interactions(p3, 
                                    avg_data = filter(avg_data_filt_lagged_plotting, tr_peak_doy_lag1_plot  < 220),
                                    x_col = "tmp_z_lag1", 
                                    y_col = "peak_doy", 
                                   temp_label = expression(paste("Temp. lag1 [", "z-score]")), 
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
                               #avg_data = avg_data_filt_lagged_plotting, #dplyr::filter(avg_data_sum_ips, peak_diff <200),
                               x_col = "peak_diff_lag1",
                               y_col = "peak_diff",
                               line_color = "darkgreen", 
                               x_title = "Peak swarming intensity [*10, lag1]", 
                               y_title = lab_peak_growth,
                               #my_title = paste("[d]", 'Peak swarming\nintensity', '[#*10]'),  
                               x_annotate = 90, lab_annotate = "***")

(p0.peak.diff)
p1.peak.diff <- 
  create_effect_plot(data = p1, 
                     avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4), x_col = "tmp_z_lag1", y_col = "peak_diff", 
                     line_color = "red", 
                     x_title = expression(paste("Temp. lag2 [", "z-score]")) , 
                     y_title = lab_peak_growth, 
                    # my_title = paste("[d]", 'Peak swarming\nintensity', '[#*10]'),
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "0.14")

(p1.peak.diff)
p2.peak.diff <- create_effect_plot(data = p2,
                                  # avg_data =  avg_data_filt_lagged_plotting,#dplyr::filter(avg_data_sum_ips, peak_diff <200), 
                                  avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4),
                                  x_col = "spei_lag2", y_col = "peak_diff", 
                               line_color = "blue", 
                               x_title = expression(paste("SPEI. lag2 [", "dim.]")) , 
                               y_title = lab_peak_growth, 
                               #my_title = "Effect of Year on Sum IPS", 
                               x_annotate = -0.5, 
                               lab_annotate = "0.14")
(p2.peak.diff)
p3.peak.diff <- plot_effect_interactions(p3,
                                         avg_data = filter(avg_data_filt_lagged_plotting, tmp_z_lag1  < 4),
                                         #avg_data = avg_data_filt_lagged_plotting, # avg_data_sum_ips, 
                                         x_col = "tmp_z_lag1", 
                                         y_col = "peak_diff", 
                                     temp_label = expression(paste("Temp. lag1 [", "z-score]")), 
                                     y_title = lab_peak_growth,
                                     x_annotate = 2,
                                     lab_annotate = "*") 

(p3.peak.diff)
ggarrange(p0.peak.diff,p1.peak.diff,p2.peak.diff,p3.peak.diff)


# put in the same figure: tmp, spei, interaction, previous years on teh bottom

p.previous <- ggarrange(p0.count, p0.agg, p0.peak, 
                         p0.peak.diff, 
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








# PLOT Moran's I -------------------------------
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


# test statistical differences and add to plot --------------- 
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


# Export models summaries --------------------------------------------------

# get characteristics
r2(fin.m.counts )
r2(fin.m.agg)
r2(fin.m.peak)
r2(fin.m.peak.diff)



summary(fin.m.counts )
summary(fin.m.agg)
summary(fin.m.peak)
summary(fin.m.peak.diff)
#summary(fin.m.moran)
#summary(fin.m.RS)



# print models outputs: ---------------------------------------------------------



# export all models
sjPlot::tab_model(fin.m.counts.previous.tw,    file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg.doy.gamma,       file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak.doy.gamma,      file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff.previous.tw, file = "outTable/model_peak_diff.doc")
sjPlot::tab_model(fin.m.moran,     file = "outTable/model_moran.doc")
#sjPlot::tab_model(fin.m.RS,        file = "outTable/model_RS_3.doc")

sum(dat_fin$sum_ips, na.rm = T)



# add significance ----------------------------------------------------
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




