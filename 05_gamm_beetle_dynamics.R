


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
#library(ggpmisc)  # add equation to plots smooth 


# Stats
library('here')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
library('cowplot')
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

# test for autocorrelation
library(lmtest)

# colors
library(RColorBrewer)

library(knitr)   # for table outputs


# plot labels: 
lab_popul_level       = "Population level [#]"
lab_colonization_time = "Aggregation timing [DOY]"
lab_peak_time         = "Peak swarming\ntiming [DOY]"
lab_peak_growth       = "Peak swarming\nintensity [#]"




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
  




# calculate population growth
dat_fin <- dat_fin %>% 
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(peak_diff = as.integer(round(peak_diff))) %>% 
 # mutate(sum_ips_lag1          = lag(sum_ips, n = 1, default = NA),
#         population_growth     = (sum_ips - sum_ips_lag1) / sum_ips_lag1 * 100,
#         population_growth2    = dplyr::lag(population_growth, n = 1, order_by = year)) %>%  # lag population growth by one more year
 ungroup(.) 

# Merge beetle dynamics inicators with XYs
dat_fin <- 
  dat_fin %>% 
  left_join(xy_df, by = c("trapID", 'year')) %>% 
  # remove years woith NAs:
  dplyr::filter(year %in% 2012:2021) 

nrow(dat_fin)



# Find teh most influential time lag :  0,1,2,3----------------------------------------------


# prepare tables

# specify individual datasets, to limit the NAs values and to scale them
dat_fin_no_RS <- dat_fin %>% 
  dplyr::select(-c( wind_beetle, harvest,
                    Morans_I_log, Morans_I, 
                    sum_ips_lag1, population_growth, population_growth2))

cols_skip <- c('trapID', 'pairID', 'Morans_I_log', 'Morans_I')

# export new df
dat_fin_Moran <-  dat_fin %>% 
  #  dplyr::filter(Morans_I_log > 0 & year %in% 2018:2021) %>%
  dplyr::select(-wind_beetle, - harvest)

dat_fin_Moran_scaled <-
  dat_fin_Moran %>% 
  ungroup(.) %>% 
  # na.omit() %>% 
  dplyr::select(-all_of(cols_skip )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_fin_Moran %>% dplyr::select(all_of(cols_skip)), .)




# Make first function for each model type: nb, beta, tweedie 
# make sure that they have the same number of observations!


fit_lags_negbin <- function(data, predictors, lags, dependent_var) {
  results <- list() # Initialize an empty list to store results
  
  # Create all lagged variables within the same table
  for (predictor in predictors) {
    for (lag in lags) {
      lagged_var_name <- paste0("lag", lag, "_", predictor)
      data[[lagged_var_name]] <- dplyr::lag(data[[predictor]], n = lag, default = NA)
    }
  }
  
  # Remove rows with any NA values in the newly created lag columns to ensure equal observation count across models
  data_complete_cases <- data %>% 
    tidyr::drop_na()
  
  # Initialize a variable to store AIC for lag0 model for comparison
  aic_lag0 <- numeric(length(predictors))
  
  for (predictor_idx in seq_along(predictors)) {
    print('1')
    predictor <- predictors[predictor_idx]
    
    # Fit and store AIC for lag0 model first to use as a baseline for comparison
    lagged_var_name_lag0 <- paste0("lag0_", predictor)
    formula_lag0 <- as.formula(paste(dependent_var, "~", lagged_var_name_lag0))
    model_lag0 <- glm.nb(formula_lag0, data = data_complete_cases, na.action = na.exclude)
    aic_lag0[predictor_idx] <- AIC(model_lag0)
    
    for (lag in lags) {
      print('lags')
      lagged_var_name <- paste0("lag", lag, "_", predictor)
      formula <- as.formula(paste(dependent_var, "~", lagged_var_name))
      model <- glm.nb(formula, data = data_complete_cases, na.action = na.exclude)
      aic_value <- AIC(model)
      print(model)
      print(nrow(data_complete_cases))
      df_rows <- nrow(data_complete_cases)
      
      # Calculate ΔAIC
      delta_aic <- aic_value - aic_lag0[predictor_idx]
      
      # Extract p-value for the lagged variable (adjust if necessary for your model structure)
      p_value <- coef(summary(model))[lagged_var_name, "Pr(>|z|)"]
      
      # Compile information into results list
      results[[paste(predictor, "Lag", lag)]] <- data.frame(
        Predictor = predictor,
        Dependent = dependent_var,
        Lag = lag,
        AIC = aic_value,
        P_Value = p_value,
        Delta_AIC = delta_aic,
        df_rows= df_rows
      )
    }
  }
  
  # Combine results into a single data frame
  do.call(rbind, results)
}

fit_lags_beta <- function(data, predictors, lags, dependent_var) {
  results <- list() # Initialize an empty list to store results
  
  data <- data %>%
    group_by(trapID) %>%
    arrange(year, .by_group = TRUE)
  
  # Create all lagged variables in one go to keep the observation count uniform across models
  for (predictor in predictors) {
    for (lag in lags) {
      lagged_var_name <- paste0(predictor, "_lag", lag)
      data[[lagged_var_name]] <- lag(data[[predictor]], n = lag, default = NA)
    }
  }
  
  # Filter out NA values after creating lagged variables
  data <- data %>% ungroup() %>% na.omit()
  
  # Initialize variables to store AIC for comparison
  aic_baseline <- numeric(length(predictors))
  
  for (predictor_idx in seq_along(predictors)) {
    predictor <- predictors[predictor_idx]
    
    # Fit model for lag0 and store its AIC for baseline comparison
    model_lag0_formula <- as.formula(paste(dependent_var, "~", paste0(predictor, "_lag0")))
    model_lag0 <- glmmTMB(model_lag0_formula, family = beta_family(link = "logit"), data = data, 
                          na.action = na.exclude)
    aic_baseline[predictor_idx] <- AIC(model_lag0)
    
    for (lag in lags) {
      lagged_var_name <- paste0(predictor, "_lag", lag)
      model_formula <- as.formula(paste(dependent_var, "~", lagged_var_name))
      
      model <- glmmTMB(model_formula, family = beta_family(link = "logit"), data = data, 
                       na.action = na.exclude)
      
      # Extract p-value for the lagged variable
      summary_model <- summary(model)
      p_value <- summary_model$coefficients$cond[lagged_var_name, "Pr(>|z|)"]
      
      aic_value <- AIC(model)
      delta_aic <- aic_value - aic_baseline[predictor_idx]
      
      # Compile the results
      results[[paste(predictor, "Lag", lag)]] <- data.frame(
        Predictor = predictor,
        Dependent = dependent_var,
        #Model = paste(predictor, "Lag", lag),
        
        Lag = lag,
        AIC = aic_value,
        P_Value = p_value,
        Delta_AIC = delta_aic
      )
    }
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}

fit_lags_tweedie <- function(data, predictors, lags, dependent_var) {
  results <- list() # Initialize an empty list to store results
  
  # Create all lagged variables within the same dataframe
  data <- data %>%
    group_by(trapID) %>%
    arrange(year, .by_group = TRUE)
  
  for (predictor in predictors) {
    for (lag in lags) {
      lagged_var_name <- paste0(predictor, "_lag", lag)
      data[[lagged_var_name]] <- lag(data[[predictor]], n = lag, default = NA)
    }
  }
  
  # Filter out NA values after creating lagged variables
  data <- data %>% ungroup() %>% na.omit()
  
  # Initialize a variable to store AIC for lag0 model for comparison
  aic_baseline <- numeric(length(predictors))
  
  for (predictor_idx in seq_along(predictors)) {
    predictor <- predictors[predictor_idx]
    
    # Fit model for lag0 and store its AIC for baseline comparison
    lagged_var_name_lag0 <- paste0(predictor, "_lag0")
    formula_lag0 <- as.formula(paste(dependent_var, "~", lagged_var_name_lag0))
    model_lag0 <- glmmTMB(formula_lag0, family = tweedie(link = "log"), data = data, 
                          na.action = na.exclude)
    aic_baseline[predictor_idx] <- AIC(model_lag0)
    
    for (lag in lags) {
      lagged_var_name <- paste0(predictor, "_lag", lag)
      formula <- as.formula(paste(dependent_var, "~", lagged_var_name))
      
      model <- glmmTMB(formula, family = tweedie(link = "log"), data = data, 
                       na.action = na.exclude)
      
      # Extract p-value for the lagged variable
      summary_model <- summary(model)
      p_value <- summary_model$coefficients$cond[lagged_var_name, "Pr(>|z|)"]
      
      aic_value <- AIC(model)
      delta_aic <- aic_value - aic_baseline[predictor_idx]
      
      # Compile the results
      results[[paste(predictor, "Lag", lag)]] <- data.frame(
        
        Predictor = predictor,
        Dependent = dependent_var,
        Lag = lag,
        AIC = aic_value,
        P_Value = p_value,
        Delta_AIC = delta_aic
      )
    }
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}


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

  
### Get lag models, validated by AIC --------------------------------------------
dat_fin_no_RS_no_agg <- dat_fin_no_RS %>% 
  dplyr::select(-tr_agg_doy ,-agg_doy )

model_lag_sum_ips <- fit_lags_negbin(dat_fin_no_RS_no_agg, 
                                     predictors  = predictors_beetle_dyn, 
                                     lags, 
                                     dependent_var = "sum_ips")

model_lag_sum_ips_spei <- fit_lags_negbin(dat_fin_no_RS, predictors  = predictors_spei, 
                                     lags, 
                                     dependent_var = "sum_ips")
# spei 3, lag 2 is the best among lags
model_lag_peak_diff <- fit_lags_negbin(dat_fin_no_RS, 
                                        predictors  = predictors_beetle_dyn, 
                                        lags, 
                                        dependent_var = "peak_diff")
model_lag_aggDOY <- fit_lags_beta(dat_fin_no_RS, 
                                   predictors  = predictors_beetle_dyn, 
                                   lags, 
                                   dependent_var = "tr_agg_doy")
model_lag_peakDOY <- fit_lags_beta(dat_fin_no_RS, 
                                   predictors  = predictors_beetle_dyn, 
                                   lags, 
                                   dependent_var = "tr_peak_doy")

model_lag_RS  <- fit_lags_negbin(dat_fin, 
                                 predictors  = predictors_RS, 
                                 lags, 
                                 dependent_var = "wind_beetle")


model_lag_Morans_log  <- fit_lags_tweedie(dat_fin_Moran_scaled, 
                                 predictors  =  predictors_Morans, 
                                 lags, 
                                 dependent_var = "Morans_I_log")

model_lag_Morans_old  <- fit_lags_tweedie(dat_fin_Moran_scaled, 
                                           predictors  =  predictors_Morans, 
                                           lags, 
                                           dependent_var = "Morans_I")



##### simplify lags stats & export tables ----------------------------------------------------
model_info_counts <- bind_rows(model_lag_sum_ips,
                    model_lag_peakDOY,
                    model_lag_aggDOY,
                    model_lag_peak_diff
                     ) %>% 
  mutate(P_Value = round(P_Value     ,3),
         Delta_AIC = round(Delta_AIC     ,3))

(model_info_counts)

model_lag_Morans <- model_lag_Morans %>% 
  mutate(P_Value = round(P_Value     ,3),
         Delta_AIC = round(Delta_AIC     ,3))


model_lag_RS <- model_lag_RS %>% 
  mutate(P_Value = round(P_Value     ,3),
         Delta_AIC = round(Delta_AIC     ,3))



# for beetle counts: the lag of limatic predictor by 2 years
# for RS and Morans: lag of 3 years

##### Export as a nice table in word: ----------------------------------------------
sjPlot::tab_df(model_info_counts,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/lags_ips_counts.doc",
               digits = 3) 

sjPlot::tab_df(model_lag_Morans,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/lags_morans.doc",
               digits = 3) 

sjPlot::tab_df(model_lag_RS,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/lags_RS.doc",
               digits = 3) 











# Run the models based on teh most important lags!! -------------------------------
# follow rules: test for intividual effects, poly tersm and interactions of teh most 
# imortant redictors
# for the most important model
# report


# IPS_counts test based on teh most important lags: lag 2 for spei3 and veg_tmp --------------
dat_spei_lags <-  dat_fin %>% 
  dplyr::select(c(year, pairID, trapID, tmp, spei, 
                  tmp_z, spei_z, 
                  sum_ips, tr_agg_doy, tr_peak_doy,peak_diff,
                  agg_doy, peak_doy, x, y)) %>% 
  group_by(trapID) %>%
  mutate(trapID = as.factor(trapID),
         pairID = as.factor(pairID),
         year_fact = as.factor(year)) %>% 
  arrange(year, .by_group = TRUE) %>%
  # mutate(sum_ips_log = log(sum_ips +1)) %>% 
  mutate(#sum_ips_lag1 = lag(sum_ips_log, n = 1, default = NA),
         tmp_lag1 = lag(tmp, n = 1, default = NA),
         tmp_lag2 = lag(tmp, n = 2, default = NA),
         tmp_lag3 = lag(tmp, n = 3, default = NA),
    tmp_z_lag1 = lag(tmp_z, n = 1, default = NA),
         tmp_z_lag2 = lag(tmp_z, n = 2, default = NA),
         tmp_z_lag3 = lag(tmp_z, n = 3, default = NA),
        spei_z_lag1 = lag(spei_z, n = 1, default = NA),
         spei_z_lag2 = lag(spei_z, n = 2, default = NA),
         spei_z_lag3 = lag(spei_z, n = 3, default = NA),
    spei_lag1 = lag(spei, n = 1, default = NA),
    spei_lag2 = lag(spei, n = 2, default = NA),
    spei_lag3 = lag(spei, n = 3, default = NA),
          ) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup(.) 

nrow(dat_spei_lags)

# remove too low values: peiting and Eschenbach_idOPf: identify loutlier from log_sum_ips
df_outliers <- dat_spei_lags %>% 
  group_by(year) %>% 
  mutate(is_outlier = ifelse(sum_ips_log > quantile(sum_ips_log, 0.75,na.rm=T) + 1.5 * IQR(sum_ips_log,na.rm=T) |
                               sum_ips_log < quantile(sum_ips_log, 0.25,na.rm=T) - 1.5 * IQR(sum_ips_log,na.rm=T), TRUE, FALSE)) # %>% 

# check which ones are outliers and how often?

pair_outliers <- df_outliers %>%
  dplyr::filter(is_outlier) %>%
  group_by(trapID, pairID) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  dplyr::filter(n > 3) %>% 
  pull(pairID)

# if outlier > 4 time, remove the whole pair; keep the correct traps
# otherwise i need to remove 35 pairs, which is too much


# export table to merge it with the spatial data: for variograms:
dat_dynamics <- dat_fin %>% 
  dplyr::select(c(year, trapID, pairID, spring_tmp, veg_tmp, veg_prcp,   
                  #spei3, 
                  spei, # spei for veg season
                  tmp_z, prcp_z, spei_z, # z score fr veg season, sd from teh mean reference conditions
                sum_ips,  peak_doy, peak_diff, tr_agg_doy, tr_peak_doy, agg_doy, peak_doy
                ))
fwrite(dat_dynamics, 'outTable/beetle_dynamic_indicators.csv')

# table for ips counts, peak diff
dat_fin_counts_m <- 
  dat_spei_lags %>%  
  dplyr::select( -c(agg_doy, tr_agg_doy)) %>% 
  na.omit()
 
nrow(dat_fin_counts_m)  # 948 rws, not sure why????

# create new table with only agg values
dat_fin_agg_m <- 
  dat_spei_lags %>% 
   ungroup() %>%
  na.omit()
nrow(dat_fin_agg_m)

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

# Print the structure of the scaled data
str(dat_fin_counts_m_scaled)
summary(dat_fin_counts_m_scaled)



# Scale and center the remaining numeric columns
dat_fin_agg_m_scaled <- dat_fin_agg_m %>%
  mutate(across(
    .cols = where(is.numeric) & !all_of(columns_to_skip_agg),
    .fns = ~ scale(.) %>% as.vector()
  ))

# Print the structure of the scaled data
str(dat_fin_agg_m_scaled)
summary(dat_fin_agg_m_scaled)


####GAM: Find the best predictor lags  (simple) --------------------------------------------------------------- 
# maybe use???
# 
dependent_vars_counts <-  c("sum_ips", "peak_diff")
dependent_vars_doy    <-  c("tr_agg_doy", "tr_peak_doy")


# Initialize a data frame to store AIC values and deviance explained
model_metrics_count <- data.frame(Predictor = character(), 
                            Dependent = character(), 
                            AIC = numeric(), DevianceExplained = numeric())

# List of dependent variables and predictors:
# keep only spei3 and spei 12 - for the short term vs long term effect
selected_predictors <- c(#"spring_tmp", "spring_tmp_lag1", "spring_tmp_lag2","spring_tmp_lag3",
  #"veg_tmp", "veg_tmp_lag1", "veg_tmp_lag2","veg_tmp_lag3",
  "tmp_z", "tmp_z_lag1", "tmp_z_lag2","tmp_z_lag3" ,
  "spei_z", "spei_z_lag1", "spei_z_lag2","spei_z_lag3") 
  #"prcp_z", "prcp_z_lag1", "prcp_z_lag2","prcp_z_lag3",
  #'tmp', 'tmp_1yr', 'tmp_2yr','tmp_3yr',  # averages of tmp
  #'tmp_z', 'tmp_z_1yr', 'tmp_z_2yr','tmp_z_3yr',  # averages of tmp
 # 'prcp', 'prcp_1yr', 'prcp_2yr','prcp_3yr',  # averages of tmp
 # 'spei', 'spei_1yr', 'spei_2yr','spei_3yr',  # averages of tmp
  #                       "spei", "spei_lag1", "spei_lag2","spei_lag3"
#  )



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
#print(model_metrics_count)

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
         predictor = str_replace_all(predictor, "spei_z", "speiz")) %>% 
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
  facet_wrap(lag~type) +
  theme_bw()



# identify the 'hot and dry conditions'
plot(dat_fin_counts_m_scaled$veg_tmp,dat_fin_counts_m_scaled$spei)
plot(dat_fin_counts_m_scaled$tmp_z,dat_fin_counts_m_scaled$spei_z)

pairs(sum_ips ~tmp_z_lag2 + spei_z_lag2,   data  = dat_fin_counts_m_scaled, panel=panel.smooth)

# identify howt and ry conditions visually
# Load your data
dat <- dat_fin_counts_m_scaled_no_outliers

# Identify the points in the lower right cluster
cluster_points <- dat[dat$tmp_z > 1.5 & dat$spei_z < 0, ]

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
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs", k = 4)) + 
  facet_wrap(.~year, scale = 'free')



#data_filtered <- data %>% 
 # na.omit() # This removes rows with any NAs across all columns

# check outliers in beetle counts
boxplot(dat_fin_counts_m$sum_ips)
boxplot(log(dat_fin_counts_m$sum_ips))



# TEST: GAMM beetle counst with temp autocorr or with previous years beetle counts -------------------------------

# account to temp autocorrelation explicitely, within gamm

# try with gamm, and different specification of the temp autocorrelation
m7 <- gamm(sum_ips ~ s(year, k =4) + s(tmp_z_lag1, k = 4) + #s(spei) +
            s(x, y, bs = "gp") + s(trapID, bs = "re"),
          data = dat_fin_counts_m_scaled, 
          family = nb,
          #rho = 0.15, 
          correlation = corAR1(form = ~ year | trapID))

# add spei
m8 <- gamm(sum_ips ~ s(year, k =4) + s(tmp_z_lag1, k = 4) + s(spei_z_lag2, k = 4) +
             s(x, y, bs = "gp") + s(trapID, bs = "re"),
           data = dat_fin_counts_m_scaled, 
           family = nb,
           correlation = corAR1(form = ~ year | trapID))


# do not add previous years: instead use groupping on pairID
# add previous years beetle counts
# try with uncsaled data: also unscaled z-score
# as this can be easier to interpret in relationship to climate change??
m10_unsc <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 8) +
              s(x, y, bs = "gp") + 
                s(pairID, bs = "re"),
            data = dat_fin_counts_m, 
            family = nb,
            correlation = corAR1(form = ~ year | trapID))


# try less correlated predictors:
m10_unsc2 <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 8) +
                   s(x, y, bs = "gp") + 
                   s(pairID, bs = "re"),
                 data = dat_fin_counts_m, 
                 family = nb,
                 correlation = corAR1(form = ~ year | trapID))


# try less correlated predictors:
m10_unsc2_int <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 8) +
                        te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                    s(x, y, bs = "gp") + 
                    s(pairID, bs = "re"),
                  data = dat_fin_counts_m, 
                  family = nb,
                  correlation = corAR1(form = ~ year | trapID))

AIC(m10_unsc2, m10_unsc) # the tmp_z1 and spei_z2 are the least correlated!

plot(m10_unsc2$gam, page = 1)

# use in interaction different values of tmp then in teh main terms
# does not converge every time, maybe just skipp the interaction
m10_unsc_int <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag2, k = 15) + s(spei_z_lag2, k = 15) +
                       te(tmp_z, spei_z_lag2, k = 30) +
                   s(x, y, bs = "gp", k =30) + 
                     s(pairID, bs = "re"),
                 data = dat_fin_counts_m, 
                 family = nb,
                 correlation = corAR1(form = ~ year | trapID))
cor(dat_fin_counts_m$tmp_z_lag2, dat_fin_counts_m$spei_z_lag2, method = 'pearson' )
cor(dat_fin_counts_m$tmp_z_lag1, dat_fin_counts_m$spei_z_lag2, method = 'pearson' )



# the best!!!
mm <- m10_unsc2_int$gam
# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(mm, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(mm, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(mm, terms = c("tmp_z_lag1", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE) #[-1, 0, 1]

p_df <- as.data.frame(p2)

# test simple plot:
ggplot(p_df, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  #geom_point(data = dat_fin_counts_m, aes(x = spei_z_lag2, y = sum_ips, color=factor(year))) +
 # ylim(0,75000) +
  #labs(x = "SPEI Lag 2", y = "Sum of Beetle Counts") +
  theme_classic2()


# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m, terms = "tmp_z_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = c("tmp_z", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

p_df <- as.data.frame(p3)


# test simple plot:
ggplot(p3, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  ylim(0,75000) +
  #labs(x = "SPEI Lag 2", y = "Sum of Beetle Counts") +
  theme_classic2()




m10_unsc_no_year <- gamm(sum_ips ~  #s(year, k =6) + 
                           s(tmp_z_lag2, k = 8) + s(spei_z_lag2, k = 8) +
                   s(x, y, bs = "gp") + s(pairID, bs = "re"),
                 data = dat_fin_counts_m, 
                 family = nb,
                 correlation = corAR1(form = ~ year | trapID))

# treat year as random effect
m10_re <- gamm(sum_ips ~ s(tmp_z_lag2, k = 8) + s(spei_z_lag2, k = 8) +
                 s(x, y, bs = "gp") + s(pairID, bs = "re"),
               data = dat_fin_counts_m, 
               family = nb,
               correlation = corAR1(form = ~ year | trapID),
               random = list(year = ~1))



AIC(m10_unsc,m10_unsc_no_year,m10_re,m10_unsc_int)
summary(m10_unsc_no_year$gam)
gam.check(m10_unsc_no_year$gam)
plot(m10_unsc$gam, page = 1, shade = T)

#dat_fin_counts_m$AR
m10 <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag2, k = 8) + s(spei_z_lag2, k = 8) +
             s(x, y, bs = "gp") + s(pairID, bs = "re"),
           data = dat_fin_counts_m_scaled, 
           family = nb,
           correlation = corAR1(form = ~ year | trapID))

# Calculate average counts for each trap pair
avg_data <- dat_spei_lags %>%
  group_by(pairID, year) %>%
  summarise(sum_ips     = mean(sum_ips, na.rm = TRUE),
            tr_agg_doy  = mean(tr_agg_doy, na.rm = TRUE),
            tr_peak_doy = mean(tr_peak_doy, na.rm = TRUE), 
            tr_agg_doy  = mean(tr_agg_doy, na.rm = TRUE),
            agg_doy     = mean(agg_doy, na.rm = TRUE),
            peak_doy    = mean(peak_doy, na.rm = TRUE),
            spei_z_lag2 = mean(spei_z_lag2, na.rm = TRUE),
            tmp_z_lag1  = mean(tmp_z_lag1, na.rm = TRUE),
            x           = mean(x, na.rm = TRUE),
            y           = mean(y, na.rm = TRUE),
            ) %>%
  ungroup(.) %>% 
  na.omit()



# test on averaged design
# try less correlated predictors:
m1 <- gamm(sum_ips ~  s(year, k =6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) +
                        te(tmp_z_lag1, spei_z_lag2, k = 15) + 
                        s(x, y, bs = "gp") + 
                        s(pairID, bs = "re"),
                      data = avg_data, 
                      family = nb,
                      correlation = corAR1(form = ~ year | pairID))
# remove xy??
m2 <- gamm(sum_ips ~  s(year, k =6) +  s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) +
             te(tmp_z_lag1, spei_z_lag2, k = 20) + 
            # s(x, y, bs = "gp") + 
             s(pairID, bs = "re"),
           data = avg_data, 
           family = nb,
           correlation = corAR1(form = ~ year | pairID))

# add spatial autocorrelation explicitely, try simple
# Define the model with spatial autocorrelation
m3 <- gamm(sum_ips ~ s(year, k = 6) + s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 4) +
             te(tmp_z_lag1, spei_z_lag2, k = 15) + s(pairID, bs = "re"),
           data = avg_data,
           family = nb,
           correlation = corAR1(form = ~ year | pairID) +
             corSpatial(form = ~ x + y, type = "exponential"))


m<- m2$gam
fin.m.count <- m2

AIC(m1$lme, m2$lme)
AIC(m1$gam, m2$gam)

summary(m)
plot.gam(m)
plot.lme( m1$lme )

acf(residuals(m2$gam))
pacf(residuals(m1$gam))


# chec for residuals:

library(gstat)

# Extract residuals from m2
residuals_m2 <- residuals(m2$lme, type = "normalized")

# Create a spatial data frame
spatial_data <- data.frame(x = avg_data$x, y = avg_data$y, residuals = residuals_m2)

# Plot residuals to check for spatial patterns
plot(spatial_data$x, spatial_data$y, col = residuals_m2, pch = 16, 
     xlab = "X Coordinate", ylab = "Y Coordinate", 
     main = "Spatial Plot of Residuals (Model m2)")
colorRampPalette(c("blue", "white", "red"))(100)

# Variogram of residuals
vgm_res <- variogram(residuals ~ 1, data = spatial_data, locations = ~ x + y)
plot(vgm_res, main = "Variogram of Residuals (Model m2)")
# residual seems not to be spatialy autocorrelated!

# Automate approach??? SUM_IPS ----------------------



# AGG DOY : test for betar family-----------------------------------
# test if it is better to keep pairs as random effect or just the xy coordinates? 
m.re <- gamm(tr_agg_doy ~ s(year, k = 6) + s(tmp_z_lag1, k = 3) + s(spei_z_lag2, k = 3) +
             te(tmp_z_lag1, spei_z_lag2, k = 15) + 
             #s(x, y, bs = "gp") + 
             s(pairID, bs = "re"),
           data = avg_data,
           family = betar,
           correlation = corAR1(form = ~ year | pairID))

m.xy <- gamm(tr_agg_doy ~ s(year, k = 6) +  s(tmp_z_lag1, k = 3) + s(spei_z_lag2, k = 3) +
             te(tmp_z_lag1, spei_z_lag2, k = 15) + 
             s(x, y, bs = "gp") ,
           data = avg_data,
           family = betar,
           correlation = corAR1(form = ~ year | pairID))

# both
m.both.spat <- gamm(tr_agg_doy ~ s(year, k = 6) +  s(tmp_z_lag1, k = 8) + s(spei_z_lag2, k = 3) +
               te(tmp_z_lag1, spei_z_lag2, k = 15) + 
               s(x, y, bs = "gp") +
                 s(pairID, bs = "re"),
             data = avg_data,
             family = betar,
             correlation = corAR1(form = ~ year | pairID))



m1 <- m.re$gam
m2 <- m.xy$gam
m3 <- m.both.spat$gam

r2(m3)


avg_data %>% 
  ggplot(aes(x = spei_z_lag2,
             y = tr_agg_doy)) +
  geom_smooth()



# quick plotting
fin.m <- m3
#  m.spei12_no_re_slope47# fin.m.peak.diff
#!!! --------------------------------------------

# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m, terms = "year [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m, terms = c("tmp_z_lag1", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

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

ggarrange(plot1,plot2,plot3,plot4, ncol = 4, nrow = 1)


AIC(m.re$lme, m.xy$lme )
# compare the two: add xy or add pairs as random effects? what is better?

m<-m1$gam
summary(m)
plot(m, page = 1)

# automate for DOY as well:  -----------------------------------------
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
#dependent_vars  <- c("tr_agg_doy")
# Run the function on your dataset
result <- compare_models(avg_data, dependent_vars)

m <- result$tr_agg_doy$models$random_effect$gam
summary(m)
plot.gam(m)

# Plot interaction term
draw(gam_model, select = "te(tmp_z_lag1,spei_z_lag2)")
draw(gam_model, select = "s(tmp_z_lag1)")
library(gratia)
gratia::draw(m)
appraise(m)
draw(m1, residuals = TRUE)



# explore model with linea terms
m2 <- gamm(tr_agg_doy ~ s(year, k = 6) + s(tmp_z_lag1,3), # +# s(spei_z_lag2,3) + te(tmp_z_lag1, spei_z_lag2, k = 8) + s(x, y, bs = 'gp') + s(pairID, bs = 're'),
           family = betar(link = "logit"),data = avg_data,  corAR1(form = ~ year | pairID))




# how many observation s i have per year?

yearly_counts <- dat_fin_counts_m_scaled %>%
  group_by(year, trapID) %>%
  summarise(count = n())

yearly_counts %>% 
  dplyr::filter(count > 1)
### glmm: SUM_IPS vegt_tmp + spei3_lag2  -----------------------------------------

m1 <- glmmTMB(sum_ips ~ veg_tmp + spei3_lag2,
              family = nbinom2,
              data = dat_fin_counts_m_scaled )
m1.nb <- glm.nb(sum_ips ~ veg_tmp + spei3_lag2,
              family = nbinom2,
              data = dat_fin_counts_m_scaled)

# veg_tmp s slightly better then veg_tmp_lag2, proceed with veg_tmp
m2 <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2,
              family = nbinom2,
              data = dat_fin_counts_m_scaled)
# add random effects
m3 <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m_scaled)

m3.time <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2 + 
                     (1|year) ,# random intercept for year
              family = nbinom2,
              data = dat_fin_counts_m_scaled)

# by Marc!
m4.time <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2 + 
                      (1 + year|pairID) ,
                    # (year|pairID) ,# random intercept for year
                   family = nbinom2,
                   data = dat_fin_counts_m_scaled)

m4_z <- glmmTMB(sum_ips ~ tmp_z*spei_z + 
                     (year|pairID) ,# random intercept for year
                   family = nbinom2,
                   data = dat_fin_counts_m_scaled)



AIC(m3, m3.time, m4.time,m4_z)
check_collinearity(m4_z)
summary(m4_z)

# exclue random effect
m4.poly <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + spei3_lag2 + veg_tmp:spei3_lag2,# + (1|pairID),
             # family = nbinom2,
              data = dat_fin_counts_m)
check_collinearity(m4.poly)


#AIC(m1,m1.nb, m2,m3,m3.time, m4.poly, m4)

# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + spei3_lag2 + veg_tmp:spei3_lag2 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)
AIC(m1,m1.nb, m2,m3,m3.time, m4.poly, m4)

check_collinearity(m4)
m5 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + poly(spei3_lag2,2) + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)




# SUM_IPS: gam: spei_z_lag1  , tmp_z ------------------------------------------- 










dat_fin_counts_m_scaled_no_outliers <- dat_fin_counts_m_scaled %>% 
  dplyr::filter(!pairID %in% pair_outliers) %>% 
  mutate(year_fact = as.factor(year))


# Fitting the model
m_no_year <- gam(sum_ips ~ s(tmp_z, k = 4) + 
           s(spei_z_lag1, k = 4) + 
           s(trapID, bs = 're'), #+ 
           #s( year_fact, bs = 're'), 
         family = nb(), 
         method = 'REML', 
         data = dat_fin_counts_m_scaled_no_outliers)


m <- gam(sum_ips ~ s(tmp_z, k = 4) + # also, increasing k from 4 to 15 increase multicoll
           s(spei_z_lag1, k = 4) + 
           s(trapID, bs = 're') + 
           s( year_fact, bs = 're'), 
         family = nb(), 
         method = 'REML', 
         data = dat_fin_counts_m_scaled_no_outliers)

# including interaction included high multicollinearity
m1 <- gam(sum_ips ~ s(tmp_z, k = 4) + 
           s(spei_z_lag1, k = 4) + 
            ti(tmp_z, spei_z_lag1, k = 4) +
           s(trapID, bs = 're') + 
           s( year_fact, bs = 're'), 
         family = nb(), 
         method = 'REML', 
         data = dat_fin_counts_m_scaled_no_outliers)


m2 <- gam(sum_ips ~ #s(tmp_z, k = 4) + 
            s(spei_z_lag1, k = 4) + 
            ti(tmp_z, spei_z_lag1, k = 4) +
            s(trapID, bs = 're') + 
            s( year_fact, bs = 're'), 
          family = nb(), 
          method = 'REML', 
          data = dat_fin_counts_m_scaled_no_outliers)


m3 <- gam(sum_ips ~ #s(tmp_z, k = 4) + 
            #s(spei_z_lag1, k = 4) + 
            ti(tmp_z, spei_z_lag1, k = 4) +
            s(trapID, bs = 're') + 
            s( year_fact, bs = 're'), 
          family = nb(), 
          method = 'REML', 
          data = dat_fin_counts_m_scaled_no_outliers)
AIC(m, m1, m2, m3)
plot(m3)
concurvity(m)
concurvity(m2)

# check how my data looks like: tmp vs drought: 
plot(dat_fin_counts_m_scaled_no_outliers$tmp_z, 
    dat_fin_counts_m_scaled_no_outliers$spei_z)



# try: random intercept for locations, random slope for year 
m1 <- gam(sum_ips ~ s(spei_z_lag2, k = 4) +
            s(trapID, bs = 're') +
            s(year, spei_z_lag2, bs = 're'),
            #s(tmp, )
            #, 
                   family = nb(),
                   method ='REML',  
                   data = dat_fin_counts_m_scaled_no_outliers)


m2 <- gam(sum_ips ~ s(spei_z_lag1, k = 4) +
            s(tmp_z_lag2, k = 4) +
            s(trapID, bs = 're') +
            s(tmp_z, by = year, bs = 'fs'),
          #s(tmp, )
          #, 
          family = nb(),
          method ='REML',  
          data = dat_fin_counts_m_scaled_no_outliers)

AIC(m1, m2)
summary(m2)
plot(m2, page = 1)

m.spei_lag1 <- gam(sum_ips ~ s(spei_lag1, k = 4) , 
                   family = nb(),
                   method ='REML',  
                   data = dat_fin_counts_m_scaled_no_outliers)

m.spei_lag2 <- gam(sum_ips ~ s(spei_lag1, k = 4), 
                   family = nb(),
                   method ='REML',
                   data = dat_fin_counts_m_scaled_no_outliers)

m.tmp <- gam(sum_ips ~  s(tmp_z, k = 4), 
             family = nb(),
             method ='REML',
             data = dat_fin_counts_m_scaled_no_outliers)

m.tmp_lag1 <- gam(sum_ips ~  s(tmp_z, k = 4), 
             family = nb(),
             method ='REML',  
             data = dat_fin_counts_m_scaled_no_outliers)


m.add <- gam(sum_ips ~ s(tmp_z, k = 4) +  s(spei_z_lag1, k = 4),# +
             #  ti(tmp_z_lag1, spei_z_lag2, k =7), 
             family = nb(),
             method ='REML', 
             data = dat_fin_counts_m_scaled_no_outliers)



m.int <- gam(sum_ips ~ s(tmp_z, k = 4) +  s(spei_z_lag1, k = 4) +
               ti(tmp_z, spei_z_lag1, k =7) +
               pairID, 
             family = nb(),
             method ='REML', 
             data = dat_fin_counts_m_scaled_no_outliers)


m.int2 <- gam(sum_ips ~ s(tmp_z_lag1, k = 4) +  s(spei_z_lag2, k = 4) +
               ti(tmp_z_lag1, spei_z_lag2, k =7) +
               pairID, 
             family = nb(),
             method ='REML', 
             data = dat_fin_counts_m_scaled_no_outliers)


AIC(m.int, m.int2)

summary(m.int2)
check_collinearity(m.int2)
concurvity(m.int2)
plot(m.int2, page = 1)








# identify significant traps:
# Extract significant traps from the model summary
summary_m.int <- summary(m.int)
# Extract coefficients
coefficients <- summary_m.int$p.coeff

# Extract p-values
p_values <- summary_m.int$p.pv


significant_traps <- 
  data.frame(coefficients = coefficients,
                                p_values  = p_values) %>% 
  
   rownames_to_column("pairID") %>%
    as.data.frame() %>% 
   # head()
   dplyr::filter(p_values < 0.05) %>% 
   mutate(pairID = gsub("pairID", '', pairID )) 


(significant_traps)


# split data into significant and non significant locations
df_sign <- dat_fin_counts_m_scaled_no_outliers %>%
  mutate(significance = ifelse(pairID %in% significant_traps$pairID, "sig", "ref")) %>% 
  dplyr::filter(significance == 'sig')


df_no_sign <- dat_fin_counts_m_scaled_no_outliers %>%
  mutate(significance = ifelse(pairID %in% significant_traps$pairID, "sig", "ref")) %>% 
  dplyr::filter(significance == 'ref')


# use different location as its own factor:
dat_fin_counts_m_scaled_no_outliers <- dat_fin_counts_m_scaled_no_outliers %>%
  mutate(significance = as.factor(ifelse(pairID %in% significant_traps$pairID, "sig", "ref"))) 

m <- gam(sum_ips ~ s(tmp_z, k = 4) + s(spei_z, k = 4) + 
           significance + 
           trapID,  # this induces collinearity/concurvity
         family = nb(),
         method ='REML', 
         dat_fin_counts_m_scaled_no_outliers)

summary(m)
check_concurvity(m)
plot(m, page = 1)

# try model just for different location stypes: signigicantly ifferent or similar to the ref
m.sign <- gam(sum_ips ~ s(tmp_z, k = 4) + s(spei_z, k = 4),# + 
          # significance + 
          # trapID,  # this induces collinearity/concurvity
         family = nb(),
         method ='REML', 
         df_sign)

summary(m.sign)
m.no_sign <- gam(sum_ips ~ s(tmp_z, k = 4) + s(spei_z, k = 4),# + 
              # significance + 
              # trapID,  # this induces collinearity/concurvity
              family = nb(),
              method ='REML', 
              df_no_sign)

plot(m.no_sign, page = 1)
summary(m.no_sign)

# no meaningful to check for them separately

# get list of outliers and remove the traps --------------------------------------

dat_fin_counts_m %>% 
  ggplot(aes(x = year,
             group = year,
             y = sum_ips_log)) +
  geom_boxplot() +
  geom_text(data = df_outliers %>% dplyr::filter(is_outlier),
            aes(label = trapID),
            position = position_jitter(width = 0.2, height = 0),
            hjust = -0.1, vjust = -0.5, size = 3, color = "red")


boxplot(dat_fin_counts_m$sum_ips_log)

# Combine data for visualization
climate_data <- dat_fin_counts_m %>%
  mutate(significance = ifelse(pairID %in% significant_traps$pairID, "sig", "ref"))

# Plot temperature comparison
p_temp <- ggplot(climate_data, aes(x = significance, y = tmp_z, fill = significance)) +
  geom_boxplot() +
  stat_compare_means(aes(group = significance), label = "p.signif") +
  labs(title = "Comparison of Temperature (tmp_z) between Significant and Non-Significant Traps",
       x = "Trap Significance", y = "Temperature (tmp_z)") +
  theme_minimal() + 
  facet_grid(.~year, scales = 'free')

# Plot spei comparison
p_spei <- ggplot(climate_data, aes(x = significance, y = spei_z, fill = significance)) +
  geom_boxplot() +
  stat_compare_means(aes(group = significance), label = "p.signif") +
  labs(title = "Comparison of SPEI (spei_z) between Significant and Non-Significant Traps",
       x = "Trap Significance", y = "SPEI (spei_z)") +
  theme_minimal() + 
  facet_grid(.~year, scales = 'free')


p_ips_counts <- ggplot(climate_data, aes(x = significance, y = sum_ips, 
                                         fill = significance)) +
  geom_boxplot() +
  stat_compare_means(aes(group = significance), label = "p.signif") +
  labs(title = "Comparison of beet counts between Significant and Non-Significant Traps",
       x = "Trap Significance", y = "Beetle pop level") +
  theme_minimal() + 
  facet_grid(.~year, scales = 'free')

p_ips_counts_log <- ggplot(climate_data, aes(x = significance, y = log(sum_ips), 
                                         fill = significance)) +
  geom_boxplot() +
  labs(title = "Comparison of beetle counts (log) between Significant and Non-Significant Traps",
       x = "Trap Significance", y = "Beetle pop level") +
  theme_minimal() + 
  facet_grid(.~year, scales = 'free')


windows()
ggarrange(p_temp, p_spei,p_ips_counts, nrow = 3) # p_ips_counts_log,



# Perform t-tests to compare climate variables
t_test_tmp <- t.test(tmp_z_lag1 ~ significance, data = climate_data)
t_test_spei <- t.test(spei_z_lag2 ~ significance, data = climate_data)

print(t_test_tmp)
print(t_test_spei)


# run t test for years


# Assume climate_data is your data frame and it has columns year, tmp_z_lag1, spei_z_lag2, and significance

# Initialize lists to store results
t_test_tmp_results <- list()
t_test_spei_results <- list()

# Get unique years
years <- unique(climate_data$year)

# Loop through each year and perform t-tests
for (year in years) {
  # Filter data for the current year
  data_year <- climate_data %>% filter(year == !!year)
  
  # Perform t-test for tmp_z_lag1
  t_test_tmp <- t.test(tmp_z_lag1 ~ significance, data = data_year)
  t_test_tmp_results[[as.character(year)]] <- round(t_test_tmp$p.value, 2)
  
  # Perform t-test for spei_z_lag2
  t_test_spei <- t.test(spei_z_lag2 ~ significance, data = data_year)
  t_test_spei_results[[as.character(year)]] <- round(t_test_spei$p.value,2)
}

# Convert results to data frames
t_test_tmp_results_df <- data.frame(year = names(t_test_tmp_results), p_value_tmp = unlist(t_test_tmp_results))
t_test_spei_results_df <- data.frame(year = names(t_test_spei_results), p_value_spei = unlist(t_test_spei_results))

# Print the results
print(t_test_tmp_results_df)
print(t_test_spei_results_df)



# run for current tmps
# Initialize lists to store results
t_test_tmp_results <- list()
t_test_spei_results <- list()

# Get unique years
years <- unique(climate_data$year)

# Loop through each year and perform t-tests
for (year in years) {
  # Filter data for the current year
  data_year <- climate_data %>% filter(year == !!year)
  
  # Perform t-test for tmp_z_lag1
  t_test_tmp <- t.test(tmp_z ~ significance, data = data_year)
  t_test_tmp_results[[as.character(year)]] <- round(t_test_tmp$p.value,2)
  
  # Perform t-test for spei_z_lag2
  t_test_spei <- t.test(spei_z ~ significance, data = data_year)
  t_test_spei_results[[as.character(year)]] <- round(t_test_spei$p.value,2)
}

# Convert results to data frames
t_test_tmp_results_df <- data.frame(year = names(t_test_tmp_results), p_value_tmp = unlist(t_test_tmp_results))
t_test_spei_results_df <- data.frame(year = names(t_test_spei_results), p_value_spei = unlist(t_test_spei_results))

# Print the results
print(t_test_tmp_results_df)
print(t_test_spei_results_df)



# PCA merge spei and tmp into one variable? ---------------------------------------------
# see the clusters;
dat_fin_counts_m_scaled_no_outliers %>% 
  ggplot(aes(x = tmp_z,
             y = spei_z,
             color = year_fact)) +
  geom_point() + 
  #geom_smooth(
  #            aes(data =  dat_fin_counts_m_scaled_no_outliers, 
  #                x = tmp_z,
   #               y = spei_z)) +
  theme_classic() + 
  theme(aspect.ratio = 1)


  
m.lm <- lm(spei_z ~ tmp_z, dat_fin_counts_m_scaled_no_outliers)
plot(m.lm, page = 1)

pc <- prcomp(dat_fin_counts_m_scaled_no_outliers[,c('tmp_z', 'spei_z')],
             center = TRUE,
             scale. = TRUE)
attributes(pc)
summary(pc)


# ad axis into data: keep only PC1
dat_fin_counts_m_scaled_no_outliers$PC1 <-  pc$x[,1]


dat_fin_counts_m_scaled_no_outliers %>% 
  ggplot(aes(x = 1:nrow(dat_fin_counts_m_scaled_no_outliers), y = PC1, color = factor(year_fact))) +
  geom_point() +
  labs(x = "Index", y = "PC1", color = "Year") +
  ggtitle("Combined Index of tmp_z and spei_z")



# check what my values mean:
dat_fin_counts_m_scaled_no_outliers %>% 
  dplyr::select(PC1, tmp_z, spei_z) %>% 
  View()

# interpret the axis:
pairs(sum_ips ~ PC1+spei_z+tmp_z,data = dat_fin_counts_m_scaled_no_outliers)

hist(dat_fin_counts_m_scaled_no_outliers$PC1)

dat_fin_counts_m_scaled_no_outliers %>% 
  ggplot(aes(y = sum_ips,
         x = PC1,
         color = year_fact)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(.~year, scales = 'free')



# allow different slope per year
m <- gam(sum_ips ~ s(PC1, k = 20) +  s(trapID, bs = 're'), 
         family = nb(), 
         dat_fin_counts_m_scaled_no_outliers)

summary(m)
plot(m, page = 1)
gam.check(m)
check_concurvity(m)


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







# SUM_IPS: gam: spei3_lag2  , veg_tmp ------------------------------------------- 
# what is more important: lagged affect of droung and current tmp!
# and what spei lag to use? 1 or 2? test!
m.spei_lag1 <- gam(sum_ips ~ s(spei3_lag1, k = 4) , 
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

m.spei_lag2 <- gam(sum_ips ~ s(spei3_lag2, k = 4), 
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

m.tmp <- gam(sum_ips ~  s(veg_tmp, k = 4), 
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

m.int <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 7) +
               ti(veg_tmp, spei3_lag2, k =7), 
             family = nb(),
             method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
             data = dat_fin_counts_m_scaled)


AIC(m.tmp, m.spei_lag2,m.spei_lag1,m.int)

# introducing random effects improves model a lot, but inceases VIF to 11 and 17 ! (spei = 11, veg_tmp = 17)
m.int.pair <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 7) +
               ti(veg_tmp, spei3_lag2, k = 7) +
                 s(pairID, bs= 're'), 
             family = nb(),
             method ='REML',  
             data = dat_fin_counts_m_scaled)


m.int.re.trap <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 7) +
                  ti(veg_tmp, spei3_lag2, k = 7) +
                  s(trapID, bs= 're'), 
                family = nb(),
                method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                data = dat_fin_counts_m_scaled)

# increase k
m.int.re.trap.pair <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 5) +
                  ti(veg_tmp, spei3_lag2, k = 5) +
                  s(pairID, trapID, bs= 're'), #+
                 #   s(trapID, bs= 're'), 
                family = nb(),
                method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                data = dat_fin_counts_m_scaled)


# Fit the model using gam with regularization
m.int.re.trap.pair.nest.reg <- gam(sum_ips ~ s(spei3_lag2, k = 6) +
                                 s(veg_tmp, k = 10) +
                                 ti(veg_tmp, spei3_lag2, k = 12) +
                                 s(pairID, bs= 're') +
                                 s(trapID, bs= 're'),  # nested design
                               family = nb(),
                               method = 'REML',
                               data = dat_fin_counts_m_scaled)





m.int.slope.spei <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                  ti(veg_tmp, spei3_lag2, k = 4) +
                  #s(pairID, bs= 're') +
                  s(spei3_lag2, pairID, bs = "re"), 
                family = nb(),
                method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                data = dat_fin_counts_m_scaled)

m.int.slope.tmp <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                             ti(veg_tmp, spei3_lag2, k = 4) +
                             #s(pairID, bs= 're') +
                             s(veg_tmp, pairID, bs = "re"), 
                           family = nb(),
                           method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                           data = dat_fin_counts_m_scaled)

# keep random slope only for temperature: seems as a better predictor
m.int.both.slopes <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                            ti(veg_tmp, spei3_lag2, k = 4) +
                            #s(pairID, bs= 're') +
                            s(veg_tmp, pairID, bs = "re") +
                              s(spei3_lag2, pairID, bs = "re"), 
                          family = nb(),
                          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                          data = dat_fin_counts_m_scaled)

# add random effect back:
m.int.re.both.slopes <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                              ti(veg_tmp, spei3_lag2, k = 4) +
                              s(pairID, bs= 're') +
                              s(veg_tmp, pairID, bs = "re") +
                              s(spei3_lag2, pairID, bs = "re"), 
                            family = nb(),
                            method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                            data = dat_fin_counts_m_scaled)

# remove tmp as has high VIF (30)
m.int.re.both.slopes.no.tmp <- gam(sum_ips ~ s(spei3_lag2, k = 4) + # s(veg_tmp, k = 4) +
                              ti(veg_tmp, spei3_lag2, k = 4) +
                              s(pairID, bs= 're') +
                              s(veg_tmp, pairID, bs = "re") +
                              s(spei3_lag2, pairID, bs = "re"), 
                            family = nb(),
                            method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                            data = dat_fin_counts_m_scaled)


# remove int term, as now both terms are correlate
m.re.both.slopes <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                                    # ti(veg_tmp, spei3_lag2, k = 4) +
                                     s(pairID, bs= 're') +
                                     s(veg_tmp, pairID, bs = "re") +
                                     s(spei3_lag2, pairID, bs = "re"), 
                                   family = nb(),
                                   method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                   data = dat_fin_counts_m_scaled)



# remove int term, as now both terms are correlate. remove random slope for veg_tmp
m.re.spei.slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                          # ti(veg_tmp, spei3_lag2, k = 4) +
                          s(pairID, bs= 're') +
                         # s(veg_tmp, pairID, bs = "re") +
                          s(spei3_lag2, pairID, bs = "re"), 
                        family = nb(),
                        method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                        data = dat_fin_counts_m_scaled)


m.re <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                         # ti(veg_tmp, spei3_lag2, k = 4) +
                         s(pairID, bs= 're') ,#+
                         # s(veg_tmp, pairID, bs = "re") +
                         #s(spei3_lag2, pairID, bs = "re"), 
                       family = nb(),
                       method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                       data = dat_fin_counts_m_scaled)

# re on trap level
m.re.trap <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
              # ti(veg_tmp, spei3_lag2, k = 4) +
              s(trapID, bs= 're') ,#+
            # s(veg_tmp, pairID, bs = "re") +
            #s(spei3_lag2, pairID, bs = "re"), 
            family = nb(),
            method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
            data = dat_fin_counts_m_scaled)

# re on trap level
m.re.trap <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 4) +
                    ti(veg_tmp, spei3_lag2, k = 4) +
                   s(trapID, bs= 're') ,#+
                 # s(veg_tmp, pairID, bs = "re") +
                 #s(spei3_lag2, pairID, bs = "re"), 
                 family = nb(),
                 method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                 data = dat_fin_counts_m_scaled)





# test gamm
# Fit a GAMM
m.gamm <- gamm(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 4),
               random = list(pairID = ~1),
               family = nb())



m.re <- gam(sum_ips ~ s(spei3_lag2, k = 4, bs = "tp", select = TRUE) +
              s(veg_tmp, k = 4, bs = "tp", select = TRUE) +
              s(pairID, bs = "re"),
            family = nb(),
            method = "REML")

AIC(m.re,m.re.trap)

AIC(m.int.re, m.int.slope.spei, m.int.slope.tmp, m.int.both.slopes, m.int.re.both.slopes)

m.re.slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp, k = 4) +
                       # ti(veg_tmp, spei3_lag2, k = 4) +
                        s(pairID, bs= 're') +
                        s(veg_tmp, pairID, bs = "re"), 
                      family = nb(),
                      method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                      data = dat_fin_counts_m_scaled)

m.re.slopetmp <- gam(sum_ips ~ #s(spei3_lag2, k = 4) +  
                    s(veg_tmp, k = 4) +
                    #ti(veg_tmp, spei3_lag2, k = 4) +
                    s(pairID, bs= 're'),# +
                    #s(veg_tmp, pairID, bs = "re"), 
                  family = nb(),
                  method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                  data = dat_fin_counts_m_scaled)

AIC(m.re.slope,m.int.re.slope,m.re.slopetmp)


AIC(m.tmp, m.spei_lag2,m.spei_lag1,m.int,m.int.re,m.int.re.slope) # best is tmp and spei_lag2


check_collinearity(m.re.slope)


# how to interpret random effects?? ---------------
# Extract random effects for 'pairID'
random_effects_intercept <- coef(m.int.re.slope)$random["pairID"]
random_effects_slope <- coef(m.int.re.slope)$random["spei3_lag2:pairID"]

# Convert to data frame for plotting
random_effects_df <- data.frame(
  pairID = levels(dat_fin_counts_m_scaled$pairID),
  intercept = random_effects_intercept,
  slope = random_effects_slope
)

# Plot random intercepts
ggplot(random_effects_df, aes(x = pairID, y = intercept)) +
  geom_bar(stat = "identity") +
  labs(title = "Random Intercepts for Trap Pairs", x = "Trap Pair ID", y = "Random Intercept") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot random slopes
ggplot(random_effects_df, aes(x = pairID, y = slope)) +
  geom_bar(stat = "identity") +
  labs(title = "Random Slopes for spei3_lag2 by Trap Pairs", x = "Trap Pair ID", y = "Random Slope") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# SUM_IPS: gam: spei3_lag2  , veg_tmp_lag2      ----------------------------------
# using both values lagged invert my relationship! I have more beetles at higher wetness!

m1 <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4), 
                  family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                  data = dat_fin_counts_m_scaled)
check_collinearity(m2.re)
#cor(dat_fin_counts_m_scaled$spei3_lag2,dat_fin_counts_m_scaled$veg_tmp_lag2, method = 'spearman')

m2.re <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
               s(trapID, bs= 're'),
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)
# adding re has high multicollinearity!

m3 <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
               ti(spei3_lag2, veg_tmp_lag2, k = 4),
             family = nb(),
             method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
             data = dat_fin_counts_m_scaled)

# add random slope
m3.sl.spei <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
            ti(spei3_lag2, veg_tmp_lag2, k = 4) +
              s(spei3_lag2, pairID, bs = "re"), # ,
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

m3.sl.tmp <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                    ti(spei3_lag2, veg_tmp_lag2, k = 4) +
                   s(veg_tmp_lag2, pairID, bs = "re"), #,
                  family = nb(),
                  method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                  data = dat_fin_counts_m_scaled)

# both random slopes:
m3.sl <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                   ti(spei3_lag2, veg_tmp_lag2, k = 4) +
               s(spei3_lag2, pairID, bs = "re") +
                   s(veg_tmp_lag2, pairID, bs = "re"), #,
                 family = nb(),
                 method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                 data = dat_fin_counts_m_scaled)

# add back random effect
# both random slopes:
m3.sl.re <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
               ti(spei3_lag2, veg_tmp_lag2, k = 4) +
                 s(pairID, bs = "re") +
               s(spei3_lag2, pairID, bs = "re") +
               s(veg_tmp_lag2, pairID, bs = "re"), #,
             family = nb(),
             method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
             data = dat_fin_counts_m_scaled)


m3.int.re <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                  ti(spei3_lag2, veg_tmp_lag2, k = 4) +
                  s(pairID, bs = "re"), #+
                  #s(spei3_lag2, pairID, bs = "re") +
                  #s(veg_tmp_lag2, pairID, bs = "re"), #,
                family = nb(),
                method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                data = dat_fin_counts_m_scaled)

m3.int.re445 <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                   ti(spei3_lag2, veg_tmp_lag2, k = 5) +
                   s(pairID, bs = "re"), #+
                 #s(spei3_lag2, pairID, bs = "re") +
                 #s(veg_tmp_lag2, pairID, bs = "re"), #,
                 family = nb(),
                 method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                 data = dat_fin_counts_m_scaled)

m3.int.re447 <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                      ti(spei3_lag2, veg_tmp_lag2, k = 7) +
                      s(pairID, bs = "re"), #+
                    #s(spei3_lag2, pairID, bs = "re") +
                    #s(veg_tmp_lag2, pairID, bs = "re"), #,
                    family = nb(),
                    method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                    data = dat_fin_counts_m_scaled)


m3.int.re4410 <- gam(sum_ips ~ s(spei3_lag2, k = 4) +  s(veg_tmp_lag2, k = 4) +
                      ti(spei3_lag2, veg_tmp_lag2, k = 10) +
                      s(pairID, bs = "re"), #+
                    #s(spei3_lag2, pairID, bs = "re") +
                    #s(veg_tmp_lag2, pairID, bs = "re"), #,
                    family = nb(),
                    method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                    data = dat_fin_counts_m_scaled)



AIC(m3.sl.tmp, m3.sl.spei, m1, m2.re, m3,m3.sl,m3.sl.re, m3.int.re445,m3.int.re447,m3.int.re4410)
summary(m3.int.re4410)
plot(m3.int.re445, page = 1)
check_collinearity(m3.int.re4410)

# add random slopes:
m.gam.spei_re_slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                             #s(pairID, bs= 're') +
                             s(spei3_lag2, pairID, bs = "re"), #+ ,
                           family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                           data = dat_fin_counts_m_scaled)

# add temp - if spei is by use as random effect, multicoll is low
m.gam.spei_tmp_re_slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 4) +
                                 # s(pairID, bs= 're') +
                                 s(spei3_lag2, pairID, bs = "re"), #+ ,
                               family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                               data = dat_fin_counts_m_scaled)

# try different k if it will help to explain deviance?
m.gam.spei_tmp_re_slope77 <- gam(sum_ips ~ s(spei3_lag2, k = 7) + s(veg_tmp, k = 7) +
                                   # s(pairID, bs= 're') +
                                   s(spei3_lag2, pairID, bs = "re"), #+ ,
                                 family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                 data = dat_fin_counts_m_scaled)






# test gamm: - replace spei3_lag2 by spei12_lag2 as it is the least correl;ated and the explains teh most deviance
m.gam.spei <- gam(sum_ips ~ s(spei3_lag2, k = 4), 
              family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
              data = dat_fin_counts_m_scaled)

m.gam.spei_re <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                       s(pairID, bs= 're'),
                  family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                  data = dat_fin_counts_m_scaled)

# add random slopes:
m.gam.spei_re_slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                       #s(pairID, bs= 're') +
                         s(spei3_lag2, pairID, bs = "re"), #+ ,
                     family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                     data = dat_fin_counts_m_scaled)

# add temp - if spei is by use as random effect, multicoll is low
m.gam.spei_tmp_re_slope <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 4) +
                            # s(pairID, bs= 're') +
                             s(spei3_lag2, pairID, bs = "re"), #+ ,
                           family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                           data = dat_fin_counts_m_scaled)

# try different k if it will help to explain deviance?
m.gam.spei_tmp_re_slope77 <- gam(sum_ips ~ s(spei3_lag2, k = 7) + s(veg_tmp, k = 7) +
                                 # s(pairID, bs= 're') +
                                 s(spei3_lag2, pairID, bs = "re"), #+ ,
                               family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                               data = dat_fin_counts_m_scaled)

# lower k for tmp as moederate correlation
m.gam.spei_tmp_re_slope75 <- gam(sum_ips ~ s(spei3_lag2, k = 7) + s(veg_tmp, k = 5) +
                                   # s(pairID, bs= 're') +
                                   s(spei3_lag2, pairID, bs = "re"), #+ ,
                                 family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                 data = dat_fin_counts_m_scaled)

# lower k for tmp and spei as moederate correlation
m.gam.spei_tmp_re_slope55 <- gam(sum_ips ~ s(spei3_lag2, k = 5) + s(veg_tmp, k = 5) +
                                   # s(pairID, bs= 're') +
                                   s(spei3_lag2, pairID, bs = "re"), #+ ,
                                 family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                 data = dat_fin_counts_m_scaled)


m.gam.spei_tmp_re_slope45 <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 5) +
                                   # s(pairID, bs= 're') +
                                   s(spei3_lag2, pairID, bs = "re"), #+ ,
                                 family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                 data = dat_fin_counts_m_scaled)

# add interaction? - it is working!
m.gam.spei_tmp_re_slope45_int <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 5) +
                                       ti(spei3_lag2, veg_tmp, k = 4) +
                                   # s(pairID, bs= 're') +
                                   s(spei3_lag2, pairID, bs = "re"), #+ ,
                                 family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                 data = dat_fin_counts_m_scaled)

# add random slope again for temp

m.gam.spei_tmp_re_slope45_temp_int <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 5) +
                                       ti(spei3_lag2, veg_tmp, k = 4) +
                                       # s(pairID, bs= 're') +
                                       s(spei3_lag2, pairID, bs = "re") +
                                         s(veg_tmp, pairID, bs = "re"), #+ ,
                                     family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                     data = dat_fin_counts_m_scaled)



# add random interacept (pairs as ranndom effects)
# # the best one!!
m.gam.spei_tmp_re_slope45_temp_int2 <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                                             s(veg_tmp, k = 5) +
                                            ti(spei3_lag2, veg_tmp, k = 4) +
                                             s(pairID, bs= 're') +
                                            s(spei3_lag2, pairID, bs = "re") +
                                            s(veg_tmp, pairID, bs = "re"), #+ ,
                                          family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                          data = dat_fin_counts_m_scaled)


AIC(m.gam.spei_tmp_re_slope, 
    m.gam.spei_tmp_re_slope45_int, 
    m.gam.spei_tmp_re_slope45_temp_int, 
    m.gam.spei_tmp_re_slope45_temp_int2) # # This one is teh best and proceeed further with veg_tmp and spei3_lag2


plot(m.gam.spei_tmp_re_slope45_temp_int2, page = 1)




########## Test based on Dominik;s examples: ---------------------
# get teh feeling about how stable the fixed effects area - if they are not hindered by random effects

#The random effect is really dominating your variance explained… I am wondering if the random effect somehow clouds the fixed effects, although I think what you are doing makes a lot of sense. I would test the following alternatives:
  
#  -random effect  of trapspairs (eg allowing for different intercept - difernt starting point) - improves model a lot! Omit the random effect
#-	Change the random effect structure (random slope and intercept)
#-	Add site as fixed effect
#-	Use x and y coordinates instead of random effect

# Test random slopes


# test gams using final predictors: spei6_lag1 vs veg_tmp_lag2

m1 <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4), 
                  family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                  data = dat_fin_counts_m_scaled)

m2 <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4) + 
            ti(spei6_lag1,veg_tmp_lag2, k = 4), 
          family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

m2.rs <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4) + 
            ti(spei6_lag1,veg_tmp_lag2, k = 4) +
              s(pairID, bs= 're'), 
          family = nb(),
          method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)
# yes, adding random effects does imporve the AIC!

m3.rs.re.tmp <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4) + 
               ti(spei6_lag1,veg_tmp_lag2, k = 4) +
               s(pairID, bs= 're') +
                 s(veg_tmp_lag2,pairID, bs= 're'), 
             family = nb(),
             method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
             data = dat_fin_counts_m_scaled)

# random slope for tmp& spei: /> high collinearity!
m3.rs.re.tmp.spei <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4) + 
                      ti(spei6_lag1,veg_tmp_lag2, k = 4) +
                      s(pairID, bs= 're') +
                      s(veg_tmp_lag2,pairID, bs= 're') +
                        s(spei6_lag1,pairID, bs= 're'), 
                    family = nb(),
                    method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                    data = dat_fin_counts_m_scaled)

# random slope for only spei: - still high coll for spei, go for random slope for tmp
m3.rs.re.spei <- gam(sum_ips ~ s(spei6_lag1, k = 4) +  s(veg_tmp_lag2, k = 4) + 
                           ti(spei6_lag1,veg_tmp_lag2, k = 4) +
                           s(pairID, bs= 're') +
                          # s(veg_tmp_lag2,pairID, bs= 're') +
                           s(spei6_lag1,pairID, bs= 're'), 
                         family = nb(),
                         method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                         data = dat_fin_counts_m_scaled)


m3.rs.re.tmp545 <- gam(sum_ips ~ s(spei6_lag1, k = 5) +  s(veg_tmp_lag2, k = 4) + 
                      ti(spei6_lag1,veg_tmp_lag2, k = 5) +
                      s(pairID, bs= 're') +
                      s(veg_tmp_lag2,pairID, bs= 're'), 
                    family = nb(),
                    method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                    data = dat_fin_counts_m_scaled)


m3.rs.545 <- gam(sum_ips ~ s(spei6_lag1, k = 5) +  s(veg_tmp_lag2, k = 4) + 
                         ti(spei6_lag1,veg_tmp_lag2, k = 5) +
                         s(pairID, bs= 're'), 
                       family = nb(),
                       method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                       data = dat_fin_counts_m_scaled)



k.check(m3.rs.re.tmp)
check_collinearity(m3.rs.re.tmp)
AIC(m2.rs, m3.rs.re.tmp, m3.rs.re.tmp.spei,m3.rs.re.spei) # # This one is teh best and proceeed further with veg_tmp and spei3_lag2


plot(m3.rs.re.tmp, page = 1)



# try glmm????

library(nlme)

# Fit the GAMM model
m2.rs.gamm <- gamm(sum_ips ~ s(spei6_lag1, k = 4) + s(veg_tmp_lag2, k = 4) + 
                     ti(spei6_lag1, veg_tmp_lag2, k = 4), 
                   random = list(pairID = ~1, trapID = ~1, year = ~1), 
                   family = nb(), 
                   method = 'REML', 
                   data = dat_fin_counts_m_scaled)

library(gamm4)
m2.rs.gamm4 <- gamm4(sum_ips ~ s(spei6_lag1, k = 4) + s(veg_tmp_lag2, k = 4) + 
                       t2(spei6_lag1, veg_tmp_lag2, k = 4), 
                     random = ~(1|pairID/trapID) + (1|year), 
                     family = nb(), 
                     data = dat_fin_counts_m_scaled)


# Fit the model with glmmTMB
m2.rs.glmmTMB <- glmmTMB(sum_ips ~ splines::ns(spei6_lag1, df = 4) + splines::ns(veg_tmp_lag2, df = 4) + 
                           splines::ns(spei6_lag1, df = 4) * splines::ns(veg_tmp_lag2, df = 4) + 
                           (1 | pairID/trapID) + (1 | year), 
                         family = nbinom2, 
                         data = dat_fin_counts_m_scaled)


###### How much variance is explained by individual terms? -----------
# Full model
m.gam.spei_tmp_re_slope45_temp_int2 <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                                             s(veg_tmp, k = 5) +
                                             ti(spei3_lag2, veg_tmp, k = 4) +
                                             s(pairID, bs= 're') +
                                             s(spei3_lag2, pairID, bs = "re") +
                                             s(veg_tmp, pairID, bs = "re"), #+ ,
                                           family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                           data = dat_fin_counts_m_scaled)



#fin.m.counts <- m.gam.spei_tmp_re_slope45_temp_int2
# Fit the reduced model without veg_tmp
m_no_veg_tmp <- gam(sum_ips ~ s(spei3_lag2, k = 4) + 
                      ti(spei3_lag2, veg_tmp, k = 4) +
                      s(pairID, bs= 're') +
                      s(spei3_lag2, pairID, bs = "re"),
                    family = nb(), method = 'REML', data = dat_fin_counts_m_scaled)

# Fit the reduced model without spei3_lag2
m_no_spei3_lag2 <- gam(sum_ips ~ s(veg_tmp, k = 5) + 
                         s(pairID, bs= 're') +
                         ti(spei3_lag2, veg_tmp, k = 4) +
                         s(veg_tmp, pairID, bs = "re"),
                       family = nb(), method = 'REML', data = dat_fin_counts_m_scaled)

# Fit the reduced model without interaction term
m_no_interaction <- gam(sum_ips ~ s(spei3_lag2, k = 4) + s(veg_tmp, k = 5) +
                          s(spei3_lag2, pairID, bs = "re") +
                          s(veg_tmp, pairID, bs = "re") +
                          s(pairID, bs= 're'),
                        family = nb(), method = 'REML', data = dat_fin_counts_m_scaled)

# Calculate deviance explained for the full model
dev_full <- summary(fin.m.counts)$dev.expl

# Calculate deviance explained for the reduced models
dev_no_veg_tmp <- summary(m_no_veg_tmp)$dev.expl
dev_no_spei3_lag2 <- summary(m_no_spei3_lag2)$dev.expl
dev_no_interaction <- summary(m_no_interaction)$dev.expl

# Calculate deviance explained by each predictor
deviance_explained_veg_tmp <- dev_full - dev_no_veg_tmp
deviance_explained_spei3_lag2 <- dev_full - dev_no_spei3_lag2
deviance_explained_interaction <- dev_full - dev_no_interaction

# Print results
cat("Deviance explained by veg_tmp: ", deviance_explained_veg_tmp, "\n")
cat("Deviance explained by spei3_lag2: ", deviance_explained_spei3_lag2, "\n")
cat("Deviance explained by interaction (spei3_lag2 * veg_tmp): ", deviance_explained_interaction, "\n")





### Peak difference --------------------------------------------------

# test gam: use the save variables for spei and veg_tmp as for sum_ips
m1 <- gam(peak_diff ~ s(spei3_lag2, k = 4) + 
                                             s(veg_tmp, k = 5) +
                                             ti(spei3_lag2, veg_tmp, k = 4) +
                                             s(pairID, bs= 're') +
                                             s(spei3_lag2, pairID, bs = "re") +
                                             s(veg_tmp, pairID, bs = "re"), #+ ,
                                           family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
                                           data = dat_fin_counts_m_scaled)

# test gam: use the save variables for spei and veg_tmp as for sum_ips
m2 <- gam(peak_diff ~ s(spei3_lag2, k = 4) + 
            s(veg_tmp, k = 7) +
            ti(spei3_lag2, veg_tmp, k = 4) +
            s(pairID, bs= 're') +
            s(spei3_lag2, pairID, bs = "re") +
            s(veg_tmp, pairID, bs = "re"), #+ ,
          family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)


# remove non significant terms:

# test gam: use the save variables for spei and veg_tmp as for sum_ips
m2_1 <- gam(peak_diff ~ s(spei3_lag2, k = 4) + 
            s(veg_tmp, k = 7) +
            ti(spei3_lag2, veg_tmp, k = 4) +
            s(pairID, bs= 're'), # +
            #s(spei3_lag2, pairID, bs = "re") +
            #s(veg_tmp, pairID, bs = "re"), #+ ,
          family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)

# increased k in spei just makes collinearity explode!!! if spei, k = 5!
# keep k low for spei
m3 <- gam(peak_diff ~ s(spei3_lag2, k = 4) + 
            s(veg_tmp, k = 7) +
            ti(spei3_lag2, veg_tmp, k = 4) +
            #s(pairID, bs= 're') +
            s(spei3_lag2, pairID, bs = "re") +
            s(veg_tmp, pairID, bs = "re"), #+ ,
          family = nb(),method ='REML', # (use ML methos if wish to compare model with glmmTMB by AIC)  "REML", 
          data = dat_fin_counts_m_scaled)


check_collinearity(m3)
summary(m3)
plot(m3, page = 1, shade = T)
AIC(m1, m2, m2_1, m3)


fin.m.peak.diff <- m3


#### Eample: quick plotting! -----------------------------------------
fin.m <- m
#  m.spei12_no_re_slope47# fin.m.peak.diff
#!!! --------------------------------------------

# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = c("tmp_z_lag1", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

p_df <- as.data.frame(p3)
# test simple plot:
ggplot(p_df, aes(x = x, y = predicted)) +
  # geom_point(data = dat_fin_counts_m, aes(x = tmp_z_lag2, y = sum_ips, color = year_fact), 
  #            alpha = 0.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  ylim(0,75000) +
  #labs(x = "SPEI Lag 2", y = "Sum of Beetle Counts") +
  theme_classic2()


# test simple plot:
ggplot(p3, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
 # ylim(0,75000) +
  #labs(x = "SPEI Lag 2", y = "Sum of Beetle Counts") +
  theme_classic2()





fin.m.lagged <- m.spei12_no_re_slope47


# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m.lagged, terms = "veg_tmp_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.lagged, terms = "spei12_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.lagged, terms = c("veg_tmp_lag2", "spei12_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)


# test simple plot:
p.lagged.tmp <- ggplot(p3, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  ggtitle('Lagged temp')


ggarrange(p.current_tmp, p.lagged.tmp)

# try spei12 and veg_tmp
fin.m. <- m.spei12



# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m., terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m., terms = "spei12_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m., terms = c("veg_tmp", "spei12_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)


# test simple plot:
 ggplot(p1, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  ggtitle('current tmp spei 12')



#OLD:  best lags: veg_tmp_lag3 + spei3_lag3





m1 <- glmmTMB(peak_diff  ~ veg_tmp_lag3 + spei3_lag3,
              family = nbinom2,
              data = dat_fin_counts_m)

m2 <- glmmTMB(peak_diff ~ veg_tmp_lag3*spei3_lag3,
              family = nbinom2,
              data = dat_fin_counts_m)
# add raodom effects
m3 <- glmmTMB(peak_diff ~ veg_tmp_lag3*spei3_lag3 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

#allEffects(m3)
# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(peak_diff ~ poly(veg_tmp_lag3, 2) + poly(spei3_lag3,2) + veg_tmp_lag2:spei3_lag3 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)
# interaction of poly terms - tto crazy, difficult to interpret, stay with simple interaction
m5 <- glmmTMB(peak_diff ~ poly(veg_tmp_lag3, 2)*poly(spei3_lag3,2)+ (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

m6 <- glmmTMB(peak_diff ~ veg_tmp_lag3 + spei3_lag3 + veg_tmp_lag2:spei3_lag3 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

plot(dat_fin_counts_m$veg_tmp_lag3, dat_fin_counts_m$peak_diff)

plot(dat_fin_counts_m$spei3_lag3, dat_fin_counts_m$peak_diff)

simulateResiduals(m4, plot = T)

AIC(m1, m2, m3, m4, m5, m6)
summary(m6)




#simulateResiduals(m4, plot = T)

summary(m4)
AIC(m1, m2, m3, m4, m5)
#fin.m.peak.diff <- #m4


# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))
(vif_values <- performance::check_collinearity(m2))
(vif_values <- performance::check_collinearity(m3))


(vif_values <- performance::check_collinearity(m.gam77_int))
##### GAM Peak difference --------------------
# baseic plot
m.gam_0 <- gam(peak_diff ~ s(veg_tmp_lag3, k = 4) + 
                 s(spei3_lag3, k = 4) + 
                 ti(veg_tmp_lag3, spei3_lag3, k = 4),# +
              #   s(pairID, bs = "re") +  # random intercept: ifferent baseline for pairs
              #   s(veg_tmp_lag3, pairID, bs = "re") +  # random slope
              #   s(spei3_lag3, pairID, bs = "re"), # random slope
               family = nb(),
               method = 'REML',
               data = dat_fin_counts_m)

# add random elocations
m.gam_re <- gam(peak_diff ~ s(veg_tmp_lag3, k = 4) + 
                 s(spei3_lag3, k = 4) + 
                 ti(veg_tmp_lag3, spei3_lag3, k = 4) +
                    s(pairID, bs = "re"), # +  # random intercept: ifferent baseline for pairs
                 #   s(veg_tmp_lag3, pairID, bs = "re") +  # random slope
                 #   s(spei3_lag3, pairID, bs = "re"), # random slope
                 family = nb(),
               method = 'REML',
               data = dat_fin_counts_m)

# add random slopes
m.gam44 <- gam(peak_diff ~ s(veg_tmp_lag3, k = 4) + 
                s(spei3_lag3, k = 4) + 
                ti(veg_tmp_lag3, spei3_lag3, k = 4) +
                s(pairID, bs = "re") +  # random intercept: ifferent baseline for pairs
                s(veg_tmp_lag3, pairID, bs = "re") +  # random slope
                s(spei3_lag3, pairID, bs = "re"), # random slope
              family = nb(),
              method = 'REML',
              data = dat_fin_counts_m)

m.gam77 <- gam(peak_diff ~ s(veg_tmp_lag3, k = 7) + 
                                   s(spei3_lag3, k = 7) + 
                                   ti(veg_tmp_lag3, spei3_lag3, k = 4) +
                                   s(pairID, bs = "re") +  # random intercept: ifferent baseline for pairs
                                   s(veg_tmp_lag3, pairID, bs = "re") +  # random slope
                                   s(spei3_lag3, pairID, bs = "re"), # random slope
                                 family = nb(),
                                 method = 'REML',
                                 data = dat_fin_counts_m)

# add random slopes for interaction
m.gam77_int <- gam(peak_diff ~ s(veg_tmp_lag3, k = 7) + 
                 s(spei3_lag3, k = 7) + 
                 ti(veg_tmp_lag3, spei3_lag3, k = 4) +
                 s(pairID, bs = "re") +  # random intercept: ifferent baseline for pairs
                 s(veg_tmp_lag3, pairID, bs = "re") +  # random slope
                 s(spei3_lag3, pairID, bs = "re") + # random slope
                 s(veg_tmp_lag3, spei3_lag3, pairID, bs = "re"), # random slope
               family = nb(),
               method = 'REML',
               data = dat_fin_counts_m)


# High multicollinearity! 
# will scaling help?

# Centering and scaling the predictors
dat_fin_counts_m$veg_tmp_lag3_sc <- scale(dat_fin_counts_m$veg_tmp_lag3, center = TRUE, scale = TRUE)
dat_fin_counts_m$spei3_lag3_sc <- scale(dat_fin_counts_m$spei3_lag3, center = TRUE, scale = TRUE)

# Refit the model with centered and scaled predictors
m.gam77_int_sc <- gam(peak_diff ~ s(veg_tmp_lag3_sc, k = 7) + 
                     s(spei3_lag3_sc, k = 7) + 
                     ti(veg_tmp_lag3_sc, spei3_lag3_sc, k = 4) +
                     s(pairID, bs = "re") + 
                     s(veg_tmp_lag3_sc, pairID, bs = "re") + 
                     s(spei3_lag3_sc, pairID, bs = "re") + 
                     s(veg_tmp_lag3_sc, spei3_lag3_sc, pairID, bs = "re"),
                   family = nb(), 
                   method = 'REML',
                   data = dat_fin_counts_m)






# test multicoollinearity bwetween predictors ------------------------------


m1<- gam(peak_diff ~ dat_fin_counts_m)

cor(dat_fin_counts_m$veg_tmp, dat_fin_counts_m$spring_tmp, method = 'spearman') # cor ~ 0.9


cor(dat_fin_counts_m$veg_tmp_lag3, dat_fin_counts_m$spei3_lag3, method = 'spearman' )
# Check VIF for the predictors
vif_model <- gam(peak_diff ~ s(veg_tmp_lag3, k = 7) + 
                   s(spei3_lag3, k = 7), 
                 data = dat_fin_counts_m)

vif(vif_model)

#heck multicollinerarity from teh all predictors
(vif_values <- performance::check_collinearity(m.gam77_int))

AIC(m.gam_0, m.gam_re, m.gam44, m.gam77, m.gam77_int)
summary(m.gam77_int)
plot(m.gam77_int, shade = T, page = 1)

dat_fin_counts_m %>%
  ggplot(aes(x = veg_tmp_lag3,
             y = peak_diff)) + 
  geom_smooth() +
  geom_point()

fin.m <- m.gam77_int
#!!! --------------------------------------------
# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m, terms = "veg_tmp_lag3 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "spei3_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = c("veg_tmp_lag3", "spei3_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)
#p3 <- ggpredict(m.gam3, terms = c("veg_tmp", "spei3_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)



# test simple plot:
ggplot(p1, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  geom_point( data = dat_fin_counts_m,
              aes(
                  x = veg_tmp_lag3,
                  y = peak_diff ))




### Colonization DOY ------------------------------------------------------------




# lowest AIC: veg_tmp, spei3_lag3
m1 <- glmmTMB(tr_agg_doy  ~ veg_tmp + spei3_lag3,
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

m2 <- glmmTMB(tr_agg_doy ~ veg_tmp*spei3_lag3,
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)
# add raodom effects
m3 <- glmmTMB(tr_agg_doy ~ veg_tmp*spei3_lag3 + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(tr_agg_doy ~ poly(veg_tmp, 2) + poly(spei3_lag3,2) + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

# add poly terms, no interaction (issuees in model fitting otherwise)
m5 <- glmmTMB(tr_agg_doy ~ veg_tmp*poly(spei3_lag3,2) + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

# add poly terms, no interaction (issuees in model fitting otherwise)
m6 <- glmmTMB(tr_agg_doy ~ veg_tmp + poly(spei3_lag3,2) +veg_tmp:spei3_lag3 + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)


# add poly terms,with interaction
m7 <- glmmTMB(tr_agg_doy ~ poly(veg_tmp, 2)*poly(spei3_lag3,2) + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)



simulateResiduals(m3, plot = T)

summary(m7)
AIC(m1, m2, m3, m4, m5,m6, m7)
fin.m.agg <- m3

# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))

### Peak DOY ------------------------------------------------------------
# best lags: veg_tmp_lag2, spei3_lag2
m1 <- glmmTMB(tr_peak_doy  ~ veg_tmp_lag2 + spei3_lag2,
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

m2 <- glmmTMB(tr_peak_doy ~ veg_tmp_lag2*spei3_lag2,
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)
# add raodom effects
m3 <- glmmTMB(tr_peak_doy ~ veg_tmp_lag2*spei3_lag2 + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(tr_peak_doy ~ poly(veg_tmp_lag2, 2) + poly(spei3_lag2,2) + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

m5 <- glmmTMB(tr_peak_doy ~ veg_tmp_lag2 + poly(spei3_lag2,2) + (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

m6 <- glmmTMB(tr_peak_doy ~ veg_tmp_lag2 + poly(spei3_lag2,2) + veg_tmp_lag2:spei3_lag2+ (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

m7 <- glmmTMB(tr_peak_doy ~ veg_tmp_lag2 + spei3_lag2 + veg_tmp_lag2:spei3_lag2+ (1|pairID),
              family = beta_family(link = "logit"),
              data = dat_fin_counts_m)

AIC(m6, m7)
fin.m.peak <- m7

simulateResiduals(m5, plot = T)

summary(m5)
AIC(m1, m2, m3, m4, m5, m6)

# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))



### MOran's I LOG values (from log(beetle counts)): ----------------------------------
# export new df
# best lags: 
# veg_tmp0, spei3_2, sum_ips0, peak_doy0,peak_diff1, agg_doy0 

cols_skip <- c('trapID', 'pairID', 'year','Morans_I',  'Morans_I_log', 'Morans_I_log_cap')

dat_fin_Moran <-  
  dat_fin %>% 
  #dplyr::filter(Morans_I_log > 0 ) %>%
  dplyr::mutate(Morans_I_log_cap  = case_when(Morans_I_log >= 1 ~ 1,
                                              Morans_I_log < 1 & Morans_I_log >= -1 ~ Morans_I_log,
                                              Morans_I_log <= -1 ~ -1)) %>%
  dplyr::select(-wind_beetle) %>% 
  # data %>%
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei3_lag2   = lag(spei3, n = 2, default = NA),
         veg_tmp_lag2 = lag(veg_tmp, n = 2, default = NA),
         agg_doy_lag1   = lag(agg_doy, n = 2, default = NA),
         peak_doy_lag2   = lag(peak_doy, n = 2, default = NA),
         sum_ips_lag2   = lag(sum_ips, n = 2, default = NA),
         peak_diff_lag1   = lag(peak_diff, n = 1, default = NA)) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup() %>% 
  na.omit()
  
dat_fin_Moran_scaled <-
  dat_fin_Moran %>% 
  ungroup(.) %>% 
  # na.omit() %>% 
  dplyr::select(-all_of(cols_skip )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_fin_Moran %>% dplyr::select(all_of(cols_skip)), .)


summary(dat_fin_Moran_scaled)



# test with capped values (-1, 0, 1 range) -------------------------------------
# remove negative Morans I vals
dat_fin_Moran_scaled_positive <- dat_fin_Moran_scaled %>% 
  dplyr::filter(Morans_I_log_cap>0)


# put all predictors 
m1 <- glmmTMB(Morans_I_log_cap ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff_lag1 + population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)


#heck multicollinerarity from teh all predictors
(vif_values <- performance::check_collinearity(m1))
summary(m1)
# add climate intearction
m2 <- glmmTMB(Morans_I_log_cap ~ veg_tmp*spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff_lag1 +population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)

summary(m2)
# remove non significan terms and climate interaction,  add interaction with colonization
m3 <- glmmTMB(Morans_I_log_cap ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + 
                sum_ips + #peak_diff_lag1 +population_growth2 +
                (1|pairID)              ,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)



AIC(m1,m2, m3)
simulateResiduals(m3, plot = T)
plot(allEffects(m3))
#fin.m.moran <- m2

# remove all intercations
m4 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + population_growth2 +
                sum_ips,# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)

m5 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + population_growth2 + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)

# improved model with random effects: remove peak doy - non significant
m6 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + population_growth2 + #peak_doy + 
                sum_ips + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)


# try poly terms
m7 <- glmmTMB(Morans_I_log ~ poly(veg_tmp,2) + poly(spei3_lag2,2) + 
                poly(agg_doy,2) + poly(population_growth2, 2) + #peak_doy + 
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)

# simplify polyterms: remove non-significant
m8 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)


# interaction climate
m9 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled_positive)

# interaction climate
m10 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                 agg_doy + peak_doy + #population_growth2 +
                 #sum_ips + 
                 (1|pairID),# + 
               #peak_diff_lag1,
               family = tweedie,
               data = dat_fin_Moran_scaled_positive)



m11 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                 agg_doy + #peak_doy + #population_growth2 +
                 #sum_ips + 
                 (1|pairID),# + 
               #peak_diff_lag1,
               family = tweedie,
               data = dat_fin_Moran_scaled_positive)

AIC(#m1,
  #m2, m3, m4, m5, m6, m7, m8, m9, 
  #m9,  
  m10, m11)

plot(allEffects(m8))
fin.m.moran <- m11

simulateResiduals(m8, plot = T)

summary(m8)

# Morans I log (raw)----------------------------







# put all predictors 
m1 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff_lag1 + population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled)


#heck multicollinerarity from teh all predictors
(vif_values <- performance::check_collinearity(m1))
summary(m1)
# add climate intearction
m2 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff_lag1 +population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled)


# remove non significan terms and climate interaction,  add interaction with colonization
m3 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                veg_tmp:agg_doy + peak_doy + population_growth2 +
                sum_ips,# + 
                #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)
# remove all intercations
m4 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + population_growth2 +
                sum_ips,# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

m5 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + population_growth2 + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

# improved model with random effects: remove peak doy - non significant
m6 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + population_growth2 + #peak_doy + 
                sum_ips + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


# try poly terms
m7 <- glmmTMB(Morans_I_log ~ poly(veg_tmp,2) + poly(spei3_lag2,2) + 
                poly(agg_doy,2) + poly(population_growth2, 2) + #peak_doy + 
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

# simplify polyterms: remove non-significant
m8 <- glmmTMB(Morans_I_log ~ veg_tmp + spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


# interaction climate
m9 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

# interaction climate
m10 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                sum_ips + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)



m11 <- glmmTMB(Morans_I_log ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                #poly(sum_ips,2) + 
                 (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


AIC(#m1,
    #m2, m3, m4, m5, m6, m7, m8, m9, 
  m9,  
  m10, m11)

plot(allEffects(m8))
#fin.m.moran <- m10

simulateResiduals(m8, plot = T)

summary(m8)


#### predict MOrans_I (from raw betle counts) --------------------------------
dat_fin_Moran_scaled_trim <- dat_fin_Moran_scaled %>% 
  dplyr::filter(Morans_I > 0)

# put all predictors 
m1 <- glmmTMB(Morans_I ~ veg_tmp +veg_tmp_lag2 + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff + population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled_trim)


#heck multicollinerarity from teh all predictors
(vif_values <- performance::check_collinearity(m1))
summary(m1)
# add climate intearction
m2 <- glmmTMB(Morans_I ~ veg_tmp + veg_tmp_lag2*spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + peak_diff +population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled_trim)


m3 <- glmmTMB(Morans_I ~ veg_tmp_lag2*spei3_lag2,# + 
                #agg_doy + peak_doy + 
                #sum_ips + peak_diff +population_growth2,
              family = tweedie,
              data = dat_fin_Moran_scaled_trim)


# remove non significan terms and climate interaction,  add interaction with colonization
m3 <- glmmTMB(Morans_I ~ veg_tmp + spei3_lag2 + 
                veg_tmp:agg_doy + peak_doy + population_growth2 +
                sum_ips,# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)
# remove all intercations
m4 <- glmmTMB(Morans_I ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + population_growth2 +
                sum_ips,# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

m5 <- glmmTMB(Morans_I ~ veg_tmp + spei3_lag2 + 
                agg_doy + peak_doy + 
                sum_ips + population_growth2 + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

# improved model with random effects: remove peak doy - non significant
m6 <- glmmTMB(Morans_I ~ veg_tmp + spei3_lag2 + 
                agg_doy + population_growth2 + #peak_doy + 
                sum_ips + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


# try poly terms
m7 <- glmmTMB(Morans_I ~ poly(veg_tmp,2) + poly(spei3_lag2,2) + 
                poly(agg_doy,2) + poly(population_growth2, 2) + #peak_doy + 
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)

# simplify polyterms: remove non-significant
m8 <- glmmTMB(Morans_I ~ veg_tmp + spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


# interaction climate
m9 <- glmmTMB(Morans_I ~ veg_tmp*spei3_lag2 + 
                agg_doy + #peak_doy + population_growth2
                poly(sum_ips,2) + (1|pairID),# + 
              #peak_diff_lag1,
              family = tweedie,
              data = dat_fin_Moran_scaled)


AIC(m1,
  m2, m3, m4, m5, m6, m7, m8, m9)

plot(allEffects(m8))
#fin.m.moran <- m9

simulateResiduals(m8, plot = T)

summary(m9)






### RS beetle damage: 2 years lag of beetle data! ----------------------------------
## always use lag of 2 years, as little difference between the models, and we wish t have more observations

cols_skip <- c('trapID', 'pairID', 'year', 'wind_beetle')

dat_fin_RS_m <-  
  dat_fin %>% 
  
  #dplyr::select(-wind_beetle) %>% 
  # data %>%
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei3_lag2   = lag(spei3, n = 2, default = NA),
         veg_tmp_lag2 = lag(veg_tmp, n = 2, default = NA),
         agg_doy_lag2   = lag(agg_doy, n = 2, default = NA),
         peak_doy_lag2   = lag(peak_doy, n = 2, default = NA),
         sum_ips_lag2   = lag(sum_ips, n = 2, default = NA),
         peak_diff_lag2   = lag(peak_diff, n = 2, default = NA),
         Morans_I_log_lag2   = lag(Morans_I_log, n = 2, default = NA)) %>%
  #sum_ips_lag3   = lag(sum_ips, n = 3, default = NA)
  dplyr::mutate(population_growth     = (sum_ips - sum_ips_lag1) / sum_ips_lag1 * 100,
                population_growth2    = dplyr::lag(population_growth, n = 2, order_by = year)) %>%  # lag population growth by one more year
  
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup() %>% 
  na.omit()

dat_fin_RS_scaled <-
  dat_fin_RS_m %>% 
  ungroup(.) %>% 
  # na.omit() %>% 
  dplyr::select(-all_of(cols_skip )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_fin_RS_m %>% dplyr::select(all_of(cols_skip)), .)


summary(dat_fin_RS_scaled)
dim(dat_fin_RS_scaled)

# put all predictors 
m1 <- glmmTMB(wind_beetle ~ veg_tmp_lag2 + spei3_lag2 + agg_doy_lag2 + peak_doy_lag2 + 
                sum_ips_lag2 + peak_diff_lag2 + Morans_I_log_lag2 + population_growth2,
              family = nbinom2,
              data = dat_fin_RS_scaled)
summary(m1)
# Assuming your model is an 'lme4' or 'glmmTMB' object, check from additive model
(vif_values <- performance::check_collinearity(m1))



# only spei is significant here, add interaction with climate
m1.2 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + 
                agg_doy_lag2 + peak_doy_lag2 + 
                sum_ips_lag2 + peak_diff_lag2 + Morans_I_log_lag2 + population_growth2,
              family = nbinom2,
              data = dat_fin_RS_scaled)

# only spei is significant here, add interaction with climate
m1.2.1 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + 
                  agg_doy_lag2 + peak_doy_lag2, # + 
                 # sum_ips_lag2 + peak_diff_lag2 + 
                  # Morans_I_log_lag2 + population_growth2,
                family = nbinom2,
                data = dat_fin_RS_scaled)


m1.2.1.re <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + 
                    agg_doy_lag2 + peak_doy_lag2  +  (1|pairID),
                  # sum_ips_lag2 + peak_diff_lag2 + 
                  # Morans_I_log_lag2 + population_growth2,
                  family = nbinom2,
                  data = dat_fin_RS_scaled)

AIC(m3,m1.2.1,m1.2.1.re)

# keep only beetle pop predictors 
m2 <- glmmTMB(wind_beetle ~ #veg_tmp_lag2*spei3_lag2 + 
                agg_doy_lag2 + peak_doy_lag2 + 
                #peak_diff_lag2 + 
                Morans_I_log_lag2 +population_growth2+ 
                sum_ips_lag2
              ,# + ,
              family = nbinom2,
              data = dat_fin_RS_scaled)

m3 <- glmmTMB(wind_beetle ~ #veg_tmp_lag2*spei3_lag2 + 
                agg_doy_lag2 + peak_doy_lag2 + 
                #peak_diff_lag2 + 
                Morans_I_log_lag2 + #population_growth2+ 
                sum_ips_lag2 + (1|pairID) 
              ,# + ,
              family = nbinom2,
              data = dat_fin_RS_scaled)

# remove random effects
m4 <- glm.nb(wind_beetle ~ #veg_tmp_lag2*spei3_lag2 + 
                agg_doy_lag2, #+ # + peak_doy_lag2 + 
                #peak_diff_lag2 + 
                #Morans_I_log_lag2 + #population_growth2+ 
                #sum_ips_lag2,# + (1|pairID) 
             # ,# + ,
            #  family = nbinom2,
              data = dat_fin_RS_scaled)

AIC(m1, m1.2, m2, m3, m4)
r2(m4)

fin.m.RS <- m3
summary(fin.m.RS)

plot(allEffects(m4))
r2(m4)


























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
         Aggregation_DOY       = stringr::str_glue("{round(mean_agg_doy,1)}±{round(sd_agg_doy,1)}"),
         Peak_growth_DOY        = stringr::str_glue("{round(mean_peak_doy,1)}±{round(sd_peak_doy,1)}"),
         Peak_growth            = stringr::str_glue("{round(mean_peak_diff,1)}±{round(sd_peak_diff,1)}")) %>% 
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
      linewidth = 1  # Make the average line slightly thicker than individual lines
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
p_spagett_ips       <- plot_data_with_average(dat_fin, "sum_ips", lab_popul_level, '[a]')
p_spagett_agg_doy   <- plot_data_with_average(dat_fin, "agg_doy", lab_colonization_time, '[b]')
p_spagett_peak_doy  <- plot_data_with_average(dat_fin, "peak_doy", lab_peak_time, '[c]')
p_spagett_peak_diff <- plot_data_with_average(dat_fin, "peak_diff", lab_peak_growth, '[d]')

windows(7,6)
p_spagetti <- ggarrange(p_spagett_ips, 
          p_spagett_agg_doy, 
          p_spagett_peak_doy, 
          p_spagett_peak_diff, ncol = 2,nrow = 2, align = 'hv')

ggsave(filename = 'outFigs/Fig1.png', plot = p_spagetti, width = 7, height = 6, dpi = 300, bg = 'white')


# save models -------------------------------------------------------------

# 
# 
save( fin.m.agg,
      fin.m.counts,
      fin.m.peak,
      fin.m.peak.diff,
      file="outData/fin_models.Rdata")
# 
# 
# 


# Effect plots ------------------------------------------------------------


# define plotiing functions


# Create effect plot function with an additional argument to control y-axis labels
create_effect_plot <- function(data, line_color = "blue", x_title = "X-axis", 
                               y_title = "Y-axis",  my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  p <- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title,
         title = my_title,
         y = y_title) +
   # ylim(y_lim) +
  #  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.7, 2.7)) + # Set x-axis breaks and limits
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black")) +
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)

  return(p)
}




# change y-axis into days (from 0-1)
create_effect_Y_trans <- function(data, line_color = "blue", x_title = "X-axis", y_title = "Y-axis", #y_lim = c(0, 1), 
                                  my_title = '',
                                  x_annotate = 0, lab_annotate = "lab ann") {
  data$predicted <- (data$predicted * (doy.end - doy.start)) + doy.start  # Reverse transformation for y-axis
  data$conf.low  <- (data$conf.low * (doy.end - doy.start)) + doy.start
  data$conf.high <- (data$conf.high * (doy.end - doy.start)) + doy.start
  
  p<- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title, y = y_title, title = my_title) +
   # ylim(y_lim) +
  #  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2.7, 2.7)) + # Set x-axis breaks and limits
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black")) +
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
  
 }




# plot interactions
plot_effect_interactions <- function(data, temp_label, y_title,x_annotate = 0, lab_annotate = "lab ann") {
  #library(ggplot2)
  
 p<- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
    geom_line(aes(color = group, linetype = group), linewidth = 1) +
    labs(x = temp_label,
         y = y_title,
         fill = "SPEI\nlevels",
         color = "SPEI\nlevels",
         linetype = "SPEI\nlevels") +  # Fixed "y_title" to "y" for correct y-axis label argument
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          legend.position = 'bottom',
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(ncol = 1), 
           fill = guide_legend(ncol = 1),
           linetype = guide_legend(ncol = 1)) +
    annotate("text", x = x_annotate, y = Inf, 
             label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}


### Beetle counts -----------------------------------------------------------
# plot only significant terms

temp_label <- expression(paste("Temperature [", degree, "C]", sep=""))
spei_label <- 'SPEI'



# Assuming 'model' is your glm.nb model
summary(fin.m.counts)
p1 <- ggpredict(fin.m.counts, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts, terms = "spei3_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.counts, terms = c("veg_tmp", "spei3_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)
#p3 <- ggpredict(m.gam3, terms = c("veg_tmp", "spei3_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)


# change the values to allow y lables to fit 
p1$predicted <- p1$predicted/100
p1$conf.low  <- p1$conf.low/100
p1$conf.high <- p1$conf.high/100

p2$predicted <- p2$predicted/100
p2$conf.low  <- p2$conf.low/100
p2$conf.high <- p2$conf.high/100

p3$predicted <- p3$predicted/100
p3$conf.low  <- p3$conf.low/100
p3$conf.high <- p3$conf.high/100



p1.count <- 
  create_effect_plot(p1, line_color = "red", 
                               x_title = temp_label, 
                               y_title = "Counts [#*100]",
                     x_annotate = 13,
                     lab_annotate = "***"
                     #,  y_lim = c(80,800)
                               ) 
 
p2.count <- create_effect_plot(p2, line_color = "blue", 
                               x_title = "SPEI [dim.]", 
                               y_title = "Counts [#*100]", #,  y_lim = c(80,800)
                               x_annotate = 0,
                               lab_annotate = "***") 
p3.count <- plot_effect_interactions(p3, 
                                     temp_label = "temp",#''temp_label, 
                                     y_title = "Counts [#*100]",
                                     x_annotate = 13,
                                     lab_annotate = "n.s.") 

ggarrange(p1.count,p2.count, p3.count)



# test simple plot:
ggplot(p3, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) 




### DOY aggregation ---------------------------------------------------------
summary(fin.m.agg)
p1 <- ggpredict(fin.m.agg, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg, terms = "spei3_lag3 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.agg, 
                terms = c("veg_tmp", "spei3_lag3 [-1, 0, 1]"), 
                #type = "re", 
                allow.new.levels = TRUE)



p3$predicted <- (p3$predicted * (doy.end - doy.start)) + doy.start  # Reverse transformation for y-axis
p3$conf.low  <- (p3$conf.low * (doy.end - doy.start)) + doy.start
p3$conf.high <- (p3$conf.high * (doy.end - doy.start)) + doy.start


p1.agg <- create_effect_Y_trans(p1, 
                                line_color = "red", 
                                x_title = temp_label, 
                                y_title = "DOY",
                                x_annotate = 13,
                                lab_annotate = "***"#, y_lim = c(80,250)
                                )
p2.agg <- create_effect_Y_trans(p2, 
                                line_color = "blue", 
                                x_title = "SPEI [dim.]", 
                                y_title = "DOY",
                                x_annotate = 0,
                                lab_annotate = "*"#, y_lim = c(80,250)
                                )
p3.agg <- plot_effect_interactions(p3, temp_label = temp_label, 
                                         y_title = "DOY",
                                   x_annotate = 13,
                                   lab_annotate = "*")

p3.agg
###### PEak Effect plots ------------------------------------------------------------
summary(fin.m.peak)
# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak, terms = "veg_tmp_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak, terms = "spei3_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak, terms = c("veg_tmp_lag2" ,"spei3_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

p3$predicted <- (p3$predicted * (doy.end - doy.start)) + doy.start  # Reverse transformation for y-axis
p3$conf.low  <- (p3$conf.low * (doy.end - doy.start)) + doy.start
p3$conf.high <- (p3$conf.high * (doy.end - doy.start)) + doy.start



p1.peak <- create_effect_Y_trans(p1, line_color = "red", x_title = temp_label, 
                                 y_title = "DOY",
                                 x_annotate = 13,
                                 lab_annotate = "*"#, 
                                 #y_lim = c(130,200) 
                                 )

p2.peak <-create_effect_Y_trans(p2, 
                                line_color = "blue", 
                                x_title = "SPEI [dim.]", 
                                y_title = "DOY",
                                x_annotate = 0,
                                lab_annotate = "."#, 
                               #y_lim = c(130,200)
                               )
p3.peak <- plot_effect_interactions(p3, temp_label = temp_label, 
                                         y_title = "DOY",
                                    x_annotate = 13,
                                    lab_annotate = ".")




### effect plot peak dif ----------------------------------------------------
# Assuming 'model' is your glm.nb model
summary(fin.m.peak.diff)
p1 <- ggpredict(fin.m.peak.diff, terms = "veg_tmp_lag3 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff, terms = "spei3_lag3 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.diff, terms = c("veg_tmp_lag3", "spei3_lag3 [-1, 0, 1]"), allow.new.levels = T)

# change the values to allow y lables to fit 
p1$predicted <- p1$predicted/10
p1$conf.low  <- p1$conf.low/10
p1$conf.high <- p1$conf.high/10

p2$predicted <- p2$predicted/10
p2$conf.low  <- p2$conf.low/10
p2$conf.high <- p2$conf.high/10

p3$predicted <- p3$predicted/10
p3$conf.low  <- p3$conf.low/10
p3$conf.high <- p3$conf.high/10




p1.peak.diff <- create_effect_plot(p1, line_color = "red", x_title = temp_label, 
                                   y_title = "Counts [#*10]",
                                   x_annotate = 13,
                                   lab_annotate = "**"#, y_lim = c(220,3500)
                                   )
p2.peak.diff <-create_effect_plot(p2, line_color = "blue", x_title = "SPEI [dim.]", 
                                  y_title = "Counts [#*10]",
                                  x_annotate = 0,
                                  lab_annotate = "**"#, y_lim = c(220,1500)
                                  )
p3.peak.diff <- plot_effect_interactions(p3, temp_label = temp_label, 
                                         y_title = "Counts [#]",
                                         x_annotate = 13,
                                         lab_annotate = "n.s.")


### all Effect plots: ips vs climate --------------------------------------------------------

# put together individual plots for temp spei
p.count     <- ggarrange(p1.count, p2.count, ncol = 1)
p.agg       <- ggarrange(p1.agg, p2.agg, ncol = 1)
p.peak      <- ggarrange(p1.peak, p2.peak, ncol = 1)
p.peak.diff <- ggarrange(p1.peak.diff, p2.peak.diff, ncol = 1)

p.out.clim <- ggarrange(p.count, p.agg, p.peak, p.peak.diff, 
          ncol=4, nrow = 1 , align = 'hv', 
          font.label = list(size = 8, color = "black", face = "plain", family = NULL),
          labels = c( "[a] Population level",
                      "[b] Colonization timing",
                      "[c] Peak timing",
                      "[d] Peak growth"))
                      

windows(7,4)
(p.out.clim)

ggsave(filename = 'outFigs/Fig3.png', plot = p.out.clim, 
       width = 7, height = 4, dpi = 300, bg = 'white')


p.clim.int <- ggarrange(p3.count, p3.agg, p3.peak, p3.peak.diff, 
          ncol=2, nrow = 2 , align = 'hv',common.legend = TRUE, legend = 'right',
          font.label = list(size = 8, color = "black", face = "plain", family = NULL),
          labels = c(  "[a] Population level",
                       "[b] Colonization timing",
                       "[c] Peak timing",
                       "[d] Peak growth"))
windows(6,6)
p.clim.int
ggsave(filename = 'outFigs/Fig4.png', plot = p.clim.int, 
       width = 6, height = 6, dpi = 300, bg = 'white')




#### Effect plots RS ----------------------------------------------------------
# does not wokr!!!
# Assuming 'model' is your glm.nb model
summary(fin.m.RS)
p1 <- ggpredict(fin.m.RS, terms = "veg_tmp_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.RS, terms = "spei3_lag1    [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.RS, terms = c("veg_tmp_lag1", "spei3_lag1 [-1, 0,1]"), allow.new.levels = TRUE)
#p4 <- ggpredict(fin.m.RS, terms = "population_growth2 [all]", allow.new.levels = TRUE)


p1.RS <- create_effect_plot(p1, line_color = "red", 
                            x_title = "Temperature [dim.]", 
                            y_title = "Tree mortality\n[# pixels]", y_lim = c(0,20), 
                            show_y_axis = TRUE)
p2.RS <- create_effect_plot(p2, line_color = "blue", 
                            x_title = "SPEI [dim.]", 
                            y_title = "Tree mortality\n[# pixels]", y_lim = c(0,20), 
                            show_y_axis = TRUE)

p3.RS <- plot_effect_interactions(p3, temp_label = temp_label, 
                                         y_title = "Tree mortality\n[# pixels]")


p1.RS

### effect plots Sum IPS ----------------------------------------------------
p.effect.RS <- ggarrange(p1.RS,p2.RS,p3.RS,
                         #p4.RS, 
                         ncol = 2, nrow = 2, align = 'hv')

windows(7,7)
p.effect.RS

# effect plot MOran's I -------------------------------
summary(fin.m.moran)

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.moran , terms = "veg_tmp  [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.moran, terms = "spei3_lag2  [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.moran, terms = "agg_doy [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m.moran, terms = c("veg_tmp", "spei3_lag2 [-1,0,1]"), allow.new.levels = TRUE) 


p1.moran <- create_effect_plot(p1, line_color = "red", 
                               x_title = "Temperature [dim.]", 
                               y_title = "Local Moran's I", 
                               x_annotate = 0,
                               lab_annotate = "***"
                               #y_lim = c(0,1)
                               )
p2.moran <- create_effect_plot(p2, line_color = "blue", 
                               x_title = "SPEI [dim.]", 
                               y_title = "Local Moran's I",  
                               #y_lim = c(0,1)
                               x_annotate = 0,
                               lab_annotate = "*"
                               )
p3.moran <- create_effect_plot(p2, line_color = "grey50", 
                               x_title = "Colonization timing [dim.]", 
                               y_title = "Local Moran's I", # y_lim = c(0,1),
                               x_annotate = 0,
                               lab_annotate = "***")
p4.moran <- plot_effect_interactions(p4, temp_label = "Temperature [dim.]", 
                                     y_title = "Local Moran's I",
                                     x_annotate = 0,
                                     lab_annotate = "n.s.") #+
p4.moran.no.leg <- p4.moran + theme(legend.position = 'none')

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
                            p1.moran, p2.moran,
                            p3.moran, 
                            p4.moran.no.leg, my_legend, 
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
                                        "   ",
                                        "[e]"))


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
sjPlot::tab_model(fin.m.counts,    file = "outTable/gam_counts2.doc")
sjPlot::tab_model(fin.m.agg,       file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak,      file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff, file = "outTable/model_peak_diff.doc")
sjPlot::tab_model(fin.m.moran,     file = "outTable/model_moran.doc")
sjPlot::tab_model(fin.m.RS,        file = "outTable/model_RS_3.doc")

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




