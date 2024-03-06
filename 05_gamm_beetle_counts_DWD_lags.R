


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

# add MOran's I values, transform the agg doy and peak doy between 0-1
dat_fin <-   dat_fin %>% 
  left_join(lisa_merged_df, 
                       by = join_by(trapID, year, sum_ips)) %>%
  mutate(tr_agg_doy   = (agg_doy - 60) / (304 - 60),
         tr_peak_doy  = (peak_doy - 60) / (304 - 60))


RS_simple <- df_RS_out %>% 
  dplyr::select(trapID, year, wind_beetle)

dat_fin <-   dat_fin %>% 
  left_join(RS_simple, 
            by = join_by(trapID, year)) 
  


# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dat_fin$tr_agg_doy  <- pmin(pmax(dat_fin$tr_agg_doy, 1e-4),  1 - 1e-4)
dat_fin$tr_peak_doy <- pmin(pmax(dat_fin$tr_peak_doy, 1e-4), 1 - 1e-4)

  


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
#  left_join(lisa_merged_df, 
#            by = join_by(trapID, year, sum_ips)) %>%
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

head(dat_fin)



# Find teh most influential time lag :  0,1,2,3----------------------------------------------



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
    predictor <- predictors[predictor_idx]
    
    # Fit and store AIC for lag0 model first to use as a baseline for comparison
    lagged_var_name_lag0 <- paste0("lag0_", predictor)
    formula_lag0 <- as.formula(paste(dependent_var, "~", lagged_var_name_lag0))
    model_lag0 <- glm.nb(formula_lag0, data = data_complete_cases, na.action = na.exclude)
    aic_lag0[predictor_idx] <- AIC(model_lag0)
    
    for (lag in lags) {
      lagged_var_name <- paste0("lag", lag, "_", predictor)
      formula <- as.formula(paste(dependent_var, "~", lagged_var_name))
      model <- glm.nb(formula, data = data_complete_cases, na.action = na.exclude)
      aic_value <- AIC(model)
      
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
        Delta_AIC = delta_aic
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


# Specify predictors and lags
# List of predictors
predictors_beetle_dyn <- c("veg_tmp",  "spei3")
predictors_Morans     <- c("veg_tmp",  "spei3", "sum_ips", "peak_doy", "peak_diff", "agg_doy")
predictors_RS         <- c("sum_ips", "peak_doy", "peak_diff", "agg_doy", "Morans_I", "Morans_I_log")
#dependent_var <- 'sum_ips'
lags <- 0:3
#data <- dat_fin # Assuming dat_fin is your dataset

dat_fin <- dat_fin %>% 
    group_by(trapID) %>%
      arrange(year, .by_group = TRUE) %>%
      mutate(sum_ips_lag1   = lag(sum_ips, n = 3, default = NA),
             population_growth     = (sum_ips - sum_ips_lag1) / sum_ips_lag1 * 100,
             population_growth2    = dplyr::lag(population_growth, n = 1, order_by = year))# %>%  # lag population growth by one more year
  

# specify individual datasets, to limit the NAs values and to scale them
dat_fin_no_RS <- dat_fin %>% 
  dplyr::select(!wind_beetle)

cols_skip <- c('trapID', 'pairID', 'Morans_I_log')

# export new df
dat_fin_Moran <-  dat_fin %>% 
#  dplyr::filter(Morans_I_log > 0 & year %in% 2018:2021) %>%
  dplyr::select(-wind_beetle)

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


  
# Get lag models, validated by AIC --------------------------------------------
model_lag_sum_ips <- fit_lags_negbin(dat_fin_no_RS, predictors  = predictors_beetle_dyn, 
                                                 lags, 
                                     dependent_var = "sum_ips")
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


model_lag_Morans  <- fit_lags_tweedie_test(dat_fin_Moran_scaled, 
                                 predictors  =  predictors_Morans, 
                                 lags, 
                                 dependent_var = "Morans_I_log")



# simplify lags stats & export tables ----------------------------------------------------
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






# rerun models based on teh most important lags!! -------------------------------
# follow rules: test for intividual effects, poly tersm and interactions of teh most 
# imortant redictors
# for the most important model
# report


# IPS_counts test based on teh most important lags: lag 2 for spei3 and veg_tmp --------------
dat_fin_counts_m <- dat_fin %>% 
  dplyr::select(-wind_beetle) %>% 
  # data %>%
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei3_lag2 = lag(spei3, n = 2, default = NA),
         spei3_lag3 = lag(spei3, n = 3, default = NA),
         veg_tmp_lag2 = lag(veg_tmp, n = 2, default = NA),
         veg_tmp_lag3 = lag(veg_tmp, n = 3, default = NA)) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>% 
  ungroup() %>% 
  na.omit()

#data_filtered <- data %>% 
 # na.omit() # This removes rows with any NAs across all columns

### SUM_IPS -----------------------------------------

m1 <- glmmTMB(sum_ips ~ veg_tmp + spei3_lag2,
              family = nbinom2,
              data = dat_fin_counts_m)
# veg_tmp s slightly better then veg_tmp_lag2, proceed with veg_tmp
m2 <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2,
              family = nbinom2,
              data = dat_fin_counts_m)
# add raodom effects
m3 <- glmmTMB(sum_ips ~ veg_tmp*spei3_lag2 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

allEffects(m4)
fin.m.counts <- m4
# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + poly(spei3_lag2,2) + veg_tmp:spei3_lag2 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

m5 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + poly(spei3_lag2,2) + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)


simulateResiduals(m3, plot = T)

summary(m4)
AIC(m1, m2, m3, m4, m5)

# Chek for multicollinearity
library(performance)

# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))


### Peak difference --------------------------------------------------
# best lags: veg_tmp_lag3 + spei3_lag3

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
# interaction of poly terms
m5 <- glmmTMB(peak_diff ~ poly(veg_tmp_lag3, 2)*poly(spei3_lag3,2)+ (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)





simulateResiduals(m5, plot = T)

summary(m4)
AIC(m1, m2, m3, m4, m5)
fin.m.peak.diff <- m5


# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))


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

summary(m6)
AIC(m1, m2, m3, m4, m5,m6, m7)
fin.m.agg <- m6

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


fin.m.peak <- m5

simulateResiduals(m5, plot = T)

summary(m5)
AIC(m1, m2, m3, m4, m5, m6)

# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m1))



### MOran's I values: 3 years lag of beetle data! ----------------------------------
# export new df
# best lags: 
# veg_tmp0, spei3_2, sum_ips0, peak_doy0,peak_diff1, agg_doy0 

cols_skip <- c('trapID', 'pairID', 'year', 'Morans_I_log')

dat_fin_Moran <-  
  dat_fin %>% 
  #dplyr::filter(Morans_I_log > 0 ) %>%
  dplyr::select(-wind_beetle) %>% 
  # data %>%
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei3_lag2   = lag(spei3, n = 2, default = NA),
         #veg_tmp_lag2 = lag(veg_tmp, n = 2, default = NA),
         # agg_doy_lag1   = lag(agg_doy, n = 2, default = NA),
         #peak_doy_lag2   = lag(peak_doy, n = 2, default = NA),
        # sum_ips_lag2   = lag(sum_ips, n = 2, default = NA),
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


AIC(m1,m2, m3, m4, m5, m6, m7, m8)

plot(allEffects(m8))
fin.m.moran <- m8

simulateResiduals(m8, plot = T)

summary(m8)





### RS beetle damage: 3 years lag of beetle data! ----------------------------------
## always use lag of 3 years, as little difference between the models, and we wish t have more observations

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
         sum_ips_lag1   = lag(sum_ips, n = 1, default = NA),
         peak_diff_lag2   = lag(peak_diff, n = 2, default = NA),
         Morans_I_log_lag2   = lag(Morans_I_log, n = 2, default = NA)) %>%
  #sum_ips_lag3   = lag(sum_ips, n = 3, default = NA)
  dplyr::mutate(population_growth     = (sum_ips - sum_ips_lag1) / sum_ips_lag1 * 100,
                population_growth2    = dplyr::lag(population_growth, n = 1, order_by = year)) %>%  # lag population growth by one more year
  
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

# add intercation between climate
m2 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + agg_doy_lag2 + peak_doy_lag2 + 
                peak_diff_lag2 + Morans_I_log_lag2 +population_growth2+ 
                sum_ips_lag2,# + ,
              family = nbinom2,
              data = dat_fin_RS_scaled)
AIC(m1, m2, m3)
# remove teh least important predicors (based on estimate)
m3 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + agg_doy_lag2 + #peak_doy_lag2 + 
                  peak_diff_lag2 + Morans_I_log_lag2 +population_growth2 +
                sum_ips_lag2,# + ,
                family = nbinom2,
                data = dat_fin_RS_scaled)

# test interation between temp and beetle population
m4 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + agg_doy_lag2 + #peak_doy_lag2 + 
                peak_diff_lag2 + Morans_I_log_lag2 +population_growth2 +
                sum_ips_lag2 + veg_tmp_lag2:sum_ips_lag2,
              family = nbinom2,
              data = dat_fin_RS_scaled)
# interaction between  temp and agg
m5 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + agg_doy_lag2 + #peak_doy_lag2 + 
                peak_diff_lag2 + Morans_I_log_lag2 +population_growth2 +
                sum_ips_lag2 + veg_tmp_lag2:agg_doy_lag2,
              family = nbinom2,
              data = dat_fin_RS_scaled)

# no interactions are important (besides the climate). add random effect)

m6 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + agg_doy_lag2 + #peak_doy_lag2 + 
                peak_diff_lag2 + Morans_I_log_lag2 +population_growth2 +
                sum_ips_lag2,# + (1|pairID),
              family = nbinom2,
              data = dat_fin_RS_scaled)
# do not use random effects




AIC(m1, m2, m3, m4, m5, m6)
summary(m6)

plot(allEffects(m5))
summary(m6)
fin.m.RS <- m6

simulateResiduals(m6, plot = T)



























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






dw_result <- dwtest(residuals(m9.tmb) ~ fitted(m9.tmb))
(dw_result)


acf(residuals(m.agg12.3))

(dw_result <- dwtest(residuals(m.agg12.3) ~ fitted(m.agg12.3)))


summary(m.agg12.3)  # is teh winner!! m.agg12.3



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
pred_interaction <- ggpredict(m3, terms = c("veg_tmp_lag2", "spei3_lag2"))
pred_spei <- ggpredict(m3, terms = c("spei3_lag2"))
pred_veg_tmp <- ggpredict(m3, terms = c("veg_tmp_lag2"))


p1 <- ggplot(pred_spei, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(linewidth = 1.2) +
  ylim(0, 80000) +
  labs(x = "SPEI [z score]", y = "Beetle population levels [counts]") +
  theme(legend.position = 'none') +
  theme_classic()

p2 <- ggplot(pred_veg_tmp, aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(linewidth = 1.2) +
  ylim(0, 80000) +
  labs(x = "Seasonal Temperature [c]", y = "Beetle population levels [counts]") +
  theme(legend.position = 'none') +
  theme_classic()


p3 <- ggplot(pred_interaction, aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1.2) +
  ylim(0, 80000) +
  labs(x = "Seasonal temperature [C]", y = "Beetle population levels [counts]") +
  theme(legend.position = 'bottom') +
  theme_classic()

# Arrange the plots side by side with the same size
p.effect.counts <- ggarrange(p1, p2, p3, ncol = 3, align = "hv")
(p.effect.counts)













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




