


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
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
library('cowplot')
library('tidyr')
library("knitr")

library('readr')

library(ggeffects)

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

# test for autocorrelation
library(lmtest)

# colors
library(RColorBrewer)

library(knitr)   # for table outputs


# Define time variables:  Vars
#doy.start  =  91 # April 1st# 60  # March 1st, 
#doy.end    = 273  # SEpt 30  304 # Oct 30
#veg.period = doy.start:doy.end



# plot labels: 
lab_popul_level       = "Population level [#]"
lab_colonization_time = "Colonization timing [DOY]"
lab_peak_time         = "Peak timing [DOY]"
lab_peak_growth       = "Peak growth [#]"




##### read data ----------------------------------------------------------------
load(file=   "outData/final_table.Rdata")
load(file =  "outData/buffers.Rdata")  # df_RS_out
load(file =  "outData/lisa.Rdata")     # read LISA Moran's I stats



# prepare data RS --------------------------------------------
df_RS_out <- df_RS_out %>% 
  dplyr::select(-c(globalid, id)) %>% 
  dplyr::rename(trapID = falsto_name)

# change column names
lisa_merged_df <- lisa_merged_df %>% 
  dplyr::rename(trapID = falsto_name,
                sum_ips = sum_beetle) %>% 
  dplyr::select(c(year, trapID, sum_ips, 
                  Morans_I, # calculated from beetles sums 
                  Morans_I_log)) # calculated from log(beetles sum)

unique(lisa_merged_df$trapID)
unique(df_RS_out$trapID)

# add MOran's I values, transform the agg doy and peak doy between 0-1
dat_fin <-   dat_fin %>% 
  left_join(lisa_merged_df, 
            by = join_by(trapID, year, sum_ips)) %>%
  mutate(tr_agg_doy   = (agg_doy - doy.start) / (doy.end - doy.start),
         tr_peak_doy  = (peak_doy - doy.start) / (doy.end - doy.start))

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dat_fin$tr_agg_doy  <- pmin(pmax(dat_fin$tr_agg_doy, 1e-4),  1 - 1e-4)
dat_fin$tr_peak_doy <- pmin(pmax(dat_fin$tr_peak_doy, 1e-4), 1 - 1e-4)



RS_simple <- df_RS_out %>% 
  dplyr::select(trapID, year, wind_beetle)

dat_fin <-   dat_fin %>% 
  left_join(RS_simple, 
            by = join_by(trapID, year)) 
  




# calculate population growth
dat_fin <- dat_fin %>% 
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(sum_ips_lag1          = lag(sum_ips, n = 1, default = NA),
         population_growth     = (sum_ips - sum_ips_lag1) / sum_ips_lag1 * 100,
         population_growth2    = dplyr::lag(population_growth, n = 1, order_by = year))# %>%  # lag population growth by one more year


View(dat_fin2)

# !!!what is teh problem?? and why is teh lag = 3??

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
predictors_spei       <- c("spei1", "spei3","spei6", "spei12", "spei24")
predictors_beetle_dyn <- c("spring_tmp", "veg_tmp",  "spei3")
predictors_Morans     <- c("spring_tmp", "veg_tmp",  "spei3", "sum_ips", "peak_doy", "peak_diff", "agg_doy")
predictors_RS         <- c("spring_tmp", "veg_tmp", "spei3", "sum_ips", "peak_doy", "peak_diff", "agg_doy", "Morans_I", "Morans_I_log")
#dependent_var <- 'sum_ips'
lags <- 0:3
#data <- dat_fin # Assuming dat_fin is your dataset

# specify individual datasets, to limit the NAs values and to scale them
dat_fin_no_RS <- dat_fin %>% 
  dplyr::select(!wind_beetle)

cols_skip <- c('trapID', 'pairID', 'Morans_I_log', 'Morans_I')

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


  
### Get lag models, validated by AIC --------------------------------------------
model_lag_sum_ips <- fit_lags_negbin(dat_fin_no_RS, predictors  = predictors_beetle_dyn, 
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






# rerun models based on teh most important lags!! -------------------------------
# follow rules: test for intividual effects, poly tersm and interactions of teh most 
# imortant redictors
# for the most important model
# report


# IPS_counts test based on teh most important lags: lag 2 for spei3 and veg_tmp --------------
dat_fin_counts_m <- 
  dat_fin %>% 
  dplyr::select(-wind_beetle) %>% 
  # data %>%
  group_by(trapID) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(spei3_lag2 = lag(spei3, n = 2, default = NA),
         spei3_lag3 = lag(spei3, n = 3, default = NA),
         veg_tmp_lag2 = lag(veg_tmp, n = 2, default = NA),
         veg_tmp_lag3 = lag(veg_tmp, n = 3, default = NA)) %>%
  mutate(peak_diff  = as.integer(peak_diff )) %>%
 # View()
  ungroup() %>%
  #View()
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

# add poly terms, no interaction (issuees in model fitting otherwise)
m4 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + spei3_lag2 + veg_tmp:spei3_lag2 + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

m5 <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) + poly(spei3_lag2,2) + (1|pairID),
              family = nbinom2,
              data = dat_fin_counts_m)

#allEffects(m4)
fin.m.counts <- m4

simulateResiduals(m4, plot = T)

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
fin.m.peak.diff <- m4


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
summary(m2)
# Assuming your model is an 'lme4' or 'glmmTMB' object, check from additive model
(vif_values <- performance::check_collinearity(m1))



# only spei is significant here, add interaction with climate
m1.2 <- glmmTMB(wind_beetle ~ veg_tmp_lag2*spei3_lag2 + 
                agg_doy_lag2 + peak_doy_lag2 + 
                sum_ips_lag2 + peak_diff_lag2 + Morans_I_log_lag2 + population_growth2,
              family = nbinom2,
              data = dat_fin_RS_scaled)


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
AIC(m1, m2, m3)




fin.m.RS <- m3
summary(m3)

plot(allEffects(m2))



























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
         Colonization_DOY       = stringr::str_glue("{round(mean_agg_doy,1)}±{round(sd_agg_doy,1)}"),
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
                     lab_annotate = "**"
                     #,  y_lim = c(80,800)
                               ) 
 
p2.count <- create_effect_plot(p2, line_color = "blue", 
                               x_title = "SPEI [dim.]", 
                               y_title = "Counts [#*100]", #,  y_lim = c(80,800)
                               x_annotate = 0,
                               lab_annotate = ".") 
p3.count <- plot_effect_interactions(p3, 
                                     temp_label = temp_label, 
                                     y_title = "Counts [#*100]",
                                     x_annotate = 13,
                                     lab_annotate = "*") 

(p3.count)




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
sjPlot::tab_model(fin.m.counts,    file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg,       file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak,      file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff, file = "outTable/model_peak_diff.doc")
sjPlot::tab_model(fin.m.moran,     file = "outTable/model_moran.doc")
sjPlot::tab_model(fin.m.RS,        file = "outTable/model_RS.doc")

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




