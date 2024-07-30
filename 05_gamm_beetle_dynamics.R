


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
  na.omit()

fwrite(avg_data, 'outTable/fin_tab_avg.csv')


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





AIC(m.re$lme, m.xy$lme )
# compare the two: add xy or add pairs as random effects? what is better?

m<-m1$gam
summary(m)
plot(m, page = 1)

# automate models over all dependent variables and over increasing model complexity  -----------------------------------------
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
avg_data_filt <- avg_data %>% 
  dplyr::filter(!pairID %in% pair_outliers )

# test year as factor 
avg_data_filt$f_year <- as.factor(avg_data_filt$year)



nrow(avg_data_filt) # 518
nrow(avg_data) # 549

m.counts.tw <- gamm(sum_ips ~ s(year, k = 6) +
                         s(tmp_z_lag1, k = 5) +
                         s(spei_lag2, k = 8) + 
                         te(tmp_z_lag1, spei_lag2, k = 20) + 
                        # s(x, y, bs = 'gp', k = 70) + 
                         s(pairID, bs = 're'),
                       data = avg_data_filt, 
                       family = tw,
                       correlation = corAR1(form = ~ year | pairID))



m.counts.tw.test.f <- gamm(sum_ips ~ f_year +
                           s(tmp_lag1, k = 5) +
                           s(spei_lag2, k = 8) + 
                           te(tmp_lag1, spei_lag2, k = 20) + 
                           # s(x, y, bs = 'gp', k = 70) + 
                           s(pairID, bs = 're'),
                         data = avg_data_filt, 
                         family = tw,
                         correlation = corAR1(form = ~ year | pairID))

# add year as random
m.counts.tw.test.f.rndm <- gamm(sum_ips ~ s(f_year,bs = 're')  +
                             s(tmp_z_lag1, k = 5) +
                             s(spei_lag2, k = 8) + 
                             te(tmp_z_lag1, spei_lag2, k = 10) + 
                             # s(x, y, bs = 'gp', k = 70) + 
                             s(pairID, bs = 're'),
                           data = avg_data_filt, 
                           family = tw,
                           correlation = corAR1(form = ~ year | pairID))


# add year as random
m.counts.tw.test.f.rndm2 <- gamm(sum_ips ~ s(f_year,bs = 're')  +
                                  s(tmp_lag1, k = 5) +
                                  s(spei_lag2, k = 8) + 
                                  te(tmp_lag1, spei_lag2, k = 10) + 
                                  # s(x, y, bs = 'gp', k = 70) + 
                                  s(pairID, bs = 're'),
                                data = avg_data_filt, 
                                family = tw,
                                correlation = corAR1(form = ~ year | pairID))

# does not converge!! if adding both year as continuous and as factor
m.counts.tw.test.f.rndm3 <- gamm(sum_ips ~ s(f_year,bs = 're')  +
                                   s(year, k = 6) +
                                   s(tmp_lag1, k = 5) +
                                   s(spei_lag2, k = 8) + 
                                   te(tmp_lag1, spei_lag2, k = 10) + 
                                   # s(x, y, bs = 'gp', k = 70) + 
                                   s(pairID, bs = 're'),
                                 data = avg_data_filt, 
                                 family = tw,
                                 correlation = corAR1(form = ~ year | pairID))



AIC(m.counts.tw, m.counts.tw.test,m.counts.tw.test.f,m.counts.tw.test.f.rndm, m.counts.tw.test.f.rndm2)





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


summary(fin.m)









###### MORANS'I TW  ------------------------------------
#avg_data_filt <- avg_data %>% 
#  dplyr::filter(!pairID %in% pair_outliers )
avg_data_moran <- avg_data %>% 
 # mutate(Morans_I_log = ifelse(Morans_I_log > 1, 1, Morans_I_log)) %>% 
  group_by(pairID) %>% 
  # get lags of beetle indicators
  arrange(year, .by_group = TRUE) %>%
  mutate(sum_ips_lag1 = lag(sum_ips , n = 1, default = NA),
         sum_ips_lag2 = lag(sum_ips , n = 2, default = NA),
         peak_diff_lag1 = lag(peak_diff , n = 1, default = NA),
         peak_diff_lag2 = lag(peak_diff , n = 2, default = NA),
         tr_agg_doy_lag1 = lag(tr_agg_doy , n = 1, default = NA),
         tr_agg_doy_lag2 = lag(tr_agg_doy , n = 2, default = NA),
         tr_peak_doy_lag1 = lag(tr_peak_doy , n = 1, default = NA),
         tr_peak_doy_lag2 = lag(tr_peak_doy , n = 2, default = NA)) %>% 
  na.omit() %>% 
  # convert to log values for facter calculation
  mutate(sum_ips = log(sum_ips),
         sum_ips_lag1 = log(sum_ips_lag1),
         sum_ips_lag2 = log(sum_ips_lag2),
         peak_diff    = log(peak_diff),
         peak_diff_lag1 = log(peak_diff_lag1),
         peak_diff_lag2 = log(peak_diff_lag2)) %>% 
  mutate(f_year = factor(year))


# remove teh crazy values, and outliers for boths
avg_data_moran_sub <- avg_data_moran %>% 
  dplyr::filter(Morans_I_log < 1.5 & Morans_I_log > 0) %>% 
  dplyr::filter(sum_ips > 8.5)

summary(avg_data_moran_sub$sum_ips)

fwrite(avg_data_moran_sub, 'outTable/input_morans_model.csv')


####### Find predicors &  lags: MOrans  -------------
dependent_moran <-  c("Morans_I_log")

# list predictors to test
selected_predictors <- c('sum_ips', 'sum_ips_lag1','sum_ips_lag2',
                         #'log_sum_ips', 'log_sum_ips_lag1','sum_ips_lag2',
                         'tr_agg_doy'   , 'tr_agg_doy_lag1', 'tr_agg_doy_lag2' ,
                         'tr_peak_doy', 'tr_peak_doy_lag1','tr_peak_doy_lag2',
                         'peak_diff', 'peak_diff_lag1',  'peak_diff_lag2'
) 

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
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 8)"))
    print(formula)
    #print(formula)
    model <- gamm(formula,  
                  family = tw,
                  method = 'REML',  
                  data = avg_data_moran_sub,
                  correlation = corAR1(form = ~ year | pairID))
    
    # Extract model summary
    model_summary <- summary(model$gam)
    
    # Store the AIC value and deviance explained
    model_metrics_moran <- rbind(model_metrics_moran, data.frame(Predictor = pred, 
                                                           Dependent = dep, 
                                                           AIC = AIC(model$lme), 
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


fin.m.moran


sjPlot::model(fin.m.moran,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_lag_predictors_moran.doc",
               digits = 1) 

### END ---------------------------

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

m.moran.tw1 <- gamm(Morans_I_log  ~ #s(year, k = 5) +
                      #s(tmp_z_lag1, k = 5) +
                      #s(spei_z_lag2, k = 4) + 
                     s(sum_ips, k =5) +
                      #te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      # s(x, y, bs = 'gp', k = 70) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))



m.moran.tw2 <- gamm(Morans_I_log  ~ #s(year, k = 5) +
                      #s(tmp_z_lag1, k = 5) +
                      #s(spei_z_lag2, k = 4) + 
                      s(sum_ips, by = f_year, k =5) +
                      #te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      # s(x, y, bs = 'gp', k = 70) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))

m.moran.tw3 <- gamm(Morans_I_log  ~ s(year, k = 5) +
                      #s(tmp_z_lag1, k = 5) +
                      #s(spei_z_lag2, k = 4) + 
                      s(sum_ips, k =5) +
                      #te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      # s(x, y, bs = 'gp', k = 70) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))

# add climate 

m.moran.tw4 <- gamm(Morans_I_log  ~ s(year, k = 5) +
                      #s(tmp_z_lag1, k = 5) +
                      s(spei_z_lag2, k = 4) + 
                      s(sum_ips, k =10) +
                      #te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                       s(x, y, bs = 'gp', k = 50) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))


m.moran.tw5 <- gamm(Morans_I_log  ~ s(year, k = 5) +
                      s(tmp_z_lag1, k = 5) +
                      s(spei_z_lag2, k = 4) + 
                      s(sum_ips, k =10) +
                      #te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      s(x, y, bs = 'gp', k = 50) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))

m.moran.tw6 <- gamm(Morans_I_log  ~ s(year, k = 5) +
                      s(tmp_z_lag1, k = 5) +
                      s(spei_z_lag2, k = 4) + 
                      s(sum_ips, k =10) +
                      te(tmp_z_lag1, spei_z_lag2, k = 20) + 
                      s(x, y, bs = 'tp', k = 50) + 
                      s(pairID, bs = 're'),
                    data = avg_data_moran_sub, 
                    family = tw,
                    correlation = corAR1(form = ~ year | pairID))

AIC(m.moran.tw2$lme, m.moran.tw1$lme, m.moran.tw3$lme,m.moran.tw4$lme,
    m.moran.tw5$lme,m.moran.tw6$lme)

summary(avg_data$Morans_I_log)
appraise(m.moran.tw6$gam)
summary(m.moran.tw6$gam)
k.check(m.moran.tw6$gam)
gam.check(m.moran.tw6$gam)
plot(m.moran.tw6$gam, page = 1)


# final model MOran's I
fin.m.moran <- m.moran.tw6$gam


ggplot(avg_data_moran_sub, aes(x = sum_ips, #peak_diff,
                           y =  Morans_I_log)) +
  geom_point() +
  geom_smooth()











##### quick plotting -----------------------------------------------
fin.m <- m.counts.tw.test.f.rndm#  fin.m.agg#$gam

# plot in easy way how tdoes the k value affect interaction
#p1 <- ggpredict(fin.m, terms = "year [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = "spei_lag2 [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m, terms = c("tmp_lag1", "spei_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

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

# explore results for individual dependent variables: select the onse with the lowest AIC/with the all random effects
fin.m.counts     <- m.counts.tw$gam #result$sum_ips$models$random_effect$gam  # lowest AIC
fin.m.agg        <- m2.agg.gamma$gam #m1.gamma$gam#result$tr_agg_doy$models$random_effect$gam  # lowest AIC
fin.m.peak       <- m.peak.gamma$gam #m.peak.diff.tw$gam #result$tr_peak_doy$models$random_effect$gam  # lowest AIC
fin.m.peak.diff  <- m.peak.diff.tw$gam #result_peak_diff$peak_diff$models$random_effect$gam  # lowest AIC


save(fin.m.counts, 
     fin.m.agg, 
     fin.m.peak,
     fin.m.peak.diff,
     m.counts.tw,
     m1.gamma,
     m.peak.gamma,
     m.peak.diff.tw,
     avg_data,
     file = "outData/fin_models.RData")









# identify significant traps: ---------------------------------------------
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
create_effect_year <- function(data, avg_data, line_color = "grey40", 
                               x_col = "year", 
                               y_col = "sum_ips", 
                               x_title = "X-axis", 
                               y_title = "Y-axis", my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() + 
    geom_line(data = avg_data, aes(x = !!x_col, y = !!y_col, group = pairID), col = "gray70", alpha = 0.4) +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color) +
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}


# Create effect plot function with additional arguments to select columns from avg_data
create_effect_plot <- function(data, avg_data, 
                               x_col = "tmp_z_lag1", 
                               y_col = "sum_ips", 
                               line_color = "blue", x_title = "X-axis", 
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
plot_effect_interactions <- function(data, temp_label = 'Temp', y_title = 'y_title',x_annotate = 0, lab_annotate = "lab ann") {
  #library(ggplot2)
  
 p<- 
   ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
    geom_line(aes(color = group, linetype = group), linewidth = 1) +
    labs(x = temp_label,
         y = y_title,
         fill = "SPEI\nlevels",
         color = "SPEI\nlevels",
         linetype = "SPEI\nlevels") +  # Fixed "y_title" to "y" for correct y-axis label argument
  
   guides(color = guide_legend(nrow = 1), 
           fill = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1)) +
    annotate("text", x = x_annotate, y = Inf, 
             label = lab_annotate, hjust = 0.5, vjust = 1.5) +
   my_theme_square() +
   theme(legend.position = "bottom")
  
  return(p)
}

plot_effect_interactions(p3)

##### PLOT: IPS SUM counts -----------------------------------------------------------

temp_label <- "Temp. [z-score]" #expression(paste("Temperature [", degree, "C]", sep=""))
spei_label <- 'SPEI [z-score]'



# Define the transformation function
adjust_predictions_counts <- function(df, divisor) {
  df <- as.data.frame(df)
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

# Output the result
transformed_value

transform_predictions_DOY(0.2)
# Assuming 'model' is your glm.nb model
summary(fin.m.counts)
p0 <- ggpredict(fin.m.counts, terms = "year [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.counts, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.counts, terms = c("tmp_z_lag1", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)

divisor_population <- 100
p0 <- adjust_predictions_counts(p0, divisor_population)
p1 <- adjust_predictions_counts(p1, divisor_population)
p2 <- adjust_predictions_counts(p2, divisor_population)
p3 <- adjust_predictions_counts(p3, divisor_population)

# get update table for the line an pint plotting
avg_data_sum_ips <- avg_data %>% 
  mutate(sum_ips = sum_ips/100,
         peak_diff = peak_diff/10)

p0.count <- create_effect_year(data = p0, 
                              avg_data = avg_data_sum_ips,
                              x_col = "year",
                              y_col = "sum_ips",
                              line_color = "darkgreen", 
                              x_title = "Year", 
                              y_title = paste(lab_popul_level, '*100'),# "Population level\n[# beetle*100]",
                              my_title = paste("[a]", 'Population level', '\n[#*100]'), 
                              x_annotate = 2018, lab_annotate = "***")

#(p0.count)
p1.count <- 
  create_effect_plot(data = p1, avg_data = avg_data_sum_ips, 
                     x_col = "tmp_z_lag1", y_col = "sum_ips", 
                     line_color = "red", 
                     x_title = temp_label , 
                     y_title = paste(lab_popul_level, '*100'), 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "n.s.")

#(p1.count)
p2.count <- create_effect_plot(data = p2, avg_data = avg_data_sum_ips, 
                               x_col = "spei_z_lag2", y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = spei_label , 
                               y_title = paste(lab_popul_level, '*100'), 
                               #my_title = "Effect of Year on Sum IPS", 
                               x_annotate = -1, 
                               lab_annotate = "**")
#(p2.count)
p3.count <- plot_effect_interactions(p3, 
                                     temp_label = temp_label, 
                                     y_title = paste(lab_popul_level, '*100'),
                                     x_annotate = 2,
                                     lab_annotate = "***") 

#(p3.count)
#ggarrange(p0.count,p1.count,p2.count, p3.count)




##### PLOT DOY aggregation ---------------------------------------------------------
summary(fin.m.agg)
p0 <- ggpredict(fin.m.agg, terms = "year [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.agg, terms = "tmp_z_lag2 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg, terms = "spei_z_lag1 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.agg, terms = c("tmp_z_lag2", "spei_z_lag1 [-1, 0, 1]"), allow.new.levels = TRUE)

# Apply transformation to each prediction data frame
p0 <- transform_predictions_DOY(p0, doy.start, doy.end)
p1 <- transform_predictions_DOY(p1, doy.start, doy.end)
p2 <- transform_predictions_DOY(p2, doy.start, doy.end)
p3 <- transform_predictions_DOY(p3, doy.start, doy.end)


# 
p0.agg <- create_effect_year(data = p0, 
                               avg_data = avg_data_sum_ips,
                               x_col = "year",
                               y_col = "agg_doy",
                               line_color = "darkgreen", 
                               x_title = "Year", 
                               y_title = lab_colonization_time,
                               my_title = paste("[b]", 'Aggregation timing', '\n[DOY]'),  
                               x_annotate = 2018, lab_annotate = "***")

(p0.agg)
p1.agg <- 
  create_effect_plot(data = p1, avg_data = avg_data_sum_ips, 
                     x_col = "tmp_z_lag1", y_col = "agg_doy", 
                     line_color = "red", 
                     x_title = temp_label , 
                     y_title = lab_colonization_time, 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "***")

(p1.agg)
p2.agg <- create_effect_plot(data = p2, avg_data = avg_data_sum_ips, 
                               x_col = "spei_z_lag2", y_col = "agg_doy", 
                               line_color = "blue", 
                               x_title = spei_label , 
                               y_title = lab_colonization_time, 
                               x_annotate = -1, 
                               lab_annotate = "n.s.")
(p2.agg)
p3.agg <- plot_effect_interactions(p3, 
                                     temp_label = temp_label, 
                                     y_title = lab_colonization_time,
                                     x_annotate = 2,
                                     lab_annotate = "***") 

(p3.agg)


#### PLOT: Peak DOY  ------------------------------------------------------------
summary(fin.m.peak)

p0 <- ggpredict(fin.m.peak, terms = "year [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak, terms = c("tmp_z_lag1" ,"spei_z_lag2 [-1, 0, 1]"), allow.new.levels = TRUE)


# transform values back to DOY (from 0-1)
p0 <- transform_predictions_DOY(p0, doy.start, doy.end)
p1 <- transform_predictions_DOY(p1, doy.start, doy.end)
p2 <- transform_predictions_DOY(p2, doy.start, doy.end)
p3 <- transform_predictions_DOY(p3, doy.start, doy.end)

# 
p0.peak <- create_effect_year(data = p0, 
                             avg_data = avg_data_sum_ips,
                             x_col = "year",
                             y_col = "peak_doy",
                             line_color = "darkgreen", 
                             x_title = "Year", 
                             y_title = lab_peak_time,
                             my_title = paste("[c]", 'Peak swarming\ntiming', '[DOY]'),  
                             x_annotate = 2018, lab_annotate = "***")

(p0.peak)
p1.peak <- 
  create_effect_plot(data = p1, avg_data = avg_data_sum_ips, 
                     x_col = "tmp_z_lag1", y_col = "peak_doy", 
                     line_color = "red", 
                     x_title = temp_label , 
                     y_title = lab_peak_time, 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = "n.s.")

(p1.peak)
p2.peak <- create_effect_plot(data = p2, avg_data = avg_data_sum_ips, 
                             x_col = "spei_z_lag2", y_col = "peak_doy", 
                             line_color = "blue", 
                             x_title = spei_label , 
                             y_title = lab_peak_time, 
                             x_annotate = -1, 
                             lab_annotate = "***")
(p2.peak)
p3.peak <- plot_effect_interactions(p3, 
                                   temp_label = temp_label, 
                                   y_title = lab_peak_time,
                                   x_annotate = 2,
                                   lab_annotate = "***") 

(p3.peak)



##### PLOT: Peak diff  ----------------------------------------------------
fin.m.peak.diff #<- fin.m.peak.diff$gam
summary(fin.m.peak.diff)
p0 <- ggpredict(fin.m.peak.diff, terms = "year [all]", allow.new.levels = TRUE)
p1 <- ggpredict(fin.m.peak.diff, terms = "tmp_z_lag1 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff, terms = "spei_z_lag2 [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.peak.diff, terms = c("tmp_z_lag1", "spei_z_lag2 [-1, 0, 1]"), allow.new.levels = T)


divisor_diff <- 10
p0 <- adjust_predictions_counts(p0, divisor_diff)
p1 <- adjust_predictions_counts(p1, divisor_diff)
p2 <- adjust_predictions_counts(p2, divisor_diff)
p3 <- adjust_predictions_counts(p3, divisor_diff)





p0.peak.diff <- create_effect_year(data = p0, 
                               avg_data = dplyr::filter(avg_data_sum_ips, peak_diff <200),
                               x_col = "year",
                               y_col = "peak_diff",
                               line_color = "darkgreen", 
                               x_title = "Year", 
                               y_title = lab_peak_growth,
                               my_title = paste("[d]", 'Peak swarming\nintensity', '[#*10]'),  
                               x_annotate = 2018, lab_annotate = "***")

#(p0.peak.diff)
p1.peak.diff <- 
  create_effect_plot(data = p1, 
                     avg_data =  dplyr::filter(avg_data_sum_ips, peak_diff <200), 
                     x_col = "tmp_z_lag1", y_col = "peak_diff", 
                     line_color = "red", 
                     x_title = temp_label , 
                     y_title = lab_peak_growth, 
                     #my_title = "Effect of Year on Sum IPS", 
                     x_annotate = 2, 
                     lab_annotate = ".")

#(p1.peak.diff)
p2.peak.diff <- create_effect_plot(data = p2,
                                   avg_data =  dplyr::filter(avg_data_sum_ips, peak_diff <200), 
                               x_col = "spei_z_lag2", y_col = "peak_diff", 
                               line_color = "blue", 
                               x_title = spei_label , 
                               y_title = lab_peak_growth, 
                               #my_title = "Effect of Year on Sum IPS", 
                               x_annotate = -1, 
                               lab_annotate = "n.s.")
#(p2.peak.diff)
p3.peak.diff <- plot_effect_interactions(p3, 
                                     temp_label = temp_label, 
                                     y_title = lab_peak_growth,
                                     x_annotate = 2,
                                     lab_annotate = "***") 

#(p3.peak.diff)
ggarrange(p0.peak.diff,p1.peak.diff,p2.peak.diff,p3.peak.diff)




# put in teh same figure? --------------------------------------------------------
# remove labels, keep only labels on top


#lab_popul_level <- 
#lab_colonization_time
#lab_peak_time
#lab_peak_growth

p.full.year <- ggarrange(p0.count, p0.agg, p0.peak, p0.peak.diff, 
                         ncol=4, nrow = 1 , #align = 'hv', 
                         font.label = list(size = 8, color = "black", face = "plain", family = NULL) )

(p.full.year)

p.temp <-  ggarrange(p1.count, p1.agg, p1.peak, p1.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv', 
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.spei <-  ggarrange(p2.count, p2.agg, p2.peak, p2.peak.diff, 
                     ncol=4, nrow = 1 , align = 'hv', 
                     font.label = list(size = 8, color = "black", face = "plain", family = NULL))


p.int <- ggarrange(p3.count, p3.agg, p3.peak, p3.peak.diff, 
                   ncol=4, nrow = 1 , align = 'hv',common.legend = TRUE, legend = 'bottom',
                   font.label = list(size = 8, color = "black", face = "plain", family = NULL))
windows(7,8)
full_preds <- ggarrange(p.full.year,p.temp, p.spei, p.int, 
                        ncol = 1, nrow= 4, 
                        align = 'hv', heights = c(1, 1, 1, 1.1),
                        widths = c(1, 1, 1, 1))

(full_preds)
ggsave(filename = 'outFigs/Fig_full_eff.png', plot = full_preds, 
       width = 7.5, height = 8, dpi = 300, bg = 'white')
















### all Effect plots: ips vs climate --------------------------------------------------------



# put together individual plots for climate
p.count     <- ggarrange(p1.count, p2.count, ncol = 1)
p.agg       <- ggarrange(p1.agg, p2.agg, ncol = 1)
p.peak      <- ggarrange(p1.peak, p2.peak, ncol = 1)
p.peak.diff <- ggarrange(p1.peak.diff, p2.peak.diff, ncol = 1)













# PLOT: Temporal development ----------------------------------
p.full.year <- ggarrange(p0.count, p0.agg, p0.peak, p0.peak.diff, 
                        ncol=4, nrow = 1 , align = 'hv', 
                        font.label = list(size = 8, color = "black", face = "plain", family = NULL),
                        labels = c( paste("[a] ", lab_popul_level),
                                    paste("[b] ", lab_colonization_time ),
                                    paste("[c] ", lab_peak_time ),
                                    paste("[d] ", lab_peak_growth )))


windows(7,5)
(p.full.year)


p.out.clim <- ggarrange(p.count, p.agg, p.peak, p.peak.diff, 
          ncol=4, nrow = 1 , align = 'hv', 
          font.label = list(size = 8, color = "black", face = "plain", family = NULL),
          labels = c( paste("[a] ", lab_popul_level),
                      paste("[b] ", lab_colonization_time ),
                      paste("[c] ", lab_peak_time ),
                      paste("[d] ", lab_peak_growth )))
                      

windows(7,4)
(p.out.clim)

ggsave(filename = 'outFigs/Fig3.png', plot = p.out.clim, 
       width = 7, height = 4, dpi = 300, bg = 'white')


p.clim.int <- ggarrange(p3.count, p3.agg, p3.peak, p3.peak.diff, 
          ncol=4, nrow = 1 , align = 'hv',common.legend = TRUE, legend = 'right',
          font.label = list(size = 8, color = "black", face = "plain", family = NULL),
          labels = c( paste("[a] ", lab_popul_level),
                      paste("[b] ", lab_colonization_time ),
                      paste("[c] ", lab_peak_time ),
                      paste("[d] ", lab_peak_growth )))
windows(7,2)
p.clim.int
ggsave(filename = 'outFigs/Fig4.png', plot = p.clim.int, 
       width = 7, height = 4, dpi = 300, bg = 'white')









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

# PLOT Moran's I -------------------------------
summary(fin.m.moran)

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.moran , terms = "tmp_z_lag1  [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.moran, terms = "spei_z_lag2  [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m.moran, terms = "sum_ips [all]", allow.new.levels = TRUE)
p4 <- ggpredict(fin.m.moran, terms = c("tmp_z_lag1", "spei_z_lag2 [-1,0,1]"), allow.new.levels = TRUE) 


p1.moran <- create_effect_plot(p1,
                               avg_data =  avg_data_moran_sub, 
                               x_col = "tmp_z_lag1", 
                               y_col = "Morans_I_log", 
                               line_color = "red", 
                               x_title = "Temp. [z-score]", 
                               y_title = "Local Moran's I", 
                               x_annotate = 2.5,
                               lab_annotate = "**"
                               #y_lim = c(0,1)
                               )

p1.moran
p2.moran <- create_effect_plot(p2,
                               avg_data =  avg_data_moran_sub, 
                               x_col = "spei_z_lag2", 
                               y_col = "Morans_I_log", 
                               line_color = "blue", 
                               x_title = "SPEI [z-score]", 
                               y_title = "Local Moran's I",  
                               #y_lim = c(0,1)
                               x_annotate = -1.5,
                               lab_annotate = "**"
                               )
p3.moran <- create_effect_plot(p3, 
                               avg_data =  avg_data_moran_sub, 
                               x_col = "sum_ips", 
                               y_col = "Morans_I_log", 
                               line_color = "grey50", 
                               x_title = "Population level [log(#)]", 
                               y_title = "Local Moran's I", # y_lim = c(0,1),
                               x_annotate = 10,
                               lab_annotate = "***")
p4.moran <- plot_effect_interactions(p4,
                                    # avg_data =  avg_data_moran_sub, 
                                    # x_col = "tmp_z_lag1", 
                                     #y_col = "Morans_I_log", 
                                     temp_label = "Temperature [dim.]", 
                                     y_title = "Local Moran's I",
                                     x_annotate = 2.5,
                                     lab_annotate = "*")# +
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



summary(fin.m.counts )
summary(fin.m.agg)
summary(fin.m.peak)
summary(fin.m.peak.diff)
#summary(fin.m.moran)
#summary(fin.m.RS)



# print models outputs: ---------------------------------------------------------



# export all models
sjPlot::tab_model(fin.m.counts,    file = "outTable/model_counts.doc")
sjPlot::tab_model(fin.m.agg,       file = "outTable/model_agg.doc")
sjPlot::tab_model(fin.m.peak,      file = "outTable/model_peak.doc")
sjPlot::tab_model(fin.m.peak.diff, file = "outTable/model_peak_diff.doc")
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




