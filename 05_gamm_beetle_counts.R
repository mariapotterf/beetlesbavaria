



# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(ggplot2)
library(ggpubr)
library(ggpmisc)  # add equation to plots smooth 

# Stats
library('here')
library('mgcv')
library('gratia')
library('gamair')
library('purrr')
library('mvnfast')
library("tibble")
library('gganimate')
library('cowplot')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('itsadug')

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





# read data ----------------------------------------------------------------

load(file="outData/final.Rdata")
load(file =  "outData/buffers.Rdata")  # df_RS_out

# get unique falsto_names of traps that remained stable (only one globalid)
# over time
stable_traps <- 
  df_RS_out %>% 
  group_by(falsto_name, globalid) %>% 
  dplyr::summarize(n = n()) %>% 
  dplyr::filter(n == 7) %>% 
  distinct(falsto_name) %>% 
  pull() %>% 
  as.vector()





# clean up data from remoe sensing --------------------------------------------
df_RS_out <- df_RS_out %>% 
  dplyr::select(-c(globalid, id)) %>% 
  dplyr::rename(trap_name = falsto_name)

# add lag: previous year counts, previous year temperature
dat_lag <-   dat_fin %>%
    ungroup(.) %>% 
   group_by(trap_name) %>%
    dplyr::mutate(previous_sum_ips =  dplyr::lag(sum_ips, order_by = year),
                 # previous_sum_ips2 =  dplyr::lag(sum_ips, order_by = year, n = 2),
                  previous_tmp     =  dplyr::lag(tmp, order_by = year),
                  previous_spei     =  dplyr::lag(spei, order_by = year)) %>% 
  dplyr::mutate(population_growth = (sum_ips - previous_sum_ips) / previous_sum_ips * 100,
                population_growth2 = dplyr::lag(population_growth, order_by = year)) %>%  # lag population growth by one more year
  left_join(df_RS_out, by = join_by(trap_name, year, sum_ips)) %>% 
  dplyr::mutate(trap_name = as.factor(trap_name)) %>% 
  dplyr::mutate(next1_RS_beetle  = dplyr::lag(wind_beetle, n = 1, order_by = year),
                next1_RS_harvest = dplyr::lag(harvest,     n = 1, order_by = year),
                next2_RS_beetle  = dplyr::lag(wind_beetle, n = 2, order_by = year),
                next2_RS_harvest = dplyr::lag(harvest,     n = 2, order_by = year)) %>%
  dplyr::mutate(temp_fact = case_when(year== 2018 ~ 'extreme',
                                      year %in% c(2015,2019, 2020) ~ 'high',
                                      year %in% c(2016,2017) ~ 'moderate',
                                      year== 2021 ~ 'low'),
                temp_fact = as.factor(temp_fact)) 
  
  
    #dplyr::select(c(year, trap_name, 
    #                previous_tmp, tmp,
    #                next2_RS_harvest, harvest )) #%>% 
  
# check if results are correct
# dat_lag %>%
#   filter(trap_name== 'Zusmarshausen_1') %>% 
#   dplyr::select(c(trap_name, year, sum_ips, previous_sum_ips,population_growth, population_growth2)) %>% 
#   View()



# split df into individual files: predist ips_sum, DOY agg, doy _peak, 
# predict RS damage

# [1] "trap_name"         "year"              "conif_prop"        "elev"              "x"                
# [6] "y"                 "sm"                "vpd"               "sm_z"              "vpd_z"            
# [11] "tmp"               "tmp_z"             "spei"              "spei_z"            "sum_ips"          
# [16] "peak_doy"          "peak_diff"         "agg_doy"           "location"          "previous_sum_ips" 
# [21] "previous_tmp"      "previous_spei"     "population_growth" "wind_beetle"       "harvest"          
# [26] "spruce_freq"       "all_dist_sum"      "cum_removed"       "remained_forest"   "next1_RS_beetle"  
# [31] "next1_RS_harvest"  "next2_RS_beetle"   "next2_RS_harvest" 

# predict aggregation date: ----------------------------------------------------
df_fin_agg <- dat_lag %>% 
  dplyr::select(-c("wind_beetle",       "harvest",          
                     "all_dist_sum",      "cum_removed",       "remained_forest" ,  "next1_RS_beetle" , 
                   "next1_RS_harvest",  "next2_RS_beetle",   "next2_RS_harvest" ))

df_fin_peak <- dat_lag %>% 
  dplyr::select(-c("wind_beetle",       "harvest",          
                   "all_dist_sum",      "cum_removed",       "remained_forest" ,  "next1_RS_beetle" , 
                   "next1_RS_harvest",  "next2_RS_beetle",   "next2_RS_harvest" ))

# df_RS_fin <- dat_lag %>% 
#   dplyr::select(c("sum_ips", 
#                   "trap_name",         "year" ,             "conif_prop" ,       "elev"  ,            "x",                
#                    "y" ,
#                    "location",
#                   "conif_prop", "elev", "spruce_freq",
#                   "peak_doy", 
#                   "peak_diff" , 
#                   "agg_doy",  "wind_beetle",  "harvest",          
#                    "all_dist_sum",      "cum_removed",       "remained_forest" ,  "next1_RS_beetle" , 
#                    "next1_RS_harvest",  "next2_RS_beetle",   "next2_RS_harvest" ))

  

# Scale predictors ================================================================================

dat_lag_scaled                <- dat_lag
dat_lag_scaled$tmp            <- scale(dat_lag_scaled$tmp)
dat_lag_scaled$vpd            <- scale(dat_lag_scaled$vpd)
dat_lag_scaled$spei           <- scale(dat_lag_scaled$spei)
dat_lag_scaled$spruce_freq    <- scale(dat_lag_scaled$spruce_freq)
# code 'year' as numeric to account for temp autocorrelation
dat_lag_scaled$year_num       <- dat_lag_scaled$year
dat_lag_scaled$year           <- factor(dat_lag_scaled$year )
dat_lag_scaled$sum_ips_scaled <- scale(dat_lag_scaled$sum_ips)
dat_lag_scaled$previous_sum_ips_scal <- scale(dat_lag_scaled$previous_sum_ips)
dat_lag_scaled$id             <- 1:nrow(dat_lag_scaled)  # Adding an observation ID

#simplify, keep all columns as vectors

dat_lag_scaled_simpl <- dat_lag_scaled %>% 
  ungroup(.) %>% 
  dplyr::mutate(across(where(~is.numeric(.) && is.matrix(.)), as.vector))

# consider year as random factor - non linear relationship betweeen years

# how to improve model:
# capp values:

cap_at_percentiles <- function(x, lower = 0.01, upper = 0.99) {
  bounds <- quantile(x, c(lower, upper))
  pmin(pmax(x, bounds[1]), bounds[2])
}

# Apply to your data
dat_lag_scaled$capped_sum_ips <- cap_at_percentiles(dat_lag_scaled$sum_ips)


# Exploratory analyses IPS_SUM --------------------------------------------------
windows()

pairs(sum_ips ~ tmp + spei + vpd, dat_lag, panel = panel.smooth)

pairs(sum_ips ~ tmp + spei + vpd + elev + spruce_freq + previous_sum_ips+population_growth, dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+population_growth + tmp + previous_tmp , dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+ population_growth + tmp + previous_tmp + previous_spei , dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+ population_growth + tmp + previous_tmp + population_growth2, dat_lag, panel = panel.smooth)

hist(dat_lag$population_growth)

# RS:  Link beetle data with observed RS damage ------------------------------------------------------------------------
# remove 0 wind-beetle damage 
# run analyses with only stable traps over time?? - split df


## [1] Only stable traps : Beetle damage  ---------------------------------------
dat_lag_stable <- dat_lag %>% 
  filter(trap_name %in% stable_traps) #

dat_lag_stable %>% 
  filter(wind_beetle > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = wind_beetle)) +
  geom_point() +
  geom_smooth() + ggtitle('Current year')

dat_lag_stable %>% 
  filter(next1_RS_beetle > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = next1_RS_beetle)) +
  geom_point() +
  geom_smooth() + ggtitle('Next year')

dat_lag_stable %>% 
  filter(next2_RS_beetle > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = next2_RS_beetle)) +
  geom_point() +
  geom_smooth() + ggtitle('Next 2 years')


## [2] Only stable traps : Harvest  ------------------------------------------
p_harv1 <- dat_lag_stable %>% 
  filter(harvest > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = harvest)) +
  geom_point() +
  geom_smooth() + ggtitle('Current year')

p_harv2 <- dat_lag_stable %>% 
  filter(next1_RS_harvest > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = next1_RS_harvest)) +
  geom_point() +
  geom_smooth() + ggtitle('Next year')

p_harv3 <- dat_lag_stable %>% 
  filter(next2_RS_harvest > 0) %>% 
  ggplot(aes(x = sum_ips, 
             y = next2_RS_harvest)) +
  geom_point() +
  geom_smooth() + ggtitle('Next 2 years')

ggarrange(p_harv1, p_harv2, p_harv3)

# !!!

pairs(next1_RS_beetle ~  previous_sum_ips+ sum_ips+ population_growth + tmp + previous_tmp + previous_spei , dat_lag_stable, panel = panel.smooth)

pairs(next1_RS_beetle ~ previous_sum_ips+ sum_ips +population_growth + population_growth2, #+
      #  tmp + previous_tmp + previous_spei , 
      dat_lag_stable, panel = panel.smooth)


# Predict beetle sums/year by environmantel predictors --------------------
#

# make glm - beetle sums by temp and spei
# test with the raw and scaled predictors;
# glm, - define the polynomial relationship

p.tmp <- dat_lag_scaled %>% 
ggplot(aes(x = tmp,
           y = sum_ips)) +
  geom_point()  +
  geom_smooth()

p.sm <- dat_lag_scaled %>% 
  ggplot(aes(x = sm,
             y = sum_ips)) +
  geom_point()  +
  geom_smooth()


p.vpd <-dat_lag_scaled %>% 
  ggplot(aes(x = vpd,
             y = sum_ips)) +
  geom_point()  +
  geom_smooth()


ggarrange(p.sm, p.tmp, p.vpd, nrow = 1)

## test GLM (no random effects) on raw data ---------------------------------------------------------------------
m1.poly <- glm(sum_ips ~ I(tmp^15), dat_lag, family = 'poisson')
m1.exp <- glm(sum_ips ~ log(tmp), dat_lag, family = 'poisson')



# see in plot
# New data for prediction
newdata_poly <- data.frame(tmp = seq(min(dat_lag$tmp), max(dat_lag$tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.poly, newdata = newdata_poly, type = "response")


ggplot(dat_lag, aes(x = tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_poly, aes(x = tmp, y = predicted), color = "red",  lwd = 1.5) +
  geom_line(data = newdata_exp, aes(x = tmp, y = predicted), color = "green", lwd = 1.5) +
  geom_smooth()



# New data for prediction: exponential
newdata_exp <- data.frame(tmp = seq(min(dat_lag$tmp), max(dat_lag$tmp), length.out = 100))

# Generate predictions
newdata_exp$predicted <- predict(m1.exp, newdata = newdata_exp, type = "response")
ggplot(dat_lag, aes(x = tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_exp, aes(x = tmp, y = predicted), color = "red")


## Fit the GLMM with negative binomial distribution & random effects ------------
m1.exp <- glm.nb(sum_ips ~ tmp, data = dat_lag, link = log)
newdata_exp <- data.frame(tmp = seq(min(dat_lag$tmp), max(dat_lag$tmp), length.out = 100))
newdata_exp$predicted <- predict(m1.exp, newdata = newdata_exp, type = "response")

ggplot(dat_lag, aes(x = tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_exp, aes(x = tmp, y = predicted), color = "red")


m1.poly2 <- glm.nb(sum_ips ~ poly(tmp, 2), data = dat_lag)
newdata_poly <- data.frame(tmp = seq(min(dat_lag$tmp), max(dat_lag$tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.poly, newdata = newdata_poly, type = "response")


m1.poly3 <- glm.nb(sum_ips ~ poly(tmp, 3), data = dat_lag)
m1.sc.poly2 <- glm.nb(adjusted_variable  ~ poly(tmp, 2), data = dat_lag_scaled)
newdata_poly <- data.frame(tmp = seq(min(dat_lag_scaled$tmp), max(dat_lag_scaled$tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.sc.poly2, newdata = newdata_poly, type = "response")


ggplot(dat_lag_scaled, aes(x = tmp, y = adjusted_variable)) +
  geom_point() +
  geom_line(data = newdata_poly, aes(x = tmp, y = predicted), color = "blue")


AIC(m1, m1.exp, m1.poly2, m1.poly3)
simulationOutput <- simulateResiduals(fittedModel = m1.sc.poly2, plot = T)
testOutliers(m1.sc.poly2, type = 'bootstrap') 


m2 <- glm.nb(sum_ips ~ tmp + spei, dat_lag)
m3 <- glm.nb(sum_ips ~ tmp + tmp_z, dat_lag)
m4 <- glm.nb(sum_ips ~ I(tmp) + previous_sum_ips   , dat_lag)
summary(m4)
plot(m4, page = 1)




# DREDGE Variable selection ------------------------------------------------------

# Fit a global model with all potential predictors
global_model <- glm.nb(sum_ips ~ tmp + spei + sm + vpd + previous_sum_ips, 
                      data = dat_lag_scaled, na.action = "na.fail")

# Run dredge on the global model
model_set <- dredge(global_model)

# View the results
print(model_set)

# Sort by AIC
model_set_sorted <- model_set[order(model_set$AIC),]

# View the top models with the lowest AIC
head(model_set_sorted)


### GLMM - with random effect and nb ----------------------------------------


m3 <- glmer.nb(sum_ips ~ tmp + spei + (1|location) + (1|trap_name), data = dat_lag)

m3.1 <- glmer.nb(sum_ips ~ tmp  + (1|trap_name), data = dat_lag)

m3.2 <- glmer.nb(sum_ips ~ tmp  + (1|location), data = dat_lag)
m3.3 <- glmer.nb(sum_ips ~ tmp + (1| year), data = dat_lag)

m3.4 <- glmer.nb(sum_ips ~ tmp + year + (1| location), data = dat_lag)

m3.5 <- glmer.nb(sum_ips ~ tmp + (1| year), data = dat_lag_scaled)
# consider as polynomial
m3.5.1 <- glmer.nb(sum_ips ~ poly(tmp,2) + (1| year), data = dat_lag_scaled)

# add locatin as random effect
m3.5.2 <- glmer.nb(sum_ips ~ poly(tmp,2) + (1| year) + (1| location), data = dat_lag_scaled)


# simplify year effect" replace by drought ref
m3.5.3 <- glmer.nb(sum_ips ~ poly(tmp,2) + (1| temp_fact) + (1| location), data = dat_lag_scaled)


summary(m3.5.2)
AICc(m3.5.2, m3.5.3, m3.7)  # m3.7 is better

m3.6 <- glmer.nb(sum_ips ~ tmp + (1| year) + (1| location), data = dat_lag_scaled)

m3.7 <- glmer.nb(sum_ips ~ tmp + (1| year) + (1| location/trap_name), data = dat_lag_scaled)

# add poly function, year is consider as random to account for diversity over years
m3.7.1 <- glmer.nb(sum_ips ~ poly(tmp,2) + (1| year) + (1| location/trap_name), data = dat_lag_scaled)


# account for temp autocorrelation: use lagged values to capture it
m3.7.1.aut <- glmer(sum_ips ~ poly(tmp, 2) + previous_sum_ips + (1 | year) + (1 | location/trap_name), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))


# scale previous sum ips
m3.7.1.aut.sc1 <- glmer(sum_ips ~ poly(tmp, 2) + previous_sum_ips_scal + (1 | year) + (1 | location/trap_name), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))

# remove 'year' as random
m3.7.1.aut.sc2 <- glmer(sum_ips ~ poly(tmp, 2) + previous_sum_ips_scal+ (1 | location/trap_name), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))


# scale the lagged predictor; previous beetle sums
m3.7.1.aut2 <- glmer(sum_ips ~ tmp + previous_sum_ips_scal + (1 | year) + (1 | location/trap_name), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))

# year as numeric, random
m3.7.1.aut2.1 <- glmer(sum_ips ~ tmp + previous_sum_ips_scal + (1 | year_num) + (1 | location/trap_name), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))

# remove 'year' as random effect
m3.7.1.aut3 <- glmer(sum_ips ~ tmp + previous_sum_ips_scal +  (1 | location/trap_name), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))

# add year as nested random effect
m3.7.1.aut4 <- glmer(sum_ips ~ tmp + previous_sum_ips_scal +  (1 | location/trap_name/year), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))


# add poly function, year is consider as random to account for diversity over years
m3.7.1.quadr <- glmer.nb(sum_ips ~ I(tmp^2) + (1| year) + (1| location/trap_name), data = dat_lag_scaled)

AICc(m3.7.1.aut, #t5his one is the best! (11/07/2023)
     m3.7.1.aut.sc1,
     m3.7.1.aut.sc2,
     m3.7.1.aut1,
     m3.7.1.aut3,
     m3.7.1.aut2.1,
     m3.7.1.aut2, # 
     m3.7.1.aut4) 

BIC(m3.7.1.aut, #t5his one is the best! (11/07/2023)
     m3.7.1.aut.sc1,
     m3.7.1.aut.sc2,
     m3.7.1.aut1,
     m3.7.1.aut3,
     m3.7.1.aut2.1,
     m3.7.1.aut2, # # the best!!!
     m3.7.1.aut4) 

acf(residuals(m3.7.1.aut.sc1))  # seems relatively ok
pacf(residuals(m3.7.1.aut.sc1))

plot(allEffects(m3.7.1.aut.sc1))
dw_result <- dwtest(residuals(m3.7.1.aut.sc1) ~ fitted(m3.7.1.aut.sc1))
(dw_result)

partialR2(m3.7.1.aut2)
# Durbin-Watson test
# 
# data:  residuals(m3.7.1.aut2) ~ fitted(m3.7.1.aut2)
# DW = 1.6613, p-value = 8.204e-08
# alternative hypothesis: true autocorrelation is greater than 0

simulationOutput <- simulateResiduals(fittedModel = m3.7.1.aut2, plot = T)
simulationOutput <- simulateResiduals(fittedModel = m3.7.1.aut, plot = T) # this one is little bit better from output

# goodness of fit!
r2(m3.7.1.aut2) # Nagelkerke s R2

# account for autocorrelation: use GAM and MER model - does not work!!
library(gamm4)
gamm4_model <- gamm4(sum_ips ~ poly(tmp, 2), 
                     random = ~(1 | year_num) + (1 | location/trap_name), 
                     data = dat_lag_scaled, 
                     #family = nb(2.8588),
                     correlation = corAR1(form = ~ 1 | year_num))


# simplify years by heat cetegories
m3.7.2 <- glmer.nb(sum_ips ~ poly(tmp,2) + (1| temp_fact) + (1| location/trap_name), data = dat_lag_scaled)

# add interaction heat vs ips_sums
m3.7.3.poly <- glmer.nb(sum_ips ~ poly(tmp,2)*temp_fact + (1| location/trap_name), data = dat_lag_scaled)

# keep it linear
m3.7.4.lin <- glmer.nb(sum_ips ~ tmp*temp_fact + (1| location/trap_name), data = dat_lag_scaled)

# keep it exp
m3.7.4.exp <- glmer.nb(sum_ips ~ exp(tmp)*temp_fact + (1| location/trap_name), data = dat_lag_scaled)


AICc(m3.7.1, m3.7.3.poly, m3.7.4.lin, m3.7.4.exp)

# try as smooths??
m3.7.3.smooth <- glmer.nb(sum_ips ~ splines::bs(tmp, df = 4)*temp_fact + (1| location/trap_name), data = dat_lag_scaled)



#splines::bs(x, degrees = 3) * factor(category)

# groups years by factor: temp factor by years???

AICc(m3.7, m3.7.1, m3.7.2) # the best: m3.7.1

m3.8 <- glmer.nb(sum_ips ~ tmp  + (1| year) + (1| location/trap_name), data = dat_lag_scaled)

m3.5.1 <- glmer.nb(sum_ips ~ tmp  + (1| year), data = dat_lag_scaled)


summary(m3.7.2)
AICc(m3.5, m3.6, m3.7)

r2(m3.7.2)

simulationOutput <- simulateResiduals(fittedModel = m3.7.3.poly, plot = T)


# Predicting beetle counts
AIC(m3.5, m3.7.1, m3.7.3.poly)

### plot GLM.nb  ===================================================
unique_trap_names <- levels(dat_lag_scaled$trap_name)
unique_location   <- levels(dat_lag_scaled$location)
unique_temp       <- levels(dat_lag_scaled$temp_fact)
temp_range        <- seq(min(dat_lag_scaled$tmp), 
                         max(dat_lag_scaled$tmp), 
                         length.out = 100)
new_data <- expand.grid(tmp = temp_range, 
                        temp_fact = levels(dat_lag_scaled$temp_fact),
                        location = unique_location[1],
                        trap_name = unique_trap_names[1])

new_data$predicted_counts <- predict(m3.7.4.exp, newdata = new_data, re.form = NULL)

# Plotting
ggplot(new_data, aes(x = tmp, y = predicted_counts, color = temp_fact)) +
  geom_line() +
  labs(title = "Predicted Beetle Counts by Temperature and Year",
       x = "Temperature (z-score)",
       y = "Predicted Beetle Counts") +
  theme_minimal()

#
r2(m3.7.4.exp)
r2(m3.7.3.poly)
r2(m3.7.4.lin)



# ceck different seasons??
ggplot(dat_lag_scaled, aes(x = tmp, y = sum_ips, color = temp_fact)) +
  geom_point() +
  geom_smooth()+
  theme_minimal()



m4 <- glmer.nb(sum_ips ~ tmp + spei + (1|location) + (1|trap_name), 
               data = dat_lag_scaled)

m5 <- glmer.nb(sum_ips ~ tmp + spei + (1 | location/trap_name), # account for nested design
               data = dat_lag_scaled)

m6 <- glmer.nb(sum_ips ~ tmp + (1 | location/trap_name), # account for nested design
               data = dat_lag_scaled)

m7 <- glmer.nb(sum_ips ~ tmp + elev + (1 | location/trap_name), # account for nested design
               data = dat_lag_scaled)  # error

m8 <- glmer.nb(sum_ips ~ tmp + spruce_freq + (1 | location/trap_name), # account for nested design
               data = dat_lag_scaled)  

# account for nteraction between temp and year account for spatial and temporal autocorrelation
m9 <- glmer.nb(sum_ips ~ tmp * year + spruce_freq + (1 | location/trap_name), 
               data = dat_lag_scaled)

# visualize teh effect of temperature over years
# Assuming 'tmp' ranges from min_tmp to max_tmp in your data
tmp_range <- seq(min(dat_lag_scaled$tmp), 
                 max(dat_lag_scaled$tmp), 
                 length.out = 100)

# Create a data frame for predictions
new_data <- expand.grid(tmp = tmp_range, 
                        year = levels(dat_lag_scaled$year), 
                        spruce_freq = mean(dat_lag_scaled$spruce_freq))  # Adjust other variables as needed

# Predict using the model
new_data$predicted <- predict(m9, newdata = new_data, re.form = NA, type = "response")

# Plotting
ggplot(new_data, 
       aes(x = tmp, y = predicted, color = year)) +
  geom_line() +
  labs(title = "Effect of Temperature across Years",
       x = "Temperature",
       y = "Predicted sum_ips") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")



# allow temparature interacting with year, but not year as individual preddictor
m10 <- glmer.nb(sum_ips ~ tmp:year + tmp + spruce_freq + (1 | location/trap_name), 
                data = dat_lag_scaled)



### glmmTMB account for temporal correlation and overdispersion in data -------------


# not meaningfull!!!
m11 <- glmmTMB(sum_ips ~ tmp + spruce_freq + 
                (1 | location/trap_name) + 
                (0 + year | location/trap_name) + 
                 ar1(year + 0|trap_name), # account for autocorrelation, as here: https://stats.stackexchange.com/questions/393734/covariance-structures-in-glmmtmb-for-temporal-autocorrelation
              data = dat_lag_scaled,
              family = nbinom2,
              ziformula = ~ 1,  # Zero-inflation formula, adjust as necessary
              dispformula = ~ 1#,  # Dispersion formula, adjust as necessary
              #corStruct = corAR1(form = ~ year | location/trap_name)
              )

# m11 does not work: as the glmmTMB does not work with the temporal autocorrelation
# create the year as random intercept to account for systematic error in year:

m12 <- glmmTMB(sum_ips ~ tmp + spruce_freq + (1 + year | location/trap_name),
              data = dat_lag_scaled,
              family = nbinom2(),
              ziformula = ~ 1,  # Zero-inflation part, adjust if needed
              dispformula = ~ 1)  # Dispersion part, adjust if needed


AICc( m8, m9, m10)
r2(m3.5)       # get R squared
confint(m5)  # get confidence interval


plot(allEffects(m3.5))

AICc(m5)
dredge(m5)



### GLM check residuals, outlieras and different tests --------------------------- 
windows()
simulationOutput <- simulateResiduals(fittedModel = m3.5, plot = T)
testOutliers(m3.5, type = 'bootstrap') # p_val is n.s.: outliers do not occurs more often then by chance
testDispersion(m3.5)  # P_val is significant, data are overdispersed
testZeroInflation(m3.5)


plotResiduals(simulationOutput)

summary(m5)
plot(m2)

AICc(m2)


r2(m2)



### GAM tested models: ips_sum --------------------------------------------------

dat_lag %>% 
  ggplot(aes(y = sum_ips,
             x = tmp,
             color = factor(year))) +
  geom_smooth()


dat_lag %>% 
  ggplot(aes(y = sum_ips,
             x = tmp,
             color = factor(temp_fact))) +
  geom_point(alpha = 0.3)+
  geom_smooth() +
  facet_grid(.~temp_fact) +
  theme_bw()

# decide about which years are warm and which not??
# check the z-score
dat_lag %>% 
  ggplot(aes(y = tmp_z,
             x = year,
             group = year)) +
  geom_boxplot(notch = T)


m.gam0 <- bam(sum_ips ~ 
                s(tmp, k = 4),# +
              #s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
              #s(location, bs = 're') +   # trap pair
              #s(x,y, k = 50)
              # ,
              dat_lag_scaled, 
              family = tw)


m.gam1 <- bam(sum_ips ~ 
                s(tmp, by=temp_fact, k = 4),# +
            #s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            #s(location, bs = 're') +   # trap pair
            #s(x,y, k = 50)
           # ,
          dat_lag_scaled, 
          family = tw)


m.gam2 <- bam(sum_ips ~ 
                s(tmp, by=temp_fact, k = 4) +
              s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
              s(location, bs = 're') +   # trap pair
              s(x,y),
              dat_lag_scaled, 
              family = tw)

# add base 'ds' for coordinates
m.gam3 <- bam(sum_ips ~ 
                s(tmp, by=temp_fact, k = 4) +
                s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(location, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)

# add base 'ds' for coordinates
m.gam4 <- bam(sum_ips ~ 
                s(tmp, by=temp_fact, k = 4) +
                s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(location, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)

# add year variable to account for autocorrelation
m.gam5 <- bam(sum_ips ~ 
                s(tmp, by=temp_fact, k = 4) +
                s(year,bs = 're', k = 3) +
                s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(location, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)


# add lagged variables to account for temp autocorrelation
m.gam6 <- bam(sum_ips ~ 
                s(previous_sum_ips, k = 3) +
                s(tmp, by=temp_fact, k = 4) +
                s(year,bs = 're', k = 3) +
                s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(location, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)


### Evaluated GAM -------------------------------------------------------------
AICc(m.gam1, m.gam2, m.gam3, m.gam4, m.gam5,m.gam6 ) 
windows()
plot(m.gam6, page = 1, shade = T)
appraise(m.gam6)
summary(m.gam6)

## Compare GAM, GLM: include anomalies categories or not? -----------

# -> conclusion: go for simpler model, do not split the data in periods
AIC(m3.6, m3.7.1, m3.7.2, m3.7.3.poly, m.gam0, m.gam5)
AICc(m3.6, m3.7.1, m3.7.2, m3.7.3.poly, m.gam0, m.gam5)
BIC(m3.6, m3.7.1, m3.7.2, m3.7.3.poly, m.gam0, m.gam5)

# the best models:
m.gam5 # lowest AIC, AICc, BIC, df 130 = very complex, r2 = 
m3.7.1 # less complex (df = 7)
m3.6   # less complex (df = 5)



r2(m.gam5) # 0.50
r2(m3.7.1) # 0.58 
r2(m3.6)   # 0.54 less complex (df = 5)


### Test autocorrelation of residuals ------------------------------------------------
acf(residuals(m3.7.1))  # seems relatively ok
pacf(residuals(m3.7.1))

plot(allEffects(m3.7.1))

# different ways to handle temporal autocorrelation --------------------------
# wrong!! low AIC, but issueas with family, and not corect residuals
m.gamm <- gamm(sum_ips ~ s(tmp, by = temp_fact, k = 4) +
                 s(year_num, by = temp_fact, bs = "re") +
                 s(location, bs = "re") +
                 s(trap_name, bs = "re"),
               data = dat_lag_scaled, 
               family = nb())
appraise(m.gamm$gam)
plot(m.gamm$gam, page = 1)

acf(residuals(m.gamm$gam))
pacf(residuals(m.gamm$gam))

# 1. if data are overdispersed, need to use neg bin family

m.glmmTMB <- glmmTMB(sum_ips ~ poly(tmp, 2) * temp_fact + 
                       (1 | year_num:temp_fact) + 
                       (1 | location) + 
                       (1 | trap_name),
                     data = dat_lag_scaled, 
                     family = nbinom2)

# 2. handle overdispersion by add observation level ID as random effect
m.gam1 <- gam(sum_ips ~ s(tmp, by = temp_fact, k = 4) +
               s(location, bs = "re") +
               s(trap_name, bs = "re") +
               s(id, bs = "re"),  # Observation-level random effect
             data = dat_lag_scaled, 
             family = nb())

m.gam2 <- gam(sum_ips ~ s(tmp, by = temp_fact, k = 4) +
                s(location, bs = "re") +
                s(trap_name, bs = "re") +
                s(id, bs = "re"),  # Observation-level random effect
              data = dat_lag_scaled, 
              family = poisson())

# the best!
m.gam3 <- gam(sum_ips ~ s(tmp, by = temp_fact, k = 4) +
                s(location, bs = "re") +
                s(trap_name, bs = "re") +
                s(id, bs = "re"),  # Observation-level random effect
              data = dat_lag_scaled, 
              family = tw())


  plot_smooth(m.gam3, 
                          view="tmp", 
                          rm.ranef =TRUE, 
                          #transform=function(x) 100*inv.logit(x),  
                          plot_all="temp_fact")$fv




###  GAM shown in a plot --------------------------------------------------------

# Select specific levels for id, location, and trap_name
selected_id <- unique(dat_lag_scaled$id)[5]      # Replace 1 with your chosen level
selected_location <- unique(dat_lag_scaled$location)[5]
selected_trap <- unique(dat_lag_scaled$trap_name)[5]

# Generate Predicted Values
new_data <- expand.grid(temp_fact = unique(dat_lag_scaled$temp_fact),
                        tmp = seq(min(dat_lag_scaled$tmp), max(dat_lag_scaled$tmp), length.out = 100),
                        location = selected_location,
                        trap_name = selected_trap,
                        id = selected_id)

# Predict using your model
new_data$predicted <- predict(m.gam3, newdata = new_data, type = "response")

# Create the Plot
ggplot(new_data, aes(x = tmp, y = predicted, color = temp_fact)) +
  geom_line() +
  labs(title = "Predicted sum_ips by Temp Factor",
       x = "Temperature",
       y = "Predicted sum_ips") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

## GAMM test model run, before improving autocorrelation - simpler - does not work!! -----------------
m.gamm1 <- gamm(
  sum_ips ~ s(tmp, by = temp_fact, k = 4) +
    s(location, bs = "re") +
    s(trap_name, bs = "re") +
    s(id, bs = "re"),
  data = dat_lag_scaled,
  family = nb()
)

# include autocorrelation; make sure that 
m.gamm2 <- gamm(
  sum_ips ~ s(tmp, by = temp_fact, k = 4) +
    s(location, bs = "re") +
    s(trap_name, bs = "re") +
    s(id, bs = "re"),
  data = dat_lag_scaled,
  family = nb(),
  correlation = corAR1(form = ~ year_num)
)

AICc(m.gam1,m.gam2,m.gam3)
appraise(m.gam3)

acf(residuals(m.gam3))
pacf(residuals(m.gam3))


AICc(m.gamm,m.glmmTMB, m.gam)


# test only for individual years, year needs to be treated as a factor for AR structure!
# does not work
m.gamm.ar <- gamm(
  sum_ips ~ s(tmp, k = 4),# +
   # s(location, bs = "re") +
  #  s(trap_name, bs = "re") +
   # s(id, bs = "re"),
  data = dat_lag_scaled,
  #family = nb(),
  correlation = corAR1(form = ~ year_num)
)



# using mixed effects GAMM with explicit temporal autocorrelation -------------
# Example of using gamm for temporal autocorrelation

# to account for autocorrelation using corAR, I need to have one year per category
m.gamm1 <- gamm(sum_ips ~ s(tmp, by=year, k = 4), 
                  s(location, bs = "re") +
                  s(trap_name, bs = "re"),
                data = dat_lag_scaled, 
                #family = nb(), 
                correlation = corAR1(form = ~ year_num))







m.gamm1 <- gamm(sum_ips ~ s(tmp, by=temp_fact, k = 4), 
               s(year_num, by = temp_fact, bs = "re") +
                 s(location, bs = "re") +
                 s(trap_name, bs = "re"),
               data = dat_lag_scaled, 
               #family = nb(), 
               correlation = corAR1(form = ~ year_num))

# Perform the Durbin-Watson test - better for linear models? ACF,pACF are better for gams, GAMM (gams with mixed effects)
dw_result <- dwtest(residuals(m.gam5) ~ fitted(m.gam5))
(dw_result)




# GAM first try, unscaled----------------------------------------------------------------

m1 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            tmp +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)

m2 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            previous_tmp + 
            tmp +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)

m3 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)


m4 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
           # s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)


m5 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)

m6 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)


m7 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            #s(elev) +
            s(spruce_freq) +
                        tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, tmp),
          dat_lag, family = tw)

m8 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_freq) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, tmp),
          method="fREML", 
          dat_lag, family = tw)


m9 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_freq) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
             s(x,y, k = 25, bs = 'ds') +
            s(year, tmp),
          method="fREML", 
          dat_lag, family = tw)

m9.nb <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_freq) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x,y, k = 25, bs = 'ds') +
            s(year, tmp),
          method="fREML", 
          dat_lag, family = nb)

m9.poisson <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_freq) +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're') +   # trap pair
               s(x,y, k = 25, bs = 'ds') +
               s(year, tmp),
             method="fREML", 
             dat_lag, family = poisson)

m9.poisson.sc <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
                    s(year, k = 6) +
                    s(elev) +
                    s(spruce_freq) +
                    tmp +
                    spei +
                    s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                    s(location, bs = 're') +   # trap pair
                    s(x,y, k = 25, bs = 'ds') +
                    s(year, tmp),
                  method="fREML", 
                  dat_lag, family = poisson, scale = -1)

# add population growth
m10 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_freq) +
            s(population_growth) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, tmp),
          dat_lag,method="fREML",  family = tw)

m10.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(tmp, year),
           dat_lag,method="fREML",  family = tw)


# remove s(tmp, year)
m10.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_freq) +
               s(population_growth) +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
               # s(x,y, k = 50) +
              # s(tmp, year),
             dat_lag,method="fREML",  family = tw)

# interaction using tensor# s(tmp, year)
m10.4 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_freq) +
               s(population_growth) +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
              ti(tmp, year),
             dat_lag,method="fREML",  family = tw)


# remove previous beetle sum
m11 <- bam(sum_ips ~ #s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(year, tmp),
           dat_lag,method="fREML",  family = tw)

# elevation as linear
m12 <- bam(sum_ips ~ #s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             elev +
             s(spruce_freq) +
             s(population_growth) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(year, tmp),
           dat_lag,method="fREML",  family = tw)



# add temperature from previous year
m13 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_freq) +
               s(population_growth) +
               s(previous_tmp) +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(tmp, year),
             dat_lag,method="fREML",  family = tw)

#previous year pop growth
m13.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 5) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth2, k = 6) +
             s(previous_tmp) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)


# set terms to linear

# include autocorrelation: A, and check for autocorrelation by XXXX test
m13.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 5) +
               elev +
               s(spruce_freq) +
               population_growth2 +
               previous_tmp +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(tmp, year),
             dat_lag,method="fREML",  family = tw)


m13.3.1 <- bam(sum_ips ~ s(previous_sum_ips, k = 4) +
               s(year, k = 5) +
               elev +
               s(spruce_freq) +
               population_growth2 +
               previous_tmp +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(tmp, year),
             dat_lag,method="fREML",  family = tw)



m13.3.1.temp.ac <- gamm(sum_ips ~ s(previous_sum_ips, k = 4) +
              s(year, k = 5) +
                elev +
                s(spruce_freq) +
                population_growth2 +
                previous_tmp +
                tmp +
                spei +
                s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(location, bs = 're'),
            data = dat_lag,
            family = tw,
            correlation = corAR1(form = ~ year | trap_name))


m1.temp.ac <- gamm(sum_ips ~ s(previous_sum_ips, k = 4) +
                          s(year, k = 5) +
                          #elev +
                          s(spruce_freq) +
                          population_growth2 +
                          previous_tmp +
                          s(tmp) +
                          #spei +
                          s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                          s(location, bs = 're'),
                        data = dat_lag,
                        family = tw,
                        correlation = corAR1(form = ~ year | trap_name))


m2.temp.ac <- gamm(sum_ips ~ #s(previous_sum_ips, k = 4) +
                     #s(year, k = 5) +
                     #elev +
                     #s(spruce_freq) +
                     #population_growth2 +
                     #previous_tmp +
                     s(tmp) +
                     #spei +
                     s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                     s(location, bs = 're'),
                   data = dat_lag,
                   family = tw,
                   correlation = corAR1(form = ~ year | trap_name))



# check for autocorrelation in residuals:
acf(residuals(m1.temp.ac$lme))  # Final model!!!

appraise(m1.temp.ac$gam)
plot(m1.temp.ac$gam, page = 1, shade = T)
summary(m1.temp.ac$gam)
AIC(m13.3.1.temp.ac,m13.3.1, m1.temp.ac )


#!!!! glm, model averaging???



# add previous year spei 
m14 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth) +
     s(previous_tmp)+
             s(previous_spei) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)


# add temperature from previous year
m13.nb <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth) +
             s(previous_tmp) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = nb(theta = 1.8))


# use as linear terms if the edf = 1
m15 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             year +
             elev +
             spruce_freq +
             s(population_growth) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


# use as linear terms if the edf = 1
m16 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             year +
             elev +
             spruce_freq +
             s(population_growth, k =80) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


# use automatic k - does not work well!
m17 <- bam(sum_ips ~ s(previous_sum_ips, k = TRUE) +
             year +
             elev +
             spruce_freq +
             s(population_growth, k =TRUE) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))




#increase k values
m16.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 130) +
             year +
             elev +
             spruce_freq +
             s(population_growth, k =130) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


#reduce tweedie
m16.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 130) +
               year +
               elev +
               spruce_freq +
               s(population_growth, k =130) +
               previous_tmp +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(tmp, year),
             dat_lag,method="fREML",  family = tw)



# increase k values for m13 model

# work on m 13: inceaase k 
m17 <- bam(sum_ips ~ s(previous_sum_ips, k = 100) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_freq) +
             s(population_growth, k = 50) +
             s(previous_tmp, k = 20) +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)

# set linear vars as linear, not smooths
m18 <- bam(sum_ips ~ s(previous_sum_ips, k = 100) +
             s(year, k = 6) +
             elev +
             spruce_freq +
             s(population_growth, k = 50) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)

# reduce 'k' consistently
m19 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
             s(year, k = 6) +
             elev +
             spruce_freq +
             s(population_growth, k = 50) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)




# add remained forest as smooth
m18.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
             s(year, k = 6) +
             elev +
             #spruce_freq +
             s(remained_forest) +
             s(population_growth, k = 20) +
             previous_tmp +
             tmp +
             spei +
             s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(location, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(tmp, year),
           dat_lag,method="fREML",  family = tw)


# remove non significant and s termsadd remained forest as smooth
m18.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
               s(year, k = 6) +
               elev +
               #spruce_freq +
               remained_forest +
               s(population_growth, k = 20) +
               previous_tmp +
               tmp +
               spei +
               s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(location, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(tmp, year),
             dat_lag,method="fREML",  family = tw)



# put linear trends as linear

summary(m18.2)
gam.check(m17)
k.check(m18)
appraise(m17)
plot(m18, page = 1, shade = T, scale = 0)
AIC(#m1, m2, m3, m4, m5, m7, 
  m8,  # intercept is NaN, which is not correct!
  m9,
  #m9.nb,
  #m9.poisson, 
  #m9.poisson.sc,
  m10, 
  #m10.2, m10.3, m10.4, 
  m13,
  m13.2,
 m13.3 # # !!!! this is the winner! having reasonable error distributions and explaining 50% od deviance
  # m13.nb
  #m14
  #m11, m12
 #m15, 
 #m16, m16.2,
 #m17,
 #m18,
# m18.2
# m19
  )



# Plot the smoothed term with an exponentiated y-axis
plot(m13.3, pages = 1, scale = 0, shade= T)
plot(m8, pages = 1, scale = 0, shade= T)

# Plot explained variability m13.3 -------------------------------------------------


# Extract the F or t statistics as measures of effect size
effect_size <- c(
  previous_sum_ips = 20.257,
  year = 5.850,
  spruce_freq = 2.037,
  trap_name = 0.000, 
  location = 2.572,
  elev = abs(1.800),                # t-value
  population_growth2 = abs(1.325),  # t-value
  previous_tmp = abs(2.787),        # t-value
  tmp = abs(0.840),                 # t-value
  spei = abs(2.952)                 # t-value
)

# Create a data frame for plotting
df <- data.frame(term = names(effect_size), effect_size = unname(effect_size))
df <- df[order(df$effect_size, decreasing = TRUE),]

# Create the bar plot
ggplot(df, aes(x = reorder(term, effect_size), y = effect_size)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  labs(x = "Term", y = "Effect Size (F or |t| value)", title = "Effect Sizes of Terms in GAM Model")







# predict DOY aggregation ---------------------------------------------------

# 
p_agg <- dat_lag_scaled %>% 
ggplot(aes(x = tmp,
           y = agg_doy,
           color = temp_fact)) +
  geom_smooth()


# peak 
p_peak <- dat_lag_scaled %>% 
  ggplot(aes(x = tmp,
             y = peak_doy,
             color = temp_fact)) +
  geom_smooth()

library(ggpubr)
ggarrange(p_agg, p_peak, common.legend = TRUE, legend="bottom", ncol = 2)



# try models
# Normalize DOY to be between 0 and 1
dat_lag_scaled$normalized_DOY_agg <- (dat_lag_scaled$agg_doy - 92) / (300 - 92)
hist(dat_lag_scaled$normalized_DOY_agg)

# try GLM:  -----------------
m.agg1 <- glm(normalized_DOY_agg ~ tmp, 
              #family = betar(link = "logit"),
              data = dat_lag_scaled)


# If this works, gradually add complexity
m.agg2 <- glm(normalized_DOY_agg ~ poly(tmp,2), 
              #family = betar(link = "logit"),
              data = dat_lag_scaled)

m.agg3 <- glmer(normalized_DOY_agg ~ tmp + previous_sum_ips_scal + (1 | year) + (1 | location/trap_name), 
                  #family = 'be,
                  data = dat_lag_scaled)

AIC(m.agg1, m.agg2,m.agg3)



plot(allEffects(m.agg2))






m.agg1 <- gam(agg_doy  ~ s(tmp, by=temp_fact), 
               family = betar(link = "logit"),  # Choosing beta family with logit link
               data = dat_lag_scaled)










# standardize the dependent and predictors vals  -------------------------------------------------------
# Standardize beetle counts for each trap over years
dat_lag_standardized <- 
  dat_lag %>%
    ungroup(.) %>% 
  dplyr::select(-c(trap_name, location, year, x, y)) %>%
  # Apply the scale function to standardize each column
  mutate_all(scale)

# move back removed columns
dat_lag_standardized$trap_name <- dat_lag$trap_name
dat_lag_standardized$location  <- dat_lag$location
dat_lag_standardized$year      <- dat_lag$year
dat_lag_standardized$x         <- dat_lag$x
dat_lag_standardized$y         <- dat_lag$y

hist(dat_lag_standardized$sum_ips)

dat_lag_standardized2 <- dat_lag_standardized %>% drop_na() %>% as.data.frame()


m8_st <- bam(sum_ips ~ s(previous_sum_ips) +
            s(year, k = 3) +
            s(elev) +
            s(spruce_freq) +
            tmp +
            spei +
            s(trap_name, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(year, tmp),
            method="fREML",
            data = dat_lag_standardized2)


m13.3_st <- bam(sum_ips ~ s(previous_sum_ips, k = 70) + 
               s(year, k = 5) + 
               elev + 
               s(spruce_freq) + 
               population_growth2 + 
               previous_tmp + 
               tmp + 
               spei + 
               s(trap_name, bs = "re") + 
               s(location, bs = "re"),
             method="fREML",
             data = dat_lag_standardized2)



AIC(m8_st, m13.3, m13.3_st)
appraise(m13.3_st)
plot(m13.3_st, page = 1, shade = T)

# 10/17/2023 - after standardizing, the AIC is way lower and calculation is faster, but teh reults are the same
# calculate beetle population growth (compare teh pup. from 
#ne year to another)
# identify transition points: how does this work???
#
# Install and load the changepoint package if not already installed
# install.packages("changepoint")
#library(changepoint)

# Sample beetle counts data (replace with your data)
beetle_counts <- c(100, 120, 150, 80, 100, 500, 800, 900, 1200, 1100, 1150)

# Perform changepoint analysis to detect transitions
cpt <- cpt.meanvar(beetle_counts)

# Plot the data with changepoints
plot(cpt)


plot_smooth(m4, view = 'tmp')
pXXX<- plot_smooth(m1, 
                          view="tmp", 
                          rm.ranef =TRUE)$fv





windows()
plot(m1, page = 1)


dat_fin %>% 
  ggplot(aes(y = sum_ips,
             x = previous_tmp)) + 
  geom_point()



m1 <- gam(sum_ips ~s(previous_sum_ips, k =50) + 
            s(elev) +
            s(conif_prop  ) +
            s(previous_tmp) +
            s(year, k = 5) +
            s(trap_name, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x,y),
          dat_fin, family = nb )

appraise(m1)
summary(m1)
plot(m1, page = 1, shade = T)

# Inspect data ----------------------------------------------------------
# VIF - variance inflation factor
# depepndence between predictors - keep only independent ones
# spatial and temporal autocorrelation

# get korelogram: if line overpass teh dashed line, there is a corelation
acf(dat_fin$sum_ips) # super correlated!





# correlation between predictors ---------------------------------------
dat_fin_cor <- dat_fin %>% 
  dplyr::select(-c(trap_name, location, year))


# Check out which predictors are highly correlated to remove them from teh VIF
preds.res <- cor(dat_fin_cor, use = "complete.obs")
round(preds.res, 2)  # simplify visualization

# p-values are missing: get them from different package:
library("Hmisc")
res2 <- Hmisc::rcorr(as.matrix(dat_fin_cor))
res2


library(corrplot)

# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.05, insig = "blank")



# Check for scatter plots between precistors:
#install.packages("PerformanceAnalytics")

library("PerformanceAnalytics")
chart.Correlation(dat_fin_cor, histogram=TRUE, pch=19)
#chart.Correlation(ips2_predictors_years[[2]], histogram=TRUE, pch=19)








# # variance inflation factor (VIF) ---------------------------------------------


# Compute the VIF values: remove multicollinearity between predictors -------------
# https://www.statology.org/variance-inflation-factor-r/
vif_model_ips <- lm(sum_ips ~ conif_prop + elev +
                    #sm +
                     # vpd +
                      sm_z + 
                      #vpd_z + 
                      #tmp + 
                      tmp_z# +
                      #spei + 
                    #  spei_z
                    , 
                    data = dat_fin)



# get a vector of vif values
(vif_model_ips <- car::vif(vif_model_ips))

# DOY of aggregation
vif_model_agg <- lm(agg_doy ~ conif_prop + elev +
                      sm +
                      vpd +
                      #sm_z + 
                      #vpd_z + 
                      #tmp + 
                      tmp_z# +
                    #spei + 
                    #  spei_z
                    , 
                    data = dat_fin)



# get a vector of vif values
(vif_model_agg <- car::vif(vif_model_agg))



# DOY of peak
vif_model_peak <- lm(peak_doy ~ conif_prop + elev +
                      sm +
                      #vpd +
                      sm_z + 
                      #vpd_z + 
                      #tmp + 
                      tmp_z# +
                   # spei + 
                     # spei_z
                    , 
                    data = dat_fin)



# get a vector of vif values
(vif_model_peak <- car::vif(vif_model_peak))



# find distribution:-Negative binomial -----------------------------------------


# Families in bam: --------------------------------------------------------------
# ocat for ordered categorical data.
# tw for Tweedie distributed data, when the power parameter relating the variance to the mean is to be estimated.
# nb for negative binomial data when the theta parameter is to be estimated.
# betar for proportions data on (0,1) when the binomial is not appropriate.
# scat scaled t for heavy tailed data that would otherwise be modelled as Gaussian.
# ziP for zero inflated Poisson data, when the zero inflation rate depends simply on the Poisson mean.
m <- gam(sum_ips ~ 1,dat_fin, family = nb )  
appraise(m)

# m <- gam(sum_ips ~ 1,dat_fin, family = nb(link = 'log' ))  
# appraise(m)
# 
 m <- gam(sum_ips ~ 1,dat_fin, family = nb(link = 'log', theta = 1.82 ))  
 appraise(m)

 
 m <- gam(sum_ips ~ 1,dat_fin, family = nb(link = 'log', theta = NULL ))  
 appraise(m)
 
 # tested several thetas, the 1.8 is teh best 






# Do drivers importance change over time? --------------
# ChatGPT: 
# To test if the importance of drivers explaining variability changed over time: 
# perform temporal interaction analysis - 
# allows to assess whether the relationships between the drivers 
# and the response variable differ across different time periods. 

# check hist of conts
windows()
hist(dat_fin$sum_ips)
median(dat_fin$sum_ips)
mean(dat_fin$sum_ips)
sd(dat_fin$sum_ips)


ggplot(dat_fin, aes(x = sum_ips)) + 
  geom_histogram(colour = 4, fill = "white", 
                 bins = 2000)

# check for 0
dat_fin %>% 
  filter(sum_ips == 0) %>% # 0 no!! 
  distinct(falsto_name)
# no 0

fitdistr(dat_fin$sum_ips, 'Poisson')


# poisson: mean and variance are equal
# negative-binomial - allows for overdispersion (variance excees the mean) - not my case

glm1 <- glm(sum_ips ~ conif_prop + elev + sm_z + tmp_z +(1 | trap_name),
            data = dat_fin,
            family = "poisson",
            na.action = "na.fail")

# Explore teh model summary:
summary(glm1)
anova(glm1)
MuMIn::dredge(glm1)

# explore residuals
windows()
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm1, 
                                              plot = T)
# Yay! the QQ plot shows that we have some issues, residuals are not independent

# test if the residuals are the same as random?
testDispersion(simulationOutput)















# try the prediction using the raw data --------------------------------

M <- list(c(1, 0.5), NA)
m_counts <- bam(fangmenge ~
                  s(spei1, k = 80)+ # + # drought
                  s(TMED, k = 30) +  # temperature
                  #s(PRCP, k = 20) +         # precip 
                  s(freq, k = 50), #+         # spruce %
                #s(month, k =7) +  # months
                #s(year, k = 8) +  # year
                #s(x, y, k = 10, bs = 'ds', m = c(1, 0.5)) + # 2D smooth
                #ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
                #ti(x, y, TMED, d = c(2,1), bs = c('ds','tp'), 
                #   m = M, k = c(25, 10)) +
                #ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
                #  m = M, k = c(25, 15)) +
                #ti(x, y, spei, d = c(2,1), bs = c('ds','tp'),
                #  m = M, k = c(25, 15)),
                data = ips_sum2, 
                method = 'fREML',
                #select = TRUE,
                #family = Tweedie(p=1.1, link = power(0)),
                tw(link = "log"),
                # knots = knots,
                nthreads = 4, 
                discrete = TRUE)

# check for k: edf and k should be balanced ~ 1
m_counts
summary(m_counts)
k.check(m_counts)
plot(m_counts, plot.all = T)


windows()
appraise(m_counts, method = 'simulate')
plot(m5, page = 1, shade = T)

gam.check(m5)
k.check(m5)
summary(m5)






# Example Chat GPT: how to analyze data?


dat_fin_sub <- dat_fin %>%
  filter(location %in% c('Weismain','Pressig', 'Weidenberg' )) #%>% 



dat_fin_sub %>%   print()
  


# Check correlation -------------------------------------------------------

cor(dat_fin[,c('sm', 'vpd', 'tmp', 'spei')], use = "complete.obs")




















# Test XY relationship pattern: observed vs modelled data ------------------
#


# add temporal autocorrelation  
dat_fin$AR.START <-dat_fin$year==2015


# What distribution? ---------------------------------------------
# https://towardsdatascience.com/insurance-risk-pricing-tweedie-approach-1d71207268fc
# poisson - only for counts
# gamma - does not take zero values
# tweedie - can handle zeros values


# useful to scale values?
scale(c(1:5))



# Get model with gam: ------------------------------------------------

# get predictors following VIF:
# > vif_model_ips
# conif_prop       elev       sm_z      tmp_z 
# 1.117820   1.230662   3.025474   3.025275 


M <- list(c(1, 0.5), NA)

m0 <- gam(sum_ips ~ 1
          ,
          method = 'fREML',
          data = dat_fin,
          family = nb,
          nthreads = 4
)

appraise(m0)


dat_fin$sum_ips

# Work on this one!!!!!

# add predictors: counst from previous year

m1 <- bam(sum_ips ~ s(year, k = 5) +
            s(tmp) +
            te(year, tmp),
          data = dat_fin,
          family = nb(theta = NULL))


m2 <- bam(sum_ips ~ s(year, k = 5) +
            s(tmp),# +
  #          te(year, tmp),
          data = dat_fin,
          family = nb(theta = NULL))

m3 <- bam(sum_ips ~ s(year, k = 5) +
            s(tmp) +
          s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m4 <- bam(sum_ips ~ s(year, k = 5) +
            s(tmp, k = 15) +
            s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m5 <- bam(sum_ips ~ s(year, k = 5) +
            tmp +
            s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m6 <- bam(sum_ips ~ s(year, k = 5) +
            tmp +
            s(trap_name, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 5))

m7 <- bam(sum_ips ~ s(year, k = 5) +
            tmp +
            s(trap_name, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m8 <- bam(sum_ips ~ s(year, k = 5) +
            tmp + 
            s(elev) +
            s(trap_name, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(location, bs = 're') +   # trap pair
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m9 <- bam(sum_ips ~ s(year, k = 5) +
            tmp + 
            elev +
            s(trap_name, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(location, bs = 're') +   # random effect for trap pair 
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m10 <- bam(sum_ips ~ s(year, k = 5) +
            tmp + 
            elev +
             spei +# spei does not improve the model
            s(trap_name, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(location, bs = 're') +   # random effect for trap pair 
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m11 <- bam(sum_ips ~ s(year, k = 5) +
             conif_prop +
             tmp + 
             elev +
             s(trap_name, bs = 're') +  
             s(location, bs = 're') +    
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))


m12 <- bam(sum_ips ~ s(year, k = 5) +
            s(conif_prop) +              # does not improve, neither as s(), or linear term
             tmp + 
             elev +
             s(trap_name, bs = 're') +  
             s(location, bs = 're') +   
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))

m13 <- bam(sum_ips ~ s(year, k = 5) +
             s(conif_prop) +              # does not improve, neither as s(), or linear term
             tmp +
             s(tmp, by = year, k = 7) +   # changes in temperature over year
             elev +
             s(trap_name, bs = 're') +  
             s(location, bs = 're') +   
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))


m14 <- bam(sum_ips ~ s(year, k = 5) +
             s(conif_prop) +              # does not improve, neither as s(), or linear term
             tmp +
             s(year, by = tmp, k = 7) +   # changes year over temperature 
             elev +
             s(trap_name, bs = 're') +  
             s(location, bs = 're') +   
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))



appraise(m14)
AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13,m14)
summary(m14)
plot(m14, page = 1, shade= T)




# cap outliers! ----------------------------------------------------------------
lower_cap <- quantile(dat_fin$sum_ips, 0.01)
upper_cap <- quantile(dat_fin$sum_ips, 0.99)


dat_fin_cap <- dat_fin 

dat_fin_cap$sum_ips <- ifelse(dat_fin_cap$sum_ips < lower_cap, lower_cap, 
                   ifelse(dat_fin_cap$sum_ips > upper_cap, upper_cap, dat_fin_cap$sum_ips))


hist(dat_fin_cap$sum_ips)
hist(dat_fin$sum_ips)


# try with capped values
m1 <- bam(sum_ips ~ s(year, k = 7) +
            tmp + 
            elev +
            s(trap_name, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(location, bs = 're') +   # random effect for trap pair 
            s(x, y, k =120),  
          data = dat_fin_cap,
          family = nb(theta = 3))

appraise(m1)
summary(m1)
k.check(m1)








m1 <- bam(sum_ips ~ 
           s(conif_prop, k = 9) +
            s(elev, k = 9) +
           
            s(spei, k = 8) +
            s(vpd, k = 8) +
           # s(vpd_z, k = 8) +
            s(spei_z, k = 8) +

            s(tmp_z, k = 8) +
            s(tmp, k = 8) +
           # s(sm_z, k = 8) +
            s(sm, k = 8) +
            #s(tmp, k = 20) +
            s(year, k = 5) +

            s(location, bs = "re") +
            s(trap_name, bs = "re") +
            s(x,y, k = 20, bs = 'ds') +
            ti(x,y, year, d = c(2,1), 
               bs = c('ds','tp'), m = M,
               k = c(20, 5)) +
            ti(x,y, tmp, d = c(2,1), 
               bs = c('ds','tp'), m = M,
               k = c(20, 5)),
          method = 'fREML',
          #AR.start=AR.START, 
          #rho=0.15,
         data = dat_fin,
         family = nb(theta = NULL),  # mean = 22000, sd = 19000
         nthreads = 4
)
appraise(m1)
summary(m1)

plot(m1, pages = 1, scheme = 2, shade = TRUE, scale = 0)



# use step-wise model selection ------------------------------------------------
library(MASS)

dat.small <- droplevels(dat_fin[dat_fin$trap_name %in% 
                                  levels(dat_fin$trap_name)[10:30],]) #len prvych par kombinacii site & zone


full.model <- bam(sum_ips ~ 
                          s(conif_prop, k = 9) +
                          s(elev, k = 9) +
                          
                          s(sm, k = 8) +
                          s(sm_z, k = 8),## +
                          
                  #  s(spei, k = 8) +
                   # s(spei_z, k = 8) +
                    
                     #     s(vpd, k = 8) +
                    #      s(vpd_z, k = 8) +
                         
                   #       s(tmp, k = 8) +
                   #       s(tmp_z, k = 8) +
                    
                   #       s(year, k = 6) +
                  
                   #       s(location, bs = "re") +
                   #       s(trap_name, bs = "re") +
                  #        s(x,y, k = 20, bs = 'ds') +
                   #       ti(x,y, year, d = c(2,1), 
                  #           bs = c('ds','tp'), m = M,
                  #           k = c(20, 5))          ,
                        method = 'fREML',
                  AR.start=AR.START, 
                  rho=0.15,
                        data = dat.small,
                        family = negative.binomial(theta = 1.8),  # mean = 22000, sd = 19000
                        nthreads = 4
)


step.model <- stepAIC(full.model, 
                      scope = list(upper = full.model, lower = ~1), 
                      direction = "both", 
                      trace = TRUE)










# Example stacklk
df <- data.frame(my_sum = rnbinom(24, 1, 0.1),
                 location = rep(c('a','b', 'c', 'd'),each =6),
                # trap_n = rep(c(1,2),12),
                # trap_ID = paste(location, trap_n),
                 year = 1:6,
                 x = c(1, 9,110,119,210,219,310,319)*100,
                 y = c(1,9,1,9,1,9,1,9)*100)

traps <- matrix(c(
  100, 100,
  900, 900,
  11000, 100,
  11900, 900,
  21000, 100,
  21900, 900,
  31000, 100,
  31900, 900
), ncol = 2, byrow = TRUE)
(df)

m1 <- gam(my_sum ~ s(vars1) + s(vars2) + 
            s(year) +
            s(location, bs = 're') +
            s(location, trap_n, bs = 're'),
            data = rats, method = 'REML')







m3 <- bam(sum_ips ~ 
            s(conif_prop, k = 9) +
            s(elev, k = 9) +
            s(sm, k = 8) +
            s(spei, k = 8) +
            s(trap_name, bs = 're') +
          
            #s(trap_name, location, bs = 're') +
            #s(vpd, k = 8) +
            # s(vpd_z, k = 8) +
            # s(spei, k = 8) +
            #s(tmp_z, k = 8) +
            # s(tmp, k = 8) +
            # s(sm_z, k = 8) +
            #s(tmp, k = 20) +
            s(year, k = 5) +
            s(x,y, k = 20, bs = 'ds') +
            ti(x,y, year, d = c(2,1), 
               bs = c('ds','tp'), m = M,
               k = c(20, 5))
          ,
          AR.start=AR.START, 
          rho=0.15,
          method = 'fREML',
          data = dat_fin,
          family = nb,
          nthreads = 4
)






interaction(c('a','b','d'), c('x', 'y', 'z'))







m2 <- gam(sum_ips ~ 
            s(conif_prop, k = 9) +
            s(elev, k = 9),
          data = dat_fin,
          family = nb
)

m3 <- gam(sum_ips ~ 
            s(conif_prop, k = 9) +
            s(elev, k = 9) +
            s(sm_z),
          data = dat_fin,
          family = nb
)
appraise(m1)
gam.check(m3)
AIC(m, m1)



M <- list(c(1, 0.5), NA)
m2 <- bam(
  sum_ips ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(conif_prop       , k = 5) +
    s(elev, k = 5) +  # elevation
    s(sm_z, k = 5)  +         # soils moisture
    s(tmp_z, k = 5) +  # temperature
    s(year, k = 5)+  # year
    s(trap_name, k = 5, bs = 're'),# + #random effect of trap
    # s(location, k = 5, bs = 're') + #random effect of trap pairs
    # s(
    #   x,
    #   y,
    #   k = 5,
    #   bs = 'ds',
    #   m = c(1, 0.5)
    # ) + # 2D smooth
    # #ti(tmp_z, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    # ti(
    #   x,
    #   y,
    #   tmp_z,
    #   d = c(2, 1),
    #   bs = c('ds', 'tp'),
    #   m = M,
    #   k = c(20, 10)
    # ) +
    # ti(
    #   x,
    #   y,
    #   year,
    #   d = c(2, 1),
    #   bs = c('ds', 'tp'),
    #   m = M,
    #   k = c(20, 5)
    # ),#  +
  # ti(
  #   x,
  #   y,
  #   spei1,
  #   d = c(2, 1),
  #   bs = c('ds', 'tp'),
  #   m = M,
  #   k = c(20, 15)
  # )    ,
  data = dat_fin,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = nb,
 # knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)




appraise(m2)













# standardize teh data: Z score: mean = 0, sd = 1 ---------------
# will ths lead to better fit? no!
ips_standard <- ips_sum2 %>% 
  mutate(ips_sum = as.vector(scale(ips_sum)),
         PRCP = as.vector(scale(PRCP)),
         TMED = as.vector(scale(TMED)),
         spei1 = as.vector(scale(spei1)),
         spei3 = as.vector(scale(spei3)),
         spei6 = as.vector(scale(spei6)),
         spei12 = as.vector(scale(spei12)),
         freq = as.vector(scale(freq) )) # %>% 
# filter(!is.na) %>% 
#complete.cases(.) %>% 
as.data.frame()



# Test gams with standardized data: ----------------------
m1 <- gam(ips_sum ~ s(spei3),
          data = ips_standard,
          family = tw(link = 'log') )

### GAM standardized y value: -------------
m1 <- gam(ips_sum_stan ~ s(spei3), 
          data = ips_sum2)

length(fitted(m1))
length(ips_standard$spei3)

appraise(m1)

ggplot(ips_sum2, aes(x = spei3, y = ips_sum_stan)) +
  geom_point(alpha=0.7)+
  geom_line(colour = "red", size = 1.2,
            aes(y = fitted(m1))) #+
geom_line(colour = "blue", size = 1.2,
          aes(y = fitted(m1))) +
  theme_bw()


# GAM: raw sums: monthly sums -------------------------------------------
knots_month <- list(doy = c(4.2, 10.5))

M <- list(c(1, 0.5), NA)
m7 <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 5) +  # temperature
    s(PRCP, k = 5)  +         # precip
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(month, k = 7, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    # s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      year,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 5)
    ),#  +
  # ti(
  #   x,
  #   y,
  #   spei1,
  #   d = c(2, 1),
  #   bs = c('ds', 'tp'),
  #   m = M,
  #   k = c(20, 15)
  # )    ,
  data = ips_sum2,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)







ips_sum2_sept <- filter(ips_sum2, month < 10)

M <- list(c(1, 0.5), NA)
m_rem10 <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 20) +  # temperature
    s(PRCP, k = 20)  +         # precip
    s(month, k = 6, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      year,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 5)
    ),#  +
  # ti(
  #   x,
  #   y,
  #   spei1,
  #   d = c(2, 1),
  #   bs = c('ds', 'tp'),
  #   m = M,
  #   k = c(20, 15)
  # )    ,
  data = ips_sum2_sept,
  AR.start=AR.START, 
  rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)


plot(m_rem10, page = 1, scale = 0, shade = T)
summary(m_rem10)
appraise(m_rem10)

draw(m_rem10, select = 13)  # xy



# m7.noAR --------------------------------------------------------
m7.noAR <- bam(
  ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(TMED, k = 5) +  # temperature
    s(PRCP, k = 5)  +         # precip
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(month, k = 7, bs = 'cc') +  # months
    s(year, k = 5)+  # year
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(
      x,
      y,
      k = 5,
      bs = 'ds',
      m = c(1, 0.5)
    ) + # 2D smooth
    ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    ti(
      x,
      y,
      TMED,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 10)
    ) +
    ti(
      x,
      y,
      PRCP,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 15)
    )  +
    ti(
      x,
      y,
      spei1,
      d = c(2, 1),
      bs = c('ds', 'tp'),
      m = M,
      k = c(20, 15)
    ),
  data = ips_sum2,
  #AR.start=AR.START, 
  # rho=0.15,
  method = 'fREML',
  family = tw(),
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)

# AIC --------------------------------------------
AIC(m7, m7.noAR )


# m7.sub <- bam(
#   ips_sum ~ #spei_cat + s(year, by=spei_cat, k=6)  +
#     #s(spei1, k = 5)  + # drought
#     s(spei3, k = 5) +
#     s(spei6, k = 5) +
#     s(spei12, k = 5) +
#     s(TMED, k = 5) +  # temperature
#     #s(PRCP, k = 5)  +         # precip
#     s(freq, k = 5) +         # conif %
#    # s(sp_prop) +             # spruce proportion
#     s(month, k = 7, bs = 'cc') +  # months
#     s(year, k = 5)+  # year
#     s(globalid, k = 5, bs = 're') + #random effect of trap
#     s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
#     s(
#       x,
#       y,
#       k = 5,
#       bs = 'ds',
#       m = c(1, 0.5)
#     ) + # 2D smooth
#     ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
#     ti(
#       x,
#       y,
#       TMED,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 10)
#     ) +
#     ti(
#       x,
#       y,
#       PRCP,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     )  +
#     ti(
#       x,
#       y,
#       spei1,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     ) +
#     ti(
#       x,
#       y,
#       ,
#       d = c(2, 1),
#       bs = c('ds', 'tp'),
#       m = M,
#       k = c(20, 15)
#     )
#   data = ips_sum2,
#   AR.start=AR.START, 
#   rho=0.15,
#   method = 'fREML',
#   family = tw(),
#   knots = knots_month ,
#   nthreads = 4,
#   discrete = TRUE
# )
# 
# 
# 



# compare more complex with less complex model: which one is better?
AIC(m7, m7.sub)


windows()
summary(m7)
plot(m7, page = 1, shade = T, scale = 0)
appraise(m7)
gratia::draw(m7, select = 13)

windows()
gratia::draw(m7, select = 16)


# Make some spatial predictions: !!! does not work!!!
ips_sum3 <- ungroup(ips_sum2)
exp_data <- with(ungroup(ips_sum3),
                 expand.grid(month = 6,
                             #DoY = 180,
                             year = seq(min(year), max(year), by = 1),
                             x  = seq(min(x), max(x), length = 100),
                             y  = seq(min(y), max(y), length = 100)))
head(exp_data)
fit <- predict(m7, exp_data)

hist(predict(m7$fitted.values))
save.image("models.RData")




# Test small model:
m8 <- gam(ips_sum ~ 
            s(month, k = 6) +
            s(year, k = 7) + 
            s(TMED) + 
            s(x,y),
          data = ips_sum2,
          family = tw()
)

summary(m8)
appraise(m8)


m.logy <- bam(
  log_sum  ~
    s(month, k = 6, bs = 'cc') +
    #s(year, k = 7) +
    s(TMED, k = 80) +
    s(PRCP, k = 10) +
    s(spei1, k = 5)  + # drought
    s(spei3, k = 5) +
    s(spei6, k = 5) +
    s(spei12, k = 5) +
    s(freq, k = 5) +         # conif %
    s(sp_prop) +             # spruce proportion
    s(globalid, k = 5, bs = 're') + #random effect of trap
    s(monsto_name, k = 5, bs = 're') + #random effect of trap pairs
    s(x, y), # +
  data = ips_sum2,
  AR.start = AR.START,
  rho = 0.15,
  family = scat,
  method = 'fREML',
  knots = knots_month ,
  nthreads = 4,
  discrete = TRUE
)

appraise(m.logy)
summary(m.logy)
plot(m.logy, page = 1, )

AIC(m7, m.logy)



# 
exp_data <- with(ips_sum2,
                 expand.grid(month = 6,
                             #DoY = 180,
                             year = seq(min(year), max(year), by = 1),
                             x  = seq(min(x), max(x), length = 100),
                             y  = seq(min(y), max(y), length = 100)))

exp_data<- exp_data %>% 
  left_join(ips_sum2)


exp_data %>% 
  filter(!is.na(globalid))

head(exp_data)
fit <- predict(m8, exp_data)



# with function -----------------------------------------------------
str(mtcars)

a2 <- with(mtcars, mpg[cyl == 8  &  disp > 350])




# Make maps: --------------------------------------------------------------
library(tidyverse)
# ips_sum2 <- ips_sum2 %>% 
#   mutate(x = as.numeric(round(x/100000, 2)),
#          y = as.numeric(round(y/100000, 2)))#

ips_sum2 <- ips_sum2  %>% 
  mutate(across(c(x, y), round, digits = 4))


ggplot(ips_sum2, 
       aes(x = x, y = y)) +
  geom_raster(aes(fill = ips_sum/10000))  + 
  facet_wrap(~ year_month, ncol = 3) # +
scale_fill_viridis(name = expression(degree*C), option = 'plasma',
                   na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'top', legend.key.width = unit(2, 'cm'))


# Galveston
ggplot(pred.g, 
       aes(x = LONGITUDE, 
           y = LATITUDE)) +
  geom_raster(aes(fill = Fitted)) + facet_wrap(~ YEAR, ncol = 12) +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma',
                     na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'top', legend.key.width = unit(2, 'cm'))

# test for Tweedie family ---------------------------- --------------------

library(mgcv)
set.seed(3)
n<-400
## Simulate data...
dat <- gamSim(1,n=n,dist="poisson",scale=.2)
dat$y <- rTweedie(exp(dat$f),p=1.3,phi=.5) ## Tweedie response

## Fit a fixed p Tweedie, with wrong link ...
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
         data=dat)
plot(b,pages=1)
print(b)
gam.check(b)
gratia::appraise(b, method = 'simulate')

## Same by approximate REML...
b1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=Tweedie(1.25,power(.1)),
          data=dat,method="REML")
plot(b1,pages=1)
print(b1)
gam.check(b1)
gratia::appraise(b1, method = 'simulate')


## estimate p as part of fitting

b2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=tw(),
          data=dat,method="REML")
plot(b2,pages=1)
print(b2)

gam.check(b2)
gratia::appraise(b2, method = 'simulate')


rm(dat)


# check the k-values: bases functions:
gam.check(m6)
summary(m6)
# plot to have all 4 plots on one page

plot(m5, pages = 1, scheme = 2, shade = TRUE) # , scale = 0

#summary(m4)
gratia::appraise(m5, method = 'simulate')

# check the k value:

gratia::draw(m5, scales = 'free')

plot(m5, pages = 1, scheme = 2, shade = TRUE)




anova(m1, m2, m3)
AIC(m1, m2, m3)




# ---------------------------------------------------------
# test gam models on full catch dataset 
# --------------------------------------------------------

m1 <- gam(fangmenge ~ s(doy, k = 100, bs = 'cc') +
            s(year, k = 8), #+
          #s(globalid, k = 300),
          data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          methods = 'REML')
windows()
plot(m1, shade = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m1, method = 'simulate')

# check the k value:
gam.check(m1)


knots <- list(doy = c(0.5, 365.5))
m2 <- bam(fangmenge ~ s(doy, k = 110, bs = 'cc') +
            s(year, k = 8), # +
          data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)



m3 <- bam(fangmenge ~ s(doy, k = 110, bs = 'cc') +
            s(year, k = 8) +
            s(globalid, bs = 're'),  # 're' = random effect
          data = ips,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m3, method = 'simulate')

# check the k value:
gam.check(m3)
plot(m3)
gratia::draw(m3, scales = 'free')

plot(m3, pages = 1, scheme = 2, shade = TRUE)

anova(m1, m2, m3)
AIC(m1, m2, m3)




# GAM with medians counts per month ---------------------------------------

# still have a lot of wiggliness over the ips catch data:
# maybe I can use the monthly medians?
# will fit with the temperature data as well

# 
ips_med <- ips %>% 
  group_by(year, month, globalid) %>% 
  summarize(fang_med = median(fangmenge, na.rm = T)) %>% 
  filter(year> 2014) %>% 
  as.data.frame()


str(ips_med)

# 
knots_month <- list(doy = c(3.5, 10.5))
m_med1 <- gam(fang_med ~ s(month, k = 7, bs = 'cc'),#+
              # s(year, k = 7), #+
              #s(globalid, k = 300),
              knots = knots_month,    # start and end of the cyclic term 'cc'
              data = ips_med,
              family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
              methods = 'REML')

gam.check(m_med1)




windows()
plot(m_med1, shade = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m_med1, method = 'simulate')

# check the k value:
gam.check(m_med1)



m_med2 <- bam(fang_med ~ s(month, k = 7, bs = 'cc')+
                s(year, k = 7), # +
              data = ips_med,
              family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
              method = 'fREML', # fast REL 
              knots = knots_month,    # start and end of the cyclic term 'cc'
              nthreads = c(4, 1), 
              discrete = TRUE)



m_med3 <- bam(fang_med ~ s(month, k = 7, bs = 'cc')+
                s(year, k = 7) +
                s(globalid, k = 500, bs = 're'),  # 're' = random effect
              data = ips_med,
              family = nb, #,nb,#tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
              method = 'fREML', # fast REL 
              knots = knots_month,    # start and end of the cyclic term 'cc'
              nthreads = c(4, 1), 
              discrete = TRUE)

# plot to have all 4 plots n one page
gratia::appraise(m_med3, method = 'simulate')
gam.check(m_med3, pages = 1)

summary(m_med3)

# check the k value:
gratia::draw(m_med3, scales = 'free')

plot(m3, pages = 1, scheme = 2, shade = TRUE)

anova(m1, m2, m3)
AIC(m1, m2, m3)



