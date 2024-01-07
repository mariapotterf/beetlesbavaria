


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


# # check if results are correct
# dat_lag %>%
#   filter(trapID== 'Zusmarshausen_1') %>%
#   dplyr::select(c(trapID, year, sum_ips, previous_sum_ips,previous_sum_ips2,
#                   veg_tmp,
#                   previous_Moran,previous_Moran2,
#                   Morans_I,
#                   population_growth, population_growth2,
#                   previous_veg_tmp)) #%>%
#   #View()
# 

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





##### IPS sparmans -----------  
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




## Scale predictors ================================================================================

# skip columns if not for scaling
#spei_cols <- grepl("spei", names(dat_lag))

additional_skip_cols <- c('trapID', 'pairID', "x", "y", 'year', 'agg_doy', 'peak_doy', 'peak_diff', 
              'sum_ips', 'previous_sum_ips', 'previous_sum_ips2')

# Combine spei columns with additional columns to skip
skip_col <- additional_skip_cols #c(names(dat_lag)[spei_cols], additional_skip_cols)

# scale also SPEIs

# export new df
dat_lag_scaled <-
  dat_lag %>%
  ungroup(.) %>% 
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






# Analyses =====================================================================

# beetle population sum vs climate drivers
# DOY aggregation vs climate drivers
# DOY peak vs climate drivers
# peak difference vs climate drivers

#### prepare individual tables for each analyses -----------------------------------------------------

windows()
hist(dat_lag_scaled)








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

# Fit a linear regression model to explain Moran's I
m_moran1 <- lm(Morans_I ~ poly(veg_tmp,2) + previous_sum_ips + poly(spei3,2), data = dat_lag_scaled_complete)

plot(allEffects(m_moran1))

summary(m_moran1)
r2(m_moran1)

# lmer: account for the random effect of teh trap: the effect of the trap can vary between years
m_moran2 <- lmer(Morans_I ~ veg_tmp + veg_prcp + (1|trapID), data = dat_lag_scaled)

# remove prec as it is correlated with veg_tmp
m_moran3 <- lmer(Morans_I ~ veg_tmp + (1|trapID), data = dat_lag_scaled)

#add pairID as random
m_moran4 <- lmer(Morans_I ~ veg_tmp + (1|pairID), data = dat_lag_scaled)

# if temp and prec are correlated, maybe there is an interaction effect?
m_moran5 <- lmer(Morans_I ~ veg_tmp * veg_prcp + (1 | trapID), data = dat_lag_scaled)


# try different family: tweedie, gaussian and linear model does not work well

m_moran6 <- glmmTMB(Morans_I ~ veg_tmp + veg_prcp,# + #(1 | trapID), 
                 data = dat_lag_scaled, 
                 family = tweedie)

# add randm effect
m_moran7 <- glmmTMB(Morans_I ~ veg_tmp + veg_prcp+ (1 | pairID), 
                    data = dat_lag_scaled, 
                    family = tweedie)

m_moran8 <- glmmTMB(Morans_I ~ veg_tmp + veg_prcp+ (1 | year), 
                    data = dat_lag_scaled, 
                    family = tweedie)

# add previous year
m_moran9 <- glmmTMB(Morans_I ~ previous_veg_tmp + veg_prcp+ (1 | year), 
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


# ips vs current years:
# inspect extreme values? wind_beetle > 200 - pixels or ha??
dat_lag_scaled %>% 
  filter(wind_beetle > 100)

dat_lag_scaled %>% 
  filter(wind_beetle < 5) %>% 
  ggplot(aes(x = sum_ips,
             y = wind_beetle )) +
  geom_smooth() +
  geom_point()

dat_lag_scaled %>% 
filter(wind_beetle > 100) %>% 
  ggplot(aes(x = sum_ips,
             y = harvest )) +
  geom_smooth() 


# RS explore lagged values 
pairs(wind_beetle ~ sum_ips + previous_sum_ips + previous_sum_ips2 + previous_agg1 + previous_agg2 + population_growth2, dat_lag_scaled)



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
simple_model <- glm.nb(wind_beetle ~ previous_sum_ips +
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
### expore plots ------------------------------------------------------------------------

## test GLM (no random effects) on raw data ---------------------------------------------------------------------

##apparoces: GLM - no random effects (glm.nb)
# random effects : (glmer.nb) - not working really well
# final: use glmmTMB -allow for specifying random structure and dependencies in data collection

# select only columns of interests, to not remove more data as necessary 
# eg lagged values, that were NA
dd <- dat_lag_scaled %>% 
  dplyr::select(c(sum_ips, veg_tmp,
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




# FINAL ONE! start simple with glmmTMB, then add random effects ============================
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



AIC(m1, m2, m3, m4, m5, m6, m7,m8, m9, m9.1, m10, m11, m12, m13)


m9.tmb <- m9

summary(m9)
r2(m9)
simRs <- simulateResiduals(m9, plot = T)
plot(allEffects(m9))
testOutliers(m9)
acf(residuals(m9.tmb))

dw_result <- dwtest(residuals(m9.tmb) ~ fitted(m9.tmb))
(dw_result)



# remove additionsl NAs
dd_complete <- dd %>%
  na.omit()



#### PLOT GLM.NB IPS vs climate veg_tmp -----------------------
# Generate a sequence of values for each predictor
dat_veg_tmp <- data.frame(veg_tmp = seq(min(dd_simpl$veg_tmp), max(dd_simpl$veg_tmp), length.out = 100),
                          previous_spei3_2 = mean(dd_simpl$previous_spei3_2, na.rm = TRUE))

# Predictions for veg_tmp
dat_veg_tmp$predicted <- predict(m9, newdata = dat_veg_tmp, type = "response")
ci_veg_tmp <- predict(m9, newdata = dat_veg_tmp, type = "response", se.fit = TRUE)
dat_veg_tmp$upper <- ci_veg_tmp$fit + 1.96 * ci_veg_tmp$se.fit
dat_veg_tmp$lower <- ci_veg_tmp$fit - 1.96 * ci_veg_tmp$se.fit


# Prepare data for previous_spei3_2
dat_previous_spei3_2 <- data.frame(previous_spei3_2 = seq(min(dd_simpl$previous_spei3_2), 
                                                          max(dd_simpl$previous_spei3_2), length.out = 100),
                                   veg_tmp = mean(dd_simpl$veg_tmp, na.rm = TRUE))

# Predictions for previous_spei3_2
dat_previous_spei3_2$predicted <- predict(m9, newdata = dat_previous_spei3_2, type = "response")
ci_previous_spei3_2         <- predict(m9, newdata = dat_previous_spei3_2, type = "response", se.fit = TRUE)
dat_previous_spei3_2$upper <- ci_previous_spei3_2$fit + 1.96 * ci_previous_spei3_2$se.fit
dat_previous_spei3_2$lower <- ci_previous_spei3_2$fit - 1.96 * ci_previous_spei3_2$se.fit

# Prepare data for previous_spei12_2
dat_previous_spei12_2 <- data.frame(previous_spei12_2 = seq(min(dd_simpl$previous_spei12_2), 
                                                            max(dd_simpl$previous_spei12_2), length.out = 100),
                                    veg_tmp = mean(dd_simpl$veg_tmp, na.rm = TRUE),
                                    previous_spei3_2 = mean(dd_simpl$previous_spei3_2, na.rm = TRUE))

# Predictions for previous_spei12_2
dat_previous_spei12_2$predicted <- predict(m9, newdata = dat_previous_spei12_2, type = "response")
ci_previous_spei12_2 <- predict(m9, newdata = dat_previous_spei12_2, type = "response", se.fit = TRUE)
dat_previous_spei12_2$upper <- ci_previous_spei12_2$fit + 1.96 * ci_previous_spei12_2$se.fit
dat_previous_spei12_2$lower <- ci_previous_spei12_2$fit - 1.96 * ci_previous_spei12_2$se.fit



# Define a common theme
common_theme <- theme_minimal(base_size = 11) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(color = "black"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5),
        aspect.ratio = 1 )

# Plot for veg_tmp with orange color
p1 <- ggplot(dat_veg_tmp, aes(x = veg_tmp, y = predicted)) +
  geom_line(color = "#E69F00") +
  ylim(10000, 35000) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#E69F00", alpha = 0.2) +
  labs(title = "", x = "Temperature (veg_tmp)", y = "Sum Beetles/year") +
  common_theme



# Plot for previous_spei3_2 with blue color
p2 <- ggplot(dat_previous_spei3_2, aes(x = previous_spei3_2, y = predicted)) +
  geom_line(color = "#0072B2") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
  ylim(10000, 35000) +
  labs(title = "", x = "SPEI (previous_spei3_2)", y = "Sum Beetles/year") +
  common_theme

# Plot for previous_spei12_2 with green color
p3 <- ggplot(dat_previous_spei12_2, aes(x = previous_spei12_2, y = predicted)) +
  geom_line(color = "#009E73") +
  ylim(10000, 35000) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#009E73", alpha = 0.2) +
  labs(title = "", x = "SPEI (previous_spei12_2)", y = "Sum Beetles/year") +
  common_theme


ggarrange(p1, p2, p3, ncol = 3, nrow = 1)


require(ggiraph)
require(ggiraphExtra)
ggPredict(m2,se=TRUE,interactive=TRUE,digits=3)





summary(m2)


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

dd <- dat_lag_scaled %>% 
  dplyr::select(c(agg_doy, 
                  peak_doy,
                  spring_tmp,
                  veg_tmp,
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
dd$tr_agg_doy  <- (dd$agg_doy - 60) / (300 - 60)
dd$tr_peak_doy <- (dd$peak_doy - 60) / (300 - 60)

# Fit a beta regression model

library(glmmTMB)

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dd$tr_agg_doy  <- pmin(pmax(dd$tr_agg_doy, 1e-4),  1 - 1e-4)
dd$tr_peak_doy <- pmin(pmax(dd$tr_peak_doy, 1e-4), 1 - 1e-4)




# Transform your response variable to fit within 0 to 1
dd_simpl$tr_agg_doy  <- (dd_simpl$agg_doy - 60) / (300 - 60)
dd_simpl$tr_peak_doy <- (dd_simpl$peak_doy - 60) / (300 - 60)

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dd_simpl$tr_agg_doy  <- pmin(pmax(dd_simpl$tr_agg_doy, 1e-4),  1 - 1e-4)
dd_simpl$tr_peak_doy <- pmin(pmax(dd_simpl$tr_peak_doy, 1e-4), 1 - 1e-4)


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




summary(m.agg02)
AIC(m.agg1, m.agg2, m.agg3, m.agg4, m.agg5, m.agg6, m.agg7, m.agg8, m.agg9, m.agg10, m.agg11,m.agg11_1, m.agg12)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.agg12.3,#m.agg12.3,  #m.agg01 
                                              plot = T)
AICc( m.agg01,m.agg13, m.agg14, m.agg15, m.agg16, m.agg12, m.agg12.2, m.agg12.3, m.agg12.4,m.agg12.5, m.agg12.3.int, m.agg12.3.12, m.agg0)

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




# what is whatL example for agg and peak:

dat_lag_scaled_complete %>% 
  dplyr::select(c(agg_doy, tr_agg_doy, peak_doy, tr_peak_doy )) %>% 
  arrange(agg_doy) %>% 
  View()







