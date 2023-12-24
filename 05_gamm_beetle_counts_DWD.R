


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





##### read data ----------------------------------------------------------------

load(file=   "outData/final_table.Rdata")
load(file =  "outData/buffers.Rdata")  # df_RS_out
load(file =  "outData/lisa.Rdata")     # read LISA Moran's I stats

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





##### prepare data RS --------------------------------------------
df_RS_out <- df_RS_out %>% 
  dplyr::select(-c(globalid, id)) %>% 
  dplyr::rename(trapID = falsto_name)

# change column names
lisa_merged_df <- lisa_merged_df %>% 
  dplyr::rename(trapID = falsto_name,
                sum_ips = sum_beetle) %>% 
  dplyr::select(c(year, trapID, sum_ips, Morans_I))

# add lag: previous year counts, previous year temperature
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
                  previous_spring_veg_tmp =  dplyr::lag(spring_veg_tmp, order_by = year),
                  previous_veg_veg_prcp   =  dplyr::lag(veg_veg_prcp,   order_by = year),
                  previous_spei1      =  dplyr::lag(spei1,      order_by = year),
                  previous_spei3      =  dplyr::lag(spei3,      order_by = year),
                  previous_spei12     =  dplyr::lag(spei12,     order_by = year),
                  previous_spei24     =  dplyr::lag(spei24,     order_by = year)) %>% 
  dplyr::mutate(population_growth     = (sum_ips - previous_sum_ips) / previous_sum_ips * 100,
                population_growth2    = dplyr::lag(population_growth, order_by = year)) %>%  # lag population growth by one more year
  left_join(df_RS_out, 
            by = join_by(trapID, year, sum_ips)) %>%
  left_join(lisa_merged_df, 
            by = join_by(trapID, year, sum_ips)) %>%
  dplyr::mutate(previous_Moran       =  dplyr::lag(Morans_I , order_by = year),
                previous_Moran2      =  dplyr::lag(Morans_I , order_by = year, n = 2),
                ) %>% 
  dplyr::mutate(trapID = as.factor(trapID)) #%>% 


# check if results are correct
dat_lag %>%
  filter(trapID== 'Zusmarshausen_1') %>%
  dplyr::select(c(trapID, year, sum_ips, previous_sum_ips,previous_sum_ips2,
                  veg_tmp,
                  previous_Moran,previous_Moran2,
                  Morans_I,
                  population_growth, population_growth2,
                  previous_veg_tmp)) #%>%
  #View()


## Scale predictors ================================================================================


dat_lag_scaled <- dat_lag #%>%
 
# skip columns if not for scaling
skip_col <- c('trapID', 'pairID', "x", "y", 'year', 'spei1','spei3','spei12','spei24','previous_spei1', 'previous_spei3', 'previous_spei12', 'previous_spei24' )

dat_lag_scaled <-
  dat_lag %>%
  ungroup(.) %>% 
  # Select columns to scale, skipping trapID and year
  dplyr::select(-all_of(skip_col )) %>%
  # Apply the scale function
  scale(center = TRUE, scale = TRUE) %>%
  # Convert back to a data frame
  as.data.frame() %>%
  # Bind the unscaled columns back
  bind_cols(dat_lag %>% dplyr::select(all_of(skip_col)), .)

#simplify, keep all columns as vectors 
dat_lag_scaled <- dat_lag_scaled %>% 
  ungroup(.) %>% 
  dplyr::mutate(across(where(~is.numeric(.) && is.matrix(.)), as.vector))

# consider year as random factor - non linear relationship betweeen years

# remove additionsl NAs
dat_lag_scaled_complete <- dat_lag_scaled %>%
  dplyr::select(-anom_wind_beetle) %>% 
  na.omit()
#  filter(across(everything(), ~ !is.na(.)))



# spatial synchronization with veg_tmp? ---------------------------------------

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

# RS:  Link beetle data with observed RS damage ------------------------------------------------------------------------
# remove 0 wind-beetle damage 
# run analyses with only stable traps over time?? - split df



# does beetle sums corresponds to beetle anomaly? YES! if I am using only smooths, not points... ---------


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



# get glm RS vs ips sums: DREDGE  ----------------------------------------------------------------------
# investigate all variables & dredge

# keep only complete cases - can alter this to keep more rows? can remove unnecessary columns
# remove unnecessare predictors, also lag values only by 1 year
#col_remove <- c( "x",  "y", "sm" ,"vpd", "sm_z","vpd_z", "previous_sum_ips2","previous_agg2", "population_growth2" )

# remove extra columns !!!! - not necessary to remove those columns, as those are complete!
#dat_lag_scaled <- dat_lag_scaled %>%
#  dplyr::select(-all_of(col_remove)) #%>% 
 # filter_all(all_vars(!is.na(.))) %>% 
#  filter(across(everything(), ~ !is.na(.)))

### get RS with beetle population predictors: ----------------------------------
# Fit a global model with all beetle dynamic predictors---------------------------------
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

## Manually create models -------------------------------------------------------

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






# RS harvest =================================================================
### get RS with beetle population predictors: ----------------------------------
# Fit a global model with all beetle dynamic predictors---------------------------------
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



# RS explore both climatic and beetles  variables -----------------------------------------------------



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




# examine best models from dredge -----------------------------------------------

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

pairs(sum_ips ~ veg_tmp + veg_prcp + spei3 + vpd, dat_lag, panel = panel.smooth)

pairs(sum_ips ~ veg_tmp + spei3 + vpd + elev + spruce_1986 + previous_sum_ips+population_growth, dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+population_growth + veg_tmp + previous_veg_tmp , dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+ population_growth + veg_tmp + previous_veg_tmp + previous_spei3 , dat_lag, panel = panel.smooth)

pairs(sum_ips ~ previous_sum_ips+ population_growth + veg_tmp + previous_veg_tmp + population_growth2, dat_lag, panel = panel.smooth)

hist(dat_lag$population_growth)


# make glm - beetle sums by temp and spei
# test with the raw and scaled predictors;
# glm, - define the polynomial relationship
### expore plots ------------------------------------------------------------------------
p.veg_tmp <- dat_lag_scaled %>% 
ggplot(aes(x = veg_tmp,
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


ggarrange(p.sm, p.veg_tmp, p.vpd, nrow = 1)

## test GLM (no random effects) on raw data ---------------------------------------------------------------------
m1.poly <- glm(sum_ips ~ I(veg_tmp^15), dat_lag, family = 'poisson')
m1.exp <- glm(sum_ips ~ log(veg_tmp), dat_lag, family = 'poisson')



# see in plot
# New data for prediction
newdata_poly <- data.frame(veg_tmp = seq(min(dat_lag$veg_tmp), max(dat_lag$veg_tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.poly, newdata = newdata_poly, type = "response")


ggplot(dat_lag, aes(x = veg_tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_poly, aes(x = veg_tmp, y = predicted), color = "red",  lwd = 1.5) +
  geom_line(data = newdata_exp, aes(x = veg_tmp, y = predicted), color = "green", lwd = 1.5) +
  geom_smooth()



# New data for prediction: exponential
newdata_exp <- data.frame(veg_tmp = seq(min(dat_lag$veg_tmp), max(dat_lag$veg_tmp), length.out = 100))

# Generate predictions
newdata_exp$predicted <- predict(m1.exp, newdata = newdata_exp, type = "response")
ggplot(dat_lag, aes(x = veg_tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_exp, aes(x = veg_tmp, y = predicted), color = "red")


## Fit the GLMM with negative binomial distribution & random effects ------------
m1.exp <- glm.nb(sum_ips ~ veg_tmp, data = dat_lag, link = log)
newdata_exp <- data.frame(veg_tmp = seq(min(dat_lag$veg_tmp), max(dat_lag$veg_tmp), length.out = 100))
newdata_exp$predicted <- predict(m1.exp, newdata = newdata_exp, type = "response")

ggplot(dat_lag, aes(x = veg_tmp, y = sum_ips)) +
  geom_point() +
  geom_line(data = newdata_exp, aes(x = veg_tmp, y = predicted), color = "red")


m1.poly2 <- glm.nb(sum_ips ~ poly(veg_tmp, 2), data = dat_lag)
newdata_poly <- data.frame(veg_tmp = seq(min(dat_lag$veg_tmp), max(dat_lag$veg_tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.poly, newdata = newdata_poly, type = "response")


m1.poly3 <- glm.nb(sum_ips ~ poly(veg_tmp, 3), data = dat_lag)
m1.sc.poly2 <- glm.nb(adjusted_variable  ~ poly(veg_tmp, 2), data = dat_lag_scaled)
newdata_poly <- data.frame(veg_tmp = seq(min(dat_lag_scaled$veg_tmp), max(dat_lag_scaled$veg_tmp), length.out = 100))
newdata_poly$predicted <- predict(m1.sc.poly2, newdata = newdata_poly, type = "response")


ggplot(dat_lag_scaled, aes(x = veg_tmp, y = adjusted_variable)) +
  geom_point() +
  geom_line(data = newdata_poly, aes(x = veg_tmp, y = predicted), color = "blue")


AIC(m1, m1.exp, m1.poly2, m1.poly3)
simulationOutput <- simulateResiduals(fittedModel = m1.sc.poly2, plot = T)
testOutliers(m1.sc.poly2, type = 'bootstrap') 


m2 <- glm.nb(sum_ips ~ veg_tmp + spei, dat_lag)
m3 <- glm.nb(sum_ips ~ veg_tmp + veg_tmp_z, dat_lag)
m4 <- glm.nb(sum_ips ~ I(veg_tmp) + previous_sum_ips   , dat_lag)
summary(m4)
plot(m4, page = 1)




# DREDGE Variable selection, IPS vs environment  ------------------------------------------------------

# > (vif_model_ips <- car::vif(global_model))
# GVIF Df GVIF^(1/(2*Df))
# poly(veg_tmp, 2)          20.319964  2        2.123150
# veg_prcp                   4.010201  1        2.002549
# spei                   2.112911  1        1.453586
# sm                     1.463907  1        1.209920
# vpd                   16.692437  1        4.085638
# previous_sum_ips  1.092418  1        1.045188
# spruce_1986           13.205598  1        3.633951
# elev                   4.470761  1        2.114417
# remained_spruce       12.135460  1        3.483599

test <- lm(elev ~ poly(veg_tmp,3), dat_lag_scaled_complete)
plot(allEffects(test))

plot(elev ~ veg_tmp, dat_lag_scaled_complete)

t1 <- lm(elev ~ veg_tmp, dat_lag_scaled_complete)
t2 <- lm(elev ~ poly(veg_tmp,2), dat_lag_scaled_complete)
t3 <- lm(elev ~ poly(veg_tmp,3), dat_lag_scaled_complete)


# Adjusted R-squared
summary(t1)$adj.r.squared
summary(t2)$adj.r.squared
summary(t3)$adj.r.squared

# AIC and BIC
AIC(t1, t2, t3)
BIC(t1, t2, t3)

# F-test for model comparison
anova(t1, t2, t3)




# keep only complete cases 
dat_lag_scaled_complete <- dat_lag_scaled %>%
  na.omit()
  #filter_all(all_vars(!is.na(.)))
  


# Fit a global model with all potential predictors
global_model <- glmer(sum_ips ~ poly(veg_tmp,2) + veg_prcp + spei + sm + 
                         #vpd + 
                         previous_sum_ips + 
                         #spruce_1986 + 
                         elev + remained_spruce + (1 | pairID/trapID), 
                      data = dat_lag_scaled_complete,
                      na.action = 'na.fail',
                    
                      family = negative.binomial(2.8588)
                      #, na.action = "na.omit"  "na.fail"
                      )

# get a vector of vif values
(vif_model_ips <- car::vif(global_model))


# Run dredge on the global model
model_set <- dredge(global_model)

# View the results
print(model_set)

# Sort by AIC
model_set_sorted <- model_set[order(model_set$AIC),]

# View the top models with the lowest AIC
head(model_set_sorted)


# examine models:

# Assuming you want to examine model 207 from your model set
i = 6
selected_model <- get.models(model_set, subset = i)  # subset = 1 for the best model; adjust if needed

# Get a summary of the selected model
summary(selected_model)  # selected_model is a list, [[1]] accesses the first element

simulationOutput <- simulateResiduals(fittedModel = selected_model  , plot = T) # m3.7.1.aut2

testOutliers(selected_model , type = 'bootstrap') # i have outliers dectected, need to simplify the model! 
# Calculate pseudo-R^2 for the selected model
r.squared <- r.squaredGLMM(selected_model[[i]])
print(r.squared)  # prints both marginal and conditional R^2; for GLM, marginal R^2 is more relevant

# my best variables are:
#(Int)        elv ply(veg_tmp,2)      prc prv_sum_ips_scl rmn_spr       sm     spe df    logLik     AICc
#95  9.882                     + -0.01892          0.6325 0.01873          0.05847  8 -51722.85 103461.7
# does not really good, the models fit is very poor and variables are not meaningfull


# Plot the effects
windows()
plot(selected_model[[1]], page = 1)

# Get the coefficients
coef(selected_model[[1]])

plot(allEffects(selected_model[[1]]))


### GLMM - ips counts with random effect and nb ----------------------------------------


m3 <- glmer.nb(sum_ips ~ veg_tmp + spei + (1|pairID) + (1|trapID), data = dat_lag)

m3.1 <- glmer.nb(sum_ips ~ veg_tmp  + (1|trapID), data = dat_lag)

m3.2 <- glmer.nb(sum_ips ~ veg_tmp  + (1|pairID), data = dat_lag)
m3.3 <- glmer.nb(sum_ips ~ veg_tmp + (1| year), data = dat_lag)

m3.4 <- glmer.nb(sum_ips ~ veg_tmp + year + (1| pairID), data = dat_lag)

m3.5 <- glmer.nb(sum_ips ~ veg_tmp + (1| year), data = dat_lag_scaled)
# consider as polynomial
m3.5.1 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + (1| year), data = dat_lag_scaled)

# add locatin as random effect
m3.5.2 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + (1| year) + (1| pairID), data = dat_lag_scaled)


# simplify year effect" replace by drought ref
m3.5.3 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + (1| temp_fact) + (1| pairID), data = dat_lag_scaled)


summary(m3.5.2)
AICc(m3.5.2, m3.5.3, m3.7)  # m3.7 is better

m3.6 <- glmer.nb(sum_ips ~ veg_tmp + (1| year) + (1| pairID), data = dat_lag_scaled)

m3.7 <- glmer.nb(sum_ips ~ veg_tmp + (1| year) + (1| pairID/trapID), data = dat_lag_scaled)

# add poly function, year is consider as random to account for diversity over years
m3.7.1 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + (1| year) + (1| pairID/trapID), data = dat_lag_scaled)


# account for temp autocorrelation: use lagged values to capture it
m3.7.1.aut <- glmer(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + (1 | year) + (1 | pairID/trapID), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))


#add precipitaion
m3.7.1.aut2 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp ,  #+ (1 | pairID/trapID)
                    data = dat_lag_scaled)
# add other enviro predistors: remained spruce
m3.7.1.aut3 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce,  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled)

# add spei
m3.7.1.aut4 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce + spei,  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled)

# add spei as poly
m3.7.1.aut5 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce + poly(spei,2),  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled)

# add previous temp
m3.7.1.aut6 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce + previous_veg_tmp + poly(spei,2),  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled_complete)

m3.7.1.aut6.1 <- glm.nb(sum_ips ~ previous_sum_ips + poly(previous_veg_tmp,2) + poly(spei,2),  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled_complete)


# remove previous temp as highly correlared
m3.7.1.aut6.2 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce + poly(spei,2),  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled_complete)


m3.7.1.aut6.3 <- glm.nb(sum_ips ~  previous_sum_ips + veg_prcp + remained_spruce + previous_veg_tmp + poly(spei,2),  #+ (1 | pairID/trapID)
                        data = dat_lag_scaled_complete)

# simplify model by removing non-significant predictors
m3.7.1.aut6.4 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          #veg_prcp + 
                          #remained_spruce + 
                          previous_veg_tmp + spei,  
                        data = dat_lag_scaled_complete)

# get previous veg_tmp as poly
m3.7.1.aut6.5 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          #veg_prcp + 
                          #remained_spruce + 
                          poly(previous_veg_tmp,2) + spei,  
                        data = dat_lag_scaled_complete)


# get spei as poly
m3.7.1.aut6.6 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          #veg_prcp + 
                          #remained_spruce + 
                          poly(previous_veg_tmp,2) + poly(spei,2),  
                        data = dat_lag_scaled_complete)


# add interaction - not improve the terms
m3.7.1.aut6.7 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          #veg_prcp + 
                          #remained_spruce + 
                          poly(previous_veg_tmp,2) + poly(spei,2) +
                          previous_sum_ips:veg_tmp,  
                        data = dat_lag_scaled_complete)



# add interaction - not improve the terms
m3.7.1.aut6.8 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          veg_tmp +
                          #veg_prcp + 
                          #remained_spruce + 
                          poly(previous_veg_tmp,2) + poly(spei,2),  
                        data = dat_lag_scaled_complete)


m3.7.1.aut6.9 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          poly(veg_tmp,2) +
                          #veg_prcp + 
                          #remained_spruce + 
                          poly(previous_veg_tmp,2) + spei,  
                        data = dat_lag_scaled_complete)

m3.7.1.aut6.10 <- glm.nb(sum_ips ~  previous_sum_ips + 
                          poly(veg_tmp,2) +
                          #veg_prcp + 
                          #remained_spruce + 
                          previous_veg_tmp + spei,  
                        data = dat_lag_scaled_complete)


AIC(m3.7.1.aut6.8, m3.7.1.aut6.9, m3.7.1.aut6.4)

# add previous temp as poly
m3.7.1.aut7 <- glm.nb(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + veg_prcp + remained_spruce + poly(previous_veg_tmp,2) + poly(spei,2),  #+ (1 | pairID/trapID)
                      data = dat_lag_scaled_complete)


AIC(m3.7.1.aut6, m3.7.1.aut6.3, m3.7.1.aut6.4, m3.7.1.aut6.5)  # the best: m3.7.1.aut6
plot(allEffects(m3.7.1.aut6.1))

summary(m3.7.1.aut6.1)
# r2(m3.7.1.aut6.10)

dat_lag_scaled_complete %>% 
  ggplot(aes(x = veg_tmp,
             y = sum_ips)) +
  geom_smooth()



dat_lag %>% 
  ggplot(aes(x = spei,
             y = sum_ips)) +
  geom_smooth()

AICc(m3.7.1.aut2, m3.7.1.aut3, m3.7.1.aut4, m3.7.1.aut5, m3.7.1.aut6)




# scale previous sum ips
m3.7.1.aut.sc1 <- glmer(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips + (1 | year) + (1 | pairID/trapID), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))

# remove 'year' as random
m3.7.1.aut.sc2 <- glmer(sum_ips ~ poly(veg_tmp, 2) + previous_sum_ips+ (1 | pairID/trapID), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))


# scale the lagged predictor; previous beetle sums
m3.7.1.aut2 <- glmer(sum_ips ~ veg_tmp + previous_sum_ips + (1 | year) + (1 | pairID/trapID), 
                    data = dat_lag_scaled, 
                    family = negative.binomial(2.8588))

# year as numeric, random
m3.7.1.aut2.1 <- glmer(sum_ips ~ veg_tmp + previous_sum_ips + (1 | year_num) + (1 | pairID/trapID), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))

# remove 'year' as random effect
m3.7.1.aut3 <- glmer(sum_ips ~ veg_tmp + previous_sum_ips +  (1 | pairID/trapID), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))

# add year as nested random effect
m3.7.1.aut4 <- glmer(sum_ips ~ veg_tmp + previous_sum_ips +  (1 | pairID/trapID/year), 
                     data = dat_lag_scaled, 
                     family = negative.binomial(2.8588))


# add poly function, year is consider as random to account for diversity over years
m3.7.1.quadr <- glmer.nb(sum_ips ~ I(veg_tmp^2) + (1| year) + (1| pairID/trapID), data = dat_lag_scaled)


# interaction between veg_tmp and beetle counts? 
m_counts1 <- glm.nb(sum_ips ~ poly(veg_tmp,2) + previous_sum_ips, dat_lag_scaled )
m_counts2 <- glm.nb(sum_ips ~ veg_tmp + previous_sum_ips + previous_sum_ips:veg_tmp, dat_lag_scaled )

plot(allEffects(m_counts1))


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

acf(residuals(m3.7.1.aut6))  # seems relatively ok
pacf(residuals(m3.7.1.aut6))


cor(dat_lag_scaled_complete$veg_tmp, dat_lag_scaled_complete$previous_veg_tmp)
cor(dat_lag_scaled_complete$spei, dat_lag_scaled_complete$veg_tmp)



plot(allEffects(m3.7.1.aut6))
dw_result <- dwtest(residuals(m3.7.1.aut6) ~ fitted(m3.7.1.aut6))
(dw_result)

summary(m3.7.1.aut6)
vif(m3.7.1.aut6)

r2(m3.7.1.aut6)
# Durbin-Watson test
# 
# data:  residuals(m3.7.1.aut2) ~ fitted(m3.7.1.aut2)
# DW = 1.6613, p-value = 8.204e-08
# alternative hypothesis: true autocorrelation is greater than 0

simulationOutput <- simulateResiduals(fittedModel = m3.7.1.aut2, plot = T)
simulationOutput <- simulateResiduals(fittedModel = m3.7.1.aut, plot = T) # this one is little bit better from output

# goodness of fit!
r2(m3.7.1.aut6.4) # Nagelkerke s R2

# account for autocorrelation: use GAM and MER model - does not work!!
library(gamm4)
gamm4_model <- gamm4(sum_ips ~ poly(veg_tmp, 2), 
                     random = ~(1 | year_num) + (1 | pairID/trapID), 
                     data = dat_lag_scaled, 
                     #family = nb(2.8588),
                     correlation = corAR1(form = ~ 1 | year_num))


# simplify years by heat cetegories
m3.7.2 <- glmer.nb(sum_ips ~ poly(veg_tmp,2) + (1| temp_fact) + (1| pairID/trapID), data = dat_lag_scaled)

# add interaction heat vs ips_sums
m3.7.3.poly <- glmer.nb(sum_ips ~ poly(veg_tmp,2)*temp_fact + (1| pairID/trapID), data = dat_lag_scaled)

# keep it linear
m3.7.4.lin <- glmer.nb(sum_ips ~ veg_tmp*temp_fact + (1| pairID/trapID), data = dat_lag_scaled)

# keep it exp
m3.7.4.exp <- glmer.nb(sum_ips ~ exp(veg_tmp)*temp_fact + (1| pairID/trapID), data = dat_lag_scaled)


AICc(m3.7.1, m3.7.3.poly, m3.7.4.lin, m3.7.4.exp)

# try as smooths??
m3.7.3.smooth <- glmer.nb(sum_ips ~ splines::bs(veg_tmp, df = 4)*temp_fact + (1| pairID/trapID), data = dat_lag_scaled)



#splines::bs(x, degrees = 3) * factor(category)

# groups years by factor: temp factor by years???

AICc(m3.7, m3.7.1, m3.7.2) # the best: m3.7.1

m3.8 <- glmer.nb(sum_ips ~ veg_tmp  + (1| year) + (1| pairID/trapID), data = dat_lag_scaled)

m3.5.1 <- glmer.nb(sum_ips ~ veg_tmp  + (1| year), data = dat_lag_scaled)


summary(m3.7.2)
AICc(m3.5, m3.6, m3.7)

r2(m3.7.2)

simulationOutput <- simulateResiduals(fittedModel = m3.7.3.poly, plot = T)


# Predicting beetle counts
AIC(m3.5, m3.7.1, m3.7.3.poly)

### plot GLM.nb  ===================================================
unique_trapIDs <- levels(dat_lag_scaled$trapID)
unique_pairID   <- levels(dat_lag_scaled$pairID)
unique_temp       <- levels(dat_lag_scaled$temp_fact)
temp_range        <- seq(min(dat_lag_scaled$veg_tmp), 
                         max(dat_lag_scaled$veg_tmp), 
                         length.out = 100)
new_data <- expand.grid(veg_tmp = temp_range, 
                        temp_fact = levels(dat_lag_scaled$temp_fact),
                        pairID = unique_pairID[1],
                        trapID = unique_trapIDs[1])

new_data$predicted_counts <- predict(m3.7.4.exp, newdata = new_data, re.form = NULL)

# Plotting
ggplot(new_data, aes(x = veg_tmp, y = predicted_counts, color = temp_fact)) +
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
ggplot(dat_lag_scaled, aes(x = veg_tmp, y = sum_ips, color = temp_fact)) +
  geom_point() +
  geom_smooth()+
  theme_minimal()



m4 <- glmer.nb(sum_ips ~ veg_tmp + spei + (1|pairID) + (1|trapID), 
               data = dat_lag_scaled)

m5 <- glmer.nb(sum_ips ~ veg_tmp + spei + (1 | pairID/trapID), # account for nested design
               data = dat_lag_scaled)

m6 <- glmer.nb(sum_ips ~ veg_tmp + (1 | pairID/trapID), # account for nested design
               data = dat_lag_scaled)

m7 <- glmer.nb(sum_ips ~ veg_tmp + elev + (1 | pairID/trapID), # account for nested design
               data = dat_lag_scaled)  # error

m8 <- glmer.nb(sum_ips ~ veg_tmp + spruce_1986 + (1 | pairID/trapID), # account for nested design
               data = dat_lag_scaled)  

# account for nteraction between temp and year account for spatial and temporal autocorrelation
m9 <- glmer.nb(sum_ips ~ veg_tmp * year + spruce_1986 + (1 | pairID/trapID), 
               data = dat_lag_scaled)

# visualize teh effect of temperature over years
# Assuming 'veg_tmp' ranges from min_veg_tmp to max_veg_tmp in your data
veg_tmp_range <- seq(min(dat_lag_scaled$veg_tmp), 
                 max(dat_lag_scaled$veg_tmp), 
                 length.out = 100)

# Create a data frame for predictions
new_data <- expand.grid(veg_tmp = veg_tmp_range, 
                        year = levels(dat_lag_scaled$year), 
                        spruce_1986 = mean(dat_lag_scaled$spruce_1986))  # Adjust other variables as needed

# Predict using the model
new_data$predicted <- predict(m9, newdata = new_data, re.form = NA, type = "response")

# Plotting
ggplot(new_data, 
       aes(x = veg_tmp, y = predicted, color = year)) +
  geom_line() +
  labs(title = "Effect of Temperature across Years",
       x = "Temperature",
       y = "Predicted sum_ips") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")



# allow temparature interacting with year, but not year as individual preddictor
m10 <- glmer.nb(sum_ips ~ veg_tmp:year + veg_tmp + spruce_1986 + (1 | pairID/trapID), 
                data = dat_lag_scaled)



### glmmTMB account for temporal correlation and overdispersion in data -------------


# not meaningfull!!!
m11 <- glmmTMB(sum_ips ~ veg_tmp + spruce_1986 + 
                (1 | pairID/trapID) + 
                (0 + year | pairID/trapID) + 
                 ar1(year + 0|trapID), # account for autocorrelation, as here: https://stats.stackexchange.com/questions/393734/covariance-structures-in-glmmtmb-for-temporal-autocorrelation
              data = dat_lag_scaled,
              family = nbinom2,
              ziformula = ~ 1,  # Zero-inflation formula, adjust as necessary
              dispformula = ~ 1#,  # Dispersion formula, adjust as necessary
              #corStruct = corAR1(form = ~ year | pairID/trapID)
              )

# m11 does not work: as the glmmTMB does not work with the temporal autocorrelation
# create the year as random intercept to account for systematic error in year:

m12 <- glmmTMB(sum_ips ~ veg_tmp + spruce_1986 + (1 + year | pairID/trapID),
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
             x = veg_tmp,
             color = factor(year))) +
  geom_smooth()


dat_lag %>% 
  ggplot(aes(y = sum_ips,
             x = veg_tmp,
             color = factor(temp_fact))) +
  geom_point(alpha = 0.3)+
  geom_smooth() +
  facet_grid(.~temp_fact) +
  theme_bw()

# decide about which years are warm and which not??
# check the z-score
dat_lag %>% 
  ggplot(aes(y = veg_tmp_z,
             x = year,
             group = year)) +
  geom_boxplot(notch = T)


m.gam0 <- bam(sum_ips ~ 
                s(veg_tmp, k = 4),# +
              #s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
              #s(pairID, bs = 're') +   # trap pair
              #s(x,y, k = 50)
              # ,
              dat_lag_scaled, 
              family = tw)


m.gam1 <- bam(sum_ips ~ 
                s(veg_tmp, by=temp_fact, k = 4),# +
            #s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            #s(pairID, bs = 're') +   # trap pair
            #s(x,y, k = 50)
           # ,
          dat_lag_scaled, 
          family = tw)


m.gam2 <- bam(sum_ips ~ 
                s(veg_tmp, by=temp_fact, k = 4) +
              s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
              s(pairID, bs = 're') +   # trap pair
              s(x,y),
              dat_lag_scaled, 
              family = tw)

# add base 'ds' for coordinates
m.gam3 <- bam(sum_ips ~ 
                s(veg_tmp, by=temp_fact, k = 4) +
                s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(pairID, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)

# add base 'ds' for coordinates
m.gam4 <- bam(sum_ips ~ 
                s(veg_tmp, by=temp_fact, k = 4) +
                s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(pairID, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)

# add year variable to account for autocorrelation
m.gam5 <- bam(sum_ips ~ 
                s(veg_tmp, by=temp_fact, k = 4) +
                s(year,bs = 're', k = 3) +
                s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(pairID, bs = 're') +   # trap pair
                s(x,y, bs = 'ds'),
              dat_lag_scaled, 
              family = tw)


# add lagged variables to account for temp autocorrelation
m.gam6 <- bam(sum_ips ~ 
                s(previous_sum_ips, k = 3) +
                s(veg_tmp, by=temp_fact, k = 4) +
                s(year,bs = 're', k = 3) +
                s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(pairID, bs = 're') +   # trap pair
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
m.gamm <- gamm(sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
                 s(year_num, by = temp_fact, bs = "re") +
                 s(pairID, bs = "re") +
                 s(trapID, bs = "re"),
               data = dat_lag_scaled, 
               family = nb())
appraise(m.gamm$gam)
plot(m.gamm$gam, page = 1)

acf(residuals(m.gamm$gam))
pacf(residuals(m.gamm$gam))

# 1. if data are overdispersed, need to use neg bin family

m.glmmTMB <- glmmTMB(sum_ips ~ poly(veg_tmp, 2) * temp_fact + 
                       (1 | year_num:temp_fact) + 
                       (1 | pairID) + 
                       (1 | trapID),
                     data = dat_lag_scaled, 
                     family = nbinom2)

# 2. handle overdispersion by add observation level ID as random effect
m.gam1 <- gam(sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
               s(pairID, bs = "re") +
               s(trapID, bs = "re") +
               s(id, bs = "re"),  # Observation-level random effect
             data = dat_lag_scaled, 
             family = nb())

m.gam2 <- gam(sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
                s(pairID, bs = "re") +
                s(trapID, bs = "re") +
                s(id, bs = "re"),  # Observation-level random effect
              data = dat_lag_scaled, 
              family = poisson())

# the best!
m.gam3 <- gam(sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
                s(pairID, bs = "re") +
                s(trapID, bs = "re") +
                s(id, bs = "re"),  # Observation-level random effect
              data = dat_lag_scaled, 
              family = tw())


  plot_smooth(m.gam3, 
                          view="veg_tmp", 
                          rm.ranef =TRUE, 
                          #transform=function(x) 100*inv.logit(x),  
                          plot_all="temp_fact")$fv




###  GAM shown in a plot --------------------------------------------------------

# Select specific levels for id, pairID, and trapID
selected_id <- unique(dat_lag_scaled$id)[5]      # Replace 1 with your chosen level
selected_pairID <- unique(dat_lag_scaled$pairID)[5]
selected_trap <- unique(dat_lag_scaled$trapID)[5]

# Generate Predicted Values
new_data <- expand.grid(temp_fact = unique(dat_lag_scaled$temp_fact),
                        veg_tmp = seq(min(dat_lag_scaled$veg_tmp), max(dat_lag_scaled$veg_tmp), length.out = 100),
                        pairID = selected_pairID,
                        trapID = selected_trap,
                        id = selected_id)

# Predict using your model
new_data$predicted <- predict(m.gam3, newdata = new_data, type = "response")

# Create the Plot
ggplot(new_data, aes(x = veg_tmp, y = predicted, color = temp_fact)) +
  geom_line() +
  labs(title = "Predicted sum_ips by Temp Factor",
       x = "Temperature",
       y = "Predicted sum_ips") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

## GAMM test model run, before improving autocorrelation - simpler - does not work!! -----------------
m.gamm1 <- gamm(
  sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
    s(pairID, bs = "re") +
    s(trapID, bs = "re") +
    s(id, bs = "re"),
  data = dat_lag_scaled,
  family = nb()
)

# include autocorrelation; make sure that 
m.gamm2 <- gamm(
  sum_ips ~ s(veg_tmp, by = temp_fact, k = 4) +
    s(pairID, bs = "re") +
    s(trapID, bs = "re") +
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
  sum_ips ~ s(veg_tmp, k = 4),# +
   # s(pairID, bs = "re") +
  #  s(trapID, bs = "re") +
   # s(id, bs = "re"),
  data = dat_lag_scaled,
  #family = nb(),
  correlation = corAR1(form = ~ year_num)
)



# using mixed effects GAMM with explicit temporal autocorrelation -------------
# Example of using gamm for temporal autocorrelation

# to account for autocorrelation using corAR, I need to have one year per category
m.gamm1 <- gamm(sum_ips ~ s(veg_tmp, by=year, k = 4), 
                  s(pairID, bs = "re") +
                  s(trapID, bs = "re"),
                data = dat_lag_scaled, 
                #family = nb(), 
                correlation = corAR1(form = ~ year_num))







m.gamm1 <- gamm(sum_ips ~ s(veg_tmp, by=temp_fact, k = 4), 
               s(year_num, by = temp_fact, bs = "re") +
                 s(pairID, bs = "re") +
                 s(trapID, bs = "re"),
               data = dat_lag_scaled, 
               #family = nb(), 
               correlation = corAR1(form = ~ year_num))

# Perform the Durbin-Watson test - better for linear models? ACF,pACF are better for gams, GAMM (gams with mixed effects)
dw_result <- dwtest(residuals(m.gam5) ~ fitted(m.gam5))
(dw_result)




# GAM first try, unscaled----------------------------------------------------------------

m1 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            veg_tmp +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)

m2 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            previous_veg_tmp + 
            veg_tmp +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)

m3 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)


m4 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
           # s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)


m5 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)

m6 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)


m7 <- gam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            #s(elev) +
            s(spruce_1986) +
                        veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag, family = tw)

m8 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_1986) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, veg_tmp),
          method="fREML", 
          dat_lag, family = tw)


m9 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_1986) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
             s(x,y, k = 25, bs = 'ds') +
            s(year, veg_tmp),
          method="fREML", 
          dat_lag, family = tw)

m9.nb <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_1986) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x,y, k = 25, bs = 'ds') +
            s(year, veg_tmp),
          method="fREML", 
          dat_lag, family = nb)

m9.poisson <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_1986) +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're') +   # trap pair
               s(x,y, k = 25, bs = 'ds') +
               s(year, veg_tmp),
             method="fREML", 
             dat_lag, family = poisson)

m9.poisson.sc <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
                    s(year, k = 6) +
                    s(elev) +
                    s(spruce_1986) +
                    veg_tmp +
                    spei +
                    s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                    s(pairID, bs = 're') +   # trap pair
                    s(x,y, k = 25, bs = 'ds') +
                    s(year, veg_tmp),
                  method="fREML", 
                  dat_lag, family = poisson, scale = -1)

# add population growth
m10 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
            s(year, k = 6) +
            s(elev) +
            s(spruce_1986) +
            s(population_growth) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            # s(x,y, k = 50) +
            s(year, veg_tmp),
          dat_lag,method="fREML",  family = tw)

m10.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_1986) +
             s(population_growth) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)


# remove s(veg_tmp, year)
m10.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_1986) +
               s(population_growth) +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
               # s(x,y, k = 50) +
              # s(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)

# interaction using tensor# s(veg_tmp, year)
m10.4 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_1986) +
               s(population_growth) +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
              ti(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)


# remove previous beetle sum
m11 <- bam(sum_ips ~ #s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_1986) +
             s(population_growth) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(year, veg_tmp),
           dat_lag,method="fREML",  family = tw)

# elevation as linear
m12 <- bam(sum_ips ~ #s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             elev +
             s(spruce_1986) +
             s(population_growth) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're') +   # trap pair
             # s(x,y, k = 50) +
             s(year, veg_tmp),
           dat_lag,method="fREML",  family = tw)



# add temperature from previous year
m13 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 6) +
               s(elev) +
               s(spruce_1986) +
               s(population_growth) +
               s(previous_veg_tmp) +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)

#previous year pop growth
m13.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 5) +
             s(elev) +
             s(spruce_1986) +
             s(population_growth2, k = 6) +
             s(previous_veg_tmp) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)


# set terms to linear

# include autocorrelation: A, and check for autocorrelation by XXXX test
m13.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
               s(year, k = 5) +
               elev +
               s(spruce_1986) +
               population_growth2 +
               previous_veg_tmp +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)


m13.3.1 <- bam(sum_ips ~ s(previous_sum_ips, k = 4) +
               s(year, k = 5) +
               elev +
               s(spruce_1986) +
               population_growth2 +
               previous_veg_tmp +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)



m13.3.1.temp.ac <- gamm(sum_ips ~ s(previous_sum_ips, k = 4) +
              s(year, k = 5) +
                elev +
                s(spruce_1986) +
                population_growth2 +
                previous_veg_tmp +
                veg_tmp +
                spei +
                s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                s(pairID, bs = 're'),
            data = dat_lag,
            family = tw,
            correlation = corAR1(form = ~ year | trapID))


m1.temp.ac <- gamm(sum_ips ~ s(previous_sum_ips, k = 4) +
                          s(year, k = 5) +
                          #elev +
                          s(spruce_1986) +
                          population_growth2 +
                          previous_veg_tmp +
                          s(veg_tmp) +
                          #spei +
                          s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                          s(pairID, bs = 're'),
                        data = dat_lag,
                        family = tw,
                        correlation = corAR1(form = ~ year | trapID))


m2.temp.ac <- gamm(sum_ips ~ #s(previous_sum_ips, k = 4) +
                     #s(year, k = 5) +
                     #elev +
                     #s(spruce_1986) +
                     #population_growth2 +
                     #previous_veg_tmp +
                     s(veg_tmp) +
                     #spei +
                     s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
                     s(pairID, bs = 're'),
                   data = dat_lag,
                   family = tw,
                   correlation = corAR1(form = ~ year | trapID))



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
             s(spruce_1986) +
             s(population_growth) +
     s(previous_veg_tmp)+
             s(previous_spei) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)


# add temperature from previous year
m13.nb <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_1986) +
             s(population_growth) +
             s(previous_veg_tmp) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = nb(theta = 1.8))


# use as linear terms if the edf = 1
m15 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             year +
             elev +
             spruce_1986 +
             s(population_growth) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


# use as linear terms if the edf = 1
m16 <- bam(sum_ips ~ s(previous_sum_ips, k = 70) +
             year +
             elev +
             spruce_1986 +
             s(population_growth, k =80) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


# use automatic k - does not work well!
m17 <- bam(sum_ips ~ s(previous_sum_ips, k = TRUE) +
             year +
             elev +
             spruce_1986 +
             s(population_growth, k =TRUE) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))




#increase k values
m16.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 130) +
             year +
             elev +
             spruce_1986 +
             s(population_growth, k =130) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw(1.5))


#reduce tweedie
m16.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 130) +
               year +
               elev +
               spruce_1986 +
               s(population_growth, k =130) +
               previous_veg_tmp +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(veg_tmp, year),
             dat_lag,method="fREML",  family = tw)



# increase k values for m13 model

# work on m 13: inceaase k 
m17 <- bam(sum_ips ~ s(previous_sum_ips, k = 100) +
             s(year, k = 6) +
             s(elev) +
             s(spruce_1986) +
             s(population_growth, k = 50) +
             s(previous_veg_tmp, k = 20) +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)

# set linear vars as linear, not smooths
m18 <- bam(sum_ips ~ s(previous_sum_ips, k = 100) +
             s(year, k = 6) +
             elev +
             spruce_1986 +
             s(population_growth, k = 50) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)

# reduce 'k' consistently
m19 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
             s(year, k = 6) +
             elev +
             spruce_1986 +
             s(population_growth, k = 50) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)




# add remained forest as smooth
m18.2 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
             s(year, k = 6) +
             elev +
             #spruce_1986 +
             s(remained_spruce) +
             s(population_growth, k = 20) +
             previous_veg_tmp +
             veg_tmp +
             spei +
             s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
             s(pairID, bs = 're'),# +   # trap pair
           # s(x,y, k = 50) +
           # s(veg_tmp, year),
           dat_lag,method="fREML",  family = tw)


# remove non significant and s termsadd remained forest as smooth
m18.3 <- bam(sum_ips ~ s(previous_sum_ips, k = 50) +
               s(year, k = 6) +
               elev +
               #spruce_1986 +
               remained_spruce +
               s(population_growth, k = 20) +
               previous_veg_tmp +
               veg_tmp +
               spei +
               s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
               s(pairID, bs = 're'),# +   # trap pair
             # s(x,y, k = 50) +
             # s(veg_tmp, year),
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
  spruce_1986 = 2.037,
  trapID = 0.000, 
  pairID = 2.572,
  elev = abs(1.800),                # t-value
  population_growth2 = abs(1.325),  # t-value
  previous_veg_tmp = abs(2.787),        # t-value
  veg_tmp = abs(0.840),                 # t-value
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
ggplot(aes(x = veg_tmp,
           y = agg_doy)) +
  geom_smooth(method = 'glm') +
 # geom_point(alpha = 0.5) +
  ylim(c(105, 190)) +
  theme_bw()


# peak 
p_peak <- dat_lag_scaled %>% 
  ggplot(aes(x = veg_tmp,
             y = peak_doy)) +
  geom_smooth(method = 'glm') +
  #geom_point(alpha = 0.5) +
  ylim(c(105, 190)) +
  theme_bw()

library(ggpubr)
ggarrange(p_agg, p_peak, common.legend = TRUE, legend="bottom", ncol = 2)

hist(dat_lag_scaled_complete$agg_doy)

# try modelsL # if bounded values, skewed: potentially beta regression with tranformation
# Transform your response variable to fit within 0 to 1
dat_lag_scaled_complete$transformed_agg_doy <- (dat_lag_scaled_complete$agg_doy - 92) / (300 - 92)
dat_lag_scaled_complete$tr_peak_doy <- (dat_lag_scaled_complete$peak_doy - 92) / (300 - 92)

# Fit a beta regression model

library(glmmTMB)

# Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
dat_lag_scaled_complete$transformed_agg_doy <- pmin(pmax(dat_lag_scaled_complete$transformed_agg_doy, 1e-4), 1 - 1e-4)
dat_lag_scaled_complete$tr_peak_doy <- pmin(pmax(dat_lag_scaled_complete$tr_peak_doy, 1e-4), 1 - 1e-4)


# try GLM aggregation  -----------------

# add random effects
m.agg4 <- glmmTMB(transformed_agg_doy ~ veg_tmp + previous_sum_ips + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

m.agg5 <- glmmTMB(transformed_agg_doy ~ poly(veg_tmp,2) + previous_sum_ips + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

# add nested trap within pairID - does not improve model
m.agg6 <- glmmTMB(transformed_agg_doy ~ poly(veg_tmp,2) + previous_sum_ips + (1 | pairID/trapID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

# add spei
m.agg7 <- glmmTMB(transformed_agg_doy ~ veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)


# add previous temperature
m.agg8 <- glmmTMB(transformed_agg_doy ~ veg_tmp + previous_veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)


# remove current temp
m.agg9 <- glmmTMB(transformed_agg_doy ~  previous_veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

AIC(m.agg4, m.agg5, m.agg6, m.agg7, m.agg8, m.agg9)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.agg9, 
                                              plot = T)


windows()
plot(allEffects(m.agg8 ))
r2(m.agg8)





### test fr peak population ----------------------------------------------------
# add previous temperature

# add random effects
m.peak4 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_sum_ips + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

m.peak5 <- glmmTMB(tr_peak_doy ~ poly(veg_tmp,2) + previous_sum_ips + (1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

# add nested trap within pairID - does not improve model
m.peak6 <- glmmTMB(tr_peak_doy ~ poly(veg_tmp,2) + previous_sum_ips + (1 | pairID/trapID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

# add spei
m.peak7 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)


# add previous temperature
m.peak8 <- glmmTMB(tr_peak_doy ~ veg_tmp + previous_veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)


# remove current temp
m.peak9 <- glmmTMB(tr_peak_doy ~  previous_veg_tmp + previous_sum_ips + spei +(1 | pairID),
                  family = beta_family(link = "logit"),
                  data = dat_lag_scaled_complete)

AIC(m.peak4, m.peak5, m.peak6, m.peak7, m.peak8, m.peak9)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.peak8, 
                                              plot = T)


windows()
plot(allEffects(m.peak8 ))
r2(m.peak8)


summary(m.peak8)


# what is whatL example for agg and peak:

dat_lag_scaled_complete %>% 
  dplyr::select(c(agg_doy, transformed_agg_doy, peak_doy, tr_peak_doy )) %>% 
  arrange(agg_doy) %>% 
  View()









# standardize the dependent and predictors vals  -------------------------------------------------------
# Standardize beetle counts for each trap over years
dat_lag_standardized <- 
  dat_lag %>%
    ungroup(.) %>% 
  dplyr::select(-c(trapID, pairID, year, x, y)) %>%
  # Apply the scale function to standardize each column
  mutate_all(scale)

# move back removed columns
dat_lag_standardized$trapID <- dat_lag$trapID
dat_lag_standardized$pairID  <- dat_lag$pairID
dat_lag_standardized$year      <- dat_lag$year
dat_lag_standardized$x         <- dat_lag$x
dat_lag_standardized$y         <- dat_lag$y

hist(dat_lag_standardized$sum_ips)

dat_lag_standardized2 <- dat_lag_standardized %>% drop_na() %>% as.data.frame()


m8_st <- bam(sum_ips ~ s(previous_sum_ips) +
            s(year, k = 3) +
            s(elev) +
            s(spruce_1986) +
            veg_tmp +
            spei +
            s(trapID, bs = 're') +  # categorical, account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(year, veg_tmp),
            method="fREML",
            data = dat_lag_standardized2)


m13.3_st <- bam(sum_ips ~ s(previous_sum_ips, k = 70) + 
               s(year, k = 5) + 
               elev + 
               s(spruce_1986) + 
               population_growth2 + 
               previous_veg_tmp + 
               veg_tmp + 
               spei + 
               s(trapID, bs = "re") + 
               s(pairID, bs = "re"),
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


plot_smooth(m4, view = 'veg_tmp')
pXXX<- plot_smooth(m1, 
                          view="veg_tmp", 
                          rm.ranef =TRUE)$fv





windows()
plot(m1, page = 1)


dat_fin %>% 
  ggplot(aes(y = sum_ips,
             x = previous_veg_tmp)) + 
  geom_point()



m1 <- gam(sum_ips ~s(previous_sum_ips, k =50) + 
            s(elev) +
            s(conif_prop  ) +
            s(previous_veg_tmp) +
            s(year, k = 5) +
            s(trapID, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
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
  dplyr::select(-c(trapID, pairID, year))


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
                      #veg_tmp + 
                      veg_tmp_z# +
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
                      #veg_tmp + 
                      veg_tmp_z# +
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
                      #veg_tmp + 
                      veg_tmp_z# +
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

glm1 <- glm(sum_ips ~ conif_prop + elev + sm_z + veg_tmp_z +(1 | trapID),
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
                  #s(veg_prcp, k = 20) +         # precip 
                  s(freq, k = 50), #+         # spruce %
                #s(month, k =7) +  # months
                #s(year, k = 8) +  # year
                #s(x, y, k = 10, bs = 'ds', m = c(1, 0.5)) + # 2D smooth
                #ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) +
                #ti(x, y, TMED, d = c(2,1), bs = c('ds','tp'), 
                #   m = M, k = c(25, 10)) +
                #ti(x, y, veg_prcp, d = c(2,1), bs = c('ds','tp'),
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
  filter(pairID %in% c('Weismain','Pressig', 'Weidenberg' )) #%>% 



dat_fin_sub %>%   print()
  


# Check correlation -------------------------------------------------------

cor(dat_fin[,c('sm', 'vpd', 'veg_tmp', 'spei')], use = "complete.obs")




















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
# conif_prop       elev       sm_z      veg_tmp_z 
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
            s(veg_tmp) +
            te(year, veg_tmp),
          data = dat_fin,
          family = nb(theta = NULL))


m2 <- bam(sum_ips ~ s(year, k = 5) +
            s(veg_tmp),# +
  #          te(year, veg_tmp),
          data = dat_fin,
          family = nb(theta = NULL))

m3 <- bam(sum_ips ~ s(year, k = 5) +
            s(veg_tmp) +
          s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m4 <- bam(sum_ips ~ s(year, k = 5) +
            s(veg_tmp, k = 15) +
            s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m5 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp +
            s(x, y),  
          data = dat_fin,
          family = nb(theta = NULL))

m6 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp +
            s(trapID, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 5))

m7 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp +
            s(trapID, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m8 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp + 
            s(elev) +
            s(trapID, bs = 're') +  # account for how much individual trap contributes to the total variation
            s(pairID, bs = 're') +   # trap pair
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m9 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp + 
            elev +
            s(trapID, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(pairID, bs = 're') +   # random effect for trap pair 
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m10 <- bam(sum_ips ~ s(year, k = 5) +
            veg_tmp + 
            elev +
             spei +# spei does not improve the model
            s(trapID, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(pairID, bs = 're') +   # random effect for trap pair 
            s(x, y),  
          data = dat_fin,
          family = nb(theta = 3))

m11 <- bam(sum_ips ~ s(year, k = 5) +
             conif_prop +
             veg_tmp + 
             elev +
             s(trapID, bs = 're') +  
             s(pairID, bs = 're') +    
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))


m12 <- bam(sum_ips ~ s(year, k = 5) +
            s(conif_prop) +              # does not improve, neither as s(), or linear term
             veg_tmp + 
             elev +
             s(trapID, bs = 're') +  
             s(pairID, bs = 're') +   
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))

m13 <- bam(sum_ips ~ s(year, k = 5) +
             s(conif_prop) +              # does not improve, neither as s(), or linear term
             veg_tmp +
             s(veg_tmp, by = year, k = 7) +   # changes in temperature over year
             elev +
             s(trapID, bs = 're') +  
             s(pairID, bs = 're') +   
             s(x, y),  
           data = dat_fin,
           family = nb(theta = 3))


m14 <- bam(sum_ips ~ s(year, k = 5) +
             s(conif_prop) +              # does not improve, neither as s(), or linear term
             veg_tmp +
             s(year, by = veg_tmp, k = 7) +   # changes year over temperature 
             elev +
             s(trapID, bs = 're') +  
             s(pairID, bs = 're') +   
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
            veg_tmp + 
            elev +
            s(trapID, bs = 're') +  # random effect for individual traps, repeatedly measured
            s(pairID, bs = 're') +   # random effect for trap pair 
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

            s(veg_tmp_z, k = 8) +
            s(veg_tmp, k = 8) +
           # s(sm_z, k = 8) +
            s(sm, k = 8) +
            #s(veg_tmp, k = 20) +
            s(year, k = 5) +

            s(pairID, bs = "re") +
            s(trapID, bs = "re") +
            s(x,y, k = 20, bs = 'ds') +
            ti(x,y, year, d = c(2,1), 
               bs = c('ds','tp'), m = M,
               k = c(20, 5)) +
            ti(x,y, veg_tmp, d = c(2,1), 
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

dat.small <- droplevels(dat_fin[dat_fin$trapID %in% 
                                  levels(dat_fin$trapID)[10:30],]) #len prvych par kombinacii site & zone


full.model <- bam(sum_ips ~ 
                          s(conif_prop, k = 9) +
                          s(elev, k = 9) +
                          
                          s(sm, k = 8) +
                          s(sm_z, k = 8),## +
                          
                  #  s(spei, k = 8) +
                   # s(spei_z, k = 8) +
                    
                     #     s(vpd, k = 8) +
                    #      s(vpd_z, k = 8) +
                         
                   #       s(veg_tmp, k = 8) +
                   #       s(veg_tmp_z, k = 8) +
                    
                   #       s(year, k = 6) +
                  
                   #       s(pairID, bs = "re") +
                   #       s(trapID, bs = "re") +
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
                 pairID = rep(c('a','b', 'c', 'd'),each =6),
                # trap_n = rep(c(1,2),12),
                # trap_ID = paste(pairID, trap_n),
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
            s(pairID, bs = 're') +
            s(pairID, trap_n, bs = 're'),
            data = rats, method = 'REML')







m3 <- bam(sum_ips ~ 
            s(conif_prop, k = 9) +
            s(elev, k = 9) +
            s(sm, k = 8) +
            s(spei, k = 8) +
            s(trapID, bs = 're') +
          
            #s(trapID, pairID, bs = 're') +
            #s(vpd, k = 8) +
            # s(vpd_z, k = 8) +
            # s(spei, k = 8) +
            #s(veg_tmp_z, k = 8) +
            # s(veg_tmp, k = 8) +
            # s(sm_z, k = 8) +
            #s(veg_tmp, k = 20) +
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
    s(veg_tmp_z, k = 5) +  # temperature
    s(year, k = 5)+  # year
    s(trapID, k = 5, bs = 're'),# + #random effect of trap
    # s(pairID, k = 5, bs = 're') + #random effect of trap pairs
    # s(
    #   x,
    #   y,
    #   k = 5,
    #   bs = 'ds',
    #   m = c(1, 0.5)
    # ) + # 2D smooth
    # #ti(veg_tmp_z, year, bs = c('cc', 'tp'), k = c(15, 6)) +
    # ti(
    #   x,
    #   y,
    #   veg_tmp_z,
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
         veg_prcp = as.vector(scale(veg_prcp)),
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
    s(veg_prcp, k = 5)  +         # precip
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
    s(veg_prcp, k = 20)  +         # precip
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
    s(veg_prcp, k = 5)  +         # precip
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
      veg_prcp,
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
#     #s(veg_prcp, k = 5)  +         # precip
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
#       veg_prcp,
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
    s(veg_prcp, k = 10) +
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



