

# Investigate the Beetle datasets:

# get some overview:
# https://anamtk.github.io/GLM_tutorials/tidyTuesday_fall2020/GLM_inR_Tutorial.html#A_Common_Distributions

# what is the ideal distribution to model my data???

# load the beetle data
# what is my distribution?
# what is my response variable?
# test several: occurences
#     predicted daily datasets
#     peaks?
#     counts until July?


# look from Cornelius data: what distribution did he used? Why?

# how to model my data and what is the best approach?
# What am I looking for?


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')


#library(sf)
library(dplyr)
library(data.table)
library(tidyr)
#library(raster)
library(rgdal)
library(tidyverse)
library(lubridate)
library(patchwork)
library(fasterize)
library(ggpubr)
#library(terra)
library(ggplot2)
library(ggpubr)


# Get beetle counts
dat      <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))



# Get to know my distribution 
# ---------------------------------------------------------

# data: counts 
# counts are discrete data
# test families:
# Poisson - count or aboundance data
# negative binomial - count or abundance with overdispersion
# quasi poisson
# trunkated poisson - Count or abundance data that won't have zeros (e.g. only counted when they were present)
# trunkated negative binomial - Count or abundance data that won't have zeros and that is overdispersed

# Try on example:

# modelled: count data (Poisson)
# predictors: continuous 

# why GLM and not LM? Many ecological data do not fit the assumptions of linear models: e,.g. normal data distribuition
# still assumptions: no weird patterns in residuals, residuals and data do not have a weird skew


# test example: salamander count data
# copmpare the abounmdance if site was mined/not mined

library(glmmTMB)
data(Salamanders)
hist(Salamanders$count)

# histogram is not normal, so GLM is more appropriate than LM

library(MASS)           # glm.nb for negative binomial models
library(glmmTMB)        # Salamanders dataset and lots of families
library(lme4)           # pseudo-R2 and and mixed model functionality
library(MuMIn)          # dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc)         # pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa)         # model diagnostics
library(effects)        # what do my marginal effects look like?
library(performance)    # binomial model diagnostics
library(emmeans)        # post hoc for categorical predictors


# The GLM Process:
# --------------------------------
# 1 .determine model structure: what is X? what is Y?
# 2. model selection process: model comparisons, (AIC)
# 3. do my data fits model assumptions? check distributions of residuals, try new model if npot randomly distributed

# check out the Data:
data(Salamanders) #from glmmTMB

str(Salamanders)


# Try first model: does the aboundance of salamanders depends on the cover/water temperature?
model <- glm(count ~ cover + Wtemp,
             data = Salamanders,
             family = "poisson")

# Explore teh model summary:
summary(model)

# find the best model: by AIC - lower AIC, better model
# AIC compared how model fits, and penalize for more parameters used
# find AIC of teh model:
AIC(model)
MuMIn::AICc(model)  # for this model, AIc and AICc (smaller data) are quite similar

# I can compare only models that come from teh same datasets!! (e.g. not two datasets: one from one stream, one from the other)

# dredge is other function to compare models
# needs to have species na.action in the 'model' string
model <- glm(count ~ cover + Wtemp,
             data = Salamanders,
             family = "poisson",
             na.action = "na.fail")  # 

MuMIn::dredge(model)

# Global model call: glm(formula = count ~ cover + Wtemp, family = "poisson", 
#data = Salamanders, na.action = "na.fail")
#---
#  Model selection table 
#(Intrc)  cover    Wtemp df    logLik   AICc delta weight
#4  0.2407 0.2708 -0.05746  3 -1388.552 2783.1  0.00  0.549
#2  0.2424 0.2780           2 -1389.758 2783.5  0.39  0.451
#3  0.2758        -0.09193  2 -1418.500 2841.0 57.88  0.000
#1  0.2799                  1 -1421.967 2845.9 62.80  0.000
#Models ranked by AICc(x) 

# Dredge evaluated several predictors as individual models:
# and compared them using AICc and Delta
# delta is in log number: magnitude of 2 is a cut-off - I can remove all
# models with delta > 2
# the top performing modesl are first 2: from those, it is better to choose teh simple one: only iuse cover 
# as a predictor


# Do my data fit model assumptions?
best_model <- glm(count ~ cover,
                  data = Salamanders,
                  family = "poisson")

# get model diagnostics: from DHARMa package, can simulate residuals
windows()
simulationOutput <- simulateResiduals(fittedModel = best_model, 
                                      plot = T)
# Yay! the QQ plot shows that we have some issues, residuals are not independent

# test if the residuals are the same as random?
testDispersion(simulationOutput)
# p-val is very low, meaning that residuals are overdispersad!!


# how to fix for overdispersion in the model:

# what to report in paper:

# get estimates: 
summary(best_model)

# get r2:
r2(best_model)

# get confidence interval
MASS::confint(best_model)


# fit my model:
# --------------------------------------------------------------
# is teh salamandra counts dependent on the presence/absence of the 
# mining site?
m0 <- glm(count ~ mined,
          data = Salamanders,
          family = "poisson",
          na.action = "na.fail")
summary(m0)

simulationOutput0 <- simulateResiduals(fittedModel = m0, plot = T)

# test dispersion
testDispersion(simulationOutput0)

# test if I have too many zeors?
testZeroInflation(simulationOutput0)

# simply plot effects: 
plot(allEffects(m0))


# my investigations:
# ------------------------------------------------------
m1 <- glm(count ~ mined + cover + Wtemp,
                  data = Salamanders,
                  family = "poisson",
          na.action = "na.fail")
summary(m1)


# visualize teh data:
Salamanders %>%
  mutate(mine_line = ifelse(mined == "no", 1, 2)) %>%
  ggplot(aes(x = mine_line, y = count)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# check model:
MuMIn::dredge(m1)

# keep only cover and mined
m2 <- glm(count ~ mined + cover,
          data = Salamanders,
          family = "poisson",
          na.action = "na.fail")
MuMIn::dredge(m2)

windows()
simulationOutput <- simulateResiduals(fittedModel = m2, plot = T)


# Need to further fix model for remove teh dependencies in teh data:
# keep only cover and mined
m3 <- glm(count ~ mined + cover,
          data = Salamanders,
          family = "poisson",
          na.action = "na.fail")
MuMIn::dredge(m3)

# try different family:
m.nb <- glm(count ~ mined + cover,
         data = Salamanders,
         family = "nbinom2",
         na.action = "na.fail")
MuMIn::dredge(m.nb)


# try random effects:
m.1 <- glm(count ~ 1,
            data = Salamanders,
            family = "poisson",
            na.action = "na.fail")


# compare poisson model with the nb (negative binomial)
MuMIn::dredge(m.1)

summary(m.1)

AIC(m3, m.nb)



# ------------------------------------------------------------------------
# Follow Kristin's example and evaluate my distributions:
# beetle counts:
# ------------------------------------------------------------------------

library(bbmle)  # for AIC 


# Read my paths -----------------------------------------------------------
source('myPaths.R')

# Get climate data for traps:
xy_clim <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = "/"))


# ----------------------------------------------------------------
# LInk beetle counts data to climate XY
# ----------------------------------------------------------------

head(dat)
head(xy_clim)

# convert data from long to wide to have variables in columns
xy_clim2 <- xy_clim %>%
  dplyr::select(ID, var, value, year, month)# %>% 
  



# seems challenging, as the climate data are by months?
unique(xy_clim$time)  # recorded once a month, from April to October (7 months), years: 2015 - 2021 (7 years) ~ 49 records in 
# need to sum beetle data: sums by months and year! then link them together

ips.sum <- ips %>%
  group_by(year, month, objectid) %>% 
  summarise(ips_month_sum = sum(fangmenge)) %>% 
  filter(year %in% 2015:2021) %>% 
  filter(month %in% 4:10)
  

# Do the same with climate variables:
unique(xy_clim$var)
# [1] "swv"  "u10"  "v10"  "t2m"  "sshf" "tp" 


# Select swm: can be mean or the max??
swm.sum <- xy_clim %>% 
  filter(var == 'swv') %>% 
  group_by(year, month, ID) %>% 
  summarize(swm_sum = sum(value, na.rm = T)) %>% 
  dplyr::rename(objectid = ID ) #%>% 
  #dplyr::select(-c(day, doy))


# temperature
t2m.sum <- xy_clim %>% 
  filter(var == 't2m') %>% 
  group_by(year, month, ID) %>% 
  summarize(t2m_sum = sum(value, na.rm = T)) %>% 
  dplyr::rename(objectid = ID ) #%>% 
#dplyr::select(-c(day, doy))





# Merg data
dat_clim <- ips.sum %>% 
  right_join(swm.sum, by = c("objectid","year", "month" )) %>% 
  right_join(t2m.sum, by = c("objectid","year", "month" )) #%>% 
#  drop_na()


# Check if data correlate??
# relatioship between counts and drought?
dat_clim %>% 
  ggplot(aes(y = ips_month_sum  ,
             x = swm_sum)) +
  geom_point() 




# does drought index predict the beetle counts???
hist(dat_clim$ips_month_sum)
hist(dat_clim$swm_sum)




# try some models fits: 
fit.binom         <- glmmTMB(ips_month_sum ~ 1, data = dat_clim, family=nbinom2)    # negative binomial from glmmTMB
fit.rand.loc      <- glmmTMB(ips_month_sum ~ (1|objectid), data = dat_clim, family=nbinom2)  # negative binomial from glmmTMB
fit.rand.month    <- glmmTMB(ips_month_sum ~ (1|objectid) + month, data = dat_clim, family=nbinom2)  # negative binomial from glmmTMB
fit.rand.month.swm    <- glmmTMB(ips_month_sum ~ (1|objectid) + month + swm_sum , data = dat_clim, family=nbinom2)  # negative binomial from glmmTMB
fit.rand.swm     <- glmmTMB(ips_month_sum ~ (1|objectid) + swm_sum , data = dat_clim, family=nbinom2)  # negative binomial from glmmTMB
fit.rand.t2m.swm    <- glmmTMB(ips_month_sum ~ (1|objectid) + t2m_sum + swm_sum , data = dat_clim, family=nbinom2)  # negative binomial from glmmTMB
# 
# by / I can set a hierarchical structure in teh model! so tehre are two random parameters: 'objection' and 'trap couple'
# fit.rand.t2m.swm.int    <- glmmTMB(ips_month_sum ~ (1|objectid/trap) + t2m_sum*swm_sum , data = dat_clim, family=nbinom2, AR())  # negative binomial from glmmTMB
# by AR - autocorrelation in time
# need to check Marek's approach! 


# Compare the models by the AIC: the lower AIC, better
# also, simpler models are better than more complex
bbmle::AICtab(fit.binom, fit.rand.loc, fit.rand.month, fit.rand.month.swm, fit.rand.swm, fit.rand.t2m.swm, fit.rand.t2m.swm.int)

# scale predictors before the fitting: Z score transformation 

# --------------------------------------


head(dat)

windows()
hist(ips$fangmenge, breaks  =1000 )

ips %>% 
  filter(fangmenge == 0)

# check distribution:
dat %>% 
  filter(art == 'Buchdrucker') %>% 
  ggplot(aes(fangmenge)) + 
  geom_histogram(aes(y =..density..)) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(dat$fangmenge), 
                                         sd = sd(dat$fangmenge)), color="red")# +

# not normally distributed

# strongly left skewed data

# do I have any zeros???


# Specify categorial variables: -----------------------
str(dat)

# The grouping variables need to be declaread as factors
dat$drought_period <- as.factor(dat$drought_period)


# Options for distribution:
# poisson
# negative binomial (alt: overdispersed, )

ips <- dat %>% filter(art == 'Buchdrucker')
  
   
# 
fit.lm      <- lm(fangmenge ~ 1, data = ips) # normal
fit.lmer    <- lmer(fangmenge ~ 1|drought_period, data = ips) # normal + random effect of drought
fit.drought <- glmmTMB(fangmenge ~ (1|drought_period), data = ips, family=nbinom2)  # negative binomial from glmmTMB
fit.year    <- glmmTMB(fangmenge ~ (1|year), data = ips, family=nbinom2)  # negative binomial from glmmTMB

# the beetle data are longitudional, represents repeated measures (over time, same measures on the locations)
# and in pairs: the pairs are more similar to each other then further distant traps
fit2 <- glmmTMB(fangmenge ~           # beetle counts 
                 (1|monsto_name) +   # enter random effect: repeated measures on traps
                 drought_period,               # dependance on year
               data = ips,           # data from  
               family=nbinom2)  # negative binomial from glmmTMB



# Compare the models by the AIC: the lower AIC, better
# also, simpler models are better than more complex
bbmle::AICtab(fit.lm, fit.lmer, fit.drought, fit.year,fit2  )

windows()
simulationOutput <- simulateResiduals(fittedModel = fit2, 
                                      plot = T)

res <- simulateResiduals(fittedModel = fit2)

plot(res, rank = T)
plot(res, asFactor = T)
# Yay! the QQ plot shows that we have some issues, residuals are not independent

# test if the residuals are the same as random?
testDispersion(simulationOutput)


summary(fit)

# simulate data:
# simulate data
sims <- simulate(fit, nsim=100)

windows()
sims %>%
  as_tibble() %>%
  gather() %>%
  ggplot() +
  geom_line(stat = "count", aes(x = value, group = key), alpha = 0.1, col="#d73027") +
  geom_density(stat="count", data = ips, aes(x = fangmenge), col = "black", size = 1) +
  theme_bw() +
  labs(x = "Cover",
       y = "Count")
