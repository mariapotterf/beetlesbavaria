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
