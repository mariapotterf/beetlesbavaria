


load(file="outData/final.Rdata")



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
library(DHARMa)
library(MASS)
library(car)     # for VIF

library(plyr)

# colors
library(RColorBrewer)


# Inspect data ----------------------------------------------------------
# VIF - variance inflation factor
# depepndence between predictors - keep only independent ones
# spatial and temporal autocorrelation

# get korelogram: if line overpass teh dashed line, there is a corelation
acf(dat_fin$sum_ips) # super correlated!


lag(1:5)
head(dat_fin)


scrambled <- slice_sample(
  tibble(year = 2000:2005, value = (0:5) ^ 2),
  prop = 1
)


right <- mutate(scrambled, previous_year_value = lag(value, order_by = year))

arrange(right, year)


# add lag: previous year counts, previous year temperature
dat_fin <- 
  dat_fin %>%
    ungroup(.) %>% 
   group_by(trap_name) %>%
    dplyr::mutate(previous_sum_ips = dplyr::lag(sum_ips, order_by = year),
                  previous_tmp =  dplyr::lag(tmp, order_by = year))
  


dat_fin %>% 
  ggplot(aes(y = sum_ips,
             x = previous_sum_ips)) + 
  geom_point()


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



