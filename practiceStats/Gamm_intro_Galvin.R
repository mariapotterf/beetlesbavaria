# Learn GAMs
library('here')
library('mgcv')
library('gratia')
library('gamair')
library('ggplot2')
library('purrr')
library('mvnfast')
library("tibble")
library('gganimate')
library('cowplot')
library('tidyr')
library("knitr")
library("viridis")
library('readr')
library('dplyr')
library('gganimate')

library('readr')
library('dplyr')
URL <-  "https://bit.ly/hadcrutv4"
#gtemp <- 
  read_delim('C:/Users/ge45lep/Documents/stats/intro-gam-webinar-2020-master/intro-gam-webinar-2020-master/data/gbtemp.csv',
                    delim = ',', 
                    col_types = 'nnnnnnnnnnnn', 
                    col_names = FALSE) %>%
  select(num_range('X', 1:2)) %>% 
  setNames(nm = c('Year', 'Temperature'))

  
gtemp <- 
 read_csv('C:/Users/ge45lep/Documents/stats/intro-gam-webinar-2020-master/intro-gam-webinar-2020-master/data/gbtemp.csv') %>% 
  select(YEAR, MEASUREMENT) %>% 
    rename(Year = YEAR,
           Temperature = MEASUREMENT)
  
  


library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)





## ----hadcrut-temp-example, echo = FALSE---------------------------------------
library('readr')
library('dplyr')
URL <-  "https://bit.ly/hadcrutv4"
#gtemp <- read_delim(URL, delim = ' ', col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
 # select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))


gtemp <- 
  read_csv('C:/Users/ge45lep/Documents/stats/intro-gam-webinar-2020-master/intro-gam-webinar-2020-master/data/gbtemp.csv') %>% 
  select(YEAR, MEASUREMENT) %>% 
  rename(Year = YEAR,
         Temperature = MEASUREMENT)




## Plot
gtemp_plt <- ggplot(gtemp, aes(x = Year, y = Temperature)) +
  geom_line() + 
  geom_point() +
  labs(x = 'Year', y = expression(Temeprature ~ degree*C))
gtemp_plt




## ----hadcrut-temp-polynomial, echo = FALSE------------------------------------
p <- c(1,3,8,15)
N <- 300
newd <- with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N)))
polyFun <- function(i, data = data) {
  lm(Temperature ~ poly(Year, degree = i), data = data)
}
mods <- lapply(p, polyFun, data = gtemp)
pred <- vapply(mods, predict, numeric(N), newdata = newd)
colnames(pred) <- p
newd <- cbind(newd, pred)
polyDat <- gather(newd, Degree, Fitted, - Year)
polyDat <- mutate(polyDat, Degree = ordered(Degree, levels = p))
gtemp_plt + geom_line(data = polyDat, mapping = aes(x = Year, y = Fitted, colour = Degree),
                      size = 1.5, alpha = 0.9) +
  scale_color_brewer(name = "Degree", palette = "PuOr") +
  theme(legend.position = "right")


## ----read-hadcrut, echo = TRUE------------------------------------------------
library('readr')
library('dplyr')
URL <-  "https://bit.ly/hadcrutv4"
gtemp <- read_delim(URL, delim = ' ', col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
  select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))


## ----show-hadcrut, echo = TRUE------------------------------------------------
gtemp


## ----hadcrutemp-fitted-gam, echo = TRUE, results = 'hide'---------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrutemp-fitted-gam, echo = FALSE--------------------------------------
library('mgcv')
m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
  mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
  geom_point() +
  geom_ribbon(data = pred,
              mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
              fill = "#fdb338") +
  geom_line(data = pred,
            mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
  labs(x = 'Year', y = expression(Temeprature ~ degree*C))


## ----hadcrtemp-plot-gam, echo = FALSE-----------------------------------------
N <- 300
newd <- as_tibble(with(gtemp, data.frame(Year = seq(min(Year), max(Year), length = N))))
pred <- as_tibble(as.data.frame(predict(m, newdata = newd, se.fit = TRUE,
                                        unconditional = TRUE)))
pred <- bind_cols(newd, pred) %>%
  mutate(upr = fit + 2 * se.fit, lwr = fit - 2*se.fit)

ggplot(gtemp, aes(x = Year, y = Temperature)) +
  geom_point() +
  geom_ribbon(data = pred,
              mapping = aes(ymin = lwr, ymax = upr, x = Year), alpha = 0.4, inherit.aes = FALSE,
              fill = "#fdb338") +
  geom_line(data = pred,
            mapping = aes(y = fit, x = Year), inherit.aes = FALSE, size = 1, colour = "#025196") +
  labs(x = 'Year', y = expression(Temeprature ~ degree*C))



## ----birds-1, echo = TRUE-----------------------------------------------------
library('gamair')
data(bird)

bird <- transform(bird,
                  crestlark = factor(crestlark),
                  linnet = factor(linnet),
                  e = x / 1000,  # transform the values to just erase the zeros; help the computation
                  n = y / 1000)
head(bird)

# creastlark data are binary: occurences: 0/1 (also some NA)


## ----birds-2, fig.width = 5, fig.height = 6, echo = FALSE---------------------
ggplot(bird, aes(x = e, 
                 y = n, 
                 colour = crestlark)) + 
  geom_point(size = 0.5) + 
  coord_fixed() + 
  scale_colour_discrete(na.value = '#bbbbbb33') + 
  labs(x = NULL, 
       y = NULL)


## ----birds-gam-1, echo = TRUE-------------------------------------------------
crest <- gam(crestlark ~  # 0/1 data - we can fit binomial family
               s(e,   # two dimensional smoothing: easting and northing 
                 n,   
                    # isotrophic this plate spline=  says that the smooths is same in both directions
                 k = 100), # k is quite hight, as it seems that birds occurence is quite clustered: required more complex model
             data = bird,
             family = binomial,
             method = 'REML')


crest_1000 <- gam(crestlark ~  # 0/1 data - we can fit binomial family
               s(e,   # two dimensional smoothing: easting and northing 
                 n,   
                 # isotrophic this plate spline=  says that the smooths is same in both directions
                 k = 1000), # k is quite hight, as it seems that birds occurence is quite clustered: required more complex model
             data = bird,
             family = binomial,
             method = 'REML')

crest_cr <- gam(crestlark ~  # 0/1 data - we can fit binomial family
                    s(e,   # two dimensional smoothing: easting and northing 
                      n,   
                      # isotrophic this plate spline=  says that the smooths is same in both directions
                      k = 100), # k is quite hight, as it seems that birds occurence is quite clustered: required more complex model
                bs = 'cr',  
                data = bird,
                  family = binomial,
                  method = 'REML')

## ----birds-gam-2, echo = TRUE-------------------------------------------------
summary(crest)


## ----munge-larks, echo = TRUE-------------------------------------------------
## convert back to numeric
## aggregate teh datasets into the QUADRICULA level and then refit the data

bird <- transform(bird,
                  crestlark = as.numeric(as.character(crestlark)),
                  linnet = as.numeric(as.character(linnet)))
## some variables to help aggregation
bird <- transform(bird, 
                  tet.n = rep(1, nrow(bird)),
                  N = rep(1, nrow(bird)), 
                  stringsAsFactors = FALSE)

## set to NA if not surveyed
bird$N[is.na(as.vector(bird$crestlark))] <- NA
## aggregate
bird2 <- aggregate(data.matrix(bird), by = list(bird$QUADRICULA),
                   FUN = sum, na.rm = TRUE)
## scale by Quads aggregated
bird2 <- transform(bird2, e = e / tet.n, n = n / tet.n)

## fit binomial GAM
crest2 <- gam(cbind(crestlark, N - crestlark) ~ s(e, n, k = 100),
              data = bird2, family = binomial, method = 'REML')


## ----crest-3, echo = TRUE-----------------------------------------------------
crest3 <- gam(cbind(crestlark, 
                    N - crestlark) ~
                s(e, n, k = 100),
              data = bird2, 
              family = quasibinomial, # just to see if data is overdispersel? yse, indicated by 'fi'
              method = 'REML')


## ----gam-check-aggregated-lark, echo = TRUE, fig.width = 4.5, fig.height = 4----
ggplot(data.frame(Fitted = fitted(crest2),
                  Resid = resid(crest2)),
       aes(Fitted, Resid)) + geom_point() 


## ----ranefs-------------------------------------------------------------------

# how to deal with random effects: these two models are actually the same!
m_nlme <- lme(travel ~ 1, data = Rail, ~ 1 | Rail, method = "REML") 

m_gam  <- gam(travel ~ s(Rail, bs = "re"), data = Rail, method = "REML")


## ----misspecify, echo = FALSE-------------------------------------------------
set.seed(15)
model_list = c("right model", 
               "wrong distribution",
               "heteroskedasticity",
               "dependent data",
               "wrong functional form")
n <- 60
sigma=1
x <- seq(-1,1, length=n)
model_data <- as.data.frame(expand.grid( x=x,model=model_list))
model_data$y <- 5*model_data$x^2 + 2*model_data$x
for(i in model_list){
  if(i == "right model"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rnorm(n,0, sigma)
  } else if(i == "wrong distribution"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rt(n,df = 3)*sigma
  } else if(i == "heteroskedasticity"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+  
      rnorm(n,0, sigma*10^(model_data[model_data$model==i, "x"]))
  } else if(i == "dependent data"){
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      arima.sim(model = list(ar=c(.7)), n = n,sd=sigma) 
  } else if(i=="wrong functional form") {
    model_data[model_data$model==i, "y"] <- model_data[model_data$model==i, "y"]+ 
      rnorm(n,0, sigma) + ifelse(model_data[model_data$model==i, "x"]>0, 5,-5)
  }
}
ggplot(aes(x,y), data= model_data)+
  geom_point()+
  geom_line(color=ifelse(model_data$model=="dependent data", "black",NA))+
  facet_wrap(~model)+
  geom_smooth(method=gam, formula = y~s(x,k=12),method.args = list(method="REML"))+
  theme(strip.text = element_text(size=16))


## ----sims, include=TRUE,echo=TRUE---------------------------------------------

# Check if model has correct base (EDF), smooths terms (more unexplained variance)
# 

set.seed(2)
n <- 400
x1 <- rnorm(n)
x2 <- rnorm(n)
y_val <- 1 + 2*cos(pi*x1) + 2/(1+exp(-5*(x2)))
y_norm <- y_val + rnorm(n, 0, 0.5)
y_negbinom <- rnbinom(n, mu = exp(y_val),size=10)
y_binom <- rbinom(n,1,prob = exp(y_val)/(1+exp(y_val)))


## ----sims_plot,fig.width = 11, fig.height = 5.5, echo = FALSE-----------------
p1 <- 
  ggplot(data.frame(x = x1, y = y_norm),
             aes(x = x, 
                 y = y)) +
  geom_point()

p2 <- ggplot(data.frame(x = x2, y = y_norm),
             aes(x = x, y = y)) +
  geom_point()

p3 <- ggplot(data.frame(x = x1, y = y_negbinom),
             aes(x = x, y = y)) +
  geom_point()

p4 <- ggplot(data.frame(x = x2, y = y_negbinom),
             aes(x = x, y = y)) +
  geom_point()

p5 <- ggplot(data.frame(x = x1, y = y_binom),
             aes(x = x, y = y)) +
  geom_point()

p6 <- ggplot(data.frame(x = x2, y = y_binom),
             aes(x = x, y = y)) +
  geom_point()

plot_grid(p1, p3, p5, p2, p4, p6, ncol = 3, align = 'hv', axis = 'lrtb')


## ----gam_check_norm1, fig.keep="none", include=TRUE,echo=TRUE, fig.width=11, fig.height = 5.5, fig.align="center"----
# set k low to illustrate when is not correctly used
norm_model_1 <- gam(y_norm~s(x1, k = 4) + 
                      s(x2, k = 4), 
                    method = 'REML')
gam.check(norm_model_1)


## ----gam_check_norm2, fig.keep="none", include=TRUE, echo=TRUE, fig.width=15, fig.height = 5.5,fig.align="center"----
# k is very low (p-value is signifcant for the s-term for x1, need to increase it!)
norm_model_2 <- gam(y_norm ~ 
                      s(x1, k = 12) + 
                      s(x2, k = 4), 
                    method = 'REML')
gam.check(norm_model_2)


## ----gam_check_norm3, fig.keep="none", include=TRUE, echo=TRUE----------------
# now the k for the x1 seems ok, but the k for x2 is too low!
# need to icrease it as well
norm_model_3 <- gam(y_norm ~ s(x1, k = 12) + 
                      s(x2, k = 12),
                    method = 'REML')
gam.check(norm_model_3)


## Plot all 3 models with different k: 4, 12 (x1), 12 (x2)--------------------------------------------
p1 <- draw(norm_model_1)
p2 <- draw(norm_model_2)
p3 <- draw(norm_model_3)

plot_grid(p1, p2, p3, nrow = 3, align = 'hv', axis = 'lrtb')


## ----gam_check_plots1, include=TRUE, echo=TRUE, results="hide"----------------
norm_model <- gam(y_norm ~ s(x1, k=12) + s(x2, k=12), method = 'REML')
gam.check(norm_model, rep = 500)


## ----gam_check_plots2, include=T, echo=TRUE, results="hide"-------------------
pois_model <- gam(y_negbinom ~ 
                    s(x1, k=12) + 
                    s(x2, k=12), 
                  family=poisson, method= 'REML')
gam.check(pois_model, rep = 500)


## ----gam_check_plots3, include=T,echo=TRUE, results="hide"--------------------
negbin_model <- gam(y_negbinom ~ 
                      s(x1, k=12) + 
                      s(x2, k=12), 
                    family = nb, 
                    method = 'REML')
gam.check(negbin_model, rep = 500)


## ----appraise-gam-check-example, fig.height = 5.5-----------------------------
# explore another function that makes the same plots
appraise(negbin_model, method = 'simulate')


## ----setup-shrinkage-example--------------------------------------------------
## an example of automatic model selection via null space penalization
set.seed(3)
n <- 200
dat <- gamSim(1, n=n, scale=.15, dist='poisson', verbose = FALSE) ## simulate data
dat <- transform(dat, x4 = runif(n, 0, 1), x5 = runif(n, 0, 1),
                 f4 = rep(0, n), f5 = rep(0, n))   ## spurious

## ----shrinkage-example-model-fit, echo = TRUE---------------------------------
b <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) +
           s(x4) + s(x5),
         data = dat, family = poisson, method = 'REML',
         select = TRUE)  # this indicated double penalty: allows to select the most appropriate variables
b2 <- gam(y ~ s(x0) + s(x1) + s(x2), 
               data = dat, family = poisson, method = 'REML',
               select = TRUE)  # this indicated double penalty: allows to select the most appropriate variables


b_noPen <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3) +
           s(x4) + s(x5),
         data = dat, family = poisson, method = 'REML',
         select = F)  # this indicated double penalty: allows to select the most appropriate variables


## ----shrinkage-example-truth, echo = FALSE------------------------------------
p1 <- ggplot(dat, aes(x = x0, y = f0)) + geom_line()
p2 <- ggplot(dat, aes(x = x1, y = f1)) + geom_line()
p3 <- ggplot(dat, aes(x = x2, y = f2)) + geom_line()
p4 <- ggplot(dat, aes(x = x3, y = f3)) + geom_line()
p5 <- ggplot(dat, aes(x = x4, y = f4)) + geom_line()
p6 <- ggplot(dat, aes(x = x5, y = f5)) + geom_line()
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, align = 'vh', labels = paste0('x', 1:6))


## ----shrinkage-example-summary------------------------------------------------
summary(b)
summary(b2)
summary(b_noPen)


## ----shrinkage-example-plot---------------------------------------------------
library('gratia'); 
draw(b, scales = 'fixed')
draw(b2, scales = 'fixed')
draw(b_noPen, scales = 'fixed')

plot.gam(b)

AIC(b,b2, b_noPen)


## Get confidence intervals 
# ----setup-confint-example, fig = TRUE, fig.width = 11, fig.height = 5.5, results = "hide", echo = FALSE----
library(mgcv)
set.seed(0)
## fake some data...
f1 <- function(x) {exp(2 * x)}
f2 <- function(x) { 
  0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10 
}
f3 <- function(x) {x*0}

n<-200
sig2 <- 12
x0 <- rep(1:4,50)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
e <- rnorm(n, 0, sqrt(sig2))
y <- 2*x0 + f1(x1) + f2(x2) + f3(x3) + e
x0 <- factor(x0)

## fit and plot...
b <- gam(y ~ x0 + s(x1) + s(x2) + s(x3))

op <- par(mar = c(4,4,1,1) + 0.1)
layout(matrix(1:9, ncol = 3, byrow = TRUE))
curve(f1)
curve(f2)
curve(f3)
plot(b, shade=TRUE)
plot(b, shade = TRUE, seWithMean = TRUE) ## better coverage intervals
layout(1)
par(op)


## ----shrinkage-example-summary------------------------------------------------
summary(b)


## ----aic-example, echo = TRUE-------------------------------------------------
AIC(b1, b2)


## ----aic-chisq, echo = TRUE---------------------------------------------------
# AIC tends to selectd the more complex model (uncorrectly) about 16% of time (0.1572992)
# needs to keep it in mind and slcte the correct model by the df and AIC both! (and brain//)
pchisq(2, 1, lower.tail = FALSE)


## Examples:
# ----co2-example-1, echo = TRUE-----------------------------------------------
data(co2s)
head(co2s)


## ----co2-example-2, echo = TRUE-----------------------------------------------
ggplot(co2s, aes(x = c.month, y = co2)) + geom_line()


## ----co2-example-3, echo = TRUE-----------------------------------------------
b1 <- gam(co2 ~ s(c.month, # cumulative months (0-x)
                 k=300, # expected lot of wigliness
                 bs="cr"), # cubic spline - many data
         data = co2s, 
         method = 'REML')
summary(b1)

gam.check(b1, rep = 500) #  problem in QQdata!!

# only from summary(b) the model looks very good ad explain 100% of deviance:
# so, lets predict the next 36 months!


## ----co2-example-4, echo = TRUE-----------------------------------------------
pd <- with(co2s, 
           data.frame(c.month = 1:(nrow(co2s)+36)))
pd <- cbind(pd, 
            predict(b, pd, se = TRUE))  # predict the next 36 months
pd <- transform(pd, 
                upr = fit + (2*se.fit), # get upper interval 
                lwr = fit - (2 * se.fit)) # get lower interval


## ----co2-example-5, echo = TRUE, fig.height = 5-------------------------------
ggplot(pd, aes(x = c.month, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = co2s, aes(c.month, co2), col = 'red') +
  geom_line(alpha = 0.5)


## ----co2-example-6, echo = TRUE-----------------------------------------------
b2 <- gam(co2 ~ s(month, bs = "cc") + # cyclic smoother (januaries over years will be more similar than other montsh)
            s(c.month, bs = "cr", k = 300), # cubic spline
          data = co2s, 
          method = 'REML',
          knots = list(month = c(0.5, 12.5)))  # 0.5 - mid january, 12.5 - mid december

b3 <- gam(co2 ~ s(month, bs = "cc") + # cyclic smoother
            s(c.month, bs = "cr", k = 300), # cubic spline
          data = co2s, 
          method = 'REML')

gam.check(b2, rep = 500)
gam.check(b3, rep = 500)

## ----co2-example-7, echo = TRUE-----------------------------------------------
summary(b3)
AIC(b1, b2, b3)


## ----co2-example-8, echo = TRUE-----------------------------------------------


nr <- nrow(co2s)
pd2 <- with(co2s, data.frame(c.month = 1:(nr+36),
                             month = rep(1:12, length.out=nr+36)))
pd2 <- cbind(pd2, predict(b2, pd2, se = TRUE))
pd2 <- transform(pd2, upr = fit + (2*se.fit), lwr = fit - (2 * se.fit))


## ----co2-example-9, echo = TRUE, fig.height = 5-------------------------------
ggplot(pd2, aes(x = c.month, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = co2s, aes(c.month, co2), col = 'red') +
  geom_line(alpha = 0.5)


## ----cross-validated, echo = FALSE--------------------------------------------
knitr::include_graphics(here('resources', 'cross-validated.png'))


## galveston model:
# variation of temperatures across the bay, accouning for the 
# days and years 
#----load-galveston-----------------------------------------------------------
galveston <- #read_csv('C:/Users/ge45lep/Documents/stats/intro-gam-webinar-2020-master/intro-gam-webinar-2020-master/data/gbtemp.csv')
  #read_csv(here('data', 'gbtemp.csv')) %>%
  read_csv('C:/Users/ge45lep/Documents/knowledge_base/stats/intro-gam-webinar-2020-master/intro-gam-webinar-2020-master/data/gbtemp.csv') %>% 
  # C:\Users\ge45lep\Documents\knowledge_base\stats\intro-gam-webinar-2020-master\intro-gam-webinar-2020-master\data
  mutate(datetime = as.POSIXct(paste(DATE, TIME), 
                               format = '%m/%d/%y %H:%M', tz = "CDT"),
         STATION_ID = factor(STATION_ID), 
         DoY = as.numeric(format(datetime, format = '%j')),
         ToD = as.numeric(format(datetime, format = '%H')) +
           (as.numeric(format(datetime, format = '%M')) / 60))
galveston


## ----galveston-full-model-----------------------------------------------------
knots <- list(DoY = c(0.5, 366.5))  # knots: beginning and end time

# model: temperature ~ tempDay + tempYear 
m <- bam(MEASUREMENT ~             # temperature
           s(ToD, k = 10) +        # how does temperature vary during a day
           s(DoY, k = 12, bs = 'cc') + # day of the year (variation within the year)
           s(YEAR, k = 30) +       # effect of years: changes over years 
           s(LONGITUDE, LATITUDE,  # account for teh spatial structure 
             k = 100, 
             bs = 'ds', 
             m = c(1, 0.5)) +
           ti(DoY, YEAR, bs = c('cc', 'tp'), k = c(12, 15)) +
           ti(LONGITUDE, LATITUDE, ToD, d = c(2,1), bs = c('ds','tp'),
              m = list(c(1, 0.5), NA), k = c(20, 10)) +
           ti(LONGITUDE, LATITUDE, DoY, d = c(2,1), bs = c('ds','cc'),
              m = list(c(1, 0.5), NA), k = c(25, 12)) +
           ti(LONGITUDE, LATITUDE, YEAR, d = c(2,1), bs = c('ds','tp'),
              m = list(c(1, 0.5), NA), k = c(25, 15)),
         data = galveston, 
         method = 'fREML', # fast REL 
         knots = knots,    # start and end of the cyclic term 'cc'
         nthreads = c(4, 1), discrete = TRUE)

# m_k <- bam(MEASUREMENT ~             # temperature
#            s(ToD, k = 20) +        # how does temperature vary during a day
#            s(DoY, k = 300, bs = 'cc') + # day of the year (variation within the year)
#            s(YEAR, k = 45) +       # effect of years: changes over years 
#            s(LONGITUDE, LATITUDE,  # account for teh spatial structure 
#              k = 100, 
#              bs = 'ds', 
#              m = c(1, 0.5)) +
#            ti(DoY, YEAR, bs = c('cc', 'tp'), k = c(12, 15)) +
#            ti(LONGITUDE, LATITUDE, ToD, d = c(2,1), bs = c('ds','tp'),
#               m = list(c(1, 0.5), NA), k = c(20, 10)) +
#            ti(LONGITUDE, LATITUDE, DoY, d = c(2,1), bs = c('ds','cc'),
#               m = list(c(1, 0.5), NA), k = c(25, 12)) +
#            ti(LONGITUDE, LATITUDE, YEAR, d = c(2,1), bs = c('ds','tp'),
#               m = list(c(1, 0.5), NA), k = c(25, 15)),
#          data = galveston, 
#          method = 'fREML', # fast REL 
#          knots = knots,    # start and end of the cyclic term 'cc'
#          nthreads = c(4, 1), discrete = TRUE)

# make station ID as a smooth as well:
# m_stat <- bam(MEASUREMENT ~  
#              s(STATION_ID) +           # temperature
#              s(ToD, k = 20) +        # how does temperature vary during a day
#              s(DoY, k = 300, bs = 'cc') + # day of the year (variation within the year)
#              s(YEAR, k = 45) +       # effect of years: changes over years 
#              s(LONGITUDE, LATITUDE,  # account for teh spatial structure 
#                k = 100, 
#                bs = 'ds', 
#                m = c(1, 0.5)) +
#              ti(DoY, YEAR, bs = c('cc', 'tp'), k = c(12, 15)) +
#              ti(LONGITUDE, LATITUDE, ToD, d = c(2,1), bs = c('ds','tp'),
#                 m = list(c(1, 0.5), NA), k = c(20, 10)) +
#              ti(LONGITUDE, LATITUDE, DoY, d = c(2,1), bs = c('ds','cc'),
#                 m = list(c(1, 0.5), NA), k = c(25, 12)) +
#              ti(LONGITUDE, LATITUDE, YEAR, d = c(2,1), bs = c('ds','tp'),
#                 m = list(c(1, 0.5), NA), k = c(25, 15)),
#            data = galveston, 
#            method = 'fREML', # fast REL 
#            knots = knots,    # start and end of the cyclic term 'cc'
#            nthreads = c(4, 1), discrete = TRUE)

## ----galveston-simple-model---------------------------------------------------
m.sub <- bam(MEASUREMENT ~
               s(ToD, by=STATION_ID , k = 10) +
               s(DoY, k = 12, bs = 'cc') +
               s(YEAR, k = 30) +
               s(LONGITUDE, LATITUDE, k = 100, bs = 'ds', m = c(1, 0.5)) +
               ti(DoY, YEAR, bs = c('cc', 'tp'), k = c(12, 15)),
             data = galveston, method = 'fREML', knots = knots,
             nthreads = c(4, 1), discrete = TRUE)


## ----galveston-compare-models-aic---------------------------------------------
AIC(m, m.sub)


## ----galveston-compare-models-anova-------------------------------------------
anova(m, m.sub, test = 'F')


## ----galveston-full-model-summary---------------------------------------------
summary(m)


## ----galveston-full-model-plot, fig.height = 5.5------------------------------
plot(m, pages = 1, scheme = 2, shade = TRUE)


## ----galveston-full-model-draw------------------------------------------------
gratia::draw(m, scales = 'free')
hi

## ----galveston-full-predict---------------------------------------------------
pdata <- with(galveston,
              expand.grid(ToD = 12,
                          DoY = 180,
                          YEAR = seq(min(YEAR), max(YEAR), by = 1),
                          LONGITUDE = seq(min(LONGITUDE), max(LONGITUDE), length = 100),
                          LATITUDE  = seq(min(LATITUDE), max(LATITUDE), length = 100)))
fit.g <- predict(m, pdata)
ind.g <- exclude.too.far(pdata$LONGITUDE, pdata$LATITUDE,
                       galveston$LONGITUDE, galveston$LATITUDE, dist = 0.1)
fit.g[ind.g] <- NA
pred.g <- cbind(pdata, Fitted = fit.g)


## ----galveston-full-predict-plot, fig.show = 'hide'---------------------------
plt <- ggplot(pred.g, aes(x = LONGITUDE, y = LATITUDE)) +
  geom_raster(aes(fill = Fitted)) + facet_wrap(~ YEAR, ncol = 12) +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma', na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'right')
plt


## ----galveston-full-predict-plot, echo = FALSE--------------------------------
plt <- ggplot(pred.g, aes(x = LONGITUDE, y = LATITUDE)) +
  geom_raster(aes(fill = Fitted)) + facet_wrap(~ YEAR, ncol = 12) +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma', na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'right')
plt


## ----galveston-animation, echo = FALSE, results = 'hide'----------------------
p <- ggplot(pred, aes(x = LONGITUDE, y = LATITUDE, frame = YEAR)) +
  geom_raster(aes(fill = Fitted)) +
  scale_fill_viridis(name = expression(degree*C), option = 'plasma', na.value = 'transparent') +
  coord_quickmap() +
  theme(legend.position = 'right')+
  labs(x = 'Longitude', y = 'Latitude')

anim <- p + transition_time(YEAR) +
  ggtitle('Year {round(frame_time, 0)}')

anim <- animate(anim,
                nframes = 200, height = anim_height, width = anim_width, res = 100,
                dev = anim_dev)

anim_save('./resources/galveston-animation.gif', anim)


## ----galveston-trends-by-month, fig.show = 'hide'-----------------------------
pdata <- with(galveston,
              expand.grid(ToD = 12,
                          DoY = c(1, 90, 180, 270),
                          YEAR = seq(min(YEAR), max(YEAR), length = 500),
                          LONGITUDE = -94.8751,
                          LATITUDE  = 29.50866))

fit <- data.frame(predict(m, newdata = pdata, se.fit = TRUE))
fit <- transform(fit, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))
pred <- cbind(pdata, fit)

plt2 <- ggplot(pred, aes(x = YEAR, y = fit, group = factor(DoY))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey', alpha = 0.5) +
  geom_line() + facet_wrap(~ DoY, scales = 'free_y') +
  labs(x = NULL, y = expression(Temperature ~ (degree * C)))
plt2


## ----galveston-trends-by-month, echo = FALSE----------------------------------
pdata <- with(galveston,
              expand.grid(ToD = 12,
                          DoY = c(1, 90, 180, 270),
                          YEAR = seq(min(YEAR), max(YEAR), length = 500),
                          LONGITUDE = -94.8751,
                          LATITUDE  = 29.50866))

fit <- data.frame(predict(m, newdata = pdata, se.fit = TRUE))
fit <- transform(fit, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))
pred <- cbind(pdata, fit)

plt2 <- ggplot(pred, aes(x = YEAR, y = fit, group = factor(DoY))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey', alpha = 0.5) +
  geom_line() + facet_wrap(~ DoY, scales = 'free_y') +
  labs(x = NULL, y = expression(Temperature ~ (degree * C)))
plt2

