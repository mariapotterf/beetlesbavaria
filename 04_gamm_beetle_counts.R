

#  Try GAMs on beetle data counts


rm(list=ls()) 

# Read my paths -----------------------------------------------------------
source('myPaths.R')



# get libs ----------------------------------------------------------------

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

# Stats
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
library('gganimate')

library(terra)





# Get beetle counts
dat      <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))

ips <- dat %>% 
  filter(art == 'Buchdrucker') %>% 
  filter(doy > 90 &  doy < 270)

str(ips)


# Get climate data for traps:
xy_clim <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = "/"))

# get trap coordinates
xy        <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                  layer = 'xy_3035') # read trap location

class(geom(xy))




# Check if there is any trends in beetles numbers?
# effect: location, effect within year, between years, between locations:

head(ips)
str(ips)

ips$globalid <- factor(ips$globalid)

length(unique(ips$globalid))

ips %>% 
  ggplot(aes(y = fangmenge,
             x = year)) +
  geom_point()
# Dummy example:///////////////////////////////////
# How to interpolate values between dates?
dd <- data.frame(count = c(100, 200, 300,
                           10, 15, 480,
                           300, 50, 40),
                 doy = c(90, 150, 230,
                         87, 135, 235,
                         85, 120, 240),
                 car = rep(c('a', 'b', 'd'), each = 3))

doy_vec = c(70:270)

# ---------------------------------------------------

dd <- data.frame(cumsum = c(4, 12, 14.5,
                            8, 15, 20,
                            10, 20, 40),
                 doy = c(4, 10, 15,
                         5, 11, 14,
                         5, 10, 15),
                 gp = rep(c('a', 'b', 'd'), each = 3))

doy_df = data.frame(doy = rep(c(2:20), 3),
                    gp = rep(c('a', 'b', 'd'), each = 19))

# complete the df by wished doy
dd %>% 
  right_join(doy_df) %>% 
  arrange(gp, doy)

# Solution: 
dd$flag <- 0

res1 <- by(dd, dd$gp, \(x) {
  a <- cbind(merge(x[-3], 
                   data.frame(doy=do.call(seq.int, 
                                          as.list(range(x$doy)))), 
                   all=TRUE), 
             gp=el(x$gp))
  na <- is.na(a$cumsum)
  int <- with(a, spline(doy, 
                        cumsum, 
                        n=nrow(a)))$y
  transform(a, 
            cumsum=replace(cumsum, na, int[na]), 
            flag=replace(flag, is.na(flag), 1))
}) |> do.call(what=rbind)
res1




# from answer:
res1 <- by(dd, dd$gp, \(x) {
  a <-
    cbind(merge(x[-3], data.frame(doy = do.call(
      seq.int, as.list(range(x$doy))
    )), all = TRUE), gp = el(x$gp))
  na <- is.na(a$cumsum)
  int <- with(a, spline(doy, cumsum, n = nrow(a)))$y
  # transform(a, cumsum=replace(cumsum, na, int[na]), flag=replace(flag, is.na(flag), 1))
  df <-
    transform(a,
              cumsum = replace(cumsum, na, int[na]),
              flag = replace(flag, is.na(flag), 1))
  merge(df, data.frame(doy = 2:20, gp = el(a$gp)), all = TRUE) |> {
    \(.) .[order(.$doy),]
  }()
}) |> do.call(what = rbind)
res1







windows()
plot(cumsum ~ doy, res1, type='n')
by(res1, res1$gp, \(x) with(x, points(doy, cumsum, type='b', pch=20, col=flag + 1)))

df <- transform(a, 
                cumsum=replace(cumsum, na, int[na]), 
                flag=replace(flag, is.na(flag), 1))    
merge(df, 
      data.frame(doy=2:20, 
                 gp=el(a$gp)), 
      all=TRUE) |> {\(.) .[order(.$doy), ]}()




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
knots_month <- list(doy = c(4.5, 10.5))
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

















# ----------------------------------------------------------------
# LInk beetle counts data to XY coordinates
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

