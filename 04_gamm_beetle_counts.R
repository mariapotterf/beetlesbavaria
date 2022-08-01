

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
library('itsadug')





# Get beetle counts
dat      <- fread(paste(myPath, outFolder, "dat.csv", sep = "/"))

# Get SPEI and clim data:
df_spei   <- fread(paste(myPath, outTable, 'xy_spei.csv'))
df_prec   <- fread(paste(myPath, outTable, 'xy_precip.csv', sep = '/'))
df_temp   <- fread(paste(myPath, outTable, 'xy_temp.csv', sep = '/'))

# Get coniferous cover data: coniferous == 2! 
# spruce: share of spruce
df_tree   <- fread(paste(myPath, outTable, 'xy_treeComp.csv', sep = '/'))
df_spruce <- fread(paste(myPath, outTable, 'xy_spruce.csv', sep = '/'))


# convert spei to wide format: 
df_spei2 <- df_spei %>% 
  pivot_wider(names_from =  scale, values_from = Series.1, names_prefix = "spei")

# convert to df_spei to month, year format: 

df_spei3 <-
  df_spei2 %>%
  mutate(date = format(as.Date(date, "%d/%m/%Y"), "%m.%Y")) %>%
  separate(date, c('month', 'year')) %>% 
  mutate(month = as.numeric(month),
         year = as.numeric(year))


# rename column names to join the datasets:
df_prec <- df_prec %>% 
  rename(PRCP = vals)

df_temp <- df_temp %>% 
  rename(TMED = vals) %>% 
  mutate(TMED = TMED/10)  # as in the description

# join data
df_clim <- df_prec %>% 
  full_join(df_temp, by = c("globalid", "month", "year")) %>% 
  full_join(df_spei3, by = c("globalid", "month", "year")) #%>% 
  

# select only coniferous: == 2, 0 is no forest (eg. covers the whole 500 m buffer)
df_conif <- df_tree %>% 
  filter(species == 2)


# plot TEMP and SPEI --------------------------------------------------------
df_spei <- df_spei %>%
  mutate(date = format(as.Date(date, "%d/%m/%Y"), "%m.%Y")) %>%
  separate(date, c('month', 'year')) %>%
  mutate(month = as.numeric(month),
         year = as.numeric(year)) #%>%

df_spei %>% 
  ggplot(aes(x = factor(year),
             y = Series.1,
             fill = factor(scale))) +
  geom_boxplot(outlier.shape = NA ) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1, lty = 'dashed') +
  theme_bw()+
  facet_grid(.~scale) 


# Classify the SPEI data:
# SPEI values can be categorized into extremely dry (SPEI ≤ −2), 
# severely dry (−2 < SPEI ≤ −1.5), 
# moderately dry (−1.5 < SPEI ≤ −1), and 
# near normal conditions (−1 < SPEI < + 1) (Slette et al., 2019).
range(df_spei$Series.1, na.rm = T)
#  -3.076544  2.980012

df_spei <- df_spei %>% 
  mutate(spei_cat = case_when(Series.1 <= -2 ~ 'ext_dry',
                              Series.1 > -2 & Series.1 <= -1.5 ~ 'sev_dry',
                              Series.1 > -1.5 & Series.1 <= -1 ~ 'mod_dry',
                              Series.1 > -1 & Series.1 <= +1 ~ 'normal',
                              Series.1 > +1 & Series.1 <= +1.5 ~ 'mod_wet',
                              Series.1 > +1.5 & Series.1 <= +2 ~ 'sev_wet',
                              Series.1 > 2 ~ 'ext_wet'
                              ))


df_spei <- df_spei %>% 
  mutate(spei_cat = factor(spei_cat,
                           levels = c('ext_dry','sev_dry', 'mod_dry','normal','mod_wet','sev_wet','ext_wet' )))

# investigate NA: they orccuf ro each location in year 2000
#df_spei %>% 
#  filter(is.na(Series.1)) %>% 
#  distinct(globalid)


# Get frequency of extremely dry locations:
# 
library(RColorBrewer)
display.brewer.pal(7, "BrBG")

df_spei %>% 
  filter(scale == 1) %>% 
  group_by(year, spei_cat) %>% 
  tally() %>% 
  #filter(spei_cat == 'ext_dry'|spei_cat == 'sev_dry') %>% 
  filter(!is.na(spei_cat     )) %>% 
  ggplot(aes(x = year,
             y = n,
             fill = spei_cat)) +
  geom_col(position = "stack") + 
  scale_fill_brewer(type = 'div', palette = "RdBu", "SPEI") + #"RdBu") #BrBG
 theme_bw()


df_temp %>% 
  ggplot(aes(x = factor(year),
             y = TMED)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_hline(yintercept = mean(df_temp$TMED), col="red") +
  #geom_hline(yintercept = -1, lty = 'dashed') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))



# Filter beetle data: get only ips ---------------------------------------------


ips <- dat %>% 
  filter(year !=2014) %>% 
  filter(art == 'Buchdrucker') %>% 
  filter(doy > 92 &  doy < 305) # April 1st to Oct 30

str(ips)


# Get climate data for traps: --------------------------------------------------
# this climate data are for soild drought!! not SPEI
# xy_clim <- fread(paste(myPath, outTable, 'xy_clim.csv', sep = "/"))
#
# get trap coordinates
xy        <- vect(paste(myPath, outFolder, "xy_3035.gpkg", sep = "/"), 
                  layer = 'xy_3035') # read trap location

class(geom(xy))

xy_4326 <- terra::project(xy, "epsg:4326")
df_xy <- data.frame(x = geom(xy_4326)[,'x'],
                    y = geom(xy_4326)[,'y'],
                    globalid = xy_4326$globalid)

# Check if there is any trends in beetles numbers? -----------------------------------------
# effect: location, effect within year, between years, between locations:

head(ips)
str(ips)

ips$globalid <- factor(ips$globalid)

length(unique(ips$globalid))



# Dummy example: Extrapolate values over days SKIP!!!:----------------------------------------------
# How to interpolate values between dates?
dd <- data.frame(count = c(100, 200, 300,
                           10, 15, 480,
                           300, 50, 40),
                 doy = c(90, 150, 230,
                         87, 135, 235,
                         85, 120, 240),
                 car = rep(c('a', 'b', 'd'), each = 3))

doy_vec = c(70:270)

#  --------------------------------------------------

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




# from answer: ------------------------------------------
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






# model raw counts??? !!! -----------------------------------
ips2 <- ips %>% 
  left_join(df_clim, by = c("globalid", "year", "month")) %>% 
  left_join(df_conif , by = c("globalid")) %>% 
  left_join(df_xy, by = c("globalid"))


# try the predicton using the raw data

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


ggplot(data = ips2, aes(y = fangmenge, 
                        x = spei1)) +
  geom_point() +
  geom_line(aes(y = fitted(m_counts)),
            colour = "red", size = 1.2) +
  theme_bw()




windows()
appraise(m_counts, method = 'simulate')
plot(m5, page = 1, shade = T)

gam.check(m5)
k.check(m5)
summary(m5)



# Test XY relationship patter: observed vs modelled data ------------------
# do not run
# linear_model <- gam(Sources ~ SampleDepth, data = isit2)
# summary(linear_model)
# 
# data_plot <- ggplot(data = isit2, aes(y = Sources, x = SampleDepth)) + 
#   geom_point() +
#   geom_line(aes(y = fitted(linear_model)),
#             colour = "red", size = 1.2) + 
#   theme_bw()
# data_plot
# 
# gam_model <- gam(Sources ~ s(SampleDepth), data = isit2)
# 

# get together montly ips counts and spei: ------------------------------------------
ips_sum <- ips %>%
  filter(representativ == 'Ja') %>% 
  group_by(globalid, year, month, monsto_name ) %>% 
  summarise(ips_sum = sum(fangmenge, na.rm = T))


#  Add SPEI, temp, precip, spruce share data --------------------------------------------
ips_sum2 <- ips_sum %>% 
  left_join(df_clim, by = c("globalid", "year", "month")) %>% 
  left_join(df_conif , by = c("globalid")) %>% 
  left_join(df_xy, by = c("globalid"))


# Investigate the Y distribution -----------------------------------
# check for zeros, for NA
min(ips_sum2$ips_sum)

ips_sum2 %>% 
  filter(ips_sum == 0) %>% 
  tally() %>% 
  distinct(globalid) %>% 
  tally()

# about 400 records with 0 values, across 180 locations
# has important impact!!
# check for NA p none
ips_sum2 %>% 
  filter(is.na(ips_sum)) #%>% 
  tally() %>% 
  distinct(globalid) %>% 
  tally()
  
# What distribution?
  # https://towardsdatascience.com/insurance-risk-pricing-tweedie-approach-1d71207268fc
  # poisson - only for counts
  # gamma - does not take zero values
  # tweete - can handle zeros values
  
  # Select optimal 'p' for tweetie - varies between 1 to2
  
  
  
# Check parameters for Tweedie
  

hist(ips_sum2$ips_sum, breaks = seq(0, 52000, 500) )
hist(ips_sum2$freq, breaks = seq(0, 1, 0.1) )
hist(ips_sum2$spei12)


# convert to factors:
ips_sum2 <- ips_sum2 %>% 
  mutate(globalid = factor(globalid),
         monsto_name = factor(monsto_name))

plot(ips_sum2$spei1, ips_sum2$month )


# Get GAM with all predictors: temp, spei, conif, XY ------------------------------------
# check the dstribution:
hist(ips_sum2$ips_sum)

plot(ips_sum2$spei6, ips_sum2$ips_sum)

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

### Test gam only with standardized y value: -------------
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


# GAM: monthly sums -------------------------------------------
knots_month <- list(doy = c(3.8, 10.5))

m <- bam(ips_sum ~ s(spei1) +
           #s(spei3) +
           s(spei6) +
           s(TMED, k = 100) +
           s(PRCP) +
           s(freq) +
           s(month, k = 7, bs = 'cc') +
           s(year, k=8) +
           ti(TMED, spei6), 
         ips_sum2, 
         family = tw(link = "log"))



k.check(m)
plot(m, page = 1, shade=T)
appraise(m)
summary(m)

m1 <- bam(ips_sum ~ s(spei1, k = 500) +
            s(TMED, k =80) + 
            s(month, k = 7, bs = 'cc') +
            s(year, k = 8), # +
            #s(globalid, k = 200, bs = 're'),  # 're' = random effect
          data = ips_sum2,
          family = tw(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

m2 <- bam(ips_sum ~ s(spei1, k = 500, bs = 'cr') +
            s(TMED, k =100) + 
            s(PRCP, k =80) + 
            s(month, k = 7, bs = 'cc') +
            s(year, k = 8) +
          s(globalid,  bs = 're'),  # 're' = random effect
          data = ips_sum2,
          family = nb, #negative binomial(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

m3 <- bam(ips_sum ~ s(spei1, k = 500, bs = 'cr') +
            s(TMED, k =100) + 
            s(PRCP, k =80) + 
            s(month, k = 7, bs = 'cc') +
            s(year, k = 8) +
            s(freq) +
            s(globalid,  bs = 're'),  # 're' = random effect
          data = ips_sum2,
          family = nb, #negative binomial(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REL 
          knots = knots_month ,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

m4 <- bam(ips_sum ~ s(spei1, k =800) +
            s(TMED, k=150) + 
            s(PRCP) + 
            s(month, k =6) +
            s(year, k = 8) +
            s(freq) +
            s(x, y, 
              k = 100, 
              bs = 'ds', 
              m = c(1, 0.5)),# +
            #s(globalid,  bs = 're'),  # 're' = random effect
          data = ips_sum2,
          family = tw(), #negative binomial(),    #Tweedie(1.25,power(.1)),#nb, #twlss,
          method = 'fREML', # fast REML
          select=F,
          knots = knots_month,    # start and end of the cyclic term 'cc'
          nthreads = c(4, 1), 
          discrete = TRUE)

M <- list(c(1, 0.5), NA)
m5 <- bam(ips_sum ~
            s(spei1, k =80) + # drought
            s(TMED, k=100) +  # temperature
            s(PRCP, k = 100) +         # precip 
            s(freq, k = 100) +         # spruce %
            s(month, k =6) +  # months
            s(year, k = 8) +  # year
            s(monsto_name, k = 10, bs = 're') +
            s(x, y, k = 200, bs = 'ds', m = c(1, 0.5))  + # 2D smooth
           ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6))  +
           ti(x, y, TMED, d = c(2,1), 
              bs = c('ds','tp'), 
              m = M, k = c(50, 10)) +
           ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
             m = M, k = c(50, 15))  +
            ti(x, y, spei1, d = c(2,1), bs = c('ds','tp'),
               m = M, k = c(50, 15)),
          data = ips_sum2, 
          method = 'fREML',
          #family = Tweedie(p=1.1, link = power(0)),
          family = tw(), #negative binomial(),
         # knots = knots,
          nthreads = 4, 
          discrete = TRUE)



M <- list(c(1, 0.5), NA)
m7 <- bam(ips_sum ~
            s(spei1, k =5) + # drought
            s(TMED, k=5) +  # temperature
            s(PRCP, k = 5) +         # precip 
            s(freq, k = 5) +         # conif %
            s(month, k =2, 'cc') +  # months
            s(year, k = 2) +  # year
            s(globalid, k = 5, bs = 're') +
            s(monsto_name, k = 5, bs = 're') +
            s(x, y, k = 5, bs = 'ds', m = c(1, 0.5))  + # 2D smooth
            ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)) ,# +
            # ti(x, y, TMED, d = c(2,1), 
            #    bs = c('ds','tp'), 
            #    m = M, k = c(20, 10)) +
            # ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(20, 15))  +
            # ti(x, y, spei1, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(20, 15)),
          data = ips_sum2, 
          method = 'fREML',
          #family = Tweedie(p=1.1, link = power(0)),
          family = tw(), #negative binomial(),
          knots = knots_month ,
          nthreads = 4, 
          discrete = TRUE)



summary(m7)
plot(m7, page = 1)
appraise(m7)


m5.glm <- gam(ips_sum ~
            spei1 + # drought
            TMED,# +  # temperature
            # s(PRCP, k = 1) +         # precip 
            # s(freq, k = 1) +         # spruce %
            # s(month, k =1) +  # months
            # s(year, k = 1) +  # year
            # s(x, y, k = 1, bs = 'ds', m = c(1, 0.5))  + # 2D smooth
            # ti(TMED, year, bs = c('cc', 'tp'), k = c(1, 1)),#  +
            # # ti(x, y, TMED, d = c(2,1), 
            #    bs = c('ds','tp'), 
            #    m = M, k = c(50, 10)) +
            # ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(50, 15))  +
            # ti(x, y, spei1, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(50, 15)),
          data = ips_sum2, 
          method = 'fREML',
          #family = Tweedie(p=1.1, link = power(0)),
          family = 'poisson') #,#tw())#, #negative binomial(),
          # knots = knots,
         # nthreads = 4)

appraise(m5.glm)


m6 <- bam(ips_sum ~
            s(spei1, k =80) + # drought
            s(TMED, k=100) +  # temperature
            s(PRCP, k = 100) +         # precip 
            s(freq, k = 100) +         # spruce %
            s(month, k =6) +  # months
            s(year, k = 8) +  # year
            s(x, y, k = 200, bs = 'ds', m = c(1, 0.5))  + # 2D smooth
            ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6))  +
            ti(x, y, TMED, d = c(2,1), 
               bs = c('ds','tp'), 
               m = M, k = c(50, 10)) +
            ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
               m = M, k = c(50, 15))  +
            ti(x, y, spei1, d = c(2,1), bs = c('ds','tp'),
               m = M, k = c(50, 15)),
          data = ips_sum2, 
          method = 'fREML',
          #family = Tweedie(p=1.1, link = power(0)),
          family = tw(), #negative binomial(),
          # knots = knots,
          nthreads = 4, 
          discrete = TRUE)



windows()
appraise(m5, method = 'simulate')
plot(m5, page = 1)

m6 <- bam(ips_sum ~
            s(spei1, k =800) + # drought
            s(TMED, k=150) +  # temperature
            s(PRCP, k = 100) +         # precip 
            s(freq, k = 400) +         # spruce %
            s(month, k =6) +  # months
            s(year, k = 8) +  # year
            s(x, y, k = 200, bs = 'ds', m = c(1, 0.5)) + # 2D smooth
            ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)),# +
           # ti(x, y, TMED, d = c(2,1), bs = c('ds','tp'), 
           #    m = M, k = c(50, 10)) +
            # ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(50, 15)) +
            # ti(x, y, spei, d = c(2,1), bs = c('ds','tp'),
            #    m = M, k = c(50, 15)),
          data = ips_sum2, 
          method = 'fREML',
          family = tw(), #negative binomial(),
          knots = knots,
          nthreads = 8, 
          discrete = TRUE)

# Usa scat family to account for heavy tail distribution
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/scat.html
m7 <- bam(ips_sum ~
            s(spei1, k =800) + # drought
            s(TMED, k=150) +  # temperature
            s(PRCP, k = 100) +         # precip 
            s(freq, k = 400) +         # spruce %
            s(month, k =6) +  # months
            s(year, k = 8) +  # year
            s(x, y, k = 200, bs = 'ds', m = c(1, 0.5)) + # 2D smooth
            ti(TMED, year, bs = c('cc', 'tp'), k = c(15, 6)),# +
          # ti(x, y, TMED, d = c(2,1), bs = c('ds','tp'), 
          #    m = M, k = c(50, 10)) +
          # ti(x, y, PRCP, d = c(2,1), bs = c('ds','tp'),
          #    m = M, k = c(50, 15)) +
          # ti(x, y, spei, d = c(2,1), bs = c('ds','tp'),
          #    m = M, k = c(50, 15)),
          data = ips_sum2, 
          method = 'fREML',
          family=scat(link="identity"),
          #family = tw(), #negative binomial(),
          knots = knots,
          nthreads = 8, 
          discrete = TRUE)

# m6 is a subset of m5, excluding lowimportance interactions. CHeck with AIC if it is really better?
AIC(m5,m6,m7)

appraise(m7, method = 'simulate')



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

