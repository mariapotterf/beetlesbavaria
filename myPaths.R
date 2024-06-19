

# ------------------------
# Paths to source
# ------------------------

# Input paths
myPath     = "C:/Users/ge45lep/Documents/2022_BarkBeetles_Bavaria"
inFolder   = "rawData"
outFolder  = "outSpatial"
outTable   = "outTable"
outReach   = "outReach"



# define variables
spring.months         = 3:5
veg.months            = 4:9  # study period
study.period.extended = 2012:2021
reference_period      = 1980:2010  # for anomalies calculation
study_period          = 2015:2021  # for anomalies calculation



# Vars
doy.start  =  91 # April 1st# 60  # March 1st, 
doy.end    = 273  # SEpt 30  304 # Oct 30
veg.period = doy.start:doy.end
