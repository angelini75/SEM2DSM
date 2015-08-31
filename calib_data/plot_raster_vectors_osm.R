rm(list=ls())
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
# library(devtools)
# install.packages("latticeExtra")
# install_github("environmentalinformatics-marburg/Rsenal")
library(Rsenal)
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
# load data
d<-read.csv("calib.data-2.1.csv")
# convert to spatial point data frame
coordinates(d) = ~X+Y
proj4string(d) <- posgar98
d<- spTransform(d, wgs84)

# only one layer, all info in popups
mapView(d, burst = F)
