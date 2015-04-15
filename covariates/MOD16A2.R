############### #### ### ## # - ALTERNATIVE WAY TO DOWNLOAD MODIS DATA - # ## ### #### ###############
#### #### For MOD16 MOD16A2_MONTHLY.MERRA_GMAO_1kmALB ####
# Purpose : Download 13 years tiles MODIS product (); mosaic them; substract AOI; calculate mean and sd per period;
# Maintainer : Marcos Angelini (marcos.angelini@wur.nl);
# Contributions : 
# Status : pre-pre-alpha
# Note :
# sessionInfo(@RStudio desktop) lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# sessionInfo(@RStudio Server) http://h2335862.stratoserver.net:8787/ (12 cores)
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-unknown-linux-gnu (64-bit)


rm(list=ls())
setwd("~/Documents/SEM2DSM/covariates") # @RStudio desktop
#setwd("~/big") # @RStudio server
#library(raster)
library(rgdal)
#library(gdalUtils)
library(doParallel)
### from REVERB it is posible to get a csv file with a query result
# http://reverb.echo.nasa.gov/reverb/#utf8=%E2%9C%93&spatial_map=satellite&spatial_type=rectangle
# it can be directly read from ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/ (Curl package)
# be aware of repeated files
install.packages("RCurl")
require(RCurl)

M <- c("M01/","M02/","M03/","M04/","M05/","M06/","M07/","M08/","M09/","M10/","M11/","M12/")
Y <- 2000:2014
Y <- paste("Y",Y,"/",sep="")
Y <- as.character(Y)
M <- as.character(M)




a <- NULL
for(i in 1:length(Y)){
  a[i] <- paste("ftp://ftp.ntsg.umt.edu/pub/MODIS/NTSG_Products/MOD16/MOD16A2_MONTHLY.MERRA_GMAO_1kmALB/",Y[i], sep="")
}
URL <- paste(rep(a,each = 12),rep(M, 15),sep="")

grep(pattern = getURL(url = URL[1], ))

MOD16A2.A2000M01.h12v12.105.2013120040418.hdf

files<- NULL
for(i in 1:length(URL)){
  b <- NULL
  b <- strsplit(try(getURL(URL[i],ftp.use.epsv = FALSE, dirlistonly = TRUE)), "\n")
  b <-b[[1]]     
  files <- append(files, b[c(grep("h12v12",b),grep("h12v13",b))])
  Sys.sleep(5)
}
     
     
     
     
     
     
     
     