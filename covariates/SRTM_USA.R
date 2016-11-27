############### #### ### ## # - DOWNLOAD SRTM & mosaic it - # ## ### #### ###############
####                     ####  ####

# Purpose        : Download 15 years tiles MODIS product (MOD13Q1); mosaic them; substract AOI; calculate mean and sd per period; 
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : Tomislav Hengl; Ype van der Velde; Wei Shangguan
# Status         : pre-pre-alpha
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# sessionInfo(@RStudio Server) http://h2335862.stratoserver.net:8787/ (12 cores)
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-unknown-linux-gnu (64-bit)


rm(list=ls())
#setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/") # @RStudio desktop
setwd("~/big/USA/SRTM/") # @RStudio server
library(raster)
library(rgdal)

#library(gdalUtils)
library(doParallel)

### from REVERB it is posible to get a csv file with a query result
# http://reverb.echo.nasa.gov/reverb/
files<-read.csv("SRTM_Platte_save_results_csv.csv") 

# to be aware
paste(length(files$Granule.UR)," files here!", sep="")
paste(sum(as.numeric(gsub(",",".",files$Size)))/1000, "Gb!")

# URL of each HDF4 file
files$id <- 1:length(files[,1])
files$Online.Access.URLs <- as.character(files$Online.Access.URLs)
granule <- files[,c(2,10)]
write.table(granule, "granule.txt",
            row.names = F,col.names = F,quote = F, sep=".")

# espesifications of each file
granuleID <- read.csv("granule.txt", sep=".", header=F)
names(granuleID) <- c("tile","ver","format","file","id")

registerDoParallel(cores=8) # each core takes a i value
for(i in seq_along(files$Online.Access.URLs)){
  download.file(files$Online.Access.URLs[i],
                destfile = paste0(granuleID$tile[i],
                                  ".zip"),
                mode="wb",
                method = "wget",
                extra = "--password=Icct5mily --user=angelini75")
}
############### #### ### ## # - IMAGE PROCESSING - # ## ### #### ###############
library(rgdal)
library(gdalUtils)
library(raster)
library(doParallel)
library(maptools)
############ Create mosaic 

# unzip files by hand :(
f <- list.files(path = "/home/mangelini/big/USA/SRTM/",pattern = ".hgt")

# subset extension 
aoi <- readShapePoly("/home/mangelini/big/USA/MODIS/Platte_area.shp")

# define projections
wgs84 <- CRS("+init=epsg:4326")

# assign projection
proj4string(aoi) <- wgs84

x <- list()
for(i in seq_along(f)){
  x[[i]] <- raster(f[i])
}
# if you have a list of Raster objects, you can use do.call
names(x)[1:2] <- c('x', 'y')
x$fun <- mean
x$na.rm <- TRUE
# mosaic
y <- do.call(mosaic, x)

# write raster
writeRaster(y, "mosaic_platte.tif",)
