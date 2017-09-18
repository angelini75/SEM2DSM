rm(list=ls()[])
setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/")

library(raster)
library(sp)
library(rgdal)
library(maptools)
library(doParallel)

## Points over DEM and its derivates
files <- list.files(pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("dem", "twi", "vdchn") 

#define crs
wgs84 <- CRS("+init=epsg:4326")
#UTM14N <- CRS("+init=epsg:32614")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
NAD83.KS.N <- CRS("+init=epsg:2796")


s0 <- readOGR("s0.shp")
plot(s0, cex=0.1)
proj4string(s0) <- NAD83.KS.N
extent(s0)


doParallel::registerDoParallel(cores = 3)

xy <- foreach(i = seq_along(files), .combine = cbind) %dopar%{
  xy <- extract(x = raster(files[i]), y = s0)
  xy
}

xy <- cbind(xy, coordinates(s0))
colnames(xy) <- c(header, "X","Y")
xy <- as.data.frame(xy)
coordinates(xy) <- ~X+Y
proj4string(xy) <- NAD83.KS.N

# list of tif files (modis)
files <- list.files(pattern = ".tif$")[c(-1,-4)]
header <- gsub(".tif", "", files)
header <-  c("evisd", "lstm", "ndwi.a", "ndwi.b") 

# transform projection
xy <- spTransform(x = xy, modis)
xy@data$evisd <- extract(x = raster(files[1]), y = xy)
xy@data$lstm <- extract(x = raster(files[2]), y = xy)
xy@data$ndwi.a <- extract(x = raster(files[3]), y = xy)
xy@data$ndwi.b <- extract(x = raster(files[4]), y = xy)
summary(xy)

xy <- spTransform(x = xy, NAD83.KS.N)

spplot(xy)


#r <- raster(xy, res = 1000)
# proj4string(r)
# plot(rasterize(xy[1],r, background= NA))
# 
# distanceFromPoints(object = , xy = , filename = )

xy <- as.data.frame(xy)
xy$twi <- log10(xy$twi)
xy$vdchn <- log10(xy$vdchn+10)
xy$ndwi.a <- (xy$ndwi.a+10)^.3
xy <- xy[complete.cases(xy),]
#
write.csv(file = "~/Documents/SEM2DSM1/Paper_4/data/xy.csv",x = xy)
