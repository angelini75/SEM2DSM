rm(list=ls()[])
setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/")

library(raster)
#library(sp)

## Points over DEM and its derivates
files <- list.files(pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("dem", "twi", "vdchn") 

#define crs
wgs84 <- CRS("+init=epsg:4326")
#UTM14N <- CRS("+init=epsg:32614")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
NAD83.KS.N <- CRS("+init=epsg:2796")

library(maptools)
sa <- readShapeSpatial("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/Platte_area_extended.shp")

spplot(sa)
proj4string(sa) <- wgs84
sa <- spTransform(sa, NAD83.KS.N)
extent(sa)

xcoord <- as.vector(extent(sa))[1:2]
x <- seq(from=xcoord[1], to = xcoord[2], 1000)
ycoord <- as.vector(extent(sa))[3:4]
y <- seq(from = ycoord[1], to = ycoord[2], 1000)
coord <- expand.grid(x,y)

xy <- SpatialPointsDataFrame(coords = coord,data =  as.data.frame(rep(NA, length(coord[,1]))))
proj4string(xy) <- NAD83.KS.N

xy <- xy[rownames(over(sa, xy, returnList = TRUE)$`0`),]


for(i in seq_along(files)){
  xy@data[,i] <- extract(x = raster(files[i]), y = xy)
  names(xy@data)[i] <- header[i]
}

# list of tif files (modis)
files <- list.files(pattern = ".tif$")[c(-1,-4)]
header <- gsub(".tif", "", files)
header <-  c("evisd", "lstm", "ndwi.a", "ndwi.b") 

# transform projection
xy <- spTransform(x = xy, modis)

xy@data[,3+seq_along(files)] <- NA
names(xy@data)[4:7] <- header

# extract values from tiff files
for(i in seq_along(files)){
  xy@data[,3+i] <- extract(x = raster(files[i]), y = xy)
}
xy <- spTransform(x = xy, NAD83.KS.N)

# r <- raster(xy, res = 1000)
# proj4string(r)
# plot(rasterize(xy[1],r, background= NA))
# 
# distanceFromPoints(object = , xy = , filename = )

xy <- as.data.frame(xy)
names(xy)[8:9] <- c("X", "Y")
xy$twi <- log10(xy$twi)
xy$vdchn <- log10(xy$vdchn+10)
xy$ndwi.a <- (xy$ndwi.a+10)^.3

#ST <- read.csv("~/Documents/SEM2DSM1/Paper_4/data/STt.ks.csv")

xy <- xy[,colnames(s)[10:18]]
xy <- as.matrix(xy)
#

#
write.csv(file = "~/Documents/SEM2DSM1/Paper_4/data/xy.csv",x = xy)
