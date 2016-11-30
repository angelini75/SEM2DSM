library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(doParallel)
library(raster)

setwd("~/big/USA/SoilMaps")
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

points <- readShapePoints("points.shp")
platte <- readShapePoly("Soil_platte.shp")

platte@data <- platte@data[,c(3,4,9,10)]
platte@data$MUSYM <- as.character(platte@data$MUSYM)
platte@data$MUKEY <- as.character(platte@data$MUKEY)
platte@data$MUSYM_2 <- as.character(platte@data$MUSYM_2)
platte@data$MUKEY_2 <- as.character(platte@data$MUKEY_2)
platte@data$MUSYM[is.na(platte@data$MUSYM)] <- 
  platte@data$MUSYM_2[is.na(platte@data$MUSYM)]
platte@data$MUKEY[is.na(platte@data$MUKEY)] <- 
  platte@data$MUKEY_2[is.na(platte@data$MUKEY)]
platte@data <- platte@data[,c(1,2)]

points@data <- cbind(points@data, over(x = points, y = platte))

sp <- read.table(file = "Table.txt", header = TRUE, sep = "|")
locations <- as.data.frame(points)
names(locations)[4:5] <- c("X", "Y")
locations$MUSYM <- as.numeric(locations$MUSYM)

m <- as.data.frame(table(locations$MUSYM)[which(table(locations$MUSYM)>500)])
z <- locations

a <- z[which(!(z$MUSYM %in% m$Var1)),]
b <- z[which(z$MUSYM %in% m$Var1),]
for(i in seq_along(m$Var1)) {
  c <- b[b$MUSYM == m$Var1[i],]
  c <- c[sample(nrow(c), 500), ]
  a <- rbind(a, c)
}
hist(a$MUSYM, 10000)
plot(a$X,a$Y, pch=3, cex=0.1)

pedons <- merge(x = a,
                y = sp,
                by.x = "MUSYM",
                by.y = "musym",
                all.x = TRUE)
#

name(pedons)

# select columns
# "X", "Y","MUSYM", "cokey", "comppct_r", "pmkind",
# "hzname", "hzdept_r", "hzdepb_r", "claytotal_l", 
# "claytotal_r", "claytotal_h", "om_l", "om_r", "om_h",
# "cec7_l", "cec7_r", "cec7_h"  
p <- pedons[c(4,5,1,9:10,17,20:22,36:44)]
names(p)
head(p,20)
p <- p[with(p, order(cokey, hzdept_r)), ]

pe <- unique(p[,3:18])


over()
tif <- raster("~/big/USA/SRTM/mosaic_platte.tif")

rsaga.module.exists(module = 1,)
(in.dem = tif, out.slope = slope.tif)
