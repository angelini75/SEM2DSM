setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
rm(list=ls())
#install.packages('reshape')
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
# install.packages('soiltexture')
#install.packages('corrgram')
# install.packages("plyr")
#library(soiltexture)
# library(reshape)
# library(plyr)
# library(sp)
# library(lattice) # required for trellis.par.set():
# #trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()
# library(corrgram)
# #library(gdalUtils)

library(dplyr)

name <- function(x) { as.data.frame(names(x))} 

d <- read.csv("calib.data-4.2.csv")[,-1]
d8 <- d[d$phw>=8,]
d8 <- d8[!is.na(d8$id.p),]
d8$id.p <- as.factor(d8$id.p)
d8 <- group_by(d8, by = id.p)

#write.csv(summarise(d8, min(X), min(Y), min(top, na.rm = T)), "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/covar2015/ph8.csv")
summarise(.data = d8, n())
a8 <- summarise(d8, min(top, na.rm = T))
names(a) <- c("id.p", "top")

b8 <- unique(merge(a, d8[,c(1,21:85)], by="id.p", all.x = T))

d7 <- d[d$phw<7,]
d7 <- d7[!is.na(d7$id.p),]
d7$id.p <- as.factor(d7$id.p)
d7 <- group_by(d7, by = id.p)
#write.csv(summarise(d7, min(X), min(Y), max(bottom, na.rm = T)), "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/covar2015/ph7.csv")
#write.csv(summarise(d8, min(X), min(Y), min(top, na.rm = T)), "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/covar2015/ph8.csv")
summarise(.data = d7, n())
a7 <- summarise(d7, max(bottom, na.rm = T))
names(a7) <- c("id.p", "bottom")

b7 <- unique(merge(a7, d7[,c(1,21:85)], by="id.p", all.x = T))
corrgram(b7[,c(2,3:16)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
corrgram(b7[,c(2,63:67)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
hist(b7$bottom)

corrgram(b8[,c(2,3:16)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
corrgram(b8[,c(2,63:67)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
hist(b8$top,20)
step(lm(bottom ~ X001_mean + X017_mean + X033_mean + X049_mean + X065_mean +
          X081_mean + X097_mean + X113_mean + X129_mean + X145_mean +
          X161_mean + X177_mean + X193_mean + X209_mean + X225_mean +
          X241_mean + X257_mean + X273_mean + X289_mean + X305_mean +
          X321_mean + X337_mean + X353_mean + X001_std + X017_std +
          X033_std + X049_std + X065_std + X081_std + X097_std +
          X113_std + X129_std + X145_std + X161_std + X177_std +
          X193_std + X209_std + X225_std + X241_std + X257_std +
          X273_std + X289_std + X305_std + X321_std + X337_std +
          X353_std, d7[d$env==1,]))

################
# SEM
library(lavaan)
b7$env <- as.factor(b7$env)

model.h <- '
1:
bottom ~ a*XX2

2:
bottom ~ a*XX2 

3:
bottom ~ dem + vdchn + XX2
'

fit.h <- sem(model.h, data = b7,group = "env",group.label = c("1","2","3"),
             fixed.x = FALSE, verbose=TRUE)#, test = "bollen.stine") for low number of samples
summary(fit.h, fit.measures=T, modindices = F, standardized = F, rsquare = T)

corrgram(b7[b7$env== 3,c(2,3:16,63:67)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
################

















library(Cubist)
summary(Cubist::cubist(x = b7[,c(3:62)],y = b7$bottom[]))





step(lm(top ~ #dem + river + wdist + maxc + mrvbf + slope + 
          #twi + vdchn.y + water + lstm + lstsd + evim + evisd + 
          X001_mean + X017_mean + X033_mean + X049_mean + X065_mean + 
          X081_mean + X097_mean + X113_mean + X129_mean + X145_mean + 
          X161_mean + X177_mean + X193_mean + X209_mean + X225_mean + 
          X241_mean + X257_mean + X273_mean + X289_mean + X305_mean + 
          X321_mean + X337_mean + X353_mean + # X + Y + 
          X001_std + X017_std + X033_std + X049_std + X065_std + X081_std + X097_std +
                    X113_std + X129_std + X145_std + X161_std + X177_std +
                    X193_std + X209_std + X225_std + X241_std + X257_std +
                    X273_std + X289_std + X305_std + X321_std + X337_std +
                    X353_std,
          b[b$top<50,]),direction =  "both")

summary(lm(formula = top ~ dem + river + wdist + maxc + vdchn.y + lstm + 
            evisd + X001 + X033 + X049 + X065 + X081 + X097 + X145 + 
            X177 + X193 + X209 + X225 + X241 + X257 + X289 + X305 + X321 + 
            X337 + X353, data = b))

summary(lm(formula = top ~ dem + river + maxc + slope + vdchn.y + lstm + evim + evisd + 
             XX1 + XX2 + XX3, data = b[b$top>50,]))

summary(lm(formula = top ~ river + wdist + mrvbf + slope + water + lstm + 
             evim + X017 + X049 + X081 + X097 + X113 + X145 + X177 + X209 + 
             X289 + X305, data = b[b$top > 50, ]))


step(lm(top~ dem + river + wdist + maxc + mrvbf + slope + 
          twi + vdchn.y + water + lstm + lstsd + evim + evisd + 
#           X001 + X017 + X033 + X049 + X065 + X081 + X097 + X113 + 
#           X129 + X145 + X161 + X177 + X193 + X209 + X225 + X241 + 
#           X257 + X273 + X289 + X305 + X321 + X337 + X353 + 
          X + Y + 
          XX1 + XX2 + XX3, b),direction =  "both")

summary(lm(formula = top ~ dem + river + wdist + maxc + mrvbf + slope + 
     twi + vdchn.y + water + lstm + evisd + X001 + X017 + X033 + 
     X049 + X065 + X081 + X097 + X113 + X145 + X161 + X177 + X209 + 
     X225 + X241 + X257 + X273 + X289 + X305 + X321, data = b[b$top <= 
                                                                50, ]))


plot(lm(formula = top ~ dem + river + maxc + slope + vdchn.y + lstm + 
             lstsd + evim + evisd + X + Y + XX1 + XX2 + XX3, data = b))

plot(top ~ X273_std, b[b$top<50,])


d$ph <- "ph7"
d$ph[d$phw>=8] <- "ph8"
e <- group_by(d, id.p)
k <- unique(e$id.p[e$ph == "ph8"])
if(e$ph == "ph8") 

e$id.p != k


m <- unique(d[!is.na(d$dren),c(1,20:79)])
m <- m[m$dren == "1.1" | 
         m$dren == "1.2" | 
         m$dren == "2.1" | 
         m$dren == "2.2" | 
         m$dren == "2.3" | 
         m$dren == "3.2" | 
         m$dren == "3.3" | 
         m$dren == "4.3" | 
         m$dren == "4.4",]
m$dren <- as.character(m$dren)
m$dren <- as.factor(m$dren)
model.dren <- C5.0(m[,c(3,5,10,16:61)], m$dren)
model.dren <- C5.0(m[,c(61,31,9,54,47,10,3,5)], m$dren)
summary(model.dren)

pred.m

str(pred.m)

#### prediction ####


################# External Drivers ##############
### rasters of external drivers, conversion to data frame and standardization
# install.packages(c("raster", "maptools", "sp", "rgdal"))
#install.packages("maptools")
library(raster)
library(maptools)
library(sp)
library(rgdal)

setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/")

# sdat files (dem) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "river", "wdist","maxc","mrvbf","slope","twi","vdchn","water") 
pred <- read.csv("mask_231m2.csv")
str(pred)

# tif files (modis)
files_m <- list.files(pattern=".tif$")
header_m <- c("lstm", "lstsd", "evim", "evisd")
# for(i in 1:4){
#   print(spplot(readGDAL(files_m[i])))
# }

# other modis files
files_n <- list.files(path = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/output/",pattern=".tif$")
header_n <- gsub(pattern = ".tif",replacement = "",x = files_n)
header_n <- paste("X",header_n, sep = "")
files_n <- paste("output/", files_n, sep = "")

# pred to spatial data frame
coordinates(pred) <- ~X+Y


#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(pred) <- modis
pred <- spTransform(pred, posgar98)

# However, if you got the data from a RasterLayer
# you can avoid the above and simply do:
# library(raster)
# pts <- rasterToPoints(r, spatial=TRUE)

# extract values from files (.sdat)
stack <- list()
for(i in seq_along(files)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files[i])
  proj4string(stack[[i]]) <- posgar98
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1] # over intersects points with raster 
  stack <- list()
  names(pred@data)[length(pred@data)] <- header[i]
}  

## extract values from modis files 
stack <- list()
# reproject endo to modis projection
pred <- spTransform( pred, modis)
for(i in seq_along(files_m)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header_m[i]
}  

# extract values from the other modis files
pred <- spTransform( pred, modis)
for(i in seq_along(files_n)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_n[i])
  proj4string(stack[[i]]) <- modis # change projection
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header_n[i]
}  
pred <- spTransform( pred, modis)
#image(raster(("mod13q1_tot_mean.tif")))

# clean points out of study area
# pred <- spTransform(pred, posgar98)
pred.df <- as.data.frame(pred)
pred.df$dren <- NA
pred.m <- as.matrix(pred.df)

dren <- data.frame(X = pred.m[,61], Y = pred.m[,62])
dren$dren <- C50::predict.C5.0(object = model.dren, newdata = pred.m, type = "class")

dren.lev <- levels(dren$dren)
dren$dren <- as.numeric(dren$dren)
####rasterize results###
library(sp)

coordinates(dren) <- ~X+Y
proj4string(dren) <- modis
#dren <- spTransform( dren, posgar98)
gridded(dren) <- T
fullgrid(dren) <- T 
spplot(dren, zcol="dren")
slot(slot(dren, "grid"), "cellsize") <- 231.6564
sp::write.asciigrid(x = dren, attr = "dren", fname = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/dren.asc")
#pred.sp <- spTransform(pred.sp, modis)
#spplot(pred.sp)


y <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/mod13q1_tot_mean.tif")
y[!y %in% NA] <- 0
proj4string(y) <- modis
names(y)<-"mask"
plot(y)
res(y)
#pred.sp <- spTransform(pred.sp, modis)
r <- rasterize(x = dren,y = y,background= NA,)
spplot(r)
ADE<- readShapePoly("/media/marcos/L0135974_DATA/UserData/BaseARG/study area/ADE_MODIS.shp")
plot(ADE)
r<-mask(x = r,mask = ADE)
TWI <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/TWI250.sdat")
proj4string(TWI) <- posgar98
#rp[!rp %in% NA] <- 0

rp <- projectRaster(from = r,to = TWI,res = res(TWI), crs = posgar98, method = 'bilinear')
plot(rp[[2:8]])
s<- as.data.frame(summary(rp[[2:8]],digits=4))
names(s) <- names(rp[[2:8]])
s[7,] <- cellStats(stat = sd, rp[[2:8]])
rownames(s)[7] <- "sd"
write.csv(s,"result.statistics.csv")






library(C50)
library(Cubist)
e <- unique(e[,16:80])
e$ph <- e$ph[e$ph=="ph8"]
e$ph <- as.factor(e$ph)

summary.C5.0(C50::C5.0(e[,c(29:74)], e$ph))
C5.
cubist()
