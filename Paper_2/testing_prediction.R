rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
#library(lavaan)
library(corrgram)
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
d <-read.csv("calib.data-3.1.csv")[,-1]
name <- function(x) { as.data.frame(names(x))}
############### PRE-PROCESSING ################## 
name(d)
#names(d)[c(5,6,9,10)]<- c("tb.A","sat.A", "oc.A","bt")
#d<-d[!is.na(d$oc.A),]
#d$sat.A[is.na(d$sat.A)] <- 100 # it is assumed 100% saturation when CaCO3 is present 
#d$d.caco3[is.na(d$d.caco3)] <- 300 # it is assumed that CaCO3 is very deep when it is absent within the solum 

#summary(d$tb.A[d$is.caco3==1], omit.na=TRUE)
#d$tb.A[is.na(d$tb.A)] <- 20 # it is assumed 20 cmol+/kg (thrid quartil) in NAs

#summary(d$esp.A[d$d.caco3<15], omit.na=T)
#d$esp.A[is.na(d$esp.A)] <- 16.3 # it is assumed mean of ESP A in NAs

#summary(d$d.caco3[is.na(d$esp.B)])
#summary(d$esp.B[d$d.caco3<30], omit.na=T)
#d$esp.B[is.na(d$esp.B)] <- 30.9 # it is assumed mean of ESP B in NAs

#nas<-d[!complete.cases(d),]
#nas

# # statistics of calibration data
# D<-d[,4:10]
# d.stat<- matrix(data = NA,nrow = 6,ncol = 7,
#                 dimnames = list(c("Min","Median","Mean", "Max", "SD","SS"),names(D)))
# d.stat[1,]<- apply(X = D,FUN = min,2) # 2 means by column
# d.stat[2,]<- apply(X = D,FUN = median,2)
# d.stat[3,]<- apply(X = D,FUN = mean,2)
# d.stat[4,]<- apply(X = D,FUN = max,2)
# d.stat[5,]<- apply(X = D,FUN = sd,2)
# for(i in 1:7){
#   d.stat[6,i] <- sum((mean(D[,i])-D[,i])^2 )
# }
# write.csv(d.stat,"summary.calibdata.csv")


# transformation
#d$esp.A <- log10(d$esp.A)
#d$esp.B <- log10(d$esp.B)
#d$is.hydro <- ordered(d$is.hydro)
#d$is.E <- ordered(d$is.E)
#d$is.caco3 <- ordered(d$is.caco3)

# reduce evi covariates
corrgram(d[,48:70])
d$XX1 <- apply(X = d[,c(48:52,70)],MARGIN = 1,FUN = mean)
d$XX2 <- apply(X = d[,c(54:61)],MARGIN = 1,FUN = mean)
d$XX3 <- apply(X = d[,c(63:67)],MARGIN = 1,FUN = mean)
d <- cbind(d[,c(1:47)], d[,73:75], d[,c(71,72)])

corrgram(d)

# #### data normalization
# N<- data.frame(mean = rep(0,24),sd=rep(0,24))
# 
# # save mean and sd
# dm <- d[,c(4:10,15:31)]
# for(i in 1:24){
#   N$mean[i] <- mean(dm[,i])
#   N$sd[i] <- sd(dm[,i])
#   rownames(N)[i] <-names(dm)[i]
# }
# ## CHECK RESULT
# N[which(rownames(N) == "oc.A"),] == c(mean(d$oc.A),sd(d$oc.A))
# 
# # normalization
# n <- c(4:10,15:31)
# D <- matrix(NA,nrow = 320,ncol = length(n), dimnames = list(NULL,names(d)[n]))
# for(i in 1:ncol(D)){
#   D[,i]<- (d[,n[i]]-mean(d[,n[i]]))/sd(d[,n[i]])
# }


# D
# boxplot(d[,2:34])  
# dimnames(D)
# D<-as.data.frame(D)
covar<- "dem + river + wdist + maxc + mrvbf + slope + twi + vdchn + water + lstm + lstsd + evim + evisd + XX1 + XX2 + XX3"
nam <- names(d)[2:34]

d$arena.mf.C <- log(d$arena.mf.C+0.01)

rep <- matrix(nrow = 33, ncol = 3, dimnames = list(NULL,c("variable","adjr2", "formula")))
for(i in 1:33){
  f <- as.formula(paste(nam[i],"~", covar, sep=""))
  m<-step(lm(f, d),direction = "both")[[11]]
  k <- summary(lm(as.formula(m),d))
  rep[i,1] <- as.character(k[[2]][[2]])
  rep[i,2] <- as.numeric(k$adj.r.squared)
  rep[i,3] <- as.character(m)[2]
}
rep <-  as.data.frame(rep)
rep$adjr2 <- round(rep$adjr2,3)
write.csv(rep, "/home/marcos/Documents/Second_paper_Marcos/Paper_2/rep.csv")
summary(lm(log(arena.mf.C) ~ river + maxc + slope + lstm + X + Y, d))

install.packages("C50")
library(C50)
install.packages("Cubist")
library(Cubist)
h <- (cubist(x = d[,c(35:52)],y = d$clay.B))
summary(h)
str(h)

isE<-(C50::C5.0(x = d[,c(15,17:53)],y = d$is.E))
summary(C50::C5.0(x = d[,15:27],y = d$is.caco3))
summary(C50::C5.0(x = d[,15:27],y = d$is.hydro))
plot(C5.0(x = d[,15:53],y = d$is.E))

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


# tif files (modis)
files_m <- list.files(pattern=".tif$")
header_m <- c("lstm", "lstsd", "evim", "evisd")
# for(i in 1:4){
#   print(spplot(readGDAL(files_m[i])))
# }

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
#image(raster(("mod13q1_tot_mean.tif")))

for(i in seq_along(files_n)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_n[i])
  proj4string(stack[[i]]) <- modis # change projection
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header_n[i]
} 
# clean points out of study area
# pred <- spTransform(pred, posgar98)
pred.df <- as.data.frame(pred)


pred.df <- pred.df[complete.cases(pred.df),] #it may clean some internal pixels
pred.df <- pred.df[,c(2:37)]

for(i in 1:14) {
  pred.df[,length(pred.df)+1] <- NA
  names(pred.df)[length(pred.df)] <- names(d)[i]
}  
pred.df <- pred.df[,c(16:29,1:15)]
pred <- pred.df
pred <- pred[,c(4:10,15:29)]



name(pred)
is.E <- cbind(pred$X,pred$Y,C50::predict.C5.0(object = isE, newdata = pred[,2:]))


library(sp)
library(raster)
is.E.sp <- as.data.frame(is.E)
names(is.E.sp) <- c("X", "Y", "isE")
coordinates(is.E.sp) <- ~X+Y
proj4string(is.E.sp) <- modis

#pred.sp <- spTransform(pred.sp, modis)
#spplot(pred.sp)


y <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/mod13q1_tot_mean.tif")
y[!y %in% NA] <- 0
proj4string(y) <- modis
names(y)<-"mask"
plot(y)

#pred.sp <- spTransform(pred.sp, modis)
r <- rasterize(x = is.E.sp,y = y,background= NA)
plot(r)
ADE<- readShapePoly("/media/marcos/L0135974_DATA/UserData/BaseARG/study area/ADE_MODIS.shp")
plot(ADE)
r<-mask(x = r,mask = ADE)
TWI <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/TWI250.sdat")
proj4string(TWI) <- posgar98
#rp[!rp %in% NA] <- 0

rp <- projectRaster(from = r,to = TWI,res = res(TWI), crs = posgar98, method = 'bilinear')
plot(rp[[2]])
summary(is.E.sp@data[,1])
