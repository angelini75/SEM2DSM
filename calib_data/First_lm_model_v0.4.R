rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
d <-read.csv("calib.data-2.1.csv")[,-1]
d$sat.A[d$id.p==480] <- 88 #error in dataset
############### PRE-PROCESSING ################## 
names(d)
names(d)[c(5,6,9,10)]<- c("tb.A","sat.A", "oc.A","bt")
d<-d[!is.na(d$oc.A),]
d$sat.A[is.na(d$sat.A)] <- 100 # it is assumed 100% saturation when CaCO3 is present 
d$d.caco3[is.na(d$d.caco3)] <- 300 # it is assumed that CaCO3 is very deep when it is absent within the solum 

summary(d$tb.A[d$is.caco3==1], omit.na=TRUE)
d$tb.A[is.na(d$tb.A)] <- 20 # it is assumed 20 cmol+/kg (thrid quartil) in NAs

summary(d$esp.A[d$d.caco3<15], omit.na=T)
d$esp.A[is.na(d$esp.A)] <- 16.3 # it is assumed mean of ESP A in NAs

summary(d$d.caco3[is.na(d$esp.B)])
summary(d$esp.B[d$d.caco3<30], omit.na=T)
d$esp.B[is.na(d$esp.B)] <- 30.9 # it is assumed mean of ESP B in NAs

nas<-d[!complete.cases(d),]
nas

# statistics of calibration data
D<-d[,4:10]
d.stat<- matrix(data = NA,nrow = 6,ncol = 7,
                dimnames = list(c("Min","Median","Mean", "Max", "SD","SS"),names(D)))
d.stat[1,]<- apply(X = D,FUN = min,2) # 2 means by column
d.stat[2,]<- apply(X = D,FUN = median,2)
d.stat[3,]<- apply(X = D,FUN = mean,2)
d.stat[4,]<- apply(X = D,FUN = max,2)
d.stat[5,]<- apply(X = D,FUN = sd,2)
for(i in 1:7){
  d.stat[6,i] <- sum((mean(D[,i])-D[,i])^2 )
}
write.csv(d.stat,"summary.calibdata.csv")


# transformation
d$esp.A <- log10(d$esp.A)
d$esp.B <- log10(d$esp.B)
d$is.hydro <- ordered(d$is.hydro)
d$is.E <- ordered(d$is.E)
d$is.caco3 <- ordered(d$is.caco3)

#### data normalization
N<- data.frame(mean = rep(0,20),sd=rep(0,20))

# save mean and sd
dm <- d[,c(4:10,15:27)]
for(i in 1:20){
  N$mean[i] <- mean(dm[,i])
  N$sd[i] <- sd(dm[,i])
  rownames(N)[i] <-names(dm)[i]
}
## CHECK RESULT
N[which(rownames(N) == "oc.A"),] == c(mean(d$oc.A),sd(d$oc.A))

# normalization
n <- c(4:10,15:27)
D <- matrix(NA,nrow = 320,ncol = length(n), dimnames = list(NULL,names(d)[n]))
for(i in 1:ncol(D)){
  D[,i]<- (d[,n[i]]-mean(d[,n[i]]))/sd(d[,n[i]])
}

D <- as.data.frame(D)


# step(lm(sat.A ~ dem+wdist+maxc+mrvbf+slope+twi+vdchn+lstm+lstsd+evim+evisd+river, d),direction = "both")
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd, data = d))
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd + tb.A + d.caco3 + river, data = d))
# oc.fit<- lm(formula = oc.A ~ lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi, data = d)
# summary(oc.fit)
# x<-predict(oc.fit,data=d)
# summary(lm(formula = tb.A ~ wdist + lstm + evim + evisd + river, data = d))
# summary(lm(formula = thick.A ~ dem + maxc + evisd, data = d))
# summary(lm(formula = esp.B ~ mrvbf + twi + vdchn + lstsd + evim + evisd + 
#              river, data = d))
# summary(lm(formula = esp.A ~ dem + mrvbf + twi + vdchn + lstsd + evim + 
#              esp.B, data = d))
# summary(lm(formula = sat.A ~ maxc + lstm + lstsd + evim + evisd, data = d))
# summary(lm(formula = d.caco3 ~ dem + mrvbf + vdchn + lstm + lstsd + evim + 
#             evisd + river, data = d))


############### FITTING MODEL ######################

#################### Multiple linear regression #########

tk <- lm(thick.A ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd, D)
oc <- lm(oc.A ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi, D)
tb <- lm(tb.A ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river, D)
sa <- lm(sat.A ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river, D)
ea <- lm(esp.A ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river, D)
eb <- lm(esp.B ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river, D)
bt <- lm(bt ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf, D)

models <- list(tk, oc,tb,sa,ea,eb,bt)

#################### PREDICTION #########################

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

# clean points out of study area
# pred <- spTransform(pred, posgar98)
pred.df <- as.data.frame(pred)


pred.df <- pred.df[complete.cases(pred.df),] #it may clean some internal pixels
pred.df <- pred.df[,c(2:16)]

for(i in 1:14) {
  pred.df[,length(pred.df)+1] <- NA
  names(pred.df)[length(pred.df)] <- names(d)[i]
}  
pred.df <- pred.df[,c(16:29,1:15)]
pred <- pred.df
pred <- pred[,c(4:10,15:29)]
##@## data normalization
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
#P <-read.csv("calib.data-2.1.csv")[,-1]
n <- c(8:20)
for(i in n){
  pred[,i]<- (pred[,i]- N[i,1])/N[i,2]
}

# for(i in n){
#   boxplot(c(d[i],pred[i]))
# }



################# Running Prediction ###########################
#N <- read.csv("N.csv")
library(utils)
library(pastecs)
stat.desc(pred)
pred <- pred[,c(1,6,2:5,7:22)]
names(pred)[1:7] <- paste0(names(pred)[1:7],"r")


for(i in 1:7) {
  pred[,i] <- predict(models[[i]], pred[,8:20])
  }
name(pred)
pred <- pred[,c(1:7,21:22)] #remove external drivers
as.data.frame(rownames(N))
# from unstandardize soil properties
M <- N[c(1,6,2:5,7),] 
name(pred)
for(i in 1:7){
  pred[,i]<- pred[,i]*M$sd[i]+M$mean[i]
}

stat.desc(pred[,1:7])

## 
summary(models[[1]])

# ##### Estimation prediction interval withd ## No stapial #######
# Var.n <-IB%*%V%*%t(IB) # diagonal vaues are variance error
# Var<-diag(Var.n)
# CI <- 1.64 * (Var.n ^ 0.5)
# 
# CI.r <- matrix(0,nrow=7,ncol = 7)
# for(i in 1:7){
#   CI.r[i,i]<- CI[i,i]*M[i,2]
# }
# CI.r<-diag(CI.r)
# names(CI.r)<-rownames(M)
# CI.r[8:9] <- NA
# CI.r[8] <- CI.r[7]
# names(CI.r) <- c(names(CI.r)[1:5], "ll.esp.B", "ul.esp.B", "ll.esp.A", "ul.esp.A")
# CI.r[7] <- 10 ^(CI.r[6] + (Var[6] * M$sd[6]^2) * 0.5)
# CI.r[6] <- 10 ^(CI.r[6] - (Var[6] * M$sd[6]^2) * 0.5)
# CI.r[9] <- 10 ^(CI.r[8] + (Var[7] * M$sd[7]^2) * 0.5)
# CI.r[8] <- 10 ^(CI.r[8] - (Var[7] * M$sd[7]^2) * 0.5)
# 
# write.table(file = "PI.csv",x = CI.r, sep = "\t")


# from log10(ESP) to ESP
pred[,5]<- 10^(pred[,5]+ 0.5 * var(pred[,5])^2)
pred[,6]<- 10^(pred[,6]+ 0.5 * var(pred[,6])^2)

stat.desc(pred[,1:7])

prediction.mean<-apply(pred[,1:7],MARGIN = 2,FUN = mean)
#write.table(file = "prediction.mean.csv",x = prediction.mean, sep = ",")

####rasterize results###
library(sp)
library(raster)
pred.sp <- as.data.frame(pred)
coordinates(pred.sp) <- ~X+Y
proj4string(pred.sp) <- modis

#pred.sp <- spTransform(pred.sp, modis)
#spplot(pred.sp)


y <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/mod13q1_tot_mean.tif")
y[!y %in% NA] <- 0
proj4string(y) <- modis
names(y)<-"mask"
plot(y)
res(y)
#pred.sp <- spTransform(pred.sp, modis)
r <- rasterize(x = pred.sp,y = y,background= NA)
#plot(r[[3]])
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
write.csv(s,"result.statistics.lm.csv")


raster::NAvalue(rp)<--99999
writeRaster(x = rp[[2:8]],filename ="march.lm.tif", overwrite=T,bylayer=TRUE,suffix=names(rp)[2:8])


#################Plotting results####################
#install.packages("semPlot")
library(semPlot)
#matrix to arrange nodes at graph
layout <-as.matrix(read.csv("matrix_semplot.csv")[,-1])
# plot sem using layout=layout

layout[]
semPaths(fit3,what = "std",whatLabels = "no", layout = layout,sizeLat = 4, cut =0.35,
         sizeInt2 = 2,sizeMan =3.5, sizeLat2 = 2, sizeInt = 1,sizeMan2 = 1.5,nCharNodes=0, font=3,
         edge.width = 2,esize=1.5, asize=1,intercepts = F, reorder = F,equalizeManifests =T, residuals = F,layoutSplit=F,
         structural = F, exoCov = F, exoVar=F,cardinal = F,style = "lisrel",#color = c("orange","blue"),
         manifests = c("tb.A", "sat.A", "evim", "evisd", "lstm", "lstsd", "dem","bt", "wdist", "mrvbf", "vdchn", 
                       "twi", "river", "thick.A", "slope", "maxc", "oc.A", "esp.B", "esp.A"))

semPlotModel(fit3)
## not run # extract original matrix
# k<-as.data.frame(semPaths(fit3,what = "est",whatLabels = "no", layout = "tree",sizeLat = 4, sizeInt2 = 2,sizeMan =5, sizeLat2 = 2, sizeInt = 1,sizeMan2 = 1.5,
#          edge.width = 2,esize=1.5, asize=1,intercepts = F, reorder = F,equalizeManifests =T, residuals = F,layoutSplit=F,
#          structural = F, exoCov = F, exoVar=F,cardinal = F,levels = c(1,4,6,8),style = "OpenMX",#color = c("orange","blue"),
#          manifests = c("tb.A", "sat.A", "evim", "evisd", "lstm", "lstsd", "dem","bt", "wdist", "mrvbf", "vdchn", 
#                        "twi", "river", "thick.A", "slope", "maxc", "oc.A", "esp.B", "esp.A"))[[2]][[5]])
# 
# fix(layout)

