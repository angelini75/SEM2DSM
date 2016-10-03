rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
d <-read.csv("calib.data-2.1.csv")[,-1]

############### PRE-PROCESSING ################## 
names(d)
names(d)[c(5,6,9,10)]<- c("tb.A","sat.A", "oc.A","bt")
d$sat.A[d$id.p==480] <- 88 #error in dataset
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
d.stat[6,] <- d.stat[5,]/d.stat[3,]
dimnames(d.stat)[[1]][6] <- "CV"
write.csv(d.stat,"summary.calibdata.csv")


# transformation
d$esp.A <- log10(d$esp.A)
d$esp.B <- log10(d$esp.B)

# save mean and sd
N <- data.frame(mean = rep(0,20), sd = rep(0,20), SStot=rep(0,20))
dm <- d[,c(4:10,15:27)]
for(i in 1:20){
  N$mean[i] <- mean(dm[,i])
  N$sd[i] <- sd(dm[,i])
  N$SStot[i] <- sum((mean(dm[,i])-dm[,i])^2)
  rownames(N)[i] <- names(dm)[i]
}
N <- data.frame(name=rownames(N),mean=N$mean,sd=N$sd, SStot=N$SStot)


# normalization
n <- c(4:10,15:27)
D <- matrix(NA,nrow = 320,ncol = length(n), dimnames = list(NULL,names(d)[n]))
for(i in 1:ncol(D)){
  D[,i]<- (d[,n[i]]-mean(d[,n[i]]))/sd(d[,n[i]])
}
dimnames(D)
D <- D[,c(1,6,2:5,7:20)]


############### FITTING MODEL ######################

#### Third Model####
third_model <- '
# measurement model
thick.Ar =~ 1*thick.A
oc.Ar =~ 1*oc.A
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
esp.Ar =~ 1*esp.A
esp.Br =~ 1*esp.B
btr =~ 1*bt

# structural model
thick.Ar ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd

oc.Ar ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi +
            esp.Br + esp.Ar + btr + thick.Ar
tb.Ar ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            oc.Ar + btr
sat.Ar ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river + 
            tb.Ar + oc.Ar                            
esp.Ar ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            esp.Br 
esp.Br ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river 

btr ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf +
            esp.Br + esp.Ar

# measurement error
thick.A ~~  0.25*thick.A
oc.A ~~     0.20*oc.A
tb.A ~~     0.20*tb.A
sat.A ~~    0.20*sat.A
esp.A ~~    0.20*esp.A
esp.B ~~    0.10*esp.B
bt  ~~      0.25*bt
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))
D <- as.data.frame(D)
fit3 <- lavaan(model = third_model, data = D,
               model.type = "sem", meanstructure = "default",
               int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default",
               orthogonal = FALSE, std.lv = FALSE, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default",
               constraints = "", estimator = "ML",
               zero.cell.warn = TRUE, start = "default")
# varTable(fit3)
# summary(fit3, standardized=F, modindices = F, fit.measures=T) 
# inspect(fit3,"std.lv") # standardized model parameters
# inspect(fit3,"partable") #Observed sample statistics
# inspect(fit3, "cov.lv") #The model-implied covariance matrix of the observed variables
# inspect(fit3,"start") # starting values for all model parameters
# inspect(fit3,"rsquare") #R-square value for all endogenous variables
full_model <- lavaanify(third_model)
# modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
# modi<-modi[modi$mi>3 & !is.na(modi$mi),]
# modi
write.csv(full_model, "full_model.csv")

#################### PREDICTION #########################

################# External Drivers ##############
### rasters of external drivers, conversion to data frame and standardization
# install.packages(c("raster", "maptools", "sp", "rgdal"))
#install.packages("maptools")
library(raster)
library(maptools)
library(sp)
library(rgdal)

setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling")

# sdat files (dem) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "river", "wdist","maxc","mrvbf","slope","twi","vdchn","water") 
pred <- read.csv("mask_231m2.csv")


# tif files (modis)
files_m <- list.files(pattern=".tif$")[1:4]
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
setwd("/mnt/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
#P <-read.csv("calib.data-2.1.csv")[,-1]
for(i in 8:20){
  pred[,i]<- (pred[,i]- N$mean[i])/N$sd[i]
}
# SAVING pred
write.csv(pred, "pred.csv")
# LOAD STANDARDISED PREDICTORS 
pred <- read.csv("pred.csv")[,-1]

# for(i in n){
#   boxplot(c(d[i],pred[i]))
# }


################# Matrices ##############
##setting up matrices
B <- inspect(fit3,"est")$beta[1:7,1:7] #matrix of coeff. latent state variables
I <-diag(nrow=7, ncol=7) #Identity matrix
A <- inspect(fit3,"est")$beta[1:7,8:19] #matrix of coeff of external drivers
V <- inspect(fit3,"est")$psi[1:7,1:7] #matrix of predicted error variance
Th <- inspect(fit3, "est")$theta[1:7,1:7] # matrix of measurement error
IB<-solve(I-B)

################# Running Prediction ###########################

pred <- pred[,8:22] # only external drivers  + X,Y
pred <- as.matrix(pred) # conversion to matrix to improve processing speed
pred <- cbind(pred,matrix(NA,nrow=dim(pred)[1],ncol = 7)) #add colums for predicted values
dimnames(pred)[[2]][16:22] <-c("thick.Ar","oc.Ar","tb.Ar","sat.Ar","esp.Ar","esp.Br","btr") #names of soil properties

# (IB%*%A%*%p) product of matrices per pixel (equation 4 paper)
for(i in seq_along(pred[,1])) {
  p=matrix(pred[i,c(1,3,5,8,5,2,6,4,12,13,10,11)],nrow=12,ncol=1)
  pred[i,16:22]=t(IB%*%A%*%p) # key equation
}

t(dimnames(pred)[[2]])
pred <- pred[,14:22] #remove external drivers

# from unstandardize soil properties
M <- N[c(1,6,2:5,7),] 
for(i in 3:9){
  pred[,i]<- pred[,i]*M$sd[i-2]+M$mean[i-2]
}


##### Estimation prediction interval withd ## No stapial #######
Var.n <-IB%*%V%*%t(IB)+Th # diagonal vaues are variance error
Var<-diag(Var.n)
CI <- 1.64 * (Var ^ 0.5)

CI.r<- CI * M$sd
names(CI.r)<-M$name
CI.r[8:9] <- NA
CI.r[9] <- CI.r[7]
CI.r[7] <- CI.r[6]
names(CI.r) <- c(names(CI.r)[1:4], "ll.esp.A", "ul.esp.A", "ll.esp.B", "ul.esp.B", "bt")
CI.r[6] <- 10 ^(CI.r[5] + (Var[5] * M$sd[5]^2) * 0.5)
CI.r[5] <- 10 ^(CI.r[5] - (Var[5] * M$sd[5]^2) * 0.5)
CI.r[8] <- 10 ^(CI.r[7] + (Var[6] * M$sd[6]^2) * 0.5)
CI.r[7] <- 10 ^(CI.r[7] - (Var[6] * M$sd[6]^2) * 0.5)
CI.r
write.table(file = "PI.csv",x = CI.r, sep = "\t")
write.csv(Var,"SEM.variance.error.csv")
write.csv(M,"N.march2016.csv")
# from log10(ESP) to ESP
# pred[,8]<- 10^(pred[,8]+(Var[6]*M$sd[6]^2)*0.5)
# pred[,9]<- 10^(pred[,9]+(Var[7]*M$sd[7]^2)*0.5)
# print(summary(pred))
# 
# prediction.mean<-apply(pred[,3:9],MARGIN = 2,FUN = mean)
# write.table(file = "prediction.mean.csv",x = prediction.mean, sep = ",")

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
write.csv(s,"result.statistics.csv")


raster::NAvalue(rp)<--99999
writeRaster(x = rp[[2:8]],filename ="oktober.tif", overwrite=T,bylayer=TRUE,suffix=names(rp)[2:8])


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

