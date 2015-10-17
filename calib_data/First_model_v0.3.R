rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
original <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/OLD/calib.data-1.0.csv")
d <-read.csv("calib.data-2.1.csv")[,-1]

############### PRE-PROCESSING ################## 
names(d)
names(d)[c(5,6,9,10)]<- c("tb.A","sat.A", "oc.A","bt")
d<-d[!is.na(d$oc.A),]
d$sat.A[is.na(d$sat.A)] <- 100 # it is assumed 100% saturation when CaCO3 is present 
d$d.caco3[is.na(d$d.caco3)] <- 300 # it is assumed that CaCO3 is very deep when it is absent within the solum 
summary(d$tb.A[d$is.caco3==1], omit.na=T)
d$tb.A[is.na(d$tb.A)] <- 20

# plot(d$esp.A[d$is.caco3==1]~d$d.caco3[d$is.caco3==1], omit.na=T)
# boxplot(d$d.caco3[is.na(d$esp.A)])
summary(d$esp.A[d$d.caco3<15], omit.na=T)
d$esp.A[is.na(d$esp.A)] <- 16.3

# plot(d$esp.B[d$is.caco3==1]~d$esp.A[d$is.caco3==1])
summary(d$d.caco3[is.na(d$esp.B)])
summary(d$esp.B[d$d.caco3<30], omit.na=T)
d$esp.B[is.na(d$esp.B)] <- 30.9

nas<-d[!complete.cases(d),]

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

# hist(10^(d$esp.A),col = "lightblue")
# summary(d$thick.A)
##@## data normalization
N<- data.frame(mean = rep(0,20),sd=rep(0,20))

# save mean and sd

dm <- d[,c(4:10,15:27)]
for(i in 1:20){
  N$mean[i] <- mean(dm[,i])
  N$sd[i] <- sd(dm[,i])
  rownames(N)[i] <-names(dm)[i]
}

# normalization
n <- c(4:10,15:27)
for(i in n){
  d[,i]<- (d[,i]-mean(d[,i]))/sd(d[,i])
}

#boxplot(d[,c(4:10,15:27)])  


# step(lm(sat.A ~ dem+wdist+maxc+mrvbf+slope+twi+vdchn+lstm+lstsd+evim+evisd+river, d),direction = "both")
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd, data = d))
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd + tb.A + d.caco3 + river, data = d))
oc.fit<- lm(formula = oc.A ~ lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi, data = d)
summary(oc.fit)
x<-predict(oc.fit,data=d)

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

#### Third Model####
third_model <- '
# measurement model
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
btr =~ 1*bt
oc.Ar =~ 1*oc.A
thick.Ar =~ 1*thick.A
esp.Br =~ 1*esp.B
esp.Ar =~ 1*esp.A

# structural model
tb.Ar ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            oc.Ar + btr
sat.Ar ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river + 
            tb.Ar + oc.Ar                            
btr ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf +
            esp.Br + esp.Ar
oc.Ar ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi +
            esp.Br + esp.Ar + btr + thick.Ar
esp.Ar ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            esp.Br 
esp.Br ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river 
            
thick.Ar ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd



# measurement error
thick.A ~~  0.25*thick.A
tb.A ~~     0.20*tb.A
sat.A ~~    0.20*sat.A
bt  ~~      0.25*bt
oc.A ~~     0.20*oc.A
esp.B ~~    0.10*esp.B
esp.A ~~    0.20*esp.A
# thick.A ~~  0.25*thick.A
# tb.A ~~     0.1*tb.A
# sat.A ~~    0.1*sat.A
# bt  ~~      0.1*bt
# oc.A ~~     0.1*oc.A
# esp.B ~~    0.10*esp.B
# esp.A ~~    0.1*esp.A
# intercepts
# tb.Ar ~1
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))

fit3 <- lavaan(model = third_model, data = d,
               model.type = "sem", meanstructure = "default",
               int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default",
               orthogonal = FALSE, std.lv = FALSE, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default",
               constraints = "", estimator = "ML",
               zero.cell.warn = TRUE, start = "default")
varTable(fit3)
summary(fit3, standardized=F, modindices = F, fit.measures=T) 
inspect(fit3,"std.lv") # standardized model parameters
inspect(fit3,"partable") #Observed sample statistics
inspect(fit3, "cov.lv") #The model-implied covariance matrix of the observed variables
inspect(fit3,"start") # starting values for all model parameters
inspect(fit3,"rsquare") #R-square value for all endogenous variables
full_model <- lavaanify(third_model)
modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>3 & !is.na(modi$mi),]
modi

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


pred.df <- pred.df[complete.cases(pred.df),]
pred.df <- pred.df[,c(2:16)]

for(i in 1:14) {
  pred.df[,length(pred.df)+1] <- NA
  names(pred.df)[length(pred.df)] <- names(d)[i]
}  
pred.df <- pred.df[,c(16:29,1:15)]
pred <- pred.df
##@## data normalization
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")
D <-read.csv("calib.data-2.1.csv")[,-1]
n <- c(15:27)
for(i in n){
  pred[,i]<- (pred[,i]-mean(D[,i]))/sd(D[,i])
}
names(D)[15:27]
# for(i in n){
#   boxplot(c(d[i],pred[i]))
# }


################# Matrices ##############
##setting up matrices
B <- inspect(fit3,"est")$beta[1:7,1:7] #matrix of coeff. latent state variables
I <-diag(nrow=7, ncol=7) #Identity matrix
A <- inspect(fit3,"est")$beta[1:7,8:19] #matrix of coeff of external drivers
V <- inspect(fit3,"est")$psi[1:7,1:7] #matrix of predicted error variance
IB<-solve(I-B)

################# Running Prediction ###########################
N <- read.csv("N.csv")
library(utils)
pb = txtProgressBar(min = 0, max = length(pred[,1]), initial = 0, style = 3)

pred <- pred[,15:29] #o nly external drivers  + X,Y
pred <- as.matrix(pred) # conversion to matrix to improve processing speed
pred <- cbind(pred,matrix(NA,nrow=dim(pred)[1],ncol = 7)) #add colums for predicted values
dimnames(pred)[[2]][16:22] <-c("tb.Ar","sat.Ar","btr","oc.Ar","thick.Ar","esp.Br","esp.Ar") #names of soil properties

# (IB%*%A%*%p) product of matrices per pixel (equation 4 paper)
for(i in seq_along(pred[,1])) {
  p=matrix(pred[i,c(12,13,10,11,1,3,5,8,7,2,6,4)],nrow=12,ncol=1)
  pred[i,16:22]=t(IB%*%A%*%p) # key equation
  setTxtProgressBar(pb,i)
}

t(dimnames(pred)[[2]])
pred <- pred[,14:22] #remove external drivers

# from unstandardize soil properties
N <- N[c(2,3,7,6,1,5,4),] 
for(i in 3:9){
  pred[,i]<- pred[,i]*N[i-2,3] + N[i-2,2]
}

##### Estimation confidence interval ## No stapial #######
Var.n <-IB%*%V%*%t(IB) # diagonal vaues are variance error
Var<-diag(Var.n)
CI <- Var.n^(1/2)*1.64
CI.r <- matrix(0,nrow=7,ncol = 7)
for(i in 1:7){
  CI.r[i,i]<- CI[i,i]*N[i,3]
}
CI.r<-diag(CI.r)
names(CI.r)<-as.character(N[,1])
CI.r

# from log10(ESP) to ESP
pred[,8]<- 10^(pred[,8]+(Var[6]*N$sd[6]^2)*0.5)
pred[,9]<- 10^(pred[,9]+(Var[7]*N$sd[7]^2)*0.5)
print(summary(pred))


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
plot(r[[3]])
ADE<- readShapePoly("/media/marcos/L0135974_DATA/UserData/BaseARG/study area/ADE_MODIS.shp")
plot(ADE)
r<-mask(x = r,mask = ADE)
TWI <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/TWI250.sdat")
proj4string(TWI) <- posgar98
#rp[!rp %in% NA] <- 0
res(rp)
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

