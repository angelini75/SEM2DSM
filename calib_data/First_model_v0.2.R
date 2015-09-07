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


# intercepts
tb.Ar ~1
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
summary(fit3, standardized=F, modindices = F, fit.measures=F) 
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
files_m <- list.files(pattern=".tif")
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
for(i in 1:length(files)) {
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
for(i in 1:length(files_m)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header_m[i]
}  
#image(raster(("mod13q1_tot_mean.tif")))

# clean points out of study area
pred <- spTransform(pred, posgar98)
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


################# Prediction ##############
##setting up matrices
B <- inspect(fit3,"est")$beta[1:7,1:7] #matrix of coeff. latent state variables
I <-diag(nrow=7, ncol=7) #Identity matrix
A <- inspect(fit3,"est")$beta[1:7,8:19] #matrix of coeff of external drivers
V <- inspect(fit3,"est")$psi[1:7,1:7] #matrix of predicted error variance

save.image("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/prediction.RData")
## go to RStudio server to run next step
#################Prediction model###########################

N <- read.csv("N.csv")
library(utils)
pb = txtProgressBar(min = 0, max = length(pred[,1]), initial = 0, style = 3)
pred <- pred[,15:29]
for(i in 1:length(pred[,1])) {
  p=matrix(c(pred$evim[i], pred$evisd[i],
             pred$lstm[i], pred$lstsd[i], pred$dem[i], 
             pred$wdist[i], pred$mrvbf[i], pred$vdchn[i],
             pred$twi[i], pred$river[i], pred$slope[i],
             pred$maxc[i]),nrow=12,ncol=1)
  n=c(0,0,0,0,0,0,0)
  n=(solve(I-B))%*%((A%*%p)) # key equation
#   pred$tb.Ar[i] <-  n[1]
#   pred$sat.Ar[i] <-  n[2]
#   pred$btr[i] <-  n[3]
  pred$oc.Ar[i] <-  n[4]
#   pred$thick.Ar[i] <-  n[5]
#   pred$esp.Br[i] <-  n[6]
#   pred$esp.Ar[i] <-  n[7]
  setTxtProgressBar(pb,i)
}
# not to run
# for(i in 1:length(pred[,1])) {
#   p=matrix(c(pred$evim[i], pred$evisd[i],
#              pred$lstm[i], pred$lstsd[i], pred$dem[i], 
#              pred$wdist[i], pred$mrvbf[i], pred$vdchn[i],
#              pred$twi[i], pred$river[i], pred$slope[i],
#              pred$maxc[i]),nrow=12,ncol=1)
#   n=c(0,0,0,0,0,0,0)
#   n=(solve(I-B))%*%((A%*%p)) # key equation
#   pred$tb.Ar[i] <-  n[1]
#   pred$sat.Ar[i] <-  n[2]
#   pred$btr[i] <-  n[3]
#   pred$oc.Ar[i] <-  n[4]
#   pred$thick.Ar[i] <-  n[5]
#   pred$esp.Br[i] <-  n[6]
#   pred$esp.Ar[i] <-  n[7]
#   setTxtProgressBar(pb,i)
# }

# result of the prediction
#pred.bk <-read.csv("pred.bk.csv")
pred.bk <-pred

#pred <- pred.bk[,-1]
t(names(pred))
pred <- pred[,14:16] #X, Y and OC
#pred <- pred[,14:22]

# N <- N[c(2,3,7,6,1,5,4),]
N <- N[c(2,3,7,6,1,5,4),][4,] # mean and sd OC

# as.data.frame(names(D))
# as.data.frame(names(pred.bk))

# 
# for(i in 3:9){
#   pred[,i]<- pred[,i]*N[i-2,3] + N[i-2,2]
# }

for(i in 3:3){ # only for OC
  pred[,i]<- pred[,i]*N[i-2,3] + N[i-2,2]
}
save.image("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/prediction2.RData")
##### estimation confidence limit of OC
ci <- sqrt(V[4,4])*1.64
pred$oc.up <-pred$oc.Ar+ ci
pred$oc.lw <-pred$oc.Ar- ci
summary(pred[,3:5])

####rasterize results
library(sp)
library(raster)
pred.sp <- pred
coordinates(pred.sp) <- ~X+Y
# spplot(pred.sp)

y <- raster("mask_231m_posgar.tif")
proj4string(y) <- posgar98
proj4string(pred.sp) <- posgar98
#pred.sp <- spTransform(pred.sp, modis)
r <- rasterize(x = pred.sp,y = y,background= NA)
plot(r)

writeRaster(x = r,filename ="rusults.tif", overwrite=T,bylayer=TRUE,suffix=r@data@names)


#################Plotting results####################

#install.packages("semPlot")
library(semPlot)
semPaths(fit3,  "model","est",style ="lisrel")
semPaths (fit3, "std",sizeLat = 4, sizeInt2 = 2,sizeMan = 3, sizeLat2 = 2, sizeInt = 1,sizeMan2 = 1.5,
          edge.width = 2 , edge.label.bg=T,layout = "circle" , equalizeManifests =F, esize=2, asize=1, 
          structural = T, intercepts = F, residuals = T, thresholds = F,reorder = T)
