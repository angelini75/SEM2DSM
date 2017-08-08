setwd("~/GSIF2016/SPcomp")

library(raster)
library(GSIF)
library(aqp)
library(plyr)
library(rgdal)
library(sp)
library(gstat)

# random seed
set.seed(180515)

################################################################################
### Download calibration points and covariates ####
# 
# ## get covariates from web  ## Only once!
# download.file("http://gsif.isric.org/zipped/NL250m_covs.zip", "NL250m_covs.zip")
# library(R.utils)
# unzip("NL250m_covs.zip")

# get point data from gdocs 
# library(RCurl)
# curl <- getCurlHandle()
# options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
# curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)
# cat(getURL("https://docs.google.com/spreadsheets/d/1J1aD6N15R4E_ZIfX41szslnCBpa-RowebAx_pXtfA4s/pub?gid=0&single=true&output=csv"), file="data.csv") 
# input <- read.csv("data.csv")
# str(input)
# rm(curl)


################################################################################
### Point data preparation ####

# backup ...
# input.backup <- input
# write.csv(input, "input_backup.csv")
# restore
# input <- input.backup
input <- read.csv("input_backup.csv")

#Clean litter horizons
input <- input[!input$UHDICM<0,]

# Unique ID
input$LOCID <- as.factor(paste("ID", input$LONWGS84, input$LATWGS84, sep="_"))

# Fill NAs in BLD 
input$BLD[is.na(input$BLD)] <- 1.378


temp <- OCSKGM(ORCDRC = input$ORCDRC, BLD = input$BLD*1000,
               HSIZE = input$LHDICM-input$UHDICM)
input$OCSKGM <- temp
input$meaERROR <- attr(temp,"measurementError")


#Soil profile collection
input.aqp <- input
depths(input.aqp) <- LOCID ~ UHDICM + LHDICM
ORCDRC.spline <- mpspline(input.aqp, var.name="ORCDRC")
BLD.spline <- mpspline(input.aqp, var.name="BLD")
#TODO: add 12 profiles with only 1 horizon!
## partir en 2

ORC.spl  <- ORCDRC.spline$var.1cm
ORC.spl <- ORC.spl[1:100,]
BLD.spl  <- BLD.spline$var.1cm
BLD.spl <- BLD.spl[1:100,]

input.1cm <- data.frame(ID = rep(ORCDRC.spline$idcol, each=100), 
                       ORC = as.vector(ORC.spl),
                       BLD = as.vector(BLD.spl)*1000)

temp <- OCSKGM(ORCDRC = input.1cm$ORC, BLD = input.1cm$BLD,
              HSIZE = 1)
input.1cm$OCSKGM <- temp
input.1cm$meaERROR <- attr(temp,"measurementError")

# Generate 2d points
coords <- matrix(unlist(strsplit(ORCDRC.spline$idcol, split = "_")), 
                 ncol=3 ,byrow = TRUE)
input.2d <- data.frame(ID=ORCDRC.spline$idcol, LONWGS84=as.numeric(coords[,2]),
                       LATWGS84=as.numeric(coords[,3]), 
                       OCSKGM=as.vector(tapply(input.1cm$OCSKGM, 
                                               INDEX = input.1cm$ID,
                                     FUN=sum, na.rm=TRUE)))


### Descriptive analysis

hist(input.2d$OCSKGM)
summary(input.2d$OCSKGM)


################################################################################
#### Covariates preparation ####

## MODIS stacks
grd.lst <- list.files(pattern="(^M_liste)")
stack250 <- stack(grd.lst, "SRTM_DEM250m.tif")
stack250 <- projectRaster(from = stack250, crs=CRS("+init=epsg:28992"))

## others
grd.lst_100 <- list.files(pattern="tif$")
grd.lst_100 <- grd.lst_100[grep("^(?!M_list)", grd.lst_100, perl=TRUE)] 
grd.lst_100 <- grd.lst_100[grep("^(?!SRTM)", grd.lst_100, perl=TRUE)] 
stack100 <- stack(grd.lst_100)
stack100




################################################################################
#### Calculate DEM derivates ####
## By: tom.hengl@isric.org

## Derive some standard DEM variables of interest for soil mapping:
saga_DEM_derivatives <- function(INPUT, MASK=NULL, 
                          sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM")){
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_hydrology 15 -DEM=\"', INPUT, '\" -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) )
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    suppressWarnings( system(paste0(saga_cmd, ' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=11' ) ) )
  }
}

## test it and plot all results:
writeRaster(stack250[["SRTM_DEM250m"]], "SRTM250.sdat", "SAGA")
saga_cmd = "/usr/bin/saga_cmd"
#saga_DEM_derivatives("SRTM250.sgrd")  ## Run only once!! time consuming
dem.lst <- list.files(pattern=glob2rx("^SRTM250_*.sdat"))
dem_derivates <- stack(dem.lst)


################################################################################
#### Add covariates to data ####

# Convert to spatial points df and project
coordinates(input.2d) <- ~ LONWGS84 + LATWGS84
input.2d@proj4string <- CRS("+init=epsg:4326")
input.2d <- spTransform(input.2d, CRSobj = CRS("+init=epsg:28992"))

input.2d@data <- cbind(input.2d@data, extract(stack250, input.2d),
                       extract(dem_derivates, input.2d),
                       extract(stack100, input.2d))

#Remove points without covariate data
input.2d@coords <- input.2d@coords[!is.na(input.2d@data$soilmap),] #TODO: improve
input.2d@data <- na.omit(input.2d@data)


## Interpolation grid
halfres <- res(stack250)[1]/2
grid <- expand.grid(x=seq(from=xmin(stack250)+halfres, to=xmax(stack250)-halfres, 
                          by=res(stack250)[1]), 
                    y=seq(from=ymin(stack250)+halfres, to=ymax(stack250)-halfres, 
                          by=res(stack250)[2]))
coordinates(grid) <- ~ x+y
gridded(grid) <- TRUE

## create table for models
Models <- data.frame(model=character(), ME=numeric(), MSPE=numeric(), 
                     RMSE=numeric(),corrObsPred=numeric(), corrObsRes=numeric())

################################################################################
#### Ordinary Kriging ####
model.ok <- gstat(formula = OCSKGM~1, data = input.2d)

model.ok.v <- variogram(model.ok)
plot(model.ok.v, plot.nu=FALSE)
plot(model.ok.v, plot.nu=TRUE)

# define initial semivariogram model
model.ok.vm <- vgm(nugget = 6000, psill = 7000, range = 15000, model = "Exp")
plot(model.ok.v, model.ok.vm)
# fit semivariogram model
model.ok.vm <- fit.variogram(model.ok.v, model.ok.vm, fit.method=7)
plot(model.ok.v, model.ok.vm)
model.ok.vm

# OKC.krige <- krige(formula = OCSKGM~1, locations = input.2d, 
#                    newdata = grid, model = model.ok.vm)

OKC.ok.cross <- krige.cv(OCSKGM~1, input.2d , model.ok.vm, verbose=T)

summary(OKC.ok.cross)
plot(OKC.ok.cross$var1.pred, OKC.ok.cross$observed)
# mean error, ideally 0:
ME <- mean(OKC.ok.cross$residual, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean(OKC.ok.cross$residual^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean(OKC.ok.cross$residual^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(OKC.ok.cross$observed, OKC.ok.cross$var1.pred)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(OKC.ok.cross$var1.pred, OKC.ok.cross$residual)

Models <- rbind(Models, cbind("OK", ME, MSPE, RMSE, corrObsPred, corrObsRes))
Models

save(model.ok.vm, file = "M_OK.RData")

################################################################################
#### Kriging with external drift ####

model.ked.lm <- lm(OCSKGM ~ . - ID, data=input.2d@data)
summary(model.ked.lm)
model.ked.slm <- step(model.ked.lm)
summary(model.ked.slm)

# Multicollinearity test
library(car)
sqrt(vif(model.ked.slm)) > 2 ## TRUE is collineal

# model.ked.slm.1 <- update(model.ked.slm, . ~ . - M_listeJulAug)
# summary(model.ked.slm.1)
# model.ked.slm.2 <- update(model.ked.slm, . ~ . - M_listeSepOct)
# summary(model.ked.slm.2)

trend <- as.formula(model.ked.slm$call)

dependend_var.v <- variogram(trend, input.2d, cutoff=50000,width=1000)
plot(dependend_var.v)
dependend_var.ovgm <- fit.variogram(dependend_var.v, vgm(psill=7000, model="Exp", 
                                                         range=2000, nugget=0))   
plot(dependend_var.v, dependend_var.ovgm, plot.nu=T)
dependend_var.ovgm

OKC.ked.cross <- krige.cv(trend, input.2d , verbose=T, dependend_var.ovgm)

summary(OKC.ked.cross)
plot(OKC.ked.cross$var1.pred, OKC.ked.cross$observed)
# mean error, ideally 0:
ME <- mean(OKC.ked.cross$residual, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean(OKC.ked.cross$residual^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean(OKC.ked.cross$residual^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(OKC.ked.cross$observed, OKC.ked.cross$var1.pred)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(OKC.ked.cross$var1.pred, OKC.ked.cross$residual)

Models <- rbind(Models, cbind("KED", ME, MSPE, RMSE, corrObsPred, corrObsRes))
Models

save(model.ked.slm, file = "M_KED.RData")

################################################################################
#### Ordinary Kriging log ####
model.oklog <- gstat(formula = log(OCSKGM)~1, data = input.2d)

model.oklog.v <- variogram(model.ok)
plot(model.oklog.v, plot.nu=FALSE)
plot(model.oklog.v, plot.nu=TRUE)

# define initial semivariogram model
model.oklog.vm <- vgm(nugget = 6000, psill = 7000, range = 15000, model = "Exp")
plot(model.oklog.v, model.oklog.vm)
# fit semivariogram model
model.oklog.vm <- fit.variogram(model.oklog.v, model.oklog.vm, fit.method=7)
plot(model.oklog.v, model.oklog.vm)
model.oklog.vm

# OKC.krige <- krige(formula = OCSKGM~1, locations = input.2d, 
#                    newdata = grid, model = model.ok.vm)

OKC.oklog.cross <- krige.cv(log(OCSKGM)~1, input.2d , model.oklog.vm)

summary(OKC.oklog.cross)
plot(OKC.oklog.cross$var1.pred, OKC.oklog.cross$observed)
# mean error, ideally 0:
ME <- mean(OKC.oklog.cross$residual, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean(OKC.oklog.cross$residual^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean(OKC.oklog.cross$residual^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(OKC.oklog.cross$observed, OKC.oklog.cross$var1.pred)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(OKC.oklog.cross$var1.pred, OKC.oklog.cross$residual)

Models <- rbind(Models, cbind("OKlog", ME, MSPE, RMSE, corrObsPred, corrObsRes))
Models

save(model.oklog.vm, file = "M_OKlog.RData")

################################################################################
#### KED log ####

model.kedlog.lm <- lm(log(OCSKGM) ~ . - ID, data=input.2d@data)
summary(model.kedlog.lm)
model.kedlog.slm <- step(model.kedlog.lm)
summary(model.kedlog.slm)

# Multicollinearity test
sqrt(vif(model.kedlog.slm)) > 2 ## TRUE is collineal

trend <- as.formula(model.kedlog.slm$call)

dependend_var.v <- variogram(trend, input.2d, cutoff=100000, width=2000)
plot(dependend_var.v)
dependend_var.ovgm <- fit.variogram(dependend_var.v, vgm(psill=0.6, model="Exp", 
                                                         range=3000, nugget=0))   
plot(dependend_var.v, dependend_var.ovgm, plot.nu=T)
dependend_var.ovgm

OKC.kedlog.cross <- krige.cv(trend, input.2d , dependend_var.ovgm, verbose=T)

summary(OKC.kedlog.cross)
plot(OKC.kedlog.cross$var1.pred, OKC.kedlog.cross$observed)
# mean error, ideally 0:
ME <- mean(OKC.kedlog.cross$residual, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean(OKC.kedlog.cross$residual^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean(OKC.kedlog.cross$residual^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(OKC.kedlog.cross$observed, OKC.kedlog.cross$var1.pred)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(OKC.kedlog.cross$var1.pred, OKC.kedlog.cross$residual)

Models <- rbind(Models, cbind("KEDlog", ME, MSPE, RMSE, corrObsPred, corrObsRes))
Models

save(model.kedlog.slm, file = "M_KEDlog.RData")


################################################################################
#### random forest base ####
library(randomForest)
model.rf <- randomForest(formula=(OCSKGM) ~ M_listeJanFeb + M_listeJulAug + M_listeMarApr + 
                           M_listeMayJun + M_listeNovDec + M_listeSepOct + SRTM_DEM250m + 
                           SRTM250_cprof + SRTM250_devmean + SRTM250_openn + SRTM250_openp + 
                           SRTM250_slope + SRTM250_twi + SRTM250_vbf + geomorfology + 
                           groundwater + landcover1970 + landcover1992 + landcover2004 + 
                           relativeElevation + soilmap, data=input.2d,
                         mtry=10, ntree=1000, nodesize=10, importance=TRUE,
                         keep.forest=T, keep.inbag=TRUE)
model.rf

plot(model.rf$predicted, model.rf$y)
# mean error, ideally 0:
ME <- mean(model.rf$predicted - model.rf$y, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean((model.rf$predicted - model.rf$y)^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean((model.rf$predicted - model.rf$y)^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(model.rf$predicted, model.rf$y)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(model.rf$predicted, model.rf$predicted - model.rf$y)

Models <- rbind(Models, cbind("rf_base", ME, MSPE, RMSE, corrObsPred, corrObsRes))
Models

save(model.rf, file = "M_randomFores.RData")

################################################################################
#### Neural networks ####

# http://stackoverflow.com/questions/32077931/specify-cross-validation-folds-with-caret

library(nnet)
library(caret)
library(datasets)
library(data.table)
library(e1071)

r <- 3  # number of repeats
k <- 5 # number of folds

# Create folds and repeats here#
input.nnet <- data.table(input.2d@data)
for (i in 1:r) {
  newcol <- paste('fold.num',i,sep='')
  input.nnet <- input.nnet[,eval(newcol):=sample(1:k, size=dim(input.nnet)[1], replace=TRUE)]
}

folds.list.out <- list()
folds.list <- list()
list.counter <- 1
for (y in 1:r) {
  newcol <- paste('fold.num', y, sep='')
  for (z in 1:k) {
    folds.list.out[[list.counter]] <- which(input.nnet[,newcol,with=FALSE]==z)
    folds.list[[list.counter]] <- which(input.nnet[,newcol,with=FALSE]!=z)
    list.counter <- list.counter + 1
  }
  input.nnet <- input.nnet[,!newcol,with=FALSE]
}

tune.grid <- expand.grid(
  size = c(4,6,8,10),
  decay = c(0,0.1,0.001,0.0001)
)

train.control <- trainControl(
  index=folds.list
  , indexOut=folds.list.out
  , verboseIter = T
  , returnData = T
  , savePredictions = T
)

input.nnet <- data.frame(input.nnet)

nnet.train <- train(
  form= OCSKGM ~ M_listeJanFeb + M_listeJulAug + M_listeMarApr +
    M_listeMayJun + M_listeNovDec + M_listeSepOct + SRTM_DEM250m +
    SRTM250_cprof + SRTM250_devmean + SRTM250_openn + SRTM250_openp +
    SRTM250_slope + SRTM250_twi + SRTM250_vbf + geomorfology +
    groundwater + landcover2004 +soilmap, data=input.nnet,
    method = "nnet",
    linout = TRUE,
    maxit=500,
    preProcess = c("center","scale"),
    metric = "RMSE",
    trControl = train.control,
    tuneGrid = tune.grid
)

nnet.train
plot(nnet.train)

plot(nnet.train$finalModel$fitted.values, input.nnet$OCSKGM)
# mean error, ideally 0:
ME <- mean(nnet.train$finalModel$residuals, na.rm=TRUE)
# Mean square prediction error, ideally small
MSPE <- mean((nnet.train$finalModel$residuals)^2, na.rm=TRUE)
# RMSE
RMSE <- sqrt(mean((nnet.train$finalModel$residuals)^2, na.rm=TRUE))
# correlation observed and predicted, ideally 1
corrObsPred <- cor(nnet.train$finalModel$fitted.values, input.nnet$OCSKGM)
# correlation predicted and residual, ideally 0
corrObsRes <- cor(nnet.train$finalModel$fitted.values, nnet.train$finalModel$residuals)

## Fix results table
Models$V1 <- as.character(Models$V1)
for(i in 2:6){
  Models[,i] <- as.numeric(as.character(Models[,i]))
}
Models <- rbind(Models, c("nnet", ME, MSPE, RMSE, 
                          as.numeric(corrObsPred), as.numeric(corrObsRes)))
Models
write.csv(Models, "Models.csv")

save(nnet.train, file = "M_NNET.RData")

################################################################################
#### Validation points ####

#Load validation data, convert to spatial and add covariates
valid <- read.csv("validationpoints.csv")

coordinates(valid) <- ~ LONWGS84 + LATWGS84
valid@proj4string <- CRS("+init=epsg:4326")
valid <- spTransform(valid, CRSobj = CRS("+init=epsg:28992"))

#Plot valid and calibration data
plot(valid, pch=20, cex=0.5, col="green")
points(input.2d, cex=0.3, col="red")

valid@data <- cbind(valid@data, extract(stack250, valid),
                       extract(dem_derivates, valid),
                       extract(stack100, valid))

# #Remove points without covariate data
# valid@coords <- valid@coords[!is.na(valid@data$soilmap),] 
# valid@data <- na.omit(valid@data)


## Predict using the best models
valid$nnet.0 <- predict(nnet.train, valid@data)
valid$RF.0 <- predict( model.rf, valid@data)

write.csv(valid, "results.csv")

################################################################################
#### Write results ####

results <- predict( model.rf, valid@data)
RMSE <- sqrt(mean((model.rf$predicted - model.rf$y)^2, na.rm=TRUE))

shapefile(input.2d, "calib.shp")
shapefile(valid, "valid.shp")

################################################################################
#### Generate a map of OCSKGM 0-100cm ####

#Resample covariates to same cellsize
covariates <- stack(stack250, dem_derivates, resample(stack100, stack250))