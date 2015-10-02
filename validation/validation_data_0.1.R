
rm(list=ls())
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data")

# load data
hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("Lab_data.csv")

# clean and transform C oxidable to OC
lab<- lab[c(-367,-366),]
lab$C_Ox<-lab$C_Ox*1.3

# Horizons A1 and A2, same site, two analysis
hor$num_lab[hor$num_lab==64561 & !is.na(hor$num_lab) ][2]<-64562
lab$labid[lab$labid==64561& !is.na(lab$labid)][2] <-64562

repl <- hor[!is.na(hor$num_lab_r),]
#--- MEASUREMENT ERROR AS STANDARD DEVIATION OF THE ERROR
replic <- unique(repl[,c(23,27)])
replic <- merge(x=replic,y=lab, by.x="num_lab", by.y="labid",all.x=T)
replic <- merge(x=replic,y=lab, by.x="num_lab_r", by.y="labid",all.x=T)
replic$pHKCl.x <- as.numeric(replic$pHKCl.x)
replic$pHKCl.y <- as.numeric(replic$pHKCl.y)
replic$C_Ox.x <- replic$C_Ox.x*1.3
replic$C_Ox.y <- replic$C_Ox.y*1.3
#fix(names) # as.data.frame(.Primitive("names"))
for(i in 3:18){
print(paste("SD",names(replic)[i],
            round(sqrt(var(replic[,i]-replic[,i+16], na.rm=T)),3),
            "Mean",round(mean(replic[,i]-replic[,i+16], na.rm=T),3),
            sep= ": "))
}
sqrt(0.5*(0.283^2))

#---
## 
#-------------------------------------#
# normalization of horizon names
as.data.frame(names(hor))
hor <- hor[,c(2:5,23)]
hor$hor <- NA
as.data.frame(t(table(hor$horizonte)))
hor$hor[grep("^A$",hor$horizonte)] <- "A"
hor$hor[grep("^a$",hor$horizonte)] <- "A"
hor$hor[grep("A1",hor$horizonte)] <- "A"
hor$hor[grep("A2",hor$horizonte)] <- "A"
hor$hor[grep("^2A",hor$horizonte)] <- "A"
hor$hor[grep("^An",hor$horizonte)] <- "A"

hor$hor[grep("E",hor$horizonte)] <- "E"

hor$hor[grep("AC",hor$horizonte)] <- "AC"

hor$hor[grep("AB",hor$horizonte)] <- "AB|BA"
hor$hor[grep("BA",hor$horizonte)] <- "AB|BA"
hor$hor[grep("Ba",hor$horizonte)] <- "AB|BA"
hor$hor[grep("B/A",hor$horizonte)] <- "AB|BA"

hor$hor[grep("^Bt",hor$horizonte)] <- "B"
hor$hor[grep("^2Bt",hor$horizonte)] <- "B"

hor$hor[grep("^BC",hor$horizonte)] <- "BC"
hor$hor[grep("Bc",hor$horizonte)] <- "BC"
hor$hor[grep("^2BC",hor$horizonte)] <- "BC"

hor$hor[grep("^C",hor$horizonte)] <- "C"
hor$hor[grep("^c",hor$horizonte)] <- "C"
hor$hor[grep("^2C",hor$horizonte)] <- "C"

as.data.frame(t(table(hor$horizonte[is.na(hor$hor)])))
#-------------------------------------#
# add coordenates to hor
hor$sitio[hor$sitio=="utic2-112"]<-"udic2-112"
hor.xy<-merge(hor,site[,c(2,5,6)], all.x=T)

# thickness of standarized horizons
library(plyr)
hor.xy <- hor.xy[!is.na(hor.xy$hor),]
hor.xy$sitio.hor <- paste(hor.xy$sitio,hor.xy$hor,sep=".")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, mintop=min(prof_s))[,c(1,2)], by= "sitio.hor")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, maxbot=max(prof_i))[,c(1,2)], by= "sitio.hor")
# number of sites with A horizons
length(unique(hor.xy$sitio.hor[hor.xy$hor=="A"]))
# Which are the top horizons
#table(hor.xy$hor[hor.xy$mintop==0])
#Thickness
hor.xy$thick<-hor.xy$maxbot-hor.xy$mintop

# load sitios strata
library(maptools)
sitios_strata <- readShapeSpatial("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data/sitios_strata.shp")
sitios_strata<- as.data.frame(sitios_strata)
hor.xy <- merge(hor.xy,sitios_strata, by="sitio", all.x=T)

#merge xy + horizons + lab data
hor.lab <- merge(hor.xy,lab, by.x="num_lab", by.y = "labid", all.x = T)
#correct strata 
hor.lab$strata[hor.lab$sitio=="2acuic1"] <- "acuic1"
hor.lab$strata[hor.lab$sitio=="co1-21"] <- "Co2"


#---------------------------------------#
#        Thickness top horizon          #
#---------------------------------------#
# SP is Soil Property
SP.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,12:15)])
SP.val <- SP.val[complete.cases(SP.val),]

library(sp)
library(rgdal)

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_thick.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")


par(pty="s")
par(mfrow = c(3, 2))
plot(SP.val@data$predicted~SP.val@data$measured, main="Thickness top horizon", xlab= "measured (cm)",
     ylab = "predicted (cm)", col = "dark red",xlim=c(10, 40), ylim=c(10, 40))
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")

# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")

library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ MSE (mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report<- data.frame(Soil_property = NA, ME=NA, RMSE= NA)
report[1,1:3]<- c("Thick.A",ME,RMSE)

#----------------------------------------#
#        Organic Carbon A horizon        #
#----------------------------------------#
as.data.frame(names(hor.lab))
SP.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,24,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]
#define crs
# wgs84 <- CRS("+init=epsg:4326")
# posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_oc.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
plot(SP.val@data$predicted~SP.val@data$measured, main="Organic carbon top horizon", xlab= "measured (%)",
     ylab = "predicted (%)", col = "dark red",xlim=c(1, 2.5), ylim=c(1, 2.5))
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ MSE (mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[2,1:3]<- c("OC.A",ME,RMSE)


#---------------------------------------#
#        Total Bases top horizon        #
#---------------------------------------#
as.data.frame(names(hor.lab))
SP.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,31,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]
#define crs
 wgs84 <- CRS("+init=epsg:4326")
 posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_tb.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
lim<-c(min(SP.val@data$measured)+0.5*sd(SP.val@data$measured),max(SP.val@data$measured)-0.5*sd(SP.val@data$measured))
plot(SP.val@data$predicted~SP.val@data$measured, main="Total bases top horizon", xlab= "measured (cmol+/kg)",
     ylab = "predicted (cmol+/kg)", col = "dark red",xlim=lim, ylim=lim)
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[3,1:3]<- c("TB.A",ME,RMSE)

#-----------------------------------------#
#        Base Saturation top horizon        #
#-----------------------------------------#
as.data.frame(names(hor.lab))
SP.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,32,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_sat.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
lim<-c(min(SP.val@data$measured)+0.5*sd(SP.val@data$measured),max(SP.val@data$measured)-0.5*sd(SP.val@data$measured))
plot(SP.val@data$predicted~SP.val@data$measured, main="Base saturation top horizon", xlab= "measured (%)",
     ylab = "predicted (%)", col = "dark red",xlim=lim, ylim=lim)
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[4,1:3]<- c("Sat.A",ME,RMSE)


#-------------------------------#
#        ESP top horizon        #
#-------------------------------#
as.data.frame(names(hor.lab))
SP.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,33,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_esp.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
lim<-c(min(SP.val@data$measured)+0*sd(SP.val@data$measured),max(SP.val@data$measured)-0*sd(SP.val@data$measured))
plot(SP.val@data$predicted~SP.val@data$measured, main="ESP top horizon", xlab= "measured (%)",
     ylab = "predicted (%)", col = "dark red",xlim=lim, ylim=lim)
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[5,1:3]<- c("ESP.A",ME,RMSE)


#-----------------------------#
#        ESP B horizon        #
#-----------------------------#
as.data.frame(names(hor.lab))
SP.val <- unique(hor.lab[hor.lab$hor=="B",c(7:9,33,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_esp.Br.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data)+1] <- over(x=SP.val,y=SP.pred)
names(SP.val@data)<-c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
lim<-c(min(SP.val@data$measured)+0*sd(SP.val@data$measured),max(SP.val@data$measured)-0*sd(SP.val@data$measured))
plot(SP.val@data$predicted~SP.val@data$measured, main="ESP B horizon", xlab= "measured (%)",
     ylab = "predicted (%)", col = "dark red",xlim=lim, ylim=lim)
abline(0,1)
abline(lm(SP.val@data$predicted~SP.val@data$measured), col="red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <-SP.val$measured-SP.val$predicted
SP.val$residuals.sq <-(SP.val$measured-SP.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(SP.val$strata))[,1]),
            n=as.data.frame(table(SP.val$strata))[,2],
            h=1:12)
SP.val<- merge(SP.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals))
H  <- 12
Nh <- as.data.frame(summarise(X,n()))[,2]
zh <- as.data.frame(summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah^2)*summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll<-zSt-qt(0.975,N-1)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,N-1)* sqrt(VzSt)

ME<- paste(round(ll,3),"<",round(zSt,3),"<",round(ul,3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah*zh.s)
paste("MSE =", round(zSt.s,3))
# variance of zSt
VzSt.s <- sum((ah^2)*summarise(X, var(residuals.sq))[2])
VzSt.s
# 95% confidence using X-square distribution
# # lowwer limit
ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # upper limit
ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[6,1:3]<- c("ESP.B",ME,RMSE)

write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report2.csv")



############# comparison between datasets

library(dplyr)
as.data.frame(names(hor.lab))
# validation dataset
val <- hor.lab[,c(2,4,7,12,32,24,31,33)]
val$hor<- as.factor(val$hor)
#val<-group_by(val,sitio,hor, add=T)
val.A<-group_by(val[val$hor=="A",], sitio)
val.A<-summarise(val.A,H=first(horizonte),THICK=mean(thick),SAT=mean(Sat),
                 CO=mean(C_Ox),TB=mean(Total_bases),ESP.A=mean(ESP))
val.B<-group_by(val[val$hor=="B",], sitio)
val.B<-summarise(val.B,ESP.B=mean(ESP))
val<-merge(val.A,val.B,by="sitio", all=T)

val$TB[is.na(val$TB)&!is.na(val$THICK)] <- 20
val$SAT[is.na(val$SAT)&!is.na(val$THICK)] <- 100
val$ESP.A[is.na(val$ESP.A)&!is.na(val$THICK)] <- 16.3
val$ESP.B[is.na(val$ESP.B)&!is.na(val$THICK)] <- 30.9
names(val)
val<- val[,c(3,6,4,7,8,5)]
# calibration dataset
cal<-read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/calib.data-2.1.csv")[,-1]
as.data.frame(names(cal))
cal<- cal[,c(4:9)]
names(cal)<- c("Thick.A (cm)","TB.A (cmol+/kg)","Sat.A (%)","log(ESP.A %)","log(ESP.B %)","OC.A (%)")
names(val)<- c("Thick.A (cm)","TB.A (cmol+/kg)","Sat.A (%)","log(ESP.A %)","log(ESP.B %)","OC.A (%)")
cal$`log(ESP.A %)`<-log(cal$`log(ESP.A %)`)
cal$`log(ESP.B %)`<-log(cal$`log(ESP.B %)`)
val$`log(ESP.A %)`<-log(val$`log(ESP.A %)`)
val$`log(ESP.B %)`<-log(val$`log(ESP.B %)`)

library(reshape2)
library(ggplot2)
c<- data.frame(value=NA,L1=NA,V3=NA)
for(i in 1:6){
  a<-list(cal[,i],val[,i])
  b<- melt(a)
  b<- b[!is.na(b$value),]
  b[,3]<- names(cal)[i]
  c<- rbind(c,as.data.frame(b))
}
c$L1[c$L1==1] <- "calibration"
c$L1[c$L1==2] <- "validation"
c<-c[!(c$values<60 & c$Property== "Sat.A (%)"),]
names(c)<- c("values","Dataset","Property")
c<-c[complete.cases(c),]
c<-group_by(c,Property)

ggplot() +  facet_wrap(~Property,scales = 'free_y') + 
  geom_jitter(data = c, mapping = aes(x=Dataset, y=values, color=Dataset), alpha = I(1/4))+
  geom_boxplot(data=c, mapping=aes(x=Dataset, y=values, color=Dataset), alpha = I(1/2))