
rm(list = ls())
name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data")

# load data
hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("Lab_data_mod.csv")
lab$X0CaCO3 <- as.character(lab$X0CaCO3)
lab$X0CaCO3[lab$X0CaCO3 == ""] <- 0
lab$X0CaCO3[lab$X0CaCO3 == "*"] <- 1
lab$X0CaCO3[lab$X0CaCO3 == "**"] <- 2
lab$X0CaCO3[lab$X0CaCO3 == "***"] <- 3
lab$X0CaCO3 <- as.numeric(lab$X0CaCO3)
#hist(lab$pHKCl)

# clean and transform C oxidable to OC
lab <- lab[c(-367,-366), ]
lab$C_Ox <- lab$C_Ox * 1.3

# # Horizons A1 and A2, same site, two analysis
# #hor$num_lab[hor$num_lab==64561 & !is.na(hor$num_lab) ][2]<-64562
# #lab$labid[lab$labid==64561& !is.na(lab$labid)][2] <-64562
# hor$num_lab_r[hor$num_lab_r==64599& !is.na(hor$num_lab_r) ]<- 64559
lab[lab$labid == 65012,] <- lab[lab$labid == 64879,]
repl <- hor[!is.na(hor$num_lab_r),]
#--- MEASUREMENT ERROR AS STANDARD DEVIATION OF THE ERROR
replic <- unique(repl[,c(23,27)])
replic <- merge(x = replic,y = lab, by.x = "num_lab", by.y = "labid", all.x = T)
replic <- merge(x = replic,y = lab, by.x = "num_lab_r", by.y = "labid", all.x = T)
replic$pHKCl.x <- as.numeric(replic$pHKCl.x)
replic$pHKCl.y <- as.numeric(replic$pHKCl.y)
replic <- replic[replic$num_lab != 64553,]
replic <- replic[replic$num_lab != 65012,]


error <- replic[,1:18]
for (i in 3:18) {
  error[,i] <- replic[,i] - replic[,i + 16]
  print(paste(names(error)[i],mean(error[,i], na.rm = T) * 
                sqrt(length(error[complete.cases(error[,i]),i])) /
                sd(error[,i], na.rm = T)))
}
summary(error)
measurement.error <- data.frame(Property = NA, MEANe = NA, VARe = NA, SDe = NA, S.Var = NA)
for (i in 3:18) {
measurement.error[i - 2,] <- c(names(error)[i], mean(error[,i], na.rm = T),
            2 * var(error[,i], na.rm = T), sqrt(2 * var(error[,i], na.rm = T)),
            var(replic[,i], na.rm = T))
}
measurement.error[,2] <- as.numeric(measurement.error[,2])
measurement.error[,3] <- as.numeric(measurement.error[,3])
measurement.error[,4] <- as.numeric(measurement.error[,4])
measurement.error[,5] <- as.numeric(measurement.error[,5])
# comapre VARe/Var(property)
measurement.error$ratio <- measurement.error$VARe / measurement.error$S.Var
write.csv(measurement.error,"/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/measurement.error.csv")

e <- matrix(NA,nrow = 35,ncol = 17)
e <- as.data.frame(e)
e[,1] <- replic[,2]
names(e)[1] <- names(replic)[2]
for (i in 3:18) {
  e[,i - 1] <- replic[,i] - replic[,i + 16]
  names(e)[i - 1] <- names(replic)[i]
}
summary(e)
#hist(e$Total_bases.x)
E <- e[e$pHw.x > 1.5 | e$pHw.x < -1.5,]
E <- E[!is.na(E$Sat_w.x),]
summary(E)

m <- matrix(NA,nrow = 35,ncol = 17)
m <- as.data.frame(m)
m[,1] <- replic[,2]
names(m)[1] <- names(replic)[2]

for (i in 3:18) {
  m[,i - 1] <- (replic[,i] + replic[,i + 16]) / 2
  names(m)[i - 1] <- names(replic)[i]
}


##### replace horizons with duplo analysis (2 measurements per sample) by mean(analysis[i])
s <- m$num_lab
for (i in 1:length(s)) {
  lab[lab$labid == s[i],] <- as.vector(as.matrix(m[m$num_lab == s[i],]))
}

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
hor$sitio[hor$sitio == "utic2-112"] <- "udic2-112"
hor.xy <- merge(hor,site[,c(2,5,6)], all.x = T)

# thickness of standarized horizons
library(plyr)
hor.xy <- hor.xy[!is.na(hor.xy$hor),]
hor.xy$sitio.hor <- paste(hor.xy$sitio,hor.xy$hor,sep = ".")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, mintop = min(prof_s))[,c(1,2)], by = "sitio.hor")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, maxbot = max(prof_i))[,c(1,2)], by = "sitio.hor")
# number of sites with A horizons
length(unique(hor.xy$sitio.hor[hor.xy$hor == "A"]))
# Which are the top horizons
#table(hor.xy$hor[hor.xy$mintop==0])
#Thickness
hor.xy$thick <- hor.xy$maxbot - hor.xy$mintop

# load sitios strata
library(maptools)
sitios_strata <- readShapeSpatial("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data/sitios_strata.shp")
sitios_strata <- as.data.frame(sitios_strata)
hor.xy <- merge(hor.xy, sitios_strata, by = "sitio", all.x = T)

#merge xy + horizons + lab data
hor.lab <- merge(hor.xy, lab, by.x = "num_lab", by.y = "labid", all.x = T)
#correct strata 
hor.lab$strata[hor.lab$sitio == "2acuic1"] <- "acuic1"
#hor.lab$strata[hor.lab$sitio == "co1-21"] <- "Co2"

# statistics of validation data
as.data.frame(names(hor.lab))
# Thick.A
# OC.A
# TB.A
# Sat.A
# ESP.A
# ESP.B

d.A <- unique(hor.lab[hor.lab$hor == "A",c(12,24,31,32,33)])
d.B <- data.frame(ESP.B = unique(hor.lab[hor.lab$hor == "B",c(33)]))

d.stat <- matrix(data = NA,nrow = 6,ncol = 6,
                dimnames = list(c("Min", "Median", "Mean", "Max", "SD", "SS"),
                                c(names(d.A),names(d.B))))
d.stat[1,1:5] <- apply(X = d.A,FUN = min,2, na.rm = T)
d.stat[1,6] <- apply(X = d.B,FUN = min,2, na.rm = T)
d.stat[2,1:5] <- apply(X = d.A,FUN = median,2, na.rm = T)
d.stat[2,6] <- apply(X = d.B,FUN = median,2, na.rm = T)
d.stat[3,1:5] <- apply(X = d.A,FUN = mean,2, na.rm = T)
d.stat[3,6] <- apply(X = d.B,FUN = mean,2, na.rm = T)
d.stat[4,1:5] <- apply(X = d.A,FUN = max,2, na.rm = T)
d.stat[4,6] <- apply(X = d.B,FUN = max,2, na.rm = T)
d.stat[5,1:5] <- apply(X = d.A,FUN = sd,2, na.rm = T)
d.stat[5,6] <- apply(X = d.B,FUN = sd,2, na.rm = T)

for (i in 1:5) {
  d.stat[6,i] <- sum((mean(d.A[complete.cases(d.A[,i]),i]) -
                        d.A[complete.cases(d.A[,i]),i]) ^ 2 )
}
d.stat[6,6] <- sum((mean(d.B[complete.cases(d.B[,1]),1]) - d.B[complete.cases(d.B[,1]),1]) ^ 2)
dimnames(d.stat)[[2]] <- c("thick.A", "co.A", "tb.A", "sat.A", "esp.A", "esp.B")

#################### VALIDATION ####################
library(sp)
library(rgdal)

#-------------------------------------#
#        Thickness A horizon          #
#-------------------------------------#
# SP is Soil Property
SP.val <- unique(hor.lab[hor.lab$hor == "A",c(7:9,12:15)])
SP.val <- SP.val[complete.cases(SP.val),]

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_thick.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val) <- ~ longitud + latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

par(pty = "s")
par(mfrow = c(3, 2))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Thickness A horizon (cm)", 
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = c(10, 40), ylim = c(10, 40))
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")

# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted) ^ 2 
SS <- sum(SP.val$residuals.sq)

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "acuic3",])
library(dplyr)
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 1
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah*zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3), "<", round(zSt, 3), "<", round(ul, 3))

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

# fill report table
report <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA)
report[1,1:4] <- c("Thick.A",ME,RMSE,SS)

#----------------------------------------#
#        Organic Carbon A horizon        #
#----------------------------------------#
as.data.frame(names(hor.lab))
rm(SP.val)
rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(7:9,24,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_oc.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val) <- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor", "measured", "strata", "area", "percentage", "predicted")

# plot residuals
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Organic carbon A horizon (%)",
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = c(1, 2.5), ylim = c(1, 2.5))
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted) ^ 2
SS <- sum(SP.val$residuals.sq) 

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "acuic3",])
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "albic1",])

#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 2
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah/A
zSt <- sum(ah * zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975,N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975,N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3), "<", round(zSt, 3), "<", round(ul, 3))

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh*zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s,3))
paste("RMSE =", round(sqrt(zSt.s),3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

# fill report table
report[2,1:4] <- c("OC.A",ME,RMSE,SS)


#-------------------------------------#
#        Total Bases A horizon        #
#-------------------------------------#
as.data.frame(names(hor.lab))
rm(SP.val)
rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(7:9,31,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_tb.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

# plot residuals
lim <- c(min(SP.val@data$measured) + 0.5 * sd(SP.val@data$measured), max(SP.val@data$measured) -
           0.5 * sd(SP.val@data$measured))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Total bases A horizon (cmol+/kg)",
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted) ^ 2 
SS <- sum(SP.val$residuals.sq)

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "acuic3",])
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "albic1",])
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 2
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X, n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh * zh) / N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah / A
zSt <- sum(ah * zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3), "<", round(zSt, 3), "<", round(ul, 3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh * zh.s) / N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s,3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

# fill report table
report[3,1:4] <- c("TB.A", ME, RMSE, SS)

#-----------------------------------------#
#        Base Saturation A horizon        #
#-----------------------------------------#
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(7:9,32,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_sat.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val)<- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

# plot residuals
lim <- c(min(SP.val@data$measured) + 0.5 * sd(SP.val@data$measured), 
         max(SP.val@data$measured) - 0.5 * sd(SP.val@data$measured))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Base saturation A horizon (%)",
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured-SP.val$predicted) ^ 2 
SS <- sum(SP.val$residuals.sq)

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "acuic3",])
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "albic1",])

#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 2
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh * zh) / N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah / A
zSt <- sum(ah * zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3),"<" ,round(zSt, 3), "<", round(ul, 3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh * zh.s) / N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s,3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)
# # variance of zSt
# VzSt.s <- sum((ah^2)*dplyr::summarise(X, var(residuals.sq))[2])
# VzSt.s
# # 95% confidence using X-square distribution
# # # lowwer limit
# ll.s <- zSt.s-qchisq(p=0.025, df=N-1, ncp = 0, lower.tail = TRUE, log.p = FALSE)*sqrt(VzSt)
# # # upper limit
# ul.s <- zSt+qchisq(p=0.975, df=N-1, ncp = 0, lower.tail = T, log.p = FALSE)* sqrt(VzSt)
# RMSE <-paste(round(ll.s,1),"<",round(zSt.s,1),"<",round(ul.s,1))

# fill report table
report[4,1:4] <- c("Sat.A", ME, RMSE, SS)


#-----------------------------#
#        ESP A horizon        #
#-----------------------------#
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A",c(7:9,33,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_esp.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val) <- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

# plot residuals
lim <- c(min(SP.val@data$measured) + 0 * sd(SP.val@data$measured), 
         max(SP.val@data$measured) - 0 * sd(SP.val@data$measured))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "ESP A horizon (%)", xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted) ^ 2 
SS <- sum(SP.val$residuals.sq)

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "acuic3",])
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "albic1",])
#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 2
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh * zh) / N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah / A
zSt <- sum(ah * zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3), "<", round(zSt, 3), "<", round(ul, 3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh * zh.s) / N 
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s,3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

# fill report table
report[5,1:4] <- c("ESP.A", ME, RMSE, SS)

#-----------------------------#
#        ESP B horizon        #
#-----------------------------#
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "B",c(7:9,33,13:15)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/oktober_esp.Br.tif")
# thickness to spatial data frame
coordinates(SP.val) <- ~longitud+latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

# plot residuals
# par(pty="s")
# par(mfrow = c(2, 2))
lim <- c(min(SP.val@data$measured) + 0 * sd(SP.val@data$measured),
         max(SP.val@data$measured) - 0 * sd(SP.val@data$measured))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "ESP B horizon (%)", xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted) ^ 2 
SS <- sum(SP.val$residuals.sq)

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata = levels(SP.val$strata),
            n = as.data.frame(table(SP.val$strata))[,2],
            h = 1:12)
SP.val <- merge(SP.val, nh, by = "strata")
SP.val <- rbind(SP.val,SP.val[SP.val$strata == "Co2",])

#vignette("introduction", package = "dplyr")
X  <- group_by(SP.val, strata)
dplyr::summarise(X,n())
# Mean error of the area
N  <- as.numeric(length(X$residuals)) - 1
H  <- 12
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2]
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2]
sum(Nh * zh) / N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2]
A <- sum(Ah)
ah <- Ah / A
zSt <- sum(ah * zh)
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2])
VzSt
# 95% confidence
# lowwer limit
ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- paste(round(ll, 3), "<", round(zSt, 3), "<",round(ul, 3))

############################################ RMSE (root mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2]
sum(Nh * zh.s) / N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s)
paste("MSE =", round(zSt.s, 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

# fill report table
report[6,1:4] <- c("ESP.B",ME,RMSE,SS)

setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

report$R2 <- NA
for (i in 1:6) {
  report$R2[i] <- t(1 - (as.numeric(report[i,4]) / d.stat[6,i]))
}
write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report.validation.csv")






############# comparison between datasets
library(dplyr)
as.data.frame(names(hor.lab))
# validation dataset
val <- hor.lab[,c(2,4,7,12,32,24,31,33)]
val$hor<- as.factor(val$hor)
#val<-group_by(val,sitio,hor, add=T)
val.A<-group_by(val[val$hor=="A",], sitio)
val.A<-dplyr::summarise(val.A,H=first(horizonte),THICK=mean(thick),SAT=mean(Sat),
                 CO=mean(C_Ox),TB=mean(Total_bases),ESP.A=mean(ESP))
val.B<-group_by(val[val$hor=="B",], sitio)
val.B<-dplyr::summarise(val.B,ESP.B=mean(ESP))
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
names(c)<- c("values","Dataset","Property")
c<-c[complete.cases(c),]
c<-group_by(c,Property)
c<-c[!(c$values<60 & c$Property== "Sat.A (%)"),]

ggplot() +  facet_wrap(~Property,scales = 'free_y') + 
  #geom_jitter(data = c, mapping = aes(x=Dataset, y=values, color=Dataset), alpha = I(1/4))+
  geom_boxplot(data=c, mapping=aes(x=Dataset, y=values, color=Dataset), alpha = I(1/2))

