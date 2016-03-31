
rm(list = ls())
name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data")

# load data
hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("64232_a_65039_FINAL.csv")
lab$X0CaCO3 <- as.character(lab$X0CaCO3)
# lab$X0CaCO3[lab$X0CaCO3 == ""] <- 0
# lab$X0CaCO3[lab$X0CaCO3 == "*"] <- 1
# lab$X0CaCO3[lab$X0CaCO3 == "**"] <- 2
# lab$X0CaCO3[lab$X0CaCO3 == "***"] <- 3
lab$X0CaCO3 <- as.numeric(lab$X0CaCO3)
#hist(lab$pHKCl)

# clean and transform C oxidable to OC
lab$oc <- lab$C_Ox * 1.3
lab$ESP <- log10(lab$ESP)
# # Horizons A1 and A2, same site, two analysis
# #hor$num_lab[hor$num_lab==64561 & !is.na(hor$num_lab) ][2]<-64562
# #lab$labid[lab$labid==64561& !is.na(lab$labid)][2] <-64562
# hor$num_lab_r[hor$num_lab_r==64599& !is.na(hor$num_lab_r) ]<- 64559
#lab[lab$labid == 65012,] <- lab[lab$labid == 64879,]
repl <- hor[!is.na(hor$num_lab_r),]
#--- MEASUREMENT ERROR AS STANDARD DEVIATION OF THE ERROR
replic <- unique(repl[,c(23,27)])
replic <- merge(x = replic,y = lab, by.x = "num_lab", by.y = "labid", all.x = T)
replic <- merge(x = replic,y = lab, by.x = "num_lab_r", by.y = "labid", all.x = T)
replic$pHKCl.x <- as.numeric(replic$pHKCl.x)
replic$pHKCl.y <- as.numeric(replic$pHKCl.y)
# replic <- replic[replic$num_lab != 64553,]
# replic <- replic[replic$num_lab != 65012,]



par(pty = "s")
par(mfrow = c(3, 3))
for (i in 3:28) {
  lim <- c(min(rbind(replic[,i],replic[,i+26]), na.rm = T), 
           max(rbind(replic[,i],replic[,i+26]), na.rm = T))
  plot(replic[,i] ~ replic[,i + 26], main = names(replic)[i], 
       xlab = "measurement 1", ylab = "measurement 2", 
       col = "dark red", xlim = lim, ylim = lim)
  abline(0,1)
  abline(lm(replic[,i] ~ replic[,i + 26]), col = "blue")
}

error <- replic[,1:28]
for (i in 3:28) {
  error[,i] <- replic[,i] - replic[,i + 26]
  print(paste(names(error)[i],mean(error[,i], na.rm = T) * 
                sqrt(length(error[complete.cases(error[,i]),i])) /
                sd(error[,i], na.rm = T)))
}
summary(error)
measurement.error <- data.frame(Property = NA, MEANe = NA, VARe = NA, SDe = NA, S.Var = NA)
for (i in 3:28) {
  measurement.error[i - 2,] <- c(names(error)[i], mean(error[,i], na.rm = T),
                                 0.5 * var(error[,i], na.rm = T), sqrt(0.5 * var(error[,i], na.rm = T)),
                                 var(replic[,i], na.rm = T))
}

#####
e <- matrix(NA,nrow = 36,ncol = 27)
e <- as.data.frame(e)
e[,1] <- replic[,2]
names(e)[1] <- names(replic)[2]
for (i in 3:28) {
  e[,i - 1] <- replic[,i] - replic[,i + 26]
  names(e)[i - 1] <- names(replic)[i]
}
summary(e)
#hist(e$Total_bases.x)
E <- e[e$pHw.x > 1.5 | e$pHw.x < -1.5,]
E <- E[!is.na(E$Sat_w.x),]
summary(E)

m <- matrix(NA,nrow = 36,ncol = 26)
m <- as.data.frame(m)
m[,1] <- replic[,2]
names(m)[1] <- names(replic)[2]

for (i in 3:28) {
  m[,i - 1] <- (replic[,i] + replic[,i + 26]) / 2
  names(m)[i - 1] <- names(replic)[i]
}


#-------------------------------------#
# normalization of horizon names
as.data.frame(names(hor))
hor <- hor[,c(2:5,23,27)]
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

##### clay ratio measuremt error
# take clay B from each place and divide by clay A of same sites
# repl <- hor[!is.na(hor$num_lab_r),]
# replic <- unique(repl[,c(1,5,6,7)])
# replic <- merge(x = replic,y = lab, by.x = "num_lab", by.y = "labid", all.x = T)
# replic <- merge(x = replic,y = lab, by.x = "num_lab_r", by.y = "labid", all.x = T)
# 
# repl.clay <- replic[replic$hor == "A"|replic$hor == "B", c(3,4,22,48)]
# repl.clay <- repl.clay[-1,]
# repl.clay$sitio <- as.character(repl.clay$sitio)
# repl.clay$sitio <- as.factor(repl.clay$sitio)
# 
# 
# write.csv(repl.clay, "/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Data/repl.clay.csv")
    # lazzy: LibreOffice::Calc compute manually four possible ratios:
    # (1)sample1:B/sample1:A; (2)sample1:B/sample2:A; (3)sample2:B/sample1:A; (4)sample2:B/sample2:A, then, six possible residuals
    # absolute values of: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4

lazzy.residuals <- read.csv("repl.clay.residuals.csv")[,1]
lazzy.ratios <- read.csv("repl.clay.ratios.csv")[,1]
    #measurement.error <- data.frame(Property = NA, MEANe = NA, VARe = NA, SDe = NA, S.Var = NA)

measurement.error[27,1:5] <- c("clay.ratio", mean(lazzy.residuals),
                                 0.5 * var(lazzy.residuals), sqrt(0.5 * var(lazzy.residuals)),
                                 var(lazzy.ratios))
    # comapre VARe/Var(property)
measurement.error[,2] <- as.numeric(measurement.error[,2])
measurement.error[,3] <- as.numeric(measurement.error[,3])
measurement.error[,4] <- as.numeric(measurement.error[,4])
measurement.error[,5] <- as.numeric(measurement.error[,5])
# comapre VARe/Var(property)
measurement.error$ratio <- measurement.error$VARe / measurement.error$S.Var
write.csv(measurement.error,"/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/measurement.error.2.csv")





#---
## 
##### replace horizons with duplo analysis (2 measurements per sample) by mean(analysis[i])
s <- m$num_lab
for (i in 1:length(s)) {
  lab[lab$labid == s[i],] <- as.vector(as.matrix(m[m$num_lab == s[i],]))
}


# add coordenates to hor
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
name(hor.lab)
# Thick.A
# OC.A
# TB.A
# Sat.A
# ESP.A
# ESP.B

d.A <- unique(hor.lab[hor.lab$hor == "A",c(2,13,44,32,33,34,36)])
d.A <- d.A[-93,] # without lab data
names(d.A)[2:7] <- paste0(names(d.A)[2:7],".A")
d.B <- unique(hor.lab[hor.lab$hor == "B",c(2,34,36)])
names(d.B)[2:3] <- paste0(names(d.B)[2:3],".B")  
d.clay <- merge(d.A[,c(1,7)],d.B[,c(1,3)], all.x=T)
d.clay$bt <- d.clay$clay.B/d.clay$clay.A
d.A <- merge(d.A, d.clay[,c(1,4)])
ratio.mean <- read.csv("repl.clay.mean_ratios.csv")
ratio.mean$sitio <- as.character(ratio.mean$sitio)
d.A$sitio <- as.character(d.A$sitio)

for(i in 1:11){
d.A[d.A$sitio == ratio.mean$sitio[i],8] <- ratio.mean$mean_ratio[i]
}
#d.A <- d.A[-24,] # co1-19 is not reliable
d.A <- d.A[,c(-1,-7)]
d.B <- d.B[-103,]
d.B <- data.frame(ESP.B=d.B[,c(-1,-3)])

d.stat <- matrix(data = NA,nrow = 6,ncol = 7,
                 dimnames = list(c("Min", "Median", "Mean", "Max", "SD", "SS"),
                                 c(names(d.A),names(d.B))))
d.stat[1,1:6] <- apply(X = d.A,FUN = min,2, na.rm = T)
d.stat[1,7] <- apply(X = d.B,FUN = min,2, na.rm = T)
d.stat[2,1:6] <- apply(X = d.A,FUN = median,2, na.rm = T)
d.stat[2,7] <- apply(X = d.B,FUN = median,2, na.rm = T)
d.stat[3,1:6] <- apply(X = d.A,FUN = mean,2, na.rm = T)
d.stat[3,7] <- apply(X = d.B,FUN = mean,2, na.rm = T)
d.stat[4,1:6] <- apply(X = d.A,FUN = max,2, na.rm = T)
d.stat[4,7] <- apply(X = d.B,FUN = max,2, na.rm = T)
d.stat[5,1:6] <- apply(X = d.A,FUN = sd,2, na.rm = T)
d.stat[5,7] <- apply(X = d.B,FUN = sd,2, na.rm = T)

sum((mean(d.A$Total_bases.A, na.rm=T) - d.A$Total_bases.A)^2, na.rm = T)
sum((mean(d.A$bt, na.rm=T) - d.A$bt)^2, na.rm = T)

for (i in 1:6) {
  d.stat[6,i] <- sum((mean(d.A[complete.cases(d.A[,i]),i]) -
                        d.A[complete.cases(d.A[,i]),i]) ^ 2 )
}
d.stat[6,7] <- sum((mean(d.B[complete.cases(d.B[,1]),1]) - d.B[complete.cases(d.B[,1]),1]) ^ 2)
dimnames(d.stat)[[2]] <- c("thick.A", "oc.A", "tb.A", "sat.A", "esp.A", "bt", "esp.B")
library(pastecs)
d.stat1 <- stat.desc(d.A)
d.stat1[,7] <- NA
d.stat1[,7] <- stat.desc(d.B)
names(d.stat1)[7] <- "esp.B"

Var <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/lm.variance.error.csv")[,-1]
Ns <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/N.march2016.csv")[,-1]
#################### VALIDATION ####################
library(sp)
library(rgdal)

#-------------------------------------#
#        Thickness A horizon          #
#-------------------------------------#
# SP is Soil Property
SP.val <- unique(hor.lab[hor.lab$hor == "A",c(8:10,13:16)])
SP.val <- SP.val[complete.cases(SP.val),]

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_thick.Ar.tif")
# thickness to spatial data frame
coordinates(SP.val) <- ~ longitud + latitud
# extract values from  predicted
proj4string(SP.val) <- wgs84
SP.val <- spTransform(SP.val, posgar98)
proj4string(SP.pred) <- posgar98
SP.val@data[,length(SP.val@data) + 1] <- over(x = SP.val, y = SP.pred)
names(SP.val@data) <- c("hor","measured","strata","area","percentage","predicted")

par(pty = "s")
par(mfrow = c(3, 3))
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Thickness A horizon (cm)", 
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = c(10, 40), ylim = c(10, 40))
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")

# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted 
write.csv(SP.val, "res.thick.csv")
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
N  <- as.numeric(length(X$residuals)) - 1 #number of samples
H  <- 12 #number of strata
Nh <- as.data.frame(dplyr::summarise(X,n()))[,2] # number of samples per strata
zh <- as.data.frame(dplyr::summarise(X, mean(residuals)))[,2] # samle mean of the stratum h
sum(Nh*zh)/N

######################## ME <- mean error of top horizon thickness of the area (considering area)
Ah <- as.data.frame(dplyr::summarise(X, mean(area)))[,2] # area of the stratum h
A <- sum(Ah) # total area
ah <- Ah/A # weights 
zSt <- sum(ah * zh) # mean error
paste("ME =", round(zSt,3))
# variance of zSt
VzSt <- sum((ah ^ 2) * dplyr::summarise(X, var(residuals))[2]) # error variance
VzSt
# 95% confidence
# lowwer limit
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="thick.A"] * Ns$sd[Ns$name=="thick.A"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="thick.A"] * Ns$sd[Ns$name=="thick.A"]^2))

# fill report table
report <- data.frame(Soil_property = NA, ME.ll=NA, ME = NA, ME.ul=NA, RMSE = NA, SS = NA, theta_mean = NA, theta_median = NA)
report[1,1:8] <- c("Thick.A",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)


#----------------------------------------#
#        Organic Carbon A horizon        #
#----------------------------------------#c(8:10,13:16)
as.data.frame(names(hor.lab))
rm(SP.val)
rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(8:10,44,14:16)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_oc.Ar.tif")
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
write.csv(SP.val, "res.oc.csv")
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="oc.A"] * Ns$sd[Ns$name=="oc.A"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="oc.A"] * Ns$sd[Ns$name=="oc.A"]^2))

# fill report table
report[2,1:8] <- c("oc.A",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)


#-------------------------------------#
#        Total Bases A horizon        #
#-------------------------------------#c(8:10,13:16)
as.data.frame(names(hor.lab))
rm(SP.val)
rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(8:10,32,14:16)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_tb.Ar.tif")
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
write.csv(SP.val, "res.tb.csv")
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="tb.A"] * Ns$sd[Ns$name=="tb.A"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="tb.A"] * Ns$sd[Ns$name=="tb.A"]^2))

# fill report table
report[3,1:8] <- c("tb.A",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)

#-----------------------------------------#
#        Base Saturation A horizon        #
#-----------------------------------------#c(8:10,13:16)
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A", c(8:10,33,14:16)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_sat.Ar.tif")
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
write.csv(SP.val, "res.sat.csv")
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="sat.A"] * Ns$sd[Ns$name=="sat.A"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="sat.A"] * Ns$sd[Ns$name=="sat.A"]^2))

# fill report table
report[4,1:8] <- c("sat.A",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)


#-----------------------------#
#        ESP A horizon        #
#-----------------------------#c(8:10,13:16)
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "A",c(8:10,34,14:16)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_esp.Ar.tif")
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
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "ESP A horizon (log(%))", xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
write.csv(SP.val, "res.espA.csv")
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="esp.A"] * Ns$sd[Ns$name=="esp.A"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="esp.A"] * Ns$sd[Ns$name=="esp.A"]^2))

# fill report table
report[5,1:8] <- c("esp.A",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)

#-----------------------------#
#        ESP B horizon        #
#-----------------------------#c(8:10,13:16)
as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(hor.lab[hor.lab$hor == "B",c(8:10,34,14:16)])
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_esp.Br.tif")
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
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "ESP B horizon (log(%))", xlab = "measured",
     ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
#write.csv(SP.val, "res.espB.csv")
SP.val$residuals.sq <- (SP.val$measured - SP.val$predicted)^2 
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="esp.B"] * Ns$sd[Ns$name=="esp.B"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="esp.B"] * Ns$sd[Ns$name=="esp.B"]^2))

# fill report table
report[6,1:8] <- c("esp.B",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)


#-----------------------------------------#
#        Clay ratio B/A horizons          #
#-----------------------------------------#c(8:10,13:16)

d.A <- unique(hor.lab[hor.lab$hor == "A",c(2,36)])
d.A <- d.A[-93,] # without lab data
names(d.A)[2] <- paste0(names(d.A)[2],".A")
d.B <- unique(hor.lab[hor.lab$hor == "B",c(2,36)])
names(d.B)[2] <- paste0(names(d.B)[2],".B")  
d.clay <- merge(d.A,d.B, all.x=T)
d.clay$bt <- d.clay$clay.B/d.clay$clay.A
d.A <- merge(d.A, d.clay[,c(1,4)])
ratio.mean <- read.csv("repl.clay.mean_ratios.csv")
ratio.mean$sitio <- as.character(ratio.mean$sitio)
d.A$sitio <- as.character(d.A$sitio)

for(i in 1:11){
  d.A[d.A$sitio == ratio.mean$sitio[i],3] <- ratio.mean$mean_ratio[i]
}

as.data.frame(names(hor.lab))
rm(SP.val);rm(SP.pred)
SP.val <- unique(merge(hor.lab[hor.lab$hor == "A", c(2,8:10,14:16)],d.A[,c(1,3)], all.y=T))
SP.val <- SP.val[,c(2:4,8,5:7)]
SP.val <- SP.val[complete.cases(SP.val),]

# load predicted
SP.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/march.lm_btr.tif")
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
plot(SP.val@data$predicted ~ SP.val@data$measured, main = "Clay ratio B/A horizons",
     xlab = "measured", ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
abline(0,1)
abline(lm(SP.val@data$predicted ~ SP.val@data$measured), col = "red")
# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

SP.val <- as.data.frame(SP.val)
SP.val$residuals <- SP.val$measured - SP.val$predicted
write.csv(SP.val, "res.bt.csv")
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
ME.ll <- zSt - qt(0.975, N - 1) * sqrt(VzSt)
# upper limit
ME.ul <- zSt + qt(0.975, N - 1) * sqrt(VzSt)

ME <- zSt

############################################ MSE (mean squared error)
zh.s <- as.data.frame(dplyr::summarise(X, mean(residuals.sq)))[,2] # sample mean of squared errors
sum(Nh * zh.s)/N
# ME <- mean error of top horizon thickness of the area (considering area)
zSt.s <- sum(ah * zh.s) # mean of squared error
paste("MSE =", round(zSt.s, 3))
paste("RMSE =", round(sqrt(zSt.s), 3))
MSE <- zSt.s
RMSE <- sqrt(zSt.s)

theta_mean <- mean(SP.val$residuals.sq / (Var$variance[Var$property=="bt"] * Ns$sd[Ns$name=="bt"]^2))
theta_median <- median(SP.val$residuals.sq / (Var$variance[Var$property=="bt"] * Ns$sd[Ns$name=="bt"]^2))

# fill report table
report[7,1:8] <- c("bt",ME.ll,ME,ME.ul,RMSE,SS, theta_mean, theta_median)

report
d.stat <- d.stat[,c(1:5,7,6)]


#########################################





setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

report$R2 <- NA
for (i in 1:7) {
  report$R2[i] <- 1 - (as.numeric(report$SS[i]) / d.stat[6,i])
}

d.stat <- rbind(d.stat,NA)
d.stat[7,6] <- length(d.A$bt)
d.stat[7,1] <- length(unique(hor.lab[hor.lab$hor == "A",c(8:10,13:16)])[,4])
d.stat[7,2] <- 92
d.stat[7,3] <- 92
d.stat[7,4] <- 92
d.stat[7,5] <- 92
d.stat[7,7] <- 100

write.csv(report, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/report.validation.lm.csv")
# write.csv(d.stat, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/summary.valdata.csv")
# write.csv(d.stat1, "/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/d.stat.val.csv")
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

