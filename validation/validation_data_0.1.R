
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



#-------------------------------------#
#        Thickness A horizon          #
#-------------------------------------#
thick.A.val <- unique(hor.lab[hor.lab$mintop==0,c(7:9,12:15)])
thick.A.val <- thick.A.val[complete.cases(thick.A.val),]

library(sp)
library(rgdal)

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
#modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# load predicted
thick.A.pred <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/Thickness.sdat")

# thickness to spatial data frame
coordinates(thick.A.val)<- ~longitud+latitud

# extract values from  predicted
proj4string(thick.A.val) <- wgs84
thick.A.val <- spTransform(thick.A.val, posgar98)
proj4string(thick.A.pred) <- posgar98
thick.A.val@data[,length(thick.A.val@data)+1] <- over(x=thick.A.val,y=thick.A.pred)
names(thick.A.val@data)<-c("hor","measured","strata","area","percentage","predicted")


sd(as.vector(thick.A.val@data$measured-thick.A.val@data$predicted), na.rm=T)/sqrt(length(thick.A.val@data$measured))
par(pty="s")
plot(thick.A.val@data$predicted~thick.A.val@data$measured,col = "dark red",xlim=c(10, 40), ylim=c(10, 40))
abline(0,1)

# Statistical Inference Stratified Simple Random Sampling
# "Sampling for natural resource monitoring" de Gruiter, Brus, Knotters. pp.92
# ME + confident interval 

thick.A.val <- as.data.frame(thick.A.val)
thick.A.val$residuals <-thick.A.val$measured-thick.A.val$predicted
thick.A.val$residuals.sq <-(thick.A.val$measured-thick.A.val$predicted)^2

# Count n (samples per starta) and h (number of strata)
nh <- cbind(strata=as.vector(as.data.frame(table(thick.A.val$strata))[,1]),
            n=as.data.frame(table(thick.A.val$strata))[,2],
            h=1:12)
thick.A.val<- merge(thick.A.val, nh, by="strata")


library(dplyr)
#vignette("introduction", package = "dplyr")

X  <- group_by(thick.A.val, strata)
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
ll<-zSt-qt(0.975,104)*sqrt(VzSt)
# upper limit
ul<-zSt+qt(0.975,104)* sqrt(VzSt)
paste("ME (95%)=",round(ll,3),"<",round(zSt,3),"<",round(ul,3))

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
# ll<-zSt.s-qt(0.975,104)*sqrt(VzSt)
# # upper limit
# ul<-zSt+qt(0.975,104)* sqrt(VzSt)
# paste("ME (95%)=",round(ll,3),"<",round(zSt,3),"<",round(ul,3))





# 
# 
# 
# summary(lm(thick.A.val@data$predicted~thick.A.val@data$measured))
# var(val.A@data$C_Ox*1.30-val.A@data$band1, na.rm=T)
# 0.22^2
# 
# qt(0.975,93)
# 
# qt(0.975,1000000)
# qt(0.975,10)
# #-------------------------------------------------------------
# 
# library(sp)
# 
# strata<-readShapeSpatial("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Stratas/stratas_soil_distance_v2.shp")
# proj4string(strata) <- posgar98
# plot(strata)
# #library(Rsenal)
# #mapView(strata, burst = F)
# SE.OC <- sd((val.A@data$OC-(val.A@data$band1-2.42))[complete.cases(val.A@data$OC-(val.A@data$band1-2.42))])
# summary((val.A@data$OC-val.A@data$band1)[complete.cases(val.A@data$OC-val.A@data$band1)])
# hist(((val@data$thick.A-val@data$band1)^2)^(1/2), add=T, col = "green")
# 
# 
# 
# 
# ## recortes
# 
# # repl.lab <- merge(repl,lab, by.x = "num_lab_r", by.y = "labid")
# # hor.lab.b <- merge(hor,lab, by.x = "num_lab", by.y = "labid",all.x = T)
# # 
# # hor.lab.c <- hor.lab.b[!is.na(hor.lab.b$num_lab_r),]
# # hor.lab.c <- hor.lab.c[,c(1,3,27:43)]
# # repl.lab <- repl.lab[,c(24,3,1,28:43)]
# # hor.lab.c<- unique(hor.lab.c)
# # repl.lab<- unique(repl.lab)
# # 
# # hor.lab.c$lab.p <- paste(hor.lab.c$num_lab,hor.lab.c$sitio, sep="&")
# # repl.lab$lab.p <- paste(repl.lab$num_lab,repl.lab$sitio, sep="&")
# 
# # acc <-merge(hor.lab.c,repl.lab, by = "lab.p",all = T)
# # acc <- acc[complete.cases(acc),]
# # cor("CO_ox.x", "CO_ox.y", acc)
# # 
# # t.test(acc$CO_ox.y,acc$CO_ox.x,paired=TRUE)
# # 
# # summary(acc$CO_ox.y-acc$CO_ox.x)
# # 
# # 
# # 
# # hist(hor.lab.b$CO_ox)
# # hist(acc$CO_ox.y-acc$CO_ox.x, add =T)
# # summary(hor.lab.b$CO_ox)
# 
# 
# #extract values for validation thick.A and OC
# # hor.lab$id.h <- paste(hor.lab$sitio,hor.lab$hor, sep=".")
# # 
# # val.A <- merge(site[,c(2,5,6)], hor.lab[hor.lab$hor=="A", c(1,2,22,4,5,12,19,20,21,23)], 
# #                by = "sitio", all = T)[c(-118:-128),]
# # val.B <- merge(site[,c(2,5,6)], hor.lab[hor.lab$hor=="B", c(1,2,22,4,5,12,19,20,21,23)], 
# #                by = "sitio", all = T)[c(-214:-224),]
# #val <- val[complete.cases(val),]
