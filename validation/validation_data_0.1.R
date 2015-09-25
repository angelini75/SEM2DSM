
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
hor.lab <- hor.lab.b[,c(1,3:6,28:43)]
hor.lab$hor <- NA
as.data.frame(t(table(hor.lab$horizonte)))
hor.lab$hor[grep("^A$",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^a$",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("A1",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("A2",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^2A",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^An",hor.lab$horizonte)] <- "A"

hor.lab$hor[grep("E",hor.lab$horizonte)] <- "E"

hor.lab$hor[grep("AC",hor.lab$horizonte)] <- "AC"

hor.lab$hor[grep("AB",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("BA",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("Ba",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("B/A",hor.lab$horizonte)] <- "AB|BA"

hor.lab$hor[grep("^Bt",hor.lab$horizonte)] <- "B"
hor.lab$hor[grep("^2Bt",hor.lab$horizonte)] <- "B"

hor.lab$hor[grep("^BC",hor.lab$horizonte)] <- "BC"
hor.lab$hor[grep("Bc",hor.lab$horizonte)] <- "BC"
hor.lab$hor[grep("^2BC",hor.lab$horizonte)] <- "BC"

hor.lab$hor[grep("^C",hor.lab$horizonte)] <- "C"
hor.lab$hor[grep("^c",hor.lab$horizonte)] <- "C"
hor.lab$hor[grep("^2C",hor.lab$horizonte)] <- "C"

as.data.frame(t(table(hor.lab$horizonte[is.na(hor.lab$hor)])))

#-------------------------------------#




repl.lab <- merge(repl,lab, by.x = "num_lab_r", by.y = "labid")
hor.lab.b <- merge(hor,lab, by.x = "num_lab", by.y = "labid",all.x = T)

hor.lab.c <- hor.lab.b[!is.na(hor.lab.b$num_lab_r),]
hor.lab.c <- hor.lab.c[,c(1,3,27:43)]
repl.lab <- repl.lab[,c(24,3,1,28:43)]
hor.lab.c<- unique(hor.lab.c)
repl.lab<- unique(repl.lab)

hor.lab.c$lab.p <- paste(hor.lab.c$num_lab,hor.lab.c$sitio, sep="&")
repl.lab$lab.p <- paste(repl.lab$num_lab,repl.lab$sitio, sep="&")

# acc <-merge(hor.lab.c,repl.lab, by = "lab.p",all = T)
# acc <- acc[complete.cases(acc),]
# cor("CO_ox.x", "CO_ox.y", acc)
# 
# t.test(acc$CO_ox.y,acc$CO_ox.x,paired=TRUE)
# 
# summary(acc$CO_ox.y-acc$CO_ox.x)
# 
# 
# 
# hist(hor.lab.b$CO_ox)
# hist(acc$CO_ox.y-acc$CO_ox.x, add =T)
# summary(hor.lab.b$CO_ox)

#-------------------------------------#
# normalization of horizon names
hor.lab <- hor.lab.b[,c(1,3:6,28:43)]
hor.lab$hor <- NA
as.data.frame(t(table(hor.lab$horizonte)))
hor.lab$hor[grep("^A$",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^a$",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("A1",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("A2",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^2A",hor.lab$horizonte)] <- "A"
hor.lab$hor[grep("^An",hor.lab$horizonte)] <- "A"

hor.lab$hor[grep("E",hor.lab$horizonte)] <- "E"

hor.lab$hor[grep("AC",hor.lab$horizonte)] <- "AC"

hor.lab$hor[grep("AB",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("BA",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("Ba",hor.lab$horizonte)] <- "AB|BA"
hor.lab$hor[grep("B/A",hor.lab$horizonte)] <- "AB|BA"

hor.lab$hor[grep("^Bt",hor.lab$horizonte)] <- "B"
hor.lab$hor[grep("^2Bt",hor.lab$horizonte)] <- "B"

hor.lab$hor[grep("^BC",hor.lab$horizonte)] <- "BC"
hor.lab$hor[grep("Bc",hor.lab$horizonte)] <- "BC"
hor.lab$hor[grep("^2BC",hor.lab$horizonte)] <- "BC"

hor.lab$hor[grep("^C",hor.lab$horizonte)] <- "C"
hor.lab$hor[grep("^c",hor.lab$horizonte)] <- "C"
hor.lab$hor[grep("^2C",hor.lab$horizonte)] <- "C"

as.data.frame(t(table(hor.lab$horizonte[is.na(hor.lab$hor)])))

#-------------------------------------#
#extract values for validation thick.A and OC
hor.lab$id.h <- paste(hor.lab$sitio,hor.lab$hor, sep=".")

val.A <- merge(site[,c(2,5,6)], hor.lab[hor.lab$hor=="A", c(1,2,22,4,5,12,19,20,21,23)], 
               by = "sitio", all = T)[c(-118:-128),]
val.B <- merge(site[,c(2,5,6)], hor.lab[hor.lab$hor=="B", c(1,2,22,4,5,12,19,20,21,23)], 
               by = "sitio", all = T)[c(-214:-224),]
#val <- val[complete.cases(val),]

library(plyr) 
# create dataset for A horizon
val.A <- merge(val.A, ddply(val.A,.(sitio), summarise, mintop=min(prof_s))[,1:2], by= "sitio")
val.A <- merge(val.A, ddply(val.A,.(sitio), summarise, maxbot=max(prof_i))[,1:2], by= "sitio")
val.A$thick.A <- val.A$maxbot - val.A$mintop
summary(val.A$thick.A)
val.A<- val.A[!is.na(val.A$num_lab),]
as.data.frame(val.A$num_lab)
val.A<- val.A[-86,]
as.data.frame(names(val.A))
val.A<- val.A[,c(-6,-7)]
val.A<- unique(val.A)

library(sp)
library(rgdal)
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

OC <- readGDAL("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/OC.sdat")
val.df<-val.A
coordinates(val.A)<- ~longitud+latitud

proj4string(val.A) <- wgs84
val.A <- spTransform(val.A, posgar98)
proj4string(OC) <- posgar98
val.A@data[,length(val.A@data)+1] <- over(x=val.A,y=OC)

sd(as.vector(val.A@data$C_Ox*1.3-val.A@data$band1), na.rm=T)/sqrt(length(val.A@data$sitio))
plot(val.A@data$C_Ox*1.30~val.A@data$band1,col = "dark red")
abline(0,1)

summary(lm(val.A@data$C_Ox*1.3~val.A@data$band1))
var(val.A@data$C_Ox*1.30-val.A@data$band1, na.rm=T)
0.22^2
qt(0.975,100)
qt(0.975,1000000)
qt(0.975,10)
#-------------------------------------------------------------

library(sp)
library(maptools)
strata<-readShapeSpatial("/media/marcos/L0135974_DATA/UserData/BaseARG/1_Sampling/Stratas/stratas_soil_distance_v2.shp")
proj4string(strata) <- posgar98
plot(strata)
#library(Rsenal)
#mapView(strata, burst = F)
SE.OC <- sd((val.A@data$OC-(val.A@data$band1-2.42))[complete.cases(val.A@data$OC-(val.A@data$band1-2.42))])
summary((val.A@data$OC-val.A@data$band1)[complete.cases(val.A@data$OC-val.A@data$band1)])
hist(((val@data$thick.A-val@data$band1)^2)^(1/2), add=T, col = "green")
