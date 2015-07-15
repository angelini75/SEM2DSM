
rm(list=ls())
setwd("/media/L0135974_DATA/UserData/BaseARG/1_Sampling/Data")

hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("lab.csv")
repl <- hor[!is.na(hor$num_lab_r),]
repl.lab <- merge(repl,lab, by.x = "num_lab_r", by.y = "lab")
hor.lab.b <- merge(hor,lab, by.x = "num_lab", by.y = "lab",all.x = T)

hor.lab.c <- hor.lab.b[!is.na(hor.lab.b$num_lab_r),]
hor.lab.c <- hor.lab.c[,c(1,3,27:34)]
repl.lab <- repl.lab[,c(24,3,1,28:34)]
hor.lab.c<- unique(hor.lab.c)
repl.lab<- unique(repl.lab)

hor.lab.c$lab.p <- paste(hor.lab.c$num_lab,hor.lab.c$sitio, sep="&")
repl.lab$lab.p <- paste(repl.lab$num_lab,repl.lab$sitio, sep="&")

acc <-merge(hor.lab.c,repl.lab, by = "lab.p",all = T)
acc <- acc[complete.cases(acc),]
cor("CO_ox.x", "CO_ox.y", acc)

t.test(acc$CO_ox.y,acc$CO_ox.x,paired=TRUE)

summary(acc$CO_ox.y-acc$CO_ox.x)



hist(hor.lab.b$CO_ox)
hist(acc$CO_ox.y-acc$CO_ox.x, add =T)
summary(hor.lab.b$CO_ox)

#-------------------------------------#
# normalization of horizon names
hor.lab <- hor.lab.b[,c(1,3:6,28:34)]
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

val <- merge(site[,c(2,5,6)], hor.lab[hor.lab$hor=="A",c(2,13,4,5)], by = "sitio", all = T)
val <- val[complete.cases(val),]

library(plyr) 
val <- merge(val, ddply(val,.(sitio), summarise, mintop=min(prof_s))[,1:2], by= "sitio")
val <- merge(val, ddply(val,.(sitio), summarise, maxbot=max(prof_i))[,1:2], by= "sitio")
val$thick.A <- val$maxbot - val$mintop
hist(val$thick.A, add=T, col = "red")
summary(val$thick.A)

library(sp)
library(rgdal)
#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

thickA <- readGDAL("/media/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/rusults_thick.Ar.tif")
val.df<-val
coordinates(val)<- ~longitud+latitud

proj4string(val) <- wgs84
val <- spTransform(val, posgar98)
proj4string(thickA) <- posgar98
val@data[,length(val@data)+1] <- over(x=val,y=thickA)

SE.thickA <- sd((val@data$thick.A-(val@data$band1-2.42))[complete.cases(val@data$thick.A-(val@data$band1-2.42))])
summary((val@data$thick.A-val@data$band1)[complete.cases(val@data$thick.A-val@data$band1)])
hist(((val@data$thick.A-val@data$band1)^2)^(1/2), add=T, col = "green")
