rm(list = ls())
name <- function(x) { as.data.frame(names(x))} 
#setwd("/home/mangelini/big/SEM_2nd_paper/")
setwd("/home/marcos/Documents/SEM2DSM1/Paper_2/data")

# load data ####
hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("64232_a_65039_FINAL.csv")
lab <- lab[,c(1,8,10,19)]

# profiles with odd horizons
# run after get object C line ~110
# hor[which(hor$sitio %in% 
#             unique(C[which(C$sitio %in% 
#                              names(table(C[,1])[table(C[,1])>1])),1])),]
hor$horizonte <- as.character(hor$horizonte)
hor[59,3] <- "X"
hor[332,3] <- "A"
hor[341,3] <- "X"
# hor[365,3] <- "X"
hor[390,3] <- "X"
hor <- hor[-364,] 
# mistakes in hor
#  
# en hor 64320 (64534 64793 64698) (64837 64903 65036)
# numero no usado 64437 64846 
# usado con otro proposito 64570 64757 

# C oxidable to OC ####
names(lab)[2] <- "OC"
lab$OC <- lab$OC * 1.3

# normalization of horizon names ####
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

# replace horizons with duplo analysis ####
# (2 measurements per sample) by mean(analysis[i]) 
repl <- hor[!is.na(hor$num_lab_r),]
replic <- unique(repl[,c(5,6)])
replic <- merge(x = replic,y = lab, by.x = "num_lab", by.y = "labid", all.x = T)
replic <- merge(x = replic,y = lab, by.x = "num_lab_r", by.y = "labid", all.x = T)

# get average
replic$OC <- (replic$OC.x + replic$OC.y) /2
replic$CEC <- (replic$CEC.x + replic$CEC.y) /2
replic$clay <- (replic$clay.x + replic$clay.y) /2
replic <- replic[,c(1,2,9,10,11)]

# delete rows duplicated from lab
lab <- lab[which(!lab$labid %in% replic$num_lab),]
lab <- lab[which(!lab$labid %in% replic$num_lab_r),]

# add replic to lab
replic <- replic[,-1]
names(replic)[1] <- "labid"
lab <- rbind(lab,replic)
# remove num_lab_r from hor
hor <- hor[,-6]

# add lab to hor
hor.lab <- merge(hor, lab, by.x="num_lab", by.y="labid", all=F)
# add coordenates to hor.lab
hor.xy <- merge(hor.lab,site[,c(2,3,5,6)], all.x = T)
write.csv(file = "validation_marcos_phd.csv",x =  hor.xy)
hor.xy <- hor.xy[(hor.xy$hor == "A" | hor.xy$hor == "B" |
                   hor.xy$hor == "C") & !is.na(hor.xy$hor),]
# thickness of standarized horizons
library(plyr)
hor.xy$sitio.hor <- paste(hor.xy$sitio,hor.xy$hor,sep = ".")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, mintop = min(prof_s))[,c(1,2)], by = "sitio.hor")
hor.xy <- merge(hor.xy, ddply(hor.xy,.(sitio.hor), summarise, maxbot = max(prof_i))[,c(1,2)], by = "sitio.hor")
# number of sites with A horizons
length(unique(hor.xy$sitio.hor[hor.xy$hor == "A"]))
# Which are the top horizons
#table(hor.xy$hor[hor.xy$mintop==0])
#Thickness
hor.xy$thick <- hor.xy$maxbot - hor.xy$mintop
samples <- unique(hor.xy[,c(2,11,12,7,13,14,15,8:10)])
samples[samples$sitio=="udic1-65",]

A <- samples[samples$hor=="A",]
B <- samples[samples$hor=="B",]
C <- samples[samples$hor=="C",]

names(A)[8:10] <- paste0(names(A)[8:10],".A")
names(B)[8:10] <- paste0(names(B)[8:10],".B")
names(C)[8:10] <- paste0(names(C)[8:10],".C")
samples <- merge(merge(A,B[,c(1,8:10)]),C[,c(1,8:10)])

# extract covariates values ####
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/")
# libraries
library(raster)
library(maptools)
library(sp)
library(rgdal)

# sdat files (dem covariates) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "river", "wdist","maxc","mrvbf","slope","twi","vdchn","water") 

# tif files (modis)
files_m <- list.files(pattern=".tif$")
# set names of covariates
header_m <- c("lstm", "lstsd", "evim", "evisd", "ndwi.a", "ndwi.b", "ndwi.bsd")
# samples to spatial object
names(samples)[c(2,3)] <- c("Y", "X")
coordinates(obj = samples) <- ~X+Y

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(samples) <- wgs84

# reproject
samples <- spTransform( samples, posgar98)

# extract values from files (.sdat)
stack <- list()
for(i in seq_along(files)) {
  samples@data[,length(samples@data)+1] <- NULL
  stack[[i]] <- readGDAL(files[i])
  proj4string(stack[[i]]) <- posgar98
  samples@data[,length(samples@data)+1] <- over(samples, stack[[i]])[,1]
  stack <- list()
  names(samples@data)[length(samples@data)] <- header[i]
}  

## extract values from modis files 
stack <- list()
# reproject endo to modis projection
samples <- spTransform( samples, modis)
for(i in seq_along(files_m)) {
  samples@data[,length(samples@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  samples@data[,length(samples@data)+1] <- over(samples, stack[[i]])[,1]
  stack <- list()
  names(samples@data)[length(samples@data)] <- header_m[i]
}  
samples <- spTransform( samples, posgar98)

samples <- as.data.frame(samples)

samples <- samples[,c(-2,-3,-4,-5)]
samples <- samples[,c(1,3,6,9,2,5,8,4,7,10,11:28)]
name(samples)
setwd("~/Documents/SEM2DSM1/Paper_2/data/")
write.csv(samples,"val.data.csv")










