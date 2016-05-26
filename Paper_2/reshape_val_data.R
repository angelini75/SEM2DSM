rm(list = ls())
name <- function(x) { as.data.frame(names(x))} 
#setwd("/home/mangelini/big/SEM_2nd_paper/")
setwd("/home/marcos/Documents/SEM2DSM1/Paper_2/data")

# load data ####
hor <- read.csv("Ficha_campo_hor.csv")
site <- read.csv("Ficha_campo_sitio.csv")
lab <- read.csv("64232_a_65039_FINAL.csv")
lab <- lab[,c(1,8,10,19)]

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
hor.xy <- merge(hor.lab,site[,c(2,5,6)], all.x = T)

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
hor.xy <- unique(hor.xy[,c(2,11,12,7,13,14,15,8:10)])

samples <- hor.xy[which(!hor.xy$sitio %in%
hor.xy$sitio[hor.xy$hor== "AB|BA" | hor.xy$hor== "BC" | hor.xy$hor== "AC" ]),]

A <- samples[samples$hor=="A",]
B <- samples[samples$hor=="B",]
C <- samples[samples$hor=="C",]
names(A)[8:10] <- paste0(names(A)[8:10],".A")
names(B)[8:10] <- paste0(names(B)[8:10],".B")
names(C)[8:10] <- paste0(names(C)[8:10],".C")
samples <- merge(merge(A,B[,c(1,8:10)]),C[,c(1,8:10)])

write.csv(samples, "validation.csv")

# it needs covariates


setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/")
# X <- 
# Y <- 

# sdat files (dem) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "river", "wdist","maxc","mrvbf","slope","twi","vdchn","water") 

#files_posgar <- files[c(1:4,9,10,12,13)]
# header <- header[c(1:4,9,10,12,13)]
# files250p <- files[11]
# header250p <- header[11]

# tif files (modis)
files_m <- list.files(pattern=".tif$")
# files250m <- files[c(7,8)]
# header250m <- header[c(7,8)]
# files1km <- files[c(5,6)]
header_m <- c("lstm", "lstsd", "evim", "evisd", "ndwi.a", "ndwi.b", "ndwi.bsd")

# files_n <- list.files(path = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/output/",pattern=".tif$")
# header_n <- gsub(pattern = ".tif",replacement = "",x = files_n)
# header_n <- paste("X",header_n, sep = "")
# files_n <- paste("output/", files_n, sep = "")
# 
# files_M <- list.files(path = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/MCD43A4/",pattern=".tif$")
# header_M <- gsub(pattern = ".tif",replacement = "",x = files_M)
# files_M <- paste("MCD43A4/", files_M, sep = "")

coordinates(obj = endo) <- ~X+Y


#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(endo) <- wgs84


# However, if you got the data from a RasterLayer (it looks like it)
# you can avoid the above and simply do:
# library(raster)
# pts <- rasterToPoints(r, spatial=TRUE)

# use spTransform
endo <- spTransform( endo, posgar98)

# extract values from files (.sdat)
stack <- list()
for(i in 1:length(files)) {
  endo@data[,length(endo@data)+1] <- NULL
  stack[[i]] <- readGDAL(files[i])
  proj4string(stack[[i]]) <- posgar98
  endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
  stack <- list()
  names(endo@data)[length(endo@data)] <- header[i]
}  

## extract values from modis files 
stack <- list()
# reproject endo to modis projection
endo <- spTransform( endo, modis)
for(i in 1:length(files_m)) {
  endo@data[,length(endo@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
  stack <- list()
  names(endo@data)[length(endo@data)] <- header_m[i]
}  
# files within folder output
# 
# for(i in 1:length(files_n)) {
#   endo@data[,length(endo@data)+1] <- NULL
#   stack[[i]] <- readGDAL(files_n[i])
#   proj4string(stack[[i]]) <- modis # change projection
#   endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
#   stack <- list()
#   names(endo@data)[length(endo@data)] <- header_n[i]
# }
# 
# # files MCD43A4
# for(i in 1:length(files_M)) {
#   endo@data[,length(endo@data)+1] <- NULL
#   stack[[i]] <- readGDAL(files_M[i])
#   proj4string(stack[[i]]) <- modis # change projection
#   endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
#   stack <- list()
#   names(endo@data)[length(endo@data)] <- header_M[i]
# }
#image(raster(("mod13q1_tot_mean.tif")))















