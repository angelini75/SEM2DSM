############### #### ### ## # - SHAPING AND MERGING DATA WITH COVARIATES - # ## ### #### ###############
#### #### #### #### #### #### FOR SEM CALIBRATION
# Purpose : soil data transformation to be merged with covariates (rasters)
# Maintainer : Marcos Angelini (marcos.angelini@wur.nl);
# Contributions : 
# Status : alpha
# Note : 
#     This section concerns about grouping same type of horizons by profile, making a 
#     weighted mean of soil properties and taking min and max boundaries of horizons; 
#   **Imput data : calib.data-0.1.2.csv; rasters from /media/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling**
#   **Output data : calib.data-2.0.csv.

# sessionInfo(@RStudio desktop) lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# sessionInfo(@RStudio Server) http://h2335862.stratoserver.net:8787/ (12 cores)
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-unknown-linux-gnu (64-bit)


rm(list=ls())
#install.packages('reshape')
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
# install.packages('soiltexture')
#install.packages('corrgram')
# install.packages("plyr")
#library(soiltexture)
library(reshape)
library(plyr)
library(sp)
library(lattice) # required for trellis.par.set():
#trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()
library(corrgram)
#library(gdalUtils)

name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))



#################### RE-SHAPING DATA  ########
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")

D <-read.csv("calib.data-0.1.1.csv")

#selected soil properties
name(D)
d <- cbind(D[,c(1:4,61,98,95,94,7:16,30:33,18,17,20:27,51,42,55,87,88,90,91,77,81)])
as.numeric(d$a_ph_kcl == 73.2) 
as.numeric(d$a_ph_kcl == 0.3) 
d[1482,13]<- 7.32
d[1516,13]<- NA
# names of soil properties
names(d)[6]<-"horizon"
names(d)[7]<-"top"
names(d)[8]<-"bottom"
names(d)[9:10]<-c("a_sum_bases","a_CEC")
names(d)[23]<-"a_OC"
#d[,23] <- d[,23]*0.581  #materia organica a CO
d$horizon <- as.character(d$horizon)

# standardising horizon names
t(table(d$horizon))
d$hor <-""
d<- cbind(d[,c(1:6,42,7:41)])
#if(d$horizon==)
count(d[grep("A1",d$horizon),c(6,7)])
d$hor[grep("A1",d$horizon)] <- "A"
d$hor[grep("^A$",d$horizon)] <- "A"
d$hor[grep("Ap",d$horizon,ignore.case = T)] <- "A"
d$hor[grep("^IA$",d$horizon)] <- "A"
d$hor[grep("Enlame",d$horizon)] <- "A"
d$hor[grep("^A2",d$horizon)] <- "E"
d$hor[grep("^IIA2",d$horizon)] <- "E"
d$hor[grep("^A3",d$horizon)] <- "AB|BA"
d$hor[grep("^IA3",d$horizon)] <- "AB|BA"
d$hor[grep("^B1",d$horizon)] <- "AB|BA"
d$hor[grep("^IB1",d$horizon)] <- "AB|BA"
d$hor[grep("AB",d$horizon)] <- "AB|BA"
d$hor[grep("^B1",d$horizon)] <- "AB|BA"
d$hor[grep("A y",d$horizon)] <- "AB|BA"
d$hor[grep("ByA",d$horizon)] <- "AB|BA"
d$hor[grep("A/B",d$horizon)] <- "AB|BA"
d$hor[grep("B2",d$horizon)] <- "B"
d$hor[grep("Bx",d$horizon)] <- "B"
d$hor[grep("Bca",d$horizon)] <- "B"

d$hor[grep("B3",d$horizon)] <- "BC"
d$hor[grep("B2",d$horizon)] <- "B"
d$hor[grep("C",d$horizon)] <- "C"
d$hor[grep("^c$",d$horizon)] <- "C"
d$hor[grep("AC",d$horizon)] <- "AC"
d$hor[grep("R",d$horizon)] <- "C"
ex <- c("^I$","^II$","^III$","^IV$","^V$","^VI$","^IVg$")
for(i in 1:length(ex)){
  d$hor[grep(ex[i],d$horizon)] <- "X"
}

# delete rows without horizon information
d0 <-d[!is.na(d$top),] 
# delete rows without analysis
# d0 <- d0[!is.na(d0$a_ph_h2o) & !is.na(d0$a_OC),]

# bottom == NA <- top + 20 cm & misstiping errors
d0[is.na(d0$bottom)& d0$top <100, 3:10]
d0$top[d0$id.h==586] <- 140
d0$bottom[d0$id.h==640] <- 19
d0$top[d0$id.h==2238] <- 198
d0$bottom[d0$id.h==2220] <- 30
d0 <-d0[-985,] # replicated horizon 
d0$bottom[is.na(d0$bottom)] <-d0$top[is.na(d0$bottom)] + 20 # NA values at bootom are consider 20 cm more than top
d0 <- d0[-2019,] # Enlame from -4 to 0
#### weighted mean per id.hor
# create id.hor to identify more than 1 B horizon in the same profile, for instance
d0$id.hor <- paste(d0$id.p, d0$hor,sep="_")
d0$thick <- d0$bottom - d0$top
#(d0$bottom - d0$top)<1
#returns sum of thick by id.hor. This requires installation of the library plyr
thickness <- ddply(d0,.(id.hor), summarise,s.thick=sum(thick))
#merge thickness and s.pr
d1 <- merge(d0, thickness, by = "id.hor", all=T)
#weights
# d1$weight <- d1$thick / d1$s.thick

# drenaje

dren <- unique(d1[,c(5,38:41)])

summary(dren)
levels(dren$p_permeabilidad)[1] <- NA
levels(dren$p_escurrimiento)[1] <- NA
levels(dren$p_anegamiento)[1] <- NA
levels(dren$p_drenaje)[1] <- NA
table(dren[2:3])
levels(dren$p_permeabilidad) <- c(2,4,3,5,1,7,6)
levels(dren$p_escurrimiento) <- c(1,3,4,2,6,5)
#dren <- dren[complete.cases(dren),1:3]
levels(dren$p_permeabilidad) <- c(2,3,2,4,1,4,4)  
levels(dren$p_escurrimiento) <- c(1,2,3,1,4,4)
table(dren[,2:3])
dren$dren <- as.factor(paste(dren$p_permeabilidad,dren$p_escurrimiento,sep = "."))
summary(dren)

d1.1 <- base::merge(x = d1[,c(1:34,42:45)], y = dren, by ="id.p", all.x = TRUE)


#compute soil properties per id.hor
name(d1)
#n(d1$id.hor)
d2 <- cbind(d1[,c(1:10)],d1[,11:34],d1[,c(36:45)])
#d2[,26:33][is.na(d2[,26:33])] <- 0
#### merge horizons
## horizon boundaries

# limits <- cbind(ddply(d1,.(id.hor), summarise, mintop=min(top))[,1:2],
#             maxbot=(ddply(d1,.(id.hor), summarise, maxbot=max(bottom))[,2]))
## soil properties 
names(d2)[11:34]
####  aggregation of horizon by id.hor. Warning! = If one horizon has NA the other horizons result in NA


library(dplyr)
d2.2 <- group_by(.data = d2, id.hor)

summarise(d2.2, n())

d2.4 <- unique(merge(d2.2[,c(1:3,5,8,9,10,43)],summarise(d2.2, n()), by = "id.hor"))
# for(i in 11:33){
#   d2.3 <- d2.2[,c(1,i,45)]
#   names(d2.3)[2] <- "X"
#   d2.4[,i-3] <- summarise(.data = d2.3, weighted.mean(x = X, w = weight, na.rm = TRUE))[,2]
#   names(d2.4)[i-3] <- names(d2.2)[i]
# }
# 
# d3 <-  cbind(ddply(d2,.(id.hor), summarise, a_S=sum(a_sum_bases))[,1:2],
#              a_CEC=(ddply(d2,.(id.hor), summarise, a_CEC=sum(a_CEC))[,2]),
#              a_base_ca=(ddply(d2,.(id.hor), summarise, a_base_ca=sum(a_base_ca))[,2]),
#              a_base_mg=(ddply(d2,.(id.hor), summarise, a_base_mg=sum(a_base_mg))[,2]),
#              a_base_k=(ddply(d2,.(id.hor), summarise, a_base_k=sum(a_base_k))[,2]),
#              a_base_na=(ddply(d2,.(id.hor), summarise, a_base_na=sum(a_base_na))[,2]),
#              a_H=(ddply(d2,.(id.hor), summarise, a_H=sum(a_h))[,2]),
#              a_sat_CEC=(ddply(d2,.(id.hor), summarise, a_sat_CEC=sum(a_saturacion_t))[,2]),
#              a_sat_SH=(ddply(d2,.(id.hor), summarise, a_sat_SH=sum(a_saturacion_s_h))[,2]),
#              a_cond=(ddply(d2,.(id.hor), summarise, a_cond=sum(a_conductividad))[,2]),
#              a_res_pasta=(ddply(d2,.(id.hor), summarise, a_res_pasta=sum(a_resistencia_pasta))[,2]),
#              a_ph_pasta=(ddply(d2,.(id.hor), summarise, a_ph_pasta=sum(a_ph_pasta))[,2]),
#              a_ph_h2o=(ddply(d2,.(id.hor), summarise, a_ph_h2o=sum(a_ph_h2o))[,2]),
#              a_ph_kcl=(ddply(d2,.(id.hor), summarise, a_ph_kcl=sum(a_ph_kcl))[,2]),
#              a_OC=(ddply(d2,.(id.hor), summarise, a_OC=sum(a_OC))[,2]),
#              a_clay=(ddply(d2,.(id.hor), summarise, a_clay=sum(a_arcilla))[,2]),
#              a_silt_20=(ddply(d2,.(id.hor), summarise, a_silt_20=sum(a_limo_2_20))[,2]),
#              a_silt_50=(ddply(d2,.(id.hor), summarise, a_silt_50=sum(a_limo_2_50))[,2]),
#              a_sand_100=(ddply(d2,.(id.hor), summarise, a_sand_100=sum(a_arena_muy_fina))[,2]),
#              a_sand_250=(ddply(d2,.(id.hor), summarise, a_sand_250=sum(a_arena_fina))[,2]),
#              a_sand_500=(ddply(d2,.(id.hor), summarise, a_sand_500=sum(a_arena_media))[,2]),
#              a_sand_1k=(ddply(d2,.(id.hor), summarise, a_sand_1k=sum(a_arena_gruesa))[,2]),
#              a_sand_2k=(ddply(d2,.(id.hor), summarise, a_sand_2k=sum(a_arena_muy_gruesa))[,2]),
#              a_caco3=(ddply(d2,.(id.hor), summarise, a_caco3=sum(a_ca_co3))[,2])
#              )

## merge (((limits + d2) + d3) + d2) #Concretions and mottles remain out of this dataset

#d4 <- merge(x = d2.4, y = dren, by = "id.p")

# d4 <- merge(x= merge(x= unique(merge(x = limits,y = d2.4[,c(1:3,5,6,8)],by = "id.hor",all.x = F, all.y = T)), 
#               y= d3, by= "id.hor", all= T),
#               y= unique(d2[,c(5,35:40)]), by= "id.p", all.x=T, all.y=F)
#write.csv(d4, "d4.csv")

# to recover concretions and mottles
# d2$moteados[(d2$moteados)==""]<-NA
# d2$is.mottles<- as.numeric(!is.na(d2$moteados))
#d2$is.mottles[d2$is.mottles==0] <-9999
#d2$is.mottles[(d2$is.mottles)<9999]<- d2$top[(d2$is.mottles)<9999]

# d2$concreciones[(d2$concreciones)==""]<-NA
# d2$is.concr<- as.numeric(!is.na(d2$concreciones))
#d2$is.concr[d2$is.concr==0] <-9999
#d2$is.concr[d2$is.concr<9999]<- d2$top[d2$is.concr<9999]
# merge d4 + concretions(depth) + mottles(depth)
# d5 <-merge(x=merge(d4, ddply(d2,.(id.p), summarise, is.mottles=min(is.mottles)),by= "id.p", all=T),
#       y=ddply(d2,.(id.p), summarise, is.concr=min(is.concr)), by= "id.p")
# d5 <-merge(x=merge(d4, ddply(d2,.(id.p), summarise, is.mottles=max(is.mottles)),by= "id.p", all=T),
#            y=ddply(d2,.(id.p), summarise, is.concr=max(is.concr)), by= "id.p")

# d5$is.concr[d5$is.concr==9999] <-NA
# d5$is.mottles[d5$is.mottles==9999] <-NA
#name(d4)
# d5 <- d5[,-c(33,34)]
# d5 <- d5[,c(1:6,8,9,10,19,21:25,27)]
d2.2$id.p.h <- paste(d2.2$id.p,"_", d2.2$horizon, sep = "")
name(d2.2)
d2.5 <- d2.2[,c(45,2,3,5,8,9,10,43,11,12,14:16,25,26,27,29)]
names(d2.5) <- c("id.p.h","X", "Y", "id.p", "hor", "top", "bottom",
                 "thick", "tb", "CEC", "phw", "phkcl", "resist", 
                 "OC", "clay", "silt20", "sand.mf")
d2.6 <- merge(d2.5,dren, by= "id.p", all.x = TRUE)
# ### order variables by horizons
# A <- unique(d2.5[d2.5$hor=="A", c(3,5:15)])
# B <- unique(d2.5[d2.5$hor=="B", c(3,5:15)])
# E <- unique(d2.5[d2.5$hor=="E", c(3,5:15)])
# BC <- unique(d2.5[d2.5$hor=="BC", c(3,5:15)])
# C <-  unique(d2.5[d2.5$hor=="C", c(3,5:15)])
# 
# names(A)[2:12] <- paste(names(A)[2:12],".A",sep="")
# names(B)[2:12] <- paste(names(B)[2:12],".B",sep="")
# names(E)[2:12] <- paste(names(E)[2:12],".E",sep="")
# names(BC)[2:12] <- paste(names(BC)[2:12],".BC",sep="")
# names(C)[2:12] <- paste(names(C)[2:12],".C",sep="")
# 
# d2.6 <- unique(d2.5[,1:3])
# d2A <- merge(d2.6, d2.5[d2.5$hor == "A",3:4], by = "id.p", all.x = TRUE)
# d2A <- merge(d2.6, A, by= "id.p", all.x = TRUE)
# AB <-merge(d2A,B, by="id.p", all=T)
# ABE <- merge(AB,E, by= "id.p", all=T)
# ABEBC <- merge(ABE, BC, by= "id.p", all=T) 
# ABEBCC <- merge(ABEBC, C, by= "id.p", all=T) 
# d6 <-  ABEBCC
# #################### THE END OF RE-SHAPING
# rm(list=ls()[ls()!="d6"])
# ################## preparing endogenous variables from conceptual model #####
# name <- function(x) { as.data.frame(names(x))} # as.data.frame(.Primitive("names"))
# name(d6)
# corrgram(d6)
# # PSI 
# d6$esp.A <- d6$a_base_na.A/d6$a_S.A*100
# d6$esp.B <- d6$a_base_na.B/d6$a_S.B*100
# # Bt ratio
# d6$Bt <-d6$a_clay.B/d6$a_clay.A
# # is.caco3
# d6$a_caco3.A[is.na(d6$a_caco3.A)==T]<-0
# d6$a_caco3.B[is.na(d6$a_caco3.B)==T]<-0
# d6$a_caco3.BC[is.na(d6$a_caco3.BC)==T]<-0
# d6$is.caco3 <- as.numeric(d6$a_caco3.A>0 | d6$a_caco3.B>0 | d6$a_caco3.BC>0)
# # depth caco3
# d6$d.caco3[d6$a_caco3.BC>0] <- d6$mintop.BC[d6$a_caco3.BC>0]
# d6$d.caco3[d6$a_caco3.B>0] <- d6$mintop.B[d6$a_caco3.B>0]
# d6$d.caco3[d6$a_caco3.A>0] <- d6$mintop.A[d6$a_caco3.A>0]
# # is.E
# d6$is.E <- as.numeric(!is.na(d6$maxbot.E))
# # thick.A
# d6$thick.A <- d6$maxbot.A - d6$mintop.A 

## dataset endogenous veriables
# endo <- d6[,c(1:25,48:58)]
#endo$is.hydro <- as.numeric(endo$is.mottles>0 & endo$is.concr>0)
############################################################################################################
endo <- as.data.frame(d2.6)
#install.packages("maptools")
library(raster)
library(maptools)
library(sp)
library(rgdal)

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
header_m <- c("env", "lstm", "lstsd", "evim", "evisd")

files_n <- list.files(path = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/output/",pattern=".tif$")
header_n <- gsub(pattern = ".tif",replacement = "",x = files_n)
header_n <- paste("X",header_n, sep = "")
files_n <- paste("output/", files_n, sep = "")

files_M <- list.files(path = "/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/MCD43A4/",pattern=".tif$")
header_M <- gsub(pattern = ".tif",replacement = "",x = files_M)
files_M <- paste("MCD43A4/", files_M, sep = "")

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

for(i in 1:length(files_n)) {
  endo@data[,length(endo@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_n[i])
  proj4string(stack[[i]]) <- modis # change projection
  endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
  stack <- list()
  names(endo@data)[length(endo@data)] <- header_n[i]
}

# files MCD43A4
for(i in 1:length(files_M)) {
  endo@data[,length(endo@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_M[i])
  proj4string(stack[[i]]) <- modis # change projection
  endo@data[,length(endo@data)+1] <- over(endo, stack[[i]])[,1]
  stack <- list()
  names(endo@data)[length(endo@data)] <- header_M[i]
}
#image(raster(("mod13q1_tot_mean.tif")))

# clean points out of study area
endo <- spTransform( endo, posgar98)
calib <- as.data.frame(endo)
calib <- calib[!is.na(calib$dem),]
corrgram(calib[calib$hor=="B",c(3:15,83,88,93,98,103,108,113,118)],lower.panel=panel.pie, upper.panel=panel.conf, diag.panel=panel.minmax)
#calib <- calib[!is.na(calib$thick.A),]
#calib <- calib[!is.na(calib$Bt),]
#calib <- calib[!is.na(calib$a_OC.A),]
calib <- calib[,c(1:34,seq(35,79, by = 2), seq(36,80, by = 2),81:122)]
name(calib)
calib$XX1 <- apply(X = calib[,c(35:38)],MARGIN = 1,FUN = mean)
calib$XX2 <- apply(X = calib[,c(42:47)],MARGIN = 1,FUN = mean)
calib$XX3 <- apply(X = calib[,c(50:54)],MARGIN = 1,FUN = mean)
calib$XD1 <- apply(X = calib[,c(58:63)],MARGIN = 1,FUN = mean)
calib$XD2 <- apply(X = calib[,c(66:70)],MARGIN = 1,FUN = mean)
calib$XD3 <- apply(X = calib[,c(73:77)],MARGIN = 1,FUN = mean)

calib <- unique(calib)
calib$top[calib$id.p.h=="350_B1"]  <- 25
calib[calib$id.p==568,1:10][2,4] <- 18
calib[calib$id.p==576,1:10][5,4] <- 42
calib[calib$id.p==710,1:10][7,c(4,5)] <- c(160,180)
calib[calib$id.p==687,1:10][c(5,6),c(4,5)] <- c(145,100,190,145)
calib$hor[calib$hor== "AB|BA"] <- "AB.BA"
name(calib)
calib1 <- calib[,c(1:34,81:128)]
#calib$name <- as.factor(p$hor)

# d <- cbind(d[,c(1:47)], d[,73:75], d[,c(71,72)])
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration")
write.csv(calib1, "calib.data-4.3.csv")

# calib <- group_by(.data = calib, hor)
# summarise(calib[calib$phw<8,], mean(resist, na.rm = TRUE))
# corrgram(calib[,35:57])
# boxplot(calib[calib$phw<8,20:23])
# summary(calib[calib$phw>=8,])
#####################################checking_dataset###########################################################



