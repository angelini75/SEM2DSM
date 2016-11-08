setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/USDA/NCSS/")
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

d <- read.csv("KS_3a.csv")[,-5:-6]
names(d)
#pedon <- read.table("csv/NCSS_Pedon_Taxonomy.csv", sep = "\t", header = T)
layer <- read.csv("csv/NCSS_Layer.csv", sep = "\t")
CEC <- read.csv("csv/CEC_and_Bases.csv", sep = "\t")
carbon <- read.csv("csv/Carbon_and_Extractions.csv", sep = "\t")
texture <- read.csv("csv/PSDA_and_Rock.csv", sep = "\t")

# pedon$site_key <- as.numeric(pedon$site_key)
# pedon$latitude_decimal_degrees <- as.numeric(as.character(pedon$latitude_decimal_degrees))
# pedon$longitude_decimal_degrees <- as.numeric(as.character(pedon$longitude_decimal_degrees))
# summary(d$Y)
# summary(d$X)

layer$hzn_top <- as.numeric(as.character(layer$hzn_top))
layer$hzn_bot <- as.numeric(as.character(layer$hzn_bot))
d2 <- merge(d,layer[,c(1,2,5,11,12,14,16)], by = "pedon_key", all.x = TRUE)

d3 <- merge(d2,CEC[,c(1,16:19)], by = "labsampnum", all.x = T)
d4 <- merge(d3,carbon[,c(1,4,7)], by = "labsampnum", all.x = T)
profiles <- merge(d4,texture[,c(1,4,7,8)], by = "labsampnum", all.x = T)

profiles <- profiles[!is.na(profiles$hzn_top),]
profiles <- profiles[!(is.na(profiles$cec_nh4) & is.na(profiles$clay_tot_psa)&
                              is.na(profiles$c_tot)),]
profiles <- profiles[!(is.na(profiles$cec_nh4) & is.na(profiles$cec_sum)&
                         is.na(profiles$cec_nhcl)),]
profiles <- profiles[!(is.na(profiles$c_tot) & is.na(profiles$oc)),]
profiles <- profiles[!(is.na(profiles$clay_tot_psa) & is.na(profiles$clay_f)),]

# comparison between soil properties
plot(profiles$c_tot~profiles$oc)
plot(profiles$clay_tot_psa~profiles$clay_f)
plot(profiles$cec_nh4~profiles$cec_sum)
name(profiles)
# [1] "labsampnum"   "pedon_key"    "X"            "Y"            "hzn_master"  
# [6] "hzn_desgn"    "cec_sum"      "cec_nh4"      "c_tot"        "oc"          
# [11] "clay_tot_psa" "clay_f" 
D <- profiles[,c(1:4,10,9,7,8,11,12,15,16,17,18)]
summary(D)

# copy values from c_tot to oc where oc == NA
D$oc[is.na(D$oc)] <- D$c_tot[which(is.na(D$oc))]
D <- D[,-11]
# copy values from cec_tot to cec_nh4 where cec_nh4 == NA
D$cec_nh4[is.na(D$cec_nh4)] <- D$cec_sum[which(is.na(D$cec_nh4))]
D <- D[,c(-9,-13)]


names(D) <- c("labsampnum","idp", "X", "Y","hzn","hzn_nom","top", "bot", "cec", "oc", "clay")
D <- D[!is.na(D$labsampnum),]
summary(D)

# checking layer depth
D$hzn <- as.character(D$hzn)
D$labsampnum <- as.character(D$labsampnum)
D[which(
  D$idp %in% 
    D[
      which(
        D$bot-D$top < 1
      ),
      2]
),
]

D[D$labsampnum=="90P00484",5] <- "BC"
D <- D[!(D$idp == 10174 & D$labsampnum == "83P01445"),]
D <- D[!(D$idp == 26252 | D$idp == 26253 | D$idp == 26254 | D$idp == 26255 | 
           D$idp == 4496 | D$idp == 15608 | D$idp == 15608),]
D <- D[!(D$hzn==""),]

D[which(D$idp %in% D[D$hzn=="","idp"]),] # should be zero

summary(D)

## extarct values from rasters (covariates)
name(D)
D.sp <- unique(D[,c(2,3,4)])
library(raster)
library(maptools)
library(sp)

# list of sdat files (DEM)
files <- 
  list.files(
    path = "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
    pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("dem", "twi", "vdchn") 

coordinates(D.sp) <- ~X+Y
#define crs
wgs84 <- CRS("+init=epsg:4326")
UTM14N <- CRS("+init=epsg:32614")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(D.sp) <- wgs84
D.sp <- spTransform(x = D.sp, UTM14N)

D.sp@data[,1+seq_along(files)] <- NA
names(D.sp@data)[2:4] <- header

# extract values from sdat files
for(i in seq_along(files)){
  D.sp@data[,1+i] <-
    extract(
      x = raster(
        paste0(
          "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
          files[i]
          )
        ),
      y = D.sp)
}

# list of tif files (modis)
files <- 
  list.files(
    path = "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
    pattern = ".tif$")
header <- gsub(".sdat", "", files)
header <-  c("evisd", "lstm", "ndwi.a", "ndwi.b") 

# transform projection
D.sp <- spTransform(x = D.sp, modis)

D.sp@data[,4+seq_along(files)] <- NA
names(D.sp@data)[5:8] <- header

# extract values from tiff files
for(i in seq_along(files)){
  D.sp@data[,4+i] <-
    extract(
      x = raster(
        paste0(
          "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
          files[i]
        )
      ),
      y = D.sp)
}
##################################
setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/covar")

files <- list.files(pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("chnbl", "EVI_M_JanFeb_250", "EVI_SD_JanFeb_250", "LS", "rsp",
             "slope", "srtm250", "twi", "Valley Depth", "vdchn") 

# assign projection
D.sp <- spTransform(D.sp, UTM14N)

D.sp@data[,8+seq_along(files)] <- NA
names(D.sp@data)[9:18] <- header

for(i in seq_along(files)){
  D.sp@data[,8+i] <- extract(x = raster(files[i]), y = D.sp) 
}
#####################################
D.sp <- as.data.frame(D.sp)
name(D.sp)

D <- merge(D,D.sp[,c(-19,-20)], by = "idp")

head(D,20)


A0 <- D[D$top<15 & D$bot>15,]
B70 <- D[D$top<70 & D$bot>70,]
C150 <- D[D$top<150 & D$bot>150,]


# step(lm(oc ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, A0), direction = "both")
summary(lm(formula = oc ~ twi + evisd + ndwi.a + X + Y + chnbl + EVI_M_JanFeb_250 + 
             EVI_SD_JanFeb_250 + LS + srtm250, data = A0))

step(lm(clay ~ dem + twi + vdchn + evisd + lstm + 
          ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
          EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
          Valley.Depth + vdchn, A0), direction = "both")
summary(lm(formula = cec ~ dem + vdchn + evisd + ndwi.a + Y + chnbl + 
             EVI_M_JanFeb_250 + rsp, data = A0))

# step(lm(clay ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, A0), direction = "both")
summary(lm(formula = clay ~ evisd + ndwi.a + ndwi.b + Y + chnbl + EVI_M_JanFeb_250 + 
             rsp + srtm250, data =A0))
#------------------==============----------------------===============#

# step(lm(oc ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, B70), direction = "both")
summary(lm(formula = oc ~ chnbl + rsp + srtm250 + twi + vdchn + X, data = B70))

# step(lm(cec ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, B70), direction = "both")
summary(lm(formula = cec ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + 
             srtm250 + vdchn + X, data = B70))

# step(lm(clay ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, B70), direction = "both")
summary(lm(formula = clay ~ chnbl + EVI_SD_JanFeb_250 + rsp + srtm250 + 
              twi + vdchn + X + Y, data = B70))
#---------------===================------------===============---------#
# step(lm(oc ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, C150), direction = "both")
summary(lm(formula = oc ~ chnbl + srtm250 + twi + Valley.Depth + vdchn + 
             Y, data = C150))

# step(lm(cec ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, C150), direction = "both")
summary(lm(formula = cec ~ chnbl + rsp + srtm250 + Valley.Depth + vdchn + 
             X + Y, data = C150))
# 
# step(lm(clay ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
#           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
#           X + Y, C150), direction = "both")
(summary(lm(formula = clay ~ chnbl + srtm250 + twi + Valley.Depth + vdchn + 
             X, data = C150)))#$r.squared


## AQP ## 
#install.packages('aqp', repos="http://R-Forge.R-project.org")
library(aqp)
names(D)[4] <- "name"
s <- D
s$name[s$name== "E" |s$name==  "AB" | s$name== "BA" |s$name==  "EB"] <- "transAB"
s$name[s$name== "BC" |s$name==  "CB"] <- "transBC"
s$name[s$name== "C" |s$name==  "R"] <- "C"
# load sample data and convert into SoilProfileCollection
depths(s) <- idp ~ top + bot
# plot
par(mfrow=c(2,1), mar=c(0, 0, 0, 0))
plot(s[1:20,])
#mtext('SD = 2', side=2, line=-1.5, font=2, cex=0.75)

# aggregate horizonation of simulated data
# note: set class_prob_mode=2 as profiles were not defined to a constant depth
s$name <- as.factor(as.character(s$name))
a <- slab(s, ~ name, class_prob_mode=5)

# convert to long format for plotting simplicity
library(reshape)
a.long <- melt(a, id.vars=c('top','bottom','contributing_fraction'),
               measure.vars=c('A','B','C'))
#a.long <- melt(a, id.vars=c('top','bottom'), measure.vars=levels(s$name)[-1])

# plot horizon probabilities derived from simulated data
# dashed lines are the original horizon boundaries
library(lattice)
xyplot(top ~ value, groups=variable, data=a.long, subset=value > 0,
       ylim=c(305, 0), type=c('l','g'), asp=1.5,
       ylab='Depth (cm)', xlab='Frequency',
       auto.key=list(columns=3, lines=TRUE, points=FALSE),
       panel=panel.depth_function, 
       #prepanel=panel.depth_function,
       cf=a.long$contributing_fraction
       )
#aqp::panel.depth_function()
# # # # # ------------ ----- --- ---- ------------ ------------- --------- -----
library(lattice)
library(grid)

# load sample data, upgrade to SoilProfileCollection
data(sp1)
depths(sp1) <- id ~ top + bottom

# aggregate entire collection with two different segment sizes
a <- slab(sp1, fm = ~ prop)
b <- slab(sp1, fm = ~ prop, slab.structure=5)

# check output
str(a)

# stack into long format
ab <- make.groups(a, b)
ab$which <- factor(ab$which, levels=c('a','b'), 
                   labels=c('1-cm Interval', '5-cm Interval'))

# plot median and IQR
# custom plotting function for uncertainty viz.
xyplot(top ~ p.q50 | which, data=ab, ylab='Depth',
       xlab='median bounded by 25th and 75th percentiles',
       lower=ab$p.q25, upper=ab$p.q75, ylim=c(250,-5),
       panel=panel.depth_function, 
       prepanel=prepanel.depth_function,
       cf=ab$contributing_fraction,
       layout=c(2,1), scales=list(x=list(alternating=1))
)
