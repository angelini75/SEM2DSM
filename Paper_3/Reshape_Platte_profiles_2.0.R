# It comes from Kansas_KS2.R


setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/USDA/")
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

d <- read.csv("NCSS/pedon_platte_extended.csv")#[,-3:-4]
name(d)
D <- d
#d <- read.table("Finnell/Table.txt", header = TRUE, sep = "|")#[,-3:-4]
pedon <- read.table("NCSS/csv/NCSS_Pedon_Taxonomy.csv", sep = ",", header = T)
layer <- read.csv("NCSS/csv/NCSS_Layer.csv", sep = "\t")
CEC <- read.csv("NCSS/csv/CEC_and_Bases.csv", sep = "\t")
carbon <- read.csv("NCSS/csv/Carbon_and_Extractions.csv", sep = "\t")
texture <- read.csv("NCSS/csv/PSDA_and_Rock.csv", sep = "\t")

# pedon$site_key <- as.numeric(pedon$site_key)
# pedon$latitude_decimal_degrees <- as.numeric(as.character(pedon$latitude_decimal_degrees))
# pedon$longitude_decimal_degrees <- as.numeric(as.character(pedon$longitude_decimal_degrees))
# summary(d$Y)
# summary(d$X)

layer$hzn_top <- as.numeric(as.character(layer$hzn_top))
layer$hzn_bot <- as.numeric(as.character(layer$hzn_bot))
d <- d[,c(1:4)]
#d2 <- merge(d,layer[,c(1,2,5,11,12,14,16)], by = "pedon_key", all.x = TRUE)
#evaluation sequence sub master horizons
d2.e <- merge(d,layer[,c(1,2,5,11,12:17)], by = "pedon_key", all.x = TRUE)

# d3 <- merge(d2,CEC[,c(1,16:19)], by = "labsampnum", all.x = T)
# d4 <- merge(d3,carbon[,c(1,4,7)], by = "labsampnum", all.x = T)
# profiles <- merge(d4,texture[,c(1,4,7,8)], by = "labsampnum", all.x = T)
# 
# profiles <- profiles[!is.na(profiles$hzn_top),]
# profiles <- profiles[!(is.na(profiles$cec_nh4) & is.na(profiles$clay_tot_psa)&
#                               is.na(profiles$c_tot)),]
# profiles <- profiles[!(is.na(profiles$cec_nh4) & is.na(profiles$cec_sum)&
#                          is.na(profiles$cec_nhcl)),]
# profiles <- profiles[!(is.na(profiles$c_tot) & is.na(profiles$oc)),]
# profiles <- profiles[!(is.na(profiles$clay_tot_psa) & is.na(profiles$clay_f)),]

# to be changed when analysis of horizons is finished
d3 <- merge(d2.e,CEC[,c(1,16:19)], by = "labsampnum", all.x = T)
d4 <- merge(d3,carbon[,c(1,4,7)], by = "labsampnum", all.x = T)
profiles <- merge(d4,texture[,c(1,4,7,8)], by = "labsampnum", all.x = T)
profiles <- profiles[profiles$hzn_master != "",]
paste("How many soil profiles there are in the area? ->",
      length(unique(profiles$pedon_key)))

A <- unique(profiles$pedon_key[profiles$hzn_master == "A"])
paste("How many have A? ->", length(unique(A)))

B <- unique(profiles$pedon_key[profiles$hzn_master == "B"])
paste("How many have B? ->",length(unique(B)))

C <- unique(profiles$pedon_key[profiles$hzn_master == "C"])
paste("How many have C? ->", length(unique(C)))
abc <- as.data.frame(table(c(A,B,C)))
abc <- as.numeric(as.character(abc$Var1[abc$Freq==3]))
paste("How many have A, B and C? ->", length(abc))

hz <- c(A,B,C)[-1239]
hz <- as.data.frame(table(hz))
hz <- as.numeric(as.character(hz$hz[hz$Freq==3]))

profiles <- profiles[which(profiles$pedon_key %in% hz),]
#write.csv(profiles, "~/Documents/platte_area.csv")
################################################################################
# profiles <- profiles[!(is.na(profiles$cec_nh4) &
#                          is.na(profiles$clay_tot_psa) &
#                          is.na(profiles$c_tot)),]
# profiles <- profiles[!(is.na(profiles$cec_nh4) & 
#                         is.na(profiles$cec_sum) &
#                          is.na(profiles$cec_nhcl)),]
# profiles <- profiles[!(is.na(profiles$c_tot) & is.na(profiles$oc)),]
# profiles <- profiles[!(is.na(profiles$clay_tot_psa) & is.na(profiles$clay_f)),]

# comparison between soil properties
plot(profiles$c_tot~profiles$oc)
plot(profiles$clay_tot_psa~profiles$clay_f)
plot(profiles$cec_nh4~profiles$cec_sum)
name(profiles)
# [1] "labsampnum"   "pedon_key"    "X"            "Y"            "hzn_master"  
# [6] "hzn_desgn"    "cec_sum"      "cec_nh4"      "c_tot"        "oc"          
# [11] "clay_tot_psa" "clay_f" 
#profiles.r <- profiles[,c(1:4,10,9,7,8,11,12,15,16,17,18)]
profiles.e <- profiles[,c(1:4,9:13,7,8,14,15,18,19,20)]
summary(profiles.e)

# # copy values from c_tot to oc where oc == NA
# profiles.r$oc[is.na(profiles.r$oc)] <- profiles.r$c_tot[which(is.na(profiles.r$oc))]
# profiles.r <- profiles.r[,-11]
# copy values from c_tot to oc where oc == NA
profiles.e$oc[is.na(profiles.e$oc)] <- profiles.e$c_tot[which(is.na(profiles.e$oc))]
profiles.e <- profiles.e[,-14]

# copy values from cec_tot to cec_nh4 where cec_nh4 == NA
profiles.e$cec_nh4[is.na(profiles.e$cec_nh4)] <- 
  profiles.e$cec_sum[which(is.na(profiles.e$cec_nh4))]
profiles.e <- profiles.e[,c(-12)]

# names(profiles.r) <- c("labsampnum","idp", "X", "Y","hzn","hzn_nom","top", "bot", "cec", "oc", "clay")
# summary(profiles.r)
names(profiles.e) <- c("labsampnum","idp", "X", "Y","hzn_old", "hzn_nom", "hzn_disc",
                       "hzn", "hzn_prime","top", "bot", "cec", "oc", "clay")
summary(profiles.e)

# adding 30 cm to bottom horizons where bot == NA (2 cases)
# profiles.r[profiles.r$idp == profiles.r$idp[which(is.na(profiles.r$bot))],]
# profiles.r$bot[which(is.na(profiles.r$bot))] <- c(234, 200)

multigenetic <- unique(profiles.e$idp[profiles.e$hzn_disc != ""])
paste("How many profiles are multigenetic? ->", length(unique(multigenetic)))

#pe <- profiles.e[which(!(profiles.e$idp %in% multigenetic)),]

pe <- profiles.e[profiles.e$hzn_disc == "",] # remove discontinued horizons
pe <- pe[pe$hzn == "A" |
           pe$hzn == "B"|
           pe$hzn == "C",]

library(sp)
library(maptools)
coordinates(pe) <- ~X+Y
NE <- readShapePoly("Finnell/SoilMaps/NE_ext.shp")
KS <- readShapePoly("Finnell/SoilMaps/KS_ext.shp")

NE@data <- NE@data[,c(-1:-2,-5:-6)]
pe@data <- cbind(pe@data, over(pe,NE))
KS@data <- KS@data[,c(-1:-2,-5:-6)]
pe@data <- cbind(pe@data, over(pe,KS))
pe <- as.data.frame(pe)
pe$MUSYM <- as.numeric(as.character(pe$MUSYM))
pe$MUSYM.1 <- as.numeric(as.character(pe$MUSYM.1))
pe$MUKEY <- as.numeric(as.character(pe$MUKEY))
pe$MUKEY.1 <- as.numeric(as.character(pe$MUKEY.1))

pe$MUSYM[is.na(pe$MUSYM)] <- pe$MUSYM.1[is.na(pe$MUSYM)]
pe$MUKEY[is.na(pe$MUKEY)] <- pe$MUKEY.1[is.na(pe$MUKEY)]
name(pe)
pe <- pe[,c(-17,-18)]
sp <- read.table(file = "Finnell/Table.txt", header = F, sep = "|")
names(sp) <- names(read.table(file = "Finnell/Table_Nov.txt", 
                              header = T, sep = "|"))
name(sp)
sp <- sp[,c(1,2,6,16:18,33,36,39)]

View(per[!(complete.cases(per[,11:13])),])
pe <- pe[pe$labsampnum != "40A14293",]
pe <- pe[pe$labsampnum != "40A14461",]
pe <- pe[pe$labsampnum != "40A14707",]
pe <- pe[pe$labsampnum != "40A14716",]

# first asumtion OC C horizon == median(pe$oc[pe$hzn == "C"])
pe$oc[pe$hzn=="C" & is.na(pe$oc)] <- median(pe$oc[pe$hzn=="C"], na.rm=T)

# remove profiles with too many NA's
per <- pe[pe$idp != 1833,]
per <- per[per$idp != 1837,]
per <- per[per$idp != 1838,]
per <- per[per$idp != 1845,]
per <- per[per$idp != 1846,]
per <- per[per$idp != 2053,]
per <- per[per$idp != 2117,]
`2125` <- per[per$idp == 2125,][-10:-12,]
per <- per[per$idp != 2125,]
per <- rbind(per,`2125`)
per <- per[!(per$idp == 2126 & per$hzn_nom == "C3"),]
per <- per[per$idp != 2165,]
per <- per[per$idp != 2166,]
per <- per[!(per$idp == 2169 & (per$hzn_nom == "Ap" | 
                                  per$hzn_nom == "B" |
                                  per$hzn_nom == "Cr" )),]
per <- per[per$idp != 7661,]
per <- per[per$idp != 7662,]
per <- per[per$idp != 14853,]
per <- per[per$idp != 14854,]
per <- per[per$idp != 17106,]
per <- per[per$idp != 23819,]
per <- per[per$idp != 23903,]
per <- per[per$idp != 23904,]

paste("how many profiles are in per? ->", length(unique(per$idp)))

A <- unique(per$idp[per$hzn == "A"])
B <- unique(per$idp[per$hzn == "B"])
C <- unique(per$idp[per$hzn == "C"])

abc <- as.data.frame(table(c(A,B,C)))
abc <- as.numeric(as.character(abc$Var1[abc$Freq==3]))
paste("how many profiles of per have ABC? ->", length(unique(abc)))

p <- per[which(per$idp %in% abc),]
rownames(p) <- 1:length(p[,1])
p <- p[c(-257:-261) ,]
rownames(p) <- 1:length(p[,1])

#View(p[which(p$idp %in% unique(p$idp[is.na(p$cec)])),])
# mu <- unique(p$MUKEY[which(p$idp %in% unique(p$idp[is.na(p$cec)]))])
# str(sp)
sp <- sp[with(sp, order(musym, mukey, -comppct_r,hzdept_r)), ]
rownames(sp) <- 1:length(sp[,1]) 
# sp <- sp[which(sp$mukey %in% mu),]
# sp <- sp[with(sp, order(musym, mukey, -comppct_r,hzdept_r)), ]

#write.csv(p, "p.csv")
#write.csv(sp, "sp.csv")
p <- read.csv("~/Documents/SEM2DSM1/Paper_3/data/p.csv")[,-1]
sp <- read.csv("~/Documents/SEM2DSM1/Paper_3/data/sp.csv")[,-1]

paste(length(unique(p$idp[is.na(p$cec)])), "profiles without CEC")
paste("how many profiles are in p? ->", length(unique(p$idp)))

data <- p[which(p$idp %in% unique(p$idp[!is.na(p$cec)])),]
data <- data[which(data$idp %in% unique(data$idp[!is.na(data$oc)])),]
paste("how many profiles are in data? ->", length(unique(data$idp)))
sp[sp$mukey == 1691123,]

# one replacement
#data$cec[is.na(data$cec)] <- sp[51464,9] # CEC extracted from the map

# 158 soil profiles
# write.csv(data, "data_KSNE.csv")

################################################################################
###                                  Start here                              ###
################################################################################
setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/USDA/")
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

# d <- read.csv("data_KSNE.csv")[,-1]
d <- data
d$thick <- d$bot - d$top
name(d)
d <- d[,c(1,2,3,4,8,10,11, 17, 12:14)]
names(d)[1] <- "idh"
#compute soil properties per id.hor
library(plyr)
library(reshape2)

d1 <- data.frame(idp = unique(d$idp))
d2 <- merge(d1, dcast(data = ddply(d, 
                                   .(idp,hzn),
                                   summarise,
                                   clay = weighted.mean(x = clay,
                                                        w =  thick,
                                                        na.rm = TRUE)),
                      idp ~ hzn), by = "idp")
names(d2)[2:4] <- paste0("clay.",names(d2)[2:4])

d2 <- merge(d2, dcast(data = ddply(d, 
                                   .(idp,hzn),
                                   summarise,
                                   cec = weighted.mean(x = cec,
                                                       w =  thick,
                                                       na.rm = TRUE)),
                      idp ~ hzn), by = "idp")
names(d2)[5:7] <- c("cec.A", "cec.B", "cec.C")

d2 <- merge(d2, dcast(data = ddply(d, 
                                   .(idp,hzn),
                                   summarise,
                                   oc = weighted.mean(x = oc,
                                                      w =  thick,
                                                      na.rm = TRUE)),
                      idp ~ hzn), by = "idp")
names(d2)[8:10] <- c("oc.A", "oc.B", "oc.C")

d3 <- merge(unique(d[,c("idp", "X", "Y")]), d2, by = "idp")

############################################################
D <- d3
D <- D[complete.cases(D),]
#############################################################
setwd("/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/")
library(raster)

## Points over DEM and its derivates
files <- list.files(pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("dem", "twi", "vdchn") 

#define crs
wgs84 <- CRS("+init=epsg:4326")
#UTM14N <- CRS("+init=epsg:32614")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
NAD83.KS.N <- CRS("+init=epsg:2796")


#devtools::install_github("environmentalinformatics-marburg/mapview", ref = "develop")
library(mapview)
coordinates(D) <- ~X+Y
proj4string(D) <- wgs84
#mapview(s, burst = TRUE)
mapview::viewExtent(D)
map <- mapview(D, map = NULL,
               map.types = mapviewGetOption("basemaps"), zcol = NULL, burst = TRUE,
               color = c("#FF00FF", "#0489B1"), alpha = 0.6, alpha.regions = 0.2,
               na.color = mapviewGetOption("na.color"), at = NULL, cex = 5, lwd = 2,
               popup = popupTable(D), legend = mapviewGetOption("legend"),
               legend.opacity = 1, layer.name = deparse(
                 substitute(x, env = parent.frame())), 
               verbose = mapviewGetOption("verbose"),
               homebutton = TRUE)
map + viewExtent(D) + raster(files[1])


# assign projection
proj4string(D) <- wgs84
D <- spTransform(D, NAD83.KS.N)

# profiles.e@data[,9+seq_along(files)] <- NA
# names(profiles.r@data)[10:19] <- header
# 
# for(i in seq_along(files)){
#   profiles.r@data[,9+i] <- extract(x = raster(files[i]), y = profiles.r) 
# }

D@data[,10+seq_along(files)] <- NA
names(D@data)[11:13] <- header

for(i in seq_along(files)){
  D@data[,10+i] <- extract(x = raster(files[i]), y = D) 
}

# list of tif files (modis)
files <- list.files(pattern = ".tif$")
header <- gsub(".tif", "", files)
header <-  c("evim", "evisd", "lstm", "lstsd", "ndwi.a", "ndwi.b") 

# transform projection
D <- spTransform(x = D, modis)

D@data[,13+seq_along(files)] <- NA
names(D@data)[14:19] <- header

# extract values from tiff files
for(i in seq_along(files)){
  D@data[,13+i] <- extract(x = raster(files[i]), y = D)
}
D <- spTransform(x = D, NAD83.KS.N)
zerodist(D, zero=0.0)
D <- as.data.frame(D)
D <- D[c(-101,-102,-104),]
setwd("~/Documents/SEM2DSM1/Paper_3/data")
write.csv(D,"KS.data-0.3.csv") 

#########################
# Argentinian data
library(ggplot2)
setwd("~/Documents/SEM2DSM1/Paper_2/data/")

d <- read.csv("calib.data-5.0.csv")[,c(-1,-20)] #remove water variable 
d <- d[,c(-12:-16,-20,-21,-25)]

meltp <- melt(unique(d,id.vars = c("id.p")))

ggplot(data = meltp,
       aes(x = value, fill=variable)) + geom_histogram() + 
  facet_wrap( ~ variable, scales = "free_x")
d$H <- NA
e <- data.frame(d[,c(1,2,5,8,11:20)])
e$H <- "A"
names(e)[2:4] <- c("CEC.B","OC.B","clay.B")
e <- rbind(e,d[,c(1,3,6,9,11:20)])
e$H[is.na(e$H)] <- "B"

names(e)[2:4] <- c("CEC.C","OC.C","clay.C")
e <- rbind(e,d[,c(1,4,7,10,11:20)])
e$H[is.na(e$H)] <- "C"
names(e)[2:4] <- c("CEC","OC","Clay")
### 
# US data

# statistics
library(ggplot2)
library(reshape2)

D <- as.data.frame(D)
head(D,20)

meltp <- melt(D,id.vars = c("idp"))


ggplot(data = meltp, aes(x = value)) + 
  geom_density(alpha = 0.2) + 
  facet_wrap( ~ variable,scales = "free")


ggplot(data = unique(melt.sp[!(melt.sp$variable == "CEC" |
                                 melt.sp$variable == "OC" |
                                 melt.sp$variable == "Clay"),]),
       aes(x = value, fill = country)) + geom_density(alpha = 0.4) + 
  facet_wrap( ~ variable,scales = "free")


ggplot(data = meltp[meltp$variable == "top" |
                        meltp$variable == "bot",],
       aes(x = value, fill = variable)) + geom_density(alpha = 0.4) + 
  facet_wrap( ~ hzn_key,scales = "free")


# ggplot(data = melt.sp[melt.sp$variable == "CEC" |
#                         melt.sp$variable == "OC" |
#                         melt.sp$variable == "Clay",],
#        aes(x = value, fill = country)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable+H,scales = "free")
# ggplot(data = unique(melt.sp[!(melt.sp$variable == "CEC" |
#                         melt.sp$variable == "OC" |
#                         melt.sp$variable == "Clay"),]),
#        aes(x = value, fill = country)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable,scales = "free")
# 
# 
# length(unique(profiles.r[,c(2)]))
# summary(unique(profiles.r[,c(2,20:26)]))
# name(profiles.r)
# meltp <- melt(unique(profiles.r[,c(2,20:26)]),id.vars = c("idp"))
# 
# ggplot(data = meltp,
#        aes(x = value, fill=variable)) + geom_histogram() + 
#   facet_wrap( ~ variable, scales = "free_x")
# 
# us <- profiles.r[profiles.r$hzn == "A" |
#                    profiles.r$hzn == "B" |
#                    profiles.r$hzn == "C", c(2:3,7:9,20:28)]
# us$country <- "US"
# e$country <- "Arg"
# names(us)[c(1,2,3,4,5,7,8)] <- c("id.p","H","CEC","OC","Clay","twi","vdchn")
# s <- rbind(e,us[,c(1,3:8,10,9,11:14,2,15)])
# melt.sp <- melt(s, id.vars = c("id.p","H", "country"))
# 
# ggplot(data = melt.sp[melt.sp$variable == "CEC" |
#                         melt.sp$variable == "OC" |
#                         melt.sp$variable == "Clay",],
#        aes(x = value, fill = country)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable+H,scales = "free")
# ggplot(data = unique(melt.sp[!(melt.sp$variable == "CEC" |
#                                  melt.sp$variable == "OC" |
#                                  melt.sp$variable == "Clay"),]),
#        aes(x = value, fill = country)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable,scales = "free")


#### THIS SECTION WAS USEFUL FOR TESTING MLR APPROACH ####
# BEGIN NOT RUN
# A0 <- profiles.r[profiles.r$top<5 & profiles.r$bot>5,]
# B70 <- profiles.r[profiles.r$top<70 & profiles.r$bot>70,]
# C150 <- profiles.r[profiles.r$top<200 & profiles.r$bot>150,]
# 
# step(lm(oc ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, A0), direction = "both")
# summary(lm(formula = oc ~ vdchn + evisd + ndwi.a + ndwi.b + X + Y + chnbl + 
#              EVI_M_JanFeb_250 + LS + rsp + dem, data = A0))
# summary(lm(formula = oc ~ dem + twi + vdchn + evisd + lstm + 
#              ndwi.a + ndwi.b + X + Y, data = A0))
# 
# step(lm(cec ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, A0), direction = "both")
# summary(lm(formula = cec ~ vdchn + evisd + lstm + ndwi.a + chnbl + EVI_M_JanFeb_250 + 
#              LS + srtm250, data = A0))
# 
# step(lm(clay ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, A0), direction = "both")
# summary(lm(formula = clay ~ vdchn + lstm + ndwi.a + chnbl + EVI_M_JanFeb_250 + 
#              EVI_SD_JanFeb_250 + LS + srtm250, data = A0))
# #------------------==============----------------------===============#
# 
# step(lm(oc ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, B70), direction = "both")
# summary(lm(formula = oc ~ dem + vdchn + ndwi.a + chnbl + EVI_M_JanFeb_250 + 
#              rsp + slope, data = B70))
# 
# step(lm(cec ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, B70), direction = "both")
# summary(lm(formula = cec ~ dem + vdchn + evisd + ndwi.a + ndwi.b + X + 
#              chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + rsp + Valley.Depth, 
#            data = B70))
# 
# step(lm(clay ~ dem + twi + vdchn + evisd + lstm + 
#           ndwi.a + ndwi.b + X + Y + chnbl + EVI_M_JanFeb_250 + 
#           EVI_SD_JanFeb_250 + LS + rsp + slope + srtm250 + twi + 
#           Valley.Depth + vdchn, B70), direction = "both")
# summary(lm(formula = clay ~ dem + twi + vdchn + evisd + lstm + ndwi.a + 
#              ndwi.b + Y + chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + 
#              rsp, data = B70))
# #---------------===================------------===============---------#
# # step(lm(oc ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
# #           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
# #           X + Y, C150), direction = "both")
# summary(lm(formula = oc ~ chnbl + srtm250 + twi + Valley.Depth + vdchn + 
#              Y, data = C150))
# 
# # step(lm(cec ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
# #           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
# #           X + Y, C150), direction = "both")
# summary(lm(formula = cec ~ chnbl + rsp + srtm250 + Valley.Depth + vdchn + 
#              X + Y, data = C150))
# # 
# # step(lm(clay ~ chnbl + EVI_M_JanFeb_250 + EVI_SD_JanFeb_250 + LS + 
# #           rsp + slope + srtm250 + twi + Valley.Depth + vdchn + 
# #           X + Y, C150), direction = "both")
# (summary(lm(formula = clay ~ chnbl + srtm250 + twi + Valley.Depth + vdchn + 
#              X, data = C150)))#$r.squared
# END NOT RUN 

library(mapview)


coordinates(D) <- ~X+Y
proj4string(D) <- wgs84
D <- spTransform(D, wgs84)
#mapview(s, burst = TRUE)
mapview::viewExtent(D)
map <- mapview(D, map = NULL,
        map.types = mapviewGetOption("basemaps"), zcol = NULL, burst = TRUE,
        color = c("#FF00FF", "#0489B1"), alpha = 0.6, alpha.regions = 0.2,
        na.color = mapviewGetOption("na.color"), at = NULL, cex = 5, lwd = 2,
        popup = popupTable(D), legend = mapviewGetOption("legend"),
        legend.opacity = 1, layer.name = deparse(
          substitute(x, env = parent.frame())), 
        verbose = mapviewGetOption("verbose"),
        homebutton = TRUE)
map + viewExtent(D)
mapview::mapshot(map,file = "actual.pdf")
mapview::viewExtent(D)

# Alternative
# http://rgraphgallery.blogspot.nl/2013/04/rg68-get-google-map-and-plot-data-in-it.html
library(ggmap)
dhanmap5 = get_googlemap(location = c(lon = -95.5, lat = 38.5), zoom = 8, 
                   maptype = 'roadmap', source = "osm")
dhanmap5 = ggmap(dhanmap5)
# data
myd <- spTransform(D, wgs84)
myd <-  as.data.frame(myd)
myd
# the bubble chart
library(grid)
dhanmap5 +   
  geom_point(
    aes(x = X, y = Y, colour = oc.A),
    data = myd) 


s <- as.data.frame(s)
length(unique(s[s$hzn_c=="c",1:2])[,1])
length(unique(s[s$hzn_c=="a",1:2])[,1])

## AQP FOR PLOT HORIZON FREQUENCY ####
#install.packages('aqp', repos="http://R-Forge.R-project.org")
library(aqp)
names(profiles.e)# <- "name"
s <- profiles.e
names(s)[6] <- "name"
s$name <- as.character(s$name)
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

