setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/USDA/NCSS/")
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}

d <- read.csv("KS_2.csv")#[,-3:-4]
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
profiles.e$cec_nh4[is.na(profiles.e$cec_nh4)] <- profiles.e$cec_sum[which(is.na(profiles.e$cec_nh4))]
profiles.e <- profiles.e[,c(-12)]

# names(profiles.r) <- c("labsampnum","idp", "X", "Y","hzn","hzn_nom","top", "bot", "cec", "oc", "clay")
# summary(profiles.r)
names(profiles.e) <- c("labsampnum","idp", "X", "Y","hzn_old", "hzn_nom", "hzn_disc",
                       "hzn", "hzn_prime","top", "bot", "cec", "oc", "clay")
summary(profiles.e)

# adding 30 cm to bottom horizons where bot == NA (2 cases)
# profiles.r[profiles.r$idp == profiles.r$idp[which(is.na(profiles.r$bot))],]
# profiles.r$bot[which(is.na(profiles.r$bot))] <- c(234, 200)

setwd("/mnt/L0135974_DATA/UserData/BaseARG/study area/USA/covar")
library(raster)

files <- list.files(pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("chnbl", "EVI_M_JanFeb_250", "EVI_SD_JanFeb_250", "LS", "rsp",
             "slope", "srtm250", "twi", "Valley Depth", "vdchn") 

coordinates(profiles.e) <- ~X+Y
#define crs
wgs84 <- CRS("+init=epsg:4326")
UTM14N <- CRS("+init=epsg:32614")

# assign projection
proj4string(profiles.e) <- wgs84
profiles.e <- spTransform(profiles.e, UTM14N)

# profiles.e@data[,9+seq_along(files)] <- NA
# names(profiles.r@data)[10:19] <- header
# 
# for(i in seq_along(files)){
#   profiles.r@data[,9+i] <- extract(x = raster(files[i]), y = profiles.r) 
# }

profiles.e@data[,12+seq_along(files)] <- NA
names(profiles.e@data)[13:22] <- header

for(i in seq_along(files)){
  profiles.e@data[,12+i] <- extract(x = raster(files[i]), y = profiles.e) 
}


###################
# list of sdat files (DEM)
files <- 
  list.files(
    path = "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
    pattern = ".dat$")
header <- gsub(".sdat", "", files)
header <-  c("dem", "twi", "vdchn") 

modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
profiles.e <- spTransform(x = profiles.e, UTM14N)

profiles.e@data[,22+seq_along(files)] <- NA
names(profiles.e@data)[23:25] <- header

# extract values from sdat files
for(i in seq_along(files)){
  profiles.e@data[,22+i] <-
    extract(
      x = raster(
        paste0(
          "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
          files[i]
        )
      ),
      y = profiles.e)
}

# list of tif files (modis)
files <- 
  list.files(
    path = "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
    pattern = ".tif$")
header <- gsub(".sdat", "", files)
header <-  c("evisd", "lstm", "ndwi.a", "ndwi.b") 

# transform projection
profiles.e <- spTransform(x = profiles.e, modis)

profiles.e@data[,25+seq_along(files)] <- NA
names(profiles.e@data)[26:29] <- header

# extract values from tiff files
for(i in seq_along(files)){
  profiles.e@data[,25+i] <-
    extract(
      x = raster(
        paste0(
          "/mnt/L0135974_DATA/UserData/BaseARG/COVARIATES/USA/modelling/",
          files[i]
        )
      ),
      y = profiles.e)
}
#########################
# Argentinian data
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

profiles.e <- as.data.frame(profiles.e)
head(profiles.e,20)

profiles.e$hzn_key <- NA
profiles.e$hzn_nom <- as.character(profiles.e$hzn_nom)
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bt")] <- "BtX"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Btk")] <- "Btk"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bt$")] <- "Bt"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bt1")] <- "Bt1"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bt2")] <- "Bt2"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bt[3-9]")] <- "Bt9"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bw")] <- "Bw"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^B2.")] <- "Bw"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^Bk.")] <- "Bk"
profiles.e$hzn_key[grep(profiles.e$hzn_nom,pattern = "^C")] <- "C"

profiles.e$oc <- log10(profiles.e$oc)

B <- profiles.e$hzn_nom[grep(profiles.e$hzn_nom,pattern = "^B")] 
B <- B[grep(B,pattern = "[^BC].")]
B <- B[grep(B,pattern = "[^Bt]..")]
B <- B[grep(B,pattern = "[^Bt]...")]


length(unique(profiles.e[,c(2)]))
summary(unique(profiles.e[,c(2,20:26)]))
name(profiles.e)
meltp <- melt(unique(profiles.e[, c(2,8:12,32)]),id.vars = c("idp", "hzn_key"))


ggplot(data = meltp[meltp$hzn_key == "Bt" |
                      meltp$hzn_key == "Bt1" |
                      meltp$hzn_key == "Bt2" |
                      meltp$hzn_key == "Bt9" |
                      meltp$hzn_key == "C" |
                      meltp$hzn_key == "Btk" ,],
       aes(x = value, fill = hzn_key)) + geom_density(alpha = 0.2) + 
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

# explore spatial
library(raster)
library(sp)
#devtools::install_github("environmentalinformatics-marburg/mapview", ref = "develop")
library(mapview)

profiles.e <- profiles.e[profiles.e$hzn_disc == "",]
name(profiles.e)


s <- profiles.e[profiles.e$hzn == "C",c(2,8:12,30,31)]

profiles.e <- profiles.e[profiles.e$hzn_disc == "",]
name(profiles.e)
s <- profiles.e[profiles.e$hzn == "A" | 
                  profiles.e$hzn == "B" |
                  profiles.e$hzn == "C",
                c(2,6,8:12,30,31)]
s$hzn_c <- "a"
s$hzn_c[s$hzn=="C"] <- "c"
s[which(s$idp %in% unique(s$idp[s$hzn_c == "c"])),10] <- "c"
name(s)
s <- unique(s[,c(1,8:10)])


coordinates(s) <- ~X+Y
proj4string(s) <- modis
mercator <- CRS("+init=epsg:3785")
s <- spTransform(s, mercator)


mapview(s, burst = TRUE)
mapView(s, map = NULL,
        map.types = mapviewGetOption("basemaps"), zcol = NULL, burst = TRUE,
        color = c("#FF00FF","#0000AF", "#0489B1"), alpha = 0.5,
        col.regions = color, alpha.regions = 0.2,
        na.color = mapviewGetOption("na.color"), at = NULL, cex = 5, lwd = 2,
        popup = popupTable(s), label, legend = mapviewGetOption("legend"),
        legend.opacity = 1, layer.name = deparse(
          substitute(x, env = parent.frame())), 
        verbose = mapviewGetOption("verbose"),
        homebutton = TRUE)
s <- as.data.frame(s)
length(unique(s[s$hzn_c=="c",1:2])[,1])
length(unique(s[s$hzn_c=="a",1:2])[,1])

## AQP FOR PLOT HORIZON FREQUENCY ####
#install.packages('aqp', repos="http://R-Forge.R-project.org")
library(aqp)
names(profiles.e)# <- "name"
s <- profiles.e
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

idp <- unique(us$id.p)

idp.C <- unique(us$id.p[us$H == "C"])

idp.noC <- idp[which(!(idp %in% idp.C))]

us.noC <- layer[which(layer$pedon_key %in% idp.noC), ]









