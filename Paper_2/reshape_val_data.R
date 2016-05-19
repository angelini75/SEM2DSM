rm(list = ls())
name <- function(x) { as.data.frame(names(x))} 
setwd("/home/mangelini/big/SEM_2nd_paper/")

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

# Var is the error variance from SEM and Ns are the mean and standard deviation used to standardise the variables
Var <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/SEM.variance.error.csv")
Ns <- read.csv("/media/marcos/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model/N.march2016.csv")[,-1]
names(Var) <- c("property", "variance")