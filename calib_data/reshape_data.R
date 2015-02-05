rm(list=ls())
#install.packages('reshape')
library(reshape)
library(plyr)
library(lavaan)
library(sp)
library(lattice) # required for trellis.par.set():
trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()
library(corrgram)
#install.packages("lavaan", repos="http://www.da.ugent.be", type="source")

# The first section concerns about grouping horizons with the same clasification, making a 
# weighted mean of soil properties and taking min and max boundaries of horizons; 
# D is the original data calib.data-1.0.csv. It was splited into covD, which are the
# covariates for every location, and s.pr that refers to soil properties of every horizon.
# dat, d2 and dat2 are data.frames intermediates to rich the final dataset: data



#################### RE-SHAPING DATA  ########
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration")

D <-read.csv("calib.data-0.1.1.csv")

#selected soil properties
names(D)
d <- cbind(D[,c(1:4,61,98,95,94,7:16,30:33,18,17,20:26,51,42,55,87,88,90,91,77,81)])
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

# standardising horizon names
t(table(d$horizon))
d$hor <-""
d<- cbind(d[,c(1:6,41,7:40)])
if(d$horizon==)
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
d0 <- d0[!is.na(d0$a_ph_h2o) & !is.na(d0$a_arcilla),]

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
d0$thick <- ((d0$bottom - d0$top)^2)^(1/2)
#(d0$bottom - d0$top)<1
#returns sum of thick by id.hor. This requires installation of the library plyr
thickness <- ddply(d0,.(id.hor), summarise,s.thick=sum(thick))
#merge thickness and s.pr
d1 <- merge(d0, thickness, by = "id.hor", all=T)
#weights
d1$weight <- d1$thick / d1$s.thick

#compute soil properties per id.hor
names(d1)
count(d1$id.hor)
d2 <- cbind(d1[,c(1:10)],(d1[,11:33]* d1[,45]),d1[,c(35:45)])
d2[,26:33][is.na(d2[,26:33])] <- 0
#### merge horizons
## horizon boundaries
limits <- cbind(ddply(d2,.(id.hor), summarise, mintop=min(top))[,1:2],
            maxbot=(ddply(d2,.(id.hor), summarise, maxbot=max(bottom))[,2]))
## soil properties 
names(d2)[11:33]
####  aggregation of horizon by id.hor. Warning! = If one horizon has NA the other horizons result in NA

d3 <-  cbind(ddply(d2,.(id.hor), summarise, a_S=sum(a_sum_bases))[,1:2],
             a_CEC=(ddply(d2,.(id.hor), summarise, a_CEC=sum(a_CEC))[,2]),
             a_base_ca=(ddply(d2,.(id.hor), summarise, a_base_ca=sum(a_base_ca))[,2]),
             a_base_mg=(ddply(d2,.(id.hor), summarise, a_base_mg=sum(a_base_mg))[,2]),
             a_base_k=(ddply(d2,.(id.hor), summarise, a_base_k=sum(a_base_k))[,2]),
             a_base_na=(ddply(d2,.(id.hor), summarise, a_base_na=sum(a_base_na))[,2]),
             a_H=(ddply(d2,.(id.hor), summarise, a_H=sum(a_h))[,2]),
             a_sat_CEC=(ddply(d2,.(id.hor), summarise, a_sat_CEC=sum(a_saturacion_t))[,2]),
             a_sat_SH=(ddply(d2,.(id.hor), summarise, a_sat_SH=sum(a_saturacion_s_h))[,2]),
             a_cond=(ddply(d2,.(id.hor), summarise, a_cond=sum(a_conductividad))[,2]),
             a_res_pasta=(ddply(d2,.(id.hor), summarise, a_res_pasta=sum(a_resistencia_pasta))[,2]),
             a_ph_pasta=(ddply(d2,.(id.hor), summarise, a_ph_pasta=sum(a_ph_pasta))[,2]),
             a_ph_h2o=(ddply(d2,.(id.hor), summarise, a_ph_h2o=sum(a_ph_h2o))[,2]),
             a_ph_kcl=(ddply(d2,.(id.hor), summarise, a_ph_kcl=sum(a_ph_kcl))[,2]),
             a_OC=(ddply(d2,.(id.hor), summarise, a_OC=sum(a_OC))[,2]),
             a_clay=(ddply(d2,.(id.hor), summarise, a_clay=sum(a_arcilla))[,2]),
             a_silt_20=(ddply(d2,.(id.hor), summarise, a_silt_20=sum(a_limo_2_20))[,2]),
             a_silt_50=(ddply(d2,.(id.hor), summarise, a_silt_50=sum(a_limo_2_50))[,2]),
             a_sand_100=(ddply(d2,.(id.hor), summarise, a_sand_100=sum(a_arena_muy_fina))[,2]),
             a_sand_250=(ddply(d2,.(id.hor), summarise, a_sand_250=sum(a_arena_fina))[,2]),
             a_sand_500=(ddply(d2,.(id.hor), summarise, a_sand_500=sum(a_arena_media))[,2]),
             a_sand_1k=(ddply(d2,.(id.hor), summarise, a_sand_1k=sum(a_arena_gruesa))[,2]),
             a_sand_2k=(ddply(d2,.(id.hor), summarise, a_sand_2k=sum(a_arena_muy_gruesa))[,2])
                )

## merge (((limits + d2) + d3) + d2) #Concretions and mottles remain out of this dataset
d4 <- merge(x= merge(x= unique(merge(x = limits,y = d2[,c(1:3,5,6,8)],by = "id.hor",all.x = F, all.y = T)), 
              y= d3, by= "id.hor", all= T),
              y= unique(d2[,c(5,36:41)]), by= "id.p", all.x=T, all.y=F)


# to recover concretions and mottles
d2$moteados[(d2$moteados)==""]<-NA
d2$is.mottles<- as.numeric(!is.na(d2$moteados))
d2$is.mottles[d2$is.mottles==0] <-9999
d2$is.mottles[(d2$is.mottles)<9999]<- d2$top[(d2$is.mottles)<9999]

d2$concreciones[(d2$concreciones)==""]<-NA
d2$is.concr<- as.numeric(!is.na(d2$concreciones))
d2$is.concr[d2$is.concr==0] <-9999
d2$is.concr[d2$is.concr<9999]<- d2$top[d2$is.concr<9999]
# merge d4 + concretions(depth) + mottles(depth)
d5 <-merge(x=merge(d4, ddply(d2,.(id.p), summarise, is.mottles=min(is.mottles)),by= "id.p", all=T),
      y=ddply(d2,.(id.p), summarise, is.concr=min(is.concr)), by= "id.p")

d5$is.concr[d5$is.concr==9999] <-NA
d5$is.mottles[d5$is.mottles==9999] <-NA


### order variables by horizons
A <- d5[d5$hor=="A",c(1:4,9:31)]
B <- d5[d5$hor=="B",c(1:4,9:31)]
E <- d5[d5$hor=="E",c(1:4,9:31)]
C <- d5[d5$hor=="C",c(1:4,9:31)]
names(A)[2:27] <- paste(names(A)[2:27],".A",sep="")
names(B)[2:27] <- paste(names(B)[2:27],".B",sep="")
names(E)[2:27] <- paste(names(E)[2:27],".E",sep="")
names(C)[2:27] <- paste(names(C)[2:27],".E",sep="")

AB <-merge(A,B, by="id.p", all=T)
ABE <- merge(AB,E, by= "id.p", all=T)
ABEC <- merge(ABE, C, by= "id.p", all=T) 

d6 <- merge(unique(d5[,c(1,5:7,32:39)]),ABEC,by="id.p", all=T)
#################### THE END OF RE-SHAPING
rm(list=ls()[ls()!="d6"])
################## ESTIMATING NEW VARIABLES #####

d6$is.Bt <-d6$a_clay.B/d6$a_clay.A
boxplot(d6$is.Bt)
# coordinates(d6) <- ~X+Y
# spplot(d6,zcol ="is.Bt",edge.col="black", colorkey=T,col.regions= rainbow(1000) )


##### CEC analysis

# CEC from OC (Nyle and Weil, 2007: chapter 8)
d6$e_CEC_OM.A <- (d6$a_ph_h2o.A*50/1.5-33.33)*(d6$a_OC.A*1.72/100)
# CEC from clay + silt
d6$e_CEC_cl.sl.A <- d6$a_CEC.A - d6$e_CEC_OM.A 
# CEC from silt considering CEC of silt = 15.5 cmol/kg (range 8 to 23 from Morras, 1995)
d6$e_CEC_sl.A <- 10*(d6$a_silt_20.A/100)
# CEC from clay
d6$e_CEC_cl.A <- d6$e_CEC_cl.sl.A - d6$e_CEC_OM.A - d6$e_CEC_sl.A
# CEC from pure clay 
d6$e_CEC_pure.cl.A <- d6$e_CEC_cl.A/(d6$a_clay.A/100)

plot(d6$a_ph_kcl.A,d6$a_base_na.A)
hist(d6$e_CEC_pure.cl.A,breaks = 20, xlab = "Clay CEC A horizon")


d7 <-d6
coordinates(d7) <- ~X+Y
spplot(d7,zcol ="a_ph_kcl.A",edge.col="black", colorkey=T)

data$clayAB <- data$a_arcilla.A/data$a_arcilla.Bt
data$dist_mean <- data$dist_mean/1000
coordinates(data)=~X+Y
## coloured points plot with legend in plotting area and scales:

names(data)
spplot(data, "a_mo_c.Bt", do.log = TRUE,
       key.space=list(x=0.2,y=0.9,corner=c(0,1)),
       scales=list(draw=T))
data <- as.data.frame(data)
corrgram(data[,-c(22:31)], order=F, lower.panel=panel.shade,
         upper.panel=panel.pie)

data <- data[,c(1:35,38,57,59,60,63,64,66,70,72,75,77,82,84,85)]
data$oc_ba <- data$a_mo_c.Bt/data$a_mo_c.A
data$na_ab <- data$a_base_na.A/data$a_base_na.Bt
data$ca_ab <- data$a_base_ca.A/data$a_base_ca.Bt
data$tick.Bt <- data$mintop.Bt/data$maxbot.Bt
names(data) <- c("id.p","top.A","bot.A","ph.A","ca.A","mg.A","k.A","na.A","om.A","clay.A", "silt.A",
                  "top.Bt","bot.Bt","ph.Bt","ca.Bt","mg.Bt","k.Bt","na.Bt","om.Bt","clay.Bt","silt.Bt",
                        "top.E","bot.E","ph.E","ca.E","mg.E","k.E","na.E","om.E","clay.E",
                 "silt.E","X","Y","DEM_min","DEM_max","DEM_sd","Carea_mean","TWI_min","TWI_max","TWI_sd",
                 "LS_min","LS_sum","ChNBL_max","ChNBL_mean","VDChN_max","VDChN_mean","RSP_mean","dist_mean","clayAB","oc_BA",
                 "na_AB","ca_AB","thick.Bt")
data$ca.mg.Bt <- data$ca.Bt/data$mg.Bt
data$is.E <- as.factor(as.numeric(!is.na(data$top.E)))
corrgram(data, order=F, lower.panel=panel.shade, 
         upper.panel=panel.pie)
boxplot(DEM_min ~ is.E,data)
#write.csv(data,"data.csv")
############### STRUCTURAL EQUATION MODELLING ######################



#data <- read.csv("data.csv")
#Ecov[,"is.E"] <- lapply(Ecov[,"is.E"], ordered)
#names(Ecov)
model_full <- '            ##CREATING MODEL

# latent variables

# path analysis (regressions)
  na.Bt ~ ca.Bt + is.E
  is.E ~ DEM_sd + DEM_min
  ca.Bt ~ is.E + clay.Bt + DEM_sd + LS_sum
  clay.Bt ~ DEM_max + dist_mean + X
  ph.Bt ~ na.Bt + TWI_max + TWI_sd

#  na.A ~ 1*na.Bt + ph.Bt + is.E
#  ph.A ~ 1*na.A + ph.Bt +is.E

# (residual) variances and covariances  
  is.E ~~ is.E
  ca.Bt ~~ ca.Bt
  clay.Bt ~~ clay.Bt
  na.Bt  ~~ na.Bt
  ph.Bt ~~ ph.Bt

#  LS_sum ~~ LS_sum
#   silt.A ~~ silt.A
#  DEM_max ~~ DEM_max
# 
#   TWI_max ~~ TWI_max
#   TWI_sd ~~ TWI_sd
#   DEM_sd ~~ DEM_sd
#   dist_mean ~~ dist_mean
#   silt.Bt ~~ silt.Bt
#   ph.A ~~ ph.A
#   X ~~ X
#   na.Bt ~~ ca.Bt + is.E
#   ca.Bt ~~ is.E + clay.Bt + DEM_sd + LS_sum
#   clay.Bt ~~ silt.A + DEM_max + dist_mean + X
#   ph.Bt ~~ na.Bt + TWI_max + TWI_sd
#   na.A ~~ na.Bt + ph.Bt + is.E
#   ph.A ~~ na.A + ph.Bt +is.E  
#   TWI_max ~~ TWI_sd
'

#fitting
#fit <- lavaan(model_ca, data=data, fixed.x = T, auto.var=T, auto.fix.single=F, std.lv = T )
fit <- lavaan(model_full, data=data)

summary(fit, standardized=T, modindices = F, fit.measures=T, rsquare=T)   #summary of the fitting the model
parameterEstimates(fit, ci = FALSE, standardized = TRUE)
inspect(fit,"cov.ov")      #shows the covariance matrix with negative values
fitted(fit)
resid(fit)
lavPredict(fit,newdata = data, type = "yhat") 

piet <- lav_partable_independence(fit)
lav_partable_npar(piet)
lav_partable_ndat(piet)
lav_partable_df(piet)

##########################################plotting######################################
#install.packages("igraph")
library("igraph")
#PLOT THE MODEL:
##variables of model
attrb<-unique(c(fit@ParTable$lhs,fit@ParTable$rhs))

##model connections
el<-rbind(match(fit@ParTable$rhs,attrb),match(fit@ParTable$lhs,attrb),round(fit@Fit@est,2))
el<-el[,!(fit@ParTable$op=="~~")]
g2 <- add.edges(graph.empty(length(attrb)), el[1:2,], weight=el[3,])
g2<-set.vertex.attribute(g2,"label", index=1:length(attrb), as.vector(attrb))

## color of the variables:
colv<-rep("orange",length(attrb))
colv[which(is.na(match(1:7,el[2,])))]<-"yellow"
colv[which(is.na(match(1:7,el[1,])))]<-"brown"
### independent variables are yellow,
### dependent variables that also are used as predictors are orange,
### variables that are only predicted are brown. 

## color of the connections: 
cole<-get.edge.attribute(g2,"weight")
cole[cole<0]<-(-4)
cole[cole>0]<-2
### red is positive correlation, 
### blue is negative correlation

## reorder the positions of the variables by hand
plotnr<-tkplot(g2,vertex.color=colv,edge.color=abs(cole),vertex.label=get.vertex.attribute(g2,"label"),vertex.label=get.vertex.attribute(g2,"label"),edge.label=get.edge.attribute(g2,"weight"))
coordinates<-tkplot.getcoords(plotnr)
plot(g2,vertex.color=colv,edge.color=abs(cole),vertex.label=get.vertex.attribute(g2,"label"),vertex.label=get.vertex.attribute(g2,"label"),edge.label=get.edge.attribute(g2,"weight"),layout=coordinates)
### This graph will also open in a separate window. Here you can 'drag and drop' the variables to optimize the layout of the model.
