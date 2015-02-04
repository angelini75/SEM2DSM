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

############### PREPROCESSING DATA ###########################

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
d2 <- cbind(d1[,c(1:10)],(d1[,11:33]* d1[,45]),c(35:45))

#### merge horizons
## horizon boundaries
limits <- cbind(ddply(d2,.(id.hor), summarise, mintop=min(top))[,1:2],
            maxbot=(ddply(d2,.(id.hor), summarise, maxbot=max(bottom))[,2]))
## soil properties 
names(d2)
 d3 <-  cbind(ddply(d2,.(id.hor), summarise, a_ph_h2o=sum(a_ph_h2o))[,1:2],
                a_base_ca=(ddply(d2,.(id.hor), summarise, a_base_ca=sum(a_base_ca))[,2]),
                a_base_mg=(ddply(d2,.(id.hor), summarise, a_base_mg=sum(a_base_mg))[,2]),
                a_base_k =(ddply(d2,.(id.hor), summarise,   a_base_k=sum(a_base_k))[,2]),
                a_base_na=(ddply(d2,.(id.hor), summarise, a_base_na=sum(a_base_na))[,2]),
                a_mo_c=   (ddply(d2,.(id.hor),    summarise,    a_mo_c=sum(a_mo_c))[,2]),
                a_arcilla=(ddply(d2,.(id.hor), summarise, a_arcilla=sum(a_arcilla))[,2]),
                a_limo_2_5=(ddply(d2,.(id.hor), summarise, a_limo_2_5=sum(a_limo_2_5))[,2])
                )

## merge (((limits + d2) + dat2) + covD)
data <- merge(x= merge(x= unique(merge(x = limits,y = d2[,1:3],by = "id.hor",all.x = F, all.y = T)), 
              y= dat2, by= "id.hor", all= T),
              y= covD, by= "id.p", all.x=T, all.y=F)
# data$E <- as.numeric(data$hor == "E")
# E <-ddply(data,.(id.p), summarise,is.E=sum(E))
# 
# Ecov<- merge(E,covD, by="id.p", all.x=T)
# Ecov$is.E <- as.factor(Ecov$is.E)
# Ecov<- Ecov[order(Ecov$is.E),]


A <- data[data$hor=="A",c(1,3,4,6:13)]
Bt <- data[data$hor=="Bt",c(1,3,4,6:13)]
E <- data[data$hor=="E",c(1,3,4,6:13)]
names(E) <- c("id.p","mintop.E","maxbot.E","a_ph_h2o.E","a_base_ca.E","a_base_mg.E","a_base_k.E","a_base_na.E","a_mo_c.E","a_arcilla.E","a_limo_2_5.E")
ABt <-merge(A,Bt, by="id.p", all=T, suffixes = c(".A",".Bt"))
ABtE <- merge(ABt,E, by= "id.p", all=T)
data <- merge(ABtE, covD, by= "id.p", all.x=T) 
#names(dat3) <- tolower(names(data))
# dat3$id.hor <- as.factor(dat3$id.hor)
# dat4 <- melt(dat3, id=c("id.p", "id.hor"), na.rm=F)
# data <- cast(dat4, id.p ~ hor ~ variable)
# data[[2]]
# table(dat3)
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
