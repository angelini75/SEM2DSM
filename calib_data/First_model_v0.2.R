rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

d <-read.csv("calib.data-2.1.csv")[,-1]

############### PRE-PROCESS ################## 
names(d)
names(d)[c(5,6,9,10)]<- c("tb.A","sat.A", "oc.A","bt")
d<-d[!is.na(d$oc.A),]
d$sat.A[is.na(d$sat.A)] <- 100 # it is assumed 100% saturation when CaCO3 is present 
d$d.caco3[is.na(d$d.caco3)] <- 300 # it is assumed that CaCO3 is very deep when it is absent within the solum 
summary(d$tb.A[d$is.caco3==1], omit.na=T)
d$tb.A[is.na(d$tb.A)] <- 20

# plot(d$esp.A[d$is.caco3==1]~d$d.caco3[d$is.caco3==1], omit.na=T)
# boxplot(d$d.caco3[is.na(d$esp.A)])
summary(d$esp.A[d$d.caco3<15], omit.na=T)
d$esp.A[is.na(d$esp.A)] <- 16.3

# plot(d$esp.B[d$is.caco3==1]~d$esp.A[d$is.caco3==1])
summary(d$d.caco3[is.na(d$esp.B)])
summary(d$esp.B[d$d.caco3<30], omit.na=T)
d$esp.B[is.na(d$esp.B)] <- 30.9

nas<-d[!complete.cases(d),]
#Ecov[,"is.E"] <- lapply(Ecov[,"is.E"], ordered)
#names(Ecov)

# normality test
for(i in 1:length(d)){
  print(names(d)[i])
  print(shapiro.test(d[,i]))
}
paste(names(d),".r", " =~ ",names(d), sep="")
# transformation
d$esp.A <- log10(d$esp.A)
d$esp.B <- log10(d$esp.B)
d$is.hydro <- ordered(d$is.hydro)
d$is.E <- ordered(d$is.E)
d$is.caco3 <- ordered(d$is.caco3)
#d[,c(4:10,12,15:27)]<- scale(d[,c(4:10,12,15:27)],)
# d$river <- d$river/100000
# d$evim <- d$evim/1000
# d$evisd <- d$evisd/1000
# d$wdist <- d$wdist/1000
# d$lstm <- (d$lstm-290)/10
# d$lstsd <- d$lstsd-6
# d$vdchn <- d$vdchn/10
# d$mrvbf <- d$mr/10
# d$slope <- d$slope*10
# d$thick.A <- d$thick.A/100
# d$tb.A <- d$tb.A/100
# d$sat.A <- d$sat.A/100
# d$slope <- d$slope/10
# d$twi <- d$twi/10
# d$dem <- d$dem/100
# d$maxc <- d$maxc/100

e<- (((4-mean(d$thick.A))/sd(d$thick.A))^2)^(1/2)
##@## data normalization
n <- c(4:10,15:27)
for(i in n){
  d[,i]<- (d[,i]-mean(d[,i]))/sd(d[,i])
}

boxplot(d[,c(4:10,15:27)])  


# step(lm(sat.A ~ dem+wdist+maxc+mrvbf+slope+twi+vdchn+lstm+lstsd+evim+evisd+river, d),direction = "both")
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd, data = d))
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd + tb.A + d.caco3 + river, data = d))
summary(lm(formula = oc.A ~ mrvbf + lstm + lstsd + evisd, data = d))

summary(lm(formula = tb.A ~ wdist + lstm + evim + evisd + river, data = d))
# summary(lm(formula = thick.A ~ dem + maxc + evisd, data = d))
# summary(lm(formula = esp.B ~ mrvbf + twi + vdchn + lstsd + evim + evisd + 
#              river, data = d))
# summary(lm(formula = esp.A ~ dem + mrvbf + twi + vdchn + lstsd + evim + 
#              esp.B, data = d))
# summary(lm(formula = sat.A ~ maxc + lstm + lstsd + evim + evisd, data = d))
# summary(lm(formula = d.caco3 ~ dem + mrvbf + vdchn + lstm + lstsd + evim + 
#             evisd + river, data = d))


############### FITTING MODEL ######################

#### Third Model####
third_model <- '
# measurement model
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
btr =~ 1*bt
oc.Ar =~ 1*oc.A
thick.Ar =~ 1*thick.A
esp.Br =~ 1*esp.B
esp.Ar =~ 1*esp.A

# structural model
tb.Ar ~     evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            oc.Ar 
sat.Ar ~    evim + evisd + lstm + lstsd + dem + wdist + mrvbf + vdchn +  twi + river + 
            tb.Ar + oc.Ar                            
btr ~       lstm +  lstsd + wdist + vdchn + twi + dem + river + mrvbf +
            esp.Br + esp.Ar
oc.Ar ~     lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi +
            esp.Br + esp.Ar + btr + thick.Ar
esp.Ar ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river + 
            esp.Br 
esp.Br ~    lstm +  lstsd + dem + wdist + mrvbf + vdchn + twi + river 
            
thick.Ar ~  dem + wdist + mrvbf + vdchn + twi + river + slope + maxc + evim + evisd 


# residual (co)variance
thick.A ~~  0.25*thick.A
tb.A ~~     0.25*tb.A




# intercepts
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))

fit3 <- lavaan(model = third_model, data = d,
               model.type = "sem", meanstructure = "default",
               int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default",
               orthogonal = FALSE, std.lv = FALSE, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default",
               constraints = "", estimator = "ML",
               zero.cell.warn = TRUE, start = "default")
varTable(fit3)
summary(fit3, standardized=F, modindices = F, fit.measures=F) 
inspect(fit3,"std.lv") # standardized model parameters
inspect(fit3,"sampstat") #Observed sample statistics
inspect(fit3, "cov.lv") #The model-implied covariance matrix of the observed variables
inspect(fit3,"start") # starting values for all model parameters
inspect(fit3,"rsquare") #R-square value for all endogenous variables
full_model <- lavaanify(third_model)
modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>3 & !is.na(modi$mi),]
modi

write.csv(full_model, "full_model.csv")

#################### PREDICTION #########################

#################Prediction 1A##############
##setting up matrices
co <- as.data.frame(coef(fit3))
co$fun <- row.names(co)
row.names(co)<-NULL
{B <- inspect(fit3,"est")$beta[1:7,1:7]} #matrix of coeff. latent state variables
 
# B <-matrix(c(0,coef(fit3)["silt_1A_r~thick_1A_r"],coef(fit3)["OM_1A_r~thick_1A_r"],  
#               coef(fit3)["thick_1A_r~silt_1A_r"],0,coef(fit3)["OM_1A_r~silt_1A_r"],
#               coef(fit3)["thick_1A_r~OM_1A_r"],coef(fit3)["silt_1A_r~OM_1A_r"],0),
#             ncol=3,nrow=3)
#  B[is.na(B)] <- 0} #Matrix endogenous variables B

{I<-diag(nrow=7, ncol=7)} #Identity matrix

{A<- inspect(fit3,"est")$beta[1:7,8:19]} #matrix of coeff of external drivers
# {A<-matrix(c(coef(fit3)["thick_1A_r~age_hy"],coef(fit3)["silt_1A_r~age_hy"],         
#              coef(fit3)["OM_1A_r~age_hy"],coef(fit3)["thick_1A_r~VI"],
#              coef(fit3)["silt_r~VI"],coef(fit3)["OM_1A_r~VI"]),
#            nrow=3,ncol=2)
#  A[is.na(A)] <- 0
# } #Matrix exogenous variables


##### ...to be continued.



#Create test rasters
age <- as.vector(seq(from=0,to=140,by=140/49))
VI <- seq(from=100,to=40,by=-(100-40)/49)

age_df=matrix(nrow=50,ncol=50)
for (i in 1:50){
  age_df[i,]<-age}
VI_df=matrix(nrow=50,ncol=50)
for (i in 1:50){
  VI_df[,i]<-VI}

#################Prediction model###########################
#For test data
thick_pred<-matrix(nrow=50,ncol=50)
silt_pred<-matrix(nrow=50,ncol=50)
OM_pred<-matrix(nrow=50,ncol=50)

k=matrix(c(coef(fit3)["thick_1A_r~1"],coef(fit3)["silt_1A_r~1"],coef(fit3)["OM_1A_r~1"]),nrow=3,ncol=1)
for (i in 1:50){
  for (j in 1:50){
    p=matrix(c(age_df[i,j],VI_df[i,j]),nrow=2,ncol=1)
    n=c(0,0,0)
    n=(solve(I-B))%*%((A%*%p)+k)
    thick_pred[i,j]<-n[1]
    silt_pred[i,j] <-n[2]
    OM_pred[i,j] <-n[3]
  }
}

thick_raster<-raster(x=thick_pred)  
OM_raster<-raster(x=OM_pred)
silt_raster<-raster(silt_pred)

##On study area
thick_1A<-age_valley
silt_1A<-age_valley
OM_1A<-age_valley

for (i in 1:nrow(age_valley)){
  for (j in 1:ncol(age_valley)){
    p=matrix(c(age_valley[i,j],VI_valley[i,j]),nrow=2,ncol=1)
    n=c(0,0,0)
    k=matrix(c(coef(fit3)["thick_1A_r~1"],coef(fit3)["silt_1A_r~1"],coef(fit3)["OM_1A_r~1"]),nrow=3,ncol=1)
    n=(solve(I-B))%*%(A%*%p+k)
    thick_1A[i,j]<-n[1]
    silt_1A[i,j] <-n[2]
    OM_1A[i,j] <-n[3]
  }
}


#################Plotting results####################
##Study area
terr<-readShapePoly("E:/thesis/Spitsbergen 2014/data veldwerk/gis-bestanden/study_area_individual terraces_estimate.shp")

bmp(filename="E:/thesis/statistics/thick_1A.bmp",width=450,height=450)
plot(thick_1A,main="Thickness (cm)")
plot(terr,add=TRUE)
dev.off()

bmp(filename="E:/thesis/statistics/silt_1A.bmp",width=450,height=450)
plot(silt_1A,main="Silt (%)")
plot(terr,add=TRUE)
dev.off()

silt_1A_total<-silt_1A*thick_1A/100*1651/100
bmp(filename="E:/thesis/statistics/silt_1A_total.bmp",width=450,height=450)
plot(silt_1A_total,main="Silt total (kg)")
plot(terr,add=TRUE)
dev.off()

bmp(filename="E:/thesis/statistics/OM_1A.bmp",width=450,height=450)
plot(OM_1A,main="Organic matter (%)")
plot(terr,add=TRUE)
dev.off()

OM_1A_total<-thick_1A/100*(OM_1A/100)*1651
bmp(filename="E:/thesis/statistics/OM_1A_total.bmp",width=450,height=450)
plot(OM_1A_total,main="Organic matter (kg)")
plot(terr,add=TRUE)
dev.off()

bmp(filename="direction where the file has to go",width=450,height=450)
plot(plot van de grafiek)
plot(eventuele dingen bijplotten)
dev.off()




##test rasters
##test data
library(ggplot2)
library(raster)
library(rasterVis)
install.packages("lattice", dependencies = TRUE) 
library(lattice)

V_data$x<-V_data$age_hy/140
V_data$y<-(V_data$VI-40)/(100-40)
V<-subset(V_data,OM_1A>0,select=c(samp_nr))
plot_points<-subset(V_data,OM_1A>0,select=c(x,y))
plot_points<-as.matrix(plot_points)
plot_pts<-SpatialPoints(plot_points,proj4string=CRS(as.character(NA)))
plots<-SpatialPointsDataFrame(plot_points,V, proj4string=CRS(as.character(NA)))
xy.coords(plot_points)
pts <- sampleRandom(thick_raster, size=20, sp=TRUE)

#thickness
levelplot(thick_raster, xlab="age",ylab="vegetation index", main=paste("Thickness (cm)"),
          col.regions = rev(terrain.colors(255)),
          scales=list(x=list(at=seq(0,1,by=0.25),labels=c(0,3500,7000,10500,14000)),
                      y=list(at=c(0,1/6,2/6,3/6,4/6,5/6,6/6),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))),
          cuts=254, margin=FALSE)+layer(sp.points(plot_points,pch=1,cex=2,col="black"))

#silt
levelplot(silt_raster, xlab="age",ylab="vegetation index", main=paste("Silt (%)"),
          col.regions = rev(terrain.colors(255)),
          scales=list(x=list(at=seq(0,1,by=0.25),labels=c(0,3500,7000,10500,14000)),
                      y=list(at=c(0,1/6,2/6,3/6,4/6,5/6,6/6),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))),
          cuts=254, margin=FALSE)+layer(sp.points(plot_points,pch=1,cex=2,col="black"))

#Organic matter
levelplot(OM_raster, xlab="age",ylab="vegetation index", main=paste("Organic matter (%)"),
          col.regions = rev(terrain.colors(255)),
          scales=list(x=list(at=seq(0,1,by=0.25),labels=c(0,3500,7000,10500,14000)),
                      y=list(at=c(0,1/6,2/6,3/6,4/6,5/6,6/6),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))),
          cuts=254, margin=FALSE)+layer(sp.points(plot_points,pch=1,cex=2,col="black"))






########Plot SEM chart (made by Ype)########
library("igraph")
#PLOT THE MODEL:
##variables of model
attrb<-unique(c(fit3@ParTable$lhs,fit3@ParTable$rhs))

##model connections
el<-rbind(match(fit3@ParTable$rhs,attrb),match(fit3@ParTable$lhs,attrb),round(fit3@Fit@est,2))
el<-el[,!(fit3@ParTable$op=="~~")]
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
