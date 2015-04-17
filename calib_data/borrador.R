
# it comes from reshape file



##### CEC analysis

# # CEC from OC (Nyle and Weil, 2007: chapter 8)
# d6$e_CEC_OM.A <- (d6$a_ph_h2o.A*50/1.5-33.33)*(d6$a_OC.A*1.72/100)
# # CEC from clay + silt
# d6$e_CEC_cl.sl.A <- d6$a_CEC.A - d6$e_CEC_OM.A 
# # CEC from silt considering CEC of silt = 15.5 cmol/kg (range 8 to 23 from Morras, 1995)
# d6$e_CEC_sl.A <- 10*(d6$a_silt_20.A/100)
# # CEC from clay
# d6$e_CEC_cl.A <- d6$e_CEC_cl.sl.A - d6$e_CEC_OM.A - d6$e_CEC_sl.A
# # CEC from pure clay 
# d6$e_CEC_pure.cl.A <- d6$e_CEC_cl.A/(d6$a_clay.A/100)
######################################## PLOTS
write.csv(d6, "calib.data-1.2.csv")
##
d6 <-read.csv("calib.data-1.2.csv")[,-1]
names(d6)
# Titles and font size
boxplot(d6[,c(31,33,83,85,57,59,109,111)],main="Percentage of clay and silt by horizon",ylab="percentage", xlab="Particle size by horizon",
        cex.main=1.5,cex.lab=1.3,col=c(rep("brown",2),rep("violet",2),rep("orange",2),rep("yellow",2)))


d6$a_base_k.BC[d6$a_base_k.BC>10 & !is.na(d6$a_base_k.BC)]<-1.852

# soiltexture::TT.plot()
# 
# d6 <- cbind(d6,(d6[34]+d6[35]+d6[36]+d6[37]+d6[38]))
# texture.A<- d6[,c(31,33,123)]
# names(texture.A) <- c("CLAY","SILT","SAND")
# texture.A <- na.omit(texture.A)
# texture.A <- texture.A[c(-162,-108),]
# texture.A <- cbind(texture.A,"sum"=(texture.A[1]+texture.A[2]+texture.A[3]))
# names(texture.A)[4] <- "sum"
# texture.A$CLAY <- texture.A$CLAY/texture.A$sum*100
# texture.A$SILT <- texture.A$SILT/texture.A$sum*100
# texture.A$SAND <- texture.A$SAND/texture.A$sum*100
# TT.plot(tri.data =texture.A,  class.sys = "USDA.TT", lang = "en",
#         col = "blue", cex = 0.5)  # English, default
# d6$thickness.A <- d6$maxbot.A-d6$mintop.A
# d6$thickness.B <- d6$maxbot.B-d6$mintop.B
# d6$thickness.E <- d6$maxbot.E-d6$mintop.E
# d6$is.E <- as.factor(!is.na(d6$thickness.E>0))

plot(d6$a_base_na.BC,d6$a_silt_50.BC)
#+d6$a_base_k.B )

summary(lm(a_clay.B~a_CEC.B, d6))



#


















hist(d6$e_CEC_pure.cl.A,breaks = 20, xlab = "Clay CEC A horizon")


d7 <-d6

coordinates(d7) <- ~X+Y
spplot(d7,zcol ="a_caco3.BC",edge.col="black", colorkey=T)

#install.packages("ggmap")
library(ggmap)
ggmap(get_googlemap(c(-60,-35)))

ext <-  get_map(c(min(d6$X),min(d6$Y),-59,max(d6$Y)), color = "bw")
ext <- get_googlemap(center = c((-59+min(d6$X))/2, (max(d6$Y)+min(d6$Y))/2),color="bw",
                     zoom=8, maptype =  "terrain", scale=2)


ggmap(ext,) +
  geom_point(aes(x = X, y = Y, colour = a_silt_50.BC), data = d6, alpha = 1)



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
corrgram(d6[,120:127], order=F, lower.panel=panel.shade,   upper.panel=panel.pie)
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