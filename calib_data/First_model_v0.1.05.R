rm(list=ls())
# install.packages("lavaan")
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
setwd("D:/UserData/BaseARG/2_Calibration/simplest_model")

d <-read.csv("calib.data-2.0.csv")[,-1]

#### preprocess 
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
############### STRUCTURAL EQUATION MODELLING ######################

#### Third Model####
third_model <- '
tb.Ar =~ 1*tb.A
sat.Ar =~ 1*sat.A
btr =~ 1*bt
# is.caco3r =~ 1*is.caco3
# is.hydror =~ 1*is.hydro
# is.Er =~ 1*is.E
oc.Ar =~ 1*oc.A
thick.Ar =~ 1*thick.A
esp.Br =~ 1*esp.B
esp.Ar =~ 1*esp.A


tb.Ar ~     evim + evisd + river + lstm +  dem + wdist +          #mrvbf + vdchn +  twi +
            oc.Ar + esp.Br + esp.Ar + btr                         #is.caco3r +
sat.Ar ~    lstm + lstsd + wdist +   mrvbf +  vdchn +            #evim + evisd +twi + dem +
            tb.Ar + esp.Ar + oc.Ar                                #is.caco3r + is.hydror +
btr ~       lstm +  wdist + vdchn +   twi + dem + river +         # mrvbf + lstsd +
            esp.Br + esp.Ar                                       #+ is.hydror + is.caco3r
# is.caco3r ~  wdist + vdchn + lstm + lstsd + evim +                 #+ twi + mrvbf + dem 
#             esp.Br 
# is.hydror ~  wdist + vdchn + mrvbf + evisd + lstm +                #twi +    dem +
#             esp.Br + esp.Ar                                     # + is.Er + is.caco3r
# is.Er ~      wdist + vdchn +  dem +                                #mrvbf + twi + 
#             is.hydror + esp.Ar + esp.Br
oc.Ar ~     wdist +  mrvbf + twi +  evisd + lstm + lstsd + dem +  #evim + vdchn + 
            esp.Br + esp.Ar                                     #+ is.hydror                           #thick.A +
thick.Ar ~  wdist + maxc + dem + vdchn +                          # slope + mrvbf + lstm + lstsd + river + twi + 
            esp.Br + esp.Ar                                     #+ is.caco3r                            # + is.hydro
esp.Ar ~    lstsd +  lstm + river + evim +                        #evisd + wdist + vdchn +  mrvbf + twi +  dem +
            esp.Br 
esp.Br ~    lstm + lstsd + vdchn +  mrvbf + twi + river + evisd     # + #  evim + wdist +  dem 
                                                                      #is.caco3r 
#Measurement error
thick.A ~~  0.25*thick.A
tb.A ~~     0.25*tb.A

# tb.A ~1
# sat.A ~1
# bt ~1
# is.caco3 ~1
# is.hydro ~1
# is.E ~1
# oc.A ~1
# thick.A ~1
# esp.A ~1
# esp.B ~1
'
##### fitting ####
#fit3<- sem(third_model, d,meanstructure = T,std.lv = T, ordered = c("is.E","is.caco3","is.hydro"))

fit3 <- lavaan(model = third_model, data = d,
               model.type = "sem", meanstructure = "default",
               int.ov.free = F, int.lv.free = F, fixed.x = "default",
               orthogonal = FALSE, std.lv = T, 
               parameterization = "default", auto.fix.first = F,
               auto.fix.single = T, auto.var = T, auto.cov.lv.x = FALSE,
               auto.cov.y = T, auto.th = T, auto.delta = FALSE,
               std.ov = FALSE, missing = "default", #ordered = c("is.hydro","is.caco3", "is.E"),
               constraints = "", estimator = "ML", #  "WLSMV",
               zero.cell.warn = TRUE, start = "default")
varTable(fit3)
summary(fit3, standardized=F, modindices = F, fit.measures=F) 
inspect(fit3,"std.lv") # standardized model parameters
inspect(fit3,"sampstat") #Observed sample statistics
inspect(fit3, "cov.lv") #The model-implied covariance matrix of the observed variables
inspect(fit3,"start") # starting values for all model parameters
inspect(fit3,"rsquare") #R-square value for all endogenous variables
lavaanify(third_model)
modi<-summary(fit3, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>3 & !is.na(modi$mi),]
modi








# std.ov = If TRUE, observed variables are standardized before entering the analysis
# std.lv = If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0.
fit1 <- sem(third_model, data=d, estimator ="WLSMV", meanstructure = T, 
            ordered = c("is.E","is.hydro", "is.caco3"))#, fixed.x=F)
varTable(fit1)
# inspect(fit1,"std.lv") # standardized model parameters
# inspect(fit1,"sampstat") #Observed sample statistics
inspect(fit1, "cor.ov") #The model-implied correlation matrix of the observed variables
# inspect(fit1,"start") # starting values for all model parameters
# inspect(fit1,"vcov.std.all") #variance covariance matrix of the standardized estimated model parameters
inspect(fit1,"rsquare") #R-square value for all endogenous variables
summary(fit1, standardized=F, modindices = F, fit.measures=T) 
modi<-summary(fit1, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>3 & !is.na(modi$mi),]
modi