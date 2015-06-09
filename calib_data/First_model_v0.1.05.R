rm(list=ls())
#install.packages("lavaan")
#install.packages("lavaan", repos="http://www.da.ugent.be", type="source") 
library(lavaan)
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

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

tb.Ar ~      evim + evisd + river + lstm +  dem + wdist +          #mrvbf + vdchn +  twi +
              is.caco3 + oc.Ar + esp.Br + esp.Ar + btr
sat.Ar ~      lstm + lstsd + wdist +   mrvbf +  vdchn +            #evim + evisd +twi + dem +
              tb.Ar + is.caco3 + is.hydro + esp.Ar + oc.Ar 
btr ~        lstm +  wdist + vdchn +   twi + dem + river +         # mrvbf + lstsd +
              esp.Br + esp.Ar + is.hydro + is.caco3
is.caco3 ~  wdist + vdchn + lstm + lstsd + evim +                 #+ twi + mrvbf + dem 
              esp.Br 
is.hydro ~  wdist + vdchn + mrvbf + evisd + lstm +                #twi +    dem +
              esp.Br + esp.Ar + is.E + is.caco3
is.E ~      wdist + vdchn +  dem +                                #mrvbf + twi + 
              is.hydro + esp.Ar + esp.Br
oc.Ar ~      wdist +  mrvbf + twi +  evisd + lstm + lstsd + dem +  #evim + vdchn + 
              esp.Br + esp.Ar #+ is.hydro                           #thick.A +
thick.Ar ~   wdist + maxc + dem + vdchn +                          # slope + mrvbf + lstm + lstsd + river + twi + 
              esp.Br + esp.Ar + is.caco3                            # + is.hydro
esp.Ar ~     lstsd +  lstm + river + evim +                        #evisd + wdist + vdchn +  mrvbf + twi +  dem +
              esp.Br 
esp.Br ~     lstm + lstsd + vdchn +  mrvbf + twi + river + evisd + #  evim + wdist +  dem 
             is.caco3 

is.caco3 ~~ esp.B
thick.A ~~  thick.A
sat.A ~~       bt


# tb.A ~~ tb.A
# sat.A ~~ sat.A
# bt ~~ bt
# is.caco3 ~~ is.caco3
# is.hydro ~~ is.hydro
# is.E ~~ is.E
# oc.A ~~ oc.A
# esp.A ~~ esp.A
# esp.B ~~ esp.B

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
fit2 <- lavaan(model = third_model, data = d,
               model.type = "sem", meanstructure = T, #meanstructure requesting intercepts
               int.ov.free = F, int.lv.free = F, fixed.x = "default", #intercepts are estimated instead 0
               orthogonal = FALSE, std.lv = T, 
               parameterization = "default", auto.fix.first = F, #WHEN LATENT VARIABLES R PRESENT
               auto.fix.single = F, auto.var = T, auto.cov.lv.x = F, # auto.var If TRUE, the residual variances and the variances of exogenous latent variables are included in the model and set free.
               auto.cov.y = T, auto.th = T, auto.delta = T, # ?? auto.th if T they are free, otherwise calculated from data ??
               std.ov = FALSE, missing = "default", ordered = c("is.E","is.hydro", "is.caco3"),
               sample.cov = NULL, sample.cov.rescale = "default",
               sample.mean = NULL, sample.nobs = NULL, ridge = 1e-05,
               group = NULL, group.label = NULL, group.equal = "", group.partial = "",
               group.w.free = FALSE, cluster = NULL, 
               constraints = "", estimator = "WLSMV",
               likelihood = "default", link = "default", information = "default",
               se = "default", test = "default", bootstrap = 1000L, mimic = "default",
               representation = "default", do.fit = TRUE, control = list(),
               WLS.V = NULL, NACOV = NULL, 
               zero.add = "default", zero.keep.margins = "default", 
               zero.cell.warn = TRUE, start = "default", 
               slotOptions = NULL, slotParTable = NULL,
               slotSampleStats = NULL, slotData = NULL, slotModel = NULL,
               slotCache = NULL, verbose = FALSE, warn = TRUE, debug = FALSE)

summary(fit2, standardized=F, modindices = F, fit.measures=F) 
inspect(fit2,"std.lv") # standardized model parameters
inspect(fit2,"sampstat") #Observed sample statistics
inspect(fit2, "cor.lv") #The model-implied correlation matrix of the observed variables
inspect(fit2,"start") # starting values for all model parameters
inspect(fit2,"information") #variance covariance matrix of the standardized estimated model parameters
inspect(fit2,"rsquare") #R-square value for all endogenous variables

modi<-summary(fit2, standardized=F, modindices = T, fit.measures=F) 
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