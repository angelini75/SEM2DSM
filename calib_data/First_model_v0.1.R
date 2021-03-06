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

plot(d$esp.A[d$is.caco3==1]~d$d.caco3[d$is.caco3==1], omit.na=T)
boxplot(d$d.caco3[is.na(d$esp.A)])
summary(d$esp.A[d$d.caco3<15], omit.na=T)
d$esp.A[is.na(d$esp.A)] <- 16.3

plot(d$esp.B[d$is.caco3==1]~d$esp.A[d$is.caco3==1])
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
# d$is.hydro <- ordered(d$is.hydro)
# d$is.E <- ordered(d$is.E)


# step(lm(sat.A ~ dem+wdist+maxc+mrvbf+slope+twi+vdchn+lstm+lstsd+evim+evisd+river, d),direction = "both")
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd, data = d))
# summary(lm(formula = bt ~ dem + maxc + slope + lstm + evisd + tb.A + d.caco3 + river, data = d))
# summary(lm(formula = oc.A ~ mrvbf + lstm + lstsd + evisd, data = d))
# summary(lm(formula = tb.A ~ wdist + lstm + evim + evisd + river, data = d))
# summary(lm(formula = thick.A ~ dem + maxc + evisd, data = d))
# summary(lm(formula = esp.B ~ mrvbf + twi + vdchn + lstsd + evim + evisd + 
#              river, data = d))
# summary(lm(formula = esp.A ~ dem + mrvbf + twi + vdchn + lstsd + evim + 
#              esp.B, data = d))
# summary(lm(formula = sat.A ~ maxc + lstm + lstsd + evim + evisd, data = d))
# summary(lm(formula = d.caco3 ~ dem + mrvbf + vdchn + lstm + lstsd + evim + 
#             evisd + river, data = d))
############### STRUCTURAL EQUATION MODELLING ######################

first_model <- '            ##CREATING MODEL
# measurement model           ## Indicate in what manner the latent variables depend on the measured variables 
                              ## (arrows in conceptual model)
                              ## multiplication by 1 to set factor loading to 1 
thick.A.r =~ thick.A 
tb.A.r =~ 1*tb.A 
sat.A.r =~ sat.A
esp.A.r =~ 1*esp.A 
esp.B.r =~ 1*esp.B 
oc.A.r =~ 1*oc.A 
bt.r =~ 1*bt 
d.caco3.r =~ d.caco3 
is.E.r =~ 1*is.E 
is.hydro.r =~ 1*is.hydro 
dem.r =~ 1*dem 
river.r =~ 1*river 
wdist.r =~ 1*wdist 
#water.r =~ 1*water 
maxc.r =~ 1*maxc 
mrvbf.r =~ 1*mrvbf 
slope.r =~ 1*slope 
twi.r =~ 1*twi 
vdchn.r =~ 1*vdchn 
lstm.r =~ 1*lstm 
lstsd.r =~ 1*lstsd 
evim.r =~ 1*evim 
evisd.r =~ 1*evisd 

# path analysis (regressions)           #indicates the regressions between the different latent variables
tb.A.r ~      dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r+
              oc.A.r + d.caco3.r
sat.A.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              oc.A.r +  d.caco3.r 
bt.r ~        dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + river.r +
              esp.A.r + d.caco3.r #+ esp.B.r + 
esp.A.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              esp.B.r
esp.B.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r 
            # esp.A.r
oc.A.r ~      dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              thick.A.r
d.caco3.r ~   dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r +
              esp.B + bt.r
is.hydro.r ~  dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r +
              esp.A.r # it does not converge with esp.B.r but it does with esp.A.r 
is.E.r ~      dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r +
              esp.B.r #esp.A.r + 
thick.A.r ~   dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r + maxc.r + slope.r

river.r  ~   dem.r
#mrvbf.r  ~ slope.r
lstsd.r  ~   dem.r
vdchn.r  ~    river.r
twi.r  ~    vdchn.r
evisd.r  ~      dem.r + vdchn.r + twi.r
vdchn.r  ~ slope.r
wdist.r  ~  lstm.r

# residual variances observed variables #indicate the residual variance in the measurements, as indicator of measurement error
thick.A ~~ thick.A

#factor variances                       #Only for endogenous variables, so not for age and vegetation index
tb.A.r ~~      oc.A.r + d.caco3.r
sat.A.r ~~     oc.A.r + d.caco3.r
bt.r ~~        esp.A.r + esp.B.r + d.caco3.r 
esp.A.r ~~     esp.B.r
oc.A.r ~~      thick.A.r
is.hydro.r ~~  esp.A.r + esp.B.r
is.E.r ~~      esp.A.r + esp.B.r
dem.r ~~ river.r
mrvbf.r ~~ slope.r
thick.A.r ~~ thick.A.r
river.r ~~ lstm.r
twi.r ~~    vdchn.r
d.caco3.r ~~ d.caco3.r
maxc.r ~~   vdchn.r
vdchn.r ~~   evisd.r

#Intercepts         #This code indicates that intercepts should be taken into account with the regression
tb.A.r ~1
sat.A.r ~1
bt.r ~1   
esp.A.r ~1
esp.B.r ~1
oc.A.r ~1 
d.caco3.r ~1  
is.hydro.r ~1 
is.E.r ~1     
thick.A.r ~1
'

#fitting
# std.ov = If TRUE, observed variables are standardized before entering the analysis
# std.lv = If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0.
fit1 <- lavaan(first_model, std.ov=T, std.lv = T, data=d)#, ordered = c("is.E","is.hydro")) 
fit1
# str(fit1)    #summary of the model fit
# str(fit1@Options)

summary(fit1, standardized=F, modindices = F, fit.measures=T) 
modi<-summary(fit1, standardized=F, modindices = T, fit.measures=F) 
modi<-modi[modi$mi>20 & !is.na(modi$mi),]
sum(modi$mi[modi$mi>20])
modi
#


inspect(fit1,"cov.ov")      #shows the covariance matrix with negative values





#################Prediction 1A##############
##setting up matrices

{B <-matrix(c(0,coef(fit1)["silt_1A_r~thick_1A_r"],coef(fit1)["OM_1A_r~thick_1A_r"],  
              coef(fit1)["thick_1A_r~silt_1A_r"],0,coef(fit1)["OM_1A_r~silt_1A_r"],
              coef(fit1)["thick_1A_r~OM_1A_r"],coef(fit1)["silt_1A_r~OM_1A_r"],0),
            ncol=3,nrow=3)
 B[is.na(B)] <- 0} #Matrix endogenous variables B
{I<-diag(nrow=3,
         ncol=3)} #Identity matrix
{A<-matrix(c(coef(fit1Al)["thick_1A_r~age_hy"],coef(fit1Al)["silt_1A_r~age_hy"],         
             coef(fit1Al)["OM_1A_r~age_hy"],coef(fit1Al)["thick_1A_r~VI"],
             coef(fit1Al)["silt_r~VI"],coef(fit1Al)["OM_1A_r~VI"]),
           nrow=3,ncol=2)
 A[is.na(A)] <- 0
} #Matrix exogenous variables

#Create test rasters
age <- as.vector(seq(from=0,to=140,by=140/49))
VI <- seq(from=100,to=40,by=-(100-40)/49)

age_df=matrix(nrow=50,ncol=50)
for (i in 1:50){
  age_df[i,]<-age}
VI_df=matrix(nrow=50,ncol=50)
for (i in 1:50){
  VI_df[,i]<-VI}