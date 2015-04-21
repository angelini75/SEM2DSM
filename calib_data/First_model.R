rm(list=ls())
#install.packages("lavaan")
library(lavaan)
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration/simplest_model")

d <-read.csv("calib.data-2.0.csv")[,-1]

############### STRUCTURAL EQUATION MODELLING ######################
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

first_model <- '            ##CREATING MODEL


# measurement model           ## Indicate in what manner the latent variables depend on the measured variables 
                              ## (arrows in conceptual model)
                              ## multiplication by 1 to set factor loading to 1 
 
thick.A.r =~ 1*thick.A 
tb.A.r =~ 1*tb.A 
sat.A.r =~ 1*sat.A
esp.A.r =~ 1*esp.A 
esp.B.r =~ 1*esp.B 
oc.A.r =~ 1*oc.A 
bt.r =~ 1*bt 
is.caco3.r =~ 1*is.caco3 
d.caco3.r =~ 1*d.caco3 
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
              oc.A.r + is.caco3.r + d.caco3.r
sat.A.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              oc.A.r +  d.caco3.r #+ is.caco3.r 
bt.r ~        dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r +
              esp.A.r + esp.B.r + d.caco3.r #+is.caco3.r #
esp.A.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              esp.B.r
esp.B.r ~     dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r 
            #  esp.A.r
oc.A.r ~      dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r + evim.r + evisd.r +
              thick.A.r
is.caco3.r ~  dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r
              
d.caco3.r ~   dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + lstm.r + lstsd.r #+
            #  is.caco3.r
is.hydro.r ~  dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r #+
            #  esp.A.r + esp.B.r
is.E.r ~      dem.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r +
              esp.A.r + esp.B.r
thick.A.r ~   dem.r + river.r + wdist.r + mrvbf.r + twi.r + vdchn.r + evim.r + evisd.r + maxc.r + slope.r
              
# evim.r ~      is.hydro.r + is.E.r
# evisd.r ~     is.hydro.r + is.E.r

# residual variances observed variables #indicate the residual variance in the measurements, as indicator of measurement error
thick.A ~~ 2*thick.A

#factor variances                       #Only for endogenous variables, so not for age and vegetation index
tb.A.r ~~      oc.A.r + is.caco3.r + d.caco3.r
sat.A.r ~~     oc.A.r + is.caco3.r + d.caco3.r
bt.r ~~        esp.A.r + esp.B.r + is.caco3.r + d.caco3.r 
esp.A.r ~~     esp.B.r
oc.A.r ~~      thick.A.r
d.caco3.r ~~   is.caco3.r
is.hydro.r ~~  esp.A.r + esp.B.r
is.E.r ~~      esp.A.r + esp.B.r

#Intercepts         #This code indicates that intercepts should be taken into account with the regression
tb.A.r ~1
sat.A.r ~1
bt.r ~1   
esp.A.r ~1
esp.B.r ~1
oc.A.r ~1 
is.caco3.r ~1 
d.caco3.r ~1  
is.hydro.r ~1 
is.E.r ~1     
thick.A.r ~1

'

#fitting
fit1 <- lavaan(first_model,std.lv = T, data=d) #Fit the model with the data

summary(fit1)    #summary of the model fit

summary(fit, standardized=T, modindices = F, fit.measures=T)   #summary of the fitting the model
inspect(fit1,"cov.ov")      #shows the covariance matrix with negative values

piet <- lav_partable_independence(fit1)
lav_partable_npar(piet)
lav_partable_ndat(piet)
lav_partable_df(piet)

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