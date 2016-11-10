# Purpose        : Fit a SEM model with Argentinian data and apply in KS data
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 

#### This code come from SEM_Arg2KS_1.0.R #######
# Libraries ####
library(lavaan)

rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
# chose one


# Dictionary of elements in this script ######
# d = calibration dataset. It comes from replacement_of_NAs.Rm 
#     (different versions: 5.0 to 5.3)
# ST = original mean and standard deviation of all variables
# STt = mean and standard deviation of transformed data
# nor = normalisation funcion (x-mean)/sd
# D = transformed and normalised data
# my.model = lavaan syntax
# my.fit = model fitted
# mod = modification indices (for respecification)
#------------------------------------------------#

e <- read.csv("~/Documents/SEM2DSM1/Paper_2/data/calib.data-5.0.csv")[,c(-1,-20)] #remove water variable 
# Descriptive statistics and normality test. ####
round(stat.desc(e,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
STe <- t(stat.desc(e,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
e$wdist <- e$wdist^0.5
e$maxc <- (e$maxc+20000)^2
e$slope <- e$slope^0.25
e$twi <- log10(e$twi)
e$vdchn <- log10(e$vdchn+10)
e$ndwi.a <- (e$ndwi.a+10)^.3
# New statistics
round(stat.desc(e,norm = TRUE),3)
# New mean and sd
STte <- t(stat.desc(e,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
E <- std(e,STte)
E[,1] <- e[,1] 

# Model with latent variables ####
## Model (re)specification
my.model.lv.e <- '
# Measurement model (lamda and epsilon)
#--------------------#
CEC.Ar =~ 1*CEC.A
CEC.Br =~ 1*CEC.B
CEC.Cr =~ 1*CEC.C
OC.Ar =~ 1*OC.A
OC.Br =~ 1*OC.B
OC.Cr =~ 1*OC.C
clay.Ar =~ 1*clay.A
clay.Br =~ 1*clay.B
clay.Cr =~ 1*clay.C
## Measurement error #
CEC.A ~~ 0.1 * CEC.A
CEC.B ~~ 0.1 * CEC.B
CEC.C ~~ 0.05 * CEC.C
OC.A ~~ 0.2 * OC.A
OC.B ~~ 0.6 * OC.B
OC.C ~~ 0.6 * OC.C
clay.A ~~ 0.3 * clay.A
clay.B ~~ 0.14 * clay.B
clay.C ~~ 0.12 *clay.C

#--------------------#

# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn 
clay.Ar ~ clay.Cr + 
evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
vdchn + twi + ndwi.b

OC.Ar ~ clay.Ar +
evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ clay.Br
CEC.Cr ~ clay.Cr
#------------------#

# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + CEC.Cr
CEC.Cr ~~ CEC.Br
OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 
#------------------#

# lavaan suggestions
#------------------#
OC.Ar  ~     dem
clay.Br  ~     dem
clay.Br  ~    lstm
clay.Cr  ~    lstm
# # 
CEC.Br  ~  ndwi.a
CEC.Cr  ~     dem
clay.Cr  ~  ndwi.a
clay.Cr  ~   evisd
# # 
CEC.Ar ~~ clay.Br
CEC.Ar ~~ clay.Cr
OC.Cr ~~ clay.Cr
CEC.Cr ~~ clay.Ar
#------------------#
'
# Model calibration ####
my.fit.lv.ML.e <- sem(model = my.model.lv.e,data = E, meanstructure = FALSE, 
                    fixed.x = T)

# Model evaluation ####
summary(my.fit.lv.ML.e, fit.measures=TRUE, rsquare = F)
#------------------------------------------------------------#####

#### This code come from SEM_KS2.0.R
setwd("~/Documents/SEM2DSM1/Paper_3/data/")
d <- read.csv("calib-data.KS.0.1.csv")[,c(-1)] 
name(d)
d <- d[,c(-11:-20)]
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
# remove outlayers
d <- d[d$idp!=26058,]
d <- d[d$idp!=22961,]


# Descriptive statistics and normality test. ####
round(stat.desc(d[,-20],norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST <- t(stat.desc(d[,c(-20)],norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
d$twi.1 <- log10(d$twi.1)
d$vdchn.1 <- log10(d$vdchn.1+10)
d$ndwi.a <- (d$ndwi.a+10)^.3

# New statistics
d <- d[,-20:-21]
names(d)[12:13] <- c("twi", "vdchn")
round(stat.desc(d,norm = TRUE),3)
# New mean and sd
STt <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
D <- std(d,STt)
D[,1] <- d[,1] 

#### Prediction ####
pre <- cbind(D[1,], matrix(nrow=1,ncol= 9, data = NA,
                           dimnames = list(NULL,paste0(names(D)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)] # for mean and median
Var <- pre[-(1:nrow(pre)),c(2:10)] # model variance (constant)
pb = txtProgressBar(min = 0, max = length(d[,1]), initial = 0, style = 3)

# Loop: cal is calibration data, pre is prediction place

# Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #
# Matrix of Beta coefficients
B <- inspect(my.fit.lv.ML.e, "est")$beta[1:9,1:9]
# Identity matrix (Kappa coefficients)
I <- diag(nrow = 9, ncol = 9)
# Matrix of Gamma coefficients
A <- inspect(my.fit.lv.ML.e, "est")$beta[1:9,10:16]
# Matrix of Psi coefficients (model error variance-covariance)
V <- inspect(my.fit.lv.ML.e, "est")$psi[1:9,1:9]
# Matrix of measurement error (Epsylon)
Th <- inspect(my.fit.lv.ML.e, "est")$theta[1:9,1:9]
IB <- solve(I - B)
# Running Prediction @ i location #
# p is a matrix with the 10 external drivers
for(i in seq_along(D[,1])){ 
  pre[i,] <- D[ i,]
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 7, ncol = 1)           # by lavaan sequence
  # prediction
  pre[i,20:26] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(20:26)] # predicted values
  v <- diag(IB%*%V%*%t(IB)+Th) # error variance (it is not diagonal!)
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v # theta
  
  # Error variance #
  Var[i,] <- diag(IB %*% V %*% t(IB))
  
  #bar time
  setTxtProgressBar(pb, i)
}

# Model variance
summary(Var)
Var <- apply(Var, MARGIN = 2, FUN = mean)

# function to unstandardise the data
unstd<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}

# Accuracy measures ####
# Residuals #
Res <- cbind(pre[,1], unstd(pre[,2:10], STt[2:10,]), unstd(pre[,20:28],
                                                           STt[2:10,]))
# plot residuals
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))

lim.cec <- c(min(cbind(Res[,2:4],Res[,11:13])), max(cbind(Res[,2:4],Res[,11:13])))
lim.oc <- c(min(cbind(Res[,5:7],Res[,14:16])), max(cbind(Res[,5:7],Res[,14:16])))
lim.clay <- c(min(cbind(Res[,8:10],Res[,17:19])), max(cbind(Res[,8:10],Res[,17:19])))

lim <- data.frame(min=NA,max=NA)
lim[1:3,1] <- lim.cec[1]
lim[1:3,2] <- lim.cec[2]
lim[4:6,1] <- lim.oc[1]
lim[4:6,2] <- lim.oc[2]
lim[7:9,1] <- lim.clay[1]
lim[7:9,2] <- lim.clay[2]


for (i in 2:10) {
  limi = c(lim[i-1,1],lim[i-1,2])
  plot(Res[,i+9] ~ Res[,i], main = paste(names(Res)[i]), xlab = "measured",
       ylab = "predicted", col = "dark red", xlim = limi, ylim = limi)
  abline(0,1)
  abline(lm(Res[,i+9] ~ Res[,i]), col = "blue")
}

# create report
report <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                     mean_theta = NA, median_th = NA)
for (i in 2:10) {
  # ME <- mean error 
  ME  <-  mean(Res[,i] - Res[,i + 9])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((Res[,i] - Res[,i + 9]) ^ 2))
  MSE <- mean((Res[,i] - Res[,i + 9]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((Res[,i] - Res[,i + 9]) ^ 2)
  # fill report table
  report[i-1,1:4] <- c(names(pre)[i], ME, RMSE, SS)
}

for(i in 1:9){
  report$mean_theta[i] <- mean(theta[,i])
  report$median_th[i] <- median(theta[,i])
}
report
#d.stat <- read.csv("summary.calibdata.csv")
STt <- as.data.frame(STt)
STt$SS <- NA 
for(i in seq_along(names(d))){
  STt$SS[i] <- sum(( d[i] - STt$mean[i])^2)
}

report$R2 <- 1 - (as.numeric(report$SS) / as.numeric(STt$SS[2:10]))
report
















