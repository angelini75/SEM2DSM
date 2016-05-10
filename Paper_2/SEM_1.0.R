############### #### ### ## # - SEM for second paper - # ## ### #### ###############
# Purpose        : Load, standardise and create a SE model
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : Gerard?
# Status         : 
# Note           : 
# sessionInfo(@RStudio desktop)  lenovo ThinkPad T430 (4 cores)
# R version 3.0.2 (2013-09-25)
# Platform: x86_64-pc-linux-gnu (64-bit)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base   
# other attached packages:
#   [1] questionr_0.5

library(lavaan)
library(pastecs)
library(utils)

rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
setwd("~/big/SEM_2nd_paper/")

##### Dictionary of elements in this script ######
# d = calibration dataset. It comes from replacement_of_NAs.Rm (different versions: 5.0 to 5.3)
# ST = original mean and standard deviation of all variables
# STt = mean and standard deviation of transformed data
# nor = normalisation funcion (x-mean)/sd
# D = transformed and normalised data
# my.model = lavaan syntax
# my.fit = model fitted
# mod = modification indices (for respecification)
##################################################

d <- read.csv("calib.data-5.0.csv")[,c(-1,-20)] #remove water variable 

# Descriptive statistics and normality test.
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and sd in ST
ST <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
d$wdist <- d$wdist^0.5
d$maxc <- (d$maxc+20000)^2
d$slope <- d$slope^0.25
d$twi <- log10(d$twi)
d$vdchn <- log10(d$vdchn+10)
d$ndwi.a <- (d$ndwi.a+10)^.3
# New statistics
round(stat.desc(d,norm = TRUE),3)
# New mean and sd
STt <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# normalisation
nor <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
# standardised data set
D <- nor(d,STt)
D[,1] <- d[,1] 
# statistics
round(stat.desc(D,norm = TRUE),0)

#### Begining SEM ####
# first, a model without latent variables (no measurement error)
# names of variables
name(D)
# my.model <- '
# clay.C ~ a01*river + a02*X + a03*Y + a04*vdchn + a05*dem 
# clay.A ~ b01*clay.C + 
#          a07*evisd + a08*lstm + a09*ndwi.b 
# clay.B ~ b02*clay.A + b03*clay.C + 
#          a10*vdchn + a11*twi + a12*river + a13*Y + a14*ndwi.b
# 
# OC.A ~ b04*clay.A +
#        a15*evisd + a16*lstm + a17*ndwi.b  
# OC.B ~ b05*OC.A + b06*clay.B + 
#        a18*evisd + a19*lstm + a20*ndwi.a + a21*vdchn
# 
# CEC.A ~ b07*OC.A + b08*clay.A 
# CEC.B ~ b09*OC.B + b10*clay.B
# CEC.C ~ b11*clay.C
# 
# CEC.A ~~ CEC.B + CEC.C
# CEC.C ~~ CEC.B
# OC.C ~ OC.B + clay.C
# 
# #modifications in order of relevances
# clay.C ~~  CEC.C # 1st suggestion. CFI 0.933 and RMSEA 0.068 means that
#                  # the model is moderately good
# clay.B  ~    dem # 2nd suggestion, which is reasonable 
# clay.B ~~   OC.B #
# clay.B ~~  CEC.B
# '
# my.fit.ML <- sem(model = my.model,data = D, meanstructure = FALSE, fixed.x = T, estimator = "ML")
# summary(my.fit.ML, fit.measures=TRUE, rsquare = T)
# mod <- modindices(my.fit.ML)
# mod[mod$mi>5,] # suggestion where mi is higher than 10 (most significant mi)
# as.data.frame(lavaan::fitMeasures(my.fit.ML,fit.measures = "all"))

###### Second, a model with measurement error

my.model.lv <- '
# Measurement model
CEC.Ar =~ 1*CEC.A
CEC.Br =~ 1*CEC.B
CEC.Cr =~ 1*CEC.C
OC.Ar =~ 1*OC.A
OC.Br =~ 1*OC.B
OC.Cr =~ 1*OC.C
clay.Ar =~ 1*clay.A
clay.Br =~ 1*clay.B
clay.Cr =~ 1*clay.C

# Structural model
clay.Cr ~ dem + river + vdchn + X + Y 
clay.Ar ~ clay.Cr + 
          evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
          vdchn + twi + river + Y + ndwi.b

OC.Ar ~ clay.Ar +
        evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
        evisd + lstm + ndwi.a + vdchn

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ OC.Br + clay.Br
CEC.Cr ~ clay.Cr

CEC.Ar ~~ CEC.Br + CEC.Cr
CEC.Cr ~~ CEC.Br
OC.Cr ~ OC.Br + clay.Cr

# Measurement error
CEC.A ~~ 0.1 * CEC.A
CEC.B ~~ 0.1 * CEC.B
CEC.C ~~ 0.1 * CEC.C
OC.A ~~ 0.1 * OC.A
OC.B ~~ 0.1 * OC.B
OC.C ~~ 0.1 * OC.C
clay.A ~~ 0.1 * clay.A
clay.B ~~ 0.1 * clay.B
clay.C ~~ 0.1 * clay.C

# suggestions
#CEC.Cr ~~ clay.Cr # CFI .939 RMSEA .063 GFI .940 SMRM .039 
'
my.fit.lv.ML <- sem(model = my.model.lv,data = D, meanstructure = FALSE, fixed.x = T, estimator = "ML")
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = F)
fitMeasures(my.fit.lv.ML,fit.measures = "gfi")
fitMeasures(my.fit.lv.ML,fit.measures = "srmr")
mod <- modindices(my.fit.lv.ML,sort. = T)
mod[mod$mi>10,] # suggestion where mi is higher than 10 (most significant mi)

# reorder variables in D as they appear in fit model.
#D <- D[,c("id.p",colnames(inspect(my.fit.lv.ML, "est")$theta))]

##### cross-validation #####


#### Model ####
my.model.lv <- '
# Measurement model
CEC.Ar =~ 1*CEC.A
CEC.Br =~ 1*CEC.B
CEC.Cr =~ 1*CEC.C
OC.Ar =~ 1*OC.A
OC.Br =~ 1*OC.B
OC.Cr =~ 1*OC.C
clay.Ar =~ 1*clay.A
clay.Br =~ 1*clay.B
clay.Cr =~ 1*clay.C

# Structural model
clay.Cr ~ river + X + Y + vdchn + dem 
clay.Ar ~ clay.Cr + 
evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
vdchn + twi + river + Y + ndwi.b

OC.Ar ~ clay.Ar +
evisd + ndwi.b + lstm +  lstsd + evim + evisd + dem + wdist + mrvbf + vdchn + twi  
OC.Br ~ OC.Ar + clay.Br + 
evisd + lstm + ndwi.a + vdchn

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ OC.Br + clay.Br
CEC.Cr ~ clay.Cr

#   CEC.Ar ~~ CEC.Br + CEC.Cr
#   CEC.Cr ~~ CEC.Br
OC.Cr ~ OC.Br + clay.Cr

# Measurement error
#   CEC.A ~~ 0.1 * CEC.A
#   CEC.B ~~ 0.1 * CEC.B
#   CEC.C ~~ 0.1 * CEC.C
#   OC.A ~~ 0.1 * OC.A
#   OC.B ~~ 0.1 * OC.B
#   OC.C ~~ 0.1 * OC.C
#   clay.A ~~ 0.1 * clay.A
#   clay.B ~~ 0.1 * clay.B
#   clay.C ~~ 0.1 * clay.C
#   
# suggestions
#  CEC.Cr ~~ clay.Cr # CFI .939 RMSEA .063
'
pre <- cbind(D[1,],
             matrix(nrow=1,ncol= 9, data = NA,dimnames = list(NULL,paste0(names(D)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)]
Var <- pre[-(1:nrow(pre)),c(2:10)]
pb = txtProgressBar(min = 0, max = length(d[,1]), initial = 0, style = 3)

for(i in seq_along(D[,1])){ 
  cal <- D[-i,]
  pre[i,] <- D[ i,]
  
  #### Fiting ####
  my.fit.lv.ML <- sem(model = my.model.lv,data = D, meanstructure = FALSE, fixed.x = T, estimator = "ML")
  #pre <- pre[,c("id.p",colnames(inspect(my.fit.lv.ML, "est")$theta),names(pre)[21:29])]
  
  #### Prediction @ i ####
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9] # matrix of coeff. latent state variables
  I <- diag(nrow = 9, ncol = 9) # Identity matrix
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:19] # matrix of coeff of external drivers
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9] # matrix of predicted error variance
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9] # matrix of measurement error
  IB <- solve(I - B)
  ################# Running Prediction ###########################
  # (IB%*%A%*%p) product of matrices per pixel (equation 4 paper)
  p = as.vector(as.matrix(pre[i,colnames(A)])) #values of covariates ordered respecting A sequence
  p = matrix(p, nrow = 10, ncol = 1)
  pre[i,21:29] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(28:36)] # predicted values
  v <- diag(IB%*%V%*%t(IB)+Th) # error variance (now is not diagonal)
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v # theta 
  
  ##### Error variance #######
  Var[i,] <- diag(IB %*% V %*% t(IB))
  
  #bar time
  setTxtProgressBar(pb, i)
}

summary(Var)
Var <- apply(Var, MARGIN = 2, FUN = mean)
apply(theta, MARGIN = 2,FUN = mean)
apply(theta, MARGIN = 2,FUN = median)

# function to unstandardise the data
unnor<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}

# 
Res <- cbind(pre[,1], unnor(pre[,2:10], STt[2:10,]), unnor(pre[,28:36], STt[2:10,]))

par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))

for (i in 2:10) {
  lim = c(0, max(c(Res[,i],Res[,i+9])))
  plot(Res[,i+9] ~ Res[,i], main = paste(names(Res)[i]), xlab = "measured",
       ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
  abline(0,1)
  abline(lm(Res[,i+9] ~ Res[,i]), col = "blue")
}

report <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA, mean_theta = NA, median_th = NA)
for (i in 2:10) {
  ######################## ME <- mean error 
  ME  <-  mean(Res[,i] - Res[,i + 9])
  ############################################ RMSE (root mean squared error)
  RMSE <- sqrt(mean((Res[,i] - Res[,i + 9]) ^ 2))
  MSE <- mean((Res[,i] - Res[,i + 9]) ^ 2)
  ############################################ SS (Sum of squares)
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
report$R2 <- 1 - (as.numeric(report$SS) / STt)
report


par(mfrow = c(1, 1), pty="s",mai=rep(0.7,4))


CEC <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),as.matrix(Res[,c(4,13)]))
CEC <- as.data.frame(CEC)
plot(CEC[,2]~CEC[,1])
abline(lm(CEC[,2]~CEC[,1]),col = "red")

OC <- rbind(as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),as.matrix(Res[,c(7,16)]))
OC <- as.data.frame(OC)
plot(OC[,2]~OC[,1])
abline(lm(OC[,2]~OC[,1]),col = "red")


clay <- rbind(as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),as.matrix(Res[,c(10,19)]))
clay <- as.data.frame(clay)
plot(clay[,2]~clay[,1])
abline(lm(clay[,2]~clay[,1]),col = "red")


plot(d$CEC.C~d$clay.C)









