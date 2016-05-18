############### #### ### ## # - SEM for second paper - # ## ### #### ###########
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

# Libraries ####
library(lavaan)
library(pastecs)
library(utils)

rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
setwd("~/big/SEM_2nd_paper/")

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

d <- read.csv("calib.data-5.0.csv")[,c(-1,-20)] #remove water variable 

# Descriptive statistics and normality test. ####
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
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
# statistics
round(stat.desc(D,norm = TRUE),0)

# SEM ####
# Model without latent variables (CFA) ####
## disabled

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
# my.fit.ML <- sem(model = my.model,data = D, meanstructure = FALSE,
# fixed.x = T, estimator = "ML")
# summary(my.fit.ML, fit.measures=TRUE, rsquare = T)
# mod <- modindices(my.fit.ML)
# mod[mod$mi>5,] # suggestion where mi is higher than 10 (most significant mi)
# as.data.frame(lavaan::fitMeasures(my.fit.ML,fit.measures = "all"))

# Model with latent variables ####
## Model (re)specification
my.model.lv <- '

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
## Measurement error
CEC.A ~~ 0.1 * CEC.A
CEC.B ~~ 0.1 * CEC.B
CEC.C ~~ 0.1 * CEC.C
OC.A ~~ 0.1 * OC.A
OC.B ~~ 0.1 * OC.B
OC.C ~~ 0.1 * OC.C
clay.A ~~ 0.1 * clay.A
clay.B ~~ 0.1 * clay.B
clay.C ~~ 0.1 * clay.C
#--------------------#

# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + river + vdchn + X + Y 
clay.Ar ~ clay.Cr + 
          evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
          vdchn + twi + river + Y + ndwi.b

OC.Ar ~ clay.Ar +
        evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
        evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ OC.Br + clay.Br
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
CEC.Cr ~~ clay.Cr 
clay.Br  ~     dem
OC.Br ~~ clay.Br
CEC.Br  ~  ndwi.a
clay.Ar  ~ clay.Br
#------------------#
'
# Model calibration ####
my.fit.lv.ML <- sem(model = my.model.lv,data = D, meanstructure = FALSE, 
                    fixed.x = T, estimator = "ML", se="bootstrap")

# Model evaluation ####
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = F)

# Model respecification: modification indices ####
fitMeasures(my.fit.lv.ML,fit.measures = "gfi")
fitMeasures(my.fit.lv.ML,fit.measures = "srmr")
mod <- modindices(my.fit.lv.ML,sort. = T)
mod[mod$mi>10,] # suggestion where mi is higher than 10 (most significant mi)


# Cross-validation #####
# same model than before
my.model.lv <- '
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
## Measurement error
CEC.A ~~ 0.1 * CEC.A
CEC.B ~~ 0.1 * CEC.B
CEC.C ~~ 0.1 * CEC.C
OC.A ~~ 0.1 * OC.A
OC.B ~~ 0.1 * OC.B
OC.C ~~ 0.1 * OC.C
clay.A ~~ 0.1 * clay.A
clay.B ~~ 0.1 * clay.B
clay.C ~~ 0.1 * clay.C
#--------------------#

# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + river + vdchn + X + Y 
clay.Ar ~ clay.Cr + 
evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
vdchn + twi + river + Y + ndwi.b

OC.Ar ~ clay.Ar +
evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ OC.Br + clay.Br
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
CEC.Cr ~~ clay.Cr 
clay.Br  ~     dem
OC.Br ~~ clay.Br
CEC.Br  ~  ndwi.a
clay.Ar  ~ clay.Br
#------------------#
'
# Element definition
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
for(i in seq_along(D[,1])){ 
  cal <- D[-i,]
  pre[i,] <- D[ i,]
  
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = cal, fixed.x = T,
                      estimator = "ML")
  #pre <- pre[,c("id.p",colnames(inspect(my.fit.lv.ML, "est")$theta),
  #                              names(pre)[21:29])]
  
  # Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #
  # Matrix of Beta coefficients
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9] 
  # Identity matrix (Kappa coefficients)
  I <- diag(nrow = 9, ncol = 9)
  # Matrix of Gamma coefficients
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:19]
  # Matrix of Psi coefficients (model error variance-covariance)
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9] 
  # Matrix of measurement error (Epsylon)
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9] 
  IB <- solve(I - B)
  # Running Prediction @ i location #
  # p is a matrix with the 10 external drivers
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 10, ncol = 1)           # by lavaan sequence
  # prediction 
  pre[i,28:36] = t(IB %*% A %*% p) # key equation 
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(28:36)] # predicted values
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
Res <- cbind(pre[,1], unstd(pre[,2:10], STt[2:10,]), unstd(pre[,28:36],
                                                           STt[2:10,]))
# plot residuals
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))
for (i in 2:10) {
  lim = c(0, max(c(Res[,i],Res[,i+9])))
  plot(Res[,i+9] ~ Res[,i], main = paste(names(Res)[i]), xlab = "measured",
       ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
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


# plot mesured vs predicted combined ####
par(mfrow = c(3, 1), pty="s",mai=rep(0.7,4))

r2<- NULL
CEC <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),
             as.matrix(Res[,c(4,13)]))
colnames(CEC) <- c("CECo","CECp")
rownames(CEC) <- 1:length(rownames(CEC))
CEC <- as.data.frame(CEC)
r2[1] <- 1 - (sum((CEC$CECo - CEC$CECp)^2)/
                 sum((mean(CEC$CECo)-CEC$CECo)^2))
plot(CEC[,2]~CEC[,1])
abline(lm(CEC[,2]~CEC[,1]),col = "red")

OC <- rbind(as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),
            as.matrix(Res[,c(7,16)]))
colnames(OC) <- c("OCo","OCp")
rownames(OC) <- 1:length(rownames(OC))
OC <- as.data.frame(OC)
r2[2] <- 1 - (sum((OC$OCo - OC$OCp)^2)/
             sum((mean(OC$OCo)-OC$OCo)^2))
plot(OC[,2]~OC[,1])
abline(lm(OC[,2]~OC[,1]),col = "red")


clay <- rbind(as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),
              as.matrix(Res[,c(10,19)]))

colnames(clay) <- c("clayo","clayp")
rownames(clay) <- 1:length(rownames(clay))
clay <- as.data.frame(clay)
r2[3] <- 1 - (sum((clay$clayo - clay$clayp)^2)/
             sum((mean(clay$clayo)-clay$clayo)^2))
plot(clay[,2]~clay[,1])
abline(lm(clay[,2]~clay[,1]),col = "red")


# create report by soil property
report2 <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                     r2 = NA)
z <- cbind(CEC, OC, clay)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2[i,1:4] <- c(names(z)[i], ME, RMSE, SS)
}

report2$r2[1] <- 1 - (as.numeric(report$SS[1]) / sum((mean(z[,1])-z[,1])^2))
report2$r2[3] <- 1 - (as.numeric(report$SS[3]) / sum((mean(z[,3])-z[,3])^2))
report2$r2[5] <- 1 - (as.numeric(report$SS[5]) / sum((mean(z[,5])-z[,5])^2))

report2 <- report2[c(-4,-2),]

# Covariation assessment ####

unstd(x = t(B),st = STt[2:10,2])
apply(d[,2:10],2,mean)%*%t(apply(d[,2:10],2,mean))
round(cor(d[,2:10])-B,3)





# Validation ####
