# Purpose        : Fit a SEM model with Argentinian data and apply in KS data
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 

#==============================================================================#
             #### This code come from SEM_Arg2KS_1.0.R ####
#==============================================================================#

# Libraries ####
library(lavaan)
library(pastecs)
library(utils)

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

e <- read.csv("~/Documents/SEM2DSM1/Paper_2/data/calib.data-5.0.csv")[,c(-1,-20)]
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
# CEC.A ~~ 0.1 * CEC.A
# CEC.B ~~ 0.1 * CEC.B
# CEC.C ~~ 0.05 * CEC.C
# OC.A ~~ 0.2 * OC.A
# OC.B ~~ 0.6 * OC.B
# OC.C ~~ 0.6 * OC.C
# clay.A ~~ 0.3 * clay.A
# clay.B ~~ 0.14 * clay.B
# clay.C ~~ 0.12 *clay.C

#--------------------#

# Structural model (gamma and beta matrices)
#--------------------#
clay.Cr ~ dem + vdchn 
clay.Ar ~ clay.Cr + evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + dem + vdchn + twi + ndwi.b

OC.Ar ~ clay.Ar +
evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ OC.Br + clay.Br
CEC.Cr ~ OC.Cr + clay.Cr
#------------------#

# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + CEC.Cr
CEC.Cr ~~ CEC.Br
#OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 
#------------------#

# lavaan suggestions
#------------------#
# OC.Ar  ~     dem
# clay.Br  ~     dem
# clay.Br  ~    lstm
# clay.Cr  ~    lstm
# # # 
# CEC.Br  ~  ndwi.a
# CEC.Cr  ~     dem
# clay.Cr  ~  ndwi.a
# clay.Cr  ~   evisd
# # # 
# CEC.Ar ~~ clay.Br
# CEC.Ar ~~ clay.Cr
# OC.Cr ~~ clay.Cr
# CEC.Cr ~~ clay.Ar
#------------------#
'
# Model calibration ####
my.fit.lv.ML.e <- sem(model = my.model.lv.e,data = E, meanstructure = FALSE, 
                    fixed.x = T)

# Model evaluation ####
summary(my.fit.lv.ML.e, fit.measures=TRUE, rsquare = F)
mod.e <- modindices(my.fit.lv.ML.e,sort. = T)
mod.e[mod.e$mi>3 & (mod.e$op == "~"|mod.e$op == "~~"),] 


#==============================================================================#
                 #### This code come from SEM_KS2.0.R ####
#==============================================================================#


setwd("~/Documents/SEM2DSM1/Paper_3/data/")
d <- read.csv("KS.data-0.2.csv")[,c(-1)] 
name(d)
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
# remove outlayers
# d <- d[d$idp!=26058,]
# d <- d[d$idp!=22961,]
d <- cbind(d[1],
           d[,colnames(E)[2:10]],
           d[11:21])
# Descriptive statistics and normality test. ####
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# Based on normtest.W the following covariates need to be transformed
d$twi <- log10(d$twi)
d$vdchn <- log10(d$vdchn+10)
d$ndwi.a <- (d$ndwi.a+10)^.3
# OC as log10 of OC
# d$OC.A <- log10(d$OC.A)
# d$OC.B <- log10(d$OC.B)
# d$OC.C <- log10(d$OC.C)

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

# Model validation ####
my.fit.lv.ML.d <- sem(model = my.model.lv.e,data = D, meanstructure = FALSE, 
                      fixed.x = T)
# Model evaluation ####
summary(my.fit.lv.ML.d, fit.measures=TRUE, rsquare = F)
# mod <- modindices(my.fit.lv.ML.d,sort. = T)
# mod[mod$mi>3 & (mod$op == "~"|mod$op == "~"),] 

#### Prediction ####
pre <- cbind(D[1,], matrix(nrow=1,ncol= 9, data = NA,
                           dimnames = list(NULL,paste0(names(D)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)] # for mean and median
Var <- pre[-(1:nrow(pre)),c(2:10)] # model variance (constant)

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
  pre[i,20:28] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(20:28)] # predicted values
  v <- diag(IB%*%V%*%t(IB)+Th) # error variance (it is not diagonal!)
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v # theta
  
  # Error variance #
  Var[i,] <- diag(IB %*% V %*% t(IB))
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
# # back transform OC
# [,8]<- 10^(pred[,8]+(Var[6]*M$sd[6]^2)*0.5)
# Res[c(5:7)] <- 10^(Res[c(5:7)])
# Res[c(14)] <- 10^(Res[14] + (Var[4] * STt[4 + 1, 2]^2) * 0.5)
# Res[c(15)] <- 10^(Res[15] + (Var[5] * STt[5 + 1, 2]^2) * 0.5)
# Res[c(16)] <- 10^(Res[16] + (Var[6] * STt[6 + 1, 2]^2) * 0.5)

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
  report[i-1,1:4] <- c(names(Res)[i], ME, RMSE, SS)
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

# ST <- as.data.frame(ST)
# ST$SS <- NA 
# for(i in seq_along(names(d))){
#   ST$SS[i] <- sum(( d[i] - ST$mean[i])^2)
# }
# 
# report$R2 <- 1 - (as.numeric(report$SS) / as.numeric(ST$SS[2:10]))
# report
# 





# > report
# Soil_property                    ME              RMSE               SS mean_theta median_th         R2
# 1        clay.A -8.13906053375816e-16  11.2150962095849  24149.449534127   1.053935 0.4044461 0.04980644
# 2        clay.B -6.71491456279294e-14  14.1657025026762 38528.0884597109   1.519931 0.9465283 0.03410628
# 3        clay.C -3.72052318762931e-14  14.4363487350168  40014.367641412   1.300804 0.6880442 0.04407188
# 4         CEC.A  2.77179515414891e-14   9.1998603724659 16250.3867275907   1.192665 0.3270950 0.06196571
# 5         CEC.B -2.17980247206017e-14  9.18703482619708 16205.1089083696   1.014783 0.4078339 0.05027280
# 6         CEC.C  -1.1518563880486e-14  8.93336100797746 15322.5482685796   1.002315 0.3310788 0.01143328
# 7          OC.A -1.78329573330416e-15  1.20760452199481 279.995266856124   1.035802 0.3993229 0.04414984
# 8          OC.B -9.04658292721905e-16 0.284069138225369 15.4934928560839   1.600814 0.8762842 0.03178272
# 9          OC.C -4.65897943608216e-16 0.123405811509595 2.92396690835362   1.334580 0.6329824 0.04580599


# Statistics with log10(OC)
#   Soil_property                    ME              RMSE               SS mean_theta median_th         R2
#            OC.A    0.0737883121116024  1.25335116623517  301.61071601339  1.0557795 0.5246712 0.44353690
#            OC.B    0.0365582711336917 0.282103675337733 15.2798368586989  0.9899792 0.3855678 0.90497819
#            OC.C    0.0239350424903057  0.12762878580584 3.12750853752443  0.9782184 0.4097314 0.98839888
# Analysis by Soil Property
# plot mesured vs predicted combined ####par(mfrow = c(1,3), pty="s",mai=rep(0.7,4))
par(mfrow = c(1, 3), pty="s",mai=rep(0.7,4))
rsq<- NULL
CEC <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),
             as.matrix(Res[,c(4,13)]))
colnames(CEC) <- c("CECo","CECp")
rownames(CEC) <- 1:length(rownames(CEC))
CEC <- as.data.frame(CEC)
rsq[1] <- 1 - (sum((CEC$CECo - CEC$CECp)^2)/
                 sum((mean(CEC$CECo)-CEC$CECo)^2))
lim = round(c(min(c(CEC[,1],CEC[,2])), max(c(CEC[,1],CEC[,2]))))
plot(CEC[,2]~CEC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "CEC residuals", col = "dark red")
abline(0,1)
abline(lm(CEC[,2]~CEC[,1]),col = "blue")

OC <- rbind(as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),
            as.matrix(Res[,c(7,16)]))
colnames(OC) <- c("OCo","OCp")
rownames(OC) <- 1:length(rownames(OC))
OC <- as.data.frame(OC)
rsq[2] <- 1 - (sum((OC$OCo - OC$OCp)^2)/
                 sum((mean(OC$OCo)-OC$OCo)^2))
lim = round(c(min(c(OC[,1],OC[,2])), max(c(OC[,1],OC[,2]))))
plot(OC[,2]~OC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "OC residuals", col = "dark red")
abline(0,1)
abline(lm(OC[,2]~OC[,1]),col = "blue")

clay <- rbind(as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),
              as.matrix(Res[,c(10,19)]))

colnames(clay) <- c("clayo","clayp")
rownames(clay) <- 1:length(rownames(clay))
clay <- as.data.frame(clay)
rsq[3] <- 1 - (sum((clay$clayo - clay$clayp)^2)/
                 sum((mean(clay$clayo)-clay$clayo)^2))
lim = round(c(min(c(clay[,1],clay[,2])), max(c(clay[,1],clay[,2]))))
plot(clay[,2]~clay[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "Clay residuals", col = "dark red")
abline(0,1)
abline(lm(clay[,2]~clay[,1]),col = "blue")


# create report by soil property
report2 <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, r2 = NA)
z <- cbind(CEC, OC, clay)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2[i,1:3] <- c(names(z)[i], ME, RMSE)
}

report2$r2[1] <- rsq[1]
report2$r2[3] <- rsq[2]
report2$r2[5] <- rsq[3]

report2 <- report2[c(-4,-2),]
report2

# Soil_property                    ME             RMSE        r2
# 1          CECo -3.50606654854419e-14 13.3523240821413 0.1514665
# 3           OCo  -1.8666799153864e-15  9.1075788709089 0.0952693
# 5         clayo -1.05132460040218e-15 0.71977611128942 0.5102711

 
# report2 w/ log10(OC)
# Soil_property                    ME              RMSE        r2
#           OCo    0.0447605419118666 0.745377026859205 0.4748143




##########################################
#### Prediction ####
pre <- cbind(E[1,], matrix(nrow=1,ncol= 9, data = NA,
                           dimnames = list(NULL,paste0(names(D)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)] # for mean and median
Var <- pre[-(1:nrow(pre)),c(2:10)] # model variance (constant)
pb = txtProgressBar(min = 0, max = length(d[,1]), initial = 0, style = 3)

# Loop: cal is calibration data, pre is prediction place
for(i in seq_along(E[,1])){ 
  cal <- E[-i,]
  pre[i,] <- E[ i,]
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv.e, data = cal, fixed.x = T,
                      estimator = "ML")
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
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 7, ncol = 1)           # by lavaan sequence
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

# Accuracy measures ####
# Residuals #
Res <- cbind(pre[,1], unstd(pre[,2:10], STte[2:10,]), unstd(pre[,28:36],
                                                           STte[2:10,]))
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

# create report.e
report.e <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                       mean_theta = NA, median_th = NA)
for (i in 2:10) {
  # ME <- mean error 
  ME  <-  mean(Res[,i] - Res[,i + 9])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((Res[,i] - Res[,i + 9]) ^ 2))
  MSE <- mean((Res[,i] - Res[,i + 9]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((Res[,i] - Res[,i + 9]) ^ 2)
  # fill report.e table
  report.e[i-1,1:4] <- c(names(pre)[i], ME, RMSE, SS)
}

for(i in 1:9){
  report.e$mean_theta[i] <- mean(theta[,i])
  report.e$median_th[i] <- median(theta[,i])
}
report.e
#d.stat <- read.csv("summary.calibdata.csv")
STte <- as.data.frame(STte)
STte$SS <- NA 
for(i in seq_along(names(e))){
  STte$SS[i] <- sum(( e[i] - STte$mean[i])^2)
}

report.e$R2 <- 1 - (as.numeric(report.e$SS) / as.numeric(STte$SS[2:10]))
report.e

#### report2.e ####
rsq<- NULL
CEC <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),
             as.matrix(Res[,c(4,13)]))
colnames(CEC) <- c("CECo","CECp")
rownames(CEC) <- 1:length(rownames(CEC))
CEC <- as.data.frame(CEC)
rsq[1] <- 1 - (sum((CEC$CECo - CEC$CECp)^2)/
                 sum((mean(CEC$CECo)-CEC$CECo)^2))
lim = round(c(min(c(CEC[,1],CEC[,2])), max(c(CEC[,1],CEC[,2]))))
plot(CEC[,2]~CEC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "CEC residuals", col = "dark red")
abline(0,1)
abline(lm(CEC[,2]~CEC[,1]),col = "blue")

OC <- rbind(as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),
            as.matrix(Res[,c(7,16)]))
colnames(OC) <- c("OCo","OCp")
rownames(OC) <- 1:length(rownames(OC))
OC <- as.data.frame(OC)
rsq[2] <- 1 - (sum((OC$OCo - OC$OCp)^2)/
                 sum((mean(OC$OCo)-OC$OCo)^2))
lim = round(c(min(c(OC[,1],OC[,2])), max(c(OC[,1],OC[,2]))))
plot(OC[,2]~OC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "OC residuals", col = "dark red")
abline(0,1)
abline(lm(OC[,2]~OC[,1]),col = "blue")

clay <- rbind(as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),
              as.matrix(Res[,c(10,19)]))

colnames(clay) <- c("clayo","clayp")
rownames(clay) <- 1:length(rownames(clay))
clay <- as.data.frame(clay)
rsq[3] <- 1 - (sum((clay$clayo - clay$clayp)^2)/
                 sum((mean(clay$clayo)-clay$clayo)^2))
lim = round(c(min(c(clay[,1],clay[,2])), max(c(clay[,1],clay[,2]))))
plot(clay[,2]~clay[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "Clay residuals", col = "dark red")
abline(0,1)
abline(lm(clay[,2]~clay[,1]),col = "blue")

# create report by soil property
report2.e <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, r2 = NA)
z <- cbind(CEC, OC, clay)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2.e[i,1:3] <- c(names(z)[i], ME, RMSE)
}

report2.e$r2[1] <- rsq[1]
report2.e$r2[3] <- rsq[2]
report2.e$r2[5] <- rsq[3]

report2.e <- report2.e[c(-4,-2),]
report2.e
#------------------------------------------------------------------------------#
report
report.e
report2
report2.e

#------------------------------------------------------------------------------#
# performance of MvLR
mod.ceca <-  lm(CEC.A  ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.cecb <-  lm(CEC.B  ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.cecc <-  lm(CEC.C  ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.oca <-   lm(OC.A   ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.ocb <-   lm(OC.B   ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.occ <-   lm(OC.C   ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.claya <- lm(clay.A ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.clayb <- lm(clay.B ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)
mod.clayc <- lm(clay.C ~ dem + vdchn + evisd + lstm + ndwi.b + twi + ndwi.a, E)

mod <- list(mod.ceca, mod.cecb, mod.cecc,
            mod.oca, mod.ocb, mod.occ,
            mod.claya, mod.clayb, mod.clayc)
P <- cbind(D[1:10], D[2:10])
names(P)[11:19] <- paste0(names(P)[11:19],".p")
P[11:19] <- NA
for(i in seq_along(mod)){
  P[,10+i] <- predict.lm(mod[[i]],D)
}

Res <- P[2:10] - P[11:19]

Res <- cbind(P[,1], unstd(P[,2:10], STt[2:10,]), unstd(P[,11:19], STt[2:10,]))

# # back transform OC
# [,8]<- 10^(pred[,8]+(Var[6]*M$sd[6]^2)*0.5)
# Res[c(5:7)] <- 10^(Res[c(5:7)])
# Res[c(14)] <- 10^(Res[14] + (Var[4] * STt[4 + 1, 2]^2) * 0.5)
# Res[c(15)] <- 10^(Res[15] + (Var[5] * STt[5 + 1, 2]^2) * 0.5)
# Res[c(16)] <- 10^(Res[16] + (Var[6] * STt[6 + 1, 2]^2) * 0.5)

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
report.m <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
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
  report.m[i-1,1:4] <- c(names(Res)[i], ME, RMSE, SS)
}

for(i in 1:9){
  report.m$mean_theta[i] <- mean(theta[,i])
  report.m$median_th[i] <- median(theta[,i])
}
report.m
#d.stat <- read.csv("summary.calibdata.csv")
STt <- as.data.frame(STt)
STt$SS <- NA 
for(i in seq_along(names(d))){
  STt$SS[i] <- sum(( d[i] - STt$mean[i])^2)
}

report.m$R2 <- 1 - (as.numeric(report.m$SS) / as.numeric(STt$SS[2:10]))
report.m

# ST <- as.data.frame(ST)
# ST$SS <- NA 
# for(i in seq_along(names(d))){
#   ST$SS[i] <- sum(( d[i] - ST$mean[i])^2)
# }

par(mfrow = c(1, 3), pty="s",mai=rep(0.7,4))
rsq<- NULL
CEC <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),
             as.matrix(Res[,c(4,13)]))
colnames(CEC) <- c("CECo","CECp")
rownames(CEC) <- 1:length(rownames(CEC))
CEC <- as.data.frame(CEC)
rsq[1] <- 1 - (sum((CEC$CECo - CEC$CECp)^2)/
                 sum((mean(CEC$CECo)-CEC$CECo)^2))
lim = round(c(min(c(CEC[,1],CEC[,2])), max(c(CEC[,1],CEC[,2]))))
plot(CEC[,2]~CEC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "CEC residuals", col = "dark red")
abline(0,1)
abline(lm(CEC[,2]~CEC[,1]),col = "blue")

OC <- rbind(as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),
            as.matrix(Res[,c(7,16)]))
colnames(OC) <- c("OCo","OCp")
rownames(OC) <- 1:length(rownames(OC))
OC <- as.data.frame(OC)
rsq[2] <- 1 - (sum((OC$OCo - OC$OCp)^2)/
                 sum((mean(OC$OCo)-OC$OCo)^2))
lim = round(c(min(c(OC[,1],OC[,2])), max(c(OC[,1],OC[,2]))))
plot(OC[,2]~OC[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "OC residuals", col = "dark red")
abline(0,1)
abline(lm(OC[,2]~OC[,1]),col = "blue")

clay <- rbind(as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),
              as.matrix(Res[,c(10,19)]))

colnames(clay) <- c("clayo","clayp")
rownames(clay) <- 1:length(rownames(clay))
clay <- as.data.frame(clay)
rsq[3] <- 1 - (sum((clay$clayo - clay$clayp)^2)/
                 sum((mean(clay$clayo)-clay$clayo)^2))
lim = round(c(min(c(clay[,1],clay[,2])), max(c(clay[,1],clay[,2]))))
plot(clay[,2]~clay[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "Clay residuals", col = "dark red")
abline(0,1)
abline(lm(clay[,2]~clay[,1]),col = "blue")


# create report by soil property
report2.m <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, r2 = NA)
z <- cbind(CEC, OC, clay)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2.m[i,1:3] <- c(names(z)[i], ME, RMSE)
}

report2.m$r2[1] <- rsq[1]
report2.m$r2[3] <- rsq[2]
report2.m$r2[5] <- rsq[3]

report2.m <- report2.m[c(-4,-2),]

report2
report2.m

report
report.m

# > report
# Soil_property                    ME              RMSE               SS mean_theta  median_th           R2
# 1         CEC.A -6.96944370183068e-16  5.92575761273541 5548.10731901288  1.1242340 0.35577341 -0.025334724
# 2         CEC.B  8.09469459201117e-16  6.35290999896399 6376.79554187998  1.2644777 0.43767222  0.116836425
# 3         CEC.C  3.14951124252346e-16  7.00648334397663 7756.34779820866  1.0245463 0.42709644  0.111379459
# 4          OC.A  -1.6298749500623e-16 0.795062648356146  99.875689140152  1.6376447 0.61203458 -0.293087703
# 5          OC.B  9.87622317870477e-16 0.537465199123197 45.6412767624289  1.1764330 0.14196800 -0.103206002
# 6          OC.C  4.19432874178757e-16 0.975510384854676 150.356040731572  1.0094245 0.02771345 -0.009062377
# 7        clay.A  6.97229831008229e-16  8.54949047506805 11548.8184065581  0.9925004 0.29778703  0.076076607
# 8        clay.B  2.13561614990411e-15   9.5805213002397 14502.2493647268  1.2651786 0.41392886  0.185965908
# 9        clay.C  1.11490842489967e-15  10.1559374973902 16296.6044992417  1.1330173 0.20891365  0.143961992
# > report.m
# Soil_property                    ME              RMSE               SS mean_theta median_th           R2
# 1         CEC.A -1.92472126053567e-14  6.33891709943952 6348.73545898355  0.9781116 0.4329402 -0.173297224
# 2         CEC.B  5.04574262064301e-14  6.76741819523969 7236.07594662327  1.0036666 0.4609274 -0.002170865
# 3         CEC.C  2.46664105170629e-14  7.16494518149537 8111.15743370572  0.8848770 0.3396347  0.070730027
# 4          OC.A  -6.6887554248406e-15 0.792615602792903 99.2618400189397  1.0038438 0.4155480 -0.285140216
# 5          OC.B  4.71879970571004e-15 0.535563312381524  45.318833727913  1.0016106 0.3570525 -0.095412156
# 6          OC.C  3.11571651498876e-15 0.986367890365801  153.72161519286  0.9876415 0.3513463 -0.031649262
# 7        clay.A -1.28836677393766e-14  8.60207633905194 11691.3233401746  0.9625760 0.3098189  0.064675991
# 8        clay.B  1.24344978758018e-14  9.33735421967365 13775.4170441379  0.9836433 0.4102715  0.226764150
# 9        clay.C  8.43296714516505e-14  10.0306771577704   15897.08851046  1.0000000 0.4844208  0.164948012
# > report2
# Soil_property                   ME              RMSE        r2
# 1          CECo 1.42327382152539e-16  6.44372763725186 0.1550577
# 3           OCo 4.14661040803837e-16 0.790066234784347 0.3060031
# 5         clayo 1.31457454802148e-15  9.45204108854313 0.2429860
# > report2.m
# Soil_property                   ME              RMSE         r2
# 1          CECo 1.86259313652657e-14  6.76550708808196 0.06856314
# 3           OCo 3.81931387507722e-16 0.793303046446629 0.30030499
# 5         clayo 2.79598390761905e-14  9.34159846082658 0.26057338











