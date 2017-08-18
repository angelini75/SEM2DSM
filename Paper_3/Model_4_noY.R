# Purpose        : Fit a SEM model with Argentinian data and apply in KS data
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 

#==============================================================================#
#### This code come from SEM_Arg2KS_1.0.R ####
#==============================================================================#

# Libraries ####
library(lavaan)
library(pastecs)
library(utils)
library(lattice)

rm(list=ls()[ls()!="t"])
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
ST.arg <- t(stat.desc(e,norm = TRUE)[c(9,13),])

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
STt.arg <- t(stat.desc(e,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
Arg <- std(e,STt.arg)
Arg[,1] <- e[,1] 

setwd("~/Documents/SEM2DSM1/Paper_3/data/")
d <- read.csv("KS.data-0.3.csv")[,c(-1)] 
name(d)
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
# remove outlayers
# d <- d[d$idp!=26058,]
# d <- d[d$idp!=22961,]
d <- cbind(d[1],
           d[,colnames(Arg)[2:10]],
           d[11:21])
d <- d[c(-156,-157),]
# Descriptive statistics and normality test. ####
round(stat.desc(d,norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
ST.ks <- t(stat.desc(d,norm = TRUE)[c(9,13),])

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
STt.ks <- t(stat.desc(d,norm = TRUE)[c(9,13),])

# standardised data set ####
std <- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] - st[i,1]) / st[i,2]
  }
  y
}
ks <- std(d,STt.ks)
ks[,1] <- d[,1] 
#######################################

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
## Measurement error #
#CEC.A ~~ 0.05 * CEC.A
#CEC.B ~~ 0.05 * CEC.B
# CEC.C ~~ 0.05 * CEC.C
OC.A ~~ 0.05 * OC.A
OC.B ~~ 0.05 * OC.B
OC.C ~~ 0.05 * OC.C
clay.A ~~ 0.05 * clay.A
clay.B ~~ 0.05 * clay.B
clay.C ~~ 0.05 *clay.C

#--------------------#
# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn + X + lstm 
clay.Ar ~ clay.Cr + evisd + lstm + ndwi.b #+ Y 
clay.Br ~ clay.Ar + clay.Cr + vdchn + twi + ndwi.b + X #+ Y

OC.Ar ~ clay.Ar + evisd + lstm + ndwi.b 
OC.Br ~ OC.Ar + clay.Br + evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ clay.Br + 0*OC.Br
CEC.Cr ~ clay.Cr + 0*OC.Cr

#------------------#
# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + 0*CEC.Cr
CEC.Cr ~~ CEC.Br
#OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 

#------------------#
# lavaan suggestions
#------------------#
clay.Br  ~   lstm
OC.Br  ~      dem
clay.Ar  ~    twi
OC.Cr  ~      dem
OC.Ar  ~   ndwi.a

OC.Ar  ~  clay.Br
OC.Br  ~  clay.Ar

CEC.Ar  ~     dem
CEC.Cr  ~   evisd
CEC.Br  ~     dem
CEC.Br  ~       X
CEC.Br  ~   evisd

OC.Cr ~~ clay.Cr
CEC.Ar ~~ clay.Br
#------------------#
'


# Model calibration ####
my.fit.lv.ML <- sem(model = my.model.lv,data = ks, meanstructure = FALSE, 
                    fixed.x = T)
inspect(my.fit.lv.ML,"cov.lv")
# Model evaluation ####
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = T)
my.fit.lv.ML
mod.ks <- modindices(my.fit.lv.ML,sort. = T)
mod.ks[mod.ks$mi & (mod.ks$op == "~"|mod.ks$op == "~~"),]# & 
# mod.ks$lhs != "CEC.Ar",]# & mod.ks$lhs == "clay.Cr",] 

write.csv(partable(my.fit.lv.ML)[partable(my.fit.lv.ML)$free != 0,], 
          "~/Documents/SEM2DSM1/Paper_3/data/partable_model4.csv")

# Model validation ####
#### Prediction ####
pre <- cbind(ks[1,], matrix(nrow=1,ncol= 9, data = NA,
                            dimnames = list(NULL,paste0(names(ks)[2:10],".p"))))
a <- pre[-(1:nrow(pre)),c(2:10)] #observed
b <- pre[-(1:nrow(pre)),c(2:10)] #predicted
v <- pre[-(1:nrow(pre)),c(2:10)] #variance(s)
resids <- pre[-(1:nrow(pre)),c(2:10)] #residuals
theta <- pre[-(1:nrow(pre)),c(2:10)] # for mean and median
Var <- pre[-(1:nrow(pre)),c(2:10)] # model variance (constant)

# Loop: cal is calibration data, pre is prediction place
# Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #

# Running Prediction @ i location #
# p is a matrix with the 10 external drivers
for(i in seq_along(ks[,1])){ 
  cal <- ks[-i,]
  pre[i,] <- ks[ i,]
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = cal, fixed.x = T,
                      estimator = "ML")
  
  # Matrix of Beta coefficients
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9]
  # Identity matrix (Kappa coefficients)
  I <- diag(nrow = 9, ncol = 9)
  # Matrix of Gamma coefficients
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:17]
  # Matrix of Psi coefficients (model error variance-covariance)
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9]
  # Matrix of measurement error (Epsylon)
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9]
  IB <- solve(I - B)
  
  # Running Prediction @ i location #
  # p is a matrix with the 10 external drivers
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 8, ncol = 1)           # by lavaan sequence
  # prediction
  pre[i,22:30] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- pre[i,c(2:10)] # observed values
  b[i,] <- pre[i,c(22:30)] # predicted values
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
Res <- cbind(pre[,1], unstd(pre[,2:10], STt.ks[2:10,]), unstd(pre[,22:30],
                                                              STt.ks[2:10,]))
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
  report[i-1,1] <- c(names(Res)[i])
  report[i-1,2] <- ME
  report[i-1,3] <- RMSE
  report[i-1,4] <- SS
}

for(i in 1:9){
  report$mean_theta[i] <- mean(theta[,i])
  report$median_th[i] <- median(theta[,i])
}
report
#d.stat <- read.csv("summary.calibdata.csv")
STt.ks <- as.data.frame(STt.ks)
STt.ks$SS <- NA 
for(i in seq_along(names(d))){
  STt.ks$SS[i] <- sum(( d[i] - STt.ks$mean[i])^2)
}

report$R2 <- 1 - (as.numeric(report$SS) / as.numeric(STt.ks$SS[2:10]))
report

# > report
#   Soil_property                 ME              RMSE               SS mean_theta  median_th          R2
# 1         CEC.A  -1.34191957183381  5.96747926364206 5626.50778439567  1.2358620 0.49772528 -0.03982376
# 2         CEC.B  -3.12455178339896  7.19852058698208 8187.35438531036  1.4499497 0.56871387 -0.13391956
# 3         CEC.C  -3.17991426471198  7.65434747281699 9257.06756706994  1.0286451 0.45802684 -0.06055332
# 4          OC.A -0.249497925624309 0.770620989528276 93.8293601012431  1.6001875 0.60191998 -0.21480605
# 5          OC.B  0.161131245639884 0.540332878165909 46.1296198378747  1.1583263 0.14104585 -0.11500986
# 6          OC.C  0.172189050861731 0.986648941090494 153.809229006889  1.0414483 0.03009904 -0.03223725
# 7        clay.A   1.27834695419596  8.63119289335688 11770.6035404488  0.9966927 0.28607540  0.05833345
# 8        clay.B  -3.44645803628208  9.95957589670949  15672.518022686  1.2215032 0.40126265  0.12027689
# 9        clay.C   5.01849584447639  11.3176134997709 20237.9633021711  1.2334849 0.20654433 -0.06307211

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
  report2[i,1] <- c(names(z)[i])
  report2[i,2:3] <- c(ME, RMSE)
}

report2$r2[1] <- rsq[1]
report2$r2[3] <- rsq[2]
report2$r2[5] <- rsq[3]

report2 <- report2[c(-4,-2),]
report2

#   Soil_property          ME      RMSE          r2
# 1          CECo -2.54879521  6.976593 0.009534237
# 3           OCo  0.02794079  0.787251 0.310940080
# 5         clayo  0.95012825 10.029606 0.147645073

################################################################################
library(lattice)

# reshape measured and predicted values for plotting
res <- rbind(as.matrix(Res[,c(2,11)]), as.matrix(Res[,c(3,12)]),
             as.matrix(Res[,c(4,13)]))
res <- rbind(res, as.matrix(Res[,c(5,14)]), as.matrix(Res[,c(6,15)]),
             as.matrix(Res[,c(7,16)]))
res <- rbind(res, as.matrix(Res[,c(8,17)]), as.matrix(Res[,c(9,18)]),
             as.matrix(Res[,c(10,19)]))

rownames(res) <- 1:(147*9)
colnames(res) <- c("observed","predicted")
res <- as.data.frame(res)
res$hor <- rep(c("A", "B", "C"), each= 147,3)
res$sp <- rep(c("CEC", "OC", "clay"), each= 147*3)
res$hor <- as.factor(res$hor)