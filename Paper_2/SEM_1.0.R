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
# chose one
#setwd("~/big/SEM2DSM1/Paper_2/data/")
setwd("~/Documents/SEM2DSM1/Paper_2/data/")

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
d <- d[,c(-12:-16,-20,-21,-25)]

meltp <- melt(unique(d,id.vars = c("id.p")))

ggplot(data = meltp,
       aes(x = value, fill=variable)) + geom_histogram() + 
  facet_wrap( ~ variable, scales = "free_x")
d$H <- NA
e <- data.frame(d[,c(1,2,5,8,11:20)])
e$H <- "A"
names(e)[2:4] <- c("CEC.B","OC.B","clay.B")
e <- rbind(e,d[,c(1,3,6,9,11:20)])
e$H[is.na(e$H)] <- "B"

names(e)[2:4] <- c("CEC.C","OC.C","clay.C")
e <- rbind(e,d[,c(1,4,7,10,11:20)])
e$H[is.na(e$H)] <- "C"
names(e)[2:4] <- c("CEC","OC","Clay")

meltp <- melt(e,id.vars = c("id.p","H"))
ggplot(data = meltp[meltp$variable=="CEC" |
                      meltp$variable=="OC" |
                      meltp$variable=="Clay",],
       aes(x = value, fill = H)) + geom_density(alpha = 0.4) + 
  facet_wrap( ~ variable,scales = "free")
ggplot(data = unique(meltp[!(meltp$variable=="CEC" |
                      meltp$variable=="OC" |
                      meltp$variable=="Clay"),]),
       aes(x = value)) + geom_density(alpha = 0.4) + 
  facet_wrap( ~ variable,scales = "free")
#############
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

# WEPAL document WEPAL_ISE900.pdf ####
# search for literature
# CEC (NH4) sd = 1.68 cmol+/kg CV 10.5%
# CO sd = 0.218 % CV 11.8%
# clay sd = 3.11 % CV 14.1%

# measurement error as a constant
(1.68/STt[2:4,2])^2 # standardised measurement error for CEC
(0.218/STt[5:7,2])^2 # standardised measurement error for OC
(3.11/STt[8:10,2])^2 # standardised measurement error for clay

# measurement error as relative to the mean value
(0.105)^2 # standardised measurement error for CEC
(0.118)^2 # standardised measurement error for OC
(0.141)^2 # standardised measurement error for clay

# model error variance when measurement errors are zero
# they are constrains for measurement error
# CEC.Ar            0.389    0.032   12.124    0.000
# CEC.Br            0.362    0.030   12.124    0.000
# CEC.Cr            0.568    0.047   12.124    0.000
# OC.Ar             0.730    0.060   12.124    0.000
# OC.Br             0.922    0.076   12.124    0.000
# OC.Cr             0.946    0.078   12.124    0.000
# clay.Ar           0.841    0.069   12.124    0.000
# clay.Br           0.357    0.029   12.124    0.000
# clay.Cr           0.568    0.047   12.124    0.000

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
clay.Cr ~ dem + river + vdchn + X + 0*Y 
clay.Ar ~ clay.Cr + 
evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
vdchn + twi + river + 0*Y + ndwi.b

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
  OC.Ar  ~       Y
clay.Ar  ~       X
  OC.Ar  ~     dem
clay.Br  ~       X
clay.Br  ~     dem
clay.Br  ~    lstm
  OC.Br  ~       X
  OC.Ar  ~   river
# # 
 CEC.Cr  ~   river
 CEC.Br  ~  ndwi.a
 CEC.Cr  ~       X
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
my.fit.lv.ML <- sem(model = my.model.lv,data = D, meanstructure = FALSE, 
                    fixed.x = T)

# Model evaluation ####
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = F)
# model Rsquare

# r2 <- rbind(r2,round(inspect(my.fit.lv.ML, "rsquare")[1:9] * 
#               inspect(my.fit.lv.ML, "rsquare")[10:18],3))
# r2
# Model respecification: modification indices ####
fitMeasures(my.fit.lv.ML,fit.measures = 
              c("chisq","df","pvalue","cfi","rmsea","gfi", "srmr"))
mod <- modindices(my.fit.lv.ML,sort. = T)
mod[mod$mi>3 & (mod$op == "~~"|mod$op == "~~"),] 


# CROSS-VALIDATION #####
# same model than before

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
  lim = c(min(c(Res[,i],Res[,i+9])), max(c(Res[,i],Res[,i+9])))
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

setwd("~/Documents/SEM2DSM1/Paper_2/reports/")
#setwd("~/big/SEM2DSM1/Paper_2/reports/")
#write.csv(report, "report.byhor.csv")
# R2
# 1 0.17997897
# 2 0.49938652
# 3 0.45497635
# 4 0.23497330
# 5 0.03238344
# 6 0.01771871
# 7 0.14873541
# 8 0.60412970
# 9 0.42054296

# plot mesured vs predicted combined ####
par(mfrow = c(1,3), pty="s",mai=rep(0.7,4))

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
#write.csv(report2, "report.bysp.csv")
# Soil_property                   ME              RMSE        r2
# 1          CECo -0.00419372399590419  4.29902877759584 0.5262497
# 3           OCo 3.66189158698424e-05 0.249115654848832 0.9079135
# 5         clayo -0.00703797144384494  5.44721980811741 0.7165984


# sem plot ####
library(semPlot)
semPaths(my.fit.lv.ML,what = "est",style = "LISREL", layout = "circle")

# COVARIATION ASSESSMENT ####
# write.csv(round(residuals(my.fit.lv.ML, "raw")$cov[1:9,1:9], 3), 
#           "residual.matrix.csv")

# How to estimate Sigma.hat and residual matrix by Yves Rosseel
# Two things: you need the full matrices, including the
# x's (not just the first three rows/cols), and S is divided by N (not N-1):

attach(inspect(my.fit.lv.ML, "est"))
IB.inv <- solve(diag( nrow(beta) ) - beta)
Sigma.hat <- lambda %*% IB.inv %*% psi %*% t(IB.inv) %*% t(lambda) + 
  inspect(my.fit.lv.ML, "est")$theta
# should be the same as fitted(my.fit)$cov

# sample cov (divided by N, instead of N-1)
S <- cov(D[,lavNames(my.fit.lv.ML)]) * (nobs(my.fit.lv.ML) - 1) /
    nobs(my.fit.lv.ML)

# residuals
colors <- colorRampPalette(c('white', 'black'))(256)
DIF <- abs(S - Sigma.hat)
DIF[lower.tri(DIF)] <- NA
S[lower.tri(S)] <- NA
Sigma.hat[lower.tri(Sigma.hat)] <- NA

dif <- levelplot(DIF[1:9,9:1],
                 col.regions=colors,
                 at=seq(0,0.4,0.05),
                 xlab = NULL, ylab = NULL,                                                                                                                                                                                                                
                 scales=list(x=list(rot=90)),
                 names.attr="SEM")
s <- levelplot(S[1:9,9:1],
               col.regions=colors,
               at=seq(0,1.1,0.1),
               xlab = NULL, ylab = NULL,                                                                                                                                                                                                                
               scales=list(x=list(rot=90)))
sigmah <- levelplot(Sigma.hat[1:9,9:1],
                    col.regions=colors,
                    at=seq(0,1.1,0.1),
                    xlab = NULL, ylab = NULL,                                                                                                                                                                                                                
                    scales=list(x=list(rot=90)))

# round((Sigma.hat)[1:9,1:9],3) 
# round(resid(my.fit.lv.ML)$cov[1:9,1:9],3)

# MLR comparison ####
attach(D)
mod.ceca <- lm(CEC.A ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.cecb <- lm(CEC.B ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.cecc <- lm(CEC.C ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.oca <- lm(OC.A ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.ocb <- lm(OC.B ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.occ <- lm(OC.C ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.claya <- lm(clay.A ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.clayb <- lm(clay.B ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                         twi + ndwi.a)
mod.clayc <- lm(clay.C ~ dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                  twi + ndwi.a)

beta.lm <- matrix(data = 0, nrow = 19, ncol = 19,dimnames = 
                    list(rownames(beta), colnames(beta)))
beta.lm[1,10:19] <- mod.ceca$coefficients[2:11]
beta.lm[2,10:19] <- mod.cecb$coefficients[2:11]
beta.lm[3,10:19] <- mod.cecc$coefficients[2:11]
beta.lm[4,10:19] <- mod.oca$coefficients[2:11]
beta.lm[5,10:19] <- mod.ocb$coefficients[2:11]
beta.lm[6,10:19] <- mod.occ$coefficients[2:11]
beta.lm[7,10:19] <- mod.claya$coefficients[2:11]
beta.lm[8,10:19] <- mod.clayb$coefficients[2:11]
beta.lm[9,10:19] <- mod.clayc$coefficients[2:11]

psi.lm <- matrix(data = 0, nrow = 19, ncol = 19,dimnames = 
                    list(rownames(beta), colnames(beta)))
psi.lm[1,1] <- summary(mod.ceca)$sigma^2
psi.lm[2,2] <- summary(mod.cecb)$sigma^2
psi.lm[3,3] <- summary(mod.cecc)$sigma^2
psi.lm[4,4] <- summary(mod.oca)$sigma^2
psi.lm[5,5] <- summary(mod.ocb)$sigma^2
psi.lm[6,6] <- summary(mod.occ)$sigma^2
psi.lm[7,7] <- summary(mod.claya)$sigma^2
psi.lm[8,8] <- summary(mod.clayb)$sigma^2
psi.lm[9,9] <- summary(mod.clayc)$sigma^2
psi.lm[10:19,10:19] <- psi[10:19,10:19]
#

IB.lm.inv <- solve(diag( nrow(beta.lm) ) - beta.lm)
lambda.lm <- diag( nrow(beta.lm) )
Sigma.hat.lm <- lambda.lm %*% IB.lm.inv %*% psi.lm %*% t(IB.lm.inv) %*% t(lambda.lm) 
rownames(Sigma.hat.lm) <- colnames(Sigma.hat)
colnames(Sigma.hat.lm) <- colnames(Sigma.hat)
Sigma.hat.lm[lower.tri(Sigma.hat.lm)] <- NA

sigmah.lm <- levelplot(Sigma.hat.lm[1:9,9:1],
                       col.regions=colors,
                       at=seq(0,1.1,0.1),
                       xlab = NULL, ylab = NULL,                                                                                                                                                                                                                
                       scales=list(x=list(rot=90)))
DIF.lm <- abs(Sigma.hat.lm - S)
DIF.lm[lower.tri(DIF.lm)] <- NA

dif.lm <- levelplot(DIF.lm[1:9,9:1],
                    col.regions=colors,
                    at=seq(0,0.4,0.05),
                     xlab = NULL, ylab = NULL,                                                                                                                                                                                                                
                    scales=list(x=list(rot=90)))

# plotting ####
print(dif, split = c(1,1,2,1),)
print(dif.lm, split=c(2,1,2,1), newpage=FALSE)
print(c(dif, dif.lm))

print(c(sigmah,s, sigmah.lm), main =NULL)

# print(s, split = c(1,1,3,1))
# print(sigmah, split=c(2,1,3,1), newpage=FALSE)
# print(sigmah.lm, split=c(3,1,3,1), newpage=FALSE)

plot(pre$CEC.A.p~pre$OC.A.p)
plot(D$CEC.A~D$OC.A)
plot(predict.lm(mod.ceca)~predict.lm(mod.oca))

plot(pre$CEC.B.p~pre$CEC.C.p)
plot(D$CEC.B~D$CEC.C)
plot(predict.lm(mod.cecb)~predict.lm(mod.cecc))

# SRMR

sum(sum(()))^.5

### CROSS-VALIDATION MLR ####
# comparison with multivariate 
attach(D)
mod.sp <- lm(cbind(CEC.A, CEC.B, CEC.C,
                   OC.A, OC.B, OC.C, 
                   clay.A, clay.B, clay.C) ~
               dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
               twi + ndwi.a)

summary(mod.sp)
class(mod.sp)
summary(mod.ceca)
summary(manova(mod.sp))
library(car)
Manova(mod.sp, type = "II")
mod.ceca$residuals - mod.sp$residuals[,1] #same

# Cross-validation
# P <- predicted, C <- calibration, V <- validation, R <- residuals
P <- as.data.frame(D[1:10])
P[,11:19] <- NA
names(P)[11:19] <- paste0(names(D)[2:10],".pred")
attach(C)
for(i in seq_along(D[,1])){
  C <- D[-i,]
  V <- D[i,]
  # calibration LOO
  mod.sp <- lm(cbind(CEC.A, CEC.B, CEC.C,
                     OC.A, OC.B, OC.C, 
                     clay.A, clay.B, clay.C) ~
                 dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
                 twi + ndwi.a, C)
  
  P[i,11:19] <- predict(mod.sp, V)
}
####
unstd<- function(x, st){
  y <- x
  for(i in seq_along(names(x))){
    y[,i] <- (x[,i] * st[i,2]) + st[i,1]
  }
  y
}
####
theta.MLR.mean <- rep(NA,10)
theta.MLR.median <- rep(NA,10)
R <- P[,1:10]
R[2:10] <- P[,2:10]-P[,11:19]
for(i in 2:10){
  theta.MLR.mean[i] <- mean((R[,i]^2)/as.numeric(summary(mod.sp)[[i-1]][6])^2)
  theta.MLR.median[i] <- median((R[,i]^2)/as.numeric(summary(mod.sp)[[i-1]][6])^2)
}


P[2:10] <- unstd(P[,2:10], STt[2:10,])
P[11:19] <- unstd(P[,11:19], STt[2:10,])
R <- P[,1:10]
R[2:10] <- P[,2:10]-P[,11:19]

# Accuracy measures of MLR
# create report

reportMLR <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                     mean_theta = NA, median_th = NA)
for (i in 2:10) {
  # ME <- mean error 
  ME  <-  mean(R[,i])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((R[,i]) ^ 2))
  MSE <- mean((R[,i]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((R[,i]) ^ 2)
  # fill report table
  reportMLR[i-1,1:4] <- c(names(P)[i], ME, RMSE, SS)
}


for(i in 1:9){
  reportMLR$mean_theta[i] <- mean(theta[,i])
  reportMLR$median_th[i] <- median(theta[,i])
}
reportMLR
#d.stat <- read.csv("summary.calibdata.csv")
STt <- as.data.frame(STt)
STt$SS <- NA 
for(i in seq_along(names(d))){
  STt$SS[i] <- sum(( d[i] - STt$mean[i])^2)
}

reportMLR$R2 <- 1 - (as.numeric(reportMLR$SS) / as.numeric(STt$SS[2:10]))
reportMLR


# VALIDATION ####
setwd("~/Documents/SEM2DSM1/Paper_2/data/")
val <- read.csv("val.data.csv")[,-1]
val <- val[,-19]
# Statistics, normalization and standardisation ####
# Descriptive statistics and normality test. #
round(stat.desc(val[,-1],norm = TRUE),3)
# Soil properties does not present strong deviation from normality.
# But some covariates need to be transformed. First, we store original mean and 
# sd in ST
STv <- t(stat.desc(val,norm = TRUE)[c(9,13),])[-1,]

# Based on normtest.W the following covariates need to be transformed
val$wdist <- val$wdist^0.5
val$maxc <- (val$maxc+20000)^2
val$slope <- val$slope^0.25
val$twi <- log10(val$twi)
val$vdchn <- log10(val$vdchn+10)
val$ndwi.a <- (val$ndwi.a+10)^.3
# New statistics
round(stat.desc(val[,2:27],norm = TRUE),3)
# New mean and sd
STtv <- t(stat.desc(val[,1:27],norm = TRUE)[c(9,13),])

# standardised data set ####
VAL <- val
VAL[,2:27] <- std(val[,2:27],STt[2:27,])

# statistics
round(stat.desc(VAL[,2:27],norm = TRUE),2)

# Model validation ####
val.fit.lv.ML <- sem(model = my.model.lv,data = VAL, meanstructure = FALSE, 
                    fixed.x = T)

summary(val.fit.lv.ML, fit.measures=TRUE, rsquare = F)

fitMeasures(val.fit.lv.ML,fit.measures = 
              c("chisq","df","pvalue","cfi","rmsea","gfi", "srmr"))
mod.v <- modindices(val.fit.lv.ML,sort. = T)
mod.v[mod.v$mi>3 & (mod.v$op == "~~"|mod.v$op == "~"),] 

# model prediction @ val location ####
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

# element definition
VAL[,28:36] <- NA
names(VAL)[28:36] <- paste0(names(VAL)[2:10], ".p")
a <- VAL[-(1:nrow(VAL)),c(2:10)] #observed
b <- VAL[-(1:nrow(VAL)),c(28:36)] #predicted
v <- VAL[-(1:nrow(VAL)),c(2:10)] #variance(s)
resids <- VAL[-(1:nrow(VAL)),c(2:10)] #residuals
theta <- VAL[-(1:nrow(VAL)),c(2:10)] # for mean and median
# Running Prediction @ i location #
for(i in seq_along(VAL[,1])){ 
  # p is a matrix with the 10 external drivers
  p = as.vector(as.matrix(VAL[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 10, ncol = 1)           # by lavaan sequence
  # prediction 
  VAL[i,28:36] = t(IB %*% A %*% p) # key equation 
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  a[i,] <- VAL[i,c(2:10)] # observed values
  b[i,] <- VAL[i,c(28:36)] # predicted values
  v <- diag(IB%*%V%*%t(IB)+Th) # error variance (it is not diagonal!)
  resids[i,] <- a[i,] - b[i,] # residuals
  theta[i,] <- (resids[i,]^2)/v # theta 
}

# Accuracy measures ####
# Residuals #
Res.v <- cbind(VAL[,1], 
               unstd(VAL[,2:10],  STt[2:10,]), 
               unstd(VAL[,28:36], STt[2:10,])
               )
# plot residuals
par(mfrow = c(3, 3), pty="s",mai=rep(0.7,4))
for (i in 2:10) {
  lim = c(min(c(Res.v[,i],Res.v[,i+9])), max(c(Res.v[,i],Res.v[,i+9])))
  plot(Res.v[,i+9] ~ Res.v[,i], main = paste(names(Res.v)[i]), xlab = "measured",
       ylab = "predicted", col = "dark red", xlim = lim, ylim = lim)
  abline(0,1)
  abline(lm(Res.v[,i+9] ~ Res.v[,i]), col = "blue")
}

# create report ####
report.v <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
                     mean_theta = NA, median_th = NA)
for (i in 2:10) {
  # ME <- mean error 
  ME  <-  mean(Res.v[,i] - Res.v[,i + 9])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((Res.v[,i] - Res.v[,i + 9]) ^ 2))
  MSE <- mean((Res.v[,i] - Res.v[,i + 9]) ^ 2)
  # SS (Sum of squares)
  SS <- sum((Res.v[,i] - Res.v[,i + 9]) ^ 2)
  # fill report table
  report.v[i-1,1:4] <- c(names(VAL)[i], ME, RMSE, SS)
}

for(i in 1:9){
  report.v$mean_theta[i] <- mean(theta[,i])
  report.v$median_th[i] <- median(theta[,i])
}
report.v
STtv <- as.data.frame(STtv)
STtv$SS <- NA 
for(i in seq_along(names(val))){
  STtv$SS[i] <- sum(( val[,i] - STtv$mean[i])^2)
}

report.v$R2 <- 1 - (as.numeric(report.v$SS) / as.numeric(STtv$SS[2:10]))
report.v

# plot mesured vs predicted combined ####
par(mfrow = c(1,3), pty="s",mai=rep(0.7,4))

rsq.v<- NULL
CEC.v <- rbind(as.matrix(Res.v[,c(2,11)]), as.matrix(Res.v[,c(3,12)]),
             as.matrix(Res.v[,c(4,13)]))
colnames(CEC.v) <- c("CECo","CECp")
#rownames(CEC.v) <- 1:length(rownames(CEC.v))
CEC.v <- as.data.frame(CEC.v)
rsq.v[1] <- 1 - (sum((CEC.v$CECo - CEC.v$CECp)^2)/
                 sum((mean(CEC.v$CECo)-CEC.v$CECo)^2))
lim = round(c(min(c(CEC.v[,1],CEC.v[,2])), max(c(CEC.v[,1],CEC.v[,2]))))
plot(CEC.v[,2]~CEC.v[,1], xlim = lim, ylim= lim, xlab = "measured",
       ylab = "predicted", main = "CEC residuals", col = "dark red")
abline(0,1)
abline(lm(CEC.v[,2]~CEC.v[,1]),col = "blue")

OC.v <- rbind(as.matrix(Res.v[,c(5,14)]), as.matrix(Res.v[,c(6,15)]),
            as.matrix(Res.v[,c(7,16)]))
colnames(OC.v) <- c("OCo","OCp")
OC.v <- as.data.frame(OC.v)
rsq.v[2] <- 1 - (sum((OC.v$OCo - OC.v$OCp)^2)/
             sum((mean(OC.v$OCo)-OC.v$OCo)^2))
lim = round(c(min(c(OC.v[,1],OC.v[,2])), max(c(OC.v[,1],OC.v[,2]))))
plot(OC.v[,2]~OC.v[,1], xlim = lim, ylim= lim, xlab = "measured",
       ylab = "predicted", main = "OC residuals", col = "dark red")
abline(0,1)
abline(lm(OC.v[,2]~OC.v[,1]),col = "blue")

clay.v <- rbind(as.matrix(Res.v[,c(8,17)]), as.matrix(Res.v[,c(9,18)]),
              as.matrix(Res.v[,c(10,19)]))

colnames(clay.v) <- c("clayo","clayp")
clay.v <- as.data.frame(clay.v)
rsq.v[3] <- 1 - (sum((clay.v$clayo - clay.v$clayp)^2)/
             sum((mean(clay.v$clayo)-clay.v$clayo)^2))
lim = round(c(min(c(clay.v[,1],clay.v[,2])), max(c(clay.v[,1],clay.v[,2]))))
plot(clay.v[,2]~clay.v[,1], xlim = lim, ylim= lim, xlab = "measured",
       ylab = "predicted", main = "Clay residuals", col = "dark red")
abline(0,1)
abline(lm(clay.v[,2]~clay.v[,1]),col = "blue")


# create report by soil property ####
report2.v <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, r2 = NA)
z <- cbind(CEC.v, OC.v, clay.v)
for (i in c(1,3,5)) {
  # ME <- mean error 
  ME  <-  mean(z[,i] - z[,i + 1])
  # RMSE (root mean squared error)
  RMSE <- sqrt(mean((z[,i] - z[,i + 1]) ^ 2))
  MSE <- mean((z[,i] - z[,i + 1]) ^ 2)
  # fill report table
  report2.v[i,1:3] <- c(names(z)[i], ME, RMSE)
}

report2.v$r2[1] <- rsq.v[1]
report2.v$r2[3] <- rsq.v[2]
report2.v$r2[5] <- rsq.v[3]

report2.v <- report2.v[c(-4,-2),]
report2.v


# MAPS PREDICTION ####


# extract covariates values ####
setwd("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/")
# libraries
library(raster)
library(maptools)
library(sp)
library(rgdal)

# location of predictions
pred <- read.csv("mask_231m2.csv")
coordinates(pred) <- ~X+Y

#define crs
wgs84 <- CRS("+init=epsg:4326")
posgar98 <- CRS("+init=epsg:22175")
modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# assign projection
proj4string(pred) <- modis
pred <- spTransform(pred, posgar98)

# sdat files (dem covariates) 
files <- list.files(pattern=".sdat$")
header <- gsub(".sdat", "", files)
header <- c("dem", "river", "wdist","maxc","mrvbf","slope","twi","vdchn","water") 

# tif files (modis)
files_m <- list.files(pattern=".tif$")
# set names of covariates
header_m <- c("lstm", "lstsd", "evim", "evisd", "ndwi.a", "ndwi.b", "ndwi.bsd")

# extract values from files (.sdat)
stack <- list()
for(i in seq_along(files)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files[i])
  proj4string(stack[[i]]) <- posgar98
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header[i]
}  

## extract values from modis files 
stack <- list()
# reproject endo to modis projection
pred <- spTransform(pred, modis)
for(i in seq_along(files_m)) {
  pred@data[,length(pred@data)+1] <- NULL
  stack[[i]] <- readGDAL(files_m[i])
  proj4string(stack[[i]]) <- modis # change projection
  pred@data[,length(pred@data)+1] <- over(pred, stack[[i]])[,1]
  stack <- list()
  names(pred@data)[length(pred@data)] <- header_m[i]
}  
pred <- spTransform(pred, posgar98)
pred.df <- as.data.frame(pred)
pred.df <- pred.df[complete.cases(pred.df),] 

# transform
pred.df$wdist <- pred.df$wdist^0.5
pred.df$maxc <- (pred.df$maxc+20000)^2
pred.df$slope <- pred.df$slope^0.25
pred.df$twi <- log10(pred.df$twi)
pred.df$vdchn <- log10(pred.df$vdchn+10)
pred.df$ndwi.a <- (pred.df$ndwi.a+10)^.3

# standardise
as.data.frame(rownames(STt))
pred.df <- pred.df[,-1]
pred.df <- pred.df[,-9]
pred.st <- std(x = pred.df,st = STt[11:27,])
name(pred.st)
# prediction
pred.st[,18:26] <- NA
names(pred.st)[18:26] <- names(D)[2:10]
pred.st <- as.matrix(pred.st) # conversion to matrix to improve processing speed
pb = txtProgressBar(min = 0, max = length(pred.st[,1]), initial = 0, style = 3)

for(i in seq_along(pred.st[,1])){
  p = as.vector(as.matrix(pred.st[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 10, ncol = 1)           # by lavaan sequence
  # prediction 
  pred.st[i,18:26] = t(IB %*% A %*% p) # 
  setTxtProgressBar(pb, i)
}

# unstandardize soil properties
pred.un <- pred.st[,rownames(STt)[2:27]]
for(i in seq_along(pred.un[1,])){
  ST <- STt[2:27,]
  pred.un[,i] <- (pred.un[,i] * ST$std.dev[i]) + ST$mean[i]
}
pred.un <- as.data.frame(pred.un)

# rasterize results 
library(sp)
library(raster)
coordinates(pred.un) <- ~X+Y
proj4string(pred.un) <- posgar98
pred.un <- spTransform(pred.un, modis)
#spplot(pred.un[,2])

y <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/mod13q1_tot_mean.tif")
y[!y %in% NA] <- 0
proj4string(y) <- modis
names(y)<-"mask"
res(y)

r <- rasterize(x = pred.un[1:9],y = y,background= NA)
spplot(r[[1]])
ADE<- readShapePoly("/media/marcos/L0135974_DATA/UserData/BaseARG/study area/ADE_MODIS.shp")
plot(ADE)
r<-mask(x = r,mask = ADE)
TWI <- raster("/media/marcos/L0135974_DATA/UserData/BaseARG/COVARIATES/modelling/TWI250.sdat")
proj4string(TWI) <- posgar98

rp <- projectRaster(from = r,to = TWI,res = res(TWI), crs = posgar98, method = 'bilinear')
library(rasterVis)
levelplot(rp[[2:4]],)
levelplot(rp[[5:7]],)
levelplot(rp[[8:10]],)
s<- as.data.frame(summary(rp[[2:8]],digits=4))
names(s) <- names(rp[[2:8]])

raster::NAvalue(rp)<--99999
writeRaster(x = rp[[2:10]],filename ="paper2.tif", overwrite=T,bylayer=TRUE,suffix=names(rp)[2:10])







