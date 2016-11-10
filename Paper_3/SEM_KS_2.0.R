############### #### ### ## # - SEM for Third Paper - # ## ### #### ###########
# Purpose        : Load, standardise and create a SE model
# Maintainer     : Marcos Angelini  (marcos.angelini@wur.nl); 
# Contributions  : 
# Status         : 
# Note           : 
# SessionInfo()
# R version 3.3.2 (2016-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.1 LTS
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] tools_3.3.2

# Libraries ####
library(lavaan)
library(pastecs)
library(utils)

rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
# chose one
#setwd("~/big/SEM2DSM1/Paper_2/data/")
setwd("~/Documents/SEM2DSM1/Paper_3/data/")

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


d <- read.csv("calib-data.KS.0.1.csv")[,c(-1)] 
name(d)
d <- d[,c(-11:-20)]
names(d)[5:10] <- c("CEC.A","CEC.B","CEC.C","OC.A","OC.B","OC.C")
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

# OC as log10 of OC
# d$oc.A <- log10(d$oc.A)
# d$oc.B <- log10(d$oc.B)
# d$oc.C <- log10(d$oc.C)

# d$clay.A <- log10(d$clay.A)
# d$clay.B <- log10(d$clay.B)
# d$clay.C <- log10(d$clay.C)
# 
# d$cec.A <- log10(d$cec.A)
# d$cec.B <- log10(d$cec.B)
# d$cec.C <- log10(d$cec.C)





library(reshape)
library(ggplot2)
meltp <- unique(melt(d,id.vars = c("idp")))

ggplot(data = meltp,
       aes(x = value, fill=variable)) + geom_histogram() + 
  facet_wrap( ~ variable, scales = "free_x")
d$H <- NA
name(d)
e <- data.frame(d[,c(1,2,5,8,11:19)])
e$H <- "A"
name(e)
names(e)[2:4] <- c("clay.B","CEC.B","OC.B")
e <- rbind(e,d[,c(1,3,6,9,11:20)])
name(e)
e$H[is.na(e$H)] <- "B"

names(e)[2:4] <- c("clay.C","CEC.C","OC.C")
e <- rbind(e,d[,c(1,4,7,10,11:20)])
e$H[is.na(e$H)] <- "C"
names(e)[2:4] <- c("Clay","CEC","OC")

meltp <- melt(e,id.vars = c("idp","H"))
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
# statistics
round(stat.desc(D,norm = TRUE),0)

#### correlogram
library(psych)
corr <- d[,1:10]
par(mfrow=c(1,1),pty = "s")
pairs.panels(d[2:10], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )#, cex= 1)  # plot correlogram
## correlogram ####
library(GGally)
library(ggplot2)

ggscatmat(d[c(-192,-188),2:10], columns = 1:9,  alpha=0.15)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


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
CEC.A ~~ 0.05 * CEC.A
CEC.B ~~ 0.05 * CEC.B
CEC.C ~~ 0.05 * CEC.C
OC.A ~~ 0.05 * OC.A
OC.B ~~ 0.05 * OC.B
OC.C ~~ 0.05 * OC.C
clay.A ~~ 0.05 * clay.A
clay.B ~~ 0.05 * clay.B
clay.C ~~ 0.05 *clay.C

#--------------------#

# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn + X + 0*Y 
clay.Ar ~ clay.Cr + 
evisd + lstm + ndwi.b 
clay.Br ~ clay.Ar + clay.Cr + 
vdchn + twi + 0*Y + ndwi.b

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
# # 
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
mod[mod$mi>3 & (mod$op == "~"|mod$op == "~"),] 


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
# for(i in seq_along(D[,1])){ 
#   cal <- D[,]
#   pre[i,] <- D[ i,]
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = D, fixed.x = T,
                      estimator = "ML")
  # Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #
  # Matrix of Beta coefficients
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9]
  # Identity matrix (Kappa coefficients)
  I <- diag(nrow = 9, ncol = 9)
  # Matrix of Gamma coefficients
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:18]
  # Matrix of Psi coefficients (model error variance-covariance)
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9]
  # Matrix of measurement error (Epsylon)
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9]
  IB <- solve(I - B)
  # Running Prediction @ i location #
  # p is a matrix with the 10 external drivers
for(i in seq_along(D[,1])){ 
  pre[i,] <- D[ i,]
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 9, ncol = 1)           # by lavaan sequence
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

#########################################
### Autocorrelation in SEM residuals ####

Res[2:10] <- Res[,2:10]-Res[,11:19]
Res <- Res[,1:10]
names(Res)[1] <- "id.p"

R <- merge(Res, unique(d[,c(1,26,27)], by = "id.p"))

library(sp)
library(gstat)
# R as geo
coordinates(R) <- ~X+Y

# independent of device size
g <- list()
vg <- list()
vgm <- list()

# CEC
for(i in 2:4){
  g[[i-1]] <-  gstat(id = c(names(R@data)[i]), formula = formula(paste0(names(R@data)[i],"~1")),
            data = R)
  vg[[i-1]] <-  variogram(g[[i-1]])
  # vg = variogram(g, width = 20000, cutoff = 600000)
  # vg = variogram(g, boundaries = c(1E4,2E4,4E4,7E4,1E5,2E5,4E5,7E5,1E6))
  # print(plot(vg, plot.numbers = TRUE, main = names(R@data)[i]))
  
  # # choose initial variogram model and plot:
  vgm[[i-1]] <- vgm(nugget = var(R@data[,i]),
             psill=var(R@data[,i]),
             range=5E4,
             model = "Sph")
  #plot(vg, vgm, main = names(R@data)[i])
  # 
  # # fit variogram model:
  vgm[[i-1]] <-  fit.variogram(vg[[i-1]], vgm[[i-1]], fit.method = 2)
  print(names(R@data)[i])
  print(vgm[[i-1]])
  print(attr(vgm[[i-1]], "SSErr"))
}
# OC
for(i in 5:7){
  g[[i-1]] <-  gstat(id = c(names(R@data)[i]), formula = formula(paste0(names(R@data)[i],"~1")),
                     data = R)
  vg[[i-1]] <-  variogram(g[[i-1]])
  # vg = variogram(g, width = 20000, cutoff = 600000)
  # vg = variogram(g, boundaries = c(1E4,2E4,4E4,7E4,1E5,2E5,4E5,7E5,1E6))
  # print(plot(vg, plot.numbers = TRUE, main = names(R@data)[i]))
  
  # # choose initial variogram model and plot:
  vgm[[i-1]] <- vgm(nugget = var(R@data[,i]),
                    psill=var(R@data[,i]),
                    range=5E4,
                    model = "Sph")
  #plot(vg, vgm, main = names(R@data)[i])
  # 
  # # fit variogram model:
  vgm[[i-1]] <-  fit.variogram(vg[[i-1]], vgm[[i-1]], fit.method = 1)
  print(names(R@data)[i])
  print(vgm[[i-1]])
  print(attr(vgm[[i-1]], "SSErr"))
}
# Clay
for(i in 8:10){
  g[[i-1]] <-  gstat(id = c(names(R@data)[i]), formula = formula(paste0(names(R@data)[i],"~1")),
                     data = R)
  vg[[i-1]] <-  variogram(g[[i-1]])
  # vg = variogram(g, width = 20000, cutoff = 600000)
  # vg = variogram(g, boundaries = c(1E4,2E4,4E4,7E4,1E5,2E5,4E5,7E5,1E6))
  # print(plot(vg, plot.numbers = TRUE, main = names(R@data)[i]))
  
  # # choose initial variogram model and plot:
  vgm[[i-1]] <- vgm(nugget = var(R@data[,i]),
                    psill=var(R@data[,i]),
                    range=5E4,
                    model = "Sph")
  #plot(vg, vgm, main = names(R@data)[i])
  # 
  # # fit variogram model:
  vgm[[i-1]] <-  fit.variogram(vg[[i-1]], vgm[[i-1]], fit.method = 2)
  print(names(R@data)[i])
  print(vgm[[i-1]])
  print(attr(vgm[[i-1]], "SSErr"))
}

par(mfrow = c(1, 3), pty = "s")
plot(variogramLine(vgm[[1]], maxdist=50000), type="l", lwd=2,col="#AA0000", 
     main="CEC",
     xlab = "Distance", ylab = "Semivariance", cex.lab = 1.3, ylim=c(5,35))
legend(x= "topleft",legend = "SSErr", bty = "n")
legend(x= 12,legend = round(attr(vgm[[1]], "SSErr"),3), text.col ="#AA0000", bty = "n")
lines(variogramLine(vgm[[2]], maxdist=50000), lwd=2, col="#00AA00")
legend(x= 17,legend = round(attr(vgm[[2]], "SSErr"),3), text.col="#00AA00", bty = "n")
lines(variogramLine(vgm[[3]], maxdist=50000), lwd=2, col="#0000AA")
legend(x= 31,legend = round(attr(vgm[[3]], "SSErr"),3), text.col="#0000AA", bty = "n")

plot(variogramLine(vgm[[4]], maxdist=50000), type="l", lwd=2,col="#AA0000", 
     main="OC",
     xlab = "Distance", ylab = "Semivariance", cex.lab = 1.3, ylim=c(0,0.5)) 
legend(x= "topleft",legend = "SSErr", bty = "n")
legend(x= 0.25,legend = round(attr(vgm[[4]], "SSErr"),3), text.col ="#AA0000", bty = "n")
lines(variogramLine(vgm[[5]], maxdist=50000), lwd=2, col="#00AA00")
legend(x= 0.07,legend = round(attr(vgm[[5]], "SSErr"),3), text.col ="#00AA00", bty = "n")
lines(variogramLine(vgm[[6]], maxdist=50000), lwd=2, col="#0000AA")
legend(x= 0.016,legend = round(attr(vgm[[6]], "SSErr"),3), text.col ="#0000AA", bty = "n")

plot(variogramLine(vgm[[7]], maxdist=50000), type="l", lwd=2,col="#AA0000", 
     main="Clay",
     xlab = "Distance", ylab = "Semivariance", cex.lab = 1.3, ylim=c(10,50)) 
legend(x= "topleft",legend = "SSErr", bty = "n")
legend(x= 19,legend = round(attr(vgm[[7]], "SSErr"),3), text.col ="#AA0000", bty = "n")
lines(variogramLine(vgm[[8]], maxdist=50000), lwd=2, col="#00AA00")
legend(x= 29,legend = round(attr(vgm[[8]], "SSErr"),3), text.col ="#00AA00", bty = "n")
lines(variogramLine(vgm[[9]], maxdist=50000), lwd=2, col="#0000AA")
legend(x= 44,legend = round(attr(vgm[[9]], "SSErr"),3), text.col ="#0000AA", bty = "n")
par(mfrow = c(1, 1))

title("Semivariograms of residuals",
  sub="red: A horizon, green: B horizon, blue: C horizon")


vgm[[4]] <- vgm(nugget = 0.1,
                  psill=0.1,
                  range=5E4,
                  model = "Wav")
plot(vg[[4]], vgm[[4]],ylim=c(0,0.5))
vgm[[4]] <-  fit.variogram(vg[[4]], vgm[[4]], fit.method = 2)
plot(vg[[4]], vgm[[4]],ylim=c(0,0.5))
attr(vgm[[4]], "SSErr")
vgm[[4]]
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

library(lattice)
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

### CROSS-VALIDATION Multivariate LR ####
# comparison with multivariate 
attach(D)
mod.sp <- lm(cbind(CEC.A, CEC.B, CEC.C,
                   OC.A, OC.B, OC.C, 
                   clay.A, clay.B, clay.C) ~
               dem + river + vdchn + X + Y + evisd + lstm + ndwi.b +
               twi + ndwi.a)

var(mod.sp$residuals)
diag(psi.lm)[1:9] - diag(var(mod.sp$residuals))[1:9]

psi.lm[1:9,1:9] <- var(mod.sp$residuals)[1:9,1:9]

IB.lm.inv <- solve(diag(nrow(beta.lm)) - beta.lm)
lambda.lm <- diag( nrow(beta.lm) )
Sigma.hat.lm <- lambda.lm %*% IB.lm.inv %*% psi.lm %*% t(IB.lm.inv) %*% t(lambda.lm) 
rownames(Sigma.hat.lm) <- colnames(Sigma.hat)
colnames(Sigma.hat.lm) <- colnames(Sigma.hat)
Sigma.hat.lm[lower.tri(Sigma.hat.lm)] <- NA

library(lattice)
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


# summary(mod.sp)
# class(mod.sp)
# summary(mod.ceca)
# summary(manova(mod.sp))
# library(car)
# Manova(mod.sp, type = "II")
# mod.ceca$residuals - mod.sp$residuals[,1] #same

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


reportMLR$mean_theta[1:9] <- theta.MLR.mean[2:10]
reportMLR$median_th[1:9] <- theta.MLR.median[2:10]

round(reportMLR,3)
#d.stat <- read.csv("summary.calibdata.csv")
STt <- as.data.frame(STt)
STt$SS <- NA 
for(i in seq_along(names(d))){
  STt$SS[i] <- sum(( d[i] - STt$mean[i])^2)
}

reportMLR$R2 <- 1 - (as.numeric(reportMLR$SS) / as.numeric(STt$SS[2:10]))
reportMLR

write.csv(reportMLR, "~/Documents/SEM2DSM1/Paper_2/reports/reportMLR.csv")
################################################################################
# Predicting soil properties with mean values ####
# Cross-validation
# P <- predicted, C <- calibration, V <- validation, R <- residuals
P <- as.data.frame(D[1:10])
P[,11:19] <- NA
names(P)[11:19] <- paste0(names(D)[2:10],".pred")
attach(C)
V <- D[,2:10]
for(i in seq_along(D[,1])){
  C <- D[-i,]
  # calibration LOO
  P[i,11:19] <- sapply(C[,2:10],mean)
  V[i,] <- sapply(C[,2:10],var)
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
theta.MEAN.mean <- rep(NA,10)
theta.MEAN.median <- rep(NA,10)
R <- P[,1:10]
R[2:10] <- P[,2:10]-P[,11:19]

for(i in 2:10){
  theta.MEAN.mean[i] <- mean(R[,i]^2)
  theta.MEAN.median[i] <- median(R[,i]^2)
}


P[2:10] <- unstd(P[,2:10], STt[2:10,])
P[11:19] <- unstd(P[,11:19], STt[2:10,])
R <- P[,1:10]
R[2:10] <- P[,2:10]-P[,11:19]

# Accuracy measures of MLR
# create report

reportMEAN <- data.frame(Soil_property = NA, ME = NA, RMSE = NA, SS = NA,
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
  reportMEAN[i-1,1:4] <- c(names(P)[i], ME, RMSE, SS)
}


reportMEAN$mean_theta[1:9] <- theta.MEAN.mean[2:10]
reportMEAN$median_th[1:9] <- theta.MEAN.median[2:10]

#d.stat <- read.csv("summary.calibdata.csv")
STt <- as.data.frame(STt)
STt$SS <- NA 
for(i in seq_along(names(d))){
  STt$SS[i] <- sum(( d[i] - STt$mean[i])^2)
}

reportMEAN$R2 <- 1 - (as.numeric(reportMEAN$SS) / as.numeric(STt$SS[2:10]))
reportMEAN


CEC.v <- rbind(as.matrix(P[,c(2,2+9)]),
              as.matrix(P[,c(3,3+9)]),
              as.matrix(P[,c(4,4+9)]))
colnames(CEC.v) <- c("CECo","CECp")
CEC.v <- as.data.frame(CEC.v)
rsq.v <- data.frame(CEC=NA,OC=NA,clay=NA)
rsq.v[,1] <- 1 - (sum((CEC.v$CECo - CEC.v$CECp)^2)/
                   sum((mean(CEC.v$CECo)-CEC.v$CECo)^2))
lim = round(c(min(c(CEC.v[,1],CEC.v[,2])), max(c(CEC.v[,1],CEC.v[,2]))))
plot(CEC.v[,2]~CEC.v[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "CEC residuals", col = "dark red")
abline(0,1)
abline(lm(CEC.v[,2]~CEC.v[,1]),col = "blue")

OC.v <- rbind(as.matrix(P[,c(5,5+9)]),
              as.matrix(P[,c(6,6+9)]),
              as.matrix(P[,c(7,7+9)]))
colnames(OC.v) <- c("OCo","OCp")
OC.v <- as.data.frame(OC.v)
rsq.v[,2] <- 1 - (sum((OC.v$OCo - OC.v$OCp)^2)/
                    sum((mean(OC.v$OCo)-OC.v$OCo)^2))
lim = round(c(min(c(OC.v[,1],OC.v[,2])), max(c(OC.v[,1],OC.v[,2]))))
plot(OC.v[,2]~OC.v[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "OC residuals", col = "dark red")
abline(0,1)
abline(lm(OC.v[,2]~OC.v[,1]),col = "blue")

clay.v <- rbind(as.matrix(P[,c(8,8+9)]),
                as.matrix(P[,c(9,9+9)]),
                as.matrix(P[,c(10,10+9)]))
colnames(clay.v) <- c("clayo","clayp")
clay.v <- as.data.frame(clay.v)
rsq.v[,3] <- 1 - (sum((clay.v$clayo - clay.v$clayp)^2)/
                    sum((mean(clay.v$clayo)-clay.v$clayo)^2))
lim = round(c(min(c(clay.v[,1],clay.v[,2])), max(c(clay.v[,1],clay.v[,2]))))
plot(clay.v[,2]~clay.v[,1], xlim = lim, ylim= lim, xlab = "measured",
     ylab = "predicted", main = "clay residuals", col = "dark red")
abline(0,1)
abline(lm(clay.v[,2]~clay.v[,1]),col = "blue")

rsq.v

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







