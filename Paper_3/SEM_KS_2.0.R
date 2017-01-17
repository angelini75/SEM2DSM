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
rm(list=ls())
name <- function(x) { as.data.frame(names(x))}
# chose one
# Libraries ####
library(lavaan)
library(pastecs)
library(utils)

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



# library(reshape)
# library(ggplot2)
# meltp <- unique(melt(d,id.vars = c("idp")))
# 
# ggplot(data = meltp,
#        aes(x = value, fill=variable)) + geom_histogram() + 
#   facet_wrap( ~ variable, scales = "free_x")
# d$H <- NA
# name(d)
# e <- data.frame(d[,c(1,2,5,8,11:19)])
# e$H <- "A"
# name(e)
# names(e)[2:4] <- c("clay.B","CEC.B","OC.B")
# e <- rbind(e,d[,c(1,3,6,9,11:20)])
# name(e)
# e$H[is.na(e$H)] <- "B"
# 
# names(e)[2:4] <- c("clay.C","CEC.C","OC.C")
# e <- rbind(e,d[,c(1,4,7,10,11:20)])
# e$H[is.na(e$H)] <- "C"
# names(e)[2:4] <- c("Clay","CEC","OC")
# 
# meltp <- melt(e,id.vars = c("idp","H"))
# ggplot(data = meltp[meltp$variable=="CEC" |
#                       meltp$variable=="OC" |
#                       meltp$variable=="Clay",],
#        aes(x = value, fill = H)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable,scales = "free")
# ggplot(data = unique(meltp[!(meltp$variable=="CEC" |
#                       meltp$variable=="OC" |
#                       meltp$variable=="Clay"),]),
#        aes(x = value)) + geom_density(alpha = 0.4) + 
#   facet_wrap( ~ variable,scales = "free")
#############

# # New statistics
# d <- d[,-20:-21]
# names(d)[12:13] <- c("twi", "vdchn")
# round(stat.desc(d,norm = TRUE),3)
# # New mean and sd
# STt <- t(stat.desc(d,norm = TRUE)[c(9,13),])
# 
# # standardised data set ####
# std <- function(x, st){
#   y <- x
#   for(i in seq_along(names(x))){
#     y[,i] <- (x[,i] - st[i,1]) / st[i,2]
#   }
#   y
# }
# 
# D <- std(d,STt)
# 
# D[,1] <- d[,1] 
# # statistics
# round(stat.desc(D,norm = TRUE),0)

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
# ## correlogram ####
# library(GGally)
# library(ggplot2)
# 
# ggscatmat(d[c(-192,-188),2:10], columns = 1:9,  alpha=0.15) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


# SEM ####


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
# CEC.A ~~ 0.1 * CEC.A
# CEC.B ~~ 0.1 * CEC.B
# CEC.C ~~ 0.1 * CEC.C
# OC.A ~~ 0.1 * OC.A
# OC.B ~~ 0.1 * OC.B
# OC.C ~~ 0.1 * OC.C
# clay.A ~~ 0.05 * clay.A
# clay.B ~~ 0.05 * clay.B
# clay.C ~~ 0.05 *clay.C

#--------------------#

# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn + twi
clay.Ar ~ clay.Cr + 
          evisd + lstm + ndwi.b + twi
clay.Br ~ clay.Ar + clay.Cr + 
          vdchn + twi + ndwi.b

OC.Ar ~ clay.Ar +
        evisd + lstm + ndwi.b  
OC.Br ~ OC.Ar + clay.Br + 
        evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ clay.Br + OC.Br
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
clay.Br  ~     dem
OC.Ar  ~ clay.Br
OC.Cr  ~   OC.Ar

CEC.Ar ~~ clay.Ar
CEC.Br ~~   OC.Br

CEC.Cr  ~ ndwi.a
CEC.Br  ~ ndwi.a
CEC.Br  ~     dem
CEC.Ar  ~     dem

#------------------#
'
# Model calibration ####
my.fit.lv.ML <- sem(model = my.model.lv,data = D, meanstructure = FALSE, 
                    fixed.x = T)
# inspect(my.fit.lv.ML, "cov.lv")
# inspect(my.fit.lv.ML, "theta")
# inspect(my.fit.lv.ML, "est")
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
mod[mod$mi>6 & (mod$op == "~"|mod$op == "~~") ,]#& (mod$lhs == "CEC.Br"),] 

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
  #Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = cal, meanstructure = FALSE, 
                      fixed.x = T)
  # Matrix dedinition (Section 3.3 2nd paper and Fig. 5) #
  # Matrix of Beta coefficients
  B <- inspect(my.fit.lv.ML, "est")$beta[1:9,1:9]
  # Identity matrix (Kappa coefficients)
  I <- diag(nrow = 9, ncol = 9)
  # Matrix of Gamma coefficients
  A <- inspect(my.fit.lv.ML, "est")$beta[1:9,10:16]
  # Matrix of Psi coefficients (model error variance-covariance)
  V <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9]
  # Matrix of measurement error (Epsylon)
  Th <- inspect(my.fit.lv.ML, "est")$theta[1:9,1:9]
  IB <- solve(I - B)
  # Running Prediction @ i location #
  # p is a matrix with the 10 external drivers
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
# Soil_property     R2          R2MrLR
# 1         CEC.A   0.19144818  0.2404
# 2         CEC.B   0.08929689  0.1849
# 3         CEC.C   0.12246822  0.2003
# 4          OC.A   0.26613067  0.3061
# 5          OC.B   0.07876769  0.1002
# 6          OC.C   0.09265098  0.1224
# 7        clay.A   0.21246431  0.2649
# 8        clay.B   0.11342301  0.2023
# 9        clay.C   0.09477576  0.1823

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

#### MrLR ####
mod.sp <- lm(cbind(CEC.A, CEC.B, CEC.C,
                   OC.A, OC.B, OC.C, 
                   clay.A, clay.B, clay.C) ~
               dem + vdchn + X + Y + evisd + lstm + ndwi.b +
               twi + ndwi.a, D)

summary(mod.sp)





