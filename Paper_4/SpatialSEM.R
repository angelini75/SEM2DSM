###########
# Model 4 ####
# from Arg2Ks_6models.Rmd

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
d <- read.csv("KS.data-0.2.csv")[,c(-1)] 
name(d)
names(d)[2:10] <- c("Clay.A", "Clay.B", "Clay.C",
                    "CEC.A","CEC.B","CEC.C",
                    "OC.A","OC.B","OC.C")
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
#write.csv(STt.ks, "/home/marcos/Documents/SEM2DSM1/Paper_4/data/STt.ks.csv")

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

# remove samples with equal coordinates ####
ks[order(ks[,c(20)]),c(1,20,21)]

# X is the difference between samples and mean of samples
X = sqrt((sapply(X = ks[c(103,105,106),c(2:10)],FUN = median) - 
            ks[c(103,105,106),c(2:10)])^2)
# sum of difference
sum(X[1,])
sum(X[2,])
sum(X[3,])
# 106 is the most similar to median
ks <- ks[c(-103, -105),]

# again
ks[order(ks[,c(20)]),c(1,20,21)]
X = sqrt((sapply(X = ks[c(100,102),c(2:10)],FUN = mean) - 
            ks[c(100,102),c(2:10)])^2)
# sum of difference
sum(X[1,])
sum(X[2,])
# 106 is the most similar to median
ks <- ks[c(-102),]
#write.csv(ks,"/home/marcos/Documents/SEM2DSM1/Paper_4/data/ks.csv")
### END ###
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
# CEC.A ~~ 0.12 * CEC.A
# CEC.B ~~ 0.25 * CEC.B
# CEC.C ~~ 0.3 * CEC.C
OC.A ~~ 0.1 * OC.A
OC.B ~~ 0.1 * OC.B
OC.C ~~ 0.1 * OC.C
clay.A ~~ 0.05 * clay.A
clay.B ~~ 0.05 * clay.B
clay.C ~~ 0.05 *clay.C

#--------------------#
# Structural model (gamma and betta matrices)
#--------------------#
clay.Cr ~ dem + vdchn + X + lstm 
clay.Ar ~ clay.Cr + evisd + lstm + ndwi.b #+ Y 
clay.Br ~ clay.Ar + clay.Cr + vdchn + twi + ndwi.b + Y

OC.Ar ~ clay.Ar + evisd + lstm + ndwi.b 
OC.Br ~ OC.Ar + clay.Br + evisd + lstm + ndwi.a + vdchn
OC.Cr ~ OC.Br 

CEC.Ar ~ OC.Ar + clay.Ar 
CEC.Br ~ clay.Br + 0*OC.Br
CEC.Cr ~ clay.Cr + 0*OC.Cr

#------------------#
# Model error covariance (Psi)
#------------------#
CEC.Ar ~~ CEC.Br + CEC.Cr
CEC.Cr ~~ CEC.Br
#OC.Cr ~~ 0*CEC.Br + 0*CEC.Cr + 0*CEC.Ar 

#------------------#
# lavaan suggestions
#------------------#
clay.Br  ~    lstm
clay.Ar ~ twi

CEC.Ar ~~ clay.Br
OC.Ar ~   clay.Br
OC.Br ~   dem
OC.Cr ~~   clay.Br

CEC.Cr  ~  ndwi.a
CEC.Br  ~  ndwi.a + dem
CEC.Ar  ~    dem
#------------------#
'


# Model calibration ####
my.fit.lv.ML <- sem(model = my.model.lv,data = ks, meanstructure = FALSE, 
                    fixed.x = T)
inspect(my.fit.lv.ML,"cov.lv")
# Model evaluation ####
summary(my.fit.lv.ML, fit.measures=TRUE, rsquare = T)
mod.ks <- modindices(my.fit.lv.ML,sort. = T)
mod.ks[mod.ks$mi>7 & (mod.ks$op == "~"|mod.ks$op == "~~"),]# & mod.ks$lhs == "clay.Cr",] 


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
  cal <- ks[,] # not LOO. Set [-i,] for LOO 
  pre[i,] <- ks[i,]
  # Fiting #
  my.fit.lv.ML <- sem(model = my.model.lv,data = cal, fixed.x = T,
                      estimator = "ML")
  
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
  p = as.vector(as.matrix(pre[i,colnames(A)])) # values of covariates ordered
  p = matrix(p, nrow = 9, ncol = 1)           # by lavaan sequence
  # prediction
  pre[i,22:30] = t(IB %*% A %*% p) # key equation
  # calculate standarised squared standard error
  ## theta is standarised squared standard error
  ## theta = ((observed-predicted)^2)/error variance=(standard error^2)
  # a[i,] <- pre[i,c(2:10)] # observed values
  # b[i,] <- pre[i,c(22:30)] # predicted values
  v <- IB%*%V%*%t(IB)+Th # error variance (it is not diagonal!)
  # resids[i,] <- a[i,] - b[i,] # residuals
  # theta[i,] <- (resids[i,]^2)/v # theta
  # Error variance #
  Var[i,] <- diag(IB %*% V %*% t(IB))
}

# Model variance

# summary(Var)
# Var <- apply(Var, MARGIN = 2, FUN = mean)
# 
# # function to unstandardise the data
# unstd<- function(x, st){
#   y <- x
#   for(i in seq_along(names(x))){
#     y[,i] <- (x[,i] * st[i,2]) + st[i,1]
#   }
#   y
# }

# Accuracy measures ####
# Unstandardized residuals #
# Res <- cbind(pre[,1], unstd(pre[,2:10], STt.ks[2:10,]), unstd(pre[,22:30],
#                                                               STt.ks[2:10,]))

# Standardized residuals #
Res <-  cbind(pre[,1], pre[,2:10] - pre[,22:30])

(var.Res <- var(Res[2:10]))

(psi <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9])

(IB%*%V%*%t(IB)+Th) - var.Res

Res
#################################################
### Autocorrelation in SEM residuals ####

names(Res)[1] <- "id.p"

R <- merge(x = Res, y = unique(d[,c(1,20,21)]), by.x = "id.p", by.y="idp")




library(sp)
library(gstat)
# R as geo
coordinates(R) <- ~X+Y

#define crs
# wgs84 <- CRS("+init=epsg:4326")
# UTM14N <- CRS("+init=epsg:32614")
# modis <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
NAD83.KS.N <- CRS("+init=epsg:2796")

# Assign projection
proj4string(R) <- NAD83.KS.N

# check zero distance between profiles
zerodist(R, zero=0.0)
#should be zero
# R <- remove.duplicates(R, zero = 0.0, remove.second = TRUE)

#reduce ####
# g <- list()
# vg <- list()
# vgm <- list()
# par(mfrow = c(3, 3), pty = "s", mar=c(4,5,2,2), family="serif")
# # CEC.A
# g[[1]] <-  gstat(formula = CEC.A ~ 1, data = R)
# vg[[1]] <- variogram(g[[1]], width = 7000, cutoff = 450000, cressie = TRUE)
# plot(vg[[1]], plot.numbers = TRUE)
# 
# # # choose initial variogram model and plot:
# vgm[[1]] <- vgm(nugget = 20,
#                 psill= 1,
#                 range=100000,
#                 model = "Exp")
# vgm[[1]] <- fit.variogram(vg[[1]], vgm[[1]], fit.method = 6)
# plot(vg[[1]], vgm[[1]], main = "CEC.A")
# attr(vgm[[1]], "SSErr")
# 
# # CEC.B
# g[[2]] <-  gstat(formula = CEC.B ~ 1, data = R)
# vg[[2]] <- variogram(g[[2]], width = 7000, cutoff = 400000, cressie = TRUE)
# plot(vg[[2]], plot.numbers = TRUE)
# 
# # # choose initial variogram model and plot:
# vgm[[2]] <- vgm(nugget = 5,
#                 psill= 30,
#                 range=100000,
#                 model = "Exp")
# vgm[[2]] <- fit.variogram(vg[[2]], vgm[[2]], fit.method = 7)
# plot(vg[[2]], vgm[[2]], main = "CEC.B")
# attr(vgm[[2]], "SSErr")
# 
# # CEC.C
# g[[3]] <-  gstat(formula = CEC.C ~ 1, data = R)
# vg[[3]] <- variogram(g[[3]], width = 7000, cutoff = 400000, cressie = TRUE)
# plot(vg[[3]], plot.numbers = TRUE)
# 
# # # choose initial variogram model and plot:
# vgm[[3]] <- vgm(nugget = 5,
#                 psill= 30,
#                 range=50000,
#                 model = "Sph")
# vgm[[3]] <- fit.variogram(vg[[3]], vgm[[3]], fit.method = 7)
# plot(vg[[3]], vgm[[3]], main = "CEC.C")
# attr(vgm[[3]], "SSErr")
# vgm[[3]]
# 
# # OC.A
# g[[4]] <-  gstat(formula = OC.A ~ 1, data = R)
# vg[[4]] <- variogram(g[[4]], width = 15000, cutoff = 400000, cressie = TRUE)
# plot(vg[[4]], plot.numbers = TRUE)
# 
# vgm[[4]] <- vgm(nugget = 0.05,
#                 psill= 0.4,
#                 range=50000,
#                 model = "Exp")
# vgm[[4]] <- fit.variogram(vg[[4]], vgm[[4]], fit.method = 7)
# plot(vg[[4]], vgm[[4]], main = "OC.A")
# attr(vgm[[4]], "SSErr")
# vgm[[4]]
# 
# # OC.B
# g[[5]] <-  gstat(formula = OC.B ~ 1, data = R)
# vg[[5]] <- variogram(g[[5]], width = 15000, cutoff = 400000, cressie = TRUE)
# plot(vg[[5]], plot.numbers = TRUE)
# 
# vgm[[5]] <- vgm(nugget = 0,
#                 psill= 0.8,
#                 range=20000,
#                 model = "Exp")
# vgm[[5]] <- fit.variogram(vg[[5]], vgm[[5]], fit.method = 7)
# plot(vg[[5]], vgm[[5]], main = "OC.B")
# attr(vgm[[5]], "SSErr")
# vgm[[5]]
# 
# 
# # OC.C
# g[[6]] <-  gstat(formula = OC.C ~ 1, data = R)
# vg[[6]] <- variogram(g[[6]], width = 10000, cutoff = 400000, cressie = TRUE)
# plot(vg[[6]], plot.numbers = TRUE)
# 
# vgm[[6]] <- vgm(nugget = 0.1,
#                 psill= 0.5,
#                 range=20000,
#                 model = "Gau")
# vgm[[6]] <- fit.variogram(vg[[6]], vgm[[6]], fit.method = 2)
# plot(vg[[6]], vgm[[6]], main = "OC.C")
# attr(vgm[[6]], "SSErr")
# vgm[[6]]
# 
# 
# # Clay.A
# g[[7]] <-  gstat(formula = clay.A ~ 1, data = R)
# vg[[7]] <- variogram(g[[7]], width = 7000, cutoff = 400000, cressie = TRUE)
# plot(vg[[7]], plot.numbers = TRUE)
# 
# vgm[[7]] <- vgm(nugget = 20,
#                 psill= 80,
#                 range=100000,
#                 model = "Exp")
# vgm[[7]] <- fit.variogram(vg[[7]], vgm[[7]], fit.method = 1)
# plot(vg[[7]], vgm[[7]], main = "Clay.A")
# attr(vgm[[7]], "SSErr")
# vgm[[7]]
# 
# # Clay.B
# g[[8]] <-  gstat(formula = clay.B ~ 1, data = R)
# vg[[8]] <- variogram(g[[8]], width = 7000, cutoff = 400000, cressie = TRUE)
# plot(vg[[8]], plot.numbers = TRUE)
# 
# vgm[[8]] <- vgm(nugget = 15,
#                 psill= 120,
#                 range=100000,
#                 model = "Sph")
# vgm[[8]] <- fit.variogram(vg[[8]], vgm[[8]], fit.method = 1)
# plot(vg[[8]], vgm[[8]], main = "Clay.B")
# attr(vgm[[8]], "SSErr")
# vgm[[8]]
# 
# # Clay.C
# g[[9]] <-  gstat(formula = clay.C ~ 1, data = R)
# vg[[9]] <- variogram(g[[9]], width = 15000, cutoff = 400000, cressie = TRUE)
# plot(vg[[9]], plot.numbers = TRUE)
# 
# vgm[[9]] <- vgm(nugget = 20,
#                 psill= 100,
#                 range=100000,
#                 model = "Gau")
# vgm[[9]] <- fit.variogram(vg[[9]], vgm[[9]], fit.method = 2)
# plot(vg[[9]], vgm[[9]], main = "Clay.C")
# attr(vgm[[9]], "SSErr")
# vgm[[9]]
# 
# 
# # Three graphs (soil properties)
# # tiff(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig1.tif", 
# #       width = 2500, height = 1000, res =  350)
# par(mfrow = c(1, 3), pty = "s", mar=c(4,5,2,2), family="serif")
# ### CEC
# ## A
# plot(variogramLine(vgm[[1]], maxdist=500000), 
#      type="l", lwd=2,col="#AA0000",
#      main= "CEC", 
#      xlab = "Distance / m", 
#      ylab = expression("Semivariance"~~"/"~("cmol"[c]~~"kg"^{-1})^{2}),
#      cex.lab = 1.3, ylim=c(0,1.2))
# points(gamma ~ dist, vg[[1]], col="#770000")
# #legend(x= "topleft",legend = "SSErr", bty = "n")
# #legend(x= 12,legend = round(attr(vgm[[1]], "SSErr"),3), text.col ="#AA0000", 
# #        bty = "n")
# ## B
# lines(variogramLine(vgm[[2]], maxdist=500000), lwd=2, col="#00AA00")
# points(gamma ~ dist, vg[[2]], col="#007700")
# ## C
# lines(variogramLine(vgm[[3]], maxdist=500000), lwd=2, col="#0000AA")
# points(gamma ~ dist, vg[[3]], col="#000077")
# 
# ### OC
# ## A
# plot(variogramLine(vgm[[4]], maxdist=500000), type="l", lwd=2,col="#AA0000", 
#      main= "OC",
#      xlab = "Distance / m", 
#      ylab = expression("Semivariance"~~"/ %"^{2}),
#      cex.lab = 1.3, ylim=c(0,1.2)) 
# points(gamma ~ dist, vg[[4]], col="#770000")
# ## B
# lines(variogramLine(vgm[[5]], maxdist=500000), lwd=2, col="#00AA00")
# points(gamma ~ dist, vg[[5]], col="#007700")
# ## C
# lines(variogramLine(vgm[[6]], maxdist=500000), lwd=2, col="#0000AA")
# points(gamma ~ dist, vg[[6]], col="#000077")
# 
# ### Clay
# ## A
# plot(variogramLine(vgm[[7]], maxdist=500000), type="l", lwd=2,col="#AA0000", 
#      main="Clay",
#      xlab = "Distance / m", 
#      ylab = expression("Semivariance"~~"/ %"^{2}),
#      cex.lab = 1.3, ylim=c(0,1.2)) 
# points(gamma ~ dist, vg[[7]], col="#770000")
# ## B
# lines(variogramLine(vgm[[8]], maxdist=500000), lwd=2, col="#00AA00")
# points(gamma ~ dist, vg[[8]], col="#007700")
# ## C
# lines(variogramLine(vgm[[9]], maxdist=500000), lwd=2, col="#0000AA")
# points(gamma ~ dist, vg[[9]], col="#000077")
# 
# #dev.off()

###### Cross-variograms ######

g <- list()
vg <- list()
vgm <- list()

# # CEC.A
# g[[10]] <-  gstat(id = c("CEC.A", "CEC.B", "CEC.C"), data = R)
# vg[[1]] <- variogram(g[[1]], width = 7000, cutoff = 450000, cressie = TRUE)
# plot(vg[[1]], plot.numbers = TRUE)
# 
# # # choose initial variogram model and plot:
# vgm[[1]] <- vgm(nugget = 20,
#                 psill= 1,
#                 range=100000,
#                 model = "Exp")
# vgm[[1]] <- fit.variogram(vg[[1]], vgm[[1]], fit.method = 6)
# plot(vg[[1]], vgm[[1]], main = "CEC.A")

names(R)[8:10] <- c("Clay.A", "Clay.B", "Clay.C") 
rm(cv)
cv <- gstat(id = "CEC.A", formula = CEC.A ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "CEC.B", formula = CEC.B ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "CEC.C", formula = CEC.C ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.A", formula = OC.A ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.B", formula = OC.B ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "OC.C", formula = OC.C ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.A", formula = Clay.A ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.B", formula = Clay.B ~ 1, data = R, nmax = 10)
cv <- gstat(cv, id = "Clay.C", formula = Clay.C ~ 1, data = R, nmax = 10)
cv <- gstat(cv, 
            model = vgm(nugget = 0.20,
                        psill= 1.1,
                        range=100000,
                        model = "Exp"), 
            fill.all = T)
cv
cv.var<- variogram(object = cv, cutoff = 450000) 

cv.fit<-fit.lmc(v = cv.var,g =  cv, fit.lmc = F, fit.ranges = F) 

# png(filename = "~/Dropbox/PhD Marcos/Paper 4/Figures/Fig2.png", 
     # width = 3000, height = 3000, res =  250)
plot(cv.var, model=cv.fit, 
     main="Variograms and cross-variograms of standardized residuals",
     xlab = "Distance / m", 
     ylab = "Semivariance",
     scales=list(x = list(alternating = 1), y = list(alternating = 1)),
     par.settings=list(grid.pars=list(fontfamily="serif")))
# dev.off()

Rdist <- spDists(R)
Rdist[1:8,1:8]
f <- names(cv$model)
p <- lower.tri(matrix(data = NA, nrow = 9, ncol = 9), diag = TRUE)
up <- upper.tri(matrix(data = NA, nrow = 9, ncol = 9), diag = FALSE)
p[p==FALSE] <- NA
place <- which(p)
C <- p
C0 <- p
a <- p
alpha <- p

cv$model[[1]]
for(i in seq_along(f)){
  a[place[i]] <- cv$model[[f[i]]][2,"range"]
  C[place[i]] <- cv$model[[f[i]]][2,"psill"] - cv$model[[f[i]]][1,"psill"]
  C0[place[i]] <- cv$model[[f[i]]][1,"psill"]
  alpha[place[i]] <- C[place[i]]/(C0[place[i]] + C[place[i]])
}
# this is Sigma fitted with lavaan
IB%*%V%*%t(IB)+Th # Th is the diag matr. of measurement errors

# The starting values of Sigma zero (Sigma0) depart from
# Matrix of Beta coefficients
(B <- lavTech(my.fit.lv.ML, "start")$beta[1:9,1:9])
# Identity matrix
I <- diag(nrow = 9, ncol = 9)
# Matrix of Gamma coefficients
A <- lavTech(my.fit.lv.ML, "start")$beta[1:9,10:18]
# Matrix of Psi coefficients (model error variance-covariance)
(Psi <- inspect(my.fit.lv.ML, "est")$psi[1:9,1:9])
# Matrix of measurement error (Epsylon)
Th <- lavTech(my.fit.lv.ML, "start")$theta[1:9,1:9]
(IB <- solve(I - B))

# Sigma0
(Sigma0 <- IB%*%V%*%t(IB)+Th)
# Make Sigma0 np x np
Sigma0.e <- kronecker(Y = Sigma0, X = matrix(1,nrow = 153,ncol = 153),FUN = "*")
round(Sigma0.e[1:12,1:12],2)
library(lattice)
levelplot(Sigma0.e[1:90,1:90])

# Create exp(-h/a)
a[which(up)] <- t(a)[upper.tri(a, F)]
exp <- exp(kronecker(X = -Rdist, Y = a,FUN = "/"))
levelplot(exp[1:90,1:90])

# Make alpha np x np
alpha[which(up)] <- t(alpha)[upper.tri(alpha, F)]
alpha.e <- kronecker(X = alpha, Y = matrix(1,nrow = 153,ncol = 153),FUN = "*")
levelplot(alpha.e[1:90,1:90]) # it is constant because we choose une single value for all cross variograms

# Let us call Rho = alpha * exp[-h/a]
Rho <- exp
for(i in seq_along(Rdist)) {
  if(exp[i] == 1){
    Rho[i] <- 1
  } else {
    Rho[i] <- alpha.e[i] * exp[i]
  }
}
levelplot(Rho[1:90,1:90])

Sigma <- Rho * Sigma0.e
levelplot(Sigma[1:90,1:90])


# png(filename = "/home/marcos/Desktop/Sigma.png", width = 3200, height = 3000, res = 500)
# levelplot(Sigma)
# dev.off()
